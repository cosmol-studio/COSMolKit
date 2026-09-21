//! Typed V3000 SGroup and enhanced-stereo collection lowering.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, BondId, SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle,
    SGroupCState, SGroupConnection, StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind, TopologyBlock,
};

use crate::sdf::{
    SdfReadError, SdfWriteError, get_v3000_line, parse_rdkit_atof, parse_rdkit_double,
    parse_rdkit_int, parse_rdkit_unsigned, rdkit_substr,
};

fn sgroup_kind_from_rdkit_type(value: &str) -> SubstanceGroupKind {
    match value {
        "DAT" => SubstanceGroupKind::Data,
        "SUP" => SubstanceGroupKind::Superatom,
        "MUL" => SubstanceGroupKind::MultipleGroup,
        "SRU" => SubstanceGroupKind::StructuralRepeatUnit,
        "MON" => SubstanceGroupKind::Monomer,
        "COP" => SubstanceGroupKind::Copolymer,
        "CRO" => SubstanceGroupKind::Crosslink,
        "GRA" => SubstanceGroupKind::Graft,
        "MOD" => SubstanceGroupKind::Modification,
        "MER" => SubstanceGroupKind::Mer,
        "ANY" => SubstanceGroupKind::AnyPolymer,
        "COM" => SubstanceGroupKind::MixtureComponent,
        "MIX" => SubstanceGroupKind::Mixture,
        "FOR" => SubstanceGroupKind::Formulation,
        other => SubstanceGroupKind::Generic(other.to_owned()),
    }
}

fn is_valid_rdkit_sgroup_type(value: &str) -> bool {
    matches!(
        value,
        "SRU"
            | "MON"
            | "COP"
            | "CRO"
            | "GRA"
            | "MOD"
            | "MER"
            | "ANY"
            | "COM"
            | "MIX"
            | "FOR"
            | "SUP"
            | "MUL"
            | "DAT"
            | "GEN"
    )
}

fn is_valid_rdkit_sgroup_subtype(value: &str) -> bool {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroupChecks::isValidSubType
    // RDKit✔️✔️: const std::vector<std::string> sGroupSubtypes = {"ALT", "RAN", "BLO"};
    // RDKit✔️✔️: bool SubstanceGroupChecks::isValidSubType(const std::string &type) {
    // RDKit✔️✔️:   return std::find(SubstanceGroupChecks::sGroupSubtypes.begin(),
    // RDKit✔️✔️:                    SubstanceGroupChecks::sGroupSubtypes.end(),
    // RDKit✔️✔️:                    type) != SubstanceGroupChecks::sGroupSubtypes.end();
    // RDKit✔️✔️: }
    matches!(value, "ALT" | "RAN" | "BLO")
    // The fixed three-value match is behaviorally identical to the source
    // vector search and has constant time with no allocation.
    // END RDKIT CPP FUNCTION
}

fn sgroup_connection_from_rdkit(value: &str) -> Option<SGroupConnection> {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroupChecks::isValidConnectType
    // RDKit✔️✔️: const std::vector<std::string> sGroupConnectTypes = {"HH", "HT", "EU"};
    // RDKit✔️✔️: bool SubstanceGroupChecks::isValidConnectType(const std::string &type) {
    // RDKit✔️✔️:   return std::find(SubstanceGroupChecks::sGroupConnectTypes.begin(),
    // RDKit✔️✔️:                    SubstanceGroupChecks::sGroupConnectTypes.end(),
    // RDKit✔️✔️:                    type) != SubstanceGroupChecks::sGroupConnectTypes.end();
    // RDKit✔️✔️: }
    match value {
        "HH" => Some(SGroupConnection::HeadToHead),
        "HT" => Some(SGroupConnection::HeadToTail),
        "EU" => Some(SGroupConnection::Either),
        _ => None,
    }
    // The source validates one of three exact strings and stores that string.
    // The canonical enum stores the same three states without a second raw
    // property. Matching is constant time and allocation-free.
    // END RDKIT CPP FUNCTION
}

fn is_valid_rdkit_sgroup_class(value: &str) -> bool {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroupChecks::isValidClass
    // RDKit✔️✔️: const std::vector<std::string> sGroupClasses = {
    // RDKit✔️✔️:     "AA",        "dAA",    "DNA",     "RNA",      "SUGAR",    "BASE",
    // RDKit✔️✔️:     "PHOSPHATE", "LINKER", "CHEM",    "LGRP",     "MODAA",    "MODdAA",
    // RDKit✔️✔️:     "MODDNA",    "MODRNA", "XLINKAA", "XLINKdAA", "XLINKDNA", "XLINKRNA",
    // RDKit✔️✔️: };
    // RDKit✔️✔️: bool SubstanceGroupChecks::isValidClass(const std::string &sgroupClass) {
    // RDKit✔️✔️:   return std::find(SubstanceGroupChecks::sGroupClasses.begin(),
    // RDKit✔️✔️:                    SubstanceGroupChecks::sGroupClasses.end(),
    // RDKit✔️✔️:                    sgroupClass) != SubstanceGroupChecks::sGroupClasses.end();
    // RDKit✔️✔️: }
    matches!(
        value,
        "AA" | "dAA"
            | "DNA"
            | "RNA"
            | "SUGAR"
            | "BASE"
            | "PHOSPHATE"
            | "LINKER"
            | "CHEM"
            | "LGRP"
            | "MODAA"
            | "MODdAA"
            | "MODDNA"
            | "MODRNA"
            | "XLINKAA"
            | "XLINKdAA"
            | "XLINKDNA"
            | "XLINKRNA"
    )
    // The fixed match has the exact case-sensitive source vocabulary, does
    // not allocate, and is no worse than the source's bounded linear search.
    // END RDKIT CPP FUNCTION
}

fn split_sgroup_line(
    line: &str,
    line_number: usize,
) -> Result<(&str, &str, &str, &str), SdfReadError> {
    let bytes = line.as_bytes();
    let mut position = 0_usize;
    let mut next_field = || {
        while bytes
            .get(position)
            .is_some_and(|byte| byte.is_ascii_whitespace())
        {
            position += 1;
        }
        let start = position;
        while bytes
            .get(position)
            .is_some_and(|byte| !byte.is_ascii_whitespace())
        {
            position += 1;
        }
        &line[start..position]
    };
    let sequence = next_field();
    let kind = next_field();
    let external_id = next_field();
    if sequence.is_empty() || kind.is_empty() || external_id.is_empty() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    Ok((sequence, kind, external_id, &line[position..]))
}

fn parse_sgroup_header_unsigned(
    text: &str,
    field: &str,
    line_number: usize,
) -> Result<u32, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (formatted header fields)
    // RDKit❗✔️: unsigned int sequenceId;
    // RDKit❗✔️: unsigned int externalId;
    // RDKit❗✔️: std::string type;
    // RDKit❗✔️:
    // RDKit❗✔️: std::stringstream lineStream(tempStr);
    // RDKit❗✔️: lineStream >> sequenceId;
    // RDKit❗✔️: lineStream >> type;
    // RDKit❗✔️: lineStream >> externalId;
    let (negative, digits) = match text.as_bytes().first() {
        Some(b'+') => (false, &text[1..]),
        Some(b'-') => (true, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(SdfReadError::Parse(format!(
            "Cannot convert '{text}' to unsigned {field} on line {line_number}"
        )));
    }
    let magnitude = digits.parse::<u32>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Unsigned {field} '{text}' is out of range on line {line_number}"
        ))
    })?;
    // Formatted unsigned extraction accepts a sign. For a representable
    // magnitude, a leading minus is assigned modulo the destination width;
    // fixed RDKit 2026.03.1 therefore observes `-1` as `UINT_MAX`.
    let value = if negative {
        0_u32.wrapping_sub(magnitude)
    } else {
        magnitude
    };
    // Behavioral review: complete signed decimal tokens match the pinned
    // formatted `unsigned int` assignment, including `+` and negative wrap.
    // The source leaves these destinations uninitialized on no conversion and
    // subsequently reads them, and overflow also poisons later extraction;
    // those C++ undefined/uninitialized cases are rejected structurally here
    // rather than assigned a fabricated zero. This boundary is not marked as
    // all-input source equivalence.
    // Complexity review: one linear validation scan and one bounded decimal
    // conversion match the source formatted extraction asymptotically, with
    // no allocation beyond an error message on failure.
    Ok(value)
    // END RDKIT CPP FUNCTION
}

#[derive(Clone, Copy)]
struct SGroupLineCursor<'a> {
    text: &'a str,
    position: usize,
    failed: bool,
}

impl<'a> SGroupLineCursor<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            text,
            position: 0,
            failed: false,
        }
    }

    fn get(&mut self) -> Option<u8> {
        if self.failed {
            return None;
        }
        let byte = self.text.as_bytes().get(self.position).copied()?;
        self.position += 1;
        Some(byte)
    }

    fn peek(&self) -> Option<u8> {
        (!self.failed)
            .then(|| self.text.as_bytes().get(self.position).copied())
            .flatten()
    }

    fn skip_c_locale_whitespace(&mut self) {
        while self
            .peek()
            .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
        {
            self.position += 1;
        }
    }

    fn read_unsigned(&mut self) -> Option<u32> {
        if self.failed {
            return None;
        }
        self.skip_c_locale_whitespace();
        let negative = self.peek() == Some(b'-');
        if matches!(self.peek(), Some(b'+') | Some(b'-')) {
            self.position += 1;
        }
        let start = self.position;
        let mut value = 0_u32;
        let mut overflow = false;
        while let Some(byte) = self.peek().filter(u8::is_ascii_digit) {
            overflow |= value
                .checked_mul(10)
                .and_then(|current| current.checked_add(u32::from(byte - b'0')))
                .map_or(true, |next| {
                    value = next;
                    false
                });
            self.position += 1;
        }
        if self.position == start {
            self.failed = true;
            return Some(0);
        }
        if overflow {
            self.failed = true;
            return Some(u32::MAX);
        }
        Some(if negative {
            value.wrapping_neg()
        } else {
            value
        })
    }

    fn read_double(&mut self) -> Option<f64> {
        // BEGIN LIBSTDC++ CPP FUNCTION num_get::_M_extract_float (C locale)
        // libstdc++✔️✔️: else if ((__c == __lit[__num_base::_S_ie]
        // libstdc++✔️✔️:           || __c == __lit[__num_base::_S_iE])
        // libstdc++✔️✔️:          && !__found_sci && __found_mantissa)
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     // Scientific notation.
        // libstdc++✔️✔️:     __xtrc += 'e';
        // libstdc++✔️✔️:     __found_sci = true;
        // libstdc++✔️✔️:
        // libstdc++✔️✔️:     // Remove optional plus or minus sign, if they exist.
        // libstdc++✔️✔️:     if (++__beg != __end)
        // libstdc++✔️✔️:       {
        // libstdc++✔️✔️:         __c = *__beg;
        // libstdc++✔️✔️:         const bool __plus = __c == __lit[__num_base::_S_iplus];
        // libstdc++✔️✔️:         if (__plus || __c == __lit[__num_base::_S_iminus])
        // libstdc++✔️✔️:           __xtrc += __plus ? '+' : '-';
        // libstdc++✔️✔️:         else
        // libstdc++✔️✔️:           continue;
        // libstdc++✔️✔️:       }
        // libstdc++✔️✔️:     else
        // libstdc++✔️✔️:       {
        // libstdc++✔️✔️:         __testeof = true;
        // libstdc++✔️✔️:         break;
        // libstdc++✔️✔️:       }
        // libstdc++✔️✔️:   }
        // END LIBSTDC++ CPP FUNCTION num_get::_M_extract_float (C locale)
        // BEGIN LIBSTDC++ CPP FUNCTION num_get::do_get(double&)
        // libstdc++✔️✔️: string __xtrc;
        // libstdc++✔️✔️: __xtrc.reserve(32);
        // libstdc++✔️✔️: __beg = _M_extract_float(__beg, __end, __io, __err, __xtrc);
        // libstdc++✔️✔️: std::__convert_to_v(__xtrc.c_str(), __v, __err, _S_get_c_locale());
        // libstdc++✔️✔️: if (__beg == __end)
        // libstdc++✔️✔️:   __err |= ios_base::eofbit;
        // libstdc++✔️✔️: return __beg;
        // END LIBSTDC++ CPP FUNCTION num_get::do_get(double&)
        // BEGIN LIBSTDC++ CPP FUNCTION __convert_to_v(double&)
        // libstdc++✔️✔️: char* __sanity;
        // libstdc++✔️✔️: __v = __strtod_l(__s, &__sanity, __cloc);
        // libstdc++✔️✔️:
        // libstdc++✔️✔️: // _GLIBCXX_RESOLVE_LIB_DEFECTS
        // libstdc++✔️✔️: // 23. Num_get overflow result.
        // libstdc++✔️✔️: if (__sanity == __s || *__sanity != '\0')
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     __v = 0.0;
        // libstdc++✔️✔️:     __err = ios_base::failbit;
        // libstdc++✔️✔️:   }
        // libstdc++✔️✔️: else if (__v == numeric_limits<double>::infinity())
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     __v = numeric_limits<double>::max();
        // libstdc++✔️✔️:     __err = ios_base::failbit;
        // libstdc++✔️✔️:   }
        // libstdc++✔️✔️: else if (__v == -numeric_limits<double>::infinity())
        // libstdc++✔️✔️:   {
        // libstdc++✔️✔️:     __v = -numeric_limits<double>::max();
        // libstdc++✔️✔️:     __err = ios_base::failbit;
        // libstdc++✔️✔️:   }
        // END LIBSTDC++ CPP FUNCTION __convert_to_v(double&)
        // Behavioral review: the scanner consumes an exponent marker and its
        // optional sign even when no exponent digit follows. The fixed GCC
        // 15.2 libstdc++ conversion then assigns zero and failbit; decimal
        // overflow assigns the signed finite maximum and failbit. Once failed,
        // the cursor models the stream sentry by refusing subsequent reads.
        // The existing raw-atof owner remains responsible for supported finite
        // decimal conversion. Exceptionally long significands, subnormal
        // boundaries, hexadecimal floats, and NaN payloads retain the approved
        // non-equivalence boundary and are not claimed here.
        // Complexity review: this is one allocation-free linear byte scan plus
        // the existing linear decimal conversion, matching formatted extraction
        // asymptotically without repeated scans or new buffering.
        if self.failed {
            return None;
        }
        self.skip_c_locale_whitespace();
        let start = self.position;
        if matches!(self.peek(), Some(b'+') | Some(b'-')) {
            self.position += 1;
        }
        let integer_start = self.position;
        while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
            self.position += 1;
        }
        let mut has_digit = self.position != integer_start;
        if self.peek() == Some(b'.') {
            self.position += 1;
            let fraction_start = self.position;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.position += 1;
            }
            has_digit |= self.position != fraction_start;
        }
        if !has_digit {
            self.failed = true;
            return Some(0.0);
        }
        if matches!(self.peek(), Some(b'e') | Some(b'E')) {
            self.position += 1;
            if matches!(self.peek(), Some(b'+') | Some(b'-')) {
                self.position += 1;
            }
            let exponent_digits = self.position;
            while self.peek().is_some_and(|byte| byte.is_ascii_digit()) {
                self.position += 1;
            }
            if self.position == exponent_digits {
                self.failed = true;
                return Some(0.0);
            }
        }
        let value = parse_rdkit_atof(&self.text[start..self.position]);
        if value.is_infinite() {
            self.failed = true;
            Some(value.signum() * f64::MAX)
        } else {
            Some(value)
        }
    }

    fn read_string_token(&mut self) -> Option<&'a str> {
        if self.failed {
            return None;
        }
        self.skip_c_locale_whitespace();
        let start = self.position;
        while self
            .peek()
            .is_some_and(|byte| !matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
        {
            self.position += 1;
        }
        if self.position == start {
            self.failed = true;
            None
        } else {
            Some(&self.text[start..self.position])
        }
    }

    fn read_label(&mut self) -> Option<&'a str> {
        if self.failed || self.position == self.text.len() {
            return None;
        }
        let start = self.position;
        while let Some(byte) = self.get() {
            if byte == b'=' {
                return Some(&self.text[start..self.position - 1]);
            }
        }
        (start != self.position).then_some(&self.text[start..self.position])
    }
}

fn parse_array<T: Copy>(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    max_count: Option<usize>,
    strict_parsing: bool,
    mut parse: impl FnMut(&mut SGroupLineCursor<'_>) -> Option<T>,
) -> Result<Vec<T>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000Array
    // RDKit✔️✔️: auto paren = stream.get();  // discard parentheses
    // RDKit✔️✔️: if (paren != '(') {
    // RDKit✔️✔️:   BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:       << "WARNING: first character of V3000 array is not '('" << std::endl;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int count = 0;
    // RDKit✔️✔️: stream >> count;
    // RDKit✔️✔️: std::vector<T> values;
    // RDKit✔️✔️: if (maxV >= 0 && count > static_cast<unsigned int>(maxV)) {
    // RDKit✔️✔️:   SGroupWarnOrThrow(strictParsing, "invalid count value");
    // RDKit✔️✔️:   return values;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: values.reserve(count);
    // RDKit✔️✔️: T value;
    // RDKit✔️✔️: for (unsigned i = 0; i < count; ++i) {
    // RDKit✔️✔️:   stream >> value;
    // RDKit✔️✔️:   values.push_back(value);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: paren = stream.get();  // discard parentheses
    // RDKit✔️✔️: if (paren != ')') {
    // RDKit✔️✔️:   BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:       << "WARNING: final character of V3000 array is not ')'" << std::endl;
    // RDKit✔️✔️: }
    // Missing parentheses are warnings only. The cursor deliberately consumes
    // the same bytes as the source stream. Formatted extraction assigns zero
    // to an arithmetic destination when no characters can be converted and
    // sets failbit; later sentry failures leave that zero in place. Overflow
    // assigns the corresponding saturated boundary and also sets failbit.
    let _opening_parenthesis = cursor.get();
    let count = cursor.read_unsigned().unwrap_or(0) as usize;
    if max_count.is_some_and(|maximum| count > maximum) {
        return if strict_parsing {
            Err(SdfReadError::Parse("invalid count value".to_owned()))
        } else {
            Ok(Vec::new())
        };
    }
    let mut values = Vec::with_capacity(count);
    let mut previous = None;
    for _ in 0..count {
        let value = parse(cursor).or(previous).ok_or_else(|| {
            SdfReadError::Parse(format!(
                "V3000 array has no initialized value on line {line_number}"
            ))
        })?;
        values.push(value);
        previous = Some(value);
    }
    let _closing_parenthesis = cursor.get();
    Ok(values)
    // END RDKIT CPP FUNCTION
}

fn parse_u32_array(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    max_count: Option<usize>,
    strict_parsing: bool,
) -> Result<Vec<u32>, SdfReadError> {
    parse_array(cursor, line_number, max_count, strict_parsing, |cursor| {
        cursor.read_unsigned()
    })
}

fn parse_f64_array(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    max_count: Option<usize>,
    strict_parsing: bool,
) -> Result<Vec<f64>, SdfReadError> {
    parse_array(cursor, line_number, max_count, strict_parsing, |cursor| {
        cursor.read_double()
    })
}

fn parse_string_property(cursor: &mut SGroupLineCursor<'_>) -> String {
    // BEGIN RDKIT CPP FUNCTION ParseV3000StringPropLabel
    // RDKit✔️✔️: std::string ParseV3000StringPropLabel(std::stringstream &stream) {
    // RDKit✔️✔️:   std::string strValue;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto nextChar = stream.peek();
    // RDKit✔️✔️:   if (nextChar == ' ') {
    // RDKit✔️✔️:     // empty value, we peeked at the next field's separator
    // RDKit✔️✔️:     return strValue;
    // RDKit✔️✔️:   } else if (nextChar == '"') {
    // RDKit✔️✔️:     // skip the opening quote:
    // RDKit✔️✔️:     stream.get();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // this is a bit gross because it's legal to include a \" in a value,
    // RDKit✔️✔️:     // but the way that's done is by doubling it. So
    // RDKit✔️✔️:     // FIELDINFO=""""
    // RDKit✔️✔️:     // should assign the value \" to FIELDINFO
    // RDKit✔️✔️:     char chr;
    // RDKit✔️✔️:     while (stream.get(chr)) {
    // RDKit✔️✔️:       if (chr == '"') {
    // RDKit✔️✔️:         nextChar = stream.peek();
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if the next element in the stream is a \" then we have a quoted \".
    // RDKit✔️✔️:         // Otherwise we're done
    // RDKit✔️✔️:         if (nextChar != '"') {
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           // skip the second \"
    // RDKit✔️✔️:           stream.get();
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       strValue += chr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (nextChar == '\'') {
    // RDKit✔️✔️:     std::getline(stream, strValue, '\'');
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     stream >> strValue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::trim_right(strValue);
    // RDKit✔️✔️:   return strValue;
    // RDKit✔️✔️: }
    let mut value = Vec::new();
    match cursor.peek() {
        Some(b' ') => {}
        Some(b'"') => {
            cursor.get();
            while let Some(byte) = cursor.get() {
                if byte == b'"' {
                    if cursor.peek() != Some(b'"') {
                        break;
                    }
                    cursor.get();
                }
                value.push(byte);
            }
        }
        Some(b'\'') => {
            // `std::getline(stream, value, '\'')` starts on the opening quote,
            // so the source returns an empty value and consumes that quote.
            cursor.get();
        }
        Some(_) => {
            cursor.skip_c_locale_whitespace();
            while let Some(byte) = cursor.peek() {
                if matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c) {
                    break;
                }
                value.push(byte);
                cursor.position += 1;
            }
        }
        None => {}
    }
    while value.last().is_some_and(|byte| byte.is_ascii_whitespace()) {
        value.pop();
    }
    String::from_utf8(value)
        .expect("SGroup cursor removes only ASCII delimiters from valid UTF-8 input")
    // END RDKIT CPP FUNCTION
}

fn parse_cstate(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000CStateLabel
    // RDKit✔️✔️: stream.get();  // discard parentheses
    // RDKit✔️✔️: unsigned int count;
    // RDKit✔️✔️: unsigned int bondMark;
    // RDKit✔️✔️: stream >> count >> bondMark;
    // RDKit✔️✔️:
    // RDKit✔️✔️: std::string type = sgroup.getProp<std::string>("TYPE");
    // RDKit✔️✔️: if ((type != "SUP" && count != 1) || (type == "SUP" && count != 4)) {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:   errout << "Unexpected number of fields for CSTATE field on line " << line;
    // RDKit✔️✔️:   SGroupWarnOrThrow<>(strictParsing, errout.str());
    // RDKit✔️✔️:   sgroup.setIsValid(false);
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: Bond *bond = mol->getUniqueBondWithBookmark(bondMark);
    // RDKit✔️✔️:
    // RDKit✔️✔️: RDGeom::Point3D vector;
    // RDKit✔️✔️: if (type == "SUP") {
    // RDKit✔️✔️:   stream >> vector.x >> vector.y >> vector.z;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: try {
    // RDKit✔️✔️:   sgroup.addCState(bond->getIdx(), vector);
    // RDKit✔️✔️: } catch (const std::exception &e) {
    // RDKit✔️✔️:   SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit✔️✔️:   sgroup.setIsValid(false);
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: stream.get();  // discard final parentheses
    // END RDKIT CPP FUNCTION
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addCState
    // RDKit✔️✔️: void SubstanceGroup::addCState(unsigned int bondIdx,
    // RDKit✔️✔️:                                const RDGeom::Point3D &vector) {
    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
    // RDKit✔️✔️:   PRECONDITION(!d_bonds.empty(), "no bonds");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (getBondType(bondIdx) != SubstanceGroup::BondType::XBOND) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Bond with index " << bondIdx
    // RDKit✔️✔️:            << " is not an XBOND for current SubstanceGroup";
    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_cstates.push_back({bondIdx, vector});
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::getBondType
    // RDKit✔️✔️: SubstanceGroup::BondType SubstanceGroup::getBondType(
    // RDKit✔️✔️:     unsigned int bondIdx) const {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end(),
    // RDKit✔️✔️:       "bond is not part of the SubstanceGroup")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bond = dp_mol->getBondWithIdx(bondIdx);
    // RDKit✔️✔️:   bool begin_atom_in_sgroup =
    // RDKit✔️✔️:       std::find(d_atoms.begin(), d_atoms.end(), bond->getBeginAtomIdx()) !=
    // RDKit✔️✔️:       d_atoms.end();
    // RDKit✔️✔️:   bool end_atom_in_sgroup = std::find(d_atoms.begin(), d_atoms.end(),
    // RDKit✔️✔️:                                       bond->getEndAtomIdx()) != d_atoms.end();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (begin_atom_in_sgroup && end_atom_in_sgroup) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::CBOND;
    // RDKit✔️✔️:   } else if (begin_atom_in_sgroup || end_atom_in_sgroup) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::XBOND;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Neither beginning nor ending atoms of bond " << bond->getIdx()
    // RDKit✔️✔️:            << " is in this SubstanceGroup.";
    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION

    let _opening_parenthesis = cursor.get();
    let count = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        ))
    })?;
    let expected = if group.kind() == &SubstanceGroupKind::Superatom {
        4
    } else {
        1
    };
    if count != expected {
        return Err(SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        )));
    }
    let bookmark = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unexpected number of fields for CSTATE field on line {line_number}"
        ))
    })?;
    let bond = *bonds.get(&bookmark).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup bond index {bookmark} out of range on line {line_number}"
        ))
    })?;
    let mut vector = [0.0; 3];
    if group.kind() == &SubstanceGroupKind::Superatom {
        // Point3D is zero-initialized. Formatted extraction stops after its
        // first failed component, leaving that and all later components zero.
        for component in &mut vector {
            let Some(parsed) = cursor.read_double() else {
                break;
            };
            *component = parsed;
        }
    }
    let (begin, end) = bond_endpoints(bond).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup bond index {bookmark} out of range on line {line_number}"
        ))
    })?;
    let begin_is_member = group.atoms().contains(&begin);
    let end_is_member = group.atoms().contains(&end);
    if !group.bonds().contains(&bond) || begin_is_member == end_is_member {
        return Err(SdfReadError::Parse(format!(
            "Bond with index {} is not an XBOND for current SubstanceGroup",
            bond.index()
        )));
    }
    group.push_cstate(SGroupCState { bond, vector });
    let _closing_parenthesis = cursor.get();
    // Behavior review: the parser retains all XYZ components for SUP groups,
    // uses Point3D's zero vector for every other group kind, resolves the
    // externally unique V3000 bookmark to the canonical BondId, and admits
    // only the group's typed crossing-bond membership. Missing/invalid SUP
    // components reproduce formatted-stream fail-state propagation rather
    // than being parsed independently. A missing count/bookmark is rejected
    // structurally because the source reads an uninitialized destination in
    // that malformed case. Complexity review: this is one bounded scalar
    // scan plus one ordered-map lookup and one linear membership lookup,
    // matching the source helper chain without a second representation,
    // whole-table clone, or repeated graph scan.
    Ok(())
}

fn parse_sap(
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    atoms: &BTreeMap<u32, AtomId>,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SAPLabel
    // RDKit✔️✔️: stream.get();  // discard parentheses
    // RDKit✔️✔️: unsigned int count = 0;
    // RDKit✔️✔️: unsigned int aIdxMark = 0;
    // RDKit✔️✔️: std::string lvIdxStr;  // In V3000 this may be a string
    // RDKit✔️✔️: std::string sapIdStr;
    // RDKit✔️✔️: stream >> count >> aIdxMark >> lvIdxStr >> sapIdStr;
    // RDKit✔️✔️: sapIdStr.pop_back();
    // RDKit✔️✔️: unsigned int aIdx = mol->getAtomWithBookmark(aIdxMark)->getIdx();
    // RDKit✔️✔️: int lvIdx = -1;
    // RDKit✔️✔️:
    // RDKit✔️✔️: boost::to_upper(lvIdxStr);
    // RDKit✔️✔️: if (lvIdxStr == "AIDX") {
    // RDKit✔️✔️:   lvIdx = aIdx;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   unsigned int lvIdxTmp = FileParserUtils::toInt(lvIdxStr);
    // RDKit✔️✔️:   if (lvIdxTmp > 0) {
    // RDKit✔️✔️:     lvIdx = mol->getAtomWithBookmark(lvIdxTmp)->getIdx();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: try {
    // RDKit✔️✔️:   sgroup.addAttachPoint(aIdx, lvIdx, sapIdStr);
    // RDKit✔️✔️: } catch (const std::exception &e) {
    // RDKit✔️✔️:   SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit✔️✔️:   sgroup.setIsValid(false);
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    let _discarded_opening = cursor.get();
    let _count = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!("SGroup SAP count missing on line {line_number}"))
    })?;
    let atom_bookmark = cursor.read_unsigned().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup attach atom index missing on line {line_number}"
        ))
    })?;
    let leaving_text = cursor.read_string_token().ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup leaving atom index missing on line {line_number}"
        ))
    })?;
    let label_with_final_byte = cursor.read_string_token().ok_or_else(|| {
        SdfReadError::Parse(format!("SGroup SAP label missing on line {line_number}"))
    })?;
    let label_end = label_with_final_byte.len().checked_sub(1).ok_or_else(|| {
        SdfReadError::Parse(format!("SGroup SAP label missing on line {line_number}"))
    })?;
    let label = label_with_final_byte.get(..label_end).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup SAP label has invalid UTF-8 byte truncation on line {line_number}"
        ))
    })?;
    let atom = *atoms.get(&atom_bookmark).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "SGroup attach atom index {atom_bookmark} out of range on line {line_number}"
        ))
    })?;
    let leaving_atom = if leaving_text.eq_ignore_ascii_case("AIDX") {
        Some(atom)
    } else {
        let bookmark = parse_rdkit_int(leaving_text).map_err(|()| {
            SdfReadError::Parse(format!(
                "Cannot convert '{leaving_text}' to int on line {line_number}"
            ))
        })? as u32;
        if bookmark == 0 {
            None
        } else {
            Some(*atoms.get(&bookmark).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup leaving atom index {bookmark} out of range on line {line_number}"
                ))
            })?)
        }
    };
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addAttachPoint
    // RDKit✔️✔️: void SubstanceGroup::addAttachPoint(unsigned int aIdx, int lvIdx,
    // RDKit✔️✔️:                                     const std::string &idStr) {
    // RDKit✔️✔️:   d_saps.push_back({aIdx, lvIdx, idStr});
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SubstanceGroup::addAttachPoint
    group.push_attach_point(SGroupAttachPoint {
        atom,
        leaving_atom,
        label: Some(label.to_owned()),
        order: None,
    });
    // Behavior review: the count is source-read but deliberately not checked;
    // atom bookmarks use formatted unsigned extraction, AIDX is uppercased by
    // the source, and numeric leaving bookmarks use the distinct toInt
    // contract (so leading '+' is no-conversion zero). Both nonzero bookmarks
    // must resolve in the molecule, but addAttachPoint imposes no ATOMS-member
    // precondition. The final byte of the label token is removed regardless of
    // whether it is ')'; invalid UTF-8 byte truncation is a structured Rust
    // boundary error instead of constructing an invalid string.
    // Complexity review: this advances the existing shared cursor once, does
    // two ordered bookmark lookups at most, and appends one typed value. It is
    // linear in the four token lengths with no token vector or parallel model,
    // matching the source asymptotic and allocation shape.
    Ok(())
    // END RDKIT CPP FUNCTION
}

fn parse_label(
    label: &str,
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    parents: &mut BTreeMap<u32, u32>,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    strict_parsing: bool,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel
    // RDKit✔️✔️: if (label == "XBHEAD" || label == "XBCORR") {
    // RDKit✔️✔️:   std::vector<unsigned int> bvect = ParseV3000Array<unsigned int>(
    // RDKit✔️✔️:       lineStream, mol->getNumBonds(), strictParsing);
    // RDKit✔️✔️:   std::transform(bvect.begin(), bvect.end(), bvect.begin(),
    // RDKit✔️✔️:                  [](unsigned int v) -> unsigned int { return v - 1; });
    // RDKit✔️✔️:   sgroup.setProp(label, bvect);
    // RDKit✔️✔️: } else if (label == "ATOMS") {
    // RDKit❗✔️:   for (auto atomIdx : ParseV3000Array<unsigned int>(
    // RDKit❗✔️:            lineStream, mol->getNumAtoms(), strictParsing)) {
    // RDKit❗✔️:     sgroup.addAtomWithBookmark(atomIdx);
    // RDKit❗✔️:   }
    // RDKit❗✔️: } else if (label == "PATOMS") {
    // RDKit❗✔️:   for (auto patomIdx : ParseV3000Array<unsigned int>(
    // RDKit❗✔️:            lineStream, mol->getNumAtoms(), strictParsing)) {
    // RDKit❗✔️:     sgroup.addParentAtomWithBookmark(patomIdx);
    // RDKit❗✔️:   }
    // RDKit❗✔️: } else if (label == "CBONDS" || label == "XBONDS") {
    // RDKit❗✔️:   for (auto bondIdx : ParseV3000Array<unsigned int>(
    // RDKit❗✔️:            lineStream, mol->getNumBonds(), strictParsing)) {
    // RDKit❗✔️:     sgroup.addBondWithBookmark(bondIdx);
    // RDKit❗✔️:   }
    // RDKit✔️✔️: } else if (label == "BRKXYZ") {
    // RDKit✔️✔️:   auto coords = ParseV3000Array<double>(lineStream, 9, strictParsing);
    // RDKit✔️✔️:   if (coords.size() != 9) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Unexpected number of coordinates for BRKXYZ on line "
    // RDKit✔️✔️:            << line;
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   SubstanceGroup::Bracket bracket;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; ++i) {
    // RDKit✔️✔️:     bracket[i] = RDGeom::Point3D(*(coords.begin() + (3 * i)),
    // RDKit✔️✔️:                                      *(coords.begin() + (3 * i) + 1),
    // RDKit✔️✔️:                                      *(coords.begin() + (3 * i) + 2));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sgroup.addBracket(bracket);
    // RDKit✔️✔️: } else if (label == "CSTATE") {
    // RDKit❗✔️:   ParseV3000CStateLabel(mol, sgroup, lineStream, line, strictParsing);
    // RDKit❗✔️: } else if (label == "SAP") {
    // RDKit❗✔️:   ParseV3000SAPLabel(mol, sgroup, lineStream, strictParsing);
    // RDKit❗✔️: } else if (label == "PARENT") {
    // RDKit❗✔️:   lineStream >> parentIdx;
    // RDKit❗✔️:   sgroup.setProp<unsigned int>("PARENT", parentIdx);
    // RDKit❗✔️: } else if (label == "COMPNO") {
    // RDKit❗✔️:   lineStream >> compno;
    // RDKit❗✔️:   if (compno > 256u) {
    // RDKit❗✔️:     std::ostringstream errout;
    // RDKit❗✔️:     errout << "SGroup SNC value over 256: '" << compno << "' on line "
    // RDKit❗✔️:            << line;
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   sgroup.setProp<unsigned int>("COMPNO", compno);
    // RDKit❗✔️: } else if (label == "FIELDDATA") {
    // RDKit❗✔️:   auto strValue = ParseV3000StringPropLabel(lineStream);
    // RDKit❗✔️:   if (strictParsing) {
    // RDKit❗✔️:     strValue = strValue.substr(0, 200);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   dataFields.push_back(strValue);
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   auto strValue = ParseV3000StringPropLabel(lineStream);
    // RDKit❗✔️:   sgroup.setProp(label, strValue);
    // RDKit❗✔️: }
    match label {
        "XBHEAD" | "XBCORR" => {
            // Unlike XBONDS/CBONDS, these values are one-based bond-table
            // positions, not V3000 bond bookmarks. The canonical model stores
            // the corresponding zero-based BondId so remapping and writers use
            // the same typed reference. RDKit's unsigned `v - 1` can retain an
            // invalid wrapped/out-of-range property value; COSMolKit rejects
            // that value at the detached-model boundary because persistent
            // typed references must satisfy local topology invariants.
            for source_position in
                parse_u32_array(cursor, line_number, Some(bonds.len()), strict_parsing)?
            {
                let index = source_position.checked_sub(1).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup {label} bond-row position 0 is invalid on line {line_number}"
                    ))
                })? as usize;
                if index >= bonds.len() {
                    return Err(SdfReadError::Parse(format!(
                        "SGroup {label} bond-row position {source_position} out of range on line {line_number}"
                    )));
                }
                let bond = BondId::new(index);
                if label == "XBHEAD" {
                    group.push_head_crossing_bond(bond);
                } else {
                    group.push_crossing_bond_correspondence(bond);
                }
            }
        }
        "ATOMS" | "PATOMS" => {
            for bookmark in parse_u32_array(cursor, line_number, Some(atoms.len()), strict_parsing)?
            {
                let atom = *atoms.get(&bookmark).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup atom index {bookmark} out of range on line {line_number}"
                    ))
                })?;
                if label == "ATOMS" {
                    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addAtomWithBookmark
                    // RDKit✔️✔️: void SubstanceGroup::addAtomWithBookmark(int mark) {
                    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
                    // RDKit✔️✔️:   Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
                    // RDKit✔️✔️:   PRECONDITION(atom, "atom not found");
                    // RDKit✔️✔️:   d_atoms.push_back(atom->getIdx());
                    // RDKit✔️✔️: }
                    group.push_atom(atom);
                    // END RDKIT CPP FUNCTION
                } else {
                    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addParentAtomWithBookmark
                    // RDKit✔️✔️: void SubstanceGroup::addParentAtomWithBookmark(int mark) {
                    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
                    // RDKit✔️✔️:
                    // RDKit✔️✔️:   Atom *atom = dp_mol->getUniqueAtomWithBookmark(mark);
                    // RDKit✔️✔️:   unsigned int idx = atom->getIdx();
                    // RDKit✔️✔️:   if (std::find(d_atoms.begin(), d_atoms.end(), idx) == d_atoms.end()) {
                    // RDKit✔️✔️:     std::ostringstream errout;
                    // RDKit✔️✔️:     errout << "Atom with bookmark " << mark
                    // RDKit✔️✔️:            << " is not a member of current SubstanceGroup ";
                    // RDKit✔️✔️:     throw SubstanceGroupException(errout.str());
                    // RDKit✔️✔️:   }
                    // RDKit✔️✔️:
                    // RDKit✔️✔️:   d_patoms.push_back(idx);
                    // RDKit✔️✔️: }
                    if !group.atoms().contains(&atom) {
                        return Err(SdfReadError::Parse(format!(
                            "Atom with bookmark {bookmark} is not a member of current SubstanceGroup on line {line_number}"
                        )));
                    }
                    group.push_parent_atom(atom);
                    // END RDKIT CPP FUNCTION
                }
            }
        }
        "CBONDS" | "XBONDS" => {
            let role = if label == "CBONDS" {
                SGroupBondRole::Contained
            } else {
                SGroupBondRole::Crossing
            };
            for bookmark in parse_u32_array(cursor, line_number, Some(bonds.len()), strict_parsing)?
            {
                let bond = *bonds.get(&bookmark).ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "SGroup bond index {bookmark} out of range on line {line_number}"
                    ))
                })?;
                // BEGIN RDKIT CPP FUNCTION SubstanceGroup::addBondWithBookmark
                // RDKit✔️✔️: void SubstanceGroup::addBondWithBookmark(int mark) {
                // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
                // RDKit✔️✔️:   Bond *bond = dp_mol->getUniqueBondWithBookmark(mark);
                // RDKit✔️✔️:   d_bonds.push_back(bond->getIdx());
                // RDKit✔️✔️: }
                group.push_bond_with_role(bond, role);
                // END RDKIT CPP FUNCTION
            }
        }
        "BRKXYZ" => {
            let coordinates = parse_f64_array(cursor, line_number, Some(9), strict_parsing)?;
            if coordinates.len() != 9 {
                return Err(SdfReadError::Parse(format!(
                    "Unexpected number of coordinates for BRKXYZ on line {line_number}"
                )));
            }
            group.display_mut().brackets.push(SGroupBracket {
                points: [
                    [coordinates[0], coordinates[1], coordinates[2]],
                    [coordinates[3], coordinates[4], coordinates[5]],
                    [coordinates[6], coordinates[7], coordinates[8]],
                ],
            });
            // The fixed nine-value array is copied once into the canonical
            // three-point display value. This remains O(1), allocates only the
            // source-equivalent appended bracket slot, and never enters the
            // molecule coordinate/conformer block.
        }
        _ => {
            parse_non_array_label(
                label,
                cursor,
                line_number,
                group,
                parents,
                atoms,
                bonds,
                bond_endpoints,
                strict_parsing,
            )?;
        }
    }
    Ok(())
    // END RDKIT CPP FUNCTION
}

#[allow(clippy::too_many_arguments)]
fn parse_non_array_label(
    label: &str,
    cursor: &mut SGroupLineCursor<'_>,
    line_number: usize,
    group: &mut SubstanceGroup,
    parents: &mut BTreeMap<u32, u32>,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    strict_parsing: bool,
) -> Result<(), SdfReadError> {
    match label {
        "CSTATE" => {
            // The source passes the shared lineStream into the helper. In
            // particular, a formatted-extraction failbit inside CSTATE stops
            // the outer label loop; isolating the value in a second cursor
            // would incorrectly allow later labels to be consumed.
            parse_cstate(cursor, line_number, group, bonds, bond_endpoints)?;
        }
        "SAP" => {
            parse_sap(cursor, line_number, group, atoms)?;
        }
        "PARENT" => {
            // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (PARENT)
            // RDKit✔️✔️: } else if (label == "PARENT") {
            // RDKit✔️✔️:   // Store relationship until all SGroups have been read
            // RDKit✔️✔️:   unsigned int parentIdx;
            // RDKit✔️✔️:   if (lineStream.eof()) {
            // RDKit✔️✔️:     std::ostringstream errout;
            // RDKit✔️✔️:     errout << "PARENT label not found on line " << line;
            // RDKit✔️✔️:     throw FileParseException(errout.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   lineStream >> parentIdx;
            // RDKit✔️✔️:   if (lineStream.fail()) {
            // RDKit✔️✔️:     std::ostringstream errout;
            // RDKit✔️✔️:     errout << "Invalid PARENT label found on line " << line;
            // RDKit✔️✔️:     throw FileParseException(errout.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   sgroup.setProp<unsigned int>("PARENT", parentIdx);
            if cursor.peek().is_none() {
                return Err(SdfReadError::Parse(format!(
                    "PARENT label not found on line {line_number}"
                )));
            }
            let parent = cursor.read_unsigned().ok_or_else(|| {
                SdfReadError::Parse(format!("Invalid PARENT label found on line {line_number}"))
            })?;
            if cursor.failed {
                return Err(SdfReadError::Parse(format!(
                    "Invalid PARENT label found on line {line_number}"
                )));
            }
            if let Some(sequence) = group.rdkit_sequence_id() {
                parents.insert(sequence, parent);
            }
            // Behavioral review: the shared cursor performs the source's
            // formatted unsigned extraction, including C-locale whitespace,
            // both signs, prefix consumption and fail state. The temporary
            // sequence relation is resolved only after all surviving rows are
            // known, so no raw or dangling PARENT property enters the typed
            // canonical model. Undefined/uninitialized overflow state remains
            // a structured error, matching the documented identity boundary.
            // Complexity review: one bounded decimal scan plus one ordered-map
            // insertion is O(token length + log S), with no duplicate model or
            // additional line scan; this is comparable to formatted extraction
            // followed by source property insertion.
            // END RDKIT CPP FUNCTION
        }
        "COMPNO" => {
            // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (COMPNO)
            // RDKit✔️✔️: } else if (label == "COMPNO") {
            // RDKit✔️✔️:   unsigned int compno;
            // RDKit✔️✔️:   lineStream >> compno;
            // RDKit✔️✔️:   if (compno > 256u) {
            // RDKit✔️✔️:     std::ostringstream errout;
            // RDKit✔️✔️:     errout << "SGroup SNC value over 256: '" << compno << "' on line "
            // RDKit✔️✔️:            << line;
            // RDKit✔️✔️:     throw FileParseException(errout.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   sgroup.setProp<unsigned int>("COMPNO", compno);
            let number = cursor.read_unsigned().ok_or_else(|| {
                SdfReadError::Parse(format!("Invalid COMPNO label found on line {line_number}"))
            })?;
            if cursor.failed {
                return Err(SdfReadError::Parse(format!(
                    "Invalid COMPNO label found on line {line_number}"
                )));
            }
            if number > 256 {
                return Err(SdfReadError::Parse(format!(
                    "SGroup SNC value over 256: '{number}' on line {line_number}"
                )));
            }
            group.set_component_number(number);
            // Formatted unsigned extraction advances the shared label stream;
            // the bound check and one typed assignment are constant time and
            // allocate no intermediate token.
            // END RDKIT CPP FUNCTION
        }
        "FIELDDATA" => {
            // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (DAT fields)
            // RDKit❗✔️: } else if (label == "FIELDDATA") {
            // RDKit❗✔️:   auto strValue = ParseV3000StringPropLabel(lineStream);
            // RDKit❗✔️:   if (strictParsing) {
            // RDKit❗✔️:     strValue = strValue.substr(0, 200);
            // RDKit❗✔️:   }
            // RDKit❗✔️:   dataFields.push_back(strValue);
            let parsed = parse_string_property(cursor);
            let parsed = if strict_parsing && parsed.len() > 200 {
                parsed
                    .get(..200)
                    .ok_or_else(|| {
                        SdfReadError::Parse(format!(
                            "FIELDDATA 200-byte truncation splits UTF-8 on line {line_number}"
                        ))
                    })?
                    .to_owned()
            } else {
                parsed
            };
            group.data_mut().values.push(parsed);
            // The canonical DAT value vector preserves append order. Strict
            // mode performs the source's 200-byte prefix operation; an invalid
            // UTF-8 cut is a structured Rust text-boundary error, so behavior
            // is intentionally partial at that representation boundary.
            // Non-strict mode moves the complete parsed string without a copy.
            // END RDKIT CPP FUNCTION
        }
        "SUBTYPE" => {
            let parsed = parse_string_property(cursor);
            if !is_valid_rdkit_sgroup_subtype(&parsed) {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported SGroup subtype '{parsed}' on line {line_number}"
                )));
            }
            group.set_subtype(parsed);
        }
        "CONNECT" => {
            let parsed = parse_string_property(cursor);
            let connection = sgroup_connection_from_rdkit(&parsed).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "Unsupported SGroup connection type '{parsed}' on line {line_number}"
                ))
            })?;
            group.set_connection(connection);
        }
        "CLASS" => {
            let parsed = parse_string_property(cursor);
            if !is_valid_rdkit_sgroup_class(&parsed) {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported SGroup template class '{parsed}' on line {line_number}"
                )));
            }
            group.set_class(parsed);
        }
        "LABEL" => group.set_label(parse_string_property(cursor)),
        // BEGIN RDKIT CPP FUNCTION ParseV3000ParseLabel (string-property tail)
        // RDKit✔️✔️: } else {
        // RDKit✔️✔️:   // Parse string props
        // RDKit✔️✔️:   auto strValue = ParseV3000StringPropLabel(lineStream);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (label == "SUBTYPE" &&
        // RDKit✔️✔️:       !SubstanceGroupChecks::isValidSubType(strValue)) {
        // RDKit✔️✔️:     std::ostringstream errout;
        // RDKit✔️✔️:     errout << "Unsupported SGroup subtype '" << strValue << "' on line "
        // RDKit✔️✔️:            << line;
        // RDKit✔️✔️:     throw FileParseException(errout.str());
        // RDKit✔️✔️:   } else if (label == "CONNECT" &&
        // RDKit✔️✔️:              !SubstanceGroupChecks::isValidConnectType(strValue)) {
        // RDKit✔️✔️:     std::ostringstream errout;
        // RDKit✔️✔️:     errout << "Unsupported SGroup connection type '" << strValue
        // RDKit✔️✔️:            << "' on line " << line;
        // RDKit✔️✔️:     throw FileParseException(errout.str());
        // RDKit✔️✔️:   } else if (label == "CLASS" &&
        // RDKit✔️✔️:              !SubstanceGroupChecks::isValidClass(strValue)) {
        // RDKit✔️✔️:     std::ostringstream errout;
        // RDKit✔️✔️:     errout << "Unsupported SGroup template class '" << strValue
        // RDKit✔️✔️:            << "' on line " << line;
        // RDKit✔️✔️:     throw FileParseException(errout.str());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   // NATREPLACE is not validated nor used
        // RDKit✔️✔️:
        // RDKit✔️✔️:   sgroup.setProp(label, strValue);
        // RDKit✔️✔️: }
        "FIELDNAME" => group.data_mut().field_name = Some(parse_string_property(cursor)),
        "FIELDTYPE" => group.data_mut().field_type = Some(parse_string_property(cursor)),
        "FIELDINFO" => group.data_mut().field_info = Some(parse_string_property(cursor)),
        "FIELDDISP" => group.data_mut().field_display = Some(parse_string_property(cursor)),
        "QUERYTYPE" => group.data_mut().query_type = Some(parse_string_property(cursor)),
        "QUERYOP" => group.data_mut().query_op = Some(parse_string_property(cursor)),
        // The six DAT metadata labels use the one shared source string parser
        // and overwrite their single canonical typed slots in O(value length),
        // matching the source property assignment without a parallel raw prop.
        // END RDKIT CPP FUNCTION
        "ESTATE" => group.set_expansion_state(parse_string_property(cursor)),
        "BRKTYP" => {
            let parsed = parse_string_property(cursor);
            let style = match parsed.as_str() {
                "BRACKET" => SGroupBracketStyle::Bracket,
                "PAREN" => SGroupBracketStyle::Parenthesis,
                "" => SGroupBracketStyle::None,
                other => SGroupBracketStyle::Unknown(other.to_owned()),
            };
            group.set_bracket_style(style);
        }
        other => group.set_prop(other, parse_string_property(cursor)),
    }
    Ok(())
}

fn skip_overridden_default_value(
    cursor: &mut SGroupLineCursor<'_>,
    label: &str,
    line_number: usize,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (overridden default)
    // RDKit✔️✔️: spacer = lineStream.peek();
    // RDKit✔️✔️: if (spacer == ' ') {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:   errout << "Found unexpected whitespace at DEFAULT label " << label;
    // RDKit✔️✔️:   if (strictParsing) {
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit✔️✔️:     sgroup.setIsValid(false);
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (spacer == '(') {
    // RDKit✔️✔️:   std::getline(lineStream, label, ')');
    // RDKit✔️✔️:   lineStream.get(spacer);
    // RDKit✔️✔️: } else if (spacer == '"') {
    // RDKit✔️✔️:   lineStream.get(spacer);
    // RDKit✔️✔️:   std::getline(lineStream, label, '"');
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   std::getline(lineStream, label, ' ');
    // RDKit✔️✔️:   lineStream.putback(' ');
    // RDKit✔️✔️: }
    match cursor.peek() {
        Some(b' ') => {
            return Err(SdfReadError::Parse(format!(
                "Found unexpected whitespace at DEFAULT label {label} on line {line_number}"
            )));
        }
        Some(b'(') => {
            while cursor.get().is_some_and(|byte| byte != b')') {}
            let _ = cursor.get();
        }
        Some(b'"') => {
            let _ = cursor.get();
            while cursor.get().is_some_and(|byte| byte != b'"') {}
        }
        Some(_) => {
            while cursor.peek().is_some_and(|byte| byte != b' ') {
                cursor.position += 1;
            }
        }
        None => {}
    }
    // Behavioral review: this preserves the source's distinct skip grammar
    // instead of parsing the overridden value. A parenthesized value consumes
    // one byte after `)`, a quoted value stops at the first quote without
    // doubled-quote handling, and a scalar leaves its separating space for
    // the outer label loop. A literal-space empty value invalidates the group.
    // Complexity review: each skipped value is scanned once in place with no
    // allocation, matching the source's linear getline/stream operations.
    Ok(())
    // END RDKIT CPP FUNCTION
}

fn apply_labels(
    text: &str,
    line_number: usize,
    group: &mut SubstanceGroup,
    parents: &mut BTreeMap<u32, u32>,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    seen: &mut BTreeSet<String>,
    defaults_only: bool,
    strict_parsing: bool,
) -> Result<bool, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (label loop)
    // RDKit❗✔️: while (sgroup.getIsValid() && !lineStream.eof() && !lineStream.fail()) {
    // RDKit❗✔️:   lineStream.get(spacer);
    // RDKit❗✔️:   if (lineStream.gcount() == 0) {
    // RDKit❗✔️:     continue;
    // RDKit❗✔️:   } else if (spacer != ' ') {
    // RDKit❗✔️:     if (strictParsing) {
    // RDKit❗✔️:       throw FileParseException(errout.str());
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:       sgroup.setIsValid(false);
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
    // RDKit❗✔️:                        nSgroups, mol, strictParsing);
    // RDKit❗✔️: }
    // RDKit❗✔️: } catch (const std::exception &e) {
    // RDKit❗✔️:   SGroupWarnOrThrow<>(strictParsing, e.what());
    // RDKit❗✔️:   sgroup.setIsValid(false);
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    let mut cursor = SGroupLineCursor::new(text);
    while cursor.position < cursor.text.len() && !cursor.failed {
        let Some(spacer) = cursor.get() else {
            continue;
        };
        if spacer != b' ' {
            let error = SdfReadError::Parse(format!(
                "Found character when expecting a separator (space) on line {line_number}"
            ));
            return if strict_parsing {
                Err(error)
            } else {
                Ok(false)
            };
        }
        let Some(label) = cursor.read_label() else {
            continue;
        };
        if label.is_empty() {
            continue;
        }
        if defaults_only && seen.contains(label) {
            if let Err(error) = skip_overridden_default_value(&mut cursor, label, line_number) {
                return if strict_parsing {
                    Err(error)
                } else {
                    Ok(false)
                };
            }
        } else {
            if let Err(error) = parse_label(
                label,
                &mut cursor,
                line_number,
                group,
                parents,
                atoms,
                bonds,
                bond_endpoints,
                strict_parsing,
            ) {
                if strict_parsing {
                    return Err(error);
                }
                if matches!(&error, SdfReadError::Parse(message) if message == "invalid count value")
                {
                    seen.insert(label.to_owned());
                    continue;
                }
                return Ok(false);
            }
        }
        // The source's `parsedLabels` records only the row's explicit labels;
        // defaults do not add to it, so repeated non-overridden defaults are
        // all parsed in encounter order and later assignments may overwrite.
        if !defaults_only {
            seen.insert(label.to_owned());
        }
    }
    Ok(true)
    // END RDKIT CPP FUNCTION
}

pub(super) fn parse_v3000_sgroup_block(
    lines: &[&str],
    cursor: &mut usize,
    expected_count: usize,
    atoms: &BTreeMap<u32, AtomId>,
    bonds: &BTreeMap<u32, BondId>,
    bond_endpoints: &dyn Fn(BondId) -> Option<(AtomId, AtomId)>,
    strict_parsing: bool,
) -> Result<Vec<SubstanceGroup>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock
    // RDKit❗✔️: // SGroups may be written in unsorted ID order, according to spec, so we will
    // RDKit❗✔️: // temporarily store them in a map before adding them to the mol
    // RDKit❗✔️: IDX_TO_SGROUP_MAP sGroupMap;
    // RDKit❗✔️: auto tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️: if (tempStr.substr(0, 7) == "DEFAULT" && tempStr.length() > 8) {
    // RDKit❗✔️:   defaultString = tempStr.substr(7);
    // RDKit❗✔️:   tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️: }
    // RDKit❗✔️: for (unsigned int si = 0; si < nSgroups; ++si) {
    // RDKit❗✔️:   lineStream >> sequenceId;
    // RDKit❗✔️:   lineStream >> type;
    // RDKit❗✔️:   lineStream >> externalId;
    // RDKit❗✔️:   if (strictParsing && !SubstanceGroupChecks::isValidType(type)) {
    // RDKit❗✔️:     throw MolFileUnhandledFeatureException(errout.str());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   SubstanceGroup sgroup(mol, type);
    // RDKit❗✔️:   sgroup.setProp<unsigned int>("index", sequenceId);
    // RDKit❗✔️:   if (externalId > 0) {
    // RDKit❗✔️:     sgroup.setProp<unsigned int>("ID", externalId);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ParseV3000ParseLabel(label, lineStream, dataFields, line, sgroup,
    // RDKit❗✔️:                        nSgroups, mol, strictParsing);
    // RDKit❗✔️:   // Process defaults
    // RDKit❗✔️:   if (std::find(parsedLabels.begin(), parsedLabels.end(), label) ==
    // RDKit❗✔️:       parsedLabels.end()) {
    // RDKit❗✔️:     ParseV3000ParseLabel(label, lineStream, dataFields, defaultLineNum,
    // RDKit❗✔️:                          sgroup, nSgroups, mol, strictParsing);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   sGroupMap.emplace(sequenceId, sgroup);
    // RDKit❗✔️:   tempStr = FileParserUtils::getV3000Line(inStream, line);
    // RDKit❗✔️: }
    // RDKit❗✔️: if (sGroupMap.size() != nSgroups) {
    // RDKit❗✔️:   std::ostringstream errout;
    // RDKit❗✔️:   errout << "Found " << sGroupMap.size() << " SGroups when " << nSgroups
    // RDKit❗✔️:          << " were expected." << std::endl;
    // RDKit❗✔️:   if (strictParsing) {
    // RDKit❗✔️:     throw FileParseException(errout.str());
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: for (const auto &sg : sGroupMap) {
    // RDKit❗✔️:   if (sg.second.getIsValid()) {
    // RDKit❗✔️:     addSubstanceGroup(*mol, sg.second);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit❗✔️:                             << " is invalid and will be ignored" << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let (mut current, mut line_number) = get_v3000_line(lines, cursor)?;
    let mut defaults = String::new();
    let mut default_line = line_number;
    if current.starts_with("DEFAULT") && current.len() > 8 {
        defaults = current[7..].trim_end().to_owned();
        default_line = line_number;
        (current, line_number) = get_v3000_line(lines, cursor)?;
    }

    let mut groups = BTreeMap::<u32, SubstanceGroup>::new();
    let mut invalid_sequences = BTreeSet::new();
    let mut parents = BTreeMap::<u32, u32>::new();
    for _ in 0..expected_count {
        if current.to_ascii_uppercase().starts_with("END SGROUP") {
            break;
        }
        let (sequence_text, kind_text, external_text, labels) =
            split_sgroup_line(current.trim_end(), line_number)?;
        let sequence = parse_sgroup_header_unsigned(sequence_text, "sequence ID", line_number)?;
        if strict_parsing && !is_valid_rdkit_sgroup_type(kind_text) {
            return Err(SdfReadError::Parse(format!(
                "Unsupported SGroup type '{kind_text}' on line {line_number}"
            )));
        }
        let external_id = parse_sgroup_header_unsigned(external_text, "external ID", line_number)?;
        let mut group = SubstanceGroup::new(
            SubstanceGroupId::new(groups.len()),
            sgroup_kind_from_rdkit_type(kind_text),
        );
        group.set_rdkit_sequence_id(sequence);
        group.set_prop("TYPE", kind_text);
        if external_id != 0 {
            group.set_external_id(external_id);
        }
        let mut seen = BTreeSet::new();
        let mut candidate_parents = BTreeMap::new();
        let valid = apply_labels(
            &labels,
            line_number,
            &mut group,
            &mut candidate_parents,
            atoms,
            bonds,
            bond_endpoints,
            &mut seen,
            false,
            strict_parsing,
        )?;
        let valid = if valid {
            apply_labels(
                &defaults,
                default_line,
                &mut group,
                &mut candidate_parents,
                atoms,
                bonds,
                bond_endpoints,
                &mut seen,
                true,
                strict_parsing,
            )?
        } else {
            false
        };
        if let std::collections::btree_map::Entry::Vacant(entry) = groups.entry(sequence) {
            entry.insert(group);
            parents.extend(candidate_parents);
            if !valid {
                invalid_sequences.insert(sequence);
            }
        }
        (current, line_number) = get_v3000_line(lines, cursor)?;
    }
    if !current.to_ascii_uppercase().starts_with("END SGROUP") {
        if strict_parsing {
            return Err(SdfReadError::Parse(format!(
                "END SGROUP line not found on line {line_number}"
            )));
        }
        *cursor = cursor.saturating_sub(1);
    }
    if groups.len() != expected_count && strict_parsing {
        return Err(SdfReadError::Parse(format!(
            "Found {} SGroups when {expected_count} were expected.",
            groups.len()
        )));
    }

    groups.retain(|sequence, _| !invalid_sequences.contains(sequence));

    // BEGIN RDKIT CPP FUNCTION ParseV3000SGroupsBlock (installation)
    // RDKit❗✔️: // SGroups successfully parsed, now add them to the molecule
    // RDKit❗✔️: for (const auto &sg : sGroupMap) {
    // RDKit❗✔️:   if (sg.second.getIsValid()) {
    // RDKit❗✔️:     addSubstanceGroup(*mol, sg.second);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "SGroup " << sg.first
    // RDKit❗✔️:                             << " is invalid and will be ignored" << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit stores PARENT as an unchecked unsigned property. The canonical
    // COSMolKit model instead requires a typed in-range SubstanceGroupId. In
    // strict mode an unresolved source sequence is therefore a structured
    // model-boundary error. In non-strict mode the complete child is removed,
    // recursively, so removing an invalid parent can never leave a dangling
    // or partially retained hierarchy. Self-parent and multi-row cycles remain
    // representable because all referenced source rows survive and neither the
    // pinned source nor canonical local validation forbids those relations.
    let children_by_parent = parents.iter().fold(
        BTreeMap::<u32, Vec<u32>>::new(),
        |mut children, (sequence, parent_sequence)| {
            children
                .entry(*parent_sequence)
                .or_default()
                .push(*sequence);
            children
        },
    );
    let mut unresolved = parents
        .iter()
        .filter_map(|(sequence, parent_sequence)| {
            (groups.contains_key(sequence) && !groups.contains_key(parent_sequence))
                .then_some((*sequence, *parent_sequence))
        })
        .collect::<Vec<_>>();
    if strict_parsing && let Some((sequence, parent_sequence)) = unresolved.first().copied() {
        return Err(SdfReadError::Parse(format!(
            "SGroup {sequence} references missing parent SGroup {parent_sequence}"
        )));
    }
    while let Some((sequence, _)) = unresolved.pop() {
        if groups.remove(&sequence).is_none() {
            continue;
        }
        if let Some(children) = children_by_parent.get(&sequence) {
            unresolved.extend(
                children
                    .iter()
                    .filter_map(|child| groups.contains_key(child).then_some((*child, sequence))),
            );
        }
    }

    let ids = groups
        .keys()
        .enumerate()
        .map(|(index, sequence)| (*sequence, SubstanceGroupId::new(index)))
        .collect::<BTreeMap<_, _>>();
    for (sequence, id) in &ids {
        groups
            .get_mut(sequence)
            .expect("key came from group map")
            .set_id(*id);
    }
    for (sequence, parent_sequence) in parents {
        let Some(group) = groups.get_mut(&sequence) else {
            continue;
        };
        let parent = ids
            .get(&parent_sequence)
            .expect("unresolved parent groups were removed before compact ID assignment");
        group.set_parent(*parent);
    }
    // Behavioral review: sorted surviving source rows receive compact IDs once
    // and every retained PARENT relation is installed against that same map;
    // strict failure returns no partial record and non-strict cascading removal
    // retains no child whose parent was discarded. This is an explicit typed
    // model boundary beyond RDKit's unchecked property storage.
    // Complexity review: the reverse parent index and removal work are linear
    // apart from ordered-map/set lookups, O(S log S), with each row removed at
    // most once. This preserves the source installation scale without repeated
    // hierarchy scans or cloning group state.
    // END RDKIT CPP FUNCTION
    Ok(groups.into_values().collect())
    // END RDKIT CPP FUNCTION
}

fn parse_stereo_collection_line(
    line: &str,
    line_number: usize,
    atom_count: usize,
) -> Result<Option<StereoGroup>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION parseEnhancedStereo (line recognition)
    // RDKit✔️✔️: const regex stereo_label(
    // RDKit✔️✔️:     R"regex(MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *)regex");
    // Recognition and the matched tag/index branches are reproduced below;
    // the caller owns duplicate-ABS strictness and collection installation.
    // Marker review: behavior verified by the closed collection regressions
    // (skip/recognition, id rules, row positions, strictness); the manual
    // single-pass scan replaces the source's per-call `std::regex`
    // construction and match without changing any accept/reject outcome.
    // `regex_match` requires the whole line to match: the tag is exactly three
    // characters, the optional group id is digits only, one or more spaces
    // separate it from `ATOMS=(`, the count is `[0-9]+` followed by spaces,
    // and only spaces may follow the closing parenthesis. The block caller
    // reproduces the source's asymmetric normalization: its first logical
    // payload is uppercased, while subsequent payloads reach this matcher
    // with their original case.
    // Non-matching lines are unrecognized collection types and are skipped,
    // not parsed and not errors.
    let Some(rest) = line.strip_prefix("MDLV30/STE") else {
        return Ok(None);
    };
    if rest.len() < 3 {
        return Ok(None);
    }
    // The source regex consumes exactly three bytes before ASCII digits/spaces.
    // A split inside UTF-8 cannot satisfy that suffix; it is an unrecognized
    // collection line, not a Rust string slicing panic.
    let Some(tag) = rest.get(..3) else {
        return Ok(None);
    };
    let after_tag = &rest[3..];
    let id_digits = after_tag.bytes().take_while(u8::is_ascii_digit).count();
    let after_digits = &after_tag[id_digits..];
    let after_id_separator = after_digits.trim_start_matches(' ');
    if after_id_separator.len() == after_digits.len() {
        return Ok(None);
    }
    let Some(after_atoms) = after_id_separator.strip_prefix("ATOMS=(") else {
        return Ok(None);
    };
    let count_digits = after_atoms.bytes().take_while(u8::is_ascii_digit).count();
    if count_digits == 0 {
        return Ok(None);
    }
    let count_text = &after_atoms[..count_digits];
    let after_count = &after_atoms[count_digits..];
    let atoms_region = after_count.trim_start_matches(' ');
    if atoms_region.len() == after_count.len() {
        return Ok(None);
    }
    // `(.*)\)` is greedy: the group ends at the last ')' that is followed
    // only by the regex's trailing spaces.
    let Some(close) = atoms_region.rfind(')') else {
        return Ok(None);
    };
    if !atoms_region[close + 1..].bytes().all(|byte| byte == b' ') {
        return Ok(None);
    }
    let atoms_text = &atoms_region[..close];
    // BEGIN RDKIT CPP FUNCTION parseEnhancedStereo (tag dispatch)
    // RDKit✔️✔️:       if (match[1] == "ABS") {
    // RDKit✔️✔️:         grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️✔️:       } else if (match[1] == "REL") {
    // RDKit✔️✔️:         grouptype = RDKit::StereoGroupType::STEREO_OR;
    // RDKit✔️✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit✔️✔️:       } else if (match[1] == "RAC") {
    // RDKit✔️✔️:         grouptype = RDKit::StereoGroupType::STEREO_AND;
    // RDKit✔️✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         errout << "Unrecognized stereogroup type : '" << tempStr
    // RDKit✔️✔️:                << "' on line" << line;
    // RDKit✔️✔️:         throw FileParseException(errout.str());
    // RDKit✔️✔️:       }
    let kind = match tag {
        "ABS" => StereoGroupKind::Absolute,
        "REL" => StereoGroupKind::Or,
        "RAC" => StereoGroupKind::And,
        _ => {
            return Err(SdfReadError::Parse(format!(
                "Unrecognized stereogroup type : '{line}' on line{line_number}"
            )));
        }
    };
    // The source constructs every group with its `groupid`: ABS keeps the
    // initialized zero, while REL/RAC run their digits through
    // `FileParserUtils::toUnsigned` — which resolves empty or overflowing
    // digit strings to zero because the `std::from_chars` error code is
    // ignored (the regex already guarantees digits, so no cast throw exists).
    let group_id = if tag == "ABS" {
        0
    } else {
        parse_rdkit_unsigned(&after_tag[..id_digits]).map_err(|()| {
            SdfReadError::Parse(format!(
                "Cannot convert stereo group id on line {line_number}"
            ))
        })?
    };
    // RDKit✔️✔️:       const unsigned int count = FileParserUtils::toUnsigned(match[3], true);
    // RDKit✔️✔️:       std::vector<Atom *> atoms;
    // RDKit✔️✔️:       std::stringstream ss(match[4]);
    // RDKit✔️✔️:       unsigned int index;
    // RDKit✔️✔️:       for (size_t i = 0; i < count; ++i) {
    // RDKit✔️✔️:         ss >> index;
    // RDKit✔️✔️:         // atoms are 1 indexed in molfiles
    // RDKit✔️✔️:         atoms.push_back(mol->getAtomWithIdx(index - 1));
    // RDKit✔️✔️:       }
    // COLLECTION atom values are 1-based main CTAB atom row positions
    // (`getAtomWithIdx(index - 1)`), not V3000 atom bookmarks: a nonsequential
    // bookmark table cannot change collection membership. `ss >> index`
    // is decimal formatted extraction, not whitespace tokenization: it consumes
    // a numeric prefix, so `1+2` supplies two values and `1x` supplies one.
    // Its sentry skips C-locale whitespace; EOF leaves the previous index
    // untouched. No digits produces zero, unsigned negation wraps, and overflow
    // produces UINT_MAX/failbit. Zero/out-of-range indices fail immediately.
    // EOF before the first index would use uninitialized C++ state: reject it
    // structurally rather than fabricate an atom. Later EOF repeats the last
    // initialized index, including when the declared count exceeds the text.
    let count = parse_rdkit_unsigned(count_text).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{count_text}' to unsigned int on line {line_number}"
        ))
    })? as usize;
    let bytes = atoms_text.as_bytes();
    let mut position_cursor = 0;
    let mut previous_position = None;
    let mut extraction_failed = false;
    let group_atoms = (0..count)
        .map(|_| {
            let out_of_range = |value: String| {
                SdfReadError::Parse(format!(
                    "Stereo group atom index {value} out of range on line {line_number}"
                ))
            };
            while bytes
                .get(position_cursor)
                .is_some_and(|byte| matches!(byte, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c))
            {
                position_cursor += 1;
            }
            let position = if extraction_failed || position_cursor == bytes.len() {
                previous_position.ok_or_else(|| {
                    SdfReadError::Parse(format!(
                        "Stereo group has no initialized atom index on line {line_number}"
                    ))
                })?
            } else {
                let negative = bytes[position_cursor] == b'-';
                if matches!(bytes[position_cursor], b'+' | b'-') {
                    position_cursor += 1;
                }
                let digits_start = position_cursor;
                let mut value = Some(0_u32);
                while let Some(byte) = bytes.get(position_cursor).filter(|b| b.is_ascii_digit()) {
                    value =
                        value.and_then(|v| v.checked_mul(10)?.checked_add(u32::from(byte - b'0')));
                    position_cursor += 1;
                }
                extraction_failed = value.is_none() || position_cursor == digits_start;
                let value = match value {
                    None => u32::MAX,
                    Some(value) if negative => value.wrapping_neg(),
                    Some(value) => value,
                };
                previous_position = Some(value);
                value
            };
            let row = usize::try_from(
                position
                    .checked_sub(1)
                    .ok_or_else(|| out_of_range(position.to_string()))?,
            )
            .map_err(|_| out_of_range(position.to_string()))?;
            if row >= atom_count {
                return Err(out_of_range(position.to_string()));
            }
            Ok(AtomId::new(row))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(
        StereoGroup::new(kind, group_atoms, Vec::new()).with_id(group_id),
    ))
}

pub(super) fn parse_v3000_collection_block(
    lines: &[&str],
    cursor: &mut usize,
    atom_count: usize,
    strict_parsing: bool,
) -> Result<Vec<StereoGroup>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION parseEnhancedStereo
    // RDKit❗✔️: std::string parseEnhancedStereo(std::istream *inStream, unsigned int &line,
    // RDKit❗✔️:                                 RWMol *mol, bool strictParsing) {
    // RDKit❗✔️:   // Lines like (absolute, relative, racemic):
    // RDKit❗✔️:   // M  V30 MDLV30/STEABS ATOMS=(2 2 3)
    // RDKit❗✔️:   // M  V30 MDLV30/STEREL1 ATOMS=(1 12)
    // RDKit❗✔️:   // M  V30 MDLV30/STERAC1 ATOMS=(1 12)
    // RDKit❗✔️:   const regex stereo_label(
    // RDKit❗✔️:       R"regex(MDLV30/STE(...)([0-9]*) +ATOMS=\(([0-9]+) +(.*)\) *)regex");
    // RDKit❗✔️:
    // RDKit❗✔️:   smatch match;
    // RDKit❗✔️:   std::vector<StereoGroup> groups;
    // RDKit❗✔️:
    // RDKit❗✔️:   // Read the collection until the end
    // RDKit❗✔️:   auto tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   boost::to_upper(tempStr);
    // RDKit❗✔️:   unsigned abs_group_seen = 0;
    // RDKit❗✔️:   while (!startsWith(tempStr, "END", 3)) {
    // RDKit❗✔️:     // If this line in the collection is part of a stereo group
    // RDKit❗✔️:     if (regex_match(tempStr, match, stereo_label)) {
    // RDKit❗✔️:       StereoGroupType grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit❗✔️:       unsigned groupid = 0;
    // RDKit❗✔️:
    // RDKit❗✔️:       if (match[1] == "ABS") {
    // RDKit❗✔️:         grouptype = RDKit::StereoGroupType::STEREO_ABSOLUTE;
    // RDKit❗✔️:         // Warn only one per mol about multiple ABS groups
    // RDKit❗✔️:         if (abs_group_seen == 1) {
    // RDKit❗✔️:           std::ostringstream errout;
    // RDKit❗✔️:           errout << "Seen a second ABS stereo group on line " << line
    // RDKit❗✔️:                  << std::endl;
    // RDKit❗✔️:           if (strictParsing) {
    // RDKit❗✔️:             throw FileParseException(errout.str());
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             BOOST_LOG(rdWarningLog) << errout.str() << std::endl;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         ++abs_group_seen;
    // RDKit❗✔️:       } else if (match[1] == "REL") {
    // RDKit❗✔️:         grouptype = RDKit::StereoGroupType::STEREO_OR;
    // RDKit❗✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit❗✔️:       } else if (match[1] == "RAC") {
    // RDKit❗✔️:         grouptype = RDKit::StereoGroupType::STEREO_AND;
    // RDKit❗✔️:         groupid = FileParserUtils::toUnsigned(match[2], true);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         std::ostringstream errout;
    // RDKit❗✔️:         errout << "Unrecognized stereogroup type : '" << tempStr << "' on line"
    // RDKit❗✔️:                << line;
    // RDKit❗✔️:         throw FileParseException(errout.str());
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       const unsigned int count = FileParserUtils::toUnsigned(match[3], true);
    // RDKit❗✔️:       std::vector<Atom *> atoms;
    // RDKit❗✔️:       std::stringstream ss(match[4]);
    // RDKit❗✔️:       unsigned int index;
    // RDKit❗✔️:       for (size_t i = 0; i < count; ++i) {
    // RDKit❗✔️:         ss >> index;
    // RDKit❗✔️:         // atoms are 1 indexed in molfiles
    // RDKit❗✔️:         atoms.push_back(mol->getAtomWithIdx(index - 1));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       std::vector<Bond *> newBonds;
    // RDKit❗✔️:       groups.emplace_back(grouptype, std::move(atoms), std::move(newBonds),
    // RDKit❗✔️:                           groupid);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // skip collection types we don't know how to read. Only one documented
    // RDKit❗✔️:       // is MDLV30/HILITE
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "Skipping unrecognized collection type at "
    // RDKit❗✔️:                                  "line "
    // RDKit❗✔️:                               << line << ": " << tempStr << std::endl;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!groups.empty()) {
    // RDKit❗✔️:     mol->setStereoGroups(std::move(groups));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   tempStr = getV3000Line(inStream, line);
    // RDKit❗✔️:   return tempStr;
    // RDKit❗✔️: }
    // The recognition/index extraction is delegated to the private helper;
    // installation is owned by read_v3000_record_detached. The full source
    // remains here; this correction does not certify every parser branch.
    // Behavior review: get_v3000_line preserves case. Exactly the first
    // assembled payload is uppercased before the loop; later payloads and the
    // END-prefix check retain their source case, matching the placement of
    // boost::to_upper above. A new call resets this first-line rule. The
    // remaining partial marker is solely the documented first-extraction
    // undefined-state boundary in parse_stereo_collection_line.
    // Complexity review: each logical line is read and scanned once. Only the
    // first payload allocates an uppercase copy, reducing rather than
    // increasing the source-corresponding linear normalization work; group
    // accumulation remains linear in total recognized membership.
    let mut groups = Vec::new();
    let mut absolute_count = 0_usize;
    let (first_line, mut line_number) = get_v3000_line(lines, cursor)?;
    let mut line = first_line.to_ascii_uppercase();
    while !line.starts_with("END") {
        if let Some(group) = parse_stereo_collection_line(&line, line_number, atom_count)? {
            if group.kind() == StereoGroupKind::Absolute {
                absolute_count += 1;
                if absolute_count > 1 && strict_parsing {
                    return Err(SdfReadError::Parse(format!(
                        "Seen a second ABS stereo group on line {line_number}\n"
                    )));
                }
            }
            groups.push(group);
        }
        (line, line_number) = get_v3000_line(lines, cursor)?;
    }
    Ok(groups)
    // END RDKIT CPP FUNCTION
}

#[derive(Debug)]
pub(super) struct V2000SgroupState {
    strict_parsing: bool,
    groups: BTreeMap<u32, SubstanceGroup>,
    invalid_sequences: BTreeSet<u32>,
    parent_by_sequence: BTreeMap<u32, u32>,
    scd_counter: u32,
    last_data_group: u32,
    current_data_field: String,
    pending_attach_points: Vec<(u32, AtomId, Option<String>)>,
}

impl V2000SgroupState {
    pub(super) fn new(strict_parsing: bool) -> Self {
        Self {
            strict_parsing,
            groups: BTreeMap::new(),
            invalid_sequences: BTreeSet::new(),
            parent_by_sequence: BTreeMap::new(),
            scd_counter: 0,
            last_data_group: 0,
            current_data_field: String::new(),
            pending_attach_points: Vec::new(),
        }
    }
}

impl Default for V2000SgroupState {
    fn default() -> Self {
        Self::new(true)
    }
}

fn parse_v2000_int_field(
    line: &str,
    line_number: usize,
    position: &mut usize,
    is_counter: bool,
) -> Result<u32, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseSGroupIntField
    // RDKit❗✔️: ++pos;  // Account for separation space
    // RDKit❗✔️: size_t len = 3 - isFieldCounter;  // field counters are smaller
    // RDKit❗✔️: fieldValue = FileParserUtils::toInt(text.substr(pos, len));
    // RDKit❗✔️: pos += len;
    // RDKit❗✔️: return fieldValue;
    *position += 1;
    let length = 3 - usize::from(is_counter);
    if *position >= line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let text = rdkit_substr(line, *position, length);
    let value = parse_rdkit_int(text).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{text}' to int on line {line_number}"
        ))
    })? as u32;
    *position += length;
    Ok(value)
    // END RDKIT CPP FUNCTION
}

fn parse_v2000_double_field(
    line: &str,
    line_number: usize,
    position: &mut usize,
) -> Result<f64, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION ParseSGroupDoubleField
    // RDKit❗✔️: size_t len = 10;
    // RDKit❗✔️: fieldValue = FileParserUtils::toDouble(text.substr(pos, len));
    // RDKit❗✔️: pos += len;
    // RDKit❗✔️: return fieldValue;
    let length = 10;
    if *position >= line.len() {
        return Err(SdfReadError::Parse(format!(
            "SGroup line too short: '{line}' on line {line_number}"
        )));
    }
    let text = rdkit_substr(line, *position, length);
    let value = parse_rdkit_double(text).map_err(|()| {
        SdfReadError::Parse(format!(
            "Cannot convert '{text}' to double on line {line_number}"
        ))
    })?;
    *position += length;
    Ok(value)
    // END RDKIT CPP FUNCTION
}

impl V2000SgroupState {
    fn group_mut(&mut self, sequence: u32) -> Option<&mut SubstanceGroup> {
        // BEGIN RDKIT CPP FUNCTION FindSgIdx
        // RDKit❗✔️: auto sgIt = sGroupMap.find(sgIdx);
        // RDKit❗✔️: if (sgIt == sGroupMap.end()) {
        // RDKit❗✔️:   return nullptr;
        // RDKit❗✔️: }
        // RDKit❗✔️: return &sgIt->second;
        self.groups.get_mut(&sequence)
        // END RDKIT CPP FUNCTION
    }

    fn group_mut_if_present(&mut self, sequence: u32) -> Option<&mut SubstanceGroup> {
        self.group_mut(sequence)
    }

    fn recoverable<T>(&self, result: Result<T, SdfReadError>) -> Result<Option<T>, SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupIntField (policy overload)
        // RDKit❗✔️: try {
        // RDKit❗✔️:   res = ParseSGroupIntField(text, line, pos, isFieldCounter);
        // RDKit❗✔️: } catch (const std::exception &e) {
        // RDKit❗✔️:   if (strictParsing) {
        // RDKit❗✔️:     throw;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     ok = false;
        // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << e.what() << std::endl;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        match result {
            Ok(value) => Ok(Some(value)),
            Err(error) if self.strict_parsing => Err(error),
            Err(_) => Ok(None),
        }
        // END RDKIT CPP FUNCTION
    }

    fn invalidate_or_throw(
        &mut self,
        sequence: u32,
        error: SdfReadError,
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION SGroupWarnOrThrow
        // RDKit❗✔️: if (strictParsing) {
        // RDKit❗✔️:   throw Exc(msg);
        // RDKit❗✔️: } else {
        // RDKit❗✔️:   BOOST_LOG(rdWarningLog) << msg << std::endl;
        // RDKit❗✔️: }
        if self.strict_parsing {
            Err(error)
        } else {
            self.invalid_sequences.insert(sequence);
            Ok(())
        }
        // END RDKIT CPP FUNCTION
    }

    pub(super) fn parse_line(
        &mut self,
        line: &str,
        line_number: usize,
        atom_count: usize,
        bond_endpoints: &[(AtomId, AtomId)],
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseMolBlockProperties (SGroup dispatch)
        // RDKit❗✔️: } else if (lineBeg == "M  STY") {
        // RDKit❗✔️:   ParseSGroupV2000STYLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SST") {
        // RDKit❗✔️:   ParseSGroupV2000SSTLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SLB") {
        // RDKit❗✔️:   ParseSGroupV2000SLBLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SCN") {
        // RDKit❗✔️:   ParseSGroupV2000SCNLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SDS") {
        // RDKit❗✔️:   ParseSGroupV2000SDSLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SAL" || lineBeg == "M  SBL" ||
        // RDKit❗✔️:            lineBeg == "M  SPA") {
        // RDKit❗✔️:   ParseSGroupV2000VectorDataLine(sGroupMap, mol, tempStr, line,
        // RDKit❗✔️:                                  strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SMT") {
        // RDKit❗✔️:   ParseSGroupV2000SMTLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SDI") {
        // RDKit❗✔️:   ParseSGroupV2000SDILine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SBV") {
        // RDKit❗✔️:   ParseSGroupV2000SBVLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SDT") {
        // RDKit❗✔️:   ParseSGroupV2000SDTLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SDD") {
        // RDKit❗✔️:   ParseSGroupV2000SDDLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SCD" || lineBeg == "M  SED") {
        // RDKit❗✔️:   ParseSGroupV2000SCDSEDLine(sGroupMap, dataFieldsMap, mol, tempStr, line,
        // RDKit❗✔️:                              strictParsing, SCDcounter, lastDataSGroup,
        // RDKit❗✔️:                              currentDataField);
        // RDKit❗✔️: } else if (lineBeg == "M  SPL") {
        // RDKit❗✔️:   ParseSGroupV2000SPLLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SNC") {
        // RDKit❗✔️:   ParseSGroupV2000SNCLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SAP") {
        // RDKit❗✔️:   ParseSGroupV2000SAPLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SCL") {
        // RDKit❗✔️:   ParseSGroupV2000SCLLine(sGroupMap, mol, tempStr, line, strictParsing);
        // RDKit❗✔️: } else if (lineBeg == "M  SBT") {
        // RDKit❗✔️:   ParseSGroupV2000SBTLine(sGroupMap, mol, tempStr, line, strictParsing);
        match rdkit_substr(line, 0, 6) {
            "M  STY" => self.parse_sty(line, line_number),
            "M  SST" => self.parse_sst(line, line_number),
            "M  SLB" => self.parse_slb(line, line_number),
            "M  SCN" => self.parse_scn(line, line_number),
            "M  SDS" => self.parse_sds(line, line_number),
            "M  SAL" | "M  SBL" | "M  SPA" => {
                self.parse_vector(line, line_number, atom_count, bond_endpoints.len())
            }
            "M  SMT" => self.parse_smt(line, line_number),
            "M  SDI" => self.parse_sdi(line, line_number),
            "M  SBV" => self.parse_sbv(line, line_number, bond_endpoints),
            "M  SDT" => self.parse_sdt(line, line_number),
            "M  SDD" => self.parse_sdd(line, line_number),
            "M  SCD" | "M  SED" => self.parse_scd_sed(line, line_number),
            "M  SPL" => self.parse_spl(line, line_number),
            "M  SNC" => self.parse_snc(line, line_number),
            "M  SAP" => self.parse_sap(line, line_number, atom_count),
            "M  SCL" => self.parse_scl(line, line_number),
            "M  SBT" => self.parse_sbt(line, line_number),
            prefix => Err(SdfReadError::Parse(format!(
                "Unknown SGroup line '{prefix}' on line {line_number}"
            ))),
        }
        // END RDKIT CPP FUNCTION
    }

    fn parse_sty(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000STYLine
        // RDKit❗✔️: unsigned int pos = 6;
        // RDKit❗✔️: unsigned int nent =
        // RDKit❗✔️:     ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️: unsigned int sequenceId =
        // RDKit❗✔️:     ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: std::string typ = text.substr(pos + 1, 3);
        // RDKit❗✔️: if (SubstanceGroupChecks::isValidType(typ)) {
        // RDKit❗✔️:   auto sgroup = SubstanceGroup(mol, typ);
        // RDKit❗✔️:   sgroup.setProp<unsigned int>("index", sequenceId);
        // RDKit❗✔️:   sGroupMap.emplace(sequenceId, sgroup);
        // RDKit❗✔️: }
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup STY line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            let kind_text = rdkit_substr(line, position + 1, 3);
            if !is_valid_rdkit_sgroup_type(kind_text) {
                if self.strict_parsing {
                    return Err(SdfReadError::Parse(format!(
                        "S group {kind_text} on line {line_number}"
                    )));
                }
                position += 4;
                continue;
            }
            let mut group = SubstanceGroup::new(
                SubstanceGroupId::new(self.groups.len()),
                sgroup_kind_from_rdkit_type(kind_text),
            );
            group.set_rdkit_sequence_id(sequence);
            group.set_prop("TYPE", kind_text);
            self.groups.entry(sequence).or_insert(group);
            position += 4;
        }
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_sst(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SST line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let subtype = rdkit_substr(line, position + 1, 3);
            if !matches!(subtype, "ALT" | "RAN" | "BLO") {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "Unsupported SGroup subtype '{subtype}' on line {line_number}"
                    )),
                );
            }
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup subtype")
                .set_subtype(subtype);
            position += 4;
        }
        Ok(())
    }

    fn parse_slb(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SLB line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let Some(id) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup label ID")
                .set_external_id(id);
        }
        Ok(())
    }

    fn parse_scn(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 7 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SCN line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let text = rdkit_substr(line, position + 1, 2);
            let Some(connection) = sgroup_connection_from_rdkit(text) else {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "Unsupported SGroup connection type '{text}' on line {line_number}"
                    )),
                );
            };
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup connection")
                .set_connection(connection);
            position += 3;
        }
        Ok(())
    }

    fn parse_sds(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        if !line.starts_with("M  SDS EXP") {
            return Err(SdfReadError::Parse(format!(
                "bad SDS line on line {line_number}"
            )));
        }
        let mut position = 10;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 4 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SDS line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            let Some(group) = self.group_mut_if_present(sequence) else {
                return Ok(());
            };
            group.set_expansion_state("E");
        }
        Ok(())
    }

    fn parse_vector(
        &mut self,
        line: &str,
        line_number: usize,
        atom_count: usize,
        bond_count: usize,
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000VectorDataLine
        // RDKit❗✔️: if (typ == "SAL") {
        // RDKit❗✔️:   sGroupAddIndexedElement = &SubstanceGroup::addAtomWithBookmark;
        // RDKit❗✔️: } else if (typ == "SBL") {
        // RDKit❗✔️:   sGroupAddIndexedElement = &SubstanceGroup::addBondWithBookmark;
        // RDKit❗✔️: } else if (typ == "SPA") {
        // RDKit❗✔️:   sGroupAddIndexedElement = &SubstanceGroup::addParentAtomWithBookmark;
        // RDKit❗✔️: }
        // RDKit❗✔️: unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: unsigned int nent =
        // RDKit❗✔️:     ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️: if (!ok) {
        // RDKit❗✔️:   sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: if (text.size() < pos + 4) {
        // RDKit❗✔️:   SGroupWarnOrThrow<>(strictParsing, errout.str());
        // RDKit❗✔️:   sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: unsigned int nbr = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: if (!ok) {
        // RDKit❗✔️:   sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: (sgroup->*sGroupAddIndexedElement)(nbr);
        let kind = rdkit_substr(line, 3, 3);
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 4 {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup line too short: '{line}' on line {line_number}"
                    )),
                );
            }
            let Some(bookmark) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let Some(index) = bookmark.checked_sub(1).map(|index| index as usize) else {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup {kind} index 0 out of range on line {line_number}"
                    )),
                );
            };
            if (kind == "SBL" && index >= bond_count) || (kind != "SBL" && index >= atom_count) {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup {kind} index {bookmark} out of range on line {line_number}"
                    )),
                );
            }
            let group = self
                .group_mut_if_present(sequence)
                .expect("presence checked before parsing vector entries");
            match kind {
                "SAL" => group.push_atom(AtomId::new(index)),
                "SBL" => group.push_bond(BondId::new(index)),
                "SPA" => group.push_parent_atom(AtomId::new(index)),
                _ => unreachable!("dispatch limits V2000 vector tags"),
            }
        }
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_smt(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        if position >= line.len() {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup line too short: '{line}' on line {line_number}"
                )),
            );
        }
        let label = &line[position..];
        let group = self
            .group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup label");
        if group.kind() == &SubstanceGroupKind::MultipleGroup {
            group.set_prop("MULT", label);
        } else {
            group.set_label(label);
        }
        Ok(())
    }

    fn parse_sdi(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        if count != 4 {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "Unexpected number of coordinates for SDI on line {line_number}"
                )),
            );
        }
        let mut coordinate = [0.0; 4];
        for value in &mut coordinate {
            let Some(parsed) =
                self.recoverable(parse_v2000_double_field(line, line_number, &mut position))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            *value = parsed;
        }
        self.group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup bracket")
            .display_mut()
            .brackets
            .push(SGroupBracket {
                points: [
                    [coordinate[0], coordinate[1], 0.0],
                    [coordinate[2], coordinate[3], 0.0],
                    [0.0, 0.0, 0.0],
                ],
            });
        Ok(())
    }

    fn parse_sbv(
        &mut self,
        line: &str,
        line_number: usize,
        bond_endpoints: &[(AtomId, AtomId)],
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SBVLine
        // RDKit❗✔️: unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: if (!ok) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️: if (!sgroup) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: unsigned int bondMark =
        // RDKit❗✔️:     ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: if (!ok) {
        // RDKit❗✔️:   sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: if (sgroup->getProp<std::string>("TYPE") == "SUP") {
        // RDKit❗✔️:   vector.x = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   vector.y = ParseSGroupDoubleField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (!ok) {
        // RDKit❗✔️:     sgroup->setIsValid(false);
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: try {
        // RDKit❗✔️:   sgroup->addCState(bond->getIdx(), vector);
        // RDKit❗✔️: } catch (const std::exception &e) {
        // RDKit❗✔️:   SGroupWarnOrThrow<>(strictParsing, e.what());
        // RDKit❗✔️:   sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        let Some(group) = self.group_mut_if_present(sequence) else {
            return Ok(());
        };
        let is_superatom = group.kind() == &SubstanceGroupKind::Superatom;
        let Some(bookmark) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        let Some(index) = bookmark.checked_sub(1).map(|index| index as usize) else {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup bond index 0 out of range on line {line_number}"
                )),
            );
        };
        let Some(&(begin, end)) = bond_endpoints.get(index) else {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup bond index {bookmark} out of range on line {line_number}"
                )),
            );
        };
        let vector = if is_superatom {
            let Some(x) =
                self.recoverable(parse_v2000_double_field(line, line_number, &mut position))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let Some(y) =
                self.recoverable(parse_v2000_double_field(line, line_number, &mut position))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            [x, y, 0.0]
        } else {
            [0.0, 0.0, 0.0]
        };
        let group = self
            .group_mut_if_present(sequence)
            .expect("presence checked before parsing CState");
        let bond = BondId::new(index);
        let bond_is_listed = group.bonds().contains(&bond);
        let begin_is_member = group.atoms().contains(&begin);
        let end_is_member = group.atoms().contains(&end);
        if !bond_is_listed || begin_is_member == end_is_member {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "Bond with index {index} is not an XBOND for current SubstanceGroup"
                )),
            );
        }
        group.push_cstate(SGroupCState { bond, vector });
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_sdt(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        let field_name = rdkit_substr(line, position, 30).trim_end().to_owned();
        position += 30;
        let field_type = rdkit_substr(line, position, 2).trim_end().to_owned();
        position += 2;
        let field_info = rdkit_substr(line, position, 20).trim_end().to_owned();
        position += 20;
        let query_type = rdkit_substr(line, position, 2).trim_end().to_owned();
        position += 2;
        let query_op = rdkit_substr(line, position, line.len().saturating_sub(position))
            .trim_end()
            .to_owned();
        let data = self
            .group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup data header")
            .data_mut();
        data.field_name = (!field_name.is_empty()).then_some(field_name);
        data.field_type = (!field_type.is_empty()).then_some(field_type);
        data.field_info = (!field_info.is_empty()).then_some(field_info);
        data.query_type = (!query_type.is_empty()).then_some(query_type);
        data.query_op = (!query_op.is_empty()).then_some(query_op);
        Ok(())
    }

    fn parse_sdd(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        if position < line.len() {
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup display data")
                .data_mut()
                .field_display = Some(line[position..].to_owned());
        }
        Ok(())
    }

    fn parse_scd_sed(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let kind = rdkit_substr(line, 3, 3);
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        if self.last_data_group != 0 && self.last_data_group != sequence {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "Found a Data Field not matching the SGroup of the last Data Field at line {line_number}"
                )),
            );
        } else if self.last_data_group == 0 && kind == "SCD" {
            self.last_data_group = sequence;
        } else if kind == "SED" {
            self.last_data_group = 0;
        }
        if self.strict_parsing && kind == "SCD" && self.scd_counter > 2 {
            return Err(SdfReadError::Parse(format!(
                "Found too many consecutive SCD lines, (#{} at line {line_number}) for SGroup {sequence}",
                self.scd_counter + 1
            )));
        }
        if position + 1 < line.len() {
            self.current_data_field
                .push_str(rdkit_substr(line, position + 1, 69));
            if kind == "SED" {
                let value = rdkit_substr(self.current_data_field.trim_end(), 0, 200).to_owned();
                self.group_mut_if_present(sequence)
                    .expect("presence checked before parsing SGroup data value")
                    .data_mut()
                    .values
                    .push(value);
                self.current_data_field.clear();
                self.scd_counter = 0;
            } else {
                self.scd_counter += 1;
            }
        }
        Ok(())
    }

    fn parse_spl(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SPL line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            // RDKit intentionally uses the throwing overload for PARENT even
            // when the surrounding parse is non-strict.
            let parent = parse_v2000_int_field(line, line_number, &mut position, false)?;
            self.parent_by_sequence.insert(sequence, parent);
        }
        Ok(())
    }

    fn parse_snc(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SNC line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let Some(component) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            if component > 256 {
                return self.invalidate_or_throw(
                    sequence,
                    SdfReadError::Parse(format!(
                        "SGroup SNC value over 256: '{component}' on line {line_number}"
                    )),
                );
            }
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup component number")
                .set_component_number(component);
        }
        Ok(())
    }

    fn parse_sap(
        &mut self,
        line: &str,
        line_number: usize,
        atom_count: usize,
    ) -> Result<(), SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseSGroupV2000SAPLine
        // RDKit❗✔️: unsigned int sgIdx = ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: if (!ok) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: SubstanceGroup *sgroup = FindSgIdx(sGroupMap, sgIdx, line);
        // RDKit❗✔️: if (!sgroup) {
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: unsigned int nent =
        // RDKit❗✔️:     ParseSGroupIntField(ok, strictParsing, text, line, pos, true);
        // RDKit❗✔️: if (!ok) {
        // RDKit❗✔️:   sgroup->setIsValid(false);
        // RDKit❗✔️:   return;
        // RDKit❗✔️: }
        // RDKit❗✔️: int lvIdx = -1;
        // RDKit❗✔️: if (text.size() < pos + 11) {
        // RDKit❗✔️:   if (strictParsing) {
        // RDKit❗✔️:     throw FileParseException(errout.str());
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     if (text.size() < pos + 4) {
        // RDKit❗✔️:       sgroup->setIsValid(false);
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     lvIdx = mol->getNumAtoms();
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: std::string id = "  ";
        // RDKit❗✔️: unsigned int aIdxMark =
        // RDKit❗✔️:     ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️: if (lvIdx == -1) {
        // RDKit❗✔️:   unsigned int lvIdxMark =
        // RDKit❗✔️:       ParseSGroupIntField(ok, strictParsing, text, line, pos);
        // RDKit❗✔️:   if (lvIdxMark != 0) {
        // RDKit❗✔️:     lvIdx = mol->getAtomWithBookmark(lvIdxMark)->getIdx();
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (text.size() >= pos + 3) {
        // RDKit❗✔️:     id = text.substr(pos + 1, 2);
        // RDKit❗✔️:     pos += 3;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: sgroup->addAttachPoint(aIdx, lvIdx, id);
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            self.invalid_sequences.insert(sequence);
            return Ok(());
        };
        for _ in 0..count {
            let missing_leaving_atom = line.len() < position + 11;
            if missing_leaving_atom && self.strict_parsing {
                return Err(SdfReadError::Parse(format!(
                    "SGroup SAP line too short: '{line}' on line {line_number}"
                )));
            }
            if missing_leaving_atom && line.len() < position + 4 {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            }
            let Some(atom_bookmark) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let atom_index = atom_bookmark.checked_sub(1).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup attach atom index 0 out of range on line {line_number}"
                ))
            })? as usize;
            if atom_index >= atom_count {
                return Err(SdfReadError::Parse(format!(
                    "SGroup attach atom index out of range on line {line_number}"
                )));
            }
            let atom = AtomId::new(atom_index);
            if missing_leaving_atom {
                self.pending_attach_points
                    .push((sequence, atom, Some("  ".to_owned())));
                continue;
            }
            let Some(leaving_bookmark) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            if leaving_bookmark as usize > atom_count {
                return Err(SdfReadError::Parse(format!(
                    "SGroup attach atom index out of range on line {line_number}"
                )));
            }
            let label = if line.len() >= position + 3 {
                let label = rdkit_substr(line, position + 1, 2).to_owned();
                position += 3;
                Some(label)
            } else {
                Some("  ".to_owned())
            };
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing attachment points")
                .push_attach_point(SGroupAttachPoint {
                    atom,
                    leaving_atom: (leaving_bookmark != 0)
                        .then(|| AtomId::new(leaving_bookmark as usize - 1)),
                    label,
                    order: None,
                });
        }
        Ok(())
        // END RDKIT CPP FUNCTION
    }

    fn parse_scl(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(sequence) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            false,
        ))?
        else {
            return Ok(());
        };
        if self.group_mut_if_present(sequence).is_none() {
            return Ok(());
        }
        position += 1;
        if position >= line.len() {
            return self.invalidate_or_throw(
                sequence,
                SdfReadError::Parse(format!(
                    "SGroup SCL line too short: '{line}' on line {line_number}"
                )),
            );
        }
        self.group_mut_if_present(sequence)
            .expect("presence checked before parsing SGroup class")
            .set_class(&line[position..]);
        Ok(())
    }

    fn parse_sbt(&mut self, line: &str, line_number: usize) -> Result<(), SdfReadError> {
        let mut position = 6;
        let Some(count) = self.recoverable(parse_v2000_int_field(
            line,
            line_number,
            &mut position,
            true,
        ))?
        else {
            return Ok(());
        };
        for _ in 0..count {
            if line.len() < position + 8 {
                let error = SdfReadError::Parse(format!(
                    "SGroup SBT line too short: '{line}' on line {line_number}"
                ));
                return if self.strict_parsing {
                    Err(error)
                } else {
                    Ok(())
                };
            }
            let Some(sequence) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                return Ok(());
            };
            if self.group_mut_if_present(sequence).is_none() {
                return Ok(());
            }
            let Some(bracket_type) = self.recoverable(parse_v2000_int_field(
                line,
                line_number,
                &mut position,
                false,
            ))?
            else {
                self.invalid_sequences.insert(sequence);
                return Ok(());
            };
            let style = match bracket_type {
                0 => SGroupBracketStyle::Bracket,
                1 => SGroupBracketStyle::Parenthesis,
                _ => {
                    return self.invalidate_or_throw(
                        sequence,
                        SdfReadError::Parse(format!(
                            "Invalid SBT value '{bracket_type}' on line {line_number}"
                        )),
                    );
                }
            };
            self.group_mut_if_present(sequence)
                .expect("presence checked before parsing SGroup bracket type")
                .set_bracket_style(style);
        }
        Ok(())
    }

    pub(super) fn finish(
        mut self,
        bond_endpoints: &[(AtomId, AtomId)],
    ) -> Result<Vec<SubstanceGroup>, SdfReadError> {
        // BEGIN RDKIT CPP FUNCTION ParseMolBlockProperties (SGroup finalization)
        // RDKit❗✔️: for (auto &sgroup : sGroupMap) {
        // RDKit❗✔️:   if (sgroup.second.getIsValid()) {
        // RDKit❗✔️:     sgroup.second.setProp("DATAFIELDS", dataFieldsMap[sgroup.first]);
        // RDKit❗✔️:     sgroup.second.setIsValid(checkAttachmentPointsAreValid(mol, sgroup));
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (sgroup.second.getIsValid()) {
        // RDKit❗✔️:     addSubstanceGroup(*mol, sgroup.second);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     if (strictParsing) {
        // RDKit❗✔️:       throw FileParseException(errout.str());
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       BOOST_LOG(rdWarningLog)
        // RDKit❗✔️:           << errout.str() << " and will be ignored" << std::endl;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // BEGIN RDKIT CPP FUNCTION checkAttachmentPointsAreValid
        // RDKit❗✔️: for (auto &attachPoint : attachPoints) {
        // RDKit❗✔️:   if (attachPoint.lvIdx == nAtoms) {
        // RDKit❗✔️:     const std::vector<unsigned int> &bonds = sgroup.second.getBonds();
        // RDKit❗✔️:     if (bonds.size() == 1) {
        // RDKit❗✔️:       const auto bond = mol->getBondWithIdx(bonds.front());
        // RDKit❗✔️:       if (bond->getBeginAtomIdx() == attachPoint.aIdx ||
        // RDKit❗✔️:           bond->getEndAtomIdx() == attachPoint.aIdx) {
        // RDKit❗✔️:         attachPoint.lvIdx = bond->getOtherAtomIdx(attachPoint.aIdx);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (attachPoint.lvIdx == nAtoms) {
        // RDKit❗✔️:     res = false;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        for (sequence, atom, label) in std::mem::take(&mut self.pending_attach_points) {
            if self.invalid_sequences.contains(&sequence) {
                continue;
            }
            let leaving_atom = self.groups.get(&sequence).and_then(|group| {
                let [bond] = group.bonds() else {
                    return None;
                };
                let &(begin, end) = bond_endpoints.get(bond.index())?;
                if begin == atom {
                    Some(end)
                } else if end == atom {
                    Some(begin)
                } else {
                    None
                }
            });
            if let Some(leaving_atom) = leaving_atom {
                self.groups
                    .get_mut(&sequence)
                    .expect("pending attachment point belongs to a parsed group")
                    .push_attach_point(SGroupAttachPoint {
                        atom,
                        leaving_atom: Some(leaving_atom),
                        label,
                        order: None,
                    });
            } else {
                self.invalid_sequences.insert(sequence);
            }
        }
        // END RDKIT CPP FUNCTION
        if let Some(sequence) = self.invalid_sequences.iter().next().copied()
            && self.strict_parsing
        {
            return Err(SdfReadError::Parse(format!("SGroup {sequence} is invalid")));
        }
        self.groups
            .retain(|sequence, _| !self.invalid_sequences.contains(sequence));
        let ids = self
            .groups
            .keys()
            .enumerate()
            .map(|(index, sequence)| (*sequence, SubstanceGroupId::new(index)))
            .collect::<BTreeMap<_, _>>();
        for (sequence, id) in &ids {
            self.groups
                .get_mut(sequence)
                .expect("key came from group map")
                .set_id(*id);
        }
        for (sequence, parent_sequence) in self.parent_by_sequence {
            let parent = *ids.get(&parent_sequence).ok_or_else(|| {
                SdfReadError::Parse(format!(
                    "SGroup {parent_sequence} referenced as parent not found."
                ))
            })?;
            self.groups
                .get_mut(&sequence)
                .ok_or_else(|| SdfReadError::Parse(format!("SGroup {sequence} missing")))?
                .set_parent(parent);
        }
        Ok(self.groups.into_values().collect())
        // END RDKIT CPP FUNCTION
    }
}

fn rdkit_sgroup_type(group: &SubstanceGroup) -> &str {
    group.props().get("TYPE").map_or_else(
        || match group.kind() {
            SubstanceGroupKind::Data => "DAT",
            SubstanceGroupKind::Superatom => "SUP",
            SubstanceGroupKind::MultipleGroup => "MUL",
            SubstanceGroupKind::StructuralRepeatUnit => "SRU",
            SubstanceGroupKind::Monomer => "MON",
            SubstanceGroupKind::Copolymer => "COP",
            SubstanceGroupKind::Crosslink => "CRO",
            SubstanceGroupKind::Graft => "GRA",
            SubstanceGroupKind::Modification => "MOD",
            SubstanceGroupKind::Mer => "MER",
            SubstanceGroupKind::AnyPolymer => "ANY",
            SubstanceGroupKind::MixtureComponent => "COM",
            SubstanceGroupKind::Mixture => "MIX",
            SubstanceGroupKind::Formulation => "FOR",
            SubstanceGroupKind::Generic(value) => value,
        },
        String::as_str,
    )
}

fn v3000_index_block<T>(
    name: &str,
    values: impl IntoIterator<Item = T>,
    index: impl Fn(T) -> usize,
) -> String {
    // BEGIN RDKIT CPP FUNCTION BuildV3000IdxVectorDataBlock
    // RDKit✔️✔️: size_t size = dataVectorEnd - dataVectorBegin;
    // RDKit✔️✔️: if (size) {
    // RDKit✔️✔️:   ret << ' ' << key << "=(" << size;
    // RDKit✔️✔️:   for (auto itr = dataVectorBegin; itr < dataVectorEnd; ++itr) {
    // RDKit✔️✔️:     ret << ' ' << 1 + *itr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ret << ')';
    // RDKit✔️✔️: }
    let indices = values.into_iter().map(index).collect::<Vec<_>>();
    if indices.is_empty() {
        return String::new();
    }
    let mut output = format!(" {name}=({}", indices.len());
    for value in indices {
        output.push_str(&format!(" {}", value + 1));
    }
    output.push(')');
    output
    // END RDKIT CPP FUNCTION
}

fn v3000_string_block(name: &str, value: Option<&str>) -> String {
    let Some(value) = value.filter(|value| !value.is_empty()) else {
        return String::new();
    };
    // BEGIN RDKIT CPP FUNCTION FormatV3000StringPropertyBlock
    // RDKit✔️✔️: bool needsQuotes = propValue.find(' ') != std::string::npos ||
    // RDKit✔️✔️:                    propValue.find('"') != std::string::npos ||
    // RDKit✔️✔️:                    propValue.find('(') != std::string::npos;
    // RDKit✔️✔️: if (needsQuotes) {
    // RDKit✔️✔️:   ret << "\"";
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (auto chr : propValue) {
    // RDKit✔️✔️:   ret << chr;
    // RDKit✔️✔️:   if (chr == '"') {
    // RDKit✔️✔️:     ret << chr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (needsQuotes) {
    // RDKit✔️✔️:   ret << "\"";
    // RDKit✔️✔️: }
    let needs_quotes = value.contains([' ', '"', '(']);
    let escaped = value.replace('"', "\"\"");
    if needs_quotes {
        format!(" {name}=\"{escaped}\"")
    } else {
        format!(" {name}={escaped}")
    }
    // END RDKIT CPP FUNCTION
}

fn add_v3000_block(block: &str, current: &mut String, output: &mut String) {
    // BEGIN RDKIT CPP FUNCTION addBlockToSGroupString
    // RDKit❗✔️: if (block.empty()) {
    // RDKit❗✔️:   return;
    // RDKit❗✔️: }
    // RDKit❗✔️: if (currentLine.length() + block.length() < 78) {
    // RDKit❗✔️:   currentLine += block;
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   os << currentLine << " -\n";
    // RDKit❗✔️:   unsigned int length = block.size();
    // RDKit❗✔️:   unsigned int start = 0;
    // RDKit❗✔️:   while (length - start >= 73) {
    // RDKit❗✔️:     os << "M  V30";
    // RDKit❗✔️:     if (start) {
    // RDKit❗✔️:       os << ' ';
    // RDKit❗✔️:     }
    // RDKit❗✔️:     os << block.substr(start, 72);
    // RDKit❗✔️:     start += 72;
    // RDKit❗✔️:     if (start < length) {
    // RDKit❗✔️:       // need to write more, so add another "-"
    // RDKit❗✔️:       os << "-\n";
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    if block.is_empty() {
        return;
    }
    if current.len() + block.len() < 78 {
        current.push_str(block);
        return;
    }
    output.push_str(current);
    output.push_str(" -\n");
    // SGroup blocks are ASCII syntax plus user strings. Split only at a UTF-8
    // character boundary while retaining RDKit's 72-byte target where possible.
    let mut remainder = block;
    let mut continued = false;
    while remainder.len() >= 73 {
        let mut split = 72;
        while !remainder.is_char_boundary(split) {
            split -= 1;
        }
        output.push_str("M  V30");
        if continued {
            output.push(' ');
        }
        output.push_str(&remainder[..split]);
        remainder = &remainder[split..];
        continued = true;
        if !remainder.is_empty() {
            output.push_str("-\n");
        }
    }
    if remainder.is_empty() {
        current.clear();
    } else {
        *current = format!("M  V30{}{remainder}", if continued { " " } else { "" });
    }
    // END RDKIT CPP FUNCTION
}

fn connection_text(value: &SGroupConnection) -> &str {
    match value {
        SGroupConnection::HeadToHead => "HH",
        SGroupConnection::HeadToTail => "HT",
        SGroupConnection::Either => "EU",
        SGroupConnection::Unknown(value) => value,
    }
}

fn bracket_style_text(value: &SGroupBracketStyle) -> &str {
    match value {
        SGroupBracketStyle::Bracket => "BRACKET",
        SGroupBracketStyle::Parenthesis => "PAREN",
        SGroupBracketStyle::None => "",
        SGroupBracketStyle::Unknown(value) => value,
    }
}

fn assigned_stereo_group_ids(groups: &[StereoGroup]) -> Vec<Option<u32>> {
    let mut assigned = groups
        .iter()
        .map(|group| match group.kind() {
            StereoGroupKind::Absolute => None,
            StereoGroupKind::Or | StereoGroupKind::And => group.id().filter(|id| *id != 0),
        })
        .collect::<Vec<_>>();
    for kind in [StereoGroupKind::Or, StereoGroupKind::And] {
        let mut used = BTreeSet::new();
        for (group, id) in groups.iter().zip(&mut assigned) {
            if group.kind() == kind
                && let Some(value) = *id
                && !used.insert(value)
            {
                *id = None;
            }
        }
        let mut next = 0_u32;
        for (group, id) in groups.iter().zip(&mut assigned) {
            if group.kind() == kind && id.is_none() {
                next += 1;
                while used.contains(&next) {
                    next += 1;
                }
                *id = Some(next);
            }
        }
    }
    assigned
}

fn write_v3000_sgroup(sequence: usize, group: &SubstanceGroup) -> String {
    // BEGIN RDKIT CPP FUNCTION GetV3000MolFileSGroupLines
    // RDKit❗✔️: std::string currLine = (boost::format("M  V30 %d %s %d") % idx %
    // RDKit❗✔️:                         sgroup.getProp<std::string>("TYPE") % id).str();
    // RDKit❗✔️: addBlockToSGroupString(
    // RDKit❗✔️:     BuildV3000IdxVectorDataBlock("ATOMS", sgroup.getAtoms()), currLine, os);
    // RDKit✔️✔️: std::string BuildV3000BondsBlock(const SubstanceGroup &sgroup) {
    // RDKit✔️✔️:   std::ostringstream ret;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto isXBond = [&sgroup](unsigned int bondIdx) {
    // RDKit✔️✔️:     return SubstanceGroup::BondType::XBOND == sgroup.getBondType(bondIdx);
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bonds = sgroup.getBonds();
    // RDKit✔️✔️:   auto first_cbond = std::stable_partition(bonds.begin(), bonds.end(), isXBond);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ret << BuildV3000IdxVectorDataBlock("XBONDS", bonds.begin(), first_cbond);
    // RDKit✔️✔️:   ret << BuildV3000IdxVectorDataBlock("CBONDS", first_cbond, bonds.end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (sgroup.hasProp("XBHEAD")) {
    // RDKit✔️✔️:     auto v = sgroup.getProp<std::vector<unsigned int>>("XBHEAD");
    // RDKit✔️✔️:     ret << BuildV3000IdxVectorDataBlock("XBHEAD", v.begin(), v.end());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (sgroup.hasProp("XBCORR")) {
    // RDKit✔️✔️:     auto v = sgroup.getProp<std::vector<unsigned int>>("XBCORR");
    // RDKit✔️✔️:     ret << BuildV3000IdxVectorDataBlock("XBCORR", v.begin(), v.end());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return ret.str();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: addBlockToSGroupString(BuildV3000BondsBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(
    // RDKit❗✔️:     BuildV3000IdxVectorDataBlock("PATOMS", sgroup.getParentAtoms()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("SUBTYPE", sgroup),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("CONNECT", sgroup),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000ParentBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CompNoBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("LABEL", sgroup),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000BracketBlock(sgroup.getBrackets()),
    // RDKit❗✔️:                        currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CStateBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000FieldDataBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000AttachPointBlock(sgroup.getAttachPoints()),
    // RDKit❗✔️:                        currLine, os);
    let mut output = String::new();
    let mut current = format!(
        "M  V30 {sequence} {} {}",
        rdkit_sgroup_type(group),
        group.external_id().unwrap_or(0)
    );
    let mut blocks = vec![v3000_index_block(
        "ATOMS",
        group.atoms().iter().copied(),
        AtomId::index,
    )];
    blocks.push(v3000_index_block(
        "XBONDS",
        group
            .bonds()
            .iter()
            .copied()
            .filter(|bond| group.bond_role(*bond) == SGroupBondRole::Crossing),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "CBONDS",
        group
            .bonds()
            .iter()
            .copied()
            .filter(|bond| group.bond_role(*bond) == SGroupBondRole::Contained),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "XBHEAD",
        group.head_crossing_bonds().iter().copied(),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "XBCORR",
        group.crossing_bond_correspondence().iter().copied(),
        BondId::index,
    ));
    blocks.push(v3000_index_block(
        "PATOMS",
        group.parent_atoms().iter().copied(),
        AtomId::index,
    ));
    blocks.push(v3000_string_block("SUBTYPE", group.subtype()));
    blocks.push(v3000_string_block(
        "MULT",
        group.props().get("MULT").map(String::as_str),
    ));
    blocks.push(v3000_string_block(
        "CONNECT",
        group.connection().map(connection_text),
    ));
    if let Some(parent) = group.parent() {
        blocks.push(format!(" PARENT={}", parent.index() + 1));
    }
    if let Some(component) = group.component_number() {
        blocks.push(format!(" COMPNO={component}"));
    }
    blocks.push(v3000_string_block("LABEL", group.label()));
    if let Some(display) = group.display() {
        for bracket in &display.brackets {
            blocks.push(format!(
                " BRKXYZ=(9 {:.4} {:.4} 0 {:.4} {:.4} 0 0 0 0)",
                bracket.points[0][0],
                bracket.points[0][1],
                bracket.points[1][0],
                bracket.points[1][1]
            ));
        }
    }
    blocks.push(v3000_string_block("ESTATE", group.expansion_state()));
    for cstate in group.cstates() {
        blocks.push(if group.kind() == &SubstanceGroupKind::Superatom {
            format!(
                " CSTATE=(4 {} {:.4} {:.4} 0)",
                cstate.bond.index() + 1,
                cstate.vector[0],
                cstate.vector[1]
            )
        } else {
            format!(" CSTATE=(1 {})", cstate.bond.index() + 1)
        });
    }
    if let Some(data) = group.data() {
        blocks.push(v3000_string_block("FIELDNAME", data.field_name.as_deref()));
        blocks.push(v3000_string_block("FIELDINFO", data.field_info.as_deref()));
        blocks.push(v3000_string_block(
            "FIELDDISP",
            data.field_display.as_deref(),
        ));
        blocks.push(v3000_string_block("QUERYTYPE", data.query_type.as_deref()));
        blocks.push(v3000_string_block("QUERYOP", data.query_op.as_deref()));
        for value in &data.values {
            blocks.push(format!(" FIELDDATA=\"{}\"", value.replace('"', "\"\"")));
        }
    }
    blocks.push(v3000_string_block("CLASS", group.class()));
    for point in group.attach_points() {
        let leaving = if point.leaving_atom == Some(point.atom) {
            "aidx".to_owned()
        } else {
            point
                .leaving_atom
                .map_or_else(|| "0".to_owned(), |atom| (atom.index() + 1).to_string())
        };
        blocks.push(format!(
            " SAP=(3 {} {leaving} {})",
            point.atom.index() + 1,
            point.label.as_deref().unwrap_or_default()
        ));
    }
    blocks.push(v3000_string_block(
        "BRKTYP",
        group.bracket_style().map(bracket_style_text),
    ));
    for block in blocks {
        add_v3000_block(&block, &mut current, &mut output);
    }
    if !current.is_empty() {
        output.push_str(&current);
        output.push('\n');
    } else if output.ends_with(" -\n") {
        output.truncate(output.len() - 3);
        output.push('\n');
    }
    output
    // END RDKIT CPP FUNCTION
}

pub(super) fn write_v3000_typed_blocks(topology: &TopologyBlock) -> String {
    let mut output = String::new();
    if !topology.substance_groups.is_empty() {
        // BEGIN RDKIT CPP FUNCTION getV3000CTAB (SGroup block)
        // RDKit✔️✔️: if (nSGroups > 0) {
        // RDKit✔️✔️:   res += "M  V30 BEGIN SGROUP\n";
        // RDKit✔️✔️:   unsigned int idx = 0;
        // RDKit✔️✔️:   for (const auto &sgroup : sgroups) {
        // RDKit✔️✔️:     res += GetV3000MolFileSGroupLines(++idx, sgroup);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   res += "M  V30 END SGROUP\n";
        // RDKit✔️✔️: }
        output.push_str("M  V30 BEGIN SGROUP\n");
        for (index, group) in topology.substance_groups.iter().enumerate() {
            output.push_str(&write_v3000_sgroup(index + 1, group));
        }
        output.push_str("M  V30 END SGROUP\n");
        // END RDKIT CPP FUNCTION
    }
    if !topology.stereo_groups.is_empty() {
        // BEGIN RDKIT CPP FUNCTION appendEnhancedStereoGroups
        // RDKit❗✔️: auto stereo_groups = tmol.getStereoGroups();
        // RDKit❗✔️: assignStereoGroupIds(stereo_groups);
        // RDKit❗✔️: res += "M  V30 BEGIN COLLECTION\n";
        // RDKit❗✔️: std::string tmp;
        // RDKit❗✔️: tmp.reserve(80);
        // RDKit❗✔️: for (auto &&group : stereo_groups) {
        // RDKit❗✔️:   tmp += "M  V30 MDLV30/";
        // RDKit❗✔️:   switch (group.getGroupType()) {
        // RDKit❗✔️:     case RDKit::StereoGroupType::STEREO_ABSOLUTE:
        // RDKit❗✔️:       tmp += "STEABS";
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDKit::StereoGroupType::STEREO_OR:
        // RDKit❗✔️:       tmp += "STEREL";
        // RDKit❗✔️:       tmp += std::to_string(group.getWriteId());
        // RDKit❗✔️:       break;
        // RDKit❗✔️:     case RDKit::StereoGroupType::STEREO_AND:
        // RDKit❗✔️:       tmp += "STERAC";
        // RDKit❗✔️:       tmp += std::to_string(group.getWriteId());
        // RDKit❗✔️:       break;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   tmp += " ATOMS=(";
        // RDKit❗✔️:   tmp += std::to_string(atomIds.size());
        // RDKit❗✔️:   for (auto &&atom : atomIds) {
        // RDKit❗✔️:     tmp += ' ';
        // RDKit❗✔️:     // atoms are 1 indexed in molfiles
        // RDKit❗✔️:     auto idxStr = std::to_string(atom + 1);
        // RDKit❗✔️:     if (tmp.size() + idxStr.size() >= 78) {
        // RDKit❗✔️:       res += tmp + "-\n";
        // RDKit❗✔️:       tmp = "M  V30 ";
        // RDKit❗✔️:     }
        // RDKit❗✔️:     tmp += idxStr;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   res += tmp + ")\n";
        // RDKit❗✔️:   tmp.clear();
        // RDKit❗✔️: }
        // RDKit❗✔️: res += tmp + "M  V30 END COLLECTION\n";
        output.push_str("M  V30 BEGIN COLLECTION\n");
        let group_ids = assigned_stereo_group_ids(&topology.stereo_groups);
        for (group, assigned_id) in topology.stereo_groups.iter().zip(group_ids) {
            let label = match group.kind() {
                StereoGroupKind::Absolute => "STEABS".to_owned(),
                StereoGroupKind::Or => format!("STEREL{}", assigned_id.expect("OR ID assigned")),
                StereoGroupKind::And => {
                    format!("STERAC{}", assigned_id.expect("AND ID assigned"))
                }
            };
            let mut current = format!("M  V30 MDLV30/{label} ATOMS=({}", group.atoms().len());
            for atom in group.atoms() {
                current.push(' ');
                let index = (atom.index() + 1).to_string();
                if current.len() + index.len() >= 78 {
                    output.push_str(&current);
                    output.push_str("-\n");
                    current.clear();
                    current.push_str("M  V30 ");
                }
                current.push_str(&index);
            }
            output.push_str(&current);
            output.push_str(")\n");
        }
        output.push_str("M  V30 END COLLECTION\n");
        // END RDKIT CPP FUNCTION
    }
    output
}

fn v2000_int(value: usize) -> String {
    format!(" {value:>3}")
}

fn v2000_count(value: usize) -> String {
    format!(" {value:>2}")
}

fn append_v2000_pairs(output: &mut String, code: &str, pairs: &[(usize, String)], per_line: usize) {
    for chunk in pairs.chunks(per_line) {
        output.push_str(&format!("M  {code}{}", v2000_count(chunk.len())));
        for (index, value) in chunk {
            output.push_str(&v2000_int(*index));
            output.push(' ');
            output.push_str(value);
        }
        output.push('\n');
    }
}

fn append_v2000_indices(
    output: &mut String,
    code: &str,
    group: usize,
    indices: impl IntoIterator<Item = usize>,
) {
    let values = indices.into_iter().collect::<Vec<_>>();
    for chunk in values.chunks(15) {
        output.push_str(&format!(
            "M  {code}{}{}",
            v2000_int(group),
            v2000_count(chunk.len())
        ));
        for index in chunk {
            output.push_str(&v2000_int(index + 1));
        }
        output.push('\n');
    }
}

fn utf8_chunks_at_most(value: &str, maximum_bytes: usize) -> Vec<&str> {
    if value.is_empty() {
        return vec![value];
    }
    let mut chunks = Vec::new();
    let mut start = 0;
    while start < value.len() {
        let mut end = (start + maximum_bytes).min(value.len());
        while end > start && !value.is_char_boundary(end) {
            end -= 1;
        }
        // `maximum_bytes` is non-zero at the only call site. This guard keeps
        // the helper total if that changes and a single scalar is wider than
        // the requested chunk size.
        if end == start {
            end = value[start..]
                .char_indices()
                .nth(1)
                .map_or(value.len(), |(offset, _)| start + offset);
        }
        chunks.push(&value[start..end]);
        start = end;
    }
    chunks
}

pub(super) fn write_v2000_sgroups(topology: &TopologyBlock) -> Result<String, SdfWriteError> {
    if topology.substance_groups.is_empty() {
        return Ok(String::new());
    }
    // BEGIN RDKIT CPP FUNCTION GetMolFileSGroupInfo
    // RDKit❗✔️: ret << BuildV2000STYLines(mol);
    // RDKit❗✔️: ret << BuildV2000SLBLines(mol);
    // RDKit❗✔️: ret << BuildV2000StringPropLines(8, mol, "SUBTYPE", "SST", 3);
    // RDKit❗✔️: ret << BuildV2000StringPropLines(8, mol, "CONNECT", "SCN", 3);
    // RDKit❗✔️: ret << BuildV2000SDSLines(mol);
    // RDKit❗✔️: ret << BuildV2000SPLLines(mol);
    // RDKit❗✔️: ret << BuildV2000SNCLines(mol);
    // RDKit❗✔️: ret << BuildV2000SBTLines(mol);
    // RDKit❗✔️: for (const auto &sgroup : getSubstanceGroups(mol)) {
    // RDKit❗✔️:   ret << BuildV2000IdxVectorDataLines(15, idx, "SAL", sgroup.getAtoms());
    // RDKit❗✔️:   ret << BuildV2000IdxVectorDataLines(15, idx, "SPA", sgroup.getParentAtoms());
    // RDKit❗✔️:   ret << BuildV2000IdxVectorDataLines(15, idx, "SBL", sgroup.getBonds());
    // RDKit❗✔️:   ret << BuildV2000SDILine(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SMTLine(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SBVLine(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SDTLine(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SDDLine(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SCDSEDLines(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SAPLines(idx, sgroup);
    // RDKit❗✔️:   ret << BuildV2000SCLLine(idx, sgroup);
    // RDKit❗✔️: }
    let mut output = String::new();
    let sty = topology
        .substance_groups
        .iter()
        .enumerate()
        .map(|(index, group)| (index + 1, format!("{:<3}", rdkit_sgroup_type(group))))
        .collect::<Vec<_>>();
    append_v2000_pairs(&mut output, "STY", &sty, 8);
    let slb = topology
        .substance_groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            group
                .external_id()
                .map(|id| (index + 1, format!("{id:>3}")))
        })
        .collect::<Vec<_>>();
    append_v2000_pairs(&mut output, "SLB", &slb, 8);
    let subtype = topology
        .substance_groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            group
                .subtype()
                .map(|value| (index + 1, format!("{value:<3}")))
        })
        .collect::<Vec<_>>();
    append_v2000_pairs(&mut output, "SST", &subtype, 8);
    let connections = topology
        .substance_groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            group
                .connection()
                .map(|value| (index + 1, format!("{:<3}", connection_text(value))))
        })
        .collect::<Vec<_>>();
    append_v2000_pairs(&mut output, "SCN", &connections, 8);
    let expanded = topology
        .substance_groups
        .iter()
        .enumerate()
        .filter(|(_, group)| group.expansion_state() == Some("E"))
        .map(|(index, _)| index + 1)
        .collect::<Vec<_>>();
    for chunk in expanded.chunks(15) {
        output.push_str(&format!("M  SDS EXP{}", v2000_count(chunk.len())));
        for index in chunk {
            output.push_str(&v2000_int(*index));
        }
        output.push('\n');
    }
    let parents = topology
        .substance_groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            group
                .parent()
                .map(|parent| (index + 1, format!("{:>3}", parent.index() + 1)))
        })
        .collect::<Vec<_>>();
    append_v2000_pairs(&mut output, "SPL", &parents, 8);
    let components = topology
        .substance_groups
        .iter()
        .enumerate()
        .filter_map(|(index, group)| {
            group
                .component_number()
                .map(|value| (index + 1, format!("{value:>3}")))
        })
        .collect::<Vec<_>>();
    append_v2000_pairs(&mut output, "SNC", &components, 8);
    let mut brackets = Vec::new();
    for (index, group) in topology.substance_groups.iter().enumerate() {
        let Some(style) = group.bracket_style() else {
            continue;
        };
        let value = match style {
            SGroupBracketStyle::Bracket => 0,
            SGroupBracketStyle::Parenthesis => 1,
            other => {
                return Err(SdfWriteError::SubstanceGroup(format!(
                    "V2000 cannot encode bracket style {other:?}"
                )));
            }
        };
        brackets.push((index + 1, format!("{value:>3}")));
    }
    append_v2000_pairs(&mut output, "SBT", &brackets, 8);

    for (zero_index, group) in topology.substance_groups.iter().enumerate() {
        let index = zero_index + 1;
        append_v2000_indices(
            &mut output,
            "SAL",
            index,
            group.atoms().iter().map(|atom| atom.index()),
        );
        append_v2000_indices(
            &mut output,
            "SPA",
            index,
            group.parent_atoms().iter().map(|atom| atom.index()),
        );
        append_v2000_indices(
            &mut output,
            "SBL",
            index,
            group.bonds().iter().map(|bond| bond.index()),
        );
        if let Some(display) = group.display() {
            for bracket in &display.brackets {
                output.push_str(&format!(
                    "M  SDI{}{}{:>10.4}{:>10.4}{:>10.4}{:>10.4}\n",
                    v2000_int(index),
                    v2000_count(4),
                    bracket.points[0][0],
                    bracket.points[0][1],
                    bracket.points[1][0],
                    bracket.points[1][1]
                ));
            }
        }
        let label = if group.kind() == &SubstanceGroupKind::MultipleGroup {
            group.props().get("MULT").map(String::as_str)
        } else {
            group.label()
        };
        if let Some(label) = label {
            output.push_str(&format!("M  SMT{} {}\n", v2000_int(index), label));
        }
        for cstate in group.cstates() {
            output.push_str(&format!(
                "M  SBV{}{}",
                v2000_int(index),
                v2000_int(cstate.bond.index() + 1)
            ));
            if group.kind() == &SubstanceGroupKind::Superatom {
                output.push_str(&format!(
                    "{:>10.4}{:>10.4}",
                    cstate.vector[0], cstate.vector[1]
                ));
            }
            output.push('\n');
        }
        if let Some(data) = group.data() {
            if let Some(field_name) = &data.field_name {
                output.push_str(&format!(
                    "M  SDT{} {:<30}{:<2}{:<20}{:<2}{}\n",
                    v2000_int(index),
                    field_name,
                    data.field_type.as_deref().unwrap_or("T"),
                    data.field_info.as_deref().unwrap_or_default(),
                    data.query_type.as_deref().unwrap_or_default(),
                    data.query_op.as_deref().unwrap_or_default()
                ));
            }
            if let Some(display) = &data.field_display {
                output.push_str(&format!("M  SDD{} {}\n", v2000_int(index), display));
            }
            for value in &data.values {
                if value.len() > 200 {
                    return Err(SdfWriteError::SubstanceGroup(format!(
                        "data field in SGroup {index} is longer than 200 bytes"
                    )));
                }
                let chunks = utf8_chunks_at_most(value, 69);
                for (chunk_index, chunk) in chunks.iter().enumerate() {
                    let code = if chunk_index + 1 == chunks.len() {
                        "SED"
                    } else {
                        "SCD"
                    };
                    output.push_str(&format!("M  {code}{} {chunk}\n", v2000_int(index)));
                }
            }
        }
        for chunk in group.attach_points().chunks(6) {
            output.push_str(&format!(
                "M  SAP{}{}",
                v2000_int(index),
                v2000_count(chunk.len())
            ));
            for point in chunk {
                output.push_str(&v2000_int(point.atom.index() + 1));
                output.push_str(&v2000_int(
                    point.leaving_atom.map_or(0, |atom| atom.index() + 1),
                ));
                output.push(' ');
                output.push_str(&format!(
                    "{:<2}",
                    point.label.as_deref().unwrap_or_default()
                ));
            }
            output.push('\n');
        }
        if let Some(class) = group.class() {
            output.push_str(&format!("M  SCL{} {}\n", v2000_int(index), class));
        }
    }
    Ok(output)
    // END RDKIT CPP FUNCTION
}
