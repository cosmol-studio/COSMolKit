use crate::CxParseError;

pub(crate) fn read_number(text: &str, cursor: &mut usize) -> Result<usize, CxParseError> {
    // RDKit✔️🔝: template <typename Iterator>
    // RDKit✔️🔝: bool read_int(Iterator &first, Iterator last, unsigned int &res) {
    // RDKit✔️🔝:   std::string num = "";
    // RDKit✔️🔝:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️🔝:     num += *first;
    // RDKit✔️🔝:     ++first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   if (num.empty()) {
    // RDKit✔️🔝:     return false;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   res = boost::lexical_cast<unsigned int>(num);
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    // The Rust scanner parses the same single digit run directly from the
    // borrowed input slice, avoiding the source's incrementally grown string
    // while preserving one forward pass and unsigned-int range semantics.
    let start = *cursor;
    while text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
        *cursor += 1;
    }
    text[start..*cursor]
        .parse::<u32>()
        .map(|value| value as usize)
        .map_err(|_| CxParseError::new(start, "invalid CX integer"))
}

pub(crate) fn read_pair(
    text: &str,
    cursor: &mut usize,
    separator: u8,
) -> Result<(usize, usize), CxParseError> {
    // RDKit✔️✔️: template <typename Iterator>
    // RDKit✔️✔️: bool read_int_pair(Iterator &first, Iterator last, unsigned int &n1,
    // RDKit✔️✔️:                    unsigned int &n2, char sep = '.') {
    // RDKit✔️✔️:   if (!read_int(first, last, n1)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != sep) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return read_int(first, last, n2);
    // RDKit✔️✔️: }
    let first = read_number(text, cursor)?;
    expect_byte(text, cursor, separator)?;
    let second = read_number(text, cursor)?;
    Ok((first, second))
}

pub(crate) fn parse_delimited_number_list(
    text: &str,
    cursor: &mut usize,
    separator: u8,
) -> Result<Vec<usize>, CxParseError> {
    // RDKit✔️🔝: template <typename Iterator>
    // RDKit✔️🔝: bool read_int_list(Iterator &first, Iterator last,
    // RDKit✔️🔝:                    std::vector<unsigned int> &res, char sep = ',') {
    // RDKit✔️🔝:   while (1) {
    // RDKit✔️🔝:     std::string num = "";
    // RDKit✔️🔝:     while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️🔝:       num += *first;
    // RDKit✔️🔝:       ++first;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (!num.empty()) {
    // RDKit✔️🔝:       res.push_back(boost::lexical_cast<unsigned int>(num));
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (first >= last || *first != sep) {
    // RDKit✔️🔝:       break;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     ++first;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    // Each number reuses `read_number`, so the Rust path avoids one temporary
    // owned string per list element while retaining linear scan and order.
    let mut values = Vec::new();
    loop {
        if text.as_bytes().get(*cursor).is_some_and(u8::is_ascii_digit) {
            values.push(read_number(text, cursor)?);
        }
        if text.as_bytes().get(*cursor) != Some(&separator) {
            break;
        }
        *cursor += 1;
    }
    Ok(values)
}

pub(crate) fn read_text_to(
    text: &str,
    cursor: &mut usize,
    delimiters: &[u8],
) -> Result<String, CxParseError> {
    // RDKit✔️✔️: template <typename Iterator>
    // RDKit✔️✔️: std::string read_text_to(Iterator &first, Iterator last, std::string delims) {
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:   Iterator start = first;
    // RDKit✔️✔️:   // EFF: there are certainly faster ways to do this
    // RDKit✔️✔️:   while (first <= last && delims.find_first_of(*first) == std::string::npos) {
    // RDKit✔️✔️:     if (*first == '&' && std::distance(first, last) > 2 &&
    // RDKit✔️✔️:         *(first + 1) == '#') {
    // RDKit✔️✔️:       // escaped char
    // RDKit✔️✔️:       if (start != first) {
    // RDKit✔️✔️:         res += std::string(start, first);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Iterator next = first + 2;
    // RDKit✔️✔️:       while (next != last && *next >= '0' && *next <= '9') {
    // RDKit✔️✔️:         ++next;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (next == last || *next != ';') {
    // RDKit✔️✔️:         throw RDKit::SmilesParseException(
    // RDKit✔️✔️:             "failure parsing CXSMILES extensions: quoted block not terminated "
    // RDKit✔️✔️:             "with ';'");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (next > first + 2) {
    // RDKit✔️✔️:         std::string blk = std::string(first + 2, next);
    // RDKit✔️✔️:         res += (char)(boost::lexical_cast<int>(blk));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       first = next + 1;
    // RDKit✔️✔️:       start = first;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (start != first) {
    // RDKit✔️✔️:     res += std::string(start, first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let mut result = String::new();
    let mut segment_start = *cursor;
    while let Some(&byte) = text.as_bytes().get(*cursor) {
        if delimiters.contains(&byte) {
            break;
        }
        if byte == b'&' && text.as_bytes().get(*cursor + 1) == Some(&b'#') {
            result.push_str(&text[segment_start..*cursor]);
            let entity_start = *cursor;
            let mut next = *cursor + 2;
            while text.as_bytes().get(next).is_some_and(u8::is_ascii_digit) {
                next += 1;
            }
            if text.as_bytes().get(next) != Some(&b';') {
                return Err(CxParseError::new(
                    entity_start,
                    "failure parsing CXSMILES extensions: quoted block not terminated with ';'",
                ));
            }
            if next > entity_start + 2 {
                // The pinned target has eight-bit signed `char`: RDKit first
                // parses the full signed-int domain and its cast then retains
                // the low byte (for example, 256 -> NUL and 321 -> `A`).
                // `String` cannot hold RDKit's raw 0x80..=0xff bytes, so the
                // detached text model lifts that byte reversibly to U+0080..=
                // U+00FF instead of rejecting it or creating invalid UTF-8.
                let value = text[entity_start + 2..next]
                    .parse::<i32>()
                    .map_err(|_| CxParseError::new(entity_start, "invalid CX character code"))?;
                result.push(char::from(value as u8));
            }
            *cursor = next + 1;
            segment_start = *cursor;
        } else {
            *cursor += 1;
        }
    }
    result.push_str(&text[segment_start..*cursor]);
    Ok(result)
}

pub(crate) fn expect_byte(
    text: &str,
    cursor: &mut usize,
    expected: u8,
) -> Result<(), CxParseError> {
    if text.as_bytes().get(*cursor) == Some(&expected) {
        *cursor += 1;
        Ok(())
    } else {
        Err(CxParseError::new(
            *cursor,
            format!("expected '{}', found CX syntax mismatch", expected as char),
        ))
    }
}
