//! Private bio-gated numeric owner (BIO-CID-NUM N02): the single
//! factored home of the pinned Gemmi `fast_from_chars` port consumed by
//! the PDB reader and, later, CID numeric parsing. No cross-crate export,
//! no second parser, no generic framework.

/// Source-shaped fast_float prefix scan over the remaining field bytes.
pub(crate) fn parse_fast_float_prefix(input: &[u8]) -> Option<(f64, usize)> {
    let negative = input.first() == Some(&b'-');
    let unsigned = if negative { &input[1..] } else { input };

    // fast_float source: third_party/gemmi/include/gemmi/third_party/fast_float.h
    // fast_float❗❗: if (fastfloat_strncasecmp3(first, str_const_nan<UC>())) {
    // fast_float❗❗:   answer.ptr = (first += 3);
    // fast_float❗❗:   value = minusSign ? -std::numeric_limits<T>::quiet_NaN()
    // fast_float❗❗:                     : std::numeric_limits<T>::quiet_NaN();
    // fast_float❗❗:   if (first != last && *first == UC('(')) {
    // fast_float❗❗:     for (UC const *ptr = first + 1; ptr != last; ++ptr) {
    // fast_float❗❗:       if (*ptr == UC(')')) {
    // fast_float❗❗:         answer.ptr = ptr + 1; // valid nan(n-char-seq-opt)
    // fast_float❗❗:         break;
    // fast_float❗❗:       } else if (!((UC('a') <= *ptr && *ptr <= UC('z')) ||
    // fast_float❗❗:                  (UC('A') <= *ptr && *ptr <= UC('Z')) ||
    // fast_float❗❗:                  (UC('0') <= *ptr && *ptr <= UC('9')) || *ptr == UC('_')))
    // fast_float❗❗:         break; // forbidden char, not nan(n-char-seq-opt)
    // fast_float❗❗:     }
    // fast_float❗❗:   }
    // fast_float❗❗:   return answer;
    // fast_float❗❗: }
    // fast_float❗❗: if (fastfloat_strncasecmp3(first, str_const_inf<UC>())) {
    // fast_float❗❗:   if ((last - first >= 8) &&
    // fast_float❗❗:       fastfloat_strncasecmp5(first + 3, str_const_inf<UC>() + 3)) {
    // fast_float❗❗:     answer.ptr = first + 8;
    // fast_float❗❗:   } else {
    // fast_float❗❗:     answer.ptr = first + 3;
    // fast_float❗❗:   }
    // fast_float❗❗:   value = minusSign ? -std::numeric_limits<T>::infinity()
    // fast_float❗❗:                     : std::numeric_limits<T>::infinity();
    // fast_float❗❗:   return answer;
    // fast_float❗❗: }
    if unsigned
        .get(..3)
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"nan"))
    {
        let sign_length = usize::from(negative);
        let mut consumed = 3;
        if unsigned.get(consumed) == Some(&b'(') {
            let mut position = consumed + 1;
            while let Some(&byte) = unsigned.get(position) {
                if byte == b')' {
                    consumed = position + 1;
                    break;
                }
                if !byte.is_ascii_alphanumeric() && byte != b'_' {
                    break;
                }
                position += 1;
            }
        }
        return Some((
            if negative { -f64::NAN } else { f64::NAN },
            sign_length + consumed,
        ));
    }
    if unsigned
        .get(..3)
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"inf"))
    {
        let consumed = if unsigned
            .get(..8)
            .is_some_and(|prefix| prefix.eq_ignore_ascii_case(b"infinity"))
        {
            8
        } else {
            3
        };
        return Some((
            if negative {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            },
            consumed + usize::from(negative),
        ));
    }

    // fast_float source: third_party/gemmi/include/gemmi/third_party/fast_float.h
    // fast_float❗❗: answer.negative = (*p == UC('-'));
    // fast_float❗❗: if ((*p == UC('-')) || (uint64_t(fmt & chars_format::allow_leading_plus) &&
    // fast_float❗❗:                           !basic_json_fmt && *p == UC('+'))) {
    // fast_float❗❗:   ++p;
    // fast_float❗❗: }
    // fast_float❗❗: uint64_t i = 0; // an unsigned int avoids signed overflows (which are bad)
    // fast_float❗❗: while ((p != pend) && is_integer(*p)) {
    // fast_float❗❗:   i = 10 * i +
    // fast_float❗❗:       uint64_t(*p - UC('0')); // might overflow, we will handle the overflow later
    // fast_float❗❗:   ++p;
    // fast_float❗❗: }
    // fast_float❗❗: bool const has_decimal_point = (p != pend) && (*p == decimal_point);
    // fast_float❗❗: if (has_decimal_point) {
    // fast_float❗❗:   ++p;
    // fast_float❗❗:   while ((p != pend) && is_integer(*p)) {
    // fast_float❗❗:     uint8_t digit = uint8_t(*p - UC('0'));
    // fast_float❗❗:     ++p;
    // fast_float❗❗:     i = i * 10 + digit;
    // fast_float❗❗:   }
    // fast_float❗❗: }
    // fast_float❗❗: if ((uint64_t(fmt & chars_format::scientific) && (p != pend) &&
    // fast_float❗❗:      ((UC('e') == *p) || (UC('E') == *p))) ||
    // fast_float❗❗:     (uint64_t(fmt & detail::basic_fortran_fmt) && (p != pend) &&
    // fast_float❗❗:      ((UC('+') == *p) || (UC('-') == *p) || (UC('d') == *p) ||
    // fast_float❗❗:       (UC('D') == *p)))) {
    // fast_float❗❗:   UC const *location_of_e = p;
    // fast_float❗❗:   if ((p == pend) || !is_integer(*p)) {
    // fast_float❗❗:     if (!uint64_t(fmt & chars_format::fixed)) {
    // fast_float❗❗:       return report_parse_error<UC>(p,
    // fast_float❗❗:                                     parse_error::missing_exponential_part);
    // fast_float❗❗:     }
    // fast_float❗❗:     // Otherwise, we will be ignoring the 'e'.
    // fast_float❗❗:     p = location_of_e;
    // fast_float❗❗:   }
    // fast_float❗❗:   answer.lastmatch = p;
    // fast_float❗❗:   answer.valid = true;
    // fast_float source: third_party/gemmi/include/gemmi/third_party/fast_float.h,
    // `from_chars_advanced(parsed_number_string_t&, T&)`.
    // fast_float❗❗:   answer.ptr = pns.lastmatch;
    // Behavior review: `fast_from_chars` has already removed at most one `+`.
    // This scan accepts the source general-format decimal prefix and leaves an
    // incomplete `e/E` suffix unconsumed; `FromStr` receives only the proven
    // ASCII token and supplies correctly rounded binary64 conversion.
    // Complexity review: one bounded byte scan and one standard-library parse
    // are both linear in the fixed field width, with no owned string buffer.
    let mut index = usize::from(negative);
    let integer_start = index;
    while input.get(index).is_some_and(u8::is_ascii_digit) {
        index += 1;
    }
    let mut has_digits = index != integer_start;

    if input.get(index) == Some(&b'.') {
        index += 1;
        let fraction_start = index;
        while input.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        has_digits |= index != fraction_start;
    }
    if !has_digits {
        return None;
    }

    if matches!(input.get(index), Some(b'e' | b'E')) {
        let exponent_start = index;
        index += 1;
        if matches!(input.get(index), Some(b'+' | b'-')) {
            index += 1;
        }
        let exponent_digits_start = index;
        while input.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        if index == exponent_digits_start {
            index = exponent_start;
        }
    }

    let token = std::str::from_utf8(&input[..index])
        .expect("the source numeric grammar accepts ASCII bytes only");
    Some((
        token
            .parse::<f64>()
            .expect("the validated source decimal prefix parses as binary64"),
        index,
    ))
}

/// Gemmi `fast_from_chars` over a fixed NUL-terminated PDB field.
pub(crate) fn gemmi_fast_atof_with_end(field: &[u8]) -> (f64, usize) {
    // Gemmi source: third_party/gemmi/include/gemmi/atof.hpp.
    // Gemmi❗✔️: inline from_chars_result fast_from_chars(const char* start, const char* end, double& d) {
    // Gemmi❗✔️:   while (start < end && is_space(*start))
    // Gemmi❗✔️:     ++start;
    // Gemmi❗✔️:   if (start < end && *start == '+')
    // Gemmi❗✔️:     ++start;
    // Gemmi❗✔️:   return fast_float::from_chars(start, end, d);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline double fast_atof(const char* p, const char** endptr=nullptr) {
    // Gemmi❗✔️:   double d = 0;
    // Gemmi❗✔️:   auto result = fast_from_chars(p, d);
    // Gemmi❗✔️:   if (endptr)
    // Gemmi❗✔️:     *endptr = result.ptr;
    // Gemmi❗✔️:   return d;
    // Gemmi❗✔️: }
    // Behavior review: stop at the first source NUL, skip Gemmi C-locale
    // whitespace and at most one leading plus, and report the parser's actual
    // end position even when no number converts. Invalid input keeps +0.0;
    // numeric range status is not exposed by fast_atof's returned value.
    // Complexity review: one prefix whitespace scan plus the existing
    // source-shaped fast-float prefix scan/conversion; no allocation and
    // linear work in this PDB record's bounded line length.
    // BIO-CID-NUM N04: delegated to the typed NUL-terminated overload;
    // fast_atof's observable surface (value with +0.0 prior, actual end
    // position even without conversion) is unchanged.
    let outcome = fast_from_chars_cstring_typed(field, 0.0);
    (outcome.value, outcome.consumed)
}

/// Pinned `fast_float` error categories reachable by this port
/// (BIO-CID-NUM N03): `std::errc::invalid_argument` (no number converts)
/// and `std::errc::result_out_of_range` (magnitude outside binary64).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FastFromCharsCategory {
    InvalidArgument,
    ResultOutOfRange,
}

/// Typed `from_chars` outcome (BIO-CID-NUM N03) preserving the pinned
/// `from_chars_result_t` fields the CID consumer needs: the consumed
/// offset (`ptr`), the error category (`ec`), and whether the output
/// value was assigned. The caller's prior value is retained on source
/// nonassignment (the source leaves `d` untouched on error).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FastFromCharsOutcome {
    /// Assigned value, or the caller's prior value when `assigned` is
    /// false.
    pub(crate) value: f64,
    /// Consumed byte count (`ptr` as an offset within the input slice).
    pub(crate) consumed: usize,
    /// `ec`: `None` is `std::errc()` success.
    pub(crate) category: Option<FastFromCharsCategory>,
    /// Whether the source assigned the value.
    pub(crate) assigned: bool,
}

/// Typed `fast_float::from_chars` over the remaining field bytes
/// (BIO-CID-NUM N03).
pub(crate) fn fast_float_from_chars_typed(input: &[u8], prior: f64) -> FastFromCharsOutcome {
    // fast_float✔️✔️: template <typename UC> struct from_chars_result_t {
    // fast_float✔️✔️:   UC const *ptr;
    // fast_float✔️✔️:   std::errc ec;
    // fast_float✔️✔️:   constexpr explicit operator bool() const noexcept {
    // fast_float✔️✔️:     return ec == std::errc();
    // fast_float✔️✔️:   }
    // fast_float✔️✔️: };
    // fast_float✔️✔️: if (!pns.valid) {
    // fast_float✔️✔️:   if (uint64_t(fmt & chars_format::no_infnan)) {
    // fast_float✔️✔️:     answer.ec = std::errc::invalid_argument;
    // fast_float✔️✔️:     answer.ptr = first;
    // fast_float✔️✔️:     return answer;
    // fast_float✔️✔️:   } else {
    // fast_float✔️✔️:     return detail::parse_infnan(first, last, value, fmt);
    // fast_float✔️✔️:   }
    // fast_float✔️✔️: }
    // fast_float✔️✔️:   to_float(pns.negative, am, value);
    // fast_float✔️✔️:   // Test for over/underflow.
    // fast_float✔️✔️:   if ((pns.mantissa != 0 && am.mantissa == 0 && am.power2 == 0) ||
    // fast_float✔️✔️:       am.power2 == binary_format<T>::infinite_power()) {
    // fast_float✔️✔️:     answer.ec = std::errc::result_out_of_range;
    // fast_float✔️✔️:   }
    // fast_float✔️✔️:   return answer;
    // Behavior (BIO-NUM-RANGE correction, fast_float.h:4728-4734): the
    // decimal path calls to_float BEFORE the range test, so the value is
    // ALWAYS assigned there — decimal overflow reports
    // result_out_of_range WITH the assigned signed infinity, and a
    // nonzero mantissa rounding to zero reports result_out_of_range WITH
    // the assigned signed zero. The zero-mantissa distinction is the
    // source's own: `pns.mantissa != 0` is the parsed significant-digit
    // state, so "0e999" (mantissa exactly zero) stays an ordinary
    // success, while "1e-400" carries the range status. Detection here
    // derives from the token's digit state (any nonzero digit in the
    // integer/fraction parts) and the converted magnitude — never a
    // blanket zero-value guard. Literal inf/nan spellings remain
    // separate parse_infnan SUCCESSES (assigned, ec = 0). Only a failed
    // decimal lex (invalid_argument) leaves the value unassigned and
    // retains the caller's prior. The Rust stdlib f64 conversion supplies
    // the signed saturating magnitude (±inf / ±0), which is exactly the
    // assigned value the pinned to_float produces.
    // Complexity: one lexing pass plus one correctly rounded binary64
    // conversion of the validated ASCII token; no allocation.
    match parse_fast_float_prefix(input) {
        None => FastFromCharsOutcome {
            value: prior,
            consumed: 0,
            category: Some(FastFromCharsCategory::InvalidArgument),
            assigned: false,
        },
        Some((value, consumed)) => {
            let literal_special = {
                let unsigned = if input.first() == Some(&b'-') {
                    &input[1..]
                } else {
                    input
                };
                unsigned.get(..3).is_some_and(|prefix| {
                    prefix.eq_ignore_ascii_case(b"inf") || prefix.eq_ignore_ascii_case(b"nan")
                })
            };
            if literal_special {
                return FastFromCharsOutcome {
                    value,
                    consumed,
                    category: None,
                    assigned: true,
                };
            }
            // Source range status: overflow (infinite_power) or nonzero
            // mantissa rounded to zero. The mantissa-nonzero test reads
            // the token's significant digits, mirroring pns.mantissa.
            // pns.mantissa covers ONLY the significant digits (integer
            // and fraction parts), never the exponent digits.
            let mantissa_end = input[..consumed]
                .iter()
                .position(|byte| *byte == b'e' || *byte == b'E')
                .unwrap_or(consumed);
            let mantissa_nonzero = input[..mantissa_end]
                .iter()
                .any(|byte| byte.is_ascii_digit() && *byte != b'0');
            let out_of_range = value.is_infinite() || (mantissa_nonzero && value == 0.0);
            FastFromCharsOutcome {
                value,
                consumed,
                category: if out_of_range {
                    Some(FastFromCharsCategory::ResultOutOfRange)
                } else {
                    None
                },
                assigned: true,
            }
        }
    }
}

/// Typed Gemmi `fast_from_chars` bounded overload (BIO-CID-NUM N04):
/// skips Gemmi C-locale whitespace and at most one leading plus, then
/// applies the typed canonical owner; `consumed` is an offset within the
/// ORIGINAL input (skip included). No NUL termination: the bounded
/// overload scans exactly `[start, end)`.
pub(crate) fn fast_from_chars_bounded_typed(input: &[u8], prior: f64) -> FastFromCharsOutcome {
    // Gemmi✔️✔️: inline from_chars_result fast_from_chars(const char* start, const char* end, double& d) {
    // Gemmi✔️✔️:   while (start < end && is_space(*start))
    // Gemmi✔️✔️:     ++start;
    // Gemmi✔️✔️:   if (start < end && *start == '+')
    // Gemmi✔️✔️:     ++start;
    // Gemmi✔️✔️:   return fast_float::from_chars(start, end, d);
    // Gemmi✔️✔️: }
    // Behavior: exact source order — whitespace bytes 9-13 and 32 only
    // (the pinned is_space table), then ONE optional plus (a second plus
    // or any minus belongs to the number grammar). Errors keep the
    // source's ptr semantics expressed as original-input offsets: on
    // invalid_argument the ptr is the post-skip position (the inner
    // result's first), which this port reports as the skip length alone
    // in `consumed`.
    // Complexity: one prefix scan plus the canonical owner; no
    // allocation.
    let mut start = 0;
    while start < input.len() && crate::bio_pdb::gemmi_is_space(input[start]) {
        start += 1;
    }
    if input.get(start) == Some(&b'+') {
        start += 1;
    }
    let inner = fast_float_from_chars_typed(&input[start..], prior);
    FastFromCharsOutcome {
        value: inner.value,
        consumed: start + inner.consumed,
        category: inner.category,
        assigned: inner.assigned,
    }
}

/// Typed Gemmi `fast_from_chars` NUL-terminated overload (BIO-CID-NUM
/// N04): stops the field at the first NUL (`std::strlen`), then applies
/// the bounded semantics over the remaining bytes.
pub(crate) fn fast_from_chars_cstring_typed(field: &[u8], prior: f64) -> FastFromCharsOutcome {
    // Gemmi✔️✔️: inline from_chars_result fast_from_chars(const char* start, double& d) {
    // Gemmi✔️✔️:   while (is_space(*start))
    // Gemmi✔️✔️:     ++start;
    // Gemmi✔️✔️:   if (*start == '+')
    // Gemmi✔️✔️:     ++start;
    // Gemmi✔️✔️:   return fast_float::from_chars(start, start + std::strlen(start), d);
    // Gemmi✔️✔️: }
    // Behavior: the C-string overload first computes strlen — a NUL
    // anywhere ends the field — then shares the bounded wrapper's skip
    // and conversion. Offsets remain original-input positions.
    // Complexity: one NUL scan plus the bounded path.
    let c_end = field
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(field.len());
    fast_from_chars_bounded_typed(&field[..c_end], prior)
}

#[cfg(test)]
mod bio_cid_num_n05_tests {
    use super::parse_fast_float_prefix;

    // Expectations derived from pinned parse_number_string
    // (fast_float.h:2085-2130) under the default general format:
    // fixed+scientific, so an incomplete exponent rolls back to the
    // marker, and no hexadecimal grammar exists (chars_format::hex is
    // never passed).
    #[test]
    fn bio_cid_num_n05_integer_fraction_exponent_forms() {
        assert_eq!(parse_fast_float_prefix(b"42"), Some((42.0, 2)));
        assert_eq!(parse_fast_float_prefix(b"42."), Some((42.0, 3)));
        assert_eq!(parse_fast_float_prefix(b".5"), Some((0.5, 2)));
        assert_eq!(parse_fast_float_prefix(b"1.25"), Some((1.25, 4)));
        assert_eq!(parse_fast_float_prefix(b"1e3"), Some((1000.0, 3)));
        assert_eq!(parse_fast_float_prefix(b"1E3"), Some((1000.0, 3)));
        assert_eq!(parse_fast_float_prefix(b"1e+3"), Some((1000.0, 4)));
        assert_eq!(parse_fast_float_prefix(b"1e-3"), Some((0.001, 4)));
        assert_eq!(parse_fast_float_prefix(b"-2.5e2"), Some((-250.0, 6)));
    }

    #[test]
    fn bio_cid_num_n05_incomplete_exponent_rolls_back() {
        // "1e", "1e+", "1e-" have no exponent digits; general format
        // allows fixed, so the marker is ignored and the token ends
        // before it (p = location_of_e).
        assert_eq!(parse_fast_float_prefix(b"1e"), Some((1.0, 1)));
        assert_eq!(parse_fast_float_prefix(b"1e+"), Some((1.0, 1)));
        assert_eq!(parse_fast_float_prefix(b"1e-"), Some((1.0, 1)));
        assert_eq!(parse_fast_float_prefix(b"2.5e"), Some((2.5, 3)));
        // A bare marker with nothing before it is not a number.
        assert_eq!(parse_fast_float_prefix(b"e5"), None);
        assert_eq!(parse_fast_float_prefix(b".e5"), None);
    }

    #[test]
    fn bio_cid_num_n05_no_hex_grammar_and_trailing_characters() {
        // "0x1A": the decimal lexer consumes only the digit "0"; "x" is
        // an ordinary terminator (the hex branch requires
        // chars_format::hex, never passed by the default format).
        assert_eq!(parse_fast_float_prefix(b"0x1A"), Some((0.0, 1)));
        assert_eq!(parse_fast_float_prefix(b"0X1A"), Some((0.0, 1)));
        // Trailing non-digit characters stop the token without error.
        assert_eq!(parse_fast_float_prefix(b"1.5 "), Some((1.5, 3)));
        assert_eq!(parse_fast_float_prefix(b"1.5x"), Some((1.5, 3)));
        assert_eq!(parse_fast_float_prefix(b"12;34"), Some((12.0, 2)));
        // The exponent consumes only digits after one optional sign;
        // "1e3z" keeps the "z" unconsumed.
        assert_eq!(parse_fast_float_prefix(b"1e3z"), Some((1000.0, 3)));
    }
}

/// BIO-NUM-RANGE-MATRIX prior set: a finite value, negative zero
/// (distinct sign bit), and a NaN with a specific payload — compared by
/// exact bits, never by is_nan alone.
#[cfg(test)]
fn matrix_priors() -> [(&'static str, f64); 3] {
    [
        ("finite", -7.25),
        ("neg_zero", -0.0),
        ("nan_payload", f64::from_bits(0x7ff8_1234_5678_9abc)),
    ]
}

#[cfg(test)]
mod bio_cid_num_n06_tests {
    use super::parse_fast_float_prefix;

    // Expectations derived from pinned parse_infnan (fast_float.h:4469-4522):
    // one optional '-' ONLY under the default format ("C++17 20.19.3.(7.1)
    // explicitly forbids '+' sign here"), then case-insensitive "inf" with
    // the 8-byte "inity" extension requiring last-first >= 8.
    #[test]
    fn bio_cid_num_n06_inf_case_classes_and_spellings() {
        // All case classes of the 3-letter spelling.
        for token in [&b"inf"[..], b"Inf", b"iNf", b"INF", b"inF", b"Inf"] {
            assert_eq!(
                parse_fast_float_prefix(token),
                Some((f64::INFINITY, 3)),
                "{token:?}"
            );
        }
        // Full 8-letter spelling in mixed case.
        for token in [&b"infinity"[..], b"INFINITY", b"Infinity", b"iNfInItY"] {
            assert_eq!(
                parse_fast_float_prefix(token),
                Some((f64::INFINITY, 8)),
                "{token:?}"
            );
        }
        // Suffixes stay unconsumed; short prefixes cap at 3 bytes.
        assert_eq!(
            parse_fast_float_prefix(b"infinit"),
            Some((f64::INFINITY, 3))
        );
        assert_eq!(parse_fast_float_prefix(b"infin"), Some((f64::INFINITY, 3)));
        assert_eq!(parse_fast_float_prefix(b"inf!"), Some((f64::INFINITY, 3)));
        assert_eq!(
            parse_fast_float_prefix(b"infinityx"),
            Some((f64::INFINITY, 8))
        );
        assert_eq!(
            parse_fast_float_prefix(b"Infinity)"),
            Some((f64::INFINITY, 8))
        );
    }

    #[test]
    fn bio_cid_num_n06_signed_inf_and_plus_rejection() {
        assert_eq!(
            parse_fast_float_prefix(b"-inf"),
            Some((f64::NEG_INFINITY, 4))
        );
        assert_eq!(
            parse_fast_float_prefix(b"-INFINITY"),
            Some((f64::NEG_INFINITY, 9))
        );
        // '+' is forbidden in parse_infnan under the default format: the
        // token fails the special branches and the decimal grammar (the
        // Gemmi wrapper strips at most one '+' before this layer).
        assert_eq!(parse_fast_float_prefix(b"+inf"), None);
        assert_eq!(parse_fast_float_prefix(b"+1"), None);
        // Incomplete spellings are not specials.
        assert_eq!(parse_fast_float_prefix(b"in"), None);
        assert_eq!(parse_fast_float_prefix(b"i"), None);
    }
}

#[cfg(test)]
mod bio_cid_num_n07_tests {
    use super::parse_fast_float_prefix;

    fn nan_with_consumed(token: &[u8]) -> Option<usize> {
        match parse_fast_float_prefix(token) {
            Some((value, consumed)) if value.is_nan() => Some(consumed),
            _ => None,
        }
    }

    // Expectations derived from pinned parse_infnan (fast_float.h:4483-4500):
    // the payload scan accepts [A-Za-z0-9_] ONLY; a ')' consumes through
    // itself, any forbidden character stops the scan with ptr at +3 while
    // the NaN stays assigned with success status.
    #[test]
    fn bio_cid_num_n07_payload_classes_and_delimiters() {
        // All case classes, plain and signed.
        for token in [&b"nan"[..], b"NaN", b"NAN", b"nAn"] {
            assert_eq!(nan_with_consumed(token), Some(3), "{token:?}");
        }
        assert_eq!(nan_with_consumed(b"-nan"), Some(4));
        assert_eq!(nan_with_consumed(b"-NaN"), Some(4));
        // Empty payload closes immediately.
        assert_eq!(nan_with_consumed(b"nan()"), Some(5));
        // Valid payloads of every allowed class, closed at various
        // positions.
        assert_eq!(nan_with_consumed(b"nan(a)"), Some(6));
        assert_eq!(nan_with_consumed(b"nan(A9_)"), Some(8));
        assert_eq!(nan_with_consumed(b"nan(ind)"), Some(8));
        assert_eq!(nan_with_consumed(b"nan(snan)"), Some(9));
        // Unclosed payload: ptr stays at +3, NaN still assigned.
        assert_eq!(nan_with_consumed(b"nan(a"), Some(3));
        assert_eq!(nan_with_consumed(b"nan(a_9"), Some(3));
        assert_eq!(nan_with_consumed(b"nan("), Some(3));
        // Forbidden payload character: same +3 stop.
        assert_eq!(nan_with_consumed(b"nan(a!)"), Some(3));
        assert_eq!(nan_with_consumed(b"nan(a-b)"), Some(3));
        assert_eq!(nan_with_consumed(b"nan(.a)"), Some(3));
        // The FIRST ')' closes the payload (the source loop returns at
        // the first closing delimiter): "nan(ab)c)" consumes through it;
        // later ')' bytes are ordinary trailing text.
        assert_eq!(nan_with_consumed(b"nan(ab)c)"), Some(7));
    }

    #[test]
    fn bio_cid_num_n07_suffixes_and_sign_rules() {
        // Text after a closed payload stays unconsumed.
        assert_eq!(nan_with_consumed(b"nan(a)x"), Some(6));
        assert_eq!(nan_with_consumed(b"nan())"), Some(5));
        // '+' is forbidden by the same default-format sign rule as inf.
        assert_eq!(parse_fast_float_prefix(b"+nan"), None);
        // Incomplete spelling is not special.
        assert_eq!(parse_fast_float_prefix(b"na"), None);
        assert_eq!(parse_fast_float_prefix(b"n"), None);
    }
}

/// Test-only expression of the parse_atom_inequality numeric-boundary
/// DECISIONS (BIO-CID-NUM N11) against the canonical typed outcome —
/// this is NOT a CID parser and is never exposed beyond tests. Source:
/// select.cpp:129-141: `fast_from_chars(cid.c_str() + pos, r.value)` (the
/// NUL-terminated overload), `result.ec != std::errc()` rejects with
/// " (expected number)", `pos = result.ptr - cid.c_str()`, then trailing
/// spaces are skipped before the end check.
/// Test-only expression of the parse_atom_inequality NUMERIC TAIL
/// decisions (BIO-NUM-EVID Step 4) against the canonical typed outcome
/// — this is NOT a CID parser and is never exposed beyond tests. The
/// property letter, surrounding spaces, and relation dispatch
/// (select.cpp:119-128) happen BEFORE this tail; their complete dispatch
/// is future unit C08, not this helper. Modeled decisions, source
/// select.cpp:129-140: `fast_from_chars(cid.c_str() + pos, r.value)`
/// (the NUL-terminated overload), `result.ec != std::errc()` rejects
/// with " (expected number)" — including result_out_of_range despite
/// assignment — then `pos = result.ptr - cid.c_str()`, then only ' '
/// bytes are skipped, and finally `pos != end` rejects.
#[cfg(test)]
fn inequality_decision(field: &[u8], field_end: usize) -> Result<(f64, usize, usize), ()> {
    let outcome = fast_from_chars_cstring_typed(field, 0.0);
    if outcome.category.is_some() {
        // Any nonzero ec — including result_out_of_range, which the
        // source rejects even though the value was assigned.
        return Err(());
    }
    let number_end = outcome.consumed;
    let mut pos = number_end;
    while field.get(pos) == Some(&b' ') {
        pos += 1;
    }
    if pos != field_end {
        // select.cpp:139-140: `if (pos != end) wrong_syntax(cid, pos);`
        return Err(());
    }
    Ok((outcome.value, number_end, pos))
}

#[cfg(test)]
mod bio_cid_num_n12_tests {
    use super::{FastFromCharsCategory, fast_float_from_chars_typed};
    use crate::cif::format_cif_f64;

    // Owner-map regressions (N12): the parse owner and the %.9g
    // serialization owner round-trip the canonical profile, and the
    // documented boundary counterexamples stay pinned.
    #[test]
    fn bio_cid_num_n12_owner_map_roundtrip() {
        for token in ["1.5", "-2.25", "42", "150"] {
            let outcome = fast_float_from_chars_typed(token.as_bytes(), 0.0);
            assert_eq!(outcome.category, None);
            assert_eq!(format_cif_f64(outcome.value), token);
        }
        // Extremes cross the owners with the pinned spellings.
        let inf = fast_float_from_chars_typed(b"inf", 0.0);
        assert_eq!(format_cif_f64(inf.value), "Inf");
        let over = fast_float_from_chars_typed(b"1.8e308", 0.0);
        assert_eq!(over.category, Some(FastFromCharsCategory::ResultOutOfRange));
        assert_eq!(format_cif_f64(over.value), "Inf");
    }

    #[test]
    fn bio_cid_num_n12_preserved_counterexamples() {
        // Hex-looking prefix consumes only the leading decimal digits.
        assert_eq!(fast_float_from_chars_typed(b"0x1A", 0.0).consumed, 1);
        // Locale comma is an ordinary terminator.
        assert_eq!(fast_float_from_chars_typed(b"1,5", 0.0).consumed, 1);
        // '+' inside specials is rejected by the default format.
        assert_eq!(
            fast_float_from_chars_typed(b"+inf", 0.0).category,
            Some(FastFromCharsCategory::InvalidArgument)
        );
        // d/D exponents belong to non-default format flags only.
        assert_eq!(fast_float_from_chars_typed(b"1d5", 0.0).consumed, 1);
        // The one-ulp rounding neighborhood stays pinned. The
        // 159-ending literal is itself out of range for the compiler, so
        // it reaches the formatter through the parse owner (which is the
        // pinned path anyway).
        assert_eq!(format_cif_f64(1.7976931348623158e308), "1.79769313e+308");
        let crosses = fast_float_from_chars_typed(b"1.7976931348623159e308", 0.0);
        assert_eq!(format_cif_f64(crosses.value), "Inf");
    }
}

#[cfg(test)]
mod bio_cid_num_n11_class_matrix_tests {
    use super::inequality_decision;

    // BIO-NUM-EVID-MATRIX Step 2: fixed class-by-ending PRODUCT with
    // ZERO case-selection discretion — ten literal tokens x five literal
    // suffixes at field_end = whole length (50 rows), plus ten
    // semicolon-boundary rows at field_end = token length. Expectations
    // derive from the literal table and select.cpp:129-140 (ec gate, ptr
    // offset, space-only skip, pos == end), never from the
    // implementation. Successes are exactly the first six tokens with
    // the empty or three-space suffix; every other whole-length
    // combination is wrong_syntax. The same six succeed at the
    // semicolon boundary; 1e999/1e-400 (range) and x/empty (invalid)
    // reject under every ending.
    #[test]
    fn bio_cid_num_n11_class_ending_matrix() {
        const NEG_ZERO_BITS: u64 = 0x8000_0000_0000_0000;
        let inf = f64::INFINITY;

        // Whole-end success rows: (token, suffix, value, consumed).
        // nan is asserted via is_nan per the source-assigned quiet NaN.
        let successes: &[(&[u8], &[u8], f64, usize)] = &[
            (b"1.5", b"", 1.5, 3),
            (b"1.5", b"   ", 1.5, 3),
            (b"-0", b"", f64::from_bits(NEG_ZERO_BITS), 2),
            (b"-0", b"   ", f64::from_bits(NEG_ZERO_BITS), 2),
            (b"inf", b"", inf, 3),
            (b"inf", b"   ", inf, 3),
            (b"-Infinity", b"", -inf, 9),
            (b"-Infinity", b"   ", -inf, 9),
            (b"nan", b"", f64::NAN, 3),
            (b"nan", b"   ", f64::NAN, 3),
            (b"0e999", b"", 0.0, 5),
            (b"0e999", b"   ", 0.0, 5),
        ];
        for &(token, suffix, value, consumed) in successes {
            let field: Vec<u8> = token
                .iter()
                .copied()
                .chain(suffix.iter().copied())
                .collect();
            let whole = field.len();
            let result = inequality_decision(&field, whole);
            let Ok((got_value, got_end, got_pos)) = result else {
                panic!("expected success for {token:?}+{suffix:?}");
            };
            if token == b"nan" {
                assert!(got_value.is_nan(), "nan row must stay NaN");
            } else if got_value == 0.0 {
                assert_eq!(
                    got_value.to_bits(),
                    value.to_bits(),
                    "zero sign for {token:?}"
                );
            } else {
                assert_eq!(got_value, value, "value for {token:?}+{suffix:?}");
            }
            assert_eq!(
                got_end, consumed,
                "consumed offset for {token:?}+{suffix:?}"
            );
            assert_eq!(got_pos, whole, "final position must equal the field end");
        }

        // All 50 whole-end rows: everything not in the success table
        // rejects. Range and invalid classes reject under every ending;
        // the six good tokens reject under tab/x/;q>2 endings because
        // only ' ' is skipped before the pos == end check.
        let tokens: [&[u8]; 10] = [
            b"1.5",
            b"-0",
            b"inf",
            b"-Infinity",
            b"nan",
            b"0e999",
            b"1e999",
            b"1e-400",
            b"x",
            b"",
        ];
        let suffixes: [&[u8]; 5] = [b"", b"   ", b"\t", b"x", b";q>2"];
        for &token in &tokens {
            for &suffix in &suffixes {
                let is_success = successes
                    .iter()
                    .any(|&(t, s, _, _)| t == token && s == suffix);
                if is_success {
                    continue;
                }
                let field: Vec<u8> = token
                    .iter()
                    .copied()
                    .chain(suffix.iter().copied())
                    .collect();
                assert_eq!(
                    inequality_decision(&field, field.len()),
                    Err(()),
                    "whole-end row {token:?}+{suffix:?} must be wrong_syntax"
                );
            }
        }

        // Ten semicolon-boundary rows: token + ";q>2" with field_end =
        // token length; only the first six classes reach success.
        for &(token, suffix, value, consumed) in successes {
            debug_assert!(suffix == b"" || suffix == b"   ");
            let field: Vec<u8> = token
                .iter()
                .copied()
                .chain(b";q>2".iter().copied())
                .collect();
            let result = inequality_decision(&field, token.len());
            let Ok((got_value, got_end, got_pos)) = result else {
                panic!("expected semicolon-boundary success for {token:?}");
            };
            if token == b"nan" {
                assert!(got_value.is_nan());
            } else if got_value == 0.0 {
                assert_eq!(got_value.to_bits(), value.to_bits());
            } else {
                assert_eq!(got_value, value);
            }
            assert_eq!(got_end, consumed);
            assert_eq!(got_pos, token.len());
        }
        for token in [b"1e999" as &[u8], b"1e-400", b"x", b""] {
            let field: Vec<u8> = token
                .iter()
                .copied()
                .chain(b";q>2".iter().copied())
                .collect();
            let end = token.len();
            assert_eq!(
                inequality_decision(&field, end),
                Err(()),
                "semicolon-boundary row {token:?} must be wrong_syntax"
            );
        }
    }
}

#[cfg(test)]
mod bio_cid_num_n11_tests {
    use super::inequality_decision;

    // Original inputs from BIO-CID-NUM N11, now with the explicit
    // source field end (BIO-NUM-EVID Step 4). The former unused
    // q/b/relation loop variables are removed — actual property and
    // relation dispatch is future unit C08; expectations here derive
    // only from the numeric tail, select.cpp:129-140. Renamed from
    // bio_cid_num_n11_normal_values_and_relations for the same reason.
    #[test]
    fn bio_cid_num_n11_values_and_end_equality() {
        // Normal value, exact end, no trailing text.
        assert_eq!(inequality_decision(b"1.5", 3), Ok((1.5, 3, 3)));
        // Integer and exponent forms.
        assert_eq!(inequality_decision(b"42", 2), Ok((42.0, 2, 2)));
        assert_eq!(inequality_decision(b"2e2", 3), Ok((200.0, 3, 3)));
        // Trailing spaces are skipped after the number: the
        // number end stays exact, the final position moves past
        // the spaces and equals the field end.
        assert_eq!(inequality_decision(b"1.5   ", 6), Ok((1.5, 3, 6)));
        assert_eq!(inequality_decision(b"1.5 x", 4), Ok((1.5, 3, 4)));
        assert_eq!(inequality_decision(b"1.5 x", 5), Err(()));
        // Tabs are NOT the source's trailing-space rule: the scan stops
        // at the number end, which must BE the field end.
        assert_eq!(inequality_decision(b"1.5\t", 3), Ok((1.5, 3, 3)));
        assert_eq!(inequality_decision(b"1.5\t", 4), Err(()));
        // A field that ends mid-number or inside trailing text is
        // wrong_syntax via the same pos != end gate.
        assert_eq!(inequality_decision(b"1.5", 2), Err(()));
        assert_eq!(inequality_decision(b"1.5;q>2", 5), Err(()));
    }

    #[test]
    fn bio_cid_num_n11_specials_range_and_errors() {
        // Literal inf/nan parse successfully (ec == 0) — accepted values
        // with their exact end offsets.
        assert_eq!(inequality_decision(b"inf", 3), Ok((f64::INFINITY, 3, 3)));
        assert_eq!(inequality_decision(b"inf ", 4), Ok((f64::INFINITY, 3, 4)));
        assert_eq!(inequality_decision(b"infx", 4), Err(()));
        assert_eq!(
            inequality_decision(b"-Infinity ", 10),
            Ok((f64::NEG_INFINITY, 9, 10))
        );
        let Ok((value, end, _)) = inequality_decision(b"nan", 3) else {
            panic!("nan must be an accepted value");
        };
        assert!(value.is_nan());
        assert_eq!(end, 3);
        // Range errors REJECT despite assignment: result.ec != 0 is the
        // source's sole gate.
        assert_eq!(inequality_decision(b"1e999", 5), Err(()));
        assert_eq!(inequality_decision(b"1e999", 3), Err(()));
        assert_eq!(inequality_decision(b"1e-400", 6), Err(()));
        // Zero-mantissa extremes succeed (no range status).
        assert_eq!(inequality_decision(b"0e999", 5), Ok((0.0, 5, 5)));
        // Lexical errors reject with the same gate.
        assert_eq!(inequality_decision(b"x", 1), Err(()));
        assert_eq!(inequality_decision(b"", 0), Err(()));
        // Leading whitespace is SKIPPED by the NUL-terminated overload
        // (fast_from_chars strips is_space bytes first), so the value
        // parses with the offset advanced past the spaces.
        assert_eq!(inequality_decision(b" 1.5", 4), Ok((1.5, 4, 4)));
        // Semicolon terminates a CID field: the number ends before it.
        assert_eq!(inequality_decision(b"1.5;q>2", 3), Ok((1.5, 3, 3)));
        // Incomplete exponent inside the field still succeeds as a
        // short number with the marker left over.
        assert_eq!(inequality_decision(b"1e;q>2", 1), Ok((1.0, 1, 1)));
    }
}

#[cfg(test)]
mod bio_cid_num_n08_tests {
    use super::{FastFromCharsCategory, fast_float_from_chars_typed};

    // Boundary-neighborhood regressions (N08) around the range-reporting
    // edges; expectations derived from IEEE-754 round-to-nearest and the
    // pinned assignment/status branches, not from the implementation.
    #[test]
    fn bio_cid_num_n08_max_finite_neighborhoods() {
        // DBL_MAX exactly: success.
        let max = fast_float_from_chars_typed(b"1.7976931348623157e308", 0.0);
        assert_eq!(
            (max.category, max.assigned, max.value),
            (None, true, f64::MAX)
        );
        // Neighborhood: the 158-ending decimal still rounds to DBL_MAX
        // (success); the 159-ending decimal crosses the round-to-nearest
        // midpoint and saturates to infinity with the range status.
        let near = fast_float_from_chars_typed(b"1.7976931348623158e308", 0.0);
        assert_eq!(
            (near.category, near.assigned, near.value),
            (None, true, f64::MAX)
        );
        let crosses = fast_float_from_chars_typed(b"1.7976931348623159e308", 0.0);
        assert_eq!(
            crosses.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert_eq!(crosses.value, f64::INFINITY);
        // A clearly-past decimal overflows identically.
        let over = fast_float_from_chars_typed(b"1.8e308", 0.0);
        assert_eq!(over.category, Some(FastFromCharsCategory::ResultOutOfRange));
        assert_eq!(over.value, f64::INFINITY);
        // Negative side mirrors both.
        let neg_over = fast_float_from_chars_typed(b"-1.8e308", 0.0);
        assert_eq!(
            neg_over.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert_eq!(neg_over.value, f64::NEG_INFINITY);
    }

    #[test]
    fn bio_cid_num_n08_min_normal_subnormal_round_to_zero() {
        // Min normal exactly and one decimal below it (still normal).
        let min_normal = fast_float_from_chars_typed(b"2.2250738585072014e-308", 0.0);
        assert_eq!(
            (min_normal.category, min_normal.value),
            (None, f64::MIN_POSITIVE)
        );
        let just_below = fast_float_from_chars_typed(b"2.225073858507201e-308", 0.0);
        assert_eq!(just_below.category, None);
        assert!(just_below.value > 0.0 && just_below.value < f64::MIN_POSITIVE);
        // Smallest subnormal: success; one decimal below it rounds to
        // zero with the range status (nonzero mantissa).
        let min_sub = fast_float_from_chars_typed(b"5e-324", 0.0);
        assert_eq!((min_sub.category, min_sub.value), (None, f64::from_bits(1)));
        let rounds_zero = fast_float_from_chars_typed(b"2e-324", 0.0);
        assert_eq!(
            rounds_zero.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert!(rounds_zero.assigned);
        assert_eq!(rounds_zero.value, 0.0);
        assert!(!rounds_zero.value.is_sign_negative());
    }

    #[test]
    fn bio_cid_num_n08_signed_zeros_and_exponent_extremes() {
        // Signed zeros are exact successes in both directions.
        let pos_zero = fast_float_from_chars_typed(b"0.0", 0.0);
        assert_eq!((pos_zero.category, pos_zero.value), (None, 0.0));
        assert!(!pos_zero.value.is_sign_negative());
        let neg_zero = fast_float_from_chars_typed(b"-0.0", 0.0);
        assert_eq!(neg_zero.category, None);
        assert_eq!(neg_zero.value, 0.0);
        assert!(neg_zero.value.is_sign_negative());
        // Exponent extremes with zero mantissa stay successes; with a
        // nonzero mantissa they carry the range status and assigned
        // extremes.
        let zero_big = fast_float_from_chars_typed(b"0.0e999999", 0.0);
        assert_eq!((zero_big.category, zero_big.value), (None, 0.0));
        let zero_small = fast_float_from_chars_typed(b"-0.0e-999999", 0.0);
        assert_eq!(zero_small.category, None);
        assert!(zero_small.value.is_sign_negative());
        let one_big = fast_float_from_chars_typed(b"1.0e999999", 0.0);
        assert_eq!(
            one_big.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert_eq!(one_big.value, f64::INFINITY);
        let one_small = fast_float_from_chars_typed(b"-1.0e-999999", 0.0);
        assert_eq!(
            one_small.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert_eq!(one_small.value, 0.0);
        assert!(one_small.value.is_sign_negative());
    }
}

#[cfg(test)]
mod bio_cid_num_rangefix_tests {
    use super::{FastFromCharsCategory, fast_float_from_chars_typed};

    // BIO-NUM-RANGE fixed matrix: expectations derived from pinned
    // fast_float.h:4728-4734 (to_float BEFORE the range test) and the
    // mantissa-zero distinction. BIO-NUM-RANGE-MATRIX completion: every
    // input runs across the finite/negative-zero/NaN-payload prior set
    // with consumed/assigned/category and exact meaningful bits asserted
    // per row; every original case is preserved.
    #[test]
    fn bio_cid_num_rangefix_signed_overflow_assigns_infinity() {
        for (token, expected) in [
            (&b"1e999"[..], f64::INFINITY),
            (b"-1e999", f64::NEG_INFINITY),
            (b"2e308", f64::INFINITY),
            (b"-2.5e400", f64::NEG_INFINITY),
            (b"1e400000000000000000000", f64::INFINITY),
        ] {
            for (prior_name, prior) in super::matrix_priors() {
                let outcome = fast_float_from_chars_typed(token, prior);
                assert_eq!(
                    outcome.category,
                    Some(FastFromCharsCategory::ResultOutOfRange),
                    "{token:?} prior {prior_name}"
                );
                assert!(outcome.assigned, "{token:?} prior {prior_name}");
                assert_eq!(
                    outcome.value.to_bits(),
                    expected.to_bits(),
                    "{token:?} prior {prior_name}"
                );
                assert_eq!(
                    outcome.consumed,
                    token.len(),
                    "{token:?} prior {prior_name}"
                );
            }
        }
    }

    #[test]
    fn bio_cid_num_rangefix_nonzero_underflow_assigns_signed_zero() {
        for (token, negative, consumed) in [
            (&b"1e-400"[..], false, 6usize),
            (b"-1e-400", true, 7),
            (b"3.7e-325", false, 8),
            (b"-9e-999999", true, 10),
        ] {
            for (prior_name, prior) in super::matrix_priors() {
                let outcome = fast_float_from_chars_typed(token, prior);
                assert_eq!(
                    outcome.category,
                    Some(FastFromCharsCategory::ResultOutOfRange),
                    "{token:?} prior {prior_name}"
                );
                assert!(outcome.assigned, "{token:?} prior {prior_name}");
                assert_eq!(outcome.value, 0.0, "{token:?} prior {prior_name}");
                assert_eq!(
                    outcome.value.is_sign_negative(),
                    negative,
                    "{token:?} prior {prior_name}"
                );
                // A negative-zero prior legitimately bit-matches an
                // assigned -0.0; assignment is proven by the flags above.
                assert_eq!(outcome.consumed, consumed, "{token:?} prior {prior_name}");
            }
        }
    }

    #[test]
    fn bio_cid_num_rangefix_zero_mantissa_and_literals() {
        // Zero mantissa with any exponent: exact zero, ordinary success.
        for token in [&b"0e999"[..], b"-0e999", b"0.000e-999", b"-0.0", b"0"] {
            for (prior_name, prior) in super::matrix_priors() {
                let outcome = fast_float_from_chars_typed(token, prior);
                assert_eq!(outcome.category, None, "{token:?} prior {prior_name}");
                assert!(outcome.assigned, "{token:?} prior {prior_name}");
                assert_eq!(outcome.value, 0.0, "{token:?} prior {prior_name}");
                assert_eq!(
                    outcome.consumed,
                    token.len(),
                    "{token:?} prior {prior_name}"
                );
            }
            // Signed-zero detail asserted once per token.
            let value = fast_float_from_chars_typed(token, 1.0).value;
            assert_eq!(value.is_sign_negative(), token[0] == b'-', "{token:?}");
        }
        // Literal inf/nan: parse_infnan successes, never range errors,
        // across every prior with exact offsets.
        for (token, consumed) in [
            (&b"inf"[..], 3usize),
            (b"-Infinity", 9),
            (b"nan", 3),
            (b"-NAN(payload)", 13),
        ] {
            for (prior_name, prior) in super::matrix_priors() {
                let outcome = fast_float_from_chars_typed(token, prior);
                assert_eq!(outcome.category, None, "{token:?} prior {prior_name}");
                assert!(outcome.assigned, "{token:?} prior {prior_name}");
                assert!(
                    outcome.value.is_nan() || outcome.value.is_infinite(),
                    "{token:?} prior {prior_name}"
                );
                assert_eq!(outcome.consumed, consumed, "{token:?} prior {prior_name}");
            }
        }
    }

    #[test]
    fn bio_cid_num_rangefix_invalid_input_retains_exact_prior_bits() {
        // Only a failed lex retains the prior: every prior pattern is
        // compared by exact bits — including the NaN payload — never by
        // is_nan alone.
        let nan_payload_bits = 0x7ff8_1234_5678_9abcu64;
        for (token, consumed) in [
            (&b"zz"[..], 0usize),
            (b"x", 0),
            (b"", 0),
            (b" ", 0),
            (b"e5", 0),
        ] {
            for (prior_name, prior) in super::matrix_priors() {
                let outcome = fast_float_from_chars_typed(token, prior);
                assert_eq!(
                    outcome.category,
                    Some(FastFromCharsCategory::InvalidArgument),
                    "{token:?} prior {prior_name}"
                );
                assert!(!outcome.assigned, "{token:?} prior {prior_name}");
                assert_eq!(
                    outcome.value.to_bits(),
                    prior.to_bits(),
                    "{token:?} prior {prior_name}"
                );
                assert_eq!(outcome.consumed, consumed, "{token:?} prior {prior_name}");
            }
            // The NaN-payload prior specifically: every payload bit
            // retained, not merely nan-ness.
            let outcome = fast_float_from_chars_typed(token, f64::from_bits(nan_payload_bits));
            assert_eq!(outcome.value.to_bits(), nan_payload_bits, "{token:?}");
        }
    }

    #[test]
    fn bio_cid_num_rangefix_finite_and_max_boundary_success() {
        // DBL_MAX itself parses exactly and stays a success.
        let max = fast_float_from_chars_typed(b"1.7976931348623157e308", 0.0);
        assert_eq!(max.category, None);
        assert!(max.assigned);
        assert_eq!(max.value, f64::MAX);
        assert_eq!(max.consumed, 22);
        // Min normal and min subnormal boundaries: ordinary successes.
        let min_normal = fast_float_from_chars_typed(b"2.2250738585072014e-308", 0.0);
        assert_eq!(min_normal.category, None);
        assert_eq!(min_normal.value, f64::MIN_POSITIVE);
        assert_eq!(min_normal.consumed, 23);
        let min_sub = fast_float_from_chars_typed(b"5e-324", 0.0);
        assert_eq!(min_sub.category, None);
        assert!(min_sub.assigned);
        assert_eq!(min_sub.value, f64::from_bits(1)); // smallest subnormal
        assert_eq!(min_sub.consumed, b"5e-324".len());
        // A decimal clearly past DBL_MAX overflows with assignment. (A
        // one-ulp-past decimal such as 1.7976931348623159e308 also
        // saturates to infinity with ResultOutOfRange under round-to-
        // nearest; only the 158-ending neighbor stays at f64::MAX — see
        // bio_cid_num_n08_tests for that neighborhood.)
        let over = fast_float_from_chars_typed(b"1.8e308", 0.0);
        assert_eq!(over.category, Some(FastFromCharsCategory::ResultOutOfRange));
        assert!(over.assigned);
        assert_eq!(over.value, f64::INFINITY);
        assert_eq!(over.consumed, 7);
    }
}
#[cfg(test)]
mod bio_cid_num_n04_tests {
    use super::{
        FastFromCharsCategory, fast_from_chars_bounded_typed, fast_from_chars_cstring_typed,
    };

    // Expectations derived from pinned atof.hpp:16-27 and the is_space
    // table (bytes 9-13 and 32), independently of the implementation.
    #[test]
    fn bio_cid_num_n04_whitespace_and_plus_skip_offsets() {
        // Every whitespace byte is skipped before one optional plus; the
        // consumed offset counts the ORIGINAL input.
        for ws in [9u8, 10, 11, 12, 13, 32] {
            let field = [ws, b'1', b'.', b'5'];
            let outcome = fast_from_chars_bounded_typed(&field, 0.0);
            assert_eq!(outcome.category, None, "ws {ws}");
            assert_eq!(outcome.value, 1.5, "ws {ws}");
            assert_eq!(outcome.consumed, 4, "ws {ws}");
        }
        // Plus after whitespace; minus belongs to the number.
        assert_eq!(fast_from_chars_bounded_typed(b"\t+2", 0.0).consumed, 3);
        assert_eq!(fast_from_chars_bounded_typed(b" -2", 0.0).value, -2.0);
        assert_eq!(fast_from_chars_bounded_typed(b" -2", 0.0).consumed, 3);
        // A second plus is NOT skipped: "++" converts nothing.
        let outcome = fast_from_chars_bounded_typed(b"++2", 0.0);
        assert_eq!(
            outcome.category,
            Some(FastFromCharsCategory::InvalidArgument)
        );
        assert_eq!(outcome.consumed, 1);
        assert!(!outcome.assigned);
        // Skip-only input: invalid argument with the post-skip offset.
        let blank = fast_from_chars_bounded_typed(b" \t ", 0.0);
        assert_eq!(blank.category, Some(FastFromCharsCategory::InvalidArgument));
        assert_eq!(blank.consumed, 3);
    }

    #[test]
    fn bio_cid_num_n04_nul_terminated_overload() {
        // The C-string overload stops at the first NUL: trailing bytes
        // after it never influence the parse or the offset.
        assert_eq!(fast_from_chars_cstring_typed(b"1.5\0 2.0", 0.0).value, 1.5);
        assert_eq!(fast_from_chars_cstring_typed(b"1.5\0 2.0", 0.0).consumed, 3);
        assert_eq!(fast_from_chars_cstring_typed(b"\0x", 0.0).consumed, 0);
        assert_eq!(fast_from_chars_cstring_typed(b"\0x", 0.0).assigned, false);
        // An embedded NUL ends the token mid-number; the bounded
        // overload by contrast treats NUL as an ordinary terminator
        // AFTER the grammar stops there anyway — both report offset 3
        // for "12\0" (digits 1,2 then NUL stops).
        assert_eq!(fast_from_chars_cstring_typed(b"12\0", 0.0).value, 12.0);
        assert_eq!(fast_from_chars_cstring_typed(b"12\0", 0.0).consumed, 2);
        // Empty suffix (field only whitespace then NUL).
        let outcome = fast_from_chars_cstring_typed(b" \0", 0.0);
        assert_eq!(
            outcome.category,
            Some(FastFromCharsCategory::InvalidArgument)
        );
        assert_eq!(outcome.consumed, 1);
    }
}

#[cfg(test)]
mod bio_cid_num_n03_tests {
    use super::{FastFromCharsCategory, fast_float_from_chars_typed};

    // Expectations derived from the pinned from_chars_result_t contract
    // (fast_float.h:204-212) and its conversion branches, independently
    // of the implementation under test.
    #[test]
    fn bio_cid_num_n03_empty_and_invalid_input() {
        // No number converts: ec = invalid_argument, ptr back at first,
        // value unassigned (caller prior retained).
        for input in [&b""[..], b"x", b" ", b"-", b".", b"e5"] {
            let outcome = fast_float_from_chars_typed(input, -7.25);
            assert_eq!(
                outcome.category,
                Some(FastFromCharsCategory::InvalidArgument)
            );
            assert_eq!(outcome.consumed, 0);
            assert!(!outcome.assigned);
            assert_eq!(outcome.value, -7.25);
        }
        // Prior retention is the caller's exact value, not a zeroed
        // substitute.
        let outcome = fast_float_from_chars_typed(b"zz", 42.5);
        assert_eq!(outcome.value, 42.5);
    }

    #[test]
    fn bio_cid_num_n03_success_assignment_and_offsets() {
        let outcome = fast_float_from_chars_typed(b"1.5e2x", 0.0);
        assert_eq!(outcome.category, None);
        assert!(outcome.assigned);
        assert_eq!(outcome.value, 150.0);
        assert_eq!(outcome.consumed, 5);
        // Literal inf/nan spellings are parse_infnan SUCCESSES with
        // assignment, distinct from decimal overflow.
        let inf = fast_float_from_chars_typed(b"-inf", 0.0);
        assert_eq!(inf.category, None);
        assert!(inf.assigned);
        assert_eq!(inf.value, f64::NEG_INFINITY);
        assert_eq!(inf.consumed, 4);
        let nan = fast_float_from_chars_typed(b"nan(x", 0.0);
        assert_eq!(nan.category, None);
        assert!(nan.assigned);
        assert!(nan.value.is_nan());
        assert_eq!(nan.consumed, 3);
    }

    #[test]
    fn bio_cid_num_n03_out_of_range_and_underflow() {
        // BIO-NUM-RANGE correction (fast_float.h:4728-4734): to_float
        // runs BEFORE the range test, so range errors still ASSIGN.
        // History: this test's first green version wrongly asserted
        // prior-retention on overflow and plain success on nonzero
        // underflow; those green-but-wrong expectations are corrected
        // here and retained in the receipt, not erased.
        let outcome = fast_float_from_chars_typed(b"1e999", -3.5);
        assert_eq!(
            outcome.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert!(outcome.assigned);
        assert_eq!(outcome.value, f64::INFINITY);
        assert_eq!(outcome.consumed, 5);
        let neg = fast_float_from_chars_typed(b"-1e999", 9.0);
        assert_eq!(neg.category, Some(FastFromCharsCategory::ResultOutOfRange));
        assert!(neg.assigned);
        assert_eq!(neg.value, f64::NEG_INFINITY);
        assert_eq!(neg.consumed, 6);
        let under = fast_float_from_chars_typed(b"1e-400", 0.0);
        assert_eq!(
            under.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert!(under.assigned);
        assert_eq!(under.value, 0.0);
        assert!(!under.value.is_sign_negative());
        let neg_under = fast_float_from_chars_typed(b"-1e-400", 0.0);
        assert_eq!(
            neg_under.category,
            Some(FastFromCharsCategory::ResultOutOfRange)
        );
        assert!(neg_under.assigned);
        assert_eq!(neg_under.value, 0.0);
        assert!(neg_under.value.is_sign_negative());
        // Zero mantissa with an extreme exponent stays an ordinary
        // success (pns.mantissa == 0 never trips the underflow branch).
        let zero_exp = fast_float_from_chars_typed(b"0e999", 0.0);
        assert_eq!(zero_exp.category, None);
        assert!(zero_exp.assigned);
        assert_eq!(zero_exp.value, 0.0);
        // Subnormal magnitudes stay ordinary successes.
        let sub = fast_float_from_chars_typed(b"5e-324", 0.0);
        assert_eq!(sub.category, None);
        assert!(sub.assigned);
        assert!(sub.value > 0.0);
    }
}

#[cfg(test)]
mod bio_cid_num_n02_tests {
    use super::{gemmi_fast_atof_with_end, parse_fast_float_prefix};

    // Expectations derived from pinned atof.hpp:14-19 and pdb.cpp
    // read_double consumers, independently of the implementation.
    #[test]
    fn bio_cid_num_n02_prefix_identity_after_factoring() {
        // The canonical owner keeps the exact pre-factoring prefix
        // semantics (fast_float.h decimal + parse_infnan branches).
        assert_eq!(parse_fast_float_prefix(b"1.5e2x"), Some((150.0, 5)));
        assert_eq!(parse_fast_float_prefix(b"-.25"), Some((-0.25, 4)));
        assert_eq!(parse_fast_float_prefix(b"1e+"), Some((1.0, 1)));
        assert_eq!(parse_fast_float_prefix(b"inf"), Some((f64::INFINITY, 3)));
        assert!(matches!(parse_fast_float_prefix(b"nan(x"), Some((v, 3)) if v.is_nan()));
        assert_eq!(parse_fast_float_prefix(b"x"), None);
    }

    #[test]
    fn bio_cid_num_n02_read_double_field_semantics() {
        // gemmi_fast_atof_with_end over fixed fields: NUL stops the field,
        // leading is_space is skipped, one optional '+' consumed,
        // non-conversion keeps +0.0 with the post-skip offset.
        assert_eq!(gemmi_fast_atof_with_end(b" 1.5\0 2.0"), (1.5, 4));
        assert_eq!(gemmi_fast_atof_with_end(b"\t+2.25x"), (2.25, 6));
        assert_eq!(gemmi_fast_atof_with_end(b"  \0junk"), (0.0, 2));
        assert_eq!(gemmi_fast_atof_with_end(b"-3e2 "), (-300.0, 4));
        assert_eq!(gemmi_fast_atof_with_end(b""), (0.0, 0));
    }
}
