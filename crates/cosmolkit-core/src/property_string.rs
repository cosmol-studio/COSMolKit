//! Source-backed string projection for detached atom and bond properties.

use cosmolkit_model::{PropertyValue, PropertyValueKind};
use std::fmt::{self, Write as _};
use thiserror::Error;

const DOUBLE_SIGNIFICANT_DIGITS: usize = 17;
const DOUBLE_FRACTIONAL_DIGITS: usize = DOUBLE_SIGNIFICANT_DIGITS - 1;
// One digit, decimal point, sixteen fractional digits, `e`, exponent sign,
// and the three digits needed for binary64's -324..=308 decimal exponents.
const SCIENTIFIC_BUFFER_LEN: usize = 1 + 1 + DOUBLE_FRACTIONAL_DIGITS + 1 + 1 + 3;
// The formatted result adds at most the original value's sign.
const DOUBLE_OUTPUT_CAPACITY: usize = SCIENTIFIC_BUFFER_LEN + 1;

struct ScientificBuffer {
    bytes: [u8; SCIENTIFIC_BUFFER_LEN],
    len: usize,
}

impl ScientificBuffer {
    const fn new() -> Self {
        Self {
            bytes: [0; SCIENTIFIC_BUFFER_LEN],
            len: 0,
        }
    }

    fn as_bytes(&self) -> &[u8] {
        &self.bytes[..self.len]
    }
}

impl fmt::Write for ScientificBuffer {
    fn write_str(&mut self, value: &str) -> fmt::Result {
        let Some(end) = self.len.checked_add(value.len()) else {
            return Err(fmt::Error);
        };
        if end > self.bytes.len() {
            return Err(fmt::Error);
        }
        self.bytes[self.len..end].copy_from_slice(value.as_bytes());
        self.len = end;
        Ok(())
    }
}

fn parse_scientific_exponent(bytes: &[u8]) -> i16 {
    let (negative, digits) = match bytes.first() {
        Some(b'-') => (true, &bytes[1..]),
        Some(b'+') => (false, &bytes[1..]),
        _ => (false, bytes),
    };
    assert!(
        !digits.is_empty() && digits.len() <= 3,
        "binary64 scientific exponent has one to three digits"
    );

    let mut magnitude = 0_i16;
    for digit in digits {
        assert!(digit.is_ascii_digit(), "scientific exponent is decimal");
        magnitude = magnitude * 10 + i16::from(*digit - b'0');
    }
    if negative { -magnitude } else { magnitude }
}

fn push_ascii_digits(output: &mut String, digits: &[u8]) {
    output.push_str(std::str::from_utf8(digits).expect("formatter emits ASCII decimal digits"));
}

fn push_scientific_exponent(output: &mut String, exponent: i16) {
    output.push('e');
    output.push(if exponent < 0 { '-' } else { '+' });

    let magnitude = exponent.unsigned_abs();
    if magnitude >= 100 {
        output.push(char::from(b'0' + (magnitude / 100) as u8));
    }
    output.push(char::from(b'0' + ((magnitude / 10) % 10) as u8));
    output.push(char::from(b'0' + (magnitude % 10) as u8));
}

fn format_boost_double(value: f64) -> String {
    // BEGIN BOOST CPP FUNCTION get_inf_nan_impl
    // Boost✔️✔️: if (boost::core::isnan(value)) {
    // Boost✔️✔️:     if (boost::core::signbit(value)) {
    // Boost✔️✔️:         return lc_minus_nan;
    // Boost✔️✔️:     }
    // Boost✔️✔️:     return lc_nan;
    // Boost✔️✔️: } else if (boost::core::isinf(value)) {
    // Boost✔️✔️:     if (boost::core::signbit(value)) {
    // Boost✔️✔️:         return lc_minus_infinity;
    // Boost✔️✔️:     }
    // Boost✔️✔️:     return lc_infinity;
    // Boost✔️✔️: }
    // Boost✔️✔️: return nullptr;
    // END BOOST CPP FUNCTION get_inf_nan_impl
    let negative = value.is_sign_negative();
    if value.is_nan() {
        return if negative { "-nan" } else { "nan" }.to_owned();
    }
    if value.is_infinite() {
        return if negative { "-inf" } else { "inf" }.to_owned();
    }
    if value == 0.0 {
        return if negative { "-0" } else { "0" }.to_owned();
    }

    // BEGIN BOOST CPP FUNCTION shl_real_type(double, char*)
    // Boost✔️✔️: finish = start +
    // Boost✔️✔️:     boost::core::snprintf(begin, CharacterBufferSize,
    // Boost✔️✔️:     "%.*g", static_cast<int>(boost::detail::lcast_get_precision<double>()), val);
    // Boost✔️✔️: return finish > start;
    // END BOOST CPP FUNCTION shl_real_type(double, char*)
    // Approved independent system-formatter boundary: Rust 1.98.1 LowerExp
    // performs the sole numeric conversion and ties-to-even rounding. The
    // code below only parses those rounded digits and implements `%g` layout;
    // it does not translate the glibc formatter implementation.
    let mut scientific = ScientificBuffer::new();
    write!(&mut scientific, "{:.16e}", value.abs())
        .expect("binary64 precision-17 scientific form fits the proven stack bound");
    let scientific = scientific.as_bytes();
    let exponent_marker = scientific
        .iter()
        .position(|byte| *byte == b'e')
        .expect("LowerExp emits an exponent marker");
    assert_eq!(
        exponent_marker,
        1 + 1 + DOUBLE_FRACTIONAL_DIGITS,
        "precision-16 LowerExp emits one integer and sixteen fractional digits"
    );
    assert_eq!(scientific[1], b'.', "LowerExp emits a decimal point");

    let mut digits = [0_u8; DOUBLE_SIGNIFICANT_DIGITS];
    digits[0] = scientific[0];
    digits[1..].copy_from_slice(&scientific[2..exponent_marker]);
    assert!(
        digits.iter().all(u8::is_ascii_digit) && digits[0] != b'0',
        "finite nonzero scientific mantissa has seventeen decimal digits"
    );
    let exponent = parse_scientific_exponent(&scientific[exponent_marker + 1..]);
    let significant_end = digits
        .iter()
        .rposition(|digit| *digit != b'0')
        .expect("finite nonzero mantissa has a nonzero digit")
        + 1;

    let mut output = String::with_capacity(DOUBLE_OUTPUT_CAPACITY);
    if negative {
        output.push('-');
    }

    if (-4..17).contains(&exponent) {
        let decimal_position = exponent + 1;
        if decimal_position <= 0 {
            output.push_str("0.");
            for _ in 0..-decimal_position {
                output.push('0');
            }
            push_ascii_digits(&mut output, &digits[..significant_end]);
        } else {
            let decimal_position = decimal_position as usize;
            push_ascii_digits(&mut output, &digits[..decimal_position]);
            if significant_end > decimal_position {
                output.push('.');
                push_ascii_digits(&mut output, &digits[decimal_position..significant_end]);
            }
        }
    } else {
        output.push(char::from(digits[0]));
        if significant_end > 1 {
            output.push('.');
            push_ascii_digits(&mut output, &digits[1..significant_end]);
        }
        push_scientific_exponent(&mut output, exponent);
    }

    // Behavior evidence: fixed owner regressions plus 1,102,984 native
    // byte-exact oracle cases passed with Rust 1.98.1. That diagnostic is not
    // exhaustive and does not establish WASM or future-toolchain behavior.
    // Complexity review: one bounded stack conversion, bounded digit
    // repositioning, and one output allocation match the source path's
    // constant binary64 work and allocation shape.
    output
}

/// Error returned when a modeled property kind has no verified source string
/// projection yet.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[non_exhaustive]
pub enum PropertyStringError {
    #[error("property value kind {kind:?} has no verified source string projection")]
    UnsupportedKind { kind: PropertyValueKind },
}

impl PropertyStringError {
    #[must_use]
    pub const fn kind(self) -> PropertyValueKind {
        match self {
            Self::UnsupportedKind { kind } => kind,
        }
    }
}

/// Convert a canonical detached property value using the pinned source's
/// modeled `RDValue` string-conversion behavior.
pub fn property_value_to_string(value: &PropertyValue) -> Result<String, PropertyStringError> {
    // BEGIN RDKIT CPP FUNCTION rdvalue_tostring
    // RDKit✔️✔️: switch (val.getTag()) {
    // RDKit✔️✔️:   case RDTypeTag::StringTag:
    // RDKit✔️✔️:     res = rdvalue_cast<std::string>(val);
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case RDTypeTag::IntTag:
    // RDKit✔️✔️:     res = boost::lexical_cast<std::string>(rdvalue_cast<int>(val));
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case RDTypeTag::DoubleTag: {
    // RDKit✔️✔️:     Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit✔️✔️:     res = boost::lexical_cast<std::string>(rdvalue_cast<double>(val));
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: #ifdef RDVALUE_HASBOOL
    // RDKit✔️✔️:   case RDTypeTag::BoolTag:
    // RDKit✔️✔️:     res = boost::lexical_cast<std::string>(rdvalue_cast<bool>(val));
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️: #endif
    // END RDKIT CPP FUNCTION rdvalue_tostring
    // Behavior review: String bytes are copied without inspection, signed
    // i32 values use normalized decimal spelling, and Boost's bool stream
    // spelling is reproduced as `1` or `0`. Double uses the approved pure-Rust
    // precision-17 `%g` implementation above; fixed owner regressions and the
    // frozen native diagnostic support the modeled byte-equivalence marker.
    // Complexity review: each scalar branch performs one output allocation
    // and a constant-size conversion, matching the source assignment into its
    // output std::string without an intermediate value cache.
    match value {
        PropertyValue::String(value) => Ok(value.clone()),
        PropertyValue::Int(value) => Ok(value.to_string()),
        PropertyValue::Double(value) => Ok(format_boost_double(*value)),
        PropertyValue::Bool(value) => Ok(if *value { "1" } else { "0" }.to_owned()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_double_cases(cases: &[(u64, &str)]) {
        for &(bits, expected) in cases {
            assert_eq!(
                property_value_to_string(&PropertyValue::Double(f64::from_bits(bits))),
                Ok(expected.to_owned()),
                "binary64 bits {bits:#018x}"
            );
        }
    }

    #[test]
    fn property_string_double_covers_signs_boundaries_extrema_and_specials() {
        // Fixed oracle bytes from RDKit 2026.03.1 / Boost 1.85. These cases
        // stay local and never invoke or regenerate the upstream reference.
        let cases = [
            (0x0000_0000_0000_0000, "0"),
            (0x8000_0000_0000_0000, "-0"),
            (0x3ff0_0000_0000_0000, "1"),
            (0xbff0_0000_0000_0000, "-1"),
            (0x3fb9_9999_9999_999a, "0.10000000000000001"),
            (0x3ee4_f8b5_88e3_68f0, "9.9999999999999991e-06"),
            (0x3ee4_f8b5_88e3_68f1, "1.0000000000000001e-05"),
            (0x3ee4_f8b5_88e3_68f2, "1.0000000000000003e-05"),
            (0x3f1a_36e2_eb1c_432c, "9.9999999999999991e-05"),
            (0x3f1a_36e2_eb1c_432d, "0.0001"),
            (0x3f1a_36e2_eb1c_432e, "0.00010000000000000002"),
            (0x4341_c379_37e0_8000, "10000000000000000"),
            (0x4376_3457_85d8_9fff, "99999999999999984"),
            (0x4376_3457_85d8_a000, "1e+17"),
            (0x4376_3457_85d8_a001, "1.0000000000000002e+17"),
            (0x0000_0000_0000_0001, "4.9406564584124654e-324"),
            (0x000f_ffff_ffff_ffff, "2.2250738585072009e-308"),
            (0x0010_0000_0000_0000, "2.2250738585072014e-308"),
            (0x7fef_ffff_ffff_ffff, "1.7976931348623157e+308"),
            (0x7ff0_0000_0000_0000, "inf"),
            (0xfff0_0000_0000_0000, "-inf"),
            (0x7ff8_0000_0000_0000, "nan"),
            (0xfff8_0000_0000_0000, "-nan"),
            (0x7ff0_0000_0000_0001, "nan"),
            (0xfff0_0000_0000_0001, "-nan"),
            (0x7ff8_1234_5678_9abc, "nan"),
            (0xfff8_1234_5678_9abc, "-nan"),
        ];
        assert_double_cases(&cases);
    }

    #[test]
    fn property_string_double_recovers_all_gpoint_counterexamples() {
        // All 35 mismatches from the frozen gpoint 0.3.0 diagnostic: signed
        // zero plus the signed predecessor cases at decimal-power boundaries.
        let cases = [
            (0x8000_0000_0000_0000, "-0"),
            (0x3f1a_36e2_eb1c_432c, "9.9999999999999991e-05"),
            (0xbf1a_36e2_eb1c_432c, "-9.9999999999999991e-05"),
            (0x3f84_7ae1_47ae_147a, "0.0099999999999999985"),
            (0xbf84_7ae1_47ae_147a, "-0.0099999999999999985"),
            (0x3fb9_9999_9999_9999, "0.099999999999999992"),
            (0xbfb9_9999_9999_9999, "-0.099999999999999992"),
            (0x4058_ffff_ffff_ffff, "99.999999999999986"),
            (0xc058_ffff_ffff_ffff, "-99.999999999999986"),
            (0x408f_3fff_ffff_ffff, "999.99999999999989"),
            (0xc08f_3fff_ffff_ffff, "-999.99999999999989"),
            (0x40c3_87ff_ffff_ffff, "9999.9999999999982"),
            (0xc0c3_87ff_ffff_ffff, "-9999.9999999999982"),
            (0x40f8_69ff_ffff_ffff, "99999.999999999985"),
            (0xc0f8_69ff_ffff_ffff, "-99999.999999999985"),
            (0x412e_847f_ffff_ffff, "999999.99999999988"),
            (0xc12e_847f_ffff_ffff, "-999999.99999999988"),
            (0x4163_12cf_ffff_ffff, "9999999.9999999981"),
            (0xc163_12cf_ffff_ffff, "-9999999.9999999981"),
            (0x4197_d783_ffff_ffff, "99999999.999999985"),
            (0xc197_d783_ffff_ffff, "-99999999.999999985"),
            (0x41cd_cd64_ffff_ffff, "999999999.99999988"),
            (0xc1cd_cd64_ffff_ffff, "-999999999.99999988"),
            (0x4202_a05f_1fff_ffff, "9999999999.9999981"),
            (0xc202_a05f_1fff_ffff, "-9999999999.9999981"),
            (0x4237_4876_e7ff_ffff, "99999999999.999985"),
            (0xc237_4876_e7ff_ffff, "-99999999999.999985"),
            (0x426d_1a94_a1ff_ffff, "999999999999.99988"),
            (0xc26d_1a94_a1ff_ffff, "-999999999999.99988"),
            (0x42d6_bcc4_1e8f_ffff, "99999999999999.984"),
            (0xc2d6_bcc4_1e8f_ffff, "-99999999999999.984"),
            (0x430c_6bf5_2633_ffff, "999999999999999.88"),
            (0xc30c_6bf5_2633_ffff, "-999999999999999.88"),
            (0x4376_3457_85d8_9fff, "99999999999999984"),
            (0xc376_3457_85d8_9fff, "-99999999999999984"),
        ];
        assert_eq!(cases.len(), 35);
        assert_double_cases(&cases);
    }

    #[test]
    fn property_string_scalar_preserves_string_bytes() {
        let bytes = "leading\0middle.\u{00e9}\ntrailing ".to_owned();
        assert_eq!(
            property_value_to_string(&PropertyValue::String(bytes.clone())),
            Ok(bytes)
        );
    }

    #[test]
    fn property_string_scalar_normalizes_signed_ints_and_bool_words() {
        let cases = [
            (i32::MIN, "-2147483648"),
            (-7, "-7"),
            (0, "0"),
            (7, "7"),
            (i32::MAX, "2147483647"),
        ];
        for (value, expected) in cases {
            assert_eq!(
                property_value_to_string(&PropertyValue::Int(value)),
                Ok(expected.to_owned())
            );
        }
        assert_eq!(
            property_value_to_string(&PropertyValue::Bool(false)),
            Ok("0".to_owned())
        );
        assert_eq!(
            property_value_to_string(&PropertyValue::Bool(true)),
            Ok("1".to_owned())
        );
    }

    #[test]
    fn property_string_scalar_retains_explicit_wrong_kind_access_errors() {
        let string_as_int = PropertyValue::String("007".to_owned())
            .as_int()
            .unwrap_err();
        assert_eq!(string_as_int.expected(), PropertyValueKind::Int);
        assert_eq!(string_as_int.actual(), PropertyValueKind::String);

        let int_as_string = PropertyValue::Int(7).as_string().unwrap_err();
        assert_eq!(int_as_string.expected(), PropertyValueKind::String);
        assert_eq!(int_as_string.actual(), PropertyValueKind::Int);

        let bool_as_double = PropertyValue::Bool(true).as_double().unwrap_err();
        assert_eq!(bool_as_double.expected(), PropertyValueKind::Double);
        assert_eq!(bool_as_double.actual(), PropertyValueKind::Bool);
    }
}
