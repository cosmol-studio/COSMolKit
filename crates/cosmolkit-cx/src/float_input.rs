//! Original CX numeric boundary regressions, retaining all counted inputs.
//! Numeric behavior has one production owner in foundational CORE.
//! The former private implementation and anchors are frozen in v49 custody;
//! reached source bodies stay inside the sole CORE implementing functions.

#[cfg(test)]
mod tests {
    use cosmolkit_core::source_lexical_double;

    #[test]
    fn byte_float_special_values_preserve_counted_opaque_payload_and_sign() {
        for text in [
            b"nan".as_slice(),
            b"NaN()",
            b"+NAN(x)",
            b"nan(\xff)",
            b"NaN(\0)",
            b"nan(()))",
        ] {
            assert_eq!(
                source_lexical_double(text).unwrap().to_bits(),
                0x7ff8_0000_0000_0000,
                "{text:?}"
            );
        }
        for text in [b"-nan".as_slice(), b"-NAN(\xff\0)"] {
            assert_eq!(
                source_lexical_double(text).unwrap().to_bits(),
                0xfff8_0000_0000_0000,
                "{text:?}"
            );
        }
        for (text, bits) in [
            (b"inf".as_slice(), 0x7ff0_0000_0000_0000),
            (b"+iNf", 0x7ff0_0000_0000_0000),
            (b"-INFINITY", 0xfff0_0000_0000_0000),
        ] {
            assert_eq!(
                source_lexical_double(text).unwrap().to_bits(),
                bits,
                "{text:?}"
            );
        }
    }

    #[test]
    fn byte_float_complete_source_grammar_errors_remain_errors() {
        for text in [
            b"".as_slice(),
            b"+",
            b"-",
            b".",
            b"1e",
            b"1e+",
            b"1e-",
            b" 1",
            b"1 ",
            b"1\t",
            b"1\0",
            b"\xff",
            b"0x1p0",
            b"1_2",
            b"--1",
            b"1..2",
            b"1e2e3",
            b"nan(x)y",
            b"nan(",
            b"infinityx",
        ] {
            assert!(source_lexical_double(text).is_err(), "{text:?}");
        }
    }

    #[test]
    fn byte_float_range_underflow_signed_zero_and_rounding_match_fixed_source_bits() {
        // Fixed boundary results from the actual pinned Boost1.81 reached
        // double-input helper. The ordinary test never invokes an oracle.
        for (text, bits) in [
            (b"-0".as_slice(), 0x8000_0000_0000_0000),
            (b"+0.0", 0),
            (b"-0.000e999999999999999999999", 0x8000_0000_0000_0000),
            (b"1.7976931348623157e308", 0x7fef_ffff_ffff_ffff),
            (b"1.7976931348623158e308", 0x7fef_ffff_ffff_ffff),
            (b"5e-324", 1),
            (b"2.4703282292062327e-324", 0),
            (b"2.4703282292062328e-324", 1),
            (b"1e-999999999999999999999", 0),
            (b"-1e-999999999999999999999", 0x8000_0000_0000_0000),
            (b"2.2250738585072014e-308", 0x0010_0000_0000_0000),
            (b"2.2250738585072011e-308", 0x000f_ffff_ffff_ffff),
            (b"9007199254740993", 0x4340_0000_0000_0000),
            (b"9007199254740995", 0x4340_0000_0000_0002),
            (
                b"1.00000000000000011102230246251565404236316680908203125",
                0x3ff0_0000_0000_0000,
            ),
            (
                b"1.00000000000000033306690738754696212708950042724609375",
                0x3ff0_0000_0000_0002,
            ),
        ] {
            assert_eq!(
                source_lexical_double(text).unwrap().to_bits(),
                bits,
                "{text:?}"
            );
        }
        for text in [
            b"1.7976931348623159e308".as_slice(),
            b"1e309",
            b"-1e309",
            b"1e999999999999999999999",
        ] {
            assert!(source_lexical_double(text).is_err(), "{text:?}");
        }
    }
}
