//! Private C-locale byte-prefix float conversion with pinned source branches.
//!
//! Provenance: musl decimal conversion and pinned runtime glibc hex/NaN branches.
//! Decimal conversion reuses the source-backed musl port. Hexadecimal
//! rounding and signed NaN payloads use the pinned glibc 2.43 source branches.
//! Acceptance remains independent; no numerical input exemptions are applied.
//! This module is private to `cosmolkit-io`; no public API is added.
//!
//! Source basis: musl 1.2.5 (`third_party/musl/`), selected binary64
//! configuration (`LDBL_MANT_DIG == 53 && LDBL_MAX_EXP == 1024`), prec=1,
//! pok=1. Verbatim source anchors live inside the implementing functions.

mod cursor;
mod decimal_float;
mod float_scan;
mod hex_float;
mod scalbn;

pub(crate) use float_scan::float_scan;

/// Parse the longest numeric prefix of `bytes` as a binary64 value and return
/// it with the number of consumed bytes.
pub(crate) fn parse_rdkit_atof_prefix(bytes: &[u8]) -> (f64, usize) {
    // Pinned RDKit source calls atof. The original source-backed decimal
    // converter is reused; hexadecimal conversion and signed NaN payloads
    // use complete glibc 2.43 branches in their implementing functions.
    // C locale, FE_TONEAREST, prefix/no-conversion semantics are required.
    // No input class or known counterexample is waived from source parity.
    // Missing fields and per-format screening remain reader responsibilities.
    // This helper is not a complete general locale/fenv glibc implementation.
    // musl❗✔️: static long double strtox(const char *s, char **p, int prec)
    // musl❗✔️: {
    // musl❗✔️: 	FILE f;
    // musl❗✔️: 	sh_fromstring(&f, s);
    // musl❗✔️: 	shlim(&f, 0);
    // musl❗✔️: 	long double y = __floatscan(&f, prec, 1);
    // musl❗✔️: 	off_t cnt = shcnt(&f);
    // musl❗✔️: 	if (p) *p = cnt ? (char *)s + cnt : (char *)s;
    // musl❗✔️: 	return y;
    // musl❗✔️: }
    // Specialization: prec=1, binary64, return the offset rather than a pointer.
    let mut cursor = cursor::ScanCursor::new(bytes);
    let value = float_scan(&mut cursor, 1, true);
    (value, cursor.consumed())
}

#[cfg(test)]
mod coordinate_contract_tests;

#[cfg(test)]
mod tests {
    use super::parse_rdkit_atof_prefix;

    fn bits(value: f64) -> u64 {
        value.to_bits()
    }

    fn run(token: &[u8]) -> (u64, usize) {
        let (value, offset) = parse_rdkit_atof_prefix(token);
        (bits(value), offset)
    }

    #[test]
    fn atof_source_exact_glibc_hex_and_nan_branches() {
        for (token, expected, consumed) in [
            ("0x1.00000003p-1044", 0x0000_0000_4000_0001, 18),
            ("0x1.00000005p-1044", 0x0000_0000_4000_0001, 18),
            ("-nan", 0xfff8_0000_0000_0000, 4),
            ("nan(123)", 0x7ff8_0000_0000_007b, 8),
            ("nan(0173)", 0x7ff8_0000_0000_007b, 9),
            ("-NAN(0x7b)", 0xfff8_0000_0000_007b, 10),
            ("nan(09)", 0x7ff8_0000_0000_0000, 7),
            ("nan(18446744073709551616)", 0x7fff_ffff_ffff_ffff, 25),
            ("nan(0b1)", 0x7ff8_0000_0000_0000, 8),
            ("nan(abc)", 0x7ff8_0000_0000_0000, 8),
            ("nan(123", 0x7ff8_0000_0000_0000, 3),
            ("-nan(12 3)", 0xfff8_0000_0000_0000, 4),
            ("0x1p-1075", 0, 9),
            ("0x1.0000000000001p-1075", 1, 23),
            ("0x1.8p-1074", 2, 11),
            ("0x1.4p-1074", 1, 11),
            ("0x1.fffffffffffffp-1023", 0x0010_0000_0000_0000, 23),
            ("0x1.fffffffffffffp1023", 0x7fef_ffff_ffff_ffff, 22),
            ("0x1.fffffffffffff8p1023", 0x7ff0_0000_0000_0000, 23),
            ("0x.0008p0tail", 0x3f20_0000_0000_0000, 9),
        ] {
            assert_eq!(run(token.as_bytes()), (expected, consumed), "{token}");
        }
    }

    #[test]
    fn atof_port_hex_complete_counterexamples_and_bases() {
        let fifty = b"9".repeat(50);
        let mut overflow = b"0x10p".to_vec();
        overflow.extend_from_slice(&fifty);
        assert_eq!(run(&overflow), (0x7ff0_0000_0000_0000, 55));

        let mut underflow = b"0x0.1p-".to_vec();
        underflow.extend_from_slice(&fifty);
        assert_eq!(run(&underflow), (0x0000_0000_0000_0000, 57));

        assert_eq!(run(b"0x10"), (0x4030_0000_0000_0000, 4));
        assert_eq!(run(b"0x1p+1"), (0x4000_0000_0000_0000, 6));
        assert_eq!(run(b"0x1.8p+1"), (0x4008_0000_0000_0000, 8));
        assert_eq!(run(b"0x.8p0"), (0x3fe0_0000_0000_0000, 6));
        assert_eq!(run(b"0x1p-1"), (0x3fe0_0000_0000_0000, 6));
        assert_eq!(run(b"0x1p-1074"), (0x0000_0000_0000_0001, 9));
        assert_eq!(run(b"0x1.00000000000008p0"), (0x3ff0_0000_0000_0000, 20));
        assert_eq!(run(b"0x1.0000000000001p0"), (0x3ff0_0000_0000_0001, 19));
        assert_eq!(run(b"0x"), (0x0000_0000_0000_0000, 1));
        assert_eq!(run(b"0X1P+1"), (0x4000_0000_0000_0000, 6));
    }

    #[test]
    fn atof_port_decimal_complete_boundaries() {
        assert_eq!(run(b"1e400"), (0x7ff0_0000_0000_0000, 5));
        assert_eq!(run(b"1e-400"), (0x0000_0000_0000_0000, 6));
        assert_eq!(run(b"1e-323"), (0x0000_0000_0000_0002, 6));
        assert_eq!(run(b"1e308"), (0x7fe1_ccf3_85eb_c8a0, 5));
        assert_eq!(run(b"0"), (0x0000_0000_0000_0000, 1));
        assert_eq!(run(b"1"), (0x3ff0_0000_0000_0000, 1));
        assert_eq!(run(b"1.5"), (0x3ff8_0000_0000_0000, 3));
        assert_eq!(run(b".5"), (0x3fe0_0000_0000_0000, 2));
        assert_eq!(run(b"5."), (0x4014_0000_0000_0000, 2));
        assert_eq!(run(b"1e"), (0x3ff0_0000_0000_0000, 1));
        assert_eq!(run(b"1e+"), (0x3ff0_0000_0000_0000, 1));
        assert_eq!(run(b"  +3.25e2xyz"), (0x4074_5000_0000_0000, 9));
    }

    #[test]
    fn atof_port_special_and_signed_zero() {
        assert_eq!(run(b"inf"), (0x7ff0_0000_0000_0000, 3));
        assert_eq!(run(b"-infinity"), (0xfff0_0000_0000_0000, 9));
        assert!(f64::from_bits(run(b"nan").0).is_nan());
        assert_eq!(run(b"nan").1, 3);
        assert!(f64::from_bits(run(b"-nan").0).is_nan());
        assert_eq!(run(b"-0"), (0x8000_0000_0000_0000, 2));
        assert_eq!(run(b"+0"), (0x0000_0000_0000_0000, 2));
        assert_eq!(run(b"-0.0"), (0x8000_0000_0000_0000, 4));
    }

    #[test]
    fn atof_port_dispatch_no_conversion_and_suffix() {
        assert_eq!(run(b""), (0x0000_0000_0000_0000, 0));
        assert_eq!(run(b"notanumber"), (0x0000_0000_0000_0000, 0));
        assert_eq!(run(b"  "), (0x0000_0000_0000_0000, 0));
        assert_eq!(run(b"0x1p+1rest"), (0x4000_0000_0000_0000, 6));
        assert_eq!(run(b"12,5"), (0x4028_0000_0000_0000, 2));
        assert_eq!(run(b"\0"), (0x0000_0000_0000_0000, 0));
    }
}
