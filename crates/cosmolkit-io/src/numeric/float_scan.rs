//! musl `__floatscan` dispatcher and `scanexp` exponent scanner (1.2.5).
//!
//! Ported from `third_party/musl/src/internal/floatscan.c:36-60` and
//! `:426-507`, together with the selected binary64 `prec=1`/`pok=1` path.

use super::cursor::ScanCursor;
use super::decimal_float::decfloat;
use super::hex_float::hexfloat;

/// `(unsigned)c - '0' < 10U` for a byte-or-EOF value, matching the source's
/// unsigned digit test.
fn is_digit(c: i32) -> bool {
    (c.wrapping_sub(i32::from(b'0')) as u32) < 10
}

/// `isspace(c)` from musl `include/ctype.h`: `c == ' ' || (unsigned)c-'\t' < 5`.
fn is_space(c: i32) -> bool {
    c == i32::from(b' ') || (c.wrapping_sub(i32::from(b'\t')) as u32) < 5
}

/// musl `floatscan.c:36-60::scanexp` (1.2.5).
///
/// Returns `i64::MIN` for an absent exponent. The `x`/`y` accumulation guards
/// bound the value well inside `i32`/`i64` (see the integration report), and the
/// final digit loop only consumes remaining digits.
pub(super) fn scanexp(cursor: &mut ScanCursor, pok: bool) -> i64 {
    // musl✔️✔️: static long long scanexp(FILE *f, int pok)
    // musl✔️✔️: {
    // musl✔️✔️: 	int c;
    // musl✔️✔️: 	int x;
    // musl✔️✔️: 	long long y;
    // musl✔️✔️: 	int neg = 0;
    // musl✔️✔️:
    // musl✔️✔️: 	c = shgetc(f);
    // musl✔️✔️: 	if (c=='+' || c=='-') {
    // musl✔️✔️: 		neg = (c=='-');
    // musl✔️✔️: 		c = shgetc(f);
    // musl✔️✔️: 		if (c-'0'>=10U && pok) shunget(f);
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (c-'0'>=10U) {
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 		return LLONG_MIN;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	for (x=0; c-'0'<10U && x<INT_MAX/10; c = shgetc(f))
    // musl✔️✔️: 		x = 10*x + c-'0';
    // musl✔️✔️: 	for (y=x; c-'0'<10U && y<LLONG_MAX/100; c = shgetc(f))
    // musl✔️✔️: 		y = 10*y + c-'0';
    // musl✔️✔️: 	for (; c-'0'<10U; c = shgetc(f));
    // musl✔️✔️: 	shunget(f);
    // musl✔️✔️: 	return neg ? -y : y;
    // musl✔️✔️: }
    const INT_MAX: i32 = i32::MAX;
    const LLONG_MAX: i64 = i64::MAX;

    let mut c = cursor.getc();
    let mut neg = false;
    if c == i32::from(b'+') || c == i32::from(b'-') {
        neg = c == i32::from(b'-');
        c = cursor.getc();
        if !is_digit(c) && pok {
            cursor.ungetc();
        }
    }
    if !is_digit(c) {
        cursor.ungetc();
        return i64::MIN;
    }
    let mut x: i32 = 0;
    while is_digit(c) && x < INT_MAX / 10 {
        x = 10 * x + (c - i32::from(b'0'));
        c = cursor.getc();
    }
    let mut y: i64 = i64::from(x);
    while is_digit(c) && y < LLONG_MAX / 100 {
        y = 10 * y + i64::from(c - i32::from(b'0'));
        c = cursor.getc();
    }
    while is_digit(c) {
        c = cursor.getc();
    }
    cursor.ungetc();
    if neg { -y } else { y }
}

/// musl `floatscan.c:426-507::__floatscan` (1.2.5), prec dispatch and
/// sign/infinity/NaN/hexadecimal/zero handling.
///
/// `prec=1` selects binary64 (`bits = DBL_MANT_DIG = 53`,
/// `emin = DBL_MIN_EXP - bits = -1074`). `prec=2` selects the retained binary64
/// long-double configuration, identical to `prec=1`. `prec=0` (f32) is part of
/// the dispatch but outside this port's scope; it is not selected by the V3000
/// call sites.
pub(crate) fn float_scan(cursor: &mut ScanCursor, prec: i32, pok: bool) -> f64 {
    // musl✔️✔️: long double __floatscan(FILE *f, int prec, int pok)
    // musl✔️✔️: {
    // musl✔️✔️: 	int sign = 1;
    // musl✔️✔️: 	size_t i;
    // musl✔️✔️: 	int bits;
    // musl✔️✔️: 	int emin;
    // musl✔️✔️: 	int c;
    // musl✔️✔️:
    // musl✔️✔️: 	switch (prec) {
    // musl✔️✔️: 	case 0:
    // musl✔️✔️: 		bits = FLT_MANT_DIG;
    // musl✔️✔️: 		emin = FLT_MIN_EXP-bits;
    // musl✔️✔️: 		break;
    // musl✔️✔️: 	case 1:
    // musl✔️✔️: 		bits = DBL_MANT_DIG;
    // musl✔️✔️: 		emin = DBL_MIN_EXP-bits;
    // musl✔️✔️: 		break;
    // musl✔️✔️: 	case 2:
    // musl✔️✔️: 		bits = LDBL_MANT_DIG;
    // musl✔️✔️: 		emin = LDBL_MIN_EXP-bits;
    // musl✔️✔️: 		break;
    // musl✔️✔️: 	default:
    // musl✔️✔️: 		return 0;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	while (isspace((c=shgetc(f))));
    // musl✔️✔️:
    // musl✔️✔️: 	if (c=='+' || c=='-') {
    // musl✔️✔️: 		sign -= 2*(c=='-');
    // musl✔️✔️: 		c = shgetc(f);
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	for (i=0; i<8 && (c|32)=="infinity"[i]; i++)
    // musl✔️✔️: 		if (i<7) c = shgetc(f);
    // musl✔️✔️: 	if (i==3 || i==8 || (i>3 && pok)) {
    // musl✔️✔️: 		if (i!=8) {
    // musl✔️✔️: 			shunget(f);
    // musl✔️✔️: 			if (pok) for (; i>3; i--) shunget(f);
    // musl✔️✔️: 		}
    // musl✔️✔️: 		return sign * INFINITY;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (!i) for (i=0; i<3 && (c|32)=="nan"[i]; i++)
    // musl✔️✔️: 		if (i<2) c = shgetc(f);
    // musl✔️✔️: 	if (i==3) {
    // musl✔️✔️: 		if (shgetc(f) != '(') {
    // musl✔️✔️: 			shunget(f);
    // musl✔️✔️: 			return NAN;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		for (i=1; ; i++) {
    // musl✔️✔️: 			c = shgetc(f);
    // musl✔️✔️: 			if (c-'0'<10U || c-'A'<26U || c-'a'<26U || c=='_')
    // musl✔️✔️: 				continue;
    // musl✔️✔️: 			if (c==')') return NAN;
    // musl✔️✔️: 			shunget(f);
    // musl✔️✔️: 			if (!pok) {
    // musl✔️✔️: 				errno = EINVAL;
    // musl✔️✔️: 				shlim(f, 0);
    // musl✔️✔️: 				return 0;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			while (i--) shunget(f);
    // musl✔️✔️: 			return NAN;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		return NAN;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	if (i) {
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 		errno = EINVAL;
    // musl✔️✔️: 		shlim(f, 0);
    // musl✔️✔️: 		return 0;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	if (c=='0') {
    // musl✔️✔️: 		c = shgetc(f);
    // musl✔️✔️: 		if ((c|32) == 'x')
    // musl✔️✔️: 			return hexfloat(f, bits, emin, sign, pok);
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 		c = '0';
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	return decfloat(f, c, bits, emin, sign, pok);
    // musl✔️✔️: }
    const INFINITY: &[u8; 8] = b"infinity";
    const NAN_WORD: &[u8; 3] = b"nan";

    let (bits, emin) = match prec {
        0 => (24i32, -125 - 24),
        1 => (53i32, -1021 - 53),
        2 => (53i32, -1021 - 53),
        _ => return 0.0,
    };

    let mut sign = 1i32;
    let mut c = cursor.getc();
    while is_space(c) {
        c = cursor.getc();
    }

    if c == i32::from(b'+') || c == i32::from(b'-') {
        if c == i32::from(b'-') {
            sign = -1;
        }
        c = cursor.getc();
    }

    let mut i = 0usize;
    while i < 8 && (c | 32) == i32::from(INFINITY[i]) {
        if i < 7 {
            c = cursor.getc();
        }
        i += 1;
    }
    if i == 3 || i == 8 || (i > 3 && pok) {
        if i != 8 {
            cursor.ungetc();
            if pok {
                while i > 3 {
                    cursor.ungetc();
                    i -= 1;
                }
            }
        }
        return (sign as f64) * f64::INFINITY;
    }

    if i == 0 {
        while i < 3 && (c | 32) == i32::from(NAN_WORD[i]) {
            if i < 2 {
                c = cursor.getc();
            }
            i += 1;
        }
    }
    if i == 3 {
        if cursor.getc() != i32::from(b'(') {
            cursor.ungetc();
            return f64::NAN;
        }
        let mut i = 1usize;
        loop {
            c = cursor.getc();
            let alnum = is_digit(c)
                || ((c.wrapping_sub(i32::from(b'A')) as u32) < 26)
                || ((c.wrapping_sub(i32::from(b'a')) as u32) < 26)
                || c == i32::from(b'_');
            if alnum {
                i += 1;
                continue;
            }
            if c == i32::from(b')') {
                return f64::NAN;
            }
            cursor.ungetc();
            if !pok {
                cursor.reset_count();
                return 0.0;
            }
            while i > 0 {
                cursor.ungetc();
                i -= 1;
            }
            return f64::NAN;
        }
    }

    if i != 0 {
        cursor.ungetc();
        cursor.reset_count();
        return 0.0;
    }

    if c == i32::from(b'0') {
        c = cursor.getc();
        if (c | 32) == i32::from(b'x') {
            return hexfloat(cursor, bits, emin, sign, pok);
        }
        cursor.ungetc();
        c = i32::from(b'0');
    }

    decfloat(cursor, c, bits, emin, sign, pok)
}

#[cfg(test)]
mod tests {
    use super::{is_digit, is_space, scanexp};
    use crate::numeric::cursor::ScanCursor;

    #[test]
    fn atof_port_scanexp_digit_classification() {
        assert!(is_digit(i32::from(b'0')));
        assert!(is_digit(i32::from(b'9')));
        assert!(!is_digit(i32::from(b'e')));
        assert!(!is_digit(0));
        assert!(!is_digit(-1));
        assert!(is_space(i32::from(b' ')));
        assert!(is_space(i32::from(b'\t')));
        assert!(is_space(i32::from(b'\n')));
        assert!(!is_space(0));
        assert!(!is_space(i32::from(b'0')));
    }

    fn scan(bytes: &[u8]) -> (i64, usize) {
        let mut cursor = ScanCursor::new(bytes);
        let value = scanexp(&mut cursor, true);
        (value, cursor.consumed())
    }

    #[test]
    fn atof_port_scanexp_signs_and_missing_digits() {
        assert_eq!(scan(b"+12"), (12, 3));
        assert_eq!(scan(b"-12"), (-12, 3));
        assert_eq!(scan(b"0"), (0, 1));
        // Missing digits: the exponent is absent and the cursor rewinds.
        assert_eq!(scan(b"e"), (i64::MIN, 0));
        assert_eq!(scan(b"+"), (i64::MIN, 0));
    }

    #[test]
    fn atof_port_scanexp_long_runs_bound_without_saturation() {
        // Fifty 9s far exceed i64; the source guard stops accumulating and only
        // consumes the remaining digits, so the value stays below LLONG_MAX
        // (about 9.22e17) instead of saturating.
        let token = [b'9'; 50];
        let (value, offset) = scan(&token);
        assert_eq!(offset, 50);
        assert!(value > 0);
        assert!(value < i64::MAX);
        assert_eq!(value, 99999999999999999);
    }

    #[test]
    fn atof_port_scanexp_thousand_digit_exponent_consumes_all() {
        let token = [b'9'; 1000];
        let (value, offset) = scan(&token);
        assert_eq!(offset, 1000);
        assert_eq!(value, 99999999999999999);
    }
}
