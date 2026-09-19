//! musl `hexfloat` (1.2.5) for the selected binary64 configuration.
//!
//! Ported from `third_party/musl/src/internal/floatscan.c:314-424`. The
//! arithmetic order is preserved: the first eight hexadecimal digits go into a
//! `u32` `x`, the next six into a binary64 `y` fraction, and the remaining
//! nonzero digits become one sticky half-step, exactly as in the source.

use super::cursor::ScanCursor;
use super::scalbn::scalbn;

/// musl `floatscan.c:314-424::hexfloat` (1.2.5), binary64.
pub(super) fn hexfloat(
    cursor: &mut ScanCursor,
    bits_in: i32,
    emin: i32,
    sign: i32,
    pok: bool,
) -> f64 {
    // musl✔️✔️: static long double hexfloat(FILE *f, int bits, int emin, int sign, int pok)
    // musl✔️✔️: {
    // musl✔️✔️: 	uint32_t x = 0;
    // musl✔️✔️: 	long double y = 0;
    // musl✔️✔️: 	long double scale = 1;
    // musl✔️✔️: 	long double bias = 0;
    // musl✔️✔️: 	int gottail = 0, gotrad = 0, gotdig = 0;
    // musl✔️✔️: 	long long rp = 0;
    // musl✔️✔️: 	long long dc = 0;
    // musl✔️✔️: 	long long e2 = 0;
    // musl✔️✔️: 	int d;
    // musl✔️✔️: 	int c;
    // musl✔️✔️:
    // musl✔️✔️: 	c = shgetc(f);
    // musl✔️✔️:
    // musl✔️✔️: 	/* Skip leading zeros */
    // musl✔️✔️: 	for (; c=='0'; c = shgetc(f)) gotdig = 1;
    // musl✔️✔️:
    // musl✔️✔️: 	if (c=='.') {
    // musl✔️✔️: 		gotrad = 1;
    // musl✔️✔️: 		c = shgetc(f);
    // musl✔️✔️: 		/* Count zeros after the radix point before significand */
    // musl✔️✔️: 		for (rp=0; c=='0'; c = shgetc(f), rp--) gotdig = 1;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	for (; c-'0'<10U || (c|32)-'a'<6U || c=='.'; c = shgetc(f)) {
    // musl✔️✔️: 		if (c=='.') {
    // musl✔️✔️: 			if (gotrad) break;
    // musl✔️✔️: 			rp = dc;
    // musl✔️✔️: 			gotrad = 1;
    // musl✔️✔️: 		} else {
    // musl✔️✔️: 			gotdig = 1;
    // musl✔️✔️: 			if (c > '9') d = (c|32)+10-'a';
    // musl✔️✔️: 			else d = c-'0';
    // musl✔️✔️: 			if (dc<8) {
    // musl✔️✔️: 				x = x*16 + d;
    // musl✔️✔️: 			} else if (dc < LDBL_MANT_DIG/4+1) {
    // musl✔️✔️: 				y += d*(scale/=16);
    // musl✔️✔️: 			} else if (d && !gottail) {
    // musl✔️✔️: 				y += 0.5*scale;
    // musl✔️✔️: 				gottail = 1;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			dc++;
    // musl✔️✔️: 		}
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (!gotdig) {
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 		if (pok) {
    // musl✔️✔️: 			shunget(f);
    // musl✔️✔️: 			if (gotrad) shunget(f);
    // musl✔️✔️: 		} else {
    // musl✔️✔️: 			shlim(f, 0);
    // musl✔️✔️: 		}
    // musl✔️✔️: 		return sign * 0.0;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (!gotrad) rp = dc;
    // musl✔️✔️: 	while (dc<8) x *= 16, dc++;
    // musl✔️✔️: 	if ((c|32)=='p') {
    // musl✔️✔️: 		e2 = scanexp(f, pok);
    // musl✔️✔️: 		if (e2 == LLONG_MIN) {
    // musl✔️✔️: 			if (pok) {
    // musl✔️✔️: 				shunget(f);
    // musl✔️✔️: 			} else {
    // musl✔️✔️: 				shlim(f, 0);
    // musl✔️✔️: 				return 0;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			e2 = 0;
    // musl✔️✔️: 		}
    // musl✔️✔️: 	} else {
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 	}
    // musl✔️✔️: 	e2 += 4*rp - 32;
    // musl✔️✔️:
    // musl✔️✔️: 	if (!x) return sign * 0.0;
    // musl✔️✔️: 	if (e2 > -emin) {
    // musl✔️✔️: 		errno = ERANGE;
    // musl✔️✔️: 		return sign * LDBL_MAX * LDBL_MAX;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (e2 < emin-2*LDBL_MANT_DIG) {
    // musl✔️✔️: 		errno = ERANGE;
    // musl✔️✔️: 		return sign * LDBL_MIN * LDBL_MIN;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	while (x < 0x80000000) {
    // musl✔️✔️: 		if (y>=0.5) {
    // musl✔️✔️: 			x += x + 1;
    // musl✔️✔️: 			y += y - 1;
    // musl✔️✔️: 		} else {
    // musl✔️✔️: 			x += x;
    // musl✔️✔️: 			y += y;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		e2--;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	if (bits > 32+e2-emin) {
    // musl✔️✔️: 		bits = 32+e2-emin;
    // musl✔️✔️: 		if (bits<0) bits=0;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	if (bits < LDBL_MANT_DIG)
    // musl✔️✔️: 		bias = copysignl(scalbn(1, 32+LDBL_MANT_DIG-bits-1), sign);
    // musl✔️✔️:
    // musl✔️✔️: 	if (bits<32 && y && !(x&1)) x++, y=0;
    // musl✔️✔️:
    // musl✔️✔️: 	y = bias + sign*(long double)x + sign*y;
    // musl✔️✔️: 	y -= bias;
    // musl✔️✔️:
    // musl✔️✔️: 	if (!y) errno = ERANGE;
    // musl✔️✔️:
    // musl✔️✔️: 	return scalbnl(y, e2);
    // musl✔️✔️: }
    const LDBL_MANT_DIG: i32 = 53;

    let mut x: u32 = 0;
    let mut y: f64 = 0.0;
    let mut scale: f64 = 1.0;
    let mut bias: f64 = 0.0;
    let mut gottail = false;
    let mut gotrad = false;
    let mut gotdig = false;
    let mut rp: i64 = 0;
    let mut dc: i64 = 0;
    let mut e2: i64 = 0;

    let mut c = cursor.getc();

    while c == i32::from(b'0') {
        gotdig = true;
        c = cursor.getc();
    }

    if c == i32::from(b'.') {
        gotrad = true;
        c = cursor.getc();
        rp = 0;
        while c == i32::from(b'0') {
            gotdig = true;
            rp -= 1;
            c = cursor.getc();
        }
    }

    loop {
        let is_hex = ((c.wrapping_sub(i32::from(b'0')) as u32) < 10)
            || (((c | 32).wrapping_sub(i32::from(b'a')) as u32) < 6);
        if !(is_hex || c == i32::from(b'.')) {
            break;
        }
        if c == i32::from(b'.') {
            if gotrad {
                break;
            }
            rp = dc;
            gotrad = true;
        } else {
            gotdig = true;
            let d: i32 = if c > i32::from(b'9') {
                (c | 32) + 10 - i32::from(b'a')
            } else {
                c - i32::from(b'0')
            };
            if dc < 8 {
                x = x.wrapping_mul(16).wrapping_add(d as u32);
            } else if dc < i64::from(LDBL_MANT_DIG / 4 + 1) {
                scale /= 16.0;
                y += f64::from(d) * scale;
            } else if d != 0 && !gottail {
                y += 0.5 * scale;
                gottail = true;
            }
            dc += 1;
        }
        c = cursor.getc();
    }

    if !gotdig {
        cursor.ungetc();
        if pok {
            cursor.ungetc();
            if gotrad {
                cursor.ungetc();
            }
        } else {
            cursor.reset_count();
        }
        return (sign as f64) * 0.0;
    }
    if !gotrad {
        rp = dc;
    }
    while dc < 8 {
        x = x.wrapping_mul(16);
        dc += 1;
    }
    if (c | 32) == i32::from(b'p') {
        e2 = super::float_scan::scanexp(cursor, pok);
        if e2 == i64::MIN {
            if pok {
                cursor.ungetc();
            } else {
                cursor.reset_count();
                return 0.0;
            }
            e2 = 0;
        }
    } else {
        cursor.ungetc();
    }
    e2 += 4 * rp - 32;

    if x == 0 {
        return (sign as f64) * 0.0;
    }
    if e2 > i64::from(-emin) {
        return (sign as f64) * f64::MAX * f64::MAX;
    }
    if e2 < i64::from(emin - 2 * LDBL_MANT_DIG) {
        return (sign as f64) * f64::MIN_POSITIVE * f64::MIN_POSITIVE;
    }

    while x < 0x8000_0000 {
        if y >= 0.5 {
            x = x.wrapping_add(x).wrapping_add(1);
            y += y - 1.0;
        } else {
            x = x.wrapping_add(x);
            y += y;
        }
        e2 -= 1;
    }

    let mut bits = bits_in;
    if i64::from(bits) > 32 + e2 - i64::from(emin) {
        bits = (32 + e2 - i64::from(emin)) as i32;
        if bits < 0 {
            bits = 0;
        }
    }

    if bits < LDBL_MANT_DIG {
        bias = scalbn(1.0, 32 + LDBL_MANT_DIG - bits - 1).copysign(sign as f64);
    }

    if bits < 32 && y != 0.0 && (x & 1) == 0 {
        x += 1;
        y = 0.0;
    }

    y = bias + (sign as f64) * f64::from(x) + (sign as f64) * y;
    y -= bias;

    scalbn(y, e2 as i32)
}
