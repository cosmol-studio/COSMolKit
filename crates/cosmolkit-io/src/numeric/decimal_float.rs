//! musl `decfloat` (1.2.5) for the selected binary64 configuration.
//!
//! Ported from `third_party/musl/src/internal/floatscan.c:63-312`. The base-1e9
//! ring buffer, alignment, up/down scaling, bias rounding and scaling are
//! transliterated in source order. `KMAX=128`, `LD_B1B_DIG=2` and
//! `LD_B1B_MAX={9007199,254740991}` are the binary64 constants. The source uses
//! plain `int` indices with `& MASK`; the port keeps signed indices so the
//! two's-complement wrap behavior is identical.

use super::cursor::ScanCursor;
use super::scalbn::{fmod, scalbn};

const KMAX: i32 = 128;
const MASK: i32 = KMAX - 1;
const LD_B1B_DIG: i32 = 2;
const TH: [u32; 2] = [9007199, 254740991];
const P10S: [i32; 8] = [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000];
const LDBL_MANT_DIG: i32 = 53;

#[inline]
fn idx(value: i32) -> usize {
    (value & MASK) as usize
}

/// musl `floatscan.c:63-312::decfloat` (1.2.5), binary64.
pub(super) fn decfloat(
    cursor: &mut ScanCursor,
    mut c: i32,
    bits_in: i32,
    emin: i32,
    sign: i32,
    pok: bool,
) -> f64 {
    // musl✔️✔️: static long double decfloat(FILE *f, int c, int bits, int emin, int sign, int pok)
    // musl✔️✔️: {
    // musl✔️✔️: 	uint32_t x[KMAX];
    // musl✔️✔️: 	static const uint32_t th[] = { LD_B1B_MAX };
    // musl✔️✔️: 	int i, j, k, a, z;
    // musl✔️✔️: 	long long lrp=0, dc=0;
    // musl✔️✔️: 	long long e10=0;
    // musl✔️✔️: 	int lnz = 0;
    // musl✔️✔️: 	int gotdig = 0, gotrad = 0;
    // musl✔️✔️: 	int rp;
    // musl✔️✔️: 	int e2;
    // musl✔️✔️: 	int emax = -emin-bits+3;
    // musl✔️✔️: 	int denormal = 0;
    // musl✔️✔️: 	long double y;
    // musl✔️✔️: 	long double frac=0;
    // musl✔️✔️: 	long double bias=0;
    // musl✔️✔️: 	static const int p10s[] = { 10, 100, 1000, 10000,
    // musl✔️✔️: 		100000, 1000000, 10000000, 100000000 };
    // musl✔️✔️:
    // musl✔️✔️: 	j=0;
    // musl✔️✔️: 	k=0;
    // musl✔️✔️:
    // musl✔️✔️: 	/* Don't let leading zeros consume buffer space */
    // musl✔️✔️: 	for (; c=='0'; c = shgetc(f)) gotdig=1;
    // musl✔️✔️: 	if (c=='.') {
    // musl✔️✔️: 		gotrad = 1;
    // musl✔️✔️: 		for (c = shgetc(f); c=='0'; c = shgetc(f)) gotdig=1, lrp--;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	x[0] = 0;
    // musl✔️✔️: 	for (; c-'0'<10U || c=='.'; c = shgetc(f)) {
    // musl✔️✔️: 		if (c == '.') {
    // musl✔️✔️: 			if (gotrad) break;
    // musl✔️✔️: 			gotrad = 1;
    // musl✔️✔️: 			lrp = dc;
    // musl✔️✔️: 		} else if (k < KMAX-3) {
    // musl✔️✔️: 			dc++;
    // musl✔️✔️: 			if (c!='0') lnz = dc;
    // musl✔️✔️: 			if (j) x[k] = x[k]*10 + c-'0';
    // musl✔️✔️: 			else x[k] = c-'0';
    // musl✔️✔️: 			if (++j==9) {
    // musl✔️✔️: 				k++;
    // musl✔️✔️: 				j=0;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			gotdig=1;
    // musl✔️✔️: 		} else {
    // musl✔️✔️: 			dc++;
    // musl✔️✔️: 			if (c!='0') {
    // musl✔️✔️: 				lnz = (KMAX-4)*9;
    // musl✔️✔️: 				x[KMAX-4] |= 1;
    // musl✔️✔️: 			}
    // musl✔️✔️: 		}
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (!gotrad) lrp=dc;
    // musl✔️✔️:
    // musl✔️✔️: 	if (gotdig && (c|32)=='e') {
    // musl✔️✔️: 		e10 = scanexp(f, pok);
    // musl✔️✔️: 		if (e10 == LLONG_MIN) {
    // musl✔️✔️: 			if (pok) {
    // musl✔️✔️: 				shunget(f);
    // musl✔️✔️: 			} else {
    // musl✔️✔️: 				shlim(f, 0);
    // musl✔️✔️: 				return 0;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			e10 = 0;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		lrp += e10;
    // musl✔️✔️: 	} else if (c>=0) {
    // musl✔️✔️: 		shunget(f);
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (!gotdig) {
    // musl✔️✔️: 		errno = EINVAL;
    // musl✔️✔️: 		shlim(f, 0);
    // musl✔️✔️: 		return 0;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Handle zero specially to avoid nasty special cases later */
    // musl✔️✔️: 	if (!x[0]) return sign * 0.0;
    // musl✔️✔️:
    // musl✔️✔️: 	/* Optimize small integers (w/no exponent) and over/under-flow */
    // musl✔️✔️: 	if (lrp==dc && dc<10 && (bits>30 || x[0]>>bits==0))
    // musl✔️✔️: 		return sign * (long double)x[0];
    // musl✔️✔️: 	if (lrp > -emin/2) {
    // musl✔️✔️: 		errno = ERANGE;
    // musl✔️✔️: 		return sign * LDBL_MAX * LDBL_MAX;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (lrp < emin-2*LDBL_MANT_DIG) {
    // musl✔️✔️: 		errno = ERANGE;
    // musl✔️✔️: 		return sign * LDBL_MIN * LDBL_MIN;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Align incomplete final B1B digit */
    // musl✔️✔️: 	if (j) {
    // musl✔️✔️: 		for (; j<9; j++) x[k]*=10;
    // musl✔️✔️: 		k++;
    // musl✔️✔️: 		j=0;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	a = 0;
    // musl✔️✔️: 	z = k;
    // musl✔️✔️: 	e2 = 0;
    // musl✔️✔️: 	rp = lrp;
    // musl✔️✔️:
    // musl✔️✔️: 	/* Optimize small to mid-size integers (even in exp. notation) */
    // musl✔️✔️: 	if (lnz<9 && lnz<=rp && rp < 18) {
    // musl✔️✔️: 		if (rp == 9) return sign * (long double)x[0];
    // musl✔️✔️: 		if (rp < 9) return sign * (long double)x[0] / p10s[8-rp];
    // musl✔️✔️: 		int bitlim = bits-3*(int)(rp-9);
    // musl✔️✔️: 		if (bitlim>30 || x[0]>>bitlim==0)
    // musl✔️✔️: 			return sign * (long double)x[0] * p10s[rp-10];
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Drop trailing zeros */
    // musl✔️✔️: 	for (; !x[z-1]; z--);
    // musl✔️✔️:
    // musl✔️✔️: 	/* Align radix point to B1B digit boundary */
    // musl✔️✔️: 	if (rp % 9) {
    // musl✔️✔️: 		int rpm9 = rp>=0 ? rp%9 : rp%9+9;
    // musl✔️✔️: 		int p10 = p10s[8-rpm9];
    // musl✔️✔️: 		uint32_t carry = 0;
    // musl✔️✔️: 		for (k=a; k!=z; k++) {
    // musl✔️✔️: 			uint32_t tmp = x[k] % p10;
    // musl✔️✔️: 			x[k] = x[k]/p10 + carry;
    // musl✔️✔️: 			carry = 1000000000/p10 * tmp;
    // musl✔️✔️: 			if (k==a && !x[k]) {
    // musl✔️✔️: 				a = (a+1 & MASK);
    // musl✔️✔️: 				rp -= 9;
    // musl✔️✔️: 			}
    // musl✔️✔️: 		}
    // musl✔️✔️: 		if (carry) x[z++] = carry;
    // musl✔️✔️: 		rp += 9-rpm9;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Upscale until desired number of bits are left of radix point */
    // musl✔️✔️: 	while (rp < 9*LD_B1B_DIG || (rp == 9*LD_B1B_DIG && x[a]<th[0])) {
    // musl✔️✔️: 		uint32_t carry = 0;
    // musl✔️✔️: 		e2 -= 29;
    // musl✔️✔️: 		for (k=(z-1 & MASK); ; k=(k-1 & MASK)) {
    // musl✔️✔️: 			uint64_t tmp = ((uint64_t)x[k] << 29) + carry;
    // musl✔️✔️: 			if (tmp > 1000000000) {
    // musl✔️✔️: 				carry = tmp / 1000000000;
    // musl✔️✔️: 				x[k] = tmp % 1000000000;
    // musl✔️✔️: 			} else {
    // musl✔️✔️: 				carry = 0;
    // musl✔️✔️: 				x[k] = tmp;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			if (k==(z-1 & MASK) && k!=a && !x[k]) z = k;
    // musl✔️✔️: 			if (k==a) break;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		if (carry) {
    // musl✔️✔️: 			rp += 9;
    // musl✔️✔️: 			a = (a-1 & MASK);
    // musl✔️✔️: 			if (a == z) {
    // musl✔️✔️: 				z = (z-1 & MASK);
    // musl✔️✔️: 				x[z-1 & MASK] |= x[z];
    // musl✔️✔️: 			}
    // musl✔️✔️: 			x[a] = carry;
    // musl✔️✔️: 		}
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Downscale until exactly number of bits are left of radix point */
    // musl✔️✔️: 	for (;;) {
    // musl✔️✔️: 		uint32_t carry = 0;
    // musl✔️✔️: 		int sh = 1;
    // musl✔️✔️: 		for (i=0; i<LD_B1B_DIG; i++) {
    // musl✔️✔️: 			k = (a+i & MASK);
    // musl✔️✔️: 			if (k == z || x[k] < th[i]) {
    // musl✔️✔️: 				i=LD_B1B_DIG;
    // musl✔️✔️: 				break;
    // musl✔️✔️: 			}
    // musl✔️✔️: 			if (x[a+i & MASK] > th[i]) break;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		if (i==LD_B1B_DIG && rp==9*LD_B1B_DIG) break;
    // musl✔️✔️: 		/* FIXME: find a way to compute optimal sh */
    // musl✔️✔️: 		if (rp > 9+9*LD_B1B_DIG) sh = 9;
    // musl✔️✔️: 		e2 += sh;
    // musl✔️✔️: 		for (k=a; k!=z; k=(k+1 & MASK)) {
    // musl✔️✔️: 			uint32_t tmp = x[k] & (1<<sh)-1;
    // musl✔️✔️: 			x[k] = (x[k]>>sh) + carry;
    // musl✔️✔️: 			carry = (1000000000>>sh) * tmp;
    // musl✔️✔️: 			if (k==a && !x[k]) {
    // musl✔️✔️: 				a = (a+1 & MASK);
    // musl✔️✔️: 				i--;
    // musl✔️✔️: 				rp -= 9;
    // musl✔️✔️: 			}
    // musl✔️✔️: 		}
    // musl✔️✔️: 		if (carry) {
    // musl✔️✔️: 			if ((z+1 & MASK) != a) {
    // musl✔️✔️: 				x[z] = carry;
    // musl✔️✔️: 				z = (z+1 & MASK);
    // musl✔️✔️: 			} else x[z-1 & MASK] |= 1;
    // musl✔️✔️: 		}
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Assemble desired bits into floating point variable */
    // musl✔️✔️: 	for (y=i=0; i<LD_B1B_DIG; i++) {
    // musl✔️✔️: 		if ((a+i & MASK)==z) x[(z=(z+1 & MASK))-1] = 0;
    // musl✔️✔️: 		y = 1000000000.0L * y + x[a+i & MASK];
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	y *= sign;
    // musl✔️✔️:
    // musl✔️✔️: 	/* Limit precision for denormal results */
    // musl✔️✔️: 	if (bits > LDBL_MANT_DIG+e2-emin) {
    // musl✔️✔️: 		bits = LDBL_MANT_DIG+e2-emin;
    // musl✔️✔️: 		if (bits<0) bits=0;
    // musl✔️✔️: 		denormal = 1;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Calculate bias term to force rounding, move out lower bits */
    // musl✔️✔️: 	if (bits < LDBL_MANT_DIG) {
    // musl✔️✔️: 		bias = copysignl(scalbn(1, 2*LDBL_MANT_DIG-bits-1), y);
    // musl✔️✔️: 		frac = fmodl(y, scalbn(1, LDBL_MANT_DIG-bits));
    // musl✔️✔️: 		y -= frac;
    // musl✔️✔️: 		y += bias;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* Process tail of decimal input so it can affect rounding */
    // musl✔️✔️: 	if ((a+i & MASK) != z) {
    // musl✔️✔️: 		uint32_t t = x[a+i & MASK];
    // musl✔️✔️: 		if (t < 500000000 && (t || (a+i+1 & MASK) != z))
    // musl✔️✔️: 			frac += 0.25*sign;
    // musl✔️✔️: 		else if (t > 500000000)
    // musl✔️✔️: 			frac += 0.75*sign;
    // musl✔️✔️: 		else if (t == 500000000) {
    // musl✔️✔️: 			if ((a+i+1 & MASK) == z)
    // musl✔️✔️: 				frac += 0.5*sign;
    // musl✔️✔️: 			else
    // musl✔️✔️: 				frac += 0.75*sign;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		if (LDBL_MANT_DIG-bits >= 2 && !fmodl(frac, 1))
    // musl✔️✔️: 			frac++;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	y += frac;
    // musl✔️✔️: 	y -= bias;
    // musl✔️✔️:
    // musl✔️✔️: 	if ((e2+LDBL_MANT_DIG & INT_MAX) > emax-5) {
    // musl✔️✔️: 		if (fabsl(y) >= 2/LDBL_EPSILON) {
    // musl✔️✔️: 			if (denormal && bits==LDBL_MANT_DIG+e2-emin)
    // musl✔️✔️: 				denormal = 0;
    // musl✔️✔️: 			y *= 0.5;
    // musl✔️✔️: 			e2++;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		if (e2+LDBL_MANT_DIG>emax || (denormal && frac))
    // musl✔️✔️: 			errno = ERANGE;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	return scalbnl(y, e2);
    // musl✔️✔️: }
    let mut x = [0u32; KMAX as usize];
    let mut i: i32;
    let mut j: i32 = 0;
    let mut k: i32 = 0;
    let mut a: i32;
    let mut z: i32;
    let mut lrp: i64 = 0;
    let mut dc: i64 = 0;
    let mut e10: i64 = 0;
    let mut lnz: i32 = 0;
    let mut gotdig = false;
    let mut gotrad = false;
    let mut rp: i32;
    let mut e2: i32;
    let emax: i32 = -emin - bits_in + 3;
    let mut denormal = false;
    let mut y: f64;
    let mut frac: f64 = 0.0;
    let mut bias: f64 = 0.0;
    let sign_f = sign as f64;

    while c == i32::from(b'0') {
        gotdig = true;
        c = cursor.getc();
    }
    if c == i32::from(b'.') {
        gotrad = true;
        c = cursor.getc();
        while c == i32::from(b'0') {
            gotdig = true;
            lrp -= 1;
            c = cursor.getc();
        }
    }

    x[0] = 0;
    loop {
        let is_digit = (c.wrapping_sub(i32::from(b'0')) as u32) < 10;
        if !(is_digit || c == i32::from(b'.')) {
            break;
        }
        if c == i32::from(b'.') {
            if gotrad {
                break;
            }
            gotrad = true;
            lrp = dc;
        } else if k < KMAX - 3 {
            dc += 1;
            if c != i32::from(b'0') {
                lnz = dc as i32;
            }
            if j != 0 {
                x[idx(k)] = x[idx(k)] * 10 + (c - i32::from(b'0')) as u32;
            } else {
                x[idx(k)] = (c - i32::from(b'0')) as u32;
            }
            j += 1;
            if j == 9 {
                k += 1;
                j = 0;
            }
            gotdig = true;
        } else {
            dc += 1;
            if c != i32::from(b'0') {
                lnz = (KMAX - 4) * 9;
                x[idx(KMAX - 4)] |= 1;
            }
        }
        c = cursor.getc();
    }
    if !gotrad {
        lrp = dc;
    }

    if gotdig && (c | 32) == i32::from(b'e') {
        e10 = super::float_scan::scanexp(cursor, pok);
        if e10 == i64::MIN {
            if pok {
                cursor.ungetc();
            } else {
                cursor.reset_count();
                return 0.0;
            }
            e10 = 0;
        }
        lrp += e10;
    } else if c >= 0 {
        cursor.ungetc();
    }
    if !gotdig {
        cursor.reset_count();
        return 0.0;
    }

    if x[0] == 0 {
        return sign_f * 0.0;
    }

    if lrp == dc && dc < 10 && (bits_in > 30 || (x[0] >> bits_in) == 0) {
        return sign_f * f64::from(x[0]);
    }
    if lrp > i64::from(-emin / 2) {
        return sign_f * f64::MAX * f64::MAX;
    }
    if lrp < i64::from(emin - 2 * LDBL_MANT_DIG) {
        return sign_f * f64::MIN_POSITIVE * f64::MIN_POSITIVE;
    }

    if j != 0 {
        while j < 9 {
            x[idx(k)] = x[idx(k)].wrapping_mul(10);
            j += 1;
        }
        k += 1;
        j = 0;
    }

    a = 0;
    z = k;
    e2 = 0;
    rp = lrp as i32;

    if lnz < 9 && lnz <= rp && rp < 18 {
        if rp == 9 {
            return sign_f * f64::from(x[0]);
        }
        if rp < 9 {
            return sign_f * f64::from(x[0]) / f64::from(P10S[(8 - rp) as usize]);
        }
        let bitlim = bits_in - 3 * (rp - 9);
        if bitlim > 30 || (x[0] >> bitlim) == 0 {
            return sign_f * f64::from(x[0]) * f64::from(P10S[(rp - 10) as usize]);
        }
    }

    while x[idx(z - 1)] == 0 {
        z -= 1;
    }

    if rp % 9 != 0 {
        let rpm9 = if rp >= 0 { rp % 9 } else { rp % 9 + 9 };
        let p10 = P10S[(8 - rpm9) as usize];
        let mut carry: u32 = 0;
        k = a;
        while k != z {
            let tmp = x[idx(k)] % p10 as u32;
            x[idx(k)] = x[idx(k)] / p10 as u32 + carry;
            carry = 1000000000u32 / p10 as u32 * tmp;
            if k == a && x[idx(k)] == 0 {
                a = (a + 1) & MASK;
                rp -= 9;
            }
            k += 1;
        }
        if carry != 0 {
            x[idx(z)] = carry;
            z += 1;
        }
        rp += 9 - rpm9;
    }

    while rp < 9 * LD_B1B_DIG || (rp == 9 * LD_B1B_DIG && x[idx(a)] < TH[0]) {
        let mut carry: u32 = 0;
        e2 -= 29;
        k = (z - 1) & MASK;
        loop {
            let tmp = (u64::from(x[idx(k)]) << 29) + u64::from(carry);
            if tmp > 1000000000 {
                carry = (tmp / 1000000000) as u32;
                x[idx(k)] = (tmp % 1000000000) as u32;
            } else {
                carry = 0;
                x[idx(k)] = tmp as u32;
            }
            if k == ((z - 1) & MASK) && k != a && x[idx(k)] == 0 {
                z = k;
            }
            if k == a {
                break;
            }
            k = (k - 1) & MASK;
        }
        if carry != 0 {
            rp += 9;
            a = (a - 1) & MASK;
            if a == z {
                z = (z - 1) & MASK;
                x[idx(z - 1)] |= x[idx(z)];
            }
            x[idx(a)] = carry;
        }
    }

    loop {
        let mut carry: u32 = 0;
        let mut sh: u32 = 1;
        i = 0;
        while i < LD_B1B_DIG {
            k = (a + i) & MASK;
            if k == z || x[idx(k)] < TH[i as usize] {
                i = LD_B1B_DIG;
                break;
            }
            if x[idx(a + i)] > TH[i as usize] {
                break;
            }
            i += 1;
        }
        if i == LD_B1B_DIG && rp == 9 * LD_B1B_DIG {
            break;
        }
        if rp > 9 + 9 * LD_B1B_DIG {
            sh = 9;
        }
        e2 += sh as i32;
        k = a;
        while k != z {
            let tmp = x[idx(k)] & ((1u32 << sh) - 1);
            x[idx(k)] = (x[idx(k)] >> sh) + carry;
            carry = (1000000000u32 >> sh) * tmp;
            if k == a && x[idx(k)] == 0 {
                a = (a + 1) & MASK;
                i -= 1;
                rp -= 9;
            }
            k = (k + 1) & MASK;
        }
        if carry != 0 {
            if ((z + 1) & MASK) != a {
                x[idx(z)] = carry;
                z = (z + 1) & MASK;
            } else {
                x[idx(z - 1)] |= 1;
            }
        }
    }

    y = 0.0;
    i = 0;
    while i < LD_B1B_DIG {
        if (a + i) & MASK == z {
            z = (z + 1) & MASK;
            x[idx(z - 1)] = 0;
        }
        y = 1000000000.0 * y + f64::from(x[idx(a + i)]);
        i += 1;
    }

    y *= sign_f;

    let mut bits = bits_in;
    if bits > LDBL_MANT_DIG + e2 - emin {
        bits = LDBL_MANT_DIG + e2 - emin;
        if bits < 0 {
            bits = 0;
        }
        denormal = true;
    }

    if bits < LDBL_MANT_DIG {
        bias = scalbn(1.0, 2 * LDBL_MANT_DIG - bits - 1).copysign(y);
        frac = fmod(y, scalbn(1.0, LDBL_MANT_DIG - bits));
        y -= frac;
        y += bias;
    }

    if (a + i) & MASK != z {
        let t = x[idx(a + i)];
        if t < 500000000 && (t != 0 || ((a + i + 1) & MASK) != z) {
            frac += 0.25 * sign_f;
        } else if t > 500000000 {
            frac += 0.75 * sign_f;
        } else if t == 500000000 {
            if ((a + i + 1) & MASK) == z {
                frac += 0.5 * sign_f;
            } else {
                frac += 0.75 * sign_f;
            }
        }
        if LDBL_MANT_DIG - bits >= 2 && fmod(frac, 1.0) == 0.0 {
            frac += 1.0;
        }
    }

    y += frac;
    y -= bias;

    if ((e2 + LDBL_MANT_DIG) & i32::MAX) > emax - 5 {
        if y.abs() >= 2.0 / f64::EPSILON {
            if denormal && bits == LDBL_MANT_DIG + e2 - emin {
                denormal = false;
            }
            y *= 0.5;
            e2 += 1;
        }
    }

    scalbn(y, e2)
}
