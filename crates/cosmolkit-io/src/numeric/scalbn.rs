//! Bit-exact power-of-two scaling and remainder helpers from musl 1.2.5.
//!
//! `scalbn` is the binary64 path used by musl `scalbnl` (`scalbnl.c:3-7`) and
//! by the decimal and hexadecimal finalization. `fmod` is the binary64 path
//! used by musl `fmodl` (`fmodl.c:3-7`) for the decimal rounding construction.
//! Both retain the source floating-point operations and integer bit handling.

/// musl `src/math/scalbn.c` (1.2.5), binary64.
///
/// Performance/complexity: identical constant-count bit algorithm; no
/// allocation or extra passes.
pub(crate) fn scalbn(x: f64, mut n: i32) -> f64 {
    // musl✔️✔️: double scalbn(double x, int n)
    // musl✔️✔️: {
    // musl✔️✔️: 	union {double f; uint64_t i;} u;
    // musl✔️✔️: 	double_t y = x;
    // musl✔️✔️:
    // musl✔️✔️: 	if (n > 1023) {
    // musl✔️✔️: 		y *= 0x1p1023;
    // musl✔️✔️: 		n -= 1023;
    // musl✔️✔️: 		if (n > 1023) {
    // musl✔️✔️: 			y *= 0x1p1023;
    // musl✔️✔️: 			n -= 1023;
    // musl✔️✔️: 			if (n > 1023)
    // musl✔️✔️: 				n = 1023;
    // musl✔️✔️: 		}
    // musl✔️✔️: 	} else if (n < -1022) {
    // musl✔️✔️: 		/* make sure final n < -53 to avoid double
    // musl✔️✔️: 		   rounding in the subnormal range */
    // musl✔️✔️: 		y *= 0x1p-1022 * 0x1p53;
    // musl✔️✔️: 		n += 1022 - 53;
    // musl✔️✔️: 		if (n < -1022) {
    // musl✔️✔️: 			y *= 0x1p-1022 * 0x1p53;
    // musl✔️✔️: 			n += 1022 - 53;
    // musl✔️✔️: 			if (n < -1022)
    // musl✔️✔️: 				n = -1022;
    // musl✔️✔️: 		}
    // musl✔️✔️: 	}
    // musl✔️✔️: 	u.i = (uint64_t)(0x3ff+n)<<52;
    // musl✔️✔️: 	x = y * u.f;
    // musl✔️✔️: 	return x;
    // musl✔️✔️: }
    // Exact powers of two as verified literal bits:
    //   0x1p1023  = 0x7FE0000000000000
    //   0x1p-1022 = 0x0010000000000000
    //   0x1p53    = 0x4340000000000000
    const P1023: f64 = f64::from_bits(0x7FE0_0000_0000_0000);
    const P_M1022: f64 = f64::from_bits(0x0010_0000_0000_0000);
    const P53: f64 = f64::from_bits(0x4340_0000_0000_0000);

    let mut y = x;
    if n > 1023 {
        y *= P1023;
        n -= 1023;
        if n > 1023 {
            y *= P1023;
            n -= 1023;
            if n > 1023 {
                n = 1023;
            }
        }
    } else if n < -1022 {
        y *= P_M1022 * P53;
        n += 1022 - 53;
        if n < -1022 {
            y *= P_M1022 * P53;
            n += 1022 - 53;
            if n < -1022 {
                n = -1022;
            }
        }
    }
    let scale = f64::from_bits(((0x3ff + i64::from(n)) as u64) << 52);
    y * scale
}

/// musl `src/math/fmod.c` (1.2.5), binary64.
///
/// Performance/complexity: identical bit-serial algorithm; no allocation.
pub(crate) fn fmod(x: f64, y: f64) -> f64 {
    // musl✔️✔️: double fmod(double x, double y)
    // musl✔️✔️: {
    // musl✔️✔️: 	union {double f; uint64_t i;} ux = {x}, uy = {y};
    // musl✔️✔️: 	int ex = ux.i>>52 & 0x7ff;
    // musl✔️✔️: 	int ey = uy.i>>52 & 0x7ff;
    // musl✔️✔️: 	int sx = ux.i>>63;
    // musl✔️✔️: 	uint64_t i;
    // musl✔️✔️:
    // musl✔️✔️: 	/* in the followings uxi should be ux.i, but then gcc wrongly adds */
    // musl✔️✔️: 	/* float load/store to inner loops ruining performance and code size */
    // musl✔️✔️: 	uint64_t uxi = ux.i;
    // musl✔️✔️:
    // musl✔️✔️: 	if (uy.i<<1 == 0 || isnan(y) || ex == 0x7ff)
    // musl✔️✔️: 		return (x*y)/(x*y);
    // musl✔️✔️: 	if (uxi<<1 <= uy.i<<1) {
    // musl✔️✔️: 		if (uxi<<1 == uy.i<<1)
    // musl✔️✔️: 			return 0*x;
    // musl✔️✔️: 		return x;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* normalize x and y */
    // musl✔️✔️: 	if (!ex) {
    // musl✔️✔️: 		for (i = uxi<<12; i>>63 == 0; ex--, i <<= 1);
    // musl✔️✔️: 		uxi <<= -ex + 1;
    // musl✔️✔️: 	} else {
    // musl✔️✔️: 		uxi &= -1ULL >> 12;
    // musl✔️✔️: 		uxi |= 1ULL << 52;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	if (!ey) {
    // musl✔️✔️: 		for (i = uy.i<<12; i>>63 == 0; ey--, i <<= 1);
    // musl✔️✔️: 		uy.i <<= -ey + 1;
    // musl✔️✔️: 	} else {
    // musl✔️✔️: 		uy.i &= -1ULL >> 12;
    // musl✔️✔️: 		uy.i |= 1ULL << 52;
    // musl✔️✔️: 	}
    // musl✔️✔️:
    // musl✔️✔️: 	/* x mod y */
    // musl✔️✔️: 	for (; ex > ey; ex--) {
    // musl✔️✔️: 		i = uxi - uy.i;
    // musl✔️✔️: 		if (i >> 63 == 0) {
    // musl✔️✔️: 			if (i == 0)
    // musl✔️✔️: 				return 0*x;
    // musl✔️✔️: 			uxi = i;
    // musl✔️✔️: 		}
    // musl✔️✔️: 		uxi <<= 1;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	i = uxi - uy.i;
    // musl✔️✔️: 	if (i >> 63 == 0) {
    // musl✔️✔️: 		if (i == 0)
    // musl✔️✔️: 			return 0*x;
    // musl✔️✔️: 		uxi = i;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	for (; uxi>>52 == 0; uxi <<= 1, ex--);
    // musl✔️✔️:
    // musl✔️✔️: 	/* scale result */
    // musl✔️✔️: 	if (ex > 0) {
    // musl✔️✔️: 		uxi -= 1ULL << 52;
    // musl✔️✔️: 		uxi |= (uint64_t)ex << 52;
    // musl✔️✔️: 	} else {
    // musl✔️✔️: 		uxi >>= -ex + 1;
    // musl✔️✔️: 	}
    // musl✔️✔️: 	uxi |= (uint64_t)sx << 63;
    // musl✔️✔️: 	ux.i = uxi;
    // musl✔️✔️: 	return ux.f;
    // musl✔️✔️: }
    let ux = x.to_bits();
    let mut uy = y.to_bits();
    let mut ex = ((ux >> 52) & 0x7ff) as i32;
    let mut ey = ((uy >> 52) & 0x7ff) as i32;
    let sx = ux >> 63;

    let mut uxi = ux;

    if (uy << 1) == 0 || y.is_nan() || ex == 0x7ff {
        return (x * y) / (x * y);
    }
    if (uxi << 1) <= (uy << 1) {
        if (uxi << 1) == (uy << 1) {
            return 0.0 * x;
        }
        return x;
    }

    if ex == 0 {
        let mut i = uxi << 12;
        while (i >> 63) == 0 {
            ex -= 1;
            i <<= 1;
        }
        uxi <<= -ex + 1;
    } else {
        uxi &= u64::MAX >> 12;
        uxi |= 1u64 << 52;
    }
    if ey == 0 {
        let mut i = uy << 12;
        while (i >> 63) == 0 {
            ey -= 1;
            i <<= 1;
        }
        uy <<= -ey + 1;
    } else {
        uy &= u64::MAX >> 12;
        uy |= 1u64 << 52;
    }

    while ex > ey {
        let i = uxi.wrapping_sub(uy);
        if (i >> 63) == 0 {
            if i == 0 {
                return 0.0 * x;
            }
            uxi = i;
        }
        uxi <<= 1;
        ex -= 1;
    }
    let i = uxi.wrapping_sub(uy);
    if (i >> 63) == 0 {
        if i == 0 {
            return 0.0 * x;
        }
        uxi = i;
    }
    while (uxi >> 52) == 0 {
        uxi <<= 1;
        ex -= 1;
    }

    if ex > 0 {
        uxi -= 1u64 << 52;
        uxi |= (ex as u64) << 52;
    } else {
        uxi >>= -ex + 1;
    }
    uxi |= sx << 63;
    f64::from_bits(uxi)
}

#[cfg(test)]
mod tests {
    use super::{fmod, scalbn};

    fn bits(value: f64) -> u64 {
        value.to_bits()
    }

    #[test]
    fn atof_port_scalbn_normal_and_zero() {
        assert_eq!(bits(scalbn(1.0, 0)), bits(1.0));
        assert_eq!(bits(scalbn(1.0, 1)), bits(2.0));
        assert_eq!(bits(scalbn(1.0, -1)), bits(0.5));
        assert_eq!(bits(scalbn(-1.0, 3)), bits(-8.0));
        assert_eq!(bits(scalbn(0.0, 100)), bits(0.0));
        assert_eq!(bits(scalbn(-0.0, -100)), bits(-0.0));
    }

    #[test]
    fn atof_port_scalbn_subnormal_boundaries() {
        assert_eq!(bits(scalbn(1.0, -1074)), 0x0000_0000_0000_0001);
        assert_eq!(bits(scalbn(1.0, -1075)), 0x0000_0000_0000_0000);
        assert_eq!(bits(scalbn(1.0, -1022)), 0x0010_0000_0000_0000);
        assert_eq!(bits(scalbn(1.0, -1023)), 0x0008_0000_0000_0000);
    }

    #[test]
    fn atof_port_scalbn_extremes_and_infinities() {
        assert_eq!(bits(scalbn(1.0, 1023)), 0x7FE0_0000_0000_0000);
        assert_eq!(bits(scalbn(1.0, 2000)), bits(f64::INFINITY));
        assert_eq!(bits(scalbn(1.0, -2000)), 0);
        assert_eq!(bits(scalbn(f64::INFINITY, 0)), bits(f64::INFINITY));
        assert_eq!(bits(scalbn(1.0, 100000)), bits(f64::INFINITY));
    }

    #[test]
    fn atof_port_scalbn_rounding_boundaries() {
        // 0x1.fffffffffffffp1023 * 2 -> inf (overflow), and exact max finite.
        assert_eq!(bits(scalbn(f64::MAX, 0)), bits(f64::MAX));
        // 2^-1022 scaled by 2^-52 is the least subnormal; by 2^-53 it is zero.
        assert_eq!(
            bits(scalbn(f64::from_bits(0x0010_0000_0000_0000), -52)),
            0x0000_0000_0000_0001
        );
        assert_eq!(
            bits(scalbn(f64::from_bits(0x0010_0000_0000_0000), -53)),
            0x0000_0000_0000_0000
        );
    }

    #[test]
    fn atof_port_fmod_matches_semantics() {
        assert_eq!(bits(fmod(5.5, 2.0)), bits(1.5));
        assert_eq!(bits(fmod(-5.5, 2.0)), bits(-1.5));
        assert_eq!(bits(fmod(4.0, 2.0)), bits(0.0));
        assert_eq!(bits(fmod(-4.0, 2.0)), bits(-0.0));
        assert_eq!(bits(fmod(1.5, f64::INFINITY)), bits(1.5));
        assert!(fmod(1.0, 0.0).is_nan());
    }
}
