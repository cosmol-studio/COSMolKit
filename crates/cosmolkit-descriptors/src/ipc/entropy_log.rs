//! Binary64 logarithm for the nonzero probabilities in Ipc's entropy sum.
//!
//! Adapted from musl v1.2.5, src/math/{log.c,log_data.c}, at
//! 0784374d561435f7c787a555aeab8ede699ed298. Copyright (c) 2018, Arm Limited.
//! SPDX-License-Identifier: MIT
//! Explicit fused operations select the rounding of the pinned RDKit reference
//! on every target, including WASM targets without a hardware FMA.

pub(super) fn probability_ln(x: f64) -> f64 {
    // musl✔️✔️: if (ix * 2 == 0)
    // musl✔️✔️:     return __math_divzero(1);
    // musl✔️✔️: if (ix == asuint64(INFINITY)) /* log(inf) == inf.  */
    // musl✔️✔️:     return x;
    // musl✔️✔️: if ((top & 0x8000) || top >= 0x7ff0)
    // musl✔️✔️:     return __math_invalid(x);
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return x;
    }
    if x.is_nan() || x < 0.0 {
        let invalid = x - x;
        return invalid / invalid;
    }
    // musl✔️✔️: #define LO asuint64(1.0 - 0x1p-4)
    // musl✔️✔️: if (predict_false(ix - LO < HI - LO)) {
    if x >= 1.0 - 1.0 / 16.0 {
        return near_one(x);
    }

    // musl✔️✔️: ix = asuint64(x);
    // musl✔️✔️: ix = asuint64(x * 0x1p52);
    // musl✔️✔️: ix -= 52ULL << 52;
    let ix = if x < f64::MIN_POSITIVE {
        (x * ((1_u64 << 52) as f64))
            .to_bits()
            .wrapping_sub(52_u64 << 52)
    } else {
        x.to_bits()
    };

    // musl✔️✔️: #define OFF 0x3fe6000000000000
    // musl✔️✔️: tmp = ix - OFF;
    // musl✔️✔️: i = (tmp >> (52 - LOG_TABLE_BITS)) % N;
    // musl✔️✔️: k = (int64_t)tmp >> 52; /* arithmetic shift */
    // musl✔️✔️: iz = ix - (tmp & 0xfffULL << 52);
    // musl✔️✔️: invc = T[i].invc;
    // musl✔️✔️: logc = T[i].logc;
    // musl✔️✔️: z = asdouble(iz);
    let tmp = ix.wrapping_sub(0x3fe6_0000_0000_0000);
    let index = ((tmp >> (52 - 7)) & 127) as usize;
    let k = ((tmp as i64) >> 52) as f64;
    let z = f64::from_bits(ix.wrapping_sub(tmp & (0xfff_u64 << 52)));
    let (invc, logc) = TABLE[index];

    // musl✔️✔️: r = __builtin_fma(z, invc, -1.0);
    // musl✔️✔️: w = kd * Ln2hi + logc;
    // musl✔️✔️: hi = w + r;
    // musl✔️✔️: lo = w - hi + r + kd * Ln2lo;
    let r = z.mul_add(invc, -1.0);
    let w = k.mul_add(LN2_HI, logc);
    let hi = w + r;
    let lo = k.mul_add(LN2_LO, (w - hi) + r);

    // musl✔️✔️: r2 = r * r; /* rounding error: 0x1p-54/N^2.  */
    // musl✔️✔️: y = lo + r2 * A[0] +
    // musl✔️✔️:     r * r2 * (A[1] + r * A[2] + r2 * (A[3] + r * A[4])) + hi;
    let r2 = r * r;
    let poly = r2.mul_add(r.mul_add(A[4], A[3]), r.mul_add(A[2], A[1]));
    (r * r2).mul_add(poly, r2.mul_add(A[0], lo)) + hi
}

fn near_one(x: f64) -> f64 {
    // musl✔️✔️: if (WANT_ROUNDING && predict_false(ix == asuint64(1.0)))
    // musl✔️✔️:     return 0;
    if x == 1.0 {
        return 0.0;
    }
    // musl✔️✔️: r = x - 1.0;
    // musl✔️✔️: r2 = r * r;
    // musl✔️✔️: r3 = r * r2;
    // musl✔️✔️: y = r3 *
    // musl✔️✔️:     (B[1] + r * B[2] + r2 * B[3] +
    // musl✔️✔️:      r3 * (B[4] + r * B[5] + r2 * B[6] +
    // musl✔️✔️:            r3 * (B[7] + r * B[8] + r2 * B[9] + r3 * B[10])));
    let r = x - 1.0;
    let r2 = r * r;
    let r3 = r * r2;
    let p1 = r2.mul_add(B[3], r.mul_add(B[2], B[1]));
    let p2 = r2.mul_add(B[6], r.mul_add(B[5], B[4]));
    let p3 = r2.mul_add(B[9], r.mul_add(B[8], B[7]));
    let tail = r3.mul_add(r3.mul_add(r3.mul_add(B[10], p3), p2), p1);

    // musl✔️✔️: w = r * 0x1p27;
    // musl✔️✔️: double_t rhi = r + w - w;
    // musl✔️✔️: double_t rlo = r - rhi;
    // musl✔️✔️: w = rhi * rhi * B[0]; /* B[0] == -0.5.  */
    // musl✔️✔️: hi = r + w;
    // musl✔️✔️: lo = r - hi + w;
    // musl✔️✔️: lo += B[0] * rlo * (rhi + r);
    // musl✔️✔️: y += lo;
    // musl✔️✔️: y += hi;
    let w = r * ((1_u32 << 27) as f64);
    let rhi = (r + w) - w;
    let rlo = r - rhi;
    let square = rhi * rhi;
    let hi = square.mul_add(B[0], r);
    let lo = (B[0] * rlo).mul_add(rhi + r, square.mul_add(B[0], r - hi));
    r3.mul_add(tail, lo) + hi
}

#[cfg(test)]
mod tests {
    use super::probability_ln;

    #[test]
    fn probability_log_matches_the_reference_at_rounding_boundaries() {
        // The two probabilities on either side of 15/16 exercise different
        // reductions. The last two distinguish fused from unfused evaluation.
        for (x, logarithm) in [
            (1_u64, 0xc087_4385_446d_71c3),
            (1.0_f64.to_bits(), 0),
            (0x3fed_ffff_ffff_ffff, 0xbfb0_8598_b59e_3a0f),
            (0x3fee_0000_0000_0000, 0xbfb0_8598_b59e_3a07),
            (0x3fe9_718a_85f7_e4c4, 0xbfcd_58e0_83de_7cbe),
            (0x3fee_0877_aaa6_1fe9, 0xbfb0_3d61_7a32_28ae),
        ] {
            assert_eq!(probability_ln(f64::from_bits(x)).to_bits(), logarithm);
        }
    }
}

// musl✔️✔️: .ln2hi = 0x1.62e42fefa3800p-1,
const LN2_HI: f64 = f64::from_bits(0x3fe62e42fefa3800);
// musl✔️✔️: .ln2lo = 0x1.ef35793c76730p-45,
const LN2_LO: f64 = f64::from_bits(0x3d2ef35793c76730);

#[rustfmt::skip]
const A: [f64; 5] = [
    // musl✔️✔️: -0x1.0000000000001p-1,
    f64::from_bits(0xbfe0000000000001),
    // musl✔️✔️: 0x1.555555551305bp-2,
    f64::from_bits(0x3fd555555551305b),
    // musl✔️✔️: -0x1.fffffffeb459p-3,
    f64::from_bits(0xbfcfffffffeb4590),
    // musl✔️✔️: 0x1.999b324f10111p-3,
    f64::from_bits(0x3fc999b324f10111),
    // musl✔️✔️: -0x1.55575e506c89fp-3,
    f64::from_bits(0xbfc55575e506c89f),
];

#[rustfmt::skip]
const B: [f64; 11] = [
    // musl✔️✔️: -0x1p-1,
    f64::from_bits(0xbfe0000000000000),
    // musl✔️✔️: 0x1.5555555555577p-2,
    f64::from_bits(0x3fd5555555555577),
    // musl✔️✔️: -0x1.ffffffffffdcbp-3,
    f64::from_bits(0xbfcffffffffffdcb),
    // musl✔️✔️: 0x1.999999995dd0cp-3,
    f64::from_bits(0x3fc999999995dd0c),
    // musl✔️✔️: -0x1.55555556745a7p-3,
    f64::from_bits(0xbfc55555556745a7),
    // musl✔️✔️: 0x1.24924a344de3p-3,
    f64::from_bits(0x3fc24924a344de30),
    // musl✔️✔️: -0x1.fffffa4423d65p-4,
    f64::from_bits(0xbfbfffffa4423d65),
    // musl✔️✔️: 0x1.c7184282ad6cap-4,
    f64::from_bits(0x3fbc7184282ad6ca),
    // musl✔️✔️: -0x1.999eb43b068ffp-4,
    f64::from_bits(0xbfb999eb43b068ff),
    // musl✔️✔️: 0x1.78182f7afd085p-4,
    f64::from_bits(0x3fb78182f7afd085),
    // musl✔️✔️: -0x1.5521375d145cdp-4,
    f64::from_bits(0xbfb5521375d145cd),
];

#[rustfmt::skip]
const TABLE: [(f64, f64); 128] = [
    // musl✔️✔️: {0x1.734f0c3e0de9fp+0, -0x1.7cc7f79e69000p-2},
    (f64::from_bits(0x3ff734f0c3e0de9f), f64::from_bits(0xbfd7cc7f79e69000)),
    // musl✔️✔️: {0x1.713786a2ce91fp+0, -0x1.76feec20d0000p-2},
    (f64::from_bits(0x3ff713786a2ce91f), f64::from_bits(0xbfd76feec20d0000)),
    // musl✔️✔️: {0x1.6f26008fab5a0p+0, -0x1.713e31351e000p-2},
    (f64::from_bits(0x3ff6f26008fab5a0), f64::from_bits(0xbfd713e31351e000)),
    // musl✔️✔️: {0x1.6d1a61f138c7dp+0, -0x1.6b85b38287800p-2},
    (f64::from_bits(0x3ff6d1a61f138c7d), f64::from_bits(0xbfd6b85b38287800)),
    // musl✔️✔️: {0x1.6b1490bc5b4d1p+0, -0x1.65d5590807800p-2},
    (f64::from_bits(0x3ff6b1490bc5b4d1), f64::from_bits(0xbfd65d5590807800)),
    // musl✔️✔️: {0x1.69147332f0cbap+0, -0x1.602d076180000p-2},
    (f64::from_bits(0x3ff69147332f0cba), f64::from_bits(0xbfd602d076180000)),
    // musl✔️✔️: {0x1.6719f18224223p+0, -0x1.5a8ca86909000p-2},
    (f64::from_bits(0x3ff6719f18224223), f64::from_bits(0xbfd5a8ca86909000)),
    // musl✔️✔️: {0x1.6524f99a51ed9p+0, -0x1.54f4356035000p-2},
    (f64::from_bits(0x3ff6524f99a51ed9), f64::from_bits(0xbfd54f4356035000)),
    // musl✔️✔️: {0x1.63356aa8f24c4p+0, -0x1.4f637c36b4000p-2},
    (f64::from_bits(0x3ff63356aa8f24c4), f64::from_bits(0xbfd4f637c36b4000)),
    // musl✔️✔️: {0x1.614b36b9ddc14p+0, -0x1.49da7fda85000p-2},
    (f64::from_bits(0x3ff614b36b9ddc14), f64::from_bits(0xbfd49da7fda85000)),
    // musl✔️✔️: {0x1.5f66452c65c4cp+0, -0x1.445923989a800p-2},
    (f64::from_bits(0x3ff5f66452c65c4c), f64::from_bits(0xbfd445923989a800)),
    // musl✔️✔️: {0x1.5d867b5912c4fp+0, -0x1.3edf439b0b800p-2},
    (f64::from_bits(0x3ff5d867b5912c4f), f64::from_bits(0xbfd3edf439b0b800)),
    // musl✔️✔️: {0x1.5babccb5b90dep+0, -0x1.396ce448f7000p-2},
    (f64::from_bits(0x3ff5babccb5b90de), f64::from_bits(0xbfd396ce448f7000)),
    // musl✔️✔️: {0x1.59d61f2d91a78p+0, -0x1.3401e17bda000p-2},
    (f64::from_bits(0x3ff59d61f2d91a78), f64::from_bits(0xbfd3401e17bda000)),
    // musl✔️✔️: {0x1.5805612465687p+0, -0x1.2e9e2ef468000p-2},
    (f64::from_bits(0x3ff5805612465687), f64::from_bits(0xbfd2e9e2ef468000)),
    // musl✔️✔️: {0x1.56397cee76bd3p+0, -0x1.2941b3830e000p-2},
    (f64::from_bits(0x3ff56397cee76bd3), f64::from_bits(0xbfd2941b3830e000)),
    // musl✔️✔️: {0x1.54725e2a77f93p+0, -0x1.23ec58cda8800p-2},
    (f64::from_bits(0x3ff54725e2a77f93), f64::from_bits(0xbfd23ec58cda8800)),
    // musl✔️✔️: {0x1.52aff42064583p+0, -0x1.1e9e129279000p-2},
    (f64::from_bits(0x3ff52aff42064583), f64::from_bits(0xbfd1e9e129279000)),
    // musl✔️✔️: {0x1.50f22dbb2bddfp+0, -0x1.1956d2b48f800p-2},
    (f64::from_bits(0x3ff50f22dbb2bddf), f64::from_bits(0xbfd1956d2b48f800)),
    // musl✔️✔️: {0x1.4f38f4734ded7p+0, -0x1.141679ab9f800p-2},
    (f64::from_bits(0x3ff4f38f4734ded7), f64::from_bits(0xbfd141679ab9f800)),
    // musl✔️✔️: {0x1.4d843cfde2840p+0, -0x1.0edd094ef9800p-2},
    (f64::from_bits(0x3ff4d843cfde2840), f64::from_bits(0xbfd0edd094ef9800)),
    // musl✔️✔️: {0x1.4bd3ec078a3c8p+0, -0x1.09aa518db1000p-2},
    (f64::from_bits(0x3ff4bd3ec078a3c8), f64::from_bits(0xbfd09aa518db1000)),
    // musl✔️✔️: {0x1.4a27fc3e0258ap+0, -0x1.047e65263b800p-2},
    (f64::from_bits(0x3ff4a27fc3e0258a), f64::from_bits(0xbfd047e65263b800)),
    // musl✔️✔️: {0x1.4880524d48434p+0, -0x1.feb224586f000p-3},
    (f64::from_bits(0x3ff4880524d48434), f64::from_bits(0xbfcfeb224586f000)),
    // musl✔️✔️: {0x1.46dce1b192d0bp+0, -0x1.f474a7517b000p-3},
    (f64::from_bits(0x3ff46dce1b192d0b), f64::from_bits(0xbfcf474a7517b000)),
    // musl✔️✔️: {0x1.453d9d3391854p+0, -0x1.ea4443d103000p-3},
    (f64::from_bits(0x3ff453d9d3391854), f64::from_bits(0xbfcea4443d103000)),
    // musl✔️✔️: {0x1.43a2744b4845ap+0, -0x1.e020d44e9b000p-3},
    (f64::from_bits(0x3ff43a2744b4845a), f64::from_bits(0xbfce020d44e9b000)),
    // musl✔️✔️: {0x1.420b54115f8fbp+0, -0x1.d60a22977f000p-3},
    (f64::from_bits(0x3ff420b54115f8fb), f64::from_bits(0xbfcd60a22977f000)),
    // musl✔️✔️: {0x1.40782da3ef4b1p+0, -0x1.cc00104959000p-3},
    (f64::from_bits(0x3ff40782da3ef4b1), f64::from_bits(0xbfccc00104959000)),
    // musl✔️✔️: {0x1.3ee8f5d57fe8fp+0, -0x1.c202956891000p-3},
    (f64::from_bits(0x3ff3ee8f5d57fe8f), f64::from_bits(0xbfcc202956891000)),
    // musl✔️✔️: {0x1.3d5d9a00b4ce9p+0, -0x1.b81178d811000p-3},
    (f64::from_bits(0x3ff3d5d9a00b4ce9), f64::from_bits(0xbfcb81178d811000)),
    // musl✔️✔️: {0x1.3bd60c010c12bp+0, -0x1.ae2c9ccd3d000p-3},
    (f64::from_bits(0x3ff3bd60c010c12b), f64::from_bits(0xbfcae2c9ccd3d000)),
    // musl✔️✔️: {0x1.3a5242b75dab8p+0, -0x1.a45402e129000p-3},
    (f64::from_bits(0x3ff3a5242b75dab8), f64::from_bits(0xbfca45402e129000)),
    // musl✔️✔️: {0x1.38d22cd9fd002p+0, -0x1.9a877681df000p-3},
    (f64::from_bits(0x3ff38d22cd9fd002), f64::from_bits(0xbfc9a877681df000)),
    // musl✔️✔️: {0x1.3755bc5847a1cp+0, -0x1.90c6d69483000p-3},
    (f64::from_bits(0x3ff3755bc5847a1c), f64::from_bits(0xbfc90c6d69483000)),
    // musl✔️✔️: {0x1.35dce49ad36e2p+0, -0x1.87120a645c000p-3},
    (f64::from_bits(0x3ff35dce49ad36e2), f64::from_bits(0xbfc87120a645c000)),
    // musl✔️✔️: {0x1.34679984dd440p+0, -0x1.7d68fb4143000p-3},
    (f64::from_bits(0x3ff34679984dd440), f64::from_bits(0xbfc7d68fb4143000)),
    // musl✔️✔️: {0x1.32f5cceffcb24p+0, -0x1.73cb83c627000p-3},
    (f64::from_bits(0x3ff32f5cceffcb24), f64::from_bits(0xbfc73cb83c627000)),
    // musl✔️✔️: {0x1.3187775a10d49p+0, -0x1.6a39a9b376000p-3},
    (f64::from_bits(0x3ff3187775a10d49), f64::from_bits(0xbfc6a39a9b376000)),
    // musl✔️✔️: {0x1.301c8373e3990p+0, -0x1.60b3154b7a000p-3},
    (f64::from_bits(0x3ff301c8373e3990), f64::from_bits(0xbfc60b3154b7a000)),
    // musl✔️✔️: {0x1.2eb4ebb95f841p+0, -0x1.5737d76243000p-3},
    (f64::from_bits(0x3ff2eb4ebb95f841), f64::from_bits(0xbfc5737d76243000)),
    // musl✔️✔️: {0x1.2d50a0219a9d1p+0, -0x1.4dc7b8fc23000p-3},
    (f64::from_bits(0x3ff2d50a0219a9d1), f64::from_bits(0xbfc4dc7b8fc23000)),
    // musl✔️✔️: {0x1.2bef9a8b7fd2ap+0, -0x1.4462c51d20000p-3},
    (f64::from_bits(0x3ff2bef9a8b7fd2a), f64::from_bits(0xbfc4462c51d20000)),
    // musl✔️✔️: {0x1.2a91c7a0c1babp+0, -0x1.3b08abc830000p-3},
    (f64::from_bits(0x3ff2a91c7a0c1bab), f64::from_bits(0xbfc3b08abc830000)),
    // musl✔️✔️: {0x1.293726014b530p+0, -0x1.31b996b490000p-3},
    (f64::from_bits(0x3ff293726014b530), f64::from_bits(0xbfc31b996b490000)),
    // musl✔️✔️: {0x1.27dfa5757a1f5p+0, -0x1.2875490a44000p-3},
    (f64::from_bits(0x3ff27dfa5757a1f5), f64::from_bits(0xbfc2875490a44000)),
    // musl✔️✔️: {0x1.268b39b1d3bbfp+0, -0x1.1f3b9f879a000p-3},
    (f64::from_bits(0x3ff268b39b1d3bbf), f64::from_bits(0xbfc1f3b9f879a000)),
    // musl✔️✔️: {0x1.2539d838ff5bdp+0, -0x1.160c8252ca000p-3},
    (f64::from_bits(0x3ff2539d838ff5bd), f64::from_bits(0xbfc160c8252ca000)),
    // musl✔️✔️: {0x1.23eb7aac9083bp+0, -0x1.0ce7f57f72000p-3},
    (f64::from_bits(0x3ff23eb7aac9083b), f64::from_bits(0xbfc0ce7f57f72000)),
    // musl✔️✔️: {0x1.22a012ba940b6p+0, -0x1.03cdc49fea000p-3},
    (f64::from_bits(0x3ff22a012ba940b6), f64::from_bits(0xbfc03cdc49fea000)),
    // musl✔️✔️: {0x1.2157996cc4132p+0, -0x1.f57bdbc4b8000p-4},
    (f64::from_bits(0x3ff2157996cc4132), f64::from_bits(0xbfbf57bdbc4b8000)),
    // musl✔️✔️: {0x1.201201dd2fc9bp+0, -0x1.e370896404000p-4},
    (f64::from_bits(0x3ff201201dd2fc9b), f64::from_bits(0xbfbe370896404000)),
    // musl✔️✔️: {0x1.1ecf4494d480bp+0, -0x1.d17983ef94000p-4},
    (f64::from_bits(0x3ff1ecf4494d480b), f64::from_bits(0xbfbd17983ef94000)),
    // musl✔️✔️: {0x1.1d8f5528f6569p+0, -0x1.bf9674ed8a000p-4},
    (f64::from_bits(0x3ff1d8f5528f6569), f64::from_bits(0xbfbbf9674ed8a000)),
    // musl✔️✔️: {0x1.1c52311577e7cp+0, -0x1.adc79202f6000p-4},
    (f64::from_bits(0x3ff1c52311577e7c), f64::from_bits(0xbfbadc79202f6000)),
    // musl✔️✔️: {0x1.1b17c74cb26e9p+0, -0x1.9c0c3e7288000p-4},
    (f64::from_bits(0x3ff1b17c74cb26e9), f64::from_bits(0xbfb9c0c3e7288000)),
    // musl✔️✔️: {0x1.19e010c2c1ab6p+0, -0x1.8a646b372c000p-4},
    (f64::from_bits(0x3ff19e010c2c1ab6), f64::from_bits(0xbfb8a646b372c000)),
    // musl✔️✔️: {0x1.18ab07bb670bdp+0, -0x1.78d01b3ac0000p-4},
    (f64::from_bits(0x3ff18ab07bb670bd), f64::from_bits(0xbfb78d01b3ac0000)),
    // musl✔️✔️: {0x1.1778a25efbcb6p+0, -0x1.674f145380000p-4},
    (f64::from_bits(0x3ff1778a25efbcb6), f64::from_bits(0xbfb674f145380000)),
    // musl✔️✔️: {0x1.1648d354c31dap+0, -0x1.55e0e6d878000p-4},
    (f64::from_bits(0x3ff1648d354c31da), f64::from_bits(0xbfb55e0e6d878000)),
    // musl✔️✔️: {0x1.151b990275fddp+0, -0x1.4485cdea1e000p-4},
    (f64::from_bits(0x3ff151b990275fdd), f64::from_bits(0xbfb4485cdea1e000)),
    // musl✔️✔️: {0x1.13f0ea432d24cp+0, -0x1.333d94d6aa000p-4},
    (f64::from_bits(0x3ff13f0ea432d24c), f64::from_bits(0xbfb333d94d6aa000)),
    // musl✔️✔️: {0x1.12c8b7210f9dap+0, -0x1.22079f8c56000p-4},
    (f64::from_bits(0x3ff12c8b7210f9da), f64::from_bits(0xbfb22079f8c56000)),
    // musl✔️✔️: {0x1.11a3028ecb531p+0, -0x1.10e4698622000p-4},
    (f64::from_bits(0x3ff11a3028ecb531), f64::from_bits(0xbfb10e4698622000)),
    // musl✔️✔️: {0x1.107fbda8434afp+0, -0x1.ffa6c6ad20000p-5},
    (f64::from_bits(0x3ff107fbda8434af), f64::from_bits(0xbfaffa6c6ad20000)),
    // musl✔️✔️: {0x1.0f5ee0f4e6bb3p+0, -0x1.dda8d4a774000p-5},
    (f64::from_bits(0x3ff0f5ee0f4e6bb3), f64::from_bits(0xbfadda8d4a774000)),
    // musl✔️✔️: {0x1.0e4065d2a9fcep+0, -0x1.bbcece4850000p-5},
    (f64::from_bits(0x3ff0e4065d2a9fce), f64::from_bits(0xbfabbcece4850000)),
    // musl✔️✔️: {0x1.0d244632ca521p+0, -0x1.9a1894012c000p-5},
    (f64::from_bits(0x3ff0d244632ca521), f64::from_bits(0xbfa9a1894012c000)),
    // musl✔️✔️: {0x1.0c0a77ce2981ap+0, -0x1.788583302c000p-5},
    (f64::from_bits(0x3ff0c0a77ce2981a), f64::from_bits(0xbfa788583302c000)),
    // musl✔️✔️: {0x1.0af2f83c636d1p+0, -0x1.5715e67d68000p-5},
    (f64::from_bits(0x3ff0af2f83c636d1), f64::from_bits(0xbfa5715e67d68000)),
    // musl✔️✔️: {0x1.09ddb98a01339p+0, -0x1.35c8a49658000p-5},
    (f64::from_bits(0x3ff09ddb98a01339), f64::from_bits(0xbfa35c8a49658000)),
    // musl✔️✔️: {0x1.08cabaf52e7dfp+0, -0x1.149e364154000p-5},
    (f64::from_bits(0x3ff08cabaf52e7df), f64::from_bits(0xbfa149e364154000)),
    // musl✔️✔️: {0x1.07b9f2f4e28fbp+0, -0x1.e72c082eb8000p-6},
    (f64::from_bits(0x3ff07b9f2f4e28fb), f64::from_bits(0xbf9e72c082eb8000)),
    // musl✔️✔️: {0x1.06ab58c358f19p+0, -0x1.a55f152528000p-6},
    (f64::from_bits(0x3ff06ab58c358f19), f64::from_bits(0xbf9a55f152528000)),
    // musl✔️✔️: {0x1.059eea5ecf92cp+0, -0x1.63d62cf818000p-6},
    (f64::from_bits(0x3ff059eea5ecf92c), f64::from_bits(0xbf963d62cf818000)),
    // musl✔️✔️: {0x1.04949cdd12c90p+0, -0x1.228fb8caa0000p-6},
    (f64::from_bits(0x3ff04949cdd12c90), f64::from_bits(0xbf9228fb8caa0000)),
    // musl✔️✔️: {0x1.038c6c6f0ada9p+0, -0x1.c317b20f90000p-7},
    (f64::from_bits(0x3ff038c6c6f0ada9), f64::from_bits(0xbf8c317b20f90000)),
    // musl✔️✔️: {0x1.02865137932a9p+0, -0x1.419355daa0000p-7},
    (f64::from_bits(0x3ff02865137932a9), f64::from_bits(0xbf8419355daa0000)),
    // musl✔️✔️: {0x1.0182427ea7348p+0, -0x1.81203c2ec0000p-8},
    (f64::from_bits(0x3ff0182427ea7348), f64::from_bits(0xbf781203c2ec0000)),
    // musl✔️✔️: {0x1.008040614b195p+0, -0x1.0040979240000p-9},
    (f64::from_bits(0x3ff008040614b195), f64::from_bits(0xbf60040979240000)),
    // musl✔️✔️: {0x1.fe01ff726fa1ap-1, 0x1.feff384900000p-9},
    (f64::from_bits(0x3fefe01ff726fa1a), f64::from_bits(0x3f6feff384900000)),
    // musl✔️✔️: {0x1.fa11cc261ea74p-1, 0x1.7dc41353d0000p-7},
    (f64::from_bits(0x3fefa11cc261ea74), f64::from_bits(0x3f87dc41353d0000)),
    // musl✔️✔️: {0x1.f6310b081992ep-1, 0x1.3cea3c4c28000p-6},
    (f64::from_bits(0x3fef6310b081992e), f64::from_bits(0x3f93cea3c4c28000)),
    // musl✔️✔️: {0x1.f25f63ceeadcdp-1, 0x1.b9fc114890000p-6},
    (f64::from_bits(0x3fef25f63ceeadcd), f64::from_bits(0x3f9b9fc114890000)),
    // musl✔️✔️: {0x1.ee9c8039113e7p-1, 0x1.1b0d8ce110000p-5},
    (f64::from_bits(0x3feee9c8039113e7), f64::from_bits(0x3fa1b0d8ce110000)),
    // musl✔️✔️: {0x1.eae8078cbb1abp-1, 0x1.58a5bd001c000p-5},
    (f64::from_bits(0x3feeae8078cbb1ab), f64::from_bits(0x3fa58a5bd001c000)),
    // musl✔️✔️: {0x1.e741aa29d0c9bp-1, 0x1.95c8340d88000p-5},
    (f64::from_bits(0x3fee741aa29d0c9b), f64::from_bits(0x3fa95c8340d88000)),
    // musl✔️✔️: {0x1.e3a91830a99b5p-1, 0x1.d276aef578000p-5},
    (f64::from_bits(0x3fee3a91830a99b5), f64::from_bits(0x3fad276aef578000)),
    // musl✔️✔️: {0x1.e01e009609a56p-1, 0x1.07598e598c000p-4},
    (f64::from_bits(0x3fee01e009609a56), f64::from_bits(0x3fb07598e598c000)),
    // musl✔️✔️: {0x1.dca01e577bb98p-1, 0x1.253f5e30d2000p-4},
    (f64::from_bits(0x3fedca01e577bb98), f64::from_bits(0x3fb253f5e30d2000)),
    // musl✔️✔️: {0x1.d92f20b7c9103p-1, 0x1.42edd8b380000p-4},
    (f64::from_bits(0x3fed92f20b7c9103), f64::from_bits(0x3fb42edd8b380000)),
    // musl✔️✔️: {0x1.d5cac66fb5ccep-1, 0x1.606598757c000p-4},
    (f64::from_bits(0x3fed5cac66fb5cce), f64::from_bits(0x3fb606598757c000)),
    // musl✔️✔️: {0x1.d272caa5ede9dp-1, 0x1.7da76356a0000p-4},
    (f64::from_bits(0x3fed272caa5ede9d), f64::from_bits(0x3fb7da76356a0000)),
    // musl✔️✔️: {0x1.cf26e3e6b2ccdp-1, 0x1.9ab434e1c6000p-4},
    (f64::from_bits(0x3fecf26e3e6b2ccd), f64::from_bits(0x3fb9ab434e1c6000)),
    // musl✔️✔️: {0x1.cbe6da2a77902p-1, 0x1.b78c7bb0d6000p-4},
    (f64::from_bits(0x3fecbe6da2a77902), f64::from_bits(0x3fbb78c7bb0d6000)),
    // musl✔️✔️: {0x1.c8b266d37086dp-1, 0x1.d431332e72000p-4},
    (f64::from_bits(0x3fec8b266d37086d), f64::from_bits(0x3fbd431332e72000)),
    // musl✔️✔️: {0x1.c5894bd5d5804p-1, 0x1.f0a3171de6000p-4},
    (f64::from_bits(0x3fec5894bd5d5804), f64::from_bits(0x3fbf0a3171de6000)),
    // musl✔️✔️: {0x1.c26b533bb9f8cp-1, 0x1.067152b914000p-3},
    (f64::from_bits(0x3fec26b533bb9f8c), f64::from_bits(0x3fc067152b914000)),
    // musl✔️✔️: {0x1.bf583eeece73fp-1, 0x1.147858292b000p-3},
    (f64::from_bits(0x3febf583eeece73f), f64::from_bits(0x3fc147858292b000)),
    // musl✔️✔️: {0x1.bc4fd75db96c1p-1, 0x1.2266ecdca3000p-3},
    (f64::from_bits(0x3febc4fd75db96c1), f64::from_bits(0x3fc2266ecdca3000)),
    // musl✔️✔️: {0x1.b951e0c864a28p-1, 0x1.303d7a6c55000p-3},
    (f64::from_bits(0x3feb951e0c864a28), f64::from_bits(0x3fc303d7a6c55000)),
    // musl✔️✔️: {0x1.b65e2c5ef3e2cp-1, 0x1.3dfc33c331000p-3},
    (f64::from_bits(0x3feb65e2c5ef3e2c), f64::from_bits(0x3fc3dfc33c331000)),
    // musl✔️✔️: {0x1.b374867c9888bp-1, 0x1.4ba366b7a8000p-3},
    (f64::from_bits(0x3feb374867c9888b), f64::from_bits(0x3fc4ba366b7a8000)),
    // musl✔️✔️: {0x1.b094b211d304ap-1, 0x1.5933928d1f000p-3},
    (f64::from_bits(0x3feb094b211d304a), f64::from_bits(0x3fc5933928d1f000)),
    // musl✔️✔️: {0x1.adbe885f2ef7ep-1, 0x1.66acd2418f000p-3},
    (f64::from_bits(0x3feadbe885f2ef7e), f64::from_bits(0x3fc66acd2418f000)),
    // musl✔️✔️: {0x1.aaf1d31603da2p-1, 0x1.740f8ec669000p-3},
    (f64::from_bits(0x3feaaf1d31603da2), f64::from_bits(0x3fc740f8ec669000)),
    // musl✔️✔️: {0x1.a82e63fd358a7p-1, 0x1.815c0f51af000p-3},
    (f64::from_bits(0x3fea82e63fd358a7), f64::from_bits(0x3fc815c0f51af000)),
    // musl✔️✔️: {0x1.a5740ef09738bp-1, 0x1.8e92954f68000p-3},
    (f64::from_bits(0x3fea5740ef09738b), f64::from_bits(0x3fc8e92954f68000)),
    // musl✔️✔️: {0x1.a2c2a90ab4b27p-1, 0x1.9bb3602f84000p-3},
    (f64::from_bits(0x3fea2c2a90ab4b27), f64::from_bits(0x3fc9bb3602f84000)),
    // musl✔️✔️: {0x1.a01a01393f2d1p-1, 0x1.a8bed1c2c0000p-3},
    (f64::from_bits(0x3fea01a01393f2d1), f64::from_bits(0x3fca8bed1c2c0000)),
    // musl✔️✔️: {0x1.9d79f24db3c1bp-1, 0x1.b5b515c01d000p-3},
    (f64::from_bits(0x3fe9d79f24db3c1b), f64::from_bits(0x3fcb5b515c01d000)),
    // musl✔️✔️: {0x1.9ae2505c7b190p-1, 0x1.c2967ccbcc000p-3},
    (f64::from_bits(0x3fe9ae2505c7b190), f64::from_bits(0x3fcc2967ccbcc000)),
    // musl✔️✔️: {0x1.9852ef297ce2fp-1, 0x1.cf635d5486000p-3},
    (f64::from_bits(0x3fe9852ef297ce2f), f64::from_bits(0x3fccf635d5486000)),
    // musl✔️✔️: {0x1.95cbaeea44b75p-1, 0x1.dc1bd3446c000p-3},
    (f64::from_bits(0x3fe95cbaeea44b75), f64::from_bits(0x3fcdc1bd3446c000)),
    // musl✔️✔️: {0x1.934c69de74838p-1, 0x1.e8c01b8cfe000p-3},
    (f64::from_bits(0x3fe934c69de74838), f64::from_bits(0x3fce8c01b8cfe000)),
    // musl✔️✔️: {0x1.90d4f2f6752e6p-1, 0x1.f5509c0179000p-3},
    (f64::from_bits(0x3fe90d4f2f6752e6), f64::from_bits(0x3fcf5509c0179000)),
    // musl✔️✔️: {0x1.8e6528effd79dp-1, 0x1.00e6c121fb800p-2},
    (f64::from_bits(0x3fe8e6528effd79d), f64::from_bits(0x3fd00e6c121fb800)),
    // musl✔️✔️: {0x1.8bfce9fcc007cp-1, 0x1.071b80e93d000p-2},
    (f64::from_bits(0x3fe8bfce9fcc007c), f64::from_bits(0x3fd071b80e93d000)),
    // musl✔️✔️: {0x1.899c0dabec30ep-1, 0x1.0d46b9e867000p-2},
    (f64::from_bits(0x3fe899c0dabec30e), f64::from_bits(0x3fd0d46b9e867000)),
    // musl✔️✔️: {0x1.87427aa2317fbp-1, 0x1.13687334bd000p-2},
    (f64::from_bits(0x3fe87427aa2317fb), f64::from_bits(0x3fd13687334bd000)),
    // musl✔️✔️: {0x1.84f00acb39a08p-1, 0x1.1980d67234800p-2},
    (f64::from_bits(0x3fe84f00acb39a08), f64::from_bits(0x3fd1980d67234800)),
    // musl✔️✔️: {0x1.82a49e8653e55p-1, 0x1.1f8ffe0cc8000p-2},
    (f64::from_bits(0x3fe82a49e8653e55), f64::from_bits(0x3fd1f8ffe0cc8000)),
    // musl✔️✔️: {0x1.8060195f40260p-1, 0x1.2595fd7636800p-2},
    (f64::from_bits(0x3fe8060195f40260), f64::from_bits(0x3fd2595fd7636800)),
    // musl✔️✔️: {0x1.7e22563e0a329p-1, 0x1.2b9300914a800p-2},
    (f64::from_bits(0x3fe7e22563e0a329), f64::from_bits(0x3fd2b9300914a800)),
    // musl✔️✔️: {0x1.7beb377dcb5adp-1, 0x1.3187210436000p-2},
    (f64::from_bits(0x3fe7beb377dcb5ad), f64::from_bits(0x3fd3187210436000)),
    // musl✔️✔️: {0x1.79baa679725c2p-1, 0x1.377266dec1800p-2},
    (f64::from_bits(0x3fe79baa679725c2), f64::from_bits(0x3fd377266dec1800)),
    // musl✔️✔️: {0x1.77907f2170657p-1, 0x1.3d54ffbaf3000p-2},
    (f64::from_bits(0x3fe77907f2170657), f64::from_bits(0x3fd3d54ffbaf3000)),
    // musl✔️✔️: {0x1.756cadbd6130cp-1, 0x1.432eee32fe000p-2},
    (f64::from_bits(0x3fe756cadbd6130c), f64::from_bits(0x3fd432eee32fe000)),
];
