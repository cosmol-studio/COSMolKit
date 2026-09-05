//! Gemmi's locale-independent `%g` formatting path.

const BOT: [f64; 23] = [
    1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16,
    1e17, 1e18, 1e19, 1e20, 1e21, 1e22,
];
const NEG_BOT: [f64; 22] = [
    1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15,
    1e-16, 1e-17, 1e-18, 1e-19, 1e-20, 1e-21, 1e-22,
];
const NEG_BOT_ERR: [f64; 22] = [
    -5.551115123125783e-18,
    -2.0816681711721684e-19,
    -2.0816681711721686e-20,
    -4.79217360238593e-21,
    -8.18030539140313e-22,
    4.525188817411374e-23,
    4.525188817411374e-24,
    -2.092256083012847e-25,
    -6.228159145777985e-26,
    -3.6432197315497743e-27,
    6.050303071806019e-28,
    2.0113352370744385e-29,
    -3.037374556340037e-30,
    1.1806906454401013e-32,
    -7.770539987666108e-32,
    2.0902213275965398e-33,
    -7.154242405462192e-34,
    -7.154242405462193e-35,
    2.475407316473987e-36,
    5.484672854579043e-37,
    9.246254777210363e-38,
    -4.859677432657087e-39,
];
const TOP: [f64; 13] = [
    1e23, 1e46, 1e69, 1e92, 1e115, 1e138, 1e161, 1e184, 1e207, 1e230, 1e253, 1e276, 1e299,
];
const NEG_TOP: [f64; 13] = [
    1e-23, 1e-46, 1e-69, 1e-92, 1e-115, 1e-138, 1e-161, 1e-184, 1e-207, 1e-230, 1e-253, 1e-276,
    1e-299,
];
const TOP_ERR: [f64; 13] = [
    8_388_608.0,
    6.860180964052972e28,
    -7.253143638152922e52,
    -4.3377296974619174e75,
    -1.5559416129466825e98,
    -3.2841562489204913e121,
    -3.7745893248228135e144,
    -1.7356668416969134e167,
    -3.8893577551088374e190,
    -9.956644432600512e213,
    6.364129306223243e236,
    -5.206914080024981e259,
    -5.250476025520439e282,
];
const NEG_TOP_ERR: [f64; 13] = [
    3.9565301985100693e-40,
    -2.299904345391321e-63,
    3.6506201437945798e-86,
    1.1875228833981544e-109,
    -5.064490231692861e-132,
    -6.715683724786543e-155,
    -2.812077463003139e-178,
    -5.777891238658995e-201,
    7.499710055933453e-224,
    -4.643966891513449e-247,
    -6.369110076296214e-270,
    -9.436808465446358e-293,
    8.0970921678015e-317,
];
const POW_TEN: [u64; 20] = [
    1,
    10,
    100,
    1_000,
    10_000,
    100_000,
    1_000_000,
    10_000_000,
    100_000_000,
    1_000_000_000,
    10_000_000_000,
    100_000_000_000,
    1_000_000_000_000,
    10_000_000_000_000,
    100_000_000_000_000,
    1_000_000_000_000_000,
    10_000_000_000_000_000,
    100_000_000_000_000_000,
    1_000_000_000_000_000_000,
    10_000_000_000_000_000_000,
];

#[derive(Debug, Clone, PartialEq, Eq)]
enum RealDigits {
    Finite { digits: String, decimal_pos: i32 },
    Special(&'static str),
}

fn split_high(value: f64) -> f64 {
    f64::from_bits(value.to_bits() & (!0_u64 << 27))
}

// BEGIN GEMMI C FUNCTION stbsp__ddmulthi
// Gemmi✔️✔️: #define stbsp__ddmulthi(oh, ol, xh, yh)                            \
// Gemmi✔️✔️:    {                                                               \
// Gemmi✔️✔️:       double ahi = 0, alo, bhi = 0, blo;                           \
// Gemmi✔️✔️:       stbsp__int64 bt;                                             \
// Gemmi✔️✔️:       oh = xh * yh;                                                \
// Gemmi✔️✔️:       STBSP__COPYFP(bt, xh);                                       \
// Gemmi✔️✔️:       bt &= ((~(stbsp__uint64)0) << 27);                           \
// Gemmi✔️✔️:       STBSP__COPYFP(ahi, bt);                                      \
// Gemmi✔️✔️:       alo = xh - ahi;                                              \
// Gemmi✔️✔️:       STBSP__COPYFP(bt, yh);                                       \
// Gemmi✔️✔️:       bt &= ((~(stbsp__uint64)0) << 27);                           \
// Gemmi✔️✔️:       STBSP__COPYFP(bhi, bt);                                      \
// Gemmi✔️✔️:       blo = yh - bhi;                                              \
// Gemmi✔️✔️:       ol = ((ahi * bhi - oh) + ahi * blo + alo * bhi) + alo * blo; \
// Gemmi✔️✔️:    }
// END GEMMI C FUNCTION
fn double_double_multiply_high(xh: f64, yh: f64) -> (f64, f64) {
    let high = xh * yh;
    let ahi = split_high(xh);
    let alo = xh - ahi;
    let bhi = split_high(yh);
    let blo = yh - bhi;
    let low = ((ahi * bhi - high) + ahi * blo + alo * bhi) + alo * blo;
    (high, low)
}

fn renormalize(high: f64, low: f64) -> (f64, f64) {
    let sum = high + low;
    (sum, low - (sum - high))
}

// BEGIN GEMMI C FUNCTION stbsp__raise_to_power10
// Gemmi✔️✔️: static void stbsp__raise_to_power10(double *ohi, double *olo, double d, stbsp__int32 power) // power can be -323 to +350
// Gemmi✔️✔️: {
// Gemmi✔️✔️:    double ph, pl;
// Gemmi✔️✔️:    if ((power >= 0) && (power <= 22)) {
// Gemmi✔️✔️:       stbsp__ddmulthi(ph, pl, d, stbsp__bot[power]);
// Gemmi✔️✔️:    } else {
// Gemmi✔️✔️:       stbsp__int32 e, et, eb;
// Gemmi✔️✔️:       double p2h, p2l;
// Gemmi✔️✔️:
// Gemmi✔️✔️:       e = power;
// Gemmi✔️✔️:       if (power < 0)
// Gemmi✔️✔️:          e = -e;
// Gemmi✔️✔️:       et = (e * 0x2c9) >> 14; /* %23 */
// Gemmi✔️✔️:       if (et > 13)
// Gemmi✔️✔️:          et = 13;
// Gemmi✔️✔️:       eb = e - (et * 23);
// Gemmi✔️✔️:
// Gemmi✔️✔️:       ph = d;
// Gemmi✔️✔️:       pl = 0.0;
// Gemmi✔️✔️:       if (power < 0) {
// Gemmi✔️✔️:          if (eb) {
// Gemmi✔️✔️:             --eb;
// Gemmi✔️✔️:             stbsp__ddmulthi(ph, pl, d, stbsp__negbot[eb]);
// Gemmi✔️✔️:             stbsp__ddmultlos(ph, pl, d, stbsp__negboterr[eb]);
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          if (et) {
// Gemmi✔️✔️:             stbsp__ddrenorm(ph, pl);
// Gemmi✔️✔️:             --et;
// Gemmi✔️✔️:             stbsp__ddmulthi(p2h, p2l, ph, stbsp__negtop[et]);
// Gemmi✔️✔️:             stbsp__ddmultlo(p2h, p2l, ph, pl, stbsp__negtop[et], stbsp__negtoperr[et]);
// Gemmi✔️✔️:             ph = p2h;
// Gemmi✔️✔️:             pl = p2l;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:       } else {
// Gemmi✔️✔️:          if (eb) {
// Gemmi✔️✔️:             e = eb;
// Gemmi✔️✔️:             if (eb > 22)
// Gemmi✔️✔️:                eb = 22;
// Gemmi✔️✔️:             e -= eb;
// Gemmi✔️✔️:             stbsp__ddmulthi(ph, pl, d, stbsp__bot[eb]);
// Gemmi✔️✔️:             if (e) {
// Gemmi✔️✔️:                stbsp__ddrenorm(ph, pl);
// Gemmi✔️✔️:                stbsp__ddmulthi(p2h, p2l, ph, stbsp__bot[e]);
// Gemmi✔️✔️:                stbsp__ddmultlos(p2h, p2l, stbsp__bot[e], pl);
// Gemmi✔️✔️:                ph = p2h;
// Gemmi✔️✔️:                pl = p2l;
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          if (et) {
// Gemmi✔️✔️:             stbsp__ddrenorm(ph, pl);
// Gemmi✔️✔️:             --et;
// Gemmi✔️✔️:             stbsp__ddmulthi(p2h, p2l, ph, stbsp__top[et]);
// Gemmi✔️✔️:             stbsp__ddmultlo(p2h, p2l, ph, pl, stbsp__top[et], stbsp__toperr[et]);
// Gemmi✔️✔️:             ph = p2h;
// Gemmi✔️✔️:             pl = p2l;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:    stbsp__ddrenorm(ph, pl);
// Gemmi✔️✔️:    *ohi = ph;
// Gemmi✔️✔️:    *olo = pl;
// Gemmi✔️✔️: }
// END GEMMI C FUNCTION
fn raise_to_power10(value: f64, power: i32) -> (f64, f64) {
    let (mut high, mut low) = if (0..=22).contains(&power) {
        double_double_multiply_high(value, BOT[power as usize])
    } else {
        let exponent = power.unsigned_abs() as i32;
        let mut top_index = (exponent * 0x2c9) >> 14;
        top_index = top_index.min(13);
        let mut remainder = exponent - top_index * 23;
        let mut high = value;
        let mut low = 0.0;
        if power < 0 {
            if remainder != 0 {
                remainder -= 1;
                (high, low) = double_double_multiply_high(value, NEG_BOT[remainder as usize]);
                low += value * NEG_BOT_ERR[remainder as usize];
            }
            if top_index != 0 {
                (high, low) = renormalize(high, low);
                top_index -= 1;
                let (next_high, mut next_low) =
                    double_double_multiply_high(high, NEG_TOP[top_index as usize]);
                next_low +=
                    high * NEG_TOP_ERR[top_index as usize] + low * NEG_TOP[top_index as usize];
                high = next_high;
                low = next_low;
            }
        } else {
            if remainder != 0 {
                let remaining = remainder;
                remainder = remainder.min(22);
                (high, low) = double_double_multiply_high(value, BOT[remainder as usize]);
                let second = remaining - remainder;
                if second != 0 {
                    (high, low) = renormalize(high, low);
                    let (next_high, mut next_low) =
                        double_double_multiply_high(high, BOT[second as usize]);
                    next_low += BOT[second as usize] * low;
                    high = next_high;
                    low = next_low;
                }
            }
            if top_index != 0 {
                (high, low) = renormalize(high, low);
                top_index -= 1;
                let (next_high, mut next_low) =
                    double_double_multiply_high(high, TOP[top_index as usize]);
                next_low += high * TOP_ERR[top_index as usize] + low * TOP[top_index as usize];
                high = next_high;
                low = next_low;
            }
        }
        (high, low)
    };
    (high, low) = renormalize(high, low);
    (high, low)
}

// BEGIN GEMMI C FUNCTION stbsp__ddtoS64
// Gemmi✔️✔️: #define stbsp__ddtoS64(ob, xh, xl)          \
// Gemmi✔️✔️:    {                                        \
// Gemmi✔️✔️:       double ahi = 0, alo, vh, t;           \
// Gemmi✔️✔️:       ob = (stbsp__int64)xh;                \
// Gemmi✔️✔️:       vh = (double)ob;                      \
// Gemmi✔️✔️:       ahi = (xh - vh);                      \
// Gemmi✔️✔️:       t = (ahi - xh);                       \
// Gemmi✔️✔️:       alo = (xh - (ahi - t)) - (vh + t);    \
// Gemmi✔️✔️:       ob += (stbsp__int64)(ahi + alo + xl); \
// Gemmi✔️✔️:    }
// END GEMMI C FUNCTION
fn double_double_to_i64(high: f64, low: f64) -> i64 {
    let mut output = high as i64;
    let integer_high = output as f64;
    let residual_high = high - integer_high;
    let transfer = residual_high - high;
    let residual_low = (high - (residual_high - transfer)) - (integer_high + transfer);
    output = output.wrapping_add((residual_high + residual_low + low) as i64);
    output
}

// BEGIN GEMMI C FUNCTION stbsp__real_to_str
// Gemmi✔️✔️: static stbsp__int32 stbsp__real_to_str(char const **start, stbsp__uint32 *len, char *out, stbsp__int32 *decimal_pos, double value, stbsp__uint32 frac_digits)
// Gemmi✔️✔️: {
// Gemmi✔️✔️:    double d;
// Gemmi✔️✔️:    stbsp__int64 bits = 0;
// Gemmi✔️✔️:    stbsp__int32 expo, e, ng, tens;
// Gemmi✔️✔️:
// Gemmi✔️✔️:    d = value;
// Gemmi✔️✔️:    STBSP__COPYFP(bits, d);
// Gemmi✔️✔️:    expo = (stbsp__int32)((bits >> 52) & 2047);
// Gemmi✔️✔️:    ng = (stbsp__int32)((stbsp__uint64) bits >> 63);
// Gemmi✔️✔️:    if (ng)
// Gemmi✔️✔️:       d = -d;
// Gemmi✔️✔️:
// Gemmi✔️✔️:    if (expo == 2047) // is nan or inf?
// Gemmi✔️✔️:    {
// Gemmi✔️✔️:       *start = (bits & ((((stbsp__uint64)1) << 52) - 1)) ? "NaN" : "Inf";
// Gemmi✔️✔️:       *decimal_pos = STBSP__SPECIAL;
// Gemmi✔️✔️:       *len = 3;
// Gemmi✔️✔️:       return ng;
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:
// Gemmi✔️✔️:    if (expo == 0) // is zero or denormal
// Gemmi✔️✔️:    {
// Gemmi✔️✔️:       if (((stbsp__uint64) bits << 1) == 0) // do zero
// Gemmi✔️✔️:       {
// Gemmi✔️✔️:          *decimal_pos = 1;
// Gemmi✔️✔️:          *start = out;
// Gemmi✔️✔️:          out[0] = '0';
// Gemmi✔️✔️:          *len = 1;
// Gemmi✔️✔️:          return ng;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       // find the right expo for denormals
// Gemmi✔️✔️:       {
// Gemmi✔️✔️:          stbsp__int64 v = ((stbsp__uint64)1) << 51;
// Gemmi✔️✔️:          while ((bits & v) == 0) {
// Gemmi✔️✔️:             --expo;
// Gemmi✔️✔️:             v >>= 1;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:
// Gemmi✔️✔️:    // find the decimal exponent as well as the decimal bits of the value
// Gemmi✔️✔️:    {
// Gemmi✔️✔️:       double ph, pl;
// Gemmi✔️✔️:
// Gemmi✔️✔️:       // log10 estimate - very specifically tweaked to hit or undershoot by no more than 1 of log10 of all expos 1..2046
// Gemmi✔️✔️:       tens = expo - 1023;
// Gemmi✔️✔️:       tens = (tens < 0) ? ((tens * 617) / 2048) : (((tens * 1233) / 4096) + 1);
// Gemmi✔️✔️:
// Gemmi✔️✔️:       // move the significant bits into position and stick them into an int
// Gemmi✔️✔️:       stbsp__raise_to_power10(&ph, &pl, d, 18 - tens);
// Gemmi✔️✔️:
// Gemmi✔️✔️:       // get full as much precision from double-double as possible
// Gemmi✔️✔️:       stbsp__ddtoS64(bits, ph, pl);
// Gemmi✔️✔️:
// Gemmi✔️✔️:       // check if we undershot
// Gemmi✔️✔️:       if (((stbsp__uint64)bits) >= stbsp__tento19th)
// Gemmi✔️✔️:          ++tens;
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:
// Gemmi✔️✔️:    // now do the rounding in integer land
// Gemmi✔️✔️:    frac_digits = (frac_digits & 0x80000000) ? ((frac_digits & 0x7ffffff) + 1) : (tens + frac_digits);
// Gemmi✔️✔️:    if ((frac_digits < 24)) {
// Gemmi✔️✔️:       stbsp__uint32 dg = 1;
// Gemmi✔️✔️:       if ((stbsp__uint64)bits >= stbsp__powten[9])
// Gemmi✔️✔️:          dg = 10;
// Gemmi✔️✔️:       while ((stbsp__uint64)bits >= stbsp__powten[dg]) {
// Gemmi✔️✔️:          ++dg;
// Gemmi✔️✔️:          if (dg == 20)
// Gemmi✔️✔️:             goto noround;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       if (frac_digits < dg) {
// Gemmi✔️✔️:          stbsp__uint64 r;
// Gemmi✔️✔️:          // add 0.5 at the right position and round
// Gemmi✔️✔️:          e = dg - frac_digits;
// Gemmi✔️✔️:          if ((stbsp__uint32)e >= 24)
// Gemmi✔️✔️:             goto noround;
// Gemmi✔️✔️:          r = stbsp__powten[e];
// Gemmi✔️✔️:          bits = bits + (r / 2);
// Gemmi✔️✔️:          if ((stbsp__uint64)bits >= stbsp__powten[dg])
// Gemmi✔️✔️:             ++tens;
// Gemmi✔️✔️:          bits /= r;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:    noround:;
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:
// Gemmi✔️✔️:    // kill long trailing runs of zeros
// Gemmi✔️✔️:    if (bits) {
// Gemmi✔️✔️:       stbsp__uint32 n;
// Gemmi✔️✔️:       for (;;) {
// Gemmi✔️✔️:          if (bits <= 0xffffffff)
// Gemmi✔️✔️:             break;
// Gemmi✔️✔️:          if (bits % 1000)
// Gemmi✔️✔️:             goto donez;
// Gemmi✔️✔️:          bits /= 1000;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       n = (stbsp__uint32)bits;
// Gemmi✔️✔️:       while ((n % 1000) == 0)
// Gemmi✔️✔️:          n /= 1000;
// Gemmi✔️✔️:       bits = n;
// Gemmi✔️✔️:    donez:;
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:
// Gemmi✔️✔️:    // convert to string
// Gemmi✔️✔️:    out += 64;
// Gemmi✔️✔️:    e = 0;
// Gemmi✔️✔️:    for (;;) {
// Gemmi✔️✔️:       stbsp__uint32 n;
// Gemmi✔️✔️:       char *o = out - 8;
// Gemmi✔️✔️:       // do the conversion in chunks of U32s (avoid most 64-bit divides, worth it, constant denomiators be damned)
// Gemmi✔️✔️:       if (bits >= 100000000) {
// Gemmi✔️✔️:          n = (stbsp__uint32)(bits % 100000000);
// Gemmi✔️✔️:          bits /= 100000000;
// Gemmi✔️✔️:       } else {
// Gemmi✔️✔️:          n = (stbsp__uint32)bits;
// Gemmi✔️✔️:          bits = 0;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       while (n) {
// Gemmi✔️✔️:          out -= 2;
// Gemmi✔️✔️:          *(stbsp__uint16 *)out = *(stbsp__uint16 *)&stbsp__digitpair.pair[(n % 100) * 2];
// Gemmi✔️✔️:          n /= 100;
// Gemmi✔️✔️:          e += 2;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       if (bits == 0) {
// Gemmi✔️✔️:          if ((e) && (out[0] == '0')) {
// Gemmi✔️✔️:             ++out;
// Gemmi✔️✔️:             --e;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          break;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       while (out != o) {
// Gemmi✔️✔️:          *--out = '0';
// Gemmi✔️✔️:          ++e;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:    }
// Gemmi✔️✔️:
// Gemmi✔️✔️:    *decimal_pos = tens;
// Gemmi✔️✔️:    *start = out;
// Gemmi✔️✔️:    *len = e;
// Gemmi✔️✔️:    return ng;
// Gemmi✔️✔️: }
// END GEMMI C FUNCTION
fn real_to_digits(value: f64, requested_digits: u32) -> (bool, RealDigits) {
    let raw = value.to_bits();
    let mut exponent = ((raw >> 52) & 2047) as i32;
    let negative = raw >> 63 != 0;
    let magnitude = if negative { -value } else { value };

    if exponent == 2047 {
        let special = if raw & ((1_u64 << 52) - 1) != 0 {
            "NaN"
        } else {
            "Inf"
        };
        return (negative, RealDigits::Special(special));
    }
    if exponent == 0 {
        if raw << 1 == 0 {
            return (
                negative,
                RealDigits::Finite {
                    digits: "0".to_string(),
                    decimal_pos: 1,
                },
            );
        }
        let mut bit = 1_u64 << 51;
        while raw & bit == 0 {
            exponent -= 1;
            bit >>= 1;
        }
    }

    let mut decimal_pos = exponent - 1023;
    decimal_pos = if decimal_pos < 0 {
        decimal_pos * 617 / 2048
    } else {
        decimal_pos * 1233 / 4096 + 1
    };
    let (high, low) = raise_to_power10(magnitude, 18 - decimal_pos);
    let mut bits = double_double_to_i64(high, low) as u64;
    if bits >= 1_000_000_000_000_000_000 {
        decimal_pos += 1;
    }

    let fraction_digits = if requested_digits & 0x8000_0000 != 0 {
        (requested_digits & 0x07ff_ffff) + 1
    } else {
        (decimal_pos as u32).wrapping_add(requested_digits)
    };
    if fraction_digits < 24 {
        let mut digit_count = if bits >= POW_TEN[9] { 10 } else { 1 };
        while digit_count < 20 && bits >= POW_TEN[digit_count] {
            digit_count += 1;
        }
        if digit_count < 20 && fraction_digits < digit_count as u32 {
            let discarded = digit_count as u32 - fraction_digits;
            if discarded < 24 {
                let power = POW_TEN[discarded as usize];
                bits = bits.wrapping_add(power / 2);
                if bits >= POW_TEN[digit_count] {
                    decimal_pos += 1;
                }
                bits /= power;
            }
        }
    }

    if bits != 0 {
        while bits > u32::MAX as u64 && bits % 1_000 == 0 {
            bits /= 1_000;
        }
        if bits <= u32::MAX as u64 {
            while bits % 1_000 == 0 {
                bits /= 1_000;
            }
        }
    }
    (
        negative,
        RealDigits::Finite {
            digits: bits.to_string(),
            decimal_pos,
        },
    )
}

// BEGIN GEMMI C FUNCTION stbsp_vsprintfcb case 'g'
// Gemmi✔️✔️:       case 'G': // float
// Gemmi✔️✔️:       case 'g': // float
// Gemmi✔️✔️:          h = (f[0] == 'G') ? hexu : hex;
// Gemmi✔️✔️:          fv = va_arg(va, double);
// Gemmi✔️✔️:          if (pr == -1)
// Gemmi✔️✔️:             pr = 6;
// Gemmi✔️✔️:          else if (pr == 0)
// Gemmi✔️✔️:             pr = 1; // default is 6
// Gemmi✔️✔️:          // read the double into a string
// Gemmi✔️✔️:          if (stbsp__real_to_str(&sn, &l, num, &dp, fv, (pr - 1) | 0x80000000))
// Gemmi✔️✔️:             fl |= STBSP__NEGATIVE;
// Gemmi✔️✔️:
// Gemmi✔️✔️:          // clamp the precision and delete extra zeros after clamp
// Gemmi✔️✔️:          n = pr;
// Gemmi✔️✔️:          if (l > (stbsp__uint32)pr)
// Gemmi✔️✔️:             l = pr;
// Gemmi✔️✔️:          while ((l > 1) && (pr) && (sn[l - 1] == '0')) {
// Gemmi✔️✔️:             --pr;
// Gemmi✔️✔️:             --l;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:
// Gemmi✔️✔️:          // should we use %e
// Gemmi✔️✔️:          if ((dp <= -4) || (dp > (stbsp__int32)n)) {
// Gemmi✔️✔️:             if (pr > (stbsp__int32)l)
// Gemmi✔️✔️:                pr = l - 1;
// Gemmi✔️✔️:             else if (pr)
// Gemmi✔️✔️:                --pr; // when using %e, there is one digit before the decimal
// Gemmi✔️✔️:             goto doexpfromg;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          // this is the insane action to get the pr to match %g semantics for %f
// Gemmi✔️✔️:          if (dp > 0) {
// Gemmi✔️✔️:             pr = (dp < (stbsp__int32)l) ? l - dp : 0;
// Gemmi✔️✔️:          } else {
// Gemmi✔️✔️:             pr = -dp + ((pr > (stbsp__int32)l) ? (stbsp__int32) l : pr);
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          goto dofloatfromg;
// END GEMMI C FUNCTION
pub(super) fn format_general(value: f64, precision: usize) -> String {
    let precision = precision.max(1);
    let (negative, real) = real_to_digits(value, (precision as u32 - 1) | 0x8000_0000);
    let sign = if negative { "-" } else { "" };
    let RealDigits::Finite {
        mut digits,
        decimal_pos,
    } = real
    else {
        let RealDigits::Special(special) = real else {
            unreachable!()
        };
        return format!("{sign}{special}");
    };

    digits.truncate(precision);
    while digits.len() > 1 && digits.ends_with('0') {
        digits.pop();
    }
    let digit_count = digits.len();
    if decimal_pos <= -4 || decimal_pos > precision as i32 {
        let mut output = String::with_capacity(digit_count + 7);
        output.push_str(sign);
        output.push(digits.as_bytes()[0] as char);
        if digit_count > 1 {
            output.push('.');
            output.push_str(&digits[1..]);
        }
        let exponent = decimal_pos - 1;
        output.push('e');
        output.push(if exponent < 0 { '-' } else { '+' });
        let absolute = exponent.unsigned_abs();
        if absolute < 10 {
            output.push('0');
        }
        output.push_str(&absolute.to_string());
        return output;
    }

    let mut output = String::with_capacity(sign.len() + precision + 8);
    output.push_str(sign);
    if decimal_pos <= 0 {
        output.push('0');
        output.push('.');
        output.extend(std::iter::repeat_n('0', (-decimal_pos) as usize));
        output.push_str(&digits);
    } else if decimal_pos as usize >= digit_count {
        output.push_str(&digits);
        output.extend(std::iter::repeat_n('0', decimal_pos as usize - digit_count));
    } else {
        let split = decimal_pos as usize;
        output.push_str(&digits[..split]);
        output.push('.');
        output.push_str(&digits[split..]);
    }
    output
}

// BEGIN GEMMI C FUNCTION stbsp_vsprintfcb case 'f'
// Gemmi✔️✔️:       case 'f': // float
// Gemmi✔️✔️:          fv = va_arg(va, double);
// Gemmi✔️✔️:       doafloat:
// Gemmi✔️✔️:          // do kilos
// Gemmi✔️✔️:          if (fl & STBSP__METRIC_SUFFIX) {
// Gemmi✔️✔️:             double divisor;
// Gemmi✔️✔️:             divisor = 1000.0f;
// Gemmi✔️✔️:             if (fl & STBSP__METRIC_1024)
// Gemmi✔️✔️:                divisor = 1024.0;
// Gemmi✔️✔️:             while (fl < 0x4000000) {
// Gemmi✔️✔️:                if ((fv < divisor) && (fv > -divisor))
// Gemmi✔️✔️:                   break;
// Gemmi✔️✔️:                fv /= divisor;
// Gemmi✔️✔️:                fl += 0x1000000;
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          if (pr == -1)
// Gemmi✔️✔️:             pr = 6; // default is 6
// Gemmi✔️✔️:          // read the double into a string
// Gemmi✔️✔️:          if (stbsp__real_to_str(&sn, &l, num, &dp, fv, pr))
// Gemmi✔️✔️:             fl |= STBSP__NEGATIVE;
// Gemmi✔️✔️:       dofloatfromg:
// Gemmi✔️✔️:          tail[0] = 0;
// Gemmi✔️✔️:          stbsp__lead_sign(fl, lead);
// Gemmi✔️✔️:          if (dp == STBSP__SPECIAL) {
// Gemmi✔️✔️:             s = (char *)sn;
// Gemmi✔️✔️:             cs = 0;
// Gemmi✔️✔️:             pr = 0;
// Gemmi✔️✔️:             goto scopy;
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          s = num + 64;
// Gemmi✔️✔️:
// Gemmi✔️✔️:          // handle the three decimal varieties
// Gemmi✔️✔️:          if (dp <= 0) {
// Gemmi✔️✔️:             stbsp__int32 i;
// Gemmi✔️✔️:             // handle 0.000*000xxxx
// Gemmi✔️✔️:             *s++ = '0';
// Gemmi✔️✔️:             if (pr)
// Gemmi✔️✔️:                *s++ = stbsp__period;
// Gemmi✔️✔️:             n = -dp;
// Gemmi✔️✔️:             if ((stbsp__int32)n > pr)
// Gemmi✔️✔️:                n = pr;
// Gemmi✔️✔️:             i = n;
// Gemmi✔️✔️:             while (i) {
// Gemmi✔️✔️:                if ((((stbsp__uintptr)s) & 3) == 0)
// Gemmi✔️✔️:                   break;
// Gemmi✔️✔️:                *s++ = '0';
// Gemmi✔️✔️:                --i;
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:             while (i >= 4) {
// Gemmi✔️✔️:                *(stbsp__uint32 *)s = 0x30303030;
// Gemmi✔️✔️:                s += 4;
// Gemmi✔️✔️:                i -= 4;
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:             while (i) {
// Gemmi✔️✔️:                *s++ = '0';
// Gemmi✔️✔️:                --i;
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:             if ((stbsp__int32)(l + n) > pr)
// Gemmi✔️✔️:                l = pr - n;
// Gemmi✔️✔️:             i = l;
// Gemmi✔️✔️:             while (i) {
// Gemmi✔️✔️:                *s++ = *sn++;
// Gemmi✔️✔️:                --i;
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:             tz = pr - (n + l);
// Gemmi✔️✔️:             cs = 1 + (3 << 24); // how many tens did we write (for commas below)
// Gemmi✔️✔️:          } else {
// Gemmi✔️✔️:             cs = (fl & STBSP__TRIPLET_COMMA) ? ((600 - (stbsp__uint32)dp) % 3) : 0;
// Gemmi✔️✔️:             if ((stbsp__uint32)dp >= l) {
// Gemmi✔️✔️:                // handle xxxx000*000.0
// Gemmi✔️✔️:                n = 0;
// Gemmi✔️✔️:                for (;;) {
// Gemmi✔️✔️:                   if ((fl & STBSP__TRIPLET_COMMA) && (++cs == 4)) {
// Gemmi✔️✔️:                      cs = 0;
// Gemmi✔️✔️:                      *s++ = stbsp__comma;
// Gemmi✔️✔️:                   } else {
// Gemmi✔️✔️:                      *s++ = sn[n];
// Gemmi✔️✔️:                      ++n;
// Gemmi✔️✔️:                      if (n >= l)
// Gemmi✔️✔️:                         break;
// Gemmi✔️✔️:                   }
// Gemmi✔️✔️:                }
// Gemmi✔️✔️:                if (n < (stbsp__uint32)dp) {
// Gemmi✔️✔️:                   n = dp - n;
// Gemmi✔️✔️:                   if ((fl & STBSP__TRIPLET_COMMA) == 0) {
// Gemmi✔️✔️:                      while (n) {
// Gemmi✔️✔️:                         if ((((stbsp__uintptr)s) & 3) == 0)
// Gemmi✔️✔️:                            break;
// Gemmi✔️✔️:                         *s++ = '0';
// Gemmi✔️✔️:                         --n;
// Gemmi✔️✔️:                      }
// Gemmi✔️✔️:                      while (n >= 4) {
// Gemmi✔️✔️:                         *(stbsp__uint32 *)s = 0x30303030;
// Gemmi✔️✔️:                         s += 4;
// Gemmi✔️✔️:                         n -= 4;
// Gemmi✔️✔️:                      }
// Gemmi✔️✔️:                   }
// Gemmi✔️✔️:                   while (n) {
// Gemmi✔️✔️:                      if ((fl & STBSP__TRIPLET_COMMA) && (++cs == 4)) {
// Gemmi✔️✔️:                         cs = 0;
// Gemmi✔️✔️:                         *s++ = stbsp__comma;
// Gemmi✔️✔️:                      } else {
// Gemmi✔️✔️:                         *s++ = '0';
// Gemmi✔️✔️:                         --n;
// Gemmi✔️✔️:                      }
// Gemmi✔️✔️:                   }
// Gemmi✔️✔️:                }
// Gemmi✔️✔️:                cs = (int)(s - (num + 64)) + (3 << 24); // cs is how many tens
// Gemmi✔️✔️:                if (pr) {
// Gemmi✔️✔️:                   *s++ = stbsp__period;
// Gemmi✔️✔️:                   tz = pr;
// Gemmi✔️✔️:                }
// Gemmi✔️✔️:             } else {
// Gemmi✔️✔️:                // handle xxxxx.xxxx000*000
// Gemmi✔️✔️:                n = 0;
// Gemmi✔️✔️:                for (;;) {
// Gemmi✔️✔️:                   if ((fl & STBSP__TRIPLET_COMMA) && (++cs == 4)) {
// Gemmi✔️✔️:                      cs = 0;
// Gemmi✔️✔️:                      *s++ = stbsp__comma;
// Gemmi✔️✔️:                   } else {
// Gemmi✔️✔️:                      *s++ = sn[n];
// Gemmi✔️✔️:                      ++n;
// Gemmi✔️✔️:                      if (n >= (stbsp__uint32)dp)
// Gemmi✔️✔️:                         break;
// Gemmi✔️✔️:                   }
// Gemmi✔️✔️:                }
// Gemmi✔️✔️:                cs = (int)(s - (num + 64)) + (3 << 24); // cs is how many tens
// Gemmi✔️✔️:                if (pr)
// Gemmi✔️✔️:                   *s++ = stbsp__period;
// Gemmi✔️✔️:                if ((l - dp) > (stbsp__uint32)pr)
// Gemmi✔️✔️:                   l = pr + dp;
// Gemmi✔️✔️:                while (n < l) {
// Gemmi✔️✔️:                   *s++ = sn[n];
// Gemmi✔️✔️:                   ++n;
// Gemmi✔️✔️:                }
// Gemmi✔️✔️:                tz = pr - (l - dp);
// Gemmi✔️✔️:             }
// Gemmi✔️✔️:          }
// Gemmi✔️✔️:          pr = 0;
// END GEMMI C FUNCTION
pub(super) fn format_fixed(value: f64, precision: usize) -> String {
    let (negative, real) = real_to_digits(value, precision as u32);
    let sign = if negative { "-" } else { "" };
    let RealDigits::Finite {
        digits,
        decimal_pos,
    } = real
    else {
        let RealDigits::Special(special) = real else {
            unreachable!()
        };
        return format!("{sign}{special}");
    };

    let mut output = String::with_capacity(sign.len() + digits.len() + precision + 8);
    output.push_str(sign);
    if decimal_pos <= 0 {
        output.push('0');
        if precision != 0 {
            output.push('.');
        }
        let leading_zeros = (-decimal_pos).min(precision as i32) as usize;
        output.extend(std::iter::repeat_n('0', leading_zeros));
        let copied = digits.len().min(precision.saturating_sub(leading_zeros));
        output.push_str(&digits[..copied]);
        output.extend(std::iter::repeat_n('0', precision - leading_zeros - copied));
    } else if decimal_pos as usize >= digits.len() {
        output.push_str(&digits);
        output.extend(std::iter::repeat_n(
            '0',
            decimal_pos as usize - digits.len(),
        ));
        if precision != 0 {
            output.push('.');
            output.extend(std::iter::repeat_n('0', precision));
        }
    } else {
        let split = decimal_pos as usize;
        output.push_str(&digits[..split]);
        if precision != 0 {
            output.push('.');
            let copied_end = digits.len().min(split + precision);
            output.push_str(&digits[split..copied_end]);
            output.extend(std::iter::repeat_n('0', precision - (copied_end - split)));
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gemmi_general_format_locks_stb_halfway_and_special_values() {
        assert_eq!(
            format_general(f32::from_bits(0xc6cf_fa80) as f64, 6),
            "-26621.3"
        );
        assert_eq!(
            format_general(f32::from_bits(0x4838_f2a0) as f64, 6),
            "189387"
        );
        assert_eq!(format_general(0.0, 9), "0");
        assert_eq!(format_general(-0.0, 9), "-0");
        assert_eq!(format_general(f64::INFINITY, 9), "Inf");
        assert_eq!(format_general(f64::NEG_INFINITY, 9), "-Inf");
        assert_eq!(format_general(f64::NAN, 9), "NaN");
    }

    #[test]
    fn gemmi_fixed_format_locks_stb_precision_and_special_values() {
        assert_eq!(format_fixed(2.045, 4), "2.0450");
        assert_eq!(format_fixed(0.00004, 4), "0.0000");
        assert_eq!(format_fixed(-0.0, 4), "-0.0000");
        assert_eq!(format_fixed(123.456789, 4), "123.4568");
        assert_eq!(format_fixed(99999.99996, 4), "100000.0000");
        assert_eq!(format_fixed(0.00123456, 4), "0.0012");
        assert_eq!(format_fixed(f64::INFINITY, 4), "Inf");
    }
}
