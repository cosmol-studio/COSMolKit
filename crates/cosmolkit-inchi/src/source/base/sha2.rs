#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Sha2Context {
    pub(crate) total: [u64; 2],
    pub(crate) state: [u64; 8],
    pub(crate) buffer: [u8; 64],
}

pub(crate) fn sha2_starts(ctx: &mut Sha2Context) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:57 sha2_starts
    // INCHI✔️✔️: void sha2_starts(sha2_context *ctx)
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     ctx->total[0] = 0;
    // INCHI✔️✔️:     ctx->total[1] = 0;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ctx->state[0] = 0x6A09E667;
    // INCHI✔️✔️:     ctx->state[1] = 0xBB67AE85;
    // INCHI✔️✔️:     ctx->state[2] = 0x3C6EF372;
    // INCHI✔️✔️:     ctx->state[3] = 0xA54FF53A;
    // INCHI✔️✔️:     ctx->state[4] = 0x510E527F;
    // INCHI✔️✔️:     ctx->state[5] = 0x9B05688C;
    // INCHI✔️✔️:     ctx->state[6] = 0x1F83D9AB;
    // INCHI✔️✔️:     ctx->state[7] = 0x5BE0CD19;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: sha2_starts
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_starts
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; unsigned long is 64-bit.
    // INCHI✔️✔️: No conditional or macro-only function behavior is active in this source frame.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_starts

    ctx.total = [0, 0];
    ctx.state = [
        0x6a09_e667,
        0xbb67_ae85,
        0x3c6e_f372,
        0xa54f_f53a,
        0x510e_527f,
        0x9b05_688c,
        0x1f83_d9ab,
        0x5be0_cd19,
    ];
}

pub(crate) fn sha2_process(ctx: &mut Sha2Context, data: &[u8; 64]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:72 sha2_process
    // INCHI✔️✔️: static void sha2_process(sha2_context *ctx, unsigned char data[64])
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     unsigned long temp1, temp2, W[64];
    // INCHI✔️✔️:     unsigned long A, B, C, D, E, F, G, H;
    // INCHI✔️✔️:
    // INCHI✔️✔️:     GET_UINT32_BE(W[0], data, 0);
    // INCHI✔️✔️:     GET_UINT32_BE(W[1], data, 4);
    // INCHI✔️✔️:     GET_UINT32_BE(W[2], data, 8);
    // INCHI✔️✔️:     GET_UINT32_BE(W[3], data, 12);
    // INCHI✔️✔️:     GET_UINT32_BE(W[4], data, 16);
    // INCHI✔️✔️:     GET_UINT32_BE(W[5], data, 20);
    // INCHI✔️✔️:     GET_UINT32_BE(W[6], data, 24);
    // INCHI✔️✔️:     GET_UINT32_BE(W[7], data, 28);
    // INCHI✔️✔️:     GET_UINT32_BE(W[8], data, 32);
    // INCHI✔️✔️:     GET_UINT32_BE(W[9], data, 36);
    // INCHI✔️✔️:     GET_UINT32_BE(W[10], data, 40);
    // INCHI✔️✔️:     GET_UINT32_BE(W[11], data, 44);
    // INCHI✔️✔️:     GET_UINT32_BE(W[12], data, 48);
    // INCHI✔️✔️:     GET_UINT32_BE(W[13], data, 52);
    // INCHI✔️✔️:     GET_UINT32_BE(W[14], data, 56);
    // INCHI✔️✔️:     GET_UINT32_BE(W[15], data, 60);
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define SHR(x, n) ((x & 0xFFFFFFFF) >> n)
    // INCHI✔️✔️: #define ROTR(x, n) (SHR(x, n) | (x << (32 - n)))
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define S0(x) (ROTR(x, 7) ^ ROTR(x, 18) ^ SHR(x, 3))
    // INCHI✔️✔️: #define S1(x) (ROTR(x, 17) ^ ROTR(x, 19) ^ SHR(x, 10))
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define S2(x) (ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22))
    // INCHI✔️✔️: #define S3(x) (ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25))
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define F0(x, y, z) ((x & y) | (z & (x | y)))
    // INCHI✔️✔️: #define F1(x, y, z) (z ^ (x & (y ^ z)))
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define R(t)                             \
    // INCHI✔️✔️:     (                                    \
    // INCHI✔️✔️:         W[t] = S1(W[t - 2]) + W[t - 7] + \
    // INCHI✔️✔️:                S0(W[t - 15]) + W[t - 16])
    // INCHI✔️✔️:
    // INCHI✔️✔️: #define P(a, b, c, d, e, f, g, h, x, K)          \
    // INCHI✔️✔️:     {                                            \
    // INCHI✔️✔️:         temp1 = h + S3(e) + F1(e, f, g) + K + x; \
    // INCHI✔️✔️:         temp2 = S2(a) + F0(a, b, c);             \
    // INCHI✔️✔️:         d += temp1;                              \
    // INCHI✔️✔️:         h = temp1 + temp2;                       \
    // INCHI✔️✔️:     }
    // INCHI✔️✔️:
    // INCHI✔️✔️:     A = ctx->state[0];
    // INCHI✔️✔️:     B = ctx->state[1];
    // INCHI✔️✔️:     C = ctx->state[2];
    // INCHI✔️✔️:     D = ctx->state[3];
    // INCHI✔️✔️:     E = ctx->state[4];
    // INCHI✔️✔️:     F = ctx->state[5];
    // INCHI✔️✔️:     G = ctx->state[6];
    // INCHI✔️✔️:     H = ctx->state[7];
    // INCHI✔️✔️:
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, W[0], 0x428A2F98);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, W[1], 0x71374491);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, W[2], 0xB5C0FBCF);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, W[3], 0xE9B5DBA5);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, W[4], 0x3956C25B);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, W[5], 0x59F111F1);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, W[6], 0x923F82A4);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, W[7], 0xAB1C5ED5);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, W[8], 0xD807AA98);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, W[9], 0x12835B01);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, W[10], 0x243185BE);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, W[11], 0x550C7DC3);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, W[12], 0x72BE5D74);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, W[13], 0x80DEB1FE);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, W[14], 0x9BDC06A7);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, W[15], 0xC19BF174);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, R(16), 0xE49B69C1);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, R(17), 0xEFBE4786);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, R(18), 0x0FC19DC6);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, R(19), 0x240CA1CC);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, R(20), 0x2DE92C6F);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, R(21), 0x4A7484AA);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, R(22), 0x5CB0A9DC);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, R(23), 0x76F988DA);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, R(24), 0x983E5152);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, R(25), 0xA831C66D);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, R(26), 0xB00327C8);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, R(27), 0xBF597FC7);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, R(28), 0xC6E00BF3);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, R(29), 0xD5A79147);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, R(30), 0x06CA6351);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, R(31), 0x14292967);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, R(32), 0x27B70A85);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, R(33), 0x2E1B2138);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, R(34), 0x4D2C6DFC);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, R(35), 0x53380D13);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, R(36), 0x650A7354);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, R(37), 0x766A0ABB);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, R(38), 0x81C2C92E);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, R(39), 0x92722C85);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, R(40), 0xA2BFE8A1);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, R(41), 0xA81A664B);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, R(42), 0xC24B8B70);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, R(43), 0xC76C51A3);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, R(44), 0xD192E819);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, R(45), 0xD6990624);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, R(46), 0xF40E3585);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, R(47), 0x106AA070);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, R(48), 0x19A4C116);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, R(49), 0x1E376C08);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, R(50), 0x2748774C);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, R(51), 0x34B0BCB5);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, R(52), 0x391C0CB3);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, R(53), 0x4ED8AA4A);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, R(54), 0x5B9CCA4F);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, R(55), 0x682E6FF3);
    // INCHI✔️✔️:     P(A, B, C, D, E, F, G, H, R(56), 0x748F82EE);
    // INCHI✔️✔️:     P(H, A, B, C, D, E, F, G, R(57), 0x78A5636F);
    // INCHI✔️✔️:     P(G, H, A, B, C, D, E, F, R(58), 0x84C87814);
    // INCHI✔️✔️:     P(F, G, H, A, B, C, D, E, R(59), 0x8CC70208);
    // INCHI✔️✔️:     P(E, F, G, H, A, B, C, D, R(60), 0x90BEFFFA);
    // INCHI✔️✔️:     P(D, E, F, G, H, A, B, C, R(61), 0xA4506CEB);
    // INCHI✔️✔️:     P(C, D, E, F, G, H, A, B, R(62), 0xBEF9A3F7);
    // INCHI✔️✔️:     P(B, C, D, E, F, G, H, A, R(63), 0xC67178F2);
    // INCHI✔️✔️:
    // INCHI✔️✔️:     ctx->state[0] += A;
    // INCHI✔️✔️:     ctx->state[1] += B;
    // INCHI✔️✔️:     ctx->state[2] += C;
    // INCHI✔️✔️:     ctx->state[3] += D;
    // INCHI✔️✔️:     ctx->state[4] += E;
    // INCHI✔️✔️:     ctx->state[5] += F;
    // INCHI✔️✔️:     ctx->state[6] += G;
    // INCHI✔️✔️:     ctx->state[7] += H;
    // INCHI✔️✔️: }
    // END INCHI C FUNCTION: sha2_process
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_process
    // INCHI✔️✔️: #ifndef GET_UINT32_BE
    // INCHI✔️✔️: #define GET_UINT32_BE(n, b, i)                                                                                                                            \
    // INCHI✔️✔️:     {                                                                                                                                                     \
    // INCHI✔️✔️:         (n) = ((unsigned long)(b)[(i)] << 24) | ((unsigned long)(b)[(i) + 1] << 16) | ((unsigned long)(b)[(i) + 2] << 8) | ((unsigned long)(b)[(i) + 3]); \
    // INCHI✔️✔️:     }
    // INCHI✔️✔️: #endif
    // INCHI✔️✔️: typedef struct
    // INCHI✔️✔️: {
    // INCHI✔️✔️:     unsigned long total[2];
    // INCHI✔️✔️:     unsigned long state[8];
    // INCHI✔️✔️:     unsigned char buffer[64];
    // INCHI✔️✔️: } sha2_context;
    // INCHI✔️✔️: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; GET_UINT32_BE is defined by sha2.c.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_process

    #[inline(always)]
    fn shr(value: u64, amount: u32) -> u64 {
        (value & 0xffff_ffff) >> amount
    }

    #[inline(always)]
    fn rotr(value: u64, amount: u32) -> u64 {
        shr(value, amount) | (value << (32 - amount))
    }

    #[inline(always)]
    fn s0(value: u64) -> u64 {
        rotr(value, 7) ^ rotr(value, 18) ^ shr(value, 3)
    }

    #[inline(always)]
    fn s1(value: u64) -> u64 {
        rotr(value, 17) ^ rotr(value, 19) ^ shr(value, 10)
    }

    #[inline(always)]
    fn s2(value: u64) -> u64 {
        rotr(value, 2) ^ rotr(value, 13) ^ rotr(value, 22)
    }

    #[inline(always)]
    fn s3(value: u64) -> u64 {
        rotr(value, 6) ^ rotr(value, 11) ^ rotr(value, 25)
    }

    #[inline(always)]
    fn f0(x: u64, y: u64, z: u64) -> u64 {
        (x & y) | (z & (x | y))
    }

    #[inline(always)]
    fn f1(x: u64, y: u64, z: u64) -> u64 {
        z ^ (x & (y ^ z))
    }

    let mut words = [0_u64; 64];
    for (word, bytes) in words[..16].iter_mut().zip(data.chunks_exact(4)) {
        *word = (u64::from(bytes[0]) << 24)
            | (u64::from(bytes[1]) << 16)
            | (u64::from(bytes[2]) << 8)
            | u64::from(bytes[3]);
    }

    let mut a = ctx.state[0];
    let mut b = ctx.state[1];
    let mut c = ctx.state[2];
    let mut d = ctx.state[3];
    let mut e = ctx.state[4];
    let mut f = ctx.state[5];
    let mut g = ctx.state[6];
    let mut h = ctx.state[7];

    macro_rules! round {
        ($a:ident, $b:ident, $c:ident, $d:ident, $e:ident, $f:ident, $g:ident, $h:ident, $word:expr, $constant:expr) => {{
            let temp1 = $h
                .wrapping_add(s3($e))
                .wrapping_add(f1($e, $f, $g))
                .wrapping_add($constant)
                .wrapping_add($word);
            let temp2 = s2($a).wrapping_add(f0($a, $b, $c));
            $d = $d.wrapping_add(temp1);
            $h = temp1.wrapping_add(temp2);
        }};
    }

    macro_rules! schedule {
        ($index:expr) => {{
            words[$index] = s1(words[$index - 2])
                .wrapping_add(words[$index - 7])
                .wrapping_add(s0(words[$index - 15]))
                .wrapping_add(words[$index - 16]);
            words[$index]
        }};
    }

    round!(a, b, c, d, e, f, g, h, words[0], 0x428a_2f98);
    round!(h, a, b, c, d, e, f, g, words[1], 0x7137_4491);
    round!(g, h, a, b, c, d, e, f, words[2], 0xb5c0_fbcf);
    round!(f, g, h, a, b, c, d, e, words[3], 0xe9b5_dba5);
    round!(e, f, g, h, a, b, c, d, words[4], 0x3956_c25b);
    round!(d, e, f, g, h, a, b, c, words[5], 0x59f1_11f1);
    round!(c, d, e, f, g, h, a, b, words[6], 0x923f_82a4);
    round!(b, c, d, e, f, g, h, a, words[7], 0xab1c_5ed5);
    round!(a, b, c, d, e, f, g, h, words[8], 0xd807_aa98);
    round!(h, a, b, c, d, e, f, g, words[9], 0x1283_5b01);
    round!(g, h, a, b, c, d, e, f, words[10], 0x2431_85be);
    round!(f, g, h, a, b, c, d, e, words[11], 0x550c_7dc3);
    round!(e, f, g, h, a, b, c, d, words[12], 0x72be_5d74);
    round!(d, e, f, g, h, a, b, c, words[13], 0x80de_b1fe);
    round!(c, d, e, f, g, h, a, b, words[14], 0x9bdc_06a7);
    round!(b, c, d, e, f, g, h, a, words[15], 0xc19b_f174);
    round!(a, b, c, d, e, f, g, h, schedule!(16), 0xe49b_69c1);
    round!(h, a, b, c, d, e, f, g, schedule!(17), 0xefbe_4786);
    round!(g, h, a, b, c, d, e, f, schedule!(18), 0x0fc1_9dc6);
    round!(f, g, h, a, b, c, d, e, schedule!(19), 0x240c_a1cc);
    round!(e, f, g, h, a, b, c, d, schedule!(20), 0x2de9_2c6f);
    round!(d, e, f, g, h, a, b, c, schedule!(21), 0x4a74_84aa);
    round!(c, d, e, f, g, h, a, b, schedule!(22), 0x5cb0_a9dc);
    round!(b, c, d, e, f, g, h, a, schedule!(23), 0x76f9_88da);
    round!(a, b, c, d, e, f, g, h, schedule!(24), 0x983e_5152);
    round!(h, a, b, c, d, e, f, g, schedule!(25), 0xa831_c66d);
    round!(g, h, a, b, c, d, e, f, schedule!(26), 0xb003_27c8);
    round!(f, g, h, a, b, c, d, e, schedule!(27), 0xbf59_7fc7);
    round!(e, f, g, h, a, b, c, d, schedule!(28), 0xc6e0_0bf3);
    round!(d, e, f, g, h, a, b, c, schedule!(29), 0xd5a7_9147);
    round!(c, d, e, f, g, h, a, b, schedule!(30), 0x06ca_6351);
    round!(b, c, d, e, f, g, h, a, schedule!(31), 0x1429_2967);
    round!(a, b, c, d, e, f, g, h, schedule!(32), 0x27b7_0a85);
    round!(h, a, b, c, d, e, f, g, schedule!(33), 0x2e1b_2138);
    round!(g, h, a, b, c, d, e, f, schedule!(34), 0x4d2c_6dfc);
    round!(f, g, h, a, b, c, d, e, schedule!(35), 0x5338_0d13);
    round!(e, f, g, h, a, b, c, d, schedule!(36), 0x650a_7354);
    round!(d, e, f, g, h, a, b, c, schedule!(37), 0x766a_0abb);
    round!(c, d, e, f, g, h, a, b, schedule!(38), 0x81c2_c92e);
    round!(b, c, d, e, f, g, h, a, schedule!(39), 0x9272_2c85);
    round!(a, b, c, d, e, f, g, h, schedule!(40), 0xa2bf_e8a1);
    round!(h, a, b, c, d, e, f, g, schedule!(41), 0xa81a_664b);
    round!(g, h, a, b, c, d, e, f, schedule!(42), 0xc24b_8b70);
    round!(f, g, h, a, b, c, d, e, schedule!(43), 0xc76c_51a3);
    round!(e, f, g, h, a, b, c, d, schedule!(44), 0xd192_e819);
    round!(d, e, f, g, h, a, b, c, schedule!(45), 0xd699_0624);
    round!(c, d, e, f, g, h, a, b, schedule!(46), 0xf40e_3585);
    round!(b, c, d, e, f, g, h, a, schedule!(47), 0x106a_a070);
    round!(a, b, c, d, e, f, g, h, schedule!(48), 0x19a4_c116);
    round!(h, a, b, c, d, e, f, g, schedule!(49), 0x1e37_6c08);
    round!(g, h, a, b, c, d, e, f, schedule!(50), 0x2748_774c);
    round!(f, g, h, a, b, c, d, e, schedule!(51), 0x34b0_bcb5);
    round!(e, f, g, h, a, b, c, d, schedule!(52), 0x391c_0cb3);
    round!(d, e, f, g, h, a, b, c, schedule!(53), 0x4ed8_aa4a);
    round!(c, d, e, f, g, h, a, b, schedule!(54), 0x5b9c_ca4f);
    round!(b, c, d, e, f, g, h, a, schedule!(55), 0x682e_6ff3);
    round!(a, b, c, d, e, f, g, h, schedule!(56), 0x748f_82ee);
    round!(h, a, b, c, d, e, f, g, schedule!(57), 0x78a5_636f);
    round!(g, h, a, b, c, d, e, f, schedule!(58), 0x84c8_7814);
    round!(f, g, h, a, b, c, d, e, schedule!(59), 0x8cc7_0208);
    round!(e, f, g, h, a, b, c, d, schedule!(60), 0x90be_fffa);
    round!(d, e, f, g, h, a, b, c, schedule!(61), 0xa450_6ceb);
    round!(c, d, e, f, g, h, a, b, schedule!(62), 0xbef9_a3f7);
    round!(b, c, d, e, f, g, h, a, schedule!(63), 0xc671_78f2);

    ctx.state[0] = ctx.state[0].wrapping_add(a);
    ctx.state[1] = ctx.state[1].wrapping_add(b);
    ctx.state[2] = ctx.state[2].wrapping_add(c);
    ctx.state[3] = ctx.state[3].wrapping_add(d);
    ctx.state[4] = ctx.state[4].wrapping_add(e);
    ctx.state[5] = ctx.state[5].wrapping_add(f);
    ctx.state[6] = ctx.state[6].wrapping_add(g);
    ctx.state[7] = ctx.state[7].wrapping_add(h);
}

pub(crate) fn sha2_update(ctx: &mut Sha2Context, input: &[u8], ilen: i32) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:206 sha2_update
    // INCHI✔️❌: void sha2_update(sha2_context *ctx, unsigned char *input, int ilen)
    // INCHI✔️❌: {
    // INCHI✔️❌:     int fill;
    // INCHI✔️❌:     unsigned long left;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ilen <= 0)
    // INCHI✔️❌:         return;
    // INCHI✔️❌:
    // INCHI✔️❌:     left = ctx->total[0] & 0x3F;
    // INCHI✔️❌:     fill = 64 - left;
    // INCHI✔️❌:
    // INCHI✔️❌:     ctx->total[0] += ilen;
    // INCHI✔️❌:     ctx->total[0] &= 0xFFFFFFFF;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ctx->total[0] < (unsigned long)ilen)
    // INCHI✔️❌:         ctx->total[1]++;
    // INCHI✔️❌:
    // INCHI✔️❌:     if (left && ilen >= fill)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memcpy((void *)(ctx->buffer + left),
    // INCHI✔️❌:                (void *)input, fill);
    // INCHI✔️❌:         sha2_process(ctx, ctx->buffer);
    // INCHI✔️❌:         input += fill;
    // INCHI✔️❌:         ilen -= fill;
    // INCHI✔️❌:         left = 0;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     while (ilen >= 64)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         sha2_process(ctx, input); /* djb-rwth: ignoring LLVM warning as ilen >= 64 just in test case */
    // INCHI✔️❌:         input += 64;
    // INCHI✔️❌:         ilen -= 64;
    // INCHI✔️❌:     }
    // INCHI✔️❌:
    // INCHI✔️❌:     if (ilen > 0)
    // INCHI✔️❌:     {
    // INCHI✔️❌:         memcpy((void *)(ctx->buffer + left),
    // INCHI✔️❌:                (void *)input, ilen);
    // INCHI✔️❌:     }
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: sha2_update
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_update
    // INCHI✔️❌: #include <string.h>
    // INCHI✔️❌: unsigned long is 64-bit and int is 32-bit under the selected GCC/Linux LP64 target.
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux; libc memcpy is active.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_update

    if ilen <= 0 {
        return;
    }

    let mut remaining = ilen as usize;
    let mut input_offset = 0_usize;
    let mut left = (ctx.total[0] & 0x3f) as usize;
    let fill = 64_usize - left;

    ctx.total[0] = ctx.total[0].wrapping_add(ilen as u64);
    ctx.total[0] &= 0xffff_ffff;

    if ctx.total[0] < ilen as u64 {
        ctx.total[1] = ctx.total[1].wrapping_add(1);
    }

    if left != 0 && remaining >= fill {
        ctx.buffer[left..left + fill].copy_from_slice(&input[..fill]);
        // The source reads ctx->buffer in place. Safe Rust copies this fixed block
        // to avoid aliasing the mutable context, adding known hot-path work.
        let block = ctx.buffer;
        sha2_process(ctx, &block);
        input_offset += fill;
        remaining -= fill;
        left = 0;
    }

    while remaining >= 64 {
        let block: &[u8; 64] = input[input_offset..input_offset + 64]
            .try_into()
            .expect("the source loop requires a complete 64-byte block");
        sha2_process(ctx, block);
        input_offset += 64;
        remaining -= 64;
    }

    if remaining > 0 {
        ctx.buffer[left..left + remaining]
            .copy_from_slice(&input[input_offset..input_offset + remaining]);
    }
}

const SHA2_PADDING: [u8; 64] = [
    0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,
];

pub(crate) fn sha2_finish(ctx: &mut Sha2Context, output: &mut [u8; 32]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:257 sha2_finish
    // INCHI✔️❌: void sha2_finish(sha2_context *ctx, unsigned char output[32])
    // INCHI✔️❌: {
    // INCHI✔️❌:     unsigned long last, padn;
    // INCHI✔️❌:     unsigned long high, low;
    // INCHI✔️❌:     unsigned char msglen[8];
    // INCHI✔️❌:
    // INCHI✔️❌:     high = (ctx->total[0] >> 29) | (ctx->total[1] << 3);
    // INCHI✔️❌:     low = (ctx->total[0] << 3);
    // INCHI✔️❌:
    // INCHI✔️❌:     PUT_UINT32_BE(high, msglen, 0);
    // INCHI✔️❌:     PUT_UINT32_BE(low, msglen, 4);
    // INCHI✔️❌:
    // INCHI✔️❌:     last = ctx->total[0] & 0x3F;
    // INCHI✔️❌:     padn = (last < 56) ? (56 - last) : (120 - last);
    // INCHI✔️❌:
    // INCHI✔️❌:     sha2_update(ctx, (unsigned char *)sha2_padding, padn);
    // INCHI✔️❌:     sha2_update(ctx, msglen, 8);
    // INCHI✔️❌:
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[0], output, 0);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[1], output, 4);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[2], output, 8);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[3], output, 12);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[4], output, 16);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[5], output, 20);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[6], output, 24);
    // INCHI✔️❌:     PUT_UINT32_BE(ctx->state[7], output, 28);
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: sha2_finish
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_finish
    // INCHI✔️❌: #ifndef PUT_UINT32_BE
    // INCHI✔️❌: #define PUT_UINT32_BE(n, b, i)                     \
    // INCHI✔️❌:     {                                              \
    // INCHI✔️❌:         (b)[(i)] = (unsigned char)((n) >> 24);     \
    // INCHI✔️❌:         (b)[(i) + 1] = (unsigned char)((n) >> 16); \
    // INCHI✔️❌:         (b)[(i) + 2] = (unsigned char)((n) >> 8);  \
    // INCHI✔️❌:         (b)[(i) + 3] = (unsigned char)((n));       \
    // INCHI✔️❌:     }
    // INCHI✔️❌: #endif
    // INCHI✔️❌: static const unsigned char sha2_padding[64] =
    // INCHI✔️❌:     {
    // INCHI✔️❌:         0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    // INCHI✔️❌:         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    // INCHI✔️❌:         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    // INCHI✔️❌:         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; unsigned long is 64-bit.
    // INCHI✔️❌: PUT_UINT32_BE and sha2_padding are defined by sha2.c and active.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_finish

    #[inline(always)]
    fn put_uint32_be(value: u64, output: &mut [u8], offset: usize) {
        output[offset] = (value >> 24) as u8;
        output[offset + 1] = (value >> 16) as u8;
        output[offset + 2] = (value >> 8) as u8;
        output[offset + 3] = value as u8;
    }

    let high = (ctx.total[0] >> 29) | ctx.total[1].wrapping_shl(3);
    let low = ctx.total[0].wrapping_shl(3);
    let mut message_length = [0_u8; 8];
    put_uint32_be(high, &mut message_length, 0);
    put_uint32_be(low, &mut message_length, 4);

    let last = ctx.total[0] & 0x3f;
    let padding_length = if last < 56 { 56 - last } else { 120 - last };

    sha2_update(ctx, &SHA2_PADDING, padding_length as i32);
    sha2_update(ctx, &message_length, 8);

    for (index, state) in ctx.state.iter().copied().enumerate() {
        put_uint32_be(state, output, index * 4);
    }
}

pub(crate) fn sha2_csum(input: &[u8], ilen: i32, output: &mut [u8; 32]) {
    // BEGIN INCHI C FUNCTION: third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:312 sha2_csum
    // INCHI✔️❌: void sha2_csum(unsigned char *input, int ilen,
    // INCHI✔️❌:                unsigned char output[32])
    // INCHI✔️❌: {
    // INCHI✔️❌:     sha2_context ctx;
    // INCHI✔️❌:
    // INCHI✔️❌:     sha2_starts(&ctx);
    // INCHI✔️❌:     sha2_update(&ctx, input, ilen);
    // INCHI✔️❌:     sha2_finish(&ctx, output);
    // INCHI✔️❌: }
    // END INCHI C FUNCTION: sha2_csum
    // BEGIN INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_csum
    // INCHI✔️❌: COMPILE_ANSI_ONLY; TARGET_API_LIB; GCC/Linux LP64; int is 32-bit.
    // INCHI✔️❌: No conditional or macro-only function behavior is active in this source frame.
    // END INCHI ACTIVE HEADER/MACRO CONFIGURATION: sha2_csum

    let mut ctx = Sha2Context {
        total: [0; 2],
        state: [0; 8],
        buffer: [0; 64],
    };
    sha2_starts(&mut ctx);
    sha2_update(&mut ctx, input, ilen);
    sha2_finish(&mut ctx, output);
}

#[cfg(test)]
mod tests {
    use super::{Sha2Context, sha2_csum, sha2_finish, sha2_process, sha2_starts, sha2_update};

    #[test]
    fn source_port__sha2__sha2_starts__line_57() {
        let buffer = core::array::from_fn(|index| (index as u8).wrapping_mul(37));
        let mut ctx = Sha2Context {
            total: [u64::MAX, 0x0123_4567_89ab_cdef],
            state: [u64::MAX; 8],
            buffer,
        };

        sha2_starts(&mut ctx);

        assert_eq!(ctx.total, [0, 0]);
        assert_eq!(
            ctx.state,
            [
                0x6a09_e667,
                0xbb67_ae85,
                0x3c6e_f372,
                0xa54f_f53a,
                0x510e_527f,
                0x9b05_688c,
                0x1f83_d9ab,
                0x5be0_cd19,
            ]
        );
        assert_eq!(ctx.buffer, buffer);
    }

    #[test]
    fn source_port__sha2__sha2_csum__line_312() {
        const EMPTY_SHA256: [u8; 32] = [
            0xe3, 0xb0, 0xc4, 0x42, 0x98, 0xfc, 0x1c, 0x14, 0x9a, 0xfb, 0xf4, 0xc8, 0x99, 0x6f,
            0xb9, 0x24, 0x27, 0xae, 0x41, 0xe4, 0x64, 0x9b, 0x93, 0x4c, 0xa4, 0x95, 0x99, 0x1b,
            0x78, 0x52, 0xb8, 0x55,
        ];
        const ABC_SHA256: [u8; 32] = [
            0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea, 0x41, 0x41, 0x40, 0xde, 0x5d, 0xae,
            0x22, 0x23, 0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c, 0xb4, 0x10, 0xff, 0x61,
            0xf2, 0x00, 0x15, 0xad,
        ];

        let mut output = [0xa5_u8; 32];
        sha2_csum(&[], 0, &mut output);
        assert_eq!(output, EMPTY_SHA256);

        output.fill(0x5a);
        sha2_csum(b"ignored", -1, &mut output);
        assert_eq!(output, EMPTY_SHA256);

        let abc = *b"abc";
        let abc_before = abc;
        output.fill(0xa5);
        sha2_csum(&abc, 3, &mut output);
        assert_eq!(output, ABC_SHA256);
        assert_eq!(abc, abc_before);

        let input: Vec<u8> = (0_u8..=u8::MAX).collect();
        for ilen in [1, 55, 56, 63, 64, 65, 127, 128, 129, 255, 256] {
            let mut via_csum = [0xa5_u8; 32];
            sha2_csum(&input, ilen, &mut via_csum);

            let mut ctx = Sha2Context {
                total: [u64::MAX; 2],
                state: [u64::MAX; 8],
                buffer: [0xa5; 64],
            };
            sha2_starts(&mut ctx);
            sha2_update(&mut ctx, &input, ilen);
            let mut via_source_sequence = [0x5a_u8; 32];
            sha2_finish(&mut ctx, &mut via_source_sequence);
            assert_eq!(via_csum, via_source_sequence, "ilen={ilen}");
        }
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__sha2_csum__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--sha2-csum-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "sha2_csum");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_u8_array = |value: &Value, expected_length: usize| {
                let values = value
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: bytes must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| {
                        u8::try_from(value.as_u64().expect("byte must be unsigned"))
                            .unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input = parse_u8_array(&official["input"]["bytes"], 1024);
            let input_before = input.clone();
            let ilen = i32::try_from(
                official["input"]["ilen"]
                    .as_i64()
                    .expect("ilen must be signed"),
            )
            .unwrap_or_else(|_| panic!("{case_id}: ilen exceeds i32"));
            let input_pointer_null = official["input"]["input_pointer_null"]
                .as_bool()
                .expect("input_pointer_null must be boolean");
            let mut digest: [u8; 32] = parse_u8_array(&official["input"]["output"], 32)
                .try_into()
                .unwrap_or_else(|_| panic!("{case_id}: output length changed"));

            sha2_csum(
                if input_pointer_null { &[] } else { &input },
                ilen,
                &mut digest,
            );

            assert_eq!(
                input,
                parse_u8_array(&official["output"]["bytes"], 1024),
                "{case_id}"
            );
            assert_eq!(input, input_before, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            assert_eq!(
                digest.as_slice(),
                parse_u8_array(&official["output"]["digest"], 32),
                "{case_id}"
            );
            record_count += 1;
        }
        assert_eq!(record_count, 26);
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__sha2_starts__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--sha2-starts-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "sha2_starts");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_u64_array = |value: &Value, expected_length: usize| {
                let values = value
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: u64 values must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| value.as_u64().expect("word must be unsigned"))
                    .collect::<Vec<_>>()
            };
            let parse_u8_array = |value: &Value, expected_length: usize| {
                parse_u64_array(value, expected_length)
                    .into_iter()
                    .map(|value| {
                        u8::try_from(value).unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let mut ctx = Sha2Context {
                total: parse_u64_array(&official["input"]["total"], 2)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: total length changed")),
                state: parse_u64_array(&official["input"]["state"], 8)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: state length changed")),
                buffer: parse_u8_array(&official["input"]["buffer"], 64)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: buffer length changed")),
            };
            let buffer_before = ctx.buffer;

            sha2_starts(&mut ctx);

            assert_eq!(
                ctx.total.as_slice(),
                parse_u64_array(&official["output"]["total"], 2),
                "{case_id}"
            );
            assert_eq!(
                ctx.state.as_slice(),
                parse_u64_array(&official["output"]["state"], 8),
                "{case_id}"
            );
            assert_eq!(
                ctx.buffer.as_slice(),
                parse_u8_array(&official["output"]["buffer"], 64),
                "{case_id}"
            );
            assert_eq!(ctx.buffer, buffer_before, "{case_id}");
            assert_eq!(official["output"]["buffer_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 32);
    }

    #[test]
    fn source_port__sha2__sha2_process__line_72() {
        let total = [0x0123_4567_89ab_cdef, 0xfedc_ba98_7654_3210];
        let buffer = [0xa5_u8; 64];
        let mut ctx = Sha2Context {
            total,
            state: [
                0x6a09_e667,
                0xbb67_ae85,
                0x3c6e_f372,
                0xa54f_f53a,
                0x510e_527f,
                0x9b05_688c,
                0x1f83_d9ab,
                0x5be0_cd19,
            ],
            buffer,
        };
        let mut data = [0_u8; 64];
        data[0] = 0x80;
        let data_before = data;

        sha2_process(&mut ctx, &data);

        assert_eq!(
            ctx.state.map(|word| word as u32),
            [
                0xe3b0_c442,
                0x98fc_1c14,
                0x9afb_f4c8,
                0x996f_b924,
                0x27ae_41e4,
                0x649b_934c,
                0xa495_991b,
                0x7852_b855,
            ]
        );
        assert_eq!(ctx.total, total);
        assert_eq!(ctx.buffer, buffer);
        assert_eq!(data, data_before);
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__sha2_process__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--sha2-process-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "sha2_process");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_u64_array = |value: &Value, expected_length: usize| {
                let values = value
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: u64 values must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| value.as_u64().expect("word must be unsigned"))
                    .collect::<Vec<_>>()
            };
            let parse_u8_array = |value: &Value, expected_length: usize| {
                parse_u64_array(value, expected_length)
                    .into_iter()
                    .map(|value| {
                        u8::try_from(value).unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let input_total = parse_u64_array(&official["input"]["total"], 2);
            let input_state = parse_u64_array(&official["input"]["state"], 8);
            let input_buffer = parse_u8_array(&official["input"]["buffer"], 64);
            let input_data = parse_u8_array(&official["input"]["data"], 64);
            let mut ctx = Sha2Context {
                total: input_total
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: total length changed")),
                state: input_state
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: state length changed")),
                buffer: input_buffer
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: buffer length changed")),
            };
            let data: [u8; 64] = input_data
                .try_into()
                .unwrap_or_else(|_| panic!("{case_id}: data length changed"));

            sha2_process(&mut ctx, &data);

            assert_eq!(
                ctx.total.as_slice(),
                parse_u64_array(&official["output"]["total"], 2),
                "{case_id}"
            );
            assert_eq!(
                ctx.state.as_slice(),
                parse_u64_array(&official["output"]["state"], 8),
                "{case_id}"
            );
            assert_eq!(
                ctx.buffer.as_slice(),
                parse_u8_array(&official["output"]["buffer"], 64),
                "{case_id}"
            );
            assert_eq!(
                data.as_slice(),
                parse_u8_array(&official["output"]["data"], 64),
                "{case_id}"
            );
            assert_eq!(official["output"]["total_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["buffer_unchanged"], true, "{case_id}");
            assert_eq!(official["output"]["data_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 256);
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__sha2_update__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--sha2-update-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "sha2_update");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_u64_array = |value: &Value, expected_length: usize| {
                let values = value
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: u64 values must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| value.as_u64().expect("word must be unsigned"))
                    .collect::<Vec<_>>()
            };
            let parse_u8_array = |value: &Value, expected_length: usize| {
                parse_u64_array(value, expected_length)
                    .into_iter()
                    .map(|value| {
                        u8::try_from(value).unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let mut ctx = Sha2Context {
                total: parse_u64_array(&official["input"]["total"], 2)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: total length changed")),
                state: parse_u64_array(&official["input"]["state"], 8)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: state length changed")),
                buffer: parse_u8_array(&official["input"]["buffer"], 64)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: buffer length changed")),
            };
            let input = parse_u8_array(&official["input"]["bytes"], 256);
            let input_before = input.clone();
            let ilen = i32::try_from(
                official["input"]["ilen"]
                    .as_i64()
                    .expect("ilen must be signed"),
            )
            .unwrap_or_else(|_| panic!("{case_id}: ilen exceeds i32"));

            sha2_update(&mut ctx, &input, ilen);

            assert_eq!(
                ctx.total.as_slice(),
                parse_u64_array(&official["output"]["total"], 2),
                "{case_id}"
            );
            assert_eq!(
                ctx.state.as_slice(),
                parse_u64_array(&official["output"]["state"], 8),
                "{case_id}"
            );
            assert_eq!(
                ctx.buffer.as_slice(),
                parse_u8_array(&official["output"]["buffer"], 64),
                "{case_id}"
            );
            assert_eq!(
                input,
                parse_u8_array(&official["output"]["bytes"], 256),
                "{case_id}"
            );
            assert_eq!(input, input_before, "{case_id}");
            assert_eq!(official["output"]["input_unchanged"], true, "{case_id}");
            record_count += 1;
        }
        assert_eq!(record_count, 36);
    }

    #[test]
    #[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
    fn official_c_oracle__sha2_finish__exact() {
        use std::path::Path;
        use std::process::Command;

        use serde_json::Value;

        let repository_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .and_then(Path::parent)
            .expect("cosmolkit-inchi must be located under crates/");
        let runner = repository_root.join("tools/oracles/official_inchi/run.sh");
        let oracle = Command::new("sh")
            .arg(&runner)
            .arg("--sha2-finish-records")
            .current_dir(repository_root)
            .output()
            .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner.display()));
        assert!(
            oracle.status.success(),
            "official C oracle failed with {}:\n{}",
            oracle.status,
            String::from_utf8_lossy(&oracle.stderr)
        );
        let output =
            String::from_utf8(oracle.stdout).expect("official C oracle output must be UTF-8");
        let mut record_count = 0_usize;
        for line in output.lines() {
            let official: Value = serde_json::from_str(line).expect("oracle record must be JSON");
            assert_eq!(official["schema_version"], "cosmolkit-inchi-official-c-v1");
            assert_eq!(official["operation"], "sha2_finish");
            let case_id = official["case_id"].as_str().expect("case_id must be text");
            let parse_u64_array = |value: &Value, expected_length: usize| {
                let values = value
                    .as_array()
                    .unwrap_or_else(|| panic!("{case_id}: u64 values must be an array"));
                assert_eq!(values.len(), expected_length, "{case_id}");
                values
                    .iter()
                    .map(|value| value.as_u64().expect("word must be unsigned"))
                    .collect::<Vec<_>>()
            };
            let parse_u8_array = |value: &Value, expected_length: usize| {
                parse_u64_array(value, expected_length)
                    .into_iter()
                    .map(|value| {
                        u8::try_from(value).unwrap_or_else(|_| panic!("{case_id}: byte exceeds u8"))
                    })
                    .collect::<Vec<_>>()
            };
            let mut ctx = Sha2Context {
                total: parse_u64_array(&official["input"]["total"], 2)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: total length changed")),
                state: parse_u64_array(&official["input"]["state"], 8)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: state length changed")),
                buffer: parse_u8_array(&official["input"]["buffer"], 64)
                    .try_into()
                    .unwrap_or_else(|_| panic!("{case_id}: buffer length changed")),
            };
            let mut digest: [u8; 32] = parse_u8_array(&official["input"]["output"], 32)
                .try_into()
                .unwrap_or_else(|_| panic!("{case_id}: output length changed"));

            sha2_finish(&mut ctx, &mut digest);

            assert_eq!(
                ctx.total.as_slice(),
                parse_u64_array(&official["output"]["total"], 2),
                "{case_id}"
            );
            assert_eq!(
                ctx.state.as_slice(),
                parse_u64_array(&official["output"]["state"], 8),
                "{case_id}"
            );
            assert_eq!(
                ctx.buffer.as_slice(),
                parse_u8_array(&official["output"]["buffer"], 64),
                "{case_id}"
            );
            assert_eq!(
                digest.as_slice(),
                parse_u8_array(&official["output"]["digest"], 32),
                "{case_id}"
            );
            record_count += 1;
        }
        assert_eq!(record_count, 12);
    }

    fn test_context(total: [u64; 2]) -> Sha2Context {
        Sha2Context {
            total,
            state: [
                0x6a09_e667,
                0xbb67_ae85,
                0x3c6e_f372,
                0xa54f_f53a,
                0x510e_527f,
                0x9b05_688c,
                0x1f83_d9ab,
                0x5be0_cd19,
            ],
            buffer: core::array::from_fn(|index| (index as u8).wrapping_mul(13)),
        }
    }

    #[test]
    fn source_port__sha2__sha2_update__line_206() {
        let input: Vec<u8> = (0_u8..=u8::MAX).collect();

        for ilen in [i32::MIN, -1, 0] {
            let mut ctx = test_context([17, 23]);
            let before = ctx.clone();
            sha2_update(&mut ctx, &[], ilen);
            assert_eq!(ctx, before, "ilen={ilen}");
        }

        let mut short = test_context([0, 7]);
        let short_state = short.state;
        let short_buffer = short.buffer;
        sha2_update(&mut short, &input, 17);
        assert_eq!(short.total, [17, 7]);
        assert_eq!(short.state, short_state);
        assert_eq!(&short.buffer[..17], &input[..17]);
        assert_eq!(&short.buffer[17..], &short_buffer[17..]);

        let mut partial_tail = test_context([11, 9]);
        let partial_state = partial_tail.state;
        let partial_buffer = partial_tail.buffer;
        sha2_update(&mut partial_tail, &input, 12);
        assert_eq!(partial_tail.total, [23, 9]);
        assert_eq!(partial_tail.state, partial_state);
        assert_eq!(&partial_tail.buffer[..11], &partial_buffer[..11]);
        assert_eq!(&partial_tail.buffer[11..23], &input[..12]);
        assert_eq!(&partial_tail.buffer[23..], &partial_buffer[23..]);

        let mut fills_partial = test_context([13, 4]);
        let mut expected_fill = fills_partial.clone();
        expected_fill.total = [64, 4];
        expected_fill.buffer[13..].copy_from_slice(&input[..51]);
        let expected_block = expected_fill.buffer;
        sha2_process(&mut expected_fill, &expected_block);
        sha2_update(&mut fills_partial, &input, 51);
        assert_eq!(fills_partial, expected_fill);

        let mut blocks_and_tail = test_context([0xffff_fff0, 0x1234]);
        let mut expected_blocks = blocks_and_tail.clone();
        expected_blocks.total = [0x81, 0x1235];
        expected_blocks.buffer[48..].copy_from_slice(&input[..16]);
        let partial_block = expected_blocks.buffer;
        sha2_process(&mut expected_blocks, &partial_block);
        let first: &[u8; 64] = input[16..80].try_into().expect("64-byte fixture");
        let second: &[u8; 64] = input[80..144].try_into().expect("64-byte fixture");
        sha2_process(&mut expected_blocks, first);
        sha2_process(&mut expected_blocks, second);
        expected_blocks.buffer[0] = input[144];
        sha2_update(&mut blocks_and_tail, &input, 145);
        assert_eq!(blocks_and_tail, expected_blocks);
        assert_eq!(input, (0_u8..=u8::MAX).collect::<Vec<_>>());
    }

    #[test]
    fn source_port__sha2__sha2_finish__line_257() {
        let mut empty = test_context([0, 0]);
        empty.buffer = [0xa5; 64];
        let mut empty_output = [0xa5; 32];
        sha2_finish(&mut empty, &mut empty_output);
        assert_eq!(
            empty_output,
            [
                0xe3, 0xb0, 0xc4, 0x42, 0x98, 0xfc, 0x1c, 0x14, 0x9a, 0xfb, 0xf4, 0xc8, 0x99, 0x6f,
                0xb9, 0x24, 0x27, 0xae, 0x41, 0xe4, 0x64, 0x9b, 0x93, 0x4c, 0xa4, 0x95, 0x99, 0x1b,
                0x78, 0x52, 0xb8, 0x55,
            ]
        );
        assert_eq!(empty.total, [64, 0]);

        let mut abc = test_context([0, 0]);
        abc.buffer = [0; 64];
        sha2_update(&mut abc, b"abc", 3);
        let mut abc_output = [0; 32];
        sha2_finish(&mut abc, &mut abc_output);
        assert_eq!(
            abc_output,
            [
                0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea, 0x41, 0x41, 0x40, 0xde, 0x5d, 0xae,
                0x22, 0x23, 0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c, 0xb4, 0x10, 0xff, 0x61,
                0xf2, 0x00, 0x15, 0xad,
            ]
        );
        assert_eq!(abc.total, [64, 0]);

        for (initial_total, expected_total) in [
            ([55, 7], [64, 7]),
            ([56, 7], [128, 7]),
            ([63, u64::MAX], [128, u64::MAX]),
            ([0xffff_ffff, u64::MAX], [64, 0]),
            (
                [0x1234_5678_0000_003f, 0xfedc_ba98_7654_3210],
                [128, 0xfedc_ba98_7654_3210],
            ),
        ] {
            let mut ctx = test_context(initial_total);
            let mut output = [0xa5; 32];
            sha2_finish(&mut ctx, &mut output);
            assert_eq!(ctx.total, expected_total, "initial_total={initial_total:?}");
            for (word, bytes) in ctx.state.iter().zip(output.chunks_exact(4)) {
                assert_eq!(
                    bytes,
                    [
                        (word >> 24) as u8,
                        (word >> 16) as u8,
                        (word >> 8) as u8,
                        *word as u8,
                    ],
                    "initial_total={initial_total:?}"
                );
            }
            assert_ne!(output, [0xa5; 32], "initial_total={initial_total:?}");
        }
    }
}
