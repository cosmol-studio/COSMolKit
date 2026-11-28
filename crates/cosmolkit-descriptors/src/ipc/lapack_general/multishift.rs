//! Source aggressive early deflation and multishift bulge chase.
use super::super::{
    lapack_blas_level3 as blas3, lapack_symmetric::larfg, lapack_tridiagonal::dlamch,
};
use super::{GeneralError, ix};
use blas3::MatrixSource::{Separate, Within};
use blas3::Transpose::{NoTrans, Trans};

#[inline]
fn at(i: isize, j: isize, ld: usize) -> usize {
    debug_assert!(i >= 1 && j >= 1);
    (j as usize - 1) * ld + i as usize - 1
}

fn laqr1(n: usize, h: &[f64], ldh: usize, sr1: f64, si1: f64, sr2: f64, si2: f64) -> [f64; 3] {
    // LAPACK❗✔️:       SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   SI1, SI2, SR1, SR2
    // LAPACK❗✔️:       INTEGER            LDH, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), V( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  ================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H21S, H31S, S
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.NE.2 .AND. N.NE.3 ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.2 ) THEN
    // LAPACK❗✔️:          S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) )
    // LAPACK❗✔️:          IF( S.EQ.ZERO ) THEN
    // LAPACK❗✔️:             V( 1 ) = ZERO
    // LAPACK❗✔️:             V( 2 ) = ZERO
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             H21S = H( 2, 1 ) / S
    // LAPACK❗✔️:             V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-SR1 )*
    // LAPACK❗✔️:      $               ( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S )
    // LAPACK❗✔️:             V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) ) +
    // LAPACK❗✔️:      $       ABS( H( 3, 1 ) )
    // LAPACK❗✔️:          IF( S.EQ.ZERO ) THEN
    // LAPACK❗✔️:             V( 1 ) = ZERO
    // LAPACK❗✔️:             V( 2 ) = ZERO
    // LAPACK❗✔️:             V( 3 ) = ZERO
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             H21S = H( 2, 1 ) / S
    // LAPACK❗✔️:             H31S = H( 3, 1 ) / S
    // LAPACK❗✔️:             V( 1 ) = ( H( 1, 1 )-SR1 )*( ( H( 1, 1 )-SR2 ) / S ) -
    // LAPACK❗✔️:      $               SI1*( SI2 / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S
    // LAPACK❗✔️:             V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 ) +
    // LAPACK❗✔️:      $               H( 2, 3 )*H31S
    // LAPACK❗✔️:             V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-SR1-SR2 ) +
    // LAPACK❗✔️:      $               H21S*H( 3, 2 )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       END

    let mut v = [0.0; 3];
    if n == 2 {
        let s = (h[0] - sr2).abs() + si2.abs() + h[1].abs();
        if s != 0.0 {
            let h21s = h[1] / s;
            v[0] = h21s * h[ldh] + (h[0] - sr1) * ((h[0] - sr2) / s) - si1 * (si2 / s);
            v[1] = h21s * (h[0] + h[ldh + 1] - sr1 - sr2);
        }
    } else if n == 3 {
        let s = (h[0] - sr2).abs() + si2.abs() + h[1].abs() + h[2].abs();
        if s != 0.0 {
            let h21s = h[1] / s;
            let h31s = h[2] / s;
            v[0] = (h[0] - sr1) * ((h[0] - sr2) / s) - si1 * (si2 / s)
                + h[ldh] * h21s
                + h[2 * ldh] * h31s;
            v[1] = h21s * (h[0] + h[ldh + 1] - sr1 - sr2) + h[2 * ldh + 1] * h31s;
            v[2] = h31s * (h[0] + h[2 * ldh + 2] - sr1 - sr2) + h21s * h[ldh + 2];
        }
    }
    v
}

fn vigilant_deflation(
    h: &mut [f64],
    ldh: usize,
    k: isize,
    ktop: isize,
    kbot: isize,
    smlnum: f64,
    ulp: f64,
) {
    // The identical Ahues--Tisseur check in DLAQR5's two-row and three-row paths.
    macro_rules! h {
        ($i:expr,$j:expr) => {
            h[at($i, $j, ldh)]
        };
    }
    if k < ktop || h!(k + 1, k) == 0.0 {
        return;
    }
    let mut tst1 = h!(k, k).abs() + h!(k + 1, k + 1).abs();
    if tst1 == 0.0 {
        if k >= ktop + 1 {
            tst1 += h!(k, k - 1).abs();
        }
        if k >= ktop + 2 {
            tst1 += h!(k, k - 2).abs();
        }
        if k >= ktop + 3 {
            tst1 += h!(k, k - 3).abs();
        }
        if k <= kbot - 2 {
            tst1 += h!(k + 2, k + 1).abs();
        }
        if k <= kbot - 3 {
            tst1 += h!(k + 3, k + 1).abs();
        }
        if k <= kbot - 4 {
            tst1 += h!(k + 4, k + 1).abs();
        }
    }
    if h!(k + 1, k).abs() <= smlnum.max(ulp * tst1) {
        let h12 = h!(k + 1, k).abs().max(h!(k, k + 1).abs());
        let h21 = h!(k + 1, k).abs().min(h!(k, k + 1).abs());
        let h11 = h!(k + 1, k + 1)
            .abs()
            .max((h!(k, k) - h!(k + 1, k + 1)).abs());
        let h22 = h!(k + 1, k + 1)
            .abs()
            .min((h!(k, k) - h!(k + 1, k + 1)).abs());
        let scl = h11 + h12;
        let tst2 = h22 * (h11 / scl);
        if tst2 == 0.0 || h21 * (h12 / scl) <= smlnum.max(ulp * tst2) {
            h!(k + 1, k) = 0.0;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn laqr5(
    wantt: bool,
    wantz: bool,
    kacc22: usize,
    n: usize,
    ktop: usize,
    kbot: usize,
    nshifts: usize,
    sr: &mut [f64],
    si: &mut [f64],
    h: &mut [f64],
    ldh: usize,
    iloz: usize,
    ihiz: usize,
    z: &mut [f64],
    ldz: usize,
    v: &mut [f64],
    ldv: usize,
    uoff: usize,
    ldu: usize,
    nv: usize,
    wvoff: usize,
    ldwv: usize,
    nh: usize,
    whoff: usize,
    ldwh: usize,
) -> Result<(), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
    // LAPACK❗✔️:      $                   SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,
    // LAPACK❗✔️:      $                   LDU, NV, WV, LDWV, NH, WH, LDWH )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
    // LAPACK❗✔️:      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
    // LAPACK❗✔️:       LOGICAL            WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ),
    // LAPACK❗✔️:      $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ),
    // LAPACK❗✔️:      $                   Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  ================================================================
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM,
    // LAPACK❗✔️:      $                   SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, T1, T2,
    // LAPACK❗✔️:      $                   T3, TST1, TST2, ULP
    // LAPACK❗✔️:       INTEGER            I, I2, I4, INCOL, J, JBOT, JCOL, JLEN,
    // LAPACK❗✔️:      $                   JROW, JTOP, K, K1, KDU, KMS, KRCOL,
    // LAPACK❗✔️:      $                   M, M22, MBOT, MTOP, NBMPS, NDCOL,
    // LAPACK❗✔️:      $                   NS, NU
    // LAPACK❗✔️:       LOGICAL            ACCUM, BMP22
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️:       EXTERNAL           DLAMCH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INTRINSIC          ABS, DBLE, MAX, MIN, MOD
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       DOUBLE PRECISION   VT( 3 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DGEMM, DLABAD, DLACPY, DLAQR1, DLARFG, DLASET,
    // LAPACK❗✔️:      $                   DTRMM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== If there are no shifts, then there is nothing to do. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NSHFTS.LT.2 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== If the active block is empty or 1-by-1, then there
    // LAPACK❗✔️: *     .    is nothing to do. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( KTOP.GE.KBOT )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Shuffle shifts into pairs of real shifts and pairs
    // LAPACK❗✔️: *     .    of complex conjugate shifts assuming complex
    // LAPACK❗✔️: *     .    conjugate shifts are already adjacent to one
    // LAPACK❗✔️: *     .    another. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 I = 1, NSHFTS - 2, 2
    // LAPACK❗✔️:          IF( SI( I ).NE.-SI( I+1 ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             SWAP = SR( I )
    // LAPACK❗✔️:             SR( I ) = SR( I+1 )
    // LAPACK❗✔️:             SR( I+1 ) = SR( I+2 )
    // LAPACK❗✔️:             SR( I+2 ) = SWAP
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             SWAP = SI( I )
    // LAPACK❗✔️:             SI( I ) = SI( I+1 )
    // LAPACK❗✔️:             SI( I+1 ) = SI( I+2 )
    // LAPACK❗✔️:             SI( I+2 ) = SWAP
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== NSHFTS is supposed to be even, but if it is odd,
    // LAPACK❗✔️: *     .    then simply reduce it by one.  The shuffle above
    // LAPACK❗✔️: *     .    ensures that the dropped shift is real and that
    // LAPACK❗✔️: *     .    the remaining shifts are paired. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NS = NSHFTS - MOD( NSHFTS, 2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Machine constants for deflation ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'SAFE MINIMUM' )
    // LAPACK❗✔️:       SAFMAX = ONE / SAFMIN
    // LAPACK❗✔️:       CALL DLABAD( SAFMIN, SAFMAX )
    // LAPACK❗✔️:       ULP = DLAMCH( 'PRECISION' )
    // LAPACK❗✔️:       SMLNUM = SAFMIN*( DBLE( N ) / ULP )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Use accumulated reflections to update far-from-diagonal
    // LAPACK❗✔️: *     .    entries ? ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== clear trash ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( KTOP+2.LE.KBOT )
    // LAPACK❗✔️:      $   H( KTOP+2, KTOP ) = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== NBMPS = number of 2-shift bulges in the chain ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NBMPS = NS / 2
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== KDU = width of slab ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       KDU = 4*NBMPS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Create and chase chains of NBMPS bulges ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 180 INCOL = KTOP - 2*NBMPS + 1, KBOT - 2, 2*NBMPS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        JTOP = Index from which updates from the right start.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( ACCUM ) THEN
    // LAPACK❗✔️:             JTOP = MAX( KTOP, INCOL )
    // LAPACK❗✔️:          ELSE IF( WANTT ) THEN
    // LAPACK❗✔️:             JTOP = 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             JTOP = KTOP
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NDCOL = INCOL + KDU
    // LAPACK❗✔️:          IF( ACCUM )
    // LAPACK❗✔️:      $      CALL DLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Near-the-diagonal bulge chase.  The following loop
    // LAPACK❗✔️: *        .    performs the near-the-diagonal part of a small bulge
    // LAPACK❗✔️: *        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
    // LAPACK❗✔️: *        .    chunk extends from column INCOL to column NDCOL
    // LAPACK❗✔️: *        .    (including both column INCOL and column NDCOL). The
    // LAPACK❗✔️: *        .    following loop chases a 2*NBMPS+1 column long chain of
    // LAPACK❗✔️: *        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
    // LAPACK❗✔️: *        .    may be less than KTOP and and NDCOL may be greater than
    // LAPACK❗✔️: *        .    KBOT indicating phantom columns from which to chase
    // LAPACK❗✔️: *        .    bulges before they are actually introduced or to which
    // LAPACK❗✔️: *        .    to chase bulges beyond column KBOT.)  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 145 KRCOL = INCOL, MIN( INCOL+2*NBMPS-1, KBOT-2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Bulges number MTOP to MBOT are active double implicit
    // LAPACK❗✔️: *           .    shift bulges.  There may or may not also be small
    // LAPACK❗✔️: *           .    2-by-2 bulge, if there is room.  The inactive bulges
    // LAPACK❗✔️: *           .    (if any) must wait until the active bulges have moved
    // LAPACK❗✔️: *           .    down the diagonal to make room.  The phantom matrix
    // LAPACK❗✔️: *           .    paradigm described above helps keep track.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             MTOP = MAX( 1, ( KTOP-KRCOL ) / 2+1 )
    // LAPACK❗✔️:             MBOT = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 2 )
    // LAPACK❗✔️:             M22 = MBOT + 1
    // LAPACK❗✔️:             BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+2*( M22-1 ) ).EQ.
    // LAPACK❗✔️:      $              ( KBOT-2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Generate reflections to chase the chain right
    // LAPACK❗✔️: *           .    one column.  (The minimum value of K is KTOP-1.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF ( BMP22 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Special case: 2-by-2 reflection at bottom treated
    // LAPACK❗✔️: *              .    separately ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                K = KRCOL + 2*( M22-1 )
    // LAPACK❗✔️:                IF( K.EQ.KTOP-1 ) THEN
    // LAPACK❗✔️:                   CALL DLAQR1( 2, H( K+1, K+1 ), LDH, SR( 2*M22-1 ),
    // LAPACK❗✔️:      $                         SI( 2*M22-1 ), SR( 2*M22 ), SI( 2*M22 ),
    // LAPACK❗✔️:      $                         V( 1, M22 ) )
    // LAPACK❗✔️:                   BETA = V( 1, M22 )
    // LAPACK❗✔️:                   CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   BETA = H( K+1, K )
    // LAPACK❗✔️:                   V( 2, M22 ) = H( K+2, K )
    // LAPACK❗✔️:                   CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
    // LAPACK❗✔️:                   H( K+1, K ) = BETA
    // LAPACK❗✔️:                   H( K+2, K ) = ZERO
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Perform update from right within
    // LAPACK❗✔️: *              .    computational window. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                T1 = V( 1, M22 )
    // LAPACK❗✔️:                T2 = T1*V( 2, M22 )
    // LAPACK❗✔️:                DO 30 J = JTOP, MIN( KBOT, K+3 )
    // LAPACK❗✔️:                   REFSUM = H( J, K+1 ) + V( 2, M22 )*H( J, K+2 )
    // LAPACK❗✔️:                   H( J, K+1 ) = H( J, K+1 ) - REFSUM*T1
    // LAPACK❗✔️:                   H( J, K+2 ) = H( J, K+2 ) - REFSUM*T2
    // LAPACK❗✔️:    30          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Perform update from left within
    // LAPACK❗✔️: *              .    computational window. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( ACCUM ) THEN
    // LAPACK❗✔️:                   JBOT = MIN( NDCOL, KBOT )
    // LAPACK❗✔️:                ELSE IF( WANTT ) THEN
    // LAPACK❗✔️:                   JBOT = N
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   JBOT = KBOT
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:                T1 = V( 1, M22 )
    // LAPACK❗✔️:                T2 = T1*V( 2, M22 )
    // LAPACK❗✔️:                DO 40 J = K+1, JBOT
    // LAPACK❗✔️:                   REFSUM = H( K+1, J ) + V( 2, M22 )*H( K+2, J )
    // LAPACK❗✔️:                   H( K+1, J ) = H( K+1, J ) - REFSUM*T1
    // LAPACK❗✔️:                   H( K+2, J ) = H( K+2, J ) - REFSUM*T2
    // LAPACK❗✔️:    40          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== The following convergence test requires that
    // LAPACK❗✔️: *              .    the tradition small-compared-to-nearby-diagonals
    // LAPACK❗✔️: *              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
    // LAPACK❗✔️: *              .    criteria both be satisfied.  The latter improves
    // LAPACK❗✔️: *              .    accuracy in some examples. Falling back on an
    // LAPACK❗✔️: *              .    alternate convergence criterion when TST1 or TST2
    // LAPACK❗✔️: *              .    is zero (as done here) is traditional but probably
    // LAPACK❗✔️: *              .    unnecessary. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( K.GE.KTOP ) THEN
    // LAPACK❗✔️:                   IF( H( K+1, K ).NE.ZERO ) THEN
    // LAPACK❗✔️:                      TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) )
    // LAPACK❗✔️:                      IF( TST1.EQ.ZERO ) THEN
    // LAPACK❗✔️:                         IF( K.GE.KTOP+1 )
    // LAPACK❗✔️:      $                     TST1 = TST1 + ABS( H( K, K-1 ) )
    // LAPACK❗✔️:                         IF( K.GE.KTOP+2 )
    // LAPACK❗✔️:      $                     TST1 = TST1 + ABS( H( K, K-2 ) )
    // LAPACK❗✔️:                         IF( K.GE.KTOP+3 )
    // LAPACK❗✔️:      $                     TST1 = TST1 + ABS( H( K, K-3 ) )
    // LAPACK❗✔️:                         IF( K.LE.KBOT-2 )
    // LAPACK❗✔️:      $                     TST1 = TST1 + ABS( H( K+2, K+1 ) )
    // LAPACK❗✔️:                         IF( K.LE.KBOT-3 )
    // LAPACK❗✔️:      $                     TST1 = TST1 + ABS( H( K+3, K+1 ) )
    // LAPACK❗✔️:                         IF( K.LE.KBOT-4 )
    // LAPACK❗✔️:      $                     TST1 = TST1 + ABS( H( K+4, K+1 ) )
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                      IF( ABS( H( K+1, K ) )
    // LAPACK❗✔️:      $                   .LE.MAX( SMLNUM, ULP*TST1 ) ) THEN
    // LAPACK❗✔️:                         H12 = MAX( ABS( H( K+1, K ) ),
    // LAPACK❗✔️:      $                             ABS( H( K, K+1 ) ) )
    // LAPACK❗✔️:                         H21 = MIN( ABS( H( K+1, K ) ),
    // LAPACK❗✔️:      $                             ABS( H( K, K+1 ) ) )
    // LAPACK❗✔️:                         H11 = MAX( ABS( H( K+1, K+1 ) ),
    // LAPACK❗✔️:      $                             ABS( H( K, K )-H( K+1, K+1 ) ) )
    // LAPACK❗✔️:                         H22 = MIN( ABS( H( K+1, K+1 ) ),
    // LAPACK❗✔️:      $                        ABS( H( K, K )-H( K+1, K+1 ) ) )
    // LAPACK❗✔️:                         SCL = H11 + H12
    // LAPACK❗✔️:                         TST2 = H22*( H11 / SCL )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         IF( TST2.EQ.ZERO .OR. H21*( H12 / SCL ).LE.
    // LAPACK❗✔️:      $                      MAX( SMLNUM, ULP*TST2 ) ) THEN
    // LAPACK❗✔️:                            H( K+1, K ) = ZERO
    // LAPACK❗✔️:                         END IF
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Accumulate orthogonal transformations. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( ACCUM ) THEN
    // LAPACK❗✔️:                   KMS = K - INCOL
    // LAPACK❗✔️:                   T1 = V( 1, M22 )
    // LAPACK❗✔️:                   T2 = T1*V( 2, M22 )
    // LAPACK❗✔️:                   DO 50 J = MAX( 1, KTOP-INCOL ), KDU
    // LAPACK❗✔️:                      REFSUM = U( J, KMS+1 ) + V( 2, M22 )*U( J, KMS+2 )
    // LAPACK❗✔️:                      U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM*T1
    // LAPACK❗✔️:                      U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*T2
    // LAPACK❗✔️:   50                 CONTINUE
    // LAPACK❗✔️:                ELSE IF( WANTZ ) THEN
    // LAPACK❗✔️:                   T1 = V( 1, M22 )
    // LAPACK❗✔️:                   T2 = T1*V( 2, M22 )
    // LAPACK❗✔️:                   DO 60 J = ILOZ, IHIZ
    // LAPACK❗✔️:                      REFSUM = Z( J, K+1 )+V( 2, M22 )*Z( J, K+2 )
    // LAPACK❗✔️:                      Z( J, K+1 ) = Z( J, K+1 ) - REFSUM*T1
    // LAPACK❗✔️:                      Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*T2
    // LAPACK❗✔️:   60              CONTINUE
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Normal case: Chain of 3-by-3 reflections ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 80 M = MBOT, MTOP, -1
    // LAPACK❗✔️:                K = KRCOL + 2*( M-1 )
    // LAPACK❗✔️:                IF( K.EQ.KTOP-1 ) THEN
    // LAPACK❗✔️:                   CALL DLAQR1( 3, H( KTOP, KTOP ), LDH, SR( 2*M-1 ),
    // LAPACK❗✔️:      $                         SI( 2*M-1 ), SR( 2*M ), SI( 2*M ),
    // LAPACK❗✔️:      $                         V( 1, M ) )
    // LAPACK❗✔️:                   ALPHA = V( 1, M )
    // LAPACK❗✔️:                   CALL DLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Perform delayed transformation of row below
    // LAPACK❗✔️: *                 .    Mth bulge. Exploit fact that first two elements
    // LAPACK❗✔️: *                 .    of row are actually zero. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   T1 = V( 1, M )
    // LAPACK❗✔️:                   T2 = T1*V( 2, M )
    // LAPACK❗✔️:                   T3 = T1*V( 3, M )
    // LAPACK❗✔️:                   REFSUM = V( 3, M )*H( K+3, K+2 )
    // LAPACK❗✔️:                   H( K+3, K   ) = -REFSUM*T1
    // LAPACK❗✔️:                   H( K+3, K+1 ) = -REFSUM*T2
    // LAPACK❗✔️:                   H( K+3, K+2 ) = H( K+3, K+2 ) - REFSUM*T3
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Calculate reflection to move
    // LAPACK❗✔️: *                 .    Mth bulge one step. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   BETA      = H( K+1, K )
    // LAPACK❗✔️:                   V( 2, M ) = H( K+2, K )
    // LAPACK❗✔️:                   V( 3, M ) = H( K+3, K )
    // LAPACK❗✔️:                   CALL DLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== A Bulge may collapse because of vigilant
    // LAPACK❗✔️: *                 .    deflation or destructive underflow.  In the
    // LAPACK❗✔️: *                 .    underflow case, try the two-small-subdiagonals
    // LAPACK❗✔️: *                 .    trick to try to reinflate the bulge.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   IF( H( K+3, K ).NE.ZERO .OR. H( K+3, K+1 ).NE.
    // LAPACK❗✔️:      $                ZERO .OR. H( K+3, K+2 ).EQ.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    ==== Typical case: not collapsed (yet). ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      H( K+1, K ) = BETA
    // LAPACK❗✔️:                      H( K+2, K ) = ZERO
    // LAPACK❗✔️:                      H( K+3, K ) = ZERO
    // LAPACK❗✔️:                   ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    ==== Atypical case: collapsed.  Attempt to
    // LAPACK❗✔️: *                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
    // LAPACK❗✔️: *                    .    If the fill resulting from the new
    // LAPACK❗✔️: *                    .    reflector is too large, then abandon it.
    // LAPACK❗✔️: *                    .    Otherwise, use the new one. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      CALL DLAQR1( 3, H( K+1, K+1 ), LDH, SR( 2*M-1 ),
    // LAPACK❗✔️:      $                            SI( 2*M-1 ), SR( 2*M ), SI( 2*M ),
    // LAPACK❗✔️:      $                            VT )
    // LAPACK❗✔️:                      ALPHA = VT( 1 )
    // LAPACK❗✔️:                      CALL DLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
    // LAPACK❗✔️:                      T1 = VT( 1 )
    // LAPACK❗✔️:                      T2 = T1*VT( 2 )
    // LAPACK❗✔️:                      T3 = T1*VT( 3 )
    // LAPACK❗✔️:                      REFSUM = H( K+1, K ) + VT( 2 )*H( K+2, K )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      IF( ABS( H( K+2, K )-REFSUM*T2 )+
    // LAPACK❗✔️:      $                   ABS( REFSUM*T3 ).GT.ULP*
    // LAPACK❗✔️:      $                   ( ABS( H( K, K ) )+ABS( H( K+1,
    // LAPACK❗✔️:      $                   K+1 ) )+ABS( H( K+2, K+2 ) ) ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                       ==== Starting a new bulge here would
    // LAPACK❗✔️: *                       .    create non-negligible fill.  Use
    // LAPACK❗✔️: *                       .    the old one with trepidation. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         H( K+1, K ) = BETA
    // LAPACK❗✔️:                         H( K+2, K ) = ZERO
    // LAPACK❗✔️:                         H( K+3, K ) = ZERO
    // LAPACK❗✔️:                      ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                       ==== Starting a new bulge here would
    // LAPACK❗✔️: *                       .    create only negligible fill.
    // LAPACK❗✔️: *                       .    Replace the old reflector with
    // LAPACK❗✔️: *                       .    the new one. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         H( K+1, K ) = H( K+1, K ) - REFSUM*T1
    // LAPACK❗✔️:                         H( K+2, K ) = ZERO
    // LAPACK❗✔️:                         H( K+3, K ) = ZERO
    // LAPACK❗✔️:                         V( 1, M ) = VT( 1 )
    // LAPACK❗✔️:                         V( 2, M ) = VT( 2 )
    // LAPACK❗✔️:                         V( 3, M ) = VT( 3 )
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ====  Apply reflection from the right and
    // LAPACK❗✔️: *              .     the first column of update from the left.
    // LAPACK❗✔️: *              .     These updates are required for the vigilant
    // LAPACK❗✔️: *              .     deflation check. We still delay most of the
    // LAPACK❗✔️: *              .     updates from the left for efficiency. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                T1 = V( 1, M )
    // LAPACK❗✔️:                T2 = T1*V( 2, M )
    // LAPACK❗✔️:                T3 = T1*V( 3, M )
    // LAPACK❗✔️:                DO 70 J = JTOP, MIN( KBOT, K+3 )
    // LAPACK❗✔️:                   REFSUM = H( J, K+1 ) + V( 2, M )*H( J, K+2 )
    // LAPACK❗✔️:      $                     + V( 3, M )*H( J, K+3 )
    // LAPACK❗✔️:                   H( J, K+1 ) = H( J, K+1 ) - REFSUM*T1
    // LAPACK❗✔️:                   H( J, K+2 ) = H( J, K+2 ) - REFSUM*T2
    // LAPACK❗✔️:                   H( J, K+3 ) = H( J, K+3 ) - REFSUM*T3
    // LAPACK❗✔️:    70          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Perform update from left for subsequent
    // LAPACK❗✔️: *              .    column. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                REFSUM = H( K+1, K+1 ) + V( 2, M )*H( K+2, K+1 )
    // LAPACK❗✔️:      $                  + V( 3, M )*H( K+3, K+1 )
    // LAPACK❗✔️:                H( K+1, K+1 ) = H( K+1, K+1 ) - REFSUM*T1
    // LAPACK❗✔️:                H( K+2, K+1 ) = H( K+2, K+1 ) - REFSUM*T2
    // LAPACK❗✔️:                H( K+3, K+1 ) = H( K+3, K+1 ) - REFSUM*T3
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== The following convergence test requires that
    // LAPACK❗✔️: *              .    the tradition small-compared-to-nearby-diagonals
    // LAPACK❗✔️: *              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
    // LAPACK❗✔️: *              .    criteria both be satisfied.  The latter improves
    // LAPACK❗✔️: *              .    accuracy in some examples. Falling back on an
    // LAPACK❗✔️: *              .    alternate convergence criterion when TST1 or TST2
    // LAPACK❗✔️: *              .    is zero (as done here) is traditional but probably
    // LAPACK❗✔️: *              .    unnecessary. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( K.LT.KTOP)
    // LAPACK❗✔️:      $              CYCLE
    // LAPACK❗✔️:                IF( H( K+1, K ).NE.ZERO ) THEN
    // LAPACK❗✔️:                   TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) )
    // LAPACK❗✔️:                   IF( TST1.EQ.ZERO ) THEN
    // LAPACK❗✔️:                      IF( K.GE.KTOP+1 )
    // LAPACK❗✔️:      $                  TST1 = TST1 + ABS( H( K, K-1 ) )
    // LAPACK❗✔️:                      IF( K.GE.KTOP+2 )
    // LAPACK❗✔️:      $                  TST1 = TST1 + ABS( H( K, K-2 ) )
    // LAPACK❗✔️:                      IF( K.GE.KTOP+3 )
    // LAPACK❗✔️:      $                  TST1 = TST1 + ABS( H( K, K-3 ) )
    // LAPACK❗✔️:                      IF( K.LE.KBOT-2 )
    // LAPACK❗✔️:      $                  TST1 = TST1 + ABS( H( K+2, K+1 ) )
    // LAPACK❗✔️:                      IF( K.LE.KBOT-3 )
    // LAPACK❗✔️:      $                  TST1 = TST1 + ABS( H( K+3, K+1 ) )
    // LAPACK❗✔️:                      IF( K.LE.KBOT-4 )
    // LAPACK❗✔️:      $                  TST1 = TST1 + ABS( H( K+4, K+1 ) )
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                   IF( ABS( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) )
    // LAPACK❗✔️:      $                 THEN
    // LAPACK❗✔️:                      H12 = MAX( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
    // LAPACK❗✔️:                      H21 = MIN( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
    // LAPACK❗✔️:                      H11 = MAX( ABS( H( K+1, K+1 ) ),
    // LAPACK❗✔️:      $                     ABS( H( K, K )-H( K+1, K+1 ) ) )
    // LAPACK❗✔️:                      H22 = MIN( ABS( H( K+1, K+1 ) ),
    // LAPACK❗✔️:      $                     ABS( H( K, K )-H( K+1, K+1 ) ) )
    // LAPACK❗✔️:                      SCL = H11 + H12
    // LAPACK❗✔️:                      TST2 = H22*( H11 / SCL )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      IF( TST2.EQ.ZERO .OR. H21*( H12 / SCL ).LE.
    // LAPACK❗✔️:      $                   MAX( SMLNUM, ULP*TST2 ) ) THEN
    // LAPACK❗✔️:                         H( K+1, K ) = ZERO
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:    80       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Multiply H by reflections from the left ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( ACCUM ) THEN
    // LAPACK❗✔️:                JBOT = MIN( NDCOL, KBOT )
    // LAPACK❗✔️:             ELSE IF( WANTT ) THEN
    // LAPACK❗✔️:                JBOT = N
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                JBOT = KBOT
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 100 M = MBOT, MTOP, -1
    // LAPACK❗✔️:                K = KRCOL + 2*( M-1 )
    // LAPACK❗✔️:                T1 = V( 1, M )
    // LAPACK❗✔️:                T2 = T1*V( 2, M )
    // LAPACK❗✔️:                T3 = T1*V( 3, M )
    // LAPACK❗✔️:                DO 90 J = MAX( KTOP, KRCOL + 2*M ), JBOT
    // LAPACK❗✔️:                   REFSUM = H( K+1, J ) + V( 2, M )*H( K+2, J )
    // LAPACK❗✔️:      $                     + V( 3, M )*H( K+3, J )
    // LAPACK❗✔️:                   H( K+1, J ) = H( K+1, J ) - REFSUM*T1
    // LAPACK❗✔️:                   H( K+2, J ) = H( K+2, J ) - REFSUM*T2
    // LAPACK❗✔️:                   H( K+3, J ) = H( K+3, J ) - REFSUM*T3
    // LAPACK❗✔️:    90          CONTINUE
    // LAPACK❗✔️:   100       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Accumulate orthogonal transformations. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( ACCUM ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Accumulate U. (If needed, update Z later
    // LAPACK❗✔️: *              .    with an efficient matrix-matrix
    // LAPACK❗✔️: *              .    multiply.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 120 M = MBOT, MTOP, -1
    // LAPACK❗✔️:                   K = KRCOL + 2*( M-1 )
    // LAPACK❗✔️:                   KMS = K - INCOL
    // LAPACK❗✔️:                   I2 = MAX( 1, KTOP-INCOL )
    // LAPACK❗✔️:                   I2 = MAX( I2, KMS-(KRCOL-INCOL)+1 )
    // LAPACK❗✔️:                   I4 = MIN( KDU, KRCOL + 2*( MBOT-1 ) - INCOL + 5 )
    // LAPACK❗✔️:                   T1 = V( 1, M )
    // LAPACK❗✔️:                   T2 = T1*V( 2, M )
    // LAPACK❗✔️:                   T3 = T1*V( 3, M )
    // LAPACK❗✔️:                   DO 110 J = I2, I4
    // LAPACK❗✔️:                      REFSUM = U( J, KMS+1 ) + V( 2, M )*U( J, KMS+2 )
    // LAPACK❗✔️:      $                        + V( 3, M )*U( J, KMS+3 )
    // LAPACK❗✔️:                      U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM*T1
    // LAPACK❗✔️:                      U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*T2
    // LAPACK❗✔️:                      U( J, KMS+3 ) = U( J, KMS+3 ) - REFSUM*T3
    // LAPACK❗✔️:   110             CONTINUE
    // LAPACK❗✔️:   120          CONTINUE
    // LAPACK❗✔️:             ELSE IF( WANTZ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== U is not accumulated, so update Z
    // LAPACK❗✔️: *              .    now by multiplying by reflections
    // LAPACK❗✔️: *              .    from the right. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 140 M = MBOT, MTOP, -1
    // LAPACK❗✔️:                   K = KRCOL + 2*( M-1 )
    // LAPACK❗✔️:                   T1 = V( 1, M )
    // LAPACK❗✔️:                   T2 = T1*V( 2, M )
    // LAPACK❗✔️:                   T3 = T1*V( 3, M )
    // LAPACK❗✔️:                   DO 130 J = ILOZ, IHIZ
    // LAPACK❗✔️:                      REFSUM = Z( J, K+1 ) + V( 2, M )*Z( J, K+2 )
    // LAPACK❗✔️:      $                        + V( 3, M )*Z( J, K+3 )
    // LAPACK❗✔️:                      Z( J, K+1 ) = Z( J, K+1 ) - REFSUM*T1
    // LAPACK❗✔️:                      Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*T2
    // LAPACK❗✔️:                      Z( J, K+3 ) = Z( J, K+3 ) - REFSUM*T3
    // LAPACK❗✔️:   130             CONTINUE
    // LAPACK❗✔️:   140          CONTINUE
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== End of near-the-diagonal bulge chase. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   145    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Use U (if accumulated) to update far-from-diagonal
    // LAPACK❗✔️: *        .    entries in H.  If required, use U to update Z as
    // LAPACK❗✔️: *        .    well. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( ACCUM ) THEN
    // LAPACK❗✔️:             IF( WANTT ) THEN
    // LAPACK❗✔️:                JTOP = 1
    // LAPACK❗✔️:                JBOT = N
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                JTOP = KTOP
    // LAPACK❗✔️:                JBOT = KBOT
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             K1 = MAX( 1, KTOP-INCOL )
    // LAPACK❗✔️:             NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Horizontal Multiply ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 150 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
    // LAPACK❗✔️:                JLEN = MIN( NH, JBOT-JCOL+1 )
    // LAPACK❗✔️:                CALL DGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ),
    // LAPACK❗✔️:      $                        LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH,
    // LAPACK❗✔️:      $                        LDWH )
    // LAPACK❗✔️:                CALL DLACPY( 'ALL', NU, JLEN, WH, LDWH,
    // LAPACK❗✔️:      $                         H( INCOL+K1, JCOL ), LDH )
    // LAPACK❗✔️:   150       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Vertical multiply ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 160 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
    // LAPACK❗✔️:                JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
    // LAPACK❗✔️:                CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE,
    // LAPACK❗✔️:      $                     H( JROW, INCOL+K1 ), LDH, U( K1, K1 ),
    // LAPACK❗✔️:      $                     LDU, ZERO, WV, LDWV )
    // LAPACK❗✔️:                CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV,
    // LAPACK❗✔️:      $                      H( JROW, INCOL+K1 ), LDH )
    // LAPACK❗✔️:   160       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Z multiply (also vertical) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( WANTZ ) THEN
    // LAPACK❗✔️:                DO 170 JROW = ILOZ, IHIZ, NV
    // LAPACK❗✔️:                   JLEN = MIN( NV, IHIZ-JROW+1 )
    // LAPACK❗✔️:                   CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE,
    // LAPACK❗✔️:      $                        Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ),
    // LAPACK❗✔️:      $                        LDU, ZERO, WV, LDWV )
    // LAPACK❗✔️:                   CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV,
    // LAPACK❗✔️:      $                         Z( JROW, INCOL+K1 ), LDZ )
    // LAPACK❗✔️:   170          CONTINUE
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:   180 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== End of DLAQR5 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if nshifts < 2 || ktop >= kbot {
        return Ok(());
    }
    for i in (0..nshifts - 2).step_by(2) {
        if si[i] != -si[i + 1] {
            let swap = sr[i];
            sr[i] = sr[i + 1];
            sr[i + 1] = sr[i + 2];
            sr[i + 2] = swap;
            let swap = si[i];
            si[i] = si[i + 1];
            si[i + 1] = si[i + 2];
            si[i + 2] = swap;
        }
    }
    let ns = nshifts - nshifts % 2;
    let safmin = dlamch(b'S');
    let ulp = dlamch(b'P');
    let smlnum = safmin * (n as f64 / ulp);
    let accum = kacc22 == 1 || kacc22 == 2;
    let (ktop, kbot, n) = (ktop as isize, kbot as isize, n as isize);
    let nbmps = (ns / 2) as isize;
    let kdu = 4 * nbmps;
    macro_rules! h {
        ($i:expr,$j:expr) => {
            h[at($i, $j, ldh)]
        };
    }
    macro_rules! v {
        ($i:expr,$j:expr) => {
            v[at($i, $j, ldv)]
        };
    }
    macro_rules! u {
        ($i:expr,$j:expr) => {
            h[uoff + at($i, $j, ldu)]
        };
    }
    if ktop + 2 <= kbot {
        h!(ktop + 2, ktop) = 0.0;
    }
    let mut incol = ktop - 2 * nbmps + 1;
    while incol <= kbot - 2 {
        let mut jtop = if accum {
            ktop.max(incol)
        } else if wantt {
            1
        } else {
            ktop
        };
        let ndcol = incol + kdu;
        if accum {
            for j in 1..=kdu {
                for i in 1..=kdu {
                    u!(i, j) = if i == j { 1.0 } else { 0.0 };
                }
            }
        }
        for krcol in incol..=(incol + 2 * nbmps - 1).min(kbot - 2) {
            let mtop = 1.max((ktop - krcol) / 2 + 1);
            let mbot = nbmps.min((kbot - krcol - 1) / 2);
            let m22 = mbot + 1;
            let bmp22 = mbot < nbmps && krcol + 2 * (m22 - 1) == kbot - 2;
            if bmp22 {
                let k = krcol + 2 * (m22 - 1);
                let mut tmp = [0.0; 2];
                if k == ktop - 1 {
                    let init = laqr1(
                        2,
                        &h[at(k + 1, k + 1, ldh)..],
                        ldh,
                        sr[(2 * m22 - 2) as usize],
                        si[(2 * m22 - 2) as usize],
                        sr[(2 * m22 - 1) as usize],
                        si[(2 * m22 - 1) as usize],
                    );
                    tmp.copy_from_slice(&init[..2]);
                    v!(1, m22) = larfg(&mut tmp);
                    v!(2, m22) = tmp[1];
                } else {
                    tmp = [h!(k + 1, k), h!(k + 2, k)];
                    v!(1, m22) = larfg(&mut tmp);
                    v!(2, m22) = tmp[1];
                    h!(k + 1, k) = tmp[0];
                    h!(k + 2, k) = 0.0;
                }
                let t1 = v!(1, m22);
                let t2 = t1 * v!(2, m22);
                for j in jtop..=kbot.min(k + 3) {
                    let refsum = h!(j, k + 1) + v!(2, m22) * h!(j, k + 2);
                    h!(j, k + 1) -= refsum * t1;
                    h!(j, k + 2) -= refsum * t2;
                }
                let jbot = if accum {
                    ndcol.min(kbot)
                } else if wantt {
                    n
                } else {
                    kbot
                };
                for j in k + 1..=jbot {
                    let refsum = h!(k + 1, j) + v!(2, m22) * h!(k + 2, j);
                    h!(k + 1, j) -= refsum * t1;
                    h!(k + 2, j) -= refsum * t2;
                }
                vigilant_deflation(h, ldh, k, ktop, kbot, smlnum, ulp);
                if accum {
                    let kms = k - incol;
                    for j in 1.max(ktop - incol)..=kdu {
                        let refsum = u!(j, kms + 1) + v!(2, m22) * u!(j, kms + 2);
                        u!(j, kms + 1) -= refsum * t1;
                        u!(j, kms + 2) -= refsum * t2;
                    }
                } else if wantz {
                    for j in iloz as isize..=ihiz as isize {
                        let refsum = z[at(j, k + 1, ldz)] + v!(2, m22) * z[at(j, k + 2, ldz)];
                        z[at(j, k + 1, ldz)] -= refsum * t1;
                        z[at(j, k + 2, ldz)] -= refsum * t2;
                    }
                }
            }
            for m in (mtop..=mbot).rev() {
                let k = krcol + 2 * (m - 1);
                if k == ktop - 1 {
                    let mut tmp = laqr1(
                        3,
                        &h[at(ktop, ktop, ldh)..],
                        ldh,
                        sr[(2 * m - 2) as usize],
                        si[(2 * m - 2) as usize],
                        sr[(2 * m - 1) as usize],
                        si[(2 * m - 1) as usize],
                    );
                    v!(1, m) = larfg(&mut tmp);
                    v!(2, m) = tmp[1];
                    v!(3, m) = tmp[2];
                } else {
                    let t1 = v!(1, m);
                    let t2 = t1 * v!(2, m);
                    let t3 = t1 * v!(3, m);
                    let refsum = v!(3, m) * h!(k + 3, k + 2);
                    h!(k + 3, k) = -refsum * t1;
                    h!(k + 3, k + 1) = -refsum * t2;
                    h!(k + 3, k + 2) -= refsum * t3;
                    let mut tmp = [h!(k + 1, k), h!(k + 2, k), h!(k + 3, k)];
                    v!(1, m) = larfg(&mut tmp);
                    v!(2, m) = tmp[1];
                    v!(3, m) = tmp[2];
                    let beta = tmp[0];
                    if h!(k + 3, k) != 0.0 || h!(k + 3, k + 1) != 0.0 || h!(k + 3, k + 2) == 0.0 {
                        h!(k + 1, k) = beta;
                        h!(k + 2, k) = 0.0;
                        h!(k + 3, k) = 0.0;
                    } else {
                        let mut vt = laqr1(
                            3,
                            &h[at(k + 1, k + 1, ldh)..],
                            ldh,
                            sr[(2 * m - 2) as usize],
                            si[(2 * m - 2) as usize],
                            sr[(2 * m - 1) as usize],
                            si[(2 * m - 1) as usize],
                        );
                        let tau = larfg(&mut vt);
                        vt[0] = tau;
                        let t1 = vt[0];
                        let t2 = t1 * vt[1];
                        let t3 = t1 * vt[2];
                        let refsum = h!(k + 1, k) + vt[1] * h!(k + 2, k);
                        if (h!(k + 2, k) - refsum * t2).abs() + (refsum * t3).abs()
                            > ulp
                                * (h!(k, k).abs() + h!(k + 1, k + 1).abs() + h!(k + 2, k + 2).abs())
                        {
                            h!(k + 1, k) = beta;
                            h!(k + 2, k) = 0.0;
                            h!(k + 3, k) = 0.0;
                        } else {
                            h!(k + 1, k) -= refsum * t1;
                            h!(k + 2, k) = 0.0;
                            h!(k + 3, k) = 0.0;
                            v!(1, m) = vt[0];
                            v!(2, m) = vt[1];
                            v!(3, m) = vt[2];
                        }
                    }
                }
                let t1 = v!(1, m);
                let t2 = t1 * v!(2, m);
                let t3 = t1 * v!(3, m);
                for j in jtop..=kbot.min(k + 3) {
                    let refsum = h!(j, k + 1) + v!(2, m) * h!(j, k + 2) + v!(3, m) * h!(j, k + 3);
                    h!(j, k + 1) -= refsum * t1;
                    h!(j, k + 2) -= refsum * t2;
                    h!(j, k + 3) -= refsum * t3;
                }
                let refsum =
                    h!(k + 1, k + 1) + v!(2, m) * h!(k + 2, k + 1) + v!(3, m) * h!(k + 3, k + 1);
                h!(k + 1, k + 1) -= refsum * t1;
                h!(k + 2, k + 1) -= refsum * t2;
                h!(k + 3, k + 1) -= refsum * t3;
                vigilant_deflation(h, ldh, k, ktop, kbot, smlnum, ulp);
            }
            let jbot = if accum {
                ndcol.min(kbot)
            } else if wantt {
                n
            } else {
                kbot
            };
            for m in (mtop..=mbot).rev() {
                let k = krcol + 2 * (m - 1);
                let t1 = v!(1, m);
                let t2 = t1 * v!(2, m);
                let t3 = t1 * v!(3, m);
                for j in ktop.max(krcol + 2 * m)..=jbot {
                    let refsum = h!(k + 1, j) + v!(2, m) * h!(k + 2, j) + v!(3, m) * h!(k + 3, j);
                    h!(k + 1, j) -= refsum * t1;
                    h!(k + 2, j) -= refsum * t2;
                    h!(k + 3, j) -= refsum * t3;
                }
            }
            if accum {
                for m in (mtop..=mbot).rev() {
                    let k = krcol + 2 * (m - 1);
                    let kms = k - incol;
                    let i2 = 1.max(ktop - incol).max(kms - (krcol - incol) + 1);
                    let i4 = kdu.min(krcol + 2 * (mbot - 1) - incol + 5);
                    let t1 = v!(1, m);
                    let t2 = t1 * v!(2, m);
                    let t3 = t1 * v!(3, m);
                    for j in i2..=i4 {
                        let refsum =
                            u!(j, kms + 1) + v!(2, m) * u!(j, kms + 2) + v!(3, m) * u!(j, kms + 3);
                        u!(j, kms + 1) -= refsum * t1;
                        u!(j, kms + 2) -= refsum * t2;
                        u!(j, kms + 3) -= refsum * t3;
                    }
                }
            } else if wantz {
                for m in (mtop..=mbot).rev() {
                    let k = krcol + 2 * (m - 1);
                    let t1 = v!(1, m);
                    let t2 = t1 * v!(2, m);
                    let t3 = t1 * v!(3, m);
                    for j in iloz as isize..=ihiz as isize {
                        let refsum = z[at(j, k + 1, ldz)]
                            + v!(2, m) * z[at(j, k + 2, ldz)]
                            + v!(3, m) * z[at(j, k + 3, ldz)];
                        z[at(j, k + 1, ldz)] -= refsum * t1;
                        z[at(j, k + 2, ldz)] -= refsum * t2;
                        z[at(j, k + 3, ldz)] -= refsum * t3;
                    }
                }
            }
        }
        if accum {
            jtop = if wantt { 1 } else { ktop };
            let jbot = if wantt { n } else { kbot };
            let k1 = 1.max(ktop - incol);
            let nu = (kdu - 0.max(ndcol - kbot)) - k1 + 1;
            let mut jcol = ndcol.min(kbot) + 1;
            while jcol <= jbot {
                let jlen = (nh as isize).min(jbot - jcol + 1);
                blas3::gemm_within(
                    Trans,
                    NoTrans,
                    nu as usize,
                    jlen as usize,
                    nu as usize,
                    1.0,
                    h,
                    uoff + at(k1, k1, ldu),
                    ldu,
                    at(incol + k1, jcol, ldh),
                    ldh,
                    0.0,
                    whoff,
                    ldwh,
                )?;
                for j in 0..jlen {
                    for i in 0..nu {
                        h[at(incol + k1 + i, jcol + j, ldh)] = h[whoff + at(i + 1, j + 1, ldwh)];
                    }
                }
                jcol += nh as isize;
            }
            let mut jrow = jtop;
            while jrow < ktop.max(incol) {
                let jlen = (nv as isize).min(ktop.max(incol) - jrow);
                blas3::gemm_within(
                    NoTrans,
                    NoTrans,
                    jlen as usize,
                    nu as usize,
                    nu as usize,
                    1.0,
                    h,
                    at(jrow, incol + k1, ldh),
                    ldh,
                    uoff + at(k1, k1, ldu),
                    ldu,
                    0.0,
                    wvoff,
                    ldwv,
                )?;
                for j in 0..nu {
                    for i in 0..jlen {
                        h[at(jrow + i, incol + k1 + j, ldh)] = h[wvoff + at(i + 1, j + 1, ldwv)];
                    }
                }
                jrow += nv as isize;
            }
            if wantz {
                let mut jrow = iloz as isize;
                while jrow <= ihiz as isize {
                    let jlen = (nv as isize).min(ihiz as isize - jrow + 1);
                    blas3::gemm_into_within(
                        NoTrans,
                        NoTrans,
                        jlen as usize,
                        nu as usize,
                        nu as usize,
                        1.0,
                        h,
                        Separate(&z[at(jrow, incol + k1, ldz)..]),
                        ldz,
                        Within(uoff + at(k1, k1, ldu)),
                        ldu,
                        0.0,
                        wvoff,
                        ldwv,
                    )?;
                    for j in 0..nu {
                        for i in 0..jlen {
                            z[at(jrow + i, incol + k1 + j, ldz)] =
                                h[wvoff + at(i + 1, j + 1, ldwv)];
                        }
                    }
                    jrow += nv as isize;
                }
            }
        }
        incol += 2 * nbmps;
    }
    Ok(())
}

/// The IPARMQ constants and fixed real-LAPACK driver names used here.
fn recommended_shifts(nh: usize) -> usize {
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22, ICOST
    // LAPACK❗✔️:       PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
    // LAPACK❗✔️:      $                   ISHFTS = 15, IACC22 = 16, ICOST = 17 )
    // LAPACK❗✔️:       INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP, RCOST
    // LAPACK❗✔️:       PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
    // LAPACK❗✔️:      $                   NIBBLE = 14, KNWSWP = 500, RCOST = 10 )
    // LAPACK❗✔️:       REAL               TWO
    // LAPACK❗✔️:       PARAMETER          ( TWO = 2.0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            NH, NS
    // LAPACK❗✔️:       INTEGER            I, IC, IZ
    // LAPACK❗✔️:       CHARACTER          SUBNAM*6
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          LOG, MAX, MOD, NINT, REAL
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️:       IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
    // LAPACK❗✔️:      $    ( ISPEC.EQ.IACC22 ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Set the number simultaneous shifts ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NH = IHI - ILO + 1
    // LAPACK❗✔️:          NS = 2
    // LAPACK❗✔️:          IF( NH.GE.30 )
    // LAPACK❗✔️:      $      NS = 4
    // LAPACK❗✔️:          IF( NH.GE.60 )
    // LAPACK❗✔️:      $      NS = 10
    // LAPACK❗✔️:          IF( NH.GE.150 )
    // LAPACK❗✔️:      $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
    // LAPACK❗✔️:          IF( NH.GE.590 )
    // LAPACK❗✔️:      $      NS = 64
    // LAPACK❗✔️:          IF( NH.GE.3000 )
    // LAPACK❗✔️:      $      NS = 128
    // LAPACK❗✔️:          IF( NH.GE.6000 )
    // LAPACK❗✔️:      $      NS = 256
    // LAPACK❗✔️:          NS = MAX( 2, NS-MOD( NS, 2 ) )
    // LAPACK❗✔️:       END IF
    let mut ns = 2;
    if nh >= 30 {
        ns = 4;
    }
    if nh >= 60 {
        ns = 10;
    }
    if nh >= 150 {
        ns = 10.max(nh / ((nh as f32).ln() / 2.0_f32.ln()).round() as usize);
    }
    if nh >= 590 {
        ns = 64;
    }
    if nh >= 3000 {
        ns = 128;
    }
    if nh >= 6000 {
        ns = 256;
    }
    2.max(ns - ns % 2)
}
fn aed_work(recursive: bool, jw: usize) -> usize {
    if jw <= 2 {
        return 1;
    }
    let lwk1 = jw * 32 + 65 * 64; // DGEHRD(JW,1,JW-1), preferred workspace.
    let lwk2 = super::hessenberg::ormhr_work(jw);
    let result = jw + lwk1.max(lwk2);
    if recursive {
        result.max(qr_work(false, jw, 1, jw))
    } else {
        result
    }
}
pub(super) fn qr_work(recursive: bool, n: usize, ilo: usize, ihi: usize) -> usize {
    if n <= 15 {
        return 1;
    }
    let nh = ihi - ilo + 1;
    let ns = recommended_shifts(nh);
    let nwr = 2
        .max(if nh <= 500 { ns } else { 3 * ns / 2 })
        .min(nh)
        .min((n - 1) / 3);
    let mut nsr = ns.min((n - 3) / 6).min(ihi - ilo);
    nsr = 2.max(nsr - nsr % 2);
    (3 * nsr / 2).max(aed_work(recursive, (nwr + 1).min(nh)))
}

/// DLAQR3 when recursive=true; DLAQR2 otherwise. Their shared state updates
/// are identical; only the source-selected window QR/workspace call differs.
#[allow(clippy::too_many_arguments)]
fn aed(
    recursive: bool,
    wantt: bool,
    wantz: bool,
    n: usize,
    ktop: usize,
    kbot: usize,
    nw: usize,
    h: &mut [f64],
    ldh: usize,
    iloz: usize,
    ihiz: usize,
    z: &mut [f64],
    ldz: usize,
    sr: &mut [f64],
    si: &mut [f64],
    voff: usize,
    ldv: usize,
    nh: usize,
    toff: usize,
    ldt: usize,
    nv: usize,
    wvoff: usize,
    ldwv: usize,
    work: &mut [f64],
) -> Result<(usize, usize), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH,
    // LAPACK❗✔️:      $                   ILOZ,
    // LAPACK❗✔️:      $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
    // LAPACK❗✔️:      $                   LDT, NV, WV, LDWV, WORK, LWORK )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
    // LAPACK❗✔️:      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
    // LAPACK❗✔️:       LOGICAL            WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
    // LAPACK❗✔️:      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
    // LAPACK❗✔️:      $                   Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  ================================================================
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
    // LAPACK❗✔️:      $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
    // LAPACK❗✔️:       INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
    // LAPACK❗✔️:      $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,
    // LAPACK❗✔️:      $                   LWKOPT, NMIN
    // LAPACK❗✔️:       LOGICAL            BULGE, SORTED
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           DLAMCH, ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DCOPY, DGEHRD, DGEMM, DLACPY, DLAHQR,
    // LAPACK❗✔️:      $                   DLANV2,
    // LAPACK❗✔️:      $                   DLAQR4, DLARF1F, DLARFG, DLASET, DORMHR, DTREXC
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Estimate optimal workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       JW = MIN( NW, KBOT-KTOP+1 )
    // LAPACK❗✔️:       IF( JW.LE.2 ) THEN
    // LAPACK❗✔️:          LWKOPT = 1
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DGEHRD ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
    // LAPACK❗✔️:          LWK1 = INT( WORK( 1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DORMHR ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DORMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V,
    // LAPACK❗✔️:      $                LDV,
    // LAPACK❗✔️:      $                WORK, -1, INFO )
    // LAPACK❗✔️:          LWK2 = INT( WORK( 1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DLAQR4 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLAQR4( .true., .true., JW, 1, JW, T, LDT, SR, SI, 1,
    // LAPACK❗✔️:      $                JW,
    // LAPACK❗✔️:      $                V, LDV, WORK, -1, INFQR )
    // LAPACK❗✔️:          LWK3 = INT( WORK( 1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Optimal workspace ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Quick return in case of workspace query. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LWORK.EQ.-1 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Nothing to do ...
    // LAPACK❗✔️: *     ... for an empty active block ... ====
    // LAPACK❗✔️:       NS = 0
    // LAPACK❗✔️:       ND = 0
    // LAPACK❗✔️:       WORK( 1 ) = ONE
    // LAPACK❗✔️:       IF( KTOP.GT.KBOT )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *     ... nor for an empty deflation window. ====
    // LAPACK❗✔️:       IF( NW.LT.1 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Machine constants ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'SAFE MINIMUM' )
    // LAPACK❗✔️:       SAFMAX = ONE / SAFMIN
    // LAPACK❗✔️:       ULP = DLAMCH( 'PRECISION' )
    // LAPACK❗✔️:       SMLNUM = SAFMIN*( DBLE( N ) / ULP )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Setup deflation window ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       JW = MIN( NW, KBOT-KTOP+1 )
    // LAPACK❗✔️:       KWTOP = KBOT - JW + 1
    // LAPACK❗✔️:       IF( KWTOP.EQ.KTOP ) THEN
    // LAPACK❗✔️:          S = ZERO
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          S = H( KWTOP, KWTOP-1 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( KBOT.EQ.KWTOP ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== 1-by-1 deflation window: not much to do ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SR( KWTOP ) = H( KWTOP, KWTOP )
    // LAPACK❗✔️:          SI( KWTOP ) = ZERO
    // LAPACK❗✔️:          NS = 1
    // LAPACK❗✔️:          ND = 0
    // LAPACK❗✔️:          IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) )
    // LAPACK❗✔️:      $        THEN
    // LAPACK❗✔️:             NS = 0
    // LAPACK❗✔️:             ND = 1
    // LAPACK❗✔️:             IF( KWTOP.GT.KTOP )
    // LAPACK❗✔️:      $         H( KWTOP, KWTOP-1 ) = ZERO
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          WORK( 1 ) = ONE
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Convert to spike-triangular form.  (In case of a
    // LAPACK❗✔️: *     .    rare QR failure, this routine continues to do
    // LAPACK❗✔️: *     .    aggressive early deflation using that part of
    // LAPACK❗✔️: *     .    the deflation window that converged using INFQR
    // LAPACK❗✔️: *     .    here and there to keep track.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
    // LAPACK❗✔️:       CALL DCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ),
    // LAPACK❗✔️:      $            LDT+1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
    // LAPACK❗✔️:       NMIN = ILAENV( 12, 'DLAQR3', 'SV', JW, 1, JW, LWORK )
    // LAPACK❗✔️:       IF( JW.GT.NMIN ) THEN
    // LAPACK❗✔️:          CALL DLAQR4( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ),
    // LAPACK❗✔️:      $                SI( KWTOP ), 1, JW, V, LDV, WORK, LWORK, INFQR )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          CALL DLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ),
    // LAPACK❗✔️:      $                SI( KWTOP ), 1, JW, V, LDV, INFQR )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== DTREXC needs a clean margin near the diagonal ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 J = 1, JW - 3
    // LAPACK❗✔️:          T( J+2, J ) = ZERO
    // LAPACK❗✔️:          T( J+3, J ) = ZERO
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       IF( JW.GT.2 )
    // LAPACK❗✔️:      $   T( JW, JW-2 ) = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Deflation detection loop ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NS = JW
    // LAPACK❗✔️:       ILST = INFQR + 1
    // LAPACK❗✔️:    20 CONTINUE
    // LAPACK❗✔️:       IF( ILST.LE.NS ) THEN
    // LAPACK❗✔️:          IF( NS.EQ.1 ) THEN
    // LAPACK❗✔️:             BULGE = .FALSE.
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             BULGE = T( NS, NS-1 ).NE.ZERO
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Small spike tip test for deflation ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( .NOT. BULGE ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Real eigenvalue ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             FOO = ABS( T( NS, NS ) )
    // LAPACK❗✔️:             IF( FOO.EQ.ZERO )
    // LAPACK❗✔️:      $         FOO = ABS( S )
    // LAPACK❗✔️:             IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Deflatable ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = NS - 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Undeflatable.   Move it up out of the way.
    // LAPACK❗✔️: *              .    (DTREXC can not fail in this case.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IFST = NS
    // LAPACK❗✔️:                CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST,
    // LAPACK❗✔️:      $                      WORK,
    // LAPACK❗✔️:      $                      INFO )
    // LAPACK❗✔️:                ILST = ILST + 1
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Complex conjugate pair ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )*
    // LAPACK❗✔️:      $            SQRT( ABS( T( NS-1, NS ) ) )
    // LAPACK❗✔️:             IF( FOO.EQ.ZERO )
    // LAPACK❗✔️:      $         FOO = ABS( S )
    // LAPACK❗✔️:             IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE.
    // LAPACK❗✔️:      $          MAX( SMLNUM, ULP*FOO ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Deflatable ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = NS - 2
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Undeflatable. Move them up out of the way.
    // LAPACK❗✔️: *              .    Fortunately, DTREXC does the right thing with
    // LAPACK❗✔️: *              .    ILST in case of a rare exchange failure. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IFST = NS
    // LAPACK❗✔️:                CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST,
    // LAPACK❗✔️:      $                      WORK,
    // LAPACK❗✔️:      $                      INFO )
    // LAPACK❗✔️:                ILST = ILST + 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== End deflation detection loop ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          GO TO 20
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Return to Hessenberg form ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NS.EQ.0 )
    // LAPACK❗✔️:      $   S = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NS.LT.JW ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== sorting diagonal blocks of T improves accuracy for
    // LAPACK❗✔️: *        .    graded matrices.  Bubble sort deals well with
    // LAPACK❗✔️: *        .    exchange failures. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SORTED = .false.
    // LAPACK❗✔️:          I = NS + 1
    // LAPACK❗✔️:    30    CONTINUE
    // LAPACK❗✔️:          IF( SORTED )
    // LAPACK❗✔️:      $      GO TO 50
    // LAPACK❗✔️:          SORTED = .true.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KEND = I - 1
    // LAPACK❗✔️:          I = INFQR + 1
    // LAPACK❗✔️:          IF( I.EQ.NS ) THEN
    // LAPACK❗✔️:             K = I + 1
    // LAPACK❗✔️:          ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
    // LAPACK❗✔️:             K = I + 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             K = I + 2
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️:          IF( K.LE.KEND ) THEN
    // LAPACK❗✔️:             IF( K.EQ.I+1 ) THEN
    // LAPACK❗✔️:                EVI = ABS( T( I, I ) )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )*
    // LAPACK❗✔️:      $               SQRT( ABS( T( I, I+1 ) ) )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( K.EQ.KEND ) THEN
    // LAPACK❗✔️:                EVK = ABS( T( K, K ) )
    // LAPACK❗✔️:             ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                EVK = ABS( T( K, K ) )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )*
    // LAPACK❗✔️:      $               SQRT( ABS( T( K, K+1 ) ) )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( EVI.GE.EVK ) THEN
    // LAPACK❗✔️:                I = K
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                SORTED = .false.
    // LAPACK❗✔️:                IFST = I
    // LAPACK❗✔️:                ILST = K
    // LAPACK❗✔️:                CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST,
    // LAPACK❗✔️:      $                      WORK,
    // LAPACK❗✔️:      $                      INFO )
    // LAPACK❗✔️:                IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️:                   I = ILST
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   I = K
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( I.EQ.KEND ) THEN
    // LAPACK❗✔️:                K = I + 1
    // LAPACK❗✔️:             ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                K = I + 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                K = I + 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             GO TO 40
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          GO TO 30
    // LAPACK❗✔️:    50    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Restore shift/eigenvalue array from T ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       I = JW
    // LAPACK❗✔️:    60 CONTINUE
    // LAPACK❗✔️:       IF( I.GE.INFQR+1 ) THEN
    // LAPACK❗✔️:          IF( I.EQ.INFQR+1 ) THEN
    // LAPACK❗✔️:             SR( KWTOP+I-1 ) = T( I, I )
    // LAPACK❗✔️:             SI( KWTOP+I-1 ) = ZERO
    // LAPACK❗✔️:             I = I - 1
    // LAPACK❗✔️:          ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
    // LAPACK❗✔️:             SR( KWTOP+I-1 ) = T( I, I )
    // LAPACK❗✔️:             SI( KWTOP+I-1 ) = ZERO
    // LAPACK❗✔️:             I = I - 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             AA = T( I-1, I-1 )
    // LAPACK❗✔️:             CC = T( I, I-1 )
    // LAPACK❗✔️:             BB = T( I-1, I )
    // LAPACK❗✔️:             DD = T( I, I )
    // LAPACK❗✔️:             CALL DLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ),
    // LAPACK❗✔️:      $                   SI( KWTOP+I-2 ), SR( KWTOP+I-1 ),
    // LAPACK❗✔️:      $                   SI( KWTOP+I-1 ), CS, SN )
    // LAPACK❗✔️:             I = I - 2
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          GO TO 60
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
    // LAPACK❗✔️:          IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Reflect spike back into lower triangle ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DCOPY( NS, V, LDV, WORK, 1 )
    // LAPACK❗✔️:             BETA = WORK( 1 )
    // LAPACK❗✔️:             CALL DLARFG( NS, BETA, WORK( 2 ), 1, TAU )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ),
    // LAPACK❗✔️:      $                   LDT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARF1F( 'L', NS, JW, WORK, 1, TAU, T, LDT,
    // LAPACK❗✔️:      $                  WORK( JW+1 ) )
    // LAPACK❗✔️:             CALL DLARF1F( 'R', NS, NS, WORK, 1, TAU, T, LDT,
    // LAPACK❗✔️:      $                  WORK( JW+1 ) )
    // LAPACK❗✔️:             CALL DLARF1F( 'R', JW, NS, WORK, 1, TAU, V, LDV,
    // LAPACK❗✔️:      $                  WORK( JW+1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
    // LAPACK❗✔️:      $                   LWORK-JW, INFO )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Copy updated reduced window into place ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( KWTOP.GT.1 )
    // LAPACK❗✔️:      $      H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
    // LAPACK❗✔️:          CALL DLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
    // LAPACK❗✔️:          CALL DCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),
    // LAPACK❗✔️:      $               LDH+1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Accumulate orthogonal matrix in order update
    // LAPACK❗✔️: *        .    H and Z, if requested.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( NS.GT.1 .AND. S.NE.ZERO )
    // LAPACK❗✔️:      $      CALL DORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V,
    // LAPACK❗✔️:      $                   LDV,
    // LAPACK❗✔️:      $                   WORK( JW+1 ), LWORK-JW, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Update vertical slab in H ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️:             LTOP = 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             LTOP = KTOP
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          DO 70 KROW = LTOP, KWTOP - 1, NV
    // LAPACK❗✔️:             KLN = MIN( NV, KWTOP-KROW )
    // LAPACK❗✔️:             CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),
    // LAPACK❗✔️:      $                  LDH, V, LDV, ZERO, WV, LDWV )
    // LAPACK❗✔️:             CALL DLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ),
    // LAPACK❗✔️:      $                   LDH )
    // LAPACK❗✔️:    70    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Update horizontal slab in H ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️:             DO 80 KCOL = KBOT + 1, N, NH
    // LAPACK❗✔️:                KLN = MIN( NH, N-KCOL+1 )
    // LAPACK❗✔️:                CALL DGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,
    // LAPACK❗✔️:      $                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
    // LAPACK❗✔️:                CALL DLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),
    // LAPACK❗✔️:      $                      LDH )
    // LAPACK❗✔️:    80       CONTINUE
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Update vertical slab in Z ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTZ ) THEN
    // LAPACK❗✔️:             DO 90 KROW = ILOZ, IHIZ, NV
    // LAPACK❗✔️:                KLN = MIN( NV, IHIZ-KROW+1 )
    // LAPACK❗✔️:                CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW,
    // LAPACK❗✔️:      $                     KWTOP ),
    // LAPACK❗✔️:      $                     LDZ, V, LDV, ZERO, WV, LDWV )
    // LAPACK❗✔️:                CALL DLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),
    // LAPACK❗✔️:      $                      LDZ )
    // LAPACK❗✔️:    90       CONTINUE
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Return the number of deflations ... ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ND = JW - NS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== ... and the number of shifts. (Subtracting
    // LAPACK❗✔️: *     .    INFQR from the spike length takes care
    // LAPACK❗✔️: *     .    of the case of a rare QR failure while
    // LAPACK❗✔️: *     .    calculating eigenvalues of the deflation
    // LAPACK❗✔️: *     .    window.)  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NS = NS - INFQR
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *      ==== Return optimal workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== End of DLAQR3 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    // LAPACK❗✔️:       SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH,
    // LAPACK❗✔️:      $                   ILOZ,
    // LAPACK❗✔️:      $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
    // LAPACK❗✔️:      $                   LDT, NV, WV, LDWV, WORK, LWORK )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
    // LAPACK❗✔️:      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
    // LAPACK❗✔️:       LOGICAL            WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
    // LAPACK❗✔️:      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
    // LAPACK❗✔️:      $                   Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  ================================================================
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
    // LAPACK❗✔️:      $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
    // LAPACK❗✔️:       INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
    // LAPACK❗✔️:      $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2,
    // LAPACK❗✔️:      $                   LWKOPT
    // LAPACK❗✔️:       LOGICAL            BULGE, SORTED
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️:       EXTERNAL           DLAMCH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DCOPY, DGEHRD, DGEMM, DLACPY,
    // LAPACK❗✔️:      $                   DLAHQR,
    // LAPACK❗✔️:      $                   DLANV2, DLARF1F, DLARFG, DLASET, DORMHR, DTREXC
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Estimate optimal workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       JW = MIN( NW, KBOT-KTOP+1 )
    // LAPACK❗✔️:       IF( JW.LE.2 ) THEN
    // LAPACK❗✔️:          LWKOPT = 1
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DGEHRD ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
    // LAPACK❗✔️:          LWK1 = INT( WORK( 1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DORMHR ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DORMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V,
    // LAPACK❗✔️:      $                LDV,
    // LAPACK❗✔️:      $                WORK, -1, INFO )
    // LAPACK❗✔️:          LWK2 = INT( WORK( 1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Optimal workspace ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          LWKOPT = JW + MAX( LWK1, LWK2 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Quick return in case of workspace query. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LWORK.EQ.-1 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Nothing to do ...
    // LAPACK❗✔️: *     ... for an empty active block ... ====
    // LAPACK❗✔️:       NS = 0
    // LAPACK❗✔️:       ND = 0
    // LAPACK❗✔️:       WORK( 1 ) = ONE
    // LAPACK❗✔️:       IF( KTOP.GT.KBOT )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *     ... nor for an empty deflation window. ====
    // LAPACK❗✔️:       IF( NW.LT.1 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Machine constants ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'SAFE MINIMUM' )
    // LAPACK❗✔️:       SAFMAX = ONE / SAFMIN
    // LAPACK❗✔️:       ULP = DLAMCH( 'PRECISION' )
    // LAPACK❗✔️:       SMLNUM = SAFMIN*( DBLE( N ) / ULP )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Setup deflation window ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       JW = MIN( NW, KBOT-KTOP+1 )
    // LAPACK❗✔️:       KWTOP = KBOT - JW + 1
    // LAPACK❗✔️:       IF( KWTOP.EQ.KTOP ) THEN
    // LAPACK❗✔️:          S = ZERO
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          S = H( KWTOP, KWTOP-1 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( KBOT.EQ.KWTOP ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== 1-by-1 deflation window: not much to do ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SR( KWTOP ) = H( KWTOP, KWTOP )
    // LAPACK❗✔️:          SI( KWTOP ) = ZERO
    // LAPACK❗✔️:          NS = 1
    // LAPACK❗✔️:          ND = 0
    // LAPACK❗✔️:          IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) )
    // LAPACK❗✔️:      $        THEN
    // LAPACK❗✔️:             NS = 0
    // LAPACK❗✔️:             ND = 1
    // LAPACK❗✔️:             IF( KWTOP.GT.KTOP )
    // LAPACK❗✔️:      $         H( KWTOP, KWTOP-1 ) = ZERO
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          WORK( 1 ) = ONE
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Convert to spike-triangular form.  (In case of a
    // LAPACK❗✔️: *     .    rare QR failure, this routine continues to do
    // LAPACK❗✔️: *     .    aggressive early deflation using that part of
    // LAPACK❗✔️: *     .    the deflation window that converged using INFQR
    // LAPACK❗✔️: *     .    here and there to keep track.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
    // LAPACK❗✔️:       CALL DCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ),
    // LAPACK❗✔️:      $            LDT+1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
    // LAPACK❗✔️:       CALL DLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ),
    // LAPACK❗✔️:      $             SI( KWTOP ), 1, JW, V, LDV, INFQR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== DTREXC needs a clean margin near the diagonal ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 J = 1, JW - 3
    // LAPACK❗✔️:          T( J+2, J ) = ZERO
    // LAPACK❗✔️:          T( J+3, J ) = ZERO
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       IF( JW.GT.2 )
    // LAPACK❗✔️:      $   T( JW, JW-2 ) = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Deflation detection loop ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NS = JW
    // LAPACK❗✔️:       ILST = INFQR + 1
    // LAPACK❗✔️:    20 CONTINUE
    // LAPACK❗✔️:       IF( ILST.LE.NS ) THEN
    // LAPACK❗✔️:          IF( NS.EQ.1 ) THEN
    // LAPACK❗✔️:             BULGE = .FALSE.
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             BULGE = T( NS, NS-1 ).NE.ZERO
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Small spike tip test for deflation ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( .NOT.BULGE ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Real eigenvalue ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             FOO = ABS( T( NS, NS ) )
    // LAPACK❗✔️:             IF( FOO.EQ.ZERO )
    // LAPACK❗✔️:      $         FOO = ABS( S )
    // LAPACK❗✔️:             IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Deflatable ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = NS - 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Undeflatable.   Move it up out of the way.
    // LAPACK❗✔️: *              .    (DTREXC can not fail in this case.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IFST = NS
    // LAPACK❗✔️:                CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST,
    // LAPACK❗✔️:      $                      WORK,
    // LAPACK❗✔️:      $                      INFO )
    // LAPACK❗✔️:                ILST = ILST + 1
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Complex conjugate pair ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )*
    // LAPACK❗✔️:      $            SQRT( ABS( T( NS-1, NS ) ) )
    // LAPACK❗✔️:             IF( FOO.EQ.ZERO )
    // LAPACK❗✔️:      $         FOO = ABS( S )
    // LAPACK❗✔️:             IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE.
    // LAPACK❗✔️:      $          MAX( SMLNUM, ULP*FOO ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Deflatable ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = NS - 2
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Undeflatable. Move them up out of the way.
    // LAPACK❗✔️: *              .    Fortunately, DTREXC does the right thing with
    // LAPACK❗✔️: *              .    ILST in case of a rare exchange failure. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IFST = NS
    // LAPACK❗✔️:                CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST,
    // LAPACK❗✔️:      $                      WORK,
    // LAPACK❗✔️:      $                      INFO )
    // LAPACK❗✔️:                ILST = ILST + 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== End deflation detection loop ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          GO TO 20
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Return to Hessenberg form ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NS.EQ.0 )
    // LAPACK❗✔️:      $   S = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NS.LT.JW ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== sorting diagonal blocks of T improves accuracy for
    // LAPACK❗✔️: *        .    graded matrices.  Bubble sort deals well with
    // LAPACK❗✔️: *        .    exchange failures. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SORTED = .false.
    // LAPACK❗✔️:          I = NS + 1
    // LAPACK❗✔️:    30    CONTINUE
    // LAPACK❗✔️:          IF( SORTED )
    // LAPACK❗✔️:      $      GO TO 50
    // LAPACK❗✔️:          SORTED = .true.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KEND = I - 1
    // LAPACK❗✔️:          I = INFQR + 1
    // LAPACK❗✔️:          IF( I.EQ.NS ) THEN
    // LAPACK❗✔️:             K = I + 1
    // LAPACK❗✔️:          ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
    // LAPACK❗✔️:             K = I + 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             K = I + 2
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️:          IF( K.LE.KEND ) THEN
    // LAPACK❗✔️:             IF( K.EQ.I+1 ) THEN
    // LAPACK❗✔️:                EVI = ABS( T( I, I ) )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )*
    // LAPACK❗✔️:      $               SQRT( ABS( T( I, I+1 ) ) )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( K.EQ.KEND ) THEN
    // LAPACK❗✔️:                EVK = ABS( T( K, K ) )
    // LAPACK❗✔️:             ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                EVK = ABS( T( K, K ) )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )*
    // LAPACK❗✔️:      $               SQRT( ABS( T( K, K+1 ) ) )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( EVI.GE.EVK ) THEN
    // LAPACK❗✔️:                I = K
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                SORTED = .false.
    // LAPACK❗✔️:                IFST = I
    // LAPACK❗✔️:                ILST = K
    // LAPACK❗✔️:                CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST,
    // LAPACK❗✔️:      $                      WORK,
    // LAPACK❗✔️:      $                      INFO )
    // LAPACK❗✔️:                IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️:                   I = ILST
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   I = K
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( I.EQ.KEND ) THEN
    // LAPACK❗✔️:                K = I + 1
    // LAPACK❗✔️:             ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                K = I + 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                K = I + 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             GO TO 40
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          GO TO 30
    // LAPACK❗✔️:    50    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Restore shift/eigenvalue array from T ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       I = JW
    // LAPACK❗✔️:    60 CONTINUE
    // LAPACK❗✔️:       IF( I.GE.INFQR+1 ) THEN
    // LAPACK❗✔️:          IF( I.EQ.INFQR+1 ) THEN
    // LAPACK❗✔️:             SR( KWTOP+I-1 ) = T( I, I )
    // LAPACK❗✔️:             SI( KWTOP+I-1 ) = ZERO
    // LAPACK❗✔️:             I = I - 1
    // LAPACK❗✔️:          ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
    // LAPACK❗✔️:             SR( KWTOP+I-1 ) = T( I, I )
    // LAPACK❗✔️:             SI( KWTOP+I-1 ) = ZERO
    // LAPACK❗✔️:             I = I - 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             AA = T( I-1, I-1 )
    // LAPACK❗✔️:             CC = T( I, I-1 )
    // LAPACK❗✔️:             BB = T( I-1, I )
    // LAPACK❗✔️:             DD = T( I, I )
    // LAPACK❗✔️:             CALL DLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ),
    // LAPACK❗✔️:      $                   SI( KWTOP+I-2 ), SR( KWTOP+I-1 ),
    // LAPACK❗✔️:      $                   SI( KWTOP+I-1 ), CS, SN )
    // LAPACK❗✔️:             I = I - 2
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          GO TO 60
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
    // LAPACK❗✔️:          IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Reflect spike back into lower triangle ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DCOPY( NS, V, LDV, WORK, 1 )
    // LAPACK❗✔️:             BETA = WORK( 1 )
    // LAPACK❗✔️:             CALL DLARFG( NS, BETA, WORK( 2 ), 1, TAU )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ),
    // LAPACK❗✔️:      $                   LDT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARF1F( 'L', NS, JW, WORK, 1, TAU, T, LDT,
    // LAPACK❗✔️:      $                  WORK( JW+1 ) )
    // LAPACK❗✔️:             CALL DLARF1F( 'R', NS, NS, WORK, 1, TAU, T, LDT,
    // LAPACK❗✔️:      $                  WORK( JW+1 ) )
    // LAPACK❗✔️:             CALL DLARF1F( 'R', JW, NS, WORK, 1, TAU, V, LDV,
    // LAPACK❗✔️:      $                  WORK( JW+1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
    // LAPACK❗✔️:      $                   LWORK-JW, INFO )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Copy updated reduced window into place ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( KWTOP.GT.1 )
    // LAPACK❗✔️:      $      H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
    // LAPACK❗✔️:          CALL DLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
    // LAPACK❗✔️:          CALL DCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),
    // LAPACK❗✔️:      $               LDH+1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Accumulate orthogonal matrix in order update
    // LAPACK❗✔️: *        .    H and Z, if requested.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( NS.GT.1 .AND. S.NE.ZERO )
    // LAPACK❗✔️:      $      CALL DORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V,
    // LAPACK❗✔️:      $                   LDV,
    // LAPACK❗✔️:      $                   WORK( JW+1 ), LWORK-JW, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Update vertical slab in H ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️:             LTOP = 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             LTOP = KTOP
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          DO 70 KROW = LTOP, KWTOP - 1, NV
    // LAPACK❗✔️:             KLN = MIN( NV, KWTOP-KROW )
    // LAPACK❗✔️:             CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),
    // LAPACK❗✔️:      $                  LDH, V, LDV, ZERO, WV, LDWV )
    // LAPACK❗✔️:             CALL DLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ),
    // LAPACK❗✔️:      $                   LDH )
    // LAPACK❗✔️:    70    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Update horizontal slab in H ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️:             DO 80 KCOL = KBOT + 1, N, NH
    // LAPACK❗✔️:                KLN = MIN( NH, N-KCOL+1 )
    // LAPACK❗✔️:                CALL DGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,
    // LAPACK❗✔️:      $                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
    // LAPACK❗✔️:                CALL DLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),
    // LAPACK❗✔️:      $                      LDH )
    // LAPACK❗✔️:    80       CONTINUE
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Update vertical slab in Z ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTZ ) THEN
    // LAPACK❗✔️:             DO 90 KROW = ILOZ, IHIZ, NV
    // LAPACK❗✔️:                KLN = MIN( NV, IHIZ-KROW+1 )
    // LAPACK❗✔️:                CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW,
    // LAPACK❗✔️:      $                     KWTOP ),
    // LAPACK❗✔️:      $                     LDZ, V, LDV, ZERO, WV, LDWV )
    // LAPACK❗✔️:                CALL DLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),
    // LAPACK❗✔️:      $                      LDZ )
    // LAPACK❗✔️:    90       CONTINUE
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Return the number of deflations ... ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ND = JW - NS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== ... and the number of shifts. (Subtracting
    // LAPACK❗✔️: *     .    INFQR from the spike length takes care
    // LAPACK❗✔️: *     .    of the case of a rare QR failure while
    // LAPACK❗✔️: *     .    calculating eigenvalues of the deflation
    // LAPACK❗✔️: *     .    window.)  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NS = NS - INFQR
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *      ==== Return optimal workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== End of DLAQR2 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if ktop > kbot || nw < 1 {
        return Ok((0, 0));
    }
    let safmin = dlamch(b'S');
    let ulp = dlamch(b'P');
    let smlnum = safmin * (n as f64 / ulp);
    let jw = nw.min(kbot - ktop + 1);
    let kwtop = kbot - jw + 1;
    let mut s = if kwtop == ktop {
        0.0
    } else {
        h[ix(kwtop, kwtop - 1, ldh)]
    };
    macro_rules! t {
        ($i:expr,$j:expr) => {
            h[toff + ix($i, $j, ldt)]
        };
    }
    macro_rules! v {
        ($i:expr,$j:expr) => {
            h[voff + ix($i, $j, ldv)]
        };
    }
    if kbot == kwtop {
        sr[kwtop - 1] = h[ix(kwtop, kwtop, ldh)];
        si[kwtop - 1] = 0.0;
        if s.abs() <= smlnum.max(ulp * h[ix(kwtop, kwtop, ldh)].abs()) {
            if kwtop > ktop {
                h[ix(kwtop, kwtop - 1, ldh)] = 0.0;
            }
            return Ok((0, 1));
        }
        return Ok((1, 0));
    }
    for j in 1..=jw {
        for i in 1..=j {
            t!(i, j) = h[ix(kwtop + i - 1, kwtop + j - 1, ldh)];
        }
    }
    for i in 1..jw {
        t!(i + 1, i) = h[ix(kwtop + i, kwtop + i - 1, ldh)];
    }
    for j in 1..=jw {
        for i in 1..=jw {
            v!(i, j) = if i == j { 1.0 } else { 0.0 };
        }
    }
    // V lies entirely in earlier columns than T; split at T's first element.
    let infqr = {
        let (before, t) = h.split_at_mut(toff);
        let v = &mut before[voff..];
        if recursive && jw > 75 {
            qr(
                false,
                true,
                true,
                jw,
                1,
                jw,
                t,
                ldt,
                &mut sr[kwtop - 1..],
                &mut si[kwtop - 1..],
                1,
                jw,
                v,
                ldv,
                work,
            )?
        } else {
            super::schur::lahqr(
                true,
                true,
                jw,
                1,
                jw,
                t,
                ldt,
                &mut sr[kwtop - 1..],
                &mut si[kwtop - 1..],
                1,
                jw,
                v,
                ldv,
            )
        }
    };
    if jw >= 4 {
        for j in 1..=jw - 3 {
            t!(j + 2, j) = 0.0;
            t!(j + 3, j) = 0.0;
        }
    }
    if jw > 2 {
        t!(jw, jw - 2) = 0.0;
    }
    let mut ns = jw;
    let mut ilst = infqr + 1;
    while ilst <= ns {
        let bulge = ns > 1 && t!(ns, ns - 1) != 0.0;
        if !bulge {
            let mut foo = t!(ns, ns).abs();
            if foo == 0.0 {
                foo = s.abs();
            }
            if (s * v!(1, ns)).abs() <= smlnum.max(ulp * foo) {
                ns -= 1;
            } else {
                let mut ifst = ns;
                let (before, t) = h.split_at_mut(toff);
                super::swaps::trexc(jw, t, ldt, &mut before[voff..], ldv, &mut ifst, &mut ilst);
                ilst += 1;
            }
        } else {
            let mut foo =
                t!(ns, ns).abs() + t!(ns, ns - 1).abs().sqrt() * t!(ns - 1, ns).abs().sqrt();
            if foo == 0.0 {
                foo = s.abs();
            }
            if (s * v!(1, ns)).abs().max((s * v!(1, ns - 1)).abs()) <= smlnum.max(ulp * foo) {
                ns -= 2;
            } else {
                let mut ifst = ns;
                let (before, t) = h.split_at_mut(toff);
                super::swaps::trexc(jw, t, ldt, &mut before[voff..], ldv, &mut ifst, &mut ilst);
                ilst += 2;
            }
        }
    }
    if ns == 0 {
        s = 0.0;
    }
    if ns < jw {
        let mut sorted = false;
        let mut i = ns + 1;
        while !sorted {
            sorted = true;
            let kend = i - 1;
            i = infqr + 1;
            let mut k = if i == ns || t!(i + 1, i) == 0.0 {
                i + 1
            } else {
                i + 2
            };
            while k <= kend {
                let evi = if k == i + 1 {
                    t!(i, i).abs()
                } else {
                    t!(i, i).abs() + t!(i + 1, i).abs().sqrt() * t!(i, i + 1).abs().sqrt()
                };
                let evk = if k == kend || t!(k + 1, k) == 0.0 {
                    t!(k, k).abs()
                } else {
                    t!(k, k).abs() + t!(k + 1, k).abs().sqrt() * t!(k, k + 1).abs().sqrt()
                };
                if evi >= evk {
                    i = k;
                } else {
                    sorted = false;
                    let mut ifst = i;
                    let mut ilst = k;
                    let info = {
                        let (before, t) = h.split_at_mut(toff);
                        super::swaps::trexc(
                            jw,
                            t,
                            ldt,
                            &mut before[voff..],
                            ldv,
                            &mut ifst,
                            &mut ilst,
                        )
                    };
                    i = if info == 0 { ilst } else { k };
                }
                k = if i == kend || t!(i + 1, i) == 0.0 {
                    i + 1
                } else {
                    i + 2
                };
            }
        }
    }
    let mut i = jw;
    while i >= infqr + 1 {
        if i == infqr + 1 || t!(i, i - 1) == 0.0 {
            sr[kwtop + i - 2] = t!(i, i);
            si[kwtop + i - 2] = 0.0;
            i -= 1;
        } else {
            let pair = super::lanv2(t!(i - 1, i - 1), t!(i - 1, i), t!(i, i - 1), t!(i, i));
            sr[kwtop + i - 3] = pair.real[0];
            si[kwtop + i - 3] = pair.imag[0];
            sr[kwtop + i - 2] = pair.real[1];
            si[kwtop + i - 2] = pair.imag[1];
            i -= 2;
        }
    }
    if ns < jw || s == 0.0 {
        if ns > 1 && s != 0.0 {
            for i in 1..=ns {
                work[i - 1] = v!(1, i);
            }
            let beta = work[0];
            let tau = larfg(&mut work[..ns]);
            work[0] = beta;
            for j in 1..=jw - 2 {
                for i in j..=jw - 2 {
                    t!(i + 2, j) = 0.0;
                }
            }
            {
                let (reflector, scratch) = work.split_at_mut(jw);
                super::hessenberg::larf1f(
                    true,
                    ns,
                    jw,
                    reflector,
                    tau,
                    &mut h[toff..],
                    ldt,
                    scratch,
                );
                super::hessenberg::larf1f(
                    false,
                    ns,
                    ns,
                    reflector,
                    tau,
                    &mut h[toff..],
                    ldt,
                    scratch,
                );
                super::hessenberg::larf1f(
                    false,
                    jw,
                    ns,
                    reflector,
                    tau,
                    &mut h[voff..],
                    ldv,
                    scratch,
                );
                super::hessenberg::gehrd(jw, 1, ns, &mut h[toff..], ldt, reflector, scratch)?;
            }
        }
        if kwtop > 1 {
            h[ix(kwtop, kwtop - 1, ldh)] = s * v!(1, 1);
        }
        for j in 1..=jw {
            for i in 1..=j {
                h[ix(kwtop + i - 1, kwtop + j - 1, ldh)] = t!(i, j);
            }
        }
        for i in 1..jw {
            h[ix(kwtop + i, kwtop + i - 1, ldh)] = t!(i + 1, i);
        }
        if ns > 1 && s != 0.0 {
            let (tau, scratch) = work.split_at_mut(jw);
            let (before, t) = h.split_at_mut(toff);
            super::hessenberg::ormhr(
                jw,
                ns,
                1,
                ns,
                t,
                ldt,
                tau,
                &mut before[voff..],
                ldv,
                scratch,
            )?;
        }
        let ltop = if wantt { 1 } else { ktop };
        for krow in (ltop..kwtop).step_by(nv) {
            let kln = nv.min(kwtop - krow);
            blas3::gemm_within(
                NoTrans,
                NoTrans,
                kln,
                jw,
                jw,
                1.0,
                h,
                ix(krow, kwtop, ldh),
                ldh,
                voff,
                ldv,
                0.0,
                wvoff,
                ldwv,
            )?;
            for j in 1..=jw {
                for i in 1..=kln {
                    h[ix(krow + i - 1, kwtop + j - 1, ldh)] = h[wvoff + ix(i, j, ldwv)];
                }
            }
        }
        if wantt {
            for kcol in (kbot + 1..=n).step_by(nh) {
                let kln = nh.min(n - kcol + 1);
                blas3::gemm_within(
                    Trans,
                    NoTrans,
                    jw,
                    kln,
                    jw,
                    1.0,
                    h,
                    voff,
                    ldv,
                    ix(kwtop, kcol, ldh),
                    ldh,
                    0.0,
                    toff,
                    ldt,
                )?;
                for j in 1..=kln {
                    for i in 1..=jw {
                        h[ix(kwtop + i - 1, kcol + j - 1, ldh)] = t!(i, j);
                    }
                }
            }
        }
        if wantz {
            for krow in (iloz..=ihiz).step_by(nv) {
                let kln = nv.min(ihiz - krow + 1);
                blas3::gemm_into_within(
                    NoTrans,
                    NoTrans,
                    kln,
                    jw,
                    jw,
                    1.0,
                    h,
                    Separate(&z[ix(krow, kwtop, ldz)..]),
                    ldz,
                    Within(voff),
                    ldv,
                    0.0,
                    wvoff,
                    ldwv,
                )?;
                for j in 1..=jw {
                    for i in 1..=kln {
                        z[ix(krow + i - 1, kwtop + j - 1, ldz)] = h[wvoff + ix(i, j, ldwv)];
                    }
                }
            }
        }
    }
    Ok((ns - infqr, jw - ns))
}

/// Shared literal DLAQR0/DLAQR4 loop. The former uses DLAQR3 and may call
/// DLAQR4 for shift generation; the latter uses DLAQR2 and DLAHQR there.
#[allow(clippy::too_many_arguments)]
pub(super) fn qr(
    recursive: bool,
    wantt: bool,
    wantz: bool,
    n: usize,
    ilo: usize,
    ihi: usize,
    h: &mut [f64],
    ldh: usize,
    wr: &mut [f64],
    wi: &mut [f64],
    iloz: usize,
    ihiz: usize,
    z: &mut [f64],
    ldz: usize,
    work: &mut [f64],
) -> Result<usize, GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
    // LAPACK❗✔️:      $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
    // LAPACK❗✔️:       LOGICAL            WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
    // LAPACK❗✔️:      $                   Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  ================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Matrices of order NTINY or smaller must be processed by
    // LAPACK❗✔️: *     .    DLAHQR because of insufficient subdiagonal scratch space.
    // LAPACK❗✔️: *     .    (This is a hard limit.) ====
    // LAPACK❗✔️:       INTEGER            NTINY
    // LAPACK❗✔️:       PARAMETER          ( NTINY = 15 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Exceptional deflation windows:  try to cure rare
    // LAPACK❗✔️: *     .    slow convergence by varying the size of the
    // LAPACK❗✔️: *     .    deflation window after KEXNW iterations. ====
    // LAPACK❗✔️:       INTEGER            KEXNW
    // LAPACK❗✔️:       PARAMETER          ( KEXNW = 5 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Exceptional shifts: try to cure rare slow convergence
    // LAPACK❗✔️: *     .    with ad-hoc exceptional shifts every KEXSH iterations.
    // LAPACK❗✔️: *     .    ====
    // LAPACK❗✔️:       INTEGER            KEXSH
    // LAPACK❗✔️:       PARAMETER          ( KEXSH = 6 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== The constants WILK1 and WILK2 are used to form the
    // LAPACK❗✔️: *     .    exceptional shifts. ====
    // LAPACK❗✔️:       DOUBLE PRECISION   WILK1, WILK2
    // LAPACK❗✔️:       PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
    // LAPACK❗✔️:       INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
    // LAPACK❗✔️:      $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
    // LAPACK❗✔️:      $                   LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS,
    // LAPACK❗✔️:      $                   NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
    // LAPACK❗✔️:       LOGICAL            SORTED
    // LAPACK❗✔️:       CHARACTER          JBCMPZ*2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZDUM( 1, 1 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR3, DLAQR4, DLAQR5
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Quick return for N = 0: nothing to do. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = ONE
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.NTINY ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Tiny matrices must use DLAHQR. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          LWKOPT = 1
    // LAPACK❗✔️:          IF( LWORK.NE.-1 )
    // LAPACK❗✔️:      $      CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
    // LAPACK❗✔️:      $                   ILOZ, IHIZ, Z, LDZ, INFO )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Use small bulge multi-shift QR with aggressive early
    // LAPACK❗✔️: *        .    deflation on larger-than-tiny matrices. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Hope for the best. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Set up job flags for ILAENV. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️:             JBCMPZ( 1: 1 ) = 'S'
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             JBCMPZ( 1: 1 ) = 'E'
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( WANTZ ) THEN
    // LAPACK❗✔️:             JBCMPZ( 2: 2 ) = 'V'
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             JBCMPZ( 2: 2 ) = 'N'
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NWR = recommended deflation window size.  At this
    // LAPACK❗✔️: *        .    point,  N .GT. NTINY = 15, so there is enough
    // LAPACK❗✔️: *        .    subdiagonal workspace for NWR.GE.2 as required.
    // LAPACK❗✔️: *        .    (In fact, there is enough subdiagonal space for
    // LAPACK❗✔️: *        .    NWR.GE.4.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NWR = ILAENV( 13, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NWR = MAX( 2, NWR )
    // LAPACK❗✔️:          NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NSR = recommended number of simultaneous shifts.
    // LAPACK❗✔️: *        .    At this point N .GT. NTINY = 15, so there is at
    // LAPACK❗✔️: *        .    enough subdiagonal workspace for NSR to be even
    // LAPACK❗✔️: *        .    and greater than or equal to two as required. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NSR = ILAENV( 15, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
    // LAPACK❗✔️:          NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Estimate optimal workspace ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DLAQR3 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,
    // LAPACK❗✔️:      $                IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH,
    // LAPACK❗✔️:      $                N, H, LDH, WORK, -1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Quick return in case of workspace query. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LWORK.EQ.-1 ) THEN
    // LAPACK❗✔️:             WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️:             RETURN
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== DLAHQR/DLAQR0 crossover point ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NMIN = ILAENV( 12, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NMIN = MAX( NTINY, NMIN )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Nibble crossover point ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NIBBLE = ILAENV( 14, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NIBBLE = MAX( 0, NIBBLE )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Accumulate reflections during ttswp?  Use block
    // LAPACK❗✔️: *        .    2-by-2 structure during matrix-matrix multiply? ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KACC22 = ILAENV( 16, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          KACC22 = MAX( 0, KACC22 )
    // LAPACK❗✔️:          KACC22 = MIN( 2, KACC22 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NWMAX = the largest possible deflation window for
    // LAPACK❗✔️: *        .    which there is sufficient workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
    // LAPACK❗✔️:          NW = NWMAX
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NSMAX = the Largest number of simultaneous shifts
    // LAPACK❗✔️: *        .    for which there is sufficient workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NSMAX = MIN( ( N-3 ) / 6, 2*LWORK / 3 )
    // LAPACK❗✔️:          NSMAX = NSMAX - MOD( NSMAX, 2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NDFL: an iteration count restarted at deflation. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NDFL = 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== ITMAX = iteration limit ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Last row and column in the active block ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KBOT = IHI
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Main Loop ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 80 IT = 1, ITMAX
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Done when KBOT falls below ILO ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( KBOT.LT.ILO )
    // LAPACK❗✔️:      $         GO TO 90
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Locate active block ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 10 K = KBOT, ILO + 1, -1
    // LAPACK❗✔️:                IF( H( K, K-1 ).EQ.ZERO )
    // LAPACK❗✔️:      $            GO TO 20
    // LAPACK❗✔️:    10       CONTINUE
    // LAPACK❗✔️:             K = ILO
    // LAPACK❗✔️:    20       CONTINUE
    // LAPACK❗✔️:             KTOP = K
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Select deflation window size:
    // LAPACK❗✔️: *           .    Typical Case:
    // LAPACK❗✔️: *           .      If possible and advisable, nibble the entire
    // LAPACK❗✔️: *           .      active block.  If not, use size MIN(NWR,NWMAX)
    // LAPACK❗✔️: *           .      or MIN(NWR+1,NWMAX) depending upon which has
    // LAPACK❗✔️: *           .      the smaller corresponding subdiagonal entry
    // LAPACK❗✔️: *           .      (a heuristic).
    // LAPACK❗✔️: *           .
    // LAPACK❗✔️: *           .    Exceptional Case:
    // LAPACK❗✔️: *           .      If there have been no deflations in KEXNW or
    // LAPACK❗✔️: *           .      more iterations, then vary the deflation window
    // LAPACK❗✔️: *           .      size.   At first, because, larger windows are,
    // LAPACK❗✔️: *           .      in general, more powerful than smaller ones,
    // LAPACK❗✔️: *           .      rapidly increase the window to the maximum possible.
    // LAPACK❗✔️: *           .      Then, gradually reduce the window size. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NH = KBOT - KTOP + 1
    // LAPACK❗✔️:             NWUPBD = MIN( NH, NWMAX )
    // LAPACK❗✔️:             IF( NDFL.LT.KEXNW ) THEN
    // LAPACK❗✔️:                NW = MIN( NWUPBD, NWR )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                NW = MIN( NWUPBD, 2*NW )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( NW.LT.NWMAX ) THEN
    // LAPACK❗✔️:                IF( NW.GE.NH-1 ) THEN
    // LAPACK❗✔️:                   NW = NH
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   KWTOP = KBOT - NW + 1
    // LAPACK❗✔️:                   IF( ABS( H( KWTOP, KWTOP-1 ) ).GT.
    // LAPACK❗✔️:      $                ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( NDFL.LT.KEXNW ) THEN
    // LAPACK❗✔️:                NDEC = -1
    // LAPACK❗✔️:             ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
    // LAPACK❗✔️:                NDEC = NDEC + 1
    // LAPACK❗✔️:                IF( NW-NDEC.LT.2 )
    // LAPACK❗✔️:      $            NDEC = 0
    // LAPACK❗✔️:                NW = NW - NDEC
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Aggressive early deflation:
    // LAPACK❗✔️: *           .    split workspace under the subdiagonal into
    // LAPACK❗✔️: *           .      - an nw-by-nw work array V in the lower
    // LAPACK❗✔️: *           .        left-hand-corner,
    // LAPACK❗✔️: *           .      - an NW-by-at-least-NW-but-more-is-better
    // LAPACK❗✔️: *           .        (NW-by-NHO) horizontal work array along
    // LAPACK❗✔️: *           .        the bottom edge,
    // LAPACK❗✔️: *           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
    // LAPACK❗✔️: *           .        vertical work array along the left-hand-edge.
    // LAPACK❗✔️: *           .        ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             KV = N - NW + 1
    // LAPACK❗✔️:             KT = NW + 1
    // LAPACK❗✔️:             NHO = ( N-NW-1 ) - KT + 1
    // LAPACK❗✔️:             KWV = NW + 2
    // LAPACK❗✔️:             NVE = ( N-NW ) - KWV + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Aggressive early deflation ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
    // LAPACK❗✔️:      $                   IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH,
    // LAPACK❗✔️:      $                   NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH,
    // LAPACK❗✔️:      $                   WORK, LWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Adjust KBOT accounting for new deflations. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             KBOT = KBOT - LD
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== KS points to the shifts. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             KS = KBOT - LS + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Skip an expensive QR sweep if there is a (partly
    // LAPACK❗✔️: *           .    heuristic) reason to expect that many eigenvalues
    // LAPACK❗✔️: *           .    will deflate without it.  Here, the QR sweep is
    // LAPACK❗✔️: *           .    skipped if many eigenvalues have just been deflated
    // LAPACK❗✔️: *           .    or if the remaining active block is small.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
    // LAPACK❗✔️:      $          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== NS = nominal number of simultaneous shifts.
    // LAPACK❗✔️: *              .    This may be lowered (slightly) if DLAQR3
    // LAPACK❗✔️: *              .    did not provide that many shifts. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
    // LAPACK❗✔️:                NS = NS - MOD( NS, 2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== If there have been no deflations
    // LAPACK❗✔️: *              .    in a multiple of KEXSH iterations,
    // LAPACK❗✔️: *              .    then try exceptional shifts.
    // LAPACK❗✔️: *              .    Otherwise use shifts provided by
    // LAPACK❗✔️: *              .    DLAQR3 above or from the eigenvalues
    // LAPACK❗✔️: *              .    of a trailing principal submatrix. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
    // LAPACK❗✔️:                   KS = KBOT - NS + 1
    // LAPACK❗✔️:                   DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
    // LAPACK❗✔️:                      SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
    // LAPACK❗✔️:                      AA = WILK1*SS + H( I, I )
    // LAPACK❗✔️:                      BB = SS
    // LAPACK❗✔️:                      CC = WILK2*SS
    // LAPACK❗✔️:                      DD = AA
    // LAPACK❗✔️:                      CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ),
    // LAPACK❗✔️:      $                            WR( I ), WI( I ), CS, SN )
    // LAPACK❗✔️:    30             CONTINUE
    // LAPACK❗✔️:                   IF( KS.EQ.KTOP ) THEN
    // LAPACK❗✔️:                      WR( KS+1 ) = H( KS+1, KS+1 )
    // LAPACK❗✔️:                      WI( KS+1 ) = ZERO
    // LAPACK❗✔️:                      WR( KS ) = WR( KS+1 )
    // LAPACK❗✔️:                      WI( KS ) = WI( KS+1 )
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
    // LAPACK❗✔️: *                 .    DLAHQR on a trailing principal submatrix to
    // LAPACK❗✔️: *                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
    // LAPACK❗✔️: *                 .    there is enough space below the subdiagonal
    // LAPACK❗✔️: *                 .    to fit an NS-by-NS scratch array.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   IF( KBOT-KS+1.LE.NS / 2 ) THEN
    // LAPACK❗✔️:                      KS = KBOT - NS + 1
    // LAPACK❗✔️:                      KT = N - NS + 1
    // LAPACK❗✔️:                      CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH,
    // LAPACK❗✔️:      $                            H( KT, 1 ), LDH )
    // LAPACK❗✔️:                      IF( NS.GT.NMIN ) THEN
    // LAPACK❗✔️:                         CALL DLAQR4( .false., .false., NS, 1, NS,
    // LAPACK❗✔️:      $                               H( KT, 1 ), LDH, WR( KS ),
    // LAPACK❗✔️:      $                               WI( KS ), 1, 1, ZDUM, 1, WORK,
    // LAPACK❗✔️:      $                               LWORK, INF )
    // LAPACK❗✔️:                      ELSE
    // LAPACK❗✔️:                         CALL DLAHQR( .false., .false., NS, 1, NS,
    // LAPACK❗✔️:      $                               H( KT, 1 ), LDH, WR( KS ),
    // LAPACK❗✔️:      $                               WI( KS ), 1, 1, ZDUM, 1, INF )
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                      KS = KS + INF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    ==== In case of a rare QR failure use
    // LAPACK❗✔️: *                    .    eigenvalues of the trailing 2-by-2
    // LAPACK❗✔️: *                    .    principal submatrix.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      IF( KS.GE.KBOT ) THEN
    // LAPACK❗✔️:                         AA = H( KBOT-1, KBOT-1 )
    // LAPACK❗✔️:                         CC = H( KBOT, KBOT-1 )
    // LAPACK❗✔️:                         BB = H( KBOT-1, KBOT )
    // LAPACK❗✔️:                         DD = H( KBOT, KBOT )
    // LAPACK❗✔️:                         CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ),
    // LAPACK❗✔️:      $                               WI( KBOT-1 ), WR( KBOT ),
    // LAPACK❗✔️:      $                               WI( KBOT ), CS, SN )
    // LAPACK❗✔️:                         KS = KBOT - 1
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   IF( KBOT-KS+1.GT.NS ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    ==== Sort the shifts (Helps a little)
    // LAPACK❗✔️: *                    .    Bubble sort keeps complex conjugate
    // LAPACK❗✔️: *                    .    pairs together. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      SORTED = .false.
    // LAPACK❗✔️:                      DO 50 K = KBOT, KS + 1, -1
    // LAPACK❗✔️:                         IF( SORTED )
    // LAPACK❗✔️:      $                     GO TO 60
    // LAPACK❗✔️:                         SORTED = .true.
    // LAPACK❗✔️:                         DO 40 I = KS, K - 1
    // LAPACK❗✔️:                            IF( ABS( WR( I ) )+ABS( WI( I ) ).LT.
    // LAPACK❗✔️:      $                         ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
    // LAPACK❗✔️:                               SORTED = .false.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                               SWAP = WR( I )
    // LAPACK❗✔️:                               WR( I ) = WR( I+1 )
    // LAPACK❗✔️:                               WR( I+1 ) = SWAP
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                               SWAP = WI( I )
    // LAPACK❗✔️:                               WI( I ) = WI( I+1 )
    // LAPACK❗✔️:                               WI( I+1 ) = SWAP
    // LAPACK❗✔️:                            END IF
    // LAPACK❗✔️:    40                   CONTINUE
    // LAPACK❗✔️:    50                CONTINUE
    // LAPACK❗✔️:    60                CONTINUE
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Shuffle shifts into pairs of real shifts
    // LAPACK❗✔️: *                 .    and pairs of complex conjugate shifts
    // LAPACK❗✔️: *                 .    assuming complex conjugate shifts are
    // LAPACK❗✔️: *                 .    already adjacent to one another. (Yes,
    // LAPACK❗✔️: *                 .    they are.)  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   DO 70 I = KBOT, KS + 2, -2
    // LAPACK❗✔️:                      IF( WI( I ).NE.-WI( I-1 ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         SWAP = WR( I )
    // LAPACK❗✔️:                         WR( I ) = WR( I-1 )
    // LAPACK❗✔️:                         WR( I-1 ) = WR( I-2 )
    // LAPACK❗✔️:                         WR( I-2 ) = SWAP
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         SWAP = WI( I )
    // LAPACK❗✔️:                         WI( I ) = WI( I-1 )
    // LAPACK❗✔️:                         WI( I-1 ) = WI( I-2 )
    // LAPACK❗✔️:                         WI( I-2 ) = SWAP
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:    70             CONTINUE
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== If there are only two shifts and both are
    // LAPACK❗✔️: *              .    real, then use only one.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( KBOT-KS+1.EQ.2 ) THEN
    // LAPACK❗✔️:                   IF( WI( KBOT ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                      IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT.
    // LAPACK❗✔️:      $                   ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
    // LAPACK❗✔️:                         WR( KBOT-1 ) = WR( KBOT )
    // LAPACK❗✔️:                      ELSE
    // LAPACK❗✔️:                         WR( KBOT ) = WR( KBOT-1 )
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Use up to NS of the the smallest magnitude
    // LAPACK❗✔️: *              .    shifts.  If there aren't NS shifts available,
    // LAPACK❗✔️: *              .    then use them all, possibly dropping one to
    // LAPACK❗✔️: *              .    make the number of shifts even. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = MIN( NS, KBOT-KS+1 )
    // LAPACK❗✔️:                NS = NS - MOD( NS, 2 )
    // LAPACK❗✔️:                KS = KBOT - NS + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Small-bulge multi-shift QR sweep:
    // LAPACK❗✔️: *              .    split workspace under the subdiagonal into
    // LAPACK❗✔️: *              .    - a KDU-by-KDU work array U in the lower
    // LAPACK❗✔️: *              .      left-hand-corner,
    // LAPACK❗✔️: *              .    - a KDU-by-at-least-KDU-but-more-is-better
    // LAPACK❗✔️: *              .      (KDU-by-NHo) horizontal work array WH along
    // LAPACK❗✔️: *              .      the bottom edge,
    // LAPACK❗✔️: *              .    - and an at-least-KDU-but-more-is-better-by-KDU
    // LAPACK❗✔️: *              .      (NVE-by-KDU) vertical work WV arrow along
    // LAPACK❗✔️: *              .      the left-hand-edge. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                KDU = 2*NS
    // LAPACK❗✔️:                KU = N - KDU + 1
    // LAPACK❗✔️:                KWH = KDU + 1
    // LAPACK❗✔️:                NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
    // LAPACK❗✔️:                KWV = KDU + 4
    // LAPACK❗✔️:                NVE = N - KDU - KWV + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Small-bulge multi-shift QR sweep ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,
    // LAPACK❗✔️:      $                      WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z,
    // LAPACK❗✔️:      $                      LDZ, WORK, 3, H( KU, 1 ), LDH, NVE,
    // LAPACK❗✔️:      $                      H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Note progress (or the lack of it). ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LD.GT.0 ) THEN
    // LAPACK❗✔️:                NDFL = 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                NDFL = NDFL + 1
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== End of main loop ====
    // LAPACK❗✔️:    80    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Iteration limit exceeded.  Set INFO to show where
    // LAPACK❗✔️: *        .    the problem occurred and exit. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          INFO = KBOT
    // LAPACK❗✔️:    90    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Return the optimal value of LWORK. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== End of DLAQR0 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    // LAPACK❗✔️:       SUBROUTINE DLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
    // LAPACK❗✔️:      $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
    // LAPACK❗✔️:       LOGICAL            WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
    // LAPACK❗✔️:      $                   Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  ================================================================
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Matrices of order NTINY or smaller must be processed by
    // LAPACK❗✔️: *     .    DLAHQR because of insufficient subdiagonal scratch space.
    // LAPACK❗✔️: *     .    (This is a hard limit.) ====
    // LAPACK❗✔️:       INTEGER            NTINY
    // LAPACK❗✔️:       PARAMETER          ( NTINY = 15 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Exceptional deflation windows:  try to cure rare
    // LAPACK❗✔️: *     .    slow convergence by varying the size of the
    // LAPACK❗✔️: *     .    deflation window after KEXNW iterations. ====
    // LAPACK❗✔️:       INTEGER            KEXNW
    // LAPACK❗✔️:       PARAMETER          ( KEXNW = 5 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Exceptional shifts: try to cure rare slow convergence
    // LAPACK❗✔️: *     .    with ad-hoc exceptional shifts every KEXSH iterations.
    // LAPACK❗✔️: *     .    ====
    // LAPACK❗✔️:       INTEGER            KEXSH
    // LAPACK❗✔️:       PARAMETER          ( KEXSH = 6 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== The constants WILK1 and WILK2 are used to form the
    // LAPACK❗✔️: *     .    exceptional shifts. ====
    // LAPACK❗✔️:       DOUBLE PRECISION   WILK1, WILK2
    // LAPACK❗✔️:       PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
    // LAPACK❗✔️:       INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
    // LAPACK❗✔️:      $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
    // LAPACK❗✔️:      $                   LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS,
    // LAPACK❗✔️:      $                   NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
    // LAPACK❗✔️:       LOGICAL            SORTED
    // LAPACK❗✔️:       CHARACTER          JBCMPZ*2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZDUM( 1, 1 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR2, DLAQR5
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Quick return for N = 0: nothing to do. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = ONE
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.NTINY ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Tiny matrices must use DLAHQR. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          LWKOPT = 1
    // LAPACK❗✔️:          IF( LWORK.NE.-1 )
    // LAPACK❗✔️:      $      CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
    // LAPACK❗✔️:      $                   ILOZ, IHIZ, Z, LDZ, INFO )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Use small bulge multi-shift QR with aggressive early
    // LAPACK❗✔️: *        .    deflation on larger-than-tiny matrices. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Hope for the best. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Set up job flags for ILAENV. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️:             JBCMPZ( 1: 1 ) = 'S'
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             JBCMPZ( 1: 1 ) = 'E'
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( WANTZ ) THEN
    // LAPACK❗✔️:             JBCMPZ( 2: 2 ) = 'V'
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             JBCMPZ( 2: 2 ) = 'N'
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NWR = recommended deflation window size.  At this
    // LAPACK❗✔️: *        .    point,  N .GT. NTINY = 15, so there is enough
    // LAPACK❗✔️: *        .    subdiagonal workspace for NWR.GE.2 as required.
    // LAPACK❗✔️: *        .    (In fact, there is enough subdiagonal space for
    // LAPACK❗✔️: *        .    NWR.GE.4.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NWR = ILAENV( 13, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NWR = MAX( 2, NWR )
    // LAPACK❗✔️:          NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NSR = recommended number of simultaneous shifts.
    // LAPACK❗✔️: *        .    At this point N .GT. NTINY = 15, so there is at
    // LAPACK❗✔️: *        .    enough subdiagonal workspace for NSR to be even
    // LAPACK❗✔️: *        .    and greater than or equal to two as required. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NSR = ILAENV( 15, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
    // LAPACK❗✔️:          NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Estimate optimal workspace ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Workspace query call to DLAQR2 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLAQR2( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,
    // LAPACK❗✔️:      $                IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH,
    // LAPACK❗✔️:      $                N, H, LDH, WORK, -1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Optimal workspace = MAX(DLAQR5, DLAQR2) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Quick return in case of workspace query. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LWORK.EQ.-1 ) THEN
    // LAPACK❗✔️:             WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️:             RETURN
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== DLAHQR/DLAQR0 crossover point ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NMIN = ILAENV( 12, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NMIN = MAX( NTINY, NMIN )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Nibble crossover point ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NIBBLE = ILAENV( 14, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          NIBBLE = MAX( 0, NIBBLE )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Accumulate reflections during ttswp?  Use block
    // LAPACK❗✔️: *        .    2-by-2 structure during matrix-matrix multiply? ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KACC22 = ILAENV( 16, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
    // LAPACK❗✔️:          KACC22 = MAX( 0, KACC22 )
    // LAPACK❗✔️:          KACC22 = MIN( 2, KACC22 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NWMAX = the largest possible deflation window for
    // LAPACK❗✔️: *        .    which there is sufficient workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
    // LAPACK❗✔️:          NW = NWMAX
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NSMAX = the Largest number of simultaneous shifts
    // LAPACK❗✔️: *        .    for which there is sufficient workspace. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NSMAX = MIN( ( N-3 ) / 6, 2*LWORK / 3 )
    // LAPACK❗✔️:          NSMAX = NSMAX - MOD( NSMAX, 2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== NDFL: an iteration count restarted at deflation. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NDFL = 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== ITMAX = iteration limit ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Last row and column in the active block ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KBOT = IHI
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Main Loop ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 80 IT = 1, ITMAX
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Done when KBOT falls below ILO ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( KBOT.LT.ILO )
    // LAPACK❗✔️:      $         GO TO 90
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Locate active block ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 10 K = KBOT, ILO + 1, -1
    // LAPACK❗✔️:                IF( H( K, K-1 ).EQ.ZERO )
    // LAPACK❗✔️:      $            GO TO 20
    // LAPACK❗✔️:    10       CONTINUE
    // LAPACK❗✔️:             K = ILO
    // LAPACK❗✔️:    20       CONTINUE
    // LAPACK❗✔️:             KTOP = K
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Select deflation window size:
    // LAPACK❗✔️: *           .    Typical Case:
    // LAPACK❗✔️: *           .      If possible and advisable, nibble the entire
    // LAPACK❗✔️: *           .      active block.  If not, use size MIN(NWR,NWMAX)
    // LAPACK❗✔️: *           .      or MIN(NWR+1,NWMAX) depending upon which has
    // LAPACK❗✔️: *           .      the smaller corresponding subdiagonal entry
    // LAPACK❗✔️: *           .      (a heuristic).
    // LAPACK❗✔️: *           .
    // LAPACK❗✔️: *           .    Exceptional Case:
    // LAPACK❗✔️: *           .      If there have been no deflations in KEXNW or
    // LAPACK❗✔️: *           .      more iterations, then vary the deflation window
    // LAPACK❗✔️: *           .      size.   At first, because, larger windows are,
    // LAPACK❗✔️: *           .      in general, more powerful than smaller ones,
    // LAPACK❗✔️: *           .      rapidly increase the window to the maximum possible.
    // LAPACK❗✔️: *           .      Then, gradually reduce the window size. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NH = KBOT - KTOP + 1
    // LAPACK❗✔️:             NWUPBD = MIN( NH, NWMAX )
    // LAPACK❗✔️:             IF( NDFL.LT.KEXNW ) THEN
    // LAPACK❗✔️:                NW = MIN( NWUPBD, NWR )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                NW = MIN( NWUPBD, 2*NW )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( NW.LT.NWMAX ) THEN
    // LAPACK❗✔️:                IF( NW.GE.NH-1 ) THEN
    // LAPACK❗✔️:                   NW = NH
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   KWTOP = KBOT - NW + 1
    // LAPACK❗✔️:                   IF( ABS( H( KWTOP, KWTOP-1 ) ).GT.
    // LAPACK❗✔️:      $                ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( NDFL.LT.KEXNW ) THEN
    // LAPACK❗✔️:                NDEC = -1
    // LAPACK❗✔️:             ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
    // LAPACK❗✔️:                NDEC = NDEC + 1
    // LAPACK❗✔️:                IF( NW-NDEC.LT.2 )
    // LAPACK❗✔️:      $            NDEC = 0
    // LAPACK❗✔️:                NW = NW - NDEC
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Aggressive early deflation:
    // LAPACK❗✔️: *           .    split workspace under the subdiagonal into
    // LAPACK❗✔️: *           .      - an nw-by-nw work array V in the lower
    // LAPACK❗✔️: *           .        left-hand-corner,
    // LAPACK❗✔️: *           .      - an NW-by-at-least-NW-but-more-is-better
    // LAPACK❗✔️: *           .        (NW-by-NHO) horizontal work array along
    // LAPACK❗✔️: *           .        the bottom edge,
    // LAPACK❗✔️: *           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
    // LAPACK❗✔️: *           .        vertical work array along the left-hand-edge.
    // LAPACK❗✔️: *           .        ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             KV = N - NW + 1
    // LAPACK❗✔️:             KT = NW + 1
    // LAPACK❗✔️:             NHO = ( N-NW-1 ) - KT + 1
    // LAPACK❗✔️:             KWV = NW + 2
    // LAPACK❗✔️:             NVE = ( N-NW ) - KWV + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Aggressive early deflation ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
    // LAPACK❗✔️:      $                   IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH,
    // LAPACK❗✔️:      $                   NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH,
    // LAPACK❗✔️:      $                   WORK, LWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Adjust KBOT accounting for new deflations. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             KBOT = KBOT - LD
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== KS points to the shifts. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             KS = KBOT - LS + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Skip an expensive QR sweep if there is a (partly
    // LAPACK❗✔️: *           .    heuristic) reason to expect that many eigenvalues
    // LAPACK❗✔️: *           .    will deflate without it.  Here, the QR sweep is
    // LAPACK❗✔️: *           .    skipped if many eigenvalues have just been deflated
    // LAPACK❗✔️: *           .    or if the remaining active block is small.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
    // LAPACK❗✔️:      $          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== NS = nominal number of simultaneous shifts.
    // LAPACK❗✔️: *              .    This may be lowered (slightly) if DLAQR2
    // LAPACK❗✔️: *              .    did not provide that many shifts. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
    // LAPACK❗✔️:                NS = NS - MOD( NS, 2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== If there have been no deflations
    // LAPACK❗✔️: *              .    in a multiple of KEXSH iterations,
    // LAPACK❗✔️: *              .    then try exceptional shifts.
    // LAPACK❗✔️: *              .    Otherwise use shifts provided by
    // LAPACK❗✔️: *              .    DLAQR2 above or from the eigenvalues
    // LAPACK❗✔️: *              .    of a trailing principal submatrix. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
    // LAPACK❗✔️:                   KS = KBOT - NS + 1
    // LAPACK❗✔️:                   DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
    // LAPACK❗✔️:                      SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
    // LAPACK❗✔️:                      AA = WILK1*SS + H( I, I )
    // LAPACK❗✔️:                      BB = SS
    // LAPACK❗✔️:                      CC = WILK2*SS
    // LAPACK❗✔️:                      DD = AA
    // LAPACK❗✔️:                      CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ),
    // LAPACK❗✔️:      $                            WR( I ), WI( I ), CS, SN )
    // LAPACK❗✔️:    30             CONTINUE
    // LAPACK❗✔️:                   IF( KS.EQ.KTOP ) THEN
    // LAPACK❗✔️:                      WR( KS+1 ) = H( KS+1, KS+1 )
    // LAPACK❗✔️:                      WI( KS+1 ) = ZERO
    // LAPACK❗✔️:                      WR( KS ) = WR( KS+1 )
    // LAPACK❗✔️:                      WI( KS ) = WI( KS+1 )
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Got NS/2 or fewer shifts? Use DLAHQR
    // LAPACK❗✔️: *                 .    on a trailing principal submatrix to
    // LAPACK❗✔️: *                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
    // LAPACK❗✔️: *                 .    there is enough space below the subdiagonal
    // LAPACK❗✔️: *                 .    to fit an NS-by-NS scratch array.) ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   IF( KBOT-KS+1.LE.NS / 2 ) THEN
    // LAPACK❗✔️:                      KS = KBOT - NS + 1
    // LAPACK❗✔️:                      KT = N - NS + 1
    // LAPACK❗✔️:                      CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH,
    // LAPACK❗✔️:      $                            H( KT, 1 ), LDH )
    // LAPACK❗✔️:                      CALL DLAHQR( .false., .false., NS, 1, NS,
    // LAPACK❗✔️:      $                            H( KT, 1 ), LDH, WR( KS ), WI( KS ),
    // LAPACK❗✔️:      $                            1, 1, ZDUM, 1, INF )
    // LAPACK❗✔️:                      KS = KS + INF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    ==== In case of a rare QR failure use
    // LAPACK❗✔️: *                    .    eigenvalues of the trailing 2-by-2
    // LAPACK❗✔️: *                    .    principal submatrix.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      IF( KS.GE.KBOT ) THEN
    // LAPACK❗✔️:                         AA = H( KBOT-1, KBOT-1 )
    // LAPACK❗✔️:                         CC = H( KBOT, KBOT-1 )
    // LAPACK❗✔️:                         BB = H( KBOT-1, KBOT )
    // LAPACK❗✔️:                         DD = H( KBOT, KBOT )
    // LAPACK❗✔️:                         CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ),
    // LAPACK❗✔️:      $                               WI( KBOT-1 ), WR( KBOT ),
    // LAPACK❗✔️:      $                               WI( KBOT ), CS, SN )
    // LAPACK❗✔️:                         KS = KBOT - 1
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   IF( KBOT-KS+1.GT.NS ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    ==== Sort the shifts (Helps a little)
    // LAPACK❗✔️: *                    .    Bubble sort keeps complex conjugate
    // LAPACK❗✔️: *                    .    pairs together. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      SORTED = .false.
    // LAPACK❗✔️:                      DO 50 K = KBOT, KS + 1, -1
    // LAPACK❗✔️:                         IF( SORTED )
    // LAPACK❗✔️:      $                     GO TO 60
    // LAPACK❗✔️:                         SORTED = .true.
    // LAPACK❗✔️:                         DO 40 I = KS, K - 1
    // LAPACK❗✔️:                            IF( ABS( WR( I ) )+ABS( WI( I ) ).LT.
    // LAPACK❗✔️:      $                         ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
    // LAPACK❗✔️:                               SORTED = .false.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                               SWAP = WR( I )
    // LAPACK❗✔️:                               WR( I ) = WR( I+1 )
    // LAPACK❗✔️:                               WR( I+1 ) = SWAP
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                               SWAP = WI( I )
    // LAPACK❗✔️:                               WI( I ) = WI( I+1 )
    // LAPACK❗✔️:                               WI( I+1 ) = SWAP
    // LAPACK❗✔️:                            END IF
    // LAPACK❗✔️:    40                   CONTINUE
    // LAPACK❗✔️:    50                CONTINUE
    // LAPACK❗✔️:    60                CONTINUE
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Shuffle shifts into pairs of real shifts
    // LAPACK❗✔️: *                 .    and pairs of complex conjugate shifts
    // LAPACK❗✔️: *                 .    assuming complex conjugate shifts are
    // LAPACK❗✔️: *                 .    already adjacent to one another. (Yes,
    // LAPACK❗✔️: *                 .    they are.)  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   DO 70 I = KBOT, KS + 2, -2
    // LAPACK❗✔️:                      IF( WI( I ).NE.-WI( I-1 ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         SWAP = WR( I )
    // LAPACK❗✔️:                         WR( I ) = WR( I-1 )
    // LAPACK❗✔️:                         WR( I-1 ) = WR( I-2 )
    // LAPACK❗✔️:                         WR( I-2 ) = SWAP
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                         SWAP = WI( I )
    // LAPACK❗✔️:                         WI( I ) = WI( I-1 )
    // LAPACK❗✔️:                         WI( I-1 ) = WI( I-2 )
    // LAPACK❗✔️:                         WI( I-2 ) = SWAP
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:    70             CONTINUE
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== If there are only two shifts and both are
    // LAPACK❗✔️: *              .    real, then use only one.  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( KBOT-KS+1.EQ.2 ) THEN
    // LAPACK❗✔️:                   IF( WI( KBOT ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                      IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT.
    // LAPACK❗✔️:      $                   ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
    // LAPACK❗✔️:                         WR( KBOT-1 ) = WR( KBOT )
    // LAPACK❗✔️:                      ELSE
    // LAPACK❗✔️:                         WR( KBOT ) = WR( KBOT-1 )
    // LAPACK❗✔️:                      END IF
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Use up to NS of the the smallest magnitude
    // LAPACK❗✔️: *              .    shifts.  If there aren't NS shifts available,
    // LAPACK❗✔️: *              .    then use them all, possibly dropping one to
    // LAPACK❗✔️: *              .    make the number of shifts even. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NS = MIN( NS, KBOT-KS+1 )
    // LAPACK❗✔️:                NS = NS - MOD( NS, 2 )
    // LAPACK❗✔️:                KS = KBOT - NS + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Small-bulge multi-shift QR sweep:
    // LAPACK❗✔️: *              .    split workspace under the subdiagonal into
    // LAPACK❗✔️: *              .    - a KDU-by-KDU work array U in the lower
    // LAPACK❗✔️: *              .      left-hand-corner,
    // LAPACK❗✔️: *              .    - a KDU-by-at-least-KDU-but-more-is-better
    // LAPACK❗✔️: *              .      (KDU-by-NHo) horizontal work array WH along
    // LAPACK❗✔️: *              .      the bottom edge,
    // LAPACK❗✔️: *              .    - and an at-least-KDU-but-more-is-better-by-KDU
    // LAPACK❗✔️: *              .      (NVE-by-KDU) vertical work WV arrow along
    // LAPACK❗✔️: *              .      the left-hand-edge. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                KDU = 2*NS
    // LAPACK❗✔️:                KU = N - KDU + 1
    // LAPACK❗✔️:                KWH = KDU + 1
    // LAPACK❗✔️:                NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
    // LAPACK❗✔️:                KWV = KDU + 4
    // LAPACK❗✔️:                NVE = N - KDU - KWV + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== Small-bulge multi-shift QR sweep ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,
    // LAPACK❗✔️:      $                      WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z,
    // LAPACK❗✔️:      $                      LDZ, WORK, 3, H( KU, 1 ), LDH, NVE,
    // LAPACK❗✔️:      $                      H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Note progress (or the lack of it). ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LD.GT.0 ) THEN
    // LAPACK❗✔️:                NDFL = 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                NDFL = NDFL + 1
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== End of main loop ====
    // LAPACK❗✔️:    80    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Iteration limit exceeded.  Set INFO to show where
    // LAPACK❗✔️: *        .    the problem occurred and exit. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          INFO = KBOT
    // LAPACK❗✔️:    90    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Return the optimal value of LWORK. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = DBLE( LWKOPT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== End of DLAQR4 ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if n == 0 {
        return Ok(0);
    }
    if n <= 15 {
        return Ok(super::schur::lahqr(
            wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz,
        ));
    }
    let nh0 = ihi - ilo + 1;
    let recommended = recommended_shifts(nh0);
    let nwr = 2
        .max(if nh0 <= 500 {
            recommended
        } else {
            3 * recommended / 2
        })
        .min(nh0)
        .min((n - 1) / 3);
    let mut nsr = recommended.min((n - 3) / 6).min(ihi - ilo);
    nsr = 2.max(nsr - nsr % 2);
    let nmin = 75usize;
    let nibble = 14usize;
    let kacc22 = if recommended >= 14 { 2 } else { 0 };
    let nwmax = ((n - 1) / 3).min(work.len() / 2);
    let mut nw = nwmax;
    let mut nsmax = ((n - 3) / 6).min(2 * work.len() / 3);
    nsmax -= nsmax % 2;
    let mut ndfl = 1;
    let mut ndec: isize = -1;
    let mut kbot = ihi;
    let itmax = 30 * 10.max(nh0);
    for _it in 1..=itmax {
        if kbot < ilo {
            return Ok(0);
        }
        let mut ktop = ilo;
        for k in (ilo + 1..=kbot).rev() {
            if h[ix(k, k - 1, ldh)] == 0.0 {
                ktop = k;
                break;
            }
        }
        let nh = kbot - ktop + 1;
        let nwupbd = nh.min(nwmax);
        nw = if ndfl < 5 {
            nwupbd.min(nwr)
        } else {
            nwupbd.min(2 * nw)
        };
        if nw < nwmax {
            if nw >= nh - 1 {
                nw = nh;
            } else {
                let kwtop = kbot - nw + 1;
                if h[ix(kwtop, kwtop - 1, ldh)].abs() > h[ix(kwtop - 1, kwtop - 2, ldh)].abs() {
                    nw += 1;
                }
            }
        }
        if ndfl < 5 {
            ndec = -1;
        } else if ndec >= 0 || nw >= nwupbd {
            ndec += 1;
            if nw as isize - ndec < 2 {
                ndec = 0;
            }
            nw -= ndec as usize;
        }
        let kv = n - nw + 1;
        let kt = nw + 1;
        let nho = n - 2 * nw - 1;
        let kwv = nw + 2;
        let nve = n - 2 * nw - 1;
        let (ls, ld) = aed(
            recursive,
            wantt,
            wantz,
            n,
            ktop,
            kbot,
            nw,
            h,
            ldh,
            iloz,
            ihiz,
            z,
            ldz,
            wr,
            wi,
            ix(kv, 1, ldh),
            ldh,
            nho,
            ix(kv, kt, ldh),
            ldh,
            nve,
            ix(kwv, 1, ldh),
            ldh,
            work,
        )?;
        kbot -= ld;
        let mut ks = kbot - ls + 1;
        if ld == 0 || (100 * ld <= nw * nibble && kbot + 1 - ktop > nmin.min(nwmax)) {
            let mut ns = nsmax.min(nsr).min(2.max(kbot - ktop));
            ns -= ns % 2;
            if ndfl % 6 == 0 {
                ks = kbot - ns + 1;
                let lower = (ks + 1).max(ktop + 2);
                let mut i = kbot;
                while i >= lower {
                    let ss = h[ix(i, i - 1, ldh)].abs() + h[ix(i - 1, i - 2, ldh)].abs();
                    let aa = 0.75 * ss + h[ix(i, i, ldh)];
                    let p = super::lanv2(aa, ss, -0.4375 * ss, aa);
                    wr[i - 2] = p.real[0];
                    wi[i - 2] = p.imag[0];
                    wr[i - 1] = p.real[1];
                    wi[i - 1] = p.imag[1];
                    i -= 2;
                }
                if ks == ktop {
                    wr[ks] = h[ix(ks + 1, ks + 1, ldh)];
                    wi[ks] = 0.0;
                    wr[ks - 1] = wr[ks];
                    wi[ks - 1] = wi[ks];
                }
            } else {
                if kbot - ks + 1 <= ns / 2 {
                    ks = kbot - ns + 1;
                    let kt = n - ns + 1;
                    for j in 1..=ns {
                        for i in 1..=ns {
                            h[ix(kt + i - 1, j, ldh)] = h[ix(ks + i - 1, ks + j - 1, ldh)];
                        }
                    }
                    let inf = if recursive && ns > nmin {
                        qr(
                            false,
                            false,
                            false,
                            ns,
                            1,
                            ns,
                            &mut h[ix(kt, 1, ldh)..],
                            ldh,
                            &mut wr[ks - 1..],
                            &mut wi[ks - 1..],
                            1,
                            1,
                            &mut [0.0],
                            1,
                            work,
                        )?
                    } else {
                        super::schur::lahqr(
                            false,
                            false,
                            ns,
                            1,
                            ns,
                            &mut h[ix(kt, 1, ldh)..],
                            ldh,
                            &mut wr[ks - 1..],
                            &mut wi[ks - 1..],
                            1,
                            1,
                            &mut [0.0],
                            1,
                        )
                    };
                    ks += inf;
                    if ks >= kbot {
                        let p = super::lanv2(
                            h[ix(kbot - 1, kbot - 1, ldh)],
                            h[ix(kbot - 1, kbot, ldh)],
                            h[ix(kbot, kbot - 1, ldh)],
                            h[ix(kbot, kbot, ldh)],
                        );
                        wr[kbot - 2] = p.real[0];
                        wi[kbot - 2] = p.imag[0];
                        wr[kbot - 1] = p.real[1];
                        wi[kbot - 1] = p.imag[1];
                        ks = kbot - 1;
                    }
                }
                if kbot - ks + 1 > ns {
                    let mut sorted = false;
                    for k in (ks + 1..=kbot).rev() {
                        if sorted {
                            break;
                        }
                        sorted = true;
                        for i in ks..k {
                            if wr[i - 1].abs() + wi[i - 1].abs() < wr[i].abs() + wi[i].abs() {
                                sorted = false;
                                wr.swap(i - 1, i);
                                wi.swap(i - 1, i);
                            }
                        }
                    }
                }
                let mut i = kbot;
                while i >= ks + 2 {
                    if wi[i - 1] != -wi[i - 2] {
                        let swap = wr[i - 1];
                        wr[i - 1] = wr[i - 2];
                        wr[i - 2] = wr[i - 3];
                        wr[i - 3] = swap;
                        let swap = wi[i - 1];
                        wi[i - 1] = wi[i - 2];
                        wi[i - 2] = wi[i - 3];
                        wi[i - 3] = swap;
                    }
                    i -= 2;
                }
            }
            if kbot - ks + 1 == 2 && wi[kbot - 1] == 0.0 {
                if (wr[kbot - 1] - h[ix(kbot, kbot, ldh)]).abs()
                    < (wr[kbot - 2] - h[ix(kbot, kbot, ldh)]).abs()
                {
                    wr[kbot - 2] = wr[kbot - 1];
                } else {
                    wr[kbot - 1] = wr[kbot - 2];
                }
            }
            ns = ns.min(kbot - ks + 1);
            ns -= ns % 2;
            ks = kbot - ns + 1;
            let kdu = 2 * ns;
            let ku = n - kdu + 1;
            let kwh = kdu + 1;
            let nho = n - 2 * kdu - 3;
            let kwv = kdu + 4;
            let nve = n - 2 * kdu - 3;
            laqr5(
                wantt,
                wantz,
                kacc22,
                n,
                ktop,
                kbot,
                ns,
                &mut wr[ks - 1..],
                &mut wi[ks - 1..],
                h,
                ldh,
                iloz,
                ihiz,
                z,
                ldz,
                work,
                3,
                ix(ku, 1, ldh),
                ldh,
                nve,
                ix(kwv, 1, ldh),
                ldh,
                nho,
                ix(ku, kwh, ldh),
                ldh,
            )?;
        }
        ndfl = if ld > 0 { 1 } else { ndfl + 1 };
    }
    Ok(kbot)
}
