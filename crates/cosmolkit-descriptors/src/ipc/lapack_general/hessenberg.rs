//! Source Hessenberg reduction and internal block reflectors.
use super::super::lapack_blas as blas;
use super::super::lapack_blas_level3 as blas3;
use super::{GeneralError, ix};
use blas3::{
    Diagonal::{NonUnit, Unit},
    Side::{Left, Right},
    Transpose::{NoTrans, Trans},
    Triangle::{Lower, Upper},
};

/// Fixed INCV=1; all DGEEV/AED callers use stored column reflectors.
#[allow(clippy::too_many_arguments)]
pub(super) fn larf1f(
    left: bool,
    m: usize,
    n: usize,
    v: &[f64],
    tau: f64,
    c: &mut [f64],
    ldc: usize,
    work: &mut [f64],
) {
    // LAPACK❗✔️:       SUBROUTINE DLARF1F( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          SIDE
    // LAPACK❗✔️:       INTEGER            INCV, LDC, M, N
    // LAPACK❗✔️:       DOUBLE PRECISION   TAU
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE, ZERO
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            APPLYLEFT
    // LAPACK❗✔️:       INTEGER            I, LASTV, LASTC
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DGEMV, DGER, DAXPY, DSCAL
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       INTEGER            ILADLR, ILADLC
    // LAPACK❗✔️:       EXTERNAL           LSAME, ILADLR, ILADLC
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       APPLYLEFT = LSAME( SIDE, 'L' )
    // LAPACK❗✔️:       LASTV = 1
    // LAPACK❗✔️:       LASTC = 0
    // LAPACK❗✔️:       IF( TAU.NE.ZERO ) THEN
    // LAPACK❗✔️: !     Set up variables for scanning V.  LASTV begins pointing to the end
    // LAPACK❗✔️: !     of V.
    // LAPACK❗✔️:          IF( APPLYLEFT ) THEN
    // LAPACK❗✔️:             LASTV = M
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             LASTV = N
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( INCV.GT.0 ) THEN
    // LAPACK❗✔️:             I = 1 + (LASTV-1) * INCV
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             I = 1
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: !     Look for the last non-zero row in V.
    // LAPACK❗✔️: !        Since we are assuming that V(1) = 1, and it is not stored, so we
    // LAPACK❗✔️: !        shouldn't access it.
    // LAPACK❗✔️:          DO WHILE( LASTV.GT.1 .AND. V( I ).EQ.ZERO )
    // LAPACK❗✔️:             LASTV = LASTV - 1
    // LAPACK❗✔️:             I = I - INCV
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:          IF( APPLYLEFT ) THEN
    // LAPACK❗✔️: !     Scan for the last non-zero column in C(1:lastv,:).
    // LAPACK❗✔️:             LASTC = ILADLC(LASTV, N, C, LDC)
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: !     Scan for the last non-zero row in C(:,1:lastv).
    // LAPACK❗✔️:             LASTC = ILADLR(M, LASTV, C, LDC)
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( LASTC.EQ.0 ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( APPLYLEFT ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Form  H * C
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ! Check if lastv = 1. This means v = 1, So we just need to compute
    // LAPACK❗✔️:          ! C := HC = (1-\tau)C.
    // LAPACK❗✔️:          IF( LASTV.EQ.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           C(1,1:lastc) := ( 1 - tau ) * C(1,1:lastc)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DSCAL(LASTC, ONE - TAU, C, LDC)
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ! w(1:lastc,1) := C(2:lastv,1:lastc)**T * v(2:lastv,1)
    // LAPACK❗✔️:             CALL DGEMV( 'Transpose', LASTV-1, LASTC, ONE, C(1+1,1),
    // LAPACK❗✔️:      $                  LDC, V(1+INCV), INCV, ZERO, WORK, 1)
    // LAPACK❗✔️:             ! w(1:lastc,1) += C(1,1:lastc)**T * v(1,1) = C(1,1:lastc)**T
    // LAPACK❗✔️:             CALL DAXPY(LASTC, ONE, C, LDC, WORK, 1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ! C(1, 1:lastc)   := C(...) - tau * v(1,1) * w(1:lastc,1)**T
    // LAPACK❗✔️:             !                  = C(...) - tau * w(1:lastc,1)**T
    // LAPACK❗✔️:             CALL DAXPY(LASTC, -TAU, WORK, 1, C, LDC)
    // LAPACK❗✔️:             ! C(2:lastv,1:lastc) := C(...) - tau * v(2:lastv,1)*w(1:lastc,1)**T
    // LAPACK❗✔️:             CALL DGER(LASTV-1, LASTC, -TAU, V(1+INCV), INCV, WORK, 1,
    // LAPACK❗✔️:      $                  C(1+1,1), LDC)
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Form  C * H
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ! Check if n = 1. This means v = 1, so we just need to compute
    // LAPACK❗✔️:          ! C := CH = C(1-\tau).
    // LAPACK❗✔️:          IF( LASTV.EQ.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           C(1:lastc,1) := ( 1 - tau ) * C(1:lastc,1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DSCAL(LASTC, ONE - TAU, C, 1)
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ! w(1:lastc,1) := C(1:lastc,2:lastv) * v(2:lastv,1)
    // LAPACK❗✔️:             CALL DGEMV( 'No transpose', LASTC, LASTV-1, ONE,
    // LAPACK❗✔️:      $         C(1,1+1), LDC, V(1+INCV), INCV, ZERO, WORK, 1 )
    // LAPACK❗✔️:             ! w(1:lastc,1) += C(1:lastc,1) v(1,1) = C(1:lastc,1)
    // LAPACK❗✔️:             CALL DAXPY(LASTC, ONE, C, 1, WORK, 1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ! C(1:lastc,1)     := C(...) - tau * w(1:lastc,1) * v(1,1)**T
    // LAPACK❗✔️:             !                   = C(...) - tau * w(1:lastc,1)
    // LAPACK❗✔️:             CALL DAXPY(LASTC, -TAU, WORK, 1, C, 1)
    // LAPACK❗✔️:             ! C(1:lastc,2:lastv) := C(...) - tau * w(1:lastc,1) * v(2:lastv)**T
    // LAPACK❗✔️:             CALL DGER( LASTC, LASTV-1, -TAU, WORK, 1, V(1+INCV),
    // LAPACK❗✔️:      $                  INCV, C(1,1+1), LDC )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLARF1F
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if m == 0 || n == 0 || tau == 0.0 {
        return;
    }
    let mut lastv = if left { m } else { n };
    while lastv > 1 && v[lastv - 1] == 0.0 {
        lastv -= 1;
    }
    let lastc = if left {
        last_column(lastv, n, c, ldc)
    } else {
        last_row(m, lastv, c, ldc)
    };
    if lastc == 0 {
        return;
    }
    if left {
        if lastv == 1 {
            blas::scal(lastc, 1.0 - tau, c, ldc);
        } else {
            blas::gemv(
                true,
                lastv - 1,
                lastc,
                1.0,
                &c[1..],
                ldc,
                &v[1..],
                1,
                0.0,
                work,
                1,
            );
            blas::axpy(lastc, 1.0, c, ldc, work, 1);
            blas::axpy(lastc, -tau, work, 1, c, ldc);
            blas::ger(
                lastv - 1,
                lastc,
                -tau,
                &v[1..],
                1,
                work,
                1,
                &mut c[1..],
                ldc,
            );
        }
    } else if lastv == 1 {
        blas::scal(lastc, 1.0 - tau, c, 1);
    } else {
        blas::gemv(
            false,
            lastc,
            lastv - 1,
            1.0,
            &c[ldc..],
            ldc,
            &v[1..],
            1,
            0.0,
            work,
            1,
        );
        blas::axpy(lastc, 1.0, c, 1, work, 1);
        blas::axpy(lastc, -tau, work, 1, c, 1);
        blas::ger(
            lastc,
            lastv - 1,
            -tau,
            work,
            1,
            &v[1..],
            1,
            &mut c[ldc..],
            ldc,
        );
    }
}

fn last_column(m: usize, n: usize, a: &[f64], lda: usize) -> usize {
    // LAPACK❗✔️:       INTEGER FUNCTION ILADLC( M, N, A, LDA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            M, N, LDA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION ZERO
    // LAPACK❗✔️:       PARAMETER ( ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER I
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick test for the common case where one corner is non-zero.
    // LAPACK❗✔️:       IF( N.EQ.0 ) THEN
    // LAPACK❗✔️:          ILADLC = N
    // LAPACK❗✔️:       ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
    // LAPACK❗✔️:          ILADLC = N
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *     Now scan each column from the end, returning with the first non-zero.
    // LAPACK❗✔️:          DO ILADLC = N, 1, -1
    // LAPACK❗✔️:             DO I = 1, M
    // LAPACK❗✔️:                IF( A(I, ILADLC).NE.ZERO ) RETURN
    // LAPACK❗✔️:             END DO
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️:       END
    if n == 0 {
        return 0;
    }
    if a[ix(1, n, lda)] != 0.0 || a[ix(m, n, lda)] != 0.0 {
        return n;
    }
    for j in (1..=n).rev() {
        for i in 1..=m {
            if a[ix(i, j, lda)] != 0.0 {
                return j;
            }
        }
    }
    0
}
fn last_row(m: usize, n: usize, a: &[f64], lda: usize) -> usize {
    // LAPACK❗✔️:       INTEGER FUNCTION ILADLR( M, N, A, LDA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            M, N, LDA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION ZERO
    // LAPACK❗✔️:       PARAMETER ( ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER I, J
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick test for the common case where one corner is non-zero.
    // LAPACK❗✔️:       IF( M.EQ.0 ) THEN
    // LAPACK❗✔️:          ILADLR = M
    // LAPACK❗✔️:       ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
    // LAPACK❗✔️:          ILADLR = M
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *     Scan up each column tracking the last zero row seen.
    // LAPACK❗✔️:          ILADLR = 0
    // LAPACK❗✔️:          DO J = 1, N
    // LAPACK❗✔️:             I=M
    // LAPACK❗✔️:             DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
    // LAPACK❗✔️:                I=I-1
    // LAPACK❗✔️:             ENDDO
    // LAPACK❗✔️:             ILADLR = MAX( ILADLR, I )
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️:       END
    if m == 0 {
        return 0;
    }
    if a[ix(m, 1, lda)] != 0.0 || a[ix(m, n, lda)] != 0.0 {
        return m;
    }
    let mut result = 0;
    for j in 1..=n {
        let mut i = m;
        while i > 0 && a[ix(i, j, lda)] == 0.0 {
            i -= 1;
        }
        result = result.max(i);
    }
    result
}

fn gehd2(
    n: usize,
    ilo: usize,
    ihi: usize,
    a: &mut [f64],
    lda: usize,
    tau: &mut [f64],
    work: &mut [f64],
) {
    // LAPACK❗✔️:       SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHI, ILO, INFO, LDA, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLARF1F, DLARFG, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DGEHD2', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 I = ILO, IHI - 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
    // LAPACK❗✔️:      $                TAU( I ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Apply H(i) to A(1:ihi,i+1:ihi) from the right
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARF1F( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ),
    // LAPACK❗✔️:      $               A( 1, I+1 ), LDA, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Apply H(i) to A(i+1:ihi,i+1:n) from the left
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARF1F( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ),
    // LAPACK❗✔️:      $               A( I+1, I+1 ), LDA, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DGEHD2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    for i in ilo..ihi {
        let start = ix(i + 1, i, lda);
        tau[i - 1] = super::super::lapack_symmetric::larfg(&mut a[start..start + ihi - i]);
        // V is in column I; both updated submatrices start at column I+1.
        // Split column storage rather than create aliased strided suffixes.
        let (columns_before, columns_after) = a.split_at_mut(i * lda);
        let v = &columns_before[start..];
        larf1f(false, ihi, ihi - i, v, tau[i - 1], columns_after, lda, work);
        larf1f(
            true,
            ihi - i,
            n - i,
            v,
            tau[i - 1],
            &mut columns_after[i..],
            lda,
            work,
        );
    }
}

/// Forward, columnwise block reflector; both SIDE modes are internally reachable.
#[allow(clippy::too_many_arguments)]
fn larfb(
    left: bool,
    transpose: bool,
    m: usize,
    n: usize,
    k: usize,
    v: &[f64],
    ldv: usize,
    t: &[f64],
    ldt: usize,
    c: &mut [f64],
    ldc: usize,
    work: &mut [f64],
    ldw: usize,
) -> Result<(), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
    // LAPACK❗✔️:      $                   T, LDT, C, LDC, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          DIRECT, SIDE, STOREV, TRANS
    // LAPACK❗✔️:       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
    // LAPACK❗✔️:      $                   WORK( LDWORK, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       CHARACTER          TRANST
    // LAPACK❗✔️:       INTEGER            I, J
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DCOPY, DGEMM, DTRMM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( M.LE.0 .OR. N.LE.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( TRANS, 'N' ) ) THEN
    // LAPACK❗✔️:          TRANST = 'T'
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          TRANST = 'N'
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( STOREV, 'C' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LSAME( DIRECT, 'F' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Let  V =  ( V1 )    (first K rows)
    // LAPACK❗✔️: *                     ( V2 )
    // LAPACK❗✔️: *           where  V1  is unit lower triangular.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LSAME( SIDE, 'L' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  H * C  or  H**T * C  where  C = ( C1 )
    // LAPACK❗✔️: *                                                    ( C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 10 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
    // LAPACK❗✔️:    10          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
    // LAPACK❗✔️:      $                     K, ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C2**T * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
    // LAPACK❗✔️:      $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
    // LAPACK❗✔️:      $                        ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T**T  or  W * T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - V * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C2 := C2 - V2 * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
    // LAPACK❗✔️:      $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
    // LAPACK❗✔️:      $                        C( K+1, 1 ), LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
    // LAPACK❗✔️:      $                     ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C1 := C1 - W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 30 J = 1, K
    // LAPACK❗✔️:                   DO 20 I = 1, N
    // LAPACK❗✔️:                      C( J, I ) = C( J, I ) - WORK( I, J )
    // LAPACK❗✔️:    20             CONTINUE
    // LAPACK❗✔️:    30          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ELSE IF( LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 40 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
    // LAPACK❗✔️:    40          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
    // LAPACK❗✔️:      $                     K, ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C2 * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
    // LAPACK❗✔️:      $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
    // LAPACK❗✔️:      $                        ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T  or  W * T**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - W * V**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C2 := C2 - W * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
    // LAPACK❗✔️:      $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
    // LAPACK❗✔️:      $                        C( 1, K+1 ), LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
    // LAPACK❗✔️:      $                     ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C1 := C1 - W
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 60 J = 1, K
    // LAPACK❗✔️:                   DO 50 I = 1, M
    // LAPACK❗✔️:                      C( I, J ) = C( I, J ) - WORK( I, J )
    // LAPACK❗✔️:    50             CONTINUE
    // LAPACK❗✔️:    60          CONTINUE
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Let  V =  ( V1 )
    // LAPACK❗✔️: *                     ( V2 )    (last K rows)
    // LAPACK❗✔️: *           where  V2  is unit upper triangular.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LSAME( SIDE, 'L' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  H * C  or  H**T * C  where  C = ( C1 )
    // LAPACK❗✔️: *                                                    ( C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 70 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
    // LAPACK❗✔️:    70          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
    // LAPACK❗✔️:      $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C1**T * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
    // LAPACK❗✔️:      $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T**T  or  W * T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - V * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C1 := C1 - V1 * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
    // LAPACK❗✔️:      $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
    // LAPACK❗✔️:      $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C2 := C2 - W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 90 J = 1, K
    // LAPACK❗✔️:                   DO 80 I = 1, N
    // LAPACK❗✔️:                      C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
    // LAPACK❗✔️:    80             CONTINUE
    // LAPACK❗✔️:    90          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ELSE IF( LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 100 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
    // LAPACK❗✔️:   100          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
    // LAPACK❗✔️:      $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C1 * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K,
    // LAPACK❗✔️:      $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T  or  W * T**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - W * V**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C1 := C1 - W * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K,
    // LAPACK❗✔️:      $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
    // LAPACK❗✔️:      $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C2 := C2 - W
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 120 J = 1, K
    // LAPACK❗✔️:                   DO 110 I = 1, M
    // LAPACK❗✔️:                      C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
    // LAPACK❗✔️:   110             CONTINUE
    // LAPACK❗✔️:   120          CONTINUE
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( LSAME( STOREV, 'R' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LSAME( DIRECT, 'F' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Let  V =  ( V1  V2 )    (V1: first K columns)
    // LAPACK❗✔️: *           where  V1  is unit upper triangular.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LSAME( SIDE, 'L' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  H * C  or  H**T * C  where  C = ( C1 )
    // LAPACK❗✔️: *                                                    ( C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 130 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
    // LAPACK❗✔️:   130          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
    // LAPACK❗✔️:      $                     ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C2**T * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
    // LAPACK❗✔️:      $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
    // LAPACK❗✔️:      $                        WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T**T  or  W * T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - V**T * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C2 := C2 - V2**T * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
    // LAPACK❗✔️:      $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
    // LAPACK❗✔️:      $                        C( K+1, 1 ), LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
    // LAPACK❗✔️:      $                     K, ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C1 := C1 - W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 150 J = 1, K
    // LAPACK❗✔️:                   DO 140 I = 1, N
    // LAPACK❗✔️:                      C( J, I ) = C( J, I ) - WORK( I, J )
    // LAPACK❗✔️:   140             CONTINUE
    // LAPACK❗✔️:   150          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ELSE IF( LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 160 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
    // LAPACK❗✔️:   160          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
    // LAPACK❗✔️:      $                     ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C2 * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
    // LAPACK❗✔️:      $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
    // LAPACK❗✔️:      $                        ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T  or  W * T**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - W * V
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C2 := C2 - W * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
    // LAPACK❗✔️:      $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
    // LAPACK❗✔️:      $                        C( 1, K+1 ), LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
    // LAPACK❗✔️:      $                     K, ONE, V, LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C1 := C1 - W
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 180 J = 1, K
    // LAPACK❗✔️:                   DO 170 I = 1, M
    // LAPACK❗✔️:                      C( I, J ) = C( I, J ) - WORK( I, J )
    // LAPACK❗✔️:   170             CONTINUE
    // LAPACK❗✔️:   180          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Let  V =  ( V1  V2 )    (V2: last K columns)
    // LAPACK❗✔️: *           where  V2  is unit lower triangular.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LSAME( SIDE, 'L' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  H * C  or  H**T * C  where  C = ( C1 )
    // LAPACK❗✔️: *                                                    ( C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 190 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
    // LAPACK❗✔️:   190          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
    // LAPACK❗✔️:      $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C1**T * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
    // LAPACK❗✔️:      $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T**T  or  W * T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - V**T * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( M.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C1 := C1 - V1**T * W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
    // LAPACK❗✔️:      $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
    // LAPACK❗✔️:      $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C2 := C2 - W**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 210 J = 1, K
    // LAPACK❗✔️:                   DO 200 I = 1, N
    // LAPACK❗✔️:                      C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
    // LAPACK❗✔️:   200             CONTINUE
    // LAPACK❗✔️:   210          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             ELSE IF( LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Form  C * H  or  C * H'  where  C = ( C1  C2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := C2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 220 J = 1, K
    // LAPACK❗✔️:                   CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
    // LAPACK❗✔️:   220          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
    // LAPACK❗✔️:      $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 W := W + C1 * V1**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K,
    // LAPACK❗✔️:      $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * T  or  W * T**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
    // LAPACK❗✔️:      $                     ONE, T, LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C := C - W * V
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( N.GT.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 C1 := C1 - W * V1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K,
    // LAPACK❗✔️:      $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              W := W * V2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
    // LAPACK❗✔️:      $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              C1 := C1 - W
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 240 J = 1, K
    // LAPACK❗✔️:                   DO 230 I = 1, M
    // LAPACK❗✔️:                      C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
    // LAPACK❗✔️:   230             CONTINUE
    // LAPACK❗✔️:   240          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLARFB
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if m == 0 || n == 0 {
        return Ok(());
    }
    if left {
        for j in 1..=k {
            for i in 1..=n {
                work[ix(i, j, ldw)] = c[ix(j, i, ldc)];
            }
        }
        blas3::trmm(Right, Lower, NoTrans, Unit, n, k, 1.0, v, ldv, work, ldw)?;
        if m > k {
            blas3::gemm(
                Trans,
                NoTrans,
                n,
                k,
                m - k,
                1.0,
                &c[k..],
                ldc,
                &v[k..],
                ldv,
                1.0,
                work,
                ldw,
            )?;
        }
        blas3::trmm(
            Right,
            Upper,
            if transpose { NoTrans } else { Trans },
            NonUnit,
            n,
            k,
            1.0,
            t,
            ldt,
            work,
            ldw,
        )?;
        if m > k {
            blas3::gemm(
                NoTrans,
                Trans,
                m - k,
                n,
                k,
                -1.0,
                &v[k..],
                ldv,
                work,
                ldw,
                1.0,
                &mut c[k..],
                ldc,
            )?;
        }
        blas3::trmm(Right, Lower, Trans, Unit, n, k, 1.0, v, ldv, work, ldw)?;
        for j in 1..=k {
            for i in 1..=n {
                c[ix(j, i, ldc)] -= work[ix(i, j, ldw)];
            }
        }
    } else {
        for j in 1..=k {
            for i in 1..=m {
                work[ix(i, j, ldw)] = c[ix(i, j, ldc)];
            }
        }
        blas3::trmm(Right, Lower, NoTrans, Unit, m, k, 1.0, v, ldv, work, ldw)?;
        if n > k {
            blas3::gemm(
                NoTrans,
                NoTrans,
                m,
                k,
                n - k,
                1.0,
                &c[k * ldc..],
                ldc,
                &v[k..],
                ldv,
                1.0,
                work,
                ldw,
            )?;
        }
        blas3::trmm(
            Right,
            Upper,
            if transpose { Trans } else { NoTrans },
            NonUnit,
            m,
            k,
            1.0,
            t,
            ldt,
            work,
            ldw,
        )?;
        if n > k {
            blas3::gemm(
                NoTrans,
                Trans,
                m,
                n - k,
                k,
                -1.0,
                work,
                ldw,
                &v[k..],
                ldv,
                1.0,
                &mut c[k * ldc..],
                ldc,
            )?;
        }
        blas3::trmm(Right, Lower, Trans, Unit, m, k, 1.0, v, ldv, work, ldw)?;
        for j in 1..=k {
            for i in 1..=m {
                c[ix(i, j, ldc)] -= work[ix(i, j, ldw)];
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn lahr2(
    n: usize,
    k: usize,
    nb: usize,
    a: &mut [f64],
    lda: usize,
    tau: &mut [f64],
    t: &mut [f64],
    ldt: usize,
    y: &mut [f64],
    ldy: usize,
) -> Result<(), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            K, LDA, LDT, LDY, N, NB
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION  A( LDA, * ), T( LDT, NB ), TAU( NB ),
    // LAPACK❗✔️:      $                   Y( LDY, NB )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION  ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0,
    // LAPACK❗✔️:      $                     ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I
    // LAPACK❗✔️:       DOUBLE PRECISION  EI
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DAXPY, DCOPY, DGEMM, DGEMV, DLACPY,
    // LAPACK❗✔️:      $                   DLARFG, DSCAL, DTRMM, DTRMV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.1 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 I = 1, NB
    // LAPACK❗✔️:          IF( I.GT.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Update A(K+1:N,I)
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Update I-th column of A - Y * V**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY,
    // LAPACK❗✔️:      $                  A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply I - V * T**T * V**T to this column (call it b) from the
    // LAPACK❗✔️: *           left, using the last column of T as workspace
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
    // LAPACK❗✔️: *                    ( V2 )             ( b2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           where V1 is unit lower triangular
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           w := V1**T * b1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
    // LAPACK❗✔️:             CALL DTRMV( 'Lower', 'Transpose', 'UNIT',
    // LAPACK❗✔️:      $                  I-1, A( K+1, 1 ),
    // LAPACK❗✔️:      $                  LDA, T( 1, NB ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           w := w + V2**T * b2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DGEMV( 'Transpose', N-K-I+1, I-1,
    // LAPACK❗✔️:      $                  ONE, A( K+I, 1 ),
    // LAPACK❗✔️:      $                  LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           w := T**T * w
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DTRMV( 'Upper', 'Transpose', 'NON-UNIT',
    // LAPACK❗✔️:      $                  I-1, T, LDT,
    // LAPACK❗✔️:      $                  T( 1, NB ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           b2 := b2 - V2*w
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE,
    // LAPACK❗✔️:      $                  A( K+I, 1 ),
    // LAPACK❗✔️:      $                  LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           b1 := b1 - V1*w
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DTRMV( 'Lower', 'NO TRANSPOSE',
    // LAPACK❗✔️:      $                  'UNIT', I-1,
    // LAPACK❗✔️:      $                  A( K+1, 1 ), LDA, T( 1, NB ), 1 )
    // LAPACK❗✔️:             CALL DAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             A( K+I-1, I-1 ) = EI
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Generate the elementary reflector H(I) to annihilate
    // LAPACK❗✔️: *        A(K+I+1:N,I)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1,
    // LAPACK❗✔️:      $                TAU( I ) )
    // LAPACK❗✔️:          EI = A( K+I, I )
    // LAPACK❗✔️:          A( K+I, I ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute  Y(K+1:N,I)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DGEMV( 'NO TRANSPOSE', N-K, N-K-I+1,
    // LAPACK❗✔️:      $               ONE, A( K+1, I+1 ),
    // LAPACK❗✔️:      $               LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
    // LAPACK❗✔️:          CALL DGEMV( 'Transpose', N-K-I+1, I-1,
    // LAPACK❗✔️:      $               ONE, A( K+I, 1 ), LDA,
    // LAPACK❗✔️:      $               A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
    // LAPACK❗✔️:          CALL DGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE,
    // LAPACK❗✔️:      $               Y( K+1, 1 ), LDY,
    // LAPACK❗✔️:      $               T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
    // LAPACK❗✔️:          CALL DSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute T(1:I,I)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
    // LAPACK❗✔️:          CALL DTRMV( 'Upper', 'No Transpose', 'NON-UNIT',
    // LAPACK❗✔️:      $               I-1, T, LDT,
    // LAPACK❗✔️:      $               T( 1, I ), 1 )
    // LAPACK❗✔️:          T( I, I ) = TAU( I )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       A( K+NB, NB ) = EI
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Compute Y(1:K,1:NB)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
    // LAPACK❗✔️:       CALL DTRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE',
    // LAPACK❗✔️:      $            'UNIT', K, NB,
    // LAPACK❗✔️:      $            ONE, A( K+1, 1 ), LDA, Y, LDY )
    // LAPACK❗✔️:       IF( N.GT.K+NB )
    // LAPACK❗✔️:      $   CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K,
    // LAPACK❗✔️:      $               NB, N-K-NB, ONE,
    // LAPACK❗✔️:      $               A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y,
    // LAPACK❗✔️:      $               LDY )
    // LAPACK❗✔️:       CALL DTRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE',
    // LAPACK❗✔️:      $            'NON-UNIT', K, NB,
    // LAPACK❗✔️:      $            ONE, T, LDT, Y, LDY )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLAHR2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if n <= 1 {
        return Ok(());
    }
    let mut ei = 0.0;
    for i in 1..=nb {
        let ac = (i - 1) * lda;
        let tc = (i - 1) * ldt;
        let yc = (i - 1) * ldy;
        if i > 1 {
            {
                let (previous, current) = a.split_at_mut(ac);
                blas::gemv(
                    false,
                    n - k,
                    i - 1,
                    -1.0,
                    &y[k..],
                    ldy,
                    &previous[k + i - 2..],
                    lda,
                    1.0,
                    &mut current[k..],
                    1,
                );
            }
            let scratch = (nb - 1) * ldt;
            for j in 0..i - 1 {
                t[scratch + j] = a[ac + k + j];
            }
            blas::trmv(true, true, true, i - 1, &a[k..], lda, &mut t[scratch..], 1);
            blas::gemv(
                true,
                n - k - i + 1,
                i - 1,
                1.0,
                &a[k + i - 1..],
                lda,
                &a[ac + k + i - 1..],
                1,
                1.0,
                &mut t[scratch..],
                1,
            );
            {
                let (upper, last) = t.split_at_mut(scratch);
                blas::trmv(false, true, false, i - 1, upper, ldt, last, 1);
            }
            {
                let (previous, current) = a.split_at_mut(ac);
                blas::gemv(
                    false,
                    n - k - i + 1,
                    i - 1,
                    -1.0,
                    &previous[k + i - 1..],
                    lda,
                    &t[scratch..],
                    1,
                    1.0,
                    &mut current[k + i - 1..],
                    1,
                );
            }
            blas::trmv(true, false, true, i - 1, &a[k..], lda, &mut t[scratch..], 1);
            blas::axpy(i - 1, -1.0, &t[scratch..], 1, &mut a[ac + k..], 1);
            a[ix(k + i - 1, i - 1, lda)] = ei;
        }
        let start = ac + k + i - 1;
        tau[i - 1] = super::super::lapack_symmetric::larfg(&mut a[start..start + n - k - i + 1]);
        ei = a[start];
        a[start] = 1.0;
        blas::gemv(
            false,
            n - k,
            n - k - i + 1,
            1.0,
            &a[ix(k + 1, i + 1, lda)..],
            lda,
            &a[start..],
            1,
            0.0,
            &mut y[yc + k..],
            1,
        );
        blas::gemv(
            true,
            n - k - i + 1,
            i - 1,
            1.0,
            &a[k + i - 1..],
            lda,
            &a[start..],
            1,
            0.0,
            &mut t[tc..],
            1,
        );
        if i > 1 {
            let (previous, current) = y.split_at_mut(yc);
            blas::gemv(
                false,
                n - k,
                i - 1,
                -1.0,
                &previous[k..],
                ldy,
                &t[tc..],
                1,
                1.0,
                &mut current[k..],
                1,
            );
        }
        blas::scal(n - k, tau[i - 1], &mut y[yc + k..], 1);
        blas::scal(i - 1, -tau[i - 1], &mut t[tc..], 1);
        {
            let (previous, current) = t.split_at_mut(tc);
            blas::trmv(false, false, false, i - 1, previous, ldt, current, 1);
        }
        t[ix(i, i, ldt)] = tau[i - 1];
    }
    a[ix(k + nb, nb, lda)] = ei;
    for j in 1..=nb {
        for i in 1..=k {
            y[ix(i, j, ldy)] = a[ix(i, j + 1, lda)];
        }
    }
    blas3::trmm(
        Right,
        Lower,
        NoTrans,
        Unit,
        k,
        nb,
        1.0,
        &a[k..],
        lda,
        y,
        ldy,
    )?;
    if n > k + nb {
        blas3::gemm(
            NoTrans,
            NoTrans,
            k,
            nb,
            n - k - nb,
            1.0,
            &a[(1 + nb) * lda..],
            lda,
            &a[k + nb..],
            lda,
            1.0,
            y,
            ldy,
        )?;
    }
    blas3::trmm(Right, Upper, NoTrans, NonUnit, k, nb, 1.0, t, ldt, y, ldy)?;
    Ok(())
}

pub(super) fn gehrd(
    n: usize,
    ilo: usize,
    ihi: usize,
    a: &mut [f64],
    lda: usize,
    tau: &mut [f64],
    work: &mut [f64],
) -> Result<(), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       INTEGER            NBMAX, LDT, TSIZE
    // LAPACK❗✔️:       PARAMETER          ( NBMAX = 64, LDT = NBMAX+1,
    // LAPACK❗✔️:      $                     TSIZE = LDT*NBMAX )
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0,
    // LAPACK❗✔️:      $                     ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            LQUERY
    // LAPACK❗✔️:       INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB,
    // LAPACK❗✔️:      $                   NBMIN, NH, NX
    // LAPACK❗✔️:       DOUBLE PRECISION   EI
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DAXPY, DGEHD2, DGEMM, DLAHR2, DLARFB, DTRMM,
    // LAPACK❗✔️:      $                   XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       LQUERY = ( LWORK.EQ.-1 )
    // LAPACK❗✔️:       IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:          INFO = -8
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NH = IHI - ILO + 1
    // LAPACK❗✔️:       IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute the workspace requirements
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( NH.LE.1 ) THEN
    // LAPACK❗✔️:             LWKOPT = 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             NB = MIN( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI,
    // LAPACK❗✔️:      $                              -1 ) )
    // LAPACK❗✔️:             LWKOPT = N*NB + TSIZE
    // LAPACK❗✔️:          ENDIF
    // LAPACK❗✔️:          WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DGEHRD', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       ELSE IF( LQUERY ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 I = 1, ILO - 1
    // LAPACK❗✔️:          TAU( I ) = ZERO
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       DO 20 I = MAX( 1, IHI ), N - 1
    // LAPACK❗✔️:          TAU( I ) = ZERO
    // LAPACK❗✔️:    20 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NH.LE.1 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = 1
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Determine the block size
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NB = MIN( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
    // LAPACK❗✔️:       NBMIN = 2
    // LAPACK❗✔️:       IF( NB.GT.1 .AND. NB.LT.NH ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Determine when to cross over from blocked to unblocked code
    // LAPACK❗✔️: *        (last block is always handled by unblocked code)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NX = MAX( NB, ILAENV( 3, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
    // LAPACK❗✔️:          IF( NX.LT.NH ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Determine if workspace is large enough for blocked code
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( LWORK.LT.LWKOPT ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Not enough workspace to use optimal NB:  determine the
    // LAPACK❗✔️: *              minimum value of NB, and reduce NB or force use of
    // LAPACK❗✔️: *              unblocked code
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N, ILO, IHI,
    // LAPACK❗✔️:      $                 -1 ) )
    // LAPACK❗✔️:                IF( LWORK.GE.(N*NBMIN + TSIZE) ) THEN
    // LAPACK❗✔️:                   NB = (LWORK-TSIZE) / N
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   NB = 1
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       LDWORK = N
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Use unblocked code below
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          I = ILO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Use blocked code
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IWT = 1 + N*NB
    // LAPACK❗✔️:          DO 40 I = ILO, IHI - 1 - NX, NB
    // LAPACK❗✔️:             IB = MIN( NB, IHI-I )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Reduce columns i:i+ib-1 to Hessenberg form, returning the
    // LAPACK❗✔️: *           matrices V and T of the block reflector H = I - V*T*V**T
    // LAPACK❗✔️: *           which performs the reduction, and also the matrix Y = A*V*T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ),
    // LAPACK❗✔️:      $                   WORK( IWT ), LDT, WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
    // LAPACK❗✔️: *           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
    // LAPACK❗✔️: *           to 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             EI = A( I+IB, I+IB-1 )
    // LAPACK❗✔️:             A( I+IB, I+IB-1 ) = ONE
    // LAPACK❗✔️:             CALL DGEMM( 'No transpose', 'Transpose',
    // LAPACK❗✔️:      $                  IHI, IHI-I-IB+1,
    // LAPACK❗✔️:      $                  IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE,
    // LAPACK❗✔️:      $                  A( 1, I+IB ), LDA )
    // LAPACK❗✔️:             A( I+IB, I+IB-1 ) = EI
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
    // LAPACK❗✔️: *           right
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DTRMM( 'Right', 'Lower', 'Transpose',
    // LAPACK❗✔️:      $                  'Unit', I, IB-1,
    // LAPACK❗✔️:      $                  ONE, A( I+1, I ), LDA, WORK, LDWORK )
    // LAPACK❗✔️:             DO 30 J = 0, IB-2
    // LAPACK❗✔️:                CALL DAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1,
    // LAPACK❗✔️:      $                     A( 1, I+J+1 ), 1 )
    // LAPACK❗✔️:    30       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
    // LAPACK❗✔️: *           left
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFB( 'Left', 'Transpose', 'Forward',
    // LAPACK❗✔️:      $                   'Columnwise',
    // LAPACK❗✔️:      $                   IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA,
    // LAPACK❗✔️:      $                   WORK( IWT ), LDT, A( I+1, I+IB ), LDA,
    // LAPACK❗✔️:      $                   WORK, LDWORK )
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Use unblocked code to reduce the rest of the matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = LWKOPT
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DGEHRD
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    const LDT: usize = 65;
    const TSIZE: usize = LDT * 64;
    tau[..ilo.saturating_sub(1)].fill(0.0);
    if ihi >= 1 {
        tau[ihi - 1..n.saturating_sub(1)].fill(0.0);
    }
    let nh = ihi + 1 - ilo;
    if nh <= 1 {
        return Ok(());
    }
    let mut nb = 32usize;
    let mut nbmin = 2;
    let mut nx = 0;
    if nb > 1 && nb < nh {
        nx = 128;
        if nx < nh && work.len() < n * nb + TSIZE {
            nbmin = 2;
            if work.len() >= n * nbmin + TSIZE {
                nb = (work.len() - TSIZE) / n;
            } else {
                nb = 1;
            }
        }
    }
    let mut i = ilo;
    if nb >= nbmin && nb < nh {
        let iwt = n * nb;
        while i + nx < ihi {
            let ib = nb.min(ihi - i);
            {
                let (y, t) = work.split_at_mut(iwt);
                lahr2(
                    ihi,
                    i,
                    ib,
                    &mut a[(i - 1) * lda..],
                    lda,
                    &mut tau[i - 1..],
                    t,
                    LDT,
                    y,
                    n,
                )?;
            }
            let position = ix(i + ib, i + ib - 1, lda);
            let ei = a[position];
            a[position] = 1.0;
            blas3::gemm_into_within(
                NoTrans,
                Trans,
                ihi,
                ihi - i - ib + 1,
                ib,
                -1.0,
                a,
                blas3::MatrixSource::Separate(&work[..iwt]),
                n,
                blas3::MatrixSource::Within(ix(i + ib, i, lda)),
                lda,
                1.0,
                ix(1, i + ib, lda),
                lda,
            )?;
            a[position] = ei;
            blas3::trmm(
                Right,
                Lower,
                Trans,
                Unit,
                i,
                ib - 1,
                1.0,
                &a[ix(i + 1, i, lda)..],
                lda,
                &mut work[..iwt],
                n,
            )?;
            for j in 0..ib - 1 {
                blas::axpy(
                    i,
                    -1.0,
                    &work[n * j..],
                    1,
                    &mut a[ix(1, i + j + 1, lda)..],
                    1,
                );
            }
            {
                let (previous, current) = a.split_at_mut((i + ib - 1) * lda);
                let (scratch, t) = work.split_at_mut(iwt);
                larfb(
                    true,
                    true,
                    ihi - i,
                    n - i - ib + 1,
                    ib,
                    &previous[ix(i + 1, i, lda)..],
                    lda,
                    t,
                    LDT,
                    &mut current[i..],
                    lda,
                    scratch,
                    n,
                )?;
            }
            i += nb;
        }
    }
    gehd2(n, i, ihi, a, lda, tau, work);
    Ok(())
}

/// Fixed F/C DLARFT dispatch reached through DORMQR's NB<=32 blocks.
fn larft_qr(n: usize, k: usize, v: &[f64], ldv: usize, tau: &[f64], t: &mut [f64], ldt: usize) {
    // LAPACK❗✔️:       IF(N.EQ.0.OR.K.EQ.0) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF(N.EQ.1.OR.K.EQ.1) THEN
    // LAPACK❗✔️:          T(1,1) = TAU(1)
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       NX = ILAENV(3, "DLARFT", DIRECT // STOREV, N, K, -1, -1)
    // LAPACK❗✔️:       IF(K.LT.NX) THEN
    // LAPACK❗✔️:          CALL DLARFT_LVL2(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT)
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // ILAENV: DORMQR NB=32, DLARFT NX=64. Its private callers can never
    // reach the recursive K>=64 branch, irrespective of molecular dimension.
    if n == 0 || k == 0 {
        return;
    }
    if n == 1 || k == 1 {
        t[0] = tau[0];
        return;
    }
    larft_level2(n, k, v, ldv, tau, t, ldt);
}

#[allow(clippy::too_many_arguments)]
fn larft_level2(n: usize, k: usize, v: &[f64], ldv: usize, tau: &[f64], t: &mut [f64], ldt: usize) {
    // LAPACK❗✔️:       SUBROUTINE DLARFT_LVL2( DIRECT, STOREV, N, K, V, LDV, TAU,
    // LAPACK❗✔️:      $            T, LDT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          DIRECT, STOREV
    // LAPACK❗✔️:       INTEGER            K, LDT, LDV, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE, ZERO
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I, J, PREVLASTV, LASTV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DGEMV, DTRMV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( DIRECT, 'F' ) ) THEN
    // LAPACK❗✔️:          PREVLASTV = N
    // LAPACK❗✔️:          DO I = 1, K
    // LAPACK❗✔️:             PREVLASTV = MAX( I, PREVLASTV )
    // LAPACK❗✔️:             IF( TAU( I ).EQ.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              H(i)  =  I
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO J = 1, I
    // LAPACK❗✔️:                   T( J, I ) = ZERO
    // LAPACK❗✔️:                END DO
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              general case
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( LSAME( STOREV, 'C' ) ) THEN
    // LAPACK❗✔️: *                 Skip any trailing zeros.
    // LAPACK❗✔️:                   DO LASTV = N, I+1, -1
    // LAPACK❗✔️:                      IF( V( LASTV, I ).NE.ZERO ) EXIT
    // LAPACK❗✔️:                   END DO
    // LAPACK❗✔️:                   DO J = 1, I-1
    // LAPACK❗✔️:                      T( J, I ) = -TAU( I ) * V( I , J )
    // LAPACK❗✔️:                   END DO
    // LAPACK❗✔️:                   J = MIN( LASTV, PREVLASTV )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ),
    // LAPACK❗✔️:      $                        V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE,
    // LAPACK❗✔️:      $                        T( 1, I ), 1 )
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *                 Skip any trailing zeros.
    // LAPACK❗✔️:                   DO LASTV = N, I+1, -1
    // LAPACK❗✔️:                      IF( V( I, LASTV ).NE.ZERO ) EXIT
    // LAPACK❗✔️:                   END DO
    // LAPACK❗✔️:                   DO J = 1, I-1
    // LAPACK❗✔️:                      T( J, I ) = -TAU( I ) * V( J , I )
    // LAPACK❗✔️:                   END DO
    // LAPACK❗✔️:                   J = MIN( LASTV, PREVLASTV )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ),
    // LAPACK❗✔️:      $                        V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE,
    // LAPACK❗✔️:      $                        T( 1, I ), 1 )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1,
    // LAPACK❗✔️:      $                     T,
    // LAPACK❗✔️:      $                     LDT, T( 1, I ), 1 )
    // LAPACK❗✔️:                T( I, I ) = TAU( I )
    // LAPACK❗✔️:                IF( I.GT.1 ) THEN
    // LAPACK❗✔️:                   PREVLASTV = MAX( PREVLASTV, LASTV )
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   PREVLASTV = LASTV
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          PREVLASTV = 1
    // LAPACK❗✔️:          DO I = K, 1, -1
    // LAPACK❗✔️:             IF( TAU( I ).EQ.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              H(i)  =  I
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO J = I, K
    // LAPACK❗✔️:                   T( J, I ) = ZERO
    // LAPACK❗✔️:                END DO
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              general case
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( I.LT.K ) THEN
    // LAPACK❗✔️:                   IF( LSAME( STOREV, 'C' ) ) THEN
    // LAPACK❗✔️: *                    Skip any leading zeros.
    // LAPACK❗✔️:                      DO LASTV = 1, I-1
    // LAPACK❗✔️:                         IF( V( LASTV, I ).NE.ZERO ) EXIT
    // LAPACK❗✔️:                      END DO
    // LAPACK❗✔️:                      DO J = I+1, K
    // LAPACK❗✔️:                         T( J, I ) = -TAU( I ) * V( N-K+I , J )
    // LAPACK❗✔️:                      END DO
    // LAPACK❗✔️:                      J = MAX( LASTV, PREVLASTV )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      CALL DGEMV( 'Transpose', N-K+I-J, K-I,
    // LAPACK❗✔️:      $                           -TAU( I ),
    // LAPACK❗✔️:      $                           V( J, I+1 ), LDV, V( J, I ), 1, ONE,
    // LAPACK❗✔️:      $                           T( I+1, I ), 1 )
    // LAPACK❗✔️:                   ELSE
    // LAPACK❗✔️: *                    Skip any leading zeros.
    // LAPACK❗✔️:                      DO LASTV = 1, I-1
    // LAPACK❗✔️:                         IF( V( I, LASTV ).NE.ZERO ) EXIT
    // LAPACK❗✔️:                      END DO
    // LAPACK❗✔️:                      DO J = I+1, K
    // LAPACK❗✔️:                         T( J, I ) = -TAU( I ) * V( J, N-K+I )
    // LAPACK❗✔️:                      END DO
    // LAPACK❗✔️:                      J = MAX( LASTV, PREVLASTV )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      CALL DGEMV( 'No transpose', K-I, N-K+I-J,
    // LAPACK❗✔️:      $                    -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV,
    // LAPACK❗✔️:      $                    ONE, T( I+1, I ), 1 )
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DTRMV( 'Lower', 'No transpose', 'Non-unit',
    // LAPACK❗✔️:      $                        K-I,
    // LAPACK❗✔️:      $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
    // LAPACK❗✔️:                   IF( I.GT.1 ) THEN
    // LAPACK❗✔️:                      PREVLASTV = MIN( PREVLASTV, LASTV )
    // LAPACK❗✔️:                   ELSE
    // LAPACK❗✔️:                      PREVLASTV = LASTV
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:                T( I, I ) = TAU( I )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLARFT_LVL2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if n == 0 {
        return;
    }
    let mut prevlastv = n;
    for i in 1..=k {
        prevlastv = i.max(prevlastv);
        if tau[i - 1] == 0.0 {
            for j in 1..=i {
                t[ix(j, i, ldt)] = 0.0;
            }
        } else {
            let mut lastv = i;
            for row in (i + 1..=n).rev() {
                if v[ix(row, i, ldv)] != 0.0 {
                    lastv = row;
                    break;
                }
            }
            for j in 1..i {
                t[ix(j, i, ldt)] = -tau[i - 1] * v[ix(i, j, ldv)];
            }
            let j = lastv.min(prevlastv);
            if i > 1 {
                blas::gemv(
                    true,
                    j - i,
                    i - 1,
                    -tau[i - 1],
                    &v[ix(i + 1, 1, ldv)..],
                    ldv,
                    &v[ix(i + 1, i, ldv)..],
                    1,
                    1.0,
                    &mut t[(i - 1) * ldt..],
                    1,
                );
                let (previous, current) = t.split_at_mut((i - 1) * ldt);
                blas::trmv(false, false, false, i - 1, previous, ldt, current, 1);
            }
            t[ix(i, i, ldt)] = tau[i - 1];
            prevlastv = if i > 1 { prevlastv.max(lastv) } else { lastv };
        }
    }
}

/// SIDE=Right, TRANS=NoTranspose, as required by both AED routines.
#[allow(clippy::too_many_arguments)]
fn orm2r(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    tau: &[f64],
    c: &mut [f64],
    ldc: usize,
    work: &mut [f64],
) {
    // LAPACK❗✔️:       SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
    // LAPACK❗✔️:      $                   WORK, INFO )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          SIDE, TRANS
    // LAPACK❗✔️:       INTEGER            INFO, K, LDA, LDC, M, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            LEFT, NOTRAN
    // LAPACK❗✔️:       INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL            XERBLA, DLARF1F
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input arguments
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       LEFT = LSAME( SIDE, 'L' )
    // LAPACK❗✔️:       NOTRAN = LSAME( TRANS, 'N' )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     NQ is the order of Q
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LEFT ) THEN
    // LAPACK❗✔️:          NQ = M
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          NQ = N
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( M.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
    // LAPACK❗✔️:          INFO = -7
    // LAPACK❗✔️:       ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
    // LAPACK❗✔️:          INFO = -10
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DORM2R', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
    // LAPACK❗✔️:      $     THEN
    // LAPACK❗✔️:          I1 = 1
    // LAPACK❗✔️:          I2 = K
    // LAPACK❗✔️:          I3 = 1
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          I1 = K
    // LAPACK❗✔️:          I2 = 1
    // LAPACK❗✔️:          I3 = -1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LEFT ) THEN
    // LAPACK❗✔️:          NI = N
    // LAPACK❗✔️:          JC = 1
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          MI = M
    // LAPACK❗✔️:          IC = 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 10 I = I1, I2, I3
    // LAPACK❗✔️:          IF( LEFT ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           H(i) is applied to C(i:m,1:n)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             MI = M - I + 1
    // LAPACK❗✔️:             IC = I
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           H(i) is applied to C(1:m,i:n)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NI = N - I + 1
    // LAPACK❗✔️:             JC = I
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Apply H(i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARF1F( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC,
    // LAPACK❗✔️:      $               JC ),
    // LAPACK❗✔️:      $               LDC, WORK )
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DORM2R
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if m == 0 || n == 0 || k == 0 {
        return;
    }
    for i in 1..=k {
        larf1f(
            false,
            m,
            n - i + 1,
            &a[ix(i, i, lda)..],
            tau[i - 1],
            &mut c[(i - 1) * ldc..],
            ldc,
            work,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn ormqr(
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    tau: &[f64],
    c: &mut [f64],
    ldc: usize,
    work: &mut [f64],
) -> Result<(), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
    // LAPACK❗✔️:      $                   WORK, LWORK, INFO )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          SIDE, TRANS
    // LAPACK❗✔️:       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       INTEGER            NBMAX, LDT, TSIZE
    // LAPACK❗✔️:       PARAMETER          ( NBMAX = 64, LDT = NBMAX+1,
    // LAPACK❗✔️:      $                     TSIZE = LDT*NBMAX )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            LEFT, LQUERY, NOTRAN
    // LAPACK❗✔️:       INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,
    // LAPACK❗✔️:      $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           LSAME, ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input arguments
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       LEFT = LSAME( SIDE, 'L' )
    // LAPACK❗✔️:       NOTRAN = LSAME( TRANS, 'N' )
    // LAPACK❗✔️:       LQUERY = ( LWORK.EQ.-1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     NQ is the order of Q and NW is the minimum dimension of WORK
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LEFT ) THEN
    // LAPACK❗✔️:          NQ = M
    // LAPACK❗✔️:          NW = MAX( 1, N )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          NQ = N
    // LAPACK❗✔️:          NW = MAX( 1, M )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( M.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
    // LAPACK❗✔️:          INFO = -7
    // LAPACK❗✔️:       ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
    // LAPACK❗✔️:          INFO = -10
    // LAPACK❗✔️:       ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:          INFO = -12
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute the workspace requirements
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N,
    // LAPACK❗✔️:      $             K,
    // LAPACK❗✔️:      $        -1 ) )
    // LAPACK❗✔️:           LWKOPT = NW*NB + LDT*NB
    // LAPACK❗✔️:          WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DORMQR', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       ELSE IF( LQUERY ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = 1
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NBMIN = 2
    // LAPACK❗✔️:       LDWORK = NW
    // LAPACK❗✔️:       IF( NB.GT.1 .AND. NB.LT.K ) THEN
    // LAPACK❗✔️:          IF( LWORK.LT.LWKOPT ) THEN
    // LAPACK❗✔️:              NB = LWORK / (LDWORK + LDT)
    // LAPACK❗✔️:             NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N,
    // LAPACK❗✔️:      $                   K,
    // LAPACK❗✔️:      $              -1 ) )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Use unblocked code
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
    // LAPACK❗✔️:      $                WORK,
    // LAPACK❗✔️:      $                IINFO )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Use blocked code
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IWT = 1 + NW*NB
    // LAPACK❗✔️:          IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
    // LAPACK❗✔️:      $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
    // LAPACK❗✔️:             I1 = 1
    // LAPACK❗✔️:             I2 = K
    // LAPACK❗✔️:             I3 = NB
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             I1 = ( ( K-1 ) / NB )*NB + 1
    // LAPACK❗✔️:             I2 = 1
    // LAPACK❗✔️:             I3 = -NB
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LEFT ) THEN
    // LAPACK❗✔️:             NI = N
    // LAPACK❗✔️:             JC = 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             MI = M
    // LAPACK❗✔️:             IC = 1
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 10 I = I1, I2, I3
    // LAPACK❗✔️:             IB = MIN( NB, K-I+1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Form the triangular factor of the block reflector
    // LAPACK❗✔️: *           H = H(i) H(i+1) . . . H(i+ib-1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I,
    // LAPACK❗✔️:      $                   I ),
    // LAPACK❗✔️:      $                   LDA, TAU( I ), WORK( IWT ), LDT )
    // LAPACK❗✔️:             IF( LEFT ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              H or H**T is applied to C(i:m,1:n)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                MI = M - I + 1
    // LAPACK❗✔️:                IC = I
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              H or H**T is applied to C(1:m,i:n)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NI = N - I + 1
    // LAPACK❗✔️:                JC = I
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply H or H**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI,
    // LAPACK❗✔️:      $                   NI,
    // LAPACK❗✔️:      $                   IB, A( I, I ), LDA, WORK( IWT ), LDT,
    // LAPACK❗✔️:      $                   C( IC, JC ), LDC, WORK, LDWORK )
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DORMQR
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if m == 0 || n == 0 || k == 0 {
        return Ok(());
    }
    const LDT: usize = 65;
    let mut nb = 32usize;
    if nb > 1 && nb < k && work.len() < m * nb + LDT * nb {
        nb = work.len() / (m + LDT);
    }
    if nb < 2 || nb >= k {
        orm2r(m, n, k, a, lda, tau, c, ldc, work);
    } else {
        let iwt = m * nb;
        for i in (1..=k).step_by(nb) {
            let ib = nb.min(k - i + 1);
            let (scratch, t) = work.split_at_mut(iwt);
            larft_qr(
                n - i + 1,
                ib,
                &a[ix(i, i, lda)..],
                lda,
                &tau[i - 1..],
                t,
                LDT,
            );
            larfb(
                false,
                false,
                m,
                n - i + 1,
                ib,
                &a[ix(i, i, lda)..],
                lda,
                t,
                LDT,
                &mut c[(i - 1) * ldc..],
                ldc,
                scratch,
                m,
            )?;
        }
    }
    Ok(())
}

pub(super) fn ormhr_work(m: usize) -> usize {
    m * 32
}

#[allow(clippy::too_many_arguments)]
pub(super) fn ormhr(
    m: usize,
    n: usize,
    ilo: usize,
    ihi: usize,
    a: &[f64],
    lda: usize,
    tau: &[f64],
    c: &mut [f64],
    ldc: usize,
    work: &mut [f64],
) -> Result<(), GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
    // LAPACK❗✔️:      $                   LDC, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          SIDE, TRANS
    // LAPACK❗✔️:       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            LEFT, LQUERY
    // LAPACK❗✔️:       INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           LSAME, ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DORMQR, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input arguments
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       NH = IHI - ILO
    // LAPACK❗✔️:       LEFT = LSAME( SIDE, 'L' )
    // LAPACK❗✔️:       LQUERY = ( LWORK.EQ.-1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     NQ is the order of Q and NW is the minimum dimension of WORK
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LEFT ) THEN
    // LAPACK❗✔️:          NQ = M
    // LAPACK❗✔️:          NW = MAX( 1, N )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          NQ = N
    // LAPACK❗✔️:          NW = MAX( 1, M )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) )
    // LAPACK❗✔️:      $          THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( M.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, NQ ) ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( IHI.LT.MIN( ILO, NQ ) .OR. IHI.GT.NQ ) THEN
    // LAPACK❗✔️:          INFO = -6
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
    // LAPACK❗✔️:          INFO = -8
    // LAPACK❗✔️:       ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
    // LAPACK❗✔️:          INFO = -11
    // LAPACK❗✔️:       ELSE IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:          INFO = -13
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️:          IF( LEFT ) THEN
    // LAPACK❗✔️:             NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, NH, N, NH, -1 )
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, NH, NH, -1 )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          LWKOPT = NW*NB
    // LAPACK❗✔️:          WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DORMHR', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       ELSE IF( LQUERY ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( M.EQ.0 .OR. N.EQ.0 .OR. NH.EQ.0 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = 1
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LEFT ) THEN
    // LAPACK❗✔️:          MI = NH
    // LAPACK❗✔️:          NI = N
    // LAPACK❗✔️:          I1 = ILO + 1
    // LAPACK❗✔️:          I2 = 1
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          MI = M
    // LAPACK❗✔️:          NI = NH
    // LAPACK❗✔️:          I1 = 1
    // LAPACK❗✔️:          I2 = ILO + 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DORMQR( SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA,
    // LAPACK❗✔️:      $             TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DORMHR
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    let nh = ihi - ilo;
    if m == 0 || n == 0 || nh == 0 {
        return Ok(());
    }
    ormqr(
        m,
        nh,
        nh,
        &a[ix(ilo + 1, ilo, lda)..],
        lda,
        &tau[ilo - 1..],
        &mut c[ilo * ldc..],
        ldc,
        work,
    )
}
