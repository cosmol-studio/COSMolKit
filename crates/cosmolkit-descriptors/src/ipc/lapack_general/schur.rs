//! Source double-shift and multishift Hessenberg eigenvalue iteration.
use super::super::lapack_tridiagonal::dlamch;
use super::{ix, lanv2};

#[allow(clippy::too_many_arguments)]
pub(super) fn lahqr(
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
) -> usize {
    // LAPACK❗✔️:       SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
    // LAPACK❗✔️:      $                   ILOZ, IHIZ, Z, LDZ, INFO )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
    // LAPACK❗✔️:       LOGICAL            WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =========================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE, TWO
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   DAT1, DAT2
    // LAPACK❗✔️:       PARAMETER          ( DAT1 = 3.0d0 / 4.0d0, DAT2 = -0.4375d0 )
    // LAPACK❗✔️:       INTEGER            KEXSH
    // LAPACK❗✔️:       PARAMETER          ( KEXSH = 10 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S,
    // LAPACK❗✔️:      $                   H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX,
    // LAPACK❗✔️:      $                   SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST,
    // LAPACK❗✔️:      $                   ULP, V2, V3
    // LAPACK❗✔️:       INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ,
    // LAPACK❗✔️:      $                   KDEFL
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       DOUBLE PRECISION   V( 3 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️:       EXTERNAL           DLAMCH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DCOPY, DLABAD, DLANV2, DLARFG, DROT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️:       IF( ILO.EQ.IHI ) THEN
    // LAPACK❗✔️:          WR( ILO ) = H( ILO, ILO )
    // LAPACK❗✔️:          WI( ILO ) = ZERO
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== clear out the trash ====
    // LAPACK❗✔️:       DO 10 J = ILO, IHI - 3
    // LAPACK❗✔️:          H( J+2, J ) = ZERO
    // LAPACK❗✔️:          H( J+3, J ) = ZERO
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       IF( ILO.LE.IHI-2 )
    // LAPACK❗✔️:      $   H( IHI, IHI-2 ) = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NH = IHI - ILO + 1
    // LAPACK❗✔️:       NZ = IHIZ - ILOZ + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Set machine-dependent constants for the stopping criterion.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'SAFE MINIMUM' )
    // LAPACK❗✔️:       SAFMAX = ONE / SAFMIN
    // LAPACK❗✔️:       CALL DLABAD( SAFMIN, SAFMAX )
    // LAPACK❗✔️:       ULP = DLAMCH( 'PRECISION' )
    // LAPACK❗✔️:       SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     I1 and I2 are the indices of the first row and last column of H
    // LAPACK❗✔️: *     to which transformations must be applied. If eigenvalues only are
    // LAPACK❗✔️: *     being computed, I1 and I2 are set inside the main loop.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( WANTT ) THEN
    // LAPACK❗✔️:          I1 = 1
    // LAPACK❗✔️:          I2 = N
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ITMAX is the total number of QR iterations allowed.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ITMAX = 30 * MAX( 10, NH )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     KDEFL counts the number of iterations since a deflation
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       KDEFL = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     The main loop begins here. I is the loop index and decreases from
    // LAPACK❗✔️: *     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
    // LAPACK❗✔️: *     with the active submatrix in rows and columns L to I.
    // LAPACK❗✔️: *     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
    // LAPACK❗✔️: *     H(L,L-1) is negligible so that the matrix splits.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       I = IHI
    // LAPACK❗✔️:    20 CONTINUE
    // LAPACK❗✔️:       L = ILO
    // LAPACK❗✔️:       IF( I.LT.ILO )
    // LAPACK❗✔️:      $   GO TO 160
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Perform QR iterations on rows and columns ILO to I until a
    // LAPACK❗✔️: *     submatrix of order 1 or 2 splits off at the bottom because a
    // LAPACK❗✔️: *     subdiagonal element has become negligible.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 140 ITS = 0, ITMAX
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Look for a single small subdiagonal element.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 30 K = I, L + 1, -1
    // LAPACK❗✔️:             IF( ABS( H( K, K-1 ) ).LE.SMLNUM )
    // LAPACK❗✔️:      $         GO TO 40
    // LAPACK❗✔️:             TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
    // LAPACK❗✔️:             IF( TST.EQ.ZERO ) THEN
    // LAPACK❗✔️:                IF( K-2.GE.ILO )
    // LAPACK❗✔️:      $            TST = TST + ABS( H( K-1, K-2 ) )
    // LAPACK❗✔️:                IF( K+1.LE.IHI )
    // LAPACK❗✔️:      $            TST = TST + ABS( H( K+1, K ) )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *           ==== The following is a conservative small subdiagonal
    // LAPACK❗✔️: *           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
    // LAPACK❗✔️: *           .    1997). It has better mathematical foundation and
    // LAPACK❗✔️: *           .    improves accuracy in some cases.  ====
    // LAPACK❗✔️:             IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
    // LAPACK❗✔️:                AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
    // LAPACK❗✔️:                BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
    // LAPACK❗✔️:                AA = MAX( ABS( H( K, K ) ),
    // LAPACK❗✔️:      $              ABS( H( K-1, K-1 )-H( K, K ) ) )
    // LAPACK❗✔️:                BB = MIN( ABS( H( K, K ) ),
    // LAPACK❗✔️:      $              ABS( H( K-1, K-1 )-H( K, K ) ) )
    // LAPACK❗✔️:                S = AA + AB
    // LAPACK❗✔️:                IF( BA*( AB / S ).LE.MAX( SMLNUM,
    // LAPACK❗✔️:      $             ULP*( BB*( AA / S ) ) ) )GO TO 40
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:    30    CONTINUE
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️:          L = K
    // LAPACK❗✔️:          IF( L.GT.ILO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           H(L,L-1) is negligible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             H( L, L-1 ) = ZERO
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Exit from loop if a submatrix of order 1 or 2 has split off.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( L.GE.I-1 )
    // LAPACK❗✔️:      $      GO TO 150
    // LAPACK❗✔️:          KDEFL = KDEFL + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Now the active submatrix is in rows and columns L to I. If
    // LAPACK❗✔️: *        eigenvalues only are being computed, only the active submatrix
    // LAPACK❗✔️: *        need be transformed.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( .NOT.WANTT ) THEN
    // LAPACK❗✔️:             I1 = L
    // LAPACK❗✔️:             I2 = I
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( MOD(KDEFL,2*KEXSH).EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Exceptional shift.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
    // LAPACK❗✔️:             H11 = DAT1*S + H( I, I )
    // LAPACK❗✔️:             H12 = DAT2*S
    // LAPACK❗✔️:             H21 = S
    // LAPACK❗✔️:             H22 = H11
    // LAPACK❗✔️:          ELSE IF( MOD(KDEFL,KEXSH).EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Exceptional shift.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             S = ABS( H( L+1, L ) ) + ABS( H( L+2, L+1 ) )
    // LAPACK❗✔️:             H11 = DAT1*S + H( L, L )
    // LAPACK❗✔️:             H12 = DAT2*S
    // LAPACK❗✔️:             H21 = S
    // LAPACK❗✔️:             H22 = H11
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Prepare to use Francis' double shift
    // LAPACK❗✔️: *           (i.e. 2nd degree generalized Rayleigh quotient)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             H11 = H( I-1, I-1 )
    // LAPACK❗✔️:             H21 = H( I, I-1 )
    // LAPACK❗✔️:             H12 = H( I-1, I )
    // LAPACK❗✔️:             H22 = H( I, I )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
    // LAPACK❗✔️:          IF( S.EQ.ZERO ) THEN
    // LAPACK❗✔️:             RT1R = ZERO
    // LAPACK❗✔️:             RT1I = ZERO
    // LAPACK❗✔️:             RT2R = ZERO
    // LAPACK❗✔️:             RT2I = ZERO
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             H11 = H11 / S
    // LAPACK❗✔️:             H21 = H21 / S
    // LAPACK❗✔️:             H12 = H12 / S
    // LAPACK❗✔️:             H22 = H22 / S
    // LAPACK❗✔️:             TR = ( H11+H22 ) / TWO
    // LAPACK❗✔️:             DET = ( H11-TR )*( H22-TR ) - H12*H21
    // LAPACK❗✔️:             RTDISC = SQRT( ABS( DET ) )
    // LAPACK❗✔️:             IF( DET.GE.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== complex conjugate shifts ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                RT1R = TR*S
    // LAPACK❗✔️:                RT2R = RT1R
    // LAPACK❗✔️:                RT1I = RTDISC*S
    // LAPACK❗✔️:                RT2I = -RT1I
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== real shifts (use only one of them)  ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                RT1R = TR + RTDISC
    // LAPACK❗✔️:                RT2R = TR - RTDISC
    // LAPACK❗✔️:                IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
    // LAPACK❗✔️:                   RT1R = RT1R*S
    // LAPACK❗✔️:                   RT2R = RT1R
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   RT2R = RT2R*S
    // LAPACK❗✔️:                   RT1R = RT2R
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:                RT1I = ZERO
    // LAPACK❗✔️:                RT2I = ZERO
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Look for two consecutive small subdiagonal elements.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 50 M = I - 2, L, -1
    // LAPACK❗✔️: *           Determine the effect of starting the double-shift QR
    // LAPACK❗✔️: *           iteration at row M, and see if this would make H(M,M-1)
    // LAPACK❗✔️: *           negligible.  (The following uses scaling to avoid
    // LAPACK❗✔️: *           overflows and most underflows.)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             H21S = H( M+1, M )
    // LAPACK❗✔️:             S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
    // LAPACK❗✔️:             H21S = H( M+1, M ) / S
    // LAPACK❗✔️:             V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )*
    // LAPACK❗✔️:      $               ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
    // LAPACK❗✔️:             V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
    // LAPACK❗✔️:             V( 3 ) = H21S*H( M+2, M+1 )
    // LAPACK❗✔️:             S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
    // LAPACK❗✔️:             V( 1 ) = V( 1 ) / S
    // LAPACK❗✔️:             V( 2 ) = V( 2 ) / S
    // LAPACK❗✔️:             V( 3 ) = V( 3 ) / S
    // LAPACK❗✔️:             IF( M.EQ.L )
    // LAPACK❗✔️:      $         GO TO 60
    // LAPACK❗✔️:             IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE.
    // LAPACK❗✔️:      $          ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M,
    // LAPACK❗✔️:      $          M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
    // LAPACK❗✔️:    50    CONTINUE
    // LAPACK❗✔️:    60    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Double-shift QR step
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 130 K = M, I - 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           The first iteration of this loop determines a reflection G
    // LAPACK❗✔️: *           from the vector V and applies it from left and right to H,
    // LAPACK❗✔️: *           thus creating a nonzero bulge below the subdiagonal.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Each subsequent iteration determines a reflection G to
    // LAPACK❗✔️: *           restore the Hessenberg form in the (K-1)th column, and thus
    // LAPACK❗✔️: *           chases the bulge one step toward the bottom of the active
    // LAPACK❗✔️: *           submatrix. NR is the order of G.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NR = MIN( 3, I-K+1 )
    // LAPACK❗✔️:             IF( K.GT.M )
    // LAPACK❗✔️:      $         CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
    // LAPACK❗✔️:             CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
    // LAPACK❗✔️:             IF( K.GT.M ) THEN
    // LAPACK❗✔️:                H( K, K-1 ) = V( 1 )
    // LAPACK❗✔️:                H( K+1, K-1 ) = ZERO
    // LAPACK❗✔️:                IF( K.LT.I-1 )
    // LAPACK❗✔️:      $            H( K+2, K-1 ) = ZERO
    // LAPACK❗✔️:             ELSE IF( M.GT.L ) THEN
    // LAPACK❗✔️: *               ==== Use the following instead of
    // LAPACK❗✔️: *               .    H( K, K-1 ) = -H( K, K-1 ) to
    // LAPACK❗✔️: *               .    avoid a bug when v(2) and v(3)
    // LAPACK❗✔️: *               .    underflow. ====
    // LAPACK❗✔️:                H( K, K-1 ) = H( K, K-1 )*( ONE-T1 )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             V2 = V( 2 )
    // LAPACK❗✔️:             T2 = T1*V2
    // LAPACK❗✔️:             IF( NR.EQ.3 ) THEN
    // LAPACK❗✔️:                V3 = V( 3 )
    // LAPACK❗✔️:                T3 = T1*V3
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply G from the left to transform the rows of the matrix
    // LAPACK❗✔️: *              in columns K to I2.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 70 J = K, I2
    // LAPACK❗✔️:                   SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
    // LAPACK❗✔️:                   H( K, J ) = H( K, J ) - SUM*T1
    // LAPACK❗✔️:                   H( K+1, J ) = H( K+1, J ) - SUM*T2
    // LAPACK❗✔️:                   H( K+2, J ) = H( K+2, J ) - SUM*T3
    // LAPACK❗✔️:    70          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply G from the right to transform the columns of the
    // LAPACK❗✔️: *              matrix in rows I1 to min(K+3,I).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 80 J = I1, MIN( K+3, I )
    // LAPACK❗✔️:                   SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
    // LAPACK❗✔️:                   H( J, K ) = H( J, K ) - SUM*T1
    // LAPACK❗✔️:                   H( J, K+1 ) = H( J, K+1 ) - SUM*T2
    // LAPACK❗✔️:                   H( J, K+2 ) = H( J, K+2 ) - SUM*T3
    // LAPACK❗✔️:    80          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( WANTZ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 Accumulate transformations in the matrix Z
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   DO 90 J = ILOZ, IHIZ
    // LAPACK❗✔️:                      SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
    // LAPACK❗✔️:                      Z( J, K ) = Z( J, K ) - SUM*T1
    // LAPACK❗✔️:                      Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
    // LAPACK❗✔️:                      Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
    // LAPACK❗✔️:    90             CONTINUE
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             ELSE IF( NR.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply G from the left to transform the rows of the matrix
    // LAPACK❗✔️: *              in columns K to I2.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 100 J = K, I2
    // LAPACK❗✔️:                   SUM = H( K, J ) + V2*H( K+1, J )
    // LAPACK❗✔️:                   H( K, J ) = H( K, J ) - SUM*T1
    // LAPACK❗✔️:                   H( K+1, J ) = H( K+1, J ) - SUM*T2
    // LAPACK❗✔️:   100          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply G from the right to transform the columns of the
    // LAPACK❗✔️: *              matrix in rows I1 to min(K+3,I).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                DO 110 J = I1, I
    // LAPACK❗✔️:                   SUM = H( J, K ) + V2*H( J, K+1 )
    // LAPACK❗✔️:                   H( J, K ) = H( J, K ) - SUM*T1
    // LAPACK❗✔️:                   H( J, K+1 ) = H( J, K+1 ) - SUM*T2
    // LAPACK❗✔️:   110          CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( WANTZ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 Accumulate transformations in the matrix Z
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   DO 120 J = ILOZ, IHIZ
    // LAPACK❗✔️:                      SUM = Z( J, K ) + V2*Z( J, K+1 )
    // LAPACK❗✔️:                      Z( J, K ) = Z( J, K ) - SUM*T1
    // LAPACK❗✔️:                      Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
    // LAPACK❗✔️:   120             CONTINUE
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:   130    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   140 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Failure to converge in remaining number of iterations
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = I
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   150 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( L.EQ.I ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        H(I,I-1) is negligible: one eigenvalue has converged.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          WR( I ) = H( I, I )
    // LAPACK❗✔️:          WI( I ) = ZERO
    // LAPACK❗✔️:       ELSE IF( L.EQ.I-1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Transform the 2-by-2 submatrix to standard Schur form,
    // LAPACK❗✔️: *        and compute and store the eigenvalues.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ),
    // LAPACK❗✔️:      $                H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ),
    // LAPACK❗✔️:      $                CS, SN )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTT ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply the transformation to the rest of H.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( I2.GT.I )
    // LAPACK❗✔️:      $         CALL DROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH,
    // LAPACK❗✔️:      $                    CS, SN )
    // LAPACK❗✔️:             CALL DROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( WANTZ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Apply the transformation to Z.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *     reset deflation counter
    // LAPACK❗✔️:       KDEFL = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     return to start of the main loop with new value of I.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       I = L - 1
    // LAPACK❗✔️:       GO TO 20
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   160 CONTINUE
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLAHQR
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if n == 0 {
        return 0;
    }
    if ilo == ihi {
        wr[ilo - 1] = h[ix(ilo, ilo, ldh)];
        wi[ilo - 1] = 0.0;
        return 0;
    }
    if ihi >= ilo + 3 {
        for j in ilo..=ihi - 3 {
            h[ix(j + 2, j, ldh)] = 0.0;
            h[ix(j + 3, j, ldh)] = 0.0;
        }
    }
    if ilo + 2 <= ihi {
        h[ix(ihi, ihi - 2, ldh)] = 0.0;
    }
    let nh = ihi - ilo + 1;
    let safmin = dlamch(b'S');
    let ulp = dlamch(b'P');
    // DLABAD leaves IEEE binary64's exponent range unchanged.
    let smlnum = safmin * (nh as f64 / ulp);
    let (mut i1, mut i2) = if wantt { (1, n) } else { (ilo, ihi) };
    let itmax = 30 * 10.max(nh);
    let mut kdefl = 0;
    let mut i = ihi;
    let mut v = [0.0; 3];
    while i >= ilo {
        let mut l = ilo;
        let mut converged = false;
        for _its in 0..=itmax {
            let mut found = l;
            for k in (l + 1..=i).rev() {
                if h[ix(k, k - 1, ldh)].abs() <= smlnum {
                    found = k;
                    break;
                }
                let mut tst = h[ix(k - 1, k - 1, ldh)].abs() + h[ix(k, k, ldh)].abs();
                if tst == 0.0 {
                    if k >= ilo + 2 {
                        tst += h[ix(k - 1, k - 2, ldh)].abs();
                    }
                    if k + 1 <= ihi {
                        tst += h[ix(k + 1, k, ldh)].abs();
                    }
                }
                if h[ix(k, k - 1, ldh)].abs() <= ulp * tst {
                    let ab = h[ix(k, k - 1, ldh)].abs().max(h[ix(k - 1, k, ldh)].abs());
                    let ba = h[ix(k, k - 1, ldh)].abs().min(h[ix(k - 1, k, ldh)].abs());
                    let aa = h[ix(k, k, ldh)]
                        .abs()
                        .max((h[ix(k - 1, k - 1, ldh)] - h[ix(k, k, ldh)]).abs());
                    let bb = h[ix(k, k, ldh)]
                        .abs()
                        .min((h[ix(k - 1, k - 1, ldh)] - h[ix(k, k, ldh)]).abs());
                    let s = aa + ab;
                    if ba * (ab / s) <= smlnum.max(ulp * (bb * (aa / s))) {
                        found = k;
                        break;
                    }
                }
            }
            l = found;
            if l > ilo {
                h[ix(l, l - 1, ldh)] = 0.0;
            }
            if l >= i - 1 {
                converged = true;
                break;
            }
            kdefl += 1;
            if !wantt {
                i1 = l;
                i2 = i;
            }
            let (mut h11, mut h12, mut h21, mut h22);
            if kdefl % 20 == 0 {
                let s = h[ix(i, i - 1, ldh)].abs() + h[ix(i - 1, i - 2, ldh)].abs();
                h11 = 0.75 * s + h[ix(i, i, ldh)];
                h12 = -0.4375 * s;
                h21 = s;
                h22 = h11;
            } else if kdefl % 10 == 0 {
                let s = h[ix(l + 1, l, ldh)].abs() + h[ix(l + 2, l + 1, ldh)].abs();
                h11 = 0.75 * s + h[ix(l, l, ldh)];
                h12 = -0.4375 * s;
                h21 = s;
                h22 = h11;
            } else {
                h11 = h[ix(i - 1, i - 1, ldh)];
                h21 = h[ix(i, i - 1, ldh)];
                h12 = h[ix(i - 1, i, ldh)];
                h22 = h[ix(i, i, ldh)];
            }
            let s = h11.abs() + h12.abs() + h21.abs() + h22.abs();
            let (rt1r, rt1i, rt2r, rt2i);
            if s == 0.0 {
                rt1r = 0.0;
                rt1i = 0.0;
                rt2r = 0.0;
                rt2i = 0.0;
            } else {
                h11 /= s;
                h21 /= s;
                h12 /= s;
                h22 /= s;
                let tr = (h11 + h22) / 2.0;
                let det = (h11 - tr) * (h22 - tr) - h12 * h21;
                let rtdisc = det.abs().sqrt();
                if det >= 0.0 {
                    rt1r = tr * s;
                    rt2r = rt1r;
                    rt1i = rtdisc * s;
                    rt2i = -rt1i;
                } else {
                    let r1 = tr + rtdisc;
                    let r2 = tr - rtdisc;
                    if (r1 - h22).abs() <= (r2 - h22).abs() {
                        rt1r = r1 * s;
                        rt2r = rt1r;
                    } else {
                        rt2r = r2 * s;
                        rt1r = rt2r;
                    }
                    rt1i = 0.0;
                    rt2i = 0.0;
                }
            }
            let mut m = l;
            for row in (l..=i - 2).rev() {
                m = row;
                let s = (h[ix(m, m, ldh)] - rt2r).abs() + rt2i.abs() + h[ix(m + 1, m, ldh)].abs();
                let h21s = h[ix(m + 1, m, ldh)] / s;
                v[0] = h21s * h[ix(m, m + 1, ldh)]
                    + (h[ix(m, m, ldh)] - rt1r) * ((h[ix(m, m, ldh)] - rt2r) / s)
                    - rt1i * (rt2i / s);
                v[1] = h21s * (h[ix(m, m, ldh)] + h[ix(m + 1, m + 1, ldh)] - rt1r - rt2r);
                v[2] = h21s * h[ix(m + 2, m + 1, ldh)];
                let s = v[0].abs() + v[1].abs() + v[2].abs();
                v[0] /= s;
                v[1] /= s;
                v[2] /= s;
                if m == l {
                    break;
                }
                if h[ix(m, m - 1, ldh)].abs() * (v[1].abs() + v[2].abs())
                    <= ulp
                        * v[0].abs()
                        * (h[ix(m - 1, m - 1, ldh)].abs()
                            + h[ix(m, m, ldh)].abs()
                            + h[ix(m + 1, m + 1, ldh)].abs())
                {
                    break;
                }
            }
            for k in m..i {
                let nr = 3.min(i - k + 1);
                if k > m {
                    for j in 0..nr {
                        v[j] = h[ix(k + j, k - 1, ldh)];
                    }
                }
                let t1 = super::super::lapack_symmetric::larfg(&mut v[..nr]);
                if k > m {
                    h[ix(k, k - 1, ldh)] = v[0];
                    h[ix(k + 1, k - 1, ldh)] = 0.0;
                    if k < i - 1 {
                        h[ix(k + 2, k - 1, ldh)] = 0.0;
                    }
                } else if m > l {
                    h[ix(k, k - 1, ldh)] *= 1.0 - t1;
                }
                let v2 = v[1];
                let t2 = t1 * v2;
                if nr == 3 {
                    let v3 = v[2];
                    let t3 = t1 * v3;
                    for j in k..=i2 {
                        let sum = h[ix(k, j, ldh)]
                            + v2 * h[ix(k + 1, j, ldh)]
                            + v3 * h[ix(k + 2, j, ldh)];
                        h[ix(k, j, ldh)] -= sum * t1;
                        h[ix(k + 1, j, ldh)] -= sum * t2;
                        h[ix(k + 2, j, ldh)] -= sum * t3;
                    }
                    for j in i1..=(k + 3).min(i) {
                        let sum = h[ix(j, k, ldh)]
                            + v2 * h[ix(j, k + 1, ldh)]
                            + v3 * h[ix(j, k + 2, ldh)];
                        h[ix(j, k, ldh)] -= sum * t1;
                        h[ix(j, k + 1, ldh)] -= sum * t2;
                        h[ix(j, k + 2, ldh)] -= sum * t3;
                    }
                    if wantz {
                        for j in iloz..=ihiz {
                            let sum = z[ix(j, k, ldz)]
                                + v2 * z[ix(j, k + 1, ldz)]
                                + v3 * z[ix(j, k + 2, ldz)];
                            z[ix(j, k, ldz)] -= sum * t1;
                            z[ix(j, k + 1, ldz)] -= sum * t2;
                            z[ix(j, k + 2, ldz)] -= sum * t3;
                        }
                    }
                } else {
                    for j in k..=i2 {
                        let sum = h[ix(k, j, ldh)] + v2 * h[ix(k + 1, j, ldh)];
                        h[ix(k, j, ldh)] -= sum * t1;
                        h[ix(k + 1, j, ldh)] -= sum * t2;
                    }
                    for j in i1..=i {
                        let sum = h[ix(j, k, ldh)] + v2 * h[ix(j, k + 1, ldh)];
                        h[ix(j, k, ldh)] -= sum * t1;
                        h[ix(j, k + 1, ldh)] -= sum * t2;
                    }
                    if wantz {
                        for j in iloz..=ihiz {
                            let sum = z[ix(j, k, ldz)] + v2 * z[ix(j, k + 1, ldz)];
                            z[ix(j, k, ldz)] -= sum * t1;
                            z[ix(j, k + 1, ldz)] -= sum * t2;
                        }
                    }
                }
            }
        }
        if !converged {
            return i;
        }
        if l == i {
            wr[i - 1] = h[ix(i, i, ldh)];
            wi[i - 1] = 0.0;
        } else {
            let p = lanv2(
                h[ix(i - 1, i - 1, ldh)],
                h[ix(i - 1, i, ldh)],
                h[ix(i, i - 1, ldh)],
                h[ix(i, i, ldh)],
            );
            h[ix(i - 1, i - 1, ldh)] = p.a;
            h[ix(i - 1, i, ldh)] = p.b;
            h[ix(i, i - 1, ldh)] = p.c;
            h[ix(i, i, ldh)] = p.d;
            wr[i - 2] = p.real[0];
            wr[i - 1] = p.real[1];
            wi[i - 2] = p.imag[0];
            wi[i - 1] = p.imag[1];
            if wantt {
                if i2 > i {
                    super::super::lapack_blas::rot_within(
                        i2 - i,
                        h,
                        ix(i - 1, i + 1, ldh),
                        ldh,
                        ix(i, i + 1, ldh),
                        ldh,
                        p.cs,
                        p.sn,
                    );
                }
                super::super::lapack_blas::rot_within(
                    i - i1 - 1,
                    h,
                    ix(i1, i - 1, ldh),
                    1,
                    ix(i1, i, ldh),
                    1,
                    p.cs,
                    p.sn,
                );
            }
            if wantz {
                super::super::lapack_blas::rot_within(
                    ihiz - iloz + 1,
                    z,
                    ix(iloz, i - 1, ldz),
                    1,
                    ix(iloz, i, ldz),
                    1,
                    p.cs,
                    p.sn,
                );
            }
        }
        kdefl = 0;
        i = l - 1;
    }
    0
}

/// DHSEQR JOB='E', COMPZ='N', with the original padded rare-failure retry.
pub(super) fn hseqr(
    n: usize,
    ilo: usize,
    ihi: usize,
    h: &mut [f64],
    ldh: usize,
    wr: &mut [f64],
    wi: &mut [f64],
    work: &mut [f64],
) -> Result<usize, super::GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
    // LAPACK❗✔️:      $                   LDZ, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
    // LAPACK❗✔️:       CHARACTER          COMPZ, JOB
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
    // LAPACK❗✔️:      $                   Z( LDZ, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Matrices of order NTINY or smaller must be processed by
    // LAPACK❗✔️: *     .    DLAHQR because of insufficient subdiagonal scratch space.
    // LAPACK❗✔️: *     .    (This is a hard limit.) ====
    // LAPACK❗✔️:       INTEGER            NTINY
    // LAPACK❗✔️:       PARAMETER          ( NTINY = 15 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== NL allocates some local workspace to help small matrices
    // LAPACK❗✔️: *     .    through a rare DLAHQR failure.  NL > NTINY = 15 is
    // LAPACK❗✔️: *     .    required and NL <= NMIN = ILAENV(ISPEC=12,...) is recom-
    // LAPACK❗✔️: *     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
    // LAPACK❗✔️: *     .    allows up to six simultaneous shifts and a 16-by-16
    // LAPACK❗✔️: *     .    deflation window.  ====
    // LAPACK❗✔️:       INTEGER            NL
    // LAPACK❗✔️:       PARAMETER          ( NL = 49 )
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       DOUBLE PRECISION   HL( NL, NL ), WORKL( NL )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I, KBOT, NMIN
    // LAPACK❗✔️:       LOGICAL            INITZ, LQUERY, WANTT, WANTZ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           ILAENV, LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLACPY, DLAHQR, DLAQR0, DLASET, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          DBLE, MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== Decode and check the input parameters. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WANTT = LSAME( JOB, 'S' )
    // LAPACK❗✔️:       INITZ = LSAME( COMPZ, 'I' )
    // LAPACK❗✔️:       WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
    // LAPACK❗✔️:       WORK( 1 ) = DBLE( MAX( 1, N ) )
    // LAPACK❗✔️:       LQUERY = LWORK.EQ.-1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -7
    // LAPACK❗✔️:       ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
    // LAPACK❗✔️:          INFO = -11
    // LAPACK❗✔️:       ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:          INFO = -13
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Quick return in case of invalid argument. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL XERBLA( 'DHSEQR', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( N.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Quick return in case N = 0; nothing to do. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( LQUERY ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Quick return in case of a workspace query ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO,
    // LAPACK❗✔️:      $                IHI, Z, LDZ, WORK, LWORK, INFO )
    // LAPACK❗✔️: *        ==== Ensure reported workspace size is backward-compatible with
    // LAPACK❗✔️: *        .    previous LAPACK versions. ====
    // LAPACK❗✔️:          WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== copy eigenvalues isolated by DGEBAL ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 10 I = 1, ILO - 1
    // LAPACK❗✔️:             WR( I ) = H( I, I )
    // LAPACK❗✔️:             WI( I ) = ZERO
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️:          DO 20 I = IHI + 1, N
    // LAPACK❗✔️:             WR( I ) = H( I, I )
    // LAPACK❗✔️:             WI( I ) = ZERO
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Initialize Z, if requested ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( INITZ )
    // LAPACK❗✔️:      $      CALL DLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Quick return if possible ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( ILO.EQ.IHI ) THEN
    // LAPACK❗✔️:             WR( ILO ) = H( ILO, ILO )
    // LAPACK❗✔️:             WI( ILO ) = ZERO
    // LAPACK❗✔️:             RETURN
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== DLAHQR/DLAQR0 crossover point ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NMIN = ILAENV( 12, 'DHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N,
    // LAPACK❗✔️:      $          ILO, IHI, LWORK )
    // LAPACK❗✔️:          NMIN = MAX( NTINY, NMIN )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( N.GT.NMIN ) THEN
    // LAPACK❗✔️:             CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO,
    // LAPACK❗✔️:      $                   IHI, Z, LDZ, WORK, LWORK, INFO )
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           ==== Small matrix ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO,
    // LAPACK❗✔️:      $                   IHI, Z, LDZ, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( INFO.GT.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
    // LAPACK❗✔️: *              .    when DLAHQR fails. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                KBOT = INFO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( N.GE.NL ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Larger matrices have enough subdiagonal scratch
    // LAPACK❗✔️: *                 .    space to call DLAQR0 directly. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR,
    // LAPACK❗✔️:      $                         WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 ==== Tiny matrices don't have enough subdiagonal
    // LAPACK❗✔️: *                 .    scratch space to benefit from DLAQR0.  Hence,
    // LAPACK❗✔️: *                 .    tiny matrices must be copied into a larger
    // LAPACK❗✔️: *                 .    array before calling DLAQR0. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DLACPY( 'A', N, N, H, LDH, HL, NL )
    // LAPACK❗✔️:                   HL( N+1, N ) = ZERO
    // LAPACK❗✔️:                   CALL DLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ),
    // LAPACK❗✔️:      $                         NL )
    // LAPACK❗✔️:                   CALL DLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR,
    // LAPACK❗✔️:      $                         WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )
    // LAPACK❗✔️:                   IF( WANTT .OR. INFO.NE.0 )
    // LAPACK❗✔️:      $               CALL DLACPY( 'A', N, N, HL, NL, H, LDH )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Clear out the trash, if necessary. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 )
    // LAPACK❗✔️:      $      CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ==== Ensure reported workspace size is backward-compatible with
    // LAPACK❗✔️: *        .    previous LAPACK versions. ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     ==== End of DHSEQR ====
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if n == 0 {
        return Ok(0);
    }
    for i in 1..ilo {
        wr[i - 1] = h[ix(i, i, ldh)];
        wi[i - 1] = 0.0;
    }
    for i in ihi + 1..=n {
        wr[i - 1] = h[ix(i, i, ldh)];
        wi[i - 1] = 0.0;
    }
    if ilo == ihi {
        wr[ilo - 1] = h[ix(ilo, ilo, ldh)];
        wi[ilo - 1] = 0.0;
        return Ok(0);
    }
    let mut info;
    if n > 75 {
        info = super::multishift::qr(
            true,
            false,
            false,
            n,
            ilo,
            ihi,
            h,
            ldh,
            wr,
            wi,
            ilo,
            ihi,
            &mut [0.0],
            1,
            work,
        )?;
    } else {
        info = lahqr(
            false,
            false,
            n,
            ilo,
            ihi,
            h,
            ldh,
            wr,
            wi,
            ilo,
            ihi,
            &mut [0.0],
            1,
        );
        if info > 0 {
            let kbot = info;
            if n >= 49 {
                info = super::multishift::qr(
                    true,
                    false,
                    false,
                    n,
                    ilo,
                    kbot,
                    h,
                    ldh,
                    wr,
                    wi,
                    ilo,
                    ihi,
                    &mut [0.0],
                    1,
                    work,
                )?;
            } else {
                let mut hl = [0.0; 49 * 49];
                let mut workl = [0.0; 49];
                for j in 1..=n {
                    for i in 1..=n {
                        hl[ix(i, j, 49)] = h[ix(i, j, ldh)];
                    }
                }
                hl[ix(n + 1, n, 49)] = 0.0;
                info = super::multishift::qr(
                    true,
                    false,
                    false,
                    49,
                    ilo,
                    kbot,
                    &mut hl,
                    49,
                    wr,
                    wi,
                    ilo,
                    ihi,
                    &mut [0.0],
                    1,
                    &mut workl,
                )?;
                if info != 0 {
                    for j in 1..=n {
                        for i in 1..=n {
                            h[ix(i, j, ldh)] = hl[ix(i, j, 49)];
                        }
                    }
                }
            }
        }
    }
    if info != 0 && n > 2 {
        for j in 1..=n - 2 {
            for i in j..=n - 2 {
                h[ix(i + 2, j, ldh)] = 0.0;
            }
        }
    }
    Ok(info)
}
