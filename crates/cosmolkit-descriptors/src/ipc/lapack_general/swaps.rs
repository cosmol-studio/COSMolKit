//! Adjacent real-Schur block exchanges used by aggressive early deflation.
use super::super::{lapack_blas as blas, lapack_symmetric::larfg, lapack_tridiagonal::dlamch};
use super::{iamax, ix, lanv2, lartg, swap_within};

/// Fixed LTRANL=false, LTRANR=false, ISGN=-1 from DLAEXC.
fn lasy2(n1: usize, n2: usize, d: &[f64; 16]) -> (f64, [f64; 4]) {
    // LAPACK❗✔️:       SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
    // LAPACK❗✔️:      $                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       LOGICAL            LTRANL, LTRANR
    // LAPACK❗✔️:       INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
    // LAPACK❗✔️:       DOUBLE PRECISION   SCALE, XNORM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
    // LAPACK❗✔️:      $                   X( LDX, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: * =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   TWO, HALF, EIGHT
    // LAPACK❗✔️:       PARAMETER          ( TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            BSWAP, XSWAP
    // LAPACK❗✔️:       INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
    // LAPACK❗✔️:       DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1,
    // LAPACK❗✔️:      $                   TEMP, U11, U12, U22, XMAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
    // LAPACK❗✔️:       INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ),
    // LAPACK❗✔️:      $                   LOCU22( 4 )
    // LAPACK❗✔️:       DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       INTEGER            IDAMAX
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️:       EXTERNAL           IDAMAX, DLAMCH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DCOPY, DSWAP
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, MAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Data statements ..
    // LAPACK❗✔️:       DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / ,
    // LAPACK❗✔️:      $                   LOCU22 / 4, 3, 2, 1 /
    // LAPACK❗✔️:       DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
    // LAPACK❗✔️:       DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Do not check the input parameters for errors
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N1.EQ.0 .OR. N2.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Set constants to control overflow
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       EPS = DLAMCH( 'P' )
    // LAPACK❗✔️:       SMLNUM = DLAMCH( 'S' ) / EPS
    // LAPACK❗✔️:       SGN = ISGN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       K = N1 + N1 + N2 - 2
    // LAPACK❗✔️:       GO TO ( 10, 20, 30, 50 )K
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     1 by 1: TL11*X + SGN*X*TR11 = B11
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
    // LAPACK❗✔️:       BET = ABS( TAU1 )
    // LAPACK❗✔️:       IF( BET.LE.SMLNUM ) THEN
    // LAPACK❗✔️:          TAU1 = SMLNUM
    // LAPACK❗✔️:          BET = SMLNUM
    // LAPACK❗✔️:          INFO = 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SCALE = ONE
    // LAPACK❗✔️:       GAM = ABS( B( 1, 1 ) )
    // LAPACK❗✔️:       IF( SMLNUM*GAM.GT.BET )
    // LAPACK❗✔️:      $   SCALE = ONE / GAM
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
    // LAPACK❗✔️:       XNORM = ABS( X( 1, 1 ) )
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     1 by 2:
    // LAPACK❗✔️: *     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
    // LAPACK❗✔️: *                                       [TR21 TR22]
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    20 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ),
    // LAPACK❗✔️:      $       ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ),
    // LAPACK❗✔️:      $       SMLNUM )
    // LAPACK❗✔️:       TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
    // LAPACK❗✔️:       TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
    // LAPACK❗✔️:       IF( LTRANR ) THEN
    // LAPACK❗✔️:          TMP( 2 ) = SGN*TR( 2, 1 )
    // LAPACK❗✔️:          TMP( 3 ) = SGN*TR( 1, 2 )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          TMP( 2 ) = SGN*TR( 1, 2 )
    // LAPACK❗✔️:          TMP( 3 ) = SGN*TR( 2, 1 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       BTMP( 1 ) = B( 1, 1 )
    // LAPACK❗✔️:       BTMP( 2 ) = B( 1, 2 )
    // LAPACK❗✔️:       GO TO 40
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     2 by 1:
    // LAPACK❗✔️: *          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
    // LAPACK❗✔️: *            [TL21 TL22] [X21]         [X21]         [B21]
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    30 CONTINUE
    // LAPACK❗✔️:       SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ),
    // LAPACK❗✔️:      $       ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ),
    // LAPACK❗✔️:      $       SMLNUM )
    // LAPACK❗✔️:       TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
    // LAPACK❗✔️:       TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
    // LAPACK❗✔️:       IF( LTRANL ) THEN
    // LAPACK❗✔️:          TMP( 2 ) = TL( 1, 2 )
    // LAPACK❗✔️:          TMP( 3 ) = TL( 2, 1 )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          TMP( 2 ) = TL( 2, 1 )
    // LAPACK❗✔️:          TMP( 3 ) = TL( 1, 2 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       BTMP( 1 ) = B( 1, 1 )
    // LAPACK❗✔️:       BTMP( 2 ) = B( 2, 1 )
    // LAPACK❗✔️:    40 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Solve 2 by 2 system using complete pivoting.
    // LAPACK❗✔️: *     Set pivots less than SMIN to SMIN.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IPIV = IDAMAX( 4, TMP, 1 )
    // LAPACK❗✔️:       U11 = TMP( IPIV )
    // LAPACK❗✔️:       IF( ABS( U11 ).LE.SMIN ) THEN
    // LAPACK❗✔️:          INFO = 1
    // LAPACK❗✔️:          U11 = SMIN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       U12 = TMP( LOCU12( IPIV ) )
    // LAPACK❗✔️:       L21 = TMP( LOCL21( IPIV ) ) / U11
    // LAPACK❗✔️:       U22 = TMP( LOCU22( IPIV ) ) - U12*L21
    // LAPACK❗✔️:       XSWAP = XSWPIV( IPIV )
    // LAPACK❗✔️:       BSWAP = BSWPIV( IPIV )
    // LAPACK❗✔️:       IF( ABS( U22 ).LE.SMIN ) THEN
    // LAPACK❗✔️:          INFO = 1
    // LAPACK❗✔️:          U22 = SMIN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( BSWAP ) THEN
    // LAPACK❗✔️:          TEMP = BTMP( 2 )
    // LAPACK❗✔️:          BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
    // LAPACK❗✔️:          BTMP( 1 ) = TEMP
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       SCALE = ONE
    // LAPACK❗✔️:       IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR.
    // LAPACK❗✔️:      $    ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
    // LAPACK❗✔️:          SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
    // LAPACK❗✔️:          BTMP( 1 ) = BTMP( 1 )*SCALE
    // LAPACK❗✔️:          BTMP( 2 ) = BTMP( 2 )*SCALE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       X2( 2 ) = BTMP( 2 ) / U22
    // LAPACK❗✔️:       X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
    // LAPACK❗✔️:       IF( XSWAP ) THEN
    // LAPACK❗✔️:          TEMP = X2( 2 )
    // LAPACK❗✔️:          X2( 2 ) = X2( 1 )
    // LAPACK❗✔️:          X2( 1 ) = TEMP
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       X( 1, 1 ) = X2( 1 )
    // LAPACK❗✔️:       IF( N1.EQ.1 ) THEN
    // LAPACK❗✔️:          X( 1, 2 ) = X2( 2 )
    // LAPACK❗✔️:          XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          X( 2, 1 ) = X2( 2 )
    // LAPACK❗✔️:          XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     2 by 2:
    // LAPACK❗✔️: *     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
    // LAPACK❗✔️: *       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Solve equivalent 4 by 4 system using complete pivoting.
    // LAPACK❗✔️: *     Set pivots less than SMIN to SMIN.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    50 CONTINUE
    // LAPACK❗✔️:       SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
    // LAPACK❗✔️:      $       ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
    // LAPACK❗✔️:       SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
    // LAPACK❗✔️:      $       ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
    // LAPACK❗✔️:       SMIN = MAX( EPS*SMIN, SMLNUM )
    // LAPACK❗✔️:       BTMP( 1 ) = ZERO
    // LAPACK❗✔️:       CALL DCOPY( 16, BTMP, 0, T16, 1 )
    // LAPACK❗✔️:       T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
    // LAPACK❗✔️:       T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
    // LAPACK❗✔️:       T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
    // LAPACK❗✔️:       T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
    // LAPACK❗✔️:       IF( LTRANL ) THEN
    // LAPACK❗✔️:          T16( 1, 2 ) = TL( 2, 1 )
    // LAPACK❗✔️:          T16( 2, 1 ) = TL( 1, 2 )
    // LAPACK❗✔️:          T16( 3, 4 ) = TL( 2, 1 )
    // LAPACK❗✔️:          T16( 4, 3 ) = TL( 1, 2 )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          T16( 1, 2 ) = TL( 1, 2 )
    // LAPACK❗✔️:          T16( 2, 1 ) = TL( 2, 1 )
    // LAPACK❗✔️:          T16( 3, 4 ) = TL( 1, 2 )
    // LAPACK❗✔️:          T16( 4, 3 ) = TL( 2, 1 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( LTRANR ) THEN
    // LAPACK❗✔️:          T16( 1, 3 ) = SGN*TR( 1, 2 )
    // LAPACK❗✔️:          T16( 2, 4 ) = SGN*TR( 1, 2 )
    // LAPACK❗✔️:          T16( 3, 1 ) = SGN*TR( 2, 1 )
    // LAPACK❗✔️:          T16( 4, 2 ) = SGN*TR( 2, 1 )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          T16( 1, 3 ) = SGN*TR( 2, 1 )
    // LAPACK❗✔️:          T16( 2, 4 ) = SGN*TR( 2, 1 )
    // LAPACK❗✔️:          T16( 3, 1 ) = SGN*TR( 1, 2 )
    // LAPACK❗✔️:          T16( 4, 2 ) = SGN*TR( 1, 2 )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       BTMP( 1 ) = B( 1, 1 )
    // LAPACK❗✔️:       BTMP( 2 ) = B( 2, 1 )
    // LAPACK❗✔️:       BTMP( 3 ) = B( 1, 2 )
    // LAPACK❗✔️:       BTMP( 4 ) = B( 2, 2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Perform elimination
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 100 I = 1, 3
    // LAPACK❗✔️:          XMAX = ZERO
    // LAPACK❗✔️:          DO 70 IP = I, 4
    // LAPACK❗✔️:             DO 60 JP = I, 4
    // LAPACK❗✔️:                IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
    // LAPACK❗✔️:                   XMAX = ABS( T16( IP, JP ) )
    // LAPACK❗✔️:                   IPSV = IP
    // LAPACK❗✔️:                   JPSV = JP
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:    60       CONTINUE
    // LAPACK❗✔️:    70    CONTINUE
    // LAPACK❗✔️:          IF( IPSV.NE.I ) THEN
    // LAPACK❗✔️:             CALL DSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
    // LAPACK❗✔️:             TEMP = BTMP( I )
    // LAPACK❗✔️:             BTMP( I ) = BTMP( IPSV )
    // LAPACK❗✔️:             BTMP( IPSV ) = TEMP
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( JPSV.NE.I )
    // LAPACK❗✔️:      $      CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
    // LAPACK❗✔️:          JPIV( I ) = JPSV
    // LAPACK❗✔️:          IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
    // LAPACK❗✔️:             INFO = 1
    // LAPACK❗✔️:             T16( I, I ) = SMIN
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          DO 90 J = I + 1, 4
    // LAPACK❗✔️:             T16( J, I ) = T16( J, I ) / T16( I, I )
    // LAPACK❗✔️:             BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
    // LAPACK❗✔️:             DO 80 K = I + 1, 4
    // LAPACK❗✔️:                T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
    // LAPACK❗✔️:    80       CONTINUE
    // LAPACK❗✔️:    90    CONTINUE
    // LAPACK❗✔️:   100 CONTINUE
    // LAPACK❗✔️:       IF( ABS( T16( 4, 4 ) ).LT.SMIN ) THEN
    // LAPACK❗✔️:          INFO = 1
    // LAPACK❗✔️:          T16( 4, 4 ) = SMIN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       SCALE = ONE
    // LAPACK❗✔️:       IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR.
    // LAPACK❗✔️:      $    ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR.
    // LAPACK❗✔️:      $    ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR.
    // LAPACK❗✔️:      $    ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
    // LAPACK❗✔️:          SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ),
    // LAPACK❗✔️:      $           ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
    // LAPACK❗✔️:          BTMP( 1 ) = BTMP( 1 )*SCALE
    // LAPACK❗✔️:          BTMP( 2 ) = BTMP( 2 )*SCALE
    // LAPACK❗✔️:          BTMP( 3 ) = BTMP( 3 )*SCALE
    // LAPACK❗✔️:          BTMP( 4 ) = BTMP( 4 )*SCALE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       DO 120 I = 1, 4
    // LAPACK❗✔️:          K = 5 - I
    // LAPACK❗✔️:          TEMP = ONE / T16( K, K )
    // LAPACK❗✔️:          TMP( K ) = BTMP( K )*TEMP
    // LAPACK❗✔️:          DO 110 J = K + 1, 4
    // LAPACK❗✔️:             TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
    // LAPACK❗✔️:   110    CONTINUE
    // LAPACK❗✔️:   120 CONTINUE
    // LAPACK❗✔️:       DO 130 I = 1, 3
    // LAPACK❗✔️:          IF( JPIV( 4-I ).NE.4-I ) THEN
    // LAPACK❗✔️:             TEMP = TMP( 4-I )
    // LAPACK❗✔️:             TMP( 4-I ) = TMP( JPIV( 4-I ) )
    // LAPACK❗✔️:             TMP( JPIV( 4-I ) ) = TEMP
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:   130 CONTINUE
    // LAPACK❗✔️:       X( 1, 1 ) = TMP( 1 )
    // LAPACK❗✔️:       X( 2, 1 ) = TMP( 2 )
    // LAPACK❗✔️:       X( 1, 2 ) = TMP( 3 )
    // LAPACK❗✔️:       X( 2, 2 ) = TMP( 4 )
    // LAPACK❗✔️:       XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ),
    // LAPACK❗✔️:      $        ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLASY2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    let eps = dlamch(b'P');
    let smlnum = dlamch(b'S') / eps;
    let mut x = [0.0; 4];
    if n1 == 0 || n2 == 0 {
        return (1.0, x);
    }
    if n1 == 1 && n2 == 1 {
        let mut tau = d[0] - d[5];
        let mut bet = tau.abs();
        if bet <= smlnum {
            tau = smlnum;
            bet = smlnum;
        }
        let gam = d[4].abs();
        let scale = if smlnum * gam > bet { 1.0 / gam } else { 1.0 };
        x[0] = (d[4] * scale) / tau;
        return (scale, x);
    }
    if n1 == 1 || n2 == 1 {
        let (mut tmp, mut rhs, smin) = if n1 == 1 {
            (
                [d[0] - d[5], -d[9], -d[6], d[0] - d[10]],
                [d[4], d[8]],
                (eps * d[0]
                    .abs()
                    .max(d[5].abs())
                    .max(d[9].abs())
                    .max(d[6].abs())
                    .max(d[10].abs()))
                .max(smlnum),
            )
        } else {
            (
                [d[0] - d[10], d[1], d[4], d[5] - d[10]],
                [d[8], d[9]],
                (eps * d[10]
                    .abs()
                    .max(d[0].abs())
                    .max(d[4].abs())
                    .max(d[1].abs())
                    .max(d[5].abs()))
                .max(smlnum),
            )
        };
        let p = iamax(4, &tmp, 1) - 1;
        let mut u11 = tmp[p];
        if u11.abs() <= smin {
            u11 = smin;
        }
        let u12 = tmp[[2, 3, 0, 1][p]];
        let l21 = tmp[[1, 0, 3, 2][p]] / u11;
        let mut u22 = tmp[[3, 2, 1, 0][p]] - u12 * l21;
        if u22.abs() <= smin {
            u22 = smin;
        }
        if [false, true, false, true][p] {
            let value = rhs[1];
            rhs[1] = rhs[0] - l21 * value;
            rhs[0] = value;
        } else {
            rhs[1] -= l21 * rhs[0];
        }
        let mut scale = 1.0;
        if (2.0 * smlnum) * rhs[1].abs() > u22.abs() || (2.0 * smlnum) * rhs[0].abs() > u11.abs() {
            scale = 0.5 / rhs[0].abs().max(rhs[1].abs());
            rhs[0] *= scale;
            rhs[1] *= scale;
        }
        let mut x2 = [0.0; 2];
        x2[1] = rhs[1] / u22;
        x2[0] = rhs[0] / u11 - (u12 / u11) * x2[1];
        if [false, false, true, true][p] {
            x2.swap(0, 1);
        }
        x[0] = x2[0];
        if n1 == 1 {
            x[2] = x2[1];
        } else {
            x[1] = x2[1];
        }
        return (scale, x);
    }
    let mut smin = d[10]
        .abs()
        .max(d[14].abs())
        .max(d[11].abs())
        .max(d[15].abs());
    smin = smin
        .max(d[0].abs())
        .max(d[4].abs())
        .max(d[1].abs())
        .max(d[5].abs());
    smin = (eps * smin).max(smlnum);
    let mut t = [0.0; 16];
    t[0] = d[0] - d[10];
    t[5] = d[5] - d[10];
    t[10] = d[0] - d[15];
    t[15] = d[5] - d[15];
    t[4] = d[4];
    t[1] = d[1];
    t[14] = d[4];
    t[11] = d[1];
    t[8] = -d[11];
    t[13] = -d[11];
    t[2] = -d[14];
    t[7] = -d[14];
    let mut rhs = [d[8], d[9], d[12], d[13]];
    let mut jpiv = [0; 3];
    for i in 0..3 {
        let (mut ip, mut jp, mut xmax) = (i, i, 0.0);
        for row in i..4 {
            for col in i..4 {
                if t[row + col * 4].abs() >= xmax {
                    xmax = t[row + col * 4].abs();
                    ip = row;
                    jp = col;
                }
            }
        }
        if ip != i {
            swap_within(&mut t, 4, ip, 4, i, 4);
            rhs.swap(i, ip);
        }
        if jp != i {
            swap_within(&mut t, 4, jp * 4, 1, i * 4, 1);
        }
        jpiv[i] = jp;
        if t[i + i * 4].abs() < smin {
            t[i + i * 4] = smin;
        }
        for j in i + 1..4 {
            t[j + i * 4] /= t[i + i * 4];
            rhs[j] -= t[j + i * 4] * rhs[i];
            for k in i + 1..4 {
                t[j + k * 4] -= t[j + i * 4] * t[i + k * 4];
            }
        }
    }
    if t[15].abs() < smin {
        t[15] = smin;
    }
    let mut scale = 1.0;
    if (0..4).any(|i| (8.0 * smlnum) * rhs[i].abs() > t[i + i * 4].abs()) {
        scale = (1.0 / 8.0)
            / rhs[0]
                .abs()
                .max(rhs[1].abs())
                .max(rhs[2].abs())
                .max(rhs[3].abs());
        for value in &mut rhs {
            *value *= scale;
        }
    }
    let mut tmp = [0.0; 4];
    for k in (0..4).rev() {
        let temp = 1.0 / t[k + k * 4];
        tmp[k] = rhs[k] * temp;
        for j in k + 1..4 {
            tmp[k] -= (temp * t[k + j * 4]) * tmp[j];
        }
    }
    for i in (0..3).rev() {
        if jpiv[i] != i {
            tmp.swap(i, jpiv[i]);
        }
    }
    (scale, tmp)
}

/// DLAEXC only calls the explicitly unrolled order-three DLARFX path.
fn larfx3(left: bool, m: usize, n: usize, v: &[f64; 3], tau: f64, c: &mut [f64], ldc: usize) {
    // LAPACK❗✔️:       SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          SIDE
    // LAPACK❗✔️:       INTEGER            LDC, M, N
    // LAPACK❗✔️:       DOUBLE PRECISION   TAU
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            J
    // LAPACK❗✔️:       DOUBLE PRECISION   SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9,
    // LAPACK❗✔️:      $                   V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLARF
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( TAU.EQ.ZERO )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️:       IF( LSAME( SIDE, 'L' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Form  H * C, where H has order m.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          GO TO ( 10, 30, 50, 70, 90, 110, 130, 150,
    // LAPACK❗✔️:      $           170, 190 )M
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Code for general M
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 1 x 1 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T1 = ONE - TAU*V( 1 )*V( 1 )
    // LAPACK❗✔️:          DO 20 J = 1, N
    // LAPACK❗✔️:             C( 1, J ) = T1*C( 1, J )
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:    30    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 2 x 2 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          DO 40 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:    50    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 3 x 3 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          DO 60 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:    60    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:    70    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 4 x 4 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          DO 80 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:    80    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:    90    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 5 x 5 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          DO 100 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J ) + V5*C( 5, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:             C( 5, J ) = C( 5, J ) - SUM*T5
    // LAPACK❗✔️:   100    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   110    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 6 x 6 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          DO 120 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:             C( 5, J ) = C( 5, J ) - SUM*T5
    // LAPACK❗✔️:             C( 6, J ) = C( 6, J ) - SUM*T6
    // LAPACK❗✔️:   120    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   130    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 7 x 7 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          DO 140 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
    // LAPACK❗✔️:      $            V7*C( 7, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:             C( 5, J ) = C( 5, J ) - SUM*T5
    // LAPACK❗✔️:             C( 6, J ) = C( 6, J ) - SUM*T6
    // LAPACK❗✔️:             C( 7, J ) = C( 7, J ) - SUM*T7
    // LAPACK❗✔️:   140    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   150    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 8 x 8 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          V8 = V( 8 )
    // LAPACK❗✔️:          T8 = TAU*V8
    // LAPACK❗✔️:          DO 160 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
    // LAPACK❗✔️:      $            V7*C( 7, J ) + V8*C( 8, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:             C( 5, J ) = C( 5, J ) - SUM*T5
    // LAPACK❗✔️:             C( 6, J ) = C( 6, J ) - SUM*T6
    // LAPACK❗✔️:             C( 7, J ) = C( 7, J ) - SUM*T7
    // LAPACK❗✔️:             C( 8, J ) = C( 8, J ) - SUM*T8
    // LAPACK❗✔️:   160    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   170    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 9 x 9 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          V8 = V( 8 )
    // LAPACK❗✔️:          T8 = TAU*V8
    // LAPACK❗✔️:          V9 = V( 9 )
    // LAPACK❗✔️:          T9 = TAU*V9
    // LAPACK❗✔️:          DO 180 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
    // LAPACK❗✔️:      $            V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:             C( 5, J ) = C( 5, J ) - SUM*T5
    // LAPACK❗✔️:             C( 6, J ) = C( 6, J ) - SUM*T6
    // LAPACK❗✔️:             C( 7, J ) = C( 7, J ) - SUM*T7
    // LAPACK❗✔️:             C( 8, J ) = C( 8, J ) - SUM*T8
    // LAPACK❗✔️:             C( 9, J ) = C( 9, J ) - SUM*T9
    // LAPACK❗✔️:   180    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   190    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 10 x 10 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          V8 = V( 8 )
    // LAPACK❗✔️:          T8 = TAU*V8
    // LAPACK❗✔️:          V9 = V( 9 )
    // LAPACK❗✔️:          T9 = TAU*V9
    // LAPACK❗✔️:          V10 = V( 10 )
    // LAPACK❗✔️:          T10 = TAU*V10
    // LAPACK❗✔️:          DO 200 J = 1, N
    // LAPACK❗✔️:             SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
    // LAPACK❗✔️:      $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
    // LAPACK❗✔️:      $            V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) +
    // LAPACK❗✔️:      $            V10*C( 10, J )
    // LAPACK❗✔️:             C( 1, J ) = C( 1, J ) - SUM*T1
    // LAPACK❗✔️:             C( 2, J ) = C( 2, J ) - SUM*T2
    // LAPACK❗✔️:             C( 3, J ) = C( 3, J ) - SUM*T3
    // LAPACK❗✔️:             C( 4, J ) = C( 4, J ) - SUM*T4
    // LAPACK❗✔️:             C( 5, J ) = C( 5, J ) - SUM*T5
    // LAPACK❗✔️:             C( 6, J ) = C( 6, J ) - SUM*T6
    // LAPACK❗✔️:             C( 7, J ) = C( 7, J ) - SUM*T7
    // LAPACK❗✔️:             C( 8, J ) = C( 8, J ) - SUM*T8
    // LAPACK❗✔️:             C( 9, J ) = C( 9, J ) - SUM*T9
    // LAPACK❗✔️:             C( 10, J ) = C( 10, J ) - SUM*T10
    // LAPACK❗✔️:   200    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Form  C * H, where H has order n.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          GO TO ( 210, 230, 250, 270, 290, 310, 330, 350,
    // LAPACK❗✔️:      $           370, 390 )N
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Code for general N
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   210    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 1 x 1 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T1 = ONE - TAU*V( 1 )*V( 1 )
    // LAPACK❗✔️:          DO 220 J = 1, M
    // LAPACK❗✔️:             C( J, 1 ) = T1*C( J, 1 )
    // LAPACK❗✔️:   220    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   230    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 2 x 2 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          DO 240 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:   240    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   250    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 3 x 3 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          DO 260 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:   260    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   270    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 4 x 4 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          DO 280 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:   280    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   290    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 5 x 5 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          DO 300 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 ) + V5*C( J, 5 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:             C( J, 5 ) = C( J, 5 ) - SUM*T5
    // LAPACK❗✔️:   300    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   310    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 6 x 6 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          DO 320 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:             C( J, 5 ) = C( J, 5 ) - SUM*T5
    // LAPACK❗✔️:             C( J, 6 ) = C( J, 6 ) - SUM*T6
    // LAPACK❗✔️:   320    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   330    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 7 x 7 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          DO 340 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
    // LAPACK❗✔️:      $            V7*C( J, 7 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:             C( J, 5 ) = C( J, 5 ) - SUM*T5
    // LAPACK❗✔️:             C( J, 6 ) = C( J, 6 ) - SUM*T6
    // LAPACK❗✔️:             C( J, 7 ) = C( J, 7 ) - SUM*T7
    // LAPACK❗✔️:   340    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   350    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 8 x 8 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          V8 = V( 8 )
    // LAPACK❗✔️:          T8 = TAU*V8
    // LAPACK❗✔️:          DO 360 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
    // LAPACK❗✔️:      $            V7*C( J, 7 ) + V8*C( J, 8 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:             C( J, 5 ) = C( J, 5 ) - SUM*T5
    // LAPACK❗✔️:             C( J, 6 ) = C( J, 6 ) - SUM*T6
    // LAPACK❗✔️:             C( J, 7 ) = C( J, 7 ) - SUM*T7
    // LAPACK❗✔️:             C( J, 8 ) = C( J, 8 ) - SUM*T8
    // LAPACK❗✔️:   360    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   370    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 9 x 9 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          V8 = V( 8 )
    // LAPACK❗✔️:          T8 = TAU*V8
    // LAPACK❗✔️:          V9 = V( 9 )
    // LAPACK❗✔️:          T9 = TAU*V9
    // LAPACK❗✔️:          DO 380 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
    // LAPACK❗✔️:      $            V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:             C( J, 5 ) = C( J, 5 ) - SUM*T5
    // LAPACK❗✔️:             C( J, 6 ) = C( J, 6 ) - SUM*T6
    // LAPACK❗✔️:             C( J, 7 ) = C( J, 7 ) - SUM*T7
    // LAPACK❗✔️:             C( J, 8 ) = C( J, 8 ) - SUM*T8
    // LAPACK❗✔️:             C( J, 9 ) = C( J, 9 ) - SUM*T9
    // LAPACK❗✔️:   380    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:   390    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Special code for 10 x 10 Householder
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          V1 = V( 1 )
    // LAPACK❗✔️:          T1 = TAU*V1
    // LAPACK❗✔️:          V2 = V( 2 )
    // LAPACK❗✔️:          T2 = TAU*V2
    // LAPACK❗✔️:          V3 = V( 3 )
    // LAPACK❗✔️:          T3 = TAU*V3
    // LAPACK❗✔️:          V4 = V( 4 )
    // LAPACK❗✔️:          T4 = TAU*V4
    // LAPACK❗✔️:          V5 = V( 5 )
    // LAPACK❗✔️:          T5 = TAU*V5
    // LAPACK❗✔️:          V6 = V( 6 )
    // LAPACK❗✔️:          T6 = TAU*V6
    // LAPACK❗✔️:          V7 = V( 7 )
    // LAPACK❗✔️:          T7 = TAU*V7
    // LAPACK❗✔️:          V8 = V( 8 )
    // LAPACK❗✔️:          T8 = TAU*V8
    // LAPACK❗✔️:          V9 = V( 9 )
    // LAPACK❗✔️:          T9 = TAU*V9
    // LAPACK❗✔️:          V10 = V( 10 )
    // LAPACK❗✔️:          T10 = TAU*V10
    // LAPACK❗✔️:          DO 400 J = 1, M
    // LAPACK❗✔️:             SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
    // LAPACK❗✔️:      $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
    // LAPACK❗✔️:      $            V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) +
    // LAPACK❗✔️:      $            V10*C( J, 10 )
    // LAPACK❗✔️:             C( J, 1 ) = C( J, 1 ) - SUM*T1
    // LAPACK❗✔️:             C( J, 2 ) = C( J, 2 ) - SUM*T2
    // LAPACK❗✔️:             C( J, 3 ) = C( J, 3 ) - SUM*T3
    // LAPACK❗✔️:             C( J, 4 ) = C( J, 4 ) - SUM*T4
    // LAPACK❗✔️:             C( J, 5 ) = C( J, 5 ) - SUM*T5
    // LAPACK❗✔️:             C( J, 6 ) = C( J, 6 ) - SUM*T6
    // LAPACK❗✔️:             C( J, 7 ) = C( J, 7 ) - SUM*T7
    // LAPACK❗✔️:             C( J, 8 ) = C( J, 8 ) - SUM*T8
    // LAPACK❗✔️:             C( J, 9 ) = C( J, 9 ) - SUM*T9
    // LAPACK❗✔️:             C( J, 10 ) = C( J, 10 ) - SUM*T10
    // LAPACK❗✔️:   400    CONTINUE
    // LAPACK❗✔️:          GO TO 410
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:   410 CONTINUE
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLARFX
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if tau == 0.0 {
        return;
    }
    let (v1, v2, v3) = (v[0], v[1], v[2]);
    let (t1, t2, t3) = (tau * v1, tau * v2, tau * v3);
    if left {
        for j in 1..=n {
            let sum = v1 * c[ix(1, j, ldc)] + v2 * c[ix(2, j, ldc)] + v3 * c[ix(3, j, ldc)];
            c[ix(1, j, ldc)] -= sum * t1;
            c[ix(2, j, ldc)] -= sum * t2;
            c[ix(3, j, ldc)] -= sum * t3;
        }
    } else {
        for j in 1..=m {
            let sum = v1 * c[ix(j, 1, ldc)] + v2 * c[ix(j, 2, ldc)] + v3 * c[ix(j, 3, ldc)];
            c[ix(j, 1, ldc)] -= sum * t1;
            c[ix(j, 2, ldc)] -= sum * t2;
            c[ix(j, 3, ldc)] -= sum * t3;
        }
    }
}

fn standardize(
    wantq: bool,
    n: usize,
    t: &mut [f64],
    ldt: usize,
    q: &mut [f64],
    ldq: usize,
    j: usize,
) {
    let p = lanv2(
        t[ix(j, j, ldt)],
        t[ix(j, j + 1, ldt)],
        t[ix(j + 1, j, ldt)],
        t[ix(j + 1, j + 1, ldt)],
    );
    t[ix(j, j, ldt)] = p.a;
    t[ix(j, j + 1, ldt)] = p.b;
    t[ix(j + 1, j, ldt)] = p.c;
    t[ix(j + 1, j + 1, ldt)] = p.d;
    if j + 2 <= n {
        blas::rot_within(
            n - j - 1,
            t,
            ix(j, j + 2, ldt),
            ldt,
            ix(j + 1, j + 2, ldt),
            ldt,
            p.cs,
            p.sn,
        );
    }
    blas::rot_within(j - 1, t, ix(1, j, ldt), 1, ix(1, j + 1, ldt), 1, p.cs, p.sn);
    if wantq {
        blas::rot_within(n, q, ix(1, j, ldq), 1, ix(1, j + 1, ldq), 1, p.cs, p.sn);
    }
}

#[allow(clippy::too_many_arguments)]
fn laexc(
    wantq: bool,
    n: usize,
    t: &mut [f64],
    ldt: usize,
    q: &mut [f64],
    ldq: usize,
    j1: usize,
    n1: usize,
    n2: usize,
) -> usize {
    // LAPACK❗✔️:       SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
    // LAPACK❗✔️:      $                   INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       LOGICAL            WANTQ
    // LAPACK❗✔️:       INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   TEN
    // LAPACK❗✔️:       PARAMETER          ( TEN = 1.0D+1 )
    // LAPACK❗✔️:       INTEGER            LDD, LDX
    // LAPACK❗✔️:       PARAMETER          ( LDD = 4, LDX = 2 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            IERR, J2, J3, J4, K, ND
    // LAPACK❗✔️:       DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
    // LAPACK❗✔️:      $                   T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2,
    // LAPACK❗✔️:      $                   WR1, WR2, XNORM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ),
    // LAPACK❗✔️:      $                   X( LDX, 2 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DLANGE
    // LAPACK❗✔️:       EXTERNAL           DLAMCH, DLANGE
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLACPY, DLANV2, DLARFG, DLARFX, DLARTG, DLASY2,
    // LAPACK❗✔️:      $                   DROT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, MAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️:       IF( J1+N1.GT.N )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       J2 = J1 + 1
    // LAPACK❗✔️:       J3 = J1 + 2
    // LAPACK❗✔️:       J4 = J1 + 3
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N1.EQ.1 .AND. N2.EQ.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Swap two 1-by-1 blocks.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T11 = T( J1, J1 )
    // LAPACK❗✔️:          T22 = T( J2, J2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Determine the transformation to perform the interchange.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Apply transformation to the matrix T.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( J3.LE.N )
    // LAPACK❗✔️:      $      CALL DROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS,
    // LAPACK❗✔️:      $                 SN )
    // LAPACK❗✔️:          CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T( J1, J1 ) = T22
    // LAPACK❗✔️:          T( J2, J2 ) = T11
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTQ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Accumulate transformation in the matrix Q.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Swapping involves at least one 2-by-2 block.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Copy the diagonal block of order N1+N2 to the local array D
    // LAPACK❗✔️: *        and compute its norm.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ND = N1 + N2
    // LAPACK❗✔️:          CALL DLACPY( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
    // LAPACK❗✔️:          DNORM = DLANGE( 'Max', ND, ND, D, LDD, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute machine-dependent threshold for test for accepting
    // LAPACK❗✔️: *        swap.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          EPS = DLAMCH( 'P' )
    // LAPACK❗✔️:          SMLNUM = DLAMCH( 'S' ) / EPS
    // LAPACK❗✔️:          THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Solve T11*X - X*T22 = scale*T12 for X.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD,
    // LAPACK❗✔️:      $                D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X,
    // LAPACK❗✔️:      $                LDX, XNORM, IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Swap the adjacent diagonal blocks.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          K = N1 + N1 + N2 - 3
    // LAPACK❗✔️:          GO TO ( 10, 20, 30 )K
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        N1 = 1, N2 = 2: generate elementary reflector H so that:
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        ( scale, X11, X12 ) H = ( 0, 0, * )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          U( 1 ) = SCALE
    // LAPACK❗✔️:          U( 2 ) = X( 1, 1 )
    // LAPACK❗✔️:          U( 3 ) = X( 1, 2 )
    // LAPACK❗✔️:          CALL DLARFG( 3, U( 3 ), U, 1, TAU )
    // LAPACK❗✔️:          U( 3 ) = ONE
    // LAPACK❗✔️:          T11 = T( J1, J1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Perform swap provisionally on diagonal block in D.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Test whether to reject swap.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3,
    // LAPACK❗✔️:      $       3 )-T11 ) ).GT.THRESH )GO TO 50
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Accept swap: apply transformation to the entire matrix T.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T( J3, J1 ) = ZERO
    // LAPACK❗✔️:          T( J3, J2 ) = ZERO
    // LAPACK❗✔️:          T( J3, J3 ) = T11
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTQ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Accumulate transformation in the matrix Q.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          GO TO 40
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        N1 = 2, N2 = 1: generate elementary reflector H so that:
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        H (  -X11 ) = ( * )
    // LAPACK❗✔️: *          (  -X21 ) = ( 0 )
    // LAPACK❗✔️: *          ( scale ) = ( 0 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          U( 1 ) = -X( 1, 1 )
    // LAPACK❗✔️:          U( 2 ) = -X( 2, 1 )
    // LAPACK❗✔️:          U( 3 ) = SCALE
    // LAPACK❗✔️:          CALL DLARFG( 3, U( 1 ), U( 2 ), 1, TAU )
    // LAPACK❗✔️:          U( 1 ) = ONE
    // LAPACK❗✔️:          T33 = T( J3, J3 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Perform swap provisionally on diagonal block in D.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Test whether to reject swap.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1,
    // LAPACK❗✔️:      $       1 )-T33 ) ).GT.THRESH )GO TO 50
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Accept swap: apply transformation to the entire matrix T.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFX( 'R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T( J1, J1 ) = T33
    // LAPACK❗✔️:          T( J2, J1 ) = ZERO
    // LAPACK❗✔️:          T( J3, J1 ) = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTQ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Accumulate transformation in the matrix Q.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          GO TO 40
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    30    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
    // LAPACK❗✔️: *        that:
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        H(2) H(1) (  -X11  -X12 ) = (  *  * )
    // LAPACK❗✔️: *                  (  -X21  -X22 )   (  0  * )
    // LAPACK❗✔️: *                  ( scale    0  )   (  0  0 )
    // LAPACK❗✔️: *                  (    0  scale )   (  0  0 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          U1( 1 ) = -X( 1, 1 )
    // LAPACK❗✔️:          U1( 2 ) = -X( 2, 1 )
    // LAPACK❗✔️:          U1( 3 ) = SCALE
    // LAPACK❗✔️:          CALL DLARFG( 3, U1( 1 ), U1( 2 ), 1, TAU1 )
    // LAPACK❗✔️:          U1( 1 ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          TEMP = -TAU1*( X( 1, 2 )+U1( 2 )*X( 2, 2 ) )
    // LAPACK❗✔️:          U2( 1 ) = -TEMP*U1( 2 ) - X( 2, 2 )
    // LAPACK❗✔️:          U2( 2 ) = -TEMP*U1( 3 )
    // LAPACK❗✔️:          U2( 3 ) = SCALE
    // LAPACK❗✔️:          CALL DLARFG( 3, U2( 1 ), U2( 2 ), 1, TAU2 )
    // LAPACK❗✔️:          U2( 1 ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Perform swap provisionally on diagonal block in D.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, 4, U1, TAU1, D, LDD, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', 4, 3, U1, TAU1, D, LDD, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, 4, U2, TAU2, D( 2, 1 ), LDD, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', 4, 3, U2, TAU2, D( 1, 2 ), LDD, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Test whether to reject swap.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ),
    // LAPACK❗✔️:      $       ABS( D( 4, 2 ) ) ).GT.THRESH )GO TO 50
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Accept swap: apply transformation to the entire matrix T.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, N-J1+1, U1, TAU1, T( J1, J1 ), LDT, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', J4, 3, U1, TAU1, T( 1, J1 ), LDT, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'L', 3, N-J1+1, U2, TAU2, T( J2, J1 ), LDT, WORK )
    // LAPACK❗✔️:          CALL DLARFX( 'R', J4, 3, U2, TAU2, T( 1, J2 ), LDT, WORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          T( J3, J1 ) = ZERO
    // LAPACK❗✔️:          T( J3, J2 ) = ZERO
    // LAPACK❗✔️:          T( J4, J1 ) = ZERO
    // LAPACK❗✔️:          T( J4, J2 ) = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTQ ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Accumulate transformation in the matrix Q.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFX( 'R', N, 3, U1, TAU1, Q( 1, J1 ), LDQ, WORK )
    // LAPACK❗✔️:             CALL DLARFX( 'R', N, 3, U2, TAU2, Q( 1, J2 ), LDQ, WORK )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( N2.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Standardize new 2-by-2 block T11
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ),
    // LAPACK❗✔️:      $                   T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
    // LAPACK❗✔️:             CALL DROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT,
    // LAPACK❗✔️:      $                 CS, SN )
    // LAPACK❗✔️:             CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
    // LAPACK❗✔️:             IF( WANTQ )
    // LAPACK❗✔️:      $         CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( N1.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Standardize new 2-by-2 block T22
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             J3 = J1 + N2
    // LAPACK❗✔️:             J4 = J3 + 1
    // LAPACK❗✔️:             CALL DLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ),
    // LAPACK❗✔️:      $                   T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
    // LAPACK❗✔️:             IF( J3+2.LE.N )
    // LAPACK❗✔️:      $         CALL DROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ),
    // LAPACK❗✔️:      $                    LDT, CS, SN )
    // LAPACK❗✔️:             CALL DROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
    // LAPACK❗✔️:             IF( WANTQ )
    // LAPACK❗✔️:      $         CALL DROT( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Exit with INFO = 1 if swap was rejected.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    50 CONTINUE
    // LAPACK❗✔️:       INFO = 1
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLAEXC
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if n == 0 || n1 == 0 || n2 == 0 || j1 + n1 > n {
        return 0;
    }
    let j2 = j1 + 1;
    let j3 = j1 + 2;
    let j4 = j1 + 3;
    if n1 == 1 && n2 == 1 {
        let t11 = t[ix(j1, j1, ldt)];
        let t22 = t[ix(j2, j2, ldt)];
        let (cs, sn, _) = lartg(t[ix(j1, j2, ldt)], t22 - t11);
        if j3 <= n {
            blas::rot_within(
                n - j1 - 1,
                t,
                ix(j1, j3, ldt),
                ldt,
                ix(j2, j3, ldt),
                ldt,
                cs,
                sn,
            );
        }
        blas::rot_within(j1 - 1, t, ix(1, j1, ldt), 1, ix(1, j2, ldt), 1, cs, sn);
        t[ix(j1, j1, ldt)] = t22;
        t[ix(j2, j2, ldt)] = t11;
        if wantq {
            blas::rot_within(n, q, ix(1, j1, ldq), 1, ix(1, j2, ldq), 1, cs, sn);
        }
        return 0;
    }
    let nd = n1 + n2;
    let mut d = [0.0; 16];
    for j in 0..nd {
        for i in 0..nd {
            d[i + 4 * j] = t[ix(j1 + i, j1 + j, ldt)];
        }
    }
    let mut dnorm: f64 = 0.0;
    for j in 0..nd {
        for i in 0..nd {
            dnorm = dnorm.max(d[i + 4 * j].abs());
        }
    }
    let eps = dlamch(b'P');
    let smlnum = dlamch(b'S') / eps;
    let thresh = (10.0 * eps * dnorm).max(smlnum);
    let (scale, x) = lasy2(n1, n2, &d);
    if n1 == 1 {
        let mut values = [x[2], scale, x[0]];
        let tau = larfg(&mut values);
        let u = [values[1], values[2], 1.0];
        let t11 = t[ix(j1, j1, ldt)];
        larfx3(true, 3, 3, &u, tau, &mut d, 4);
        larfx3(false, 3, 3, &u, tau, &mut d, 4);
        if d[2].abs().max(d[6].abs()).max((d[10] - t11).abs()) > thresh {
            return 1;
        }
        larfx3(true, 3, n - j1 + 1, &u, tau, &mut t[ix(j1, j1, ldt)..], ldt);
        larfx3(false, j2, 3, &u, tau, &mut t[ix(1, j1, ldt)..], ldt);
        t[ix(j3, j1, ldt)] = 0.0;
        t[ix(j3, j2, ldt)] = 0.0;
        t[ix(j3, j3, ldt)] = t11;
        if wantq {
            larfx3(false, n, 3, &u, tau, &mut q[ix(1, j1, ldq)..], ldq);
        }
    } else if n2 == 1 {
        let mut u = [-x[0], -x[1], scale];
        let tau = larfg(&mut u);
        u[0] = 1.0;
        let t33 = t[ix(j3, j3, ldt)];
        larfx3(true, 3, 3, &u, tau, &mut d, 4);
        larfx3(false, 3, 3, &u, tau, &mut d, 4);
        if d[1].abs().max(d[2].abs()).max((d[0] - t33).abs()) > thresh {
            return 1;
        }
        larfx3(false, j3, 3, &u, tau, &mut t[ix(1, j1, ldt)..], ldt);
        larfx3(true, 3, n - j1, &u, tau, &mut t[ix(j1, j2, ldt)..], ldt);
        t[ix(j1, j1, ldt)] = t33;
        t[ix(j2, j1, ldt)] = 0.0;
        t[ix(j3, j1, ldt)] = 0.0;
        if wantq {
            larfx3(false, n, 3, &u, tau, &mut q[ix(1, j1, ldq)..], ldq);
        }
    } else {
        let mut u1 = [-x[0], -x[1], scale];
        let tau1 = larfg(&mut u1);
        u1[0] = 1.0;
        let temp = -tau1 * (x[2] + u1[1] * x[3]);
        let mut u2 = [-temp * u1[1] - x[3], -temp * u1[2], scale];
        let tau2 = larfg(&mut u2);
        u2[0] = 1.0;
        larfx3(true, 3, 4, &u1, tau1, &mut d, 4);
        larfx3(false, 4, 3, &u1, tau1, &mut d, 4);
        larfx3(true, 3, 4, &u2, tau2, &mut d[1..], 4);
        larfx3(false, 4, 3, &u2, tau2, &mut d[4..], 4);
        if d[2].abs().max(d[6].abs()).max(d[3].abs()).max(d[7].abs()) > thresh {
            return 1;
        }
        larfx3(
            true,
            3,
            n - j1 + 1,
            &u1,
            tau1,
            &mut t[ix(j1, j1, ldt)..],
            ldt,
        );
        larfx3(false, j4, 3, &u1, tau1, &mut t[ix(1, j1, ldt)..], ldt);
        larfx3(
            true,
            3,
            n - j1 + 1,
            &u2,
            tau2,
            &mut t[ix(j2, j1, ldt)..],
            ldt,
        );
        larfx3(false, j4, 3, &u2, tau2, &mut t[ix(1, j2, ldt)..], ldt);
        t[ix(j3, j1, ldt)] = 0.0;
        t[ix(j3, j2, ldt)] = 0.0;
        t[ix(j4, j1, ldt)] = 0.0;
        t[ix(j4, j2, ldt)] = 0.0;
        if wantq {
            larfx3(false, n, 3, &u1, tau1, &mut q[ix(1, j1, ldq)..], ldq);
            larfx3(false, n, 3, &u2, tau2, &mut q[ix(1, j2, ldq)..], ldq);
        }
    }
    if n2 == 2 {
        standardize(wantq, n, t, ldt, q, ldq, j1);
    }
    if n1 == 2 {
        standardize(wantq, n, t, ldt, q, ldq, j1 + n2);
    }
    0
}

#[allow(clippy::too_many_arguments)]
pub(super) fn trexc(
    n: usize,
    t: &mut [f64],
    ldt: usize,
    q: &mut [f64],
    ldq: usize,
    ifst: &mut usize,
    ilst: &mut usize,
) -> usize {
    // LAPACK❗✔️:       SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
    // LAPACK❗✔️:      $                   INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          COMPQ
    // LAPACK❗✔️:       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            WANTQ
    // LAPACK❗✔️:       INTEGER            HERE, NBF, NBL, NBNEXT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLAEXC, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Decode and test the input arguments.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       WANTQ = LSAME( COMPQ, 'V' )
    // LAPACK❗✔️:       IF( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
    // LAPACK❗✔️:          INFO = -6
    // LAPACK❗✔️:       ELSE IF(( IFST.LT.1 .OR. IFST.GT.N ).AND.( N.GT.0 )) THEN
    // LAPACK❗✔️:          INFO = -7
    // LAPACK❗✔️:       ELSE IF(( ILST.LT.1 .OR. ILST.GT.N ).AND.( N.GT.0 )) THEN
    // LAPACK❗✔️:          INFO = -8
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DTREXC', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.1 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Determine the first row of specified block
    // LAPACK❗✔️: *     and find out it is 1 by 1 or 2 by 2.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( IFST.GT.1 ) THEN
    // LAPACK❗✔️:          IF( T( IFST, IFST-1 ).NE.ZERO )
    // LAPACK❗✔️:      $      IFST = IFST - 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       NBF = 1
    // LAPACK❗✔️:       IF( IFST.LT.N ) THEN
    // LAPACK❗✔️:          IF( T( IFST+1, IFST ).NE.ZERO )
    // LAPACK❗✔️:      $      NBF = 2
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Determine the first row of the final block
    // LAPACK❗✔️: *     and find out it is 1 by 1 or 2 by 2.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ILST.GT.1 ) THEN
    // LAPACK❗✔️:          IF( T( ILST, ILST-1 ).NE.ZERO )
    // LAPACK❗✔️:      $      ILST = ILST - 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       NBL = 1
    // LAPACK❗✔️:       IF( ILST.LT.N ) THEN
    // LAPACK❗✔️:          IF( T( ILST+1, ILST ).NE.ZERO )
    // LAPACK❗✔️:      $      NBL = 2
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( IFST.EQ.ILST )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( IFST.LT.ILST ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Update ILST
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( NBF.EQ.2 .AND. NBL.EQ.1 )
    // LAPACK❗✔️:      $      ILST = ILST - 1
    // LAPACK❗✔️:          IF( NBF.EQ.1 .AND. NBL.EQ.2 )
    // LAPACK❗✔️:      $      ILST = ILST + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          HERE = IFST
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Swap block with next one below
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Current block either 1 by 1 or 2 by 2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NBNEXT = 1
    // LAPACK❗✔️:             IF( HERE+NBF+1.LE.N ) THEN
    // LAPACK❗✔️:                IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO )
    // LAPACK❗✔️:      $            NBNEXT = 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT,
    // LAPACK❗✔️:      $                   WORK, INFO )
    // LAPACK❗✔️:             IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:                ILST = HERE
    // LAPACK❗✔️:                RETURN
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             HERE = HERE + NBNEXT
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Test if 2 by 2 block breaks into two 1 by 1 blocks
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( NBF.EQ.2 ) THEN
    // LAPACK❗✔️:                IF( T( HERE+1, HERE ).EQ.ZERO )
    // LAPACK❗✔️:      $            NBF = 3
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Current block consists of two 1 by 1 blocks each of which
    // LAPACK❗✔️: *           must be swapped individually
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NBNEXT = 1
    // LAPACK❗✔️:             IF( HERE+3.LE.N ) THEN
    // LAPACK❗✔️:                IF( T( HERE+3, HERE+2 ).NE.ZERO )
    // LAPACK❗✔️:      $            NBNEXT = 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT,
    // LAPACK❗✔️:      $                   WORK, INFO )
    // LAPACK❗✔️:             IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:                ILST = HERE
    // LAPACK❗✔️:                RETURN
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( NBNEXT.EQ.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Swap two 1 by 1 blocks, no problems possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT,
    // LAPACK❗✔️:      $                      WORK, INFO )
    // LAPACK❗✔️:                HERE = HERE + 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Recompute NBNEXT in case 2 by 2 split
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( T( HERE+2, HERE+1 ).EQ.ZERO )
    // LAPACK❗✔️:      $            NBNEXT = 1
    // LAPACK❗✔️:                IF( NBNEXT.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 2 by 2 Block did not split
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1,
    // LAPACK❗✔️:      $                         NBNEXT, WORK, INFO )
    // LAPACK❗✔️:                   IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:                      ILST = HERE
    // LAPACK❗✔️:                      RETURN
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                   HERE = HERE + 2
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 2 by 2 Block did split
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1,
    // LAPACK❗✔️:      $                         WORK, INFO )
    // LAPACK❗✔️:                   CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1,
    // LAPACK❗✔️:      $                         WORK, INFO )
    // LAPACK❗✔️:                   HERE = HERE + 2
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( HERE.LT.ILST )
    // LAPACK❗✔️:      $      GO TO 10
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          HERE = IFST
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Swap block with next one above
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Current block either 1 by 1 or 2 by 2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NBNEXT = 1
    // LAPACK❗✔️:             IF( HERE.GE.3 ) THEN
    // LAPACK❗✔️:                IF( T( HERE-1, HERE-2 ).NE.ZERO )
    // LAPACK❗✔️:      $            NBNEXT = 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT,
    // LAPACK❗✔️:      $                   NBF, WORK, INFO )
    // LAPACK❗✔️:             IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:                ILST = HERE
    // LAPACK❗✔️:                RETURN
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             HERE = HERE - NBNEXT
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Test if 2 by 2 block breaks into two 1 by 1 blocks
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( NBF.EQ.2 ) THEN
    // LAPACK❗✔️:                IF( T( HERE+1, HERE ).EQ.ZERO )
    // LAPACK❗✔️:      $            NBF = 3
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Current block consists of two 1 by 1 blocks each of which
    // LAPACK❗✔️: *           must be swapped individually
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NBNEXT = 1
    // LAPACK❗✔️:             IF( HERE.GE.3 ) THEN
    // LAPACK❗✔️:                IF( T( HERE-1, HERE-2 ).NE.ZERO )
    // LAPACK❗✔️:      $            NBNEXT = 2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT,
    // LAPACK❗✔️:      $                   1, WORK, INFO )
    // LAPACK❗✔️:             IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:                ILST = HERE
    // LAPACK❗✔️:                RETURN
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( NBNEXT.EQ.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Swap two 1 by 1 blocks, no problems possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1,
    // LAPACK❗✔️:      $                      WORK, INFO )
    // LAPACK❗✔️:                HERE = HERE - 1
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Recompute NBNEXT in case 2 by 2 split
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( T( HERE, HERE-1 ).EQ.ZERO )
    // LAPACK❗✔️:      $            NBNEXT = 1
    // LAPACK❗✔️:                IF( NBNEXT.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 2 by 2 Block did not split
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1,
    // LAPACK❗✔️:      $                         WORK, INFO )
    // LAPACK❗✔️:                   IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:                      ILST = HERE
    // LAPACK❗✔️:                      RETURN
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                   HERE = HERE - 2
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                 2 by 2 Block did split
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1,
    // LAPACK❗✔️:      $                         WORK, INFO )
    // LAPACK❗✔️:                   CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1,
    // LAPACK❗✔️:      $                         WORK, INFO )
    // LAPACK❗✔️:                   HERE = HERE - 2
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          IF( HERE.GT.ILST )
    // LAPACK❗✔️:      $      GO TO 20
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       ILST = HERE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DTREXC
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    if n <= 1 {
        return 0;
    }
    if *ifst > 1 && t[ix(*ifst, *ifst - 1, ldt)] != 0.0 {
        *ifst -= 1;
    }
    let mut nbf = if *ifst < n && t[ix(*ifst + 1, *ifst, ldt)] != 0.0 {
        2
    } else {
        1
    };
    if *ilst > 1 && t[ix(*ilst, *ilst - 1, ldt)] != 0.0 {
        *ilst -= 1;
    }
    let nbl = if *ilst < n && t[ix(*ilst + 1, *ilst, ldt)] != 0.0 {
        2
    } else {
        1
    };
    if *ifst == *ilst {
        return 0;
    }
    let mut here = *ifst;
    if *ifst < *ilst {
        if nbf == 2 && nbl == 1 {
            *ilst -= 1;
        }
        if nbf == 1 && nbl == 2 {
            *ilst += 1;
        }
        loop {
            if nbf == 1 || nbf == 2 {
                let nbnext = if here + nbf + 1 <= n && t[ix(here + nbf + 1, here + nbf, ldt)] != 0.0
                {
                    2
                } else {
                    1
                };
                let info = laexc(true, n, t, ldt, q, ldq, here, nbf, nbnext);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                here += nbnext;
                if nbf == 2 && t[ix(here + 1, here, ldt)] == 0.0 {
                    nbf = 3;
                }
            } else {
                let mut nbnext = if here + 3 <= n && t[ix(here + 3, here + 2, ldt)] != 0.0 {
                    2
                } else {
                    1
                };
                let info = laexc(true, n, t, ldt, q, ldq, here + 1, 1, nbnext);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                if nbnext == 1 {
                    laexc(true, n, t, ldt, q, ldq, here, 1, 1);
                    here += 1;
                } else {
                    if t[ix(here + 2, here + 1, ldt)] == 0.0 {
                        nbnext = 1;
                    }
                    if nbnext == 2 {
                        let info = laexc(true, n, t, ldt, q, ldq, here, 1, 2);
                        if info != 0 {
                            *ilst = here;
                            return info;
                        }
                        here += 2;
                    } else {
                        laexc(true, n, t, ldt, q, ldq, here, 1, 1);
                        laexc(true, n, t, ldt, q, ldq, here + 1, 1, 1);
                        here += 2;
                    }
                }
            }
            if here >= *ilst {
                break;
            }
        }
    } else {
        loop {
            if nbf == 1 || nbf == 2 {
                let nbnext = if here >= 3 && t[ix(here - 1, here - 2, ldt)] != 0.0 {
                    2
                } else {
                    1
                };
                let info = laexc(true, n, t, ldt, q, ldq, here - nbnext, nbnext, nbf);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                here -= nbnext;
                if nbf == 2 && t[ix(here + 1, here, ldt)] == 0.0 {
                    nbf = 3;
                }
            } else {
                let mut nbnext = if here >= 3 && t[ix(here - 1, here - 2, ldt)] != 0.0 {
                    2
                } else {
                    1
                };
                let info = laexc(true, n, t, ldt, q, ldq, here - nbnext, nbnext, 1);
                if info != 0 {
                    *ilst = here;
                    return info;
                }
                if nbnext == 1 {
                    laexc(true, n, t, ldt, q, ldq, here, 1, 1);
                    here -= 1;
                } else {
                    if t[ix(here, here - 1, ldt)] == 0.0 {
                        nbnext = 1;
                    }
                    if nbnext == 2 {
                        let info = laexc(true, n, t, ldt, q, ldq, here - 1, 2, 1);
                        if info != 0 {
                            *ilst = here;
                            return info;
                        }
                        here -= 2;
                    } else {
                        laexc(true, n, t, ldt, q, ldq, here, 1, 1);
                        laexc(true, n, t, ldt, q, ldq, here - 1, 1, 1);
                        here -= 2;
                    }
                }
            }
            if here <= *ilst {
                break;
            }
        }
    }
    *ilst = here;
    0
}
