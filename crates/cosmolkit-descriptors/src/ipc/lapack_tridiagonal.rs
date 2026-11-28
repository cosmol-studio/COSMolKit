//! Eigenvalue-only tridiagonal LAPACK kernels from the target NumPy backend.
//!
//! OpenBLAS 446c436e10450c348169808f1a6b3fae0925c9f7, lapack-netlib/SRC.
//! LAPACK copyright: University of Tennessee, University of California Berkeley,
//! University of Colorado Denver, and NAG Ltd.; BSD-3-Clause.
//! The arithmetic follows that source. The target wheel's DSTERF and scalar
//! helpers contain separate SSE2 operations, without FMA contraction. Complete
//! numerical equivalence remains unverified; no external BLAS is used.

// LAPACK license notice (the source distribution retains its original placeholders):
// Copyright (c) 1992-2023 The University of Tennessee and The University
//                         of Tennessee Research Foundation.  All rights
//                         reserved.
// Copyright (c) 2000-2023 The University of California Berkeley. All
//                         rights reserved.
// Copyright (c) 2006-2023 The University of Colorado Denver.  All rights
//                         reserved.
//
// $COPYRIGHT$
//
// Additional copyrights may follow
//
// $HEADER$
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// - Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//
// - Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer listed
//   in this license in the documentation and/or other materials
//   provided with the distribution.
//
// - Neither the name of the copyright holders nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// The copyright holders provide no reassurances that the source code
// provided does not infringe any patent, copyright, or any other
// intellectual property rights of third parties.  The copyright holders
// disclaim any liability to any recipient for claims brought against
// recipient by any third party for infringement of that parties
// intellectual property rights.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum DsterfError {
    InvalidOffDiagonalLength { expected: usize, actual: usize },
    DimensionOverflow,
    InvalidArgument { routine: &'static str, argument: u8 },
    SortingStackExhausted,
    NonConvergence { remaining: usize },
}
impl fmt::Display for DsterfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOffDiagonalLength { expected, actual } => write!(
                f,
                "DSTERF requires at least {expected} off-diagonal entries, got {actual}"
            ),
            Self::DimensionOverflow => {
                f.write_str("DSTERF dimension exceeds the iteration index range")
            }
            Self::InvalidArgument { routine, argument } => {
                write!(f, "{routine} argument {argument} is invalid")
            }
            Self::SortingStackExhausted => f.write_str("DLASRT fixed sorting stack exhausted"),
            Self::NonConvergence { remaining } => write!(
                f,
                "DSTERF did not converge: {remaining} off-diagonal entries remain"
            ),
        }
    }
}
impl std::error::Error for DsterfError {}

/// DLAMCH on the target IEEE binary64 environment; control characters are ASCII.
pub(super) fn dlamch(cmach: u8) -> f64 {
    // LAPACK❗✔️:       DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine (version 3.7.0) --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *     December 2016
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          CMACH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: * =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE, ZERO
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, HUGEVAL, RMACH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT,
    // LAPACK❗✔️:      $                   MINEXPONENT, RADIX, TINY
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Assume rounding, not chopping. Always.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RND = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ONE.EQ.RND ) THEN
    // LAPACK❗✔️:          EPS = EPSILON(ZERO) * 0.5
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          EPS = EPSILON(ZERO)
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( CMACH, 'E' ) ) THEN
    // LAPACK❗✔️:          RMACH = EPS
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'S' ) ) THEN
    // LAPACK❗✔️:          SFMIN = TINY(ZERO)
    // LAPACK❗✔️:          HUGEVAL = HUGE(ZERO)
    // LAPACK❗✔️:          IF( HUGEVAL*SFMIN.LE.ONE ) THEN
    // LAPACK❗✔️:             SMALL = ONE / HUGEVAL
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Use SMALL plus a bit, to avoid the possibility of rounding
    // LAPACK❗✔️: *           causing overflow when computing  1/sfmin.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             SFMIN = SMALL*( ONE+EPS )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          RMACH = SFMIN
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'B' ) ) THEN
    // LAPACK❗✔️:          RMACH = RADIX(ZERO)
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'P' ) ) THEN
    // LAPACK❗✔️:          RMACH = EPS * RADIX(ZERO)
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'N' ) ) THEN
    // LAPACK❗✔️:          RMACH = DIGITS(ZERO)
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'R' ) ) THEN
    // LAPACK❗✔️:          RMACH = RND
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'M' ) ) THEN
    // LAPACK❗✔️:          RMACH = MINEXPONENT(ZERO)
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'U' ) ) THEN
    // LAPACK❗✔️:          RMACH = TINY(ZERO)
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'L' ) ) THEN
    // LAPACK❗✔️:          RMACH = MAXEXPONENT(ZERO)
    // LAPACK❗✔️:       ELSE IF( LSAME( CMACH, 'O' ) ) THEN
    // LAPACK❗✔️:          RMACH = HUGE(ZERO)
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          RMACH = ZERO
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DLAMCH = RMACH
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLAMCH
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    let eps = f64::EPSILON * 0.5;
    match cmach.to_ascii_uppercase() {
        b'E' => eps,
        b'S' => {
            let mut sfmin = f64::MIN_POSITIVE;
            let hugeval = f64::MAX;
            if hugeval * sfmin <= 1.0 {
                let small = 1.0 / hugeval;
                sfmin = small * (1.0 + eps);
            }
            sfmin
        }
        b'B' => f64::RADIX as f64,
        b'P' => eps * f64::RADIX as f64,
        b'N' => f64::MANTISSA_DIGITS as f64,
        b'R' => 1.0,
        b'M' => f64::MIN_EXP as f64,
        b'U' => f64::MIN_POSITIVE,
        b'L' => f64::MAX_EXP as f64,
        b'O' => f64::MAX,
        _ => 0.0,
    }
}

/// DLANST('M'): the source visits the final diagonal first, then each d/e pair.
fn dlanst_max(d: &[f64], e: &[f64]) -> f64 {
    // LAPACK❗✔️:       DOUBLE PRECISION FUNCTION DLANST( NORM, N, D, E )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          NORM
    // LAPACK❗✔️:       INTEGER            N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   D( * ), E( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE, ZERO
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I
    // LAPACK❗✔️:       DOUBLE PRECISION   ANORM, SCALE, SUM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME, DISNAN
    // LAPACK❗✔️:       EXTERNAL           LSAME, DISNAN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLASSQ
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.0 ) THEN
    // LAPACK❗✔️:          ANORM = ZERO
    // LAPACK❗✔️:       ELSE IF( LSAME( NORM, 'M' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Find max(abs(A(i,j))).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ANORM = ABS( D( N ) )
    // LAPACK❗✔️:          DO 10 I = 1, N - 1
    // LAPACK❗✔️:             SUM = ABS( D( I ) )
    // LAPACK❗✔️:             IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
    // LAPACK❗✔️:             SUM = ABS( E( I ) )
    // LAPACK❗✔️:             IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️:       ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
    // LAPACK❗✔️:      $         LSAME( NORM, 'I' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Find norm1(A).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( N.EQ.1 ) THEN
    // LAPACK❗✔️:             ANORM = ABS( D( 1 ) )
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             ANORM = ABS( D( 1 ) )+ABS( E( 1 ) )
    // LAPACK❗✔️:             SUM = ABS( E( N-1 ) )+ABS( D( N ) )
    // LAPACK❗✔️:             IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
    // LAPACK❗✔️:             DO 20 I = 2, N - 1
    // LAPACK❗✔️:                SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) )
    // LAPACK❗✔️:                IF( ANORM .LT. SUM .OR. DISNAN( SUM ) ) ANORM = SUM
    // LAPACK❗✔️:    20       CONTINUE
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Find normF(A).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SCALE = ZERO
    // LAPACK❗✔️:          SUM = ONE
    // LAPACK❗✔️:          IF( N.GT.1 ) THEN
    // LAPACK❗✔️:             CALL DLASSQ( N-1, E, 1, SCALE, SUM )
    // LAPACK❗✔️:             SUM = 2*SUM
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          CALL DLASSQ( N, D, 1, SCALE, SUM )
    // LAPACK❗✔️:          ANORM = SCALE*SQRT( SUM )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DLANST = ANORM
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLANST
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    let Some(last) = d.last() else {
        return 0.0;
    };
    let mut anorm = last.abs();
    for i in 0..d.len() - 1 {
        let sum = d[i].abs();
        if anorm < sum || sum.is_nan() {
            anorm = sum;
        }
        let sum = e[i].abs();
        if anorm < sum || sum.is_nan() {
            anorm = sum;
        }
    }
    anorm
}

/// Source DLAPY2 expression; Rust's hypot may use different arithmetic.
pub(super) fn dlapy2(x: f64, y: f64) -> f64 {
    // LAPACK❗✔️:       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   X, Y
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   W, XABS, YABS, Z, HUGEVAL
    // LAPACK❗✔️:       LOGICAL            X_IS_NAN, Y_IS_NAN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            DISNAN
    // LAPACK❗✔️:       EXTERNAL           DISNAN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, MAX, MIN, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       X_IS_NAN = DISNAN( X )
    // LAPACK❗✔️:       Y_IS_NAN = DISNAN( Y )
    // LAPACK❗✔️:       IF ( X_IS_NAN ) DLAPY2 = X
    // LAPACK❗✔️:       IF ( Y_IS_NAN ) DLAPY2 = Y
    // LAPACK❗✔️:       HUGEVAL = DLAMCH( 'Overflow' )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
    // LAPACK❗✔️:          XABS = ABS( X )
    // LAPACK❗✔️:          YABS = ABS( Y )
    // LAPACK❗✔️:          W = MAX( XABS, YABS )
    // LAPACK❗✔️:          Z = MIN( XABS, YABS )
    // LAPACK❗✔️:          IF( Z.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
    // LAPACK❗✔️:             DLAPY2 = W
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLAPY2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if y.is_nan() {
        return y;
    }
    if x.is_nan() {
        return x;
    }
    let xabs = x.abs();
    let yabs = y.abs();
    let w = if xabs > yabs { xabs } else { yabs };
    let z = if xabs < yabs { xabs } else { yabs };
    if z == 0.0 || w > dlamch(b'O') {
        w
    } else {
        let ratio = z / w;
        w * (1.0 + ratio * ratio).sqrt()
    }
}

/// Eigenvalues of [[a,b],[b,c]], ordered by decreasing absolute value.
pub(super) fn dlae2(a: f64, b: f64, c: f64) -> (f64, f64) {
    // LAPACK❗✔️:       SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A, B, C, RT1, RT2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: * =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   TWO
    // LAPACK❗✔️:       PARAMETER          ( TWO = 2.0D0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   HALF
    // LAPACK❗✔️:       PARAMETER          ( HALF = 0.5D0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Compute the eigenvalues
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SM = A + C
    // LAPACK❗✔️:       DF = A - C
    // LAPACK❗✔️:       ADF = ABS( DF )
    // LAPACK❗✔️:       TB = B + B
    // LAPACK❗✔️:       AB = ABS( TB )
    // LAPACK❗✔️:       IF( ABS( A ).GT.ABS( C ) ) THEN
    // LAPACK❗✔️:          ACMX = A
    // LAPACK❗✔️:          ACMN = C
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          ACMX = C
    // LAPACK❗✔️:          ACMN = A
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( ADF.GT.AB ) THEN
    // LAPACK❗✔️:          RT = ADF*SQRT( ONE+( AB / ADF )**2 )
    // LAPACK❗✔️:       ELSE IF( ADF.LT.AB ) THEN
    // LAPACK❗✔️:          RT = AB*SQRT( ONE+( ADF / AB )**2 )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Includes case AB=ADF=0
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RT = AB*SQRT( TWO )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( SM.LT.ZERO ) THEN
    // LAPACK❗✔️:          RT1 = HALF*( SM-RT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Order of execution important.
    // LAPACK❗✔️: *        To get fully accurate smaller eigenvalue,
    // LAPACK❗✔️: *        next line needs to be executed in higher precision.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
    // LAPACK❗✔️:       ELSE IF( SM.GT.ZERO ) THEN
    // LAPACK❗✔️:          RT1 = HALF*( SM+RT )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Order of execution important.
    // LAPACK❗✔️: *        To get fully accurate smaller eigenvalue,
    // LAPACK❗✔️: *        next line needs to be executed in higher precision.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Includes case RT1 = RT2 = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RT1 = HALF*RT
    // LAPACK❗✔️:          RT2 = -HALF*RT
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLAE2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    let sm = a + c;
    let df = a - c;
    let adf = df.abs();
    let tb = b + b;
    let ab = tb.abs();
    let (acmx, acmn) = if a.abs() > c.abs() { (a, c) } else { (c, a) };
    let rt = if adf > ab {
        let ratio = ab / adf;
        adf * (1.0 + ratio * ratio).sqrt()
    } else if adf < ab {
        let ratio = adf / ab;
        ab * (1.0 + ratio * ratio).sqrt()
    } else {
        ab * 2.0_f64.sqrt()
    };
    if sm < 0.0 {
        let rt1 = 0.5 * (sm - rt);
        let rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
        (rt1, rt2)
    } else if sm > 0.0 {
        let rt1 = 0.5 * (sm + rt);
        let rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
        (rt1, rt2)
    } else {
        (0.5 * rt, -0.5 * rt)
    }
}

/// DLASCL('G') for a contiguous column/vector, preserving each scaling pass.
/// Callers control the matrix shape; this helper never visits padding/triangles.
pub(super) fn dlascl_vector(cfrom: f64, cto: f64, values: &mut [f64]) -> Result<(), DsterfError> {
    // LAPACK❗✔️:       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA,
    // LAPACK❗✔️:      $                   INFO )
    // LAPACK❗✔️:       IMPLICIT NONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          TYPE
    // LAPACK❗✔️:       INTEGER            INFO, KL, KU, LDA, M, N
    // LAPACK❗✔️:       DOUBLE PRECISION   CFROM, CTO
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            DONE
    // LAPACK❗✔️:       INTEGER            I, ITYPE, J, K1, K2, K3, K4
    // LAPACK❗✔️:       DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME, DISNAN
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH
    // LAPACK❗✔️:       EXTERNAL           LSAME, DLAMCH, DISNAN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input arguments
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( TYPE, 'G' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 0
    // LAPACK❗✔️:       ELSE IF( LSAME( TYPE, 'L' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 1
    // LAPACK❗✔️:       ELSE IF( LSAME( TYPE, 'U' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 2
    // LAPACK❗✔️:       ELSE IF( LSAME( TYPE, 'H' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 3
    // LAPACK❗✔️:       ELSE IF( LSAME( TYPE, 'B' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 4
    // LAPACK❗✔️:       ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 5
    // LAPACK❗✔️:       ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
    // LAPACK❗✔️:          ITYPE = 6
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          ITYPE = -1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ITYPE.EQ.-1 ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( DISNAN(CTO) ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( M.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -6
    // LAPACK❗✔️:       ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
    // LAPACK❗✔️:      $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
    // LAPACK❗✔️:          INFO = -7
    // LAPACK❗✔️:       ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
    // LAPACK❗✔️:          INFO = -9
    // LAPACK❗✔️:       ELSE IF( ITYPE.GE.4 ) THEN
    // LAPACK❗✔️:          IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
    // LAPACK❗✔️:             INFO = -2
    // LAPACK❗✔️:          ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
    // LAPACK❗✔️:      $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
    // LAPACK❗✔️:      $             THEN
    // LAPACK❗✔️:             INFO = -3
    // LAPACK❗✔️:          ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
    // LAPACK❗✔️:      $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
    // LAPACK❗✔️:      $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
    // LAPACK❗✔️:             INFO = -9
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DLASCL', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 .OR. M.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Get machine parameters
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SMLNUM = DLAMCH( 'S' )
    // LAPACK❗✔️:       BIGNUM = ONE / SMLNUM
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CFROMC = CFROM
    // LAPACK❗✔️:       CTOC = CTO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       CFROM1 = CFROMC*SMLNUM
    // LAPACK❗✔️:       IF( CFROM1.EQ.CFROMC ) THEN
    // LAPACK❗✔️: !        CFROMC is an inf.  Multiply by a correctly signed zero for
    // LAPACK❗✔️: !        finite CTOC, or a NaN if CTOC is infinite.
    // LAPACK❗✔️:          MUL = CTOC / CFROMC
    // LAPACK❗✔️:          DONE = .TRUE.
    // LAPACK❗✔️:          CTO1 = CTOC
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          CTO1 = CTOC / BIGNUM
    // LAPACK❗✔️:          IF( CTO1.EQ.CTOC ) THEN
    // LAPACK❗✔️: !           CTOC is either 0 or an inf.  In both cases, CTOC itself
    // LAPACK❗✔️: !           serves as the correct multiplication factor.
    // LAPACK❗✔️:             MUL = CTOC
    // LAPACK❗✔️:             DONE = .TRUE.
    // LAPACK❗✔️:             CFROMC = ONE
    // LAPACK❗✔️:          ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
    // LAPACK❗✔️:             MUL = SMLNUM
    // LAPACK❗✔️:             DONE = .FALSE.
    // LAPACK❗✔️:             CFROMC = CFROM1
    // LAPACK❗✔️:          ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
    // LAPACK❗✔️:             MUL = BIGNUM
    // LAPACK❗✔️:             DONE = .FALSE.
    // LAPACK❗✔️:             CTOC = CTO1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             MUL = CTOC / CFROMC
    // LAPACK❗✔️:             DONE = .TRUE.
    // LAPACK❗✔️:             IF (MUL .EQ. ONE)
    // LAPACK❗✔️:      $         RETURN
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ITYPE.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Full matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 30 J = 1, N
    // LAPACK❗✔️:             DO 20 I = 1, M
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:    20       CONTINUE
    // LAPACK❗✔️:    30    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ITYPE.EQ.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Lower triangular matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 50 J = 1, MIN( M, N )
    // LAPACK❗✔️:             DO 40 I = J, M
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:    40       CONTINUE
    // LAPACK❗✔️:    50    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ITYPE.EQ.2 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Upper triangular matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 70 J = 1, N
    // LAPACK❗✔️:             DO 60 I = 1, MIN( J, M )
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:    60       CONTINUE
    // LAPACK❗✔️:    70    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ITYPE.EQ.3 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Upper Hessenberg matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 90 J = 1, N
    // LAPACK❗✔️:             DO 80 I = 1, MIN( J+1, M )
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:    80       CONTINUE
    // LAPACK❗✔️:    90    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ITYPE.EQ.4 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Lower half of a symmetric band matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          K3 = KL + 1
    // LAPACK❗✔️:          K4 = N + 1
    // LAPACK❗✔️:          DO 110 J = 1, N
    // LAPACK❗✔️:             DO 100 I = 1, MIN( K3, K4-J )
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:   100       CONTINUE
    // LAPACK❗✔️:   110    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ITYPE.EQ.5 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Upper half of a symmetric band matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          K1 = KU + 2
    // LAPACK❗✔️:          K3 = KU + 1
    // LAPACK❗✔️:          DO 130 J = 1, N
    // LAPACK❗✔️:             DO 120 I = MAX( K1-J, 1 ), K3
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:   120       CONTINUE
    // LAPACK❗✔️:   130    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ITYPE.EQ.6 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Band matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          K1 = KL + KU + 2
    // LAPACK❗✔️:          K2 = KL + 1
    // LAPACK❗✔️:          K3 = 2*KL + KU + 1
    // LAPACK❗✔️:          K4 = KL + KU + 1 + M
    // LAPACK❗✔️:          DO 150 J = 1, N
    // LAPACK❗✔️:             DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
    // LAPACK❗✔️:                A( I, J ) = A( I, J )*MUL
    // LAPACK❗✔️:   140       CONTINUE
    // LAPACK❗✔️:   150    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( .NOT.DONE )
    // LAPACK❗✔️:      $   GO TO 10
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLASCL
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if cfrom == 0.0 || cfrom.is_nan() {
        return Err(DsterfError::InvalidArgument {
            routine: "DLASCL",
            argument: 4,
        });
    }
    if cto.is_nan() {
        return Err(DsterfError::InvalidArgument {
            routine: "DLASCL",
            argument: 5,
        });
    }
    if values.is_empty() {
        return Ok(());
    }
    let smlnum = dlamch(b'S');
    let bignum = 1.0 / smlnum;
    let mut cfromc = cfrom;
    let mut ctoc = cto;
    loop {
        let cfrom1 = cfromc * smlnum;
        let (mul, done);
        if cfrom1 == cfromc {
            mul = ctoc / cfromc;
            done = true;
        } else {
            let cto1 = ctoc / bignum;
            if cto1 == ctoc {
                mul = ctoc;
                done = true;
                cfromc = 1.0;
            } else if cfrom1.abs() > ctoc.abs() && ctoc != 0.0 {
                mul = smlnum;
                done = false;
                cfromc = cfrom1;
            } else if cto1.abs() > cfromc.abs() {
                mul = bignum;
                done = false;
                ctoc = cto1;
            } else {
                mul = ctoc / cfromc;
                done = true;
                if mul == 1.0 {
                    return Ok(());
                }
            }
        }
        for value in values.iter_mut() {
            *value *= mul;
        }
        if done {
            return Ok(());
        }
    }
}

/// DLASRT('I'), including the source insertion cutoff and median-of-three swaps.
fn dlasrt_increasing(d: &mut [f64]) -> Result<(), DsterfError> {
    // LAPACK❗✔️:       SUBROUTINE DLASRT( ID, N, D, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          ID
    // LAPACK❗✔️:       INTEGER            INFO, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   D( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       INTEGER            SELECT
    // LAPACK❗✔️:       PARAMETER          ( SELECT = 20 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            DIR, ENDD, I, J, START, STKPNT
    // LAPACK❗✔️:       DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       INTEGER            STACK( 2, 32 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       EXTERNAL           LSAME
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       DIR = -1
    // LAPACK❗✔️:       IF( LSAME( ID, 'D' ) ) THEN
    // LAPACK❗✔️:          DIR = 0
    // LAPACK❗✔️:       ELSE IF( LSAME( ID, 'I' ) ) THEN
    // LAPACK❗✔️:          DIR = 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( DIR.EQ.-1 ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DLASRT', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.1 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       STKPNT = 1
    // LAPACK❗✔️:       STACK( 1, 1 ) = 1
    // LAPACK❗✔️:       STACK( 2, 1 ) = N
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       START = STACK( 1, STKPNT )
    // LAPACK❗✔️:       ENDD = STACK( 2, STKPNT )
    // LAPACK❗✔️:       STKPNT = STKPNT - 1
    // LAPACK❗✔️:       IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Do Insertion sort on D( START:ENDD )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( DIR.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Sort into decreasing order
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 30 I = START + 1, ENDD
    // LAPACK❗✔️:                DO 20 J = I, START + 1, -1
    // LAPACK❗✔️:                   IF( D( J ).GT.D( J-1 ) ) THEN
    // LAPACK❗✔️:                      DMNMX = D( J )
    // LAPACK❗✔️:                      D( J ) = D( J-1 )
    // LAPACK❗✔️:                      D( J-1 ) = DMNMX
    // LAPACK❗✔️:                   ELSE
    // LAPACK❗✔️:                      GO TO 30
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:    20          CONTINUE
    // LAPACK❗✔️:    30       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Sort into increasing order
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 50 I = START + 1, ENDD
    // LAPACK❗✔️:                DO 40 J = I, START + 1, -1
    // LAPACK❗✔️:                   IF( D( J ).LT.D( J-1 ) ) THEN
    // LAPACK❗✔️:                      DMNMX = D( J )
    // LAPACK❗✔️:                      D( J ) = D( J-1 )
    // LAPACK❗✔️:                      D( J-1 ) = DMNMX
    // LAPACK❗✔️:                   ELSE
    // LAPACK❗✔️:                      GO TO 50
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:    40          CONTINUE
    // LAPACK❗✔️:    50       CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ENDD-START.GT.SELECT ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Partition D( START:ENDD ) and stack parts, largest one first
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Choose partition entry as median of 3
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          D1 = D( START )
    // LAPACK❗✔️:          D2 = D( ENDD )
    // LAPACK❗✔️:          I = ( START+ENDD ) / 2
    // LAPACK❗✔️:          D3 = D( I )
    // LAPACK❗✔️:          IF( D1.LT.D2 ) THEN
    // LAPACK❗✔️:             IF( D3.LT.D1 ) THEN
    // LAPACK❗✔️:                DMNMX = D1
    // LAPACK❗✔️:             ELSE IF( D3.LT.D2 ) THEN
    // LAPACK❗✔️:                DMNMX = D3
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                DMNMX = D2
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             IF( D3.LT.D2 ) THEN
    // LAPACK❗✔️:                DMNMX = D2
    // LAPACK❗✔️:             ELSE IF( D3.LT.D1 ) THEN
    // LAPACK❗✔️:                DMNMX = D3
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                DMNMX = D1
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( DIR.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Sort into decreasing order
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             I = START - 1
    // LAPACK❗✔️:             J = ENDD + 1
    // LAPACK❗✔️:    60       CONTINUE
    // LAPACK❗✔️:    70       CONTINUE
    // LAPACK❗✔️:             J = J - 1
    // LAPACK❗✔️:             IF( D( J ).LT.DMNMX )
    // LAPACK❗✔️:      $         GO TO 70
    // LAPACK❗✔️:    80       CONTINUE
    // LAPACK❗✔️:             I = I + 1
    // LAPACK❗✔️:             IF( D( I ).GT.DMNMX )
    // LAPACK❗✔️:      $         GO TO 80
    // LAPACK❗✔️:             IF( I.LT.J ) THEN
    // LAPACK❗✔️:                TMP = D( I )
    // LAPACK❗✔️:                D( I ) = D( J )
    // LAPACK❗✔️:                D( J ) = TMP
    // LAPACK❗✔️:                GO TO 60
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( J-START.GT.ENDD-J-1 ) THEN
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = START
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = J
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = J + 1
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = ENDD
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = J + 1
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = ENDD
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = START
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = J
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Sort into increasing order
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             I = START - 1
    // LAPACK❗✔️:             J = ENDD + 1
    // LAPACK❗✔️:    90       CONTINUE
    // LAPACK❗✔️:   100       CONTINUE
    // LAPACK❗✔️:             J = J - 1
    // LAPACK❗✔️:             IF( D( J ).GT.DMNMX )
    // LAPACK❗✔️:      $         GO TO 100
    // LAPACK❗✔️:   110       CONTINUE
    // LAPACK❗✔️:             I = I + 1
    // LAPACK❗✔️:             IF( D( I ).LT.DMNMX )
    // LAPACK❗✔️:      $         GO TO 110
    // LAPACK❗✔️:             IF( I.LT.J ) THEN
    // LAPACK❗✔️:                TMP = D( I )
    // LAPACK❗✔️:                D( I ) = D( J )
    // LAPACK❗✔️:                D( J ) = TMP
    // LAPACK❗✔️:                GO TO 90
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( J-START.GT.ENDD-J-1 ) THEN
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = START
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = J
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = J + 1
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = ENDD
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = J + 1
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = ENDD
    // LAPACK❗✔️:                STKPNT = STKPNT + 1
    // LAPACK❗✔️:                STACK( 1, STKPNT ) = START
    // LAPACK❗✔️:                STACK( 2, STKPNT ) = J
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( STKPNT.GT.0 )
    // LAPACK❗✔️:      $   GO TO 10
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLASRT
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if d.len() <= 1 {
        return Ok(());
    }
    let end = isize::try_from(d.len() - 1).map_err(|_| DsterfError::DimensionOverflow)?;
    let mut stack = [(0_isize, 0_isize); 32];
    stack[0] = (0, end);
    let mut stkpnt = 1;
    while stkpnt > 0 {
        stkpnt -= 1;
        let (start, endd) = stack[stkpnt];
        if endd - start <= 20 && endd - start > 0 {
            for i in start + 1..=endd {
                let mut j = i;
                while j > start {
                    if d[j as usize] < d[(j - 1) as usize] {
                        d.swap(j as usize, (j - 1) as usize);
                    } else {
                        break;
                    }
                    j -= 1;
                }
            }
        } else if endd - start > 20 {
            let d1 = d[start as usize];
            let d2 = d[endd as usize];
            let d3 = d[((start + endd) / 2) as usize];
            let dmnmx = if d1 < d2 {
                if d3 < d1 {
                    d1
                } else if d3 < d2 {
                    d3
                } else {
                    d2
                }
            } else if d3 < d2 {
                d2
            } else if d3 < d1 {
                d3
            } else {
                d1
            };
            let mut i = start - 1;
            let mut j = endd + 1;
            loop {
                j -= 1;
                while d[j as usize] > dmnmx {
                    j -= 1;
                }
                i += 1;
                while d[i as usize] < dmnmx {
                    i += 1;
                }
                if i < j {
                    d.swap(i as usize, j as usize);
                } else {
                    break;
                }
            }
            if stkpnt + 2 > stack.len() {
                return Err(DsterfError::SortingStackExhausted);
            }
            if j - start > endd - j - 1 {
                stack[stkpnt] = (start, j);
                stack[stkpnt + 1] = (j + 1, endd);
            } else {
                stack[stkpnt] = (j + 1, endd);
                stack[stkpnt + 1] = (start, j);
            }
            stkpnt += 2;
        }
    }
    Ok(())
}

/// Pal-Walker-Kahan QL/QR eigenvalues, in place; `e` is destroyed as in LAPACK.
pub(super) fn dsterf(d: &mut [f64], e: &mut [f64]) -> Result<(), DsterfError> {
    // LAPACK❗✔️:       SUBROUTINE DSTERF( N, D, E, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            INFO, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   D( * ), E( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE, TWO, THREE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
    // LAPACK❗✔️:      $                   THREE = 3.0D0 )
    // LAPACK❗✔️:       INTEGER            MAXIT
    // LAPACK❗✔️:       PARAMETER          ( MAXIT = 30 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M,
    // LAPACK❗✔️:      $                   NMAXIT
    // LAPACK❗✔️:       DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,
    // LAPACK❗✔️:      $                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN,
    // LAPACK❗✔️:      $                   SIGMA, SSFMAX, SSFMIN, RMAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
    // LAPACK❗✔️:       EXTERNAL           DLAMCH, DLANST, DLAPY2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, SIGN, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:          CALL XERBLA( 'DSTERF', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( N.LE.1 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Determine the unit roundoff for this environment.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       EPS = DLAMCH( 'E' )
    // LAPACK❗✔️:       EPS2 = EPS**2
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'S' )
    // LAPACK❗✔️:       SAFMAX = ONE / SAFMIN
    // LAPACK❗✔️:       SSFMAX = SQRT( SAFMAX ) / THREE
    // LAPACK❗✔️:       SSFMIN = SQRT( SAFMIN ) / EPS2
    // LAPACK❗✔️:       RMAX = DLAMCH( 'O' )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Compute the eigenvalues of the tridiagonal matrix.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NMAXIT = N*MAXIT
    // LAPACK❗✔️:       SIGMA = ZERO
    // LAPACK❗✔️:       JTOT = 0
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Determine where the matrix splits and choose QL or QR iteration
    // LAPACK❗✔️: *     for each block, according to whether top or bottom diagonal
    // LAPACK❗✔️: *     element is smaller.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       L1 = 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10 CONTINUE
    // LAPACK❗✔️:       IF( L1.GT.N )
    // LAPACK❗✔️:      $   GO TO 170
    // LAPACK❗✔️:       IF( L1.GT.1 )
    // LAPACK❗✔️:      $   E( L1-1 ) = ZERO
    // LAPACK❗✔️:       DO 20 M = L1, N - 1
    // LAPACK❗✔️:          IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
    // LAPACK❗✔️:      $       1 ) ) ) )*EPS ) THEN
    // LAPACK❗✔️:             E( M ) = ZERO
    // LAPACK❗✔️:             GO TO 30
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:    20 CONTINUE
    // LAPACK❗✔️:       M = N
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    30 CONTINUE
    // LAPACK❗✔️:       L = L1
    // LAPACK❗✔️:       LSV = L
    // LAPACK❗✔️:       LEND = M
    // LAPACK❗✔️:       LENDSV = LEND
    // LAPACK❗✔️:       L1 = M + 1
    // LAPACK❗✔️:       IF( LEND.EQ.L )
    // LAPACK❗✔️:      $   GO TO 10
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Scale submatrix in rows and columns L to LEND
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ANORM = DLANST( 'M', LEND-L+1, D( L ), E( L ) )
    // LAPACK❗✔️:       ISCALE = 0
    // LAPACK❗✔️:       IF( ANORM.EQ.ZERO )
    // LAPACK❗✔️:      $   GO TO 10
    // LAPACK❗✔️:       IF( (ANORM.GT.SSFMAX) ) THEN
    // LAPACK❗✔️:          ISCALE = 1
    // LAPACK❗✔️:          CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
    // LAPACK❗✔️:      $                INFO )
    // LAPACK❗✔️:          CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
    // LAPACK❗✔️:      $                INFO )
    // LAPACK❗✔️:       ELSE IF( ANORM.LT.SSFMIN ) THEN
    // LAPACK❗✔️:          ISCALE = 2
    // LAPACK❗✔️:          CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
    // LAPACK❗✔️:      $                INFO )
    // LAPACK❗✔️:          CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
    // LAPACK❗✔️:      $                INFO )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO 40 I = L, LEND - 1
    // LAPACK❗✔️:          E( I ) = E( I )**2
    // LAPACK❗✔️:    40 CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Choose between QL and QR iteration
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
    // LAPACK❗✔️:          LEND = LSV
    // LAPACK❗✔️:          L = LENDSV
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LEND.GE.L ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        QL Iteration
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Look for small subdiagonal element.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    50    CONTINUE
    // LAPACK❗✔️:          IF( L.NE.LEND ) THEN
    // LAPACK❗✔️:             DO 60 M = L, LEND - 1
    // LAPACK❗✔️:                IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) )
    // LAPACK❗✔️:      $            GO TO 70
    // LAPACK❗✔️:    60       CONTINUE
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          M = LEND
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    70    CONTINUE
    // LAPACK❗✔️:          IF( M.LT.LEND )
    // LAPACK❗✔️:      $      E( M ) = ZERO
    // LAPACK❗✔️:          P = D( L )
    // LAPACK❗✔️:          IF( M.EQ.L )
    // LAPACK❗✔️:      $      GO TO 90
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        If remaining matrix is 2 by 2, use DLAE2 to compute its
    // LAPACK❗✔️: *        eigenvalues.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( M.EQ.L+1 ) THEN
    // LAPACK❗✔️:             RTE = SQRT( E( L ) )
    // LAPACK❗✔️:             CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
    // LAPACK❗✔️:             D( L ) = RT1
    // LAPACK❗✔️:             D( L+1 ) = RT2
    // LAPACK❗✔️:             E( L ) = ZERO
    // LAPACK❗✔️:             L = L + 2
    // LAPACK❗✔️:             IF( L.LE.LEND )
    // LAPACK❗✔️:      $         GO TO 50
    // LAPACK❗✔️:             GO TO 150
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( JTOT.EQ.NMAXIT )
    // LAPACK❗✔️:      $      GO TO 150
    // LAPACK❗✔️:          JTOT = JTOT + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Form shift.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RTE = SQRT( E( L ) )
    // LAPACK❗✔️:          SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
    // LAPACK❗✔️:          R = DLAPY2( SIGMA, ONE )
    // LAPACK❗✔️:          SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          C = ONE
    // LAPACK❗✔️:          S = ZERO
    // LAPACK❗✔️:          GAMMA = D( M ) - SIGMA
    // LAPACK❗✔️:          P = GAMMA*GAMMA
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Inner loop
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 80 I = M - 1, L, -1
    // LAPACK❗✔️:             BB = E( I )
    // LAPACK❗✔️:             R = P + BB
    // LAPACK❗✔️:             IF( I.NE.M-1 )
    // LAPACK❗✔️:      $         E( I+1 ) = S*R
    // LAPACK❗✔️:             OLDC = C
    // LAPACK❗✔️:             C = P / R
    // LAPACK❗✔️:             S = BB / R
    // LAPACK❗✔️:             OLDGAM = GAMMA
    // LAPACK❗✔️:             ALPHA = D( I )
    // LAPACK❗✔️:             GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
    // LAPACK❗✔️:             D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
    // LAPACK❗✔️:             IF( C.NE.ZERO ) THEN
    // LAPACK❗✔️:                P = ( GAMMA*GAMMA ) / C
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                P = OLDC*BB
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:    80    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          E( L ) = S*P
    // LAPACK❗✔️:          D( L ) = SIGMA + GAMMA
    // LAPACK❗✔️:          GO TO 50
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Eigenvalue found.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    90    CONTINUE
    // LAPACK❗✔️:          D( L ) = P
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          L = L + 1
    // LAPACK❗✔️:          IF( L.LE.LEND )
    // LAPACK❗✔️:      $      GO TO 50
    // LAPACK❗✔️:          GO TO 150
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        QR Iteration
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Look for small superdiagonal element.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   100    CONTINUE
    // LAPACK❗✔️:          DO 110 M = L, LEND + 1, -1
    // LAPACK❗✔️:             IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) )
    // LAPACK❗✔️:      $         GO TO 120
    // LAPACK❗✔️:   110    CONTINUE
    // LAPACK❗✔️:          M = LEND
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   120    CONTINUE
    // LAPACK❗✔️:          IF( M.GT.LEND )
    // LAPACK❗✔️:      $      E( M-1 ) = ZERO
    // LAPACK❗✔️:          P = D( L )
    // LAPACK❗✔️:          IF( M.EQ.L )
    // LAPACK❗✔️:      $      GO TO 140
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        If remaining matrix is 2 by 2, use DLAE2 to compute its
    // LAPACK❗✔️: *        eigenvalues.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( M.EQ.L-1 ) THEN
    // LAPACK❗✔️:             RTE = SQRT( E( L-1 ) )
    // LAPACK❗✔️:             CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
    // LAPACK❗✔️:             D( L ) = RT1
    // LAPACK❗✔️:             D( L-1 ) = RT2
    // LAPACK❗✔️:             E( L-1 ) = ZERO
    // LAPACK❗✔️:             L = L - 2
    // LAPACK❗✔️:             IF( L.GE.LEND )
    // LAPACK❗✔️:      $         GO TO 100
    // LAPACK❗✔️:             GO TO 150
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( JTOT.EQ.NMAXIT )
    // LAPACK❗✔️:      $      GO TO 150
    // LAPACK❗✔️:          JTOT = JTOT + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Form shift.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          RTE = SQRT( E( L-1 ) )
    // LAPACK❗✔️:          SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
    // LAPACK❗✔️:          R = DLAPY2( SIGMA, ONE )
    // LAPACK❗✔️:          SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          C = ONE
    // LAPACK❗✔️:          S = ZERO
    // LAPACK❗✔️:          GAMMA = D( M ) - SIGMA
    // LAPACK❗✔️:          P = GAMMA*GAMMA
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Inner loop
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 130 I = M, L - 1
    // LAPACK❗✔️:             BB = E( I )
    // LAPACK❗✔️:             R = P + BB
    // LAPACK❗✔️:             IF( I.NE.M )
    // LAPACK❗✔️:      $         E( I-1 ) = S*R
    // LAPACK❗✔️:             OLDC = C
    // LAPACK❗✔️:             C = P / R
    // LAPACK❗✔️:             S = BB / R
    // LAPACK❗✔️:             OLDGAM = GAMMA
    // LAPACK❗✔️:             ALPHA = D( I+1 )
    // LAPACK❗✔️:             GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
    // LAPACK❗✔️:             D( I ) = OLDGAM + ( ALPHA-GAMMA )
    // LAPACK❗✔️:             IF( C.NE.ZERO ) THEN
    // LAPACK❗✔️:                P = ( GAMMA*GAMMA ) / C
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                P = OLDC*BB
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:   130    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          E( L-1 ) = S*P
    // LAPACK❗✔️:          D( L ) = SIGMA + GAMMA
    // LAPACK❗✔️:          GO TO 100
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Eigenvalue found.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   140    CONTINUE
    // LAPACK❗✔️:          D( L ) = P
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          L = L - 1
    // LAPACK❗✔️:          IF( L.GE.LEND )
    // LAPACK❗✔️:      $      GO TO 100
    // LAPACK❗✔️:          GO TO 150
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Undo scaling if necessary
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   150 CONTINUE
    // LAPACK❗✔️:       IF( ISCALE.EQ.1 )
    // LAPACK❗✔️:      $   CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
    // LAPACK❗✔️:      $                D( LSV ), N, INFO )
    // LAPACK❗✔️:       IF( ISCALE.EQ.2 )
    // LAPACK❗✔️:      $   CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
    // LAPACK❗✔️:      $                D( LSV ), N, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Check for no convergence to an eigenvalue after a total
    // LAPACK❗✔️: *     of N*MAXIT iterations.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( JTOT.LT.NMAXIT )
    // LAPACK❗✔️:      $   GO TO 10
    // LAPACK❗✔️:       DO 160 I = 1, N - 1
    // LAPACK❗✔️:          IF( E( I ).NE.ZERO )
    // LAPACK❗✔️:      $      INFO = INFO + 1
    // LAPACK❗✔️:   160 CONTINUE
    // LAPACK❗✔️:       GO TO 180
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Sort eigenvalues in increasing order.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   170 CONTINUE
    // LAPACK❗✔️:       CALL DLASRT( 'I', N, D, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:   180 CONTINUE
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DSTERF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    let n = d.len();
    let expected = n.saturating_sub(1);
    if e.len() < expected {
        return Err(DsterfError::InvalidOffDiagonalLength {
            expected,
            actual: e.len(),
        });
    }
    if n <= 1 {
        return Ok(());
    }
    let nmaxit = n
        .checked_mul(30)
        .filter(|&count| i64::try_from(count).is_ok())
        .ok_or(DsterfError::DimensionOverflow)?;
    let eps = dlamch(b'E');
    let eps2 = eps * eps;
    let safmin = dlamch(b'S');
    let safmax = 1.0 / safmin;
    let ssfmax = safmax.sqrt() / 3.0;
    let ssfmin = safmin.sqrt() / eps2;
    // The source queries RMAX but never consumes it in this routine.
    let _rmax = dlamch(b'O');
    let mut jtot = 0;
    let mut l1 = 1;
    while l1 <= n {
        if l1 > 1 {
            e[l1 - 2] = 0.0;
        }
        let mut m = n;
        for k in l1..n {
            if e[k - 1].abs() <= (d[k - 1].abs().sqrt() * d[k].abs().sqrt()) * eps {
                e[k - 1] = 0.0;
                m = k;
                break;
            }
        }
        let mut l = l1;
        let lsv = l;
        let mut lend = m;
        let lendsv = lend;
        l1 = m + 1;
        if lend == l {
            continue;
        }
        let anorm = dlanst_max(&d[l - 1..lend], &e[l - 1..lend - 1]);
        let mut iscale = 0;
        if anorm == 0.0 {
            continue;
        }
        if anorm > ssfmax {
            iscale = 1;
            dlascl_vector(anorm, ssfmax, &mut d[l - 1..lend])?;
            dlascl_vector(anorm, ssfmax, &mut e[l - 1..lend - 1])?;
        } else if anorm < ssfmin {
            iscale = 2;
            dlascl_vector(anorm, ssfmin, &mut d[l - 1..lend])?;
            dlascl_vector(anorm, ssfmin, &mut e[l - 1..lend - 1])?;
        }
        for i in l..lend {
            e[i - 1] = e[i - 1] * e[i - 1];
        }
        if d[lend - 1].abs() < d[l - 1].abs() {
            lend = lsv;
            l = lendsv;
        }
        if lend >= l {
            loop {
                m = lend;
                if l != lend {
                    for k in l..lend {
                        if e[k - 1].abs() <= eps2 * (d[k - 1] * d[k]).abs() {
                            m = k;
                            break;
                        }
                    }
                }
                if m < lend {
                    e[m - 1] = 0.0;
                }
                let mut p = d[l - 1];
                if m == l {
                    d[l - 1] = p;
                    l += 1;
                    if l <= lend {
                        continue;
                    }
                    break;
                }
                if m == l + 1 {
                    let rte = e[l - 1].sqrt();
                    let (rt1, rt2) = dlae2(d[l - 1], rte, d[l]);
                    d[l - 1] = rt1;
                    d[l] = rt2;
                    e[l - 1] = 0.0;
                    l += 2;
                    if l <= lend {
                        continue;
                    }
                    break;
                }
                if jtot == nmaxit {
                    break;
                }
                jtot += 1;
                let rte = e[l - 1].sqrt();
                let mut sigma = (d[l] - p) / (2.0 * rte);
                let r = dlapy2(sigma, 1.0);
                sigma = p - (rte / (sigma + r.copysign(sigma)));
                let mut c = 1.0;
                let mut s = 0.0;
                let mut gamma = d[m - 1] - sigma;
                p = gamma * gamma;
                for i in (l..m).rev() {
                    let bb = e[i - 1];
                    let r = p + bb;
                    if i != m - 1 {
                        e[i] = s * r;
                    }
                    let oldc = c;
                    c = p / r;
                    s = bb / r;
                    let oldgam = gamma;
                    let alpha = d[i - 1];
                    gamma = c * (alpha - sigma) - s * oldgam;
                    d[i] = oldgam + (alpha - gamma);
                    p = if c != 0.0 {
                        (gamma * gamma) / c
                    } else {
                        oldc * bb
                    };
                }
                e[l - 1] = s * p;
                d[l - 1] = sigma + gamma;
            }
        } else {
            loop {
                m = lend;
                for k in (lend + 1..=l).rev() {
                    if e[k - 2].abs() <= eps2 * (d[k - 1] * d[k - 2]).abs() {
                        m = k;
                        break;
                    }
                }
                if m > lend {
                    e[m - 2] = 0.0;
                }
                let mut p = d[l - 1];
                if m == l {
                    d[l - 1] = p;
                    l -= 1;
                    if l >= lend {
                        continue;
                    }
                    break;
                }
                if m == l - 1 {
                    let rte = e[l - 2].sqrt();
                    let (rt1, rt2) = dlae2(d[l - 1], rte, d[l - 2]);
                    d[l - 1] = rt1;
                    d[l - 2] = rt2;
                    e[l - 2] = 0.0;
                    l -= 2;
                    if l >= lend {
                        continue;
                    }
                    break;
                }
                if jtot == nmaxit {
                    break;
                }
                jtot += 1;
                let rte = e[l - 2].sqrt();
                let mut sigma = (d[l - 2] - p) / (2.0 * rte);
                let r = dlapy2(sigma, 1.0);
                sigma = p - (rte / (sigma + r.copysign(sigma)));
                let mut c = 1.0;
                let mut s = 0.0;
                let mut gamma = d[m - 1] - sigma;
                p = gamma * gamma;
                for i in m..l {
                    let bb = e[i - 1];
                    let r = p + bb;
                    if i != m {
                        e[i - 2] = s * r;
                    }
                    let oldc = c;
                    c = p / r;
                    s = bb / r;
                    let oldgam = gamma;
                    let alpha = d[i];
                    gamma = c * (alpha - sigma) - s * oldgam;
                    d[i - 1] = oldgam + (alpha - gamma);
                    p = if c != 0.0 {
                        (gamma * gamma) / c
                    } else {
                        oldc * bb
                    };
                }
                e[l - 2] = s * p;
                d[l - 1] = sigma + gamma;
            }
        }
        if iscale == 1 {
            dlascl_vector(ssfmax, anorm, &mut d[lsv - 1..lendsv])?;
        }
        if iscale == 2 {
            dlascl_vector(ssfmin, anorm, &mut d[lsv - 1..lendsv])?;
        }
        if jtot == nmaxit {
            let remaining = e[..n - 1].iter().filter(|&&x| x != 0.0).count();
            return if remaining == 0 {
                Ok(())
            } else {
                Err(DsterfError::NonConvergence { remaining })
            };
        }
    }
    dlasrt_increasing(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dsterf_analytic_tridiagonal_repeated_and_split_roots() {
        for (mut d, mut e, expected) in [
            (vec![2., 2.], vec![1.], vec![1., 3.]),
            (vec![3., 1., 3., 1.], vec![0., 0., 0.], vec![1., 1., 3., 3.]),
            (
                vec![2., 2., 2.],
                vec![1., 1.],
                vec![2. - 2_f64.sqrt(), 2., 2. + 2_f64.sqrt()],
            ),
            (
                vec![4., 2., 0.],
                vec![1., 1.],
                vec![2. - 6_f64.sqrt(), 2., 2. + 6_f64.sqrt()],
            ),
        ] {
            dsterf(&mut d, &mut e).unwrap();
            for (actual, expected) in d.iter().zip(expected) {
                assert!(
                    (actual - expected).abs() < 2e-14,
                    "{d:?}, expected {expected}"
                );
            }
        }
        // A split diagonal of length22 enters the source quicksort partition,
        // rather than only its small insertion-sort branch.
        let mut diagonal: Vec<_> = (0..22).rev().map(f64::from).collect();
        dsterf(&mut diagonal, &mut [0.; 21]).unwrap();
        assert_eq!(diagonal, (0..22).map(f64::from).collect::<Vec<_>>());
        // Exercise both source scaling directions without absolute tolerances
        // that would accidentally accept zero for a tiny eigenvalue.
        for scale in [2.0_f64.powi(-900), 2.0_f64.powi(900)] {
            let mut diagonal = [2.0 * scale, 2.0 * scale];
            dsterf(&mut diagonal, &mut [scale]).unwrap();
            assert!((diagonal[0] / scale - 1.0).abs() < 2e-14);
            assert!((diagonal[1] / scale - 3.0).abs() < 2e-14);
        }
        let mut d = vec![2., 2., 2.];
        let error = dsterf(&mut d, &mut [f64::NAN, f64::NAN]).unwrap_err();
        assert_eq!(error, DsterfError::NonConvergence { remaining: 2 });
        assert!(d.iter().any(|x| x.is_nan()));
        assert!(matches!(
            dsterf(&mut [0., 0.], &mut []),
            Err(DsterfError::InvalidOffDiagonalLength { .. })
        ));
    }
}
