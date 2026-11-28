//! Eigenvalues-only, lower-triangle DSYEVD path for default IPC.
//! LAPACK/OpenBLAS 446c436e10450c348169808f1a6b3fae0925c9f7.
//! LAPACK BSD notice is retained in LICENSE and lapack_tridiagonal.rs.
//! Numerical operation ordering remains under verification against the frozen build.

use super::lapack_blas as blas;
use super::lapack_tridiagonal::{dlamch, dlapy2};

/// DLARFG with unit stride: the first value is ALPHA and the rest are X.
/// The caller owns its detached column/work vector; no matrix or vector clone.
pub(super) fn larfg(vector: &mut [f64]) -> f64 {
    // LAPACK❗✔️:       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *       .. Scalar Arguments ..
    // LAPACK❗✔️: *       INTEGER            INCX, N
    // LAPACK❗✔️: *       DOUBLE PRECISION   ALPHA, TAU
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *       .. Array Arguments ..
    // LAPACK❗✔️: *       DOUBLE PRECISION   X( * )
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Purpose:
    // LAPACK❗✔️: *  =============
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> DLARFG generates a real elementary reflector H of order n, such
    // LAPACK❗✔️: *> that
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>       H * ( alpha ) = ( beta ),   H**T * H = I.
    // LAPACK❗✔️: *>           (   x   )   (   0  )
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> where alpha and beta are scalars, and x is an (n-1)-element real
    // LAPACK❗✔️: *> vector. H is represented in the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
    // LAPACK❗✔️: *>                     ( v )
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> where tau is a real scalar and v is a real (n-1)-element
    // LAPACK❗✔️: *> vector.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> If the elements of x are all zero, then tau = 0 and H is taken to be
    // LAPACK❗✔️: *> the unit matrix.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> Otherwise  1 <= tau <= 2.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Arguments:
    // LAPACK❗✔️: *  ==========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \param[in] N
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          N is INTEGER
    // LAPACK❗✔️: *>          The order of the elementary reflector.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in,out] ALPHA
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          ALPHA is DOUBLE PRECISION
    // LAPACK❗✔️: *>          On entry, the value alpha.
    // LAPACK❗✔️: *>          On exit, it is overwritten with the value beta.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in,out] X
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          X is DOUBLE PRECISION array, dimension
    // LAPACK❗✔️: *>                         (1+(N-2)*abs(INCX))
    // LAPACK❗✔️: *>          On entry, the vector x.
    // LAPACK❗✔️: *>          On exit, it is overwritten with the vector v.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] INCX
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          INCX is INTEGER
    // LAPACK❗✔️: *>          The increment between elements of X. INCX > 0.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] TAU
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          TAU is DOUBLE PRECISION
    // LAPACK❗✔️: *>          The value tau.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Authors:
    // LAPACK❗✔️: *  ========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \author Univ. of Tennessee
    // LAPACK❗✔️: *> \author Univ. of California Berkeley
    // LAPACK❗✔️: *> \author Univ. of Colorado Denver
    // LAPACK❗✔️: *> \author NAG Ltd.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \ingroup doubleOTHERauxiliary
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️:       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       INTEGER            INCX, N
    // LAPACK❗✔️:       DOUBLE PRECISION   ALPHA, TAU
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   X( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE, ZERO
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            J, KNT
    // LAPACK❗✔️:       DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
    // LAPACK❗✔️:       EXTERNAL           DLAMCH, DLAPY2, DNRM2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, SIGN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DSCAL
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.1 ) THEN
    // LAPACK❗✔️:          TAU = ZERO
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       XNORM = DNRM2( N-1, X, INCX )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( XNORM.EQ.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        H  =  I
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          TAU = ZERO
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        general case
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
    // LAPACK❗✔️:          SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
    // LAPACK❗✔️:          KNT = 0
    // LAPACK❗✔️:          IF( ABS( BETA ).LT.SAFMIN ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           XNORM, BETA may be inaccurate; scale X and recompute them
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             RSAFMN = ONE / SAFMIN
    // LAPACK❗✔️:    10       CONTINUE
    // LAPACK❗✔️:             KNT = KNT + 1
    // LAPACK❗✔️:             CALL DSCAL( N-1, RSAFMN, X, INCX )
    // LAPACK❗✔️:             BETA = BETA*RSAFMN
    // LAPACK❗✔️:             ALPHA = ALPHA*RSAFMN
    // LAPACK❗✔️:             IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) )
    // LAPACK❗✔️:      $         GO TO 10
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           New BETA is at most 1, at least SAFMIN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             XNORM = DNRM2( N-1, X, INCX )
    // LAPACK❗✔️:             BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          TAU = ( BETA-ALPHA ) / BETA
    // LAPACK❗✔️:          CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        If ALPHA is subnormal, it may lose relative accuracy
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 20 J = 1, KNT
    // LAPACK❗✔️:             BETA = BETA*SAFMIN
    // LAPACK❗✔️:  20      CONTINUE
    // LAPACK❗✔️:          ALPHA = BETA
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLARFG
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if vector.len() <= 1 {
        return 0.0;
    }
    let n = vector.len() - 1;
    let mut alpha = vector[0];
    let mut xnorm = blas::nrm2(n, &vector[1..], 1);
    if xnorm == 0.0 {
        return 0.0;
    }
    let mut beta = -dlapy2(alpha, xnorm).copysign(alpha);
    let safmin = dlamch(b'S') / dlamch(b'E');
    let mut count = 0;
    if beta.abs() < safmin {
        let reciprocal = 1.0 / safmin;
        loop {
            count += 1;
            blas::scal(n, reciprocal, &mut vector[1..], 1);
            beta *= reciprocal;
            alpha *= reciprocal;
            if !(beta.abs() < safmin && count < 20) {
                break;
            }
        }
        xnorm = blas::nrm2(n, &vector[1..], 1);
        beta = -dlapy2(alpha, xnorm).copysign(alpha);
    }
    let tau = (beta - alpha) / beta;
    blas::scal(n, 1.0 / (alpha - beta), &mut vector[1..], 1);
    for _ in 0..count {
        beta *= safmin;
    }
    vector[0] = beta;
    tau
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn larfg_source_reflector_and_subnormal_scaling() {
        let mut vector = [3.0, 4.0];
        let tau = larfg(&mut vector);
        assert_eq!(tau.to_bits(), 1.6_f64.to_bits());
        assert_eq!(
            vector.map(f64::to_bits),
            [(-5.0_f64).to_bits(), 0.5_f64.to_bits()]
        );
        let mut identity = [-0.0, 0.0];
        assert_eq!(larfg(&mut identity).to_bits(), 0.0_f64.to_bits());
        assert_eq!(identity[0].to_bits(), (-0.0_f64).to_bits());
        let mut signed = [-0.0, 3.0];
        assert_eq!(larfg(&mut signed).to_bits(), 1.0_f64.to_bits());
        assert_eq!(
            signed.map(f64::to_bits),
            [3.0_f64.to_bits(), (-1.0_f64).to_bits()]
        );
        let tiny = f64::from_bits(23_u64 << 52); // 2^-1000.
        let mut scaled = [0.0, tiny];
        assert_eq!(larfg(&mut scaled).to_bits(), 1.0_f64.to_bits());
        assert_eq!(
            scaled.map(f64::to_bits),
            [(-tiny).to_bits(), 1.0_f64.to_bits()]
        );
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SymmetricError {
    Tridiagonal(super::lapack_tridiagonal::DsterfError),
    Blas(super::lapack_blas_level3::Blas3Error),
}
impl From<super::lapack_tridiagonal::DsterfError> for SymmetricError {
    fn from(error: super::lapack_tridiagonal::DsterfError) -> Self {
        Self::Tridiagonal(error)
    }
}
impl From<super::lapack_blas_level3::Blas3Error> for SymmetricError {
    fn from(error: super::lapack_blas_level3::Blas3Error) -> Self {
        Self::Blas(error)
    }
}

/// Fixed DSYEVD JOBZ=N, UPLO=L with the queried preferred workspace.
/// Input is a checked square column-major buffer, consumed by the reduction.
pub(super) fn eigenvalues(n: usize, mut a: Vec<f64>) -> Result<Vec<f64>, SymmetricError> {
    // LAPACK❗✔️:       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
    // LAPACK❗✔️: *                          LIWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *       .. Scalar Arguments ..
    // LAPACK❗✔️: *       CHARACTER          JOBZ, UPLO
    // LAPACK❗✔️: *       INTEGER            INFO, LDA, LIWORK, LWORK, N
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *       .. Array Arguments ..
    // LAPACK❗✔️: *       INTEGER            IWORK( * )
    // LAPACK❗✔️: *       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Purpose:
    // LAPACK❗✔️: *  =============
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
    // LAPACK❗✔️: *> real symmetric matrix A. If eigenvectors are desired, it uses a
    // LAPACK❗✔️: *> divide and conquer algorithm.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> Because of large use of BLAS of level 3, DSYEVD needs N**2 more
    // LAPACK❗✔️: *> workspace than DSYEVX.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Arguments:
    // LAPACK❗✔️: *  ==========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \param[in] JOBZ
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          JOBZ is CHARACTER*1
    // LAPACK❗✔️: *>          = 'N':  Compute eigenvalues only;
    // LAPACK❗✔️: *>          = 'V':  Compute eigenvalues and eigenvectors.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] UPLO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          UPLO is CHARACTER*1
    // LAPACK❗✔️: *>          = 'U':  Upper triangle of A is stored;
    // LAPACK❗✔️: *>          = 'L':  Lower triangle of A is stored.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] N
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          N is INTEGER
    // LAPACK❗✔️: *>          The order of the matrix A.  N >= 0.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in,out] A
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          A is DOUBLE PRECISION array, dimension (LDA, N)
    // LAPACK❗✔️: *>          On entry, the symmetric matrix A.  If UPLO = 'U', the
    // LAPACK❗✔️: *>          leading N-by-N upper triangular part of A contains the
    // LAPACK❗✔️: *>          upper triangular part of the matrix A.  If UPLO = 'L',
    // LAPACK❗✔️: *>          the leading N-by-N lower triangular part of A contains
    // LAPACK❗✔️: *>          the lower triangular part of the matrix A.
    // LAPACK❗✔️: *>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
    // LAPACK❗✔️: *>          orthonormal eigenvectors of the matrix A.
    // LAPACK❗✔️: *>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
    // LAPACK❗✔️: *>          or the upper triangle (if UPLO='U') of A, including the
    // LAPACK❗✔️: *>          diagonal, is destroyed.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LDA
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LDA is INTEGER
    // LAPACK❗✔️: *>          The leading dimension of the array A.  LDA >= max(1,N).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] W
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          W is DOUBLE PRECISION array, dimension (N)
    // LAPACK❗✔️: *>          If INFO = 0, the eigenvalues in ascending order.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] WORK
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    // LAPACK❗✔️: *>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LWORK
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LWORK is INTEGER
    // LAPACK❗✔️: *>          The dimension of the array WORK.
    // LAPACK❗✔️: *>          If N <= 1,               LWORK must be at least 1.
    // LAPACK❗✔️: *>          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
    // LAPACK❗✔️: *>          If JOBZ = 'V' and N > 1, LWORK must be at least
    // LAPACK❗✔️: *>                                                1 + 6*N + 2*N**2.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>          If LWORK = -1, then a workspace query is assumed; the routine
    // LAPACK❗✔️: *>          only calculates the optimal sizes of the WORK and IWORK
    // LAPACK❗✔️: *>          arrays, returns these values as the first entries of the WORK
    // LAPACK❗✔️: *>          and IWORK arrays, and no error message related to LWORK or
    // LAPACK❗✔️: *>          LIWORK is issued by XERBLA.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] IWORK
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
    // LAPACK❗✔️: *>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LIWORK
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LIWORK is INTEGER
    // LAPACK❗✔️: *>          The dimension of the array IWORK.
    // LAPACK❗✔️: *>          If N <= 1,                LIWORK must be at least 1.
    // LAPACK❗✔️: *>          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
    // LAPACK❗✔️: *>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>          If LIWORK = -1, then a workspace query is assumed; the
    // LAPACK❗✔️: *>          routine only calculates the optimal sizes of the WORK and
    // LAPACK❗✔️: *>          IWORK arrays, returns these values as the first entries of
    // LAPACK❗✔️: *>          the WORK and IWORK arrays, and no error message related to
    // LAPACK❗✔️: *>          LWORK or LIWORK is issued by XERBLA.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] INFO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          INFO is INTEGER
    // LAPACK❗✔️: *>          = 0:  successful exit
    // LAPACK❗✔️: *>          < 0:  if INFO = -i, the i-th argument had an illegal value
    // LAPACK❗✔️: *>          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
    // LAPACK❗✔️: *>                to converge; i off-diagonal elements of an intermediate
    // LAPACK❗✔️: *>                tridiagonal form did not converge to zero;
    // LAPACK❗✔️: *>                if INFO = i and JOBZ = 'V', then the algorithm failed
    // LAPACK❗✔️: *>                to compute an eigenvalue while working on the submatrix
    // LAPACK❗✔️: *>                lying in rows and columns INFO/(N+1) through
    // LAPACK❗✔️: *>                mod(INFO,N+1).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Authors:
    // LAPACK❗✔️: *  ========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \author Univ. of Tennessee
    // LAPACK❗✔️: *> \author Univ. of California Berkeley
    // LAPACK❗✔️: *> \author Univ. of Colorado Denver
    // LAPACK❗✔️: *> \author NAG Ltd.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \ingroup heevd
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Contributors:
    // LAPACK❗✔️: *  ==================
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> Jeff Rutter, Computer Science Division, University of California
    // LAPACK❗✔️: *> at Berkeley, USA \n
    // LAPACK❗✔️: *>  Modified by Francoise Tisseur, University of Tennessee \n
    // LAPACK❗✔️: *>  Modified description of INFO. Sven, 16 Feb 05. \n
    // LAPACK❗✔️:
    // LAPACK❗✔️:
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️:       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,
    // LAPACK❗✔️:      $                   LIWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK driver routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          JOBZ, UPLO
    // LAPACK❗✔️:       INTEGER            INFO, LDA, LIWORK, LWORK, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       INTEGER            IWORK( * )
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       LOGICAL            LOWER, LQUERY, WANTZ
    // LAPACK❗✔️:       INTEGER            IINFO, INDE, INDTAU, INDWK2, INDWRK, ISCALE,
    // LAPACK❗✔️:      $                   LIOPT, LIWMIN, LLWORK, LLWRK2, LOPT, LWMIN
    // LAPACK❗✔️:       DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
    // LAPACK❗✔️:      $                   SMLNUM
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DLANSY
    // LAPACK❗✔️:       EXTERNAL           LSAME, DLAMCH, DLANSY, ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLACPY, DLASCL, DORMTR, DSCAL, DSTEDC, DSTERF,
    // LAPACK❗✔️:      $                   DSYTRD, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WANTZ = LSAME( JOBZ, 'V' )
    // LAPACK❗✔️:       LOWER = LSAME( UPLO, 'L' )
    // LAPACK❗✔️:       LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️:          IF( N.LE.1 ) THEN
    // LAPACK❗✔️:             LIWMIN = 1
    // LAPACK❗✔️:             LWMIN = 1
    // LAPACK❗✔️:             LOPT = LWMIN
    // LAPACK❗✔️:             LIOPT = LIWMIN
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             IF( WANTZ ) THEN
    // LAPACK❗✔️:                LIWMIN = 3 + 5*N
    // LAPACK❗✔️:                LWMIN = 1 + 6*N + 2*N**2
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                LIWMIN = 1
    // LAPACK❗✔️:                LWMIN = 2*N + 1
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             LOPT = MAX( LWMIN, 2*N +
    // LAPACK❗✔️:      $                  N*ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
    // LAPACK❗✔️:             LIOPT = LIWMIN
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          WORK( 1 ) = LOPT
    // LAPACK❗✔️:          IWORK( 1 ) = LIOPT
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:             INFO = -8
    // LAPACK❗✔️:          ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:             INFO = -10
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DSYEVD', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       ELSE IF( LQUERY ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.1 ) THEN
    // LAPACK❗✔️:          W( 1 ) = A( 1, 1 )
    // LAPACK❗✔️:          IF( WANTZ )
    // LAPACK❗✔️:      $      A( 1, 1 ) = ONE
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Get machine constants.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'Safe minimum' )
    // LAPACK❗✔️:       EPS = DLAMCH( 'Precision' )
    // LAPACK❗✔️:       SMLNUM = SAFMIN / EPS
    // LAPACK❗✔️:       BIGNUM = ONE / SMLNUM
    // LAPACK❗✔️:       RMIN = SQRT( SMLNUM )
    // LAPACK❗✔️:       RMAX = SQRT( BIGNUM )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Scale matrix to allowable range, if necessary.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
    // LAPACK❗✔️:       ISCALE = 0
    // LAPACK❗✔️:       IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
    // LAPACK❗✔️:          ISCALE = 1
    // LAPACK❗✔️:          SIGMA = RMIN / ANRM
    // LAPACK❗✔️:       ELSE IF( ANRM.GT.RMAX ) THEN
    // LAPACK❗✔️:          ISCALE = 1
    // LAPACK❗✔️:          SIGMA = RMAX / ANRM
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( ISCALE.EQ.1 )
    // LAPACK❗✔️:      $   CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INDE = 1
    // LAPACK❗✔️:       INDTAU = INDE + N
    // LAPACK❗✔️:       INDWRK = INDTAU + N
    // LAPACK❗✔️:       LLWORK = LWORK - INDWRK + 1
    // LAPACK❗✔️:       INDWK2 = INDWRK + N*N
    // LAPACK❗✔️:       LLWRK2 = LWORK - INDWK2 + 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
    // LAPACK❗✔️:      $             WORK( INDWRK ), LLWORK, IINFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     For eigenvalues only, call DSTERF.  For eigenvectors, first call
    // LAPACK❗✔️: *     DSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
    // LAPACK❗✔️: *     tridiagonal matrix, then call DORMTR to multiply it by the
    // LAPACK❗✔️: *     Householder transformations stored in A.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( .NOT.WANTZ ) THEN
    // LAPACK❗✔️:          CALL DSTERF( N, W, WORK( INDE ), INFO )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          CALL DSTEDC( 'I', N, W, WORK( INDE ), WORK( INDWRK ), N,
    // LAPACK❗✔️:      $                WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO )
    // LAPACK❗✔️:          CALL DORMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ),
    // LAPACK❗✔️:      $                WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
    // LAPACK❗✔️:          CALL DLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     If matrix was scaled, then rescale eigenvalues appropriately.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( ISCALE.EQ.1 )
    // LAPACK❗✔️:      $   CALL DSCAL( N, ONE / SIGMA, W, 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = LOPT
    // LAPACK❗✔️:       IWORK( 1 ) = LIOPT
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DSYEVD
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    // The mode is fixed by NumPy eigvalsh_lo. Neither eigenvector machinery
    // nor the unused upper-triangle path is part of this private interface.
    // The full matrix is the input buffer; work is O(n*32), not another n² copy.
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return Ok(vec![a[0]]);
    }
    let safe_min = dlamch(b'S');
    let epsilon = dlamch(b'P');
    let small = safe_min / epsilon;
    let large = 1.0 / small;
    let minimum = small.sqrt();
    let maximum = large.sqrt();
    let mut norm = 0.0_f64;
    for j in 0..n {
        for i in j..n {
            let value = a[i + j * n].abs();
            if norm < value || value.is_nan() {
                norm = value;
            }
        }
    }
    let scale = if norm > 0.0 && norm < minimum {
        Some(minimum / norm)
    } else if norm > maximum {
        Some(maximum / norm)
    } else {
        None
    };
    if let Some(scale) = scale {
        // DLASCL(L) applies the same staged factor to each lower column.
        // Columns are independent, so the per-column call preserves every
        // element's operation order without copying the matrix triangle.
        for j in 0..n {
            super::lapack_tridiagonal::dlascl_vector(1.0, scale, &mut a[j + j * n..n + j * n])?;
        }
    }
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    let mut tau = vec![0.0; n];
    sytrd_lower(n, &mut a, &mut d, &mut e, &mut tau)?;
    let outcome = super::lapack_tridiagonal::dsterf(&mut d, &mut e);
    // DSYEVD rescales all W entries even when DSTERF returns positive INFO.
    if let Some(scale) = scale {
        blas::scal(n, 1.0 / scale, &mut d, 1);
    }
    outcome?;
    Ok(d)
}

fn sytrd_lower(
    n: usize,
    a: &mut [f64],
    d: &mut [f64],
    e: &mut [f64],
    tau: &mut [f64],
) -> Result<(), SymmetricError> {
    // LAPACK❗✔️:       SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *       .. Scalar Arguments ..
    // LAPACK❗✔️: *       CHARACTER          UPLO
    // LAPACK❗✔️: *       INTEGER            INFO, LDA, LWORK, N
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *       .. Array Arguments ..
    // LAPACK❗✔️: *       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
    // LAPACK❗✔️: *      $                   WORK( * )
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Purpose:
    // LAPACK❗✔️: *  =============
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> DSYTRD reduces a real symmetric matrix A to real symmetric
    // LAPACK❗✔️: *> tridiagonal form T by an orthogonal similarity transformation:
    // LAPACK❗✔️: *> Q**T * A * Q = T.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Arguments:
    // LAPACK❗✔️: *  ==========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \param[in] UPLO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          UPLO is CHARACTER*1
    // LAPACK❗✔️: *>          = 'U':  Upper triangle of A is stored;
    // LAPACK❗✔️: *>          = 'L':  Lower triangle of A is stored.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] N
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          N is INTEGER
    // LAPACK❗✔️: *>          The order of the matrix A.  N >= 0.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in,out] A
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          A is DOUBLE PRECISION array, dimension (LDA,N)
    // LAPACK❗✔️: *>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
    // LAPACK❗✔️: *>          N-by-N upper triangular part of A contains the upper
    // LAPACK❗✔️: *>          triangular part of the matrix A, and the strictly lower
    // LAPACK❗✔️: *>          triangular part of A is not referenced.  If UPLO = 'L', the
    // LAPACK❗✔️: *>          leading N-by-N lower triangular part of A contains the lower
    // LAPACK❗✔️: *>          triangular part of the matrix A, and the strictly upper
    // LAPACK❗✔️: *>          triangular part of A is not referenced.
    // LAPACK❗✔️: *>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
    // LAPACK❗✔️: *>          of A are overwritten by the corresponding elements of the
    // LAPACK❗✔️: *>          tridiagonal matrix T, and the elements above the first
    // LAPACK❗✔️: *>          superdiagonal, with the array TAU, represent the orthogonal
    // LAPACK❗✔️: *>          matrix Q as a product of elementary reflectors; if UPLO
    // LAPACK❗✔️: *>          = 'L', the diagonal and first subdiagonal of A are over-
    // LAPACK❗✔️: *>          written by the corresponding elements of the tridiagonal
    // LAPACK❗✔️: *>          matrix T, and the elements below the first subdiagonal, with
    // LAPACK❗✔️: *>          the array TAU, represent the orthogonal matrix Q as a product
    // LAPACK❗✔️: *>          of elementary reflectors. See Further Details.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LDA
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LDA is INTEGER
    // LAPACK❗✔️: *>          The leading dimension of the array A.  LDA >= max(1,N).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] D
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          D is DOUBLE PRECISION array, dimension (N)
    // LAPACK❗✔️: *>          The diagonal elements of the tridiagonal matrix T:
    // LAPACK❗✔️: *>          D(i) = A(i,i).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] E
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          E is DOUBLE PRECISION array, dimension (N-1)
    // LAPACK❗✔️: *>          The off-diagonal elements of the tridiagonal matrix T:
    // LAPACK❗✔️: *>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] TAU
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          TAU is DOUBLE PRECISION array, dimension (N-1)
    // LAPACK❗✔️: *>          The scalar factors of the elementary reflectors (see Further
    // LAPACK❗✔️: *>          Details).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] WORK
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    // LAPACK❗✔️: *>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LWORK
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LWORK is INTEGER
    // LAPACK❗✔️: *>          The dimension of the array WORK.  LWORK >= 1.
    // LAPACK❗✔️: *>          For optimum performance LWORK >= N*NB, where NB is the
    // LAPACK❗✔️: *>          optimal blocksize.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>          If LWORK = -1, then a workspace query is assumed; the routine
    // LAPACK❗✔️: *>          only calculates the optimal size of the WORK array, returns
    // LAPACK❗✔️: *>          this value as the first entry of the WORK array, and no error
    // LAPACK❗✔️: *>          message related to LWORK is issued by XERBLA.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] INFO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          INFO is INTEGER
    // LAPACK❗✔️: *>          = 0:  successful exit
    // LAPACK❗✔️: *>          < 0:  if INFO = -i, the i-th argument had an illegal value
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Authors:
    // LAPACK❗✔️: *  ========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \author Univ. of Tennessee
    // LAPACK❗✔️: *> \author Univ. of California Berkeley
    // LAPACK❗✔️: *> \author Univ. of Colorado Denver
    // LAPACK❗✔️: *> \author NAG Ltd.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \ingroup hetrd
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Further Details:
    // LAPACK❗✔️: *  =====================
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  If UPLO = 'U', the matrix Q is represented as a product of elementary
    // LAPACK❗✔️: *>  reflectors
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     Q = H(n-1) . . . H(2) H(1).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  Each H(i) has the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where tau is a real scalar, and v is a real vector with
    // LAPACK❗✔️: *>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    // LAPACK❗✔️: *>  A(1:i-1,i+1), and tau in TAU(i).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  If UPLO = 'L', the matrix Q is represented as a product of elementary
    // LAPACK❗✔️: *>  reflectors
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     Q = H(1) H(2) . . . H(n-1).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  Each H(i) has the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where tau is a real scalar, and v is a real vector with
    // LAPACK❗✔️: *>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    // LAPACK❗✔️: *>  and tau in TAU(i).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  The contents of A on exit are illustrated by the following examples
    // LAPACK❗✔️: *>  with n = 5:
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  if UPLO = 'U':                       if UPLO = 'L':
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>    (  d   e   v2  v3  v4 )              (  d                  )
    // LAPACK❗✔️: *>    (      d   e   v3  v4 )              (  e   d              )
    // LAPACK❗✔️: *>    (          d   e   v4 )              (  v1  e   d          )
    // LAPACK❗✔️: *>    (              d   e  )              (  v1  v2  e   d      )
    // LAPACK❗✔️: *>    (                  d  )              (  v1  v2  v3  e   d  )
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where d and e denote diagonal and off-diagonal elements of T, and vi
    // LAPACK❗✔️: *>  denotes an element of the vector defining H(i).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️:       SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          UPLO
    // LAPACK❗✔️:       INTEGER            INFO, LDA, LWORK, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
    // LAPACK❗✔️:      $                   WORK( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            LQUERY, UPPER
    // LAPACK❗✔️:       INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
    // LAPACK❗✔️:      $                   NBMIN, NX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       INTEGER            ILAENV
    // LAPACK❗✔️:       EXTERNAL           LSAME, ILAENV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       UPPER = LSAME( UPLO, 'U' )
    // LAPACK❗✔️:       LQUERY = ( LWORK.EQ.-1 )
    // LAPACK❗✔️:       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:          INFO = -9
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Determine the block size.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
    // LAPACK❗✔️:          LWKOPT = MAX( 1, N*NB )
    // LAPACK❗✔️:          WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DSYTRD', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       ELSE IF( LQUERY ) THEN
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 ) THEN
    // LAPACK❗✔️:          WORK( 1 ) = 1
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NX = N
    // LAPACK❗✔️:       IWS = 1
    // LAPACK❗✔️:       IF( NB.GT.1 .AND. NB.LT.N ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Determine when to cross over from blocked to unblocked code
    // LAPACK❗✔️: *        (last block is always handled by unblocked code).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
    // LAPACK❗✔️:          IF( NX.LT.N ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Determine if workspace is large enough for blocked code.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             LDWORK = N
    // LAPACK❗✔️:             IWS = LDWORK*NB
    // LAPACK❗✔️:             IF( LWORK.LT.IWS ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Not enough workspace to use optimal NB:  determine the
    // LAPACK❗✔️: *              minimum value of NB, and reduce NB or force use of
    // LAPACK❗✔️: *              unblocked code by setting NX = N.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                NB = MAX( LWORK / LDWORK, 1 )
    // LAPACK❗✔️:                NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
    // LAPACK❗✔️:                IF( NB.LT.NBMIN )
    // LAPACK❗✔️:      $            NX = N
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             NX = N
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          NB = 1
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( UPPER ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Reduce the upper triangle of A.
    // LAPACK❗✔️: *        Columns 1:kk are handled by the unblocked method.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          KK = N - ( ( N-NX+NB-1 ) / NB )*NB
    // LAPACK❗✔️:          DO 20 I = N - NB + 1, KK + 1, -NB
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Reduce columns i:i+nb-1 to tridiagonal form and form the
    // LAPACK❗✔️: *           matrix W which is needed to update the unreduced part of
    // LAPACK❗✔️: *           the matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
    // LAPACK❗✔️:      $                   LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Update the unreduced submatrix A(1:i-1,1:i-1), using an
    // LAPACK❗✔️: *           update of the form:  A := A - V*W**T - W*V**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
    // LAPACK❗✔️:      $                   LDA, WORK, LDWORK, ONE, A, LDA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Copy superdiagonal elements back into A, and diagonal
    // LAPACK❗✔️: *           elements into D
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 10 J = I, I + NB - 1
    // LAPACK❗✔️:                A( J-1, J ) = E( J-1 )
    // LAPACK❗✔️:                D( J ) = A( J, J )
    // LAPACK❗✔️:    10       CONTINUE
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Use unblocked code to reduce the last or only block
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Reduce the lower triangle of A
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 40 I = 1, N - NX, NB
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Reduce columns i:i+nb-1 to tridiagonal form and form the
    // LAPACK❗✔️: *           matrix W which is needed to update the unreduced part of
    // LAPACK❗✔️: *           the matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
    // LAPACK❗✔️:      $                   TAU( I ), WORK, LDWORK )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
    // LAPACK❗✔️: *           an update of the form:  A := A - V*W**T - W*V**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
    // LAPACK❗✔️:      $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
    // LAPACK❗✔️:      $                   A( I+NB, I+NB ), LDA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Copy subdiagonal elements back into A, and diagonal
    // LAPACK❗✔️: *           elements into D
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO 30 J = I, I + NB - 1
    // LAPACK❗✔️:                A( J+1, J ) = E( J )
    // LAPACK❗✔️:                D( J ) = A( J, J )
    // LAPACK❗✔️:    30       CONTINUE
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Use unblocked code to reduce the last or only block
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
    // LAPACK❗✔️:      $                TAU( I ), IINFO )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = LWKOPT
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DSYTRD
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if n == 0 {
        return Ok(());
    }
    // ILAENV(1,DSYTRD)=32, ILAENV(3,DSYTRD)=32. DSYEVD queries
    // 2*n+n*32 entries, so the panel workspace is never artificially reduced.
    const BLOCK: usize = 32;
    let mut start = 0;
    if n > BLOCK {
        let mut work = vec![0.0; n * BLOCK];
        while start < n - BLOCK {
            let remaining = n - start;
            let sub = &mut a[start + start * n..];
            latrd_lower(
                remaining,
                BLOCK,
                sub,
                n,
                &mut e[start..],
                &mut tau[start..],
                &mut work,
                n,
            );
            let (panels, trailing) = sub.split_at_mut(BLOCK * n);
            super::lapack_blas_level3::syr2k_lower(
                remaining - BLOCK,
                BLOCK,
                -1.0,
                &panels[BLOCK..],
                n,
                &work[BLOCK..],
                n,
                1.0,
                &mut trailing[BLOCK..],
                n,
            )?;
            for j in start..start + BLOCK {
                a[j + 1 + j * n] = e[j];
                d[j] = a[j + j * n];
            }
            start += BLOCK;
        }
    }
    sytd2_lower(
        n - start,
        &mut a[start + start * n..],
        n,
        &mut d[start..],
        &mut e[start..],
        &mut tau[start..],
    );
    Ok(())
}

fn sytd2_lower(n: usize, a: &mut [f64], lda: usize, d: &mut [f64], e: &mut [f64], tau: &mut [f64]) {
    // LAPACK❗✔️:       SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *       .. Scalar Arguments ..
    // LAPACK❗✔️: *       CHARACTER          UPLO
    // LAPACK❗✔️: *       INTEGER            INFO, LDA, N
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *       .. Array Arguments ..
    // LAPACK❗✔️: *       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Purpose:
    // LAPACK❗✔️: *  =============
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
    // LAPACK❗✔️: *> form T by an orthogonal similarity transformation: Q**T * A * Q = T.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Arguments:
    // LAPACK❗✔️: *  ==========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \param[in] UPLO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          UPLO is CHARACTER*1
    // LAPACK❗✔️: *>          Specifies whether the upper or lower triangular part of the
    // LAPACK❗✔️: *>          symmetric matrix A is stored:
    // LAPACK❗✔️: *>          = 'U':  Upper triangular
    // LAPACK❗✔️: *>          = 'L':  Lower triangular
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] N
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          N is INTEGER
    // LAPACK❗✔️: *>          The order of the matrix A.  N >= 0.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in,out] A
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          A is DOUBLE PRECISION array, dimension (LDA,N)
    // LAPACK❗✔️: *>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
    // LAPACK❗✔️: *>          n-by-n upper triangular part of A contains the upper
    // LAPACK❗✔️: *>          triangular part of the matrix A, and the strictly lower
    // LAPACK❗✔️: *>          triangular part of A is not referenced.  If UPLO = 'L', the
    // LAPACK❗✔️: *>          leading n-by-n lower triangular part of A contains the lower
    // LAPACK❗✔️: *>          triangular part of the matrix A, and the strictly upper
    // LAPACK❗✔️: *>          triangular part of A is not referenced.
    // LAPACK❗✔️: *>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
    // LAPACK❗✔️: *>          of A are overwritten by the corresponding elements of the
    // LAPACK❗✔️: *>          tridiagonal matrix T, and the elements above the first
    // LAPACK❗✔️: *>          superdiagonal, with the array TAU, represent the orthogonal
    // LAPACK❗✔️: *>          matrix Q as a product of elementary reflectors; if UPLO
    // LAPACK❗✔️: *>          = 'L', the diagonal and first subdiagonal of A are over-
    // LAPACK❗✔️: *>          written by the corresponding elements of the tridiagonal
    // LAPACK❗✔️: *>          matrix T, and the elements below the first subdiagonal, with
    // LAPACK❗✔️: *>          the array TAU, represent the orthogonal matrix Q as a product
    // LAPACK❗✔️: *>          of elementary reflectors. See Further Details.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LDA
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LDA is INTEGER
    // LAPACK❗✔️: *>          The leading dimension of the array A.  LDA >= max(1,N).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] D
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          D is DOUBLE PRECISION array, dimension (N)
    // LAPACK❗✔️: *>          The diagonal elements of the tridiagonal matrix T:
    // LAPACK❗✔️: *>          D(i) = A(i,i).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] E
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          E is DOUBLE PRECISION array, dimension (N-1)
    // LAPACK❗✔️: *>          The off-diagonal elements of the tridiagonal matrix T:
    // LAPACK❗✔️: *>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] TAU
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          TAU is DOUBLE PRECISION array, dimension (N-1)
    // LAPACK❗✔️: *>          The scalar factors of the elementary reflectors (see Further
    // LAPACK❗✔️: *>          Details).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] INFO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          INFO is INTEGER
    // LAPACK❗✔️: *>          = 0:  successful exit
    // LAPACK❗✔️: *>          < 0:  if INFO = -i, the i-th argument had an illegal value.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Authors:
    // LAPACK❗✔️: *  ========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \author Univ. of Tennessee
    // LAPACK❗✔️: *> \author Univ. of California Berkeley
    // LAPACK❗✔️: *> \author Univ. of Colorado Denver
    // LAPACK❗✔️: *> \author NAG Ltd.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \ingroup doubleSYcomputational
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Further Details:
    // LAPACK❗✔️: *  =====================
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  If UPLO = 'U', the matrix Q is represented as a product of elementary
    // LAPACK❗✔️: *>  reflectors
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     Q = H(n-1) . . . H(2) H(1).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  Each H(i) has the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where tau is a real scalar, and v is a real vector with
    // LAPACK❗✔️: *>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    // LAPACK❗✔️: *>  A(1:i-1,i+1), and tau in TAU(i).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  If UPLO = 'L', the matrix Q is represented as a product of elementary
    // LAPACK❗✔️: *>  reflectors
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     Q = H(1) H(2) . . . H(n-1).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  Each H(i) has the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where tau is a real scalar, and v is a real vector with
    // LAPACK❗✔️: *>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    // LAPACK❗✔️: *>  and tau in TAU(i).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  The contents of A on exit are illustrated by the following examples
    // LAPACK❗✔️: *>  with n = 5:
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  if UPLO = 'U':                       if UPLO = 'L':
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>    (  d   e   v2  v3  v4 )              (  d                  )
    // LAPACK❗✔️: *>    (      d   e   v3  v4 )              (  e   d              )
    // LAPACK❗✔️: *>    (          d   e   v4 )              (  v1  e   d          )
    // LAPACK❗✔️: *>    (              d   e  )              (  v1  v2  e   d      )
    // LAPACK❗✔️: *>    (                  d  )              (  v1  v2  v3  e   d  )
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where d and e denote diagonal and off-diagonal elements of T, and vi
    // LAPACK❗✔️: *>  denotes an element of the vector defining H(i).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️:       SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          UPLO
    // LAPACK❗✔️:       INTEGER            INFO, LDA, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ONE, ZERO, HALF
    // LAPACK❗✔️:       PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
    // LAPACK❗✔️:      $                   HALF = 1.0D0 / 2.0D0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            UPPER
    // LAPACK❗✔️:       INTEGER            I
    // LAPACK❗✔️:       DOUBLE PRECISION   ALPHA, TAUI
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       DOUBLE PRECISION   DDOT
    // LAPACK❗✔️:       EXTERNAL           LSAME, DDOT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input parameters
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       UPPER = LSAME( UPLO, 'U' )
    // LAPACK❗✔️:       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DSYTD2', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( UPPER ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Reduce the upper triangle of A
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 10 I = N - 1, 1, -1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Generate elementary reflector H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *           to annihilate A(1:i-1,i+1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
    // LAPACK❗✔️:             E( I ) = A( I, I+1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( TAUI.NE.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply H(i) from both sides to A(1:i,1:i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                A( I, I+1 ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Compute  x := tau * A * v  storing x in TAU(1:i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
    // LAPACK❗✔️:      $                     TAU, 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Compute  w := x - 1/2 * tau * (x**T * v) * v
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
    // LAPACK❗✔️:                CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply the transformation as a rank-2 update:
    // LAPACK❗✔️: *                 A := A - v * w**T - w * v**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
    // LAPACK❗✔️:      $                     LDA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                A( I, I+1 ) = E( I )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             D( I+1 ) = A( I+1, I+1 )
    // LAPACK❗✔️:             TAU( I ) = TAUI
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️:          D( 1 ) = A( 1, 1 )
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Reduce the lower triangle of A
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 20 I = 1, N - 1
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Generate elementary reflector H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *           to annihilate A(i+2:n,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
    // LAPACK❗✔️:      $                   TAUI )
    // LAPACK❗✔️:             E( I ) = A( I+1, I )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( TAUI.NE.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply H(i) from both sides to A(i+1:n,i+1:n)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                A( I+1, I ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Compute  x := tau * A * v  storing y in TAU(i:n-1)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
    // LAPACK❗✔️:      $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Compute  w := x - 1/2 * tau * (x**T * v) * v
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
    // LAPACK❗✔️:      $                 1 )
    // LAPACK❗✔️:                CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Apply the transformation as a rank-2 update:
    // LAPACK❗✔️: *                 A := A - v * w**T - w * v**T
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
    // LAPACK❗✔️:      $                     A( I+1, I+1 ), LDA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                A( I+1, I ) = E( I )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             D( I ) = A( I, I )
    // LAPACK❗✔️:             TAU( I ) = TAUI
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️:          D( N ) = A( N, N )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DSYTD2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    for i in 0..n - 1 {
        let m = n - i - 1;
        let reflector = larfg(&mut a[i + 1 + i * lda..n + i * lda]);
        e[i] = a[i + 1 + i * lda];
        if reflector != 0.0 {
            a[i + 1 + i * lda] = 1.0;
            let trailing = i + 1 + (i + 1) * lda;
            let (previous, lower) = a.split_at_mut(trailing);
            let v = &previous[i + 1 + i * lda..n + i * lda];
            let w = &mut tau[i..i + m];
            blas::symv_lower(m, reflector, lower, lda, v, 1, 0.0, w, 1);
            let alpha = (-0.5 * reflector) * blas::dot(m, w, 1, v, 1);
            blas::axpy(m, alpha, v, 1, w, 1);
            blas::syr2_lower(m, -1.0, v, 1, w, 1, lower, lda);
            a[i + 1 + i * lda] = e[i];
        }
        d[i] = a[i + i * lda];
        tau[i] = reflector;
    }
    d[n - 1] = a[n - 1 + (n - 1) * lda];
}

fn latrd_lower(
    n: usize,
    nb: usize,
    a: &mut [f64],
    lda: usize,
    e: &mut [f64],
    tau: &mut [f64],
    w: &mut [f64],
    ldw: usize,
) {
    // LAPACK❗✔️:       SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *       .. Scalar Arguments ..
    // LAPACK❗✔️: *       CHARACTER          UPLO
    // LAPACK❗✔️: *       INTEGER            LDA, LDW, N, NB
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *       .. Array Arguments ..
    // LAPACK❗✔️: *       DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
    // LAPACK❗✔️: *       ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Purpose:
    // LAPACK❗✔️: *  =============
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> DLATRD reduces NB rows and columns of a real symmetric matrix A to
    // LAPACK❗✔️: *> symmetric tridiagonal form by an orthogonal similarity
    // LAPACK❗✔️: *> transformation Q**T * A * Q, and returns the matrices V and W which are
    // LAPACK❗✔️: *> needed to apply the transformation to the unreduced part of A.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
    // LAPACK❗✔️: *> matrix, of which the upper triangle is supplied;
    // LAPACK❗✔️: *> if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
    // LAPACK❗✔️: *> matrix, of which the lower triangle is supplied.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> This is an auxiliary routine called by DSYTRD.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Arguments:
    // LAPACK❗✔️: *  ==========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \param[in] UPLO
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          UPLO is CHARACTER*1
    // LAPACK❗✔️: *>          Specifies whether the upper or lower triangular part of the
    // LAPACK❗✔️: *>          symmetric matrix A is stored:
    // LAPACK❗✔️: *>          = 'U': Upper triangular
    // LAPACK❗✔️: *>          = 'L': Lower triangular
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] N
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          N is INTEGER
    // LAPACK❗✔️: *>          The order of the matrix A.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] NB
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          NB is INTEGER
    // LAPACK❗✔️: *>          The number of rows and columns to be reduced.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in,out] A
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          A is DOUBLE PRECISION array, dimension (LDA,N)
    // LAPACK❗✔️: *>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
    // LAPACK❗✔️: *>          n-by-n upper triangular part of A contains the upper
    // LAPACK❗✔️: *>          triangular part of the matrix A, and the strictly lower
    // LAPACK❗✔️: *>          triangular part of A is not referenced.  If UPLO = 'L', the
    // LAPACK❗✔️: *>          leading n-by-n lower triangular part of A contains the lower
    // LAPACK❗✔️: *>          triangular part of the matrix A, and the strictly upper
    // LAPACK❗✔️: *>          triangular part of A is not referenced.
    // LAPACK❗✔️: *>          On exit:
    // LAPACK❗✔️: *>          if UPLO = 'U', the last NB columns have been reduced to
    // LAPACK❗✔️: *>            tridiagonal form, with the diagonal elements overwriting
    // LAPACK❗✔️: *>            the diagonal elements of A; the elements above the diagonal
    // LAPACK❗✔️: *>            with the array TAU, represent the orthogonal matrix Q as a
    // LAPACK❗✔️: *>            product of elementary reflectors;
    // LAPACK❗✔️: *>          if UPLO = 'L', the first NB columns have been reduced to
    // LAPACK❗✔️: *>            tridiagonal form, with the diagonal elements overwriting
    // LAPACK❗✔️: *>            the diagonal elements of A; the elements below the diagonal
    // LAPACK❗✔️: *>            with the array TAU, represent the  orthogonal matrix Q as a
    // LAPACK❗✔️: *>            product of elementary reflectors.
    // LAPACK❗✔️: *>          See Further Details.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LDA
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LDA is INTEGER
    // LAPACK❗✔️: *>          The leading dimension of the array A.  LDA >= (1,N).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] E
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          E is DOUBLE PRECISION array, dimension (N-1)
    // LAPACK❗✔️: *>          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
    // LAPACK❗✔️: *>          elements of the last NB columns of the reduced matrix;
    // LAPACK❗✔️: *>          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
    // LAPACK❗✔️: *>          the first NB columns of the reduced matrix.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] TAU
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          TAU is DOUBLE PRECISION array, dimension (N-1)
    // LAPACK❗✔️: *>          The scalar factors of the elementary reflectors, stored in
    // LAPACK❗✔️: *>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
    // LAPACK❗✔️: *>          See Further Details.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[out] W
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          W is DOUBLE PRECISION array, dimension (LDW,NB)
    // LAPACK❗✔️: *>          The n-by-nb matrix W required to update the unreduced part
    // LAPACK❗✔️: *>          of A.
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \param[in] LDW
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>          LDW is INTEGER
    // LAPACK❗✔️: *>          The leading dimension of the array W. LDW >= max(1,N).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  Authors:
    // LAPACK❗✔️: *  ========
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \author Univ. of Tennessee
    // LAPACK❗✔️: *> \author Univ. of California Berkeley
    // LAPACK❗✔️: *> \author Univ. of Colorado Denver
    // LAPACK❗✔️: *> \author NAG Ltd.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \ingroup doubleOTHERauxiliary
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *> \par Further Details:
    // LAPACK❗✔️: *  =====================
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *> \verbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  If UPLO = 'U', the matrix Q is represented as a product of elementary
    // LAPACK❗✔️: *>  reflectors
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     Q = H(n) H(n-1) . . . H(n-nb+1).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  Each H(i) has the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where tau is a real scalar, and v is a real vector with
    // LAPACK❗✔️: *>  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
    // LAPACK❗✔️: *>  and tau in TAU(i-1).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  If UPLO = 'L', the matrix Q is represented as a product of elementary
    // LAPACK❗✔️: *>  reflectors
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     Q = H(1) H(2) . . . H(nb).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  Each H(i) has the form
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>     H(i) = I - tau * v * v**T
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where tau is a real scalar, and v is a real vector with
    // LAPACK❗✔️: *>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
    // LAPACK❗✔️: *>  and tau in TAU(i).
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  The elements of the vectors v together form the n-by-nb matrix V
    // LAPACK❗✔️: *>  which is needed, with W, to apply the transformation to the unreduced
    // LAPACK❗✔️: *>  part of the matrix, using a symmetric rank-2k update of the form:
    // LAPACK❗✔️: *>  A := A - V*W**T - W*V**T.
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  The contents of A on exit are illustrated by the following examples
    // LAPACK❗✔️: *>  with n = 5 and nb = 2:
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  if UPLO = 'U':                       if UPLO = 'L':
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>    (  a   a   a   v4  v5 )              (  d                  )
    // LAPACK❗✔️: *>    (      a   a   v4  v5 )              (  1   d              )
    // LAPACK❗✔️: *>    (          a   1   v5 )              (  v1  1   a          )
    // LAPACK❗✔️: *>    (              d   1  )              (  v1  v2  a   a      )
    // LAPACK❗✔️: *>    (                  d  )              (  v1  v2  a   a   a  )
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *>  where d denotes a diagonal element of the reduced matrix, a denotes
    // LAPACK❗✔️: *>  an element of the original matrix that is unchanged, and vi denotes
    // LAPACK❗✔️: *>  an element of the vector defining H(i).
    // LAPACK❗✔️: *> \endverbatim
    // LAPACK❗✔️: *>
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️:       SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          UPLO
    // LAPACK❗✔️:       INTEGER            LDA, LDW, N, NB
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE, HALF
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       INTEGER            I, IW
    // LAPACK❗✔️:       DOUBLE PRECISION   ALPHA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME
    // LAPACK❗✔️:       DOUBLE PRECISION   DDOT
    // LAPACK❗✔️:       EXTERNAL           LSAME, DDOT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick return if possible
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.LE.0 )
    // LAPACK❗✔️:      $   RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( UPLO, 'U' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Reduce last NB columns of upper triangle
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 10 I = N, N - NB + 1, -1
    // LAPACK❗✔️:             IW = I - N + NB
    // LAPACK❗✔️:             IF( I.LT.N ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Update A(1:i,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
    // LAPACK❗✔️:      $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
    // LAPACK❗✔️:      $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( I.GT.1 ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Generate elementary reflector H(i) to annihilate
    // LAPACK❗✔️: *              A(1:i-2,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
    // LAPACK❗✔️:                E( I-1 ) = A( I-1, I )
    // LAPACK❗✔️:                A( I-1, I ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Compute W(1:i-1,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
    // LAPACK❗✔️:      $                     ZERO, W( 1, IW ), 1 )
    // LAPACK❗✔️:                IF( I.LT.N ) THEN
    // LAPACK❗✔️:                   CALL DGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
    // LAPACK❗✔️:      $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
    // LAPACK❗✔️:                   CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
    // LAPACK❗✔️:      $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
    // LAPACK❗✔️:      $                        W( 1, IW ), 1 )
    // LAPACK❗✔️:                   CALL DGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
    // LAPACK❗✔️:      $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
    // LAPACK❗✔️:                   CALL DGEMV( 'No transpose', I-1, N-I, -ONE,
    // LAPACK❗✔️:      $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
    // LAPACK❗✔️:      $                        W( 1, IW ), 1 )
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:                CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
    // LAPACK❗✔️:                ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
    // LAPACK❗✔️:      $                 A( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    10    CONTINUE
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Reduce first NB columns of lower triangle
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 20 I = 1, NB
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Update A(i:n,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
    // LAPACK❗✔️:      $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
    // LAPACK❗✔️:             CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
    // LAPACK❗✔️:      $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
    // LAPACK❗✔️:             IF( I.LT.N ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Generate elementary reflector H(i) to annihilate
    // LAPACK❗✔️: *              A(i+2:n,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
    // LAPACK❗✔️:      $                      TAU( I ) )
    // LAPACK❗✔️:                E( I ) = A( I+1, I )
    // LAPACK❗✔️:                A( I+1, I ) = ONE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *              Compute W(i+1:n,i)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
    // LAPACK❗✔️:      $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
    // LAPACK❗✔️:                CALL DGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
    // LAPACK❗✔️:      $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
    // LAPACK❗✔️:      $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
    // LAPACK❗✔️:                CALL DGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
    // LAPACK❗✔️:      $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
    // LAPACK❗✔️:      $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
    // LAPACK❗✔️:                CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
    // LAPACK❗✔️:                ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
    // LAPACK❗✔️:      $                 A( I+1, I ), 1 )
    // LAPACK❗✔️:                CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLATRD
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    for i in 0..nb {
        // Prior panel columns and the current column are disjoint. The
        // source's strided views need no full-matrix copy or aliasing slices.
        let (prior, current) = a.split_at_mut(i * lda);
        if i != 0 {
            blas::gemv(
                false,
                n - i,
                i,
                -1.0,
                &prior[i..],
                lda,
                &w[i..],
                ldw,
                1.0,
                &mut current[i..n],
                1,
            );
            blas::gemv(
                false,
                n - i,
                i,
                -1.0,
                &w[i..],
                ldw,
                &prior[i..],
                lda,
                1.0,
                &mut current[i..n],
                1,
            );
        }
        if i + 1 < n {
            let m = n - i - 1;
            tau[i] = larfg(&mut current[i + 1..n]);
            e[i] = current[i + 1];
            current[i + 1] = 1.0;
            let v = &current[i + 1..n];
            let future = &current[i + 1 + lda..];
            let (prior_w, current_w) = w.split_at_mut(i * ldw);
            let (head, tail) = current_w.split_at_mut(i + 1);
            blas::symv_lower(m, 1.0, future, lda, v, 1, 0.0, &mut tail[..m], 1);
            if i != 0 {
                blas::gemv(
                    true,
                    m,
                    i,
                    1.0,
                    &prior_w[i + 1..],
                    ldw,
                    v,
                    1,
                    0.0,
                    &mut head[..i],
                    1,
                );
                blas::gemv(
                    false,
                    m,
                    i,
                    -1.0,
                    &prior[i + 1..],
                    lda,
                    &head[..i],
                    1,
                    1.0,
                    &mut tail[..m],
                    1,
                );
                blas::gemv(
                    true,
                    m,
                    i,
                    1.0,
                    &prior[i + 1..],
                    lda,
                    v,
                    1,
                    0.0,
                    &mut head[..i],
                    1,
                );
                blas::gemv(
                    false,
                    m,
                    i,
                    -1.0,
                    &prior_w[i + 1..],
                    ldw,
                    &head[..i],
                    1,
                    1.0,
                    &mut tail[..m],
                    1,
                );
            }
            blas::scal(m, tau[i], &mut tail[..m], 1);
            let alpha = (-0.5 * tau[i]) * blas::dot(m, &tail[..m], 1, v, 1);
            blas::axpy(m, alpha, v, 1, &mut tail[..m], 1);
        }
    }
}

#[cfg(test)]
mod driver_tests {
    use super::*;

    #[test]
    fn dsyevd_lower_driver_small_and_preferred_blocked_workspace() {
        assert!(eigenvalues(0, Vec::new()).unwrap().is_empty());
        assert_eq!(
            eigenvalues(1, vec![-0.]).unwrap()[0].to_bits(),
            (-0.0_f64).to_bits()
        );
        // Analytic roots of [[3,4],[4,-3]] are exactly -5 and +5.
        assert_eq!(
            eigenvalues(2, vec![3., 4., 4., -3.]).unwrap(),
            vec![-5., 5.]
        );
        // Exercise the source's NB=NX=32 crossover without another eigensolver
        // as an oracle: a diagonal matrix's eigenvalues are its diagonal.
        for n in [32, 33, 65] {
            let mut a = vec![0.; n * n];
            for i in 0..n {
                a[i + i * n] = (n - 1 - i) as f64;
            }
            let expected: Vec<_> = (0..n).map(|i| i as f64).collect();
            assert_eq!(eigenvalues(n, a).unwrap(), expected);
        }
    }
}
