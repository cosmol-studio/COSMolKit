//! Real general eigenvalue kernels used by the default IPC polynomial.
//! Source: OpenBLAS 446c436e10450c348169808f1a6b3fae0925c9f7,
//! lapack-netlib/SRC (LAPACK BSD-3-Clause; notice in lapack_tridiagonal.rs).
//! Only eigenvalues are requested by the driver; internal Schur transformations
//! remain part of aggressive early deflation. Floating contraction/reduction
//! equivalence stays explicitly unverified until compared with the frozen build.

use super::lapack_tridiagonal::{dlamch, dlapy2};
mod hessenberg;
mod multishift;
mod schur;
mod swaps;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum GeneralError {
    InvalidArgument { routine: &'static str, argument: u8 },
    NonConvergence { info: usize },
    Blas(super::lapack_blas_level3::Blas3Error),
    Scalar(super::lapack_tridiagonal::DsterfError),
    DimensionOverflow,
}
impl From<super::lapack_blas_level3::Blas3Error> for GeneralError {
    fn from(value: super::lapack_blas_level3::Blas3Error) -> Self {
        Self::Blas(value)
    }
}

impl From<super::lapack_tridiagonal::DsterfError> for GeneralError {
    fn from(value: super::lapack_tridiagonal::DsterfError) -> Self {
        Self::Scalar(value)
    }
}
impl std::fmt::Display for GeneralError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidArgument { routine, argument } => {
                write!(f, "{routine} argument {argument} is invalid")
            }
            Self::NonConvergence { info } => {
                write!(f, "DGEEV failed to converge at eigenvalue {info}")
            }
            Self::Blas(error) => std::fmt::Display::fmt(error, f),
            Self::Scalar(error) => std::fmt::Display::fmt(error, f),
            Self::DimensionOverflow => f.write_str("DGEEV dimensions overflow the index range"),
        }
    }
}
impl std::error::Error for GeneralError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Blas(error) => Some(error),
            Self::Scalar(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Schur2 {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    real: [f64; 2],
    imag: [f64; 2],
    cs: f64,
    sn: f64,
}

fn lanv2(mut a: f64, mut b: f64, mut c: f64, mut d: f64) -> Schur2 {
    // LAPACK❗✔️:       SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, HALF, ONE, TWO
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0,
    // LAPACK❗✔️:      $                     TWO = 2.0D0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   MULTPL
    // LAPACK❗✔️:       PARAMETER          ( MULTPL = 4.0D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       DOUBLE PRECISION   AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB,
    // LAPACK❗✔️:      $                   SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z, SAFMIN,
    // LAPACK❗✔️:      $                   SAFMN2, SAFMX2
    // LAPACK❗✔️:       INTEGER            COUNT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DLAPY2
    // LAPACK❗✔️:       EXTERNAL           DLAMCH, DLAPY2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SAFMIN = DLAMCH( 'S' )
    // LAPACK❗✔️:       EPS = DLAMCH( 'P' )
    // LAPACK❗✔️:       SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
    // LAPACK❗✔️:      $            LOG( DLAMCH( 'B' ) ) / TWO )
    // LAPACK❗✔️:       SAFMX2 = ONE / SAFMN2
    // LAPACK❗✔️:       IF( C.EQ.ZERO ) THEN
    // LAPACK❗✔️:          CS = ONE
    // LAPACK❗✔️:          SN = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( B.EQ.ZERO ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Swap rows and columns
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CS = ZERO
    // LAPACK❗✔️:          SN = ONE
    // LAPACK❗✔️:          TEMP = D
    // LAPACK❗✔️:          D = A
    // LAPACK❗✔️:          A = TEMP
    // LAPACK❗✔️:          B = -C
    // LAPACK❗✔️:          C = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( ( A-D ).EQ.ZERO .AND. SIGN( ONE, B ).NE.SIGN( ONE, C ) )
    // LAPACK❗✔️:      $          THEN
    // LAPACK❗✔️:          CS = ONE
    // LAPACK❗✔️:          SN = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          TEMP = A - D
    // LAPACK❗✔️:          P = HALF*TEMP
    // LAPACK❗✔️:          BCMAX = MAX( ABS( B ), ABS( C ) )
    // LAPACK❗✔️:          BCMIS = MIN( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C )
    // LAPACK❗✔️:          SCALE = MAX( ABS( P ), BCMAX )
    // LAPACK❗✔️:          Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        If Z is of the order of the machine accuracy, postpone the
    // LAPACK❗✔️: *        decision on the nature of eigenvalues
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( Z.GE.MULTPL*EPS ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Real eigenvalues. Compute A and D.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             Z = P + SIGN( SQRT( SCALE )*SQRT( Z ), P )
    // LAPACK❗✔️:             A = D + Z
    // LAPACK❗✔️:             D = D - ( BCMAX / Z )*BCMIS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Compute B and the rotation matrix
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             TAU = DLAPY2( C, Z )
    // LAPACK❗✔️:             CS = Z / TAU
    // LAPACK❗✔️:             SN = C / TAU
    // LAPACK❗✔️:             B = B - C
    // LAPACK❗✔️:             C = ZERO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Complex eigenvalues, or real (almost) equal eigenvalues.
    // LAPACK❗✔️: *           Make diagonal elements equal.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             COUNT = 0
    // LAPACK❗✔️:             SIGMA = B + C
    // LAPACK❗✔️:    10       CONTINUE
    // LAPACK❗✔️:             COUNT = COUNT + 1
    // LAPACK❗✔️:             SCALE = MAX( ABS(TEMP), ABS(SIGMA) )
    // LAPACK❗✔️:             IF( SCALE.GE.SAFMX2 ) THEN
    // LAPACK❗✔️:                SIGMA = SIGMA * SAFMN2
    // LAPACK❗✔️:                TEMP = TEMP * SAFMN2
    // LAPACK❗✔️:                IF (COUNT .LE. 20)
    // LAPACK❗✔️:      $            GOTO 10
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( SCALE.LE.SAFMN2 ) THEN
    // LAPACK❗✔️:                SIGMA = SIGMA * SAFMX2
    // LAPACK❗✔️:                TEMP = TEMP * SAFMX2
    // LAPACK❗✔️:                IF (COUNT .LE. 20)
    // LAPACK❗✔️:      $            GOTO 10
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             P = HALF*TEMP
    // LAPACK❗✔️:             TAU = DLAPY2( SIGMA, TEMP )
    // LAPACK❗✔️:             CS = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
    // LAPACK❗✔️:             SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
    // LAPACK❗✔️: *                   [ CC  DD ]   [ C  D ] [ SN  CS ]
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             AA = A*CS + B*SN
    // LAPACK❗✔️:             BB = -A*SN + B*CS
    // LAPACK❗✔️:             CC = C*CS + D*SN
    // LAPACK❗✔️:             DD = -C*SN + D*CS
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
    // LAPACK❗✔️: *                   [ C  D ]   [-SN  CS ] [ CC  DD ]
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Note: Some of the multiplications are wrapped in parentheses to
    // LAPACK❗✔️: *                 prevent compilers from using FMA instructions. See
    // LAPACK❗✔️: *                 https://github.com/Reference-LAPACK/lapack/issues/1031.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             A = AA*CS + CC*SN
    // LAPACK❗✔️:             B = ( BB*CS ) + ( DD*SN )
    // LAPACK❗✔️:             C = -( AA*SN ) + ( CC*CS )
    // LAPACK❗✔️:             D = -BB*SN + DD*CS
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             TEMP = HALF*( A+D )
    // LAPACK❗✔️:             A = TEMP
    // LAPACK❗✔️:             D = TEMP
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( C.NE.ZERO ) THEN
    // LAPACK❗✔️:                IF( B.NE.ZERO ) THEN
    // LAPACK❗✔️:                   IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *                    Real eigenvalues: reduce to upper triangular form
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                      SAB = SQRT( ABS( B ) )
    // LAPACK❗✔️:                      SAC = SQRT( ABS( C ) )
    // LAPACK❗✔️:                      P = SIGN( SAB*SAC, C )
    // LAPACK❗✔️:                      TAU = ONE / SQRT( ABS( B+C ) )
    // LAPACK❗✔️:                      A = TEMP + P
    // LAPACK❗✔️:                      D = TEMP - P
    // LAPACK❗✔️:                      B = B - C
    // LAPACK❗✔️:                      C = ZERO
    // LAPACK❗✔️:                      CS1 = SAB*TAU
    // LAPACK❗✔️:                      SN1 = SAC*TAU
    // LAPACK❗✔️:                      TEMP = CS*CS1 - SN*SN1
    // LAPACK❗✔️:                      SN = CS*SN1 + SN*CS1
    // LAPACK❗✔️:                      CS = TEMP
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                ELSE
    // LAPACK❗✔️:                   B = -C
    // LAPACK❗✔️:                   C = ZERO
    // LAPACK❗✔️:                   TEMP = CS
    // LAPACK❗✔️:                   CS = -SN
    // LAPACK❗✔️:                   SN = TEMP
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RT1R = A
    // LAPACK❗✔️:       RT2R = D
    // LAPACK❗✔️:       IF( C.EQ.ZERO ) THEN
    // LAPACK❗✔️:          RT1I = ZERO
    // LAPACK❗✔️:          RT2I = ZERO
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️:          RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
    // LAPACK❗✔️:          RT2I = -RT1I
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DLANV2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    let eps = dlamch(b'P');
    let safmin = dlamch(b'S');
    let base = dlamch(b'B');
    let exponent = ((safmin / eps).ln() / base.ln() / 2.0) as i32;
    let safmn2 = base.powi(exponent);
    let safmx2 = 1.0 / safmn2;
    let (mut cs, mut sn);
    if c == 0.0 {
        cs = 1.0;
        sn = 0.0;
    } else if b == 0.0 {
        cs = 0.0;
        sn = 1.0;
        std::mem::swap(&mut a, &mut d);
        b = -c;
        c = 0.0;
    } else if a - d == 0.0 && 1.0_f64.copysign(b) != 1.0_f64.copysign(c) {
        cs = 1.0;
        sn = 0.0;
    } else {
        let mut temp = a - d;
        let mut p = 0.5 * temp;
        let bcmax = b.abs().max(c.abs());
        let bcmis = b.abs().min(c.abs()) * 1.0_f64.copysign(b) * 1.0_f64.copysign(c);
        let mut scale = p.abs().max(bcmax);
        let mut z = (p / scale) * p + (bcmax / scale) * bcmis;
        if z >= 4.0 * eps {
            z = p + (scale.sqrt() * z.sqrt()).copysign(p);
            a = d + z;
            d -= (bcmax / z) * bcmis;
            let tau = dlapy2(c, z);
            cs = z / tau;
            sn = c / tau;
            b -= c;
            c = 0.0;
        } else {
            let mut count = 0;
            let mut sigma = b + c;
            loop {
                count += 1;
                scale = temp.abs().max(sigma.abs());
                if scale >= safmx2 {
                    sigma *= safmn2;
                    temp *= safmn2;
                    if count <= 20 {
                        continue;
                    }
                }
                if scale <= safmn2 {
                    sigma *= safmx2;
                    temp *= safmx2;
                    if count <= 20 {
                        continue;
                    }
                }
                break;
            }
            p = 0.5 * temp;
            let tau = dlapy2(sigma, temp);
            cs = (0.5 * (1.0 + sigma.abs() / tau)).sqrt();
            sn = -(p / (tau * cs)) * 1.0_f64.copysign(sigma);
            let aa = a * cs + b * sn;
            let bb = -a * sn + b * cs;
            let cc = c * cs + d * sn;
            let dd = -c * sn + d * cs;
            a = aa * cs + cc * sn;
            // Source parentheses explicitly inhibit contraction here (LAPACK #1031).
            b = (bb * cs) + (dd * sn);
            c = -(aa * sn) + (cc * cs);
            d = -bb * sn + dd * cs;
            temp = 0.5 * (a + d);
            a = temp;
            d = temp;
            if c != 0.0 {
                if b != 0.0 {
                    if 1.0_f64.copysign(b) == 1.0_f64.copysign(c) {
                        let sab = b.abs().sqrt();
                        let sac = c.abs().sqrt();
                        p = (sab * sac).copysign(c);
                        let tau = 1.0 / (b + c).abs().sqrt();
                        a = temp + p;
                        d = temp - p;
                        b -= c;
                        c = 0.0;
                        let cs1 = sab * tau;
                        let sn1 = sac * tau;
                        temp = cs * cs1 - sn * sn1;
                        sn = cs * sn1 + sn * cs1;
                        cs = temp;
                    }
                } else {
                    b = -c;
                    c = 0.0;
                    temp = cs;
                    cs = -sn;
                    sn = temp;
                }
            }
        }
    }
    let imag = if c == 0.0 {
        [0.0, 0.0]
    } else {
        let value = b.abs().sqrt() * c.abs().sqrt();
        [value, -value]
    };
    Schur2 {
        a,
        b,
        c,
        d,
        real: [a, d],
        imag,
        cs,
        sn,
    }
}

pub(super) fn lartg(f: f64, g: f64) -> (f64, f64, f64) {
    // LAPACK❗✔️: subroutine DLARTG( f, g, c, s, r )
    // LAPACK❗✔️:    use LA_CONSTANTS, &
    // LAPACK❗✔️:    only: wp=>dp, zero=>dzero, half=>dhalf, one=>done, &
    // LAPACK❗✔️:          safmin=>dsafmin, safmax=>dsafmax
    // LAPACK❗✔️: !
    // LAPACK❗✔️: !  -- LAPACK auxiliary routine --
    // LAPACK❗✔️: !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: !     February 2021
    // LAPACK❗✔️: !
    // LAPACK❗✔️: !  .. Scalar Arguments ..
    // LAPACK❗✔️:    real(wp) :: c, f, g, r, s
    // LAPACK❗✔️: !  ..
    // LAPACK❗✔️: !  .. Local Scalars ..
    // LAPACK❗✔️:    real(wp) :: d, f1, fs, g1, gs, u, rtmin, rtmax
    // LAPACK❗✔️: !  ..
    // LAPACK❗✔️: !  .. Intrinsic Functions ..
    // LAPACK❗✔️:    intrinsic :: abs, sign, sqrt
    // LAPACK❗✔️: !  ..
    // LAPACK❗✔️: !  .. Constants ..
    // LAPACK❗✔️:    rtmin = sqrt( safmin )
    // LAPACK❗✔️:    rtmax = sqrt( safmax/2 )
    // LAPACK❗✔️: !  ..
    // LAPACK❗✔️: !  .. Executable Statements ..
    // LAPACK❗✔️: !
    // LAPACK❗✔️:    f1 = abs( f )
    // LAPACK❗✔️:    g1 = abs( g )
    // LAPACK❗✔️:    if( g == zero ) then
    // LAPACK❗✔️:       c = one
    // LAPACK❗✔️:       s = zero
    // LAPACK❗✔️:       r = f
    // LAPACK❗✔️:    else if( f == zero ) then
    // LAPACK❗✔️:       c = zero
    // LAPACK❗✔️:       s = sign( one, g )
    // LAPACK❗✔️:       r = g1
    // LAPACK❗✔️:    else if( f1 > rtmin .and. f1 < rtmax .and. &
    // LAPACK❗✔️:             g1 > rtmin .and. g1 < rtmax ) then
    // LAPACK❗✔️:       d = sqrt( f*f + g*g )
    // LAPACK❗✔️:       c = f1 / d
    // LAPACK❗✔️:       r = sign( d, f )
    // LAPACK❗✔️:       s = g / r
    // LAPACK❗✔️:    else
    // LAPACK❗✔️:       u = min( safmax, max( safmin, f1, g1 ) )
    // LAPACK❗✔️:       fs = f / u
    // LAPACK❗✔️:       gs = g / u
    // LAPACK❗✔️:       d = sqrt( fs*fs + gs*gs )
    // LAPACK❗✔️:       c = abs( fs ) / d
    // LAPACK❗✔️:       r = sign( d, f )
    // LAPACK❗✔️:       s = gs / r
    // LAPACK❗✔️:       r = r*u
    // LAPACK❗✔️:    end if
    // LAPACK❗✔️:    return
    // LAPACK❗✔️: end subroutine
    let safmin = dlamch(b'S');
    let safmax = 1.0 / safmin;
    let rtmin = safmin.sqrt();
    let rtmax = (safmax / 2.0).sqrt();
    let f1 = f.abs();
    let g1 = g.abs();
    if g == 0.0 {
        (1.0, 0.0, f)
    } else if f == 0.0 {
        (0.0, 1.0_f64.copysign(g), g1)
    } else if f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax {
        let d = (f * f + g * g).sqrt();
        let r = d.copysign(f);
        (f1 / d, g / r, r)
    } else {
        let u = safmax.min(safmin.max(f1).max(g1));
        let fs = f / u;
        let gs = g / u;
        let d = (fs * fs + gs * gs).sqrt();
        let r = d.copysign(f);
        (fs.abs() / d, gs / r, r * u)
    }
}

#[inline]
fn ix(i: usize, j: usize, ld: usize) -> usize {
    (j - 1) * ld + i - 1
}

fn swap_within(a: &mut [f64], n: usize, x: usize, incx: usize, y: usize, incy: usize) {
    // OpenBLAS swap_sse2.S/SWAP: values are exchanged without floating arithmetic.
    // The two ranges may share their identical first element, as in Fortran.
    for i in 0..n {
        a.swap(x + i * incx, y + i * incy);
    }
}

fn iamax(n: usize, x: &[f64], incx: usize) -> usize {
    // OpenBLAS interface/imax.c and iamax_sse2.S: return the first index
    // attaining the largest absolute value (one-based LAPACK convention).
    if n == 0 || incx == 0 {
        return 0;
    }
    let mut result = 0;
    let mut best = x[0].abs();
    for i in 1..n {
        let value = x[i * incx].abs();
        if value > best {
            best = value;
            result = i;
        }
    }
    result + 1
}

/// The JOB='B' path used by DGEEV, retaining one-based active block bounds.
fn gebal(n: usize, a: &mut [f64], scale: &mut [f64]) -> Result<(usize, usize), isize> {
    // LAPACK❗✔️:       SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK computational routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          JOB
    // LAPACK❗✔️:       INTEGER            IHI, ILO, INFO, LDA, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), SCALE( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   SCLFAC
    // LAPACK❗✔️:       PARAMETER          ( SCLFAC = 2.0D+0 )
    // LAPACK❗✔️:       DOUBLE PRECISION   FACTOR
    // LAPACK❗✔️:       PARAMETER          ( FACTOR = 0.95D+0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            NOCONV, CANSWAP
    // LAPACK❗✔️:       INTEGER            I, ICA, IRA, J, K, L
    // LAPACK❗✔️:       DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1,
    // LAPACK❗✔️:      $                   SFMIN2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            DISNAN, LSAME
    // LAPACK❗✔️:       INTEGER            IDAMAX
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DNRM2
    // LAPACK❗✔️:       EXTERNAL           DISNAN, LSAME, IDAMAX, DLAMCH, DNRM2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DSCAL, DSWAP, XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          ABS, MAX, MIN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     Test the input parameters
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
    // LAPACK❗✔️:      $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DGEBAL', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Quick returns.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( N.EQ.0 ) THEN
    // LAPACK❗✔️:          ILO = 1
    // LAPACK❗✔️:          IHI = 0
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( JOB, 'N' ) ) THEN
    // LAPACK❗✔️:          DO I = 1, N
    // LAPACK❗✔️:             SCALE( I ) = ONE
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:          ILO = 1
    // LAPACK❗✔️:          IHI = N
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Permutation to isolate eigenvalues if possible.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       K = 1
    // LAPACK❗✔️:       L = N
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( .NOT.LSAME( JOB, 'S' ) ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Row and column exchange.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          NOCONV = .TRUE.
    // LAPACK❗✔️:          DO WHILE( NOCONV )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Search for rows isolating an eigenvalue and push them down.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NOCONV = .FALSE.
    // LAPACK❗✔️:             DO I = L, 1, -1
    // LAPACK❗✔️:                CANSWAP = .TRUE.
    // LAPACK❗✔️:                DO J = 1, L
    // LAPACK❗✔️:                   IF( I.NE.J .AND. A( I, J ).NE.ZERO ) THEN
    // LAPACK❗✔️:                      CANSWAP = .FALSE.
    // LAPACK❗✔️:                      EXIT
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( CANSWAP ) THEN
    // LAPACK❗✔️:                   SCALE( L ) = I
    // LAPACK❗✔️:                   IF( I.NE.L ) THEN
    // LAPACK❗✔️:                      CALL DSWAP( L, A( 1, I ), 1, A( 1, L ), 1 )
    // LAPACK❗✔️:                      CALL DSWAP( N-K+1, A( I, K ), LDA, A( L, K ), LDA )
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                   NOCONV = .TRUE.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   IF( L.EQ.1 ) THEN
    // LAPACK❗✔️:                      ILO = 1
    // LAPACK❗✔️:                      IHI = 1
    // LAPACK❗✔️:                      RETURN
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   L = L - 1
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️:
    // LAPACK❗✔️:          NOCONV = .TRUE.
    // LAPACK❗✔️:          DO WHILE( NOCONV )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Search for columns isolating an eigenvalue and push them left.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             NOCONV = .FALSE.
    // LAPACK❗✔️:             DO J = K, L
    // LAPACK❗✔️:                CANSWAP = .TRUE.
    // LAPACK❗✔️:                DO I = K, L
    // LAPACK❗✔️:                   IF( I.NE.J .AND. A( I, J ).NE.ZERO ) THEN
    // LAPACK❗✔️:                      CANSWAP = .FALSE.
    // LAPACK❗✔️:                      EXIT
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                IF( CANSWAP ) THEN
    // LAPACK❗✔️:                   SCALE( K ) = J
    // LAPACK❗✔️:                   IF( J.NE.K ) THEN
    // LAPACK❗✔️:                      CALL DSWAP( L, A( 1, J ), 1, A( 1, K ), 1 )
    // LAPACK❗✔️:                      CALL DSWAP( N-K+1, A( J, K ), LDA, A( K, K ), LDA )
    // LAPACK❗✔️:                   END IF
    // LAPACK❗✔️:                   NOCONV = .TRUE.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:                   K = K + 1
    // LAPACK❗✔️:                END IF
    // LAPACK❗✔️:             END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Initialize SCALE for non-permuted submatrix.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       DO I = K, L
    // LAPACK❗✔️:          SCALE( I ) = ONE
    // LAPACK❗✔️:       END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     If we only had to permute, we are done.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( LSAME( JOB, 'P' ) ) THEN
    // LAPACK❗✔️:          ILO = K
    // LAPACK❗✔️:          IHI = L
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Balance the submatrix in rows K to L.
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Iterative loop for norm reduction.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
    // LAPACK❗✔️:       SFMAX1 = ONE / SFMIN1
    // LAPACK❗✔️:       SFMIN2 = SFMIN1*SCLFAC
    // LAPACK❗✔️:       SFMAX2 = ONE / SFMIN2
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       NOCONV = .TRUE.
    // LAPACK❗✔️:       DO WHILE( NOCONV )
    // LAPACK❗✔️:          NOCONV = .FALSE.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO I = K, L
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             C = DNRM2( L-K+1, A( K, I ), 1 )
    // LAPACK❗✔️:             R = DNRM2( L-K+1, A( I, K ), LDA )
    // LAPACK❗✔️:             ICA = IDAMAX( L, A( 1, I ), 1 )
    // LAPACK❗✔️:             CA = ABS( A( ICA, I ) )
    // LAPACK❗✔️:             IRA = IDAMAX( N-K+1, A( I, K ), LDA )
    // LAPACK❗✔️:             RA = ABS( A( I, IRA+K-1 ) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Guard against zero C or R due to underflow.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( C.EQ.ZERO .OR. R.EQ.ZERO ) CYCLE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Exit if NaN to avoid infinite loop
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( DISNAN( C+CA+R+RA ) ) THEN
    // LAPACK❗✔️:                INFO = -3
    // LAPACK❗✔️:                CALL XERBLA( 'DGEBAL', -INFO )
    // LAPACK❗✔️:                RETURN
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             G = R / SCLFAC
    // LAPACK❗✔️:             F = ONE
    // LAPACK❗✔️:             S = C + R
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO WHILE( C.LT.G .AND. MAX( F, C, CA ).LT.SFMAX2 .AND.
    // LAPACK❗✔️:      $                MIN( R, G, RA ).GT.SFMIN2 )
    // LAPACK❗✔️:                F = F*SCLFAC
    // LAPACK❗✔️:                C = C*SCLFAC
    // LAPACK❗✔️:                CA = CA*SCLFAC
    // LAPACK❗✔️:                R = R / SCLFAC
    // LAPACK❗✔️:                G = G / SCLFAC
    // LAPACK❗✔️:                RA = RA / SCLFAC
    // LAPACK❗✔️:             END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             G = C / SCLFAC
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             DO WHILE( G.GE.R .AND. MAX( R, RA ).LT.SFMAX2 .AND.
    // LAPACK❗✔️:      $                MIN( F, C, G, CA ).GT.SFMIN2 )
    // LAPACK❗✔️:                F = F / SCLFAC
    // LAPACK❗✔️:                C = C / SCLFAC
    // LAPACK❗✔️:                G = G / SCLFAC
    // LAPACK❗✔️:                CA = CA / SCLFAC
    // LAPACK❗✔️:                R = R*SCLFAC
    // LAPACK❗✔️:                RA = RA*SCLFAC
    // LAPACK❗✔️:             END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Now balance.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             IF( ( C+R ).GE.FACTOR*S ) CYCLE
    // LAPACK❗✔️:             IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
    // LAPACK❗✔️:                IF( F*SCALE( I ).LE.SFMIN1 ) CYCLE
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
    // LAPACK❗✔️:                IF( SCALE( I ).GE.SFMAX1 / F ) CYCLE
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             G = ONE / F
    // LAPACK❗✔️:             SCALE( I ) = SCALE( I )*F
    // LAPACK❗✔️:             NOCONV = .TRUE.
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             CALL DSCAL( N-K+1, G, A( I, K ), LDA )
    // LAPACK❗✔️:             CALL DSCAL( L, F, A( 1, I ), 1 )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END DO
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ILO = K
    // LAPACK❗✔️:       IHI = L
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DGEBAL
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END
    if n == 0 {
        return Ok((1, 0));
    }
    let mut k = 1;
    let mut l = n;
    loop {
        let mut changed = false;
        for i in (1..=l).rev() {
            let can_swap = (1..=l).all(|j| i == j || a[ix(i, j, n)] == 0.0);
            if can_swap {
                scale[l - 1] = i as f64;
                if i != l {
                    swap_within(a, l, ix(1, i, n), 1, ix(1, l, n), 1);
                    swap_within(a, n - k + 1, ix(i, k, n), n, ix(l, k, n), n);
                }
                changed = true;
                if l == 1 {
                    return Ok((1, 1));
                }
                l -= 1;
            }
        }
        if !changed {
            break;
        }
    }
    loop {
        let mut changed = false;
        for j in k..=l {
            let can_swap = (k..=l).all(|i| i == j || a[ix(i, j, n)] == 0.0);
            if can_swap {
                scale[k - 1] = j as f64;
                if j != k {
                    swap_within(a, l, ix(1, j, n), 1, ix(1, k, n), 1);
                    swap_within(a, n - k + 1, ix(j, k, n), n, ix(k, k, n), n);
                }
                changed = true;
                k += 1;
            }
        }
        if !changed {
            break;
        }
    }
    scale[k - 1..l].fill(1.0);
    let sfmin1 = dlamch(b'S') / dlamch(b'P');
    let sfmax1 = 1.0 / sfmin1;
    let sfmin2 = sfmin1 * 2.0;
    let sfmax2 = 1.0 / sfmin2;
    loop {
        let mut changed = false;
        for i in k..=l {
            let mut c = super::lapack_blas::nrm2(l - k + 1, &a[ix(k, i, n)..], 1);
            let mut r = super::lapack_blas::nrm2(l - k + 1, &a[ix(i, k, n)..], n);
            let ica = iamax(l, &a[ix(1, i, n)..], 1);
            let mut ca = a[ix(ica, i, n)].abs();
            let ira = iamax(n - k + 1, &a[ix(i, k, n)..], n);
            let mut ra = a[ix(i, ira + k - 1, n)].abs();
            if c == 0.0 || r == 0.0 {
                continue;
            }
            if (c + ca + r + ra).is_nan() {
                return Err(-3);
            }
            let mut g = r / 2.0;
            let mut f: f64 = 1.0;
            let sum = c + r;
            while c < g && f.max(c).max(ca) < sfmax2 && r.min(g).min(ra) > sfmin2 {
                f *= 2.0;
                c *= 2.0;
                ca *= 2.0;
                r /= 2.0;
                g /= 2.0;
                ra /= 2.0;
            }
            g = c / 2.0;
            while g >= r && r.max(ra) < sfmax2 && f.min(c).min(g).min(ca) > sfmin2 {
                f /= 2.0;
                c /= 2.0;
                g /= 2.0;
                ca /= 2.0;
                r *= 2.0;
                ra *= 2.0;
            }
            if c + r >= 0.95 * sum {
                continue;
            }
            if f < 1.0 && scale[i - 1] < 1.0 && f * scale[i - 1] <= sfmin1 {
                continue;
            }
            if f > 1.0 && scale[i - 1] > 1.0 && scale[i - 1] >= sfmax1 / f {
                continue;
            }
            g = 1.0 / f;
            scale[i - 1] *= f;
            changed = true;
            super::lapack_blas::scal(n - k + 1, g, &mut a[ix(i, k, n)..], n);
            super::lapack_blas::scal(l, f, &mut a[ix(1, i, n)..], 1);
        }
        if !changed {
            break;
        }
    }
    Ok((k, l))
}

/// DGEEV JOBVL='N', JOBVR='N'; accepts the caller-owned column-major input.
/// Roots retain LAPACK's order and positive/negative conjugate-pair order.
pub(super) fn eigenvalues(n: usize, mut a: Vec<f64>) -> Result<Vec<[f64; 2]>, GeneralError> {
    // LAPACK❗✔️:       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
    // LAPACK❗✔️:      $                  LDVR, WORK, LWORK, INFO )
    // LAPACK❗✔️:       implicit none
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  -- LAPACK driver routine --
    // LAPACK❗✔️: *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    // LAPACK❗✔️: *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Scalar Arguments ..
    // LAPACK❗✔️:       CHARACTER          JOBVL, JOBVR
    // LAPACK❗✔️:       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Array Arguments ..
    // LAPACK❗✔️:       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
    // LAPACK❗✔️:      $                   WI( * ), WORK( * ), WR( * )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *  =====================================================================
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     .. Parameters ..
    // LAPACK❗✔️:       DOUBLE PRECISION   ZERO, ONE
    // LAPACK❗✔️:       PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Scalars ..
    // LAPACK❗✔️:       LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
    // LAPACK❗✔️:       CHARACTER          SIDE
    // LAPACK❗✔️:       INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K,
    // LAPACK❗✔️:      $                   LWORK_TREVC, MAXWRK, MINWRK, NOUT
    // LAPACK❗✔️:       DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM,
    // LAPACK❗✔️:      $                   SN
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Local Arrays ..
    // LAPACK❗✔️:       LOGICAL            SELECT( 1 )
    // LAPACK❗✔️:       DOUBLE PRECISION   DUM( 1 )
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Subroutines ..
    // LAPACK❗✔️:       EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLABAD, DLACPY,
    // LAPACK❗✔️:      $                   DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC3,
    // LAPACK❗✔️:      $                   XERBLA
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. External Functions ..
    // LAPACK❗✔️:       LOGICAL            LSAME, DISNAN
    // LAPACK❗✔️:       INTEGER            IDAMAX, ILAENV
    // LAPACK❗✔️:       DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
    // LAPACK❗✔️:       EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2,
    // LAPACK❗✔️:      $                   DNRM2
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Intrinsic Functions ..
    // LAPACK❗✔️:       INTRINSIC          MAX, SQRT
    // LAPACK❗✔️: *     ..
    // LAPACK❗✔️: *     .. Executable Statements ..
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Test the input arguments
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       INFO = 0
    // LAPACK❗✔️:       LQUERY = ( LWORK.EQ.-1 )
    // LAPACK❗✔️:       WANTVL = LSAME( JOBVL, 'V' )
    // LAPACK❗✔️:       WANTVR = LSAME( JOBVR, 'V' )
    // LAPACK❗✔️:       IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
    // LAPACK❗✔️:          INFO = -1
    // LAPACK❗✔️:       ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
    // LAPACK❗✔️:          INFO = -2
    // LAPACK❗✔️:       ELSE IF( N.LT.0 ) THEN
    // LAPACK❗✔️:          INFO = -3
    // LAPACK❗✔️:       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
    // LAPACK❗✔️:          INFO = -5
    // LAPACK❗✔️:       ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
    // LAPACK❗✔️:          INFO = -9
    // LAPACK❗✔️:       ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
    // LAPACK❗✔️:          INFO = -11
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Compute workspace
    // LAPACK❗✔️: *      (Note: Comments in the code beginning "Workspace:" describe the
    // LAPACK❗✔️: *       minimal amount of workspace needed at that point in the code,
    // LAPACK❗✔️: *       as well as the preferred amount for good performance.
    // LAPACK❗✔️: *       NB refers to the optimal block size for the immediately
    // LAPACK❗✔️: *       following subroutine, as returned by ILAENV.
    // LAPACK❗✔️: *       HSWORK refers to the workspace preferred by DHSEQR, as
    // LAPACK❗✔️: *       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
    // LAPACK❗✔️: *       the worst case.)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.EQ.0 ) THEN
    // LAPACK❗✔️:          IF( N.EQ.0 ) THEN
    // LAPACK❗✔️:             MINWRK = 1
    // LAPACK❗✔️:             MAXWRK = 1
    // LAPACK❗✔️:          ELSE
    // LAPACK❗✔️:             MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
    // LAPACK❗✔️:             IF( WANTVL ) THEN
    // LAPACK❗✔️:                MINWRK = 4*N
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1,
    // LAPACK❗✔️:      $                       'DORGHR', ' ', N, 1, N, -1 ) )
    // LAPACK❗✔️:                CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL,
    // LAPACK❗✔️:      $                      WORK, -1, INFO )
    // LAPACK❗✔️:                HSWORK = INT( WORK(1) )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
    // LAPACK❗✔️:                CALL DTREVC3( 'L', 'B', SELECT, N, A, LDA,
    // LAPACK❗✔️:      $                       VL, LDVL, VR, LDVR, N, NOUT,
    // LAPACK❗✔️:      $                       WORK, -1, IERR )
    // LAPACK❗✔️:                LWORK_TREVC = INT( WORK(1) )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, 4*N )
    // LAPACK❗✔️:             ELSE IF( WANTVR ) THEN
    // LAPACK❗✔️:                MINWRK = 4*N
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1,
    // LAPACK❗✔️:      $                       'DORGHR', ' ', N, 1, N, -1 ) )
    // LAPACK❗✔️:                CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR,
    // LAPACK❗✔️:      $                      WORK, -1, INFO )
    // LAPACK❗✔️:                HSWORK = INT( WORK(1) )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
    // LAPACK❗✔️:                CALL DTREVC3( 'R', 'B', SELECT, N, A, LDA,
    // LAPACK❗✔️:      $                       VL, LDVL, VR, LDVR, N, NOUT,
    // LAPACK❗✔️:      $                       WORK, -1, IERR )
    // LAPACK❗✔️:                LWORK_TREVC = INT( WORK(1) )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, 4*N )
    // LAPACK❗✔️:             ELSE
    // LAPACK❗✔️:                MINWRK = 3*N
    // LAPACK❗✔️:                CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR,
    // LAPACK❗✔️:      $                      WORK, -1, INFO )
    // LAPACK❗✔️:                HSWORK = INT( WORK(1) )
    // LAPACK❗✔️:                MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:             MAXWRK = MAX( MAXWRK, MINWRK )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:          WORK( 1 ) = MAXWRK
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
    // LAPACK❗✔️:             INFO = -13
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 ) THEN
    // LAPACK❗✔️:          CALL XERBLA( 'DGEEV ', -INFO )
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
    // LAPACK❗✔️: *     Get machine constants
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       EPS = DLAMCH( 'P' )
    // LAPACK❗✔️:       SMLNUM = DLAMCH( 'S' )
    // LAPACK❗✔️:       BIGNUM = ONE / SMLNUM
    // LAPACK❗✔️:       CALL DLABAD( SMLNUM, BIGNUM )
    // LAPACK❗✔️:       SMLNUM = SQRT( SMLNUM ) / EPS
    // LAPACK❗✔️:       BIGNUM = ONE / SMLNUM
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Scale A if max element outside range [SMLNUM,BIGNUM]
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
    // LAPACK❗✔️:       SCALEA = .FALSE.
    // LAPACK❗✔️:       IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
    // LAPACK❗✔️:          SCALEA = .TRUE.
    // LAPACK❗✔️:          CSCALE = SMLNUM
    // LAPACK❗✔️:       ELSE IF( ANRM.GT.BIGNUM ) THEN
    // LAPACK❗✔️:          SCALEA = .TRUE.
    // LAPACK❗✔️:          CSCALE = BIGNUM
    // LAPACK❗✔️:       ELSE IF ( DISNAN (ANRM ) ) THEN
    // LAPACK❗✔️:          INFO = -4
    // LAPACK❗✔️:          CALL XERBLA ( 'DGEEV', -INFO )
    // LAPACK❗✔️:          RETURN
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️:       IF( SCALEA )
    // LAPACK❗✔️:      $   CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Balance the matrix
    // LAPACK❗✔️: *     (Workspace: need N)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IBAL = 1
    // LAPACK❗✔️:       CALL DGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Reduce to upper Hessenberg form
    // LAPACK❗✔️: *     (Workspace: need 3*N, prefer 2*N+N*NB)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ITAU = IBAL + N
    // LAPACK❗✔️:       IWRK = ITAU + N
    // LAPACK❗✔️:       CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
    // LAPACK❗✔️:      $             LWORK-IWRK+1, IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( WANTVL ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Want left eigenvectors
    // LAPACK❗✔️: *        Copy Householder vectors to VL
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SIDE = 'L'
    // LAPACK❗✔️:          CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Generate orthogonal matrix in VL
    // LAPACK❗✔️: *        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ),
    // LAPACK❗✔️:      $                LWORK-IWRK+1, IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Perform QR iteration, accumulating Schur vectors in VL
    // LAPACK❗✔️: *        (Workspace: need N+1, prefer N+HSWORK (see comments) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IWRK = ITAU
    // LAPACK❗✔️:          CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL,
    // LAPACK❗✔️:      $                WORK( IWRK ), LWORK-IWRK+1, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IF( WANTVR ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *           Want left and right eigenvectors
    // LAPACK❗✔️: *           Copy Schur vectors to VR
    // LAPACK❗✔️: *
    // LAPACK❗✔️:             SIDE = 'B'
    // LAPACK❗✔️:             CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE IF( WANTVR ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Want right eigenvectors
    // LAPACK❗✔️: *        Copy Householder vectors to VR
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          SIDE = 'R'
    // LAPACK❗✔️:          CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Generate orthogonal matrix in VR
    // LAPACK❗✔️: *        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ),
    // LAPACK❗✔️:      $                LWORK-IWRK+1, IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Perform QR iteration, accumulating Schur vectors in VR
    // LAPACK❗✔️: *        (Workspace: need N+1, prefer N+HSWORK (see comments) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IWRK = ITAU
    // LAPACK❗✔️:          CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
    // LAPACK❗✔️:      $                WORK( IWRK ), LWORK-IWRK+1, INFO )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       ELSE
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute eigenvalues only
    // LAPACK❗✔️: *        (Workspace: need N+1, prefer N+HSWORK (see comments) )
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          IWRK = ITAU
    // LAPACK❗✔️:          CALL DHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
    // LAPACK❗✔️:      $                WORK( IWRK ), LWORK-IWRK+1, INFO )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     If INFO .NE. 0 from DHSEQR, then quit
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( INFO.NE.0 )
    // LAPACK❗✔️:      $   GO TO 50
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( WANTVL .OR. WANTVR ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Compute left and/or right eigenvectors
    // LAPACK❗✔️: *        (Workspace: need 4*N, prefer N + N + 2*N*NB)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DTREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
    // LAPACK❗✔️:      $                 N, NOUT, WORK( IWRK ), LWORK-IWRK+1, IERR )
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( WANTVL ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Undo balancing of left eigenvectors
    // LAPACK❗✔️: *        (Workspace: need N)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL,
    // LAPACK❗✔️:      $                IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Normalize left eigenvectors and make largest component real
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 20 I = 1, N
    // LAPACK❗✔️:             IF( WI( I ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DSCAL( N, SCL, VL( 1, I ), 1 )
    // LAPACK❗✔️:             ELSE IF( WI( I ).GT.ZERO ) THEN
    // LAPACK❗✔️:                SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ),
    // LAPACK❗✔️:      $               DNRM2( N, VL( 1, I+1 ), 1 ) )
    // LAPACK❗✔️:                CALL DSCAL( N, SCL, VL( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
    // LAPACK❗✔️:                DO 10 K = 1, N
    // LAPACK❗✔️:                   WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
    // LAPACK❗✔️:    10          CONTINUE
    // LAPACK❗✔️:                K = IDAMAX( N, WORK( IWRK ), 1 )
    // LAPACK❗✔️:                CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
    // LAPACK❗✔️:                CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
    // LAPACK❗✔️:                VL( K, I+1 ) = ZERO
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:    20    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       IF( WANTVR ) THEN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Undo balancing of right eigenvectors
    // LAPACK❗✔️: *        (Workspace: need N)
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          CALL DGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR,
    // LAPACK❗✔️:      $                IERR )
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *        Normalize right eigenvectors and make largest component real
    // LAPACK❗✔️: *
    // LAPACK❗✔️:          DO 40 I = 1, N
    // LAPACK❗✔️:             IF( WI( I ).EQ.ZERO ) THEN
    // LAPACK❗✔️:                SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DSCAL( N, SCL, VR( 1, I ), 1 )
    // LAPACK❗✔️:             ELSE IF( WI( I ).GT.ZERO ) THEN
    // LAPACK❗✔️:                SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ),
    // LAPACK❗✔️:      $               DNRM2( N, VR( 1, I+1 ), 1 ) )
    // LAPACK❗✔️:                CALL DSCAL( N, SCL, VR( 1, I ), 1 )
    // LAPACK❗✔️:                CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
    // LAPACK❗✔️:                DO 30 K = 1, N
    // LAPACK❗✔️:                   WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
    // LAPACK❗✔️:    30          CONTINUE
    // LAPACK❗✔️:                K = IDAMAX( N, WORK( IWRK ), 1 )
    // LAPACK❗✔️:                CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
    // LAPACK❗✔️:                CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
    // LAPACK❗✔️:                VR( K, I+1 ) = ZERO
    // LAPACK❗✔️:             END IF
    // LAPACK❗✔️:    40    CONTINUE
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     Undo scaling if necessary
    // LAPACK❗✔️: *
    // LAPACK❗✔️:    50 CONTINUE
    // LAPACK❗✔️:       IF( SCALEA ) THEN
    // LAPACK❗✔️:          CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ),
    // LAPACK❗✔️:      $                MAX( N-INFO, 1 ), IERR )
    // LAPACK❗✔️:          CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ),
    // LAPACK❗✔️:      $                MAX( N-INFO, 1 ), IERR )
    // LAPACK❗✔️:          IF( INFO.GT.0 ) THEN
    // LAPACK❗✔️:             CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N,
    // LAPACK❗✔️:      $                   IERR )
    // LAPACK❗✔️:             CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N,
    // LAPACK❗✔️:      $                   IERR )
    // LAPACK❗✔️:          END IF
    // LAPACK❗✔️:       END IF
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       WORK( 1 ) = MAXWRK
    // LAPACK❗✔️:       RETURN
    // LAPACK❗✔️: *
    // LAPACK❗✔️: *     End of DGEEV
    // LAPACK❗✔️: *
    // LAPACK❗✔️:       END

    let square = n.checked_mul(n).ok_or(GeneralError::DimensionOverflow)?;
    if a.len() != square {
        return Err(GeneralError::InvalidArgument {
            routine: "DGEEV",
            argument: 4,
        });
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    // Reproduce NumPy's driver workspace query, including DHSEQR's query.
    // Giving DGEHRD additional private workspace changes its selected NB.
    let hswork = n.max(multishift::qr_work(true, n, 1, n));
    let maxwrk = (2 * n + n * 32).max(n + 1).max(n + hswork).max(3 * n);
    let mut work = vec![0.0; maxwrk];
    let mut real = vec![0.0; n];
    let mut imag = vec![0.0; n];
    let eps = dlamch(b'P');
    let mut smlnum = dlamch(b'S');
    // DLABAD's exponent-range test is false on IEEE binary64.
    smlnum = smlnum.sqrt() / eps;
    let bignum = 1.0 / smlnum;
    let mut anrm: f64 = 0.0;
    for &value in &a {
        if value.is_nan() {
            return Err(GeneralError::InvalidArgument {
                routine: "DGEEV",
                argument: 4,
            });
        }
        anrm = anrm.max(value.abs());
    }
    let cscale = if anrm > 0.0 && anrm < smlnum {
        Some(smlnum)
    } else if anrm > bignum {
        Some(bignum)
    } else {
        None
    };
    if let Some(cscale) = cscale {
        super::lapack_tridiagonal::dlascl_vector(anrm, cscale, &mut a)?;
    }
    let (scale, remaining) = work.split_at_mut(n);
    let (ilo, ihi) = gebal(n, &mut a, scale).map_err(|info| GeneralError::InvalidArgument {
        routine: "DGEBAL",
        argument: (-info) as u8,
    })?;
    let (tau, scratch) = remaining.split_at_mut(n);
    hessenberg::gehrd(n, ilo, ihi, &mut a, n, tau, scratch)?;
    let info = schur::hseqr(n, ilo, ihi, &mut a, n, &mut real, &mut imag, remaining)?;
    if let Some(cscale) = cscale {
        super::lapack_tridiagonal::dlascl_vector(cscale, anrm, &mut real[info..])?;
        super::lapack_tridiagonal::dlascl_vector(cscale, anrm, &mut imag[info..])?;
        if info > 0 {
            super::lapack_tridiagonal::dlascl_vector(cscale, anrm, &mut real[..ilo - 1])?;
            super::lapack_tridiagonal::dlascl_vector(cscale, anrm, &mut imag[..ilo - 1])?;
        }
    }
    if info > 0 {
        return Err(GeneralError::NonConvergence { info });
    }
    Ok(real
        .into_iter()
        .zip(imag)
        .map(|(real, imag)| [real, imag])
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn general_schur_pair_and_safe_rotation() {
        let balanced = vec![2.0, 0.0, 0.0, 0.0, 0.0, -0.125, 0.0, 8.0, 0.0];
        assert_eq!(
            eigenvalues(3, balanced).unwrap(),
            [[0.0, 1.0], [0.0, -1.0], [2.0, 0.0]]
        );
        // Cross the fixed DGEHRD and DHSEQR blocked-driver boundaries. Each
        // irreducible 2x2 block has the exactly representable roots +/-i.
        let n = 150;
        let mut block_diagonal = vec![0.0; n * n];
        for i in (0..n).step_by(2) {
            block_diagonal[i + (i + 1) * n] = 1.0;
            block_diagonal[i + 1 + i * n] = -1.0;
        }
        let roots = eigenvalues(n, block_diagonal).unwrap();
        assert_eq!(roots.len(), n);
        for pair in roots.chunks_exact(2) {
            assert_eq!(pair, [[0.0, 1.0], [0.0, -1.0]]);
        }
        let real = lanv2(3.0, 2.0, 0.0, 1.0);
        assert_eq!(real.real, [3.0, 1.0]);
        assert_eq!(real.imag, [0.0, 0.0]);
        let complex = lanv2(0.0, 1.0, -1.0, 0.0);
        assert_eq!(complex.real, [0.0, 0.0]);
        assert_eq!(complex.imag, [1.0, -1.0]);
        let mut h = [2.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0];
        let mut real = [0.0; 3];
        let mut imag = [0.0; 3];
        assert_eq!(
            schur::lahqr(
                false,
                false,
                3,
                1,
                3,
                &mut h,
                3,
                &mut real,
                &mut imag,
                1,
                0,
                &mut [],
                1
            ),
            0
        );
        assert_eq!(real, [2.0, 0.0, 0.0]);
        assert_eq!(imag, [0.0, 1.0, -1.0]);
        let mut q = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let (mut from, mut to) = (1, 3);
        assert_eq!(swaps::trexc(3, &mut h, 3, &mut q, 3, &mut from, &mut to), 0);
        assert_eq!((h[0], h[4], h[8]), (0.0, 0.0, 2.0));
        from = 3;
        to = 1;
        assert_eq!(swaps::trexc(3, &mut h, 3, &mut q, 3, &mut from, &mut to), 0);
        assert_eq!((h[0], h[4], h[8]), (2.0, 0.0, 0.0));
        assert_eq!(lartg(3.0, 4.0), (0.6, 0.8, 5.0));
        assert_eq!(lartg(-0.0, 0.0).2.to_bits(), (-0.0_f64).to_bits());
        let (cs, sn, r) = lartg(f64::MIN_POSITIVE, f64::MIN_POSITIVE);
        assert!(cs.is_finite() && sn.is_finite() && r.is_normal());
    }
}
