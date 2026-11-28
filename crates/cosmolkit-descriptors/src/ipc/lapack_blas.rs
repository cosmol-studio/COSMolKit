//! Arithmetic kernels for the pinned single-threaded SkylakeX LAPACK path.
//! Matrices are column-major; vector slices start at their logical first element.
//! Explicit lane reductions and fused operations reproduce the reference binary.

// Source: OpenBLAS 446c436e10450c348169808f1a6b3fae0925c9f7.
// /*********************************************************************/
// /* Copyright 2009, 2010 The University of Texas at Austin.           */
// /* All rights reserved.                                              */
// /*                                                                   */
// /* Redistribution and use in source and binary forms, with or        */
// /* without modification, are permitted provided that the following   */
// /* conditions are met:                                               */
// /*                                                                   */
// /*   1. Redistributions of source code must retain the above         */
// /*      copyright notice, this list of conditions and the following  */
// /*      disclaimer.                                                  */
// /*                                                                   */
// /*   2. Redistributions in binary form must reproduce the above      */
// /*      copyright notice, this list of conditions and the following  */
// /*      disclaimer in the documentation and/or other materials       */
// /*      provided with the distribution.                              */
// /*                                                                   */
// /*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
// /*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
// /*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
// /*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
// /*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
// /*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
// /*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
// /*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
// /*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
// /*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
// /*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
// /*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
// /*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
// /*    POSSIBILITY OF SUCH DAMAGE.                                    */
// /*                                                                   */
// /* The views and conclusions contained in the software and           */
// /* documentation are those of the authors and should not be          */
// /* interpreted as representing official policies, either expressed   */
// /* or implied, of The University of Texas at Austin.                 */
// /*********************************************************************/
//
// #define ASSEMBLER

// /***************************************************************************
// Copyright (c) 2014, The OpenBLAS Project
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in
// the documentation and/or other materials provided with the
// distribution.
// 3. Neither the name of the OpenBLAS project nor the names of
// its contributors may be used to endorse or promote products
// derived from this software without specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// *****************************************************************************/
//

pub(super) use super::lapack_blas_level3::{gemm, syr2k_lower, trmm};

#[inline]
fn reduce4(a: [f64; 4]) -> f64 {
    (a[0] + a[2]) + (a[1] + a[3])
}

pub(super) fn dot(n: usize, x: &[f64], incx: usize, y: &[f64], incy: usize) -> f64 {
    // Source: kernel/x86_64/ddot.c + ddot_microk_skylakex-2.c.
    // OpenBLAS✔️❌: 	BLASLONG ix=0,iy=0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLOAT  dot = 0.0 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( n <= 0 )  return(dot);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( (inc_x == 1) && (inc_y == 1) )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		BLASLONG n1 = n & -16;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( n1 )
    // OpenBLAS✔️❌: 			ddot_kernel_8(n1, x, y , &dot );
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		i = n1;
    // OpenBLAS✔️❌: 		while(i < n)
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			dot += y[i] * x[i] ;
    // OpenBLAS✔️❌: 			i++ ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		return(dot);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLOAT temp1 = 0.0;
    // OpenBLAS✔️❌: 	FLOAT temp2 = 0.0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         BLASLONG n1 = n & -4;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	while(i < n1)
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		FLOAT m1 = y[iy]       * x[ix] ;
    // OpenBLAS✔️❌: 		FLOAT m2 = y[iy+inc_y] * x[ix+inc_x] ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		FLOAT m3 = y[iy+2*inc_y] * x[ix+2*inc_x] ;
    // OpenBLAS✔️❌: 		FLOAT m4 = y[iy+3*inc_y] * x[ix+3*inc_x] ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		ix  += inc_x*4 ;
    // OpenBLAS✔️❌: 		iy  += inc_y*4 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		temp1 += m1+m3;
    // OpenBLAS✔️❌: 		temp2 += m2+m4;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		i+=4 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	while(i < n)
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		temp1 += y[iy] * x[ix] ;
    // OpenBLAS✔️❌: 		ix  += inc_x ;
    // OpenBLAS✔️❌: 		iy  += inc_y ;
    // OpenBLAS✔️❌: 		i++ ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌: 	dot = temp1 + temp2;
    // OpenBLAS✔️❌: 	return(dot);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #if defined(SMP)
    // OpenBLAS✔️❌: static int dot_thread_function(BLASLONG n, BLASLONG dummy0,
    // OpenBLAS✔️❌:         BLASLONG dummy1, FLOAT dummy2, FLOAT *x, BLASLONG inc_x, FLOAT *y,
    // OpenBLAS✔️❌:         BLASLONG inc_y, RETURN_TYPE *result, BLASLONG dummy3)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:         *(RETURN_TYPE *)result = dot_compute(n, x, inc_x, y, inc_y);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         return 0;
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: extern int blas_level1_thread_with_return_value(int mode, BLASLONG m, BLASLONG n,
    // OpenBLAS✔️❌:         BLASLONG k, void *alpha, void *a, BLASLONG lda, void *b, BLASLONG ldb,
    // OpenBLAS✔️❌:         void *c, BLASLONG ldc, int (*function)(void), int nthreads);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: FLOAT CNAME(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌: 	accum_3 = _mm256_setzero_pd();
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef __AVX512CD__
    // OpenBLAS✔️❌: 	__m512d accum_05, accum_15, accum_25, accum_35;
    // OpenBLAS✔️❌: 	int n32;
    // OpenBLAS✔️❌: 	n32 = n & (~31);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	accum_05 = _mm512_setzero_pd();
    // OpenBLAS✔️❌: 	accum_15 = _mm512_setzero_pd();
    // OpenBLAS✔️❌: 	accum_25 = _mm512_setzero_pd();
    // OpenBLAS✔️❌: 	accum_35 = _mm512_setzero_pd();
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (; i < n32; i += 32) {
    // OpenBLAS✔️❌: 		accum_05 += _mm512_loadu_pd(&x[i+ 0]) * _mm512_loadu_pd(&y[i+ 0]);
    // OpenBLAS✔️❌: 		accum_15 += _mm512_loadu_pd(&x[i+ 8]) * _mm512_loadu_pd(&y[i+ 8]);
    // OpenBLAS✔️❌: 		accum_25 += _mm512_loadu_pd(&x[i+16]) * _mm512_loadu_pd(&y[i+16]);
    // OpenBLAS✔️❌: 		accum_35 += _mm512_loadu_pd(&x[i+24]) * _mm512_loadu_pd(&y[i+24]);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/*
    // OpenBLAS✔️❌: 	 * we need to fold our 512 bit wide accumulator vectors into 256 bit wide vectors so that the AVX2 code
    // OpenBLAS✔️❌: 	 * below can continue using the intermediate results in its loop
    // OpenBLAS✔️❌: 	 */
    // OpenBLAS✔️❌: 	accum_0 = _mm512_extractf64x4_pd(accum_05, 0) + _mm512_extractf64x4_pd(accum_05, 1);
    // OpenBLAS✔️❌: 	accum_1 = _mm512_extractf64x4_pd(accum_15, 0) + _mm512_extractf64x4_pd(accum_15, 1);
    // OpenBLAS✔️❌: 	accum_2 = _mm512_extractf64x4_pd(accum_25, 0) + _mm512_extractf64x4_pd(accum_25, 1);
    // OpenBLAS✔️❌: 	accum_3 = _mm512_extractf64x4_pd(accum_35, 0) + _mm512_extractf64x4_pd(accum_35, 1);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌: 	for (; i < n; i += 16) {
    // OpenBLAS✔️❌: 		accum_0 += _mm256_loadu_pd(&x[i+ 0]) * _mm256_loadu_pd(&y[i+ 0]);
    // OpenBLAS✔️❌: 		accum_1 += _mm256_loadu_pd(&x[i+ 4]) * _mm256_loadu_pd(&y[i+ 4]);
    // OpenBLAS✔️❌: 		accum_2 += _mm256_loadu_pd(&x[i+ 8]) * _mm256_loadu_pd(&y[i+ 8]);
    // OpenBLAS✔️❌: 		accum_3 += _mm256_loadu_pd(&x[i+12]) * _mm256_loadu_pd(&y[i+12]);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* we now have the partial sums of the dot product in the 4 accumulation vectors, time to consolidate */
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	accum_0 = accum_0 + accum_1 + accum_2 + accum_3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__m128d half_accum0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* Add upper half to lower half of each of the 256 bit vector to get a 128 bit vector */
    // OpenBLAS✔️❌: 	half_accum0 = _mm_add_pd(_mm256_extractf128_pd(accum_0, 0), _mm256_extractf128_pd(accum_0, 1));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* in 128 bit land there is a hadd operation to do the rest of the element-wise sum in one go */
    // OpenBLAS✔️❌: 	half_accum0 = _mm_hadd_pd(half_accum0, half_accum0);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	*dot = half_accum0[0];
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #else
    // OpenBLAS✔️❌: #include "ddot_microk_haswell-2.c"
    // OpenBLAS✔️❌: #endif

    // Scalar lanes preserve the AVX512/AVX2 dependency graph on every target.
    let mut result = 0.0;
    if incx == 1 && incy == 1 {
        let n16 = n & !15;
        let n32 = n16 & !31;
        let mut wide = [[0.0_f64; 8]; 4];
        for i in (0..n32).step_by(32) {
            for block in 0..4 {
                for lane in 0..8 {
                    let j = i + block * 8 + lane;
                    wide[block][lane] = x[j].mul_add(y[j], wide[block][lane]);
                }
            }
        }
        let mut narrow = [[0.0_f64; 4]; 4];
        for block in 0..4 {
            for lane in 0..4 {
                narrow[block][lane] = wide[block][lane] + wide[block][lane + 4];
            }
        }
        for i in (n32..n16).step_by(16) {
            for block in 0..4 {
                for lane in 0..4 {
                    let j = i + block * 4 + lane;
                    narrow[block][lane] = x[j].mul_add(y[j], narrow[block][lane]);
                }
            }
        }
        if n16 != 0 {
            result = reduce4(core::array::from_fn(|lane| {
                ((narrow[0][lane] + narrow[1][lane]) + narrow[2][lane]) + narrow[3][lane]
            }));
        }
        for i in n16..n {
            result = x[i].mul_add(y[i], result);
        }
    } else {
        let mut second = 0.0;
        let n4 = n & !3;
        for i in (0..n4).step_by(4) {
            result += x[i * incx].mul_add(y[i * incy], x[(i + 2) * incx] * y[(i + 2) * incy]);
            second +=
                x[(i + 1) * incx].mul_add(y[(i + 1) * incy], x[(i + 3) * incx] * y[(i + 3) * incy]);
        }
        for i in n4..n {
            result = x[i * incx].mul_add(y[i * incy], result);
        }
        result += second;
    }
    result
}

// DNRM2 needs only nonnegative extended arithmetic. The explicit significand
// has the x87 64-bit precision; exponent denotes the significand's lowest bit.
// No heap allocation or host floating-point environment changes are required.
#[derive(Clone, Copy, Default)]
struct NormExtended {
    significand: u64,
    exponent: i32,
}

impl NormExtended {
    fn rounded(value: u128, exponent: i32) -> Self {
        if value == 0 {
            return Self::default();
        }
        let bits = 128 - value.leading_zeros();
        if bits <= 64 {
            let shift = 64 - bits;
            return Self {
                significand: (value << shift) as u64,
                exponent: exponent - shift as i32,
            };
        }
        let shift = bits - 64;
        let retained = value >> shift;
        let discarded = value & ((1_u128 << shift) - 1);
        let half = 1_u128 << (shift - 1);
        let rounded =
            retained + u128::from(discarded > half || (discarded == half && retained & 1 != 0));
        if rounded == 1_u128 << 64 {
            Self {
                significand: 1_u64 << 63,
                exponent: exponent + shift as i32 + 1,
            }
        } else {
            Self {
                significand: rounded as u64,
                exponent: exponent + shift as i32,
            }
        }
    }
    fn square(value: f64) -> Self {
        let bits = value.to_bits() & !(1_u64 << 63);
        let field = (bits >> 52) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        let (significand, exponent) = if field == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1_u64 << 52), field - 1075)
        };
        Self::rounded(
            u128::from(significand) * u128::from(significand),
            exponent * 2,
        )
    }
    fn add(self, other: Self) -> Self {
        if self.significand == 0 {
            return other;
        }
        if other.significand == 0 {
            return self;
        }
        let (a, b) = if self.exponent >= other.exponent {
            (self, other)
        } else {
            (other, self)
        };
        let shift = (a.exponent - b.exponent) as u32;
        let bits = u128::from(b.significand) << 3;
        let aligned = if shift >= 128 {
            1
        } else if shift == 0 {
            bits
        } else {
            (bits >> shift) | u128::from(bits & ((1_u128 << shift) - 1) != 0)
        };
        Self::rounded((u128::from(a.significand) << 3) + aligned, a.exponent - 3)
    }
    fn sqrt(self) -> Self {
        if self.significand == 0 {
            return self;
        }
        // The shift makes the integer square root a normalized 64-bit value,
        // while leaving an even power of two outside the square root.
        let shift = if self.exponent & 1 == 0 { 64 } else { 63 };
        let radicand = u128::from(self.significand) << shift;
        let root = radicand.isqrt();
        let remainder = radicand - root * root;
        // (root+1/2)^2 = root^2+root+1/4; the integer radicand cannot tie.
        let rounded = root + u128::from(remainder > root);
        Self::rounded(rounded, (self.exponent - shift as i32) / 2)
    }
    fn to_f64(self) -> f64 {
        if self.significand == 0 {
            return 0.0;
        }
        let highest = self.exponent + 63;
        let shift = if highest >= -1022 {
            11
        } else {
            (-1074 - self.exponent) as u32
        };
        if shift > 64 {
            return 0.0;
        }
        let value = u128::from(self.significand);
        let retained = value >> shift;
        let discarded = value & ((1_u128 << shift) - 1);
        let half = 1_u128 << (shift - 1);
        let rounded =
            retained + u128::from(discarded > half || (discarded == half && retained & 1 != 0));
        if highest < -1022 {
            return f64::from_bits(rounded as u64);
        }
        let mut exponent = highest + 1023;
        let mut mantissa = rounded as u64;
        if mantissa == 1_u64 << 53 {
            exponent += 1;
            mantissa >>= 1;
        }
        if exponent >= 2047 {
            return f64::INFINITY;
        }
        f64::from_bits(((exponent as u64) << 52) | (mantissa & ((1_u64 << 52) - 1)))
    }
}

pub(super) fn nrm2(n: usize, x: &[f64], incx: usize) -> f64 {
    // Source: OpenBLAS interface/nrm2.c; these precede its kernel dispatch.
    // OpenBLAS✔️✔️:   if (n <= 0) return 0.;
    // OpenBLAS✔️✔️:   if (n == 1)
    // OpenBLAS✔️✔️:     return fabs(x[0]);
    // OpenBLAS✔️✔️:   if (incx == 0)
    // OpenBLAS✔️✔️:   return (sqrt((double)n)*fabs(x[0]));
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return x[0].abs();
    }
    if incx == 0 {
        return (n as f64).sqrt() * x[0].abs();
    }
    // The safe slice contract proves every native load, including strided tails.
    // Checking the maximum element once also prevents address arithmetic overflow.
    let last = (n - 1)
        .checked_mul(incx)
        .expect("NRM2 vector index overflow");
    assert!(last < x.len(), "NRM2 input shorter than its strided layout");
    #[cfg(target_arch = "x86_64")]
    {
        let mut control = 0_u16;
        // SAFETY: FNSTCW writes exactly this local u16. It only reads the x87
        // mode, leaves the floating-point stack unchanged, and does not alter CW.
        unsafe {
            core::arch::asm!("fnstcw word ptr [{control}]",control=in(reg)&mut control,options(nostack,preserves_flags));
        }
        // Do not cache: callers may change their thread's rounding/precision.
        // Other modes retain the pinned 64-significand/RN-even portable result.
        if control & 0x0f00 == 0x0300 {
            let stride = incx
                .checked_mul(core::mem::size_of::<f64>())
                .expect("NRM2 byte stride overflow");
            // SAFETY: n>=2, positive stride and the maximum-element check above
            // prove all loads remain within x. x87 is mandatory on x86_64;
            // the local CW check proves the frozen kernel's precision/rounding.
            return unsafe { nrm2_x87(n, x.as_ptr(), stride) };
        }
    }
    nrm2_portable(n, x, incx)
}

#[cfg(target_arch = "x86_64")]
unsafe fn nrm2_x87(n: usize, x: *const f64, stride: usize) -> f64 {
    // Source: kernel/x86_64/nrm2.S (frozen DOUBLE x86_64 Linux kernel).
    // OpenBLAS✔️✔️: 	testq	INCX, INCX
    // OpenBLAS✔️✔️: 	je	.L999
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	salq	$BASE_SHIFT, INCX
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	fldz
    // OpenBLAS✔️✔️: 	fldz
    // OpenBLAS✔️✔️: 	fldz
    // OpenBLAS✔️✔️: 	cmpq	$SIZE, INCX
    // OpenBLAS✔️✔️: 	jne	.L40
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	movq	M,  I
    // OpenBLAS✔️✔️: 	sarq	$3, I
    // OpenBLAS✔️✔️: 	jle	.L20
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L10:
    // OpenBLAS✔️✔️: #ifdef PREFETCH
    // OpenBLAS✔️✔️: 	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
    // OpenBLAS✔️✔️: #endif
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	FLD	0 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	1 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	2 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	3 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	faddp	%st, %st(7)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(5)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(3)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(1)
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	FLD	4 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	5 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	6 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	7 * SIZE(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	addq	$8 * SIZE, X
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	faddp	%st, %st(7)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(5)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(3)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(1)
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	decq	I
    // OpenBLAS✔️✔️: 	jg	.L10
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L20:
    // OpenBLAS✔️✔️: 	andq	$7,  M
    // OpenBLAS✔️✔️: 	jle	.L998
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L21:
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	faddp	%st,%st(1)
    // OpenBLAS✔️✔️: 	addq	$1 * SIZE, X
    // OpenBLAS✔️✔️: 	decq	M
    // OpenBLAS✔️✔️: 	jg	.L21
    // OpenBLAS✔️✔️: 	jmp	.L998
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L40:
    // OpenBLAS✔️✔️: 	movq	M, I
    // OpenBLAS✔️✔️: 	sarq	$3,   I
    // OpenBLAS✔️✔️: 	jle	.L60
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L50:
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	faddp	%st, %st(7)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(5)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(3)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(1)
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	faddp	%st, %st(7)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(5)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(3)
    // OpenBLAS✔️✔️: 	faddp	%st, %st(1)
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	decq	I
    // OpenBLAS✔️✔️: 	jg	.L50
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L60:
    // OpenBLAS✔️✔️: 	andq	$7,  M
    // OpenBLAS✔️✔️: 	jle	.L998
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L61:
    // OpenBLAS✔️✔️: 	FLD	(X)
    // OpenBLAS✔️✔️: 	addq	INCX, X
    // OpenBLAS✔️✔️: 	fmul	%st(0), %st
    // OpenBLAS✔️✔️: 	faddp	%st,%st(1)
    // OpenBLAS✔️✔️: 	decq	M
    // OpenBLAS✔️✔️: 	jg	.L61
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L998:
    // OpenBLAS✔️✔️: 	faddp	%st,%st(2)
    // OpenBLAS✔️✔️: 	faddp	%st,%st(1)
    // OpenBLAS✔️✔️: 	faddp	%st,%st(1)
    // OpenBLAS✔️✔️: 	ALIGN_4
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: .L999:
    // OpenBLAS✔️✔️: 	fsqrt
    // OpenBLAS✔️✔️: #ifndef XDOUBLE
    // OpenBLAS✔️✔️: 	sub	$2 * SIZE, %rsp
    // OpenBLAS✔️✔️: 	FST	(%rsp)
    // OpenBLAS✔️✔️: 	MOVSD	(%rsp), %xmm0
    // OpenBLAS✔️✔️: 	add	$2 * SIZE, %rsp
    // OpenBLAS✔️✔️: #endif
    // OpenBLAS✔️✔️: 	ret
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	EPILOGUE
    // OpenBLAS✔️✔️:

    let mut result = 0.0_f64;
    // SAFETY: the caller has checked all source addresses and the current x87
    // precision/rounding. Rust's x86_64 ABI uses XMM for f64 and requires an
    // empty x87 stack at this boundary. Four accumulators plus four products
    // use at most all eight slots; each faddp and the final fstpl balance them.
    // All eight x87 registers and modified integer operands are clobbered.
    // The asm touches only x's checked input region and the local f64 result;
    // it does not write a control word, MXCSR, global state, or machine stack.
    unsafe {
        core::arch::asm!(
            "fldz",
            "fldz",
            "fldz",
            "fldz",
            "cmpq $8, {stride}",
            "jne 5f",
            "movq {remaining}, {blocks}",
            "sarq $3, {blocks}",
            "jle 3f",
            ".p2align 4",
            "2:",
            "fldl 0({input})",
            "fmul %st(0), %st",
            "fldl 8({input})",
            "fmul %st(0), %st",
            "fldl 16({input})",
            "fmul %st(0), %st",
            "fldl 24({input})",
            "fmul %st(0), %st",
            "faddp %st, %st(7)",
            "faddp %st, %st(5)",
            "faddp %st, %st(3)",
            "faddp %st, %st(1)",
            "fldl 32({input})",
            "fmul %st(0), %st",
            "fldl 40({input})",
            "fmul %st(0), %st",
            "fldl 48({input})",
            "fmul %st(0), %st",
            "fldl 56({input})",
            "fmul %st(0), %st",
            "addq $64, {input}",
            "faddp %st, %st(7)",
            "faddp %st, %st(5)",
            "faddp %st, %st(3)",
            "faddp %st, %st(1)",
            "decq {blocks}",
            "jg 2b",
            "3:",
            "andq $7, {remaining}",
            "jle 9f",
            ".p2align 4",
            "4:",
            "fldl ({input})",
            "fmul %st(0), %st",
            "faddp %st, %st(1)",
            "addq $8, {input}",
            "decq {remaining}",
            "jg 4b",
            "jmp 9f",
            "5:",
            "movq {remaining}, {blocks}",
            "sarq $3, {blocks}",
            "jle 7f",
            ".p2align 4",
            "6:",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "faddp %st, %st(7)",
            "faddp %st, %st(5)",
            "faddp %st, %st(3)",
            "faddp %st, %st(1)",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "faddp %st, %st(7)",
            "faddp %st, %st(5)",
            "faddp %st, %st(3)",
            "faddp %st, %st(1)",
            "decq {blocks}",
            "jg 6b",
            "7:",
            "andq $7, {remaining}",
            "jle 9f",
            ".p2align 4",
            "8:",
            "fldl ({input})",
            "addq {stride}, {input}",
            "fmul %st(0), %st",
            "faddp %st, %st(1)",
            "decq {remaining}",
            "jg 8b",
            "9:",
            "faddp %st, %st(2)",
            "faddp %st, %st(1)",
            "faddp %st, %st(1)",
            "fsqrt",
            "fstpl ({output})",
            input=inout(reg)x=>_,remaining=inout(reg)n=>_,stride=in(reg)stride,
            blocks=out(reg)_,output=in(reg)&mut result,
            out("st(0)") _,out("st(1)") _,out("st(2)") _,out("st(3)") _,
            out("st(4)") _,out("st(5)") _,out("st(6)") _,out("st(7)") _,
            options(att_syntax,nostack)
        );
    }
    result
}

fn nrm2_portable(n: usize, x: &[f64], incx: usize) -> f64 {
    // Source: interface/nrm2.c + kernel/x86_64/nrm2.S.
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: FLOATRET NAME(blasint *N, FLOAT *x, blasint *INCX){
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   BLASLONG n    = *N;
    // OpenBLAS✔️❌:   BLASLONG incx = *INCX;
    // OpenBLAS✔️❌:   FLOATRET ret;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   PRINT_DEBUG_NAME;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (n <= 0) return 0.;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef COMPLEX
    // OpenBLAS✔️❌:   if (n == 1)
    // OpenBLAS✔️❌: #ifdef DOUBLE
    // OpenBLAS✔️❌:     return fabs(x[0]);
    // OpenBLAS✔️❌: #else
    // OpenBLAS✔️❌:     return fabsf(x[0]);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incx == 0)
    // OpenBLAS✔️❌: #ifndef COMPLEX
    // OpenBLAS✔️❌: #ifdef DOUBLE
    // OpenBLAS✔️❌:   return (sqrt((double)n)*fabs(x[0]));
    // OpenBLAS✔️❌: 	emms
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	fldz
    // OpenBLAS✔️❌: 	testq	M, M
    // OpenBLAS✔️❌: 	jle	.L999
    // OpenBLAS✔️❌: 	testq	INCX, INCX
    // OpenBLAS✔️❌: 	je	.L999
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	salq	$BASE_SHIFT, INCX
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	fldz
    // OpenBLAS✔️❌: 	fldz
    // OpenBLAS✔️❌: 	fldz
    // OpenBLAS✔️❌: 	cmpq	$SIZE, INCX
    // OpenBLAS✔️❌: 	jne	.L40
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	movq	M,  I
    // OpenBLAS✔️❌: 	sarq	$3, I
    // OpenBLAS✔️❌: 	jle	.L20
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L10:
    // OpenBLAS✔️❌: #ifdef PREFETCH
    // OpenBLAS✔️❌: 	PREFETCH (PREFETCHSIZE +  0) - PREOFFSET(X)
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLD	0 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	1 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	2 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	3 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	faddp	%st, %st(7)
    // OpenBLAS✔️❌: 	faddp	%st, %st(5)
    // OpenBLAS✔️❌: 	faddp	%st, %st(3)
    // OpenBLAS✔️❌: 	faddp	%st, %st(1)
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLD	4 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	5 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	6 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	7 * SIZE(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	addq	$8 * SIZE, X
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	faddp	%st, %st(7)
    // OpenBLAS✔️❌: 	faddp	%st, %st(5)
    // OpenBLAS✔️❌: 	faddp	%st, %st(3)
    // OpenBLAS✔️❌: 	faddp	%st, %st(1)
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	decq	I
    // OpenBLAS✔️❌: 	jg	.L10
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L20:
    // OpenBLAS✔️❌: 	andq	$7,  M
    // OpenBLAS✔️❌: 	jle	.L998
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L21:
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	faddp	%st,%st(1)
    // OpenBLAS✔️❌: 	addq	$1 * SIZE, X
    // OpenBLAS✔️❌: 	decq	M
    // OpenBLAS✔️❌: 	jg	.L21
    // OpenBLAS✔️❌: 	jmp	.L998
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L40:
    // OpenBLAS✔️❌: 	movq	M, I
    // OpenBLAS✔️❌: 	sarq	$3,   I
    // OpenBLAS✔️❌: 	jle	.L60
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L50:
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	faddp	%st, %st(7)
    // OpenBLAS✔️❌: 	faddp	%st, %st(5)
    // OpenBLAS✔️❌: 	faddp	%st, %st(3)
    // OpenBLAS✔️❌: 	faddp	%st, %st(1)
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	faddp	%st, %st(7)
    // OpenBLAS✔️❌: 	faddp	%st, %st(5)
    // OpenBLAS✔️❌: 	faddp	%st, %st(3)
    // OpenBLAS✔️❌: 	faddp	%st, %st(1)
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	decq	I
    // OpenBLAS✔️❌: 	jg	.L50
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L60:
    // OpenBLAS✔️❌: 	andq	$7,  M
    // OpenBLAS✔️❌: 	jle	.L998
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L61:
    // OpenBLAS✔️❌: 	FLD	(X)
    // OpenBLAS✔️❌: 	addq	INCX, X
    // OpenBLAS✔️❌: 	fmul	%st(0), %st
    // OpenBLAS✔️❌: 	faddp	%st,%st(1)
    // OpenBLAS✔️❌: 	decq	M
    // OpenBLAS✔️❌: 	jg	.L61
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L998:
    // OpenBLAS✔️❌: 	faddp	%st,%st(2)
    // OpenBLAS✔️❌: 	faddp	%st,%st(1)
    // OpenBLAS✔️❌: 	faddp	%st,%st(1)
    // OpenBLAS✔️❌: 	ALIGN_4
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: .L999:
    // OpenBLAS✔️❌: 	fsqrt
    // OpenBLAS✔️❌: #ifndef XDOUBLE
    // OpenBLAS✔️❌: 	sub	$2 * SIZE, %rsp
    // OpenBLAS✔️❌: 	FST	(%rsp)
    // OpenBLAS✔️❌: 	MOVSD	(%rsp), %xmm0
    // OpenBLAS✔️❌: 	add	$2 * SIZE, %rsp
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌: 	ret
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	EPILOGUE
    // OpenBLAS✔️❌:

    // The source's x87 PC=64, RN-even operations are emulated exactly.
    // Integer emulation keeps O(n)/O(1) space, but costs more than native x87.
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return x[0].abs();
    }
    if incx == 0 {
        return (n as f64).sqrt() * x[0].abs();
    }
    let mut accumulators = [NormExtended::default(); 4];
    let mut infinity = false;
    let mut nan = false;
    let end = n & !7;
    for i in 0..n {
        let value = x[i * incx];
        if value.is_nan() {
            nan = true;
            continue;
        }
        if value.is_infinite() {
            infinity = true;
            continue;
        }
        let lane = if i < end { i & 3 } else { 0 };
        accumulators[lane] = accumulators[lane].add(NormExtended::square(value));
    }
    if nan {
        return f64::NAN;
    }
    if infinity {
        return f64::INFINITY;
    }
    accumulators[3]
        .add(accumulators[2].add(accumulators[0]).add(accumulators[1]))
        .sqrt()
        .to_f64()
}

pub(super) fn scal(n: usize, alpha: f64, x: &mut [f64], incx: usize) {
    // Source: kernel/x86_64/dscal.c (public dummy2=1).
    // OpenBLAS✔️✔️: }
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: int CNAME(BLASLONG n, BLASLONG dummy0, BLASLONG dummy1, FLOAT da, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *dummy, BLASLONG dummy2)
    // OpenBLAS✔️✔️: {
    // OpenBLAS✔️✔️:     BLASLONG i = 0, j = 0;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️:     // Resolved issue 4728 when the caller is dscal
    // OpenBLAS✔️✔️:     if (dummy2 == 1 && da == 0.0)
    // OpenBLAS✔️✔️:     {
    // OpenBLAS✔️✔️:         if ( inc_x != 1 )
    // OpenBLAS✔️✔️:         {
    // OpenBLAS✔️✔️:             BLASLONG n1 = n & -8;
    // OpenBLAS✔️✔️:             if ( n1 > 0 )
    // OpenBLAS✔️✔️:             {
    // OpenBLAS✔️✔️:                 dscal_kernel_inc_8(n1, &da, x, inc_x);
    // OpenBLAS✔️✔️:                 i = n1 * inc_x;
    // OpenBLAS✔️✔️:                 j = n1;
    // OpenBLAS✔️✔️:             }
    // OpenBLAS✔️✔️:             while(j < n)
    // OpenBLAS✔️✔️:             {
    // OpenBLAS✔️✔️:                 x[i] *= da;
    // OpenBLAS✔️✔️:                 i += inc_x ;
    // OpenBLAS✔️✔️:                 j++;
    // OpenBLAS✔️✔️:             }
    // OpenBLAS✔️✔️:         }
    // OpenBLAS✔️✔️:         else
    // OpenBLAS✔️✔️:         {
    // OpenBLAS✔️✔️:             BLASLONG n1 = n & -8;
    // OpenBLAS✔️✔️:             if ( n1 > 0)
    // OpenBLAS✔️✔️:                 dscal_kernel_8(n1 , &da , x);
    // OpenBLAS✔️✔️:             for ( i = n1 ; i < n; i++ )
    // OpenBLAS✔️✔️:                 x[i] *= da;
    // OpenBLAS✔️✔️:         }
    // OpenBLAS✔️✔️:     }
    // OpenBLAS✔️✔️:     else
    // OpenBLAS✔️✔️:     {
    // OpenBLAS✔️✔️:         if ( inc_x != 1 )

    // Public DSCAL preserves IEEE multiplication even for alpha == zero.
    for i in 0..n {
        x[i * incx] *= alpha;
    }
}

fn scale_output(n: usize, beta: f64, y: &mut [f64], incy: usize) {
    // OpenBLAS✔️✔️:             {
    // OpenBLAS✔️✔️:                 BLASLONG n1 = n & -2;
    // OpenBLAS✔️✔️:                 while(j < n1)
    // OpenBLAS✔️✔️:                 {
    // OpenBLAS✔️✔️:                     x[i] = 0.0;
    // OpenBLAS✔️✔️:                     x[i+inc_x] = 0.0;
    // OpenBLAS✔️✔️:                     i += 2 * inc_x ;
    // OpenBLAS✔️✔️:                     j += 2;
    // OpenBLAS✔️✔️:                 }
    // OpenBLAS✔️✔️:                 while(j < n)
    // OpenBLAS✔️✔️:                 {
    // OpenBLAS✔️✔️:                     x[i] = 0.0;
    // OpenBLAS✔️✔️:                     i += inc_x ;
    // OpenBLAS✔️✔️:                     j++;
    // OpenBLAS✔️✔️:                 }
    // OpenBLAS✔️✔️:             }
    // OpenBLAS✔️✔️:             else
    // OpenBLAS✔️✔️:             {
    // OpenBLAS✔️✔️:                 BLASLONG n1 = n & -8;
    // OpenBLAS✔️✔️:                 if ( n1 > 0 )
    // OpenBLAS✔️✔️:                 {
    // OpenBLAS✔️✔️:                     dscal_kernel_inc_8(n1, &da, x, inc_x);
    // OpenBLAS✔️✔️:                     i = n1 * inc_x;
    // OpenBLAS✔️✔️:                     j = n1;
    // OpenBLAS✔️✔️:                 }
    // OpenBLAS✔️✔️:                 while(j < n)
    // OpenBLAS✔️✔️:                 {
    // OpenBLAS✔️✔️:                     x[i] *= da;
    // OpenBLAS✔️✔️:                     i += inc_x ;
    // OpenBLAS✔️✔️:                     j++;
    // OpenBLAS✔️✔️:                 }
    // OpenBLAS✔️✔️:             }
    // OpenBLAS✔️✔️:         }
    // OpenBLAS✔️✔️:         else
    // OpenBLAS✔️✔️:         {
    // OpenBLAS✔️✔️:             if ( da == 0.0 )
    // OpenBLAS✔️✔️:             {
    // OpenBLAS✔️✔️:                 BLASLONG n1 = n & -8;
    // OpenBLAS✔️✔️:                 if ( n1 > 0)
    // OpenBLAS✔️✔️:                    dscal_kernel_8_zero(n1, &da, x);
    // OpenBLAS✔️✔️:                 for ( i = n1 ; i < n; i++ )
    // OpenBLAS✔️✔️:                     x[i] = 0.0;
    // OpenBLAS✔️✔️:             }
    // OpenBLAS✔️✔️:             else
    // OpenBLAS✔️✔️:             {
    // OpenBLAS✔️✔️:                 BLASLONG n1 = n & -8;
    // OpenBLAS✔️✔️:                 if ( n1 > 0)
    // OpenBLAS✔️✔️:                     dscal_kernel_8(n1 , &da , x);
    // OpenBLAS✔️✔️:                 for ( i = n1 ; i < n; i++ )
    // OpenBLAS✔️✔️:                     x[i] *= da;
    // OpenBLAS✔️✔️:             }
    // OpenBLAS✔️✔️:         }
    // OpenBLAS✔️✔️:     }
    // OpenBLAS✔️✔️:     return(0);
    // OpenBLAS✔️✔️: }

    // The BLAS2 wrappers call SCAL_K with dummy2=0, unlike public DSCAL.
    if beta == 0.0 {
        for i in 0..n {
            y[i * incy] = 0.0;
        }
    } else if beta != 1.0 {
        scal(n, beta, y, incy);
    }
}

pub(super) fn axpy(n: usize, alpha: f64, x: &[f64], incx: usize, y: &mut [f64], incy: usize) {
    // Source: kernel/x86_64/daxpy.c + daxpy_microk_skylakex-2.c.
    // OpenBLAS✔️✔️: 	BLASLONG ix=0,iy=0;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	if ( n <= 0 )  return(0);
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	if ( (inc_x == 1) && (inc_y == 1) )
    // OpenBLAS✔️✔️: 	{
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		BLASLONG n1 = n & -16;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		if ( n1 )
    // OpenBLAS✔️✔️: 			daxpy_kernel_8(n1, x, y , &da );
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		i = n1;
    // OpenBLAS✔️✔️: 		while(i < n)
    // OpenBLAS✔️✔️: 		{
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 			y[i] += da * x[i] ;
    // OpenBLAS✔️✔️: 			i++ ;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		}
    // OpenBLAS✔️✔️: 		return(0);
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	}
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	BLASLONG n1 = n & -4;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	while(i < n1)
    // OpenBLAS✔️✔️: 	{
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		FLOAT m1      = da * x[ix] ;
    // OpenBLAS✔️✔️: 		FLOAT m2      = da * x[ix+inc_x] ;
    // OpenBLAS✔️✔️: 		FLOAT m3      = da * x[ix+2*inc_x] ;
    // OpenBLAS✔️✔️: 		FLOAT m4      = da * x[ix+3*inc_x] ;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		y[iy]         += m1 ;
    // OpenBLAS✔️✔️: 		y[iy+inc_y]   += m2 ;
    // OpenBLAS✔️✔️: 		y[iy+2*inc_y] += m3 ;
    // OpenBLAS✔️✔️: 		y[iy+3*inc_y] += m4 ;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		ix  += inc_x*4 ;
    // OpenBLAS✔️✔️: 		iy  += inc_y*4 ;
    // OpenBLAS✔️✔️: 		i+=4 ;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	}
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	while(i < n)
    // OpenBLAS✔️✔️: 	{
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 		y[iy] += da * x[ix] ;
    // OpenBLAS✔️✔️: 		ix  += inc_x ;
    // OpenBLAS✔️✔️: 		iy  += inc_y ;
    // OpenBLAS✔️✔️: 		i++ ;
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: 	}
    // OpenBLAS✔️✔️: 	return(0);
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️: }
    // OpenBLAS✔️✔️:
    // OpenBLAS✔️✔️:

    if alpha == 0.0 {
        return;
    }
    axpy_kernel(n, alpha, x, incx, y, incy);
}

fn axpy_kernel(n: usize, alpha: f64, x: &[f64], incx: usize, y: &mut [f64], incy: usize) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: runtime feature detection verifies CPU and OS register support.
        // This shares the exact raw-kernel body, including zero-alpha NaN behavior.
        unsafe {
            axpy_kernel_native(n, alpha, x, incx, y, incy);
        }
        return;
    }
    axpy_kernel_impl(n, alpha, x, incx, y, incy);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn axpy_kernel_native(
    n: usize,
    alpha: f64,
    x: &[f64],
    incx: usize,
    y: &mut [f64],
    incy: usize,
) {
    axpy_kernel_impl(n, alpha, x, incx, y, incy);
}

#[inline(always)]
fn axpy_kernel_impl(n: usize, alpha: f64, x: &[f64], incx: usize, y: &mut [f64], incy: usize) {
    // Internal AXPYU_K calls do not apply the public zero-alpha early return.
    for i in 0..n {
        y[i * incy] = alpha.mul_add(x[i * incx], y[i * incy]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_skylakex_arithmetic_and_extended_norm() {
        // OpenBLAS 446c436e, single-threaded SkylakeX DNRM2: exercise
        // extended overflow/underflow range, the four-lane tail and RN-even.
        let cases: &[(&[f64], u64)] = &[
            (&[3.0, 4.0], 0x4014000000000000),
            (&[1e300, 2e300, 3e300], 0x7e5659371f6da5b0),
            (&[1e-300, 2e-300], 0x01b7f5ac55394353),
            (&[0.0, f64::from_bits(1), f64::from_bits(1)], 1),
        ];
        for &(values, expected) in cases {
            assert_eq!(nrm2(values.len(), values, 1).to_bits(), expected);
        }
        let x: Vec<_> = (0..17)
            .map(|i| {
                if i % 2 == 0 {
                    (i + 1) as f64 / 7.0
                } else {
                    -((i + 1) as f64) / 7.0
                }
            })
            .collect();
        assert_eq!(nrm2(x.len(), &x, 1).to_bits(), 0x40182476961e6b67);
        assert_eq!(nrm2(0, &[], 1), 0.0);
        assert!(nrm2(2, &[f64::INFINITY, f64::NAN], 1).is_nan());
        assert_eq!(nrm2(2, &[f64::MAX, f64::MAX], 1), f64::INFINITY);
        assert_eq!(dot(3, &[1., 2., 3.], 1, &[4., 5., 6.], 1), 32.0);
        let mut y = [1., 2., 3.];
        axpy(3, 2., &[4., 5., 6.], 1, &mut y, 1);
        assert_eq!(y, [9., 12., 15.]);
        scal(3, 2., &mut y, 1);
        assert_eq!(y, [18., 24., 30.]);
        // Shared active-disjoint matrix storage must not clone padding or
        // change the source GEMV arithmetic when output occupies that padding.
        let mut backing = [1.0, 3.0, 7.0, 8.0, 2.0, 4.0];
        gemv_shared(
            false,
            2,
            2,
            1.0,
            MatrixSource::Within(0),
            4,
            MatrixSource::Separate(&[5.0, 6.0]),
            1,
            1.0,
            &mut backing,
            2,
            1,
        );
        assert_eq!(backing, [1.0, 3.0, 24.0, 47.0, 2.0, 4.0]);
        // Internal rank updates use AXPYU_K even when a column multiplier is
        // zero; the public DAXPY zero-alpha shortcut must not hide NaN.
        let mut zero_update = [1.0];
        syr2_lower(1, 1.0, &[0.0], 1, &[f64::INFINITY], 1, &mut zero_update, 1);
        assert!(zero_update[0].is_nan());
    }
}

pub(super) fn symv_lower(
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    x: &[f64],
    incx: usize,
    beta: f64,
    y: &mut [f64],
    incy: usize,
) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: runtime detection includes OS support for AVX register state.
        // Only code generation differs: the exact same source arithmetic body
        // runs with hardware FMA, preserving all original reduction boundaries.
        unsafe {
            symv_lower_native(n, alpha, a, lda, x, incx, beta, y, incy);
        }
        return;
    }
    symv_lower_impl(n, alpha, a, lda, x, incx, beta, y, incy);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn symv_lower_native(
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    x: &[f64],
    incx: usize,
    beta: f64,
    y: &mut [f64],
    incy: usize,
) {
    symv_lower_impl(n, alpha, a, lda, x, incx, beta, y, incy);
}

#[inline(always)]
fn symv_lower_impl(
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    x: &[f64],
    incx: usize,
    beta: f64,
    y: &mut [f64],
    incy: usize,
) {
    // Source: kernel/x86_64/dsymv_L.c + dsymv_L_microk_skylakex-2.c.
    // OpenBLAS✔️❌: int CNAME(BLASLONG m, BLASLONG offset, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌: 	BLASLONG i;
    // OpenBLAS✔️❌: 	BLASLONG ix,iy;
    // OpenBLAS✔️❌: 	BLASLONG jx,jy;
    // OpenBLAS✔️❌: 	BLASLONG j;
    // OpenBLAS✔️❌: 	FLOAT temp1;
    // OpenBLAS✔️❌: 	FLOAT temp2;
    // OpenBLAS✔️❌: 	FLOAT tmp1[4];
    // OpenBLAS✔️❌: 	FLOAT tmp2[4];
    // OpenBLAS✔️❌: 	FLOAT *ap[4];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #if 0
    // OpenBLAS✔️❌: 	if ( m != offset )
    // OpenBLAS✔️❌: 		printf("Symv_L: m=%d offset=%d\n",m,offset);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( (inc_x != 1) || (inc_y != 1) )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		jx = 0;
    // OpenBLAS✔️❌: 		jy = 0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		for (j=0; j<offset; j++)
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			temp1 = alpha * x[jx];
    // OpenBLAS✔️❌: 			temp2 = 0.0;
    // OpenBLAS✔️❌: 			y[jy] += temp1 * a[j*lda+j];
    // OpenBLAS✔️❌: 			iy = jy;
    // OpenBLAS✔️❌: 			ix = jx;
    // OpenBLAS✔️❌: 			for (i=j+1; i<m; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				ix += inc_x;
    // OpenBLAS✔️❌: 				iy += inc_y;
    // OpenBLAS✔️❌: 				y[iy] += temp1 * a[j*lda+i];
    // OpenBLAS✔️❌: 				temp2 += a[j*lda+i] * x[ix];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌: 			y[jy] += alpha * temp2;
    // OpenBLAS✔️❌: 			jx    += inc_x;
    // OpenBLAS✔️❌: 			jy    += inc_y;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		return(0);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	BLASLONG offset1 = (offset/4)*4;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (j=0; j<offset1; j+=4)
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		tmp1[0] = alpha * x[j];
    // OpenBLAS✔️❌: 		tmp1[1] = alpha * x[j+1];
    // OpenBLAS✔️❌: 		tmp1[2] = alpha * x[j+2];
    // OpenBLAS✔️❌: 		tmp1[3] = alpha * x[j+3];
    // OpenBLAS✔️❌: 		tmp2[0] = 0.0;
    // OpenBLAS✔️❌: 		tmp2[1] = 0.0;
    // OpenBLAS✔️❌: 		tmp2[2] = 0.0;
    // OpenBLAS✔️❌: 		tmp2[3] = 0.0;
    // OpenBLAS✔️❌: 		ap[0]   = &a[j*lda];
    // OpenBLAS✔️❌: 		ap[1]   = ap[0] + lda;
    // OpenBLAS✔️❌: 		ap[2]   = ap[1] + lda;
    // OpenBLAS✔️❌: 		ap[3]   = ap[2] + lda;
    // OpenBLAS✔️❌: 		y[j]   += tmp1[0] * ap[0][j];
    // OpenBLAS✔️❌: 		y[j+1] += tmp1[1] * ap[1][j+1];
    // OpenBLAS✔️❌: 		y[j+2] += tmp1[2] * ap[2][j+2];
    // OpenBLAS✔️❌: 		y[j+3] += tmp1[3] * ap[3][j+3];
    // OpenBLAS✔️❌: 		BLASLONG from = j+1;
    // OpenBLAS✔️❌: 		if ( m - from >=12 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			BLASLONG m2 = (m/4)*4;
    // OpenBLAS✔️❌: 			for (i=j+1; i<j+4; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[0] * ap[0][i];
    // OpenBLAS✔️❌: 				tmp2[0] += ap[0][i] * x[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j+2; i<j+4; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[1] * ap[1][i];
    // OpenBLAS✔️❌: 				tmp2[1] += ap[1][i] * x[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j+3; i<j+4; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[2] * ap[2][i];
    // OpenBLAS✔️❌: 				tmp2[2] += ap[2][i] * x[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			if ( m2 > j+4 )
    // OpenBLAS✔️❌: 				dsymv_kernel_4x4(j+4,m2,ap,x,y,tmp1,tmp2);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=m2; i<m; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[0] * ap[0][i];
    // OpenBLAS✔️❌: 				tmp2[0] += ap[0][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y[i] += tmp1[1] * ap[1][i];
    // OpenBLAS✔️❌: 				tmp2[1] += ap[1][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y[i] += tmp1[2] * ap[2][i];
    // OpenBLAS✔️❌: 				tmp2[2] += ap[2][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y[i] += tmp1[3] * ap[3][i];
    // OpenBLAS✔️❌: 				tmp2[3] += ap[3][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j+1; i<j+4; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[0] * ap[0][i];
    // OpenBLAS✔️❌: 				tmp2[0] += ap[0][i] * x[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j+2; i<j+4; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[1] * ap[1][i];
    // OpenBLAS✔️❌: 				tmp2[1] += ap[1][i] * x[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j+3; i<j+4; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[2] * ap[2][i];
    // OpenBLAS✔️❌: 				tmp2[2] += ap[2][i] * x[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j+4; i<m; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += tmp1[0] * ap[0][i];
    // OpenBLAS✔️❌: 				tmp2[0] += ap[0][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y[i] += tmp1[1] * ap[1][i];
    // OpenBLAS✔️❌: 				tmp2[1] += ap[1][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y[i] += tmp1[2] * ap[2][i];
    // OpenBLAS✔️❌: 				tmp2[2] += ap[2][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y[i] += tmp1[3] * ap[3][i];
    // OpenBLAS✔️❌: 				tmp2[3] += ap[3][i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		y[j]   += alpha * tmp2[0];
    // OpenBLAS✔️❌: 		y[j+1] += alpha * tmp2[1];
    // OpenBLAS✔️❌: 		y[j+2] += alpha * tmp2[2];
    // OpenBLAS✔️❌: 		y[j+3] += alpha * tmp2[3];
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (j=offset1; j<offset; j++)
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		temp1 = alpha * x[j];
    // OpenBLAS✔️❌: 		temp2 = 0.0;
    // OpenBLAS✔️❌: 		y[j] += temp1 * a[j*lda+j];
    // OpenBLAS✔️❌: 		BLASLONG from = j+1;
    // OpenBLAS✔️❌: 		if ( m - from >=8 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			BLASLONG j1 = ((from + 4)/4)*4;
    // OpenBLAS✔️❌: 			BLASLONG j2 = (m/4)*4;
    // OpenBLAS✔️❌: 			for (i=from; i<j1; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += temp1 * a[j*lda+i];
    // OpenBLAS✔️❌: 				temp2 += a[j*lda+i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j1; i<j2; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += temp1 * a[j*lda+i];
    // OpenBLAS✔️❌: 				temp2 += a[j*lda+i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for (i=j2; i<m; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += temp1 * a[j*lda+i];
    // OpenBLAS✔️❌: 				temp2 += a[j*lda+i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			for (i=from; i<m; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y[i] += temp1 * a[j*lda+i];
    // OpenBLAS✔️❌: 				temp2 += a[j*lda+i] * x[i];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		y[j] += alpha * temp2;
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌: 	return(0);
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__m256d accum_0, accum_1, accum_2, accum_3;
    // OpenBLAS✔️❌: 	__m256d temp1_0, temp1_1, temp1_2, temp1_3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* the 256 bit wide acculmulator vectors start out as zero */
    // OpenBLAS✔️❌: 	accum_0 = _mm256_setzero_pd();
    // OpenBLAS✔️❌: 	accum_1 = _mm256_setzero_pd();
    // OpenBLAS✔️❌: 	accum_2 = _mm256_setzero_pd();
    // OpenBLAS✔️❌: 	accum_3 = _mm256_setzero_pd();
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	temp1_0 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[0]));
    // OpenBLAS✔️❌: 	temp1_1 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[1]));
    // OpenBLAS✔️❌: 	temp1_2 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[2]));
    // OpenBLAS✔️❌: 	temp1_3 = _mm256_broadcastsd_pd(_mm_load_sd(&temp1[3]));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef __AVX512CD__
    // OpenBLAS✔️❌: 	__m512d accum_05, accum_15, accum_25, accum_35;
    // OpenBLAS✔️❌: 	__m512d temp1_05, temp1_15, temp1_25, temp1_35;
    // OpenBLAS✔️❌: 	BLASLONG to2;
    // OpenBLAS✔️❌: 	int delta;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* the 512 bit wide accumulator vectors start out as zero */
    // OpenBLAS✔️❌: 	accum_05 = _mm512_setzero_pd();
    // OpenBLAS✔️❌: 	accum_15 = _mm512_setzero_pd();
    // OpenBLAS✔️❌: 	accum_25 = _mm512_setzero_pd();
    // OpenBLAS✔️❌: 	accum_35 = _mm512_setzero_pd();
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	temp1_05 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[0]));
    // OpenBLAS✔️❌: 	temp1_15 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[1]));
    // OpenBLAS✔️❌: 	temp1_25 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[2]));
    // OpenBLAS✔️❌: 	temp1_35 = _mm512_broadcastsd_pd(_mm_load_sd(&temp1[3]));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	delta = (to - from) & ~7;
    // OpenBLAS✔️❌: 	to2 = from + delta;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (; from < to2; from += 8) {
    // OpenBLAS✔️❌: 		__m512d _x, _y;
    // OpenBLAS✔️❌: 		__m512d a0, a1, a2, a3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		_y = _mm512_loadu_pd(&y[from]);
    // OpenBLAS✔️❌: 		_x = _mm512_loadu_pd(&x[from]);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		a0 = _mm512_loadu_pd(&a[0][from]);
    // OpenBLAS✔️❌: 		a1 = _mm512_loadu_pd(&a[1][from]);
    // OpenBLAS✔️❌: 		a2 = _mm512_loadu_pd(&a[2][from]);
    // OpenBLAS✔️❌: 		a3 = _mm512_loadu_pd(&a[3][from]);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		_y += temp1_05 * a0 + temp1_15 * a1 + temp1_25 * a2 + temp1_35 * a3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		accum_05 += _x * a0;
    // OpenBLAS✔️❌: 		accum_15 += _x * a1;
    // OpenBLAS✔️❌: 		accum_25 += _x * a2;
    // OpenBLAS✔️❌: 		accum_35 += _x * a3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		_mm512_storeu_pd(&y[from], _y);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	};
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/*
    // OpenBLAS✔️❌: 	 * we need to fold our 512 bit wide accumulator vectors into 256 bit wide vectors so that the AVX2 code
    // OpenBLAS✔️❌: 	 * below can continue using the intermediate results in its loop
    // OpenBLAS✔️❌: 	 */
    // OpenBLAS✔️❌: 	accum_0 = _mm256_add_pd(_mm512_extractf64x4_pd(accum_05, 0), _mm512_extractf64x4_pd(accum_05, 1));
    // OpenBLAS✔️❌: 	accum_1 = _mm256_add_pd(_mm512_extractf64x4_pd(accum_15, 0), _mm512_extractf64x4_pd(accum_15, 1));
    // OpenBLAS✔️❌: 	accum_2 = _mm256_add_pd(_mm512_extractf64x4_pd(accum_25, 0), _mm512_extractf64x4_pd(accum_25, 1));
    // OpenBLAS✔️❌: 	accum_3 = _mm256_add_pd(_mm512_extractf64x4_pd(accum_35, 0), _mm512_extractf64x4_pd(accum_35, 1));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (; from != to; from += 4) {
    // OpenBLAS✔️❌: 		__m256d _x, _y;
    // OpenBLAS✔️❌: 		__m256d a0, a1, a2, a3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		_y = _mm256_loadu_pd(&y[from]);
    // OpenBLAS✔️❌: 		_x = _mm256_loadu_pd(&x[from]);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		/* load 4 rows of matrix data */
    // OpenBLAS✔️❌: 		a0 = _mm256_loadu_pd(&a[0][from]);
    // OpenBLAS✔️❌: 		a1 = _mm256_loadu_pd(&a[1][from]);
    // OpenBLAS✔️❌: 		a2 = _mm256_loadu_pd(&a[2][from]);
    // OpenBLAS✔️❌: 		a3 = _mm256_loadu_pd(&a[3][from]);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		_y += temp1_0 * a0 + temp1_1 * a1 + temp1_2 * a2 + temp1_3 * a3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		accum_0 += _x * a0;
    // OpenBLAS✔️❌: 		accum_1 += _x * a1;
    // OpenBLAS✔️❌: 		accum_2 += _x * a2;
    // OpenBLAS✔️❌: 		accum_3 += _x * a3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		_mm256_storeu_pd(&y[from], _y);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	};
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/*
    // OpenBLAS✔️❌: 	 * we now have 4 accumulator vectors. Each vector needs to be summed up element wise and stored in the temp2
    // OpenBLAS✔️❌: 	 * output array. There is no direct instruction for this in 256 bit space, only in 128 space.
    // OpenBLAS✔️❌: 	 */
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__m128d half_accum0, half_accum1, half_accum2, half_accum3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* Add upper half to lower half of each of the four 256 bit vectors to get to four 128 bit vectors */
    // OpenBLAS✔️❌: 	half_accum0 = _mm_add_pd(_mm256_extractf128_pd(accum_0, 0), _mm256_extractf128_pd(accum_0, 1));
    // OpenBLAS✔️❌: 	half_accum1 = _mm_add_pd(_mm256_extractf128_pd(accum_1, 0), _mm256_extractf128_pd(accum_1, 1));
    // OpenBLAS✔️❌: 	half_accum2 = _mm_add_pd(_mm256_extractf128_pd(accum_2, 0), _mm256_extractf128_pd(accum_2, 1));
    // OpenBLAS✔️❌: 	half_accum3 = _mm_add_pd(_mm256_extractf128_pd(accum_3, 0), _mm256_extractf128_pd(accum_3, 1));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* in 128 bit land there is a hadd operation to do the rest of the element-wise sum in one go */
    // OpenBLAS✔️❌: 	half_accum0 = _mm_hadd_pd(half_accum0, half_accum0);
    // OpenBLAS✔️❌: 	half_accum1 = _mm_hadd_pd(half_accum1, half_accum1);
    // OpenBLAS✔️❌: 	half_accum2 = _mm_hadd_pd(half_accum2, half_accum2);
    // OpenBLAS✔️❌: 	half_accum3 = _mm_hadd_pd(half_accum3, half_accum3);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	/* and store the lowest double value from each of these vectors in the temp2 output */
    // OpenBLAS✔️❌: 	temp2[0] += half_accum0[0];
    // OpenBLAS✔️❌: 	temp2[1] += half_accum1[0];
    // OpenBLAS✔️❌: 	temp2[2] += half_accum2[0];

    if n == 0 {
        return;
    }
    scale_output(n, beta, y, incy);
    if alpha == 0.0 {
        return;
    }
    if incx != 1 || incy != 1 {
        for j in 0..n {
            let t = alpha * x[j * incx];
            let mut dot = 0.0;
            y[j * incy] = t.mul_add(a[j * lda + j], y[j * incy]);
            for i in j + 1..n {
                y[i * incy] = t.mul_add(a[j * lda + i], y[i * incy]);
                dot = a[j * lda + i].mul_add(x[i * incx], dot);
            }
            y[j * incy] = alpha.mul_add(dot, y[j * incy]);
        }
        return;
    }
    let n4 = n & !3;
    for j in (0..n4).step_by(4) {
        let t: [f64; 4] = core::array::from_fn(|c| alpha * x[j + c]);
        let mut sums = [0.0_f64; 4];
        for c in 0..4 {
            y[j + c] = t[c].mul_add(a[(j + c) * lda + j + c], y[j + c]);
        }
        for c in 0..3 {
            for i in j + c + 1..j + 4 {
                let v = a[(j + c) * lda + i];
                y[i] = t[c].mul_add(v, y[i]);
                sums[c] = v.mul_add(x[i], sums[c]);
            }
        }
        let mut i = j + 4;
        if n - (j + 1) >= 12 {
            let to = n & !3;
            let to8 = i + ((to - i) & !7);
            let mut wide = [[0.0_f64; 8]; 4];
            while i < to8 {
                for lane in 0..8 {
                    let row = i + lane;
                    let v: [f64; 4] = core::array::from_fn(|c| a[(j + c) * lda + row]);
                    let total =
                        t[3].mul_add(v[3], t[2].mul_add(v[2], t[0].mul_add(v[0], t[1] * v[1])));
                    y[row] += total;
                    for c in 0..4 {
                        wide[c][lane] = v[c].mul_add(x[row], wide[c][lane]);
                    }
                }
                i += 8;
            }
            let mut narrow: [[f64; 4]; 4] =
                core::array::from_fn(|c| core::array::from_fn(|l| wide[c][l] + wide[c][l + 4]));
            while i < to {
                for lane in 0..4 {
                    let row = i + lane;
                    let v: [f64; 4] = core::array::from_fn(|c| a[(j + c) * lda + row]);
                    y[row] +=
                        t[3].mul_add(v[3], t[2].mul_add(v[2], t[0].mul_add(v[0], t[1] * v[1])));
                    for c in 0..4 {
                        narrow[c][lane] = v[c].mul_add(x[row], narrow[c][lane]);
                    }
                }
                i += 4;
            }
            for c in 0..4 {
                sums[c] += reduce4(narrow[c]);
            }
        }
        for row in i..n {
            for c in 0..4 {
                let v = a[(j + c) * lda + row];
                y[row] = t[c].mul_add(v, y[row]);
                sums[c] = v.mul_add(x[row], sums[c]);
            }
        }
        for c in 0..4 {
            y[j + c] = alpha.mul_add(sums[c], y[j + c]);
        }
    }
    for j in n4..n {
        let t = alpha * x[j];
        let mut sum = 0.0;
        y[j] = t.mul_add(a[j * lda + j], y[j]);
        for i in j + 1..n {
            let v = a[j * lda + i];
            y[i] = t.mul_add(v, y[i]);
            sum = v.mul_add(x[i], sum);
        }
        y[j] = alpha.mul_add(sum, y[j]);
    }
}

#[derive(Clone, Copy)]
pub(super) enum MatrixSource<'a> {
    Separate(&'a [f64]),
    Within(usize),
}

struct GemvStorage<'a, 'b> {
    a: MatrixSource<'a>,
    x: MatrixSource<'a>,
    buffer: &'b mut [f64],
    y_offset: usize,
    incy: usize,
}
impl GemvStorage<'_, '_> {
    #[inline]
    fn matrix(&self, index: usize) -> f64 {
        match self.a {
            MatrixSource::Separate(a) => a[index],
            MatrixSource::Within(offset) => self.buffer[offset + index],
        }
    }
    #[inline]
    fn input(&self, index: usize) -> f64 {
        match self.x {
            MatrixSource::Separate(x) => x[index],
            MatrixSource::Within(offset) => self.buffer[offset + index],
        }
    }
    #[inline]
    fn output(&self, index: usize) -> f64 {
        self.buffer[self.y_offset + index * self.incy]
    }
    #[inline]
    fn store(&mut self, index: usize, value: f64) {
        self.buffer[self.y_offset + index * self.incy] = value;
    }
    #[inline]
    fn block_read(&self, row: usize, i: usize, scratch: &[f64]) -> f64 {
        if self.incy == 1 {
            self.output(row + i)
        } else {
            scratch[i]
        }
    }
    #[inline]
    fn block_store(&mut self, row: usize, i: usize, scratch: &mut [f64], value: f64) {
        if self.incy == 1 {
            self.store(row + i, value);
        } else {
            scratch[i] = value;
        }
    }
}

#[cfg(target_arch = "x86_64")]
struct GemvPointers {
    a: *const f64,
    x: *const f64,
    y: *mut f64,
}
#[cfg(target_arch = "x86_64")]
impl GemvStorage<'_, '_> {
    fn pointers(
        &mut self,
        transpose: bool,
        m: usize,
        n: usize,
        lda: usize,
        incx: usize,
    ) -> GemvPointers {
        // All lengths are validated before entering a native microkernel. The
        // addresses remain local to this driver call and cannot outlive buffer.
        let matrix_len = (n - 1)
            .checked_mul(lda)
            .and_then(|v| v.checked_add(m))
            .expect("GEMV matrix extent");
        let x_len = (if transpose { m } else { n } - 1)
            .checked_mul(incx)
            .and_then(|v| v.checked_add(1))
            .expect("GEMV input extent");
        let y_len = (if transpose { n } else { m } - 1)
            .checked_mul(self.incy)
            .and_then(|v| v.checked_add(1))
            .expect("GEMV output extent");
        let buffer_len = self.buffer.len();
        let buffer = self.buffer.as_mut_ptr();
        fn source(
            source: MatrixSource<'_>,
            buffer: *mut f64,
            buffer_len: usize,
            required: usize,
        ) -> *const f64 {
            match source {
                MatrixSource::Separate(slice) => {
                    assert!(required <= slice.len(), "GEMV input storage");
                    slice.as_ptr()
                }
                MatrixSource::Within(offset) => {
                    assert!(
                        offset <= buffer_len && required <= buffer_len - offset,
                        "GEMV input storage"
                    );
                    buffer.wrapping_add(offset)
                }
            }
        }
        let a = source(self.a, buffer, buffer_len, matrix_len);
        let x = source(self.x, buffer, buffer_len, x_len);
        assert!(
            self.y_offset <= buffer_len && y_len <= buffer_len - self.y_offset,
            "GEMV output storage"
        );
        // All Within pointers descend from the same mutable buffer pointer.
        // No overlapping shared/mutable slices are materialized, including
        // padded suffixes whose active matrix/output cells are disjoint.
        let y = buffer.wrapping_add(self.y_offset);
        GemvPointers { a, x, y }
    }
}

pub(super) fn gemv(
    transpose: bool,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    x: &[f64],
    incx: usize,
    beta: f64,
    y: &mut [f64],
    incy: usize,
) {
    gemv_shared(
        transpose,
        m,
        n,
        alpha,
        MatrixSource::Separate(a),
        lda,
        MatrixSource::Separate(x),
        incx,
        beta,
        y,
        0,
        incy,
    );
}

pub(super) fn gemv_shared(
    transpose: bool,
    m: usize,
    n: usize,
    alpha: f64,
    a: MatrixSource<'_>,
    lda: usize,
    x: MatrixSource<'_>,
    incx: usize,
    beta: f64,
    buffer: &mut [f64],
    y_offset: usize,
    incy: usize,
) {
    // One arithmetic kernel serves separate and active-disjoint shared storage;
    // offsets replace raw source pointers without cloning the backing matrix.
    if m == 0 || n == 0 {
        return;
    }
    let mut storage = GemvStorage {
        a,
        x,
        buffer,
        y_offset,
        incy,
    };
    for i in 0..if transpose { n } else { m } {
        if beta == 0.0 {
            storage.store(i, 0.0);
        } else if beta != 1.0 {
            storage.store(i, beta * storage.output(i));
        }
    }
    if alpha == 0.0 {
        return;
    }
    if transpose {
        gemv_t(m, n, alpha, lda, incx, &mut storage);
    } else {
        gemv_n(m, n, alpha, lda, incx, &mut storage);
    }
}

// The source microkernels below replace only SIMD arithmetic inside the common
// driver. They retain its 2048-row blocks, source tails and shared-buffer rules.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn gemv_n_wide(n: usize, a: *const f64, lda: usize, x: [f64; 4], alpha: f64, y: *mut f64) {
    // OpenBLAS❗❗: 	n5 = n & ~7;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	x05 = _mm512_broadcastsd_pd(_mm_load_sd(&x[0]));
    // OpenBLAS❗❗: 	x15 = _mm512_broadcastsd_pd(_mm_load_sd(&x[1]));
    // OpenBLAS❗❗: 	x25 = _mm512_broadcastsd_pd(_mm_load_sd(&x[2]));
    // OpenBLAS❗❗: 	x35 = _mm512_broadcastsd_pd(_mm_load_sd(&x[3]));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	__alpha5 = _mm512_broadcastsd_pd(_mm_load_sd(alpha));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	for (; i < n5; i+= 8) {
    // OpenBLAS❗❗: 		__m512d tempY;
    // OpenBLAS❗❗: 		__m512d sum;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 		sum = _mm512_loadu_pd(&ap[0][i]) * x05 +
    // OpenBLAS❗❗: 		      _mm512_loadu_pd(&ap[1][i]) * x15 +
    // OpenBLAS❗❗: 		      _mm512_loadu_pd(&ap[2][i]) * x25 +
    // OpenBLAS❗❗: 		      _mm512_loadu_pd(&ap[3][i]) * x35;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 		tempY = _mm512_loadu_pd(&y[i]);
    // OpenBLAS❗❗: 		tempY += sum *  __alpha5;
    // OpenBLAS❗❗: 		_mm512_storeu_pd(&y[i], tempY);
    // OpenBLAS❗❗: 	}
    use core::arch::x86_64::*;
    // SAFETY: caller validates all four column extents and n output rows,
    // n is a multiple of eight, and AVX512F is runtime guarded. Each lane
    // evaluates the same expression; loads/stores never form Rust references.
    unsafe {
        let x0 = _mm512_set1_pd(x[0]);
        let x1 = _mm512_set1_pd(x[1]);
        let x2 = _mm512_set1_pd(x[2]);
        let x3 = _mm512_set1_pd(x[3]);
        let alpha = _mm512_set1_pd(alpha);
        for i in (0..n).step_by(8) {
            let v0 = _mm512_loadu_pd(a.add(i));
            let v1 = _mm512_loadu_pd(a.add(lda + i));
            let v2 = _mm512_loadu_pd(a.add(2 * lda + i));
            let v3 = _mm512_loadu_pd(a.add(3 * lda + i));
            let sum = _mm512_fmadd_pd(
                x3,
                v3,
                _mm512_fmadd_pd(x2, v2, _mm512_fmadd_pd(x0, v0, _mm512_mul_pd(x1, v1))),
            );
            _mm512_storeu_pd(
                y.add(i),
                _mm512_fmadd_pd(alpha, sum, _mm512_loadu_pd(y.add(i))),
            );
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn gemv_n_micro<const W: usize>(
    n: usize,
    a: *const f64,
    lda: usize,
    x: *const f64,
    incx: usize,
    alpha: f64,
    y: *mut f64,
) {
    // OpenBLAS❗❗: 	for (; i < n; i+= 4) {
    // OpenBLAS❗❗: 		__m256d tempY;
    // OpenBLAS❗❗: 		__m256d sum;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 		sum = _mm256_loadu_pd(&ap[0][i]) * x0 +
    // OpenBLAS❗❗: 		      _mm256_loadu_pd(&ap[1][i]) * x1 +
    // OpenBLAS❗❗: 		      _mm256_loadu_pd(&ap[2][i]) * x2 +
    // OpenBLAS❗❗: 		      _mm256_loadu_pd(&ap[3][i]) * x3;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 		tempY = _mm256_loadu_pd(&y[i]);
    // OpenBLAS❗❗: 		tempY += sum *  __alpha;
    // OpenBLAS❗❗: 		_mm256_storeu_pd(&y[i], tempY);
    // OpenBLAS❗❗: 	}
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: }
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: #define HAVE_KERNEL_4x2
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: static void dgemv_kernel_4x2( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
    // OpenBLAS❗❗: {
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	int i = 0;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	__m256d x0, x1;
    // OpenBLAS❗❗: 	__m256d __alpha;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	x0 = _mm256_broadcastsd_pd(_mm_load_sd(&x[0]));
    // OpenBLAS❗❗: 	x1 = _mm256_broadcastsd_pd(_mm_load_sd(&x[1]));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	__alpha = _mm256_broadcastsd_pd(_mm_load_sd(alpha));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	for (i = 0; i < n; i+= 4) {
    // OpenBLAS❗❗: 		__m256d tempY;
    // OpenBLAS❗❗: 		__m256d sum;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 		sum = _mm256_loadu_pd(&ap[0][i]) * x0 + _mm256_loadu_pd(&ap[1][i]) * x1;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 		tempY = _mm256_loadu_pd(&y[i]);
    // OpenBLAS❗❗: 		tempY +=  sum *  __alpha;
    // OpenBLAS❗❗: 		_mm256_storeu_pd(&y[i], tempY);
    // OpenBLAS❗❗: 	}
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: }
    // The W=1 source remainder is kernel/x86_64/dgemv_n_4.c:
    // OpenBLAS❗❗: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
    // OpenBLAS❗❗: {
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         BLASLONG register i = 0;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         __asm__  __volatile__
    // OpenBLAS❗❗:         (
    // OpenBLAS❗❗:         "movsd          (%2), %%xmm12            \n\t"  // x0
    // OpenBLAS❗❗:         "mulsd          (%5), %%xmm12            \n\t"  // alpha
    // OpenBLAS❗❗:         "shufpd $0,  %%xmm12, %%xmm12            \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	//        ".align 16                               \n\t"
    // OpenBLAS❗❗:         "1:                             \n\t"
    // OpenBLAS❗❗:         "movups       (%4,%0,8), %%xmm8          \n\t"  // 2 * a
    // OpenBLAS❗❗:         "movups     16(%4,%0,8), %%xmm9          \n\t"  // 2 * a
    // OpenBLAS❗❗:         "movups       (%3,%0,8), %%xmm4          \n\t"  // 2 * y
    // OpenBLAS❗❗:         "movups     16(%3,%0,8), %%xmm5          \n\t"  // 2 * y
    // OpenBLAS❗❗: 	"mulpd          %%xmm12, %%xmm8          \n\t"
    // OpenBLAS❗❗: 	"mulpd          %%xmm12, %%xmm9          \n\t"
    // OpenBLAS❗❗:         "addpd          %%xmm8 , %%xmm4          \n\t"
    // OpenBLAS❗❗:         "addpd          %%xmm9 , %%xmm5          \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movups  %%xmm4 ,    (%3,%0,8)           \n\t"    // 2 * y
    // OpenBLAS❗❗: 	"movups  %%xmm5 ,  16(%3,%0,8)           \n\t"    // 2 * y
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "addq           $4 , %0                  \n\t"
    // OpenBLAS❗❗:         "subq           $4 , %1                  \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "jnz            1b               \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         :
    // OpenBLAS❗❗:           "+r" (i),     // 0
    // OpenBLAS❗❗:           "+r" (n)      // 1
    // OpenBLAS❗❗:         :
    // OpenBLAS❗❗:           "r" (x),      // 2
    // OpenBLAS❗❗:           "r" (y),      // 3
    // OpenBLAS❗❗:           "r" (ap),     // 4
    // OpenBLAS❗❗:           "r" (alpha)   // 5
    // OpenBLAS❗❗:         : "cc",
    // OpenBLAS❗❗:           "%xmm4", "%xmm5",
    // OpenBLAS❗❗:           "%xmm6", "%xmm7",
    // OpenBLAS❗❗:           "%xmm8", "%xmm9", "%xmm10", "%xmm11",
    // OpenBLAS❗❗:           "%xmm12", "%xmm13", "%xmm14", "%xmm15",
    // OpenBLAS❗❗:           "memory"
    // OpenBLAS❗❗:         );
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: }
    use core::arch::x86_64::*;
    // SAFETY: driver validates n multiple-of-four rows and W input columns;
    // W is instantiated only as 1,2,4. AVX2/FMA guard is paired with entry.
    unsafe {
        let xx = [
            *x,
            if W > 1 { *x.add(incx) } else { 0. },
            if W > 2 { *x.add(2 * incx) } else { 0. },
            if W > 3 { *x.add(3 * incx) } else { 0. },
        ];
        let mut i = 0;
        if W == 4 && std::is_x86_feature_detected!("avx512f") {
            let bulk = n & !7;
            gemv_n_wide(bulk, a, lda, xx, alpha, y);
            i = bulk;
        }
        let x0 = _mm256_set1_pd(xx[0]);
        let x1 = _mm256_set1_pd(xx[1]);
        let x2 = _mm256_set1_pd(xx[2]);
        let x3 = _mm256_set1_pd(xx[3]);
        let aa = _mm256_set1_pd(alpha);
        let single = _mm256_set1_pd(alpha * xx[0]);
        while i < n {
            let previous = _mm256_loadu_pd(y.add(i));
            let result = if W == 4 {
                let sum = _mm256_fmadd_pd(
                    x3,
                    _mm256_loadu_pd(a.add(3 * lda + i)),
                    _mm256_fmadd_pd(
                        x2,
                        _mm256_loadu_pd(a.add(2 * lda + i)),
                        _mm256_fmadd_pd(
                            x0,
                            _mm256_loadu_pd(a.add(i)),
                            _mm256_mul_pd(x1, _mm256_loadu_pd(a.add(lda + i))),
                        ),
                    ),
                );
                _mm256_fmadd_pd(aa, sum, previous)
            } else if W == 2 {
                let sum = _mm256_fmadd_pd(
                    x0,
                    _mm256_loadu_pd(a.add(i)),
                    _mm256_mul_pd(x1, _mm256_loadu_pd(a.add(lda + i))),
                );
                _mm256_fmadd_pd(aa, sum, previous)
            } else {
                _mm256_add_pd(previous, _mm256_mul_pd(_mm256_loadu_pd(a.add(i)), single))
            };
            _mm256_storeu_pd(y.add(i), result);
            i += 4;
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn gemv_t_micro(
    n: usize,
    a: *const f64,
    lda: usize,
    x: *const f64,
    incx: usize,
) -> [f64; 4] {
    // OpenBLAS❗❗: static void dgemv_kernel_4x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y) __attribute__ ((noinline));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: static void dgemv_kernel_4x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y)
    // OpenBLAS❗❗: {
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	BLASLONG register i = 0;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	__asm__  __volatile__
    // OpenBLAS❗❗: 	(
    // OpenBLAS❗❗: 	"vzeroupper			         \n\t"
    // OpenBLAS❗❗: 	"vxorpd		%%ymm4 , %%ymm4, %%ymm4  \n\t"
    // OpenBLAS❗❗: 	"vxorpd		%%ymm5 , %%ymm5, %%ymm5  \n\t"
    // OpenBLAS❗❗: 	"vxorpd		%%ymm6 , %%ymm6, %%ymm6  \n\t"
    // OpenBLAS❗❗: 	"vxorpd		%%ymm7 , %%ymm7, %%ymm7  \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "testq          $0x04, %1                      \n\t"
    // OpenBLAS❗❗:         "jz             2f                    \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"vmovups	(%2,%0,8), %%ymm12       \n\t"	// 4 * x
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"vfmadd231pd   (%4,%0,8), %%ymm12, %%ymm4      \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%5,%0,8), %%ymm12, %%ymm5      \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%6,%0,8), %%ymm12, %%ymm6      \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%7,%0,8), %%ymm12, %%ymm7      \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "addq		$4 , %0	  	 	      \n\t"
    // OpenBLAS❗❗: 	"subq	        $4 , %1			      \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "2:                                  \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "cmpq           $0, %1                         \n\t"
    // OpenBLAS❗❗:         "je             3f                      \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	// ".align 16				 \n\t"
    // OpenBLAS❗❗: 	"1:				 \n\t"
    // OpenBLAS❗❗: 	// "prefetcht0	 384(%2,%0,8)		 \n\t"
    // OpenBLAS❗❗: 	"vmovups	(%2,%0,8), %%ymm12       \n\t"	// 4 * x
    // OpenBLAS❗❗: 	"vmovups      32(%2,%0,8), %%ymm13       \n\t"	// 4 * x
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	// "prefetcht0	 384(%4,%0,8)		       \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%4,%0,8), %%ymm12, %%ymm4      \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%5,%0,8), %%ymm12, %%ymm5      \n\t"
    // OpenBLAS❗❗: 	// "prefetcht0	 384(%5,%0,8)		       \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%6,%0,8), %%ymm12, %%ymm6      \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd   (%7,%0,8), %%ymm12, %%ymm7      \n\t"
    // OpenBLAS❗❗: 	// "prefetcht0	 384(%6,%0,8)		       \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd 32(%4,%0,8), %%ymm13, %%ymm4      \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd 32(%5,%0,8), %%ymm13, %%ymm5      \n\t"
    // OpenBLAS❗❗:         "addq		$8 , %0	  	 	       \n\t"
    // OpenBLAS❗❗: 	// "prefetcht0	 384(%7,%0,8)		       \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd -32(%6,%0,8), %%ymm13, %%ymm6     \n\t"
    // OpenBLAS❗❗: 	"subq	        $8 , %1			       \n\t"
    // OpenBLAS❗❗: 	"vfmadd231pd -32(%7,%0,8), %%ymm13, %%ymm7     \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"jnz		1b		      \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "3:                                   \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"vextractf128   $1 , %%ymm4, %%xmm12	      \n\t"
    // OpenBLAS❗❗: 	"vextractf128   $1 , %%ymm5, %%xmm13	      \n\t"
    // OpenBLAS❗❗: 	"vextractf128   $1 , %%ymm6, %%xmm14	      \n\t"
    // OpenBLAS❗❗: 	"vextractf128   $1 , %%ymm7, %%xmm15	      \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"vaddpd		%%xmm4, %%xmm12, %%xmm4       \n\t"
    // OpenBLAS❗❗: 	"vaddpd		%%xmm5, %%xmm13, %%xmm5       \n\t"
    // OpenBLAS❗❗: 	"vaddpd		%%xmm6, %%xmm14, %%xmm6       \n\t"
    // OpenBLAS❗❗: 	"vaddpd		%%xmm7, %%xmm15, %%xmm7       \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "vhaddpd        %%xmm4, %%xmm4, %%xmm4  \n\t"
    // OpenBLAS❗❗:         "vhaddpd        %%xmm5, %%xmm5, %%xmm5  \n\t"
    // OpenBLAS❗❗:         "vhaddpd        %%xmm6, %%xmm6, %%xmm6  \n\t"
    // OpenBLAS❗❗:         "vhaddpd        %%xmm7, %%xmm7, %%xmm7  \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "vmovsd         %%xmm4,    (%3)         \n\t"
    // OpenBLAS❗❗:         "vmovsd         %%xmm5,   8(%3)         \n\t"
    // OpenBLAS❗❗:         "vmovsd         %%xmm6,  16(%3)         \n\t"
    // OpenBLAS❗❗:         "vmovsd         %%xmm7,  24(%3)         \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"vzeroupper			 \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	:
    // OpenBLAS❗❗:           "+r" (i),	// 0
    // OpenBLAS❗❗: 	  "+r" (n)  	// 1
    // OpenBLAS❗❗: 	:
    // OpenBLAS❗❗:           "r" (x),      // 2
    // OpenBLAS❗❗:           "r" (y),      // 3
    // OpenBLAS❗❗:           "r" (ap[0]),  // 4
    // OpenBLAS❗❗:           "r" (ap[1]),  // 5
    // OpenBLAS❗❗:           "r" (ap[2]),  // 6
    // OpenBLAS❗❗:           "r" (ap[3])   // 7
    // OpenBLAS❗❗: 	: "cc",
    // OpenBLAS❗❗: 	  "%xmm0", "%xmm1", "%xmm2", "%xmm3",
    // OpenBLAS❗❗: 	  "%xmm4", "%xmm5", "%xmm6", "%xmm7",
    // OpenBLAS❗❗: 	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
    // OpenBLAS❗❗: 	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
    // OpenBLAS❗❗: 	  "memory"
    // OpenBLAS❗❗: 	);
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: }
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:
    use core::arch::x86_64::*;
    // SAFETY: driver validates four complete columns and n source rows; n is
    // a multiple of four. Only active source loads occur, including positive
    // or zero strided x. Four independent YMM accumulators reuse each X load,
    // matching source per-lane order and the (lane0+2)+(lane1+3) fold.
    unsafe {
        let mut s0 = _mm256_setzero_pd();
        let mut s1 = _mm256_setzero_pd();
        let mut s2 = _mm256_setzero_pd();
        let mut s3 = _mm256_setzero_pd();
        for i in (0..n).step_by(4) {
            let xx = if incx == 1 {
                _mm256_loadu_pd(x.add(i))
            } else {
                _mm256_set_pd(
                    *x.add((i + 3) * incx),
                    *x.add((i + 2) * incx),
                    *x.add((i + 1) * incx),
                    *x.add(i * incx),
                )
            };
            s0 = _mm256_fmadd_pd(_mm256_loadu_pd(a.add(i)), xx, s0);
            s1 = _mm256_fmadd_pd(_mm256_loadu_pd(a.add(lda + i)), xx, s1);
            s2 = _mm256_fmadd_pd(_mm256_loadu_pd(a.add(2 * lda + i)), xx, s2);
            s3 = _mm256_fmadd_pd(_mm256_loadu_pd(a.add(3 * lda + i)), xx, s3);
        }
        let fold = |v| {
            let pair = _mm_add_pd(_mm256_castpd256_pd128(v), _mm256_extractf128_pd::<1>(v));
            _mm_cvtsd_f64(_mm_hadd_pd(pair, pair))
        };
        [fold(s0), fold(s1), fold(s2), fold(s3)]
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn gemv_t_tail<const W: usize>(
    n: usize,
    a: *const f64,
    lda: usize,
    x: *const f64,
    incx: usize,
) -> [f64; 2] {
    // OpenBLAS❗❗: static void dgemv_kernel_4x2(BLASLONG n, FLOAT *ap0, FLOAT *ap1, FLOAT *x, FLOAT *y)  __attribute__ ((noinline));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: static void dgemv_kernel_4x2(BLASLONG n, FLOAT *ap0, FLOAT *ap1, FLOAT *x, FLOAT *y)
    // OpenBLAS❗❗: {
    // OpenBLAS❗❗: 	BLASLONG i;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	i=0;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         __asm__  __volatile__
    // OpenBLAS❗❗:  	(
    // OpenBLAS❗❗: 	"xorpd %%xmm10 , %%xmm10		\n\t"
    // OpenBLAS❗❗: 	"xorpd %%xmm11 , %%xmm11		\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"testq	$2 , %1				\n\t"
    // OpenBLAS❗❗: 	"jz	2f			\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movups  (%5,%0,8) , %%xmm14		\n\t" // x
    // OpenBLAS❗❗: 	"movups  (%3,%0,8) , %%xmm12		\n\t" // ap0
    // OpenBLAS❗❗: 	"movups  (%4,%0,8) , %%xmm13		\n\t" // ap1
    // OpenBLAS❗❗: 	"mulpd   %%xmm14   , %%xmm12 		\n\t"
    // OpenBLAS❗❗: 	"mulpd   %%xmm14   , %%xmm13 		\n\t"
    // OpenBLAS❗❗:         "addq           $2 , %0                 \n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS❗❗:         "subq           $2 , %1                 \n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm13   , %%xmm11		\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "2:                           \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"cmpq	$0, %1				\n\t"
    // OpenBLAS❗❗: 	"je	3f			\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	//        ".align 16                              \n\t"
    // OpenBLAS❗❗:         "1:                            \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movups  (%5,%0,8) , %%xmm14		\n\t" // x
    // OpenBLAS❗❗: 	"movups  (%3,%0,8) , %%xmm12		\n\t" // ap0
    // OpenBLAS❗❗: 	"movups  (%4,%0,8) , %%xmm13		\n\t" // ap1
    // OpenBLAS❗❗: 	"mulpd   %%xmm14   , %%xmm12 		\n\t"
    // OpenBLAS❗❗: 	"mulpd   %%xmm14   , %%xmm13 		\n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm13   , %%xmm11		\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movups  16(%5,%0,8) , %%xmm14		\n\t" // x
    // OpenBLAS❗❗: 	"movups  16(%3,%0,8) , %%xmm12		\n\t" // ap0
    // OpenBLAS❗❗: 	"movups  16(%4,%0,8) , %%xmm13		\n\t" // ap1
    // OpenBLAS❗❗: 	"mulpd   %%xmm14   , %%xmm12 		\n\t"
    // OpenBLAS❗❗: 	"mulpd   %%xmm14   , %%xmm13 		\n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm13   , %%xmm11		\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "addq           $4 , %0                 \n\t"
    // OpenBLAS❗❗:         "subq           $4 , %1                 \n\t"
    // OpenBLAS❗❗:         "jnz            1b              \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "3:                             \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"haddpd        %%xmm10, %%xmm10         \n\t"
    // OpenBLAS❗❗: 	"haddpd        %%xmm11, %%xmm11         \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movsd	       %%xmm10, (%2)	        \n\t"
    // OpenBLAS❗❗: 	"movsd	       %%xmm11,8(%2)	        \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         :
    // OpenBLAS❗❗: 	"+r" (i),	 // 0
    // OpenBLAS❗❗: 	"+r" (n)	 // 1
    // OpenBLAS❗❗: 	:
    // OpenBLAS❗❗:         "r" (y),         // 2
    // OpenBLAS❗❗:         "r" (ap0),       // 3
    // OpenBLAS❗❗:         "r" (ap1),       // 4
    // OpenBLAS❗❗:         "r" (x)          // 5
    // OpenBLAS❗❗:         : "cc",
    // OpenBLAS❗❗:        	"%xmm4", "%xmm5", "%xmm10", "%xmm11",
    // OpenBLAS❗❗:        	"%xmm12", "%xmm13", "%xmm14", "%xmm15",
    // OpenBLAS❗❗:        	"memory"
    // OpenBLAS❗❗:        	);
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: }
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y)  __attribute__ ((noinline));
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y)
    // OpenBLAS❗❗: {
    // OpenBLAS❗❗: 	BLASLONG i;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	i=0;
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         __asm__  __volatile__
    // OpenBLAS❗❗:  	(
    // OpenBLAS❗❗: 	"xorpd %%xmm9  , %%xmm9 		\n\t"
    // OpenBLAS❗❗: 	"xorpd %%xmm10 , %%xmm10		\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"testq	$2 , %1				\n\t"
    // OpenBLAS❗❗: 	"jz	2f			\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movups  (%3,%0,8) , %%xmm12		\n\t"
    // OpenBLAS❗❗: 	"movups  (%4,%0,8) , %%xmm11		\n\t"
    // OpenBLAS❗❗: 	"mulpd   %%xmm11   , %%xmm12 		\n\t"
    // OpenBLAS❗❗:         "addq           $2 , %0                 \n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS❗❗:         "subq           $2 , %1                 \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "2:                           \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"cmpq	$0, %1				\n\t"
    // OpenBLAS❗❗: 	"je	3f			\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	//        ".align 16                              \n\t"
    // OpenBLAS❗❗:         "1:                            \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movups    (%3,%0,8) , %%xmm12		\n\t"
    // OpenBLAS❗❗: 	"movups  16(%3,%0,8) , %%xmm14		\n\t"
    // OpenBLAS❗❗: 	"movups    (%4,%0,8) , %%xmm11		\n\t"
    // OpenBLAS❗❗: 	"movups  16(%4,%0,8) , %%xmm13		\n\t"
    // OpenBLAS❗❗: 	"mulpd   %%xmm11   , %%xmm12 		\n\t"
    // OpenBLAS❗❗: 	"mulpd   %%xmm13   , %%xmm14 		\n\t"
    // OpenBLAS❗❗:         "addq           $4 , %0                 \n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS❗❗:         "subq           $4 , %1                 \n\t"
    // OpenBLAS❗❗: 	"addpd   %%xmm14   , %%xmm9 		\n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "jnz            1b              \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         "3:                             \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"addpd	       %%xmm9 , %%xmm10         \n\t"
    // OpenBLAS❗❗: 	"haddpd        %%xmm10, %%xmm10         \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: 	"movsd	       %%xmm10, (%2)	        \n\t"
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:         :
    // OpenBLAS❗❗: 	"+r" (i),	 // 0
    // OpenBLAS❗❗: 	"+r" (n)	 // 1
    // OpenBLAS❗❗: 	:
    // OpenBLAS❗❗:         "r" (y),         // 2
    // OpenBLAS❗❗:         "r" (ap),        // 3
    // OpenBLAS❗❗:         "r" (x)          // 4
    // OpenBLAS❗❗:         : "cc",
    // OpenBLAS❗❗:        	"%xmm9", "%xmm10" ,
    // OpenBLAS❗❗:        	"%xmm11", "%xmm12", "%xmm13", "%xmm14",
    // OpenBLAS❗❗:        	"memory"
    // OpenBLAS❗❗:        	);
    // OpenBLAS❗❗:
    // OpenBLAS❗❗:
    // OpenBLAS❗❗: }
    // OpenBLAS❗❗:
    use core::arch::x86_64::*;
    // SAFETY: the common driver validates W complete columns and n rows; n
    // is a multiple of four, W=1/2. These source SSE2 tails use separate mul/add,
    // two-lane accumulators, and the exact original final horizontal fold.
    unsafe {
        let mut first = _mm_setzero_pd();
        let mut second = _mm_setzero_pd();
        let mut other = _mm_setzero_pd();
        for i in (0..n).step_by(4) {
            let x01 = if incx == 1 {
                _mm_loadu_pd(x.add(i))
            } else {
                _mm_set_pd(*x.add((i + 1) * incx), *x.add(i * incx))
            };
            let x23 = if incx == 1 {
                _mm_loadu_pd(x.add(i + 2))
            } else {
                _mm_set_pd(*x.add((i + 3) * incx), *x.add((i + 2) * incx))
            };
            first = _mm_add_pd(first, _mm_mul_pd(_mm_loadu_pd(a.add(i)), x01));
            if W == 2 {
                first = _mm_add_pd(first, _mm_mul_pd(_mm_loadu_pd(a.add(i + 2)), x23));
                other = _mm_add_pd(other, _mm_mul_pd(_mm_loadu_pd(a.add(lda + i)), x01));
                other = _mm_add_pd(other, _mm_mul_pd(_mm_loadu_pd(a.add(lda + i + 2)), x23));
            } else {
                second = _mm_add_pd(second, _mm_mul_pd(_mm_loadu_pd(a.add(i + 2)), x23));
            }
        }
        if W == 1 {
            first = _mm_add_pd(first, second);
        }
        [
            _mm_cvtsd_f64(_mm_hadd_pd(first, first)),
            _mm_cvtsd_f64(_mm_hadd_pd(other, other)),
        ]
    }
}

fn gemv_n(
    m: usize,
    n: usize,
    alpha: f64,
    lda: usize,
    incx: usize,
    storage: &mut GemvStorage<'_, '_>,
) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: feature detection verifies the CPU and OS register support.
        // The target-feature entry only inlines the same source arithmetic body.
        unsafe {
            gemv_n_native(m, n, alpha, lda, incx, storage);
        }
        return;
    }
    gemv_n_impl::<false>(m, n, alpha, lda, incx, storage);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn gemv_n_native(
    m: usize,
    n: usize,
    alpha: f64,
    lda: usize,
    incx: usize,
    storage: &mut GemvStorage<'_, '_>,
) {
    gemv_n_impl::<true>(m, n, alpha, lda, incx, storage);
}

#[inline(always)]
fn gemv_n_impl<const NATIVE: bool>(
    m: usize,
    n: usize,
    alpha: f64,
    lda: usize,
    incx: usize,
    storage: &mut GemvStorage<'_, '_>,
) {
    // Source: kernel/x86_64/dgemv_n_4.c + dgemv_n_microk_skylakex-4.c.
    // OpenBLAS✔️❌: 	FLOAT *a_ptr;
    // OpenBLAS✔️❌: 	FLOAT *x_ptr;
    // OpenBLAS✔️❌: 	FLOAT *y_ptr;
    // OpenBLAS✔️❌: 	FLOAT *ap[4];
    // OpenBLAS✔️❌: 	BLASLONG n1;
    // OpenBLAS✔️❌: 	BLASLONG m1;
    // OpenBLAS✔️❌: 	BLASLONG m2;
    // OpenBLAS✔️❌: 	BLASLONG m3;
    // OpenBLAS✔️❌: 	BLASLONG n2;
    // OpenBLAS✔️❌: 	BLASLONG lda4 =  lda << 2;
    // OpenBLAS✔️❌: 	FLOAT xbuffer[8],*ybuffer;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         if ( m < 1 ) return(0);
    // OpenBLAS✔️❌:         if ( n < 1 ) return(0);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	ybuffer = buffer;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	n1 = n >> 2 ;
    // OpenBLAS✔️❌: 	n2 = n &  3 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         m3 = m & 3  ;
    // OpenBLAS✔️❌:         m1 = m & -4 ;
    // OpenBLAS✔️❌:         m2 = (m & (NBMAX-1)) - m3 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	y_ptr = y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	BLASLONG NB = NBMAX;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	while ( NB == NBMAX )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		m1 -= NB;
    // OpenBLAS✔️❌: 		if ( m1 < 0)
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			if ( m2 == 0 ) break;
    // OpenBLAS✔️❌: 			NB = m2;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		a_ptr = a;
    // OpenBLAS✔️❌: 		x_ptr = x;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		ap[0] = a_ptr;
    // OpenBLAS✔️❌: 		ap[1] = a_ptr + lda;
    // OpenBLAS✔️❌: 		ap[2] = ap[1] + lda;
    // OpenBLAS✔️❌: 		ap[3] = ap[2] + lda;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( inc_y != 1 )
    // OpenBLAS✔️❌: 			memset(ybuffer,0,NB*8);
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 			ybuffer = y_ptr;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( inc_x == 1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < n1 ; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				dgemv_kernel_4x4(NB,ap,x_ptr,ybuffer,&alpha);
    // OpenBLAS✔️❌: 				ap[0] += lda4;
    // OpenBLAS✔️❌: 				ap[1] += lda4;
    // OpenBLAS✔️❌: 				ap[2] += lda4;
    // OpenBLAS✔️❌: 				ap[3] += lda4;
    // OpenBLAS✔️❌: 				a_ptr += lda4;
    // OpenBLAS✔️❌: 				x_ptr += 4;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			if ( n2 & 2 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				dgemv_kernel_4x2(NB,ap,x_ptr,ybuffer,&alpha);
    // OpenBLAS✔️❌: 				a_ptr += lda*2;
    // OpenBLAS✔️❌: 				x_ptr += 2;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			if ( n2 & 1 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				dgemv_kernel_4x1(NB,a_ptr,x_ptr,ybuffer,&alpha);
    // OpenBLAS✔️❌: 				/* a_ptr += lda;
    // OpenBLAS✔️❌: 				x_ptr += 1; */
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < n1 ; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				xbuffer[0] = x_ptr[0];
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌: 				xbuffer[1] =  x_ptr[0];
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌: 				xbuffer[2] =  x_ptr[0];
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌: 				xbuffer[3] = x_ptr[0];
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌: 				dgemv_kernel_4x4(NB,ap,xbuffer,ybuffer,&alpha);
    // OpenBLAS✔️❌: 				ap[0] += lda4;
    // OpenBLAS✔️❌: 				ap[1] += lda4;
    // OpenBLAS✔️❌: 				ap[2] += lda4;
    // OpenBLAS✔️❌: 				ap[3] += lda4;
    // OpenBLAS✔️❌: 				a_ptr += lda4;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < n2 ; i++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				xbuffer[0] = x_ptr[0];
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌: 				dgemv_kernel_4x1(NB,a_ptr,xbuffer,ybuffer,&alpha);
    // OpenBLAS✔️❌: 				a_ptr += lda;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		a     += NB;
    // OpenBLAS✔️❌: 		if ( inc_y != 1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			add_y(NB,ybuffer,y_ptr,inc_y);
    // OpenBLAS✔️❌: 			y_ptr += NB * inc_y;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 			y_ptr += NB ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	int i = 0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__m256d x0, x1, x2, x3;
    // OpenBLAS✔️❌: 	__m256d __alpha;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	x0 = _mm256_broadcastsd_pd(_mm_load_sd(&x[0]));
    // OpenBLAS✔️❌: 	x1 = _mm256_broadcastsd_pd(_mm_load_sd(&x[1]));
    // OpenBLAS✔️❌: 	x2 = _mm256_broadcastsd_pd(_mm_load_sd(&x[2]));
    // OpenBLAS✔️❌: 	x3 = _mm256_broadcastsd_pd(_mm_load_sd(&x[3]));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__alpha = _mm256_broadcastsd_pd(_mm_load_sd(alpha));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef __AVX512CD__
    // OpenBLAS✔️❌: 	int n5;
    // OpenBLAS✔️❌: 	__m512d x05, x15, x25, x35;
    // OpenBLAS✔️❌: 	__m512d __alpha5;
    // OpenBLAS✔️❌: 	n5 = n & ~7;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	x05 = _mm512_broadcastsd_pd(_mm_load_sd(&x[0]));
    // OpenBLAS✔️❌: 	x15 = _mm512_broadcastsd_pd(_mm_load_sd(&x[1]));
    // OpenBLAS✔️❌: 	x25 = _mm512_broadcastsd_pd(_mm_load_sd(&x[2]));
    // OpenBLAS✔️❌: 	x35 = _mm512_broadcastsd_pd(_mm_load_sd(&x[3]));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__alpha5 = _mm512_broadcastsd_pd(_mm_load_sd(alpha));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (; i < n5; i+= 8) {
    // OpenBLAS✔️❌: 		__m512d tempY;
    // OpenBLAS✔️❌: 		__m512d sum;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		sum = _mm512_loadu_pd(&ap[0][i]) * x05 +
    // OpenBLAS✔️❌: 		      _mm512_loadu_pd(&ap[1][i]) * x15 +
    // OpenBLAS✔️❌: 		      _mm512_loadu_pd(&ap[2][i]) * x25 +
    // OpenBLAS✔️❌: 		      _mm512_loadu_pd(&ap[3][i]) * x35;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		tempY = _mm512_loadu_pd(&y[i]);
    // OpenBLAS✔️❌: 		tempY += sum *  __alpha5;
    // OpenBLAS✔️❌: 		_mm512_storeu_pd(&y[i], tempY);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (; i < n; i+= 4) {
    // OpenBLAS✔️❌: 		__m256d tempY;
    // OpenBLAS✔️❌: 		__m256d sum;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		sum = _mm256_loadu_pd(&ap[0][i]) * x0 +
    // OpenBLAS✔️❌: 		      _mm256_loadu_pd(&ap[1][i]) * x1 +
    // OpenBLAS✔️❌: 		      _mm256_loadu_pd(&ap[2][i]) * x2 +
    // OpenBLAS✔️❌: 		      _mm256_loadu_pd(&ap[3][i]) * x3;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		tempY = _mm256_loadu_pd(&y[i]);
    // OpenBLAS✔️❌: 		tempY += sum *  __alpha;
    // OpenBLAS✔️❌: 		_mm256_storeu_pd(&y[i], tempY);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #define HAVE_KERNEL_4x2
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void dgemv_kernel_4x2( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	int i = 0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__m256d x0, x1;
    // OpenBLAS✔️❌: 	__m256d __alpha;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	x0 = _mm256_broadcastsd_pd(_mm_load_sd(&x[0]));
    // OpenBLAS✔️❌: 	x1 = _mm256_broadcastsd_pd(_mm_load_sd(&x[1]));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__alpha = _mm256_broadcastsd_pd(_mm_load_sd(alpha));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	for (i = 0; i < n; i+= 4) {
    // OpenBLAS✔️❌: 		__m256d tempY;
    // OpenBLAS✔️❌: 		__m256d sum;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		sum = _mm256_loadu_pd(&ap[0][i]) * x0 + _mm256_loadu_pd(&ap[1][i]) * x1;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		tempY = _mm256_loadu_pd(&y[i]);
    // OpenBLAS✔️❌: 		tempY +=  sum *  __alpha;
    // OpenBLAS✔️❌: 		_mm256_storeu_pd(&y[i], tempY);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #else
    // OpenBLAS✔️❌: #include "dgemv_n_microk_haswell-4.c"
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT *alpha) __attribute__ ((noinline));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y, FLOAT *alpha)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         BLASLONG register i = 0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         __asm__  __volatile__
    // OpenBLAS✔️❌:         (
    // OpenBLAS✔️❌:         "movsd          (%2), %%xmm12            \n\t"  // x0
    // OpenBLAS✔️❌:         "mulsd          (%5), %%xmm12            \n\t"  // alpha
    // OpenBLAS✔️❌:         "shufpd $0,  %%xmm12, %%xmm12            \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	//        ".align 16                               \n\t"
    // OpenBLAS✔️❌:         "1:                             \n\t"
    // OpenBLAS✔️❌:         "movups       (%4,%0,8), %%xmm8          \n\t"  // 2 * a
    // OpenBLAS✔️❌:         "movups     16(%4,%0,8), %%xmm9          \n\t"  // 2 * a
    // OpenBLAS✔️❌:         "movups       (%3,%0,8), %%xmm4          \n\t"  // 2 * y
    // OpenBLAS✔️❌:         "movups     16(%3,%0,8), %%xmm5          \n\t"  // 2 * y
    // OpenBLAS✔️❌: 	"mulpd          %%xmm12, %%xmm8          \n\t"
    // OpenBLAS✔️❌: 	"mulpd          %%xmm12, %%xmm9          \n\t"
    // OpenBLAS✔️❌:         "addpd          %%xmm8 , %%xmm4          \n\t"
    // OpenBLAS✔️❌:         "addpd          %%xmm9 , %%xmm5          \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movups  %%xmm4 ,    (%3,%0,8)           \n\t"    // 2 * y
    // OpenBLAS✔️❌: 	"movups  %%xmm5 ,  16(%3,%0,8)           \n\t"    // 2 * y
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "addq           $4 , %0                  \n\t"
    // OpenBLAS✔️❌:         "subq           $4 , %1                  \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "jnz            1b               \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         :
    // OpenBLAS✔️❌:           "+r" (i),     // 0
    // OpenBLAS✔️❌:           "+r" (n)      // 1
    // OpenBLAS✔️❌:         :
    // OpenBLAS✔️❌:           "r" (x),      // 2
    // OpenBLAS✔️❌:           "r" (y),      // 3
    // OpenBLAS✔️❌:           "r" (ap),     // 4
    // OpenBLAS✔️❌:           "r" (alpha)   // 5
    // OpenBLAS✔️❌:         : "cc",
    // OpenBLAS✔️❌:           "%xmm4", "%xmm5",
    // OpenBLAS✔️❌: 	if ( m3 == 3 )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		a_ptr = a;
    // OpenBLAS✔️❌: 		x_ptr = x;
    // OpenBLAS✔️❌: 		FLOAT temp0 = 0.0;
    // OpenBLAS✔️❌: 		FLOAT temp1 = 0.0;
    // OpenBLAS✔️❌: 		FLOAT temp2 = 0.0;
    // OpenBLAS✔️❌: 		if ( lda == 3 && inc_x ==1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < ( n & -4 ); i+=4 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				temp0 += a_ptr[0] * x_ptr[0] + a_ptr[3] * x_ptr[1];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[1] * x_ptr[0] + a_ptr[4] * x_ptr[1];
    // OpenBLAS✔️❌: 				temp2 += a_ptr[2] * x_ptr[0] + a_ptr[5] * x_ptr[1];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				temp0 += a_ptr[6] * x_ptr[2] + a_ptr[9]  * x_ptr[3];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[7] * x_ptr[2] + a_ptr[10] * x_ptr[3];
    // OpenBLAS✔️❌: 				temp2 += a_ptr[8] * x_ptr[2] + a_ptr[11] * x_ptr[3];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				a_ptr += 12;
    // OpenBLAS✔️❌: 				x_ptr += 4;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( ; i < n; i++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp0 += a_ptr[0] * x_ptr[0];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[1] * x_ptr[0];
    // OpenBLAS✔️❌: 				temp2 += a_ptr[2] * x_ptr[0];
    // OpenBLAS✔️❌: 				a_ptr += 3;
    // OpenBLAS✔️❌: 				x_ptr ++;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < n; i++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp0 += a_ptr[0] * x_ptr[0];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[1] * x_ptr[0];
    // OpenBLAS✔️❌: 				temp2 += a_ptr[2] * x_ptr[0];
    // OpenBLAS✔️❌: 				a_ptr += lda;
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		y_ptr[0] += alpha * temp0;
    // OpenBLAS✔️❌: 		y_ptr += inc_y;
    // OpenBLAS✔️❌: 		y_ptr[0] += alpha * temp1;
    // OpenBLAS✔️❌: 		y_ptr += inc_y;
    // OpenBLAS✔️❌: 		y_ptr[0] += alpha * temp2;
    // OpenBLAS✔️❌: 		return(0);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( m3 == 2 )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		a_ptr = a;
    // OpenBLAS✔️❌: 		x_ptr = x;
    // OpenBLAS✔️❌: 		FLOAT temp0 = 0.0;
    // OpenBLAS✔️❌: 		FLOAT temp1 = 0.0;
    // OpenBLAS✔️❌: 		if ( lda == 2 && inc_x ==1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < (n & -4) ; i+=4 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp0 += a_ptr[0] * x_ptr[0] + a_ptr[2] * x_ptr[1];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[1] * x_ptr[0] + a_ptr[3] * x_ptr[1];
    // OpenBLAS✔️❌: 				temp0 += a_ptr[4] * x_ptr[2] + a_ptr[6] * x_ptr[3];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[5] * x_ptr[2] + a_ptr[7] * x_ptr[3];
    // OpenBLAS✔️❌: 				a_ptr += 8;
    // OpenBLAS✔️❌: 				x_ptr += 4;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( ; i < n; i++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp0 += a_ptr[0]   * x_ptr[0];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[1]   * x_ptr[0];
    // OpenBLAS✔️❌: 				a_ptr += 2;
    // OpenBLAS✔️❌: 				x_ptr ++;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < n; i++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp0 += a_ptr[0] * x_ptr[0];
    // OpenBLAS✔️❌: 				temp1 += a_ptr[1] * x_ptr[0];
    // OpenBLAS✔️❌: 				a_ptr += lda;
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		y_ptr[0] += alpha * temp0;
    // OpenBLAS✔️❌: 		y_ptr += inc_y;
    // OpenBLAS✔️❌: 		y_ptr[0] += alpha * temp1;
    // OpenBLAS✔️❌: 		return(0);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( m3 == 1 )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		a_ptr = a;
    // OpenBLAS✔️❌: 		x_ptr = x;
    // OpenBLAS✔️❌: 		FLOAT temp = 0.0;
    // OpenBLAS✔️❌: 		if ( lda == 1 && inc_x ==1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < (n & -4); i+=4 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp += a_ptr[i] * x_ptr[i] + a_ptr[i+1] * x_ptr[i+1] + a_ptr[i+2] * x_ptr[i+2] + a_ptr[i+3] * x_ptr[i+3];
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( ; i < n; i++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp += a_ptr[i] * x_ptr[i];
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for( i = 0; i < n; i++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				temp += a_ptr[0] * x_ptr[0];
    // OpenBLAS✔️❌: 				a_ptr += lda;
    // OpenBLAS✔️❌: 				x_ptr += inc_x;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		y_ptr[0] += alpha * temp;
    // OpenBLAS✔️❌: 		return(0);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	return(0);
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:

    #[cfg(target_arch = "x86_64")]
    let pointers = NATIVE.then(|| storage.pointers(false, m, n, lda, incx));
    let m4 = m & !3;
    let mut workspace = [core::mem::MaybeUninit::<f64>::uninit(); 2048];
    for row in (0..m4).step_by(2048) {
        let count = (m4 - row).min(2048);
        let scratch: &mut [f64] = if storage.incy != 1 {
            for cell in &mut workspace[..count] {
                cell.write(0.0);
            }
            // SAFETY: exactly count cells were initialized above, the workspace
            // is aligned for f64 and this unique slice ends with this block.
            // Unit-stride output never reads or initializes this source buffer.
            unsafe { core::slice::from_raw_parts_mut(workspace.as_mut_ptr().cast::<f64>(), count) }
        } else {
            &mut []
        };
        let n4 = n & !3;
        for col in (0..n4).step_by(4) {
            #[cfg(target_arch = "x86_64")]
            if let Some(p) = &pointers {
                // SAFETY: pointers() validates the active extents. count is a
                // multiple of four; the native entry is guarded by AVX2/FMA.
                unsafe {
                    let y = if storage.incy == 1 {
                        p.y.add(row)
                    } else {
                        scratch.as_mut_ptr()
                    };
                    gemv_n_micro::<4>(
                        count,
                        p.a.add(col * lda + row),
                        lda,
                        p.x.add(col * incx),
                        incx,
                        alpha,
                        y,
                    );
                }
                continue;
            }
            for i in 0..count {
                let r = row + i;
                let total = storage.input((col + 3) * incx).mul_add(
                    storage.matrix((col + 3) * lda + r),
                    storage.input((col + 2) * incx).mul_add(
                        storage.matrix((col + 2) * lda + r),
                        storage.input(col * incx).mul_add(
                            storage.matrix(col * lda + r),
                            storage.input((col + 1) * incx) * storage.matrix((col + 1) * lda + r),
                        ),
                    ),
                );
                let result = alpha.mul_add(total, storage.block_read(row, i, scratch));
                storage.block_store(row, i, scratch, result);
            }
        }
        let mut col = n4;
        if incx == 1 && n & 2 != 0 {
            #[cfg(target_arch = "x86_64")]
            let handled = if let Some(p) = &pointers {
                // SAFETY: same checked extents and four-row block as above.
                unsafe {
                    let y = if storage.incy == 1 {
                        p.y.add(row)
                    } else {
                        scratch.as_mut_ptr()
                    };
                    gemv_n_micro::<2>(
                        count,
                        p.a.add(col * lda + row),
                        lda,
                        p.x.add(col),
                        1,
                        alpha,
                        y,
                    );
                }
                true
            } else {
                false
            };
            #[cfg(not(target_arch = "x86_64"))]
            let handled = false;
            if !handled {
                for i in 0..count {
                    let r = row + i;
                    let total = storage.input(col).mul_add(
                        storage.matrix(col * lda + r),
                        storage.input(col + 1) * storage.matrix((col + 1) * lda + r),
                    );
                    let result = alpha.mul_add(total, storage.block_read(row, i, scratch));
                    storage.block_store(row, i, scratch, result);
                }
            }
            col += 2;
        }
        for col in col..n {
            #[cfg(target_arch = "x86_64")]
            if let Some(p) = &pointers {
                // SAFETY: checked extents; microkernel preserves the source's
                // separate multiply/add for this one-column remainder.
                unsafe {
                    let y = if storage.incy == 1 {
                        p.y.add(row)
                    } else {
                        scratch.as_mut_ptr()
                    };
                    gemv_n_micro::<1>(
                        count,
                        p.a.add(col * lda + row),
                        lda,
                        p.x.add(col * incx),
                        incx,
                        alpha,
                        y,
                    );
                }
                continue;
            }
            let value = alpha * storage.input(col * incx);
            for i in 0..count {
                let result = storage.block_read(row, i, scratch)
                    + storage.matrix(col * lda + row + i) * value;
                storage.block_store(row, i, scratch, result);
            }
        }
        if storage.incy != 1 {
            #[cfg(target_arch = "x86_64")]
            if let Some(p) = &pointers {
                // SAFETY: validated strided output cells; pointers remain
                // derived from one buffer borrow throughout native blocks.
                unsafe {
                    for i in 0..count {
                        let output = p.y.add((row + i) * storage.incy);
                        *output += scratch[i];
                    }
                }
                continue;
            }
            for i in 0..count {
                storage.store(row + i, storage.output(row + i) + scratch[i]);
            }
        }
    }
    let rem = m - m4;
    for row in m4..m {
        let mut sum = 0.0;
        let mut col = 0;
        if incx == 1 && lda == rem {
            for c in (0..n & !3).step_by(4) {
                if rem == 1 {
                    sum += storage.input(c + 3).mul_add(
                        storage.matrix((c + 3) * lda + row),
                        storage.input(c + 2).mul_add(
                            storage.matrix((c + 2) * lda + row),
                            storage.input(c).mul_add(
                                storage.matrix(c * lda + row),
                                storage.input(c + 1) * storage.matrix((c + 1) * lda + row),
                            ),
                        ),
                    );
                } else {
                    sum += storage.input(c).mul_add(
                        storage.matrix(c * lda + row),
                        storage.input(c + 1) * storage.matrix((c + 1) * lda + row),
                    );
                    sum += storage.input(c + 2).mul_add(
                        storage.matrix((c + 2) * lda + row),
                        storage.input(c + 3) * storage.matrix((c + 3) * lda + row),
                    );
                }
                col = c + 4;
            }
        }
        for c in col..n {
            sum = storage
                .input(c * incx)
                .mul_add(storage.matrix(c * lda + row), sum);
        }
        storage.store(row, alpha.mul_add(sum, storage.output(row)));
    }
}

fn gemv_t(
    m: usize,
    n: usize,
    alpha: f64,
    lda: usize,
    incx: usize,
    storage: &mut GemvStorage<'_, '_>,
) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
        // SAFETY: feature detection verifies the CPU and OS register support.
        // The target-feature entry only inlines the same source arithmetic body.
        unsafe {
            gemv_t_native(m, n, alpha, lda, incx, storage);
        }
        return;
    }
    gemv_t_impl::<false>(m, n, alpha, lda, incx, storage);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
unsafe fn gemv_t_native(
    m: usize,
    n: usize,
    alpha: f64,
    lda: usize,
    incx: usize,
    storage: &mut GemvStorage<'_, '_>,
) {
    gemv_t_impl::<true>(m, n, alpha, lda, incx, storage);
}

#[inline(always)]
fn gemv_t_impl<const NATIVE: bool>(
    m: usize,
    n: usize,
    alpha: f64,
    lda: usize,
    incx: usize,
    storage: &mut GemvStorage<'_, '_>,
) {
    // Source: kernel/x86_64/dgemv_t_4.c + dgemv_t_microk_haswell-4.c.
    // OpenBLAS✔️❌: int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT alpha, FLOAT *a, BLASLONG lda, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT *buffer)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌: 	BLASLONG register i;
    // OpenBLAS✔️❌: 	BLASLONG register j;
    // OpenBLAS✔️❌: 	FLOAT *a_ptr;
    // OpenBLAS✔️❌: 	FLOAT *x_ptr;
    // OpenBLAS✔️❌: 	FLOAT *y_ptr;
    // OpenBLAS✔️❌: 	BLASLONG n0;
    // OpenBLAS✔️❌: 	BLASLONG n1;
    // OpenBLAS✔️❌: 	BLASLONG m1;
    // OpenBLAS✔️❌: 	BLASLONG m2;
    // OpenBLAS✔️❌: 	BLASLONG m3;
    // OpenBLAS✔️❌: 	BLASLONG n2;
    // OpenBLAS✔️❌: 	FLOAT ybuffer[4],*xbuffer;
    // OpenBLAS✔️❌: 	FLOAT *ytemp;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         if ( m < 1 ) return(0);
    // OpenBLAS✔️❌:         if ( n < 1 ) return(0);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	xbuffer = buffer;
    // OpenBLAS✔️❌: 	ytemp   = buffer + (m < NBMAX ? m : NBMAX);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	n0 = n / NBMAX;
    // OpenBLAS✔️❌:         n1 = (n % NBMAX)  >> 2 ;
    // OpenBLAS✔️❌:         n2 = n & 3  ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	m3 = m & 3  ;
    // OpenBLAS✔️❌:         m1 = m & -4 ;
    // OpenBLAS✔️❌:         m2 = (m & (NBMAX-1)) - m3 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	BLASLONG NB = NBMAX;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	while ( NB == NBMAX )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		m1 -= NB;
    // OpenBLAS✔️❌: 		if ( m1 < 0)
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			if ( m2 == 0 ) break;
    // OpenBLAS✔️❌: 			NB = m2;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		y_ptr = y;
    // OpenBLAS✔️❌: 		a_ptr = a;
    // OpenBLAS✔️❌: 		x_ptr = x;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( inc_x == 1 )
    // OpenBLAS✔️❌: 			xbuffer = x_ptr;
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 			copy_x(NB,x_ptr,xbuffer,inc_x);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		FLOAT *ap[4];
    // OpenBLAS✔️❌: 		FLOAT *yp;
    // OpenBLAS✔️❌: 		BLASLONG register lda4 = 4 * lda;
    // OpenBLAS✔️❌: 		ap[0] = a_ptr;
    // OpenBLAS✔️❌: 		ap[1] = a_ptr + lda;
    // OpenBLAS✔️❌: 		ap[2] = ap[1] + lda;
    // OpenBLAS✔️❌: 		ap[3] = ap[2] + lda;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( n0 > 0 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			BLASLONG nb1 = NBMAX / 4;
    // OpenBLAS✔️❌: 			for( j=0; j<n0; j++)
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				yp = ytemp;
    // OpenBLAS✔️❌: 				for( i = 0; i < nb1  ; i++)
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌: 					dgemv_kernel_4x4(NB,ap,xbuffer,yp);
    // OpenBLAS✔️❌: 					ap[0] += lda4 ;
    // OpenBLAS✔️❌: 					ap[1] += lda4 ;
    // OpenBLAS✔️❌: 					ap[2] += lda4 ;
    // OpenBLAS✔️❌: 					ap[3] += lda4 ;
    // OpenBLAS✔️❌: 					yp += 4;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌: 				add_y(nb1*4, alpha, ytemp, y_ptr, inc_y );
    // OpenBLAS✔️❌: 				y_ptr += nb1 * inc_y * 4;
    // OpenBLAS✔️❌: 				a_ptr += nb1 * lda4 ;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		yp = ytemp;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		for( i = 0; i < n1 ; i++)
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			dgemv_kernel_4x4(NB,ap,xbuffer,yp);
    // OpenBLAS✔️❌: 			ap[0] += lda4 ;
    // OpenBLAS✔️❌: 			ap[1] += lda4 ;
    // OpenBLAS✔️❌: 			ap[2] += lda4 ;
    // OpenBLAS✔️❌: 			ap[3] += lda4 ;
    // OpenBLAS✔️❌: 			yp += 4;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		if ( n1 > 0 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			add_y(n1*4, alpha, ytemp, y_ptr, inc_y );
    // OpenBLAS✔️❌: 			y_ptr += n1 * inc_y * 4;
    // OpenBLAS✔️❌: 			a_ptr += n1 * lda4 ;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( n2 & 2 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			dgemv_kernel_4x2(NB,ap[0],ap[1],xbuffer,ybuffer);
    // OpenBLAS✔️❌: 			a_ptr  += lda * 2;
    // OpenBLAS✔️❌: 			*y_ptr += ybuffer[0] * alpha;
    // OpenBLAS✔️❌: 			y_ptr  += inc_y;
    // OpenBLAS✔️❌: 			*y_ptr += ybuffer[1] * alpha;
    // OpenBLAS✔️❌: 			y_ptr  += inc_y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( n2 & 1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			dgemv_kernel_4x1(NB,a_ptr,xbuffer,ybuffer);
    // OpenBLAS✔️❌: 			// a_ptr  += lda;
    // OpenBLAS✔️❌: 			*y_ptr += ybuffer[0] * alpha;
    // OpenBLAS✔️❌: 			// y_ptr  += inc_y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		a += NB;
    // OpenBLAS✔️❌: 		x += NB * inc_x;
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( m3 == 0 ) return(0);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	x_ptr = x;
    // OpenBLAS✔️❌: 	a_ptr = a;
    // OpenBLAS✔️❌: 	if ( m3 == 3 )
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void dgemv_kernel_4x4( BLASLONG n, FLOAT **ap, FLOAT *x, FLOAT *y)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	BLASLONG register i = 0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	__asm__  __volatile__
    // OpenBLAS✔️❌: 	(
    // OpenBLAS✔️❌: 	"vzeroupper			         \n\t"
    // OpenBLAS✔️❌: 	"vxorpd		%%ymm4 , %%ymm4, %%ymm4  \n\t"
    // OpenBLAS✔️❌: 	"vxorpd		%%ymm5 , %%ymm5, %%ymm5  \n\t"
    // OpenBLAS✔️❌: 	"vxorpd		%%ymm6 , %%ymm6, %%ymm6  \n\t"
    // OpenBLAS✔️❌: 	"vxorpd		%%ymm7 , %%ymm7, %%ymm7  \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "testq          $0x04, %1                      \n\t"
    // OpenBLAS✔️❌:         "jz             2f                    \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"vmovups	(%2,%0,8), %%ymm12       \n\t"	// 4 * x
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%4,%0,8), %%ymm12, %%ymm4      \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%5,%0,8), %%ymm12, %%ymm5      \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%6,%0,8), %%ymm12, %%ymm6      \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%7,%0,8), %%ymm12, %%ymm7      \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "addq		$4 , %0	  	 	      \n\t"
    // OpenBLAS✔️❌: 	"subq	        $4 , %1			      \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "2:                                  \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "cmpq           $0, %1                         \n\t"
    // OpenBLAS✔️❌:         "je             3f                      \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	// ".align 16				 \n\t"
    // OpenBLAS✔️❌: 	"1:				 \n\t"
    // OpenBLAS✔️❌: 	// "prefetcht0	 384(%2,%0,8)		 \n\t"
    // OpenBLAS✔️❌: 	"vmovups	(%2,%0,8), %%ymm12       \n\t"	// 4 * x
    // OpenBLAS✔️❌: 	"vmovups      32(%2,%0,8), %%ymm13       \n\t"	// 4 * x
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	// "prefetcht0	 384(%4,%0,8)		       \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%4,%0,8), %%ymm12, %%ymm4      \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%5,%0,8), %%ymm12, %%ymm5      \n\t"
    // OpenBLAS✔️❌: 	// "prefetcht0	 384(%5,%0,8)		       \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%6,%0,8), %%ymm12, %%ymm6      \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd   (%7,%0,8), %%ymm12, %%ymm7      \n\t"
    // OpenBLAS✔️❌: 	// "prefetcht0	 384(%6,%0,8)		       \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd 32(%4,%0,8), %%ymm13, %%ymm4      \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd 32(%5,%0,8), %%ymm13, %%ymm5      \n\t"
    // OpenBLAS✔️❌:         "addq		$8 , %0	  	 	       \n\t"
    // OpenBLAS✔️❌: 	// "prefetcht0	 384(%7,%0,8)		       \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd -32(%6,%0,8), %%ymm13, %%ymm6     \n\t"
    // OpenBLAS✔️❌: 	"subq	        $8 , %1			       \n\t"
    // OpenBLAS✔️❌: 	"vfmadd231pd -32(%7,%0,8), %%ymm13, %%ymm7     \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"jnz		1b		      \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "3:                                   \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"vextractf128   $1 , %%ymm4, %%xmm12	      \n\t"
    // OpenBLAS✔️❌: 	"vextractf128   $1 , %%ymm5, %%xmm13	      \n\t"
    // OpenBLAS✔️❌: 	"vextractf128   $1 , %%ymm6, %%xmm14	      \n\t"
    // OpenBLAS✔️❌: 	"vextractf128   $1 , %%ymm7, %%xmm15	      \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"vaddpd		%%xmm4, %%xmm12, %%xmm4       \n\t"
    // OpenBLAS✔️❌: 	"vaddpd		%%xmm5, %%xmm13, %%xmm5       \n\t"
    // OpenBLAS✔️❌: 	"vaddpd		%%xmm6, %%xmm14, %%xmm6       \n\t"
    // OpenBLAS✔️❌: 	"vaddpd		%%xmm7, %%xmm15, %%xmm7       \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "vhaddpd        %%xmm4, %%xmm4, %%xmm4  \n\t"
    // OpenBLAS✔️❌:         "vhaddpd        %%xmm5, %%xmm5, %%xmm5  \n\t"
    // OpenBLAS✔️❌:         "vhaddpd        %%xmm6, %%xmm6, %%xmm6  \n\t"
    // OpenBLAS✔️❌:         "vhaddpd        %%xmm7, %%xmm7, %%xmm7  \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "vmovsd         %%xmm4,    (%3)         \n\t"
    // OpenBLAS✔️❌:         "vmovsd         %%xmm5,   8(%3)         \n\t"
    // OpenBLAS✔️❌:         "vmovsd         %%xmm6,  16(%3)         \n\t"
    // OpenBLAS✔️❌:         "vmovsd         %%xmm7,  24(%3)         \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"vzeroupper			 \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	:
    // OpenBLAS✔️❌:           "+r" (i),	// 0
    // OpenBLAS✔️❌: 	  "+r" (n)  	// 1
    // OpenBLAS✔️❌: 	:
    // OpenBLAS✔️❌:           "r" (x),      // 2
    // OpenBLAS✔️❌:           "r" (y),      // 3
    // OpenBLAS✔️❌:           "r" (ap[0]),  // 4
    // OpenBLAS✔️❌:           "r" (ap[1]),  // 5
    // OpenBLAS✔️❌:           "r" (ap[2]),  // 6
    // OpenBLAS✔️❌:           "r" (ap[3])   // 7
    // OpenBLAS✔️❌: 	: "cc",
    // OpenBLAS✔️❌: 	  "%xmm0", "%xmm1", "%xmm2", "%xmm3",
    // OpenBLAS✔️❌: 	  "%xmm4", "%xmm5", "%xmm6", "%xmm7",
    // OpenBLAS✔️❌: 	  "%xmm8", "%xmm9", "%xmm10", "%xmm11",
    // OpenBLAS✔️❌: 	  "%xmm12", "%xmm13", "%xmm14", "%xmm15",
    // OpenBLAS✔️❌: 	  "memory"
    // OpenBLAS✔️❌: static void dgemv_kernel_4x2(BLASLONG n, FLOAT *ap0, FLOAT *ap1, FLOAT *x, FLOAT *y)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌: 	BLASLONG i;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	i=0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         __asm__  __volatile__
    // OpenBLAS✔️❌:  	(
    // OpenBLAS✔️❌: 	"xorpd %%xmm10 , %%xmm10		\n\t"
    // OpenBLAS✔️❌: 	"xorpd %%xmm11 , %%xmm11		\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"testq	$2 , %1				\n\t"
    // OpenBLAS✔️❌: 	"jz	2f			\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movups  (%5,%0,8) , %%xmm14		\n\t" // x
    // OpenBLAS✔️❌: 	"movups  (%3,%0,8) , %%xmm12		\n\t" // ap0
    // OpenBLAS✔️❌: 	"movups  (%4,%0,8) , %%xmm13		\n\t" // ap1
    // OpenBLAS✔️❌: 	"mulpd   %%xmm14   , %%xmm12 		\n\t"
    // OpenBLAS✔️❌: 	"mulpd   %%xmm14   , %%xmm13 		\n\t"
    // OpenBLAS✔️❌:         "addq           $2 , %0                 \n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS✔️❌:         "subq           $2 , %1                 \n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm13   , %%xmm11		\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "2:                           \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"cmpq	$0, %1				\n\t"
    // OpenBLAS✔️❌: 	"je	3f			\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	//        ".align 16                              \n\t"
    // OpenBLAS✔️❌:         "1:                            \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movups  (%5,%0,8) , %%xmm14		\n\t" // x
    // OpenBLAS✔️❌: 	"movups  (%3,%0,8) , %%xmm12		\n\t" // ap0
    // OpenBLAS✔️❌: 	"movups  (%4,%0,8) , %%xmm13		\n\t" // ap1
    // OpenBLAS✔️❌: 	"mulpd   %%xmm14   , %%xmm12 		\n\t"
    // OpenBLAS✔️❌: 	"mulpd   %%xmm14   , %%xmm13 		\n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm13   , %%xmm11		\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movups  16(%5,%0,8) , %%xmm14		\n\t" // x
    // OpenBLAS✔️❌: 	"movups  16(%3,%0,8) , %%xmm12		\n\t" // ap0
    // OpenBLAS✔️❌: 	"movups  16(%4,%0,8) , %%xmm13		\n\t" // ap1
    // OpenBLAS✔️❌: 	"mulpd   %%xmm14   , %%xmm12 		\n\t"
    // OpenBLAS✔️❌: 	"mulpd   %%xmm14   , %%xmm13 		\n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm13   , %%xmm11		\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "addq           $4 , %0                 \n\t"
    // OpenBLAS✔️❌:         "subq           $4 , %1                 \n\t"
    // OpenBLAS✔️❌:         "jnz            1b              \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "3:                             \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"haddpd        %%xmm10, %%xmm10         \n\t"
    // OpenBLAS✔️❌: 	"haddpd        %%xmm11, %%xmm11         \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movsd	       %%xmm10, (%2)	        \n\t"
    // OpenBLAS✔️❌: 	"movsd	       %%xmm11,8(%2)	        \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         :
    // OpenBLAS✔️❌: 	"+r" (i),	 // 0
    // OpenBLAS✔️❌: 	"+r" (n)	 // 1
    // OpenBLAS✔️❌: 	:
    // OpenBLAS✔️❌:         "r" (y),         // 2
    // OpenBLAS✔️❌:         "r" (ap0),       // 3
    // OpenBLAS✔️❌:         "r" (ap1),       // 4
    // OpenBLAS✔️❌:         "r" (x)          // 5
    // OpenBLAS✔️❌:         : "cc",
    // OpenBLAS✔️❌:        	"%xmm4", "%xmm5", "%xmm10", "%xmm11",
    // OpenBLAS✔️❌:        	"%xmm12", "%xmm13", "%xmm14", "%xmm15",
    // OpenBLAS✔️❌:        	"memory"
    // OpenBLAS✔️❌:        	);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y)  __attribute__ ((noinline));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void dgemv_kernel_4x1(BLASLONG n, FLOAT *ap, FLOAT *x, FLOAT *y)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌: 	BLASLONG i;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	i=0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         __asm__  __volatile__
    // OpenBLAS✔️❌:  	(
    // OpenBLAS✔️❌: 	"xorpd %%xmm9  , %%xmm9 		\n\t"
    // OpenBLAS✔️❌: 	"xorpd %%xmm10 , %%xmm10		\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"testq	$2 , %1				\n\t"
    // OpenBLAS✔️❌: 	"jz	2f			\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movups  (%3,%0,8) , %%xmm12		\n\t"
    // OpenBLAS✔️❌: 	"movups  (%4,%0,8) , %%xmm11		\n\t"
    // OpenBLAS✔️❌: 	"mulpd   %%xmm11   , %%xmm12 		\n\t"
    // OpenBLAS✔️❌:         "addq           $2 , %0                 \n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS✔️❌:         "subq           $2 , %1                 \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "2:                           \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"cmpq	$0, %1				\n\t"
    // OpenBLAS✔️❌: 	"je	3f			\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	//        ".align 16                              \n\t"
    // OpenBLAS✔️❌:         "1:                            \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movups    (%3,%0,8) , %%xmm12		\n\t"
    // OpenBLAS✔️❌: 	"movups  16(%3,%0,8) , %%xmm14		\n\t"
    // OpenBLAS✔️❌: 	"movups    (%4,%0,8) , %%xmm11		\n\t"
    // OpenBLAS✔️❌: 	"movups  16(%4,%0,8) , %%xmm13		\n\t"
    // OpenBLAS✔️❌: 	"mulpd   %%xmm11   , %%xmm12 		\n\t"
    // OpenBLAS✔️❌: 	"mulpd   %%xmm13   , %%xmm14 		\n\t"
    // OpenBLAS✔️❌:         "addq           $4 , %0                 \n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm12   , %%xmm10		\n\t"
    // OpenBLAS✔️❌:         "subq           $4 , %1                 \n\t"
    // OpenBLAS✔️❌: 	"addpd   %%xmm14   , %%xmm9 		\n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "jnz            1b              \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         "3:                             \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"addpd	       %%xmm9 , %%xmm10         \n\t"
    // OpenBLAS✔️❌: 	"haddpd        %%xmm10, %%xmm10         \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	"movsd	       %%xmm10, (%2)	        \n\t"
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         :
    // OpenBLAS✔️❌: 	"+r" (i),	 // 0
    // OpenBLAS✔️❌: 	"+r" (n)	 // 1
    // OpenBLAS✔️❌: 	:
    // OpenBLAS✔️❌:         "r" (y),         // 2
    // OpenBLAS✔️❌:         "r" (ap),        // 3
    // OpenBLAS✔️❌:         "r" (x)          // 4
    // OpenBLAS✔️❌:         : "cc",
    // OpenBLAS✔️❌:        	"%xmm9", "%xmm10" ,
    // OpenBLAS✔️❌:        	"%xmm11", "%xmm12", "%xmm13", "%xmm14",
    // OpenBLAS✔️❌:        	"memory"
    // OpenBLAS✔️❌:        	);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void copy_x(BLASLONG n, FLOAT *src, FLOAT *dest, BLASLONG inc_src)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:         BLASLONG i;
    // OpenBLAS✔️❌:         for ( i=0; i<n; i++ )
    // OpenBLAS✔️❌:         {
    // OpenBLAS✔️❌:                 *dest = *src;
    // OpenBLAS✔️❌:                 dest++;
    // OpenBLAS✔️❌:                 src += inc_src;
    // OpenBLAS✔️❌:         }
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void add_y(BLASLONG n, FLOAT da , FLOAT *src, FLOAT *dest, BLASLONG inc_dest) __attribute__ ((noinline));
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: static void add_y(BLASLONG n, FLOAT da , FLOAT *src, FLOAT *dest, BLASLONG inc_dest)
    // OpenBLAS✔️❌: {
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         BLASLONG i;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( inc_dest != 1 )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:         	for ( i=0; i<n; i++ )
    // OpenBLAS✔️❌:         	{
    // OpenBLAS✔️❌:                 	*dest += src[i]  * da;
    // OpenBLAS✔️❌:                 	dest  += inc_dest;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		return;
    // OpenBLAS✔️❌:         }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	i=0;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:         __asm__  __volatile__
    // OpenBLAS✔️❌:  	(
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		FLOAT xtemp0 = *x_ptr * alpha;
    // OpenBLAS✔️❌: 		x_ptr += inc_x;
    // OpenBLAS✔️❌: 		FLOAT xtemp1 = *x_ptr * alpha;
    // OpenBLAS✔️❌: 		x_ptr += inc_x;
    // OpenBLAS✔️❌: 		FLOAT xtemp2 = *x_ptr * alpha;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		FLOAT *aj = a_ptr;
    // OpenBLAS✔️❌: 		y_ptr = y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( lda == 3 && inc_y == 1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for ( j=0; j< ( n & -4) ; j+=4 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				y_ptr[j]   += aj[0] * xtemp0 + aj[1]  * xtemp1 + aj[2]  * xtemp2;
    // OpenBLAS✔️❌: 				y_ptr[j+1] += aj[3] * xtemp0 + aj[4]  * xtemp1 + aj[5]  * xtemp2;
    // OpenBLAS✔️❌: 				y_ptr[j+2] += aj[6] * xtemp0 + aj[7]  * xtemp1 + aj[8]  * xtemp2;
    // OpenBLAS✔️❌: 				y_ptr[j+3] += aj[9] * xtemp0 + aj[10] * xtemp1 + aj[11] * xtemp2;
    // OpenBLAS✔️❌: 			 	aj        += 12;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for ( ; j<n; j++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y_ptr[j]  += aj[0] * xtemp0 + aj[1] * xtemp1 + aj[2] * xtemp2;
    // OpenBLAS✔️❌: 			 	aj        += 3;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			if ( inc_y == 1 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				BLASLONG register lda2 = lda << 1;
    // OpenBLAS✔️❌: 				BLASLONG register lda4 = lda << 2;
    // OpenBLAS✔️❌: 				BLASLONG register lda3 = lda2 + lda;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				for ( j=0; j< ( n & -4 ); j+=4 )
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 					y_ptr[j]    += *aj        * xtemp0 + *(aj+1)      * xtemp1 + *(aj+2)      * xtemp2;
    // OpenBLAS✔️❌: 					y_ptr[j+1]  += *(aj+lda)  * xtemp0 + *(aj+lda+1)  * xtemp1 + *(aj+lda+2)  * xtemp2;
    // OpenBLAS✔️❌: 					y_ptr[j+2]  += *(aj+lda2) * xtemp0 + *(aj+lda2+1) * xtemp1 + *(aj+lda2+2) * xtemp2;
    // OpenBLAS✔️❌: 					y_ptr[j+3]  += *(aj+lda3) * xtemp0 + *(aj+lda3+1) * xtemp1 + *(aj+lda3+2) * xtemp2;
    // OpenBLAS✔️❌: 			 		aj          += lda4;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				for ( ; j< n ; j++ )
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 					y_ptr[j]    += *aj * xtemp0 + *(aj+1) * xtemp1 + *(aj+2) * xtemp2 ;
    // OpenBLAS✔️❌: 			 		aj          += lda;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌: 			else
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				for ( j=0; j<n; j++ )
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌: 					*y_ptr += *aj * xtemp0 + *(aj+1) * xtemp1 + *(aj+2) * xtemp2;
    // OpenBLAS✔️❌: 				 	y_ptr += inc_y;
    // OpenBLAS✔️❌: 			 		aj    += lda;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		return(0);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if ( m3 == 2 )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		FLOAT xtemp0 = *x_ptr * alpha;
    // OpenBLAS✔️❌: 		x_ptr += inc_x;
    // OpenBLAS✔️❌: 		FLOAT xtemp1 = *x_ptr * alpha;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		FLOAT *aj = a_ptr;
    // OpenBLAS✔️❌: 		y_ptr = y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( lda == 2 && inc_y == 1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for ( j=0; j< ( n & -4) ; j+=4 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y_ptr[j]   += aj[0] * xtemp0 + aj[1] * xtemp1 ;
    // OpenBLAS✔️❌: 				y_ptr[j+1] += aj[2] * xtemp0 + aj[3] * xtemp1 ;
    // OpenBLAS✔️❌: 				y_ptr[j+2] += aj[4] * xtemp0 + aj[5] * xtemp1 ;
    // OpenBLAS✔️❌: 				y_ptr[j+3] += aj[6] * xtemp0 + aj[7] * xtemp1 ;
    // OpenBLAS✔️❌: 			 	aj         += 8;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for ( ; j<n; j++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y_ptr[j] += aj[0] * xtemp0 + aj[1] * xtemp1 ;
    // OpenBLAS✔️❌: 			 	aj       += 2;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			if ( inc_y == 1 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				BLASLONG register lda2 = lda << 1;
    // OpenBLAS✔️❌: 				BLASLONG register lda4 = lda << 2;
    // OpenBLAS✔️❌: 				BLASLONG register lda3 = lda2 + lda;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				for ( j=0; j< ( n & -4 ); j+=4 )
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 					y_ptr[j]    += *aj        * xtemp0 + *(aj+1)      * xtemp1 ;
    // OpenBLAS✔️❌: 					y_ptr[j+1]  += *(aj+lda)  * xtemp0 + *(aj+lda+1)  * xtemp1 ;
    // OpenBLAS✔️❌: 					y_ptr[j+2]  += *(aj+lda2) * xtemp0 + *(aj+lda2+1) * xtemp1 ;
    // OpenBLAS✔️❌: 					y_ptr[j+3]  += *(aj+lda3) * xtemp0 + *(aj+lda3+1) * xtemp1 ;
    // OpenBLAS✔️❌: 			 		aj          += lda4;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 				for ( ; j< n ; j++ )
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 					y_ptr[j]    += *aj * xtemp0 + *(aj+1) * xtemp1 ;
    // OpenBLAS✔️❌: 			 		aj          += lda;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌: 			else
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				for ( j=0; j<n; j++ )
    // OpenBLAS✔️❌: 				{
    // OpenBLAS✔️❌: 					*y_ptr += *aj * xtemp0 + *(aj+1) * xtemp1 ;
    // OpenBLAS✔️❌: 			 		y_ptr += inc_y;
    // OpenBLAS✔️❌: 			 		aj    += lda;
    // OpenBLAS✔️❌: 				}
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		return(0);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLOAT xtemp = *x_ptr * alpha;
    // OpenBLAS✔️❌: 	FLOAT *aj = a_ptr;
    // OpenBLAS✔️❌: 	y_ptr = y;
    // OpenBLAS✔️❌: 	if ( lda == 1 && inc_y == 1 )
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		for ( j=0; j< ( n & -4) ; j+=4 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			y_ptr[j]   += aj[j]   * xtemp;
    // OpenBLAS✔️❌: 			y_ptr[j+1] += aj[j+1] * xtemp;
    // OpenBLAS✔️❌: 			y_ptr[j+2] += aj[j+2] * xtemp;
    // OpenBLAS✔️❌: 			y_ptr[j+3] += aj[j+3] * xtemp;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		for ( ; j<n   ; j++ )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			y_ptr[j] += aj[j] * xtemp;
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌: 	else
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌: 		if ( inc_y == 1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			BLASLONG register lda2 = lda << 1;
    // OpenBLAS✔️❌: 			BLASLONG register lda4 = lda << 2;
    // OpenBLAS✔️❌: 			BLASLONG register lda3 = lda2 + lda;
    // OpenBLAS✔️❌: 			for ( j=0; j< ( n & -4 ); j+=4 )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y_ptr[j]    += *aj        * xtemp;
    // OpenBLAS✔️❌: 				y_ptr[j+1]  += *(aj+lda)  * xtemp;
    // OpenBLAS✔️❌: 				y_ptr[j+2]  += *(aj+lda2) * xtemp;
    // OpenBLAS✔️❌: 				y_ptr[j+3]  += *(aj+lda3) * xtemp;
    // OpenBLAS✔️❌: 		 		aj          += lda4  ;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 			for ( ; j<n; j++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				y_ptr[j]  += *aj * xtemp;
    // OpenBLAS✔️❌: 		 		aj        += lda;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 		else
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			for ( j=0; j<n; j++ )
    // OpenBLAS✔️❌: 			{
    // OpenBLAS✔️❌: 				*y_ptr += *aj * xtemp;
    // OpenBLAS✔️❌: 		 		y_ptr += inc_y;
    // OpenBLAS✔️❌: 		 		aj    += lda;
    // OpenBLAS✔️❌: 			}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:

    let m4 = m & !3;
    for row in (0..m4).step_by(2048) {
        #[cfg(target_arch = "x86_64")]
        let pointers = NATIVE.then(|| storage.pointers(true, m, n, lda, incx));
        let count = (m4 - row).min(2048);
        let n4 = n & !3;
        let mut col = 0;
        while col < n {
            #[cfg(target_arch = "x86_64")]
            if col < n4
                && let Some(p) = &pointers
            {
                // SAFETY: four complete columns, count a multiple of four,
                // source/input extents prevalidated, AVX2/FMA guard paired.
                let sums = unsafe {
                    gemv_t_micro(
                        count,
                        p.a.add(col * lda + row),
                        lda,
                        p.x.add(row * incx),
                        incx,
                    )
                };
                for c in 0..4 {
                    // SAFETY: validated output extent, no new slice borrow
                    // invalidates same-buffer source pointers between calls.
                    unsafe {
                        let output = p.y.add((col + c) * storage.incy);
                        let previous = *output;
                        *output = if storage.incy == 1 {
                            previous + alpha * sums[c]
                        } else {
                            alpha.mul_add(sums[c], previous)
                        };
                    }
                }
                col += 4;
                continue;
            }
            #[cfg(target_arch = "x86_64")]
            if let Some(p) = &pointers {
                let width = if n & 2 != 0 && col < n4 + 2 { 2 } else { 1 };
                // SAFETY: prevalidated source extents and complete SSE2 tails.
                // Separate mul/add in these source tails is intentionally not
                // changed to FMA; only the driver output update is fused.
                unsafe {
                    let a = p.a.add(col * lda + row);
                    let x = p.x.add(row * incx);
                    let sums = if width == 2 {
                        gemv_t_tail::<2>(count, a, lda, x, incx)
                    } else {
                        gemv_t_tail::<1>(count, a, lda, x, incx)
                    };
                    for c in 0..width {
                        let output = p.y.add((col + c) * storage.incy);
                        *output = alpha.mul_add(sums[c], *output);
                    }
                }
                col += width;
                continue;
            }
            let sum = if col < n4 {
                let mut lanes = [0.0_f64; 4];
                for i in 0..count {
                    let l = i & 3;
                    lanes[l] = storage
                        .matrix(col * lda + row + i)
                        .mul_add(storage.input((row + i) * incx), lanes[l]);
                }
                reduce4(lanes)
            } else if n & 2 != 0 && col < n4 + 2 {
                let mut lanes = [0.0_f64; 2];
                for i in 0..count {
                    lanes[i & 1] +=
                        storage.matrix(col * lda + row + i) * storage.input((row + i) * incx);
                }
                lanes[0] + lanes[1]
            } else {
                let mut first = [0.0_f64; 2];
                let mut second = [0.0_f64; 2];
                for i in (0..count).step_by(4) {
                    for l in 0..2 {
                        first[l] += storage.matrix(col * lda + row + i + l)
                            * storage.input((row + i + l) * incx);
                        second[l] += storage.matrix(col * lda + row + i + 2 + l)
                            * storage.input((row + i + 2 + l) * incx);
                    }
                }
                (first[0] + second[0]) + (first[1] + second[1])
            };
            let result = if col < n4 && storage.incy == 1 {
                storage.output(col) + alpha * sum
            } else {
                alpha.mul_add(sum, storage.output(col))
            };
            storage.store(col, result);
            col += 1;
        }
    }
    let rem = m - m4;
    if rem != 0 {
        let values: [f64; 3] = core::array::from_fn(|i| {
            if i < rem {
                alpha * storage.input((m4 + i) * incx)
            } else {
                0.0
            }
        });
        for col in 0..n {
            let mut sum = storage.matrix(col * lda + m4) * values[0];
            if rem > 1 {
                sum = storage
                    .matrix(col * lda + m4)
                    .mul_add(values[0], storage.matrix(col * lda + m4 + 1) * values[1]);
            }
            if rem > 2 {
                sum = storage.matrix(col * lda + m4 + 2).mul_add(values[2], sum);
            }
            let result = if rem == 1 {
                storage
                    .matrix(col * lda + m4)
                    .mul_add(values[0], storage.output(col))
            } else {
                storage.output(col) + sum
            };
            storage.store(col, result);
        }
    }
}

pub(super) fn syr2_lower(
    n: usize,
    alpha: f64,
    x: &[f64],
    incx: usize,
    y: &[f64],
    incy: usize,
    a: &mut [f64],
    lda: usize,
) {
    // Source: driver/level2/syr2_k.c (LOWER).
    // OpenBLAS✔️❌: int CNAME(BLASLONG m, FLOAT alpha_r, FLOAT *x, BLASLONG incx,
    // OpenBLAS✔️❌: 		       FLOAT *y, BLASLONG incy, FLOAT *a, BLASLONG lda, FLOAT *buffer){
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   BLASLONG i;
    // OpenBLAS✔️❌:   FLOAT *X, *Y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   X = x;
    // OpenBLAS✔️❌:   Y = y;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incx != 1) {
    // OpenBLAS✔️❌:     COPY_K(m, x, incx, buffer, 1);
    // OpenBLAS✔️❌:     X = buffer;
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incy != 1) {
    // OpenBLAS✔️❌:     COPY_K(m, y, incy, (FLOAT *)((BLASLONG)buffer + (BUFFER_SIZE / 2)), 1);
    // OpenBLAS✔️❌:     Y = (FLOAT *)((BLASLONG)buffer + (BUFFER_SIZE / 2));
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   for (i = 0; i < m; i++){
    // OpenBLAS✔️❌: #ifndef LOWER
    // OpenBLAS✔️❌:     AXPYU_K(i + 1, 0, 0, alpha_r * X[i], Y,     1, a, 1, NULL, 0);
    // OpenBLAS✔️❌:     AXPYU_K(i + 1, 0, 0, alpha_r * Y[i], X,     1, a, 1, NULL, 0);
    // OpenBLAS✔️❌:     a += lda;
    // OpenBLAS✔️❌: #else
    // OpenBLAS✔️❌:     AXPYU_K(m - i, 0, 0, alpha_r * X[i], Y + i, 1, a, 1, NULL, 0);
    // OpenBLAS✔️❌:     AXPYU_K(m - i, 0, 0, alpha_r * Y[i], X + i, 1, a, 1, NULL, 0);
    // OpenBLAS✔️❌:     a += 1 + lda;
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:     }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   return 0;
    // OpenBLAS✔️❌: }

    if alpha == 0.0 {
        return;
    }
    // Source packs strided inputs to unit stride before the two AXPY calls.
    let xp;
    let yp;
    let xx = if incx == 1 {
        x
    } else {
        xp = (0..n).map(|i| x[i * incx]).collect::<Vec<_>>();
        &xp
    };
    let yy = if incy == 1 {
        y
    } else {
        yp = (0..n).map(|i| y[i * incy]).collect::<Vec<_>>();
        &yp
    };
    for i in 0..n {
        axpy_kernel(n - i, alpha * xx[i], &yy[i..], 1, &mut a[i * lda + i..], 1);
        axpy_kernel(n - i, alpha * yy[i], &xx[i..], 1, &mut a[i * lda + i..], 1);
    }
}

pub(super) fn ger(
    m: usize,
    n: usize,
    alpha: f64,
    x: &[f64],
    incx: usize,
    y: &[f64],
    incy: usize,
    a: &mut [f64],
    lda: usize,
) {
    // Source: kernel/x86_64/dger.c.
    // OpenBLAS✔️❌: 	 FLOAT *a, BLASLONG lda, FLOAT *buffer){
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	FLOAT *X    = x;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	if (incx != 1) {
    // OpenBLAS✔️❌:     		X = buffer;
    // OpenBLAS✔️❌:     		COPY_K(m, x, incx, X, 1);
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	BLASLONG m1 = m & -16;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	while (n > 0)
    // OpenBLAS✔️❌: 	{
    // OpenBLAS✔️❌:    		FLOAT y0 = alpha * *y;
    // OpenBLAS✔️❌: 		if ( m1 > 0 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌: 			#ifdef HAVE_KERNEL_16
    // OpenBLAS✔️❌: 				dger_kernel_16(m1, X, a, &y0);
    // OpenBLAS✔️❌: 			#else
    // OpenBLAS✔️❌:     				AXPYU_K(m1, 0, 0, y0, X, 1, a, 1, NULL, 0);
    // OpenBLAS✔️❌: 			#endif
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 		if ( m > m1 )
    // OpenBLAS✔️❌: 		{
    // OpenBLAS✔️❌:     			AXPYU_K(m-m1, 0, 0, y0, X+m1 , 1, a+m1, 1, NULL, 0);
    // OpenBLAS✔️❌: 		}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:     		a += lda;
    // OpenBLAS✔️❌:     		y += incy;
    // OpenBLAS✔️❌:     		n --;
    // OpenBLAS✔️❌: 	}
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: 	return 0;
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:

    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    for col in 0..n {
        let scale = alpha * y[col * incy];
        for row in 0..m {
            a[col * lda + row] = scale.mul_add(x[row * incx], a[col * lda + row]);
        }
    }
}

#[inline]
fn rotated(x: f64, y: f64, c: f64, s: f64) -> (f64, f64) {
    // OpenBLAS✔️✔️:         unsigned char tail_mask8 = (((unsigned char) 0xff) >> (8 -(n1&7)));
    // OpenBLAS✔️✔️: 	__m512d tail_x = _mm512_maskz_loadu_pd(*((__mmask8*) &tail_mask8), &x[tail_index_8]);
    // OpenBLAS✔️✔️: 	__m512d tail_y = _mm512_maskz_loadu_pd(*((__mmask8*) &tail_mask8), &y[tail_index_8]);
    // OpenBLAS✔️✔️: 	__m512d temp = _mm512_mul_pd(s_512, tail_y);
    // OpenBLAS✔️✔️: 	temp = _mm512_fmadd_pd(c_512, tail_x, temp);
    // OpenBLAS✔️✔️: 	_mm512_mask_storeu_pd(&x[tail_index_8],*((__mmask8*)&tail_mask8), temp);
    // OpenBLAS✔️✔️:         temp = _mm512_mul_pd(s_512, tail_x);
    // OpenBLAS✔️✔️:         temp = _mm512_fmsub_pd(c_512, tail_y, temp);
    // OpenBLAS✔️✔️:         _mm512_mask_storeu_pd(&y[tail_index_8], *((__mmask8*)&tail_mask8), temp);
    // OpenBLAS✔️✔️:     }
    // OpenBLAS✔️✔️: }
    // OpenBLAS✔️✔️: #endif

    (c.mul_add(x, s * y), c.mul_add(y, -(s * x)))
}

pub(super) fn rot_within(
    n: usize,
    a: &mut [f64],
    x_offset: usize,
    incx: usize,
    y_offset: usize,
    incy: usize,
    c: f64,
    s: f64,
) {
    // Same source kernel; offsets express disjoint strided vectors within one
    // owned backing slice without overlapping mutable suffix borrows or copies.
    for i in 0..n {
        let xi = x_offset + i * incx;
        let yi = y_offset + i * incy;
        let (x, y) = rotated(a[xi], a[yi], c, s);
        a[xi] = x;
        a[yi] = y;
    }
}

pub(super) fn trmv(
    lower: bool,
    transpose: bool,
    unit: bool,
    n: usize,
    a: &[f64],
    lda: usize,
    x: &mut [f64],
    incx: usize,
) {
    // Source: driver/level2/trmv_U.c and trmv_L.c; DTB_ENTRIES=64.
    // OpenBLAS✔️❌:   FLOAT *gemvbuffer = (FLOAT *)buffer;
    // OpenBLAS✔️❌:   FLOAT *B = b;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incb != 1) {
    // OpenBLAS✔️❌:     B = buffer;
    // OpenBLAS✔️❌:     gemvbuffer = (FLOAT *)(((BLASLONG)buffer + m * sizeof(FLOAT) + 4095) & ~4095);
    // OpenBLAS✔️❌:     COPY_K(m, b, incb, buffer, 1);
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   for (is = 0; is < m; is += DTB_ENTRIES){
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:     min_i = MIN(m - is, DTB_ENTRIES);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef TRANSA
    // OpenBLAS✔️❌:     if (is > 0){
    // OpenBLAS✔️❌:       GEMV_N(is, min_i, 0, dp1,
    // OpenBLAS✔️❌: 	     a + is * lda,  lda,
    // OpenBLAS✔️❌: 	     B + is, 1,
    // OpenBLAS✔️❌: 	     B,      1, gemvbuffer);
    // OpenBLAS✔️❌:     }
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:     for (i = 0; i < min_i; i++) {
    // OpenBLAS✔️❌:       FLOAT *AA = a + is + (i + is) * lda;
    // OpenBLAS✔️❌:       FLOAT *BB = B + is;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef TRANSA
    // OpenBLAS✔️❌:       if (i > 0) AXPYU_K(i, 0, 0, BB[i], AA, 1, BB, 1, NULL, 0);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef UNIT
    // OpenBLAS✔️❌:       BB[i] *= AA[i];
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef TRANSA
    // OpenBLAS✔️❌:       if (i < min_i - 1) BB[i] += DOTU_K(min_i - i - 1, AA + i + 1, 1, BB + i + 1, 1);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:     }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef TRANSA
    // OpenBLAS✔️❌:     if (m - is > min_i){
    // OpenBLAS✔️❌:       GEMV_T(m - is - min_i, min_i, 0, dp1,
    // OpenBLAS✔️❌: 	     a + is + min_i + is * lda,  lda,
    // OpenBLAS✔️❌: 	     B + is + min_i, 1,
    // OpenBLAS✔️❌: 	     B + is,         1, gemvbuffer);
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incb != 1) {
    // OpenBLAS✔️❌:     COPY_K(m, buffer, 1, b, incb);
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   return 0;
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   FLOAT *gemvbuffer = (FLOAT *)buffer;
    // OpenBLAS✔️❌:   FLOAT *B = b;
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incb != 1) {
    // OpenBLAS✔️❌:     B = buffer;
    // OpenBLAS✔️❌:     gemvbuffer = (FLOAT *)(((BLASLONG)buffer + m * sizeof(FLOAT) + 4095) & ~4095);
    // OpenBLAS✔️❌:     COPY_K(m, b, incb, buffer, 1);
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   for (is = m; is > 0; is -= DTB_ENTRIES){
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:     min_i = MIN(is, DTB_ENTRIES);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef TRANSA
    // OpenBLAS✔️❌:     if (m - is > 0){
    // OpenBLAS✔️❌:       GEMV_N(m - is, min_i, 0, dp1,
    // OpenBLAS✔️❌: 	     a + is + (is - min_i) * lda, lda,
    // OpenBLAS✔️❌: 	     B + is - min_i, 1,
    // OpenBLAS✔️❌: 	     B + is, 1, gemvbuffer);
    // OpenBLAS✔️❌:     }
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:     for (i = 0; i < min_i; i++) {
    // OpenBLAS✔️❌:       FLOAT *AA = a + (is - i - 1) + (is - i - 1) * lda;
    // OpenBLAS✔️❌:       FLOAT *BB = B + (is - i - 1);
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef TRANSA
    // OpenBLAS✔️❌:       if (i > 0) AXPYU_K(i, 0, 0, BB[0], AA + 1, 1, BB + 1, 1, NULL, 0);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifndef UNIT
    // OpenBLAS✔️❌:       BB[0] *= AA[0];
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef TRANSA
    // OpenBLAS✔️❌:       if (i < min_i - 1) BB[0] += DOTU_K(min_i - i - 1, AA - (min_i - i - 1), 1, BB - (min_i - i - 1), 1);
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:     }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌: #ifdef TRANSA
    // OpenBLAS✔️❌:     if (is - min_i > 0){
    // OpenBLAS✔️❌:       GEMV_T(is - min_i, min_i, 0, dp1,
    // OpenBLAS✔️❌: 	     a +  (is - min_i) * lda, lda,
    // OpenBLAS✔️❌: 	     B,              1,
    // OpenBLAS✔️❌: 	     B + is - min_i, 1, gemvbuffer);
    // OpenBLAS✔️❌:     }
    // OpenBLAS✔️❌: #endif
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   if (incb != 1) {
    // OpenBLAS✔️❌:     COPY_K(m, buffer, 1, b, incb);
    // OpenBLAS✔️❌:   }
    // OpenBLAS✔️❌:
    // OpenBLAS✔️❌:   return 0;
    // OpenBLAS✔️❌: }
    // OpenBLAS✔️❌:

    if n == 0 {
        return;
    }
    if incx != 1 {
        let mut packed: Vec<_> = (0..n).map(|i| x[i * incx]).collect();
        trmv(lower, transpose, unit, n, a, lda, &mut packed, 1);
        for i in 0..n {
            x[i * incx] = packed[i];
        }
        return;
    }
    // common.h: DTB_DEFAULT_ENTRIES = 64; the pinned dynamic table inherits it.
    const BLOCK: usize = 64;
    if lower != transpose {
        let mut end = n;
        while end != 0 {
            let start = end.saturating_sub(BLOCK);
            let count = end - start;
            if !transpose {
                if end < n {
                    let (left, right) = x.split_at_mut(end);
                    gemv(
                        false,
                        n - end,
                        count,
                        1.0,
                        &a[start * lda + end..],
                        lda,
                        &left[start..],
                        1,
                        1.0,
                        right,
                        1,
                    );
                }
                for idx in (start..end).rev() {
                    let value = x[idx];
                    if idx + 1 < end {
                        axpy_kernel(
                            end - idx - 1,
                            value,
                            &a[idx * lda + idx + 1..],
                            1,
                            &mut x[idx + 1..],
                            1,
                        );
                    }
                    if !unit {
                        x[idx] *= a[idx * lda + idx];
                    }
                }
            } else {
                for idx in (start..end).rev() {
                    if !unit {
                        x[idx] *= a[idx * lda + idx];
                    }
                    if idx > start {
                        x[idx] += dot(idx - start, &a[idx * lda + start..], 1, &x[start..], 1);
                    }
                }
                if start > 0 {
                    let (left, right) = x.split_at_mut(start);
                    gemv(
                        true,
                        start,
                        count,
                        1.0,
                        &a[start * lda..],
                        lda,
                        left,
                        1,
                        1.0,
                        right,
                        1,
                    );
                }
            }
            end = start;
        }
    } else {
        for start in (0..n).step_by(BLOCK) {
            let end = (start + BLOCK).min(n);
            let count = end - start;
            if !transpose {
                if start > 0 {
                    let (left, right) = x.split_at_mut(start);
                    gemv(
                        false,
                        start,
                        count,
                        1.0,
                        &a[start * lda..],
                        lda,
                        right,
                        1,
                        1.0,
                        left,
                        1,
                    );
                }
                for idx in start..end {
                    let value = x[idx];
                    if idx > start {
                        axpy_kernel(
                            idx - start,
                            value,
                            &a[idx * lda + start..],
                            1,
                            &mut x[start..],
                            1,
                        );
                    }
                    if !unit {
                        x[idx] *= a[idx * lda + idx];
                    }
                }
            } else {
                for idx in start..end {
                    if !unit {
                        x[idx] *= a[idx * lda + idx];
                    }
                    if idx + 1 < end {
                        x[idx] += dot(
                            end - idx - 1,
                            &a[idx * lda + idx + 1..],
                            1,
                            &x[idx + 1..],
                            1,
                        );
                    }
                }
                if end < n {
                    let (left, right) = x.split_at_mut(end);
                    gemv(
                        true,
                        n - end,
                        count,
                        1.0,
                        &a[start * lda + end..],
                        lda,
                        right,
                        1,
                        1.0,
                        &mut left[start..],
                        1,
                    );
                }
            }
        }
    }
}
