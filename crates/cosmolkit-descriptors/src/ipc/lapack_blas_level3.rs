// Copyright (c) 2011-2014, The OpenBLAS Project
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//    1. Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in
//       the documentation and/or other materials provided with the
//       distribution.
//    3. Neither the name of the OpenBLAS project nor the names of
//       its contributors may be used to endorse or promote products
//       derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//! Single-thread SkylakeX BLAS3 numerical paths from OpenBLAS
//! 446c436e10450c348169808f1a6b3fae0925c9f7. Storage is column major.
//! FMA lanes, source reductions and packed workspaces are explicit. The source
//! NN tail buffer is retained; inputs otherwise use their original storage.
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Transpose {
    NoTrans,
    Trans,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Side {
    Left,
    Right,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Triangle {
    Upper,
    Lower,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Diagonal {
    Unit,
    NonUnit,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum Blas3Error {
    InvalidArgument {
        routine: &'static str,
        argument: u8,
    },
    StorageTooShort {
        routine: &'static str,
        matrix: &'static str,
        required: usize,
        actual: usize,
    },
    DimensionOverflow,
}
impl fmt::Display for Blas3Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidArgument { routine, argument } => {
                write!(f, "{routine} argument {argument} is invalid")
            }
            Self::StorageTooShort {
                routine,
                matrix,
                required,
                actual,
            } => write!(
                f,
                "{routine} {matrix} storage requires {required} entries, got {actual}"
            ),
            Self::DimensionOverflow => {
                f.write_str("BLAS3 matrix dimensions overflow the index range")
            }
        }
    }
}
impl std::error::Error for Blas3Error {}
const P: usize = 192;
const Q: usize = 384;
const R: usize = 10704;

type RowFma = fn(&[f64], f64, &mut [f64]);
fn row_fma() -> RowFma {
    row_fma_portable
}
fn row_fma_portable(a: &[f64], b: f64, c: &mut [f64]) {
    for (av, cv) in a.iter().zip(c) {
        *cv = av.mul_add(b, *cv);
    }
}
pub(super) use super::lapack_blas::MatrixSource;
impl<'a> MatrixSource<'a> {
    fn matrix(self, ld: usize, trans: Transpose) -> Matrix<'a> {
        match self {
            Self::Separate(data) => Matrix {
                separate: Some(data),
                offset: 0,
                ld,
                trans,
            },
            Self::Within(offset) => Matrix {
                separate: None,
                offset,
                ld,
                trans,
            },
        }
    }
}
#[derive(Clone, Copy)]
struct Matrix<'a> {
    separate: Option<&'a [f64]>,
    offset: usize,
    ld: usize,
    trans: Transpose,
}
impl<'a> Matrix<'a> {
    #[inline]
    fn get(self, shared: &[f64], row: usize, col: usize) -> f64 {
        let (row, col) = if self.trans == Transpose::NoTrans {
            (row, col)
        } else {
            (col, row)
        };
        self.separate.unwrap_or(shared)[self.offset + row + col * self.ld]
    }
}
fn leading(routine: &'static str, rows: usize, ld: usize, argument: u8) -> Result<(), Blas3Error> {
    if ld < rows.max(1) {
        Err(Blas3Error::InvalidArgument { routine, argument })
    } else {
        Ok(())
    }
}
fn storage(
    routine: &'static str,
    name: &'static str,
    rows: usize,
    cols: usize,
    ld: usize,
    offset: usize,
    len: usize,
    arg: u8,
) -> Result<(), Blas3Error> {
    leading(routine, rows, ld, arg)?;
    if rows == 0 || cols == 0 {
        return Ok(());
    }
    let required = (cols - 1)
        .checked_mul(ld)
        .and_then(|x| x.checked_add(rows))
        .and_then(|x| x.checked_add(offset))
        .filter(|&x| x <= isize::MAX as usize)
        .ok_or(Blas3Error::DimensionOverflow)?;
    if required > len {
        return Err(Blas3Error::StorageTooShort {
            routine,
            matrix: name,
            required,
            actual: len,
        });
    }
    Ok(())
}
fn op_storage(
    routine: &'static str,
    name: &'static str,
    rows: usize,
    cols: usize,
    a: Matrix<'_>,
    shared_len: usize,
    arg: u8,
) -> Result<(), Blas3Error> {
    let (rows, cols) = if a.trans == Transpose::NoTrans {
        (rows, cols)
    } else {
        (cols, rows)
    };
    storage(
        routine,
        name,
        rows,
        cols,
        a.ld,
        a.offset,
        a.separate.map_or(shared_len, <[f64]>::len),
        arg,
    )
}
fn independent(
    rows_a: usize,
    cols_a: usize,
    a: Matrix<'_>,
    rows_c: usize,
    cols_c: usize,
    c_offset: usize,
    ldc: usize,
    routine: &'static str,
    arg: u8,
) -> Result<(), Blas3Error> {
    if rows_a == 0 || cols_a == 0 || rows_c == 0 || cols_c == 0 {
        return Ok(());
    }
    if a.separate.is_some() {
        return Ok(());
    }
    let (rows_a, cols_a) = if a.trans == Transpose::NoTrans {
        (rows_a, cols_a)
    } else {
        (cols_a, rows_a)
    };
    let (mut i, mut j) = (0, 0);
    while i < cols_a && j < cols_c {
        let aa = a.offset + i * a.ld;
        let cc = c_offset + j * ldc;
        if aa + rows_a <= cc {
            i += 1;
        } else if cc + rows_c <= aa {
            j += 1;
        } else {
            return Err(Blas3Error::InvalidArgument {
                routine,
                argument: arg,
            });
        }
    }
    Ok(())
}
#[inline]
fn row_tile(m: usize) -> usize {
    if m >= 16 {
        16
    } else if m >= 8 {
        8
    } else if m >= 4 {
        4
    } else if m >= 2 {
        2
    } else {
        1
    }
}
#[inline]
fn col_tile(n: usize) -> usize {
    if n >= 12 {
        12
    } else if n >= 10 {
        10
    } else if n >= 8 {
        8
    } else if n >= 6 {
        6
    } else if n >= 4 {
        4
    } else if n >= 2 {
        2
    } else {
        1
    }
}
// The actual DYNAMIC_ARCH driver is compiled without SKYLAKEX defined:
// its source #else dispatch uses 3 * GEMM_UNROLL_N (6), then 2, then 1.
fn driver_columns(n: usize) -> usize {
    if n >= 6 {
        6
    } else if n > 2 {
        2
    } else {
        n
    }
}
fn driver_chunks(mut start: usize, end: usize) -> impl Iterator<Item = (usize, usize)> {
    std::iter::from_fn(move || {
        if start >= end {
            None
        } else {
            let w = driver_columns(end - start);
            let old = start;
            start += w;
            Some((old, w))
        }
    })
}
fn panel_size(m: usize, align: usize) -> usize {
    if m >= 2 * P {
        P
    } else if m > P {
        (m / 2).div_ceil(align) * align
    } else {
        m
    }
}
fn scale(m: usize, n: usize, beta: f64, c: &mut [f64], offset: usize, ldc: usize) {
    // OpenBLAS❗✔️: /*********************************************************************/
    // OpenBLAS❗✔️: /* Copyright 2009, 2010 The University of Texas at Austin.           */
    // OpenBLAS❗✔️: /* All rights reserved.                                              */
    // OpenBLAS❗✔️: /*                                                                   */
    // OpenBLAS❗✔️: /* Redistribution and use in source and binary forms, with or        */
    // OpenBLAS❗✔️: /* without modification, are permitted provided that the following   */
    // OpenBLAS❗✔️: /* conditions are met:                                               */
    // OpenBLAS❗✔️: /*                                                                   */
    // OpenBLAS❗✔️: /*   1. Redistributions of source code must retain the above         */
    // OpenBLAS❗✔️: /*      copyright notice, this list of conditions and the following  */
    // OpenBLAS❗✔️: /*      disclaimer.                                                  */
    // OpenBLAS❗✔️: /*                                                                   */
    // OpenBLAS❗✔️: /*   2. Redistributions in binary form must reproduce the above      */
    // OpenBLAS❗✔️: /*      copyright notice, this list of conditions and the following  */
    // OpenBLAS❗✔️: /*      disclaimer in the documentation and/or other materials       */
    // OpenBLAS❗✔️: /*      provided with the distribution.                              */
    // OpenBLAS❗✔️: /*                                                                   */
    // OpenBLAS❗✔️: /*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
    // OpenBLAS❗✔️: /*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
    // OpenBLAS❗✔️: /*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
    // OpenBLAS❗✔️: /*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
    // OpenBLAS❗✔️: /*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
    // OpenBLAS❗✔️: /*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
    // OpenBLAS❗✔️: /*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
    // OpenBLAS❗✔️: /*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
    // OpenBLAS❗✔️: /*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
    // OpenBLAS❗✔️: /*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
    // OpenBLAS❗✔️: /*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
    // OpenBLAS❗✔️: /*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
    // OpenBLAS❗✔️: /*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
    // OpenBLAS❗✔️: /*    POSSIBILITY OF SUCH DAMAGE.                                    */
    // OpenBLAS❗✔️: /*                                                                   */
    // OpenBLAS❗✔️: /* The views and conclusions contained in the software and           */
    // OpenBLAS❗✔️: /* documentation are those of the authors and should not be          */
    // OpenBLAS❗✔️: /* interpreted as representing official policies, either expressed   */
    // OpenBLAS❗✔️: /* or implied, of The University of Texas at Austin.                 */
    // OpenBLAS❗✔️: /*********************************************************************/
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #include "common.h"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #include <immintrin.h>
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: int CNAME(BLASLONG m, BLASLONG n, BLASLONG dummy1, FLOAT beta,
    // OpenBLAS❗✔️: 	  FLOAT *dummy2, BLASLONG dummy3, FLOAT *dummy4, BLASLONG dummy5,
    // OpenBLAS❗✔️: 	  FLOAT *c, BLASLONG ldc){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   BLASLONG i, j;
    // OpenBLAS❗✔️:   FLOAT *c_offset1, *c_offset;
    // OpenBLAS❗✔️:   FLOAT ctemp1, ctemp2, ctemp3, ctemp4;
    // OpenBLAS❗✔️:   FLOAT ctemp5, ctemp6, ctemp7, ctemp8;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   /* fast path.. just zero the whole matrix */
    // OpenBLAS❗✔️:   if (m == ldc && beta == ZERO) {
    // OpenBLAS❗✔️: 	memset(c, 0, m * n * sizeof(FLOAT));
    // OpenBLAS❗✔️: 	return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (m == 0 || n == 0)
    // OpenBLAS❗✔️: 	return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   c_offset = c;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (beta == ZERO){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     j = n;
    // OpenBLAS❗✔️:     do {
    // OpenBLAS❗✔️:       c_offset1 = c_offset;
    // OpenBLAS❗✔️:       c_offset += ldc;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = m;
    // OpenBLAS❗✔️: #ifdef __AVX2__
    // OpenBLAS❗✔️: #ifdef __AVX512CD__
    // OpenBLAS❗✔️:       while (i >= 32) {
    // OpenBLAS❗✔️: 	  __m512d z_zero = _mm512_setzero_pd();
    // OpenBLAS❗✔️: 	  _mm512_storeu_pd(c_offset1, z_zero);
    // OpenBLAS❗✔️: 	  _mm512_storeu_pd(c_offset1 + 8, z_zero);
    // OpenBLAS❗✔️: 	  _mm512_storeu_pd(c_offset1 + 16, z_zero);
    // OpenBLAS❗✔️: 	  _mm512_storeu_pd(c_offset1 + 24 , z_zero);
    // OpenBLAS❗✔️: 	  c_offset1 += 32;
    // OpenBLAS❗✔️: 	  i -= 32;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:       while (i >= 8) {
    // OpenBLAS❗✔️: #ifdef __AVX512CD__
    // OpenBLAS❗✔️: 	  __m512d z_zero = _mm512_setzero_pd();
    // OpenBLAS❗✔️: 	  _mm512_storeu_pd(c_offset1, z_zero);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	 __m256d y_zero = _mm256_setzero_pd();
    // OpenBLAS❗✔️: 	 _mm256_storeu_pd(c_offset1, y_zero);
    // OpenBLAS❗✔️: 	 _mm256_storeu_pd(c_offset1 + 4, y_zero);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	  c_offset1 += 8;
    // OpenBLAS❗✔️: 	  i -= 8;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:       while (i > 0) {
    // OpenBLAS❗✔️: 	  *c_offset1 = ZERO;
    // OpenBLAS❗✔️: 	  c_offset1 ++;
    // OpenBLAS❗✔️: 	  i --;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       j --;
    // OpenBLAS❗✔️:     } while (j > 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   } else {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     j = n;
    // OpenBLAS❗✔️:     do {
    // OpenBLAS❗✔️:       c_offset1 = c_offset;
    // OpenBLAS❗✔️:       c_offset += ldc;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m >> 3);
    // OpenBLAS❗✔️:       if (i > 0){
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  ctemp1 = *(c_offset1 + 0);
    // OpenBLAS❗✔️: 	  ctemp2 = *(c_offset1 + 1);
    // OpenBLAS❗✔️: 	  ctemp3 = *(c_offset1 + 2);
    // OpenBLAS❗✔️: 	  ctemp4 = *(c_offset1 + 3);
    // OpenBLAS❗✔️: 	  ctemp5 = *(c_offset1 + 4);
    // OpenBLAS❗✔️: 	  ctemp6 = *(c_offset1 + 5);
    // OpenBLAS❗✔️: 	  ctemp7 = *(c_offset1 + 6);
    // OpenBLAS❗✔️: 	  ctemp8 = *(c_offset1 + 7);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp1 *= beta;
    // OpenBLAS❗✔️: 	  ctemp2 *= beta;
    // OpenBLAS❗✔️: 	  ctemp3 *= beta;
    // OpenBLAS❗✔️: 	  ctemp4 *= beta;
    // OpenBLAS❗✔️: 	  ctemp5 *= beta;
    // OpenBLAS❗✔️: 	  ctemp6 *= beta;
    // OpenBLAS❗✔️: 	  ctemp7 *= beta;
    // OpenBLAS❗✔️: 	  ctemp8 *= beta;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  *(c_offset1 + 0) = ctemp1;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 1) = ctemp2;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 2) = ctemp3;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 3) = ctemp4;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 4) = ctemp5;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 5) = ctemp6;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 6) = ctemp7;
    // OpenBLAS❗✔️: 	  *(c_offset1 + 7) = ctemp8;
    // OpenBLAS❗✔️: 	  c_offset1 += 8;
    // OpenBLAS❗✔️: 	  i --;
    // OpenBLAS❗✔️: 	} while (i > 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m & 7);
    // OpenBLAS❗✔️:       if (i > 0){
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  ctemp1 = *c_offset1;
    // OpenBLAS❗✔️: 	  ctemp1 *= beta;
    // OpenBLAS❗✔️: 	  *c_offset1 = ctemp1;
    // OpenBLAS❗✔️: 	  c_offset1 ++;
    // OpenBLAS❗✔️: 	  i --;
    // OpenBLAS❗✔️: 	} while (i > 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       j --;
    // OpenBLAS❗✔️:     } while (j > 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:   return 0;
    // OpenBLAS❗✔️: };
    if beta == 1. {
        return;
    }
    for j in 0..n {
        for i in 0..m {
            let x = &mut c[offset + i + j * ldc];
            *x = if beta == 0. { 0. } else { *x * beta };
        }
    }
}
fn pack_a(
    m: usize,
    k: usize,
    a: Matrix<'_>,
    shared: &[f64],
    row: usize,
    col: usize,
    out: &mut [f64],
) {
    // OpenBLAS❗✔️:   IFLOAT *aoffset;
    // OpenBLAS❗✔️:   IFLOAT *aoffset1, *aoffset2, *aoffset3, *aoffset4;
    // OpenBLAS❗✔️:   IFLOAT *aoffset5, *aoffset6, *aoffset7, *aoffset8;
    // OpenBLAS❗✔️:   IFLOAT *aoffset9, *aoffset10, *aoffset11, *aoffset12;
    // OpenBLAS❗✔️:   IFLOAT *aoffset13, *aoffset14, *aoffset15, *aoffset16;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   IFLOAT *boffset;
    // OpenBLAS❗✔️:   IFLOAT ctemp01, ctemp02, ctemp03, ctemp04;
    // OpenBLAS❗✔️:   IFLOAT ctemp05, ctemp06, ctemp07, ctemp08;
    // OpenBLAS❗✔️:   IFLOAT ctemp09, ctemp10, ctemp11, ctemp12;
    // OpenBLAS❗✔️:   IFLOAT ctemp13, ctemp14, ctemp15, ctemp16;
    // OpenBLAS❗✔️:   IFLOAT ctemp17, ctemp18, ctemp19, ctemp20;
    // OpenBLAS❗✔️:   IFLOAT ctemp21, ctemp22, ctemp23, ctemp24;
    // OpenBLAS❗✔️:   IFLOAT ctemp25, ctemp26, ctemp27, ctemp28;
    // OpenBLAS❗✔️:   IFLOAT ctemp29, ctemp30, ctemp31, ctemp32;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   aoffset = a;
    // OpenBLAS❗✔️:   boffset = b;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   j = (n >> 4);
    // OpenBLAS❗✔️:   if (j > 0){
    // OpenBLAS❗✔️:     do{
    // OpenBLAS❗✔️:       aoffset1  = aoffset;
    // OpenBLAS❗✔️:       aoffset2  = aoffset1  + lda;
    // OpenBLAS❗✔️:       aoffset3  = aoffset2  + lda;
    // OpenBLAS❗✔️:       aoffset4  = aoffset3  + lda;
    // OpenBLAS❗✔️:       aoffset5  = aoffset4  + lda;
    // OpenBLAS❗✔️:       aoffset6  = aoffset5  + lda;
    // OpenBLAS❗✔️:       aoffset7  = aoffset6  + lda;
    // OpenBLAS❗✔️:       aoffset8  = aoffset7  + lda;
    // OpenBLAS❗✔️:       aoffset9  = aoffset8  + lda;
    // OpenBLAS❗✔️:       aoffset10 = aoffset9  + lda;
    // OpenBLAS❗✔️:       aoffset11 = aoffset10 + lda;
    // OpenBLAS❗✔️:       aoffset12 = aoffset11 + lda;
    // OpenBLAS❗✔️:       aoffset13 = aoffset12 + lda;
    // OpenBLAS❗✔️:       aoffset14 = aoffset13 + lda;
    // OpenBLAS❗✔️:       aoffset15 = aoffset14 + lda;
    // OpenBLAS❗✔️:       aoffset16 = aoffset15 + lda;
    // OpenBLAS❗✔️:       aoffset += 16 * lda;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m >> 1);
    // OpenBLAS❗✔️:       if (i > 0){
    // OpenBLAS❗✔️: 	do{
    // OpenBLAS❗✔️: 	  ctemp01 = *(aoffset1 +  0);
    // OpenBLAS❗✔️: 	  ctemp02 = *(aoffset1 +  1);
    // OpenBLAS❗✔️: 	  ctemp03 = *(aoffset2 +  0);
    // OpenBLAS❗✔️: 	  ctemp04 = *(aoffset2 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp05 = *(aoffset3 +  0);
    // OpenBLAS❗✔️: 	  ctemp06 = *(aoffset3 +  1);
    // OpenBLAS❗✔️: 	  ctemp07 = *(aoffset4 +  0);
    // OpenBLAS❗✔️: 	  ctemp08 = *(aoffset4 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp09 = *(aoffset5 +  0);
    // OpenBLAS❗✔️: 	  ctemp10 = *(aoffset5 +  1);
    // OpenBLAS❗✔️: 	  ctemp11 = *(aoffset6 +  0);
    // OpenBLAS❗✔️: 	  ctemp12 = *(aoffset6 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp13 = *(aoffset7 +  0);
    // OpenBLAS❗✔️: 	  ctemp14 = *(aoffset7 +  1);
    // OpenBLAS❗✔️: 	  ctemp15 = *(aoffset8 +  0);
    // OpenBLAS❗✔️: 	  ctemp16 = *(aoffset8 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp17 = *(aoffset9 +  0);
    // OpenBLAS❗✔️: 	  ctemp18 = *(aoffset9 +  1);
    // OpenBLAS❗✔️: 	  ctemp19 = *(aoffset10 +  0);
    // OpenBLAS❗✔️: 	  ctemp20 = *(aoffset10 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp21 = *(aoffset11 +  0);
    // OpenBLAS❗✔️: 	  ctemp22 = *(aoffset11 +  1);
    // OpenBLAS❗✔️: 	  ctemp23 = *(aoffset12 +  0);
    // OpenBLAS❗✔️: 	  ctemp24 = *(aoffset12 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp25 = *(aoffset13 +  0);
    // OpenBLAS❗✔️: 	  ctemp26 = *(aoffset13 +  1);
    // OpenBLAS❗✔️: 	  ctemp27 = *(aoffset14 +  0);
    // OpenBLAS❗✔️: 	  ctemp28 = *(aoffset14 +  1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ctemp29 = *(aoffset15 +  0);
    // OpenBLAS❗✔️: 	  ctemp30 = *(aoffset15 +  1);
    // OpenBLAS❗✔️: 	  ctemp31 = *(aoffset16 +  0);
    // OpenBLAS❗✔️:       _mm512_storeu_pd(dst1+16,z3); _mm512_storeu_pd(dst1+24,z4); dst1 += 8 * dim_second;
    // OpenBLAS❗✔️:       dim1_count -= 8;
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:     dst1 -= 4 * (dim_second - dim2_count);
    // OpenBLAS❗✔️:     if(dim1_count>3){
    // OpenBLAS❗✔️:       y1 = _mm256_loadu_pd(src1); src1 += 4;
    // OpenBLAS❗✔️:       y2 = _mm256_loadu_pd(src2); src2 += 4;
    // OpenBLAS❗✔️:       y3 = _mm256_loadu_pd(src3); src3 += 4;
    // OpenBLAS❗✔️:       y4 = _mm256_loadu_pd(src4); src4 += 4;
    // OpenBLAS❗✔️:       _mm256_storeu_pd(dst1+ 0,y1); _mm256_storeu_pd(dst1+ 4,y2);
    // OpenBLAS❗✔️:       _mm256_storeu_pd(dst1+ 8,y3); _mm256_storeu_pd(dst1+12,y4); dst1 += 4 * dim_second;
    // OpenBLAS❗✔️:       dim1_count -= 4;
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:     dst1 -= 2 * (dim_second - dim2_count);
    // OpenBLAS❗✔️:     if(dim1_count>1){
    // OpenBLAS❗✔️:       x1 = _mm_loadu_pd(src1); src1 += 2;
    // OpenBLAS❗✔️:       x2 = _mm_loadu_pd(src2); src2 += 2;
    // OpenBLAS❗✔️:       x3 = _mm_loadu_pd(src3); src3 += 2;
    // OpenBLAS❗✔️:       x4 = _mm_loadu_pd(src4); src4 += 2;
    // OpenBLAS❗✔️:       _mm_storeu_pd(dst1+0,x1); _mm_storeu_pd(dst1+2,x2);
    // OpenBLAS❗✔️:       _mm_storeu_pd(dst1+4,x3); _mm_storeu_pd(dst1+6,x4); dst1 += 2 * dim_second;
    // OpenBLAS❗✔️:       dim1_count -= 2;
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:     dst1 -= dim_second - dim2_count;
    // OpenBLAS❗✔️:     if(dim1_count>0){
    // OpenBLAS❗✔️:       s1 = *src1; src1++; s2 = *src2; src2++; s3 = *src3; src3++; s4 = *src4; src4++;
    // OpenBLAS❗✔️:       dst1[0] = s1; dst1[1] = s2; dst1[2] = s3; dst1[3] = s4;
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:     src1 += src_inc; src2 += src_inc; src3 += src_inc; src4 += src_inc;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:   src_inc -= 2 * lead_dim;
    // OpenBLAS❗✔️:   for(; dim2_count>1; dim2_count-=2){
    // OpenBLAS❗✔️:     dst1 = dst + 16 * (dim_second - dim2_count);
    // OpenBLAS❗✔️:     for(dim1_count=dim_first; dim1_count>15; dim1_count-=16){
    // OpenBLAS❗✔️:       z1 = _mm512_loadu_pd(src1); z2 = _mm512_loadu_pd(src1+8); src1 += 16;
    // OpenBLAS❗✔️:       z3 = _mm512_loadu_pd(src2); z4 = _mm512_loadu_pd(src2+8); src2 += 16;
    // OpenBLAS❗✔️:       _mm512_storeu_pd(dst1+ 0,z1); _mm512_storeu_pd(dst1+ 8,z2);
    // OpenBLAS❗✔️:       _mm512_storeu_pd(dst1+16,z3); _mm512_storeu_pd(dst1+24,z4); dst1 += 16 * dim_second;
    // OpenBLAS❗✔️:     }
    let (mut i, mut p) = (0, 0);
    while i < m {
        let h = row_tile(m - i);
        if a.trans == Transpose::NoTrans {
            let source = a.separate.unwrap_or(shared);
            for l in 0..k {
                let begin = a.offset + row + i + (col + l) * a.ld;
                out[p..p + h].copy_from_slice(&source[begin..begin + h]);
                p += h;
            }
        } else {
            for l in 0..k {
                for ii in 0..h {
                    out[p] = a.get(shared, row + i + ii, col + l);
                    p += 1;
                }
            }
        }
        i += h;
    }
}
fn pack_b(
    k: usize,
    n: usize,
    b: Matrix<'_>,
    shared: &[f64],
    row: usize,
    col: usize,
    out: &mut [f64],
) {
    // OpenBLAS❗✔️: #include "common.h"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b){
    // OpenBLAS❗✔️:   BLASLONG i, j;
    // OpenBLAS❗✔️:   IFLOAT *a_offset, *a_offset1, *a_offset2;
    // OpenBLAS❗✔️:   IFLOAT *b_offset;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   a_offset = a;
    // OpenBLAS❗✔️:   b_offset = b;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   j = (n >> 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (j > 0){
    // OpenBLAS❗✔️:     do {
    // OpenBLAS❗✔️:       a_offset1 = a_offset;
    // OpenBLAS❗✔️:       a_offset2 = a_offset + lda;
    // OpenBLAS❗✔️:       a_offset += 2 * lda;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m >> 2);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (i > 0){
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  *(b_offset + 0) = *(a_offset1 + 0);
    // OpenBLAS❗✔️: 	  *(b_offset + 1) = *(a_offset2 + 0);
    // OpenBLAS❗✔️: 	  *(b_offset + 2) = *(a_offset1 + 1);
    // OpenBLAS❗✔️: 	  *(b_offset + 3) = *(a_offset2 + 1);
    // OpenBLAS❗✔️: 	  *(b_offset + 4) = *(a_offset1 + 2);
    // OpenBLAS❗✔️: 	  *(b_offset + 5) = *(a_offset2 + 2);
    // OpenBLAS❗✔️: 	  *(b_offset + 6) = *(a_offset1 + 3);
    // OpenBLAS❗✔️: 	  *(b_offset + 7) = *(a_offset2 + 3);
    // OpenBLAS❗✔️: 	  a_offset1 += 4;
    // OpenBLAS❗✔️: 	  a_offset2 += 4;
    // OpenBLAS❗✔️: 	  b_offset += 8;
    // OpenBLAS❗✔️: 	  i --;
    // OpenBLAS❗✔️: 	} while (i > 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m & 3);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (i > 0){
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  *(b_offset + 0) = *(a_offset1 + 0);
    // OpenBLAS❗✔️: 	  *(b_offset + 1) = *(a_offset2 + 0);
    // OpenBLAS❗✔️: 	  a_offset1 ++;
    // OpenBLAS❗✔️: 	  a_offset2 ++;
    // OpenBLAS❗✔️: 	  b_offset += 2;
    // OpenBLAS❗✔️: 	  i --;
    // OpenBLAS❗✔️: 	} while (i > 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       j --;
    // OpenBLAS❗✔️:     } while (j > 0);
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (n & 1){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     i = (m >> 3);
    // OpenBLAS❗✔️:     if (i > 0){
    // OpenBLAS❗✔️:       do {
    // OpenBLAS❗✔️: 	*(b_offset + 0) = *(a_offset + 0);
    // OpenBLAS❗✔️: 	*(b_offset + 1) = *(a_offset + 1);
    // OpenBLAS❗✔️: 	*(b_offset + 2) = *(a_offset + 2);
    // OpenBLAS❗✔️: 	*(b_offset + 3) = *(a_offset + 3);
    // OpenBLAS❗✔️: 	*(b_offset + 4) = *(a_offset + 4);
    // OpenBLAS❗✔️: 	*(b_offset + 5) = *(a_offset + 5);
    // OpenBLAS❗✔️: 	*(b_offset + 6) = *(a_offset + 6);
    // OpenBLAS❗✔️: 	*(b_offset + 7) = *(a_offset + 7);
    // OpenBLAS❗✔️: 	a_offset += 8;
    // OpenBLAS❗✔️: 	b_offset += 8;
    // OpenBLAS❗✔️: 	i --;
    // OpenBLAS❗✔️:       } while (i > 0);
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     i = (m & 7);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     if (i > 0){
    // OpenBLAS❗✔️:       do {
    // OpenBLAS❗✔️: 	*(b_offset + 0) = *(a_offset + 0);
    // OpenBLAS❗✔️: 	a_offset ++;
    // OpenBLAS❗✔️: 	b_offset ++;
    // OpenBLAS❗✔️: 	i --;
    // OpenBLAS❗✔️:       } while (i > 0);
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   return 0;
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️: #include "common.h"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: int CNAME(BLASLONG m, BLASLONG n, IFLOAT *a, BLASLONG lda, IFLOAT *b){
    // OpenBLAS❗✔️:   BLASLONG i, j;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   IFLOAT *a_offset, *a_offset1, *a_offset2;
    // OpenBLAS❗✔️:   IFLOAT *b_offset, *b_offset1, *b_offset2;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   a_offset = a;
    // OpenBLAS❗✔️:   b_offset = b;
    // OpenBLAS❗✔️:   b_offset2 = b + m * (n & ~1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   i = (m >> 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (i > 0) {
    // OpenBLAS❗✔️:     do {
    // OpenBLAS❗✔️:       a_offset1 = a_offset;
    // OpenBLAS❗✔️:       a_offset2 = a_offset + lda;
    // OpenBLAS❗✔️:       a_offset += 2 * lda;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       b_offset1 = b_offset;
    // OpenBLAS❗✔️:       b_offset += 4;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       j = (n >> 1);
    // OpenBLAS❗✔️:       if (j > 0){
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  *(b_offset1 + 0) = *(a_offset1 + 0);
    // OpenBLAS❗✔️: 	  *(b_offset1 + 1) = *(a_offset1 + 1);
    // OpenBLAS❗✔️: 	  *(b_offset1 + 2) = *(a_offset2 + 0);
    // OpenBLAS❗✔️: 	  *(b_offset1 + 3) = *(a_offset2 + 1);
    // OpenBLAS❗✔️: 	  a_offset1 += 2;
    // OpenBLAS❗✔️: 	  a_offset2 += 2;
    // OpenBLAS❗✔️: 	  b_offset1 += m * 2;
    // OpenBLAS❗✔️: 	  j--;
    // OpenBLAS❗✔️: 	} while (j > 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (n & 1){
    // OpenBLAS❗✔️: 	  *(b_offset2 + 0) = *(a_offset1 + 0);
    // OpenBLAS❗✔️: 	  *(b_offset2 + 1) = *(a_offset2 + 0);
    // OpenBLAS❗✔️: 	  b_offset2 += 2;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       i --;
    // OpenBLAS❗✔️:     } while (i > 0);
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (m & 1) {
    // OpenBLAS❗✔️:     j = (n >> 1);
    // OpenBLAS❗✔️:     if (j > 0){
    // OpenBLAS❗✔️:       do {
    // OpenBLAS❗✔️: 	*(b_offset + 0) = *(a_offset + 0);
    // OpenBLAS❗✔️: 	*(b_offset + 1) = *(a_offset + 1);
    // OpenBLAS❗✔️: 	a_offset += 2;
    // OpenBLAS❗✔️: 	b_offset += m * 2;
    // OpenBLAS❗✔️: 	j--;
    // OpenBLAS❗✔️:       } while (j > 0);
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     if (n & 1){
    // OpenBLAS❗✔️:       *(b_offset2 + 0) = *(a_offset + 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   return 0;
    // OpenBLAS❗✔️: }
    let (mut j, mut p) = (0, 0);
    while j < n {
        let w = if n - j >= 2 { 2 } else { 1 };
        let source = b.separate.unwrap_or(shared);
        if b.trans == Transpose::Trans {
            for l in 0..k {
                let begin = b.offset + col + j + (row + l) * b.ld;
                out[p..p + w].copy_from_slice(&source[begin..begin + w]);
                p += w;
            }
        } else if w == 1 {
            let begin = b.offset + row + (col + j) * b.ld;
            out[p..p + k].copy_from_slice(&source[begin..begin + k]);
            p += k;
        } else {
            let first = b.offset + row + (col + j) * b.ld;
            let second = first + b.ld;
            for (pair, (&x, &y)) in out[p..p + 2 * k].chunks_exact_mut(2).zip(
                source[first..first + k]
                    .iter()
                    .zip(&source[second..second + k]),
            ) {
                pair[0] = x;
                pair[1] = y;
            }
            p += 2 * k;
        }
        j += w;
    }
}
#[inline]
fn source_simd() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        return std::is_x86_feature_detected!("avx512f")
            && std::is_x86_feature_detected!("avx2")
            && std::is_x86_feature_detected!("fma");
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_wide_packed<const H: usize, const G: usize, const W: usize>(
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    co: usize,
    ldc: usize,
    begin: usize,
    end: usize,
    overwrite: bool,
) {
    // OpenBLAS❗❗: "vmovupd (%0),%%zmm1; vmovupd 64(%0),%%zmm2; addq $128,%0;"
    // OpenBLAS❗❗: "vbroadcastsd (%1),%%zmm3; vfmadd231pd %%zmm1,%%zmm3,%%zmm8; vfmadd231pd %%zmm2,%%zmm3,%%zmm9;"
    // OpenBLAS❗❗: "vfmadd213pd (%2),%%zmm0,%%zmm8; vmovupd %%zmm8,(%2); vfmadd213pd 64(%2),%%zmm0,%%zmm9; vmovupd %%zmm9,64(%2); addq $128,%2;"
    // The source's m16/m8 paths have one accumulator per output, with no
    // k-lane reduction. A column pair's triangular prefix/main/suffix steps
    // form one contiguous range, so their increasing-k recurrence is shared.
    use std::arch::x86_64::*;
    let mut sums = [[_mm512_setzero_pd(); G]; W];
    for l in begin..end {
        let mut av = [_mm512_setzero_pd(); G];
        for g in 0..G {
            // SAFETY: H=8*G (the source m16/m8 row tiles), l<k and the
            // packed panel contains H*k entries; CPU guard is source_simd.
            av[g] = unsafe { _mm512_loadu_pd(a.as_ptr().add(l * H + g * 8)) };
        }
        for col in 0..W {
            let bp = if W == 1 {
                l
            } else {
                col / 2 * 2 * k + l * 2 + col % 2
            };
            // SAFETY: paired B packing has W*k entries; the source loop
            // stays within this panel. Broadcast reads exactly one f64.
            let bv = unsafe { _mm512_broadcastsd_pd(_mm_load_sd(b.as_ptr().add(bp))) };
            for g in 0..G {
                sums[col][g] = _mm512_fmadd_pd(av[g], bv, sums[col][g]);
            }
        }
    }
    let alpha = _mm512_set1_pd(alpha);
    for col in 0..W {
        for g in 0..G {
            let ci = co + g * 8 + col * ldc;
            // SAFETY: the outer BLAS entry checked active output storage; source
            // row tiles contain eight complete rows and do not touch padding.
            unsafe {
                let result = if overwrite {
                    _mm512_mul_pd(alpha, sums[col][g])
                } else {
                    _mm512_fmadd_pd(alpha, sums[col][g], _mm512_loadu_pd(c.as_ptr().add(ci)))
                };
                _mm512_storeu_pd(c.as_mut_ptr().add(ci), result);
            }
        }
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_wide_packed_dispatch<const H: usize, const G: usize>(
    w: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    co: usize,
    ldc: usize,
    begin: usize,
    end: usize,
    overwrite: bool,
) {
    // Closed physical source tile sizes, not extensible business dispatch.
    macro_rules! run {
        ($w:expr) => {
            unsafe {
                native_wide_packed::<H, G, $w>(k, alpha, a, b, c, co, ldc, begin, end, overwrite)
            }
        };
    }
    match w {
        12 => run!(12),
        10 => run!(10),
        8 => run!(8),
        6 => run!(6),
        4 => run!(4),
        2 => run!(2),
        _ => run!(1),
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_wide_direct<
    const H: usize,
    const G: usize,
    const W: usize,
    const BT: bool,
    const MASKED: bool,
>(
    active: usize,
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: Matrix<'_>,
    ar: usize,
    bc: usize,
    beta: f64,
    c: &mut [f64],
    co: usize,
    row_stride: usize,
    col_stride: usize,
) {
    // OpenBLAS❗❗: #define LOAD_A_512(M, N) __m512d Aval##M = _mm512_loadu_pd(&A[lda * k + i + (M*8)])
    // OpenBLAS❗❗: #define MASK_LOAD_A_512(M, N) __m512d Aval##M = _mm512_maskz_loadu_pd(mask, &A[lda * k + i + (M*8)])
    // OpenBLAS❗❗: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // OpenBLAS❗❗: result##M##N = _mm512_mul_pd(result##M##N, alpha_512);
    // OpenBLAS❗❗: result##M##N = _mm512_fmadd_pd(tmp##M##N, beta_512, result##M##N);
    use std::arch::x86_64::*;
    let mut sums = [[_mm512_setzero_pd(); G]; W];
    let mask = if MASKED {
        ((1_u16 << active) - 1) as u8
    } else {
        255
    };
    let ap = a.separate.unwrap_or(c).as_ptr();
    let bp = b.separate.unwrap_or(c).as_ptr();
    for l in 0..k {
        let mut av = [_mm512_setzero_pd(); G];
        for g in 0..G {
            // SAFETY: the outer entry checked source rows and k. Complete
            // groups read eight values; masked residuals read only active rows.
            av[g] = unsafe {
                if MASKED {
                    _mm512_maskz_loadu_pd(mask, ap.add(a.offset + ar + l * a.ld))
                } else {
                    _mm512_loadu_pd(ap.add(a.offset + ar + g * 8 + l * a.ld))
                }
            };
        }
        for col in 0..W {
            let index = b.offset
                + if BT {
                    bc + col + l * b.ld
                } else {
                    l + (bc + col) * b.ld
                };
            // SAFETY: BT is the source's NN/NT dispatch, fixed before its k
            // loop; the checked input matrix contains this scalar broadcast.
            let bv = unsafe { _mm512_broadcastsd_pd(_mm_load_sd(bp.add(index))) };
            for g in 0..G {
                sums[col][g] = _mm512_fmadd_pd(av[g], bv, sums[col][g]);
            }
        }
    }
    let alpha = _mm512_set1_pd(alpha);
    let beta_v = _mm512_set1_pd(beta);
    for col in 0..W {
        for g in 0..G {
            let ci = co + g * 8 * row_stride + col * col_stride;
            let count = if MASKED { active } else { 8 };
            let mut tail = [0.; 8];
            let mut out = _mm512_mul_pd(alpha, sums[col][g]);
            if beta != 0. {
                let cv = if row_stride == 1 {
                    unsafe { _mm512_maskz_loadu_pd(mask, c.as_ptr().add(ci)) }
                } else {
                    for v in 0..count {
                        tail[v] = c[ci + v * row_stride];
                    }
                    unsafe { _mm512_loadu_pd(tail.as_ptr()) }
                };
                out = _mm512_fmadd_pd(cv, beta_v, out);
            }
            // SAFETY: output masks and strided staging cover only checked active
            // rows. The TT/NT-tail projection swaps independent output axes only.
            unsafe {
                if row_stride == 1 {
                    _mm512_mask_storeu_pd(c.as_mut_ptr().add(ci), mask, out);
                } else {
                    _mm512_storeu_pd(tail.as_mut_ptr(), out);
                    for v in 0..count {
                        c[ci + v * row_stride] = tail[v];
                    }
                }
            }
        }
    }
    let _ = H;
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_wide_direct_dispatch<const H: usize, const G: usize, const MASKED: bool>(
    active: usize,
    w: usize,
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: Matrix<'_>,
    ar: usize,
    bc: usize,
    beta: f64,
    c: &mut [f64],
    co: usize,
    row_stride: usize,
    col_stride: usize,
) {
    macro_rules! run {
        ($w:expr,$bt:expr) => {
            unsafe {
                native_wide_direct::<H, G, $w, $bt, MASKED>(
                    active, k, alpha, a, b, ar, bc, beta, c, co, row_stride, col_stride,
                )
            }
        };
    }
    if b.trans == Transpose::Trans {
        match w {
            6 => run!(6, true),
            4 => run!(4, true),
            2 => run!(2, true),
            _ => run!(1, true),
        }
    } else {
        match w {
            6 => run!(6, false),
            4 => run!(4, false),
            2 => run!(2, false),
            _ => run!(1, false),
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_tail_packed<const H: usize, const L: usize, const W: usize>(
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    co: usize,
    ldc: usize,
    prefix_begin: usize,
    main_begin: usize,
    main_end: usize,
    suffix_end: usize,
    overwrite: bool,
) {
    // OpenBLAS❗❗: "vmovddup (%0),%%zmm1; vmovddup 8(%0),%%zmm2; addq $64,%0;"
    // OpenBLAS❗❗: "vfmadd231pd %%zmm1,%%zmm3,%%zmm"#c1_no"; vfmadd231pd %%zmm2,%%zmm3,%%zmm"#c2_no";"
    // OpenBLAS❗❗: "vextractf64x4 $1,%%zmm"#c1_no",%%ymm30; vaddpd %%ymm"#c1_no",%%ymm30,%%ymm"#c1_no";"
    // OpenBLAS❗❗: "vaddpd %%xmm20,%%xmm21,%%xmm20; vaddpd %%xmm22,%%xmm23,%%xmm22; vaddpd %%xmm20,%%xmm22,%%xmm"#c1_no";"
    // Pack (k-lane,row,column) into the same eight physical FMA lanes. Source
    // m4n12 uses two k lanes; m2n12/m1n12 use four. No arithmetic reduction is
    // introduced for the other source kernels.
    use std::arch::x86_64::*;
    macro_rules! ai {
        ($s:expr) => {
            ($s / (H * W) * H + ($s / W) % H) as i64
        };
    }
    macro_rules! bi {
        ($s:expr) => {
            ($s / (H * W) * W + $s % W) as i64
        };
    }
    let avars = [
        _mm512_set_epi64(
            ai!(7),
            ai!(6),
            ai!(5),
            ai!(4),
            ai!(3),
            ai!(2),
            ai!(1),
            ai!(0),
        ),
        _mm512_set_epi64(
            ai!(15),
            ai!(14),
            ai!(13),
            ai!(12),
            ai!(11),
            ai!(10),
            ai!(9),
            ai!(8),
        ),
    ];
    let bvars = [
        _mm512_set_epi64(
            bi!(7),
            bi!(6),
            bi!(5),
            bi!(4),
            bi!(3),
            bi!(2),
            bi!(1),
            bi!(0),
        ),
        _mm512_set_epi64(
            bi!(15),
            bi!(14),
            bi!(13),
            bi!(12),
            bi!(11),
            bi!(10),
            bi!(9),
            bi!(8),
        ),
    ];
    let count = (H * W * L).div_ceil(8);
    let mut sums = [_mm512_setzero_pd(); 2];
    let cells = H * W;
    let cells_mask = ((1_u16 << cells) - 1) as u8;
    macro_rules! scalar_step {
        ($l:expr) => {{
            let l = $l;
            // SAFETY: prefix/tail k stays within packed storage; mask loads read
            // exactly H and W valid elements, and mask3 preserves other k lanes.
            unsafe {
                let av = _mm512_maskz_loadu_pd(((1_u16 << H) - 1) as u8, a.as_ptr().add(l * H));
                let bv = _mm512_maskz_loadu_pd(((1_u16 << W) - 1) as u8, b.as_ptr().add(l * W));
                sums[0] = _mm512_mask3_fmadd_pd(
                    _mm512_permutexvar_pd(avars[0], av),
                    _mm512_permutexvar_pd(bvars[0], bv),
                    sums[0],
                    cells_mask,
                );
            }
        }};
    }
    for l in prefix_begin..main_begin {
        scalar_step!(l);
    }
    let bulk = main_begin + (main_end - main_begin) / L * L;
    for l in (main_begin..bulk).step_by(L) {
        // SAFETY: a full source k block contains H*L and W*L elements; both
        // are at most eight. The driver validates packed dimensions/ranges.
        unsafe {
            let av = _mm512_maskz_loadu_pd(((1_u16 << (H * L)) - 1) as u8, a.as_ptr().add(l * H));
            let bv = _mm512_maskz_loadu_pd(((1_u16 << (W * L)) - 1) as u8, b.as_ptr().add(l * W));
            for v in 0..count {
                sums[v] = _mm512_fmadd_pd(
                    _mm512_permutexvar_pd(avars[v], av),
                    _mm512_permutexvar_pd(bvars[v], bv),
                    sums[v],
                );
            }
        }
    }
    if L > 1 && bulk > main_begin {
        let mut flat = [0_f64; 16];
        // SAFETY: each local store has eight valid f64 slots.
        unsafe {
            _mm512_storeu_pd(flat.as_mut_ptr(), sums[0]);
            _mm512_storeu_pd(flat.as_mut_ptr().add(8), sums[1]);
        }
        for cell in 0..cells {
            flat[cell] = if L == 2 {
                flat[cell] + flat[cell + cells]
            } else {
                (flat[cell] + flat[cell + cells])
                    + (flat[cell + cells * 2] + flat[cell + cells * 3])
            };
        }
        // SAFETY: local first eight lanes contain the reduced source outputs.
        sums[0] = unsafe { _mm512_loadu_pd(flat.as_ptr()) };
    }
    for l in bulk..suffix_end {
        scalar_step!(l);
    }
    let mut flat = [0_f64; 8];
    unsafe {
        _mm512_storeu_pd(flat.as_mut_ptr(), sums[0]);
    }
    for row in 0..H {
        for col in 0..W {
            let ci = co + row + col * ldc;
            let sum = flat[row * W + col];
            c[ci] = if overwrite {
                alpha * sum
            } else {
                alpha.mul_add(sum, c[ci])
            };
        }
    }
    let _ = k;
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_tail_dispatch<const H: usize, const L: usize>(
    bw: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    co: usize,
    ldc: usize,
    prefix_begin: usize,
    main_begin: usize,
    main_end: usize,
    suffix_end: usize,
    overwrite: bool,
) {
    // SAFETY: paired source_simd guard and source row/column tile dimensions.
    unsafe {
        if bw == 2 {
            native_tail_packed::<H, L, 2>(
                k,
                alpha,
                a,
                b,
                c,
                co,
                ldc,
                prefix_begin,
                main_begin,
                main_end,
                suffix_end,
                overwrite,
            );
        } else {
            native_tail_packed::<H, L, 1>(
                k,
                alpha,
                a,
                b,
                c,
                co,
                ldc,
                prefix_begin,
                main_begin,
                main_end,
                suffix_end,
                overwrite,
            );
        }
    }
}

#[derive(Clone, Copy)]
struct TriKernel {
    side: Side,
    upper: bool,
    diagonal: isize,
}
fn packed_kernel(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    offset: usize,
    ldc: usize,
) {
    packed_product(m, n, k, alpha, a, b, c, offset, ldc, None)
}
fn packed_product(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    c: &mut [f64],
    offset: usize,
    ldc: usize,
    tri: Option<TriKernel>,
) {
    // OpenBLAS❗✔️: /* %0 = a_ptr, %1 = b_ptr, %2 = c_ptr, %3 = c_tmp, %4 = ldc(bytes), %5 = k_counter, %6 = b_pref */
    // OpenBLAS❗✔️: /* r10 = tmp, r11 = m_counter, r12 = size_of_1_tile_in_b, r13 = k, r14 = b_head, r15 = %1+3*r12 */
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if (defined (LEFT) && !defined(TRANSA)) || (!defined (LEFT) && defined(TRANSA))
    // OpenBLAS❗✔️:   #define BACKWARDS 1
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define BACKWARDS 0
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define GEMM_SET_PB "movq %%r14,%1; leaq (%%r14,%%r12,2),%%r15; addq %%r12,%%r15;"
    // OpenBLAS❗✔️: #define set_p_copy1(ptr) "sarq $1,%%r12; addq %%r12,"#ptr"; salq $1,%%r12; salq $3,%%r13; subq %%r13,"#ptr"; sarq $3,%%r13;"
    // OpenBLAS❗✔️: #define set_p_copy2(ptr) "addq %%r12,"#ptr"; salq $4,%%r13; subq %%r13,"#ptr"; sarq $4,%%r13;"
    // OpenBLAS❗✔️: #define set_p_copy4(ptr) "leaq ("#ptr",%%r12,2),"#ptr"; salq $5,%%r13; subq %%r13,"#ptr"; sarq $5,%%r13;"
    // OpenBLAS❗✔️: #define set_p_copy8(ptr) "leaq ("#ptr",%%r12,4),"#ptr"; salq $6,%%r13; subq %%r13,"#ptr"; sarq $6,%%r13;"
    // OpenBLAS❗✔️: #define set_p_copy16(ptr) "leaq ("#ptr",%%r12,8),"#ptr"; salq $7,%%r13; subq %%r13,"#ptr"; sarq $7,%%r13;"
    // OpenBLAS❗✔️: #define set_p_b_dim1(ptr) set_p_copy1(ptr)
    // OpenBLAS❗✔️: #define set_p_b_dim2(ptr) set_p_copy2(ptr)
    // OpenBLAS❗✔️: #define set_p_b_dim4(ptr) set_p_copy2(ptr)
    // OpenBLAS❗✔️: #define set_p_b_dim6(ptr) set_p_copy2(ptr)
    // OpenBLAS❗✔️: #define set_p_b_dim8(ptr) set_p_copy2(ptr)
    // OpenBLAS❗✔️: #define set_p_b_dim10(ptr) set_p_copy2(ptr)
    // OpenBLAS❗✔️: #define set_p_b_dim12(ptr) set_p_copy2(ptr)
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #if BACKWARDS == 1
    // OpenBLAS❗✔️:     #define INIT_set_papb(mdim,ndim) GEMM_SET_PB set_p_copy##mdim(%0) set_p_b_dim##ndim(%1) set_p_b_dim##ndim(%%r15)
    // OpenBLAS❗✔️:     #define SAVE_set_pa(mdim) ""
    // OpenBLAS❗✔️:   #else
    // OpenBLAS❗✔️:     #define INIT_set_papb(mdim,ndim) GEMM_SET_PB
    // OpenBLAS❗✔️:     #define SAVE_set_pa(mdim) set_p_copy##mdim(%0)
    // OpenBLAS❗✔️:   #endif
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define INIT_set_papb(mdim,ndim) GEMM_SET_PB
    // OpenBLAS❗✔️:   #define SAVE_set_pa(mdim) ""
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT)
    // OpenBLAS❗✔️:   #if BACKWARDS == 1
    // OpenBLAS❗✔️:     #define HEAD_SET_OFF(ndim) {}
    // OpenBLAS❗✔️:     #define TAIL_SET_OFF(ndim) {off += ndim;}
    // OpenBLAS❗✔️:     #define kernel_kstart_n4(mdim,updk) KERNEL_k1m##mdim##n2 KERNEL_k1m##mdim##n2 "addq $32,%%r15; "#updk" $2,%5;"
    // OpenBLAS❗✔️:     #define kernel_kstart_n6(mdim,updk) kernel_kstart_n4(mdim,updk) KERNEL_k1m##mdim##n4 KERNEL_k1m##mdim##n4 "addq $32,%%r15; "#updk" $2,%5;"
    // OpenBLAS❗✔️:     #define kernel_kstart_n8(mdim,updk) kernel_kstart_n6(mdim,updk) KERNEL_k1m##mdim##n6 KERNEL_k1m##mdim##n6 "addq $32,%%r15; "#updk" $2,%5;"
    // OpenBLAS❗✔️:     #define kernel_kstart_n10(mdim,updk) kernel_kstart_n8(mdim,updk) KERNEL_k1m##mdim##n8 KERNEL_k1m##mdim##n8 #updk" $2,%5;"
    // OpenBLAS❗✔️:     #define kernel_kstart_n12(mdim,updk) kernel_kstart_n10(mdim,updk) KERNEL_k1m##mdim##n10 KERNEL_k1m##mdim##n10 #updk" $2,%5;"
    // OpenBLAS❗✔️:     #define kernel_kend_n4(mdim) ""
    // OpenBLAS❗✔️:     #define kernel_kend_n6(mdim) ""
    // OpenBLAS❗✔️:     #define kernel_kend_n8(mdim) ""
    // OpenBLAS❗✔️:     #define kernel_kend_n10(mdim) ""
    // OpenBLAS❗✔️:     #define kernel_kend_n12(mdim) ""
    // OpenBLAS❗✔️:   #else
    // OpenBLAS❗✔️:     #define HEAD_SET_OFF(ndim) {off += (ndim > 2 ? 2 : ndim);}
    // OpenBLAS❗✔️:     #define TAIL_SET_OFF(ndim) {off += (ndim > 2 ? (ndim-2) : 0);}
    // OpenBLAS❗✔️:     #define kernel_kstart_n4(mdim,updk) ""
    // OpenBLAS❗✔️:     #define kernel_kstart_n6(mdim,updk) ""
    // OpenBLAS❗✔️:     #define kernel_kstart_n8(mdim,updk) ""
    // OpenBLAS❗✔️:     #define kernel_kstart_n10(mdim,updk) ""
    // OpenBLAS❗✔️:     #define kernel_kstart_n12(mdim,updk) ""
    // OpenBLAS❗✔️:     #define kernel_kend_n4(mdim) "xorq %3,%3;"\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(0)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(16)
    // OpenBLAS❗✔️:     #define kernel_kend_n6(mdim) "xorq %3,%3;"\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(0) acc_kend_nc3_k1m##mdim(0)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(16) acc_kend_nc3_k1m##mdim(16)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(32)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(48)
    // OpenBLAS❗✔️:     #define kernel_kend_n8(mdim) "xorq %3,%3;"\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(0) acc_kend_nc3_k1m##mdim(0) acc_kend_nc4_k1m##mdim(0)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(16) acc_kend_nc3_k1m##mdim(16) acc_kend_nc4_k1m##mdim(16)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(32) acc_kend_nc4_k1m##mdim(32)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(48) acc_kend_nc4_k1m##mdim(48)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc4_k1m##mdim(64)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc4_k1m##mdim(80)
    // OpenBLAS❗✔️:     #define kernel_kend_n10(mdim) "xorq %3,%3;"\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(0) acc_kend_nc3_k1m##mdim(0) acc_kend_nc4_k1m##mdim(0) acc_kend_nc5_k1m##mdim(0)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(16) acc_kend_nc3_k1m##mdim(16) acc_kend_nc4_k1m##mdim(16) acc_kend_nc5_k1m##mdim(16)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(32) acc_kend_nc4_k1m##mdim(32) acc_kend_nc5_k1m##mdim(32)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(48) acc_kend_nc4_k1m##mdim(48) acc_kend_nc5_k1m##mdim(48)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc4_k1m##mdim(64) acc_kend_nc5_k1m##mdim(64)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc4_k1m##mdim(80) acc_kend_nc5_k1m##mdim(80)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc5_k1m##mdim(96)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc5_k1m##mdim(112)
    // OpenBLAS❗✔️:     #define kernel_kend_n12(mdim) "xorq %3,%3;"\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(0) acc_kend_nc3_k1m##mdim(0) acc_kend_nc4_k1m##mdim(0) acc_kend_nc5_k1m##mdim(0) acc_kend_nc6_k1m##mdim(0)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc2_k1m##mdim(16) acc_kend_nc3_k1m##mdim(16) acc_kend_nc4_k1m##mdim(16) acc_kend_nc5_k1m##mdim(16) acc_kend_nc6_k1m##mdim(16)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(32) acc_kend_nc4_k1m##mdim(32) acc_kend_nc5_k1m##mdim(32) acc_kend_nc6_k1m##mdim(32)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc3_k1m##mdim(48) acc_kend_nc4_k1m##mdim(48) acc_kend_nc5_k1m##mdim(48) acc_kend_nc6_k1m##mdim(48)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc4_k1m##mdim(64) acc_kend_nc5_k1m##mdim(64) acc_kend_nc6_k1m##mdim(64)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc4_k1m##mdim(80) acc_kend_nc5_k1m##mdim(80) acc_kend_nc6_k1m##mdim(80)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc5_k1m##mdim(96) acc_kend_nc6_k1m##mdim(96)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc5_k1m##mdim(112) acc_kend_nc6_k1m##mdim(112)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc6_k1m##mdim(128)\
    // OpenBLAS❗✔️:       loada_kend_k1m##mdim acc_kend_nc6_k1m##mdim(144)
    // OpenBLAS❗✔️:   #endif
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define HEAD_SET_OFF(ndim) {}
    // OpenBLAS❗✔️:   #define TAIL_SET_OFF(ndim) {}
    // OpenBLAS❗✔️:   #define kernel_kstart_n4(mdim,updk) ""
    // OpenBLAS❗✔️:   #define kernel_kstart_n6(mdim,updk) ""
    // OpenBLAS❗✔️:   #define kernel_kstart_n8(mdim,updk) ""
    // OpenBLAS❗✔️:   #define kernel_kstart_n10(mdim,updk) ""
    // OpenBLAS❗✔️:   #define kernel_kstart_n12(mdim,updk) ""
    // OpenBLAS❗✔️:   #define kernel_kend_n4(mdim) ""
    // OpenBLAS❗✔️:   #define kernel_kend_n6(mdim) ""
    // OpenBLAS❗✔️:   #define kernel_kend_n8(mdim) ""
    // OpenBLAS❗✔️:   #define kernel_kend_n10(mdim) ""
    // OpenBLAS❗✔️:   #define kernel_kend_n12(mdim) ""
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define kernel_kstart_n1(mdim,updk) ""
    // OpenBLAS❗✔️: #define kernel_kstart_n2(mdim,updk) ""
    // OpenBLAS❗✔️: #define kernel_kend_n1(mdim) ""
    // OpenBLAS❗✔️: #define kernel_kend_n2(mdim) ""
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #if BACKWARDS == 1
    // OpenBLAS❗✔️:     #define INITASM_SET_K "movq %10,%%r13; subq %9,%%r13;"
    // OpenBLAS❗✔️:   #else
    // OpenBLAS❗✔️:     #define INITASM_SET_K "movq %9,%%r13;"
    // OpenBLAS❗✔️:   #endif
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define INITASM_SET_K "movq %10,%%r13;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && defined(LEFT)
    // OpenBLAS❗✔️:   #if BACKWARDS==1
    // OpenBLAS❗✔️:     #define init_update_k(mdim) ""
    // OpenBLAS❗✔️:     #define save_update_k(mdim) "subq $"#mdim",%%r13;"
    // OpenBLAS❗✔️:   #else
    // OpenBLAS❗✔️:     #define init_update_k(mdim) "addq $"#mdim",%%r13;"
    // OpenBLAS❗✔️:     #define save_update_k(mdim) ""
    // OpenBLAS❗✔️:   #endif
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define init_update_k(mdim) ""
    // OpenBLAS❗✔️:   #define save_update_k(mdim) ""
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_h_k1m16n1 \
    // OpenBLAS❗✔️:   "vmovupd (%0),%%zmm1; vmovupd 64(%0),%%zmm2; addq $128,%0;"\
    // OpenBLAS❗✔️:   "vbroadcastsd (%1),%%zmm3; vfmadd231pd %%zmm1,%%zmm3,%%zmm8; vfmadd231pd %%zmm2,%%zmm3,%%zmm9;"
    // OpenBLAS❗✔️: #define KERNEL_k1m16n1 KERNEL_h_k1m16n1 "addq $8,%1;"
    // OpenBLAS❗✔️: #ifdef BROADCAST_KERNEL
    // OpenBLAS❗✔️:  #define KERNEL_h_k1m16n2 KERNEL_h_k1m16n1\
    // OpenBLAS❗✔️:   "vbroadcastsd 8(%1),%%zmm4; vfmadd231pd %%zmm1,%%zmm4,%%zmm10; vfmadd231pd %%zmm2,%%zmm4,%%zmm11;"
    // OpenBLAS❗✔️:  #define unit_acc_gen_m16n2(c1_no,c2_no,c3_no,c4_no,boff1,...)\
    // OpenBLAS❗✔️:   "vbroadcastsd "#boff1"("#__VA_ARGS__"),%%zmm3; vfmadd231pd %%zmm1,%%zmm3,%%zmm"#c1_no"; vfmadd231pd %%zmm2,%%zmm3,%%zmm"#c2_no";"\
    // OpenBLAS❗✔️:   "vbroadcastsd "#boff1"+8("#__VA_ARGS__"),%%zmm4; vfmadd231pd %%zmm1,%%zmm4,%%zmm"#c3_no"; vfmadd231pd %%zmm2,%%zmm4,%%zmm"#c4_no";"
    // OpenBLAS❗✔️:  #define unit_acc_m16n2(c1_no,c2_no,c3_no,c4_no,...) unit_acc_gen_m16n2(c1_no,c2_no,c3_no,c4_no,0,__VA_ARGS__)
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:  #define unit_acc_gen_m16n2(c1_no,c2_no,c3_no,c4_no,boff1,...)\
    // OpenBLAS❗✔️:   "vbroadcastf32x4 "#boff1"("#__VA_ARGS__"),%%zmm5; vfmadd231pd %%zmm1,%%zmm5,%%zmm"#c1_no"; vfmadd231pd %%zmm2,%%zmm5,%%zmm"#c2_no";"\
    // OpenBLAS❗✔️:   "vfmadd231pd %%zmm3,%%zmm5,%%zmm"#c3_no"; vfmadd231pd %%zmm4,%%zmm5,%%zmm"#c4_no";"
    // OpenBLAS❗✔️:  #define unit_acc_m16n2(c1_no,c2_no,c3_no,c4_no,...) unit_acc_gen_m16n2(c1_no,c2_no,c3_no,c4_no,0,__VA_ARGS__)
    // OpenBLAS❗✔️:  #define KERNEL_h_k1m16n2 \
    // OpenBLAS❗✔️:   "vmovddup (%0),%%zmm1; vmovddup 8(%0),%%zmm2; vmovddup 64(%0),%%zmm3; vmovddup 72(%0),%%zmm4; addq $128,%0;"\
    // OpenBLAS❗✔️:   unit_acc_m16n2(8,9,10,11,%1)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define KERNEL_k1m16n2 KERNEL_h_k1m16n2 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m16n4 KERNEL_h_k1m16n2 "prefetcht0 384(%0);" unit_acc_m16n2(12,13,14,15,%1,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m16n4 KERNEL_h_k1m16n4 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_k1m16n6 KERNEL_h_k1m16n4 unit_acc_m16n2(16,17,18,19,%1,%%r12,2) "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m16n8 KERNEL_k1m16n6 "prefetcht0 448(%0);" unit_acc_m16n2(20,21,22,23,%%r15)
    // OpenBLAS❗✔️: #define KERNEL_k1m16n8 KERNEL_h_k1m16n8 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m16n10 KERNEL_h_k1m16n8 unit_acc_m16n2(24,25,26,27,%%r15,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m16n10 KERNEL_h_k1m16n10 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m16n12 KERNEL_h_k1m16n10 unit_acc_m16n2(28,29,30,31,%%r15,%%r12,2)
    // OpenBLAS❗✔️: #define KERNEL_k1m16n12 KERNEL_h_k1m16n12 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT) && (BACKWARDS == 0)
    // OpenBLAS❗✔️:  #ifdef BROADCAST_KERNEL
    // OpenBLAS❗✔️:   #define loada_kend_k1m16 "vmovupd (%0,%3,1),%%zmm1; vmovupd 64(%0,%3,1),%%zmm2; addq $128,%3;"
    // OpenBLAS❗✔️:  #else
    // OpenBLAS❗✔️:   #define loada_kend_k1m16 "vmovddup (%0,%3,1),%%zmm1; vmovddup 8(%0,%3,1),%%zmm2; vmovddup 64(%0,%3,1),%%zmm3; vmovddup 72(%0,%3,1),%%zmm4; addq $128,%3;"
    // OpenBLAS❗✔️:  #endif
    // OpenBLAS❗✔️:  #define acc_kend_nc2_k1m16(boff1) unit_acc_gen_m16n2(12,13,14,15,boff1,%1,%%r12,1)
    // OpenBLAS❗✔️:  #define acc_kend_nc3_k1m16(boff1) unit_acc_gen_m16n2(16,17,18,19,boff1,%1,%%r12,2)
    // OpenBLAS❗✔️:  #define acc_kend_nc4_k1m16(boff1) unit_acc_gen_m16n2(20,21,22,23,boff1,%%r15)
    // OpenBLAS❗✔️:  #define acc_kend_nc5_k1m16(boff1) unit_acc_gen_m16n2(24,25,26,27,boff1,%%r15,%%r12,1)
    // OpenBLAS❗✔️:  #define acc_kend_nc6_k1m16(boff1) unit_acc_gen_m16n2(28,29,30,31,boff1,%%r15,%%r12,2)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define save_init_m16 "movq %2,%3; addq $128,%2;"
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #define SAVE_m16n1 "vmulpd %%zmm8,%%zmm0,%%zmm8; vmovupd %%zmm8,(%2); vmulpd %%zmm9,%%zmm0,%%zmm9; vmovupd %%zmm9,64(%2); addq $128,%2;"
    // OpenBLAS❗✔️:  #ifdef BROADCAST_KERNEL
    // OpenBLAS❗✔️:   #define unit_save_m16n2(c1_no,c2_no,c3_no,c4_no)\
    // OpenBLAS❗✔️:     "vmulpd %%zmm"#c1_no",%%zmm0,%%zmm"#c1_no"; vmovupd %%zmm"#c1_no",(%3); vmulpd %%zmm"#c2_no",%%zmm0,%%zmm"#c2_no"; vmovupd %%zmm"#c2_no",64(%3);"\
    // OpenBLAS❗✔️:     "vmulpd %%zmm"#c3_no",%%zmm0,%%zmm"#c3_no"; vmovupd %%zmm"#c3_no",(%3,%4,1); vmulpd %%zmm"#c4_no",%%zmm0,%%zmm"#c4_no"; vmovupd %%zmm"#c4_no",64(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️:  #else
    // OpenBLAS❗✔️:   #define unit_save_m16n2(c1_no,c2_no,c3_no,c4_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm1; vunpcklpd %%zmm"#c4_no",%%zmm"#c3_no",%%zmm2; vunpckhpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm3; vunpckhpd %%zmm"#c4_no",%%zmm"#c3_no",%%zmm4;"\
    // OpenBLAS❗✔️:     "vmulpd %%zmm1,%%zmm0,%%zmm1; vmovupd %%zmm1,(%3); vmulpd %%zmm2,%%zmm0,%%zmm2; vmovupd %%zmm2,64(%3);"\
    // OpenBLAS❗✔️:     "vmulpd %%zmm3,%%zmm0,%%zmm3; vmovupd %%zmm3,(%3,%4,1); vmulpd %%zmm4,%%zmm0,%%zmm4; vmovupd %%zmm4,64(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️:  #endif
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define SAVE_m16n1 "vfmadd213pd (%2),%%zmm0,%%zmm8; vmovupd %%zmm8,(%2); vfmadd213pd 64(%2),%%zmm0,%%zmm9; vmovupd %%zmm9,64(%2); addq $128,%2;"
    // OpenBLAS❗✔️:  #ifdef BROADCAST_KERNEL
    // OpenBLAS❗✔️:   #define unit_save_m16n2(c1_no,c2_no,c3_no,c4_no)\
    // OpenBLAS❗✔️:     "vfmadd213pd (%3),%%zmm0,%%zmm"#c1_no"; vmovupd %%zmm"#c1_no",(%3); vfmadd213pd 64(%3),%%zmm0,%%zmm"#c2_no"; vmovupd %%zmm"#c2_no",64(%3);"\
    // OpenBLAS❗✔️:     "vfmadd213pd (%3,%4,1),%%zmm0,%%zmm"#c3_no"; vmovupd %%zmm"#c3_no",(%3,%4,1); vfmadd213pd 64(%3,%4,1),%%zmm0,%%zmm"#c4_no"; vmovupd %%zmm"#c4_no",64(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️:  #else
    // OpenBLAS❗✔️:   #define unit_save_m16n2(c1_no,c2_no,c3_no,c4_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm1; vunpcklpd %%zmm"#c4_no",%%zmm"#c3_no",%%zmm2; vunpckhpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm3; vunpckhpd %%zmm"#c4_no",%%zmm"#c3_no",%%zmm4;"\
    // OpenBLAS❗✔️:     "vfmadd213pd (%3),%%zmm0,%%zmm1; vmovupd %%zmm1,(%3); vfmadd213pd 64(%3),%%zmm0,%%zmm2; vmovupd %%zmm2,64(%3);"\
    // OpenBLAS❗✔️:     "vfmadd213pd (%3,%4,1),%%zmm0,%%zmm3; vmovupd %%zmm3,(%3,%4,1); vfmadd213pd 64(%3,%4,1),%%zmm0,%%zmm4; vmovupd %%zmm4,64(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️:  #endif
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define SAVE_m16n2 save_init_m16 unit_save_m16n2(8,9,10,11)
    // OpenBLAS❗✔️: #define SAVE_m16n4 SAVE_m16n2 unit_save_m16n2(12,13,14,15)
    // OpenBLAS❗✔️: #define SAVE_m16n6 SAVE_m16n4 unit_save_m16n2(16,17,18,19)
    // OpenBLAS❗✔️: #define SAVE_m16n8 SAVE_m16n6 unit_save_m16n2(20,21,22,23)
    // OpenBLAS❗✔️: #define SAVE_m16n10 SAVE_m16n8 unit_save_m16n2(24,25,26,27)
    // OpenBLAS❗✔️: #define SAVE_m16n12 SAVE_m16n10 unit_save_m16n2(28,29,30,31)
    // OpenBLAS❗✔️: #define unit_init_2zmm(c1_no,c2_no) "vpxorq %%zmm"#c1_no",%%zmm"#c1_no",%%zmm"#c1_no"; vpxorq %%zmm"#c2_no",%%zmm"#c2_no",%%zmm"#c2_no";"
    // OpenBLAS❗✔️: #define unit_init_4zmm(c1_no,c2_no,c3_no,c4_no) unit_init_2zmm(c1_no,c2_no) unit_init_2zmm(c3_no,c4_no)
    // OpenBLAS❗✔️: #define INIT_m16n1 unit_init_2zmm(8,9)
    // OpenBLAS❗✔️: #define INIT_m16n2 unit_init_4zmm(8,9,10,11)
    // OpenBLAS❗✔️: #define INIT_m16n4 INIT_m16n2 unit_init_4zmm(12,13,14,15)
    // OpenBLAS❗✔️: #define INIT_m16n6 INIT_m16n4 unit_init_4zmm(16,17,18,19)
    // OpenBLAS❗✔️: #define INIT_m16n8 INIT_m16n6 unit_init_4zmm(20,21,22,23)
    // OpenBLAS❗✔️: #define INIT_m16n10 INIT_m16n8 unit_init_4zmm(24,25,26,27)
    // OpenBLAS❗✔️: #define INIT_m16n12 INIT_m16n10 unit_init_4zmm(28,29,30,31)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_k1m8n1 \
    // OpenBLAS❗✔️:   "vbroadcastsd (%1),%%zmm1; addq $8,%1;"\
    // OpenBLAS❗✔️:   "vfmadd231pd (%0),%%zmm1,%%zmm8; addq $64,%0;"
    // OpenBLAS❗✔️: #define unit_acc_gen_m8n2(c1_no,c2_no,boff,...)\
    // OpenBLAS❗✔️:   "vbroadcastf32x4 "#boff"("#__VA_ARGS__"),%%zmm3; vfmadd231pd %%zmm1,%%zmm3,%%zmm"#c1_no"; vfmadd231pd %%zmm2,%%zmm3,%%zmm"#c2_no";"
    // OpenBLAS❗✔️: #define unit_acc_m8n2(c1_no,c2_no,...) unit_acc_gen_m8n2(c1_no,c2_no,0,__VA_ARGS__)
    // OpenBLAS❗✔️: #define KERNEL_h_k1m8n2 \
    // OpenBLAS❗✔️:   "vmovddup (%0),%%zmm1; vmovddup 8(%0),%%zmm2; addq $64,%0;" unit_acc_m8n2(8,9,%1)
    // OpenBLAS❗✔️: #define KERNEL_k1m8n2 KERNEL_h_k1m8n2 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m8n4 KERNEL_h_k1m8n2 unit_acc_m8n2(10,11,%1,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m8n4 KERNEL_h_k1m8n4 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_k1m8n6 KERNEL_h_k1m8n4 unit_acc_m8n2(12,13,%1,%%r12,2) "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m8n8 KERNEL_k1m8n6 unit_acc_m8n2(14,15,%%r15)
    // OpenBLAS❗✔️: #define KERNEL_k1m8n8 KERNEL_h_k1m8n8 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m8n10 KERNEL_h_k1m8n8 unit_acc_m8n2(16,17,%%r15,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m8n10 KERNEL_h_k1m8n10 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m8n12 KERNEL_h_k1m8n10 unit_acc_m8n2(18,19,%%r15,%%r12,2)
    // OpenBLAS❗✔️: #define KERNEL_k1m8n12 KERNEL_h_k1m8n12 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT) && (BACKWARDS == 0)
    // OpenBLAS❗✔️:   #define loada_kend_k1m8 "vmovddup (%0,%3,1),%%zmm1; vmovddup 8(%0,%3,1),%%zmm2; addq $64,%3;"
    // OpenBLAS❗✔️:   #define acc_kend_nc2_k1m8(boff1) unit_acc_gen_m8n2(10,11,boff1,%1,%%r12,1)
    // OpenBLAS❗✔️:   #define acc_kend_nc3_k1m8(boff1) unit_acc_gen_m8n2(12,13,boff1,%1,%%r12,2)
    // OpenBLAS❗✔️:   #define acc_kend_nc4_k1m8(boff1) unit_acc_gen_m8n2(14,15,boff1,%%r15)
    // OpenBLAS❗✔️:   #define acc_kend_nc5_k1m8(boff1) unit_acc_gen_m8n2(16,17,boff1,%%r15,%%r12,1)
    // OpenBLAS❗✔️:   #define acc_kend_nc6_k1m8(boff1) unit_acc_gen_m8n2(18,19,boff1,%%r15,%%r12,2)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define save_init_m8 "movq %2,%3; addq $64,%2;"
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #define SAVE_m8n1 "vmulpd %%zmm8,%%zmm0,%%zmm8; vmovupd %%zmm8,(%2); addq $64,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m8n2(c1_no,c2_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm1; vmulpd %%zmm1,%%zmm0,%%zmm1; vmovupd %%zmm1,(%3);"\
    // OpenBLAS❗✔️:     "vunpckhpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm2; vmulpd %%zmm2,%%zmm0,%%zmm2; vmovupd %%zmm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define SAVE_m8n1 "vfmadd213pd (%2),%%zmm0,%%zmm8; vmovupd %%zmm8,(%2); addq $64,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m8n2(c1_no,c2_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm1; vfmadd213pd (%3),%%zmm0,%%zmm1; vmovupd %%zmm1,(%3);"\
    // OpenBLAS❗✔️:     "vunpckhpd %%zmm"#c2_no",%%zmm"#c1_no",%%zmm2; vfmadd213pd (%3,%4,1),%%zmm0,%%zmm2; vmovupd %%zmm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define SAVE_m8n2 save_init_m8 unit_save_m8n2(8,9)
    // OpenBLAS❗✔️: #define SAVE_m8n4 SAVE_m8n2 unit_save_m8n2(10,11)
    // OpenBLAS❗✔️: #define SAVE_m8n6 SAVE_m8n4 unit_save_m8n2(12,13)
    // OpenBLAS❗✔️: #define SAVE_m8n8 SAVE_m8n6 unit_save_m8n2(14,15)
    // OpenBLAS❗✔️: #define SAVE_m8n10 SAVE_m8n8 unit_save_m8n2(16,17)
    // OpenBLAS❗✔️: #define SAVE_m8n12 SAVE_m8n10 unit_save_m8n2(18,19)
    // OpenBLAS❗✔️: #define INIT_m8n1 "vpxorq %%zmm8,%%zmm8,%%zmm8;"
    // OpenBLAS❗✔️: #define INIT_m8n2 unit_init_2zmm(8,9)
    // OpenBLAS❗✔️: #define INIT_m8n4 INIT_m8n2 unit_init_2zmm(10,11)
    // OpenBLAS❗✔️: #define INIT_m8n6 INIT_m8n4 unit_init_2zmm(12,13)
    // OpenBLAS❗✔️: #define INIT_m8n8 INIT_m8n6 unit_init_2zmm(14,15)
    // OpenBLAS❗✔️: #define INIT_m8n10 INIT_m8n8 unit_init_2zmm(16,17)
    // OpenBLAS❗✔️: #define INIT_m8n12 INIT_m8n10 unit_init_2zmm(18,19)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_k1m4n1 \
    // OpenBLAS❗✔️:   "vbroadcastsd (%1),%%ymm1; addq $8,%1;"\
    // OpenBLAS❗✔️:   "vfmadd231pd (%0),%%ymm1,%%ymm4; addq $32,%0;"
    // OpenBLAS❗✔️: #define unit_acc_gen_m4n2(c1_no,c2_no,boff,...)\
    // OpenBLAS❗✔️:   "vbroadcastf128 "#boff"("#__VA_ARGS__"),%%ymm3; vfmadd231pd %%ymm1,%%ymm3,%%ymm"#c1_no"; vfmadd231pd %%ymm2,%%ymm3,%%ymm"#c2_no";"
    // OpenBLAS❗✔️: #define unit_acc_m4n2(c1_no,c2_no,...) unit_acc_gen_m4n2(c1_no,c2_no,0,__VA_ARGS__)
    // OpenBLAS❗✔️: #define KERNEL_h_k1m4n2 \
    // OpenBLAS❗✔️:   "vmovddup (%0),%%ymm1; vmovddup 8(%0),%%ymm2; addq $32,%0;" unit_acc_m4n2(4,5,%1)
    // OpenBLAS❗✔️: #define KERNEL_k1m4n2 KERNEL_h_k1m4n2 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m4n4 KERNEL_h_k1m4n2 unit_acc_m4n2(6,7,%1,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m4n4 KERNEL_h_k1m4n4 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_k1m4n6 KERNEL_h_k1m4n4 unit_acc_m4n2(8,9,%1,%%r12,2) "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m4n8 KERNEL_k1m4n6 unit_acc_m4n2(10,11,%%r15)
    // OpenBLAS❗✔️: #define KERNEL_k1m4n8 KERNEL_h_k1m4n8 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m4n10 KERNEL_h_k1m4n8 unit_acc_m4n2(12,13,%%r15,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m4n10 KERNEL_h_k1m4n10 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m4n12 KERNEL_h_k1m4n10 unit_acc_m4n2(14,15,%%r15,%%r12,2)
    // OpenBLAS❗✔️: //#define KERNEL_k1m4n12 KERNEL_h_k1m4n12 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define unit_acc_k2m4n2(c1_no,c2_no,...)\
    // OpenBLAS❗✔️:   "vbroadcastf64x4 ("#__VA_ARGS__"),%%zmm3; vpermpd %%zmm3,%%zmm30,%%zmm3;"\
    // OpenBLAS❗✔️:   "vfmadd231pd %%zmm1,%%zmm3,%%zmm"#c1_no"; vfmadd231pd %%zmm2,%%zmm3,%%zmm"#c2_no";"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define unit_merge_to_ymm(c1_no) \
    // OpenBLAS❗✔️:   "vextractf64x4 $1,%%zmm"#c1_no",%%ymm30; vaddpd %%ymm"#c1_no",%%ymm30,%%ymm"#c1_no";"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_k1m4n12 \
    // OpenBLAS❗✔️:   "cmpq $2, %5; jb 104912f;"\
    // OpenBLAS❗✔️:   "vmovupd 64+%11,%%zmm30;"\
    // OpenBLAS❗✔️:   "\n204912:"\
    // OpenBLAS❗✔️:   "vmovddup (%0),%%zmm1; vmovddup 8(%0),%%zmm2; addq $64,%0;" \
    // OpenBLAS❗✔️:   unit_acc_k2m4n2(4,5,%1) unit_acc_k2m4n2(6,7,%1,%%r12,1) unit_acc_k2m4n2(8, 9, %1, %%r12, 2) "addq $32,%1;" \
    // OpenBLAS❗✔️:   unit_acc_k2m4n2(10,11,%%r15) unit_acc_k2m4n2(12,13,%%r15,%%r12,1) unit_acc_k2m4n2(14,15,%%r15,%%r12,2) "addq $32,%%r15;" \
    // OpenBLAS❗✔️:   "subq $2, %5; cmpq $2, %5; jnb 204912b;"\
    // OpenBLAS❗✔️:   unit_merge_to_ymm(4) unit_merge_to_ymm(5) unit_merge_to_ymm(6) unit_merge_to_ymm(7) \
    // OpenBLAS❗✔️:   unit_merge_to_ymm(8) unit_merge_to_ymm(9) unit_merge_to_ymm(10) unit_merge_to_ymm(11) \
    // OpenBLAS❗✔️:   unit_merge_to_ymm(12) unit_merge_to_ymm(13) unit_merge_to_ymm(14) unit_merge_to_ymm(15) \
    // OpenBLAS❗✔️:   "testq %5, %5; jz 1004912f;"\
    // OpenBLAS❗✔️:   "\n104912:"\
    // OpenBLAS❗✔️:   KERNEL_h_k1m4n12 "addq $16,%%r15;"\
    // OpenBLAS❗✔️:   "decq %5; jnz 104912b;"\
    // OpenBLAS❗✔️:   "\n1004912:"\
    // OpenBLAS❗✔️:   "incq %5;"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT) && (BACKWARDS == 0)
    // OpenBLAS❗✔️:   #define loada_kend_k1m4 "vmovddup (%0,%3,1),%%ymm1; vmovddup 8(%0,%3,1),%%ymm2; addq $32,%3;"
    // OpenBLAS❗✔️:   #define acc_kend_nc2_k1m4(boff1) unit_acc_gen_m4n2(6,7,boff1,%1,%%r12,1)
    // OpenBLAS❗✔️:   #define acc_kend_nc3_k1m4(boff1) unit_acc_gen_m4n2(8,9,boff1,%1,%%r12,2)
    // OpenBLAS❗✔️:   #define acc_kend_nc4_k1m4(boff1) unit_acc_gen_m4n2(10,11,boff1,%%r15)
    // OpenBLAS❗✔️:   #define acc_kend_nc5_k1m4(boff1) unit_acc_gen_m4n2(12,13,boff1,%%r15,%%r12,1)
    // OpenBLAS❗✔️:   #define acc_kend_nc6_k1m4(boff1) unit_acc_gen_m4n2(14,15,boff1,%%r15,%%r12,2)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define save_init_m4 "movq %2,%3; addq $32,%2;"
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #define SAVE_m4n1 "vmulpd %%ymm4,%%ymm0,%%ymm4; vmovupd %%ymm4,(%2); addq $32,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m4n2(c1_no,c2_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%ymm"#c2_no",%%ymm"#c1_no",%%ymm1; vmulpd %%ymm1,%%ymm0,%%ymm1; vmovupd %%ymm1,(%3);"\
    // OpenBLAS❗✔️:     "vunpckhpd %%ymm"#c2_no",%%ymm"#c1_no",%%ymm2; vmulpd %%ymm2,%%ymm0,%%ymm2; vmovupd %%ymm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define SAVE_m4n1 "vfmadd213pd (%2),%%ymm0,%%ymm4; vmovupd %%ymm4,(%2); addq $32,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m4n2(c1_no,c2_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%ymm"#c2_no",%%ymm"#c1_no",%%ymm1; vfmadd213pd (%3),%%ymm0,%%ymm1; vmovupd %%ymm1,(%3);"\
    // OpenBLAS❗✔️:     "vunpckhpd %%ymm"#c2_no",%%ymm"#c1_no",%%ymm2; vfmadd213pd (%3,%4,1),%%ymm0,%%ymm2; vmovupd %%ymm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define SAVE_m4n2 save_init_m4 unit_save_m4n2(4,5)
    // OpenBLAS❗✔️: #define SAVE_m4n4 SAVE_m4n2 unit_save_m4n2(6,7)
    // OpenBLAS❗✔️: #define SAVE_m4n6 SAVE_m4n4 unit_save_m4n2(8,9)
    // OpenBLAS❗✔️: #define SAVE_m4n8 SAVE_m4n6 unit_save_m4n2(10,11)
    // OpenBLAS❗✔️: #define SAVE_m4n10 SAVE_m4n8 unit_save_m4n2(12,13)
    // OpenBLAS❗✔️: #define SAVE_m4n12 SAVE_m4n10 unit_save_m4n2(14,15)
    // OpenBLAS❗✔️: #define INIT_m4n1 "vpxor %%ymm4,%%ymm4,%%ymm4;"
    // OpenBLAS❗✔️: #define unit_init_2ymm(c1_no,c2_no) "vpxor %%ymm"#c1_no",%%ymm"#c1_no",%%ymm"#c1_no"; vpxor %%ymm"#c2_no",%%ymm"#c2_no",%%ymm"#c2_no";"
    // OpenBLAS❗✔️: #define INIT_m4n2 unit_init_2ymm(4,5)
    // OpenBLAS❗✔️: #define INIT_m4n4 INIT_m4n2 unit_init_2ymm(6,7)
    // OpenBLAS❗✔️: #define INIT_m4n6 INIT_m4n4 unit_init_2ymm(8,9)
    // OpenBLAS❗✔️: #define INIT_m4n8 INIT_m4n6 unit_init_2ymm(10,11)
    // OpenBLAS❗✔️: #define INIT_m4n10 INIT_m4n8 unit_init_2ymm(12,13)
    // OpenBLAS❗✔️: #define INIT_m4n12 INIT_m4n10 unit_init_2ymm(14,15)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_k1m2n1 \
    // OpenBLAS❗✔️:   "vmovddup (%1),%%xmm1; addq $8,%1;"\
    // OpenBLAS❗✔️:   "vfmadd231pd (%0),%%xmm1,%%xmm4; addq $16,%0;"
    // OpenBLAS❗✔️: #define unit_acc_gen_m2n2(c1_no,c2_no,boff,...)\
    // OpenBLAS❗✔️:   "vmovupd "#boff"("#__VA_ARGS__"),%%xmm3; vfmadd231pd %%xmm1,%%xmm3,%%xmm"#c1_no"; vfmadd231pd %%xmm2,%%xmm3,%%xmm"#c2_no";"
    // OpenBLAS❗✔️: #define unit_acc_m2n2(c1_no,c2_no,...) unit_acc_gen_m2n2(c1_no,c2_no,0,__VA_ARGS__)
    // OpenBLAS❗✔️: #define KERNEL_h_k1m2n2 \
    // OpenBLAS❗✔️:   "vmovddup (%0),%%xmm1; vmovddup 8(%0),%%xmm2; addq $16,%0;" unit_acc_m2n2(4,5,%1)
    // OpenBLAS❗✔️: #define KERNEL_k1m2n2 KERNEL_h_k1m2n2 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m2n4 KERNEL_h_k1m2n2 unit_acc_m2n2(6,7,%1,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m2n4 KERNEL_h_k1m2n4 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_k1m2n6 KERNEL_h_k1m2n4 unit_acc_m2n2(8,9,%1,%%r12,2) "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m2n8 KERNEL_k1m2n6 unit_acc_m2n2(10,11,%%r15)
    // OpenBLAS❗✔️: #define KERNEL_k1m2n8 KERNEL_h_k1m2n8 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m2n10 KERNEL_h_k1m2n8 unit_acc_m2n2(12,13,%%r15,%%r12,1)
    // OpenBLAS❗✔️: #define KERNEL_k1m2n10 KERNEL_h_k1m2n10 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m2n12 KERNEL_h_k1m2n10 unit_acc_m2n2(14,15,%%r15,%%r12,2)
    // OpenBLAS❗✔️: //#define KERNEL_k1m2n12 KERNEL_h_k1m2n12 "addq $16,%%r15;"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define unit_acc_k4m2n2(c1_no,c2_no,...) \
    // OpenBLAS❗✔️:   "vmovupd ("#__VA_ARGS__"),%%zmm3; vfmadd231pd %%zmm1,%%zmm3,%%zmm"#c1_no"; vfmadd231pd %%zmm2,%%zmm3,%%zmm"#c2_no";"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define unit_merge_to_xmm(c1_no) \
    // OpenBLAS❗✔️:   "vextractf64x2 $0,%%zmm"#c1_no",%%xmm20; vextractf64x2 $1,%%zmm"#c1_no",%%xmm21; vextractf64x2 $2,%%zmm"#c1_no",%%xmm22; vextractf64x2 $3,%%zmm"#c1_no",%%xmm23;"\
    // OpenBLAS❗✔️:   "vaddpd %%xmm20,%%xmm21,%%xmm20; vaddpd %%xmm22,%%xmm23,%%xmm22; vaddpd %%xmm20,%%xmm22,%%xmm"#c1_no";"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_k1m2n12 \
    // OpenBLAS❗✔️:   "cmpq $4,%5; jb 102912f;"\
    // OpenBLAS❗✔️:   "\n402912:"\
    // OpenBLAS❗✔️:   "vmovddup (%0),%%zmm1; vmovddup 8(%0),%%zmm2; addq $64,%0;" \
    // OpenBLAS❗✔️:   unit_acc_k4m2n2(4,5,%1) unit_acc_k4m2n2(6,7,%1,%%r12,1) unit_acc_k4m2n2(8,9,%1,%%r12,2) "addq $64,%1;" \
    // OpenBLAS❗✔️:   unit_acc_k4m2n2(10,11,%%r15) unit_acc_k4m2n2(12,13,%%r15,%%r12,1) unit_acc_k4m2n2(14,15,%%r15,%%r12,2) "addq $64,%%r15;" \
    // OpenBLAS❗✔️:   "subq $4,%5; cmpq $4,%5; jnb 402912b;"\
    // OpenBLAS❗✔️:   unit_merge_to_xmm(4) unit_merge_to_xmm(5) unit_merge_to_xmm(6) unit_merge_to_xmm(7) unit_merge_to_xmm(8) unit_merge_to_xmm(9) \
    // OpenBLAS❗✔️:   unit_merge_to_xmm(10) unit_merge_to_xmm(11) unit_merge_to_xmm(12) unit_merge_to_xmm(13) unit_merge_to_xmm(14) unit_merge_to_xmm(15) \
    // OpenBLAS❗✔️:   "testq %5,%5; jz 1002912f;"\
    // OpenBLAS❗✔️:   "\n102912:"\
    // OpenBLAS❗✔️:   KERNEL_h_k1m2n12 "addq $16,%%r15;" \
    // OpenBLAS❗✔️:   "decq %5; jnz 102912b;" \
    // OpenBLAS❗✔️:   "\n1002912:"\
    // OpenBLAS❗✔️:   "incq %5;"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT) && (BACKWARDS == 0)
    // OpenBLAS❗✔️:   #define loada_kend_k1m2 "vmovddup (%0,%3,1),%%xmm1; vmovddup 8(%0,%3,1),%%xmm2; addq $16,%3;"
    // OpenBLAS❗✔️:   #define acc_kend_nc2_k1m2(boff1) unit_acc_gen_m2n2(6,7,boff1,%1,%%r12,1)
    // OpenBLAS❗✔️:   #define acc_kend_nc3_k1m2(boff1) unit_acc_gen_m2n2(8,9,boff1,%1,%%r12,2)
    // OpenBLAS❗✔️:   #define acc_kend_nc4_k1m2(boff1) unit_acc_gen_m2n2(10,11,boff1,%%r15)
    // OpenBLAS❗✔️:   #define acc_kend_nc5_k1m2(boff1) unit_acc_gen_m2n2(12,13,boff1,%%r15,%%r12,1)
    // OpenBLAS❗✔️:   #define acc_kend_nc6_k1m2(boff1) unit_acc_gen_m2n2(14,15,boff1,%%r15,%%r12,2)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define save_init_m2 "movq %2,%3; addq $16,%2;"
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #define SAVE_m2n1 "vmulpd %%xmm4,%%xmm0,%%xmm4; vmovupd %%xmm4,(%2); addq $16,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m2n2(c1_no,c2_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%xmm"#c2_no",%%xmm"#c1_no",%%xmm1; vmulpd %%xmm1,%%xmm0,%%xmm1; vmovupd %%xmm1,(%3);"\
    // OpenBLAS❗✔️:     "vunpckhpd %%xmm"#c2_no",%%xmm"#c1_no",%%xmm2; vmulpd %%xmm2,%%xmm0,%%xmm2; vmovupd %%xmm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define SAVE_m2n1 "vfmadd213pd (%2),%%xmm0,%%xmm4; vmovupd %%xmm4,(%2); addq $16,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m2n2(c1_no,c2_no)\
    // OpenBLAS❗✔️:     "vunpcklpd %%xmm"#c2_no",%%xmm"#c1_no",%%xmm1; vfmadd213pd (%3),%%xmm0,%%xmm1; vmovupd %%xmm1,(%3);"\
    // OpenBLAS❗✔️:     "vunpckhpd %%xmm"#c2_no",%%xmm"#c1_no",%%xmm2; vfmadd213pd (%3,%4,1),%%xmm0,%%xmm2; vmovupd %%xmm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define SAVE_m2n2 save_init_m2 unit_save_m2n2(4,5)
    // OpenBLAS❗✔️: #define SAVE_m2n4 SAVE_m2n2 unit_save_m2n2(6,7)
    // OpenBLAS❗✔️: #define SAVE_m2n6 SAVE_m2n4 unit_save_m2n2(8,9)
    // OpenBLAS❗✔️: #define SAVE_m2n8 SAVE_m2n6 unit_save_m2n2(10,11)
    // OpenBLAS❗✔️: #define SAVE_m2n10 SAVE_m2n8 unit_save_m2n2(12,13)
    // OpenBLAS❗✔️: #define SAVE_m2n12 SAVE_m2n10 unit_save_m2n2(14,15)
    // OpenBLAS❗✔️: #define INIT_m2n1 "vpxor %%xmm4,%%xmm4,%%xmm4;"
    // OpenBLAS❗✔️: #define unit_init_2xmm(c1_no,c2_no) "vpxor %%xmm"#c1_no",%%xmm"#c1_no",%%xmm"#c1_no"; vpxor %%xmm"#c2_no",%%xmm"#c2_no",%%xmm"#c2_no";"
    // OpenBLAS❗✔️: #define INIT_m2n2 unit_init_2xmm(4,5)
    // OpenBLAS❗✔️: #define INIT_m2n4 INIT_m2n2 unit_init_2xmm(6,7)
    // OpenBLAS❗✔️: #define INIT_m2n6 INIT_m2n4 unit_init_2xmm(8,9)
    // OpenBLAS❗✔️: #define INIT_m2n8 INIT_m2n6 unit_init_2xmm(10,11)
    // OpenBLAS❗✔️: #define INIT_m2n10 INIT_m2n8 unit_init_2xmm(12,13)
    // OpenBLAS❗✔️: #define INIT_m2n12 INIT_m2n10 unit_init_2xmm(14,15)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define KERNEL_k1m1n1 \
    // OpenBLAS❗✔️:   "vmovsd (%1),%%xmm1; addq $8,%1;"\
    // OpenBLAS❗✔️:   "vfmadd231sd (%0),%%xmm1,%%xmm4; addq $8,%0;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m1n2 \
    // OpenBLAS❗✔️:   "vmovddup (%0),%%xmm1; addq $8,%0;"\
    // OpenBLAS❗✔️:   "vfmadd231pd (%1),%%xmm1,%%xmm4;"
    // OpenBLAS❗✔️: #define KERNEL_k1m1n2 KERNEL_h_k1m1n2 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m1n4 KERNEL_h_k1m1n2 "vfmadd231pd (%1,%%r12,1),%%xmm1,%%xmm5;"
    // OpenBLAS❗✔️: #define KERNEL_k1m1n4 KERNEL_h_k1m1n4 "addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_k1m1n6 KERNEL_h_k1m1n4 "vfmadd231pd (%1,%%r12,2),%%xmm1,%%xmm6; addq $16,%1;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m1n8 KERNEL_k1m1n6 "vfmadd231pd (%%r15),%%xmm1,%%xmm7;"
    // OpenBLAS❗✔️: #define KERNEL_k1m1n8 KERNEL_h_k1m1n8 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m1n10 KERNEL_h_k1m1n8 "vfmadd231pd (%%r15,%%r12,1),%%xmm1,%%xmm8;"
    // OpenBLAS❗✔️: #define KERNEL_k1m1n10 KERNEL_h_k1m1n10 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_h_k1m1n12 KERNEL_h_k1m1n10 "vfmadd231pd (%%r15,%%r12,2),%%xmm1,%%xmm9;"
    // OpenBLAS❗✔️: //#define KERNEL_k1m1n12 KERNEL_h_k1m1n12 "addq $16,%%r15;"
    // OpenBLAS❗✔️: #define KERNEL_k1m1n12 \
    // OpenBLAS❗✔️:   "cmpq $4,%5; jb 101912f;" \
    // OpenBLAS❗✔️:   "vmovupd %11,%%zmm2;"\
    // OpenBLAS❗✔️:   "\n401912:"\
    // OpenBLAS❗✔️:   "vmovupd (%0),%%ymm1; vpermpd %%zmm1,%%zmm2,%%zmm1; addq $32,%0;" \
    // OpenBLAS❗✔️:   "vfmadd231pd (%1),%%zmm1,%%zmm4; vfmadd231pd (%1,%%r12,1),%%zmm1,%%zmm5; vfmadd231pd (%1,%%r12,2),%%zmm1,%%zmm6; addq $64,%1;"\
    // OpenBLAS❗✔️:   "vfmadd231pd (%%r15),%%zmm1,%%zmm7; vfmadd231pd (%%r15,%%r12,1),%%zmm1,%%zmm8; vfmadd231pd (%%r15,%%r12,2),%%zmm1,%%zmm9; addq $64,%%r15;"\
    // OpenBLAS❗✔️:   "subq $4,%5; cmpq $4,%5; jnb 401912b;"\
    // OpenBLAS❗✔️:   unit_merge_to_xmm(4) unit_merge_to_xmm(5) unit_merge_to_xmm(6) \
    // OpenBLAS❗✔️:   unit_merge_to_xmm(7) unit_merge_to_xmm(8) unit_merge_to_xmm(9) \
    // OpenBLAS❗✔️:   "testq %5,%5; jz 1001912f;"\
    // OpenBLAS❗✔️:   "\n101912:"\
    // OpenBLAS❗✔️:   KERNEL_h_k1m1n12 "addq $16,%%r15;" \
    // OpenBLAS❗✔️:   "decq %5; jnz 101912b;" \
    // OpenBLAS❗✔️:   "\n1001912:"\
    // OpenBLAS❗✔️:   "incq %5;"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT) && (BACKWARDS == 0)
    // OpenBLAS❗✔️:   #define loada_kend_k1m1 "vmovddup (%0,%3,1),%%xmm1; addq $8,%3;"
    // OpenBLAS❗✔️:   #define acc_kend_nc2_k1m1(boff1) "vfmadd231pd "#boff1"(%1,%%r12,1),%%xmm1,%%xmm5;"
    // OpenBLAS❗✔️:   #define acc_kend_nc3_k1m1(boff1) "vfmadd231pd "#boff1"(%1,%%r12,2),%%xmm1,%%xmm6;"
    // OpenBLAS❗✔️:   #define acc_kend_nc4_k1m1(boff1) "vfmadd231pd "#boff1"(%%r15),%%xmm1,%%xmm7;"
    // OpenBLAS❗✔️:   #define acc_kend_nc5_k1m1(boff1) "vfmadd231pd "#boff1"(%%r15,%%r12,1),%%xmm1,%%xmm8;"
    // OpenBLAS❗✔️:   #define acc_kend_nc6_k1m1(boff1) "vfmadd231pd "#boff1"(%%r15,%%r12,2),%%xmm1,%%xmm9;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define save_init_m1 "movq %2,%3; addq $8,%2;"
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #define SAVE_m1n1 "vmulsd %%xmm4,%%xmm0,%%xmm4; vmovsd %%xmm4,(%2); addq $8,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m1n2(c1_no)\
    // OpenBLAS❗✔️:     "vmulpd %%xmm"#c1_no",%%xmm0,%%xmm2; vmovsd %%xmm2,(%3); vmovhpd %%xmm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define SAVE_m1n1 "vfmadd213sd (%2),%%xmm0,%%xmm4; vmovsd %%xmm4,(%2); addq $8,%2;"
    // OpenBLAS❗✔️:   #define unit_save_m1n2(c1_no)\
    // OpenBLAS❗✔️:     "vmovsd (%3),%%xmm2; vmovhpd (%3,%4,1),%%xmm2,%%xmm2; vfmadd231pd %%xmm"#c1_no",%%xmm0,%%xmm2; vmovsd %%xmm2,(%3); vmovhpd %%xmm2,(%3,%4,1); leaq (%3,%4,2),%3;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define SAVE_m1n2 save_init_m1 unit_save_m1n2(4)
    // OpenBLAS❗✔️: #define SAVE_m1n4 SAVE_m1n2 unit_save_m1n2(5)
    // OpenBLAS❗✔️: #define SAVE_m1n6 SAVE_m1n4 unit_save_m1n2(6)
    // OpenBLAS❗✔️: #define SAVE_m1n8 SAVE_m1n6 unit_save_m1n2(7)
    // OpenBLAS❗✔️: #define SAVE_m1n10 SAVE_m1n8 unit_save_m1n2(8)
    // OpenBLAS❗✔️: #define SAVE_m1n12 SAVE_m1n10 unit_save_m1n2(9)
    // OpenBLAS❗✔️: #define INIT_m1n1 "vpxor %%xmm4,%%xmm4,%%xmm4;"
    // OpenBLAS❗✔️: #define INIT_m1n2 INIT_m1n1
    // OpenBLAS❗✔️: #define INIT_m1n4 INIT_m1n2 "vpxor %%xmm5,%%xmm5,%%xmm5;"
    // OpenBLAS❗✔️: #define INIT_m1n6 INIT_m1n4 "vpxor %%xmm6,%%xmm6,%%xmm6;"
    // OpenBLAS❗✔️: #define INIT_m1n8 INIT_m1n6 "vpxor %%xmm7,%%xmm7,%%xmm7;"
    // OpenBLAS❗✔️: #define INIT_m1n10 INIT_m1n8 "vpxor %%xmm8,%%xmm8,%%xmm8;"
    // OpenBLAS❗✔️: #define INIT_m1n12 INIT_m1n10 "vpxor %%xmm9,%%xmm9,%%xmm9;"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define COMPUTE_SIMPLE(mdim,ndim)\
    // OpenBLAS❗✔️:   init_update_k(mdim) INIT_m##mdim##n##ndim\
    // OpenBLAS❗✔️:   "movq %%r13,%5;" INIT_set_papb(mdim,ndim)\
    // OpenBLAS❗✔️:   kernel_kstart_n##ndim(mdim,subq)\
    // OpenBLAS❗✔️:   "testq %5,%5; jz 7"#mdim"7"#ndim"9f;"\
    // OpenBLAS❗✔️:   "7"#mdim"7"#ndim"1:\n\t"\
    // OpenBLAS❗✔️:   KERNEL_k1m##mdim##n##ndim "decq %5; jnz 7"#mdim"7"#ndim"1b;"\
    // OpenBLAS❗✔️:   "7"#mdim"7"#ndim"9:\n\t"\
    // OpenBLAS❗✔️:   kernel_kend_n##ndim(mdim)\
    // OpenBLAS❗✔️:   SAVE_set_pa(mdim) SAVE_m##mdim##n##ndim save_update_k(mdim)
    // OpenBLAS❗✔️: #define COMPUTE_m16n1 COMPUTE_SIMPLE(16,1)
    // OpenBLAS❗✔️: #define COMPUTE_m16n2 COMPUTE_SIMPLE(16,2)
    // OpenBLAS❗✔️: #define COMPUTE_m16n4 COMPUTE_SIMPLE(16,4)
    // OpenBLAS❗✔️: #define COMPUTE_m16n6 COMPUTE_SIMPLE(16,6)
    // OpenBLAS❗✔️: #define COMPUTE_m16n8 COMPUTE_SIMPLE(16,8)
    // OpenBLAS❗✔️: #define COMPUTE_m16n10 COMPUTE_SIMPLE(16,10)
    // OpenBLAS❗✔️: #if defined(TRMMKERNEL) && !defined(LEFT) && defined(TRANSA)
    // OpenBLAS❗✔️:   #define INVERSE_K_MID "negq %5; leaq 6(%%r13,%5,1),%5;"
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   #define INVERSE_K_MID "negq %5; leaq 16(%%r13,%5,1),%5;"
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define COMPUTE_m16n12 \
    // OpenBLAS❗✔️:   init_update_k(16) INIT_m16n12 "movq %%r13,%5;" INIT_set_papb(16,12) "movq %2,%3;"\
    // OpenBLAS❗✔️:   kernel_kstart_n12(16,subq)\
    // OpenBLAS❗✔️:   "cmpq $16,%5; jb 7167123f; movq $16,%5;"\
    // OpenBLAS❗✔️:   "7167121:\n\t"\
    // OpenBLAS❗✔️:   KERNEL_k1m16n12 "addq $4,%5; testq $12,%5; movq $172,%%r10; cmovz %4,%%r10;"\
    // OpenBLAS❗✔️:   KERNEL_k1m16n12 "prefetcht1 (%3); subq $129,%3; addq %%r10,%3;"\
    // OpenBLAS❗✔️:   KERNEL_k1m16n12 "prefetcht1 (%6); addq $32,%6; cmpq $208,%5; cmoveq %2,%3;"\
    // OpenBLAS❗✔️:   KERNEL_k1m16n12 "cmpq %5,%%r13; jnb 7167121b;"\
    // OpenBLAS❗✔️:   "movq %2,%3;" INVERSE_K_MID\
    // OpenBLAS❗✔️:   "7167123:\n\t"\
    // OpenBLAS❗✔️:   "testq %5,%5; jz 7167129f;"\
    // OpenBLAS❗✔️:   "7167125:\n\t"\
    // OpenBLAS❗✔️:   "prefetcht0 (%3); prefetcht0 64(%3); prefetcht0 127(%3);"\
    // OpenBLAS❗✔️:   KERNEL_k1m16n12 "addq %4,%3; decq %5;jnz 7167125b;"\
    // OpenBLAS❗✔️:   "7167129:\n\t"\
    // OpenBLAS❗✔️:   kernel_kend_n12(16)\
    // OpenBLAS❗✔️:   "prefetcht0 (%%r14);" SAVE_set_pa(16) SAVE_m16n12 save_update_k(16)
    // OpenBLAS❗✔️: #define COMPUTE(ndim) {\
    // OpenBLAS❗✔️:   b_pref = b_ptr + ndim * K; HEAD_SET_OFF(ndim)\
    // OpenBLAS❗✔️:   __asm__ __volatile__(\
    // OpenBLAS❗✔️:     "vbroadcastsd %8,%%zmm0; movq %7,%%r11; movq %1,%%r14; movq %10,%%r12; salq $4,%%r12;" INITASM_SET_K\
    // OpenBLAS❗✔️:     "cmpq $16,%%r11; jb "#ndim"33102f;"\
    // OpenBLAS❗✔️:     #ndim"33101:\n\t"\
    // OpenBLAS❗✔️:     COMPUTE_m16n##ndim "subq $16,%%r11; cmpq $16,%%r11; jnb "#ndim"33101b;"\
    // OpenBLAS❗✔️:     #ndim"33102:\n\t"\
    // OpenBLAS❗✔️:     "cmpq $8,%%r11; jb "#ndim"33103f;"\
    // OpenBLAS❗✔️:     COMPUTE_SIMPLE(8,ndim) "subq $8,%%r11;"\
    // OpenBLAS❗✔️:     #ndim"33103:\n\t"\
    // OpenBLAS❗✔️:     "cmpq $4,%%r11; jb "#ndim"33104f;"\
    // OpenBLAS❗✔️:     COMPUTE_SIMPLE(4,ndim) "subq $4,%%r11;"\
    // OpenBLAS❗✔️:     #ndim"33104:\n\t"\
    // OpenBLAS❗✔️:     "cmpq $2,%%r11; jb "#ndim"33105f;"\
    // OpenBLAS❗✔️:     COMPUTE_SIMPLE(2,ndim) "subq $2,%%r11;"\
    // OpenBLAS❗✔️:     #ndim"33105:\n\t"\
    // OpenBLAS❗✔️:     "testq %%r11,%%r11; jz "#ndim"33106f;"\
    // OpenBLAS❗✔️:     COMPUTE_SIMPLE(1,ndim) "subq $1,%%r11;"\
    // OpenBLAS❗✔️:     #ndim"33106:\n\t"\
    // OpenBLAS❗✔️:     "movq %%r14,%1;"\
    // OpenBLAS❗✔️:   :"+r"(a_ptr),"+r"(b_ptr),"+r"(c_ptr),"+r"(c_tmp),"+r"(ldc_in_bytes),"+r"(k_count),"+r"(b_pref):"m"(M),"m"(ALPHA),"m"(off),"m"(K), "o"(permute_table):"r10","r11","r12","r13","r14","r15","cc","memory",\
    // OpenBLAS❗✔️:     "zmm0","zmm1","zmm2","zmm3","zmm4","zmm5","zmm6","zmm7","zmm8","zmm9","zmm10","zmm11","zmm12","zmm13","zmm14","zmm15",\
    // OpenBLAS❗✔️:     "zmm16","zmm17","zmm18","zmm19","zmm20","zmm21","zmm22","zmm23","zmm24","zmm25","zmm26","zmm27","zmm28","zmm29","zmm30","zmm31");\
    // OpenBLAS❗✔️:   a_ptr -= M * K; b_ptr += ndim * K; c_ptr += ndim * ldc - M; TAIL_SET_OFF(ndim)\
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #include "common.h"
    // OpenBLAS❗✔️: #include <stdint.h>
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: int __attribute__ ((noinline))
    // OpenBLAS❗✔️: CNAME(BLASLONG m, BLASLONG n, BLASLONG k, double alpha, double * __restrict__ A, double * __restrict__ B, double * __restrict__ C, BLASLONG ldc
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   , BLASLONG offset
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: )
    // OpenBLAS❗✔️: {
    // OpenBLAS❗✔️:     if(m==0||n==0) return 0;
    // OpenBLAS❗✔️:     int64_t ldc_in_bytes = (int64_t)ldc * sizeof(double); double ALPHA = alpha;
    // OpenBLAS❗✔️:     int64_t M = (int64_t)m, K = (int64_t)k, k_count = 0;
    // OpenBLAS❗✔️:     BLASLONG n_count = n, off = 0;
    // OpenBLAS❗✔️:     double *a_ptr = A,*b_ptr = B,*c_ptr = C,*c_tmp = C,*b_pref = B;
    // OpenBLAS❗✔️:     int64_t permute_table[] = {
    // OpenBLAS❗✔️: 	    0, 0, 1, 1, 2, 2, 3, 3,  // abcdxxxx -> aabbccdd
    // OpenBLAS❗✔️: 	    0, 1, 0, 1, 2, 3, 2, 3,  // abcdxxxx -> ababcdcd
    // OpenBLAS❗✔️:     };
    // OpenBLAS❗✔️: #ifdef TRMMKERNEL
    // OpenBLAS❗✔️:   #ifdef LEFT
    // OpenBLAS❗✔️:     off = offset;
    // OpenBLAS❗✔️:   #else
    // OpenBLAS❗✔️:     off = -offset;
    // OpenBLAS❗✔️:   #endif
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     for(;n_count>11;n_count-=12) COMPUTE(12)
    // OpenBLAS❗✔️:     for(;n_count>9;n_count-=10) COMPUTE(10)
    // OpenBLAS❗✔️:     for(;n_count>7;n_count-=8) COMPUTE(8)
    // OpenBLAS❗✔️:     for(;n_count>5;n_count-=6) COMPUTE(6)
    // OpenBLAS❗✔️:     for(;n_count>3;n_count-=4) COMPUTE(4)
    // OpenBLAS❗✔️:     for(;n_count>1;n_count-=2) COMPUTE(2)
    // OpenBLAS❗✔️:     if(n_count>0) COMPUTE(1)
    // OpenBLAS❗✔️:     return 0;
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️:
    #[cfg(target_arch = "x86_64")]
    let native = source_simd();
    let fma_rows = row_fma();
    let mut j = 0;
    while j < n {
        let w = col_tile(n - j);
        let mut i = 0;
        while i < m {
            let h = row_tile(m - i);
            let lanes = if w == 12 && h == 4 {
                2
            } else if w == 12 && h <= 2 {
                4
            } else {
                1
            };
            let mut begin = 0;
            let mut end = k;
            if let Some(t) = tri {
                match (t.side, t.upper) {
                    (Side::Left, true) => begin = (t.diagonal + i as isize) as usize,
                    (Side::Left, false) => end = (t.diagonal + i as isize + h as isize) as usize,
                    (Side::Right, true) => {
                        end = (t.diagonal + j as isize + w.min(2) as isize) as usize
                    }
                    (Side::Right, false) => begin = (t.diagonal + j as isize) as usize,
                }
            }
            #[cfg(target_arch = "x86_64")]
            if native {
                if h >= 8 && tri.is_none() {
                    let pa = &a[i * k..i * k + h * k];
                    let pb = &b[j * k..(j + w) * k];
                    let ci = offset + i + j * ldc;
                    // SAFETY: paired source_simd guard, complete eight-row
                    // groups, checked packed inputs and active output region.
                    unsafe {
                        if h == 16 {
                            native_wide_packed_dispatch::<16, 2>(
                                w, k, alpha, pa, pb, c, ci, ldc, begin, end, false,
                            );
                        } else {
                            native_wide_packed_dispatch::<8, 1>(
                                w, k, alpha, pa, pb, c, ci, ldc, begin, end, false,
                            );
                        }
                    }
                    i += h;
                    continue;
                }
                for pair in (0..w).step_by(2) {
                    let bw = (w - pair).min(2);
                    let right_lower = matches!(
                        tri,
                        Some(TriKernel {
                            side: Side::Right,
                            upper: false,
                            ..
                        })
                    );
                    let right_upper = matches!(
                        tri,
                        Some(TriKernel {
                            side: Side::Right,
                            upper: true,
                            ..
                        })
                    );
                    let prefix = begin + if right_lower { pair } else { 0 };
                    let main = begin + if right_lower { w.saturating_sub(2) } else { 0 };
                    let suffix = end + if right_upper { pair } else { 0 };
                    let pa = &a[i * k..i * k + h * k];
                    let pb = &b[(j + pair) * k..(j + pair + bw) * k];
                    let ci = offset + i + (j + pair) * ldc;
                    // SAFETY: all source k ranges lie inside the checked pack;
                    // output coordinates and CPU feature guard were validated.
                    unsafe {
                        match (h, lanes) {
                            (16, 1) => native_wide_packed_dispatch::<16, 2>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                suffix,
                                tri.is_some(),
                            ),
                            (8, 1) => native_wide_packed_dispatch::<8, 1>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                suffix,
                                tri.is_some(),
                            ),
                            (4, 2) => native_tail_dispatch::<4, 2>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                main,
                                end,
                                suffix,
                                tri.is_some(),
                            ),
                            (2, 4) => native_tail_dispatch::<2, 4>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                main,
                                end,
                                suffix,
                                tri.is_some(),
                            ),
                            (1, 4) => native_tail_dispatch::<1, 4>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                main,
                                end,
                                suffix,
                                tri.is_some(),
                            ),
                            (4, 1) => native_tail_dispatch::<4, 1>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                main,
                                end,
                                suffix,
                                tri.is_some(),
                            ),
                            (2, 1) => native_tail_dispatch::<2, 1>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                main,
                                end,
                                suffix,
                                tri.is_some(),
                            ),
                            _ => native_tail_dispatch::<1, 1>(
                                bw,
                                k,
                                alpha,
                                pa,
                                pb,
                                c,
                                ci,
                                ldc,
                                prefix,
                                main,
                                end,
                                suffix,
                                tri.is_some(),
                            ),
                        }
                    }
                }
                i += h;
                continue;
            }
            let mut sums = [[[0_f64; 16]; 4]; 12];
            // Right/lower source expands the leading staircase two k values
            // at a time before the shared full-width recurrence.
            if matches!(
                tri,
                Some(TriKernel {
                    side: Side::Right,
                    upper: false,
                    ..
                })
            ) {
                for prefix in (0..w.saturating_sub(2)).step_by(2) {
                    for l in begin..begin + 2 {
                        accumulate(
                            fma_rows,
                            &mut sums,
                            h,
                            w,
                            k,
                            i,
                            j,
                            l,
                            0,
                            0,
                            prefix + 2,
                            a,
                            b,
                            n,
                        );
                    }
                    begin += 2;
                }
            }
            let bulk = (end - begin) / lanes * lanes + begin;
            for l in begin..bulk {
                accumulate(
                    fma_rows,
                    &mut sums,
                    h,
                    w,
                    k,
                    i,
                    j,
                    l,
                    (l - begin) % lanes,
                    0,
                    w,
                    a,
                    b,
                    n,
                );
            }
            if lanes > 1 && bulk > begin {
                for jj in 0..w {
                    for ii in 0..h {
                        let s = [
                            sums[jj][0][ii],
                            sums[jj][1][ii],
                            sums[jj][2][ii],
                            sums[jj][3][ii],
                        ];
                        sums[jj][0][ii] = if lanes == 2 {
                            s[0] + s[1]
                        } else {
                            (s[0] + s[1]) + (s[2] + s[3])
                        };
                    }
                }
            }
            for l in bulk..end {
                accumulate(fma_rows, &mut sums, h, w, k, i, j, l, 0, 0, w, a, b, n);
            }
            // Right/upper source fills the remaining column-pair staircase
            // after the common first-pair recurrence.
            if matches!(
                tri,
                Some(TriKernel {
                    side: Side::Right,
                    upper: true,
                    ..
                })
            ) {
                for suffix in (2..w).step_by(2) {
                    for l in end..end + 2 {
                        accumulate(fma_rows, &mut sums, h, w, k, i, j, l, 0, suffix, w, a, b, n);
                    }
                    end += 2;
                }
            }
            for jj in 0..w {
                for ii in 0..h {
                    let ci = offset + i + ii + (j + jj) * ldc;
                    c[ci] = if tri.is_some() {
                        alpha * sums[jj][0][ii]
                    } else {
                        alpha.mul_add(sums[jj][0][ii], c[ci])
                    };
                }
            }
            i += h;
        }
        j += w;
    }
}
#[inline]
fn accumulate(
    fma_rows: RowFma,
    sums: &mut [[[f64; 16]; 4]; 12],
    h: usize,
    _w: usize,
    k: usize,
    i: usize,
    j: usize,
    l: usize,
    lane: usize,
    from: usize,
    to: usize,
    a: &[f64],
    b: &[f64],
    n: usize,
) {
    for jj in from..to {
        let pair = jj / 2 * 2;
        let bw = (n - j - pair).min(2);
        let bv = b[(j + pair) * k + l * bw + jj - pair];
        fma_rows(
            &a[i * k + l * h..i * k + l * h + h],
            bv,
            &mut sums[jj][lane][..h],
        );
    }
}
#[inline]
fn reduce8(v: [f64; 8], four: bool) -> f64 {
    // OpenBLAS❗✔️: #define REDUCE_4(rr0, rr1, rr2, rr3) \
    // OpenBLAS❗✔️: 	__m512d r0, r1, r2, r3, t0, t1, t2, t3;\
    // OpenBLAS❗✔️: 	r0 = _mm512_unpacklo_pd(rr0, rr1); r1 = _mm512_unpackhi_pd(rr0, rr1); \
    // OpenBLAS❗✔️: 	r2 = _mm512_unpacklo_pd(rr2, rr3); r3 = _mm512_unpackhi_pd(rr2, rr3); \
    // OpenBLAS❗✔️: 	t0 = _mm512_permutex2var_pd(r0, idx_lo, r2); t1 = _mm512_permutex2var_pd(r1, idx_lo, r3); \
    // OpenBLAS❗✔️: 	t2 = _mm512_permutex2var_pd(r0, idx_hi, r2); t3 = _mm512_permutex2var_pd(r1, idx_hi, r3); \
    // OpenBLAS❗✔️: 	r0 = _mm512_add_pd(t0, t1); r1 = _mm512_add_pd(t2, t3); t0 = _mm512_add_pd(r0, r1); \
    // OpenBLAS❗✔️: 	__m256d s0, s1; \
    // OpenBLAS❗✔️: 	s0 = _mm512_extractf64x4_pd(t0, 0); s1 = _mm512_extractf64x4_pd(t0, 1); \
    // OpenBLAS❗✔️: 	s0 = _mm256_add_pd(s0, s1); s0 = _mm256_mul_pd(alpha_256, s0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define REDUCE_M4(N) REDUCE_4(result0##N, result1##N, result2##N, result3##N)
    // OpenBLAS❗✔️: #define REDUCE_N4(M) REDUCE_4(result##M##0, result##M##1, result##M##2, result##M##3)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: #define STORE_REDUCE(M, N) C[(j+N)*ldc + i + M] = alpha * _mm512_reduce_add_pd(result##M##N)
    // OpenBLAS❗✔️: #define STORE_M4(N, s0) _mm256_storeu_pd(&C[(j + N)*ldc + i], s0);
    // OpenBLAS❗✔️: #define STORE_N4(M, s0) _mm256_i64scatter_pd(&C[j*ldc + i + M], vindex_n, s0, 8);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: #define STORE_REDUCE(M, N) C[(j+N)*ldc + i + M] = alpha * _mm512_reduce_add_pd(result##M##N) + beta * C[(j+N)*ldc + i + M]
    // OpenBLAS❗✔️: #define STORE_M4(N, s0) \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0": "+v"(s0):"r"(&C[(j + N)*ldc + i]), "v"(beta_256)); \
    // OpenBLAS❗✔️: 	_mm256_storeu_pd(&C[(j + N)*ldc + i], s0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define STORE_N4(M, s0) \
    // OpenBLAS❗✔️: 	s0 = _mm256_fmadd_pd(_mm256_i64gather_pd(&C[j*ldc + i + M], vindex_n, 8), beta_256, s0); \
    // OpenBLAS❗✔️: 	_mm256_i64scatter_pd(&C[j*ldc + i + M], vindex_n, s0, 8);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #define STORE_REDUCE_M4(N) {\
    // OpenBLAS❗✔️: 	REDUCE_M4(N) \
    // OpenBLAS❗✔️: 	STORE_M4(N, s0) \
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️: #define STORE_REDUCE_N4(M) {\
    // OpenBLAS❗✔️: 	REDUCE_N4(M) \
    // OpenBLAS❗✔️: 	STORE_N4(M, s0) \
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT * C, BLASLONG ldc)
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT beta, FLOAT * C, BLASLONG ldc)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: {
    // OpenBLAS❗✔️: 	// column major
    // OpenBLAS❗✔️: 	BLASLONG i, j, k;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	BLASLONG m4 = M & ~3;
    // OpenBLAS❗✔️: 	BLASLONG m2 = M & ~1;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	BLASLONG n4 = N & ~3;
    if four {
        ((v[0] + v[1]) + (v[2] + v[3])) + ((v[4] + v[5]) + (v[6] + v[7]))
    } else {
        ((v[6] + v[2]) + (v[4] + v[0])) + ((v[7] + v[3]) + (v[5] + v[1]))
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_small_direct<const SWAP_TAIL: bool>(
    n: usize,
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: Matrix<'_>,
    beta: f64,
    c: &mut [f64],
    offset: usize,
    row_stride: usize,
    col_stride: usize,
    begin: usize,
    end: usize,
) {
    let mut i = begin;
    while i < end {
        let remaining = end - i;
        if SWAP_TAIL && remaining < 8 {
            // Source NT residuals broadcast A and vectorize contiguous B.
            // Reuse this same executor with its output axes exchanged once.
            let av = Matrix {
                trans: Transpose::NoTrans,
                ..b
            };
            let bv = Matrix {
                offset: a.offset + i,
                trans: Transpose::Trans,
                ..a
            };
            unsafe {
                native_small_direct::<false>(
                    remaining,
                    k,
                    alpha,
                    av,
                    bv,
                    beta,
                    c,
                    offset + i * row_stride,
                    col_stride,
                    row_stride,
                    0,
                    n,
                );
            }
            return;
        }
        let h = if remaining >= 32 {
            32
        } else if remaining >= 16 {
            16
        } else if remaining >= 8 {
            8
        } else {
            remaining
        };
        let mut j = 0;
        while j < n {
            let limit = if h == 32 && b.trans == Transpose::NoTrans {
                4
            } else {
                6
            };
            let left = n - j;
            let w = if left >= limit {
                limit
            } else if left >= 4 {
                4
            } else if left >= 2 {
                2
            } else {
                1
            };
            let ci = offset + i * row_stride + j * col_stride;
            // SAFETY: source_simd guard and active source/output dimensions
            // are checked by the outer entry. H/G are closed source row tiles.
            unsafe {
                match h {
                    32 => native_wide_direct_dispatch::<32, 4, false>(
                        32, w, k, alpha, a, b, i, j, beta, c, ci, row_stride, col_stride,
                    ),
                    16 => native_wide_direct_dispatch::<16, 2, false>(
                        16, w, k, alpha, a, b, i, j, beta, c, ci, row_stride, col_stride,
                    ),
                    8 => native_wide_direct_dispatch::<8, 1, false>(
                        8, w, k, alpha, a, b, i, j, beta, c, ci, row_stride, col_stride,
                    ),
                    _ => native_wide_direct_dispatch::<8, 1, true>(
                        h, w, k, alpha, a, b, i, j, beta, c, ci, row_stride, col_stride,
                    ),
                }
            }
            j += w;
        }
        i += h;
    }
}

fn small_direct(
    _m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: Matrix<'_>,
    beta: f64,
    c: &mut [f64],
    offset: usize,
    ldc: usize,
    begin: usize,
    end: usize,
) {
    // OpenBLAS❗❗: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // OpenBLAS❗❗: result##M##N = _mm512_mul_pd(result##M##N, alpha_512);
    // OpenBLAS❗❗: asm("vfmadd231pd (%1), %2, %0": "+v"(result##M##N):"r"(&C[(j+N)*ldc + i + (M*8)]), "v"(beta_512));
    // All output lanes keep increasing-k FMA order. Independent rows/columns
    // are tiled together so A's source vector is reused across the columns.
    #[cfg(target_arch = "x86_64")]
    if source_simd() {
        // SAFETY: the feature guard is paired with the helper's target_feature;
        // all source/output dimensions and active independence were validated.
        unsafe {
            if a.trans == Transpose::NoTrans {
                if b.trans == Transpose::Trans {
                    native_small_direct::<true>(
                        n, k, alpha, a, b, beta, c, offset, 1, ldc, begin, end,
                    );
                } else {
                    native_small_direct::<false>(
                        n, k, alpha, a, b, beta, c, offset, 1, ldc, begin, end,
                    );
                }
            } else {
                // Source TT vectorizes contiguous columns of B and scatters C.
                // This swaps only independent output axes, never summation order.
                let av = Matrix {
                    trans: Transpose::NoTrans,
                    ..b
                };
                let bv = Matrix {
                    trans: Transpose::NoTrans,
                    ..a
                };
                native_small_direct::<false>(
                    end - begin,
                    k,
                    alpha,
                    av,
                    bv,
                    beta,
                    c,
                    offset + begin,
                    ldc,
                    1,
                    0,
                    n,
                );
            }
        }
        return;
    }
    let fma = row_fma();
    for i in (begin..end).step_by(32) {
        let h = (end - i).min(32);
        for j in (0..n).step_by(4) {
            let w = (n - j).min(4);
            let mut sums = [[0_f64; 32]; 4];
            for l in 0..k {
                let mut gathered = [0_f64; 32];
                let av = if a.trans == Transpose::NoTrans {
                    let p = a.offset + i + l * a.ld;
                    &a.separate.unwrap_or(c)[p..p + h]
                } else {
                    for (ii, x) in gathered[..h].iter_mut().enumerate() {
                        *x = a.get(c, i + ii, l);
                    }
                    &gathered[..h]
                };
                for jj in 0..w {
                    fma(av, b.get(c, l, j + jj), &mut sums[jj][..h]);
                }
            }
            for jj in 0..w {
                let ci = offset + i + (j + jj) * ldc;
                for value in &mut sums[jj][..h] {
                    *value *= alpha;
                }
                if beta != 0. {
                    fma(&c[ci..ci + h], beta, &mut sums[jj][..h]);
                }
                c[ci..ci + h].copy_from_slice(&sums[jj][..h]);
            }
        }
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx2,fma")]
unsafe fn native_dot_pair<const H: usize, const W: usize>(
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    ar: usize,
    b: Matrix<'_>,
    bc: usize,
    beta: f64,
    c: &mut [f64],
    co: usize,
    ldc: usize,
    four: bool,
    vector_store: bool,
) {
    // OpenBLAS❗❗: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // OpenBLAS❗❗: #define LOAD_KA_512(M, N) __m512d Aval##M = _mm512_loadu_pd(&A[(i + M)*lda + k]);
    // OpenBLAS❗❗: #define LOAD_KB_512(M, N) __m512d Bval##N = _mm512_loadu_pd(&B[(j + N)*ldb + k])
    // OpenBLAS❗❗: r0 = _mm512_add_pd(t0, t1); r1 = _mm512_add_pd(t2, t3); t0 = _mm512_add_pd(r0, r1);
    // OpenBLAS❗❗: s0 = _mm256_add_pd(s0, s1); s0 = _mm256_mul_pd(alpha_256, s0);
    use std::arch::x86_64::*;
    let mut sums = [[_mm512_setzero_pd(); W]; H];
    let bulk = k / 8 * 8;
    for l in (0..bulk).step_by(8) {
        let mut av = [_mm512_setzero_pd(); H];
        let mut bv = [_mm512_setzero_pd(); W];
        // SAFETY: each source dot input is contiguous along k; the entry
        // checked all active dimensions. These complete eight-lane blocks
        // stay inside their rows. Feature availability is source_simd.
        unsafe {
            for row in 0..H {
                av[row] = _mm512_loadu_pd(
                    a.separate
                        .unwrap_or(c)
                        .as_ptr()
                        .add(a.offset + (ar + row) * a.ld + l),
                );
            }
            for col in 0..W {
                bv[col] = _mm512_loadu_pd(
                    b.separate
                        .unwrap_or(c)
                        .as_ptr()
                        .add(b.offset + (bc + col) * b.ld + l),
                );
            }
        }
        for row in 0..H {
            for col in 0..W {
                sums[row][col] = _mm512_fmadd_pd(av[row], bv[col], sums[row][col]);
            }
        }
    }
    if bulk < k {
        let mask = ((1_u16 << (k - bulk)) - 1) as u8;
        let mut av = [_mm512_setzero_pd(); H];
        let mut bv = [_mm512_setzero_pd(); W];
        // SAFETY: mask loads read only the remaining valid k elements, at a
        // pointer inside each checked input row; other lanes become zero.
        unsafe {
            for row in 0..H {
                av[row] = _mm512_maskz_loadu_pd(
                    mask,
                    a.separate
                        .unwrap_or(c)
                        .as_ptr()
                        .add(a.offset + (ar + row) * a.ld + bulk),
                );
            }
            for col in 0..W {
                bv[col] = _mm512_maskz_loadu_pd(
                    mask,
                    b.separate
                        .unwrap_or(c)
                        .as_ptr()
                        .add(b.offset + (bc + col) * b.ld + bulk),
                );
            }
        }
        for row in 0..H {
            for col in 0..W {
                sums[row][col] = _mm512_fmadd_pd(av[row], bv[col], sums[row][col]);
            }
        }
    }
    for row in 0..H {
        for col in 0..W {
            let mut lanes = [0_f64; 8];
            // SAFETY: the local unaligned store has eight valid f64 lanes.
            unsafe {
                _mm512_storeu_pd(lanes.as_mut_ptr(), sums[row][col]);
            }
            let sum = reduce8(lanes, four);
            let ci = co + row + col * ldc;
            let value = alpha * sum;
            c[ci] = if beta == 0. {
                value
            } else if vector_store {
                beta.mul_add(c[ci], value)
            } else {
                alpha.mul_add(sum, beta * c[ci])
            };
        }
    }
}

fn small_gemm(
    ta: Transpose,
    tb: Transpose,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: Matrix<'_>,
    beta: f64,
    c: &mut [f64],
    offset: usize,
    ldc: usize,
) {
    // OpenBLAS❗✔️: #define DECLARE_RESULT_512(M, N) __m512d result##M##N = _mm512_setzero_pd()
    // OpenBLAS❗✔️: #define LOAD_A_512(M, N) __m512d Aval##M = _mm512_loadu_pd(&A[lda * k + i + (M*8)])
    // OpenBLAS❗✔️: #define MASK_LOAD_A_512(M, N) __m512d Aval##M = _mm512_maskz_loadu_pd(mask, &A[lda * k + i + (M*8)])
    // OpenBLAS❗✔️: #define BROADCAST_LOAD_B_512(M, N) __m512d Bval##N = _mm512_broadcastsd_pd(_mm_load_pd1(&B[k + ldb * (j+N)]))
    // OpenBLAS❗✔️: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: #define STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 			_mm512_storeu_pd(&C[(j+N)*ldc + i + (M*8)], result##M##N)
    // OpenBLAS❗✔️: #define MASK_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 			_mm512_mask_storeu_pd(&C[(j+N)*ldc + i + (M*8)], mask, result##M##N)
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: #define STORE_512(M, N) \
    // OpenBLAS❗✔️: 	result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0": "+v"(result##M##N):"r"(&C[(j+N)*ldc + i + (M*8)]), "v"(beta_512)); \
    // OpenBLAS❗✔️: 	_mm512_storeu_pd(&C[(j+N)*ldc + i + (M*8)], result##M##N)
    // OpenBLAS❗✔️: #define MASK_STORE_512(M, N) \
    // OpenBLAS❗✔️: 	result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0 %{%3%}": "+v"(result##M##N):"r"(&C[(j+N)*ldc + i + (M*8)]), "v"(beta_512), "Yk"(mask)); \
    // OpenBLAS❗✔️: 	_mm512_mask_storeu_pd(&C[(j+N)*ldc + i + (M*8)], mask, result##M##N)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define LOAD_KA_512(M, N) __m512d Aval##M = _mm512_loadu_pd(&mbuf[(mi + M)*K + k]);
    // OpenBLAS❗✔️: #define LOAD_KB_512(M, N) __m512d Bval##N = _mm512_loadu_pd(&B[(j + N)*ldb + k])
    // OpenBLAS❗✔️: #define MASK_LOAD_KA_512(M, N) __m512d Aval##M = _mm512_maskz_loadu_pd(mask, &mbuf[(mi + M)*K + k])
    // OpenBLAS❗✔️: #define MASK_LOAD_KB_512(M, N) __m512d Bval##N = _mm512_maskz_loadu_pd(mask, &B[(j + N)*ldb + k])
    // OpenBLAS❗✔️: #define REDUCE_4(rr0, rr1, rr2, rr3) \
    // OpenBLAS❗✔️: 	__m512d r0, r1, r2, r3, t0, t1, t2, t3;\
    // OpenBLAS❗✔️: 	r0 = _mm512_unpacklo_pd(rr0, rr1); r1 = _mm512_unpackhi_pd(rr0, rr1); \
    // OpenBLAS❗✔️: 	r2 = _mm512_unpacklo_pd(rr2, rr3); r3 = _mm512_unpackhi_pd(rr2, rr3); \
    // OpenBLAS❗✔️: 	t0 = _mm512_permutex2var_pd(r0, idx_lo, r2); t1 = _mm512_permutex2var_pd(r1, idx_lo, r3); \
    // OpenBLAS❗✔️: 	t2 = _mm512_permutex2var_pd(r0, idx_hi, r2); t3 = _mm512_permutex2var_pd(r1, idx_hi, r3); \
    // OpenBLAS❗✔️: 	r0 = _mm512_add_pd(t0, t1); r1 = _mm512_add_pd(t2, t3); t0 = _mm512_add_pd(r0, r1); \
    // OpenBLAS❗✔️: 	__m256d s0, s1; \
    // OpenBLAS❗✔️: 	s0 = _mm512_extractf64x4_pd(t0, 0); s1 = _mm512_extractf64x4_pd(t0, 1); \
    // OpenBLAS❗✔️: 	s0 = _mm256_add_pd(s0, s1); s0 = _mm256_mul_pd(alpha_256, s0);
    // OpenBLAS❗✔️: #define REDUCE_M4(N) REDUCE_4(result0##N, result1##N, result2##N, result3##N)
    // OpenBLAS❗✔️: #define REDUCE_N4(M) REDUCE_4(result##M##0, result##M##1, result##M##2, result##M##3)
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: #define STORE_REDUCE(M, N) C[(j+N)*ldc + i + M] = alpha * _mm512_reduce_add_pd(result##M##N);
    // OpenBLAS❗✔️: #define STORE_REDUCE_M4(N) {\
    // OpenBLAS❗✔️: 	REDUCE_M4(N) \
    // OpenBLAS❗✔️: 	_mm256_storeu_pd(&C[(j + N)*ldc + i], s0); \
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️: #define STORE_REDUCE_N4(M) {\
    // OpenBLAS❗✔️: 	REDUCE_N4(M) \
    // OpenBLAS❗✔️: 	_mm256_i64scatter_pd(&C[j*ldc + i + M], vindex_n, s0, 8); \
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: #define STORE_REDUCE(M, N) C[(j+N)*ldc + i + M] = alpha * _mm512_reduce_add_pd(result##M##N) + beta * C[(j+N)*ldc + i + M];
    // OpenBLAS❗✔️: #define STORE_REDUCE_M4(N) {\
    // OpenBLAS❗✔️: 	REDUCE_M4(N) \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0": "+v"(s0):"r"(&C[(j + N)*ldc + i]), "v"(beta_256)); \
    // OpenBLAS❗✔️: 	_mm256_storeu_pd(&C[(j + N)*ldc + i], s0); \
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️: #define STORE_REDUCE_N4(M) {\
    // OpenBLAS❗✔️: 	REDUCE_N4(M) \
    // OpenBLAS❗✔️: 	s1 = _mm256_i64gather_pd(&C[j*ldc + i + M], vindex_n, 8); \
    // OpenBLAS❗✔️: 	s0 = _mm256_fmadd_pd(s1, beta_256, s0); \
    // OpenBLAS❗✔️: 	_mm256_i64scatter_pd(&C[j*ldc + i + M], vindex_n, s0, 8); \
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT * C, BLASLONG ldc)
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT beta, FLOAT * C, BLASLONG ldc)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: {
    // OpenBLAS❗✔️: 	// column major
    // OpenBLAS❗✔️: 	BLASLONG i, j, k;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	BLASLONG m32 = M & ~31;
    // OpenBLAS❗✔️: 	BLASLONG m16 = M & ~15;
    // OpenBLAS❗✔️: 	BLASLONG m8 = M & ~7;
    // OpenBLAS❗✔️: 	BLASLONG m4 = M & ~3;
    // OpenBLAS❗✔️: 	BLASLONG m2 = M & ~1;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	BLASLONG n6 = N - (N % 6);
    // OpenBLAS❗✔️: 	BLASLONG n4 = N & ~3;
    // OpenBLAS❗✔️: 	BLASLONG n2 = N & ~1;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	__m512d alpha_512 = _mm512_broadcastsd_pd(_mm_load_pd1(&alpha));
    // OpenBLAS❗✔️: #if !defined(B0)
    // OpenBLAS❗✔️: 	__m512d beta_512 = _mm512_broadcastsd_pd(_mm_load_pd1(&beta));
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	for (i = 0; i < m32; i += 32) {
    // OpenBLAS❗✔️: 		for (j = 0; j < n4; j += 4) {
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 1);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 2);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 3);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 4);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 5);
    // OpenBLAS❗✔️: 			for (k = 0; k < K; k++) {
    // OpenBLAS❗✔️: 				MASK_LOAD_A_512(0, x);
    // OpenBLAS❗✔️: 				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
    // OpenBLAS❗✔️: 				BROADCAST_LOAD_B_512(x, 2); BROADCAST_LOAD_B_512(x, 3);
    // OpenBLAS❗✔️: 				BROADCAST_LOAD_B_512(x, 4); BROADCAST_LOAD_B_512(x, 5);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 				MATMUL_512(0, 0);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 1);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 2);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 3);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 4);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 5);
    // OpenBLAS❗✔️: 			}
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 0);
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 1);
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 2);
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 3);
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 4);
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 5);
    // OpenBLAS❗✔️: 		}
    // OpenBLAS❗✔️: 		for (; j < n2; j += 2) {
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 0);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 1);
    // OpenBLAS❗✔️: 			for (k = 0; k < K; k++) {
    // OpenBLAS❗✔️: 				MASK_LOAD_A_512(0, x);
    // OpenBLAS❗✔️: 				BROADCAST_LOAD_B_512(x, 0); BROADCAST_LOAD_B_512(x, 1);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 0);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 1);
    // OpenBLAS❗✔️: 			}
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 0);
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 1);
    // OpenBLAS❗✔️: 		}
    // OpenBLAS❗✔️: 		for (; j < N; j++) {
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 0);
    // OpenBLAS❗✔️: 			for (k = 0; k < K; k++) {
    // OpenBLAS❗✔️: 				MASK_LOAD_A_512(0, x);
    // OpenBLAS❗✔️: 				BROADCAST_LOAD_B_512(x, 0);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 0);
    // OpenBLAS❗✔️: 			}
    // OpenBLAS❗✔️: 			MASK_STORE_512(0, 0);
    // OpenBLAS❗✔️: 		}
    // OpenBLAS❗✔️: 	} else {
    // OpenBLAS❗✔️: 		/* M => [1, 4]
    // OpenBLAS❗✔️: 		 *
    // OpenBLAS❗✔️: 		 * This kernel use dot-like style to calc a value - C(x, y):
    // OpenBLAS❗✔️: 		 * C(x, y) = A(x, 0)*B(0, y) + A(x, 1)*B(1, y) +....+ A(x, K)*B(K, y)
    // OpenBLAS❗✔️: 		 *
    // OpenBLAS❗✔️: 		 * Alloc a buf to copy rest of A as row major,
    // OpenBLAS❗✔️: 		 * so memory access from 0 to K is continuous for both A & B.
    // OpenBLAS❗✔️: 		 *
    // OpenBLAS❗✔️: 		 * Loading to zmm and FMA 8 of k at one loop,
    // OpenBLAS❗✔️: 		 * finally reduce_add zmm to a single float result in C(x, y).
    // OpenBLAS❗✔️: 		 *
    // OpenBLAS❗✔️: 		 * Note: performance is bad when K is small.
    // OpenBLAS❗✔️: 		 */
    // OpenBLAS❗✔️: 		FLOAT *mbuf = (FLOAT *) malloc(sizeof(FLOAT)*mm*K);
    // OpenBLAS❗✔️: 		__mmask8 mask = (1UL << mm) - 1;
    // OpenBLAS❗✔️: 		BLASLONG k8 = K & ~7;
    // OpenBLAS❗✔️: 		BLASLONG k4 = K & ~3;
    // OpenBLAS❗✔️: 		for (k = 0; k < k4; k += 4) {
    // OpenBLAS❗✔️: 			__m256d  r0, r1, r2, r3;
    // OpenBLAS❗✔️: 			__m256d  t0, t1, t2, t3;
    // OpenBLAS❗✔️: 			r0 = _mm256_maskz_loadu_pd(mask, &A[i + lda*(0 + k)]);
    // OpenBLAS❗✔️: 	BLASLONG n2 = N & ~1;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	BLASLONG k8 = K & ~7;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	__mmask8 mask;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	__m256i vindex_n = _mm256_set_epi64x(ldc*3, ldc*2, ldc, 0);
    // OpenBLAS❗✔️: 	__m256d alpha_256 = _mm256_broadcast_sd(&alpha);
    // OpenBLAS❗✔️: #if !defined(B0)
    // OpenBLAS❗✔️: 	__m256d beta_256 = _mm256_broadcast_sd(&beta);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	long long permute_table[] = {
    // OpenBLAS❗✔️: 		0, 1, 0|8, 1|8, 4, 5, 4|8, 5|8,
    // OpenBLAS❗✔️: 		2, 3, 2|8, 3|8, 6, 7, 6|8, 7|8,
    // OpenBLAS❗✔️: 	};
    // OpenBLAS❗✔️: 	__m512i idx_lo = _mm512_loadu_si512(permute_table);
    // OpenBLAS❗✔️: 	__m512i idx_hi = _mm512_loadu_si512(permute_table + 8);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	for (i = 0; i < m4; i += 4) {
    // OpenBLAS❗✔️: 		for (j = 0; j < n4; j += 4) {
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1); DECLARE_RESULT_512(2, 1); DECLARE_RESULT_512(3, 1);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 2); DECLARE_RESULT_512(1, 2); DECLARE_RESULT_512(2, 2); DECLARE_RESULT_512(3, 2);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 3); DECLARE_RESULT_512(1, 3); DECLARE_RESULT_512(2, 3); DECLARE_RESULT_512(3, 3);
    // OpenBLAS❗✔️: 			for (k = 0; k < k8; k += 8) {
    // OpenBLAS❗✔️: 				LOAD_KA_512(0, x); LOAD_KA_512(1, x); LOAD_KA_512(2, x); LOAD_KA_512(3, x);
    // OpenBLAS❗✔️: 				LOAD_KB_512(x, 0); LOAD_KB_512(x, 1); LOAD_KB_512(x, 2); LOAD_KB_512(x, 3);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 				MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 2); MATMUL_512(1, 2); MATMUL_512(2, 2); MATMUL_512(3, 2);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 3); MATMUL_512(1, 3); MATMUL_512(2, 3); MATMUL_512(3, 3);
    // OpenBLAS❗✔️: 			}
    // OpenBLAS❗✔️: 			int remains = K - k;
    // OpenBLAS❗✔️: 			if (remains) {
    // OpenBLAS❗✔️: 				mask = (1UL << remains) - 1;
    // OpenBLAS❗✔️: 				MASK_LOAD_KA_512(0, x); MASK_LOAD_KA_512(1, x); MASK_LOAD_KA_512(2, x); MASK_LOAD_KA_512(3, x);
    // OpenBLAS❗✔️: 				MASK_LOAD_KB_512(x, 0); MASK_LOAD_KB_512(x, 1); MASK_LOAD_KB_512(x, 2); MASK_LOAD_KB_512(x, 3);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 				MATMUL_512(0, 0); MATMUL_512(1, 0); MATMUL_512(2, 0); MATMUL_512(3, 0);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 1); MATMUL_512(1, 1); MATMUL_512(2, 1); MATMUL_512(3, 1);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 2); MATMUL_512(1, 2); MATMUL_512(2, 2); MATMUL_512(3, 2);
    // OpenBLAS❗✔️: 				MATMUL_512(0, 3); MATMUL_512(1, 3); MATMUL_512(2, 3); MATMUL_512(3, 3);
    // OpenBLAS❗✔️: 			}
    // OpenBLAS❗✔️: 			STORE_REDUCE_M4(0); STORE_REDUCE_M4(1); STORE_REDUCE_M4(2); STORE_REDUCE_M4(3);
    // OpenBLAS❗✔️: 		}
    // OpenBLAS❗✔️: 		for (; j < n2; j += 2) {
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 0); DECLARE_RESULT_512(1, 0); DECLARE_RESULT_512(2, 0); DECLARE_RESULT_512(3, 0);
    // OpenBLAS❗✔️: 			DECLARE_RESULT_512(0, 1); DECLARE_RESULT_512(1, 1); DECLARE_RESULT_512(2, 1); DECLARE_RESULT_512(3, 1);
    // OpenBLAS❗✔️: 			for (k = 0; k < k8; k += 8) {
    // OpenBLAS❗✔️: 				LOAD_KA_512(0, x); LOAD_KA_512(1, x); LOAD_KA_512(2, x); LOAD_KA_512(3, x);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define DECLARE_RESULT_512(M, N) __m512d result##M##N = _mm512_setzero_pd()
    // OpenBLAS❗✔️: #define LOAD_A_512(M, N) __m512d Aval##M = _mm512_loadu_pd(&A[lda * k + i + (M*8)])
    // OpenBLAS❗✔️: #define MASK_LOAD_A_512(M, N) __m512d Aval##M = _mm512_maskz_loadu_pd(mask, &A[lda * k + i + (M*8)])
    // OpenBLAS❗✔️: #define BROADCAST_LOAD_B_512(M, N) __m512d Bval##N = _mm512_broadcastsd_pd(_mm_load_sd(&B[ldb * k + j + N]))
    // OpenBLAS❗✔️: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define BROADCAST_LOAD_A_512(M, N) __m512d Aval##M = _mm512_broadcastsd_pd(_mm_load_sd(&A[lda * k + i + M]))
    // OpenBLAS❗✔️: #define LOAD_B_512(M, N) __m512d Bval##N = _mm512_loadu_pd(&B[ldb * k + j + (N*8)])
    // OpenBLAS❗✔️: #define MASK_LOAD_B_512(M, N) __m512d Bval##N = _mm512_maskz_loadu_pd(mask, &B[ldb * k + j + (N*8)])
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: #define STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 			_mm512_storeu_pd(&C[(j+N)*ldc + i + (M*8)], result##M##N)
    // OpenBLAS❗✔️: #define MASK_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 			_mm512_mask_storeu_pd(&C[(j+N)*ldc + i + (M*8)], mask, result##M##N)
    // OpenBLAS❗✔️: #define SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				_mm512_i64scatter_pd(&C[(j + N*8)*ldc + i + M], vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #define MASK_SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				_mm512_mask_i64scatter_pd(&C[(j + N*8)*ldc + i + M], mask, vindex_n, result##M##N, 8)
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: #define STORE_512(M, N) \
    // OpenBLAS❗✔️: 	result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0": "+v"(result##M##N):"r"(&C[(j+N)*ldc + i + (M*8)]), "v"(beta_512)); \
    // OpenBLAS❗✔️: 	_mm512_storeu_pd(&C[(j+N)*ldc + i + (M*8)], result##M##N)
    // OpenBLAS❗✔️: #define MASK_STORE_512(M, N) \
    // OpenBLAS❗✔️: 	result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0 %{%3%}": "+v"(result##M##N):"r"(&C[(j+N)*ldc + i + (M*8)]), "v"(beta_512), "Yk"(mask)); \
    // OpenBLAS❗✔️: 	_mm512_mask_storeu_pd(&C[(j+N)*ldc + i + (M*8)], mask, result##M##N)
    // OpenBLAS❗✔️: #define SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				__m512d tmp##M##N = _mm512_i64gather_pd(vindex_n, &C[(j + N*8)*ldc + i + M], 8); \
    // OpenBLAS❗✔️: 				result##M##N = _mm512_fmadd_pd(tmp##M##N, beta_512, result##M##N); \
    // OpenBLAS❗✔️: 				_mm512_i64scatter_pd(&C[(j + N*8)*ldc + i + M], vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #define MASK_SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				__m512d tmp##M##N = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), mask, vindex_n, &C[(j + N*8)*ldc + i + M], 8); \
    // OpenBLAS❗✔️: 				result##M##N = _mm512_fmadd_pd(tmp##M##N, beta_512, result##M##N); \
    // OpenBLAS❗✔️: 				_mm512_mask_i64scatter_pd(&C[(j + N*8)*ldc + i + M], mask, vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT * C, BLASLONG ldc)
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: int CNAME(BLASLONG M, BLASLONG N, BLASLONG K, FLOAT * A, BLASLONG lda, FLOAT alpha, FLOAT * B, BLASLONG ldb, FLOAT beta, FLOAT * C, BLASLONG ldc)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: {
    // OpenBLAS❗✔️: 	// column major
    // OpenBLAS❗✔️: 	BLASLONG i, j, k;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	BLASLONG m32 = M & ~31;
    // OpenBLAS❗✔️: 	BLASLONG m16 = M & ~15;
    // OpenBLAS❗✔️: 	BLASLONG m8 = M & ~7;
    // OpenBLAS❗✔️: 	BLASLONG m4 = M & ~3;
    // OpenBLAS❗✔️: 	BLASLONG m2 = M & ~1;
    // OpenBLAS❗✔️: #define DECLARE_RESULT_512(M, N) __m512d result##M##N = _mm512_setzero_pd()
    // OpenBLAS❗✔️: #define BROADCAST_LOAD_A_512(M, N) __m512d Aval##M = _mm512_broadcastsd_pd(_mm_load_sd(&A[k  + lda * (i+M)]))
    // OpenBLAS❗✔️: #define LOAD_B_512(M,N)  __m512d Bval##N = _mm512_loadu_pd(&B[ldb * k + j + (N*8)])
    // OpenBLAS❗✔️: #define MASK_LOAD_B_512(M, N) __m512d Bval##N = _mm512_maskz_loadu_pd(mask, &B[ldb * k + j + (N*8)])
    // OpenBLAS❗✔️: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(B0)
    // OpenBLAS❗✔️: #define STORE_8xy(v, N, x, y) _mm512_storeu_pd(&C[(j + N*8 + x + y*8)*ldc + i], v)
    // OpenBLAS❗✔️: #define STORE_4xy(v, N, x, y) _mm256_storeu_pd(&C[(j + N*8 + x + y*4)*ldc + i], v)
    // OpenBLAS❗✔️: #define SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				_mm512_i64scatter_pd(&C[(j + N*8)*ldc + i + M], vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #define MASK_SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				    _mm512_mask_i64scatter_pd(&C[(j + N*8)*ldc + i + M], mask, vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: #define STORE_8xy(v, N, x, y) \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0": "+v"(v): "r"(&C[(j + N*8 + x + y*8)*ldc + i]), "v"(beta_512)); \
    // OpenBLAS❗✔️: 	_mm512_storeu_pd(&C[(j + N*8 + x + y*8)*ldc + i], v)
    // OpenBLAS❗✔️: #define STORE_4xy(v, N, x, y) \
    // OpenBLAS❗✔️: 	asm("vfmadd231pd (%1), %2, %0": "+v"(v): "r"(&C[(j + N*8 + x + y*4)*ldc + i]), "v"(beta_256)); \
    // OpenBLAS❗✔️: 	_mm256_storeu_pd(&C[(j + N*8 + x + y*4)*ldc + i], v)
    // OpenBLAS❗✔️: #define SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				__m512d tmp##M##N = _mm512_i64gather_pd(vindex_n, &C[(j + N*8)*ldc + i + M], 8); \
    // OpenBLAS❗✔️: 				result##M##N = _mm512_fmadd_pd(tmp##M##N, beta_512, result##M##N); \
    // OpenBLAS❗✔️: 				_mm512_i64scatter_pd(&C[(j + N*8)*ldc + i + M], vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #define MASK_SCATTER_STORE_512(M, N) result##M##N = _mm512_mul_pd(result##M##N, alpha_512); \
    // OpenBLAS❗✔️: 				__m512d tmp##M##N = _mm512_mask_i64gather_pd(_mm512_setzero_pd(), mask, vindex_n, &C[(j + N*8)*ldc + i + M], 8); \
    // OpenBLAS❗✔️: 				result##M##N = _mm512_fmadd_pd(tmp##M##N, beta_512, result##M##N); \
    // OpenBLAS❗✔️: 				_mm512_mask_i64scatter_pd(&C[(j + N*8)*ldc + i + M], mask, vindex_n, result##M##N, 8);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #define REORDER_8x8(r0, r1, r2, r3, r4, r5, r6, r7) \
    // OpenBLAS❗✔️: 	__m512d t0, t1, t2, t3, t4, t5, t6, t7; \
    // OpenBLAS❗✔️: 	t0 = _mm512_unpacklo_pd(r0, r1); \
    // OpenBLAS❗✔️: 	t1 = _mm512_unpackhi_pd(r0, r1); \
    // OpenBLAS❗✔️: 	t2 = _mm512_unpacklo_pd(r2, r3); \
    // OpenBLAS❗✔️: 	t3 = _mm512_unpackhi_pd(r2, r3); \
    // OpenBLAS❗✔️: 	t4 = _mm512_unpacklo_pd(r4, r5); \
    // OpenBLAS❗✔️: 	t5 = _mm512_unpackhi_pd(r4, r5); \
    // OpenBLAS❗✔️: 	t6 = _mm512_unpacklo_pd(r6, r7); \
    let direct_rows = if ta == Transpose::NoTrans && tb == Transpose::NoTrans {
        m / 8 * 8
    } else {
        m
    };
    let tail_dot = ta == Transpose::Trans && tb == Transpose::NoTrans;
    if !tail_dot {
        small_direct(m, n, k, alpha, a, b, beta, c, offset, ldc, 0, direct_rows);
    }
    let start = if tail_dot { 0 } else { direct_rows };
    let mm = m - start;
    // OpenBLAS❗✔️: if (mm > 4 || K < 16) {
    if !tail_dot && (mm > 4 || k < 16) {
        small_direct(m, n, k, alpha, a, b, beta, c, offset, ldc, start, m);
        return;
    }
    // NN's last <=4 rows are transposed into its source row-major mbuf.
    let mut mbuf = if !tail_dot {
        vec![0.; mm * k]
    } else {
        Vec::new()
    };
    if !tail_dot {
        for i in 0..mm {
            for l in 0..k {
                mbuf[i * k + l] = a.get(c, start + i, l);
            }
        }
    }
    #[cfg(target_arch = "x86_64")]
    if source_simd() {
        let dot_a = if tail_dot {
            a
        } else {
            Matrix {
                separate: Some(&mbuf),
                offset: 0,
                ld: k,
                trans: Transpose::Trans,
            }
        };
        let mut i = start;
        while i < m {
            let h = (m - i).min(4);
            let h = if h >= 4 {
                4
            } else if h >= 2 {
                2
            } else {
                1
            };
            let mut j = 0;
            while j < n {
                let w = if n - j >= 4 {
                    4
                } else if n - j >= 2 {
                    2
                } else {
                    1
                };
                let four = h == 4 || j < n / 4 * 4;
                let ar = if tail_dot { i } else { i - start };
                let ci = offset + i + j * ldc;
                // SAFETY: guarded source eight-lane dot panels and checked
                // independent output storage; these are the source's 9 tiles.
                macro_rules! run {
                    ($h:expr,$w:expr) => {
                        unsafe {
                            native_dot_pair::<$h, $w>(
                                k,
                                alpha,
                                dot_a,
                                ar,
                                b,
                                j,
                                beta,
                                c,
                                ci,
                                ldc,
                                four,
                                h > 1 || four,
                            )
                        }
                    };
                }
                match (h, w) {
                    (4, 4) => run!(4, 4),
                    (4, 2) => run!(4, 2),
                    (4, 1) => run!(4, 1),
                    (2, 4) => run!(2, 4),
                    (2, 2) => run!(2, 2),
                    (2, 1) => run!(2, 1),
                    (1, 4) => run!(1, 4),
                    (1, 2) => run!(1, 2),
                    _ => run!(1, 1),
                }
                j += w;
            }
            i += h;
        }
        return;
    }
    let mut i = start;
    while i < m {
        let h = (m - i).min(4);
        let h = if h >= 4 {
            4
        } else if h >= 2 {
            2
        } else {
            1
        };
        for j in 0..n {
            let four = h == 4 || j < n / 4 * 4;
            for ii in 0..h {
                let mut sums = [0_f64; 8];
                for l in 0..k {
                    let av = if tail_dot {
                        a.get(c, i + ii, l)
                    } else {
                        mbuf[(i + ii - start) * k + l]
                    };
                    sums[l % 8] = av.mul_add(b.get(c, l, j), sums[l % 8]);
                }
                // The final masked FMA uses zeroed inactive lanes, even if a
                // lane already contains signed zero.
                if k % 8 != 0 {
                    for lane in k % 8..8 {
                        sums[lane] = 0_f64.mul_add(0., sums[lane]);
                    }
                }
                let sum = reduce8(sums, four);
                let ci = offset + i + ii + j * ldc;
                let value = alpha * sum;
                c[ci] = if beta == 0. {
                    value
                } else {
                    if four || h > 1 {
                        beta.mul_add(c[ci], value)
                    } else {
                        alpha.mul_add(sum, beta * c[ci])
                    }
                };
            }
        }
        i += h;
    }
}
fn gemm_inner(
    ta: Transpose,
    tb: Transpose,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: Matrix<'_>,
    beta: f64,
    c: &mut [f64],
    offset: usize,
    ldc: usize,
) -> Result<(), Blas3Error> {
    // OpenBLAS❗✔️: {
    // OpenBLAS❗✔️: 	double MNK = (double) M * (double) N * (double) K;
    // OpenBLAS❗✔️: 	if (MNK > 100.0*100.0*100.0)  // disable for big size matrix
    // OpenBLAS❗✔️: 		return 0;
    // OpenBLAS❗✔️: 	if (transa && !transb) {
    // OpenBLAS❗✔️: 		/* TN kernel perform not good when:
    // OpenBLAS❗✔️: 		 * 1. C matrix is too big
    // OpenBLAS❗✔️: 		 * 2. K is too small
    // OpenBLAS❗✔️: 		 */
    // OpenBLAS❗✔️: 		if (M * N > 1200 || K < 32)
    // OpenBLAS❗✔️: 			return 0;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️: 	return 1;
    // OpenBLAS❗✔️: }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(GEMM_GEMV_FORWARD) && !defined(GEMM3M) && !defined(COMPLEX) && HFLOAT16_GEMM_GEMV_FORWARD && BFLOAT16_GEMM_GEMV_FORWARD
    // OpenBLAS❗✔️: #if defined(ARCH_ARM64)
    // OpenBLAS❗✔️:   // The gemv kernels in arm64/{gemv_n.S,gemv_n_sve.c,gemv_t.S,gemv_t_sve.c}
    // OpenBLAS❗✔️:   // perform poorly in certain circumstances. We use the following boolean
    // OpenBLAS❗✔️:   // variable along with the gemv argument values to avoid these inefficient
    // OpenBLAS❗✔️:   // gemv cases, see github issue#4951.
    // OpenBLAS❗✔️:   bool have_tuned_gemv = false;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   bool have_tuned_gemv = true;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:   // Check if we can convert GEMM -> GEMV
    // OpenBLAS❗✔️:   if (args.k != 0) {
    // OpenBLAS❗✔️:     if (args.n == 1) {
    // OpenBLAS❗✔️:       blasint inc_x = 1;
    // OpenBLAS❗✔️:       blasint inc_y = 1;
    // OpenBLAS❗✔️:       // These were passed in as blasint, but the struct translates them to blaslong
    // OpenBLAS❗✔️:       blasint m = args.m;
    // OpenBLAS❗✔️:       blasint n = args.k;
    // OpenBLAS❗✔️:       blasint lda = args.lda;
    // OpenBLAS❗✔️:       // Create new transpose parameters
    // OpenBLAS❗✔️:       char NT = 'N';
    // OpenBLAS❗✔️:       if (transa & 1) {
    // OpenBLAS❗✔️:         NT = 'T';
    // OpenBLAS❗✔️:         m = args.k;
    // OpenBLAS❗✔️:         n = args.m;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       if (transb & 1) {
    // OpenBLAS❗✔️:         inc_x = args.ldb;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       bool is_efficient_gemv = have_tuned_gemv || ((NT == 'N') || (NT == 'T' && inc_x == 1));
    // OpenBLAS❗✔️:       if (is_efficient_gemv) {
    // OpenBLAS❗✔️:         GEMV(&NT, &m, &n, args.alpha, args.a, &lda, args.b, &inc_x, args.beta, args.c, &inc_y);
    // OpenBLAS❗✔️:         return;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:     if (args.m == 1) {
    // OpenBLAS❗✔️:       blasint inc_x = args.lda;
    // OpenBLAS❗✔️:       blasint inc_y = args.ldc;
    // OpenBLAS❗✔️:       // These were passed in as blasint, but the struct translates them to blaslong
    // OpenBLAS❗✔️:       blasint m = args.k;
    // OpenBLAS❗✔️:       blasint n = args.n;
    // OpenBLAS❗✔️:       blasint ldb = args.ldb;
    // OpenBLAS❗✔️:       // Create new transpose parameters
    // OpenBLAS❗✔️:       char NT = 'T';
    // OpenBLAS❗✔️:       if (transa & 1) {
    // OpenBLAS❗✔️:         inc_x = 1;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       if (transb & 1) {
    // OpenBLAS❗✔️:         NT = 'N';
    // OpenBLAS❗✔️:         m = args.n;
    // OpenBLAS❗✔️:         n = args.k;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:       bool is_efficient_gemv = have_tuned_gemv || ((NT == 'N' && inc_y == 1) || (NT == 'T' && inc_x == 1));
    // OpenBLAS❗✔️:       if (is_efficient_gemv) {
    // OpenBLAS❗✔️:         GEMV(&NT, &m, &n, args.alpha, args.b, &ldb, args.a, &inc_x, args.beta, args.c, &inc_y);
    // OpenBLAS❗✔️:         return;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   IDEBUG_START;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   FUNCTION_PROFILE_START();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if USE_SMALL_MATRIX_OPT
    // OpenBLAS❗✔️: #if !defined(COMPLEX)
    // OpenBLAS❗✔️:   if(GEMM_SMALL_MATRIX_PERMIT(transa, transb, args.m, args.n, args.k, *(FLOAT *)(args.alpha), *(FLOAT *)(args.beta))){
    // OpenBLAS❗✔️: 	  if(*(FLOAT *)(args.beta) == 0.0){
    // OpenBLAS❗✔️: 		(GEMM_SMALL_KERNEL_B0((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, args.c, args.ldc);
    // OpenBLAS❗✔️: 	  }else{
    // OpenBLAS❗✔️: 		(GEMM_SMALL_KERNEL((transb << 2) | transa))(args.m, args.n, args.k, args.a, args.lda, *(FLOAT *)(args.alpha), args.b, args.ldb, *(FLOAT *)(args.beta), args.c, args.ldc);
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:   if (range_n) {
    // OpenBLAS❗✔️:     n_from = *(((BLASLONG *)range_n) + 0);
    // OpenBLAS❗✔️:     n_to   = *(((BLASLONG *)range_n) + 1);
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (beta) {
    // OpenBLAS❗✔️: #if !defined(XDOUBLE) || !defined(QUAD_PRECISION)
    // OpenBLAS❗✔️: #ifndef COMPLEX
    // OpenBLAS❗✔️:     if (beta[0] != ONE
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     if ((beta[0] != ONE) || (beta[1] != ZERO)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     if (((beta[0].x[1] != 0x3fff000000000000UL) || beta[0].x[0] != 0)
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 	&&(((beta[1].x[0] | beta[1].x[1]) << 1) != 0)
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	) {
    // OpenBLAS❗✔️: #if defined(XDOUBLE) && defined(QUAD_PRECISION)
    // OpenBLAS❗✔️: 	  xidouble xbeta;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  qtox(&xbeta, beta);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	  BETA_OPERATION(m_from, m_to, n_from, n_to, beta, c, ldc);
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if ((k == 0) || (alpha == NULL)) return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if !defined(XDOUBLE) || !defined(QUAD_PRECISION)
    // OpenBLAS❗✔️:   if ( alpha[0] == ZERO
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️:       && alpha[1] == ZERO
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	 ) return 0;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   if (((alpha[0].x[0] | alpha[0].x[1]
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️:        | alpha[1].x[0] | alpha[1].x[1]
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:        ) << 1) == 0) return 0;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(XDOUBLE)  && defined(QUAD_PRECISION)
    // OpenBLAS❗✔️:   qtox(&xalpha, alpha);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   l2size = GEMM_P * GEMM_Q;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if 0
    // OpenBLAS❗✔️:   fprintf(stderr, "GEMM(Single): M_from : %ld  M_to : %ld  N_from : %ld  N_to : %ld  k : %ld\n", m_from, m_to, n_from, n_to, k);
    // OpenBLAS❗✔️:   fprintf(stderr, "GEMM(Single):: P = %4ld  Q = %4ld  R = %4ld\n", (BLASLONG)GEMM_P, (BLASLONG)GEMM_Q, (BLASLONG)GEMM_R);
    // OpenBLAS❗✔️: 	//  fprintf(stderr, "GEMM: SA .. %p  SB .. %p\n", sa, sb);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	//  fprintf(stderr, "A = %p  B = %p  C = %p\n\tlda = %ld  ldb = %ld ldc = %ld\n", a, b, c, lda, ldb, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifdef TIMING
    // OpenBLAS❗✔️:   innercost = 0;
    // OpenBLAS❗✔️:   outercost = 0;
    // OpenBLAS❗✔️:   kernelcost = 0;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   for(js = n_from; js < n_to; js += GEMM_R){
    // OpenBLAS❗✔️:     if (openblas_cancel_poll(cancel_slot, cancel_gen)) return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     min_j = n_to - js;
    // OpenBLAS❗✔️:     if (min_j > GEMM_R) min_j = GEMM_R;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = 0; ls < k; ls += min_l){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (openblas_cancel_poll(cancel_slot, cancel_gen)) return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       min_l = k - ls;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (min_l >= GEMM_Q * 2) {
    // OpenBLAS❗✔️: 	// gemm_p = GEMM_P;
    // OpenBLAS❗✔️: 	min_l  = GEMM_Q;
    // OpenBLAS❗✔️:       } else {
    // OpenBLAS❗✔️: 	if (min_l > GEMM_Q) {
    // OpenBLAS❗✔️: 	  min_l = ((min_l / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️: 	gemm_p = ((l2size / min_l + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	while (gemm_p * min_l > l2size) gemm_p -= GEMM_UNROLL_M;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       BLASLONG pad_min_l = min_l;
    // OpenBLAS❗✔️: #if defined(BFLOAT16)
    // OpenBLAS❗✔️:     pad_min_l = (min_l + BFLOAT16_ALIGN_K - 1) & ~(BFLOAT16_ALIGN_K - 1);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       /* First, we have to move data A to L2 cache */
    // OpenBLAS❗✔️:       min_i = m_to - m_from;
    // OpenBLAS❗✔️:       l1stride = 1;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (min_i >= GEMM_P * 2) {
    // OpenBLAS❗✔️: 	min_i = GEMM_P;
    // OpenBLAS❗✔️:       } else {
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) {
    // OpenBLAS❗✔️: 	  min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	} else {
    // OpenBLAS❗✔️: 	  l1stride = 0;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_from, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if defined(FUSED_GEMM) && !defined(TIMING)
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       FUSED_KERNEL_OPERATION(min_i, min_j, min_l, alpha,
    // OpenBLAS❗✔️: 			     sa, sb, b, ldb, c, ldc, m_from, js, ls);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️: 	if (openblas_cancel_poll(cancel_slot, cancel_gen)) return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:         if (min_jj >= 3*GEMM_UNROLL_N) min_jj = 3*GEMM_UNROLL_N;
    // OpenBLAS❗✔️:         else
    // OpenBLAS❗✔️: /*
    // OpenBLAS❗✔️: 		if (min_jj >= 2*GEMM_UNROLL_N) min_jj = 2*GEMM_UNROLL_N;
    // OpenBLAS❗✔️:         	else
    // OpenBLAS❗✔️: */
    // OpenBLAS❗✔️:           		if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	OCOPY_OPERATION(min_l, min_jj, b, ldb, ls, jjs,
    // OpenBLAS❗✔️: 			sb + pad_min_l * (jjs - js) * COMPSIZE * l1stride);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(outercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if !defined(XDOUBLE)  || !defined(QUAD_PRECISION)
    // OpenBLAS❗✔️: 	KERNEL_OPERATION(min_i, min_jj, min_l, alpha,
    // OpenBLAS❗✔️: 			 sa, sb + pad_min_l * (jjs - js)  * COMPSIZE * l1stride, c, ldc, m_from, jjs);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	KERNEL_OPERATION(min_i, min_jj, min_l, (void *)&xalpha,
    // OpenBLAS❗✔️: 			 sa, sb + pad_min_l * (jjs - js)  * COMPSIZE * l1stride, c, ldc, m_from, jjs);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(kernelcost);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = m_from + min_i; is < m_to; is += min_i){
    // OpenBLAS❗✔️: 	if (openblas_cancel_poll(cancel_slot, cancel_gen)) return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	min_i = m_to - is;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (min_i >= GEMM_P * 2) {
    // OpenBLAS❗✔️: 	  min_i = GEMM_P;
    // OpenBLAS❗✔️: 	} else
    // OpenBLAS❗✔️: 	  if (min_i > GEMM_P) {
    // OpenBLAS❗✔️: 	    min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	ICOPY_OPERATION(min_l, min_i, a, lda, ls, is, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if !defined(XDOUBLE)  || !defined(QUAD_PRECISION)
    // OpenBLAS❗✔️: 	KERNEL_OPERATION(min_i, min_j, min_l, alpha, sa, sb, c, ldc, is, js);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	KERNEL_OPERATION(min_i, min_j, min_l, (void *)&xalpha, sa, sb, c, ldc, is, js);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    if m == 0 || n == 0 {
        leading(
            "DGEMM",
            if ta == Transpose::NoTrans { m } else { k },
            a.ld,
            8,
        )?;
        leading(
            "DGEMM",
            if tb == Transpose::NoTrans { k } else { n },
            b.ld,
            10,
        )?;
        leading("DGEMM", m, ldc, 13)?;
        return Ok(());
    }
    op_storage("DGEMM", "A", m, k, a, c.len(), 8)?;
    op_storage("DGEMM", "B", k, n, b, c.len(), 10)?;
    storage("DGEMM", "C", m, n, ldc, offset, c.len(), 13)?;
    independent(m, k, a, m, n, offset, ldc, "DGEMM", 7)?;
    independent(k, n, b, m, n, offset, ldc, "DGEMM", 9)?;
    if m == 0 || n == 0 {
        return Ok(());
    }
    // Makefile.system defines GEMM_GEMV_FORWARD for arm64/riscv64/power,
    // not the frozen x86_64 wheel: M=1/N=1 still use its small GEMM kernels.
    if (m as f64) * (n as f64) * (k as f64) <= 1_000_000.
        && !(ta == Transpose::Trans && tb == Transpose::NoTrans && (m * n > 1200 || k < 32))
    {
        small_gemm(ta, tb, m, n, k, alpha, a, b, beta, c, offset, ldc);
        return Ok(());
    }
    scale(m, n, beta, c, offset, ldc);
    if k == 0 || alpha == 0. {
        return Ok(());
    }
    let mut sa = vec![0.; m.min(P) * k.min(Q)];
    let mut sb = vec![0.; k.min(Q) * if m > P { n.min(R) } else { n.min(6) }];
    for js in (0..n).step_by(R) {
        let nj = (n - js).min(R);
        let mut ls = 0;
        while ls < k {
            let remaining = k - ls;
            let kl = if remaining >= Q * 2 {
                Q
            } else if remaining > Q {
                (remaining / 2).div_ceil(16) * 16
            } else {
                remaining
            };
            let mi = panel_size(m, 16);
            let keep_b = m > P;
            pack_a(mi, kl, a, c, 0, ls, &mut sa[..mi * kl]);
            let mut j = 0;
            while j < nj {
                let w = driver_columns(nj - j);
                let bp = if keep_b { j * kl } else { 0 };
                pack_b(kl, w, b, c, ls, js + j, &mut sb[bp..bp + kl * w]);
                packed_kernel(
                    mi,
                    w,
                    kl,
                    alpha,
                    &sa[..mi * kl],
                    &sb[bp..bp + kl * w],
                    c,
                    offset + (js + j) * ldc,
                    ldc,
                );
                j += w;
            }
            let mut i = mi;
            while i < m {
                let h = panel_size(m - i, 16);
                pack_a(h, kl, a, c, i, ls, &mut sa[..h * kl]);
                packed_kernel(
                    h,
                    nj,
                    kl,
                    alpha,
                    &sa[..h * kl],
                    &sb[..nj * kl],
                    c,
                    offset + i + js * ldc,
                    ldc,
                );
                i += h;
            }
            ls += kl;
        }
    }
    Ok(())
}
#[allow(clippy::too_many_arguments)]
pub(super) fn gemm(
    ta: Transpose,
    tb: Transpose,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    beta: f64,
    c: &mut [f64],
    ldc: usize,
) -> Result<(), Blas3Error> {
    gemm_inner(
        ta,
        tb,
        m,
        n,
        k,
        alpha,
        Matrix {
            separate: Some(a),
            offset: 0,
            ld: lda,
            trans: ta,
        },
        Matrix {
            separate: Some(b),
            offset: 0,
            ld: ldb,
            trans: tb,
        },
        beta,
        c,
        0,
        ldc,
    )
}
#[allow(clippy::too_many_arguments)]
pub(super) fn gemm_within(
    ta: Transpose,
    tb: Transpose,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    data: &mut [f64],
    a_offset: usize,
    lda: usize,
    b_offset: usize,
    ldb: usize,
    beta: f64,
    c_offset: usize,
    ldc: usize,
) -> Result<(), Blas3Error> {
    gemm_inner(
        ta,
        tb,
        m,
        n,
        k,
        alpha,
        Matrix {
            separate: None,
            offset: a_offset,
            ld: lda,
            trans: ta,
        },
        Matrix {
            separate: None,
            offset: b_offset,
            ld: ldb,
            trans: tb,
        },
        beta,
        data,
        c_offset,
        ldc,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn gemm_into_within(
    ta: Transpose,
    tb: Transpose,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    data: &mut [f64],
    a: MatrixSource<'_>,
    lda: usize,
    b: MatrixSource<'_>,
    ldb: usize,
    beta: f64,
    c_offset: usize,
    ldc: usize,
) -> Result<(), Blas3Error> {
    gemm_inner(
        ta,
        tb,
        m,
        n,
        k,
        alpha,
        a.matrix(lda, ta),
        b.matrix(ldb, tb),
        beta,
        data,
        c_offset,
        ldc,
    )
}

fn syr_kernel(
    mut m: usize,
    mut n: usize,
    k: usize,
    alpha: f64,
    mut a: &[f64],
    mut b: &[f64],
    c: &mut [f64],
    mut ci: usize,
    ldc: usize,
    mut offset: isize,
    flag: bool,
) {
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	   FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG offset, int flag){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   BLASLONG i, j;
    // OpenBLAS❗✔️:   BLASLONG loop;
    // OpenBLAS❗✔️:   FLOAT subbuffer[GEMM_UNROLL_MN * GEMM_UNROLL_MN * COMPSIZE];
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (m + offset < 0) {
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(m, n, k,
    // OpenBLAS❗✔️: 		alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		a, b, c, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (n < offset) {
    // OpenBLAS❗✔️: #ifdef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(m, n, k,
    // OpenBLAS❗✔️: 		alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		a, b, c, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (offset > 0) {
    // OpenBLAS❗✔️: #ifdef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(m, offset, k,
    // OpenBLAS❗✔️: 		alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		a, b, c, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     b += offset * k   * COMPSIZE;
    // OpenBLAS❗✔️:     c += offset * ldc * COMPSIZE;
    // OpenBLAS❗✔️:     n -= offset;
    // OpenBLAS❗✔️:     offset = 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     if (n <= 0) return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (n > m + offset) {
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:       GEMM_KERNEL_N(m, n - m - offset, k,
    // OpenBLAS❗✔️: 		  alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		  alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		  a,
    // OpenBLAS❗✔️: 		  b + (m + offset) * k   * COMPSIZE,
    // OpenBLAS❗✔️: 		  c + (m + offset) * ldc * COMPSIZE, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     n = m + offset;
    // OpenBLAS❗✔️:     if (n <= 0) return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (offset < 0) {
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(-offset, n, k,
    // OpenBLAS❗✔️: 		alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		a, b, c, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     a -= offset * k   * COMPSIZE;
    // OpenBLAS❗✔️:     c -= offset       * COMPSIZE;
    // OpenBLAS❗✔️:     m += offset;
    // OpenBLAS❗✔️:     offset = 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (m <= 0) return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (m > n - offset) {
    // OpenBLAS❗✔️: #ifdef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(m - n + offset, n, k,
    // OpenBLAS❗✔️: 		alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		a + (n - offset) * k * COMPSIZE,
    // OpenBLAS❗✔️: 		b,
    // OpenBLAS❗✔️: 		c + (n - offset)     * COMPSIZE, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     m = n + offset;
    // OpenBLAS❗✔️:   if (m <= 0) return 0;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   for (loop = 0; loop < n; loop += GEMM_UNROLL_MN) {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     int mm, nn;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     mm = (loop & ~(GEMM_UNROLL_MN - 1));
    // OpenBLAS❗✔️:     nn = MIN(GEMM_UNROLL_MN, n - loop);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(mm, nn, k,
    // OpenBLAS❗✔️: 		  alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		  alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		  a, b + loop * k * COMPSIZE, c + loop * ldc * COMPSIZE, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     if (flag) {
    // OpenBLAS❗✔️:       GEMM_BETA(nn, nn, 0, ZERO,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		NULL, 0, NULL, 0, subbuffer, nn);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_KERNEL_N(nn, nn, k,
    // OpenBLAS❗✔️: 		    alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    a + loop * k * COMPSIZE, b + loop * k * COMPSIZE, subbuffer, nn);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for (j = 0; j < nn; j ++) {
    // OpenBLAS❗✔️: 	for (i = 0; i <= j; i ++) {
    // OpenBLAS❗✔️: #ifndef COMPLEX
    // OpenBLAS❗✔️: 	  c[i + loop + (j + loop) * ldc] +=
    // OpenBLAS❗✔️: 	    subbuffer[i + j * nn] + subbuffer[j + i * nn];
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	  c[(i + loop + (j + loop) * ldc) * 2 + 0] +=
    // OpenBLAS❗✔️: 	    subbuffer[(i + j * nn) * 2 + 0] + subbuffer[(j + i * nn) * 2 + 0];
    // OpenBLAS❗✔️: 	  c[(i + loop + (j + loop) * ldc) * 2 + 1] +=
    // OpenBLAS❗✔️: 	    subbuffer[(i + j * nn) * 2 + 1] + subbuffer[(j + i * nn) * 2 + 1];
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       for (j = 0; j < nn; j ++) {
    // OpenBLAS❗✔️: 	for (i = j; i < nn; i ++) {
    // OpenBLAS❗✔️: #ifndef COMPLEX
    // OpenBLAS❗✔️: 	  c[i + loop + (j + loop) * ldc] +=
    // OpenBLAS❗✔️: 	    subbuffer[i + j * nn] + subbuffer[j + i * nn];
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	  c[(i + loop + (j + loop) * ldc) * 2 + 0] +=
    // OpenBLAS❗✔️: 	    subbuffer[(i + j * nn) * 2 + 0] + subbuffer[(j + i * nn) * 2 + 0];
    // OpenBLAS❗✔️: 	  c[(i + loop + (j + loop) * ldc) * 2 + 1] +=
    // OpenBLAS❗✔️: 	    subbuffer[(i + j * nn) * 2 + 1] + subbuffer[(j + i * nn) * 2 + 1];
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifdef LOWER
    // OpenBLAS❗✔️:     GEMM_KERNEL_N(m - mm - nn, nn, k,
    // OpenBLAS❗✔️: 		  alpha_r,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		  alpha_i,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		  a + (mm + nn) * k * COMPSIZE, b + loop * k * COMPSIZE,
    // OpenBLAS❗✔️: 		  c + (mm + nn + loop * ldc) * COMPSIZE, ldc);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   return 0;
    // OpenBLAS❗✔️: }
    if m as isize + offset < 0 {
        return;
    }
    if (n as isize) < offset {
        packed_kernel(m, n, k, alpha, a, b, c, ci, ldc);
        return;
    }
    if offset > 0 {
        let q = offset as usize;
        packed_kernel(m, q, k, alpha, a, b, c, ci, ldc);
        b = &b[q * k..];
        ci += q * ldc;
        n -= q;
        offset = 0;
        if n == 0 {
            return;
        }
    }
    if n as isize > m as isize + offset {
        n = (m as isize + offset) as usize;
        if n == 0 {
            return;
        }
    }
    if offset < 0 {
        let q = (-offset) as usize;
        a = &a[q * k..];
        ci += q;
        m -= q;
        offset = 0;
        if m == 0 {
            return;
        }
    }
    if m > n - offset as usize {
        packed_kernel(m - n, n, k, alpha, &a[n * k..], b, c, ci + n, ldc);
        m = n;
        if m == 0 {
            return;
        }
    }
    let mut scratch = [0_f64; 32 * 32];
    for l in (0..n).step_by(32) {
        let nn = (n - l).min(32);
        if flag {
            scratch[..nn * nn].fill(0.);
            packed_kernel(
                nn,
                nn,
                k,
                alpha,
                &a[l * k..],
                &b[l * k..],
                &mut scratch[..nn * nn],
                0,
                nn,
            );
            for j in 0..nn {
                for i in j..nn {
                    c[ci + i + l + (j + l) * ldc] += scratch[i + j * nn] + scratch[j + i * nn];
                }
            }
        }
        packed_kernel(
            m - l - nn,
            nn,
            k,
            alpha,
            &a[(l + nn) * k..],
            &b[l * k..],
            c,
            ci + l + nn + l * ldc,
            ldc,
        );
    }
}
#[allow(clippy::too_many_arguments)]
pub(super) fn syr2k_lower(
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    beta: f64,
    c: &mut [f64],
    ldc: usize,
) -> Result<(), Blas3Error> {
    // OpenBLAS❗✔️:       ) return 0;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   for(js = n_from; js < n_to; js += GEMM_R){
    // OpenBLAS❗✔️:     min_j = n_to - js;
    // OpenBLAS❗✔️:     if (min_j > GEMM_R) min_j = GEMM_R;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:     m_start = m_from;
    // OpenBLAS❗✔️:     m_end   = js + min_j;
    // OpenBLAS❗✔️:     if (m_end > m_to) m_end = m_to;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     m_start = m_from;
    // OpenBLAS❗✔️:     m_end   = m_to;
    // OpenBLAS❗✔️:     if (m_start < js) m_start = js;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = 0; ls < k; ls += min_l){
    // OpenBLAS❗✔️:       min_l = k - ls;
    // OpenBLAS❗✔️:       if (min_l >= GEMM_Q * 2) {
    // OpenBLAS❗✔️: 	min_l = GEMM_Q;
    // OpenBLAS❗✔️:       } else
    // OpenBLAS❗✔️: 	if (min_l > GEMM_Q) {
    // OpenBLAS❗✔️: 	  min_l = (min_l + 1) / 2;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       min_i = m_end - m_start;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (min_i >= GEMM_P * 2) {
    // OpenBLAS❗✔️: 	min_i = GEMM_P;
    // OpenBLAS❗✔️:       } else
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) {
    // OpenBLAS❗✔️: 	  min_i = ((min_i / 2 + GEMM_UNROLL_MN - 1)/GEMM_UNROLL_MN) * GEMM_UNROLL_MN;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef LOWER
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (m_start >= js) {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_start, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	aa = sb + min_l * (m_start - js)  * COMPSIZE;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	OCOPY_OPERATION(min_l, min_i, b, ldb, ls, m_start, aa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	KERNEL_OPERATION(min_i, min_i, min_l, alpha, sa, aa, c, ldc, m_start, m_start, 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	jjs = m_start + min_i;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       } else {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_start, sa);
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	ICOPY_OPERATION(min_l, min_i, b, ldb, ls, is, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	KERNEL_OPERATION_C(min_i, min_j, min_l, alpha, sa, sb, c, ldc, is, js, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       aa = sb + min_l * (m_start - js) * COMPSIZE;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       ICOPY_OPERATION(min_l, min_i, a, lda, ls, m_start, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       OCOPY_OPERATION(min_l, min_i, b, ldb, ls, m_start, aa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       KERNEL_OPERATION(min_i, MIN(min_i, min_j + js - m_start), min_l, alpha,
    // OpenBLAS❗✔️: 		       sa, aa, c, ldc, m_start, m_start, 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = js; jjs < m_start; jjs += GEMM_UNROLL_MN){
    // OpenBLAS❗✔️: 	min_jj = m_start - jjs;
    // OpenBLAS❗✔️: 	if (min_jj > GEMM_UNROLL_MN) min_jj = GEMM_UNROLL_MN;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	OCOPY_OPERATION(min_l, min_jj, b, ldb, ls, jjs, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	KERNEL_OPERATION(min_i, min_jj, min_l, alpha,
    // OpenBLAS❗✔️: 			 sa, sb + min_l * (jjs - js)  * COMPSIZE, c, ldc, m_start, jjs, 1);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = m_start + min_i; is < m_end; is += min_i){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	min_i = m_end - is;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (min_i >= GEMM_P * 2) {
    // OpenBLAS❗✔️: 	  min_i = GEMM_P;
    // OpenBLAS❗✔️: 	} else
    // OpenBLAS❗✔️: 	  if (min_i > GEMM_P) {
    // OpenBLAS❗✔️: 	    min_i = ((min_i / 2 + GEMM_UNROLL_MN - 1)/GEMM_UNROLL_MN) * GEMM_UNROLL_MN;
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	aa = sb + min_l * (is - js) * COMPSIZE;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (is  < js + min_j) {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ICOPY_OPERATION(min_l, min_i, a, lda, ls, is, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  OCOPY_OPERATION(min_l, min_i, b, ldb, ls, is, aa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  KERNEL_OPERATION(min_i, MIN(min_i, min_j - is + js), min_l, alpha,  sa, aa,  c, ldc, is, is, 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  KERNEL_OPERATION(min_i, is - js, min_l, alpha, sa, sb,  c, ldc, is, js, 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  } else {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	    ICOPY_OPERATION(min_l, min_i, a, lda, ls, is, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	    KERNEL_OPERATION(min_i, min_j, min_l, alpha, sa, sb,  c, ldc, is, js, 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       min_i = m_end - m_start;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (min_i >= GEMM_P * 2) {
    // OpenBLAS❗✔️: 	min_i = GEMM_P;
    // OpenBLAS❗✔️:       } else
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) {
    // OpenBLAS❗✔️: 	  min_i = ((min_i / 2 + GEMM_UNROLL_MN - 1)/GEMM_UNROLL_MN) * GEMM_UNROLL_MN;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       aa = sb + min_l * (m_start - js) * COMPSIZE;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       ICOPY_OPERATION(min_l, min_i, b, ldb, ls, m_start, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       OCOPY_OPERATION(min_l, min_i, a, lda, ls, m_start, aa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       KERNEL_OPERATION_C(min_i, MIN(min_i, min_j + js - m_start), min_l, alpha,
    // OpenBLAS❗✔️: 		       sa, aa, c, ldc, m_start, m_start, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = js; jjs < m_start; jjs += GEMM_UNROLL_MN){
    // OpenBLAS❗✔️: 	min_jj = m_start - jjs;
    // OpenBLAS❗✔️: 	if (min_jj > GEMM_UNROLL_MN) min_jj = GEMM_UNROLL_MN;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	OCOPY_OPERATION(min_l, min_jj, a, lda, ls, jjs, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	KERNEL_OPERATION_C(min_i, min_jj, min_l, alpha,
    // OpenBLAS❗✔️: 			 sa, sb + min_l * (jjs - js)  * COMPSIZE, c, ldc, m_start, jjs, 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = m_start + min_i; is < m_end; is += min_i){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	min_i = m_end - is;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (min_i >= GEMM_P * 2) {
    // OpenBLAS❗✔️: 	  min_i = GEMM_P;
    // OpenBLAS❗✔️: 	} else
    // OpenBLAS❗✔️: 	  if (min_i > GEMM_P) {
    // OpenBLAS❗✔️: 	    min_i = ((min_i / 2 + GEMM_UNROLL_MN - 1)/GEMM_UNROLL_MN) * GEMM_UNROLL_MN;
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	aa = sb + min_l * (is - js) * COMPSIZE;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (is  < js + min_j) {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  ICOPY_OPERATION(min_l, min_i, b, ldb, ls, is, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  OCOPY_OPERATION(min_l, min_i, a, lda, ls, is, aa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  KERNEL_OPERATION_C(min_i, MIN(min_i, min_j - is + js), min_l, alpha,  sa, aa,  c, ldc, is, is, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  KERNEL_OPERATION_C(min_i, is - js, min_l, alpha, sa, sb,  c, ldc, is, js, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  } else {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	    ICOPY_OPERATION(min_l, min_i, b, ldb, ls, is, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	    KERNEL_OPERATION_C(min_i, min_j, min_l, alpha, sa, sb,  c, ldc, is, js, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   return 0;
    // OpenBLAS❗✔️: }
    storage("DSYR2K", "A", n, k, lda, 0, a.len(), 7)?;
    storage("DSYR2K", "B", n, k, ldb, 0, b.len(), 9)?;
    storage("DSYR2K", "C", n, n, ldc, 0, c.len(), 12)?;
    if beta != 1. {
        for j in 0..n {
            for i in j..n {
                let ci = i + j * ldc;
                c[ci] = if beta == 0. { 0. } else { c[ci] * beta };
            }
        }
    }
    if k == 0 || alpha == 0. || n == 0 {
        return Ok(());
    }
    let mut sa = vec![0.; n.min(P) * k.min(Q)];
    let mut sb = vec![0.; k.min(Q) * n.min(n.min(R) + P)];
    for js in (0..n).step_by(R) {
        let nj = (n - js).min(R);
        let mut ls = 0;
        while ls < k {
            let remain = k - ls;
            let kl = if remain >= Q * 2 {
                Q
            } else if remain > Q {
                remain.div_ceil(2)
            } else {
                remain
            };
            for (aa, lda, bb, ldb, flag) in [(a, lda, b, ldb, true), (b, ldb, a, lda, false)] {
                let am = Matrix {
                    separate: Some(aa),
                    offset: 0,
                    ld: lda,
                    trans: Transpose::NoTrans,
                };
                let bm = Matrix {
                    separate: Some(bb),
                    offset: 0,
                    ld: ldb,
                    trans: Transpose::Trans,
                };
                let first = panel_size(n - js, 32);
                pack_a(first, kl, am, c, js, ls, &mut sa[..first * kl]);
                pack_b(kl, first, bm, c, ls, js, &mut sb[..first * kl]);
                syr_kernel(
                    first,
                    first.min(nj),
                    kl,
                    alpha,
                    &sa[..first * kl],
                    &sb[..first * kl],
                    c,
                    js + js * ldc,
                    ldc,
                    0,
                    flag,
                );
                let mut i = js + first;
                while i < n {
                    let h = panel_size(n - i, 32);
                    pack_a(h, kl, am, c, i, ls, &mut sa[..h * kl]);
                    if i < js + nj {
                        let bp = (i - js) * kl;
                        pack_b(kl, h, bm, c, ls, i, &mut sb[bp..bp + h * kl]);
                        syr_kernel(
                            h,
                            h.min(nj + js - i),
                            kl,
                            alpha,
                            &sa[..h * kl],
                            &sb[bp..bp + h * kl],
                            c,
                            i + i * ldc,
                            ldc,
                            0,
                            flag,
                        );
                        syr_kernel(
                            h,
                            i - js,
                            kl,
                            alpha,
                            &sa[..h * kl],
                            &sb,
                            c,
                            i + js * ldc,
                            ldc,
                            (i - js) as isize,
                            flag,
                        );
                    } else {
                        syr_kernel(
                            h,
                            nj,
                            kl,
                            alpha,
                            &sa[..h * kl],
                            &sb,
                            c,
                            i + js * ldc,
                            ldc,
                            (i - js) as isize,
                            flag,
                        );
                    }
                    i += h;
                }
            }
            ls += kl;
        }
    }
    Ok(())
}

#[inline]
fn triangular(
    a: Matrix<'_>,
    shared: &[f64],
    row: usize,
    col: usize,
    upper: bool,
    unit: bool,
) -> f64 {
    if (upper && row > col) || (!upper && row < col) {
        0.
    } else if unit && row == col {
        1.
    } else {
        a.get(shared, row, col)
    }
}
fn pack_tri_a(
    m: usize,
    k: usize,
    a: Matrix<'_>,
    shared: &[f64],
    row: usize,
    col: usize,
    upper: bool,
    unit: bool,
    out: &mut [f64],
) {
    // OpenBLAS❗✔️: #include "common.h"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   BLASLONG i, js, ii;
    // OpenBLAS❗✔️:   BLASLONG X;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   FLOAT *a01, *a02, *a03 ,*a04, *a05, *a06, *a07, *a08;
    // OpenBLAS❗✔️:   FLOAT *a09, *a10, *a11, *a12, *a13, *a14, *a15, *a16;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   js = (n >> 4);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (js > 0){
    // OpenBLAS❗✔️:     do {
    // OpenBLAS❗✔️:       X = posX;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (posX <= posY) {
    // OpenBLAS❗✔️: 	a01 = a + posX + (posY +  0) * lda;
    // OpenBLAS❗✔️: 	a02 = a + posX + (posY +  1) * lda;
    // OpenBLAS❗✔️: 	a03 = a + posX + (posY +  2) * lda;
    // OpenBLAS❗✔️: 	a04 = a + posX + (posY +  3) * lda;
    // OpenBLAS❗✔️: 	a05 = a + posX + (posY +  4) * lda;
    // OpenBLAS❗✔️: 	a06 = a + posX + (posY +  5) * lda;
    // OpenBLAS❗✔️: 	a07 = a + posX + (posY +  6) * lda;
    // OpenBLAS❗✔️: 	a08 = a + posX + (posY +  7) * lda;
    // OpenBLAS❗✔️: 	a09 = a + posX + (posY +  8) * lda;
    // OpenBLAS❗✔️: 	a10 = a + posX + (posY +  9) * lda;
    // OpenBLAS❗✔️: 	a11 = a + posX + (posY + 10) * lda;
    // OpenBLAS❗✔️: 	a12 = a + posX + (posY + 11) * lda;
    // OpenBLAS❗✔️: 	a13 = a + posX + (posY + 12) * lda;
    // OpenBLAS❗✔️: 	a14 = a + posX + (posY + 13) * lda;
    // OpenBLAS❗✔️: 	a15 = a + posX + (posY + 14) * lda;
    // OpenBLAS❗✔️: 	a16 = a + posX + (posY + 15) * lda;
    // OpenBLAS❗✔️:       } else {
    // OpenBLAS❗✔️: 	a01 = a + posY + (posX +  0) * lda;
    // OpenBLAS❗✔️: 	a02 = a + posY + (posX +  1) * lda;
    // OpenBLAS❗✔️: 	a03 = a + posY + (posX +  2) * lda;
    // OpenBLAS❗✔️: 	a04 = a + posY + (posX +  3) * lda;
    // OpenBLAS❗✔️: 	a05 = a + posY + (posX +  4) * lda;
    // OpenBLAS❗✔️: 	a06 = a + posY + (posX +  5) * lda;
    // OpenBLAS❗✔️: 	a07 = a + posY + (posX +  6) * lda;
    // OpenBLAS❗✔️: 	a08 = a + posY + (posX +  7) * lda;
    // OpenBLAS❗✔️: 	a09 = a + posY + (posX +  8) * lda;
    // OpenBLAS❗✔️: 	a10 = a + posY + (posX +  9) * lda;
    // OpenBLAS❗✔️: 	a11 = a + posY + (posX + 10) * lda;
    // OpenBLAS❗✔️: 	a12 = a + posY + (posX + 11) * lda;
    // OpenBLAS❗✔️: 	a13 = a + posY + (posX + 12) * lda;
    // OpenBLAS❗✔️: 	a14 = a + posY + (posX + 13) * lda;
    // OpenBLAS❗✔️: 	a15 = a + posY + (posX + 14) * lda;
    // OpenBLAS❗✔️: 	a16 = a + posY + (posX + 15) * lda;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m >> 4);
    // OpenBLAS❗✔️:       if (i > 0) {
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  if (X < posY) {
    // OpenBLAS❗✔️: 	    a01 += 16;
    // OpenBLAS❗✔️: 	    a02 += 16;
    // OpenBLAS❗✔️: 	    a03 += 16;
    // OpenBLAS❗✔️: 	    a04 += 16;
    // OpenBLAS❗✔️: 	    a05 += 16;
    // OpenBLAS❗✔️: 	    a06 += 16;
    // OpenBLAS❗✔️: 	    a07 += 16;
    // OpenBLAS❗✔️: 	    a08 += 16;
    // OpenBLAS❗✔️: 	    a09 += 16;
    // OpenBLAS❗✔️: 	    a10 += 16;
    let (mut i, mut p) = (0, 0);
    while i < m {
        let h = row_tile(m - i);
        for l in 0..k {
            for ii in 0..h {
                out[p] = triangular(a, shared, row + i + ii, col + l, upper, unit);
                p += 1;
            }
        }
        i += h;
    }
}
fn pack_tri_b(
    k: usize,
    n: usize,
    b: Matrix<'_>,
    shared: &[f64],
    row: usize,
    col: usize,
    upper: bool,
    unit: bool,
    out: &mut [f64],
) {
    // OpenBLAS❗✔️: #include "common.h"
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   BLASLONG i, js;
    // OpenBLAS❗✔️:   BLASLONG X;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   FLOAT data01, data02, data03, data04;
    // OpenBLAS❗✔️:   FLOAT *ao1, *ao2;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   js = (n >> 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (js > 0){
    // OpenBLAS❗✔️:     do {
    // OpenBLAS❗✔️:       X = posX;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (posX <= posY) {
    // OpenBLAS❗✔️: 	ao1 = a + posX + (posY + 0) * lda;
    // OpenBLAS❗✔️: 	ao2 = a + posX + (posY + 1) * lda;
    // OpenBLAS❗✔️:       } else {
    // OpenBLAS❗✔️: 	ao1 = a + posY + (posX + 0) * lda;
    // OpenBLAS❗✔️: 	ao2 = a + posY + (posX + 1) * lda;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       i = (m >> 1);
    // OpenBLAS❗✔️:       if (i > 0) {
    // OpenBLAS❗✔️: 	do {
    // OpenBLAS❗✔️: 	  if (X < posY) {
    // OpenBLAS❗✔️: 	    data01 = *(ao1 + 0);
    // OpenBLAS❗✔️: 	    data02 = *(ao1 + 1);
    // OpenBLAS❗✔️: 	    data03 = *(ao2 + 0);
    // OpenBLAS❗✔️: 	    data04 = *(ao2 + 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	    b[ 0] = data01;
    // OpenBLAS❗✔️: 	    b[ 1] = data03;
    // OpenBLAS❗✔️: 	    b[ 2] = data02;
    // OpenBLAS❗✔️: 	    b[ 3] = data04;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	    ao1 += 2;
    // OpenBLAS❗✔️: 	    ao2 += 2;
    // OpenBLAS❗✔️: 	    b   += 4;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  } else
    // OpenBLAS❗✔️: 	    if (X > posY) {
    // OpenBLAS❗✔️: 	      ao1 += 2 * lda;
    // OpenBLAS❗✔️: 	      ao2 += 2 * lda;
    // OpenBLAS❗✔️: 	      b += 4;
    // OpenBLAS❗✔️: 	    } else {
    // OpenBLAS❗✔️: #ifdef UNIT
    // OpenBLAS❗✔️: 	      data03 = *(ao2 + 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	      b[ 0] = ONE;
    // OpenBLAS❗✔️: 	      b[ 1] = data03;
    // OpenBLAS❗✔️: 	      b[ 2] = ZERO;
    // OpenBLAS❗✔️: 	      b[ 3] = ONE;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	      data01 = *(ao1 + 0);
    // OpenBLAS❗✔️: 	      data03 = *(ao2 + 0);
    // OpenBLAS❗✔️: 	      data04 = *(ao2 + 1);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	      b[ 0] = data01;
    // OpenBLAS❗✔️: 	      b[ 1] = data03;
    // OpenBLAS❗✔️: 	      b[ 2] = ZERO;
    // OpenBLAS❗✔️: 	      b[ 3] = data04;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	      ao1 += 2 * lda;
    // OpenBLAS❗✔️: 	      ao2 += 2 * lda;
    // OpenBLAS❗✔️: 	      b += 4;
    // OpenBLAS❗✔️: 	    }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  X += 2;
    // OpenBLAS❗✔️: 	  i --;
    // OpenBLAS❗✔️: 	} while (i > 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       if (m & 1) {
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (X < posY) {
    // OpenBLAS❗✔️: 	  data01 = *(ao1 + 0);
    // OpenBLAS❗✔️: 	  data03 = *(ao2 + 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  b[ 0] = data01;
    // OpenBLAS❗✔️: 	  b[ 1] = data03;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	  /* ao1 += 1;
    let (mut j, mut p) = (0, 0);
    while j < n {
        let w = (n - j).min(2);
        for l in 0..k {
            for jj in 0..w {
                out[p] = triangular(b, shared, row + l, col + j + jj, upper, unit);
                p += 1;
            }
        }
        j += w;
    }
}
fn trmm_right(
    upper: bool,
    unit: bool,
    m: usize,
    n: usize,
    a: Matrix<'_>,
    data: &mut [f64],
    bo: usize,
    ldb: usize,
) {
    // OpenBLAS❗✔️:     if ((beta[0] != ONE) || (beta[1] != ZERO))
    // OpenBLAS❗✔️:       GEMM_BETA(m, n, 0, beta[0], beta[1], NULL, 0, NULL, 0, b, ldb);
    // OpenBLAS❗✔️:     if ((beta[0] == ZERO) && (beta[1] == ZERO)) return 0;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if (!defined(UPPER) && !defined(TRANSA)) || (defined(UPPER) && defined(TRANSA))
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   for(js = 0; js < n; js += GEMM_R){
    // OpenBLAS❗✔️:     min_j = n - js;
    // OpenBLAS❗✔️:     if (min_j > GEMM_R) min_j = GEMM_R;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = js; ls < js + min_j; ls += GEMM_Q){
    // OpenBLAS❗✔️:       min_l = js + min_j - ls;
    // OpenBLAS❗✔️:       if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:       min_i = m;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_ITCOPY(min_l, min_i, b + (ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = 0; jjs < ls - js; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = ls - js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	GEMM_ONCOPY(min_l, min_jj, a + (ls + (js + jjs) * lda) * COMPSIZE, lda, sb + min_l * jjs * COMPSIZE);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	GEMM_OTCOPY(min_l, min_jj, a + ((js + jjs) + ls * lda) * COMPSIZE, lda, sb + min_l * jjs * COMPSIZE);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb + min_l * jjs * COMPSIZE,
    // OpenBLAS❗✔️: 		    b + ((js + jjs) * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = 0; jjs < min_l; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = min_l - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	TRMM_OLNCOPY(min_l, min_jj, a, lda, ls, ls + jjs, sb + min_l * (ls - js + jjs) * COMPSIZE);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	TRMM_OUTCOPY(min_l, min_jj, a, lda, ls, ls + jjs, sb + min_l * (ls - js + jjs) * COMPSIZE);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_T(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa,
    // OpenBLAS❗✔️: 		      sb + (ls - js + jjs) * min_l * COMPSIZE,
    // OpenBLAS❗✔️: 		      b + ((ls + jjs) * ldb) * COMPSIZE, ldb, -jjs);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = min_i; is < m; is += GEMM_P){
    // OpenBLAS❗✔️: 	min_i = m - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_ITCOPY(min_l, min_i, b + (is + ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, ls - js, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb,
    // OpenBLAS❗✔️: 		    b + (is + js * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_T(min_i, min_l, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa,
    // OpenBLAS❗✔️: 		      sb + (ls - js) * min_l * COMPSIZE,
    // OpenBLAS❗✔️: 		      b + (is + ls * ldb) * COMPSIZE, ldb, 0);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = js + min_j; ls < n; ls += GEMM_Q){
    // OpenBLAS❗✔️:       min_l = n - ls;
    // OpenBLAS❗✔️:       if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:       min_i = m;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_ITCOPY(min_l, min_i, b + (ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	GEMM_ONCOPY(min_l, min_jj, a + (ls + jjs * lda) * COMPSIZE, lda, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	GEMM_OTCOPY(min_l, min_jj, a + (jjs + ls * lda) * COMPSIZE, lda, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb + min_l * (jjs - js) * COMPSIZE,
    // OpenBLAS❗✔️: 		    b + (jjs * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = min_i; is < m; is += GEMM_P){
    // OpenBLAS❗✔️: 	min_i = m - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_ITCOPY(min_l, min_i, b + (is + ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb, b + (is + js * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:   for(js = n; js > 0; js -= GEMM_R){
    // OpenBLAS❗✔️:     min_j = js;
    // OpenBLAS❗✔️:     if (min_j > GEMM_R) min_j = GEMM_R;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     start_ls = js - min_j;
    // OpenBLAS❗✔️:     while (start_ls + GEMM_Q < js) start_ls += GEMM_Q;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = start_ls; ls >= js - min_j; ls -= GEMM_Q){
    // OpenBLAS❗✔️:       min_l = js - ls;
    // OpenBLAS❗✔️:       if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:       min_i = m;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_ITCOPY(min_l, min_i, b + (ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = 0; jjs < min_l; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = min_l - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	TRMM_OUNCOPY(min_l, min_jj, a, lda, ls, ls + jjs, sb + min_l * jjs * COMPSIZE);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	TRMM_OLTCOPY(min_l, min_jj, a, lda, ls, ls + jjs, sb + min_l * jjs * COMPSIZE);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_N(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa,
    // OpenBLAS❗✔️: 		      sb + min_l * jjs * COMPSIZE,
    // OpenBLAS❗✔️: 		      b + ((ls + jjs) * ldb) * COMPSIZE, ldb, -jjs);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = 0; jjs < js - ls - min_l; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = js - ls - min_l - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	GEMM_ONCOPY(min_l, min_jj, a + (ls + (ls + min_l + jjs) * lda) * COMPSIZE, lda,
    // OpenBLAS❗✔️: 		    sb + min_l * (min_l + jjs) * COMPSIZE);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	GEMM_OTCOPY(min_l, min_jj, a + ((ls + min_l + jjs) + ls * lda) * COMPSIZE, lda,
    // OpenBLAS❗✔️: 		    sb + min_l * (min_l + jjs) * COMPSIZE);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa,
    // OpenBLAS❗✔️: 		    sb + min_l * (min_l + jjs) * COMPSIZE,
    // OpenBLAS❗✔️: 		    b + ((ls + min_l + jjs) * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = min_i; is < m; is += GEMM_P){
    // OpenBLAS❗✔️: 	min_i = m - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_ITCOPY(min_l, min_i, b + (is + ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_N(min_i, min_l, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa,
    // OpenBLAS❗✔️: 		      sb,
    // OpenBLAS❗✔️: 		      b + (is + ls * ldb) * COMPSIZE, ldb, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if (js - ls - min_l > 0) {
    // OpenBLAS❗✔️: 	  GEMM_KERNEL(min_i, js - ls - min_l, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa,
    // OpenBLAS❗✔️: 		      sb + min_l * min_l * COMPSIZE,
    // OpenBLAS❗✔️: 		      b + (is + (ls + min_l) * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = 0; ls < js - min_j; ls += GEMM_Q){
    // OpenBLAS❗✔️:       min_l = js - min_j - ls;
    // OpenBLAS❗✔️:       if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:       min_i = m;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_ITCOPY(min_l, min_i, b + (ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	GEMM_ONCOPY(min_l, min_jj, a + (ls + (jjs - min_j) * lda) * COMPSIZE, lda, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	GEMM_OTCOPY(min_l, min_jj, a + ((jjs - min_j) + ls * lda) * COMPSIZE, lda, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb + min_l * (jjs - js) * COMPSIZE,
    // OpenBLAS❗✔️: 		    b + ((jjs - min_j) * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = min_i; is < m; is += GEMM_P){
    // OpenBLAS❗✔️: 	min_i = m - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_ITCOPY(min_l, min_i, b + (is + ls * ldb) * COMPSIZE, ldb, sa);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb, b + (is + (js - min_j) * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   return 0;
    // OpenBLAS❗✔️: }
    let bm = Matrix {
        separate: None,
        offset: bo,
        ld: ldb,
        trans: Transpose::NoTrans,
    };
    let mut sa = vec![0.; P * Q];
    let mut sb = vec![0.; Q * n.min(R)];
    let mut edge = if upper { n } else { 0 };
    while if upper { edge > 0 } else { edge < n } {
        let nj = if upper {
            edge.min(R)
        } else {
            (n - edge).min(R)
        };
        let js = if upper { edge - nj } else { edge };
        let je = js + nj;
        let mut cursor = if upper { js + (nj - 1) / Q * Q } else { js };
        loop {
            let ls = cursor;
            let kl = (je - ls).min(Q);
            let first = m.min(P);
            pack_a(first, kl, bm, data, 0, ls, &mut sa[..first * kl]);
            let (offstart, offend) = if upper { (ls + kl, je) } else { (js, ls) };
            // With a lower multiplier, earlier columns are updated before
            // overwriting this diagonal block; with upper, source reverses it.
            if !upper {
                for (j, w) in driver_chunks(offstart, offend) {
                    let bp = (j - js) * kl;
                    pack_b(kl, w, a, data, ls, j, &mut sb[bp..bp + kl * w]);
                    packed_kernel(
                        first,
                        w,
                        kl,
                        1.,
                        &sa[..first * kl],
                        &sb[bp..bp + kl * w],
                        data,
                        bo + j * ldb,
                        ldb,
                    );
                }
            }
            for (j, w) in driver_chunks(ls, ls + kl) {
                let bp = (j - js) * kl;
                pack_tri_b(kl, w, a, data, ls, j, upper, unit, &mut sb[bp..bp + kl * w]);
                packed_product(
                    first,
                    w,
                    kl,
                    1.,
                    &sa[..first * kl],
                    &sb[bp..bp + kl * w],
                    data,
                    bo + j * ldb,
                    ldb,
                    Some(TriKernel {
                        side: Side::Right,
                        upper,
                        diagonal: (j - ls) as isize,
                    }),
                );
            }
            if upper {
                for (j, w) in driver_chunks(offstart, offend) {
                    let bp = (j - js) * kl;
                    pack_b(kl, w, a, data, ls, j, &mut sb[bp..bp + kl * w]);
                    packed_kernel(
                        first,
                        w,
                        kl,
                        1.,
                        &sa[..first * kl],
                        &sb[bp..bp + kl * w],
                        data,
                        bo + j * ldb,
                        ldb,
                    );
                }
            }
            for i in (first..m).step_by(P) {
                let h = (m - i).min(P);
                pack_a(h, kl, bm, data, i, ls, &mut sa[..h * kl]);
                if !upper && offstart < offend {
                    packed_kernel(
                        h,
                        offend - offstart,
                        kl,
                        1.,
                        &sa[..h * kl],
                        &sb[..(offend - offstart) * kl],
                        data,
                        bo + i + js * ldb,
                        ldb,
                    );
                }
                let bp = (ls - js) * kl;
                packed_product(
                    h,
                    kl,
                    kl,
                    1.,
                    &sa[..h * kl],
                    &sb[bp..bp + kl * kl],
                    data,
                    bo + i + ls * ldb,
                    ldb,
                    Some(TriKernel {
                        side: Side::Right,
                        upper,
                        diagonal: 0,
                    }),
                );
                if upper && offstart < offend {
                    let bp = (offstart - js) * kl;
                    packed_kernel(
                        h,
                        offend - offstart,
                        kl,
                        1.,
                        &sa[..h * kl],
                        &sb[bp..bp + (offend - offstart) * kl],
                        data,
                        bo + i + offstart * ldb,
                        ldb,
                    );
                }
            }
            if upper {
                if cursor == js {
                    break;
                }
                cursor -= Q;
            } else {
                cursor += kl;
                if cursor >= je {
                    break;
                }
            }
        }
        let (external_start, external_end) = if upper { (0, js) } else { (je, n) };
        for ls in (external_start..external_end).step_by(Q) {
            let kl = (external_end - ls).min(Q);
            let first = m.min(P);
            pack_a(first, kl, bm, data, 0, ls, &mut sa[..first * kl]);
            for (j, w) in driver_chunks(js, je) {
                let bp = (j - js) * kl;
                pack_b(kl, w, a, data, ls, j, &mut sb[bp..bp + kl * w]);
                packed_kernel(
                    first,
                    w,
                    kl,
                    1.,
                    &sa[..first * kl],
                    &sb[bp..bp + kl * w],
                    data,
                    bo + j * ldb,
                    ldb,
                );
            }
            for i in (first..m).step_by(P) {
                let h = (m - i).min(P);
                pack_a(h, kl, bm, data, i, ls, &mut sa[..h * kl]);
                packed_kernel(
                    h,
                    nj,
                    kl,
                    1.,
                    &sa[..h * kl],
                    &sb[..nj * kl],
                    data,
                    bo + i + js * ldb,
                    ldb,
                );
            }
        }
        edge = if upper { js } else { je };
    }
}
fn left_rows(m: usize) -> usize {
    let h = m.min(P);
    if h > 16 { h / 16 * 16 } else { h }
}
fn trmm_left(
    upper: bool,
    unit: bool,
    m: usize,
    n: usize,
    a: Matrix<'_>,
    data: &mut [f64],
    bo: usize,
    ldb: usize,
) {
    // OpenBLAS❗✔️:     if ((beta[0] != ONE) || (beta[1] != ZERO))
    // OpenBLAS❗✔️:       GEMM_BETA(m, n, 0, beta[0], beta[1], NULL, 0, NULL, 0, b, ldb);
    // OpenBLAS❗✔️:     if ((beta[0] == ZERO) && (beta[1] == ZERO)) return 0;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   for(js = 0; js < n; js += GEMM_R){
    // OpenBLAS❗✔️:     min_j = n - js;
    // OpenBLAS❗✔️:     if (min_j > GEMM_R) min_j = GEMM_R;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if (defined(UPPER) && !defined(TRANSA)) || (!defined(UPPER) && defined(TRANSA))
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     min_l = m;
    // OpenBLAS❗✔️:     if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:     min_i = min_l;
    // OpenBLAS❗✔️:     if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:     if( min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️:         min_i =  (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️:     TRMM_IUTCOPY(min_l, min_i, a, lda, 0, 0, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     TRMM_ILNCOPY(min_l, min_i, a, lda, 0, 0, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️:       min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️:       /* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️:       if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️:       else
    // OpenBLAS❗✔️:       	if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_ONCOPY(min_l, min_jj, b + (jjs * ldb) * COMPSIZE, ldb, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(outercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       TRMM_KERNEL_N(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb + min_l * (jjs - js) * COMPSIZE, b + (jjs * ldb) * COMPSIZE, ldb, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(is = min_i; is < min_l; is += min_i){
    // OpenBLAS❗✔️:       min_i = min_l - is;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:       if( min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️:           min_i =  (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:         }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️:       TRMM_IUTCOPY(min_l, min_i, a, lda, 0, is, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       TRMM_ILNCOPY(min_l, min_i, a, lda, 0, is, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       TRMM_KERNEL_N(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb, b + (is + js * ldb) * COMPSIZE, ldb, is);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = min_l; ls < m; ls += GEMM_Q){
    // OpenBLAS❗✔️:       min_l = m - ls;
    // OpenBLAS❗✔️:       if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:       min_i = ls;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:       if( min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️:           min_i =  (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:           }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️:       GEMM_ITCOPY(min_l, min_i, a + (ls * lda) * COMPSIZE, lda, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       GEMM_INCOPY(min_l, min_i, a + (ls      ) * COMPSIZE, lda, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:         if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️:         else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_ONCOPY(min_l, min_jj, b + (ls + jjs * ldb) * COMPSIZE, ldb, sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(gemmcost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb + min_l * (jjs - js) * COMPSIZE,
    // OpenBLAS❗✔️: 		    b + (jjs * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(gemmcost);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = min_i; is < ls; is += min_i){
    // OpenBLAS❗✔️: 	min_i = ls - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️: 	if( min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️: 	    min_i =  (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:           }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	GEMM_ITCOPY(min_l, min_i, a + (is + ls * lda) * COMPSIZE, lda, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	GEMM_INCOPY(min_l, min_i, a + (ls + is * lda) * COMPSIZE, lda, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb, b + (is + js * ldb) * COMPSIZE, ldb);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(gemmcost);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = ls; is < ls + min_l; is += min_i){
    // OpenBLAS❗✔️: 	min_i = ls + min_l - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️: 	if( min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️: 	    min_i =  (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	TRMM_IUTCOPY(min_l, min_i, a, lda, ls, is, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	TRMM_ILNCOPY(min_l, min_i, a, lda, ls, is, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_N(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa, sb, b + (is + js * ldb) * COMPSIZE, ldb, is - ls);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     min_l = m;
    // OpenBLAS❗✔️:     if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:     min_i = min_l;
    // OpenBLAS❗✔️:     if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:     if (min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️:         min_i = (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️:     TRMM_ILTCOPY(min_l, min_i, a, lda, m - min_l, m - min_l, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     TRMM_IUNCOPY(min_l, min_i, a, lda, m - min_l, m - min_l, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️:       min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️:       /* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️:       if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️:       else
    // OpenBLAS❗✔️:         if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       GEMM_ONCOPY(min_l, min_jj, b + (m - min_l + jjs * ldb) * COMPSIZE, ldb,
    // OpenBLAS❗✔️: 		  sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(outercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       TRMM_KERNEL_T(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb + min_l * (jjs - js) * COMPSIZE,
    // OpenBLAS❗✔️: 		    b + (m - min_l + jjs * ldb) * COMPSIZE, ldb, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(is = m - min_l + min_i; is < m; is += min_i){
    // OpenBLAS❗✔️:       min_i = m - is;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:       if (min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️:           min_i = (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️:       TRMM_ILTCOPY(min_l, min_i, a, lda, m - min_l, is, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       TRMM_IUNCOPY(min_l, min_i, a, lda, m - min_l, is, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       TRMM_KERNEL_T(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		    ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		    sa, sb, b + (is + js * ldb) * COMPSIZE, ldb, is - m + min_l);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:     }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     for(ls = m - min_l; ls > 0; ls -= GEMM_Q){
    // OpenBLAS❗✔️:       min_l = ls;
    // OpenBLAS❗✔️:       if (min_l > GEMM_Q) min_l = GEMM_Q;
    // OpenBLAS❗✔️:       min_i = min_l;
    // OpenBLAS❗✔️:       if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️:       if (min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️:           min_i = (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️:       TRMM_ILTCOPY(min_l, min_i, a, lda, ls - min_l, ls - min_l, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:       TRMM_IUNCOPY(min_l, min_i, a, lda, ls - min_l, ls - min_l, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(jjs = js; jjs < js + min_j; jjs += min_jj){
    // OpenBLAS❗✔️: 	min_jj = min_j + js - jjs;
    // OpenBLAS❗✔️: #if defined(SKYLAKEX) || defined(COOPERLAKE) || defined(SAPPHIRERAPIDS)
    // OpenBLAS❗✔️: 	/* the current AVX512 s/d/c/z GEMM kernel requires n>=6*GEMM_UNROLL_N to achieve the best performance */
    // OpenBLAS❗✔️: 	if (min_jj >= 6*GEMM_UNROLL_N) min_jj = 6*GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:         if (min_jj >= GEMM_UNROLL_N*3) min_jj = GEMM_UNROLL_N*3;
    // OpenBLAS❗✔️:         else
    // OpenBLAS❗✔️: 	  if (min_jj > GEMM_UNROLL_N) min_jj = GEMM_UNROLL_N;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_ONCOPY(min_l, min_jj, b + (ls - min_l + jjs * ldb) * COMPSIZE, ldb,
    // OpenBLAS❗✔️: 		    sb + min_l * (jjs - js) * COMPSIZE);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(outercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_T(min_i, min_jj, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa, sb + min_l * (jjs - js) * COMPSIZE,
    // OpenBLAS❗✔️: 		      b + (ls - min_l + jjs * ldb) * COMPSIZE, ldb, 0);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = ls - min_l + min_i; is < ls; is += min_i){
    // OpenBLAS❗✔️: 	min_i = ls - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️: 	    min_i = (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	TRMM_ILTCOPY(min_l, min_i, a, lda, ls - min_l, is, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	TRMM_IUNCOPY(min_l, min_i, a, lda, ls - min_l, is, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	TRMM_KERNEL_T(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    // OpenBLAS❗✔️: 		      ZERO,
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️: 		      sa, sb, b + (is + js * ldb) * COMPSIZE, ldb, is - ls + min_l);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(trmmcost);
    // OpenBLAS❗✔️:       }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:       for(is = ls; is < m; is += min_i){
    // OpenBLAS❗✔️: 	min_i = m - is;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_P) min_i = GEMM_P;
    // OpenBLAS❗✔️: 	if (min_i > GEMM_UNROLL_M){
    // OpenBLAS❗✔️: 	    min_i = (min_i / GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #ifndef TRANSA
    // OpenBLAS❗✔️: 	GEMM_ITCOPY(min_l, min_i, a + (is + (ls - min_l) * lda) * COMPSIZE, lda, sa);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	GEMM_INCOPY(min_l, min_i, a + ((ls - min_l) + is * lda) * COMPSIZE, lda, sa);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	STOP_RPCC(innercost);
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	START_RPCC();
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	GEMM_KERNEL(min_i, min_j, min_l, dp1,
    // OpenBLAS❗✔️: #ifdef COMPLEX
    let bm = Matrix {
        separate: None,
        offset: bo,
        ld: ldb,
        trans: Transpose::NoTrans,
    };
    let mut sa = vec![0.; P * Q];
    let mut sb = vec![0.; Q * n.min(R)];
    for js in (0..n).step_by(R) {
        let nj = (n - js).min(R);
        let mut edge = if upper { 0 } else { m };
        while if upper { edge < m } else { edge > 0 } {
            let kl = if upper {
                (m - edge).min(Q)
            } else {
                edge.min(Q)
            };
            let ls = if upper { edge } else { edge - kl };
            let le = ls + kl;
            let (offstart, offend) = if upper { (0, ls) } else { (le, m) };
            // Source upper path packs the whole B panel while updating the
            // first off-diagonal A panel; lower packs it with the first
            // triangular A panel. The data are retained for remaining rows.
            if upper && offend > 0 {
                let first = left_rows(offend);
                pack_a(first, kl, a, data, 0, ls, &mut sa[..first * kl]);
                for (j, w) in driver_chunks(js, js + nj) {
                    let bp = (j - js) * kl;
                    pack_b(kl, w, bm, data, ls, j, &mut sb[bp..bp + kl * w]);
                    packed_kernel(
                        first,
                        w,
                        kl,
                        1.,
                        &sa[..first * kl],
                        &sb[bp..bp + kl * w],
                        data,
                        bo + j * ldb,
                        ldb,
                    );
                }
                let mut i = first;
                while i < offend {
                    let h = left_rows(offend - i);
                    pack_a(h, kl, a, data, i, ls, &mut sa[..h * kl]);
                    packed_kernel(
                        h,
                        nj,
                        kl,
                        1.,
                        &sa[..h * kl],
                        &sb[..nj * kl],
                        data,
                        bo + i + js * ldb,
                        ldb,
                    );
                    i += h;
                }
            }
            let first = left_rows(kl);
            pack_tri_a(
                first,
                kl,
                a,
                data,
                ls,
                ls,
                upper,
                unit,
                &mut sa[..first * kl],
            );
            if !upper || ls == 0 {
                for (j, w) in driver_chunks(js, js + nj) {
                    let bp = (j - js) * kl;
                    pack_b(kl, w, bm, data, ls, j, &mut sb[bp..bp + kl * w]);
                    packed_product(
                        first,
                        w,
                        kl,
                        1.,
                        &sa[..first * kl],
                        &sb[bp..bp + kl * w],
                        data,
                        bo + ls + j * ldb,
                        ldb,
                        Some(TriKernel {
                            side: Side::Left,
                            upper,
                            diagonal: 0,
                        }),
                    );
                }
            } else {
                packed_product(
                    first,
                    nj,
                    kl,
                    1.,
                    &sa[..first * kl],
                    &sb[..nj * kl],
                    data,
                    bo + ls + js * ldb,
                    ldb,
                    Some(TriKernel {
                        side: Side::Left,
                        upper,
                        diagonal: 0,
                    }),
                );
            }
            let mut i = ls + first;
            while i < le {
                let h = left_rows(le - i);
                pack_tri_a(h, kl, a, data, i, ls, upper, unit, &mut sa[..h * kl]);
                packed_product(
                    h,
                    nj,
                    kl,
                    1.,
                    &sa[..h * kl],
                    &sb[..nj * kl],
                    data,
                    bo + i + js * ldb,
                    ldb,
                    Some(TriKernel {
                        side: Side::Left,
                        upper,
                        diagonal: (i - ls) as isize,
                    }),
                );
                i += h;
            }
            if !upper {
                let mut i = offstart;
                while i < offend {
                    let h = left_rows(offend - i);
                    pack_a(h, kl, a, data, i, ls, &mut sa[..h * kl]);
                    packed_kernel(
                        h,
                        nj,
                        kl,
                        1.,
                        &sa[..h * kl],
                        &sb[..nj * kl],
                        data,
                        bo + i + js * ldb,
                        ldb,
                    );
                    i += h;
                }
            }
            edge = if upper { le } else { ls };
        }
    }
}
fn trmm_inner(
    side: Side,
    triangle: Triangle,
    trans: Transpose,
    diag: Diagonal,
    m: usize,
    n: usize,
    alpha: f64,
    a: Matrix<'_>,
    b: &mut [f64],
    bo: usize,
    ldb: usize,
) -> Result<(), Blas3Error> {
    // OpenBLAS❗✔️:     m = m_to - m_from;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:     b += m_from * COMPSIZE;
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   if (beta) {
    // OpenBLAS❗✔️: #ifndef COMPLEX
    // OpenBLAS❗✔️:     if (beta[0] != ONE)
    // OpenBLAS❗✔️:       GEMM_BETA(m, n, 0, beta[0], NULL, 0, NULL, 0, b, ldb);
    // OpenBLAS❗✔️:     if (beta[0] == ZERO) return 0;
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️:     if ((beta[0] != ONE) || (beta[1] != ZERO))
    // OpenBLAS❗✔️:       GEMM_BETA(m, n, 0, beta[0], beta[1], NULL, 0, NULL, 0, b, ldb);
    // OpenBLAS❗✔️:     if ((beta[0] == ZERO) && (beta[1] == ZERO)) return 0;
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:   }
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if (!defined(UPPER) && !defined(TRANSA)) || (defined(UPPER) && defined(TRANSA))
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:   for(js = 0; js < n; js += GEMM_R){
    // OpenBLAS❗✔️:     min_j = n - js;
    let dim = if side == Side::Left { m } else { n };
    op_storage("DTRMM", "A", dim, dim, a, b.len(), 9)?;
    storage("DTRMM", "B", m, n, ldb, bo, b.len(), 11)?;
    independent(dim, dim, a, m, n, bo, ldb, "DTRMM", 8)?;
    if m == 0 || n == 0 {
        return Ok(());
    }
    scale(m, n, alpha, b, bo, ldb);
    if alpha == 0. {
        return Ok(());
    }
    let upper = (triangle == Triangle::Upper) != (trans == Transpose::Trans);
    if side == Side::Right {
        trmm_right(upper, diag == Diagonal::Unit, m, n, a, b, bo, ldb);
    } else {
        trmm_left(upper, diag == Diagonal::Unit, m, n, a, b, bo, ldb);
    }
    Ok(())
}
#[allow(clippy::too_many_arguments)]
pub(super) fn trmm(
    side: Side,
    triangle: Triangle,
    trans: Transpose,
    diag: Diagonal,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
) -> Result<(), Blas3Error> {
    trmm_inner(
        side,
        triangle,
        trans,
        diag,
        m,
        n,
        alpha,
        Matrix {
            separate: Some(a),
            offset: 0,
            ld: lda,
            trans,
        },
        b,
        0,
        ldb,
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    fn data(n: usize, seed: usize) -> Vec<f64> {
        (0..n)
            .map(|i| (((i * 37 + seed * 13) % 127) as f64 - 63.) / 17.)
            .collect()
    }
    fn hash(v: &[f64]) -> u64 {
        v.iter().fold(0xcbf29ce484222325_u64, |h, x| {
            (h ^ x.to_bits()).wrapping_mul(0x100000001b3)
        })
    }
    #[test]
    fn openblas3_packing_block_tails_and_transposes() {
        // Actual OpenBLAS 446c single-thread SkylakeX output bits, including
        // untouched column padding. One owner regression exercises dispatch,
        // packing, M/N/K blocks, transpose and triangular kernel boundaries.
        use Transpose::{NoTrans as N, Trans as T};
        let mut failures = Vec::new();
        let gemm_cases: &[(Transpose, Transpose, usize, usize, usize, u64)] = &[
            (N, N, 3, 7, 35, 6583727534588862214),
            (N, N, 7, 5, 13, 9341543982069039717),
            (T, N, 7, 5, 35, 18301563737336003168),
            (N, T, 7, 5, 13, 14707173111442767769),
            (T, T, 7, 5, 13, 17601126659619605538),
            (N, N, 23, 59, 769, 15115784566659338684),
            (N, T, 23, 59, 769, 12960827240900065200),
            (T, N, 23, 59, 769, 12683977528994847788),
            (T, T, 23, 59, 769, 6134089932201126934),
            (N, N, 389, 5, 769, 13397492776733283728),
            (N, N, 2, 10705, 65, 14112443786020398363),
        ];
        for &(ta, tb, m, n, k, expected) in gemm_cases {
            let lda = if ta == N { m + 1 } else { k + 1 };
            let ldb = if tb == N { k + 1 } else { n + 1 };
            let ldc = m + 2;
            let a = data(lda * if ta == N { k } else { m }, 1);
            let b = data(ldb * if tb == N { n } else { k }, 2);
            let mut c = data(ldc * n, 3);
            gemm(ta, tb, m, n, k, 0.7, &a, lda, &b, ldb, -0.3, &mut c, ldc).unwrap();
            let actual = hash(&c);
            if actual != expected {
                failures.push(format!(
                    "GEMM {ta:?}/{tb:?} {m}/{n}/{k}: {actual} != {expected}"
                ));
            }
        }
        // x86_64 does not enable GEMM_GEMV_FORWARD. These axes still use
        // small/packed GEMM, including the NN K=15/K=16 tail crossover.
        let b0_cases: &[(Transpose, Transpose, usize, usize, usize, u64)] = &[
            (N, N, 1, 10, 10, 14230941537548679026),
            (N, T, 1, 10, 10, 16072053539473217814),
            (T, N, 1, 10, 10, 13854560869759696902),
            (T, T, 1, 10, 10, 15998136326706971664),
            (N, N, 11, 1, 35, 4875885622875864571),
            (N, T, 11, 1, 35, 9462155637298142013),
            (T, N, 11, 1, 35, 16167909369190688562),
            (T, T, 11, 1, 35, 5376703158766575206),
            (N, N, 11, 10, 15, 11042005802762824451),
            (N, N, 11, 10, 16, 3491336554642472481),
        ];
        for &(ta, tb, m, n, k, expected) in b0_cases {
            let lda = if ta == N { m + 1 } else { k + 1 };
            let ldb = if tb == N { k + 1 } else { n + 1 };
            let ldc = m + 2;
            let a = data(lda * if ta == N { k } else { m }, 1);
            let b = data(ldb * if tb == N { n } else { k }, 2);
            let mut c = data(ldc * n, 3);
            gemm(ta, tb, m, n, k, 1., &a, lda, &b, ldb, 0., &mut c, ldc).unwrap();
            let actual = hash(&c);
            if actual != expected {
                failures.push(format!(
                    "B0 GEMM {ta:?}/{tb:?} {m}/{n}/{k}: {actual} != {expected}"
                ));
            }
        }
        let syr_cases: &[(usize, usize, u64)] = &[
            (35, 13, 13143543308237184550),
            (199, 385, 18318609062981588617),
        ];
        for &(n, k, expected) in syr_cases {
            let a = data((n + 1) * k, 1);
            let b = data((n + 2) * k, 2);
            let mut c = data((n + 3) * n, 3);
            syr2k_lower(n, k, 0.7, &a, n + 1, &b, n + 2, -0.3, &mut c, n + 3).unwrap();
            let actual = hash(&c);
            if actual != expected {
                failures.push(format!("SYR2K {n}/{k}: {actual} != {expected}"));
            }
        }
        use Diagonal::{NonUnit as Full, Unit as One};
        use Side::{Left as L, Right as R};
        use Triangle::{Lower as D, Upper as U};
        let tri_cases: &[(Side, Triangle, Transpose, Diagonal, usize, usize, u64)] = &[
            (L, U, N, One, 23, 13, 8833656859852383639),
            (L, U, N, Full, 23, 13, 3244722745600105478),
            (L, U, T, One, 23, 13, 14558589841839089153),
            (L, U, T, Full, 23, 13, 14686811153650557173),
            (L, D, N, One, 23, 13, 8238931179595325329),
            (L, D, N, Full, 23, 13, 17251781386062409302),
            (L, D, T, One, 23, 13, 476847641466253459),
            (L, D, T, Full, 23, 13, 12341980689419102786),
            (R, U, N, One, 23, 13, 11531341460376497713),
            (R, U, N, Full, 23, 13, 13486171366216007501),
            (R, U, T, One, 23, 13, 270026691668154310),
            (R, U, T, Full, 23, 13, 93710698049720575),
            (R, D, N, One, 23, 13, 17239598992857270912),
            (R, D, N, Full, 23, 13, 1772897595565387234),
            (R, D, T, One, 23, 13, 18389555409070902198),
            (R, D, T, Full, 23, 13, 1176540359500558824),
            (L, U, N, One, 389, 13, 6828093926903151287),
            (L, D, N, One, 389, 13, 9674655553529421790),
            (R, U, N, One, 23, 389, 17775131083870818462),
            (R, D, N, One, 23, 389, 9342470205600653733),
        ];
        for &(side, tri, trans, diag, m, n, expected) in tri_cases {
            let dim = if side == L { m } else { n };
            let a = data((dim + 1) * dim, 1);
            let mut b = data((m + 2) * n, 2);
            trmm(
                side,
                tri,
                trans,
                diag,
                m,
                n,
                0.7,
                &a,
                dim + 1,
                &mut b,
                m + 2,
            )
            .unwrap();
            let actual = hash(&b);
            if actual != expected {
                failures.push(format!(
                    "TRMM {side:?}/{tri:?}/{trans:?}/{diag:?} {m}/{n}: {actual} != {expected}"
                ));
            }
        }
        {
            let a = data(6, 1);
            let mut shared = data(16, 2);
            let before = shared.clone();
            let mut expected = shared[..2]
                .iter()
                .copied()
                .chain(shared[8..10].iter().copied())
                .collect::<Vec<_>>();
            gemm(
                N,
                N,
                2,
                2,
                3,
                0.7,
                &a,
                2,
                &before[3..],
                8,
                -0.3,
                &mut expected,
                2,
            )
            .unwrap();
            gemm_into_within(
                N,
                N,
                2,
                2,
                3,
                0.7,
                &mut shared,
                MatrixSource::Separate(&a),
                2,
                MatrixSource::Within(3),
                8,
                -0.3,
                0,
                8,
            )
            .unwrap();
            assert_eq!(&shared[..2], &expected[..2]);
            assert_eq!(&shared[8..10], &expected[2..]);
            assert_eq!(&shared[3..6], &before[3..6]);
            assert_eq!(&shared[11..14], &before[11..14]);
            gemm_within(
                N,
                N,
                0,
                2,
                3,
                0.,
                &mut [],
                usize::MAX,
                1,
                usize::MAX,
                3,
                1.,
                usize::MAX,
                1,
            )
            .unwrap();
        }
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}
