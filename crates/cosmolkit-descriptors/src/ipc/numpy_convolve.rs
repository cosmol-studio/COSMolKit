//! NumPy2.5.3 polynomial-factor convolution and X86_V4 complex magnitude.
//! NumPy/OpenBLAS source notices are retained in LICENSE.
//! Only contiguous owned coefficient buffers and the [1,-root] kernel enter
//! these helpers, exactly as in the pinned Graphs.CharacteristicPolynomial.

/// DOUBLE_dot of one element followed by NumPy's outer accumulation.
fn dot_one(a: f64, b: f64) -> f64 {
    // NumPy❗✔️: @name@_dot(char *ip1, npy_intp is1, char *ip2, npy_intp is2, char *op,
    // NumPy❗✔️:            npy_intp n, void *NPY_UNUSED(ignore))
    // NumPy❗✔️: {
    // NumPy❗✔️: #if defined(HAVE_CBLAS)
    // NumPy❗✔️:     CBLAS_INT is1b = blas_stride(is1, sizeof(@type@));
    // NumPy❗✔️:     CBLAS_INT is2b = blas_stride(is2, sizeof(@type@));
    // NumPy❗✔️:
    // NumPy❗✔️:     if (is1b && is2b)
    // NumPy❗✔️:     {
    // NumPy❗✔️:         double sum = 0.;  /* double for stability */
    // NumPy❗✔️:
    // NumPy❗✔️:         while (n > 0) {
    // NumPy❗✔️:             CBLAS_INT chunk = n < NPY_CBLAS_CHUNK ? n : NPY_CBLAS_CHUNK;
    // NumPy❗✔️:
    // NumPy❗✔️:             sum += CBLAS_FUNC(cblas_@prefix@dot)(chunk,
    // NumPy❗✔️:                                      (@type@ *) ip1, is1b,
    // NumPy❗✔️:                                      (@type@ *) ip2, is2b);
    // NumPy❗✔️:             /* use char strides here */
    // NumPy❗✔️:             ip1 += chunk * is1;
    // NumPy❗✔️:             ip2 += chunk * is2;
    // NumPy❗✔️:             n -= chunk;
    // NumPy❗✔️:         }
    // NumPy❗✔️:         *((@type@ *)op) = (@type@)sum;
    // NumPy❗✔️:     }
    // NumPy❗✔️:     else
    // NumPy❗✔️: #endif
    // NumPy❗✔️:     {
    // NumPy❗✔️:         @type@ sum = (@type@)0;  /* could make this double */
    // NumPy❗✔️:         npy_intp i;
    // NumPy❗✔️:
    // NumPy❗✔️:         for (i = 0; i < n; i++, ip1 += is1, ip2 += is2) {
    // NumPy❗✔️:             const @type@ ip1r = *((@type@ *)ip1);
    // NumPy❗✔️:             const @type@ ip2r = *((@type@ *)ip2);
    // NumPy❗✔️:
    // NumPy❗✔️:             sum += ip1r * ip2r;
    // NumPy❗✔️:         }
    // NumPy❗✔️:         *((@type@ *)op) = sum;
    // NumPy❗✔️:     }
    // NumPy❗✔️: }
    // PyArray_Correlate requests CARRAY for both arrays: the Python reversed
    // kernel is copied to a positive-stride contiguous array before this call.
    // NumPy's double sum starts at +0 even for the single CBLAS chunk.
    0.0 + super::lapack_blas::dot(1, &[a], 1, &[b], 1)
}

/// Multiply a real polynomial by [1,-root], in the source's correlate order.
pub(super) fn real_factor(polynomial: &[f64], root: f64) -> Vec<f64> {
    // NumPy❗✔️: small_correlate(const char * d_, npy_intp dstride,
    // NumPy❗✔️:                 npy_intp nd, enum NPY_TYPES dtype,
    // NumPy❗✔️:                 const char * k_, npy_intp kstride,
    // NumPy❗✔️:                 npy_intp nk, enum NPY_TYPES ktype,
    // NumPy❗✔️:                 char * out_, npy_intp ostride)
    // NumPy❗✔️: {
    // NumPy❗✔️:     /* only handle small kernels and uniform types */
    // NumPy❗✔️:     if (nk > 11 || dtype != ktype) {
    // NumPy❗✔️:         return 0;
    // NumPy❗✔️:     }
    // NumPy❗✔️:
    // NumPy❗✔️:     switch (dtype) {
    // NumPy❗✔️: /**begin repeat
    // NumPy❗✔️:  * Float types
    // NumPy❗✔️:  *  #type = npy_float, npy_double#
    // NumPy❗✔️:  *  #TYPE = NPY_FLOAT, NPY_DOUBLE#
    // NumPy❗✔️:  */
    // NumPy❗✔️:         case @TYPE@:
    // NumPy❗✔️:             {
    // NumPy❗✔️:                 npy_intp i;
    // NumPy❗✔️:                 const @type@ * d = (@type@*)d_;
    // NumPy❗✔️:                 const @type@ * k = (@type@*)k_;
    // NumPy❗✔️:                 @type@ * out = (@type@*)out_;
    // NumPy❗✔️:                 dstride /= sizeof(@type@);
    // NumPy❗✔️:                 kstride /= sizeof(@type@);
    // NumPy❗✔️:                 ostride /= sizeof(@type@);
    // NumPy❗✔️:                 /* unroll inner loop to optimize register usage of the kernel*/
    // NumPy❗✔️:                 switch (nk) {
    // NumPy❗✔️: /**begin repeat1
    // NumPy❗✔️:  *  #ksz_outer = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11# */
    // NumPy❗✔️:                     case @ksz_outer@:
    // NumPy❗✔️:                     {
    // NumPy❗✔️: /**begin repeat2
    // NumPy❗✔️:  *  #ksz = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11# */
    // NumPy❗✔️: #if @ksz@ <= @ksz_outer@
    // NumPy❗✔️:                         /* load kernel */
    // NumPy❗✔️:                         const @type@ k@ksz@ = k[(@ksz@ - 1) * kstride];
    // NumPy❗✔️: #endif
    // NumPy❗✔️: /**end repeat2**/
    // NumPy❗✔️:                         for (i = 0; i < nd; i++) {
    // NumPy❗✔️:                             @type@ s = 0;
    // NumPy❗✔️: /**begin repeat2
    // NumPy❗✔️:  *  #ksz = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11# */
    // NumPy❗✔️: #if @ksz@ <= @ksz_outer@
    // NumPy❗✔️:                             s += d[(i + @ksz@ - 1) * dstride] * k@ksz@;
    // NumPy❗✔️: #endif
    // NumPy❗✔️: /**end repeat2**/
    // NumPy❗✔️:                             out[i * ostride] = s;
    // NumPy❗✔️:                         }
    // NumPy❗✔️:                         return 1;
    // NumPy❗✔️:                     }
    // NumPy❗✔️: /**end repeat1**/
    // NumPy❗✔️:                     default:
    // NumPy❗✔️:                         return 0;
    // NumPy❗✔️:                 }
    // NumPy❗✔️:             }
    // NumPy❗✔️: /**end repeat**/
    // NumPy❗✔️:         default:
    // NumPy❗✔️:             return 0;
    // NumPy❗✔️:     }
    // NumPy❗✔️: }
    assert!(
        !polynomial.is_empty(),
        "source polynomial starts with one coefficient"
    );
    let mut out = vec![0.; polynomial.len() + 1];
    let negative = -root;
    if polynomial.len() == 1 {
        // convolve swaps arguments when its kernel is longer than the input.
        out[0] = dot_one(1.0, polynomial[0]);
        out[1] = dot_one(negative, polynomial[0]);
    } else {
        out[0] = dot_one(polynomial[0], 1.0);
        for i in 1..polynomial.len() {
            // small_correlate's DOUBLE kernel2 uses separate multiply/adds,
            // unlike BLAS dot's FMA. Do not simplify the 0 or unit products:
            // their signed-zero/Inf/NaN operations are source behavior.
            let mut sum = 0.0;
            sum += polynomial[i - 1] * negative;
            sum += polynomial[i] * 1.0;
            out[i] = sum;
        }
        out[polynomial.len()] = dot_one(polynomial[polynomial.len() - 1], negative);
    }
    out
}

/// CDOUBLE_dot with length one/two and unit positive strides.
fn complex_dot<const N: usize>(a: [[f64; 2]; N], b: [[f64; 2]; N]) -> [f64; 2] {
    // OpenBLAS❗✔️: static void zdot_compute (BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y,OPENBLAS_COMPLEX_FLOAT *result)
    // OpenBLAS❗✔️: {
    // OpenBLAS❗✔️: 	BLASLONG i;
    // OpenBLAS❗✔️: 	BLASLONG ix,iy;
    // OpenBLAS❗✔️: 	FLOAT  dot[4] = { 0.0, 0.0, 0.0 , 0.0 } ;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if ( n <= 0 )
    // OpenBLAS❗✔️: 	{
    // OpenBLAS❗✔️: 		OPENBLAS_COMPLEX_FLOAT res=OPENBLAS_MAKE_COMPLEX_FLOAT(0.0,0.0);
    // OpenBLAS❗✔️: 		*result=res;
    // OpenBLAS❗✔️: 		return;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	if ( (inc_x == 1) && (inc_y == 1) )
    // OpenBLAS❗✔️: 	{
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 		BLASLONG n1 = n & -8;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 		if ( n1 )
    // OpenBLAS❗✔️: 			zdot_kernel_8(n1, x, y , dot );
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 		i = n1;
    // OpenBLAS❗✔️: 		BLASLONG j = i * 2;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 		while( i < n )
    // OpenBLAS❗✔️: 		{
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 			dot[0] += x[j]   * y[j]   ;
    // OpenBLAS❗✔️: 			dot[1] += x[j+1] * y[j+1] ;
    // OpenBLAS❗✔️: 			dot[2] += x[j]   * y[j+1] ;
    // OpenBLAS❗✔️: 			dot[3] += x[j+1] * y[j]   ;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 			j+=2;
    // OpenBLAS❗✔️: 			i++ ;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 		}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️: 	else
    // OpenBLAS❗✔️: 	{
    // OpenBLAS❗✔️: 		i=0;
    // OpenBLAS❗✔️: 		ix=0;
    // OpenBLAS❗✔️: 		iy=0;
    // OpenBLAS❗✔️: 		inc_x *= 2;
    // OpenBLAS❗✔️: 		inc_y *= 2;
    // OpenBLAS❗✔️: 		while(i < n)
    // OpenBLAS❗✔️: 		{
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 			dot[0] += x[ix]   * y[iy]   ;
    // OpenBLAS❗✔️: 			dot[1] += x[ix+1] * y[iy+1] ;
    // OpenBLAS❗✔️: 			dot[2] += x[ix]   * y[iy+1] ;
    // OpenBLAS❗✔️: 			dot[3] += x[ix+1] * y[iy]   ;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 			ix  += inc_x ;
    // OpenBLAS❗✔️: 			iy  += inc_y ;
    // OpenBLAS❗✔️: 			i++ ;
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: 		}
    // OpenBLAS❗✔️: 	}
    // OpenBLAS❗✔️:
    // OpenBLAS❗✔️: #if !defined(CONJ)
    // OpenBLAS❗✔️: 	OPENBLAS_COMPLEX_FLOAT res=OPENBLAS_MAKE_COMPLEX_FLOAT(dot[0]-dot[1],dot[2]+dot[3]);
    // OpenBLAS❗✔️: #else
    // OpenBLAS❗✔️: 	OPENBLAS_COMPLEX_FLOAT res=OPENBLAS_MAKE_COMPLEX_FLOAT(dot[0]+dot[1],dot[2]-dot[3]);
    // OpenBLAS❗✔️: #endif
    // OpenBLAS❗✔️:         *result=res;
    // OpenBLAS❗✔️: 	return;
    // OpenBLAS❗✔️: }
    // NumPy❗✔️: @name@_dot(char *ip1, npy_intp is1, char *ip2, npy_intp is2,
    // NumPy❗✔️:            char *op, npy_intp n, void *NPY_UNUSED(ignore))
    // NumPy❗✔️: {
    // NumPy❗✔️: #if defined(HAVE_CBLAS)
    // NumPy❗✔️:     CBLAS_INT is1b = blas_stride(is1, sizeof(@ctype@));
    // NumPy❗✔️:     CBLAS_INT is2b = blas_stride(is2, sizeof(@ctype@));
    // NumPy❗✔️:
    // NumPy❗✔️:     if (is1b && is2b) {
    // NumPy❗✔️:         double sum[2] = {0., 0.};  /* double for stability */
    // NumPy❗✔️:
    // NumPy❗✔️:         while (n > 0) {
    // NumPy❗✔️:             CBLAS_INT chunk = n < NPY_CBLAS_CHUNK ? n : NPY_CBLAS_CHUNK;
    // NumPy❗✔️:             @type@ tmp[2];
    // NumPy❗✔️:
    // NumPy❗✔️:             CBLAS_FUNC(cblas_@prefix@dotu_sub)(
    // NumPy❗✔️:                     (CBLAS_INT)chunk, ip1, is1b, ip2, is2b, tmp);
    // NumPy❗✔️:             sum[0] += (double)tmp[0];
    // NumPy❗✔️:             sum[1] += (double)tmp[1];
    // NumPy❗✔️:             /* use char strides here */
    // NumPy❗✔️:             ip1 += chunk * is1;
    // NumPy❗✔️:             ip2 += chunk * is2;
    // NumPy❗✔️:             n -= chunk;
    // NumPy❗✔️:         }
    // NumPy❗✔️:         ((@type@ *)op)[0] = (@type@)sum[0];
    // NumPy❗✔️:         ((@type@ *)op)[1] = (@type@)sum[1];
    // NumPy❗✔️:     }
    // NumPy❗✔️:     else
    // NumPy❗✔️: #endif
    // NumPy❗✔️:     {
    // NumPy❗✔️:         @type@ sumr = (@type@)0.0;
    // NumPy❗✔️:         @type@ sumi = (@type@)0.0;
    // NumPy❗✔️:         npy_intp i;
    // NumPy❗✔️:
    // NumPy❗✔️:         for (i = 0; i < n; i++, ip1 += is1, ip2 += is2) {
    // NumPy❗✔️:             const @type@ ip1r = ((@type@ *)ip1)[0];
    // NumPy❗✔️:             const @type@ ip1i = ((@type@ *)ip1)[1];
    // NumPy❗✔️:             const @type@ ip2r = ((@type@ *)ip2)[0];
    // NumPy❗✔️:             const @type@ ip2i = ((@type@ *)ip2)[1];
    // NumPy❗✔️:
    // NumPy❗✔️:             sumr += ip1r * ip2r - ip1i * ip2i;
    // NumPy❗✔️:             sumi += ip1r * ip2i + ip1i * ip2r;
    // NumPy❗✔️:         }
    // NumPy❗✔️:         ((@type@ *)op)[0] = sumr;
    // NumPy❗✔️:         ((@type@ *)op)[1] = sumi;
    // NumPy❗✔️:     }
    // NumPy❗✔️: }
    // Only the unit-stride scalar tail is reached for this fixed kernel2.
    // The actual wheel contracts its four independent += products to FMA;
    // both the real/imag merge and NumPy's chunk accumulation stay separate.
    let mut dot = [0.0_f64; 4];
    for i in 0..N {
        dot[0] = a[i][0].mul_add(b[i][0], dot[0]);
        dot[1] = a[i][1].mul_add(b[i][1], dot[1]);
        dot[2] = a[i][0].mul_add(b[i][1], dot[2]);
        dot[3] = a[i][1].mul_add(b[i][0], dot[3]);
    }
    let re = dot[0] - dot[1];
    let im = dot[2] + dot[3];
    // OPENBLAS_MAKE_COMPLEX_FLOAT expands to real + I*imag on this build.
    // Actual zdot_compute uses an FMA for the real component; retain even
    // imag*0, which turns an infinite imaginary component into real NaN.
    [0.0 + im.mul_add(0.0, re), 0.0 + im]
}

/// Multiply a complex polynomial by [1,-root]; no early real projection.
pub(super) fn complex_factor(
    polynomial: &[[f64; 2]],
    root: [f64; 2],
    real_roots: bool,
) -> Vec<[f64; 2]> {
    // NumPy❗✔️: for r in roots[numpy.argsort(numpy.abs(roots))]:
    // NumPy❗✔️:   res = numpy.convolve(res, numpy.array([1.0, -r]))
    // NumPy❗✔️: return numpy.real(res)
    // small_correlate has no complex dtype case, so all points use CDOUBLE_dot.
    assert!(
        !polynomial.is_empty(),
        "source polynomial starts with one coefficient"
    );
    let mut out = vec![[0.; 2]; polynomial.len() + 1];
    let one = [1.0, 0.0];
    // When eigvals cast all-real roots to float64, the real [1,-r] factor
    // acquires +0 imaginary parts during common-dtype conversion. Complex
    // roots instead negate both components before convolution.
    let negative = [-root[0], if real_roots { 0.0 } else { -root[1] }];
    if polynomial.len() == 1 {
        out[0] = complex_dot([one], [polynomial[0]]);
        out[1] = complex_dot([negative], [polynomial[0]]);
    } else {
        out[0] = complex_dot([polynomial[0]], [one]);
        for i in 1..polynomial.len() {
            out[i] = complex_dot([polynomial[i - 1], polynomial[i]], [negative, one]);
        }
        out[polynomial.len()] = complex_dot([polynomial[polynomial.len() - 1]], [negative]);
    }
    out
}

/// One lane of the actual non-overlapping X86_V4 complex absolute ufunc.
pub(super) fn complex_absolute(value: [f64; 2]) -> f64 {
    // NumPy❗✔️: simd_cabsolute_@sfx@(npyv_@sfx@ re, npyv_@sfx@ im)
    // NumPy❗✔️: {
    // NumPy❗✔️:     const npyv_@sfx@ inf = npyv_setall_@sfx@(@INF@);
    // NumPy❗✔️:     const npyv_@sfx@ nan = npyv_setall_@sfx@(@NAN@);
    // NumPy❗✔️:
    // NumPy❗✔️:     re = npyv_abs_@sfx@(re);
    // NumPy❗✔️:     im = npyv_abs_@sfx@(im);
    // NumPy❗✔️:     /*
    // NumPy❗✔️:      * If real or imag = INF, then convert it to inf + j*inf
    // NumPy❗✔️:      * Handles: inf + j*nan, nan + j*inf
    // NumPy❗✔️:      */
    // NumPy❗✔️:     npyv_@bsfx@ re_infmask = npyv_cmpeq_@sfx@(re, inf);
    // NumPy❗✔️:     npyv_@bsfx@ im_infmask = npyv_cmpeq_@sfx@(im, inf);
    // NumPy❗✔️:     im = npyv_select_@sfx@(re_infmask, inf, im);
    // NumPy❗✔️:     re = npyv_select_@sfx@(im_infmask, inf, re);
    // NumPy❗✔️:     /*
    // NumPy❗✔️:      * If real or imag = NAN, then convert it to nan + j*nan
    // NumPy❗✔️:      * Handles: x + j*nan, nan + j*x
    // NumPy❗✔️:      */
    // NumPy❗✔️:     npyv_@bsfx@ re_nnanmask = npyv_notnan_@sfx@(re);
    // NumPy❗✔️:     npyv_@bsfx@ im_nnanmask = npyv_notnan_@sfx@(im);
    // NumPy❗✔️:     im = npyv_select_@sfx@(re_nnanmask, im, nan);
    // NumPy❗✔️:     re = npyv_select_@sfx@(im_nnanmask, re, nan);
    // NumPy❗✔️:
    // NumPy❗✔️:     npyv_@sfx@ larger  = npyv_max_@sfx@(re, im);
    // NumPy❗✔️:     npyv_@sfx@ smaller = npyv_min_@sfx@(im, re);
    // NumPy❗✔️:     /*
    // NumPy❗✔️:      * Calculate div_mask to prevent 0./0. and inf/inf operations in div
    // NumPy❗✔️:      */
    // NumPy❗✔️:     npyv_@bsfx@ zeromask = npyv_cmpeq_@sfx@(larger, npyv_zero_@sfx@());
    // NumPy❗✔️:     npyv_@bsfx@ infmask = npyv_cmpeq_@sfx@(smaller, inf);
    // NumPy❗✔️:     npyv_@bsfx@ div_mask = npyv_not_@bsfx@(npyv_or_@bsfx@(zeromask, infmask));
    // NumPy❗✔️:
    // NumPy❗✔️:     npyv_@sfx@ ratio = npyv_ifdivz_@sfx@(div_mask, smaller, larger);
    // NumPy❗✔️:     npyv_@sfx@ hypot = npyv_sqrt_@sfx@(
    // NumPy❗✔️:         npyv_muladd_@sfx@(ratio, ratio, npyv_setall_@sfx@(1.0@c@)
    // NumPy❗✔️:     ));
    // NumPy❗✔️:     return npyv_mul_@sfx@(hypot, larger);
    // NumPy❗✔️: }
    // The eigenvalue array and magnitude output are contiguous and disjoint.
    // Even a short tail uses the masked SIMD path, never scalar libc hypot.
    let [a, b] = value.map(f64::abs);
    if a.is_infinite() || b.is_infinite() {
        return f64::INFINITY;
    }
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    let hi = a.max(b);
    let lo = a.min(b);
    let ratio = if hi == 0.0 { 0.0 } else { lo / hi };
    ratio.mul_add(ratio, 1.0).sqrt() * hi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numpy_factor_convolution_retains_complex_and_nonfinite_arithmetic() {
        assert_eq!(real_factor(&[1., 2.], 3.), vec![1., -1., -6.]);
        assert_eq!(
            real_factor(&[1.], 0.)
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            vec![1.0_f64.to_bits(), 0.0_f64.to_bits()]
        );
        let first = complex_factor(&[[1., 0.]], [0., 1.], false);
        let second = complex_factor(&first, [0., -1.], false);
        assert_eq!(second, vec![[1., 0.], [0., 0.], [1., 0.]]);
        assert_eq!(complex_absolute([3., 4.]).to_bits(), 5.0_f64.to_bits());
        assert_eq!(complex_absolute([f64::NAN, f64::INFINITY]), f64::INFINITY);
        assert!(complex_absolute([f64::NAN, 1.]).is_nan());
        let overflow = complex_dot([[0., f64::MAX], [0., f64::MAX]], [[1., 0.], [1., 0.]]);
        assert!(overflow[0].is_nan());
        assert_eq!(overflow[1], f64::INFINITY);
    }
}
