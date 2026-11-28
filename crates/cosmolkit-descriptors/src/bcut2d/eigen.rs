// Source-derived portions from Eigen 3.3.4:
// Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010 Jitse Niesen <jitse@maths.leeds.ac.uk>
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2010 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
// Copyright (C) 2014 Pedro Gonnet (pedro.gonnet@gmail.com)
// Copyright (C) 2016 Gael Guennebaud <gael.guennebaud@inria.fr>
// SPDX-License-Identifier: MPL-2.0
//! Eigen 3.3.4's dynamic, real, column-major/lower self-adjoint solver path.
//!
//! Packet2d reductions are expressed as ordered scalar lanes on every target.
//! This is the Eigen version installed by the Linux RDKit 2026.09.1 wheel build.
//! This preserves the fixed SSE2 reference path rather than choosing a new
//! reduction order from the host's SIMD width. No Eigen ABI or dependency is used.
//! Only eigenvalues are observable: construction and rotation of Q are omitted.
//! The upstream diagonal/subdiagonal updates do not depend on Q; omitting it
//! preserves their arithmetic and avoids unused cubic work and matrix storage.

use super::burden::{allocate_dense_values, allocate_vector};
use crate::DescriptorError;

pub(super) struct Workspace {
    dimension: usize,
    matrix: Vec<f64>,
    diagonal: Vec<f64>,
    off_diagonal: Vec<f64>,
    reflector: Vec<f64>,
    update: Vec<f64>,
}

impl Workspace {
    pub(super) fn new(n: usize) -> Result<Self, DescriptorError> {
        Ok(Self {
            dimension: n,
            matrix: allocate_dense_values(n, 0.0)?,
            diagonal: allocate_vector(n)?,
            off_diagonal: allocate_vector(n.saturating_sub(1))?,
            reflector: allocate_vector(n)?,
            update: allocate_vector(n)?,
        })
    }
}

fn reduce_packet2(
    len: usize,
    f: impl Fn(usize) -> f64,
    op: impl Fn(f64, f64) -> f64,
    predux: impl Fn(f64, f64) -> f64,
) -> f64 {
    // Eigen❗✔️: struct redux_impl<Func, Derived, LinearVectorizedTraversal, NoUnrolling>
    // Eigen❗✔️: {
    // Eigen❗✔️:   typedef typename Derived::Scalar Scalar;
    // Eigen❗✔️:   typedef typename redux_traits<Func, Derived>::PacketType PacketScalar;
    // Eigen❗✔️:
    // Eigen❗✔️:   static Scalar run(const Derived &mat, const Func& func)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     const Index size = mat.size();
    // Eigen❗✔️:
    // Eigen❗✔️:     const Index packetSize = redux_traits<Func, Derived>::PacketSize;
    // Eigen❗✔️:     const int packetAlignment = unpacket_traits<PacketScalar>::alignment;
    // Eigen❗✔️:     enum {
    // Eigen❗✔️:       alignment0 = (bool(Derived::Flags & DirectAccessBit) && bool(packet_traits<Scalar>::AlignedOnScalar)) ? int(packetAlignment) : int(Unaligned),
    // Eigen❗✔️:       alignment = EIGEN_PLAIN_ENUM_MAX(alignment0, Derived::Alignment)
    // Eigen❗✔️:     };
    // Eigen❗✔️:     const Index alignedStart = internal::first_default_aligned(mat.nestedExpression());
    // Eigen❗✔️:     const Index alignedSize2 = ((size-alignedStart)/(2*packetSize))*(2*packetSize);
    // Eigen❗✔️:     const Index alignedSize = ((size-alignedStart)/(packetSize))*(packetSize);
    // Eigen❗✔️:     const Index alignedEnd2 = alignedStart + alignedSize2;
    // Eigen❗✔️:     const Index alignedEnd  = alignedStart + alignedSize;
    // Eigen❗✔️:     Scalar res;
    // Eigen❗✔️:     if(alignedSize)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       PacketScalar packet_res0 = mat.template packet<alignment,PacketScalar>(alignedStart);
    // Eigen❗✔️:       if(alignedSize>packetSize) // we have at least two packets to partly unroll the loop
    // Eigen❗✔️:       {
    // Eigen❗✔️:         PacketScalar packet_res1 = mat.template packet<alignment,PacketScalar>(alignedStart+packetSize);
    // Eigen❗✔️:         for(Index index = alignedStart + 2*packetSize; index < alignedEnd2; index += 2*packetSize)
    // Eigen❗✔️:         {
    // Eigen❗✔️:           packet_res0 = func.packetOp(packet_res0, mat.template packet<alignment,PacketScalar>(index));
    // Eigen❗✔️:           packet_res1 = func.packetOp(packet_res1, mat.template packet<alignment,PacketScalar>(index+packetSize));
    // Eigen❗✔️:         }
    // Eigen❗✔️:
    // Eigen❗✔️:         packet_res0 = func.packetOp(packet_res0,packet_res1);
    // Eigen❗✔️:         if(alignedEnd>alignedEnd2)
    // Eigen❗✔️:           packet_res0 = func.packetOp(packet_res0, mat.template packet<alignment,PacketScalar>(alignedEnd2));
    // Eigen❗✔️:       }
    // Eigen❗✔️:       res = func.predux(packet_res0);
    // Eigen❗✔️:
    // Eigen❗✔️:       for(Index index = 0; index < alignedStart; ++index)
    // Eigen❗✔️:         res = func(res,mat.coeff(index));
    // Eigen❗✔️:
    // Eigen❗✔️:       for(Index index = alignedEnd; index < size; ++index)
    // Eigen❗✔️:         res = func(res,mat.coeff(index));
    // Eigen❗✔️:     }
    // Eigen❗✔️:     else // too small to vectorize anything.
    // Eigen❗✔️:          // since this is dynamic-size hence inefficient anyway for such small sizes, don't try to optimize.
    // Eigen❗✔️:     {
    // Eigen❗✔️:       res = mat.coeff(0);
    // Eigen❗✔️:       for(Index index = 1; index < size; ++index)
    // Eigen❗✔️:         res = func(res,mat.coeff(index));
    // Eigen❗✔️:     }
    // Eigen❗✔️:
    // Eigen❗✔️:     return res;
    // Eigen❗✔️:   }
    // Eigen❗✔️: };
    // Eigen❗✔️:
    // Eigen❗✔️: // NOTE: for SliceVectorizedTraversal we simply bypass unrolling

    if len == 0 {
        return 0.0;
    }
    let end = len / 2 * 2;
    if end == 0 {
        return f(0);
    }
    let mut a = [f(0), f(1)];
    if end > 2 {
        let mut b = [f(2), f(3)];
        let end4 = len / 4 * 4;
        let mut i = 4;
        while i < end4 {
            a[0] = op(a[0], f(i));
            a[1] = op(a[1], f(i + 1));
            b[0] = op(b[0], f(i + 2));
            b[1] = op(b[1], f(i + 3));
            i += 4;
        }
        a[0] = op(a[0], b[0]);
        a[1] = op(a[1], b[1]);
        if end > end4 {
            a[0] = op(a[0], f(end4));
            a[1] = op(a[1], f(end4 + 1));
        }
    }
    let mut result = predux(a[0], a[1]);
    for i in end..len {
        result = op(result, f(i));
    }
    result
}

fn sum_packet2(len: usize, f: impl Fn(usize) -> f64) -> f64 {
    reduce_packet2(len, f, |a, b| a + b, |a, b| a + b)
}

fn symv(a: &[f64], n: usize, base: usize, v: &[f64], h: f64, w: &mut [f64], offset: usize) {
    // Eigen❗✔️: void selfadjoint_matrix_vector_product<Scalar,Index,StorageOrder,UpLo,ConjugateLhs,ConjugateRhs,Version>::run(
    // Eigen❗✔️:   Index size,
    // Eigen❗✔️:   const Scalar*  lhs, Index lhsStride,
    // Eigen❗✔️:   const Scalar*  rhs,
    // Eigen❗✔️:   Scalar* res,
    // Eigen❗✔️:   Scalar alpha)
    // Eigen❗✔️: {
    // Eigen❗✔️:   typedef typename packet_traits<Scalar>::type Packet;
    // Eigen❗✔️:   typedef typename NumTraits<Scalar>::Real RealScalar;
    // Eigen❗✔️:   const Index PacketSize = sizeof(Packet)/sizeof(Scalar);
    // Eigen❗✔️:
    // Eigen❗✔️:   enum {
    // Eigen❗✔️:     IsRowMajor = StorageOrder==RowMajor ? 1 : 0,
    // Eigen❗✔️:     IsLower = UpLo == Lower ? 1 : 0,
    // Eigen❗✔️:     FirstTriangular = IsRowMajor == IsLower
    // Eigen❗✔️:   };
    // Eigen❗✔️:
    // Eigen❗✔️:   conj_helper<Scalar,Scalar,NumTraits<Scalar>::IsComplex && EIGEN_LOGICAL_XOR(ConjugateLhs,  IsRowMajor), ConjugateRhs> cj0;
    // Eigen❗✔️:   conj_helper<Scalar,Scalar,NumTraits<Scalar>::IsComplex && EIGEN_LOGICAL_XOR(ConjugateLhs, !IsRowMajor), ConjugateRhs> cj1;
    // Eigen❗✔️:   conj_helper<RealScalar,Scalar,false, ConjugateRhs> cjd;
    // Eigen❗✔️:
    // Eigen❗✔️:   conj_helper<Packet,Packet,NumTraits<Scalar>::IsComplex && EIGEN_LOGICAL_XOR(ConjugateLhs,  IsRowMajor), ConjugateRhs> pcj0;
    // Eigen❗✔️:   conj_helper<Packet,Packet,NumTraits<Scalar>::IsComplex && EIGEN_LOGICAL_XOR(ConjugateLhs, !IsRowMajor), ConjugateRhs> pcj1;
    // Eigen❗✔️:
    // Eigen❗✔️:   Scalar cjAlpha = ConjugateRhs ? numext::conj(alpha) : alpha;
    // Eigen❗✔️:
    // Eigen❗✔️:   Index bound = numext::maxi(Index(0), size-8) & 0xfffffffe;
    // Eigen❗✔️:   if (FirstTriangular)
    // Eigen❗✔️:     bound = size - bound;
    // Eigen❗✔️:
    // Eigen❗✔️:   for (Index j=FirstTriangular ? bound : 0;
    // Eigen❗✔️:        j<(FirstTriangular ? size : bound);j+=2)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     const Scalar* EIGEN_RESTRICT A0 = lhs + j*lhsStride;
    // Eigen❗✔️:     const Scalar* EIGEN_RESTRICT A1 = lhs + (j+1)*lhsStride;
    // Eigen❗✔️:
    // Eigen❗✔️:     Scalar t0 = cjAlpha * rhs[j];
    // Eigen❗✔️:     Packet ptmp0 = pset1<Packet>(t0);
    // Eigen❗✔️:     Scalar t1 = cjAlpha * rhs[j+1];
    // Eigen❗✔️:     Packet ptmp1 = pset1<Packet>(t1);
    // Eigen❗✔️:
    // Eigen❗✔️:     Scalar t2(0);
    // Eigen❗✔️:     Packet ptmp2 = pset1<Packet>(t2);
    // Eigen❗✔️:     Scalar t3(0);
    // Eigen❗✔️:     Packet ptmp3 = pset1<Packet>(t3);
    // Eigen❗✔️:
    // Eigen❗✔️:     Index starti = FirstTriangular ? 0 : j+2;
    // Eigen❗✔️:     Index endi   = FirstTriangular ? j : size;
    // Eigen❗✔️:     Index alignedStart = (starti) + internal::first_default_aligned(&res[starti], endi-starti);
    // Eigen❗✔️:     Index alignedEnd = alignedStart + ((endi-alignedStart)/(PacketSize))*(PacketSize);
    // Eigen❗✔️:
    // Eigen❗✔️:     res[j]   += cjd.pmul(numext::real(A0[j]), t0);
    // Eigen❗✔️:     res[j+1] += cjd.pmul(numext::real(A1[j+1]), t1);
    // Eigen❗✔️:     if(FirstTriangular)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       res[j]   += cj0.pmul(A1[j],   t1);
    // Eigen❗✔️:       t3       += cj1.pmul(A1[j],   rhs[j]);
    // Eigen❗✔️:     }
    // Eigen❗✔️:     else
    // Eigen❗✔️:     {
    // Eigen❗✔️:       res[j+1] += cj0.pmul(A0[j+1],t0);
    // Eigen❗✔️:       t2 += cj1.pmul(A0[j+1], rhs[j+1]);
    // Eigen❗✔️:     }
    // Eigen❗✔️:
    // Eigen❗✔️:     for (Index i=starti; i<alignedStart; ++i)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       res[i] += cj0.pmul(A0[i], t0) + cj0.pmul(A1[i],t1);
    // Eigen❗✔️:       t2 += cj1.pmul(A0[i], rhs[i]);
    // Eigen❗✔️:       t3 += cj1.pmul(A1[i], rhs[i]);
    // Eigen❗✔️:     }
    // Eigen❗✔️:     // Yes this an optimization for gcc 4.3 and 4.4 (=> huge speed up)
    // Eigen❗✔️:     // gcc 4.2 does this optimization automatically.
    // Eigen❗✔️:     const Scalar* EIGEN_RESTRICT a0It  = A0  + alignedStart;
    // Eigen❗✔️:     const Scalar* EIGEN_RESTRICT a1It  = A1  + alignedStart;
    // Eigen❗✔️:     const Scalar* EIGEN_RESTRICT rhsIt = rhs + alignedStart;
    // Eigen❗✔️:           Scalar* EIGEN_RESTRICT resIt = res + alignedStart;
    // Eigen❗✔️:     for (Index i=alignedStart; i<alignedEnd; i+=PacketSize)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       Packet A0i = ploadu<Packet>(a0It);  a0It  += PacketSize;
    // Eigen❗✔️:       Packet A1i = ploadu<Packet>(a1It);  a1It  += PacketSize;
    // Eigen❗✔️:       Packet Bi  = ploadu<Packet>(rhsIt); rhsIt += PacketSize; // FIXME should be aligned in most cases
    // Eigen❗✔️:       Packet Xi  = pload <Packet>(resIt);
    // Eigen❗✔️:
    // Eigen❗✔️:       Xi    = pcj0.pmadd(A0i,ptmp0, pcj0.pmadd(A1i,ptmp1,Xi));
    // Eigen❗✔️:       ptmp2 = pcj1.pmadd(A0i,  Bi, ptmp2);
    // Eigen❗✔️:       ptmp3 = pcj1.pmadd(A1i,  Bi, ptmp3);
    // Eigen❗✔️:       pstore(resIt,Xi); resIt += PacketSize;
    // Eigen❗✔️:     }
    // Eigen❗✔️:     for (Index i=alignedEnd; i<endi; i++)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       res[i] += cj0.pmul(A0[i], t0) + cj0.pmul(A1[i],t1);
    // Eigen❗✔️:       t2 += cj1.pmul(A0[i], rhs[i]);
    // Eigen❗✔️:       t3 += cj1.pmul(A1[i], rhs[i]);
    // Eigen❗✔️:     }
    // Eigen❗✔️:
    // Eigen❗✔️:     res[j]   += alpha * (t2 + predux(ptmp2));
    // Eigen❗✔️:     res[j+1] += alpha * (t3 + predux(ptmp3));
    // Eigen❗✔️:   }
    // Eigen❗✔️:   for (Index j=FirstTriangular ? 0 : bound;j<(FirstTriangular ? bound : size);j++)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     const Scalar* EIGEN_RESTRICT A0 = lhs + j*lhsStride;
    // Eigen❗✔️:
    // Eigen❗✔️:     Scalar t1 = cjAlpha * rhs[j];
    // Eigen❗✔️:     Scalar t2(0);
    // Eigen❗✔️:     res[j] += cjd.pmul(numext::real(A0[j]), t1);
    // Eigen❗✔️:     for (Index i=FirstTriangular ? 0 : j+1; i<(FirstTriangular ? j : size); i++)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       res[i] += cj0.pmul(A0[i], t1);
    // Eigen❗✔️:       t2 += cj1.pmul(A0[i], rhs[i]);
    // Eigen❗✔️:     }
    // Eigen❗✔️:     res[j] += alpha * t2;
    // Eigen❗✔️:   }
    // Eigen❗✔️: }

    let size = v.len();
    w.fill(0.0);
    let bound = size.saturating_sub(8) & !1;
    let at = |r: usize, c: usize| a[(base + c) * n + base + r];
    for j in (0..bound).step_by(2) {
        let t0 = h * v[j];
        let t1 = h * v[j + 1];
        let mut t2 = 0.0;
        let mut t3 = 0.0;
        let mut p2 = [0.0; 2];
        let mut p3 = [0.0; 2];
        let start = j + 2;
        let aligned = (start + ((offset + start) % 2)).min(size);
        let end = aligned + (size - aligned) / 2 * 2;
        w[j] += at(j, j) * t0;
        w[j + 1] += at(j + 1, j + 1) * t1;
        w[j + 1] += at(j + 1, j) * t0;
        t2 += at(j + 1, j) * v[j + 1];
        for i in start..aligned {
            w[i] += at(i, j) * t0 + at(i, j + 1) * t1;
            t2 += at(i, j) * v[i];
            t3 += at(i, j + 1) * v[i];
        }
        for i in (aligned..end).step_by(2) {
            for lane in 0..2 {
                let row = i + lane;
                w[row] = at(row, j) * t0 + (at(row, j + 1) * t1 + w[row]);
                p2[lane] = at(row, j) * v[row] + p2[lane];
                p3[lane] = at(row, j + 1) * v[row] + p3[lane];
            }
        }
        for i in end..size {
            w[i] += at(i, j) * t0 + at(i, j + 1) * t1;
            t2 += at(i, j) * v[i];
            t3 += at(i, j + 1) * v[i];
        }
        w[j] += h * (t2 + (p2[0] + p2[1]));
        w[j + 1] += h * (t3 + (p3[0] + p3[1]));
    }
    for j in bound..size {
        let t1 = h * v[j];
        let mut t2 = 0.0;
        w[j] += at(j, j) * t1;
        for i in j + 1..size {
            w[i] += at(i, j) * t1;
            t2 += at(i, j) * v[i];
        }
        w[j] += h * t2;
    }
}

fn tridiagonalize(w: &mut Workspace) {
    // Eigen❗✔️: void MatrixBase<Derived>::makeHouseholder(
    // Eigen❗✔️:   EssentialPart& essential,
    // Eigen❗✔️:   Scalar& tau,
    // Eigen❗✔️:   RealScalar& beta) const
    // Eigen❗✔️: {
    // Eigen❗✔️:   using std::sqrt;
    // Eigen❗✔️:   using numext::conj;
    // Eigen❗✔️:
    // Eigen❗✔️:   EIGEN_STATIC_ASSERT_VECTOR_ONLY(EssentialPart)
    // Eigen❗✔️:   VectorBlock<const Derived, EssentialPart::SizeAtCompileTime> tail(derived(), 1, size()-1);
    // Eigen❗✔️:
    // Eigen❗✔️:   RealScalar tailSqNorm = size()==1 ? RealScalar(0) : tail.squaredNorm();
    // Eigen❗✔️:   Scalar c0 = coeff(0);
    // Eigen❗✔️:   const RealScalar tol = (std::numeric_limits<RealScalar>::min)();
    // Eigen❗✔️:
    // Eigen❗✔️:   if(tailSqNorm <= tol && numext::abs2(numext::imag(c0))<=tol)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     tau = RealScalar(0);
    // Eigen❗✔️:     beta = numext::real(c0);
    // Eigen❗✔️:     essential.setZero();
    // Eigen❗✔️:   }
    // Eigen❗✔️:   else
    // Eigen❗✔️:   {
    // Eigen❗✔️:     beta = sqrt(numext::abs2(c0) + tailSqNorm);
    // Eigen❗✔️:     if (numext::real(c0)>=RealScalar(0))
    // Eigen❗✔️:       beta = -beta;
    // Eigen❗✔️:     essential = tail / (c0 - beta);
    // Eigen❗✔️:     tau = conj((beta - c0) / beta);
    // Eigen❗✔️:   }
    // Eigen❗✔️: }
    // Eigen❗✔️: void tridiagonalization_inplace(MatrixType& matA, CoeffVectorType& hCoeffs);
    // Eigen❗✔️: }
    // Eigen❗✔️: struct selfadjoint_rank2_update_selector<Scalar,Index,UType,VType,Lower>
    // Eigen❗✔️: {
    // Eigen❗✔️:   static EIGEN_DEVICE_FUNC
    // Eigen❗✔️:   void run(Scalar* mat, Index stride, const UType& u, const VType& v, const Scalar& alpha)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     const Index size = u.size();
    // Eigen❗✔️:     for (Index i=0; i<size; ++i)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       Map<Matrix<Scalar,Dynamic,1> >(mat+stride*i+i, size-i) +=
    // Eigen❗✔️:                         (numext::conj(alpha) * numext::conj(u.coeff(i))) * v.tail(size-i)
    // Eigen❗✔️:                       + (alpha * numext::conj(v.coeff(i))) * u.tail(size-i);
    // Eigen❗✔️:     }
    // Eigen❗✔️:   }
    // Eigen❗✔️: };

    let n = w.dimension;
    for column in 0..n.saturating_sub(1) {
        let base = column + 1;
        let size = n - base;
        let begin = column * n + base;
        let c0 = w.matrix[begin];
        let tail = sum_packet2(size - 1, |k| {
            let x = w.matrix[begin + k + 1];
            x * x
        });
        let (h, beta) = if tail <= f64::MIN_POSITIVE {
            for k in 1..size {
                w.matrix[begin + k] = 0.0;
            }
            (0.0, c0)
        } else {
            let mut beta = (c0 * c0 + tail).sqrt();
            if c0 >= 0.0 {
                beta = -beta;
            }
            let denominator = c0 - beta;
            for k in 1..size {
                w.matrix[begin + k] /= denominator;
            }
            ((beta - c0) / beta, beta)
        };
        w.matrix[begin] = 1.0;
        w.reflector[..size].copy_from_slice(&w.matrix[begin..begin + size]);
        symv(
            &w.matrix,
            n,
            base,
            &w.reflector[..size],
            h,
            &mut w.update[..size],
            column,
        );
        let dot = sum_packet2(size, |k| w.update[k] * w.reflector[k]);
        let correction = h * -0.5 * dot;
        for k in 0..size {
            w.update[k] += correction * w.reflector[k];
        }
        for c in 0..size {
            let left = -1.0 * w.reflector[c];
            let right = -1.0 * w.update[c];
            for r in c..size {
                w.matrix[(base + c) * n + base + r] += left * w.update[r] + right * w.reflector[r];
            }
        }
        w.matrix[begin] = beta;
    }
    for i in 0..n {
        w.diagonal[i] = w.matrix[i * n + i];
    }
    for i in 0..n.saturating_sub(1) {
        w.off_diagonal[i] = w.matrix[i * n + i + 1];
    }
}

fn eigen_hypot(x: f64, y: f64) -> f64 {
    // Eigen✔️✔️:   static inline RealScalar run(const Scalar& x, const Scalar& y)
    // Eigen✔️✔️:   {
    // Eigen✔️✔️:     EIGEN_USING_STD_MATH(abs);
    // Eigen✔️✔️:     EIGEN_USING_STD_MATH(sqrt);
    // Eigen✔️✔️:     RealScalar _x = abs(x);
    // Eigen✔️✔️:     RealScalar _y = abs(y);
    // Eigen✔️✔️:     Scalar p, qp;
    // Eigen✔️✔️:     if(_x>_y)
    // Eigen✔️✔️:     {
    // Eigen✔️✔️:       p = _x;
    // Eigen✔️✔️:       qp = _y / p;
    // Eigen✔️✔️:     }
    // Eigen✔️✔️:     else
    // Eigen✔️✔️:     {
    // Eigen✔️✔️:       p = _y;
    // Eigen✔️✔️:       qp = _x / p;
    // Eigen✔️✔️:     }
    // Eigen✔️✔️:     if(p==RealScalar(0)) return RealScalar(0);
    // Eigen✔️✔️:     return p * sqrt(RealScalar(1) + qp*qp);
    // Eigen✔️✔️:   }

    let (x, y) = (x.abs(), y.abs());
    let (p, qp) = if x > y { (x, y / x) } else { (y, x / y) };
    if p == 0.0 {
        return 0.0;
    }
    p * (1.0 + qp * qp).sqrt()
}

fn givens(p: f64, q: f64) -> (f64, f64) {
    // Eigen❗✔️: void JacobiRotation<Scalar>::makeGivens(const Scalar& p, const Scalar& q, Scalar* r, internal::false_type)
    // Eigen❗✔️: {
    // Eigen❗✔️:   using std::sqrt;
    // Eigen❗✔️:   using std::abs;
    // Eigen❗✔️:   if(q==Scalar(0))
    // Eigen❗✔️:   {
    // Eigen❗✔️:     m_c = p<Scalar(0) ? Scalar(-1) : Scalar(1);
    // Eigen❗✔️:     m_s = Scalar(0);
    // Eigen❗✔️:     if(r) *r = abs(p);
    // Eigen❗✔️:   }
    // Eigen❗✔️:   else if(p==Scalar(0))
    // Eigen❗✔️:   {
    // Eigen❗✔️:     m_c = Scalar(0);
    // Eigen❗✔️:     m_s = q<Scalar(0) ? Scalar(1) : Scalar(-1);
    // Eigen❗✔️:     if(r) *r = abs(q);
    // Eigen❗✔️:   }
    // Eigen❗✔️:   else if(abs(p) > abs(q))
    // Eigen❗✔️:   {
    // Eigen❗✔️:     Scalar t = q/p;
    // Eigen❗✔️:     Scalar u = sqrt(Scalar(1) + numext::abs2(t));
    // Eigen❗✔️:     if(p<Scalar(0))
    // Eigen❗✔️:       u = -u;
    // Eigen❗✔️:     m_c = Scalar(1)/u;
    // Eigen❗✔️:     m_s = -t * m_c;
    // Eigen❗✔️:     if(r) *r = p * u;
    // Eigen❗✔️:   }
    // Eigen❗✔️:   else
    // Eigen❗✔️:   {
    // Eigen❗✔️:     Scalar t = p/q;
    // Eigen❗✔️:     Scalar u = sqrt(Scalar(1) + numext::abs2(t));
    // Eigen❗✔️:     if(q<Scalar(0))
    // Eigen❗✔️:       u = -u;
    // Eigen❗✔️:     m_s = -Scalar(1)/u;
    // Eigen❗✔️:     m_c = -t * m_s;
    // Eigen❗✔️:     if(r) *r = q * u;
    // Eigen❗✔️:   }
    // Eigen❗✔️:
    // Eigen❗✔️: }

    if q == 0.0 {
        (if p < 0.0 { -1.0 } else { 1.0 }, 0.0)
    } else if p == 0.0 {
        (0.0, if q < 0.0 { 1.0 } else { -1.0 })
    } else if p.abs() > q.abs() {
        let t = q / p;
        let mut u = (1.0 + t * t).sqrt();
        if p < 0.0 {
            u = -u;
        }
        let c = 1.0 / u;
        (c, -t * c)
    } else {
        let t = p / q;
        let mut u = (1.0 + t * t).sqrt();
        if q < 0.0 {
            u = -u;
        }
        let s = -1.0 / u;
        (-t * s, s)
    }
}

fn qr_step(d: &mut [f64], sd: &mut [f64], start: usize, end: usize) {
    // Eigen❗✔️: static void tridiagonal_qr_step(RealScalar* diag, RealScalar* subdiag, Index start, Index end, Scalar* matrixQ, Index n)
    // Eigen❗✔️: {
    // Eigen❗✔️:   using std::abs;
    // Eigen❗✔️:   RealScalar td = (diag[end-1] - diag[end])*RealScalar(0.5);
    // Eigen❗✔️:   RealScalar e = subdiag[end-1];
    // Eigen❗✔️:   // Note that thanks to scaling, e^2 or td^2 cannot overflow, however they can still
    // Eigen❗✔️:   // underflow thus leading to inf/NaN values when using the following commented code:
    // Eigen❗✔️: //   RealScalar e2 = numext::abs2(subdiag[end-1]);
    // Eigen❗✔️: //   RealScalar mu = diag[end] - e2 / (td + (td>0 ? 1 : -1) * sqrt(td*td + e2));
    // Eigen❗✔️:   // This explain the following, somewhat more complicated, version:
    // Eigen❗✔️:   RealScalar mu = diag[end];
    // Eigen❗✔️:   if(td==RealScalar(0))
    // Eigen❗✔️:     mu -= abs(e);
    // Eigen❗✔️:   else
    // Eigen❗✔️:   {
    // Eigen❗✔️:     RealScalar e2 = numext::abs2(subdiag[end-1]);
    // Eigen❗✔️:     RealScalar h = numext::hypot(td,e);
    // Eigen❗✔️:     if(e2==RealScalar(0)) mu -= (e / (td + (td>RealScalar(0) ? RealScalar(1) : RealScalar(-1)))) * (e / h);
    // Eigen❗✔️:     else                  mu -= e2 / (td + (td>RealScalar(0) ? h : -h));
    // Eigen❗✔️:   }
    // Eigen❗✔️:
    // Eigen❗✔️:   RealScalar x = diag[start] - mu;
    // Eigen❗✔️:   RealScalar z = subdiag[start];
    // Eigen❗✔️:   for (Index k = start; k < end; ++k)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     JacobiRotation<RealScalar> rot;
    // Eigen❗✔️:     rot.makeGivens(x, z);
    // Eigen❗✔️:
    // Eigen❗✔️:     // do T = G' T G
    // Eigen❗✔️:     RealScalar sdk = rot.s() * diag[k] + rot.c() * subdiag[k];
    // Eigen❗✔️:     RealScalar dkp1 = rot.s() * subdiag[k] + rot.c() * diag[k+1];
    // Eigen❗✔️:
    // Eigen❗✔️:     diag[k] = rot.c() * (rot.c() * diag[k] - rot.s() * subdiag[k]) - rot.s() * (rot.c() * subdiag[k] - rot.s() * diag[k+1]);
    // Eigen❗✔️:     diag[k+1] = rot.s() * sdk + rot.c() * dkp1;
    // Eigen❗✔️:     subdiag[k] = rot.c() * sdk - rot.s() * dkp1;
    // Eigen❗✔️:
    // Eigen❗✔️:
    // Eigen❗✔️:     if (k > start)
    // Eigen❗✔️:       subdiag[k - 1] = rot.c() * subdiag[k-1] - rot.s() * z;
    // Eigen❗✔️:
    // Eigen❗✔️:     x = subdiag[k];
    // Eigen❗✔️:
    // Eigen❗✔️:     if (k < end - 1)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       z = -rot.s() * subdiag[k+1];
    // Eigen❗✔️:       subdiag[k + 1] = rot.c() * subdiag[k+1];
    // Eigen❗✔️:     }
    // Eigen❗✔️:
    // Eigen❗✔️:     // apply the givens rotation to the unit matrix Q = Q * G
    // Eigen❗✔️:     if (matrixQ)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       // FIXME if StorageOrder == RowMajor this operation is not very efficient
    // Eigen❗✔️:       Map<Matrix<Scalar,Dynamic,Dynamic,StorageOrder> > q(matrixQ,n,n);
    // Eigen❗✔️:       q.applyOnTheRight(k,k+1,rot);
    // Eigen❗✔️:     }
    // Eigen❗✔️:   }
    // Eigen❗✔️: }

    let td = (d[end - 1] - d[end]) * 0.5;
    let e = sd[end - 1];
    let mut mu = d[end];
    if td == 0.0 {
        mu -= e.abs();
    } else {
        let e2 = e * e;
        let h = eigen_hypot(td, e);
        if e2 == 0.0 {
            mu -= (e / (td + if td > 0.0 { 1.0 } else { -1.0 })) * (e / h);
        } else {
            mu -= e2 / (td + if td > 0.0 { h } else { -h });
        }
    }
    let mut x = d[start] - mu;
    let mut z = sd[start];
    for k in start..end {
        let (c, s) = givens(x, z);
        let sdk = s * d[k] + c * sd[k];
        let next = s * sd[k] + c * d[k + 1];
        d[k] = c * (c * d[k] - s * sd[k]) - s * (c * sd[k] - s * d[k + 1]);
        d[k + 1] = s * sdk + c * next;
        sd[k] = c * sdk - s * next;
        if k > start {
            sd[k - 1] = c * sd[k - 1] - s * z;
        }
        x = sd[k];
        if k < end - 1 {
            z = -s * sd[k + 1];
            sd[k + 1] *= c;
        }
    }
}

fn diagonalize(w: &mut Workspace) -> bool {
    // Eigen❗✔️: ComputationInfo computeFromTridiagonal_impl(DiagType& diag, SubDiagType& subdiag, const Index maxIterations, bool computeEigenvectors, MatrixType& eivec)
    // Eigen❗✔️: {
    // Eigen❗✔️:   using std::abs;
    // Eigen❗✔️:
    // Eigen❗✔️:   ComputationInfo info;
    // Eigen❗✔️:   typedef typename MatrixType::Scalar Scalar;
    // Eigen❗✔️:
    // Eigen❗✔️:   Index n = diag.size();
    // Eigen❗✔️:   Index end = n-1;
    // Eigen❗✔️:   Index start = 0;
    // Eigen❗✔️:   Index iter = 0; // total number of iterations
    // Eigen❗✔️:
    // Eigen❗✔️:   typedef typename DiagType::RealScalar RealScalar;
    // Eigen❗✔️:   const RealScalar considerAsZero = (std::numeric_limits<RealScalar>::min)();
    // Eigen❗✔️:   const RealScalar precision = RealScalar(2)*NumTraits<RealScalar>::epsilon();
    // Eigen❗✔️:
    // Eigen❗✔️:   while (end>0)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     for (Index i = start; i<end; ++i)
    // Eigen❗✔️:       if (internal::isMuchSmallerThan(abs(subdiag[i]),(abs(diag[i])+abs(diag[i+1])),precision) || abs(subdiag[i]) <= considerAsZero)
    // Eigen❗✔️:         subdiag[i] = 0;
    // Eigen❗✔️:
    // Eigen❗✔️:     // find the largest unreduced block
    // Eigen❗✔️:     while (end>0 && subdiag[end-1]==RealScalar(0))
    // Eigen❗✔️:     {
    // Eigen❗✔️:       end--;
    // Eigen❗✔️:     }
    // Eigen❗✔️:     if (end<=0)
    // Eigen❗✔️:       break;
    // Eigen❗✔️:
    // Eigen❗✔️:     // if we spent too many iterations, we give up
    // Eigen❗✔️:     iter++;
    // Eigen❗✔️:     if(iter > maxIterations * n) break;
    // Eigen❗✔️:
    // Eigen❗✔️:     start = end - 1;
    // Eigen❗✔️:     while (start>0 && subdiag[start-1]!=0)
    // Eigen❗✔️:       start--;
    // Eigen❗✔️:
    // Eigen❗✔️:     internal::tridiagonal_qr_step<MatrixType::Flags&RowMajorBit ? RowMajor : ColMajor>(diag.data(), subdiag.data(), start, end, computeEigenvectors ? eivec.data() : (Scalar*)0, n);
    // Eigen❗✔️:   }
    // Eigen❗✔️:   if (iter <= maxIterations * n)
    // Eigen❗✔️:     info = Success;
    // Eigen❗✔️:   else
    // Eigen❗✔️:     info = NoConvergence;
    // Eigen❗✔️:
    // Eigen❗✔️:   // Sort eigenvalues and corresponding vectors.
    // Eigen❗✔️:   // TODO make the sort optional ?
    // Eigen❗✔️:   // TODO use a better sort algorithm !!
    // Eigen❗✔️:   if (info == Success)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     for (Index i = 0; i < n-1; ++i)
    // Eigen❗✔️:     {
    // Eigen❗✔️:       Index k;
    // Eigen❗✔️:       diag.segment(i,n-i).minCoeff(&k);
    // Eigen❗✔️:       if (k > 0)
    // Eigen❗✔️:       {
    // Eigen❗✔️:         std::swap(diag[i], diag[k+i]);
    // Eigen❗✔️:         if(computeEigenvectors)
    // Eigen❗✔️:           eivec.col(i).swap(eivec.col(k+i));
    // Eigen❗✔️:       }
    // Eigen❗✔️:     }
    // Eigen❗✔️:   }
    // Eigen❗✔️:   return info;
    // Eigen❗✔️: }

    let n = w.dimension;
    let mut end = n - 1;
    let mut start = 0;
    let mut iter = 0;
    let precision = 2.0 * f64::EPSILON;
    while end > 0 {
        for i in start..end {
            if w.off_diagonal[i].abs()
                <= (w.diagonal[i].abs() + w.diagonal[i + 1].abs()) * precision
                || w.off_diagonal[i].abs() <= f64::MIN_POSITIVE
            {
                w.off_diagonal[i] = 0.0;
            }
        }
        while end > 0 && w.off_diagonal[end - 1] == 0.0 {
            end -= 1;
        }
        if end == 0 {
            break;
        }
        iter += 1;
        if iter > n.saturating_mul(30) {
            return false;
        }
        start = end - 1;
        while start > 0 && w.off_diagonal[start - 1] != 0.0 {
            start -= 1;
        }
        qr_step(&mut w.diagonal, &mut w.off_diagonal, start, end);
    }
    for i in 0..n - 1 {
        let mut k = i;
        for j in i + 1..n {
            if w.diagonal[j] < w.diagonal[k] {
                k = j;
            }
        }
        if k > i {
            w.diagonal.swap(i, k);
        }
    }
    true
}

pub(super) fn extrema(matrix: &[f64], w: &mut Workspace) -> Result<(f64, f64), DescriptorError> {
    // Eigen❗✔️: ::compute(const EigenBase<InputType>& a_matrix, int options)
    // Eigen❗✔️: {
    // Eigen❗✔️:   check_template_parameters();
    // Eigen❗✔️:
    // Eigen❗✔️:   const InputType &matrix(a_matrix.derived());
    // Eigen❗✔️:
    // Eigen❗✔️:   using std::abs;
    // Eigen❗✔️:   eigen_assert(matrix.cols() == matrix.rows());
    // Eigen❗✔️:   eigen_assert((options&~(EigVecMask|GenEigMask))==0
    // Eigen❗✔️:           && (options&EigVecMask)!=EigVecMask
    // Eigen❗✔️:           && "invalid option parameter");
    // Eigen❗✔️:   bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;
    // Eigen❗✔️:   Index n = matrix.cols();
    // Eigen❗✔️:   m_eivalues.resize(n,1);
    // Eigen❗✔️:
    // Eigen❗✔️:   if(n==1)
    // Eigen❗✔️:   {
    // Eigen❗✔️:     m_eivec = matrix;
    // Eigen❗✔️:     m_eivalues.coeffRef(0,0) = numext::real(m_eivec.coeff(0,0));
    // Eigen❗✔️:     if(computeEigenvectors)
    // Eigen❗✔️:       m_eivec.setOnes(n,n);
    // Eigen❗✔️:     m_info = Success;
    // Eigen❗✔️:     m_isInitialized = true;
    // Eigen❗✔️:     m_eigenvectorsOk = computeEigenvectors;
    // Eigen❗✔️:     return *this;
    // Eigen❗✔️:   }
    // Eigen❗✔️:
    // Eigen❗✔️:   // declare some aliases
    // Eigen❗✔️:   RealVectorType& diag = m_eivalues;
    // Eigen❗✔️:   EigenvectorsType& mat = m_eivec;
    // Eigen❗✔️:
    // Eigen❗✔️:   // map the matrix coefficients to [-1:1] to avoid over- and underflow.
    // Eigen❗✔️:   mat = matrix.template triangularView<Lower>();
    // Eigen❗✔️:   RealScalar scale = mat.cwiseAbs().maxCoeff();
    // Eigen❗✔️:   if(scale==RealScalar(0)) scale = RealScalar(1);
    // Eigen❗✔️:   mat.template triangularView<Lower>() /= scale;
    // Eigen❗✔️:   m_subdiag.resize(n-1);
    // Eigen❗✔️:   internal::tridiagonalization_inplace(mat, diag, m_subdiag, computeEigenvectors);
    // Eigen❗✔️:
    // Eigen❗✔️:   m_info = internal::computeFromTridiagonal_impl(diag, m_subdiag, m_maxIterations, computeEigenvectors, m_eivec);
    // Eigen❗✔️:
    // Eigen❗✔️:   // scale back the eigen values
    // Eigen❗✔️:   m_eivalues *= scale;
    // Eigen❗✔️:
    // Eigen❗✔️:   m_isInitialized = true;
    // Eigen❗✔️:   m_eigenvectorsOk = computeEigenvectors;
    // Eigen❗✔️:   return *this;
    // Eigen❗✔️: }

    let n = w.dimension;
    if n == 0 {
        return Ok((0.0, 0.0));
    }
    if n == 1 {
        return Ok((matrix[0], matrix[0]));
    }
    for c in 0..n {
        for r in 0..n {
            let v = if r >= c { matrix[r * n + c] } else { 0.0 };
            w.matrix[c * n + r] = v;
        }
    }
    let mut scale = reduce_packet2(
        n * n,
        |k| w.matrix[k].abs(),
        |a, b| if a < b { b } else { a },
        |a, b| if a > b { a } else { b },
    );
    if scale == 0.0 {
        scale = 1.0;
    }
    for c in 0..n {
        for r in c..n {
            w.matrix[c * n + r] /= scale;
        }
    }
    tridiagonalize(w);
    diagonalize(w);
    for value in &mut w.diagonal {
        *value *= scale;
    }
    Ok((w.diagonal[n - 1], w.diagonal[0]))
}
