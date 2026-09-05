//! Source-backed RDKit quaternion alignment kernel.

const TOLERANCE: f64 = 1.0e-6;

pub(crate) type Transform3D = [[f64; 4]; 4];

pub(crate) fn identity() -> Transform3D {
    // BEGIN RDKIT CPP FUNCTION RDGeom::Transform3D::Transform3D (Transform3D.h:43-47)
    // RDKit✔️✔️:   Transform3D() : RDNumeric::SquareMatrix<double>(DIM_3D, 0.0) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < DIM_3D; i++) {
    // RDKit✔️✔️:       unsigned int id = i * (DIM_3D + 1);
    // RDKit✔️✔️:       d_data[id] = 1.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDGeom::Transform3D::Transform3D
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

pub(crate) fn transform_point(trans: &Transform3D, point: [f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP FUNCTION RDGeom::Transform3D::TransformPoint (Transform3D.cpp:29-36)
    // RDKit✔️✔️:   double x = data[0] * pt.x + data[1] * pt.y + data[2] * pt.z + data[3];
    // RDKit✔️✔️:   double y = data[4] * pt.x + data[5] * pt.y + data[6] * pt.z + data[7];
    // RDKit✔️✔️:   double z = data[8] * pt.x + data[9] * pt.y + data[10] * pt.z + data[11];
    // END RDKIT CPP FUNCTION RDGeom::Transform3D::TransformPoint
    [
        trans[0][0] * point[0] + trans[0][1] * point[1] + trans[0][2] * point[2] + trans[0][3],
        trans[1][0] * point[0] + trans[1][1] * point[1] + trans[1][2] * point[2] + trans[1][3],
        trans[2][0] * point[0] + trans[2][1] * point[1] + trans[2][2] * point[2] + trans[2][3],
    ]
}

pub(crate) fn set_translation(trans: &mut Transform3D, value: [f64; 3]) {
    // BEGIN RDKIT CPP FUNCTION RDGeom::Transform3D::SetTranslation (Transform3D.cpp:49-59)
    // RDKit✔️✔️:   data[i] = move.x;
    // RDKit✔️✔️:   data[i] = move.y;
    // RDKit✔️✔️:   data[i] = move.z;
    // END RDKIT CPP FUNCTION RDGeom::Transform3D::SetTranslation
    trans[0][3] = value[0];
    trans[1][3] = value[1];
    trans[2][3] = value[2];
}

pub(crate) fn mul(lhs: &Transform3D, rhs: &Transform3D) -> Transform3D {
    // BEGIN RDKIT CPP FUNCTION RDGeom::operator* (Transform3D.h:103-109)
    // RDKit✔️✔️:  t3 = t1*t2;
    // RDKit✔️✔️:  t3(point) = t1(t2(point));
    // END RDKIT CPP FUNCTION RDGeom::operator*
    let mut out = [[0.0; 4]; 4];
    for row in 0..4 {
        for col in 0..4 {
            out[row][col] = (0..4).map(|k| lhs[row][k] * rhs[k][col]).sum();
        }
    }
    out
}

fn set_rotation_from_quaternion(trans: &mut Transform3D, q: [f64; 4]) {
    // BEGIN RDKIT CPP FUNCTION RDGeom::Transform3D::SetRotationFromQuaternion (Transform3D.cpp:114-139)
    // RDKit✔️✔️:   double q00 = quaternion[0] * quaternion[0];
    // RDKit✔️✔️:   double q11 = quaternion[1] * quaternion[1];
    // RDKit✔️✔️:   double q22 = quaternion[2] * quaternion[2];
    // RDKit✔️✔️:   double q33 = quaternion[3] * quaternion[3];
    // RDKit✔️✔️:   double sumSq = q00 + q11 + q22 + q33;
    // END RDKIT CPP FUNCTION RDGeom::Transform3D::SetRotationFromQuaternion
    let q00 = q[0] * q[0];
    let q11 = q[1] * q[1];
    let q22 = q[2] * q[2];
    let q33 = q[3] * q[3];
    let sum_sq = q00 + q11 + q22 + q33;
    let q01 = 2.0 * q[0] * q[1];
    let q02 = 2.0 * q[0] * q[2];
    let q03 = 2.0 * q[0] * q[3];
    let q12 = 2.0 * q[1] * q[2];
    let q13 = 2.0 * q[1] * q[3];
    let q23 = 2.0 * q[2] * q[3];
    *trans = identity();
    trans[0][0] = (q00 + q11 - q22 - q33) / sum_sq;
    trans[0][1] = (q12 + q03) / sum_sq;
    trans[0][2] = (q13 - q02) / sum_sq;
    trans[1][0] = (q12 - q03) / sum_sq;
    trans[1][1] = (q00 - q11 + q22 - q33) / sum_sq;
    trans[1][2] = (q23 + q01) / sum_sq;
    trans[2][0] = (q13 + q02) / sum_sq;
    trans[2][1] = (q23 - q01) / sum_sq;
    trans[2][2] = (q00 - q11 - q22 + q33) / sum_sq;
}

fn reflect(trans: &mut Transform3D) {
    // BEGIN RDKIT CPP FUNCTION RDGeom::Transform3D::Reflect (Transform3D.cpp:141-146)
    // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j < 3; j++) {
    // RDKit✔️✔️:       d_data[i * DIM_3D + j] *= -1.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDGeom::Transform3D::Reflect
    for row in trans.iter_mut().take(3) {
        for cell in row.iter_mut().take(3) {
            *cell = -*cell;
        }
    }
}

fn weighted_sum(points: &[[f64; 3]], weights: Option<&[f64]>) -> [f64; 3] {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfPoints (AlignPoints.cpp:22-34)
    // RDKit✔️✔️:   for (pti = points.begin(); pti != points.end(); pti++) {
    // RDKit✔️✔️:     tmpPt = (*(*pti));
    // RDKit✔️✔️:     if (weights) { tmpPt *= wData[i]; }
    // RDKit✔️✔️:     res += tmpPt;
    // RDKit✔️✔️:     i++;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfPoints
    points.iter().enumerate().fold([0.0; 3], |mut sum, (i, p)| {
        let w = weights.map_or(1.0, |weights| weights[i]);
        sum[0] += w * p[0];
        sum[1] += w * p[1];
        sum[2] += w * p[2];
        sum
    })
}

fn weighted_len_sq(points: &[[f64; 3]], weights: Option<&[f64]>) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfLenSq (AlignPoints.cpp:43-49)
    // RDKit✔️✔️:     auto l = pti->lengthSq();
    // RDKit✔️✔️:     if (weights) { l *= wData[i]; }
    // RDKit✔️✔️:     res += l;
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfLenSq
    points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let w = weights.map_or(1.0, |weights| weights[i]);
            w * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
        })
        .sum()
}

fn covariance(ref_points: &[[f64; 3]], probe_points: &[[f64; 3]], weights: Option<&[f64]>) -> [[f64; 3]; 3] {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_computeCovarianceMat (AlignPoints.cpp:67-88)
    // RDKit✔️✔️:     double w = weights ? wData[i] : 1.0;
    // RDKit✔️✔️:     covMat[0][0] += w * (ppt->x) * (rpt->x);
    // RDKit✔️✔️:     covMat[0][1] += w * (ppt->x) * (rpt->y);
    // RDKit✔️✔️:     covMat[0][2] += w * (ppt->x) * (rpt->z);
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_computeCovarianceMat
    let mut out = [[0.0; 3]; 3];
    for (i, (r, p)) in ref_points.iter().zip(probe_points).enumerate() {
        let w = weights.map_or(1.0, |weights| weights[i]);
        for a in 0..3 {
            for b in 0..3 {
                out[a][b] += w * p[a] * r[b];
            }
        }
    }
    out
}

fn quad(c: [[f64; 3]; 3], r: [f64; 3], p: [f64; 3], w: f64) -> [[f64; 4]; 4] {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_covertCovMatToQuad (AlignPoints.cpp:97-129)
    // RDKit✔️✔️:   temp = pptSum.x / wtsSum;
    // RDKit✔️✔️:   PxRx = covMat[0][0] - temp * rptSum.x;
    // RDKit✔️✔️:   quad[0][0] = -2.0 * (PxRx + PyRy + PzRz);
    // RDKit✔️✔️:   quad[0][1] = quad[1][0] = 2.0 * (PyRz - PzRy);
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_covertCovMatToQuad
    let px_rx = c[0][0] - p[0] / w * r[0];
    let px_ry = c[0][1] - p[0] / w * r[1];
    let px_rz = c[0][2] - p[0] / w * r[2];
    let py_rx = c[1][0] - p[1] / w * r[0];
    let py_ry = c[1][1] - p[1] / w * r[1];
    let py_rz = c[1][2] - p[1] / w * r[2];
    let pz_rx = c[2][0] - p[2] / w * r[0];
    let pz_ry = c[2][1] - p[2] / w * r[1];
    let pz_rz = c[2][2] - p[2] / w * r[2];
    let mut q = [[0.0; 4]; 4];
    q[0][0] = -2.0 * (px_rx + py_ry + pz_rz);
    q[1][1] = -2.0 * (px_rx - py_ry - pz_rz);
    q[2][2] = -2.0 * (py_ry - pz_rz - px_rx);
    q[3][3] = -2.0 * (pz_rz - px_rx - py_ry);
    q[0][1] = 2.0 * (py_rz - pz_ry);
    q[1][0] = q[0][1];
    q[0][2] = 2.0 * (pz_rx - px_rz);
    q[2][0] = q[0][2];
    q[0][3] = 2.0 * (px_ry - py_rx);
    q[3][0] = q[0][3];
    q[1][2] = -2.0 * (px_ry + py_rx);
    q[2][1] = q[1][2];
    q[1][3] = -2.0 * (pz_rx + px_rz);
    q[3][1] = q[1][3];
    q[2][3] = -2.0 * (py_rz + pz_ry);
    q[3][2] = q[2][3];
    q
}

fn jacobi(mut a: [[f64; 4]; 4], max_iter: usize) -> ([f64; 4], [[f64; 4]; 4]) {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::jacobi (AlignPoints.cpp:163-270)
    // RDKit✔️✔️:   for (l = 0; l < maxIter; l++) {
    // RDKit✔️✔️:     if (fabs(diagNorm) > 1.e-16 && (offDiagNorm / diagNorm) <= TOLERANCE) {
    // RDKit✔️✔️:       goto Exit_now;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::jacobi
    let mut v = [[0.0; 4]; 4];
    let mut d = [0.0; 4];
    for j in 0..4 {
        v[j][j] = 1.0;
        d[j] = a[j][j];
    }
    for _ in 0..max_iter {
        let diag: f64 = d.iter().map(|x| x.abs()).sum();
        let mut off = 0.0;
        for j in 0..4 {
            for i in 0..j {
                off += a[i][j].abs();
            }
        }
        if diag.abs() > 1.0e-16 && off / diag <= TOLERANCE {
            break;
        }
        for j in 1..4 {
            for i in 0..j {
                let b = a[i][j];
                if b.abs() <= 0.0 {
                    continue;
                }
                let dma = d[j] - d[i];
                let t = if (dma.abs() + b.abs()) <= dma.abs() {
                    b / dma
                } else {
                    let q = 0.5 * dma / b;
                    let mut t = 1.0 / (q.abs() + (1.0 + q * q).sqrt());
                    if q < 0.0 {
                        t = -t;
                    }
                    t
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                a[i][j] = 0.0;
                for k in 0..i {
                    let x = c * a[k][i] - s * a[k][j];
                    a[k][j] = s * a[k][i] + c * a[k][j];
                    a[k][i] = x;
                }
                for k in (i + 1)..j {
                    let x = c * a[i][k] - s * a[k][j];
                    a[k][j] = s * a[i][k] + c * a[k][j];
                    a[i][k] = x;
                }
                for k in (j + 1)..4 {
                    let x = c * a[i][k] - s * a[j][k];
                    a[j][k] = s * a[i][k] + c * a[j][k];
                    a[i][k] = x;
                }
                for row in &mut v {
                    let x = c * row[i] - s * row[j];
                    row[j] = s * row[i] + c * row[j];
                    row[i] = x;
                }
                let x = c * c * d[i] + s * s * d[j] - 2.0 * c * s * b;
                d[j] = s * s * d[i] + c * c * d[j] + 2.0 * c * s * b;
                d[i] = x;
            }
        }
    }
    for j in 0..3 {
        let mut k = j;
        for i in (j + 1)..4 {
            if d[i] < d[k] {
                k = i;
            }
        }
        if k != j {
            d.swap(k, j);
            for row in &mut v {
                row.swap(k, j);
            }
        }
    }
    (d, v)
}

/// RDKit RDNumeric::Alignments::AlignPoints, including weighted and reflection paths.
pub(crate) fn align_points(
    ref_points: &[[f64; 3]],
    probe_points: &[[f64; 3]],
    weights: Option<&[f64]>,
    reflect_input: bool,
    max_iterations: usize,
) -> Result<(f64, Transform3D), &'static str> {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::AlignPoints
    // RDKit✔️✔️:   unsigned int npt = refPoints.size();
    // RDKit✔️✔️:   PRECONDITION(npt == probePoints.size(), "Mismatch in number of points");
    // RDKit✔️✔️:   trans.setToIdentity();
    // RDKit✔️✔️:   if (weights) { wtsSum = _sumOfWeights(*weights); }
    // RDKit✔️✔️:   else { wtsSum = static_cast<double>(npt); }
    // RDKit✔️✔️:   _computeCovarianceMat(refPoints, probePoints, weights, covMat);
    // RDKit✔️✔️:   if (reflect) { rptSum *= -1.0; reflectCovMat(covMat); }
    // RDKit✔️✔️:   _covertCovMatToQuad(covMat, rptSum, pptSum, wtsSum, quad);
    // RDKit✔️✔️:   jacobi(quad, eigenVals, eigenVecs, maxIterations);
    // RDKit✔️✔️:   trans.SetRotationFromQuaternion(quater);
    // RDKit✔️✔️:   if (reflect) { trans.Reflect(); }
    // RDKit✔️✔️:   double ssr = eigenVals[0] - (pptSum.lengthSq() + rptSum.lengthSq()) / wtsSum +
    // RDKit✔️✔️:                rptSumLenSq + pptSumLenSq;
    // RDKit✔️✔️:   if ((ssr < 0.0) && (fabs(ssr) < TOLERANCE)) { ssr = 0.0; }
    // RDKit✔️✔️:   trans.TransformPoint(pptSum);
    // RDKit✔️✔️:   move = rptSum; move -= pptSum; move /= wtsSum;
    // RDKit✔️✔️:   trans.SetTranslation(move);
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::AlignPoints
    if ref_points.len() != probe_points.len() {
        return Err("Mismatch in number of points");
    }
    if ref_points.is_empty() {
        return Err("alignment requires at least one point");
    }
    if let Some(weights) = weights {
        if weights.len() != ref_points.len() {
            return Err("Mismatch in number of points");
        }
        if weights.iter().any(|w| !(*w > 0.0)) {
            return Err("Negative weight specified for a point");
        }
    }
    let wsum = weights.map_or(ref_points.len() as f64, |w| w.iter().sum());
    let mut rsum = weighted_sum(ref_points, weights);
    let psum = weighted_sum(probe_points, weights);
    let rlen = weighted_len_sq(ref_points, weights);
    let plen = weighted_len_sq(probe_points, weights);
    let mut cov = covariance(ref_points, probe_points, weights);
    if reflect_input {
        for x in &mut rsum {
            *x = -*x;
        }
        for row in &mut cov {
            for x in row {
                *x = -*x;
            }
        }
    }
    let (evals, evecs) = jacobi(quad(cov, rsum, psum, wsum), max_iterations);
    let mut trans = identity();
    set_rotation_from_quaternion(&mut trans, [evecs[0][0], evecs[1][0], evecs[2][0], evecs[3][0]]);
    if reflect_input {
        reflect(&mut trans);
    }
    let mut ssr = evals[0] - (psum.iter().map(|x| x * x).sum::<f64>() + rsum.iter().map(|x| x * x).sum::<f64>()) / wsum
        + rlen
        + plen;
    if ssr < 0.0 && ssr.abs() < TOLERANCE {
        ssr = 0.0;
    }
    if reflect_input {
        for x in &mut rsum {
            *x = -*x;
        }
    }
    let moved = transform_point(&trans, psum);
    set_translation(
        &mut trans,
        [
            (rsum[0] - moved[0]) / wsum,
            (rsum[1] - moved[1]) / wsum,
            (rsum[2] - moved[2]) / wsum,
        ],
    );
    Ok((ssr, trans))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rigid_translation_has_zero_ssr() {
        let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let probe = [[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]];
        let (ssr, transform) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!(ssr.abs() < 1.0e-10, "ssr={ssr}");
        for (expected, point) in reference.iter().zip(probe) {
            let actual = transform_point(&transform, point);
            for axis in 0..3 {
                assert!((actual[axis] - expected[axis]).abs() < 1.0e-8);
            }
        }
    }

    #[test]
    fn rotation_uses_rdkit_transform_orientation() {
        let reference = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let probe = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let (ssr, transform) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!(ssr.abs() < 1.0e-10, "ssr={ssr}");
        for (expected, point) in reference.iter().zip(probe) {
            let actual = transform_point(&transform, point);
            for axis in 0..3 {
                assert!((actual[axis] - expected[axis]).abs() < 1.0e-8);
            }
        }
    }

    #[test]
    fn weighted_alignment_rejects_invalid_weights() {
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert_eq!(
            align_points(&points, &points, Some(&[1.0]), false, 50),
            Err("Mismatch in number of points")
        );
        assert_eq!(
            align_points(&points, &points, Some(&[1.0, 0.0]), false, 50),
            Err("Negative weight specified for a point")
        );
    }

    #[test]
    fn weighted_triangle_matches_rdkit_source_regression() {
        let c30 = std::f64::consts::FRAC_PI_6.cos();
        let s30 = std::f64::consts::FRAC_PI_6.sin();
        let reference = [[-c30, -s30, 0.0], [c30, -s30, 0.0], [0.0, 1.0, 0.0]];
        let probe = [
            [-2.0 * s30 + 3.0, 2.0 * c30, 4.0],
            [-2.0 * s30 + 3.0, -2.0 * c30, 4.0],
            [5.0, 0.0, 4.0],
        ];
        let (unweighted, _) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!((unweighted - 3.0).abs() < 1.0e-4);
        let (weighted, _) = align_points(&reference, &probe, Some(&[1.0, 1.0, 2.0]), false, 50).unwrap();
        assert!((weighted - 3.75).abs() < 1.0e-4);
        let (weighted, _) = align_points(&reference, &probe, Some(&[1.0, 2.0, 2.0]), false, 50).unwrap();
        assert!((weighted - 4.8).abs() < 1.0e-4);
    }

    #[test]
    fn reflection_matches_rdkit_source_regression() {
        let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let probe = [[2.0, 2.0, 3.0], [3.0, 2.0, 3.0], [2.0, 2.0, 4.0], [2.0, 3.0, 3.0]];
        let (without_reflection, _) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!((without_reflection - 1.0).abs() < 1.0e-4);
        let (with_reflection, transform) = align_points(&reference, &probe, None, true, 50).unwrap();
        assert!(with_reflection.abs() < 1.0e-4);
        for (expected, point) in reference.iter().zip(probe) {
            let actual = transform_point(&transform, point);
            for axis in 0..3 {
                assert!((actual[axis] - expected[axis]).abs() < 1.0e-4);
            }
        }
    }

    #[test]
    fn alignment_rejects_empty_and_mismatched_inputs() {
        assert_eq!(
            align_points(&[], &[], None, false, 50),
            Err("alignment requires at least one point")
        );
        assert_eq!(
            align_points(&[[0.0, 0.0, 0.0]], &[], None, false, 50),
            Err("Mismatch in number of points")
        );
    }
}
