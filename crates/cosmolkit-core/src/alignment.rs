//! Source-backed RDKit quaternion alignment kernel.

const TOLERANCE: f64 = 1.0e-6;

type Transform3D = [[f64; 4]; 4];

fn identity() -> Transform3D {
    let values = crate::Transform3D::identity();
    std::array::from_fn(|row| std::array::from_fn(|column| values.values()[4 * row + column]))
}

/// Apply a detached alignment matrix through the existing transform owner.
pub fn alignment_transform_point(trans: &Transform3D, point: [f64; 3]) -> [f64; 3] {
    crate::Transform3D::from_values(std::array::from_fn(|i| trans[i / 4][i % 4]))
        .transform_point(point)
}

pub(crate) fn set_translation(trans: &mut Transform3D, value: [f64; 3]) {
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: void Transform3D::SetTranslation(const Point3D &move) {
    // RDKit✔️✔️:   unsigned int i = DIM_3D - 1;
    // RDKit✔️✔️:   double *data = d_data.get();
    // RDKit✔️✔️:   data[i] = move.x;
    // RDKit✔️✔️:   i += DIM_3D;
    // RDKit✔️✔️:   data[i] = move.y;
    // RDKit✔️✔️:   i += DIM_3D;
    // RDKit✔️✔️:   data[i] = move.z;
    // RDKit✔️✔️:   i += DIM_3D;
    // RDKit✔️✔️:   data[i] = 1.0;
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

    trans[0][3] = value[0];
    trans[1][3] = value[1];
    trans[2][3] = value[2];
}

fn set_rotation_from_quaternion(trans: &mut Transform3D, q: [f64; 4]) {
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: void Transform3D::SetRotationFromQuaternion(const double quaternion[4]) {
    // RDKit✔️✔️:   double q00 = quaternion[0] * quaternion[0];
    // RDKit✔️✔️:   double q11 = quaternion[1] * quaternion[1];
    // RDKit✔️✔️:   double q22 = quaternion[2] * quaternion[2];
    // RDKit✔️✔️:   double q33 = quaternion[3] * quaternion[3];
    // RDKit✔️✔️:   double sumSq = q00 + q11 + q22 + q33;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double q01 = 2 * quaternion[0] * quaternion[1];
    // RDKit✔️✔️:   double q02 = 2 * quaternion[0] * quaternion[2];
    // RDKit✔️✔️:   double q03 = 2 * quaternion[0] * quaternion[3];
    // RDKit✔️✔️:   double q12 = 2 * quaternion[1] * quaternion[2];
    // RDKit✔️✔️:   double q13 = 2 * quaternion[1] * quaternion[3];
    // RDKit✔️✔️:   double q23 = 2 * quaternion[2] * quaternion[3];
    // RDKit✔️✔️:   double *data = d_data.get();
    // RDKit✔️✔️:   data[0] = (q00 + q11 - q22 - q33) / sumSq;
    // RDKit✔️✔️:   data[1] = (q12 + q03) / sumSq;
    // RDKit✔️✔️:   data[2] = (q13 - q02) / sumSq;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   data[4] = (q12 - q03) / sumSq;
    // RDKit✔️✔️:   data[5] = (q00 - q11 + q22 - q33) / sumSq;
    // RDKit✔️✔️:   data[6] = (q23 + q01) / sumSq;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   data[8] = (q13 + q02) / sumSq;
    // RDKit✔️✔️:   data[9] = (q23 - q01) / sumSq;
    // RDKit✔️✔️:   data[10] = (q00 - q11 - q22 + q33) / sumSq;
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

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
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: void Transform3D::Reflect() {
    // RDKit✔️✔️:   double *data = d_data.get();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < DIM_3D - 1; i++) {
    // RDKit✔️✔️:     unsigned int id = i * DIM_3D;
    // RDKit✔️✔️:     for (unsigned int j = 0; j < DIM_3D - 1; j++) {
    // RDKit✔️✔️:       data[id + j] *= -1.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

    for row in trans.iter_mut().take(3) {
        for cell in row.iter_mut().take(3) {
            *cell = -*cell;
        }
    }
}

fn weighted_sum(points: &[[f64; 3]], weights: Option<&[f64]>) -> [f64; 3] {
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: RDGeom::Point3D _weightedSumOfPoints(const RDGeom::Point3DConstPtrVect &points,
    // RDKit✔️✔️:                                      const DoubleVector *weights) {
    // RDKit✔️✔️:   PRECONDITION(!weights || points.size() == weights->size(), "");
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect_CI pti;
    // RDKit✔️✔️:   RDGeom::Point3D tmpPt, res;
    // RDKit✔️✔️:   const double *wData = weights ? weights->getData() : nullptr;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   for (pti = points.begin(); pti != points.end(); pti++) {
    // RDKit✔️✔️:     tmpPt = (*(*pti));
    // RDKit✔️✔️:     if (weights) {
    // RDKit✔️✔️:       tmpPt *= wData[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += tmpPt;
    // RDKit✔️✔️:     i++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

    points.iter().enumerate().fold([0.0; 3], |mut sum, (i, p)| {
        let w = weights.map_or(1.0, |weights| weights[i]);
        sum[0] += w * p[0];
        sum[1] += w * p[1];
        sum[2] += w * p[2];
        sum
    })
}

fn weighted_len_sq(points: &[[f64; 3]], weights: Option<&[f64]>) -> f64 {
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: double _weightedSumOfLenSq(const RDGeom::Point3DConstPtrVect &points,
    // RDKit✔️✔️:                            const DoubleVector *weights) {
    // RDKit✔️✔️:   PRECONDITION(!weights || (points.size() == weights->size()), "");
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   const double *wData = weights ? weights->getData() : nullptr;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   for (const auto &pti : points) {
    // RDKit✔️✔️:     auto l = pti->lengthSq();
    // RDKit✔️✔️:     if (weights) {
    // RDKit✔️✔️:       l *= wData[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += l;
    // RDKit✔️✔️:     i++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

    points
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let w = weights.map_or(1.0, |weights| weights[i]);
            w * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
        })
        .sum()
}

fn covariance(
    ref_points: &[[f64; 3]],
    probe_points: &[[f64; 3]],
    weights: Option<&[f64]>,
) -> [[f64; 3]; 3] {
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: void _computeCovarianceMat(const RDGeom::Point3DConstPtrVect &refPoints,
    // RDKit✔️✔️:                            const RDGeom::Point3DConstPtrVect &probePoints,
    // RDKit✔️✔️:                            const DoubleVector *weights, double covMat[3][3]) {
    // RDKit✔️✔️:   memset(static_cast<void *>(covMat), 0, 9 * sizeof(double));
    // RDKit✔️✔️:   unsigned int npt = refPoints.size();
    // RDKit✔️✔️:   CHECK_INVARIANT(npt == probePoints.size(), "Number of points mismatch");
    // RDKit✔️✔️:   CHECK_INVARIANT(!weights || (npt == weights->size()),
    // RDKit✔️✔️:                   "Number of points and number of weights do not match");
    // RDKit✔️✔️:   const double *wData = weights ? weights->getData() : nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const RDGeom::Point3D *rpt, *ppt;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < npt; i++) {
    // RDKit✔️✔️:     rpt = refPoints[i];
    // RDKit✔️✔️:     ppt = probePoints[i];
    // RDKit✔️✔️:     double w = weights ? wData[i] : 1.0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     covMat[0][0] += w * (ppt->x) * (rpt->x);
    // RDKit✔️✔️:     covMat[0][1] += w * (ppt->x) * (rpt->y);
    // RDKit✔️✔️:     covMat[0][2] += w * (ppt->x) * (rpt->z);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     covMat[1][0] += w * (ppt->y) * (rpt->x);
    // RDKit✔️✔️:     covMat[1][1] += w * (ppt->y) * (rpt->y);
    // RDKit✔️✔️:     covMat[1][2] += w * (ppt->y) * (rpt->z);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     covMat[2][0] += w * (ppt->z) * (rpt->x);
    // RDKit✔️✔️:     covMat[2][1] += w * (ppt->z) * (rpt->y);
    // RDKit✔️✔️:     covMat[2][2] += w * (ppt->z) * (rpt->z);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

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
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: void _covertCovMatToQuad(const double covMat[3][3],
    // RDKit✔️✔️:                          const RDGeom::Point3D &rptSum,
    // RDKit✔️✔️:                          const RDGeom::Point3D &pptSum, double wtsSum,
    // RDKit✔️✔️:                          double quad[4][4]) {
    // RDKit✔️✔️:   double PxRx, PxRy, PxRz;
    // RDKit✔️✔️:   double PyRx, PyRy, PyRz;
    // RDKit✔️✔️:   double PzRx, PzRy, PzRz;
    // RDKit✔️✔️:   double temp;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   temp = pptSum.x / wtsSum;
    // RDKit✔️✔️:   PxRx = covMat[0][0] - temp * rptSum.x;
    // RDKit✔️✔️:   PxRy = covMat[0][1] - temp * rptSum.y;
    // RDKit✔️✔️:   PxRz = covMat[0][2] - temp * rptSum.z;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   temp = pptSum.y / wtsSum;
    // RDKit✔️✔️:   PyRx = covMat[1][0] - temp * rptSum.x;
    // RDKit✔️✔️:   PyRy = covMat[1][1] - temp * rptSum.y;
    // RDKit✔️✔️:   PyRz = covMat[1][2] - temp * rptSum.z;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   temp = pptSum.z / wtsSum;
    // RDKit✔️✔️:   PzRx = covMat[2][0] - temp * rptSum.x;
    // RDKit✔️✔️:   PzRy = covMat[2][1] - temp * rptSum.y;
    // RDKit✔️✔️:   PzRz = covMat[2][2] - temp * rptSum.z;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   quad[0][0] = -2.0 * (PxRx + PyRy + PzRz);
    // RDKit✔️✔️:   quad[1][1] = -2.0 * (PxRx - PyRy - PzRz);
    // RDKit✔️✔️:   quad[2][2] = -2.0 * (PyRy - PzRz - PxRx);
    // RDKit✔️✔️:   quad[3][3] = -2.0 * (PzRz - PxRx - PyRy);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   quad[0][1] = quad[1][0] = 2.0 * (PyRz - PzRy);
    // RDKit✔️✔️:   quad[0][2] = quad[2][0] = 2.0 * (PzRx - PxRz);
    // RDKit✔️✔️:   quad[0][3] = quad[3][0] = 2.0 * (PxRy - PyRx);
    // RDKit✔️✔️:   quad[1][2] = quad[2][1] = -2.0 * (PxRy + PyRx);
    // RDKit✔️✔️:   quad[1][3] = quad[3][1] = -2.0 * (PzRx + PxRz);
    // RDKit✔️✔️:   quad[2][3] = quad[3][2] = -2.0 * (PyRz + PzRy);
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

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
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: unsigned int jacobi(double quad[4][4], double eigenVals[4],
    // RDKit✔️✔️:                     double eigenVecs[4][4], unsigned int maxIter) {
    // RDKit✔️✔️:   double offDiagNorm, diagNorm;
    // RDKit✔️✔️:   double b, dma, q, t, c, s;
    // RDKit✔️✔️:   double atemp, vtemp, dtemp;
    // RDKit✔️✔️:   int i, j, k;
    // RDKit✔️✔️:   unsigned int l;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // initialize the eigen vector to Identity
    // RDKit✔️✔️:   for (j = 0; j <= 3; j++) {
    // RDKit✔️✔️:     for (i = 0; i <= 3; i++) {
    // RDKit✔️✔️:       eigenVecs[i][j] = 0.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     eigenVecs[j][j] = 1.0;
    // RDKit✔️✔️:     eigenVals[j] = quad[j][j];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (l = 0; l < maxIter; l++) {
    // RDKit✔️✔️:     diagNorm = 0.0;
    // RDKit✔️✔️:     offDiagNorm = 0.0;
    // RDKit✔️✔️:     for (j = 0; j <= 3; j++) {
    // RDKit✔️✔️:       diagNorm += fabs(eigenVals[j]);
    // RDKit✔️✔️:       for (i = 0; i <= j - 1; i++) {
    // RDKit✔️✔️:         offDiagNorm += fabs(quad[i][j]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (fabs(diagNorm) > 1.e-16 && (offDiagNorm / diagNorm) <= TOLERANCE) {
    // RDKit✔️✔️:       goto Exit_now;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (j = 1; j <= 3; j++) {
    // RDKit✔️✔️:       for (i = 0; i <= j - 1; i++) {
    // RDKit✔️✔️:         b = quad[i][j];
    // RDKit✔️✔️:         if (fabs(b) > 0.0) {
    // RDKit✔️✔️:           dma = eigenVals[j] - eigenVals[i];
    // RDKit✔️✔️:           if ((fabs(dma) + fabs(b)) <= fabs(dma)) {
    // RDKit✔️✔️:             t = b / dma;
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             q = 0.5 * dma / b;
    // RDKit✔️✔️:             t = 1.0 / (fabs(q) + sqrt(1.0 + q * q));
    // RDKit✔️✔️:             if (q < 0.0) {
    // RDKit✔️✔️:               t = -t;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           c = 1.0 / sqrt(t * t + 1.0);
    // RDKit✔️✔️:           s = t * c;
    // RDKit✔️✔️:           quad[i][j] = 0.0;
    // RDKit✔️✔️:           for (k = 0; k <= i - 1; k++) {
    // RDKit✔️✔️:             atemp = c * quad[k][i] - s * quad[k][j];
    // RDKit✔️✔️:             quad[k][j] = s * quad[k][i] + c * quad[k][j];
    // RDKit✔️✔️:             quad[k][i] = atemp;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (k = i + 1; k <= j - 1; k++) {
    // RDKit✔️✔️:             atemp = c * quad[i][k] - s * quad[k][j];
    // RDKit✔️✔️:             quad[k][j] = s * quad[i][k] + c * quad[k][j];
    // RDKit✔️✔️:             quad[i][k] = atemp;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (k = j + 1; k <= 3; k++) {
    // RDKit✔️✔️:             atemp = c * quad[i][k] - s * quad[j][k];
    // RDKit✔️✔️:             quad[j][k] = s * quad[i][k] + c * quad[j][k];
    // RDKit✔️✔️:             quad[i][k] = atemp;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (k = 0; k <= 3; k++) {
    // RDKit✔️✔️:             vtemp = c * eigenVecs[k][i] - s * eigenVecs[k][j];
    // RDKit✔️✔️:             eigenVecs[k][j] = s * eigenVecs[k][i] + c * eigenVecs[k][j];
    // RDKit✔️✔️:             eigenVecs[k][i] = vtemp;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           dtemp = c * c * eigenVals[i] + s * s * eigenVals[j] - 2.0 * c * s * b;
    // RDKit✔️✔️:           eigenVals[j] =
    // RDKit✔️✔️:               s * s * eigenVals[i] + c * c * eigenVals[j] + 2.0 * c * s * b;
    // RDKit✔️✔️:           eigenVals[i] = dtemp;
    // RDKit✔️✔️:         } /* end if */
    // RDKit✔️✔️:       } /* end for i */
    // RDKit✔️✔️:     } /* end for j */
    // RDKit✔️✔️:   } /* end for l */
    // RDKit✔️✔️:
    // RDKit✔️✔️: Exit_now:
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (j = 0; j <= 2; j++) {
    // RDKit✔️✔️:     k = j;
    // RDKit✔️✔️:     dtemp = eigenVals[k];
    // RDKit✔️✔️:     for (i = j + 1; i <= 3; i++) {
    // RDKit✔️✔️:       if (eigenVals[i] < dtemp) {
    // RDKit✔️✔️:         k = i;
    // RDKit✔️✔️:         dtemp = eigenVals[k];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (k > j) {
    // RDKit✔️✔️:       eigenVals[k] = eigenVals[j];
    // RDKit✔️✔️:       eigenVals[j] = dtemp;
    // RDKit✔️✔️:       for (i = 0; i <= 3; i++) {
    // RDKit✔️✔️:         dtemp = eigenVecs[i][k];
    // RDKit✔️✔️:         eigenVecs[i][k] = eigenVecs[i][j];
    // RDKit✔️✔️:         eigenVecs[i][j] = dtemp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return l + 1;
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

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
                // Negate the source positive comparison itself: NaN does not
                // enter `if (fabs(b) > 0.0)` and must leave eigenvectors intact.
                if !(b.abs() > 0.0) {
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
pub fn align_points(
    ref_points: &[[f64; 3]],
    probe_points: &[[f64; 3]],
    weights: Option<&[f64]>,
    reflect_input: bool,
    max_iterations: usize,
) -> Result<(f64, Transform3D), &'static str> {
    // Fixed RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, complete source anchor.
    // RDKit✔️✔️: double AlignPoints(const RDGeom::Point3DConstPtrVect &refPoints,
    // RDKit✔️✔️:                    const RDGeom::Point3DConstPtrVect &probePoints,
    // RDKit✔️✔️:                    RDGeom::Transform3D &trans, const DoubleVector *weights,
    // RDKit✔️✔️:                    bool reflect, unsigned int maxIterations) {
    // RDKit✔️✔️:   unsigned int npt = refPoints.size();
    // RDKit✔️✔️:   PRECONDITION(npt == probePoints.size(), "Mismatch in number of points");
    // RDKit✔️✔️:   trans.setToIdentity();
    // RDKit✔️✔️:   double wtsSum = 0.0;
    // RDKit✔️✔️:   if (weights) {
    // RDKit✔️✔️:     PRECONDITION(npt == weights->size(), "Mismatch in number of points");
    // RDKit✔️✔️:     wtsSum = _sumOfWeights(*weights);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     wtsSum = static_cast<double>(npt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D rptSum = _weightedSumOfPoints(refPoints, weights);
    // RDKit✔️✔️:   RDGeom::Point3D pptSum = _weightedSumOfPoints(probePoints, weights);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double rptSumLenSq = _weightedSumOfLenSq(refPoints, weights);
    // RDKit✔️✔️:   double pptSumLenSq = _weightedSumOfLenSq(probePoints, weights);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double covMat[3][3];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // compute the co-variance matrix
    // RDKit✔️✔️:   _computeCovarianceMat(refPoints, probePoints, weights, covMat);
    // RDKit✔️✔️:   if (reflect) {
    // RDKit✔️✔️:     rptSum *= -1.0;
    // RDKit✔️✔️:     reflectCovMat(covMat);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // convert the covariance matrix to a 4x4 matrix that needs to be diagonalized
    // RDKit✔️✔️:   double quad[4][4];
    // RDKit✔️✔️:   _covertCovMatToQuad(covMat, rptSum, pptSum, wtsSum, quad);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // get the eigenVecs and eigenVals for the matrix
    // RDKit✔️✔️:   double eigenVecs[4][4], eigenVals[4];
    // RDKit✔️✔️:   jacobi(quad, eigenVals, eigenVecs, maxIterations);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // get the quaternion
    // RDKit✔️✔️:   double quater[4];
    // RDKit✔️✔️:   quater[0] = eigenVecs[0][0];
    // RDKit✔️✔️:   quater[1] = eigenVecs[1][0];
    // RDKit✔️✔️:   quater[2] = eigenVecs[2][0];
    // RDKit✔️✔️:   quater[3] = eigenVecs[3][0];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   trans.SetRotationFromQuaternion(quater);
    // RDKit✔️✔️:   if (reflect) {
    // RDKit✔️✔️:     // put the flip in the rotation matrix
    // RDKit✔️✔️:     trans.Reflect();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // compute the SSR value
    // RDKit✔️✔️:   double ssr = eigenVals[0] - (pptSum.lengthSq() + rptSum.lengthSq()) / wtsSum +
    // RDKit✔️✔️:                rptSumLenSq + pptSumLenSq;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((ssr < 0.0) && (fabs(ssr) < TOLERANCE)) {
    // RDKit✔️✔️:     ssr = 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (reflect) {
    // RDKit✔️✔️:     rptSum *= -1.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // set the translation
    // RDKit✔️✔️:   trans.TransformPoint(pptSum);
    // RDKit✔️✔️:   RDGeom::Point3D move = rptSum;
    // RDKit✔️✔️:   move -= pptSum;
    // RDKit✔️✔️:   move /= wtsSum;
    // RDKit✔️✔️:   trans.SetTranslation(move);
    // RDKit✔️✔️:   return ssr;
    // RDKit✔️✔️: }
    // Behavior verified by original native profiles and kernel/source-condition regressions; fixed-size buffers and source loop complexity retained.

    // Verbatim source helpers implemented by the weight and reflection loops below.
    // RDKit✔️✔️: double _sumOfWeights(const DoubleVector &weights) {
    // RDKit✔️✔️:   const double *wData = weights.getData();
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < weights.size(); i++) {
    // RDKit✔️✔️:     CHECK_INVARIANT(wData[i] > 0.0, "Negative weight specified for a point");
    // RDKit✔️✔️:     res += wData[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void reflectCovMat(double covMat[3][3]) {
    // RDKit✔️✔️:   unsigned int i, j;
    // RDKit✔️✔️:   for (i = 0; i < 3; i++) {
    // RDKit✔️✔️:     for (j = 0; j < 3; j++) {
    // RDKit✔️✔️:       covMat[i][j] = -covMat[i][j];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
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
    set_rotation_from_quaternion(
        &mut trans,
        [evecs[0][0], evecs[1][0], evecs[2][0], evecs[3][0]],
    );
    if reflect_input {
        reflect(&mut trans);
    }
    let mut ssr = evals[0]
        - (psum.iter().map(|x| x * x).sum::<f64>() + rsum.iter().map(|x| x * x).sum::<f64>())
            / wsum
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
    let moved = alignment_transform_point(&trans, psum);
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
    fn positive_infinite_weight_keeps_source_nan_rotation_branch() {
        // Pinned AlignPoints.cpp:195 enters rotation only for fabs(b)>0.
        // The source-built native probe accepts positive infinity and returns
        // NaN SSR/translation with the identity rotation, rather than rejecting it.
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let (ssr, matrix) =
            align_points(&points, &points, Some(&[f64::INFINITY, 1.0]), false, 50).unwrap();
        assert!(ssr.is_nan());
        for row in 0..3 {
            for column in 0..3 {
                assert_eq!(matrix[row][column], if row == column { 1.0 } else { 0.0 });
            }
            assert!(matrix[row][3].is_nan());
        }
        assert_eq!(matrix[3], [0.0, 0.0, 0.0, 1.0]);
    }

    #[test]
    fn rigid_translation_has_zero_ssr() {
        let reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let probe = [[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]];
        let (ssr, transform) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!(ssr.abs() < 1.0e-10, "ssr={ssr}");
        for (expected, point) in reference.iter().zip(probe) {
            let actual = alignment_transform_point(&transform, point);
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
            let actual = alignment_transform_point(&transform, point);
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
        let (weighted, _) =
            align_points(&reference, &probe, Some(&[1.0, 1.0, 2.0]), false, 50).unwrap();
        assert!((weighted - 3.75).abs() < 1.0e-4);
        let (weighted, _) =
            align_points(&reference, &probe, Some(&[1.0, 2.0, 2.0]), false, 50).unwrap();
        assert!((weighted - 4.8).abs() < 1.0e-4);
    }

    #[test]
    fn reflection_matches_rdkit_source_regression() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let probe = [
            [2.0, 2.0, 3.0],
            [3.0, 2.0, 3.0],
            [2.0, 2.0, 4.0],
            [2.0, 3.0, 3.0],
        ];
        let (without_reflection, _) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!((without_reflection - 1.0).abs() < 1.0e-4);
        let (with_reflection, transform) =
            align_points(&reference, &probe, None, true, 50).unwrap();
        assert!(with_reflection.abs() < 1.0e-4);
        for (expected, point) in reference.iter().zip(probe) {
            let actual = alignment_transform_point(&transform, point);
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
