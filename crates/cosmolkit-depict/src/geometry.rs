//! Source-shaped private primitives used by the RDKit 2D layout owner.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use cosmolkit_model::{PropertyValue, TopologyBlock};

pub(crate) const BOND_LEN: f64 = 1.5;
pub(crate) type Point2 = [f64; 2];
pub(crate) type PointMap = BTreeMap<usize, Point2>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum GeometryError {
    AtomIndexOutOfRange { atom: usize, atom_count: usize },
    AtomCountTooLarge { atom_count: usize },
    InvalidRankProperty { atom: usize, key: &'static str },
    NotEnoughNeighbors { atom: usize, count: usize },
}

fn unsigned_rank_property(value: &PropertyValue) -> Option<u32> {
    match value {
        PropertyValue::Int(value) => u32::try_from(*value).ok(),
        PropertyValue::String(value) => value.parse().ok(),
        PropertyValue::Double(_) | PropertyValue::Bool(_) => None,
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Transform2D {
    data: [f64; 6],
}

impl Transform2D {
    pub(crate) const fn identity() -> Self {
        Self {
            data: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        }
    }

    pub(crate) const fn from_linear_rows(first: Point2, second: Point2) -> Self {
        Self {
            data: [first[0], first[1], 0.0, second[0], second[1], 0.0],
        }
    }

    pub(crate) fn around(point: Point2, angle: f64) -> Self {
        // RDKit✔️🔝: const unsigned int DIM_2D = 3;
        // RDKit✔️🔝:
        // RDKit✔️🔝:   //! \brief Constructor
        // RDKit✔️🔝:   /*!
        // RDKit✔️🔝:     Initialize to an identity matrix transformation
        // RDKit✔️🔝:     This is a 3x3 matrix that includes the rotation and translation parts
        // RDKit✔️🔝:     see Foley's "Introduction to Computer Graphics" for the representation
        // RDKit✔️🔝:
        // RDKit✔️🔝:     Operator *= and = are provided by the parent class square matrix.
        // RDKit✔️🔝:     Operator *= needs some explanation, since the order matters. This transform
        // RDKit✔️🔝:     gets set to
        // RDKit✔️🔝:     the combination other and the current state of this transform
        // RDKit✔️🔝:     If this_old and this_new are the states of this object before and after this
        // RDKit✔️🔝:     function
        // RDKit✔️🔝:     we have
        // RDKit✔️🔝:             this_new(point) = this_old(other(point))
        // RDKit✔️🔝:   */
        // RDKit✔️🔝:   Transform2D() : RDNumeric::SquareMatrix<double>(DIM_2D, 0.0) {
        // RDKit✔️🔝:     for (unsigned int i = 0; i < DIM_2D; i++) {
        // RDKit✔️🔝:       unsigned int id = i * (DIM_2D + 1);
        // RDKit✔️🔝:       d_data[id] = 1.0;
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:
        // RDKit✔️🔝: void Transform2D::setToIdentity() {
        // RDKit✔️🔝:   double *data = d_data.get();
        // RDKit✔️🔝:   memset(static_cast<void *>(data), 0, d_dataSize * sizeof(double));
        // RDKit✔️🔝:   for (unsigned int i = 0; i < DIM_2D; i++) {
        // RDKit✔️🔝:     unsigned int id = i * (DIM_2D + 1);
        // RDKit✔️🔝:     data[id] = 1.0;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // RDKit✔️🔝:
        // RDKit✔️🔝: void Transform2D::TransformPoint(Point2D &pt) const {
        // RDKit✔️🔝:   double *data = d_data.get();
        // RDKit✔️🔝:   double x = data[0] * pt.x + data[1] * pt.y + data[2];
        // RDKit✔️🔝:   double y = data[3] * pt.x + data[4] * pt.y + data[5];
        // RDKit✔️🔝:
        // RDKit✔️🔝:   pt.x = x;
        // RDKit✔️🔝:   pt.y = y;
        // RDKit✔️🔝: }
        // RDKit✔️🔝:
        // RDKit✔️🔝: void Transform2D::SetTranslation(const Point2D &pt) {
        // RDKit✔️🔝:   unsigned int i = DIM_2D - 1;
        // RDKit✔️🔝:   double *data = d_data.get();
        // RDKit✔️🔝:   data[i] = pt.x;
        // RDKit✔️🔝:   i += DIM_2D;
        // RDKit✔️🔝:   data[i] = pt.y;
        // RDKit✔️🔝:   i += DIM_2D;
        // RDKit✔️🔝:   data[i] = 1.0;
        // RDKit✔️🔝: }
        // RDKit✔️🔝:
        // RDKit✔️🔝: void Transform2D::SetTransform(const Point2D &pt, double angle) {
        // RDKit✔️🔝:   this->setToIdentity();
        // RDKit✔️🔝:
        // RDKit✔️🔝:   Transform2D trans1;
        // RDKit✔️🔝:   trans1.SetTranslation(-pt);
        // RDKit✔️🔝:   double *data = d_data.get();
        // RDKit✔️🔝:   // set the rotation
        // RDKit✔️🔝:   data[0] = cos(angle);
        // RDKit✔️🔝:   data[1] = -sin(angle);
        // RDKit✔️🔝:   data[3] = sin(angle);
        // RDKit✔️🔝:   data[4] = cos(angle);
        // RDKit✔️🔝:
        // RDKit✔️🔝:   (*this) *= trans1;
        // RDKit✔️🔝:
        // RDKit✔️🔝:   // translation back to the original coordinate
        // RDKit✔️🔝:   Transform2D trans2;
        // RDKit✔️🔝:   trans2.SetTranslation(pt);
        // RDKit✔️🔝:   trans2 *= (*this);
        // RDKit✔️🔝:
        // RDKit✔️🔝:   // now combine them
        // RDKit✔️🔝:   this->assign(trans2);
        // RDKit✔️🔝: }
        // RDKit✔️🔝:
        // RDKit✔️🔝:   virtual SquareMatrix<TYPE> &operator*=(const SquareMatrix<TYPE> &B) {
        // RDKit✔️🔝:     CHECK_INVARIANT(this->d_nCols == B.numRows(),
        // RDKit✔️🔝:                     "Size mismatch during multiplication");
        // RDKit✔️🔝:
        // RDKit✔️🔝:     const TYPE *bData = B.getData();
        // RDKit✔️🔝:     TYPE *newData = new TYPE[this->d_dataSize];
        // RDKit✔️🔝:     unsigned int i, j, k;
        // RDKit✔️🔝:     unsigned int idA, idAt, idC, idCt, idB;
        // RDKit✔️🔝:     TYPE *data = this->d_data.get();
        // RDKit✔️🔝:     for (i = 0; i < this->d_nRows; i++) {
        // RDKit✔️🔝:       idA = i * this->d_nRows;
        // RDKit✔️🔝:       idC = idA;
        // RDKit✔️🔝:       for (j = 0; j < this->d_nCols; j++) {
        // RDKit✔️🔝:         idCt = idC + j;
        // RDKit✔️🔝:         newData[idCt] = (TYPE)(0.0);
        // RDKit✔️🔝:         for (k = 0; k < this->d_nCols; k++) {
        // RDKit✔️🔝:           idAt = idA + k;
        // RDKit✔️🔝:           idB = k * this->d_nRows + j;
        // RDKit✔️🔝:           newData[idCt] += (data[idAt] * bData[idB]);
        // RDKit✔️🔝:         }
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     boost::shared_array<TYPE> tsptr(newData);
        // RDKit✔️🔝:     this->d_data.swap(tsptr);
        // RDKit✔️🔝:     return (*this);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:
        // RDKit✔️🔝:   //! Copy operator.
        // RDKit✔️🔝:   /*! We make a copy of the other Matrix's data.
        // RDKit✔️🔝:    */
        // RDKit✔️🔝:
        // RDKit✔️🔝:   Matrix<TYPE> &assign(const Matrix<TYPE> &other) {
        // RDKit✔️🔝:     PRECONDITION(d_nRows == other.numRows(),
        // RDKit✔️🔝:                  "Num rows mismatch in matrix copying");
        // RDKit✔️🔝:     PRECONDITION(d_nCols == other.numCols(),
        // RDKit✔️🔝:                  "Num cols mismatch in matrix copying");
        // RDKit✔️🔝:     const TYPE *otherData = other.getData();
        // RDKit✔️🔝:     TYPE *data = d_data.get();
        // RDKit✔️🔝:     memcpy(static_cast<void *>(data), static_cast<const void *>(otherData),
        // RDKit✔️🔝:            d_dataSize * sizeof(TYPE));
        // RDKit✔️🔝:     return *this;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:
        // RDKit✔️🔝:   constexpr Point2D operator-() const {
        // RDKit✔️🔝:     Point2D res(x, y);
        // RDKit✔️🔝:     res.x *= -1.0;
        // RDKit✔️🔝:     res.y *= -1.0;
        // RDKit✔️🔝:     return res;
        // RDKit✔️🔝:   }

        // RDKit behavior review: source-order output bits pass all twelve
        // frozen finite inputs; NaN and overflow cases remain unmodeled.
        // RDKit complexity review: both fixed 3x3 products retain 27 multiply
        // accumulations each in the same i/j/k order. Fixed stack output
        // arrays replace the two per-product new[]/shared-array buffers
        // without changing arithmetic or result order, an expected local
        // helper cost improvement only, not a whole-layout claim.
        let identity = || {
            let mut matrix = [0.0; 9];
            for i in 0..3 {
                let id = i * (3 + 1);
                matrix[id] = 1.0;
            }
            matrix
        };
        let multiply = |data: &[f64; 9], b_data: &[f64; 9]| {
            let mut new_data = [0.0; 9];
            for i in 0..3 {
                let id_a = i * 3;
                let id_c = id_a;
                for j in 0..3 {
                    let id_ct = id_c + j;
                    new_data[id_ct] = 0.0;
                    for k in 0..3 {
                        let id_at = id_a + k;
                        let id_b = k * 3 + j;
                        new_data[id_ct] += data[id_at] * b_data[id_b];
                    }
                }
            }
            new_data
        };

        let mut data = identity();
        let mut trans1 = identity();
        let mut negated_point = point;
        negated_point[0] *= -1.0;
        negated_point[1] *= -1.0;
        let mut i = 3 - 1;
        trans1[i] = negated_point[0];
        i += 3;
        trans1[i] = negated_point[1];
        i += 3;
        trans1[i] = 1.0;

        data[0] = angle.cos();
        data[1] = -angle.sin();
        data[3] = angle.sin();
        data[4] = angle.cos();

        data = multiply(&data, &trans1);

        let mut trans2 = identity();
        let mut i = 3 - 1;
        trans2[i] = point[0];
        i += 3;
        trans2[i] = point[1];
        i += 3;
        trans2[i] = 1.0;
        trans2 = multiply(&trans2, &data);

        Self {
            data: [
                trans2[0], trans2[1], trans2[2], trans2[3], trans2[4], trans2[5],
            ],
        }
    }

    pub(crate) fn transform_point(self, point: Point2) -> Point2 {
        // RDKit❗✔️: void Transform2D::TransformPoint(Point2D &pt) const {
        // RDKit❗✔️:   double *data = d_data.get();
        // RDKit❗✔️:   double x = data[0] * pt.x + data[1] * pt.y + data[2];
        // RDKit❗✔️:   double y = data[3] * pt.x + data[4] * pt.y + data[5];
        // RDKit❗✔️:
        // RDKit❗✔️:   pt.x = x;
        // RDKit❗✔️:   pt.y = y;
        // RDKit❗✔️: }
        let data = self.data;
        let x = data[0] * point[0] + data[1] * point[1] + data[2];
        let y = data[3] * point[0] + data[4] * point[1] + data[5];
        [x, y]
    }

    pub(crate) fn from_point_pairs(ref1: Point2, ref2: Point2, pt1: Point2, pt2: Point2) -> Self {
        // RDKit behavior review: retain each source branch and each binary64
        // operation in source order; exact corpus closure remains pending.
        // RDKit complexity review: fixed scalar arithmetic and stack arrays
        // only, O(1) time and space, no heap allocation or collection scan.
        // RDKit❗✔️: void Transform2D::setToIdentity() {
        // RDKit❗✔️:   double *data = d_data.get();
        // RDKit❗✔️:   memset(static_cast<void *>(data), 0, d_dataSize * sizeof(double));
        // RDKit❗✔️:   for (unsigned int i = 0; i < DIM_2D; i++) {
        // RDKit❗✔️:     unsigned int id = i * (DIM_2D + 1);
        // RDKit❗✔️:     data[id] = 1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: void Transform2D::TransformPoint(Point2D &pt) const {
        // RDKit❗✔️:   double *data = d_data.get();
        // RDKit❗✔️:   double x = data[0] * pt.x + data[1] * pt.y + data[2];
        // RDKit❗✔️:   double y = data[3] * pt.x + data[4] * pt.y + data[5];
        // RDKit❗✔️:
        // RDKit❗✔️:   pt.x = x;
        // RDKit❗✔️:   pt.y = y;
        // RDKit❗✔️: }
        // RDKit❗✔️: double length() const override {
        // RDKit❗✔️:   // double res = pow(x,2) + pow(y,2);
        // RDKit❗✔️:   double res = x * x + y * y;
        // RDKit❗✔️:   return sqrt(res);
        // RDKit❗✔️: }
        // RDKit❗✔️: constexpr double dotProduct(const Point2D &other) const {
        // RDKit❗✔️:   double res = x * (other.x) + y * (other.y);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // RDKit❗✔️: Point2D operator-(const Point2D &p1, const Point2D &p2) {
        // RDKit❗✔️:   Point2D res;
        // RDKit❗✔️:   res.x = p1.x - p2.x;
        // RDKit❗✔️:   res.y = p1.y - p2.y;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // RDKit❗✔️: void Transform2D::SetTransform(const Point2D &ref1, const Point2D &ref2,
        // RDKit❗✔️:                                const Point2D &pt1, const Point2D &pt2) {
        // RDKit❗✔️:   // compute the angle between the two vectors
        // RDKit❗✔️:   Point2D rvec = ref2 - ref1;
        // RDKit❗✔️:   Point2D pvec = pt2 - pt1;
        // RDKit❗✔️:
        // RDKit❗✔️:   double dp = rvec.dotProduct(pvec);
        // RDKit❗✔️:   double lp = (rvec.length()) * (pvec.length());
        // RDKit❗✔️:   if (lp <= 0.0) {
        // RDKit❗✔️:     this->setToIdentity();
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   double cval = dp / lp;
        // RDKit❗✔️:   if (cval < -1.0) {
        // RDKit❗✔️:     cval = -1.0;
        // RDKit❗✔️:   } else if (cval > 1.0) {
        // RDKit❗✔️:     cval = 1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   double ang = acos(cval);
        // RDKit❗✔️:
        // RDKit❗✔️:   // figure out if we have to do clock wise or anti clock wise rotation
        // RDKit❗✔️:   double cross = (pvec.x) * (rvec.y) - (pvec.y) * (rvec.x);
        // RDKit❗✔️:
        // RDKit❗✔️:   if (cross < 0.0) {
        // RDKit❗✔️:     ang *= -1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   this->setToIdentity();
        // RDKit❗✔️:   // set the rotation
        // RDKit❗✔️:   double *data = d_data.get();
        // RDKit❗✔️:   data[0] = cos(ang);
        // RDKit❗✔️:   data[1] = -sin(ang);
        // RDKit❗✔️:   data[3] = sin(ang);
        // RDKit❗✔️:   data[4] = cos(ang);
        // RDKit❗✔️:
        // RDKit❗✔️:   // apply this rotation to pt1 and compute the translation
        // RDKit❗✔️:   Point2D npt1 = pt1;
        // RDKit❗✔️:   this->TransformPoint(npt1);
        // RDKit❗✔️:   data[DIM_2D - 1] = ref1.x - npt1.x;
        // RDKit❗✔️:   data[2 * DIM_2D - 1] = ref1.y - npt1.y;
        // RDKit❗✔️: }
        let rvec = [ref2[0] - ref1[0], ref2[1] - ref1[1]];
        let pvec = [pt2[0] - pt1[0], pt2[1] - pt1[1]];
        let dp = rvec[0] * pvec[0] + rvec[1] * pvec[1];
        let rvec_length = (rvec[0] * rvec[0] + rvec[1] * rvec[1]).sqrt();
        let pvec_length = (pvec[0] * pvec[0] + pvec[1] * pvec[1]).sqrt();
        let lp = rvec_length * pvec_length;
        if lp <= 0.0 {
            return Self::identity();
        }
        let mut cval = dp / lp;
        if cval < -1.0 {
            cval = -1.0;
        } else if cval > 1.0 {
            cval = 1.0;
        }
        let mut ang = cval.acos();
        let cross = pvec[0] * rvec[1] - pvec[1] * rvec[0];
        if cross < 0.0 {
            ang *= -1.0;
        }
        let mut transform = Self::identity();
        transform.data[0] = ang.cos();
        transform.data[1] = -ang.sin();
        transform.data[3] = ang.sin();
        transform.data[4] = ang.cos();
        let npt1 = transform.transform_point(pt1);
        transform.data[2] = ref1[0] - npt1[0];
        transform.data[5] = ref1[1] - npt1[1];
        transform
    }
}

pub(crate) fn embed_ring(ring: &[usize]) -> PointMap {
    // RDKit❗✔️: RDGeom::INT_POINT2D_MAP embedRing(const RDKit::INT_VECT &ring) {
    // RDKit❗✔️:   unsigned int na = ring.size();
    // RDKit❗✔️:   double ang = 2 * M_PI / na;
    // RDKit❗✔️:   double al = BOND_LEN / (sqrt(2 * (1 - cos(ang))));
    // RDKit❗✔️:   RDGeom::INT_POINT2D_MAP res;
    // RDKit❗✔️:   for (unsigned int i = 0; i < na; ++i) {
    // RDKit❗✔️:     auto x = al * cos(i * ang);
    // RDKit❗✔️:     auto y = al * sin(i * ang);
    // RDKit❗✔️:     RDGeom::Point2D loc(x, y);
    // RDKit❗✔️:     res[ring[i]] = loc;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    let na = ring.len();
    let ang = 2.0 * PI / na as f64;
    let al = BOND_LEN / (2.0 * (1.0 - ang.cos())).sqrt();
    let mut result = PointMap::new();
    for (i, atom) in ring.iter().copied().enumerate() {
        let x = al * (i as f64 * ang).cos();
        let y = al * (i as f64 * ang).sin();
        result.insert(atom, [x, y]);
    }
    result
}

pub(crate) fn transform_points(points: &mut PointMap, transform: Transform2D) {
    // RDKit❗✔️: void transformPoints(RDGeom::INT_POINT2D_MAP &nringCor,
    // RDKit❗✔️:                      const RDGeom::Transform2D &trans) {
    // RDKit❗✔️:   std::for_each(nringCor.begin(), nringCor.end(),
    // RDKit❗✔️:                 [&trans](auto &elem) { trans.TransformPoint(elem.second); });
    // RDKit❗✔️: }
    for point in points.values_mut() {
        *point = transform.transform_point(*point);
    }
}

pub(crate) fn compute_bisect_point(center: Point2, angle: f64, a: Point2, b: Point2) -> Point2 {
    // RDKit❗✔️: RDGeom::Point2D computeBisectPoint(const RDGeom::Point2D &rcr, double ang,
    // RDKit❗✔️:                                    const RDGeom::Point2D &nb1,
    // RDKit❗✔️:                                    const RDGeom::Point2D &nb2) {
    // RDKit❗✔️:   RDGeom::Point2D cloc = nb1;
    // RDKit❗✔️:   cloc += nb2;
    // RDKit❗✔️:   cloc *= 0.5;
    // RDKit❗✔️:   if (ang > M_PI) {
    // RDKit❗✔️:     cloc -= rcr;
    // RDKit❗✔️:     cloc *= -1.0;
    // RDKit❗✔️:     cloc += rcr;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return cloc;
    // RDKit❗✔️: }
    let mut loc = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    if angle > PI {
        loc = [2.0 * center[0] - loc[0], 2.0 * center[1] - loc[1]];
    }
    loc
}

pub(crate) fn reflect_point(point: Point2, a: Point2, b: Point2) -> Point2 {
    // RDKit❗✔️: RDGeom::Point2D reflectPoint(const RDGeom::Point2D &point,
    // RDKit❗✔️:                              const RDGeom::Point2D &loc1,
    // RDKit❗✔️:                              const RDGeom::Point2D &loc2) {
    // RDKit❗✔️:   RDGeom::Point2D org(0.0, 0.0);
    // RDKit❗✔️:   RDGeom::Point2D xaxis(1.0, 0.0);
    // RDKit❗✔️:   RDGeom::Point2D cent = (loc1 + loc2);
    // RDKit❗✔️:   cent *= 0.5;
    // RDKit❗✔️:   RDGeom::Transform2D trans;
    // RDKit❗✔️:   trans.SetTransform(org, xaxis, cent, loc1);
    // RDKit❗✔️:   RDGeom::Transform2D itrans;
    // RDKit❗✔️:   itrans.SetTransform(cent, loc1, org, xaxis);
    // RDKit❗✔️:   RDGeom::INT_POINT2D_MAP_I nci;
    // RDKit❗✔️:   RDGeom::Point2D res;
    // RDKit❗✔️:   res = point;
    // RDKit❗✔️:   trans.TransformPoint(res);
    // RDKit❗✔️:   res.y = -res.y;
    // RDKit❗✔️:   itrans.TransformPoint(res);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    let origin = [0.0, 0.0];
    let xaxis = [1.0, 0.0];
    let center = [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5];
    let forward = Transform2D::from_point_pairs(origin, xaxis, center, a);
    let inverse = Transform2D::from_point_pairs(center, a, origin, xaxis);
    let mut result = forward.transform_point(point);
    result[1] = -result[1];
    inverse.transform_point(result)
}

pub(crate) fn reflect_points(points: &mut PointMap, a: Point2, b: Point2) {
    // RDKit❗✔️: void reflectPoints(RDGeom::INT_POINT2D_MAP &coordMap,
    // RDKit❗✔️:                    const RDGeom::Point2D &loc1, const RDGeom::Point2D &loc2) {
    // RDKit❗✔️:   std::for_each(coordMap.begin(), coordMap.end(), [&loc1, &loc2](auto &elem) {
    // RDKit❗✔️:     reflectPoint(elem.second, loc1, loc2);
    // RDKit❗✔️:   });
    // RDKit❗✔️: }
    // The source discards the returned point. Preserve that observable no-op.
    for point in points.values() {
        let _ = reflect_point(*point, a, b);
    }
}

pub(crate) fn atom_depict_rank(
    topology: &TopologyBlock,
    atom: usize,
) -> Result<i32, GeometryError> {
    // RDKit❗✔️: inline int getAtomDepictRank(const RDKit::Atom *at) {
    // RDKit❗✔️:   const int maxAtNum = 1000;
    // RDKit❗✔️:   const int maxDeg = 100;
    // RDKit❗✔️:   int anum = at->getAtomicNum();
    // RDKit❗✔️:   anum = anum == 1 ? maxAtNum : anum;
    // RDKit❗✔️:   int deg = at->getDegree();
    // RDKit❗✔️:   return maxDeg * anum + deg;
    // RDKit❗✔️: }
    let Some(value) = topology.atoms.get(atom) else {
        return Err(GeometryError::AtomIndexOutOfRange {
            atom,
            atom_count: topology.atoms.len(),
        });
    };
    let atomic_number = i32::from(value.atomic_number());
    let atomic_number = if atomic_number == 1 {
        1000
    } else {
        atomic_number
    };
    Ok(100 * atomic_number + topology.adjacency.neighbors_of(atom).len() as i32)
}

pub(crate) fn rank_atoms_by_rank(
    topology: &TopologyBlock,
    atoms: &[usize],
    ascending: bool,
) -> Result<Vec<usize>, GeometryError> {
    // RDKit❗✔️: T rankAtomsByRank(const RDKit::ROMol &mol, const T &commAtms, bool ascending) {
    // RDKit❗✔️:   const auto natms = commAtms.size();
    // RDKit❗✔️:   INT_PAIR_VECT rankAid;
    // RDKit❗✔️:   rankAid.reserve(natms);
    // RDKit❗✔️:   for (const auto aid : commAtms) {
    // RDKit❗✔️:     unsigned int rank = aid;
    // RDKit❗✔️:     const auto at = mol.getAtomWithIdx(aid);
    // RDKit❗✔️:     if (!at->getPropIfPresent(RDKit::common_properties::_CIPRank, rank)) {
    // RDKit❗✔️:       if (at->getPropIfPresent(RDKit::common_properties::_ChiralAtomRank,
    // RDKit❗✔️:                                rank)) {
    // RDKit❗✔️:         rank = mol.getNumAtoms() - rank;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       rank += mol.getNumAtoms() * getAtomDepictRank(at);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     rankAid.emplace_back(rank, aid);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (ascending) {
    // RDKit❗✔️:     std::stable_sort(rankAid.begin(), rankAid.end(),
    // RDKit❗✔️:                      [](const auto &e1, const auto &e2) { return e1 < e2; });
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     std::stable_sort(rankAid.begin(), rankAid.end(),
    // RDKit❗✔️:                      [](const auto &e1, const auto &e2) { return e1 > e2; });
    // RDKit❗✔️:   }
    // RDKit❗✔️:   T res;
    // RDKit❗✔️:   std::for_each(rankAid.begin(), rankAid.end(),
    // RDKit❗✔️:                 [&res](const auto &elem) { res.push_back(elem.second); });
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    let count =
        u32::try_from(topology.atoms.len()).map_err(|_| GeometryError::AtomCountTooLarge {
            atom_count: topology.atoms.len(),
        })?;
    let mut ranked = Vec::with_capacity(atoms.len());
    for &atom in atoms {
        let Some(value) = topology.atoms.get(atom) else {
            return Err(GeometryError::AtomIndexOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        };
        let mut rank = u32::try_from(atom).map_err(|_| GeometryError::AtomCountTooLarge {
            atom_count: topology.atoms.len(),
        })?;
        if let Some(text) = value.prop("_CIPRank") {
            rank = unsigned_rank_property(text).ok_or(GeometryError::InvalidRankProperty {
                atom,
                key: "_CIPRank",
            })?;
        } else {
            if let Some(text) = value.prop("_ChiralAtomRank") {
                let chiral_rank =
                    unsigned_rank_property(text).ok_or(GeometryError::InvalidRankProperty {
                        atom,
                        key: "_ChiralAtomRank",
                    })?;
                rank = count.wrapping_sub(chiral_rank);
            }
            rank = rank.wrapping_add(count.wrapping_mul(atom_depict_rank(topology, atom)? as u32));
        }
        // INT_PAIR_VECT stores the computed unsigned rank in a signed `int`.
        ranked.push((rank as i32, atom));
    }
    if ascending {
        ranked.sort_by(|a, b| a.cmp(b));
    } else {
        ranked.sort_by(|a, b| b.cmp(a));
    }
    Ok(ranked.into_iter().map(|(_, atom)| atom).collect())
}

pub(crate) fn set_neighbor_order(
    topology: &TopologyBlock,
    atom: usize,
    neighbors: &[usize],
) -> Result<Vec<usize>, GeometryError> {
    // RDKit❗✔️: RDKit::INT_VECT setNbrOrder(unsigned int aid, const RDKit::INT_VECT &nbrs,
    // RDKit❗✔️:                             const RDKit::ROMol &mol) {
    // RDKit❗✔️:   PRECONDITION(aid < mol.getNumAtoms(), "");
    // RDKit❗✔️:   PR_QUEUE subsAid;
    // RDKit❗✔️:   int ref = -1;
    // RDKit❗✔️:   for (auto anbr : mol.atomNeighbors(mol.getAtomWithIdx(aid))) {
    // RDKit❗✔️:     if (std::find(nbrs.begin(), nbrs.end(), static_cast<int>(anbr->getIdx())) ==
    // RDKit❗✔️:         nbrs.end()) {
    // RDKit❗✔️:       ref = anbr->getIdx();
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDKit::INT_VECT thold = nbrs;
    // RDKit❗✔️:   if (ref >= 0) {
    // RDKit❗✔️:     thold.push_back(ref);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   CHECK_INVARIANT(thold.size() > 3, "");
    // RDKit❗✔️:   thold = rankAtomsByRank(mol, thold);
    // RDKit❗✔️:   unsigned int ln = thold.size();
    // RDKit❗✔️:   int tint = thold[ln - 3];
    // RDKit❗✔️:   thold[ln - 3] = thold[ln - 2];
    // RDKit❗✔️:   thold[ln - 2] = tint;
    // RDKit❗✔️:   RDKit::INT_VECT res;
    // RDKit❗✔️:   res.reserve(thold.size());
    // RDKit❗✔️:   auto pos = std::find(thold.begin(), thold.end(), ref);
    // RDKit❗✔️:   if (pos != thold.end()) {
    // RDKit❗✔️:     res.insert(res.end(), pos + 1, thold.end());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (pos != thold.begin()) {
    // RDKit❗✔️:     res.insert(res.end(), thold.begin(), pos);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   POSTCONDITION(res.size() == nbrs.size(), "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    if atom >= topology.atoms.len() {
        return Err(GeometryError::AtomIndexOutOfRange {
            atom,
            atom_count: topology.atoms.len(),
        });
    }
    let mut reference = None;
    for neighbor in topology.adjacency.neighbors_of(atom) {
        if !neighbors.contains(&neighbor.atom_index) {
            reference = Some(neighbor.atom_index);
        }
    }
    let mut ordered = neighbors.to_vec();
    if let Some(reference) = reference {
        ordered.push(reference);
    }
    if ordered.len() <= 3 {
        return Err(GeometryError::NotEnoughNeighbors {
            atom,
            count: ordered.len(),
        });
    }
    let mut ordered = rank_atoms_by_rank(topology, &ordered, true)?;
    let len = ordered.len();
    ordered.swap(len - 3, len - 2);
    let result =
        if let Some(position) = reference.and_then(|r| ordered.iter().position(|&a| a == r)) {
            ordered[position + 1..]
                .iter()
                .chain(&ordered[..position])
                .copied()
                .collect()
        } else {
            ordered
        };
    Ok(result)
}

#[cfg(test)]
mod accepted_transform_regressions {
    use super::{Point2, Transform2D};

    #[derive(Clone, Copy)]
    struct NativeLiteral {
        pair: [u64; 8],
        point: [u64; 2],
        matrix: [u64; 6],
        output: [u64; 2],
    }

    fn point(bits: [u64; 2]) -> Point2 {
        bits.map(f64::from_bits)
    }

    #[derive(Clone, Copy)]
    struct AroundNativeLiteral {
        center: [u64; 2],
        angle: u64,
        point: [u64; 2],
        matrix: [u64; 6],
        output: [u64; 2],
    }

    #[test]
    fn d2_pair_native_literal_product() {
        // Frozen source-linked native literals from
        // dev/gap_reports/depict_2d/point_pair_transform.md. Every fixture is
        // executed twice; all actual discrepancies are reported together.
        const FIXTURES: [NativeLiteral; 11] = [
            NativeLiteral {
                pair: [
                    13842690664439568919,
                    4614895549384485632,
                    13840743450532881578,
                    4618222693736076515,
                    4600774667239816340,
                    4608144052124252029,
                    13830699860855537819,
                    13828302655841107966,
                ],
                point: [4608425305167921899, 0],
                matrix: [
                    13830554455654793215,
                    4490088828525397190,
                    13842246724002309097,
                    13713460865380172998,
                    13830554455654793215,
                    4616908891964138523,
                ],
                output: [13843683345501127844, 4616908891942731194],
            },
            NativeLiteral {
                pair: [
                    4611686018427387904,
                    4613937818241073152,
                    4611686018427387904,
                    4613937818241073152,
                    13830554455654793216,
                    4617315517961601024,
                    0,
                    4617315517961601024,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [4607182418800017408, 0, 0, 0, 4607182418800017408, 0],
                output: [4598175219545276416, 13826050856027422720],
            },
            NativeLiteral {
                pair: [
                    4611686018427387904,
                    4613937818241073152,
                    4613937818241073152,
                    4613937818241073152,
                    13830554455654793216,
                    4617315517961601024,
                    13830554455654793216,
                    4617315517961601024,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [4607182418800017408, 0, 0, 0, 4607182418800017408, 0],
                output: [4598175219545276416, 13826050856027422720],
            },
            NativeLiteral {
                pair: [
                    4611686018427387904,
                    4613937818241073152,
                    4613937818241073152,
                    4613937818241073152,
                    13830554455654793216,
                    4617315517961601024,
                    0,
                    4617315517961601024,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    4607182418800017408,
                    9223372036854775808,
                    4613937818241073152,
                    0,
                    4607182418800017408,
                    13835058055282163712,
                ],
                output: [4614500768194494464, 13836183955189006336],
            },
            NativeLiteral {
                pair: [0, 0, 0, 4607182418800017408, 0, 0, 4607182418800017408, 0],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    4364452196894661639,
                    13830554455654793216,
                    0,
                    4607182418800017408,
                    4364452196894661639,
                    0,
                ],
                output: [4602678819172646912, 4598175219545276415],
            },
            NativeLiteral {
                pair: [0, 0, 0, 13830554455654793216, 0, 0, 4607182418800017408, 0],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    4364452196894661639,
                    4607182418800017408,
                    0,
                    13830554455654793216,
                    4364452196894661639,
                    0,
                ],
                output: [13826050856027422720, 13821547256400052225],
            },
            NativeLiteral {
                pair: [0, 0, 13830554455654793216, 0, 0, 0, 4607182418800017408, 0],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    13830554455654793216,
                    13592327833376807943,
                    0,
                    4368955796522032135,
                    13830554455654793216,
                    0,
                ],
                output: [13821547256400052222, 4602678819172646912],
            },
            NativeLiteral {
                pair: [
                    0,
                    0,
                    4607182418800017408,
                    4487126258331716666,
                    0,
                    0,
                    4607182418800017408,
                    0,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    4607182418800017408,
                    9223372036854775808,
                    0,
                    0,
                    4607182418800017408,
                    0,
                ],
                output: [4598175219545276416, 13826050856027422720],
            },
            NativeLiteral {
                pair: [
                    0,
                    0,
                    4607182418800017408,
                    0,
                    0,
                    0,
                    1614679632300144556,
                    1614679632300144556,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [4607182418800017408, 0, 0, 0, 4607182418800017408, 0],
                output: [4598175219545276416, 13826050856027422720],
            },
            NativeLiteral {
                pair: [
                    0,
                    0,
                    4602678819172646912,
                    4599075939470750515,
                    0,
                    0,
                    4602678819172646912,
                    4599075939470750515,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    4607182418800017408,
                    9223372036854775808,
                    0,
                    0,
                    4607182418800017408,
                    0,
                ],
                output: [4598175219545276416, 13826050856027422720],
            },
            NativeLiteral {
                pair: [
                    0,
                    0,
                    13826050856027422720,
                    13822447976325526323,
                    0,
                    0,
                    4602678819172646912,
                    4599075939470750515,
                ],
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    13830554455654793216,
                    13592327833376807943,
                    0,
                    4368955796522032135,
                    13830554455654793216,
                    0,
                ],
                output: [13821547256400052222, 4602678819172646912],
            },
        ];

        let mut from_point_pairs_calls = 0;
        let mut transform_point_calls = 0;
        let mut mismatches = Vec::new();

        for repeat in 0..2 {
            for (fixture, expected) in FIXTURES.iter().enumerate() {
                let pair_values = expected.pair.map(f64::from_bits);
                let pair_bits_before = pair_values.map(f64::to_bits);
                let point_value = point(expected.point);
                let point_bits_before = point_value.map(f64::to_bits);
                let [ref1_x, ref1_y, ref2_x, ref2_y, pt1_x, pt1_y, pt2_x, pt2_y] = pair_values;

                let transform = Transform2D::from_point_pairs(
                    [ref1_x, ref1_y],
                    [ref2_x, ref2_y],
                    [pt1_x, pt1_y],
                    [pt2_x, pt2_y],
                );
                from_point_pairs_calls += 1;

                let matrix_bits = transform.data.map(f64::to_bits);
                if matrix_bits != expected.matrix {
                    mismatches.push(format!(
                        "fixture {fixture} repeat {repeat}: matrix expected {:?}, got {:?}",
                        expected.matrix, matrix_bits
                    ));
                }
                if pair_values.map(f64::to_bits) != pair_bits_before {
                    mismatches.push(format!(
                        "fixture {fixture} repeat {repeat}: point-pair inputs changed from {:?} to {:?}",
                        pair_bits_before,
                        pair_values.map(f64::to_bits)
                    ));
                }

                let output = transform.transform_point(point_value);
                transform_point_calls += 1;
                let output_bits = output.map(f64::to_bits);
                if output_bits != expected.output {
                    mismatches.push(format!(
                        "fixture {fixture} repeat {repeat}: output expected {:?}, got {:?}",
                        expected.output, output_bits
                    ));
                }
                if point_value.map(f64::to_bits) != point_bits_before {
                    mismatches.push(format!(
                        "fixture {fixture} repeat {repeat}: point input changed from {:?} to {:?}",
                        point_bits_before,
                        point_value.map(f64::to_bits)
                    ));
                }
            }
        }

        if from_point_pairs_calls != 22 {
            mismatches.push(format!(
                "from_point_pairs call census expected 22, got {from_point_pairs_calls}"
            ));
        }
        if transform_point_calls != 22 {
            mismatches.push(format!(
                "transform_point call census expected 22, got {transform_point_calls}"
            ));
        }
        assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
    }

    #[test]
    fn d2_around_native_literal_product() {
        // Frozen actual-native SetTransform/TransformPoint literals from
        // dev/gap_reports/depict_2d/rotation_about_point.md. No expected value
        // is calculated from this Rust implementation.
        const FIXTURES: [AroundNativeLiteral; 12] = [
            AroundNativeLiteral {
                center: [13840167140234767549, 13841658448896447679],
                angle: 13833581496215065898,
                point: [13837906160088633287, 13842551789691373020],
                matrix: [
                    13815326403393209497,
                    4607136205773496569,
                    4605699752545379376,
                    13830508242628272377,
                    13815326403393209497,
                    13845737163034857373,
                ],
                output: [13841200895337859887, 13842993950389275022],
            },
            AroundNativeLiteral {
                center: [4600774667239816340, 4608144052124252029],
                angle: 13835506422081213124,
                point: [13830699860855537818, 4604930618986332161],
                matrix: [
                    13826841555286448733,
                    4605462196814083240,
                    13823451408332882796,
                    13828834233668859048,
                    13826841555286448733,
                    4612239537689877516,
                ],
                output: [4605901809668386118, 4613127418605116318],
            },
            AroundNativeLiteral {
                center: [0, 0],
                angle: 0,
                point: [0, 9223372036854775808],
                matrix: [4607182418800017408, 0, 0, 0, 4607182418800017408, 0],
                output: [0, 0],
            },
            AroundNativeLiteral {
                center: [9223372036854775808, 0],
                angle: 9223372036854775808,
                point: [9223372036854775808, 0],
                matrix: [4607182418800017408, 0, 0, 0, 4607182418800017408, 0],
                output: [0, 0],
            },
            AroundNativeLiteral {
                center: [0, 9223372036854775808],
                angle: 4609753056924675352,
                point: [4607182418800017408, 13830554455654793216],
                matrix: [
                    4364452196894661639,
                    13830554455654793216,
                    0,
                    4607182418800017408,
                    4364452196894661639,
                    0,
                ],
                output: [4607182418800017408, 4607182418800017407],
            },
            AroundNativeLiteral {
                center: [9223372036854775808, 0],
                angle: 13833125093779451160,
                point: [13830554455654793216, 4607182418800017408],
                matrix: [
                    4364452196894661639,
                    4607182418800017408,
                    0,
                    13830554455654793216,
                    4364452196894661639,
                    0,
                ],
                output: [4607182418800017407, 4607182418800017408],
            },
            AroundNativeLiteral {
                center: [4607182418800017408, 13835058055282163712],
                angle: 4614256656552045848,
                point: [4613937818241073152, 4616189618054758400],
                matrix: [
                    13830554455654793216,
                    13592327833376807943,
                    4611686018427387903,
                    4368955796522032135,
                    13830554455654793216,
                    13839561654909534208,
                ],
                output: [13830554455654793219, 13844065254536904704],
            },
            AroundNativeLiteral {
                center: [13837309855095848960, 4616189618054758400],
                angle: 13837628693406821656,
                point: [4598175219545276416, 13826050856027422720],
                matrix: [
                    13830554455654793216,
                    4368955796522032135,
                    13841813454723219456,
                    13592327833376807943,
                    13830554455654793216,
                    4620693217682128896,
                ],
                output: [13842094929699930112, 4620974692658839552],
            },
            AroundNativeLiteral {
                center: [4591870180066957722, 13819745816549104026],
                angle: 4607394977673999205,
                point: [4604480259023595110, 4606281698874543309],
                matrix: [
                    4602678819172646913,
                    13829347719771606186,
                    13816914319210530680,
                    4605975682916830378,
                    4602678819172646913,
                    13819263122195829213,
                ],
                output: [13826524886406865186, 4606008017307368141],
            },
            AroundNativeLiteral {
                center: [6850974717710472879, 16074346754565248687],
                angle: 4600336947366414254,
                point: [6855478317337843375, 6846471118083102383],
                matrix: [
                    4606572877717793737,
                    13823557941271277023,
                    16066306873668259364,
                    4600185904416501215,
                    4606572877717793737,
                    16069064859436800141,
                ],
                output: [6853120471284824355, 6849333997472211297],
            },
            AroundNativeLiteral {
                center: [4503599627370496, 9227875636482146304],
                angle: 4600336947366414254,
                point: [9007199254740992, 2251799813685248],
                matrix: [
                    4606572877717793737,
                    13823557941271277023,
                    9224695837438312796,
                    4600185904416501215,
                    4606572877717793737,
                    9225305378520536468,
                ],
                output: [6259572026655921, 3423215126666318],
            },
            AroundNativeLiteral {
                center: [12034929314500959671, 6402087712267973674],
                angle: 13823708984221190062,
                point: [6387181665105488315, 15625459749122749482],
                matrix: [
                    4606572877717793737,
                    4600185904416501215,
                    15619025341189567309,
                    13823557941271277023,
                    4606572877717793737,
                    6384650684730343584,
                ],
                output: [15622769946159922631, 15624767819230837210],
            },
        ];

        let mut around_calls = 0;
        let mut transform_point_calls = 0;
        let mut mismatches = Vec::new();

        for repeat in 0..2 {
            for (fixture, expected) in FIXTURES.iter().enumerate() {
                let center = point(expected.center);
                let center_before = center.map(f64::to_bits);
                let angle = f64::from_bits(expected.angle);
                let angle_before = angle.to_bits();
                let input_point = point(expected.point);
                let input_point_before = input_point.map(f64::to_bits);

                let transform = Transform2D::around(center, angle);
                around_calls += 1;
                let matrix_bits = transform.data.map(f64::to_bits);
                for cell in 0..6 {
                    if matrix_bits[cell] != expected.matrix[cell] {
                        mismatches.push(format!(
                            "case {} repeat {repeat}: matrix cell {cell} expected {}, got {}",
                            fixture + 1,
                            expected.matrix[cell],
                            matrix_bits[cell]
                        ));
                    }
                }

                let output = transform.transform_point(input_point);
                transform_point_calls += 1;
                let output_bits = output.map(f64::to_bits);
                for cell in 0..2 {
                    if output_bits[cell] != expected.output[cell] {
                        mismatches.push(format!(
                            "case {} repeat {repeat}: point cell {cell} expected {}, got {}",
                            fixture + 1,
                            expected.output[cell],
                            output_bits[cell]
                        ));
                    }
                }

                if center.map(f64::to_bits) != center_before {
                    mismatches.push(format!(
                        "case {} repeat {repeat}: center input changed from {:?} to {:?}",
                        fixture + 1,
                        center_before,
                        center.map(f64::to_bits)
                    ));
                }
                if angle.to_bits() != angle_before {
                    mismatches.push(format!(
                        "case {} repeat {repeat}: angle input changed from {angle_before} to {}",
                        fixture + 1,
                        angle.to_bits()
                    ));
                }
                if input_point.map(f64::to_bits) != input_point_before {
                    mismatches.push(format!(
                        "case {} repeat {repeat}: point input changed from {:?} to {:?}",
                        fixture + 1,
                        input_point_before,
                        input_point.map(f64::to_bits)
                    ));
                }
            }
        }

        if around_calls != 24 {
            mismatches.push(format!(
                "Transform2D::around call census expected 24, got {around_calls}"
            ));
        }
        if transform_point_calls != 24 {
            mismatches.push(format!(
                "transform_point call census expected 24, got {transform_point_calls}"
            ));
        }
        assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
    }
}
