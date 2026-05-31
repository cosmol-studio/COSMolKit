use std::cell::RefCell;
use std::ops::{Add, Div, DivAssign, Index, IndexMut, Mul, Neg, Sub};

const FUNCTOL: f64 = 1.0e-4;
const MOVETOL: f64 = 1.0e-7;
const EPS: f64 = 3.0e-8;
const TOLX: f64 = 4.0 * EPS;
const MAXSTEP: f64 = 100.0;
const ROW1_ITER1_CHECKPOINT_BITS: [u64; 18] = [
    0x3fe475e01a824896,
    0xbfdedf178888fe31,
    0xbfbd3692c7d4e18a,
    0xbfe5525f9ea5a52d,
    0x3fe89829272a4cde,
    0x3fde4fb92f4fe23d,
    0x3ff650c540ee9981,
    0xbfe16bc86d346ecd,
    0xbfe1912ff74f169f,
    0x3ff05cc01120f8d1,
    0x3fed71cfd920a5eb,
    0x3fe19bf52f5c05e2,
    0xbff066948d8f75d9,
    0xbff226a6db7818d0,
    0x3fdb6ef8a69179f8,
    0xbff5d8b1026e6e2b,
    0x3fde3d51d03c598c,
    0xbfe9434bca030132,
];

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForceFieldVec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl ForceFieldVec3 {
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 0.0 }
    }

    #[must_use]
    pub const fn new4(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }

    #[must_use]
    pub fn length(self) -> f64 {
        // RDKit✔️✔️: double res = x * x + y * y + z * z;
        // RDKit✔️✔️: return sqrt(res);
        self.length_sq().sqrt()
    }

    #[must_use]
    pub const fn length_sq(self) -> f64 {
        // RDKit✔️✔️: double res = x * x + y * y + z * z;
        // RDKit✔️✔️: return res;
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    #[must_use]
    pub const fn dot_product(self, other: Self) -> f64 {
        // RDKit✔️✔️: double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit✔️✔️: return res;
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    #[must_use]
    pub const fn cross_product(self, other: Self) -> Self {
        // RDKit✔️✔️: res.x = y * (other.z) - z * (other.y);
        // RDKit✔️✔️: res.y = -x * (other.z) + z * (other.x);
        // RDKit✔️✔️: res.z = x * (other.y) - y * (other.x);
        Self {
            x: self.y * other.z - self.z * other.y,
            y: -self.x * other.z + self.z * other.x,
            z: self.x * other.y - self.y * other.x,
            w: 0.0,
        }
    }
}

impl Index<usize> for ForceFieldVec3 {
    type Output = f64;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            3 => &self.w,
            _ => panic!("Invalid index on Point3D"),
        }
    }
}

impl IndexMut<usize> for ForceFieldVec3 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            3 => &mut self.w,
            _ => panic!("Invalid index on Point3D"),
        }
    }
}

impl Add for ForceFieldVec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        // RDKit✔️✔️: res.x = p1.x + p2.x;
        // RDKit✔️✔️: res.y = p1.y + p2.y;
        // RDKit✔️✔️: res.z = p1.z + p2.z;
        Self::new4(
            self.x + rhs.x,
            self.y + rhs.y,
            self.z + rhs.z,
            self.w + rhs.w,
        )
    }
}

impl Sub for ForceFieldVec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        // RDKit✔️✔️: res.x = p1.x - p2.x;
        // RDKit✔️✔️: res.y = p1.y - p2.y;
        // RDKit✔️✔️: res.z = p1.z - p2.z;
        Self::new4(
            self.x - rhs.x,
            self.y - rhs.y,
            self.z - rhs.z,
            self.w - rhs.w,
        )
    }
}

impl Neg for ForceFieldVec3 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        // RDKit✔️✔️: res.x *= -1.0;
        // RDKit✔️✔️: res.y *= -1.0;
        // RDKit✔️✔️: res.z *= -1.0;
        Self::new4(-self.x, -self.y, -self.z, -self.w)
    }
}

impl Div<f64> for ForceFieldVec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        // RDKit✔️✔️: res.x = p1.x / v;
        // RDKit✔️✔️: res.y = p1.y / v;
        // RDKit✔️✔️: res.z = p1.z / v;
        Self::new4(self.x / rhs, self.y / rhs, self.z / rhs, self.w / rhs)
    }
}

impl DivAssign<f64> for ForceFieldVec3 {
    fn div_assign(&mut self, rhs: f64) {
        // RDKit✔️✔️: x /= scale;
        // RDKit✔️✔️: y /= scale;
        // RDKit✔️✔️: z /= scale;
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
        self.w /= rhs;
    }
}

impl Mul<f64> for ForceFieldVec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::new4(self.x * rhs, self.y * rhs, self.z * rhs, self.w * rhs)
    }
}

pub fn normalize_angle_deg(angle_deg: &mut f64) {
    // RDKit✔️✔️: angleDeg = fmod(angleDeg, 360.0);
    // RDKit✔️✔️: if (angleDeg < -180.0) {
    // RDKit✔️✔️:   angleDeg += 360.0;
    // RDKit✔️✔️: } else if (angleDeg > 180.0) {
    // RDKit✔️✔️:   angleDeg -= 360.0;
    // RDKit✔️✔️: }
    *angle_deg %= 360.0;
    if *angle_deg < -180.0 {
        *angle_deg += 360.0;
    } else if *angle_deg > 180.0 {
        *angle_deg -= 360.0;
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DihedralOutput {
    pub dihedral: Option<f64>,
    pub cos_phi: f64,
    pub r: [ForceFieldVec3; 4],
    pub t: [ForceFieldVec3; 2],
    pub d: [f64; 2],
}

#[must_use]
pub fn compute_dihedral_from_points(
    p1: &ForceFieldVec3,
    p2: &ForceFieldVec3,
    p3: &ForceFieldVec3,
    p4: &ForceFieldVec3,
    want_dihedral: bool,
) -> DihedralOutput {
    // RDKit✔️✔️: PRECONDITION(p1, "p1 must not be null");
    // RDKit✔️✔️: PRECONDITION(p2, "p2 must not be null");
    // RDKit✔️✔️: PRECONDITION(p3, "p3 must not be null");
    // RDKit✔️✔️: PRECONDITION(p4, "p4 must not be null");
    let mut r = [ForceFieldVec3::default(); 4];
    let mut t = [ForceFieldVec3::default(); 2];
    let mut d = [0.0; 2];
    // RDKit✔️✔️: r[0] = *p1 - *p2;
    // RDKit✔️✔️: r[1] = *p3 - *p2;
    // RDKit✔️✔️: r[2] = -r[1];
    // RDKit✔️✔️: r[3] = *p4 - *p3;
    r[0] = *p1 - *p2;
    r[1] = *p3 - *p2;
    r[2] = -r[1];
    r[3] = *p4 - *p3;
    // RDKit✔️✔️: t[0] = r[0].crossProduct(r[1]);
    // RDKit✔️✔️: d[0] = (std::max)(t[0].length(), 1.0e-5);
    // RDKit✔️✔️: t[0] /= d[0];
    // RDKit✔️✔️: t[1] = r[2].crossProduct(r[3]);
    // RDKit✔️✔️: d[1] = (std::max)(t[1].length(), 1.0e-5);
    // RDKit✔️✔️: t[1] /= d[1];
    t[0] = r[0].cross_product(r[1]);
    d[0] = t[0].length().max(1.0e-5);
    t[0] /= d[0];
    t[1] = r[2].cross_product(r[3]);
    d[1] = t[1].length().max(1.0e-5);
    t[1] /= d[1];
    // RDKit✔️✔️: *cosPhi = (std::max)(-1.0, (std::min)(t[0].dotProduct(t[1]), 1.0));
    let cos_phi = t[0].dot_product(t[1]).clamp(-1.0, 1.0);
    let dihedral = if want_dihedral {
        // RDKit✔️✔️: RDGeom::Point3D m = t[0].crossProduct(r[1]);
        // RDKit✔️✔️: double mLength = (std::max)(m.length(), 1.0e-5);
        // RDKit✔️✔️: *dihedral = -atan2(m.dotProduct(t[1]) / mLength, *cosPhi);
        let m = t[0].cross_product(r[1]);
        let m_length = m.length().max(1.0e-5);
        Some(-(m.dot_product(t[1]) / m_length).atan2(cos_phi))
    } else {
        None
    };
    DihedralOutput {
        dihedral,
        cos_phi,
        r,
        t,
        d,
    }
}

#[must_use]
pub fn compute_dihedral_from_flat(
    pos: &[f64],
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    want_dihedral: bool,
) -> DihedralOutput {
    // RDKit✔️✔️: RDGeom::Point3D p1(pos[3 * idx1], pos[3 * idx1 + 1], pos[3 * idx1 + 2]);
    // RDKit✔️✔️: RDGeom::Point3D p2(pos[3 * idx2], pos[3 * idx2 + 1], pos[3 * idx2 + 2]);
    // RDKit✔️✔️: RDGeom::Point3D p3(pos[3 * idx3], pos[3 * idx3 + 1], pos[3 * idx3 + 2]);
    // RDKit✔️✔️: RDGeom::Point3D p4(pos[3 * idx4], pos[3 * idx4 + 1], pos[3 * idx4 + 2]);
    let p1 = ForceFieldVec3::new(pos[3 * idx1], pos[3 * idx1 + 1], pos[3 * idx1 + 2]);
    let p2 = ForceFieldVec3::new(pos[3 * idx2], pos[3 * idx2 + 1], pos[3 * idx2 + 2]);
    let p3 = ForceFieldVec3::new(pos[3 * idx3], pos[3 * idx3 + 1], pos[3 * idx3 + 2]);
    let p4 = ForceFieldVec3::new(pos[3 * idx4], pos[3 * idx4 + 1], pos[3 * idx4 + 2]);
    compute_dihedral_from_points(&p1, &p2, &p3, &p4, want_dihedral)
}

#[must_use]
pub fn compute_dihedral_from_position_vec(
    pos: &[ForceFieldVec3],
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    want_dihedral: bool,
) -> DihedralOutput {
    // RDKit✔️✔️: computeDihedral(static_cast<RDGeom::Point3D *>(pos[idx1]),
    // RDKit✔️✔️:                 static_cast<RDGeom::Point3D *>(pos[idx2]),
    // RDKit✔️✔️:                 static_cast<RDGeom::Point3D *>(pos[idx3]),
    // RDKit✔️✔️:                 static_cast<RDGeom::Point3D *>(pos[idx4]), dihedral, cosPhi,
    // RDKit✔️✔️:                 r, t, d);
    compute_dihedral_from_points(
        &pos[idx1],
        &pos[idx2],
        &pos[idx3],
        &pos[idx4],
        want_dihedral,
    )
}

pub trait ForceFieldContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib>;
    fn set_force_field(&mut self, _owner: *const ForceField) {}
    fn get_energy(&self, pos: &[f64]) -> f64;
    fn get_grad(&self, pos: &[f64], grad: &mut [f64]);
}

#[derive(Clone, Debug, PartialEq)]
pub struct ForceFieldSnapshot {
    pub positions: Vec<f64>,
    pub energy: f64,
}

pub struct ForceField {
    dimension: usize,
    initialized: bool,
    num_points: usize,
    dist_mat: RefCell<Vec<f64>>,
    positions: Vec<ForceFieldVec3>,
    contribs: Vec<Box<dyn ForceFieldContrib>>,
    fixed_points: Vec<usize>,
    mat_size: usize,
}

impl ForceField {
    #[must_use]
    pub fn new(dimension: usize) -> Self {
        // RDKit✔️✔️: ForceField(unsigned int dimension = 3) : d_dimension(dimension) {}
        Self {
            dimension,
            initialized: false,
            num_points: 0,
            dist_mat: RefCell::new(Vec::new()),
            positions: Vec::new(),
            contribs: Vec::new(),
            fixed_points: Vec::new(),
            mat_size: 0,
        }
    }

    #[must_use]
    pub fn dimension(&self) -> usize {
        // RDKit✔️✔️: unsigned int dimension() const { return d_dimension; }
        self.dimension
    }

    #[must_use]
    pub fn num_points(&self) -> usize {
        // RDKit✔️✔️: unsigned int numPoints() const { return d_numPoints; }
        self.num_points
    }

    pub fn positions(&self) -> &[ForceFieldVec3] {
        // RDKit✔️✔️: const RDGeom::PointPtrVect &positions() const { return d_positions; }
        &self.positions
    }

    pub fn positions_mut(&mut self) -> &mut Vec<ForceFieldVec3> {
        // RDKit✔️✔️: RDGeom::PointPtrVect &positions() { return d_positions; }
        &mut self.positions
    }

    pub fn contribs(&self) -> &[Box<dyn ForceFieldContrib>] {
        // RDKit✔️✔️: const ContribPtrVect &contribs() const { return d_contribs; }
        &self.contribs
    }

    pub fn contribs_mut(&mut self) -> &mut Vec<Box<dyn ForceFieldContrib>> {
        // RDKit✔️✔️: ContribPtrVect &contribs() { return d_contribs; }
        &mut self.contribs
    }

    pub fn add_contrib(&mut self, mut contrib: Box<dyn ForceFieldContrib>) {
        contrib.set_force_field(self as *const ForceField);
        self.contribs.push(contrib);
    }

    fn refresh_contrib_owners(&mut self) {
        let owner = self as *const ForceField;
        for contrib in &mut self.contribs {
            contrib.set_force_field(owner);
        }
    }

    pub fn fixed_points(&self) -> &[usize] {
        &self.fixed_points
    }

    pub fn fixed_points_mut(&mut self) -> &mut Vec<usize> {
        &mut self.fixed_points
    }

    pub fn initialize(&mut self) {
        // RDKit✔️✔️: df_init = false;
        // RDKit✔️✔️: delete[] dp_distMat;
        // RDKit✔️✔️: dp_distMat = nullptr;
        // Rust ForceField values can move after contrib construction; refresh
        // RDKit-style owner pointers before any energy or gradient evaluation.
        self.refresh_contrib_owners();
        self.initialized = false;
        self.dist_mat.get_mut().clear();
        // RDKit✔️✔️: d_numPoints = d_positions.size();
        // RDKit✔️✔️: d_matSize = d_numPoints * (d_numPoints + 1) / 2;
        // RDKit✔️✔️: dp_distMat = new double[d_matSize];
        self.num_points = self.positions.len();
        self.mat_size = self.num_points * (self.num_points + 1) / 2;
        *self.dist_mat.get_mut() = vec![0.0; self.mat_size];
        self.init_distance_matrix();
        // RDKit✔️✔️: df_init = true;
        self.initialized = true;
    }

    pub fn distance(&self, mut i: usize, mut j: usize, pos: Option<&[f64]>) -> f64 {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: URANGE_CHECK(i, d_numPoints);
        // RDKit✔️✔️: URANGE_CHECK(j, d_numPoints);
        assert!(self.initialized, "not initialized");
        assert!(i < self.num_points);
        assert!(j < self.num_points);
        // RDKit✔️✔️: if (j < i) {
        // RDKit✔️✔️:   int tmp = j;
        // RDKit✔️✔️:   j = i;
        // RDKit✔️✔️:   i = tmp;
        // RDKit✔️✔️: }
        if j < i {
            std::mem::swap(&mut i, &mut j);
        }
        // RDKit✔️✔️: unsigned int idx = i + j * (j + 1) / 2;
        // RDKit✔️✔️: CHECK_INVARIANT(idx < d_matSize, "Bad index");
        let mat_idx = i + j * (j + 1) / 2;
        assert!(mat_idx < self.mat_size, "Bad index");
        // RDKit✔️✔️: double &res = dp_distMat[idx];
        // RDKit✔️✔️: if (res < 0.0) {
        let mut dist_mat = self.dist_mat.borrow_mut();
        if dist_mat[mat_idx] < 0.0 {
            let res = self.distance2_no_cache(i, j, pos).sqrt();
            dist_mat[mat_idx] = res;
        }
        dist_mat[mat_idx]
    }

    pub fn distance_const(&self, i: usize, j: usize, pos: Option<&[f64]>) -> f64 {
        // RDKit✔️✔️: auto res = sqrt(distance2(i, j, pos));
        // RDKit✔️✔️: return res;
        self.distance2(i, j, pos).sqrt()
    }

    pub fn distance2(&self, mut i: usize, mut j: usize, pos: Option<&[f64]>) -> f64 {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: URANGE_CHECK(i, d_numPoints);
        // RDKit✔️✔️: URANGE_CHECK(j, d_numPoints);
        assert!(self.initialized, "not initialized");
        assert!(i < self.num_points);
        assert!(j < self.num_points);
        // RDKit✔️✔️: if (j < i) {
        if j < i {
            std::mem::swap(&mut i, &mut j);
        }
        self.distance2_no_cache(i, j, pos)
    }

    fn distance2_no_cache(&self, i: usize, j: usize, pos: Option<&[f64]>) -> f64 {
        let mut res = 0.0;
        if let Some(pos) = pos {
            // RDKit✔️✔️: double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
            // RDKit✔️✔️: for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
            // RDKit✔️✔️:   double tmp = *pi - *pj;
            // RDKit✔️✔️:   res += tmp * tmp;
            // RDKit✔️✔️: }
            for idx in 0..self.dimension {
                let tmp = pos[self.dimension * i + idx] - pos[self.dimension * j + idx];
                res += tmp * tmp;
            }
        } else {
            // RDKit✔️✔️: double tmp =
            // RDKit✔️✔️:     (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
            // RDKit✔️✔️: res += tmp * tmp;
            for idx in 0..self.dimension {
                let tmp = self.positions[i][idx] - self.positions[j][idx];
                res += tmp * tmp;
            }
        }
        res
    }

    pub fn calc_energy_current(&self, contrib_energies: Option<&mut Vec<f64>>) -> f64 {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        assert!(self.initialized, "not initialized");
        let mut res = 0.0;
        // RDKit✔️✔️: if (d_contribs.empty()) {
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        if self.contribs.is_empty() {
            return res;
        }
        let mut contrib_energies = contrib_energies;
        if let Some(contrib_energies) = &mut contrib_energies {
            // RDKit✔️✔️: contribs->clear();
            // RDKit✔️✔️: contribs->reserve(d_contribs.size());
            contrib_energies.clear();
            contrib_energies.reserve(self.contribs.len());
        }
        // RDKit✔️✔️: auto *pos = new double[d_dimension * N];
        // RDKit✔️✔️: this->scatter(pos);
        let mut pos = vec![0.0; self.dimension * self.positions.len()];
        self.scatter(&mut pos);
        // RDKit✔️✔️: for (const auto &d_contrib : d_contribs) {
        // RDKit✔️✔️:   double e = d_contrib->getEnergy(pos);
        // RDKit✔️✔️:   res += e;
        // RDKit✔️✔️:   if (contribs) {
        // RDKit✔️✔️:     contribs->push_back(e);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        for contrib in &self.contribs {
            let energy = contrib.get_energy(&pos);
            res += energy;
            if let Some(contrib_energies) = &mut contrib_energies {
                contrib_energies.push(energy);
            }
        }
        res
    }

    pub fn calc_energy(&mut self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: PRECONDITION(pos, "bad position vector");
        assert!(self.initialized, "not initialized");
        assert!(!pos.is_empty(), "bad position vector");
        // RDKit✔️✔️: this->initDistanceMatrix();
        self.init_distance_matrix();
        let mut res = 0.0;
        // RDKit✔️✔️: if (d_contribs.empty()) {
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        if self.contribs.is_empty() {
            return res;
        }
        // RDKit✔️✔️: for (ContribPtrVect::const_iterator contrib = d_contribs.begin();
        // RDKit✔️✔️:      contrib != d_contribs.end(); contrib++) {
        // RDKit✔️✔️:   double E = (*contrib)->getEnergy(pos);
        // RDKit✔️✔️:   res += E;
        // RDKit✔️✔️: }
        for contrib in &self.contribs {
            res += contrib.get_energy(pos);
        }
        res
    }

    pub fn calc_grad_current(&self, grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: PRECONDITION(grad, "bad gradient vector");
        assert!(self.initialized, "not initialized");
        assert!(!grad.is_empty(), "bad gradient vector");
        // RDKit✔️✔️: if (d_contribs.empty()) {
        // RDKit✔️✔️:   return;
        // RDKit✔️✔️: }
        if self.contribs.is_empty() {
            return;
        }
        // RDKit✔️✔️: auto *pos = new double[d_dimension * N];
        // RDKit✔️✔️: this->scatter(pos);
        let mut pos = vec![0.0; self.dimension * self.positions.len()];
        self.scatter(&mut pos);
        // RDKit✔️✔️: for (const auto &d_contrib : d_contribs) {
        // RDKit✔️✔️:   d_contrib->getGrad(pos, grad);
        // RDKit✔️✔️: }
        for contrib in &self.contribs {
            contrib.get_grad(&pos, grad);
        }
        self.zero_fixed_point_grads(grad);
    }

    pub fn calc_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: PRECONDITION(pos, "bad position vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad gradient vector");
        assert!(self.initialized, "not initialized");
        assert!(!pos.is_empty(), "bad position vector");
        assert!(!grad.is_empty(), "bad gradient vector");
        if self.contribs.is_empty() {
            return;
        }
        // RDKit✔️✔️: (*contrib)->getGrad(pos, grad);
        let trace_row1_iter1 =
            row1_bfgs_trace_enabled(pos.len()) && row1_iter1_checkpoint_pos_matches(pos);
        for (idx, contrib) in self.contribs.iter().enumerate() {
            contrib.get_grad(pos, grad);
            if trace_row1_iter1 {
                let bits: Vec<String> = grad
                    .iter()
                    .map(|value| format!("{:#018x}", value.to_bits()))
                    .collect();
                println!(
                    "[row1-calcgrad-cumulative] contrib_idx={idx} values=[{}]",
                    bits.join(",")
                );
            }
        }
        self.zero_fixed_point_grads(grad);
    }

    fn zero_fixed_point_grads(&self, grad: &mut [f64]) {
        // RDKit✔️✔️: for (int d_fixedPoint : d_fixedPoints) {
        // RDKit✔️✔️:   CHECK_INVARIANT(static_cast<unsigned int>(d_fixedPoint) < d_numPoints,
        // RDKit✔️✔️:                   "bad fixed point index");
        // RDKit✔️✔️:   unsigned int idx = d_dimension * d_fixedPoint;
        // RDKit✔️✔️:   for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:     grad[idx + di] = 0.0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        for &fixed_point in &self.fixed_points {
            assert!(fixed_point < self.num_points, "bad fixed point index");
            let idx = self.dimension * fixed_point;
            for di in 0..self.dimension {
                grad[idx + di] = 0.0;
            }
        }
    }

    pub fn scatter(&self, pos: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: PRECONDITION(pos, "bad position vector");
        assert!(self.initialized, "not initialized");
        assert!(!pos.is_empty(), "bad position vector");
        // RDKit✔️✔️: unsigned int tab = 0;
        // RDKit✔️✔️: for (auto d_position : d_positions) {
        // RDKit✔️✔️:   for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:     pos[tab + di] = (*d_position)[di];  //->x;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   tab += this->dimension();
        // RDKit✔️✔️: }
        let mut tab = 0;
        for position in &self.positions {
            for di in 0..self.dimension {
                pos[tab + di] = position[di];
            }
            tab += self.dimension;
        }
        // RDKit✔️✔️: POSTCONDITION(tab == this->dimension() * d_positions.size(), "bad index");
        assert_eq!(tab, self.dimension * self.positions.len(), "bad index");
    }

    pub fn gather(&mut self, pos: &[f64]) {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: PRECONDITION(pos, "bad position vector");
        assert!(self.initialized, "not initialized");
        assert!(!pos.is_empty(), "bad position vector");
        // RDKit✔️✔️: unsigned int tab = 0;
        // RDKit✔️✔️: for (auto &d_position : d_positions) {
        // RDKit✔️✔️:   for (unsigned int di = 0; di < this->dimension(); ++di) {
        // RDKit✔️✔️:     (*d_position)[di] = pos[tab + di];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   tab += this->dimension();
        // RDKit✔️✔️: }
        let mut tab = 0;
        for position in &mut self.positions {
            for di in 0..self.dimension {
                position[di] = pos[tab + di];
            }
            tab += self.dimension;
        }
    }

    pub fn init_distance_matrix(&mut self) {
        // RDKit✔️✔️: PRECONDITION(d_numPoints, "no points");
        // RDKit✔️✔️: PRECONDITION(dp_distMat, "no distance matrix");
        // RDKit✔️✔️: PRECONDITION(static_cast<unsigned int>(d_numPoints * (d_numPoints + 1) / 2) <=
        // RDKit✔️✔️:                  d_matSize,
        // RDKit✔️✔️:              "matrix size mismatch");
        assert!(self.num_points != 0, "no points");
        assert!(!self.dist_mat.borrow().is_empty(), "no distance matrix");
        assert!(
            self.num_points * (self.num_points + 1) / 2 <= self.mat_size,
            "matrix size mismatch"
        );
        // RDKit✔️✔️: for (unsigned int i = 0; i < d_numPoints * (d_numPoints + 1) / 2; i++) {
        // RDKit✔️✔️:   dp_distMat[i] = -1.0;
        // RDKit✔️✔️: }
        for value in self
            .dist_mat
            .borrow_mut()
            .iter_mut()
            .take(self.num_points * (self.num_points + 1) / 2)
        {
            *value = -1.0;
        }
    }

    pub fn minimize(&mut self, max_its: usize, force_tol: f64, energy_tol: f64) -> i32 {
        // RDKit✔️✔️: return minimize(0, nullptr, maxIts, forceTol, energyTol);
        self.minimize_with_snapshots(0, None, max_its, force_tol, energy_tol)
    }

    pub fn minimize_with_snapshots(
        &mut self,
        snapshot_freq: usize,
        mut snapshot_vect: Option<&mut Vec<ForceFieldSnapshot>>,
        max_its: usize,
        force_tol: f64,
        energy_tol: f64,
    ) -> i32 {
        // RDKit✔️✔️: PRECONDITION(df_init, "not initialized");
        // RDKit✔️✔️: PRECONDITION(static_cast<unsigned int>(d_numPoints) == d_positions.size(),
        // RDKit✔️✔️:              "size mismatch");
        assert!(self.initialized, "not initialized");
        assert_eq!(self.num_points, self.positions.len(), "size mismatch");
        // RDKit✔️✔️: if (d_contribs.empty()) {
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
        if self.contribs.is_empty() {
            return 0;
        }
        // RDKit✔️✔️: unsigned int dim = this->d_numPoints * d_dimension;
        // RDKit✔️✔️: std::vector<double> points(dim);
        // RDKit✔️✔️: this->scatter(points.data());
        let dim = self.num_points * self.dimension;
        let mut points = vec![0.0; dim];
        self.scatter(&mut points);
        let (res, final_points) = {
            let ff_cell = std::cell::RefCell::new(&mut *self);
            bfgs_minimize(
                points,
                force_tol,
                |pos| ff_cell.borrow_mut().calc_energy(pos),
                |pos, grad| {
                    let ff_ref = ff_cell.borrow();
                    calc_gradient_wrapper(&ff_ref, pos, grad)
                },
                snapshot_freq,
                snapshot_vect.as_deref_mut(),
                energy_tol,
                max_its,
            )
        };
        // RDKit✔️✔️: this->gather(points.data());
        self.gather(&final_points);
        res
    }

    #[must_use]
    pub fn dist_mat_value(&self, i: usize) -> f64 {
        self.dist_mat.borrow()[i]
    }
}

impl Clone for ForceField {
    fn clone(&self) -> Self {
        // RDKit✔️✔️: : d_dimension(other.d_dimension),
        // RDKit✔️✔️:   df_init(false),
        // RDKit✔️✔️:   d_numPoints(other.d_numPoints),
        // RDKit✔️✔️:   dp_distMat(nullptr) {
        let mut cloned = Self {
            dimension: self.dimension,
            initialized: false,
            num_points: self.num_points,
            dist_mat: RefCell::new(Vec::new()),
            positions: self.positions.clone(),
            contribs: Vec::new(),
            fixed_points: self.fixed_points.clone(),
            mat_size: 0,
        };
        // RDKit✔️✔️: for (const auto &contrib : other.d_contribs) {
        // RDKit✔️✔️:   ForceFieldContrib *ncontrib = contrib->copy();
        // RDKit✔️✔️:   ncontrib->dp_forceField = this;
        // RDKit✔️✔️:   d_contribs.push_back(ContribPtr(ncontrib));
        // RDKit✔️✔️: }
        for contrib in &self.contribs {
            let mut copied = contrib.copy();
            copied.set_force_field(&cloned as *const ForceField);
            cloned.contribs.push(copied);
        }
        cloned
    }
}

impl Drop for ForceField {
    fn drop(&mut self) {
        // RDKit✔️✔️: d_numPoints = 0;
        // RDKit✔️✔️: d_positions.clear();
        // RDKit✔️✔️: d_contribs.clear();
        // RDKit✔️✔️: delete[] dp_distMat;
        // RDKit✔️✔️: dp_distMat = nullptr;
        self.num_points = 0;
        self.positions.clear();
        self.contribs.clear();
        self.dist_mat.get_mut().clear();
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DistanceConstraintContrib {
    owner: *const ForceField,
    end1_idx: usize,
    end2_idx: usize,
    min_len: f64,
    max_len: f64,
    force_constant: f64,
}

impl DistanceConstraintContrib {
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, owner->positions().size());
        // RDKit✔️✔️: PRECONDITION(maxLen >= minLen, "bad bounds");
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(max_len >= min_len, "bad bounds");
        // RDKit✔️✔️: dp_forceField = owner;
        // RDKit✔️✔️: d_end1Idx = idx1;
        // RDKit✔️✔️: d_end2Idx = idx2;
        // RDKit✔️✔️: d_minLen = minLen;
        // RDKit✔️✔️: d_maxLen = maxLen;
        // RDKit✔️✔️: d_forceConstant = forceConst;
        Self {
            owner,
            end1_idx: idx1,
            end2_idx: idx2,
            min_len,
            max_len,
            force_constant,
        }
    }

    pub fn new_relative(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        relative: bool,
        mut min_len: f64,
        mut max_len: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: const RDGeom::PointPtrVect &pos = owner->positions();
        // RDKit✔️✔️: URANGE_CHECK(idx1, pos.size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, pos.size());
        // RDKit✔️✔️: PRECONDITION(maxLen >= minLen, "bad bounds");
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(max_len >= min_len, "bad bounds");
        // RDKit✔️✔️: if (relative) {
        if relative {
            // RDKit✔️✔️: RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
            // RDKit✔️✔️: RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
            // RDKit✔️✔️: const auto dist = (p1 - p2).length();
            // RDKit✔️✔️: minLen = std::max(dist + minLen, 0.0);
            // RDKit✔️✔️: maxLen = std::max(dist + maxLen, 0.0);
            // RDKit✔️✔️: }
            let dist = (owner.positions()[idx1] - owner.positions()[idx2]).length();
            min_len = (dist + min_len).max(0.0);
            max_len = (dist + max_len).max(0.0);
        }
        // RDKit✔️✔️: d_minLen = minLen;
        // RDKit✔️✔️: d_maxLen = maxLen;
        // RDKit✔️✔️: dp_forceField = owner;
        // RDKit✔️✔️: d_end1Idx = idx1;
        // RDKit✔️✔️: d_end2Idx = idx2;
        // RDKit✔️✔️: d_forceConstant = forceConst;
        Self {
            owner,
            end1_idx: idx1,
            end2_idx: idx2,
            min_len,
            max_len,
            force_constant,
        }
    }

    #[must_use]
    pub fn min_len(&self) -> f64 {
        self.min_len
    }

    #[must_use]
    pub fn max_len(&self) -> f64 {
        self.max_len
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: ForceFieldContrib objects store the owning ForceField pointer just
        // like RDKit's dp_forceField. Constructors and set_force_field maintain it.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for DistanceConstraintContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️❌: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️: double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
        //
        // COSMolKit stores cached distances from the initialized point set in
        // `dist_mat`. For contrib evaluation against an explicit `pos` buffer,
        // match RDKit's live-position behavior by recomputing the current pair
        // distance instead of reusing the initialized cache entry.
        let dist = self
            .owner()
            .distance2(self.end1_idx, self.end2_idx, Some(pos))
            .sqrt();
        // RDKit✔️✔️: double distTerm = 0.0;
        let mut dist_term = 0.0;
        // RDKit✔️✔️: if (dist < d_minLen) {
        // RDKit✔️✔️:   distTerm = d_minLen - dist;
        // RDKit✔️✔️: } else if (dist > d_maxLen) {
        // RDKit✔️✔️:   distTerm = dist - d_maxLen;
        // RDKit✔️✔️: }
        if dist < self.min_len {
            dist_term = self.min_len - dist;
        } else if dist > self.max_len {
            dist_term = dist - self.max_len;
        }
        // RDKit✔️✔️: double res = 0.5 * d_forceConstant * distTerm * distTerm;
        // RDKit✔️✔️: return res;
        0.5 * self.force_constant * dist_term * dist_term
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️❌: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️: double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
        //
        // See `get_energy()`: gradients must use the current trial coordinates,
        // not the initialization-time cached distance.
        let dist = self
            .owner()
            .distance2(self.end1_idx, self.end2_idx, Some(pos))
            .sqrt();
        // RDKit✔️✔️: double preFactor = 0.0;
        // RDKit✔️✔️: if (dist < d_minLen) {
        // RDKit✔️✔️:   preFactor = dist - d_minLen;
        // RDKit✔️✔️: } else if (dist > d_maxLen) {
        // RDKit✔️✔️:   preFactor = dist - d_maxLen;
        // RDKit✔️✔️: } else {
        // RDKit✔️✔️:   return;
        // RDKit✔️✔️: }
        let mut pre_factor = if dist < self.min_len {
            dist - self.min_len
        } else if dist > self.max_len {
            dist - self.max_len
        } else {
            return;
        };
        // RDKit✔️✔️: preFactor *= d_forceConstant;
        pre_factor *= self.force_constant;
        // RDKit✔️✔️: double *end1Coords = &(pos[3 * d_end1Idx]);
        // RDKit✔️✔️: double *end2Coords = &(pos[3 * d_end2Idx]);
        // RDKit✔️✔️: for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:   double dGrad =
        // RDKit✔️✔️:       preFactor * (end1Coords[i] - end2Coords[i]) / std::max(dist, 1.0e-8);
        // RDKit✔️✔️:   grad[3 * d_end1Idx + i] += dGrad;
        // RDKit✔️✔️:   grad[3 * d_end2Idx + i] -= dGrad;
        // RDKit✔️✔️: }
        let end1 = 3 * self.end1_idx;
        let end2 = 3 * self.end2_idx;
        for i in 0..3 {
            let d_grad = pre_factor * (pos[end1 + i] - pos[end2 + i]) / dist.max(1.0e-8);
            grad[end1 + i] += d_grad;
            grad[end2 + i] -= d_grad;
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistanceConstraintContribsParams {
    idx1: usize,
    idx2: usize,
    min_len: f64,
    max_len: f64,
    force_constant: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DistanceConstraintContribs {
    owner: *const ForceField,
    contribs: Vec<DistanceConstraintContribsParams>,
}

impl DistanceConstraintContribs {
    pub fn new(owner: &ForceField) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: dp_forceField = owner;
        Self {
            owner,
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(
        &mut self,
        idx1: usize,
        idx2: usize,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) {
        // RDKit✔️✔️: URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️: PRECONDITION(maxLen >= minLen, "bad bounds");
        let owner = self.owner();
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(max_len >= min_len, "bad bounds");
        // RDKit✔️✔️: d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
        self.contribs.push(DistanceConstraintContribsParams {
            idx1,
            idx2,
            min_len,
            max_len,
            force_constant,
        });
    }

    pub fn add_contrib_relative(
        &mut self,
        idx1: usize,
        idx2: usize,
        relative: bool,
        mut min_len: f64,
        mut max_len: f64,
        force_constant: f64,
    ) {
        // RDKit✔️✔️: const RDGeom::PointPtrVect &pos = dp_forceField->positions();
        // RDKit✔️✔️: URANGE_CHECK(idx1, pos.size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, pos.size());
        // RDKit✔️✔️: PRECONDITION(maxLen >= minLen, "bad bounds");
        let owner = self.owner();
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(max_len >= min_len, "bad bounds");
        // RDKit✔️✔️: if (relative) {
        if relative {
            // RDKit✔️✔️: const RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
            // RDKit✔️✔️: const RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
            // RDKit✔️✔️: const auto distance = (p1 - p2).length();
            // RDKit✔️✔️: minLen = std::max(minLen + distance, 0.0);
            // RDKit✔️✔️: maxLen = std::max(maxLen + distance, 0.0);
            // RDKit✔️✔️: }
            let distance = (owner.positions()[idx1] - owner.positions()[idx2]).length();
            min_len = (min_len + distance).max(0.0);
            max_len = (max_len + distance).max(0.0);
        }
        // RDKit✔️✔️: d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
        self.contribs.push(DistanceConstraintContribsParams {
            idx1,
            idx2,
            min_len,
            max_len,
            force_constant,
        });
    }

    #[must_use]
    pub fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the RDKit owner pointer lifetime convention for contribs.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for DistanceConstraintContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️: double accum = 0.0;
        let mut accum = 0.0;
        // RDKit✔️✔️: for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            // RDKit✔️✔️: const auto distance2 =
            // RDKit✔️✔️:     dp_forceField->distance2(contrib.idx1, contrib.idx2, pos);
            let distance2 = self
                .owner()
                .distance2(contrib.idx1, contrib.idx2, Some(pos));
            // RDKit✔️✔️: double difference = 0.0;
            // RDKit✔️✔️: if (distance2 < contrib.minLen * contrib.minLen) {
            // RDKit✔️✔️:   difference = contrib.minLen - std::sqrt(distance2);
            // RDKit✔️✔️: } else if (distance2 > contrib.maxLen * contrib.maxLen) {
            // RDKit✔️✔️:   difference = std::sqrt(distance2) - contrib.maxLen;
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   continue;
            // RDKit✔️✔️: }
            let difference = if distance2 < contrib.min_len * contrib.min_len {
                contrib.min_len - distance2.sqrt()
            } else if distance2 > contrib.max_len * contrib.max_len {
                distance2.sqrt() - contrib.max_len
            } else {
                continue;
            };
            // RDKit✔️✔️: accum += 0.5 * contrib.forceConstant * difference * difference;
            accum += 0.5 * contrib.force_constant * difference * difference;
        }
        // RDKit✔️✔️: return accum;
        accum
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️: for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            // RDKit✔️✔️: double preFactor = 0.0;
            // RDKit✔️✔️: double distance = 0.0;
            // RDKit✔️✔️: const auto distance2 =
            // RDKit✔️✔️:     dp_forceField->distance2(contrib.idx1, contrib.idx2, pos);
            let distance2 = self
                .owner()
                .distance2(contrib.idx1, contrib.idx2, Some(pos));
            // RDKit✔️✔️: if (distance2 < contrib.minLen * contrib.minLen) {
            // RDKit✔️✔️:   distance = std::sqrt(distance2);
            // RDKit✔️✔️:   preFactor = distance - contrib.minLen;
            // RDKit✔️✔️: } else if (distance2 > contrib.maxLen * contrib.maxLen) {
            // RDKit✔️✔️:   distance = std::sqrt(distance2);
            // RDKit✔️✔️:   preFactor = distance - contrib.maxLen;
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   continue;
            // RDKit✔️✔️: }
            let (distance, mut pre_factor) = if distance2 < contrib.min_len * contrib.min_len {
                let distance = distance2.sqrt();
                (distance, distance - contrib.min_len)
            } else if distance2 > contrib.max_len * contrib.max_len {
                let distance = distance2.sqrt();
                (distance, distance - contrib.max_len)
            } else {
                continue;
            };
            // RDKit✔️✔️: preFactor *= contrib.forceConstant;
            // RDKit✔️✔️: preFactor /= std::max(1.0e-8, distance);
            pre_factor *= contrib.force_constant;
            pre_factor /= distance.max(1.0e-8);
            // RDKit✔️✔️: const double *atom1Coords = &(pos[3 * contrib.idx1]);
            // RDKit✔️✔️: const double *atom2Coords = &(pos[3 * contrib.idx2]);
            // RDKit✔️✔️: for (unsigned int i = 0; i < 3; i++) {
            // RDKit✔️✔️:   const double dGrad = preFactor * (atom1Coords[i] - atom2Coords[i]);
            // RDKit✔️✔️:   grad[3 * contrib.idx1 + i] += dGrad;
            // RDKit✔️✔️:   grad[3 * contrib.idx2 + i] -= dGrad;
            // RDKit✔️✔️: }
            let atom1 = 3 * contrib.idx1;
            let atom2 = 3 * contrib.idx2;
            for i in 0..3 {
                let d_grad = pre_factor * (pos[atom1 + i] - pos[atom2 + i]);
                grad[atom1 + i] += d_grad;
                grad[atom2 + i] -= d_grad;
            }
        }
    }
}

const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

fn point_from_flat(pos: &[f64], idx: usize) -> ForceFieldVec3 {
    ForceFieldVec3::new(pos[3 * idx], pos[3 * idx + 1], pos[3 * idx + 2])
}

fn angle_degrees(p1: ForceFieldVec3, p2: ForceFieldVec3, p3: ForceFieldVec3) -> f64 {
    let r = [p1 - p2, p3 - p2];
    let r_length_sq = [r[0].length_sq().max(1.0e-5), r[1].length_sq().max(1.0e-5)];
    let cos_theta =
        (r[0].dot_product(r[1]) / (r_length_sq[0] * r_length_sq[1]).sqrt()).clamp(-1.0, 1.0);
    RAD2DEG * cos_theta.acos()
}

#[derive(Clone, Debug, PartialEq)]
pub struct AngleConstraintContrib {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    at3_idx: usize,
    min_angle_deg: f64,
    max_angle_deg: f64,
    force_constant: f64,
}

impl AngleConstraintContrib {
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit✔️✔️: RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        // RDKit✔️✔️: URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, owner->positions().size());
        // RDKit✔️✔️: URANGE_CHECK(idx3, owner->positions().size());
        // RDKit✔️✔️: PRECONDITION(!(minAngleDeg > maxAngleDeg),
        // RDKit✔️✔️:              "minAngleDeg must be <= maxAngleDeg");
        assert!((0.0..=180.0).contains(&min_angle_deg));
        assert!((0.0..=180.0).contains(&max_angle_deg));
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(
            min_angle_deg <= max_angle_deg,
            "minAngleDeg must be <= maxAngleDeg"
        );
        // RDKit✔️✔️: dp_forceField = owner;
        // RDKit✔️✔️: d_at1Idx = idx1;
        // RDKit✔️✔️: d_at2Idx = idx2;
        // RDKit✔️✔️: d_at3Idx = idx3;
        // RDKit✔️✔️: d_minAngleDeg = minAngleDeg;
        // RDKit✔️✔️: d_maxAngleDeg = maxAngleDeg;
        // RDKit✔️✔️: d_forceConstant = forceConst;
        Self {
            owner,
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        }
    }

    pub fn new_relative(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        relative: bool,
        mut min_angle_deg: f64,
        mut max_angle_deg: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: const RDGeom::PointPtrVect &pos = owner->positions();
        // RDKit✔️✔️: URANGE_CHECK(idx1, pos.size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, pos.size());
        // RDKit✔️✔️: URANGE_CHECK(idx3, pos.size());
        // RDKit✔️✔️: PRECONDITION(!(minAngleDeg > maxAngleDeg),
        // RDKit✔️✔️:              "minAngleDeg must be <= maxAngleDeg");
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(
            min_angle_deg <= max_angle_deg,
            "minAngleDeg must be <= maxAngleDeg"
        );
        // RDKit✔️✔️: if (relative) {
        if relative {
            // RDKit✔️✔️: const RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
            // RDKit✔️✔️: const RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
            // RDKit✔️✔️: const RDGeom::Point3D &p3 = *((RDGeom::Point3D *)pos[idx3]);
            // RDKit✔️✔️: const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
            // RDKit✔️✔️: const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
            // RDKit✔️✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
            // RDKit✔️✔️: double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
            // RDKit✔️✔️: cosTheta = std::clamp(cosTheta, -1.0, 1.0);
            // RDKit✔️✔️: const double angle = RAD2DEG * acos(cosTheta);
            // RDKit✔️✔️: minAngleDeg += angle;
            // RDKit✔️✔️: maxAngleDeg += angle;
            // RDKit✔️✔️: }
            let angle = angle_degrees(
                owner.positions()[idx1],
                owner.positions()[idx2],
                owner.positions()[idx3],
            );
            min_angle_deg += angle;
            max_angle_deg += angle;
        }
        // RDKit✔️✔️: RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit✔️✔️: RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        assert!((0.0..=180.0).contains(&min_angle_deg));
        assert!((0.0..=180.0).contains(&max_angle_deg));
        Self {
            owner,
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        }
    }

    #[must_use]
    pub fn compute_angle_term(&self, angle: f64) -> f64 {
        // RDKit✔️✔️: double angleTerm = 0.0;
        let mut angle_term = 0.0;
        // RDKit✔️✔️: if (angle < d_minAngleDeg) {
        // RDKit✔️✔️:   angleTerm = angle - d_minAngleDeg;
        // RDKit✔️✔️: } else if (angle > d_maxAngleDeg) {
        // RDKit✔️✔️:   angleTerm = angle - d_maxAngleDeg;
        // RDKit✔️✔️: }
        if angle < self.min_angle_deg {
            angle_term = angle - self.min_angle_deg;
        } else if angle > self.max_angle_deg {
            angle_term = angle - self.max_angle_deg;
        }
        // RDKit✔️✔️: return angleTerm;
        angle_term
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the RDKit owner pointer lifetime convention for contribs.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for AngleConstraintContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️: const RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️: const RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️: const RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at3Idx + 2]);
        let p1 = point_from_flat(pos, self.at1_idx);
        let p2 = point_from_flat(pos, self.at2_idx);
        let p3 = point_from_flat(pos, self.at3_idx);
        // RDKit✔️✔️: const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit✔️✔️: const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit✔️✔️:                                std::max(1.0e-5, r[1].lengthSq())};
        // RDKit✔️✔️: double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit✔️✔️: cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit✔️✔️: const double angle = RAD2DEG * acos(cosTheta);
        let angle = angle_degrees(p1, p2, p3);
        // RDKit✔️✔️: const double angleTerm = computeAngleTerm(angle);
        let angle_term = self.compute_angle_term(angle);
        // RDKit✔️✔️: return d_forceConstant * angleTerm * angleTerm;
        self.force_constant * angle_term * angle_term
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️: const RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️: const RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️: const RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at3Idx + 2]);
        let p1 = point_from_flat(pos, self.at1_idx);
        let p2 = point_from_flat(pos, self.at2_idx);
        let p3 = point_from_flat(pos, self.at3_idx);
        // RDKit✔️✔️: const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        let r = [p1 - p2, p3 - p2];
        // RDKit✔️✔️: const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit✔️✔️:                                std::max(1.0e-5, r[1].lengthSq())};
        let r_length_sq = [r[0].length_sq().max(1.0e-5), r[1].length_sq().max(1.0e-5)];
        // RDKit✔️✔️: double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit✔️✔️: cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit✔️✔️: const double angle = RAD2DEG * acos(cosTheta);
        let cos_theta =
            (r[0].dot_product(r[1]) / (r_length_sq[0] * r_length_sq[1]).sqrt()).clamp(-1.0, 1.0);
        let angle = RAD2DEG * cos_theta.acos();
        // RDKit✔️✔️: const double angleTerm = computeAngleTerm(angle);
        let angle_term = self.compute_angle_term(angle);
        // RDKit✔️✔️: double dE_dTheta = 2.0 * RAD2DEG * d_forceConstant * angleTerm;
        let d_e_d_theta = 2.0 * RAD2DEG * self.force_constant * angle_term;
        // RDKit✔️✔️: RDGeom::Point3D rp = r[1].crossProduct(r[0]);
        // RDKit✔️✔️: double prefactor = dE_dTheta / std::max(1.0e-5, rp.length());
        // RDKit✔️✔️: double t[2] = {-prefactor / rLengthSq[0], prefactor / rLengthSq[1]};
        let rp = r[1].cross_product(r[0]);
        let prefactor = d_e_d_theta / rp.length().max(1.0e-5);
        let t = [-prefactor / r_length_sq[0], prefactor / r_length_sq[1]];
        // RDKit✔️✔️: RDGeom::Point3D dedp[3];
        // RDKit✔️✔️: dedp[0] = r[0].crossProduct(rp) * t[0];
        // RDKit✔️✔️: dedp[2] = r[1].crossProduct(rp) * t[1];
        // RDKit✔️✔️: dedp[1] = -dedp[0] - dedp[2];
        let dedp0 = r[0].cross_product(rp) * t[0];
        let dedp2 = r[1].cross_product(rp) * t[1];
        let dedp = [dedp0, -dedp0 - dedp2, dedp2];
        // RDKit✔️✔️: double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit✔️✔️:                 &(grad[3 * d_at3Idx])};
        // RDKit✔️✔️: for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:   g[i][0] += dedp[i].x;
        // RDKit✔️✔️:   g[i][1] += dedp[i].y;
        // RDKit✔️✔️:   g[i][2] += dedp[i].z;
        // RDKit✔️✔️: }
        for (atom_idx, delta) in [
            (self.at1_idx, dedp[0]),
            (self.at2_idx, dedp[1]),
            (self.at3_idx, dedp[2]),
        ] {
            grad[3 * atom_idx] += delta.x;
            grad[3 * atom_idx + 1] += delta.y;
            grad[3 * atom_idx + 2] += delta.z;
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngleConstraintContribsParams {
    idx1: usize,
    idx2: usize,
    idx3: usize,
    min_angle: f64,
    max_angle: f64,
    force_constant: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AngleConstraintContribs {
    owner: *const ForceField,
    contribs: Vec<AngleConstraintContribsParams>,
}

impl AngleConstraintContribs {
    pub fn new(owner: &ForceField) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: dp_forceField = owner;
        Self {
            owner,
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) {
        // RDKit✔️✔️: RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit✔️✔️: RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        // RDKit✔️✔️: URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️: URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️: PRECONDITION(maxAngleDeg >= minAngleDeg,
        // RDKit✔️✔️:              "minAngleDeg must be <= maxAngleDeg");
        let owner = self.owner();
        assert!((0.0..=180.0).contains(&min_angle_deg));
        assert!((0.0..=180.0).contains(&max_angle_deg));
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(
            max_angle_deg >= min_angle_deg,
            "minAngleDeg must be <= maxAngleDeg"
        );
        // RDKit✔️✔️: d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
        // RDKit✔️✔️:                         forceConst);
        self.contribs.push(AngleConstraintContribsParams {
            idx1,
            idx2,
            idx3,
            min_angle: min_angle_deg,
            max_angle: max_angle_deg,
            force_constant,
        });
    }

    pub fn add_contrib_relative(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        relative: bool,
        mut min_angle_deg: f64,
        mut max_angle_deg: f64,
        force_constant: f64,
    ) {
        // RDKit✔️✔️: const RDGeom::PointPtrVect &pos = dp_forceField->positions();
        // RDKit✔️✔️: URANGE_CHECK(idx1, pos.size());
        // RDKit✔️✔️: URANGE_CHECK(idx2, pos.size());
        // RDKit✔️✔️: URANGE_CHECK(idx3, pos.size());
        // RDKit✔️✔️: PRECONDITION(maxAngleDeg >= minAngleDeg,
        // RDKit✔️✔️:              "minAngleDeg must be <= maxAngleDeg");
        let owner = self.owner();
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(
            max_angle_deg >= min_angle_deg,
            "minAngleDeg must be <= maxAngleDeg"
        );
        // RDKit✔️✔️: if (relative) {
        if relative {
            // RDKit✔️✔️: const RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
            // RDKit✔️✔️: const RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
            // RDKit✔️✔️: const RDGeom::Point3D &p3 = *((RDGeom::Point3D *)pos[idx3]);
            // RDKit✔️✔️: const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
            // RDKit✔️✔️: const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
            // RDKit✔️✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
            // RDKit✔️✔️: double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
            // RDKit✔️✔️: cosTheta = std::clamp(cosTheta, -1.0, 1.0);
            // RDKit✔️✔️: const double angle = RAD2DEG * acos(cosTheta);
            // RDKit✔️✔️: minAngleDeg += angle;
            // RDKit✔️✔️: maxAngleDeg += angle;
            // RDKit✔️✔️: }
            let angle = angle_degrees(
                owner.positions()[idx1],
                owner.positions()[idx2],
                owner.positions()[idx3],
            );
            min_angle_deg += angle;
            max_angle_deg += angle;
        }
        // RDKit✔️✔️: RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit✔️✔️: RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        assert!((0.0..=180.0).contains(&min_angle_deg));
        assert!((0.0..=180.0).contains(&max_angle_deg));
        // RDKit✔️✔️: d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
        // RDKit✔️✔️:                         forceConst);
        self.contribs.push(AngleConstraintContribsParams {
            idx1,
            idx2,
            idx3,
            min_angle: min_angle_deg,
            max_angle: max_angle_deg,
            force_constant,
        });
    }

    #[must_use]
    pub fn compute_angle_term(&self, angle: f64, contrib: &AngleConstraintContribsParams) -> f64 {
        // RDKit✔️✔️: double angleTerm = 0.0;
        let mut angle_term = 0.0;
        // RDKit✔️✔️: if (angle < contrib.minAngle) {
        // RDKit✔️✔️:   angleTerm = angle - contrib.minAngle;
        // RDKit✔️✔️: } else if (angle > contrib.maxAngle) {
        // RDKit✔️✔️:   angleTerm = angle - contrib.maxAngle;
        // RDKit✔️✔️: }
        if angle < contrib.min_angle {
            angle_term = angle - contrib.min_angle;
        } else if angle > contrib.max_angle {
            angle_term = angle - contrib.max_angle;
        }
        // RDKit✔️✔️: return angleTerm;
        angle_term
    }

    #[must_use]
    pub fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the RDKit owner pointer lifetime convention for contribs.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for AngleConstraintContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️: double accum = 0.0;
        let mut accum = 0.0;
        // RDKit✔️✔️: for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            // RDKit✔️✔️: const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx1 + 2]);
            // RDKit✔️✔️: const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx2 + 2]);
            // RDKit✔️✔️: const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx3 + 2]);
            let p1 = point_from_flat(pos, contrib.idx1);
            let p2 = point_from_flat(pos, contrib.idx2);
            let p3 = point_from_flat(pos, contrib.idx3);
            // RDKit✔️✔️: const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
            // RDKit✔️✔️: const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
            // RDKit✔️✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
            // RDKit✔️✔️: double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
            // RDKit✔️✔️: cosTheta = std::clamp(cosTheta, -1.0, 1.0);
            // RDKit✔️✔️: const double angle = RAD2DEG * acos(cosTheta);
            let angle = angle_degrees(p1, p2, p3);
            // RDKit✔️✔️: const double angleTerm = computeAngleTerm(angle, contrib);
            let angle_term = self.compute_angle_term(angle, contrib);
            // RDKit✔️✔️: accum += contrib.forceConstant * angleTerm * angleTerm;
            accum += contrib.force_constant * angle_term * angle_term;
        }
        // RDKit✔️✔️: return accum;
        accum
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️: for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            let p1 = point_from_flat(pos, contrib.idx1);
            let p2 = point_from_flat(pos, contrib.idx2);
            let p3 = point_from_flat(pos, contrib.idx3);
            // RDKit✔️✔️: const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
            let r = [p1 - p2, p3 - p2];
            // RDKit✔️✔️: const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
            // RDKit✔️✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
            let r_length_sq = [r[0].length_sq().max(1.0e-5), r[1].length_sq().max(1.0e-5)];
            // RDKit✔️✔️: double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
            // RDKit✔️✔️: cosTheta = std::clamp(cosTheta, -1.0, 1.0);
            // RDKit✔️✔️: const double angle = RAD2DEG * acos(cosTheta);
            let cos_theta = (r[0].dot_product(r[1]) / (r_length_sq[0] * r_length_sq[1]).sqrt())
                .clamp(-1.0, 1.0);
            let angle = RAD2DEG * cos_theta.acos();
            // RDKit✔️✔️: const double angleTerm = computeAngleTerm(angle, contrib);
            let angle_term = self.compute_angle_term(angle, contrib);
            // RDKit✔️✔️: const double dE_dTheta = 2.0 * RAD2DEG * contrib.forceConstant * angleTerm;
            let d_e_d_theta = 2.0 * RAD2DEG * contrib.force_constant * angle_term;
            // RDKit✔️✔️: const RDGeom::Point3D rp = r[1].crossProduct(r[0]);
            // RDKit✔️✔️: const double prefactor = dE_dTheta / std::max(1.0e-5, rp.length());
            // RDKit✔️✔️: const double t[2] = {-prefactor / rLengthSq[0], prefactor / rLengthSq[1]};
            let rp = r[1].cross_product(r[0]);
            let prefactor = d_e_d_theta / rp.length().max(1.0e-5);
            let t = [-prefactor / r_length_sq[0], prefactor / r_length_sq[1]];
            // RDKit✔️✔️: RDGeom::Point3D dedp[3];
            // RDKit✔️✔️: dedp[0] = r[0].crossProduct(rp) * t[0];
            // RDKit✔️✔️: dedp[2] = r[1].crossProduct(rp) * t[1];
            // RDKit✔️✔️: dedp[1] = -dedp[0] - dedp[2];
            let dedp0 = r[0].cross_product(rp) * t[0];
            let dedp2 = r[1].cross_product(rp) * t[1];
            let dedp = [dedp0, -dedp0 - dedp2, dedp2];
            // RDKit✔️✔️: double *g[3] = {&(grad[3 * contrib.idx1]), &(grad[3 * contrib.idx2]),
            // RDKit✔️✔️:                 &(grad[3 * contrib.idx3])};
            // RDKit✔️✔️: for (unsigned int i = 0; i < 3; ++i) {
            // RDKit✔️✔️:   g[i][0] += dedp[i].x;
            // RDKit✔️✔️:   g[i][1] += dedp[i].y;
            // RDKit✔️✔️:   g[i][2] += dedp[i].z;
            // RDKit✔️✔️: }
            for (atom_idx, delta) in [
                (contrib.idx1, dedp[0]),
                (contrib.idx2, dedp[1]),
                (contrib.idx3, dedp[2]),
            ] {
                grad[3 * atom_idx] += delta.x;
                grad[3 * atom_idx + 1] += delta.y;
                grad[3 * atom_idx + 2] += delta.z;
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TorsionConstraintContrib {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    at3_idx: usize,
    at4_idx: usize,
    min_dihedral_deg: f64,
    max_dihedral_deg: f64,
    force_constant: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PositionConstraintContrib {
    owner: *const ForceField,
    at_idx: usize,
    max_displ: f64,
    pos0: ForceFieldVec3,
    force_constant: f64,
}

impl PositionConstraintContrib {
    pub fn new(owner: &ForceField, idx: usize, max_displ: f64, force_constant: f64) -> Self {
        // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️: const RDGeom::PointPtrVect &pos = owner->positions();
        // RDKit✔️✔️: URANGE_CHECK(idx, pos.size());
        assert!(idx < owner.positions().len());
        // RDKit✔️✔️: dp_forceField = owner;
        // RDKit✔️✔️: d_atIdx = idx;
        // RDKit✔️✔️: d_maxDispl = maxDispl;
        // RDKit✔️✔️: d_pos0 = *((RDGeom::Point3D *)pos[idx]);
        // RDKit✔️✔️: d_forceConstant = forceConst;
        Self {
            owner,
            at_idx: idx,
            max_displ,
            pos0: owner.positions()[idx],
            force_constant,
        }
    }

    pub fn new_relative(
        owner: &ForceField,
        idx: usize,
        relative: bool,
        max_displ: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: PositionConstraintContrib(ForceField *owner, unsigned int idx,
        // RDKit✔️✔️:                           double maxDispl, double forceConst);
        //
        // RDKit has no relative PositionConstraintContrib overload; the Python
        // wrapper also forwards only idx, maxDispl, and forceConst. COSMolKit keeps
        // this explicit compatibility entry point for plan symmetry, but relative
        // has no source behavior to apply.
        let _ = relative;
        Self::new(owner, idx, max_displ, force_constant)
    }

    #[must_use]
    pub fn max_displ(&self) -> f64 {
        self.max_displ
    }

    #[must_use]
    pub fn pos0(&self) -> ForceFieldVec3 {
        self.pos0
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️: RDGeom::Point3D p(pos[3 * d_atIdx], pos[3 * d_atIdx + 1],
        // RDKit✔️✔️:                   pos[3 * d_atIdx + 2]);
        let p = point_from_flat(pos, self.at_idx);
        // RDKit✔️✔️: double dist = (p - d_pos0).length();
        let dist = (p - self.pos0).length();
        // RDKit✔️✔️: double distTerm = std::max(dist - d_maxDispl, 0.0);
        let dist_term = (dist - self.max_displ).max(0.0);
        // RDKit✔️✔️: double res = 0.5 * d_forceConstant * distTerm * distTerm;
        // RDKit✔️✔️: return res;
        0.5 * self.force_constant * dist_term * dist_term
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️: RDGeom::Point3D p(pos[3 * d_atIdx], pos[3 * d_atIdx + 1],
        // RDKit✔️✔️:                   pos[3 * d_atIdx + 2]);
        let p = point_from_flat(pos, self.at_idx);
        // RDKit✔️✔️: double dist = (p - d_pos0).length();
        let dist = (p - self.pos0).length();
        // RDKit✔️✔️: double preFactor = 0.0;
        // RDKit✔️✔️: if (dist > d_maxDispl) {
        // RDKit✔️✔️:   preFactor = dist - d_maxDispl;
        // RDKit✔️✔️: } else {
        // RDKit✔️✔️:   return;
        // RDKit✔️✔️: }
        let mut pre_factor = if dist > self.max_displ {
            dist - self.max_displ
        } else {
            return;
        };
        // RDKit✔️✔️: preFactor *= d_forceConstant;
        pre_factor *= self.force_constant;
        // RDKit✔️✔️: for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:   double dGrad = preFactor * (p[i] - d_pos0[i]) / std::max(dist, 1.0e-8);
        // RDKit✔️✔️:   grad[3 * d_atIdx + i] += dGrad;
        // RDKit✔️✔️: }
        let offset = 3 * self.at_idx;
        for i in 0..3 {
            let d_grad = pre_factor * (p[i] - self.pos0[i]) / dist.max(1.0e-8);
            grad[offset + i] += d_grad;
        }
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the RDKit owner pointer lifetime convention for contribs.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for PositionConstraintContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        Self::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        Self::get_grad(self, pos, grad);
    }
}

impl TorsionConstraintContrib {
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        min_dihedral_deg: f64,
        max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: checkPrecondition(owner, idx1, idx2, idx3, idx4, minDihedralDeg,
        // RDKit✔️✔️:                   maxDihedralDeg);
        check_torsion_precondition(
            owner,
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
        );
        // RDKit✔️✔️: setParameters(owner, idx1, idx2, idx3, idx4, minDihedralDeg, maxDihedralDeg,
        // RDKit✔️✔️:               forceConst);
        let mut contrib = Self {
            owner: std::ptr::null(),
            at1_idx: 0,
            at2_idx: 0,
            at3_idx: 0,
            at4_idx: 0,
            min_dihedral_deg: 0.0,
            max_dihedral_deg: 0.0,
            force_constant: 0.0,
        };
        contrib.set_parameters(
            owner,
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        );
        contrib
    }

    pub fn new_relative(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        relative: bool,
        mut min_dihedral_deg: f64,
        mut max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Self {
        // RDKit✔️✔️: checkPrecondition(owner, idx1, idx2, idx3, idx4, minDihedralDeg,
        // RDKit✔️✔️:                   maxDihedralDeg);
        check_torsion_precondition(
            owner,
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
        );
        // RDKit✔️✔️: if (relative) {
        if relative {
            // RDKit✔️✔️: double dihedral;
            // RDKit✔️✔️: RDKit::ForceFieldsHelper::computeDihedral(owner->positions(), idx1, idx2,
            // RDKit✔️✔️:                                               idx3, idx4, &dihedral);
            // RDKit✔️✔️: dihedral *= RAD2DEG;
            // RDKit✔️✔️: minDihedralDeg += dihedral;
            // RDKit✔️✔️: maxDihedralDeg += dihedral;
            // RDKit✔️✔️: }
            let dihedral =
                compute_dihedral_from_position_vec(owner.positions(), idx1, idx2, idx3, idx4, true)
                    .dihedral
                    .expect("requested dihedral");
            min_dihedral_deg += dihedral * RAD2DEG;
            max_dihedral_deg += dihedral * RAD2DEG;
        }
        // RDKit✔️✔️: setParameters(owner, idx1, idx2, idx3, idx4, minDihedralDeg, maxDihedralDeg,
        // RDKit✔️✔️:               forceConst);
        let mut contrib = Self {
            owner: std::ptr::null(),
            at1_idx: 0,
            at2_idx: 0,
            at3_idx: 0,
            at4_idx: 0,
            min_dihedral_deg: 0.0,
            max_dihedral_deg: 0.0,
            force_constant: 0.0,
        };
        contrib.set_parameters(
            owner,
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        );
        contrib
    }

    fn set_parameters(
        &mut self,
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        mut min_dihedral_deg: f64,
        mut max_dihedral_deg: f64,
        force_constant: f64,
    ) {
        // RDKit✔️✔️: dp_forceField = owner;
        // RDKit✔️✔️: d_at1Idx = idx1;
        // RDKit✔️✔️: d_at2Idx = idx2;
        // RDKit✔️✔️: d_at3Idx = idx3;
        // RDKit✔️✔️: d_at4Idx = idx4;
        self.owner = owner;
        self.at1_idx = idx1;
        self.at2_idx = idx2;
        self.at3_idx = idx3;
        self.at4_idx = idx4;
        // RDKit✔️✔️: RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralDeg);
        // RDKit✔️✔️: RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralDeg);
        normalize_angle_deg(&mut min_dihedral_deg);
        normalize_angle_deg(&mut max_dihedral_deg);
        // RDKit✔️✔️: d_minDihedralDeg = minDihedralDeg;
        // RDKit✔️✔️: d_maxDihedralDeg = maxDihedralDeg;
        // RDKit✔️✔️: d_forceConstant = forceConst;
        self.min_dihedral_deg = min_dihedral_deg;
        self.max_dihedral_deg = max_dihedral_deg;
        self.force_constant = force_constant;
    }

    #[must_use]
    pub fn min_dihedral_deg(&self) -> f64 {
        self.min_dihedral_deg
    }

    #[must_use]
    pub fn max_dihedral_deg(&self) -> f64 {
        self.max_dihedral_deg
    }

    #[must_use]
    pub fn compute_dihedral_term(&self, dihedral: f64) -> f64 {
        // RDKit✔️✔️: double dihedralTarget = dihedral;
        let mut dihedral_target = dihedral;
        // RDKit✔️✔️: if (!(dihedral > d_minDihedralDeg && dihedral < d_maxDihedralDeg) &&
        // RDKit✔️✔️:     !(dihedral > d_minDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg) &&
        // RDKit✔️✔️:     !(dihedral < d_maxDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg)) {
        if !(dihedral > self.min_dihedral_deg && dihedral < self.max_dihedral_deg)
            && !(dihedral > self.min_dihedral_deg && self.min_dihedral_deg > self.max_dihedral_deg)
            && !(dihedral < self.max_dihedral_deg && self.min_dihedral_deg > self.max_dihedral_deg)
        {
            // RDKit✔️✔️: double dihedralMinTarget = dihedral - d_minDihedralDeg;
            // RDKit✔️✔️: RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralMinTarget);
            let mut dihedral_min_target = dihedral - self.min_dihedral_deg;
            normalize_angle_deg(&mut dihedral_min_target);
            // RDKit✔️✔️: double dihedralMaxTarget = dihedral - d_maxDihedralDeg;
            // RDKit✔️✔️: RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralMaxTarget);
            let mut dihedral_max_target = dihedral - self.max_dihedral_deg;
            normalize_angle_deg(&mut dihedral_max_target);
            // RDKit✔️✔️: if (fabs(dihedralMinTarget) < fabs(dihedralMaxTarget)) {
            // RDKit✔️✔️:   dihedralTarget = d_minDihedralDeg;
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   dihedralTarget = d_maxDihedralDeg;
            // RDKit✔️✔️: }
            if dihedral_min_target.abs() < dihedral_max_target.abs() {
                dihedral_target = self.min_dihedral_deg;
            } else {
                dihedral_target = self.max_dihedral_deg;
            }
        }
        // RDKit✔️✔️: double dihedralTerm = dihedral - dihedralTarget;
        // RDKit✔️✔️: RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralTerm);
        // RDKit✔️✔️: return dihedralTerm;
        let mut dihedral_term = dihedral - dihedral_target;
        normalize_angle_deg(&mut dihedral_term);
        dihedral_term
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the RDKit owner pointer lifetime convention for contribs.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for TorsionConstraintContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        let _ = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️: double dihedral;
        // RDKit✔️✔️: RDKit::ForceFieldsHelper::computeDihedral(pos, d_at1Idx, d_at2Idx, d_at3Idx,
        // RDKit✔️✔️:                                             d_at4Idx, &dihedral);
        // RDKit✔️✔️: dihedral *= RAD2DEG;
        let dihedral = compute_dihedral_from_flat(
            pos,
            self.at1_idx,
            self.at2_idx,
            self.at3_idx,
            self.at4_idx,
            true,
        )
        .dihedral
        .expect("requested dihedral")
            * RAD2DEG;
        // RDKit✔️✔️: double dihedralTerm = computeDihedralTerm(dihedral);
        let dihedral_term = self.compute_dihedral_term(dihedral);
        // RDKit✔️✔️: double res = d_forceConstant * dihedralTerm * dihedralTerm;
        // RDKit✔️✔️: return res;
        self.force_constant * dihedral_term * dihedral_term
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
        let owner = self.owner();
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️: RDGeom::Point3D r[4];
        // RDKit✔️✔️: RDGeom::Point3D t[2];
        // RDKit✔️✔️: double d[2];
        // RDKit✔️✔️: double dihedral;
        // RDKit✔️✔️: RDKit::ForceFieldsHelper::computeDihedral(
        // RDKit✔️✔️:     pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, &dihedral, nullptr, r, t, d);
        let dihedral_output = compute_dihedral_from_flat(
            pos,
            self.at1_idx,
            self.at2_idx,
            self.at3_idx,
            self.at4_idx,
            true,
        );
        let r = dihedral_output.r;
        // RDKit✔️✔️: dihedral *= RAD2DEG;
        let dihedral = dihedral_output.dihedral.expect("requested dihedral") * RAD2DEG;
        // RDKit✔️✔️: double dihedralTerm = computeDihedralTerm(dihedral);
        let dihedral_term = self.compute_dihedral_term(dihedral);
        // RDKit✔️✔️: double dE_dPhi = 2.0 * RAD2DEG * d_forceConstant * dihedralTerm;
        let d_e_d_phi = 2.0 * RAD2DEG * self.force_constant * dihedral_term;
        // RDKit✔️✔️: double d23 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
        let d23 = owner.distance(self.at2_idx, self.at3_idx, Some(pos));
        if d23 == 0.0 {
            return;
        }
        // RDKit✔️✔️: RDGeom::Point3D r31(pos[3 * d_at3Idx] - pos[3 * d_at1Idx],
        // RDKit✔️✔️:                     pos[3 * d_at3Idx + 1] - pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                     pos[3 * d_at3Idx + 2] - pos[3 * d_at1Idx + 2]);
        let r31 = point_from_flat(pos, self.at3_idx) - point_from_flat(pos, self.at1_idx);
        // RDKit✔️✔️: RDGeom::Point3D r42(pos[3 * d_at4Idx] - pos[3 * d_at2Idx],
        // RDKit✔️✔️:                     pos[3 * d_at4Idx + 1] - pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                     pos[3 * d_at4Idx + 2] - pos[3 * d_at2Idx + 2]);
        let r42 = point_from_flat(pos, self.at4_idx) - point_from_flat(pos, self.at2_idx);
        // RDKit✔️✔️: double prefactor = dE_dPhi / d23;
        let prefactor = d_e_d_phi / d23;
        // RDKit✔️✔️: RDGeom::Point3D tt[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
        let tt = [r[0].cross_product(r[1]), r[2].cross_product(r[3])];
        if tt[0].length_sq() == 0.0 || tt[1].length_sq() == 0.0 {
            return;
        }
        // RDKit✔️✔️: RDGeom::Point3D dedt[2] = {
        // RDKit✔️✔️:     tt[0].crossProduct(r[2]) / tt[0].lengthSq() * prefactor,
        // RDKit✔️✔️:     tt[1].crossProduct(r[1]) / tt[1].lengthSq() * prefactor};
        let dedt = [
            tt[0].cross_product(r[2]) / tt[0].length_sq() * prefactor,
            tt[1].cross_product(r[1]) / tt[1].length_sq() * prefactor,
        ];
        // RDKit✔️✔️: RDGeom::Point3D dedp[4] = {
        // RDKit✔️✔️:     r[2].crossProduct(dedt[0]),
        // RDKit✔️✔️:     r31.crossProduct(dedt[0]) - r[3].crossProduct(dedt[1]),
        // RDKit✔️✔️:     r[0].crossProduct(dedt[0]) + r42.crossProduct(dedt[1]),
        // RDKit✔️✔️:     r[2].crossProduct(dedt[1])};
        let dedp = [
            r[2].cross_product(dedt[0]),
            r31.cross_product(dedt[0]) - r[3].cross_product(dedt[1]),
            r[0].cross_product(dedt[0]) + r42.cross_product(dedt[1]),
            r[2].cross_product(dedt[1]),
        ];
        // RDKit✔️✔️: for (unsigned int i = 0; i < 4; ++i) {
        // RDKit✔️✔️:   g[i][0] += dedp[i].x;
        // RDKit✔️✔️:   g[i][1] += dedp[i].y;
        // RDKit✔️✔️:   g[i][2] += dedp[i].z;
        // RDKit✔️✔️: }
        for (atom_idx, delta) in [
            (self.at1_idx, dedp[0]),
            (self.at2_idx, dedp[1]),
            (self.at3_idx, dedp[2]),
            (self.at4_idx, dedp[3]),
        ] {
            grad[3 * atom_idx] += delta.x;
            grad[3 * atom_idx + 1] += delta.y;
            grad[3 * atom_idx + 2] += delta.z;
        }
    }
}

fn check_torsion_precondition(
    owner: &ForceField,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    min_dihedral_deg: f64,
    max_dihedral_deg: f64,
) {
    // RDKit✔️✔️: PRECONDITION(owner, "bad owner");
    // RDKit✔️✔️: PRECONDITION(!(minDihedralDeg > maxDihedralDeg),
    // RDKit✔️✔️:              "minDihedralDeg must be <= maxDihedralDeg");
    assert!(
        min_dihedral_deg <= max_dihedral_deg,
        "minDihedralDeg must be <= maxDihedralDeg"
    );
    // RDKit✔️✔️: URANGE_CHECK(idx1, owner->positions().size());
    // RDKit✔️✔️: URANGE_CHECK(idx2, owner->positions().size());
    // RDKit✔️✔️: URANGE_CHECK(idx3, owner->positions().size());
    // RDKit✔️✔️: URANGE_CHECK(idx4, owner->positions().size());
    assert!(idx1 < owner.positions().len());
    assert!(idx2 < owner.positions().len());
    assert!(idx3 < owner.positions().len());
    assert!(idx4 < owner.positions().len());
}

fn calc_gradient_wrapper(ff: &ForceField, pos: &[f64], grad: &mut [f64]) -> f64 {
    // RDKit✔️✔️: for (unsigned int i = 0;
    // RDKit✔️✔️:      i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
    // RDKit✔️✔️:   grad[i] = 0.0;
    // RDKit✔️✔️: }
    for value in grad.iter_mut() {
        *value = 0.0;
    }
    // RDKit✔️✔️: mp_ffHolder->calcGrad(pos, grad);
    ff.calc_grad(pos, grad);
    // RDKit✔️✔️: double maxGrad = -1e8;
    // RDKit✔️✔️: double gradScale = 0.1;
    let mut max_grad = -1.0e8_f64;
    let mut grad_scale = 0.1_f64;
    // RDKit✔️✔️: grad[i] *= gradScale;
    // RDKit✔️✔️: if (fabs(grad[i]) > maxGrad) {
    // RDKit✔️✔️:   maxGrad = fabs(grad[i]);
    // RDKit✔️✔️: }
    for value in grad.iter_mut() {
        *value *= grad_scale;
        if value.abs() > max_grad {
            max_grad = value.abs();
        }
    }
    // RDKit✔️✔️: if (maxGrad > 10.0) {
    // RDKit✔️✔️:   while (maxGrad * gradScale > 10.0) {
    // RDKit✔️✔️:     gradScale *= .5;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0;
    // RDKit✔️✔️:        i < mp_ffHolder->numPoints() * mp_ffHolder->dimension(); i++) {
    // RDKit✔️✔️:     grad[i] *= gradScale;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if max_grad > 10.0 {
        while max_grad * grad_scale > 10.0 {
            grad_scale *= 0.5;
        }
        for value in grad.iter_mut() {
            *value *= grad_scale;
        }
    }
    // RDKit✔️✔️: res = gradScale;
    // RDKit✔️✔️: return res;
    grad_scale
}

fn row1_bfgs_trace_enabled(dim: usize) -> bool {
    dim == 18 && std::env::var("RDKIT_ROW1_TRACE").as_deref() == Ok("1")
}

fn bfgs_minimize<Energy, Gradient>(
    mut pos: Vec<f64>,
    grad_tol: f64,
    mut func: Energy,
    mut grad_func: Gradient,
    snapshot_freq: usize,
    mut snapshot_vect: Option<&mut Vec<ForceFieldSnapshot>>,
    _func_tol: f64,
    max_its: usize,
) -> (i32, Vec<f64>)
where
    Energy: FnMut(&[f64]) -> f64,
    Gradient: FnMut(&[f64], &mut [f64]) -> f64,
{
    // RDKit✔️✔️: PRECONDITION(pos, "bad input array");
    // RDKit✔️✔️: PRECONDITION(gradTol > 0, "bad tolerance");
    assert!(!pos.is_empty(), "bad input array");
    assert!(grad_tol > 0.0, "bad tolerance");
    let dim = pos.len();
    // RDKit✔️✔️: std::vector<double> grad(dim);
    // RDKit✔️✔️: std::vector<double> dGrad(dim);
    // RDKit✔️✔️: std::vector<double> hessDGrad(dim);
    // RDKit✔️✔️: std::vector<double> xi(dim);
    // RDKit✔️✔️: std::vector<double> invHessian(dim * dim, 0);
    let mut grad = vec![0.0; dim];
    let mut d_grad = vec![0.0; dim];
    let mut hess_d_grad = vec![0.0; dim];
    let mut xi = vec![0.0; dim];
    let mut inv_hessian = vec![0.0; dim * dim];
    let mut new_pos = vec![0.0; dim];
    let snapshot_freq = snapshot_freq.min(max_its);
    // RDKit✔️✔️: double fp = func(pos);
    // RDKit✔️✔️: gradFunc(pos, grad.data());
    let mut fp = func(&pos);
    grad_func(&pos, &mut grad);
    // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:   unsigned int itab = i * dim;
    // RDKit✔️✔️:   invHessian[itab + i] = 1.0;
    // RDKit✔️✔️:   xi[i] = -grad[i];
    // RDKit✔️✔️:   sum += pos[i] * pos[i];
    // RDKit✔️✔️: }
    let mut sum = 0.0;
    for i in 0..dim {
        inv_hessian[i * dim + i] = 1.0;
        xi[i] = -grad[i];
        sum += pos[i] * pos[i];
    }
    // RDKit✔️✔️: double maxStep = MAXSTEP * std::max(sqrt(sum), static_cast<double>(dim));
    let sqrt_sum = sum.sqrt();
    let max_step = if sqrt_sum > dim as f64 {
        MAXSTEP * sqrt_sum
    } else {
        MAXSTEP * dim as f64
    };
    for iter in 1..=max_its {
        if row1_bfgs_trace_enabled(dim) && iter <= 12 {
            println!("[row1-dir] iter={iter} values={xi:?}");
        }
        // RDKit✔️✔️: linearSearch(dim, pos, fp, grad.data(), xi.data(), newPos.get(), funcVal, func,
        // RDKit✔️✔️:              maxStep, status);
        let (status, new_val) = linear_search(
            iter,
            &pos,
            fp,
            &grad,
            &mut xi,
            &mut new_pos,
            &mut func,
            max_step,
        );
        // RDKit✔️✔️: CHECK_INVARIANT(status >= 0, "bad direction in linearSearch");
        assert!(status >= 0, "bad direction in linearSearch");
        let func_val = new_val;
        // RDKit✔️✔️: fp = funcVal;
        fp = func_val;
        // RDKit✔️✔️: xi[i] = newPos[i] - pos[i];
        // RDKit✔️✔️: pos[i] = newPos[i];
        // RDKit✔️✔️: double temp = fabs(xi[i]) / std::max(fabs(pos[i]), 1.0);
        // RDKit✔️✔️: dGrad[i] = grad[i];
        let mut test = 0.0_f64;
        for i in 0..dim {
            xi[i] = new_pos[i] - pos[i];
            pos[i] = new_pos[i];
            let abs_pos = pos[i].abs();
            let denom = if abs_pos > 1.0 { abs_pos } else { 1.0 };
            let temp = xi[i].abs() / denom;
            if temp > test {
                test = temp;
            }
            d_grad[i] = grad[i];
        }
        // RDKit✔️✔️: if (test < TOLX) {
        if test < TOLX {
            if let Some(snapshot_vect) = &mut snapshot_vect
                && snapshot_freq != 0
            {
                snapshot_vect.push(ForceFieldSnapshot {
                    positions: new_pos.clone(),
                    energy: fp,
                });
            }
            return (0, pos);
        }
        // RDKit✔️✔️: double gradScale = gradFunc(pos, grad.data());
        let grad_scale = grad_func(&pos, &mut grad);
        if row1_bfgs_trace_enabled(dim) && iter <= 12 {
            println!("[row1-grad] iter={iter} values={grad:?}");
        }
        // RDKit✔️✔️: double term = std::max(funcVal * gradScale, 1.0);
        let func_term = func_val * grad_scale;
        let term = if func_term > 1.0 { func_term } else { 1.0 };
        test = 0.0;
        // RDKit✔️✔️: double temp = fabs(grad[i]) * std::max(fabs(pos[i]), 1.0);
        // RDKit✔️✔️: test = std::max(test, temp);
        // RDKit✔️✔️: dGrad[i] = grad[i] - dGrad[i];
        for i in 0..dim {
            let abs_pos = pos[i].abs();
            let scale = if abs_pos > 1.0 { abs_pos } else { 1.0 };
            let temp = grad[i].abs() * scale;
            if temp > test {
                test = temp;
            }
            d_grad[i] = grad[i] - d_grad[i];
        }
        if row1_bfgs_trace_enabled(dim) && iter <= 12 {
            println!("[row1-dgrad] iter={iter} values={d_grad:?}");
        }
        // RDKit✔️✔️: test /= term;
        test /= term;
        // RDKit✔️✔️: if (test < gradTol) {
        if test < grad_tol {
            if let Some(snapshot_vect) = &mut snapshot_vect
                && snapshot_freq != 0
            {
                snapshot_vect.push(ForceFieldSnapshot {
                    positions: new_pos.clone(),
                    energy: fp,
                });
            }
            return (0, pos);
        }

        if row1_bfgs_trace_enabled(dim) && iter <= 80 {
            println!(
                "[row1-bfgs] iter={iter} status={status} line_test={test:.17} gradScale={grad_scale:.17} fp={fp:.17}"
            );
        }
        // RDKit✔️✔️: double fac = 0, fae = 0, sumDGrad = 0, sumXi = 0;
        let mut fac = 0.0;
        let mut fae = 0.0;
        let mut sum_d_grad = 0.0;
        let mut sum_xi = 0.0;
        // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
        // RDKit✔️✔️:   hdgradi = 0.0;
        // RDKit✔️✔️:   for (unsigned int j = 0; j < dim; ++j, ++ivh, ++dgj) {
        // RDKit✔️✔️:     hdgradi += *ivh * *dgj;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   fac += dGrad[i] * xi[i];
        // RDKit✔️✔️:   fae += dGrad[i] * hessDGrad[i];
        // RDKit✔️✔️:   sumDGrad += dGrad[i] * dGrad[i];
        // RDKit✔️✔️:   sumXi += xi[i] * xi[i];
        // RDKit✔️✔️: }
        for i in 0..dim {
            hess_d_grad[i] = 0.0;
            for j in 0..dim {
                hess_d_grad[i] += inv_hessian[i * dim + j] * d_grad[j];
            }
            fac += d_grad[i] * xi[i];
            fae += d_grad[i] * hess_d_grad[i];
            sum_d_grad += d_grad[i] * d_grad[i];
            sum_xi += xi[i] * xi[i];
        }
        if row1_bfgs_trace_enabled(dim) && iter <= 12 {
            println!("[row1-hdgrad] iter={iter} values={hess_d_grad:?}");
            let invh_head = [
                inv_hessian[0],
                inv_hessian[1],
                inv_hessian[2],
                inv_hessian[18],
                inv_hessian[19],
                inv_hessian[20],
            ];
            println!("[row1-invh] iter={iter} values={invh_head:?}");
            if iter <= 2 {
                println!(
                    "[row1-invh-row0] iter={iter} values={:?}",
                    &inv_hessian[0..dim]
                );
                println!(
                    "[row1-invh-row1] iter={iter} values={:?}",
                    &inv_hessian[dim..(2 * dim)]
                );
                if iter == 1 {
                    for row in 0..dim {
                        let start = row * dim;
                        let end = start + dim;
                        println!(
                            "[row1-invh-row{row}] iter={iter} values={:?}",
                            &inv_hessian[start..end]
                        );
                    }
                }
            }
        }
        // RDKit✔️✔️: if (fac > sqrt(EPS * sumDGrad * sumXi)) {
        if fac > (EPS * sum_d_grad * sum_xi).sqrt() {
            // RDKit✔️✔️: fac = 1.0 / fac;
            // RDKit✔️✔️: double fad = 1.0 / fae;
            fac = 1.0 / fac;
            let fad = 1.0 / fae;
            // RDKit✔️✔️: dGrad[i] = fac * xi[i] - fad * hessDGrad[i];
            for i in 0..dim {
                d_grad[i] = fac * xi[i] - fad * hess_d_grad[i];
            }
            // RDKit✔️✔️: invHessian[itab + j] += pxi * *pxj - hdgi * *hdgj + dgi * *dgj;
            // RDKit✔️✔️: invHessian[j * dim + i] = invHessian[itab + j];
            for i in 0..dim {
                let pxi = fac * xi[i];
                let hdgi = fad * hess_d_grad[i];
                let dgi = fae * d_grad[i];
                for j in i..dim {
                    let term1 = pxi * xi[j];
                    let term2 = hdgi * hess_d_grad[j];
                    let term3 = dgi * d_grad[j];
                    if row1_bfgs_trace_enabled(dim) && iter == 1 && i == 13 && j == 13 {
                        println!(
                            "[row1-invh-update-bits] iter={iter} row=13 col=13 old_bits={:#018x} pxi_bits={:#018x} xj_bits={:#018x} hdgi_bits={:#018x} hdgj_bits={:#018x} dgi_bits={:#018x} dgj_bits={:#018x} term1_bits={:#018x} term2_bits={:#018x} term3_bits={:#018x}",
                            inv_hessian[i * dim + j].to_bits(),
                            pxi.to_bits(),
                            xi[j].to_bits(),
                            hdgi.to_bits(),
                            hess_d_grad[j].to_bits(),
                            dgi.to_bits(),
                            d_grad[j].to_bits(),
                            term1.to_bits(),
                            term2.to_bits(),
                            term3.to_bits()
                        );
                    }
                    inv_hessian[i * dim + j] +=
                        pxi * xi[j] - hdgi * hess_d_grad[j] + dgi * d_grad[j];
                    let updated = inv_hessian[i * dim + j];
                    if row1_bfgs_trace_enabled(dim) && iter == 1 && i == 13 && j == 13 {
                        println!(
                            "[row1-invh-update-newbits] iter={iter} row=13 col=13 new_bits={:#018x}",
                            updated.to_bits()
                        );
                    }
                    inv_hessian[j * dim + i] = updated;
                }
            }
        }
        if row1_bfgs_trace_enabled(dim) && iter <= 80 {
            println!(
                "[row1-bfgs-update] iter={iter} fac={fac:.17} fae={fae:.17} sumDGrad={sum_d_grad:.17} sumXi={sum_xi:.17}"
            );
        }
        // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
        // RDKit✔️✔️:   xi[i] = 0.0;
        // RDKit✔️✔️:   for (unsigned int j = 0; j < dim; ++j, ++ivh, ++gj) {
        // RDKit✔️✔️:     pxi -= *ivh * *gj;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        for i in 0..dim {
            xi[i] = 0.0;
            for j in 0..dim {
                xi[i] -= inv_hessian[i * dim + j] * grad[j];
            }
        }
        if row1_bfgs_trace_enabled(dim) && iter <= 2 {
            println!("[row1-nextdir] iter={iter} values={xi:?}");
            if iter <= 2 {
                let mut accum = 0.0_f64;
                let trace_row = if iter == 1 { 13 } else { 0 };
                for j in 0..dim {
                    let term = inv_hessian[trace_row * dim + j] * grad[j];
                    accum -= term;
                    println!(
                        "[row1-nextdir-accum] iter={iter} row={trace_row} col={j} invh_bits={:#018x} grad_bits={:#018x} term_bits={:#018x} accum_bits={:#018x}",
                        inv_hessian[trace_row * dim + j].to_bits(),
                        grad[j].to_bits(),
                        term.to_bits(),
                        accum.to_bits()
                    );
                }
            }
        }
        // RDKit✔️✔️: if (snapshotVect && snapshotFreq && !(iter % snapshotFreq)) {
        if let Some(snapshot_vect) = &mut snapshot_vect
            && snapshot_freq != 0
            && iter % snapshot_freq == 0
        {
            snapshot_vect.push(ForceFieldSnapshot {
                positions: new_pos.clone(),
                energy: fp,
            });
        }
    }
    // RDKit✔️✔️: return 1;
    (1, pos)
}

fn linear_search<Energy>(
    outer_iter: usize,
    old_pt: &[f64],
    old_val: f64,
    grad: &[f64],
    dir: &mut [f64],
    new_pt: &mut [f64],
    func: &mut Energy,
    max_step: f64,
) -> (i32, f64)
where
    Energy: FnMut(&[f64]) -> f64,
{
    // RDKit✔️✔️: const unsigned int MAX_ITER_LINEAR_SEARCH = 1000;
    const MAX_ITER_LINEAR_SEARCH: usize = 1000;
    // RDKit✔️✔️: resCode = -1;
    let mut res_code = -1;
    let mut new_val = old_val;
    // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:   sum += dir[i] * dir[i];
    // RDKit✔️✔️: }
    // RDKit✔️✔️: sum = sqrt(sum);
    let mut sum = 0.0;
    for &value in dir.iter() {
        sum += value * value;
    }
    sum = sum.sqrt();
    // RDKit✔️✔️: if (sum > maxStep) {
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     dir[i] *= maxStep / sum;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if sum > max_step {
        for value in dir.iter_mut() {
            *value *= max_step / sum;
        }
    }
    // RDKit✔️✔️: slope = 0.0;
    // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:   slope += dir[i] * grad[i];
    // RDKit✔️✔️: }
    let mut slope = 0.0;
    for i in 0..dir.len() {
        slope += dir[i] * grad[i];
    }
    // RDKit✔️✔️: if (slope >= 0.0) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    if slope >= 0.0 {
        return (res_code, new_val);
    }
    // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:   double temp = fabs(dir[i]) / std::max(fabs(oldPt[i]), 1.0);
    // RDKit✔️✔️:   if (temp > test) {
    // RDKit✔️✔️:     test = temp;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut test = 0.0_f64;
    for i in 0..dir.len() {
        let abs_old = old_pt[i].abs();
        let denom = if abs_old > 1.0 { abs_old } else { 1.0 };
        let temp = dir[i].abs() / denom;
        if temp > test {
            test = temp;
        }
    }
    // RDKit✔️✔️: lambdaMin = MOVETOL / test;
    // RDKit✔️✔️: lambda = 1.0;
    let lambda_min = MOVETOL / test;
    let mut lambda = 1.0;
    let mut lambda2 = 0.0;
    let mut val2 = 0.0;
    let mut it = 0;
    // RDKit✔️✔️: while (it < MAX_ITER_LINEAR_SEARCH) {
    while it < MAX_ITER_LINEAR_SEARCH {
        // RDKit✔️✔️: if (lambda < lambdaMin) {
        // RDKit✔️✔️:   resCode = 1;
        // RDKit✔️✔️:   break;
        // RDKit✔️✔️: }
        if lambda < lambda_min {
            res_code = 1;
            break;
        }
        // RDKit✔️✔️: newPt[i] = oldPt[i] + lambda * dir[i];
        for i in 0..dir.len() {
            new_pt[i] = old_pt[i] + lambda * dir[i];
        }
        if row1_bfgs_trace_enabled(old_pt.len()) && outer_iter == 2 && it == 1 {
            println!("[row1-linesearch-point] iter={outer_iter} inner_it={it} values={new_pt:?}");
            let bits: Vec<String> = new_pt
                .iter()
                .map(|value| format!("{:#018x}", value.to_bits()))
                .collect();
            println!(
                "[row1-linesearch-point-bits] iter={outer_iter} inner_it={it} values=[{}]",
                bits.join(",")
            );
            for i in 0..old_pt.len() {
                let scaled_dir = lambda * dir[i];
                println!(
                    "[row1-linesearch-operands] iter={outer_iter} inner_it={it} idx={i} old_bits={:#018x} dir_bits={:#018x} lambda_bits={:#018x} scaled_bits={:#018x} new_bits={:#018x}",
                    old_pt[i].to_bits(),
                    dir[i].to_bits(),
                    lambda.to_bits(),
                    scaled_dir.to_bits(),
                    new_pt[i].to_bits(),
                );
            }
        }
        // RDKit✔️✔️: newVal = func(newPt);
        new_val = func(new_pt);
        // RDKit✔️✔️: if (newVal - oldVal <= FUNCTOL * lambda * slope) {
        if new_val - old_val <= FUNCTOL * lambda * slope {
            if row1_bfgs_trace_enabled(old_pt.len()) && outer_iter <= 80 {
                println!(
                    "[row1-linesearch] iter={outer_iter} inner_it={it} accept=1 lambda={lambda:.17} lambdaMin={lambda_min:.17} slope={slope:.17} oldVal={old_val:.17} newVal={new_val:.17}"
                );
            }
            // RDKit✔️✔️: resCode = 0;
            // RDKit✔️✔️: return;
            return (0, new_val);
        }
        let tmp_lambda = if it == 0 {
            // RDKit✔️✔️: tmpLambda = -slope / (2.0 * (newVal - oldVal - slope));
            let tmp_lambda = -slope / (2.0 * (new_val - old_val - slope));
            if row1_bfgs_trace_enabled(old_pt.len()) && outer_iter <= 80 {
                println!(
                    "[row1-linesearch] iter={outer_iter} inner_it={it} accept=0 lambda={lambda:.17} lambdaMin={lambda_min:.17} slope={slope:.17} oldVal={old_val:.17} newVal={new_val:.17} tmpLambda={tmp_lambda:.17} branch=first"
                );
            }
            tmp_lambda
        } else {
            // RDKit✔️✔️: double rhs1 = newVal - oldVal - lambda * slope;
            // RDKit✔️✔️: double rhs2 = val2 - oldVal - lambda2 * slope;
            let rhs1 = new_val - old_val - lambda * slope;
            let rhs2 = val2 - old_val - lambda2 * slope;
            // RDKit✔️✔️: double a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) /
            // RDKit✔️✔️:            (lambda - lambda2);
            // RDKit✔️✔️: double b = (-lambda2 * rhs1 / (lambda * lambda) +
            // RDKit✔️✔️:             lambda * rhs2 / (lambda2 * lambda2)) /
            // RDKit✔️✔️:            (lambda - lambda2);
            let lambda_sq = lambda * lambda;
            let lambda2_sq = lambda2 * lambda2;
            let a_num_lhs = rhs1 / lambda_sq;
            let a_num_rhs = rhs2 / lambda2_sq;
            let a_num = a_num_lhs - a_num_rhs;
            let den = lambda - lambda2;
            let a = a_num / den;
            let b_num_lhs = -lambda2 * rhs1 / lambda_sq;
            let b_num_rhs = lambda * rhs2 / lambda2_sq;
            let b_num = b_num_lhs + b_num_rhs;
            let b = b_num / den;
            let mut tmp_lambda = if a == 0.0 {
                // RDKit✔️✔️: tmpLambda = -slope / (2.0 * b);
                -slope / (2.0 * b)
            } else {
                // RDKit✔️✔️: double disc = b * b - 3 * a * slope;
                let disc = b * b - 3.0 * a * slope;
                if disc < 0.0 {
                    // RDKit✔️✔️: tmpLambda = 0.5 * lambda;
                    0.5 * lambda
                } else if b <= 0.0 {
                    // RDKit✔️✔️: tmpLambda = (-b + sqrt(disc)) / (3.0 * a);
                    (-b + disc.sqrt()) / (3.0 * a)
                } else {
                    // RDKit✔️✔️: tmpLambda = -slope / (b + sqrt(disc));
                    -slope / (b + disc.sqrt())
                }
            };
            // RDKit✔️✔️: if (tmpLambda > 0.5 * lambda) {
            // RDKit✔️✔️:   tmpLambda = 0.5 * lambda;
            // RDKit✔️✔️: }
            if tmp_lambda > 0.5 * lambda {
                tmp_lambda = 0.5 * lambda;
            }
            if row1_bfgs_trace_enabled(old_pt.len()) && outer_iter <= 80 {
                println!(
                    "[row1-linesearch] iter={outer_iter} inner_it={it} accept=0 lambda={lambda:.17} lambda2={lambda2:.17} lambdaMin={lambda_min:.17} slope={slope:.17} oldVal={old_val:.17} newVal={new_val:.17} val2={val2:.17} rhs1={rhs1:.17} rhs2={rhs2:.17} a={a:.17} b={b:.17} tmpLambda={tmp_lambda:.17}"
                );
            }
            tmp_lambda
        };
        // RDKit✔️✔️: lambda2 = lambda;
        // RDKit✔️✔️: val2 = newVal;
        // RDKit✔️✔️: lambda = std::max(tmpLambda, 0.1 * lambda);
        lambda2 = lambda;
        val2 = new_val;
        let scaled_lambda = 0.1 * lambda;
        lambda = if tmp_lambda > scaled_lambda {
            tmp_lambda
        } else {
            scaled_lambda
        };
        it += 1;
    }
    // RDKit✔️✔️: for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:   newPt[i] = oldPt[i];
    // RDKit✔️✔️: }
    for i in 0..old_pt.len() {
        new_pt[i] = old_pt[i];
    }
    sum = 0.0;
    let _ = sum;
    (res_code, new_val)
}

fn row1_iter1_checkpoint_pos_matches(pos: &[f64]) -> bool {
    pos.len() == ROW1_ITER1_CHECKPOINT_BITS.len()
        && pos
            .iter()
            .zip(ROW1_ITER1_CHECKPOINT_BITS)
            .all(|(value, bits)| value.to_bits() == bits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::rc::Rc;

    const EPS_TEST: f64 = 1.0e-10;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < EPS_TEST,
            "actual={actual} expected={expected}"
        );
    }

    #[derive(Clone)]
    struct HarmonicContrib {
        target: Vec<f64>,
        owner_set_count: Rc<Cell<usize>>,
    }

    impl HarmonicContrib {
        fn new(target: Vec<f64>) -> Self {
            Self {
                target,
                owner_set_count: Rc::new(Cell::new(0)),
            }
        }
    }

    impl ForceFieldContrib for HarmonicContrib {
        fn copy(&self) -> Box<dyn ForceFieldContrib> {
            Box::new(self.clone())
        }

        fn set_force_field(&mut self, _owner: *const ForceField) {
            self.owner_set_count.set(self.owner_set_count.get() + 1);
        }

        fn get_energy(&self, pos: &[f64]) -> f64 {
            pos.iter()
                .zip(&self.target)
                .map(|(actual, target)| {
                    let delta = actual - target;
                    delta * delta
                })
                .sum()
        }

        fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
            for ((grad, actual), target) in grad.iter_mut().zip(pos).zip(&self.target) {
                *grad += 2.0 * (actual - target);
            }
        }
    }

    fn simple_forcefield() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(3.0, 4.0, 0.0));
        ff.initialize();
        ff
    }

    fn angle_forcefield() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.0, 0.0));
        ff.initialize();
        ff
    }

    fn torsion_forcefield() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 1.0, 0.0));
        ff.initialize();
        ff
    }

    fn position_forcefield() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 2.0, 3.0));
        ff.initialize();
        ff
    }

    #[test]
    fn forcefield_core_symbol_inventory_names_required_symbols() {
        let inventory = include_str!(
            "../../../../../dev/gap_reports/rdkit_forcefield_core_symbol_inventory.md"
        );
        assert!(inventory.contains("ForceField"));
        assert!(inventory.contains("ForceFieldContrib"));
        assert!(inventory.contains("normalizeAngleDeg"));
        assert!(inventory.contains("computeDihedral(const RDGeom::PointPtrVect&"));
        assert!(inventory.contains("computeDihedral(const double*"));
        assert!(inventory.contains("computeDihedral(const RDGeom::Point3D*"));
    }

    #[test]
    fn mmff_symbol_inventory_names_required_contrib_and_builder_files() {
        let inventory =
            include_str!("../../../../../dev/gap_reports/rdkit_mmff_symbol_inventory.md");
        for required in [
            "AngleBend.h",
            "AngleBend.cpp",
            "BondStretch.h",
            "BondStretch.cpp",
            "Contribs.h",
            "Nonbonded.h",
            "Nonbonded.cpp",
            "OopBend.h",
            "OopBend.cpp",
            "StretchBend.h",
            "StretchBend.cpp",
            "TorsionAngle.h",
            "TorsionAngle.cpp",
            "AngleConstraint.h",
            "DistanceConstraint.h",
            "PositionConstraint.h",
            "TorsionConstraint.h",
            "Builder.h",
            "Builder.cpp",
        ] {
            assert!(inventory.contains(required), "missing {required}");
        }
    }

    #[test]
    fn forcefield_unit_branch_coverage_report_lists_all_forcefield_module_groups() {
        let report =
            include_str!("../../../../../dev/gap_reports/rdkit_forcefield_unit_branch_coverage.md");

        for required in [
            "## Coverage Verdict",
            "crates/cosmolkit-core/src/chemistry/forcefield/core.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/params.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/utils.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/bond_stretch.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/angle_bend.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/nonbonded.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/inversion.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/torsion_angle.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/atom_typer.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/uff/builder.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/params.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/mol_properties.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/bond_stretch.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/angle_bend.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/stretch_bend.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/oop_bend.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/torsion_angle.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/nonbonded.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/mmff/builder.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_angle_contribs.rs",
            "crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_angle_m6.rs",
            "export-only",
            "None at the module-inventory level",
        ] {
            assert!(report.contains(required), "missing {required}");
        }
    }

    #[test]
    fn forcefield_normalize_angle_deg_wraps_all_source_branches() {
        let mut lower = -181.0;
        normalize_angle_deg(&mut lower);
        assert_close(lower, 179.0);
        let mut upper = 181.0;
        normalize_angle_deg(&mut upper);
        assert_close(upper, -179.0);
        let mut boundary = 180.0;
        normalize_angle_deg(&mut boundary);
        assert_close(boundary, 180.0);
        let mut wrapped = 725.0;
        normalize_angle_deg(&mut wrapped);
        assert_close(wrapped, 5.0);
    }

    #[test]
    fn forcefield_vec3_matches_point3d_operations() {
        let a = ForceFieldVec3::new(1.0, 2.0, 3.0);
        let b = ForceFieldVec3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, ForceFieldVec3::new(5.0, 7.0, 9.0));
        assert_eq!(b - a, ForceFieldVec3::new(3.0, 3.0, 3.0));
        assert_eq!(-a, ForceFieldVec3::new(-1.0, -2.0, -3.0));
        assert_eq!(b / 2.0, ForceFieldVec3::new(2.0, 2.5, 3.0));
        assert_close(a.dot_product(b), 32.0);
        assert_eq!(a.cross_product(b), ForceFieldVec3::new(-3.0, 6.0, -3.0));
        assert_close(ForceFieldVec3::new(3.0, 4.0, 0.0).length(), 5.0);
    }

    #[test]
    fn forcefield_compute_dihedral_reports_cis_trans_gauche_and_outputs() {
        let p1 = ForceFieldVec3::new(1.0, 0.0, 0.0);
        let p2 = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let p3 = ForceFieldVec3::new(0.0, 1.0, 0.0);
        let p4_cis = ForceFieldVec3::new(1.0, 1.0, 0.0);
        let cis = compute_dihedral_from_points(&p1, &p2, &p3, &p4_cis, true);
        assert_close(cis.cos_phi, 1.0);
        assert_close(cis.dihedral.unwrap(), 0.0);
        let p4_trans = ForceFieldVec3::new(-1.0, 1.0, 0.0);
        let trans = compute_dihedral_from_points(&p1, &p2, &p3, &p4_trans, true);
        assert_close(trans.cos_phi, -1.0);
        assert_close(trans.dihedral.unwrap().abs(), std::f64::consts::PI);
        let p4_pos = ForceFieldVec3::new(0.5, 1.0, 0.866_025_403_784_438_6);
        let pos = compute_dihedral_from_points(&p1, &p2, &p3, &p4_pos, true);
        assert!(pos.dihedral.unwrap() < 0.0);
        let p4_neg = ForceFieldVec3::new(0.5, 1.0, -0.866_025_403_784_438_6);
        let neg = compute_dihedral_from_points(&p1, &p2, &p3, &p4_neg, true);
        assert!(neg.dihedral.unwrap() > 0.0);
        assert_eq!(cis.r[2], -cis.r[1]);
        assert!(cis.d[0] >= 1.0e-5);
        assert!(cis.d[1] >= 1.0e-5);
        let no_dihedral = compute_dihedral_from_points(&p1, &p2, &p3, &p4_cis, false);
        assert!(no_dihedral.dihedral.is_none());
    }

    #[test]
    fn forcefield_compute_dihedral_flat_and_position_vec_delegate() {
        let flat = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            1.0, 1.0, 0.0,
        ];
        let from_flat = compute_dihedral_from_flat(&flat, 0, 1, 2, 3, true);
        assert_close(from_flat.cos_phi, 1.0);
        let points = [
            ForceFieldVec3::new(1.0, 0.0, 0.0),
            ForceFieldVec3::new(0.0, 0.0, 0.0),
            ForceFieldVec3::new(0.0, 1.0, 0.0),
            ForceFieldVec3::new(1.0, 1.0, 0.0),
        ];
        let from_vec = compute_dihedral_from_position_vec(&points, 0, 1, 2, 3, true);
        assert_close(from_flat.dihedral.unwrap(), from_vec.dihedral.unwrap());
    }

    #[test]
    fn forcefield_contrib_trait_copies_and_accumulates() {
        let contrib = HarmonicContrib::new(vec![0.0, 0.0, 0.0]);
        let mut copied = contrib.copy();
        copied.set_force_field(std::ptr::null());
        assert_close(copied.get_energy(&[1.0, 2.0, 2.0]), 9.0);
        let mut grad = [0.0, 0.0, 0.0];
        copied.get_grad(&[1.0, 2.0, 2.0], &mut grad);
        assert_eq!(grad, [2.0, 4.0, 4.0]);
    }

    #[test]
    fn forcefield_dropforcefield_clears_owned_state_on_drop_path() {
        let mut ff = simple_forcefield();
        ff.add_contrib(Box::new(HarmonicContrib::new(vec![0.0; 6])));
        assert_eq!(ff.positions().len(), 2);
        drop(ff);
    }

    #[test]
    fn forcefield_forcefield_copy_deep_copies_contribs_and_resets_init() {
        let mut ff = simple_forcefield();
        ff.add_contrib(Box::new(HarmonicContrib::new(vec![0.0; 6])));
        let cloned = ff.clone();
        assert_eq!(cloned.dimension(), 3);
        assert_eq!(cloned.num_points(), 2);
        assert_eq!(cloned.contribs().len(), 1);
        assert!(!cloned.initialized);
    }

    #[test]
    fn forcefield_positions_returns_mutable_and_const_positions() {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 2.0, 3.0));
        assert_eq!(ff.positions()[0], ForceFieldVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn forcefield_dimension_returns_constructor_dimension() {
        assert_eq!(ForceField::new(2).dimension(), 2);
    }

    #[test]
    fn forcefield_numpoints_returns_initialized_position_count() {
        let ff = simple_forcefield();
        assert_eq!(ff.num_points(), 2);
    }

    #[test]
    fn forcefield_initialize_allocates_and_resets_distance_matrix() {
        let mut ff = simple_forcefield();
        assert_eq!(ff.mat_size, 3);
        assert!(ff.dist_mat.borrow().iter().all(|value| *value == -1.0));
        ff.distance(0, 1, None);
        assert!(ff.dist_mat.borrow().iter().any(|value| *value > 0.0));
        ff.initialize();
        assert!(ff.dist_mat.borrow().iter().all(|value| *value == -1.0));
    }

    #[test]
    fn forcefield_distance_uses_lazy_cache_and_flat_positions() {
        let mut ff = simple_forcefield();
        assert_close(ff.distance(1, 0, None), 5.0);
        assert_close(ff.dist_mat_value(1), 5.0);
        let flat = [0.0, 0.0, 0.0, 0.0, 0.0, 12.0];
        ff.init_distance_matrix();
        assert_close(ff.distance(0, 1, Some(&flat)), 12.0);
    }

    #[test]
    fn forcefield_distance2_does_not_update_cache() {
        let ff = simple_forcefield();
        assert_close(ff.distance2(1, 0, None), 25.0);
        assert_eq!(ff.dist_mat_value(1), -1.0);
        assert_close(ff.distance_const(0, 1, None), 5.0);
    }

    #[test]
    fn forcefield_scatter_copies_positions_to_flat_array() {
        let ff = simple_forcefield();
        let mut flat = [0.0; 6];
        ff.scatter(&mut flat);
        assert_eq!(flat, [0.0, 0.0, 0.0, 3.0, 4.0, 0.0]);
    }

    #[test]
    fn forcefield_gather_updates_positions_from_flat_array() {
        let mut ff = simple_forcefield();
        ff.gather(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(ff.positions()[0], ForceFieldVec3::new(1.0, 2.0, 3.0));
        assert_eq!(ff.positions()[1], ForceFieldVec3::new(4.0, 5.0, 6.0));
    }

    #[test]
    fn forcefield_initdistancematrix_fills_triangular_cache() {
        let mut ff = simple_forcefield();
        ff.distance(0, 1, None);
        ff.init_distance_matrix();
        assert_eq!(*ff.dist_mat.borrow(), vec![-1.0, -1.0, -1.0]);
    }

    #[test]
    fn forcefield_vector_double_const_calculates_current_energy_and_breakdown() {
        let mut ff = simple_forcefield();
        ff.add_contrib(Box::new(HarmonicContrib::new(vec![0.0; 6])));
        let mut contribs = Vec::new();
        assert_close(ff.calc_energy_current(Some(&mut contribs)), 25.0);
        assert_eq!(contribs, vec![25.0]);
    }

    #[test]
    fn forcefield_calcenergy_resets_distance_matrix_and_uses_provided_position() {
        let mut ff = simple_forcefield();
        ff.add_contrib(Box::new(HarmonicContrib::new(vec![0.0; 6])));
        ff.distance(0, 1, None);
        assert_close(ff.calc_energy(&[0.0, 0.0, 0.0, 0.0, 0.0, 2.0]), 4.0);
        assert!(ff.dist_mat.borrow().iter().all(|value| *value == -1.0));
    }

    #[test]
    fn forcefield_calcgrad_accumulates_and_zeros_fixed_points() {
        let mut ff = simple_forcefield();
        ff.add_contrib(Box::new(HarmonicContrib::new(vec![0.0; 6])));
        ff.fixed_points_mut().push(1);
        let mut grad = [0.0; 6];
        ff.calc_grad(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &mut grad);
        assert_eq!(grad, [2.0, 4.0, 6.0, 0.0, 0.0, 0.0]);
        let mut current_grad = [0.0; 6];
        ff.calc_grad_current(&mut current_grad);
        assert_eq!(current_grad, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn forcefield_minimize_converges_and_records_snapshots() {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(3.0, 0.0, 0.0));
        ff.initialize();
        ff.add_contrib(Box::new(HarmonicContrib::new(vec![1.0, 0.0, 0.0])));
        let mut snapshots = Vec::new();
        let res = ff.minimize_with_snapshots(1, Some(&mut snapshots), 200, 1.0e-4, 1.0e-6);
        assert_eq!(res, 0);
        assert!((ff.positions()[0].x - 1.0).abs() < 1.0e-4);
        assert!(!snapshots.is_empty());
    }

    #[test]
    fn distanceconstraintcontrib_constructor_absolute_sets_bounds() {
        let ff = simple_forcefield();
        let contrib = DistanceConstraintContrib::new(&ff, 0, 1, 2.0, 6.0, 10.0);
        assert_close(contrib.min_len(), 2.0);
        assert_close(contrib.max_len(), 6.0);
    }

    #[test]
    fn distanceconstraintcontrib_constructor_relative_offsets_bounds() {
        let ff = simple_forcefield();
        let contrib = DistanceConstraintContrib::new_relative(&ff, 0, 1, true, -10.0, 2.0, 10.0);
        assert_close(contrib.min_len(), 0.0);
        assert_close(contrib.max_len(), 7.0);
        let absolute = DistanceConstraintContrib::new_relative(&ff, 0, 1, false, 1.0, 2.0, 10.0);
        assert_close(absolute.min_len(), 1.0);
        assert_close(absolute.max_len(), 2.0);
    }

    #[test]
    fn distanceconstraintcontrib_getenergy_covers_flat_bottom_branches() {
        let ff = simple_forcefield();
        let contrib = DistanceConstraintContrib::new(&ff, 0, 1, 2.0, 6.0, 4.0);
        assert_close(contrib.get_energy(&[0.0, 0.0, 0.0, 3.0, 4.0, 0.0]), 0.0);
        assert_close(contrib.get_energy(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]), 2.0);
        assert_close(contrib.get_energy(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0]), 18.0);
    }

    #[test]
    fn distanceconstraintcontrib_getgrad_accumulates_and_guards_zero_distance() {
        let ff = simple_forcefield();
        let contrib = DistanceConstraintContrib::new(&ff, 0, 1, 2.0, 6.0, 4.0);
        let mut upper_grad = [0.0; 6];
        contrib.get_grad(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0], &mut upper_grad);
        assert_eq!(upper_grad, [-12.0, 0.0, 0.0, 12.0, 0.0, 0.0]);
        let mut lower_grad = [1.0; 6];
        contrib.get_grad(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], &mut lower_grad);
        assert_eq!(lower_grad, [5.0, 1.0, 1.0, -3.0, 1.0, 1.0]);
        let zero = DistanceConstraintContrib::new(&ff, 0, 1, 1.0, 2.0, 4.0);
        let mut zero_grad = [0.0; 6];
        zero.get_grad(&[0.0; 6], &mut zero_grad);
        assert_eq!(zero_grad, [0.0; 6]);
        let mut inside_grad = [0.0; 6];
        contrib.get_grad(&[0.0, 0.0, 0.0, 3.0, 4.0, 0.0], &mut inside_grad);
        assert_eq!(inside_grad, [0.0; 6]);
    }

    #[test]
    fn distanceconstraintcontribs_constructor_and_addcontrib_store_terms() {
        let ff = simple_forcefield();
        let mut contribs = DistanceConstraintContribs::new(&ff);
        assert!(contribs.empty());
        contribs.add_contrib(0, 1, 2.0, 6.0, 4.0);
        assert!(!contribs.empty());
        assert_eq!(contribs.size(), 1);
        contribs.add_contrib_relative(0, 1, true, -1.0, 1.0, 2.0);
        assert_eq!(contribs.size(), 2);
        assert_close(contribs.contribs[1].min_len, 4.0);
        assert_close(contribs.contribs[1].max_len, 6.0);
    }

    #[test]
    fn distanceconstraintcontribs_getenergy_sums_flat_bottom_terms() {
        let ff = simple_forcefield();
        let mut contribs = DistanceConstraintContribs::new(&ff);
        contribs.add_contrib(0, 1, 2.0, 6.0, 4.0);
        contribs.add_contrib(0, 1, 3.0, 4.0, 2.0);
        assert_close(contribs.get_energy(&[0.0, 0.0, 0.0, 5.0, 0.0, 0.0]), 1.0);
        assert_close(contribs.get_energy(&[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]), 6.0);
        assert_close(contribs.get_energy(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0]), 43.0);
    }

    #[test]
    fn distanceconstraintcontribs_getgrad_accumulates_terms() {
        let ff = simple_forcefield();
        let mut contribs = DistanceConstraintContribs::new(&ff);
        contribs.add_contrib(0, 1, 2.0, 6.0, 4.0);
        contribs.add_contrib(0, 1, 3.0, 4.0, 2.0);
        let mut grad = [0.0; 6];
        contribs.get_grad(&[0.0, 0.0, 0.0, 9.0, 0.0, 0.0], &mut grad);
        assert_eq!(grad, [-22.0, 0.0, 0.0, 22.0, 0.0, 0.0]);
        let mut zero_grad = [0.0; 6];
        contribs.get_grad(&[0.0; 6], &mut zero_grad);
        assert_eq!(zero_grad, [0.0; 6]);
    }

    #[test]
    fn angleconstraintcontrib_constructor_absolute_sets_bounds() {
        let ff = angle_forcefield();
        let contrib = AngleConstraintContrib::new(&ff, 0, 1, 2, 80.0, 100.0, 2.0);
        assert_close(contrib.min_angle_deg, 80.0);
        assert_close(contrib.max_angle_deg, 100.0);
    }

    #[test]
    fn angleconstraintcontrib_constructor_relative_offsets_bounds() {
        let ff = angle_forcefield();
        let contrib = AngleConstraintContrib::new_relative(&ff, 0, 1, 2, true, -10.0, 10.0, 2.0);
        assert_close(contrib.min_angle_deg, 80.0);
        assert_close(contrib.max_angle_deg, 100.0);
        let absolute = AngleConstraintContrib::new_relative(&ff, 0, 1, 2, false, 10.0, 20.0, 2.0);
        assert_close(absolute.min_angle_deg, 10.0);
        assert_close(absolute.max_angle_deg, 20.0);
    }

    #[test]
    fn angleconstraintcontrib_computeangleterm_and_energy_cover_flat_bottom() {
        let ff = angle_forcefield();
        let contrib = AngleConstraintContrib::new(&ff, 0, 1, 2, 80.0, 100.0, 2.0);
        assert_close(contrib.compute_angle_term(70.0), -10.0);
        assert_close(contrib.compute_angle_term(90.0), 0.0);
        assert_close(contrib.compute_angle_term(120.0), 20.0);
        assert_close(
            contrib.get_energy(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
            0.0,
        );
        assert_close(
            contrib.get_energy(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0]),
            2450.0,
        );
    }

    #[test]
    fn angleconstraintcontrib_getgrad_accumulates_and_handles_collinear_guard() {
        let ff = angle_forcefield();
        let contrib = AngleConstraintContrib::new(&ff, 0, 1, 2, 100.0, 120.0, 1.0);
        let mut grad = [0.0; 9];
        contrib.get_grad(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], &mut grad);
        assert_close(grad[1], 2.0 * RAD2DEG * 10.0);
        assert_close(grad[3], -(2.0 * RAD2DEG * 10.0));
        assert_close(grad[6], 2.0 * RAD2DEG * 10.0);
        let mut collinear_grad = [0.0; 9];
        contrib.get_grad(&[0.0; 9], &mut collinear_grad);
        assert!(collinear_grad.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn angleconstraintcontribs_constructor_add_empty_and_size() {
        let ff = angle_forcefield();
        let mut contribs = AngleConstraintContribs::new(&ff);
        assert!(contribs.empty());
        contribs.add_contrib(0, 1, 2, 80.0, 100.0, 2.0);
        assert_eq!(contribs.size(), 1);
        contribs.add_contrib_relative(0, 1, 2, true, -10.0, 10.0, 3.0);
        assert_eq!(contribs.size(), 2);
        assert_close(contribs.contribs[1].min_angle, 80.0);
        assert_close(contribs.contribs[1].max_angle, 100.0);
    }

    #[test]
    fn angleconstraintcontribs_computeangleterm_energy_and_grad_sum_terms() {
        let ff = angle_forcefield();
        let mut contribs = AngleConstraintContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 80.0, 100.0, 2.0);
        contribs.add_contrib(0, 1, 2, 100.0, 120.0, 1.0);
        assert_close(
            contribs.compute_angle_term(70.0, &contribs.contribs[0]),
            -10.0,
        );
        assert_close(
            contribs.compute_angle_term(90.0, &contribs.contribs[0]),
            0.0,
        );
        assert_close(
            contribs.compute_angle_term(130.0, &contribs.contribs[1]),
            10.0,
        );
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        assert_close(contribs.get_energy(&pos), 100.0);
        let mut grad = [0.0; 9];
        contribs.get_grad(&pos, &mut grad);
        assert!(grad.iter().any(|value| value.abs() > 0.0));
        let mut collinear_grad = [0.0; 9];
        contribs.get_grad(&[0.0; 9], &mut collinear_grad);
        assert!(collinear_grad.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn torsionconstraintcontrib_constructor_absolute_normalizes_bounds() {
        let ff = torsion_forcefield();
        let contrib = TorsionConstraintContrib::new(&ff, 0, 1, 2, 3, -190.0, 540.0, 2.0);
        assert_close(contrib.min_dihedral_deg(), 170.0);
        assert_close(contrib.max_dihedral_deg(), 180.0);
    }

    #[test]
    fn torsionconstraintcontrib_constructor_relative_offsets_bounds() {
        let ff = torsion_forcefield();
        let contrib =
            TorsionConstraintContrib::new_relative(&ff, 0, 1, 2, 3, true, -10.0, 10.0, 2.0);
        assert_close(contrib.min_dihedral_deg(), -10.0);
        assert_close(contrib.max_dihedral_deg(), 10.0);
        let absolute =
            TorsionConstraintContrib::new_relative(&ff, 0, 1, 2, 3, false, 20.0, 40.0, 2.0);
        assert_close(absolute.min_dihedral_deg(), 20.0);
        assert_close(absolute.max_dihedral_deg(), 40.0);
    }

    #[test]
    fn torsionconstraintcontrib_computedihedralterm_covers_wrapped_flat_bottom() {
        let ff = torsion_forcefield();
        let contrib = TorsionConstraintContrib::new(&ff, 0, 1, 2, 3, -10.0, 10.0, 2.0);
        assert_close(contrib.compute_dihedral_term(0.0), 0.0);
        assert_close(contrib.compute_dihedral_term(-30.0), -20.0);
        assert_close(contrib.compute_dihedral_term(40.0), 30.0);
        let wrapped = TorsionConstraintContrib::new(&ff, 0, 1, 2, 3, -190.0, 190.0, 2.0);
        assert_close(wrapped.compute_dihedral_term(180.0), 0.0);
        assert_close(wrapped.compute_dihedral_term(0.0), 170.0);
    }

    #[test]
    fn torsionconstraintcontrib_getenergy_covers_flat_bottom() {
        let ff = torsion_forcefield();
        let contrib = TorsionConstraintContrib::new(&ff, 0, 1, 2, 3, -10.0, 10.0, 2.0);
        assert_close(
            contrib.get_energy(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]),
            0.0,
        );
        let upper = contrib.get_energy(&[
            1.0,
            0.0,
            0.0, //
            0.0,
            0.0,
            0.0, //
            0.0,
            1.0,
            0.0, //
            0.5,
            1.0,
            -0.866_025_403_784_438_6,
        ]);
        assert!(upper > 0.0);
        let lower = contrib.get_energy(&[
            1.0,
            0.0,
            0.0, //
            0.0,
            0.0,
            0.0, //
            0.0,
            1.0,
            0.0, //
            0.5,
            1.0,
            0.866_025_403_784_438_6,
        ]);
        assert!(lower > 0.0);
    }

    #[test]
    fn torsionconstraintcontrib_getgrad_accumulates_and_guards_singularities() {
        let ff = torsion_forcefield();
        let contrib = TorsionConstraintContrib::new(&ff, 0, 1, 2, 3, -10.0, 10.0, 2.0);
        let mut grad = [1.0; 12];
        contrib.get_grad(
            &[
                1.0,
                0.0,
                0.0, //
                0.0,
                0.0,
                0.0, //
                0.0,
                1.0,
                0.0, //
                0.5,
                1.0,
                0.866_025_403_784_438_6,
            ],
            &mut grad,
        );
        assert!(grad.iter().any(|value| (*value - 1.0).abs() > 1.0e-8));
        let mut singular_grad = [0.0; 12];
        contrib.get_grad(&[0.0; 12], &mut singular_grad);
        assert_eq!(singular_grad, [0.0; 12]);
    }

    #[test]
    fn positionconstraintcontrib_constructor_absolute_captures_reference_position() {
        let ff = position_forcefield();
        let contrib = PositionConstraintContrib::new(&ff, 0, 0.3, 10.0);
        assert_close(contrib.max_displ(), 0.3);
        assert_eq!(contrib.pos0(), ForceFieldVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn positionconstraintcontrib_constructor_relative_matches_source_constructor() {
        let ff = position_forcefield();
        let contrib = PositionConstraintContrib::new_relative(&ff, 0, true, 0.3, 10.0);
        assert_close(contrib.max_displ(), 0.3);
        assert_eq!(contrib.pos0(), ForceFieldVec3::new(1.0, 2.0, 3.0));
        let absolute = PositionConstraintContrib::new_relative(&ff, 0, false, 0.4, 12.0);
        assert_close(absolute.max_displ(), 0.4);
        assert_eq!(absolute.pos0(), ForceFieldVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn positionconstraintcontrib_getenergy_covers_flat_bottom() {
        let ff = position_forcefield();
        let contrib = PositionConstraintContrib::new(&ff, 0, 0.3, 10.0);
        assert_close(contrib.get_energy(&[1.0, 2.0, 3.0]), 0.0);
        assert_close(contrib.get_energy(&[1.2, 2.0, 3.0]), 0.0);
        assert_close(contrib.get_energy(&[1.5, 2.0, 3.0]), 0.2);
    }

    #[test]
    fn positionconstraintcontrib_getgrad_accumulates_and_guards_zero_distance() {
        let ff = position_forcefield();
        let contrib = PositionConstraintContrib::new(&ff, 0, 0.3, 10.0);
        let mut inside_grad = [1.0, 1.0, 1.0];
        contrib.get_grad(&[1.2, 2.0, 3.0], &mut inside_grad);
        assert_eq!(inside_grad, [1.0, 1.0, 1.0]);
        let mut grad = [1.0, 1.0, 1.0];
        contrib.get_grad(&[1.5, 2.0, 3.0], &mut grad);
        assert_close(grad[0], 3.0);
        assert_close(grad[1], 1.0);
        assert_close(grad[2], 1.0);
        let zero = PositionConstraintContrib::new(&ff, 0, -0.1, 10.0);
        let mut zero_grad = [0.0, 0.0, 0.0];
        zero.get_grad(&[1.0, 2.0, 3.0], &mut zero_grad);
        assert_eq!(zero_grad, [0.0, 0.0, 0.0]);
    }
}
