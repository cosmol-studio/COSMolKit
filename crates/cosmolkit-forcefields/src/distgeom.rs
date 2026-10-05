//! Complete RDKit distance/chiral/fourth-dimension force-field contributions.
use crate::geometry::Point3;
use crate::kernel::{EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError};
use std::sync::Arc;
fn point3(x: f64, y: f64, z: f64) -> Point3 {
    Point3 { x, y, z }
}
// BEGIN RDKIT CPP ENUM DistGeom::ChiralSetStructureFlags (ChiralSet.h:18-21)
// RDKit✔️✔️: enum class ChiralSetStructureFlags : std::uint64_t {
// RDKit✔️✔️:   IN_FUSED_SMALL_RINGS =
// RDKit✔️✔️:       1 << 0,  // a chiral center involved in fusing 2 or more small rings
// RDKit✔️✔️: };
// END RDKIT CPP ENUM DistGeom::ChiralSetStructureFlags
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u64)]
pub enum ChiralSetStructureFlags {
    InFusedSmallRings = 1 << 0,
}

/// Class used to store a quartet of points and chiral volume bounds on them.
#[derive(Debug, Clone, PartialEq)]
pub struct ChiralSet {
    pub idx0: usize,
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub volume_lower_bound: f64,
    pub volume_upper_bound: f64,
    pub structure_flags: u64,
}

impl ChiralSet {
    #[must_use]
    pub fn new(
        pid0: usize,
        pid1: usize,
        pid2: usize,
        pid3: usize,
        pid4: usize,
        lower_vol_bound: f64,
        upper_vol_bound: f64,
        structure_flags: u64,
    ) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::ChiralSet::ChiralSet (ChiralSet.h:40-57)
        // RDKit✔️✔️: ChiralSet(unsigned int pid0, unsigned int pid1, unsigned int pid2,
        // RDKit✔️✔️:           unsigned int pid3, unsigned int pid4, double lowerVolBound,
        // RDKit✔️✔️:           double upperVolBound, std::uint64_t structureFlags = 0)
        // RDKit✔️✔️:     : d_idx0(pid0),
        // RDKit✔️✔️:       d_idx1(pid1),
        // RDKit✔️✔️:       d_idx2(pid2),
        // RDKit✔️✔️:       d_idx3(pid3),
        // RDKit✔️✔️:       d_idx4(pid4),
        // RDKit✔️✔️:       d_volumeLowerBound(lowerVolBound),
        // RDKit✔️✔️:       d_volumeUpperBound(upperVolBound),
        // RDKit✔️✔️:       d_structureFlags(structureFlags) {
        // RDKit✔️✔️:   CHECK_INVARIANT(lowerVolBound <= upperVolBound, "Inconsistent bounds\n");
        // RDKit✔️✔️:   d_volumeLowerBound = lowerVolBound;
        // RDKit✔️✔️:   d_volumeUpperBound = upperVolBound;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::ChiralSet::ChiralSet
        assert!(lower_vol_bound <= upper_vol_bound, "Inconsistent bounds\n");
        Self {
            idx0: pid0,
            idx1: pid1,
            idx2: pid2,
            idx3: pid3,
            idx4: pid4,
            volume_lower_bound: lower_vol_bound,
            volume_upper_bound: upper_vol_bound,
            structure_flags,
        }
    }

    #[must_use]
    pub fn with_default_structure_flags(
        pid0: usize,
        pid1: usize,
        pid2: usize,
        pid3: usize,
        pid4: usize,
        lower_vol_bound: f64,
        upper_vol_bound: f64,
    ) -> Self {
        Self::new(
            pid0,
            pid1,
            pid2,
            pid3,
            pid4,
            lower_vol_bound,
            upper_vol_bound,
            0,
        )
    }

    #[must_use]
    pub fn get_upper_volume_bound(&self) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralSet::getUpperVolumeBound (ChiralSet.h:59)
        // RDKit✔️✔️: inline double getUpperVolumeBound() const { return d_volumeUpperBound; }
        // END RDKIT CPP METHOD DistGeom::ChiralSet::getUpperVolumeBound
        self.volume_upper_bound
    }

    #[must_use]
    pub fn get_lower_volume_bound(&self) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralSet::getLowerVolumeBound (ChiralSet.h:61)
        // RDKit✔️✔️: inline double getLowerVolumeBound() const { return d_volumeLowerBound; }
        // END RDKIT CPP METHOD DistGeom::ChiralSet::getLowerVolumeBound
        self.volume_lower_bound
    }
}

// BEGIN RDKIT CPP TYPEDEFS DistGeom::ChiralSetPtr/VECT_CHIRALSET (ChiralSet.h:64-65)
// RDKit✔️✔️: typedef boost::shared_ptr<ChiralSet> ChiralSetPtr;
// RDKit✔️✔️: typedef std::vector<ChiralSetPtr> VECT_CHIRALSET;
// END RDKIT CPP TYPEDEFS DistGeom::ChiralSetPtr/VECT_CHIRALSET
pub type ChiralSetPtr = Arc<ChiralSet>;
pub type VectChiralSet = Vec<ChiralSetPtr>;

#[must_use]
pub fn calc_chiral_volume_flat(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    pos: &[f64],
    dim: usize,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::calcChiralVolume flat overload (ChiralViolationContribs.cpp:15-35)
    // RDKit✔️✔️: double calcChiralVolume(const unsigned int idx1, const unsigned int idx2,
    // RDKit✔️✔️:                         const unsigned int idx3, const unsigned int idx4,
    // RDKit✔️✔️:                         const double *pos, const unsigned int dim) {
    // RDKit✔️✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit✔️✔️:   // calculated using only the first 3 dimensions
    // RDKit✔️✔️:   RDGeom::Point3D v1(pos[idx1 * dim] - pos[idx4 * dim],
    // RDKit✔️✔️:                      pos[idx1 * dim + 1] - pos[idx4 * dim + 1],
    // RDKit✔️✔️:                      pos[idx1 * dim + 2] - pos[idx4 * dim + 2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2(pos[idx2 * dim] - pos[idx4 * dim],
    // RDKit✔️✔️:                      pos[idx2 * dim + 1] - pos[idx4 * dim + 1],
    // RDKit✔️✔️:                      pos[idx2 * dim + 2] - pos[idx4 * dim + 2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v3(pos[idx3 * dim] - pos[idx4 * dim],
    // RDKit✔️✔️:                      pos[idx3 * dim + 1] - pos[idx4 * dim + 1],
    // RDKit✔️✔️:                      pos[idx3 * dim + 2] - pos[idx4 * dim + 2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double vol = v1.dotProduct(v2xv3);
    // RDKit✔️✔️:   return vol;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::calcChiralVolume flat overload
    let v1 = point3(
        pos[idx1 * dim] - pos[idx4 * dim],
        pos[idx1 * dim + 1] - pos[idx4 * dim + 1],
        pos[idx1 * dim + 2] - pos[idx4 * dim + 2],
    );
    let v2 = point3(
        pos[idx2 * dim] - pos[idx4 * dim],
        pos[idx2 * dim + 1] - pos[idx4 * dim + 1],
        pos[idx2 * dim + 2] - pos[idx4 * dim + 2],
    );
    let v3 = point3(
        pos[idx3 * dim] - pos[idx4 * dim],
        pos[idx3 * dim + 1] - pos[idx4 * dim + 1],
        pos[idx3 * dim + 2] - pos[idx4 * dim + 2],
    );

    v1.dot_product(&v2.cross_product(&v3))
}

#[must_use]
pub fn calc_chiral_volume_points(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    pts: &[Point3],
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::calcChiralVolume PointPtrVect overload (ChiralViolationContribs.cpp:36-56)
    // RDKit✔️✔️: double calcChiralVolume(const unsigned int idx1, const unsigned int idx2,
    // RDKit✔️✔️:                         const unsigned int idx3, const unsigned int idx4,
    // RDKit✔️✔️:                         const RDGeom::PointPtrVect &pts) {
    // RDKit✔️✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit✔️✔️:   // calculated using only the first 3 dimensions
    // RDKit✔️✔️:   RDGeom::Point3D v1((*pts[idx1])[0] - (*pts[idx4])[0],
    // RDKit✔️✔️:                      (*pts[idx1])[1] - (*pts[idx4])[1],
    // RDKit✔️✔️:                      (*pts[idx1])[2] - (*pts[idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2((*pts[idx2])[0] - (*pts[idx4])[0],
    // RDKit✔️✔️:                      (*pts[idx2])[1] - (*pts[idx4])[1],
    // RDKit✔️✔️:                      (*pts[idx2])[2] - (*pts[idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v3((*pts[idx3])[0] - (*pts[idx4])[0],
    // RDKit✔️✔️:                      (*pts[idx3])[1] - (*pts[idx4])[1],
    // RDKit✔️✔️:                      (*pts[idx3])[2] - (*pts[idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double vol = v1.dotProduct(v2xv3);
    // RDKit✔️✔️:   return vol;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::calcChiralVolume PointPtrVect overload
    let v1 = Point3::difference(&pts[idx1], &pts[idx4]);
    let v2 = Point3::difference(&pts[idx2], &pts[idx4]);
    let v3 = Point3::difference(&pts[idx3], &pts[idx4]);

    v1.dot_product(&v2.cross_product(&v3))
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChiralViolationContribsParams {
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub vol_upper: f64,
    pub vol_lower: f64,
    pub weight: f64,
}

impl ChiralViolationContribsParams {
    #[must_use]
    pub fn new(
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        vol_upper: f64,
        vol_lower: f64,
        weight: f64,
    ) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribsParams (ChiralViolationContribs.h:25-37)
        // RDKit✔️✔️: ChiralViolationContribsParams(unsigned int i1, unsigned int i2,
        // RDKit✔️✔️:                               unsigned int i3, unsigned int i4, double u,
        // RDKit✔️✔️:                               double l, double w = 1.0)
        // RDKit✔️✔️:     : idx1(i1),
        // RDKit✔️✔️:       idx2(i2),
        // RDKit✔️✔️:       idx3(i3),
        // RDKit✔️✔️:       idx4(i4),
        // RDKit✔️✔️:       volUpper(u),
        // RDKit✔️✔️:       volLower(l),
        // RDKit✔️✔️:       weight(w) {};
        // END RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribsParams
        Self {
            idx1,
            idx2,
            idx3,
            idx4,
            vol_upper,
            vol_lower,
            weight,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChiralViolationContribs {
    owner_dimension: Option<usize>,
    owner_points: Option<usize>,
    contribs: Vec<ChiralViolationContribsParams>,
}

impl Default for ChiralViolationContribs {
    fn default() -> Self {
        Self {
            owner_dimension: None,
            owner_points: None,
            contribs: Vec::new(),
        }
    }
}

impl ChiralViolationContribs {
    #[must_use]
    pub fn new(owner: &ForceField<'_>) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribs::ChiralViolationContribs (ChiralViolationContribs.cpp:57-62)
        // RDKit✔️✔️: ChiralViolationContribs::ChiralViolationContribs(
        // RDKit✔️✔️:     ForceFields::ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribs::ChiralViolationContribs
        Self {
            owner_dimension: Some(owner.dimension() as usize),
            owner_points: Some(owner.positions().len()),
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, cset: &ChiralSet, weight: f64) {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralViolationContribs::addContrib (ChiralViolationContribs.cpp:63-76)
        // RDKit✔️✔️: void ChiralViolationContribs::addContrib(const ChiralSet *cset, double weight) {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(cset, "bad chiral set");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx4, dp_forceField->positions().size());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   d_contribs.emplace_back(cset->d_idx1, cset->d_idx2, cset->d_idx3,
        // RDKit✔️✔️:                           cset->d_idx4, cset->getUpperVolumeBound(),
        // RDKit✔️✔️:                           cset->getLowerVolumeBound(), weight);
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::ChiralViolationContribs::addContrib
        let num_points = self.owner_points.expect("no owner");
        assert!(cset.idx1 < num_points);
        assert!(cset.idx2 < num_points);
        assert!(cset.idx3 < num_points);
        assert!(cset.idx4 < num_points);
        self.contribs.push(ChiralViolationContribsParams::new(
            cset.idx1,
            cset.idx2,
            cset.idx3,
            cset.idx4,
            cset.get_upper_volume_bound(),
            cset.get_lower_volume_bound(),
            weight,
        ));
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

    #[must_use]
    pub fn contribs(&self) -> &[ChiralViolationContribsParams] {
        &self.contribs
    }
}

impl ForceFieldContribution for ChiralViolationContribs {
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️: ChiralViolationContribs *copy() const override {
        // RDKit✔️✔️:   return new ChiralViolationContribs(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        let pos = context.coordinates();
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getEnergy (ChiralViolationContribs.cpp:78-94)
        // RDKit✔️✔️: double ChiralViolationContribs::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const unsigned int dim = dp_forceField->dimension();
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   for (const auto &c : d_contribs) {
        // RDKit✔️✔️:     double vol = calcChiralVolume(c.idx1, c.idx2, c.idx3, c.idx4, pos, dim);
        // RDKit✔️✔️:     if (vol < c.volLower) {
        // RDKit✔️✔️:       res += c.weight * (vol - c.volLower) * (vol - c.volLower);
        // RDKit✔️✔️:     } else if (vol > c.volUpper) {
        // RDKit✔️✔️:       res += c.weight * (vol - c.volUpper) * (vol - c.volUpper);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getEnergy
        assert!(!pos.is_empty(), "bad vector");
        let dim = self.owner_dimension.expect("no owner");
        assert_eq!(
            dim,
            context.dimension() as usize,
            "force field has wrong dimension"
        );
        let mut res = 0.0;
        for c in &self.contribs {
            let vol = calc_chiral_volume_flat(c.idx1, c.idx2, c.idx3, c.idx4, pos, dim);
            if vol < c.vol_lower {
                res += c.weight * (vol - c.vol_lower) * (vol - c.vol_lower);
            } else if vol > c.vol_upper {
                res += c.weight * (vol - c.vol_upper) * (vol - c.vol_upper);
            }
        }
        Ok(res)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        grad: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        let pos = context.coordinates();
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getGrad (ChiralViolationContribs.cpp:96-174)
        // RDKit✔️✔️: void ChiralViolationContribs::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const unsigned int dim = dp_forceField->dimension();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (const auto &c : d_contribs) {
        // RDKit✔️✔️:     // even if we are minimizing in higher dimension the chiral volume is
        // RDKit✔️✔️:     // calculated using only the first 3 dimensions
        // RDKit✔️✔️:     RDGeom::Point3D v1(pos[c.idx1 * dim] - pos[c.idx4 * dim],
        // RDKit✔️✔️:                        pos[c.idx1 * dim + 1] - pos[c.idx4 * dim + 1],
        // RDKit✔️✔️:                        pos[c.idx1 * dim + 2] - pos[c.idx4 * dim + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D v2(pos[c.idx2 * dim] - pos[c.idx4 * dim],
        // RDKit✔️✔️:                        pos[c.idx2 * dim + 1] - pos[c.idx4 * dim + 1],
        // RDKit✔️✔️:                        pos[c.idx2 * dim + 2] - pos[c.idx4 * dim + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D v3(pos[c.idx3 * dim] - pos[c.idx4 * dim],
        // RDKit✔️✔️:                        pos[c.idx3 * dim + 1] - pos[c.idx4 * dim + 1],
        // RDKit✔️✔️:                        pos[c.idx3 * dim + 2] - pos[c.idx4 * dim + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
        // RDKit✔️✔️:     double vol = v1.dotProduct(v2xv3);
        // RDKit✔️✔️:     double preFactor;
        // RDKit✔️✔️:     if (vol < c.volLower) {
        // RDKit✔️✔️:       preFactor = c.weight * (vol - c.volLower);
        // RDKit✔️✔️:     } else if (vol > c.volUpper) {
        // RDKit✔️✔️:       preFactor = c.weight * (vol - c.volUpper);
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        // END RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getGrad
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        let dim = self.owner_dimension.expect("no owner");
        assert_eq!(
            dim,
            context.dimension() as usize,
            "force field has wrong dimension"
        );

        for c in &self.contribs {
            let v1 = point3(
                pos[c.idx1 * dim] - pos[c.idx4 * dim],
                pos[c.idx1 * dim + 1] - pos[c.idx4 * dim + 1],
                pos[c.idx1 * dim + 2] - pos[c.idx4 * dim + 2],
            );
            let v2 = point3(
                pos[c.idx2 * dim] - pos[c.idx4 * dim],
                pos[c.idx2 * dim + 1] - pos[c.idx4 * dim + 1],
                pos[c.idx2 * dim + 2] - pos[c.idx4 * dim + 2],
            );
            let v3 = point3(
                pos[c.idx3 * dim] - pos[c.idx4 * dim],
                pos[c.idx3 * dim + 1] - pos[c.idx4 * dim + 1],
                pos[c.idx3 * dim + 2] - pos[c.idx4 * dim + 2],
            );

            let vol = v1.dot_product(&v2.cross_product(&v3));
            let pre_factor = if vol < c.vol_lower {
                c.weight * (vol - c.vol_lower)
            } else if vol > c.vol_upper {
                c.weight * (vol - c.vol_upper)
            } else {
                continue;
            };

            // RDKit✔️✔️:     grad[dim * c.idx1] += preFactor * ((v2.y) * (v3.z) - (v3.y) * (v2.z));
            // RDKit✔️✔️:     grad[dim * c.idx1 + 1] += preFactor * ((v3.x) * (v2.z) - (v2.x) * (v3.z));
            // RDKit✔️✔️:     grad[dim * c.idx1 + 2] += preFactor * ((v2.x) * (v3.y) - (v3.x) * (v2.y));
            grad[dim * c.idx1] += pre_factor * (v2.y * v3.z - v3.y * v2.z);
            grad[dim * c.idx1 + 1] += pre_factor * (v3.x * v2.z - v2.x * v3.z);
            grad[dim * c.idx1 + 2] += pre_factor * (v2.x * v3.y - v3.x * v2.y);

            // RDKit✔️✔️:     grad[dim * c.idx2] += preFactor * ((v3.y) * (v1.z) - (v3.z) * (v1.y));
            // RDKit✔️✔️:     grad[dim * c.idx2 + 1] += preFactor * ((v3.z) * (v1.x) - (v3.x) * (v1.z));
            // RDKit✔️✔️:     grad[dim * c.idx2 + 2] += preFactor * ((v3.x) * (v1.y) - (v3.y) * (v1.x));
            grad[dim * c.idx2] += pre_factor * (v3.y * v1.z - v3.z * v1.y);
            grad[dim * c.idx2 + 1] += pre_factor * (v3.z * v1.x - v3.x * v1.z);
            grad[dim * c.idx2 + 2] += pre_factor * (v3.x * v1.y - v3.y * v1.x);

            // RDKit✔️✔️:     grad[dim * c.idx3] += preFactor * ((v2.z) * (v1.y) - (v2.y) * (v1.z));
            // RDKit✔️✔️:     grad[dim * c.idx3 + 1] += preFactor * ((v2.x) * (v1.z) - (v2.z) * (v1.x));
            // RDKit✔️✔️:     grad[dim * c.idx3 + 2] += preFactor * ((v2.y) * (v1.x) - (v2.x) * (v1.y));
            grad[dim * c.idx3] += pre_factor * (v2.z * v1.y - v2.y * v1.z);
            grad[dim * c.idx3 + 1] += pre_factor * (v2.x * v1.z - v2.z * v1.x);
            grad[dim * c.idx3 + 2] += pre_factor * (v2.y * v1.x - v2.x * v1.y);

            // RDKit✔️✔️:     grad[dim * c.idx4] +=
            // RDKit✔️✔️:         preFactor * (pos[c.idx1 * dim + 2] *
            // RDKit✔️✔️:                          (pos[c.idx2 * dim + 1] - pos[c.idx3 * dim + 1]) +
            // RDKit✔️✔️:                      pos[c.idx2 * dim + 2] *
            // RDKit✔️✔️:                          (pos[c.idx3 * dim + 1] - pos[c.idx1 * dim + 1]) +
            // RDKit✔️✔️:                      pos[c.idx3 * dim + 2] *
            // RDKit✔️✔️:                          (pos[c.idx1 * dim + 1] - pos[c.idx2 * dim + 1]));
            grad[dim * c.idx4] += pre_factor
                * (pos[c.idx1 * dim + 2] * (pos[c.idx2 * dim + 1] - pos[c.idx3 * dim + 1])
                    + pos[c.idx2 * dim + 2] * (pos[c.idx3 * dim + 1] - pos[c.idx1 * dim + 1])
                    + pos[c.idx3 * dim + 2] * (pos[c.idx1 * dim + 1] - pos[c.idx2 * dim + 1]));

            // RDKit✔️✔️:     grad[dim * c.idx4 + 1] +=
            // RDKit✔️✔️:         preFactor *
            // RDKit✔️✔️:         (pos[c.idx1 * dim] * (pos[c.idx2 * dim + 2] - pos[c.idx3 * dim + 2]) +
            // RDKit✔️✔️:          pos[c.idx2 * dim] * (pos[c.idx3 * dim + 2] - pos[c.idx1 * dim + 2]) +
            // RDKit✔️✔️:          pos[c.idx3 * dim] * (pos[c.idx1 * dim + 2] - pos[c.idx2 * dim + 2]));
            grad[dim * c.idx4 + 1] += pre_factor
                * (pos[c.idx1 * dim] * (pos[c.idx2 * dim + 2] - pos[c.idx3 * dim + 2])
                    + pos[c.idx2 * dim] * (pos[c.idx3 * dim + 2] - pos[c.idx1 * dim + 2])
                    + pos[c.idx3 * dim] * (pos[c.idx1 * dim + 2] - pos[c.idx2 * dim + 2]));

            // RDKit✔️✔️:     grad[dim * c.idx4 + 2] +=
            // RDKit✔️✔️:         preFactor *
            // RDKit✔️✔️:         (pos[c.idx1 * dim + 1] * (pos[c.idx2 * dim] - pos[c.idx3 * dim]) +
            // RDKit✔️✔️:          pos[c.idx2 * dim + 1] * (pos[c.idx3 * dim] - pos[c.idx1 * dim]) +
            // RDKit✔️✔️:          pos[c.idx3 * dim + 1] * (pos[c.idx1 * dim] - pos[c.idx2 * dim]));
            grad[dim * c.idx4 + 2] += pre_factor
                * (pos[c.idx1 * dim + 1] * (pos[c.idx2 * dim] - pos[c.idx3 * dim])
                    + pos[c.idx2 * dim + 1] * (pos[c.idx3 * dim] - pos[c.idx1 * dim])
                    + pos[c.idx3 * dim + 1] * (pos[c.idx1 * dim] - pos[c.idx2 * dim]));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DistViolationContribsParams {
    pub idx1: usize,
    pub idx2: usize,
    pub ub: f64,
    pub lb: f64,
    pub ub2: f64,
    pub lb2: f64,
    pub weight: f64,
}

impl DistViolationContribsParams {
    #[must_use]
    pub fn new(idx1: usize, idx2: usize, ub: f64, lb: f64, weight: f64) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribsParams (DistViolationContribs.h:18-30)
        // RDKit✔️✔️: DistViolationContribsParams(unsigned int i1, unsigned int i2, double u,
        // RDKit✔️✔️:                             double l, double w = 1.0)
        // RDKit✔️✔️:     : idx1(i1), idx2(i2), ub(u), lb(l), ub2(u * u), lb2(l * l), weight(w) {};
        // END RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribsParams
        Self {
            idx1,
            idx2,
            ub,
            lb,
            ub2: ub * ub,
            lb2: lb * lb,
            weight,
        }
    }
}

// BEGIN RDKIT CPP LOCAL HELPER DistGeom::distance2 (DistViolationContribs.cpp:21-31)
// RDKit✔️✔️: inline double distance2(const unsigned int idx1, const unsigned int idx2,
// RDKit✔️✔️:                         const double *pos, const unsigned int dim) {
// RDKit✔️✔️:   const auto *end1Coords = &(pos[dim * idx1]);
// RDKit✔️✔️:   const auto *end2Coords = &(pos[dim * idx2]);
// RDKit✔️✔️:   double d2 = 0.0;
// RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
// RDKit✔️✔️:     double d = end1Coords[i] - end2Coords[i];
// RDKit✔️✔️:     d2 += d * d;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return d2;
// RDKit✔️✔️: }
// END RDKIT CPP LOCAL HELPER DistGeom::distance2
fn dist_violation_distance2(idx1: usize, idx2: usize, pos: &[f64], dim: usize) -> f64 {
    let mut d2 = 0.0;
    for i in 0..dim {
        let d = pos[dim * idx1 + i] - pos[dim * idx2 + i];
        d2 += d * d;
    }
    d2
}

// BEGIN RDKIT CPP LOCAL HELPER DistGeom::distance (DistViolationContribs.cpp:33-36)
// RDKit✔️✔️: inline double distance(const unsigned int idx1, const unsigned int idx2,
// RDKit✔️✔️:                        const double *pos, const unsigned int dim) {
// RDKit✔️✔️:   return sqrt(distance2(idx1, idx2, pos, dim));
// RDKit✔️✔️: }
// END RDKIT CPP LOCAL HELPER DistGeom::distance
fn dist_violation_distance(idx1: usize, idx2: usize, pos: &[f64], dim: usize) -> f64 {
    dist_violation_distance2(idx1, idx2, pos, dim).sqrt()
}

#[derive(Debug, Clone, PartialEq)]
pub struct DistViolationContribs {
    owner_dimension: Option<usize>,
    owner_points: Option<usize>,
    contribs: Vec<DistViolationContribsParams>,
}

impl Default for DistViolationContribs {
    fn default() -> Self {
        Self {
            owner_dimension: None,
            owner_points: None,
            contribs: Vec::new(),
        }
    }
}

impl DistViolationContribs {
    #[must_use]
    pub fn new(owner: &ForceField<'_>) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribs::DistViolationContribs (DistViolationContribs.cpp:16-19)
        // RDKit✔️✔️: DistViolationContribs::DistViolationContribs(ForceFields::ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribs::DistViolationContribs
        Self {
            owner_dimension: Some(owner.dimension() as usize),
            owner_points: Some(owner.positions().len()),
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, idx1: usize, idx2: usize, ub: f64, lb: f64, weight: f64) {
        // BEGIN RDKIT CPP METHOD DistGeom::DistViolationContribs::addContrib (DistViolationContribs.h:49-52)
        // RDKit✔️✔️: void addContrib(unsigned int idx1, unsigned int idx2, double ub, double lb,
        // RDKit✔️✔️:                 double weight = 1.0) {
        // RDKit✔️✔️:   d_contribs.emplace_back(idx1, idx2, ub, lb, weight);
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::DistViolationContribs::addContrib
        self.contribs
            .push(DistViolationContribsParams::new(idx1, idx2, ub, lb, weight));
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

    #[must_use]
    pub fn contribs(&self) -> &[DistViolationContribsParams] {
        &self.contribs
    }
}

impl ForceFieldContribution for DistViolationContribs {
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️: DistViolationContribs *copy() const override {
        // RDKit✔️✔️:   return new DistViolationContribs(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        let pos = context.coordinates();
        // BEGIN RDKIT CPP METHOD DistGeom::DistViolationContribs::getEnergy (DistViolationContribs.cpp:38-59)
        // RDKit✔️✔️: double DistViolationContribs::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   double accum = 0.0;
        // RDKit✔️✔️:   auto contrib = [&](const auto &c) {
        // RDKit✔️✔️:     double d2 = distance2(c.idx1, c.idx2, pos, dp_forceField->dimension());
        // RDKit✔️✔️:     double val = 0.0;
        // RDKit✔️✔️:     if (d2 > c.ub2) {
        // RDKit✔️✔️:       val = (d2 / (c.ub2)) - 1.0;
        // RDKit✔️✔️:     } else if (d2 < c.lb2) {
        // RDKit✔️✔️:       val = ((2 * c.lb2) / (c.lb2 + d2)) - 1.0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (val > 0.0) {
        // RDKit✔️✔️:       accum += c.weight * val * val;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   };
        // RDKit✔️✔️:   for (const auto &c : d_contribs) {
        // RDKit✔️✔️:     contrib(c);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return accum;
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::DistViolationContribs::getEnergy
        assert!(!pos.is_empty(), "bad vector");
        let dim = self.owner_dimension.expect("no owner");
        assert_eq!(
            dim,
            context.dimension() as usize,
            "force field has wrong dimension"
        );
        let mut accum = 0.0;
        for c in &self.contribs {
            let d2 = dist_violation_distance2(c.idx1, c.idx2, pos, dim);
            let mut val = 0.0;
            if d2 > c.ub2 {
                val = d2 / c.ub2 - 1.0;
            } else if d2 < c.lb2 {
                val = (2.0 * c.lb2) / (c.lb2 + d2) - 1.0;
            }
            if val > 0.0 {
                accum += c.weight * val * val;
            }
        }
        Ok(accum)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        grad: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        let pos = context.coordinates();
        // BEGIN RDKIT CPP METHOD DistGeom::DistViolationContribs::getGrad (DistViolationContribs.cpp:61-96)
        // RDKit✔️✔️: void DistViolationContribs::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:   const unsigned int dim = this->dp_forceField->dimension();
        // RDKit✔️✔️:   auto contrib = [&](const auto &c) {
        // RDKit✔️✔️:     double d2 = distance2(c.idx1, c.idx2, pos, dp_forceField->dimension());
        // RDKit✔️✔️:     double d;
        // RDKit✔️✔️:     double preFactor = 0.0;
        // RDKit✔️✔️:     if (d2 > c.ub2) {
        // RDKit✔️✔️:       d = sqrt(d2);
        // RDKit✔️✔️:       preFactor = 4. * (((d * d) / c.ub2) - 1.0) * (d / c.ub2);
        // RDKit✔️✔️:     } else if (d2 < c.lb2) {
        // RDKit✔️✔️:       d = sqrt(d2);
        // RDKit✔️✔️:       double l2d2 = d2 + c.lb2;
        // RDKit✔️✔️:       preFactor = 8. * c.lb2 * d * (1. - 2 * c.lb2 / l2d2) / (l2d2 * l2d2);
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // END RDKIT CPP METHOD DistGeom::DistViolationContribs::getGrad
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        let dim = self.owner_dimension.expect("no owner");
        assert_eq!(
            dim,
            context.dimension() as usize,
            "force field has wrong dimension"
        );
        for c in &self.contribs {
            let d2 = dist_violation_distance2(c.idx1, c.idx2, pos, dim);
            let d;
            let pre_factor;
            if d2 > c.ub2 {
                d = dist_violation_distance(c.idx1, c.idx2, pos, dim);
                pre_factor = 4.0 * ((d * d) / c.ub2 - 1.0) * (d / c.ub2);
            } else if d2 < c.lb2 {
                d = dist_violation_distance(c.idx1, c.idx2, pos, dim);
                let l2d2 = d2 + c.lb2;
                pre_factor = 8.0 * c.lb2 * d * (1.0 - 2.0 * c.lb2 / l2d2) / (l2d2 * l2d2);
            } else {
                continue;
            }
            // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
            // RDKit✔️✔️:       const auto p1 = dim * c.idx1 + i;
            // RDKit✔️✔️:       const auto p2 = dim * c.idx2 + i;
            // RDKit✔️✔️:       double dGrad;
            // RDKit✔️✔️:       if (d > 0.0) {
            // RDKit✔️✔️:         dGrad = c.weight * preFactor * (pos[p1] - pos[p2]) / d;
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         // FIX: this likely isn't right
            // RDKit✔️✔️:         dGrad = c.weight * preFactor * (pos[p1] - pos[p2]);
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       grad[p1] += dGrad;
            // RDKit✔️✔️:       grad[p2] -= dGrad;
            // RDKit✔️✔️:     }
            for i in 0..dim {
                let p1 = dim * c.idx1 + i;
                let p2 = dim * c.idx2 + i;
                let d_grad = if d > 0.0 {
                    c.weight * pre_factor * (pos[p1] - pos[p2]) / d
                } else {
                    c.weight * pre_factor * (pos[p1] - pos[p2])
                };
                grad[p1] += d_grad;
                grad[p2] -= d_grad;
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FourthDimContribsParams {
    pub idx: usize,
    pub weight: f64,
}

impl FourthDimContribsParams {
    #[must_use]
    pub fn new(idx: usize, weight: f64) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribsParams (FourthDimContribs.h:19-23)
        // RDKit✔️✔️: struct FourthDimContribsParams {
        // RDKit✔️✔️:   unsigned int idx{0};
        // RDKit✔️✔️:   double weight{0.0};
        // RDKit✔️✔️:   FourthDimContribsParams(unsigned int idx, double w) : idx(idx), weight(w) {};
        // RDKit✔️✔️: };
        // END RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribsParams
        Self { idx, weight }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FourthDimContribs {
    owner_dimension: Option<usize>,
    owner_points: Option<usize>,
    contribs: Vec<FourthDimContribsParams>,
}

impl Default for FourthDimContribs {
    fn default() -> Self {
        // RDKit✔️✔️: FourthDimContribs() = default;
        Self {
            owner_dimension: None,
            owner_points: None,
            contribs: Vec::new(),
        }
    }
}

impl FourthDimContribs {
    #[must_use]
    pub fn new(owner: &ForceField<'_>) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribs::FourthDimContribs (FourthDimContribs.h:36-40)
        // RDKit✔️✔️: FourthDimContribs(ForceFields::ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad force field");
        // RDKit✔️✔️:   PRECONDITION(owner->dimension() == 4, "force field has wrong dimension");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribs::FourthDimContribs
        assert_eq!(owner.dimension(), 4, "force field has wrong dimension");
        Self {
            owner_dimension: Some(owner.dimension() as usize),
            owner_points: Some(owner.positions().len()),
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, idx: usize, weight: f64) {
        // BEGIN RDKIT CPP METHOD DistGeom::FourthDimContribs::addContrib (FourthDimContribs.h:42-44)
        // RDKit✔️✔️: void addContrib(unsigned int idx, double weight) {
        // RDKit✔️✔️:   d_contribs.emplace_back(idx, weight);
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::FourthDimContribs::addContrib
        self.contribs
            .push(FourthDimContribsParams::new(idx, weight));
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

    #[must_use]
    pub fn contribs(&self) -> &[FourthDimContribsParams] {
        &self.contribs
    }
}

impl ForceFieldContribution for FourthDimContribs {
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️: FourthDimContribs *copy() const override {
        // RDKit✔️✔️:   return new FourthDimContribs(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        let pos = context.coordinates();
        // BEGIN RDKIT CPP METHOD DistGeom::FourthDimContribs::getEnergy (FourthDimContribs.h:47-57)
        // RDKit✔️✔️: double getEnergy(double *pos) const override {
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   constexpr unsigned int ffdim = 4;
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        // RDKit✔️✔️:     unsigned int pid = contrib.idx * ffdim + 3;
        // RDKit✔️✔️:     res += contrib.weight * pos[pid] * pos[pid];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::FourthDimContribs::getEnergy
        assert!(!pos.is_empty(), "bad vector");
        assert_eq!(
            self.owner_dimension.expect("no owner"),
            context.dimension() as usize,
            "force field has wrong dimension"
        );
        let ffdim = 4;
        let mut res = 0.0;
        for contrib in &self.contribs {
            let pid = contrib.idx * ffdim + 3;
            res += contrib.weight * pos[pid] * pos[pid];
        }
        Ok(res)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        grad: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        let pos = context.coordinates();
        // BEGIN RDKIT CPP METHOD DistGeom::FourthDimContribs::getGrad (FourthDimContribs.h:61-70)
        // RDKit✔️✔️: void getGrad(double *pos, double *grad) const override {
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   constexpr unsigned int ffdim = 4;
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        // RDKit✔️✔️:     unsigned int pid = contrib.idx * ffdim + 3;
        // RDKit✔️✔️:     grad[pid] += contrib.weight * pos[pid];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::FourthDimContribs::getGrad
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        assert_eq!(
            self.owner_dimension.expect("no owner"),
            context.dimension() as usize,
            "force field has wrong dimension"
        );
        let ffdim = 4;
        for contrib in &self.contribs {
            let pid = contrib.idx * ffdim + 3;
            grad[pid] += contrib.weight * pos[pid];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_field<'a>(dimension: u32, rows: &'a mut [Vec<f64>]) -> ForceField<'a> {
        let mut field = ForceField::new(dimension);
        field
            .positions_mut()
            .extend(rows.iter_mut().map(Vec::as_mut_slice));
        field
    }
    fn evaluate_energy(
        term: &dyn ForceFieldContribution,
        coordinates: &[f64],
        dimension: u32,
    ) -> f64 {
        let n = coordinates.len() / dimension as usize;
        let mut cache = vec![-1.; n * (n + 1) / 2];
        let mut context =
            EvaluationContext::for_distgeom_test(coordinates, &mut cache, n as u32, dimension);
        term.get_energy(&mut context).unwrap()
    }
    fn evaluate_grad(
        term: &dyn ForceFieldContribution,
        coordinates: &[f64],
        gradient: &mut [f64],
        dimension: u32,
    ) {
        let n = coordinates.len() / dimension as usize;
        let mut cache = vec![-1.; n * (n + 1) / 2];
        let mut context =
            EvaluationContext::for_distgeom_test(coordinates, &mut cache, n as u32, dimension);
        term.get_grad(&mut context, gradient).unwrap();
    }
    #[test]
    fn distgeom_chiral_set_constructor_stores_indices_bounds_and_flags() {
        let chiral_set = ChiralSet::new(
            10,
            11,
            12,
            13,
            14,
            -2.5,
            3.5,
            ChiralSetStructureFlags::InFusedSmallRings as u64,
        );

        assert_eq!(chiral_set.idx0, 10);
        assert_eq!(chiral_set.idx1, 11);
        assert_eq!(chiral_set.idx2, 12);
        assert_eq!(chiral_set.idx3, 13);
        assert_eq!(chiral_set.idx4, 14);
        assert_eq!(chiral_set.volume_lower_bound, -2.5);
        assert_eq!(chiral_set.volume_upper_bound, 3.5);
        assert_eq!(
            chiral_set.structure_flags,
            ChiralSetStructureFlags::InFusedSmallRings as u64
        );
    }

    #[test]
    fn distgeom_chiral_set_default_structure_flags_match_rdkit_default_argument() {
        let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

        assert_eq!(chiral_set.structure_flags, 0);
    }

    #[test]
    fn distgeom_chiral_set_getters_return_volume_bounds() {
        let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -7.0, -3.0);

        assert_eq!(chiral_set.get_lower_volume_bound(), -7.0);
        assert_eq!(chiral_set.get_upper_volume_bound(), -3.0);
    }

    #[test]
    fn distgeom_chiral_set_allows_equal_volume_bounds() {
        let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, 2.0, 2.0);

        assert_eq!(chiral_set.get_lower_volume_bound(), 2.0);
        assert_eq!(chiral_set.get_upper_volume_bound(), 2.0);
    }

    #[test]
    #[should_panic(expected = "Inconsistent bounds")]
    fn distgeom_chiral_set_rejects_lower_bound_above_upper_bound() {
        let _ = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, 2.0, 1.0);
    }

    #[test]
    fn distgeom_chiral_set_aliases_model_shared_pointer_vector() {
        let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, -1.0, 1.0,
        ));
        let chiral_sets: VectChiralSet = vec![Arc::clone(&chiral_set)];

        assert_eq!(chiral_sets.len(), 1);
        assert!(Arc::ptr_eq(&chiral_set, &chiral_sets[0]));
    }

    #[test]
    fn chiral_volume_flat_returns_signed_triple_product() {
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 0.0,
        ];

        assert_eq!(calc_chiral_volume_flat(0, 1, 2, 3, &pos, 3), 1.0);
        assert_eq!(calc_chiral_volume_flat(0, 2, 1, 3, &pos, 3), -1.0);
    }

    #[test]
    fn chiral_volume_flat_uses_idx4_as_reference_point() {
        let pos = [
            2.0, 2.0, 2.0, //
            3.0, 2.0, 2.0, //
            2.0, 3.0, 2.0, //
            2.0, 2.0, 3.0,
        ];

        assert_eq!(calc_chiral_volume_flat(1, 2, 3, 0, &pos, 3), 1.0);
    }

    #[test]
    fn chiral_volume_flat_ignores_dimensions_after_first_three() {
        let pos = [
            1.0, 0.0, 0.0, 100.0, //
            0.0, 1.0, 0.0, 200.0, //
            0.0, 0.0, 1.0, 300.0, //
            0.0, 0.0, 0.0, 400.0,
        ];

        assert_eq!(calc_chiral_volume_flat(0, 1, 2, 3, &pos, 4), 1.0);
    }

    #[test]
    fn chiral_volume_points_returns_signed_triple_product() {
        let pts = [
            point3(1.0, 0.0, 0.0),
            point3(0.0, 1.0, 0.0),
            point3(0.0, 0.0, 1.0),
            point3(0.0, 0.0, 0.0),
        ];

        assert_eq!(calc_chiral_volume_points(0, 1, 2, 3, &pts), 1.0);
        assert_eq!(calc_chiral_volume_points(0, 2, 1, 3, &pts), -1.0);
    }

    fn chiral_violation_pos() -> Vec<f64> {
        vec![
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            0.0, 0.0, 0.0,
        ]
    }

    #[test]
    fn chiral_violation_contribs_constructor_starts_empty() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let contribs = ChiralViolationContribs::new(&ff);

        assert!(contribs.empty());
        assert_eq!(contribs.size(), 0);
    }

    #[test]
    fn chiral_violation_contribs_add_contrib_copies_chiral_set_bounds_and_weight() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = ChiralViolationContribs::new(&ff);
        let cset = ChiralSet::with_default_structure_flags(99, 0, 1, 2, 3, -0.5, 0.5);

        contribs.add_contrib(&cset, 2.5);

        assert!(!contribs.empty());
        assert_eq!(contribs.size(), 1);
        assert_eq!(
            contribs.contribs()[0],
            ChiralViolationContribsParams::new(0, 1, 2, 3, 0.5, -0.5, 2.5)
        );
    }

    #[test]
    #[should_panic]
    fn chiral_violation_contribs_add_contrib_rejects_out_of_range_indices() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = ChiralViolationContribs::new(&ff);
        let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 4, -0.5, 0.5);

        contribs.add_contrib(&cset, 1.0);
    }

    #[test]
    fn chiral_violation_contribs_get_energy_returns_zero_inside_bounds() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = ChiralViolationContribs::new(&ff);
        let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, 0.5, 1.5);
        contribs.add_contrib(&cset, 2.0);

        assert_eq!(evaluate_energy(&contribs, &chiral_violation_pos(), 3), 0.0);
    }

    #[test]
    fn chiral_violation_contribs_get_energy_accumulates_lower_and_upper_violations() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = ChiralViolationContribs::new(&ff);
        let upper = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, -0.5, 0.5);
        let lower = ChiralSet::with_default_structure_flags(0, 0, 2, 1, 3, -0.5, 0.5);
        contribs.add_contrib(&upper, 2.0);
        contribs.add_contrib(&lower, 3.0);

        assert_eq!(
            evaluate_energy(&contribs, &chiral_violation_pos(), 3),
            2.0 * 0.5 * 0.5 + 3.0 * 0.5 * 0.5
        );
    }

    #[test]
    fn chiral_violation_contribs_get_grad_returns_early_inside_bounds() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = ChiralViolationContribs::new(&ff);
        let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, 0.5, 1.5);
        contribs.add_contrib(&cset, 2.0);
        let mut grad = vec![10.0; 12];

        evaluate_grad(&contribs, &chiral_violation_pos(), &mut grad, 3);

        assert_eq!(grad, vec![10.0; 12]);
    }

    #[test]
    fn chiral_violation_contribs_get_grad_matches_source_formula_for_upper_violation() {
        let mut fixture_rows = vec![
            vec![1., 0., 0.],
            vec![0., 1., 0.],
            vec![0., 0., 1.],
            vec![0., 0., 0.],
        ];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = ChiralViolationContribs::new(&ff);
        let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, -1.0, 0.0);
        contribs.add_contrib(&cset, 2.0);
        let mut grad = vec![0.0; 12];

        evaluate_grad(&contribs, &chiral_violation_pos(), &mut grad, 3);

        assert_eq!(
            grad,
            vec![
                2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, -2.0, -2.0, -2.0
            ]
        );
    }

    fn dist_violation_pos(distance: f64) -> Vec<f64> {
        vec![0.0, 0.0, 0.0, distance, 0.0, 0.0]
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "actual={actual} expected={expected}"
        );
    }

    #[test]
    fn dist_violation_contribs_constructor_and_add_contrib_store_squared_bounds() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        assert!(contribs.empty());

        contribs.add_contrib(0, 1, 3.0, 1.5, 2.0);

        assert_eq!(contribs.size(), 1);
        assert_eq!(
            contribs.contribs()[0],
            DistViolationContribsParams::new(0, 1, 3.0, 1.5, 2.0)
        );
        assert_eq!(contribs.contribs()[0].ub2, 9.0);
        assert_eq!(contribs.contribs()[0].lb2, 2.25);
    }

    #[test]
    fn dist_violation_contribs_distance_helpers_follow_source_dim_loop() {
        let pos = [0.0, 0.0, 0.0, 5.0, 1.0, 2.0, 2.0, 9.0];

        assert_eq!(dist_violation_distance2(0, 1, &pos, 4), 25.0);
        assert_eq!(dist_violation_distance(0, 1, &pos, 4), 25.0_f64.sqrt());
    }

    #[test]
    fn dist_violation_contribs_get_energy_returns_zero_inside_bounds() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        contribs.add_contrib(0, 1, 3.0, 1.0, 2.0);

        assert_eq!(evaluate_energy(&contribs, &dist_violation_pos(2.0), 3), 0.0);
    }

    #[test]
    fn dist_violation_contribs_get_energy_matches_upper_violation_formula() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        contribs.add_contrib(0, 1, 1.0, 0.0, 2.0);

        assert_eq!(
            evaluate_energy(&contribs, &dist_violation_pos(2.0), 3),
            18.0
        );
    }

    #[test]
    fn dist_violation_contribs_get_energy_matches_lower_violation_formula() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        contribs.add_contrib(0, 1, 10.0, 2.0, 3.0);

        assert_close(
            evaluate_energy(&contribs, &dist_violation_pos(1.0), 3),
            1.08,
        );
    }

    #[test]
    fn dist_violation_contribs_get_grad_returns_early_inside_bounds() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        contribs.add_contrib(0, 1, 3.0, 1.0, 2.0);
        let mut grad = vec![7.0; 6];

        evaluate_grad(&contribs, &dist_violation_pos(2.0), &mut grad, 3);

        assert_eq!(grad, vec![7.0; 6]);
    }

    #[test]
    fn dist_violation_contribs_get_grad_matches_upper_violation_formula() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        contribs.add_contrib(0, 1, 1.0, 0.0, 2.0);
        let mut grad = vec![0.0; 6];

        evaluate_grad(&contribs, &dist_violation_pos(2.0), &mut grad, 3);

        assert_eq!(grad, vec![-48.0, 0.0, 0.0, 48.0, 0.0, 0.0]);
    }

    #[test]
    fn dist_violation_contribs_get_grad_matches_lower_violation_formula() {
        let mut fixture_rows = vec![vec![0., 0., 0.], vec![2., 0., 0.]];
        let ff = fixture_field(3, &mut fixture_rows);
        let mut contribs = DistViolationContribs::new(&ff);
        contribs.add_contrib(0, 1, 10.0, 2.0, 3.0);
        let mut grad = vec![0.0; 6];

        evaluate_grad(&contribs, &dist_violation_pos(1.0), &mut grad, 3);

        assert_close(grad[0], 2.304);
        assert_eq!(grad[1], 0.0);
        assert_eq!(grad[2], 0.0);
        assert_close(grad[3], -2.304);
        assert_eq!(grad[4], 0.0);
        assert_eq!(grad[5], 0.0);
    }

    fn fourth_dim_pos() -> Vec<f64> {
        vec![
            0.0, 0.0, 0.0, 2.0, //
            1.0, 2.0, 3.0, -3.0,
        ]
    }

    #[test]
    fn fourth_dim_contribs_default_starts_without_owner_and_empty() {
        let contribs = FourthDimContribs::default();

        assert!(contribs.empty());
        assert_eq!(contribs.size(), 0);
    }

    #[test]
    fn fourth_dim_contribs_constructor_requires_four_dimensional_forcefield() {
        let mut fixture_rows = vec![vec![0., 0., 0., 0.], vec![1., 0., 0., 0.]];
        let ff = fixture_field(4, &mut fixture_rows);
        let contribs = FourthDimContribs::new(&ff);

        assert!(contribs.empty());
        assert_eq!(contribs.size(), 0);
    }

    #[test]
    #[should_panic(expected = "force field has wrong dimension")]
    fn fourth_dim_contribs_constructor_rejects_non_four_dimensional_forcefield() {
        let ff = ForceField::new(3);

        let _ = FourthDimContribs::new(&ff);
    }

    #[test]
    fn fourth_dim_contribs_add_contrib_appends_index_and_weight() {
        let mut fixture_rows = vec![vec![0., 0., 0., 0.], vec![1., 0., 0., 0.]];
        let ff = fixture_field(4, &mut fixture_rows);
        let mut contribs = FourthDimContribs::new(&ff);

        contribs.add_contrib(1, 2.5);

        assert!(!contribs.empty());
        assert_eq!(contribs.size(), 1);
        assert_eq!(contribs.contribs()[0], FourthDimContribsParams::new(1, 2.5));
    }

    #[test]
    fn fourth_dim_contribs_get_energy_accumulates_weighted_fourth_coordinate_squares() {
        let mut fixture_rows = vec![vec![0., 0., 0., 0.], vec![1., 0., 0., 0.]];
        let ff = fixture_field(4, &mut fixture_rows);
        let mut contribs = FourthDimContribs::new(&ff);
        contribs.add_contrib(0, 2.0);
        contribs.add_contrib(1, 3.0);

        assert_eq!(evaluate_energy(&contribs, &fourth_dim_pos(), 4), 35.0);
    }

    #[test]
    fn fourth_dim_contribs_get_grad_adds_source_weighted_fourth_coordinate_terms() {
        let mut fixture_rows = vec![vec![0., 0., 0., 0.], vec![1., 0., 0., 0.]];
        let ff = fixture_field(4, &mut fixture_rows);
        let mut contribs = FourthDimContribs::new(&ff);
        contribs.add_contrib(0, 2.0);
        contribs.add_contrib(1, 3.0);
        let mut grad = vec![10.0; 8];

        evaluate_grad(&contribs, &fourth_dim_pos(), &mut grad, 4);

        assert_eq!(grad, vec![10.0, 10.0, 10.0, 14.0, 10.0, 10.0, 10.0, 1.0]);
    }

    #[test]
    fn fourth_dim_contribs_copy_preserves_contribs_and_behavior() {
        let mut fixture_rows = vec![vec![0., 0., 0., 0.], vec![1., 0., 0., 0.]];
        let ff = fixture_field(4, &mut fixture_rows);
        let mut contribs = FourthDimContribs::new(&ff);
        contribs.add_contrib(1, 3.0);
        let copied = contribs.copy();

        assert_eq!(evaluate_energy(&*copied, &fourth_dim_pos(), 4), 27.0);
    }
}

mod construction;

pub use construction::{
    ConformerOptimizer, ConformerOptimizerError, DistanceBoundsRead,
    DistanceGeometryForceFieldParams, calc_chiral_volume_rows,
};
