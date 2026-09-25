//! Source-shaped private primitives used by the RDKit 2D layout owner.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use cosmolkit_model::TopologyBlock;

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
        // RDKit❗✔️: void Transform2D::SetTransform(const Point2D &pt, double angle) {
        // RDKit❗✔️:   this->setToIdentity();
        // RDKit❗✔️:   Transform2D trans1;
        // RDKit❗✔️:   trans1.SetTranslation(-pt);
        // RDKit❗✔️:   double *data = d_data.get();
        // RDKit❗✔️:   data[0] = cos(angle);
        // RDKit❗✔️:   data[1] = -sin(angle);
        // RDKit❗✔️:   data[3] = sin(angle);
        // RDKit❗✔️:   data[4] = cos(angle);
        // RDKit❗✔️:   (*this) *= trans1;
        // RDKit❗✔️:   Transform2D trans2;
        // RDKit❗✔️:   trans2.SetTranslation(pt);
        // RDKit❗✔️:   trans2 *= (*this);
        // RDKit❗✔️:   this->assign(trans2);
        // RDKit❗✔️: }
        let (sin, cos) = angle.sin_cos();
        Self {
            data: [
                cos,
                -sin,
                point[0] - cos * point[0] + sin * point[1],
                sin,
                cos,
                point[1] - sin * point[0] - cos * point[1],
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
        // RDKit❗✔️: void Transform2D::SetTransform(const Point2D &ref1, const Point2D &ref2,
        // RDKit❗✔️:                                const Point2D &pt1, const Point2D &pt2) {
        // RDKit❗✔️:   Point2D rvec = ref2 - ref1;
        // RDKit❗✔️:   Point2D pvec = pt2 - pt1;
        // RDKit❗✔️:   double dp = rvec.dotProduct(pvec);
        // RDKit❗✔️:   double lp = (rvec.length()) * (pvec.length());
        // RDKit❗✔️:   if (lp <= 0.0) {
        // RDKit❗✔️:     this->setToIdentity();
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   double cval = dp / lp;
        // RDKit❗✔️:   if (cval < -1.0) {
        // RDKit❗✔️:     cval = -1.0;
        // RDKit❗✔️:   } else if (cval > 1.0) {
        // RDKit❗✔️:     cval = 1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   double ang = acos(cval);
        // RDKit❗✔️:   double cross = (pvec.x) * (rvec.y) - (pvec.y) * (rvec.x);
        // RDKit❗✔️:   if (cross < 0.0) {
        // RDKit❗✔️:     ang *= -1.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   this->setToIdentity();
        // RDKit❗✔️:   double *data = d_data.get();
        // RDKit❗✔️:   data[0] = cos(ang);
        // RDKit❗✔️:   data[1] = -sin(ang);
        // RDKit❗✔️:   data[3] = sin(ang);
        // RDKit❗✔️:   data[4] = cos(ang);
        // RDKit❗✔️:   Point2D npt1 = pt1;
        // RDKit❗✔️:   this->TransformPoint(npt1);
        // RDKit❗✔️:   data[DIM_2D - 1] = ref1.x - npt1.x;
        // RDKit❗✔️:   data[2 * DIM_2D - 1] = ref1.y - npt1.y;
        // RDKit❗✔️: }
        let rvec = [ref2[0] - ref1[0], ref2[1] - ref1[1]];
        let pvec = [pt2[0] - pt1[0], pt2[1] - pt1[1]];
        let dp = rvec[0] * pvec[0] + rvec[1] * pvec[1];
        let lp = rvec[0].hypot(rvec[1]) * pvec[0].hypot(pvec[1]);
        if lp <= 0.0 {
            return Self::identity();
        }
        let cval = (dp / lp).clamp(-1.0, 1.0);
        let mut ang = cval.acos();
        let cross = pvec[0] * rvec[1] - pvec[1] * rvec[0];
        if cross < 0.0 {
            ang *= -1.0;
        }
        let mut transform = Self {
            data: [ang.cos(), -ang.sin(), 0.0, ang.sin(), ang.cos(), 0.0],
        };
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
            rank = text
                .parse()
                .map_err(|_| GeometryError::InvalidRankProperty {
                    atom,
                    key: "_CIPRank",
                })?;
        } else {
            if let Some(text) = value.prop("_ChiralAtomRank") {
                let chiral_rank: u32 =
                    text.parse()
                        .map_err(|_| GeometryError::InvalidRankProperty {
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
