//! RDKit-aligned atom chiral-tag assignment from detached 3D coordinates.
//!
//! This module owns only the geometry kernel. Live molecule checkout, cache
//! effects, `_StereochemDone` removal, and commit remain in `cosmolkit`.

use std::f64::consts::PI;

use cosmolkit_model::{
    AtomId, AtomPropertyError, Bond, BondDirection, BondOrder, ChiralTag, CoordinateBlock,
    TopologyBlock, TopologyValidationError,
};

use crate::{StereoOrderError, ValenceAssignment, bond_affects_atom_chirality};

// RDKit✔️✔️: static constexpr double zero_tolerance = 1.e-16;
const ZERO_VECTOR_TOLERANCE: f64 = 1.0e-16;
const ZERO_VOLUME_TOLERANCE: f64 = 0.1;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StructureTagParams {
    pub conformer_id: i32,
    pub replace_existing_tags: bool,
}

impl Default for StructureTagParams {
    fn default() -> Self {
        Self {
            conformer_id: -1,
            replace_existing_tags: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StructureTagAssignment {
    pub topology: TopologyBlock,
    pub selected_conformer_id: Option<usize>,
    pub clear_stereochem_done: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("duplicate 3D conformer id {id}")]
    DuplicateConformerId { id: usize },
    #[error("cannot find conformer with id {requested}")]
    ConformerNotFound { requested: i32 },
    #[error("3D conformer {conformer} has {rows} coordinate rows, expected {atom_count}")]
    ConformerAtomCountMismatch {
        conformer: usize,
        rows: usize,
        atom_count: usize,
    },
    #[error("valence field {field} has {actual} rows, expected {atom_count}")]
    InvalidValence {
        field: &'static str,
        actual: usize,
        atom_count: usize,
    },
    #[error("valence field {field} has invalid value {value} for atom {atom}")]
    InvalidValenceValue {
        field: &'static str,
        atom: AtomId,
        value: i32,
    },
    #[error("cannot normalize zero-length vector from atom {center} to atom {neighbor}")]
    ZeroLengthVector { center: AtomId, neighbor: AtomId },
    #[error("stereo-order input is invalid: {0}")]
    StereoOrder(#[from] StereoOrderError),
    #[error("atom property update failed: {0}")]
    AtomProperty(#[from] AtomPropertyError),
}

#[derive(Clone, Copy)]
struct Vec3([f64; 3]);

impl Vec3 {
    fn between(from: [f64; 3], to: [f64; 3]) -> Self {
        // BEGIN RDKIT CPP FUNCTION Point3D operator-
        // RDKit✔️✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = p1.x - p2.x;
        // RDKit✔️✔️:   res.y = p1.y - p2.y;
        // RDKit✔️✔️:   res.z = p1.z - p2.z;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Point3D operator-
        Self([to[0] - from[0], to[1] - from[1], to[2] - from[2]])
    }

    fn normalized_between(
        from: [f64; 3],
        to: [f64; 3],
        center: AtomId,
        neighbor: AtomId,
    ) -> Result<Self, StereoError> {
        // BEGIN RDKIT CPP FUNCTION Point3D::directionVector/normalize
        // RDKit✔️✔️:   Point3D directionVector(const Point3D &other) const {
        // RDKit✔️✔️:     Point3D res;
        // RDKit✔️✔️:     res.x = other.x - x;
        // RDKit✔️✔️:     res.y = other.y - y;
        // RDKit✔️✔️:     res.z = other.z - z;
        // RDKit✔️✔️:     res.normalize();
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   constexpr void normalize() override {
        // RDKit✔️✔️:     double l = this->length();
        // RDKit✔️✔️:     if (l < zero_tolerance) {
        // RDKit✔️✔️:       throw std::runtime_error("Cannot normalize a zero length vector");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     x /= l;
        // RDKit✔️✔️:     y /= l;
        // RDKit✔️✔️:     z /= l;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   double length() const override {
        // RDKit✔️✔️:     double res = x * x + y * y + z * z;
        // RDKit✔️✔️:     return sqrt(res);
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION Point3D::directionVector/normalize
        let mut value = Self::between(from, to);
        let length =
            (value.0[0] * value.0[0] + value.0[1] * value.0[1] + value.0[2] * value.0[2]).sqrt();
        if length < ZERO_VECTOR_TOLERANCE {
            return Err(StereoError::ZeroLengthVector { center, neighbor });
        }
        value.0[0] /= length;
        value.0[1] /= length;
        value.0[2] /= length;
        Ok(value)
    }

    fn dot(self, other: Self) -> f64 {
        // BEGIN RDKIT CPP FUNCTION Point3D::dotProduct
        // RDKit✔️✔️:   constexpr double dotProduct(const Point3D &other) const {
        // RDKit✔️✔️:     double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION Point3D::dotProduct
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    fn cross(self, other: Self) -> Self {
        // BEGIN RDKIT CPP FUNCTION Point3D::crossProduct
        // RDKit✔️✔️:   constexpr Point3D crossProduct(const Point3D &other) const {
        // RDKit✔️✔️:     Point3D res;
        // RDKit✔️✔️:     res.x = y * (other.z) - z * (other.y);
        // RDKit✔️✔️:     res.y = -x * (other.z) + z * (other.x);
        // RDKit✔️✔️:     res.z = x * (other.y) - y * (other.x);
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION Point3D::crossProduct
        Self([
            self.0[1] * other.0[2] - self.0[2] * other.0[1],
            -self.0[0] * other.0[2] + self.0[2] * other.0[0],
            self.0[0] * other.0[1] - self.0[1] * other.0[0],
        ])
    }

    fn angle_to(self, other: Self) -> f64 {
        // BEGIN RDKIT CPP FUNCTION Point3D::angleTo
        // RDKit✔️✔️:   double angleTo(const Point3D &other) const {
        // RDKit✔️✔️:     double lsq = lengthSq() * other.lengthSq();
        // RDKit✔️✔️:     double dotProd = dotProduct(other);
        // RDKit✔️✔️:     dotProd /= sqrt(lsq);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // watch for roundoff error:
        // RDKit✔️✔️:     if (dotProd <= -1.0) {
        // RDKit✔️✔️:       return M_PI;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (dotProd >= 1.0) {
        // RDKit✔️✔️:       return 0.0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     return acos(dotProd);
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION Point3D::angleTo
        let self_sq = self.dot(self);
        let other_sq = other.dot(other);
        let dot = self.dot(other) / (self_sq * other_sq).sqrt();
        if dot <= -1.0 {
            PI
        } else if dot >= 1.0 {
            0.0
        } else {
            dot.acos()
        }
    }
}

fn volume_test(vectors: &[Vec3; 6], x: usize, y: usize, z: usize) -> bool {
    // RDKit✔️✔️: #define VOLTEST(X, Y, Z) (v[X].dotProduct(v[Y].crossProduct(v[Z])) >= 0.0)
    vectors[x].dot(vectors[y].cross(vectors[z])) >= 0.0
}

fn octahedral_permutation(pair: &[u8; 6], vectors: &[Vec3; 6]) -> u32 {
    // BEGIN RDKIT CPP FUNCTION OctahedralPermFrom3D
    // RDKit✔️✔️: static unsigned int OctahedralPermFrom3D(unsigned char *pair,
    // RDKit✔️✔️:                                          const RDGeom::Point3D *v) {
    // RDKit✔️✔️:   switch (pair[0]) {
    // RDKit✔️✔️:     case 2:  // a-b
    // RDKit✔️✔️:       switch (pair[2]) {
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           return VOLTEST(0, 3, 4) ? 28 : 27;
    // RDKit✔️✔️:         case 5:
    // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 25 : 30;
    // RDKit✔️✔️:         default:  // 0 or 6
    // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 26 : 29;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 3:  // a-c
    // RDKit✔️✔️:       switch (pair[1]) {
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           return VOLTEST(0, 3, 4) ? 22 : 21;
    // RDKit✔️✔️:         case 5:
    // RDKit✔️✔️:           return VOLTEST(0, 1, 3) ? 19 : 24;
    // RDKit✔️✔️:         default:  // 0 or 6
    // RDKit✔️✔️:           return VOLTEST(0, 1, 3) ? 20 : 23;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 4:  // a-d
    // RDKit✔️✔️:       switch (pair[1]) {
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           return VOLTEST(0, 2, 4) ? 13 : 12;
    // RDKit✔️✔️:         case 5:
    // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 6 : 18;
    // RDKit✔️✔️:         default:  // 0 or 6
    // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 7 : 17;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 5:  // a-e
    // RDKit✔️✔️:       switch (pair[1]) {
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 11 : 9;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 3 : 16;
    // RDKit✔️✔️:         default:  // 0 or 6
    // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 5 : 15;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:  // 0 or 6  a-f
    // RDKit✔️✔️:       switch (pair[1]) {
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 10 : 8;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 1 : 2;
    // RDKit✔️✔️:         default:  // 5
    // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 4 : 14;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // unreachable
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION OctahedralPermFrom3D
    match pair[0] {
        2 => match pair[2] {
            4 => {
                if volume_test(vectors, 0, 3, 4) {
                    28
                } else {
                    27
                }
            }
            5 => {
                if volume_test(vectors, 0, 2, 3) {
                    25
                } else {
                    30
                }
            }
            _ => {
                if volume_test(vectors, 0, 2, 3) {
                    26
                } else {
                    29
                }
            }
        },
        3 => match pair[1] {
            4 => {
                if volume_test(vectors, 0, 3, 4) {
                    22
                } else {
                    21
                }
            }
            5 => {
                if volume_test(vectors, 0, 1, 3) {
                    19
                } else {
                    24
                }
            }
            _ => {
                if volume_test(vectors, 0, 1, 3) {
                    20
                } else {
                    23
                }
            }
        },
        4 => match pair[1] {
            3 => {
                if volume_test(vectors, 0, 2, 4) {
                    13
                } else {
                    12
                }
            }
            5 => {
                if volume_test(vectors, 0, 1, 2) {
                    6
                } else {
                    18
                }
            }
            _ => {
                if volume_test(vectors, 0, 1, 2) {
                    7
                } else {
                    17
                }
            }
        },
        5 => match pair[1] {
            3 => {
                if volume_test(vectors, 0, 2, 3) {
                    11
                } else {
                    9
                }
            }
            4 => {
                if volume_test(vectors, 0, 1, 2) {
                    3
                } else {
                    16
                }
            }
            _ => {
                if volume_test(vectors, 0, 1, 2) {
                    5
                } else {
                    15
                }
            }
        },
        _ => match pair[1] {
            3 => {
                if volume_test(vectors, 0, 2, 3) {
                    10
                } else {
                    8
                }
            }
            4 => {
                if volume_test(vectors, 0, 1, 2) {
                    1
                } else {
                    2
                }
            }
            _ => {
                if volume_test(vectors, 0, 1, 2) {
                    4
                } else {
                    14
                }
            }
        },
    }
}

fn is_wiggly_bond(bond: &Bond, center: AtomId) -> bool {
    // BEGIN RDKIT CPP FUNCTION isWigglyBond
    // RDKit✔️✔️: bool isWigglyBond(const Bond *bond, const Atom *atom) {
    // RDKit✔️✔️:   int hasWigglyBond = 0;
    // RDKit✔️✔️:   if (bond->getBeginAtomIdx() == atom->getIdx() &&
    // RDKit✔️✔️:       bond->getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:       (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
    // RDKit✔️✔️:        (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
    // RDKit✔️✔️:                                     hasWigglyBond) &&
    // RDKit✔️✔️:         hasWigglyBond))) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isWigglyBond
    bond.begin() == center
        && bond.order() == BondOrder::Single
        && (bond.direction() == BondDirection::Unknown || bond.unknown_stereo())
}

fn non_tetrahedral_assignment(
    topology: &TopologyBlock,
    positions: &[[f64; 3]],
    center: AtomId,
) -> Result<Option<(ChiralTag, u32)>, StereoError> {
    // BEGIN RDKIT CPP FUNCTION assignNontetrahedralChiralTypeFrom3D
    // RDKit✔️✔️: static bool assignNontetrahedralChiralTypeFrom3D(ROMol &mol,
    // RDKit✔️✔️:                                                  const Conformer &conf,
    // RDKit✔️✔️:                                                  Atom *atom,
    // RDKit✔️✔️:                                                  double tolerance = 0.1) {
    // RDKit✔️✔️:   // FIX: add tests for dative and zero order bonds
    // RDKit✔️✔️:   // Fail fast check for non-tetrahedral elements
    // RDKit✔️✔️:   if (atom->getAtomicNum() < 15) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check for wiggly bonds
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (isWigglyBond(bond, atom)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D cen = conf.getAtomPos(atom->getIdx());
    // RDKit✔️✔️:   RDGeom::Point3D v[6];
    // RDKit✔️✔️:   unsigned int count = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     if (count == 6) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     RDGeom::Point3D p = conf.getAtomPos(*nbrIdx);
    // RDKit✔️✔️:     v[count] = cen.directionVector(p);
    // RDKit✔️✔️:     ++count;
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (count < 3) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned char pair[6];
    // RDKit✔️✔️:   memset(pair, 0, 6);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int pairs = 0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < count; j++) {
    // RDKit✔️✔️:       if (v[i].dotProduct(v[j]) < -(1 - tolerance)) {
    // RDKit✔️✔️:         if (pair[i] || pair[j]) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         pair[i] = j + 1;
    // RDKit✔️✔️:         pair[j] = i + 1;
    // RDKit✔️✔️:         pairs++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Atom::ChiralType tag;
    // RDKit✔️✔️:   unsigned int perm;
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   switch (pairs) {
    // RDKit✔️✔️:     case 0:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 1:
    // RDKit✔️✔️:       switch (count) {
    // RDKit✔️✔️:         case 3: /* T-shape */
    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           if (pair[0] == 0) {
    // RDKit✔️✔️:             perm = 3;  // Z
    // RDKit✔️✔️:           } else if (pair[0] == 2) {
    // RDKit✔️✔️:             perm = 2;  // 4
    // RDKit✔️✔️:           } else /* pair[0] == 3 */ {
    // RDKit✔️✔️:             perm = 1;  // U
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:                /* See-saw */
    // RDKit✔️✔️:           if (pair[0] == 2) {  // a b
    // RDKit✔️✔️:             if (v[2].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 2, 3) ? 25 : 29;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 2, 3) ? 7 : 8;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[0] == 3) {  // a c
    // RDKit✔️✔️:             if (v[1].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 19 : 23;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 5 : 6;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[0] == 4) {  // a d
    // RDKit✔️✔️:             if (v[1].angleTo(v[2]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 2) ? 6 : 17;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 2) ? 3 : 4;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[1] == 3) {  // b c
    // RDKit✔️✔️:             if (v[0].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 10 : 8;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(1, 0, 3) ? 13 : 14;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[1] == 4) {  // b d
    // RDKit✔️✔️:             if (v[0].angleTo(v[2]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 1 : 2;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(1, 0, 2) ? 10 : 12;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else /* pair[2] == 4 */ {  // c d
    // RDKit✔️✔️:             if (v[0].angleTo(v[1]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 4 : 14;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(3, 0, 1) ? 16 : 19;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setChiralTag(tag);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 5: /* Trigonal bipyramidal */
    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           if (pair[0] == 2) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 2, 3) ? 7 : 8;  // a b
    // RDKit✔️✔️:           } else if (pair[0] == 3) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 3) ? 5 : 6;  // a c
    // RDKit✔️✔️:           } else if (pair[0] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 2) ? 3 : 4;  // a d
    // RDKit✔️✔️:           } else if (pair[0] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 2) ? 1 : 2;  // a e
    // RDKit✔️✔️:           } else if (pair[1] == 3) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 3) ? 13 : 14;  // b c
    // RDKit✔️✔️:           } else if (pair[1] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 2) ? 10 : 12;  // b d
    // RDKit✔️✔️:           } else if (pair[1] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 2) ? 9 : 11;  // b e
    // RDKit✔️✔️:           } else if (pair[2] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(2, 0, 1) ? 16 : 19;  // c d
    // RDKit✔️✔️:           } else if (pair[2] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(2, 0, 1) ? 15 : 20;  // c e
    // RDKit✔️✔️:           } else /* pair[2] == 4 */ {
    // RDKit✔️✔️:             perm = VOLTEST(3, 0, 1) ? 17 : 18;  // d e
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 2:
    // RDKit✔️✔️:       if (count == 4) {
    // RDKit✔️✔️:         /* Square planar */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         if (pair[0] == 2) {
    // RDKit✔️✔️:           perm = 2;  // 4
    // RDKit✔️✔️:         } else if (pair[0] == 3) {
    // RDKit✔️✔️:           perm = 1;  // U
    // RDKit✔️✔️:         } else /* pair[1] == 4 */ {
    // RDKit✔️✔️:           perm = 3;  // Z
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       } else if (count == 5) {
    // RDKit✔️✔️:         /* Square pyramidal */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         perm = OctahedralPermFrom3D(pair, v);
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 3:
    // RDKit✔️✔️:       if (count == 6) {
    // RDKit✔️✔️:         /* Octahedral */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         perm = OctahedralPermFrom3D(pair, v);
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION assignNontetrahedralChiralTypeFrom3D
    let atom = &topology.atoms[center.index()];
    if atom.atomic_number() < 15 {
        return Ok(None);
    }
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        if is_wiggly_bond(&topology.bonds[neighbor.bond.index()], center) {
            return Ok(None);
        }
    }
    let center_position = positions[center.index()];
    let neighbors = topology.adjacency.neighbors_of(center.index());
    if neighbors.len() > 6 || neighbors.len() < 3 {
        return Ok(None);
    }
    let mut vectors = [Vec3([0.0; 3]); 6];
    for (index, neighbor) in neighbors.iter().enumerate() {
        let neighbor_id = AtomId::new(neighbor.atom_index);
        vectors[index] = Vec3::normalized_between(
            center_position,
            positions[neighbor.atom_index],
            center,
            neighbor_id,
        )?;
    }
    let count = neighbors.len();
    let mut pair = [0u8; 6];
    let mut pairs = 0;
    for i in 0..count {
        for j in i + 1..count {
            if vectors[i].dot(vectors[j]) < -0.9 {
                if pair[i] != 0 || pair[j] != 0 {
                    return Ok(None);
                }
                pair[i] = (j + 1) as u8;
                pair[j] = (i + 1) as u8;
                pairs += 1;
            }
        }
    }

    let result = match pairs {
        1 if count == 3 => {
            let permutation = if pair[0] == 0 {
                3
            } else if pair[0] == 2 {
                2
            } else {
                1
            };
            Some((ChiralTag::SquarePlanar, permutation))
        }
        1 if count == 4 => {
            let threshold = 100.0 * PI / 180.0;
            let (tag, permutation) = if pair[0] == 2 {
                if vectors[2].angle_to(vectors[3]) < threshold {
                    (
                        ChiralTag::Octahedral,
                        if volume_test(&vectors, 0, 2, 3) {
                            25
                        } else {
                            29
                        },
                    )
                } else {
                    (
                        ChiralTag::TrigonalBipyramidal,
                        if volume_test(&vectors, 0, 2, 3) { 7 } else { 8 },
                    )
                }
            } else if pair[0] == 3 {
                if vectors[1].angle_to(vectors[3]) < threshold {
                    (
                        ChiralTag::Octahedral,
                        if volume_test(&vectors, 0, 1, 3) {
                            19
                        } else {
                            23
                        },
                    )
                } else {
                    (
                        ChiralTag::TrigonalBipyramidal,
                        if volume_test(&vectors, 0, 1, 3) { 5 } else { 6 },
                    )
                }
            } else if pair[0] == 4 {
                if vectors[1].angle_to(vectors[2]) < threshold {
                    (
                        ChiralTag::Octahedral,
                        if volume_test(&vectors, 0, 1, 2) {
                            6
                        } else {
                            17
                        },
                    )
                } else {
                    (
                        ChiralTag::TrigonalBipyramidal,
                        if volume_test(&vectors, 0, 1, 2) { 3 } else { 4 },
                    )
                }
            } else if pair[1] == 3 {
                if vectors[0].angle_to(vectors[3]) < threshold {
                    (
                        ChiralTag::Octahedral,
                        if volume_test(&vectors, 0, 1, 3) {
                            10
                        } else {
                            8
                        },
                    )
                } else {
                    (
                        ChiralTag::TrigonalBipyramidal,
                        if volume_test(&vectors, 1, 0, 3) {
                            13
                        } else {
                            14
                        },
                    )
                }
            } else if pair[1] == 4 {
                if vectors[0].angle_to(vectors[2]) < threshold {
                    (
                        ChiralTag::Octahedral,
                        if volume_test(&vectors, 0, 1, 3) { 1 } else { 2 },
                    )
                } else {
                    (
                        ChiralTag::TrigonalBipyramidal,
                        if volume_test(&vectors, 1, 0, 2) {
                            10
                        } else {
                            12
                        },
                    )
                }
            } else if vectors[0].angle_to(vectors[1]) < threshold {
                (
                    ChiralTag::Octahedral,
                    if volume_test(&vectors, 0, 1, 3) {
                        4
                    } else {
                        14
                    },
                )
            } else {
                (
                    ChiralTag::TrigonalBipyramidal,
                    if volume_test(&vectors, 3, 0, 1) {
                        16
                    } else {
                        19
                    },
                )
            };
            Some((tag, permutation))
        }
        1 if count == 5 => {
            let permutation = if pair[0] == 2 {
                if volume_test(&vectors, 0, 2, 3) { 7 } else { 8 }
            } else if pair[0] == 3 {
                if volume_test(&vectors, 0, 1, 3) { 5 } else { 6 }
            } else if pair[0] == 4 {
                if volume_test(&vectors, 0, 1, 2) { 3 } else { 4 }
            } else if pair[0] == 5 {
                if volume_test(&vectors, 0, 1, 2) { 1 } else { 2 }
            } else if pair[1] == 3 {
                if volume_test(&vectors, 1, 0, 3) {
                    13
                } else {
                    14
                }
            } else if pair[1] == 4 {
                if volume_test(&vectors, 1, 0, 2) {
                    10
                } else {
                    12
                }
            } else if pair[1] == 5 {
                if volume_test(&vectors, 1, 0, 2) {
                    9
                } else {
                    11
                }
            } else if pair[2] == 4 {
                if volume_test(&vectors, 2, 0, 1) {
                    16
                } else {
                    19
                }
            } else if pair[2] == 5 {
                if volume_test(&vectors, 2, 0, 1) {
                    15
                } else {
                    20
                }
            } else if volume_test(&vectors, 3, 0, 1) {
                17
            } else {
                18
            };
            Some((ChiralTag::TrigonalBipyramidal, permutation))
        }
        2 if count == 4 => {
            let permutation = if pair[0] == 2 {
                2
            } else if pair[0] == 3 {
                1
            } else {
                3
            };
            Some((ChiralTag::SquarePlanar, permutation))
        }
        2 if count == 5 => Some((
            ChiralTag::Octahedral,
            octahedral_permutation(&pair, &vectors),
        )),
        3 if count == 6 => Some((
            ChiralTag::Octahedral,
            octahedral_permutation(&pair, &vectors),
        )),
        _ => None,
    };
    Ok(result)
}

fn validate_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), StereoError> {
    let atom_count = topology.atoms.len();
    for (field, values) in [
        ("explicit_valence", &valence.explicit_valence),
        ("implicit_hydrogens", &valence.implicit_hydrogens),
    ] {
        if values.len() != atom_count {
            return Err(StereoError::InvalidValence {
                field,
                actual: values.len(),
                atom_count,
            });
        }
        if let Some((index, value)) = values
            .iter()
            .copied()
            .enumerate()
            .find(|(_, value)| *value < 0)
        {
            return Err(StereoError::InvalidValenceValue {
                field,
                atom: AtomId::new(index),
                value,
            });
        }
    }
    Ok(())
}

fn total_hydrogens(topology: &TopologyBlock, valence: &ValenceAssignment, center: AtomId) -> usize {
    // BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalNumHs
    let atom = &topology.atoms[center.index()];
    let implicit = if atom.no_implicit() {
        0
    } else {
        valence.implicit_hydrogens[center.index()] as usize
    };
    usize::from(atom.explicit_hydrogens()) + implicit
}

fn nonzero_degree(topology: &TopologyBlock, center: AtomId) -> Result<usize, StereoError> {
    // BEGIN RDKIT CPP FUNCTION ROMol adjacency accessors
    // RDKit✔️✔️:   CXXAtomIterator<const MolGraph, Atom *const, MolGraph::adjacency_iterator>
    // RDKit✔️✔️:   atomNeighbors(Atom const *at) const {
    // RDKit✔️✔️:     auto pr = getAtomNeighbors(at);
    // RDKit✔️✔️:     return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CXXAtomIterator<MolGraph, Atom *, MolGraph::adjacency_iterator> atomNeighbors(
    // RDKit✔️✔️:       Atom const *at) {
    // RDKit✔️✔️:     auto pr = getAtomNeighbors(at);
    // RDKit✔️✔️:     return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CXXBondIterator<const MolGraph, Bond *const, MolGraph::out_edge_iterator>
    // RDKit✔️✔️:   atomBonds(Atom const *at) const {
    // RDKit✔️✔️:     auto pr = getAtomBonds(at);
    // RDKit✔️✔️:     return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CXXBondIterator<MolGraph, Bond *, MolGraph::out_edge_iterator> atomBonds(
    // RDKit✔️✔️:       Atom const *at) {
    // RDKit✔️✔️:     auto pr = getAtomBonds(at);
    // RDKit✔️✔️:     return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
    // RDKit✔️✔️: };
    // RDKit✔️✔️:
    // RDKit✔️✔️: ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return boost::out_edges(at->getIdx(), d_graph);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol adjacency accessors
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::getAtomNonzeroDegree
    // RDKit✔️✔️: unsigned int getAtomNonzeroDegree(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad pointer");
    // RDKit✔️✔️:   PRECONDITION(atom->hasOwningMol(), "no owning molecule");
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (auto bond : atom->getOwningMol().atomBonds(atom)) {
    // RDKit✔️✔️:     if (!bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // END RDKIT CPP FUNCTION Chirality::detail::getAtomNonzeroDegree
    let mut degree = 0;
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        if bond_affects_atom_chirality(&topology.bonds[neighbor.bond.index()], center)? {
            degree += 1;
        }
    }
    Ok(degree)
}

fn nontetrahedral_enabled() -> bool {
    // BEGIN RDKIT CPP FUNCTION getValFromEnvironment/getAllowNontetrahedralChirality
    // RDKit✔️❌: constexpr auto nonTetrahedralStereoEnvVar = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";
    // RDKit✔️❌: constexpr bool nonTetrahedralStereoDefaultVal =
    // RDKit✔️❌:     true;  //!< whether or not nontetrahedral stereo is perceived by default
    // RDKit✔️❌: bool getValFromEnvironment(const char *var, bool defVal) {
    // RDKit✔️❌:   auto evar = std::getenv(var);
    // RDKit✔️❌:   if (evar != nullptr) {
    // RDKit✔️❌:     if (!strcmp(evar, "0")) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       return true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return defVal;
    // RDKit✔️❌: }
    // RDKit✔️❌: bool getAllowNontetrahedralChirality() {
    // RDKit✔️❌:   return getValFromEnvironment(nonTetrahedralStereoEnvVar,
    // RDKit✔️❌:                                nonTetrahedralStereoDefaultVal);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION getValFromEnvironment/getAllowNontetrahedralChirality
    std::env::var_os("RDK_ENABLE_NONTETRAHEDRAL_STEREO").is_none_or(|value| value != "0")
}

/// Assign atom chiral tags from the selected detached 3D conformer.
pub fn assign_chiral_tags_from_structure(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    valence: &ValenceAssignment,
    params: &StructureTagParams,
) -> Result<StructureTagAssignment, StereoError> {
    // BEGIN RDKIT CPP FUNCTION assignChiralTypesFrom3D
    // RDKit✔️❌: void assignChiralTypesFrom3D(ROMol &mol, int confId, bool replaceExistingTags) {
    // RDKit✔️❌:   const double ZERO_VOLUME_TOL = 0.1;
    // RDKit✔️❌:   if (!mol.getNumConformers()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️❌:   if (!conf.is3D()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // if the molecule already has stereochemistry
    // RDKit✔️❌:   // perceived, remove the flags that indicate
    // RDKit✔️❌:   // this... what we're about to do will require
    // RDKit✔️❌:   // that we go again.
    // RDKit✔️❌:   if (mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️❌:     mol.clearProp(common_properties::_StereochemDone);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   auto allowNontetrahedralStereo = Chirality::getAllowNontetrahedralChirality();
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> explicitAtoms;
    // RDKit✔️❌:   explicitAtoms.resize(mol.getNumAtoms(), 0);
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     auto bondDir = bond->getBondDir();
    // RDKit✔️❌:     if (bondDir == Bond::BondDir::BEGINWEDGE ||
    // RDKit✔️❌:         bondDir == Bond::BondDir::BEGINDASH) {
    // RDKit✔️❌:       explicitAtoms[bond->getBeginAtom()->getIdx()] = 1;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       explicitAtoms[atom->getIdx()] = 1;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     // if we aren't replacing existing tags and the atom is already tagged,
    // RDKit✔️❌:     // punt:
    // RDKit✔️❌:     if (!replaceExistingTags && atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:     // additional reasons to skip the atom:
    // RDKit✔️❌:     auto nzDegree = Chirality::detail::getAtomNonzeroDegree(atom);
    // RDKit✔️❌:     auto tnzDegree = nzDegree + atom->getTotalNumHs();
    // RDKit✔️❌:     if (nzDegree < 3 || tnzDegree > 6) {
    // RDKit✔️❌:       // not enough explicit neighbors or too many total neighbors
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (allowNontetrahedralStereo &&
    // RDKit✔️❌:         assignNontetrahedralChiralTypeFrom3D(mol, conf, atom)) {
    // RDKit✔️❌:       if (explicitAtoms[atom->getIdx()] == 0) {
    // RDKit✔️❌:         atom->setProp(common_properties::_NonExplicit3DChirality, 1);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     /* We're only doing tetrahedral cases here */
    // RDKit✔️❌:     if (tnzDegree > 4) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     int anum = atom->getAtomicNum();
    // RDKit✔️❌:     if (anum != 16 && anum != 34 &&  // S or Se are special
    // RDKit✔️❌:                                      // (just using the InChI list for now)
    // RDKit✔️❌:         tnzDegree != 4               // not enough total neighbors
    // RDKit✔️❌:     ) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     const auto &p0 = conf.getAtomPos(atom->getIdx());
    // RDKit✔️❌:     const RDGeom::Point3D *nbrs[4];
    // RDKit✔️❌:     unsigned int nbrIdx = 0;
    // RDKit✔️❌:     int hasWigglyBond = 0;
    // RDKit✔️❌:     for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️❌:       hasWigglyBond = isWigglyBond(bond, atom);
    // RDKit✔️❌:       if (hasWigglyBond) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (!Chirality::detail::bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       nbrs[nbrIdx++] = &conf.getAtomPos(bond->getOtherAtomIdx(atom->getIdx()));
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (hasWigglyBond) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     auto v1 = *nbrs[0] - p0;
    // RDKit✔️❌:     auto v2 = *nbrs[1] - p0;
    // RDKit✔️❌:     auto v3 = *nbrs[2] - p0;
    // RDKit✔️❌:
    // RDKit✔️❌:     double chiralVol = v1.dotProduct(v2.crossProduct(v3));
    // RDKit✔️❌:     bool chiralitySet = false;
    // RDKit✔️❌:     if (chiralVol < -ZERO_VOLUME_TOL) {
    // RDKit✔️❌:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️❌:       chiralitySet = true;
    // RDKit✔️❌:     } else if (chiralVol > ZERO_VOLUME_TOL) {
    // RDKit✔️❌:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️❌:       chiralitySet = true;
    // RDKit✔️❌:     } else if (nbrIdx == 4) {
    // RDKit✔️❌:       // The first three neighbors are on the same plane as the chiral atom (or
    // RDKit✔️❌:       // very close to it). If a 4th neighbor is present, let's see if this one
    // RDKit✔️❌:       // determines a chiral volume
    // RDKit✔️❌:
    // RDKit✔️❌:       auto v4 = *nbrs[3] - p0;
    // RDKit✔️❌:       // v4 would be in the opposite direction to v3
    // RDKit✔️❌:       chiralVol = -v1.dotProduct(v2.crossProduct(v4));
    // RDKit✔️❌:       if (chiralVol < -ZERO_VOLUME_TOL) {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️❌:         chiralitySet = true;
    // RDKit✔️❌:       } else if (chiralVol > ZERO_VOLUME_TOL) {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️❌:         chiralitySet = true;
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (chiralitySet && explicitAtoms[atom->getIdx()] == 0) {
    // RDKit✔️❌:       atom->setProp<int>(common_properties::_NonExplicit3DChirality, 1);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION assignChiralTypesFrom3D
    // BEGIN RDKIT CPP FUNCTION ROMol::getConformer
    // RDKit✔️✔️: const Conformer &ROMol::getConformer(int id) const {
    // RDKit✔️✔️:   // make sure we have more than one conformation
    // RDKit✔️✔️:   if (d_confs.size() == 0) {
    // RDKit✔️✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (id < 0) {
    // RDKit✔️✔️:     return *(d_confs.front());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto cid = (unsigned int)id;
    // RDKit✔️✔️:   for (auto conf : d_confs) {
    // RDKit✔️✔️:     if (conf->getId() == cid) {
    // RDKit✔️✔️:       return *conf;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // we did not find a conformation with the specified ID
    // RDKit✔️✔️:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit✔️✔️:   mesg += id;
    // RDKit✔️✔️:   throw ConformerException(mesg);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::getConformer
    // The detached assignment clones topology before computing changes so an
    // error has no committable partial result. That adds O(atoms + bonds +
    // properties) work versus RDKit's in-place mutation, hence ✔️❌.
    topology.validate()?;
    if coordinates.conformers_3d.is_empty() {
        return Ok(StructureTagAssignment {
            topology: topology.clone(),
            selected_conformer_id: None,
            clear_stereochem_done: false,
        });
    }
    let mut ids = std::collections::BTreeSet::new();
    for conformer in &coordinates.conformers_3d {
        if !ids.insert(conformer.id()) {
            return Err(StereoError::DuplicateConformerId { id: conformer.id() });
        }
    }
    let conformer = if params.conformer_id < 0 {
        &coordinates.conformers_3d[0]
    } else {
        coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == params.conformer_id as usize)
            .ok_or(StereoError::ConformerNotFound {
                requested: params.conformer_id,
            })?
    };
    if !conformer.is_3d() {
        return Ok(StructureTagAssignment {
            topology: topology.clone(),
            selected_conformer_id: Some(conformer.id()),
            clear_stereochem_done: false,
        });
    }
    if conformer.coordinates().len() != topology.atoms.len() {
        return Err(StereoError::ConformerAtomCountMismatch {
            conformer: conformer.id(),
            rows: conformer.coordinates().len(),
            atom_count: topology.atoms.len(),
        });
    }
    validate_valence(topology, valence)?;

    let positions = conformer.coordinates();
    let allow_nontetrahedral = nontetrahedral_enabled();
    let mut explicit = vec![false; topology.atoms.len()];
    for bond in &topology.bonds {
        if matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) {
            explicit[bond.begin().index()] = true;
        }
    }
    for atom in &topology.atoms {
        if atom.chiral_tag() != ChiralTag::Unspecified {
            explicit[atom.id().index()] = true;
        }
    }

    let mut working = topology.clone();
    for index in 0..working.atoms.len() {
        let center = AtomId::new(index);
        if !params.replace_existing_tags
            && working.atoms[index].chiral_tag() != ChiralTag::Unspecified
        {
            continue;
        }
        working.atoms[index].set_chiral_tag(ChiralTag::Unspecified);
        let degree = nonzero_degree(&working, center)?;
        let total_degree = degree + total_hydrogens(&working, valence, center);
        if degree < 3 || total_degree > 6 {
            continue;
        }
        if allow_nontetrahedral
            && let Some((tag, permutation)) =
                non_tetrahedral_assignment(&working, positions, center)?
        {
            working.atoms[index].set_chiral_tag(tag);
            working.atoms[index].set_chiral_permutation(Some(permutation));
            working.atoms[index].set_prop("_chiralPermutation", permutation.to_string())?;
            if !explicit[index] {
                working.atoms[index].set_prop("_NonExplicit3DChirality", "1")?;
            }
            continue;
        }
        if total_degree > 4 {
            continue;
        }
        let atomic_number = working.atoms[index].atomic_number();
        if atomic_number != 16 && atomic_number != 34 && total_degree != 4 {
            continue;
        }
        let mut neighbor_positions = [[0.0; 3]; 4];
        let mut neighbor_count = 0;
        let mut wiggly = false;
        for neighbor in working.adjacency.neighbors_of(index) {
            let bond = &working.bonds[neighbor.bond.index()];
            if is_wiggly_bond(bond, center) {
                wiggly = true;
                break;
            }
            if !bond_affects_atom_chirality(bond, center)? {
                continue;
            }
            neighbor_positions[neighbor_count] = positions[neighbor.atom_index];
            neighbor_count += 1;
        }
        if wiggly {
            continue;
        }
        let center_position = positions[index];
        let v1 = Vec3::between(center_position, neighbor_positions[0]);
        let v2 = Vec3::between(center_position, neighbor_positions[1]);
        let v3 = Vec3::between(center_position, neighbor_positions[2]);
        let mut volume = v1.dot(v2.cross(v3));
        let mut tag = if volume < -ZERO_VOLUME_TOLERANCE {
            Some(ChiralTag::TetrahedralCw)
        } else if volume > ZERO_VOLUME_TOLERANCE {
            Some(ChiralTag::TetrahedralCcw)
        } else {
            None
        };
        if tag.is_none() && neighbor_count == 4 {
            let v4 = Vec3::between(center_position, neighbor_positions[3]);
            volume = -v1.dot(v2.cross(v4));
            tag = if volume < -ZERO_VOLUME_TOLERANCE {
                Some(ChiralTag::TetrahedralCw)
            } else if volume > ZERO_VOLUME_TOLERANCE {
                Some(ChiralTag::TetrahedralCcw)
            } else {
                None
            };
        }
        working.atoms[index].set_chiral_tag(tag.unwrap_or(ChiralTag::Unspecified));
        if tag.is_some() && !explicit[index] {
            working.atoms[index].set_prop("_NonExplicit3DChirality", "1")?;
        }
    }

    Ok(StructureTagAssignment {
        topology: working,
        selected_conformer_id: Some(conformer.id()),
        clear_stereochem_done: true,
    })
}
