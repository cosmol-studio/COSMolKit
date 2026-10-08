use cosmolkit_model::{AtomId, BondId, Conformer3D, TopologyBlock};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag};

use crate::ValenceError;
use crate::valence::{
    source_atom_implicit_hydrogens, source_atom_needs_cache_update, update_source_atom_cache,
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BondDirectionStereoError {
    #[error("invalid stereochemistry state: {0}")]
    InvalidState(String),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

#[derive(Clone, Copy)]
struct Vector3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector3 {
    fn between(from: [f64; 3], to: [f64; 3]) -> Self {
        Self {
            x: to[0] - from[0],
            y: to[1] - from[1],
            z: to[2] - from[2],
        }
    }

    fn length_squared(self) -> f64 {
        self.dot(self)
    }

    fn length(self) -> f64 {
        self.length_squared().sqrt()
    }

    fn normalized(
        self,
        center: AtomId,
        neighbor: AtomId,
    ) -> Result<Self, BondDirectionStereoError> {
        let [x, y, z] = crate::structure_tags::normalize_vector_components(
            [self.x, self.y, self.z],
            center,
            neighbor,
        )
        .map_err(|error| BondDirectionStereoError::InvalidState(error.to_string()))?;
        Ok(Self { x, y, z })
    }

    fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    fn difference(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

fn pseudo_3d_chiral_tag(
    topology: &TopologyBlock,
    bond_id: BondId,
    conformer: &Conformer3D,
) -> Result<Option<ChiralTag>, BondDirectionStereoError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::atomChiralTypeFromBondDirPseudo3D complete source
    // RDKit❗❌: std::optional<Atom::ChiralType> atomChiralTypeFromBondDirPseudo3D(
    // RDKit❗❌:     const ROMol &mol, const Bond *bond, const Conformer *conf) {
    // RDKit❗❌:   PRECONDITION(bond, "no bond");
    // RDKit❗❌:   PRECONDITION(conf, "no conformer");
    // RDKit❗❌:   auto bondDir = bond->getBondDir();
    // RDKit❗❌:   PRECONDITION(bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH,
    // RDKit❗❌:                "bad bond direction");
    // RDKit❗❌:   constexpr double coordZeroTol = 1e-4;
    // RDKit❗❌:   constexpr double zeroTol = 1e-3;
    // RDKit❗❌:   constexpr double tShapeTol =
    // RDKit❗❌:       0.00031;  // used to recognize T-shaped arrangements
    // RDKit❗❌:   // corresponds to an angle between the two vectors of just under 178 degrees
    // RDKit❗❌:   // degree
    // RDKit❗❌:
    // RDKit❗❌:   constexpr double pseudo3DOffset = 0.1;  // z-displacement for wedged bonds
    // RDKit❗❌:
    // RDKit❗❌:   constexpr double volumeTolerance =
    // RDKit❗❌:       0.00174;  // used to recognize zero chiral volume
    // RDKit❗❌:   // This is what we get for a T-shaped arrangement with just over 178 degrees
    // RDKit❗❌:
    // RDKit❗❌:   // NOTE that according to the CT file spec, wedging assigns chirality
    // RDKit❗❌:   // to the atom at the point of the wedge, (atom 1 in the bond).
    // RDKit❗❌:   const auto atom = bond->getBeginAtom();
    // RDKit❗❌:   PRECONDITION(atom, "no atom");
    // RDKit❗❌:
    // RDKit❗❌:   // we can't do anything with atoms that have more than 4 neighbors:
    // RDKit❗❌:   if (atom->getDegree() > 4) {
    // RDKit❗❌:     return Atom::CHI_UNSPECIFIED;
    // RDKit❗❌:   }
    // RDKit❗❌:   const auto bondAtom = bond->getEndAtom();
    // RDKit❗❌:
    // RDKit❗❌:   Atom::ChiralType res = Atom::CHI_UNSPECIFIED;
    // RDKit❗❌:
    // RDKit❗❌:   auto centerLoc = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:   centerLoc.z = 0.0;
    // RDKit❗❌:   auto refPt = conf->getAtomPos(bondAtom->getIdx());
    // RDKit❗❌:
    // RDKit❗❌:   // Github #7305: in some odd cases, we get conformers with
    // RDKit❗❌:   // weird scalings. In these, we need to scale the 3d offset
    // RDKit❗❌:   // or it might be irrelevant or dominate over the coordinates.
    // RDKit❗❌:   auto refLength = (centerLoc - refPt).length();
    // RDKit❗❌:   refPt.z =
    // RDKit❗❌:       bondDir == Bond::BondDir::BEGINWEDGE ? pseudo3DOffset : -pseudo3DOffset;
    // RDKit❗❌:   if (refLength) {
    // RDKit❗❌:     refPt.z *= refLength;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  collect indices and bond vectors of neighbors and track whether or
    // RDKit❗❌:   //  not there's an H neighbor and if all bonds are single
    // RDKit❗❌:   //
    // RDKit❗❌:   //  at the end of this process bond 0 is the input wedged bond
    // RDKit❗❌:   //
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:
    // RDKit❗❌:   INT_VECT neighborBondIndices;
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int refIdx = mol.getNumBonds() + 1;
    // RDKit❗❌:   std::vector<RDGeom::Point3D> bondVects;
    // RDKit❗❌:   bool allSingle = true;
    // RDKit❗❌:   unsigned int nbrIdx = 0;
    // RDKit❗❌:   for (const auto nbrBond : mol.atomBonds(atom)) {
    // RDKit❗❌:     const auto oAtom = nbrBond->getOtherAtom(atom);
    // RDKit❗❌:     auto tmpPt = conf->getAtomPos(oAtom->getIdx());
    // RDKit❗❌:     if (nbrBond == bond) {
    // RDKit❗❌:       refIdx = nbrIdx;
    // RDKit❗❌:       tmpPt = refPt;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // theoretically we could confirm that this is a single bond,
    // RDKit❗❌:       // but it's not impossible that at some point in the future we
    // RDKit❗❌:       // could allow wedged multiple bonds for things like atropisomers
    // RDKit❗❌:       if (nbrBond->getBeginAtomIdx() == atom->getIdx() &&
    // RDKit❗❌:           (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit❗❌:            nbrBond->getBondDir() == Bond::BondDir::BEGINDASH)) {
    // RDKit❗❌:         // scale the 3d offset based on the reference bond here too
    // RDKit❗❌:         tmpPt.z = nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE
    // RDKit❗❌:                       ? pseudo3DOffset
    // RDKit❗❌:                       : -pseudo3DOffset;
    // RDKit❗❌:         if (refLength) {
    // RDKit❗❌:           tmpPt.z *= refLength;
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:       } else {
    // RDKit❗❌:         tmpPt.z = 0;
    // RDKit❗❌:       }
    // RDKit❗❌:       // check for overly short bonds. Note that we're doing this check *after*
    // RDKit❗❌:       // adjusting the z coordinate.
    // RDKit❗❌:       //    We want to allow atoms to overlap in x-y space if they are connected
    // RDKit❗❌:       //    via a wedged bond.
    // RDKit❗❌:       if ((centerLoc - tmpPt).lengthSq() < zeroTol) {
    // RDKit❗❌:         BOOST_LOG(rdWarningLog)
    // RDKit❗❌:             << "Warning: ambiguous stereochemistry - zero-length (or near zero-length) bond - at atom "
    // RDKit❗❌:             << atom->getIdx() << " ignored." << std::endl;
    // RDKit❗❌:         return std::nullopt;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     ++nbrIdx;
    // RDKit❗❌:     if (nbrBond->getBondType() != Bond::SINGLE) {
    // RDKit❗❌:       allSingle = false;
    // RDKit❗❌:     }
    // RDKit❗❌:     bondVects.push_back(centerLoc.directionVector(tmpPt));
    // RDKit❗❌:     neighborBondIndices.push_back(nbrBond->getIdx());
    // RDKit❗❌:   }
    // RDKit❗❌:   CHECK_INVARIANT(refIdx < mol.getNumBonds(),
    // RDKit❗❌:                   "could not find reference bond in neighbors");
    // RDKit❗❌:
    // RDKit❗❌:   auto nNbrs = bondVects.size();
    // RDKit❗❌:
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Return now if there aren't at least 3 bonds to the atom.
    // RDKit❗❌:   //  (we can implicitly add a single H to 3 coordinate atoms).
    // RDKit❗❌:   //
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   if (nNbrs < 3 || nNbrs > 4) {
    // RDKit❗❌:     return std::nullopt;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   //  Check for neighbor atoms which overlap
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   for (auto i = 0u; i < nNbrs; ++i) {
    // RDKit❗❌:     for (auto j = 0u; j < i; ++j) {
    // RDKit❗❌:       if ((bondVects[i] - bondVects[j]).lengthSq() < zeroTol) {
    // RDKit❗❌:         BOOST_LOG(rdWarningLog)
    // RDKit❗❌:             << "Warning: ambiguous stereochemistry - overlapping neighbors  - at atom "
    // RDKit❗❌:             << atom->getIdx() << " ignored" << std::endl;
    // RDKit❗❌:         return std::nullopt;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Continue if there are all single bonds or if we're considering
    // RDKit❗❌:   //  4-coordinate P or S
    // RDKit❗❌:   //
    // RDKit❗❌:   //----------------------------------------------------------
    // RDKit❗❌:   if (allSingle || atom->getAtomicNum() == 15 || atom->getAtomicNum() == 16) {
    // RDKit❗❌:     double vol;
    // RDKit❗❌:     unsigned int order[4] = {0, 1, 2, 3};
    // RDKit❗❌:     double prefactor = 1;
    // RDKit❗❌:     if (refIdx != 0) {
    // RDKit❗❌:       // bring the wedged bond to the front so that we always consider it
    // RDKit❗❌:       std::swap(order[0], order[refIdx]);
    // RDKit❗❌:       prefactor *= -1;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // check for the case that bonds 1 and 2 are co-linear but 1 and 0 are
    // RDKit❗❌:     // not:
    // RDKit❗❌:     if (nNbrs > 3 &&
    // RDKit❗❌:         bondVects[order[1]].crossProduct(bondVects[order[2]]).lengthSq() <
    // RDKit❗❌:             10 * zeroTol &&
    // RDKit❗❌:         bondVects[order[1]].crossProduct(bondVects[order[0]]).lengthSq() >
    // RDKit❗❌:             10 * zeroTol) {
    // RDKit❗❌:       bondVects[order[1]].z = bondVects[order[0]].z * -1;
    // RDKit❗❌:       // that bondVect is no longer normalized, but this hopefully won't break
    // RDKit❗❌:       // anything
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     //----------------------------------------------------------
    // RDKit❗❌:     //
    // RDKit❗❌:     // order the bonds so that the rotation order is:
    // RDKit❗❌:     //   0 - 1 - 2        for three coordinate
    // RDKit❗❌:     // or
    // RDKit❗❌:     //   0 - 1 - 2 - 3    for four coordinate
    // RDKit❗❌:     //
    // RDKit❗❌:     // this makes the rest of the code a lot simpler
    // RDKit❗❌:     //
    // RDKit❗❌:     //----------------------------------------------------------
    // RDKit❗❌:
    // RDKit❗❌:     // checks to see if the vectors 1 and 2 need to have their order
    // RDKit❗❌:     //    relative to vector 0 swapped.
    // RDKit❗❌:     // we don't actually pass the vectors in, but use their cross products
    // RDKit❗❌:     // and dot products to vector 0 to figure out if they need to be swapped
    // RDKit❗❌: #if defined(__clang__)
    // RDKit❗❌: // Clang apparently doesn't need to capture the constexpr zeroTol, and complains
    // RDKit❗❌: // about it being specified, but MSVC does need it, and removing it will break
    // RDKit❗❌: // the build
    // RDKit❗❌: #pragma GCC diagnostic push
    // RDKit❗❌: #pragma GCC diagnostic ignored "-Wunused-lambda-capture"
    // RDKit❗❌: #endif
    // RDKit❗❌:     auto needsSwap = [&zeroTol](const RDGeom::Point3D &cp01,
    // RDKit❗❌:                                 const RDGeom::Point3D &cp02, double dp01,
    // RDKit❗❌:                                 double dp02) -> bool {
    // RDKit❗❌:       if (fabs(dp01) - 1 > -zeroTol) {
    // RDKit❗❌:         if (cp02.z < 0) {
    // RDKit❗❌:           return true;
    // RDKit❗❌:         }
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (fabs(dp02) - 1 > -zeroTol) {
    // RDKit❗❌:         if (cp01.z < 0) {
    // RDKit❗❌:           return true;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       if ((cp01.z * cp02.z) < -zeroTol) {
    // RDKit❗❌:         if (cp01.z < cp02.z) {
    // RDKit❗❌:           return true;
    // RDKit❗❌:         }
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (dp01 * dp02 < -zeroTol) {
    // RDKit❗❌:         if (dp01 < dp02) {
    // RDKit❗❌:           return true;
    // RDKit❗❌:         }
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       return fabs(dp01) > fabs(dp02);
    // RDKit❗❌:     };
    // RDKit❗❌: #if defined(__clang__)
    // RDKit❗❌: #pragma GCC diagnostic pop
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:     if (nNbrs == 3) {
    // RDKit❗❌:       // this case is simple, we either need to swap vectors 1 and 2 or we
    // RDKit❗❌:       // don't:
    // RDKit❗❌:       auto cp01 = bondVects[order[0]].crossProduct(bondVects[order[1]]);
    // RDKit❗❌:       auto cp02 = bondVects[order[0]].crossProduct(bondVects[order[2]]);
    // RDKit❗❌:       auto dp01 = bondVects[order[0]].dotProduct(bondVects[order[1]]);
    // RDKit❗❌:       auto dp02 = bondVects[order[0]].dotProduct(bondVects[order[2]]);
    // RDKit❗❌:       if (needsSwap(cp01, cp02, dp01, dp02)) {
    // RDKit❗❌:         std::swap(order[1], order[2]);
    // RDKit❗❌:         prefactor *= -1;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else if (nNbrs > 3) {
    // RDKit❗❌:       // here there are more permutations. Rather than hand-coding all of them
    // RDKit❗❌:       // we'll just sort bonds 1, 2, and 3 based on their cross- and dot-
    // RDKit❗❌:       // products to bond 0
    // RDKit❗❌:       std::vector<std::tuple<double, double, unsigned>> orderedBonds(3);
    // RDKit❗❌:       for (auto i = 1u; i < 4; ++i) {
    // RDKit❗❌:         auto cp0i = bondVects[order[0]].crossProduct(bondVects[order[i]]);
    // RDKit❗❌:         auto sgn = cp0i.z < -zeroTol ? -1 : 1;
    // RDKit❗❌:         auto dp0i = bondVects[order[0]].dotProduct(bondVects[order[i]]);
    // RDKit❗❌:         orderedBonds[i - 1] = std::make_tuple(sgn, sgn * dp0i, order[i]);
    // RDKit❗❌:       }
    // RDKit❗❌:       std::sort(orderedBonds.rbegin(), orderedBonds.rend());
    // RDKit❗❌:
    // RDKit❗❌:       // update the order array and figure out whether or not we've done a
    // RDKit❗❌:       // cyclic permutation
    // RDKit❗❌:       auto nChanged = 0;
    // RDKit❗❌:       for (auto i = 1u; i < 4; ++i) {
    // RDKit❗❌:         auto ni = std::get<2>(orderedBonds[i - 1]);
    // RDKit❗❌:         if (order[i] != ni) {
    // RDKit❗❌:           order[i] = ni;
    // RDKit❗❌:           ++nChanged;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (nChanged == 2) {
    // RDKit❗❌:         // this is always an acyclic permutation
    // RDKit❗❌:         prefactor *= -1;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // std::cerr<<"ORDER "<<neighborBondIndices[order[0]]<<"
    // RDKit❗❌:     // "<<neighborBondIndices[order[1]]<<" "<<neighborBondIndices[order[2]]<<"
    // RDKit❗❌:     // "<<neighborBondIndices[order[3]]<<std::endl;
    // RDKit❗❌:
    // RDKit❗❌:     // check for opposing bonds with opposite wedging
    // RDKit❗❌:     for (auto i = 0u; i < nNbrs; ++i) {
    // RDKit❗❌:       for (auto j = i + 1; j < nNbrs; ++j) {
    // RDKit❗❌:         if (bondVects[order[i]].z * bondVects[order[j]].z < -zeroTol) {
    // RDKit❗❌:           auto cp =
    // RDKit❗❌:               bondVects[order[i]].crossProduct(bondVects[order[j]]).lengthSq();
    // RDKit❗❌:           if (cp < 0.01) {
    // RDKit❗❌:             // exception to our rejection of these structures: in some horrible
    // RDKit❗❌:             // pseudo-3D drawings of things like sugars the ring substituents
    // RDKit❗❌:             // are drawn 180 degrees apart and with opposite wedging. Let that
    // RDKit❗❌:             // one pass.
    // RDKit❗❌:             if (nNbrs == 4 &&
    // RDKit❗❌:                 fabs(bondVects[order[i]].dotProduct(bondVects[order[j]]) + 1) <
    // RDKit❗❌:                     zeroTol) {
    // RDKit❗❌:               // this is allowed for neighboring bonds
    // RDKit❗❌:               if (j - i == 1 || (i == 0 && j == 3)) {
    // RDKit❗❌:                 // std::cerr << " skip it " << std::endl;
    // RDKit❗❌:                 bondVects[order[j]].z = 0.0;
    // RDKit❗❌:                 continue;
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             BOOST_LOG(rdWarningLog)
    // RDKit❗❌:                 << "Warning: ambiguous stereochemistry - opposing bonds have opposite wedging - at atom "
    // RDKit❗❌:                 << atom->getIdx() << " ignored." << std::endl;
    // RDKit❗❌:             return std::nullopt;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // three-coordinate special cases where chirality cannot be determined
    // RDKit❗❌:     //
    // RDKit❗❌:     //  Case 1:
    // RDKit❗❌:     //  this one is never allowed with different directions for the bonds to 1
    // RDKit❗❌:     //  and 2
    // RDKit❗❌:     //     0   2
    // RDKit❗❌:     //      \ /
    // RDKit❗❌:     //       C
    // RDKit❗❌:     //       *
    // RDKit❗❌:     //       1
    // RDKit❗❌:     //   This is ST-1.2.10 in the IUPAC guidelines
    // RDKit❗❌:     //
    // RDKit❗❌:     //  Case 2: all bonds are wedged in the same direction
    // RDKit❗❌:     if (nNbrs == 3) {
    // RDKit❗❌:       bool conflict = false;
    // RDKit❗❌:       if (bondVects[order[1]].z * bondVects[order[0]].z < -coordZeroTol &&
    // RDKit❗❌:           fabs(bondVects[order[2]].z) < coordZeroTol) {
    // RDKit❗❌:         conflict = bondVects[order[2]].crossProduct(bondVects[order[0]]).z *
    // RDKit❗❌:                        bondVects[order[2]].crossProduct(bondVects[order[1]]).z <
    // RDKit❗❌:                    -1e-4;
    // RDKit❗❌:       } else if (bondVects[order[2]].z * bondVects[order[0]].z <
    // RDKit❗❌:                      -coordZeroTol &&
    // RDKit❗❌:                  fabs(bondVects[order[1]].z) < coordZeroTol) {
    // RDKit❗❌:         conflict = bondVects[order[1]].crossProduct(bondVects[order[0]]).z *
    // RDKit❗❌:                        bondVects[order[1]].crossProduct(bondVects[order[2]]).z <
    // RDKit❗❌:                    -coordZeroTol;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (conflict) {
    // RDKit❗❌:         BOOST_LOG(rdWarningLog)
    // RDKit❗❌:             << "Warning: conflicting stereochemistry - bond wedging contradiction - at atom "
    // RDKit❗❌:             << atom->getIdx() << " ignored" << std::endl;
    // RDKit❗❌:         return std::nullopt;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // for the purposes of the cross products we ignore any pseudo-3D
    // RDKit❗❌:     // coordinates
    // RDKit❗❌:     auto bv1 = bondVects[order[1]];
    // RDKit❗❌:     bv1.z = 0;
    // RDKit❗❌:     auto bv2 = bondVects[order[2]];
    // RDKit❗❌:     bv2.z = 0;
    // RDKit❗❌:     auto crossp1 = bv1.crossProduct(bv2);
    // RDKit❗❌:     // catch linear arrangements
    // RDKit❗❌:     if (nNbrs == 3) {
    // RDKit❗❌:       if (crossp1.lengthSq() < tShapeTol) {
    // RDKit❗❌:         // in a linear relationship with three neighbors we assume that the
    // RDKit❗❌:         // two perpendicular bonds are wedged in the other direction from the
    // RDKit❗❌:         // one that was provided.
    // RDKit❗❌:         // that's this situation:
    // RDKit❗❌:         //
    // RDKit❗❌:         //              0
    // RDKit❗❌:         //              |   <- wedged up
    // RDKit❗❌:         //           1--C--2
    // RDKit❗❌:         //
    // RDKit❗❌:         //  here we assume that bonds C-1 and C-2 are wedged down
    // RDKit❗❌:         //
    // RDKit❗❌:         // ST-1.2.12 of the IUPAC guidelines says that this form is wrong since
    // RDKit❗❌:         // it's for a "T-shaped" configuration instead of a tetrahedron, but it
    // RDKit❗❌:         // shows up fairly frequently, particularly with fused ring systems
    // RDKit❗❌:         bv1.z = -bondVects[order[0]].z;
    // RDKit❗❌:         bv2.z = -bondVects[order[0]].z;
    // RDKit❗❌:         crossp1 = bv1.crossProduct(bv2);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else if (crossp1.lengthSq() < 10 * zeroTol) {
    // RDKit❗❌:       // if the other bond is flat:
    // RDKit❗❌:       if (fabs(bondVects[order[3]].z) < coordZeroTol) {
    // RDKit❗❌:         // By construction this is a neighboring bond, so make it the opposite
    // RDKit❗❌:         // wedging from us.
    // RDKit❗❌:         bondVects[order[3]].z = -1 * bondVects[order[0]].z;
    // RDKit❗❌:         // that bondVect is no longer normalized, but this hopefully won't break
    // RDKit❗❌:         // anything
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     vol = crossp1.dotProduct(bondVects[order[0]]);
    // RDKit❗❌:     if (nNbrs == 4) {
    // RDKit❗❌:       const auto dotp1 = bondVects[order[1]].dotProduct(bondVects[order[2]]);
    // RDKit❗❌:       // for the purposes of the cross products we ignore any pseudo-3D
    // RDKit❗❌:       // coordinates
    // RDKit❗❌:       auto bv3 = bondVects[order[3]];
    // RDKit❗❌:       bv3.z = 0;
    // RDKit❗❌:       const auto crossp2 = bv1.crossProduct(bv3);
    // RDKit❗❌:       const auto dotp2 = bondVects[order[1]].dotProduct(bondVects[order[3]]);
    // RDKit❗❌:       auto vol2 = crossp2.dotProduct(bondVects[order[0]]);
    // RDKit❗❌:
    // RDKit❗❌:       // detect the case where there's no chiral volume for the default
    // RDKit❗❌:       // evaluation
    // RDKit❗❌:       if (fabs(vol) < zeroTol) {
    // RDKit❗❌:         // and check the other evaluation:
    // RDKit❗❌:         if (fabs(vol2) < zeroTol) {
    // RDKit❗❌:           BOOST_LOG(rdWarningLog)
    // RDKit❗❌:               << "Warning: ambiguous stereochemistry - no chiral volume - at atom "
    // RDKit❗❌:               << atom->getIdx() << " ignored" << std::endl;
    // RDKit❗❌:           return std::nullopt;
    // RDKit❗❌:         }
    // RDKit❗❌:         vol = vol2;
    // RDKit❗❌:         prefactor *= -1;
    // RDKit❗❌:       } else if (vol * vol2 > 0 && fabs(vol2) > volumeTolerance &&
    // RDKit❗❌:                  dotp1 < dotp2) {
    // RDKit❗❌:         // both volumes give the same answer, but in the second case the cross
    // RDKit❗❌:         // product is between two bonds with a better dot product
    // RDKit❗❌:         vol = vol2;
    // RDKit❗❌:         prefactor *= -1;
    // RDKit❗❌:       } else if (fabs(vol) < volumeTolerance && fabs(vol2) > volumeTolerance) {
    // RDKit❗❌:         // if the first volume is too small, but the second isn't, take the
    // RDKit❗❌:         // second
    // RDKit❗❌:         if (vol * vol2 < 0) {
    // RDKit❗❌:           prefactor *= -1;
    // RDKit❗❌:         }
    // RDKit❗❌:         vol = vol2;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     vol *= prefactor;
    // RDKit❗❌:     // std::cerr << " final " << vol << std::endl;
    // RDKit❗❌:
    // RDKit❗❌:     // at this point we can assign our atomic stereo based on the sign of the
    // RDKit❗❌:     // chiral volume
    // RDKit❗❌:     if (vol > volumeTolerance) {
    // RDKit❗❌:       res = Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    // RDKit❗❌:     } else if (vol < -volumeTolerance) {
    // RDKit❗❌:       res = Atom::ChiralType::CHI_TETRAHEDRAL_CW;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog)
    // RDKit❗❌:           << "Warning: ambiguous stereochemistry - zero final chiral volume - at atom "
    // RDKit❗❌:           << atom->getIdx() << " ignored" << std::endl;
    // RDKit❗❌:       return std::nullopt;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION Chirality::atomChiralTypeFromBondDirPseudo3D complete source
    // The native body is preserved verbatim. The extra neighbor copy/sort
    // restores physical bond insertion order from detached adjacency; warning
    // emission and nonfinite tuple ordering remain unproven, not source parity.
    const COORD_ZERO_TOL: f64 = 1e-4;
    const ZERO_TOL: f64 = 1e-3;
    const T_SHAPE_TOL: f64 = 0.00031;
    const PSEUDO_3D_OFFSET: f64 = 0.1;
    const VOLUME_TOLERANCE: f64 = 0.00174;

    let bond = &topology.bonds[bond_id.index()];
    let direction = bond.direction();
    if !matches!(
        direction,
        BondDirection::BeginWedge | BondDirection::BeginDash
    ) {
        return Ok(None);
    }
    let center = bond.begin();
    let center_neighbors = topology.adjacency.neighbors_of(center.index());
    if center_neighbors.len() > 4 {
        return Ok(Some(ChiralTag::Unspecified));
    }
    let coords = conformer.coordinates();
    let mut center_point = coords[center.index()];
    center_point[2] = 0.0;
    let mut reference_point = coords[bond.end().index()];
    let reference_length = Vector3::between(reference_point, center_point).length();
    reference_point[2] = if direction == BondDirection::BeginWedge {
        PSEUDO_3D_OFFSET
    } else {
        -PSEUDO_3D_OFFSET
    };
    if reference_length != 0.0 {
        reference_point[2] *= reference_length;
    }

    let mut neighbors = center_neighbors.to_vec();
    neighbors.sort_by_key(|neighbor| neighbor.bond.index());
    let mut reference_index = topology.bonds.len() + 1;
    let mut vectors = Vec::with_capacity(neighbors.len());
    let mut all_single = true;
    for (neighbor_index, neighbor) in neighbors.iter().enumerate() {
        let neighbor_bond = &topology.bonds[neighbor.bond.index()];
        let mut point = coords[neighbor.atom_index];
        if neighbor.bond == bond_id {
            reference_index = neighbor_index;
            point = reference_point;
        } else {
            if neighbor_bond.begin() == center
                && matches!(
                    neighbor_bond.direction(),
                    BondDirection::BeginWedge | BondDirection::BeginDash
                )
            {
                point[2] = if neighbor_bond.direction() == BondDirection::BeginWedge {
                    PSEUDO_3D_OFFSET
                } else {
                    -PSEUDO_3D_OFFSET
                };
                if reference_length != 0.0 {
                    point[2] *= reference_length;
                }
            } else {
                point[2] = 0.0;
            }
            if Vector3::between(point, center_point).length_squared() < ZERO_TOL {
                return Ok(None);
            }
        }
        if neighbor_bond.order() != BondOrder::Single {
            all_single = false;
        }
        vectors.push(
            Vector3::between(center_point, point)
                .normalized(center, AtomId::new(neighbor.atom_index))?,
        );
    }
    if !(3..=4).contains(&vectors.len()) || reference_index >= vectors.len() {
        return Ok(None);
    }
    for index in 0..vectors.len() {
        for previous in 0..index {
            if vectors[index]
                .difference(vectors[previous])
                .length_squared()
                < ZERO_TOL
            {
                return Ok(None);
            }
        }
    }
    if !all_single && !matches!(topology.atoms[center.index()].atomic_number(), 15 | 16) {
        return Ok(Some(ChiralTag::Unspecified));
    }

    let mut order = [0, 1, 2, 3];
    let mut prefactor = 1.0;
    if reference_index != 0 {
        order.swap(0, reference_index);
        prefactor *= -1.0;
    }
    if vectors.len() > 3
        && vectors[order[1]].cross(vectors[order[2]]).length_squared() < 10.0 * ZERO_TOL
        && vectors[order[1]].cross(vectors[order[0]]).length_squared() > 10.0 * ZERO_TOL
    {
        vectors[order[1]].z = -vectors[order[0]].z;
    }

    let needs_swap = |cp01: Vector3, cp02: Vector3, dp01: f64, dp02: f64| {
        if dp01.abs() - 1.0 > -ZERO_TOL {
            return cp02.z < 0.0;
        }
        if dp02.abs() - 1.0 > -ZERO_TOL && cp01.z < 0.0 {
            return true;
        }
        if cp01.z * cp02.z < -ZERO_TOL {
            return cp01.z < cp02.z;
        }
        if dp01 * dp02 < -ZERO_TOL {
            return dp01 < dp02;
        }
        dp01.abs() > dp02.abs()
    };
    if vectors.len() == 3 {
        let cp01 = vectors[order[0]].cross(vectors[order[1]]);
        let cp02 = vectors[order[0]].cross(vectors[order[2]]);
        let dp01 = vectors[order[0]].dot(vectors[order[1]]);
        let dp02 = vectors[order[0]].dot(vectors[order[2]]);
        if needs_swap(cp01, cp02, dp01, dp02) {
            order.swap(1, 2);
            prefactor *= -1.0;
        }
    } else {
        let mut ordered = (1..4)
            .map(|index| {
                let cross = vectors[order[0]].cross(vectors[order[index]]);
                let sign = if cross.z < -ZERO_TOL { -1.0 } else { 1.0 };
                (
                    sign,
                    sign * vectors[order[0]].dot(vectors[order[index]]),
                    order[index],
                )
            })
            .collect::<Vec<_>>();
        ordered.sort_by(|left, right| {
            right
                .0
                .total_cmp(&left.0)
                .then_with(|| right.1.total_cmp(&left.1))
                .then_with(|| right.2.cmp(&left.2))
        });
        let mut changed = 0;
        for index in 1..4 {
            if order[index] != ordered[index - 1].2 {
                order[index] = ordered[index - 1].2;
                changed += 1;
            }
        }
        if changed == 2 {
            prefactor *= -1.0;
        }
    }

    for index in 0..vectors.len() {
        for next in index + 1..vectors.len() {
            if vectors[order[index]].z * vectors[order[next]].z < -ZERO_TOL
                && vectors[order[index]]
                    .cross(vectors[order[next]])
                    .length_squared()
                    < 0.01
            {
                if vectors.len() == 4
                    && (vectors[order[index]].dot(vectors[order[next]]) + 1.0).abs() < ZERO_TOL
                    && (next - index == 1 || (index == 0 && next == 3))
                {
                    vectors[order[next]].z = 0.0;
                    continue;
                }
                return Ok(None);
            }
        }
    }

    if vectors.len() == 3 {
        let mut conflict = false;
        if vectors[order[1]].z * vectors[order[0]].z < -COORD_ZERO_TOL
            && vectors[order[2]].z.abs() < COORD_ZERO_TOL
        {
            conflict = vectors[order[2]].cross(vectors[order[0]]).z
                * vectors[order[2]].cross(vectors[order[1]]).z
                < -1e-4;
        } else if vectors[order[2]].z * vectors[order[0]].z < -COORD_ZERO_TOL
            && vectors[order[1]].z.abs() < COORD_ZERO_TOL
        {
            conflict = vectors[order[1]].cross(vectors[order[0]]).z
                * vectors[order[1]].cross(vectors[order[2]]).z
                < -COORD_ZERO_TOL;
        }
        if conflict {
            return Ok(None);
        }
    }

    let mut vector1 = vectors[order[1]];
    vector1.z = 0.0;
    let mut vector2 = vectors[order[2]];
    vector2.z = 0.0;
    let mut cross1 = vector1.cross(vector2);
    if vectors.len() == 3 && cross1.length_squared() < T_SHAPE_TOL {
        vector1.z = -vectors[order[0]].z;
        vector2.z = -vectors[order[0]].z;
        cross1 = vector1.cross(vector2);
    } else if vectors.len() == 4
        && cross1.length_squared() < 10.0 * ZERO_TOL
        && vectors[order[3]].z.abs() < COORD_ZERO_TOL
    {
        vectors[order[3]].z = -vectors[order[0]].z;
    }
    let mut volume = cross1.dot(vectors[order[0]]);
    if vectors.len() == 4 {
        let dot1 = vectors[order[1]].dot(vectors[order[2]]);
        let mut vector3 = vectors[order[3]];
        vector3.z = 0.0;
        let cross2 = vector1.cross(vector3);
        let dot2 = vectors[order[1]].dot(vectors[order[3]]);
        let volume2 = cross2.dot(vectors[order[0]]);
        if volume.abs() < ZERO_TOL {
            if volume2.abs() < ZERO_TOL {
                return Ok(None);
            }
            volume = volume2;
            prefactor *= -1.0;
        } else if volume * volume2 > 0.0 && volume2.abs() > VOLUME_TOLERANCE && dot1 < dot2 {
            volume = volume2;
            prefactor *= -1.0;
        } else if volume.abs() < VOLUME_TOLERANCE && volume2.abs() > VOLUME_TOLERANCE {
            if volume * volume2 < 0.0 {
                prefactor *= -1.0;
            }
            volume = volume2;
        }
    }
    volume *= prefactor;
    if volume > VOLUME_TOLERANCE {
        Ok(Some(ChiralTag::TetrahedralCcw))
    } else if volume < -VOLUME_TOLERANCE {
        Ok(Some(ChiralTag::TetrahedralCw))
    } else {
        Ok(None)
    }
}

/// Assign tetrahedral tags from molfile-style wedge/dash bonds and a 2D
/// conformer, matching RDKit's CXSMILES post-parser path.
pub fn assign_chiral_types_from_bond_dirs(
    topology: &mut TopologyBlock,
    conformer: &Conformer3D,
    replace_existing_tags: bool,
) -> Result<(), BondDirectionStereoError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::assignChiralTypesFromBondDirs complete source
    // RDKit✔️❌: void assignChiralTypesFromBondDirs(ROMol &mol, const int confId,
    // RDKit✔️❌:                                    const bool replaceExistingTags) {
    // RDKit✔️❌:   if (!mol.getNumConformers()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto conf = mol.getConformer(confId);
    // RDKit✔️❌:   boost::dynamic_bitset<> atomsSet(mol.getNumAtoms(), 0);
    // RDKit✔️❌:   for (auto &bond : mol.bonds()) {
    // RDKit✔️❌:     const Bond::BondDir dir = bond->getBondDir();
    // RDKit✔️❌:     Atom *atom = bond->getBeginAtom();
    // RDKit✔️❌:     if (dir == Bond::UNKNOWN) {
    // RDKit✔️❌:       if (atomsSet[atom->getIdx()] || replaceExistingTags) {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:         atomsSet.set(atom->getIdx());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // the bond is marked as chiral:
    // RDKit✔️❌:       if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    // RDKit✔️❌:         if (atomsSet[atom->getIdx()] ||
    // RDKit✔️❌:             (!replaceExistingTags &&
    // RDKit✔️❌:              atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
    // RDKit✔️❌:           continue;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (atom->needsUpdatePropertyCache()) {
    // RDKit✔️❌:           atom->updatePropertyCache(false);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         Atom::ChiralType code =
    // RDKit✔️❌:             Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
    // RDKit✔️❌:                 .value_or(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️❌:         if (code != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit✔️❌:           atomsSet.set(atom->getIdx());
    // RDKit✔️❌:           //   std::cerr << "atom " << atom->getIdx() << " code " << code
    // RDKit✔️❌:           //             << " from bond " << bond->getIdx() << std::endl;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         atom->setChiralTag(code);
    // RDKit✔️❌:
    // RDKit✔️❌:         // within the RD representation, if a three-coordinate atom
    // RDKit✔️❌:         // is chiral and has an implicit H, that H needs to be made explicit:
    // RDKit✔️❌:         if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
    // RDKit✔️❌:             atom->getNumImplicitHs() == 1) {
    // RDKit✔️❌:           atom->setNumExplicitHs(1);
    // RDKit✔️❌:           // recalculated number of implicit Hs:
    // RDKit✔️❌:           atom->updatePropertyCache();
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::assignChiralTypesFromBondDirs complete source
    // The caller supplies the selected, present native conformer. Borrowing
    // preserves its rows without the native auto conf copy. Structural input
    // validation adds O(V+E) work; no chemical cache is eagerly recomputed.
    conformer
        .validate_for_atom_count(topology.atoms.len())
        .map_err(|error| BondDirectionStereoError::InvalidState(error.to_string()))?;
    topology
        .validate()
        .map_err(|error| BondDirectionStereoError::InvalidState(error.to_string()))?;
    let mut assigned = vec![false; topology.atoms.len()];
    for bond_index in 0..topology.bonds.len() {
        let bond_id = BondId::new(bond_index);
        let direction = topology.bonds[bond_index].direction();
        let atom = topology.bonds[bond_index].begin();
        if direction == BondDirection::Unknown {
            if assigned[atom.index()] || replace_existing_tags {
                topology.atoms[atom.index()].set_chiral_tag(ChiralTag::Unspecified);
                assigned[atom.index()] = true;
            }
        } else if matches!(
            direction,
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) {
            if assigned[atom.index()]
                || (!replace_existing_tags
                    && topology.atoms[atom.index()].chiral_tag() != ChiralTag::Unspecified)
            {
                continue;
            }
            if source_atom_needs_cache_update(&topology.atoms[atom.index()]) {
                update_source_atom_cache(topology, atom, false)?;
            }
            let tag = pseudo_3d_chiral_tag(topology, bond_id, conformer)?
                .unwrap_or(ChiralTag::Unspecified);
            if tag != ChiralTag::Unspecified {
                assigned[atom.index()] = true;
            }
            topology.atoms[atom.index()].set_chiral_tag(tag);
            if topology.adjacency.neighbors_of(atom.index()).len() == 3
                && topology.atoms[atom.index()].explicit_hydrogens() == 0
                && source_atom_implicit_hydrogens(&topology.atoms[atom.index()])
                    .map_err(|message| BondDirectionStereoError::InvalidState(message.into()))?
                    == 1
            {
                topology.atoms[atom.index()].set_explicit_hydrogens(1);
                update_source_atom_cache(topology, atom, true)?;
            }
        }
    }
    Ok(())
}
