use cosmolkit_core::{ValenceModel, assign_valence_with_options_for_topology};
use cosmolkit_model::{BondId, Conformer3D, TopologyBlock};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag};

use crate::StereoError;

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

    fn normalized(self) -> Self {
        let length = self.length();
        if length == 0.0 {
            self
        } else {
            Self {
                x: self.x / length,
                y: self.y / length,
                z: self.z / length,
            }
        }
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
) -> Option<ChiralTag> {
    // BEGIN RDKIT CPP FUNCTION Chirality::atomChiralTypeFromBondDirPseudo3D
    // RDKit✔️✔️: auto bondDir = bond->getBondDir();
    // RDKit✔️✔️: PRECONDITION(bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH,
    // RDKit✔️✔️:              "bad bond direction");
    // RDKit✔️✔️: constexpr double coordZeroTol = 1e-4;
    // RDKit✔️✔️: constexpr double zeroTol = 1e-3;
    // RDKit✔️✔️: constexpr double tShapeTol = 0.00031;
    // RDKit✔️✔️: constexpr double pseudo3DOffset = 0.1;
    // RDKit✔️✔️: constexpr double volumeTolerance = 0.00174;
    // RDKit✔️✔️: const auto atom = bond->getBeginAtom();
    // RDKit✔️✔️: if (atom->getDegree() > 4) { return Atom::CHI_UNSPECIFIED; }
    // RDKit✔️✔️: const auto bondAtom = bond->getEndAtom();
    // RDKit✔️✔️: auto centerLoc = conf->getAtomPos(atom->getIdx());
    // RDKit✔️✔️: centerLoc.z = 0.0;
    // RDKit✔️✔️: auto refPt = conf->getAtomPos(bondAtom->getIdx());
    // RDKit✔️✔️: auto refLength = (centerLoc - refPt).length();
    // RDKit✔️✔️: refPt.z = bondDir == Bond::BondDir::BEGINWEDGE ? pseudo3DOffset
    // RDKit✔️✔️:                                                     : -pseudo3DOffset;
    // RDKit✔️✔️: if (refLength) { refPt.z *= refLength; }
    // RDKit✔️✔️: for (const auto nbrBond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:   const auto oAtom = nbrBond->getOtherAtom(atom);
    // RDKit✔️✔️:   auto tmpPt = conf->getAtomPos(oAtom->getIdx());
    // RDKit✔️✔️:   if (nbrBond == bond) { refIdx = nbrIdx; tmpPt = refPt; }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (nbrBond->getBeginAtomIdx() == atom->getIdx() &&
    // RDKit✔️✔️:         (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit✔️✔️:          nbrBond->getBondDir() == Bond::BondDir::BEGINDASH)) {
    // RDKit✔️✔️:       tmpPt.z = nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE
    // RDKit✔️✔️:                     ? pseudo3DOffset
    // RDKit✔️✔️:                     : -pseudo3DOffset;
    // RDKit✔️✔️:       if (refLength) {
    // RDKit✔️✔️:         tmpPt.z *= refLength;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       tmpPt.z = 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   if ((centerLoc - tmpPt).lengthSq() < zeroTol) { return std::nullopt; }
    // RDKit✔️✔️:   if (nbrBond->getBondType() != Bond::SINGLE) { allSingle = false; }
    // RDKit✔️✔️:   bondVects.push_back(centerLoc.directionVector(tmpPt));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (nNbrs < 3 || nNbrs > 4) { return std::nullopt; }
    // RDKit✔️✔️: for (auto i = 0u; i < nNbrs; ++i) {
    // RDKit✔️✔️:   for (auto j = 0u; j < i; ++j) {
    // RDKit✔️✔️:     if ((bondVects[i] - bondVects[j]).lengthSq() < zeroTol) return std::nullopt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (allSingle || atom->getAtomicNum() == 15 || atom->getAtomicNum() == 16) {
    // RDKit✔️✔️:   unsigned int order[4] = {0, 1, 2, 3};
    // RDKit✔️✔️:   double prefactor = 1;
    // RDKit✔️✔️:   if (refIdx != 0) { std::swap(order[0], order[refIdx]); prefactor *= -1; }
    // RDKit✔️✔️:   if (nNbrs > 3 && cross(vectors 1,2).lengthSq() < 10*zeroTol &&
    // RDKit✔️✔️:       cross(vectors 1,0).lengthSq() > 10*zeroTol) {
    // RDKit✔️✔️:     bondVects[order[1]].z = bondVects[order[0]].z * -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // order the bonds so that the rotation order is:
    // RDKit✔️✔️:   //   0 - 1 - 2        for three coordinate
    // RDKit✔️✔️:   // or
    // RDKit✔️✔️:   //   0 - 1 - 2 - 3    for four coordinate
    // RDKit✔️✔️:   auto bv1 = bondVects[order[1]]; bv1.z = 0;
    // RDKit✔️✔️:   auto bv2 = bondVects[order[2]]; bv2.z = 0;
    // RDKit✔️✔️:   auto crossp1 = bv1.crossProduct(bv2);
    // RDKit✔️✔️:   // for the purposes of the cross products we ignore any pseudo-3D
    // RDKit✔️✔️:   // coordinates
    // RDKit✔️✔️:   vol = crossp1.dotProduct(bondVects[order[0]]);
    // RDKit✔️✔️:   if (nNbrs == 4) {
    // RDKit✔️✔️:     const auto dotp1 = bondVects[order[1]].dotProduct(bondVects[order[2]]);
    // RDKit✔️✔️:   vol *= prefactor;
    // RDKit✔️✔️:   if (vol > volumeTolerance) res = Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️:   else if (vol < -volumeTolerance) res = Atom::CHI_TETRAHEDRAL_CW;
    // RDKit✔️✔️:   else return std::nullopt;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return res;
    // END RDKIT CPP FUNCTION Chirality::atomChiralTypeFromBondDirPseudo3D
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
        return None;
    }
    let center = bond.begin();
    let center_neighbors = topology.adjacency.neighbors_of(center.index());
    if center_neighbors.len() > 4 {
        return Some(ChiralTag::Unspecified);
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
                return None;
            }
        }
        if neighbor_bond.order() != BondOrder::Single {
            all_single = false;
        }
        vectors.push(Vector3::between(center_point, point).normalized());
    }
    if !(3..=4).contains(&vectors.len()) || reference_index >= vectors.len() {
        return None;
    }
    for index in 0..vectors.len() {
        for previous in 0..index {
            if vectors[index]
                .difference(vectors[previous])
                .length_squared()
                < ZERO_TOL
            {
                return None;
            }
        }
    }
    if !all_single && !matches!(topology.atoms[center.index()].atomic_number(), 15 | 16) {
        return Some(ChiralTag::Unspecified);
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
                return None;
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
            return None;
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
                return None;
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
        Some(ChiralTag::TetrahedralCcw)
    } else if volume < -VOLUME_TOLERANCE {
        Some(ChiralTag::TetrahedralCw)
    } else {
        None
    }
}

/// Assign tetrahedral tags from molfile-style wedge/dash bonds and a 2D
/// conformer, matching RDKit's CXSMILES post-parser path.
pub fn assign_chiral_types_from_bond_dirs(
    topology: &mut TopologyBlock,
    conformer: &Conformer3D,
    replace_existing_tags: bool,
) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::assignChiralTypesFromBondDirs
    // RDKit✔️✔️: if (!mol.getNumConformers()) { return; }
    // RDKit✔️✔️: boost::dynamic_bitset<> atomsSet(mol.getNumAtoms(), 0);
    // RDKit✔️✔️: for (auto &bond : mol.bonds()) {
    // RDKit✔️✔️:   const Bond::BondDir dir = bond->getBondDir();
    // RDKit✔️✔️:   Atom *atom = bond->getBeginAtom();
    // RDKit✔️✔️:   if (dir == Bond::UNKNOWN) {
    // RDKit✔️✔️:     if (atomsSet[atom->getIdx()] || replaceExistingTags) {
    // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:       atomsSet.set(atom->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    // RDKit✔️✔️:     if (atomsSet[atom->getIdx()] || (!replaceExistingTags &&
    // RDKit✔️✔️:         atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) { continue; }
    // RDKit✔️✔️:     Atom::ChiralType code =
    // RDKit✔️✔️:         Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
    // RDKit✔️✔️:             .value_or(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:     if (code != Atom::CHI_UNSPECIFIED) { atomsSet.set(atom->getIdx()); }
    // RDKit✔️✔️:     atom->setChiralTag(code);
    // RDKit✔️✔️:     if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
    // RDKit✔️✔️:         atom->getNumImplicitHs() == 1) {
    // RDKit✔️✔️:       atom->setNumExplicitHs(1);
    // RDKit✔️✔️:       atom->updatePropertyCache();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::assignChiralTypesFromBondDirs
    conformer
        .validate_for_atom_count(topology.atoms.len())
        .map_err(|error| StereoError::InvalidState(error.to_string()))?;
    let valence =
        assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, false)?;
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
            let tag = pseudo_3d_chiral_tag(topology, bond_id, conformer)
                .unwrap_or(ChiralTag::Unspecified);
            if tag != ChiralTag::Unspecified {
                assigned[atom.index()] = true;
            }
            topology.atoms[atom.index()].set_chiral_tag(tag);
            if topology.adjacency.neighbors_of(atom.index()).len() == 3
                && topology.atoms[atom.index()].explicit_hydrogens() == 0
                && valence.implicit_hydrogens[atom.index()] == 1
            {
                topology.atoms[atom.index()].set_explicit_hydrogens(1);
            }
        }
    }
    Ok(())
}
