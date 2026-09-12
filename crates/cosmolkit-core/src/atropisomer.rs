//! Detached, source-backed atropisomer perception and wedge assignment.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, Bond, BondId, Conformer2D, Conformer3D, CoordinateValidationError, StereoGroup,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Hybridization};

use crate::RingInfo;

const REALLY_SMALL_BOND_LEN: f64 = 0.000_000_1;

#[derive(Debug, Clone, Copy)]
pub enum AtropisomerConformer<'a> {
    TwoD(&'a Conformer2D),
    ThreeD(&'a Conformer3D),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtropisomerRejectionKind {
    MissingCarrier,
    UnknownCarrierDirection,
    InconsistentDirections,
    ZeroLengthAxis,
    CollinearCarrier,
    SameSideCarriers,
    CoplanarCarriers,
    NoUsableWedgeBond,
    DirectionConflict,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtropisomerDiagnostic {
    pub bond: BondId,
    pub kind: AtropisomerRejectionKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtropisomerBondUpdate {
    pub bond: BondId,
    pub stereo: BondStereo,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AtropisomerAssignment {
    pub bond_updates: Vec<AtropisomerBondUpdate>,
    pub diagnostics: Vec<AtropisomerDiagnostic>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtropisomerWedgeUpdate {
    pub bond: BondId,
    pub begin: AtomId,
    pub end: AtomId,
    pub direction: BondDirection,
    pub atropisomer_bond: BondId,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AtropisomerWedgeAssignment {
    pub bond_updates: Vec<AtropisomerWedgeUpdate>,
    pub diagnostics: Vec<AtropisomerDiagnostic>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct StereoGroupAssignment {
    pub groups: Vec<StereoGroup>,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AtropisomerError {
    #[error("invalid topology: {source}")]
    InvalidTopology { source: TopologyValidationError },
    #[error("invalid coordinates: {source}")]
    InvalidCoordinates { source: CoordinateValidationError },
    #[error("3D conformer {conformer} is marked as non-3D")]
    ConformerNotThreeDimensional { conformer: usize },
    #[error("ring information must be SSSR or better")]
    RingInfoNotSssr,
    #[error("ring atom {atom} is out of range for {atom_count} atoms")]
    RingAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("ring bond {bond} is out of range for {bond_count} bonds")]
    RingBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("stereo group atom {atom} is out of range for {atom_count} atoms")]
    StereoGroupAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("stereo group bond {bond} is out of range for {bond_count} bonds")]
    StereoGroupBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("atropisomer assignment bond {bond} is out of range for {bond_count} bonds")]
    AssignmentBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("bond {bond} has invalid atropisomer stereo {stereo:?}")]
    InvalidAtropisomerStereo { bond: BondId, stereo: BondStereo },
}

#[derive(Debug, Clone)]
struct AtropEnd {
    atom: AtomId,
    bonds: Vec<BondId>,
}

#[derive(Debug, Clone, Copy)]
struct Frame {
    y: [f64; 3],
    z: [f64; 3],
}

fn can_have_direction(bond: &Bond) -> bool {
    // Complete pinned source: Bond.h::canHaveDirection.
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
    matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
}

fn other_atom(bond: &Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom {
        bond.end()
    } else {
        bond.begin()
    }
}

fn total_degree(topology: &TopologyBlock, atom: AtomId) -> usize {
    let value = &topology.atoms[atom.index()];
    topology.adjacency.neighbors_of(atom.index()).len()
        + usize::from(value.explicit_hydrogens())
        + usize::from(value.implicit_hydrogen())
}

fn atropisomer_ends(topology: &TopologyBlock, bond: &Bond) -> Option<[AtropEnd; 2]> {
    // Complete pinned source: getAtropisomerAtomsAndBonds.
    // RDKit✔️✔️: bool getAtropisomerAtomsAndBonds(const Bond *bond,
    // RDKit✔️✔️:                                  AtropAtomAndBondVec atomsAndBondVects[2],
    // RDKit✔️✔️:                                  const ROMol &mol) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    // RDKit✔️✔️:   atomsAndBondVects[0].first = bond->getBeginAtom();
    // RDKit✔️✔️:   atomsAndBondVects[1].first = bond->getEndAtom();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // get the one or two bonds on each end
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit✔️✔️:     for (const auto nbrBond :
    // RDKit✔️✔️:          mol.atomBonds(atomsAndBondVects[bondAtomIndex].first)) {
    // RDKit✔️✔️:       if (nbrBond == bond) {
    // RDKit✔️✔️:         continue;  // a bond is NOT its own neighbor
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atomsAndBondVects[bondAtomIndex].second.push_back(nbrBond);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 0) {
    // RDKit✔️✔️:       return false;  // no neighbor bonds found
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // make sure the bond with this lowest atom is is first
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 2 &&
    // RDKit✔️✔️:         atomsAndBondVects[bondAtomIndex]
    // RDKit✔️✔️:                 .second[1]
    // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
    // RDKit✔️✔️:                 ->getIdx() <
    // RDKit✔️✔️:             atomsAndBondVects[bondAtomIndex]
    // RDKit✔️✔️:                 .second[0]
    // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
    // RDKit✔️✔️:                 ->getIdx()) {
    // RDKit✔️✔️:       std::swap(atomsAndBondVects[bondAtomIndex].second[0],
    // RDKit✔️✔️:                 atomsAndBondVects[bondAtomIndex].second[1]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    let mut result = [
        AtropEnd {
            atom: bond.begin(),
            bonds: Vec::new(),
        },
        AtropEnd {
            atom: bond.end(),
            bonds: Vec::new(),
        },
    ];
    for end in &mut result {
        end.bonds = topology
            .adjacency
            .neighbors_of(end.atom.index())
            .iter()
            .filter_map(|neighbor| (neighbor.bond != bond.id()).then_some(neighbor.bond))
            .collect();
        if end.bonds.is_empty() || end.bonds.len() > 2 {
            return None;
        }
        end.bonds
            .sort_by_key(|id| other_atom(&topology.bonds[id.index()], end.atom).index());
    }
    Some(result)
}

fn point(conformer: AtropisomerConformer<'_>, atom: AtomId) -> [f64; 3] {
    match conformer {
        AtropisomerConformer::TwoD(value) => {
            let p = value.coordinates()[atom.index()];
            [p[0], p[1], 0.0]
        }
        AtropisomerConformer::ThreeD(value) => value.coordinates()[atom.index()],
    }
}

fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn neg(value: [f64; 3]) -> [f64; 3] {
    [-value[0], -value[1], -value[2]]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn length(value: [f64; 3]) -> f64 {
    dot(value, value).sqrt()
}

fn normalized(value: [f64; 3]) -> [f64; 3] {
    let scale = 1.0 / length(value);
    [value[0] * scale, value[1] * scale, value[2] * scale]
}

fn frame_of_reference(bond: &Bond, conformer: AtropisomerConformer<'_>) -> Option<Frame> {
    // Complete pinned source: getBondFrameOfReference.
    // RDKit✔️✔️: bool getBondFrameOfReference(const Bond *bond, const Conformer *conf,
    // RDKit✔️✔️:                              RDGeom::Point3D &xAxis, RDGeom::Point3D &yAxis,
    // RDKit✔️✔️:                              RDGeom::Point3D &zAxis) {
    // RDKit✔️✔️:   // create a frame of reference that has its X-axis along the atrop bond
    // RDKit✔️✔️:   // for 2D confs, the yAxis is in the 2D plane and the zAxis is perpendicular
    // RDKit✔️✔️:   // to that plane) for 3D confs  the yAxis and the zAxis are arbitrary.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   xAxis = conf->getAtomPos(bond->getEndAtom()->getIdx()) -
    // RDKit✔️✔️:           conf->getAtomPos(bond->getBeginAtom()->getIdx());
    // RDKit✔️✔️:   if (xAxis.length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     return false;  // bond len is xero
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   xAxis.normalize();
    // RDKit✔️✔️:   if (!conf->is3D()) {
    // RDKit✔️✔️:     yAxis = RDGeom::Point3D(-xAxis.y, xAxis.x, 0);
    // RDKit✔️✔️:     yAxis.normalize();
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(0.0, 0.0, 1.0);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // here for 3D conf
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (fabs(xAxis.x) > REALLY_SMALL_BOND_LEN ||
    // RDKit✔️✔️:       fabs(xAxis.y) > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(
    // RDKit✔️✔️:         0, 0, 1);  // since X or Y value of the new x xaxis is NOT 0, this
    // RDKit✔️✔️:                    // new temp z axis cannnot be colinear with it
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(1, 0, 0);  // since the new x axis is exactly along
    // RDKit✔️✔️:                                        // the (old) z axis, this new temp z axis
    // RDKit✔️✔️:                                        // is NOT colinear with it
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   yAxis = zAxis.crossProduct(xAxis);
    // RDKit✔️✔️:   zAxis = xAxis.crossProduct(yAxis);
    // RDKit✔️✔️:   yAxis.normalize();
    // RDKit✔️✔️:   zAxis.normalize();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    let axis = sub(point(conformer, bond.end()), point(conformer, bond.begin()));
    if length(axis) < REALLY_SMALL_BOND_LEN {
        return None;
    }
    let x = normalized(axis);
    if matches!(conformer, AtropisomerConformer::TwoD(_)) {
        return Some(Frame {
            y: normalized([-x[1], x[0], 0.0]),
            z: [0.0, 0.0, 1.0],
        });
    }
    let initial_z = if x[0].abs() > REALLY_SMALL_BOND_LEN || x[1].abs() > REALLY_SMALL_BOND_LEN {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let y = normalized(cross(initial_z, x));
    let z = normalized(cross(x, y));
    Some(Frame { y, z })
}

fn end_vector(
    topology: &TopologyBlock,
    end: &AtropEnd,
    frame: Frame,
    conformer: AtropisomerConformer<'_>,
) -> Result<[f64; 3], AtropisomerRejectionKind> {
    // Complete pinned source: getAtropIsomerEndVect.
    // RDKit✔️✔️: bool getAtropIsomerEndVect(const AtropAtomAndBondVec &atomAndBondVec,
    // RDKit✔️✔️:                            const RDGeom::Point3D &yAxis,
    // RDKit✔️✔️:                            const RDGeom::Point3D &zAxis, const Conformer *conf,
    // RDKit✔️✔️:                            RDGeom::Point3D &bondVec) {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       atomAndBondVec.second.size() > 0 && atomAndBondVec.second.size() < 3,
    // RDKit✔️✔️:       "bad bond size");
    // RDKit✔️✔️:   PRECONDITION(atomAndBondVec.second[0], "bad first bond");
    // RDKit✔️✔️:   PRECONDITION(atomAndBondVec.second.size() == 1 || atomAndBondVec.second[1],
    // RDKit✔️✔️:                "bad second bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec = conf->getAtomPos(atomAndBondVec.second[0]
    // RDKit✔️✔️:                                  ->getOtherAtom(atomAndBondVec.first)
    // RDKit✔️✔️:                                  ->getIdx()) -
    // RDKit✔️✔️:             conf->getAtomPos(
    // RDKit✔️✔️:                 atomAndBondVec.first->getIdx());  // in old frame of reference
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec = RDGeom::Point3D(0.0, bondVec.dotProduct(yAxis),
    // RDKit✔️✔️:                             bondVec.dotProduct(zAxis));  // in new frame
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // make sure the other atom is on the other side
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atomAndBondVec.second.size() == 2) {
    // RDKit✔️✔️:     RDGeom::Point3D otherVec =
    // RDKit✔️✔️:         conf->getAtomPos(atomAndBondVec.second[1]
    // RDKit✔️✔️:                              ->getOtherAtom(atomAndBondVec.first)
    // RDKit✔️✔️:                              ->getIdx()) -
    // RDKit✔️✔️:         conf->getAtomPos(
    // RDKit✔️✔️:             atomAndBondVec.first->getIdx());  // in old frame of reference
    // RDKit✔️✔️:     otherVec = RDGeom::Point3D(0.0, otherVec.dotProduct(yAxis),
    // RDKit✔️✔️:                                otherVec.dotProduct(zAxis));  // in new frame
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:       bondVec = -otherVec;  // put it on the other side of otherVec
    // RDKit✔️✔️:     } else if (bondVec.dotProduct(otherVec) > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:       // the product of dotproducts (y-values) should be
    // RDKit✔️✔️:       // negative (or at least zero)
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Both bonds on one end of an atropisomer are on the same side - atoms is : "
    // RDKit✔️✔️:           << atomAndBondVec.first->getIdx() << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Could not find a bond on one end of an atropisomer that is not co-linear - atoms are : "
    // RDKit✔️✔️:         << atomAndBondVec.first->getIdx() << std::endl;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec.normalize();
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    let projected = |bond_id: BondId| {
        let carrier = &topology.bonds[bond_id.index()];
        let vector = sub(
            point(conformer, other_atom(carrier, end.atom)),
            point(conformer, end.atom),
        );
        [0.0, dot(vector, frame.y), dot(vector, frame.z)]
    };
    let mut result = projected(end.bonds[0]);
    if end.bonds.len() == 2 {
        let other = projected(end.bonds[1]);
        if length(result) < REALLY_SMALL_BOND_LEN {
            result = neg(other);
        } else if dot(result, other) > REALLY_SMALL_BOND_LEN {
            return Err(AtropisomerRejectionKind::SameSideCarriers);
        }
    }
    if length(result) < REALLY_SMALL_BOND_LEN {
        return Err(AtropisomerRejectionKind::CollinearCarrier);
    }
    Ok(normalized(result))
}

fn effective_direction(
    topology: &TopologyBlock,
    updates: &BTreeMap<BondId, AtropisomerWedgeUpdate>,
    bond: BondId,
) -> BondDirection {
    updates
        .get(&bond)
        .map_or(topology.bonds[bond.index()].direction(), |value| {
            value.direction
        })
}

fn interpreted_end_direction(
    topology: &TopologyBlock,
    end: &AtropEnd,
    updates: &BTreeMap<BondId, AtropisomerWedgeUpdate>,
) -> Result<BondDirection, AtropisomerRejectionKind> {
    // Complete pinned source: getBondDir.
    // RDKit✔️✔️: std::pair<bool, Bond::BondDir> getBondDir(
    // RDKit✔️✔️:     const Bond *bond, const AtropAtomAndBondVec &atomAndBondVec) {
    // RDKit✔️✔️:   // get the wedge dir for this end of the bond
    // RDKit✔️✔️:   // if the first bond 1 has a bondDir, use it
    // RDKit✔️✔️:   // if the second bond has a bond dir use the opposite of if
    // RDKit✔️✔️:   // if both bonds have a dir, make sure they are different
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bond1Dir = atomAndBondVec.second[0]->getBondDir();
    // RDKit✔️✔️:   if (bond1Dir != Bond::BEGINWEDGE && bond1Dir != Bond::BEGINDASH) {
    // RDKit✔️✔️:     bond1Dir = Bond::NONE;  //  we dont care if it any thing else
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto bond2Dir = atomAndBondVec.second.size() == 2
    // RDKit✔️✔️:                       ? atomAndBondVec.second[1]->getBondDir()
    // RDKit✔️✔️:                       : Bond::NONE;
    // RDKit✔️✔️:   if (bond2Dir != Bond::BEGINWEDGE && bond2Dir != Bond::BEGINDASH) {
    // RDKit✔️✔️:     bond2Dir = Bond::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if both are set to a direction, they must NOT be the same - one
    // RDKit✔️✔️:   // must be a dash and the other a hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bond1Dir != Bond::NONE && bond2Dir != Bond::NONE &&
    // RDKit✔️✔️:       bond1Dir == bond2Dir) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "The bonds on one end of an atropisomer are both UP or both DOWN - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return {false, Bond::BondDir::NONE};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bond1Dir == Bond::BEGINWEDGE || bond2Dir == Bond::BEGINDASH) {
    // RDKit✔️✔️:     return {true, Bond::BondDir::BEGINWEDGE};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bond1Dir == Bond::BEGINDASH || bond2Dir == Bond::BEGINWEDGE) {
    // RDKit✔️✔️:     return {true, Bond::BondDir::BEGINDASH};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return {true, Bond::BondDir::NONE};
    // RDKit✔️✔️: }
    let normalize = |value| match value {
        BondDirection::BeginWedge | BondDirection::BeginDash => value,
        _ => BondDirection::None,
    };
    let first = normalize(effective_direction(topology, updates, end.bonds[0]));
    let second = if end.bonds.len() == 2 {
        normalize(effective_direction(topology, updates, end.bonds[1]))
    } else {
        BondDirection::None
    };
    if first != BondDirection::None && first == second {
        return Err(AtropisomerRejectionKind::InconsistentDirections);
    }
    if first == BondDirection::BeginWedge || second == BondDirection::BeginDash {
        Ok(BondDirection::BeginWedge)
    } else if first == BondDirection::BeginDash || second == BondDirection::BeginWedge {
        Ok(BondDirection::BeginDash)
    } else {
        Ok(BondDirection::None)
    }
}

fn validate_conformer(
    conformer: AtropisomerConformer<'_>,
    atom_count: usize,
) -> Result<(), AtropisomerError> {
    match conformer {
        AtropisomerConformer::TwoD(value) => value
            .validate_for_atom_count(atom_count)
            .map_err(|source| AtropisomerError::InvalidCoordinates { source }),
        AtropisomerConformer::ThreeD(value) => {
            value
                .validate_for_atom_count(atom_count)
                .map_err(|source| AtropisomerError::InvalidCoordinates { source })?;
            if !value.is_3d() {
                return Err(AtropisomerError::ConformerNotThreeDimensional {
                    conformer: value.id(),
                });
            }
            Ok(())
        }
    }
}

fn detect_one(
    topology: &TopologyBlock,
    bond: &Bond,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<BondStereo, AtropisomerRejectionKind> {
    // Complete pinned source: DetectAtropisomerChiralityOneBond.
    // RDKit✔️✔️: void DetectAtropisomerChiralityOneBond(Bond *bond, ROMol &mol,
    // RDKit✔️✔️:                                        const Conformer *conf) {
    // RDKit✔️✔️:   // the approach is this:
    // RDKit✔️✔️:   // we will view the system along the line from the potential atropisomer
    // RDKit✔️✔️:   // bond, from atom1 to atom 2 and we do a coordinate transformation to
    // RDKit✔️✔️:   // the plane of reference where that vector, from a1 to a2, is the x-AXIS.
    // RDKit✔️✔️:   // For 2D, the y axis is in the 2D plane, and the zaxis is perpendicaul to
    // RDKit✔️✔️:   // the 2D plane For 3D, the Y and Z axes are taken arbitrarily to form
    // RDKit✔️✔️:   // a right-handed system with the X-axis.
    // RDKit✔️✔️:   //  atoms 1 and 2 each have one or two bonds out from the main potential
    // RDKit✔️✔️:   //  atrop bond. for each end of the main bond, we find a vector to reprent
    // RDKit✔️✔️:   //  the neighbor atom with the smallest index as its projection onto the
    // RDKit✔️✔️:   //  x=0 plane.
    // RDKit✔️✔️:   // (In 2d, this projection is on the y-AXIS for the end that does NOT have
    // RDKit✔️✔️:   // a wedge/hash bond, and  on the z axis - out of the plane - for the end
    // RDKit✔️✔️:   // that does have a wedge/hash). The chirality is recorded as the
    // RDKit✔️✔️:   // direction we rotate from, atom 1's projection to atom2's proejection -
    // RDKit✔️✔️:   // either clockwise or counter clockwise
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // one vector for each end - each one - should end up with 1 or 2 entries
    // RDKit✔️✔️:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit✔️✔️:     return;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // make sure we do not have wiggle bonds
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atomAndBondVec : atomAndBondVecs) {
    // RDKit✔️✔️:     for (auto endBond : atomAndBondVec.second) {
    // RDKit✔️✔️:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         return;  // not an atropisomer
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // the convention is that in the absence of coords, the coordiates are choosen
    // RDKit✔️✔️:   // with the lowest numbered atom of the atrop bond down, and the other atom
    // RDKit✔️✔️:   // straight up.
    // RDKit✔️✔️:   // On each end, the lowest numbered connecting atom is on the left
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //              a      b
    // RDKit✔️✔️:   //               \   /
    // RDKit✔️✔️:   //                 c
    // RDKit✔️✔️:   //                 |
    // RDKit✔️✔️:   //                 d
    // RDKit✔️✔️:   //               /   \     aaa
    // RDKit✔️✔️:   //              e      f
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // where  c > d
    // RDKit✔️✔️:   //        a < b
    // RDKit✔️✔️:   //        e < f
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (conf == nullptr) {
    // RDKit✔️✔️:     std::pair<bool, Bond::BondDir> bond1DirResult;
    // RDKit✔️✔️:     bond1DirResult = getBondDir(bond, atomAndBondVecs[0]);
    // RDKit✔️✔️:     if (!bond1DirResult.first) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::pair<bool, Bond::BondDir> bond2DirResult;
    // RDKit✔️✔️:     bond2DirResult = getBondDir(bond, atomAndBondVecs[1]);
    // RDKit✔️✔️:     if (!bond2DirResult.first) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (bond1DirResult.second == bond2DirResult.second) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "inconsistent bond wedging for an atropisomer.  Atoms are: "
    // RDKit✔️✔️:           << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bond1DirResult.second == Bond::BEGINWEDGE ||
    // RDKit✔️✔️:         bond2DirResult.second == Bond::BEGINDASH) {
    // RDKit✔️✔️:       bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
    // RDKit✔️✔️:     } else if (bond1DirResult.second == Bond::BEGINDASH ||
    // RDKit✔️✔️:                bond2DirResult.second == Bond::BEGINWEDGE) {
    // RDKit✔️✔️:       bond->setStereo(Bond::BondStereo::STEREOATROPCW);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // create a frame of reference that has its X-axis along the atrop bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D xAxis, yAxis, zAxis;
    // RDKit✔️✔️:   if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // RDKit✔️✔️:     // connot percieve atroisomer
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Failed to get a frame of reference along an atropisomer bond - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
    // RDKit✔️✔️:                                 // potential atropisomer bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit✔️✔️:     // if the conf is 2D, we use the wedge bonds to set the coords for the
    // RDKit✔️✔️:     // projected vector onto the xAxis perpendicular plane (looking down
    // RDKit✔️✔️:     // the atrop bond )
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!conf->is3D()) {
    // RDKit✔️✔️:       // get the wedge dir for this end of the bond
    // RDKit✔️✔️:       // if the first bond 1 has a bondDir, use it
    // RDKit✔️✔️:       // if the second bond has a bond dir use the opposite of if
    // RDKit✔️✔️:       // if both bonds have a dir, make sure they are different
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::pair<bool, Bond::BondDir> bondDirResult;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       bondDirResult = getBondDir(bond, atomAndBondVecs[bondAtomIndex]);
    // RDKit✔️✔️:       if (!bondDirResult.first) {
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!getAtropIsomerEndVect(atomAndBondVecs[bondAtomIndex], yAxis, zAxis,
    // RDKit✔️✔️:                                  conf, bondVecs[bondAtomIndex])) {
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondDirResult.second == Bond::BEGINWEDGE) {
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].y *= 0.707;
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].z = fabs(bondVecs[bondAtomIndex].y);
    // RDKit✔️✔️:       } else if (bondDirResult.second == Bond::BEGINDASH) {
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].y *= 0.707;
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].z = -fabs(bondVecs[bondAtomIndex].y);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {  // the conf is 3D
    // RDKit✔️✔️:       // to be considered, one or more neighbor bonds must have a wedge or
    // RDKit✔️✔️:       // hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // find the projection of the bond(s) on this end in the frame of
    // RDKit✔️✔️:       // reference's  x=0  plane
    // RDKit✔️✔️:       RDGeom::Point3D tempBondVec =
    // RDKit✔️✔️:           conf->getAtomPos(
    // RDKit✔️✔️:               atomAndBondVecs[bondAtomIndex]
    // RDKit✔️✔️:                   .second[0]
    // RDKit✔️✔️:                   ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit✔️✔️:                   ->getIdx()) -
    // RDKit✔️✔️:           conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
    // RDKit✔️✔️:       bondVecs[bondAtomIndex] = RDGeom::Point3D(
    // RDKit✔️✔️:           0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (atomAndBondVecs[bondAtomIndex].second.size() == 2) {
    // RDKit✔️✔️:         tempBondVec =
    // RDKit✔️✔️:             conf->getAtomPos(
    // RDKit✔️✔️:                 atomAndBondVecs[bondAtomIndex]
    // RDKit✔️✔️:                     .second[1]
    // RDKit✔️✔️:                     ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit✔️✔️:                     ->getIdx()) -
    // RDKit✔️✔️:             conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // get the projection of the 2nd bond on the x=0 plane
    // RDKit✔️✔️:
    // RDKit✔️✔️:         RDGeom::Point3D otherBondVec = RDGeom::Point3D(
    // RDKit✔️✔️:             0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if the first atom is co-linear with the main atrop bond, use
    // RDKit✔️✔️:         // the opposite of the 2nd atom
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:           bondVecs[bondAtomIndex] =
    // RDKit✔️✔️:               -otherBondVec;  // note - it might still be co-linear-
    // RDKit✔️✔️:                               // this is checked below
    // RDKit✔️✔️:         } else if (bondVecs[bondAtomIndex].dotProduct(otherBondVec) >
    // RDKit✔️✔️:                    REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "Both bonds on one end of an atropisomer are on the same side - atoms are: "
    // RDKit✔️✔️:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:           return;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Failed to find a bond on one end of an atropisomer that is NOT co-linear - atoms are: "
    // RDKit✔️✔️:             << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto crossProduct = bondVecs[1].crossProduct(bondVecs[0]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (crossProduct.x > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
    // RDKit✔️✔️:   } else if (crossProduct.x < -REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREOATROPCW);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "The 2 defining bonds for an atropisomer are co-planar - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let ends = atropisomer_ends(topology, bond).ok_or(AtropisomerRejectionKind::MissingCarrier)?;
    if ends
        .iter()
        .flat_map(|end| &end.bonds)
        .any(|id| topology.bonds[id.index()].direction() == BondDirection::Unknown)
    {
        return Err(AtropisomerRejectionKind::UnknownCarrierDirection);
    }
    let no_updates = BTreeMap::new();
    if conformer.is_none() {
        let first = interpreted_end_direction(topology, &ends[0], &no_updates)?;
        let second = interpreted_end_direction(topology, &ends[1], &no_updates)?;
        if first == second {
            return Err(AtropisomerRejectionKind::InconsistentDirections);
        }
        return if first == BondDirection::BeginWedge || second == BondDirection::BeginDash {
            Ok(BondStereo::AtropCcw)
        } else if first == BondDirection::BeginDash || second == BondDirection::BeginWedge {
            Ok(BondStereo::AtropCw)
        } else {
            Err(AtropisomerRejectionKind::InconsistentDirections)
        };
    }
    let conformer = conformer.expect("checked above");
    let frame =
        frame_of_reference(bond, conformer).ok_or(AtropisomerRejectionKind::ZeroLengthAxis)?;
    let mut vectors = [
        end_vector(topology, &ends[0], frame, conformer)?,
        end_vector(topology, &ends[1], frame, conformer)?,
    ];
    if matches!(conformer, AtropisomerConformer::TwoD(_)) {
        for (index, end) in ends.iter().enumerate() {
            match interpreted_end_direction(topology, end, &no_updates)? {
                BondDirection::BeginWedge => {
                    vectors[index][1] *= 0.707;
                    vectors[index][2] = vectors[index][1].abs();
                }
                BondDirection::BeginDash => {
                    vectors[index][1] *= 0.707;
                    vectors[index][2] = -vectors[index][1].abs();
                }
                _ => {}
            }
        }
    }
    let x = cross(vectors[1], vectors[0])[0];
    if x > REALLY_SMALL_BOND_LEN {
        Ok(BondStereo::AtropCcw)
    } else if x < -REALLY_SMALL_BOND_LEN {
        Ok(BondStereo::AtropCw)
    } else {
        Err(AtropisomerRejectionKind::CoplanarCarriers)
    }
}

pub fn detect_atropisomer_chirality(
    topology: &TopologyBlock,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<AtropisomerAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    // Complete pinned source: detectAtropisomerChirality.
    // RDKit✔️✔️: void detectAtropisomerChirality(ROMol &mol, const Conformer *conf) {
    // RDKit✔️✔️:   PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
    // RDKit✔️✔️:                "conformer does not belong to molecule");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::set<Bond *> bondsToTry;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (canHaveDirection(*bond) &&
    // RDKit✔️✔️:         (bond->getBondDir() == Bond::BondDir::BEGINDASH ||
    // RDKit✔️✔️:          bond->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
    // RDKit✔️✔️:       for (const auto &nbrBond : mol.atomBonds(bond->getBeginAtom())) {
    // RDKit✔️✔️:         if (nbrBond == bond) {
    // RDKit✔️✔️:           continue;  // a bond is NOT its own neighbor
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         bondsToTry.insert(nbrBond);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bondsToTry.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // First, do a simple check with TotalDegree to see if any bonds might be
    // RDKit✔️✔️:   // candidates before doing the expensive hybridization calculation.
    // RDKit✔️✔️:   bool anyBondPassesDegreeCheck = false;
    // RDKit✔️✔️:   for (auto bondToTry : bondsToTry) {
    // RDKit✔️✔️:     if (bondToTry->getBeginAtom()->needsUpdatePropertyCache()) {
    // RDKit✔️✔️:       bondToTry->getBeginAtom()->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (bondToTry->getEndAtom()->needsUpdatePropertyCache()) {
    // RDKit✔️✔️:       bondToTry->getEndAtom()->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bondToTry->getBondType() == Bond::SINGLE &&
    // RDKit✔️✔️:         bondToTry->getStereo() != Bond::BondStereo::STEREOANY &&
    // RDKit✔️✔️:         bondToTry->getBeginAtom()->getTotalDegree() >= 2 &&
    // RDKit✔️✔️:         bondToTry->getBeginAtom()->getTotalDegree() <= 3 &&
    // RDKit✔️✔️:         bondToTry->getEndAtom()->getTotalDegree() >= 2 &&
    // RDKit✔️✔️:         bondToTry->getEndAtom()->getTotalDegree() <= 3) {
    // RDKit✔️✔️:       anyBondPassesDegreeCheck = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!anyBondPassesDegreeCheck) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // defer cache update on the whole mol unless we actually have bonds to try
    // RDKit✔️✔️:   // we need to do an update on the whole mol and not just incident atoms
    // RDKit✔️✔️:   // because we need to calculate hybridization, which is non-local
    // RDKit✔️✔️:   bool needsUpdate =
    // RDKit✔️✔️:       mol.needsUpdatePropertyCache() ||
    // RDKit✔️✔️:       std::any_of(mol.atoms().begin(), mol.atoms().end(), [](const auto atom) {
    // RDKit✔️✔️:         return atom->getAtomicNum() != 0 &&
    // RDKit✔️✔️:                atom->getHybridization() == Atom::HybridizationType::UNSPECIFIED;
    // RDKit✔️✔️:       });
    // RDKit✔️✔️:   if (needsUpdate) {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:     MolOps::setConjugation(mol);
    // RDKit✔️✔️:     MolOps::setHybridization(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto bondToTry : bondsToTry) {
    // RDKit✔️✔️:     if (bondToTry->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:         bondToTry->getStereo() == Bond::BondStereo::STEREOANY ||
    // RDKit✔️✔️:         // before, we checked only on totalDegree = 2 or 3,
    // RDKit✔️✔️:         // but this causes false positives for something like a chiral sulfoxide
    // RDKit✔️✔️:         // since the S is tetrahedral (sp3) but has only 3 substituents.
    // RDKit✔️✔️:         // the hybridization code relies on totalDegree,
    // RDKit✔️✔️:         // but modified to include and making sure to include conjugation
    // RDKit✔️✔️:         // so while this is more expensive per molecule, it is closer to intent
    // RDKit✔️✔️:         bondToTry->getBeginAtom()->getHybridization() != Atom::SP2 ||
    // RDKit✔️✔️:         bondToTry->getEndAtom()->getHybridization() != Atom::SP2) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     DetectAtropisomerChiralityOneBond(bondToTry, mol, conf);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut candidates = BTreeSet::new();
    for marker in &topology.bonds {
        if can_have_direction(marker)
            && matches!(
                marker.direction(),
                BondDirection::BeginDash | BondDirection::BeginWedge
            )
        {
            for neighbor in topology.adjacency.neighbors_of(marker.begin().index()) {
                if neighbor.bond != marker.id() {
                    candidates.insert(neighbor.bond);
                }
            }
        }
    }
    let any_degree_candidate = candidates.iter().copied().any(|id| {
        let bond = &topology.bonds[id.index()];
        bond.order() == BondOrder::Single
            && bond.stereo() != BondStereo::Any
            && (2..=3).contains(&total_degree(topology, bond.begin()))
            && (2..=3).contains(&total_degree(topology, bond.end()))
    });
    if !any_degree_candidate {
        return Ok(AtropisomerAssignment::default());
    }
    let mut assignment = AtropisomerAssignment::default();
    for id in candidates {
        let bond = &topology.bonds[id.index()];
        if bond.order() != BondOrder::Single
            || bond.stereo() == BondStereo::Any
            || topology.atoms[bond.begin().index()].hybridization() != Hybridization::Sp2
            || topology.atoms[bond.end().index()].hybridization() != Hybridization::Sp2
        {
            continue;
        }
        match detect_one(topology, bond, conformer) {
            Ok(stereo) => assignment
                .bond_updates
                .push(AtropisomerBondUpdate { bond: id, stereo }),
            Err(kind) => assignment
                .diagnostics
                .push(AtropisomerDiagnostic { bond: id, kind }),
        }
    }
    Ok(assignment)
}

fn effective_stereo(
    topology: &TopologyBlock,
    assignment: &AtropisomerAssignment,
    bond: BondId,
) -> BondStereo {
    assignment
        .bond_updates
        .iter()
        .find(|update| update.bond == bond)
        .map_or(topology.bonds[bond.index()].stereo(), |update| {
            update.stereo
        })
}

pub fn cleanup_atropisomer_stereo_groups(
    topology: &TopologyBlock,
    detected: &AtropisomerAssignment,
) -> Result<StereoGroupAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    for update in &detected.bond_updates {
        if update.bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::AssignmentBondOutOfRange {
                bond: update.bond,
                bond_count: topology.bonds.len(),
            });
        }
        if !matches!(update.stereo, BondStereo::AtropCw | BondStereo::AtropCcw) {
            return Err(AtropisomerError::InvalidAtropisomerStereo {
                bond: update.bond,
                stereo: update.stereo,
            });
        }
    }
    // Complete pinned source: cleanupAtropisomerStereoGroups.
    // RDKit✔️✔️: void cleanupAtropisomerStereoGroups(ROMol &mol) {
    // RDKit✔️✔️:   std::vector<StereoGroup> newsgs;
    // RDKit✔️✔️:   for (auto sg : mol.getStereoGroups()) {
    // RDKit✔️✔️:     std::vector<Atom *> okatoms;
    // RDKit✔️✔️:     std::vector<Bond *> okbonds;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (auto atom : sg.getAtoms()) {
    // RDKit✔️✔️:       bool foundAtrop = false;
    // RDKit✔️✔️:       for (auto bndI : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:         auto bond = (mol)[bndI];
    // RDKit✔️✔️:         if (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
    // RDKit✔️✔️:             bond->getStereo() == Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️✔️:           foundAtrop = true;
    // RDKit✔️✔️:           if (std::find(okbonds.begin(), okbonds.end(), bond) ==
    // RDKit✔️✔️:               okbonds.end()) {
    // RDKit✔️✔️:             okbonds.push_back(bond);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!foundAtrop) {
    // RDKit✔️✔️:         okatoms.push_back(atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (okbonds.empty()) {
    // RDKit✔️✔️:       newsgs.push_back(sg);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       newsgs.emplace_back(sg.getGroupType(), std::move(okatoms),
    // RDKit✔️✔️:                           std::move(okbonds));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setStereoGroups(std::move(newsgs));
    // RDKit✔️✔️: }
    let mut groups = Vec::with_capacity(topology.stereo_groups.len());
    for group in &topology.stereo_groups {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        for atom in group.atoms() {
            let mut found = false;
            for neighbor in topology.adjacency.neighbors_of(atom.index()) {
                if matches!(
                    effective_stereo(topology, detected, neighbor.bond),
                    BondStereo::AtropCw | BondStereo::AtropCcw
                ) {
                    found = true;
                    if !bonds.contains(&neighbor.bond) {
                        bonds.push(neighbor.bond);
                    }
                }
            }
            if !found {
                atoms.push(*atom);
            }
        }
        if bonds.is_empty() {
            groups.push(group.clone());
        } else {
            let mut replacement = StereoGroup::new(group.kind(), atoms, bonds);
            if let Some(id) = group.id() {
                replacement = replacement.with_id(id);
            }
            groups.push(replacement);
        }
    }
    Ok(StereoGroupAssignment { groups })
}

pub fn stereo_group_atom_ids(
    topology: &TopologyBlock,
    group: &StereoGroup,
    wedges: &AtropisomerWedgeAssignment,
) -> Result<Vec<AtomId>, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    for atom in group.atoms() {
        if atom.index() >= topology.atoms.len() {
            return Err(AtropisomerError::StereoGroupAtomOutOfRange {
                atom: *atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    for bond in group.bonds() {
        if bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::StereoGroupBondOutOfRange {
                bond: *bond,
                bond_count: topology.bonds.len(),
            });
        }
    }
    // Complete pinned source: getAllAtomIdsForStereoGroup.
    // RDKit✔️❌: void getAllAtomIdsForStereoGroup(
    // RDKit✔️❌:     const ROMol &mol, const StereoGroup &group,
    // RDKit✔️❌:     std::vector<unsigned int> &atomIds,
    // RDKit✔️❌:     const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit✔️❌:         &wedgeBonds) {
    // RDKit✔️❌:   atomIds.clear();
    // RDKit✔️❌:   for (auto &&atom : group.getAtoms()) {
    // RDKit✔️❌:     atomIds.push_back(atom->getIdx());
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto &&bond : group.getBonds()) {
    // RDKit✔️❌:     // figure out which atoms of the bond get wedge/hash indications
    // RDKit✔️❌:     // mark the atom with the wedge/hash
    // RDKit✔️❌:
    // RDKit✔️❌:     for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
    // RDKit✔️❌:       for (const auto atomBond : mol.atomBonds(atom)) {
    // RDKit✔️❌:         if (atomBond->getIdx() == bond->getIdx()) {
    // RDKit✔️❌:           continue;
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         if (atomBond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit✔️❌:             atomBond->getBondDir() == Bond::BEGINDASH ||
    // RDKit✔️❌:             (wedgeBonds.find(atomBond->getIdx()) != wedgeBonds.end() &&
    // RDKit✔️❌:              (wedgeBonds.at(atomBond->getIdx())->getType()) ==
    // RDKit✔️❌:                  Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer)) {
    // RDKit✔️❌:           if (std::find(atomIds.begin(), atomIds.end(), atom->getIdx()) ==
    // RDKit✔️❌:               atomIds.end()) {
    // RDKit✔️❌:             atomIds.push_back(atom->getIdx());
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Rust preserves the modeled behavior but linearly scans wedge updates for each
    // adjacent bond instead of using the source map lookup, so the local complexity
    // review records the known performance gap.
    let mut ids = group.atoms().to_vec();
    for axial_id in group.bonds() {
        let axial = &topology.bonds[axial_id.index()];
        for atom in [axial.begin(), axial.end()] {
            let marked = topology
                .adjacency
                .neighbors_of(atom.index())
                .iter()
                .filter(|neighbor| neighbor.bond != *axial_id)
                .any(|neighbor| {
                    matches!(
                        topology.bonds[neighbor.bond.index()].direction(),
                        BondDirection::BeginWedge | BondDirection::BeginDash
                    ) || wedges.bond_updates.iter().any(|update| {
                        update.bond == neighbor.bond && update.atropisomer_bond == *axial_id
                    })
                });
            if marked && !ids.contains(&atom) {
                ids.push(atom);
            }
        }
    }
    Ok(ids)
}

fn no_conf_direction(stereo: BondStereo, which_end: usize, which_bond: usize) -> BondDirection {
    // Complete pinned source: getBondDirForAtropisomerNoConf.
    // RDKit✔️✔️: Bond::BondDir getBondDirForAtropisomerNoConf(Bond::BondStereo bondStereo,
    // RDKit✔️✔️:                                              unsigned int whichEnd,
    // RDKit✔️✔️:                                              unsigned int whichBond) {
    // RDKit✔️✔️:   // the convention is that in the absence of coords, the coordiates are choosen
    // RDKit✔️✔️:   // with the lowest numbered atom of the atrop bond down, and the other atom
    // RDKit✔️✔️:   // straight up.
    // RDKit✔️✔️:   // On each end, the lowest numbered connecting atom is on the left
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //              a      b
    // RDKit✔️✔️:   //               \   /
    // RDKit✔️✔️:   //                 c
    // RDKit✔️✔️:   //                 |
    // RDKit✔️✔️:   //                 d
    // RDKit✔️✔️:   //               /   \     aaa
    // RDKit✔️✔️:   //              e      f
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // where  c > d
    // RDKit✔️✔️:   //        a < b
    // RDKit✔️✔️:   //        e < f
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
    // RDKit✔️✔️:                    bondStereo == Bond::BondStereo::STEREOATROPCCW,
    // RDKit✔️✔️:                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int flips = 0;
    // RDKit✔️✔️:   if (bondStereo == Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichBond == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichEnd == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return flips % 2 ? Bond::BEGINDASH : Bond::BEGINWEDGE;
    // RDKit✔️✔️: }
    let flips = usize::from(stereo == BondStereo::AtropCw) + which_bond + which_end;
    if flips % 2 == 1 {
        BondDirection::BeginDash
    } else {
        BondDirection::BeginWedge
    }
}

fn two_d_direction(
    vectors: [[f64; 3]; 2],
    stereo: BondStereo,
    which_end: usize,
    which_bond: usize,
) -> BondDirection {
    // Complete pinned source: getBondDirForAtropisomer2d.
    // RDKit✔️✔️: Bond::BondDir getBondDirForAtropisomer2d(RDGeom::Point3D bondVecs[2],
    // RDKit✔️✔️:                                          Bond::BondStereo bondStereo,
    // RDKit✔️✔️:                                          unsigned int whichEnd,
    // RDKit✔️✔️:                                          unsigned int whichBond) {
    // RDKit✔️✔️:   PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
    // RDKit✔️✔️:                    bondStereo == Bond::BondStereo::STEREOATROPCCW,
    // RDKit✔️✔️:                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int flips = 0;
    // RDKit✔️✔️:   if (bondStereo == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichBond == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichEnd == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bondVecs[1 - whichEnd].y < 0) {
    // RDKit✔️✔️:     ++flips;  // if the OTHER end is negative for the low index bond vec, it
    // RDKit✔️✔️:               // is a flip
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return flips % 2 ? Bond::BEGINWEDGE : Bond::BEGINDASH;
    // RDKit✔️✔️: }
    let flips = usize::from(stereo == BondStereo::AtropCcw)
        + which_bond
        + which_end
        + usize::from(vectors[1 - which_end][1] < 0.0);
    if flips % 2 == 1 {
        BondDirection::BeginWedge
    } else {
        BondDirection::BeginDash
    }
}

fn three_d_direction(
    topology: &TopologyBlock,
    bond: BondId,
    center: AtomId,
    conformer: AtropisomerConformer<'_>,
) -> BondDirection {
    // Complete pinned source: getBondDirForAtropisomer3d.
    // RDKit✔️✔️: Bond::BondDir getBondDirForAtropisomer3d(Bond *whichBond,
    // RDKit✔️✔️:                                          const Conformer *conf) {
    // RDKit✔️✔️:   // for 3D we mark it as wedge or hash depending on the z-value of the bond
    // RDKit✔️✔️:   // vector
    // RDKit✔️✔️:   //  IT really doesn't matter since we ignore these except as MARKERS for
    // RDKit✔️✔️:   //  which bonds are atropisomer bonds
    // RDKit✔️✔️:   if ((conf->getAtomPos(whichBond->getEndAtom()->getIdx()).z -
    // RDKit✔️✔️:        conf->getAtomPos(whichBond->getBeginAtom()->getIdx()).z) >
    // RDKit✔️✔️:       REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     return Bond::BondDir::BEGINWEDGE;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return Bond::BondDir::BEGINDASH;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let value = &topology.bonds[bond.index()];
    let other = other_atom(value, center);
    if point(conformer, other)[2] - point(conformer, center)[2] > REALLY_SMALL_BOND_LEN {
        BondDirection::BeginWedge
    } else {
        BondDirection::BeginDash
    }
}

#[derive(Clone, Copy)]
enum WedgeMode<'a> {
    NoConformer,
    TwoD { vectors: [[f64; 3]; 2] },
    ThreeD(AtropisomerConformer<'a>),
}

fn desired_direction(
    topology: &TopologyBlock,
    axial: &Bond,
    mode: WedgeMode<'_>,
    end: usize,
    carrier: usize,
    carrier_id: BondId,
    center: AtomId,
) -> BondDirection {
    match mode {
        WedgeMode::NoConformer => no_conf_direction(axial.stereo(), end, carrier),
        WedgeMode::TwoD { vectors, .. } => two_d_direction(vectors, axial.stereo(), end, carrier),
        WedgeMode::ThreeD(conformer) => three_d_direction(topology, carrier_id, center, conformer),
    }
}

fn record_update(
    topology: &TopologyBlock,
    updates: &mut BTreeMap<BondId, AtropisomerWedgeUpdate>,
    carrier: BondId,
    center: AtomId,
    direction: BondDirection,
    axial: BondId,
) {
    let bond = &topology.bonds[carrier.index()];
    let other = other_atom(bond, center);
    updates.insert(
        carrier,
        AtropisomerWedgeUpdate {
            bond: carrier,
            begin: center,
            end: other,
            direction,
            atropisomer_bond: axial,
        },
    );
}

fn wedge_one(
    topology: &TopologyBlock,
    rings: &RingInfo,
    axial: &Bond,
    conformer: Option<AtropisomerConformer<'_>>,
    occupied: &BTreeSet<BondId>,
    updates: &mut BTreeMap<BondId, AtropisomerWedgeUpdate>,
) -> Result<(), AtropisomerRejectionKind> {
    let ends = atropisomer_ends(topology, axial).ok_or(AtropisomerRejectionKind::MissingCarrier)?;
    if ends
        .iter()
        .flat_map(|end| &end.bonds)
        .any(|id| effective_direction(topology, updates, *id) == BondDirection::Unknown)
    {
        return Err(AtropisomerRejectionKind::UnknownCarrierDirection);
    }
    let mode = match conformer {
        None => WedgeMode::NoConformer,
        Some(value @ AtropisomerConformer::TwoD(_)) => {
            let frame =
                frame_of_reference(axial, value).ok_or(AtropisomerRejectionKind::ZeroLengthAxis)?;
            WedgeMode::TwoD {
                vectors: [
                    end_vector(topology, &ends[0], frame, value)?,
                    end_vector(topology, &ends[1], frame, value)?,
                ],
            }
        }
        Some(value @ AtropisomerConformer::ThreeD(_)) => WedgeMode::ThreeD(value),
    };
    // The three source functions first reuse every wedge/hash carrier whose
    // narrow end is the axial endpoint.
    let mut existing = Vec::new();
    for (end_index, end) in ends.iter().enumerate() {
        for (carrier_index, carrier) in end.bonds.iter().copied().enumerate() {
            let bond = &topology.bonds[carrier.index()];
            if matches!(
                effective_direction(topology, updates, carrier),
                BondDirection::BeginWedge | BondDirection::BeginDash
            ) && updates
                .get(&carrier)
                .map_or(bond.begin(), |update| update.begin)
                == end.atom
                && if matches!(mode, WedgeMode::ThreeD(_)) {
                    can_have_direction(bond)
                } else {
                    // The no-conformer and 2D source paths call
                    // canHaveDirection() on the axial bond at this point.
                    can_have_direction(axial)
                }
            {
                existing.push((end_index, carrier_index, carrier));
            }
        }
    }
    if !existing.is_empty() {
        for (end_index, carrier_index, carrier) in existing {
            let direction = desired_direction(
                topology,
                axial,
                mode,
                end_index,
                carrier_index,
                carrier,
                ends[end_index].atom,
            );
            record_update(
                topology,
                updates,
                carrier,
                ends[end_index].atom,
                direction,
                axial.id(),
            );
        }
        return Ok(());
    }

    #[derive(Clone, Copy)]
    struct Best {
        end: usize,
        carrier_index: usize,
        bond: BondId,
        ring_count: usize,
        ring_size: usize,
        single: bool,
        direction: BondDirection,
    }
    let mut best: Option<Best> = None;
    let mut largest_ring_size = 0;
    for (end_index, end) in ends.iter().enumerate() {
        for (carrier_index, carrier) in end.bonds.iter().copied().enumerate() {
            let candidate = &topology.bonds[carrier.index()];
            let source_occupied = if matches!(mode, WedgeMode::ThreeD(_)) {
                occupied.contains(&axial.id()) || updates.contains_key(&axial.id())
            } else {
                occupied.contains(&carrier) || updates.contains_key(&carrier)
            };
            if !can_have_direction(candidate) || source_occupied {
                continue;
            }
            let direction_now = effective_direction(topology, updates, carrier);
            if direction_now != BondDirection::None {
                let begin = updates
                    .get(&carrier)
                    .map_or(candidate.begin(), |value| value.begin);
                match mode {
                    WedgeMode::NoConformer if begin == end.atom => {
                        return Err(AtropisomerRejectionKind::DirectionConflict);
                    }
                    WedgeMode::TwoD { .. }
                        if begin == end.atom
                            && matches!(
                                direction_now,
                                BondDirection::BeginWedge | BondDirection::BeginDash
                            ) =>
                    {
                        return Err(AtropisomerRejectionKind::DirectionConflict);
                    }
                    WedgeMode::ThreeD(_) => {
                        // The source reorients the candidate to the axial end
                        // before this check, so any remaining direction is a
                        // conflict in this branch.
                        return Err(AtropisomerRejectionKind::DirectionConflict);
                    }
                    _ => {}
                }
                continue;
            }
            let source_ring_count = rings.num_bond_rings(carrier);
            let (ring_count, ring_size) = if matches!(mode, WedgeMode::NoConformer) {
                (source_ring_count, 0)
            } else if source_ring_count == 0 {
                (10, 0)
            } else {
                let size = rings.min_bond_ring_size(carrier);
                (source_ring_count, if size > 8 { 0 } else { size })
            };
            let direction = desired_direction(
                topology,
                axial,
                mode,
                end_index,
                carrier_index,
                carrier,
                end.atom,
            );
            let row = Best {
                end: end_index,
                carrier_index,
                bond: carrier,
                ring_count,
                ring_size,
                single: candidate.order() == BondOrder::Single,
                direction,
            };
            let (replace, update_largest_ring_size) = match best {
                None => (true, true),
                Some(current) if row.ring_count < current.ring_count => (true, true),
                Some(current)
                    if row.ring_count == current.ring_count
                        && row.ring_size > largest_ring_size =>
                {
                    (true, true)
                }
                Some(current)
                    if row.ring_count == current.ring_count && row.single != current.single =>
                {
                    (row.single, false)
                }
                Some(current)
                    if row.ring_count == current.ring_count && row.single == current.single =>
                {
                    (
                        current.direction == BondDirection::BeginDash
                            && row.direction == BondDirection::BeginWedge,
                        false,
                    )
                }
                _ => (false, false),
            };
            if replace {
                if update_largest_ring_size {
                    largest_ring_size = row.ring_size;
                }
                best = Some(row);
            }
        }
    }
    // Complete pinned source: WedgeBondFromAtropisomerOneBondNoConf/2d/3d.
    // RDKit✔️✔️: bool WedgeBondFromAtropisomerOneBondNoConf(
    // RDKit✔️✔️:     Bond *bond, const ROMol &mol,
    // RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit✔️✔️:         &wedgeBonds) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit✔️✔️:     return false;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   //  make sure we do not have wiggle bonds
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atomAndBondVec : atomAndBondVecs) {
    // RDKit✔️✔️:     for (auto endBond : atomAndBondVec.second) {
    // RDKit✔️✔️:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         return false;  // not an atropisomer)
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // first see if any candidate bond is already set to a wedge or hash
    // RDKit✔️✔️:   // if so, we will use that bond as a wedge or hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<int> useBondsAtEnd[2];
    // RDKit✔️✔️:   bool foundBondDir = false;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:     for (unsigned int whichBond = 0;
    // RDKit✔️✔️:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit✔️✔️:       auto bondDir = atomAndBondVecs[whichEnd].second[whichBond]->getBondDir();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // see if it is a wedge or hash and its origin is the atom in the
    // RDKit✔️✔️:       // main bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
    // RDKit✔️✔️:           atomAndBondVecs[whichEnd].second[whichBond]->getBeginAtom() ==
    // RDKit✔️✔️:               atomAndBondVecs[whichEnd].first &&
    // RDKit✔️✔️:           canHaveDirection(*bond)) {
    // RDKit✔️✔️:         useBondsAtEnd[whichEnd].push_back(whichBond);
    // RDKit✔️✔️:         foundBondDir = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (foundBondDir) {
    // RDKit✔️✔️:     for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:       for (unsigned int whichBondIndex = 0;
    // RDKit✔️✔️:            whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
    // RDKit✔️✔️:         atomAndBondVecs[whichEnd]
    // RDKit✔️✔️:             .second[useBondsAtEnd[whichEnd][whichBondIndex]]
    // RDKit✔️✔️:             ->setBondDir(getBondDirForAtropisomerNoConf(
    // RDKit✔️✔️:                 bond->getStereo(), whichEnd,
    // RDKit✔️✔️:                 useBondsAtEnd[whichEnd][whichBondIndex]));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // did not find a good bond dir - pick one to use
    // RDKit✔️✔️:   // we would like to have one that is not in a ring, and will be a wedge
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int bestBondEnd = -1, bestBondNumber = -1;
    // RDKit✔️✔️:   bool bestBondIsSingle = false;
    // RDKit✔️✔️:   unsigned int bestRingCount = INT_MAX;
    // RDKit✔️✔️:   Bond::BondDir bestBondDir = Bond::BondDir::NONE;
    // RDKit✔️✔️:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:     for (unsigned int whichBond = 0;
    // RDKit✔️✔️:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit✔️✔️:       auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!canHaveDirection(*bondToTry) ||
    // RDKit✔️✔️:           wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
    // RDKit✔️✔️:         continue;  // must be a single OR aromatic bond and not already
    // RDKit✔️✔️:                    // spoken for by a chiral center
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
    // RDKit✔️✔️:         if (bondToTry->getBeginAtom()->getIdx() ==
    // RDKit✔️✔️:             atomAndBondVecs[whichEnd].first->getIdx()) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
    // RDKit✔️✔️:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           continue;  // wedge or hash bond affecting the OTHER atom
    // RDKit✔️✔️:                      // = perhaps a chiral center
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto ringCount = ri->numBondRings(bondToTry->getIdx());
    // RDKit✔️✔️:       if (ringCount > bestRingCount) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       else if (ringCount < bestRingCount) {
    // RDKit✔️✔️:         bestBondEnd = whichEnd;
    // RDKit✔️✔️:         bestBondNumber = whichBond;
    // RDKit✔️✔️:         bestRingCount = ringCount;
    // RDKit✔️✔️:         bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit✔️✔️:         bestBondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
    // RDKit✔️✔️:                                                      whichEnd, whichBond);
    // RDKit✔️✔️:       } else if (bestBondIsSingle &&
    // RDKit✔️✔️:                  bondToTry->getBondType() != Bond::BondType::SINGLE) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       } else if (!bestBondIsSingle &&
    // RDKit✔️✔️:                  bondToTry->getBondType() == Bond::BondType::SINGLE) {
    // RDKit✔️✔️:         bestBondEnd = whichEnd;
    // RDKit✔️✔️:         bestBondNumber = whichBond;
    // RDKit✔️✔️:         bestRingCount = ringCount;
    // RDKit✔️✔️:         bestBondIsSingle = true;
    // RDKit✔️✔️:         bestBondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
    // RDKit✔️✔️:                                                      whichEnd, whichBond);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto bondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
    // RDKit✔️✔️:                                                       whichEnd, whichBond);
    // RDKit✔️✔️:         if (bestBondDir == Bond::BondDir::NONE ||
    // RDKit✔️✔️:             (bestBondDir == Bond::BondDir::BEGINDASH &&
    // RDKit✔️✔️:              bondDir == Bond::BondDir::BEGINWEDGE)) {
    // RDKit✔️✔️:           bestBondEnd = whichEnd;
    // RDKit✔️✔️:           bestBondNumber = whichBond;
    // RDKit✔️✔️:           bestRingCount = ringCount;
    // RDKit✔️✔️:           bestBondIsSingle =
    // RDKit✔️✔️:               (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit✔️✔️:           bestBondDir = bondDir;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bestBondEnd >= 0)  // we found a good one
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     // make sure the atoms on the bond are in the right order for the
    // RDKit✔️✔️:     // wedge/hash the atom on the end of the main bond must be listed
    // RDKit✔️✔️:     // first for the wedge/has bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
    // RDKit✔️✔️:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
    // RDKit✔️✔️:       bestBond->setEndAtom(bestBond->getBeginAtom());
    // RDKit✔️✔️:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     bestBond->setBondDir(bestBondDir);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
    // RDKit✔️✔️:         new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
    // RDKit✔️✔️:                                                    bestBondDir));
    // RDKit✔️✔️:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: bool WedgeBondFromAtropisomerOneBond2d(
    // RDKit✔️✔️:     Bond *bond, const ROMol &mol, const Conformer *conf,
    // RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit✔️✔️:         &wedgeBonds) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit✔️✔️:     return false;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   //  make sure we do not have wiggle bonds
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atomAndBondVec : atomAndBondVecs) {
    // RDKit✔️✔️:     for (auto endBond : atomAndBondVec.second) {
    // RDKit✔️✔️:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         return false;  // not an atropisomer)
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // create a frame of reference that has its X-axis along the atrop bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D xAxis, yAxis, zAxis;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // RDKit✔️✔️:     // connot percieve atroisomer bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Cound not get a frame of reference for an atropisomer bond - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
    // RDKit✔️✔️:                                 // potential atropisome bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit✔️✔️:     // find a vector to represent the lowest numbered atom on each end
    // RDKit✔️✔️:     // this vector is NOT the bond vector, but is y-value in the bond
    // RDKit✔️✔️:     // frame or reference
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!getAtropIsomerEndVect(atomAndBondVecs[bondAtomIndex], yAxis, zAxis,
    // RDKit✔️✔️:                                conf, bondVecs[bondAtomIndex])) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:       // did not find a non-colinear bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Failed to get a representative vector for the defining bond of an atropisomer - atoms are: "
    // RDKit✔️✔️:           << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // first see if any candidate bond is already set to a wedge or hash
    // RDKit✔️✔️:   // if so, we will use that bond as a wedge or hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<int> useBondsAtEnd[2];
    // RDKit✔️✔️:   bool foundBondDir = false;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:     for (unsigned int whichBond = 0;
    // RDKit✔️✔️:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit✔️✔️:       auto bondDir = atomAndBondVecs[whichEnd].second[whichBond]->getBondDir();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // see if it is a wedge or hash and its origin is the atom in the
    // RDKit✔️✔️:       // main bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
    // RDKit✔️✔️:           atomAndBondVecs[whichEnd].second[whichBond]->getBeginAtom() ==
    // RDKit✔️✔️:               atomAndBondVecs[whichEnd].first &&
    // RDKit✔️✔️:           canHaveDirection(*bond)) {
    // RDKit✔️✔️:         useBondsAtEnd[whichEnd].push_back(whichBond);
    // RDKit✔️✔️:         foundBondDir = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (foundBondDir) {
    // RDKit✔️✔️:     for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:       for (unsigned int whichBondIndex = 0;
    // RDKit✔️✔️:            whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
    // RDKit✔️✔️:         atomAndBondVecs[whichEnd]
    // RDKit✔️✔️:             .second[useBondsAtEnd[whichEnd][whichBondIndex]]
    // RDKit✔️✔️:             ->setBondDir(getBondDirForAtropisomer2d(
    // RDKit✔️✔️:                 bondVecs, bond->getStereo(), whichEnd,
    // RDKit✔️✔️:                 useBondsAtEnd[whichEnd][whichBondIndex]));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // did not find a good bond dir - pick one to use
    // RDKit✔️✔️:   // we would like to have one that is in a ring, and will favor it being a
    // RDKit✔️✔️:   // wedge
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // We favor rings here because wedging non-ring bonds makes it too likely that
    // RDKit✔️✔️:   // we'll end up accidentally creating new atropisomeric bonds. This was github
    // RDKit✔️✔️:   // issue 7371
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int bestBondEnd = -1, bestBondNumber = -1;
    // RDKit✔️✔️:   bool bestBondIsSingle = false;
    // RDKit✔️✔️:   unsigned int bestRingCount = INT_MAX;
    // RDKit✔️✔️:   unsigned int largestRingSize = 0;
    // RDKit✔️✔️:   Bond::BondDir bestBondDir = Bond::BondDir::NONE;
    // RDKit✔️✔️:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:     for (unsigned int whichBond = 0;
    // RDKit✔️✔️:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit✔️✔️:       auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!canHaveDirection(*bondToTry) ||
    // RDKit✔️✔️:           wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
    // RDKit✔️✔️:         continue;  // must be a single OR aromatic bond and not already
    // RDKit✔️✔️:                    // spoken for by a chiral center
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
    // RDKit✔️✔️:         if (bondToTry->getBeginAtom()->getIdx() ==
    // RDKit✔️✔️:             atomAndBondVecs[whichEnd].first->getIdx()) {
    // RDKit✔️✔️:           if (bondToTry->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit✔️✔️:               bondToTry->getBondDir() == Bond::BEGINDASH) {
    // RDKit✔️✔️:             BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:                 << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
    // RDKit✔️✔️:                 << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:                 << std::endl;
    // RDKit✔️✔️:             return false;
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             continue;  // probably a slash up or down for a double bond
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           continue;  // wedge or hash bond affecting the OTHER atom
    // RDKit✔️✔️:                      // = perhaps a chiral center
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto ringCount = ri->numBondRings(bondToTry->getIdx());
    // RDKit✔️✔️:       unsigned int ringSize = 0;
    // RDKit✔️✔️:       if (!ringCount) {
    // RDKit✔️✔️:         ringCount = 10;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // we're going to prefer to put wedges in larger rings, but don't want
    // RDKit✔️✔️:         // to end up wedging macrocyles if it's avoidable.
    // RDKit✔️✔️:         ringSize = ri->minBondRingSize(bondToTry->getIdx());
    // RDKit✔️✔️:         if (ringSize > 8) {
    // RDKit✔️✔️:           ringSize = 0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (ringCount > bestRingCount) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       } else if (ringCount < bestRingCount || ringSize > largestRingSize) {
    // RDKit✔️✔️:         bestBondEnd = whichEnd;
    // RDKit✔️✔️:         bestBondNumber = whichBond;
    // RDKit✔️✔️:         bestRingCount = ringCount;
    // RDKit✔️✔️:         largestRingSize = ringSize;
    // RDKit✔️✔️:         bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit✔️✔️:         bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
    // RDKit✔️✔️:                                                  whichEnd, whichBond);
    // RDKit✔️✔️:       } else if (bestBondIsSingle &&
    // RDKit✔️✔️:                  bondToTry->getBondType() != Bond::BondType::SINGLE) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       } else if (!bestBondIsSingle &&
    // RDKit✔️✔️:                  bondToTry->getBondType() == Bond::BondType::SINGLE) {
    // RDKit✔️✔️:         bestBondEnd = whichEnd;
    // RDKit✔️✔️:         bestBondNumber = whichBond;
    // RDKit✔️✔️:         bestRingCount = ringCount;
    // RDKit✔️✔️:         bestBondIsSingle = true;
    // RDKit✔️✔️:         bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
    // RDKit✔️✔️:                                                  whichEnd, whichBond);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto bondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
    // RDKit✔️✔️:                                                   whichEnd, whichBond);
    // RDKit✔️✔️:         if (bestBondDir == Bond::BondDir::NONE ||
    // RDKit✔️✔️:             (bestBondDir == Bond::BondDir::BEGINDASH &&
    // RDKit✔️✔️:              bondDir == Bond::BondDir::BEGINWEDGE)) {
    // RDKit✔️✔️:           bestBondEnd = whichEnd;
    // RDKit✔️✔️:           bestBondNumber = whichBond;
    // RDKit✔️✔️:           bestRingCount = ringCount;
    // RDKit✔️✔️:           bestBondIsSingle =
    // RDKit✔️✔️:               (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit✔️✔️:           bestBondDir = bondDir;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bestBondEnd >= 0) {
    // RDKit✔️✔️:     // we found a good one
    // RDKit✔️✔️:     // make sure the atoms on the bond are in the right order for the
    // RDKit✔️✔️:     // wedge/hash the atom on the end of the main bond must be listed
    // RDKit✔️✔️:     // first for the wedge/has bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
    // RDKit✔️✔️:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
    // RDKit✔️✔️:       bestBond->setEndAtom(bestBond->getBeginAtom());
    // RDKit✔️✔️:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bestBond->setBondDir(bestBondDir);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
    // RDKit✔️✔️:         new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
    // RDKit✔️✔️:                                                    bestBondDir));
    // RDKit✔️✔️:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: bool WedgeBondFromAtropisomerOneBond3d(
    // RDKit✔️✔️:     Bond *bond, const ROMol &mol, const Conformer *conf,
    // RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit✔️✔️:         &wedgeBonds) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit✔️✔️:     return false;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   //  make sure we do not have wiggle bonds
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atomAndBondVecs : atomAndBondVecs) {
    // RDKit✔️✔️:     for (auto endBond : atomAndBondVecs.second) {
    // RDKit✔️✔️:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         return false;  // not an atropisomer)
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // first see if any candidate bond is already set to a wedge or hash
    // RDKit✔️✔️:   // if so, we will use that bond as a wedge or hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Bond *> useBonds;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:     for (unsigned int whichBond = 0;
    // RDKit✔️✔️:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit✔️✔️:       auto bond = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit✔️✔️:       auto bondDir = bond->getBondDir();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // see if it is a wedge or hash and its origin is the atom in the
    // RDKit✔️✔️:       // main bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
    // RDKit✔️✔️:           bond->getBeginAtom() == atomAndBondVecs[whichEnd].first &&
    // RDKit✔️✔️:           canHaveDirection(*bond)) {
    // RDKit✔️✔️:         useBonds.push_back(bond);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // the following may seem redundant, since we just found the useBonds
    // RDKit✔️✔️:   // based on their bond dir PRESENCE, but this endures that the values are
    // RDKit✔️✔️:   // correct.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (useBonds.size() > 0) {
    // RDKit✔️✔️:     for (auto useBond : useBonds) {
    // RDKit✔️✔️:       useBond->setBondDir(getBondDirForAtropisomer3d(useBond, conf));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // did not find a used bond dir - pick one to use
    // RDKit✔️✔️:   // we would like to have one that is not in a ring, and will be a dash
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Bond *bestBond = nullptr;
    // RDKit✔️✔️:   int bestBondEnd = -1;
    // RDKit✔️✔️:   unsigned int bestRingCount = UINT_MAX;
    // RDKit✔️✔️:   unsigned int largestRingSize = 0;
    // RDKit✔️✔️:   Bond::BondDir bestBondDir = Bond::BondDir::NONE;
    // RDKit✔️✔️:   bool bestBondIsSingle = false;
    // RDKit✔️✔️:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit✔️✔️:     for (unsigned int whichBond = 0;
    // RDKit✔️✔️:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit✔️✔️:       auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // cannot use a bond that is not single, nor if it is already slated
    // RDKit✔️✔️:       // to be used for a chiral center
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!canHaveDirection(*bondToTry) ||
    // RDKit✔️✔️:           wedgeBonds.find(bond->getIdx()) != wedgeBonds.end()) {
    // RDKit✔️✔️:         continue;  // must be a single bond and not already spoken
    // RDKit✔️✔️:                    // for by a chiral center
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // make sure the atoms on the bond are in the right order for the
    // RDKit✔️✔️:       // wedge/hash the atom on the end of the main bond must be listed
    // RDKit✔️✔️:       // first
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondToTry->getBeginAtom() != atomAndBondVecs[whichEnd].first) {
    // RDKit✔️✔️:         bondToTry->setEndAtom(bondToTry->getBeginAtom());
    // RDKit✔️✔️:         bondToTry->setBeginAtom(atomAndBondVecs[whichEnd].first);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
    // RDKit✔️✔️:         if (bondToTry->getBeginAtom()->getIdx() ==
    // RDKit✔️✔️:             atomAndBondVecs[whichEnd].first->getIdx()) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
    // RDKit✔️✔️:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           continue;  // wedge or hash bond affecting the OTHER atom
    // RDKit✔️✔️:                      // = perhaps a chiral center
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto ringCount = ri->numBondRings(bondToTry->getIdx());
    // RDKit✔️✔️:       unsigned int ringSize = 0;
    // RDKit✔️✔️:       if (!ringCount) {
    // RDKit✔️✔️:         ringCount = 10;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // we're going to prefer to put wedges in larger rings, but don't want
    // RDKit✔️✔️:         // to end up wedging macrocyles if it's avoidable.
    // RDKit✔️✔️:         ringSize = ri->minBondRingSize(bondToTry->getIdx());
    // RDKit✔️✔️:         if (ringSize > 8) {
    // RDKit✔️✔️:           ringSize = 0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (ringCount > bestRingCount) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       } else if (ringCount < bestRingCount || ringSize > largestRingSize) {
    // RDKit✔️✔️:         bestBond = bondToTry;
    // RDKit✔️✔️:         bestBondEnd = whichEnd;
    // RDKit✔️✔️:         bestRingCount = ringCount;
    // RDKit✔️✔️:         largestRingSize = ringSize;
    // RDKit✔️✔️:         bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit✔️✔️:         bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
    // RDKit✔️✔️:       } else if (bestBondIsSingle &&
    // RDKit✔️✔️:                  bondToTry->getBondType() != Bond::BondType::SINGLE) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       } else if (!bestBondIsSingle &&
    // RDKit✔️✔️:                  bondToTry->getBondType() == Bond::BondType::SINGLE) {
    // RDKit✔️✔️:         bestBondEnd = whichEnd;
    // RDKit✔️✔️:         bestBond = bondToTry;
    // RDKit✔️✔️:         bestRingCount = ringCount;
    // RDKit✔️✔️:         bestBondIsSingle = true;
    // RDKit✔️✔️:         bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto bondDir = getBondDirForAtropisomer3d(bondToTry, conf);
    // RDKit✔️✔️:         if (bestBondDir == Bond::BondDir::NONE ||
    // RDKit✔️✔️:             (bestBondDir == Bond::BondDir::BEGINDASH &&
    // RDKit✔️✔️:              bondDir == Bond::BondDir::BEGINWEDGE)) {
    // RDKit✔️✔️:           bestBond = bondToTry;
    // RDKit✔️✔️:           bestBondEnd = whichEnd;
    // RDKit✔️✔️:           bestRingCount = ringCount;
    // RDKit✔️✔️:           bestBondIsSingle =
    // RDKit✔️✔️:               (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:           bestBondDir = bondDir;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bestBond != nullptr) {
    // RDKit✔️✔️:     // we found a good one
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // make sure the atoms on the bond are in the right order for the
    // RDKit✔️✔️:     // wedge/hash the atom on the end of the main bond must be listed
    // RDKit✔️✔️:     // first for the wedge/has bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
    // RDKit✔️✔️:       bestBond->setEndAtom(bestBond->getBeginAtom());
    // RDKit✔️✔️:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bestBond->setBondDir(bestBondDir);
    // RDKit✔️✔️:     auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
    // RDKit✔️✔️:         new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
    // RDKit✔️✔️:                                                    bestBondDir));
    // RDKit✔️✔️:
    // RDKit✔️✔️:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    let best = best.ok_or(AtropisomerRejectionKind::NoUsableWedgeBond)?;
    debug_assert_eq!(best.bond, ends[best.end].bonds[best.carrier_index]);
    record_update(
        topology,
        updates,
        best.bond,
        ends[best.end].atom,
        best.direction,
        axial.id(),
    );
    Ok(())
}

fn validate_rings(topology: &TopologyBlock, rings: &RingInfo) -> Result<(), AtropisomerError> {
    if !rings.is_sssr_or_better() {
        return Err(AtropisomerError::RingInfoNotSssr);
    }
    for ring in rings.atom_rings() {
        for atom in ring {
            if atom.index() >= topology.atoms.len() {
                return Err(AtropisomerError::RingAtomOutOfRange {
                    atom: *atom,
                    atom_count: topology.atoms.len(),
                });
            }
        }
    }
    for ring in rings.bond_rings() {
        for bond in ring {
            if bond.index() >= topology.bonds.len() {
                return Err(AtropisomerError::RingBondOutOfRange {
                    bond: *bond,
                    bond_count: topology.bonds.len(),
                });
            }
        }
    }
    Ok(())
}

pub fn wedge_bonds_from_atropisomers(
    topology: &TopologyBlock,
    rings: &RingInfo,
    conformer: Option<AtropisomerConformer<'_>>,
    occupied_bonds: &BTreeSet<BondId>,
) -> Result<AtropisomerWedgeAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    validate_rings(topology, rings)?;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    for bond in occupied_bonds {
        if bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::AssignmentBondOutOfRange {
                bond: *bond,
                bond_count: topology.bonds.len(),
            });
        }
    }
    // Complete pinned source: wedgeBondsFromAtropisomers.
    // RDKit✔️✔️: void wedgeBondsFromAtropisomers(
    // RDKit✔️✔️:     const ROMol &mol, const Conformer *conf,
    // RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit✔️✔️:         &wedgeBonds) {
    // RDKit✔️✔️:   PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
    // RDKit✔️✔️:                "conformer does not belong to molecule");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // WedgeBondFromAtropisomerOneBond 2d/3d requires ring bond counts
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:     RDKit::MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto bondStereo = bond->getStereo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bond->getBondType() != Bond::BondType::SINGLE ||
    // RDKit✔️✔️:         (bondStereo != Bond::BondStereo::STEREOATROPCW &&
    // RDKit✔️✔️:          bondStereo != Bond::BondStereo::STEREOATROPCCW) ||
    // RDKit✔️✔️:         bond->getBeginAtom()->getTotalDegree() < 2 ||
    // RDKit✔️✔️:         bond->getEndAtom()->getTotalDegree() < 2 ||
    // RDKit✔️✔️:         bond->getBeginAtom()->getTotalDegree() > 3 ||
    // RDKit✔️✔️:         bond->getEndAtom()->getTotalDegree() > 3) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (conf) {
    // RDKit✔️✔️:       if (conf->is3D()) {
    // RDKit✔️✔️:         WedgeBondFromAtropisomerOneBond3d(bond, mol, conf, wedgeBonds);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         WedgeBondFromAtropisomerOneBond2d(bond, mol, conf, wedgeBonds);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {  // no conformer
    // RDKit✔️✔️:       WedgeBondFromAtropisomerOneBondNoConf(bond, mol, wedgeBonds);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut updates = BTreeMap::new();
    let mut diagnostics = Vec::new();
    for axial in &topology.bonds {
        if axial.order() != BondOrder::Single
            || !matches!(axial.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw)
            || !(2..=3).contains(&total_degree(topology, axial.begin()))
            || !(2..=3).contains(&total_degree(topology, axial.end()))
        {
            continue;
        }
        if let Err(kind) = wedge_one(
            topology,
            rings,
            axial,
            conformer,
            occupied_bonds,
            &mut updates,
        ) {
            diagnostics.push(AtropisomerDiagnostic {
                bond: axial.id(),
                kind,
            });
        }
    }
    Ok(AtropisomerWedgeAssignment {
        bond_updates: updates.into_values().collect(),
        diagnostics,
    })
}

#[must_use]
pub fn does_topology_have_atropisomers(topology: &TopologyBlock) -> bool {
    // Complete pinned source: doesMolHaveAtropisomers.
    // RDKit✔️✔️: bool doesMolHaveAtropisomers(const ROMol &mol) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto bondStereo = bond->getStereo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bondStereo == Bond::BondStereo::STEREOATROPCW ||
    // RDKit✔️✔️:         bondStereo == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    topology
        .bonds
        .iter()
        .any(|bond| matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
}
