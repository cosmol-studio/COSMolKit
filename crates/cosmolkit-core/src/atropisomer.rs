//! Detached, source-backed atropisomer perception and wedge assignment.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, Bond, BondId, BondValueError, Conformer2D, Conformer3D, CoordinateValidationError,
    StereoGroup, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Hybridization};

use crate::{HybridizationAssignment, RingInfo, WedgeAssignments, WedgeInfo};

const REALLY_SMALL_BOND_LEN: f64 = 0.000_000_1;

#[derive(Debug, Clone, Copy)]
pub enum AtropisomerConformer<'a> {
    TwoD(&'a Conformer2D),
    ThreeD(&'a Conformer3D),
}

fn conformer_is_3d(conformer: AtropisomerConformer<'_>) -> bool {
    // BEGIN RDKIT CPP FUNCTION getBondFrameOfReference
    // RDKit✔️✔️:   if (!conf->is3D()) {
    // RDKit✔️✔️:     yAxis = RDGeom::Point3D(-xAxis.y, xAxis.x, 0);
    // RDKit✔️✔️:     yAxis.normalize();
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(0.0, 0.0, 1.0);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Behavior review: a Conformer3D is a lossless XYZ storage carrier, not
    // proof that the source conformer flag is 3D. The independent is_3d bit
    // selects the same source branch while all XYZ bits remain available.
    // Complexity review: one enum match and one flag read are constant-time
    // and introduce no allocation or coordinate projection.
    match conformer {
        AtropisomerConformer::TwoD(_) => false,
        AtropisomerConformer::ThreeD(value) => value.is_3d(),
    }
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtropisomerCarrierEnd {
    focus: AtomId,
    carrier_bonds: Vec<BondId>,
}

impl AtropisomerCarrierEnd {
    #[must_use]
    pub const fn focus(&self) -> AtomId {
        self.focus
    }

    #[must_use]
    pub fn carrier_bonds(&self) -> &[BondId] {
        &self.carrier_bonds
    }
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
    #[error("ring atom membership has {actual} rows; expected {expected}")]
    RingAtomRowCount { actual: usize, expected: usize },
    #[error("ring bond membership has {actual} rows; expected {expected}")]
    RingBondRowCount { actual: usize, expected: usize },
    #[error("hybridization assignment has {actual} rows; expected {expected}")]
    HybridizationAssignmentLength { actual: usize, expected: usize },
    #[error("failed to clear atropisomer stereo from bond {bond}: {source}")]
    BondUpdate {
        bond: BondId,
        source: BondValueError,
    },
    #[error("stereo group atom {atom} is out of range for {atom_count} atoms")]
    StereoGroupAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("stereo group bond {bond} is out of range for {bond_count} bonds")]
    StereoGroupBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("atropisomer assignment bond {bond} is out of range for {bond_count} bonds")]
    AssignmentBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("atropisomer axial bond {bond} is out of range for {bond_count} bonds")]
    AxialBondOutOfRange { bond: BondId, bond_count: usize },
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
        if end.bonds.is_empty() {
            return None;
        }
        if end.bonds.len() == 2
            && other_atom(&topology.bonds[end.bonds[1].index()], end.atom).index()
                < other_atom(&topology.bonds[end.bonds[0].index()], end.atom).index()
        {
            end.bonds.swap(0, 1);
        }
    }
    Some(result)
}

pub fn atropisomer_carriers(
    topology: &TopologyBlock,
    axial_bond: BondId,
) -> Result<Option<[AtropisomerCarrierEnd; 2]>, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    let bond =
        topology
            .bonds
            .get(axial_bond.index())
            .ok_or(AtropisomerError::AxialBondOutOfRange {
                bond: axial_bond,
                bond_count: topology.bonds.len(),
            })?;
    Ok(atropisomer_ends(topology, bond).map(|ends| {
        ends.map(|end| AtropisomerCarrierEnd {
            focus: end.atom,
            carrier_bonds: end.bonds,
        })
    }))
}

/// Returns carrier ends for requested axial bonds after validating topology once.
pub fn atropisomer_carriers_for_bonds(
    topology: &TopologyBlock,
    axial_bonds: &[BondId],
) -> Result<Vec<(BondId, Option<[AtropisomerCarrierEnd; 2]>)>, AtropisomerError> {
    // Validate once, then reuse the source-shaped `atropisomer_ends` owner
    // helper for each requested bond instead of repeating a whole-topology
    // scan for every atropisomer neighbor in the CX bond-config writer.
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    axial_bonds
        .iter()
        .copied()
        .map(|axial_bond| {
            let bond = topology.bonds.get(axial_bond.index()).ok_or(
                AtropisomerError::AxialBondOutOfRange {
                    bond: axial_bond,
                    bond_count: topology.bonds.len(),
                },
            )?;
            let carriers = atropisomer_ends(topology, bond).map(|ends| {
                ends.map(|end| AtropisomerCarrierEnd {
                    focus: end.atom,
                    carrier_bonds: end.bonds,
                })
            });
            Ok((axial_bond, carriers))
        })
        .collect()
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
    if !conformer_is_3d(conformer) {
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
        AtropisomerConformer::ThreeD(value) => value
            .validate_for_atom_count(atom_count)
            .map_err(|source| AtropisomerError::InvalidCoordinates { source }),
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
    if !conformer_is_3d(conformer) {
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
            // The pinned three-argument constructor propagates neither the
            // source read ID nor its write ID; this fresh value resets both.
            groups.push(StereoGroup::new(group.kind(), atoms, bonds));
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

pub fn get_all_atom_ids_for_stereo_group(
    topology: &TopologyBlock,
    group: &StereoGroup,
    wedge_bonds: &WedgeAssignments,
) -> Result<Vec<AtomId>, AtropisomerError> {
    validate_stereo_group_topology(topology)?;
    validate_stereo_group_members(topology, group)?;
    Ok(collect_stereo_group_atom_ids(topology, group, wedge_bonds))
}

/// Collects source-ordered atom IDs for each checked stereo group.
///
/// The topology is validated once for the batch. Group membership is checked
/// in input order before each source collection pass.
pub fn get_all_atom_ids_for_stereo_groups(
    topology: &TopologyBlock,
    groups: &[StereoGroup],
    wedge_bonds: &WedgeAssignments,
) -> Result<Vec<Vec<AtomId>>, AtropisomerError> {
    validate_stereo_group_topology(topology)?;

    let mut atom_ids_by_group = Vec::with_capacity(groups.len());
    for group in groups {
        validate_stereo_group_members(topology, group)?;
        atom_ids_by_group.push(collect_stereo_group_atom_ids(topology, group, wedge_bonds));
    }
    Ok(atom_ids_by_group)
}

#[cfg(test)]
std::thread_local! {
    static STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS: std::cell::Cell<usize> = const {
        std::cell::Cell::new(0)
    };
}

fn validate_stereo_group_topology(topology: &TopologyBlock) -> Result<(), AtropisomerError> {
    #[cfg(test)]
    STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS.with(|calls| calls.set(calls.get() + 1));

    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })
}

fn validate_stereo_group_members(
    topology: &TopologyBlock,
    group: &StereoGroup,
) -> Result<(), AtropisomerError> {
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
    Ok(())
}

fn collect_stereo_group_atom_ids(
    topology: &TopologyBlock,
    group: &StereoGroup,
    wedge_bonds: &WedgeAssignments,
) -> Vec<AtomId> {
    // Complete pinned source: Atropisomers::getAllAtomIdsForStereoGroup.
    // RDKit❗❌: void getAllAtomIdsForStereoGroup(
    // RDKit❗❌:     const ROMol &mol, const StereoGroup &group,
    // RDKit❗❌:     std::vector<unsigned int> &atomIds,
    // RDKit❗❌:     const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗❌:         &wedgeBonds) {
    // RDKit❗❌:   atomIds.clear();
    // RDKit❗❌:   for (auto &&atom : group.getAtoms()) {
    // RDKit❗❌:     atomIds.push_back(atom->getIdx());
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto &&bond : group.getBonds()) {
    // RDKit❗❌:     // figure out which atoms of the bond get wedge/hash indications
    // RDKit❗❌:     // mark the atom with the wedge/hash
    // RDKit❗❌:
    // RDKit❗❌:     for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
    // RDKit❗❌:       for (const auto atomBond : mol.atomBonds(atom)) {
    // RDKit❗❌:         if (atomBond->getIdx() == bond->getIdx()) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         if (atomBond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit❗❌:             atomBond->getBondDir() == Bond::BEGINDASH ||
    // RDKit❗❌:             (wedgeBonds.find(atomBond->getIdx()) != wedgeBonds.end() &&
    // RDKit❗❌:              (wedgeBonds.at(atomBond->getIdx())->getType()) ==
    // RDKit❗❌:                  Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer)) {
    // RDKit❗❌:           if (std::find(atomIds.begin(), atomIds.end(), atom->getIdx()) ==
    // RDKit❗❌:               atomIds.end()) {
    // RDKit❗❌:             atomIds.push_back(atom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Behavior review: the shared wedge map contributes only Atropisomer entries;
    // existing BEGINWEDGE/BEGINDASH directions, group member order and adjacency
    // source order are preserved, and newly found endpoint IDs are appended once.
    // Complexity review: adjacency and BTreeMap lookups avoid rescanning wedge
    // updates; linear Vec membership matches source std::find. Full detached
    // topology validation adds an O(V+E) pass absent from the pointer-valid source.
    let mut atom_ids = group.atoms().to_vec();
    for group_bond_id in group.bonds() {
        let group_bond = &topology.bonds[group_bond_id.index()];
        for atom in [group_bond.begin(), group_bond.end()] {
            for adjacent in topology.adjacency.neighbors_of(atom.index()) {
                if adjacent.bond == *group_bond_id {
                    continue;
                }

                let adjacent_bond = &topology.bonds[adjacent.bond.index()];
                if matches!(
                    adjacent_bond.direction(),
                    BondDirection::BeginWedge | BondDirection::BeginDash
                ) || matches!(
                    wedge_bonds.get(adjacent.bond),
                    Some(WedgeInfo::Atropisomer { .. })
                ) {
                    if !atom_ids.contains(&atom) {
                        atom_ids.push(atom);
                    }
                }
            }
        }
    }
    atom_ids
}

#[cfg(test)]
mod stereo_group_batch_tests {
    use super::*;
    use crate::{WedgeAssignments, pick_bonds_to_wedge};
    use cosmolkit_model::{Atom, AtomSpec, BondSpec, StereoGroupKind, TopologyValidationError};
    use cosmolkit_types::{BondOrder, ChiralTag, Element};

    fn group_topology(directions: [BondDirection; 3]) -> TopologyBlock {
        let atoms = (0..4)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::Unspecified),
                )
            })
            .collect();
        let edges = [(0, 1), (1, 2), (1, 3)];
        let bonds = edges
            .into_iter()
            .zip(directions)
            .enumerate()
            .map(|(index, ((begin, end), direction))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
                        .with_direction(direction),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("valid stereo-group test topology")
    }

    fn tetrahedral_wedge_topology() -> TopologyBlock {
        let atoms = (0..5)
            .map(|index| {
                let tag = if index == 0 {
                    ChiralTag::TetrahedralCw
                } else {
                    ChiralTag::Unspecified
                };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_chiral_tag(tag),
                )
            })
            .collect();
        let bonds = (0..4)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(0), AtomId::new(index + 1), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("valid tetrahedral wedge topology")
    }

    fn reset_validation_count() {
        STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS.with(|calls| calls.set(0));
    }

    fn validation_count() -> usize {
        STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS.with(std::cell::Cell::get)
    }

    #[test]
    fn cf_smi_groups_batch_preserves_member_order_duplicates_and_overlaps() {
        let topology = group_topology([
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
        ]);
        let groups = vec![
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(2), AtomId::new(0), AtomId::new(2)],
                Vec::new(),
            ),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(1)],
                Vec::new(),
            ),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(3)],
                vec![BondId::new(0)],
            ),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(1)],
                vec![BondId::new(0)],
            ),
        ];

        let expected = vec![
            vec![AtomId::new(2), AtomId::new(0), AtomId::new(2)],
            vec![AtomId::new(0), AtomId::new(1)],
            vec![AtomId::new(3), AtomId::new(1)],
            vec![AtomId::new(1)],
        ];
        let batch =
            get_all_atom_ids_for_stereo_groups(&topology, &groups, &WedgeAssignments::default())
                .expect("collect all valid groups");
        let independent = groups
            .iter()
            .map(|group| {
                get_all_atom_ids_for_stereo_group(&topology, group, &WedgeAssignments::default())
                    .expect("collect one valid group")
            })
            .collect::<Vec<_>>();

        assert_eq!(batch, expected);
        assert_eq!(batch, independent);
    }

    #[test]
    fn cf_smi_groups_batch_handles_empty_and_single_group() {
        let topology = group_topology([BondDirection::None; 3]);
        let empty =
            get_all_atom_ids_for_stereo_groups(&topology, &[], &WedgeAssignments::default())
                .expect("empty checked batch");
        assert!(empty.is_empty());

        let group = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(3), AtomId::new(0)],
            Vec::new(),
        );
        assert_eq!(
            get_all_atom_ids_for_stereo_groups(
                &topology,
                std::slice::from_ref(&group),
                &WedgeAssignments::default(),
            )
            .expect("one checked group"),
            vec![vec![AtomId::new(3), AtomId::new(0)]]
        );
    }

    #[test]
    fn cf_smi_groups_wedge_map_accepts_only_atropisomer_entries() {
        let topology = group_topology([BondDirection::None; 3]);
        let axial_group =
            StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), vec![BondId::new(0)]);
        let update = AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(1),
            end: AtomId::new(2),
            direction: BondDirection::None,
            atropisomer_bond: BondId::new(0),
        };
        let atrop_map =
            WedgeAssignments::from_atropisomer_wedge_assignment(AtropisomerWedgeAssignment {
                bond_updates: vec![update],
                diagnostics: Vec::new(),
            });
        assert_eq!(
            get_all_atom_ids_for_stereo_group(&topology, &axial_group, &atrop_map)
                .expect("Atropisomer map entry expands its endpoint"),
            vec![AtomId::new(1)]
        );

        let tetrahedral = tetrahedral_wedge_topology();
        let chiral_map = pick_bonds_to_wedge(&tetrahedral, None)
            .expect("valid tetrahedral topology receives a Chiral wedge map entry");
        let (chiral_bond, entry) = chiral_map
            .iter()
            .find(|(_, entry)| matches!(entry, WedgeInfo::Chiral { .. }))
            .expect("a tetrahedral center receives a Chiral map entry");
        assert!(matches!(entry, WedgeInfo::Chiral { .. }));
        let group_bond = (0..4)
            .map(BondId::new)
            .find(|bond| *bond != chiral_bond)
            .expect("a second central bond remains for the group");
        let chiral_group =
            StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), vec![group_bond]);
        assert_eq!(
            get_all_atom_ids_for_stereo_group(&tetrahedral, &chiral_group, &chiral_map)
                .expect("Chiral map entries do not expand an enhanced group"),
            Vec::<AtomId>::new()
        );
    }

    #[test]
    fn cf_smi_groups_batch_preserves_typed_errors_and_precedence() {
        let topology = group_topology([BondDirection::None; 3]);
        let invalid_members = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(99)],
            vec![BondId::new(99)],
        );
        assert_eq!(
            get_all_atom_ids_for_stereo_groups(
                &topology,
                std::slice::from_ref(&invalid_members),
                &WedgeAssignments::default(),
            ),
            Err(AtropisomerError::StereoGroupAtomOutOfRange {
                atom: AtomId::new(99),
                atom_count: 4,
            })
        );

        let invalid_bond = StereoGroup::new(StereoGroupKind::Or, Vec::new(), vec![BondId::new(99)]);
        assert_eq!(
            get_all_atom_ids_for_stereo_groups(
                &topology,
                &[
                    StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(0)], Vec::new(),),
                    invalid_bond
                ],
                &WedgeAssignments::default(),
            ),
            Err(AtropisomerError::StereoGroupBondOutOfRange {
                bond: BondId::new(99),
                bond_count: 3,
            })
        );

        let mut invalid_topology = group_topology([BondDirection::None; 3]);
        invalid_topology.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(99), AtomId::new(1), BondOrder::Single),
        );
        assert!(matches!(
            get_all_atom_ids_for_stereo_groups(
                &invalid_topology,
                &[StereoGroup::new(
                    StereoGroupKind::Absolute,
                    vec![AtomId::new(0)],
                    Vec::new(),
                )],
                &WedgeAssignments::default(),
            ),
            Err(AtropisomerError::InvalidTopology {
                source: TopologyValidationError::BondEndpointOutOfRange {
                    bond,
                    ..
                }
            })
                if bond == BondId::new(0)
        ));
    }

    #[test]
    fn cf_smi_groups_batch_validates_topology_once_per_call() {
        let topology = group_topology([BondDirection::None; 3]);
        let groups = (0..12)
            .map(|index| {
                StereoGroup::new(
                    StereoGroupKind::Or,
                    vec![AtomId::new(index % 4)],
                    Vec::new(),
                )
            })
            .collect::<Vec<_>>();

        reset_validation_count();
        get_all_atom_ids_for_stereo_groups(&topology, &groups, &WedgeAssignments::default())
            .expect("batch of twelve groups");
        assert_eq!(validation_count(), 1);

        reset_validation_count();
        get_all_atom_ids_for_stereo_group(&topology, &groups[0], &WedgeAssignments::default())
            .expect("checked singular group");
        assert_eq!(validation_count(), 1);

        reset_validation_count();
        get_all_atom_ids_for_stereo_groups(&topology, &[], &WedgeAssignments::default())
            .expect("empty checked batch still validates the topology");
        assert_eq!(validation_count(), 1);
    }
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
        Some(value) if !conformer_is_3d(value) => {
            let frame =
                frame_of_reference(axial, value).ok_or(AtropisomerRejectionKind::ZeroLengthAxis)?;
            WedgeMode::TwoD {
                vectors: [
                    end_vector(topology, &ends[0], frame, value)?,
                    end_vector(topology, &ends[1], frame, value)?,
                ],
            }
        }
        Some(value) => WedgeMode::ThreeD(value),
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
    if rings.atom_row_count() != topology.atoms.len() {
        return Err(AtropisomerError::RingAtomRowCount {
            actual: rings.atom_row_count(),
            expected: topology.atoms.len(),
        });
    }
    if rings.bond_row_count() != topology.bonds.len() {
        return Err(AtropisomerError::RingBondRowCount {
            actual: rings.bond_row_count(),
            expected: topology.bonds.len(),
        });
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

fn check_invalid_atrop_bond(
    bond: &mut Bond,
    hybridization: &HybridizationAssignment,
    rings: &RingInfo,
) -> Result<bool, AtropisomerError> {
    // Complete pinned source: MolOps.cpp anonymous-namespace checkBond.
    // RDKit✔️✔️: bool checkBond(RWMol &mol, Bond *bond, MolOps::Hybridizations &hybs) {
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:     RDKit::MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const RingInfo *ri = mol.getRingInfo();
    // RDKit✔️✔️:   if (hybs[bond->getBeginAtomIdx()] != Atom::SP2 ||
    // RDKit✔️✔️:       hybs[bond->getEndAtomIdx()] != Atom::SP2 ||
    // RDKit✔️✔️:       // do not clear bonds that part of a macrocycle
    // RDKit✔️✔️:       // because they can be linking actual atropisomeric portions
    // RDKit✔️✔️:       (ri->numBondRings(bond->getIdx()) > 0 &&
    // RDKit✔️✔️:        ri->minBondRingSize(bond->getIdx()) < 8)) {
    // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Ring discovery is performed by the canonical ring owner before this
    // detached primitive is called; `validate_rings()` enforces the same
    // SSSR-or-better postcondition before any bond can be changed.
    let invalid = hybridization.values[bond.begin().index()] != Hybridization::Sp2
        || hybridization.values[bond.end().index()] != Hybridization::Sp2
        || (rings.num_bond_rings(bond.id()) > 0 && rings.min_bond_ring_size(bond.id()) < 8);
    if invalid {
        let bond_id = bond.id();
        bond.set_stereo(BondStereo::None)
            .map_err(|source| AtropisomerError::BondUpdate {
                bond: bond_id,
                source,
            })?;
        return Ok(true);
    }
    Ok(false)
}

/// Clear source-invalid atropisomer bond stereo over detached topology values.
///
/// The caller supplies canonical hybridization and ring assignments from the
/// immediately preceding sanitize stages. The borrowed topology is never
/// modified; failure returns no partial result.
pub fn cleanup_invalid_atropisomers(
    topology: &TopologyBlock,
    hybridization: &HybridizationAssignment,
    rings: &RingInfo,
) -> Result<TopologyBlock, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    validate_rings(topology, rings)?;
    if hybridization.values.len() != topology.atoms.len() {
        return Err(AtropisomerError::HybridizationAssignmentLength {
            actual: hybridization.values.len(),
            expected: topology.atoms.len(),
        });
    }

    // The one private engine below owns the clone/tag-loop/group-cleanup/
    // final-validation body; this wrapper keeps the source wrapper-before-
    // loop prevalidation and supplies the existing detached check primitive.
    cleanup_invalid_atropisomers_engine(topology, |bond| {
        check_invalid_atrop_bond(bond, hybridization, rings)
    })
}

/// One private error-generic cleanup engine shared by the public detached
/// wrapper and the owned-ring-state sanitize wrapper.
///
/// `check` is called ONLY for AtropCw/AtropCcw bonds; it reports whether it
/// cleared the tag. Stereo-group cleanup inspects the RESULT topology after
/// tag updates and runs only if a tag was cleared, exactly as the source
/// does after its bond loop.
fn cleanup_invalid_atropisomers_engine<E>(
    topology: &TopologyBlock,
    mut check: impl FnMut(&mut Bond) -> Result<bool, E>,
) -> Result<TopologyBlock, E>
where
    E: From<AtropisomerError>,
{
    // Complete pinned source: MolOps::cleanupAtropisomers(RWMol &, Hybridizations &).
    // RDKit✔️❌: void cleanupAtropisomers(RWMol &mol, MolOps::Hybridizations &hybs) {
    // RDKit✔️❌:   // make sure that ring info is available
    // RDKit✔️❌:   // (defensive, current calls have it available)
    // RDKit✔️❌:   bool needCleanupAtropisomerStereoGroups = false;
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     switch (bond->getStereo()) {
    // RDKit✔️❌:       case Bond::BondStereo::STEREOATROPCW:
    // RDKit✔️❌:       case Bond::BondStereo::STEREOATROPCCW:
    // RDKit✔️❌:         if (checkBond(mol, bond, hybs)) {
    // RDKit✔️❌:           needCleanupAtropisomerStereoGroups = true;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       default:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (needCleanupAtropisomerStereoGroups) {
    // RDKit✔️❌:     Atropisomers::cleanupAtropisomerStereoGroups(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Returning a detached owned value requires one O(atoms+bonds) clone that
    // the source in-place operation does not perform; loop and lookup costs
    // after that clone remain linear/direct-indexed.
    let mut result = topology.clone();
    let mut need_stereo_group_cleanup = false;
    for bond in &mut result.bonds {
        if matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw) && check(bond)? {
            need_stereo_group_cleanup = true;
        }
    }
    if need_stereo_group_cleanup {
        #[cfg(test)]
        cleanup_ring_state_probe::record_group_cleanup();
        result.stereo_groups =
            cleanup_atropisomer_stereo_groups(&result, &AtropisomerAssignment::default())
                .map_err(E::from)?
                .groups;
    }
    result
        .validate()
        .map_err(|source| E::from(AtropisomerError::InvalidTopology { source }))?;
    Ok(result)
}

/// Private two-cause error for the owned-ring-state cleanup wrapper.
///
/// `Rings` carries the canonical finder failure raised while lazily
/// acquiring SSSR for a tagged bond; `Algorithm` carries the existing
/// detached atropisomer causes. `source()` borrows the inner error.
#[derive(Debug, thiserror::Error)]
pub(crate) enum AtropisomerCleanupError {
    #[error("ring finding failed during atropisomer cleanup: {0}")]
    Rings(#[from] crate::RingFindingError),
    #[error("atropisomer cleanup failed: {0}")]
    Algorithm(#[from] AtropisomerError),
}

/// Test-only observation points at the actual wrapper acquisition site.
/// Counters are never reset by production code; production builds contain
/// no counters or hooks.
#[cfg(test)]
pub(crate) mod cleanup_ring_state_probe {
    use std::cell::Cell;

    thread_local! {
        static FIND_CALLS: Cell<u64> = const { Cell::new(0) };
        static CONSUMED_CALLS: Cell<u64> = const { Cell::new(0) };
        static GROUP_CLEANUP_CALLS: Cell<u64> = const { Cell::new(0) };
        static ACQUIRED_ATOM_ROW_PTRS: Cell<Vec<usize>> = const { Cell::new(Vec::new()) };
    }

    pub(crate) fn record_find() {
        FIND_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_consumed() {
        CONSUMED_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_group_cleanup() {
        GROUP_CLEANUP_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_acquired_atom_rows(rows: *const Vec<cosmolkit_model::AtomId>) {
        ACQUIRED_ATOM_ROW_PTRS.with(|buffer| {
            let mut buffer = buffer.take();
            buffer.push(rows as usize);
            ACQUIRED_ATOM_ROW_PTRS.with(|cell| cell.set(buffer));
        });
    }

    pub(crate) fn find_calls() -> u64 {
        FIND_CALLS.with(Cell::get)
    }

    pub(crate) fn consumed_calls() -> u64 {
        CONSUMED_CALLS.with(Cell::get)
    }

    pub(crate) fn group_cleanup_calls() -> u64 {
        GROUP_CLEANUP_CALLS.with(Cell::get)
    }

    pub(crate) fn acquired_atom_row_pointers() -> Vec<usize> {
        ACQUIRED_ATOM_ROW_PTRS.with(Cell::take)
    }
}

/// Owned-ring-state cleanup used by the sanitize CLEANUP_ATROPISOMERS stage.
///
/// The supplied `Option<RingInfo>` is the one state carrier: it is moved in
/// and moved back out; no full RingInfo clone serves this path. With no
/// AtropCw/AtropCcw bonds the callback never runs: NO ring acquisition or
/// validation happens and the supplied state is preserved exactly (a
/// malformed no-tag state remains unconsumed). For a tagged bond the source
/// checkBond guard runs first: a state below SSSR-or-better acquires ONE
/// canonical find_sssr from the ORIGINAL topology (cleanup changes only
/// stereo tags in its clone, never graph/bond orders/endpoints/IDs) and
/// replaces the carrier by move BEFORE any Sp2 endpoint short-circuit; an
/// initialized-empty SSSR/Symm satisfies the guard and is never re-found.
/// The consumed state is then dimension/index validated and the existing
/// detached check primitive decides the tag.
pub(crate) fn cleanup_invalid_atropisomers_with_ring_state(
    topology: &TopologyBlock,
    hybridization: &HybridizationAssignment,
    rings: Option<RingInfo>,
) -> Result<(TopologyBlock, Option<RingInfo>), AtropisomerCleanupError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if hybridization.values.len() != topology.atoms.len() {
        return Err(AtropisomerError::HybridizationAssignmentLength {
            actual: hybridization.values.len(),
            expected: topology.atoms.len(),
        }
        .into());
    }
    let mut state = rings;
    let output =
        cleanup_invalid_atropisomers_engine::<AtropisomerCleanupError>(topology, |bond| {
            #[cfg(test)]
            cleanup_ring_state_probe::record_consumed();
            // Complete pinned source: MolOps.cpp anonymous-namespace checkBond
            // ring-state guard, applied to the owned carrier.
            // RDKit✔️✔️:   if (!mol.getRingInfo()->isSssrOrBetter()) {
            // RDKit✔️✔️:     RDKit::MolOps::findSSSR(mol);
            // RDKit✔️✔️:   }
            let below_quality = !state
                .as_ref()
                .is_some_and(|carrier| carrier.is_sssr_or_better());
            if below_quality {
                let fresh = crate::rings::find_sssr(topology, &crate::RingSearchParams::default())?;
                #[cfg(test)]
                cleanup_ring_state_probe::record_find();
                #[cfg(test)]
                cleanup_ring_state_probe::record_acquired_atom_rows(fresh.atom_rings().as_ptr());
                state = Some(fresh);
            }
            let carrier = state
                .as_ref()
                .expect("ring state present after the guard branch");
            validate_rings(topology, carrier)?;
            check_invalid_atrop_bond(bond, hybridization, carrier)
                .map_err(AtropisomerCleanupError::from)
        })?;
    Ok((output, state))
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

// L05: the owned-ring-state cleanup wrapper products. Counters are
// thread-local and never reset; all observations are per-call deltas.
// Probe pointers prove MOVE preservation of the one carrier (supplied
// buffers survive with the same heap address; the acquired SSSR address
// recorded at the real acquisition site is the address returned).
#[cfg(test)]
mod cleanup_ring_state_tests {
    use super::cleanup_invalid_atropisomers_with_ring_state;
    use super::cleanup_ring_state_probe as probe;
    use crate::atropisomer::AtropisomerError;
    use crate::{HybridizationAssignment, RingFindType, RingInfo};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, BondStereo, Element, Hybridization};

    fn cycle(n: usize, stereo: BondStereo) -> TopologyBlock {
        let atoms = (0..n)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        let bonds = (0..n)
            .map(|id| {
                let mut spec = BondSpec::new(
                    AtomId::new(id),
                    AtomId::new((id + 1) % n),
                    BondOrder::Single,
                );
                if id == 0 {
                    spec = spec.with_stereo(stereo);
                }
                Bond::from_spec(BondId::new(id), spec)
            })
            .collect::<Vec<_>>();
        // Source-shaped sentinel: one ordered stereo group over all cycle
        // atoms; complete ordered group output is asserted every call.
        let sentinel = StereoGroup::new(
            StereoGroupKind::Or,
            (0..n).map(AtomId::new).collect(),
            Vec::new(),
        )
        .with_id(7);
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), vec![sentinel]).unwrap()
    }

    fn hybs(values: Vec<Hybridization>) -> HybridizationAssignment {
        HybridizationAssignment { values }
    }

    fn all(n: usize, value: Hybridization) -> Vec<Hybridization> {
        vec![value; n]
    }

    fn full_ring(find_type: RingFindType, n: usize) -> RingInfo {
        // Frozen full-row convention: reversed atom/bond order, e.g. the
        // six-cycle rows atoms[0,5,4,3,2,1], bonds[5,4,3,2,1,0].
        let mut info = RingInfo::new(find_type, n, n);
        let rows: Vec<usize> = (0..n).rev().collect();
        info.add_ring(&rows, &rows).unwrap();
        info
    }

    fn reset_state() -> RingInfo {
        let mut info = RingInfo::new(RingFindType::Sssr, 6, 6);
        info.reset();
        info
    }

    fn row_indices(rings: &RingInfo) -> Vec<Vec<usize>> {
        rings
            .atom_rings()
            .iter()
            .map(|row| row.iter().map(|atom| atom.index()).collect())
            .collect()
    }

    fn bond_row_indices(rings: &RingInfo) -> Vec<Vec<usize>> {
        rings
            .bond_rings()
            .iter()
            .map(|row| row.iter().map(|bond| bond.index()).collect())
            .collect()
    }

    #[test]
    fn sanitize_ring_l05_lazy_quality_sixty_call_product() {
        let states: Vec<(&'static str, Option<RingInfo>)> = vec![
            ("none", None),
            ("reset", Some(reset_state())),
            (
                "other-empty",
                Some(RingInfo::new(RingFindType::OtherOrUnknown, 6, 6)),
            ),
            ("fast-empty", Some(RingInfo::new(RingFindType::Fast, 6, 6))),
            ("sssr-empty", Some(RingInfo::new(RingFindType::Sssr, 6, 6))),
            (
                "symm-empty",
                Some(RingInfo::new(RingFindType::SymmSssr, 6, 6)),
            ),
            (
                "other-full6",
                Some(full_ring(RingFindType::OtherOrUnknown, 6)),
            ),
            ("fast-full6", Some(full_ring(RingFindType::Fast, 6))),
            ("sssr-full6", Some(full_ring(RingFindType::Sssr, 6))),
            ("symm-full6", Some(full_ring(RingFindType::SymmSssr, 6))),
        ];
        let mut calls = 0usize;
        for (state_name, state) in states {
            for stereo in [BondStereo::None, BondStereo::AtropCw, BondStereo::AtropCcw] {
                for endpoint in [Hybridization::Sp2, Hybridization::Sp3] {
                    let label = format!("{state_name}/{stereo:?}/{endpoint:?}");
                    let topology = cycle(6, stereo);
                    let groups_snapshot = topology.stereo_groups.clone();
                    let topology_snapshot = topology.clone();
                    let hybridization = hybs(all(6, endpoint));
                    let hybs_snapshot = hybridization.clone();
                    let state_snapshot = state.clone();
                    // The pointer must be captured from the value actually
                    // PASSED (the moved-in carrier), not the loop-local
                    // original, for a real move-preservation proof.
                    let supplied = state.clone();
                    let supplied_ptr = supplied.as_ref().map(|rings| rings.atom_rings().as_ptr());
                    let find_before = probe::find_calls();
                    let consumed_before = probe::consumed_calls();
                    let group_before = probe::group_cleanup_calls();
                    let acquired_before = probe::acquired_atom_row_pointers().len();
                    let (output, returned) = cleanup_invalid_atropisomers_with_ring_state(
                        &topology,
                        &hybridization,
                        supplied,
                    )
                    .unwrap_or_else(|error| panic!("{label}: unexpected error {error:?}"));
                    calls += 1;
                    let find_delta = probe::find_calls() - find_before;
                    let consumed_delta = probe::consumed_calls() - consumed_before;
                    let group_delta = probe::group_cleanup_calls() - group_before;
                    let acquired = probe::acquired_atom_row_pointers();
                    assert_eq!(topology, topology_snapshot, "{label}: input mutated");
                    assert_eq!(hybridization, hybs_snapshot, "{label}: hybs mutated");
                    let tagged = !matches!(stereo, BondStereo::None);
                    let quality = state_snapshot
                        .as_ref()
                        .is_some_and(|rings| rings.is_sssr_or_better());
                    if !tagged {
                        // No tags: NO acquisition, NO validation, NO group
                        // cleanup; the supplied state is preserved exactly.
                        assert_eq!(consumed_delta, 0, "{label}: consumed");
                        assert_eq!(find_delta, 0, "{label}: find");
                        assert_eq!(group_delta, 0, "{label}: group cleanup");
                        assert_eq!(returned, state_snapshot, "{label}: state changed");
                        assert_eq!(output.bonds[0].stereo(), BondStereo::None, "{label}");
                    } else {
                        assert!(consumed_delta >= 1, "{label}: tagged bond not consumed");
                        // Disposition on the six-cycle: Sp3 always clears;
                        // Sp2 clears whenever b0 sits in a size-6 ring
                        // (acquired full6 or supplied full6) and retains only
                        // for supplied initialized-EMPTY SSSR/Symm.
                        let supplied_empty = quality
                            && state_snapshot
                                .as_ref()
                                .is_some_and(|r| r.atom_rings().is_empty());
                        let retains = endpoint == Hybridization::Sp2 && supplied_empty;
                        if quality {
                            assert_eq!(find_delta, 0, "{label}: re-found supplied quality state");
                            assert_eq!(returned, state_snapshot, "{label}: borrowed state changed");
                            assert_eq!(
                                returned.as_ref().map(|r| r.atom_rings().as_ptr()),
                                supplied_ptr,
                                "{label}: supplied buffer not move-preserved"
                            );
                        } else {
                            // Below SSSR: ONE canonical full6 SSSR acquisition
                            // replaces the carrier by move BEFORE the Sp2
                            // short-circuit, whatever the endpoint value.
                            assert_eq!(find_delta, 1, "{label}: acquisition count");
                            let acquired_state = returned
                                .as_ref()
                                .unwrap_or_else(|| panic!("{label}: acquired state missing"));
                            assert_eq!(acquired_state.find_type(), RingFindType::Sssr, "{label}");
                            assert_eq!(
                                row_indices(acquired_state),
                                vec![vec![0, 5, 4, 3, 2, 1]],
                                "{label}"
                            );
                            assert_eq!(
                                bond_row_indices(acquired_state),
                                vec![vec![5, 4, 3, 2, 1, 0]],
                                "{label}"
                            );
                            assert_eq!(acquired.len(), acquired_before + 1, "{label}");
                            assert_eq!(
                                acquired_state.atom_rings().as_ptr() as usize,
                                acquired[acquired_before],
                                "{label}: returned buffer differs from real acquisition site"
                            );
                        }
                        if retains {
                            assert_eq!(output.bonds[0].stereo(), stereo, "{label}");
                            assert_eq!(group_delta, 0, "{label}: group cleanup ran");
                        } else {
                            assert_eq!(output.bonds[0].stereo(), BondStereo::None, "{label}");
                            // Cleared tags alone trigger group cleanup; with
                            // no remaining atrop bonds the ordered group
                            // output is unchanged.
                            assert_eq!(group_delta, 1, "{label}: group cleanup missing");
                            assert_eq!(output.stereo_groups, groups_snapshot, "{label}: groups");
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 60, "exact census");
    }

    #[test]
    fn sanitize_ring_l05_macrocycle_endpoint_forty_eight_call_product() {
        let mk_state = |name: &str, n: usize| -> Option<RingInfo> {
            match name {
                "none" => None,
                "fast-full" => Some(full_ring(RingFindType::Fast, n)),
                "sssr-full" => Some(full_ring(RingFindType::Sssr, n)),
                _ => Some(full_ring(RingFindType::SymmSssr, n)),
            }
        };
        let mut calls = 0usize;
        for (state_name, quality) in [
            ("none", false),
            ("fast-full", false),
            ("sssr-full", true),
            ("symm-full", true),
        ] {
            for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
                for n in [6usize, 8usize] {
                    for endpoint in [0usize, 1usize, 2usize] {
                        let label = format!("{state_name}/{stereo:?}/c{n}/ep{endpoint}");
                        let mut values = all(n, Hybridization::Sp2);
                        if endpoint == 1 {
                            values[0] = Hybridization::Sp3;
                        }
                        if endpoint == 2 {
                            values[1] = Hybridization::Sp3;
                        }
                        let topology = cycle(n, stereo);
                        let groups_snapshot = topology.stereo_groups.clone();
                        let topology_snapshot = topology.clone();
                        let hybridization = hybs(values);
                        let hybs_snapshot = hybridization.clone();
                        let state = mk_state(state_name, n);
                        let state_snapshot = state.clone();
                        let supplied = state.clone();
                        let supplied_ptr = supplied.as_ref().map(|r| r.atom_rings().as_ptr());
                        let find_before = probe::find_calls();
                        let group_before = probe::group_cleanup_calls();
                        let acquired_before = probe::acquired_atom_row_pointers().len();
                        let (output, returned) = cleanup_invalid_atropisomers_with_ring_state(
                            &topology,
                            &hybridization,
                            supplied,
                        )
                        .unwrap_or_else(|error| panic!("{label}: unexpected error {error:?}"));
                        calls += 1;
                        let find_delta = probe::find_calls() - find_before;
                        let group_delta = probe::group_cleanup_calls() - group_before;
                        let acquired = probe::acquired_atom_row_pointers();
                        assert_eq!(topology, topology_snapshot, "{label}: input mutated");
                        assert_eq!(hybridization, hybs_snapshot, "{label}: hybs mutated");
                        if quality {
                            assert_eq!(find_delta, 0, "{label}: re-found");
                            assert_eq!(returned, state_snapshot, "{label}: state changed");
                            assert_eq!(
                                returned.as_ref().map(|r| r.atom_rings().as_ptr()),
                                supplied_ptr,
                                "{label}: buffer not move-preserved"
                            );
                        } else {
                            assert_eq!(find_delta, 1, "{label}: acquisition count");
                            let acquired_state = returned.as_ref().unwrap();
                            assert_eq!(acquired_state.find_type(), RingFindType::Sssr, "{label}");
                            // Frozen row convention: atoms start at 0 and
                            // run backwards to 1; bonds run n-1..0.
                            let expected_atoms: Vec<usize> =
                                std::iter::once(0).chain((1..n).rev()).collect();
                            let expected_bonds: Vec<usize> = (0..n).rev().collect();
                            assert_eq!(
                                row_indices(acquired_state),
                                vec![expected_atoms],
                                "{label}"
                            );
                            assert_eq!(
                                bond_row_indices(acquired_state),
                                vec![expected_bonds],
                                "{label}"
                            );
                            assert_eq!(acquired.len(), acquired_before + 1, "{label}");
                            assert_eq!(
                                acquired_state.atom_rings().as_ptr() as usize,
                                acquired[acquired_before],
                                "{label}: acquired buffer mismatch"
                            );
                        }
                        // Six-cycle always clears; eight-cycle retains only
                        // with BOTH endpoints Sp2.
                        let retains = n == 8 && endpoint == 0;
                        if retains {
                            assert_eq!(output.bonds[0].stereo(), stereo, "{label}");
                            assert_eq!(group_delta, 0, "{label}: group cleanup ran");
                        } else {
                            assert_eq!(output.bonds[0].stereo(), BondStereo::None, "{label}");
                            assert_eq!(group_delta, 1, "{label}: group cleanup missing");
                            assert_eq!(output.stereo_groups, groups_snapshot, "{label}: groups");
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 48, "exact census");
    }

    #[test]
    fn sanitize_ring_l05_malformed_inputs_and_vocabulary_controls() {
        // Bad hybridization length is rejected before the engine runs.
        let tagged = cycle(6, BondStereo::AtropCw);
        let short = hybs(all(5, Hybridization::Sp2));
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(&tagged, &short, None),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::HybridizationAssignmentLength {
                    actual: 5,
                    expected: 6
                }
            ))
        ));
        // Invalid topology is rejected first in the existing order.
        let invalid = TopologyBlock {
            adjacency: cosmolkit_model::AdjacencyList::from_topology(0, &[]),
            ..tagged.clone()
        };
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(
                &invalid,
                &hybs(all(6, Hybridization::Sp2)),
                None
            ),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::InvalidTopology { .. }
            ))
        ));
        // A CONSUMED SSSR with wrong membership dimensions is rejected by
        // the existing validate_rings causes.
        let atom_rows_wrong = RingInfo::new(RingFindType::Sssr, 7, 6);
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(
                &tagged,
                &hybs(all(6, Hybridization::Sp2)),
                Some(atom_rows_wrong)
            ),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::RingAtomRowCount {
                    actual: 7,
                    expected: 6
                }
            ))
        ));
        let bond_rows_wrong = RingInfo::new(RingFindType::Sssr, 6, 5);
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(
                &tagged,
                &hybs(all(6, Hybridization::Sp2)),
                Some(bond_rows_wrong)
            ),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::RingBondRowCount {
                    actual: 5,
                    expected: 6
                }
            ))
        ));
        // Out-of-range membership cannot reach the wrapper from valid
        // RingInfo values: add_ring GROWS membership rows, and the
        // persisted-components constructor rejects out-of-range rows at
        // construction, so the dimension check above always fires first.
        // That structural fact is recorded here instead of an unreachable
        // assertion.
        // No-tag malformed supplied state stays unconsumed and is returned
        // exactly.
        let untagged = cycle(6, BondStereo::None);
        let malformed = RingInfo::new(RingFindType::Sssr, 7, 6);
        let malformed_snapshot = malformed.clone();
        let find_before = probe::find_calls();
        let (output, returned) = cleanup_invalid_atropisomers_with_ring_state(
            &untagged,
            &hybs(all(6, Hybridization::Sp2)),
            Some(malformed),
        )
        .unwrap();
        assert_eq!(returned, Some(malformed_snapshot), "no-tag state consumed");
        assert_eq!(probe::find_calls() - find_before, 0, "no-tag find");
        assert_eq!(output.bonds[0].stereo(), BondStereo::None);
        // Vocabulary control ONLY: constructed errors prove Display and the
        // borrowed Error::source chain, never a real acquisition failure.
        // A canonical find failure is unreachable from these valid inputs.
        let rings_error = crate::RingFindingError::Value { message: "probe" };
        let cleanup_error = super::AtropisomerCleanupError::Rings(rings_error);
        assert_eq!(
            cleanup_error.to_string(),
            "ring finding failed during atropisomer cleanup: probe"
        );
        assert!(std::error::Error::source(&cleanup_error).is_some());
        let algorithm_error =
            super::AtropisomerCleanupError::Algorithm(AtropisomerError::RingInfoNotSssr);
        assert_eq!(
            algorithm_error.to_string(),
            "atropisomer cleanup failed: ring information must be SSSR or better"
        );
        assert!(std::error::Error::source(&algorithm_error).is_some());
    }
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
