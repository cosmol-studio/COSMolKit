// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{collections::BTreeSet, f64::consts::PI};

use cosmolkit_model::{
    AtomId, Bond, BondId, BondValueError, Conformer3D, CoordinateValidationError, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

use crate::RingInfo;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondControl {
    Atom(AtomId),
    Implicit,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondStereoSpecified {
    Unspecified,
    Specified,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondStereoDescriptor {
    Cis,
    Trans,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereoInfo {
    pub bond: BondId,
    pub controlling_atoms: [DoubleBondControl; 4],
    pub specified: DoubleBondStereoSpecified,
    pub descriptor: Option<DoubleBondStereoDescriptor>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DoubleBondStereoAssignment {
    pub topology: TopologyBlock,
    pub has_unassigned: bool,
    pub assigned_any: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DoubleBondStereoError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("invalid bond value: {0}")]
    BondValue(#[from] BondValueError),
    #[error("bond {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("bond {bond} has order {order:?}, expected Double")]
    NotDoubleBond { bond: BondId, order: BondOrder },
    #[error("bond {bond} has undefined stereo {stereo:?}")]
    UndefinedStereo { bond: BondId, stereo: BondStereo },
    #[error("bond {bond} has unsupported stereo {stereo:?}")]
    UnsupportedStereo { bond: BondId, stereo: BondStereo },
    #[error("bond {bond} {endpoint} endpoint has invalid degree {degree}")]
    InvalidEndpointDegree {
        bond: BondId,
        endpoint: &'static str,
        degree: usize,
    },
    #[error("bond direction {direction:?} is not slash/backslash stereo")]
    InvalidDirection { direction: BondDirection },
    #[error("rank row has {actual} entries, expected {atom_count}")]
    RankCount { actual: usize, atom_count: usize },
    #[error("ring information is not initialized")]
    RingInfoNotInitialized,
    #[error("ring {dimension} rows have length {actual}, expected {expected}")]
    RingRowCount {
        dimension: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error("invalid conformer: {0}")]
    InvalidConformer(#[from] CoordinateValidationError),
    #[error("bond {bond} stereo {stereo:?} requires two stereo references")]
    StereoReferenceRequired { bond: BondId, stereo: BondStereo },
    #[error("bond {bond} {endpoint} reference {reference} is not an endpoint neighbor")]
    StereoReferenceNotNeighbor {
        bond: BondId,
        endpoint: &'static str,
        reference: AtomId,
    },
    #[error("bond {bond} {endpoint} reference {reference} is the opposite endpoint")]
    StereoReferenceIsOppositeEndpoint {
        bond: BondId,
        endpoint: &'static str,
        reference: AtomId,
    },
    #[error("bond {bond} {endpoint} controlling atom {reference} does not match either slot")]
    ControllingAtomMismatch {
        bond: BondId,
        endpoint: &'static str,
        reference: AtomId,
    },
    #[error("bond {bond} {endpoint} endpoint degree {degree} exceeds two controlling slots")]
    ControllingSlotOverflow {
        bond: BondId,
        endpoint: &'static str,
        degree: usize,
    },
}

#[must_use]
pub const fn has_stereo_bond_direction(direction: BondDirection) -> bool {
    // BEGIN RDKIT CPP FUNCTION hasStereoBondDir
    // RDKit✔️✔️: bool hasStereoBondDir(const Bond *bond) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    // RDKit✔️✔️:   return bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
    // RDKit✔️✔️:          bond->getBondDir() == Bond::BondDir::ENDUPRIGHT;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hasStereoBondDir
    matches!(
        direction,
        BondDirection::EndDownRight | BondDirection::EndUpRight
    )
}

pub fn opposite_stereo_bond_direction(
    direction: BondDirection,
) -> Result<BondDirection, DoubleBondStereoError> {
    // BEGIN RDKIT CPP FUNCTION getOppositeBondDir
    // RDKit✔️✔️: Bond::BondDir getOppositeBondDir(Bond::BondDir dir) {
    // RDKit✔️✔️:   PRECONDITION(dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT,
    // RDKit✔️✔️:                "bad bond direction");
    // RDKit✔️✔️:   switch (dir) {
    // RDKit✔️✔️:     case Bond::ENDDOWNRIGHT:
    // RDKit✔️✔️:       return Bond::ENDUPRIGHT;
    // RDKit✔️✔️:     case Bond::ENDUPRIGHT:
    // RDKit✔️✔️:       return Bond::ENDDOWNRIGHT;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return Bond::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getOppositeBondDir
    match direction {
        BondDirection::EndDownRight => Ok(BondDirection::EndUpRight),
        BondDirection::EndUpRight => Ok(BondDirection::EndDownRight),
        direction => Err(DoubleBondStereoError::InvalidDirection { direction }),
    }
}

pub fn should_detect_double_bond_stereo(
    topology: &TopologyBlock,
    rings: &RingInfo,
    bond: BondId,
) -> Result<bool, DoubleBondStereoError> {
    validate_ring_inputs(topology, rings)?;
    let bond_value = checked_bond(topology, bond)?;
    // BEGIN RDKIT CPP FUNCTION shouldDetectDoubleBondStereo
    // RDKit✔️✔️: bool shouldDetectDoubleBondStereo(const Bond *bond) {
    // RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:   return (!ri->numBondRings(bond->getIdx()) ||
    // RDKit✔️✔️:           ri->minBondRingSize(bond->getIdx()) >=
    // RDKit✔️✔️:               Chirality::minRingSizeForDoubleBondStereo);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION shouldDetectDoubleBondStereo
    Ok(bond_value.order() == BondOrder::Double
        && (rings.num_bond_rings(bond) == 0 || rings.min_bond_ring_size(bond) >= 8))
}

pub fn is_double_bond_stereo_candidate(
    topology: &TopologyBlock,
    rings: &RingInfo,
    bond: BondId,
) -> Result<bool, DoubleBondStereoError> {
    validate_ring_inputs(topology, rings)?;
    let bond_value = checked_bond(topology, bond)?;
    // BEGIN RDKIT CPP FUNCTION isBondCandidateForStereo
    // RDKit✔️✔️: bool isBondCandidateForStereo(const Bond *bond) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    // RDKit✔️✔️:   return bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:          bond->getStereo() != Bond::STEREOANY &&
    // RDKit✔️✔️:          bond->getBondDir() != Bond::EITHERDOUBLE &&
    // RDKit✔️✔️:          bond->getBeginAtom()->getDegree() > 1u &&
    // RDKit✔️✔️:          bond->getEndAtom()->getDegree() > 1u &&
    // RDKit✔️✔️:          shouldDetectDoubleBondStereo(bond);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isBondCandidateForStereo
    Ok(bond_value.order() == BondOrder::Double
        && bond_value.stereo() != BondStereo::Any
        && bond_value.direction() != BondDirection::EitherDouble
        && degree(topology, bond_value.begin()) > 1
        && degree(topology, bond_value.end()) > 1
        && (rings.num_bond_rings(bond) == 0 || rings.min_bond_ring_size(bond) >= 8))
}

pub fn neighboring_directed_bond(
    topology: &TopologyBlock,
    atom: AtomId,
) -> Result<Option<BondId>, DoubleBondStereoError> {
    topology.validate()?;
    check_atom(topology, atom)?;
    // BEGIN RDKIT CPP FUNCTION getNeighboringDirectedBond
    // RDKit✔️✔️: const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "no atom");
    // RDKit✔️✔️:   for (const auto &bondIdx :
    // RDKit✔️✔️:        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:     const Bond *bond = mol[bondIdx];
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:         hasStereoBondDir(bond)) {
    // RDKit✔️✔️:       return bond;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getNeighboringDirectedBond
    Ok(topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .find_map(|neighbor| {
            let bond = &topology.bonds[neighbor.bond.index()];
            (bond.order() != BondOrder::Double && has_stereo_bond_direction(bond.direction()))
                .then_some(bond.id())
        }))
}

#[must_use]
pub const fn translate_ez_to_cis_trans(stereo: BondStereo) -> BondStereo {
    // BEGIN RDKIT CPP FUNCTION translateEZLabelToCisTrans
    // RDKit✔️✔️: Bond::BondStereo translateEZLabelToCisTrans(Bond::BondStereo label) {
    // RDKit✔️✔️:   switch (label) {
    // RDKit✔️✔️:     case Bond::STEREOE:
    // RDKit✔️✔️:       return Bond::STEREOTRANS;
    // RDKit✔️✔️:     case Bond::STEREOZ:
    // RDKit✔️✔️:       return Bond::STEREOCIS;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return label;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION translateEZLabelToCisTrans
    match stereo {
        BondStereo::E => BondStereo::Trans,
        BondStereo::Z => BondStereo::Cis,
        stereo => stereo,
    }
}

pub fn find_double_bond_stereo_atoms(
    topology: &TopologyBlock,
    bond: BondId,
    ranks: &[u32],
) -> Result<Option<[AtomId; 2]>, DoubleBondStereoError> {
    topology.validate()?;
    validate_ranks(topology, ranks)?;
    let bond_value = require_double_bond(topology, bond)?;
    // BEGIN RDKIT CPP FUNCTION findStereoAtoms
    // RDKit✔️✔️: INT_VECT findStereoAtoms(const Bond *bond) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(bond->hasOwningMol(), "no mol");
    // RDKit✔️✔️:   PRECONDITION(bond->getBondType() == Bond::DOUBLE, "not double bond");
    // RDKit✔️✔️:   PRECONDITION(bond->getStereo() > Bond::BondStereo::STEREOANY,
    // RDKit✔️✔️:                "no defined stereo");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!bond->getStereoAtoms().empty()) {
    // RDKit✔️✔️:     return bond->getStereoAtoms();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bond->getStereo() == Bond::BondStereo::STEREOE ||
    // RDKit✔️✔️:       bond->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit✔️✔️:     const Atom *startStereoAtom =
    // RDKit✔️✔️:         findHighestCIPNeighbor(bond->getBeginAtom(), bond->getEndAtom());
    // RDKit✔️✔️:     const Atom *endStereoAtom =
    // RDKit✔️✔️:         findHighestCIPNeighbor(bond->getEndAtom(), bond->getBeginAtom());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (startStereoAtom == nullptr || endStereoAtom == nullptr) {
    // RDKit✔️✔️:       return {};
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int startStereoAtomIdx = static_cast<int>(startStereoAtom->getIdx());
    // RDKit✔️✔️:     int endStereoAtomIdx = static_cast<int>(endStereoAtom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return {startStereoAtomIdx, endStereoAtomIdx};
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "Unable to assign stereo atoms for bond "
    // RDKit✔️✔️:                             << bond->getIdx() << std::endl;
    // RDKit✔️✔️:     return {};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findStereoAtoms
    if let Some(references) = bond_value.stereo_atoms() {
        validate_stereo_references(topology, bond_value, references)?;
        return Ok(Some(references));
    }
    match bond_value.stereo() {
        BondStereo::E | BondStereo::Z => {
            let begin = unique_highest_ranked_neighbor(
                topology,
                bond_value.begin(),
                bond_value.end(),
                ranks,
            );
            let end = unique_highest_ranked_neighbor(
                topology,
                bond_value.end(),
                bond_value.begin(),
                ranks,
            );
            Ok(begin.zip(end).map(|(begin, end)| [begin, end]))
        }
        BondStereo::Cis | BondStereo::Trans => {
            Err(DoubleBondStereoError::StereoReferenceRequired {
                bond,
                stereo: bond_value.stereo(),
            })
        }
        BondStereo::None | BondStereo::Any => Err(DoubleBondStereoError::UndefinedStereo {
            bond,
            stereo: bond_value.stereo(),
        }),
        stereo => Err(DoubleBondStereoError::UnsupportedStereo { bond, stereo }),
    }
}

pub fn double_bond_stereo_info(
    topology: &TopologyBlock,
    bond: BondId,
) -> Result<DoubleBondStereoInfo, DoubleBondStereoError> {
    // BEGIN RDKIT CPP FUNCTION getStereoInfo double-bond branch
    // RDKit❗✔️: if (bond->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗✔️:     if (beginAtom->getDegree() < 1 || endAtom->getDegree() < 1 ||
    // RDKit❗✔️:         beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
    // RDKit❗✔️:       throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     sinfo.type = StereoType::Bond_Double;
    // RDKit❗✔️:     sinfo.centeredOn = bond->getIdx();
    // RDKit❗✔️:     sinfo.controllingAtoms.reserve(4);
    // RDKit❗✔️:
    // RDKit❗✔️:     bool seenSquiggleBond = false;
    // RDKit❗✔️:     const auto &mol = bond->getOwningMol();
    // RDKit❗✔️:
    // RDKit❗✔️:     auto explore_bond_end = [&mol, &bond, &sinfo,
    // RDKit❗✔️:                              &seenSquiggleBond](const Atom *atom) {
    // RDKit❗✔️:       for (const auto nbr : mol.atomBonds(atom)) {
    // RDKit❗✔️:         if (nbr->getIdx() != bond->getIdx()) {
    // RDKit❗✔️:           if (nbr->getBondDir() == Bond::BondDir::UNKNOWN) {
    // RDKit❗✔️:             seenSquiggleBond = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           sinfo.controllingAtoms.push_back(
    // RDKit❗✔️:               nbr->getOtherAtomIdx(atom->getIdx()));
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       for (unsigned i = atom->getDegree(); i < 3; ++i) {
    // RDKit❗✔️:         sinfo.controllingAtoms.push_back(Atom::NOATOM);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     };
    // RDKit❗✔️:
    // RDKit❗✔️:     explore_bond_end(beginAtom);
    // RDKit❗✔️:     explore_bond_end(endAtom);
    // RDKit❗✔️:
    // RDKit❗✔️:     if (!seenSquiggleBond) {
    // RDKit❗✔️:       // check to see if either the begin or end atoms has the _UnknownStereo
    // RDKit❗✔️:       // property set. This happens if there was a squiggle bond to an H
    // RDKit❗✔️:       int explicitUnknownStereo = 0;
    // RDKit❗✔️:       if ((bond->getBeginAtom()->getPropIfPresent<int>(
    // RDKit❗✔️:                common_properties::_UnknownStereo, explicitUnknownStereo) &&
    // RDKit❗✔️:            explicitUnknownStereo) ||
    // RDKit❗✔️:           (bond->getEndAtom()->getPropIfPresent<int>(
    // RDKit❗✔️:                common_properties::_UnknownStereo, explicitUnknownStereo) &&
    // RDKit❗✔️:            explicitUnknownStereo)) {
    // RDKit❗✔️:         seenSquiggleBond = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     Bond::BondStereo stereo = bond->getStereo();
    // RDKit❗✔️:     if (stereo == Bond::BondStereo::STEREOANY ||
    // RDKit❗✔️:         bond->getBondDir() == Bond::BondDir::EITHERDOUBLE || seenSquiggleBond) {
    // RDKit❗✔️:       sinfo.specified = Chirality::StereoSpecified::Unknown;
    // RDKit❗✔️:     } else if (stereo != Bond::BondStereo::STEREONONE) {
    // RDKit❗✔️:       if (stereo == Bond::BondStereo::STEREOE ||
    // RDKit❗✔️:           stereo == Bond::BondStereo::STEREOZ) {
    // RDKit❗✔️:         stereo = Chirality::translateEZLabelToCisTrans(stereo);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       sinfo.specified = Chirality::StereoSpecified::Specified;
    // RDKit❗✔️:       const auto satoms = bond->getStereoAtoms();
    // RDKit❗✔️:       if (satoms.size() != 2) {
    // RDKit❗✔️:         throw ValueErrorException("only can support 2 stereo neighbors");
    // RDKit❗✔️:       }
    // RDKit❗✔️:       bool firstAtBegin;
    // RDKit❗✔️:       if (satoms[0] == static_cast<int>(sinfo.controllingAtoms[0])) {
    // RDKit❗✔️:         firstAtBegin = true;
    // RDKit❗✔️:       } else if (satoms[0] == static_cast<int>(sinfo.controllingAtoms[1])) {
    // RDKit❗✔️:         firstAtBegin = false;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         throw ValueErrorException("controlling atom mismatch at begin");
    // RDKit❗✔️:       }
    // RDKit❗✔️:       bool firstAtEnd;
    // RDKit❗✔️:       if (satoms[1] == static_cast<int>(sinfo.controllingAtoms[2])) {
    // RDKit❗✔️:         firstAtEnd = true;
    // RDKit❗✔️:       } else if (satoms[1] == static_cast<int>(sinfo.controllingAtoms[3])) {
    // RDKit❗✔️:         firstAtEnd = false;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         throw ValueErrorException("controlling atom mismatch at end");
    // RDKit❗✔️:       }
    // RDKit❗✔️:       auto mismatch = firstAtBegin ^ firstAtEnd;
    // RDKit❗✔️:       if (mismatch) {
    // RDKit❗✔️:         stereo = (stereo == Bond::BondStereo::STEREOCIS
    // RDKit❗✔️:                       ? Bond::BondStereo::STEREOTRANS
    // RDKit❗✔️:                       : Bond::BondStereo::STEREOCIS);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       switch (stereo) {
    // RDKit❗✔️:         case Bond::BondStereo::STEREOCIS:
    // RDKit❗✔️:           sinfo.descriptor = Chirality::StereoDescriptor::Bond_Cis;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case Bond::BondStereo::STEREOTRANS:
    // RDKit❗✔️:           sinfo.descriptor = Chirality::StereoDescriptor::Bond_Trans;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           UNDER_CONSTRUCTION("unrecognized bond stereo type");
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       sinfo.specified = Chirality::StereoSpecified::Unspecified;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION getStereoInfo double-bond branch
    topology.validate()?;
    let bond_value = require_double_bond(topology, bond)?;
    let begin_controls = controlling_slots(topology, bond_value, bond_value.begin(), "begin")?;
    let end_controls = controlling_slots(topology, bond_value, bond_value.end(), "end")?;
    let controlling_atoms = [
        begin_controls[0],
        begin_controls[1],
        end_controls[0],
        end_controls[1],
    ];
    let seen_unknown = atom_has_unknown_stereo(topology, bond_value.begin())
        || atom_has_unknown_stereo(topology, bond_value.end())
        || incident_bonds(topology, bond_value.begin())
            .chain(incident_bonds(topology, bond_value.end()))
            .filter(|candidate| candidate.id() != bond)
            .any(|candidate| candidate.direction() == BondDirection::Unknown);
    if seen_unknown
        || bond_value.stereo() == BondStereo::Any
        || bond_value.direction() == BondDirection::EitherDouble
    {
        return Ok(DoubleBondStereoInfo {
            bond,
            controlling_atoms,
            specified: DoubleBondStereoSpecified::Unknown,
            descriptor: None,
        });
    }
    if bond_value.stereo() == BondStereo::None {
        return Ok(DoubleBondStereoInfo {
            bond,
            controlling_atoms,
            specified: DoubleBondStereoSpecified::Unspecified,
            descriptor: None,
        });
    }
    let stereo = translate_ez_to_cis_trans(bond_value.stereo());
    if !matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
        return Err(DoubleBondStereoError::UnsupportedStereo {
            bond,
            stereo: bond_value.stereo(),
        });
    }
    let references =
        bond_value
            .stereo_atoms()
            .ok_or(DoubleBondStereoError::StereoReferenceRequired {
                bond,
                stereo: bond_value.stereo(),
            })?;
    validate_stereo_references(topology, bond_value, references)?;
    let begin_second = control_position(begin_controls, references[0], bond, "begin")? == 1;
    let end_second = control_position(end_controls, references[1], bond, "end")? == 1;
    let flipped = begin_second ^ end_second;
    let descriptor = match (stereo, flipped) {
        (BondStereo::Cis, false) | (BondStereo::Trans, true) => DoubleBondStereoDescriptor::Cis,
        (BondStereo::Trans, false) | (BondStereo::Cis, true) => DoubleBondStereoDescriptor::Trans,
        _ => unreachable!(),
    };
    Ok(DoubleBondStereoInfo {
        bond,
        controlling_atoms,
        specified: DoubleBondStereoSpecified::Specified,
        descriptor: Some(descriptor),
    })
}

pub fn with_double_bond_stereo_reference(
    mut topology: TopologyBlock,
    bond: BondId,
    stereo: BondStereo,
    use_cx_ordering: bool,
) -> Result<TopologyBlock, DoubleBondStereoError> {
    topology.validate()?;
    require_double_bond(&topology, bond)?;
    set_stereo_for_bond(&mut topology, bond, stereo, use_cx_ordering)?;
    topology.validate()?;
    Ok(topology)
}

pub fn assign_double_bond_stereo_from_directions(
    mut topology: TopologyBlock,
) -> Result<TopologyBlock, DoubleBondStereoError> {
    topology.validate()?;
    let assignments = topology
        .bonds
        .iter()
        .filter(|bond| bond.order() == BondOrder::Double && bond.stereo() != BondStereo::Any)
        .filter_map(|bond| {
            let begin_directed = neighboring_directed_bond_unchecked(&topology, bond.begin())?;
            let end_directed = neighboring_directed_bond_unchecked(&topology, bond.end())?;
            let begin_reference = other_atom(begin_directed, bond.begin());
            let end_reference = other_atom(end_directed, bond.end());
            let mut begin_direction = begin_directed.direction();
            if begin_directed.begin() == bond.begin() {
                begin_direction = opposite_unchecked(begin_direction);
            }
            let mut end_direction = end_directed.direction();
            if end_directed.end() == bond.end() {
                end_direction = opposite_unchecked(end_direction);
            }
            Some((
                bond.id(),
                [begin_reference, end_reference],
                if begin_direction == end_direction {
                    BondStereo::Trans
                } else {
                    BondStereo::Cis
                },
            ))
        })
        .collect::<Vec<_>>();
    // BEGIN RDKIT CPP FUNCTION setBondStereoFromDirections
    // RDKit✔️✔️: void setBondStereoFromDirections(ROMol &mol) {
    // RDKit✔️✔️:   mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️:   for (Bond *bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bond->getStereo() != Bond::STEREOANY) {
    // RDKit✔️✔️:       const Atom *stereoBondBeginAtom = bond->getBeginAtom();
    // RDKit✔️✔️:       const Atom *stereoBondEndAtom = bond->getEndAtom();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       const Bond *directedBondAtBegin =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondBeginAtom);
    // RDKit✔️✔️:       const Bond *directedBondAtEnd =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondEndAtom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (directedBondAtBegin != nullptr && directedBondAtEnd != nullptr) {
    // RDKit✔️✔️:         unsigned beginSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtBegin->getOtherAtomIdx(stereoBondBeginAtom->getIdx());
    // RDKit✔️✔️:         unsigned endSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtEnd->getOtherAtomIdx(stereoBondEndAtom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:         bond->setStereoAtoms(beginSideStereoAtom, endSideStereoAtom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         auto beginSideBondDirection = directedBondAtBegin->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtBegin->getBeginAtom() == stereoBondBeginAtom) {
    // RDKit✔️✔️:           beginSideBondDirection = getOppositeBondDir(beginSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         auto endSideBondDirection = directedBondAtEnd->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtEnd->getEndAtom() == stereoBondEndAtom) {
    // RDKit✔️✔️:           endSideBondDirection = getOppositeBondDir(endSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (beginSideBondDirection == endSideBondDirection) {
    // RDKit✔️✔️:           bond->setStereo(Bond::STEREOTRANS);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bond->setStereo(Bond::STEREOCIS);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION setBondStereoFromDirections
    for (bond, references, stereo) in assignments {
        topology.bonds[bond.index()].set_stereo_atoms(Some(references));
        topology.bonds[bond.index()].set_stereo(stereo)?;
    }
    topology.validate()?;
    Ok(topology)
}

pub fn assign_directional_double_bond_stereo(
    mut topology: TopologyBlock,
    ranks: &[u32],
    rings: &RingInfo,
) -> Result<DoubleBondStereoAssignment, DoubleBondStereoError> {
    validate_ring_inputs(&topology, rings)?;
    validate_ranks(&topology, ranks)?;
    let mut assignments = Vec::new();
    let mut clear_directions = BTreeSet::new();
    let mut unassigned_bonds = 0usize;
    let mut assigned_any = false;
    for bond in &mut topology.bonds {
        if bond.order() == BondOrder::Double && bond.stereo() == BondStereo::None {
            bond.set_stereo_atoms(None);
        }
    }
    // BEGIN RDKIT CPP FUNCTION assignBondStereoCodes
    // RDKit❗✔️: std::pair<bool, bool> assignBondStereoCodes(ROMol &mol, UINT_VECT &ranks) {
    // RDKit❗✔️:   PRECONDITION((!ranks.size() || ranks.size() == mol.getNumAtoms()),
    // RDKit❗✔️:                "bad rank vector size");
    // RDKit❗✔️:   bool assignedABond = false;
    // RDKit❗✔️:   unsigned int unassignedBonds = 0;
    // RDKit❗✔️:   boost::dynamic_bitset<> bondsToClear(mol.getNumBonds());
    // RDKit❗✔️:   // find the double bonds:
    // RDKit❗✔️:   for (auto dblBond : mol.bonds()) {
    // RDKit❗✔️:     if (dblBond->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗✔️:       if (dblBond->getStereo() != Bond::BondStereo::STEREONONE) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (!ranks.size()) {
    // RDKit❗✔️:         assignAtomCIPRanks(mol, ranks);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       dblBond->getStereoAtoms().clear();
    // RDKit❗✔️:
    // RDKit❗✔️:       // at the moment we are ignoring stereochem on ring bonds with less than
    // RDKit❗✔️:       // 8 members.
    // RDKit❗✔️:       if (shouldDetectDoubleBondStereo(dblBond)) {
    // RDKit❗✔️:         const Atom *begAtom = dblBond->getBeginAtom();
    // RDKit❗✔️:         const Atom *endAtom = dblBond->getEndAtom();
    // RDKit❗✔️:         // we're only going to handle 2 or three coordinate atoms:
    // RDKit❗✔️:         if ((begAtom->getDegree() == 2 || begAtom->getDegree() == 3) &&
    // RDKit❗✔️:             (endAtom->getDegree() == 2 || endAtom->getDegree() == 3)) {
    // RDKit❗✔️:           ++unassignedBonds;
    // RDKit❗✔️:
    // RDKit❗✔️:           // look around each atom and see if it has at least one bond with
    // RDKit❗✔️:           // direction marked:
    // RDKit❗✔️:
    // RDKit❗✔️:           // the pairs here are: atomIdx,bonddir
    // RDKit❗✔️:           Chirality::INT_PAIR_VECT begAtomNeighbors, endAtomNeighbors;
    // RDKit❗✔️:           bool hasExplicitUnknownStereo = false;
    // RDKit❗✔️:           int bgn_stereo = false, end_stereo = false;
    // RDKit❗✔️:           if ((dblBond->getBeginAtom()->getPropIfPresent(
    // RDKit❗✔️:                    common_properties::_UnknownStereo, bgn_stereo) &&
    // RDKit❗✔️:                bgn_stereo) ||
    // RDKit❗✔️:               (dblBond->getEndAtom()->getPropIfPresent(
    // RDKit❗✔️:                    common_properties::_UnknownStereo, end_stereo) &&
    // RDKit❗✔️:                end_stereo)) {
    // RDKit❗✔️:             hasExplicitUnknownStereo = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           Chirality::findAtomNeighborDirHelper(mol, begAtom, dblBond, ranks,
    // RDKit❗✔️:                                                begAtomNeighbors,
    // RDKit❗✔️:                                                hasExplicitUnknownStereo);
    // RDKit❗✔️:           Chirality::findAtomNeighborDirHelper(mol, endAtom, dblBond, ranks,
    // RDKit❗✔️:                                                endAtomNeighbors,
    // RDKit❗✔️:                                                hasExplicitUnknownStereo);
    // RDKit❗✔️:
    // RDKit❗✔️:           if (begAtomNeighbors.size() && endAtomNeighbors.size()) {
    // RDKit❗✔️:             // Each atom has at least one neighboring bond with marked
    // RDKit❗✔️:             // directionality.  Find the highest-ranked directionality
    // RDKit❗✔️:             // on each side:
    // RDKit❗✔️:
    // RDKit❗✔️:             int begDir, endDir, endNbrAid, begNbrAid;
    // RDKit❗✔️:             if (begAtomNeighbors.size() == 1 ||
    // RDKit❗✔️:                 ranks[begAtomNeighbors[0].first] >
    // RDKit❗✔️:                     ranks[begAtomNeighbors[1].first]) {
    // RDKit❗✔️:               begDir = begAtomNeighbors[0].second;
    // RDKit❗✔️:               begNbrAid = begAtomNeighbors[0].first;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               begDir = begAtomNeighbors[1].second;
    // RDKit❗✔️:               begNbrAid = begAtomNeighbors[1].first;
    // RDKit❗✔️:             }
    // RDKit❗✔️:             if (endAtomNeighbors.size() == 1 ||
    // RDKit❗✔️:                 ranks[endAtomNeighbors[0].first] >
    // RDKit❗✔️:                     ranks[endAtomNeighbors[1].first]) {
    // RDKit❗✔️:               endDir = endAtomNeighbors[0].second;
    // RDKit❗✔️:               endNbrAid = endAtomNeighbors[0].first;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               endDir = endAtomNeighbors[1].second;
    // RDKit❗✔️:               endNbrAid = endAtomNeighbors[1].first;
    // RDKit❗✔️:             }
    // RDKit❗✔️:
    // RDKit❗✔️:             bool conflictingBegin =
    // RDKit❗✔️:                 (begAtomNeighbors.size() == 2 &&
    // RDKit❗✔️:                  begAtomNeighbors[0].second == begAtomNeighbors[1].second);
    // RDKit❗✔️:             bool conflictingEnd =
    // RDKit❗✔️:                 (endAtomNeighbors.size() == 2 &&
    // RDKit❗✔️:                  endAtomNeighbors[0].second == endAtomNeighbors[1].second);
    // RDKit❗✔️:             if (conflictingBegin || conflictingEnd) {
    // RDKit❗✔️:               dblBond->setStereo(Bond::STEREONONE);
    // RDKit❗✔️:               BOOST_LOG(rdWarningLog) << "Conflicting single bond directions "
    // RDKit❗✔️:                                          "around double bond at index "
    // RDKit❗✔️:                                       << dblBond->getIdx() << "." << std::endl;
    // RDKit❗✔️:               BOOST_LOG(rdWarningLog) << "  BondStereo set to STEREONONE and "
    // RDKit❗✔️:                                          "single bond directions set to NONE."
    // RDKit❗✔️:                                       << std::endl;
    // RDKit❗✔️:               assignedABond = true;
    // RDKit❗✔️:               if (conflictingBegin) {
    // RDKit❗✔️:                 bondsToClear[mol.getBondBetweenAtoms(begAtomNeighbors[0].first,
    // RDKit❗✔️:                                                      begAtom->getIdx())
    // RDKit❗✔️:                                  ->getIdx()] = 1;
    // RDKit❗✔️:                 bondsToClear[mol.getBondBetweenAtoms(begAtomNeighbors[1].first,
    // RDKit❗✔️:                                                      begAtom->getIdx())
    // RDKit❗✔️:                                  ->getIdx()] = 1;
    // RDKit❗✔️:               }
    // RDKit❗✔️:               if (conflictingEnd) {
    // RDKit❗✔️:                 bondsToClear[mol.getBondBetweenAtoms(endAtomNeighbors[0].first,
    // RDKit❗✔️:                                                      endAtom->getIdx())
    // RDKit❗✔️:                                  ->getIdx()] = 1;
    // RDKit❗✔️:                 bondsToClear[mol.getBondBetweenAtoms(endAtomNeighbors[1].first,
    // RDKit❗✔️:                                                      endAtom->getIdx())
    // RDKit❗✔️:                                  ->getIdx()] = 1;
    // RDKit❗✔️:               }
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               dblBond->getStereoAtoms().push_back(begNbrAid);
    // RDKit❗✔️:               dblBond->getStereoAtoms().push_back(endNbrAid);
    // RDKit❗✔️:               if (hasExplicitUnknownStereo) {
    // RDKit❗✔️:                 dblBond->setStereo(Bond::STEREOANY);
    // RDKit❗✔️:                 assignedABond = true;
    // RDKit❗✔️:               } else if (begDir == endDir) {
    // RDKit❗✔️:                 // In findAtomNeighborDirHelper, we've set up the
    // RDKit❗✔️:                 // bond directions here so that they correspond to
    // RDKit❗✔️:                 // having both single bonds START at the double bond.
    // RDKit❗✔️:                 // This means that if the single bonds point in the same
    // RDKit❗✔️:                 // direction, the bond is cis, "Z"
    // RDKit❗✔️:                 dblBond->setStereo(Bond::STEREOZ);
    // RDKit❗✔️:                 assignedABond = true;
    // RDKit❗✔️:               } else {
    // RDKit❗✔️:                 dblBond->setStereo(Bond::STEREOE);
    // RDKit❗✔️:                 assignedABond = true;
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:             --unassignedBonds;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    // RDKit❗✔️:     if (bondsToClear[i]) {
    // RDKit❗✔️:       mol.getBondWithIdx(i)->setBondDir(Bond::NONE);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return std::make_pair(unassignedBonds > 0, assignedABond);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION assignBondStereoCodes
    for bond in &topology.bonds {
        if bond.order() != BondOrder::Double || bond.stereo() != BondStereo::None {
            continue;
        }
        if rings.num_bond_rings(bond.id()) != 0 && rings.min_bond_ring_size(bond.id()) < 8 {
            continue;
        }
        if !matches!(degree(&topology, bond.begin()), 2 | 3)
            || !matches!(degree(&topology, bond.end()), 2 | 3)
        {
            continue;
        }
        unassigned_bonds += 1;
        let mut explicit_unknown = atom_has_unknown_stereo(&topology, bond.begin())
            || atom_has_unknown_stereo(&topology, bond.end());
        let begin_neighbors = neighbor_directions(
            &topology,
            bond.begin(),
            bond.id(),
            ranks,
            &mut explicit_unknown,
        );
        let end_neighbors = neighbor_directions(
            &topology,
            bond.end(),
            bond.id(),
            ranks,
            &mut explicit_unknown,
        );
        if begin_neighbors.is_empty() || end_neighbors.is_empty() {
            continue;
        }
        let (begin_atom, begin_direction) = highest_ranked_direction(&begin_neighbors, ranks);
        let (end_atom, end_direction) = highest_ranked_direction(&end_neighbors, ranks);
        let conflicting_begin =
            begin_neighbors.len() == 2 && begin_neighbors[0].1 == begin_neighbors[1].1;
        let conflicting_end = end_neighbors.len() == 2 && end_neighbors[0].1 == end_neighbors[1].1;
        if conflicting_begin || conflicting_end {
            assigned_any = true;
            if conflicting_begin {
                clear_directions.extend(direction_bonds(&topology, bond.begin(), bond.id()));
            }
            if conflicting_end {
                clear_directions.extend(direction_bonds(&topology, bond.end(), bond.id()));
            }
        } else {
            assignments.push((
                bond.id(),
                [begin_atom, end_atom],
                if explicit_unknown {
                    BondStereo::Any
                } else if begin_direction == end_direction {
                    BondStereo::Z
                } else {
                    BondStereo::E
                },
            ));
            assigned_any = true;
        }
        unassigned_bonds -= 1;
    }
    for bond in clear_directions {
        topology.bonds[bond.index()].set_direction(BondDirection::None);
    }
    for (bond, references, stereo) in assignments {
        topology.bonds[bond.index()].set_stereo_atoms(Some(references));
        topology.bonds[bond.index()].set_stereo(stereo)?;
    }
    topology.validate()?;
    Ok(DoubleBondStereoAssignment {
        topology,
        has_unassigned: unassigned_bonds > 0,
        assigned_any,
    })
}

pub fn set_double_bond_neighbor_directions(
    mut topology: TopologyBlock,
    rings: &RingInfo,
    conformer: Option<&Conformer3D>,
) -> Result<TopologyBlock, DoubleBondStereoError> {
    validate_ring_inputs(&topology, rings)?;
    if let Some(conformer) = conformer {
        conformer.validate_for_atom_count(topology.atoms.len())?;
    }
    let bond_count = topology.bonds.len();
    let mut single_bond_counts = vec![0usize; bond_count];
    let mut double_bond_neighbors = vec![Vec::<BondId>::new(); bond_count];
    let mut single_bond_neighbors = vec![Vec::<BondId>::new(); bond_count];
    let mut needs_direction = vec![false; bond_count];
    let mut bonds_in_play = Vec::new();
    // BEGIN RDKIT CPP FUNCTION setDoubleBondNeighborDirections
    // RDKit❗✔️: void setDoubleBondNeighborDirections(ROMol &mol, const Conformer *conf) {
    // RDKit❗✔️:   // used to store the number of single bonds a given
    // RDKit❗✔️:   // single bond is adjacent to
    // RDKit❗✔️:   std::vector<unsigned int> singleBondCounts(mol.getNumBonds(), 0);
    // RDKit❗✔️:   std::vector<Bond *> bondsInPlay;
    // RDKit❗✔️:   // keeps track of which single bonds are adjacent to each double bond:
    // RDKit❗✔️:   VECT_INT_VECT dblBondNbrs(mol.getNumBonds());
    // RDKit❗✔️:   // keeps track of which double bonds are adjacent to each single bond:
    // RDKit❗✔️:   VECT_INT_VECT singleBondNbrs(mol.getNumBonds());
    // RDKit❗✔️:   // keeps track of which single bonds need a dir set and which double bonds
    // RDKit❗✔️:   // need to have their neighbors' dirs set
    // RDKit❗✔️:   boost::dynamic_bitset<> needsDir(mol.getNumBonds());
    // RDKit❗✔️:
    // RDKit❗✔️:   // find double bonds that should be considered for
    // RDKit❗✔️:   // stereochemistry
    // RDKit❗✔️:   // NOTE that we are explicitly excluding double bonds in rings
    // RDKit❗✔️:   // with this test.
    // RDKit❗✔️:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit❗✔️:     RDKit::MolOps::symmetrizeSSSR(mol);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (isBondCandidateForStereo(bond)) {
    // RDKit❗✔️:       bool isCandidate = true;
    // RDKit❗✔️:       for (const auto bondAtom : {bond->getBeginAtom(), bond->getEndAtom()}) {
    // RDKit❗✔️:         for (const auto nbrBond : mol.atomBonds(bondAtom)) {
    // RDKit❗✔️:           if (nbrBond->getBondType() == Bond::SINGLE ||
    // RDKit❗✔️:               nbrBond->getBondType() == Bond::AROMATIC) {
    // RDKit❗✔️:             singleBondCounts[nbrBond->getIdx()] += 1;
    // RDKit❗✔️:             auto nbrDir = nbrBond->getBondDir();
    // RDKit❗✔️:             int hasUnknownStereo = 0;
    // RDKit❗✔️:             if (nbrBond->getBeginAtom() == bondAtom &&
    // RDKit❗✔️:                 nbrDir == Bond::BondDir::UNKNOWN &&
    // RDKit❗✔️:                 nbrBond->getPropIfPresent(common_properties::_UnknownStereo,
    // RDKit❗✔️:                                           hasUnknownStereo) &&
    // RDKit❗✔️:                 hasUnknownStereo) {
    // RDKit❗✔️:               // if there's a wiggly bond starting here, then we're not a
    // RDKit❗✔️:               // candidate for stereo
    // RDKit❗✔️:               isCandidate = false;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               needsDir[bond->getIdx()] = 1;
    // RDKit❗✔️:               if (nbrDir == Bond::BondDir::NONE ||
    // RDKit❗✔️:                   nbrDir == Bond::BondDir::ENDDOWNRIGHT ||
    // RDKit❗✔️:                   nbrDir == Bond::BondDir::ENDUPRIGHT) {
    // RDKit❗✔️:                 needsDir[nbrBond->getIdx()] = 1;
    // RDKit❗✔️:                 dblBondNbrs[bond->getIdx()].push_back(nbrBond->getIdx());
    // RDKit❗✔️:                 // the search may seem inefficient, but these vectors are
    // RDKit❗✔️:                 // going to be at most 2 long (with very few exceptions). It's
    // RDKit❗✔️:                 // just not worth using a different data structure
    // RDKit❗✔️:                 if (std::find(singleBondNbrs[nbrBond->getIdx()].begin(),
    // RDKit❗✔️:                               singleBondNbrs[nbrBond->getIdx()].end(),
    // RDKit❗✔️:                               bond->getIdx()) ==
    // RDKit❗✔️:                     singleBondNbrs[nbrBond->getIdx()].end()) {
    // RDKit❗✔️:                   singleBondNbrs[nbrBond->getIdx()].push_back(bond->getIdx());
    // RDKit❗✔️:                 }
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:           if (!isCandidate) {
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (!isCandidate) {
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (isCandidate) {
    // RDKit❗✔️:         bondsInPlay.push_back(bond);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!bondsInPlay.size()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // order the double bonds based on the singleBondCounts of their neighbors:
    // RDKit❗✔️:   std::vector<std::pair<unsigned int, Bond *>> orderedBondsInPlay;
    // RDKit❗✔️:   for (auto dblBond : bondsInPlay) {
    // RDKit❗✔️:     unsigned int countHere =
    // RDKit❗✔️:         std::accumulate(dblBondNbrs[dblBond->getIdx()].begin(),
    // RDKit❗✔️:                         dblBondNbrs[dblBond->getIdx()].end(), 0);
    // RDKit❗✔️:     // and favor double bonds that are *not* in rings. The combination of
    // RDKit❗✔️:     // using the sum above (instead of the max) and this ring-membershipt test
    // RDKit❗✔️:     // seem to fix sf.net issue 3009836
    // RDKit❗✔️:     if (!(mol.getRingInfo()->numBondRings(dblBond->getIdx()))) {
    // RDKit❗✔️:       countHere *= 10;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     orderedBondsInPlay.push_back(std::make_pair(countHere, dblBond));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::sort(orderedBondsInPlay.begin(), orderedBondsInPlay.end());
    // RDKit❗✔️:
    // RDKit❗✔️:   // oof, now loop over the double bonds in that order and
    // RDKit❗✔️:   // update their neighbor directionalities:
    // RDKit❗✔️:   std::vector<std::pair<unsigned int, Bond *>>::reverse_iterator pairIter;
    // RDKit❗✔️:   for (pairIter = orderedBondsInPlay.rbegin();
    // RDKit❗✔️:        pairIter != orderedBondsInPlay.rend(); ++pairIter) {
    // RDKit❗✔️:     // std::cerr << "RESET?: " << pairIter->second->getIdx() << " "
    // RDKit❗✔️:     //           << pairIter->second->getStereo() << std::endl;
    // RDKit❗✔️:     updateDoubleBondNeighbors(mol, pairIter->second, conf, needsDir,
    // RDKit❗✔️:                               singleBondCounts, singleBondNbrs);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION setDoubleBondNeighborDirections
    // BEGIN RDKIT CPP FUNCTION detectBondStereochemistry
    // RDKit✔️✔️: void detectBondStereochemistry(ROMol &mol, int confId) {
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   setDoubleBondNeighborDirections(mol, &conf);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detectBondStereochemistry
    for double_bond in &topology.bonds {
        if !candidate_unchecked(&topology, rings, double_bond) {
            continue;
        }
        let mut candidate = true;
        for endpoint in [double_bond.begin(), double_bond.end()] {
            for neighbor in topology.adjacency.neighbors_of(endpoint.index()) {
                let neighbor_bond = &topology.bonds[neighbor.bond.index()];
                if matches!(
                    neighbor_bond.order(),
                    BondOrder::Single | BondOrder::Aromatic
                ) {
                    single_bond_counts[neighbor.bond.index()] += 1;
                    if neighbor_bond.begin() == endpoint
                        && neighbor_bond.direction() == BondDirection::Unknown
                        && property_is_true(neighbor_bond.prop("_UnknownStereo"))
                    {
                        candidate = false;
                    } else {
                        needs_direction[double_bond.id().index()] = true;
                        if matches!(
                            neighbor_bond.direction(),
                            BondDirection::None
                                | BondDirection::EndDownRight
                                | BondDirection::EndUpRight
                        ) {
                            needs_direction[neighbor.bond.index()] = true;
                            double_bond_neighbors[double_bond.id().index()].push(neighbor.bond);
                            if !single_bond_neighbors[neighbor.bond.index()]
                                .contains(&double_bond.id())
                            {
                                single_bond_neighbors[neighbor.bond.index()].push(double_bond.id());
                            }
                        }
                    }
                }
                if !candidate {
                    break;
                }
            }
            if !candidate {
                break;
            }
        }
        if candidate {
            bonds_in_play.push(double_bond.id());
        }
    }
    let mut ordered = bonds_in_play
        .into_iter()
        .map(|bond| {
            let mut score = double_bond_neighbors[bond.index()]
                .iter()
                .map(|neighbor| neighbor.index())
                .sum::<usize>();
            if rings.num_bond_rings(bond) == 0 {
                score *= 10;
            }
            (score, bond)
        })
        .collect::<Vec<_>>();
    ordered.sort_by_key(|(score, bond)| (*score, bond.index()));
    for (_, bond) in ordered.into_iter().rev() {
        update_double_bond_neighbors(
            &mut topology,
            bond,
            conformer,
            &mut needs_direction,
            &single_bond_counts,
            &single_bond_neighbors,
        )?;
    }
    topology.validate()?;
    Ok(topology)
}

pub fn clear_single_bond_directions(
    mut topology: TopologyBlock,
    only_wedge_flags: bool,
) -> Result<TopologyBlock, DoubleBondStereoError> {
    topology.validate()?;
    // BEGIN RDKIT CPP FUNCTION clearSingleBondDirFlags
    // RDKit✔️✔️: void clearSingleBondDirFlags(ROMol &mol, bool onlyWedgeFlags) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:       if (bond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         bond->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!onlyWedgeFlags ||
    // RDKit✔️✔️:           (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
    // RDKit✔️✔️:            bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:         bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION clearSingleBondDirFlags
    for bond in &mut topology.bonds {
        if bond.order() != BondOrder::Single {
            continue;
        }
        if bond.direction() == BondDirection::Unknown {
            bond.set_prop("_UnknownStereo", "1")?;
        }
        if !only_wedge_flags || !has_stereo_bond_direction(bond.direction()) {
            bond.set_direction(BondDirection::None);
        }
    }
    Ok(topology)
}

pub fn clear_bond_directions(
    mut topology: TopologyBlock,
    only_wedge_type_bond_directions: bool,
) -> Result<TopologyBlock, DoubleBondStereoError> {
    topology.validate()?;
    // BEGIN RDKIT CPP FUNCTION clearDirFlags
    // RDKit✔️✔️: void clearDirFlags(ROMol &mol, bool onlyWedgeTypeBondDirs) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondDir() == Bond::UNKNOWN ||
    // RDKit✔️✔️:         bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit✔️✔️:       bond->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (onlyWedgeTypeBondDirs == false ||
    // RDKit✔️✔️:         (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
    // RDKit✔️✔️:          bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:       bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION clearDirFlags
    // BEGIN RDKIT CPP FUNCTION clearAllBondDirFlags
    // RDKit✔️✔️: void clearAllBondDirFlags(ROMol &mol) { clearDirFlags(mol, false); }
    // END RDKIT CPP FUNCTION clearAllBondDirFlags
    for bond in &mut topology.bonds {
        if matches!(
            bond.direction(),
            BondDirection::Unknown | BondDirection::EitherDouble
        ) {
            bond.set_prop("_UnknownStereo", "1")?;
        }
        if !only_wedge_type_bond_directions || !has_stereo_bond_direction(bond.direction()) {
            bond.set_direction(BondDirection::None);
        }
    }
    Ok(topology)
}

fn validate_ring_inputs(
    topology: &TopologyBlock,
    rings: &RingInfo,
) -> Result<(), DoubleBondStereoError> {
    topology.validate()?;
    if !rings.is_initialized() {
        return Err(DoubleBondStereoError::RingInfoNotInitialized);
    }
    if rings.atom_row_count() != topology.atoms.len() {
        return Err(DoubleBondStereoError::RingRowCount {
            dimension: "atom",
            actual: rings.atom_row_count(),
            expected: topology.atoms.len(),
        });
    }
    if rings.bond_row_count() != topology.bonds.len() {
        return Err(DoubleBondStereoError::RingRowCount {
            dimension: "bond",
            actual: rings.bond_row_count(),
            expected: topology.bonds.len(),
        });
    }
    Ok(())
}

fn validate_ranks(topology: &TopologyBlock, ranks: &[u32]) -> Result<(), DoubleBondStereoError> {
    if ranks.len() != topology.atoms.len() {
        Err(DoubleBondStereoError::RankCount {
            actual: ranks.len(),
            atom_count: topology.atoms.len(),
        })
    } else {
        Ok(())
    }
}

fn checked_bond(topology: &TopologyBlock, bond: BondId) -> Result<&Bond, DoubleBondStereoError> {
    topology
        .bonds
        .get(bond.index())
        .ok_or(DoubleBondStereoError::BondOutOfRange {
            bond,
            bond_count: topology.bonds.len(),
        })
}

fn require_double_bond(
    topology: &TopologyBlock,
    bond: BondId,
) -> Result<&Bond, DoubleBondStereoError> {
    let bond_value = checked_bond(topology, bond)?;
    if bond_value.order() != BondOrder::Double {
        return Err(DoubleBondStereoError::NotDoubleBond {
            bond,
            order: bond_value.order(),
        });
    }
    Ok(bond_value)
}

fn check_atom(topology: &TopologyBlock, atom: AtomId) -> Result<(), DoubleBondStereoError> {
    if atom.index() >= topology.atoms.len() {
        Err(DoubleBondStereoError::AtomOutOfRange {
            atom,
            atom_count: topology.atoms.len(),
        })
    } else {
        Ok(())
    }
}

fn degree(topology: &TopologyBlock, atom: AtomId) -> usize {
    topology.adjacency.neighbors_of(atom.index()).len()
}

fn incident_bonds(topology: &TopologyBlock, atom: AtomId) -> impl Iterator<Item = &Bond> {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| &topology.bonds[neighbor.bond.index()])
}

fn other_atom(bond: &Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom {
        bond.end()
    } else {
        bond.begin()
    }
}

fn atom_has_unknown_stereo(topology: &TopologyBlock, atom: AtomId) -> bool {
    topology.atoms[atom.index()].unknown_stereo()
        || property_is_true(topology.atoms[atom.index()].prop("_UnknownStereo"))
}

fn property_is_true(value: Option<&str>) -> bool {
    value
        .and_then(|value| value.parse::<i32>().ok())
        .is_some_and(|value| value != 0)
}

fn opposite_unchecked(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        _ => unreachable!("caller proved slash/backslash direction"),
    }
}

fn unique_highest_ranked_neighbor(
    topology: &TopologyBlock,
    atom: AtomId,
    skip: AtomId,
    ranks: &[u32],
) -> Option<AtomId> {
    // BEGIN RDKIT CPP FUNCTION findHighestCIPNeighbor
    // RDKit✔️✔️: const Atom *findHighestCIPNeighbor(const Atom *atom, const Atom *skipAtom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned bestCipRank = 0;
    // RDKit✔️✔️:   const Atom *bestCipRankedAtom = nullptr;
    // RDKit✔️✔️:   const auto &mol = atom->getOwningMol();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto neighbor : mol.atomNeighbors(atom)) {
    // RDKit✔️✔️:     if (neighbor == skipAtom) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned cip = 0;
    // RDKit✔️✔️:     if (!neighbor->getPropIfPresent(common_properties::_CIPRank, cip)) {
    // RDKit✔️✔️:       // If at least one of the atoms doesn't have a CIP rank, the highest rank
    // RDKit✔️✔️:       // does not make sense, so return a nullptr.
    // RDKit✔️✔️:       return nullptr;
    // RDKit✔️✔️:     } else if (cip > bestCipRank || bestCipRankedAtom == nullptr) {
    // RDKit✔️✔️:       bestCipRank = cip;
    // RDKit✔️✔️:       bestCipRankedAtom = neighbor;
    // RDKit✔️✔️:     } else if (cip == bestCipRank) {
    // RDKit✔️✔️:       // This also doesn't make sense if there is a tie (if that's possible).
    // RDKit✔️✔️:       // We still keep the best CIP rank in case something better comes around
    // RDKit✔️✔️:       // (also not sure if that's possible).
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       bestCipRankedAtom = nullptr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return bestCipRankedAtom;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findHighestCIPNeighbor
    let mut best = None;
    let mut best_rank = 0;
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let candidate = AtomId::new(neighbor.atom_index);
        if candidate == skip {
            continue;
        }
        let rank = ranks[candidate.index()];
        if best.is_none() || rank > best_rank {
            best = Some(candidate);
            best_rank = rank;
        } else if rank == best_rank {
            best = None;
        }
    }
    best
}

fn validate_stereo_references(
    topology: &TopologyBlock,
    bond: &Bond,
    references: [AtomId; 2],
) -> Result<(), DoubleBondStereoError> {
    for (endpoint, center, opposite, reference) in [
        ("begin", bond.begin(), bond.end(), references[0]),
        ("end", bond.end(), bond.begin(), references[1]),
    ] {
        if reference == opposite {
            return Err(DoubleBondStereoError::StereoReferenceIsOppositeEndpoint {
                bond: bond.id(),
                endpoint,
                reference,
            });
        }
        if reference.index() >= topology.atoms.len()
            || !topology
                .adjacency
                .neighbors_of(center.index())
                .iter()
                .any(|neighbor| neighbor.atom_index == reference.index())
        {
            return Err(DoubleBondStereoError::StereoReferenceNotNeighbor {
                bond: bond.id(),
                endpoint,
                reference,
            });
        }
    }
    Ok(())
}

fn controlling_slots(
    topology: &TopologyBlock,
    bond: &Bond,
    endpoint_atom: AtomId,
    endpoint: &'static str,
) -> Result<[DoubleBondControl; 2], DoubleBondStereoError> {
    let endpoint_degree = degree(topology, endpoint_atom);
    if !(1..=3).contains(&endpoint_degree) {
        return Err(DoubleBondStereoError::InvalidEndpointDegree {
            bond: bond.id(),
            endpoint,
            degree: endpoint_degree,
        });
    }
    let controls = topology
        .adjacency
        .neighbors_of(endpoint_atom.index())
        .iter()
        .filter(|neighbor| neighbor.bond != bond.id())
        .map(|neighbor| DoubleBondControl::Atom(AtomId::new(neighbor.atom_index)))
        .collect::<Vec<_>>();
    if controls.len() > 2 {
        return Err(DoubleBondStereoError::ControllingSlotOverflow {
            bond: bond.id(),
            endpoint,
            degree: endpoint_degree,
        });
    }
    Ok([
        controls
            .first()
            .copied()
            .unwrap_or(DoubleBondControl::Implicit),
        controls
            .get(1)
            .copied()
            .unwrap_or(DoubleBondControl::Implicit),
    ])
}

fn control_position(
    slots: [DoubleBondControl; 2],
    reference: AtomId,
    bond: BondId,
    endpoint: &'static str,
) -> Result<usize, DoubleBondStereoError> {
    slots
        .iter()
        .position(|control| *control == DoubleBondControl::Atom(reference))
        .ok_or(DoubleBondStereoError::ControllingAtomMismatch {
            bond,
            endpoint,
            reference,
        })
}

fn set_stereo_for_bond(
    topology: &mut TopologyBlock,
    bond: BondId,
    stereo: BondStereo,
    use_cx_ordering: bool,
) -> Result<(), DoubleBondStereoError> {
    // BEGIN RDKIT CPP FUNCTION setStereoForBond
    // RDKit✔️✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit✔️✔️:                       bool useCXSmilesOrdering) {
    // RDKit✔️✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit✔️✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit✔️✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit✔️✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit✔️✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit✔️✔️:   // atom
    // RDKit✔️✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit✔️✔️:   auto endAtom = bond->getEndAtom();
    // RDKit✔️✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit✔️✔️:     std::swap(begAtom, endAtom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit✔️✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit✔️✔️:       if (nbr == endAtom) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit✔️✔️:       if (nbr == begAtom) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit✔️✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit✔️✔️:       std::swap(begControl, endControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit✔️✔️:     bond->setStereo(stereo);
    // RDKit✔️✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION setStereoForBond
    let snapshot = topology.bonds[bond.index()].clone();
    let (low, high, reversed) = if snapshot.begin().index() <= snapshot.end().index() {
        (snapshot.begin(), snapshot.end(), false)
    } else {
        (snapshot.end(), snapshot.begin(), true)
    };
    if degree(topology, low) <= 1 || degree(topology, high) <= 1 {
        return Ok(());
    }
    let low_control = topology
        .adjacency
        .neighbors_of(low.index())
        .iter()
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .filter(|atom| *atom != high)
        .min_by_key(|atom| atom.index())
        .expect("degree gate proves a control");
    let high_controls = topology
        .adjacency
        .neighbors_of(high.index())
        .iter()
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .filter(|atom| *atom != low);
    let high_control = if use_cx_ordering {
        high_controls.min_by_key(|atom| atom.index())
    } else {
        high_controls.max_by_key(|atom| atom.index())
    }
    .expect("degree gate proves a control");
    let references = if reversed {
        [high_control, low_control]
    } else {
        [low_control, high_control]
    };
    topology.bonds[bond.index()].set_stereo_atoms(Some(references));
    topology.bonds[bond.index()].set_stereo(stereo)?;
    Ok(())
}

fn neighboring_directed_bond_unchecked(topology: &TopologyBlock, atom: AtomId) -> Option<&Bond> {
    incident_bonds(topology, atom).find(|bond| {
        bond.order() != BondOrder::Double && has_stereo_bond_direction(bond.direction())
    })
}

fn neighbor_directions(
    topology: &TopologyBlock,
    atom: AtomId,
    reference: BondId,
    ranks: &[u32],
    explicit_unknown: &mut bool,
) -> Vec<(AtomId, BondDirection)> {
    // BEGIN RDKIT CPP FUNCTION findAtomNeighborDirHelper
    // RDKit❗✔️: void findAtomNeighborDirHelper(const ROMol &mol, const Atom *atom,
    // RDKit❗✔️:                                const Bond *refBond, UINT_VECT &ranks,
    // RDKit❗✔️:                                INT_PAIR_VECT &neighbors,
    // RDKit❗✔️:                                bool &hasExplicitUnknownStereo) {
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   PRECONDITION(refBond, "bad bond");
    // RDKit❗✔️:
    // RDKit❗✔️:   bool seenDir = false;
    // RDKit❗✔️:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit❗✔️:     // check whether this bond is explicitly set to have unknown stereo
    // RDKit❗✔️:     if (!hasExplicitUnknownStereo) {
    // RDKit❗✔️:       int explicit_unknown_stereo;
    // RDKit❗✔️:       if (bond->getBondDir() == Bond::UNKNOWN  // there's a squiggle bond
    // RDKit❗✔️:           || (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
    // RDKit❗✔️:                                           explicit_unknown_stereo) &&
    // RDKit❗✔️:               explicit_unknown_stereo)) {
    // RDKit❗✔️:         hasExplicitUnknownStereo = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     Bond::BondDir dir = bond->getBondDir();
    // RDKit❗✔️:     if (bond->getIdx() != refBond->getIdx()) {
    // RDKit❗✔️:       if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
    // RDKit❗✔️:         seenDir = true;
    // RDKit❗✔️:         // If we're considering the bond "backwards", (i.e. from end
    // RDKit❗✔️:         // to beginning, reverse the effective direction:
    // RDKit❗✔️:         if (atom != bond->getBeginAtom()) {
    // RDKit❗✔️:           if (dir == Bond::ENDDOWNRIGHT) {
    // RDKit❗✔️:             dir = Bond::ENDUPRIGHT;
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             dir = Bond::ENDDOWNRIGHT;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       Atom *nbrAtom = bond->getOtherAtom(atom);
    // RDKit❗✔️:       neighbors.push_back(std::make_pair(nbrAtom->getIdx(), dir));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!seenDir) {
    // RDKit❗✔️:     neighbors.clear();
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     if (neighbors.size() == 2 &&
    // RDKit❗✔️:         ranks[neighbors[0].first] == ranks[neighbors[1].first]) {
    // RDKit❗✔️:       // the two substituents are identical, no stereochemistry here:
    // RDKit❗✔️:       neighbors.clear();
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // it's possible that direction was set only one of the bonds, set the
    // RDKit❗✔️:       // other
    // RDKit❗✔️:       // bond's direction to be reversed:
    // RDKit❗✔️:       if (neighbors[0].second != Bond::ENDDOWNRIGHT &&
    // RDKit❗✔️:           neighbors[0].second != Bond::ENDUPRIGHT) {
    // RDKit❗✔️:         CHECK_INVARIANT(neighbors.size() > 1, "too few neighbors");
    // RDKit❗✔️:         neighbors[0].second = neighbors[1].second == Bond::ENDDOWNRIGHT
    // RDKit❗✔️:                                   ? Bond::ENDUPRIGHT
    // RDKit❗✔️:                                   : Bond::ENDDOWNRIGHT;
    // RDKit❗✔️:       } else if (neighbors.size() > 1 &&
    // RDKit❗✔️:                  neighbors[1].second != Bond::ENDDOWNRIGHT &&
    // RDKit❗✔️:                  neighbors[1].second != Bond::ENDUPRIGHT) {
    // RDKit❗✔️:         neighbors[1].second = neighbors[0].second == Bond::ENDDOWNRIGHT
    // RDKit❗✔️:                                   ? Bond::ENDUPRIGHT
    // RDKit❗✔️:                                   : Bond::ENDDOWNRIGHT;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION findAtomNeighborDirHelper
    let mut result = Vec::new();
    let mut saw_direction = false;
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let bond = &topology.bonds[neighbor.bond.index()];
        if !*explicit_unknown
            && (bond.direction() == BondDirection::Unknown
                || bond.unknown_stereo()
                || property_is_true(bond.prop("_UnknownStereo")))
        {
            *explicit_unknown = true;
        }
        if neighbor.bond == reference {
            continue;
        }
        let mut direction = bond.direction();
        if has_stereo_bond_direction(direction) {
            saw_direction = true;
            if atom != bond.begin() {
                direction = opposite_unchecked(direction);
            }
        }
        result.push((AtomId::new(neighbor.atom_index), direction));
    }
    if !saw_direction
        || result.len() == 2 && ranks[result[0].0.index()] == ranks[result[1].0.index()]
    {
        return Vec::new();
    }
    if !has_stereo_bond_direction(result[0].1) {
        result[0].1 = opposite_unchecked(result[1].1);
    } else if result.len() > 1 && !has_stereo_bond_direction(result[1].1) {
        result[1].1 = opposite_unchecked(result[0].1);
    }
    result
}

fn highest_ranked_direction(
    neighbors: &[(AtomId, BondDirection)],
    ranks: &[u32],
) -> (AtomId, BondDirection) {
    if neighbors.len() == 1 || ranks[neighbors[0].0.index()] > ranks[neighbors[1].0.index()] {
        neighbors[0]
    } else {
        neighbors[1]
    }
}

fn direction_bonds(topology: &TopologyBlock, atom: AtomId, reference: BondId) -> Vec<BondId> {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| neighbor.bond != reference)
        .map(|neighbor| neighbor.bond)
        .collect()
}

fn candidate_unchecked(topology: &TopologyBlock, rings: &RingInfo, bond: &Bond) -> bool {
    bond.order() == BondOrder::Double
        && bond.stereo() != BondStereo::Any
        && bond.direction() != BondDirection::EitherDouble
        && degree(topology, bond.begin()) > 1
        && degree(topology, bond.end()) > 1
        && (rings.num_bond_rings(bond.id()) == 0 || rings.min_bond_ring_size(bond.id()) >= 8)
}

#[derive(Clone, Copy)]
struct Controls {
    primary: Option<BondId>,
    secondary: Option<BondId>,
    squiggle: bool,
}

fn controlling_bonds(
    topology: &TopologyBlock,
    needs_direction: &[bool],
    counts: &[usize],
    double_bond: BondId,
    atom: AtomId,
) -> Controls {
    // BEGIN RDKIT CPP FUNCTION controllingBondFromAtom
    // RDKit❗✔️: void controllingBondFromAtom(const ROMol &mol,
    // RDKit❗✔️:                              const boost::dynamic_bitset<> &needsDir,
    // RDKit❗✔️:                              const std::vector<unsigned int> &singleBondCounts,
    // RDKit❗✔️:                              const Bond *dblBond, const Atom *atom, Bond *&bond,
    // RDKit❗✔️:                              Bond *&obond, bool &squiggleBondSeen,
    // RDKit❗✔️:                              bool &doubleBondSeen) {
    // RDKit❗✔️:   bond = nullptr;
    // RDKit❗✔️:   obond = nullptr;
    // RDKit❗✔️:   for (const auto tBond : mol.atomBonds(atom)) {
    // RDKit❗✔️:     if (tBond == dblBond) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if ((tBond->getBondType() == Bond::SINGLE ||
    // RDKit❗✔️:          tBond->getBondType() == Bond::AROMATIC) &&
    // RDKit❗✔️:         (tBond->getBondDir() == Bond::BondDir::NONE ||
    // RDKit❗✔️:          tBond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
    // RDKit❗✔️:          tBond->getBondDir() == Bond::BondDir::ENDUPRIGHT)) {
    // RDKit❗✔️:       // prefer bonds that already have their directionality set
    // RDKit❗✔️:       // or that are adjacent to more double bonds:
    // RDKit❗✔️:       if (!bond) {
    // RDKit❗✔️:         bond = tBond;
    // RDKit❗✔️:       } else if (needsDir[tBond->getIdx()]) {
    // RDKit❗✔️:         if (singleBondCounts[tBond->getIdx()] >
    // RDKit❗✔️:             singleBondCounts[bond->getIdx()]) {
    // RDKit❗✔️:           obond = bond;
    // RDKit❗✔️:           bond = tBond;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           obond = tBond;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         obond = bond;
    // RDKit❗✔️:         bond = tBond;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (tBond->getBondType() == Bond::DOUBLE) {
    // RDKit❗✔️:       doubleBondSeen = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     int explicit_unknown_stereo;
    // RDKit❗✔️:     if ((tBond->getBondType() == Bond::SINGLE ||
    // RDKit❗✔️:          tBond->getBondType() == Bond::AROMATIC) &&
    // RDKit❗✔️:         (tBond->getBondDir() == Bond::UNKNOWN ||
    // RDKit❗✔️:          ((tBond->getPropIfPresent<int>(common_properties::_UnknownStereo,
    // RDKit❗✔️:                                         explicit_unknown_stereo) &&
    // RDKit❗✔️:            explicit_unknown_stereo)))) {
    // RDKit❗✔️:       squiggleBondSeen = true;
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION controllingBondFromAtom
    let mut primary: Option<BondId> = None;
    let mut secondary: Option<BondId> = None;
    let mut squiggle = false;
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        if neighbor.bond == double_bond {
            continue;
        }
        let bond = &topology.bonds[neighbor.bond.index()];
        if matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
            && matches!(
                bond.direction(),
                BondDirection::None | BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            if let Some(current) = primary {
                if needs_direction[neighbor.bond.index()] {
                    if counts[neighbor.bond.index()] > counts[current.index()] {
                        secondary = primary;
                        primary = Some(neighbor.bond);
                    } else {
                        secondary = Some(neighbor.bond);
                    }
                } else {
                    secondary = primary;
                    primary = Some(neighbor.bond);
                }
            } else {
                primary = Some(neighbor.bond);
            }
        }
        if matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
            && (bond.direction() == BondDirection::Unknown
                || bond.unknown_stereo()
                || property_is_true(bond.prop("_UnknownStereo")))
        {
            squiggle = true;
            break;
        }
    }
    Controls {
        primary,
        secondary,
        squiggle,
    }
}

fn update_double_bond_neighbors(
    topology: &mut TopologyBlock,
    double_bond: BondId,
    conformer: Option<&Conformer3D>,
    needs_direction: &mut [bool],
    counts: &[usize],
    single_bond_neighbors: &[Vec<BondId>],
) -> Result<(), DoubleBondStereoError> {
    // BEGIN RDKIT CPP FUNCTION updateDoubleBondNeighbors
    // RDKit❗✔️: void updateDoubleBondNeighbors(ROMol &mol, Bond *dblBond, const Conformer *conf,
    // RDKit❗✔️:                                boost::dynamic_bitset<> &needsDir,
    // RDKit❗✔️:                                std::vector<unsigned int> &singleBondCounts,
    // RDKit❗✔️:                                const VECT_INT_VECT &singleBondNbrs) {
    // RDKit❗✔️:   // we want to deal only with double bonds:
    // RDKit❗✔️:   PRECONDITION(dblBond, "bad bond");
    // RDKit❗✔️:   PRECONDITION(dblBond->getBondType() == Bond::DOUBLE, "not a double bond");
    // RDKit❗✔️:   if (!needsDir[dblBond->getIdx()]) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   needsDir.set(dblBond->getIdx(), 0);
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<Bond *> followupBonds;
    // RDKit❗✔️:
    // RDKit❗✔️:   Bond *bond1 = nullptr, *obond1 = nullptr;
    // RDKit❗✔️:   bool squiggleBondSeen = false;
    // RDKit❗✔️:   bool doubleBondSeen = false;
    // RDKit❗✔️:
    // RDKit❗✔️:   controllingBondFromAtom(mol, needsDir, singleBondCounts, dblBond,
    // RDKit❗✔️:                           dblBond->getBeginAtom(), bond1, obond1,
    // RDKit❗✔️:                           squiggleBondSeen, doubleBondSeen);
    // RDKit❗✔️:
    // RDKit❗✔️:   // Don't do any direction setting if we've seen a squiggle bond, but do mark
    // RDKit❗✔️:   // the double bond as a crossed bond and return
    // RDKit❗✔️:   if (squiggleBondSeen) {
    // RDKit❗✔️:     Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!bond1) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   Bond *bond2 = nullptr, *obond2 = nullptr;
    // RDKit❗✔️:   controllingBondFromAtom(mol, needsDir, singleBondCounts, dblBond,
    // RDKit❗✔️:                           dblBond->getEndAtom(), bond2, obond2,
    // RDKit❗✔️:                           squiggleBondSeen, doubleBondSeen);
    // RDKit❗✔️:
    // RDKit❗✔️:   // Don't do any direction setting if we've seen a squiggle bond, but do mark
    // RDKit❗✔️:   // the double bond as a crossed bond and return
    // RDKit❗✔️:   if (squiggleBondSeen) {
    // RDKit❗✔️:     Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!bond2) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   CHECK_INVARIANT(bond1 && bond2, "no bonds found");
    // RDKit❗✔️:
    // RDKit❗✔️:   bool sameTorsionDir = false;
    // RDKit❗✔️:   if (conf) {
    // RDKit❗✔️:     RDGeom::Point3D beginP = conf->getAtomPos(dblBond->getBeginAtomIdx());
    // RDKit❗✔️:     RDGeom::Point3D endP = conf->getAtomPos(dblBond->getEndAtomIdx());
    // RDKit❗✔️:     RDGeom::Point3D bond1P =
    // RDKit❗✔️:         conf->getAtomPos(bond1->getOtherAtomIdx(dblBond->getBeginAtomIdx()));
    // RDKit❗✔️:     RDGeom::Point3D bond2P =
    // RDKit❗✔️:         conf->getAtomPos(bond2->getOtherAtomIdx(dblBond->getEndAtomIdx()));
    // RDKit❗✔️:     // check for a linear arrangement of atoms on either end:
    // RDKit❗✔️:     bool linear = false;
    // RDKit❗✔️:     RDGeom::Point3D p1;
    // RDKit❗✔️:     RDGeom::Point3D p2;
    // RDKit❗✔️:     p1 = bond1P - beginP;
    // RDKit❗✔️:     p2 = endP - beginP;
    // RDKit❗✔️:     if (isLinearArrangement(p1, p2)) {
    // RDKit❗✔️:       if (!obond1) {
    // RDKit❗✔️:         linear = true;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // one of the bonds was linear; what about the other one?
    // RDKit❗✔️:         Bond *tBond = bond1;
    // RDKit❗✔️:         bond1 = obond1;
    // RDKit❗✔️:         obond1 = tBond;
    // RDKit❗✔️:         bond1P = conf->getAtomPos(
    // RDKit❗✔️:             bond1->getOtherAtomIdx(dblBond->getBeginAtomIdx()));
    // RDKit❗✔️:         p1 = bond1P - beginP;
    // RDKit❗✔️:         if (isLinearArrangement(p1, p2)) {
    // RDKit❗✔️:           linear = true;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!linear) {
    // RDKit❗✔️:       p1 = bond2P - endP;
    // RDKit❗✔️:       p2 = beginP - endP;
    // RDKit❗✔️:       if (isLinearArrangement(p1, p2)) {
    // RDKit❗✔️:         if (!obond2) {
    // RDKit❗✔️:           linear = true;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           Bond *tBond = bond2;
    // RDKit❗✔️:           bond2 = obond2;
    // RDKit❗✔️:           obond2 = tBond;
    // RDKit❗✔️:           bond2P = conf->getAtomPos(
    // RDKit❗✔️:               bond2->getOtherAtomIdx(dblBond->getEndAtomIdx()));
    // RDKit❗✔️:           p1 = bond2P - beginP;
    // RDKit❗✔️:           if (isLinearArrangement(p1, p2)) {
    // RDKit❗✔️:             linear = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (linear) {
    // RDKit❗✔️:       Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     double ang = RDGeom::computeDihedralAngle(bond1P, beginP, endP, bond2P);
    // RDKit❗✔️:     sameTorsionDir = ang >= M_PI / 2;
    // RDKit❗✔️:     // std::cerr << "   angle: " << ang << " sameTorsionDir: " << sameTorsionDir
    // RDKit❗✔️:     // << "\n";
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     if (dblBond->getStereo() == Bond::STEREOCIS ||
    // RDKit❗✔️:         dblBond->getStereo() == Bond::STEREOZ) {
    // RDKit❗✔️:       sameTorsionDir = false;
    // RDKit❗✔️:     } else if (dblBond->getStereo() == Bond::STEREOTRANS ||
    // RDKit❗✔️:                dblBond->getStereo() == Bond::STEREOE) {
    // RDKit❗✔️:       sameTorsionDir = true;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // if bond1 or bond2 are not to the stereo-controlling atoms, flip
    // RDKit❗✔️:     // our expections of the torsion dir
    // RDKit❗✔️:     int bond1AtomIdx = bond1->getOtherAtomIdx(dblBond->getBeginAtomIdx());
    // RDKit❗✔️:     if (bond1AtomIdx != dblBond->getStereoAtoms()[0] &&
    // RDKit❗✔️:         bond1AtomIdx != dblBond->getStereoAtoms()[1]) {
    // RDKit❗✔️:       sameTorsionDir = !sameTorsionDir;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     int bond2AtomIdx = bond2->getOtherAtomIdx(dblBond->getEndAtomIdx());
    // RDKit❗✔️:     if (bond2AtomIdx != dblBond->getStereoAtoms()[0] &&
    // RDKit❗✔️:         bond2AtomIdx != dblBond->getStereoAtoms()[1]) {
    // RDKit❗✔️:       sameTorsionDir = !sameTorsionDir;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   /*
    // RDKit❗✔️:      Time for some clarificatory text, because this gets really
    // RDKit❗✔️:      confusing really fast.
    // RDKit❗✔️:
    // RDKit❗✔️:      The dihedral angle analysis above is based on viewing things
    // RDKit❗✔️:      with an atom order as follows:
    // RDKit❗✔️:
    // RDKit❗✔️:      1
    // RDKit❗✔️:       \
    // RDKit❗✔️:        2 = 3
    // RDKit❗✔️:             \
    // RDKit❗✔️:              4
    // RDKit❗✔️:
    // RDKit❗✔️:      so dihedrals > 90 correspond to sameDir=true
    // RDKit❗✔️:
    // RDKit❗✔️:      however, the stereochemistry representation is
    // RDKit❗✔️:      based on something more like this:
    // RDKit❗✔️:
    // RDKit❗✔️:      2
    // RDKit❗✔️:       \
    // RDKit❗✔️:        1 = 3
    // RDKit❗✔️:             \
    // RDKit❗✔️:              4
    // RDKit❗✔️:      (i.e. we consider the direction-setting single bonds to be
    // RDKit❗✔️:       starting at the double-bonded atom)
    // RDKit❗✔️:
    // RDKit❗✔️:   */
    // RDKit❗✔️:   bool reverseBondDir = sameTorsionDir;
    // RDKit❗✔️:
    // RDKit❗✔️:   Atom *atom1 = dblBond->getBeginAtom(), *atom2 = dblBond->getEndAtom();
    // RDKit❗✔️:   if (needsDir[bond1->getIdx()]) {
    // RDKit❗✔️:     for (auto bidx : singleBondNbrs[bond1->getIdx()]) {
    // RDKit❗✔️:       // std::cerr << "       neighbor from: " << bond1->getIdx() << " " << bidx
    // RDKit❗✔️:       //           << ": " << needsDir[bidx] << std::endl;
    // RDKit❗✔️:       if (needsDir[bidx]) {
    // RDKit❗✔️:         followupBonds.push_back(mol.getBondWithIdx(bidx));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (needsDir[bond2->getIdx()]) {
    // RDKit❗✔️:     for (auto bidx : singleBondNbrs[bond2->getIdx()]) {
    // RDKit❗✔️:       // std::cerr << "       neighbor from: " << bond2->getIdx() << " " << bidx
    // RDKit❗✔️:       //           << ": " << needsDir[bidx] << std::endl;
    // RDKit❗✔️:       if (needsDir[bidx]) {
    // RDKit❗✔️:         followupBonds.push_back(mol.getBondWithIdx(bidx));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!needsDir[bond1->getIdx()]) {
    // RDKit❗✔️:     if (!needsDir[bond2->getIdx()]) {
    // RDKit❗✔️:       // check that we agree
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       if (bond1->getBeginAtom() != atom1) {
    // RDKit❗✔️:         reverseBondDir = !reverseBondDir;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       setBondDirRelativeToAtom(bond2, atom2, bond1->getBondDir(),
    // RDKit❗✔️:                                reverseBondDir, needsDir);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (!needsDir[bond2->getIdx()]) {
    // RDKit❗✔️:     if (bond2->getBeginAtom() != atom2) {
    // RDKit❗✔️:       reverseBondDir = !reverseBondDir;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     setBondDirRelativeToAtom(bond1, atom1, bond2->getBondDir(), reverseBondDir,
    // RDKit❗✔️:                              needsDir);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     setBondDirRelativeToAtom(bond1, atom1, Bond::ENDDOWNRIGHT, false, needsDir);
    // RDKit❗✔️:     setBondDirRelativeToAtom(bond2, atom2, Bond::ENDDOWNRIGHT, reverseBondDir,
    // RDKit❗✔️:                              needsDir);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   needsDir[bond1->getIdx()] = 0;
    // RDKit❗✔️:   needsDir[bond2->getIdx()] = 0;
    // RDKit❗✔️:   if (obond1 && needsDir[obond1->getIdx()]) {
    // RDKit❗✔️:     setBondDirRelativeToAtom(obond1, atom1, bond1->getBondDir(),
    // RDKit❗✔️:                              bond1->getBeginAtom() == atom1, needsDir);
    // RDKit❗✔️:     needsDir[obond1->getIdx()] = 0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (obond2 && needsDir[obond2->getIdx()]) {
    // RDKit❗✔️:     setBondDirRelativeToAtom(obond2, atom2, bond2->getBondDir(),
    // RDKit❗✔️:                              bond2->getBeginAtom() == atom2, needsDir);
    // RDKit❗✔️:     needsDir[obond2->getIdx()] = 0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (Bond *oDblBond : followupBonds) {
    // RDKit❗✔️:     updateDoubleBondNeighbors(mol, oDblBond, conf, needsDir, singleBondCounts,
    // RDKit❗✔️:                               singleBondNbrs);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION updateDoubleBondNeighbors
    if !needs_direction[double_bond.index()] {
        return Ok(());
    }
    needs_direction[double_bond.index()] = false;
    let double = topology.bonds[double_bond.index()].clone();
    let mut begin_controls = controlling_bonds(
        topology,
        needs_direction,
        counts,
        double_bond,
        double.begin(),
    );
    if begin_controls.squiggle {
        set_stereo_for_bond(topology, double_bond, BondStereo::Any, false)?;
        return Ok(());
    }
    let Some(mut begin_bond) = begin_controls.primary else {
        return Ok(());
    };
    let mut end_controls =
        controlling_bonds(topology, needs_direction, counts, double_bond, double.end());
    if end_controls.squiggle {
        set_stereo_for_bond(topology, double_bond, BondStereo::Any, false)?;
        return Ok(());
    }
    let Some(mut end_bond) = end_controls.primary else {
        return Ok(());
    };
    let mut same_torsion_direction;
    if let Some(conformer) = conformer {
        let coordinates = conformer.coordinates();
        let begin_point = coordinates[double.begin().index()];
        let end_point = coordinates[double.end().index()];
        let mut begin_neighbor = other_atom(&topology.bonds[begin_bond.index()], double.begin());
        let mut end_neighbor = other_atom(&topology.bonds[end_bond.index()], double.end());
        let mut begin_neighbor_point = coordinates[begin_neighbor.index()];
        let mut end_neighbor_point = coordinates[end_neighbor.index()];
        let mut linear = is_linear(
            sub(begin_neighbor_point, begin_point),
            sub(end_point, begin_point),
        );
        if linear {
            if let Some(alternate) = begin_controls.secondary {
                begin_controls.secondary = Some(begin_bond);
                begin_bond = alternate;
                begin_neighbor = other_atom(&topology.bonds[begin_bond.index()], double.begin());
                begin_neighbor_point = coordinates[begin_neighbor.index()];
                linear = is_linear(
                    sub(begin_neighbor_point, begin_point),
                    sub(end_point, begin_point),
                );
            }
        }
        if !linear {
            linear = is_linear(
                sub(end_neighbor_point, end_point),
                sub(begin_point, end_point),
            );
            if linear {
                if let Some(alternate) = end_controls.secondary {
                    end_controls.secondary = Some(end_bond);
                    end_bond = alternate;
                    end_neighbor = other_atom(&topology.bonds[end_bond.index()], double.end());
                    end_neighbor_point = coordinates[end_neighbor.index()];
                    linear = is_linear(
                        sub(end_neighbor_point, begin_point),
                        sub(begin_point, end_point),
                    );
                }
            }
        }
        if linear {
            set_stereo_for_bond(topology, double_bond, BondStereo::Any, false)?;
            return Ok(());
        }
        same_torsion_direction = dihedral(
            begin_neighbor_point,
            begin_point,
            end_point,
            end_neighbor_point,
        ) >= PI / 2.0;
    } else {
        same_torsion_direction = match double.stereo() {
            BondStereo::Cis | BondStereo::Z => false,
            BondStereo::Trans | BondStereo::E => true,
            _ => return Ok(()),
        };
        let references =
            double
                .stereo_atoms()
                .ok_or(DoubleBondStereoError::StereoReferenceRequired {
                    bond: double_bond,
                    stereo: double.stereo(),
                })?;
        validate_stereo_references(topology, &double, references)?;
        let begin_atom = other_atom(&topology.bonds[begin_bond.index()], double.begin());
        if !references.contains(&begin_atom) {
            same_torsion_direction = !same_torsion_direction;
        }
        let end_atom = other_atom(&topology.bonds[end_bond.index()], double.end());
        if !references.contains(&end_atom) {
            same_torsion_direction = !same_torsion_direction;
        }
    }
    let mut reverse = same_torsion_direction;
    let mut followups = Vec::new();
    if needs_direction[begin_bond.index()] {
        followups.extend(
            single_bond_neighbors[begin_bond.index()]
                .iter()
                .copied()
                .filter(|bond| needs_direction[bond.index()]),
        );
    }
    if needs_direction[end_bond.index()] {
        followups.extend(
            single_bond_neighbors[end_bond.index()]
                .iter()
                .copied()
                .filter(|bond| needs_direction[bond.index()]),
        );
    }
    match (
        needs_direction[begin_bond.index()],
        needs_direction[end_bond.index()],
    ) {
        (false, true) => {
            if topology.bonds[begin_bond.index()].begin() != double.begin() {
                reverse = !reverse;
            }
            let direction = topology.bonds[begin_bond.index()].direction();
            set_direction_relative(topology, end_bond, double.end(), direction, reverse)?;
        }
        (true, false) => {
            if topology.bonds[end_bond.index()].begin() != double.end() {
                reverse = !reverse;
            }
            let direction = topology.bonds[end_bond.index()].direction();
            set_direction_relative(topology, begin_bond, double.begin(), direction, reverse)?;
        }
        (true, true) => {
            set_direction_relative(
                topology,
                begin_bond,
                double.begin(),
                BondDirection::EndDownRight,
                false,
            )?;
            set_direction_relative(
                topology,
                end_bond,
                double.end(),
                BondDirection::EndDownRight,
                reverse,
            )?;
        }
        (false, false) => {}
    }
    needs_direction[begin_bond.index()] = false;
    needs_direction[end_bond.index()] = false;
    if let Some(secondary) = begin_controls.secondary
        && needs_direction[secondary.index()]
    {
        let direction = topology.bonds[begin_bond.index()].direction();
        let reverse = topology.bonds[begin_bond.index()].begin() == double.begin();
        set_direction_relative(topology, secondary, double.begin(), direction, reverse)?;
        needs_direction[secondary.index()] = false;
    }
    if let Some(secondary) = end_controls.secondary
        && needs_direction[secondary.index()]
    {
        let direction = topology.bonds[end_bond.index()].direction();
        let reverse = topology.bonds[end_bond.index()].begin() == double.end();
        set_direction_relative(topology, secondary, double.end(), direction, reverse)?;
        needs_direction[secondary.index()] = false;
    }
    for followup in followups {
        update_double_bond_neighbors(
            topology,
            followup,
            conformer,
            needs_direction,
            counts,
            single_bond_neighbors,
        )?;
    }
    Ok(())
}

fn set_direction_relative(
    topology: &mut TopologyBlock,
    bond: BondId,
    atom: AtomId,
    mut direction: BondDirection,
    mut reverse: bool,
) -> Result<(), DoubleBondStereoError> {
    // BEGIN RDKIT CPP FUNCTION setBondDirRelativeToAtom
    // RDKit✔️✔️: void setBondDirRelativeToAtom(Bond *bond, Atom *atom, Bond::BondDir dir,
    // RDKit✔️✔️:                               bool reverse, boost::dynamic_bitset<> &) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(dir == Bond::ENDUPRIGHT || dir == Bond::ENDDOWNRIGHT, "bad dir");
    // RDKit✔️✔️:   PRECONDITION(atom == bond->getBeginAtom() || atom == bond->getEndAtom(),
    // RDKit✔️✔️:                "atom doesn't belong to bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bond->getBeginAtom() != atom) {
    // RDKit✔️✔️:     reverse = !reverse;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (reverse) {
    // RDKit✔️✔️:     dir = (dir == Bond::ENDUPRIGHT ? Bond::ENDDOWNRIGHT : Bond::ENDUPRIGHT);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // to ensure maximum compatibility, even when a bond has unknown stereo (set
    // RDKit✔️✔️:   // explicitly and recorded in _UnknownStereo property), I will still let a
    // RDKit✔️✔️:   // direction to be computed. You must check the _UnknownStereo property to
    // RDKit✔️✔️:   // make sure whether this bond is explicitly set to have no direction info.
    // RDKit✔️✔️:   // This makes sense because the direction info are all derived from
    // RDKit✔️✔️:   // coordinates, the _UnknownStereo property is like extra metadata to be
    // RDKit✔️✔️:   // used with the direction info.
    // RDKit✔️✔️:   bond->setBondDir(dir);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION setBondDirRelativeToAtom
    if !has_stereo_bond_direction(direction) {
        return Err(DoubleBondStereoError::InvalidDirection { direction });
    }
    if topology.bonds[bond.index()].begin() != atom {
        reverse = !reverse;
    }
    if reverse {
        direction = opposite_unchecked(direction);
    }
    topology.bonds[bond.index()].set_direction(direction);
    Ok(())
}

fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
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

fn norm_squared(vector: [f64; 3]) -> f64 {
    dot(vector, vector)
}

fn is_linear(left: [f64; 3], right: [f64; 3]) -> bool {
    // BEGIN RDKIT CPP FUNCTION isLinearArrangement
    // RDKit✔️✔️: bool isLinearArrangement(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2) {
    // RDKit✔️✔️:   double lsq = v1.lengthSq() * v2.lengthSq();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // treat zero length vectors as linear
    // RDKit✔️✔️:   if (lsq < 1.0e-6) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double dotProd = v1.dotProduct(v2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double cos178 =
    // RDKit✔️✔️:       -0.999388;  // == cos(M_PI-0.035), corresponds to a tolerance of 2 degrees
    // RDKit✔️✔️:   return dotProd < cos178 * sqrt(lsq);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isLinearArrangement
    let product = norm_squared(left) * norm_squared(right);
    product < 1.0e-6 || dot(left, right) < -0.999_388 * product.sqrt()
}

fn dihedral(i: [f64; 3], j: [f64; 3], k: [f64; 3], l: [f64; 3]) -> f64 {
    // BEGIN RDKIT CPP FUNCTION computeDihedralAngle
    // RDKit✔️✔️: Point3D begEndVec = pt3 - pt2;
    // RDKit✔️✔️: Point3D begNbrVec = pt1 - pt2;
    // RDKit✔️✔️: Point3D crs1 = begNbrVec.crossProduct(begEndVec);
    // RDKit✔️✔️: Point3D endNbrVec = pt4 - pt3;
    // RDKit✔️✔️: Point3D crs2 = endNbrVec.crossProduct(begEndVec);
    // RDKit✔️✔️: double ang = crs1.angleTo(crs2);
    // RDKit✔️✔️: return ang;
    // END RDKIT CPP FUNCTION computeDihedralAngle
    // Behavior review: unlike computeSignedDihedralAngle, this source helper
    // returns the unsigned angle in [0, pi]; clamping only contains floating
    // roundoff before acos and does not alter an in-range cosine.
    // Complexity review: two cross products, one dot product and one acos are
    // the same constant-time arithmetic shape as Point3D::angleTo.
    let begin_end = sub(k, j);
    let begin_neighbor = sub(i, j);
    let first = cross(begin_neighbor, begin_end);
    let end_neighbor = sub(l, k);
    let second = cross(end_neighbor, begin_end);
    let cosine = dot(first, second) / (norm_squared(first) * norm_squared(second)).sqrt();
    cosine.clamp(-1.0, 1.0).acos()
}
