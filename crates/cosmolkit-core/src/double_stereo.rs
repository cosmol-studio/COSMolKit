// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{collections::BTreeSet, f64::consts::PI};

use cosmolkit_model::{
    AtomId, Bond, BondId, BondValueError, Conformer3D, CoordinateValidationError, PropertyValue,
    TopologyBlock, TopologyValidationError,
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

/// Detached result of source bond stereo/direction updates.
/// `needs_detect_bond_stereo` records an executed ordinary Int(1) property
/// write. False means no write, including retention of any preexisting flag.
/// Consumers apply this event to their existing molecule property carrier.
#[derive(Debug, Clone, PartialEq)]
pub struct DoubleBondStereoUpdate {
    pub topology: TopologyBlock,
    pub needs_detect_bond_stereo: bool,
    /// Native insufficient ring state is replaced before direction dispatch.
    /// None retains the supplied already-Symm carrier.
    pub ring_update: Option<RingInfo>,
}

/// Source diagnostic events from stereo-reference discovery.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoAtomSearchWarning {
    DuplicateCipRank { center: AtomId },
    UnableToAssign { bond: BondId },
}

impl std::fmt::Display for StereoAtomSearchWarning {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // RDKit❗✔️:           << "Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"
        // RDKit❗✔️:           << std::endl;
        // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Unable to assign stereo atoms for bond "
        // RDKit❗✔️:                             << bond->getIdx() << std::endl;
        // Exact message payloads in the modeled index range; callers deliver
        // the newline immediately. Global logger state is not modeled here.
        // Constant payload or one decimal integer, no temporary String.
        match self {
            Self::DuplicateCipRank { .. } => formatter
                .write_str("Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"),
            Self::UnableToAssign { bond } => write!(
                formatter,
                "Unable to assign stereo atoms for bond {}",
                bond.index()
            ),
        }
    }
}

/// Found references plus ordered source warnings.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoAtomSearch {
    pub atoms: Option<[AtomId; 2]>,
    pub warnings: Vec<StereoAtomSearchWarning>,
}

/// Discover source stereo references with a lazy reader of already-present
/// CIP ranks. The reader returns None for absence and propagates source
/// property conversion errors. No rank assignment is performed.
#[doc(hidden)]
pub fn find_double_bond_stereo_atoms_with_rank_reader<E>(
    topology: &TopologyBlock,
    bond: BondId,
    mut rank_reader: impl FnMut(AtomId) -> Result<Option<u32>, E>,
    mut emit_warning: impl FnMut(StereoAtomSearchWarning),
) -> Result<StereoAtomSearch, E>
where
    E: From<DoubleBondStereoError>,
{
    // BEGIN RDKIT CPP FUNCTION findStereoAtoms
    // RDKit❗❌: INT_VECT findStereoAtoms(const Bond *bond) {
    // RDKit❗❌:   PRECONDITION(bond, "bad bond");
    // RDKit❗❌:   PRECONDITION(bond->hasOwningMol(), "no mol");
    // RDKit❗❌:   PRECONDITION(bond->getBondType() == Bond::DOUBLE, "not double bond");
    // RDKit❗❌:   PRECONDITION(bond->getStereo() > Bond::BondStereo::STEREOANY,
    // RDKit❗❌:                "no defined stereo");
    // RDKit❗❌:
    // RDKit❗❌:   if (!bond->getStereoAtoms().empty()) {
    // RDKit❗❌:     return bond->getStereoAtoms();
    // RDKit❗❌:   }
    // RDKit❗❌:   if (bond->getStereo() == Bond::BondStereo::STEREOE ||
    // RDKit❗❌:       bond->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit❗❌:     const Atom *startStereoAtom =
    // RDKit❗❌:         findHighestCIPNeighbor(bond->getBeginAtom(), bond->getEndAtom());
    // RDKit❗❌:     const Atom *endStereoAtom =
    // RDKit❗❌:         findHighestCIPNeighbor(bond->getEndAtom(), bond->getBeginAtom());
    // RDKit❗❌:
    // RDKit❗❌:     if (startStereoAtom == nullptr || endStereoAtom == nullptr) {
    // RDKit❗❌:       return {};
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     int startStereoAtomIdx = static_cast<int>(startStereoAtom->getIdx());
    // RDKit❗❌:     int endStereoAtomIdx = static_cast<int>(endStereoAtom->getIdx());
    // RDKit❗❌:
    // RDKit❗❌:     return {startStereoAtomIdx, endStereoAtomIdx};
    // RDKit❗❌:   } else {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog) << "Unable to assign stereo atoms for bond "
    // RDKit❗❌:                             << bond->getIdx() << std::endl;
    // RDKit❗❌:     return {};
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION findStereoAtoms
    // Detached validation and fixed two-ID references retain additional input
    // constraints versus native owning pointers and INT_VECT. Native logger
    // global filtering/stream state and signed source index width are separate
    // unmodeled boundaries, not reasons to infer rank or suppress errors.
    // Cost: source endpoint scans are O(d), but whole topology validation adds
    // O(V+E) and the retained diagnostic report can allocate. Warning delivery
    // occurs at the source statement, including before a later reader error.
    topology.validate().map_err(DoubleBondStereoError::from)?;
    let bond_value = require_double_bond(topology, bond)?;
    if matches!(bond_value.stereo(), BondStereo::None | BondStereo::Any) {
        return Err(DoubleBondStereoError::UndefinedStereo {
            bond,
            stereo: bond_value.stereo(),
        }
        .into());
    }
    let mut result = StereoAtomSearch {
        atoms: bond_value.stereo_atoms(),
        warnings: Vec::new(),
    };
    if result.atoms.is_some() {
        return Ok(result);
    }
    if matches!(bond_value.stereo(), BondStereo::E | BondStereo::Z) {
        let mut warn = |center| {
            let warning = StereoAtomSearchWarning::DuplicateCipRank { center };
            emit_warning(warning.clone());
            result.warnings.push(warning);
        };
        let begin = highest_ranked_neighbor_from_reader(
            topology,
            bond_value.begin(),
            bond_value.end(),
            &mut rank_reader,
            &mut warn,
        )?;
        // Source evaluates both endpoints even when the first is absent.
        let end = highest_ranked_neighbor_from_reader(
            topology,
            bond_value.end(),
            bond_value.begin(),
            &mut rank_reader,
            &mut warn,
        )?;
        result.atoms = begin.zip(end).map(|(begin, end)| [begin, end]);
    } else {
        let warning = StereoAtomSearchWarning::UnableToAssign { bond };
        emit_warning(warning.clone());
        result.warnings.push(warning);
    }
    Ok(result)
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DoubleBondStereoError {
    #[error("source SymmSSSR preparation failed: {0}")]
    Rings(#[from] crate::RingFindingError),
    #[error(
        "property {property} unsigned value {value} causes positive_overflow converting UInt to signed int at atom {atom:?} bond {bond:?}"
    )]
    UnsignedPropertyOverflow {
        atom: Option<AtomId>,
        bond: Option<BondId>,
        property: &'static str,
        value: u32,
    },
    #[error("property {property} has invalid kind {kind:?} at atom {atom:?} bond {bond:?}")]
    InvalidPropertyKind {
        atom: Option<AtomId>,
        bond: Option<BondId>,
        property: &'static str,
        kind: cosmolkit_model::PropertyValueKind,
    },
    #[error("property {property} numeric read failed at atom {atom:?} bond {bond:?}: {source}")]
    NumericPropertyRead {
        atom: Option<AtomId>,
        bond: Option<BondId>,
        property: &'static str,
        #[source]
        source: crate::PropertyIntReadError,
    },
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
    checked_bond(topology, bond)?;
    // BEGIN RDKIT CPP FUNCTION shouldDetectDoubleBondStereo
    // RDKit✔️✔️: bool shouldDetectDoubleBondStereo(const Bond *bond) {
    // RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:   return (!ri->numBondRings(bond->getIdx()) ||
    // RDKit✔️✔️:           ri->minBondRingSize(bond->getIdx()) >=
    // RDKit✔️✔️:               Chirality::minRingSizeForDoubleBondStereo);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION shouldDetectDoubleBondStereo
    // The native helper does not test bond order or reject short ring rows:
    // numBondRings explicitly returns zero outside the member table. Preserve
    // its initialization precondition only when the caller reaches this read.
    if !rings.is_initialized() {
        return Err(DoubleBondStereoError::RingInfoNotInitialized);
    }
    Ok(rings.num_bond_rings(bond) == 0 || rings.min_bond_ring_size(bond) >= 8)
}

pub fn is_double_bond_stereo_candidate(
    topology: &TopologyBlock,
    rings: &RingInfo,
    bond: BondId,
) -> Result<bool, DoubleBondStereoError> {
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
    // Native && order is observable: the ring precondition is reached only
    // after both degree guards. Reuse the single canonical source helper.
    // O(1) bond/degree guards, then O(number of this bond's ring memberships);
    // no unrelated graph validation, row-length gate, allocation or clone.
    if bond_value.order() != BondOrder::Double
        || bond_value.stereo() == BondStereo::Any
        || bond_value.direction() == BondDirection::EitherDouble
    {
        return Ok(false);
    }
    check_atom(topology, bond_value.begin())?;
    if degree(topology, bond_value.begin()) <= 1 {
        return Ok(false);
    }
    check_atom(topology, bond_value.end())?;
    if degree(topology, bond_value.end()) <= 1 {
        return Ok(false);
    }
    should_detect_double_bond_stereo(topology, rings, bond)
}

/// Canonical detached source traversal over borrowed incident bonds.
#[doc(hidden)]
pub fn neighboring_directed_bond_from_incident<'a, E>(
    bonds: impl IntoIterator<Item = Result<&'a Bond, E>>,
) -> Result<Option<&'a Bond>, E> {
    // RDKit❗✔️: const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(atom, "no atom");
    // RDKit❗✔️:   for (const auto &bondIdx :
    // RDKit❗✔️:        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit❗✔️:     const Bond *bond = mol[bondIdx];
    // RDKit❗✔️:
    // RDKit❗✔️:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit❗✔️:         hasStereoBondDir(bond)) {
    // RDKit❗✔️:       return bond;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nullptr;
    // RDKit❗✔️: }
    // Same O(degree) physical encounter order and first-return/error boundary;
    // the iterator adapts source graph access without ownership or allocations.
    for bond in bonds {
        let bond = bond?;
        if bond.order() != BondOrder::Double && has_stereo_bond_direction(bond.direction()) {
            return Ok(Some(bond));
        }
    }
    Ok(None)
}

pub fn neighboring_directed_bond(
    topology: &TopologyBlock,
    atom: AtomId,
) -> Result<Option<BondId>, DoubleBondStereoError> {
    // Existing detached entry validates all model rows before native traversal.
    // This adds O(V+E) work and a non-native validation/error-order boundary to
    // the source O(degree) neighbor helper, so neither source axis is promoted.
    // Preserve this existing local invariant contract/coverage during the first
    // complete port; its source-facing boundary is explicitly deferred along
    // with other differences until after all source functions are processed.
    topology.validate()?;
    check_atom(topology, atom)?;
    // BEGIN RDKIT CPP FUNCTION getNeighboringDirectedBond
    // RDKit❗❌: const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
    // RDKit❗❌:   PRECONDITION(atom, "no atom");
    // RDKit❗❌:   for (const auto &bondIdx :
    // RDKit❗❌:        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit❗❌:     const Bond *bond = mol[bondIdx];
    // RDKit❗❌:
    // RDKit❗❌:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit❗❌:         hasStereoBondDir(bond)) {
    // RDKit❗❌:       return bond;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return nullptr;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION getNeighboringDirectedBond
    // Complete native traversal: physical incident order, skip actual DOUBLE
    // before the exact canonical two-direction predicate, return first match
    // or None. No sorting, chemical heuristic, extra bond-type filter, clone
    // or allocation occurs in this source traversal itself.
    neighboring_directed_bond_from_incident(
        topology
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .map(|neighbor| Ok::<_, DoubleBondStereoError>(&topology.bonds[neighbor.bond.index()])),
    )
    .map(|bond| bond.map(Bond::id))
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
    // Existing explicit rank-table input adapter retains its structural and
    // reference-neighbor validation contract. The source discovery behavior
    // has one owner: the lazy reader above, with immediate diagnostics.
    topology.validate()?;
    validate_ranks(topology, ranks)?;
    let bond_value = require_double_bond(topology, bond)?;
    if let Some(references) = bond_value.stereo_atoms() {
        validate_stereo_references(topology, bond_value, references)?;
    }
    find_double_bond_stereo_atoms_with_rank_reader(
        topology,
        bond,
        |candidate| Ok::<_, DoubleBondStereoError>(Some(ranks[candidate.index()])),
        |warning| eprintln!("{warning}"),
    )
    .map(|found| found.atoms)
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
    let seen_unknown = incident_bonds(topology, bond_value.begin())
        .chain(incident_bonds(topology, bond_value.end()))
        .filter(|candidate| candidate.id() != bond)
        .any(|candidate| candidate.direction() == BondDirection::Unknown)
        || atom_has_unknown_stereo(topology, bond_value.begin())?
        || atom_has_unknown_stereo(topology, bond_value.end())?;
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
) -> Result<DoubleBondStereoUpdate, DoubleBondStereoError> {
    topology.validate()?;
    require_double_bond(&topology, bond)?;
    let needs_detect_bond_stereo =
        set_stereo_for_bond(&mut topology, bond, stereo, use_cx_ordering)?;
    topology.validate()?;
    Ok(DoubleBondStereoUpdate {
        topology,
        needs_detect_bond_stereo,
        ring_update: None,
    })
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
    let (has_unassigned, assigned_any) =
        assign_directional_double_bond_stereo_source(&mut topology, ranks, rings)?;
    Ok(DoubleBondStereoAssignment {
        topology,
        has_unassigned,
        assigned_any,
    })
}

// The same canonical algorithm borrows actual detached state for native callers;
// ownership-returning callers retain their original signatures without cloning.
#[doc(hidden)]
pub fn assign_directional_double_bond_stereo_source(
    mut topology: &mut TopologyBlock,
    ranks: &[u32],
    rings: &RingInfo,
) -> Result<(bool, bool), DoubleBondStereoError> {
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
        let mut explicit_unknown = atom_has_unknown_stereo(&topology, bond.begin())?
            || atom_has_unknown_stereo(&topology, bond.end())?;
        let begin_neighbors = neighbor_directions(
            &topology,
            bond.begin(),
            bond.id(),
            ranks,
            &mut explicit_unknown,
        )?;
        let end_neighbors = neighbor_directions(
            &topology,
            bond.end(),
            bond.id(),
            ranks,
            &mut explicit_unknown,
        )?;
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
    Ok((unassigned_bonds > 0, assigned_any))
}

pub fn set_double_bond_neighbor_directions(
    mut topology: TopologyBlock,
    rings: &RingInfo,
    conformer: Option<&Conformer3D>,
) -> Result<DoubleBondStereoUpdate, DoubleBondStereoError> {
    validate_ring_inputs(&topology, rings)?;
    if let Some(conformer) = conformer {
        conformer.validate_for_atom_count(topology.atoms.len())?;
    }
    let mut needs_detect_bond_stereo = false;
    let bond_count = topology.bonds.len();
    let mut single_bond_counts = vec![0u32; bond_count];
    let mut double_bond_neighbors = vec![Vec::<BondId>::new(); bond_count];
    let mut single_bond_neighbors = vec![Vec::<BondId>::new(); bond_count];
    let mut needs_direction = vec![false; bond_count];
    let mut bonds_in_play = Vec::new();
    // BEGIN RECOVERY CHEM-14 SOURCE setDoubleBondNeighborDirections
    // RDKit❗❌: void setDoubleBondNeighborDirections(ROMol &mol, const Conformer *conf) {
    // RDKit❗❌:   // used to store the number of single bonds a given
    // RDKit❗❌:   // single bond is adjacent to
    // RDKit❗❌:   std::vector<unsigned int> singleBondCounts(mol.getNumBonds(), 0);
    // RDKit❗❌:   std::vector<Bond *> bondsInPlay;
    // RDKit❗❌:   // keeps track of which single bonds are adjacent to each double bond:
    // RDKit❗❌:   VECT_INT_VECT dblBondNbrs(mol.getNumBonds());
    // RDKit❗❌:   // keeps track of which double bonds are adjacent to each single bond:
    // RDKit❗❌:   VECT_INT_VECT singleBondNbrs(mol.getNumBonds());
    // RDKit❗❌:   // keeps track of which single bonds need a dir set and which double bonds
    // RDKit❗❌:   // need to have their neighbors' dirs set
    // RDKit❗❌:   boost::dynamic_bitset<> needsDir(mol.getNumBonds());
    // RDKit❗❌:
    // RDKit❗❌:   // find double bonds that should be considered for
    // RDKit❗❌:   // stereochemistry
    // RDKit❗❌:   // NOTE that we are explicitly excluding double bonds in rings
    // RDKit❗❌:   // with this test.
    // RDKit❗❌:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit❗❌:     RDKit::MolOps::symmetrizeSSSR(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bond : mol.bonds()) {
    // RDKit❗❌:     if (isBondCandidateForStereo(bond)) {
    // RDKit❗❌:       bool isCandidate = true;
    // RDKit❗❌:       for (const auto bondAtom : {bond->getBeginAtom(), bond->getEndAtom()}) {
    // RDKit❗❌:         for (const auto nbrBond : mol.atomBonds(bondAtom)) {
    // RDKit❗❌:           if (nbrBond->getBondType() == Bond::SINGLE ||
    // RDKit❗❌:               nbrBond->getBondType() == Bond::AROMATIC) {
    // RDKit❗❌:             singleBondCounts[nbrBond->getIdx()] += 1;
    // RDKit❗❌:             auto nbrDir = nbrBond->getBondDir();
    // RDKit❗❌:             int hasUnknownStereo = 0;
    // RDKit❗❌:             if (nbrBond->getBeginAtom() == bondAtom &&
    // RDKit❗❌:                 nbrDir == Bond::BondDir::UNKNOWN &&
    // RDKit❗❌:                 nbrBond->getPropIfPresent(common_properties::_UnknownStereo,
    // RDKit❗❌:                                           hasUnknownStereo) &&
    // RDKit❗❌:                 hasUnknownStereo) {
    // RDKit❗❌:               // if there's a wiggly bond starting here, then we're not a
    // RDKit❗❌:               // candidate for stereo
    // RDKit❗❌:               isCandidate = false;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               needsDir[bond->getIdx()] = 1;
    // RDKit❗❌:               if (nbrDir == Bond::BondDir::NONE ||
    // RDKit❗❌:                   nbrDir == Bond::BondDir::ENDDOWNRIGHT ||
    // RDKit❗❌:                   nbrDir == Bond::BondDir::ENDUPRIGHT) {
    // RDKit❗❌:                 needsDir[nbrBond->getIdx()] = 1;
    // RDKit❗❌:                 dblBondNbrs[bond->getIdx()].push_back(nbrBond->getIdx());
    // RDKit❗❌:                 // the search may seem inefficient, but these vectors are
    // RDKit❗❌:                 // going to be at most 2 long (with very few exceptions). It's
    // RDKit❗❌:                 // just not worth using a different data structure
    // RDKit❗❌:                 if (std::find(singleBondNbrs[nbrBond->getIdx()].begin(),
    // RDKit❗❌:                               singleBondNbrs[nbrBond->getIdx()].end(),
    // RDKit❗❌:                               bond->getIdx()) ==
    // RDKit❗❌:                     singleBondNbrs[nbrBond->getIdx()].end()) {
    // RDKit❗❌:                   singleBondNbrs[nbrBond->getIdx()].push_back(bond->getIdx());
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           if (!isCandidate) {
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         if (!isCandidate) {
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (isCandidate) {
    // RDKit❗❌:         bondsInPlay.push_back(bond);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!bondsInPlay.size()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // order the double bonds based on the singleBondCounts of their neighbors:
    // RDKit❗❌:   std::vector<std::pair<unsigned int, Bond *>> orderedBondsInPlay;
    // RDKit❗❌:   for (auto dblBond : bondsInPlay) {
    // RDKit❗❌:     unsigned int countHere =
    // RDKit❗❌:         std::accumulate(dblBondNbrs[dblBond->getIdx()].begin(),
    // RDKit❗❌:                         dblBondNbrs[dblBond->getIdx()].end(), 0);
    // RDKit❗❌:     // and favor double bonds that are *not* in rings. The combination of
    // RDKit❗❌:     // using the sum above (instead of the max) and this ring-membershipt test
    // RDKit❗❌:     // seem to fix sf.net issue 3009836
    // RDKit❗❌:     if (!(mol.getRingInfo()->numBondRings(dblBond->getIdx()))) {
    // RDKit❗❌:       countHere *= 10;
    // RDKit❗❌:     }
    // RDKit❗❌:     orderedBondsInPlay.push_back(std::make_pair(countHere, dblBond));
    // RDKit❗❌:   }
    // RDKit❗❌:   std::ranges::sort(orderedBondsInPlay, [](const auto &a, const auto &b) {
    // RDKit❗❌:     // sort in decreasing order of priority
    // RDKit❗❌:     if (a.first != b.first) {
    // RDKit❗❌:       return a.first > b.first;
    // RDKit❗❌:     }
    // RDKit❗❌:     // in case of ties, use the bond index to decide the order
    // RDKit❗❌:     return a.second->getIdx() < b.second->getIdx();
    // RDKit❗❌:   });
    // RDKit❗❌:
    // RDKit❗❌:   // oof, now loop over the double bonds in that order and
    // RDKit❗❌:   // update their neighbor directionalities:
    // RDKit❗❌:   for (const auto &pairIter : orderedBondsInPlay) {
    // RDKit❗❌:     updateDoubleBondNeighbors(mol, pairIter.second, conf, needsDir,
    // RDKit❗❌:                               singleBondCounts, singleBondNbrs);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RECOVERY CHEM-14 SOURCE setDoubleBondNeighborDirections
    // BEGIN RDKIT CPP FUNCTION detectBondStereochemistry
    // RDKit✔️✔️: void detectBondStereochemistry(ROMol &mol, int confId) {
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   setDoubleBondNeighborDirections(mol, &conf);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detectBondStereochemistry
    // Complete detached source ring effect, even when no candidate exists.
    // Actual record/product callers prepare their supplied properties via the
    // owned carrier adapter first; already-Symm input causes no extra find.
    let ring_update = if rings.is_symm_sssr() {
        None
    } else {
        Some(crate::symmetrized_sssr(
            &topology,
            &crate::RingSearchParams::default(),
        )?)
    };
    let rings = ring_update.as_ref().unwrap_or(rings);
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
                    single_bond_counts[neighbor.bond.index()] =
                        single_bond_counts[neighbor.bond.index()].wrapping_add(1);
                    if neighbor_bond.begin() == endpoint
                        && neighbor_bond.direction() == BondDirection::Unknown
                        && property_is_true(
                            neighbor_bond.prop("_UnknownStereo"),
                            None,
                            Some(neighbor_bond.id()),
                        )?
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
                .fold(0u32, |sum, neighbor| {
                    sum.wrapping_add(neighbor.index() as u32)
                });
            score = direction_priority_source_score(score, rings.num_bond_rings(bond) != 0);
            (score, bond)
        })
        .collect::<Vec<_>>();
    // Exact changed ordering: unsigned priority descending, bond index ascending,
    // then forward traversal. Existing candidate/cache/coordinate/error behavior
    // and adjacency allocation remain qualified by the full-function markers.
    ordered.sort_unstable_by(|(left_score, left), (right_score, right)| {
        right_score.cmp(left_score).then_with(|| left.cmp(right))
    });
    for (_, bond) in ordered {
        update_double_bond_neighbors(
            &mut topology,
            bond,
            conformer,
            &mut needs_direction,
            &single_bond_counts,
            &single_bond_neighbors,
            &mut needs_detect_bond_stereo,
        )?;
    }
    topology.validate()?;
    Ok(DoubleBondStereoUpdate {
        topology,
        needs_detect_bond_stereo,
        ring_update,
    })
}

fn direction_priority_source_score(count_here: u32, is_ring_bond: bool) -> u32 {
    // RDKit✔️✔️: if (!(mol.getRingInfo()->numBondRings(dblBond->getIdx()))) {
    // RDKit✔️✔️:   countHere *= 10;
    // RDKit✔️✔️: }
    // Source countHere is unsigned int. The sum is kept in the caller; only
    // this multiplication has defined modulo-2^32 overflow semantics.
    if is_ring_bond {
        count_here
    } else {
        count_here.wrapping_mul(10)
    }
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
            bond.set_prop("_UnknownStereo", 1_i32)?;
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
            bond.set_prop("_UnknownStereo", 1_i32)?;
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
    // SmartsWrite.cpp::FragmentSmartsConstruct intentionally supplies this
    // initialized, zero-row cache before Canon's stereo perception.
    // RDKit✔️✔️:   mol.getRingInfo()->reset();
    // RDKit✔️✔️:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);
    // RingInfo.cpp defines reads beyond its membership rows as non-ring:
    // RDKit✔️✔️: unsigned int RingInfo::numBondRings(unsigned int idx) const {
    // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (idx < d_bondMembers.size()) {
    // RDKit✔️✔️:     return rdcast<unsigned int>(d_bondMembers[idx].size());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // Behavior: accept only the complete empty SymmSSSR sentinel here; keep
    // rejecting uninitialized and partially populated mismatched carriers.
    // Complexity: three O(1) metadata reads; no ring finding or allocation.
    if rings.is_symm_sssr() && rings.atom_row_count() == 0 && rings.bond_row_count() == 0 {
        return Ok(());
    }
    let preserved_prefix = (rings.atom_row_count() != topology.atoms.len()
        || rings.bond_row_count() != topology.bonds.len())
        && crate::rings::preserves_appended_terminal_hydrogen_ring_prefix(topology, rings);
    if rings.atom_row_count() != topology.atoms.len() && !preserved_prefix {
        return Err(DoubleBondStereoError::RingRowCount {
            dimension: "atom",
            actual: rings.atom_row_count(),
            expected: topology.atoms.len(),
        });
    }
    if rings.bond_row_count() != topology.bonds.len() && !preserved_prefix {
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

fn atom_has_unknown_stereo(
    topology: &TopologyBlock,
    atom: AtomId,
) -> Result<bool, DoubleBondStereoError> {
    Ok(topology.atoms[atom.index()].unknown_stereo()
        || property_is_true(
            topology.atoms[atom.index()].prop("_UnknownStereo"),
            Some(atom),
            None,
        )?)
}

fn property_is_true(
    value: Option<&PropertyValue>,
    atom: Option<AtomId>,
    bond: Option<BondId>,
) -> Result<bool, DoubleBondStereoError> {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
    // RDKit✔️✔️:     for (const auto &data : _data) {
    // RDKit✔️✔️:       if (data.key == what) {
    // RDKit✔️✔️:         res = from_rdvalue<T>(data.val);
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450
    // RDKit❗✔️: template <>
    // RDKit❗✔️: inline int rdvalue_cast<int>(RDValue_cast_t v) {
    // RDKit❗✔️:   if (rdvalue_is<int>(v)) {
    // RDKit❗✔️:     return v.value.i;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit❗✔️:     return boost::numeric_cast<int>(v.value.u);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   throw std::bad_any_cast();
    // RDKit❗✔️: }
    // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450

    // RDKit❗✔️: int explicitUnknownStereo = 0;
    // RDKit❗✔️: common_properties::_UnknownStereo, explicitUnknownStereo) &&
    // A vector cannot be cast to int; errors propagate at the guarded read.
    // Present Bool/Double are also incompatible with int and must propagate
    // the same structural wrong-tag failure, not become a false stereo flag.
    // Numeric/String branches reuse source arithmetic behavior; no copy.
    let Some(value) = value else {
        return Ok(false);
    };
    crate::property_value_to_int(value)
        .map(|value| value != 0)
        .map_err(|source| match source {
            crate::PropertyIntReadError::UnsignedOverflow { value } => {
                DoubleBondStereoError::UnsignedPropertyOverflow {
                    atom,
                    bond,
                    property: "_UnknownStereo",
                    value,
                }
            }
            crate::PropertyIntReadError::InvalidKind { kind } => {
                DoubleBondStereoError::InvalidPropertyKind {
                    atom,
                    bond,
                    property: "_UnknownStereo",
                    kind,
                }
            }
            source => DoubleBondStereoError::NumericPropertyRead {
                atom,
                bond,
                property: "_UnknownStereo",
                source,
            },
        })
}

fn opposite_unchecked(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        _ => unreachable!("caller proved slash/backslash direction"),
    }
}

fn highest_ranked_neighbor_from_reader<E>(
    topology: &TopologyBlock,
    atom: AtomId,
    skip: AtomId,
    rank_reader: &mut impl FnMut(AtomId) -> Result<Option<u32>, E>,
    warn: &mut impl FnMut(AtomId),
) -> Result<Option<AtomId>, E> {
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
    // Source physical neighbor order is retained. A tie clears only the best
    // pointer; the next present rank is selected even if lower than best_rank,
    // exactly because the source condition includes bestCipRankedAtom==nullptr.
    // Skip precedes the lazy rank read; a missing rank returns immediately,
    // reader errors propagate and duplicate-rank diagnostics remain ordered.
    // Cost: one O(degree) scan, one rank read per non-skipped neighbor, constant
    // local storage, no sorting, allocation, graph clone or eager rank pass.
    // Reader/diagnostic callbacks transport the existing source property/logger
    // state; their canonical callers own conversions and warning delivery.
    let mut best = None;
    let mut best_rank = 0;
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let candidate = AtomId::new(neighbor.atom_index);
        if candidate == skip {
            continue;
        }
        let Some(rank) = rank_reader(candidate)? else {
            return Ok(None);
        };
        if best.is_none() || rank > best_rank {
            best = Some(candidate);
            best_rank = rank;
        } else if rank == best_rank {
            warn(atom);
            best = None;
        }
    }
    Ok(best)
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

/// Source control selection on the caller's validated ordered adjacency.
/// This internal domain boundary returns references without graph ownership,
/// property writes, runtime authority, or a second stereo algorithm.
#[doc(hidden)]
pub fn double_bond_stereo_reference_atoms<D, F, I>(
    bond: BondId,
    atom_count: usize,
    begin: AtomId,
    end: AtomId,
    use_cx_ordering: bool,
    mut degree_of: D,
    mut neighbors_of: F,
) -> Result<Option<[AtomId; 2]>, DoubleBondStereoError>
where
    D: FnMut(AtomId) -> usize,
    F: FnMut(AtomId) -> I,
    I: IntoIterator<Item = AtomId>,
{
    // BEGIN COMPLETE existing CORE setStereoForBond selector projection
    // RDKit❗✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit❗✔️:                       bool useCXSmilesOrdering) {
    // RDKit❗✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit❗✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit❗✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit❗✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit❗✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit❗✔️:   // atom
    // RDKit❗✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit❗✔️:   auto endAtom = bond->getEndAtom();
    // RDKit❗✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit❗✔️:     std::swap(begAtom, endAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit❗✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit❗✔️:       if (nbr == endAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit❗✔️:       if (nbr == begAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit❗✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit❗✔️:       std::swap(begControl, endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit❗✔️:     bond->setStereo(stereo);
    // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END COMPLETE existing CORE setStereoForBond selector projection
    // BEGIN COMPLETE source reference preconditions
    // RDKit❗✔️: void Bond::setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx) {
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       getOwningMol().getBondBetweenAtoms(getBeginAtomIdx(), bgnIdx) != nullptr,
    // RDKit❗✔️:       "bgnIdx not connected to begin atom of bond");
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       getOwningMol().getBondBetweenAtoms(getEndAtomIdx(), endIdx) != nullptr,
    // RDKit❗✔️:       "endIdx not connected to end atom of bond");
    // RDKit❗✔️:
    // RDKit❗✔️:   auto &atoms = getStereoAtoms();
    // RDKit❗✔️:   atoms.clear();
    // RDKit❗✔️:   atoms.push_back(bgnIdx);
    // RDKit❗✔️:   atoms.push_back(endIdx);
    // RDKit❗✔️: }
    // END COMPLETE source reference preconditions
    // Behavior: degree checks follow sorted endpoint order and short-circuit.
    // Source count/zero initializers remain explicit sentinels, not heuristic
    // defaults. All reference validation precedes caller mutation. The caller
    // owns source setStereoAtoms, setStereo and Int flag side effects.
    // Complexity: two ordered degree-bounded scans, O(1) scratch, no Bond or
    // topology clone. Only CORE contains this control-selection algorithm.
    let (low, high) = if begin.index() > end.index() {
        (end, begin)
    } else {
        (begin, end)
    };
    if degree_of(low) <= 1 || degree_of(high) <= 1 {
        return Ok(None);
    }
    let mut low_control = AtomId::new(atom_count);
    for neighbor in neighbors_of(low) {
        if neighbor != high && neighbor.index() < low_control.index() {
            low_control = neighbor;
        }
    }
    let mut high_control = AtomId::new(if use_cx_ordering { atom_count } else { 0 });
    for neighbor in neighbors_of(high) {
        if neighbor == low {
            continue;
        }
        if (use_cx_ordering && neighbor.index() < high_control.index())
            || (!use_cx_ordering && neighbor.index() > high_control.index())
        {
            high_control = neighbor;
        }
    }
    let references = if low != begin {
        [high_control, low_control]
    } else {
        [low_control, high_control]
    };
    for (endpoint, reference) in [("begin", references[0]), ("end", references[1])] {
        if reference.index() >= atom_count {
            return Err(DoubleBondStereoError::StereoReferenceNotNeighbor {
                bond,
                endpoint,
                reference,
            });
        }
    }
    Ok(Some(references))
}

fn set_stereo_for_bond(
    topology: &mut TopologyBlock,
    bond: BondId,
    stereo: BondStereo,
    use_cx_ordering: bool,
) -> Result<bool, DoubleBondStereoError> {
    // BEGIN COMPLETE CORE source stereo write adapter
    // RDKit❗✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit❗✔️:                       bool useCXSmilesOrdering) {
    // RDKit❗✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit❗✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit❗✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit❗✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit❗✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit❗✔️:   // atom
    // RDKit❗✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit❗✔️:   auto endAtom = bond->getEndAtom();
    // RDKit❗✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit❗✔️:     std::swap(begAtom, endAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit❗✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit❗✔️:       if (nbr == endAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit❗✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit❗✔️:       if (nbr == begAtom) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit❗✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit❗✔️:       std::swap(begControl, endControl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit❗✔️:     bond->setStereo(stereo);
    // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END COMPLETE CORE source stereo write adapter
    let value = &topology.bonds[bond.index()];
    let references = double_bond_stereo_reference_atoms(
        bond,
        topology.atoms.len(),
        value.begin(),
        value.end(),
        use_cx_ordering,
        |atom| degree(topology, atom),
        |atom| {
            topology
                .adjacency
                .neighbors_of(atom.index())
                .iter()
                .map(|neighbor| AtomId::new(neighbor.atom_index))
        },
    )?;
    let Some(references) = references else {
        return Ok(false);
    };
    topology.bonds[bond.index()].set_stereo_atoms(Some(references));
    topology.bonds[bond.index()].set_stereo(stereo)?;
    // The source property write is represented by the detached result.
    Ok(true)
}

fn neighboring_directed_bond_unchecked(topology: &TopologyBlock, atom: AtomId) -> Option<&Bond> {
    match neighboring_directed_bond_from_incident(
        incident_bonds(topology, atom).map(Ok::<_, std::convert::Infallible>),
    ) {
        Ok(bond) => bond,
        Err(never) => match never {},
    }
}

fn neighbor_directions(
    topology: &TopologyBlock,
    atom: AtomId,
    reference: BondId,
    ranks: &[u32],
    explicit_unknown: &mut bool,
) -> Result<Vec<(AtomId, BondDirection)>, DoubleBondStereoError> {
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
                || property_is_true(bond.prop("_UnknownStereo"), None, Some(bond.id()))?)
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
        return Ok(Vec::new());
    }
    if !has_stereo_bond_direction(result[0].1) {
        result[0].1 = opposite_unchecked(result[1].1);
    } else if result.len() > 1 && !has_stereo_bond_direction(result[1].1) {
        result[1].1 = opposite_unchecked(result[0].1);
    }
    Ok(result)
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
    counts: &[u32],
    double_bond: BondId,
    atom: AtomId,
    double_bond_seen: &mut bool,
) -> Result<Controls, DoubleBondStereoError> {
    // BEGIN RDKIT CPP FUNCTION controllingBondFromAtom
    // RDKit✔️✔️: void controllingBondFromAtom(const ROMol &mol,
    // RDKit✔️✔️:                              const boost::dynamic_bitset<> &needsDir,
    // RDKit✔️✔️:                              const std::vector<unsigned int> &singleBondCounts,
    // RDKit✔️✔️:                              const Bond *dblBond, const Atom *atom, Bond *&bond,
    // RDKit✔️✔️:                              Bond *&obond, bool &squiggleBondSeen,
    // RDKit✔️✔️:                              bool &doubleBondSeen) {
    // RDKit✔️✔️:   bond = nullptr;
    // RDKit✔️✔️:   obond = nullptr;
    // RDKit✔️✔️:   for (const auto tBond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (tBond == dblBond) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if ((tBond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:          tBond->getBondType() == Bond::AROMATIC) &&
    // RDKit✔️✔️:         (tBond->getBondDir() == Bond::BondDir::NONE ||
    // RDKit✔️✔️:          tBond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
    // RDKit✔️✔️:          tBond->getBondDir() == Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:       // prefer bonds that already have their directionality set
    // RDKit✔️✔️:       // or that are adjacent to more double bonds:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         bond = tBond;
    // RDKit✔️✔️:       } else if (needsDir[tBond->getIdx()]) {
    // RDKit✔️✔️:         if (singleBondCounts[tBond->getIdx()] >
    // RDKit✔️✔️:             singleBondCounts[bond->getIdx()]) {
    // RDKit✔️✔️:           obond = bond;
    // RDKit✔️✔️:           bond = tBond;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           obond = tBond;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         obond = bond;
    // RDKit✔️✔️:         bond = tBond;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (tBond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       doubleBondSeen = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     int explicit_unknown_stereo;
    // RDKit✔️✔️:     if ((tBond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:          tBond->getBondType() == Bond::AROMATIC) &&
    // RDKit✔️✔️:         (tBond->getBondDir() == Bond::UNKNOWN ||
    // RDKit✔️✔️:          ((tBond->getPropIfPresent<int>(common_properties::_UnknownStereo,
    // RDKit✔️✔️:                                         explicit_unknown_stereo) &&
    // RDKit✔️✔️:            explicit_unknown_stereo)))) {
    // RDKit✔️✔️:       squiggleBondSeen = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION controllingBondFromAtom
    // Source adjacency order, first candidate, strict greater-count tie rule,
    // and already-directed preference remain exact. The passed double flag
    // accumulates across the caller's two scans; a squiggle stops immediately.
    // Direction UNKNOWN short-circuits the int property getter. The separate
    // detached unknown-stereo metadata bit does not bypass that native getter
    // or hide a wrong present property tag. No rank sorting or fallback.
    // Cost: one source-ordered adjacency pass, indexed counts/direction reads,
    // fixed-key dictionary lookup only on reached SINGLE/AROMATIC guards,
    // constant temporary state, no allocation or molecule/bond clone.
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
        } else if bond.order() == BondOrder::Double {
            *double_bond_seen = true;
        }
        if matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
            && (bond.direction() == BondDirection::Unknown
                || property_is_true(bond.prop("_UnknownStereo"), None, Some(bond.id()))?)
        {
            squiggle = true;
            break;
        }
    }
    Ok(Controls {
        primary,
        secondary,
        squiggle,
    })
}

fn update_double_bond_neighbors(
    topology: &mut TopologyBlock,
    double_bond: BondId,
    conformer: Option<&Conformer3D>,
    needs_direction: &mut [bool],
    counts: &[u32],
    single_bond_neighbors: &[Vec<BondId>],
    needs_detect_bond_stereo: &mut bool,
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
    // Ordinary Int(1) property writes are accumulated as an exact event;
    // callers apply it to their existing carrier after the detached kernel.
    // No molecule property is read by this kernel after those source writes.
    // No completion marker is promoted by this uninstalled proposal.
    // Locally: one source neighbor scan per endpoint and one owned followup
    // list; only scalar bond fields are copied, never a property/query carrier.
    let (double_begin, double_end, double_stereo, double_references) = {
        let value = require_double_bond(topology, double_bond)?;
        (
            value.begin(),
            value.end(),
            value.stereo(),
            value.stereo_atoms(),
        )
    };
    if !needs_direction[double_bond.index()] {
        return Ok(());
    }
    needs_direction[double_bond.index()] = false;
    let mut double_bond_seen = false;
    let mut begin_controls = controlling_bonds(
        topology,
        needs_direction,
        counts,
        double_bond,
        double_begin,
        &mut double_bond_seen,
    )?;
    if begin_controls.squiggle {
        *needs_detect_bond_stereo |=
            set_stereo_for_bond(topology, double_bond, BondStereo::Any, false)?;
        return Ok(());
    }
    let Some(mut begin_bond) = begin_controls.primary else {
        return Ok(());
    };
    let mut end_controls = controlling_bonds(
        topology,
        needs_direction,
        counts,
        double_bond,
        double_end,
        &mut double_bond_seen,
    )?;
    if end_controls.squiggle {
        *needs_detect_bond_stereo |=
            set_stereo_for_bond(topology, double_bond, BondStereo::Any, false)?;
        return Ok(());
    }
    let Some(mut end_bond) = end_controls.primary else {
        return Ok(());
    };
    let mut same_torsion_direction;
    if let Some(conformer) = conformer {
        let coordinates = conformer.coordinates();
        let begin_point = coordinates[double_begin.index()];
        let end_point = coordinates[double_end.index()];
        let mut begin_neighbor = other_atom(&topology.bonds[begin_bond.index()], double_begin);
        let mut end_neighbor = other_atom(&topology.bonds[end_bond.index()], double_end);
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
                begin_neighbor = other_atom(&topology.bonds[begin_bond.index()], double_begin);
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
                    end_neighbor = other_atom(&topology.bonds[end_bond.index()], double_end);
                    end_neighbor_point = coordinates[end_neighbor.index()];
                    linear = is_linear(
                        sub(end_neighbor_point, begin_point),
                        sub(begin_point, end_point),
                    );
                }
            }
        }
        if linear {
            *needs_detect_bond_stereo |=
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
        same_torsion_direction = match double_stereo {
            BondStereo::Cis | BondStereo::Z => false,
            BondStereo::Trans | BondStereo::E => true,
            _ => return Ok(()),
        };
        let references =
            double_references.ok_or(DoubleBondStereoError::StereoReferenceRequired {
                bond: double_bond,
                stereo: double_stereo,
            })?;
        let begin_atom = other_atom(&topology.bonds[begin_bond.index()], double_begin);
        if !references.contains(&begin_atom) {
            same_torsion_direction = !same_torsion_direction;
        }
        let end_atom = other_atom(&topology.bonds[end_bond.index()], double_end);
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
            if topology.bonds[begin_bond.index()].begin() != double_begin {
                reverse = !reverse;
            }
            let direction = topology.bonds[begin_bond.index()].direction();
            set_direction_relative(topology, end_bond, double_end, direction, reverse)?;
        }
        (true, false) => {
            if topology.bonds[end_bond.index()].begin() != double_end {
                reverse = !reverse;
            }
            let direction = topology.bonds[end_bond.index()].direction();
            set_direction_relative(topology, begin_bond, double_begin, direction, reverse)?;
        }
        (true, true) => {
            set_direction_relative(
                topology,
                begin_bond,
                double_begin,
                BondDirection::EndDownRight,
                false,
            )?;
            set_direction_relative(
                topology,
                end_bond,
                double_end,
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
        let reverse = topology.bonds[begin_bond.index()].begin() == double_begin;
        set_direction_relative(topology, secondary, double_begin, direction, reverse)?;
        needs_direction[secondary.index()] = false;
    }
    if let Some(secondary) = end_controls.secondary
        && needs_direction[secondary.index()]
    {
        let direction = topology.bonds[end_bond.index()].direction();
        let reverse = topology.bonds[end_bond.index()].begin() == double_end;
        set_direction_relative(topology, secondary, double_end, direction, reverse)?;
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
            needs_detect_bond_stereo,
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
    // RDKit✔️✔️: double computeDihedralAngle(const Point3D &pt1, const Point3D &pt2,
    // RDKit✔️✔️:                             const Point3D &pt3, const Point3D &pt4) {
    // RDKit✔️✔️:   Point3D begEndVec = pt3 - pt2;
    // RDKit✔️✔️:   Point3D begNbrVec = pt1 - pt2;
    // RDKit✔️✔️:   Point3D crs1 = begNbrVec.crossProduct(begEndVec);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Point3D endNbrVec = pt4 - pt3;
    // RDKit✔️✔️:   Point3D crs2 = endNbrVec.crossProduct(begEndVec);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double ang = crs1.angleTo(crs2);
    // RDKit✔️✔️:   return ang;
    // RDKit✔️✔️: }
    // Canonical geometry owner already implements the exact source cross
    // expression and angleTo ordered roundoff checks. This private adapter
    // forwards fixed-size points once, with no allocation or duplicate math.
    crate::unsigned_dihedral_radians([i, j, k, l])
}

#[cfg(test)]
mod uint_unknown_proposed_tests {
    use super::*;
    #[test]
    fn proposed_uint_unknown_signed_getter_and_native_flag_guard() {
        for (value, truth) in [
            (0_u32, false),
            (1, true),
            (2147483646, true),
            (2147483647, true),
        ] {
            assert_eq!(
                property_is_true(
                    Some(&cosmolkit_model::PropertyValue::UInt(value)),
                    Some(AtomId::new(0)),
                    None
                ),
                Ok(truth)
            );
        }
        for value in [2147483648_u32, 4294967295] {
            assert_eq!(
                property_is_true(
                    Some(&cosmolkit_model::PropertyValue::UInt(value)),
                    Some(AtomId::new(0)),
                    None
                ),
                Err(DoubleBondStereoError::UnsignedPropertyOverflow {
                    atom: Some(AtomId::new(0)),
                    bond: None,
                    property: "_UnknownStereo",
                    value
                })
            );
            let atom = cosmolkit_model::Atom::from_spec(
                AtomId::new(0),
                cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C)
                    .with_unknown_stereo(true)
                    .with_prop(
                        "_UnknownStereo",
                        cosmolkit_model::PropertyValue::UInt(value),
                    )
                    .unwrap(),
            );
            let graph = TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
            assert_eq!(atom_has_unknown_stereo(&graph, AtomId::new(0)), Ok(true));
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, PropertyValue};
    use cosmolkit_types::Element;

    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/double_stereo__UnknownStereo_0
    #[test]
    fn uint_cell_signed_consumer_core_double_stereo__unknownstereo_0() {
        let v = PropertyValue::UInt(0_u32);
        assert_eq!(
            property_is_true(Some(&v), Some(AtomId::new(0)), None),
            Ok(false)
        );
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/double_stereo__UnknownStereo_1
    #[test]
    fn uint_cell_signed_consumer_core_double_stereo__unknownstereo_1() {
        let v = PropertyValue::UInt(1_u32);
        assert_eq!(
            property_is_true(Some(&v), Some(AtomId::new(0)), None),
            Ok(true)
        );
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/double_stereo__UnknownStereo_2147483646
    #[test]
    fn uint_cell_signed_consumer_core_double_stereo__unknownstereo_2147483646() {
        let v = PropertyValue::UInt(2147483646_u32);
        assert_eq!(
            property_is_true(Some(&v), Some(AtomId::new(0)), None),
            Ok(true)
        );
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/double_stereo__UnknownStereo_2147483647
    #[test]
    fn uint_cell_signed_consumer_core_double_stereo__unknownstereo_2147483647() {
        let v = PropertyValue::UInt(2147483647_u32);
        assert_eq!(
            property_is_true(Some(&v), Some(AtomId::new(0)), None),
            Ok(true)
        );
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/double_stereo__UnknownStereo_2147483648
    #[test]
    fn uint_cell_signed_consumer_core_double_stereo__unknownstereo_2147483648() {
        let v = PropertyValue::UInt(2147483648_u32);
        assert_eq!(
            property_is_true(Some(&v), Some(AtomId::new(0)), None),
            Err(DoubleBondStereoError::UnsignedPropertyOverflow {
                atom: Some(AtomId::new(0)),
                bond: None,
                property: "_UnknownStereo",
                value: 2147483648_u32
            })
        );
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/double_stereo__UnknownStereo_4294967295
    #[test]
    fn uint_cell_signed_consumer_core_double_stereo__unknownstereo_4294967295() {
        let v = PropertyValue::UInt(4294967295_u32);
        assert_eq!(
            property_is_true(Some(&v), Some(AtomId::new(0)), None),
            Err(DoubleBondStereoError::UnsignedPropertyOverflow {
                atom: Some(AtomId::new(0)),
                bond: None,
                property: "_UnknownStereo",
                value: 4294967295_u32
            })
        );
    }
}

#[cfg(test)]
mod direction_priority_source_boundaries {
    use super::direction_priority_source_score;

    #[test]
    fn non_ring_priority_wraps_at_the_native_unsigned_boundary() {
        // All source sums here fit signed int: no signed accumulate overflow
        // or implementation-defined Bond-ID conversion is assumed.
        for (score, expected) in [
            (0, 0),
            (429_496_729, 4_294_967_290),
            (429_496_730, 4),
            (2_147_483_647, 4_294_967_286),
        ] {
            assert_eq!(direction_priority_source_score(score, false), expected);
        }
    }

    #[test]
    fn ring_priority_preserves_the_unmultiplied_native_sum() {
        for score in [0, 429_496_729, 429_496_730, 2_147_483_647] {
            assert_eq!(direction_priority_source_score(score, true), score);
        }
    }
}

#[cfg(test)]
mod source_highest_cip_neighbor_complete_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;

    fn fixture(neighbors: &[usize]) -> TopologyBlock {
        let atoms = (0..5)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = neighbors
            .iter()
            .enumerate()
            .map(|(i, &other)| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(0), AtomId::new(other), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    #[test]
    fn source_cip_tie_then_lower_rank_keeps_native_pointer_reset_and_physical_order() {
        let topology = fixture(&[1, 3, 2, 4]);
        let mut reads = Vec::new();
        let mut warnings = Vec::new();
        let result = highest_ranked_neighbor_from_reader(
            &topology,
            AtomId::new(0),
            AtomId::new(1),
            &mut |candidate| {
                reads.push(candidate.index());
                Ok::<_, &'static str>(Some(if candidate.index() == 4 { 1 } else { 9 }))
            },
            &mut |center| warnings.push(center),
        )
        .unwrap();
        assert_eq!(reads, [3, 2, 4]);
        assert_eq!(warnings, [AtomId::new(0)]);
        assert_eq!(
            result,
            Some(AtomId::new(4)),
            "native null best pointer admits a lower following rank"
        );
    }

    #[test]
    fn source_cip_missing_rank_returns_before_later_neighbor_reads() {
        let topology = fixture(&[1, 3, 2, 4]);
        let mut reads = Vec::new();
        let mut warnings = Vec::new();
        let result = highest_ranked_neighbor_from_reader(
            &topology,
            AtomId::new(0),
            AtomId::new(1),
            &mut |candidate| {
                reads.push(candidate.index());
                Ok::<_, &'static str>(if candidate.index() == 2 {
                    None
                } else {
                    Some(9)
                })
            },
            &mut |center| warnings.push(center),
        )
        .unwrap();
        assert_eq!(result, None);
        assert_eq!(reads, [3, 2]);
        assert!(warnings.is_empty());
    }

    #[test]
    fn source_cip_rank_read_failure_propagates_without_later_reads_or_diagnostics() {
        let topology = fixture(&[1, 3, 2, 4]);
        let mut reads = Vec::new();
        let mut warnings = Vec::new();
        let result = highest_ranked_neighbor_from_reader(
            &topology,
            AtomId::new(0),
            AtomId::new(1),
            &mut |candidate| {
                reads.push(candidate.index());
                if candidate.index() == 2 {
                    Err("source rank conversion error")
                } else {
                    Ok(Some(9))
                }
            },
            &mut |center| warnings.push(center),
        );
        assert_eq!(result, Err("source rank conversion error"));
        assert_eq!(reads, [3, 2]);
        assert!(warnings.is_empty());
    }

    #[test]
    fn source_cip_zero_unsigned_max_and_empty_neighbor_boundaries() {
        for (neighbors, ranks, expected) in [
            (vec![1], [0_u32; 5], None),
            (vec![1, 3], [0_u32; 5], Some(AtomId::new(3))),
            (
                vec![1, 3, 2, 4],
                [0, 0, 1, u32::MAX, 2],
                Some(AtomId::new(3)),
            ),
            (Vec::new(), [0_u32; 5], None),
        ] {
            let topology = fixture(&neighbors);
            let mut reads = Vec::new();
            let result = highest_ranked_neighbor_from_reader(
                &topology,
                AtomId::new(0),
                AtomId::new(1),
                &mut |candidate| {
                    reads.push(candidate.index());
                    Ok::<_, &'static str>(Some(ranks[candidate.index()]))
                },
                &mut |_| panic!("none of these source boundary rows contains a rank tie"),
            )
            .unwrap();
            assert_eq!(result, expected);
            assert!(!reads.contains(&1));
        }
    }
}

#[cfg(test)]
mod source_find_stereo_atoms_complete_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;

    fn fixture(stereo: BondStereo, references: bool, fork: bool) -> (TopologyBlock, BondId) {
        let count = if fork { 6 } else { 4 };
        let atoms = (0..count)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let edges = if fork {
            vec![(2, 0), (2, 1), (2, 3), (3, 4), (3, 5)]
        } else {
            vec![(0, 1), (1, 2), (2, 3)]
        };
        let double = if fork { 2 } else { 1 };
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (begin, end))| {
                let mut spec = BondSpec::new(
                    AtomId::new(begin),
                    AtomId::new(end),
                    if i == double {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                );
                if i == double {
                    spec = spec.with_stereo(stereo);
                    if references {
                        spec = spec.with_stereo_atoms(AtomId::new(0), AtomId::new(3));
                    }
                }
                Bond::from_spec(BondId::new(i), spec)
            })
            .collect();
        (
            TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap(),
            BondId::new(double),
        )
    }

    #[test]
    fn stored_references_skip_rank_reads_but_follow_defined_stereo_precondition() {
        for stereo in [
            BondStereo::E,
            BondStereo::Z,
            BondStereo::Cis,
            BondStereo::Trans,
        ] {
            let (topology, bond) = fixture(stereo, true, false);
            let found = find_double_bond_stereo_atoms_with_rank_reader::<DoubleBondStereoError>(
                &topology,
                bond,
                |_| panic!("stored references precede lazy rank reads"),
                |_| panic!("stored references produce no diagnostic"),
            )
            .unwrap();
            assert_eq!(found.atoms, Some([AtomId::new(0), AtomId::new(3)]));
            assert!(found.warnings.is_empty());
        }
        for stereo in [BondStereo::None, BondStereo::Any] {
            let (topology, bond) = fixture(stereo, true, false);
            assert_eq!(
                find_double_bond_stereo_atoms_with_rank_reader::<DoubleBondStereoError>(
                    &topology,
                    bond,
                    |_| panic!("source precondition precedes rank read"),
                    |_| panic!("source precondition precedes diagnostics")
                ),
                Err(DoubleBondStereoError::UndefinedStereo { bond, stereo })
            );
        }
    }

    #[test]
    fn missing_begin_cip_still_evaluates_end_then_returns_empty() {
        for stereo in [BondStereo::E, BondStereo::Z] {
            let (topology, bond) = fixture(stereo, false, false);
            let mut reads = Vec::new();
            let found = find_double_bond_stereo_atoms_with_rank_reader::<DoubleBondStereoError>(
                &topology,
                bond,
                |atom| {
                    reads.push(atom.index());
                    Ok(if atom.index() == 0 { None } else { Some(0) })
                },
                |_| panic!("missing ranks are not tie warnings"),
            )
            .unwrap();
            assert_eq!(reads, [0, 3]);
            assert_eq!(found.atoms, None);
            assert!(found.warnings.is_empty());
        }
    }

    #[test]
    fn specified_non_ez_without_references_warns_and_returns_empty_without_rank_reads() {
        for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
            let (topology, bond) = fixture(stereo, false, false);
            let mut emitted = Vec::new();
            let found = find_double_bond_stereo_atoms_with_rank_reader::<DoubleBondStereoError>(
                &topology,
                bond,
                |_| panic!("source non-E/Z branch has no CIP read"),
                |warning| emitted.push(warning),
            )
            .unwrap();
            assert_eq!(found.atoms, None);
            assert_eq!(emitted, [StereoAtomSearchWarning::UnableToAssign { bond }]);
            assert_eq!(found.warnings, emitted);
            assert_eq!(
                emitted[0].to_string(),
                "Unable to assign stereo atoms for bond 1"
            );
            assert_eq!(
                find_double_bond_stereo_atoms(&topology, bond, &[0; 4]).unwrap(),
                None
            );
        }
    }

    #[derive(Debug, PartialEq)]
    enum ReadError {
        Core(DoubleBondStereoError),
        UInt(crate::PropertyUIntReadError),
    }
    impl From<DoubleBondStereoError> for ReadError {
        fn from(error: DoubleBondStereoError) -> Self {
            Self::Core(error)
        }
    }

    #[test]
    fn begin_tie_diagnostic_survives_later_actual_uint_conversion_error() {
        let (mut topology, bond) = fixture(BondStereo::E, false, true);
        topology.atoms[0]
            .set_prop("_CIPRank", PropertyValue::UInt(9))
            .unwrap();
        topology.atoms[1]
            .set_prop("_CIPRank", PropertyValue::UInt(9))
            .unwrap();
        topology.atoms[4]
            .set_prop("_CIPRank", PropertyValue::String("bad".into()))
            .unwrap();
        let mut reads = Vec::new();
        let mut emitted = Vec::new();
        let result = find_double_bond_stereo_atoms_with_rank_reader(
            &topology,
            bond,
            |atom| {
                reads.push(atom.index());
                topology.atoms[atom.index()]
                    .prop("_CIPRank")
                    .map(crate::property_value_to_uint)
                    .transpose()
                    .map_err(ReadError::UInt)
            },
            |warning| emitted.push(warning),
        );
        assert!(matches!(
            result,
            Err(ReadError::UInt(
                crate::PropertyUIntReadError::Lexical { .. }
            ))
        ));
        assert_eq!(reads, [0, 1, 4]);
        assert_eq!(
            emitted,
            [StereoAtomSearchWarning::DuplicateCipRank {
                center: AtomId::new(2)
            }]
        );
        assert_eq!(
            emitted[0].to_string(),
            "Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"
        );
    }
}

#[cfg(test)]
mod complete_incident_directed_bond_source_tests {
    use super::*;
    use cosmolkit_model::BondSpec;

    fn bond(row: usize, order: BondOrder, direction: BondDirection) -> Bond {
        Bond::from_spec(
            BondId::new(row),
            BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_direction(direction),
        )
    }

    #[test]
    fn exact_source_direction_predicate_skips_double_and_returns_first_physical_match() {
        let bonds = [
            bond(0, BondOrder::Double, BondDirection::EndUpRight),
            bond(1, BondOrder::Single, BondDirection::Unknown),
            bond(2, BondOrder::Single, BondDirection::EndDownRight),
            bond(3, BondOrder::Single, BondDirection::EndUpRight),
        ];
        let found =
            neighboring_directed_bond_from_incident(bonds.iter().map(Ok::<_, &'static str>))
                .unwrap()
                .unwrap();
        assert_eq!(found.id(), BondId::new(2));
        assert!(std::ptr::eq(found, &bonds[2]));
    }

    #[test]
    fn source_accepts_any_non_double_order_and_empty_incident_range_returns_none() {
        let zero = bond(0, BondOrder::Zero, BondDirection::EndUpRight);
        assert!(
            neighboring_directed_bond_from_incident([Ok::<_, &'static str>(&zero)])
                .unwrap()
                .is_some()
        );
        assert!(
            neighboring_directed_bond_from_incident(
                std::iter::empty::<Result<&Bond, &'static str>>()
            )
            .unwrap()
            .is_none()
        );
    }

    #[test]
    fn first_return_and_reached_error_stop_iterator_without_reading_later_bonds() {
        let first = bond(0, BondOrder::Single, BondDirection::EndUpRight);
        assert!(
            neighboring_directed_bond_from_incident([Ok(&first), Err("unreached")])
                .unwrap()
                .is_some()
        );
        assert_eq!(
            neighboring_directed_bond_from_incident([Err("reached"), Ok(&first)]).unwrap_err(),
            "reached"
        );
    }
}

#[cfg(test)]
mod recovery_chem14 {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn fixture(first: BondStereo, weighted: bool, reversed: bool) -> TopologyBlock {
        let elements = [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::CL,
            Element::F,
            Element::C,
        ];
        let ids: [usize; 8] = std::array::from_fn(|i| if reversed { 7 - i } else { i });
        let mut ordered = [Element::C; 8];
        for i in 0..8 {
            ordered[ids[i]] = elements[i];
        }
        let atoms = ordered
            .into_iter()
            .enumerate()
            .map(|(i, e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
            .collect();
        let mut edges = [
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (4, 7, BondOrder::Single),
            (4, 6, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (1, 5, BondOrder::Single),
            (3, 4, BondOrder::Double),
        ];
        if weighted {
            edges.swap(2, 5);
        }
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(ids[a]), AtomId::new(ids[b]), o),
                )
            })
            .collect();
        let mut g = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        for (b, stereo, pair) in [(1, first, [0, 3]), (6, BondStereo::Trans, [2, 7])] {
            g.bonds[b].set_stereo_atoms(Some(pair.map(|a| AtomId::new(ids[a]))));
            g.bonds[b].set_stereo(stereo).unwrap();
        }
        g
    }
    fn run(g: TopologyBlock) -> TopologyBlock {
        let rings = crate::symmetrized_sssr(&g, &crate::RingSearchParams::default()).unwrap();
        let update = set_double_bond_neighbor_directions(g, &rings, None).unwrap();
        assert!(!update.needs_detect_bond_stereo);
        assert!(update.ring_update.is_none());
        update.topology
    }
    #[test]
    fn equal_priorities_use_lower_actual_bond_id_first() {
        // Actual adjacent SINGLE IDs [0,5,4] and [4,2,3] each sum9, score90.
        // Two real connected double bonds; changing first assignment is visible
        // in final physical single directions, not merely in a mock sort trace.
        for reversed in [false, true] {
            let g = run(fixture(BondStereo::Cis, false, reversed));
            assert_eq!(
                [0, 2, 3, 4, 5].map(|i| g.bonds[i].direction()),
                [
                    BondDirection::EndUpRight,
                    BondDirection::EndDownRight,
                    BondDirection::EndUpRight,
                    BondDirection::EndDownRight,
                    BondDirection::EndUpRight
                ]
            );
            assert_eq!(g.bonds[1].stereo(), BondStereo::Cis);
            assert_eq!(g.bonds[6].stereo(), BondStereo::Trans);
        }
    }
    #[test]
    fn priority_is_sum_of_actual_ids_and_descending_before_tie() {
        // Still three single neighbors each: now sums6 and12, scores60/120.
        // Bond6 runs first. Neighbor counts or ascending priority both fail.
        for reversed in [false, true] {
            let g = run(fixture(BondStereo::Cis, true, reversed));
            assert_eq!(
                [0, 2, 3, 4, 5].map(|i| g.bonds[i].direction()),
                [
                    BondDirection::EndDownRight,
                    BondDirection::EndDownRight,
                    BondDirection::EndDownRight,
                    BondDirection::EndUpRight,
                    BondDirection::EndUpRight
                ]
            );
            assert_eq!(g.bonds[1].stereo(), BondStereo::Cis);
            assert_eq!(g.bonds[6].stereo(), BondStereo::Trans);
        }
    }
    #[test]
    fn compatible_trans_control_retains_directions() {
        for weighted in [false, true] {
            let g = run(fixture(BondStereo::Trans, weighted, false));
            let order = if weighted {
                [0, 5, 3, 4, 2]
            } else {
                [0, 2, 3, 4, 5]
            };
            assert_eq!(
                order.map(|i| g.bonds[i].direction()),
                [
                    BondDirection::EndUpRight,
                    BondDirection::EndUpRight,
                    BondDirection::EndDownRight,
                    BondDirection::EndUpRight,
                    BondDirection::EndUpRight
                ]
            );
        }
    }
}
