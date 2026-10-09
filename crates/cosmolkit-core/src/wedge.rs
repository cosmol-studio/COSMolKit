use crate::stereo_graph::StereoGraphMut;
use crate::stereo_graph::{StereoAtomAccess, StereoGraphAccess};
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

use crate::potential_stereo::{self, PotentialStereoError};
use crate::rings::{RingFindType, RingFindingError, RingInfo, RingSearchParams, find_sssr};
use crate::stereo_order::{StereoOrderError, count_swaps_to_interconvert};
use crate::structure_tags::{StereoError, Vec3};
use crate::{
    AtropisomerConformer, AtropisomerDiagnostic, AtropisomerError, AtropisomerWedgeAssignment,
    AtropisomerWedgeUpdate, ValenceAssignment, wedge_bonds_from_atropisomers,
};
use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, CoordinateValidationError, PropertyValue, PropertyValueKind,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

const ATTACHMENT_POINT_PROPERTY: &str = "_fromAttchpt";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WedgeInfo {
    Chiral { center: AtomId },
    Atropisomer { update: AtropisomerWedgeUpdate },
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct WedgeAssignments {
    by_bond: BTreeMap<BondId, WedgeInfo>,
    diagnostics: Vec<AtropisomerDiagnostic>,
}

impl WedgeAssignments {
    pub(crate) fn insert_atropisomer_source(&mut self, update: AtropisomerWedgeUpdate) {
        // RDKit❗✔️:   WedgeInfoBase(int idxInit) : idx(idxInit){};
        // RDKit❗✔️:   WedgeInfoAtropisomer(int bondId, RDKit::Bond::BondDir dirInit)
        // RDKit❗✔️:       : WedgeInfoBase(bondId) {
        // RDKit❗✔️:     dir = dirInit;
        // RDKit❗✔️:   };
        // Field-copy constructors preserve the source axial ID and direction;
        // unsigned BondId/signed native idx iteration is a final transport
        // comparison topic. No domain calculation or heap carrier clone.
        // RDKit✔️✔️:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
        self.by_bond
            .insert(update.bond, WedgeInfo::Atropisomer { update });
    }

    pub(crate) fn append_source_atropisomer_diagnostics(
        &mut self,
        diagnostics: Vec<AtropisomerDiagnostic>,
    ) {
        self.diagnostics.extend(diagnostics);
    }

    /// Returns the source assignment associated with a bond, if one was made.
    #[must_use]
    pub fn get(&self, bond: BondId) -> Option<&WedgeInfo> {
        self.by_bond.get(&bond)
    }

    /// Iterates assignments in source bond-index order.
    pub fn iter(&self) -> impl Iterator<Item = (BondId, &WedgeInfo)> {
        self.by_bond
            .iter()
            .map(|(bond, assignment)| (*bond, assignment))
    }

    /// Returns diagnostics produced while assigning atropisomer wedges.
    #[must_use]
    pub fn diagnostics(&self) -> &[AtropisomerDiagnostic] {
        &self.diagnostics
    }

    /// Converts source-ordered atropisomer updates into the shared bond map.
    #[must_use]
    pub fn from_atropisomer_wedge_assignment(assignment: AtropisomerWedgeAssignment) -> Self {
        let mut result = Self::from_atropisomer_wedge_parts_source(
            &assignment.bond_updates,
            &assignment.source_map_writes,
        );
        result.diagnostics = assignment.diagnostics;
        result
    }

    pub(crate) fn from_atropisomer_wedge_parts_source(
        bond_updates: &[crate::AtropisomerWedgeUpdate],
        source_map_writes: &[BondId],
    ) -> Self {
        // RDKit❗✔️:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
        // Project only actual info writes, retaining the final scalar update
        // for each key. The same projection serves owned and borrowed inputs;
        // graph-only endpoint/direction updates never fabricate an info entry.
        // Two ordered maps match the previous adapter's O(U log U + W log U)
        // transport cost. Values are borrowed and copied only at real inserts;
        // the canonical source info constructor owns the field copies.
        let mut updates: BTreeMap<_, _> = bond_updates
            .iter()
            .map(|update| (update.bond, update))
            .collect();
        let mut result = Self::default();
        for bond in source_map_writes {
            if let Some(update) = updates.remove(bond) {
                result.insert_atropisomer_source(*update);
            }
        }
        result
    }
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum WedgeError {
    #[error("source wedge precondition: {message}")]
    SourcePrecondition { message: &'static str },
    #[error("source wedge {state} index {index} is outside {count} entries")]
    SourceStateIndex {
        state: &'static str,
        index: usize,
        count: usize,
    },
    #[error("source wedge {state} index {index} exceeds unsigned32 transport")]
    SourceIndexWidth { state: &'static str, index: usize },
    #[error("source signed wedge score at atom {atom} overflows during {operation}")]
    SourceSignedScoreOverflow {
        atom: AtomId,
        operation: &'static str,
    },

    #[error(
        "atom {atom} property {property} unsigned value {value} causes positive_overflow converting UInt to signed int"
    )]
    UnsignedRankOverflow {
        atom: AtomId,
        property: &'static str,
        value: u32,
    },
    #[error("invalid topology: {0}")]
    Topology(#[from] TopologyValidationError),
    #[error("invalid coordinates: {0}")]
    Coordinates(#[from] CoordinateValidationError),
    #[error("bond {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("bond {bond} is not single and cannot be wedged")]
    NonSingleBond { bond: BondId },
    #[error("bond {bond} is not incident to wedge center {center}")]
    CenterNotIncident { bond: BondId, center: AtomId },
    #[error("wedge center {center} has unsupported chiral tag {tag:?}")]
    UnsupportedCenterTag { center: AtomId, tag: ChiralTag },
    #[error(transparent)]
    Geometry(#[from] StereoError),
    #[error(transparent)]
    StereoOrder(#[from] StereoOrderError),
    #[error(transparent)]
    PotentialStereo(#[from] PotentialStereoError),
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Atropisomer(#[from] AtropisomerError),
    #[error("atom {atom} property {property} is not an integer stereo rank: {value:?}")]
    InvalidStereoRank {
        atom: AtomId,
        property: &'static str,
        value: cosmolkit_model::PropertyText,
    },
    #[error("atom {atom} property {property} has {kind:?} value, expected an integer stereo rank")]
    InvalidStereoRankType {
        atom: AtomId,
        property: &'static str,
        kind: PropertyValueKind,
    },
}

#[derive(Debug, Clone, Copy)]
pub struct CrossedBondContext<'a> {
    topology: &'a TopologyBlock,
    valence: &'a ValenceAssignment,
    rings: &'a RingInfo,
    use_legacy_stereo_perception: bool,
}

impl<'a> CrossedBondContext<'a> {
    /// Creates a reusable crossed-bond view over one validated detached input.
    pub fn new(
        topology: &'a TopologyBlock,
        valence: &'a ValenceAssignment,
        rings: &'a RingInfo,
        use_legacy_stereo_perception: bool,
    ) -> Result<Self, WedgeError> {
        // Validate once when creating this reusable view; callers probing many
        // bonds should share it instead of repeating the topology walk.
        topology.validate()?;
        Ok(Self {
            topology,
            valence,
            rings,
            use_legacy_stereo_perception,
        })
    }

    fn should_be_crossed_bond(self, bond_id: BondId) -> Result<bool, WedgeError> {
        // RDKit❗✔️: bool shouldBeACrossedBond(const Bond *bond) {
        // RDKit❗✔️:   PRECONDITION(bond, "");
        // RDKit❗✔️:
        // RDKit❗✔️:   // double bond stereochemistry -
        // RDKit❗✔️:   // if the bond isn't specified, then it should go in the mol block
        // RDKit❗✔️:   // as "any", this was sf.net issue 2963522.
        // RDKit❗✔️:   // two caveats to this:
        // RDKit❗✔️:   // 1) if it's a ring bond, we'll only put the "any"
        // RDKit❗✔️:   //    in the mol block if the user specifically asked for it.
        // RDKit❗✔️:   //    Constantly seeing crossed bonds in rings, though maybe
        // RDKit❗✔️:   //    technically correct, is irritating.
        // RDKit❗✔️:   // 2) if it's a terminal bond (where there's no chance of
        // RDKit❗✔️:   //    stereochemistry anyway), we also skip the any.
        // RDKit❗✔️:   //    this was sf.net issue 3009756
        // RDKit❗✔️:
        // RDKit❗✔️:   if (bond->getStereo() == Bond::STEREOANY) {
        // RDKit❗✔️:     // see if any of the neighbors have a wiggle bond - if so, do NOT make this
        // RDKit❗✔️:     // one a cross bond
        // RDKit❗✔️:     for (auto nbrBond : bond->getOwningMol().atomBonds(bond->getBeginAtom())) {
        // RDKit❗✔️:       if (nbrBond->getBondDir() == Bond::UNKNOWN &&
        // RDKit❗✔️:           nbrBond->getBeginAtom()->getIdx() == bond->getBeginAtom()->getIdx()) {
        // RDKit❗✔️:         return false;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     for (auto nbrBond : bond->getOwningMol().atomBonds(bond->getEndAtom())) {
        // RDKit❗✔️:       if (nbrBond->getBondDir() == Bond::UNKNOWN &&
        // RDKit❗✔️:           nbrBond->getBeginAtom()->getIdx() == bond->getEndAtom()->getIdx()) {
        // RDKit❗✔️:         return false;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     return true;  // crossed double bond
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (bond->getStereo() != Bond::BondStereo::STEREONONE) {
        // RDKit❗✔️:     return false;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // if it is in a ring it is not makred as stereo.
        // RDKit❗✔️:   // If either end is terminal, it is not stereo
        // RDKit❗✔️:
        // RDKit❗✔️:   if (!Chirality::detail::isBondPotentialStereoBond(bond)) {
        // RDKit❗✔️:     return false;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // we don't know that it's explicitly unspecified (covered above with
        // RDKit❗✔️:   // the ==STEREOANY check)
        // RDKit❗✔️:
        // RDKit❗✔️:   if (bond->getBondDir() == Bond::EITHERDOUBLE) {
        // RDKit❗✔️:     return true;  // crossed double bond
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   const auto beginAtom = bond->getBeginAtom();
        // RDKit❗✔️:   const auto endAtom = bond->getEndAtom();
        // RDKit❗✔️:   if (beginAtom->getDegree() > 1 && endAtom->getDegree() > 1 &&
        // RDKit❗✔️:       (beginAtom->getTotalValence() - beginAtom->getTotalDegree()) == 1 &&
        // RDKit❗✔️:       (endAtom->getTotalValence() - endAtom->getTotalDegree()) == 1) {
        // RDKit❗✔️:     // we only do this if each atom only has one unsaturation
        // RDKit❗✔️:     // FIX: this is the fix for github #2649, but we will need to
        // RDKit❗✔️:     // change it once we start handling allenes properly
        // RDKit❗✔️:
        // RDKit❗✔️:     if (canBeStereoBond(bond)) {
        // RDKit❗✔️:       return true;  // crossed double bond
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return false;  // NOT crossed double bond
        // RDKit❗✔️: }
        let bond = self
            .topology
            .bonds
            .get(bond_id.index())
            .ok_or(WedgeError::BondOutOfRange {
                bond: bond_id,
                bond_count: self.topology.bonds.len(),
            })?;

        if bond.stereo() == BondStereo::Any {
            for endpoint in [bond.begin(), bond.end()] {
                for neighbor in self.topology.adjacency.neighbors_of(endpoint.index()) {
                    let neighbor_bond = &self.topology.bonds[neighbor.bond.index()];
                    if neighbor_bond.direction() == BondDirection::Unknown
                        && neighbor_bond.begin() == endpoint
                    {
                        return Ok(false);
                    }
                }
            }
            return Ok(true);
        }
        if bond.stereo() != BondStereo::None {
            return Ok(false);
        }

        if !potential_stereo::is_potential_bond_source(
            self.topology,
            self.valence,
            self.rings,
            bond,
        )? {
            return Ok(false);
        }
        if bond.direction() == BondDirection::EitherDouble {
            return Ok(true);
        }

        if self
            .topology
            .adjacency
            .neighbors_of(bond.begin().index())
            .len()
            <= 1
            || self
                .topology
                .adjacency
                .neighbors_of(bond.end().index())
                .len()
                <= 1
        {
            return Ok(false);
        }

        // Each false comparison short-circuits before the next cache read.
        let begin_total_valence = source_total_valence(self.valence, bond.begin(), self.topology)?;
        let begin_total_degree =
            source_u32_total_degree(self.topology, self.valence, bond.begin())?;
        if begin_total_valence.wrapping_sub(begin_total_degree) != 1 {
            return Ok(false);
        }
        let end_total_valence = source_total_valence(self.valence, bond.end(), self.topology)?;
        let end_total_degree = source_u32_total_degree(self.topology, self.valence, bond.end())?;
        if end_total_valence.wrapping_sub(end_total_degree) != 1 {
            return Ok(false);
        }
        if can_be_stereo_bond(self.topology, bond, self.use_legacy_stereo_perception)? {
            return Ok(true);
        }
        Ok(false)
        // Source degree/cache helpers are shared with their canonical owners.
        // Neighbor scans/rank vector are source-shaped, with no ring-vector
        // allocation or speculative endpoint cache read.
    }
}

/// Molfile-style direction and endpoint-order information for one bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolFileBondStereoInfo {
    pub direction: BondDirection,
    pub direction_code: i32,
    pub reverse: bool,
}

/// Reproduces RDKit's `GetMolFileBondStereoInfo` for a detached bond using a
/// shared wedge assignment and a reusable crossed-bond context.
pub fn get_molfile_bond_stereo_info(
    crossed_bonds: &CrossedBondContext<'_>,
    wedge_assignments: &WedgeAssignments,
    bond_id: BondId,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<MolFileBondStereoInfo, WedgeError> {
    let mut direction_code = 0;
    let mut reverse = false;
    let direction = get_molfile_bond_stereo_code_source(
        crossed_bonds,
        wedge_assignments,
        bond_id,
        conformer,
        &mut direction_code,
        &mut reverse,
    )?;
    Ok(MolFileBondStereoInfo {
        direction,
        direction_code,
        reverse,
    })
}

fn get_molfile_bond_stereo_code_source(
    crossed_bonds: &CrossedBondContext<'_>,
    wedge_assignments: &WedgeAssignments,
    bond_id: BondId,
    conformer: Option<AtropisomerConformer<'_>>,
    direction_code: &mut i32,
    reverse: &mut bool,
) -> Result<BondDirection, WedgeError> {
    // RDKit❗✔️: void GetMolFileBondStereoInfo(
    // RDKit❗✔️:     const Bond *bond,
    // RDKit❗✔️:     const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗✔️:         &wedgeBonds,
    // RDKit❗✔️:     const Conformer *conf, int &dirCode, bool &reverse) {
    // RDKit❗✔️:   Bond::BondDir dir;
    // RDKit❗✔️:   GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dir, reverse);
    // RDKit❗✔️:   dirCode = BondGetDirCode(dir);
    // RDKit❗✔️: }
    // The native local is assigned by the successful child before use. Rust's
    // enum initializer never escapes on error and does not supply a missing
    // source result. Code assignment occurs only after the direction child
    // succeeds; reverse preserves that child's actual mutation/error prefix.
    // Return the same scalar direction for the existing detached value API,
    // without a second direction algorithm or child call. O(1) adapter work
    // and stack storage plus the sole source child/finite code switch.
    let mut direction = BondDirection::None;
    get_molfile_bond_stereo_direction_source(
        crossed_bonds,
        wedge_assignments,
        bond_id,
        conformer,
        &mut direction,
        reverse,
    )?;
    *direction_code = bond_get_dir_code(direction);
    Ok(direction)
}

fn get_molfile_bond_stereo_direction_source(
    crossed_bonds: &CrossedBondContext<'_>,
    wedge_assignments: &WedgeAssignments,
    bond_id: BondId,
    conformer: Option<AtropisomerConformer<'_>>,
    direction: &mut BondDirection,
    reverse: &mut bool,
) -> Result<(), WedgeError> {
    // RDKit❗✔️: void GetMolFileBondStereoInfo(
    // RDKit❗✔️:     const Bond *bond,
    // RDKit❗✔️:     const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗✔️:         &wedgeBonds,
    // RDKit❗✔️:     const Conformer *conf, Bond::BondDir &dir, bool &reverse) {
    // RDKit❗✔️:   PRECONDITION(bond, "");
    // RDKit❗✔️:   reverse = false;
    // RDKit❗✔️:   dir = Bond::NONE;
    // RDKit❗✔️:   if (canHaveDirection(*bond)) {
    // RDKit❗✔️:     // single bond stereo chemistry
    // RDKit❗✔️:
    // RDKit❗✔️:     dir = Chirality::detail::determineBondWedgeState(bond, wedgeBonds, conf);
    // RDKit❗✔️:
    // RDKit❗✔️:     // if this bond needs to be wedged it is possible that this
    // RDKit❗✔️:     // wedging was determined by a chiral atom at the end of the
    // RDKit❗✔️:     // bond (instead of at the beginning). In this case we need to
    // RDKit❗✔️:     // reverse the begin and end atoms for the bond when we write
    // RDKit❗✔️:     // the mol file
    // RDKit❗✔️:     if ((dir == Bond::BEGINDASH) ||
    // RDKit❗✔️:         (dir == Bond::BEGINWEDGE || dir == Bond::UNKNOWN)) {
    // RDKit❗✔️:       auto wbi = wedgeBonds.find(bond->getIdx());
    // RDKit❗✔️:       if (wbi != wedgeBonds.end() &&
    // RDKit❗✔️:           wbi->second->getType() ==
    // RDKit❗✔️:               Chirality::WedgeInfoType::WedgeInfoTypeChiral &&
    // RDKit❗✔️:           static_cast<unsigned int>(wbi->second->getIdx()) !=
    // RDKit❗✔️:               bond->getBeginAtomIdx()) {
    // RDKit❗✔️:         reverse = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit❗✔️:     if (Chirality::shouldBeACrossedBond(bond)) {
    // RDKit❗✔️:       dir = Bond::BondDir::EITHERDOUBLE;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit❗✔️:   auto bondType = bond.getBondType();
    // RDKit❗✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit❗✔️: }
    // PRECONDITION precedes both output resets. Later child errors expose
    // the already-reset output prefix, matching native reference parameters.
    // O(log W) lookups plus selected source child, O(1) local storage. No
    // graph/assignment/conformer clone or eager unsupported-branch work.
    let topology = crossed_bonds.topology;
    let bond = topology
        .bonds
        .get(bond_id.index())
        .ok_or(WedgeError::BondOutOfRange {
            bond: bond_id,
            bond_count: topology.bonds.len(),
        })?;
    *reverse = false;
    *direction = BondDirection::None;
    if matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic) {
        get_directional_bond_stereo_source(
            topology,
            wedge_assignments,
            bond_id,
            conformer,
            direction,
            reverse,
        )?;
    } else if bond.order() == BondOrder::Double {
        if crossed_bonds.should_be_crossed_bond(bond_id)? {
            *direction = BondDirection::EitherDouble;
        }
    }
    Ok(())
}

fn determine_bond_wedge_state_from_assignments<G: StereoGraphAccess>(
    topology: &G,
    bond_id: BondId,
    wedge_assignments: &WedgeAssignments,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<BondDirection, WedgeError> {
    // RDKit❗✔️: Bond::BondDir determineBondWedgeState(
    // RDKit❗✔️:     const Bond *bond,
    // RDKit❗✔️:     const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗✔️:         &wedgeBonds,
    // RDKit❗✔️:     const Conformer *conf) {
    // RDKit❗✔️:   PRECONDITION(bond, "no bond");
    // RDKit❗✔️:   int bid = bond->getIdx();
    // RDKit❗✔️:   auto wbi = wedgeBonds.find(bid);
    // RDKit❗✔️:   if (wbi == wedgeBonds.end()) {
    // RDKit❗✔️:     return bond->getBondDir();
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (wbi->second->getType() ==
    // RDKit❗✔️:       Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer) {
    // RDKit❗✔️:     return wbi->second->getDir();
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return determineBondWedgeState(bond, wbi->second->getIdx(), conf);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:   int getIdx() const { return idx; }
    // RDKit❗✔️:   Bond::BondDir getDir() const override { return dir; }
    // RDKit❗✔️:   WedgeInfoType getType() const override {
    // RDKit❗✔️:     return Chirality::WedgeInfoType::WedgeInfoTypeChiral;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   WedgeInfoType getType() const override {
    // RDKit❗✔️:     return Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int getIdx() const { return d_index; }
    // RDKit❗✔️:   BondDir getBondDir() const { return static_cast<BondDir>(d_dirTag); }
    // One map find by actual object ID, then only the selected source branch.
    // No bond-order/tag/endpoint/conformer validation in missing-map or Atrop
    // branches. Chiral entries delegate to the sole full coordinate kernel.
    // O(log W) lookup plus selected child cost, O(1) local storage, no clone.
    let bond = topology
        .bonds()
        .get(bond_id.index())
        .ok_or(WedgeError::BondOutOfRange {
            bond: bond_id,
            bond_count: topology.bonds().len(),
        })?;
    match wedge_assignments.by_bond.get(&bond.id()) {
        None => Ok(bond.direction()),
        Some(WedgeInfo::Atropisomer { update }) => Ok(update.direction),
        Some(WedgeInfo::Chiral { center }) => {
            determine_bond_wedge_state_graph(topology, bond_id, *center, conformer)
        }
    }
}

fn bond_get_dir_code(direction: BondDirection) -> i32 {
    // BEGIN RDKIT CPP FUNCTION BondGetDirCode
    // RDKit✔️✔️: int BondGetDirCode(const Bond::BondDir dir) {
    // RDKit✔️✔️:   int res = 0;
    // RDKit✔️✔️:   switch (dir) {
    // RDKit✔️✔️:     case Bond::NONE:
    // RDKit✔️✔️:       res = 0;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::BEGINWEDGE:
    // RDKit✔️✔️:       res = 1;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::BEGINDASH:
    // RDKit✔️✔️:       res = 6;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::UNKNOWN:
    // RDKit✔️✔️:       res = 4;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::BondDir::EITHERDOUBLE:
    // RDKit✔️✔️:       res = 3;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BondGetDirCode
    match direction {
        BondDirection::None => 0,
        BondDirection::BeginWedge => 1,
        BondDirection::BeginDash => 6,
        BondDirection::Unknown => 4,
        BondDirection::EitherDouble => 3,
        BondDirection::EndDownRight | BondDirection::EndUpRight => 0,
    }
    // This is a finite enum dispatch with no allocation, matching the source
    // switch's constant-time behavior.
}

fn source_total_valence(
    valence: &ValenceAssignment,
    atom: AtomId,
    topology: &TopologyBlock,
) -> Result<u32, WedgeError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getTotalValence
    // RDKit❗✔️: unsigned int Atom::getTotalValence() const {
    // RDKit❗✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalValence
    // RDKit❗✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit❗✔️:   if (!dp_mol) {
    // RDKit❗✔️:     return 0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit❗✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit❗✔️:        d_implicitValence > -1),
    // RDKit❗✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit❗✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit❗✔️:     return d_explicitValence;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let explicit = valence.explicit_valence.get(atom.index()).copied().ok_or(
        PotentialStereoError::InvalidValence {
            field: "explicit_valence",
            actual: valence.explicit_valence.len(),
            atom_count: topology.atoms.len(),
        },
    )?;
    let explicit =
        u32::try_from(explicit).map_err(|_| PotentialStereoError::InvalidValenceValue {
            field: "explicit_valence",
            atom,
            value: explicit,
        })?;
    let source_atom = topology
        .atoms
        .get(atom.index())
        .ok_or(WedgeError::SourceStateIndex {
            state: "atom",
            index: atom.index(),
            count: topology.atoms.len(),
        })?;
    let implicit = crate::hcount::implicit_hydrogen_count(source_atom, valence)
        .map_err(PotentialStereoError::from)?;
    Ok(explicit.wrapping_add(implicit))
    // Behavior remains provisional until Step 28 raw-valence cases run.
    // The two indexed reads and wrapping addition are O(1), matching the
    // source getters; no detached state is cloned.
}

fn source_u32_total_degree(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
) -> Result<u32, WedgeError> {
    let degree = potential_stereo::total_degree(topology, valence, atom)?;
    u32::try_from(degree)
        .map_err(|_| PotentialStereoError::InvalidAtomDegree { atom, degree }.into())
}

fn can_be_stereo_bond(
    topology: &TopologyBlock,
    bond: &Bond,
    use_legacy_stereo_perception: bool,
) -> Result<bool, WedgeError> {
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
    // RDKit❗✔️: bool canBeStereoBond(const Bond *bond) {
    // RDKit❗✔️:   PRECONDITION(bond, "no bond");
    // RDKit❗✔️:   if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit❗✔️:       bond->getBondType() != Bond::BondType::AROMATIC) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto beginAtom = bond->getBeginAtom();
    // RDKit❗✔️:   auto endAtom = bond->getEndAtom();
    // RDKit❗✔️:   for (const auto atom : {beginAtom, endAtom}) {
    // RDKit❗✔️:     std::vector<int> nbrRanks;
    // RDKit❗✔️:     for (auto nbrBond : bond->getOwningMol().atomBonds(atom)) {
    // RDKit❗✔️:       if (nbrBond == bond) {
    // RDKit❗✔️:         continue;  // a bond is NOT its own neighbor
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       if (nbrBond->getBondType() == Bond::SINGLE) {
    // RDKit❗✔️:         // if a neighbor has a wedge or hash bond, do NOT mark it as double
    // RDKit❗✔️:         // crossed
    // RDKit❗✔️:         if (nbrBond->getBondDir() == Bond::ENDUPRIGHT ||
    // RDKit❗✔️:             nbrBond->getBondDir() == Bond::ENDDOWNRIGHT) {
    // RDKit❗✔️:           return false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         // if a neighbor has a wiggle bond, do NOT mark it as crossed (although
    // RDKit❗✔️:         // it is unknown
    // RDKit❗✔️:         if (nbrBond->getBondDir() == Bond::BondDir::UNKNOWN &&
    // RDKit❗✔️:             nbrBond->getBeginAtom() == atom) {
    // RDKit❗✔️:           return false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         // if two neighbors havr the same CIP ranking, this is not stereo
    // RDKit❗✔️:         const auto otherAtom = nbrBond->getOtherAtom(atom);
    // RDKit❗✔️:         int rank;
    // RDKit❗✔️:         if (RDKit::Chirality::getUseLegacyStereoPerception()) {
    // RDKit❗✔️:           if (!otherAtom->getPropIfPresent(common_properties::_CIPRank, rank)) {
    // RDKit❗✔️:             rank = -1;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else {  // NOT legacy stereo
    // RDKit❗✔️:           if (!otherAtom->getPropIfPresent(common_properties::_ChiralAtomRank,
    // RDKit❗✔️:                                            rank)) {
    // RDKit❗✔️:             rank = -1;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (rank >= 0) {
    // RDKit❗✔️:           if (std::find(nbrRanks.begin(), nbrRanks.end(), rank) !=
    // RDKit❗✔️:               nbrRanks.end()) {
    // RDKit❗✔️:             return false;
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             nbrRanks.push_back(rank);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    if bond.order() != BondOrder::Double && bond.order() != BondOrder::Aromatic {
        return Ok(false);
    }

    let rank_property = if use_legacy_stereo_perception {
        "_CIPRank"
    } else {
        "_ChiralAtomRank"
    };
    for atom in [bond.begin(), bond.end()] {
        let mut neighbor_ranks = Vec::new();
        for neighbor in topology.adjacency.neighbors_of(atom.index()) {
            if neighbor.bond == bond.id() {
                continue;
            }
            let neighbor_bond = &topology.bonds[neighbor.bond.index()];
            if neighbor_bond.order() != BondOrder::Single {
                continue;
            }
            if matches!(
                neighbor_bond.direction(),
                BondDirection::EndUpRight | BondDirection::EndDownRight
            ) {
                return Ok(false);
            }
            if neighbor_bond.direction() == BondDirection::Unknown && neighbor_bond.begin() == atom
            {
                return Ok(false);
            }
            let neighbor_atom = &topology.atoms[neighbor.atom_index];
            let Some(rank_text) = neighbor_atom.prop(rank_property) else {
                continue;
            };
            let rank = match rank_text {
                PropertyValue::Int(rank) => *rank,
                PropertyValue::UInt(rank) => {
                    i32::try_from(*rank).map_err(|_| WedgeError::UnsignedRankOverflow {
                        atom: neighbor_atom.id(),
                        property: rank_property,
                        value: *rank,
                    })?
                }
                value @ PropertyValue::String(rank_text) => crate::property_value_to_int(value)
                    .map_err(|_| WedgeError::InvalidStereoRank {
                        atom: neighbor_atom.id(),
                        property: rank_property,
                        value: rank_text.to_owned(),
                    })?,
                PropertyValue::StringVector(_)
                | PropertyValue::IntVector(_)
                | PropertyValue::Double(_)
                | PropertyValue::Bool(_) => {
                    return Err(WedgeError::InvalidStereoRankType {
                        atom: neighbor_atom.id(),
                        property: rank_property,
                        kind: rank_text.kind(),
                    });
                }
            };
            if rank >= 0 {
                if neighbor_ranks.contains(&rank) {
                    return Ok(false);
                }
                neighbor_ranks.push(rank);
            }
        }
    }
    Ok(true)
    // Behavior remains provisional until Step 28 rank/profile cases run.
    // Complexity review: endpoint scans and rank-vector comparisons match the
    // source O(degree^2) shape, but parsing the model's string-backed rank
    // properties adds work compared with RDKit's typed property lookup.
}

pub fn determine_bond_wedge_state(
    topology: &TopologyBlock,
    bond_id: BondId,
    from_atom: AtomId,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<BondDirection, WedgeError> {
    determine_bond_wedge_state_graph(topology, bond_id, from_atom, conformer)
}

fn determine_bond_wedge_state_graph<G: StereoGraphAccess>(
    topology: &G,
    bond_id: BondId,
    from_atom: AtomId,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<BondDirection, WedgeError> {
    // BEGIN RDKIT CPP FUNCTION detail::determineBondWedgeState
    // RDKit❗✔️: Bond::BondDir determineBondWedgeState(const Bond *bond,
    // RDKit❗✔️:                                       unsigned int fromAtomIdx,
    // RDKit❗✔️:                                       const Conformer *conf) {
    // RDKit❗✔️:   PRECONDITION(bond, "no bond");
    // RDKit❗✔️:   PRECONDITION(bond->getBondType() == Bond::SINGLE,
    // RDKit❗✔️:                "bad bond order for wedging");
    // RDKit❗✔️:   const auto mol = &(bond->getOwningMol());
    // RDKit❗✔️:   PRECONDITION(mol, "no mol");
    // RDKit❗✔️:
    // RDKit❗✔️:   auto res = bond->getBondDir();
    // RDKit❗✔️:   if (!conf) {
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   Atom *atom;
    // RDKit❗✔️:   Atom *bondAtom;
    // RDKit❗✔️:   if (bond->getBeginAtom()->getIdx() == fromAtomIdx) {
    // RDKit❗✔️:     atom = bond->getBeginAtom();
    // RDKit❗✔️:     bondAtom = bond->getEndAtom();
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     atom = bond->getEndAtom();
    // RDKit❗✔️:     bondAtom = bond->getBeginAtom();
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   auto chiralType = atom->getChiralTag();
    // RDKit❗✔️:   TEST_ASSERT(chiralType == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:               chiralType == Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we got this far, we really need to think about it:
    // RDKit❗✔️:   std::list<int> neighborBondIndices;
    // RDKit❗✔️:   std::list<double> neighborBondAngles;
    // RDKit❗✔️:   auto centerLoc = conf->getAtomPos(atom->getIdx());
    // RDKit❗✔️:   auto tmpPt = conf->getAtomPos(bondAtom->getIdx());
    // RDKit❗✔️:   centerLoc.z = 0.0;
    // RDKit❗✔️:   tmpPt.z = 0.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   RDGeom::Point3D refVect;
    // RDKit❗✔️:   try {
    // RDKit❗✔️:     refVect = centerLoc.directionVector(tmpPt);
    // RDKit❗✔️:   } catch (const std::runtime_error &) {
    // RDKit❗✔️:     // we have a problem with the reference bond;
    // RDKit❗✔️:     // it's probably that the center and the tmp atom overlap
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   neighborBondIndices.push_back(bond->getIdx());
    // RDKit❗✔️:   neighborBondAngles.push_back(0.0);
    // RDKit❗✔️:   for (const auto nbrBond : mol->atomBonds(atom)) {
    // RDKit❗✔️:     const auto otherAtom = nbrBond->getOtherAtom(atom);
    // RDKit❗✔️:     if (nbrBond != bond) {
    // RDKit❗✔️:       tmpPt = conf->getAtomPos(otherAtom->getIdx());
    // RDKit❗✔️:       tmpPt.z = 0.0;
    // RDKit❗✔️:       RDGeom::Point3D tmpVect;
    // RDKit❗✔️:       try {
    // RDKit❗✔️:         tmpVect = centerLoc.directionVector(tmpPt);
    // RDKit❗✔️:       } catch (const std::runtime_error &) {
    // RDKit❗✔️:         // we have a problem with the tmp bond;
    // RDKit❗✔️:         // it's probably that the atoms overlap
    // RDKit❗✔️:         return res;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       auto angle = refVect.signedAngleTo(tmpVect);
    // RDKit❗✔️:       if (angle < 0.0) {
    // RDKit❗✔️:         angle += 2. * M_PI;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       auto nbrIt = neighborBondIndices.begin();
    // RDKit❗✔️:       auto angleIt = neighborBondAngles.begin();
    // RDKit❗✔️:       // find the location of this neighbor in our angle-sorted list
    // RDKit❗✔️:       // of neighbors:
    // RDKit❗✔️:       while (angleIt != neighborBondAngles.end() && angle > (*angleIt)) {
    // RDKit❗✔️:         ++angleIt;
    // RDKit❗✔️:         ++nbrIt;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       neighborBondAngles.insert(angleIt, angle);
    // RDKit❗✔️:       neighborBondIndices.insert(nbrIt, nbrBond->getIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // at this point, neighborBondIndices contains a list of bond
    // RDKit❗✔️:   // indices from the central atom.  They are arranged starting
    // RDKit❗✔️:   // at the reference bond in CCW order (based on the current
    // RDKit❗✔️:   // depiction).
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we already have one bond with direction set, then we can use it to
    // RDKit❗✔️:   // decide what the direction of this one is
    // RDKit❗✔️:
    // RDKit❗✔️:   // we're starting from scratch... do the work!
    // RDKit❗✔️:   int nSwaps = atom->getPerturbationOrder(neighborBondIndices);
    // RDKit❗✔️:
    // RDKit❗✔️:   // in the case of three-coordinated atoms we may have to worry about
    // RDKit❗✔️:   // the location of the implicit hydrogen - Issue 209
    // RDKit❗✔️:   // Check if we have one of these situation
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //      0        1 0 2
    // RDKit❗✔️:   //      *         \*/
    // RDKit❗✔️:   //  1 - C - 2      C
    // RDKit❗✔️:   //
    // RDKit❗✔️:   // here the hydrogen will be between 1 and 2 and we need to add an
    // RDKit❗✔️:   // additional swap
    // RDKit❗✔️:   if (neighborBondAngles.size() == 3) {
    // RDKit❗✔️:     // three coordinated
    // RDKit❗✔️:     auto angleIt = neighborBondAngles.begin();
    // RDKit❗✔️:     ++angleIt;  // the first is the 0 (or reference bond - we will ignore
    // RDKit❗✔️:                 // that
    // RDKit❗✔️:     double angle1 = (*angleIt);
    // RDKit❗✔️:     ++angleIt;
    // RDKit❗✔️:     double angle2 = (*angleIt);
    // RDKit❗✔️:     constexpr double angleTol =
    // RDKit❗✔️:         M_PI * 1.9 / 180.;  // just under 2 degrees tolerance, which is what we
    // RDKit❗✔️:                             // use when perceiving T-shaped geometries
    // RDKit❗✔️:     if (angle2 - angle1 >= (M_PI - angleTol)) {
    // RDKit❗✔️:       // we have the above situation
    // RDKit❗✔️:       nSwaps++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️: #ifdef VERBOSE_STEREOCHEM
    // RDKit❗✔️:   BOOST_LOG(rdDebugLog) << "--------- " << nSwaps << std::endl;
    // RDKit❗✔️:   std::copy(neighborBondIndices.begin(), neighborBondIndices.end(),
    // RDKit❗✔️:             std::ostream_iterator<int>(BOOST_LOG(rdDebugLog), " "));
    // RDKit❗✔️:   BOOST_LOG(rdDebugLog) << std::endl;
    // RDKit❗✔️:   std::copy(neighborBondAngles.begin(), neighborBondAngles.end(),
    // RDKit❗✔️:             std::ostream_iterator<double>(BOOST_LOG(rdDebugLog), " "));
    // RDKit❗✔️:   BOOST_LOG(rdDebugLog) << std::endl;
    // RDKit❗✔️: #endif
    // RDKit❗✔️:   if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:     if (nSwaps % 2 == 1) {
    // RDKit❗✔️:       res = Bond::BEGINDASH;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       res = Bond::BEGINWEDGE;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     if (nSwaps % 2 == 1) {
    // RDKit❗✔️:       res = Bond::BEGINWEDGE;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       res = Bond::BEGINDASH;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION detail::determineBondWedgeState
    // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder
    // RDKit❗✔️: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
    // RDKit❗✔️:   INT_LIST ref;
    // RDKit❗✔️:   for (const auto bnd : getOwningMol().atomBonds(this)) {
    // RDKit❗✔️:     ref.push_back(bnd->getIdx());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Atom::getPerturbationOrder
    // BEGIN RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit❗✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit❗✔️:   typename T::const_iterator refIt = ref.begin();
    // RDKit❗✔️:   typename T::iterator probeIt = probe.begin();
    // RDKit❗✔️:   typename T::iterator probeIt2;
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int nSwaps = 0;
    // RDKit❗✔️:   while (refIt != ref.end()) {
    // RDKit❗✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit❗✔️:       bool foundIt = false;
    // RDKit❗✔️:       probeIt2 = probeIt;
    // RDKit❗✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit❗✔️:         ++probeIt2;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (probeIt2 != probe.end()) {
    // RDKit❗✔️:         foundIt = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit❗✔️:
    // RDKit❗✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit❗✔️:       nSwaps++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     ++probeIt;
    // RDKit❗✔️:     ++refIt;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nSwaps;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDGeneral::countSwapsToInterconvert
    // RDKit❗✔️: Atom *Bond::getBeginAtom() const {
    // RDKit❗✔️:   PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    // RDKit❗✔️:   return dp_mol->getAtomWithIdx(d_beginAtomIdx);
    // RDKit❗✔️: };
    // RDKit❗✔️: Atom *Bond::getEndAtom() const {
    // RDKit❗✔️:   PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    // RDKit❗✔️:   return dp_mol->getAtomWithIdx(d_endAtomIdx);
    // RDKit❗✔️: };
    // RDKit❗✔️: Atom *Bond::getOtherAtom(Atom const *what) const {
    // RDKit❗✔️:   PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    // RDKit❗✔️:
    // RDKit❗✔️:   return dp_mol->getAtomWithIdx(getOtherAtomIdx(what->getIdx()));
    // RDKit❗✔️: };
    // RDKit❗✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit❗✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_endAtomIdx;
    // RDKit❗✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_beginAtomIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit❗✔️:   // above, but no need to be redundant, so just throw.
    // RDKit❗✔️:   POSTCONDITION(false, "bad index");
    // RDKit❗✔️: }
    // RDKit❗✔️: const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
    // RDKit❗✔️:   if (dp_mol) {
    // RDKit❗✔️:     PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   URANGE_CHECK(atomId, d_positions.size());
    // RDKit❗✔️:   return d_positions.at(atomId);
    // RDKit❗✔️: }
    // RDKit❗✔️: Atom *ROMol::getAtomWithIdx(unsigned int idx) {
    // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit❗✔️:   auto res = d_graph[vd];
    // RDKit❗✔️:   POSTCONDITION(res, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️:   BondType getBondType() const { return static_cast<BondType>(d_bondType); }
    // RDKit❗✔️:   unsigned int getIdx() const { return d_index; }
    // RDKit❗✔️:   BondDir getBondDir() const { return static_cast<BondDir>(d_dirTag); }
    // RDKit❗✔️:   ROMol &getOwningMol() const {
    // RDKit❗✔️:     PRECONDITION(dp_mol, "no owner");
    // RDKit❗✔️:     return *dp_mol;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int getIdx() const { return d_index; }
    // RDKit❗✔️:   ChiralType getChiralTag() const {
    // RDKit❗✔️:     return static_cast<ChiralType>(d_chiralTag);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ROMol &getOwningMol() const {
    // RDKit❗✔️:     PRECONDITION(dp_mol, "no owner");
    // RDKit❗✔️:     return *dp_mol;
    // RDKit❗✔️:   }
    // RDKit❗✔️: ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
    // RDKit❗✔️:   PRECONDITION(at, "no atom");
    // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit❗✔️:                "atom not associated with this molecule");
    // RDKit❗✔️:   return boost::out_edges(at->getIdx(), d_graph);
    // RDKit❗✔️: }
    // RDKit❗✔️:   CXXBondIterator<MolGraph, Bond *, MolGraph::out_edge_iterator> atomBonds(
    // RDKit❗✔️:       Atom const *at) {
    // RDKit❗✔️:     auto pr = getAtomBonds(at);
    // RDKit❗✔️:     return {&d_graph, pr.first, pr.second};
    // RDKit❗✔️:   }
    // Source-local access order is intentional: no whole graph/conformer
    // scan precedes the null-conformer return or the overlap return. Detached
    // conformers model native conformers without a molecule owner; positions
    // are range checked only when getAtomPos is reached. The source owner's
    // optional cardinality precondition requires explicit provenance, never a
    // guessed owner. Native pointer identity is represented by bond row index.
    // Cost: O(d^2) insertion/permutation work and O(d) scratch as in the source;
    // contiguous pairs replace two lists without an additional graph scan.
    let bond = topology
        .bonds()
        .get(bond_id.index())
        .ok_or(WedgeError::BondOutOfRange {
            bond: bond_id,
            bond_count: topology.bonds().len(),
        })?;
    if bond.order() != BondOrder::Single {
        return Err(WedgeError::NonSingleBond { bond: bond_id });
    }
    let mut result = bond.direction();
    let Some(conformer) = conformer else {
        return Ok(result);
    };
    let atom = |id: AtomId| {
        topology
            .atoms()
            .get(id.index())
            .ok_or(WedgeError::SourceStateIndex {
                state: "atom",
                index: id.index(),
                count: topology.atoms().len(),
            })
    };
    let begin = atom(bond.begin())?;
    let (center_atom, bond_atom) = if begin.id() == from_atom {
        (begin, atom(bond.end())?)
    } else {
        (atom(bond.end())?, begin)
    };
    let center = center_atom.id();
    let chiral_tag = center_atom.chiral_tag();
    if !matches!(
        chiral_tag,
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    ) {
        return Err(WedgeError::UnsupportedCenterTag {
            center,
            tag: chiral_tag,
        });
    }
    // Native always zeros z, including NaN/Inf z, independent of is3D.
    let coordinate = |id: AtomId| -> Result<[f64; 3], WedgeError> {
        let (position, count) = match conformer {
            AtropisomerConformer::TwoD(conf) => (
                conf.coordinates()
                    .get(id.index())
                    .map(|p| [p[0], p[1], 0.0]),
                conf.coordinates().len(),
            ),
            AtropisomerConformer::ThreeD(conf) => (
                conf.coordinates()
                    .get(id.index())
                    .map(|p| [p[0], p[1], 0.0]),
                conf.coordinates().len(),
            ),
        };
        position.ok_or(WedgeError::SourceStateIndex {
            state: "conformer position",
            index: id.index(),
            count,
        })
    };
    let center_location = coordinate(center)?;
    let reference_point = coordinate(bond_atom.id())?;
    let reference_vector =
        match Vec3::normalized_between(center_location, reference_point, center, bond_atom.id()) {
            Ok(vector) => vector,
            Err(StereoError::ZeroLengthVector { .. }) => return Ok(result),
            Err(error) => return Err(WedgeError::Geometry(error)),
        };
    let mut neighbor_angles = vec![(bond.id(), 0.0_f64)];
    let neighbors = topology
        .adjacency()
        .try_neighbors_of(center.index())
        .ok_or(WedgeError::SourcePrecondition {
            message: "source wedge atom adjacency row is absent",
        })?;
    for neighbor in neighbors.iter() {
        let neighbor_bond =
            topology
                .bonds()
                .get(neighbor.bond.index())
                .ok_or(WedgeError::BondOutOfRange {
                    bond: neighbor.bond,
                    bond_count: topology.bonds().len(),
                })?;
        // getOtherAtom is evaluated even for the reference pointer. Use actual
        // bond endpoints/object IDs, never the cached neighbor atom index.
        let other_id = if neighbor_bond.begin() == center {
            neighbor_bond.end()
        } else if neighbor_bond.end() == center {
            neighbor_bond.begin()
        } else {
            return Err(WedgeError::CenterNotIncident {
                bond: neighbor.bond,
                center,
            });
        };
        let other_atom = atom(other_id)?;
        if neighbor.bond == bond_id {
            continue;
        }
        let neighbor_point = coordinate(other_atom.id())?;
        let neighbor_vector = match Vec3::normalized_between(
            center_location,
            neighbor_point,
            center,
            other_atom.id(),
        ) {
            Ok(vector) => vector,
            Err(StereoError::ZeroLengthVector { .. }) => return Ok(result),
            Err(error) => return Err(WedgeError::Geometry(error)),
        };
        let mut angle = reference_vector.signed_projected_angle_to(neighbor_vector);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        let mut insertion = 0;
        while insertion < neighbor_angles.len() && angle > neighbor_angles[insertion].1 {
            insertion += 1;
        }
        neighbor_angles.insert(insertion, (neighbor_bond.id(), angle));
    }
    let source_bond_order: Vec<_> = neighbors
        .iter()
        .map(|neighbor| {
            topology
                .bonds()
                .get(neighbor.bond.index())
                .map(|bond| bond.id())
                .ok_or(WedgeError::BondOutOfRange {
                    bond: neighbor.bond,
                    bond_count: topology.bonds().len(),
                })
        })
        .collect::<Result<_, _>>()?;
    let angle_bond_order: Vec<_> = neighbor_angles.iter().map(|(bond, _)| *bond).collect();
    let mut swaps =
        count_swaps_to_interconvert(&angle_bond_order, &source_bond_order)? as u32 as i32;
    if neighbor_angles.len() == 3 {
        let angle1 = neighbor_angles[1].1;
        let angle2 = neighbor_angles[2].1;
        let angle_tolerance = PI * 1.9 / 180.0;
        if angle2 - angle1 >= PI - angle_tolerance {
            swaps = swaps
                .checked_add(1)
                .ok_or(WedgeError::SourceSignedScoreOverflow {
                    atom: center,
                    operation: "perturbation-order implicit-H increment",
                })?;
        }
    }
    // Preserve the literal signed predicate after native unsigned-to-int cast;
    // negative odd remainders take the source else branch, too.
    result = match (chiral_tag, swaps % 2 == 1) {
        (ChiralTag::TetrahedralCcw, true) | (ChiralTag::TetrahedralCw, false) => {
            BondDirection::BeginDash
        }
        (ChiralTag::TetrahedralCcw, false) | (ChiralTag::TetrahedralCw, true) => {
            BondDirection::BeginWedge
        }
        _ => unreachable!("the source branch above validates tetrahedral chirality"),
    };
    Ok(result)
}

fn get_double_bond_presence<G: StereoGraphAccess>(topology: &G, atom: &G::Atom) -> (u32, u32, u32) {
    // BEGIN RDKIT CPP FUNCTION getDoubleBondPresence
    // RDKit❗✔️: std::tuple<unsigned int, unsigned int, unsigned int> getDoubleBondPresence(
    // RDKit❗✔️:     const ROMol &mol, const Atom &atom) {
    // RDKit❗✔️:   unsigned int hasDouble = 0;
    // RDKit❗✔️:   unsigned int hasKnownDouble = 0;
    // RDKit❗✔️:   unsigned int hasAnyDouble = 0;
    // RDKit❗✔️:   for (const auto bond : mol.atomBonds(&atom)) {
    // RDKit❗✔️:     if (bond->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗✔️:       ++hasDouble;
    // RDKit❗✔️:       if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
    // RDKit❗✔️:         ++hasAnyDouble;
    // RDKit❗✔️:       } else if (bond->getStereo() > Bond::BondStereo::STEREOANY) {
    // RDKit❗✔️:         ++hasKnownDouble;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return std::make_tuple(hasDouble, hasKnownDouble, hasAnyDouble);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION getDoubleBondPresence
    // Behavior remains provisional until the fixed Step 4 branch tests. Complexity
    // matches source: scan this atom's incident bonds with no extra allocation.
    let mut has_double = 0_u32;
    let mut has_known_double = 0_u32;
    let mut has_any_double = 0_u32;

    for neighbor in topology.adjacency().neighbors_of(atom.id().index()).iter() {
        let bond = &topology.bonds()[neighbor.bond.index()];
        if bond.order() != BondOrder::Double {
            continue;
        }

        has_double = has_double.wrapping_add(1);
        if bond.stereo() == BondStereo::Any {
            has_any_double = has_any_double.wrapping_add(1);
        } else if bond.stereo().rdkit_code() > BondStereo::Any.rdkit_code() {
            has_known_double = has_known_double.wrapping_add(1);
        }
    }

    (has_double, has_known_double, has_any_double)
}

fn count_chiral_neighbors<G: StereoGraphAccess>(
    topology: &G,
    no_neighbors: i32,
) -> (bool, Vec<i32>) {
    // BEGIN RDKIT CPP FUNCTION countChiralNbrs
    // RDKit❗✔️: std::pair<bool, INT_VECT> countChiralNbrs(const ROMol &mol, int noNbrs) {
    // RDKit❗✔️:   INT_VECT nChiralNbrs(mol.getNumAtoms(), noNbrs);
    // RDKit❗✔️:
    // RDKit❗✔️:   // start by looking for bonds that are already wedged
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (bond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit❗✔️:         bond->getBondDir() == Bond::BEGINDASH ||
    // RDKit❗✔️:         bond->getBondDir() == Bond::UNKNOWN) {
    // RDKit❗✔️:       if (bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:           bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:         nChiralNbrs[bond->getBeginAtomIdx()] = noNbrs + 1;
    // RDKit❗✔️:       } else if (bond->getEndAtom()->getChiralTag() ==
    // RDKit❗✔️:                      Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:                  bond->getEndAtom()->getChiralTag() ==
    // RDKit❗✔️:                      Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:         nChiralNbrs[bond->getEndAtomIdx()] = noNbrs + 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // now rank atoms by the number of chiral neighbors or Hs they have:
    // RDKit❗✔️:   bool chiNbrs = false;
    // RDKit❗✔️:   for (const auto at : mol.atoms()) {
    // RDKit❗✔️:     if (nChiralNbrs[at->getIdx()] > noNbrs) {
    // RDKit❗✔️:       // std::cerr << " SKIPPING1: " << at->getIdx() << std::endl;
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     auto type = at->getChiralTag();
    // RDKit❗✔️:     if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     nChiralNbrs[at->getIdx()] = 0;
    // RDKit❗✔️:     chiNbrs = true;
    // RDKit❗✔️:     for (const auto nat : mol.atomNeighbors(at)) {
    // RDKit❗✔️:       if (nat->getAtomicNum() == 1) {
    // RDKit❗✔️:         // special case: it's an H... we weight these especially high:
    // RDKit❗✔️:         nChiralNbrs[at->getIdx()] -= 10;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       type = nat->getChiralTag();
    // RDKit❗✔️:       if (type != Atom::CHI_TETRAHEDRAL_CW &&
    // RDKit❗✔️:           type != Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       nChiralNbrs[at->getIdx()] -= 1;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return std::make_pair(chiNbrs, nChiralNbrs);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION countChiralNbrs
    // Behavior remains provisional until the fixed Step 4 branch tests. Complexity
    // matches source: one atom-sized score vector and source-order bond/atom/neighbor
    // scans, O(atoms + bonds), with no nested whole-graph searches.
    let mut chiral_neighbor_counts = vec![no_neighbors; topology.atoms().len()];

    for bond in topology.bonds() {
        if !matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash | BondDirection::Unknown
        ) {
            continue;
        }

        let begin = &topology.atoms()[bond.begin().index()];
        let end = &topology.atoms()[bond.end().index()];
        if matches!(
            begin.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            chiral_neighbor_counts[begin.id().index()] = no_neighbors + 1;
        } else if matches!(
            end.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            chiral_neighbor_counts[end.id().index()] = no_neighbors + 1;
        }
    }

    let mut has_chiral_neighbors = false;
    for atom in topology.atoms() {
        let atom_index = atom.id().index();
        if chiral_neighbor_counts[atom_index] > no_neighbors {
            continue;
        }
        if !matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            continue;
        }

        chiral_neighbor_counts[atom_index] = 0;
        has_chiral_neighbors = true;
        for neighbor in topology.adjacency().neighbors_of(atom_index).iter() {
            let neighbor_atom = &topology.atoms()[neighbor.atom_index];
            if neighbor_atom.atomic_number() == 1 {
                chiral_neighbor_counts[atom_index] -= 10;
                continue;
            }
            if matches!(
                neighbor_atom.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                chiral_neighbor_counts[atom_index] -= 1;
            }
        }
    }

    (has_chiral_neighbors, chiral_neighbor_counts)
}

fn pick_bond_to_wedge(
    topology: &TopologyBlock,
    rings: &mut RingInfo,
    center: AtomId,
    counts: &[i32],
    assignments: &WedgeAssignments,
    no_neighbors: i32,
) -> Result<Option<BondId>, WedgeError> {
    // Existing graph-only interface carries the source's known-empty molecule
    // dictionary. A supplied actual dictionary is passed by the source entry.
    // Native pickBondsToWedge consumes only a nonnegative returned int.
    let id = pick_bond_to_wedge_with_source_properties(
        topology,
        rings,
        center,
        counts,
        assignments,
        no_neighbors,
        None,
    )?;
    Ok((id >= 0).then(|| BondId::new(id as usize)))
}

fn source_wedge_candidate(score: i32, bond_index: usize) -> Result<(i32, i32), WedgeError> {
    // RDKit❗✔️: int bid = bond->getIdx();
    // RDKit❗✔️: nbrScores.emplace_back(nbrScore, bid);
    // Source uint32 ID converts to pinned signed32 int before std::pair's
    // lexicographic minimum. Preserve high-bit IDs; wider detached IDs are
    // a structured transport failure, not a chemistry fallback.
    let bits = u32::try_from(bond_index).map_err(|_| WedgeError::SourceIndexWidth {
        state: "bond",
        index: bond_index,
    })?;
    Ok((score, bits as i32))
}

/// Source scalar over actual topology, ring cache, score/map carriers and an
/// optional actual molecule dictionary. Cache/property errors retain source
/// prefix state; the return is the native signed int, including -1 for empty.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn pick_bond_to_wedge_with_source_properties(
    topology: &TopologyBlock,
    rings: &mut RingInfo,
    center: AtomId,
    counts: &[i32],
    assignments: &WedgeAssignments,
    no_neighbors: i32,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
) -> Result<i32, WedgeError> {
    pick_bond_to_wedge_graph_with_source_properties(
        topology,
        rings,
        center,
        counts,
        assignments,
        no_neighbors,
        properties,
    )
}

fn pick_bond_to_wedge_graph_with_source_properties<G: StereoGraphAccess>(
    topology: &G,
    rings: &mut RingInfo,
    center: AtomId,
    counts: &[i32],
    assignments: &WedgeAssignments,
    no_neighbors: i32,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
) -> Result<i32, WedgeError> {
    // RDKit❗❌: int pickBondToWedge(
    // RDKit❗❌:     const Atom *atom, const ROMol &mol, const INT_VECT &nChiralNbrs,
    // RDKit❗❌:     const std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds,
    // RDKit❗❌:     int noNbrs) {
    // RDKit❗❌:   // here is what we are going to do
    // RDKit❗❌:   // - at each chiral center look for a bond that is begins at the atom and
    // RDKit❗❌:   //   is not yet picked to be wedged for a different chiral center, preferring
    // RDKit❗❌:   //   bonds to Hs
    // RDKit❗❌:   // - if we do not find a bond that begins at the chiral center - we will take
    // RDKit❗❌:   //   the first bond that is not yet picked by any other chiral centers
    // RDKit❗❌:   // we use the orders calculated above to determine which order to do the
    // RDKit❗❌:   // wedging
    // RDKit❗❌:
    // RDKit❗❌:   // we need ring information; make sure findSSSR has been called before
    // RDKit❗❌:   // if not call now
    // RDKit❗❌:   if (!mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit❗❌:     MolOps::findSSSR(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<std::pair<int, int>> nbrScores;
    // RDKit❗❌:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit❗❌:     // can only wedge single bonds:
    // RDKit❗❌:     if (bond->getBondType() != Bond::SINGLE) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     int bid = bond->getIdx();
    // RDKit❗❌:     if (wedgeBonds.find(bid) == wedgeBonds.end()) {
    // RDKit❗❌:       // very strong preference for Hs:
    // RDKit❗❌:       auto *oatom = bond->getOtherAtom(atom);
    // RDKit❗❌:       if (oatom->getAtomicNum() == 1) {
    // RDKit❗❌:         nbrScores.emplace_back(-1000000,
    // RDKit❗❌:                                bid);  // lower than anything else can be
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       // prefer lower atomic numbers with lower degrees and no specified
    // RDKit❗❌:       // chirality:
    // RDKit❗❌:       int nbrScore = oatom->getAtomicNum() + 100 * oatom->getDegree() +
    // RDKit❗❌:                      1000 * ((oatom->getChiralTag() != Atom::CHI_UNSPECIFIED));
    // RDKit❗❌:       // prefer neighbors that are nonchiral or have as few chiral neighbors
    // RDKit❗❌:       // as possible:
    // RDKit❗❌:       int oIdx = oatom->getIdx();
    // RDKit❗❌:       if (nChiralNbrs[oIdx] < noNbrs) {
    // RDKit❗❌:         // the counts are negative, so we have to subtract them off
    // RDKit❗❌:         nbrScore -= 100000 * nChiralNbrs[oIdx];
    // RDKit❗❌:       }
    // RDKit❗❌:       // prefer bonds to non-ring atoms:
    // RDKit❗❌:       nbrScore += 10000 * mol.getRingInfo()->numAtomRings(oIdx);
    // RDKit❗❌:       // prefer non-ring bonds;
    // RDKit❗❌:       nbrScore += 20000 * mol.getRingInfo()->numBondRings(bid);
    // RDKit❗❌:       // prefer bonds to atoms which don't have a double bond from them
    // RDKit❗❌:       auto [hasDoubleBond, hasKnownDoubleBond, hasAnyDoubleBond] =
    // RDKit❗❌:           getDoubleBondPresence(mol, *oatom);
    // RDKit❗❌:       nbrScore += 11000 * hasDoubleBond;
    // RDKit❗❌:       nbrScore += 12000 * hasKnownDoubleBond;
    // RDKit❗❌:       nbrScore += 23000 * hasAnyDoubleBond;
    // RDKit❗❌:
    // RDKit❗❌:       // if at all possible, do not go to marked attachment points
    // RDKit❗❌:       // since they may well be removed when we write a mol block
    // RDKit❗❌:       if (oatom->hasProp(common_properties::_fromAttachPoint)) {
    // RDKit❗❌:         nbrScore += 500000;
    // RDKit❗❌:       }
    // RDKit❗❌:       // std::cerr << "    nrbScore: " << idx << " - " << oIdx << " : "
    // RDKit❗❌:       //           << nbrScore << " nChiralNbrs: " << nChiralNbrs[oIdx]
    // RDKit❗❌:       //           << std::endl;
    // RDKit❗❌:       nbrScores.emplace_back(nbrScore, bid);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // There's still one situation where this whole thing can fail: an unlucky
    // RDKit❗❌:   // situation where all neighbors of all neighbors of an atom are chiral
    // RDKit❗❌:   // and that atom ends up being the last one picked for stereochem
    // RDKit❗❌:   // assignment. This also happens in cases where the chiral atom doesn't
    // RDKit❗❌:   // have all of its neighbors (like when working with partially sanitized
    // RDKit❗❌:   // fragments)
    // RDKit❗❌:   //
    // RDKit❗❌:   // We'll bail here by returning -1
    // RDKit❗❌:   if (nbrScores.empty()) {
    // RDKit❗❌:     return -1;
    // RDKit❗❌:   }
    // RDKit❗❌:   auto minPr = std::min_element(nbrScores.begin(), nbrScores.end());
    // RDKit❗❌:   return minPr->second;
    // RDKit❗❌: }
    // Same candidate vector, one actual incident scan, log map membership and
    // signed pair minimum. Existing SSSR packed-state/output projection costs
    // remain ❌; no graph/property clone or alternate ring/wedge algorithm.
    // Mixed unsigned C++ operations below explicitly use uint32 arithmetic.
    // Pure signed source overflow is undefined: report a structured width
    // error when reached, rather than guessing a wrap/score/selected bond.
    if !rings.is_sssr_or_better() {
        #[cfg(test)]
        drawing_ring_probe::acquire();
        crate::rings::find_sssr_with_source_outputs_from_graph(
            topology, rings, properties, None, false, false,
        )?;
    }
    #[cfg(test)]
    drawing_ring_probe::scoring(rings);
    if center.index() >= topology.atoms().len() {
        return Err(WedgeError::SourceStateIndex {
            state: "center",
            index: center.index(),
            count: topology.atoms().len(),
        });
    }
    let mut scores = Vec::<(i32, i32)>::new();
    for nb in topology.adjacency().neighbors_of(center.index()).iter() {
        let bond = topology
            .bonds()
            .get(nb.bond.index())
            .ok_or(WedgeError::SourceStateIndex {
                state: "bond",
                index: nb.bond.index(),
                count: topology.bonds().len(),
            })?;
        if bond.order() != BondOrder::Single {
            continue;
        }
        // Compute source int bid at the source assignment point, before map
        // membership and hydrogen shortcut, with no full input prevalidation.
        let (_, bid) = source_wedge_candidate(0, bond.id().index())?;
        if assignments.by_bond.contains_key(&bond.id()) {
            continue;
        }
        let atom = topology
            .atoms()
            .get(nb.atom_index)
            .ok_or(WedgeError::SourceStateIndex {
                state: "otherAtom",
                index: nb.atom_index,
                count: topology.atoms().len(),
            })?;
        if atom.atomic_number() == 1 {
            scores.push((-1_000_000, bid));
            continue;
        }
        let degree = topology.adjacency().neighbors_of(atom.id().index()).len() as u32;
        let initial = u32::from(atom.atomic_number())
            .wrapping_add(100u32.wrapping_mul(degree))
            .wrapping_add(1000 * u32::from(atom.chiral_tag() != ChiralTag::Unspecified));
        let mut score = initial as i32;
        let other = atom.id().index();
        let count = *counts.get(other).ok_or(WedgeError::SourceStateIndex {
            state: "nChiralNbrs",
            index: other,
            count: counts.len(),
        })?;
        if count < no_neighbors {
            let penalty =
                100_000i32
                    .checked_mul(count)
                    .ok_or(WedgeError::SourceSignedScoreOverflow {
                        atom: atom.id(),
                        operation: "100000 * nChiralNbrs",
                    })?;
            score = score
                .checked_sub(penalty)
                .ok_or(WedgeError::SourceSignedScoreOverflow {
                    atom: atom.id(),
                    operation: "nbrScore -= chiral neighbor penalty",
                })?;
        }
        score = (score as u32)
            .wrapping_add(10_000u32.wrapping_mul(rings.num_atom_rings(atom.id()) as u32))
            as i32;
        score = (score as u32)
            .wrapping_add(20_000u32.wrapping_mul(rings.num_bond_rings(bond.id()) as u32))
            as i32;
        let (total, known, any) = get_double_bond_presence(topology, atom);
        score = (score as u32).wrapping_add(11_000u32.wrapping_mul(total)) as i32;
        score = (score as u32).wrapping_add(12_000u32.wrapping_mul(known)) as i32;
        score = (score as u32).wrapping_add(23_000u32.wrapping_mul(any)) as i32;
        if atom.prop(ATTACHMENT_POINT_PROPERTY).is_some() {
            score = score
                .checked_add(500_000)
                .ok_or(WedgeError::SourceSignedScoreOverflow {
                    atom: atom.id(),
                    operation: "nbrScore += attachment point penalty",
                })?;
        }
        scores.push((score, bid));
    }
    Ok(scores.into_iter().min().map_or(-1, |(_, bid)| bid))
}

/// Runs native default-conformer picking using actual source append order.
#[doc(hidden)]
pub fn pick_bonds_to_wedge_default_source(
    topology: &mut TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
) -> Result<WedgeAssignments, WedgeError> {
    // RDKit❗✔️: std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> pickBondsToWedge(
    // RDKit❗✔️:     const ROMol &mol, const BondWedgingParameters *params) {
    // RDKit❗✔️:   const Conformer *conf = nullptr;
    // RDKit❗✔️:   if (mol.getNumConformers()) {
    // RDKit❗✔️:     conf = &mol.getConformer();
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return pickBondsToWedge(mol, params, conf);
    // RDKit❗✔️: }
    // RDKit✔️✔️:   inline unsigned int getNumConformers() const {
    // RDKit✔️✔️:     return rdcast<unsigned int>(d_confs.size());
    // RDKit✔️✔️:   }
    // Default native getConformer() always means front(), independent of IDs,
    // dimensional flags, import provenance, or numerical coordinate shape.
    // The model's unique first-source getter owns the actual interleaving;
    // absent mixed-dimensional order is a structural input error, not a guess.
    // Only the selected conformer reaches the explicit source caller.
    // Complexity review: two stored vector lengths, one borrowed front getter,
    // and direct delegation. No conformer/property/topology clone or scan.
    let count = coordinates
        .conformers_2d
        .len()
        .wrapping_add(coordinates.conformers_3d.len()) as u32;
    let conformer = if count != 0 {
        Some(
            match coordinates
                .first_source_conformer()?
                .ok_or(CoordinateValidationError::MissingSourceConformerOrder)?
            {
                cosmolkit_model::CoordinateSourceConformer::TwoD(value) => {
                    AtropisomerConformer::TwoD(value)
                }
                cosmolkit_model::CoordinateSourceConformer::ThreeD(value) => {
                    AtropisomerConformer::ThreeD(value)
                }
            },
        )
    } else {
        None
    };
    pick_bonds_to_wedge_source(topology, conformer, rings, properties)
}

/// Assigns source-default chiral and atropisomer wedges for one detached topology.
pub fn pick_bonds_to_wedge(
    topology: &TopologyBlock,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<WedgeAssignments, WedgeError> {
    pick_bonds_to_wedge_with_ring_info(topology, conformer).map(|(assignments, _)| assignments)
}

/// Assigns source-default wedges and returns the resulting SSSR ring state.
///
/// The CX writer consumes the same ring cache that `pickBondsToWedge` causes
/// the source molecule to retain before its later extension writers run.
pub fn pick_bonds_to_wedge_with_ring_info(
    topology: &TopologyBlock,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<(WedgeAssignments, RingInfo), WedgeError> {
    pick_bonds_to_wedge_with_existing_ring_info(topology, conformer, None)
}

/// Consume a detached ring carrier, preserving trusted SSSR/SymmSSSR rows.
/// The carrier must correspond to the final topology supplied by the caller.
enum SourceWedgeGraph<'a, G: StereoGraphMut> {
    Mutable(&'a mut G),
    Projection(&'a G),
}
impl<G: StereoGraphMut> SourceWedgeGraph<'_, G> {
    fn topology(&self) -> &G {
        match self {
            Self::Mutable(topology) => topology,
            Self::Projection(topology) => topology,
        }
    }
    fn atropisomers(
        &mut self,
        rings: &mut RingInfo,
        properties: Option<&mut cosmolkit_model::MoleculeProperties>,
        conformer: Option<AtropisomerConformer<'_>>,
        wedges: &mut WedgeAssignments,
    ) -> Result<(), WedgeError> {
        match self {
            Self::Mutable(topology) => {
                crate::atropisomer::wedge_bonds_from_atropisomers_graph_source(
                    &mut **topology,
                    rings,
                    properties,
                    conformer,
                    wedges,
                )
                .map_err(WedgeError::from)
            }
            Self::Projection(topology) => {
                let occupied = wedges.by_bond.keys().copied().collect();
                let result =
                    crate::atropisomer::wedge_bonds_from_atropisomers_projected_graph_source(
                        *topology, rings, properties, conformer, &occupied,
                    )?;
                // Only actual native map insertions/replacements join this map.
                // This query returns the native map; graph mutations are observable in
                // the actual mutable source entry.
                let atrop = WedgeAssignments::from_atropisomer_wedge_assignment(result);
                wedges.by_bond.extend(atrop.by_bond);
                wedges.diagnostics.extend(atrop.diagnostics);
                Ok(())
            }
        }
    }
}

fn pick_bonds_to_wedge_kernel<G: StereoGraphMut>(
    graph: &mut SourceWedgeGraph<'_, G>,
    conformer: Option<AtropisomerConformer<'_>>,
    rings: &mut RingInfo,
    mut properties: Option<&mut cosmolkit_model::MoleculeProperties>,
) -> Result<WedgeAssignments, WedgeError> {
    // RDKit❗❌: std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> pickBondsToWedge(
    // RDKit❗❌:     const ROMol &mol, const BondWedgingParameters *params,
    // RDKit❗❌:     const Conformer *conf) {
    // RDKit❗❌:   if (!params) {
    // RDKit❗❌:     params = &defaultWedgingParams;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<unsigned int> indices(mol.getNumAtoms());
    // RDKit❗❌:   std::iota(indices.begin(), indices.end(), 0);
    // RDKit❗❌:   static int noNbrs = 100;
    // RDKit❗❌:   auto [chiNbrs, nChiralNbrs] = detail::countChiralNbrs(mol, noNbrs);
    // RDKit❗❌:   if (chiNbrs) {
    // RDKit❗❌:     std::sort(indices.begin(), indices.end(),
    // RDKit❗❌:               [&nChiralNbrs = nChiralNbrs](auto i1, auto i2) {
    // RDKit❗❌:                 return nChiralNbrs[i1] < nChiralNbrs[i2];
    // RDKit❗❌:               });
    // RDKit❗❌:   }
    // RDKit❗❌:   std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> wedgeInfo;
    // RDKit❗❌:   for (auto idx : indices) {
    // RDKit❗❌:     if (nChiralNbrs[idx] > noNbrs) {
    // RDKit❗❌:       // std::cerr << " SKIPPING2: " << idx << std::endl;
    // RDKit❗❌:       continue;  // already have a wedged bond here
    // RDKit❗❌:     }
    // RDKit❗❌:     auto atom = mol.getAtomWithIdx(idx);
    // RDKit❗❌:     auto type = atom->getChiralTag();
    // RDKit❗❌:     // the indices are ordered such that all chiral atoms come first. If
    // RDKit❗❌:     // this has no chiral flag, we can stop the whole loop:
    // RDKit❗❌:     if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:     auto bnd1 =
    // RDKit❗❌:         detail::pickBondToWedge(atom, mol, nChiralNbrs, wedgeInfo, noNbrs);
    // RDKit❗❌:     if (bnd1 >= 0) {
    // RDKit❗❌:       auto wi = std::unique_ptr<RDKit::Chirality::WedgeInfoChiral>(
    // RDKit❗❌:           new RDKit::Chirality::WedgeInfoChiral(idx));
    // RDKit❗❌:       wedgeInfo[bnd1] = std::move(wi);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   RDKit::Atropisomers::wedgeBondsFromAtropisomers(mol, conf, wedgeInfo);
    // RDKit❗❌:
    // RDKit❗❌:   return wedgeInfo;
    // RDKit❗❌: }
    // Native params is set to the default address when null, then never read
    // again in this complete function. Erasing that unused pointer preserves
    // every parameter-value case; no two-wedge policy is inferred here.
    // Behavior review: exact native score-only sorting, threshold-before-tag
    // loop exit, actual mutable lazy cache/property calls and map writes, then
    // the unique complete source atrop owner. The immutable query uses the
    // same kernel and staged delta without topology/property/cache cloning.
    // Complexity review: native-sized index vector, O(V log V) score sort,
    // degree scans and ordered map. Projection adds occupied set/update/map
    // conversion allocations; actual native mutable carrier avoids those.
    const NO_NEIGHBORS: i32 = 100;
    let mut indices: Vec<_> = (0..graph.topology().atoms().len()).collect();
    let (has_chiral, counts) = count_chiral_neighbors(graph.topology(), NO_NEIGHBORS);
    if has_chiral {
        crate::source_sort::sort_by(&mut indices, |left, right| counts[*left] < counts[*right]);
    }
    let mut wedges = WedgeAssignments::default();
    for index in indices {
        if counts[index] > NO_NEIGHBORS {
            continue;
        }
        let atom = &graph.topology().atoms()[index];
        if !matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            break;
        }
        let center = atom.id();
        let id = pick_bond_to_wedge_graph_with_source_properties(
            graph.topology(),
            rings,
            center,
            &counts,
            &wedges,
            NO_NEIGHBORS,
            properties.as_deref_mut(),
        )?;
        if id >= 0 {
            // RDKit❗✔️:   WedgeInfoBase(int idxInit) : idx(idxInit){};
            // RDKit❗✔️:   WedgeInfoChiral(int atomId) : WedgeInfoBase(atomId){};
            // The source constructor copies the actual center ID only.
            wedges
                .by_bond
                .insert(BondId::new(id as usize), WedgeInfo::Chiral { center });
        }
    }
    // Keep instrumentation at the real native acquisition/use boundaries.
    #[cfg(test)]
    if !rings.is_sssr_or_better() {
        drawing_ring_probe::acquire();
    }
    // The whole source owner promotes here even with no chiral candidate.
    graph.atropisomers(rings, properties, conformer, &mut wedges)?;
    #[cfg(test)]
    drawing_ring_probe::atrop(rings);
    Ok(wedges)
}

/// Runs complete native explicit-conformer picking on actual detached state.
#[doc(hidden)]
pub fn pick_bonds_to_wedge_source(
    topology: &mut TopologyBlock,
    conformer: Option<AtropisomerConformer<'_>>,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
) -> Result<WedgeAssignments, WedgeError> {
    topology.validate()?;
    pick_bonds_to_wedge_kernel(
        &mut SourceWedgeGraph::Mutable(topology),
        conformer,
        rings,
        properties,
    )
}

pub fn pick_bonds_to_wedge_with_existing_ring_info(
    topology: &TopologyBlock,
    conformer: Option<AtropisomerConformer<'_>>,
    rings: Option<RingInfo>,
) -> Result<(WedgeAssignments, RingInfo), WedgeError> {
    topology.validate()?;
    let mut rings = rings.unwrap_or_else(|| {
        RingInfo::new(
            RingFindType::Fast,
            topology.atoms.len(),
            topology.bonds.len(),
        )
    });
    // Preserve this existing checked input's row guards. The separate source
    // entry accepts native sparse rows; internally acquired source rows are
    // passed directly to their source owner without padding or extra guards.
    if rings.is_sssr_or_better() {
        if rings.atom_row_count() != topology.atoms.len() {
            return Err(AtropisomerError::RingAtomRowCount {
                actual: rings.atom_row_count(),
                expected: topology.atoms.len(),
            }
            .into());
        }
        if rings.bond_row_count() != topology.bonds.len() {
            return Err(AtropisomerError::RingBondRowCount {
                actual: rings.bond_row_count(),
                expected: topology.bonds.len(),
            }
            .into());
        }
    }
    let assignments = pick_bonds_to_wedge_kernel(
        &mut SourceWedgeGraph::Projection(topology),
        conformer,
        &mut rings,
        None,
    )?;
    Ok((assignments, rings))
}

#[cfg(test)]
mod drawing_ring_probe {
    use super::*;
    use std::cell::RefCell;

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub(super) struct Observation {
        pub atoms_ptr: usize,
        pub bonds_ptr: usize,
        pub quality: RingFindType,
        pub atoms: Vec<Vec<AtomId>>,
        pub bonds: Vec<Vec<BondId>>,
        pub atom_counts: Vec<usize>,
        pub bond_counts: Vec<usize>,
    }
    pub(super) fn observe(rings: &RingInfo) -> Observation {
        Observation {
            atoms_ptr: rings.atom_rings().as_ptr() as usize,
            bonds_ptr: rings.bond_rings().as_ptr() as usize,
            quality: rings.find_type(),
            atoms: rings.atom_rings().to_vec(),
            bonds: rings.bond_rings().to_vec(),
            atom_counts: (0..rings.atom_row_count())
                .map(|i| rings.num_atom_rings(AtomId::new(i)))
                .collect(),
            bond_counts: (0..rings.bond_row_count())
                .map(|i| rings.num_bond_rings(BondId::new(i)))
                .collect(),
        }
    }
    #[derive(Default, Clone)]
    pub(super) struct State {
        pub acquisitions: usize,
        pub scoring: Vec<Observation>,
        pub atrop: Vec<Observation>,
    }
    thread_local! {
        static STATE: RefCell<State> = RefCell::new(State::default());
    }
    pub(super) fn acquire() {
        STATE.with(|state| state.borrow_mut().acquisitions += 1);
    }
    pub(super) fn scoring(rings: &RingInfo) {
        STATE.with(|state| state.borrow_mut().scoring.push(observe(rings)));
    }
    pub(super) fn atrop(rings: &RingInfo) {
        STATE.with(|state| state.borrow_mut().atrop.push(observe(rings)));
    }
    pub(super) fn state() -> State {
        STATE.with(|state| state.borrow().clone())
    }
}

#[cfg(test)]
mod tests {
    mod drawing_ring_input_tests {
        use super::*;
        use crate::AtropisomerError;
        use crate::wedge::{
            drawing_ring_probe as probe, pick_bonds_to_wedge_with_existing_ring_info,
        };

        const XY: [[f64; 2]; 7] = [
            [0., 0.],
            [1., 0.],
            [-1., 0.],
            [2., 1.],
            [2., -1.],
            [-2., 1.],
            [-2., -1.],
        ];
        const EDGES: [(usize, usize); 7] = [(0, 1), (0, 2), (1, 3), (3, 4), (4, 1), (2, 5), (2, 6)];

        fn graph() -> TopologyBlock {
            topology_from_specs(
                (0..7)
                    .map(|i| {
                        atom(
                            6,
                            if i == 0 {
                                ChiralTag::TetrahedralCw
                            } else {
                                ChiralTag::Unspecified
                            },
                        )
                    })
                    .collect(),
                EDGES
                    .iter()
                    .map(|&(a, b)| edge(a, b, BondOrder::Single))
                    .collect(),
            )
        }

        fn carrier(state: usize) -> Option<RingInfo> {
            let quality = match state {
                4 | 6 => RingFindType::Sssr,
                5 | 7 => RingFindType::SymmSssr,
                3 => RingFindType::Fast,
                _ => RingFindType::OtherOrUnknown,
            };
            if state == 0 {
                return None;
            }
            let mut rings = RingInfo::new(quality, 7, 7);
            if state == 1 {
                rings.reset();
            }
            if state >= 6 {
                rings.add_ring(&[1, 3, 4], &[2, 3, 4]).unwrap();
            }
            Some(rings)
        }

        fn prerequisites(topology: &TopologyBlock, coords: &Conformer2D) {
            topology.validate().unwrap();
            coords.validate_for_atom_count(7).unwrap();
            assert_eq!(coords.id(), 17);
            assert_eq!(coords.coordinates(), &XY);
            for (i, atom) in topology.atoms.iter().enumerate() {
                assert_eq!(atom.id(), AtomId::new(i));
                assert_eq!(atom.element(), Element::C);
                assert_eq!(
                    atom.chiral_tag(),
                    if i == 0 {
                        ChiralTag::TetrahedralCw
                    } else {
                        ChiralTag::Unspecified
                    }
                );
                assert!(!atom.is_aromatic());
                assert!(atom.props().is_empty());
                assert_eq!(
                    topology.adjacency.neighbors_of(i).len(),
                    [2, 3, 3, 2, 2, 1, 1][i]
                );
            }
            for (i, bond) in topology.bonds.iter().enumerate() {
                assert_eq!(bond.id(), BondId::new(i));
                assert_eq!((bond.begin().index(), bond.end().index()), EDGES[i]);
                assert_eq!(bond.order(), BondOrder::Single);
                assert_eq!(bond.direction(), BondDirection::None);
                assert_eq!(bond.stereo(), BondStereo::None);
                assert!(!bond.is_aromatic());
                assert!(bond.props().is_empty());
            }
        }

        #[test]
        fn drawing_ring_wedge_eighteen_actual_calls() {
            let mut calls = 0;
            for state in 0..8 {
                for with_coords in [false, true] {
                    let topology = graph();
                    let coords = Conformer2D::new(17, XY.to_vec());
                    prerequisites(&topology, &coords);
                    let topology_snapshot = topology.clone();
                    let coords_snapshot = coords.clone();
                    let coordinate_bits = coords
                        .coordinates()
                        .iter()
                        .map(|p| p.map(f64::to_bits))
                        .collect::<Vec<_>>();
                    let rings = carrier(state);
                    let ring_snapshot = rings.clone();
                    let move_baseline = if state >= 4 {
                        Some(probe::observe(rings.as_ref().unwrap()))
                    } else {
                        None
                    };
                    if let Some(rings) = &rings {
                        assert_eq!(rings.is_initialized(), state != 1);
                        if state != 1 {
                            assert_eq!(rings.atom_row_count(), 7);
                            assert_eq!(rings.bond_row_count(), 7);
                            assert_eq!(
                                rings.atom_rings(),
                                if state >= 6 {
                                    vec![vec![AtomId::new(1), AtomId::new(3), AtomId::new(4)]]
                                } else {
                                    vec![]
                                }
                            );
                            assert_eq!(
                                rings.bond_rings(),
                                if state >= 6 {
                                    vec![vec![BondId::new(2), BondId::new(3), BondId::new(4)]]
                                } else {
                                    vec![]
                                }
                            );
                            for i in 0..7 {
                                assert_eq!(
                                    rings.num_atom_rings(AtomId::new(i)),
                                    usize::from(state >= 6 && [1, 3, 4].contains(&i))
                                );
                                assert_eq!(
                                    rings.num_bond_rings(BondId::new(i)),
                                    usize::from(state >= 6 && [2, 3, 4].contains(&i))
                                );
                            }
                        }
                    }
                    let baseline = probe::state();
                    let result = pick_bonds_to_wedge_with_existing_ring_info(
                        &topology,
                        with_coords.then_some(AtropisomerConformer::TwoD(&coords)),
                        rings,
                    );
                    calls += 1;
                    let observed = probe::state();
                    assert_eq!(topology, topology_snapshot);
                    assert_eq!(coords, coords_snapshot);
                    assert_eq!(
                        coords
                            .coordinates()
                            .iter()
                            .map(|p| p.map(f64::to_bits))
                            .collect::<Vec<_>>(),
                        coordinate_bits
                    );
                    let (wedges, returned) = result.unwrap();
                    let expected_bond = if state == 4 || state == 5 { 0 } else { 1 };
                    assert_eq!(
                        wedges.iter().collect::<Vec<_>>(),
                        vec![(
                            BondId::new(expected_bond),
                            &WedgeInfo::Chiral {
                                center: AtomId::new(0)
                            }
                        )]
                    );
                    assert!(wedges.diagnostics().is_empty());
                    assert_eq!(
                        observed.acquisitions - baseline.acquisitions,
                        usize::from(state < 4)
                    );
                    assert_eq!(observed.scoring.len() - baseline.scoring.len(), 1);
                    assert_eq!(observed.atrop.len() - baseline.atrop.len(), 1);
                    let final_observation = probe::observe(&returned);
                    assert_eq!(observed.scoring[baseline.scoring.len()], final_observation);
                    assert_eq!(observed.atrop[baseline.atrop.len()], final_observation);
                    if let Some(move_baseline) = move_baseline {
                        assert_eq!(final_observation, move_baseline);
                        assert_eq!(Some(returned), ring_snapshot);
                    } else {
                        assert_eq!(returned.find_type(), RingFindType::Sssr);
                        assert_eq!(
                            returned.atom_rings(),
                            &[vec![AtomId::new(1), AtomId::new(3), AtomId::new(4)]]
                        );
                        assert_eq!(
                            returned.bond_rings(),
                            &[vec![BondId::new(2), BondId::new(3), BondId::new(4)]]
                        );
                        for i in 0..7 {
                            assert_eq!(
                                returned.num_atom_rings(AtomId::new(i)),
                                usize::from([1, 3, 4].contains(&i))
                            );
                            assert_eq!(
                                returned.num_bond_rings(BondId::new(i)),
                                usize::from([2, 3, 4].contains(&i))
                            );
                        }
                    }
                }
            }
            assert_eq!(calls, 16);
            for (atoms, bonds, expected) in [
                (
                    8,
                    7,
                    AtropisomerError::RingAtomRowCount {
                        actual: 8,
                        expected: 7,
                    },
                ),
                (
                    7,
                    8,
                    AtropisomerError::RingBondRowCount {
                        actual: 8,
                        expected: 7,
                    },
                ),
            ] {
                let topology = graph();
                let coords = Conformer2D::new(17, XY.to_vec());
                prerequisites(&topology, &coords);
                let snapshot = (topology.clone(), coords.clone());
                let rings = RingInfo::new(RingFindType::Sssr, atoms, bonds);
                assert_eq!(rings.atom_row_count(), atoms);
                assert_eq!(rings.bond_row_count(), bonds);
                let baseline = probe::state();
                let result = pick_bonds_to_wedge_with_existing_ring_info(
                    &topology,
                    Some(AtropisomerConformer::TwoD(&coords)),
                    Some(rings),
                );
                calls += 1;
                assert_eq!(result, Err(WedgeError::Atropisomer(expected)));
                let observed = probe::state();
                assert_eq!(observed.acquisitions, baseline.acquisitions);
                assert_eq!(observed.scoring, baseline.scoring);
                assert_eq!(observed.atrop, baseline.atrop);
                assert_eq!((topology, coords), snapshot);
            }
            assert_eq!(calls, 18);
        }
    }
    use super::{
        AtropisomerConformer, AtropisomerWedgeUpdate, CrossedBondContext, MolFileBondStereoInfo,
        WedgeAssignments, WedgeError, WedgeInfo, bond_get_dir_code, can_be_stereo_bond,
        count_chiral_neighbors, determine_bond_wedge_state, get_double_bond_presence,
        get_molfile_bond_stereo_info, pick_bond_to_wedge, pick_bonds_to_wedge,
    };
    use crate::ValenceAssignment;
    use crate::get_all_atom_ids_for_stereo_group;
    use crate::potential_stereo;
    use crate::rings::{RingFindType, RingInfo, RingSearchParams, find_sssr};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D, StereoGroup,
        StereoGroupKind, TopologyBlock,
    };
    use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element};

    fn atom(atomic_number: u8, chiral_tag: ChiralTag) -> AtomSpec {
        AtomSpec::new(Element::from_atomic_number(atomic_number).expect("fixture element"))
            .with_chiral_tag(chiral_tag)
    }

    fn topology_from_specs(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = bond_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("valid fixed wedge fixture")
    }

    fn topology(atom_specs: &[(u8, ChiralTag)], bond_specs: Vec<BondSpec>) -> TopologyBlock {
        let atoms = atom_specs
            .iter()
            .map(|(atomic_number, chiral_tag)| atom(*atomic_number, *chiral_tag))
            .collect();
        topology_from_specs(atoms, bond_specs)
    }

    fn edge(begin: usize, end: usize, order: BondOrder) -> BondSpec {
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
    }

    fn pick(
        topology: &TopologyBlock,
        chiral_neighbor_counts: &[i32],
        occupied_bonds: &[usize],
    ) -> (Option<BondId>, RingInfo) {
        let no_neighbors = 99;
        let mut rings = RingInfo::new(
            RingFindType::Fast,
            topology.atoms.len(),
            topology.bonds.len(),
        );
        let wedge_assignments = WedgeAssignments {
            by_bond: occupied_bonds
                .iter()
                .copied()
                .map(|bond| {
                    (
                        BondId::new(bond),
                        WedgeInfo::Chiral {
                            center: AtomId::new(0),
                        },
                    )
                })
                .collect(),
            diagnostics: Vec::new(),
        };
        let selected = pick_bond_to_wedge(
            topology,
            &mut rings,
            AtomId::new(0),
            chiral_neighbor_counts,
            &wedge_assignments,
            no_neighbors,
        )
        .expect("source SSSR selection succeeds on valid topology");
        (selected, rings)
    }

    fn unweighted_neighbor_counts(topology: &TopologyBlock) -> Vec<i32> {
        vec![99; topology.atoms.len()]
    }

    fn crossed_valence(explicit_valence: &[i32], implicit_hydrogens: &[i32]) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit_valence.to_vec(),
            implicit_hydrogens: implicit_hydrogens.to_vec(),
        }
    }

    fn ordinary_crossed_topology() -> TopologyBlock {
        topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Double),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        )
    }

    fn should_cross_first_bond(
        topology: &TopologyBlock,
        valence: &ValenceAssignment,
        use_legacy_stereo_perception: bool,
    ) -> Result<bool, WedgeError> {
        let rings = find_sssr(topology, &RingSearchParams::default())
            .expect("valid fixed crossed-bond topology has SSSR ring information");
        CrossedBondContext::new(topology, valence, &rings, use_legacy_stereo_perception)?
            .should_be_crossed_bond(BondId::new(0))
    }

    fn potential_first_bond(topology: &TopologyBlock, valence: &ValenceAssignment) -> bool {
        let rings = find_sssr(topology, &RingSearchParams::default())
            .expect("valid fixed potential-bond topology has SSSR ring information");
        potential_stereo::is_potential_bond(topology, valence, &rings, &topology.bonds[0])
            .expect("fixed potential-bond inputs have source-shaped valence and rings")
    }

    fn ranked_atom(cip_rank: Option<&str>, chiral_atom_rank: Option<&str>) -> AtomSpec {
        let mut atom = atom(6, ChiralTag::Unspecified);
        if let Some(rank) = cip_rank {
            atom = atom
                .with_prop("_CIPRank", rank)
                .expect("fixed CIP rank property is valid");
        }
        if let Some(rank) = chiral_atom_rank {
            atom = atom
                .with_prop("_ChiralAtomRank", rank)
                .expect("fixed modern rank property is valid");
        }
        atom
    }

    fn ranked_crossed_topology(
        cip_ranks: [Option<&str>; 2],
        chiral_atom_ranks: [Option<&str>; 2],
        target_direction: BondDirection,
        neighbor_bonds: Vec<BondSpec>,
    ) -> TopologyBlock {
        let mut bonds = vec![edge(0, 1, BondOrder::Double).with_direction(target_direction)];
        bonds.extend(neighbor_bonds);
        topology_from_specs(
            vec![
                atom(6, ChiralTag::Unspecified),
                atom(6, ChiralTag::Unspecified),
                ranked_atom(cip_ranks[0], chiral_atom_ranks[0]),
                ranked_atom(cip_ranks[1], chiral_atom_ranks[1]),
                atom(6, ChiralTag::Unspecified),
            ],
            bonds,
        )
    }

    fn ordinary_rank_neighbors() -> Vec<BondSpec> {
        vec![
            edge(0, 2, BondOrder::Single),
            edge(0, 3, BondOrder::Single),
            edge(1, 4, BondOrder::Single),
        ]
    }

    fn ring_crossed_topology(size: usize) -> TopologyBlock {
        let atoms = (0..size).map(|_| atom(6, ChiralTag::Unspecified)).collect();
        let mut bonds = vec![edge(0, 1, BondOrder::Double)];
        for begin in 1..size {
            let end = if begin + 1 == size { 0 } else { begin + 1 };
            bonds.push(edge(begin, end, BondOrder::Single));
        }
        topology_from_specs(atoms, bonds)
    }

    fn ring_crossed_valence(size: usize) -> ValenceAssignment {
        let mut explicit_valence = vec![2; size];
        explicit_valence[0] = 3;
        explicit_valence[1] = 3;
        let mut implicit_hydrogens = vec![2; size];
        implicit_hydrogens[0] = 1;
        implicit_hydrogens[1] = 1;
        crossed_valence(&explicit_valence, &implicit_hydrogens)
    }

    fn wedge_geometry_topology(chiral_tag: ChiralTag, direction: BondDirection) -> TopologyBlock {
        topology(
            &[
                (6, chiral_tag),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single).with_direction(direction),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
            ],
        )
    }

    fn wedge_geometry_2d(topology: &TopologyBlock, coordinates: Vec<[f64; 2]>) -> BondDirection {
        determine_bond_wedge_state(
            topology,
            BondId::new(0),
            AtomId::new(0),
            Some(AtropisomerConformer::TwoD(&Conformer2D::new(
                0,
                coordinates,
            ))),
        )
        .expect("fixed wedge geometry is valid")
    }

    fn molfile_info(
        topology: &TopologyBlock,
        valence: &ValenceAssignment,
        assignments: &WedgeAssignments,
        bond_id: BondId,
        conformer: Option<AtropisomerConformer<'_>>,
    ) -> MolFileBondStereoInfo {
        let rings = find_sssr(topology, &RingSearchParams::default())
            .expect("fixed molfile fixture has valid ring information");
        let crossed_bonds = CrossedBondContext::new(topology, valence, &rings, true)
            .expect("fixed molfile topology is valid");
        get_molfile_bond_stereo_info(&crossed_bonds, assignments, bond_id, conformer)
            .expect("fixed molfile bond stereo input is valid")
    }

    fn assignment_map(entries: Vec<(usize, WedgeInfo)>) -> WedgeAssignments {
        WedgeAssignments {
            by_bond: entries
                .into_iter()
                .map(|(bond, info)| (BondId::new(bond), info))
                .collect(),
            diagnostics: Vec::new(),
        }
    }

    #[test]
    fn wedge_preexisting_wedge_and_unknown_directions_exclude_chiral_centers() {
        // WedgeBonds.cpp::countChiralNbrs first marks BEGINWEDGE, BEGINDASH,
        // and UNKNOWN endpoints, preferring the bond's begin atom.
        for direction in [
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::Unknown,
        ] {
            let molecule = topology(
                &[(6, ChiralTag::TetrahedralCw), (6, ChiralTag::Unspecified)],
                vec![edge(0, 1, BondOrder::Single).with_direction(direction)],
            );
            assert_eq!(
                count_chiral_neighbors(&molecule, 4),
                (false, vec![5, 4]),
                "direction {direction:?}"
            );
        }

        let end_center = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::TetrahedralCcw)],
            vec![edge(0, 1, BondOrder::Single).with_direction(BondDirection::Unknown)],
        );
        assert_eq!(count_chiral_neighbors(&end_center, 4), (false, vec![4, 5]));

        let both_centers = topology(
            &[
                (6, ChiralTag::TetrahedralCw),
                (6, ChiralTag::TetrahedralCcw),
            ],
            vec![edge(0, 1, BondOrder::Single).with_direction(BondDirection::Unknown)],
        );
        assert_eq!(
            count_chiral_neighbors(&both_centers, 4),
            (true, vec![5, -1])
        );
    }

    #[test]
    fn wedge_non_wedge_directions_leave_chiral_centers_eligible() {
        // Cis/trans directions and EITHERDOUBLE are not in the source's
        // preexisting-wedge exclusion set.
        for direction in [
            BondDirection::None,
            BondDirection::EndDownRight,
            BondDirection::EndUpRight,
            BondDirection::EitherDouble,
        ] {
            let molecule = topology(
                &[(6, ChiralTag::TetrahedralCw), (6, ChiralTag::Unspecified)],
                vec![edge(0, 1, BondOrder::Single).with_direction(direction)],
            );
            assert_eq!(
                count_chiral_neighbors(&molecule, 9),
                (true, vec![0, 9]),
                "direction {direction:?}"
            );
        }
    }

    #[test]
    fn wedge_neighbor_scores_ignore_ordinary_neighbors_and_weight_chiral_atoms_and_hydrogens() {
        let ordinary_neighbor = topology(
            &[(6, ChiralTag::TetrahedralCw), (6, ChiralTag::Tetrahedral)],
            vec![edge(0, 1, BondOrder::Single)],
        );
        assert_eq!(
            count_chiral_neighbors(&ordinary_neighbor, 7),
            (true, vec![0, 7])
        );

        let chiral_hydrogen_and_ordinary_neighbors = topology(
            &[
                (6, ChiralTag::TetrahedralCw),
                (6, ChiralTag::TetrahedralCcw),
                (1, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
            ],
        );
        assert_eq!(
            count_chiral_neighbors(&chiral_hydrogen_and_ordinary_neighbors, 7),
            (true, vec![-11, -1, 7, 7])
        );
    }

    #[test]
    fn source598_count_chiral_neighbors_weights_chiral_hydrogen_before_tag_without_order_filter() {
        // countChiralNbrs weighs atomic number 1 before checking the neighbor's
        // tag, and its atomNeighbors pass contains no zero/dative exclusion.
        let molecule = topology(
            &[
                (6, ChiralTag::TetrahedralCw),
                (1, ChiralTag::TetrahedralCcw),
                (6, ChiralTag::SquarePlanar),
            ],
            vec![edge(0, 1, BondOrder::Zero), edge(0, 2, BondOrder::Dative)],
        );
        assert_eq!(
            count_chiral_neighbors(&molecule, 100),
            (true, vec![-10, -1, 100])
        );
    }

    #[test]
    fn wedge_double_bond_presence_counts_every_pinned_stereo_category() {
        // The pinned helper counts STEREOANY separately and every enum value
        // above it as known; STEREONONE contributes only to hasDouble.
        let categories = [
            (BondStereo::None, (1, 0, 0)),
            (BondStereo::Any, (1, 0, 1)),
            (BondStereo::Z, (1, 1, 0)),
            (BondStereo::E, (1, 1, 0)),
            (BondStereo::Cis, (1, 1, 0)),
            (BondStereo::Trans, (1, 1, 0)),
            (BondStereo::AtropCw, (1, 1, 0)),
            (BondStereo::AtropCcw, (1, 1, 0)),
        ];

        for (stereo, expected) in categories {
            let mut double_bond = edge(0, 1, BondOrder::Double).with_stereo(stereo);
            if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
                double_bond = double_bond.with_stereo_atoms(AtomId::new(2), AtomId::new(3));
            }
            let molecule = topology(
                &[
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                ],
                vec![
                    double_bond,
                    edge(0, 2, BondOrder::Single),
                    edge(1, 3, BondOrder::Single),
                ],
            );
            assert_eq!(
                get_double_bond_presence(&molecule, &molecule.atoms[0]),
                expected,
                "stereo {stereo:?}"
            );
        }
    }

    #[test]
    fn source606_candidate_ties_compare_signed_native_int_bond_ids() {
        let a = super::source_wedge_candidate(7, 0x8000_0000).unwrap();
        let b = super::source_wedge_candidate(7, 1).unwrap();
        assert_eq!(a, (7, i32::MIN));
        assert_eq!([b, a].into_iter().min(), Some(a));
        assert_eq!(
            super::source_wedge_candidate(-1, u32::MAX as usize).unwrap(),
            (-1, -1)
        );
        #[cfg(target_pointer_width = "64")]
        assert!(matches!(
            super::source_wedge_candidate(0, (u32::MAX as usize) + 1),
            Err(WedgeError::SourceIndexWidth { state: "bond", .. })
        ));
    }
    #[test]
    fn source606_real_ring_cache_reset_and_property_error_precede_scoring() {
        let g = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Single)],
        );
        let mut rings = RingInfo::new(RingFindType::Fast, 2, 1);
        let mut props = cosmolkit_model::MoleculeProperties::default();
        props.set_prop("__computedProps", false).unwrap();
        let error = super::pick_bond_to_wedge_with_source_properties(
            &g,
            &mut rings,
            AtomId::new(0),
            &[100; 2],
            &WedgeAssignments::default(),
            100,
            Some(&mut props),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            WedgeError::RingFinding(crate::RingFindingError::MoleculeProperty(_))
        ));
        assert_eq!(rings.find_type(), RingFindType::Sssr);
        assert!(rings.is_initialized());
        assert_eq!(rings.atom_row_count(), 0);
        assert_eq!(
            props.prop("__computedProps"),
            Some(&cosmolkit_model::PropertyValue::Bool(false))
        );
    }
    #[test]
    fn source606_trusted_sssr_skips_dictionary_and_hydrogen_skips_missing_score_rows() {
        let g = topology(
            &[(6, ChiralTag::Unspecified), (1, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Single)],
        );
        let mut rings = RingInfo::new(RingFindType::SymmSssr, 2, 1);
        let original = rings.clone();
        let mut props = cosmolkit_model::MoleculeProperties::default();
        props.set_prop("__computedProps", false).unwrap();
        assert_eq!(
            super::pick_bond_to_wedge_with_source_properties(
                &g,
                &mut rings,
                AtomId::new(0),
                &[],
                &WedgeAssignments::default(),
                100,
                Some(&mut props)
            )
            .unwrap(),
            0
        );
        assert_eq!(rings, original);
        assert_eq!(
            props.prop("__computedProps"),
            Some(&cosmolkit_model::PropertyValue::Bool(false))
        );
    }
    #[test]
    fn source606_actual_properties_clear_only_source_computed_extra_rings_when_preparing() {
        let g = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Single)],
        );
        let mut rings = RingInfo::new(RingFindType::Fast, 2, 1);
        let mut props = cosmolkit_model::MoleculeProperties::default();
        props.set_computed_prop("extraRings", 1i32).unwrap();
        props.set_prop("ordinary", "kept").unwrap();
        assert_eq!(
            super::pick_bond_to_wedge_with_source_properties(
                &g,
                &mut rings,
                AtomId::new(0),
                &[100; 2],
                &WedgeAssignments::default(),
                100,
                Some(&mut props)
            )
            .unwrap(),
            0
        );
        assert!(props.prop("extraRings").is_none());
        assert_eq!(
            props.prop("ordinary"),
            Some(&cosmolkit_model::PropertyValue::String("kept".into()))
        );
        assert_eq!(rings.find_type(), RingFindType::Sssr);
    }
    #[test]
    fn source606_empty_candidates_return_native_minus_one_after_cache_preparation() {
        let g = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Double)],
        );
        let mut rings = RingInfo::new(RingFindType::Fast, 2, 1);
        assert_eq!(
            super::pick_bond_to_wedge_with_source_properties(
                &g,
                &mut rings,
                AtomId::new(0),
                &[],
                &WedgeAssignments::default(),
                100,
                None
            )
            .unwrap(),
            -1
        );
        assert_eq!(rings.find_type(), RingFindType::Sssr);
    }
    #[test]
    fn wedge_pick_prefers_hydrogen_then_lower_atomic_number_degree_and_unspecified_tag() {
        let hydrogen = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (1, ChiralTag::Unspecified),
            ],
            vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
        );
        assert_eq!(
            pick(&hydrogen, &unweighted_neighbor_counts(&hydrogen), &[]).0,
            Some(BondId::new(1)),
            "source gives hydrogen the fixed -1000000 score before other terms"
        );

        let atomic_number = topology(
            &[
                (6, ChiralTag::Unspecified),
                (9, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
        );
        assert_eq!(
            pick(
                &atomic_number,
                &unweighted_neighbor_counts(&atomic_number),
                &[]
            )
            .0,
            Some(BondId::new(1)),
            "lower atomic number wins when other score terms tie"
        );

        let degree = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(2, 3, BondOrder::Single),
            ],
        );
        assert_eq!(
            pick(&degree, &unweighted_neighbor_counts(&degree), &[]).0,
            Some(BondId::new(0)),
            "lower graph degree wins when element and tag tie"
        );

        let specified_tag = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Tetrahedral),
            ],
            vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
        );
        assert_eq!(
            pick(
                &specified_tag,
                &unweighted_neighbor_counts(&specified_tag),
                &[]
            )
            .0,
            Some(BondId::new(0)),
            "any non-unspecified chiral tag incurs the source 1000-point penalty"
        );
    }

    #[test]
    fn wedge_pick_penalizes_neighbors_with_chiral_neighbors() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
        );
        let mut counts = unweighted_neighbor_counts(&molecule);
        counts[1] = -1;
        assert_eq!(
            pick(&molecule, &counts, &[]).0,
            Some(BondId::new(1)),
            "negative chiral-neighbor counts are subtracted as a positive penalty"
        );
    }

    #[test]
    fn wedge_pick_prefers_non_ring_atoms_and_promotes_ring_information() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
                edge(3, 4, BondOrder::Single),
                edge(4, 1, BondOrder::Single),
                edge(2, 5, BondOrder::Single),
                edge(2, 6, BondOrder::Single),
            ],
        );
        let (selected, rings) = pick(&molecule, &unweighted_neighbor_counts(&molecule), &[]);
        assert_eq!(selected, Some(BondId::new(1)));
        assert!(
            rings.is_sssr_or_better(),
            "Fast ring info is promoted by findSSSR"
        );
        assert_eq!(rings.num_atom_rings(AtomId::new(1)), 1);
        assert_eq!(rings.num_atom_rings(AtomId::new(2)), 0);
        assert_eq!(rings.num_bond_rings(BondId::new(0)), 0);
        assert_eq!(rings.num_bond_rings(BondId::new(1)), 0);
    }

    #[test]
    fn wedge_pick_penalizes_ring_bonds_even_when_both_neighbors_are_ring_atoms() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Tetrahedral),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
                edge(1, 4, BondOrder::Single),
                edge(2, 5, BondOrder::Single),
                edge(5, 6, BondOrder::Single),
                edge(6, 2, BondOrder::Single),
            ],
        );
        let (selected, rings) = pick(&molecule, &unweighted_neighbor_counts(&molecule), &[]);
        assert_eq!(selected, Some(BondId::new(1)));
        assert_eq!(rings.num_atom_rings(AtomId::new(1)), 1);
        assert_eq!(rings.num_atom_rings(AtomId::new(2)), 1);
        assert_eq!(rings.num_bond_rings(BondId::new(0)), 1);
        assert_eq!(rings.num_bond_rings(BondId::new(1)), 0);
    }

    #[test]
    fn wedge_pick_applies_double_presence_known_and_any_penalties() {
        let double_presence = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Double),
                edge(1, 4, BondOrder::Single),
                edge(2, 5, BondOrder::Single),
                edge(5, 6, BondOrder::Single),
                edge(6, 2, BondOrder::Single),
            ],
        );
        assert_eq!(
            pick(
                &double_presence,
                &unweighted_neighbor_counts(&double_presence),
                &[]
            )
            .0,
            Some(BondId::new(1)),
            "one double bond adds 11000, outweighing one ring-atom penalty of 10000"
        );

        for stereo in [BondStereo::Z, BondStereo::Any] {
            let known_or_any = topology(
                &[
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Tetrahedral),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                ],
                vec![
                    edge(0, 1, BondOrder::Single),
                    edge(0, 2, BondOrder::Single),
                    edge(0, 3, BondOrder::Single),
                    edge(2, 3, BondOrder::Single),
                    edge(1, 4, BondOrder::Single),
                    edge(4, 5, BondOrder::Single),
                    edge(5, 1, BondOrder::Single),
                    edge(1, 6, BondOrder::Double).with_stereo(stereo),
                    edge(2, 7, BondOrder::Single),
                ],
            );
            let (selected, rings) = pick(
                &known_or_any,
                &unweighted_neighbor_counts(&known_or_any),
                &[],
            );
            assert_eq!(
                selected,
                Some(BondId::new(1)),
                "known/any double score {stereo:?} exceeds the 20000 ring-bond penalty"
            );
            assert_eq!(rings.num_atom_rings(AtomId::new(1)), 1);
            assert_eq!(rings.num_atom_rings(AtomId::new(2)), 1);
            assert_eq!(rings.num_bond_rings(BondId::new(0)), 0);
            assert_eq!(rings.num_bond_rings(BondId::new(1)), 1);
        }
    }

    #[test]
    fn wedge_pick_uses_attachment_property_presence_as_a_penalty() {
        // WedgeBonds.cpp uses the symbol whose pinned types.h value is
        // _fromAttchpt. Keep the original ordinary-key input as a no-penalty
        // regression and add the actual reserved property including empty value.
        for (property_key, expected_bond) in [("_fromAttachPoint", 0), ("_fromAttchpt", 1)] {
            let marked = atom(6, ChiralTag::Unspecified)
                .with_prop(property_key, "")
                .expect("empty raw attachment property is representable");
            let molecule = topology_from_specs(
                vec![
                    atom(6, ChiralTag::Unspecified),
                    marked,
                    atom(6, ChiralTag::Unspecified),
                ],
                vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
            );
            assert_eq!(
                pick(&molecule, &unweighted_neighbor_counts(&molecule), &[]).0,
                Some(BondId::new(expected_bond)),
                "property presence penalizes the marked atom even when its value is empty"
            );
        }
    }

    #[test]
    fn wedge_pick_skips_occupied_and_non_single_bonds_and_returns_none_when_empty() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (1, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
        );
        assert_eq!(
            pick(&molecule, &unweighted_neighbor_counts(&molecule), &[0]).0,
            Some(BondId::new(1)),
            "an occupied hydrogen bond is excluded before the hydrogen preference"
        );
        assert_eq!(
            pick(&molecule, &unweighted_neighbor_counts(&molecule), &[0, 1]).0,
            None,
            "all single-bond candidates occupied returns the source no-candidate value"
        );

        let only_double = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Double)],
        );
        assert_eq!(
            pick(&only_double, &unweighted_neighbor_counts(&only_double), &[]).0,
            None,
            "non-single bonds are not candidates"
        );
    }

    #[test]
    fn wedge_pick_breaks_equal_scores_by_source_bond_id() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![edge(0, 2, BondOrder::Single), edge(0, 1, BondOrder::Single)],
        );
        assert_eq!(
            pick(&molecule, &unweighted_neighbor_counts(&molecule), &[]).0,
            Some(BondId::new(0)),
            "equal scores use std::pair ordering, score then source bond index"
        );
    }

    #[test]
    fn wedge_geometry_maps_cw_and_ccw_for_even_source_order_parity() {
        let coordinates = vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [-0.5, 3.0_f64.sqrt() / 2.0],
            [-0.5, -3.0_f64.sqrt() / 2.0],
        ];

        for (tag, expected) in [
            (ChiralTag::TetrahedralCcw, BondDirection::BeginWedge),
            (ChiralTag::TetrahedralCw, BondDirection::BeginDash),
        ] {
            let topology = wedge_geometry_topology(tag, BondDirection::None);
            assert_eq!(
                wedge_geometry_2d(&topology, coordinates.clone()),
                expected,
                "{tag:?}"
            );
        }
    }

    #[test]
    fn wedge_geometry_projects_both_3d_is3d_flags_to_xy() {
        let topology = wedge_geometry_topology(ChiralTag::TetrahedralCcw, BondDirection::None);
        let expected = BondDirection::BeginWedge;
        let xy = [
            [0.0, 0.0],
            [1.0, 0.0],
            [-0.5, 3.0_f64.sqrt() / 2.0],
            [-0.5, -3.0_f64.sqrt() / 2.0],
        ];
        for is_3d in [false, true] {
            let conformer = Conformer3D::new(
                2,
                xy.into_iter()
                    .enumerate()
                    .map(|(atom, [x, y])| [x, y, (atom as f64 + 1.0) * 1000.0])
                    .collect(),
                is_3d,
            );
            assert_eq!(
                determine_bond_wedge_state(
                    &topology,
                    BondId::new(0),
                    AtomId::new(0),
                    Some(AtropisomerConformer::ThreeD(&conformer)),
                )
                .expect("fixed 3D coordinate rows are valid"),
                expected,
                "source zeros z even when is3D={is_3d}"
            );
        }
    }

    #[test]
    fn wedge_geometry_returns_existing_direction_without_coordinates_or_on_overlap() {
        let topology =
            wedge_geometry_topology(ChiralTag::TetrahedralCcw, BondDirection::EndDownRight);
        assert_eq!(
            determine_bond_wedge_state(&topology, BondId::new(0), AtomId::new(0), None)
                .expect("no conformer returns the original bond direction"),
            BondDirection::EndDownRight
        );

        for coordinates in [
            vec![[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            vec![[0.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
        ] {
            assert_eq!(
                wedge_geometry_2d(&topology, coordinates),
                BondDirection::EndDownRight,
                "overlapping reference or neighboring atoms retain the source direction"
            );
        }
    }

    #[test]
    fn wedge_geometry_equal_angles_insert_before_prior_neighbors_and_change_parity() {
        let coordinates = vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]];
        for (tag, expected) in [
            (ChiralTag::TetrahedralCcw, BondDirection::BeginDash),
            (ChiralTag::TetrahedralCw, BondDirection::BeginWedge),
        ] {
            let topology = wedge_geometry_topology(tag, BondDirection::None);
            assert_eq!(
                wedge_geometry_2d(&topology, coordinates.clone()),
                expected,
                "equal 90-degree angles reverse the source bond-row tie order for {tag:?}"
            );
        }
    }

    #[test]
    fn wedge_geometry_three_coordinate_implicit_h_correction_straddles_178_1_degrees() {
        let topology = wedge_geometry_topology(ChiralTag::TetrahedralCcw, BondDirection::None);
        for (separation_degrees, expected) in [
            (178.09_f64, BondDirection::BeginWedge),
            (178.11_f64, BondDirection::BeginDash),
        ] {
            let first_angle = 1.0_f64.to_radians();
            let second_angle = first_angle + separation_degrees.to_radians();
            let coordinates = vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [first_angle.cos(), first_angle.sin()],
                [second_angle.cos(), second_angle.sin()],
            ];
            assert_eq!(
                wedge_geometry_2d(&topology, coordinates),
                expected,
                "the source adds the implicit-H swap at >= 178.1 degrees; separation={separation_degrees}"
            );
        }
    }

    #[test]
    fn wedge_assignment_resolves_competing_chiral_centers_by_neighbor_score_and_occupancy() {
        let molecule = topology(
            &[
                (6, ChiralTag::TetrahedralCw),
                (6, ChiralTag::TetrahedralCcw),
                (6, ChiralTag::TetrahedralCw),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Double),
                edge(1, 3, BondOrder::Double),
            ],
        );

        let assignments = pick_bonds_to_wedge(&molecule, None)
            .expect("pinned selection handles this valid detached topology");

        assert_eq!(assignments.by_bond.len(), 1);
        assert_eq!(
            assignments.by_bond.get(&BondId::new(0)),
            Some(&WedgeInfo::Chiral {
                center: AtomId::new(0),
            }),
            "center 0 has two chiral neighbors and is processed before center 1, taking their only shared single bond"
        );
        assert!(assignments.diagnostics.is_empty());
    }

    #[test]
    fn wedge_assignment_handles_equal_chiral_scores_without_a_center_tie_break() {
        let molecule = topology(
            &[
                (6, ChiralTag::TetrahedralCw),
                (9, ChiralTag::Unspecified),
                (8, ChiralTag::Unspecified),
                (6, ChiralTag::TetrahedralCcw),
                (9, ChiralTag::Unspecified),
                (8, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(3, 4, BondOrder::Single),
                edge(3, 5, BondOrder::Single),
            ],
        );

        let assignments = pick_bonds_to_wedge(&molecule, None)
            .expect("independent equal-score centers can be assigned in either source sort order");

        assert_eq!(assignments.by_bond.len(), 2);
        assert_eq!(
            assignments.by_bond.get(&BondId::new(1)),
            Some(&WedgeInfo::Chiral {
                center: AtomId::new(0),
            })
        );
        assert_eq!(
            assignments.by_bond.get(&BondId::new(3)),
            Some(&WedgeInfo::Chiral {
                center: AtomId::new(3),
            })
        );
        assert!(assignments.diagnostics.is_empty());
    }

    #[test]
    fn wedge_assignment_skips_existing_wedge_and_preserves_other_source_directions() {
        let molecule = topology(
            &[
                (6, ChiralTag::TetrahedralCw),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::TetrahedralCcw),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single).with_direction(BondDirection::BeginWedge),
                edge(2, 3, BondOrder::Single).with_direction(BondDirection::EndUpRight),
            ],
        );

        let assignments =
            pick_bonds_to_wedge(&molecule, None).expect("existing source directions are valid");

        assert_eq!(assignments.by_bond.len(), 1);
        assert_eq!(
            assignments.by_bond.get(&BondId::new(1)),
            Some(&WedgeInfo::Chiral {
                center: AtomId::new(2),
            }),
            "ENDUPRIGHT does not count as an already wedged center"
        );
        assert_eq!(
            determine_bond_wedge_state(&molecule, BondId::new(1), AtomId::new(2), None)
                .expect("the no-conformer source branch retains the bond direction"),
            BondDirection::EndUpRight
        );
        assert!(!assignments.by_bond.contains_key(&BondId::new(0)));
    }

    #[test]
    fn wedge_assignment_wedges_atrop_only_topology_without_chiral_centers() {
        let mut molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single).with_stereo(BondStereo::AtropCw),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        );

        // This fixture models zero implicit H, not an uninitialized native cache.
        // Atom::getNumImplicitHs returns zero immediately for noImplicit atoms.
        molecule.atoms[0].set_no_implicit(true);
        molecule.atoms[1].set_no_implicit(true);
        let assignments = pick_bonds_to_wedge(&molecule, None)
            .expect("the source atrop stage runs without tetrahedral centers");

        assert_eq!(assignments.by_bond.len(), 1);
        assert_eq!(
            assignments.by_bond.get(&BondId::new(2)),
            Some(&WedgeInfo::Atropisomer {
                update: AtropisomerWedgeUpdate {
                    bond: BondId::new(2),
                    begin: AtomId::new(1),
                    end: AtomId::new(3),
                    direction: BondDirection::BeginWedge,
                    atropisomer_bond: BondId::new(0),
                },
            })
        );
        assert!(assignments.diagnostics.is_empty());
    }

    #[test]
    fn wedge_assignment_passes_chiral_occupancy_to_atrop_owner() {
        // Original ordinary-key input remains: oxygen scores108 versus
        // axial carbon206, so its bond3 receives the chiral wedge. With actual
        // _fromAttchpt, O/F gain500000 and bond0 wins. Atrop's unoccupied bond2
        // is preferred to the no-conf BeginDash alternative in both source cases.
        for (property_key, expected_chiral_bond) in [("_fromAttachPoint", 3), ("_fromAttchpt", 0)] {
            let oxygen = atom(8, ChiralTag::Unspecified)
                .with_prop(property_key, "")
                .expect("source attachment property can be present with an empty value");
            let fluorine = atom(9, ChiralTag::Unspecified)
                .with_prop(property_key, "")
                .expect("source attachment property can be present with an empty value");
            let mut molecule = topology_from_specs(
                vec![
                    atom(6, ChiralTag::TetrahedralCw),
                    atom(6, ChiralTag::Unspecified),
                    atom(6, ChiralTag::Unspecified),
                    atom(6, ChiralTag::Unspecified),
                    oxygen,
                    fluorine,
                ],
                vec![
                    edge(0, 1, BondOrder::Single),
                    edge(1, 2, BondOrder::Single).with_stereo(BondStereo::AtropCw),
                    edge(2, 3, BondOrder::Single),
                    edge(0, 4, BondOrder::Single),
                    edge(0, 5, BondOrder::Single),
                ],
            );

            // This fixture models zero implicit H, not an uninitialized native cache.
            // Atom::getNumImplicitHs returns zero immediately for noImplicit atoms.
            molecule.atoms[1].set_no_implicit(true);
            molecule.atoms[2].set_no_implicit(true);
            let assignments = pick_bonds_to_wedge(&molecule, None)
                .expect("the atrop owner composes after tetrahedral selection");

            assert_eq!(assignments.by_bond.len(), 2);
            assert_eq!(
                assignments.by_bond.get(&BondId::new(expected_chiral_bond)),
                Some(&WedgeInfo::Chiral {
                    center: AtomId::new(0),
                }),
                "source score selects the bond to the axial endpoint before attachment-marked alternatives"
            );
            assert_eq!(
                assignments.by_bond.get(&BondId::new(2)),
                Some(&WedgeInfo::Atropisomer {
                    update: AtropisomerWedgeUpdate {
                        bond: BondId::new(2),
                        begin: AtomId::new(2),
                        end: AtomId::new(3),
                        direction: BondDirection::BeginWedge,
                        atropisomer_bond: BondId::new(1),
                    },
                }),
                "the already selected chiral bond is occupied; the existing atrop owner selects the remaining carrier"
            );
            assert!(assignments.diagnostics.is_empty());
        }
    }

    #[test]
    fn wedge_molfile_direction_codes_match_the_pinned_switch_and_default() {
        let assignments = WedgeAssignments::default();
        let valence = crossed_valence(&[1, 1], &[3, 3]);
        for (direction, expected_code) in [
            (BondDirection::None, 0),
            (BondDirection::BeginWedge, 1),
            (BondDirection::BeginDash, 6),
            (BondDirection::Unknown, 4),
            (BondDirection::EitherDouble, 3),
            (BondDirection::EndDownRight, 0),
            (BondDirection::EndUpRight, 0),
        ] {
            let molecule = topology(
                &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
                vec![edge(0, 1, BondOrder::Single).with_direction(direction)],
            );
            let info = molfile_info(&molecule, &valence, &assignments, BondId::new(0), None);
            assert_eq!(info.direction, direction, "direction {direction:?}");
            assert_eq!(
                info.direction_code, expected_code,
                "source molfile code for {direction:?}"
            );
            assert!(!info.reverse, "an unassigned bond never reverses");
        }
    }

    #[test]
    fn wedge_molfile_crossed_double_uses_eitherdouble_and_code_three() {
        let molecule = ordinary_crossed_topology();
        let valence = crossed_valence(&[3, 3, 1, 1], &[1, 1, 3, 3]);
        let info = molfile_info(
            &molecule,
            &valence,
            &WedgeAssignments::default(),
            BondId::new(0),
            None,
        );
        assert_eq!(info.direction, BondDirection::EitherDouble);
        assert_eq!(info.direction_code, 3);
        assert!(!info.reverse);
    }

    #[test]
    fn wedge_molfile_reverses_end_chiral_assignments_but_never_atrop_assignments() {
        for center in [0, 1] {
            let mut atom_specs = vec![(6, ChiralTag::Unspecified); 4];
            atom_specs[center].1 = ChiralTag::TetrahedralCw;
            let molecule = topology(
                &atom_specs,
                vec![
                    edge(0, 1, BondOrder::Single),
                    edge(center, 2, BondOrder::Single),
                    edge(center, 3, BondOrder::Single),
                ],
            );
            let coordinates = if center == 0 {
                vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
            } else {
                vec![[1.0, 0.0], [0.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
            };
            let conformer = Conformer2D::new(0, coordinates);
            let source_conformer = AtropisomerConformer::TwoD(&conformer);
            let assignments = pick_bonds_to_wedge(&molecule, Some(source_conformer))
                .expect("the source picker assigns this valid chiral center");
            assert_eq!(
                assignments.get(BondId::new(0)),
                Some(&WedgeInfo::Chiral {
                    center: AtomId::new(center),
                })
            );
            let mut explicit_valence = vec![1; 4];
            explicit_valence[center] = 3;
            let mut implicit_hydrogens = vec![3; 4];
            implicit_hydrogens[center] = 1;
            let valence = crossed_valence(&explicit_valence, &implicit_hydrogens);
            let info = molfile_info(
                &molecule,
                &valence,
                &assignments,
                BondId::new(0),
                Some(source_conformer),
            );
            assert!(matches!(
                info.direction,
                BondDirection::BeginWedge | BondDirection::BeginDash
            ));
            assert_eq!(info.direction_code, bond_get_dir_code(info.direction));
            assert_eq!(info.reverse, center == 1);
        }

        let mut atrop = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Single).with_stereo(BondStereo::AtropCw),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        );
        // This fixture models zero implicit H, not an uninitialized native cache.
        // Atom::getNumImplicitHs returns zero immediately for noImplicit atoms.
        atrop.atoms[0].set_no_implicit(true);
        atrop.atoms[1].set_no_implicit(true);
        let atrop_assignments = pick_bonds_to_wedge(&atrop, None)
            .expect("source atrop assignment exists without a conformer");
        let Some(WedgeInfo::Atropisomer { update }) = atrop_assignments.get(BondId::new(2)) else {
            panic!("pinned assignment must be atropisomeric");
        };
        assert_eq!(update.direction, BondDirection::BeginWedge);
        let atrop_valence = crossed_valence(&[2, 2, 1, 1], &[2, 2, 3, 3]);
        let atrop_info = molfile_info(
            &atrop,
            &atrop_valence,
            &atrop_assignments,
            BondId::new(2),
            None,
        );
        assert_eq!(atrop_info.direction, update.direction);
        assert_eq!(atrop_info.direction_code, 1);
        assert!(!atrop_info.reverse);
    }

    #[test]
    fn wedge_stereo_any_respects_source_unknown_orientation_at_both_endpoints() {
        let no_unknown = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Single).with_stereo(BondStereo::Any)],
        );
        let no_valence = crossed_valence(&[], &[]);
        assert!(should_cross_first_bond(&no_unknown, &no_valence, true).unwrap());

        for (unknown, expected, case) in [
            (edge(0, 2, BondOrder::Single), false, "outgoing at begin"),
            (edge(2, 0, BondOrder::Single), true, "incoming at begin"),
            (edge(1, 2, BondOrder::Single), false, "outgoing at end"),
            (edge(2, 1, BondOrder::Single), true, "incoming at end"),
        ] {
            let molecule = topology(
                &[
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                    (6, ChiralTag::Unspecified),
                ],
                vec![
                    edge(0, 1, BondOrder::Single).with_stereo(BondStereo::Any),
                    unknown.with_direction(BondDirection::Unknown),
                ],
            );
            assert_eq!(
                should_cross_first_bond(&molecule, &no_valence, false).unwrap(),
                expected,
                "pinned STEREOANY UNKNOWN orientation: {case}"
            );
        }

        let target_unknown = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![
                edge(0, 1, BondOrder::Single)
                    .with_stereo(BondStereo::Any)
                    .with_direction(BondDirection::Unknown),
            ],
        );
        assert!(
            !should_cross_first_bond(&target_unknown, &no_valence, true).unwrap(),
            "the source endpoint adjacency scan includes the target bond itself"
        );

        let specified = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Double).with_stereo(BondStereo::Z),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        );
        assert!(
            !should_cross_first_bond(&specified, &no_valence, true).unwrap(),
            "all non-None/non-Any stereo suppresses automatic crossing"
        );

        let non_double = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Single)],
        );
        assert!(
            !should_cross_first_bond(&non_double, &no_valence, true).unwrap(),
            "an unspecified non-double bond is not a crossed bond"
        );
    }

    #[test]
    fn wedge_potential_double_bond_keeps_degree_hydrogen_and_ring_thresholds() {
        let low_degree = topology(
            &[(6, ChiralTag::Unspecified), (6, ChiralTag::Unspecified)],
            vec![edge(0, 1, BondOrder::Double)],
        );
        assert!(!potential_first_bond(
            &low_degree,
            &crossed_valence(&[2, 2], &[0, 0])
        ));

        let high_degree = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Double),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
                edge(1, 4, BondOrder::Single),
            ],
        );
        assert!(!potential_first_bond(
            &high_degree,
            &crossed_valence(&[3, 3, 1, 1, 1], &[1, 1, 3, 3, 3])
        ));

        for hydrogen_endpoint in [0, 1] {
            let mut atom_specs = vec![atom(6, ChiralTag::Unspecified); 5];
            atom_specs[2] = atom(1, ChiralTag::Unspecified);
            atom_specs[3] = atom(1, ChiralTag::Unspecified);
            let bonds = if hydrogen_endpoint == 0 {
                vec![
                    edge(0, 1, BondOrder::Double),
                    edge(0, 2, BondOrder::Single),
                    edge(0, 3, BondOrder::Single),
                    edge(1, 4, BondOrder::Single),
                ]
            } else {
                vec![
                    edge(0, 1, BondOrder::Double),
                    edge(0, 4, BondOrder::Single),
                    edge(1, 2, BondOrder::Single),
                    edge(1, 3, BondOrder::Single),
                ]
            };
            let molecule = topology_from_specs(atom_specs, bonds);
            let implicit = if hydrogen_endpoint == 0 {
                [0, 1, 0, 0, 3]
            } else {
                [1, 0, 0, 0, 3]
            };
            assert!(
                !potential_first_bond(&molecule, &crossed_valence(&[3, 3, 1, 1, 1], &implicit)),
                "getTotalNumHs(true) counts explicit neighboring hydrogen atoms at endpoint {hydrogen_endpoint}"
            );
        }

        for (size, expected) in [(7, false), (8, true)] {
            let molecule = ring_crossed_topology(size);
            let valence = ring_crossed_valence(size);
            assert_eq!(
                potential_first_bond(&molecule, &valence),
                expected,
                "ring {size}"
            );
            assert_eq!(
                should_cross_first_bond(&molecule, &valence, true).unwrap(),
                expected,
                "pinned minimum ring size is eight"
            );
        }

        let acyclic = ordinary_crossed_topology();
        let fast_ring_info =
            RingInfo::new(RingFindType::Fast, acyclic.atoms.len(), acyclic.bonds.len());
        let valence = crossed_valence(&[3, 3, 1, 1], &[1, 1, 3, 3]);
        assert!(fast_ring_info.is_initialized());
        assert!(
            CrossedBondContext::new(&acyclic, &valence, &fast_ring_info, true)
                .unwrap()
                .should_be_crossed_bond(BondId::new(0))
                .unwrap(),
            "the pinned helper needs initialized ring lists, not the separate public symmetric-SSSR contract"
        );
    }

    #[test]
    fn wedge_crossed_bond_either_double_precedes_raw_final_valence_and_both_endpoints_gate() {
        let either_double = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Double).with_direction(BondDirection::EitherDouble),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        );
        assert!(
            should_cross_first_bond(
                &either_double,
                &crossed_valence(&[0, 0, 0, 0], &[0, 0, 0, 0]),
                false
            )
            .unwrap()
        );

        let ordinary = ordinary_crossed_topology();
        for (endpoint_valences, expected) in [
            ([3, 3], true),
            ([2, 3], false),
            ([3, 2], false),
            ([2, 2], false),
        ] {
            let valence = crossed_valence(
                &[endpoint_valences[0], endpoint_valences[1], 1, 1],
                &[1, 1, 3, 3],
            );
            assert_eq!(
                should_cross_first_bond(&ordinary, &valence, true).unwrap(),
                expected,
                "source total valence minus total degree must equal one at both endpoints: {endpoint_valences:?}"
            );
        }

        let terminal = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![edge(0, 1, BondOrder::Double), edge(1, 2, BondOrder::Single)],
        );
        let terminal_valence = crossed_valence(&[2, 3, 1], &[1, 1, 3]);
        assert!(potential_first_bond(&terminal, &terminal_valence));
        assert!(
            !should_cross_first_bond(&terminal, &terminal_valence, true).unwrap(),
            "the later source getDegree() > 1 gate excludes a terminal endpoint"
        );
    }

    #[test]
    fn wedge_can_be_stereo_bond_selects_only_the_active_profile_and_nonnegative_ranks() {
        let legacy_duplicate = ranked_crossed_topology(
            [Some("7"), Some("7")],
            [Some("8"), Some("9")],
            BondDirection::None,
            ordinary_rank_neighbors(),
        );
        assert!(!can_be_stereo_bond(&legacy_duplicate, &legacy_duplicate.bonds[0], true).unwrap());
        assert!(can_be_stereo_bond(&legacy_duplicate, &legacy_duplicate.bonds[0], false).unwrap());
        let duplicate_valence = crossed_valence(&[4, 3, 1, 1, 1], &[0, 1, 3, 3, 3]);
        assert!(!should_cross_first_bond(&legacy_duplicate, &duplicate_valence, true).unwrap());
        assert!(should_cross_first_bond(&legacy_duplicate, &duplicate_valence, false).unwrap());

        let modern_duplicate = ranked_crossed_topology(
            [Some("8"), Some("9")],
            [Some("7"), Some("7")],
            BondDirection::None,
            ordinary_rank_neighbors(),
        );
        assert!(can_be_stereo_bond(&modern_duplicate, &modern_duplicate.bonds[0], true).unwrap());
        assert!(!can_be_stereo_bond(&modern_duplicate, &modern_duplicate.bonds[0], false).unwrap());

        let absent = ranked_crossed_topology(
            [None, None],
            [None, None],
            BondDirection::None,
            ordinary_rank_neighbors(),
        );
        assert!(can_be_stereo_bond(&absent, &absent.bonds[0], true).unwrap());
        let negative_duplicate = ranked_crossed_topology(
            [Some("-1"), Some("-1")],
            [None, None],
            BondDirection::None,
            ordinary_rank_neighbors(),
        );
        assert!(
            can_be_stereo_bond(&negative_duplicate, &negative_duplicate.bonds[0], true).unwrap()
        );

        let malformed_in_unselected_profile = ranked_crossed_topology(
            [Some("not-an-integer"), None],
            [Some("3"), Some("4")],
            BondDirection::None,
            ordinary_rank_neighbors(),
        );
        assert!(
            can_be_stereo_bond(
                &malformed_in_unselected_profile,
                &malformed_in_unselected_profile.bonds[0],
                false
            )
            .unwrap()
        );
        assert!(matches!(
            can_be_stereo_bond(
                &malformed_in_unselected_profile,
                &malformed_in_unselected_profile.bonds[0],
                true
            ),
            Err(WedgeError::InvalidStereoRank {
                atom,
                property,
                ..
            }) if atom == AtomId::new(2) && property == "_CIPRank"
        ));
    }

    #[test]
    fn wedge_can_be_stereo_bond_applies_only_source_single_neighbor_directions() {
        for direction in [BondDirection::EndUpRight, BondDirection::EndDownRight] {
            let molecule = ranked_crossed_topology(
                [None, None],
                [None, None],
                BondDirection::None,
                vec![
                    edge(0, 2, BondOrder::Single).with_direction(direction),
                    edge(0, 3, BondOrder::Single),
                    edge(1, 4, BondOrder::Single),
                ],
            );
            assert!(
                !can_be_stereo_bond(&molecule, &molecule.bonds[0], true).unwrap(),
                "single neighbor {direction:?} suppresses crossing"
            );
        }

        for (neighbor, expected, case) in [
            (
                edge(0, 2, BondOrder::Single),
                false,
                "begin-endpoint outgoing",
            ),
            (
                edge(2, 0, BondOrder::Single),
                true,
                "begin-endpoint incoming",
            ),
            (
                edge(1, 4, BondOrder::Single),
                false,
                "end-endpoint outgoing",
            ),
            (edge(4, 1, BondOrder::Single), true, "end-endpoint incoming"),
        ] {
            let mut neighbors = ordinary_rank_neighbors();
            let replace_index =
                if neighbor.begin() == AtomId::new(0) || neighbor.end() == AtomId::new(0) {
                    0
                } else {
                    2
                };
            neighbors[replace_index] = neighbor.with_direction(BondDirection::Unknown);
            let molecule =
                ranked_crossed_topology([None, None], [None, None], BondDirection::None, neighbors);
            assert_eq!(
                can_be_stereo_bond(&molecule, &molecule.bonds[0], true).unwrap(),
                expected,
                "pinned UNKNOWN direction rule: {case}"
            );
        }

        let ignored_nonsingle_direction = ranked_crossed_topology(
            [None, None],
            [None, None],
            BondDirection::None,
            vec![
                edge(0, 2, BondOrder::Triple).with_direction(BondDirection::EndUpRight),
                edge(0, 3, BondOrder::Single),
                edge(1, 4, BondOrder::Single),
            ],
        );
        assert!(
            can_be_stereo_bond(
                &ignored_nonsingle_direction,
                &ignored_nonsingle_direction.bonds[0],
                true
            )
            .unwrap()
        );

        let ignored_target_direction = ranked_crossed_topology(
            [None, None],
            [None, None],
            BondDirection::EndDownRight,
            ordinary_rank_neighbors(),
        );
        assert!(
            can_be_stereo_bond(
                &ignored_target_direction,
                &ignored_target_direction.bonds[0],
                true
            )
            .unwrap()
        );

        let aromatic_target = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Aromatic),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        );
        assert!(can_be_stereo_bond(&aromatic_target, &aromatic_target.bonds[0], true).unwrap());

        let triple_target = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 1, BondOrder::Triple),
                edge(0, 2, BondOrder::Single),
                edge(1, 3, BondOrder::Single),
            ],
        );
        assert!(!can_be_stereo_bond(&triple_target, &triple_target.bonds[0], true).unwrap());
    }

    #[test]
    fn wedge_stereo_group_atom_only_preserves_source_member_order() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            Vec::new(),
        );
        let group = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(5), AtomId::new(2)],
            Vec::new(),
        )
        .expect("valid distinct stereo members");

        assert_eq!(
            get_all_atom_ids_for_stereo_group(&molecule, &group, &WedgeAssignments::default())
                .unwrap(),
            vec![AtomId::new(5), AtomId::new(2)]
        );
    }

    #[test]
    fn wedge_stereo_group_bond_only_uses_assigned_atrop_map_entries() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 2, BondOrder::Single).with_stereo(BondStereo::AtropCw),
                edge(2, 1, BondOrder::Single),
                edge(2, 3, BondOrder::Single),
                edge(3, 4, BondOrder::Single).with_stereo(BondStereo::AtropCcw),
                edge(3, 5, BondOrder::Single).with_direction(BondDirection::Unknown),
            ],
        );
        let group = StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), vec![BondId::new(2)])
            .expect("valid distinct stereo members");
        let assignments = assignment_map(vec![(
            1,
            WedgeInfo::Atropisomer {
                update: AtropisomerWedgeUpdate {
                    bond: BondId::new(1),
                    begin: AtomId::new(2),
                    end: AtomId::new(1),
                    direction: BondDirection::BeginWedge,
                    atropisomer_bond: BondId::new(2),
                },
            },
        )]);

        assert_eq!(
            get_all_atom_ids_for_stereo_group(&molecule, &group, &assignments).unwrap(),
            vec![AtomId::new(2)]
        );
    }

    #[test]
    fn wedge_stereo_group_mixed_members_keep_source_endpoint_and_adjacency_order() {
        let molecule = topology(
            &[
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
                (6, ChiralTag::Unspecified),
            ],
            vec![
                edge(0, 2, BondOrder::Single),
                edge(1, 2, BondOrder::Single).with_direction(BondDirection::BeginWedge),
                edge(2, 3, BondOrder::Single),
                edge(3, 4, BondOrder::Single),
                edge(3, 5, BondOrder::Single).with_stereo(BondStereo::AtropCw),
                edge(4, 6, BondOrder::Single).with_stereo(BondStereo::AtropCcw),
                edge(4, 7, BondOrder::Single).with_direction(BondDirection::Unknown),
                edge(4, 8, BondOrder::Single).with_direction(BondDirection::BeginDash),
            ],
        );
        let group = StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(9)],
            vec![BondId::new(2), BondId::new(3)],
        )
        .expect("valid distinct stereo members");
        let assignments = assignment_map(vec![
            (
                4,
                WedgeInfo::Atropisomer {
                    update: AtropisomerWedgeUpdate {
                        bond: BondId::new(4),
                        begin: AtomId::new(3),
                        end: AtomId::new(5),
                        direction: BondDirection::BeginDash,
                        atropisomer_bond: BondId::new(2),
                    },
                },
            ),
            (
                6,
                WedgeInfo::Chiral {
                    center: AtomId::new(4),
                },
            ),
        ]);

        assert_eq!(
            get_all_atom_ids_for_stereo_group(&molecule, &group, &assignments).unwrap(),
            vec![
                AtomId::new(9),
                AtomId::new(2),
                AtomId::new(3),
                AtomId::new(4)
            ]
        );
    }
}

#[cfg(test)]
mod q01_b1_wedge_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec, PropertyValueKind};
    use cosmolkit_types::Element;
    #[test]
    fn q01_b1_wedge_vector_rank_is_named_wrong_kind() {
        let mut topology = TopologyBlock::try_from_parts(
            (0..4)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            [
                (0, 1, BondOrder::Double),
                (0, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
            ]
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                )
            })
            .collect(),
            vec![],
            vec![],
        )
        .unwrap();
        for (key, legacy) in [("_CIPRank", true), ("_ChiralAtomRank", false)] {
            assert!(can_be_stereo_bond(&topology, &topology.bonds[0], legacy).unwrap());
            topology.atoms[2].set_prop(key, vec![1_i32]).unwrap();
            assert_eq!(
                can_be_stereo_bond(&topology, &topology.bonds[0], legacy),
                Err(WedgeError::InvalidStereoRankType {
                    atom: AtomId::new(2),
                    property: key,
                    kind: PropertyValueKind::IntVector
                })
            );
            assert!(can_be_stereo_bond(&topology, &topology.bonds[0], !legacy).unwrap());
            topology.atoms[2].clear_prop(key);
        }
    }
}

#[cfg(test)]
mod uint_wedge_proposed_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    #[test]
    fn proposed_uint_wedge_reads_selected_profile_and_single_neighbors_lazily() {
        for value in [0_u32, 1, 2147483646, 2147483647, 2147483648, 4294967295] {
            let atoms = (0..4)
                .map(|index| {
                    Atom::from_spec(
                        AtomId::new(index),
                        if index == 2 {
                            AtomSpec::new(Element::C)
                                .with_prop("_CIPRank", PropertyValue::UInt(value))
                                .unwrap()
                        } else {
                            AtomSpec::new(Element::C)
                        },
                    )
                })
                .collect();
            let mut graph = TopologyBlock::try_from_parts(
                atoms,
                vec![
                    Bond::from_spec(
                        BondId::new(0),
                        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
                    ),
                    Bond::from_spec(
                        BondId::new(1),
                        BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
                    ),
                    Bond::from_spec(
                        BondId::new(2),
                        BondSpec::new(AtomId::new(1), AtomId::new(3), BondOrder::Single),
                    ),
                ],
                vec![],
                vec![],
            )
            .unwrap();
            let expected = if value <= 2147483647 {
                Ok(true)
            } else {
                Err(WedgeError::UnsignedRankOverflow {
                    atom: AtomId::new(2),
                    property: "_CIPRank",
                    value,
                })
            };
            let original = graph.clone();
            assert_eq!(can_be_stereo_bond(&graph, &graph.bonds[0], true), expected);
            assert_eq!(graph, original);
            assert_eq!(can_be_stereo_bond(&graph, &graph.bonds[0], false), Ok(true));
            graph.bonds[1].set_direction(BondDirection::EndUpRight);
            assert_eq!(can_be_stereo_bond(&graph, &graph.bonds[0], true), Ok(false));
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;

    fn graph(key: &str, v: u32) -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    if i == 2 {
                        AtomSpec::new(Element::C)
                            .with_prop(key, PropertyValue::UInt(v))
                            .unwrap()
                    } else {
                        AtomSpec::new(Element::C)
                    },
                )
            })
            .collect();
        let bonds = [
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, (a, b, o))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), o),
            )
        })
        .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }

    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__CIPRank_0
    #[test]
    fn uint_cell_signed_consumer_core_wedge__ciprank_0_wedge() {
        let g = graph("_CIPRank", 0_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], true), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_CIPRank", PropertyValue::Int(0))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(false));
        tied.atoms[4]
            .set_prop("_CIPRank", PropertyValue::Int(1))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__ChiralAtomRank_0
    #[test]
    fn uint_cell_signed_consumer_core_wedge__chiralatomrank_0_wedge() {
        let g = graph("_ChiralAtomRank", 0_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], false), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_ChiralAtomRank", PropertyValue::Int(0))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(false));
        tied.atoms[4]
            .set_prop("_ChiralAtomRank", PropertyValue::Int(1))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__CIPRank_1
    #[test]
    fn uint_cell_signed_consumer_core_wedge__ciprank_1_wedge() {
        let g = graph("_CIPRank", 1_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], true), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_CIPRank", PropertyValue::Int(1))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(false));
        tied.atoms[4]
            .set_prop("_CIPRank", PropertyValue::Int(0))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__ChiralAtomRank_1
    #[test]
    fn uint_cell_signed_consumer_core_wedge__chiralatomrank_1_wedge() {
        let g = graph("_ChiralAtomRank", 1_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], false), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_ChiralAtomRank", PropertyValue::Int(1))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(false));
        tied.atoms[4]
            .set_prop("_ChiralAtomRank", PropertyValue::Int(0))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__CIPRank_2147483646
    #[test]
    fn uint_cell_signed_consumer_core_wedge__ciprank_2147483646_wedge() {
        let g = graph("_CIPRank", 2147483646_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], true), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_CIPRank", PropertyValue::Int(2147483646))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(false));
        tied.atoms[4]
            .set_prop("_CIPRank", PropertyValue::Int(2147483645))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__ChiralAtomRank_2147483646
    #[test]
    fn uint_cell_signed_consumer_core_wedge__chiralatomrank_2147483646_wedge() {
        let g = graph("_ChiralAtomRank", 2147483646_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], false), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_ChiralAtomRank", PropertyValue::Int(2147483646))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(false));
        tied.atoms[4]
            .set_prop("_ChiralAtomRank", PropertyValue::Int(2147483645))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__CIPRank_2147483647
    #[test]
    fn uint_cell_signed_consumer_core_wedge__ciprank_2147483647_wedge() {
        let g = graph("_CIPRank", 2147483647_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], true), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_CIPRank", PropertyValue::Int(2147483647))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(false));
        tied.atoms[4]
            .set_prop("_CIPRank", PropertyValue::Int(2147483646))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], true), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__ChiralAtomRank_2147483647
    #[test]
    fn uint_cell_signed_consumer_core_wedge__chiralatomrank_2147483647_wedge() {
        let g = graph("_ChiralAtomRank", 2147483647_u32);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], false), Ok(true));
        assert_eq!(g, before);
        let mut tied = g.clone();
        tied.atoms.push(Atom::from_spec(
            AtomId::new(4),
            AtomSpec::new(Element::C)
                .with_prop("_ChiralAtomRank", PropertyValue::Int(2147483647))
                .unwrap(),
        ));
        tied.bonds.push(Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ));
        let mut tied =
            TopologyBlock::try_from_parts(tied.atoms, tied.bonds, vec![], vec![]).unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(false));
        tied.atoms[4]
            .set_prop("_ChiralAtomRank", PropertyValue::Int(2147483646))
            .unwrap();
        assert_eq!(can_be_stereo_bond(&tied, &tied.bonds[0], false), Ok(true));
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__CIPRank_2147483648
    #[test]
    fn uint_cell_signed_consumer_core_wedge__ciprank_2147483648_wedge() {
        let g = graph("_CIPRank", 2147483648_u32);
        let before = g.clone();
        assert_eq!(
            can_be_stereo_bond(&g, &g.bonds[0], true),
            Err(WedgeError::UnsignedRankOverflow {
                atom: AtomId::new(2),
                property: "_CIPRank",
                value: 2147483648_u32
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__ChiralAtomRank_2147483648
    #[test]
    fn uint_cell_signed_consumer_core_wedge__chiralatomrank_2147483648_wedge() {
        let g = graph("_ChiralAtomRank", 2147483648_u32);
        let before = g.clone();
        assert_eq!(
            can_be_stereo_bond(&g, &g.bonds[0], false),
            Err(WedgeError::UnsignedRankOverflow {
                atom: AtomId::new(2),
                property: "_ChiralAtomRank",
                value: 2147483648_u32
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__CIPRank_4294967295
    #[test]
    fn uint_cell_signed_consumer_core_wedge__ciprank_4294967295_wedge() {
        let g = graph("_CIPRank", 4294967295_u32);
        let before = g.clone();
        assert_eq!(
            can_be_stereo_bond(&g, &g.bonds[0], true),
            Err(WedgeError::UnsignedRankOverflow {
                atom: AtomId::new(2),
                property: "_CIPRank",
                value: 4294967295_u32
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: SIGNED_CONSUMER_core/wedge__ChiralAtomRank_4294967295
    #[test]
    fn uint_cell_signed_consumer_core_wedge__chiralatomrank_4294967295_wedge() {
        let g = graph("_ChiralAtomRank", 4294967295_u32);
        let before = g.clone();
        assert_eq!(
            can_be_stereo_bond(&g, &g.bonds[0], false),
            Err(WedgeError::UnsignedRankOverflow {
                atom: AtomId::new(2),
                property: "_ChiralAtomRank",
                value: 4294967295_u32
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: WEDGE_PROFILE_GUARD
    #[test]
    fn uint_cell_wedge_profile_guard_wedge() {
        let g = graph("_CIPRank", 4294967295);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], false), Ok(true));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: WEDGE_DIRECTION_GUARD
    #[test]
    fn uint_cell_wedge_direction_guard_wedge() {
        let mut g = graph("_CIPRank", 4294967295);
        g.bonds[1].set_direction(BondDirection::EndUpRight);
        let before = g.clone();
        assert_eq!(can_be_stereo_bond(&g, &g.bonds[0], true), Ok(false));
        assert_eq!(g, before);
    }
}

#[cfg(test)]
mod source662_explicit_picker_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, Conformer3D, SourceAtomValenceFacts};
    use cosmolkit_types::Element;
    fn graph(stereo: BondStereo, chiral: bool, dir: BondDirection) -> TopologyBlock {
        let mut atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(if i == 2 { Element::H } else { Element::C })
                        .with_no_implicit(true),
                )
            })
            .collect::<Vec<_>>();
        if chiral {
            atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        }
        let bonds = [(0, 1), (2, 0), (3, 1)]
            .into_iter()
            .enumerate()
            .map(|(i, (a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single)
                        .with_stereo(if i == 0 { stereo } else { BondStereo::None })
                        .with_direction(if i == 1 { dir } else { BondDirection::None }),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn sparse() -> RingInfo {
        RingInfo::new(RingFindType::Sssr, 0, 0)
    }
    fn weak() -> RingInfo {
        RingInfo::new(RingFindType::Fast, 4, 3)
    }
    fn conf() -> Conformer3D {
        Conformer3D::new(
            0,
            vec![[0.0; 3], [0.0; 3], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
            true,
        )
    }
    #[test]
    fn source662_no_chiral_centers_still_promote_actual_sparse_cache_and_run_atrop_stage() {
        let mut t = graph(BondStereo::AtropCcw, false, BondDirection::None);
        let mut r = weak();
        let m = pick_bonds_to_wedge_source(&mut t, None, &mut r, None).unwrap();
        assert!(r.is_sssr_or_better());
        assert_eq!(r.atom_row_count(), 0);
        assert_eq!(r.bond_row_count(), 0);
        assert!(
            matches!(m.get(BondId::new(1)),Some(WedgeInfo::Atropisomer {update}) if update.atropisomer_bond==BondId::new(0))
        );
        assert_eq!(t.bonds[1].begin(), AtomId::new(0));
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginWedge);
    }
    #[test]
    fn source662_only_chiral_source_map_entry_has_native_center_and_no_graph_wedge_write() {
        let mut t = graph(BondStereo::None, true, BondDirection::None);
        let before = t.clone();
        let mut r = weak();
        let m = pick_bonds_to_wedge_source(&mut t, None, &mut r, None).unwrap();
        assert_eq!(m.iter().count(), 1);
        assert_eq!(
            m.get(BondId::new(1)),
            Some(&WedgeInfo::Chiral {
                center: AtomId::new(0)
            })
        );
        assert_eq!(t, before);
        assert!(r.is_sssr_or_better());
    }
    #[test]
    fn source662_no_conformer_preserves_chiral_occupancy_and_adds_only_selected_atrop_info() {
        let mut t = graph(BondStereo::AtropCcw, true, BondDirection::None);
        let mut r = sparse();
        let m = pick_bonds_to_wedge_source(&mut t, None, &mut r, None).unwrap();
        assert_eq!(
            m.get(BondId::new(1)),
            Some(&WedgeInfo::Chiral {
                center: AtomId::new(0)
            })
        );
        assert!(
            matches!(m.get(BondId::new(2)),Some(WedgeInfo::Atropisomer {update}) if update.atropisomer_bond==BondId::new(0))
        );
        assert_eq!(t.bonds[1].direction(), BondDirection::None);
        assert_eq!(t.bonds[2].direction(), BondDirection::BeginDash);
    }
    #[test]
    fn source662_three_d_native_axial_key_quirk_overwrites_previously_chiral_carrier_entry() {
        let mut t = graph(BondStereo::AtropCcw, true, BondDirection::None);
        let mut r = sparse();
        let c = conf();
        let m = pick_bonds_to_wedge_source(
            &mut t,
            Some(AtropisomerConformer::ThreeD(&c)),
            &mut r,
            None,
        )
        .unwrap();
        assert_eq!(m.iter().count(), 1);
        assert!(
            matches!(m.get(BondId::new(1)),Some(WedgeInfo::Atropisomer {update}) if update.atropisomer_bond==BondId::new(0)&&update.direction==BondDirection::BeginWedge)
        );
        assert_eq!(t.bonds[2].begin(), AtomId::new(1));
        assert_eq!(t.bonds[2].direction(), BondDirection::None);
    }
    #[test]
    fn source662_existing_graph_direction_writes_do_not_become_native_map_entries() {
        let mut t = graph(BondStereo::AtropCcw, false, BondDirection::BeginDash);
        t.bonds[1].set_endpoints(AtomId::new(0), AtomId::new(2));
        let original = t.clone();
        let mut r = sparse();
        let m = pick_bonds_to_wedge_source(&mut t, None, &mut r, None).unwrap();
        assert_eq!(m.iter().count(), 0);
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginWedge);
        let projected = pick_bonds_to_wedge_with_existing_ring_info(&original, None, None).unwrap();
        assert_eq!(projected.0.iter().count(), 0);
    }
    #[test]
    fn source662_real_property_error_precedes_chiral_scoring_and_atrop_reads() {
        let mut t = graph(BondStereo::AtropCcw, true, BondDirection::None);
        t.atoms[0].set_no_implicit(false);
        t.atoms[0].set_source_valence_facts(SourceAtomValenceFacts::UNINITIALIZED);
        let before = t.clone();
        let mut r = weak();
        let mut props = cosmolkit_model::MoleculeProperties::default();
        props.set_prop("__computedProps", false).unwrap();
        assert!(matches!(
            pick_bonds_to_wedge_source(&mut t, None, &mut r, Some(&mut props)),
            Err(WedgeError::RingFinding(
                crate::RingFindingError::MoleculeProperty(_)
            ))
        ));
        assert!(r.is_sssr_or_better());
        assert_eq!(t, before);
        assert_eq!(
            props.prop("__computedProps"),
            Some(&PropertyValue::Bool(false))
        );
    }
    #[test]
    fn source662_mutable_and_query_projection_return_identical_native_maps_without_graph_clone() {
        for chiral in [false, true] {
            for use_3d in [false, true] {
                let original = graph(BondStereo::AtropCcw, chiral, BondDirection::None);
                let mut t = original.clone();
                let mut r = weak();
                let c = conf();
                let c = if use_3d {
                    Some(AtropisomerConformer::ThreeD(&c))
                } else {
                    None
                };
                let m = pick_bonds_to_wedge_source(&mut t, c, &mut r, None).unwrap();
                let (projected, pr) =
                    pick_bonds_to_wedge_with_existing_ring_info(&original, c, None).unwrap();
                assert_eq!(m, projected);
                assert_eq!(r, pr);
            }
        }
    }
    #[test]
    fn source662_checked_input_row_guard_and_native_sparse_source_entry_are_distinct_boundaries() {
        let original = graph(BondStereo::AtropCcw, false, BondDirection::None);
        assert!(matches!(
            pick_bonds_to_wedge_with_existing_ring_info(&original, None, Some(sparse())),
            Err(WedgeError::Atropisomer(
                AtropisomerError::RingAtomRowCount {
                    actual: 0,
                    expected: 4
                }
            ))
        ));
        let mut t = original;
        let m = pick_bonds_to_wedge_source(&mut t, None, &mut sparse(), None).unwrap();
        assert_eq!(m.iter().count(), 1);
    }
}

#[cfg(test)]
mod source666_default_picker_tests {
    use super::*;
    use cosmolkit_model::{
        AtomSpec, BondSpec, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension,
    };
    use cosmolkit_types::Element;
    fn graph() -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(Element::C).with_no_implicit(true),
                )
            })
            .collect();
        let bonds = [(0, 1), (2, 0), (3, 1)]
            .into_iter()
            .enumerate()
            .map(|(i, (a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single).with_stereo(
                        if i == 0 {
                            BondStereo::AtropCcw
                        } else {
                            BondStereo::None
                        },
                    ),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn weak() -> RingInfo {
        RingInfo::new(RingFindType::Fast, 4, 3)
    }
    fn two(id: usize) -> Conformer2D {
        Conformer2D::new(id, vec![[0.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [-1.0, 1.0]])
    }
    fn three(id: usize, z2: f64, z3: f64) -> Conformer3D {
        Conformer3D::new(
            id,
            vec![[0.0; 3], [0.0; 3], [0.0, 0.0, z2], [0.0, 0.0, z3]],
            true,
        )
    }
    fn pick(t: &mut TopologyBlock, c: &CoordinateBlock) -> WedgeAssignments {
        pick_bonds_to_wedge_default_source(t, c, &mut weak(), None).unwrap()
    }
    #[test]
    fn source666_empty_native_count_selects_no_conformer_without_provenance_guess() {
        let mut t = graph();
        let c = CoordinateBlock {
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            source_conformer_order: Some(Vec::new()),
            ..Default::default()
        };
        let m = pick(&mut t, &c);
        assert!(m.get(BondId::new(1)).is_some());
        assert_eq!(t.bonds[2].begin(), AtomId::new(3));
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginWedge);
    }
    #[test]
    fn source666_nonempty_default_uses_actual_front_not_smallest_or_zero_or_unique_id() {
        for ids in [[91, 0], [17, 17]] {
            let mut t = graph();
            let c = CoordinateBlock {
                conformers_3d: vec![three(ids[0], -1.0, 1.0), three(ids[1], 1.0, -1.0)],
                ..Default::default()
            };
            let m = pick(&mut t, &c);
            assert!(m.get(BondId::new(2)).is_some());
            assert_eq!(t.bonds[1].begin(), AtomId::new(0));
            assert_eq!(t.bonds[2].direction(), BondDirection::BeginWedge);
        }
    }
    #[test]
    fn source666_single_dimension_front_does_not_read_later_geometry() {
        let mut t = graph();
        let c = CoordinateBlock {
            conformers_2d: vec![two(29), Conformer2D::new(0, vec![[f64::NAN; 2]; 1])],
            ..Default::default()
        };
        let m = pick(&mut t, &c);
        assert!(m.get(BondId::new(1)).is_some());
        assert_eq!(t.bonds[2].begin(), AtomId::new(3));
    }
    #[test]
    fn source666_actual_mixed_append_order_overrides_ids_and_import_dimension_hint() {
        for first in [CoordinateDimension::TwoD, CoordinateDimension::ThreeD] {
            let mut t = graph();
            let second = if first == CoordinateDimension::TwoD {
                CoordinateDimension::ThreeD
            } else {
                CoordinateDimension::TwoD
            };
            let c = CoordinateBlock {
                conformers_2d: vec![two(0)],
                conformers_3d: vec![three(91, -1.0, 1.0)],
                source_conformer_order: Some(vec![first, second]),
                source_coordinate_dim: Some(second),
            };
            let m = pick(&mut t, &c);
            let id = if first == CoordinateDimension::TwoD {
                1
            } else {
                2
            };
            assert!(m.get(BondId::new(id)).is_some());
        }
    }
    #[test]
    fn source666_selected_source_is3d_flag_controls_native_dispatch_independent_of_storage() {
        for flag in [false, true] {
            let mut t = graph();
            let c = CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    47,
                    vec![
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [-1.0, 0.0, 10.0],
                        [-1.0, 1.0, -10.0],
                    ],
                    flag,
                )],
                ..Default::default()
            };
            let m = pick(&mut t, &c);
            assert!(m.get(BondId::new(1)).is_some());
            assert_eq!(t.bonds[2].begin(), AtomId::new(if flag { 1 } else { 3 }));
        }
    }
    #[test]
    fn source666_absent_mixed_order_is_structural_error_before_cache_or_graph_changes() {
        let mut t = graph();
        let before = t.clone();
        let c = CoordinateBlock {
            conformers_2d: vec![two(0)],
            conformers_3d: vec![three(91, -1.0, 1.0)],
            ..Default::default()
        };
        let mut r = weak();
        let before_r = r.clone();
        assert!(matches!(
            pick_bonds_to_wedge_default_source(&mut t, &c, &mut r, None),
            Err(WedgeError::Coordinates(
                CoordinateValidationError::MissingSourceConformerOrder
            ))
        ));
        assert_eq!(t, before);
        assert_eq!(r, before_r);
    }
}

#[cfg(test)]
mod complete_determine_wedge_state_tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, AtomSpec, BondSpec, Conformer2D, Conformer3D};
    use cosmolkit_types::Element;

    fn graph() -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(Element::C).with_chiral_tag(if i == 0 {
                        ChiralTag::TetrahedralCcw
                    } else {
                        ChiralTag::Unspecified
                    }),
                )
            })
            .collect();
        let bonds = (1..4)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(i), BondOrder::Single)
                        .with_direction(BondDirection::EndUpRight),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn xy() -> Vec<[f64; 2]> {
        vec![
            [0., 0.],
            [1., 0.],
            [-0.5, 3.0_f64.sqrt() / 2.],
            [-0.5, -3.0_f64.sqrt() / 2.],
        ]
    }
    fn run(
        topology: &TopologyBlock,
        positions: Vec<[f64; 2]>,
        from: usize,
    ) -> Result<BondDirection, WedgeError> {
        let conf = Conformer2D::new(9, positions);
        determine_bond_wedge_state(
            topology,
            BondId::new(0),
            AtomId::new(from),
            Some(AtropisomerConformer::TwoD(&conf)),
        )
    }
    #[test]
    fn null_conformer_returns_direction_before_endpoint_tag_and_adjacency_access() {
        let mut topology = graph();
        topology.atoms.clear();
        topology.adjacency = AdjacencyList::default();
        assert_eq!(
            determine_bond_wedge_state(&topology, BondId::new(0), AtomId::new(usize::MAX), None),
            Ok(BondDirection::EndUpRight)
        );
        topology.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        );
        assert_eq!(
            determine_bond_wedge_state(&topology, BondId::new(0), AtomId::new(0), None),
            Err(WedgeError::NonSingleBond {
                bond: BondId::new(0)
            })
        );
    }
    #[test]
    fn any_non_begin_from_index_selects_end_center_without_an_incidence_guard() {
        let mut topology = graph();
        topology.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single),
        );
        for from in [0, 2, usize::MAX] {
            assert_eq!(run(&topology, xy(), from), Ok(BondDirection::BeginWedge));
        }
    }
    #[test]
    fn source_zeros_nonfinite_z_and_ignores_conformer_dimension_flag() {
        let topology = graph();
        for flag in [false, true] {
            let conf = Conformer3D::new(
                9,
                xy().into_iter()
                    .enumerate()
                    .map(|(i, p)| {
                        [
                            p[0],
                            p[1],
                            if i % 2 == 0 { f64::NAN } else { f64::INFINITY },
                        ]
                    })
                    .collect(),
                flag,
            );
            assert_eq!(
                determine_bond_wedge_state(
                    &topology,
                    BondId::new(0),
                    AtomId::new(0),
                    Some(AtropisomerConformer::ThreeD(&conf))
                ),
                Ok(BondDirection::BeginWedge)
            );
        }
    }
    #[test]
    fn detached_unowned_conformer_and_unused_graph_rows_do_not_trigger_global_validation() {
        let mut topology = graph();
        topology
            .atoms
            .push(Atom::from_spec(AtomId::new(99), AtomSpec::new(Element::C)));
        topology.bonds.push(Bond::from_spec(
            BondId::new(99),
            BondSpec::new(AtomId::new(99), AtomId::new(100), BondOrder::Single),
        ));
        assert_eq!(run(&topology, xy(), 0), Ok(BondDirection::BeginWedge));
    }
    #[test]
    fn overlap_returns_original_before_missing_adjacency_or_later_coordinate_reads() {
        let mut topology = graph();
        topology.adjacency = AdjacencyList::default();
        assert_eq!(
            run(&topology, vec![[0., 0.], [0., 0.]], 0),
            Ok(BondDirection::EndUpRight)
        );
        assert_eq!(
            run(&topology, vec![[0., 0.], [1., 0.]], 0),
            Err(WedgeError::SourcePrecondition {
                message: "source wedge atom adjacency row is absent"
            })
        );
    }
    #[test]
    fn conformer_range_error_occurs_only_at_first_reached_position() {
        let topology = graph();
        assert_eq!(
            run(&topology, vec![], 0),
            Err(WedgeError::SourceStateIndex {
                state: "conformer position",
                index: 0,
                count: 0
            })
        );
        assert_eq!(
            run(&topology, vec![[0., 0.], [1., 0.]], 0),
            Err(WedgeError::SourceStateIndex {
                state: "conformer position",
                index: 2,
                count: 2
            })
        );
        assert_eq!(
            run(&topology, vec![[0., 0.], [1., 0.], [0., 0.]], 0),
            Ok(BondDirection::EndUpRight)
        );
    }
    #[test]
    fn real_bond_endpoints_replace_stale_adjacency_atom_cache() {
        let mut topology = graph();
        topology.bonds[1] = Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
        );
        topology.bonds[2] = Bond::from_spec(
            BondId::new(2),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
        );
        assert_eq!(run(&topology, xy(), 0), Ok(BondDirection::BeginDash));
        topology.bonds[1] = Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        );
        assert_eq!(
            run(&topology, xy(), 0),
            Err(WedgeError::CenterNotIncident {
                bond: BondId::new(1),
                center: AtomId::new(0)
            })
        );
    }
    #[test]
    fn nonfinite_xy_follows_native_nan_angle_insertion_instead_of_being_rejected() {
        let topology = graph();
        let mut positions = xy();
        positions[2][0] = f64::NAN;
        assert_eq!(run(&topology, positions, 0), Ok(BondDirection::BeginDash));
    }
}

#[cfg(test)]
mod complete_wedge_map_overload_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, Conformer2D};
    use cosmolkit_types::Element;
    fn raw(order: BondOrder, direction: BondDirection, id: usize) -> TopologyBlock {
        let mut topology = TopologyBlock::default();
        topology.bonds.push(Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(88), AtomId::new(99), order).with_direction(direction),
        ));
        topology
    }
    fn atrop(direction: BondDirection, key: usize) -> WedgeAssignments {
        let mut result = WedgeAssignments::default();
        result.by_bond.insert(
            BondId::new(key),
            WedgeInfo::Atropisomer {
                update: AtropisomerWedgeUpdate {
                    bond: BondId::new(999),
                    begin: AtomId::new(77),
                    end: AtomId::new(66),
                    atropisomer_bond: BondId::new(555),
                    direction,
                },
            },
        );
        result
    }
    #[test]
    fn missing_assignment_preserves_direction_without_accessing_graph_or_conformer() {
        let conf = Conformer2D::new(1, vec![]);
        for order in [BondOrder::Single, BondOrder::Double, BondOrder::Aromatic] {
            let topology = raw(order, BondDirection::EndDownRight, 0);
            assert_eq!(
                determine_bond_wedge_state_from_assignments(
                    &topology,
                    BondId::new(0),
                    &WedgeAssignments::default(),
                    Some(AtropisomerConformer::TwoD(&conf))
                ),
                Ok(BondDirection::EndDownRight)
            );
        }
    }
    #[test]
    fn atrop_entry_returns_every_stored_direction_without_reading_its_other_fields() {
        let topology = raw(BondOrder::Double, BondDirection::None, 0);
        let conf = Conformer2D::new(1, vec![]);
        for direction in [
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::EndDownRight,
            BondDirection::EndUpRight,
            BondDirection::EitherDouble,
            BondDirection::Unknown,
        ] {
            assert_eq!(
                determine_bond_wedge_state_from_assignments(
                    &topology,
                    BondId::new(0),
                    &atrop(direction, 0),
                    Some(AtropisomerConformer::TwoD(&conf))
                ),
                Ok(direction)
            );
        }
    }
    #[test]
    fn lookup_uses_actual_bond_metadata_id_instead_of_input_row_index() {
        let topology = raw(BondOrder::Double, BondDirection::EndUpRight, 7);
        assert_eq!(
            determine_bond_wedge_state_from_assignments(
                &topology,
                BondId::new(0),
                &atrop(BondDirection::Unknown, 7),
                None
            ),
            Ok(BondDirection::Unknown)
        );
        assert_eq!(
            determine_bond_wedge_state_from_assignments(
                &topology,
                BondId::new(0),
                &atrop(BondDirection::Unknown, 0),
                None
            ),
            Ok(BondDirection::EndUpRight)
        );
    }
    #[test]
    fn chiral_entry_delegates_order_precondition_and_null_conformer_short_circuit() {
        let mut assignments = WedgeAssignments::default();
        assignments.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Chiral {
                center: AtomId::new(999),
            },
        );
        assert_eq!(
            determine_bond_wedge_state_from_assignments(
                &raw(BondOrder::Single, BondDirection::Unknown, 0),
                BondId::new(0),
                &assignments,
                None
            ),
            Ok(BondDirection::Unknown)
        );
        assert_eq!(
            determine_bond_wedge_state_from_assignments(
                &raw(BondOrder::Double, BondDirection::Unknown, 0),
                BondId::new(0),
                &assignments,
                None
            ),
            Err(WedgeError::NonSingleBond {
                bond: BondId::new(0)
            })
        );
    }
    #[test]
    fn chiral_entry_delegates_non_begin_from_index_to_source_end_center() {
        let atoms = (0..2)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(Element::C).with_chiral_tag(if i == 1 {
                        ChiralTag::TetrahedralCcw
                    } else {
                        ChiralTag::Unspecified
                    }),
                )
            })
            .collect();
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let mut assignments = WedgeAssignments::default();
        assignments.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Chiral {
                center: AtomId::new(999),
            },
        );
        let conf = Conformer2D::new(1, vec![[0., 0.], [1., 0.]]);
        assert_eq!(
            determine_bond_wedge_state_from_assignments(
                &topology,
                BondId::new(0),
                &assignments,
                Some(AtropisomerConformer::TwoD(&conf))
            ),
            Ok(BondDirection::BeginWedge)
        );
        assert_eq!(
            determine_bond_wedge_state_from_assignments(
                &topology,
                BondId::new(99),
                &assignments,
                None
            ),
            Err(WedgeError::BondOutOfRange {
                bond: BondId::new(99),
                bond_count: 1
            })
        );
    }
}

#[cfg(test)]
mod complete_molfile_direction_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(order: BondOrder, direction: BondDirection, no_implicit: bool) -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(Element::C).with_no_implicit(no_implicit),
                )
            })
            .collect();
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_direction(direction),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(1), AtomId::new(3), BondOrder::Single),
            ),
        ];
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn empty_rings() -> RingInfo {
        let mut rings = RingInfo::new(RingFindType::OtherOrUnknown, 0, 0);
        rings.reset();
        rings
    }
    fn valence() -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: vec![3, 3, 1, 1],
            implicit_hydrogens: vec![0; 4],
        }
    }
    fn chiral(key: usize, center: usize) -> WedgeAssignments {
        let mut assignments = WedgeAssignments::default();
        assignments.by_bond.insert(
            BondId::new(key),
            WedgeInfo::Chiral {
                center: AtomId::new(center),
            },
        );
        assignments
    }
    fn call(
        topology: &TopologyBlock,
        valence: &ValenceAssignment,
        rings: &RingInfo,
        map: &WedgeAssignments,
        id: usize,
        direction: &mut BondDirection,
        reverse: &mut bool,
    ) -> Result<(), WedgeError> {
        let ctx = CrossedBondContext {
            topology,
            valence,
            rings,
            use_legacy_stereo_perception: true,
        };
        get_molfile_bond_stereo_direction_source(
            &ctx,
            map,
            BondId::new(id),
            None,
            direction,
            reverse,
        )
    }
    #[test]
    fn bond_precondition_precedes_both_caller_output_resets() {
        let topology = graph(BondOrder::Single, BondDirection::None, false);
        let mut dir = BondDirection::BeginDash;
        let mut reverse = true;
        assert_eq!(
            call(
                &topology,
                &valence(),
                &empty_rings(),
                &WedgeAssignments::default(),
                99,
                &mut dir,
                &mut reverse
            ),
            Err(WedgeError::BondOutOfRange {
                bond: BondId::new(99),
                bond_count: 3
            })
        );
        assert_eq!((dir, reverse), (BondDirection::BeginDash, true));
    }
    #[test]
    fn child_error_exposes_already_reset_direction_and_reverse() {
        let topology = graph(BondOrder::Aromatic, BondDirection::Unknown, false);
        let mut dir = BondDirection::BeginDash;
        let mut reverse = true;
        assert_eq!(
            call(
                &topology,
                &valence(),
                &empty_rings(),
                &chiral(0, 1),
                0,
                &mut dir,
                &mut reverse
            ),
            Err(WedgeError::NonSingleBond {
                bond: BondId::new(0)
            })
        );
        assert_eq!((dir, reverse), (BondDirection::None, false));
    }
    #[test]
    fn reverse_is_limited_to_three_direction_tags_and_chiral_info_type() {
        for direction in [
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::EndDownRight,
            BondDirection::EndUpRight,
            BondDirection::EitherDouble,
            BondDirection::Unknown,
        ] {
            let topology = graph(BondOrder::Single, direction, false);
            let mut dir = BondDirection::None;
            let mut reverse = false;
            call(
                &topology,
                &valence(),
                &empty_rings(),
                &chiral(0, 1),
                0,
                &mut dir,
                &mut reverse,
            )
            .unwrap();
            assert_eq!(dir, direction);
            assert_eq!(
                reverse,
                matches!(
                    direction,
                    BondDirection::BeginWedge | BondDirection::BeginDash | BondDirection::Unknown
                )
            );
            let mut map = WedgeAssignments::default();
            map.by_bond.insert(
                BondId::new(0),
                WedgeInfo::Atropisomer {
                    update: AtropisomerWedgeUpdate {
                        bond: BondId::new(99),
                        begin: AtomId::new(1),
                        end: AtomId::new(0),
                        atropisomer_bond: BondId::new(77),
                        direction,
                    },
                },
            );
            call(
                &topology,
                &valence(),
                &empty_rings(),
                &map,
                0,
                &mut dir,
                &mut reverse,
            )
            .unwrap();
            assert_eq!((dir, reverse), (direction, false));
        }
    }
    #[test]
    fn reverse_lookup_uses_the_same_actual_object_id_as_direction_lookup() {
        let mut topology = graph(BondOrder::Single, BondDirection::Unknown, false);
        topology.bonds[0] = Bond::from_spec(
            BondId::new(7),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_direction(BondDirection::Unknown),
        );
        let mut dir = BondDirection::None;
        let mut reverse = false;
        call(
            &topology,
            &valence(),
            &empty_rings(),
            &chiral(7, 1),
            0,
            &mut dir,
            &mut reverse,
        )
        .unwrap();
        assert_eq!((dir, reverse), (BondDirection::Unknown, true));
        call(
            &topology,
            &valence(),
            &empty_rings(),
            &chiral(0, 1),
            0,
            &mut dir,
            &mut reverse,
        )
        .unwrap();
        assert_eq!((dir, reverse), (BondDirection::Unknown, false));
    }
    #[test]
    fn unrelated_orders_ignore_assignment_and_return_reset_outputs() {
        let mut dir = BondDirection::BeginDash;
        let mut reverse = true;
        for order in [BondOrder::Triple, BondOrder::Quadruple] {
            let topology = graph(order, BondDirection::Unknown, false);
            call(
                &topology,
                &ValenceAssignment {
                    explicit_valence: vec![],
                    implicit_hydrogens: vec![],
                },
                &empty_rings(),
                &chiral(0, 1),
                0,
                &mut dir,
                &mut reverse,
            )
            .unwrap();
            assert_eq!((dir, reverse), (BondDirection::None, false));
        }
    }
    #[test]
    fn potential_bond_reads_actual_ring_vectors_without_member_row_or_initialization_guards() {
        let topology = graph(BondOrder::Double, BondDirection::None, false);
        let mut dir = BondDirection::None;
        let mut reverse = true;
        for rings in [
            empty_rings(),
            RingInfo::new(RingFindType::Fast, 0, 0),
            RingInfo::new(RingFindType::Sssr, 99, 1),
        ] {
            assert!(
                potential_stereo::is_potential_bond_source(
                    &topology,
                    &valence(),
                    &rings,
                    &topology.bonds[0]
                )
                .unwrap()
            );
            call(
                &topology,
                &valence(),
                &rings,
                &WedgeAssignments::default(),
                0,
                &mut dir,
                &mut reverse,
            )
            .unwrap();
            assert_eq!((dir, reverse), (BondDirection::EitherDouble, false));
        }
    }
    #[test]
    fn begin_unsaturation_false_short_circuits_missing_end_explicit_cache() {
        let topology = graph(BondOrder::Double, BondDirection::None, false);
        let mut cache = valence();
        cache.explicit_valence = vec![4];
        let mut dir = BondDirection::Unknown;
        let mut reverse = true;
        call(
            &topology,
            &cache,
            &empty_rings(),
            &WedgeAssignments::default(),
            0,
            &mut dir,
            &mut reverse,
        )
        .unwrap();
        assert_eq!((dir, reverse), (BondDirection::None, false));
        cache.explicit_valence[0] = 3;
        assert!(matches!(
            call(
                &topology,
                &cache,
                &empty_rings(),
                &WedgeAssignments::default(),
                0,
                &mut dir,
                &mut reverse
            ),
            Err(WedgeError::PotentialStereo(
                PotentialStereoError::InvalidValence {
                    field: "explicit_valence",
                    ..
                }
            ))
        ));
        assert_eq!((dir, reverse), (BondDirection::None, false));
    }
    #[test]
    fn no_implicit_short_circuits_absent_cache_in_degree_hydrogens_and_total_valence() {
        let topology = graph(BondOrder::Double, BondDirection::None, true);
        let mut cache = valence();
        cache.implicit_hydrogens.clear();
        let mut dir = BondDirection::None;
        let mut reverse = true;
        assert_eq!(
            potential_stereo::total_degree(&topology, &cache, AtomId::new(0)).unwrap(),
            2
        );
        assert_eq!(
            source_total_valence(&cache, AtomId::new(0), &topology).unwrap(),
            3
        );
        call(
            &topology,
            &cache,
            &empty_rings(),
            &WedgeAssignments::default(),
            0,
            &mut dir,
            &mut reverse,
        )
        .unwrap();
        assert_eq!((dir, reverse), (BondDirection::EitherDouble, false));
    }
}

#[cfg(test)]
mod complete_molfile_code_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(order: BondOrder) -> TopologyBlock {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), order)
                    .with_direction(BondDirection::Unknown),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(1), AtomId::new(3), BondOrder::Single),
            ),
        ];
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn call(
        topology: &TopologyBlock,
        cache: &ValenceAssignment,
        map: &WedgeAssignments,
        id: usize,
        code: &mut i32,
        reverse: &mut bool,
    ) -> Result<BondDirection, WedgeError> {
        let rings = RingInfo::new(RingFindType::Fast, 0, 0);
        let ctx = CrossedBondContext {
            topology,
            valence: cache,
            rings: &rings,
            use_legacy_stereo_perception: true,
        };
        get_molfile_bond_stereo_code_source(&ctx, map, BondId::new(id), None, code, reverse)
    }
    fn cache() -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: vec![3, 3, 1, 1],
            implicit_hydrogens: vec![0; 4],
        }
    }
    fn chiral() -> WedgeAssignments {
        let mut map = WedgeAssignments::default();
        map.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Chiral {
                center: AtomId::new(1),
            },
        );
        map
    }
    #[test]
    fn bond_precondition_failure_preserves_both_external_output_parameters() {
        let mut code = 99;
        let mut reverse = true;
        assert_eq!(
            call(
                &graph(BondOrder::Single),
                &cache(),
                &chiral(),
                99,
                &mut code,
                &mut reverse
            ),
            Err(WedgeError::BondOutOfRange {
                bond: BondId::new(99),
                bond_count: 3
            })
        );
        assert_eq!((code, reverse), (99, true));
    }
    #[test]
    fn direction_child_error_preserves_code_but_exposes_reverse_reset() {
        let mut code = 99;
        let mut reverse = true;
        assert_eq!(
            call(
                &graph(BondOrder::Aromatic),
                &cache(),
                &chiral(),
                0,
                &mut code,
                &mut reverse
            ),
            Err(WedgeError::NonSingleBond {
                bond: BondId::new(0)
            })
        );
        assert_eq!((code, reverse), (99, false));
    }
    #[test]
    fn crossed_child_late_error_never_converts_a_partial_direction_into_a_code() {
        let mut code = 99;
        let mut reverse = true;
        let mut data = cache();
        data.explicit_valence = vec![3];
        assert!(matches!(
            call(
                &graph(BondOrder::Double),
                &data,
                &WedgeAssignments::default(),
                0,
                &mut code,
                &mut reverse
            ),
            Err(WedgeError::PotentialStereo(
                PotentialStereoError::InvalidValence {
                    field: "explicit_valence",
                    ..
                }
            ))
        ));
        assert_eq!((code, reverse), (99, false));
    }
    #[test]
    fn successful_chiral_unknown_returns_direction_and_assigns_code_after_reverse_selection() {
        let mut code = 99;
        let mut reverse = false;
        assert_eq!(
            call(
                &graph(BondOrder::Single),
                &cache(),
                &chiral(),
                0,
                &mut code,
                &mut reverse
            ),
            Ok(BondDirection::Unknown)
        );
        assert_eq!((code, reverse), (4, true));
    }
    #[test]
    fn crossed_and_unrelated_order_successes_replace_both_prior_outputs() {
        for (order, expected, expected_code) in [
            (BondOrder::Double, BondDirection::EitherDouble, 3),
            (BondOrder::Triple, BondDirection::None, 0),
        ] {
            let mut code = 99;
            let mut reverse = true;
            assert_eq!(
                call(
                    &graph(order),
                    &cache(),
                    &WedgeAssignments::default(),
                    0,
                    &mut code,
                    &mut reverse
                ),
                Ok(expected)
            );
            assert_eq!((code, reverse), (expected_code, false));
        }
    }
}

fn get_directional_bond_stereo_source<G: StereoGraphAccess>(
    topology: &G,
    wedge_assignments: &WedgeAssignments,
    bond_id: BondId,
    conformer: Option<AtropisomerConformer<'_>>,
    direction: &mut BondDirection,
    reverse: &mut bool,
) -> Result<(), WedgeError> {
    // RDKit❗✔️:   if (canHaveDirection(*bond)) {
    // RDKit❗✔️:     // single bond stereo chemistry
    // RDKit❗✔️:
    // RDKit❗✔️:     dir = Chirality::detail::determineBondWedgeState(bond, wedgeBonds, conf);
    // RDKit❗✔️:
    // RDKit❗✔️:     // if this bond needs to be wedged it is possible that this
    // RDKit❗✔️:     // wedging was determined by a chiral atom at the end of the
    // RDKit❗✔️:     // bond (instead of at the beginning). In this case we need to
    // RDKit❗✔️:     // reverse the begin and end atoms for the bond when we write
    // RDKit❗✔️:     // the mol file
    // RDKit❗✔️:     if ((dir == Bond::BEGINDASH) ||
    // RDKit❗✔️:         (dir == Bond::BEGINWEDGE || dir == Bond::UNKNOWN)) {
    // RDKit❗✔️:       auto wbi = wedgeBonds.find(bond->getIdx());
    // RDKit❗✔️:       if (wbi != wedgeBonds.end() &&
    // RDKit❗✔️:           wbi->second->getType() ==
    // RDKit❗✔️:               Chirality::WedgeInfoType::WedgeInfoTypeChiral &&
    // RDKit❗✔️:           static_cast<unsigned int>(wbi->second->getIdx()) !=
    // RDKit❗✔️:               bond->getBeginAtomIdx()) {
    // RDKit❗✔️:         reverse = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    let bond = topology
        .bonds()
        .get(bond_id.index())
        .ok_or(WedgeError::BondOutOfRange {
            bond: bond_id,
            bond_count: topology.bonds().len(),
        })?;
    *direction = determine_bond_wedge_state_from_assignments(
        topology,
        bond_id,
        wedge_assignments,
        conformer,
    )?;
    if matches!(
        *direction,
        BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown
    ) {
        if let Some(WedgeInfo::Chiral { center }) = wedge_assignments.by_bond.get(&bond.id()) {
            if center.index() as u32 != bond.begin().index() as u32 {
                *reverse = true;
            }
        }
    }
    Ok(())
}

/// The source directional branch selected by CX canHaveDirection. Double-bond
/// crossing is deliberately not part of this narrow borrowed query entry.
#[doc(hidden)]
pub fn get_query_directional_bond_stereo_info_source(
    query: &cosmolkit_model::QueryGraph,
    wedge_assignments: &WedgeAssignments,
    bond_id: BondId,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<MolFileBondStereoInfo, WedgeError> {
    // RDKit❗✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit❗✔️:   auto bondType = bond.getBondType();
    // RDKit❗✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit❗✔️: }
    // The selected source call reaches only single/aromatic bonds. Keep this
    // precondition explicit instead of manufacturing a crossed-bond context.
    let bond = query
        .bonds()
        .get(bond_id.index())
        .ok_or(WedgeError::BondOutOfRange {
            bond: bond_id,
            bond_count: query.num_bonds(),
        })?
        .bond();
    if !matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic) {
        return Err(WedgeError::NonSingleBond { bond: bond_id });
    }
    let mut direction = BondDirection::None;
    let mut reverse = false;
    get_directional_bond_stereo_source(
        query,
        wedge_assignments,
        bond_id,
        conformer,
        &mut direction,
        &mut reverse,
    )?;
    Ok(MolFileBondStereoInfo {
        direction,
        direction_code: bond_get_dir_code(direction),
        reverse,
    })
}

#[cfg(test)]
mod cx_query_bond_config_dependency_source_tests {
    use super::*;

    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, Conformer2D, QueryAtom, QueryAtomIdentity, QueryBond,
        QueryGraph, QueryNode,
    };
    fn graph(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn query_direction_missing_map_reads_actual_direction_without_conformer() {
        let mut q = graph(2, &[(0, 1)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::BeginWedge);
        let before = q.clone();
        let info = get_query_directional_bond_stereo_info_source(
            &q,
            &WedgeAssignments::default(),
            BondId::new(0),
            None,
        )
        .unwrap();
        assert_eq!(
            info,
            MolFileBondStereoInfo {
                direction: BondDirection::BeginWedge,
                direction_code: 1,
                reverse: false
            }
        );
        assert_eq!(q, before);
    }
    #[test]
    fn query_atrop_map_branch_does_not_validate_chiral_center_or_read_coordinates() {
        let q = graph(2, &[(0, 1)]);
        let mut w = WedgeAssignments::default();
        w.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Atropisomer {
                update: AtropisomerWedgeUpdate {
                    bond: BondId::new(0),
                    begin: AtomId::new(1),
                    end: AtomId::new(0),
                    direction: BondDirection::BeginDash,
                    atropisomer_bond: BondId::new(7),
                },
            },
        );
        let info =
            get_query_directional_bond_stereo_info_source(&q, &w, BondId::new(0), None).unwrap();
        assert_eq!(
            info,
            MolFileBondStereoInfo {
                direction: BondDirection::BeginDash,
                direction_code: 6,
                reverse: false
            }
        );
    }
    #[test]
    fn query_chiral_center_at_end_runs_same_geometry_kernel_and_reverses_endpoints() {
        let mut q = graph(4, &[(0, 1), (1, 2), (1, 3)]);
        q.atoms_mut()[1].set_chiral_tag(ChiralTag::TetrahedralCcw);
        let mut w = WedgeAssignments::default();
        w.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Chiral {
                center: AtomId::new(1),
            },
        );
        let conf = Conformer2D::new(
            0,
            vec![[1.0, 0.0], [0.0, 0.0], [-0.5, 0.866], [-0.5, -0.866]],
        );
        let before = q.clone();
        let info = get_query_directional_bond_stereo_info_source(
            &q,
            &w,
            BondId::new(0),
            Some(AtropisomerConformer::TwoD(&conf)),
        )
        .unwrap();
        assert_eq!(
            info,
            MolFileBondStereoInfo {
                direction: BondDirection::BeginWedge,
                direction_code: 1,
                reverse: true
            }
        );
        assert_eq!(q, before);
    }
    #[test]
    fn query_chiral_null_conformer_return_precedes_center_and_coordinate_reads() {
        let mut q = graph(2, &[(0, 1)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::Unknown);
        let mut w = WedgeAssignments::default();
        w.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Chiral {
                center: AtomId::new(99),
            },
        );
        let info =
            get_query_directional_bond_stereo_info_source(&q, &w, BondId::new(0), None).unwrap();
        assert_eq!(info.direction_code, 4);
        assert!(info.reverse);
    }
    #[test]
    fn query_coordinate_and_center_errors_propagate_from_canonical_kernel() {
        let mut q = graph(2, &[(0, 1)]);
        let mut w = WedgeAssignments::default();
        w.by_bond.insert(
            BondId::new(0),
            WedgeInfo::Chiral {
                center: AtomId::new(0),
            },
        );
        let conf = Conformer2D::new(0, vec![]);
        assert!(matches!(
            get_query_directional_bond_stereo_info_source(
                &q,
                &w,
                BondId::new(0),
                Some(AtropisomerConformer::TwoD(&conf))
            ),
            Err(WedgeError::UnsupportedCenterTag { .. })
        ));
        q.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCcw);
        assert!(matches!(
            get_query_directional_bond_stereo_info_source(
                &q,
                &w,
                BondId::new(0),
                Some(AtropisomerConformer::TwoD(&conf))
            ),
            Err(WedgeError::SourceStateIndex {
                state: "conformer position",
                ..
            })
        ));
    }
}

/// Complete default-parameter source wedge selection over actual QueryAtom rows.
#[doc(hidden)]
pub fn pick_query_bonds_to_wedge_source(
    query: &mut cosmolkit_model::QueryGraph,
    conformer: Option<AtropisomerConformer<'_>>,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
) -> Result<WedgeAssignments, WedgeError> {
    pick_bonds_to_wedge_kernel(
        &mut SourceWedgeGraph::Mutable(query),
        conformer,
        rings,
        properties,
    )
}

#[cfg(test)]
mod query_cx_producers_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, MoleculeProperties, PropertyText, QueryAtom,
        QueryAtomIdentity, QueryBond, QueryGraph, QueryNode, SourceAtomValenceFacts,
    };
    fn graph(count: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..count)
                .map(|i| {
                    let mut a = QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    );
                    a.set_source_valence_facts(SourceAtomValenceFacts {
                        explicit_valence: 0,
                        implicit_valence: 0,
                    });
                    a
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn rings() -> RingInfo {
        RingInfo::from_source_snapshot(&cosmolkit_model::SourceRingInfo::default()).unwrap()
    }
    #[test]
    fn empty_query_still_acquires_source_sssr_cache() {
        let mut q = graph(0, &[]);
        let mut r = rings();
        assert!(
            pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None)
                .unwrap()
                .by_bond
                .is_empty()
        );
        assert!(r.is_sssr_or_better());
        assert_eq!(r.num_rings(), 0);
    }
    #[test]
    fn uninitialized_cycle_is_found_by_the_same_borrowed_ring_kernel() {
        let mut q = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let mut r = rings();
        let before = q.clone();
        pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).unwrap();
        assert_eq!(r.num_rings(), 1);
        assert_eq!(r.min_bond_ring_size(BondId::new(0)), 3);
        assert_eq!(q, before);
    }
    #[test]
    fn trusted_empty_sssr_is_not_replaced_by_inferred_cycle() {
        let mut q = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let mut r = RingInfo::new(RingFindType::Sssr, 0, 0);
        pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).unwrap();
        assert_eq!(r.num_rings(), 0);
        assert_eq!(r.atom_row_count(), 0);
    }
    #[test]
    fn source_atom_numbers_zero_and_119_drive_chiral_candidate_scores() {
        let mut q = graph(4, &[(1, 0), (1, 2), (1, 3)]);
        q.atoms_mut()[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        let mut r = rings();
        let before = q.clone();
        let w = pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).unwrap();
        assert!(
            matches!(w.get(BondId::new(0)),Some(WedgeInfo::Chiral{center})if *center==AtomId::new(1))
        );
        assert_eq!(q, before);
    }
    #[test]
    fn already_wedged_chiral_atom_skips_source_candidate_assignment() {
        let mut q = graph(4, &[(1, 0), (1, 2), (1, 3)]);
        q.atoms_mut()[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        q.bonds_mut()[0]
            .bond_mut()
            .set_direction(BondDirection::Unknown);
        let mut r = rings();
        let w = pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).unwrap();
        assert!(w.by_bond.is_empty());
        assert!(r.is_sssr_or_better());
    }
    #[test]
    fn source_ring_property_clearing_uses_actual_computed_dictionary() {
        let mut q = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let mut r = rings();
        let mut p = MoleculeProperties::default();
        p.set_prop(
            "extraRings",
            PropertyValue::String(PropertyText::from("stale")),
        )
        .unwrap();
        p.set_prop("keep", PropertyValue::Int(7)).unwrap();
        pick_query_bonds_to_wedge_source(&mut q, None, &mut r, Some(&mut p)).unwrap();
        assert!(p.prop("extraRings").is_none());
        assert_eq!(p.prop("keep"), Some(&PropertyValue::Int(7)));
    }
    #[test]
    fn no_conformer_atrop_writes_real_query_bond_direction_and_orientation() {
        let mut q = graph(6, &[(0, 1), (2, 0), (0, 3), (4, 1), (1, 5)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        let atoms_before = q.atoms().to_vec();
        let mut r = rings();
        let w = pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).unwrap();
        assert!(!w.by_bond.is_empty());
        for (id, wi) in &w.by_bond {
            if let WedgeInfo::Atropisomer { update } = wi {
                let b = q.bonds()[id.index()].bond();
                assert_eq!(b.begin(), update.begin);
                assert_eq!(b.end(), update.end);
                assert_eq!(b.direction(), update.direction);
            }
        }
        assert_eq!(q.atoms(), atoms_before);
        assert!(r.is_sssr_or_better());
    }
    #[test]
    fn direct_atrop_query_entry_and_complete_picker_share_source_kernel() {
        let mut a = graph(6, &[(0, 1), (2, 0), (0, 3), (4, 1), (1, 5)]);
        a.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCcw)
            .unwrap();
        let mut b = a.clone();
        let mut ar = rings();
        let mut br = rings();
        let aw = pick_query_bonds_to_wedge_source(&mut a, None, &mut ar, None).unwrap();
        let mut bw = WedgeAssignments::default();
        crate::atropisomer::wedge_query_bonds_from_atropisomers_source(
            &mut b, &mut br, None, None, &mut bw,
        )
        .unwrap();
        assert_eq!(a, b);
        assert_eq!(aw, bw);
        assert_eq!(ar.into_source_snapshot(), br.into_source_snapshot());
    }
    #[test]
    fn actual_implicit_cache_error_is_propagated_after_source_ring_acquisition() {
        let mut q = graph(6, &[(0, 1), (2, 0), (0, 3), (4, 1), (1, 5)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        q.atoms_mut()[0].set_source_valence_facts(SourceAtomValenceFacts::UNINITIALIZED);
        let mut r = rings();
        assert!(pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).is_err());
        assert!(r.is_sssr_or_better());
    }
    #[test]
    fn no_implicit_short_circuits_uninitialized_actual_cache() {
        let mut q = graph(6, &[(0, 1), (2, 0), (0, 3), (4, 1), (1, 5)]);
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        for a in q.atoms_mut() {
            a.set_source_valence_facts(SourceAtomValenceFacts::UNINITIALIZED);
            a.set_no_implicit(true);
        }
        let mut r = rings();
        assert!(pick_query_bonds_to_wedge_source(&mut q, None, &mut r, None).is_ok());
    }
}
