use crate::{
    AdjacencyList, Atom, AtomId, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag,
    Molecule, NeighborRef, SubstanceGroup,
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AddHydrogensError {
    #[error("unsupported atom for hydrogen addition at {atom} with atomic number {atomic_number}")]
    UnsupportedAtom {
        atom: crate::AtomId,
        atomic_number: u8,
    },
    #[error("hydrogen addition is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RemoveHydrogensError {
    #[error("unsupported hydrogen removal at atom {atom}")]
    UnsupportedHydrogen { atom: crate::AtomId },
    #[error("removeHs atom index {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("removeHs bond index {bond} is out of range for {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
    #[error("substance-group bond {bond} is neither crossing nor contained")]
    InvalidSubstanceGroupBondRole { bond: BondId },
    #[error("incomplete RDKit removeHs port for {branch}: {reason}")]
    ProtocolDebt {
        branch: &'static str,
        reason: &'static str,
    },
    #[error(transparent)]
    Sanitize(#[from] crate::SanitizeError),
    #[error("hydrogen removal is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

// BEGIN RDKIT CPP STRUCT MolOps::RemoveHsParameters
// RDKit✔️✔️: struct RDKIT_GRAPHMOL_EXPORT RemoveHsParameters {
// RDKit✔️✔️:   bool removeDegreeZero = false;    /**< hydrogens that have no bonds */
// RDKit✔️✔️:   bool removeHigherDegrees = false; /**< hydrogens with two (or more) bonds */
// RDKit✔️✔️:   bool removeOnlyHNeighbors =
// RDKit✔️✔️:       false; /**< hydrogens with bonds only to other hydrogens */
// RDKit✔️✔️:   bool removeIsotopes = false; /**< hydrogens with non-default isotopes */
// RDKit✔️✔️:   bool removeAndTrackIsotopes = false; /**< removes hydrogens with non-default
// RDKit✔️✔️:    isotopes and keeps track of the heavy atom the isotopes were attached to in
// RDKit✔️✔️:    the private _isotopicHs atom property, so they are re-added by AddHs() as
// RDKit✔️✔️:    the original isotopes if possible*/
// RDKit✔️✔️:   bool removeDummyNeighbors =
// RDKit✔️✔️:       false; /**< hydrogens with at least one dummy-atom neighbor */
// RDKit✔️✔️:   bool removeDefiningBondStereo =
// RDKit✔️✔️:       false; /**< hydrogens defining bond stereochemistry */
// RDKit✔️✔️:   bool removeWithWedgedBond = true; /**< hydrogens with wedged bonds to them */
// RDKit✔️✔️:   bool removeWithQuery = false;     /**< hydrogens with queries defined */
// RDKit✔️✔️:   bool removeMapped = true;         /**< mapped hydrogens */
// RDKit✔️✔️:   bool removeInSGroups = true;      /**< part of a SubstanceGroup.
// RDKit✔️✔️:     An H atom will only be removed if it doesn't cause any SGroup to become empty,
// RDKit✔️✔️:     and if it doesn't play a special role in the SGroup (XBOND, attach point
// RDKit✔️✔️:     or a CState) */
// RDKit✔️✔️:   bool showWarnings = true; /**< display warnings for Hs that are not removed */
// RDKit✔️✔️:   bool removeNonimplicit = true; /**< DEPRECATED equivalent of !implicitOnly */
// RDKit✔️✔️:   bool updateExplicitCount =
// RDKit✔️✔️:       false; /**< DEPRECATED equivalent of updateExplicitCount */
// RDKit✔️✔️:   bool removeHydrides = false; /**< Removing Hydrides */
// RDKit✔️✔️:   bool removeNontetrahedralNeighbors =
// RDKit✔️✔️:       false; /**<  remove Hs which are bonded to atoms with specified
// RDKit✔️✔️:                 non-tetrahedral stereochemistry */
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT MolOps::RemoveHsParameters
#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct RemoveHsParams {
    pub(crate) remove_degree_zero: bool,
    pub(crate) remove_higher_degrees: bool,
    pub(crate) remove_only_h_neighbors: bool,
    pub(crate) remove_isotopes: bool,
    pub(crate) remove_and_track_isotopes: bool,
    pub(crate) remove_dummy_neighbors: bool,
    pub(crate) remove_defining_bond_stereo: bool,
    pub(crate) remove_with_wedged_bond: bool,
    pub(crate) remove_with_query: bool,
    pub(crate) remove_mapped: bool,
    pub(crate) remove_in_sgroups: bool,
    pub(crate) show_warnings: bool,
    pub(crate) remove_nonimplicit: bool,
    pub(crate) update_explicit_count: bool,
    pub(crate) remove_hydrides: bool,
    pub(crate) remove_nontetrahedral_neighbors: bool,
}

impl Default for RemoveHsParams {
    fn default() -> Self {
        Self {
            remove_degree_zero: false,
            remove_higher_degrees: false,
            remove_only_h_neighbors: false,
            remove_isotopes: false,
            remove_and_track_isotopes: false,
            remove_dummy_neighbors: false,
            remove_defining_bond_stereo: false,
            remove_with_wedged_bond: true,
            remove_with_query: false,
            remove_mapped: true,
            remove_in_sgroups: true,
            show_warnings: true,
            remove_nonimplicit: true,
            update_explicit_count: false,
            remove_hydrides: false,
            remove_nontetrahedral_neighbors: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub(crate) struct RemoveHsAssignment {
    pub(crate) atoms_to_remove: Vec<AtomId>,
    pub(crate) atom_explicit_hydrogen_updates: Vec<AtomExplicitHydrogenUpdate>,
    pub(crate) atom_chirality_inversions: Vec<AtomId>,
    pub(crate) atom_property_updates: Vec<AtomPropertyUpdate>,
    pub(crate) bond_direction_updates: Vec<BondDirectionUpdate>,
    pub(crate) bond_stereo_updates: Vec<BondStereoUpdate>,
    pub(crate) bond_stereo_atom_replacements: Vec<BondStereoAtomReplacement>,
    pub(crate) sgroup_updates: Vec<SGroupRemoveHsUpdate>,
    pub(crate) sanitize_after_removal: bool,
    pub(crate) clear_computed_properties: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AtomExplicitHydrogenUpdate {
    pub(crate) atom: AtomId,
    pub(crate) explicit_hydrogens: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomPropertyUpdate {
    SetUnknownStereo { atom: AtomId },
    SetIsotopicHydrogens { atom: AtomId, isotopes: Vec<u16> },
    ClearExcessChiralExplicitHydrogens { atom: AtomId },
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct BondDirectionUpdate {
    pub(crate) bond: BondId,
    pub(crate) direction: BondDirection,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct BondStereoUpdate {
    pub(crate) bond: BondId,
    pub(crate) stereo: BondStereo,
    pub(crate) stereo_atoms: Option<[AtomId; 2]>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct BondStereoAtomReplacement {
    pub(crate) bond: BondId,
    pub(crate) old_atom: AtomId,
    pub(crate) new_atom: AtomId,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum SGroupRemoveHsUpdate {
    RemoveBond { bond: BondId },
    RemoveAtom { atom: AtomId },
    RemoveParentAtom { atom: AtomId },
    ClearAttachPointLeavingAtom { atom: AtomId },
}

pub fn with_hydrogens(_molecule: &Molecule) -> Result<Molecule, AddHydrogensError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::HYDROGENS_FEATURE).into())
}

pub fn without_hydrogens(_molecule: &Molecule) -> Result<Molecule, RemoveHydrogensError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::HYDROGENS_FEATURE).into())
}

pub fn without_hydrogens_with_sanitize(
    _molecule: &Molecule,
    _sanitize: bool,
) -> Result<Molecule, RemoveHydrogensError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::HYDROGENS_FEATURE).into())
}

#[allow(dead_code)]
pub(crate) fn remove_hs_assignment(
    molecule: &Molecule,
    params: &RemoveHsParams,
    sanitize: bool,
) -> Result<RemoveHsAssignment, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs(RWMol&, const RemoveHsParameters&, bool)
    // RDKit❌❌: void removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize) {
    // RDKit❌❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit❌❌:     // if there are any non-isotopic Hs remove them first
    // RDKit❌❌:     // to make sure chirality is preserved
    // RDKit❌❌:     bool needRemoveHs = false;
    // RDKit❌❌:     for (auto atom : mol.atoms()) {
    // RDKit❌❌:       if (atom->getAtomicNum() == 1 && atom->getIsotope() == 0) {
    // RDKit❌❌:         needRemoveHs = true;
    // RDKit❌❌:         break;
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:     if (needRemoveHs) {
    // RDKit❌❌:       RemoveHsParameters psCopy(ps);
    // RDKit❌❌:       psCopy.removeAndTrackIsotopes = false;
    // RDKit❌❌:       psCopy.removeIsotopes = false;
    // RDKit❌❌:       removeHs(mol, psCopy, false);
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌:   for (auto atom : mol.atoms()) {
    // RDKit❌❌:     atom->updatePropertyCache(false);
    // RDKit❌❌:   }
    // RDKit❌❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit❌❌:     for (const auto &pair : getIsoMap(mol)) {
    // RDKit❌❌:       mol.getAtomWithIdx(pair.first)
    // RDKit❌❌:           ->setProp(common_properties::_isotopicHs, pair.second);
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌:   boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
    // RDKit❌❌:
    // RDKit❌❌:   for (auto atom : mol.atoms()) {
    // RDKit❌❌:     if (shouldRemoveH(mol, atom, ps)) {
    // RDKit❌❌:       atomsToRemove.set(atom->getIdx());
    // RDKit❌❌:     }
    // RDKit❌❌:   }  // end of the loop over atoms
    // RDKit❌❌:
    // RDKit❌❌:   // Once we know which H atoms would be removed, filter out those that
    // RDKit❌❌:   // would cause any SGroups to become empty
    // RDKit❌❌:   if (ps.removeInSGroups) {
    // RDKit❌❌:     filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌:   // now that we know which atoms need to be removed, go ahead and remove them
    // RDKit❌❌:   // NOTE: there's too much complexity around stereochemistry here
    // RDKit❌❌:   // to be able to safely use batch editing.
    // RDKit❌❌:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    // RDKit❌❌:     if (atomsToRemove[idx]) {
    // RDKit❌❌:       molRemoveH(mol, idx, ps.updateExplicitCount);
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌:   mol.clearComputedProps(true);
    // RDKit❌❌:   //
    // RDKit❌❌:   //  If we didn't only remove implicit Hs, which are guaranteed to
    // RDKit❌❌:   //  be the highest numbered atoms, we may have altered atom indices.
    // RDKit❌❌:   //  This can screw up derived properties (such as ring members), so
    // RDKit❌❌:   //  do some checks:
    // RDKit❌❌:   //
    // RDKit❌❌:   if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    // RDKit❌❌:     sanitizeMol(mol);
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌:   // if we removed Hs and any chiral atoms now have more than 1 explict H,
    // RDKit❌❌:   // remove those
    // RDKit❌❌:   if (!atomsToRemove.empty()) {
    // RDKit❌❌:     for (auto atom : mol.atoms()) {
    // RDKit❌❌:       if (!atom->getNoImplicit() &&
    // RDKit❌❌:           atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❌❌:         unsigned int numExplicitHs = atom->getNumExplicitHs();
    // RDKit❌❌:         if (numExplicitHs > 1) {
    // RDKit❌❌:           atom->setNumExplicitHs(0);
    // RDKit❌❌:           atom->updatePropertyCache(false);
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌: }
    // END RDKIT CPP FUNCTION MolOps::removeHs(RWMol&, const RemoveHsParameters&, bool)
    let _ = (molecule, params, sanitize);
    Err(RemoveHydrogensError::ProtocolDebt {
        branch: "MolOps::removeHs",
        reason: "source-alignment frame is present; remove-H assignment semantics are not yet ported",
    })
}

#[allow(dead_code)]
struct RemoveHsPredicateContext<'a> {
    molecule: &'a Molecule,
    adjacency: AdjacencyList,
}

#[allow(dead_code)]
impl<'a> RemoveHsPredicateContext<'a> {
    fn new(molecule: &'a Molecule) -> Self {
        Self {
            molecule,
            adjacency: AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds()),
        }
    }

    fn atom(&self, atom: AtomId) -> Result<&'a Atom, RemoveHydrogensError> {
        self.molecule
            .atom(atom)
            .ok_or(RemoveHydrogensError::AtomOutOfRange {
                atom,
                atom_count: self.molecule.num_atoms(),
            })
    }

    fn bond(&self, bond: BondId) -> Result<&'a Bond, RemoveHydrogensError> {
        self.molecule
            .bonds()
            .get(bond.index())
            .ok_or(RemoveHydrogensError::BondOutOfRange {
                bond,
                bond_count: self.molecule.num_bonds(),
            })
    }

    fn degree(&self, atom: AtomId) -> usize {
        self.adjacency.neighbors_of(atom.index()).len()
    }

    fn neighbors(&self, atom: AtomId) -> &[NeighborRef] {
        self.adjacency.neighbors_of(atom.index())
    }

    fn sgroup_bond_is_xbond(
        &self,
        substance_group: &SubstanceGroup,
        bond: BondId,
    ) -> Result<bool, RemoveHydrogensError> {
        let bond_data = self.bond(bond)?;
        let begin_in_group = substance_group.atoms().contains(&bond_data.begin());
        let end_in_group = substance_group.atoms().contains(&bond_data.end());
        match (begin_in_group, end_in_group) {
            (true, true) => Ok(false),
            (true, false) | (false, true) => Ok(true),
            (false, false) => Err(RemoveHydrogensError::InvalidSubstanceGroupBondRole { bond }),
        }
    }
}

fn bond_has_endpoint(bond: &Bond, atom: AtomId) -> bool {
    bond.begin() == atom || bond.end() == atom
}

fn sgroup_includes_atom(substance_group: &SubstanceGroup, atom: AtomId) -> bool {
    substance_group.atoms().contains(&atom)
        || substance_group.parent_atoms().contains(&atom)
        || substance_group.attach_points().iter().any(|attach_point| {
            attach_point.atom == atom || attach_point.leaving_atom == Some(atom)
        })
}

fn has_non_tetrahedral_stereo(atom: &Atom) -> bool {
    matches!(
        atom.chiral_tag(),
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
    )
}

fn bond_stereo_greater_than_any(stereo: BondStereo) -> bool {
    !matches!(stereo, BondStereo::None | BondStereo::Any)
}

#[allow(dead_code)]
fn should_remove_h(
    context: &RemoveHsPredicateContext<'_>,
    atom: AtomId,
    params: &RemoveHsParams,
) -> Result<bool, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION shouldRemoveH
    // RDKit✔️✔️: bool shouldRemoveH(const RWMol &mol, const Atom *atom,
    // RDKit✔️✔️:                    const RemoveHsParameters &ps) {
    // RDKit✔️✔️:   if (atom->getAtomicNum() != 1) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    let atom_data = context.atom(atom)?;
    if atom_data.atomic_number() != 1 {
        return Ok(false);
    }

    // RDKit✔️✔️:   if (!ps.removeWithQuery && atom->hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_with_query && atom_data.query().is_some() {
        return Ok(false);
    }

    let degree = context.degree(atom);

    // RDKit✔️✔️:   if (!ps.removeDegreeZero && !atom->getDegree()) {
    // RDKit✔️✔️:     if (ps.showWarnings) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "WARNING: not removing hydrogen atom without neighbors"
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_degree_zero && degree == 0 {
        if params.show_warnings {
            log::warn!("WARNING: not removing hydrogen atom without neighbors");
        }
        return Ok(false);
    }

    // RDKit✔️✔️:   if (!ps.removeHigherDegrees && atom->getDegree() > 1) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_higher_degrees && degree > 1 {
        return Ok(false);
    }

    // RDKit✔️✔️:   if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes && atom->getIsotope()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_isotopes && !params.remove_and_track_isotopes && atom_data.isotope().is_some()
    {
        return Ok(false);
    }

    // RDKit✔️✔️:   if (!ps.removeNonimplicit && !atom->hasProp(common_properties::isImplicit)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_nonimplicit && atom_data.prop("isImplicit").is_none() {
        return Ok(false);
    }

    // RDKit✔️✔️:   if (!ps.removeMapped && atom->getAtomMapNum()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_mapped && atom_data.atom_map().is_some() {
        return Ok(false);
    }

    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (ps.removeInSGroups) {
    // RDKit✔️✔️:     // If removing H in SGroups, do not remove H atoms in special
    // RDKit✔️✔️:     // roles in the SGroup
    // RDKit✔️✔️:     for (const auto &sg : getSubstanceGroups(mol)) {
    if params.remove_in_sgroups {
        for substance_group in context.molecule.substance_groups() {
            // RDKit✔️✔️:       // The H atom is one of the "caps" of the SGroup. Technically,
            // RDKit✔️✔️:       // it's not part of the group, but it defines its boundaries.
            // RDKit✔️✔️:       for (const auto &bond_idx : sg.getBonds()) {
            // RDKit✔️✔️:         if (sg.getBondType(bond_idx) == SubstanceGroup::BondType::XBOND) {
            // RDKit✔️✔️:           auto bond = mol.getBondWithIdx(bond_idx);
            // RDKit✔️✔️:           if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
            // RDKit✔️✔️:             return false;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for bond in substance_group.bonds() {
                if context.sgroup_bond_is_xbond(substance_group, *bond)? {
                    let bond_data = context.bond(*bond)?;
                    if bond_has_endpoint(bond_data, atom) {
                        return Ok(false);
                    }
                }
            }

            // RDKit✔️✔️:
            // RDKit✔️✔️:       for (const auto &sap : sg.getAttachPoints()) {
            // RDKit✔️✔️:         // The H atoms is an attach point. This would be weird, but is possible.
            // RDKit✔️✔️:         // (if it is a 'leaving atom' we don't care, though)
            // RDKit✔️✔️:         if (sap.aIdx == atom->getIdx()) {
            // RDKit✔️✔️:           return false;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for attach_point in substance_group.attach_points() {
                if attach_point.atom == atom {
                    return Ok(false);
                }
            }

            // RDKit✔️✔️:
            // RDKit✔️✔️:       for (const auto &cs : sg.getCStates()) {
            // RDKit✔️✔️:         // The bond to the H atom defines a CState
            // RDKit✔️✔️:         auto bond = mol.getBondWithIdx(cs.bondIdx);
            // RDKit✔️✔️:         if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
            // RDKit✔️✔️:           return false;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for cstate in substance_group.cstates() {
                let bond_data = context.bond(cstate.bond)?;
                if bond_has_endpoint(bond_data, atom) {
                    return Ok(false);
                }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   } else {
    } else {
        // RDKit✔️✔️:     for (const auto &sg : getSubstanceGroups(mol)) {
        // RDKit✔️✔️:       if (sg.includesAtom(atom->getIdx())) {
        // RDKit✔️✔️:         return false;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for substance_group in context.molecule.substance_groups() {
            if sgroup_includes_atom(substance_group, atom) {
                return Ok(false);
            }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   if (!ps.removeHydrides && atom->getFormalCharge() == -1) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_hydrides && atom_data.formal_charge() == -1 {
        return Ok(false);
    }

    // RDKit✔️✔️:   bool removeIt = true;
    let remove_it = true;

    // RDKit✔️✔️:   if (atom->getDegree() &&
    // RDKit✔️✔️:       (!ps.removeDummyNeighbors || !ps.removeDefiningBondStereo ||
    // RDKit✔️✔️:        !ps.removeOnlyHNeighbors || !ps.removeNontetrahedralNeighbors ||
    // RDKit✔️✔️:        !ps.removeWithWedgedBond)) {
    if degree != 0
        && (!params.remove_dummy_neighbors
            || !params.remove_defining_bond_stereo
            || !params.remove_only_h_neighbors
            || !params.remove_nontetrahedral_neighbors
            || !params.remove_with_wedged_bond)
    {
        // RDKit✔️✔️:     bool onlyHNeighbors = true;
        let mut only_h_neighbors = true;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        for neighbor_ref in context.neighbors(atom) {
            let neighbor = AtomId::new(neighbor_ref.atom_index);
            let neighbor_data = context.atom(neighbor)?;
            // RDKit✔️✔️:       // is it a dummy?
            // RDKit✔️✔️:       if (!ps.removeDummyNeighbors && nbr->getAtomicNum() < 1) {
            // RDKit✔️✔️:         if (ps.showWarnings) {
            // RDKit✔️✔️:           BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
            // RDKit✔️✔️:                                      "with dummy atom neighbors"
            // RDKit✔️✔️:                                   << std::endl;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            if !params.remove_dummy_neighbors && neighbor_data.atomic_number() < 1 {
                if params.show_warnings {
                    log::warn!("WARNING: not removing hydrogen atom with dummy atom neighbors");
                }
                return Ok(false);
            }

            // RDKit✔️✔️:       // does it have non-tetrahedral stereo:
            // RDKit✔️✔️:       if (!ps.removeNontetrahedralNeighbors &&
            // RDKit✔️✔️:           Chirality::hasNonTetrahedralStereo(nbr)) {
            // RDKit✔️✔️:         if (ps.showWarnings) {
            // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
            // RDKit✔️✔️:               << "WARNING: not removing hydrogen atom "
            // RDKit✔️✔️:                  "with neighbor that has non-tetrahedral stereochemistry"
            // RDKit✔️✔️:               << std::endl;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            if !params.remove_nontetrahedral_neighbors && has_non_tetrahedral_stereo(neighbor_data)
            {
                if params.show_warnings {
                    log::warn!(
                        "WARNING: not removing hydrogen atom with neighbor that has non-tetrahedral stereochemistry"
                    );
                }
                return Ok(false);
            }

            // RDKit✔️✔️:       if (!ps.removeOnlyHNeighbors && nbr->getAtomicNum() != 1) {
            // RDKit✔️✔️:         onlyHNeighbors = false;
            // RDKit✔️✔️:       }
            if !params.remove_only_h_neighbors && neighbor_data.atomic_number() != 1 {
                only_h_neighbors = false;
            }

            // RDKit✔️✔️:       if (!ps.removeWithWedgedBond) {
            // RDKit✔️✔️:         const auto bnd = mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
            // RDKit✔️✔️:         if (bnd->getBondDir() == Bond::BEGINDASH ||
            // RDKit✔️✔️:             bnd->getBondDir() == Bond::BEGINWEDGE) {
            // RDKit✔️✔️:           if (ps.showWarnings) {
            // RDKit✔️✔️:             BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
            // RDKit✔️✔️:                                        "with wedged bond"
            // RDKit✔️✔️:                                     << std::endl;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           return false;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            if !params.remove_with_wedged_bond {
                let hydrogen_bond = context.bond(neighbor_ref.bond)?;
                if matches!(
                    hydrogen_bond.direction(),
                    BondDirection::BeginDash | BondDirection::BeginWedge
                ) {
                    if params.show_warnings {
                        log::warn!("WARNING: not removing hydrogen atom with wedged bond");
                    }
                    return Ok(false);
                }
            }

            // RDKit✔️✔️:       // Check to see if the neighbor has a double bond and we're the only
            // RDKit✔️✔️:       // neighbor at this end.  This was part of github #1810
            // RDKit✔️✔️:       if (!ps.removeDefiningBondStereo && nbr->getDegree() == 2) {
            // RDKit✔️✔️:         for (const auto bnd : mol.atomBonds(nbr)) {
            // RDKit✔️✔️:           if (bnd->getBondType() == Bond::DOUBLE &&
            // RDKit✔️✔️:               (bnd->getStereo() > Bond::STEREOANY ||
            // RDKit✔️✔️:                mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
            // RDKit✔️✔️:                        ->getBondDir() > Bond::NONE)) {
            // RDKit✔️✔️:             return false;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            if !params.remove_defining_bond_stereo && context.degree(neighbor) == 2 {
                for neighbor_bond_ref in context.neighbors(neighbor) {
                    let neighbor_bond = context.bond(neighbor_bond_ref.bond)?;
                    let hydrogen_bond = context.bond(neighbor_ref.bond)?;
                    if neighbor_bond.order() == BondOrder::Double
                        && (bond_stereo_greater_than_any(neighbor_bond.stereo())
                            || hydrogen_bond.direction() != BondDirection::None)
                    {
                        return Ok(false);
                    }
                }
            }
            // RDKit✔️✔️:     }
        }

        // RDKit✔️✔️:     if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if remove_it && (!params.remove_only_h_neighbors && only_h_neighbors) {
            return Ok(false);
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   return removeIt;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION shouldRemoveH
    Ok(remove_it)
}

#[allow(dead_code)]
fn mol_remove_h_assignment(
    molecule: &Molecule,
    atom: AtomId,
    params: &RemoveHsParams,
    assignment: &mut RemoveHsAssignment,
) -> Result<(), RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION molRemoveH
    // RDKit❌❌: void molRemoveH(RWMol &mol, unsigned int idx, bool updateExplicitCount) {
    // RDKit❌❌:   auto atom = mol.getAtomWithIdx(idx);
    // RDKit❌❌:   PRECONDITION(atom->getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");
    // RDKit❌❌:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit❌❌:     Atom *heavyAtom = bond->getOtherAtom(atom);
    // RDKit❌❌:     int heavyAtomNum = heavyAtom->getAtomicNum();
    // RDKit❌❌:
    // RDKit❌❌:     // we'll update the neighbor's explicit H count if we were told to
    // RDKit❌❌:     // *or* if the neighbor is chiral, in which case the H is needed
    // RDKit❌❌:     // in order to complete the coordination
    // RDKit❌❌:     // *or* if the neighbor has the noImplicit flag set:
    // RDKit❌❌:     if (updateExplicitCount || heavyAtom->getNoImplicit() ||
    // RDKit❌❌:         heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❌❌:       heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
    // RDKit❌❌:     } else {
    // RDKit❌❌:       // this is a special case related to Issue 228 and the
    // RDKit❌❌:       // "disappearing Hydrogen" problem discussed in MolOps::adjustHs
    // RDKit❌❌:       //
    // RDKit❌❌:       // If we remove a hydrogen from an aromatic N or P, or if
    // RDKit❌❌:       // the heavy atom it is connected to is not in its default
    // RDKit❌❌:       // valence state, we need to be *sure* to increment the
    // RDKit❌❌:       // explicit count, even if the H itself isn't marked as explicit
    // RDKit❌❌:       const INT_VECT &defaultVs =
    // RDKit❌❌:           PeriodicTable::getTable()->getValenceList(heavyAtomNum);
    // RDKit❌❌:       if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
    // RDKit❌❌:             may_need_extra_H(mol, heavyAtom)) &&
    // RDKit❌❌:            isAromaticAtom(*heavyAtom)) ||
    // RDKit❌❌:           (std::find(defaultVs.begin() + 1, defaultVs.end(),
    // RDKit❌❌:                      heavyAtom->getTotalValence()) != defaultVs.end())) {
    // RDKit❌❌:         heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:
    // RDKit❌❌:     // One other consequence of removing the H from the graph is
    // RDKit❌❌:     // that we may change the ordering of the bonds about a
    // RDKit❌❌:     // chiral center.  This may change the chiral label at that
    // RDKit❌❌:     // atom.  We deal with that by explicitly checking here:
    // RDKit❌❌:     if (heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❌❌:       INT_LIST neighborIndices;
    // RDKit❌❌:       for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit❌❌:         if (nbnd->getIdx() != bond->getIdx()) {
    // RDKit❌❌:           neighborIndices.push_back(nbnd->getIdx());
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:       neighborIndices.push_back(bond->getIdx());
    // RDKit❌❌:
    // RDKit❌❌:       int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
    // RDKit❌❌:       // std::cerr << "H: "<<atom->getIdx()<<" hvy:
    // RDKit❌❌:       // "<<heavyAtom->getIdx()<<" swaps: " << nSwaps<<std::endl;
    // RDKit❌❌:       if (nSwaps % 2) {
    // RDKit❌❌:         heavyAtom->invertChirality();
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:
    // RDKit❌❌:     // If we are removing a H atom that defines bond stereo (e.g. imines),
    // RDKit❌❌:     // Then also remove the bond stereo information, as it is no longer valid.
    // RDKit❌❌:     if (heavyAtom->getDegree() == 2) {
    // RDKit❌❌:       for (auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit❌❌:         if (nbnd != bond) {
    // RDKit❌❌:           if (nbnd->getStereo() > Bond::STEREOANY) {
    // RDKit❌❌:             nbnd->setStereo(Bond::STEREONONE);
    // RDKit❌❌:             nbnd->getStereoAtoms().clear();
    // RDKit❌❌:           }
    // RDKit❌❌:           break;
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:
    // RDKit❌❌:     // if it's a wavy bond, then we need to
    // RDKit❌❌:     // mark the beginning atom with the _UnknownStereo tag.
    // RDKit❌❌:     // so that we know later that something was affecting its
    // RDKit❌❌:     // stereochem
    // RDKit❌❌:     if (bond->getBondDir() == Bond::UNKNOWN &&
    // RDKit❌❌:         bond->getBeginAtomIdx() == heavyAtom->getIdx()) {
    // RDKit❌❌:       heavyAtom->setProp(common_properties::_UnknownStereo, 1);
    // RDKit❌❌:     } else if (bond->getBondDir() == Bond::ENDDOWNRIGHT ||
    // RDKit❌❌:                bond->getBondDir() == Bond::ENDUPRIGHT) {
    // RDKit❌❌:       // if the direction is set on this bond and the atom it's connected to
    // RDKit❌❌:       // has no other single bonds with directions set, then we need to set
    // RDKit❌❌:       // direction on one of the other neighbors in order to avoid double
    // RDKit❌❌:       // bond stereochemistry possibly being lost. This was github #754
    // RDKit❌❌:       bool foundADir = false;
    // RDKit❌❌:       Bond *oBond = nullptr;
    // RDKit❌❌:       for (const auto &nbri :
    // RDKit❌❌:            boost::make_iterator_range(mol.getAtomBonds(heavyAtom))) {
    // RDKit❌❌:         Bond *nbnd = mol[nbri];
    // RDKit❌❌:         if (nbnd->getIdx() != bond->getIdx() &&
    // RDKit❌❌:             nbnd->getBondType() == Bond::SINGLE) {
    // RDKit❌❌:           if (nbnd->getBondDir() == Bond::NONE) {
    // RDKit❌❌:             oBond = nbnd;
    // RDKit❌❌:           } else {
    // RDKit❌❌:             foundADir = true;
    // RDKit❌❌:           }
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:       if (!foundADir && oBond != nullptr) {
    // RDKit❌❌:         bool flipIt = (oBond->getBeginAtom() == heavyAtom) &&
    // RDKit❌❌:                       (bond->getBeginAtom() == heavyAtom);
    // RDKit❌❌:         if (flipIt) {
    // RDKit❌❌:           oBond->setBondDir(bond->getBondDir() == Bond::ENDDOWNRIGHT
    // RDKit❌❌:                                 ? Bond::ENDUPRIGHT
    // RDKit❌❌:                                 : Bond::ENDDOWNRIGHT);
    // RDKit❌❌:         } else {
    // RDKit❌❌:           oBond->setBondDir(bond->getBondDir());
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:       // if this atom is one of the stereoatoms for a double bond we need
    // RDKit❌❌:       // to switch the stereo atom on this end to be the other neighbor
    // RDKit❌❌:       // This was part of github #1810
    // RDKit❌❌:       adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
    // RDKit❌❌:     } else {
    // RDKit❌❌:       // if this atom is one of the stereoatoms for a double bond we need
    // RDKit❌❌:       // to switch the stereo atom on this end to be the other neighbor
    // RDKit❌❌:       // This was part of github #1810
    // RDKit❌❌:       adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
    // RDKit❌❌:     }
    // RDKit❌❌:
    // RDKit❌❌:     // remove the bond from any SGroups that might include it.
    // RDKit❌❌:     for (auto &sg : getSubstanceGroups(mol)) {
    // RDKit❌❌:       sg.removeBondWithIdx(bond->getIdx());
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌:   // Finally, remove the atom from any SGroups that might include it, so that
    // RDKit❌❌:   // the SGroups don't get removed in removeAtom(). Since we allow removing
    // RDKit❌❌:   // SGroup SAP lvidx H atoms, we need to check for those and update them.
    // RDKit❌❌:   for (auto &sg : getSubstanceGroups(mol)) {
    // RDKit❌❌:     sg.removeAtomWithIdx(idx);
    // RDKit❌❌:     sg.removeParentAtomWithIdx(idx);
    // RDKit❌❌:
    // RDKit❌❌:     for (auto &sap : sg.getAttachPoints()) {
    // RDKit❌❌:       if (sap.lvIdx == static_cast<int>(idx)) {
    // RDKit❌❌:         sap.lvIdx = -1;
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❌❌:   // computed properties will be cleared after all hydrogens are removed
    // RDKit❌❌:   bool clearProps = false;
    // RDKit❌❌:   mol.removeAtom(atom, clearProps);
    // RDKit❌❌: }
    // END RDKIT CPP FUNCTION molRemoveH
    let _ = (molecule, atom, params, assignment);
    Err(RemoveHydrogensError::ProtocolDebt {
        branch: "molRemoveH",
        reason: "source-alignment frame is present; per-hydrogen removal assignment semantics are not yet ported",
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomQueryPredicate, AtomSpec, BondSpec, Element, MoleculeBuilder, QueryNode,
        SGroupAttachPoint, SGroupCState, SubstanceGroupId, SubstanceGroupKind,
    };

    #[test]
    fn remove_hs_assignment_exposes_protocol_debt() {
        let molecule = Molecule::default();
        let error = remove_hs_assignment(&molecule, &RemoveHsParams::default(), true).unwrap_err();

        assert!(matches!(
            error,
            RemoveHydrogensError::ProtocolDebt {
                branch: "MolOps::removeHs",
                ..
            }
        ));
    }

    #[test]
    fn should_remove_h_matches_basic_rdkit_guards() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        assert!(!should(&molecule, carbon, &RemoveHsParams::default()));
        assert!(should(&molecule, hydrogen, &RemoveHsParams::default()));

        let mut builder = MoleculeBuilder::new();
        let isolated_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let molecule = builder.build().unwrap();
        assert!(!should(
            &molecule,
            isolated_hydrogen,
            &RemoveHsParams::default()
        ));
        let mut params = RemoveHsParams::default();
        params.remove_degree_zero = true;
        assert!(should(&molecule, isolated_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let query_hydrogen = builder.add_atom(
            AtomSpec::new(Element::H).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let molecule = builder.build().unwrap();
        let mut params = RemoveHsParams {
            remove_degree_zero: true,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, query_hydrogen, &params));
        params.remove_with_query = true;
        assert!(should(&molecule, query_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let isotopic_hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        let molecule = builder.build().unwrap();
        let mut params = RemoveHsParams {
            remove_degree_zero: true,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, isotopic_hydrogen, &params));
        params.remove_isotopes = true;
        assert!(should(&molecule, isotopic_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let mapped_hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_atom_map(7));
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_degree_zero: true,
            remove_mapped: false,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, mapped_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let hydride = builder.add_atom(AtomSpec::new(Element::H).with_formal_charge(-1));
        let molecule = builder.build().unwrap();
        let mut params = RemoveHsParams {
            remove_degree_zero: true,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, hydride, &params));
        params.remove_hydrides = true;
        assert!(should(&molecule, hydride, &params));
    }

    #[test]
    fn should_remove_h_matches_neighbor_guard_branches() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let dummy = builder.add_atom(AtomSpec::new(Element::DUMMY));
        builder
            .add_bond(BondSpec::new(hydrogen, dummy, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));
        let mut params = RemoveHsParams::default();
        params.remove_dummy_neighbors = true;
        assert!(should(&molecule, hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let left_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let right_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(
                left_hydrogen,
                right_hydrogen,
                BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(
            &molecule,
            left_hydrogen,
            &RemoveHsParams::default()
        ));
        let mut params = RemoveHsParams::default();
        params.remove_only_h_neighbors = true;
        assert!(should(&molecule, left_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(
                BondSpec::new(carbon, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_with_wedged_bond: false,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let center =
            builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::SquarePlanar));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));
        let mut params = RemoveHsParams::default();
        params.remove_nontetrahedral_neighbors = true;
        assert!(should(&molecule, hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let left = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(left, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(hydrogen, right, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));
        let mut params = RemoveHsParams::default();
        params.remove_higher_degrees = true;
        assert!(should(&molecule, hydrogen, &params));
    }

    #[test]
    fn should_remove_h_keeps_hydrogen_defining_double_bond_stereo() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(hydrogen, nitrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(nitrogen, carbon, BondOrder::Double)
                    .with_stereo_atoms(hydrogen, carbon)
                    .with_stereo(BondStereo::Cis),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));
        let mut params = RemoveHsParams::default();
        params.remove_defining_bond_stereo = true;
        assert!(should(&molecule, hydrogen, &params));
    }

    #[test]
    fn should_remove_h_matches_sgroup_membership_and_role_guards() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![hydrogen]),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_degree_zero: true,
            remove_in_sgroups: false,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let crossing_bond = builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![carbon])
                    .with_bonds(vec![crossing_bond]),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));

        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![carbon])
                    .with_attach_points(vec![SGroupAttachPoint {
                        atom: hydrogen,
                        leaving_atom: None,
                        label: None,
                        order: None,
                    }]),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));

        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let cstate_bond = builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![carbon])
                    .with_cstates(vec![SGroupCState {
                        bond: cstate_bond,
                        vector: [0.0, 0.0],
                    }]),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, hydrogen, &RemoveHsParams::default()));
    }

    fn should(molecule: &Molecule, atom: AtomId, params: &RemoveHsParams) -> bool {
        let context = RemoveHsPredicateContext::new(molecule);
        should_remove_h(&context, atom, params).unwrap()
    }
}
