// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::collections::{BTreeMap, BTreeSet};

use crate::{
    AdjacencyList, Atom, AtomId, AtomPdbResidueInfo, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag,
    Molecule, NeighborRef, SubstanceGroup, ValenceModel, ops::MoleculeReadParts,
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AddHydrogensError {
    #[error("unsupported atom for hydrogen addition at {atom} with atomic number {atomic_number}")]
    UnsupportedAtom { atom: crate::AtomId, atomic_number: u8 },
    #[error("incomplete RDKit addHs port for {branch}: {reason}")]
    ProtocolDebt { branch: &'static str, reason: &'static str },
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
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
    ProtocolDebt { branch: &'static str, reason: &'static str },
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
pub struct RemoveHsParams {
    pub remove_degree_zero: bool,
    pub remove_higher_degrees: bool,
    pub remove_only_h_neighbors: bool,
    pub remove_isotopes: bool,
    pub remove_and_track_isotopes: bool,
    pub remove_dummy_neighbors: bool,
    pub remove_defining_bond_stereo: bool,
    pub remove_with_wedged_bond: bool,
    pub remove_with_query: bool,
    pub remove_mapped: bool,
    pub remove_in_sgroups: bool,
    pub show_warnings: bool,
    pub remove_nonimplicit: bool,
    pub update_explicit_count: bool,
    pub remove_hydrides: bool,
    pub remove_nontetrahedral_neighbors: bool,
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

// BEGIN RDKIT CPP STRUCT MolOps::AddHsParameters
// RDKit✔️✔️: struct RDKIT_GRAPHMOL_EXPORT AddHsParameters {
// RDKit✔️✔️:   bool explicitOnly = false;   /**< only add explicit Hs */
// RDKit✔️✔️:   bool addCoords = false;      /**< add coordinates for the Hs */
// RDKit✔️✔️:   bool addResidueInfo = false; /**< add residue info to the Hs */
// RDKit✔️✔️:   bool skipQueries =
// RDKit✔️✔️:       false; /**< do not add Hs to query atoms or atoms with query bonds */
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT MolOps::AddHsParameters
#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub struct AddHsParams {
    pub explicit_only: bool,
    pub add_coords: bool,
    pub add_residue_info: bool,
    pub skip_queries: bool,
    pub only_on_atoms: Option<Vec<AtomId>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub(crate) struct AddHsAssignment {
    pub(crate) hydrogens_to_add: Vec<AddHydrogen>,
    pub(crate) atoms_to_update_property_cache: Vec<AtomId>,
    pub(crate) valence_before_add_hs: Option<crate::ValenceAssignment>,
    pub(crate) atom_explicit_hydrogen_updates: Vec<AtomExplicitHydrogenUpdate>,
    pub(crate) atom_pdb_residue_info_updates: Vec<AtomPdbResidueInfoUpdate>,
    pub(crate) clear_isotopic_hydrogen_properties: Vec<AtomId>,
    pub(crate) clear_computed_properties: bool,
    pub(crate) add_terminal_coordinates: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AddHydrogen {
    pub(crate) heavy_atom: AtomId,
    pub(crate) isotope: Option<u16>,
    pub(crate) is_implicit: bool,
    pub(crate) props: BTreeMap<String, String>,
    pub(crate) pdb_residue_info: Option<AtomPdbResidueInfo>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub(crate) struct RemoveHsAssignment {
    pub(crate) atoms_to_remove: Vec<AtomId>,
    pub(crate) atom_valence_cache_updates: Vec<AtomValenceCacheUpdate>,
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
pub(crate) struct AtomValenceCacheUpdate {
    pub(crate) atom: AtomId,
    pub(crate) explicit_valence: i32,
    pub(crate) implicit_hydrogens: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AtomExplicitHydrogenUpdate {
    pub(crate) atom: AtomId,
    pub(crate) explicit_hydrogens: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AtomPdbResidueInfoUpdate {
    pub(crate) atom: AtomId,
    pub(crate) pdb_residue_info: AtomPdbResidueInfo,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomPropertyUpdate {
    SetUnknownStereo { atom: AtomId },
    SetIsotopicHydrogens { atom: AtomId, isotopes: Vec<u16> },
    ClearIsotopicHydrogens { atom: AtomId },
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

pub(crate) fn add_hs_assignment(
    read_parts: MoleculeReadParts<'_>,
    params: &AddHsParams,
) -> Result<AddHsAssignment, AddHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::addHs(RWMol&, const AddHsParameters&, const UINT_VECT*)
    // RDKit✔️❌: void addHs(RWMol &mol, const AddHsParameters &params,
    // RDKit✔️❌:            const UINT_VECT *onlyOnAtoms) {
    // RDKit✔️❌:   mol.clearComputedProps(false);
    // RDKit✔️❌:   unsigned int numAddHyds = 0;
    // RDKit✔️❌:   boost::dynamic_bitset<> onAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   if (onlyOnAtoms) {
    // RDKit✔️✔️:     for (auto atIdx : *onlyOnAtoms) {
    // RDKit✔️✔️:       onAtoms.set(atIdx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     onAtoms.set();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned int> numExplicitHs(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   std::vector<unsigned int> numImplicitHs(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   for (auto at : mol.atoms()) {
    // RDKit✔️✔️:     numExplicitHs[at->getIdx()] = at->getNumExplicitHs();
    // RDKit✔️✔️:     numImplicitHs[at->getIdx()] = at->getNumImplicitHs();
    // RDKit✔️✔️:     if (onAtoms[at->getIdx()]) {
    // RDKit✔️✔️:       if (params.skipQueries && isQueryAtom(mol, *at)) {
    // RDKit✔️✔️:         onAtoms.set(at->getIdx(), 0);
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       numAddHyds += at->getNumExplicitHs();
    // RDKit✔️✔️:       if (!params.explicitOnly) {
    // RDKit✔️✔️:         numAddHyds += at->getNumImplicitHs();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️❌:   unsigned int nSize = mol.getNumAtoms() + numAddHyds;
    // RDKit✔️❌:   for (auto cfi = mol.beginConformers(); cfi != mol.endConformers(); ++cfi) {
    // RDKit✔️❌:     (*cfi)->reserve(nSize);
    // RDKit✔️❌:   }
    // RDKit✔️✔️:   unsigned int stopIdx = mol.getNumAtoms();
    // RDKit✔️✔️:   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    // RDKit✔️✔️:     if (!onAtoms[aidx]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Atom *newAt = mol.getAtomWithIdx(aidx);
    // RDKit✔️❌:     std::vector<unsigned int> isoHs;
    // RDKit✔️❌:     if (newAt->getPropIfPresent(common_properties::_isotopicHs, isoHs)) {
    // RDKit✔️❌:       newAt->clearProp(common_properties::_isotopicHs);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<unsigned int>::const_iterator isoH = isoHs.begin();
    // RDKit✔️✔️:     unsigned int newIdx;
    // RDKit✔️❌:     newAt->clearComputedProps();
    // RDKit✔️✔️:     unsigned int onumexpl = numExplicitHs[aidx];
    // RDKit✔️✔️:     for (unsigned int i = 0; i < onumexpl; i++) {
    // RDKit✔️✔️:       newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️✔️:       mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RDKit✔️✔️:       auto hAtom = mol.getAtomWithIdx(newIdx);
    // RDKit✔️❌:       hAtom->updatePropertyCache();
    // RDKit✔️❌:       if (params.addCoords) {
    // RDKit✔️❌:         setTerminalAtomCoords(mol, newIdx, aidx);
    // RDKit✔️❌:       }
    // RDKit✔️✔️:       if (isoH != isoHs.end()) {
    // RDKit✔️✔️:         hAtom->setIsotope(*isoH);
    // RDKit✔️✔️:         ++isoH;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     newAt->setNumExplicitHs(0);
    // RDKit✔️✔️:     if (!params.explicitOnly) {
    // RDKit✔️✔️:       for (unsigned int i = 0; i < numImplicitHs[aidx]; i++) {
    // RDKit✔️✔️:         newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️✔️:         mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RDKit✔️✔️:         auto hAtom = mol.getAtomWithIdx(newIdx);
    // RDKit✔️✔️:         hAtom->setProp(common_properties::isImplicit, 1);
    // RDKit✔️❌:         hAtom->updatePropertyCache();
    // RDKit✔️❌:         if (params.addCoords) {
    // RDKit✔️❌:           setTerminalAtomCoords(mol, newIdx, aidx);
    // RDKit✔️❌:         }
    // RDKit✔️✔️:         if (isoH != isoHs.end()) {
    // RDKit✔️✔️:           hAtom->setIsotope(*isoH);
    // RDKit✔️✔️:           ++isoH;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️❌:     newAt->updatePropertyCache(false);
    // RDKit❗❌:     if (isoH != isoHs.end()) { ... warn ... }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.addResidueInfo) {
    // RDKit✔️✔️:     AssignHsResidueInfo(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::addHs(RWMol&, const AddHsParameters&, const UINT_VECT*)
    // BEGIN RDKIT CPP BODY MolOps::addHs assignment planner
    if let Some(only_on_atoms) = &params.only_on_atoms {
        for atom in only_on_atoms {
            if read_parts.atom(*atom).is_none() {
                return Err(AddHydrogensError::UnsupportedAtom {
                    atom: *atom,
                    atomic_number: 0,
                });
            }
        }
    }
    let atom_count = read_parts.num_atoms();
    let valence = match read_parts.derived_cache().valence.as_ref() {
        Some(valence)
            if valence.explicit_valence.len() == atom_count && valence.implicit_hydrogens.len() == atom_count =>
        {
            valence.clone()
        }
        Some(_) => {
            return Err(crate::ValenceError::UnsupportedBranch {
                reason: "AddHs input valence cache length does not match the atom count",
            }
            .into());
        }
        None => {
            // COSMolKit builders store valence in a molecule-level typed cache,
            // while RDKit stores it on each Atom. Materialize the same
            // Atom::updatePropertyCache(false) values when that typed cache has
            // not yet been initialized, then follow AddHs.cpp incrementally.
            read_parts.assign_valence_with_options(ValenceModel::RdkitLike, false)?
        }
    };
    let context = RemoveHsPredicateContext::new(read_parts);
    let mut on_atoms = vec![params.only_on_atoms.is_none(); atom_count];
    if let Some(only_on_atoms) = &params.only_on_atoms {
        for atom in only_on_atoms {
            on_atoms[atom.index()] = true;
        }
    }
    let mut assignment = AddHsAssignment {
        valence_before_add_hs: Some(valence.clone()),
        clear_computed_properties: true,
        add_terminal_coordinates: params.add_coords,
        ..AddHsAssignment::default()
    };
    for atom in read_parts.atoms() {
        if !on_atoms[atom.id().index()] {
            continue;
        }
        if params.skip_queries && is_add_hs_query_atom(&context, atom.id())? {
            on_atoms[atom.id().index()] = false;
            continue;
        }
        assignment.atoms_to_update_property_cache.push(atom.id());
        let isotopic_hydrogens = isotopic_hydrogen_property(atom);
        if !atom.tracked_isotopic_hydrogens().is_empty() {
            assignment.clear_isotopic_hydrogen_properties.push(atom.id());
        }
        let mut isotopes = isotopic_hydrogens.iter().copied();

        for _ in 0..atom.explicit_hydrogens() {
            assignment.hydrogens_to_add.push(AddHydrogen {
                heavy_atom: atom.id(),
                isotope: isotopes.next(),
                is_implicit: false,
                props: BTreeMap::new(),
                pdb_residue_info: None,
            });
        }
        if !params.explicit_only {
            let implicit_count = valence.implicit_hydrogens[atom.id().index()].max(0);
            for _ in 0..implicit_count {
                assignment.hydrogens_to_add.push(AddHydrogen {
                    heavy_atom: atom.id(),
                    isotope: isotopes.next(),
                    is_implicit: true,
                    props: BTreeMap::new(),
                    pdb_residue_info: None,
                });
            }
        }
        if atom.explicit_hydrogens() != 0 {
            assignment
                .atom_explicit_hydrogen_updates
                .push(AtomExplicitHydrogenUpdate {
                    atom: atom.id(),
                    explicit_hydrogens: 0,
                });
        }
    }
    if params.add_residue_info {
        add_hs_residue_info_assignment(read_parts, &mut assignment);
    }
    // END RDKIT CPP BODY MolOps::addHs assignment planner
    Ok(assignment)
}

pub(crate) fn add_hs_valence_assignment_from_parts(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    old_atom_count: usize,
    assignment: &AddHsAssignment,
) -> Result<crate::ValenceAssignment, AddHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::addHs(RWMol&, const AddHsParameters&, const UINT_VECT*)
    // RDKit✔️❌:     newAt->clearComputedProps();
    // RDKit✔️❌:     unsigned int onumexpl = numExplicitHs[aidx];
    // RDKit✔️❌:     for (unsigned int i = 0; i < onumexpl; i++) {
    // RDKit✔️❌:       newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️❌:       mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RDKit✔️❌:       auto hAtom = mol.getAtomWithIdx(newIdx);
    // RDKit✔️❌:       hAtom->updatePropertyCache();
    // RDKit✔️❌:     }
    // RDKit✔️❌:     newAt->setNumExplicitHs(0);
    // RDKit✔️❌:     if (!params.explicitOnly) {
    // RDKit✔️❌:       for (unsigned int i = 0; i < numImplicitHs[aidx]; i++) {
    // RDKit✔️❌:         newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️❌:         mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RDKit✔️❌:         auto hAtom = mol.getAtomWithIdx(newIdx);
    // RDKit✔️❌:         hAtom->setProp(common_properties::isImplicit, 1);
    // RDKit✔️❌:         hAtom->updatePropertyCache();
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     newAt->updatePropertyCache(false);
    // END RDKIT CPP FUNCTION MolOps::addHs(RWMol&, const AddHsParameters&, const UINT_VECT*)
    // The molecule-level vectors require one O(numAtoms) clone that RDKit's
    // per-Atom cache layout does not. The per-atom calculations and source loop
    // ordering below are otherwise the same.
    let valence_before_add_hs =
        assignment
            .valence_before_add_hs
            .as_ref()
            .ok_or(crate::ValenceError::UnsupportedBranch {
                reason: "AddHs assignment is missing its original valence cache",
            })?;
    if valence_before_add_hs.explicit_valence.len() != old_atom_count
        || valence_before_add_hs.implicit_hydrogens.len() != old_atom_count
    {
        return Err(crate::ValenceError::UnsupportedBranch {
            reason: "AddHs assignment valence cache length does not match the original atom count",
        }
        .into());
    }
    if atoms.len() != old_atom_count + assignment.hydrogens_to_add.len() {
        return Err(crate::ValenceError::UnsupportedBranch {
            reason: "AddHs assignment hydrogen count does not match the appended topology",
        }
        .into());
    }

    let mut valence = valence_before_add_hs.clone();
    valence.explicit_valence.resize(atoms.len(), -1);
    valence.implicit_hydrogens.resize(atoms.len(), -1);

    let mut hydrogen_offset = 0usize;
    for heavy_atom in &assignment.atoms_to_update_property_cache {
        while assignment
            .hydrogens_to_add
            .get(hydrogen_offset)
            .is_some_and(|hydrogen| hydrogen.heavy_atom == *heavy_atom)
        {
            let hydrogen = AtomId::new(old_atom_count + hydrogen_offset);
            let (explicit, implicit) =
                crate::valence::assign_valence_state_for_atom_from_parts(atoms, bonds, adjacency, hydrogen, true)?;
            valence.explicit_valence[hydrogen.index()] = explicit;
            valence.implicit_hydrogens[hydrogen.index()] = implicit;
            hydrogen_offset += 1;
        }

        let (explicit, implicit) =
            crate::valence::assign_valence_state_for_atom_from_parts(atoms, bonds, adjacency, *heavy_atom, false)?;
        valence.explicit_valence[heavy_atom.index()] = explicit;
        valence.implicit_hydrogens[heavy_atom.index()] = implicit;
    }
    if hydrogen_offset != assignment.hydrogens_to_add.len() {
        return Err(crate::ValenceError::UnsupportedBranch {
            reason: "AddHs assignment hydrogens are not grouped by processed heavy atom",
        }
        .into());
    }

    Ok(valence)
}

#[derive(Debug, Clone, Default)]
struct AddHsResidueInfoAssignment {
    max_serial: i32,
    current_residue: Option<(String, String)>,
    current_h_id: usize,
}

impl AddHsResidueInfoAssignment {
    fn info_for_new_hydrogen(&mut self, info: &AtomPdbResidueInfo) -> AtomPdbResidueInfo {
        let residue = (info.residue_number().to_string(), info.chain_id().to_string());
        self.current_h_id += 1;
        if self.current_residue.as_ref() != Some(&residue) {
            self.current_h_id = 1;
            self.current_residue = Some(residue);
        }

        let label = hydrogen_residue_label(self.current_h_id);
        let new_info = AtomPdbResidueInfo::new(
            label,
            self.max_serial,
            info.residue_name(),
            info.residue_number(),
            info.chain_id(),
            info.is_hetero_atom(),
        );
        self.max_serial = self.max_serial.saturating_add(1);
        new_info
    }
}

fn add_hs_residue_info_assignment(read_parts: MoleculeReadParts<'_>, assignment: &mut AddHsAssignment) {
    // BEGIN RDKIT CPP FUNCTION AssignHsResidueInfo
    // RDKit✔️✔️: void AssignHsResidueInfo(RWMol &mol) {
    // RDKit✔️✔️:   int max_serial = 0;
    // RDKit✔️✔️:   unsigned int stopIdx = mol.getNumAtoms();
    // RDKit✔️✔️:   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) { ... max_serial ... }
    let mut residue_assignment = AddHsResidueInfoAssignment::default();
    for atom in read_parts.atoms() {
        if let Some(info) = atom.pdb_residue_info() {
            residue_assignment.max_serial = residue_assignment.max_serial.max(info.serial_number());
        }
    }

    let context = RemoveHsPredicateContext::new(read_parts);
    let mut new_hydrogens_by_heavy = BTreeMap::<AtomId, Vec<usize>>::new();
    for (idx, hydrogen) in assignment.hydrogens_to_add.iter().enumerate() {
        new_hydrogens_by_heavy.entry(hydrogen.heavy_atom).or_default().push(idx);
    }
    let mut planned_existing_hydrogens = BTreeSet::<AtomId>::new();

    // RDKit✔️✔️:   AtomPDBResidueInfo *current_info = nullptr;
    // RDKit✔️✔️:   int current_h_id = 0;
    // RDKit✔️✔️:   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    // RDKit✔️✔️:     Atom *newAt = mol.getAtomWithIdx(aidx);
    // RDKit✔️✔️:     auto *info = (AtomPDBResidueInfo *)(newAt->getMonomerInfo());
    for atom in read_parts.atoms() {
        let Some(info) = atom.pdb_residue_info() else {
            continue;
        };

        // RDKit✔️✔️:       while (begin != end) {
        // RDKit✔️✔️:         if (mol.getAtomWithIdx(*begin)->getAtomicNum() == 1) {
        for neighbor in context.neighbors(atom.id()) {
            let neighbor_atom = AtomId::new(neighbor.atom_index);
            let Ok(neighbor_data) = context.atom(neighbor_atom) else {
                continue;
            };
            if neighbor_data.atomic_number() != 1 {
                continue;
            }
            // RDKit✔️✔️:           ++current_h_id;
            // COSMolKit increments inside info_for_new_hydrogen when it records a typed update.
            // RDKit✔️✔️:           if (h_info && h_info->getMonomerType() == PDBRESIDUE) { continue; }
            if neighbor_data.pdb_residue_info().is_some() || planned_existing_hydrogens.contains(&neighbor_atom) {
                residue_assignment.current_h_id += 1;
                continue;
            }
            // RDKit✔️✔️:           AtomPDBResidueInfo *newInfo = new AtomPDBResidueInfo(...);
            // RDKit✔️✔️:           mol.getAtomWithIdx(*begin)->setMonomerInfo(newInfo);
            // COSMolKit records the same typed atom-state update for the operation body.
            assignment.atom_pdb_residue_info_updates.push(AtomPdbResidueInfoUpdate {
                atom: neighbor_atom,
                pdb_residue_info: residue_assignment.info_for_new_hydrogen(info),
            });
            planned_existing_hydrogens.insert(neighbor_atom);
        }

        if let Some(new_hydrogen_indices) = new_hydrogens_by_heavy.get(&atom.id()) {
            for hydrogen_index in new_hydrogen_indices {
                // RDKit✔️✔️:           AtomPDBResidueInfo *newInfo = new AtomPDBResidueInfo(...);
                // COSMolKit stores the modeled PDB residue info subset as typed atom state.
                assignment.hydrogens_to_add[*hydrogen_index].pdb_residue_info =
                    Some(residue_assignment.info_for_new_hydrogen(info));
            }
        }
    }
    // END RDKIT CPP FUNCTION AssignHsResidueInfo
}

fn hydrogen_residue_label(ordinal: usize) -> String {
    let mut label = ordinal.to_string();
    if label.len() > 3 {
        label = label[label.len() - 3..].to_string();
    }
    while label.len() < 3 {
        label.push(' ');
    }
    let label = format!("H{label}");
    format!("{}{}", &label[3..4], &label[0..3])
}

fn isotopic_hydrogen_property(atom: &Atom) -> Vec<u16> {
    // BEGIN RDKIT CPP BODY MolOps::addHs _isotopicHs replay
    // RDKit✔️✔️:     std::vector<unsigned int> isoHs;
    // RDKit✔️❌:     if (newAt->getPropIfPresent(common_properties::_isotopicHs, isoHs)) {
    // RDKit✔️❌:       newAt->clearProp(common_properties::_isotopicHs);
    // RDKit✔️✔️:     }
    atom.tracked_isotopic_hydrogens().to_vec()
    // END RDKIT CPP BODY MolOps::addHs _isotopicHs replay
}

fn is_add_hs_query_atom(context: &RemoveHsPredicateContext<'_>, atom: AtomId) -> Result<bool, AddHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION isQueryAtom
    // RDKit✔️✔️: bool isQueryAtom(const RWMol &mol, const Atom &atom) {
    // RDKit✔️✔️:   if (atom.hasQuery()) { return true; }
    // RDKit✔️✔️:   for (const auto bnd : mol.atomBonds(&atom)) {
    // RDKit✔️✔️:     if (bnd->hasQuery()) { return true; }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    let atom_data = context.atom(atom).map_err(|err| match err {
        RemoveHydrogensError::AtomOutOfRange { atom, .. } => {
            AddHydrogensError::UnsupportedAtom { atom, atomic_number: 0 }
        }
        _ => AddHydrogensError::ProtocolDebt {
            branch: "MolOps::addHs/isQueryAtom",
            reason: "unexpected remove-H context error while checking AddHs query state",
        },
    })?;
    if atom_data.query().is_some() {
        return Ok(true);
    }
    for neighbor in context.neighbors(atom) {
        let bond = context
            .bond(neighbor.bond)
            .map_err(|_| AddHydrogensError::ProtocolDebt {
                branch: "MolOps::addHs/isQueryAtom",
                reason: "bond lookup failed while checking AddHs query state",
            })?;
        if bond.query().is_some() {
            return Ok(true);
        }
    }
    Ok(false)
    // END RDKIT CPP FUNCTION isQueryAtom
}

pub(crate) fn remove_hs_assignment(
    read_parts: MoleculeReadParts<'_>,
    params: &RemoveHsParams,
    sanitize: bool,
) -> Result<RemoveHsAssignment, RemoveHydrogensError> {
    remove_hs_assignment_inner(read_parts, params, sanitize, true)
}

fn remove_hs_assignment_inner(
    read_parts: MoleculeReadParts<'_>,
    params: &RemoveHsParams,
    sanitize: bool,
    allow_isotope_prepass: bool,
) -> Result<RemoveHsAssignment, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs(RWMol&, const RemoveHsParameters&, bool)
    // RDKit✔️❌: void removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize) {
    // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit✔️✔️:     // if there are any non-isotopic Hs remove them first
    // RDKit✔️✔️:     // to make sure chirality is preserved
    // RDKit✔️✔️:     bool needRemoveHs = false;
    // RDKit✔️✔️:     for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:       if (atom->getAtomicNum() == 1 && atom->getIsotope() == 0) {
    // RDKit✔️✔️:         needRemoveHs = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (needRemoveHs) {
    // RDKit✔️✔️:       RemoveHsParameters psCopy(ps);
    // RDKit✔️✔️:       psCopy.removeAndTrackIsotopes = false;
    // RDKit✔️✔️:       psCopy.removeIsotopes = false;
    // RDKit✔️❌:       removeHs(mol, psCopy, false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit✔️❌:     for (const auto &pair : getIsoMap(mol)) {
    // RDKit✔️❌:       mol.getAtomWithIdx(pair.first)
    // RDKit✔️❌:           ->setProp(common_properties::_isotopicHs, pair.second);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (shouldRemoveH(mol, atom, ps)) {
    // RDKit✔️✔️:       atomsToRemove.set(atom->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }  // end of the loop over atoms
    // RDKit✔️✔️:
    // RDKit✔️❌:   // Once we know which H atoms would be removed, filter out those that
    // RDKit✔️❌:   // would cause any SGroups to become empty
    // RDKit✔️❌:   if (ps.removeInSGroups) {
    // RDKit✔️❌:     filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // now that we know which atoms need to be removed, go ahead and remove them
    // RDKit✔️❌:   // NOTE: there's too much complexity around stereochemistry here
    // RDKit✔️❌:   // to be able to safely use batch editing.
    // RDKit✔️❌:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    // RDKit✔️❌:     if (atomsToRemove[idx]) {
    // RDKit✔️❌:       molRemoveH(mol, idx, ps.updateExplicitCount);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol.clearComputedProps(true);
    // RDKit✔️❌:   //
    // RDKit✔️❌:   //  If we didn't only remove implicit Hs, which are guaranteed to
    // RDKit✔️❌:   //  be the highest numbered atoms, we may have altered atom indices.
    // RDKit✔️❌:   //  This can screw up derived properties (such as ring members), so
    // RDKit✔️❌:   //  do some checks:
    // RDKit✔️❌:   //
    // RDKit❗✔️:   if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    // RDKit❗✔️:     sanitizeMol(mol);
    // RDKit❗✔️:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // if we removed Hs and any chiral atoms now have more than 1 explict H,
    // RDKit✔️❌:   // remove those
    // RDKit✔️❌:   if (!atomsToRemove.empty()) {
    // RDKit✔️❌:     for (auto atom : mol.atoms()) {
    // RDKit✔️❌:       if (!atom->getNoImplicit() &&
    // RDKit✔️❌:           atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:         unsigned int numExplicitHs = atom->getNumExplicitHs();
    // RDKit✔️❌:         if (numExplicitHs > 1) {
    // RDKit✔️❌:           atom->setNumExplicitHs(0);
    // RDKit✔️❌:           atom->updatePropertyCache(false);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::removeHs(RWMol&, const RemoveHsParameters&, bool)
    // BEGIN RDKIT CPP BODY MolOps::removeHs assignment planner
    if params.remove_and_track_isotopes && allow_isotope_prepass {
        // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) {
        // RDKit✔️✔️:     bool needRemoveHs = false;
        let need_remove_hs = read_parts
            .atoms()
            .iter()
            .any(|atom| atom.atomic_number() == 1 && atom.isotope().is_none());
        if need_remove_hs {
            // RDKit✔️✔️:       RemoveHsParameters psCopy(ps);
            // RDKit✔️✔️:       psCopy.removeAndTrackIsotopes = false;
            // RDKit✔️✔️:       psCopy.removeIsotopes = false;
            // RDKit✔️✔️:       removeHs(mol, psCopy, false);
            let first_pass_params = RemoveHsParams {
                remove_and_track_isotopes: false,
                remove_isotopes: false,
                ..params.clone()
            };
            let first_pass = remove_hs_assignment_inner(read_parts, &first_pass_params, false, false)?;
            let (first_pass_molecule, first_pass_mapping) =
                remove_hs_recursive_planning_molecule(read_parts, &first_pass)?;
            let second_pass = remove_hs_assignment_inner(
                MoleculeReadParts::from_molecule(&first_pass_molecule),
                params,
                sanitize,
                false,
            )?;
            let second_pass = remap_remove_hs_assignment_to_original(
                second_pass,
                &first_pass_mapping,
                "MolOps::removeHs/removeAndTrackIsotopes",
            )?;
            return Ok(merge_remove_hs_assignments(first_pass, second_pass));
        }
    }

    let atoms_with_isotopic_hydrogen_property = if params.remove_and_track_isotopes {
        read_parts
            .atoms()
            .iter()
            .filter(|atom| !atom.tracked_isotopic_hydrogens().is_empty())
            .map(Atom::id)
            .collect()
    } else {
        Vec::new()
    };
    let isotopic_hydrogens = if params.remove_and_track_isotopes {
        tracked_isotopic_hydrogens(read_parts)?
    } else {
        Vec::new()
    };

    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (shouldRemoveH(mol, atom, ps)) {
    // RDKit✔️✔️:       atomsToRemove.set(atom->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let context = RemoveHsPredicateContext::new(read_parts);
    let mut atoms_to_remove = Vec::new();
    for atom in read_parts.atoms() {
        if should_remove_h(&context, atom.id(), params)? {
            atoms_to_remove.push(atom.id());
        }
    }

    // RDKit✔️❌:   if (ps.removeInSGroups) { filter_sgroup_emptying_hydrogens(mol, atomsToRemove); }
    // COSMolKit records the filtered atom-removal plan here; mutation is applied later through OpParts.
    if params.remove_in_sgroups {
        filter_sgroup_emptying_hydrogens_for_assignment(read_parts, &mut atoms_to_remove);
    }

    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(false);
    // RDKit✔️❌:   }
    // The assignment carries those per-atom cache values across the later
    // compacting topology mapping. This mirrors Atom's value fields; they are
    // not computed RDProps and therefore survive mol.clearComputedProps(true).
    let (atom_valence_cache_updates, cached_total_valence) = remove_hs_property_cache_update(read_parts)?;
    let mut assignment = RemoveHsAssignment {
        atoms_to_remove,
        atom_valence_cache_updates,
        sanitize_after_removal: false,
        clear_computed_properties: true,
        ..RemoveHsAssignment::default()
    };

    // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) { ... getIsoMap(mol) ... setProp(_isotopicHs, pair.second); }
    // COSMolKit records the typed isotope property update and applies it through OpParts.
    for atom in atoms_with_isotopic_hydrogen_property {
        assignment
            .atom_property_updates
            .push(AtomPropertyUpdate::ClearIsotopicHydrogens { atom });
    }
    for (atom, isotopes) in isotopic_hydrogens {
        assignment
            .atom_property_updates
            .push(AtomPropertyUpdate::SetIsotopicHydrogens { atom, isotopes });
    }

    // RDKit✔️❌:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    // RDKit✔️❌:     if (atomsToRemove[idx]) { molRemoveH(mol, idx, ps.updateExplicitCount); }
    // RDKit✔️❌:   }
    // COSMolKit does not remove atoms here; it records per-H side effects for the operation body.
    let atoms_to_remove = assignment.atoms_to_remove.clone();
    for atom in atoms_to_remove {
        mol_remove_h_assignment(read_parts, atom, params, Some(&cached_total_valence), &mut assignment)?;
    }

    // RDKit✔️✔️: ROMol *removeHs(const ROMol &mol, const RemoveHsParameters &ps, bool sanitize) {
    // RDKit✔️✔️:   auto *res = new RWMol(mol);
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     removeHs(*res, ps, sanitize);
    // RDKit✔️✔️:   } catch (const MolSanitizeException &) {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     throw;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<ROMol *>(res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: python::def("RemoveHs",
    // RDKit✔️✔️:             (ROMol * (*)(const ROMol &, const MolOps::RemoveHsParameters &, bool)) &
    // RDKit✔️✔️:                 MolOps::removeHs,
    // RDKit✔️✔️:             (python::arg("mol"), python::arg("params"),
    // RDKit✔️✔️:              python::arg("sanitize") = true),
    // RDKit✔️❌:   if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) { sanitizeMol(mol); }
    //
    // `atomsToRemove` is a boost::dynamic_bitset sized to mol.getNumAtoms();
    // empty() is false for every non-empty molecule, even when no hydrogens
    // were selected for removal.
    if sanitize && params.remove_nonimplicit && !read_parts.atoms().is_empty() {
        assignment.sanitize_after_removal = true;
    }

    // RDKit✔️❌:   if (!atomsToRemove.empty()) { ... chiral atom explicit H cleanup ... }
    // COSMolKit records cleanup as an atom update; the operation body applies it through OpParts.
    if !assignment.atoms_to_remove.is_empty() {
        for atom in read_parts.atoms() {
            if !assignment.atoms_to_remove.contains(&atom.id())
                && !atom.no_implicit()
                && atom.chiral_tag() != ChiralTag::Unspecified
                && atom.explicit_hydrogens() > 1
            {
                assignment
                    .atom_property_updates
                    .push(AtomPropertyUpdate::ClearExcessChiralExplicitHydrogens { atom: atom.id() });
            }
        }
    }

    // END RDKIT CPP BODY MolOps::removeHs assignment planner
    Ok(assignment)
}

fn merge_remove_hs_assignments(mut first: RemoveHsAssignment, second: RemoveHsAssignment) -> RemoveHsAssignment {
    first.atoms_to_remove.extend(second.atoms_to_remove);
    first.atom_valence_cache_updates = second.atom_valence_cache_updates;
    first
        .atom_explicit_hydrogen_updates
        .extend(second.atom_explicit_hydrogen_updates);
    first.atom_chirality_inversions.extend(second.atom_chirality_inversions);
    first.atom_property_updates.extend(second.atom_property_updates);
    first.bond_direction_updates.extend(second.bond_direction_updates);
    first.bond_stereo_updates.extend(second.bond_stereo_updates);
    first
        .bond_stereo_atom_replacements
        .extend(second.bond_stereo_atom_replacements);
    first.sgroup_updates.extend(second.sgroup_updates);
    first.sanitize_after_removal = second.sanitize_after_removal;
    first.clear_computed_properties |= second.clear_computed_properties;
    first
}

fn remove_hs_recursive_planning_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &RemoveHsAssignment,
) -> Result<(Molecule, crate::molecule::TopologyMapping), RemoveHydrogensError> {
    // BEGIN RDKIT CPP BODY MolOps::removeHs recursive prepass mutation model
    // RDKit✔️❌:       removeHs(mol, psCopy, false);
    // COSMolKit does not mutate the source molecule here. This creates the
    // post-prepass value used only to plan the second removeAndTrackIsotopes phase.
    let mut topology = read_parts.topology().clone();
    for update in &assignment.atom_explicit_hydrogen_updates {
        let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
            return Err(RemoveHydrogensError::AtomOutOfRange {
                atom: update.atom,
                atom_count: topology.atoms.len(),
            });
        };
        atom.set_explicit_hydrogens(update.explicit_hydrogens);
    }
    for atom_id in &assignment.atom_chirality_inversions {
        let Some(atom) = topology.atoms.get_mut(atom_id.index()) else {
            return Err(RemoveHydrogensError::AtomOutOfRange {
                atom: *atom_id,
                atom_count: topology.atoms.len(),
            });
        };
        atom.set_chiral_tag(match atom.chiral_tag() {
            ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
            ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
            other => other,
        });
    }
    for update in &assignment.atom_property_updates {
        match update {
            AtomPropertyUpdate::SetUnknownStereo { atom } => {
                let Some(atom_data) = topology.atoms.get_mut(atom.index()) else {
                    return Err(RemoveHydrogensError::AtomOutOfRange {
                        atom: *atom,
                        atom_count: topology.atoms.len(),
                    });
                };
                atom_data.set_unknown_stereo(true);
            }
            AtomPropertyUpdate::SetIsotopicHydrogens { atom, isotopes } => {
                let Some(atom_data) = topology.atoms.get_mut(atom.index()) else {
                    return Err(RemoveHydrogensError::AtomOutOfRange {
                        atom: *atom,
                        atom_count: topology.atoms.len(),
                    });
                };
                atom_data.set_tracked_isotopic_hydrogens(isotopes.clone());
            }
            AtomPropertyUpdate::ClearIsotopicHydrogens { atom } => {
                let Some(atom_data) = topology.atoms.get_mut(atom.index()) else {
                    return Err(RemoveHydrogensError::AtomOutOfRange {
                        atom: *atom,
                        atom_count: topology.atoms.len(),
                    });
                };
                atom_data.set_tracked_isotopic_hydrogens(Vec::new());
            }
            AtomPropertyUpdate::ClearExcessChiralExplicitHydrogens { atom } => {
                let Some(atom_data) = topology.atoms.get_mut(atom.index()) else {
                    return Err(RemoveHydrogensError::AtomOutOfRange {
                        atom: *atom,
                        atom_count: topology.atoms.len(),
                    });
                };
                atom_data.set_explicit_hydrogens(0);
            }
        }
    }
    for update in &assignment.bond_direction_updates {
        let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
            return Err(RemoveHydrogensError::BondOutOfRange {
                bond: update.bond,
                bond_count: topology.bonds.len(),
            });
        };
        bond.set_direction(update.direction);
    }
    for update in &assignment.bond_stereo_updates {
        let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
            return Err(RemoveHydrogensError::BondOutOfRange {
                bond: update.bond,
                bond_count: topology.bonds.len(),
            });
        };
        bond.set_stereo_atoms(update.stereo_atoms);
        bond.set_stereo(update.stereo);
    }
    for update in &assignment.bond_stereo_atom_replacements {
        let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
            return Err(RemoveHydrogensError::BondOutOfRange {
                bond: update.bond,
                bond_count: topology.bonds.len(),
            });
        };
        let Some([left, right]) = bond.stereo_atoms() else {
            continue;
        };
        bond.set_stereo_atoms(Some([
            if left == update.old_atom { update.new_atom } else { left },
            if right == update.old_atom {
                update.new_atom
            } else {
                right
            },
        ]));
    }
    for update in &assignment.sgroup_updates {
        for substance_group in &mut topology.substance_groups {
            match update {
                SGroupRemoveHsUpdate::RemoveBond { bond } => {
                    substance_group.remove_bond(*bond);
                }
                SGroupRemoveHsUpdate::RemoveAtom { atom } => {
                    substance_group.remove_atom(*atom);
                }
                SGroupRemoveHsUpdate::RemoveParentAtom { atom } => {
                    substance_group.remove_parent_atom(*atom);
                }
                SGroupRemoveHsUpdate::ClearAttachPointLeavingAtom { atom } => {
                    substance_group.clear_attach_point_leaving_atom(*atom);
                }
            }
        }
    }

    let mapping = topology.remove_atoms_with_mapping(&assignment.atoms_to_remove);
    let mut coordinates = read_parts.coordinates().clone();
    coordinates.remap_topology(&mapping.retained_atom_indices());
    let planned = Molecule::from_blocks(topology, coordinates, read_parts.properties().clone()).map_err(|_| {
        RemoveHydrogensError::ProtocolDebt {
            branch: "MolOps::removeHs/removeAndTrackIsotopes",
            reason: "recursive prepass produced an invalid temporary molecule",
        }
    })?;
    Ok((planned, mapping))
    // END RDKIT CPP BODY MolOps::removeHs recursive prepass mutation model
}

fn remap_remove_hs_assignment_to_original(
    mut assignment: RemoveHsAssignment,
    mapping: &crate::molecule::TopologyMapping,
    branch: &'static str,
) -> Result<RemoveHsAssignment, RemoveHydrogensError> {
    for update in &mut assignment.atom_valence_cache_updates {
        update.atom = original_atom_from_prepass(mapping, update.atom, branch)?;
    }
    for atom in &mut assignment.atoms_to_remove {
        *atom = original_atom_from_prepass(mapping, *atom, branch)?;
    }
    for update in &mut assignment.atom_explicit_hydrogen_updates {
        update.atom = original_atom_from_prepass(mapping, update.atom, branch)?;
    }
    for atom in &mut assignment.atom_chirality_inversions {
        *atom = original_atom_from_prepass(mapping, *atom, branch)?;
    }
    for update in &mut assignment.atom_property_updates {
        match update {
            AtomPropertyUpdate::SetUnknownStereo { atom }
            | AtomPropertyUpdate::ClearIsotopicHydrogens { atom }
            | AtomPropertyUpdate::ClearExcessChiralExplicitHydrogens { atom } => {
                *atom = original_atom_from_prepass(mapping, *atom, branch)?;
            }
            AtomPropertyUpdate::SetIsotopicHydrogens { atom, .. } => {
                *atom = original_atom_from_prepass(mapping, *atom, branch)?;
            }
        }
    }
    for update in &mut assignment.bond_direction_updates {
        update.bond = original_bond_from_prepass(mapping, update.bond, branch)?;
    }
    for update in &mut assignment.bond_stereo_updates {
        update.bond = original_bond_from_prepass(mapping, update.bond, branch)?;
        update.stereo_atoms = match update.stereo_atoms {
            Some([left, right]) => Some([
                original_atom_from_prepass(mapping, left, branch)?,
                original_atom_from_prepass(mapping, right, branch)?,
            ]),
            None => None,
        };
    }
    for update in &mut assignment.bond_stereo_atom_replacements {
        update.bond = original_bond_from_prepass(mapping, update.bond, branch)?;
        update.old_atom = original_atom_from_prepass(mapping, update.old_atom, branch)?;
        update.new_atom = original_atom_from_prepass(mapping, update.new_atom, branch)?;
    }
    for update in &mut assignment.sgroup_updates {
        match update {
            SGroupRemoveHsUpdate::RemoveBond { bond } => {
                *bond = original_bond_from_prepass(mapping, *bond, branch)?;
            }
            SGroupRemoveHsUpdate::RemoveAtom { atom }
            | SGroupRemoveHsUpdate::RemoveParentAtom { atom }
            | SGroupRemoveHsUpdate::ClearAttachPointLeavingAtom { atom } => {
                *atom = original_atom_from_prepass(mapping, *atom, branch)?;
            }
        }
    }
    Ok(assignment)
}

fn original_atom_from_prepass(
    mapping: &crate::molecule::TopologyMapping,
    atom: AtomId,
    branch: &'static str,
) -> Result<AtomId, RemoveHydrogensError> {
    mapping
        .atoms()
        .new_to_old()
        .get(atom.index())
        .and_then(|mapped| *mapped)
        .ok_or(RemoveHydrogensError::ProtocolDebt {
            branch,
            reason: "recursive remove-H prepass atom mapping did not preserve a required atom",
        })
}

fn original_bond_from_prepass(
    mapping: &crate::molecule::TopologyMapping,
    bond: BondId,
    branch: &'static str,
) -> Result<BondId, RemoveHydrogensError> {
    mapping
        .bonds()
        .new_to_old()
        .get(bond.index())
        .and_then(|mapped| *mapped)
        .ok_or(RemoveHydrogensError::ProtocolDebt {
            branch,
            reason: "recursive remove-H prepass bond mapping did not preserve a required bond",
        })
}

#[allow(dead_code)]
struct RemoveHsPredicateContext<'a> {
    read_parts: MoleculeReadParts<'a>,
    adjacency: AdjacencyList,
}

#[allow(dead_code)]
impl<'a> RemoveHsPredicateContext<'a> {
    fn new(read_parts: MoleculeReadParts<'a>) -> Self {
        Self {
            read_parts,
            adjacency: AdjacencyList::from_topology(read_parts.num_atoms(), read_parts.bonds()),
        }
    }

    fn atom(&self, atom: AtomId) -> Result<&'a Atom, RemoveHydrogensError> {
        self.read_parts.atom(atom).ok_or(RemoveHydrogensError::AtomOutOfRange {
            atom,
            atom_count: self.read_parts.num_atoms(),
        })
    }

    fn bond(&self, bond: BondId) -> Result<&'a Bond, RemoveHydrogensError> {
        self.read_parts
            .bonds()
            .get(bond.index())
            .ok_or(RemoveHydrogensError::BondOutOfRange {
                bond,
                bond_count: self.read_parts.num_bonds(),
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
        || substance_group
            .attach_points()
            .iter()
            .any(|attach_point| attach_point.atom == atom || attach_point.leaving_atom == Some(atom))
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

fn filter_sgroup_emptying_hydrogens_for_assignment(
    read_parts: MoleculeReadParts<'_>,
    atoms_to_remove: &mut Vec<AtomId>,
) {
    // BEGIN RDKIT CPP FUNCTION filter_sgroup_emptying_hydrogens
    // RDKit✔️❌: Once candidate H atoms are known, filter out removals that would empty an SGroup.
    // COSMolKit updates only the assignment bitset/vector here; SGroup mutation remains OpParts-owned.
    let original_atoms_to_remove = atoms_to_remove.clone();
    atoms_to_remove.retain(|atom| {
        read_parts.substance_groups().iter().all(|sgroup| {
            let atoms = sgroup.atoms();
            atoms.is_empty()
                || atoms
                    .iter()
                    .any(|member| member != atom && !original_atoms_to_remove.contains(member))
        })
    });
    // END RDKIT CPP FUNCTION filter_sgroup_emptying_hydrogens
}

fn remove_hs_property_cache_update(
    read_parts: MoleculeReadParts<'_>,
) -> Result<(Vec<AtomValenceCacheUpdate>, Vec<i32>), RemoveHydrogensError> {
    // BEGIN RDKIT CPP BODY MolOps::removeHs property-cache preparation
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(false);
    // RDKit✔️❌:   }
    let assignment =
        read_parts
            .assign_valence_with_options(ValenceModel::RdkitLike, false)
            .map_err(|_| RemoveHydrogensError::ProtocolDebt {
                branch: "MolOps::removeHs/updatePropertyCache",
                reason: "remove-H requires a modeled RDKit-like property-cache update before heavy-atom explicit H decisions",
            })?;
    let mut updates = Vec::with_capacity(read_parts.num_atoms());
    let mut total_valence = Vec::with_capacity(read_parts.num_atoms());
    for (index, (explicit_valence, implicit_hydrogens)) in assignment
        .explicit_valence
        .into_iter()
        .zip(assignment.implicit_hydrogens)
        .enumerate()
    {
        updates.push(AtomValenceCacheUpdate {
            atom: AtomId::new(index),
            explicit_valence,
            implicit_hydrogens,
        });
        total_valence.push(explicit_valence + implicit_hydrogens);
    }
    Ok((updates, total_valence))
    // END RDKIT CPP BODY MolOps::removeHs property-cache preparation
}

fn tracked_isotopic_hydrogens(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<(AtomId, Vec<u16>)>, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/AddHs.cpp :: getIsoMap
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     if (atom->hasProp(common_properties::_isotopicHs)) {
    // RDKit✔️❌:       atom->clearProp(common_properties::_isotopicHs);
    // RDKit✔️❌:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/AddHs.cpp :: getIsoMap
    // The operation adapter records clears before sets in `RemoveHsAssignment`.
    crate::source_port_helpers::rdkit_get_iso_map(
        read_parts.bonds().iter().map(|bond| (bond.begin(), bond.end())),
        |atom_id| {
            read_parts
                .atom(atom_id)
                .map(|atom| (atom.atomic_number(), atom.isotope()))
                .ok_or(RemoveHydrogensError::AtomOutOfRange {
                    atom: atom_id,
                    atom_count: read_parts.num_atoms(),
                })
        },
    )
}

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
    if !params.remove_isotopes && !params.remove_and_track_isotopes && atom_data.isotope().is_some() {
        return Ok(false);
    }

    // RDKit✔️✔️:   if (!ps.removeNonimplicit && !atom->hasProp(common_properties::isImplicit)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if !params.remove_nonimplicit && !atom_data.implicit_hydrogen() {
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
        for substance_group in context.read_parts.substance_groups() {
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
        for substance_group in context.read_parts.substance_groups() {
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
            if !params.remove_nontetrahedral_neighbors && has_non_tetrahedral_stereo(neighbor_data) {
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

fn mol_remove_h_assignment(
    read_parts: MoleculeReadParts<'_>,
    atom: AtomId,
    params: &RemoveHsParams,
    cached_total_valence: Option<&[i32]>,
    assignment: &mut RemoveHsAssignment,
) -> Result<(), RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION molRemoveH
    // RDKit✔️❌: void molRemoveH(RWMol &mol, unsigned int idx, bool updateExplicitCount) {
    // RDKit✔️✔️:   auto atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:   PRECONDITION(atom->getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");
    // RDKit✔️❌:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     Atom *heavyAtom = bond->getOtherAtom(atom);
    // RDKit✔️✔️:     int heavyAtomNum = heavyAtom->getAtomicNum();
    // RDKit✔️✔️:
    // RDKit✔️❌:     if (updateExplicitCount || heavyAtom->getNoImplicit() ||
    // RDKit✔️❌:         heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) { ... }
    // RDKit✔️❌:     else { ... aromatic/default-valence explicit-H preservation ... }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       ... getPerturbationOrder(...); invertChirality(); ...
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (heavyAtom->getDegree() == 2) { ... clear neighboring bond stereo ... }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (bond->getBondDir() == Bond::UNKNOWN &&
    // RDKit✔️❌:         bond->getBeginAtomIdx() == heavyAtom->getIdx()) { ... }
    // RDKit✔️❌:     else if (bond->getBondDir() == Bond::ENDDOWNRIGHT ||
    // RDKit✔️❌:                bond->getBondDir() == Bond::ENDUPRIGHT) { ... transfer dir ... }
    // RDKit✔️❌:     else { adjustStereoAtomsIfRequired(mol, atom, heavyAtom); }
    // RDKit✔️❌:
    // RDKit✔️❌:     for (auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️❌:       sg.removeBondWithIdx(bond->getIdx());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️❌:     sg.removeAtomWithIdx(idx);
    // RDKit✔️❌:     sg.removeParentAtomWithIdx(idx);
    // RDKit✔️❌:     for (auto &sap : sg.getAttachPoints()) {
    // RDKit✔️❌:       if (sap.lvIdx == static_cast<int>(idx)) { sap.lvIdx = -1; }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   bool clearProps = false;
    // RDKit✔️❌:   mol.removeAtom(atom, clearProps);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION molRemoveH
    // BEGIN RDKIT CPP BODY molRemoveH assignment planner
    // RDKit✔️✔️:   auto atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:   PRECONDITION(atom->getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");
    let context = RemoveHsPredicateContext::new(read_parts);
    let atom_data = context.atom(atom)?;
    if atom_data.atomic_number() != 1 {
        return Err(RemoveHydrogensError::UnsupportedHydrogen { atom });
    }

    // RDKit✔️❌:   for (const auto bond : mol.atomBonds(atom)) { ... }
    // COSMolKit records each incident-bond side effect as assignment data. Degree-zero
    // hydrogens naturally skip this loop and still fall through to atom/SGroup removal.
    let neighbors = context.neighbors(atom);
    for neighbor in neighbors {
        // RDKit✔️✔️:     Atom *heavyAtom = bond->getOtherAtom(atom);
        // RDKit✔️✔️:     int heavyAtomNum = heavyAtom->getAtomicNum();
        let hydrogen_bond = context.bond(neighbor.bond)?;
        let heavy_atom = if hydrogen_bond.begin() == atom {
            hydrogen_bond.end()
        } else {
            hydrogen_bond.begin()
        };
        let heavy_atom_data = context.atom(heavy_atom)?;

        // RDKit✔️❌:     if (heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) { ... invertChirality(); }
        // COSMolKit records the inversion; the operation body applies it through OpParts.
        if heavy_atom_data.chiral_tag() != ChiralTag::Unspecified {
            if should_invert_heavy_atom_chirality(&context, heavy_atom, hydrogen_bond.id())? {
                assignment.atom_chirality_inversions.push(heavy_atom);
            }
        }

        // RDKit✔️❌:     if (heavyAtom->getDegree() == 2) { ... clear neighboring bond stereo ... }
        // COSMolKit records this as BondStereoUpdate.
        if context.degree(heavy_atom) == 2
            && let Some(update) = clear_neighbor_stereo_assignment(&context, heavy_atom, hydrogen_bond.id())?
        {
            assignment.bond_stereo_updates.push(update);
        }

        // RDKit✔️❌:     if (bond->getBondDir() == Bond::UNKNOWN) { heavyAtom->setProp(...); }
        // RDKit✔️❌:     else if (bond->getBondDir() == ENDDOWNRIGHT || ENDUPRIGHT) { ... transfer dir ... }
        // RDKit✔️✔️:     else { adjustStereoAtomsIfRequired(mol, atom, heavyAtom); }
        match hydrogen_bond.direction() {
            BondDirection::None => {}
            BondDirection::Unknown => {
                if hydrogen_bond.begin() == heavy_atom {
                    assignment
                        .atom_property_updates
                        .push(AtomPropertyUpdate::SetUnknownStereo { atom: heavy_atom });
                }
            }
            BondDirection::EndDownRight | BondDirection::EndUpRight => {
                if let Some(update) =
                    copy_hydrogen_bond_direction_to_neighbor_assignment(&context, heavy_atom, hydrogen_bond)?
                {
                    assignment.bond_direction_updates.push(update);
                }
            }
            BondDirection::BeginWedge | BondDirection::BeginDash | BondDirection::EitherDouble => {}
        }

        // RDKit✔️❌:     adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
        // COSMolKit records stereo atom replacement and cis/trans flip as assignment updates.
        if let Some((replacement, stereo_update)) =
            adjust_stereo_atoms_if_required_assignment(&context, atom, heavy_atom)?
        {
            assignment.bond_stereo_updates.push(stereo_update);
            assignment.bond_stereo_atom_replacements.push(replacement);
        }

        // RDKit✔️❌:     if (updateExplicitCount || heavyAtom->getNoImplicit() || ...) {
        // RDKit✔️❌:       heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
        // RDKit✔️❌:     }
        // COSMolKit records this as AtomExplicitHydrogenUpdate.
        if should_increment_explicit_h_count(
            read_parts,
            heavy_atom,
            params.update_explicit_count,
            cached_total_valence,
        )? {
            increment_explicit_hydrogen_update(assignment, heavy_atom, heavy_atom_data.explicit_hydrogens());
        }

        // RDKit✔️❌:     for (auto &sg : getSubstanceGroups(mol)) { sg.removeBondWithIdx(...); }
        // COSMolKit records SGroup bond removal for the operation body.
        for substance_group in read_parts.substance_groups() {
            if substance_group.bonds().contains(&hydrogen_bond.id()) {
                assignment.sgroup_updates.push(SGroupRemoveHsUpdate::RemoveBond {
                    bond: hydrogen_bond.id(),
                });
            }
        }
    }

    // RDKit✔️❌:   for (auto &sg : getSubstanceGroups(mol)) { sg.removeAtomWithIdx(...); ... }
    // COSMolKit records SGroup changes for the operation body; topology compaction remaps SGroups through OpParts.
    for substance_group in read_parts.substance_groups() {
        if substance_group.atoms().contains(&atom) {
            assignment
                .sgroup_updates
                .push(SGroupRemoveHsUpdate::RemoveAtom { atom });
        }
        if substance_group.parent_atoms().contains(&atom) {
            assignment
                .sgroup_updates
                .push(SGroupRemoveHsUpdate::RemoveParentAtom { atom });
        }
        if substance_group
            .attach_points()
            .iter()
            .any(|attach_point| attach_point.leaving_atom == Some(atom))
        {
            assignment
                .sgroup_updates
                .push(SGroupRemoveHsUpdate::ClearAttachPointLeavingAtom { atom });
        }
    }

    // RDKit✔️❌:   mol.removeAtom(atom, clearProps);
    // COSMolKit records atom removal in RemoveHsAssignment; OpParts::remove_atoms applies it atomically.
    // END RDKIT CPP BODY molRemoveH assignment planner
    Ok(())
}

fn increment_explicit_hydrogen_update(
    assignment: &mut RemoveHsAssignment,
    atom: AtomId,
    original_explicit_hydrogens: u8,
) {
    if let Some(update) = assignment
        .atom_explicit_hydrogen_updates
        .iter_mut()
        .find(|update| update.atom == atom)
    {
        update.explicit_hydrogens = update.explicit_hydrogens.saturating_add(1);
        return;
    }
    assignment
        .atom_explicit_hydrogen_updates
        .push(AtomExplicitHydrogenUpdate {
            atom,
            explicit_hydrogens: original_explicit_hydrogens.saturating_add(1),
        });
}

fn should_invert_heavy_atom_chirality(
    context: &RemoveHsPredicateContext<'_>,
    heavy_atom: AtomId,
    hydrogen_bond: BondId,
) -> Result<bool, RemoveHydrogensError> {
    // BEGIN RDKIT CPP BODY Atom::getPerturbationOrder for molRemoveH
    // RDKit✔️✔️:       INT_LIST neighborIndices;
    // RDKit✔️✔️:       for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit✔️✔️:         if (nbnd->getIdx() != bond->getIdx()) { neighborIndices.push_back(nbnd->getIdx()); }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       neighborIndices.push_back(bond->getIdx());
    // RDKit✔️✔️:       int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
    let reference = incident_bond_ids(context, heavy_atom);
    let mut probe = reference
        .iter()
        .copied()
        .filter(|bond| *bond != hydrogen_bond)
        .collect::<Vec<_>>();
    probe.push(hydrogen_bond);
    Ok(count_swaps_to_interconvert(&probe, &reference)? % 2 == 1)
    // END RDKIT CPP BODY Atom::getPerturbationOrder for molRemoveH
}

fn clear_neighbor_stereo_assignment(
    context: &RemoveHsPredicateContext<'_>,
    heavy_atom: AtomId,
    hydrogen_bond: BondId,
) -> Result<Option<BondStereoUpdate>, RemoveHydrogensError> {
    // BEGIN RDKIT CPP BODY molRemoveH clear neighboring double-bond stereo
    // RDKit✔️✔️:     if (heavyAtom->getDegree() == 2) {
    // RDKit✔️✔️:       for (auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit✔️✔️:         if (nbnd != bond) {
    // RDKit✔️✔️:           if (nbnd->getStereo() > Bond::STEREOANY) {
    // RDKit✔️✔️:             nbnd->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️:             nbnd->getStereoAtoms().clear();
    // RDKit✔️✔️:           }
    for neighbor in context.neighbors(heavy_atom) {
        if neighbor.bond == hydrogen_bond {
            continue;
        }
        let bond = context.bond(neighbor.bond)?;
        if bond_stereo_greater_than_any(bond.stereo()) {
            return Ok(Some(BondStereoUpdate {
                bond: bond.id(),
                stereo: BondStereo::None,
                stereo_atoms: None,
            }));
        }
        break;
    }
    Ok(None)
    // END RDKIT CPP BODY molRemoveH clear neighboring double-bond stereo
}

fn copy_hydrogen_bond_direction_to_neighbor_assignment(
    context: &RemoveHsPredicateContext<'_>,
    heavy_atom: AtomId,
    hydrogen_bond: &Bond,
) -> Result<Option<BondDirectionUpdate>, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION copy_hydrogen_bond_direction_to_neighbor
    // RDKit✔️✔️:       bool foundADir = false;
    // RDKit✔️✔️:       Bond *oBond = nullptr;
    // RDKit✔️✔️:       for (...) {
    // RDKit✔️✔️:         if (nbnd->getIdx() != bond->getIdx() && nbnd->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:           if (nbnd->getBondDir() == Bond::NONE) { oBond = nbnd; } else { foundADir = true; }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!foundADir && oBond != nullptr) { ... setBondDir(...); }
    let mut found_direction = false;
    let mut other_single_bond = None;
    for neighbor in context.neighbors(heavy_atom) {
        if neighbor.bond == hydrogen_bond.id() {
            continue;
        }
        let bond = context.bond(neighbor.bond)?;
        if bond.order() != BondOrder::Single {
            continue;
        }
        if bond.direction() == BondDirection::None {
            other_single_bond = Some(bond);
        } else {
            found_direction = true;
        }
    }
    if found_direction {
        return Ok(None);
    }
    let Some(other_bond) = other_single_bond else {
        return Ok(None);
    };
    let mut direction = hydrogen_bond.direction();
    let flip = other_bond.begin() == heavy_atom && hydrogen_bond.begin() == heavy_atom;
    if flip {
        direction = opposite_bond_direction(direction);
    }
    Ok(Some(BondDirectionUpdate {
        bond: other_bond.id(),
        direction,
    }))
    // END RDKIT CPP FUNCTION copy_hydrogen_bond_direction_to_neighbor
}

fn adjust_stereo_atoms_if_required_assignment(
    context: &RemoveHsPredicateContext<'_>,
    atom: AtomId,
    heavy_atom: AtomId,
) -> Result<Option<(BondStereoAtomReplacement, BondStereoUpdate)>, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION adjustStereoAtomsIfRequired
    // RDKit✔️✔️:   if (heavyAtom->getDegree() == 2) { return false; }
    if context.degree(heavy_atom) == 2 {
        return Ok(None);
    }

    // RDKit✔️✔️:   for (bond around heavyAtom) {
    // RDKit✔️✔️:     if (bond is double and CIS/TRANS and atom is a stereo atom) {
    // RDKit✔️✔️:       replace atom with another neighbor and flip CIS/TRANS;
    for neighbor_ref in context.neighbors(heavy_atom) {
        let bond = context.bond(neighbor_ref.bond)?;
        if bond.order() != BondOrder::Double || !bond_stereo_greater_than_any(bond.stereo()) {
            continue;
        }
        let Some(stereo_atoms) = bond.stereo_atoms() else {
            continue;
        };
        if stereo_atoms[0] != atom && stereo_atoms[1] != atom {
            continue;
        }
        let double_neighbor = if bond.begin() == heavy_atom {
            bond.end()
        } else {
            bond.begin()
        };
        for replacement_ref in context.neighbors(heavy_atom) {
            let replacement = AtomId::new(replacement_ref.atom_index);
            if replacement == double_neighbor || replacement == atom {
                continue;
            }
            let stereo = match bond.stereo() {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                other => other,
            };
            return Ok(Some((
                BondStereoAtomReplacement {
                    bond: bond.id(),
                    old_atom: atom,
                    new_atom: replacement,
                },
                BondStereoUpdate {
                    bond: bond.id(),
                    stereo,
                    stereo_atoms: bond.stereo_atoms(),
                },
            )));
        }
    }
    Ok(None)
    // END RDKIT CPP FUNCTION adjustStereoAtomsIfRequired
}

fn incident_bond_ids(context: &RemoveHsPredicateContext<'_>, atom: AtomId) -> Vec<BondId> {
    context
        .read_parts
        .bonds()
        .iter()
        .filter_map(|bond| (bond.begin() == atom || bond.end() == atom).then_some(bond.id()))
        .collect()
}

fn count_swaps_to_interconvert(probe: &[BondId], reference: &[BondId]) -> Result<usize, RemoveHydrogensError> {
    crate::source_port_helpers::count_swaps_to_interconvert(reference, probe).map_err(|error| {
        let (branch, reason) = match error {
            crate::source_port_helpers::CountSwapsError::SizeMismatch => (
                "Atom::getPerturbationOrder/size-mismatch",
                "chirality perturbation order received mismatched bond-ordering lengths",
            ),
            crate::source_port_helpers::CountSwapsError::MissingProbeElement => (
                "Atom::getPerturbationOrder/missing-probe-element",
                "chirality perturbation order could not match RDKit bond-ordering elements",
            ),
        };
        RemoveHydrogensError::ProtocolDebt { branch, reason }
    })
}

fn opposite_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        other => other,
    }
}

fn should_increment_explicit_h_count(
    read_parts: MoleculeReadParts<'_>,
    heavy_atom: AtomId,
    update_explicit_count: bool,
    cached_total_valence: Option<&[i32]>,
) -> Result<bool, RemoveHydrogensError> {
    // BEGIN RDKIT CPP BODY molRemoveH explicit-H-count decision
    // RDKit✔️✔️:     if (updateExplicitCount || heavyAtom->getNoImplicit() ||
    // RDKit✔️✔️:         heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) { ... }
    let heavy_atom_data = read_parts
        .atom(heavy_atom)
        .ok_or(RemoveHydrogensError::AtomOutOfRange {
            atom: heavy_atom,
            atom_count: read_parts.num_atoms(),
        })?;
    if update_explicit_count || heavy_atom_data.no_implicit() || heavy_atom_data.chiral_tag() != ChiralTag::Unspecified
    {
        return Ok(true);
    }

    // RDKit✔️✔️:       const INT_VECT &defaultVs =
    // RDKit✔️✔️:           PeriodicTable::getTable()->getValenceList(heavyAtomNum);
    // RDKit✔️✔️:       if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
    // RDKit✔️✔️:             may_need_extra_H(mol, heavyAtom)) &&
    // RDKit✔️✔️:            isAromaticAtom(*heavyAtom)) ||
    // RDKit✔️✔️:           (std::find(defaultVs.begin() + 1, defaultVs.end(),
    // RDKit✔️✔️:                      heavyAtom->getTotalValence()) != defaultVs.end())) { ... }
    let total_valence = cached_total_valence
        .and_then(|values| values.get(heavy_atom.index()))
        .copied()
        .ok_or(RemoveHydrogensError::ProtocolDebt {
            branch: "molRemoveH/explicit-H-count",
            reason: "cached total valence is required before deciding whether to preserve removed H as explicit count",
        })?;
    let default_valences =
        crate::rdkit_valence_list(heavy_atom_data.atomic_number()).map_err(|_| RemoveHydrogensError::ProtocolDebt {
            branch: "molRemoveH/default-valence-list",
            reason: "RDKit valence table lookup failed while deciding removed-H explicit count",
        })?;
    let non_default_valence =
        default_valences.is_some_and(|values| values.iter().skip(1).any(|value| *value == total_valence));
    Ok(((heavy_atom_data.atomic_number() == 7
        || heavy_atom_data.atomic_number() == 15
        || may_need_extra_h(read_parts, heavy_atom, cached_total_valence)?)
        && heavy_atom_data.is_aromatic())
        || non_default_valence)
    // END RDKIT CPP BODY molRemoveH explicit-H-count decision
}

fn may_need_extra_h(
    read_parts: MoleculeReadParts<'_>,
    atom: AtomId,
    cached_total_valence: Option<&[i32]>,
) -> Result<bool, RemoveHydrogensError> {
    // BEGIN RDKIT CPP FUNCTION may_need_extra_H
    // RDKit✔️✔️: aromatic N/P and special aromatic-valence branch used by molRemoveH.
    let context = RemoveHsPredicateContext::new(read_parts);
    let mut single_bonds = 0usize;
    let mut aromatic_bonds = 0usize;
    for neighbor in context.neighbors(atom) {
        match context.bond(neighbor.bond)?.order() {
            BondOrder::Single => single_bonds += 1,
            BondOrder::Aromatic => aromatic_bonds += 1,
            _ => return Ok(false),
        }
    }
    let total_valence = cached_total_valence
        .and_then(|values| values.get(atom.index()))
        .copied()
        .ok_or(RemoveHydrogensError::ProtocolDebt {
            branch: "molRemoveH/may_need_extra_H",
            reason: "cached total valence is required before applying aromatic extra-H logic",
        })?;
    Ok(single_bonds == 1 && aromatic_bonds == 2 && total_valence == 3)
    // END RDKIT CPP FUNCTION may_need_extra_H
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomQueryPredicate, AtomSpec, BondSpec, Element, MoleculeBuilder, QueryNode, SGroupAttachPoint, SGroupCState,
        SubstanceGroupId, SubstanceGroupKind,
    };

    fn read(molecule: &Molecule) -> MoleculeReadParts<'_> {
        MoleculeReadParts::from_molecule(molecule)
    }

    #[test]
    fn shared_get_iso_map_remove_hs_adapter_preserves_atom_range_error() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let deuterium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        let bond = builder
            .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let missing = AtomId::new(7);
        molecule.topology_block_mut().bonds[bond.index()].set_endpoints(carbon, missing);

        let error = tracked_isotopic_hydrogens(read(&molecule)).unwrap_err();

        assert_eq!(
            error,
            RemoveHydrogensError::AtomOutOfRange {
                atom: missing,
                atom_count: 2,
            }
        );
    }

    #[test]
    fn shared_count_swaps_remove_hs_preserves_protocol_debt_mapping() {
        let size_error = count_swaps_to_interconvert(&[BondId::new(0)], &[BondId::new(0), BondId::new(1)]).unwrap_err();
        assert!(matches!(
            size_error,
            RemoveHydrogensError::ProtocolDebt {
                branch: "Atom::getPerturbationOrder/size-mismatch",
                ..
            }
        ));

        let missing_error =
            count_swaps_to_interconvert(&[BondId::new(0), BondId::new(2)], &[BondId::new(0), BondId::new(1)])
                .unwrap_err();
        assert!(matches!(
            missing_error,
            RemoveHydrogensError::ProtocolDebt {
                branch: "Atom::getPerturbationOrder/missing-probe-element",
                ..
            }
        ));
    }

    #[test]
    fn remove_hs_assignment_returns_noop_for_empty_molecule() {
        let molecule = Molecule::default();
        let assignment = remove_hs_assignment(read(&molecule), &RemoveHsParams::default(), true).unwrap();

        assert!(assignment.atoms_to_remove.is_empty());
        assert!(!assignment.sanitize_after_removal);
        assert!(assignment.clear_computed_properties);
    }

    #[test]
    fn remove_hs_assignment_collects_basic_single_neighbor_hydrogen() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = remove_hs_assignment(read(&molecule), &RemoveHsParams::default(), false).unwrap();

        assert_eq!(assignment.atoms_to_remove, vec![hydrogen]);
        assert!(!assignment.sanitize_after_removal);
    }

    #[test]
    fn remove_hs_assignment_collects_degree_zero_hydrogen_when_enabled() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_degree_zero: true,
            ..RemoveHsParams::default()
        };

        let assignment = remove_hs_assignment(read(&molecule), &params, false).unwrap();

        assert_eq!(assignment.atoms_to_remove, vec![hydrogen]);
        assert!(assignment.atom_explicit_hydrogen_updates.is_empty());
    }

    #[test]
    fn remove_hs_assignment_collects_multi_neighbor_hydrogen_when_enabled() {
        let mut builder = MoleculeBuilder::new();
        let left = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let right = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        builder
            .add_bond(BondSpec::new(left, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(hydrogen, right, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_higher_degrees: true,
            ..RemoveHsParams::default()
        };

        let assignment = remove_hs_assignment(read(&molecule), &params, false).unwrap();

        assert_eq!(assignment.atoms_to_remove, vec![hydrogen]);
        assert!(
            assignment
                .atom_explicit_hydrogen_updates
                .contains(&AtomExplicitHydrogenUpdate {
                    atom: left,
                    explicit_hydrogens: 1,
                })
        );
        assert!(
            assignment
                .atom_explicit_hydrogen_updates
                .contains(&AtomExplicitHydrogenUpdate {
                    atom: right,
                    explicit_hydrogens: 1,
                })
        );
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
        assert!(!should(&molecule, isolated_hydrogen, &RemoveHsParams::default()));
        let mut params = RemoveHsParams::default();
        params.remove_degree_zero = true;
        assert!(should(&molecule, isolated_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let query_hydrogen =
            builder.add_atom(AtomSpec::new(Element::H).with_query(QueryNode::predicate(AtomQueryPredicate::Any)));
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
            .add_bond(BondSpec::new(left_hydrogen, right_hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        assert!(!should(&molecule, left_hydrogen, &RemoveHsParams::default()));
        let mut params = RemoveHsParams::default();
        params.remove_only_h_neighbors = true;
        assert!(should(&molecule, left_hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single).with_direction(BondDirection::BeginWedge))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_with_wedged_bond: false,
            ..RemoveHsParams::default()
        };
        assert!(!should(&molecule, hydrogen, &params));

        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::SquarePlanar));
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
    fn remove_hs_assignment_records_chiral_neighbor_inversion() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, oxygen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, nitrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, carbon, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = remove_hs_assignment(read(&molecule), &RemoveHsParams::default(), false).unwrap();

        assert_eq!(assignment.atoms_to_remove, vec![hydrogen]);
        assert_eq!(assignment.atom_chirality_inversions, vec![center]);
    }

    #[test]
    fn remove_hs_assignment_clears_neighbor_double_bond_stereo_for_degree_two_heavy_atom() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(hydrogen, nitrogen, BondOrder::Single))
            .unwrap();
        let stereo_bond = builder
            .add_bond(
                BondSpec::new(nitrogen, carbon, BondOrder::Double)
                    .with_stereo_atoms(hydrogen, carbon)
                    .with_stereo(BondStereo::Cis),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_defining_bond_stereo: true,
            ..RemoveHsParams::default()
        };

        let assignment = remove_hs_assignment(read(&molecule), &params, false).unwrap();

        assert_eq!(
            assignment.bond_stereo_updates,
            vec![BondStereoUpdate {
                bond: stereo_bond,
                stereo: BondStereo::None,
                stereo_atoms: None,
            }]
        );
    }

    #[test]
    fn remove_hs_assignment_moves_end_direction_from_removed_hydrogen_to_single_neighbor() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single).with_direction(BondDirection::EndUpRight))
            .unwrap();
        let neighbor_bond = builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = remove_hs_assignment(read(&molecule), &RemoveHsParams::default(), false).unwrap();

        assert_eq!(
            assignment.bond_direction_updates,
            vec![BondDirectionUpdate {
                bond: neighbor_bond,
                direction: BondDirection::EndDownRight,
            }]
        );
    }

    #[test]
    fn remove_hs_assignment_adjusts_double_bond_stereo_atom_references() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, oxygen, BondOrder::Single))
            .unwrap();
        let double_bond = builder
            .add_bond(
                BondSpec::new(center, carbon, BondOrder::Double)
                    .with_stereo_atoms(hydrogen, carbon)
                    .with_stereo(BondStereo::Cis),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_defining_bond_stereo: true,
            ..RemoveHsParams::default()
        };

        let assignment = remove_hs_assignment(read(&molecule), &params, false).unwrap();

        assert_eq!(
            assignment.bond_stereo_atom_replacements,
            vec![BondStereoAtomReplacement {
                bond: double_bond,
                old_atom: hydrogen,
                new_atom: oxygen,
            }]
        );
        assert_eq!(
            assignment.bond_stereo_updates,
            vec![BondStereoUpdate {
                bond: double_bond,
                stereo: BondStereo::Trans,
                stereo_atoms: Some([hydrogen, carbon]),
            }]
        );
    }

    #[test]
    fn should_remove_h_matches_sgroup_membership_and_role_guards() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom).with_atoms(vec![hydrogen]),
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

    #[test]
    fn remove_hs_assignment_tracks_isotopic_hydrogens_and_clears_stale_property() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_tracked_isotopic_hydrogens(vec![9]));
        let deuterium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        builder
            .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_and_track_isotopes: true,
            ..RemoveHsParams::default()
        };

        let assignment = remove_hs_assignment(read(&molecule), &params, false).unwrap();

        assert_eq!(assignment.atoms_to_remove, vec![deuterium]);
        assert!(
            assignment
                .atom_property_updates
                .contains(&AtomPropertyUpdate::ClearIsotopicHydrogens { atom: carbon })
        );
        assert!(
            assignment
                .atom_property_updates
                .contains(&AtomPropertyUpdate::SetIsotopicHydrogens {
                    atom: carbon,
                    isotopes: vec![2],
                })
        );
    }

    #[test]
    fn remove_hs_assignment_tracks_isotopes_after_nonisotopic_prepass() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let protium = builder.add_atom(AtomSpec::new(Element::H));
        let deuterium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        builder
            .add_bond(BondSpec::new(carbon, protium, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = RemoveHsParams {
            remove_and_track_isotopes: true,
            ..RemoveHsParams::default()
        };

        let assignment = remove_hs_assignment(read(&molecule), &params, false).unwrap();

        assert_eq!(assignment.atoms_to_remove, vec![protium, deuterium]);
        assert!(
            assignment
                .atom_property_updates
                .contains(&AtomPropertyUpdate::SetIsotopicHydrogens {
                    atom: carbon,
                    isotopes: vec![2],
                })
        );
        assert_eq!(
            assignment.atom_explicit_hydrogen_updates.last(),
            Some(&AtomExplicitHydrogenUpdate {
                atom: carbon,
                explicit_hydrogens: 2,
            })
        );
    }

    #[test]
    fn add_hs_assignment_replays_tracked_isotopes_and_clears_tracking_property() {
        let mut builder = MoleculeBuilder::new();
        let nitrogen = builder.add_atom(
            AtomSpec::new(Element::N)
                .with_explicit_hydrogens(2)
                .with_tracked_isotopic_hydrogens(vec![2, 3]),
        );
        let molecule = builder.build().unwrap();

        let assignment = add_hs_assignment(read(&molecule), &AddHsParams::default()).unwrap();

        assert_eq!(assignment.clear_isotopic_hydrogen_properties, vec![nitrogen]);
        assert_eq!(
            assignment.hydrogens_to_add[..2]
                .iter()
                .map(|hydrogen| hydrogen.isotope)
                .collect::<Vec<_>>(),
            vec![Some(2), Some(3)]
        );
    }

    #[test]
    fn add_hs_assignment_copies_residue_info_for_new_hydrogens() {
        let mut builder = MoleculeBuilder::new();
        let nitrogen = builder.add_atom(
            AtomSpec::new(Element::N)
                .with_explicit_hydrogens(1)
                .with_pdb_residue_info(AtomPdbResidueInfo::new(" N  ", 12, "GLY", 3, "A", false)),
        );
        let molecule = builder.build().unwrap();
        let params = AddHsParams {
            explicit_only: true,
            add_residue_info: true,
            ..AddHsParams::default()
        };

        let assignment = add_hs_assignment(read(&molecule), &params).unwrap();

        assert_eq!(assignment.hydrogens_to_add.len(), 1);
        let hydrogen = &assignment.hydrogens_to_add[0];
        assert_eq!(hydrogen.heavy_atom, nitrogen);
        let info = hydrogen.pdb_residue_info.as_ref().unwrap();
        assert_eq!(info.serial_number(), 12);
        assert_eq!(info.residue_name(), "GLY");
        assert_eq!(info.residue_number(), 3);
        assert_eq!(info.chain_id(), "A");
        assert_eq!(info.atom_name(), " H1 ");
    }

    #[test]
    fn add_hs_assignment_records_residue_info_for_existing_hydrogens() {
        let mut builder = MoleculeBuilder::new();
        let nitrogen = builder.add_atom(
            AtomSpec::new(Element::N).with_pdb_residue_info(AtomPdbResidueInfo::new(" N  ", 12, "GLY", 3, "A", false)),
        );
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(crate::BondSpec::new(nitrogen, hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = AddHsParams {
            explicit_only: true,
            add_residue_info: true,
            ..AddHsParams::default()
        };

        let assignment = add_hs_assignment(read(&molecule), &params).unwrap();

        assert!(assignment.hydrogens_to_add.is_empty());
        assert_eq!(assignment.atom_pdb_residue_info_updates.len(), 1);
        let update = &assignment.atom_pdb_residue_info_updates[0];
        assert_eq!(update.atom, hydrogen);
        assert_eq!(update.pdb_residue_info.serial_number(), 12);
        assert_eq!(update.pdb_residue_info.atom_name(), " H1 ");
        assert_eq!(update.pdb_residue_info.residue_name(), "GLY");
        assert_eq!(update.pdb_residue_info.residue_number(), 3);
        assert_eq!(update.pdb_residue_info.chain_id(), "A");
    }

    #[test]
    fn add_hs_assignment_honors_only_on_atoms() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let molecule = builder.build().unwrap();
        let params = AddHsParams {
            only_on_atoms: Some(vec![nitrogen]),
            ..AddHsParams::default()
        };

        let assignment = add_hs_assignment(read(&molecule), &params).unwrap();

        assert!(
            assignment
                .hydrogens_to_add
                .iter()
                .all(|hydrogen| hydrogen.heavy_atom == nitrogen)
        );
        assert_eq!(
            assignment
                .hydrogens_to_add
                .iter()
                .filter(|hydrogen| hydrogen.heavy_atom == carbon)
                .count(),
            0
        );
    }

    #[test]
    fn add_hs_assignment_materializes_implicit_hydrogens_on_degree_zero_hydrogens_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        let first_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let second_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let molecule = builder.build().unwrap();

        let assignment = add_hs_assignment(read(&molecule), &AddHsParams::default()).unwrap();

        assert_eq!(assignment.hydrogens_to_add.len(), 4);
        assert_eq!(
            assignment
                .hydrogens_to_add
                .iter()
                .map(|hydrogen| hydrogen.heavy_atom)
                .collect::<Vec<_>>(),
            vec![oxygen, oxygen, first_hydrogen, second_hydrogen]
        );
        assert!(assignment.hydrogens_to_add.iter().all(|hydrogen| hydrogen.is_implicit));
    }

    #[test]
    fn add_hs_assignment_materializes_explicit_hydrogen_on_hydrogen_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_explicit_hydrogens(1));
        let molecule = builder.build().unwrap();

        let assignment = add_hs_assignment(read(&molecule), &AddHsParams::default()).unwrap();

        assert_eq!(assignment.hydrogens_to_add.len(), 1);
        assert_eq!(assignment.hydrogens_to_add[0].heavy_atom, hydrogen);
        assert!(!assignment.hydrogens_to_add[0].is_implicit);
        assert_eq!(
            assignment.atom_explicit_hydrogen_updates,
            vec![AtomExplicitHydrogenUpdate {
                atom: hydrogen,
                explicit_hydrogens: 0,
            }]
        );
    }

    #[test]
    fn add_hs_assignment_allows_bondless_heavy_atom_without_explicit_hydrogen_atoms() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let molecule = builder.build().unwrap();

        let assignment = add_hs_assignment(read(&molecule), &AddHsParams::default()).unwrap();

        assert_eq!(assignment.hydrogens_to_add.len(), 4);
        assert!(
            assignment
                .hydrogens_to_add
                .iter()
                .all(|hydrogen| hydrogen.heavy_atom == carbon)
        );
    }

    #[test]
    fn remove_hs_assignment_accumulates_explicit_h_count_for_multiple_removed_hydrogens() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let first_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let second_hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, first_hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(carbon, second_hydrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = remove_hs_assignment(read(&molecule), &RemoveHsParams::default(), false).unwrap();

        assert_eq!(
            assignment.atom_explicit_hydrogen_updates,
            vec![AtomExplicitHydrogenUpdate {
                atom: carbon,
                explicit_hydrogens: 2,
            }]
        );
    }

    #[test]
    fn remove_hs_assignment_allows_import_facing_wedge_directions_to_fall_through() {
        for direction in [
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::EitherDouble,
        ] {
            let mut builder = MoleculeBuilder::new();
            let carbon = builder.add_atom(AtomSpec::new(Element::C));
            let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
            builder
                .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single).with_direction(direction))
                .unwrap();
            let molecule = builder.build().unwrap();

            let assignment = remove_hs_assignment(read(&molecule), &RemoveHsParams::default(), false).unwrap();

            assert_eq!(assignment.atoms_to_remove, vec![hydrogen]);
        }
    }

    fn should(molecule: &Molecule, atom: AtomId, params: &RemoveHsParams) -> bool {
        let context = RemoveHsPredicateContext::new(read(molecule));
        should_remove_h(&context, atom, params).unwrap()
    }
}
