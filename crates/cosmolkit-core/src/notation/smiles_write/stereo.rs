use super::*;

// RDKit✔️✔️: std::string getAtomChiralityInfo(const Atom *atom) {
// RDKit✔️✔️:   auto allowNontet = Chirality::getAllowNontetrahedralChirality();
// RDKit✔️✔️:   std::string atString;
// RDKit✔️✔️:   switch (atom->getChiralTag()) {
// RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CW: atString = "@@"; break;
// RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CCW: atString = "@"; break;
// RDKit✔️✔️:     default: break;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (atString.empty() && allowNontet) {
// RDKit✔️✔️:     switch (atom->getChiralTag()) {
// RDKit✔️✔️:       case Atom::CHI_SQUAREPLANAR: atString = "@SP"; break;
// RDKit✔️✔️:       case Atom::CHI_TRIGONALBIPYRAMIDAL: atString = "@TB"; break;
// RDKit✔️✔️:       case Atom::CHI_OCTAHEDRAL: atString = "@OH"; break;
// RDKit✔️✔️:       default: break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (!atString.empty()) {
// RDKit✔️✔️:       int permutation = 0;
// RDKit✔️✔️:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER &&
// RDKit✔️✔️:           atom->getPropIfPresent(common_properties::_chiralPermutation,
// RDKit✔️✔️:                                  permutation) &&
// RDKit✔️✔️:           !SmilesParseOps::checkChiralPermutation(atom->getChiralTag(),
// RDKit✔️✔️:                                                   permutation)) {
// RDKit✔️✔️:         throw ValueErrorException("bad chirality spec");
// RDKit✔️✔️:       } else if (permutation) {
// RDKit✔️✔️:         atString += std::to_string(permutation);
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return atString;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION getAtomChiralityInfo
pub(super) fn get_atom_chirality_info_with_inversion(
    molecule: &Molecule,
    atom: AtomId,
    chiral_tag_override: Option<ChiralTag>,
    invert: bool,
    permutation_override: Option<u32>,
) -> Result<String, SmilesWriteError> {
    let atom = &molecule.atoms()[atom.index()];
    let base_chiral_tag = chiral_tag_override.unwrap_or_else(|| atom.chiral_tag());
    let chiral_tag = if invert {
        match base_chiral_tag {
            ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
            ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
            other => other,
        }
    } else {
        base_chiral_tag
    };
    match chiral_tag {
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
            let perm = permutation_override
                .or_else(|| atom.chiral_permutation())
                .unwrap_or(0);
            let res = if perm % 2 == 0 {
                match chiral_tag {
                    ChiralTag::TetrahedralCw => "@@",
                    _ => "@",
                }
            } else {
                match chiral_tag {
                    ChiralTag::TetrahedralCw => "@",
                    _ => "@@",
                }
            };
            Ok(res.to_string())
        }
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
            let mut res = match chiral_tag {
                ChiralTag::SquarePlanar => "@SP".to_string(),
                ChiralTag::TrigonalBipyramidal => "@TB".to_string(),
                ChiralTag::Octahedral => "@OH".to_string(),
                _ => unreachable!(),
            };
            if let Some(permutation) = permutation_override.or_else(|| atom.chiral_permutation()) {
                validate_writer_chiral_permutation(chiral_tag, permutation)?;
                if permutation != 0 {
                    res.push_str(&permutation.to_string());
                }
            }
            Ok(res)
        }
        _ => Ok(String::new()),
    }
}

pub(super) fn validate_writer_chiral_permutation(
    chiral_tag: ChiralTag,
    permutation: u32,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesParseOps::checkChiralPermutation
    // RDKit✔️✔️: bool checkChiralPermutation(int chiralTag, int permutation) {
    // RDKit✔️✔️:   if (chiralTag > RDKit::Atom::ChiralType::CHI_OTHER &&
    // RDKit✔️✔️:       permutationLimits.find(chiralTag) != permutationLimits.end() &&
    // RDKit✔️✔️:       (permutation < 0 || permutation > permutationLimits.at(chiralTag))) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SmilesParseOps::checkChiralPermutation
    let limit = match chiral_tag {
        ChiralTag::Allene => Some(2),
        ChiralTag::SquarePlanar => Some(3),
        ChiralTag::TrigonalBipyramidal => Some(20),
        ChiralTag::Octahedral => Some(30),
        _ => None,
    };
    if let Some(limit) = limit
        && permutation > limit
    {
        return Err(SmilesWriteError::InvalidChiralPermutation {
            chiral_tag,
            permutation,
            limit,
        });
    }
    Ok(())
}

pub(super) fn atom_needs_bracket(
    molecule: &Molecule,
    atom: AtomId,
    at_string: &str,
    params: &SmilesWriteParams,
) -> Result<bool, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION atomNeedsBracket
    // RDKit✔️✔️: bool atomNeedsBracket(const Atom *atom, const std::string &atString,
    // RDKit✔️✔️:                       const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(atom, "null atom");
    // RDKit✔️✔️:   auto num = atom->getAtomicNum();
    // RDKit✔️✔️:   if (!inOrganicSubset(num)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atom->getFormalCharge()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.doIsomericSmiles && (atom->getIsotope() || !atString.empty())) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom->hasProp(common_properties::molAtomMapNumber)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit✔️✔️:   int totalValence = atom->getTotalValence();
    // RDKit✔️✔️:   bool nonStandard = false;
    // RDKit✔️✔️:   if (atom->getNumRadicalElectrons()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit✔️✔️:              atom->getNumExplicitHs()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nonStandard) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check for bonds to a metal
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION atomNeedsBracket
    let atom_id = atom;
    let atom = &molecule.atoms()[atom.index()];
    if !in_organic_subset(atom.atomic_number())? {
        return Ok(true);
    }
    if atom.formal_charge() != 0 || atom.atom_map().is_some() {
        return Ok(true);
    }
    if params.do_isomeric_smiles && (atom.isotope().is_some() || !at_string.is_empty()) {
        return Ok(true);
    }
    // RDKit✔️✔️: if (atom->getNumRadicalElectrons()) { nonStandard = true; }
    if atom.radical_electrons() != 0 {
        return Ok(true);
    }
    // RDKit✔️✔️: else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit✔️✔️:            atom->getNumExplicitHs()) {
    // RDKit✔️✔️:   nonStandard = true;
    // RDKit✔️✔️: }
    if matches!(atom.atomic_number(), 7 | 15) && atom.is_aromatic() && atom.explicit_hydrogens() > 0
    {
        return Ok(true);
    }
    // RDKit✔️✔️:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit✔️✔️:   int totalValence = atom->getTotalValence();
    // RDKit✔️✔️:   bool nonStandard = false;
    // RDKit✔️✔️:   if (atom->getNumRadicalElectrons()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit✔️✔️:              atom->getNumExplicitHs()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit✔️✔️:   }
    if let Ok(Some(valence_list)) = crate::valence::rdkit_valence_list(atom.atomic_number()) {
        if let Some(&default_valence) = valence_list.first() {
            let total_valence = total_valence_for_writer(molecule, atom_id).unwrap_or_else(|| {
                molecule
                    .bonds()
                    .iter()
                    .fold(i32::from(atom.explicit_hydrogens()), |acc, bond| {
                        if bond.begin() == atom_id || bond.end() == atom_id {
                            acc + crate::valence::bond_valence_contrib(bond, atom_id)
                                .unwrap_or(0.0)
                                .round() as i32
                        } else {
                            acc
                        }
                    })
            });
            let total_hs = total_num_hydrogens_for_writer(molecule, atom_id) as i32;
            if total_valence != default_valence && total_hs > 0 {
                return Ok(true);
            }
        }
    }
    if molecule.bonds().iter().any(|bond| {
        bond_other_atom(bond, atom_id)
            .map(|other| rdkit_query_ops_is_metal(molecule.atoms()[other.index()].atomic_number()))
            .unwrap_or(false)
    }) {
        return Ok(true);
    }
    Ok(false)
}

pub(super) fn rdkit_query_ops_is_metal(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    // RDKit✔️✔️: // !#0!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️: bool isMetal(const Atom &atom) { ... makeMAtomQuery()->Match(&atom); }
    // END RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    !matches!(
        atomic_number,
        0 | 1
            | 2
            | 5
            | 6
            | 7
            | 8
            | 9
            | 10
            | 14
            | 15
            | 16
            | 17
            | 18
            | 33
            | 34
            | 35
            | 36
            | 52
            | 53
            | 54
            | 85
            | 86
    )
}

pub(super) fn update_property_cache_for_smiles(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles writer property-cache prep
    // RDKit✔️✔️: for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:   atom->updatePropertyCache(false);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles writer property-cache prep
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles ring-info lifetime
    // RDKit✔️✔️: // writer prep here updates atom property cache only; it does not
    // RDKit✔️✔️: // proactively recompute ring info. Traversal/ranking code consults the
    // RDKit✔️✔️: // molecule's existing RingInfo and only computes fast rings on demand in
    // RDKit✔️✔️: // specific canonical-ranking paths.
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles ring-info lifetime
    *molecule = molecule.with_assigned_valence()?;
    Ok(())
}

pub(super) fn clear_fragment_temp_molecule_computed_stereo_props_for_writer(
    molecule: &mut Molecule,
) {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags fragment temporary molecule path
    // RDKit❗✔️: if (comp.size() == 1 || ...) {
    // RDKit❗✔️:   SubsetOptions opts{.sanitize = sanitizeFrags,
    // RDKit❗✔️:                      .clearComputedProps = true,
    // RDKit❗✔️:                      .copyCoordinates = copyConformers,
    // RDKit❗✔️:                      .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
    // RDKit❗✔️:   auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   frag->beginBatchEdit();
    // RDKit❗✔️:   ...
    // RDKit❗✔️:   frag->commitBatchEdit();
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags fragment temporary molecule path
    // BEGIN RDKIT CPP FUNCTION RWMol::removeAtom / commitBatchEdit computed prop reset
    // RDKit❗✔️: // clear computed properties and reset our ring info structure
    // RDKit❗✔️: // they are pretty likely to be wrong now:
    // RDKit❗✔️: if (clearProps) {
    // RDKit❗✔️:   clearComputedProps(true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RWMol::removeAtom / commitBatchEdit computed prop reset
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗✔️: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit❗✔️:   if (includeRings) {
    // RDKit❗✔️:     this->dp_ringInfo->reset();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDProps::clearComputedProps();
    // RDKit❗✔️:   for (auto atom : atoms()) {
    // RDKit❗✔️:     atom->clearComputedProps();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto bond : bonds()) {
    // RDKit❗✔️:     bond->clearComputedProps();
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    let fragment_count = crate::notation::fragment::get_fragment_atom_mapping(molecule)
        .iter()
        .copied()
        .max()
        .map_or(0, |max_fragment| max_fragment + 1);
    if fragment_count <= 1 {
        return;
    }

    molecule.properties_mut().clear_prop("_StereochemDone");
    for atom in &mut molecule.topology_block_mut().atoms {
        atom.clear_prop("_CIPCode");
        atom.clear_prop("_CIPRank");
        atom.clear_prop("_ChiralityPossible");
        atom.clear_prop("_ringStereochemCand");
        atom.clear_prop("_ringStereoAtoms");
    }
    for bond in &mut molecule.topology_block_mut().bonds {
        bond.clear_prop("_CIPCode");
    }
}

pub(super) fn assign_stereochemistry_for_smiles(
    molecule: &mut Molecule,
    clean_stereo: bool,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::assignStereochemistry dispatcher tail
    // RDKit❗✔️: void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
    // RDKit❗✔️:                            bool flagPossibleStereoCenters) {
    // RDKit❗✔️:   if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ...
    // RDKit❗✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Chirality::assignStereochemistry dispatcher tail
    // RDKit calls assignStereochemistry() on the temporary writer molecule
    // unless stereochemistry was already prepared. COSMolKit's current stereo
    // assignment reads typed state and validates/detects writer-visible stereo;
    // future ports can make this mutating without changing the writer contract.
    // BEGIN RDKIT CPP FUNCTION Chirality::stereoPerception cleanIt prelude
    // RDKit❗✔️: void stereoPerception(ROMol &mol, bool cleanIt,
    // RDKit❗✔️:                       bool flagPossibleStereoCenters) {
    // RDKit❗✔️:   if (cleanIt) {
    // RDKit❗✔️:     for (auto atom : mol.atoms()) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_CIPCode);
    // RDKit❗✔️:       atom->clearProp(common_properties::_ChiralityPossible);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto bond : mol.bonds()) {
    // RDKit❗✔️:       bond->clearProp(common_properties::_CIPCode);
    // RDKit❗✔️:       if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit❗✔️:         bond->setStereo(Bond::BondStereo::STEREOANY);
    // RDKit❗✔️:         bond->getStereoAtoms().clear();
    // RDKit❗✔️:         bond->setBondDir(Bond::BondDir::NONE);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // we need cis/trans markers on the double bonds... set those now:
    // RDKit❗✔️:   MolOps::setBondStereoFromDirections(mol);
    // RDKit❗✔️:   ...
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Chirality::stereoPerception cleanIt prelude
    //
    // This ports the writer-visible `cleanIt` property reset and EITHERDOUBLE
    // normalization exactly. The downstream `findPotentialStereo()` /
    // `updateDoubleBondStereo()` pipeline is still only partially reproduced by
    // the current typed-state helpers below, so the behavior marker remains `❗`.
    if clean_stereo {
        for atom in &mut molecule.topology_block_mut().atoms {
            atom.clear_prop("_CIPCode");
            atom.clear_prop("_ChiralityPossible");
        }
        for bond in &mut molecule.topology_block_mut().bonds {
            bond.clear_prop("_CIPCode");
            if bond.direction() == BondDirection::EitherDouble {
                bond.set_stereo(BondStereo::Any);
                bond.set_stereo_atoms(None);
                bond.set_direction(BondDirection::None);
            }
        }
    }
    ensure_fast_rings_for_writer_stereo_perception(molecule)?;
    crate::stereo::assign_stereochemistry(molecule)?;
    assign_double_bond_stereo_for_writer_working_copy(molecule)?;
    let ranks = crate::stereo::assign_atom_cip_ranks(molecule)?;
    let (_, atom_labels, _) = crate::stereo::assign_atom_chiral_codes(molecule, &ranks)?;
    for (atom_idx, cip_code) in atom_labels {
        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(atom_idx) {
            atom_mut.set_prop("_CIPCode", cip_code);
        }
    }
    for (atom_idx, rank) in ranks.iter().copied().enumerate() {
        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(atom_idx) {
            atom_mut.set_prop("_CIPRank", rank.to_string());
        }
    }
    if clean_stereo {
        apply_clean_stereo_ring_special_cases_for_writer(molecule, &ranks)?;
    }
    molecule.properties_mut().set_prop("_StereochemDone", "1");
    Ok(())
}

pub(super) fn ensure_fast_rings_for_writer_stereo_perception(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::legacyStereoPerception ring prelude
    // RDKit✔️✔️: void legacyStereoPerception(ROMol &mol, bool cleanIt,
    // RDKit✔️✔️:                             bool flagPossibleStereoCenters) {
    // RDKit✔️✔️:   mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // later we're going to need ring information, get it now if we don't
    // RDKit✔️✔️:   // have it already:
    // RDKit✔️✔️:   // NOTE, if called from the SMART code, the ring info will be DUMMY, and
    // RDKit✔️✔️:   // contains no information
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::legacyStereoPerception ring prelude
    if molecule
        .derived_cache()
        .rings
        .as_ref()
        .is_some_and(crate::RingInfo::is_find_fast_or_better)
    {
        return Ok(());
    }
    let rings = crate::fast_find_rings(molecule)?;
    molecule.derived_cache_mut().rings = Some(rings);
    Ok(())
}

pub(super) fn apply_clean_stereo_ring_special_cases_for_writer(
    molecule: &mut Molecule,
    ranks: &[u32],
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception cleanIt ring special cases
    // RDKit❗✔️: if (cleanIt) {
    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->hasProp(common_properties::_ringStereochemCand)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_ringStereochemCand);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (atom->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_ringStereoAtoms);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    // RDKit❗✔️:   Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases, atomRanks);
    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:         !Chirality::hasNonTetrahedralStereo(atom) &&
    // RDKit❗✔️:         !atom->hasProp(common_properties::_CIPCode) &&
    // RDKit❗✔️:         (!possibleSpecialCases[atom->getIdx()] ||
    // RDKit❗✔️:          !atom->hasProp(common_properties::_ringStereoAtoms))) {
    // RDKit❗✔️:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:       ...
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION legacyStereoPerception cleanIt ring special cases
    for atom in &mut molecule.topology_block_mut().atoms {
        atom.clear_prop("_ringStereochemCand");
        atom.clear_prop("_ringStereoAtoms");
    }

    let special_cases = crate::stereo::find_chiral_atom_special_cases(molecule, ranks)?;
    let special_case_atoms = special_cases
        .iter()
        .map(|case| case.atom_idx)
        .collect::<BTreeSet<_>>();
    for case in &special_cases {
        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(case.atom_idx) {
            atom_mut.set_prop("_ringStereochemCand", "1");
            atom_mut.set_prop(
                "_ringStereoAtoms",
                serialize_ring_stereo_atoms(&case.ring_stereo_atoms),
            );
        }
    }

    let atom_ids = molecule
        .atoms()
        .iter()
        .map(|atom| atom.id())
        .collect::<Vec<_>>();
    for atom_id in atom_ids {
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.chiral_tag() == ChiralTag::Unspecified
            || crate::stereo::has_non_tetrahedral_stereo(atom)
            || atom.prop("_CIPCode").is_some()
            || (special_case_atoms.contains(&atom_id.index())
                && atom.prop("_ringStereoAtoms").is_some())
        {
            continue;
        }

        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(atom_id.index()) {
            atom_mut.set_chiral_tag(ChiralTag::Unspecified);
            if atom_mut.explicit_hydrogens() == 1
                && atom_mut.formal_charge() == 0
                && !atom_mut.is_aromatic()
            {
                atom_mut.set_explicit_hydrogens(0);
                atom_mut.set_no_implicit(false);
            }
        }
    }
    *molecule = molecule.with_assigned_valence()?;
    Ok(())
}

pub(crate) fn serialize_ring_stereo_atoms(ring_stereo_atoms: &[(bool, usize)]) -> String {
    ring_stereo_atoms
        .iter()
        .map(|(same_orientation, atom_idx)| {
            let signed = if *same_orientation {
                *atom_idx as i32 + 1
            } else {
                -(*atom_idx as i32 + 1)
            };
            signed.to_string()
        })
        .collect::<Vec<_>>()
        .join(",")
}

pub(super) fn parse_ring_stereo_atoms_prop(encoded: &str) -> Option<Vec<(bool, usize)>> {
    let mut result = Vec::new();
    for token in encoded.split(',').filter(|token| !token.is_empty()) {
        let entry = token.parse::<i32>().ok()?;
        if entry == 0 {
            return None;
        }
        let atom_idx = entry.unsigned_abs() as usize - 1;
        result.push((entry > 0, atom_idx));
    }
    Some(result)
}

pub(super) fn writer_ring_stereo_atoms(
    atom: &crate::Atom,
    atom_id: AtomId,
) -> Result<Option<Vec<(bool, usize)>>, SmilesWriteError> {
    let Some(encoded) = atom.prop("_ringStereoAtoms") else {
        if atom.prop("_ringStereochemCand").is_some() {
            return Err(SmilesWriteError::InvalidRingStereoState {
                atom: atom_id.index(),
                requirement: "`_ringStereochemCand` requires `_ringStereoAtoms`",
            });
        }
        return Ok(None);
    };
    parse_ring_stereo_atoms_prop(encoded)
        .ok_or(SmilesWriteError::InvalidRingStereoState {
            atom: atom_id.index(),
            requirement: "`_ringStereoAtoms` must be a valid encoded ring-neighbor list",
        })
        .map(Some)
}

pub(super) fn assign_double_bond_stereo_for_writer_working_copy(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::assignStereochemistry double-bond section
    // RDKit✔️✔️: assignBondStereoCodes(mol, atomRanks);
    // RDKit✔️✔️: if (bond->getStereo() == Bond::STEREONONE) {
    // RDKit✔️✔️:   ... find directed neighbor bonds around double bond ...
    // RDKit✔️✔️:   bond->setStereoAtoms(beginControl, endControl);
    // RDKit✔️✔️:   bond->setStereo(Bond::STEREOE or Bond::STEREOZ);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::assignStereochemistry double-bond section
    let ranks = crate::stereo::assign_atom_cip_ranks(molecule)?;
    let (_, assignments, changed) = crate::stereo::assign_bond_stereo_codes(molecule, &ranks);
    if !changed {
        return Ok(());
    }
    for (bond_idx, stereo, begin_control, end_control) in assignments {
        let stereo = match stereo {
            // RDKit✔️✔️: bond->setStereo(Bond::STEREOE);
            crate::stereo::DoubleBondStereo::E => BondStereo::E,
            // RDKit✔️✔️: bond->setStereo(Bond::STEREOZ);
            crate::stereo::DoubleBondStereo::Z => BondStereo::Z,
            crate::stereo::DoubleBondStereo::Unknown => BondStereo::Any,
        };
        let bond = &mut molecule.topology_block_mut().bonds[bond_idx];
        if bond.stereo() == BondStereo::None {
            bond.set_stereo_atoms(Some([AtomId::new(begin_control), AtomId::new(end_control)]));
            bond.set_stereo(stereo);
        }
    }
    molecule
        .derived_cache_mut()
        .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
    Ok(())
}

pub(super) fn canonicalize_enhanced_stereo_for_smiles(
    _molecule: &Molecule,
) -> Result<(), SmilesWriteError> {
    // Enhanced stereo group canonicalization is only needed for CX SMILES.
    // For plain SMILES, stereo groups are already in typed state.
    Ok(())
}

pub(super) fn cleanup_stereo_groups_for_cx_smiles(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // RDKit cleanupStereoGroups() also applies atropisomer-specific cleanup
    // before CX serialization, moving atropisomer participation from atom
    // members onto bond members in each stereo group.
    let current_groups = molecule.stereo_groups().to_vec();
    if current_groups.is_empty() {
        return Ok(());
    }

    let mut cleaned_groups = Vec::with_capacity(current_groups.len());
    for group in current_groups {
        let mut kept_atoms = Vec::with_capacity(group.atoms().len());
        let mut kept_bonds = group.bonds().to_vec();
        for atom in group.atoms().iter().copied() {
            let mut found_atrop = false;
            for bond_id in incident_bonds(molecule, atom) {
                let bond = &molecule.bonds()[bond_id.index()];
                if matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw) {
                    found_atrop = true;
                    if !kept_bonds.contains(&bond_id) {
                        kept_bonds.push(bond_id);
                    }
                }
            }
            if !found_atrop {
                kept_atoms.push(atom);
            }
        }
        let mut cleaned_group = crate::StereoGroup::new(group.kind(), kept_atoms, kept_bonds);
        if let Some(id) = group.id() {
            cleaned_group = cleaned_group.with_id(id);
        }
        cleaned_groups.push(cleaned_group);
    }

    molecule.topology_block_mut().stereo_groups = cleaned_groups;
    Ok(())
}

// RDKit✔️✔️: void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
// RDKit✔️✔️:               unsigned int maxBackTracks) {
// RDKit✔️✔️:   boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
// RDKit✔️✔️:   atomsToUse.set();
// RDKit✔️✔️:   boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
// RDKit✔️✔️:   bondsToUse.set();
// RDKit✔️✔️:   details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
// RDKit✔️✔️:                             canonical, maxBackTracks);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolOps::Kekulize
/// Run operation-routed kekulization on the molecule and return a new
/// `Molecule` with resolved bond orders and cleared aromatic flags.
///
/// The writer's `do_kekule=true` path is source-aligned to RDKit's
/// `Kekulize(mol, true, true, 100)` behavior: it uses canonical fragment
/// ranking when resolving aromatic systems before SMILES traversal.
pub(super) fn kekulize_for_smiles(molecule: &Molecule) -> Result<Molecule, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION details::KekulizeFragment
    // RDKit✔️✔️:     VECT_INT_VECT allringsSSSR;
    // RDKit✔️✔️:     if (!mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:       MolOps::findSSSR(mol, allringsSSSR);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     const VECT_INT_VECT &allrings =
    // RDKit✔️✔️:         allringsSSSR.empty() ? mol.getRingInfo()->atomRings() : allringsSSSR;
    // END RDKIT CPP FUNCTION details::KekulizeFragment
    let assignment = crate::kekulize::kekulize_assignment(
        molecule,
        molecule.derived_cache().rings.as_ref(),
        true,
        true,
        100,
    )?;
    let mut kekulized = molecule.clone();
    crate::kekulize::apply_kekulize_assignment(kekulized.topology_block_mut(), &assignment);
    // RDKit✔️✔️:         if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 15) &&
    // RDKit✔️✔️:             atom->getFormalCharge() == 0 && atom->getNumExplicitHs() == 1) {
    // RDKit✔️✔️:           atom->setNoImplicit(false);
    // RDKit✔️✔️:           atom->setNumExplicitHs(0);
    // RDKit✔️✔️:           atom->updatePropertyCache(false);
    // RDKit✔️✔️:         }
    //
    // COSMolKit stores RDKit's property-cache-derived implicit hydrogen state
    // in `derived_cache.valence` instead of per-atom mutable cache fields.
    // After applying the kekulize assignment, refresh that cache so
    // `GetAtomSmiles()`-equivalent total-H reads see the post-kekulization
    // atom state.
    kekulized.derived_cache_mut().valence = Some(crate::valence::assign_valence_with_options(
        &kekulized,
        crate::ValenceModel::RdkitLike,
        false,
    )?);
    Ok(kekulized)
}

pub(super) fn normalize_dative_bonds_for_plain_smiles(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // RDKit✔️✔️:     if (doingCXSmiles || !params.includeDativeBonds) {
    // RDKit✔️✔️:       for (auto bond : tmol->bonds()) {
    // RDKit✔️✔️:         if (bond->getBondType() == Bond::DATIVE) {
    // RDKit✔️✔️:           // we are intentionally only handling DATIVE here. The other weird
    // RDKit✔️✔️:           // RDKit dative alternatives really shouldn't ever show up.
    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:           // update the explicit valence of the begin atom since the implicit
    // RDKit✔️✔️:           // valence will no longer be properly perceived
    // RDKit✔️✔️:           bond->getBeginAtom()->calcExplicitValence(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    for bond in &mut molecule.topology_block_mut().bonds {
        if bond.order() == crate::BondOrder::Dative {
            bond.set_order(crate::BondOrder::Single);
        }
    }
    // COSMolKit's writer currently stores RDKit property-cache-derived valence
    // state in `derived_cache.valence`. Recomputing that cache from the
    // topology-mutated working copy changes the observable writer behavior for
    // former dative donors (`N->O` becomes `[NH3]O` instead of RDKit's `NO`).
    //
    // RDKit updates only the begin atom's internal explicit-valence cache here;
    // it does not rebuild a whole-molecule valence assignment for the writer
    // working copy. Preserving the existing derived valence cache matches the
    // current RDKit writer state boundary for the modeled input space.
    Ok(())
}

pub(super) fn normalize_dative_bonds_for_cx_smiles(
    _molecule: &Molecule,
) -> Result<(), SmilesWriteError> {
    // In CX SMILES mode, dative bonds are preserved and written as `_Z:2:...`
    // entries in the CX extension section. No molecule mutation needed.
    Ok(())
}

pub(super) fn normalize_hydrogen_bonds_for_cx_smiles(
    _molecule: &Molecule,
) -> Result<(), SmilesWriteError> {
    // In CX SMILES mode, hydrogen bonds are preserved and written as `_Z:1:...`
    // entries in the CX extension section. No molecule mutation needed.
    Ok(())
}

pub(super) fn apply_cx_bond_direction_policy(
    molecule: &mut Molecule,
    restore_bond_dirs: RestoreBondDirOption,
) -> Result<(), SmilesWriteError> {
    // In CX SMILES mode, bond directions are kept for stereo output.
    // `RestoreBondDirOption::Clear` clears unknown/either directions
    // for plain SMILES, but in CX mode we preserve them since the
    // bond direction state is already stored in typed atom/bond state.
    // For `RestoreBondDirOption::True`, RDKit restores from the molblock;
    // this is a no-op here since our directions are already in typed state.
    match restore_bond_dirs {
        RestoreBondDirOption::Clear => {
            let mut changed = false;
            for bond in &mut molecule.topology_block_mut().bonds {
                if matches!(
                    bond.direction(),
                    BondDirection::Unknown | BondDirection::EitherDouble
                ) {
                    bond.set_direction(BondDirection::None);
                    changed = true;
                }
            }
            if changed {
                molecule
                    .derived_cache_mut()
                    .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
            }
        }
        RestoreBondDirOption::True | RestoreBondDirOption::None => {
            // Keep existing direction state as-is.
        }
    }
    Ok(())
}

pub(super) fn remove_plain_smiles_only_cx_state(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    let mut changed = false;
    for bond in &mut molecule.topology_block_mut().bonds {
        if matches!(
            bond.direction(),
            BondDirection::Unknown | BondDirection::EitherDouble
        ) {
            bond.set_direction(BondDirection::None);
            changed = true;
        }
        if bond.stereo() == BondStereo::Any {
            bond.set_stereo_atoms(None);
            bond.set_stereo(BondStereo::None);
            changed = true;
        }
    }
    if changed {
        molecule
            .derived_cache_mut()
            .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
    }
    Ok(())
}

pub(super) fn validate_cx_extension_plan(fields: CxSmilesFields) -> Result<(), SmilesWriteError> {
    // All CX field types are now supported through write_cx_* functions.
    // If a specific field is requested, the corresponding writer will
    // produce output or return empty string if no data exists.
    let _ = fields;
    Ok(())
}

pub(super) fn is_minimal_plain_smiles_path(params: &SmilesWriteParams) -> bool {
    !params.do_isomeric_smiles
        && !params.do_kekule
        && !params.canonical
        && !params.clean_stereo
        && !params.all_hydrogens_explicit
        && !params.do_random
        && params.include_dative_bonds
        && !params.ignore_atom_map_numbers
}

pub(super) fn validate_minimal_plain_smiles_molecule(molecule: &Molecule) -> bool {
    for atom in molecule.atoms() {
        if atom.query().is_some()
            || atom.radical_electrons() != 0
            || atom.chiral_tag() != ChiralTag::Unspecified
        {
            // [deferred] Minimal plain SMILES path doesn't handle query atoms,
            // radical-bearing atoms, or chiral atoms. Falls through to the
            // standard path which supports these features via get_atom_smiles()
            // and isomeric SMILES / stereochemistry handling.
            return false;
        }
        if atom.props().keys().any(|key| {
            key != "dummyLabel" && key != "_SmilesStart" && key != "_supplementalSmilesLabel"
        }) {
            // [deferred] Non-whitelisted atom properties (e.g. map numbers,
            // custom data) are not supported by the minimal fast path.
            // Falls through to the standard path which writes atoms with
            // their full property set.
            return false;
        }
    }
    for bond in molecule.bonds() {
        if bond.direction() != BondDirection::None
            || bond.stereo() != BondStereo::None
            || bond.query().is_some()
            || !matches!(
                bond.order(),
                BondOrder::Single
                    | BondOrder::Double
                    | BondOrder::Triple
                    | BondOrder::Quadruple
                    | BondOrder::Aromatic
                    | BondOrder::Dative
            )
        {
            // [deferred] Bonds with direction, stereo, query, or non-standard
            // orders (e.g. Hydrogen, Zero) are not supported by the minimal
            // fast path. Falls through to the standard path.
            return false;
        }
    }
    true
}

pub(super) fn canonical_dfs_traversal(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    do_isomeric_smiles: bool,
    clean_stereo: bool,
    do_random: bool,
    bond_symbols: Option<&[String]>,
) -> Result<CanonicalTraversalResult, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalDFSTraversal
    // RDKit✔️✔️: void canonicalDFSTraversal(ROMol &mol, int atomIdx, int inBondIdx,
    // RDKit✔️✔️:                            std::vector<AtomColors> &colors,
    // RDKit✔️✔️:                            const UINT_VECT &ranks, MolStack &molStack,
    // RDKit✔️✔️:                            VECT_INT_VECT &atomRingClosures,
    // RDKit✔️✔️:                            std::vector<INT_LIST> &atomTraversalBondOrder,
    // RDKit✔️✔️:                            const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit✔️✔️:                            const std::vector<std::string> *bondSymbols,
    // RDKit✔️✔️:                            bool doRandom) {
    // RDKit✔️✔️:   dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomRingClosures,
    // RDKit✔️✔️:                 bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:   boost::dynamic_bitset<> cyclesAvailable(MAX_CYCLES);
    // RDKit✔️✔️:   cyclesAvailable.set();
    // RDKit✔️✔️:   dfsBuildStack(mol, atomIdx, inBondIdx, colors, ranks, cyclesAvailable,
    // RDKit✔️✔️:                 molStack, atomRingClosures, atomTraversalBondOrder, bondsInPlay,
    // RDKit✔️✔️:                 bondSymbols, doRandom);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::canonicalDFSTraversal
    //
    // BEGIN RDKIT CPP FUNCTION Canon::dfsFindCycles / dfsBuildStack traversal sections
    // RDKit✔️✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit✔️✔️:   std::vector<PossibleType> possibles;
    // RDKit✔️✔️:   std::sort(possibles.begin(), possibles.end(), _possibleCompare);
    // RDKit✔️✔️:   for (auto &possible : possibles) {
    // RDKit✔️✔️:     int possibleIdx = std::get<1>(possible);
    // RDKit✔️✔️:     Bond *bond = std::get<2>(possible);
    // RDKit✔️✔️:     switch (colors[possibleIdx]) {
    // RDKit✔️✔️:       case WHITE_NODE:
    // RDKit✔️✔️:         dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit✔️✔️:                       atomRingClosures, bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case GREY_NODE:
    // RDKit✔️✔️:         atomRingClosures[possibleIdx].push_back(bond->getIdx());
    // RDKit✔️✔️:         atomRingClosures[atomIdx].push_back(bond->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   colors[atomIdx] = BLACK_NODE;
    // RDKit✔️✔️:   molStack.push_back(MolStackElem(atom));
    // RDKit✔️✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit✔️✔️:   if (!atomRingClosures[atomIdx].empty()) {
    // RDKit✔️✔️:     for (auto bIdx : atomRingClosures[atomIdx]) {
    // RDKit✔️✔️:       Bond *bond = mol.getBondWithIdx(bIdx);
    // RDKit✔️✔️:       if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                              ringIdx)) {
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(ringIdx));
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         bond->setProp(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                       static_cast<unsigned int>(lowestRingIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(lowestRingIdx));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
    // RDKit✔️✔️:        ++possiblesIt) {
    // RDKit✔️✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit✔️✔️:       molStack.push_back(MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:     dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit✔️✔️:                   cyclesAvailable, molStack, atomRingClosures,
    // RDKit✔️✔️:                   atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit✔️✔️:       molStack.push_back(MolStackElem(")", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION Canon::dfsFindCycles / dfsBuildStack traversal sections
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Color {
        White,
        Grey,
        Black,
    }

    fn sorted_incident_bonds(
        molecule: &Molecule,
        atom: AtomId,
        rank_by_atom: &[usize],
        colors: &[Color],
        ring_info: &crate::RingInfo,
        do_random: bool,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) -> Vec<(BondId, AtomId)> {
        let mut incident = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .filter(|neighbor| bonds_in_play[neighbor.bond.index()])
            .map(|neighbor| (neighbor.bond, AtomId::new(neighbor.atom_index)))
            .map(|(bond, other)| {
                let rank = if do_random {
                    next_random_smiles_u64() as i64
                } else {
                    traversal_possible_rank(
                        molecule,
                        bond,
                        other,
                        rank_by_atom,
                        colors,
                        ring_info,
                        bond_symbols,
                    )
                };
                (rank, bond, other)
            })
            .collect::<Vec<_>>();
        incident.sort_by_key(|(rank, _, _)| *rank);
        incident
            .into_iter()
            .map(|(_, bond, other)| (bond, other))
            .collect()
    }

    fn traversal_possible_rank(
        molecule: &Molecule,
        bond: BondId,
        other: AtomId,
        rank_by_atom: &[usize],
        colors: &[Color],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
    ) -> i64 {
        let mut rank = rank_by_atom[other.index()] as i64;
        let bond_order_rank = rdkit_bond_order_rank(molecule.bonds()[bond.index()].order());
        if colors[other.index()] == Color::Grey {
            rank -= (CANON_MAX_BONDTYPE + 1) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
            if let Some(symbols) = bond_symbols {
                rank += i64::from(gboost_hash_range(symbols[bond.index()].as_bytes()) % 5000)
                    * CANON_MAX_NATOMS;
            } else {
                rank += (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS;
            }
        } else if ring_info.num_bond_rings(bond) > 0 {
            if let Some(symbols) = bond_symbols {
                rank += i64::from(gboost_hash_range(symbols[bond.index()].as_bytes()) % 5000)
                    * CANON_MAX_NATOMS
                    * CANON_MAX_NATOMS;
            } else {
                rank +=
                    (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
            }
        }
        rank
    }

    fn gboost_hash_range(bytes: &[u8]) -> u32 {
        let mut seed = 0u32;
        for byte in bytes {
            seed ^= u32::from(*byte)
                .wrapping_add(0x9e37_79b9)
                .wrapping_add(seed << 6)
                .wrapping_add(seed >> 2);
        }
        seed
    }

    fn rdkit_bond_order_rank(order: BondOrder) -> i64 {
        match order {
            BondOrder::Null | BondOrder::Unspecified => 0,
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::OneAndHalf => 7,
            BondOrder::TwoAndHalf => 8,
            BondOrder::ThreeAndHalf => 9,
            BondOrder::FourAndHalf => 10,
            BondOrder::FiveAndHalf => 11,
            BondOrder::Aromatic => 12,
            BondOrder::Ionic => 13,
            BondOrder::Hydrogen => 14,
            BondOrder::ThreeCenter => 15,
            BondOrder::DativeOne => 16,
            BondOrder::Dative => 17,
            BondOrder::DativeLeft => 18,
            BondOrder::DativeRight => 19,
            BondOrder::Other => 20,
            BondOrder::Zero => 21,
        }
    }

    fn dfs_find_cycles(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        parent_bond: Option<BondId>,
        do_random: bool,
        colors: &mut [Color],
        atom_ring_closures: &mut [Vec<BondId>],
        rank_by_atom: &[usize],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) {
        colors[atom.index()] = Color::Grey;
        for (bond, other) in sorted_incident_bonds(
            molecule,
            atom,
            rank_by_atom,
            colors,
            ring_info,
            do_random,
            bond_symbols,
            bonds_in_play,
        ) {
            if Some(bond) == parent_bond {
                continue;
            }
            match colors[other.index()] {
                Color::White => {
                    dfs_find_cycles(
                        molecule,
                        plan,
                        other,
                        Some(bond),
                        do_random,
                        colors,
                        atom_ring_closures,
                        rank_by_atom,
                        ring_info,
                        bond_symbols,
                        bonds_in_play,
                    );
                }
                Color::Grey => {
                    atom_ring_closures[other.index()].push(bond);
                    atom_ring_closures[atom.index()].push(bond);
                }
                Color::Black => {}
            }
        }
        colors[atom.index()] = Color::Black;
    }

    fn dfs_build_stack(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        parent_bond: Option<BondId>,
        do_random: bool,
        colors: &mut [Color],
        atom_ring_closures: &[Vec<BondId>],
        cycles_available: &mut [bool],
        traversal_ring_indices: &mut [Option<usize>],
        traversal_ring_closure_bonds: &mut [bool],
        rank_by_atom: &[usize],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        stack: &mut Vec<MolStackElem>,
        atom_traversal_bond_order: &mut [Vec<BondId>],
        bonds_in_play: &[bool],
    ) -> Result<(), SmilesWriteError> {
        let mut seen_from_here = vec![false; molecule.num_atoms()];
        seen_from_here[atom.index()] = true;
        stack.push(MolStackElem::Atom(atom));
        colors[atom.index()] = Color::Grey;
        let mut trav_list = Vec::new();
        if let Some(parent_bond) = parent_bond {
            trav_list.push(parent_bond);
        }

        let mut rings_closed = Vec::new();
        for &bond in &atom_ring_closures[atom.index()] {
            trav_list.push(bond);
            if let Some(other) = bond_other_atom(&molecule.bonds()[bond.index()], atom) {
                seen_from_here[other.index()] = true;
            }
            if let Some(ring_idx) = traversal_ring_indices[bond.index()] {
                stack.push(MolStackElem::Bond(bond, atom));
                stack.push(MolStackElem::Ring { bond, ring_idx });
                rings_closed.push(ring_idx - 1);
            } else {
                let Some(lowest_ring_idx) =
                    cycles_available.iter().position(|available| *available)
                else {
                    return invariant_stage_error(
                        SmilesPlanStage::ShortTermBondWriter,
                        "write_ring_closure() could not allocate a free ring index",
                    );
                };
                cycles_available[lowest_ring_idx] = false;
                let ring_idx = lowest_ring_idx + 1;
                traversal_ring_indices[bond.index()] = Some(ring_idx);
                traversal_ring_closure_bonds[bond.index()] = true;
                stack.push(MolStackElem::Ring { bond, ring_idx });
            }
        }
        for ring_idx in rings_closed {
            cycles_available[ring_idx] = true;
        }

        let children = sorted_incident_bonds(
            molecule,
            atom,
            rank_by_atom,
            colors,
            ring_info,
            do_random,
            bond_symbols,
            bonds_in_play,
        )
        .into_iter()
        .filter(|(bond, other)| {
            Some(*bond) != parent_bond
                && colors[other.index()] == Color::White
                && !seen_from_here[other.index()]
        })
        .collect::<Vec<_>>();

        for (idx, (bond, child)) in children.iter().enumerate() {
            let is_branch = idx + 1 != children.len();
            if is_branch {
                stack.push(MolStackElem::BranchOpen);
            }
            stack.push(MolStackElem::Bond(*bond, atom));
            trav_list.push(*bond);
            dfs_build_stack(
                molecule,
                plan,
                *child,
                Some(*bond),
                do_random,
                colors,
                atom_ring_closures,
                cycles_available,
                traversal_ring_indices,
                traversal_ring_closure_bonds,
                rank_by_atom,
                ring_info,
                bond_symbols,
                stack,
                atom_traversal_bond_order,
                bonds_in_play,
            )?;
            if is_branch {
                stack.push(MolStackElem::BranchClose);
            }
        }
        atom_traversal_bond_order[atom.index()] = trav_list;
        colors[atom.index()] = Color::Black;
        Ok(())
    }

    let mut cycle_colors = vec![Color::White; molecule.num_atoms()];
    let mut atom_ring_closures = vec![Vec::new(); molecule.num_atoms()];
    let mut stack = Vec::new();
    let mut atom_traversal_bond_order = vec![Vec::new(); molecule.num_atoms()];
    let mut rank_by_atom = vec![usize::MAX; molecule.num_atoms()];
    let mut bonds_in_play = vec![false; molecule.num_bonds()];
    let mut traversal_ring_closure_bonds = vec![false; molecule.num_bonds()];
    for (idx, atom) in plan.atoms.iter().enumerate() {
        rank_by_atom[atom.index()] = ranks[idx];
    }
    for bond in &plan.bonds {
        bonds_in_play[bond.index()] = true;
    }
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalDFSTraversal ring-info requirement
    // RDKit❗✔️: traversal only queries `numBondRings()` / `numAtomRings()` from
    // RDKit❗✔️: a `fastFindRings()`-or-better RingInfo. It does not require
    // RDKit❗✔️: symmetrized SSSR membership here.
    // END RDKIT CPP FUNCTION Canon::canonicalDFSTraversal ring-info requirement
    let ring_info = molecule
        .derived_cache()
        .rings
        .clone()
        .filter(|rings| rings.is_find_fast_or_better())
        .unwrap_or(crate::fast_find_rings(molecule)?);

    dfs_find_cycles(
        molecule,
        plan,
        start_atom,
        None,
        do_random,
        &mut cycle_colors,
        &mut atom_ring_closures,
        &rank_by_atom,
        &ring_info,
        bond_symbols,
        &bonds_in_play,
    );
    let mut stack_colors = vec![Color::White; molecule.num_atoms()];
    let mut cycles_available = vec![true; molecule.num_bonds().max(1)];
    let mut traversal_ring_indices = vec![None; molecule.num_bonds()];
    dfs_build_stack(
        molecule,
        plan,
        start_atom,
        None,
        do_random,
        &mut stack_colors,
        &atom_ring_closures,
        &mut cycles_available,
        &mut traversal_ring_indices,
        &mut traversal_ring_closure_bonds,
        &rank_by_atom,
        &ring_info,
        bond_symbols,
        &mut stack,
        &mut atom_traversal_bond_order,
        &bonds_in_play,
    )?;
    let chiral_adjustments = if do_isomeric_smiles {
        compute_writer_chiral_adjustments(
            molecule,
            plan,
            start_atom,
            &atom_ring_closures,
            &atom_traversal_bond_order,
            &stack,
            clean_stereo,
        )?
    } else {
        WriterChiralAdjustments::default()
    };

    if plan
        .atoms
        .iter()
        .any(|atom_id| stack_colors[atom_id.index()] != Color::Black)
    {
        return invariant_stage_error(
            SmilesPlanStage::ShortTermAtomWriter,
            "canonical_dfs_traversal() did not visit every atom in the fragment plan",
        );
    }

    Ok(CanonicalTraversalResult {
        stack,
        traversal_ring_closure_bonds,
        chiral_tag_overrides: chiral_adjustments.chiral_tag_overrides,
        chiral_inversions: chiral_adjustments.chiral_inversions,
        chiral_permutations: chiral_adjustments.chiral_permutations,
        broken_chiral_atoms: chiral_adjustments.broken_chiral_atoms,
    })
}

#[cfg(test)]
pub(super) fn debug_atom_ring_closures_for_writer(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    bond_symbols: Option<&[String]>,
) -> Result<(crate::RingInfo, Vec<Vec<BondId>>), SmilesWriteError> {
    fn gboost_hash_range_debug(bytes: &[u8]) -> u32 {
        let mut seed = 0u32;
        for byte in bytes {
            seed ^= u32::from(*byte)
                .wrapping_add(0x9e37_79b9)
                .wrapping_add(seed << 6)
                .wrapping_add(seed >> 2);
        }
        seed
    }

    fn rdkit_bond_order_rank_debug(order: BondOrder) -> i64 {
        match order {
            BondOrder::Null | BondOrder::Unspecified => 0,
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::OneAndHalf => 7,
            BondOrder::TwoAndHalf => 8,
            BondOrder::ThreeAndHalf => 9,
            BondOrder::FourAndHalf => 10,
            BondOrder::FiveAndHalf => 11,
            BondOrder::Aromatic => 12,
            BondOrder::Ionic => 13,
            BondOrder::Hydrogen => 14,
            BondOrder::ThreeCenter => 15,
            BondOrder::DativeOne => 16,
            BondOrder::Dative => 17,
            BondOrder::DativeLeft => 18,
            BondOrder::DativeRight => 19,
            BondOrder::Other => 20,
            BondOrder::Zero => 21,
        }
    }

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum DebugColor {
        White,
        Grey,
        Black,
    }

    fn sorted_incident_bonds_debug(
        molecule: &Molecule,
        atom: AtomId,
        rank_by_atom: &[usize],
        colors: &[DebugColor],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) -> Vec<(BondId, AtomId)> {
        let mut incident = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .filter(|neighbor| bonds_in_play[neighbor.bond.index()])
            .map(|neighbor| (neighbor.bond, AtomId::new(neighbor.atom_index)))
            .map(|(bond, other)| {
                let mut rank = rank_by_atom[other.index()] as i64;
                let bond_order_rank =
                    rdkit_bond_order_rank_debug(molecule.bonds()[bond.index()].order());
                if colors[other.index()] == DebugColor::Grey {
                    rank -= (CANON_MAX_BONDTYPE + 1) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
                    if let Some(symbols) = bond_symbols {
                        rank += i64::from(
                            gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000,
                        ) * CANON_MAX_NATOMS;
                    } else {
                        rank += (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS;
                    }
                } else if ring_info.num_bond_rings(bond) > 0 {
                    if let Some(symbols) = bond_symbols {
                        rank += i64::from(
                            gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000,
                        ) * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                    } else {
                        rank += (CANON_MAX_BONDTYPE - bond_order_rank)
                            * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                    }
                }
                (rank, bond, other)
            })
            .collect::<Vec<_>>();
        incident.sort_by_key(|(rank, _, _)| *rank);
        incident
            .into_iter()
            .map(|(_, bond, other)| (bond, other))
            .collect()
    }

    fn dfs_find_cycles_debug(
        molecule: &Molecule,
        atom: AtomId,
        parent_bond: Option<BondId>,
        colors: &mut [DebugColor],
        atom_ring_closures: &mut [Vec<BondId>],
        rank_by_atom: &[usize],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) {
        colors[atom.index()] = DebugColor::Grey;
        let possibles = sorted_incident_bonds_debug(
            molecule,
            atom,
            rank_by_atom,
            colors,
            ring_info,
            bond_symbols,
            bonds_in_play,
        );
        eprint!("debug_dfsFindCycles atom={} possibles", atom.index());
        for (bond, other) in &possibles {
            let mut rank = rank_by_atom[other.index()] as i64;
            let bond_order_rank =
                rdkit_bond_order_rank_debug(molecule.bonds()[bond.index()].order());
            if colors[other.index()] == DebugColor::Grey {
                rank -= (CANON_MAX_BONDTYPE + 1) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
                if let Some(symbols) = bond_symbols {
                    rank +=
                        i64::from(gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000)
                            * CANON_MAX_NATOMS;
                } else {
                    rank += (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS;
                }
            } else if ring_info.num_bond_rings(*bond) > 0 {
                if let Some(symbols) = bond_symbols {
                    rank +=
                        i64::from(gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000)
                            * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                } else {
                    rank += (CANON_MAX_BONDTYPE - bond_order_rank)
                        * CANON_MAX_NATOMS
                        * CANON_MAX_NATOMS;
                }
            }
            eprint!(" ({},r={},b={})", other.index(), rank, bond.index());
        }
        eprintln!();
        for (bond, other) in possibles {
            if Some(bond) == parent_bond {
                continue;
            }
            match colors[other.index()] {
                DebugColor::White => dfs_find_cycles_debug(
                    molecule,
                    other,
                    Some(bond),
                    colors,
                    atom_ring_closures,
                    rank_by_atom,
                    ring_info,
                    bond_symbols,
                    bonds_in_play,
                ),
                DebugColor::Grey => {
                    atom_ring_closures[other.index()].push(bond);
                    atom_ring_closures[atom.index()].push(bond);
                }
                DebugColor::Black => {}
            }
        }
        colors[atom.index()] = DebugColor::Black;
    }

    let mut rank_by_atom = vec![usize::MAX; molecule.num_atoms()];
    for (idx, atom) in plan.atoms.iter().enumerate() {
        rank_by_atom[atom.index()] = ranks[idx];
    }
    let bonds_in_play = {
        let mut value = vec![false; molecule.num_bonds()];
        for bond in &plan.bonds {
            value[bond.index()] = true;
        }
        value
    };
    let ring_info = crate::fast_find_rings(molecule)?;
    let mut colors = vec![DebugColor::White; molecule.num_atoms()];
    let mut atom_ring_closures = vec![Vec::new(); molecule.num_atoms()];
    dfs_find_cycles_debug(
        molecule,
        start_atom,
        None,
        &mut colors,
        &mut atom_ring_closures,
        &rank_by_atom,
        &ring_info,
        bond_symbols,
        &bonds_in_play,
    );
    Ok((ring_info, atom_ring_closures))
}

#[derive(Debug, Default)]
pub(super) struct WriterChiralAdjustments {
    pub(super) chiral_tag_overrides: BTreeMap<AtomId, ChiralTag>,
    pub(super) chiral_inversions: BTreeSet<AtomId>,
    pub(super) chiral_permutations: BTreeMap<AtomId, u32>,
    pub(super) broken_chiral_atoms: BTreeSet<AtomId>,
}

pub(super) fn compute_writer_chiral_adjustments(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    atom_ring_closures: &[Vec<BondId>],
    atom_traversal_bond_order: &[Vec<BondId>],
    stack: &[MolStackElem],
    _clean_stereo: bool,
) -> Result<WriterChiralAdjustments, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral traversal section
    // RDKit✔️✔️: std::vector<INT_LIST> atomTraversalBondOrder(nAtoms);
    // RDKit✔️✔️: std::vector<int> atomPermutationIndices(nAtoms, 0);
    // RDKit✔️✔️: const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];
    // RDKit✔️✔️: if (trueOrder.size() < atom->getDegree()) { ... append missing bonds ... }
    // RDKit✔️✔️: if (!perm) { nSwaps = atom->getPerturbationOrder(tOrder); }
    // RDKit✔️✔️: else { insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit✔️✔️:        perm = Chirality::getChiralPermutation(atom, tOrder); }
    // RDKit✔️✔️: if (doChiralInversions && chiralAtomNeedsTagInversion(...)) { ++nSwaps; }
    // RDKit✔️✔️: if (nSwaps % 2) { numSwapsChiralAtoms.set(atom->getIdx()); }
    // RDKit✔️✔️: atomPermutationIndices[atom->getIdx()] = perm;
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral traversal section
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    // RDKit✔️✔️: std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
    // RDKit✔️✔️: for (const auto &msI : molStack) { if (msI.type == MOL_STACK_ATOM) atomVisitOrders[...] = pos; ++pos; }
    // RDKit✔️✔️: boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    // RDKit✔️✔️: if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) { ... setChiralTag(CCW) ... }
    // RDKit✔️✔️: else if (msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) && mol.getStereoGroups().size() > sgidx) { ... }
    // RDKit✔️✔️: else { if (numSwapsChiralAtoms[idx]) atom->invertChirality(); else if (atomPermutationIndices[idx]) atom->setProp(_chiralPermutation, ...); }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    let plan_atoms = plan.atoms.iter().copied().collect::<BTreeSet<_>>();
    let plan_bonds = plan.bonds.iter().copied().collect::<BTreeSet<_>>();
    let mut adjustments = WriterChiralAdjustments::default();
    let mut num_swaps_chiral_atoms = vec![false; molecule.num_atoms()];
    let mut atom_permutation_indices = vec![0_u32; molecule.num_atoms()];

    for atom_id in &plan.atoms {
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.chiral_tag() == ChiralTag::Unspecified {
            continue;
        }

        let incident = incident_bonds(molecule, *atom_id);
        if incident.iter().any(|bond| !plan_bonds.contains(bond)) {
            adjustments.broken_chiral_atoms.insert(*atom_id);
            continue;
        }

        let mut true_order = atom_traversal_bond_order[atom_id.index()].clone();
        if true_order.len() < incident.len() {
            for bond in &incident {
                if !true_order.contains(bond) {
                    true_order.push(*bond);
                }
            }
        }

        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                let mut swaps = count_swaps_to_interconvert(&true_order, &incident)?;
                if chiral_atom_needs_tag_inversion_for_writer(
                    molecule,
                    *atom_id,
                    *atom_id == start_atom,
                    atom_ring_closures[atom_id.index()].len(),
                ) {
                    swaps = swaps.saturating_add(1);
                }
                if swaps % 2 == 1 {
                    num_swaps_chiral_atoms[atom_id.index()] = true;
                }
            }
            ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
                let mut probe = true_order.into_iter().map(Some).collect::<Vec<_>>();
                insert_implicit_nontetrahedral_neighbors_for_writer(
                    &mut probe,
                    atom.chiral_tag(),
                    *atom_id == start_atom,
                );
                let permutation =
                    crate::notation::smiles::nontetrahedral_chiral_permutation_for_probe(
                        molecule, *atom_id, &probe, false,
                    )
                    .map_err(|_| {
                        crate::UnsupportedFeatureError::from_spec(&crate::SMILES_WRITE_FEATURE)
                    })?;
                if permutation != 0 {
                    atom_permutation_indices[atom_id.index()] = permutation;
                }
            }
            _ => {}
        }
    }

    let stereo_group_by_atom = molecule
        .stereo_groups()
        .iter()
        .enumerate()
        .flat_map(|(group_idx, group)| {
            group
                .atoms()
                .iter()
                .copied()
                .map(move |atom_id| (atom_id, group_idx))
        })
        .collect::<BTreeMap<_, _>>();
    let mut atom_visit_orders = vec![usize::MAX; molecule.num_atoms()];
    for (pos, element) in stack.iter().enumerate() {
        if let MolStackElem::Atom(atom_id) = element {
            atom_visit_orders[atom_id.index()] = pos;
        }
    }
    let mut ring_stereo_chem_adjusted = vec![false; molecule.num_atoms()];
    for element in stack {
        let MolStackElem::Atom(atom_id) = element else {
            continue;
        };
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.chiral_tag() == ChiralTag::Unspecified
            || adjustments.broken_chiral_atoms.contains(atom_id)
        {
            continue;
        }

        let ring_stereo_atoms = writer_ring_stereo_atoms(atom, *atom_id)?;
        if let Some(ring_stereo_atoms) = ring_stereo_atoms {
            if !ring_stereo_chem_adjusted[atom_id.index()] {
                adjustments
                    .chiral_tag_overrides
                    .insert(*atom_id, ChiralTag::TetrahedralCcw);
                ring_stereo_chem_adjusted[atom_id.index()] = true;
            }
            let source_tag = adjustments
                .chiral_tag_overrides
                .get(atom_id)
                .copied()
                .unwrap_or_else(|| atom.chiral_tag());
            for (same_orientation, neighbor_idx) in ring_stereo_atoms {
                let neighbor = AtomId::new(neighbor_idx);
                if !plan_atoms.contains(&neighbor)
                    || ring_stereo_chem_adjusted[neighbor_idx]
                    || atom_visit_orders[neighbor_idx] <= atom_visit_orders[atom_id.index()]
                {
                    continue;
                }
                let mut neighbor_tag = source_tag;
                if !same_orientation {
                    neighbor_tag = invert_tetrahedral_chiral_tag(neighbor_tag);
                }
                if num_swaps_chiral_atoms[atom_id.index()] {
                    if !num_swaps_chiral_atoms[neighbor_idx] {
                        neighbor_tag = invert_tetrahedral_chiral_tag(neighbor_tag);
                    }
                } else if num_swaps_chiral_atoms[neighbor_idx] {
                    neighbor_tag = invert_tetrahedral_chiral_tag(neighbor_tag);
                }
                adjustments
                    .chiral_tag_overrides
                    .insert(neighbor, neighbor_tag);
                ring_stereo_chem_adjusted[neighbor_idx] = true;
            }
        } else if let Some(group_idx) = stereo_group_by_atom.get(atom_id).copied() {
            if let Some(group) = molecule.stereo_groups().get(group_idx) {
                let swap_it = matches!(
                    adjustments
                        .chiral_tag_overrides
                        .get(atom_id)
                        .copied()
                        .unwrap_or_else(|| atom.chiral_tag()),
                    ChiralTag::TetrahedralCw
                );
                if swap_it {
                    let current = adjustments
                        .chiral_tag_overrides
                        .get(atom_id)
                        .copied()
                        .unwrap_or_else(|| atom.chiral_tag());
                    adjustments
                        .chiral_tag_overrides
                        .insert(*atom_id, invert_tetrahedral_chiral_tag(current));
                }
                if swap_it || num_swaps_chiral_atoms[atom_id.index()] {
                    for other_atom in group.atoms() {
                        if *other_atom == *atom_id || !plan_atoms.contains(other_atom) {
                            continue;
                        }
                        let current = adjustments
                            .chiral_tag_overrides
                            .get(other_atom)
                            .copied()
                            .unwrap_or_else(|| molecule.atoms()[other_atom.index()].chiral_tag());
                        adjustments
                            .chiral_tag_overrides
                            .insert(*other_atom, invert_tetrahedral_chiral_tag(current));
                    }
                }
            }
        } else {
            match atom.chiral_tag() {
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                    if num_swaps_chiral_atoms[atom_id.index()] {
                        adjustments.chiral_inversions.insert(*atom_id);
                    }
                }
                ChiralTag::SquarePlanar
                | ChiralTag::TrigonalBipyramidal
                | ChiralTag::Octahedral => {
                    let permutation = atom_permutation_indices[atom_id.index()];
                    if permutation != 0 {
                        adjustments
                            .chiral_permutations
                            .insert(*atom_id, permutation);
                    }
                }
                _ => {}
            }
        }
    }
    Ok(adjustments)
}

pub(super) fn invert_tetrahedral_chiral_tag(chiral_tag: ChiralTag) -> ChiralTag {
    match chiral_tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        other => other,
    }
}

pub(super) fn insert_implicit_nontetrahedral_neighbors_for_writer(
    bonds: &mut Vec<Option<BondId>>,
    chiral_tag: ChiralTag,
    first_atom: bool,
) {
    // BEGIN RDKIT CPP FUNCTION Canon::insertImplicitNbors
    // RDKit✔️✔️: static void insertImplicitNbors(INT_LIST &bonds, const Atom::ChiralType tag,
    // RDKit✔️✔️:                                 const bool firstAtom) {
    // RDKit✔️✔️:   unsigned int ref_max = Chirality::getMaxNbors(tag);
    // RDKit✔️✔️:   if (bonds.size() < ref_max) {
    // RDKit✔️✔️:     if (firstAtom) {
    // RDKit✔️✔️:       bonds.insert(bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bonds.insert(++bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::insertImplicitNbors
    let Some(ref_max) = crate::notation::smiles::nontetrahedral_max_neighbors(chiral_tag) else {
        return;
    };
    if bonds.len() < ref_max {
        let missing = ref_max - bonds.len();
        let insert_at = if first_atom { 0 } else { 1.min(bonds.len()) };
        bonds.splice(insert_at..insert_at, std::iter::repeat_n(None, missing));
    }
}

pub(super) fn incident_bonds(molecule: &Molecule, atom: AtomId) -> Vec<BondId> {
    // BEGIN RDKIT CPP FUNCTION ROMol::atomBonds
    // RDKit✔️✔️: CXXBondIterator<const MolGraph, Bond *const, MolGraph::out_edge_iterator>
    // RDKit✔️✔️: atomBonds(Atom const *at) const {
    // RDKit✔️✔️:   auto pr = getAtomBonds(at);
    // RDKit✔️✔️:   return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::atomBonds
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| neighbor.bond)
        .collect()
}

pub(super) fn count_swaps_to_interconvert(
    probe: &[BondId],
    reference: &[BondId],
) -> Result<usize, SmilesWriteError> {
    if probe.len() != reference.len() {
        return invariant_stage_error(
            SmilesPlanStage::ShortTermAtomWriter,
            "count_swaps_to_interconvert() probe/reference length mismatch",
        );
    }
    let mut work = probe.to_vec();
    let mut swaps = 0usize;
    for (idx, expected) in reference.iter().copied().enumerate() {
        if work[idx] == expected {
            continue;
        }
        let Some(found_idx) = work[idx..]
            .iter()
            .position(|bond| *bond == expected)
            .map(|offset| idx + offset)
        else {
            return invariant_stage_error(
                SmilesPlanStage::ShortTermAtomWriter,
                "count_swaps_to_interconvert() expected bond missing from probe order",
            );
        };
        work.swap(idx, found_idx);
        swaps = swaps.saturating_add(1);
    }
    Ok(swaps)
}

pub(super) fn chiral_atom_needs_tag_inversion_for_writer(
    molecule: &Molecule,
    atom_id: AtomId,
    is_atom_first: bool,
    num_closures: usize,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
    // RDKit✔️✔️: bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
    // RDKit✔️✔️:                                  const RDKit::Atom *atom, bool isAtomFirst,
    // RDKit✔️✔️:                                  size_t numClosures) {
    // RDKit✔️✔️:   return atom->getDegree() == 3 &&
    // RDKit✔️✔️:          ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
    // RDKit✔️✔️:           (!details::atomHasFourthValence(atom) && numClosures == 1 &&
    // RDKit✔️✔️:            !details::isUnsaturated(atom, mol)));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
    let atom = &molecule.atoms()[atom_id.index()];
    let degree = incident_bonds(molecule, atom_id).len();
    degree == 3
        && ((is_atom_first && atom.explicit_hydrogens() == 1)
            || (!atom_has_fourth_valence_for_writer(molecule, atom_id)
                && num_closures == 1
                && !atom_is_unsaturated_for_writer(molecule, atom_id)))
}

pub(super) fn atom_has_fourth_valence_for_writer(molecule: &Molecule, atom_id: AtomId) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::details::atomHasFourthValence
    // RDKit✔️✔️: bool atomHasFourthValence(const Atom *atom) {
    // RDKit✔️✔️:   if (atom->getNumExplicitHs() == 1 ||
    // RDKit✔️✔️:       (!atom->needsUpdatePropertyCache() &&
    // RDKit✔️✔️:        atom->getValence(Atom::ValenceType::IMPLICIT) == 1)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom->hasQuery()) {
    // RDKit✔️✔️:     return hasSingleHQuery(atom->getQuery());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::details::atomHasFourthValence
    //
    // COSMolKit does not model RDKit's per-atom property-cache dirty bit.
    // The implicit-H branch therefore relies on the current derived valence
    // cache when present, matching the writer's current modeled state space.
    let atom = &molecule.atoms()[atom_id.index()];
    atom.explicit_hydrogens() == 1
        || molecule
            .derived_cache()
            .valence
            .as_ref()
            .and_then(|valence| valence.implicit_hydrogens.get(atom_id.index()))
            .is_some_and(|implicit_hydrogens| *implicit_hydrogens == 1)
        || atom
            .query()
            .is_some_and(atom_query_has_single_h_count_for_writer)
}

pub(super) fn atom_is_unsaturated_for_writer(molecule: &Molecule, atom_id: AtomId) -> bool {
    molecule
        .bonds()
        .iter()
        .filter(|bond| bond.begin() == atom_id || bond.end() == atom_id)
        .any(|bond| bond_order_as_double_for_writer(bond.order()) > 1.0)
}

pub(super) fn bond_order_as_double_for_writer(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single
        | BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::OneAndHalf | BondOrder::Aromatic => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Null
        | BondOrder::Ionic
        | BondOrder::Hydrogen
        | BondOrder::ThreeCenter
        | BondOrder::Other
        | BondOrder::Zero
        | BondOrder::Unspecified => 0.0,
    }
}

pub(super) fn atom_query_has_single_h_count_for_writer(
    query: &QueryNode<AtomQueryPredicate>,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::details::hasSingleHQuery
    // RDKit✔️✔️: bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q) {
    // RDKit✔️✔️:   PRECONDITION(q, "bad query");
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   const auto &descr = q->getDescription();
    // RDKit✔️✔️:   if (descr == "AtomAnd") {
    // RDKit✔️✔️:     for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
    // RDKit✔️✔️:       const auto &cDescr = (*cIt)->getDescription();
    // RDKit✔️✔️:       if (cDescr == "AtomHCount") {
    // RDKit✔️✔️:         return !(*cIt)->getNegation() &&
    // RDKit✔️✔️:                ((ATOM_EQUALS_QUERY *)(*cIt).get())->getVal() == 1;
    // RDKit✔️✔️:       } else if (cDescr == "AtomAnd") {
    // RDKit✔️✔️:         res = hasSingleHQuery((*cIt).get());
    // RDKit✔️✔️:         if (res) {
    // RDKit✔️✔️:           return true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::details::hasSingleHQuery
    match query {
        QueryNode::And(children) => {
            let mut found_nested = false;
            for child in children {
                match child {
                    QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(count))
                        if *count == 1 =>
                    {
                        return true;
                    }
                    QueryNode::And(_) if atom_query_has_single_h_count_for_writer(child) => {
                        found_nested = true;
                        break;
                    }
                    _ => {}
                }
            }
            found_nested
        }
        _ => false,
    }
}
