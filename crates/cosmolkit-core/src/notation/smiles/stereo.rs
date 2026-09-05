use super::*;

pub(crate) fn assign_stereochemistry_cleanup_subset(mol: &mut Molecule, clean_it: bool) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION assignStereochemistry cleanup subset
    // RDKit❗✔️: void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
    // RDKit❗✔️:                            bool flagPossibleStereoCenters) {
    // RDKit❗✔️:   ...
    // RDKit❗✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION assignStereochemistry cleanup subset
    ensure_valence_for_stereo(mol)?;
    assign_double_bond_stereo_after_smiles_parse(mol, clean_it)?;
    mol.properties_mut().set_computed_prop("_StereochemDone", "1");
    Ok(())
}

pub(crate) fn assign_double_bond_stereo_after_smiles_parse(
    mol: &mut Molecule,
    clean_it: bool,
) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION MolFromSmiles final stereochemistry assignment subset
    // RDKit❗❌: if (res && (params.sanitize || params.removeHs)) {
    // RDKit❗❌:   if (res->hasProp(detail::_needsDetectBondStereo)) {
    // RDKit❗❌:     if (conf || conf3d) {
    // RDKit❗❌:       MolOps::clearSingleBondDirFlags(*res);
    // RDKit❗❌:     }
    // RDKit❗❌:     MolOps::setDoubleBondNeighborDirections(*res, conf ? conf : conf3d);
    // RDKit❗❌:   }
    // RDKit❗❌:   res->clearProp(detail::_needsDetectBondStereo);
    // RDKit❗❌:   MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION MolFromSmiles final stereochemistry assignment subset
    // RDKit✔️✔️: mol.clearProp("_needsDetectBondStereo");
    mol.properties_mut().clear_prop("_needsDetectBondStereo");
    // RDKit✔️✔️: if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:   MolOps::fastFindRings(mol);
    // RDKit✔️✔️: }
    ensure_fast_find_rings_for_legacy_stereo(mol)?;
    //
    // `setBondStereoFromDirections()` is not part of RDKit's
    // legacyStereoPerception() parser finalization path for coordinate-free
    // SMILES. Double-bond stereo remains unassigned here until the
    // fixed-point `assignBondStereoCodes()` loop below consumes the raw single-
    // bond directions together with CIP ranks.
    //
    // COSMolKit does not yet port full assignStereochemistry(). For the
    // no-coordinate SMILES parse path, this subset reproduces:
    // 1. cleanIt small-ring double-bond cleanup
    // 2. legacy fixed-point assignAtomChiralCodes()/assignBondStereoCodes()
    //    assignment for direction-marked double bonds and dependent atom CIP labels
    //
    // We intentionally do not assign double-bond stereo directly from raw
    // neighboring directions on every double bond. RDKit gates final
    // assignment through assignBondStereoCodes(), which also checks the
    // substituent distinguishability implied by CIP ranks.
    //
    // BEGIN RDKIT CPP FUNCTION assignStereochemistry cleanIt small-ring double-bond section
    // RDKit❗✔️: if (cleanIt) {
    // RDKit❗✔️:   // enforce no stereo on small rings
    // RDKit❗✔️:   if ((bond->getBondType() == Bond::DOUBLE ||
    // RDKit❗✔️:        bond->getBondType() == Bond::AROMATIC) &&
    // RDKit❗✔️:       !shouldDetectDoubleBondStereo(bond)) {
    // RDKit❗✔️:     if (bond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit❗✔️:       bond->setBondDir(Bond::NONE);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (bond->getStereo() != Bond::STEREONONE) {
    // RDKit❗✔️:       bond->setStereo(Bond::STEREONONE);
    // RDKit❗✔️:       bond->getStereoAtoms().clear();
    // RDKit❗✔️:     }
    // RDKit❗✔️:     continue;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION assignStereochemistry cleanIt small-ring double-bond section
    //
    // RDKit's MolFromSmiles() finalization enters assignStereochemistry()
    // with `cleanIt=true` whenever the `(params.sanitize || params.removeHs)`
    // branch runs. The public SMILES parser defaults `removeHs=true`, so
    // sanitize=false inputs still take this cleanup path.
    if clean_it {
        // This local subset only ports the small-ring double-bond cleanup
        // branch plus the no-coordinate assignBondCisTrans path. Full
        // assignStereochemistry() still remains broader than the currently
        // modeled SMILES parse finalization.
        let bond_ids: Vec<BondId> = mol.bonds().iter().map(|bond| bond.id()).collect();
        for bond_id in bond_ids.iter().copied() {
            let Some(snapshot) = mol.bonds().get(bond_id.index()).cloned() else {
                continue;
            };
            if !matches!(snapshot.order(), BondOrder::Double | BondOrder::Aromatic) {
                continue;
            }
            let should_detect = crate::stereo::should_detect_double_bond_stereo(mol, bond_id)?;
            if let Some(bond_mut) = mol.topology_block_mut().bonds.get_mut(bond_id.index()) {
                if !should_detect {
                    if bond_mut.direction() == BondDirection::EitherDouble {
                        bond_mut.set_direction(BondDirection::None);
                    }
                    if bond_mut.stereo() != BondStereo::None {
                        bond_mut.set_stereo(BondStereo::None);
                        bond_mut.set_stereo_atoms(None);
                    }
                    continue;
                } else if snapshot.order() == BondOrder::Double {
                    // BEGIN RDKIT CPP FUNCTION assignStereochemistry cleanIt
                    // BEGIN RDKIT CPP FUNCTION double-bond reset branch
                    // RDKit✔️✔️: } else if (bond->getBondType() == Bond::DOUBLE) {
                    // RDKit✔️✔️:   if (bond->getBondDir() == Bond::EITHERDOUBLE) {
                    // RDKit✔️✔️:     bond->setStereo(Bond::STEREOANY);
                    // RDKit✔️✔️:     bond->getStereoAtoms().clear();
                    // RDKit✔️✔️:     bond->setBondDir(Bond::NONE);
                    // RDKit✔️✔️:   } else if (bond->getStereo() != Bond::STEREOANY) {
                    // RDKit✔️✔️:     bond->setStereo(Bond::STEREONONE);
                    // RDKit✔️✔️:     bond->getStereoAtoms().clear();
                    // RDKit✔️✔️:   }
                    // RDKit✔️✔️: }
                    // END RDKIT CPP FUNCTION assignStereochemistry cleanIt
                    // END RDKIT CPP FUNCTION double-bond reset branch
                    if bond_mut.direction() == BondDirection::EitherDouble {
                        bond_mut.set_stereo(BondStereo::Any);
                        bond_mut.set_stereo_atoms(None);
                        bond_mut.set_direction(BondDirection::None);
                    } else if bond_mut.stereo() != BondStereo::Any {
                        bond_mut.set_stereo(BondStereo::None);
                        bond_mut.set_stereo_atoms(None);
                    }
                }
            }
        }
    }
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception keepGoing gate subset
    // RDKit❗✔️: bool hasStereoAtoms = false;
    // RDKit❗✔️: bool hasPotentialStereoAtoms = false;
    // RDKit❗✔️: for (auto atom : mol.atoms()) {
    // RDKit❗✔️:   if (cleanIt) {
    // RDKit❗✔️:     if (atom->hasProp(common_properties::_CIPCode)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_CIPCode);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (atom->hasProp(common_properties::_ChiralityPossible)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_ChiralityPossible);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!hasStereoAtoms && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:       atom->getChiralTag() != Atom::CHI_OTHER) {
    // RDKit❗✔️:     hasStereoAtoms = true;
    // RDKit❗✔️:   } else if (!hasPotentialStereoAtoms) {
    // RDKit❗✔️:     UINT_VECT ranks;
    // RDKit❗✔️:     Chirality::INT_PAIR_VECT nbrs;
    // RDKit❗✔️:     hasPotentialStereoAtoms =
    // RDKit❗✔️:         isAtomPotentialChiralCenter(atom, mol, ranks, nbrs).first;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: bool hasStereoBonds = false;
    // RDKit❗✔️: bool hasPotentialStereoBonds = false;
    // RDKit❗✔️: for (auto bond : mol.bonds()) {
    // RDKit❗✔️:   ...
    // RDKit❗✔️:   if (!hasStereoBonds && bond->getBondType() == Bond::DOUBLE) {
    // RDKit❗✔️:     bool isSpecified = false;
    // RDKit❗✔️:     ...
    // RDKit❗✔️:     if (!hasPotentialStereoBonds && !isSpecified &&
    // RDKit❗✔️:         shouldDetectDoubleBondStereo(bond)) {
    // RDKit❗✔️:       hasPotentialStereoBonds = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!cleanIt && hasStereoBonds && hasPotentialStereoBonds) {
    // RDKit❗✔️:     break;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: bool keepGoing = hasStereoAtoms | hasStereoBonds;
    // RDKit❗✔️: if (!keepGoing) {
    // RDKit❗✔️:   keepGoing = flagPossibleStereoCenters &&
    // RDKit❗✔️:               (hasPotentialStereoAtoms || hasPotentialStereoBonds);
    // RDKit❗✔️: }
    let atom_ids: Vec<AtomId> = mol.atoms().iter().map(|atom| atom.id()).collect();
    let mut has_stereo_atoms = false;
    let mut has_potential_stereo_atoms = false;
    for atom_id in atom_ids.iter().copied() {
        if clean_it && let Some(atom_mut) = mol.topology_block_mut().atoms.get_mut(atom_id.index()) {
            atom_mut.clear_prop("_CIPCode");
            atom_mut.clear_prop("_ChiralityPossible");
        }
        let atom = &mol.atoms()[atom_id.index()];
        if !has_stereo_atoms && !matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) {
            has_stereo_atoms = true;
        } else if !has_potential_stereo_atoms {
            has_potential_stereo_atoms = crate::stereo::is_atom_potential_chiral_center(mol, atom_id.index(), &[]).0;
        }
    }
    let mut has_stereo_bonds = false;
    let mut has_potential_stereo_bonds = false;
    let bond_ids: Vec<BondId> = mol.bonds().iter().map(|bond| bond.id()).collect();
    for bond_id in bond_ids.iter().copied() {
        let Some(bond) = mol.bonds().get(bond_id.index()) else {
            continue;
        };
        if !has_stereo_bonds && bond.order() == BondOrder::Double {
            let mut is_specified = false;
            for batom in [bond.begin(), bond.end()] {
                for nbr in mol.topology_block().adjacency.neighbors_of(batom.index()) {
                    let Some(nbond) = mol.bonds().get(nbr.bond.index()) else {
                        continue;
                    };
                    if matches!(
                        nbond.direction(),
                        BondDirection::EndDownRight | BondDirection::EndUpRight
                    ) {
                        has_stereo_bonds = true;
                        is_specified = true;
                        break;
                    }
                }
                if has_stereo_bonds {
                    break;
                }
            }
            if !has_potential_stereo_bonds
                && !is_specified
                && crate::stereo::should_detect_double_bond_stereo(mol, bond_id)?
            {
                has_potential_stereo_bonds = true;
            }
        }
        if !clean_it && has_stereo_bonds && has_potential_stereo_bonds {
            break;
        }
    }
    let flag_possible_stereo_centers = true;
    let mut keep_going = has_stereo_atoms || has_stereo_bonds;
    if !keep_going {
        keep_going = flag_possible_stereo_centers && (has_potential_stereo_atoms || has_potential_stereo_bonds);
    }

    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception fixed-point assignment loop subset
    // RDKit❗✔️: UINT_VECT atomRanks;
    // RDKit❗✔️: bool keepGoing = hasStereoAtoms | hasStereoBonds;
    // RDKit❗✔️: bool changedStereoAtoms, changedStereoBonds;
    // RDKit❗✔️: while (keepGoing) {
    // RDKit❗✔️:   if (hasStereoAtoms || hasPotentialStereoAtoms) {
    // RDKit❗✔️:     std::tie(hasStereoAtoms, changedStereoAtoms) =
    // RDKit❗✔️:         Chirality::assignAtomChiralCodes(mol, atomRanks, flagPossibleStereoCenters);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     changedStereoAtoms = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (hasStereoBonds || hasPotentialStereoBonds) {
    // RDKit❗✔️:     std::tie(hasStereoBonds, changedStereoBonds) =
    // RDKit❗✔️:         Chirality::assignBondStereoCodes(mol, atomRanks);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     changedStereoBonds = false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   keepGoing = (hasStereoAtoms || hasStereoBonds) &&
    // RDKit❗✔️:               (changedStereoAtoms || changedStereoBonds);
    // RDKit❗✔️:   if (keepGoing) {
    // RDKit❗✔️:     Chirality::rerankAtoms(mol, atomRanks);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION legacyStereoPerception fixed-point assignment loop subset
    let mut ranks = Vec::new();
    while keep_going {
        if ranks.is_empty()
            && (has_stereo_atoms || has_potential_stereo_atoms || has_stereo_bonds || has_potential_stereo_bonds)
        {
            ranks = crate::stereo::assign_atom_cip_ranks_in_place(mol)?;
        }
        // BEGIN RDKIT CPP FUNCTION assignAtomChiralCodes subset
        // RDKit❗✔️: std::pair<bool, bool> assignAtomChiralCodes(ROMol &mol, UINT_VECT &ranks,
        // RDKit❗✔️:                                             bool flagPossibleStereoCenters) {
        // RDKit❗✔️:   for (auto atom : mol.atoms()) {
        // RDKit❗✔️:     if (atom->hasProp(common_properties::_CIPCode)) { continue; }
        // RDKit❗✔️:     if (!ranks.size()) { assignAtomCIPRanks(mol, ranks); }
        // RDKit❗✔️:     auto [legalCenter, hasDupes] = isAtomPotentialChiralCenter(atom, mol, ranks, nbrs);
        // RDKit❗✔️:     if (legalCenter && !hasDupes && tag != CHI_UNSPECIFIED && tag != CHI_OTHER) {
        // RDKit❗✔️:       int nSwaps = atom->getPerturbationOrder(nbrIndices);
        // RDKit❗✔️:       if (nbrIndices.size() == 3 && atom->getTotalNumHs() == 1) { ++nSwaps; }
        // RDKit❗✔️:       if (nSwaps % 2) { ... flip tag ... }
        // RDKit❗✔️:       atom->setProp(common_properties::_CIPCode, cipCode);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION assignAtomChiralCodes subset
        let atom_changed = if has_stereo_atoms || has_potential_stereo_atoms {
            let (unassigned_atoms, atom_assignments, possible_atoms, atom_changed) =
                crate::stereo::assign_atom_chiral_codes_with_possible(mol, &ranks, flag_possible_stereo_centers)?;
            for atom_idx in possible_atoms {
                if let Some(atom_mut) = mol.topology_block_mut().atoms.get_mut(atom_idx) {
                    atom_mut.set_prop("_ChiralityPossible", "1");
                }
            }
            for (atom_idx, cip_code) in atom_assignments {
                if let Some(atom_mut) = mol.topology_block_mut().atoms.get_mut(atom_idx) {
                    atom_mut.set_prop("_CIPCode", cip_code);
                }
            }
            has_stereo_atoms = unassigned_atoms;
            atom_changed
        } else {
            has_stereo_atoms = false;
            false
        };

        // BEGIN RDKIT CPP FUNCTION assignBondStereoCodes subset
        // RDKit❗✔️: std::pair<bool, bool> assignBondStereoCodes(ROMol &mol, UINT_VECT &ranks) {
        // RDKit❗✔️:   for (auto dblBond : mol.bonds()) {
        // RDKit❗✔️:     if (dblBond->getBondType() == Bond::BondType::DOUBLE) {
        // RDKit❗✔️:       if (dblBond->getStereo() != Bond::BondStereo::STEREONONE) {
        // RDKit❗✔️:         continue;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (!ranks.size()) {
        // RDKit❗✔️:         assignAtomCIPRanks(mol, ranks);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       ... find highest-ranked directionality on each side ...
        // RDKit❗✔️:       ... only assign stereo when neighbor ranking distinguishes the bond ...
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION assignBondStereoCodes subset
        let bond_changed = if has_stereo_bonds || has_potential_stereo_bonds {
            let (unassigned_bonds, assignments, bond_changed) = crate::stereo::assign_bond_stereo_codes(mol, &ranks);
            for (bond_idx, stereo, begin_control, end_control) in assignments {
                let Some(bond_mut) = mol.topology_block_mut().bonds.get_mut(bond_idx) else {
                    continue;
                };
                if bond_mut.stereo() != BondStereo::None {
                    continue;
                }
                bond_mut.set_stereo_atoms(Some([AtomId::new(begin_control), AtomId::new(end_control)]));
                bond_mut.set_stereo(match stereo {
                    crate::stereo::DoubleBondStereo::E => BondStereo::E,
                    crate::stereo::DoubleBondStereo::Z => BondStereo::Z,
                    crate::stereo::DoubleBondStereo::Unknown => BondStereo::Any,
                });
            }
            has_stereo_bonds = unassigned_bonds;
            bond_changed
        } else {
            has_stereo_bonds = false;
            false
        };
        keep_going = (has_stereo_atoms || has_stereo_bonds) && (atom_changed || bond_changed);
        if keep_going {
            ranks = crate::stereo::rerank_atoms_in_place(mol, &ranks)?;
        }
    }
    // RDKit✔️✔️: if (cleanIt) {
    // Every remaining transition in legacyStereoPerception belongs to this
    // source branch. The fixed-point assignment above still runs for
    // cleanIt=false, but ring-special-case and terminal cleanup do not.
    if !clean_it {
        return Ok(());
    }

    // BEGIN RDKIT CPP FUNCTION assignStereochemistry cleanIt atom cleanup subset
    // RDKit❗✔️:   boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    // RDKit❗✔️:   Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases, atomRanks);
    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:         !Chirality::hasNonTetrahedralStereo(atom) &&
    // RDKit❗✔️:         !atom->hasProp(common_properties::_CIPCode) &&
    // RDKit❗✔️:         (!possibleSpecialCases[atom->getIdx()] ||
    // RDKit❗✔️:          !atom->hasProp(common_properties::_ringStereoAtoms))) {
    // RDKit❗✔️:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:
    // RDKit❗✔️:       // If the atom has an explicit hydrogen and no charge, that H
    // RDKit❗✔️:       // was probably put there solely because of the chirality.
    // RDKit❗✔️:       // So we'll go ahead and remove it.
    // RDKit❗✔️:       // This was Issue 194
    // RDKit❗✔️:       if (atom->getNumExplicitHs() == 1 && atom->getFormalCharge() == 0 &&
    // RDKit❗✔️:           !atom->getIsAromatic()) {
    // RDKit❗✔️:         atom->setNumExplicitHs(0);
    // RDKit❗✔️:         atom->setNoImplicit(false);
    // RDKit❗✔️:         atom->calcExplicitValence(false);
    // RDKit❗✔️:         atom->calcImplicitValence(false);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION assignStereochemistry cleanIt atom cleanup subset
    //
    // This local MolFromSmiles subset now also preserves RDKit's
    // ring-stereochemistry special cases via `findChiralAtomSpecialCases()`,
    // but it still does not port the full `findPotentialStereo()` cleanup
    // pipeline.
    // RDKit✔️✔️: if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️✔️:   VECT_INT_VECT sssrs;
    // RDKit✔️✔️:   MolOps::symmetrizeSSSR(mol, sssrs);
    // RDKit✔️✔️: }
    // `findChiralAtomSpecialCases()` mutates the source molecule's RingInfo.
    // Materialize that transition before the read-only shared helper so later
    // stereo and SMARTS consumers observe the same SymmSSSR state.
    ensure_symm_sssr_for_stereo(mol)?;
    let special_case_atoms: BTreeSet<usize> = crate::stereo::find_chiral_atom_special_cases(mol, &ranks)?
        .into_iter()
        .map(|case| case.atom_idx)
        .collect();
    let atom_ids: Vec<AtomId> = mol.atoms().iter().map(|atom| atom.id()).collect();
    let mut explicit_h_to_implicit_updates = false;
    for atom_id in atom_ids {
        let atom = &mol.atoms()[atom_id.index()];
        if matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other)
            || crate::chemistry::stereo::has_non_tetrahedral_stereo(atom)
            // RDKit✔️✔️:         !atom->hasProp(common_properties::_CIPCode) &&
            // RDKit tests the label assigned against the original SSSR basis;
            // it does not reconsider center legality after symmetrizeSSSR().
            || atom.prop("_CIPCode").is_some()
            // RDKit✔️✔️:         (!possibleSpecialCases[atom->getIdx()] ||
            // RDKit✔️✔️:          !atom->hasProp(common_properties::_ringStereoAtoms))) {
            || (special_case_atoms.contains(&atom_id.index())
                && atom.prop("_ringStereoAtoms").is_some())
        {
            continue;
        }
        let atom_mut = &mut mol.topology_block_mut().atoms[atom_id.index()];
        atom_mut.set_chiral_tag(ChiralTag::Unspecified);
        if atom_mut.explicit_hydrogens() == 1 && atom_mut.formal_charge() == 0 && !atom_mut.is_aromatic() {
            atom_mut.set_explicit_hydrogens(0);
            atom_mut.set_no_implicit(false);
            explicit_h_to_implicit_updates = true;
        }
    }
    if explicit_h_to_implicit_updates {
        // RDKit✔️✔️:         atom->calcExplicitValence(false);
        // RDKit✔️✔️:         atom->calcImplicitValence(false);
        let valence =
            crate::assign_valence_with_options(mol, crate::ValenceModel::RdkitLike, false).map_err(|_error| {
                StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
                    feature: "SMILES_STEREOCHEMISTRY_CLEANUP",
                    reason: "valence assignment failed after explicit-H cleanup",
                })
            })?;
        mol.derived_cache_mut().valence = Some(valence);
    }
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception cleanIt terminal double-bond cleanup
    // RDKit✔️✔️:       // if we have either a double bond that involves a degree one atom and
    // RDKit✔️✔️:       // that is either crossed or STEREOANY, we need to clear the stereo and,
    // RDKit✔️✔️:       // possibly, direction. This can happen with imines that just have an
    // RDKit✔️✔️:       // implicit H on the N
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:           (bond->getBondDir() == Bond::EITHERDOUBLE ||
    // RDKit✔️✔️:            bond->getStereo() == Bond::STEREOANY) &&
    // RDKit✔️✔️:           (bond->getBeginAtom()->getDegree() == 1 ||
    // RDKit✔️✔️:            bond->getEndAtom()->getDegree() == 1)) {
    // RDKit✔️✔️:         if (bond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit✔️✔️:           bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         bond->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️:       }
    // END RDKIT CPP FUNCTION legacyStereoPerception cleanIt terminal double-bond cleanup
    if clean_it {
        let bond_ids = mol.bonds().iter().map(crate::Bond::id).collect::<Vec<_>>();
        for bond_id in bond_ids {
            let Some(bond) = mol.bonds().get(bond_id.index()).cloned() else {
                continue;
            };
            let has_terminal_endpoint = mol.topology_block().adjacency.neighbors_of(bond.begin().index()).len() == 1
                || mol.topology_block().adjacency.neighbors_of(bond.end().index()).len() == 1;
            if bond.order() != BondOrder::Double
                || !matches!(bond.direction(), BondDirection::EitherDouble) && bond.stereo() != BondStereo::Any
                || !has_terminal_endpoint
            {
                continue;
            }
            let bond_mut = &mut mol.topology_block_mut().bonds[bond_id.index()];
            if bond_mut.direction() == BondDirection::EitherDouble {
                bond_mut.set_direction(BondDirection::None);
            }
            bond_mut.set_stereo(BondStereo::None);
        }
    }
    // BEGIN RDKIT CPP FUNCTION assignStereochemistry cleanIt bond-dir cleanup subset
    // RDKit✔️✔️:       // check for directionality on single bonds around
    // RDKit✔️✔️:       // double bonds without stereo. This was github #2422
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:           (bond->getStereo() == Bond::STEREOANY ||
    // RDKit✔️✔️:            bond->getStereo() == Bond::STEREONONE)) {
    // RDKit✔️✔️:         std::vector<Atom *> batoms = {bond->getBeginAtom(), bond->getEndAtom()};
    // RDKit✔️✔️:         for (auto batom : batoms) {
    // RDKit✔️✔️:           for (const auto nbrBndI : mol.atomBonds(batom)) {
    // RDKit✔️✔️:             if (nbrBndI == bond) {
    // RDKit✔️✔️:               continue;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             if ((nbrBndI->getBondDir() == Bond::ENDDOWNRIGHT ||
    // RDKit✔️✔️:                  nbrBndI->getBondDir() == Bond::ENDUPRIGHT) &&
    // RDKit✔️✔️:                 (nbrBndI->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:                  nbrBndI->getBondType() == Bond::AROMATIC)) {
    // RDKit✔️✔️:               bool okToClear = true;
    // RDKit✔️✔️:               for (const auto nbrBndJ :
    // RDKit✔️✔️:                    mol.atomBonds(nbrBndI->getOtherAtom(batom))) {
    // RDKit✔️✔️:                 if (nbrBndJ->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:                     nbrBndJ->getStereo() != Bond::STEREOANY &&
    // RDKit✔️✔️:                     nbrBndJ->getStereo() != Bond::STEREONONE) {
    // RDKit✔️✔️:                   okToClear = false;
    // RDKit✔️✔️:                   break;
    // RDKit✔️✔️:                 }
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:               if (okToClear) {
    // RDKit✔️✔️:                 nbrBndI->setBondDir(Bond::NONE);
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // END RDKIT CPP FUNCTION assignStereochemistry cleanIt bond-dir cleanup subset
    let adjacency = AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
    let double_bond_ids: Vec<BondId> = mol
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Double && matches!(bond.stereo(), BondStereo::Any | BondStereo::None))
        .map(crate::Bond::id)
        .collect();
    for bond_id in double_bond_ids {
        let Some(bond) = mol.bonds().get(bond_id.index()).cloned() else {
            continue;
        };
        for batom in [bond.begin(), bond.end()] {
            let neighbor_bond_ids: Vec<BondId> = adjacency
                .neighbors_of(batom.index())
                .iter()
                .map(|neighbor| neighbor.bond)
                .filter(|nbr_bond_id| *nbr_bond_id != bond_id)
                .collect();
            for nbr_bond_id in neighbor_bond_ids {
                let Some(nbr_bond) = mol.bonds().get(nbr_bond_id.index()) else {
                    continue;
                };
                if !matches!(
                    nbr_bond.direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                ) || !matches!(nbr_bond.order(), BondOrder::Single | BondOrder::Aromatic)
                {
                    continue;
                }
                let other_atom = bond_other_endpoint(nbr_bond, batom);
                let mut ok_to_clear = true;
                for second_nbr in adjacency.neighbors_of(other_atom.index()) {
                    let Some(second_bond) = mol.bonds().get(second_nbr.bond.index()) else {
                        continue;
                    };
                    if second_bond.order() == BondOrder::Double
                        && !matches!(second_bond.stereo(), BondStereo::Any | BondStereo::None)
                    {
                        ok_to_clear = false;
                        break;
                    }
                }
                if ok_to_clear && let Some(nbr_bond_mut) = mol.topology_block_mut().bonds.get_mut(nbr_bond_id.index()) {
                    nbr_bond_mut.set_direction(BondDirection::None);
                }
            }
        }
    }
    Ok(())
}

pub(super) fn first_2d_and_3d_conformer_ids(mol: &Molecule) -> (Option<usize>, Option<usize>) {
    // RDKit✔️✔️:   const Conformer *conf = nullptr, *conf3d = nullptr;
    // RDKit✔️✔️:   if (res && res->getNumConformers() > 0) {
    // RDKit✔️✔️:     for (unsigned int confId = 0; confId < res->getNumConformers(); ++confId) {
    // RDKit✔️✔️:       auto *testConf = &res->getConformer(confId);
    // RDKit✔️✔️:       if (!testConf->is3D()) {
    // RDKit✔️✔️:         if (conf == nullptr) {  // only take the first 2d conf
    // RDKit✔️✔️:           conf = testConf;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (conf3d == nullptr) {  // only take the first 3d conf
    // RDKit✔️✔️:           conf3d = testConf;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (conf != nullptr && conf3d != nullptr) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut first_2d = None;
    let mut first_3d = None;
    for conformer in mol.conformers_3d() {
        if conformer.is_3d() {
            if first_3d.is_none() {
                first_3d = Some(conformer.id());
            }
        } else if first_2d.is_none() {
            first_2d = Some(conformer.id());
        }

        if first_2d.is_some() && first_3d.is_some() {
            break;
        }
    }
    (first_2d, first_3d)
}

pub(super) fn apply_coordinate_free_atropisomer_assignments(
    mol: &mut Molecule,
    assignments: Vec<(BondId, BondStereo)>,
) {
    apply_atropisomer_stereo_assignments(mol, assignments);
}

pub(super) fn apply_atropisomer_stereo_assignments(mol: &mut Molecule, assignments: Vec<(BondId, BondStereo)>) {
    for (bond_id, stereo) in assignments {
        if let Some(bond) = mol.topology_block_mut().bonds.get_mut(bond_id.index()) {
            bond.set_stereo(stereo);
        }
    }
}

pub(super) fn opposite_dir(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        other => other,
    }
}

pub(super) fn swap_bond_direction_if_needed(
    target_direction: BondDirection,
    source_direction: BondDirection,
    target_begin: AtomId,
    source_begin: AtomId,
) -> BondDirection {
    // BEGIN RDKIT CPP FUNCTION swapBondDirIfNeeded
    // RDKit✔️✔️: void swapBondDirIfNeeded(Bond *bond1, const Bond *bond2) {
    // RDKit✔️✔️:   if (bond1->getBondDir() == Bond::NONE && bond2->getBondDir() != Bond::NONE) {
    // RDKit✔️✔️:     bond1->setBondDir(bond2->getBondDir());
    // RDKit✔️✔️:     if (bond1->getBeginAtom() != bond2->getBeginAtom()) {
    // RDKit✔️✔️:       switch (bond1->getBondDir()) {
    // RDKit✔️✔️:         case Bond::ENDDOWNRIGHT:
    // RDKit✔️✔️:           bond1->setBondDir(Bond::ENDUPRIGHT);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::ENDUPRIGHT:
    // RDKit✔️✔️:           bond1->setBondDir(Bond::ENDDOWNRIGHT);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION swapBondDirIfNeeded
    if target_direction != BondDirection::None || source_direction == BondDirection::None {
        return target_direction;
    }
    let mut direction = source_direction;
    if target_begin != source_begin {
        direction = opposite_dir(direction);
    }
    direction
}

#[derive(Clone, Copy)]
struct SmilesControllingBondResult {
    bond: Option<BondId>,
    obond: Option<BondId>,
    squiggle_bond_seen: bool,
    double_bond_seen: bool,
}

pub(super) fn set_bond_dir_relative_to_atom(
    mol: &mut Molecule,
    bond_id: BondId,
    center: AtomId,
    mut dir: BondDirection,
    mut reverse: bool,
) {
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
    // RDKit✔️✔️:   bond->setBondDir(dir);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION setBondDirRelativeToAtom
    let Some(bond) = mol.bonds().get(bond_id.index()) else {
        return;
    };
    if bond.begin() != center {
        reverse = !reverse;
    }
    if reverse {
        dir = opposite_dir(dir);
    }
    if let Some(bond) = mol.topology_block_mut().bonds.get_mut(bond_id.index()) {
        bond.set_direction(dir);
    }
}

pub(super) fn point_sub(a: [f64; 3], b: [f64; 3]) -> (f64, f64, f64) {
    (a[0] - b[0], a[1] - b[1], a[2] - b[2])
}

pub(super) fn point_dot(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

pub(super) fn point_cross(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (a.1 * b.2 - a.2 * b.1, a.2 * b.0 - a.0 * b.2, a.0 * b.1 - a.1 * b.0)
}

pub(super) fn point_len_sq(v: (f64, f64, f64)) -> f64 {
    point_dot(v, v)
}

pub(super) fn compute_dihedral_angle_points(p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], p4: [f64; 3]) -> f64 {
    let r12 = point_sub(p2, p1);
    let r23 = point_sub(p3, p2);
    let r34 = point_sub(p4, p3);
    let n123 = point_cross(r12, r23);
    let n234 = point_cross(r23, r34);
    let n123_len_sq = point_len_sq(n123);
    let n234_len_sq = point_len_sq(n234);
    if n123_len_sq <= 1.0e-16 || n234_len_sq <= 1.0e-16 {
        return 0.0;
    }
    let cosine = (point_dot(n123, n234) / (n123_len_sq * n234_len_sq).sqrt()).clamp(-1.0, 1.0);
    cosine.acos()
}

pub(super) fn bond_other_endpoint(bond: &crate::Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom { bond.end() } else { bond.begin() }
}

pub(super) fn controlling_bond_from_atom_for_smiles(
    mol: &Molecule,
    adjacency: &AdjacencyList,
    needs_dir: &[bool],
    single_bond_counts: &[u32],
    dbl_bond_id: BondId,
    atom: AtomId,
) -> SmilesControllingBondResult {
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
    // RDKit✔️✔️:     ...
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION controllingBondFromAtom
    let mut bond = None;
    let mut obond = None;
    let mut squiggle_bond_seen = false;
    let mut double_bond_seen = false;
    for neighbor in adjacency.neighbors_of(atom.index()) {
        let t_bond_id = neighbor.bond;
        if t_bond_id == dbl_bond_id {
            continue;
        }
        let Some(t_bond) = mol.bonds().get(t_bond_id.index()) else {
            continue;
        };
        if matches!(t_bond.order(), BondOrder::Single | BondOrder::Aromatic)
            && matches!(
                t_bond.direction(),
                BondDirection::None | BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            if bond.is_none() {
                bond = Some(t_bond_id);
            } else if needs_dir[t_bond_id.index()] {
                let current = bond.expect("bond exists");
                if single_bond_counts[t_bond_id.index()] > single_bond_counts[current.index()] {
                    obond = bond;
                    bond = Some(t_bond_id);
                } else {
                    obond = Some(t_bond_id);
                }
            } else {
                obond = bond;
                bond = Some(t_bond_id);
            }
        } else if t_bond.order() == BondOrder::Double {
            double_bond_seen = true;
        }
        if matches!(t_bond.order(), BondOrder::Single | BondOrder::Aromatic)
            && (t_bond.direction() == BondDirection::Unknown || t_bond.unknown_stereo())
        {
            squiggle_bond_seen = true;
            break;
        }
    }
    SmilesControllingBondResult {
        bond,
        obond,
        squiggle_bond_seen,
        double_bond_seen,
    }
}

pub(super) fn ensure_valence_for_stereo(mol: &mut Molecule) -> Result<(), StereoError> {
    if mol.derived_cache().valence.is_some() {
        return Ok(());
    }
    let valence = crate::assign_valence_with_options(mol, crate::ValenceModel::RdkitLike, false).map_err(|_error| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "CIP_RANKING",
            reason: "valence assignment failed before stereo perception",
        })
    })?;
    mol.derived_cache_mut().valence = Some(valence);
    Ok(())
}

pub(super) fn ensure_symm_sssr_for_stereo(mol: &mut Molecule) -> Result<(), StereoError> {
    if mol
        .derived_cache()
        .rings
        .as_ref()
        .is_some_and(crate::RingInfo::is_symm_sssr)
    {
        return Ok(());
    }
    *mol = mol.with_assigned_rings().map_err(|_error| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "RING_INFO",
            reason: "ring assignment failed before double-bond stereo perception",
        })
    })?;
    Ok(())
}

fn ensure_fast_find_rings_for_legacy_stereo(mol: &mut Molecule) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception ring-info preparation
    // RDKit✔️✔️: if (!mol.getRingInfo()->isFindFastOrBetter()) {
    if !mol
        .derived_cache()
        .rings
        .as_ref()
        .is_some_and(crate::RingInfo::is_find_fast_or_better)
    {
        // RDKit✔️✔️:   MolOps::fastFindRings(mol);
        let rings = crate::fast_find_rings(mol).map_err(|_error| {
            StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
                feature: "RING_INFO",
                reason: "fast ring assignment failed before legacy stereo perception",
            })
        })?;
        mol.derived_cache_mut().rings = Some(rings);
        // RDKit✔️✔️: }
    }
    // END RDKIT CPP FUNCTION legacyStereoPerception ring-info preparation
    Ok(())
}

pub(super) fn set_double_bond_neighbor_directions_impl(
    mol: &mut Molecule,
    conf_id: Option<usize>,
) -> Result<(), StereoError> {
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️✔️:     RDKit::MolOps::symmetrizeSSSR(mol);
    // RDKit✔️✔️:   }
    ensure_symm_sssr_for_stereo(mol)?;
    let coords_ptr = conf_id.and_then(|wanted_id| {
        mol.conformers_3d()
            .iter()
            .find(|conf| conf.id() == wanted_id)
            .map(|conf| conf.coordinates() as *const [[f64; 3]])
    });
    if conf_id.is_some() && coords_ptr.is_none() {
        return Ok(());
    }
    let adjacency = AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
    let mut single_bond_counts = vec![0_u32; mol.num_bonds()];
    let mut bonds_in_play = Vec::new();
    let mut dbl_bond_nbrs = vec![Vec::<BondId>::new(); mol.num_bonds()];
    let mut single_bond_nbrs = vec![Vec::<BondId>::new(); mol.num_bonds()];
    let mut needs_dir = vec![false; mol.num_bonds()];

    for bond in mol.bonds() {
        if !crate::stereo::is_bond_candidate_for_stereo(mol, bond.id().index()) {
            continue;
        }
        let mut is_candidate = true;
        for bond_atom in [bond.begin(), bond.end()] {
            for neighbor in adjacency.neighbors_of(bond_atom.index()) {
                let nbr_bond_id = neighbor.bond;
                let Some(nbr_bond) = mol.bonds().get(nbr_bond_id.index()) else {
                    continue;
                };
                if matches!(nbr_bond.order(), BondOrder::Single | BondOrder::Aromatic) {
                    single_bond_counts[nbr_bond_id.index()] += 1;
                    if nbr_bond.begin() == bond_atom
                        && nbr_bond.direction() == BondDirection::Unknown
                        && nbr_bond.unknown_stereo()
                    {
                        is_candidate = false;
                    } else {
                        needs_dir[bond.id().index()] = true;
                        if matches!(
                            nbr_bond.direction(),
                            BondDirection::None | BondDirection::EndDownRight | BondDirection::EndUpRight
                        ) {
                            needs_dir[nbr_bond_id.index()] = true;
                            dbl_bond_nbrs[bond.id().index()].push(nbr_bond_id);
                            if !single_bond_nbrs[nbr_bond_id.index()].contains(&bond.id()) {
                                single_bond_nbrs[nbr_bond_id.index()].push(bond.id());
                            }
                        }
                    }
                }
                if !is_candidate {
                    break;
                }
            }
            if !is_candidate {
                break;
            }
        }
        if is_candidate {
            bonds_in_play.push(bond.id());
        }
    }

    if bonds_in_play.is_empty() {
        return Ok(());
    }

    // BEGIN RDKIT CPP FUNCTION setDoubleBondNeighborDirections ordering subset
    // RDKit✔️✔️:   std::vector<std::pair<unsigned int, Bond *>> orderedBondsInPlay;
    // RDKit✔️✔️:   for (auto dblBond : bondsInPlay) {
    // RDKit✔️✔️:     unsigned int countHere =
    // RDKit✔️✔️:         std::accumulate(dblBondNbrs[dblBond->getIdx()].begin(),
    // RDKit✔️✔️:                         dblBondNbrs[dblBond->getIdx()].end(), 0);
    // RDKit✔️✔️:     if (!(mol.getRingInfo()->numBondRings(dblBond->getIdx()))) {
    // RDKit✔️✔️:       countHere *= 10;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     orderedBondsInPlay.push_back(std::make_pair(countHere, dblBond));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::sort(orderedBondsInPlay.begin(), orderedBondsInPlay.end());
    // RDKit✔️✔️:   for (pairIter = orderedBondsInPlay.rbegin();
    // RDKit✔️✔️:        pairIter != orderedBondsInPlay.rend(); ++pairIter) {
    // RDKit✔️✔️:     updateDoubleBondNeighbors(mol, pairIter->second, conf, needsDir,
    // RDKit✔️✔️:                               singleBondCounts, singleBondNbrs);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION setDoubleBondNeighborDirections ordering subset
    let ring_info = mol.derived_cache().rings.as_ref();
    let mut ordered_bonds_in_play = Vec::with_capacity(bonds_in_play.len());
    for dbl_bond_id in bonds_in_play {
        let mut count_here = dbl_bond_nbrs[dbl_bond_id.index()]
            .iter()
            .map(|bond_id| bond_id.index() as u32)
            .sum::<u32>();
        let is_ring_bond = ring_info.is_some_and(|ri| ri.num_bond_rings(dbl_bond_id) > 0);
        if !is_ring_bond {
            count_here *= 10;
        }
        ordered_bonds_in_play.push((count_here, dbl_bond_id));
    }
    ordered_bonds_in_play.sort_by(|left, right| left.cmp(right));

    for (_, dbl_bond_id) in ordered_bonds_in_play.into_iter().rev() {
        update_double_bond_neighbors(
            mol,
            dbl_bond_id,
            coords_ptr,
            &adjacency,
            &mut needs_dir,
            &single_bond_counts,
            &single_bond_nbrs,
        )?;
    }
    Ok(())
}

fn set_stereo_for_bond(
    mol: &mut Molecule,
    adjacency: &AdjacencyList,
    bond_id: BondId,
    stereo: BondStereo,
    use_cx_smiles_ordering: bool,
) {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    // RDKit✔️✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit✔️✔️:                       bool useCXSmilesOrdering) {
    // RDKit✔️✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit✔️✔️:   auto endAtom = bond->getEndAtom();
    let Some(bond) = mol.bonds().get(bond_id.index()) else {
        return;
    };
    let original_begin = bond.begin();
    let mut begin_atom = original_begin;
    let mut end_atom = bond.end();
    // RDKit✔️✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    if begin_atom > end_atom {
        // RDKit✔️✔️:     std::swap(begAtom, endAtom);
        std::mem::swap(&mut begin_atom, &mut end_atom);
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    if adjacency.neighbors_of(begin_atom.index()).len() > 1 && adjacency.neighbors_of(end_atom.index()).len() > 1 {
        // RDKit✔️✔️:     unsigned int begControl = mol.getNumAtoms();
        let mut begin_control = mol.num_atoms();
        // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
        for neighbor in adjacency.neighbors_of(begin_atom.index()) {
            // RDKit✔️✔️:       if (nbr == endAtom) {
            if neighbor.atom_index == end_atom.index() {
                // RDKit✔️✔️:         continue;
                continue;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       begControl = std::min(nbr->getIdx(), begControl);
            begin_control = begin_control.min(neighbor.atom_index);
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
        let mut end_control = if use_cx_smiles_ordering { mol.num_atoms() } else { 0 };
        // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
        for neighbor in adjacency.neighbors_of(end_atom.index()) {
            // RDKit✔️✔️:       if (nbr == begAtom) {
            if neighbor.atom_index == begin_atom.index() {
                // RDKit✔️✔️:         continue;
                continue;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
            // RDKit✔️✔️:                                        : std::max(nbr->getIdx(), endControl);
            end_control = if use_cx_smiles_ordering {
                end_control.min(neighbor.atom_index)
            } else {
                end_control.max(neighbor.atom_index)
            };
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (begAtom != bond->getBeginAtom()) {
        if begin_atom != original_begin {
            // RDKit✔️✔️:       std::swap(begControl, endControl);
            std::mem::swap(&mut begin_control, &mut end_control);
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     bond->setStereoAtoms(begControl, endControl);
        // RDKit✔️✔️:     bond->setStereo(stereo);
        let bond = &mut mol.topology_block_mut().bonds[bond_id.index()];
        bond.set_stereo_atoms(Some([AtomId::new(begin_control), AtomId::new(end_control)]));
        bond.set_stereo(stereo);
        // RDKit✔️✔️:     mol.setProp("_needsDetectBondStereo", 1);
        mol.properties_mut().set_prop("_needsDetectBondStereo", "1");
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
}

pub(super) fn update_double_bond_neighbors(
    mol: &mut Molecule,
    dbl_bond_id: BondId,
    coords_ptr: Option<*const [[f64; 3]]>,
    adjacency: &AdjacencyList,
    needs_dir: &mut [bool],
    single_bond_counts: &[u32],
    single_bond_nbrs: &[Vec<BondId>],
) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION updateDoubleBondNeighbors
    // RDKit✔️✔️: void updateDoubleBondNeighbors(ROMol &mol, Bond *dblBond, const Conformer *conf,
    // RDKit✔️✔️:                                boost::dynamic_bitset<> &needsDir,
    // RDKit✔️✔️:                                std::vector<unsigned int> &singleBondCounts,
    // RDKit✔️✔️:                                const VECT_INT_VECT &singleBondNbrs) {
    // RDKit✔️✔️:   // we want to deal only with double bonds:
    // RDKit✔️✔️:   PRECONDITION(dblBond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(dblBond->getBondType() == Bond::DOUBLE, "not a double bond");
    // RDKit✔️✔️:   if (!needsDir[dblBond->getIdx()]) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   needsDir.set(dblBond->getIdx(), 0);
    // RDKit✔️✔️:   std::vector<Bond *> followupBonds;
    // RDKit✔️✔️:   Bond *bond1 = nullptr, *obond1 = nullptr;
    // RDKit✔️✔️:   bool squiggleBondSeen = false;
    // RDKit✔️✔️:   bool doubleBondSeen = false;
    // RDKit✔️✔️:   controllingBondFromAtom(mol, needsDir, singleBondCounts, dblBond,
    // RDKit✔️✔️:                           dblBond->getBeginAtom(), bond1, obond1,
    // RDKit✔️✔️:                           squiggleBondSeen, doubleBondSeen);
    // RDKit✔️✔️:   if (squiggleBondSeen) {
    // RDKit✔️✔️:     Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!bond1) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Bond *bond2 = nullptr, *obond2 = nullptr;
    // RDKit✔️✔️:   controllingBondFromAtom(mol, needsDir, singleBondCounts, dblBond,
    // RDKit✔️✔️:                           dblBond->getEndAtom(), bond2, obond2,
    // RDKit✔️✔️:                           squiggleBondSeen, doubleBondSeen);
    // RDKit✔️✔️:   if (squiggleBondSeen) {
    // RDKit✔️✔️:     Chirality::detail::setStereoForBond(mol, dblBond, Bond::STEREOANY);
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!bond2) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ... sameTorsionDir calculation ...
    // RDKit✔️✔️:   ... direction propagation ...
    // RDKit✔️✔️:   for (Bond *oDblBond : followupBonds) {
    // RDKit✔️✔️:     updateDoubleBondNeighbors(mol, oDblBond, conf, needsDir, singleBondCounts,
    // RDKit✔️✔️:                               singleBondNbrs);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION updateDoubleBondNeighbors
    if !needs_dir[dbl_bond_id.index()] {
        return Ok(());
    }
    needs_dir[dbl_bond_id.index()] = false;
    let Some(dbl_bond) = mol.bonds().get(dbl_bond_id.index()).cloned() else {
        return Ok(());
    };
    let atom1 = dbl_bond.begin();
    let atom2 = dbl_bond.end();

    let begin_result =
        controlling_bond_from_atom_for_smiles(mol, adjacency, needs_dir, single_bond_counts, dbl_bond_id, atom1);
    if begin_result.squiggle_bond_seen {
        set_stereo_for_bond(mol, adjacency, dbl_bond_id, BondStereo::Any, false);
        return Ok(());
    }
    let Some(mut bond1) = begin_result.bond else {
        return Ok(());
    };
    let mut obond1 = begin_result.obond;

    let end_result =
        controlling_bond_from_atom_for_smiles(mol, adjacency, needs_dir, single_bond_counts, dbl_bond_id, atom2);
    if end_result.squiggle_bond_seen {
        set_stereo_for_bond(mol, adjacency, dbl_bond_id, BondStereo::Any, false);
        return Ok(());
    }
    let Some(mut bond2) = end_result.bond else {
        return Ok(());
    };
    let mut obond2 = end_result.obond;

    let same_torsion_dir = if let Some(coords_ptr) = coords_ptr {
        let coords = unsafe { &*coords_ptr };
        let begin_point = coords[atom1.index()];
        let end_point = coords[atom2.index()];
        let mut bond1_point = coords[bond_other_endpoint(&mol.bonds()[bond1.index()], atom1).index()];
        let mut bond2_point = coords[bond_other_endpoint(&mol.bonds()[bond2.index()], atom2).index()];
        let mut linear = false;
        let mut p1 = point_sub(bond1_point, begin_point);
        let mut p2 = point_sub(end_point, begin_point);
        if crate::stereo::is_linear_arrangement(p1, p2) {
            if let Some(other_bond) = obond1 {
                let swap = bond1;
                bond1 = other_bond;
                obond1 = Some(swap);
                bond1_point = coords[bond_other_endpoint(&mol.bonds()[bond1.index()], atom1).index()];
                p1 = point_sub(bond1_point, begin_point);
                if crate::stereo::is_linear_arrangement(p1, p2) {
                    linear = true;
                }
            } else {
                linear = true;
            }
        }
        if !linear {
            p1 = point_sub(bond2_point, end_point);
            p2 = point_sub(begin_point, end_point);
            if crate::stereo::is_linear_arrangement(p1, p2) {
                if let Some(other_bond) = obond2 {
                    let swap = bond2;
                    bond2 = other_bond;
                    obond2 = Some(swap);
                    bond2_point = coords[bond_other_endpoint(&mol.bonds()[bond2.index()], atom2).index()];
                    p1 = point_sub(bond2_point, begin_point);
                    if crate::stereo::is_linear_arrangement(p1, p2) {
                        linear = true;
                    }
                } else {
                    linear = true;
                }
            }
        }
        if linear {
            set_stereo_for_bond(mol, adjacency, dbl_bond_id, BondStereo::Any, false);
            return Ok(());
        }
        compute_dihedral_angle_points(bond1_point, begin_point, end_point, bond2_point) >= PI / 2.0
    } else {
        let Some(stereo_atoms) = dbl_bond.stereo_atoms() else {
            return Ok(());
        };
        let mut same_torsion_dir = match dbl_bond.stereo() {
            BondStereo::Z | BondStereo::Cis => false,
            BondStereo::E | BondStereo::Trans => true,
            _ => return Ok(()),
        };
        let bond1_atom = bond_other_endpoint(&mol.bonds()[bond1.index()], atom1);
        if bond1_atom != stereo_atoms[0] && bond1_atom != stereo_atoms[1] {
            same_torsion_dir = !same_torsion_dir;
        }
        let bond2_atom = bond_other_endpoint(&mol.bonds()[bond2.index()], atom2);
        if bond2_atom != stereo_atoms[0] && bond2_atom != stereo_atoms[1] {
            same_torsion_dir = !same_torsion_dir;
        }
        same_torsion_dir
    };

    let mut reverse_bond_dir = same_torsion_dir;
    let mut followup_bonds = Vec::new();
    if needs_dir[bond1.index()] {
        for neighbor_dbl_bond in &single_bond_nbrs[bond1.index()] {
            if needs_dir[neighbor_dbl_bond.index()] {
                followup_bonds.push(*neighbor_dbl_bond);
            }
        }
    }
    if needs_dir[bond2.index()] {
        for neighbor_dbl_bond in &single_bond_nbrs[bond2.index()] {
            if needs_dir[neighbor_dbl_bond.index()] {
                followup_bonds.push(*neighbor_dbl_bond);
            }
        }
    }
    if !needs_dir[bond1.index()] {
        if needs_dir[bond2.index()] {
            if mol.bonds()[bond1.index()].begin() != atom1 {
                reverse_bond_dir = !reverse_bond_dir;
            }
            let reference_dir = mol.bonds()[bond1.index()].direction();
            set_bond_dir_relative_to_atom(mol, bond2, atom2, reference_dir, reverse_bond_dir);
        }
    } else if !needs_dir[bond2.index()] {
        if mol.bonds()[bond2.index()].begin() != atom2 {
            reverse_bond_dir = !reverse_bond_dir;
        }
        let reference_dir = mol.bonds()[bond2.index()].direction();
        set_bond_dir_relative_to_atom(mol, bond1, atom1, reference_dir, reverse_bond_dir);
    } else {
        set_bond_dir_relative_to_atom(mol, bond1, atom1, BondDirection::EndDownRight, false);
        set_bond_dir_relative_to_atom(mol, bond2, atom2, BondDirection::EndDownRight, reverse_bond_dir);
    }
    needs_dir[bond1.index()] = false;
    needs_dir[bond2.index()] = false;
    if let Some(other_bond) = obond1
        && needs_dir[other_bond.index()]
    {
        let reference_dir = mol.bonds()[bond1.index()].direction();
        let reverse = mol.bonds()[bond1.index()].begin() == atom1;
        set_bond_dir_relative_to_atom(mol, other_bond, atom1, reference_dir, reverse);
        needs_dir[other_bond.index()] = false;
    }
    if let Some(other_bond) = obond2
        && needs_dir[other_bond.index()]
    {
        let reference_dir = mol.bonds()[bond2.index()].direction();
        let reverse = mol.bonds()[bond2.index()].begin() == atom2;
        set_bond_dir_relative_to_atom(mol, other_bond, atom2, reference_dir, reverse);
        needs_dir[other_bond.index()] = false;
    }
    for followup_bond in followup_bonds {
        update_double_bond_neighbors(
            mol,
            followup_bond,
            coords_ptr,
            adjacency,
            needs_dir,
            single_bond_counts,
            single_bond_nbrs,
        )?;
    }
    Ok(())
}

pub(crate) fn set_double_bond_neighbor_directions(mol: &mut Molecule, conf_id: usize) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION setDoubleBondNeighborDirections
    // RDKit✔️✔️: void setDoubleBondNeighborDirections(ROMol &mol, const Conformer *conf) {
    // RDKit✔️✔️:   // used to store the number of single bonds a given
    // RDKit✔️✔️:   // single bond is adjacent to
    // RDKit✔️✔️:   std::vector<unsigned int> singleBondCounts(mol.getNumBonds(), 0);
    // RDKit✔️✔️:   std::vector<Bond *> bondsInPlay;
    // RDKit✔️✔️:   // keeps track of which single bonds are adjacent to each double bond:
    // RDKit✔️✔️:   VECT_INT_VECT dblBondNbrs(mol.getNumBonds());
    // RDKit✔️✔️:   // keeps track of which double bonds are adjacent to each single bond:
    // RDKit✔️✔️:   VECT_INT_VECT singleBondNbrs(mol.getNumBonds());
    // RDKit✔️✔️:   // keeps track of which single bonds need a dir set and which double bonds
    // RDKit✔️✔️:   // need to have their neighbors' dirs set
    // RDKit✔️✔️:   boost::dynamic_bitset<> needsDir(mol.getNumBonds());
    // RDKit✔️✔️:   // find double bonds that should be considered for stereochemistry
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) { ... }
    // RDKit✔️✔️:   if (!bondsInPlay.size()) { return; }
    // RDKit✔️✔️:   std::vector<std::pair<unsigned int, Bond *>> orderedBondsInPlay;
    // RDKit✔️✔️:   for (auto dblBond : bondsInPlay) { ... }
    // RDKit✔️✔️:   std::sort(orderedBondsInPlay.begin(), orderedBondsInPlay.end());
    // RDKit✔️✔️:   for (...) {
    // RDKit✔️✔️:     updateDoubleBondNeighbors(mol, pairIter->second, conf, needsDir,
    // RDKit✔️✔️:                               singleBondCounts, singleBondNbrs);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION setDoubleBondNeighborDirections
    set_double_bond_neighbor_directions_impl(mol, Some(conf_id))
}

pub(crate) fn set_double_bond_neighbor_directions_from_stereo(mol: &mut Molecule) -> Result<(), StereoError> {
    set_double_bond_neighbor_directions_impl(mol, None)
}

pub(crate) fn clear_dir_flags(mol: &mut Molecule, only_wedge_type_bond_dirs: bool) {
    // BEGIN RDKIT CPP FUNCTION clearDirFlags
    // RDKit✔️✔️: void clearDirFlags(ROMol &mol, bool onlyWedgeTypeBondDirs) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondDir() == Bond::UNKNOWN ||
    // RDKit✔️✔️:         bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit✔️✔️:       bond->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (onlyWedgeTypeBondDirs == false ||
    // RDKit✔️✔️:         (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
    // RDKit✔️✔️:          bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:       bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION clearDirFlags
    let topology = mol.topology_block_mut();
    for bond in &mut topology.bonds {
        if matches!(bond.direction(), BondDirection::Unknown | BondDirection::EitherDouble) {
            bond.set_unknown_stereo(true);
        }
        if !only_wedge_type_bond_dirs
            || !matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            bond.set_direction(BondDirection::None);
        }
    }
}

pub(crate) fn clear_all_bond_dir_flags(mol: &mut Molecule) {
    // BEGIN RDKIT CPP FUNCTION clearAllBondDirFlags
    // RDKit✔️✔️: void clearAllBondDirFlags(ROMol &mol) { clearDirFlags(mol, false); }
    // END RDKIT CPP FUNCTION clearAllBondDirFlags
    clear_dir_flags(mol, false);
}

pub(super) fn has_stereo_bond_dir(direction: BondDirection) -> bool {
    matches!(direction, BondDirection::EndDownRight | BondDirection::EndUpRight)
}

pub(super) fn neighboring_directed_bond(mol: &Molecule, atom: AtomId) -> Option<BondId> {
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
    // Local complexity review: the topology CSR adjacency gives the same
    // O(degree(atom)) lookup and allocation-free iteration as RDKit's
    // `getAtomBonds()` range. Bond lookup remains direct by stable BondId.
    for neighbor in mol.topology_block().adjacency.neighbors_of(atom.index()) {
        let bond = &mol.bonds()[neighbor.bond.index()];
        if bond.order() != BondOrder::Double && has_stereo_bond_dir(bond.direction()) {
            return Some(bond.id());
        }
    }
    None
}

pub(crate) fn set_bond_stereo_from_directions(mol: &mut Molecule) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION setBondStereoFromDirections
    // RDKit✔️✔️: void setBondStereoFromDirections(ROMol &mol) {
    // RDKit✔️✔️:   mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️:   for (Bond *bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bond->getStereo() != Bond::STEREOANY) {
    // RDKit✔️✔️:       const Atom *stereoBondBeginAtom = bond->getBeginAtom();
    // RDKit✔️✔️:       const Atom *stereoBondEndAtom = bond->getEndAtom();
    // RDKit✔️✔️:       const Bond *directedBondAtBegin =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondBeginAtom);
    // RDKit✔️✔️:       const Bond *directedBondAtEnd =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondEndAtom);
    // RDKit✔️✔️:       if (directedBondAtBegin != nullptr && directedBondAtEnd != nullptr) {
    // RDKit✔️✔️:         unsigned beginSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtBegin->getOtherAtomIdx(stereoBondBeginAtom->getIdx());
    // RDKit✔️✔️:         unsigned endSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtEnd->getOtherAtomIdx(stereoBondEndAtom->getIdx());
    // RDKit✔️✔️:         bond->setStereoAtoms(beginSideStereoAtom, endSideStereoAtom);
    // RDKit✔️✔️:         auto beginSideBondDirection = directedBondAtBegin->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtBegin->getBeginAtom() == stereoBondBeginAtom) {
    // RDKit✔️✔️:           beginSideBondDirection = getOppositeBondDir(beginSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         auto endSideBondDirection = directedBondAtEnd->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtEnd->getEndAtom() == stereoBondEndAtom) {
    // RDKit✔️✔️:           endSideBondDirection = getOppositeBondDir(endSideBondDirection);
    // RDKit✔️✔️:         }
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
    // Local complexity review: each double bond performs two CSR adjacency
    // scans matching RDKit's `getAtomBonds()` calls. The update list is O(D)
    // auxiliary storage, where D is the number of eligible double bonds; it
    // avoids cloning bonds and keeps immutable graph reads separate from the
    // final in-place stereo writes. Overall work is O(E + sum endpoint degree),
    // matching the source traversal for the modeled topology representation.
    mol.properties_mut().clear_prop("_needsDetectBondStereo");
    let mut updates = Vec::new();
    for bond in mol.bonds() {
        if bond.order() != BondOrder::Double || bond.stereo() == BondStereo::Any {
            continue;
        }
        let begin = bond.begin();
        let end = bond.end();
        let Some(begin_dir_bond) = neighboring_directed_bond(mol, begin) else {
            continue;
        };
        let Some(end_dir_bond) = neighboring_directed_bond(mol, end) else {
            continue;
        };
        let Some(begin_dir_snapshot) = mol.bonds().get(begin_dir_bond.index()) else {
            continue;
        };
        let Some(end_dir_snapshot) = mol.bonds().get(end_dir_bond.index()) else {
            continue;
        };
        let begin_side_stereo_atom = if begin_dir_snapshot.begin() == begin {
            begin_dir_snapshot.end()
        } else {
            begin_dir_snapshot.begin()
        };
        let end_side_stereo_atom = if end_dir_snapshot.begin() == end {
            end_dir_snapshot.end()
        } else {
            end_dir_snapshot.begin()
        };
        let mut begin_side_bond_direction = begin_dir_snapshot.direction();
        if begin_dir_snapshot.begin() == begin {
            begin_side_bond_direction = opposite_dir(begin_side_bond_direction);
        }
        let mut end_side_bond_direction = end_dir_snapshot.direction();
        if end_dir_snapshot.end() == end {
            end_side_bond_direction = opposite_dir(end_side_bond_direction);
        }
        let stereo = if begin_side_bond_direction == end_side_bond_direction {
            BondStereo::Trans
        } else {
            BondStereo::Cis
        };
        updates.push((bond.id(), [begin_side_stereo_atom, end_side_stereo_atom], stereo));
    }
    let topology = mol.topology_block_mut();
    for (bond_id, stereo_atoms, stereo) in updates {
        if let Some(bond) = topology.bonds.get_mut(bond_id.index()) {
            bond.set_stereo_atoms(Some(stereo_atoms));
            bond.set_stereo(stereo);
        }
    }
    Ok(())
}

pub(crate) fn assign_stereochemistry_from_3d(mol: &mut Molecule, conf_id: usize) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION assignStereochemistryFrom3D
    // RDKit✔️✔️: void assignStereochemistryFrom3D(ROMol &mol, int confId,
    // RDKit✔️✔️:                                  bool replaceExistingTags) {
    // RDKit✔️✔️:   if (!mol.getNumConformers() || !mol.getConformer(confId).is3D()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   detectBondStereochemistry(mol, confId);
    // RDKit✔️✔️:   assignChiralTypesFrom3D(mol, confId, replaceExistingTags);
    // RDKit✔️✔️:   assignStereochemistry(mol, replaceExistingTags, true, true);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION assignStereochemistryFrom3D
    let Some(is_3d) = mol
        .conformers_3d()
        .iter()
        .find(|conf| conf.id() == conf_id)
        .map(|conf| conf.is_3d())
    else {
        return Ok(());
    };
    if !is_3d {
        return Ok(());
    }
    ensure_valence_for_stereo(mol)?;
    set_double_bond_neighbor_directions(mol, conf_id)?;
    assign_chiral_types_from_3d(mol, conf_id)?;
    set_bond_stereo_from_directions(mol)?;
    Ok(())
}

pub(super) fn clear_single_bond_dir_flags(mol: &mut Molecule, only_wedge_flags: bool) {
    // RDKit✔️✔️: void clearSingleBondDirFlags(ROMol &mol, bool onlyWedgeFlags) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:       if (bond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         bond->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!onlyWedgeFlags ||
    // RDKit✔️✔️:           (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
    // RDKit✔️✔️:            bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:         bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let topology = mol.topology_block_mut();
    for bond in &mut topology.bonds {
        if bond.order() == BondOrder::Single {
            if bond.direction() == BondDirection::Unknown {
                bond.set_unknown_stereo(true);
            }
            if !only_wedge_flags
                || !matches!(
                    bond.direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                )
            {
                bond.set_direction(BondDirection::None);
            }
        }
    }
}

pub(super) fn atrop_can_have_direction(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

pub(super) fn atrop_other_atom(mol: &Molecule, bond_id: BondId, atom: AtomId) -> Option<AtomId> {
    let bond = mol.bonds().get(bond_id.index())?;
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

pub(super) fn atrop_neighbor_bonds(mol: &Molecule, focus_atom: AtomId, atrop_bond: BondId) -> Option<Vec<BondId>> {
    let mut nbr_bonds: Vec<BondId> = mol
        .bonds()
        .iter()
        .filter(|b| b.id() != atrop_bond && (b.begin() == focus_atom || b.end() == focus_atom))
        .map(|b| b.id())
        .collect();
    if nbr_bonds.is_empty() {
        return None;
    }
    if nbr_bonds.len() == 2 {
        let other0 = atrop_other_atom(mol, nbr_bonds[0], focus_atom)?;
        let other1 = atrop_other_atom(mol, nbr_bonds[1], focus_atom)?;
        if other1.index() < other0.index() {
            nbr_bonds.swap(0, 1);
        }
    }
    Some(nbr_bonds)
}

pub(super) fn atrop_end_wedge_direction(mol: &Molecule, nbr_bonds: &[BondId]) -> (bool, BondDirection) {
    let Some(bond0) = mol.bonds().get(nbr_bonds[0].index()) else {
        return (false, BondDirection::None);
    };
    let bond1 = if nbr_bonds.len() > 1 {
        mol.bonds().get(nbr_bonds[1].index())
    } else {
        None
    };

    let dir0 = match bond0.direction() {
        BondDirection::BeginWedge | BondDirection::BeginDash => bond0.direction(),
        _ => BondDirection::None,
    };
    let dir1 = match bond1.map(|bond| bond.direction()).unwrap_or(BondDirection::None) {
        BondDirection::BeginWedge | BondDirection::BeginDash => {
            bond1.map(|bond| bond.direction()).unwrap_or(BondDirection::None)
        }
        _ => BondDirection::None,
    };

    if dir0 != BondDirection::None && dir1 != BondDirection::None && dir0 == dir1 {
        return (false, BondDirection::None);
    }
    if dir0 == BondDirection::BeginWedge || dir1 == BondDirection::BeginDash {
        return (true, BondDirection::BeginWedge);
    }
    if dir0 == BondDirection::BeginDash || dir1 == BondDirection::BeginWedge {
        return (true, BondDirection::BeginDash);
    }
    (true, BondDirection::None)
}

pub(super) fn atrop_bond_frame_of_reference(
    mol: &Molecule,
    bond_id: BondId,
    conformer: &Conformer3D,
) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let bond = mol.bonds().get(bond_id.index())?;
    // RDKit✔️✔️: xAxis = conf->getAtomPos(bond->getEndAtom()->getIdx()) -
    // RDKit✔️✔️:         conf->getAtomPos(bond->getBeginAtom()->getIdx());
    let mut x_axis = vec3_sub(
        conformer.coordinates()[bond.end().index()],
        conformer.coordinates()[bond.begin().index()],
    );
    let x_len = vec3_len(x_axis);
    if x_len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    x_axis = [x_axis[0] / x_len, x_axis[1] / x_len, x_axis[2] / x_len];
    if !conformer.is_3d() {
        let y_len = (x_axis[0] * x_axis[0] + x_axis[1] * x_axis[1]).sqrt();
        if y_len < REALLY_SMALL_BOND_LEN {
            return None;
        }
        let y_axis = [-x_axis[1] / y_len, x_axis[0] / y_len, 0.0];
        let z_axis = [0.0, 0.0, 1.0];
        return Some((x_axis, y_axis, z_axis));
    }

    let mut z_axis = if x_axis[0].abs() > REALLY_SMALL_BOND_LEN || x_axis[1].abs() > REALLY_SMALL_BOND_LEN {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let mut y_axis = vec3_cross(z_axis, x_axis);
    let y_len = vec3_len(y_axis);
    if y_len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    y_axis = [y_axis[0] / y_len, y_axis[1] / y_len, y_axis[2] / y_len];
    z_axis = vec3_cross(x_axis, y_axis);
    let z_len = vec3_len(z_axis);
    if z_len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    z_axis = [z_axis[0] / z_len, z_axis[1] / z_len, z_axis[2] / z_len];
    Some((x_axis, y_axis, z_axis))
}

pub(super) fn atrop_projected_end_vector(
    mol: &Molecule,
    conformer: &Conformer3D,
    focus_atom: AtomId,
    nbr_bonds: &[BondId],
    y_axis: [f64; 3],
    z_axis: [f64; 3],
    normalize: bool,
) -> Option<[f64; 3]> {
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let other0 = atrop_other_atom(mol, nbr_bonds[0], focus_atom)?;
    let mut bond_vec = vec3_sub(
        conformer.coordinates()[other0.index()],
        conformer.coordinates()[focus_atom.index()],
    );
    bond_vec = [0.0, vec3_dot(bond_vec, y_axis), vec3_dot(bond_vec, z_axis)];

    if nbr_bonds.len() == 2 {
        let other1 = atrop_other_atom(mol, nbr_bonds[1], focus_atom)?;
        let mut other_vec = vec3_sub(
            conformer.coordinates()[other1.index()],
            conformer.coordinates()[focus_atom.index()],
        );
        other_vec = [0.0, vec3_dot(other_vec, y_axis), vec3_dot(other_vec, z_axis)];
        if vec3_len(bond_vec) < REALLY_SMALL_BOND_LEN {
            bond_vec = [-other_vec[0], -other_vec[1], -other_vec[2]];
        } else if vec3_dot(bond_vec, other_vec) > REALLY_SMALL_BOND_LEN {
            return None;
        }
    }

    let len = vec3_len(bond_vec);
    if len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    if normalize {
        Some([bond_vec[0] / len, bond_vec[1] / len, bond_vec[2] / len])
    } else {
        Some(bond_vec)
    }
}

pub(super) fn atropisomer_stereo_without_conformer(mol: &Molecule) -> Vec<(BondId, BondStereo)> {
    let Some(rings) = mol.derived_cache().rings.as_ref() else {
        return Vec::new();
    };

    let mut candidate_bonds: Vec<BondId> = Vec::new();
    for bond in mol.bonds() {
        if !atrop_can_have_direction(bond.order()) {
            continue;
        }
        if !matches!(bond.direction(), BondDirection::BeginWedge | BondDirection::BeginDash) {
            continue;
        }
        let begin = bond.begin();
        for nbr_bond in mol.bonds() {
            if nbr_bond.id() == bond.id() {
                continue;
            }
            if (nbr_bond.begin() == begin || nbr_bond.end() == begin) && !candidate_bonds.contains(&nbr_bond.id()) {
                candidate_bonds.push(nbr_bond.id());
            }
        }
    }

    let num_atoms = mol.num_atoms();
    let mut degree = vec![0usize; num_atoms];
    for bond in mol.bonds() {
        degree[bond.begin().index()] += 1;
        degree[bond.end().index()] += 1;
    }
    let hybridization: Vec<crate::Hybridization> = mol.atoms().iter().map(|atom| atom.hybridization()).collect();

    let mut assignments = Vec::new();
    for candidate_id in candidate_bonds {
        let Some(candidate) = mol.bonds().get(candidate_id.index()) else {
            continue;
        };
        if candidate.order() != BondOrder::Single || candidate.stereo() == BondStereo::Any {
            continue;
        }
        let begin_idx = candidate.begin().index();
        let end_idx = candidate.end().index();
        let deg_begin = degree.get(begin_idx).copied().unwrap_or(0);
        let deg_end = degree.get(end_idx).copied().unwrap_or(0);
        if !(2..=3).contains(&deg_begin) || !(2..=3).contains(&deg_end) {
            continue;
        }
        if hybridization
            .get(begin_idx)
            .copied()
            .unwrap_or(crate::Hybridization::Unspecified)
            != crate::Hybridization::Sp2
            || hybridization
                .get(end_idx)
                .copied()
                .unwrap_or(crate::Hybridization::Unspecified)
                != crate::Hybridization::Sp2
        {
            continue;
        }
        if rings.num_bond_rings(candidate_id) > 0 && rings.min_bond_ring_size(candidate_id) < 8 {
            continue;
        }

        let Some(nbr0) = atrop_neighbor_bonds(mol, candidate.begin(), candidate_id) else {
            continue;
        };
        let Some(nbr1) = atrop_neighbor_bonds(mol, candidate.end(), candidate_id) else {
            continue;
        };
        for end_bond in nbr0.iter().chain(nbr1.iter()) {
            let Some(bond) = mol.bonds().get(end_bond.index()) else {
                continue;
            };
            if bond.direction() == BondDirection::Unknown {
                continue;
            }
        }

        let (has_dir0, wedge_dir0) = atrop_end_wedge_direction(mol, &nbr0);
        let (has_dir1, wedge_dir1) = atrop_end_wedge_direction(mol, &nbr1);
        if !has_dir0 || !has_dir1 || wedge_dir0 == wedge_dir1 {
            continue;
        }

        if wedge_dir0 == BondDirection::BeginWedge || wedge_dir1 == BondDirection::BeginDash {
            assignments.push((candidate_id, BondStereo::AtropCcw));
        } else if wedge_dir0 == BondDirection::BeginDash || wedge_dir1 == BondDirection::BeginWedge {
            assignments.push((candidate_id, BondStereo::AtropCw));
        }
    }

    assignments
}

pub(super) fn atropisomer_stereo_from_conformer(mol: &Molecule, conf_id: usize) -> Vec<(BondId, BondStereo)> {
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let Some(conformer) = mol.conformers_3d().iter().find(|conf| conf.id() == conf_id) else {
        return Vec::new();
    };
    let Some(rings) = mol.derived_cache().rings.as_ref() else {
        return Vec::new();
    };

    let mut candidate_bonds: Vec<BondId> = Vec::new();
    for bond in mol.bonds() {
        if !atrop_can_have_direction(bond.order()) {
            continue;
        }
        if !matches!(bond.direction(), BondDirection::BeginWedge | BondDirection::BeginDash) {
            continue;
        }
        let begin = bond.begin();
        for nbr_bond in mol.bonds() {
            if nbr_bond.id() == bond.id() {
                continue;
            }
            if (nbr_bond.begin() == begin || nbr_bond.end() == begin) && !candidate_bonds.contains(&nbr_bond.id()) {
                candidate_bonds.push(nbr_bond.id());
            }
        }
    }

    let num_atoms = mol.num_atoms();
    let mut degree = vec![0usize; num_atoms];
    for bond in mol.bonds() {
        degree[bond.begin().index()] += 1;
        degree[bond.end().index()] += 1;
    }
    let hybridization: Vec<crate::Hybridization> = mol.atoms().iter().map(|atom| atom.hybridization()).collect();

    let mut assignments = Vec::new();
    'candidate: for candidate_id in candidate_bonds {
        let Some(candidate) = mol.bonds().get(candidate_id.index()) else {
            continue;
        };
        if candidate.order() != BondOrder::Single || candidate.stereo() == BondStereo::Any {
            continue;
        }
        let begin_idx = candidate.begin().index();
        let end_idx = candidate.end().index();
        let deg_begin = degree.get(begin_idx).copied().unwrap_or(0);
        let deg_end = degree.get(end_idx).copied().unwrap_or(0);
        if !(2..=3).contains(&deg_begin) || !(2..=3).contains(&deg_end) {
            continue;
        }
        if hybridization
            .get(begin_idx)
            .copied()
            .unwrap_or(crate::Hybridization::Unspecified)
            != crate::Hybridization::Sp2
            || hybridization
                .get(end_idx)
                .copied()
                .unwrap_or(crate::Hybridization::Unspecified)
                != crate::Hybridization::Sp2
        {
            continue;
        }
        if rings.num_bond_rings(candidate_id) > 0 && rings.min_bond_ring_size(candidate_id) < 8 {
            continue;
        }

        let Some(nbr0) = atrop_neighbor_bonds(mol, candidate.begin(), candidate_id) else {
            continue;
        };
        let Some(nbr1) = atrop_neighbor_bonds(mol, candidate.end(), candidate_id) else {
            continue;
        };
        for end_bond in nbr0.iter().chain(nbr1.iter()) {
            let Some(bond) = mol.bonds().get(end_bond.index()) else {
                continue 'candidate;
            };
            if bond.direction() == BondDirection::Unknown {
                continue 'candidate;
            }
        }

        let Some((_x_axis, y_axis, z_axis)) = atrop_bond_frame_of_reference(mol, candidate_id, conformer) else {
            continue;
        };
        let mut bond_vecs = [[0.0; 3]; 2];
        for (bond_atom_index, (focus_atom, nbr_bonds)) in [(candidate.begin(), &nbr0), (candidate.end(), &nbr1)]
            .into_iter()
            .enumerate()
        {
            if !conformer.is_3d() {
                let (has_dir, bond_dir) = atrop_end_wedge_direction(mol, nbr_bonds);
                if !has_dir {
                    continue 'candidate;
                }
                let Some(mut bond_vec) =
                    atrop_projected_end_vector(mol, conformer, focus_atom, nbr_bonds, y_axis, z_axis, true)
                else {
                    continue 'candidate;
                };
                if bond_dir == BondDirection::BeginWedge {
                    bond_vec[1] *= 0.707;
                    bond_vec[2] = bond_vec[1].abs();
                } else if bond_dir == BondDirection::BeginDash {
                    bond_vec[1] *= 0.707;
                    bond_vec[2] = -bond_vec[1].abs();
                }
                bond_vecs[bond_atom_index] = bond_vec;
            } else {
                let Some(bond_vec) =
                    atrop_projected_end_vector(mol, conformer, focus_atom, nbr_bonds, y_axis, z_axis, false)
                else {
                    continue 'candidate;
                };
                if vec3_len(bond_vec) < REALLY_SMALL_BOND_LEN {
                    continue 'candidate;
                }
                bond_vecs[bond_atom_index] = bond_vec;
            }
        }

        let cross_product = vec3_cross(bond_vecs[1], bond_vecs[0]);
        if cross_product[0] > REALLY_SMALL_BOND_LEN {
            assignments.push((candidate_id, BondStereo::AtropCcw));
        } else if cross_product[0] < -REALLY_SMALL_BOND_LEN {
            assignments.push((candidate_id, BondStereo::AtropCw));
        }
    }
    assignments
}

/// Port of `MolOps::assignChiralTypesFromBondDirs` — processes wedge/dash bonds
/// to assign tetrahedral chiral types to their begin atoms using pseudo-3D
/// coordinates from the given conformer.
///
/// C++ source: `third_party/rdkit/Code/GraphMol/Chirality.cpp:3765-3812`
pub(crate) fn assign_chiral_types_from_bond_dirs(mol: &mut Molecule, conf_id: usize) {
    // RDKit✔️✔️: void assignChiralTypesFromBondDirs(ROMol &mol, const int confId,
    // RDKit✔️✔️:                                    const bool replaceExistingTags) {
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto conf = mol.getConformer(confId);
    let conf_idx = mol.conformers_3d().iter().position(|c| c.id() == conf_id);
    let Some(conf_idx) = conf_idx else {
        return;
    };
    let n_atoms = mol.num_atoms();
    let mut atoms_set = vec![false; n_atoms];
    let replace_existing_tags = true; // Called from MolFromSmiles with default true
    let bond_range = 0..mol.num_bonds();
    let implicit_hydrogens = crate::assign_valence_with_options(mol, crate::ValenceModel::RdkitLike, false)
        .ok()
        .map(|valence| valence.implicit_hydrogens);

    // Phase 1: collect bonds that might need chiral assignment.
    // We collect indices first to avoid borrow conflicts between reading
    // molecule data for pseudo-3D and mutating atoms.
    let mut chiral_assignments: Vec<(usize, Option<ChiralTag>)> = Vec::new();
    for bond_idx in bond_range {
        // RDKit✔️✔️:   for (auto &bond : mol.bonds()) {
        // RDKit✔️✔️:     const Bond::BondDir dir = bond->getBondDir();
        // RDKit✔️✔️:     Atom *atom = bond->getBeginAtom();
        let bond = &mol.topology_block().bonds[bond_idx];
        let dir = bond.direction();
        let begin_idx = bond.begin().index();
        let begin_chiral = mol.topology_block().atoms[begin_idx].chiral_tag();

        // RDKit✔️✔️:     if (dir == Bond::UNKNOWN) {
        // RDKit✔️✔️:       if (atomsSet[atom->getIdx()] || replaceExistingTags) {
        // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        // RDKit✔️✔️:         atomsSet.set(atom->getIdx());
        // RDKit✔️✔️:       }
        if dir == crate::BondDirection::Unknown {
            if atoms_set[begin_idx] || replace_existing_tags {
                chiral_assignments.push((begin_idx, Some(ChiralTag::Unspecified)));
                atoms_set[begin_idx] = true;
            }
            continue;
        }

        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        if dir != crate::BondDirection::BeginWedge && dir != crate::BondDirection::BeginDash {
            continue;
        }
        // RDKit✔️✔️:         if (atomsSet[atom->getIdx()] ||
        // RDKit✔️✔️:             (!replaceExistingTags &&
        // RDKit✔️✔️:              atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
        // RDKit✔️✔️:           continue;
        // RDKit✔️✔️:         }
        // COSMolKit: we use replaceExistingTags=true from MolFromSmiles,
        // so we only skip if already set by a prior wedge bond.
        if atoms_set[begin_idx] {
            continue;
        }

        if !replace_existing_tags && begin_chiral != ChiralTag::Unspecified {
            continue;
        }

        // RDKit✔️✔️:         Atom::ChiralType code =
        // RDKit✔️✔️:             Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
        // RDKit✔️✔️:                 .value_or(Atom::ChiralType::CHI_UNSPECIFIED);
        let conformer = &mol.conformers_3d()[conf_idx];
        let code = crate::chemistry::stereo::atom_chiral_type_from_bond_dir_pseudo_3d(mol, bond_idx, conformer)
            .unwrap_or(ChiralTag::Unspecified);

        // RDKit✔️✔️:         if (code != Atom::ChiralType::CHI_UNSPECIFIED) {
        // RDKit✔️✔️:           atomsSet.set(atom->getIdx());
        // RDKit✔️✔️:         }
        if code != ChiralTag::Unspecified {
            atoms_set[begin_idx] = true;
        }

        // RDKit✔️✔️:         atom->setChiralTag(code);
        chiral_assignments.push((begin_idx, Some(code)));

        // RDKit✔️✔️:         if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
        // RDKit✔️✔️:             atom->getNumImplicitHs() == 1) {
        // RDKit✔️✔️:           atom->setNumExplicitHs(1);
        // RDKit✔️✔️:           atom->updatePropertyCache();
        // RDKit✔️✔️:         }
        // COSMolKit: compute degree as bond count to this atom.
        let deg = mol
            .bonds()
            .iter()
            .filter(|b| b.begin().index() == begin_idx || b.end().index() == begin_idx)
            .count();
        if deg == 3 && code != ChiralTag::Unspecified {
            let has_explicit = {
                let atom = &mol.topology_block().atoms[begin_idx];
                atom.explicit_hydrogens() > 0
            };
            let has_one_implicit_h = implicit_hydrogens
                .as_ref()
                .and_then(|hydrogens| hydrogens.get(begin_idx))
                .is_some_and(|hydrogens| *hydrogens == 1);
            if !has_explicit && has_one_implicit_h {
                chiral_assignments.push((begin_idx, None)); // marker for explicit H update
            }
        }
    }

    // Phase 2: apply collected assignments
    let atoms = &mut mol.topology_block_mut().atoms;
    for (idx, assignment) in chiral_assignments {
        if let Some(tag) = assignment {
            atoms[idx].set_chiral_tag(tag);
        } else {
            // marker for explicit H update
            if atoms[idx].explicit_hydrogens() == 0 {
                atoms[idx].set_explicit_hydrogens(1);
            }
        }
    }
}

/// Compute the signed volume for a tetrahedral chiral center from 3D conformer
/// coordinates. The sign of the signed volume determines CW vs CCW chirality.
///
/// Returns `ChiralTag::Unspecified` if the geometry is ambiguous (zero volume).
pub(super) fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub(super) fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(super) fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub(super) fn vec3_len(v: [f64; 3]) -> f64 {
    vec3_dot(v, v).sqrt()
}

pub(crate) fn assign_chiral_types_from_3d(molecule: &mut Molecule, conformer_id: usize) -> Result<(), StereoError> {
    crate::chemistry::stereo::assign_chiral_types_from_3d_molecule(molecule, conformer_id as i32, true)
}
