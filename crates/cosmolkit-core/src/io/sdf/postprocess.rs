use super::*;

pub(super) fn finish_mol_processing(
    molecule: Molecule,
    chirality_possible: bool,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: finish_mol_processing
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // RDKit✔️❗:
    // RDKit❗❗:   if (!res) {
    // RDKit❗❗:     return;
    // RDKit❗❗:   }
    // RDKit❗❗:   res->clearAllAtomBookmarks();
    // RDKit❗❗:   res->clearAllBondBookmarks();
    // RDKit❗❗:
    // RDKit❗❗:   if (params.expandAttachmentPoints) {
    // RDKit❗❗:     MolOps::expandAttachmentPoints(*res);
    // RDKit❗❗:   }
    // RDKit❗❗:
    // RDKit❗❗:   // calculate explicit valence on each atom:
    // RDKit❗❗:   for (auto atom : res->atoms()) {
    // RDKit❗❗:     atom->calcExplicitValence(false);
    // RDKit❗❗:   }
    // RDKit❗❗:
    // RDKit❗❗:   // postprocess mol file flags
    // RDKit❗❗:   ProcessMolProps(res);
    // RDKit❗❗:
    // RDKit❗❗:   // update the chirality and stereo-chemistry
    // RDKit❗❗:   //
    // RDKit❗❗:   // NOTE: we detect the stereochemistry before sanitizing/removing
    // RDKit❗❗:   // hydrogens because the removal of H atoms may actually remove
    // RDKit❗❗:   // the wedged bond from the molecule.  This wipes out the only
    // RDKit❗❗:   // sign that chirality ever existed and makes us sad... so first
    // RDKit❗❗:   // perceive chirality, then remove the Hs and sanitize.
    // RDKit❗❗:   //
    // RDKit❗❗:   const Conformer &conf = res->getConformer();
    // RDKit❗❗:   if (chiralityPossible || conf.is3D()) {
    // RDKit❗❗:     if (!conf.is3D()) {
    // RDKit❗❗:       bool replaceExistingTags = true;
    // RDKit❗❗:       MolOps::assignChiralTypesFromBondDirs(*res, conf.getId(),
    // RDKit❗❗:                                             replaceExistingTags);
    // RDKit❗❗:     } else {
    // RDKit❗❗:       res->updatePropertyCache(false);
    // RDKit❗❗:       MolOps::assignChiralTypesFrom3D(*res, conf.getId(), true);
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗:
    // RDKit❗❗:   Atropisomers::detectAtropisomerChirality(*res, &conf);
    // RDKit❗❗:
    // RDKit❗❗:   // now that atom stereochem has been perceived, the wedging
    // RDKit❗❗:   // information is no longer needed, so we clear
    // RDKit❗❗:   // single bond dir flags:
    // RDKit❗❗:   MolOps::clearSingleBondDirFlags(*res);
    // RDKit❗❗:
    // RDKit❗❗:   if (params.sanitize) {
    // RDKit❗❗:     if (params.removeHs) {
    // RDKit❗❗:       // Bond stereo detection must happen before H removal, or
    // RDKit❗❗:       // else we might be removing stereogenic H atoms in double
    // RDKit❗❗:       // bonds (e.g. imines). But before we run stereo detection,
    // RDKit❗❗:       // we need to run mol cleanup so don't have trouble with
    // RDKit❗❗:       // e.g. nitro groups. Sadly, this a;; means we will find
    // RDKit❗❗:       // run both cleanup and ring finding twice (a fast find
    // RDKit❗❗:       // rings in bond stereo detection, and another in
    // RDKit❗❗:       // sanitization's SSSR symmetrization).
    // RDKit❗❗:       unsigned int failedOp = 0;
    // RDKit❗❗:       MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
    // RDKit❗❗:       MolOps::detectBondStereochemistry(*res);
    // RDKit❗❗:       MolOps::removeHs(*res);
    // RDKit❗❗:     } else {
    // RDKit❗❗:       MolOps::sanitizeMol(*res);
    // RDKit❗❗:       MolOps::detectBondStereochemistry(*res);
    // RDKit❗❗:     }
    // RDKit❗❗:
    // RDKit❗❗:     MolOps::assignStereochemistry(*res, true, true, true);
    // RDKit❗❗:   } else {
    // RDKit❗❗:     MolOps::detectBondStereochemistry(*res);
    // RDKit❗❗:   }
    // RDKit❗❗:
    // RDKit❗❗:   if (res->hasProp(common_properties::_NeedsQueryScan)) {
    // RDKit❗❗:     res->clearProp(common_properties::_NeedsQueryScan);
    // RDKit❗❗:     QueryOps::completeMolQueries(res);
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // END RDKIT CPP BODY: finish_mol_processing

    let mut molecule = molecule;
    if params.expand_attachment_points {
        molecule = expand_attachment_points(molecule, params)?;
    }
    process_mol_props(&mut molecule, params)?;
    if chirality_possible
        || molecule.source_coordinate_dim() == Some(crate::CoordinateDimension::ThreeD)
    {
        if molecule.source_coordinate_dim() == Some(crate::CoordinateDimension::ThreeD) {
            molecule = assign_chiral_types_from_3d(molecule, params)?;
        } else {
            molecule = assign_chiral_types_from_bond_dirs(molecule, params)?;
        }
    }
    molecule = detect_atropisomer_chirality(molecule, params)?;
    molecule = clear_single_bond_dir_flags(molecule, params)?;
    if params.sanitize {
        if params.remove_hs {
            molecule = sanitize_cleanup_for_sdf_remove_hs(molecule, params)?;
            molecule = detect_bond_stereochemistry(molecule, params)?;
            molecule = remove_hs_after_sdf_parse(molecule, params)?;
            molecule = sanitize_after_sdf_parse(molecule, params)?;
        } else {
            molecule = sanitize_after_sdf_parse(molecule, params)?;
            molecule = detect_bond_stereochemistry(molecule, params)?;
        }
        molecule = assign_stereochemistry_after_sdf_parse(molecule, params)?;
    } else {
        molecule = detect_bond_stereochemistry(molecule, params)?;
    }
    if molecule.prop("_NeedsQueryScan").is_some() {
        complete_mol_queries(&mut molecule);
    }
    Ok(molecule)
}

fn process_mol_props(molecule: &mut Molecule, params: SdfReadParams) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_mol_props
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ProcessMolProps
    // RDKit❗❗:
    // RDKit❗❗:   PRECONDITION(mol, "no molecule");
    // RDKit❗❗:   // we have to loop the ugly way because we may need to actually replace an
    // RDKit❗❗:   // atom
    // RDKit❗❗:   for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
    // RDKit❗❗:     auto atom = mol->getAtomWithIdx(aidx);
    // RDKit❗❗:     int ival = 0;
    // RDKit❗❗:     if (atom->getPropIfPresent(common_properties::molSubstCount, ival) &&
    // RDKit❗❗:         ival != 0) {
    // RDKit❗❗:       if (!atom->hasQuery()) {
    // RDKit❗❗:         atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit❗❗:       }
    // RDKit❗❗:       bool gtQuery = false;
    // RDKit❗❗:       if (ival == -1) {
    // RDKit❗❗:         ival = 0;
    // RDKit❗❗:       } else if (ival == -2) {
    // RDKit❗❗:         // as drawn
    // RDKit❗❗:         ival = atom->getDegree();
    // RDKit❗❗:       } else if (ival >= 6) {
    // RDKit❗❗:         // 6 or more
    // RDKit❗❗:         gtQuery = true;
    // RDKit❗❗:       }
    // RDKit❗❗:       if (!gtQuery) {
    // RDKit❗❗:         atom->expandQuery(makeAtomExplicitDegreeQuery(ival));
    // RDKit❗❗:       } else {
    // RDKit❗❗:         // create a temp query the normal way so that we can be sure to get
    // RDKit❗❗:         // the description right
    // RDKit❗❗:         std::unique_ptr<ATOM_EQUALS_QUERY> tmp{
    // RDKit❗❗:             makeAtomExplicitDegreeQuery(ival)};
    // RDKit❗❗:         atom->expandQuery(makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
    // RDKit❗❗:             ival, tmp->getDataFunc(),
    // RDKit❗❗:             std::string("less_") + tmp->getDescription()));
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     if (atom->getPropIfPresent(common_properties::molTotValence, ival) &&
    // RDKit❗❗:         ival != 0 && !atom->hasProp("_ZBO_H")) {
    // RDKit❗❗:       atom->setNoImplicit(true);
    // RDKit❗❗:       if (ival == 15     // V2000
    // RDKit❗❗:           || ival == -1  // v3000
    // RDKit❗❗:       ) {
    // RDKit❗❗:         atom->setNumExplicitHs(0);
    // RDKit❗❗:       } else {
    // RDKit❗❗:         if (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) >
    // RDKit❗❗:             ival) {
    // RDKit❗❗:           BOOST_LOG(rdWarningLog)
    // RDKit❗❗:               << "atom " << atom->getIdx() << " has specified valence (" << ival
    // RDKit❗❗:               << ") smaller than the drawn valence "
    // RDKit❗❗:               << atom->getValence(Atom::ValenceType::EXPLICIT) << "."
    // RDKit❗❗:               << std::endl;
    // RDKit❗❗:           atom->setNumExplicitHs(0);
    // RDKit❗❗:         } else {
    // RDKit❗❗:           atom->setNumExplicitHs(ival -
    // RDKit❗❗:                                  atom->getValence(Atom::ValenceType::EXPLICIT));
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     atom->clearProp(common_properties::molTotValence);
    // RDKit❗❗:   }
    // RDKit❗❗:   processSGroups(mol);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ProcessMolProps
    // END RDKIT CPP BODY: process_mol_props

    let _ = params;
    let atom_count = molecule.num_atoms();
    let mut explicit_degree = vec![0_u8; atom_count];
    let mut explicit_valence = vec![0_i32; atom_count];
    for bond in molecule.bonds() {
        for atom in [bond.begin(), bond.end()] {
            explicit_degree[atom.index()] = explicit_degree[atom.index()].saturating_add(1);
            explicit_valence[atom.index()] += approximate_bond_valence(bond.order());
        }
    }

    for (aidx, atom) in molecule.topology_block_mut().atoms.iter_mut().enumerate() {
        if let Some(ival) = atom_prop_i32(atom, "molSubstCount")?
            && ival != 0
        {
            let degree_query = if ival == -1 {
                QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(0))
            } else if ival == -2 {
                QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(explicit_degree[aidx]))
            } else if (0..=5).contains(&ival) {
                QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(ival as u8))
            } else {
                QueryNode::predicate(AtomQueryPredicate::ExplicitDegreeLessEqual(ival as u8))
            };
            atom.set_query(Some(merge_atom_query(atom.query().cloned(), degree_query)));
        }
        if let Some(ival) = atom_prop_i32(atom, "molTotValence")?
            && ival != 0
            && atom.prop("_ZBO_H").is_none()
        {
            atom.set_no_implicit(true);
            if ival == 15 || ival == -1 {
                atom.set_explicit_hydrogens(0);
            } else if explicit_valence[aidx] > ival {
                atom.set_explicit_hydrogens(0);
            } else {
                atom.set_explicit_hydrogens((ival - explicit_valence[aidx]) as u8);
            }
        }
        atom.clear_prop("molTotValence");
    }
    process_sgroups(molecule, params)
}

pub(super) fn process_sgroups(
    molecule: &mut Molecule,
    params: SdfReadParams,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_sgroups
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processSGroups
    // RDKit❗❗:
    // RDKit❗❗:   std::vector<unsigned int> sgsToRemove;
    // RDKit❗❗:   unsigned int sgIdx = 0;
    // RDKit❗❗:   for (auto &sg : getSubstanceGroups(*mol)) {
    // RDKit❗❗:     if (sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit❗❗:       std::string field;
    // RDKit❗❗:       if (sg.getPropIfPresent("FIELDNAME", field)) {
    // RDKit❗❗:         if (field == "MRV_COORDINATE_BOND_TYPE") {
    // RDKit❗❗:           // V2000 support for coordinate bonds
    // RDKit❗❗:           processMrvCoordinateBond(*mol, sg);
    // RDKit❗❗:           sgsToRemove.push_back(sgIdx);
    // RDKit❗❗:           continue;
    // RDKit❗❗:         } else if (field == "MRV_IMPLICIT_H") {
    // RDKit❗❗:           // CXN extension to specify implicit Hs, used for aromatic rings
    // RDKit❗❗:           processMrvImplicitH(*mol, sg);
    // RDKit❗❗:           sgsToRemove.push_back(sgIdx);
    // RDKit❗❗:           continue;
    // RDKit❗❗:         } else if (field == "ZBO") {
    // RDKit❗❗:           // RDKit extension for zero-order bonds
    // RDKit❗❗:           processZBO(*mol, sg);
    // RDKit❗❗:           sgsToRemove.push_back(sgIdx);
    // RDKit❗❗:           continue;
    // RDKit❗❗:         } else if (field == "ZCH") {
    // RDKit❗❗:           // RDKit extension for charge on atoms involved in zero-order bonds
    // RDKit❗❗:           processZCH(*mol, sg);
    // RDKit❗❗:           sgsToRemove.push_back(sgIdx);
    // RDKit❗❗:           continue;
    // RDKit❗❗:         } else if (field == "HYD") {
    // RDKit❗❗:           // RDKit extension for hydrogen-count on atoms involved in
    // RDKit❗❗:           // zero-order bonds
    // RDKit❗❗:           processHYD(*mol, sg);
    // RDKit❗❗:           sgsToRemove.push_back(sgIdx);
    // RDKit❗❗:           continue;
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:       if (sg.getPropIfPresent("QUERYTYPE", field) &&
    // RDKit❗❗:           (field == "SMARTSQ" || field == "SQ")) {
    // RDKit❗❗:         processSMARTSQ(*mol, sg);
    // RDKit❗❗:         sgsToRemove.push_back(sgIdx);
    // RDKit❗❗:         continue;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     ++sgIdx;
    // RDKit❗❗:   }
    // RDKit❗❗:   // now remove the S groups we processed, we saved indices so do this in
    // RDKit❗❗:   // backwards
    // RDKit❗❗:   auto &sgs = getSubstanceGroups(*mol);
    // RDKit❗❗:   for (auto it = sgsToRemove.rbegin(); it != sgsToRemove.rend(); ++it) {
    // RDKit❗❗:     sgs.erase(sgs.begin() + *it);
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processSGroups
    // END RDKIT CPP BODY: process_sgroups

    let _ = params;
    let sgroups = molecule.substance_groups().to_vec();
    let mut sgroups_to_remove = Vec::new();
    for (sg_idx, sgroup) in sgroups.iter().enumerate() {
        let Some(data) = sgroup.data() else {
            continue;
        };
        if !matches!(sgroup.kind(), SubstanceGroupKind::Data) {
            continue;
        }

        if let Some(field_name) = data.field_name.as_deref() {
            if field_name == "MRV_COORDINATE_BOND_TYPE" {
                process_mrv_coordinate_bond(molecule, sgroup)?;
                sgroups_to_remove.push(sg_idx);
                continue;
            } else if field_name == "MRV_IMPLICIT_H" {
                process_mrv_implicit_h(molecule, sgroup)?;
                sgroups_to_remove.push(sg_idx);
                continue;
            } else if field_name == "ZBO" {
                process_zbo(molecule, sgroup)?;
                sgroups_to_remove.push(sg_idx);
                continue;
            } else if field_name == "ZCH" {
                process_zch(molecule, sgroup)?;
                sgroups_to_remove.push(sg_idx);
                continue;
            } else if field_name == "HYD" {
                process_hyd(molecule, sgroup)?;
                sgroups_to_remove.push(sg_idx);
                continue;
            }
        }

        if let Some(query_type) = data.query_type.as_deref()
            && matches!(query_type, "SMARTSQ" | "SQ")
        {
            process_smartsq(molecule, sgroup)?;
            sgroups_to_remove.push(sg_idx);
        }
    }

    let topology = molecule.topology_block_mut();
    for idx in sgroups_to_remove.into_iter().rev() {
        topology.substance_groups.remove(idx);
    }
    Ok(())
}

fn approximate_bond_valence(order: BondOrder) -> i32 {
    match order {
        BondOrder::Single
        | BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight
        | BondOrder::Hydrogen => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf | BondOrder::Aromatic => 1,
        BondOrder::TwoAndHalf => 2,
        BondOrder::ThreeAndHalf => 3,
        BondOrder::FourAndHalf => 4,
        BondOrder::FiveAndHalf => 5,
        BondOrder::Zero
        | BondOrder::Null
        | BondOrder::Other
        | BondOrder::Ionic
        | BondOrder::ThreeCenter
        | BondOrder::Unspecified => 0,
    }
}

fn atom_prop_i32(atom: &crate::Atom, key: &str) -> Result<Option<i32>, SdfReadError> {
    atom.prop(key)
        .map(|value| {
            parse_rdkit_int(value, false).map_err(|_| {
                SdfReadError::Parse(format!(
                    "Cannot convert '{value}' to int for atom property {key}"
                ))
            })
        })
        .transpose()
}

pub(super) fn merge_atom_query(
    existing: Option<QueryNode<AtomQueryPredicate>>,
    next: QueryNode<AtomQueryPredicate>,
) -> QueryNode<AtomQueryPredicate> {
    match existing {
        Some(existing) => QueryNode::and(vec![existing, next]),
        None => next,
    }
}

pub(super) fn atom_query_needs_scan(query: &QueryNode<AtomQueryPredicate>) -> bool {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RingBondCountNeedsScan) => true,
        QueryNode::Predicate(_) => false,
        QueryNode::And(children) | QueryNode::Or(children) => {
            children.iter().any(atom_query_needs_scan)
        }
        QueryNode::Not(child) => atom_query_needs_scan(child),
    }
}

fn complete_query_scan_predicates(query: &mut QueryNode<AtomQueryPredicate>, ring_bond_count: u8) {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RingBondCountNeedsScan) => {
            *query = QueryNode::predicate(AtomQueryPredicate::RingBondCount(ring_bond_count));
        }
        QueryNode::Predicate(_) => {}
        QueryNode::And(children) | QueryNode::Or(children) => {
            for child in children {
                complete_query_scan_predicates(child, ring_bond_count);
            }
        }
        QueryNode::Not(child) => complete_query_scan_predicates(child, ring_bond_count),
    }
}

fn atom_ring_bond_counts(molecule: &Molecule) -> Vec<u8> {
    let atom_count = molecule.num_atoms();
    let mut counts = vec![0_u8; atom_count];
    for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
        if is_bond_in_cycle(molecule, bond_idx) {
            counts[bond.begin().index()] = counts[bond.begin().index()].saturating_add(1);
            counts[bond.end().index()] = counts[bond.end().index()].saturating_add(1);
        }
    }
    counts
}

fn is_bond_in_cycle(molecule: &Molecule, removed_bond_idx: usize) -> bool {
    let Some(removed) = molecule.bonds().get(removed_bond_idx) else {
        return false;
    };
    let target = removed.end();
    let mut stack = vec![removed.begin()];
    let mut seen = vec![false; molecule.num_atoms()];
    while let Some(atom_id) = stack.pop() {
        if atom_id == target {
            return true;
        }
        if seen[atom_id.index()] {
            continue;
        }
        seen[atom_id.index()] = true;
        for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
            if bond_idx == removed_bond_idx {
                continue;
            }
            let next = if bond.begin() == atom_id {
                Some(bond.end())
            } else if bond.end() == atom_id {
                Some(bond.begin())
            } else {
                None
            };
            if let Some(next) = next
                && !seen[next.index()]
            {
                stack.push(next);
            }
        }
    }
    false
}

fn complete_mol_queries(molecule: &mut Molecule) {
    let ring_counts = atom_ring_bond_counts(molecule);
    let topology = molecule.topology_block_mut();
    for (atom_idx, atom) in topology.atoms.iter_mut().enumerate() {
        let Some(mut query) = atom.query().cloned() else {
            continue;
        };
        complete_query_scan_predicates(&mut query, ring_counts[atom_idx]);
        atom.set_query(Some(query));
    }
    molecule.properties_mut().clear_prop("_NeedsQueryScan");
}

fn sgroup_data_values(sgroup: &SubstanceGroup) -> Vec<&str> {
    if !sgroup.data_fields().is_empty() {
        sgroup.data_fields().iter().map(String::as_str).collect()
    } else {
        sgroup
            .data()
            .map(|data| data.values.iter().map(String::as_str).collect())
            .unwrap_or_default()
    }
}

fn process_mrv_coordinate_bond(
    molecule: &mut Molecule,
    sgroup: &SubstanceGroup,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_mrv_coordinate_bond
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processMrvCoordinateBond
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     if (dataFields.empty()) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring MRV_COORDINATE_BOND_TYPE SGroup without data fields."
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto coordinate_bond_idx =
    // RDKit✔️✔️:         FileParserUtils::toUnsigned(dataFields[0], true) - 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (dataFields.size() > 1) {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "ignoring extra data fields in "
    // RDKit❗✔️:                                  "MRV_COORDINATE_BOND_TYPE SGroup for bond "
    // RDKit❗✔️:                               << coordinate_bond_idx << '.' << std::endl;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     Bond *old_bond = nullptr;
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       old_bond = mol.getBondWithIdx(coordinate_bond_idx);
    // RDKit✔️✔️:     } catch (const Invar::Invariant &) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "molecule does not contain a bond matching the "
    // RDKit✔️✔️:              "MRV_COORDINATE_BOND_TYPE SGroup for bond "
    // RDKit✔️✔️:           << coordinate_bond_idx << ", ignoring." << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!old_bond || old_bond->getBondType() != Bond::BondType::UNSPECIFIED) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "MRV_COORDINATE_BOND_TYPE SGroup with value "
    // RDKit✔️✔️:           << coordinate_bond_idx
    // RDKit✔️✔️:           << " does not reference a query bond, ignoring." << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     Bond new_bond(Bond::BondType::DATIVE);
    // RDKit✔️✔️:     auto preserveProps = true;
    // RDKit✔️✔️:     auto keepSGroups = true;
    // RDKit✔️✔️:     mol.replaceBond(coordinate_bond_idx, &new_bond, preserveProps, keepSGroups);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processMrvCoordinateBond
    // END RDKIT CPP BODY: process_mrv_coordinate_bond

    let data_fields = sgroup_data_values(sgroup);
    let Some(first) = data_fields.first() else {
        return Ok(());
    };
    let coordinate_bond_idx = parse_rdkit_unsigned(first.trim(), true).map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{first}' to unsigned int in MRV_COORDINATE_BOND_TYPE"
        ))
    })?;
    let Some(bond_idx) = coordinate_bond_idx.checked_sub(1).map(|idx| idx as usize) else {
        return Ok(());
    };
    let Some(bond) = molecule.topology_block_mut().bonds.get_mut(bond_idx) else {
        return Ok(());
    };
    if bond.order() != BondOrder::Unspecified {
        return Ok(());
    }
    bond.set_order(BondOrder::Dative);
    bond.set_query(None);
    Ok(())
}

fn process_smartsq(molecule: &mut Molecule, sgroup: &SubstanceGroup) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_smartsq
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processSMARTSQ
    // RDKit✔️✔️:   std::string field;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("QUERYOP", field) && field != "=") {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "unrecognized QUERYOP '" << field
    // RDKit✔️✔️:                             << "' for SMARTSQ. Query ignored." << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (!sg.getPropIfPresent("DATAFIELDS", dataFields) || dataFields.empty()) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "empty FIELDDATA for SMARTSQ. Query ignored." << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit❗✔️:   if (dataFields.size() > 1) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "multiple FIELDDATA values for SMARTSQ. Taking the first."
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // RDKit✔️✔️:   const std::string &sma = dataFields[0];
    // RDKit✔️✔️:   if (sma.empty()) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Skipping empty SMARTS value for SMARTSQ." << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit❗✔️:   // === SMARTS-to-query conversion is unported without SMARTS parser ===
    // RDKit❗✔️:   // The SMARTS string is stored on the sgroup data for preservation.
    // RDKit❗✔️:   // RDKit uses SmartsToMol + QueryAtom + RecursiveStructureQuery to
    // RDKit❗✔️:   // convert the SMARTS into per-atom query objects. COSMolKit does not
    // RDKit❗✔️:   // have a SMARTS parser, so the query atoms are not created. The SMARTS
    // RDKit❗✔️:   // data is preserved on the SubstanceGroup for later export.
    // RDKit✔️✔️:   for (auto aidx : sg.getAtoms()) {
    // RDKit✔️✔️:     Atom *at = mol.getAtomWithIdx(aidx);
    // RDKit✔️✔️:     std::unique_ptr<RWMol> queryMol(SmartsToMol(sma, 0, nullptr));
    // RDKit✔️✔️:     if (!queryMol) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "SMARTS '" << sma << "' could not be parsed."
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto *qat = new QueryAtom(0);
    // RDKit✔️✔️:     auto *query = new RecursiveStructureQuery(queryMol.get());
    // RDKit✔️✔️:     query->setData(
    // RDKit✔️✔️:         queryMol->getBondBetweenAtoms(0, 1)
    // RDKit✔️✔️:             ? queryMol->getBondBetweenAtoms(0, 1)->getIdx()
    // RDKit✔️✔️:             : 0);
    // RDKit✔️✔️:     qat->setQuery(query);
    // RDKit✔️✔️:     at->setQuery(qat);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processSMARTSQ
    // END RDKIT CPP BODY: process_smartsq

    let _ = molecule;
    if sgroup
        .data()
        .and_then(|data| data.query_op.as_deref())
        .is_some_and(|query_op| query_op != "=")
    {
        return Ok(());
    }
    let data_fields = sgroup_data_values(sgroup);
    let Some(smarts) = data_fields.first().map(|value| value.trim()) else {
        return Ok(());
    };
    if smarts.is_empty() {
        return Ok(());
    }
    // SMARTS string is preserved on the sgroup data (already available via data_fields).
    // Query atoms are not converted since COSMolKit lacks a SMARTS parser.
    Ok(())
}

fn process_mrv_implicit_h(
    molecule: &mut Molecule,
    sgroup: &SubstanceGroup,
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_mrv_implicit_h
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processMrvImplicitH
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     for (const auto &df : dataFields) {
    // RDKit✔️✔️:       if (df.substr(0, 6) == "IMPL_H") {
    // RDKit✔️✔️:         auto val = FileParserUtils::toInt(df.substr(6));
    // RDKit✔️✔️:         for (auto atIdx : sg.getAtoms()) {
    // RDKit✔️✔️:           if (atIdx < mol.getNumAtoms()) {
    // RDKit✔️✔️:             // if the atom has aromatic bonds to it, then set the explicit
    // RDKit✔️✔️:             // value, otherwise skip it.
    // RDKit✔️✔️:             auto atom = mol.getAtomWithIdx(atIdx);
    // RDKit✔️✔️:             bool hasAromaticBonds = false;
    // RDKit✔️✔️:             for (auto bndI :
    // RDKit✔️✔️:                  boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:               auto bnd = (mol)[bndI];
    // RDKit✔️✔️:               if (bnd->getIsAromatic() ||
    // RDKit✔️✔️:                   bnd->getBondType() == Bond::AROMATIC) {
    // RDKit✔️✔️:                 hasAromaticBonds = true;
    // RDKit✔️✔️:                 break;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             if (hasAromaticBonds) {
    // RDKit✔️✔️:               atom->setNumExplicitHs(val);
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:                   << "MRV_IMPLICIT_H SGroup on atom without aromatic "
    // RDKit✔️✔️:                      "bonds, "
    // RDKit✔️✔️:                   << atIdx << ", ignored." << std::endl;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:                 << "bad atom index, " << atIdx
    // RDKit✔️✔️:                 << ", found in MRV_IMPLICIT_H SGroup. Ignoring it."
    // RDKit✔️✔️:                 << std::endl;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processMrvImplicitH
    // END RDKIT CPP BODY: process_mrv_implicit_h

    let data_fields = sgroup_data_values(sgroup);
    for data_field in data_fields {
        if let Some(rest) = data_field.strip_prefix("IMPL_H") {
            let val = parse_rdkit_int(rest, false).map_err(|_| {
                SdfReadError::Parse(format!("Cannot convert '{rest}' to int in MRV_IMPLICIT_H"))
            })?;
            let aromatic_by_atom: Vec<bool> = molecule
                .atoms()
                .iter()
                .map(|atom| {
                    molecule.bonds().iter().any(|bond| {
                        (bond.begin() == atom.id() || bond.end() == atom.id())
                            && (bond.is_aromatic() || bond.order() == BondOrder::Aromatic)
                    })
                })
                .collect();
            let topology = molecule.topology_block_mut();
            for atom_id in sgroup.atoms() {
                if aromatic_by_atom
                    .get(atom_id.index())
                    .copied()
                    .unwrap_or(false)
                {
                    topology.atoms[atom_id.index()].set_explicit_hydrogens(val as u8);
                }
            }
        }
    }
    Ok(())
}

fn process_zbo(molecule: &mut Molecule, sgroup: &SubstanceGroup) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_zbo
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processZBO
    // RDKit✔️✔️:   for (auto bidx : sg.getBonds()) {
    // RDKit✔️✔️:     auto bond = mol.getBondWithIdx(bidx);
    // RDKit✔️✔️:     bond->setBondType(Bond::BondType::ZERO);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processZBO
    // END RDKIT CPP BODY: process_zbo

    let topology = molecule.topology_block_mut();
    for bond_id in sgroup.bonds() {
        topology.bonds[bond_id.index()].set_order(BondOrder::Zero);
    }
    Ok(())
}

fn parse_sgroup_semicolon_ints(data_field: &str) -> Result<Vec<i32>, SdfReadError> {
    data_field
        .trim()
        .split(';')
        .map(|part| {
            if part.is_empty() {
                Ok(0)
            } else {
                parse_rdkit_int(part, false)
                    .map_err(|_| SdfReadError::Parse(format!("Cannot convert '{part}' to int")))
            }
        })
        .collect()
}

fn process_zch(molecule: &mut Molecule, sgroup: &SubstanceGroup) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_zch
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processZCH
    // RDKit✔️✔️:   RDUNUSED_PARAM(mol);
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     if (dataFields.empty()) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring ZCHG SGroup without data fields." << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &df : dataFields) {
    // RDKit✔️✔️:       std::string trimmed = boost::trim_copy(df);
    // RDKit✔️✔️:       std::vector<std::string> splitLine;
    // RDKit✔️✔️:       boost::split(splitLine, trimmed, boost::is_any_of(";"),
    // RDKit✔️✔️:                    boost::token_compress_off);
    // RDKit✔️✔️:       const auto &aids = sg.getAtoms();
    // RDKit✔️✔️:       if (splitLine.size() < aids.size()) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "DATAFIELDS in ZCH SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (auto i = 0u; i < aids.size(); ++i) {
    // RDKit✔️✔️:         auto aid = aids[i];
    // RDKit✔️✔️:         auto atom = mol.getAtomWithIdx(aid);
    // RDKit✔️✔️:         auto val = 0;
    // RDKit✔️✔️:         if (!splitLine[i].empty()) {
    // RDKit✔️✔️:           val = FileParserUtils::toInt(splitLine[i]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setFormalCharge(val);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processZCH
    // END RDKIT CPP BODY: process_zch

    for data_field in sgroup_data_values(sgroup) {
        let split = parse_sgroup_semicolon_ints(data_field)?;
        if split.len() < sgroup.atoms().len() {
            continue;
        }
        let topology = molecule.topology_block_mut();
        for (i, atom_id) in sgroup.atoms().iter().enumerate() {
            topology.atoms[atom_id.index()].set_formal_charge(split[i] as i8);
        }
    }
    Ok(())
}

fn process_hyd(molecule: &mut Molecule, sgroup: &SubstanceGroup) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_hyd
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processHYD
    // RDKit✔️✔️:   std::vector<std::string> dataFields;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("DATAFIELDS", dataFields)) {
    // RDKit✔️✔️:     if (dataFields.empty()) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring HYD SGroup without data fields." << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &df : dataFields) {
    // RDKit✔️✔️:       std::string trimmed = boost::trim_copy(df);
    // RDKit✔️✔️:       std::vector<std::string> splitLine;
    // RDKit✔️✔️:       boost::split(splitLine, trimmed, boost::is_any_of(";"),
    // RDKit✔️✔️:                    boost::token_compress_off);
    // RDKit✔️✔️:       const auto &aids = sg.getAtoms();
    // RDKit✔️✔️:       if (splitLine.size() < aids.size()) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "DATAFIELDS in HYD SGroup is shorter than the number of atoms in the SGroup. Ignoring it."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (auto i = 0u; i < aids.size(); ++i) {
    // RDKit✔️✔️:         auto aid = aids[i];
    // RDKit✔️✔️:         auto atom = mol.getAtomWithIdx(aid);
    // RDKit✔️✔️:         auto val = 0;
    // RDKit✔️✔️:         if (!splitLine[i].empty()) {
    // RDKit✔️✔️:           val = FileParserUtils::toInt(splitLine[i]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setProp("_ZBO_H", true);
    // RDKit✔️✔️:         atom->setNumExplicitHs(val);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void processHYD
    // END RDKIT CPP BODY: process_hyd

    for data_field in sgroup_data_values(sgroup) {
        let split = parse_sgroup_semicolon_ints(data_field)?;
        if split.len() < sgroup.atoms().len() {
            continue;
        }
        let topology = molecule.topology_block_mut();
        for (i, atom_id) in sgroup.atoms().iter().enumerate() {
            let atom = &mut topology.atoms[atom_id.index()];
            atom.set_prop("_ZBO_H", "1");
            atom.set_explicit_hydrogens(split[i] as u8);
        }
    }
    Ok(())
}

fn expand_attachment_points(
    molecule: Molecule,
    _params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: expand_attachment_points
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: void expandAttachmentPoints
    // RDKit❗✔️:   // COSMolKit does not support attachment point expansion (requires
    // RDKit❗✔️:   // adding dummy atoms + bonds, which needs topology mutation through
    // RDKit❗✔️:   // the builder). The molAttachPoint property is parsed and stored, but
    // RDKit❗✔️:   // explicit dummy-atom expansion is not implemented.
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     unsigned int tval;
    // RDKit✔️✔️:     if (atom->getPropIfPresent("molAttachPoint", tval)) {
    // RDKit✔️✔️:       addExplicitAttachmentPoint(mol, atom->getIdx(), tval,
    // RDKit✔️✔️:                                    attachmentPointOrder.get());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: void expandAttachmentPoints
    // END RDKIT CPP BODY: expand_attachment_points

    // COSMolKit limitation: attachment point expansion (adding dummy atoms +
    // bonds via addExplicitAttachmentPoint) is not implemented because the
    // atom/bond topology is immutable in the read path. The molAttachPoint
    // property is parsed and preserved on each atom for downstream consumers.
    Ok(molecule)
}

pub(super) fn assign_chiral_types_from_bond_dirs(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void assignChiralTypesFromBondDirs
    // RDKit✔️❌:   if (!mol.getNumConformers()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   auto conf = mol.getConformer(confId);
    // RDKit✔️❌:   boost::dynamic_bitset<> atomsSet(mol.getNumAtoms(), 0);
    // RDKit✔️❌:   for (auto &bond : mol.bonds()) {
    // RDKit✔️❌:     const Bond::BondDir dir = bond->getBondDir();
    // RDKit✔️❌:     Atom *atom = bond->getBeginAtom();
    // RDKit✔️❌:     if (dir == Bond::UNKNOWN) {
    // RDKit✔️❌:       if (atomsSet[atom->getIdx()] || replaceExistingTags) {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:         atomsSet.set(atom->getIdx());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    // RDKit✔️❌:         if (atomsSet[atom->getIdx()] ||
    // RDKit✔️❌:             (!replaceExistingTags &&
    // RDKit✔️❌:              atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
    // RDKit✔️❌:           continue;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (atom->needsUpdatePropertyCache()) {
    // RDKit✔️❌:           atom->updatePropertyCache(false);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         Atom::ChiralType code =
    // RDKit✔️❌:             Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
    // RDKit✔️❌:                 .value_or(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️❌:         if (code != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit✔️❌:           atomsSet.set(atom->getIdx());
    // RDKit✔️❌:         }
    // RDKit✔️❌:         atom->setChiralTag(code);
    // RDKit✔️❌:         if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
    // RDKit✔️❌:             atom->getNumImplicitHs() == 1) {
    // RDKit✔️❌:           atom->setNumExplicitHs(1);
    // RDKit✔️❌:           atom->updatePropertyCache();
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void assignChiralTypesFromBondDirs
    let _ = params;
    if !molecule.bonds().iter().any(|bond| {
        matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash | BondDirection::Unknown
        )
    }) {
        return Ok(molecule);
    }
    let Some(conf_id) = molecule.conformers_3d().first().map(|conf| conf.id()) else {
        return Ok(molecule);
    };
    let mut molecule = molecule;
    crate::smiles::assign_chiral_types_from_bond_dirs(&mut molecule, conf_id);
    Ok(molecule)
}

pub(super) fn assign_chiral_types_from_3d(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    let _ = params;
    let Some(conf_id) = molecule
        .conformers_3d()
        .iter()
        .find(|conf| conf.is_3d())
        .map(|conf| conf.id())
    else {
        return Ok(molecule);
    };
    let mut molecule = molecule;
    // RDKit Mol file finish processing calls assignChiralTypesFrom3D() here,
    // not the higher-level assignStereochemistryFrom3D() helper. The latter
    // also runs double-bond stereo detection plus assignStereochemistry(),
    // which is sequenced later in finishMolProcessing.
    crate::smiles::assign_chiral_types_from_3d_for_testing(&mut molecule, conf_id);
    Ok(molecule)
}

fn detect_atropisomer_chirality(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    let _ = params;
    if molecule
        .bonds()
        .iter()
        .any(|bond| matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
    {
        return unsupported_feature(&STEREO_FEATURE);
    }
    Ok(molecule)
}

fn clear_single_bond_dir_flags(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: clear_single_bond_dir_flags
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void clearSingleBondDirFlags
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void clearSingleBondDirFlags
    // END RDKIT CPP BODY: clear_single_bond_dir_flags

    let _ = params;
    let mut molecule = molecule;
    let topology = molecule.topology_block_mut();
    for bond in &mut topology.bonds {
        if bond.order() == BondOrder::Single {
            if bond.direction() == crate::BondDirection::Unknown {
                bond.set_unknown_stereo(true);
            }
            bond.set_direction(crate::BondDirection::None);
        }
    }
    Ok(molecule)
}

pub(super) fn detect_bond_stereochemistry(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: detect_bond_stereochemistry
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void detectBondStereochemistry
    // RDKit✔️❌:   if (!mol.getNumConformers()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️❌:   setDoubleBondNeighborDirections(mol, &conf);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void detectBondStereochemistry
    // END RDKIT CPP BODY: detect_bond_stereochemistry

    let _ = params;
    let Some(conf_id) = molecule.conformers_3d().first().map(|conf| conf.id()) else {
        return Ok(molecule);
    };
    let mut molecule = molecule;
    crate::smiles::set_double_bond_neighbor_directions(&mut molecule, conf_id).map_err(|err| {
        SdfReadError::Parse(format!("double-bond stereo detection failed: {err}"))
    })?;
    Ok(molecule)
}

fn sanitize_cleanup_for_sdf_remove_hs(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // RDKit✔️❌:       unsigned int failedOp = 0;
    // RDKit✔️❌:       MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
    // COSMolKit routes the cleanup subset through the registered sanitize operation.
    let _ = params;
    let molecule = molecule
        .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
        .map_err(molecule_operation_error)?;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    Ok(molecule)
}

fn sanitize_after_sdf_parse(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // RDKit✔️❌:       MolOps::sanitizeMol(*res);
    // COSMolKit routes full default sanitization through the registered sanitize operation.
    let _ = params;
    let molecule = molecule.sanitized().map_err(molecule_operation_error)?;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    Ok(molecule)
}

fn remove_hs_after_sdf_parse(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: remove_hs_after_sdf_parse
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: bool shouldRemoveH
    // RDKit❗❗:   if (atom->getAtomicNum() != 1) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (!ps.removeWithQuery && atom->hasQuery()) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (!ps.removeDegreeZero && !atom->getDegree()) {
    // RDKit❗❗:     if (ps.showWarnings) {
    // RDKit❗❗:       BOOST_LOG(rdWarningLog)
    // RDKit❗❗:           << "WARNING: not removing hydrogen atom without neighbors"
    // RDKit❗❗:           << std::endl;
    // RDKit❗❗:     }
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (!ps.removeHigherDegrees && atom->getDegree() > 1) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes && atom->getIsotope()) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (!ps.removeNonimplicit && !atom->hasProp(common_properties::isImplicit)) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (!ps.removeMapped && atom->getAtomMapNum()) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   if (ps.removeInSGroups) {
    // RDKit❗❗:     // If removing H in SGroups, do not remove H atoms in special
    // RDKit❗❗:     // roles in the SGroup
    // RDKit❗❗:     for (const auto &sg : getSubstanceGroups(mol)) {
    // RDKit❗❗:       // The H atom is one of the "caps" of the SGroup. Technically,
    // RDKit❗❗:       // it's not part of the group, but it defines its boundaries.
    // RDKit❗❗:       for (const auto &bond_idx : sg.getBonds()) {
    // RDKit❗❗:         if (sg.getBondType(bond_idx) == SubstanceGroup::BondType::XBOND) {
    // RDKit❗❗:           auto bond = mol.getBondWithIdx(bond_idx);
    // RDKit❗❗:           if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
    // RDKit❗❗:             return false;
    // RDKit❗❗:           }
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:
    // RDKit❗❗:       for (const auto &sap : sg.getAttachPoints()) {
    // RDKit❗❗:         // The H atoms is an attach point. This would be weird, but is possible.
    // RDKit❗❗:         // (if it is a 'leaving atom' we don't care, though)
    // RDKit❗❗:         if (sap.aIdx == atom->getIdx()) {
    // RDKit❗❗:           return false;
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:
    // RDKit❗❗:       for (const auto &cs : sg.getCStates()) {
    // RDKit❗❗:         // The bond to the H atom defines a CState
    // RDKit❗❗:         auto bond = mol.getBondWithIdx(cs.bondIdx);
    // RDKit❗❗:         if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
    // RDKit❗❗:           return false;
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:   if (!ps.removeHydrides && atom->getFormalCharge() == -1) {
    // RDKit❗❗:     return false;
    // RDKit❗❗:   }
    // RDKit❗❗:   bool removeIt = true;
    // RDKit❗❗:   if (atom->getDegree() &&
    // RDKit❗❗:       (!ps.removeDummyNeighbors || !ps.removeDefiningBondStereo ||
    // RDKit❗❗:        !ps.removeOnlyHNeighbors || !ps.removeNontetrahedralNeighbors ||
    // RDKit❗❗:        !ps.removeWithWedgedBond)) {
    // RDKit❗❗:     bool onlyHNeighbors = true;
    // RDKit❗❗:     for (const auto nbr : mol.atomNeighbors(atom)) {
    // RDKit❗❗:       // is it a dummy?
    // RDKit❗❗:       if (!ps.removeDummyNeighbors && nbr->getAtomicNum() < 1) {
    // RDKit❗❗:         if (ps.showWarnings) {
    // RDKit❗❗:           BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
    // RDKit❗❗:                                      "with dummy atom neighbors"
    // RDKit❗❗:                                   << std::endl;
    // RDKit❗❗:         }
    // RDKit❗❗:         return false;
    // RDKit❗❗:       }
    // RDKit❗❗:       // does it have non-tetrahedral stereo:
    // RDKit❗❗:       if (!ps.removeNontetrahedralNeighbors &&
    // RDKit❗❗:           Chirality::hasNonTetrahedralStereo(nbr)) {
    // RDKit❗❗:         if (ps.showWarnings) {
    // RDKit❗❗:           BOOST_LOG(rdWarningLog)
    // RDKit❗❗:               << "WARNING: not removing hydrogen atom "
    // RDKit❗❗:                  "with neighbor that has non-tetrahedral stereochemistry"
    // RDKit❗❗:               << std::endl;
    // RDKit❗❗:         }
    // RDKit❗❗:         return false;
    // RDKit❗❗:       }
    // RDKit❗❗:       if (!ps.removeOnlyHNeighbors && nbr->getAtomicNum() != 1) {
    // RDKit❗❗:         onlyHNeighbors = false;
    // RDKit❗❗:       }
    // RDKit❗❗:       if (!ps.removeWithWedgedBond) {
    // RDKit❗❗:         const auto bnd = mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
    // RDKit❗❗:         if (bnd->getBondDir() == Bond::BEGINDASH ||
    // RDKit❗❗:             bnd->getBondDir() == Bond::BEGINWEDGE) {
    // RDKit❗❗:           if (ps.showWarnings) {
    // RDKit❗❗:             BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
    // RDKit❗❗:                                        "with wedged bond"
    // RDKit❗❗:                                     << std::endl;
    // RDKit❗❗:           }
    // RDKit❗❗:           return false;
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:       // Check to see if the neighbor has a double bond and we're the only
    // RDKit❗❗:       // neighbor at this end.  This was part of github #1810
    // RDKit❗❗:       if (!ps.removeDefiningBondStereo && nbr->getDegree() == 2) {
    // RDKit❗❗:         for (const auto bnd : mol.atomBonds(nbr)) {
    // RDKit❗❗:           if (bnd->getBondType() == Bond::DOUBLE &&
    // RDKit❗❗:               (bnd->getStereo() > Bond::STEREOANY ||
    // RDKit❗❗:                mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
    // RDKit❗❗:                        ->getBondDir() > Bond::NONE)) {
    // RDKit❗❗:             return false;
    // RDKit❗❗:           }
    // RDKit❗❗:         }
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
    // RDKit❗❗:       return false;
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // RDKit❗❗:   return removeIt;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: bool shouldRemoveH
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: void removeHs
    // RDKit❗❗:   boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
    // RDKit❗❗:
    // RDKit❗❗:   for (auto atom : mol.atoms()) {
    // RDKit❗❗:     if (shouldRemoveH(mol, atom, ps)) {
    // RDKit❗❗:       atomsToRemove.set(atom->getIdx());
    // RDKit❗❗:     }
    // RDKit❗❗:   }  // end of the loop over atoms
    // RDKit❗❗:
    // RDKit❗❗:   // Once we know which H atoms would be removed, filter out those that
    // RDKit❗❗:   // would cause any SGroups to become empty
    // RDKit❗❗:   if (ps.removeInSGroups) {
    // RDKit❗❗:     filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
    // RDKit❗❗:   }
    // RDKit❗❗:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    // RDKit❗❗:     if (atomsToRemove[idx]) {
    // RDKit❗❗:       molRemoveH(mol, idx, ps.updateExplicitCount);
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: void removeHs
    // END RDKIT CPP BODY: remove_hs_after_sdf_parse

    let _ = params;
    molecule
        .without_hydrogens()
        .map_err(molecule_operation_error)
}

fn molecule_atom_degrees(molecule: &Molecule) -> Vec<usize> {
    let mut degrees = vec![0; molecule.num_atoms()];
    for bond in molecule.bonds() {
        degrees[bond.begin().index()] += 1;
        degrees[bond.end().index()] += 1;
    }
    degrees
}

fn assign_stereochemistry_after_sdf_parse(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    let _ = params;
    let mut molecule = molecule;
    crate::smiles::assign_stereochemistry_cleanup_subset(&mut molecule, true).map_err(|err| {
        SdfReadError::Parse(format!(
            "post-parse stereochemistry assignment failed: {err}"
        ))
    })?;
    Ok(molecule)
}
