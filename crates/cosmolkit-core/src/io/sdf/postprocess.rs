use std::collections::BTreeSet;

use super::*;

pub(super) fn finish_mol_processing(
    molecule: Molecule,
    chirality_possible: bool,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: finish_mol_processing
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // RDKit✔️✔️:   if (!res) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res->clearAllAtomBookmarks();
    // RDKit✔️✔️:   res->clearAllBondBookmarks();
    // RDKit✔️✔️:
    // RDKit✔️❌:   if (params.expandAttachmentPoints) {
    // RDKit✔️❌:     MolOps::expandAttachmentPoints(*res);
    // RDKit✔️❌:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // calculate explicit valence on each atom:
    // RDKit✔️✔️:   for (auto atom : res->atoms()) {
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // postprocess mol file flags
    // RDKit✔️✔️:   ProcessMolProps(res);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // update the chirality and stereo-chemistry
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // NOTE: we detect the stereochemistry before sanitizing/removing
    // RDKit✔️✔️:   // hydrogens because the removal of H atoms may actually remove
    // RDKit✔️✔️:   // the wedged bond from the molecule.  This wipes out the only
    // RDKit✔️✔️:   // sign that chirality ever existed and makes us sad... so first
    // RDKit✔️✔️:   // perceive chirality, then remove the Hs and sanitize.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   const Conformer &conf = res->getConformer();
    // RDKit✔️✔️:   if (chiralityPossible || conf.is3D()) {
    // RDKit✔️✔️:     if (!conf.is3D()) {
    // RDKit✔️✔️:       bool replaceExistingTags = true;
    // RDKit✔️❌:       MolOps::assignChiralTypesFromBondDirs(*res, conf.getId(),
    // RDKit✔️❌:                                             replaceExistingTags);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       res->updatePropertyCache(false);
    // RDKit✔️❌:       MolOps::assignChiralTypesFrom3D(*res, conf.getId(), true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️❌:   Atropisomers::detectAtropisomerChirality(*res, &conf);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now that atom stereochem has been perceived, the wedging
    // RDKit✔️✔️:   // information is no longer needed, so we clear
    // RDKit✔️✔️:   // single bond dir flags:
    // RDKit✔️✔️:   MolOps::clearSingleBondDirFlags(*res);
    // RDKit❗❗:
    // RDKit✔️❌:   if (params.sanitize) {
    // RDKit✔️❌:     if (params.removeHs) {
    // RDKit✔️❌:       // Bond stereo detection must happen before H removal, or
    // RDKit✔️❌:       // else we might be removing stereogenic H atoms in double
    // RDKit✔️❌:       // bonds (e.g. imines). But before we run stereo detection,
    // RDKit✔️❌:       // we need to run mol cleanup so don't have trouble with
    // RDKit✔️❌:       // e.g. nitro groups. Sadly, this a;; means we will find
    // RDKit✔️❌:       // run both cleanup and ring finding twice (a fast find
    // RDKit✔️❌:       // rings in bond stereo detection, and another in
    // RDKit✔️❌:       // sanitization's SSSR symmetrization).
    // RDKit✔️❌:       unsigned int failedOp = 0;
    // RDKit✔️❌:       MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
    // RDKit✔️❌:       MolOps::detectBondStereochemistry(*res);
    // RDKit✔️❌:       MolOps::removeHs(*res);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       MolOps::sanitizeMol(*res);
    // RDKit✔️❌:       MolOps::detectBondStereochemistry(*res);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit❗❌:     MolOps::assignStereochemistry(*res, true, true, true);
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     MolOps::detectBondStereochemistry(*res);
    // RDKit✔️❌:   }
    // RDKit❗❗:
    // RDKit✔️❌:   if (res->hasProp(common_properties::_NeedsQueryScan)) {
    // RDKit✔️❌:     res->clearProp(common_properties::_NeedsQueryScan);
    // RDKit✔️❌:     QueryOps::completeMolQueries(res);
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // END RDKIT CPP BODY: finish_mol_processing

    let mut molecule = molecule;
    // RDKit's null RWMol branch and bookmark clearing are represented by the
    // Rust type boundary: this function cannot receive a null Molecule, and
    // COSMolKit has no atom/bond bookmark state in the current model.
    if params.expand_attachment_points {
        molecule = expand_attachment_points(molecule, params)?;
    }
    let explicit_valence = calculate_explicit_valence_before_mol_props(&molecule)?;
    process_mol_props(&mut molecule, params, &explicit_valence)?;
    let default_conformer_is_3d = molecule
        .conformers_3d()
        .first()
        .is_some_and(crate::Conformer3D::is_3d);
    if chirality_possible || default_conformer_is_3d {
        if default_conformer_is_3d {
            let valence = crate::assign_valence_with_options(
                &molecule,
                crate::ValenceModel::RdkitLike,
                false,
            )
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
            molecule.derived_cache_mut().valence = Some(valence);
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

fn calculate_explicit_valence_before_mol_props(
    molecule: &Molecule,
) -> Result<Vec<i32>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // RDKit✔️✔️:   // calculate explicit valence on each atom:
    // RDKit✔️✔️:   for (auto atom : res->atoms()) {
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    let atom_count = molecule.num_atoms();
    let adjacency = crate::AdjacencyList::from_topology(atom_count, molecule.bonds());
    let mut explicit_valence = Vec::with_capacity(atom_count);
    for atom_idx in 0..atom_count {
        explicit_valence.push(
            crate::valence::assign_explicit_valence_for_atom_from_parts(
                molecule.atoms(),
                molecule.bonds(),
                &adjacency,
                AtomId::new(atom_idx),
                false,
            )
            .map_err(|err| SdfReadError::Parse(err.to_string()))?,
        );
    }
    Ok(explicit_valence)
}

fn process_mol_props(
    molecule: &mut Molecule,
    params: SdfReadParams,
    explicit_valence: &[i32],
) -> Result<(), SdfReadError> {
    // BEGIN RDKIT CPP BODY: process_mol_props
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ProcessMolProps
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
    // RDKit✔️✔️:   // we have to loop the ugly way because we may need to actually replace an
    // RDKit✔️✔️:   // atom
    // RDKit✔️✔️:   for (unsigned int aidx = 0; aidx < mol->getNumAtoms(); ++aidx) {
    // RDKit✔️✔️:     auto atom = mol->getAtomWithIdx(aidx);
    // RDKit✔️✔️:     int ival = 0;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molSubstCount, ival) &&
    // RDKit✔️✔️:         ival != 0) {
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bool gtQuery = false;
    // RDKit✔️✔️:       if (ival == -1) {
    // RDKit✔️✔️:         ival = 0;
    // RDKit✔️✔️:       } else if (ival == -2) {
    // RDKit✔️✔️:         // as drawn
    // RDKit✔️✔️:         ival = atom->getDegree();
    // RDKit✔️✔️:       } else if (ival >= 6) {
    // RDKit✔️✔️:         // 6 or more
    // RDKit✔️✔️:         gtQuery = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!gtQuery) {
    // RDKit✔️✔️:         atom->expandQuery(makeAtomExplicitDegreeQuery(ival));
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // create a temp query the normal way so that we can be sure to get
    // RDKit✔️✔️:         // the description right
    // RDKit✔️✔️:         std::unique_ptr<ATOM_EQUALS_QUERY> tmp{
    // RDKit✔️✔️:             makeAtomExplicitDegreeQuery(ival)};
    // RDKit✔️✔️:         atom->expandQuery(makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(
    // RDKit✔️✔️:             ival, tmp->getDataFunc(),
    // RDKit✔️✔️:             std::string("less_") + tmp->getDescription()));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molTotValence, ival) &&
    // RDKit✔️✔️:         ival != 0 && !atom->hasProp("_ZBO_H")) {
    // RDKit✔️✔️:       atom->setNoImplicit(true);
    // RDKit✔️✔️:       if (ival == 15     // V2000
    // RDKit✔️✔️:           || ival == -1  // v3000
    // RDKit✔️✔️:       ) {
    // RDKit✔️✔️:         atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) >
    // RDKit✔️✔️:             ival) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "atom " << atom->getIdx() << " has specified valence (" << ival
    // RDKit✔️✔️:               << ") smaller than the drawn valence "
    // RDKit✔️✔️:               << atom->getValence(Atom::ValenceType::EXPLICIT) << "."
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:           atom->setNumExplicitHs(0);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           atom->setNumExplicitHs(ival -
    // RDKit✔️✔️:                                  atom->getValence(Atom::ValenceType::EXPLICIT));
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->clearProp(common_properties::molTotValence);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   processSGroups(mol);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void ProcessMolProps
    // END RDKIT CPP BODY: process_mol_props

    let _ = params;
    let atom_count = molecule.num_atoms();
    if explicit_valence.len() != atom_count {
        return Err(SdfReadError::Parse(
            "mol property explicit valence row count mismatch".into(),
        ));
    }
    let mut explicit_degree = vec![0_u8; atom_count];
    for bond in molecule.bonds() {
        for atom in [bond.begin(), bond.end()] {
            explicit_degree[atom.index()] = explicit_degree[atom.index()].saturating_add(1);
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
            } else if ival < 0 || ival > i32::from(u8::MAX) {
                QueryNode::not(QueryNode::predicate(AtomQueryPredicate::Any))
            } else {
                // RDKit LessEqualQuery stores N but matches N <= atom degree.
                QueryNode::not(QueryNode::predicate(
                    AtomQueryPredicate::ExplicitDegreeLessEqual((ival - 1) as u8),
                ))
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
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: completeQueryAndChildren
    // RDKit✔️❌: void completeQueryAndChildren(Atom::QUERYATOM_QUERY *query, Atom *tgt,
    // RDKit✔️❌:                               unsigned int magicVal) {
    // RDKit✔️❌:   PRECONDITION(query, "no query");
    // RDKit✔️❌:   PRECONDITION(tgt, "no atom");
    // RDKit✔️❌:   auto eqQuery = dynamic_cast<ATOM_EQUALS_QUERY *>(query);
    // RDKit✔️❌:   if (eqQuery) {
    // RDKit✔️❌:     if (static_cast<unsigned int>(eqQuery->getVal()) == magicVal) {
    // RDKit✔️❌:       int tgtVal = eqQuery->getDataFunc()(tgt);
    // RDKit✔️❌:       eqQuery->setVal(tgtVal);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto childIt = query->beginChildren(); childIt != query->endChildren();
    // RDKit✔️❌:        ++childIt) {
    // RDKit✔️❌:     completeQueryAndChildren(childIt->get(), tgt, magicVal);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: completeQueryAndChildren
    //
    // COSMolKit does not model RDKit query data-function pointers. The only
    // current magic-value query produced by MolBlock/SDF parsing is
    // RingBondCountNeedsScan, corresponding to RDKit's AtomRingBondCount query
    // with value 0xDEADBEEF.
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
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: queryAtomRingBondCount
    // RDKit✔️❌: int queryAtomRingBondCount(const Atom *at) {
    // RDKit✔️❌:   PRECONDITION(at, "bad atom");
    // RDKit✔️❌:   PRECONDITION(at->getOwningMol(), "no owning molecule");
    // RDKit✔️❌:   int res = 0;
    // RDKit✔️❌:   ROMol::OEDGE_ITER beg, end;
    // RDKit✔️❌:   boost::tie(beg, end) = at->getOwningMol().getAtomBonds(at);
    // RDKit✔️❌:   while (beg != end) {
    // RDKit✔️❌:     if (at->getOwningMol().getRingInfo()->numBondRings((*beg)->getIdx())) {
    // RDKit✔️❌:       ++res;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++beg;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: };
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: queryAtomRingBondCount
    let rings = molecule
        .derived_cache()
        .rings
        .clone()
        .or_else(|| molecule.derived_cache().ring_families.clone())
        .or_else(|| crate::rings::find_sssr(molecule).ok());
    let mut counts = vec![0_u8; molecule.num_atoms()];
    let Some(rings) = rings else {
        return counts;
    };
    for bond in molecule.bonds() {
        if rings.num_bond_rings(bond.id()) > 0 {
            let begin = bond.begin().index();
            let end = bond.end().index();
            counts[begin] = counts[begin].saturating_add(1);
            counts[end] = counts[end].saturating_add(1);
        }
    }
    counts
}

fn complete_mol_queries(molecule: &mut Molecule) {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: completeMolQueries
    // RDKit✔️❌: void completeMolQueries(RWMol *mol, unsigned int magicVal) {
    // RDKit✔️❌:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️❌:   for (auto atom : mol->atoms()) {
    // RDKit✔️❌:     if (atom->hasQuery()) {
    // RDKit✔️❌:       completeQueryAndChildren(atom->getQuery(), atom, magicVal);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: completeMolQueries
    //
    // RDKit clears _NeedsQueryScan before calling completeMolQueries from
    // finishMolProcessing. COSMolKit keeps the clear here because this helper
    // owns the mutable molecule property map and is called only from that branch.
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
    mut molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: expand_attachment_points
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: unsigned int addExplicitAttachmentPoint
    // RDKit✔️❌:   Atom *newAtom = nullptr;
    // RDKit✔️❌:   if (addAsQuery) {
    // RDKit✔️❌:     newAtom = new QueryAtom(0);
    // RDKit✔️❌:     newAtom->setQuery(RDKit::makeAtomNullQuery());
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     newAtom = new Atom(0);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   newAtom->setProp(common_properties::_fromAttachPoint, val);
    // RDKit✔️❌:   bool updateLabel = false;
    // RDKit✔️❌:   bool takeOwnership = true;
    // RDKit✔️❌:   auto idx = mol.addAtom(newAtom, updateLabel, takeOwnership);
    // RDKit✔️❌:   mol.addBond(atomIdx, idx, Bond::SINGLE);
    // RDKit✔️❌:   mol.getAtomWithIdx(idx)->updatePropertyCache(false);
    // RDKit✔️❌:   if (addCoords) {
    // RDKit✔️❌:     setTerminalAtomCoords(mol, idx, atomIdx);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return idx;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: unsigned int addExplicitAttachmentPoint
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: void expandAttachmentPoints
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     int value;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAttachPoint, value)) {
    // RDKit✔️✔️:       std::vector<int> tgtVals;
    // RDKit✔️✔️:       if (value == 1 || value == -1) {
    // RDKit✔️✔️:         tgtVals.push_back(1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (value == 2 || value == -1) {
    // RDKit✔️✔️:         tgtVals.push_back(2);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (tgtVals.empty()) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Invalid value for molAttachPoint: " << value << " on atom "
    // RDKit✔️✔️:             << atom->getIdx() << ". Not expanding this atttachment point."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (auto tval : tgtVals) {
    // RDKit✔️✔️:         atom->clearProp(common_properties::molAttachPoint);
    // RDKit✔️✔️:         details::addExplicitAttachmentPoint(mol, atom->getIdx(), tval,
    // RDKit✔️✔️:                                             addAsQueries, addCoords);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: void expandAttachmentPoints
    // END RDKIT CPP BODY: expand_attachment_points

    let _ = params;
    let expansion_targets = molecule
        .atoms()
        .iter()
        .enumerate()
        .map(|(atom_idx, atom)| {
            let targets = match atom_prop_i32(atom, "molAttachPoint")? {
                Some(1) => vec![1_u32],
                Some(2) => vec![2_u32],
                Some(-1) => vec![1_u32, 2_u32],
                Some(_) | None => Vec::new(),
            };
            Ok((atom_idx, targets))
        })
        .collect::<Result<Vec<_>, SdfReadError>>()?;
    for (atom_idx, targets) in expansion_targets {
        if targets.is_empty() {
            continue;
        }
        for target in targets {
            molecule.topology_block_mut().atoms[atom_idx].clear_prop("molAttachPoint");
            add_explicit_attachment_point(&mut molecule, atom_idx, target, true, true)?;
        }
    }
    Ok(molecule)
}

fn add_explicit_attachment_point(
    molecule: &mut Molecule,
    atom_idx: usize,
    val: u32,
    add_as_query: bool,
    add_coords: bool,
) -> Result<usize, SdfReadError> {
    if atom_idx >= molecule.num_atoms() {
        return Err(SdfReadError::Parse(format!(
            "attachment point source atom index {atom_idx} is out of range"
        )));
    }

    let new_idx = {
        let topology = molecule.topology_block_mut();
        let atom_id = AtomId::new(topology.atoms.len());
        let mut atom_spec =
            AtomSpec::new(Element::DUMMY).with_prop("_fromAttachPoint", val.to_string());
        if add_as_query {
            atom_spec = atom_spec.with_query(QueryNode::predicate(AtomQueryPredicate::Any));
        }
        topology
            .atoms
            .push(crate::Atom::from_spec(atom_id, atom_spec));
        let bond_id = BondId::new(topology.bonds.len());
        topology.bonds.push(crate::Bond::from_spec(
            bond_id,
            BondSpec::new(AtomId::new(atom_idx), atom_id, BondOrder::Single),
        ));
        topology.adjacency =
            crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
        atom_id.index()
    };

    if add_coords {
        add_explicit_attachment_point_coords(molecule, new_idx, atom_idx)?;
    }
    molecule.derived_cache_mut().invalidate(
        crate::DerivedState::RINGS
            | crate::DerivedState::RING_FAMILIES
            | crate::DerivedState::VALENCE
            | crate::DerivedState::AROMATICITY
            | crate::DerivedState::STEREO,
    );
    Ok(new_idx)
}

fn add_explicit_attachment_point_coords(
    molecule: &mut Molecule,
    atom_idx: usize,
    other_idx: usize,
) -> Result<(), SdfReadError> {
    let adjacency = attachment_virtual_adjacency(molecule.bonds(), molecule.num_atoms());
    let mut coordinates = molecule.take_coordinate_block_or_clone();
    {
        let read_parts = crate::read_parts::MoleculeReadParts::from_molecule(&*molecule);
        for conformer in &mut coordinates.conformers_3d {
            let coord = crate::operations::ops::hydrogens::add_hs_set_terminal_atom_coord(
                read_parts,
                &adjacency,
                conformer.coordinates(),
                atom_idx,
                other_idx,
                conformer.is_3d(),
            )
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
            conformer.push_coord(coord);
        }
        for conformer in &mut coordinates.conformers_2d {
            let coords_3d = conformer
                .coordinates()
                .iter()
                .map(|coord| [coord[0], coord[1], 0.0])
                .collect::<Vec<_>>();
            let coord = crate::operations::ops::hydrogens::add_hs_set_terminal_atom_coord(
                read_parts, &adjacency, &coords_3d, atom_idx, other_idx, false,
            )
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
            conformer.push_coord([coord[0], coord[1]]);
        }
    }
    molecule.replace_coordinate_block(coordinates);
    Ok(())
}

fn attachment_virtual_adjacency(
    bonds: &[crate::Bond],
    atom_count: usize,
) -> Vec<Vec<(usize, Option<BondId>)>> {
    let mut adjacency = vec![Vec::new(); atom_count];
    for bond in bonds {
        adjacency[bond.begin().index()].push((bond.end().index(), Some(bond.id())));
        adjacency[bond.end().index()].push((bond.begin().index(), Some(bond.id())));
    }
    adjacency
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
    let fallback_2d_conformer = molecule.conformers_2d().first().map(|conformer| {
        crate::Conformer3D::new(
            conformer.id(),
            conformer
                .coordinates()
                .iter()
                .map(|coord| [coord[0], coord[1], 0.0])
                .collect(),
            false,
        )
    });
    let conformer = molecule
        .conformers_3d()
        .first()
        .or(fallback_2d_conformer.as_ref());
    let Some(conformer) = conformer else {
        return Ok(molecule);
    };
    let replace_existing_tags = true;
    let atom_count = molecule.num_atoms();
    let mut atoms_set = vec![false; atom_count];
    let mut atom_degree = vec![0_usize; atom_count];
    for bond in molecule.bonds() {
        atom_degree[bond.begin().index()] += 1;
        atom_degree[bond.end().index()] += 1;
    }
    let valence =
        crate::assign_valence_with_options(&molecule, crate::ValenceModel::RdkitLike, false)
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;
    let mut chiral_tag_updates = vec![None; atom_count];
    let mut promote_implicit_h = vec![false; atom_count];

    for bond_idx in 0..molecule.num_bonds() {
        let bond = &molecule.bonds()[bond_idx];
        let dir = bond.direction();
        let atom_idx = bond.begin().index();

        if dir == BondDirection::Unknown {
            if atoms_set[atom_idx] || replace_existing_tags {
                chiral_tag_updates[atom_idx] = Some(crate::ChiralTag::Unspecified);
                atoms_set[atom_idx] = true;
            }
            continue;
        }

        if !matches!(dir, BondDirection::BeginWedge | BondDirection::BeginDash) {
            continue;
        }
        if atoms_set[atom_idx]
            || (!replace_existing_tags
                && molecule.atoms()[atom_idx].chiral_tag() != crate::ChiralTag::Unspecified)
        {
            continue;
        }

        let code = crate::chemistry::stereo::atom_chiral_type_from_bond_dir_pseudo_3d(
            &molecule, bond_idx, conformer,
        )
        .unwrap_or(crate::ChiralTag::Unspecified);
        if code != crate::ChiralTag::Unspecified {
            atoms_set[atom_idx] = true;
        }
        chiral_tag_updates[atom_idx] = Some(code);

        if atom_degree[atom_idx] == 3
            && molecule.atoms()[atom_idx].explicit_hydrogens() == 0
            && valence.implicit_hydrogens[atom_idx] == 1
        {
            promote_implicit_h[atom_idx] = true;
        }
    }

    let mut molecule = molecule;
    let mut changed = false;
    {
        let atoms = &mut molecule.topology_block_mut().atoms;
        for atom_idx in 0..atom_count {
            if let Some(tag) = chiral_tag_updates[atom_idx] {
                atoms[atom_idx].set_chiral_tag(tag);
                changed = true;
            }
            if promote_implicit_h[atom_idx] && atoms[atom_idx].explicit_hydrogens() == 0 {
                atoms[atom_idx].set_explicit_hydrogens(1);
                changed = true;
            }
        }
    }
    if changed {
        molecule.derived_cache_mut().invalidate(
            crate::DerivedState::VALENCE
                | crate::DerivedState::STEREO
                | crate::DerivedState::DRAWING
                | crate::DerivedState::FINGERPRINT,
        );
    }
    Ok(molecule)
}

pub(super) fn assign_chiral_types_from_3d(
    mut molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    let _ = params;
    crate::chemistry::stereo::assign_chiral_types_from_3d_molecule(&mut molecule, -1, true)
        .map_err(|error| SdfReadError::Parse(error.to_string()))?;
    Ok(molecule)
}
fn detect_atropisomer_chirality(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    // BEGIN RDKIT CPP BODY: detect_atropisomer_chirality
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: void detectAtropisomerChirality
    // RDKit✔️❌:   PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
    // RDKit✔️❌:                "conformer does not belong to molecule");
    // RDKit✔️❌:
    // RDKit✔️❌:   std::set<Bond *> bondsToTry;
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     if (canHaveDirection(*bond) &&
    // RDKit✔️❌:         (bond->getBondDir() == Bond::BondDir::BEGINDASH ||
    // RDKit✔️❌:          bond->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
    // RDKit✔️❌:       for (const auto &nbrBond : mol.atomBonds(bond->getBeginAtom())) {
    // RDKit✔️❌:         if (nbrBond == bond) {
    // RDKit✔️❌:           continue;  // a bond is NOT its own neighbor
    // RDKit✔️❌:         }
    // RDKit✔️❌:         bondsToTry.insert(nbrBond);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (bondsToTry.empty()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // First, do a simple check with TotalDegree to see if any bonds might be
    // RDKit✔️❌:   // candidates before doing the expensive hybridization calculation.
    // RDKit✔️❌:   bool anyBondPassesDegreeCheck = false;
    // RDKit✔️❌:   for (auto bondToTry : bondsToTry) {
    // RDKit✔️❌:     if (bondToTry->getBeginAtom()->needsUpdatePropertyCache()) {
    // RDKit✔️❌:       bondToTry->getBeginAtom()->updatePropertyCache(false);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (bondToTry->getEndAtom()->needsUpdatePropertyCache()) {
    // RDKit✔️❌:       bondToTry->getEndAtom()->updatePropertyCache(false);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (bondToTry->getBondType() == Bond::SINGLE &&
    // RDKit✔️❌:         bondToTry->getStereo() != Bond::BondStereo::STEREOANY &&
    // RDKit✔️❌:         bondToTry->getBeginAtom()->getTotalDegree() >= 2 &&
    // RDKit✔️❌:         bondToTry->getBeginAtom()->getTotalDegree() <= 3 &&
    // RDKit✔️❌:         bondToTry->getEndAtom()->getTotalDegree() >= 2 &&
    // RDKit✔️❌:         bondToTry->getEndAtom()->getTotalDegree() <= 3) {
    // RDKit✔️❌:       anyBondPassesDegreeCheck = true;
    // RDKit✔️❌:       break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (!anyBondPassesDegreeCheck) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // defer cache update on the whole mol unless we actually have bonds to try
    // RDKit✔️❌:   // we need to do an update on the whole mol and not just incident atoms
    // RDKit✔️❌:   // because we need to calculate hybridization, which is non-local
    // RDKit✔️❌:   bool needsUpdate =
    // RDKit✔️❌:       mol.needsUpdatePropertyCache() ||
    // RDKit✔️❌:       std::any_of(mol.atoms().begin(), mol.atoms().end(), [](const auto atom) {
    // RDKit✔️❌:         return atom->getAtomicNum() != 0 &&
    // RDKit✔️❌:                atom->getHybridization() == Atom::HybridizationType::UNSPECIFIED;
    // RDKit✔️❌:       });
    // RDKit✔️❌:   if (needsUpdate) {
    // RDKit✔️❌:     mol.updatePropertyCache(false);
    // RDKit✔️❌:     MolOps::setConjugation(mol);
    // RDKit✔️❌:     MolOps::setHybridization(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto bondToTry : bondsToTry) {
    // RDKit✔️❌:     if (bondToTry->getBondType() != Bond::SINGLE ||
    // RDKit✔️❌:         bondToTry->getStereo() == Bond::BondStereo::STEREOANY ||
    // RDKit✔️❌:         // before, we checked only on totalDegree = 2 or 3,
    // RDKit✔️❌:         // but this causes false positives for something like a chiral sulfoxide
    // RDKit✔️❌:         // since the S is tetrahedral (sp3) but has only 3 substituents.
    // RDKit✔️❌:         // the hybridization code relies on totalDegree,
    // RDKit✔️❌:         // but modified to include and making sure to include conjugation
    // RDKit✔️❌:         // so while this is more expensive per molecule, it is closer to intent
    // RDKit✔️❌:         bondToTry->getBeginAtom()->getHybridization() != Atom::SP2 ||
    // RDKit✔️❌:         bondToTry->getEndAtom()->getHybridization() != Atom::SP2) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     DetectAtropisomerChiralityOneBond(bondToTry, mol, conf);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: void detectAtropisomerChirality
    // END RDKIT CPP BODY: detect_atropisomer_chirality

    // Performance note: this reproduces the source behavior, but computes
    // valence for the whole molecule before the degree precheck. RDKit lazily
    // updates only candidate endpoint atom caches for that precheck.
    let _ = params;
    let mut molecule = molecule;
    let candidates = atrop_bonds_to_try(&molecule)?;
    if candidates.is_empty() {
        return Ok(molecule);
    }

    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| {
            SdfReadError::Parse("atropisomer topology bond atom index out of range".into())
        })?;
    let valence =
        crate::assign_valence_with_options(&molecule, crate::ValenceModel::RdkitLike, false)
            .map_err(|err| SdfReadError::Parse(err.to_string()))?;

    let mut degree_check_passed = false;
    for candidate_id in &candidates {
        let Some(candidate) = molecule.bonds().get(candidate_id.index()) else {
            continue;
        };
        if candidate.order() == BondOrder::Single
            && candidate.stereo() != BondStereo::Any
            && (2..=3).contains(&atrop_total_degree(
                &molecule,
                &adjacency,
                &valence,
                candidate.begin(),
            )?)
            && (2..=3).contains(&atrop_total_degree(
                &molecule,
                &adjacency,
                &valence,
                candidate.end(),
            )?)
        {
            degree_check_passed = true;
            break;
        }
    }
    if !degree_check_passed {
        return Ok(molecule);
    }

    if molecule.derived_cache().valence.is_none()
        || molecule.atoms().iter().any(|atom| {
            atom.atomic_number() != 0 && atom.hybridization() == crate::Hybridization::Unspecified
        })
    {
        molecule.derived_cache_mut().valence = Some(valence);
        let conjugation = crate::operations::ops::sanitize_conjugation_assignment(
            crate::read_parts::MoleculeReadParts::from_molecule(&molecule),
        )
        .map_err(|err| SdfReadError::Parse(err.to_string()))?;
        for (bond, is_conjugated) in molecule
            .topology_block_mut()
            .bonds
            .iter_mut()
            .zip(conjugation)
        {
            bond.set_conjugated(is_conjugated);
        }
        let hybridization = crate::operations::ops::sanitize_hybridization_assignment(
            crate::read_parts::MoleculeReadParts::from_molecule(&molecule),
        )
        .map_err(|err| SdfReadError::Parse(err.to_string()))?;
        for (atom, hybridization) in molecule
            .topology_block_mut()
            .atoms
            .iter_mut()
            .zip(hybridization)
        {
            atom.set_hybridization(hybridization);
        }
    }

    let conformer = molecule.conformers_3d().first();
    let mut assignments = Vec::new();
    for candidate_id in candidates {
        let Some(candidate) = molecule.bonds().get(candidate_id.index()) else {
            continue;
        };
        if candidate.order() != BondOrder::Single
            || candidate.stereo() == BondStereo::Any
            || molecule.atoms()[candidate.begin().index()].hybridization()
                != crate::Hybridization::Sp2
            || molecule.atoms()[candidate.end().index()].hybridization()
                != crate::Hybridization::Sp2
        {
            continue;
        }
        if let Some(stereo) =
            detect_atropisomer_chirality_one_bond(&molecule, candidate_id, conformer)?
        {
            assignments.push((candidate_id, stereo));
        }
    }
    if !assignments.is_empty() {
        for (bond_id, stereo) in assignments {
            if let Some(bond) = molecule.topology_block_mut().bonds.get_mut(bond_id.index()) {
                bond.set_stereo(stereo);
            }
        }
        molecule.derived_cache_mut().invalidate(
            crate::DerivedState::STEREO
                | crate::DerivedState::DRAWING
                | crate::DerivedState::FINGERPRINT,
        );
    }
    Ok(molecule)
}

fn atrop_bonds_to_try(molecule: &Molecule) -> Result<Vec<BondId>, SdfReadError> {
    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| {
            SdfReadError::Parse("atropisomer topology bond atom index out of range".into())
        })?;
    let mut candidates = BTreeSet::new();
    for bond in molecule.bonds() {
        if !atrop_can_have_direction(bond.order())
            || !matches!(
                bond.direction(),
                BondDirection::BeginDash | BondDirection::BeginWedge
            )
        {
            continue;
        }
        for neighbor in adjacency.neighbors_of(bond.begin().index()) {
            if neighbor.bond == bond.id() {
                continue;
            }
            candidates.insert(neighbor.bond);
        }
    }
    Ok(candidates.into_iter().collect())
}

fn atrop_total_degree(
    molecule: &Molecule,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, SdfReadError> {
    let atom_data = molecule.atom(atom).ok_or_else(|| {
        SdfReadError::Parse("atropisomer total-degree atom index out of range".into())
    })?;
    let graph_degree = i32::try_from(adjacency.neighbors_of(atom.index()).len())
        .map_err(|_| SdfReadError::Parse("atropisomer atom degree out of range".into()))?;
    Ok(graph_degree
        + i32::from(atom_data.explicit_hydrogens())
        + valence
            .implicit_hydrogens
            .get(atom.index())
            .copied()
            .unwrap_or(0)
            .max(0))
}

fn atrop_can_have_direction(order: BondOrder) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Bond.h :: bool canHaveDirection
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Bond.h :: bool canHaveDirection
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

fn detect_atropisomer_chirality_one_bond(
    molecule: &Molecule,
    bond_id: BondId,
    conformer: Option<&Conformer3D>,
) -> Result<Option<BondStereo>, SdfReadError> {
    // BEGIN RDKIT CPP BODY: detect_atropisomer_chirality_one_bond
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: void DetectAtropisomerChiralityOneBond
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
    // RDKit✔️✔️:   RDGeom::Point3D xAxis, yAxis, zAxis;
    // RDKit✔️✔️:   if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D bondVecs[2];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit✔️✔️:     if (!conf->is3D()) {
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
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: void DetectAtropisomerChiralityOneBond
    // END RDKIT CPP BODY: detect_atropisomer_chirality_one_bond

    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let Some(candidate) = molecule.bonds().get(bond_id.index()) else {
        return Ok(None);
    };
    let Some(nbr0) = atrop_neighbor_bonds(molecule, candidate.begin(), bond_id)? else {
        return Ok(None);
    };
    let Some(nbr1) = atrop_neighbor_bonds(molecule, candidate.end(), bond_id)? else {
        return Ok(None);
    };
    for end_bond in nbr0.iter().chain(nbr1.iter()) {
        let Some(bond) = molecule.bonds().get(end_bond.index()) else {
            return Ok(None);
        };
        if bond.direction() == BondDirection::Unknown {
            return Ok(None);
        }
    }

    let Some(conformer) = conformer else {
        let (ok0, dir0) = atrop_end_wedge_direction(molecule, &nbr0);
        if !ok0 {
            return Ok(None);
        }
        let (ok1, dir1) = atrop_end_wedge_direction(molecule, &nbr1);
        if !ok1 || dir0 == dir1 {
            return Ok(None);
        }
        return Ok(
            if dir0 == BondDirection::BeginWedge || dir1 == BondDirection::BeginDash {
                Some(BondStereo::AtropCcw)
            } else if dir0 == BondDirection::BeginDash || dir1 == BondDirection::BeginWedge {
                Some(BondStereo::AtropCw)
            } else {
                None
            },
        );
    };

    if conformer.coordinates().len() < molecule.num_atoms() {
        return Err(SdfReadError::Parse(
            "atropisomer conformer coordinate count is smaller than atom count".into(),
        ));
    }
    let Some((_x_axis, y_axis, z_axis)) =
        atrop_bond_frame_of_reference(molecule, bond_id, conformer)
    else {
        return Ok(None);
    };
    let mut bond_vecs = [[0.0; 3]; 2];
    for (bond_atom_index, (focus_atom, nbr_bonds)) in
        [(candidate.begin(), &nbr0), (candidate.end(), &nbr1)]
            .into_iter()
            .enumerate()
    {
        if !conformer.is_3d() {
            let (ok, dir) = atrop_end_wedge_direction(molecule, nbr_bonds);
            if !ok {
                return Ok(None);
            }
            let Some(mut bond_vec) = atrop_projected_end_vector(
                molecule, conformer, focus_atom, nbr_bonds, y_axis, z_axis, true,
            ) else {
                return Ok(None);
            };
            if dir == BondDirection::BeginWedge {
                bond_vec[1] *= 0.707;
                bond_vec[2] = bond_vec[1].abs();
            } else if dir == BondDirection::BeginDash {
                bond_vec[1] *= 0.707;
                bond_vec[2] = -bond_vec[1].abs();
            }
            bond_vecs[bond_atom_index] = bond_vec;
        } else {
            let Some(bond_vec) = atrop_projected_end_vector(
                molecule, conformer, focus_atom, nbr_bonds, y_axis, z_axis, false,
            ) else {
                return Ok(None);
            };
            if vec3_len(bond_vec) < REALLY_SMALL_BOND_LEN {
                return Ok(None);
            }
            bond_vecs[bond_atom_index] = bond_vec;
        }
    }

    let cross_product = vec3_cross(bond_vecs[1], bond_vecs[0]);
    Ok(if cross_product[0] > REALLY_SMALL_BOND_LEN {
        Some(BondStereo::AtropCcw)
    } else if cross_product[0] < -REALLY_SMALL_BOND_LEN {
        Some(BondStereo::AtropCw)
    } else {
        None
    })
}

fn atrop_neighbor_bonds(
    molecule: &Molecule,
    focus_atom: AtomId,
    atrop_bond: BondId,
) -> Result<Option<Vec<BondId>>, SdfReadError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: bool getAtropisomerAtomsAndBonds
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: bool getAtropisomerAtomsAndBonds
    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| {
            SdfReadError::Parse("atropisomer topology bond atom index out of range".into())
        })?;
    let mut bonds: Vec<BondId> = adjacency
        .neighbors_of(focus_atom.index())
        .iter()
        .map(|neighbor| neighbor.bond)
        .filter(|bond| *bond != atrop_bond)
        .collect();
    if bonds.is_empty() {
        return Ok(None);
    }
    if bonds.len() == 2 {
        let other0 = atrop_other_atom(molecule, bonds[0], focus_atom)?;
        let other1 = atrop_other_atom(molecule, bonds[1], focus_atom)?;
        if other1.index() < other0.index() {
            bonds.swap(0, 1);
        }
    }
    Ok(Some(bonds))
}

fn atrop_other_atom(
    molecule: &Molecule,
    bond_id: BondId,
    atom: AtomId,
) -> Result<AtomId, SdfReadError> {
    let bond = molecule.bonds().get(bond_id.index()).ok_or_else(|| {
        SdfReadError::Parse("atropisomer neighbor bond index out of range".into())
    })?;
    if bond.begin() == atom {
        Ok(bond.end())
    } else if bond.end() == atom {
        Ok(bond.begin())
    } else {
        Err(SdfReadError::Parse(
            "atropisomer neighbor bond does not include focus atom".into(),
        ))
    }
}

fn atrop_end_wedge_direction(molecule: &Molecule, nbr_bonds: &[BondId]) -> (bool, BondDirection) {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: std::pair<bool, Bond::BondDir> getBondDir
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
    // RDKit✔️✔️:   if (bond1Dir != Bond::NONE && bond2Dir != Bond::NONE &&
    // RDKit✔️✔️:       bond1Dir == bond2Dir) {
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: std::pair<bool, Bond::BondDir> getBondDir
    let Some(bond0) = molecule.bonds().get(nbr_bonds[0].index()) else {
        return (false, BondDirection::None);
    };
    let bond1 = if nbr_bonds.len() == 2 {
        molecule.bonds().get(nbr_bonds[1].index())
    } else {
        None
    };
    let dir0 = match bond0.direction() {
        BondDirection::BeginWedge | BondDirection::BeginDash => bond0.direction(),
        _ => BondDirection::None,
    };
    let dir1 = match bond1.map(crate::Bond::direction) {
        Some(BondDirection::BeginWedge | BondDirection::BeginDash) => bond1
            .map(crate::Bond::direction)
            .unwrap_or(BondDirection::None),
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

fn atrop_bond_frame_of_reference(
    molecule: &Molecule,
    bond_id: BondId,
    conformer: &Conformer3D,
) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: bool getBondFrameOfReference
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
    // RDKit✔️✔️:   if (fabs(xAxis.x) > REALLY_SMALL_BOND_LEN ||
    // RDKit✔️✔️:       fabs(xAxis.y) > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(0, 0, 1);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(1, 0, 0);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   yAxis = zAxis.crossProduct(xAxis);
    // RDKit✔️✔️:   zAxis = xAxis.crossProduct(yAxis);
    // RDKit✔️✔️:   yAxis.normalize();
    // RDKit✔️✔️:   zAxis.normalize();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: bool getBondFrameOfReference
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let bond = molecule.bonds().get(bond_id.index())?;
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
        return Some((
            x_axis,
            [-x_axis[1] / y_len, x_axis[0] / y_len, 0.0],
            [0.0, 0.0, 1.0],
        ));
    }
    let mut z_axis =
        if x_axis[0].abs() > REALLY_SMALL_BOND_LEN || x_axis[1].abs() > REALLY_SMALL_BOND_LEN {
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

fn atrop_projected_end_vector(
    molecule: &Molecule,
    conformer: &Conformer3D,
    focus_atom: AtomId,
    nbr_bonds: &[BondId],
    y_axis: [f64; 3],
    z_axis: [f64; 3],
    normalize: bool,
) -> Option<[f64; 3]> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: bool getAtropIsomerEndVect
    // RDKit✔️✔️:   bondVec = conf->getAtomPos(atomAndBondVec.second[0]
    // RDKit✔️✔️:                                  ->getOtherAtom(atomAndBondVec.first)
    // RDKit✔️✔️:                                  ->getIdx()) -
    // RDKit✔️✔️:             conf->getAtomPos(
    // RDKit✔️✔️:                 atomAndBondVec.first->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec = RDGeom::Point3D(0.0, bondVec.dotProduct(yAxis),
    // RDKit✔️✔️:                             bondVec.dotProduct(zAxis));
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atropisomers.cpp :: bool getAtropIsomerEndVect
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let other0 = atrop_other_atom(molecule, nbr_bonds[0], focus_atom).ok()?;
    let mut bond_vec = vec3_sub(
        conformer.coordinates()[other0.index()],
        conformer.coordinates()[focus_atom.index()],
    );
    bond_vec = [0.0, vec3_dot(bond_vec, y_axis), vec3_dot(bond_vec, z_axis)];
    if nbr_bonds.len() == 2 {
        let other1 = atrop_other_atom(molecule, nbr_bonds[1], focus_atom).ok()?;
        let mut other_vec = vec3_sub(
            conformer.coordinates()[other1.index()],
            conformer.coordinates()[focus_atom.index()],
        );
        other_vec = [
            0.0,
            vec3_dot(other_vec, y_axis),
            vec3_dot(other_vec, z_axis),
        ];
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

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec3_len(v: [f64; 3]) -> f64 {
    vec3_dot(v, v).sqrt()
}

fn clear_single_bond_dir_flags(
    molecule: Molecule,
    params: SdfReadParams,
) -> Result<Molecule, SdfReadError> {
    let _ = params;
    clear_single_bond_dir_flags_with_mode(molecule, false)
}

fn clear_single_bond_dir_flags_with_mode(
    molecule: Molecule,
    only_wedge_flags: bool,
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

    let mut molecule = molecule;
    let topology = molecule.topology_block_mut();
    for bond in &mut topology.bonds {
        if bond.order() == BondOrder::Single {
            if bond.direction() == crate::BondDirection::Unknown {
                bond.set_unknown_stereo(true);
            }
            if !only_wedge_flags
                || !matches!(
                    bond.direction(),
                    crate::BondDirection::EndDownRight | crate::BondDirection::EndUpRight
                )
            {
                bond.set_direction(crate::BondDirection::None);
            }
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
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   setDoubleBondNeighborDirections(mol, &conf);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Chirality.cpp :: void detectBondStereochemistry
    // END RDKIT CPP BODY: detect_bond_stereochemistry

    let _ = params;
    let Some(conf_id) = molecule.conformers_3d().first().map(|conf| conf.id()) else {
        return Ok(molecule);
    };
    let mut molecule = molecule;
    // The shared stereo implementation mirrors RDKit's
    // setDoubleBondNeighborDirections data flow: one adjacency pass, per-bond
    // direction flags, neighbor vectors, ordered double-bond processing, and
    // recursive follow-up propagation.
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
        .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
    let molecule = molecule.sanitize().map_err(molecule_operation_error)?;
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
        .without_hydrogens_with_sanitize(true)
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
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    // RDKit❗❌:     MolOps::assignStereochemistry(*res, true, true, true);
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp :: void finishMolProcessing
    let _ = params;
    let mut molecule = molecule;
    crate::smiles::assign_stereochemistry_cleanup_subset(&mut molecule, true).map_err(|err| {
        SdfReadError::Parse(format!(
            "post-parse stereochemistry assignment failed: {err}"
        ))
    })?;
    crate::smiles::assign_double_bond_stereo_after_smiles_parse(&mut molecule, true).map_err(
        |err| {
            SdfReadError::Parse(format!(
                "post-parse double-bond stereochemistry assignment failed: {err}"
            ))
        },
    )?;
    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn process_mol_props_like_finish(molecule: &mut Molecule) -> Result<(), SdfReadError> {
        let explicit_valence = calculate_explicit_valence_before_mol_props(molecule)?;
        process_mol_props(molecule, SdfReadParams::default(), &explicit_valence)
    }

    #[test]
    fn process_mol_props_converts_mol_subst_count_queries_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let atom_zero = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_prop("molSubstCount", "-1"),
        );
        let atom_as_drawn = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_prop("molSubstCount", "-2"),
        );
        let atom_six_or_more = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_prop("molSubstCount", "6"),
        );
        builder
            .add_bond(BondSpec::new(atom_zero, atom_as_drawn, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                atom_as_drawn,
                atom_six_or_more,
                BondOrder::Single,
            ))
            .unwrap();
        let mut molecule = builder.build().unwrap();

        process_mol_props_like_finish(&mut molecule).unwrap();

        assert_eq!(
            molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(0)))
        );
        assert_eq!(
            molecule.atoms()[1].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(2)))
        );
        assert_eq!(
            molecule.atoms()[2].query(),
            Some(&QueryNode::not(QueryNode::predicate(
                AtomQueryPredicate::ExplicitDegreeLessEqual(5)
            )))
        );
    }

    #[test]
    fn process_mol_props_applies_mol_tot_valence_with_explicit_h_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1)
                .with_prop("molTotValence", "4"),
        );
        let oxygen = builder.add_atom(AtomSpec::new(Element::O).with_no_implicit(true));
        let nitrogen = builder.add_atom(
            AtomSpec::new(Element::N)
                .with_no_implicit(true)
                .with_prop("molTotValence", "15"),
        );
        builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single))
            .unwrap();
        let mut molecule = builder.build().unwrap();

        process_mol_props_like_finish(&mut molecule).unwrap();

        assert!(molecule.atoms()[0].no_implicit());
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 2);
        assert_eq!(molecule.atoms()[0].prop("molTotValence"), None);
        assert!(molecule.atoms()[2].no_implicit());
        assert_eq!(molecule.atoms()[2].explicit_hydrogens(), 0);
        assert_eq!(molecule.atoms()[2].id(), nitrogen);
    }

    #[test]
    fn process_mol_props_skips_zbo_h_total_valence_update_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(1)
                .with_prop("molTotValence", "4")
                .with_prop("_ZBO_H", "1"),
        );
        let mut molecule = builder.build().unwrap();

        process_mol_props_like_finish(&mut molecule).unwrap();

        assert!(!molecule.atoms()[0].no_implicit());
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(molecule.atoms()[0].prop("molTotValence"), None);
        assert_eq!(molecule.atoms()[0].prop("_ZBO_H"), Some("1"));
    }

    #[test]
    fn process_mol_props_ignores_molecule_data_fields_like_rdkit() {
        let mut builder = MoleculeBuilder::new()
            .with_property("_MolFileChiralFlag", "1")
            .with_property("molTotValence", "4");
        builder.add_atom(AtomSpec::new(Element::C));
        let mut molecule = builder.build().unwrap();

        process_mol_props_like_finish(&mut molecule).unwrap();

        assert_eq!(molecule.prop("_MolFileChiralFlag"), Some("1"));
        assert_eq!(molecule.prop("molTotValence"), Some("4"));
        assert!(!molecule.atoms()[0].no_implicit());
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 0);
    }

    fn single_bond_dir_fixture(direction: BondDirection, order: BondOrder) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(begin, end, order).with_direction(direction))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn clear_single_bond_dir_flags_clears_wedge_like_rdkit() {
        let molecule = single_bond_dir_fixture(BondDirection::EndUpRight, BondOrder::Single);

        let molecule = clear_single_bond_dir_flags(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert!(!molecule.bonds()[0].unknown_stereo());
    }

    #[test]
    fn clear_single_bond_dir_flags_clears_dash_like_rdkit() {
        let molecule = single_bond_dir_fixture(BondDirection::EndDownRight, BondOrder::Single);

        let molecule = clear_single_bond_dir_flags(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert!(!molecule.bonds()[0].unknown_stereo());
    }

    #[test]
    fn clear_single_bond_dir_flags_retains_non_single_bond_like_rdkit() {
        let molecule = single_bond_dir_fixture(BondDirection::EitherDouble, BondOrder::Double);

        let molecule = clear_single_bond_dir_flags(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EitherDouble);
        assert!(!molecule.bonds()[0].unknown_stereo());
    }

    #[test]
    fn clear_single_bond_dir_flags_only_wedge_false_records_unknown_like_rdkit() {
        let molecule = single_bond_dir_fixture(BondDirection::Unknown, BondOrder::Single);

        let molecule = clear_single_bond_dir_flags_with_mode(molecule, false).unwrap();

        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert!(molecule.bonds()[0].unknown_stereo());
    }

    #[test]
    fn clear_single_bond_dir_flags_only_wedge_true_preserves_wedge_like_rdkit() {
        let molecule = single_bond_dir_fixture(BondDirection::EndUpRight, BondOrder::Single);

        let molecule = clear_single_bond_dir_flags_with_mode(molecule, true).unwrap();

        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
        assert!(!molecule.bonds()[0].unknown_stereo());
    }

    #[test]
    fn clear_single_bond_dir_flags_runs_after_atom_stereo_perception_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, fluorine, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(
            molecule,
            true,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    }

    fn wedged_hydrogen_chiral_center_fixture() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn finish_mol_processing_perceives_wedged_hydrogen_before_remove_hs_like_rdkit() {
        let molecule = wedged_hydrogen_chiral_center_fixture();

        let molecule = finish_mol_processing(
            molecule,
            true,
            SdfReadParams {
                sanitize: false,
                remove_hs: true,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(molecule.num_atoms(), 5);
        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    }

    #[test]
    fn finish_mol_processing_assigns_3d_chirality_before_sanitize_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        for ligand in [hydrogen, fluorine, chlorine, bromine] {
            builder
                .add_bond(BondSpec::new(center, ligand, BondOrder::Single))
                .unwrap();
        }
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-1.0, -1.0, -1.0],
                ],
                true,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(
            molecule,
            false,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert!(
            molecule.atoms()[0]
                .prop("_NonExplicit3DChirality")
                .is_some()
        );
    }

    #[test]
    fn finish_mol_processing_preserves_unmodeled_bookmark_absence_like_rdkit_clearing() {
        let mut builder = MoleculeBuilder::new().with_property("user_prop", "kept");
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_prop("atom_prop", "kept"));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(
                BondSpec::new(carbon, oxygen, BondOrder::Single).with_prop("bond_prop", "kept"),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(
            molecule,
            false,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(
            molecule
                .properties()
                .props()
                .get("user_prop")
                .map(String::as_str),
            Some("kept")
        );
        assert_eq!(molecule.atoms()[0].prop("atom_prop"), Some("kept"));
        assert_eq!(molecule.bonds()[0].prop("bond_prop"), Some("kept"));
        assert!(
            !molecule
                .properties()
                .props()
                .keys()
                .any(|key| key.contains("Bookmark") || key.contains("bookmark"))
        );
        assert!(molecule.atoms().iter().all(|atom| {
            atom.props()
                .keys()
                .all(|key| !key.contains("Bookmark") && !key.contains("bookmark"))
        }));
        assert!(molecule.bonds().iter().all(|bond| {
            bond.props()
                .keys()
                .all(|key| !key.contains("Bookmark") && !key.contains("bookmark"))
        }));
    }

    #[test]
    fn finish_mol_processing_runs_process_mol_props_before_bond_dir_chirality_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_prop("molTotValence", "4"),
        );
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, fluorine, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(
            molecule,
            true,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(molecule.atoms()[0].prop("molTotValence"), None);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    }

    #[test]
    fn finish_mol_processing_uses_default_conformer_is_3d_for_chirality_branch_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        for ligand in [hydrogen, fluorine, chlorine, bromine] {
            builder
                .add_bond(BondSpec::new(center, ligand, BondOrder::Single))
                .unwrap();
        }
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-1.0, -1.0, -1.0],
                ],
                false,
            ))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                1,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-1.0, -1.0, -1.0],
                ],
                true,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(
            molecule,
            false,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(
            molecule.source_coordinate_dim(),
            Some(crate::CoordinateDimension::ThreeD)
        );
        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert!(
            molecule.atoms()[0]
                .prop("_NonExplicit3DChirality")
                .is_none()
        );
    }

    fn difluoroethene_with_optional_hydrogen(use_hydrogen: bool, is_3d: bool) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let left = builder.add_atom(AtomSpec::new(if use_hydrogen {
            Element::H
        } else {
            Element::from_atomic_number(9).unwrap()
        }));
        let right = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        builder
            .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c0, left, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, right, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 1.0, if is_3d { 0.2 } else { 0.0 }],
                    [1.0, -1.0, if is_3d { -0.2 } else { 0.0 }],
                ],
                is_3d,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn detect_bond_stereochemistry_assigns_2d_neighbor_dirs_like_rdkit() {
        let molecule = difluoroethene_with_optional_hydrogen(false, false);

        let molecule = detect_bond_stereochemistry(molecule, SdfReadParams::default()).unwrap();

        assert!(matches!(
            molecule.bonds()[1].direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ));
        assert!(matches!(
            molecule.bonds()[2].direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ));
        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
    }

    #[test]
    fn detect_bond_stereochemistry_assigns_3d_neighbor_dirs_like_rdkit() {
        let molecule = difluoroethene_with_optional_hydrogen(false, true);

        let molecule = detect_bond_stereochemistry(molecule, SdfReadParams::default()).unwrap();

        assert!(matches!(
            molecule.bonds()[1].direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ));
        assert!(matches!(
            molecule.bonds()[2].direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ));
    }

    #[test]
    fn detect_bond_stereochemistry_marks_crossed_double_bond_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        builder
            .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(f, c0, BondOrder::Single)
                    .with_direction(BondDirection::Unknown)
                    .with_unknown_stereo(true),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = detect_bond_stereochemistry(molecule, SdfReadParams::default()).unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Any);
        assert_eq!(molecule.bonds()[0].stereo_atoms(), None);
        assert_eq!(molecule.bonds()[1].direction(), BondDirection::Unknown);
    }

    #[test]
    fn detect_bond_stereochemistry_sets_stereo_atoms_after_direction_assignment_like_rdkit() {
        let mut molecule = difluoroethene_with_optional_hydrogen(false, false);

        molecule = detect_bond_stereochemistry(molecule, SdfReadParams::default()).unwrap();
        crate::smiles::assign_double_bond_stereo_from_directions(&mut molecule).unwrap();

        assert!(matches!(
            molecule.bonds()[0].stereo(),
            BondStereo::Cis | BondStereo::Trans
        ));
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn detect_bond_stereochemistry_uses_imine_hydrogen_before_removal_like_rdkit() {
        let mut molecule = difluoroethene_with_optional_hydrogen(true, false);

        molecule = detect_bond_stereochemistry(molecule, SdfReadParams::default()).unwrap();
        crate::smiles::assign_double_bond_stereo_from_directions(&mut molecule).unwrap();

        assert_eq!(molecule.num_atoms(), 4);
        assert!(matches!(
            molecule.bonds()[0].stereo(),
            BondStereo::Cis | BondStereo::Trans
        ));
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );

        let molecule = molecule.without_hydrogens().unwrap();
        assert_eq!(molecule.num_atoms(), 4);
    }

    fn sanitize_remove_hs_params() -> SdfReadParams {
        SdfReadParams {
            sanitize: true,
            remove_hs: true,
            ..Default::default()
        }
    }

    fn sanitize_keep_hs_params() -> SdfReadParams {
        SdfReadParams {
            sanitize: true,
            remove_hs: false,
            ..Default::default()
        }
    }

    fn unsanitized_params() -> SdfReadParams {
        SdfReadParams {
            sanitize: false,
            remove_hs: true,
            ..Default::default()
        }
    }

    fn v2000_atom_line_at(symbol: &str, x: f64, y: f64, z: f64) -> String {
        format!(
            "{x:>10.4}{y:>10.4}{z:>10.4} {symbol:<3}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}",
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        )
    }

    fn v2000_bond_line(begin: u32, end: u32, bond_type: u32, stereo: u32, topology: u32) -> String {
        format!(
            "{begin:>3}{end:>3}{bond_type:>3}{stereo:>3}{:>3}{topology:>3}",
            0
        )
    }

    fn nitro_group_fixture() -> (Molecule, AtomId, AtomId, AtomId) {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let oxygen_single = builder.add_atom(AtomSpec::new(Element::O));
        let oxygen_double = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(carbon, nitrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(nitrogen, oxygen_single, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(nitrogen, oxygen_double, BondOrder::Double))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [2.0, 1.0, 0.0],
                    [2.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        (
            builder.build().unwrap(),
            nitrogen,
            oxygen_single,
            oxygen_double,
        )
    }

    #[test]
    fn complete_mol_queries_completes_v2000_atom_query_and_clears_scan_flag_like_rdkit() {
        let input = format!(
            "v2000-rbcnt-as-drawn\n  COSMolKit          2D\ncomment\n  3  3  0  0  0  0            999 V2000\n{}\n{}\n{}\n{}\n{}\n{}\nM  RBC  1   1  -2\nM  END\n",
            v2000_atom_line_at("C", 0.0, 0.0, 0.0),
            v2000_atom_line_at("C", 1.0, 0.0, 0.0),
            v2000_atom_line_at("C", 0.5, 1.0, 0.0),
            v2000_bond_line(1, 2, 1, 0, 0),
            v2000_bond_line(2, 3, 1, 0, 0),
            v2000_bond_line(3, 1, 1, 0, 0),
        );

        let record = crate::io::molfile::read_mol_record_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.prop("_NeedsQueryScan"), None);
        assert_eq!(
            record.molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
    }

    #[test]
    fn complete_mol_queries_preserves_v2000_bond_query_while_completing_atom_query_like_rdkit() {
        let input = format!(
            "v2000-rbcnt-and-bond-query\n  COSMolKit          2D\ncomment\n  3  3  0  0  0  0            999 V2000\n{}\n{}\n{}\n{}\n{}\n{}\nM  RBC  1   1  -2\nM  END\n",
            v2000_atom_line_at("C", 0.0, 0.0, 0.0),
            v2000_atom_line_at("C", 1.0, 0.0, 0.0),
            v2000_atom_line_at("C", 0.5, 1.0, 0.0),
            v2000_bond_line(1, 2, 5, 0, 0),
            v2000_bond_line(2, 3, 1, 0, 0),
            v2000_bond_line(3, 1, 1, 0, 0),
        );

        let record = crate::io::molfile::read_mol_record_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.prop("_NeedsQueryScan"), None);
        assert_eq!(
            record.molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
        assert_eq!(
            record.molecule.bonds()[0].query(),
            Some(&QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])))
        );
    }

    #[test]
    fn complete_mol_queries_completes_v3000_atom_query_and_preserves_bond_query_like_rdkit() {
        let mut counts_line = "  0  0  0  0  0            999".to_string();
        while counts_line.len() < 34 {
            counts_line.push(' ');
        }
        counts_line.push_str("V3000");
        let input = format!(
            "\
v3000-rbcnt-and-bond-query
  COSMolKit

{counts_line}
M  V30 BEGIN CTAB
M  V30 COUNTS 3 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0.0 0 RBCNT=-2
M  V30 2 C 1.0 0.0 0.0 0
M  V30 3 C 0.5 1.0 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 5 1 2
M  V30 2 1 2 3
M  V30 3 1 3 1
M  V30 END BOND
M  V30 END CTAB
M  END
"
        );

        let record = crate::io::molfile::read_mol_record_from_str_with_params(
            &input,
            SdfReadParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(record.molecule.prop("_NeedsQueryScan"), None);
        assert_eq!(
            record.molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
        assert_eq!(
            record.molecule.bonds()[0].query(),
            Some(&QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])))
        );
    }

    #[test]
    fn molblock_sanitize_remove_hs_runs_cleanup_before_stereo_like_rdkit() {
        let (molecule, nitrogen, oxygen_single, _oxygen_double) = nitro_group_fixture();

        let molecule = finish_mol_processing(molecule, false, sanitize_remove_hs_params()).unwrap();

        assert_eq!(molecule.atoms()[nitrogen.index()].formal_charge(), 1);
        assert_eq!(molecule.atoms()[oxygen_single.index()].formal_charge(), -1);
        assert_eq!(molecule.bonds()[1].order(), BondOrder::Single);
        assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
    }

    #[test]
    fn molblock_sanitize_remove_hs_detects_double_bond_stereogenic_h_like_rdkit() {
        let molecule = difluoroethene_with_optional_hydrogen(true, false);

        let molecule = finish_mol_processing(molecule, false, sanitize_remove_hs_params()).unwrap();

        assert_eq!(molecule.num_atoms(), 4);
        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::E);
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn molblock_sanitize_remove_hs_preserves_wedged_hydrogen_chirality_like_rdkit() {
        let molecule = wedged_hydrogen_chiral_center_fixture();

        let molecule = finish_mol_processing(molecule, true, sanitize_remove_hs_params()).unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    }

    #[test]
    fn molblock_sanitize_remove_hs_handles_explicit_h_on_aromatic_nitrogen_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let nitrogen = builder.add_atom(
            AtomSpec::new(Element::N)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c2 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c3 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let c4 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        builder
            .add_bond(BondSpec::new(nitrogen, c1, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, c3, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c3, c4, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c4, nitrogen, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        builder
            .add_bond(BondSpec::new(nitrogen, hydrogen, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.5, 1.0, 0.0],
                    [0.5, 1.8, 0.0],
                    [-0.5, 1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(molecule, false, sanitize_remove_hs_params()).unwrap();

        assert_eq!(molecule.num_atoms(), 5);
        assert_eq!(molecule.atoms()[nitrogen.index()].explicit_hydrogens(), 1);
        assert!(molecule.atoms()[nitrogen.index()].is_aromatic());
    }

    #[test]
    fn molblock_sanitize_remove_hs_assigns_final_atom_and_bond_stereo_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        let alkene_left = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let alkene_right = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let alkene_f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let alkene_cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        builder
            .add_bond(
                BondSpec::new(center, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(alkene_left, alkene_right, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(alkene_left, alkene_f, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(alkene_right, alkene_cl, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [3.0, 0.0, 0.0],
                    [5.0, 0.0, 0.0],
                    [3.0, 1.0, 0.0],
                    [5.0, -1.0, 0.0],
                ],
                false,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let molecule = finish_mol_processing(molecule, true, sanitize_remove_hs_params()).unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.num_atoms(), 8);
        assert_eq!(molecule.bonds()[3].stereo(), BondStereo::E);
        assert_eq!(
            molecule.bonds()[3].stereo_atoms(),
            Some([AtomId::new(6), AtomId::new(7)])
        );
    }

    #[test]
    fn molblock_sanitize_keep_hs_preserves_explicit_hydrogen_like_rdkit() {
        let molecule = wedged_hydrogen_chiral_center_fixture();

        let molecule = finish_mol_processing(molecule, true, sanitize_keep_hs_params()).unwrap();

        assert_eq!(molecule.num_atoms(), 5);
        assert_eq!(molecule.atoms()[1].atomic_number(), 1);
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
    }

    #[test]
    fn molblock_sanitize_keep_hs_detects_double_bond_stereo_after_sanitize_like_rdkit() {
        let molecule = difluoroethene_with_optional_hydrogen(false, false);

        let molecule = finish_mol_processing(molecule, false, sanitize_keep_hs_params()).unwrap();

        assert_eq!(molecule.num_atoms(), 4);
        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::E);
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn molblock_sanitize_keep_hs_assigns_final_stereochemistry_like_rdkit() {
        let molecule = wedged_hydrogen_chiral_center_fixture();

        let molecule = finish_mol_processing(molecule, true, sanitize_keep_hs_params()).unwrap();

        assert_ne!(
            molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert_eq!(molecule.prop("_StereochemDone"), Some("1"));
    }

    #[test]
    fn molblock_sanitize_keep_hs_sets_property_cache_state_like_rdkit() {
        let molecule = difluoroethene_with_optional_hydrogen(false, false);

        let molecule = finish_mol_processing(molecule, false, sanitize_keep_hs_params()).unwrap();

        assert!(molecule.derived_cache().valence.is_some());
        assert!(molecule.derived_cache().rings.is_some());
        assert_eq!(molecule.prop("_StereochemDone"), Some("1"));
    }

    #[test]
    fn molblock_unsanitized_preserves_raw_valence_and_delays_cleanup_like_rdkit() {
        let (molecule, nitrogen, oxygen_single, oxygen_double) = nitro_group_fixture();

        let molecule = finish_mol_processing(molecule, false, unsanitized_params()).unwrap();

        assert_eq!(molecule.atoms()[nitrogen.index()].formal_charge(), 0);
        assert_eq!(molecule.atoms()[oxygen_single.index()].formal_charge(), 0);
        assert_eq!(molecule.bonds()[1].order(), BondOrder::Double);
        assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
        assert!(molecule.derived_cache().valence.is_none());
        assert!(molecule.derived_cache().rings.is_some());
        assert_eq!(molecule.prop("_StereochemDone"), None);

        let sanitized = molecule.sanitize().unwrap();
        assert_eq!(sanitized.atoms()[nitrogen.index()].formal_charge(), 1);
        assert_eq!(sanitized.atoms()[oxygen_single.index()].formal_charge(), -1);
        assert_eq!(sanitized.bonds()[1].order(), BondOrder::Single);
        assert_eq!(sanitized.bonds()[2].order(), BondOrder::Double);
        assert_eq!(sanitized.atoms()[oxygen_double.index()].formal_charge(), 0);
    }

    #[test]
    fn molblock_unsanitized_preserves_aromaticity_state_and_hydrogens_like_rdkit() {
        let molecule = wedged_hydrogen_chiral_center_fixture();

        let molecule = finish_mol_processing(molecule, true, unsanitized_params()).unwrap();

        assert_eq!(molecule.num_atoms(), 5);
        assert_eq!(molecule.atoms()[1].atomic_number(), 1);
        assert!(molecule.derived_cache().valence.is_none());
        assert!(molecule.derived_cache().rings.is_some());
        assert_eq!(molecule.prop("_StereochemDone"), None);
    }

    #[test]
    fn molblock_unsanitized_detects_bond_stereo_without_final_cip_like_rdkit() {
        let molecule = difluoroethene_with_optional_hydrogen(false, false);

        let molecule = finish_mol_processing(molecule, false, unsanitized_params()).unwrap();

        assert!(matches!(
            molecule.bonds()[1].direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ));
        assert!(matches!(
            molecule.bonds()[2].direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ));
        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
        assert_eq!(molecule.bonds()[0].stereo_atoms(), None);
        assert_eq!(molecule.prop("_StereochemDone"), None);
    }
}
