//! Molfile-specific detached postprocessing.

use cosmolkit_core::{
    AtropisomerConformer, RemoveHsParams, RingSearchParams, SanitizeOperations, SanitizeParams,
    StructureTagParams, ValenceModel, assign_chiral_tags_from_structure,
    assign_chiral_types_from_bond_dirs, assign_legacy_stereochemistry_with_query_state,
    assign_valence_for_topology, calculate_explicit_valence_for_topology,
    clear_single_bond_directions, detect_atropisomer_chirality, remove_hydrogens_with_query_state,
    sanitize_topology_with_query_state, set_double_bond_neighbor_directions, symmetrized_sssr,
};
use cosmolkit_model::{
    AdjacencyList, AtomId, AtomQueryPredicate, BondQueryPredicate, Conformer3D, CoordinateBlock,
    QueryAtom, QueryBond, QueryGraph, QueryNode, QueryStateRef, RecursiveStructureQuery,
    SubstanceGroup, SubstanceGroupId, TopologyBlock, TopologyMapping, remap_query_rows,
};
use cosmolkit_types::BondOrder;

use crate::sdf::{
    MolBlockRecord, QueryMolBlockRecord, parse_rdkit_int, query_from_concrete_atom_value,
};

/// Source options applied after Molfile syntax parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolPostParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub expand_attachment_points: bool,
}

impl Default for MolPostParams {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MolFileParserParams
        // RDKit✔️✔️:   bool sanitize = true;      /**< sanitize the molecule after building it */
        // RDKit✔️✔️:   bool removeHs = true;      /**< remove Hs after constructing the molecule */
        // RDKit✔️✔️:   bool expandAttachmentPoints =
        // RDKit✔️✔️:       false; /**< toggle conversion of attachment points into dummy atoms */
        // END RDKIT CPP FUNCTION
        Self {
            sanitize: true,
            remove_hs: true,
            expand_attachment_points: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MolPostError {
    #[error("attachment-point expansion is not representable at the current detached boundary")]
    AttachmentPointExpansion,
    #[error("Molfile postprocessing property is outside the detached model: {0}")]
    Representation(&'static str),
    #[error("Molfile postprocessing failed: {0}")]
    Processing(String),
}

fn data_values(group: &SubstanceGroup) -> &[String] {
    group
        .data()
        .map_or(group.data_fields(), |data| data.values.as_slice())
}

fn retain_substance_groups(
    groups: Vec<SubstanceGroup>,
    remove: &[bool],
) -> Result<Vec<SubstanceGroup>, MolPostError> {
    let map_len = groups
        .iter()
        .map(|group| group.id().index())
        .max()
        .map_or(0, |maximum| maximum + 1);
    let mut old_to_new = vec![None; map_len];
    let mut retained =
        Vec::with_capacity(groups.len() - remove.iter().filter(|flag| **flag).count());
    for (position, group) in groups.into_iter().enumerate() {
        if !remove[position] {
            old_to_new[group.id().index()] = Some(SubstanceGroupId::new(retained.len()));
            retained.push(group);
        }
    }
    for group in &mut retained {
        let new_id = old_to_new[group.id().index()].ok_or(MolPostError::Representation(
            "retained SGroup id was not mapped",
        ))?;
        group.set_id(new_id);
        if let Some(old_parent) = group.parent() {
            let new_parent = old_to_new
                .get(old_parent.index())
                .and_then(|mapped| *mapped)
                .ok_or(MolPostError::Representation(
                    "retained SGroup references a consumed parent",
                ))?;
            group.set_parent(new_parent);
        }
    }
    Ok(retained)
}

fn process_groups_on_topology(topology: &mut TopologyBlock) -> Result<(), MolPostError> {
    // BEGIN RDKIT CPP FUNCTION processSGroups
    // RDKit✔️✔️:   for (auto &sg : getSubstanceGroups(*mol)) {
    // RDKit✔️✔️:     if (sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit✔️✔️:       std::string field;
    // RDKit✔️✔️:       if (sg.getPropIfPresent("FIELDNAME", field)) {
    // RDKit✔️✔️:         if (field == "MRV_COORDINATE_BOND_TYPE") {
    // RDKit✔️✔️:           processMrvCoordinateBond(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "MRV_IMPLICIT_H") {
    // RDKit✔️✔️:           processMrvImplicitH(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "ZBO") {
    // RDKit✔️✔️:           processZBO(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "ZCH") {
    // RDKit✔️✔️:           processZCH(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         } else if (field == "HYD") {
    // RDKit✔️✔️:           processHYD(*mol, sg);
    // RDKit✔️✔️:           sgsToRemove.push_back(sgIdx);
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto it = sgsToRemove.rbegin(); it != sgsToRemove.rend(); ++it) {
    // RDKit✔️✔️:     sgs.erase(sgs.begin() + *it);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Behavior review: recognized DAT groups are processed in input order and
    // removed only after their source action. Invalid coordinate-bond targets
    // remain source-style no-ops. Checked Rust model widths are explicit errors.
    // Complexity review: one linear group pass plus membership-local work and
    // one retain pass matches the source's linear traversal/removal shape.
    let groups = topology.substance_groups.clone();
    let mut remove = vec![false; groups.len()];
    for (index, group) in groups.iter().enumerate() {
        let Some(data) = group.data() else { continue };
        let Some(field) = data.field_name.as_deref() else {
            continue;
        };
        match field {
            "MRV_COORDINATE_BOND_TYPE" => {
                if let Some(value) = data_values(group).first()
                    && let Ok(raw) = value.trim().parse::<usize>()
                    && let Some(bond) = raw
                        .checked_sub(1)
                        .and_then(|idx| topology.bonds.get_mut(idx))
                    && bond.order() == BondOrder::Unspecified
                {
                    bond.set_order(BondOrder::Dative);
                }
                remove[index] = true;
            }
            "MRV_IMPLICIT_H" => {
                for value in data_values(group) {
                    let Some(raw) = value.strip_prefix("IMPL_H") else {
                        continue;
                    };
                    let count = parse_rdkit_int(raw).unwrap_or(0);
                    let count = u8::try_from(count).map_err(|_| {
                        MolPostError::Representation("MRV_IMPLICIT_H count outside u8")
                    })?;
                    for atom_id in group.atoms() {
                        if atom_id.index() >= topology.atoms.len() {
                            continue;
                        }
                        let aromatic = topology.adjacency.neighbors_of(atom_id.index()).iter().any(
                            |neighbor| {
                                let bond = &topology.bonds[neighbor.bond.index()];
                                bond.is_aromatic() || bond.order() == BondOrder::Aromatic
                            },
                        );
                        if aromatic {
                            topology.atoms[atom_id.index()].set_explicit_hydrogens(count);
                        }
                    }
                }
                remove[index] = true;
            }
            "ZBO" => {
                for bond_id in group.bonds() {
                    if let Some(bond) = topology.bonds.get_mut(bond_id.index()) {
                        bond.set_order(BondOrder::Zero);
                    }
                }
                remove[index] = true;
            }
            "ZCH" | "HYD" => {
                for value in data_values(group) {
                    let values = value.trim().split(';').collect::<Vec<_>>();
                    if values.len() < group.atoms().len() {
                        continue;
                    }
                    for (atom_id, text) in group.atoms().iter().zip(values) {
                        let parsed = if text.is_empty() {
                            0
                        } else {
                            parse_rdkit_int(text).unwrap_or(0)
                        };
                        let Some(atom) = topology.atoms.get_mut(atom_id.index()) else {
                            continue;
                        };
                        if field == "ZCH" {
                            atom.set_formal_charge(i8::try_from(parsed).map_err(|_| {
                                MolPostError::Representation("ZCH charge outside i8")
                            })?);
                        } else {
                            atom.set_prop("_ZBO_H", "1")
                                .map_err(|error| MolPostError::Processing(error.to_string()))?;
                            atom.set_explicit_hydrogens(u8::try_from(parsed).map_err(|_| {
                                MolPostError::Representation("HYD count outside u8")
                            })?);
                        }
                    }
                }
                remove[index] = true;
            }
            _ => {}
        }
    }
    topology.substance_groups =
        retain_substance_groups(std::mem::take(&mut topology.substance_groups), &remove)?;
    Ok(())
}

fn process_atom_properties(
    topology: &mut TopologyBlock,
    mut query_atoms: Option<&mut [QueryAtom]>,
) -> Result<(), MolPostError> {
    // BEGIN RDKIT CPP FUNCTION ProcessMolProps
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
    // RDKit✔️✔️:         ival = atom->getDegree();
    // RDKit✔️✔️:       } else if (ival >= 6) {
    // RDKit✔️✔️:         gtQuery = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!gtQuery) {
    // RDKit✔️✔️:         atom->expandQuery(makeAtomExplicitDegreeQuery(ival));
    // RDKit✔️✔️:       } else {
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
    // RDKit✔️✔️:       if (ival == 15 || ival == -1) {
    // RDKit✔️✔️:         atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       } else if (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) > ival) {
    // RDKit✔️✔️:         atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         atom->setNumExplicitHs(ival - atom->getValence(Atom::ValenceType::EXPLICIT));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->clearProp(common_properties::molTotValence);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Behavior review: the optional query slice is present exactly after the
    // Molfile owner has promoted a record containing SUBST/SMARTSQ/query
    // syntax. SUBST expands the existing predicate in source order and marks
    // the carrier as a real source QueryAtom for downstream RemoveHs. The
    // total-valence branch observes the same current carrier topology.
    // Complexity review: one atom pass with O(1) predicate construction and
    // one valence calculation per atom matches the source traversal shape.
    for index in 0..topology.atoms.len() {
        let substitution = topology.atoms[index]
            .prop("molSubstCount")
            .and_then(|value| parse_rdkit_int(value).ok())
            .unwrap_or(0);
        if substitution != 0 {
            let atoms = query_atoms
                .as_deref_mut()
                .ok_or(MolPostError::Representation(
                    "molSubstCount requires query record promotion",
                ))?;
            let degree = if substitution == -1 {
                0
            } else if substitution == -2 {
                i64::try_from(topology.adjacency.neighbors_of(index).len())
                    .map_err(|_| MolPostError::Representation("molSubstCount degree outside i64"))?
            } else {
                i64::from(substitution)
            };
            let degree = u8::try_from(degree).map_err(|_| {
                MolPostError::Representation("molSubstCount query target outside u8")
            })?;
            let predicate = if substitution >= 6 {
                AtomQueryPredicate::ExplicitDegreeLessEqual(degree)
            } else {
                AtomQueryPredicate::ExplicitDegree(degree)
            };
            let current = atoms[index].predicate().clone();
            atoms[index].set_predicate(QueryNode::and(vec![
                current,
                QueryNode::predicate(predicate),
            ]));
            topology.atoms[index]
                .set_prop("_MolFileAtomQuery", "1")
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
            atoms[index]
                .atom_mut()
                .set_prop("_MolFileAtomQuery", "1")
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        }
        let value = topology.atoms[index]
            .prop("molTotValence")
            .and_then(|value| parse_rdkit_int(value).ok())
            .unwrap_or(0);
        if value != 0 && topology.atoms[index].prop("_ZBO_H").is_none() {
            let explicit =
                calculate_explicit_valence_for_topology(topology, AtomId::new(index), false, false)
                    .map_err(|error| MolPostError::Processing(error.to_string()))?;
            let hydrogens = if value == 15 || value == -1 || explicit > value {
                0
            } else {
                u8::try_from(value - explicit)
                    .map_err(|_| MolPostError::Representation("molTotValence H count outside u8"))?
            };
            topology.atoms[index].set_no_implicit(true);
            topology.atoms[index].set_explicit_hydrogens(hydrogens);
        }
        topology.atoms[index].clear_prop("molTotValence");
    }
    Ok(())
}

fn apply_stereo_and_sanitize(
    mut topology: TopologyBlock,
    mut coordinates: CoordinateBlock,
    mut properties: cosmolkit_model::MoleculeProperties,
    chirality_possible: bool,
    params: MolPostParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<
    (
        TopologyBlock,
        CoordinateBlock,
        cosmolkit_model::MoleculeProperties,
        TopologyMapping,
        Option<(Vec<QueryAtom>, Vec<QueryBond>)>,
    ),
    MolPostError,
> {
    // BEGIN RDKIT CPP FUNCTION finishMolProcessing applicable stereo/sanitize closure
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
    // RDKit✔️✔️:       MolOps::assignChiralTypesFromBondDirs(*res, conf.getId(),
    // RDKit✔️✔️:                                           replaceExistingTags);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       res->updatePropertyCache(false);
    // RDKit✔️✔️:       MolOps::assignChiralTypesFrom3D(*res, conf.getId(), true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Atropisomers::detectAtropisomerChirality(*res, &conf);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now that atom stereochem has been perceived, the wedging
    // RDKit✔️✔️:   // information is no longer needed, so we clear
    // RDKit✔️✔️:   // single bond dir flags:
    // RDKit✔️✔️:   MolOps::clearSingleBondDirFlags(*res);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (params.sanitize) {
    // RDKit✔️✔️:     if (params.removeHs) {
    // RDKit✔️✔️:       // Bond stereo detection must happen before H removal, or
    // RDKit✔️✔️:       // else we might be removing stereogenic H atoms in double
    // RDKit✔️✔️:       // bonds (e.g. imines). But before we run stereo detection,
    // RDKit✔️✔️:       // we need to run mol cleanup so don't have trouble with
    // RDKit✔️✔️:       // e.g. nitro groups. Sadly, this a;; means we will find
    // RDKit✔️✔️:       // run both cleanup and ring finding twice (a fast find
    // RDKit✔️✔️:       // rings in bond stereo detection, and another in
    // RDKit✔️✔️:       // sanitization's SSSR symmetrization).
    // RDKit✔️✔️:       unsigned int failedOp = 0;
    // RDKit✔️✔️:       MolOps::sanitizeMol(*res, failedOp, MolOps::SANITIZE_CLEANUP);
    // RDKit✔️✔️:       MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:       MolOps::removeHs(*res);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       MolOps::sanitizeMol(*res);
    // RDKit✔️✔️:       MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*res, true, true, true);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION finishMolProcessing applicable stereo/sanitize closure
    // Behavior review: source RWMol dynamic QueryAtom/QueryBond identity remains
    // available to every called owner. The detached adaptation carries the
    // validated typed query rows through sanitize, hydrogen removal and legacy
    // stereo, and remaps them with the same authoritative topology mapping.
    // Coordinate-driven double-bond directions retain source ordering. The
    // preceding bookmark, attachment, explicit-valence and ProcessMolProps
    // statements have distinct owners or the intentionally open attachment
    // gate and are not claimed by this helper.
    // Complexity review: each owner call is linear or owner-defined; this
    // orchestration introduces no repeated whole-graph clone beyond the owned
    // source-equivalent transform results.
    let original_atom_count = topology.atoms.len();
    let original_bond_count = topology.bonds.len();
    let mut mapping = TopologyMapping::identity(original_atom_count, original_bond_count);
    let first_3d = coordinates.conformers_3d.first();
    if let Some(conformer) = first_3d {
        if conformer.is_3d() {
            let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
            topology = assign_chiral_tags_from_structure(
                &topology,
                &coordinates,
                &valence,
                &StructureTagParams::default(),
            )
            .map_err(|error| MolPostError::Processing(error.to_string()))?
            .topology;
        } else if chirality_possible {
            assign_chiral_types_from_bond_dirs(&mut topology, conformer, true)
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        }
        let assignment =
            detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(conformer)))
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        for update in assignment.bond_updates {
            topology.bonds[update.bond.index()]
                .set_stereo(update.stereo)
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        }
    } else if let Some(conformer) = coordinates.conformers_2d.first() {
        if chirality_possible {
            let pseudo = Conformer3D::new(
                conformer.id(),
                conformer
                    .coordinates()
                    .iter()
                    .map(|xy| [xy[0], xy[1], 0.0])
                    .collect(),
                false,
            );
            assign_chiral_types_from_bond_dirs(&mut topology, &pseudo, true)
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        }
        let assignment =
            detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::TwoD(conformer)))
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        for update in assignment.bond_updates {
            topology.bonds[update.bond.index()]
                .set_stereo(update.stereo)
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
        }
    }
    topology = clear_single_bond_directions(topology, false)
        .map_err(|error| MolPostError::Processing(error.to_string()))?;
    if params.sanitize {
        if params.remove_hs {
            topology = sanitize_topology_with_query_state(
                &topology,
                &SanitizeParams {
                    operations: SanitizeOperations::CLEANUP,
                },
                query_state,
            )
            .map_err(|error| MolPostError::Processing(error.to_string()))?
            .topology;
            topology = detect_double_bond_stereochemistry(topology, &coordinates)?;
            let removed = remove_hydrogens_with_query_state(
                topology,
                coordinates,
                properties,
                &RemoveHsParams::default(),
                query_state,
            )
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
            topology = removed.topology;
            coordinates = removed.coordinates;
            properties = removed.properties;
            mapping = removed.mapping;
        } else {
            topology = sanitize_topology_with_query_state(
                &topology,
                &SanitizeParams::default(),
                query_state,
            )
            .map_err(|error| MolPostError::Processing(error.to_string()))?
            .topology;
            topology = detect_double_bond_stereochemistry(topology, &coordinates)?;
        }
        let remapped_query_rows = query_state
            .map(|state| remap_query_rows(state, &topology, &mapping))
            .transpose()
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
        let final_query_state = remapped_query_rows
            .as_ref()
            .map(|(atoms, bonds)| QueryStateRef::try_for_topology(atoms, bonds, &topology))
            .transpose()
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
        let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
        let rings = symmetrized_sssr(&topology, &RingSearchParams::default())
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
        topology = assign_legacy_stereochemistry_with_query_state(
            topology,
            &valence,
            &rings,
            final_query_state,
        )
        .map_err(|error| MolPostError::Processing(error.to_string()))?;
        // RDKit✔️✔️: mol.setProp(common_properties::_StereochemDone, 1, true);
        // Behavior review: molecule-level computed properties are carried by
        // the detached `MoleculeProperties` block, so this is the implementing
        // location for the source wrapper's final property write.
        // Complexity review: one ordered-map property update matches the
        // source computed-property write and adds no graph traversal.
        properties
            .set_computed_prop("_StereochemDone", "1")
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
    } else {
        topology = detect_double_bond_stereochemistry(topology, &coordinates)?;
    }
    let query_rows = query_state
        .map(|state| remap_query_rows(state, &topology, &mapping))
        .transpose()
        .map_err(|error| MolPostError::Processing(error.to_string()))?;
    Ok((topology, coordinates, properties, mapping, query_rows))
}

fn detect_double_bond_stereochemistry(
    topology: TopologyBlock,
    coordinates: &CoordinateBlock,
) -> Result<TopologyBlock, MolPostError> {
    // BEGIN RDKIT CPP FUNCTION detectBondStereochemistry
    // RDKit✔️✔️: void detectBondStereochemistry(ROMol &mol, int confId) {
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   setDoubleBondNeighborDirections(mol, &conf);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detectBondStereochemistry
    // Behavior review: the first stored conformer is used regardless of its
    // independent is_3d flag; XY is lifted with exact positive-zero Z only for
    // the detached geometry kernel. This prepares directions without assigning
    // final bond stereo, matching the source phase boundary.
    // Complexity review: ring perception plus the source-shaped direction
    // owner performs the same graph-scale work without an additional clone.
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default())
        .map_err(|error| MolPostError::Processing(error.to_string()))?;
    if let Some(conformer) = coordinates.conformers_3d.first() {
        return set_double_bond_neighbor_directions(topology, &rings, Some(conformer))
            .map_err(|error| MolPostError::Processing(error.to_string()));
    }
    let Some(conformer) = coordinates.conformers_2d.first() else {
        return Ok(topology);
    };
    let lifted = Conformer3D::new(
        conformer.id(),
        conformer
            .coordinates()
            .iter()
            .map(|xy| [xy[0], xy[1], 0.0])
            .collect(),
        false,
    );
    set_double_bond_neighbor_directions(topology, &rings, Some(&lifted))
        .map_err(|error| MolPostError::Processing(error.to_string()))
}

fn concrete_to_query(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: cosmolkit_model::MoleculeProperties,
) -> Result<QueryMolBlockRecord, MolPostError> {
    let atoms = topology
        .atoms
        .into_iter()
        .map(|atom| {
            let predicate = query_from_concrete_atom_value(&atom);
            QueryAtom::from_carrier_parts(atom, predicate)
        })
        .collect();
    let bonds = topology
        .bonds
        .into_iter()
        .map(|bond| {
            let predicate = if bond.order() == BondOrder::Unspecified {
                QueryNode::predicate(BondQueryPredicate::Any)
            } else {
                QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
            };
            QueryBond::from_carrier_parts(bond, predicate)
        })
        .collect();
    let mut props = properties.props().clone();
    if let Some(name) = properties.name() {
        props.insert("_Name".to_owned(), name.to_owned());
    }
    let query = QueryGraph::from_parts(
        atoms,
        bonds,
        props,
        coordinates.conformers_2d,
        coordinates.conformers_3d,
        topology.stereo_groups,
    )
    .map_err(|error| MolPostError::Processing(error.to_string()))?;
    Ok(QueryMolBlockRecord {
        query,
        substance_groups: topology.substance_groups,
        properties,
        source_coordinate_dim: coordinates.source_coordinate_dim,
    })
}

fn is_source_query_atom(atom: &QueryAtom) -> bool {
    !atom.predicate_is_carrier_derived()
}

fn is_source_query_bond(bond: &QueryBond) -> bool {
    !bond.predicate_is_carrier_derived()
}

fn synchronize_query_atom(mut source: QueryAtom, carrier: cosmolkit_model::Atom) -> QueryAtom {
    // BEGIN RDKIT CPP FUNCTION QueryOps::replaceAtomWithQueryAtom
    // RDKit✔️✔️: if (!atom->hasQuery()) {
    // RDKit✔️✔️:   auto *newAt = new QueryAtom(*atom);
    // RDKit✔️✔️:   res = static_cast<QueryAtom *>(mol.replaceAtom(atom->getIdx(), newAt));
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   res = static_cast<QueryAtom *>(atom);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryOps::replaceAtomWithQueryAtom
    // Behavior review: typed construction provenance distinguishes an
    // explicit detached query from Molfile IO's uniform wrapper around an
    // ordinary atom. Only the latter receives the constructor snapshot of the
    // final carrier; optional Molfile properties are not provenance.
    // Complexity review: one enum test and, for synthesized carriers, the same
    // bounded constructor predicate allocation as the source conversion.
    if is_source_query_atom(&source) {
        *source.atom_mut() = carrier;
        source
    } else {
        let predicate = query_from_concrete_atom_value(&carrier);
        QueryAtom::from_carrier_parts(carrier, predicate)
    }
}

fn synchronize_query_bond(mut source: QueryBond, carrier: cosmolkit_model::Bond) -> QueryBond {
    // BEGIN RDKIT CPP FUNCTION QueryBond::QueryBond(const Bond &)
    // RDKit✔️✔️: explicit QueryBond(const Bond &other)
    // RDKit✔️✔️:     : Bond(other), dp_query(makeBondOrderEqualsQuery(other.getBondType())) {}
    // END RDKIT CPP FUNCTION QueryBond::QueryBond(const Bond &)
    // Behavior review: explicit detached query types retain their predicate;
    // a Molfile-owned uniform wrapper for an ordinary source Bond is
    // reconstructed from the final processed carrier, exactly where RDKit
    // still owns an ordinary Bond. Optional properties are not provenance.
    // Complexity review: classification is one enum test and reconstruction
    // allocates one leaf predicate.
    if is_source_query_bond(&source) {
        *source.bond_mut() = carrier;
        source
    } else {
        let predicate = if carrier.order() == BondOrder::Unspecified {
            QueryNode::predicate(BondQueryPredicate::Any)
        } else {
            QueryNode::predicate(BondQueryPredicate::Order(carrier.order()))
        };
        QueryBond::from_carrier_parts(carrier, predicate)
    }
}

fn record_requires_query(record: &MolBlockRecord) -> bool {
    match record {
        MolBlockRecord::Concrete { topology, .. } => {
            topology.atoms.iter().any(|atom| {
                atom.prop("molSubstCount")
                    .and_then(|value| parse_rdkit_int(value).ok())
                    .is_some_and(|value| value != 0)
            }) || topology.substance_groups.iter().any(|group| {
                group.data().is_some_and(|data| {
                    matches!(data.query_type.as_deref(), Some("SMARTSQ" | "SQ"))
                })
            })
        }
        MolBlockRecord::Query(_) => false,
    }
}

fn promote_record_to_query(record: &mut MolBlockRecord) -> Result<(), MolPostError> {
    if record_requires_query(record) {
        let old = std::mem::replace(
            record,
            MolBlockRecord::Concrete {
                topology: TopologyBlock::default(),
                coordinates: CoordinateBlock::default(),
                properties: cosmolkit_model::MoleculeProperties::default(),
            },
        );
        let MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } = old
        else {
            unreachable!()
        };
        *record = MolBlockRecord::Query(concrete_to_query(topology, coordinates, properties)?);
    }
    Ok(())
}

fn process_smarts_groups(record: &mut MolBlockRecord) -> Result<(), MolPostError> {
    let MolBlockRecord::Query(query_record) = record else {
        return Ok(());
    };
    let mut remove = vec![false; query_record.substance_groups.len()];
    for (index, group) in query_record.substance_groups.iter().enumerate() {
        let Some(data) = group.data() else { continue };
        if !matches!(data.query_type.as_deref(), Some("SMARTSQ" | "SQ")) {
            continue;
        }
        remove[index] = true;
        if data.query_op.as_deref().is_some_and(|op| op != "=") {
            continue;
        }
        let Some(smarts) = data_values(group).first().filter(|value| !value.is_empty()) else {
            continue;
        };
        let Ok(parsed) =
            cosmolkit_search::parse_smarts(smarts, &cosmolkit_search::SmartsParseParams::default())
        else {
            continue;
        };
        if parsed.num_atoms() == 0 {
            continue;
        }
        let predicate = if parsed.num_atoms() == 1 {
            parsed.atoms()[0].predicate().clone()
        } else {
            QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
                RecursiveStructureQuery::from_query_graph(parsed, 0)
                    .with_source_smarts(smarts.clone()),
            ))
        };
        for atom_id in group.atoms() {
            if let Some(atom) = query_record.query.atom_mut(atom_id.index()) {
                atom.set_predicate(predicate.clone());
                atom.atom_mut()
                    .set_prop("_MolFileAtomQuery", "1")
                    .map_err(|error| MolPostError::Processing(error.to_string()))?;
                atom.atom_mut()
                    .set_prop("MRV SMA", smarts)
                    .map_err(|error| MolPostError::Processing(error.to_string()))?;
                atom.atom_mut()
                    .set_prop("_MolFileAtomQuery", "1")
                    .map_err(|error| MolPostError::Processing(error.to_string()))?;
            }
        }
    }
    query_record.substance_groups =
        retain_substance_groups(std::mem::take(&mut query_record.substance_groups), &remove)?;
    Ok(())
}

/// Apply the ordered Molfile postprocessing closure to a detached record.
pub fn finish_mol_block_record(
    mut record: MolBlockRecord,
    chirality_possible: bool,
    params: MolPostParams,
) -> Result<MolBlockRecord, MolPostError> {
    // `expandAttachmentPoints` is a strong topology edit whose source also
    // creates null-query atoms and terminal coordinates. Returning a typed
    // boundary error is the only honest behavior until that complete detached
    // transform is installed; no plausible partial graph is emitted.
    if params.expand_attachment_points {
        return Err(MolPostError::AttachmentPointExpansion);
    }
    promote_record_to_query(&mut record)?;
    match record {
        MolBlockRecord::Concrete {
            mut topology,
            coordinates,
            properties,
        } => {
            process_atom_properties(&mut topology, None)?;
            process_groups_on_topology(&mut topology)?;
            let (topology, coordinates, properties, _, _) = apply_stereo_and_sanitize(
                topology,
                coordinates,
                properties,
                chirality_possible,
                params,
                None,
            )?;
            Ok(MolBlockRecord::Concrete {
                topology,
                coordinates,
                properties,
            })
        }
        MolBlockRecord::Query(mut query_record) => {
            let mut topology = TopologyBlock {
                atoms: query_record
                    .query
                    .atoms()
                    .iter()
                    .map(|atom| atom.atom().clone())
                    .collect(),
                bonds: query_record
                    .query
                    .bonds()
                    .iter()
                    .map(|bond| bond.bond().clone())
                    .collect(),
                adjacency: AdjacencyList::from_topology(
                    query_record.query.num_atoms(),
                    &query_record
                        .query
                        .bonds()
                        .iter()
                        .map(|bond| bond.bond().clone())
                        .collect::<Vec<_>>(),
                ),
                substance_groups: query_record.substance_groups.clone(),
                stereo_groups: query_record.query.stereo_groups().to_vec(),
            };
            process_atom_properties(&mut topology, Some(query_record.query.atoms_mut()))?;
            process_groups_on_topology(&mut topology)?;
            for (query_atom, atom) in query_record
                .query
                .atoms_mut()
                .iter_mut()
                .zip(&topology.atoms)
            {
                *query_atom.atom_mut() = atom.clone();
            }
            for (query_bond, bond) in query_record
                .query
                .bonds_mut()
                .iter_mut()
                .zip(&topology.bonds)
            {
                *query_bond.bond_mut() = bond.clone();
            }
            query_record.substance_groups = topology.substance_groups;
            let mut wrapped = MolBlockRecord::Query(query_record);
            process_smarts_groups(&mut wrapped)?;
            let MolBlockRecord::Query(mut query_record) = wrapped else {
                unreachable!()
            };
            for (query_atom, atom) in query_record
                .query
                .atoms()
                .iter()
                .map(QueryAtom::atom)
                .cloned()
                .enumerate()
            {
                topology.atoms[query_atom] = atom;
            }
            topology.substance_groups = query_record.substance_groups.clone();

            let coordinates = query_record
                .query
                .coordinate_block(query_record.source_coordinate_dim);
            let old_query = query_record.query;
            let query_props = old_query.props().clone();
            let old_atoms = old_query.atoms().to_vec();
            let old_bonds = old_query.bonds().to_vec();
            let query_state = QueryStateRef::try_for_topology(&old_atoms, &old_bonds, &topology)
                .map_err(|error| MolPostError::Processing(error.to_string()))?;
            let (topology, coordinates, properties, _, query_rows) = apply_stereo_and_sanitize(
                topology,
                coordinates,
                query_record.properties,
                chirality_possible,
                params,
                Some(query_state),
            )?;
            let (query_atoms, query_bonds) = query_rows.ok_or(MolPostError::Representation(
                "query state missing after mol-post finalization",
            ))?;
            let query_atoms = query_atoms
                .into_iter()
                .zip(&topology.atoms)
                .map(|(query, carrier)| synchronize_query_atom(query, carrier.clone()))
                .collect();
            let query_bonds = query_bonds
                .into_iter()
                .zip(&topology.bonds)
                .map(|(query, carrier)| synchronize_query_bond(query, carrier.clone()))
                .collect();
            let source_coordinate_dim = coordinates.source_coordinate_dim;
            query_record.query = QueryGraph::from_parts(
                query_atoms,
                query_bonds,
                query_props,
                coordinates.conformers_2d,
                coordinates.conformers_3d,
                topology.stereo_groups,
            )
            .map_err(|error| MolPostError::Processing(error.to_string()))?;
            query_record.substance_groups = topology.substance_groups;
            query_record.properties = properties;
            query_record.source_coordinate_dim = source_coordinate_dim;
            if query_record.query.prop("_NeedsQueryScan").is_some() {
                query_record.query.clear_prop("_NeedsQueryScan");
                cosmolkit_search::complete_mol_queries(
                    &mut query_record.query,
                    cosmolkit_search::QUERY_SCAN_MAGIC_VALUE,
                );
            }
            Ok(MolBlockRecord::Query(query_record))
        }
    }
}
