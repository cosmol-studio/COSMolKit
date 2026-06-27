use super::*;

#[mol_op_body(sanitize, parts)]
pub(super) fn sanitize_impl(ops: crate::SanitizeOps) -> Result<OpOutcome, OperationError> {
    let mut parts = parts;
    let changed = run_sanitize_pipeline(&mut parts, ops, &SANITIZE_SPEC)?;
    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "requested sanitize steps only updated derived state",
        }
    })
}

fn sanitize_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    message: String,
    unsupported: Option<crate::UnsupportedFeatureError>,
) -> OperationError {
    OperationError::Sanitize {
        operation,
        source: crate::SanitizeError {
            step,
            message,
            unsupported,
        },
    }
}

fn sanitize_invalid_input_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    message: &'static str,
) -> OperationError {
    sanitize_error(operation, step, message.to_string(), None)
}

pub(super) fn sanitize_valence_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::ValenceError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_ring_finding_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::RingFindingError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_kekulize_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::KekulizeError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_aromaticity_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::AromaticityError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

pub(super) fn sanitize_stage<T, E>(
    step: crate::SanitizeStep,
    stage: impl FnOnce() -> Result<T, E>,
    map_error: impl FnOnce(crate::SanitizeStep, E) -> OperationError,
) -> Result<T, OperationError> {
    stage().map_err(|source| map_error(step, source))
}

fn sanitize_properties_stage(
    parts: &mut OpParts<'_>,
    topology: &TopologyBlock,
    coordinates: Option<&CoordinateBlock>,
    properties: Option<&MoleculeProperties>,
    strict: bool,
    operation: &'static MoleculeOpSpec,
) -> Result<crate::ValenceAssignment, OperationError> {
    sanitize_stage(
        crate::SanitizeStep::Properties,
        || {
            sanitize_recompute_property_cache(
                parts,
                topology,
                coordinates,
                properties,
                strict,
                crate::SanitizeStep::Properties,
                operation,
            )
        },
        |_step, error| error,
    )
}

fn sanitize_recompute_property_cache(
    parts: &mut OpParts<'_>,
    topology: &TopologyBlock,
    coordinates: Option<&CoordinateBlock>,
    properties: Option<&MoleculeProperties>,
    strict: bool,
    step: crate::SanitizeStep,
    operation: &'static MoleculeOpSpec,
) -> Result<crate::ValenceAssignment, OperationError> {
    let view = parts.read_parts_for_optional_blocks(topology.clone(), coordinates, properties)?;
    let valence = MoleculeReadParts::from_molecule(&view)
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, strict)
        .map_err(|source| sanitize_valence_error(operation, step, source))?;
    parts.set_valence_cache(valence.clone());
    Ok(valence)
}

fn run_sanitize_pipeline(
    parts: &mut OpParts<'_>,
    ops: crate::SanitizeOps,
    operation: &'static MoleculeOpSpec,
) -> Result<bool, OperationError> {
    let mut topology = parts.begin_topology_mut()?;
    let changed =
        run_sanitize_pipeline_on_topology(parts, &mut topology, None, None, ops, operation)?;
    parts.commit_topology(topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    Ok(changed)
}

pub(super) fn run_sanitize_pipeline_on_topology(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    coordinates: Option<&CoordinateBlock>,
    properties: Option<&MoleculeProperties>,
    ops: crate::SanitizeOps,
    operation: &'static MoleculeOpSpec,
) -> Result<bool, OperationError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::sanitizeMol
    // RDKit✔️✔️: void sanitizeMol(RWMol &mol) {
    // RDKit✔️✔️:   unsigned int failedOp = 0;
    // RDKit✔️✔️:   sanitizeMol(mol, failedOp, SANITIZE_ALL);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
    // RDKit✔️✔️:                  unsigned int sanitizeOps) {
    // RDKit✔️✔️:   // clear out any cached properties
    // RDKit✔️✔️:   mol.clearComputedProps();
    parts.clear_cache(SANITIZE_SPEC.needs_update());
    macro_rules! sanitize_read {
        ($read:ident => $body:expr) => {{
            let view = parts
                .read_parts_for_optional_blocks((*topology).clone(), coordinates, properties)
                .expect("sanitize working topology should satisfy molecule invariants");
            let $read = MoleculeReadParts::from_molecule(&view);
            $body
        }};
    }

    let mut changed = false;
    let property_cache_strict = ops.contains(crate::SanitizeOps::PROPERTIES);

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUP;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     // clean up things like nitro groups
    // RDKit✔️✔️:     cleanUp(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP) {
        let cleanup = sanitize_stage(
            crate::SanitizeStep::Cleanup,
            || sanitize_read!(read => sanitize_cleanup_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let cleanup_changed = sanitize_read!(read => cleanup_changes_molecule(read, &cleanup));
        if cleanup_changed {
            apply_sanitize_cleanup_assignment(topology, cleanup);
        }
        changed |= cleanup_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     cleanUpOrganometallics(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS) {
        let organometallics = sanitize_stage(
            crate::SanitizeStep::CleanupOrganometallics,
            || sanitize_read!(read => sanitize_organometallic_cleanup_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let organometallics_changed =
            sanitize_read!(read => organometallic_cleanup_changes_molecule(read, &organometallics));
        if organometallics_changed {
            apply_sanitize_organometallic_cleanup_assignment(topology, organometallics);
        }
        changed |= organometallics_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     mol.updatePropertyCache(true);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    let valence = sanitize_properties_stage(
        parts,
        topology,
        coordinates,
        properties,
        property_cache_strict,
        operation,
    )?;
    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SYMMRINGS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     VECT_INT_VECT arings;
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol, arings);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SYMMRINGS) {
        let rings = sanitize_stage(
            crate::SanitizeStep::SymmRings,
            || sanitize_read!(read => read.symmetrize_sssr()),
            |step, source| sanitize_ring_finding_error(operation, step, source),
        )?;
        parts.set_rings_cache(rings);
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_KEKULIZE;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     Kekulize(mol, true, false);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::KEKULIZE) {
        if sanitize_read!(read => read.derived_cache().rings.is_none()) {
            let rings = sanitize_stage(
                crate::SanitizeStep::Kekulize,
                || sanitize_read!(read => read.symmetrize_sssr()),
                |step, source| sanitize_ring_finding_error(operation, step, source),
            )?;
            parts.set_rings_cache(rings);
        }
        let ring_info = sanitize_read!(read => {
            read.derived_cache()
                .rings
                .as_ref()
                .expect("rings were recomputed immediately above")
                .clone()
        });
        let assignment = sanitize_stage(
            crate::SanitizeStep::Kekulize,
            || sanitize_read!(read => read.kekulize_assignment(Some(&ring_info), true, false, 100)),
            |step, source| sanitize_kekulize_error(operation, step, source),
        )?;
        let kekulize_changed = crate::kekulize::apply_kekulize_assignment(topology, &assignment);
        if kekulize_changed {
            sanitize_stage(
                crate::SanitizeStep::Kekulize,
                || {
                    sanitize_recompute_property_cache(
                        parts,
                        topology,
                        coordinates,
                        properties,
                        property_cache_strict,
                        crate::SanitizeStep::Kekulize,
                        operation,
                    )
                },
                |_step, error| error,
            )?;
            parts.clear_cache(
                DerivedState::AROMATICITY | DerivedState::DRAWING | DerivedState::FINGERPRINT,
            );
        }
        changed |= kekulize_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_FINDRADICALS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     assignRadicals(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::FIND_RADICALS) {
        let radicals = sanitize_stage(
            crate::SanitizeStep::FindRadicals,
            || sanitize_read!(read => read.assign_radicals()),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let radicals_changed = sanitize_read!(read => read
            .atoms()
            .iter()
            .zip(radicals.iter().copied())
            .any(|(atom, radical)| atom.radical_electrons() != radical));
        if radicals_changed {
            for (atom, radical) in topology.atoms.iter_mut().zip(radicals) {
                atom.set_radical_electrons(radical);
            }
        }
        changed |= radicals_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SETAROMATICITY;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     setAromaticity(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SET_AROMATICITY) {
        if sanitize_read!(read => read.derived_cache().rings.is_none()) {
            let rings = sanitize_stage(
                crate::SanitizeStep::SetAromaticity,
                || sanitize_read!(read => read.symmetrize_sssr()),
                |step, source| sanitize_ring_finding_error(operation, step, source),
            )?;
            parts.set_rings_cache(rings);
        }
        let assignment = sanitize_stage(
            crate::SanitizeStep::SetAromaticity,
            || sanitize_read!(read => read.set_aromaticity(crate::AromaticityModel::Default)),
            |step, source| sanitize_aromaticity_error(operation, step, source),
        )?;
        for (atom, is_aromatic) in topology
            .atoms
            .iter_mut()
            .zip(assignment.atom_aromatic.iter().copied())
        {
            atom.set_aromatic(is_aromatic);
        }
        for (bond, is_aromatic) in topology
            .bonds
            .iter_mut()
            .zip(assignment.bond_aromatic.iter().copied())
        {
            bond.set_aromatic(is_aromatic);
            if is_aromatic
                && matches!(
                    bond.order(),
                    crate::BondOrder::Single | crate::BondOrder::Double
                )
            {
                bond.set_order(crate::BondOrder::Aromatic);
            }
        }
        parts.mark_aromaticity_valid();
        parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        changed = true;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SETCONJUGATION;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     setConjugation(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SET_CONJUGATION) {
        let conjugation = sanitize_stage(
            crate::SanitizeStep::SetConjugation,
            || sanitize_read!(read => sanitize_conjugation_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let conjugation_changed = sanitize_read!(read => read
            .bonds()
            .iter()
            .zip(conjugation.iter().copied())
            .any(|(bond, is_conjugated)| bond.is_conjugated() != is_conjugated));
        if conjugation_changed {
            for (bond, is_conjugated) in topology.bonds.iter_mut().zip(conjugation) {
                bond.set_conjugated(is_conjugated);
            }
            parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        }
        changed |= conjugation_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SETHYBRIDIZATION;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     setHybridization(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SET_HYBRIDIZATION) {
        let hybridization = sanitize_stage(
            crate::SanitizeStep::SetHybridization,
            || sanitize_read!(read => sanitize_hybridization_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let hybridization_changed = sanitize_read!(read => read
            .atoms()
            .iter()
            .zip(hybridization.iter().copied())
            .any(|(atom, hybridization)| atom.hybridization() != hybridization));
        if hybridization_changed {
            for (atom, hybridization) in topology.atoms.iter_mut().zip(hybridization) {
                atom.set_hybridization(hybridization);
            }
            parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        }
        changed |= hybridization_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUPATROPISOMERS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     cleanupAtropisomers(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_ATROPISOMERS) {
        if sanitize_read!(read => read.derived_cache().rings.is_none()) {
            let rings = sanitize_stage(
                crate::SanitizeStep::CleanupAtropisomers,
                || sanitize_read!(read => read.symmetrize_sssr()),
                |step, source| sanitize_ring_finding_error(operation, step, source),
            )?;
            parts.set_rings_cache(rings);
        }
        let ring_info = sanitize_read!(read => {
            read.derived_cache()
                .rings
                .as_ref()
                .expect("rings were recomputed immediately above")
                .clone()
        });
        let atropisomer_cleanup =
            sanitize_read!(read => sanitize_cleanup_atropisomers_assignment(read, &ring_info));
        let atropisomer_cleanup_changed = sanitize_read!(read => {
            atropisomer_cleanup
                .iter()
                .any(|(bond, _)| read.bonds()[bond.index()].stereo() != crate::BondStereo::None)
        });
        if atropisomer_cleanup_changed {
            apply_sanitize_cleanup_atropisomers_assignment(topology, &atropisomer_cleanup);
            parts.clear_cache(
                DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT,
            );
        }
        changed |= atropisomer_cleanup_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     cleanupChirality(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_CHIRALITY) {
        let chirality_cleanup = sanitize_stage(
            crate::SanitizeStep::CleanupChirality,
            || sanitize_read!(read => sanitize_cleanup_chirality_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let chirality_cleanup_changed =
            sanitize_read!(read => cleanup_chirality_changes_molecule(read, &chirality_cleanup));
        if chirality_cleanup_changed {
            apply_sanitize_cleanup_chirality_assignment(topology, &chirality_cleanup);
            parts.clear_cache(
                DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT,
            );
        }
        changed |= chirality_cleanup_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_ADJUSTHS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     adjustHs(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::ADJUST_HYDROGENS) {
        let hydrogen_adjustment = sanitize_stage(
            crate::SanitizeStep::AdjustHydrogens,
            || {
                sanitize_read!(read => {
                    sanitize_adjust_hydrogens_assignment(read)
                })
            },
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let hydrogen_adjustment_changed =
            sanitize_read!(read => adjust_hydrogens_changes_molecule(read, &hydrogen_adjustment));
        if hydrogen_adjustment_changed {
            apply_sanitize_adjust_hydrogens_assignment(topology, &hydrogen_adjustment);
            parts.clear_cache(DerivedState::VALENCE);
            parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        }
        changed |= hydrogen_adjustment_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     mol.updatePropertyCache(true);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::PROPERTIES) {
        sanitize_properties_stage(parts, topology, coordinates, properties, true, operation)?;
    }
    // RDKit✔️✔️:   operationThatFailed = 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::sanitizeMol
    Ok(changed)
}

pub(crate) fn sanitize_conjugation_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<bool>, crate::ValenceError> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION MolOps::setConjugation
    // RDKit✔️✔️: void setConjugation(ROMol &mol) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     bond->setIsConjugated(bond->getIsAromatic());
    // RDKit✔️✔️:   }
    let adjacency = sanitize_adjacency(molecule)?;
    let valence = sanitize_valence_facts(molecule)?;
    let mut conjugated = molecule
        .bonds()
        .iter()
        .map(crate::Bond::is_aromatic)
        .collect::<Vec<_>>();

    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     markConjAtomBonds(atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for atom in molecule.atoms() {
        sanitize_mark_conjugated_atom_bonds(
            molecule,
            &adjacency,
            &valence,
            atom.id(),
            &mut conjugated,
        )?;
    }
    Ok(conjugated)
    // END RDKIT CPP FUNCTION MolOps::setConjugation
}

fn sanitize_mark_conjugated_atom_bonds<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
    conjugated: &mut [bool],
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION markConjAtomBonds
    // RDKit✔️✔️: void markConjAtomBonds(Atom *at) {
    // RDKit✔️✔️:   if (!isAtomConjugCand(at)) { return; }
    if !sanitize_is_atom_conjugation_candidate(molecule, adjacency, valence, atom)? {
        return Ok(());
    }

    // RDKit✔️✔️:   int sbo = at->getDegree() + at->getTotalNumHs();
    // RDKit✔️✔️:   if ((sbo < 2) || (sbo > 3)) { return; }
    let sbo = sanitize_total_degree(molecule, adjacency, valence, atom)?;
    if !(2..=3).contains(&sbo) {
        return Ok(());
    }

    // RDKit✔️✔️:   for (const auto bnd1 : mol.atomBonds(at)) {
    for neighbor_1 in adjacency.neighbors_of(atom.index()) {
        let bond_1 = molecule.bonds().get(neighbor_1.bond.index()).ok_or(
            crate::ValenceError::UnsupportedBranch {
                reason: "conjugation adjacency bond index out of range",
            },
        )?;
        let other_1 = AtomId::new(neighbor_1.atom_index);
        // RDKit✔️✔️:     if (bnd1->getValenceContrib(at) < 1.5 ||
        // RDKit✔️✔️:         !isAtomConjugCand(bnd1->getOtherAtom(at))) { continue; }
        if crate::valence::bond_valence_contrib(bond_1, atom)? < 1.5
            || !sanitize_is_atom_conjugation_candidate(molecule, adjacency, valence, other_1)?
        {
            continue;
        }

        // RDKit✔️✔️:     for (const auto bnd2 : mol.atomBonds(at)) {
        for neighbor_2 in adjacency.neighbors_of(atom.index()) {
            if neighbor_1.bond == neighbor_2.bond {
                continue;
            }
            let atom_2 = AtomId::new(neighbor_2.atom_index);
            // RDKit✔️✔️:       sbo = at2->getDegree() + at2->getTotalNumHs();
            // RDKit✔️✔️:       if (sbo > 3) { continue; }
            if sanitize_total_degree(molecule, adjacency, valence, atom_2)? > 3 {
                continue;
            }
            // RDKit✔️✔️:       if (isAtomConjugCand(at2)) {
            // RDKit✔️✔️:         bnd1->setIsConjugated(true);
            // RDKit✔️✔️:         bnd2->setIsConjugated(true);
            // RDKit✔️✔️:       }
            if sanitize_is_atom_conjugation_candidate(molecule, adjacency, valence, atom_2)? {
                conjugated[neighbor_1.bond.index()] = true;
                conjugated[neighbor_2.bond.index()] = true;
            }
        }
    }
    Ok(())
    // END RDKIT CPP FUNCTION markConjAtomBonds
}

fn sanitize_is_atom_conjugation_candidate<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<bool, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION isAtomConjugCand
    // RDKit✔️✔️: bool isAtomConjugCand(const Atom *at) {
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "conjugation atom index out of range",
        })?;
    // RDKit✔️✔️:   const auto &vals =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
    let valence_list = crate::rdkit_valence_list(atom_data.atomic_number())?.ok_or(
        crate::ValenceError::UnsupportedBranch {
            reason: "conjugation valence-list atomic number out of range",
        },
    )?;
    let total_valence = sanitize_total_valence(valence, atom);
    // RDKit✔️✔️:   if (!at->getFormalCharge() && vals.front() >= 0 &&
    // RDKit✔️✔️:       at->getTotalValence() > static_cast<unsigned int>(vals.front())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if atom_data.formal_charge() == 0
        && valence_list.first().is_some_and(|value| *value >= 0)
        && total_valence > valence_list[0]
    {
        return Ok(false);
    }
    // RDKit✔️✔️:   int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:   auto res = ((at->getAtomicNum() <= 10) || (nouter != 5 && nouter != 6) ||
    // RDKit✔️✔️:               (nouter == 6 && at->getTotalDegree() < 2u)) &&
    // RDKit✔️✔️:              MolOps::countAtomElec(at) > 0;
    let n_outer = crate::valence::periodic_table_outer_electrons(atom_data.atomic_number())?;
    Ok(((atom_data.atomic_number() <= 10)
        || (n_outer != 5 && n_outer != 6)
        || (n_outer == 6 && sanitize_total_degree(molecule, adjacency, valence, atom)? < 2))
        && sanitize_count_atom_electrons(molecule, adjacency, valence, atom)? > 0)
    // END RDKIT CPP FUNCTION isAtomConjugCand
}

pub(crate) fn sanitize_hybridization_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<crate::Hybridization>, crate::ValenceError> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION MolOps::setHybridization
    // RDKit✔️✔️: void setHybridization(ROMol &mol) {
    let adjacency = sanitize_adjacency(molecule)?;
    let valence = sanitize_valence_facts(molecule)?;
    molecule
        .atoms()
        .iter()
        .map(|atom| {
            // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
            // RDKit✔️✔️:     if (atom->getAtomicNum() == 0) {
            // RDKit✔️✔️:       atom->setHybridization(Atom::UNSPECIFIED);
            if atom.atomic_number() == 0 {
                return Ok(crate::Hybridization::Unspecified);
            }

            let total_degree = sanitize_total_degree(molecule, &adjacency, &valence, atom.id())?;
            // RDKit✔️✔️:       switch (atom->getChiralTag()) { ... stereo spec matches coordination number ... }
            match atom.chiral_tag() {
                crate::ChiralTag::Tetrahedral
                | crate::ChiralTag::TetrahedralCw
                | crate::ChiralTag::TetrahedralCcw
                    if total_degree == 4 =>
                {
                    return Ok(crate::Hybridization::Sp3);
                }
                crate::ChiralTag::SquarePlanar if (2..=4).contains(&total_degree) => {
                    return Ok(crate::Hybridization::Sp2d);
                }
                crate::ChiralTag::TrigonalBipyramidal if (2..=5).contains(&total_degree) => {
                    return Ok(crate::Hybridization::Sp3d);
                }
                crate::ChiralTag::Octahedral if (2..=6).contains(&total_degree) => {
                    return Ok(crate::Hybridization::Sp3d2);
                }
                _ => {}
            }

            // RDKit✔️✔️:       if (atom->getAtomicNum() < 89) {
            // RDKit✔️✔️:         norbs = numBondsPlusLonePairs(atom);
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         norbs = atom->getTotalDegree();
            // RDKit✔️✔️:       }
            let norbs = if atom.atomic_number() < 89 {
                sanitize_num_bonds_plus_lone_pairs(molecule, &adjacency, &valence, atom.id())?
            } else {
                total_degree
            };
            // RDKit✔️✔️:       switch (norbs) { ... }
            Ok(match norbs {
                0 | 1 => crate::Hybridization::S,
                2 => crate::Hybridization::Sp,
                3 => crate::Hybridization::Sp2,
                4 => {
                    if total_degree > 3 || !sanitize_atom_has_conjugated_bond(molecule, atom.id()) {
                        crate::Hybridization::Sp3
                    } else {
                        crate::Hybridization::Sp2
                    }
                }
                5 => crate::Hybridization::Sp3d,
                6 => crate::Hybridization::Sp3d2,
                _ => crate::Hybridization::Unspecified,
            })
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
        })
        .collect()
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::setHybridization
}

fn sanitize_num_bonds_plus_lone_pairs<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION numBondsPlusLonePairs
    // RDKit✔️✔️: int numBondsPlusLonePairs(Atom *at) {
    // RDKit✔️✔️:   int deg = at->getTotalDegree();
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "hybridization atom index out of range",
        })?;
    let mut degree = sanitize_total_degree(molecule, adjacency, valence, atom)?;
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::ZERO ||
    // RDKit✔️✔️:         (isDative(*bond) && at->getIdx() != bond->getEndAtomIdx())) { --deg; }
    // RDKit✔️✔️:   }
    for neighbor in adjacency.neighbors_of(atom.index()) {
        let bond = molecule.bonds().get(neighbor.bond.index()).ok_or(
            crate::ValenceError::UnsupportedBranch {
                reason: "hybridization adjacency bond index out of range",
            },
        )?;
        if bond.order() == crate::BondOrder::Zero
            || (sanitize_is_dative_bond(bond.order()) && atom != bond.end())
        {
            degree -= 1;
        }
    }
    // RDKit✔️✔️:   if (at->getAtomicNum() <= 1) { return deg; }
    if atom_data.atomic_number() <= 1 {
        return Ok(degree);
    }
    // RDKit✔️✔️:   int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:   int totalValence = at->getTotalValence();
    // RDKit✔️✔️:   int chg = at->getFormalCharge();
    let n_outer = crate::valence::periodic_table_outer_electrons(atom_data.atomic_number())?;
    let total_valence = sanitize_total_valence(valence, atom);
    let charge = i32::from(atom_data.formal_charge());
    // RDKit✔️✔️:   int numFreeElectrons = nouter - (totalValence + chg);
    let free_electrons = n_outer - (total_valence + charge);
    // RDKit✔️✔️:   if (totalValence + nouter - chg < 8) { ... } else { ... }
    if total_valence + n_outer - charge < 8 {
        let radicals = i32::from(atom_data.radical_electrons());
        let lone_pairs = (free_electrons - radicals) / 2;
        Ok(degree + lone_pairs + radicals)
    } else {
        Ok(degree + free_electrons / 2)
    }
    // END RDKIT CPP FUNCTION numBondsPlusLonePairs
}

fn sanitize_count_atom_electrons<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "countAtomElec atom index out of range",
        })?;
    let default_valence = crate::valence::rdkit_default_valence(atom_data.atomic_number())?;
    if default_valence <= 1 {
        return Ok(-1);
    }
    let mut degree = sanitize_total_degree(molecule, adjacency, valence, atom)?;
    for neighbor in adjacency.neighbors_of(atom.index()) {
        let bond = molecule.bonds().get(neighbor.bond.index()).ok_or(
            crate::ValenceError::UnsupportedBranch {
                reason: "countAtomElec adjacency bond index out of range",
            },
        )?;
        if crate::valence::bond_valence_contrib(bond, atom)? == 0.0
            && !sanitize_is_bond_order_query(bond.query())
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return Ok(-1);
    }
    let n_outer = crate::valence::periodic_table_outer_electrons(atom_data.atomic_number())?;
    let lone_pairs = (n_outer - default_valence - i32::from(atom_data.formal_charge())).max(0);
    let radicals = i32::from(atom_data.radical_electrons());
    let mut result = (default_valence - degree) + lone_pairs - radicals;
    if result > 1 {
        let graph_degree =
            i32::try_from(adjacency.neighbors_of(atom.index()).len()).map_err(|_| {
                crate::ValenceError::UnsupportedBranch {
                    reason: "countAtomElec atom degree out of range",
                }
            })?;
        if valence.explicit_valence[atom.index()] - graph_degree > 1 {
            result = 1;
        }
    }
    Ok(result)
}

fn sanitize_valence_facts<'a>(
    molecule: impl MoleculeReadAccess<'a>,
) -> Result<crate::ValenceAssignment, crate::ValenceError> {
    if let Some(valence) = &molecule.derived_cache().valence {
        Ok(valence.clone())
    } else {
        molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
    }
}

pub(super) fn sanitize_adjacency<'a>(
    molecule: impl MoleculeReadAccess<'a>,
) -> Result<crate::AdjacencyList, crate::ValenceError> {
    crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds()).map_err(|_| {
        crate::ValenceError::UnsupportedBranch {
            reason: "sanitize topology bond atom index out of range",
        }
    })
}

fn sanitize_total_degree<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "sanitize total-degree atom index out of range",
        })?;
    let graph_degree = i32::try_from(adjacency.neighbors_of(atom.index()).len()).map_err(|_| {
        crate::ValenceError::UnsupportedBranch {
            reason: "sanitize atom degree out of range",
        }
    })?;
    Ok(graph_degree
        + i32::from(atom_data.explicit_hydrogens())
        + valence.implicit_hydrogens[atom.index()].max(0))
}

fn sanitize_total_valence(valence: &crate::ValenceAssignment, atom: AtomId) -> i32 {
    valence.explicit_valence[atom.index()] + valence.implicit_hydrogens[atom.index()].max(0)
}

fn sanitize_atom_has_conjugated_bond<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    atom: AtomId,
) -> bool {
    molecule
        .bonds()
        .iter()
        .any(|bond| bond.is_conjugated() && (bond.begin() == atom || bond.end() == atom))
}

fn sanitize_is_dative_bond(order: crate::BondOrder) -> bool {
    matches!(
        order,
        crate::BondOrder::Dative
            | crate::BondOrder::DativeOne
            | crate::BondOrder::DativeLeft
            | crate::BondOrder::DativeRight
    )
}

pub(super) fn sanitize_is_bond_order_query(
    query: Option<&crate::QueryNode<crate::BondQueryPredicate>>,
) -> bool {
    match query {
        Some(query) => {
            crate::valence::has_bond_type_query(query)
                && !crate::valence::has_complex_bond_type_query(query)
        }
        None => false,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SanitizeCleanupAssignment {
    pub(super) atom_formal_charges: Vec<i8>,
    pub(super) bond_orders: Vec<crate::BondOrder>,
}

fn sanitize_cleanup_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<SanitizeCleanupAssignment, crate::ValenceError> {
    let molecule = read_parts;
    let adjacency = sanitize_adjacency(molecule)?;
    let mut assignment = SanitizeCleanupAssignment {
        atom_formal_charges: molecule
            .atoms()
            .iter()
            .map(crate::Atom::formal_charge)
            .collect(),
        bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
    };

    sanitize_nitrogens_cleanup_assignment(molecule, &adjacency, &mut assignment)?;
    for atom in molecule.atoms() {
        match atom.atomic_number() {
            15 => sanitize_phosphorus_cleanup_assignment(
                molecule,
                &adjacency,
                atom.id(),
                &mut assignment,
            )?,
            17 | 35 | 53 => sanitize_halogen_cleanup_assignment(
                molecule,
                &adjacency,
                atom.id(),
                &mut assignment,
            )?,
            _ => {}
        }
    }

    Ok(assignment)
}

fn sanitize_nitrogens_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    assignment: &mut SanitizeCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION nitrogensCleanup
    // RDKit✔️✔️: void nitrogensCleanup(RWMol &mol) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> nitrogensToConsider(mol.getNumAtoms());
    let mut nitrogens_to_consider = Vec::new();

    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->getAtomicNum() != 7) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getFormalCharge()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nitrogensToConsider.set(atom->getIdx());
    for atom in molecule.atoms() {
        if atom.atomic_number() != 7 || assignment.atom_formal_charges[atom.id().index()] != 0 {
            continue;
        }
        if sanitize_cleanup_explicit_valence(molecule, adjacency, assignment, atom.id())? != 5 {
            continue;
        }
        nitrogens_to_consider.push(atom.id());

        // RDKit✔️✔️:     bool updateNeeded = false;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
        // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
        // RDKit✔️✔️:            Bond::DOUBLE)) {
        // RDKit✔️✔️:         auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        // RDKit✔️✔️:         b->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:         atom->setFormalCharge(1);
        // RDKit✔️✔️:         nbr->setFormalCharge(-1);
        // RDKit✔️✔️:         updateNeeded = true;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom.id()) {
            let bond = &molecule.bonds()[bond_index];
            let neighbor = sanitize_cleanup_other_atom(bond, atom.id());
            let neighbor_atom = &molecule.atoms()[neighbor.index()];
            if neighbor_atom.atomic_number() == 8
                && assignment.atom_formal_charges[neighbor.index()] == 0
                && assignment.bond_orders[bond_index] == crate::BondOrder::Double
            {
                assignment.bond_orders[bond_index] = crate::BondOrder::Single;
                assignment.atom_formal_charges[atom.id().index()] = 1;
                assignment.atom_formal_charges[neighbor.index()] = -1;
                break;
            }
        }
    }

    // RDKit✔️✔️:   for (auto aid = nitrogensToConsider.find_first();
    // RDKit✔️✔️:        aid != boost::dynamic_bitset<>::npos;
    // RDKit✔️✔️:        aid = nitrogensToConsider.find_next(aid)) {
    for atom in nitrogens_to_consider {
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 7) && (nbr->getFormalCharge() == 0) &&
        // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
        // RDKit✔️✔️:            Bond::TRIPLE)) {
        // RDKit✔️✔️:         auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        // RDKit✔️✔️:         b->setBondType(Bond::DOUBLE);
        // RDKit✔️✔️:         atom->setFormalCharge(1);
        // RDKit✔️✔️:         nbr->setFormalCharge(-1);
        // RDKit✔️✔️:         updateNeeded = true;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
            let bond = &molecule.bonds()[bond_index];
            let neighbor = sanitize_cleanup_other_atom(bond, atom);
            let neighbor_atom = &molecule.atoms()[neighbor.index()];
            if neighbor_atom.atomic_number() == 7
                && assignment.atom_formal_charges[neighbor.index()] == 0
                && assignment.bond_orders[bond_index] == crate::BondOrder::Triple
            {
                assignment.bond_orders[bond_index] = crate::BondOrder::Double;
                assignment.atom_formal_charges[atom.index()] = 1;
                assignment.atom_formal_charges[neighbor.index()] = -1;
                break;
            }
        }
    }
    // END RDKIT CPP FUNCTION nitrogensCleanup
    Ok(())
}

fn sanitize_phosphorus_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    atom: AtomId,
    assignment: &mut SanitizeCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION phosphorusCleanup
    // RDKit✔️✔️: void phosphorusCleanup(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->getFormalCharge()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if assignment.atom_formal_charges[atom.index()] != 0 {
        return Ok(());
    }
    // RDKit✔️✔️:   if (atom->calcExplicitValence(false) == 5 && atom->getDegree() == 3) {
    if sanitize_cleanup_explicit_valence(molecule, adjacency, assignment, atom)? != 5
        || sanitize_cleanup_degree(adjacency, atom) != 3
    {
        return Ok(());
    }

    // RDKit✔️✔️:     Bond *dbl_to_O = nullptr;
    // RDKit✔️✔️:     Atom *O_atom = nullptr;
    // RDKit✔️✔️:     bool hasDoubleToCorN = false;
    let mut double_to_oxygen = None;
    let mut oxygen_atom = None;
    let mut has_double_to_carbon_or_nitrogen = false;
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
    for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
        let bond = &molecule.bonds()[bond_index];
        let neighbor = sanitize_cleanup_other_atom(bond, atom);
        let neighbor_atom = &molecule.atoms()[neighbor.index()];
        if neighbor_atom.atomic_number() == 8
            && assignment.atom_formal_charges[neighbor.index()] == 0
            && assignment.bond_orders[bond_index] == crate::BondOrder::Double
        {
            // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
            // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:            Bond::DOUBLE)) {
            // RDKit✔️✔️:         dbl_to_O = mol.getBondBetweenAtoms(aid, nbr->getIdx());
            // RDKit✔️✔️:         O_atom = nbr;
            double_to_oxygen = Some(bond_index);
            oxygen_atom = Some(neighbor);
        } else if matches!(neighbor_atom.atomic_number(), 6 | 7)
            && sanitize_cleanup_degree(adjacency, neighbor) >= 2
            && assignment.bond_orders[bond_index] == crate::BondOrder::Double
        {
            // RDKit✔️✔️:       } else if ((nbr->getAtomicNum() == 6 || nbr->getAtomicNum() == 7) &&
            // RDKit✔️✔️:                  (nbr->getDegree() >= 2) &&
            // RDKit✔️✔️:                  (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:                   Bond::DOUBLE)) {
            // RDKit✔️✔️:         hasDoubleToCorN = true;
            has_double_to_carbon_or_nitrogen = true;
        }
    }

    // RDKit✔️✔️:     if (hasDoubleToCorN && dbl_to_O != nullptr) {
    // RDKit✔️✔️:       O_atom->setFormalCharge(-1);
    // RDKit✔️✔️:       dbl_to_O->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:       atom->setFormalCharge(1);
    // RDKit✔️✔️:     }
    if has_double_to_carbon_or_nitrogen {
        if let (Some(bond_index), Some(oxygen)) = (double_to_oxygen, oxygen_atom) {
            assignment.atom_formal_charges[oxygen.index()] = -1;
            assignment.bond_orders[bond_index] = crate::BondOrder::Single;
            assignment.atom_formal_charges[atom.index()] = 1;
        }
    }
    // END RDKIT CPP FUNCTION phosphorusCleanup
    Ok(())
}

fn sanitize_halogen_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    atom: AtomId,
    assignment: &mut SanitizeCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION halogenCleanup
    // RDKit✔️✔️: void halogenCleanup(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   int ev = atom->calcExplicitValence(false);
    // RDKit✔️✔️:   if (atom->getFormalCharge() == 0 && (ev == 7 || ev == 5 || ev == 3)) {
    let explicit_valence =
        sanitize_cleanup_explicit_valence(molecule, adjacency, assignment, atom)?;
    if assignment.atom_formal_charges[atom.index()] != 0 || !matches!(explicit_valence, 3 | 5 | 7) {
        return Ok(());
    }

    // RDKit✔️✔️:     bool neighborsAllO = true;
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
    // RDKit✔️✔️:       if (nbr->getAtomicNum() != 8) {
    // RDKit✔️✔️:         neighborsAllO = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    if !sanitize_cleanup_incident_bonds(adjacency, atom)
        .into_iter()
        .all(|bond_index| {
            let neighbor = sanitize_cleanup_other_atom(&molecule.bonds()[bond_index], atom);
            molecule.atoms()[neighbor.index()].atomic_number() == 8
        })
    {
        return Ok(());
    }

    // RDKit✔️✔️:     if (neighborsAllO) {
    // RDKit✔️✔️:       int formalCharge = 0;
    // RDKit✔️✔️:       for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:         if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:           auto otherAtom = bond->getOtherAtom(atom);
    // RDKit✔️✔️:           formalCharge++;
    // RDKit✔️✔️:           otherAtom->setFormalCharge(-1);
    // RDKit✔️✔️:           otherAtom->calcExplicitValence(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->setFormalCharge(formalCharge);
    let mut formal_charge = 0;
    for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
        if assignment.bond_orders[bond_index] == crate::BondOrder::Double {
            assignment.bond_orders[bond_index] = crate::BondOrder::Single;
            let other_atom = sanitize_cleanup_other_atom(&molecule.bonds()[bond_index], atom);
            formal_charge += 1;
            assignment.atom_formal_charges[other_atom.index()] = -1;
        }
    }
    assignment.atom_formal_charges[atom.index()] = formal_charge;
    // END RDKIT CPP FUNCTION halogenCleanup
    Ok(())
}

pub(super) fn sanitize_cleanup_explicit_valence<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    assignment: &SanitizeCleanupAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    let mut accum = f64::from(molecule.atoms()[atom.index()].explicit_hydrogens());
    for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
        let mut bond = molecule.bonds()[bond_index].clone();
        bond.set_order(assignment.bond_orders[bond_index]);
        accum += crate::valence::bond_valence_contrib(&bond, atom)?;
    }
    Ok((accum + 0.1) as i32)
}

fn sanitize_cleanup_degree(adjacency: &crate::AdjacencyList, atom: AtomId) -> usize {
    adjacency.neighbors_of(atom.index()).len()
}

pub(super) fn sanitize_cleanup_incident_bonds(
    adjacency: &crate::AdjacencyList,
    atom: AtomId,
) -> Vec<usize> {
    // RDKit✔️✔️: cleanUp() neighbors/bonds helper equivalent
    adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| neighbor.bond.index())
        .collect()
}

fn sanitize_cleanup_other_atom(bond: &crate::Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom {
        bond.end()
    } else {
        bond.begin()
    }
}

fn cleanup_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &SanitizeCleanupAssignment,
) -> bool {
    let molecule = read_parts;
    molecule
        .atoms()
        .iter()
        .zip(assignment.atom_formal_charges.iter().copied())
        .any(|(atom, formal_charge)| atom.formal_charge() != formal_charge)
        || molecule
            .bonds()
            .iter()
            .zip(assignment.bond_orders.iter().copied())
            .any(|(bond, order)| bond.order() != order)
}

fn apply_sanitize_cleanup_assignment(
    topology: &mut TopologyBlock,
    assignment: SanitizeCleanupAssignment,
) {
    for (atom, formal_charge) in topology
        .atoms
        .iter_mut()
        .zip(assignment.atom_formal_charges)
    {
        atom.set_formal_charge(formal_charge);
    }
    for (bond, order) in topology.bonds.iter_mut().zip(assignment.bond_orders) {
        bond.set_order(order);
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SanitizeOrganometallicCleanupAssignment {
    pub(super) bond_orders: Vec<crate::BondOrder>,
    pub(super) bond_endpoints: Vec<(AtomId, AtomId)>,
}

fn sanitize_organometallic_cleanup_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<SanitizeOrganometallicCleanupAssignment, crate::ValenceError> {
    let molecule = read_parts;
    let mut assignment = SanitizeOrganometallicCleanupAssignment {
        bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
        bond_endpoints: molecule
            .bonds()
            .iter()
            .map(|bond| (bond.begin(), bond.end()))
            .collect(),
    };

    // BEGIN RDKIT CPP FUNCTION MolOps::cleanUpOrganometallics
    // RDKit✔️✔️: void cleanUpOrganometallics(RWMol &mol) {
    // RDKit✔️✔️:   bool needsFixing = false;
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (isHypervalentNonMetal(atom) && !noDative(atom)) {
    let valence = molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?;
    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| crate::ValenceError::UnsupportedBranch {
            reason: "organometallic cleanup adjacency could not be built",
        })?;

    let mut needs_fixing = false;
    for atom in molecule.atoms() {
        if sanitize_is_hypervalent_nonmetal(molecule, &adjacency, &valence, atom.id())?
            && !sanitize_no_dative(atom)
        {
            // RDKit✔️✔️:       for (auto bond : mol.atomBonds(atom)) {
            // RDKit✔️✔️:         if (bond->getBondType() == Bond::BondType::SINGLE &&
            // RDKit✔️✔️:             QueryOps::isMetal(*bond->getOtherAtom(atom))) {
            // RDKit✔️✔️:           needsFixing = true;
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            needs_fixing = !sanitize_organometallic_single_bonded_metals(
                molecule,
                &adjacency,
                &assignment,
                atom.id(),
            )
            .is_empty();
        }
        if needs_fixing {
            break;
        }
    }
    // RDKit✔️✔️:   if (!needsFixing) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if !needs_fixing {
        return Ok(assignment);
    }

    // RDKit✔️✔️:   mol.updatePropertyCache(false);
    // RDKit✔️✔️:   std::vector<unsigned int> ranks(mol.getNumAtoms());
    // RDKit✔️✔️:   RDKit::Canon::rankMolAtoms(mol, ranks);
    let ranks = MoleculeReadAccess::rank_mol_atoms(molecule).map_err(|_| {
        crate::ValenceError::UnsupportedBranch {
            reason: "organometallic cleanup canonical ranking failed",
        }
    })?;
    // RDKit✔️✔️:   std::sort(atom_ranks.begin(), atom_ranks.end(),
    // RDKit✔️✔️:             [](const std::pair<int, int> &p1, std::pair<int, int> &p2) -> bool {
    // RDKit✔️✔️:               return p1.second < p2.second;
    // RDKit✔️✔️:             });
    let mut atom_order: Vec<usize> = (0..molecule.num_atoms()).collect();
    atom_order.sort_by_key(|atom| ranks[*atom]);
    // RDKit✔️✔️:   for (auto ar : atom_ranks) {
    // RDKit✔️✔️:     auto atom = mol.getAtomWithIdx(ar.first);
    // RDKit✔️✔️:     metalBondCleanup(mol, atom, ranks);
    // RDKit✔️✔️:   }
    for atom_index in atom_order {
        sanitize_metal_bond_cleanup_assignment(
            molecule,
            &adjacency,
            &valence,
            &ranks,
            AtomId::new(atom_index),
            &mut assignment,
        )?;
    }
    // END RDKIT CPP FUNCTION MolOps::cleanUpOrganometallics

    Ok(assignment)
}

pub(super) fn sanitize_metal_bond_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    ranks: &[usize],
    atom: AtomId,
    assignment: &mut SanitizeOrganometallicCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION metalBondCleanup
    // RDKit✔️✔️: void metalBondCleanup(RWMol &mol, Atom *atom,
    // RDKit✔️✔️:                       const std::vector<unsigned int> &ranks) {
    // RDKit✔️✔️:   if (isHypervalentNonMetal(atom) && !noDative(atom)) {
    if !sanitize_is_hypervalent_nonmetal(molecule, adjacency, valence, atom)?
        || sanitize_no_dative(&molecule.atoms()[atom.index()])
    {
        return Ok(());
    }

    // RDKit✔️✔️:     std::vector<Atom *> metals;
    // RDKit✔️✔️:     for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:           QueryOps::isMetal(*bond->getOtherAtom(atom))) {
    // RDKit✔️✔️:         metals.push_back(bond->getOtherAtom(atom));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    let mut metals =
        sanitize_organometallic_single_bonded_metals(molecule, adjacency, assignment, atom);
    if metals.is_empty() {
        return Ok(());
    }

    // RDKit✔️✔️:       std::sort(metals.begin(), metals.end(),
    // RDKit✔️✔️:                 [&](const Atom *a1, const Atom *a2) -> bool {
    // RDKit✔️✔️:                   int nda1 = numDativeBonds(a1);
    // RDKit✔️✔️:                   int nda2 = numDativeBonds(a2);
    // RDKit✔️✔️:                   if (nda1 == nda2) {
    // RDKit✔️✔️:                     return ranks[a1->getIdx()] > ranks[a2->getIdx()];
    // RDKit✔️✔️:                   } else {
    // RDKit✔️✔️:                     return nda1 < nda2;
    // RDKit✔️✔️:                   }
    // RDKit✔️✔️:                 });
    metals.sort_by(|left, right| {
        let left_dative = sanitize_num_dative_bonds(adjacency, assignment, *left);
        let right_dative = sanitize_num_dative_bonds(adjacency, assignment, *right);
        left_dative
            .cmp(&right_dative)
            .then_with(|| ranks[right.index()].cmp(&ranks[left.index()]))
    });
    let metal = metals[0];

    // RDKit✔️✔️:       auto bond =
    // RDKit✔️✔️:           mol.getBondBetweenAtoms(atom->getIdx(), metals.front()->getIdx());
    // RDKit✔️✔️:       if (bond) {
    // RDKit✔️✔️:         bond->setBondType(RDKit::Bond::BondType::DATIVE);
    // RDKit✔️✔️:         bond->setBeginAtom(atom);
    // RDKit✔️✔️:         bond->setEndAtom(metals.front());
    // RDKit✔️✔️:       }
    if let Some(bond_index) = adjacency
        .neighbors_of(atom.index())
        .iter()
        .find(|neighbor| neighbor.atom_index == metal.index())
        .map(|neighbor| neighbor.bond.index())
    {
        assignment.bond_orders[bond_index] = crate::BondOrder::Dative;
        assignment.bond_endpoints[bond_index] = (atom, metal);
    }
    // END RDKIT CPP FUNCTION metalBondCleanup
    Ok(())
}

pub(super) fn sanitize_is_hypervalent_nonmetal<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<bool, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION isHypervalentNonMetal
    // RDKit✔️✔️: bool isHypervalentNonMetal(Atom *atom) {
    // RDKit✔️✔️:   if (QueryOps::isMetal(*atom)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    let atom_data = &molecule.atoms()[atom.index()];
    if is_metal_atomic_number(atom_data.atomic_number()) {
        return Ok(false);
    }
    // RDKit✔️✔️:   int ev = atom->getValence(Atom::ValenceType::EXPLICIT);
    let explicit_valence = valence.explicit_valence[atom.index()];
    // RDKit✔️✔️:   int effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
    let effective_atomic_number =
        i16::from(atom_data.atomic_number()) - i16::from(atom_data.formal_charge());
    // RDKit✔️✔️:   if (effAtomicNum <= 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if effective_atomic_number <= 0 {
        return Ok(false);
    }
    // RDKit✔️✔️:   const auto &otherValens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effAtomicNum);
    let Some(valences) = crate::rdkit_valence_list(effective_atomic_number as u8)? else {
        return Ok(false);
    };
    let Some(&max_valence) = valences.last() else {
        return Ok(false);
    };
    // RDKit✔️✔️:   if (maxV > 0 && (ev > maxV || (ev == maxV && atom->getIsAromatic() &&
    // RDKit✔️✔️:                                  atom->getTotalDegree() == 4))) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if max_valence > 0
        && (explicit_valence > max_valence
            || (explicit_valence == max_valence
                && atom_data.is_aromatic()
                && sanitize_total_degree(molecule, adjacency, valence, atom)? == 4))
    {
        return Ok(true);
    }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isHypervalentNonMetal
    Ok(false)
}

fn sanitize_num_dative_bonds(
    adjacency: &crate::AdjacencyList,
    assignment: &SanitizeOrganometallicCleanupAssignment,
    atom: AtomId,
) -> usize {
    adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| sanitize_is_dative_bond(assignment.bond_orders[neighbor.bond.index()]))
        .count()
}

fn sanitize_no_dative(atom: &crate::Atom) -> bool {
    matches!(atom.atomic_number(), 1 | 2 | 9 | 10)
}

pub(super) fn sanitize_organometallic_single_bonded_metals<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    assignment: &SanitizeOrganometallicCleanupAssignment,
    atom: AtomId,
) -> Vec<AtomId> {
    // RDKit✔️✔️: cleanUpOrganometallics()/metalBondCleanup() single-bonded metal helper
    adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter_map(|neighbor| {
            if assignment.bond_orders[neighbor.bond.index()] != crate::BondOrder::Single {
                return None;
            }
            let other = AtomId::new(neighbor.atom_index);
            is_metal_atomic_number(molecule.atoms()[other.index()].atomic_number()).then_some(other)
        })
        .collect()
}

fn organometallic_cleanup_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &SanitizeOrganometallicCleanupAssignment,
) -> bool {
    let molecule = read_parts;
    molecule
        .bonds()
        .iter()
        .zip(assignment.bond_orders.iter().copied())
        .any(|(bond, order)| bond.order() != order)
        || molecule
            .bonds()
            .iter()
            .zip(assignment.bond_endpoints.iter().copied())
            .any(|(bond, (begin, end))| bond.begin() != begin || bond.end() != end)
}

fn apply_sanitize_organometallic_cleanup_assignment(
    topology: &mut TopologyBlock,
    assignment: SanitizeOrganometallicCleanupAssignment,
) {
    for (bond, (order, (begin, end))) in topology.bonds.iter_mut().zip(
        assignment
            .bond_orders
            .into_iter()
            .zip(assignment.bond_endpoints),
    ) {
        bond.set_order(order);
        bond.set_endpoints(begin, end);
    }
}

fn sanitize_cleanup_atropisomers_assignment(
    read_parts: MoleculeReadParts<'_>,
    rings: &crate::RingInfo,
) -> Vec<(crate::BondId, crate::BondStereo)> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION MolOps::cleanupAtropisomers
    // RDKit✔️✔️: void cleanupAtropisomers(RWMol &mol) {
    // RDKit✔️✔️:   auto hybs = MolOps::Hybridizations(mol);
    // RDKit✔️✔️:   MolOps::cleanupAtropisomers(mol, hybs);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void cleanupAtropisomers(RWMol &mol, MolOps::Hybridizations &hybs) {
    let hybridization = molecule
        .atoms()
        .iter()
        .map(crate::Atom::hybridization)
        .collect::<Vec<_>>();
    let mut assignment = Vec::new();

    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     switch (bond->getStereo()) {
    // RDKit✔️✔️:       case Bond::BondStereo::STEREOATROPCW:
    // RDKit✔️✔️:       case Bond::BondStereo::STEREOATROPCCW:
    for bond in molecule.bonds() {
        if !matches!(
            bond.stereo(),
            crate::BondStereo::AtropCw | crate::BondStereo::AtropCcw
        ) {
            continue;
        }
        // BEGIN RDKIT CPP FUNCTION checkBond
        // RDKit✔️✔️:   if (hybs[bond->getBeginAtomIdx()] != Atom::SP2 ||
        // RDKit✔️✔️:       hybs[bond->getEndAtomIdx()] != Atom::SP2 ||
        // RDKit✔️✔️:       (ri->numBondRings(bond->getIdx()) > 0 &&
        // RDKit✔️✔️:        ri->minBondRingSize(bond->getIdx()) < 8)) {
        // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREONONE);
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION checkBond
        if hybridization[bond.begin().index()] != crate::Hybridization::Sp2
            || hybridization[bond.end().index()] != crate::Hybridization::Sp2
            || (rings.num_bond_rings(bond.id()) > 0 && rings.min_bond_ring_size(bond.id()) < 8)
        {
            assignment.push((bond.id(), crate::BondStereo::None));
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::cleanupAtropisomers
    assignment
}

fn apply_sanitize_cleanup_atropisomers_assignment(
    topology: &mut TopologyBlock,
    assignment: &[(crate::BondId, crate::BondStereo)],
) {
    for &(bond_id, stereo) in assignment {
        let bond = &mut topology.bonds[bond_id.index()];
        bond.set_stereo(stereo);
        if stereo == crate::BondStereo::None {
            bond.set_stereo_atoms(None);
            for stereo_group in &mut topology.stereo_groups {
                stereo_group.remove_bond(bond_id);
            }
        }
    }
    topology.stereo_groups.retain(|group| !group.is_empty());
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SanitizeCleanupChiralityAssignment {
    atom_chiral_tags: Vec<crate::ChiralTag>,
    chiral_permutations: Vec<Option<u32>>,
    cleanup_stereo_groups: bool,
}

fn sanitize_cleanup_chirality_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<SanitizeCleanupChiralityAssignment, crate::ValenceError> {
    let molecule = read_parts;
    let valence = molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?;
    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| crate::ValenceError::UnsupportedBranch {
            reason: "chirality cleanup adjacency could not be built",
        })?;
    let mut assignment = SanitizeCleanupChiralityAssignment {
        atom_chiral_tags: molecule
            .atoms()
            .iter()
            .map(crate::Atom::chiral_tag)
            .collect(),
        chiral_permutations: molecule
            .atoms()
            .iter()
            .map(crate::Atom::chiral_permutation)
            .collect(),
        cleanup_stereo_groups: false,
    };

    // BEGIN RDKIT CPP FUNCTION cleanupChirality
    // RDKit✔️✔️: void cleanupChirality(RWMol &mol) {
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     switch (atom->getChiralTag()) {
    for atom in molecule.atoms() {
        let atom_id = atom.id();
        match atom.chiral_tag() {
            // RDKit✔️✔️:       case Atom::CHI_TETRAHEDRAL_CW:
            // RDKit✔️✔️:       case Atom::CHI_TETRAHEDRAL_CCW:
            // RDKit✔️✔️:         if (atom->getHybridization() != Atom::SP3) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::TetrahedralCw | crate::ChiralTag::TetrahedralCcw => {
                if atom.hybridization() != crate::Hybridization::Sp3 {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                    assignment.cleanup_stereo_groups = true;
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_TETRAHEDRAL:
            // RDKit✔️✔️:         if (atom->getHybridization() != Atom::SP3) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           perm = 0;
            // RDKit✔️✔️:           atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
            // RDKit✔️✔️:           if (perm > 2) {
            // RDKit✔️✔️:             perm = 0;
            // RDKit✔️✔️:             atom->setProp(common_properties::_chiralPermutation, perm);
            crate::ChiralTag::Tetrahedral => {
                if atom.hybridization() != crate::Hybridization::Sp3 {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                    assignment.cleanup_stereo_groups = true;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 2);
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_SQUAREPLANAR:
            // RDKit✔️✔️:         degree = atom->getTotalDegree();
            // RDKit✔️✔️:         if (degree < 2 || degree > 4) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::SquarePlanar => {
                let degree = sanitize_total_degree(molecule, &adjacency, &valence, atom_id)?;
                if !(2..=4).contains(&degree) {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 3);
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_TRIGONALBIPYRAMIDAL:
            // RDKit✔️✔️:         degree = atom->getTotalDegree();
            // RDKit✔️✔️:         if (degree < 2 || degree > 5) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::TrigonalBipyramidal => {
                let degree = sanitize_total_degree(molecule, &adjacency, &valence, atom_id)?;
                if !(2..=5).contains(&degree) {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 20);
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_OCTAHEDRAL:
            // RDKit✔️✔️:         degree = atom->getTotalDegree();
            // RDKit✔️✔️:         if (degree < 2 || degree > 6) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::Octahedral => {
                let degree = sanitize_total_degree(molecule, &adjacency, &valence, atom_id)?;
                if !(2..=6).contains(&degree) {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 30);
                }
            }
            // RDKit✔️✔️:       default:
            // RDKit✔️✔️:         break;
            crate::ChiralTag::Unspecified | crate::ChiralTag::Other | crate::ChiralTag::Allene => {}
        }
    }
    // RDKit✔️✔️:   if (needCleanupStereoGroups) {
    // RDKit✔️✔️:     Chirality::cleanupStereoGroups(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION cleanupChirality
    Ok(assignment)
}

fn sanitize_reset_chiral_permutation_above(
    assignment: &mut SanitizeCleanupChiralityAssignment,
    atom: AtomId,
    limit: u32,
) {
    if assignment.chiral_permutations[atom.index()].is_some_and(|value| value > limit) {
        assignment.chiral_permutations[atom.index()] = Some(0);
    }
}

fn cleanup_chirality_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &SanitizeCleanupChiralityAssignment,
) -> bool {
    let molecule = read_parts;
    molecule
        .atoms()
        .iter()
        .zip(assignment.atom_chiral_tags.iter().copied())
        .any(|(atom, tag)| atom.chiral_tag() != tag)
        || molecule
            .atoms()
            .iter()
            .zip(assignment.chiral_permutations.iter().copied())
            .any(|(atom, permutation)| atom.chiral_permutation() != permutation)
}

fn apply_sanitize_cleanup_chirality_assignment(
    topology: &mut TopologyBlock,
    assignment: &SanitizeCleanupChiralityAssignment,
) {
    let mut cleared_atoms = Vec::new();
    for (idx, (atom, tag)) in topology
        .atoms
        .iter_mut()
        .zip(assignment.atom_chiral_tags.iter().copied())
        .enumerate()
    {
        if atom.chiral_tag() != tag && tag == crate::ChiralTag::Unspecified {
            cleared_atoms.push(AtomId::new(idx));
        }
        atom.set_chiral_tag(tag);
    }
    for (atom, permutation) in topology
        .atoms
        .iter_mut()
        .zip(assignment.chiral_permutations.iter().copied())
    {
        atom.set_chiral_permutation(permutation);
    }
    if assignment.cleanup_stereo_groups {
        for atom in cleared_atoms {
            for stereo_group in &mut topology.stereo_groups {
                stereo_group.remove_atom(atom);
            }
        }
        topology.stereo_groups.retain(|group| !group.is_empty());
    }
}

fn sanitize_adjust_hydrogens_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<u8>, crate::ValenceError> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION adjustHs
    // RDKit✔️✔️: void adjustHs(RWMol &mol) {
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     int origImplicitV = atom->getValence(Atom::ValenceType::IMPLICIT);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     int origExplicitV = atom->getNumExplicitHs();
    // RDKit✔️✔️:     int newImplicitV = atom->calcImplicitValence(false);
    let original_valence = sanitize_valence_facts(molecule)?;
    let current_valence =
        molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?;
    let mut explicit_hydrogens = molecule
        .atoms()
        .iter()
        .map(crate::Atom::explicit_hydrogens)
        .collect::<Vec<_>>();

    for atom in molecule.atoms() {
        let original_implicit = original_valence.implicit_hydrogens[atom.id().index()];
        let new_implicit = current_valence.implicit_hydrogens[atom.id().index()];
        // RDKit✔️✔️:     if (newImplicitV < origImplicitV) {
        // RDKit✔️✔️:       atom->setNumExplicitHs(origExplicitV + (origImplicitV - newImplicitV));
        // RDKit✔️✔️:       atom->calcExplicitValence(false);
        // RDKit✔️✔️:     }
        if new_implicit < original_implicit {
            let delta = u8::try_from(original_implicit - new_implicit).map_err(|_| {
                crate::ValenceError::UnsupportedBranch {
                    reason: "adjustHs implicit hydrogen delta out of range",
                }
            })?;
            explicit_hydrogens[atom.id().index()] =
                explicit_hydrogens[atom.id().index()].saturating_add(delta);
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION adjustHs
    Ok(explicit_hydrogens)
}

fn adjust_hydrogens_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    explicit_hydrogens: &[u8],
) -> bool {
    let molecule = read_parts;
    molecule
        .atoms()
        .iter()
        .zip(explicit_hydrogens.iter().copied())
        .any(|(atom, explicit_hydrogens)| atom.explicit_hydrogens() != explicit_hydrogens)
}

fn apply_sanitize_adjust_hydrogens_assignment(
    topology: &mut TopologyBlock,
    explicit_hydrogens: &[u8],
) {
    for (atom, explicit_hydrogens) in topology.atoms.iter_mut().zip(explicit_hydrogens.iter()) {
        atom.set_explicit_hydrogens(*explicit_hydrogens);
    }
}

fn is_metal_atomic_number(atomic_number: u8) -> bool {
    matches!(
        atomic_number,
        3 | 4
            | 11
            | 12
            | 13
            | 19..=31
            | 37..=50
            | 55..=83
            | 87..=118
    )
}
