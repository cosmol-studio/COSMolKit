//! Registered multi-molecule operation body for tautomer enumeration.

use std::{cell::RefCell, collections::BTreeSet};

use super::*;
use crate::chemistry::tautomer::{
    TautomerEnumeration, TautomerEnumerationStatus, TautomerEnumerator, TautomerExpandedProduct,
    TautomerExpansionAttempt, TautomerExpansionError, TautomerExpansionState, TautomerProductPlan,
    TautomerPruningError, TautomerRunError, TautomerStereoUpdatePlan, TautomerTransformAttempt,
    apply_tautomer_transform_match, evaluate_tautomer_callback,
    expand_tautomer_candidates_in_source_order, find_tautomer_transform_matches,
    materialize_tautomer_candidates_in_source_order, plan_tautomer_initialization,
    plan_tautomer_stereo_update, prune_and_rekey_tautomer_candidates_in_source_order,
};

pub(super) struct TautomerEnumerationMetadata {
    canonical_smiles: Vec<String>,
    status: TautomerEnumerationStatus,
    modified_atoms: BTreeSet<AtomId>,
    modified_bonds: BTreeSet<crate::BondId>,
}

fn tautomer_operation_error(source: TautomerRunError) -> OperationError {
    OperationError::Tautomer {
        operation: &ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC,
        source: Box::new(source),
    }
}

fn apply_derived_cache_snapshot(
    branch: &mut OpParts<'_>,
    cache: &crate::molecule::DerivedCacheBlock,
) {
    branch.clear_cache(ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC.needs_update());
    if let Some(rings) = &cache.rings {
        branch.set_rings_cache(rings.clone());
    }
    if let Some(ring_families) = &cache.ring_families {
        branch.set_ring_families_cache(ring_families.clone());
    }
    if let Some(valence) = &cache.valence {
        branch.set_valence_cache(valence.clone());
    }
    if cache.aromaticity_valid {
        branch.mark_aromaticity_valid();
    }
    if cache.stereo_valid {
        branch.mark_stereo_handled();
    }
}

fn apply_topology_and_property_snapshot(
    branch: &mut OpParts<'_>,
    topology: crate::molecule::TopologyBlock,
    properties: MoleculeProperties,
    derived_cache: &crate::molecule::DerivedCacheBlock,
) -> Result<(), OperationError> {
    branch.with_topology_and_properties_mut(|_branch, target_topology, target_properties| {
        *target_topology = topology;
        *target_properties = properties;
        Ok(())
    })?;
    branch.record_topology_edit(TopologyEditKind::Local)?;
    apply_derived_cache_snapshot(branch, derived_cache);
    Ok(())
}

fn apply_kekulize_snapshot(
    branch: &mut OpParts<'_>,
    assignment: &crate::kekulize::KekulizeAssignment,
    derived_cache: &crate::molecule::DerivedCacheBlock,
) -> Result<(), OperationError> {
    branch.with_topology_mut(|_branch, topology| {
        crate::kekulize::apply_kekulize_assignment(topology, assignment);
        Ok(())
    })?;
    branch.record_topology_edit(TopologyEditKind::Local)?;
    apply_derived_cache_snapshot(branch, derived_cache);
    Ok(())
}

fn apply_product_plan(
    branch: &mut OpParts<'_>,
    plan: &TautomerProductPlan,
) -> Result<(), OperationError> {
    apply_topology_and_property_snapshot(
        branch,
        plan.topology.clone(),
        plan.properties.clone(),
        &plan.derived_cache,
    )
}

fn apply_stereo_update_plan(
    branch: &mut OpParts<'_>,
    plan: &TautomerStereoUpdatePlan,
) -> Result<(), OperationError> {
    apply_topology_and_property_snapshot(
        branch,
        plan.topology.clone(),
        plan.properties.clone(),
        &plan.derived_cache,
    )
}

fn map_expansion_error(error: TautomerExpansionError<OperationError>) -> OperationError {
    match error {
        TautomerExpansionError::MissingKekulizedBranch { canonical_smiles } => {
            tautomer_operation_error(TautomerRunError::MissingKekulizedBranch { canonical_smiles })
        }
        TautomerExpansionError::DuplicateProductKey { canonical_smiles } => {
            tautomer_operation_error(TautomerRunError::DuplicateProductKey { canonical_smiles })
        }
        TautomerExpansionError::Backend(source) => source,
    }
}

fn map_pruning_error(error: TautomerPruningError<OperationError>) -> OperationError {
    match error {
        TautomerPruningError::MissingTautomerBranch { canonical_smiles } => {
            tautomer_operation_error(TautomerRunError::MissingTautomerBranch { canonical_smiles })
        }
        TautomerPruningError::Backend(source) => source,
    }
}

#[cosmolkit_macros::mol_multi_op_body(enumerate_tautomers_with_options, parts)]
pub(super) fn enumerate_tautomers_with_options_impl(
    enumerator: &TautomerEnumerator<'_>,
) -> Result<TautomerEnumerationMetadata, OperationError> {
    // RDKit✔️❌: TautomerEnumeratorResult TautomerEnumerator::enumerate(const ROMol &mol) const {
    // RDKit✔️❌:   PRECONDITION(dp_catalog, "no catalog!");
    // RDKit✔️❌:   const TautomerCatalogParams *tautparams = dp_catalog->getCatalogParams();
    // RDKit✔️❌:   PRECONDITION(tautparams, "");
    // RDKit✔️❌:
    // RDKit✔️❌:   TautomerEnumeratorResult res;
    // RDKit✔️❌:
    // RDKit✔️❌:   const std::vector<TautomerTransform> &transforms =
    // RDKit✔️❌:       tautparams->getTransforms();
    let initialization = parts.with_source_read_parts(|read| {
        plan_tautomer_initialization(read)
            .map_err(TautomerRunError::from)
            .map_err(tautomer_operation_error)
    })?;

    let initial_tautomer = parts.derive_from_source(|branch| {
        if initialization.valence_update.is_some() || initialization.rings_update.is_some() {
            apply_derived_cache_snapshot(branch, &initialization.target_derived_cache);
        }
        Ok(())
    })?;
    let initial_kekulized = parts.derive_from_branch(initial_tautomer, |branch| {
        apply_kekulize_snapshot(
            branch,
            &initialization.kekulize_assignment,
            &initialization.target_derived_cache,
        )
    })?;

    let mut state = TautomerExpansionState {
        candidates: initialization.into_candidate_map(initial_tautomer, initial_kekulized),
        modified_atoms: BTreeSet::new(),
        modified_bonds: BTreeSet::new(),
        status: TautomerEnumerationStatus::Completed,
        num_transforms: 0,
    };
    let options = enumerator.options();
    let transforms = enumerator.catalog().transforms();
    let executor = RefCell::new(parts);

    // RDKit✔️❌:   bool completed = false;
    // RDKit✔️❌:   bool bailOut = false;
    // RDKit✔️❌:   unsigned int nTransforms = 0;
    // RDKit✔️❌:   while (!completed && !bailOut) {
    let mut completed = false;
    let mut bailed_out = false;
    while !completed && !bailed_out {
        let expansion = expand_tautomer_candidates_in_source_order(
            &mut state,
            transforms,
            options,
            |kekulized, transform| {
                executor
                    .borrow()
                    .with_branch_read_parts(*kekulized, |read| {
                        find_tautomer_transform_matches(read, transform)
                            .map_err(tautomer_operation_error)
                    })
            },
            |kekulized, transform, matched, modified_atoms, modified_bonds, existing_smiles| {
                let attempt = executor.borrow().with_source_and_branch_read_parts(
                    *kekulized,
                    |source_read, candidate_read| {
                        apply_tautomer_transform_match(
                            source_read,
                            candidate_read,
                            transform,
                            matched,
                            modified_atoms,
                            modified_bonds,
                            existing_smiles,
                            options,
                        )
                        .map_err(TautomerRunError::from)
                        .map_err(tautomer_operation_error)
                    },
                )?;
                match attempt {
                    TautomerTransformAttempt::RecoverableKekulizeFailure {
                        modified_atoms,
                        modified_bonds,
                    } => Ok(TautomerExpansionAttempt::RecoverableKekulizeFailure {
                        modified_atoms,
                        modified_bonds,
                    }),
                    TautomerTransformAttempt::Duplicate {
                        canonical_smiles,
                        modified_atoms,
                        modified_bonds,
                    } => Ok(TautomerExpansionAttempt::Duplicate {
                        canonical_smiles,
                        modified_atoms,
                        modified_bonds,
                    }),
                    TautomerTransformAttempt::Product(plan) => {
                        let canonical_smiles = plan.canonical_smiles.clone();
                        let modified_atoms = plan.modified_atoms.clone();
                        let modified_bonds = plan.modified_bonds.clone();
                        let kekulize_assignment = plan.kekulize_assignment.clone();
                        let derived_cache = plan.derived_cache.clone();
                        let mut executor = executor.borrow_mut();
                        let tautomer = executor.derive_from_branch(*kekulized, |branch| {
                            apply_product_plan(branch, &plan)
                        })?;
                        let kekulized = executor.derive_from_branch(tautomer, |branch| {
                            apply_kekulize_snapshot(branch, &kekulize_assignment, &derived_cache)
                        })?;
                        Ok(TautomerExpansionAttempt::Product(TautomerExpandedProduct {
                            tautomer,
                            kekulized,
                            canonical_smiles,
                            modified_atoms,
                            modified_bonds,
                        }))
                    }
                }
            },
            |current| {
                let Some(callback) = enumerator.callback() else {
                    return Ok(true);
                };
                let ordered = current
                    .candidates
                    .iter()
                    .map(|(canonical_smiles, candidate)| {
                        candidate
                            .tautomer
                            .map(|branch| (canonical_smiles.clone(), branch))
                            .ok_or_else(|| {
                                tautomer_operation_error(TautomerRunError::MissingTautomerBranch {
                                    canonical_smiles: canonical_smiles.clone(),
                                })
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let branches = ordered
                    .iter()
                    .map(|(_, branch)| *branch)
                    .collect::<Vec<_>>();
                executor.borrow().with_source_and_branches_read_parts(
                    &branches,
                    |source_read, branch_reads| {
                        let candidate_reads = ordered
                            .iter()
                            .map(|(canonical_smiles, _)| canonical_smiles.clone())
                            .zip(branch_reads)
                            .collect();
                        evaluate_tautomer_callback(
                            callback,
                            source_read,
                            candidate_reads,
                            current.status,
                            &current.modified_atoms,
                            &current.modified_bonds,
                        )
                        .map_err(tautomer_operation_error)
                    },
                )
            },
        )
        .map_err(map_expansion_error)?;

        let pruning = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            options,
            expansion.bailed_out,
            |tautomer, modified_atoms, modified_bonds| {
                let plan = executor.borrow().with_source_and_branch_read_parts(
                    *tautomer,
                    |source_read, tautomer_read| {
                        plan_tautomer_stereo_update(
                            source_read,
                            tautomer_read,
                            modified_atoms,
                            modified_bonds,
                            options,
                        )
                        .map_err(tautomer_operation_error)
                    },
                )?;
                let changed = plan.changed;
                let replacement = executor
                    .borrow_mut()
                    .derive_from_branch(*tautomer, |branch| {
                        apply_stereo_update_plan(branch, &plan)
                    })?;
                *tautomer = replacement;
                Ok(changed)
            },
            |tautomer| {
                executor.borrow().with_branch_read_parts(*tautomer, |read| {
                    read.canonical_isomeric_smiles()
                        .map_err(TautomerRunError::from)
                        .map_err(tautomer_operation_error)
                })
            },
        )
        .map_err(map_pruning_error)?;
        completed = pruning.completed;
        bailed_out = pruning.bailed_out;
    }

    // RDKit✔️❌:   res.fillTautomersItVec();
    // RDKit✔️❌:   return res;
    let ordered = materialize_tautomer_candidates_in_source_order(state.candidates)
        .map_err(TautomerRunError::from)
        .map_err(tautomer_operation_error)?;
    let mut canonical_smiles = Vec::with_capacity(ordered.len());
    let executor = executor.into_inner();
    for (smiles, branch) in ordered {
        canonical_smiles.push(smiles);
        executor.emit(branch)?;
    }
    Ok(TautomerEnumerationMetadata {
        canonical_smiles,
        status: state.status,
        modified_atoms: state.modified_atoms,
        modified_bonds: state.modified_bonds,
    })
}

pub(super) fn assemble_tautomer_enumeration(
    molecules: Vec<Molecule>,
    metadata: TautomerEnumerationMetadata,
) -> Result<TautomerEnumeration, OperationError> {
    if molecules.len() != metadata.canonical_smiles.len() {
        return Err(OperationError::InvalidInput {
            operation: &ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC,
            message: "tautomer assembler received inconsistent molecule and canonical-key counts",
        });
    }
    Ok(TautomerEnumeration::from_ordered_entries(
        metadata
            .canonical_smiles
            .into_iter()
            .zip(molecules)
            .collect(),
        metadata.status,
        metadata.modified_atoms,
        metadata.modified_bonds,
    ))
}
