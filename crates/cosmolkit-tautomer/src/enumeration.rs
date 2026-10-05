//! Ordered detached enumeration and canonical selection.
use crate::engine::{
    apply_tautomer_transform_match, assign_stereo, canonical_smiles, kekulized, prepared,
    set_tautomer_stereo_and_isotopic_hydrogens, transform_matches,
};
use crate::ordered::*;
use crate::{
    TautomerCatalog, TautomerEnumerationStatus, TautomerParams, TautomerRecord, TautomerRecordView,
    TautomerRunError,
};
use cosmolkit_model::{AtomId, BondId, CoordinateBlock};
use std::{collections::BTreeSet, sync::Arc};

/// Immutable progress at the source's pre-application callback boundary.
pub struct TautomerProgress<'a> {
    state: &'a TautomerExpansionState<Arc<TautomerRecord>>,
}
impl TautomerProgress<'_> {
    pub fn len(&self) -> usize {
        self.state.candidates.len()
    }
    pub fn is_empty(&self) -> bool {
        self.state.candidates.is_empty()
    }
    pub fn status(&self) -> TautomerEnumerationStatus {
        self.state.status
    }
    pub fn num_transforms(&self) -> u32 {
        self.state.num_transforms
    }
    pub fn modified_atoms(&self) -> &BTreeSet<AtomId> {
        &self.state.modified_atoms
    }
    pub fn modified_bonds(&self) -> &BTreeSet<BondId> {
        &self.state.modified_bonds
    }
    pub fn entries(&self) -> impl ExactSizeIterator<Item = (&str, &TautomerRecord)> {
        self.state.candidates.iter().map(|(key, candidate)| {
            (
                key.as_str(),
                candidate
                    .tautomer
                    .as_deref()
                    .expect("enumeration owns materialized candidates"),
            )
        })
    }
}
/// Borrowed cancellation hook. Errors abort the detached run atomically.
pub trait TautomerEnumerationCallback: Send + Sync {
    fn should_continue(
        &self,
        source: TautomerRecordView<'_>,
        progress: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError>;
}

/// Final detached candidates and source metadata, sorted by retained SMILES.
#[derive(Debug, Clone, PartialEq)]
pub struct TautomerEnumerationOutput {
    pub entries: Vec<(String, TautomerRecord)>,
    pub status: TautomerEnumerationStatus,
    pub modified_atoms: BTreeSet<AtomId>,
    pub modified_bonds: BTreeSet<BondId>,
}
impl Default for TautomerEnumerationOutput {
    fn default() -> Self {
        Self {
            entries: Vec::new(),
            status: TautomerEnumerationStatus::Completed,
            modified_atoms: BTreeSet::new(),
            modified_bonds: BTreeSet::new(),
        }
    }
}

fn initialize_candidates_with_key(
    source: TautomerRecordView<'_>,
    key: String,
) -> Result<(TautomerRecord, TautomerExpansionState<Arc<TautomerRecord>>), TautomerRunError> {
    // RDKit✔️❌:   ROMOL_SPTR taut(new ROMol(mol));
    // RDKit✔️❌:   if (taut->needsUpdatePropertyCache()) {
    // RDKit✔️❌:     taut->updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!taut->getRingInfo()->isSymmSssr()) {
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(*taut);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   RWMOL_SPTR kekulized(new RWMol(*taut));
    // RDKit✔️❌:   MolOps::Kekulize(*kekulized, false, true);
    // Owned detached topology/properties replace source ROMol copies; immutable
    // coordinates are borrowed throughout. Existing core owns all preparation.
    let initial = prepared(source)?;
    let initial_kekulized = kekulized(&initial)?;
    let state = TautomerExpansionState {
        candidates: std::collections::BTreeMap::from([(
            key,
            TautomerCandidate {
                tautomer: Some(Arc::new(initial.clone())),
                kekulized: Some(Arc::new(initial_kekulized)),
                num_modified_atoms: 0,
                num_modified_bonds: 0,
                done: false,
            },
        )]),
        modified_atoms: BTreeSet::new(),
        modified_bonds: BTreeSet::new(),
        status: TautomerEnumerationStatus::Completed,
        num_transforms: 0,
    };
    Ok((initial, state))
}
pub fn enumerate_with_catalog(
    source: TautomerRecordView<'_>,
    catalog: &TautomerCatalog,
    params: TautomerParams,
    callback: Option<&dyn TautomerEnumerationCallback>,
) -> Result<TautomerEnumerationOutput, TautomerRunError> {
    // RDKit✔️❌:   std::string smi = MolToSmiles(mol, true);
    let key = canonical_smiles(source)?;
    let (initial, mut state) = initialize_candidates_with_key(source, key)?;
    // Preserve original tags/CIP independently of candidate computed properties.
    let original = &initial;
    // RDKit✔️❌:   while (!completed && !bailOut) {
    loop {
        let pass = expand_tautomer_candidates_in_source_order(
            &mut state,
            catalog.transforms(),
            params,
            |candidate, transform| transform_matches(candidate, source.coordinates, transform),
            |candidate, transform, matched, atoms, bonds, existing| {
                apply_tautomer_transform_match(
                    original,
                    candidate,
                    source.coordinates,
                    transform,
                    matched,
                    atoms,
                    bonds,
                    existing,
                    params,
                )
            },
            |progress| match callback {
                Some(callback) => {
                    callback.should_continue(source, TautomerProgress { state: progress })
                }
                None => Ok(true),
            },
        )
        .map_err(control_expansion_error)?;
        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            params,
            pass.bailed_out,
            |candidate, atoms, bonds| {
                set_tautomer_stereo_and_isotopic_hydrogens(
                    original,
                    Arc::make_mut(candidate),
                    atoms,
                    bonds,
                    params,
                    source.coordinates,
                )
            },
            |candidate| canonical_smiles(candidate.view(source.coordinates)),
        )
        .map_err(control_pruning_error)?;
        if pass.completed || pass.bailed_out {
            break;
        }
    }
    // RDKit✔️❌:   res.fillTautomersItVec();
    // Arc moves preserve source branch sharing during traversal; unwrapping
    // takes final detached ownership without another topology clone.
    let entries = materialize_tautomer_candidates_in_source_order(state.candidates)
        .map_err(|error| TautomerRunError::Control(error.to_string()))?
        .into_iter()
        .map(|(key, value)| (key, Arc::unwrap_or_clone(value)))
        .collect();
    Ok(TautomerEnumerationOutput {
        entries,
        status: state.status,
        modified_atoms: state.modified_atoms,
        modified_bonds: state.modified_bonds,
    })
}
fn control_expansion_error(error: TautomerExpansionError<TautomerRunError>) -> TautomerRunError {
    match error {
        TautomerExpansionError::Backend(error) => error,
        other => TautomerRunError::Control(other.to_string()),
    }
}
fn control_pruning_error(error: TautomerPruningError<TautomerRunError>) -> TautomerRunError {
    match error {
        TautomerPruningError::Backend(error) => error,
        other => TautomerRunError::Control(other.to_string()),
    }
}

/// Select from retained source keys without recomputing SMILES or finalizing values.
pub fn select_canonical_index_with<'a>(
    candidates: impl ExactSizeIterator<Item = (&'a str, TautomerRecordView<'a>)>,
    scorer: impl FnMut(TautomerRecordView<'a>) -> Result<i32, TautomerRunError>,
) -> Result<usize, TautomerRunError> {
    select_canonical_index_by(candidates, scorer, |key, _| {
        Ok(std::borrow::Cow::Borrowed(key))
    })
}
/// Source iterable selection computes canonical keys only for a winning score or tie.
pub fn select_canonical_index_from_iterable_with<'a>(
    candidates: impl ExactSizeIterator<Item = TautomerRecordView<'a>>,
    scorer: impl FnMut(TautomerRecordView<'a>) -> Result<i32, TautomerRunError>,
) -> Result<usize, TautomerRunError> {
    select_canonical_index_by(candidates.map(|value| ((), value)), scorer, |(), view| {
        canonical_smiles(view).map(std::borrow::Cow::Owned)
    })
}
// Retained results and source iterables share exactly one signed-score selection
// loop. Keys are borrowed in the result path and computed lazily in the iterable
// path, preserving scoring order, duplicate inputs and error timing.
fn select_canonical_index_by<'a, E>(
    candidates: impl ExactSizeIterator<Item = (E, TautomerRecordView<'a>)>,
    mut scorer: impl FnMut(TautomerRecordView<'a>) -> Result<i32, TautomerRunError>,
    mut key: impl FnMut(
        E,
        TautomerRecordView<'a>,
    ) -> Result<std::borrow::Cow<'a, str>, TautomerRunError>,
) -> Result<usize, TautomerRunError> {
    // RDKit✔️✔️:   ROMOL_SPTR bestMol;
    // RDKit✔️✔️:   if (tautRes.d_tautomers.size() == 1) {
    // RDKit✔️✔️:     bestMol = tautRes.d_tautomers.begin()->second.tautomer;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // Calculate score for each tautomer
    // RDKit✔️✔️:     int bestScore = std::numeric_limits<int>::min();
    // RDKit✔️✔️:     std::string bestSmiles = "";
    // RDKit✔️✔️:     for (const auto &t : tautRes.d_tautomers) {
    // RDKit✔️✔️:       auto score = scoreFunc(*t.second.tautomer);
    // RDKit✔️✔️:       if (score > bestScore) {
    // RDKit✔️✔️:         bestScore = score;
    // RDKit✔️✔️:         bestSmiles = t.first;
    // RDKit✔️✔️:         bestMol = t.second.tautomer;
    // RDKit✔️✔️:       } else if (score == bestScore) {
    // RDKit✔️✔️:         if (t.first < bestSmiles) {
    // RDKit✔️✔️:           bestSmiles = t.first;
    // RDKit✔️✔️:           bestMol = t.second.tautomer;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:     if (tautomers.size() == 1) {
    // RDKit✔️✔️:       bestMol = *tautomers.begin();
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       for (const auto &t : tautomers) {
    // RDKit✔️✔️:         auto score = scoreFunc(*t);
    // RDKit✔️✔️:         if (score > bestScore) {
    // RDKit✔️✔️:           bestScore = score;
    // RDKit✔️✔️:           bestSmiles = MolToSmiles(*t);
    // RDKit✔️✔️:           bestMol = t;
    // RDKit✔️✔️:         } else if (score == bestScore) {
    // RDKit✔️✔️:           auto smiles = MolToSmiles(*t);
    // RDKit✔️✔️:           if (smiles < bestSmiles) {
    // RDKit✔️✔️:             bestSmiles = smiles;
    // RDKit✔️✔️:             bestMol = t;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    if candidates.len() == 1 {
        return Ok(0);
    }
    let mut best_score = i32::MIN;
    let mut best_smiles = std::borrow::Cow::Borrowed("");
    let mut best = None;
    for (index, (metadata, candidate)) in candidates.enumerate() {
        let score = scorer(candidate)?;
        if score >= best_score {
            let smiles = key(metadata, candidate)?;
            if score > best_score || smiles < best_smiles {
                best_score = score;
                best_smiles = smiles;
                best = Some(index);
            }
        }
    }
    best.ok_or(TautomerRunError::NoCanonicalTautomer)
}
/// Source pickCanonical final copy and forced legacy assignment, without enumeration.
pub fn finalize_canonical_candidate(
    selected: TautomerRecordView<'_>,
) -> Result<TautomerRecord, TautomerRunError> {
    // RDKit✔️❌:   ROMol *res = new ROMol(*bestMol);
    // RDKit✔️❌:   static const bool cleanIt = true;
    // RDKit✔️❌:   static const bool force = true;
    // RDKit✔️❌:   MolOps::assignStereochemistry(*res, cleanIt, force);
    // Legacy assignment owns source-defined missing-cache/ring preparation.
    let mut output = crate::engine::copy_for_canonical_assignment(selected)?;
    assign_stereo(&mut output)?;
    Ok(output)
}
pub fn pick_canonical_with(
    result: &TautomerEnumerationOutput,
    coordinates: &CoordinateBlock,
    scorer: impl FnMut(TautomerRecordView<'_>) -> Result<i32, TautomerRunError>,
) -> Result<TautomerRecord, TautomerRunError> {
    let index = select_canonical_index_with(
        result
            .entries
            .iter()
            .map(|(key, value)| (key.as_str(), value.view(coordinates))),
        scorer,
    )?;
    finalize_canonical_candidate(result.entries[index].1.view(coordinates))
}
pub fn canonicalize_with_catalog(
    source: TautomerRecordView<'_>,
    catalog: &TautomerCatalog,
    params: TautomerParams,
    callback: Option<&dyn TautomerEnumerationCallback>,
    scorer: impl FnMut(TautomerRecordView<'_>) -> Result<i32, TautomerRunError>,
) -> Result<TautomerRecord, TautomerRunError> {
    // RDKit✔️❌:   auto thisCopy = TautomerEnumerator(*this);
    // RDKit✔️❌:   thisCopy.setReassignStereo(false);
    // RDKit✔️❌:   auto res = thisCopy.enumerate(mol);
    let result = enumerate_with_catalog(
        source,
        catalog,
        params.with_reassign_stereo(false),
        callback,
    )?;
    // RDKit✔️❌:   if (res.empty()) {
    // RDKit✔️❌:     return new ROMol(mol);
    // RDKit✔️❌:   }
    if result.entries.is_empty() {
        return prepared(source);
    }
    // RDKit✔️❌:   return pickCanonical(res, scoreFunc);
    pick_canonical_with(&result, source.coordinates, scorer)
}

#[cfg(test)]
mod initialization_tests;
