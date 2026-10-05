//! Source-ordered expansion, pruning and materialization control.
use crate::{TautomerEnumerationStatus, TautomerParams, TautomerTransform};
use cosmolkit_model::{AtomId, BondId};
use cosmolkit_search::SubstructMatchResult;
use std::collections::{BTreeMap, BTreeSet};
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerCandidate<M> {
    pub(crate) tautomer: Option<M>,
    pub(crate) kekulized: Option<M>,
    pub(crate) num_modified_atoms: usize,
    pub(crate) num_modified_bonds: usize,
    pub(crate) done: bool,
}

pub(crate) type SmilesTautomerMap<M> = BTreeMap<String, TautomerCandidate<M>>;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerExpandedProduct<M> {
    pub(crate) tautomer: M,
    pub(crate) kekulized: M,
    pub(crate) canonical_smiles: String,
    pub(crate) modified_atoms: BTreeSet<AtomId>,
    pub(crate) modified_bonds: BTreeSet<BondId>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum TautomerExpansionAttempt<M> {
    RecoverableKekulizeFailure {
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
    Duplicate {
        canonical_smiles: String,
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
    Product(TautomerExpandedProduct<M>),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TautomerExpansionState<M> {
    pub(crate) candidates: SmilesTautomerMap<M>,
    pub(crate) modified_atoms: BTreeSet<AtomId>,
    pub(crate) modified_bonds: BTreeSet<BondId>,
    pub(crate) status: TautomerEnumerationStatus,
    pub(crate) num_transforms: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TautomerExpansionPass {
    pub(crate) bailed_out: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TautomerPruningPass {
    pub(crate) completed: bool,
    pub(crate) bailed_out: bool,
}

#[derive(Debug, thiserror::Error)]
pub(crate) enum TautomerExpansionError<E> {
    #[error("tautomer candidate {canonical_smiles} has no kekulized branch")]
    MissingKekulizedBranch { canonical_smiles: String },
    #[error("tautomer expansion attempted to replace existing canonical key {canonical_smiles}")]
    DuplicateProductKey { canonical_smiles: String },
    #[error(transparent)]
    Backend(E),
}

#[derive(Debug, thiserror::Error)]
pub(crate) enum TautomerPruningError<E> {
    #[error("tautomer candidate {canonical_smiles} has no materialized branch")]
    MissingTautomerBranch { canonical_smiles: String },
    #[error(transparent)]
    Backend(E),
}

pub(crate) fn expand_tautomer_candidates_in_source_order<M: Clone, E>(
    state: &mut TautomerExpansionState<M>,
    transforms: &[TautomerTransform],
    options: TautomerParams,
    mut find_matches: impl FnMut(&M, &TautomerTransform) -> Result<Vec<SubstructMatchResult>, E>,
    mut apply_match: impl FnMut(
        &M,
        &TautomerTransform,
        &SubstructMatchResult,
        &BTreeSet<AtomId>,
        &BTreeSet<BondId>,
        &dyn Fn(&str) -> bool,
    ) -> Result<TautomerExpansionAttempt<M>, E>,
    mut callback: impl FnMut(&TautomerExpansionState<M>) -> Result<bool, E>,
) -> Result<TautomerExpansionPass, TautomerExpansionError<E>> {
    // RDKit✔️✔️:   bool completed = false;
    // RDKit✔️✔️:   bool bailOut = false;
    // RDKit✔️✔️:   unsigned int nTransforms = 0;
    // `completed` belongs to the post-expansion pruning/rekeying stage. This
    // helper owns one exact ordered expansion pass and retains the other two
    // source variables as its return flag and cumulative state counter.
    let mut bail_out = false;

    // RDKit✔️✔️:   while (!completed && !bailOut) {
    // RDKit✔️✔️:     // std::map automatically sorts res.d_tautomers into alphabetical order
    // RDKit✔️✔️:     // (SMILES)
    // RDKit✔️❌:     for (auto &smilesTautomerPair : res.d_tautomers) {
    // A BTreeMap range cursor, rather than a snapshot of keys, reproduces
    // std::map iterator behavior when a transform inserts during traversal:
    // later keys are visible in this pass and earlier keys wait for the next.
    // Each BTreeMap cursor step seeks in O(log T), unlike the source
    // std::map iterator advance; that traversal overhead remains explicit.
    let mut previous_key: Option<String> = None;
    loop {
        let next_key = match previous_key.as_deref() {
            Some(key) => state
                .candidates
                .range::<str, _>((std::ops::Bound::Excluded(key), std::ops::Bound::Unbounded))
                .next()
                .map(|(key, _)| key.clone()),
            None => state.candidates.keys().next().cloned(),
        };
        let Some(current_key) = next_key else {
            break;
        };
        previous_key = Some(current_key.clone());

        // RDKit✔️✔️:       if (smilesTautomerPair.second.d_done) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        let (done, kekulized) = {
            let candidate = &state.candidates[&current_key];
            (candidate.done, candidate.kekulized.clone())
        };
        if done {
            continue;
        }
        let kekulized =
            kekulized.ok_or_else(|| TautomerExpansionError::MissingKekulizedBranch {
                canonical_smiles: current_key.clone(),
            })?;

        // RDKit✔️✔️:       // tautomer not yet done
        // RDKit✔️✔️:       for (const auto &transform : transforms) {
        for transform in transforms {
            // RDKit✔️✔️:         if (bailOut) {
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            if bail_out {
                break;
            }

            // RDKit✔️✔️:         // kmol is the kekulized version of the tautomer
            // RDKit✔️✔️:         const auto &kmol = smilesTautomerPair.second.kekulized;
            // RDKit✔️✔️:         std::vector<MatchVectType> matches;
            // RDKit✔️✔️:         unsigned int matched =
            // RDKit✔️✔️:             SubstructMatch(*kmol, *(transform.Mol), matches);
            let matches =
                find_matches(&kekulized, transform).map_err(TautomerExpansionError::Backend)?;

            // RDKit✔️✔️:         if (!matched) {
            // RDKit✔️✔️:           continue;
            // RDKit✔️✔️:         }
            if matches.is_empty() {
                continue;
            }

            // RDKit✔️✔️:         ++nTransforms;
            state.num_transforms = state.num_transforms.wrapping_add(1);

            // RDKit✔️✔️:         // loop over transform matches
            // RDKit✔️✔️:         for (const auto &match : matches) {
            for matched in &matches {
                // RDKit✔️✔️:           if (nTransforms >= d_maxTransforms) {
                // RDKit✔️✔️:             res.d_status =
                // RDKit✔️✔️:                 TautomerEnumeratorStatus::MaxTransformsReached;
                // RDKit✔️✔️:             bailOut = true;
                // RDKit✔️✔️:           } else if (res.d_tautomers.size() >= d_maxTautomers) {
                // RDKit✔️✔️:             res.d_status =
                // RDKit✔️✔️:                 TautomerEnumeratorStatus::MaxTautomersReached;
                // RDKit✔️✔️:             bailOut = true;
                // RDKit✔️✔️:           } else if (d_callback.get() &&
                // RDKit✔️✔️:                      !(*d_callback)(mol, res)) {
                // RDKit✔️✔️:             res.d_status = TautomerEnumeratorStatus::Canceled;
                // RDKit✔️✔️:             bailOut = true;
                // RDKit✔️✔️:           }
                if state.num_transforms >= options.max_transforms() {
                    state.status = TautomerEnumerationStatus::MaxTransformsReached;
                    bail_out = true;
                } else if state.candidates.len() >= options.max_tautomers() as usize {
                    state.status = TautomerEnumerationStatus::MaxTautomersReached;
                    bail_out = true;
                } else if !callback(state).map_err(TautomerExpansionError::Backend)? {
                    state.status = TautomerEnumerationStatus::Canceled;
                    bail_out = true;
                }

                // RDKit✔️✔️:           if (bailOut) {
                // RDKit✔️✔️:             break;
                // RDKit✔️✔️:           }
                if bail_out {
                    break;
                }

                // Source res.d_tautomers.find(tsmiles) is one O(log T) lookup.
                // Borrow the retained map through this lookup, without copying keys.
                let contains_smiles = |key: &str| state.candidates.contains_key(key);
                let attempt = apply_match(
                    &kekulized,
                    transform,
                    matched,
                    &state.modified_atoms,
                    &state.modified_bonds,
                    &contains_smiles,
                )
                .map_err(TautomerExpansionError::Backend)?;

                match attempt {
                    TautomerExpansionAttempt::RecoverableKekulizeFailure {
                        modified_atoms,
                        modified_bonds,
                    } => {
                        // RDKit sets the endpoint and directly edited bond bits
                        // before either source `continue` branch.
                        state.modified_atoms = modified_atoms;
                        state.modified_bonds = modified_bonds;
                    }
                    TautomerExpansionAttempt::Duplicate {
                        canonical_smiles: _,
                        modified_atoms,
                        modified_bonds,
                    } => {
                        // RDKit sets the endpoint and directly edited bond bits
                        // before either source `continue` branch.
                        state.modified_atoms = modified_atoms;
                        state.modified_bonds = modified_bonds;
                    }
                    TautomerExpansionAttempt::Product(product) => {
                        state.modified_atoms = product.modified_atoms;
                        state.modified_bonds = product.modified_bonds;
                        let canonical_smiles = product.canonical_smiles;

                        // RDKit✔️✔️:           res.d_tautomers[tsmiles] = Tautomer(
                        // RDKit✔️✔️:               std::move(product),
                        // RDKit✔️✔️:               std::move(kekulized_product),
                        // RDKit✔️✔️:               res.d_modifiedAtoms.count(),
                        // RDKit✔️✔️:               res.d_modifiedBonds.count());
                        if state.candidates.contains_key(&canonical_smiles) {
                            return Err(TautomerExpansionError::DuplicateProductKey {
                                canonical_smiles,
                            });
                        }
                        state.candidates.insert(
                            canonical_smiles,
                            TautomerCandidate {
                                tautomer: Some(product.tautomer),
                                kekulized: Some(product.kekulized),
                                num_modified_atoms: state.modified_atoms.len(),
                                num_modified_bonds: state.modified_bonds.len(),
                                done: false,
                            },
                        );
                    }
                }
            }
        }

        // RDKit✔️✔️:       smilesTautomerPair.second.d_done = true;
        // RDKit✔️✔️:     }
        state
            .candidates
            .get_mut(&current_key)
            .expect("the current ordered-map entry cannot be removed during expansion")
            .mark_done();
    }

    Ok(TautomerExpansionPass {
        bailed_out: bail_out,
    })
}

pub(crate) fn prune_and_rekey_tautomer_candidates_in_source_order<M, E>(
    state: &mut TautomerExpansionState<M>,
    options: TautomerParams,
    mut bail_out: bool,
    mut set_stereo_and_isotopic_hydrogens: impl FnMut(
        &mut M,
        &BTreeSet<AtomId>,
        &BTreeSet<BondId>,
    ) -> Result<bool, E>,
    mut canonical_isomeric_smiles: impl FnMut(&M) -> Result<String, E>,
) -> Result<TautomerPruningPass, TautomerPruningError<E>> {
    // RDKit✔️✔️:     completed = true;
    // RDKit✔️✔️:     size_t maxNumModifiedAtoms = res.d_modifiedAtoms.count();
    // RDKit✔️✔️:     size_t maxNumModifiedBonds = res.d_modifiedBonds.count();
    let mut completed = true;
    let max_num_modified_atoms = state.modified_atoms.len();
    let max_num_modified_bonds = state.modified_bonds.len();

    // RDKit✔️✔️:     for (auto it = res.d_tautomers.begin(); it != res.d_tautomers.end();) {
    let mut current_key = state.candidates.keys().next().cloned();
    while let Some(key) = current_key {
        let (done, num_modified_atoms, num_modified_bonds) = {
            let candidate = &state.candidates[&key];
            (
                candidate.done,
                candidate.num_modified_atoms,
                candidate.num_modified_bonds,
            )
        };

        // RDKit✔️✔️:       auto &taut = it->second;
        // RDKit✔️✔️:       if (!taut.d_done) {
        // RDKit✔️✔️:         completed = false;
        // RDKit✔️✔️:       }
        if !done {
            completed = false;
        }

        // RDKit✔️✔️:       if ((taut.d_numModifiedAtoms < maxNumModifiedAtoms ||
        // RDKit✔️✔️:            taut.d_numModifiedBonds < maxNumModifiedBonds) &&
        // RDKit✔️✔️:           setTautomerStereoAndIsoHs(mol, *taut.tautomer, res)) {
        let needs_stereo_update = num_modified_atoms < max_num_modified_atoms
            || num_modified_bonds < max_num_modified_bonds;
        let stereo_changed = if needs_stereo_update {
            let candidate = state
                .candidates
                .get_mut(&key)
                .expect("the current ordered-map entry exists while pruning");
            let tautomer = candidate.tautomer.as_mut().ok_or_else(|| {
                TautomerPruningError::MissingTautomerBranch {
                    canonical_smiles: key.clone(),
                }
            })?;
            set_stereo_and_isotopic_hydrogens(
                tautomer,
                &state.modified_atoms,
                &state.modified_bonds,
            )
            .map_err(TautomerPruningError::Backend)?
        } else {
            false
        };

        if stereo_changed {
            let new_key = {
                let candidate = &state.candidates[&key];
                let tautomer = candidate.tautomer.as_ref().ok_or_else(|| {
                    TautomerPruningError::MissingTautomerBranch {
                        canonical_smiles: key.clone(),
                    }
                })?;
                canonical_isomeric_smiles(tautomer).map_err(TautomerPruningError::Backend)?
            };
            // RDKit✔️✔️:         Tautomer tautStored = std::move(taut);
            // RDKit✔️✔️:         it = res.d_tautomers.erase(it);
            let mut candidate = state
                .candidates
                .remove(&key)
                .expect("the current ordered-map entry exists while rekeying");
            let next_after_erased = state
                .candidates
                .range::<str, _>((
                    std::ops::Bound::Excluded(key.as_str()),
                    std::ops::Bound::Unbounded,
                ))
                .next()
                .map(|(next_key, _)| next_key.clone());

            // RDKit✔️✔️:         tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
            // RDKit✔️✔️:         tautStored.d_numModifiedBonds = maxNumModifiedBonds;
            candidate.update_modified_counts(max_num_modified_atoms, max_num_modified_bonds);

            // RDKit✔️✔️:         auto insertRes = res.d_tautomers.insert(std::make_pair(
            // RDKit✔️✔️:             MolToSmiles(*tautStored.tautomer), std::move(tautStored)));
            // RDKit✔️✔️:         if (insertRes.second) {
            // RDKit✔️✔️:           it = insertRes.first;
            // RDKit✔️✔️:         }
            if state.candidates.contains_key(&new_key) {
                current_key = next_after_erased;
            } else {
                state.candidates.insert(new_key.clone(), candidate);
                current_key = Some(new_key);
            }
        } else {
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         ++it;
            // RDKit✔️✔️:       }
            current_key = state
                .candidates
                .range::<str, _>((
                    std::ops::Bound::Excluded(key.as_str()),
                    std::ops::Bound::Unbounded,
                ))
                .next()
                .map(|(next_key, _)| next_key.clone());
        }
        // RDKit✔️✔️:     }
    }

    // RDKit✔️✔️:     if (bailOut && res.d_tautomers.size() < d_maxTautomers &&
    // RDKit✔️✔️:         res.d_status == TautomerEnumeratorStatus::MaxTautomersReached) {
    // RDKit✔️✔️:       res.d_status = TautomerEnumeratorStatus::Completed;
    // RDKit✔️✔️:       bailOut = false;
    // RDKit✔️✔️:     }
    if bail_out
        && state.candidates.len() < options.max_tautomers() as usize
        && state.status == TautomerEnumerationStatus::MaxTautomersReached
    {
        state.status = TautomerEnumerationStatus::Completed;
        bail_out = false;
    }

    Ok(TautomerPruningPass {
        completed,
        bailed_out: bail_out,
    })
}

pub(crate) fn materialize_tautomer_candidates_in_source_order<M>(
    candidates: SmilesTautomerMap<M>,
) -> Result<Vec<(String, M)>, TautomerEnumerationError> {
    // RDKit✔️✔️:   res.fillTautomersItVec();
    // RDKit✔️✔️:   void fillTautomersItVec() {
    // RDKit✔️✔️:     for (auto it = d_tautomers.begin(); it != d_tautomers.end(); ++it) {
    // RDKit✔️✔️:       d_tautomersItVec.push_back(it);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut entries = Vec::with_capacity(candidates.len());
    for (canonical_smiles, candidate) in candidates {
        let tautomer = candidate.tautomer.ok_or_else(|| {
            TautomerEnumerationError::MissingCandidateMolecule {
                canonical_smiles: canonical_smiles.clone(),
            }
        })?;
        entries.push((canonical_smiles, tautomer));
    }
    Ok(entries)
}
impl<M> TautomerCandidate<M> {
    fn empty() -> Self {
        // RDKit✔️✔️: Tautomer() : d_numModifiedAtoms(0), d_numModifiedBonds(0), d_done(false) {}
        Self {
            tautomer: None,
            kekulized: None,
            num_modified_atoms: 0,
            num_modified_bonds: 0,
            done: false,
        }
    }

    fn new(
        tautomer: M,
        kekulized: M,
        num_modified_atoms: usize,
        num_modified_bonds: usize,
    ) -> Self {
        // RDKit✔️✔️: Tautomer(ROMOL_SPTR t, ROMOL_SPTR k, size_t a = 0, size_t b = 0)
        // RDKit✔️✔️:     : tautomer(std::move(t)),
        // RDKit✔️✔️:       kekulized(std::move(k)),
        // RDKit✔️✔️:       d_numModifiedAtoms(a),
        // RDKit✔️✔️:       d_numModifiedBonds(b),
        // RDKit✔️✔️:       d_done(false) {}
        Self {
            tautomer: Some(tautomer),
            kekulized: Some(kekulized),
            num_modified_atoms,
            num_modified_bonds,
            done: false,
        }
    }
}

impl<M> TautomerCandidate<M> {
    fn mark_done(&mut self) {
        // RDKit✔️✔️: smilesTautomerPair.second.d_done = true;
        self.done = true;
    }

    fn update_modified_counts(&mut self, num_modified_atoms: usize, num_modified_bonds: usize) {
        // RDKit✔️✔️: tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
        // RDKit✔️✔️: tautStored.d_numModifiedBonds = maxNumModifiedBonds;
        self.num_modified_atoms = num_modified_atoms;
        self.num_modified_bonds = num_modified_bonds;
    }
}
#[derive(Debug, thiserror::Error)]
pub(crate) enum TautomerEnumerationError {
    #[error("tautomer candidate {canonical_smiles} has no materialized value")]
    MissingCandidateMolecule { canonical_smiles: String },
}

#[cfg(test)]
mod tests {
    use super::*;
    fn query(text: &str) -> cosmolkit_model::QueryGraph {
        cosmolkit_search::parse_smarts(text, &Default::default()).expect("parse fixture")
    }
    fn record(
        text: &str,
    ) -> Result<cosmolkit_smiles::SmilesRecord, cosmolkit_smiles::SmilesParseError> {
        cosmolkit_smiles::parse_smiles(text, &Default::default())
    }
    fn marked_atoms(indices: impl IntoIterator<Item = usize>) -> BTreeSet<AtomId> {
        indices.into_iter().map(AtomId::new).collect()
    }
    fn marked_bonds(indices: impl IntoIterator<Item = usize>) -> BTreeSet<BondId> {
        indices.into_iter().map(BondId::new).collect()
    }
    fn expansion_transform(name: &str) -> TautomerTransform {
        TautomerTransform::new(name, query("[C]-[C]"), Vec::new(), Vec::new())
            .expect("construct expansion-control transform")
    }

    fn expansion_candidate(handle: &str, done: bool) -> TautomerCandidate<String> {
        TautomerCandidate {
            tautomer: Some(format!("{handle}-tautomer")),
            kekulized: Some(handle.to_owned()),
            num_modified_atoms: 0,
            num_modified_bonds: 0,
            done,
        }
    }

    fn expansion_state(entries: &[(&str, &str, bool)]) -> TautomerExpansionState<String> {
        TautomerExpansionState {
            candidates: entries
                .iter()
                .map(|&(key, handle, done)| (key.to_owned(), expansion_candidate(handle, done)))
                .collect(),
            modified_atoms: BTreeSet::new(),
            modified_bonds: BTreeSet::new(),
            status: TautomerEnumerationStatus::Completed,
            num_transforms: 0,
        }
    }

    fn expansion_match(tag: usize) -> SubstructMatchResult {
        SubstructMatchResult {
            atom_mapping: vec![tag, tag],
            bond_mapping: vec![tag],
        }
    }

    fn expansion_product(key: &str, tag: usize) -> TautomerExpansionAttempt<String> {
        TautomerExpansionAttempt::Product(TautomerExpandedProduct {
            tautomer: format!("{key}-tautomer"),
            kekulized: key.to_owned(),
            canonical_smiles: key.to_owned(),
            modified_atoms: BTreeSet::from([AtomId::new(tag)]),
            modified_bonds: BTreeSet::from([BondId::new(tag)]),
        })
    }

    #[test]
    fn enumeration_expansion_visits_later_insertions_in_pass_and_defers_earlier_insertions() {
        let transforms = [expansion_transform("ordered")];
        let mut state = expansion_state(&[("m", "m", false)]);
        let mut visited = Vec::new();

        let first = expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerParams::default(),
            |handle, _| {
                visited.push(handle.clone());
                Ok::<_, &'static str>(match handle.as_str() {
                    "m" | "z" => vec![expansion_match(0)],
                    _ => Vec::new(),
                })
            },
            |handle, _, _, _, _, _| {
                Ok::<_, &'static str>(match handle.as_str() {
                    "m" => expansion_product("z", 0),
                    "z" => expansion_product("a", 1),
                    _ => unreachable!("only matching handles are applied"),
                })
            },
            |_| Ok(true),
        )
        .expect("first ordered expansion pass");

        assert!(!first.bailed_out);
        assert_eq!(visited, ["m", "z"]);
        assert_eq!(
            state
                .candidates
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["a", "m", "z"]
        );
        assert!(!state.candidates["a"].done);
        assert!(state.candidates["m"].done);
        assert!(state.candidates["z"].done);

        visited.clear();
        let second = expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerParams::default(),
            |handle, _| {
                visited.push(handle.clone());
                Ok::<_, &'static str>(Vec::new())
            },
            |_, _, _, _, _, _| -> Result<_, &'static str> {
                unreachable!("the deferred key has no match")
            },
            |_| Ok(true),
        )
        .expect("second ordered expansion pass");
        assert!(!second.bailed_out);
        assert_eq!(visited, ["a"]);
        assert!(state.candidates.values().all(|candidate| candidate.done));
    }

    #[test]
    fn enumeration_expansion_counts_one_transform_for_multiple_and_duplicate_matches() {
        let transforms = [expansion_transform("multi-match")];
        let mut state = expansion_state(&[("m", "m", false)]);
        let mut applied = 0;
        let mut callback_sizes = Vec::new();

        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerParams::default(),
            |handle, _| {
                Ok::<_, &'static str>(if handle == "m" {
                    vec![expansion_match(0), expansion_match(1)]
                } else {
                    Vec::new()
                })
            },
            |_, _, matched, modified_atoms, modified_bonds, existing| {
                applied += 1;
                if matched.bond_mapping[0] == 0 {
                    Ok::<_, &'static str>(expansion_product("z", 0))
                } else {
                    assert!(existing("z"));
                    Ok(TautomerExpansionAttempt::Duplicate {
                        canonical_smiles: "z".to_owned(),
                        modified_atoms: modified_atoms
                            .union(&BTreeSet::from([AtomId::new(1)]))
                            .copied()
                            .collect(),
                        modified_bonds: modified_bonds
                            .union(&BTreeSet::from([BondId::new(1)]))
                            .copied()
                            .collect(),
                    })
                }
            },
            |view| {
                callback_sizes.push((view.num_transforms, view.candidates.len()));
                Ok(true)
            },
        )
        .expect("expand multiple matches");

        assert_eq!(state.num_transforms, 1);
        assert_eq!(applied, 2);
        assert_eq!(callback_sizes, [(1, 1), (1, 2)]);
        assert_eq!(state.modified_atoms, marked_atoms([0, 1]));
        assert_eq!(state.modified_bonds, marked_bonds([0, 1]));
    }

    #[test]
    fn enumeration_expansion_zero_and_exact_transform_limits_use_source_increment_order() {
        for limit in [0, 1] {
            let transforms = [expansion_transform("limited")];
            let mut state = expansion_state(&[("m", "m", false)]);
            let mut callbacks = 0;
            let mut applications = 0;
            let pass = expand_tautomer_candidates_in_source_order(
                &mut state,
                &transforms,
                TautomerParams::default().with_max_transforms(limit),
                |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
                |_, _, _, _, _, _| {
                    applications += 1;
                    Ok::<_, &'static str>(expansion_product("z", 0))
                },
                |_| {
                    callbacks += 1;
                    Ok(true)
                },
            )
            .expect("apply transform limit");

            assert!(pass.bailed_out, "limit {limit}");
            assert_eq!(state.num_transforms, 1, "limit {limit}");
            assert_eq!(
                state.status,
                TautomerEnumerationStatus::MaxTransformsReached,
                "limit {limit}"
            );
            assert_eq!(callbacks, 0, "limit {limit}");
            assert_eq!(applications, 0, "limit {limit}");
            assert!(state.candidates["m"].done, "limit {limit}");
        }

        let transforms = [expansion_transform("first"), expansion_transform("second")];
        let mut state = expansion_state(&[("m", "m", false)]);
        let mut applied = Vec::new();
        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerParams::default().with_max_transforms(2),
            |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
            |_, transform, _, _, _, _| {
                applied.push(transform.name().to_owned());
                Ok::<_, &'static str>(TautomerExpansionAttempt::Duplicate {
                    canonical_smiles: "m".to_owned(),
                    modified_atoms: BTreeSet::new(),
                    modified_bonds: BTreeSet::new(),
                })
            },
            |_| Ok(true),
        )
        .expect("exact transform boundary");
        assert_eq!(applied, ["first"]);
        assert_eq!(state.num_transforms, 2);
        assert_eq!(
            state.status,
            TautomerEnumerationStatus::MaxTransformsReached
        );
    }

    #[test]
    fn enumeration_expansion_tautomer_limits_are_checked_before_callback_and_each_match() {
        let transforms = [expansion_transform("two matches")];
        for limit in [0, 1] {
            let mut state = expansion_state(&[("m", "m", false)]);
            let mut callbacks = 0;
            let pass = expand_tautomer_candidates_in_source_order(
                &mut state,
                &transforms,
                TautomerParams::default().with_max_tautomers(limit),
                |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
                |_, _, _, _, _, _| Ok::<_, &'static str>(expansion_product("z", 0)),
                |_| {
                    callbacks += 1;
                    Ok(true)
                },
            )
            .expect("apply immediate tautomer limit");
            assert!(pass.bailed_out, "limit {limit}");
            assert_eq!(
                state.status,
                TautomerEnumerationStatus::MaxTautomersReached,
                "limit {limit}"
            );
            assert_eq!(callbacks, 0, "limit {limit}");
        }

        let mut state = expansion_state(&[("m", "m", false)]);
        let mut callbacks = Vec::new();
        let mut applications = 0;
        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerParams::default().with_max_tautomers(2),
            |_, _| Ok::<_, &'static str>(vec![expansion_match(0), expansion_match(1)]),
            |_, _, _, _, _, _| {
                applications += 1;
                Ok::<_, &'static str>(expansion_product("z", 0))
            },
            |view| {
                callbacks.push(view.candidates.len());
                Ok(true)
            },
        )
        .expect("apply per-match tautomer limit");
        assert_eq!(callbacks, [1]);
        assert_eq!(applications, 1);
        assert_eq!(state.status, TautomerEnumerationStatus::MaxTautomersReached);
    }

    #[test]
    fn enumeration_expansion_callback_observes_preapplication_state_and_cancels_deterministically()
    {
        let transforms = [expansion_transform("callback")];
        let run = || {
            let mut state = expansion_state(&[("m", "m", false)]);
            let mut observations = Vec::new();
            let mut applications = 0;
            let pass = expand_tautomer_candidates_in_source_order(
                &mut state,
                &transforms,
                TautomerParams::default(),
                |_, _| Ok::<_, &'static str>(vec![expansion_match(0)]),
                |_, _, _, _, _, _| {
                    applications += 1;
                    Ok::<_, &'static str>(expansion_product("z", 0))
                },
                |view| {
                    observations.push((
                        view.num_transforms,
                        view.candidates.keys().cloned().collect::<Vec<_>>(),
                        view.status,
                    ));
                    Ok(false)
                },
            )
            .expect("callback cancellation");
            (state, pass, observations, applications)
        };

        let first = run();
        let second = run();
        assert_eq!(first, second);
        assert!(first.1.bailed_out);
        assert_eq!(first.0.status, TautomerEnumerationStatus::Canceled);
        assert_eq!(
            first.2,
            [(
                1,
                vec!["m".to_owned()],
                TautomerEnumerationStatus::Completed
            )]
        );
        assert_eq!(first.3, 0);
        assert!(first.0.candidates["m"].done);
    }

    #[test]
    fn enumeration_expansion_skips_done_candidates_and_marks_later_keys_done_after_bailout() {
        let transforms = [expansion_transform("stop")];
        let mut state = expansion_state(&[("a", "a", true), ("m", "m", false), ("z", "z", false)]);
        let mut matched_handles = Vec::new();
        expand_tautomer_candidates_in_source_order(
            &mut state,
            &transforms,
            TautomerParams::default().with_max_transforms(1),
            |handle, _| {
                matched_handles.push(handle.clone());
                Ok::<_, &'static str>(vec![expansion_match(0)])
            },
            |_, _, _, _, _, _| -> Result<_, &'static str> {
                unreachable!("limit is checked before application")
            },
            |_| Ok(true),
        )
        .expect("bail out in ordered traversal");

        assert_eq!(matched_handles, ["m"]);
        assert!(state.candidates["a"].done);
        assert!(state.candidates["m"].done);
        assert!(state.candidates["z"].done);
    }

    #[test]
    fn enumeration_pruning_rekeys_only_when_stereo_changes_and_tracks_completion() {
        let mut state = expansion_state(&[("old", "old", true), ("stable", "stable", false)]);
        state.modified_atoms = marked_atoms([0]);
        let mut stereo_calls = Vec::new();
        let mut key_calls = 0;

        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerParams::default(),
            false,
            |tautomer, modified_atoms, modified_bonds| {
                stereo_calls.push(tautomer.clone());
                assert_eq!(modified_atoms, &marked_atoms([0]));
                assert!(modified_bonds.is_empty());
                if tautomer == "old-tautomer" {
                    *tautomer = "new".to_owned();
                    Ok::<_, &'static str>(true)
                } else {
                    Ok(false)
                }
            },
            |tautomer| {
                key_calls += 1;
                Ok::<_, &'static str>(tautomer.clone())
            },
        )
        .expect("prune changed and unchanged stereo branches");

        assert!(!pass.completed);
        assert!(!pass.bailed_out);
        assert_eq!(stereo_calls, ["old-tautomer", "stable-tautomer"]);
        assert_eq!(key_calls, 1);
        assert!(!state.candidates.contains_key("old"));
        assert_eq!(state.candidates["new"].num_modified_atoms, 1);
        assert_eq!(state.candidates["new"].num_modified_bonds, 0);
        assert_eq!(state.candidates["stable"].num_modified_atoms, 0);
        assert_eq!(state.candidates["stable"].num_modified_bonds, 0);
    }

    #[test]
    fn enumeration_pruning_reproduces_source_rekey_iterator_order() {
        let mut forward = expansion_state(&[("a", "a", true), ("b", "b", true), ("c", "c", true)]);
        forward.modified_atoms = marked_atoms([0]);
        let mut forward_calls = Vec::new();
        let forward_pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut forward,
            TautomerParams::default(),
            false,
            |tautomer, _, _| {
                forward_calls.push(tautomer.clone());
                if tautomer == "a-tautomer" {
                    *tautomer = "z".to_owned();
                    Ok::<_, &'static str>(true)
                } else {
                    Ok(false)
                }
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("forward rekey traversal");
        assert!(forward_pass.completed);
        assert_eq!(forward_calls, ["a-tautomer"]);
        assert_eq!(
            forward
                .candidates
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["b", "c", "z"]
        );

        let mut backward = expansion_state(&[("a", "a", true), ("b", "b", true), ("c", "c", true)]);
        backward.modified_atoms = marked_atoms([0]);
        let mut backward_calls = Vec::new();
        let backward_pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut backward,
            TautomerParams::default(),
            false,
            |tautomer, _, _| {
                backward_calls.push(tautomer.clone());
                if tautomer == "c-tautomer" {
                    *tautomer = "aa".to_owned();
                    Ok::<_, &'static str>(true)
                } else {
                    Ok(false)
                }
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("backward rekey traversal");
        assert!(backward_pass.completed);
        assert_eq!(
            backward_calls,
            ["a-tautomer", "b-tautomer", "c-tautomer", "b-tautomer"]
        );
        assert_eq!(
            backward
                .candidates
                .keys()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["a", "aa", "b"]
        );
    }

    #[test]
    fn enumeration_pruning_collapses_duplicates_and_corrects_only_tautomer_limit_status() {
        let mut state = expansion_state(&[("a", "a", true), ("b", "b", true)]);
        state.modified_atoms = marked_atoms([0]);
        state.candidates.get_mut("a").unwrap().num_modified_atoms = 1;
        state.status = TautomerEnumerationStatus::MaxTautomersReached;

        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerParams::default().with_max_tautomers(2),
            true,
            |tautomer, _, _| {
                assert_eq!(tautomer, "b-tautomer");
                *tautomer = "a".to_owned();
                Ok::<_, &'static str>(true)
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("collapse duplicate and correct status");

        assert!(pass.completed);
        assert!(!pass.bailed_out);
        assert_eq!(state.status, TautomerEnumerationStatus::Completed);
        assert_eq!(state.candidates.len(), 1);
        assert_eq!(
            state.candidates["a"].tautomer.as_deref(),
            Some("a-tautomer")
        );

        let mut transform_limited = expansion_state(&[("a", "a", true), ("b", "b", true)]);
        transform_limited.modified_atoms = marked_atoms([0]);
        transform_limited
            .candidates
            .get_mut("a")
            .unwrap()
            .num_modified_atoms = 1;
        transform_limited.status = TautomerEnumerationStatus::MaxTransformsReached;
        let pass = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut transform_limited,
            TautomerParams::default().with_max_tautomers(2),
            true,
            |tautomer, _, _| {
                *tautomer = "a".to_owned();
                Ok::<_, &'static str>(true)
            },
            |tautomer| Ok::<_, &'static str>(tautomer.clone()),
        )
        .expect("retain non-tautomer limit status");
        assert!(pass.bailed_out);
        assert_eq!(
            transform_limited.status,
            TautomerEnumerationStatus::MaxTransformsReached
        );
    }

    #[test]
    fn enumeration_pruning_reapplies_after_modified_sets_grow_across_rounds() {
        let mut state = expansion_state(&[("m", "m", true)]);
        state.modified_atoms = marked_atoms([0]);
        state.modified_bonds = marked_bonds([0]);
        {
            let candidate = state.candidates.get_mut("m").unwrap();
            candidate.num_modified_atoms = 1;
            candidate.num_modified_bonds = 1;
        }

        let first = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerParams::default(),
            false,
            |_, _, _| -> Result<_, &'static str> {
                unreachable!("equal modified-set snapshots do not reapply stereo")
            },
            |_| -> Result<_, &'static str> {
                unreachable!("an unchanged candidate is not rekeyed")
            },
        )
        .expect("prune with unchanged modified sets");
        assert!(first.completed);

        state.modified_atoms.insert(AtomId::new(1));
        state.modified_bonds.insert(BondId::new(1));
        let mut calls = 0;
        let second = prune_and_rekey_tautomer_candidates_in_source_order(
            &mut state,
            TautomerParams::default(),
            false,
            |tautomer, modified_atoms, modified_bonds| {
                calls += 1;
                assert_eq!(modified_atoms, &marked_atoms([0, 1]));
                assert_eq!(modified_bonds, &marked_bonds([0, 1]));
                tautomer.push_str("-updated");
                Ok::<_, &'static str>(true)
            },
            |_| Ok::<_, &'static str>("m".to_owned()),
        )
        .expect("prune after modified sets grow");

        assert!(second.completed);
        assert_eq!(calls, 1);
        assert_eq!(state.candidates["m"].num_modified_atoms, 2);
        assert_eq!(state.candidates["m"].num_modified_bonds, 2);
        assert_eq!(
            state.candidates["m"].tautomer.as_deref(),
            Some("m-tautomer-updated")
        );
    }

    #[test]
    fn enumeration_pruning_materializes_final_candidates_in_key_order() {
        let state = expansion_state(&[("z", "z", true), ("a", "a", true), ("m", "m", true)]);
        let entries = materialize_tautomer_candidates_in_source_order(state.candidates)
            .expect("materialize ordered candidates");
        assert_eq!(
            entries,
            [
                ("a".to_owned(), "a-tautomer".to_owned()),
                ("m".to_owned(), "m-tautomer".to_owned()),
                ("z".to_owned(), "z-tautomer".to_owned()),
            ]
        );

        let mut missing = expansion_state(&[("missing", "missing", true)]).candidates;
        missing.get_mut("missing").unwrap().tautomer = None;
        assert!(matches!(
            materialize_tautomer_candidates_in_source_order(missing),
            Err(TautomerEnumerationError::MissingCandidateMolecule { canonical_smiles })
                if canonical_smiles == "missing"
        ));
    }

    #[test]
    fn candidate_record_default_constructor_has_source_empty_state() {
        let candidate = TautomerCandidate::<String>::empty();

        assert!(candidate.tautomer.is_none());
        assert!(candidate.kekulized.is_none());
        assert_eq!(candidate.num_modified_atoms, 0);
        assert_eq!(candidate.num_modified_bonds, 0);
        assert!(!candidate.done);
    }

    #[test]
    fn candidate_record_explicit_constructor_preserves_values_and_initial_state() {
        let tautomer = record("CCO").expect("parse tautomer");
        let kekulized = record("CCO").expect("parse kekulized tautomer");
        let candidate = TautomerCandidate::new(tautomer.clone(), kekulized.clone(), 3, 2);

        assert_eq!(candidate.tautomer.as_ref(), Some(&tautomer));
        assert_eq!(candidate.kekulized.as_ref(), Some(&kekulized));
        assert_eq!(candidate.num_modified_atoms, 3);
        assert_eq!(candidate.num_modified_bonds, 2);
        assert!(!candidate.done);
    }

    #[test]
    fn candidate_record_clone_and_state_transitions_are_independent() {
        let molecule = record("CC=O").expect("parse candidate");
        let original = TautomerCandidate::new(molecule.clone(), molecule, 1, 1);
        let mut changed = original.clone();

        changed.update_modified_counts(4, 5);
        changed.mark_done();

        assert_eq!(original.num_modified_atoms, 1);
        assert_eq!(original.num_modified_bonds, 1);
        assert!(!original.done);
        assert_eq!(changed.num_modified_atoms, 4);
        assert_eq!(changed.num_modified_bonds, 5);
        assert!(changed.done);
    }

    #[test]
    fn candidate_record_ordered_map_sorts_keys_and_replaces_one_canonical_key() {
        let molecule = record("C").expect("parse candidate");
        let candidate = |modified_atoms| {
            TautomerCandidate::new(molecule.clone(), molecule.clone(), modified_atoms, 0)
        };
        let mut candidates = SmilesTautomerMap::new();
        candidates.insert("z".to_owned(), candidate(1));
        candidates.insert("a".to_owned(), candidate(2));
        candidates.insert("m".to_owned(), candidate(3));
        let replaced = candidates.insert("a".to_owned(), candidate(9));

        assert_eq!(replaced.expect("existing key").num_modified_atoms, 2);
        assert_eq!(
            candidates.keys().map(String::as_str).collect::<Vec<_>>(),
            ["a", "m", "z"]
        );
        assert_eq!(candidates["a"].num_modified_atoms, 9);
    }
}
