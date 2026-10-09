//! Ordered detached enumeration and canonical selection.
use crate::engine::{
    apply_tautomer_transform_match, assign_stereo, canonical_smiles, get_cached_kekulized,
    prepared, set_tautomer_stereo_and_isotopic_hydrogens, transform_matches,
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
    state: &'a mut TautomerExpansionState<Arc<TautomerRecord>>,
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
    pub fn entries_readonly(
        &self,
    ) -> impl ExactSizeIterator<Item = (&cosmolkit_model::PropertyText, &TautomerRecord)> {
        self.state.candidates.iter().map(|(key, candidate)| {
            (
                key,
                candidate
                    .tautomer
                    .as_deref()
                    .expect("enumeration owns materialized candidates"),
            )
        })
    }
    pub fn entries_mut<'a>(
        &'a mut self,
        coordinates: &'a CoordinateBlock,
    ) -> impl ExactSizeIterator<
        Item = (
            &'a cosmolkit_model::PropertyText,
            crate::TautomerScoreView<'a>,
        ),
    > + 'a {
        self.state
            .candidates
            .iter_mut()
            .map(move |(key, candidate)| {
                let record = Arc::make_mut(
                    candidate
                        .tautomer
                        .as_mut()
                        .expect("enumeration owns materialized candidates"),
                );
                (&*key, record.score_view(coordinates))
            })
    }
}
/// Borrowed cancellation hook. Errors abort the detached run atomically.
pub trait TautomerEnumerationCallback {
    fn should_continue(
        &mut self,
        source: crate::TautomerScoreView<'_>,
        progress: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError>;
}

/// Final detached candidates and source metadata, sorted by retained SMILES.
#[derive(Debug, Clone, PartialEq)]
pub struct TautomerEnumerationOutput {
    pub entries: Vec<(cosmolkit_model::PropertyText, TautomerRecord)>,
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
    key: cosmolkit_model::PropertyText,
) -> Result<(TautomerRecord, TautomerExpansionState<Arc<TautomerRecord>>), TautomerRunError> {
    // RDKit✔️❌:   ROMOL_SPTR taut(new ROMol(mol));
    // RDKit✔️❌:   if (taut->needsUpdatePropertyCache()) {
    // RDKit✔️❌:     taut->updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!taut->getRingInfo()->isSymmSssr()) {
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(*taut);
    // RDKit✔️❌:   }
    // Owned detached topology/properties replace source ROMol copies; immutable
    // coordinates are borrowed throughout. Existing core owns all preparation.
    let initial = prepared(source)?;
    let state = TautomerExpansionState {
        candidates: std::collections::BTreeMap::from([(
            key,
            TautomerCandidate {
                tautomer: Some(Arc::new(initial.clone())),
                kekulized: None,
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
    mut source: crate::TautomerScoreView<'_>,
    catalog: &TautomerCatalog,
    params: TautomerParams,
    mut callback: Option<&mut dyn TautomerEnumerationCallback>,
) -> Result<TautomerEnumerationOutput, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerEnumerator::enumerate actual source cache loan
    // RDKit❗❌: TautomerEnumeratorResult TautomerEnumerator::enumerate(const ROMol &mol) const {
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:   std::cout << "**********************************" << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:   PRECONDITION(dp_catalog, "no catalog!");
    // RDKit❗❌:   const TautomerCatalogParams *tautparams = dp_catalog->getCatalogParams();
    // RDKit❗❌:   PRECONDITION(tautparams, "");
    // RDKit❗❌:
    // RDKit❗❌:   TautomerEnumeratorResult res;
    // RDKit❗❌:
    // RDKit❗❌:   const std::vector<TautomerTransform> &transforms =
    // RDKit❗❌:       tautparams->getTransforms();
    // RDKit❗❌:
    // RDKit❗❌:   // Enumerate all possible tautomers and return them as a vector.
    // RDKit❗❌:   // smi is the input molecule SMILES
    // RDKit❗❌:   std::string smi = MolToSmiles(mol, true);
    // RDKit❗❌:   // taut is a copy of the input molecule
    // RDKit❗❌:   ROMOL_SPTR taut(new ROMol(mol));
    // RDKit❗❌:   // do whatever sanitization bits are required
    // RDKit❗❌:   if (taut->needsUpdatePropertyCache()) {
    // RDKit❗❌:     taut->updatePropertyCache(false);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!taut->getRingInfo()->isSymmSssr()) {
    // RDKit❗❌:     MolOps::symmetrizeSSSR(*taut);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // Kekulized form will be created lazily when needed for transform matching.
    // RDKit❗❌:   // canonical=true is used on demand so that tautomer deduplication is
    // RDKit❗❌:   // independent of atom ordering in the molecule.
    // RDKit❗❌:   res.d_tautomers = {{smi, Tautomer(taut, 0, 0)}};
    // RDKit❗❌:   res.d_modifiedAtoms.resize(mol.getNumAtoms());
    // RDKit❗❌:   res.d_modifiedBonds.resize(mol.getNumBonds());
    // RDKit❗❌:
    // RDKit❗❌:   // Keep running counts of modified atoms/bonds.
    // RDKit❗❌:   // `boost::dynamic_bitset<>::count()` is O(n) in the number of blocks, and we
    // RDKit❗❌:   // were previously calling it once per new tautomer, which is avoidable.
    // RDKit❗❌:   size_t numModifiedAtoms = 0;
    // RDKit❗❌:   size_t numModifiedBonds = 0;
    // RDKit❗❌:   const auto markAtomModified = [&res, &numModifiedAtoms](unsigned int idx) {
    // RDKit❗❌:     if (!res.d_modifiedAtoms.test(idx)) {
    // RDKit❗❌:       res.d_modifiedAtoms.set(idx);
    // RDKit❗❌:       ++numModifiedAtoms;
    // RDKit❗❌:     }
    // RDKit❗❌:   };
    // RDKit❗❌:   const auto markBondModified = [&res, &numModifiedBonds](unsigned int idx) {
    // RDKit❗❌:     if (!res.d_modifiedBonds.test(idx)) {
    // RDKit❗❌:       res.d_modifiedBonds.set(idx);
    // RDKit❗❌:       ++numModifiedBonds;
    // RDKit❗❌:     }
    // RDKit❗❌:   };
    // RDKit❗❌:   bool completed = false;
    // RDKit❗❌:   bool bailOut = false;
    // RDKit❗❌:   unsigned int nTransforms = 0;
    // RDKit❗❌:   static const std::array<const char *, 4> statusMsg{
    // RDKit❗❌:       "completed", "max tautomers reached", "max transforms reached",
    // RDKit❗❌:       "canceled"};
    // RDKit❗❌:
    // RDKit❗❌:   while (!completed && !bailOut) {
    // RDKit❗❌:     // std::map automatically sorts res.d_tautomers into alphabetical order
    // RDKit❗❌:     // (SMILES)
    // RDKit❗❌:     for (auto &smilesTautomerPair : res.d_tautomers) {
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:       std::cout << "Current tautomers: " << std::endl;
    // RDKit❗❌:       for (const auto &smilesTautomerPair : res.d_tautomers) {
    // RDKit❗❌:         std::cout << smilesTautomerPair.first << " done "
    // RDKit❗❌:                   << smilesTautomerPair.second.d_done << std::endl;
    // RDKit❗❌:       }
    // RDKit❗❌: #endif
    // RDKit❗❌:       std::string tsmiles;
    // RDKit❗❌:       if (smilesTautomerPair.second.d_done) {
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:         std::cout << "Skipping " << smilesTautomerPair.first
    // RDKit❗❌:                   << " as already done" << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:       std::cout << "Looking at tautomer: " << smilesTautomerPair.first
    // RDKit❗❌:                 << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:       // tautomer not yet done
    // RDKit❗❌:       for (const auto &transform : transforms) {
    // RDKit❗❌:         if (bailOut) {
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:         // kmol is the kekulized version of the tautomer (created lazily)
    // RDKit❗❌:         const auto &kmol = smilesTautomerPair.second.getKekulized();
    // RDKit❗❌:         std::vector<MatchVectType> matches;
    // RDKit❗❌:         unsigned int matched = SubstructMatch(*kmol, *(transform.Mol), matches);
    // RDKit❗❌:
    // RDKit❗❌:         if (!matched) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         ++nTransforms;
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:         std::string name;
    // RDKit❗❌:         (transform.Mol)->getProp(common_properties::_Name, name);
    // RDKit❗❌:         SmilesWriteParams smilesWriteParams;
    // RDKit❗❌:         smilesWriteParams.allBondsExplicit = true;
    // RDKit❗❌:         std::cout << "kmol for " << smilesTautomerPair.first << " : "
    // RDKit❗❌:                   << MolToSmiles(*kmol, smilesWriteParams) << std::endl;
    // RDKit❗❌:         std::cout << "transform mol: " << MolToSmarts(*(transform.Mol))
    // RDKit❗❌:                   << std::endl;
    // RDKit❗❌:
    // RDKit❗❌:         std::cout << "Matched: " << name << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:         // loop over transform matches
    // RDKit❗❌:         for (const auto &match : matches) {
    // RDKit❗❌:           if (nTransforms >= d_maxTransforms) {
    // RDKit❗❌:             res.d_status = TautomerEnumeratorStatus::MaxTransformsReached;
    // RDKit❗❌:             bailOut = true;
    // RDKit❗❌:           } else if (res.d_tautomers.size() >= d_maxTautomers) {
    // RDKit❗❌:             res.d_status = TautomerEnumeratorStatus::MaxTautomersReached;
    // RDKit❗❌:             bailOut = true;
    // RDKit❗❌:           } else if (d_callback.get() && !(*d_callback)(mol, res)) {
    // RDKit❗❌:             res.d_status = TautomerEnumeratorStatus::Canceled;
    // RDKit❗❌:             bailOut = true;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (bailOut) {
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:           // Create a copy of in the input molecule so we can modify it
    // RDKit❗❌:           // Use kekule form so bonds are explicitly single/double instead of
    // RDKit❗❌:           // aromatic
    // RDKit❗❌:           RWMOL_SPTR product(new RWMol(*kmol, true));
    // RDKit❗❌:           // Remove a hydrogen from the first matched atom and add one to the
    // RDKit❗❌:           // last
    // RDKit❗❌:           int firstIdx = match.front().second;
    // RDKit❗❌:           int lastIdx = match.back().second;
    // RDKit❗❌:           Atom *first = product->getAtomWithIdx(firstIdx);
    // RDKit❗❌:           Atom *last = product->getAtomWithIdx(lastIdx);
    // RDKit❗❌:           markAtomModified(static_cast<unsigned int>(firstIdx));
    // RDKit❗❌:           markAtomModified(static_cast<unsigned int>(lastIdx));
    // RDKit❗❌:           first->setNumExplicitHs(
    // RDKit❗❌:               std::max(0, static_cast<int>(first->getTotalNumHs()) - 1));
    // RDKit❗❌:           last->setNumExplicitHs(last->getTotalNumHs() + 1);
    // RDKit❗❌:           // Remove any implicit hydrogens from the first and last atoms
    // RDKit❗❌:           // now we have set the count explicitly
    // RDKit❗❌:           first->setNoImplicit(true);
    // RDKit❗❌:           last->setNoImplicit(true);
    // RDKit❗❌:           // Adjust bond orders
    // RDKit❗❌:           unsigned int bi = 0;
    // RDKit❗❌:           for (size_t i = 0; i < transform.Mol->getNumBonds(); ++i) {
    // RDKit❗❌:             const auto tbond = transform.Mol->getBondWithIdx(i);
    // RDKit❗❌:             Bond *bond = product->getBondBetweenAtoms(
    // RDKit❗❌:                 match[tbond->getBeginAtomIdx()].second,
    // RDKit❗❌:                 match[tbond->getEndAtomIdx()].second);
    // RDKit❗❌:             ASSERT_INVARIANT(bond, "required bond not found");
    // RDKit❗❌:             // check if bonds is specified in tautomer.in file
    // RDKit❗❌:             if (!transform.BondTypes.empty()) {
    // RDKit❗❌:               bond->setBondType(transform.BondTypes[bi]);
    // RDKit❗❌:               ++bi;
    // RDKit❗❌:             } else {
    // RDKit❗❌:               Bond::BondType bondtype = bond->getBondType();
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:               std::cout << "Bond as double: " << bond->getBondTypeAsDouble()
    // RDKit❗❌:                         << std::endl;
    // RDKit❗❌:               std::cout << bondtype << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:               if (bondtype == Bond::SINGLE) {
    // RDKit❗❌:                 bond->setBondType(Bond::DOUBLE);
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:                 std::cout << "Set bond to double" << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:               }
    // RDKit❗❌:               if (bondtype == Bond::DOUBLE) {
    // RDKit❗❌:                 bond->setBondType(Bond::SINGLE);
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:                 std::cout << "Set bond to single" << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             markBondModified(bond->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:           // TODO adjust charges
    // RDKit❗❌:           if (!transform.Charges.empty()) {
    // RDKit❗❌:             unsigned int ci = 0;
    // RDKit❗❌:             for (const auto &pair : match) {
    // RDKit❗❌:               Atom *atom = product->getAtomWithIdx(pair.second);
    // RDKit❗❌:               atom->setFormalCharge(atom->getFormalCharge() +
    // RDKit❗❌:                                     transform.Charges[ci++]);
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:           {
    // RDKit❗❌:             SmilesWriteParams smilesWriteParams;
    // RDKit❗❌:             smilesWriteParams.allBondsExplicit = true;
    // RDKit❗❌:             std::cout << "pre-sanitize: "
    // RDKit❗❌:                       << MolToSmiles(*product, smilesWriteParams) << std::endl;
    // RDKit❗❌:           }
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:           try {
    // RDKit❗❌:             // We only change bond orders/H counts/charges; the molecular graph
    // RDKit❗❌:             // (and therefore ring topology) is unchanged.
    // RDKit❗❌:             // `sanitizeMol()` always calls `clearComputedProps()` which resets
    // RDKit❗❌:             // ring info and forces ring-finding for each generated tautomer.
    // RDKit❗❌:             // Avoid that by clearing computed props without touching rings,
    // RDKit❗❌:             // then running the specific sanitize steps we need.
    // RDKit❗❌:             product->clearComputedProps(false);
    // RDKit❗❌:             product->updatePropertyCache(false);
    // RDKit❗❌:             MolOps::Kekulize(*product);
    // RDKit❗❌:             MolOps::setAromaticity(*product);
    // RDKit❗❌:             MolOps::setConjugation(*product);
    // RDKit❗❌:             MolOps::setHybridization(*product);
    // RDKit❗❌:             MolOps::adjustHs(*product);
    // RDKit❗❌:           } catch (const KekulizeException &) {
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:           SmilesWriteParams smilesWriteParams;
    // RDKit❗❌:           smilesWriteParams.allBondsExplicit = true;
    // RDKit❗❌:           std::cout << "pre-setTautomerStereo: "
    // RDKit❗❌:                     << MolToSmiles(*product, smilesWriteParams) << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:           setTautomerStereoAndIsoHs(mol, *product, res);
    // RDKit❗❌:           tsmiles = MolToSmiles(*product, true);
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:           (transform.Mol)->getProp(common_properties::_Name, name);
    // RDKit❗❌:           std::cout << "Applied rule: " << name << " to "
    // RDKit❗❌:                     << smilesTautomerPair.first << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:           if (res.d_tautomers.find(tsmiles) != res.d_tautomers.end()) {
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:             std::cout << "Previous tautomer produced again: " << tsmiles
    // RDKit❗❌:                       << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // RDKit❗❌:           // in addition to the above transformations, sanitization may modify
    // RDKit❗❌:           // bonds, e.g. Cc1nc2ccccc2[nH]1
    // RDKit❗❌:           // Use parallel iteration to avoid O(n) getBondWithIdx lookups
    // RDKit❗❌:           {
    // RDKit❗❌:             auto productBondIt = product->bonds().begin();
    // RDKit❗❌:             for (const auto molBond : mol.bonds()) {
    // RDKit❗❌:               const auto productBond = *productBondIt;
    // RDKit❗❌:               ++productBondIt;
    // RDKit❗❌:               auto i = molBond->getIdx();
    // RDKit❗❌:               if (molBond->getBondType() != productBond->getBondType() &&
    // RDKit❗❌:                   !res.d_modifiedBonds.test(i)) {
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:                 std::cout << "Sanitization has modified bond " << i
    // RDKit❗❌:                           << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:                 markBondModified(static_cast<unsigned int>(i));
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           // Kekulized form will be created lazily when needed;
    // RDKit❗❌:           // canonical=true is used on demand for order-independent deduplication.
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:           auto it = res.d_tautomers.find(tsmiles);
    // RDKit❗❌:           if (it == res.d_tautomers.end()) {
    // RDKit❗❌:             std::cout << "New tautomer added as ";
    // RDKit❗❌:           } else {
    // RDKit❗❌:             std::cout << "New tautomer replaced for ";
    // RDKit❗❌:           }
    // RDKit❗❌:           std::cout << tsmiles << ", taut: " << MolToSmiles(*product)
    // RDKit❗❌:                     << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:           // BOOST_LOG(rdInfoLog)
    // RDKit❗❌:           //     << "Tautomer transform "
    // RDKit❗❌:           //     <<
    // RDKit❗❌:           //     transform.Mol->getProp<std::string>(common_properties::_Name)
    // RDKit❗❌:           //     << " produced tautomer " << tsmiles << std::endl;
    // RDKit❗❌:           res.d_tautomers[tsmiles] = Tautomer(
    // RDKit❗❌:               std::move(product),
    // RDKit❗❌:               numModifiedAtoms, numModifiedBonds);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       smilesTautomerPair.second.d_done = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     completed = true;
    // RDKit❗❌:     size_t maxNumModifiedAtoms = numModifiedAtoms;
    // RDKit❗❌:     size_t maxNumModifiedBonds = numModifiedBonds;
    // RDKit❗❌:     for (auto it = res.d_tautomers.begin(); it != res.d_tautomers.end();) {
    // RDKit❗❌:       auto &taut = it->second;
    // RDKit❗❌:       if (!taut.d_done) {
    // RDKit❗❌:         completed = false;
    // RDKit❗❌:       }
    // RDKit❗❌:       if ((taut.d_numModifiedAtoms < maxNumModifiedAtoms ||
    // RDKit❗❌:            taut.d_numModifiedBonds < maxNumModifiedBonds) &&
    // RDKit❗❌:           setTautomerStereoAndIsoHs(mol, *taut.tautomer, res)) {
    // RDKit❗❌:         Tautomer tautStored = std::move(taut);
    // RDKit❗❌:         it = res.d_tautomers.erase(it);
    // RDKit❗❌:         tautStored.d_numModifiedAtoms = maxNumModifiedAtoms;
    // RDKit❗❌:         tautStored.d_numModifiedBonds = maxNumModifiedBonds;
    // RDKit❗❌:         auto insertRes = res.d_tautomers.insert(std::make_pair(
    // RDKit❗❌:             MolToSmiles(*tautStored.tautomer), std::move(tautStored)));
    // RDKit❗❌:         if (insertRes.second) {
    // RDKit❗❌:           it = insertRes.first;
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         ++it;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (bailOut && res.d_tautomers.size() < d_maxTautomers &&
    // RDKit❗❌:         res.d_status == TautomerEnumeratorStatus::MaxTautomersReached) {
    // RDKit❗❌:       res.d_status = TautomerEnumeratorStatus::Completed;
    // RDKit❗❌:       bailOut = false;
    // RDKit❗❌:     }
    // RDKit❗❌:   }  // while
    // RDKit❗❌:   res.fillTautomersItVec();
    // RDKit❗❌:   if (!completed) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Tautomer enumeration stopped at " << res.d_tautomers.size()
    // RDKit❗❌:         << " tautomers: " << statusMsg.at(static_cast<size_t>(res.d_status))
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::TautomerEnumerator::enumerate actual source cache loan
    // The source callback receives the original mol. The exclusive input loan
    // therefore belongs to the caller; enumeration never substitutes a cloned
    // source RingInfo. Candidate preparation still copies its own record.
    // RDKit✔️❌:   std::string smi = MolToSmiles(mol, true);
    let coordinates = source.coordinates;
    let key = canonical_smiles(source.as_record_view())?;
    let (initial, mut state) = initialize_candidates_with_key(source.as_record_view(), key)?;
    // Preserve original tags/CIP independently of candidate computed properties.
    let original = &initial;
    // RDKit✔️❌:   while (!completed && !bailOut) {
    loop {
        let pass = expand_tautomer_candidates_in_source_order(
            &mut state,
            catalog.transforms(),
            params,
            get_cached_kekulized,
            |candidate, transform| transform_matches(candidate, coordinates, transform),
            |candidate, transform, matched, atoms, bonds, existing| {
                apply_tautomer_transform_match(
                    original,
                    candidate,
                    coordinates,
                    transform,
                    matched,
                    atoms,
                    bonds,
                    existing,
                    params,
                )
            },
            |progress| match callback.as_deref_mut() {
                Some(callback) => callback
                    .should_continue(source.reborrow(), TautomerProgress { state: progress }),
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
                    coordinates,
                )
            },
            |candidate| canonical_smiles(candidate.view(coordinates)),
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
    candidates: impl ExactSizeIterator<
        Item = (
            &'a cosmolkit_model::PropertyText,
            crate::TautomerScoreView<'a>,
        ),
    >,
    mut scorer: impl for<'v> FnMut(crate::TautomerScoreView<'v>) -> Result<i32, TautomerRunError>,
) -> Result<usize, TautomerRunError> {
    select_canonical_index_by(
        candidates,
        |candidate| scorer(candidate.reborrow()),
        |key, _| Ok(std::borrow::Cow::Borrowed(key)),
    )
    .and_then(|index| index.ok_or(TautomerRunError::NoCanonicalTautomer))
}
/// Source iterable selection computes canonical keys only for a winning score or tie.
pub fn select_canonical_index_from_iterable_with<'a>(
    candidates: impl ExactSizeIterator<Item = crate::TautomerScoreView<'a>>,
    mut scorer: impl for<'v> FnMut(crate::TautomerScoreView<'v>) -> Result<i32, TautomerRunError>,
) -> Result<usize, TautomerRunError> {
    select_canonical_index_by(
        candidates.map(|value| ((), value)),
        |candidate| scorer(candidate.reborrow()),
        |(), view| canonical_smiles(view.as_record_view()).map(std::borrow::Cow::Owned),
    )
    .and_then(|index| index.ok_or(TautomerRunError::NoCanonicalTautomer))
}
// Retained results and source iterables share exactly one signed-score selection
// loop. Keys are borrowed in the result path and computed lazily in the iterable
// path, preserving scoring order, duplicate inputs and error timing.
#[doc(hidden)]
pub fn select_canonical_index_by<'a, M, C, E>(
    candidates: impl ExactSizeIterator<Item = (M, C)>,
    mut scorer: impl FnMut(&mut C) -> Result<i32, E>,
    mut key: impl FnMut(M, &C) -> Result<std::borrow::Cow<'a, cosmolkit_model::PropertyText>, E>,
) -> Result<Option<usize>, E> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerEnumerator::pickCanonical shared selection
    // RDKit❗❌: ROMol *TautomerEnumerator::pickCanonical(
    // RDKit❗❌:     const TautomerEnumeratorResult &tautRes,
    // RDKit❗❌:     boost::function<int(const ROMol &mol)> scoreFunc) const {
    // RDKit❗❌:   ROMOL_SPTR bestMol;
    // RDKit❗❌:   if (tautRes.d_tautomers.size() == 1) {
    // RDKit❗❌:     bestMol = tautRes.d_tautomers.begin()->second.tautomer;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // Calculate score for each tautomer
    // RDKit❗❌:     int bestScore = std::numeric_limits<int>::min();
    // RDKit❗❌:     std::string bestSmiles = "";
    // RDKit❗❌:     for (const auto &t : tautRes.d_tautomers) {
    // RDKit❗❌:       auto score = scoreFunc(*t.second.tautomer);
    // RDKit❗❌: #ifdef VERBOSE_ENUMERATION
    // RDKit❗❌:       std::cerr << "  " << t.first << " " << score << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:       if (score > bestScore) {
    // RDKit❗❌:         bestScore = score;
    // RDKit❗❌:         bestSmiles = t.first;
    // RDKit❗❌:         bestMol = t.second.tautomer;
    // RDKit❗❌:       } else if (score == bestScore) {
    // RDKit❗❌:         if (t.first < bestSmiles) {
    // RDKit❗❌:           bestSmiles = t.first;
    // RDKit❗❌:           bestMol = t.second.tautomer;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   ROMol *res = new ROMol(*bestMol);
    // RDKit❗❌:   static const bool cleanIt = true;
    // RDKit❗❌:   static const bool force = true;
    // RDKit❗❌:   MolOps::assignStereochemistry(*res, cleanIt, force);
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::TautomerEnumerator::pickCanonical shared selection
    // The source final copy/forced stereo assignment remains in the existing
    // finalize_canonical_candidate caller. This is the sole selection loop.
    if candidates.len() == 1 {
        return Ok(Some(0));
    }
    let mut best_score = i32::MIN;
    let mut best_smiles: std::borrow::Cow<'a, cosmolkit_model::PropertyText> =
        std::borrow::Cow::Owned(cosmolkit_model::PropertyText::new());
    let mut best = None;
    for (index, (metadata, mut candidate)) in candidates.enumerate() {
        let score = scorer(&mut candidate)?;
        if score >= best_score {
            let smiles = key(metadata, &candidate)?;
            if score > best_score || smiles < best_smiles {
                best_score = score;
                best_smiles = smiles;
                best = Some(index);
            }
        }
    }
    Ok(best)
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
    result: &mut TautomerEnumerationOutput,
    coordinates: &CoordinateBlock,
    scorer: impl FnMut(crate::TautomerScoreView<'_>) -> Result<i32, TautomerRunError>,
) -> Result<TautomerRecord, TautomerRunError> {
    let index = select_canonical_index_with(
        result
            .entries
            .iter_mut()
            .map(|(key, value)| (&*key, value.score_view(coordinates))),
        scorer,
    )?;
    finalize_canonical_candidate(result.entries[index].1.view(coordinates))
}
pub fn canonicalize_with_catalog(
    mut source: crate::TautomerScoreView<'_>,
    catalog: &TautomerCatalog,
    params: TautomerParams,
    callback: Option<&mut dyn TautomerEnumerationCallback>,
    scorer: Option<&mut dyn FnMut(crate::TautomerScoreView<'_>) -> Result<i32, TautomerRunError>>,
) -> Result<TautomerRecord, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MolStandardize::TautomerEnumerator::canonicalize
    // RDKit❗❌: ROMol *TautomerEnumerator::canonicalize(
    // RDKit❗❌:     const ROMol &mol, boost::function<int(const ROMol &mol)> scoreFunc) const {
    // RDKit❗❌:   auto thisCopy = TautomerEnumerator(*this);
    // RDKit❗❌:   thisCopy.setReassignStereo(false);
    // RDKit❗❌:   auto res = thisCopy.enumerate(mol);
    // RDKit❗❌:   if (res.empty()) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "no tautomers found, returning input molecule" << std::endl;
    // RDKit❗❌:     return new ROMol(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:   // When no custom scorer provided, use optimized scoring that pre-filters
    // RDKit❗❌:   // SubstructTerm patterns once for the input molecule rather than evaluating
    // RDKit❗❌:   // all 12 substructure matches per tautomer. This is safe because
    // RDKit❗❌:   // tautomerization only moves H and changes bond orders, never creates or
    // RDKit❗❌:   // destroys heavy-atom bonds.
    // RDKit❗❌:   if (!scoreFunc) {
    // RDKit❗❌:     scoreFunc = TautomerScoringFunctions::makeOptimizedScorer(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:   ROMol *canonical = pickCanonical(res, scoreFunc);
    // RDKit❗❌:   // quickCopy during enumeration doesn't copy molecule properties or
    // RDKit❗❌:   // conformers.  Restore both from the original molecule so that
    // RDKit❗❌:   // downstream code (e.g. InChI generation) that relies on 2D/3D
    // RDKit❗❌:   // coordinates or mol-level properties works correctly.
    // RDKit❗❌:   canonical->updateProps(mol);
    // RDKit❗❌:   for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
    // RDKit❗❌:        ++confIt) {
    // RDKit❗❌:     canonical->addConformer(new Conformer(**confIt), true);
    // RDKit❗❌:   }
    // RDKit❗❌:   return canonical;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::MolStandardize::TautomerEnumerator::canonicalize

    let mut result = enumerate_with_catalog(
        source.reborrow(),
        catalog,
        params.with_reassign_stereo(false),
        callback,
    )?;
    if result.entries.is_empty() {
        return prepared(source.as_record_view());
    }
    let mut canonical = match scorer {
        Some(custom) => pick_canonical_with(&mut result, source.coordinates, custom),
        None => {
            let optimized = crate::score::make_optimized_tautomer_scorer(source.as_record_view())?;
            pick_canonical_default_with_ring_cache(&mut result, source.coordinates, &optimized)
        }
    }?;
    restore_canonical_source_metadata(&mut canonical, source.as_record_view())?;
    Ok(canonical)
}

#[cfg(test)]
mod initialization_tests;

fn restore_canonical_source_metadata(
    canonical: &mut TautomerRecord,
    source: TautomerRecordView<'_>,
) -> Result<(), TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerEnumerator::canonicalize metadata restoration
    // RDKit❗❌:   // quickCopy during enumeration doesn't copy molecule properties or
    // RDKit❗❌:   // conformers.  Restore both from the original molecule so that
    // RDKit❗❌:   // downstream code (e.g. InChI generation) that relies on 2D/3D
    // RDKit❗❌:   // coordinates or mol-level properties works correctly.
    // RDKit❗❌:   canonical->updateProps(mol);
    // RDKit❗❌:   for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
    // RDKit❗❌:        ++confIt) {
    // RDKit❗❌:     canonical->addConformer(new Conformer(**confIt), true);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION RDKit::TautomerEnumerator::canonicalize metadata restoration
    // BEGIN RDKIT CPP FUNCTION RDKit::RDProps::updateProps
    // RDKit❗❌:   void updateProps(const RDProps &source, bool preserveExisting = false) {
    // RDKit❗❌:     d_props.update(source.getDict(), preserveExisting);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION RDKit::RDProps::updateProps
    // Default preserveExisting=false REPLACES the dictionary, including the
    // computed-property vector. The modeled name/SDF projections travel with it.
    canonical.properties = source.properties.clone();
    let mut coordinates = canonical
        .coordinates
        .take()
        .unwrap_or_else(|| source.coordinates.clone());
    let order = match &source.coordinates.source_conformer_order {
        Some(order) => order.clone(),
        None if source.coordinates.conformers_2d.is_empty() => {
            vec![
                cosmolkit_model::CoordinateDimension::ThreeD;
                source.coordinates.conformers_3d.len()
            ]
        }
        None if source.coordinates.conformers_3d.is_empty() => {
            vec![cosmolkit_model::CoordinateDimension::TwoD; source.coordinates.conformers_2d.len()]
        }
        None => {
            return Err(
                cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder.into(),
            );
        }
    };
    let (mut two, mut three) = (0, 0);
    for dimension in order {
        let id = source_next_canonical_conformer_id(&coordinates)?;
        coordinates.record_source_conformer_append(dimension)?;
        match dimension {
            cosmolkit_model::CoordinateDimension::TwoD => {
                let frame = source.coordinates.conformers_2d.get(two).ok_or(
                    cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
                )?;
                coordinates.conformers_2d.push(frame.clone().with_id(id));
                two += 1;
            }
            cosmolkit_model::CoordinateDimension::ThreeD => {
                let frame = source.coordinates.conformers_3d.get(three).ok_or(
                    cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder,
                )?;
                coordinates.conformers_3d.push(frame.clone().with_id(id));
                three += 1;
            }
        }
    }
    if two != source.coordinates.conformers_2d.len()
        || three != source.coordinates.conformers_3d.len()
    {
        return Err(cosmolkit_model::CoordinateValidationError::MissingSourceConformerOrder.into());
    }
    canonical.coordinates = Some(coordinates);
    Ok(())
}
fn source_next_canonical_conformer_id(
    coordinates: &CoordinateBlock,
) -> Result<usize, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::ROMol::addConformer
    // RDKit❗❌: unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
    // RDKit❗❌:   PRECONDITION(conf, "bad conformer");
    // RDKit❗❌:   PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
    // RDKit❗❌:                "Number of atom mismatch");
    // RDKit❗❌:   if (assignId) {
    // RDKit❗❌:     int maxId = -1;
    // RDKit❗❌:     for (auto cptr : d_confs) {
    // RDKit❗❌:       maxId = std::max((int)(cptr->getId()), maxId);
    // RDKit❗❌:     }
    // RDKit❗❌:     maxId++;
    // RDKit❗❌:     conf->setId((unsigned int)maxId);
    // RDKit❗❌:   }
    // RDKit❗❌:   conf->setOwningMol(this);
    // RDKit❗❌:   CONFORMER_SPTR nConf(conf);
    // RDKit❗❌:   d_confs.push_back(nConf);
    // RDKit❗❌:   return conf->getId();
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::ROMol::addConformer
    // Source scans ALL conformers. Existing dimension-local add helpers cannot
    // implement this append. Signed overflow is undefined in C++; report it.
    let max = coordinates
        .conformers_2d
        .iter()
        .map(|c| c.id())
        .chain(coordinates.conformers_3d.iter().map(|c| c.id()))
        .map(|id| id as u32 as i32)
        .max()
        .unwrap_or(-1)
        .max(-1);
    let id = max.checked_add(1).ok_or(
        cosmolkit_model::CoordinateValidationError::ConformerIdOverflow {
            max_id: max as usize,
        },
    )?;
    Ok(id as u32 as usize)
}

#[cfg(test)]
mod search04_metadata_tests {
    use super::*;
    use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateDimension};
    fn coordinates(atoms: usize) -> CoordinateBlock {
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(9, vec![[1.0, -2.0]; atoms])],
            conformers_3d: vec![Conformer3D::new(2, vec![[3.0, 4.0, 5.0]; atoms], true)],
            source_coordinate_dim: Some(CoordinateDimension::TwoD),
            source_conformer_order: Some(vec![
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD,
            ]),
        }
    }
    #[test]
    fn search04_canonical_restores_quickcopy_properties_and_ordered_conformers() {
        let mut source = crate::engine::stereo_tests::fixture_from_smiles("CC=O").unwrap();
        source.properties.set_prop("input", 7_i32).unwrap();
        source
            .properties
            .set_computed_prop("computed-input", 9_i32)
            .unwrap();
        let coords = coordinates(source.topology.atoms.len());
        let view = source.view(&coords);
        let mut product = source.clone();
        product.coordinates = Some(CoordinateBlock::default());
        product.properties = Default::default();
        product.properties.set_prop("product-only", 3_i32).unwrap();
        restore_canonical_source_metadata(&mut product, view).unwrap();
        assert_eq!(product.properties, source.properties);
        let restored = product.coordinates.unwrap();
        assert_eq!(restored.conformers_2d[0].id(), 0);
        assert_eq!(restored.conformers_3d[0].id(), 1);
        assert_eq!(
            restored.conformers_2d[0].coordinates(),
            coords.conformers_2d[0].coordinates()
        );
        assert_eq!(
            restored.conformers_3d[0].coordinates(),
            coords.conformers_3d[0].coordinates()
        );
        assert_eq!(
            restored.source_conformer_order,
            coords.source_conformer_order
        );
        assert_eq!(coords.conformers_2d[0].id(), 9);
    }
    #[test]
    fn search04_canonical_initial_selection_appends_without_clearing_existing_conformers() {
        let source = crate::engine::stereo_tests::fixture_from_smiles("CC").unwrap();
        let coords = coordinates(source.topology.atoms.len());
        let mut selected =
            crate::engine::copy_for_canonical_assignment(source.view(&coords)).unwrap();
        restore_canonical_source_metadata(&mut selected, source.view(&coords)).unwrap();
        let restored = selected.coordinates.unwrap();
        assert_eq!(
            restored
                .conformers_2d
                .iter()
                .map(|c| c.id())
                .collect::<Vec<_>>(),
            [9, 10]
        );
        assert_eq!(
            restored
                .conformers_3d
                .iter()
                .map(|c| c.id())
                .collect::<Vec<_>>(),
            [2, 11]
        );
        assert_eq!(
            restored.source_conformer_order,
            Some(vec![
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD
            ])
        );
    }
}

fn pick_canonical_default_with_ring_cache(
    result: &mut TautomerEnumerationOutput,
    coordinates: &CoordinateBlock,
    scorer: &crate::score::OptimizedTautomerScorer,
) -> Result<TautomerRecord, TautomerRunError> {
    let index = select_canonical_index_by(
        result
            .entries
            .iter_mut()
            .map(|(key, record)| (&*key, record)),
        |record| scorer.score_record(record, coordinates),
        |key, _| Ok(std::borrow::Cow::Borrowed(key)),
    )?
    .ok_or(TautomerRunError::NoCanonicalTautomer)?;
    finalize_canonical_candidate(result.entries[index].1.view(coordinates))
}

#[cfg(test)]
mod search06_default_picker_tests {
    use super::*;
    fn cold(text: &str) -> TautomerRecord {
        let mut record = crate::engine::stereo_tests::fixture_from_smiles(text).unwrap();
        record.rings = cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            record.topology.atoms.len(),
            record.topology.bonds.len(),
        );
        record
    }
    #[test]
    fn search06_default_picker_mutates_each_actual_candidate_cache() {
        let coordinates = CoordinateBlock::default();
        let mut result = TautomerEnumerationOutput {
            entries: vec![
                ("CC".into(), cold("CC")),
                ("c1ccccc1".into(), cold("c1ccccc1")),
            ],
            ..Default::default()
        };
        let scorer =
            crate::score::make_optimized_tautomer_scorer(result.entries[1].1.view(&coordinates))
                .unwrap();
        let mut selected =
            pick_canonical_default_with_ring_cache(&mut result, &coordinates, &scorer).unwrap();
        assert!(result.entries.iter().all(|(_, r)| r.rings.is_symm_sssr()));
        assert_eq!(
            crate::score::score_tautomer_rings_(&mut selected.score_view(&coordinates)).unwrap(),
            250
        );
    }
    #[test]
    fn search06_single_default_candidate_skips_score_and_cache_write() {
        let coordinates = CoordinateBlock::default();
        let mut result = TautomerEnumerationOutput {
            entries: vec![("CC".into(), cold("CC"))],
            ..Default::default()
        };
        let scorer =
            crate::score::make_optimized_tautomer_scorer(result.entries[0].1.view(&coordinates))
                .unwrap();
        pick_canonical_default_with_ring_cache(&mut result, &coordinates, &scorer).unwrap();
        assert!(!result.entries[0].1.rings.is_symm_sssr());
    }
}

#[cfg(test)]
mod recovery_search06_domain_source_loan_tests {
    use super::*;
    use cosmolkit_core::{RingFindType, RingInfo};
    struct ScoreAndStop {
        error: bool,
    }
    impl TautomerEnumerationCallback for ScoreAndStop {
        fn should_continue(
            &mut self,
            mut source: crate::TautomerScoreView<'_>,
            mut progress: TautomerProgress<'_>,
        ) -> Result<bool, TautomerRunError> {
            crate::score_tautomer_rings_(&mut source)?;
            for (_, mut candidate) in progress.entries_mut(source.coordinates) {
                crate::score_tautomer_rings_(&mut candidate)?;
            }
            assert!(
                progress
                    .entries_readonly()
                    .all(|(_, record)| record.rings.is_symm_sssr())
            );
            if self.error {
                Err(TautomerRunError::Callback(
                    "after source-cache write".into(),
                ))
            } else {
                Ok(false)
            }
        }
    }
    struct NoScore;
    impl TautomerEnumerationCallback for NoScore {
        fn should_continue(
            &mut self,
            _: crate::TautomerScoreView<'_>,
            _: TautomerProgress<'_>,
        ) -> Result<bool, TautomerRunError> {
            Ok(false)
        }
    }
    fn cold_record() -> TautomerRecord {
        let mut value = crate::engine::stereo_tests::fixture_from_smiles("CC(C)=O").unwrap();
        value.rings = RingInfo::new(
            RingFindType::OtherOrUnknown,
            value.topology.atoms.len(),
            value.topology.bonds.len(),
        );
        value
    }
    #[test]
    fn recovery_search06_domain_enumeration_callback_error_keeps_actual_source_cache() {
        let mut source = cold_record();
        let before = source.clone();
        let coordinates = CoordinateBlock::default();
        let catalog = TautomerCatalog::current().unwrap();
        let mut callback = ScoreAndStop { error: true };
        let result = enumerate_with_catalog(
            source.score_view(&coordinates),
            &catalog,
            Default::default(),
            Some(&mut callback),
        );
        assert!(matches!(result, Err(TautomerRunError::Callback(_))));
        assert!(source.rings.is_symm_sssr());
        assert_eq!(source.topology, before.topology);
        assert_eq!(source.valence, before.valence);
        assert_eq!(source.properties, before.properties);
    }
    #[test]
    fn recovery_search06_domain_canonical_callback_error_keeps_actual_source_cache() {
        let mut source = cold_record();
        let coordinates = CoordinateBlock::default();
        let catalog = TautomerCatalog::current().unwrap();
        let mut callback = ScoreAndStop { error: true };
        let result = canonicalize_with_catalog(
            source.score_view(&coordinates),
            &catalog,
            Default::default(),
            Some(&mut callback),
            None,
        );
        assert!(matches!(result, Err(TautomerRunError::Callback(_))));
        assert!(source.rings.is_symm_sssr());
    }
    #[test]
    fn recovery_search06_domain_no_score_callback_does_not_eagerly_warm_source() {
        let mut source = cold_record();
        let coordinates = CoordinateBlock::default();
        let catalog = TautomerCatalog::current().unwrap();
        let mut callback = NoScore;
        let result = enumerate_with_catalog(
            source.score_view(&coordinates),
            &catalog,
            Default::default(),
            Some(&mut callback),
        )
        .unwrap();
        assert_eq!(result.status, TautomerEnumerationStatus::Canceled);
        assert!(!source.rings.is_symm_sssr());
    }
    #[test]
    fn recovery_search06_domain_empty_catalog_keeps_actual_source_cold() {
        let mut source = cold_record();
        let coordinates = CoordinateBlock::default();
        let catalog = TautomerCatalog::from_data(&[]).unwrap();
        let result = enumerate_with_catalog(
            source.score_view(&coordinates),
            &catalog,
            TautomerParams::default().with_reassign_stereo(false),
            None,
        )
        .unwrap();
        assert_eq!(result.entries.len(), 1);
        assert!(!source.rings.is_symm_sssr());
    }
}
