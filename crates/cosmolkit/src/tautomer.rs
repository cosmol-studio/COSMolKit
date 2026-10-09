//! Public tautomer configuration and finalized values over the unique detached owner.
use crate::{Atom, AtomId, Bond, BondId, Molecule, OperationError};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties};
pub use cosmolkit_tautomer::{
    TautomerCatalogError, TautomerEnumerationStatus, TautomerRunError, TautomerScore,
    TautomerScoreTerm,
};
use std::{
    collections::BTreeSet,
    fmt,
    sync::{Arc, OnceLock},
};

/// Source score terms, in their source order.
pub fn default_tautomer_score_terms() -> &'static [TautomerScoreTerm] {
    cosmolkit_tautomer::default_tautomer_score_terms()
}

/// Terms used by the score query or restricted callback score. None selects the source defaults.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TautomerScoreParams {
    pub terms: Option<Vec<TautomerScoreTerm>>,
}

/// Readonly graph inspection with an exclusive score-cache loan during callbacks.
pub struct TautomerMoleculeView<'a> {
    storage: TautomerViewStorage<'a>,
}
enum TautomerViewStorage<'a> {
    Borrowed(cosmolkit_tautomer::TautomerScoreView<'a>),
    Owned(Arc<TautomerViewSnapshot>),
}
#[derive(Clone)]
struct TautomerViewSnapshot {
    topology: cosmolkit_model::TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    valence: Option<cosmolkit_core::ValenceAssignment>,
    rings: Option<cosmolkit_core::RingInfo>,
}
impl<'a> TautomerMoleculeView<'a> {
    fn borrowed(view: cosmolkit_tautomer::TautomerScoreView<'a>) -> Self {
        Self {
            storage: TautomerViewStorage::Borrowed(view),
        }
    }
    fn record(&self) -> cosmolkit_tautomer::TautomerRecordView<'_> {
        match &self.storage {
            TautomerViewStorage::Borrowed(view) => view.as_record_view(),
            TautomerViewStorage::Owned(value) => cosmolkit_tautomer::TautomerRecordView {
                topology: &value.topology,
                coordinates: &value.coordinates,
                properties: &value.properties,
                valence: value.valence.as_ref(),
                rings: value.rings.as_ref(),
            },
        }
    }
    /// Retain an independent owned inspection/score snapshot after the loan.
    pub fn to_owned(&self) -> TautomerMoleculeView<'static> {
        match &self.storage {
            TautomerViewStorage::Owned(value) => TautomerMoleculeView {
                storage: TautomerViewStorage::Owned(value.clone()),
            },
            TautomerViewStorage::Borrowed(view) => Self::snapshot(view.as_record_view()),
        }
    }
    pub fn num_atoms(&self) -> usize {
        self.record().topology.atoms.len()
    }
    pub fn num_bonds(&self) -> usize {
        self.record().topology.bonds.len()
    }
    /// Canonical metadata query through its unique foundational owner.
    pub fn atom_metadata(&self) -> Result<Vec<crate::AtomMetadata>, crate::ValenceError> {
        let record = self.record();
        cosmolkit_core::atom_metadata_from_assignment(record.topology, record.valence)
    }
    /// Degree from the validated detached adjacency, independent of valence errors.
    pub fn atom_degree(&self, id: AtomId) -> Option<usize> {
        let record = self.record();
        record.topology.atoms.get(id.index())?;
        Some(record.topology.adjacency.neighbors_of(id.index()).len())
    }
    pub fn atoms(&self) -> &[Atom] {
        &self.record().topology.atoms
    }
    pub fn bonds(&self) -> &[Bond] {
        &self.record().topology.bonds
    }
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.record().topology.atoms.get(id.index())
    }
    pub fn bond(&self, id: BondId) -> Option<&Bond> {
        self.record().topology.bonds.get(id.index())
    }
    pub fn properties(&self) -> &MoleculeProperties {
        self.record().properties
    }
    pub fn to_smiles(&self) -> Result<cosmolkit_model::PropertyText, TautomerRunError> {
        let view = self.record();
        Ok(cosmolkit_smiles::write_smiles(
            cosmolkit_smiles::SmilesRecordView {
                topology: view.topology,
                coordinates: view.coordinates,
                properties: view.properties,
                rings: view.rings,
            },
        )?)
    }
}
/// Source-ordered progress before application of a matched transform.
pub struct TautomerProgress<'a> {
    storage: TautomerProgressStorage<'a>,
}
enum TautomerProgressStorage<'a> {
    Borrowed {
        inner: cosmolkit_tautomer::TautomerProgress<'a>,
        coordinates: &'a CoordinateBlock,
    },
    Owned {
        entries: Vec<(cosmolkit_model::PropertyText, TautomerMoleculeView<'static>)>,
        status: TautomerEnumerationStatus,
        num_transforms: u32,
        modified_atoms: BTreeSet<AtomId>,
        modified_bonds: BTreeSet<BondId>,
    },
}
impl TautomerProgress<'_> {
    /// Preserve this exact pre-application snapshot for a language callback.
    pub fn to_owned(&self) -> TautomerProgress<'static> {
        let entries = match &self.storage {
            TautomerProgressStorage::Borrowed { inner, coordinates } => inner
                .entries_readonly()
                .map(|(key, record)| {
                    (
                        key.to_owned(),
                        TautomerMoleculeView::snapshot(record.view(coordinates)),
                    )
                })
                .collect(),
            TautomerProgressStorage::Owned { entries, .. } => entries
                .iter()
                .map(|(key, view)| (key.to_owned(), view.to_owned()))
                .collect(),
        };
        TautomerProgress {
            storage: TautomerProgressStorage::Owned {
                entries,
                status: self.status(),
                num_transforms: self.num_transforms(),
                modified_atoms: self.modified_atoms().clone(),
                modified_bonds: self.modified_bonds().clone(),
            },
        }
    }
    pub fn len(&self) -> usize {
        match &self.storage {
            TautomerProgressStorage::Borrowed { inner, .. } => inner.len(),
            TautomerProgressStorage::Owned { entries, .. } => entries.len(),
        }
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    pub fn status(&self) -> TautomerEnumerationStatus {
        match &self.storage {
            TautomerProgressStorage::Borrowed { inner, .. } => inner.status(),
            TautomerProgressStorage::Owned { status, .. } => *status,
        }
    }
    pub fn num_transforms(&self) -> u32 {
        match &self.storage {
            TautomerProgressStorage::Borrowed { inner, .. } => inner.num_transforms(),
            TautomerProgressStorage::Owned { num_transforms, .. } => *num_transforms,
        }
    }
    pub fn modified_atoms(&self) -> &BTreeSet<AtomId> {
        match &self.storage {
            TautomerProgressStorage::Borrowed { inner, .. } => inner.modified_atoms(),
            TautomerProgressStorage::Owned { modified_atoms, .. } => modified_atoms,
        }
    }
    pub fn modified_bonds(&self) -> &BTreeSet<BondId> {
        match &self.storage {
            TautomerProgressStorage::Borrowed { inner, .. } => inner.modified_bonds(),
            TautomerProgressStorage::Owned { modified_bonds, .. } => modified_bonds,
        }
    }
    pub fn entries(
        &mut self,
    ) -> Box<
        dyn ExactSizeIterator<Item = (&cosmolkit_model::PropertyText, TautomerMoleculeView<'_>)>
            + '_,
    > {
        match &mut self.storage {
            TautomerProgressStorage::Borrowed { inner, coordinates } => Box::new(
                inner
                    .entries_mut(coordinates)
                    .map(|(key, view)| (key, TautomerMoleculeView::borrowed(view))),
            ),
            TautomerProgressStorage::Owned { entries, .. } => Box::new(
                entries
                    .iter_mut()
                    .map(|(key, view)| (&*key, view.reborrow())),
            ),
        }
    }
}
/// Returning false cancels the enumeration; an error aborts the value operation.
pub trait TautomerEnumerationCallback: Send + Sync {
    fn should_continue(
        &self,
        source: TautomerMoleculeView<'_>,
        progress: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError>;
}
/// A signed custom score. Lexical ties use the source's retained canonical keys.
pub trait TautomerScorer: Send + Sync {
    fn score(&self, molecule: TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError>;
}

/// Source policy and shared immutable catalog/callback/scorer configuration.
#[derive(Clone)]
pub struct TautomerParams {
    pub(crate) policy: cosmolkit_tautomer::TautomerParams,
    pub(crate) catalog: Arc<cosmolkit_tautomer::TautomerCatalog>,
    // Internal invocation mode for finalizing a selected validated result.
    // It is never exposed or returned as public parameter configuration.
    pub(crate) finalize_selected: bool,
    callback: Option<Arc<dyn TautomerEnumerationCallback>>,
    scorer: Option<Arc<dyn TautomerScorer>>,
    pub score_params: TautomerScoreParams,
}
impl fmt::Debug for TautomerParams {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TautomerParams")
            .field("policy", &self.policy)
            .field("transform_count", &self.catalog.transforms().len())
            .field("callback", &self.callback.is_some())
            .field("scorer", &self.scorer.is_some())
            .field("score_params", &self.score_params)
            .finish()
    }
}
impl Default for TautomerParams {
    fn default() -> Self {
        // RDKit✔️🔝:   dp_catalog.reset(new TautomerCatalog(tautParams.get()));
        // Compiled builtins are immutable: sharing the verified source rows
        // avoids recompilation while each parameter value owns its own flags.
        static CATALOG: OnceLock<Arc<cosmolkit_tautomer::TautomerCatalog>> = OnceLock::new();
        Self {
            policy: Default::default(),
            finalize_selected: false,
            catalog: CATALOG
                .get_or_init(|| {
                    Arc::new(
                        cosmolkit_tautomer::TautomerCatalog::current()
                            .expect("verified source builtin catalog"),
                    )
                })
                .clone(),
            callback: None,
            scorer: None,
            score_params: Default::default(),
        }
    }
}
impl TautomerParams {
    pub fn v1() -> Result<Self, TautomerCatalogError> {
        Ok(Self {
            catalog: Arc::new(cosmolkit_tautomer::TautomerCatalog::v1()?),
            ..Default::default()
        })
    }
    pub fn from_transform_data(
        data: &[(&str, &str, &str, &str)],
    ) -> Result<Self, TautomerCatalogError> {
        Ok(Self {
            catalog: Arc::new(cosmolkit_tautomer::TautomerCatalog::from_data(data)?),
            ..Default::default()
        })
    }
    pub fn from_transform_file(
        path: impl AsRef<std::path::Path>,
    ) -> Result<Self, TautomerCatalogError> {
        Ok(Self {
            catalog: Arc::new(cosmolkit_tautomer::TautomerCatalog::from_file(path)?),
            ..Default::default()
        })
    }
    pub fn transform_count(&self) -> usize {
        self.catalog.transforms().len()
    }
    pub fn callback(&self) -> Option<&Arc<dyn TautomerEnumerationCallback>> {
        self.callback.as_ref()
    }
    pub fn set_callback(&mut self, callback: Option<Arc<dyn TautomerEnumerationCallback>>) {
        self.callback = callback;
    }
    pub fn scorer(&self) -> Option<&Arc<dyn TautomerScorer>> {
        self.scorer.as_ref()
    }
    pub fn set_scorer(&mut self, scorer: Option<Arc<dyn TautomerScorer>>) {
        self.scorer = scorer;
    }
    pub(crate) fn score_view(
        &self,
        mut view: cosmolkit_tautomer::TautomerScoreView<'_>,
    ) -> Result<i32, TautomerRunError> {
        match &self.scorer {
            Some(scorer) => scorer.score(TautomerMoleculeView::borrowed(view)),
            None => score_view(&mut view, &self.score_params).map(TautomerScore::total),
        }
    }
}
impl TautomerParams {
    pub fn max_tautomers(&self) -> u32 {
        self.policy.max_tautomers()
    }
    pub fn set_max_tautomers(&mut self, value: u32) {
        self.policy.set_max_tautomers(value);
    }
    pub fn with_max_tautomers(mut self, value: u32) -> Self {
        self.set_max_tautomers(value);
        self
    }
}
impl TautomerParams {
    pub fn max_transforms(&self) -> u32 {
        self.policy.max_transforms()
    }
    pub fn set_max_transforms(&mut self, value: u32) {
        self.policy.set_max_transforms(value);
    }
    pub fn with_max_transforms(mut self, value: u32) -> Self {
        self.set_max_transforms(value);
        self
    }
}
impl TautomerParams {
    pub fn remove_sp3_stereo(&self) -> bool {
        self.policy.remove_sp3_stereo()
    }
    pub fn set_remove_sp3_stereo(&mut self, value: bool) {
        self.policy.set_remove_sp3_stereo(value);
    }
    pub fn with_remove_sp3_stereo(mut self, value: bool) -> Self {
        self.set_remove_sp3_stereo(value);
        self
    }
}
impl TautomerParams {
    pub fn remove_bond_stereo(&self) -> bool {
        self.policy.remove_bond_stereo()
    }
    pub fn set_remove_bond_stereo(&mut self, value: bool) {
        self.policy.set_remove_bond_stereo(value);
    }
    pub fn with_remove_bond_stereo(mut self, value: bool) -> Self {
        self.set_remove_bond_stereo(value);
        self
    }
}
impl TautomerParams {
    pub fn remove_isotopic_hydrogens(&self) -> bool {
        self.policy.remove_isotopic_hydrogens()
    }
    pub fn set_remove_isotopic_hydrogens(&mut self, value: bool) {
        self.policy.set_remove_isotopic_hydrogens(value);
    }
    pub fn with_remove_isotopic_hydrogens(mut self, value: bool) -> Self {
        self.set_remove_isotopic_hydrogens(value);
        self
    }
}
impl TautomerParams {
    pub fn reassign_stereo(&self) -> bool {
        self.policy.reassign_stereo()
    }
    pub fn set_reassign_stereo(&mut self, value: bool) {
        self.policy.set_reassign_stereo(value);
    }
    pub fn with_reassign_stereo(mut self, value: bool) -> Self {
        self.set_reassign_stereo(value);
        self
    }
}

pub(crate) struct CallbackAdapter<'a> {
    pub(crate) params: &'a TautomerParams,
}
impl cosmolkit_tautomer::TautomerEnumerationCallback for CallbackAdapter<'_> {
    fn should_continue(
        &mut self,
        mut source: cosmolkit_tautomer::TautomerScoreView<'_>,
        progress: cosmolkit_tautomer::TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError> {
        let coordinates = source.coordinates;
        let result = match &self.params.callback {
            Some(callback) => callback.should_continue(
                TautomerMoleculeView::borrowed(source.reborrow()),
                TautomerProgress {
                    storage: TautomerProgressStorage::Borrowed {
                        inner: progress,
                        coordinates,
                    },
                },
            ),
            None => Ok(true),
        };
        // The score view retains its local cache, including before a later
        // callback error. The public receiver is never a cache write target.
        result
    }
}
/// Rich result containing only runtime-validated molecule values.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TautomerEnumeration {
    entries: Vec<(cosmolkit_model::PropertyText, Molecule)>,
    status: TautomerEnumerationStatus,
    modified_atoms: BTreeSet<AtomId>,
    modified_bonds: BTreeSet<BondId>,
}
impl TautomerEnumeration {
    pub fn canonical_tautomer(&self) -> Result<Molecule, OperationError> {
        self.canonical_tautomer_with_params(&TautomerParams::default())
    }
    pub fn canonical_tautomer_with_params(
        &self,
        params: &TautomerParams,
    ) -> Result<Molecule, OperationError> {
        // COW candidates retain scorer cache transitions for selection without
        // modifying the enumeration or any molecule previously borrowed from it.
        let mut entries = self.entries.clone();
        let index = cosmolkit_tautomer::select_canonical_index_by(
            entries.iter_mut().map(|(key, molecule)| (&*key, molecule)),
            |molecule| score_live_candidate(molecule, params),
            |key, _| Ok(std::borrow::Cow::Borrowed(key)),
        )?
        .ok_or_else(|| OperationError::Tautomer(TautomerRunError::NoCanonicalTautomer))?;
        let mut invocation = params.clone();
        invocation.finalize_selected = true;
        entries[index].1.canonical_tautomer_with_params(&invocation)
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn status(&self) -> TautomerEnumerationStatus {
        self.status
    }
    pub fn modified_atoms(&self) -> &BTreeSet<AtomId> {
        &self.modified_atoms
    }
    pub fn modified_bonds(&self) -> &BTreeSet<BondId> {
        &self.modified_bonds
    }
    pub fn get(&self, index: usize) -> Option<&Molecule> {
        self.entries.get(index).map(|(_, m)| m)
    }
    pub fn canonical_smiles(&self) -> Vec<&cosmolkit_model::PropertyText> {
        self.entries.iter().map(|(k, _)| k).collect()
    }
    pub fn iter(
        &self,
    ) -> std::iter::Map<
        std::slice::Iter<'_, (cosmolkit_model::PropertyText, Molecule)>,
        fn(&(cosmolkit_model::PropertyText, Molecule)) -> &Molecule,
    > {
        self.entries.iter().map(|(_, molecule)| molecule)
    }
    pub fn entries(
        &self,
    ) -> std::iter::Map<
        std::slice::Iter<'_, (cosmolkit_model::PropertyText, Molecule)>,
        fn(
            &(cosmolkit_model::PropertyText, Molecule),
        ) -> (&cosmolkit_model::PropertyText, &Molecule),
    > {
        self.entries.iter().map(|(key, molecule)| (key, molecule))
    }
}
impl<'a> IntoIterator for &'a TautomerEnumeration {
    type Item = &'a Molecule;
    type IntoIter = std::iter::Map<
        std::slice::Iter<'a, (cosmolkit_model::PropertyText, Molecule)>,
        for<'b> fn(&'b (cosmolkit_model::PropertyText, Molecule)) -> &'b Molecule,
    >;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
impl std::ops::Index<usize> for TautomerEnumeration {
    type Output = Molecule;
    fn index(&self, index: usize) -> &Molecule {
        &self.entries[index].1
    }
}
/// Select from supplied values without enumerating, preserving their input order.
pub fn canonical_tautomer_from_molecules(
    molecules: &[Molecule],
) -> Result<Molecule, OperationError> {
    canonical_tautomer_from_molecules_with_params(molecules, &TautomerParams::default())
}
pub fn canonical_tautomer_from_molecules_with_params(
    molecules: &[Molecule],
    params: &TautomerParams,
) -> Result<Molecule, OperationError> {
    select_working_candidates(&mut molecules.to_vec(), params)
}

fn select_working_candidates(
    molecules: &mut [Molecule],
    params: &TautomerParams,
) -> Result<Molecule, OperationError> {
    let index = cosmolkit_tautomer::select_canonical_index_by(
        molecules.iter_mut().map(|molecule| ((), molecule)),
        |molecule| score_live_candidate(molecule, params),
        |(), molecule| {
            cosmolkit_tautomer::canonical_smiles(record_view(molecule))
                .map(std::borrow::Cow::Owned)
                .map_err(OperationError::Tautomer)
        },
    )?
    .ok_or_else(|| OperationError::Tautomer(TautomerRunError::NoCanonicalTautomer))?;
    let mut invocation = params.clone();
    invocation.finalize_selected = true;
    molecules[index].canonical_tautomer_with_params(&invocation)
}

/// Project cross-language ownership seam over the sole source selector.
/// Hosts lend shared candidates; one COW working value per input isolates all
/// scorer cache writes, including repeated references to the same molecule.
#[doc(hidden)]
pub fn canonical_tautomer_from_molecule_hosts_with_params(
    len: usize,
    with_candidate: impl Fn(
        usize,
        &mut dyn FnMut(&Molecule) -> Result<(), OperationError>,
    ) -> Result<(), OperationError>,
    params: &TautomerParams,
) -> Result<Molecule, OperationError> {
    let mut working = Vec::with_capacity(len);
    for index in 0..len {
        with_candidate(index, &mut |molecule| {
            working.push(molecule.clone());
            Ok(())
        })?;
    }
    select_working_candidates(&mut working, params)
}

pub(crate) struct EnumerationMetadata {
    pub(crate) keys: Vec<cosmolkit_model::PropertyText>,
    pub(crate) status: TautomerEnumerationStatus,
    pub(crate) modified_atoms: BTreeSet<AtomId>,
    pub(crate) modified_bonds: BTreeSet<BondId>,
}
pub(crate) fn assemble_enumeration(
    molecules: Vec<Molecule>,
    metadata: EnumerationMetadata,
) -> Result<TautomerEnumeration, OperationError> {
    if molecules.len() != metadata.keys.len() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "enumerate_tautomers_with_params",
            field: "retained keys",
            actual: metadata.keys.len(),
            expected: molecules.len(),
        });
    }
    Ok(TautomerEnumeration {
        entries: metadata.keys.into_iter().zip(molecules).collect(),
        status: metadata.status,
        modified_atoms: metadata.modified_atoms,
        modified_bonds: metadata.modified_bonds,
    })
}
fn score_view(
    view: &mut cosmolkit_tautomer::TautomerScoreView<'_>,
    params: &TautomerScoreParams,
) -> Result<TautomerScore, TautomerRunError> {
    match &params.terms {
        Some(terms) => cosmolkit_tautomer::score_tautomer_with_terms_(view, terms),
        None => cosmolkit_tautomer::score_tautomer_(view),
    }
}
impl Molecule {
    /// Score without changing this molecule, including its derived caches.
    pub fn tautomer_score(&self) -> Result<TautomerScore, OperationError> {
        self.tautomer_score_with_params(&Default::default())
    }
    pub fn tautomer_score_with_params(
        &self,
        params: &TautomerScoreParams,
    ) -> Result<TautomerScore, OperationError> {
        if !self
            .derived_cache_runtime()
            .valid_ring_info()
            .is_some_and(cosmolkit_core::RingInfo::is_symm_sssr)
        {
            let mut working = self.clone();
            working.assign_symm_sssr_()?;
            return cosmolkit_tautomer::score_tautomer_from_retained_cache(
                record_view(&working),
                params.terms.as_deref(),
            )
            .map_err(OperationError::Tautomer);
        }
        cosmolkit_tautomer::score_tautomer_from_retained_cache(
            record_view(self),
            params.terms.as_deref(),
        )
        .map_err(OperationError::Tautomer)
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod tests {
    fn fixed_key_text(key: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(key.as_bytes()).expect("original fixed ASCII test observation")
    }

    use super::*;
    #[test]
    #[cfg(feature = "cap-depict")]
    fn tautomer_products_drop_source_coordinates_and_canonical_restores() {
        let mut source = Molecule::from_smiles("CC(C)=O")
            .unwrap()
            .with_2d_coordinates()
            .unwrap();
        let coordinates = source.coordinates_arc_runtime();
        let outputs = source.enumerate_tautomers().unwrap();
        assert!(
            outputs
                .iter()
                .any(|output| Arc::ptr_eq(&coordinates, &output.coordinates_arc_runtime()))
        );
        assert!(
            outputs
                .iter()
                .any(|output| !Arc::ptr_eq(&coordinates, &output.coordinates_arc_runtime()))
        );
        let canonical = source.canonical_tautomer().unwrap();
        assert_eq!(canonical.coordinate_block_runtime().conformers_2d.len(), 2);
        assert_eq!(source.coordinate_block_runtime().conformers_2d.len(), 1);
    }
}

fn record_view(molecule: &Molecule) -> cosmolkit_tautomer::TautomerRecordView<'_> {
    let cache = molecule.derived_cache_runtime();
    cosmolkit_tautomer::TautomerRecordView {
        topology: molecule.topology(),
        coordinates: molecule.coordinate_block_runtime(),
        properties: molecule.properties(),
        valence: cache.valence_assignment(),
        rings: cache.valid_ring_info(),
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod canonical_selection_tests {
    fn fixed_key_text(key: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(key.as_bytes()).expect("original fixed ASCII test observation")
    }

    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    struct ScoreFn<F>(F);
    impl<F> TautomerScorer for ScoreFn<F>
    where
        F: Fn(TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> + Send + Sync,
    {
        fn score(&self, m: TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> {
            (self.0)(m)
        }
    }
    fn configured(
        f: impl Fn(TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> + Send + Sync + 'static,
    ) -> TautomerParams {
        let mut p = TautomerParams::default();
        p.set_scorer(Some(Arc::new(ScoreFn(f))));
        p
    }
    fn fixture(entries: Vec<(&str, Molecule)>) -> TautomerEnumeration {
        TautomerEnumeration {
            entries: entries.into_iter().map(|(k, m)| (k.into(), m)).collect(),
            ..Default::default()
        }
    }
    #[test]
    fn canonical_selection_rejects_empty_inputs_and_unselectable_minimum_scores() {
        assert!(matches!(
            TautomerEnumeration::default().canonical_tautomer(),
            Err(OperationError::Tautomer(
                TautomerRunError::NoCanonicalTautomer
            ))
        ));
        let mut result = fixture(vec![
            ("C", Molecule::from_smiles("C").unwrap()),
            ("CC", Molecule::from_smiles("CC").unwrap()),
        ]);
        assert!(matches!(
            result.canonical_tautomer_with_params(&configured(|_| Ok(i32::MIN))),
            Err(OperationError::Tautomer(
                TautomerRunError::NoCanonicalTautomer
            ))
        ));
    }
    #[test]
    fn canonical_selection_single_item_skips_scoring_and_returns_a_clean_copy() {
        let parsed = Molecule::from_smiles("F[C@H](Cl)Br").unwrap();
        let mut properties = parsed.properties().clone();
        properties.clear_prop("_StereochemDone");
        let mut source = parsed
            .to_builder()
            .with_properties(properties)
            .build()
            .unwrap();
        let before = source.clone();
        let mut result = fixture(vec![("retained-without-rewrite", source.clone())]);
        let calls = Arc::new(AtomicUsize::new(0));
        let seen = calls.clone();
        let selected = result
            .canonical_tautomer_with_params(&configured(move |_| {
                seen.fetch_add(1, Ordering::SeqCst);
                Ok(7)
            }))
            .unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 0);
        assert_eq!(source, before);
        assert_eq!(source.property("_StereochemDone"), None);
        assert_eq!(
            selected.property("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        assert!(
            selected
                .properties()
                .is_prop_computed("_StereochemDone")
                .unwrap()
        );
        let selected = selected
            .to_builder()
            .with_property("selected-only".into(), "yes".into())
            .unwrap()
            .build()
            .unwrap();
        assert_eq!(
            selected.property("selected-only"),
            Some(&cosmolkit_model::PropertyValue::String("yes".into()))
        );
        assert_eq!(source.property("selected-only"), None);
    }
    #[test]
    fn canonical_selection_uses_signed_maximum_and_retained_result_keys() {
        let mut result = fixture(vec![
            ("z-retained", Molecule::from_smiles("C").unwrap()),
            ("a-retained", Molecule::from_smiles("CC").unwrap()),
            ("m-retained", Molecule::from_smiles("CCC").unwrap()),
        ]);
        let negative = result
            .canonical_tautomer_with_params(&configured(|m| {
                Ok(match m.num_atoms() {
                    1 => -10,
                    2 => -1,
                    _ => -5,
                })
            }))
            .unwrap();
        assert_eq!(negative.to_smiles().unwrap().as_bytes(), b"CC");
        let positive = result
            .canonical_tautomer_with_params(&configured(|m| Ok(m.num_atoms() as i32)))
            .unwrap();
        assert_eq!(positive.to_smiles().unwrap().as_bytes(), b"CCC");
        let tied = result
            .canonical_tautomer_with_params(&configured(|_| Ok(4)))
            .unwrap();
        assert_eq!(tied.to_smiles().unwrap().as_bytes(), b"CC");
        assert_eq!(
            result
                .canonical_smiles()
                .into_iter()
                .map(fixed_key_text)
                .collect::<Vec<_>>(),
            ["z-retained", "a-retained", "m-retained"]
        );
    }
    #[test]
    fn canonical_selection_default_and_custom_scorers_share_finalization() {
        let mut source = Molecule::from_smiles("CC(C)=O").unwrap();
        let before = source.clone();
        let mut result = source.enumerate_tautomers().unwrap();
        let a = result.canonical_tautomer().unwrap();
        let b = result
            .canonical_tautomer_with_params(&configured(|mut m| Ok(m.tautomer_score()?.total())))
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(
            a.property("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        assert!(a.properties().is_prop_computed("_StereochemDone").unwrap());
        assert_eq!(source, before);
    }
    #[test]
    fn canonical_result_selection_does_not_enumerate_the_selected_candidate_again() {
        struct NeverEnumerate;
        impl TautomerEnumerationCallback for NeverEnumerate {
            fn should_continue(
                &self,
                _: TautomerMoleculeView<'_>,
                _: TautomerProgress<'_>,
            ) -> Result<bool, TautomerRunError> {
                panic!("pickCanonical does not enumerate")
            }
        }
        let mut source = Molecule::from_smiles("CC(C)=O").unwrap();
        let mut result = source.enumerate_tautomers().unwrap();
        let mut params = TautomerParams::default();
        params.set_callback(Some(Arc::new(NeverEnumerate)));
        assert_eq!(
            result
                .canonical_tautomer_with_params(&params)
                .unwrap()
                .to_smiles()
                .unwrap()
                .as_bytes(),
            b"CC(C)=O"
        );
    }
    #[test]
    fn owned_callback_views_retain_source_and_preapplication_progress_after_the_run() {
        struct Capture(
            std::sync::Mutex<Option<(TautomerMoleculeView<'static>, TautomerProgress<'static>)>>,
        );
        impl TautomerEnumerationCallback for Capture {
            fn should_continue(
                &self,
                m: TautomerMoleculeView<'_>,
                mut p: TautomerProgress<'_>,
            ) -> Result<bool, TautomerRunError> {
                *self.0.lock().unwrap() = Some((m.to_owned(), p.to_owned()));
                Ok(false)
            }
        }
        let capture = Arc::new(Capture(std::sync::Mutex::new(None)));
        let mut params = TautomerParams::default();
        params.set_callback(Some(capture.clone()));
        let mut source = Molecule::from_smiles("CC(C)=O").unwrap();
        let mut result = source.enumerate_tautomers_with_params(&params).unwrap();
        assert_eq!(result.status(), TautomerEnumerationStatus::Canceled);
        drop(result);
        drop(source);
        let mut captured = capture.0.lock().unwrap();
        let (m, p) = captured.as_mut().unwrap();
        assert_eq!(m.to_smiles().unwrap().as_bytes(), b"CC(C)=O");
        assert_eq!(m.tautomer_score().unwrap().total(), 5);
        assert_eq!(p.len(), 1);
        assert_eq!(p.status(), TautomerEnumerationStatus::Completed);
        assert_eq!(p.num_transforms(), 1);
        assert!(p.modified_atoms().is_empty());
        assert!(p.modified_bonds().is_empty());
        let entries = p
            .entries()
            .map(|(key, m)| (key.to_owned(), m.to_smiles().unwrap()))
            .collect::<Vec<_>>();
        assert_eq!(entries, [("CC(C)=O".into(), "CC(C)=O".into())]);
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
#[path = "tautomer_result_tests.rs"]
mod result_tests;

#[cfg(all(test, feature = "cap-smiles"))]
#[path = "tautomer_configuration_tests.rs"]
mod configuration_tests;

impl Clone for TautomerMoleculeView<'_> {
    fn clone(&self) -> Self {
        match &self.storage {
            TautomerViewStorage::Owned(value) => Self {
                storage: TautomerViewStorage::Owned(value.clone()),
            },
            TautomerViewStorage::Borrowed(value) => Self {
                storage: Self::snapshot(value.as_record_view()).storage,
            },
        }
    }
}
impl TautomerMoleculeView<'_> {
    fn snapshot(view: cosmolkit_tautomer::TautomerRecordView<'_>) -> TautomerMoleculeView<'static> {
        TautomerMoleculeView {
            storage: TautomerViewStorage::Owned(Arc::new(TautomerViewSnapshot {
                topology: view.topology.clone(),
                coordinates: view.coordinates.clone(),
                properties: view.properties.clone(),
                valence: view.valence.cloned(),
                rings: view.rings.cloned(),
            })),
        }
    }
    fn reborrow(&mut self) -> TautomerMoleculeView<'_> {
        match &mut self.storage {
            TautomerViewStorage::Borrowed(view) => TautomerMoleculeView::borrowed(view.reborrow()),
            TautomerViewStorage::Owned(shared) => {
                let value = Arc::make_mut(shared);
                let rings = value.rings.get_or_insert_with(|| {
                    cosmolkit_core::RingInfo::new(
                        cosmolkit_core::RingFindType::OtherOrUnknown,
                        value.topology.atoms.len(),
                        value.topology.bonds.len(),
                    )
                });
                let view = cosmolkit_tautomer::TautomerRecordView {
                    topology: &value.topology,
                    coordinates: &value.coordinates,
                    properties: &value.properties,
                    valence: value.valence.as_ref(),
                    rings: None,
                };
                TautomerMoleculeView::borrowed(cosmolkit_tautomer::TautomerScoreView::new(
                    view, rings,
                ))
            }
        }
    }
    pub fn tautomer_score(&mut self) -> Result<TautomerScore, TautomerRunError> {
        let mut loan = self.reborrow();
        match &mut loan.storage {
            TautomerViewStorage::Borrowed(view) => score_view(view, &Default::default()),
            TautomerViewStorage::Owned(_) => {
                unreachable!("reborrow is an exclusive restricted cache loan")
            }
        }
    }
    fn scored_ring_cache(&self) -> Option<&cosmolkit_core::RingInfo> {
        self.record().rings.filter(|rings| rings.is_symm_sssr())
    }
    /// Retain only the completed SymmSSSR cache of this invocation's snapshot.
    /// This is a restricted detached loan, not live Molecule storage authority.
    #[doc(hidden)]
    pub fn retain_score_cache_from(
        &mut self,
        snapshot: &TautomerMoleculeView<'_>,
    ) -> Result<(), TautomerRunError> {
        // Reject a wrong graph before an owned cold target can acquire even
        // an uninitialized placeholder. Ring rows depend on graph identities,
        // not SGroup/conformer floating metadata outside this readonly graph.
        let target = self.record();
        let saved = snapshot.record();
        if target.topology.atoms != saved.topology.atoms
            || target.topology.bonds != saved.topology.bonds
            || target.topology.adjacency != saved.topology.adjacency
        {
            return Err(TautomerRunError::ScoreCacheTargetMismatch);
        }
        let mut loan = self.reborrow();
        match &mut loan.storage {
            TautomerViewStorage::Borrowed(view) => view.retain_scored_cache_from(snapshot.record()),
            TautomerViewStorage::Owned(_) => {
                unreachable!("reborrow is an exclusive restricted cache loan")
            }
        }
    }
}

// This helper only receives COW working candidates, never public inputs.
fn score_live_candidate(
    molecule: &mut Molecule,
    params: &TautomerParams,
) -> Result<i32, OperationError> {
    if params.scorer.is_none() {
        if !molecule
            .derived_cache_runtime()
            .valid_ring_info()
            .is_some_and(cosmolkit_core::RingInfo::is_symm_sssr)
        {
            molecule.assign_symm_sssr_()?;
        }
        return cosmolkit_tautomer::score_tautomer_from_retained_cache(
            record_view(molecule),
            params.score_params.terms.as_deref(),
        )
        .map(TautomerScore::total)
        .map_err(OperationError::Tautomer);
    }
    let original = molecule.derived_cache_runtime().valid_ring_info();
    let was_symm = original.is_some_and(cosmolkit_core::RingInfo::is_symm_sssr);
    let mut rings = original.cloned().unwrap_or_else(|| {
        cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            molecule.num_atoms(),
            molecule.num_bonds(),
        )
    });
    let result = params.score_view(cosmolkit_tautomer::TautomerScoreView::new(
        record_view(molecule),
        &mut rings,
    ));
    // This happens even when the custom scorer returns a later typed error.
    if !was_symm && rings.is_symm_sssr() {
        molecule.install_tautomer_score_cache_(&rings)?;
    }
    result.map_err(OperationError::Tautomer)
}

#[cfg(all(test, feature = "cap-smiles"))]
mod recovery_search06_full_tests {
    use super::*;
    struct AfterScore {
        error: bool,
    }
    impl TautomerEnumerationCallback for AfterScore {
        fn should_continue(
            &self,
            mut source: TautomerMoleculeView<'_>,
            mut progress: TautomerProgress<'_>,
        ) -> Result<bool, TautomerRunError> {
            let _score = source.tautomer_score();
            assert!(
                source.scored_ring_cache().is_some(),
                "ring write precedes any later score component failure"
            );
            for (_, mut candidate) in progress.entries() {
                candidate.tautomer_score()?;
                assert!(candidate.scored_ring_cache().is_some());
            }
            if self.error {
                Err(TautomerRunError::Callback(
                    "controlled failure after completed source score".into(),
                ))
            } else {
                Ok(false)
            }
        }
    }
    fn cold_source() -> Molecule {
        Molecule::from_smiles_with_params(
            "CC(C)=O",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap()
    }
    fn assert_unchanged(actual: &Molecule, before: &Molecule) {
        assert_eq!(actual, before);
        assert!(Arc::ptr_eq(
            &actual.topology_arc_runtime(),
            &before.topology_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &actual.coordinates_arc_runtime(),
            &before.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &actual.properties_arc_runtime(),
            &before.properties_arc_runtime()
        ));
        assert!(Arc::ptr_eq(
            &actual.derived_cache_arc_runtime(),
            &before.derived_cache_arc_runtime()
        ));
    }
    #[test]
    fn recovery_search06_value_success_preserves_all_input_blocks() {
        let source = Molecule::from_smiles("CC(C)=O").unwrap();
        let before = source.clone();
        assert_eq!(source.tautomer_score().unwrap().total(), 5);
        assert_unchanged(&source, &before);
        let result = source.enumerate_tautomers().unwrap();
        assert_eq!(result.len(), 2);
        assert_unchanged(&source, &before);
        let canonical = source.canonical_tautomer().unwrap();
        assert_eq!(canonical.to_smiles().unwrap().as_bytes(), b"CC(C)=O");
        assert_unchanged(&source, &before);
    }
    #[test]
    fn recovery_search06_cold_score_preserves_all_input_blocks() {
        let source = Molecule::from_smiles_with_params(
            "c1ccccc1",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let before = source.clone();
        assert!(source.derived_cache_runtime().valid_ring_info().is_none());
        assert_eq!(source.tautomer_score().unwrap().ring(), 250);
        assert_unchanged(&source, &before);
    }
    #[test]
    fn recovery_search06_callback_error_preserves_source_and_shared_cache() {
        let source = cold_source();
        let peer = source.clone();
        let cache = source.derived_cache_arc_runtime();
        assert!(source.derived_cache_runtime().valid_ring_info().is_none());
        let topology = source.topology().clone();
        let properties = source.properties().clone();
        let coordinates = source.coordinates_arc_runtime();
        let mut params = TautomerParams::default();
        params.set_callback(Some(Arc::new(AfterScore { error: true })));
        assert!(matches!(
            source.enumerate_tautomers_with_params(&params),
            Err(OperationError::Tautomer(TautomerRunError::Callback(_)))
        ));
        assert!(source.derived_cache_runtime().valid_ring_info().is_none());
        assert!(Arc::ptr_eq(&cache, &source.derived_cache_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &peer.derived_cache_arc_runtime()));
        assert!(peer.derived_cache_runtime().valid_ring_info().is_none());
        assert_eq!(source.topology(), &topology);
        assert_eq!(source.properties(), &properties);
        assert!(Arc::ptr_eq(&coordinates, &source.coordinates_arc_runtime()));
        assert_unchanged(&source, &peer);
    }
    #[test]
    fn recovery_search06_callback_cancel_preserves_source_cache_and_status() {
        let source = cold_source();
        let before = source.clone();
        let cache = source.derived_cache_arc_runtime();
        let mut params = TautomerParams::default();
        params.set_callback(Some(Arc::new(AfterScore { error: false })));
        let result = source.enumerate_tautomers_with_params(&params).unwrap();
        assert_eq!(result.status(), TautomerEnumerationStatus::Canceled);
        assert_eq!(result.len(), 1);
        assert_eq!(source, before);
        assert!(source.derived_cache_runtime().valid_ring_info().is_none());
        assert!(Arc::ptr_eq(&cache, &source.derived_cache_arc_runtime()));
    }
    #[test]
    fn recovery_search06_no_score_callback_keeps_source_cold() {
        struct Cancel;
        impl TautomerEnumerationCallback for Cancel {
            fn should_continue(
                &self,
                _: TautomerMoleculeView<'_>,
                _: TautomerProgress<'_>,
            ) -> Result<bool, TautomerRunError> {
                Ok(false)
            }
        }
        let mut source = cold_source();
        let mut params = TautomerParams::default();
        params.set_callback(Some(Arc::new(Cancel)));
        source.enumerate_tautomers_with_params(&params).unwrap();
        assert!(
            source.derived_cache_runtime().valid_ring_info().is_none(),
            "no eager cache materialization"
        );
    }
    #[test]
    fn recovery_search06_owned_snapshot_score_retains_cache_and_cow_peer() {
        let source = Molecule::from_smiles_with_params(
            "c1ccccc1",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let mut snapshot = TautomerMoleculeView::snapshot(record_view(&source));
        let peer = snapshot.clone();
        assert!(snapshot.record().rings.is_none());
        assert_eq!(snapshot.tautomer_score().unwrap().ring(), 250);
        assert!(snapshot.scored_ring_cache().is_some());
        assert!(peer.record().rings.is_none());
        assert!(source.derived_cache_runtime().valid_ring_info().is_none());
    }
    #[test]
    fn recovery_search06_working_cache_rejects_invalid_ring_transport() {
        for uninitialized in [false, true] {
            let source = cold_source();
            let mut working = source.clone();
            let cache = source.derived_cache_arc_runtime();
            let mut invalid = cosmolkit_core::RingInfo::new(
                cosmolkit_core::RingFindType::OtherOrUnknown,
                source.num_atoms(),
                source.num_bonds(),
            );
            // RingInfo::new calls initialize even for OtherOrUnknown. Exercise
            // that original invalid state and a truly reset source cache.
            if uninitialized {
                invalid.reset();
            }
            assert_eq!(invalid.is_initialized(), !uninitialized);
            let expected_field = if uninitialized {
                "initialized-ring-info"
            } else {
                "symm-sssr-ring-info"
            };
            let result = working.install_tautomer_score_cache_(&invalid);
            assert!(
                matches!(
                    result,
                    Err(OperationError::InvalidAlgorithmResult {
                        operation: "with_installed_tautomer_score_cache",
                        field,
                        actual: 0,
                        expected: 1,
                    }) if field == expected_field
                ),
                "actual typed host result: {result:?}"
            );
            assert!(source.derived_cache_runtime().valid_ring_info().is_none());
            assert!(Arc::ptr_eq(&cache, &source.derived_cache_arc_runtime()));
        }
        // Controlled invalid internal cache transport, not a natural source molecule.
    }
    #[test]
    fn recovery_search06_warm_live_score_keeps_cache_arc() {
        for (text, expected) in [("c1ccccc1", 250), ("CC", 0)] {
            let mut source = Molecule::from_smiles(text).unwrap();
            let cache = source.derived_cache_arc_runtime();
            assert!(
                source
                    .derived_cache_runtime()
                    .valid_ring_info()
                    .unwrap()
                    .is_symm_sssr()
            );
            assert_eq!(source.tautomer_score().unwrap().ring(), expected);
            assert!(Arc::ptr_eq(&cache, &source.derived_cache_arc_runtime()));
        }
    }
}
#[cfg(all(test, feature = "cap-smiles"))]
mod recovery_search06_selection_cache_tests {
    use super::*;
    struct ScoreThenFail;
    impl TautomerScorer for ScoreThenFail {
        fn score(&self, mut source: TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> {
            source.tautomer_score()?;
            Err(TautomerRunError::Callback(
                "controlled error after score".into(),
            ))
        }
    }
    fn cold(text: &str) -> Molecule {
        Molecule::from_smiles_with_params(
            text,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap()
    }
    #[test]
    fn recovery_search06_iterable_custom_error_preserves_all_inputs() {
        let inputs = [cold("CC"), cold("CCC")];
        let before = inputs.clone();
        let mut params = TautomerParams::default();
        params.set_scorer(Some(Arc::new(ScoreThenFail)));
        assert!(matches!(
            canonical_tautomer_from_molecules_with_params(&inputs, &params),
            Err(OperationError::Tautomer(TautomerRunError::Callback(_)))
        ));
        assert!(
            inputs[0]
                .derived_cache_runtime()
                .valid_ring_info()
                .is_none()
        );
        assert!(
            inputs[1]
                .derived_cache_runtime()
                .valid_ring_info()
                .is_none()
        );
        for (actual, original) in inputs.iter().zip(&before) {
            assert_eq!(actual, original);
            assert!(Arc::ptr_eq(
                &actual.derived_cache_arc_runtime(),
                &original.derived_cache_arc_runtime(),
            ));
            assert_eq!(actual.topology(), original.topology());
            assert_eq!(actual.properties(), original.properties());
        }
    }
    #[test]
    fn recovery_search06_result_custom_error_preserves_all_candidates() {
        let result = TautomerEnumeration {
            entries: vec![("CC".into(), cold("CC")), ("CCC".into(), cold("CCC"))],
            ..Default::default()
        };
        let before = result.clone();
        let mut params = TautomerParams::default();
        params.set_scorer(Some(Arc::new(ScoreThenFail)));
        assert!(matches!(
            result.canonical_tautomer_with_params(&params),
            Err(OperationError::Tautomer(TautomerRunError::Callback(_)))
        ));
        assert_eq!(result, before);
        assert!(
            result.entries[0]
                .1
                .derived_cache_runtime()
                .valid_ring_info()
                .is_none()
        );
        assert!(
            result.entries[1]
                .1
                .derived_cache_runtime()
                .valid_ring_info()
                .is_none()
        );
        // Controlled retained-result state, not a naturally returned cold enumeration.
    }
    #[test]
    fn recovery_search06_single_iterable_candidate_skips_scorer_and_cache_write() {
        let mut inputs = [cold("CC")];
        let mut params = TautomerParams::default();
        params.set_scorer(Some(Arc::new(ScoreThenFail)));
        canonical_tautomer_from_molecules_with_params(&mut inputs, &params).unwrap();
        assert!(
            inputs[0]
                .derived_cache_runtime()
                .valid_ring_info()
                .is_none()
        );
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod recovery_search06_transport_host_tests {
    use super::*;
    fn cold(text: &str) -> Molecule {
        Molecule::from_smiles_with_params(
            text,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap()
    }
    #[test]
    fn recovery_search06_transport_updates_only_corresponding_loan() {
        let source = cold("c1ccccc1");
        let mut rings = cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            source.num_atoms(),
            source.num_bonds(),
        );
        let mut loan = TautomerMoleculeView::borrowed(cosmolkit_tautomer::TautomerScoreView::new(
            record_view(&source),
            &mut rings,
        ));
        let mut snapshot = loan.to_owned();
        snapshot.tautomer_score().unwrap();
        loan.retain_score_cache_from(&snapshot).unwrap();
        assert!(loan.scored_ring_cache().is_some());
        assert!(source.derived_cache_runtime().valid_ring_info().is_none());
    }
    #[test]
    fn recovery_search06_transport_rejects_wrong_graph_before_cache_write() {
        let source = cold("CC");
        let other = cold("c1ccccc1");
        let mut target = TautomerMoleculeView::snapshot(record_view(&source));
        let mut snapshot = TautomerMoleculeView::snapshot(record_view(&other));
        snapshot.tautomer_score().unwrap();
        assert!(matches!(
            target.retain_score_cache_from(&snapshot),
            Err(TautomerRunError::ScoreCacheTargetMismatch)
        ));
        assert!(target.record().rings.is_none());
    }
    #[test]
    fn recovery_search06_single_slice_control_skips_scoring() {
        let mut molecule = cold("CC");
        let params = TautomerParams::default();
        // Fn host may be called once for final copy; no simultaneous mutable
        // host loans are required. A direct mutable slice exercises the same
        // sole selector's one-candidate rule, already shared by host dispatch.
        canonical_tautomer_from_molecules_with_params(std::slice::from_mut(&mut molecule), &params)
            .unwrap();
        assert!(molecule.derived_cache_runtime().valid_ring_info().is_none());
    }
}
