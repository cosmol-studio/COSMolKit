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

/// Terms used by the read-only score query. None selects the source defaults.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TautomerScoreParams {
    pub terms: Option<Vec<TautomerScoreTerm>>,
}

/// Immutable detached inspection during a callback or custom score.
#[derive(Clone)]
pub struct TautomerMoleculeView<'a> {
    storage: TautomerViewStorage<'a>,
}
#[derive(Clone)]
enum TautomerViewStorage<'a> {
    Borrowed(cosmolkit_tautomer::TautomerRecordView<'a>),
    Owned(Arc<TautomerViewSnapshot>),
}
struct TautomerViewSnapshot {
    topology: cosmolkit_model::TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    valence: Option<cosmolkit_core::ValenceAssignment>,
    rings: Option<cosmolkit_core::RingInfo>,
}
impl<'a> TautomerMoleculeView<'a> {
    fn borrowed(view: cosmolkit_tautomer::TautomerRecordView<'a>) -> Self {
        Self {
            storage: TautomerViewStorage::Borrowed(view),
        }
    }
    fn record(&self) -> cosmolkit_tautomer::TautomerRecordView<'_> {
        match &self.storage {
            TautomerViewStorage::Borrowed(view) => *view,
            TautomerViewStorage::Owned(value) => cosmolkit_tautomer::TautomerRecordView {
                topology: &value.topology,
                coordinates: &value.coordinates,
                properties: &value.properties,
                valence: value.valence.as_ref(),
                rings: value.rings.as_ref(),
            },
        }
    }
    /// Retain an immutable view after the callback without any live runtime authority.
    pub fn to_owned(&self) -> TautomerMoleculeView<'static> {
        let storage = match &self.storage {
            TautomerViewStorage::Owned(value) => TautomerViewStorage::Owned(value.clone()),
            TautomerViewStorage::Borrowed(view) => {
                TautomerViewStorage::Owned(Arc::new(TautomerViewSnapshot {
                    topology: view.topology.clone(),
                    coordinates: view.coordinates.clone(),
                    properties: view.properties.clone(),
                    valence: view.valence.cloned(),
                    rings: view.rings.cloned(),
                }))
            }
        };
        TautomerMoleculeView { storage }
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
            },
        )?)
    }
    pub fn tautomer_score(&self) -> Result<TautomerScore, TautomerRunError> {
        cosmolkit_tautomer::score_tautomer(self.record())
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
        TautomerProgress {
            storage: TautomerProgressStorage::Owned {
                entries: self
                    .entries()
                    .map(|(key, value)| (key.to_owned(), value.to_owned()))
                    .collect(),
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
        &self,
    ) -> Box<
        dyn ExactSizeIterator<Item = (&cosmolkit_model::PropertyText, TautomerMoleculeView<'_>)>
            + '_,
    > {
        match &self.storage {
            TautomerProgressStorage::Borrowed { inner, coordinates } => {
                Box::new(inner.entries().map(|(key, value)| {
                    (key, TautomerMoleculeView::borrowed(value.view(coordinates)))
                }))
            }
            TautomerProgressStorage::Owned { entries, .. } => {
                Box::new(entries.iter().map(|(key, value)| (key, value.clone())))
            }
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
        view: cosmolkit_tautomer::TautomerRecordView<'_>,
    ) -> Result<i32, TautomerRunError> {
        match &self.scorer {
            Some(scorer) => scorer.score(TautomerMoleculeView::borrowed(view)),
            None => score_view(view, &self.score_params).map(TautomerScore::total),
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

pub(crate) struct CallbackAdapter<'a>(pub(crate) &'a TautomerParams);
impl cosmolkit_tautomer::TautomerEnumerationCallback for CallbackAdapter<'_> {
    fn should_continue(
        &self,
        source: cosmolkit_tautomer::TautomerRecordView<'_>,
        progress: cosmolkit_tautomer::TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError> {
        match &self.0.callback {
            Some(callback) => callback.should_continue(
                TautomerMoleculeView::borrowed(source),
                TautomerProgress {
                    storage: TautomerProgressStorage::Borrowed {
                        inner: progress,
                        coordinates: source.coordinates,
                    },
                },
            ),
            None => Ok(true),
        }
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
        let index = cosmolkit_tautomer::select_canonical_index_with(
            self.entries
                .iter()
                .map(|(key, molecule)| (key, record_view(molecule))),
            |view| params.score_view(view),
        )
        .map_err(OperationError::Tautomer)?;
        let mut invocation = params.clone();
        invocation.finalize_selected = true;
        self.entries[index]
            .1
            .canonical_tautomer_with_params(&invocation)
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
    let index = cosmolkit_tautomer::select_canonical_index_from_iterable_with(
        molecules.iter().map(record_view),
        |view| params.score_view(view),
    )
    .map_err(OperationError::Tautomer)?;
    let mut invocation = params.clone();
    invocation.finalize_selected = true;
    molecules[index].canonical_tautomer_with_params(&invocation)
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
    view: cosmolkit_tautomer::TautomerRecordView<'_>,
    params: &TautomerScoreParams,
) -> Result<TautomerScore, TautomerRunError> {
    match &params.terms {
        Some(terms) => cosmolkit_tautomer::score_tautomer_with_terms(view, terms),
        None => cosmolkit_tautomer::score_tautomer(view),
    }
}
impl Molecule {
    pub fn tautomer_score(&self) -> Result<TautomerScore, TautomerRunError> {
        self.tautomer_score_with_params(&Default::default())
    }
    pub fn tautomer_score_with_params(
        &self,
        params: &TautomerScoreParams,
    ) -> Result<TautomerScore, TautomerRunError> {
        let cache = self.derived_cache_runtime();
        score_view(
            cosmolkit_tautomer::TautomerRecordView {
                topology: self.topology(),
                coordinates: self.coordinate_block_runtime(),
                properties: self.properties(),
                valence: cache.valence_assignment(),
                rings: cache.valid_ring_info(),
            },
            params,
        )
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod tests {
    fn fixed_key_text(key: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(key.as_bytes()).expect("original fixed ASCII test observation")
    }

    use super::*;
    #[test]
    fn tautomer_outputs_share_the_unchanged_source_coordinate_block() {
        let source = Molecule::from_smiles("CC(C)=O").unwrap();
        let coordinates = source.coordinates_arc_runtime();
        for output in source.enumerate_tautomers().unwrap().iter() {
            assert!(Arc::ptr_eq(&coordinates, &output.coordinates_arc_runtime()));
        }
        assert!(Arc::ptr_eq(
            &coordinates,
            &source
                .canonical_tautomer()
                .unwrap()
                .coordinates_arc_runtime()
        ));
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
        let result = fixture(vec![
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
        let source = parsed
            .to_builder()
            .with_properties(properties)
            .build()
            .unwrap();
        let before = source.clone();
        let result = fixture(vec![("retained-without-rewrite", source.clone())]);
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
        let result = fixture(vec![
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
        let source = Molecule::from_smiles("CC(C)=O").unwrap();
        let before = source.clone();
        let result = source.enumerate_tautomers().unwrap();
        let a = result.canonical_tautomer().unwrap();
        let b = result
            .canonical_tautomer_with_params(&configured(|m| Ok(m.tautomer_score()?.total())))
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
        let source = Molecule::from_smiles("CC(C)=O").unwrap();
        let result = source.enumerate_tautomers().unwrap();
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
                p: TautomerProgress<'_>,
            ) -> Result<bool, TautomerRunError> {
                *self.0.lock().unwrap() = Some((m.to_owned(), p.to_owned()));
                Ok(false)
            }
        }
        let capture = Arc::new(Capture(std::sync::Mutex::new(None)));
        let mut params = TautomerParams::default();
        params.set_callback(Some(capture.clone()));
        let source = Molecule::from_smiles("CC(C)=O").unwrap();
        let result = source.enumerate_tautomers_with_params(&params).unwrap();
        assert_eq!(result.status(), TautomerEnumerationStatus::Canceled);
        drop(result);
        drop(source);
        let captured = capture.0.lock().unwrap();
        let (m, p) = captured.as_ref().unwrap();
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
