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
#[derive(Clone, Copy)]
pub struct TautomerMoleculeView<'a> {
    pub(crate) inner: cosmolkit_tautomer::TautomerRecordView<'a>,
}
impl TautomerMoleculeView<'_> {
    pub fn num_atoms(&self) -> usize {
        self.inner.topology.atoms.len()
    }
    pub fn num_bonds(&self) -> usize {
        self.inner.topology.bonds.len()
    }
    pub fn atoms(&self) -> &[Atom] {
        &self.inner.topology.atoms
    }
    pub fn bonds(&self) -> &[Bond] {
        &self.inner.topology.bonds
    }
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.inner.topology.atoms.get(id.index())
    }
    pub fn bond(&self, id: BondId) -> Option<&Bond> {
        self.inner.topology.bonds.get(id.index())
    }
    pub fn properties(&self) -> &MoleculeProperties {
        self.inner.properties
    }
    pub fn to_smiles(&self) -> Result<String, TautomerRunError> {
        Ok(cosmolkit_smiles::write_smiles(
            cosmolkit_smiles::SmilesRecordView {
                topology: self.inner.topology,
                coordinates: self.inner.coordinates,
                properties: self.inner.properties,
            },
        )?)
    }
    pub fn tautomer_score(&self) -> Result<TautomerScore, TautomerRunError> {
        cosmolkit_tautomer::score_tautomer(self.inner)
    }
}

/// Borrowed source-ordered progress before application of a matched transform.
pub struct TautomerProgress<'a> {
    inner: cosmolkit_tautomer::TautomerProgress<'a>,
    coordinates: &'a CoordinateBlock,
}
impl TautomerProgress<'_> {
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn status(&self) -> TautomerEnumerationStatus {
        self.inner.status()
    }
    pub fn num_transforms(&self) -> u32 {
        self.inner.num_transforms()
    }
    pub fn modified_atoms(&self) -> &BTreeSet<AtomId> {
        self.inner.modified_atoms()
    }
    pub fn modified_bonds(&self) -> &BTreeSet<BondId> {
        self.inner.modified_bonds()
    }
    pub fn entries(&self) -> impl ExactSizeIterator<Item = (&str, TautomerMoleculeView<'_>)> {
        self.inner.entries().map(|(key, value)| {
            (
                key,
                TautomerMoleculeView {
                    inner: value.view(self.coordinates),
                },
            )
        })
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
            Some(scorer) => scorer.score(TautomerMoleculeView { inner: view }),
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
                TautomerMoleculeView { inner: source },
                TautomerProgress {
                    inner: progress,
                    coordinates: source.coordinates,
                },
            ),
            None => Ok(true),
        }
    }
}
/// Rich result containing only runtime-validated molecule values.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TautomerEnumeration {
    entries: Vec<(String, Molecule)>,
    status: TautomerEnumerationStatus,
    modified_atoms: BTreeSet<AtomId>,
    modified_bonds: BTreeSet<BondId>,
}
impl TautomerEnumeration {
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
    pub fn canonical_smiles(&self) -> Vec<&str> {
        self.entries.iter().map(|(k, _)| k.as_str()).collect()
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &Molecule> {
        self.entries.iter().map(|(_, m)| m)
    }
    pub fn entries(&self) -> impl ExactSizeIterator<Item = (&str, &Molecule)> {
        self.entries.iter().map(|(k, m)| (k.as_str(), m))
    }
}
impl std::ops::Index<usize> for TautomerEnumeration {
    type Output = Molecule;
    fn index(&self, index: usize) -> &Molecule {
        &self.entries[index].1
    }
}
pub(crate) struct EnumerationMetadata {
    pub(crate) keys: Vec<String>,
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
