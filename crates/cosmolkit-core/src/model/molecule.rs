use std::{
    collections::BTreeMap,
    sync::{Arc, RwLock},
};

use crate::{
    Atom, AtomId, Bond, MoleculeBuilder, error::InvariantError, invariants::enforce_molecule_invariants,
    sgroup::SubstanceGroup, stereo::StereoGroup,
};

pub use cosmolkit_model::{
    Conformer2D, Conformer3D, ConformerStore, CoordinateBlock, CoordinateDimension, MoleculeProperties, PropertyStore,
    SdfPropertyList, SdfPropertyListTarget,
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmilesParseError {
    #[error("SMILES parser is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("{0}")]
    ParseError(String),
}

pub use crate::smiles_write::SmilesWriteError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopologyTrust {
    Unknown,
    TrustedGraph,
    CoordinateOnly,
}

impl Default for TopologyTrust {
    fn default() -> Self {
        Self::Unknown
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MoleculeCapabilities {
    topology_trust: TopologyTrust,
}

impl MoleculeCapabilities {
    #[must_use]
    pub const fn new(topology_trust: TopologyTrust) -> Self {
        Self { topology_trust }
    }

    #[must_use]
    pub const fn topology_trust(self) -> TopologyTrust {
        self.topology_trust
    }

    #[must_use]
    pub const fn with_topology_trust(self, topology_trust: TopologyTrust) -> Self {
        Self { topology_trust }
    }
}

pub type TopologyBlock = cosmolkit_model::TopologyBlock<
    crate::QueryNode<crate::AtomQueryPredicate>,
    crate::QueryNode<crate::BondQueryPredicate>,
>;
pub use cosmolkit_model::{AtomMapping, BondMapping, TopologyMapping};

#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct DerivedCacheBlock {
    pub(crate) rings: Option<crate::RingInfo>,
    pub(crate) ring_families: Option<crate::RingInfo>,
    pub(crate) valence: Option<crate::ValenceAssignment>,
    pub(crate) aromaticity_valid: bool,
    pub(crate) stereo_valid: bool,
}

impl DerivedCacheBlock {
    pub(crate) fn invalidate(&mut self, states: crate::DerivedState) {
        if states.contains(crate::DerivedState::VALENCE) {
            self.valence = None;
        }
        if states.contains(crate::DerivedState::RINGS) {
            self.rings = None;
        }
        if states.contains(crate::DerivedState::RING_FAMILIES) {
            self.ring_families = None;
        }
        if states.contains(crate::DerivedState::AROMATICITY) {
            self.aromaticity_valid = false;
        }
        if states.contains(crate::DerivedState::STEREO) {
            self.stereo_valid = false;
        }
    }
}

/// RDKit-compatible computed properties written through logically read-only
/// descriptor APIs. This cache is separate from the Arc-backed operation
/// cache: cloning a molecule copies the current values without coupling later
/// descriptor calls between the copies.
#[derive(Debug, Default)]
struct ComputedPropertyCache {
    crippen: RwLock<Option<[f64; 2]>>,
    crippen_atom_contributions: RwLock<Option<CrippenAtomContributionCache>>,
    connectivity_hk_deltas: RwLock<Option<Arc<[f64]>>>,
    connectivity_n_vals: RwLock<Option<Arc<[f64]>>>,
    labute: RwLock<Option<LabuteDescriptorCache>>,
    topological_distance_matrix: RwLock<Option<Arc<[f64]>>>,
    distance_matrices_3d: RwLock<BTreeMap<usize, Arc<[f64]>>>,
    potential_stereo: RwLock<Option<Arc<[crate::StereoInfo]>>>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CrippenAtomContributionCache {
    pub(crate) logp: Arc<[f64]>,
    pub(crate) mr: Arc<[f64]>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct LabuteDescriptorCache {
    pub(crate) atom_contributions: Arc<[f64]>,
    pub(crate) hydrogen_contribution: f64,
    pub(crate) asa: f64,
}

impl Clone for ComputedPropertyCache {
    fn clone(&self) -> Self {
        Self {
            crippen: RwLock::new(self.crippen()),
            crippen_atom_contributions: RwLock::new(self.crippen_atom_contributions()),
            connectivity_hk_deltas: RwLock::new(self.connectivity_hk_deltas()),
            connectivity_n_vals: RwLock::new(self.connectivity_n_vals()),
            labute: RwLock::new(self.labute()),
            topological_distance_matrix: RwLock::new(self.topological_distance_matrix()),
            distance_matrices_3d: RwLock::new(self.distance_matrices_3d()),
            potential_stereo: RwLock::new(self.potential_stereo()),
        }
    }
}

impl PartialEq for ComputedPropertyCache {
    fn eq(&self, _other: &Self) -> bool {
        // Computed properties are observational caches, not molecule state.
        true
    }
}

impl ComputedPropertyCache {
    fn crippen(&self) -> Option<[f64; 2]> {
        *self.crippen.read().unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    fn set_crippen(&self, values: [f64; 2]) {
        *self.crippen.write().unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn crippen_atom_contributions(&self) -> Option<CrippenAtomContributionCache> {
        self.crippen_atom_contributions
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn set_crippen_atom_contributions(&self, values: CrippenAtomContributionCache) {
        *self
            .crippen_atom_contributions
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn connectivity_hk_deltas(&self) -> Option<Arc<[f64]>> {
        self.connectivity_hk_deltas
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn set_connectivity_hk_deltas(&self, values: Arc<[f64]>) {
        *self
            .connectivity_hk_deltas
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn connectivity_n_vals(&self) -> Option<Arc<[f64]>> {
        self.connectivity_n_vals
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn set_connectivity_n_vals(&self, values: Arc<[f64]>) {
        *self
            .connectivity_n_vals
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn labute(&self) -> Option<LabuteDescriptorCache> {
        self.labute
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn set_labute(&self, values: LabuteDescriptorCache) {
        *self.labute.write().unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn topological_distance_matrix(&self) -> Option<Arc<[f64]>> {
        self.topological_distance_matrix
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn topological_distance_matrix_or_init(&self, initialize: impl FnOnce() -> Vec<f64>) -> Arc<[f64]> {
        if let Some(matrix) = self.topological_distance_matrix() {
            return matrix;
        }
        let mut cached = self
            .topological_distance_matrix
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(matrix) = cached.as_ref() {
            return Arc::clone(matrix);
        }
        let matrix = Arc::<[f64]>::from(initialize());
        *cached = Some(Arc::clone(&matrix));
        matrix
    }

    fn distance_matrices_3d(&self) -> BTreeMap<usize, Arc<[f64]>> {
        self.distance_matrices_3d
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn potential_stereo(&self) -> Option<Arc<[crate::StereoInfo]>> {
        self.potential_stereo
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn set_potential_stereo(&self, values: Arc<[crate::StereoInfo]>) {
        *self
            .potential_stereo
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn distance_matrix_3d_or_init(&self, conformer_id: usize, initialize: impl FnOnce() -> Vec<f64>) -> Arc<[f64]> {
        if let Some(matrix) = self
            .distance_matrices_3d
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .get(&conformer_id)
            .cloned()
        {
            return matrix;
        }
        let mut cached = self
            .distance_matrices_3d
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(matrix) = cached.get(&conformer_id) {
            return Arc::clone(matrix);
        }
        let matrix = Arc::<[f64]>::from(initialize());
        cached.insert(conformer_id, Arc::clone(&matrix));
        matrix
    }

    fn clear_coordinate_dependent(&self) {
        self.distance_matrices_3d
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clear();
    }

    fn clear(&self) {
        *self.crippen.write().unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self
            .crippen_atom_contributions
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self
            .connectivity_hk_deltas
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self
            .connectivity_n_vals
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self.labute.write().unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self
            .topological_distance_matrix
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self
            .potential_stereo
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        self.clear_coordinate_dependent();
    }
}

///
/// The only public way to create a molecule is `MoleculeBuilder`.
///
/// Existing molecules are transformed through registered operations. This type
/// intentionally does not expose mutable storage accessors.
///
/// # Examples
///
/// Build a molecule from SMILES and create a transformed value:
///
/// ```
/// use cosmolkit_core::Molecule;
///
/// let mol = Molecule::from_smiles("CCO").unwrap();
/// let named = mol.with_name("ethanol");
///
/// assert_eq!(mol.properties().name(), None);
/// assert_eq!(named.properties().name(), Some("ethanol"));
/// ```
///
/// Use explicit in-place operations when mutation is intended:
///
/// ```
/// use cosmolkit_core::Molecule;
///
/// let mut mol = Molecule::from_smiles("c1ccccc1").unwrap();
/// mol.kekulize_(true).unwrap();
/// mol.sanitize_().unwrap();
///
/// assert_eq!(mol.num_atoms(), 6);
/// ```
#[derive(Debug, PartialEq, Default)]
pub struct Molecule {
    topology: Arc<TopologyBlock>,
    coordinates: Arc<CoordinateBlock>,
    properties: Arc<MoleculeProperties>,
    derived_cache: Arc<DerivedCacheBlock>,
    capabilities: Arc<MoleculeCapabilities>,
    computed_properties: ComputedPropertyCache,
}

impl Clone for Molecule {
    fn clone(&self) -> Self {
        Self {
            topology: Arc::clone(&self.topology),
            coordinates: Arc::clone(&self.coordinates),
            properties: Arc::clone(&self.properties),
            derived_cache: Arc::clone(&self.derived_cache),
            capabilities: Arc::clone(&self.capabilities),
            computed_properties: self.computed_properties.clone(),
        }
    }
}

impl Molecule {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn builder() -> MoleculeBuilder {
        MoleculeBuilder::new()
    }

    pub(crate) fn from_blocks(
        mut topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, InvariantError> {
        topology.adjacency = crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
        let molecule = Self {
            topology: Arc::new(topology),
            coordinates: Arc::new(coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::new(DerivedCacheBlock::default()),
            capabilities: Arc::new(MoleculeCapabilities::default()),
            computed_properties: ComputedPropertyCache::default(),
        };
        enforce_molecule_invariants(&molecule)?;
        Ok(molecule)
    }

    pub(crate) fn from_blocks_with_capabilities(
        mut topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        capabilities: MoleculeCapabilities,
    ) -> Result<Self, InvariantError> {
        topology.adjacency = crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
        let molecule = Self {
            topology: Arc::new(topology),
            coordinates: Arc::new(coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::new(DerivedCacheBlock::default()),
            capabilities: Arc::new(capabilities),
            computed_properties: ComputedPropertyCache::default(),
        };
        enforce_molecule_invariants(&molecule)?;
        Ok(molecule)
    }

    pub(crate) fn from_operation_blocks(
        mut topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        derived_cache: DerivedCacheBlock,
        capabilities: MoleculeCapabilities,
    ) -> Result<Self, InvariantError> {
        topology.adjacency = crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
        let molecule = Self {
            topology: Arc::new(topology),
            coordinates: Arc::new(coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::new(derived_cache),
            capabilities: Arc::new(capabilities),
            computed_properties: ComputedPropertyCache::default(),
        };
        enforce_molecule_invariants(&molecule)?;
        Ok(molecule)
    }

    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.topology.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.topology.bonds
    }

    #[must_use]
    pub(crate) fn adjacency(&self) -> &crate::AdjacencyList {
        &self.topology.adjacency
    }

    #[must_use]
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms().get(id.index())
    }

    #[must_use]
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.atoms().iter().map(Atom::atomic_number).collect()
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms().len()
    }

    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.bonds().len()
    }

    #[must_use]
    pub fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.coordinates.conformers_2d.first().map(Conformer2D::coordinates)
    }

    #[must_use]
    pub fn conformers_2d(&self) -> &[Conformer2D] {
        &self.coordinates.conformers_2d
    }

    #[must_use]
    pub fn conformers_3d(&self) -> &[Conformer3D] {
        &self.coordinates.conformers_3d
    }

    #[must_use]
    pub fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.coordinates.source_coordinate_dim
    }

    #[must_use]
    pub fn capabilities(&self) -> MoleculeCapabilities {
        *self.capabilities
    }

    #[must_use]
    pub fn topology_trust(&self) -> TopologyTrust {
        self.capabilities.topology_trust()
    }

    #[must_use]
    pub fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.topology.substance_groups
    }

    #[must_use]
    pub fn stereo_groups(&self) -> &[StereoGroup] {
        &self.topology.stereo_groups
    }

    #[must_use]
    pub fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.properties.prop(key)
    }

    /// Returns whether a molecule property is registered as computed state.
    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.properties.is_prop_computed(key)
    }

    #[must_use]
    pub fn with_name(&self, name: impl Into<String>) -> Self {
        let mut properties = (*self.properties).clone();
        properties = properties.with_name(name);
        Self {
            topology: Arc::clone(&self.topology),
            coordinates: Arc::clone(&self.coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::clone(&self.derived_cache),
            capabilities: Arc::clone(&self.capabilities),
            computed_properties: self.computed_properties.clone(),
        }
    }

    #[must_use]
    pub fn with_prop(&self, key: impl Into<String>, value: impl Into<String>) -> Self {
        let mut properties = (*self.properties).clone();
        properties = properties.with_prop(key, value);
        Self {
            topology: Arc::clone(&self.topology),
            coordinates: Arc::clone(&self.coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::clone(&self.derived_cache),
            capabilities: Arc::clone(&self.capabilities),
            computed_properties: self.computed_properties.clone(),
        }
    }

    #[must_use]
    pub fn with_sdf_data_field(&self, key: impl Into<String>, value: impl Into<String>) -> Self {
        let mut properties = (*self.properties).clone();
        properties = properties.with_sdf_data_field(key, value);
        Self {
            topology: Arc::clone(&self.topology),
            coordinates: Arc::clone(&self.coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::clone(&self.derived_cache),
            capabilities: Arc::clone(&self.capabilities),
            computed_properties: self.computed_properties.clone(),
        }
    }

    pub fn from_smiles(smiles: &str) -> Result<Self, SmilesParseError> {
        crate::smiles::mol_from_smiles(smiles, &crate::smiles::SmilesParseParams::default())
    }

    pub fn from_smiles_with_params(smiles: &str, params: &crate::SmilesParseParams) -> Result<Self, SmilesParseError> {
        crate::smiles::mol_from_smiles(smiles, params)
    }

    pub fn from_smiles_with_sanitize(smiles: &str, sanitize: bool) -> Result<Self, SmilesParseError> {
        let params = crate::smiles::SmilesParseParams::with_sanitize(sanitize);
        crate::smiles::mol_from_smiles(smiles, &params)
    }

    #[cfg(feature = "io")]
    pub fn from_mol_block(_block: &str) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::sdf::read_sdf_from_str_with_params(_block, crate::io::sdf::SdfReadParams::default())?.molecule)
    }

    #[cfg(feature = "io")]
    pub fn from_mol_block_with_params(
        _block: &str,
        params: crate::io::sdf::SdfReadParams,
    ) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::sdf::read_sdf_from_str_with_params(_block, params)?.molecule)
    }

    #[cfg(feature = "io")]
    pub fn from_mol_file(_path: impl AsRef<std::path::Path>) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::molfile::read_mol_file(_path)?.molecule)
    }

    #[cfg(feature = "io")]
    pub fn from_mol_file_with_params(
        _path: impl AsRef<std::path::Path>,
        params: crate::io::sdf::SdfReadParams,
    ) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::molfile::read_mol_file_with_params(_path, params)?.molecule)
    }

    #[cfg(feature = "io")]
    pub fn from_xyz_block(block: &str) -> Result<Self, crate::io::xyz::XyzReadError> {
        crate::io::xyz::read_xyz_from_str(block)
    }

    pub fn to_smiles(&self, isomeric_smiles: bool) -> Result<String, SmilesWriteError> {
        let params = crate::SmilesWriteParams {
            do_isomeric_smiles: isomeric_smiles,
            ..Default::default()
        };
        self.to_smiles_with_params(&params)
    }

    pub fn to_smiles_with_params(&self, params: &crate::SmilesWriteParams) -> Result<String, SmilesWriteError> {
        crate::smiles_write::mol_to_smiles(self, params)
    }

    pub fn dg_bounds_matrix(&self) -> Result<Vec<Vec<f64>>, crate::DgBoundsError> {
        crate::distgeom::dg_bounds_matrix(self)
    }

    #[cfg(feature = "fingerprints")]
    pub fn avalon_fingerprint(
        &self,
        params: &crate::avalon_fingerprint::AvalonFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::avalon_fingerprint::avalon_fingerprint(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn morgan_fingerprint(
        &self,
        params: &crate::MorganFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::morgan_fingerprint(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn morgan_fingerprint_with_output(
        &self,
        params: &crate::MorganFingerprintParams,
    ) -> Result<crate::MorganFingerprintOutput, crate::FingerprintError> {
        crate::fingerprint::morgan_fingerprint_with_output(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn atom_pair_fingerprint(
        &self,
        params: &crate::AtomPairFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::atom_pair_fingerprint(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn atom_pair_fingerprint_with_output(
        &self,
        params: &crate::AtomPairFingerprintParams,
    ) -> Result<crate::AtomPairFingerprintOutput, crate::FingerprintError> {
        crate::fingerprint::atom_pair_fingerprint_with_output(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn pattern_fingerprint(
        &self,
        params: &crate::PatternFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::pattern_fingerprint(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn topological_fingerprint(
        &self,
        params: &crate::fingerprint::TopologicalFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::topological_fingerprint(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn topological_fingerprint_with_output(
        &self,
        params: &crate::fingerprint::TopologicalFingerprintParams,
        request: crate::fingerprint::TopologicalFingerprintOutputRequest,
    ) -> Result<crate::fingerprint::TopologicalFingerprintResult, crate::FingerprintError> {
        crate::fingerprint::topological_fingerprint_with_output(self, params, request)
    }

    #[cfg(feature = "fingerprints")]
    pub fn layered_fingerprint(
        &self,
        params: &crate::fingerprint::LayeredFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::layered_fingerprint(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn layered_fingerprint_with_output(
        &self,
        params: &crate::fingerprint::LayeredFingerprintParams,
    ) -> Result<crate::fingerprint::LayeredFingerprintResult, crate::FingerprintError> {
        crate::fingerprint::layered_fingerprint_with_output(self, params)
    }

    #[cfg(feature = "fingerprints")]
    pub fn maccs_fingerprint(
        &self,
        params: &crate::fingerprint::MaccsFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::maccs_fingerprint(self, params)
    }

    #[cfg(feature = "hashing")]
    pub fn hash(&self) -> Result<u64, crate::mol_hash::HashError> {
        crate::mol_hash::mol_hash(self)
    }

    #[cfg(feature = "hashing")]
    pub fn hash_with_ranks(&self, ranks: &[u32]) -> Result<u64, crate::mol_hash::HashError> {
        crate::mol_hash::mol_hash_with_ranks(self, ranks)
    }

    pub fn fragments(&self) -> Result<Vec<Molecule>, crate::fragment::FragmentError> {
        crate::fragment::get_mol_frags(self)
    }

    pub fn largest_fragment(&self) -> Result<Molecule, crate::fragment::FragmentError> {
        crate::fragment::get_largest_fragment(self)
    }

    #[cfg(feature = "hashing")]
    pub fn murcko_scaffold(&self) -> Result<Molecule, crate::mol_hash::HashError> {
        crate::mol_hash::mol_murcko_scaffold(self)
    }

    #[cfg(feature = "hashing")]
    pub fn net_scaffold(&self) -> Result<Molecule, crate::mol_hash::HashError> {
        crate::mol_hash::mol_net_scaffold(self)
    }

    #[cfg(feature = "depict")]
    pub fn to_svg(&self, width: u32, height: u32) -> Result<String, crate::SvgDrawError> {
        crate::draw::mol_to_svg(self, width, height)
    }

    #[cfg(feature = "depict")]
    pub fn to_png(&self, width: u32, height: u32) -> Result<Vec<u8>, crate::SvgDrawError> {
        crate::draw::mol_to_png(self, width, height)
    }

    #[cfg(feature = "io")]
    pub fn to_pdb_block(&self, conf_id: i32, flavor: u32) -> String {
        crate::pdb_writer::mol_to_pdb_block(self, conf_id, flavor)
    }

    #[cfg(feature = "depict")]
    pub(crate) fn prepared_for_drawing_parity(&self) -> Result<crate::draw::PreparedDrawMolecule, crate::SvgDrawError> {
        crate::draw::prepare_mol_for_drawing_parity(self)
    }

    pub fn tetrahedral_stereo(&self) -> Result<Vec<crate::TetrahedralStereo>, crate::StereoError> {
        crate::stereo::tetrahedral_stereo(self)
    }

    pub fn perceive_stereochemistry(&self) -> Result<(), crate::StereoError> {
        crate::stereo::perceive_stereochemistry(self)
    }

    /// Analyze potential stereochemistry without mutating this molecule.
    ///
    /// The result contains the isolated source-defined cleanup state and
    /// ordered typed potential-stereo records.
    pub fn analyze_potential_stereo(
        &self,
        options: crate::PotentialStereoOptions,
    ) -> Result<crate::PotentialStereoAnalysis, crate::PotentialStereoError> {
        crate::potential_stereo::analyze_potential_stereo(self, options)
    }

    /// Enumerate stereoisomers lazily without mutating this molecule.
    ///
    /// Candidate discovery happens when the iterator is created. Configuration
    /// application, uniqueness, optional embedding, and their errors are
    /// deferred until the iterator is consumed.
    #[cfg(feature = "stereoisomers")]
    pub fn stereoisomers(
        &self,
        options: crate::StereoisomerOptions,
    ) -> Result<crate::StereoisomerIterator, crate::EnumerationError> {
        crate::stereo_enumerate::enumerate_stereoisomers(self, options)
    }

    /// Return the source-defined arbitrary-width upper bound for stereoisomer
    /// enumeration under `options`.
    #[cfg(feature = "stereoisomers")]
    pub fn stereoisomer_count(
        &self,
        options: &crate::StereoisomerOptions,
    ) -> Result<num_bigint::BigUint, crate::EnumerationError> {
        crate::stereo_enumerate::stereoisomer_count(self, options)
    }

    #[allow(dead_code)]
    pub(crate) fn topology_block(&self) -> &TopologyBlock {
        &self.topology
    }

    // Operation-body COW accessors are reached through OpParts. Some accessors
    // are intentionally ahead of the first real operation body that uses them.
    #[allow(dead_code)]
    pub(crate) fn topology_block_mut(&mut self) -> &mut TopologyBlock {
        self.clear_computed_property_cache();
        Arc::make_mut(&mut self.topology)
    }

    /// Detach topology for private source-port workspaces whose mutation is
    /// limited to atom/bond properties that do not underlie derived state.
    ///
    /// This is not an operation-body escape hatch. Public molecule operations
    /// must continue to use `OpParts`; callers here must own an isolated value
    /// and must not change topology, chemistry state, or computed properties.
    pub(crate) fn topology_properties_mut_for_private_workspace(&mut self) -> &mut TopologyBlock {
        Arc::make_mut(&mut self.topology)
    }

    pub(crate) fn replace_topology_block(&mut self, topology: TopologyBlock) {
        self.clear_computed_property_cache();
        self.topology = Arc::new(topology);
    }

    /// Commit topology owned by the operation runtime.
    ///
    /// Operation bodies reproduce the source operation's computed-property
    /// policy explicitly through `OpParts::clear_computed_properties()`. The
    /// general mutation accessors above and below remain conservatively
    /// invalidating for non-operation callers.
    pub(crate) fn replace_topology_block_from_operation(&mut self, topology: TopologyBlock) {
        self.topology = Arc::new(topology);
    }

    pub(crate) fn take_topology_block_or_clone(&mut self) -> TopologyBlock {
        self.clear_computed_property_cache();
        self.take_topology_block_or_clone_from_operation()
    }

    /// Move topology into the operation runtime without imposing a cache
    /// transition that belongs to the source operation body.
    pub(crate) fn take_topology_block_or_clone_from_operation(&mut self) -> TopologyBlock {
        let block = std::mem::replace(&mut self.topology, Arc::new(TopologyBlock::default()));
        match Arc::try_unwrap(block) {
            Ok(block) => block,
            Err(block) => (*block).clone(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn coordinate_block(&self) -> &CoordinateBlock {
        &self.coordinates
    }

    #[allow(dead_code)]
    pub(crate) fn coordinate_block_mut(&mut self) -> &mut CoordinateBlock {
        self.clear_coordinate_dependent_computed_property_cache();
        Arc::make_mut(&mut self.coordinates)
    }

    pub(crate) fn replace_coordinate_block(&mut self, coordinates: CoordinateBlock) {
        self.clear_coordinate_dependent_computed_property_cache();
        self.coordinates = Arc::new(coordinates);
    }

    pub(crate) fn take_coordinate_block_or_clone(&mut self) -> CoordinateBlock {
        self.clear_coordinate_dependent_computed_property_cache();
        let block = std::mem::replace(&mut self.coordinates, Arc::new(CoordinateBlock::default()));
        match Arc::try_unwrap(block) {
            Ok(block) => block,
            Err(block) => (*block).clone(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn derived_cache(&self) -> &DerivedCacheBlock {
        &self.derived_cache
    }

    #[allow(dead_code)]
    pub(crate) fn derived_cache_mut(&mut self) -> &mut DerivedCacheBlock {
        Arc::make_mut(&mut self.derived_cache)
    }

    /// Restore a validated cache block while constructing a molecule at an IO
    /// boundary. Runtime operations must continue to use `OpParts` cache APIs.
    pub(crate) fn with_deserialized_derived_cache(
        mut self,
        derived_cache: DerivedCacheBlock,
    ) -> Result<Self, InvariantError> {
        self.derived_cache = Arc::new(derived_cache);
        enforce_molecule_invariants(&self)?;
        Ok(self)
    }

    pub(crate) fn crippen_descriptor_cache(&self) -> Option<[f64; 2]> {
        self.computed_properties.crippen()
    }

    pub(crate) fn set_crippen_descriptor_cache(&self, values: [f64; 2]) {
        self.computed_properties.set_crippen(values);
    }

    pub(crate) fn potential_stereo_cache(&self) -> Option<Arc<[crate::StereoInfo]>> {
        self.computed_properties.potential_stereo()
    }

    pub(crate) fn set_potential_stereo_cache(&self, values: Arc<[crate::StereoInfo]>) {
        self.computed_properties.set_potential_stereo(values);
    }

    pub(crate) fn crippen_atom_contribution_cache(&self) -> Option<CrippenAtomContributionCache> {
        self.computed_properties.crippen_atom_contributions()
    }

    pub(crate) fn set_crippen_atom_contribution_cache(&self, values: CrippenAtomContributionCache) {
        self.computed_properties.set_crippen_atom_contributions(values);
    }

    pub(crate) fn connectivity_hk_deltas_cache(&self) -> Option<Arc<[f64]>> {
        self.computed_properties.connectivity_hk_deltas()
    }

    pub(crate) fn set_connectivity_hk_deltas_cache(&self, values: Arc<[f64]>) {
        self.computed_properties.set_connectivity_hk_deltas(values);
    }

    pub(crate) fn connectivity_n_vals_cache(&self) -> Option<Arc<[f64]>> {
        self.computed_properties.connectivity_n_vals()
    }

    pub(crate) fn set_connectivity_n_vals_cache(&self, values: Arc<[f64]>) {
        self.computed_properties.set_connectivity_n_vals(values);
    }

    pub(crate) fn labute_descriptor_cache(&self) -> Option<LabuteDescriptorCache> {
        self.computed_properties.labute()
    }

    pub(crate) fn set_labute_descriptor_cache(&self, values: LabuteDescriptorCache) {
        self.computed_properties.set_labute(values);
    }

    pub(crate) fn topological_distance_matrix_cache_or_init(
        &self,
        initialize: impl FnOnce() -> Vec<f64>,
    ) -> Arc<[f64]> {
        self.computed_properties.topological_distance_matrix_or_init(initialize)
    }

    pub(crate) fn distance_matrix_3d_cache_or_init(
        &self,
        conformer_id: usize,
        initialize: impl FnOnce() -> Vec<f64>,
    ) -> Arc<[f64]> {
        self.computed_properties
            .distance_matrix_3d_or_init(conformer_id, initialize)
    }

    pub(crate) fn clear_computed_property_cache(&self) {
        self.computed_properties.clear();
    }

    fn clear_coordinate_dependent_computed_property_cache(&self) {
        self.computed_properties.clear_coordinate_dependent();
    }

    #[allow(dead_code)]
    pub(crate) fn properties_mut(&mut self) -> &mut MoleculeProperties {
        Arc::make_mut(&mut self.properties)
    }

    pub(crate) fn replace_properties(&mut self, properties: MoleculeProperties) {
        self.properties = Arc::new(properties);
    }

    pub(crate) fn take_properties_or_clone(&mut self) -> MoleculeProperties {
        let block = std::mem::replace(&mut self.properties, Arc::new(MoleculeProperties::default()));
        match Arc::try_unwrap(block) {
            Ok(block) => block,
            Err(block) => (*block).clone(),
        }
    }

    pub(crate) fn capabilities_block(&self) -> MoleculeCapabilities {
        *self.capabilities
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::{LabuteDescriptorCache, Molecule};
    use crate::avalon_fingerprint::{AvalonFingerprintParams, avalon_fingerprint};
    use crate::fingerprint::{MaccsFingerprintParams, TopologicalFingerprintParams};
    use crate::{
        AtomSpec, BondOrder, BondSpec, Conformer3D, Element, MoleculeBuilder, assign_stereochemistry, fragment,
        mol_hash, pdb_writer, perceive_stereochemistry,
    };

    #[test]
    fn molecule_fragment_helpers_match_module_functions() {
        let mol = Molecule::from_smiles("CC.C").expect("failed to parse fragments test molecule");

        let via_method = mol.fragments().expect("method fragments");
        let via_module = fragment::get_mol_frags(&mol).expect("module fragments");

        assert_eq!(via_method, via_module);
    }

    #[test]
    fn molecule_fragments_do_not_materialize_absent_isotope_or_atom_map_as_zero() {
        let mol = Molecule::from_smiles("CC.O").expect("failed to parse fragments test molecule");
        let smiles: Vec<_> = mol
            .fragments()
            .expect("method fragments")
            .into_iter()
            .map(|fragment| fragment.to_smiles(true).expect("fragment smiles"))
            .collect();

        assert_eq!(smiles, vec!["CC", "O"]);
    }

    #[test]
    fn molecule_tetrahedral_stereo_reports_smiles_chiral_center() {
        let mol = Molecule::from_smiles("F[C@H](Cl)Br").expect("failed to parse chiral molecule");
        let stereo = mol.tetrahedral_stereo().expect("tetrahedral stereo");

        assert_eq!(stereo.len(), 1);
        assert_eq!(stereo[0].center.index(), 1);
    }

    #[test]
    fn molecule_largest_fragment_helper_matches_module_function() {
        let mol = Molecule::from_smiles("CC.C").expect("failed to parse fragments test molecule");

        let via_method = mol.largest_fragment().expect("method largest fragment");
        let via_module = fragment::get_largest_fragment(&mol).expect("module largest fragment");

        assert_eq!(via_method, via_module);
    }

    #[test]
    fn molecule_scaffold_helpers_match_module_functions() {
        let mol = Molecule::from_smiles("Cc1ccccc1").expect("failed to parse scaffold molecule");

        let murcko_method = mol.murcko_scaffold().expect("method murcko scaffold");
        let murcko_module = crate::mol_hash::mol_murcko_scaffold(&mol).expect("module murcko scaffold");
        let net_method = mol.net_scaffold().expect("method net scaffold");
        let net_module = crate::mol_hash::mol_net_scaffold(&mol).expect("module net scaffold");

        assert_eq!(murcko_method, murcko_module);
        assert_eq!(net_method, net_module);
    }

    #[test]
    fn molecule_fingerprint_helpers_match_module_functions() {
        let mol = Molecule::from_smiles("CCO").expect("failed to parse fingerprint molecule");

        let avalon_params = AvalonFingerprintParams::default();
        let topological_params = TopologicalFingerprintParams::default();
        let maccs_params = MaccsFingerprintParams::default();

        assert_eq!(
            mol.avalon_fingerprint(&avalon_params)
                .expect("method Avalon fingerprint"),
            avalon_fingerprint(&mol, &avalon_params).expect("module Avalon fingerprint")
        );
        assert_eq!(
            mol.topological_fingerprint(&topological_params)
                .expect("method topological fingerprint"),
            crate::fingerprint::topological_fingerprint(&mol, &topological_params)
                .expect("module topological fingerprint")
        );
        assert_eq!(
            mol.maccs_fingerprint(&maccs_params).expect("method MACCS fingerprint"),
            crate::fingerprint::maccs_fingerprint(&mol, &maccs_params).expect("module MACCS fingerprint")
        );
    }

    #[test]
    fn molecule_hash_helpers_match_module_functions() {
        let mol = Molecule::from_smiles("Cl[C@H](Br)I").expect("failed to parse hash molecule");
        let ranks = crate::stereo::assign_atom_cip_ranks(&mol).expect("cip ranks");

        assert_eq!(
            mol.hash().expect("method hash"),
            mol_hash::mol_hash(&mol).expect("module hash")
        );
        assert_eq!(
            mol.hash_with_ranks(&ranks).expect("method hash with ranks"),
            mol_hash::mol_hash_with_ranks(&mol, &ranks).expect("module hash with ranks")
        );
    }

    #[test]
    fn molecule_pdb_helper_matches_module_function() {
        let mol = Molecule::from_smiles("CO").expect("failed to parse pdb molecule");

        assert_eq!(mol.to_pdb_block(-1, 0), pdb_writer::mol_to_pdb_block(&mol, -1, 0));
    }

    #[test]
    fn molecule_perceive_stereochemistry_matches_module_function() {
        let mol = Molecule::from_smiles("Cl[C@H](Br)I").expect("failed to parse stereo molecule");

        assert_eq!(mol.perceive_stereochemistry(), perceive_stereochemistry(&mol));
        #[allow(deprecated)]
        {
            assert_eq!(assign_stereochemistry(&mol), perceive_stereochemistry(&mol));
        }
    }

    #[test]
    fn molecule_pdb_helper_uses_conformer_index_selection() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_conformer(Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], true))
            .expect("failed to add conformer");
        let mol = builder.build().expect("failed to build conformer molecule");

        let block = mol.to_pdb_block(0, 0);
        assert_eq!(
            block,
            "HETATM    1  O1  UNL     1      +1.000  +2.000  +3.000  1.00  0.00           O  \nEND\n"
        );
    }

    #[test]
    fn molecule_perceive_stereochemistry_is_read_only() {
        let mol = Molecule::from_smiles("F/C=C/F").expect("failed to parse stereo molecule");
        let before = mol.to_smiles(true).expect("smiles before");
        let _ = mol.perceive_stereochemistry().expect("perceive stereochemistry");
        let after = mol.to_smiles(true).expect("smiles after");
        assert_eq!(before, after);
    }

    #[test]
    fn molecule_topological_fingerprint_helper_matches_module_function_with_custom_params() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        let o = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .expect("failed to add c-c bond");
        builder
            .add_bond(BondSpec::new(c2, o, BondOrder::Single))
            .expect("failed to add c-o bond");
        let mol = builder.build().expect("failed to build fingerprint molecule");
        let params = TopologicalFingerprintParams {
            min_path: 1,
            max_path: 2,
            fp_size: 512,
            num_bits_per_feature: 3,
            use_hs: true,
            target_density: 0.0,
            min_size: 128,
            branched_paths: true,
            use_bond_order: false,
            atom_invariants: None,
            from_atoms: Some(vec![0]),
        };
        assert_eq!(
            mol.topological_fingerprint(&params).expect("Molecule helper result"),
            crate::fingerprint::topological_fingerprint(&mol, &params).expect("module function result")
        );
    }

    #[test]
    fn descriptor_computed_cache_lifecycle_is_typed_independent_and_topology_invalidated() {
        let mut molecule = Molecule::from_smiles("CCO").expect("descriptor cache molecule");
        assert!(molecule.connectivity_hk_deltas_cache().is_none());
        assert!(molecule.connectivity_n_vals_cache().is_none());
        assert!(molecule.labute_descriptor_cache().is_none());

        let hk_deltas = Arc::<[f64]>::from([1.0, 2.0, 3.0]);
        let n_vals = Arc::<[f64]>::from([4.0, 5.0, 6.0]);
        let labute = LabuteDescriptorCache {
            atom_contributions: Arc::<[f64]>::from([7.0, 8.0, 9.0]),
            hydrogen_contribution: 10.0,
            asa: 11.0,
        };
        molecule.set_connectivity_hk_deltas_cache(Arc::clone(&hk_deltas));
        molecule.set_labute_descriptor_cache(labute.clone());
        molecule.set_connectivity_n_vals_cache(Arc::clone(&n_vals));

        assert_eq!(molecule.connectivity_hk_deltas_cache(), Some(hk_deltas));
        assert_eq!(molecule.connectivity_n_vals_cache(), Some(n_vals));
        assert_eq!(molecule.labute_descriptor_cache(), Some(labute));

        let clone = molecule.clone();
        molecule.set_connectivity_hk_deltas_cache(Arc::<[f64]>::from([12.0, 13.0, 14.0]));
        molecule.set_connectivity_n_vals_cache(Arc::<[f64]>::from([15.0, 16.0, 17.0]));
        molecule.set_labute_descriptor_cache(LabuteDescriptorCache {
            atom_contributions: Arc::<[f64]>::from([18.0, 19.0, 20.0]),
            hydrogen_contribution: 21.0,
            asa: 22.0,
        });
        assert_eq!(
            clone.connectivity_hk_deltas_cache().as_deref(),
            Some(&[1.0, 2.0, 3.0][..])
        );
        assert_eq!(clone.connectivity_n_vals_cache().as_deref(), Some(&[4.0, 5.0, 6.0][..]));
        assert_eq!(clone.labute_descriptor_cache().unwrap().asa, 11.0);

        let coordinates = molecule.coordinate_block().clone();
        molecule.replace_coordinate_block(coordinates);
        assert!(molecule.connectivity_hk_deltas_cache().is_some());
        assert!(molecule.connectivity_n_vals_cache().is_some());
        assert!(molecule.labute_descriptor_cache().is_some());

        let topology = molecule.topology_block().clone();
        molecule.replace_topology_block(topology);
        assert!(molecule.connectivity_hk_deltas_cache().is_none());
        assert!(molecule.connectivity_n_vals_cache().is_none());
        assert!(molecule.labute_descriptor_cache().is_none());
    }

    #[test]
    fn descriptor_computed_cache_supports_parallel_reads_without_cross_entry_aliasing() {
        let molecule = Molecule::from_smiles("CCO").expect("parallel descriptor cache molecule");
        molecule.set_connectivity_hk_deltas_cache(Arc::<[f64]>::from([1.0, 2.0, 3.0]));
        molecule.set_connectivity_n_vals_cache(Arc::<[f64]>::from([4.0, 5.0, 6.0]));
        molecule.set_labute_descriptor_cache(LabuteDescriptorCache {
            atom_contributions: Arc::<[f64]>::from([7.0, 8.0, 9.0]),
            hydrogen_contribution: 10.0,
            asa: 11.0,
        });
        let molecule = Arc::new(molecule);

        let readers = (0..16)
            .map(|_| {
                let molecule = Arc::clone(&molecule);
                std::thread::spawn(move || {
                    assert_eq!(
                        molecule.connectivity_hk_deltas_cache().as_deref(),
                        Some(&[1.0, 2.0, 3.0][..])
                    );
                    assert_eq!(
                        molecule.connectivity_n_vals_cache().as_deref(),
                        Some(&[4.0, 5.0, 6.0][..])
                    );
                    let labute = molecule.labute_descriptor_cache().unwrap();
                    assert_eq!(labute.atom_contributions.as_ref(), &[7.0, 8.0, 9.0]);
                    assert_eq!(labute.hydrogen_contribution, 10.0);
                    assert_eq!(labute.asa, 11.0);
                })
            })
            .collect::<Vec<_>>();

        for reader in readers {
            reader.join().expect("parallel descriptor cache reader");
        }
    }
}
