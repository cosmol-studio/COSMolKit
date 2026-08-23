use std::{
    collections::BTreeMap,
    sync::{Arc, RwLock},
};

use crate::{
    Atom, AtomId, Bond, MoleculeBuilder, error::InvariantError,
    invariants::enforce_molecule_invariants, sgroup::SubstanceGroup, stereo::StereoGroup,
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
pub enum CoordinateDimension {
    TwoD,
    ThreeD,
}

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

#[derive(Debug, Clone, PartialEq)]
pub struct Conformer2D {
    id: usize,
    coords: Vec<[f64; 2]>,
    props: BTreeMap<String, String>,
}

impl Conformer2D {
    #[must_use]
    pub fn new(id: usize, coords: Vec<[f64; 2]>) -> Self {
        Self {
            id,
            coords,
            props: BTreeMap::new(),
        }
    }

    #[must_use]
    pub const fn id(&self) -> usize {
        self.id
    }

    #[must_use]
    pub fn coordinates(&self) -> &[[f64; 2]] {
        &self.coords
    }

    pub fn coordinates_mut(&mut self) -> &mut [[f64; 2]] {
        &mut self.coords
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub(crate) fn with_id(mut self, id: usize) -> Self {
        self.id = id;
        self
    }

    #[allow(dead_code)]
    pub(crate) fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize], id: usize) -> Self {
        let coords = kept_old_indices
            .iter()
            .filter_map(|old_idx| self.coords.get(*old_idx).copied())
            .collect();
        Self {
            id,
            coords,
            props: self.props.clone(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn push_coord(&mut self, coord: [f64; 2]) {
        self.coords.push(coord);
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ConformerStore {
    pub conformers_2d: Vec<Conformer2D>,
    pub conformers_3d: Vec<Conformer3D>,
    pub source_coordinate_dim: Option<CoordinateDimension>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct PropertyStore {
    pub name: Option<String>,
    pub sdf_data_fields: Vec<(String, String)>,
    pub sdf_property_lists: Vec<SdfPropertyList>,
    pub props: std::collections::BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Conformer3D {
    id: usize,
    coords: Vec<[f64; 3]>,
    is_3d: bool,
    props: BTreeMap<String, String>,
}

impl Conformer3D {
    #[must_use]
    pub fn new(id: usize, coords: Vec<[f64; 3]>, is_3d: bool) -> Self {
        Self {
            id,
            coords,
            is_3d,
            props: BTreeMap::new(),
        }
    }

    #[must_use]
    pub const fn id(&self) -> usize {
        self.id
    }

    #[must_use]
    pub fn coordinates(&self) -> &[[f64; 3]] {
        &self.coords
    }

    /// Mutable access to coordinates (pub(crate) for conformer transforms).
    pub fn coordinates_mut(&mut self) -> &mut [[f64; 3]] {
        &mut self.coords
    }

    #[must_use]
    pub const fn is_3d(&self) -> bool {
        self.is_3d
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub(crate) fn with_id(mut self, id: usize) -> Self {
        self.id = id;
        self
    }

    #[allow(dead_code)]
    pub(crate) fn remapped_to_kept_atoms(&self, kept_old_indices: &[usize], id: usize) -> Self {
        let coords = kept_old_indices
            .iter()
            .filter_map(|old_idx| self.coords.get(*old_idx).copied())
            .collect();
        Self {
            id,
            coords,
            is_3d: self.is_3d,
            props: self.props.clone(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn push_coord(&mut self, coord: [f64; 3]) {
        self.coords.push(coord);
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MoleculeProperties {
    name: Option<String>,
    sdf_data_fields: Vec<(String, String)>,
    sdf_property_lists: Vec<SdfPropertyList>,
    props: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdfPropertyListTarget {
    Atom,
    Bond,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfPropertyList {
    target: SdfPropertyListTarget,
    name: String,
    values: Vec<Option<String>>,
}

impl SdfPropertyList {
    #[must_use]
    pub fn new(
        target: SdfPropertyListTarget,
        name: impl Into<String>,
        values: Vec<Option<String>>,
    ) -> Self {
        Self {
            target,
            name: name.into(),
            values,
        }
    }

    #[must_use]
    pub const fn target(&self) -> SdfPropertyListTarget {
        self.target
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn values(&self) -> &[Option<String>] {
        &self.values
    }
}

impl MoleculeProperties {
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    #[must_use]
    pub fn sdf_data_fields(&self) -> &[(String, String)] {
        &self.sdf_data_fields
    }

    #[must_use]
    pub fn sdf_property_lists(&self) -> &[SdfPropertyList] {
        &self.sdf_property_lists
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_sdf_data_field(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.sdf_data_fields.push((key.into(), value.into()));
        self
    }

    #[must_use]
    pub fn with_sdf_property_list(mut self, property_list: SdfPropertyList) -> Self {
        self.sdf_property_lists.push(property_list);
        self
    }

    #[allow(dead_code)]
    pub(crate) fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.props.insert(key.into(), value.into());
    }

    #[allow(dead_code)]
    pub(crate) fn clear_prop(&mut self, key: &str) {
        self.props.remove(key);
    }

    pub(crate) fn remap_topology(&mut self, mapping: &TopologyMapping) {
        self.sdf_property_lists = self
            .sdf_property_lists
            .iter()
            .map(|property_list| property_list.remapped_topology(mapping))
            .collect();
    }
}

impl SdfPropertyList {
    fn remapped_topology(&self, mapping: &TopologyMapping) -> Self {
        let values = match self.target {
            SdfPropertyListTarget::Atom => mapping
                .atoms
                .new_to_old
                .iter()
                .map(|old_row| {
                    old_row.and_then(|row| self.values.get(row.index()).cloned().flatten())
                })
                .collect(),
            SdfPropertyListTarget::Bond => mapping
                .bonds
                .new_to_old
                .iter()
                .map(|old_row| {
                    old_row.and_then(|row| self.values.get(row.index()).cloned().flatten())
                })
                .collect(),
        };
        Self {
            target: self.target,
            name: self.name.clone(),
            values,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TopologyBlock {
    pub(crate) atoms: Vec<Atom>,
    pub(crate) bonds: Vec<Bond>,
    pub(crate) adjacency: crate::AdjacencyList,
    pub(crate) substance_groups: Vec<SubstanceGroup>,
    pub(crate) stereo_groups: Vec<StereoGroup>,
}

impl Default for TopologyBlock {
    fn default() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: crate::AdjacencyList::from_topology(0, &[]),
            substance_groups: Vec::new(),
            stereo_groups: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct CoordinateBlock {
    /// Zero or more 2D conformers.
    ///
    /// Coordinates are stored in the same atom-index order as `TopologyBlock`.
    /// Any operation changing atom indices must remap or drop this block through
    /// a topology report. Do not mutate this block directly from operation code.
    pub(crate) conformers_2d: Vec<Conformer2D>,
    pub(crate) conformers_3d: Vec<Conformer3D>,
    pub(crate) source_coordinate_dim: Option<CoordinateDimension>,
}

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
    topological_distance_matrix: RwLock<Option<Arc<[f64]>>>,
    distance_matrices_3d: RwLock<BTreeMap<usize, Arc<[f64]>>>,
}

impl Clone for ComputedPropertyCache {
    fn clone(&self) -> Self {
        Self {
            crippen: RwLock::new(self.crippen()),
            topological_distance_matrix: RwLock::new(self.topological_distance_matrix()),
            distance_matrices_3d: RwLock::new(self.distance_matrices_3d()),
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
        *self
            .crippen
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    fn set_crippen(&self, values: [f64; 2]) {
        *self
            .crippen
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(values);
    }

    fn topological_distance_matrix(&self) -> Option<Arc<[f64]>> {
        self.topological_distance_matrix
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn topological_distance_matrix_or_init(
        &self,
        initialize: impl FnOnce() -> Vec<f64>,
    ) -> Arc<[f64]> {
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

    fn distance_matrix_3d_or_init(
        &self,
        conformer_id: usize,
        initialize: impl FnOnce() -> Vec<f64>,
    ) -> Arc<[f64]> {
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
        *self
            .crippen
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        *self
            .topological_distance_matrix
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = None;
        self.clear_coordinate_dependent();
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomMapping {
    pub(crate) old_to_new: Vec<Option<AtomId>>,
    pub(crate) new_to_old: Vec<Option<AtomId>>,
}

impl AtomMapping {
    #[must_use]
    pub fn old_to_new(&self) -> &[Option<AtomId>] {
        &self.old_to_new
    }

    #[must_use]
    pub fn new_to_old(&self) -> &[Option<AtomId>] {
        &self.new_to_old
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondMapping {
    pub(crate) old_to_new: Vec<Option<crate::BondId>>,
    pub(crate) new_to_old: Vec<Option<crate::BondId>>,
}

impl BondMapping {
    #[must_use]
    pub fn old_to_new(&self) -> &[Option<crate::BondId>] {
        &self.old_to_new
    }

    #[must_use]
    pub fn new_to_old(&self) -> &[Option<crate::BondId>] {
        &self.new_to_old
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyMapping {
    pub(crate) atoms: AtomMapping,
    pub(crate) bonds: BondMapping,
}

impl TopologyMapping {
    pub(crate) fn identity(atom_count: usize, bond_count: usize) -> Self {
        Self {
            atoms: AtomMapping {
                old_to_new: (0..atom_count).map(|idx| Some(AtomId::new(idx))).collect(),
                new_to_old: (0..atom_count).map(|idx| Some(AtomId::new(idx))).collect(),
            },
            bonds: BondMapping {
                old_to_new: (0..bond_count)
                    .map(|idx| Some(crate::BondId::new(idx)))
                    .collect(),
                new_to_old: (0..bond_count)
                    .map(|idx| Some(crate::BondId::new(idx)))
                    .collect(),
            },
        }
    }

    pub(crate) fn with_appended(
        old_atom_count: usize,
        old_bond_count: usize,
        added_atom_count: usize,
        added_bond_count: usize,
    ) -> Self {
        let mut atom_new_to_old = (0..old_atom_count)
            .map(|idx| Some(AtomId::new(idx)))
            .collect::<Vec<_>>();
        atom_new_to_old.extend((0..added_atom_count).map(|_| None));

        let mut bond_new_to_old = (0..old_bond_count)
            .map(|idx| Some(crate::BondId::new(idx)))
            .collect::<Vec<_>>();
        bond_new_to_old.extend((0..added_bond_count).map(|_| None));

        Self {
            atoms: AtomMapping {
                old_to_new: (0..old_atom_count)
                    .map(|idx| Some(AtomId::new(idx)))
                    .collect(),
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: (0..old_bond_count)
                    .map(|idx| Some(crate::BondId::new(idx)))
                    .collect(),
                new_to_old: bond_new_to_old,
            },
        }
    }

    #[must_use]
    pub fn atoms(&self) -> &AtomMapping {
        &self.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &BondMapping {
        &self.bonds
    }
}

impl TopologyBlock {
    #[allow(dead_code)]
    pub(crate) fn remove_atoms_with_mapping(
        &mut self,
        atoms_to_remove: &[AtomId],
    ) -> TopologyMapping {
        let mut remove_atom = vec![false; self.atoms.len()];
        for atom in atoms_to_remove {
            if let Some(slot) = remove_atom.get_mut(atom.index()) {
                *slot = true;
            }
        }

        let mut atom_old_to_new = vec![None; self.atoms.len()];
        let mut atom_new_to_old = Vec::new();
        let mut atoms = Vec::with_capacity(self.atoms.len().saturating_sub(atoms_to_remove.len()));
        for atom in &self.atoms {
            if remove_atom[atom.id().index()] {
                continue;
            }
            let new_id = AtomId::new(atoms.len());
            atom_old_to_new[atom.id().index()] = Some(new_id);
            atom_new_to_old.push(Some(atom.id()));
            atoms.push(atom.clone().with_id(new_id));
        }

        let mut bond_old_to_new = vec![None; self.bonds.len()];
        let mut bond_new_to_old = Vec::new();
        let mut bonds = Vec::new();
        for bond in &self.bonds {
            let Some(begin) = atom_old_to_new
                .get(bond.begin().index())
                .and_then(|mapped| *mapped)
            else {
                continue;
            };
            let Some(end) = atom_old_to_new
                .get(bond.end().index())
                .and_then(|mapped| *mapped)
            else {
                continue;
            };

            let stereo_atoms = bond.stereo_atoms().and_then(|[left, right]| {
                Some([
                    atom_old_to_new
                        .get(left.index())
                        .and_then(|mapped| *mapped)?,
                    atom_old_to_new
                        .get(right.index())
                        .and_then(|mapped| *mapped)?,
                ])
            });

            let new_id = crate::BondId::new(bonds.len());
            bond_old_to_new[bond.id().index()] = Some(new_id);
            bond_new_to_old.push(Some(bond.id()));
            bonds.push(bond.clone().remapped(new_id, begin, end, stereo_atoms));
        }

        let mut sgroup_survives = self
            .substance_groups
            .iter()
            .map(|sgroup| sgroup.can_remap_without_parent(&atom_old_to_new, &bond_old_to_new))
            .collect::<Vec<_>>();
        loop {
            let mut changed = false;
            for (idx, sgroup) in self.substance_groups.iter().enumerate() {
                if !sgroup_survives[idx] {
                    continue;
                }
                if let Some(parent) = sgroup.parent()
                    && !sgroup_survives
                        .get(parent.index())
                        .copied()
                        .unwrap_or(false)
                {
                    sgroup_survives[idx] = false;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }

        let mut sgroup_map = vec![None; self.substance_groups.len()];
        let mut next_sgroup_index = 0usize;
        for (old_index, survives) in sgroup_survives.iter().copied().enumerate() {
            if survives {
                sgroup_map[old_index] = Some(crate::SubstanceGroupId::new(next_sgroup_index));
                next_sgroup_index += 1;
            }
        }
        self.substance_groups = self
            .substance_groups
            .iter()
            .enumerate()
            .filter_map(|(old_index, sgroup)| {
                sgroup_map[old_index].and_then(|new_id| {
                    sgroup.remapped(new_id, &atom_old_to_new, &bond_old_to_new, &sgroup_map)
                })
            })
            .collect();

        self.stereo_groups = self
            .stereo_groups
            .iter()
            .filter_map(|group| group.remapped(&atom_old_to_new, &bond_old_to_new))
            .collect();

        self.atoms = atoms;
        self.bonds = bonds;
        self.adjacency = crate::AdjacencyList::from_topology(self.atoms.len(), &self.bonds);

        TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_old_to_new,
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: bond_old_to_new,
                new_to_old: bond_new_to_old,
            },
        }
    }
}

impl CoordinateBlock {
    #[allow(dead_code)]
    pub(crate) fn remap_topology(&mut self, mapping: &TopologyMapping) {
        let kept_old_indices: Vec<_> = mapping
            .atoms
            .new_to_old
            .iter()
            .filter_map(|old_atom| old_atom.map(AtomId::index))
            .collect();

        self.conformers_2d = self
            .conformers_2d
            .iter()
            .enumerate()
            .map(|(id, conformer)| conformer.remapped_to_kept_atoms(&kept_old_indices, id))
            .collect();

        self.conformers_3d = self
            .conformers_3d
            .iter()
            .enumerate()
            .map(|(id, conformer)| conformer.remapped_to_kept_atoms(&kept_old_indices, id))
            .collect();
    }
}

/// Immutable molecule value.
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
        topology.adjacency =
            crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
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
        topology.adjacency =
            crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
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
        topology.adjacency =
            crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
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
        self.coordinates
            .conformers_2d
            .first()
            .map(Conformer2D::coordinates)
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

    pub fn from_smiles_with_sanitize(
        smiles: &str,
        sanitize: bool,
    ) -> Result<Self, SmilesParseError> {
        let params = crate::smiles::SmilesParseParams::with_sanitize(sanitize);
        crate::smiles::mol_from_smiles(smiles, &params)
    }

    pub fn from_mol_block(_block: &str) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::sdf::read_sdf_from_str_with_params(
            _block,
            crate::io::sdf::SdfReadParams::default(),
        )?
        .molecule)
    }

    pub fn from_mol_block_with_params(
        _block: &str,
        params: crate::io::sdf::SdfReadParams,
    ) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::sdf::read_sdf_from_str_with_params(_block, params)?.molecule)
    }

    pub fn from_mol_file(
        _path: impl AsRef<std::path::Path>,
    ) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::molfile::read_mol_file(_path)?.molecule)
    }

    pub fn from_mol_file_with_params(
        _path: impl AsRef<std::path::Path>,
        params: crate::io::sdf::SdfReadParams,
    ) -> Result<Self, crate::io::sdf::SdfReadError> {
        Ok(crate::io::molfile::read_mol_file_with_params(_path, params)?.molecule)
    }

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

    pub fn to_smiles_with_params(
        &self,
        params: &crate::SmilesWriteParams,
    ) -> Result<String, SmilesWriteError> {
        crate::smiles_write::mol_to_smiles(self, params)
    }

    pub fn dg_bounds_matrix(&self) -> Result<Vec<Vec<f64>>, crate::DgBoundsError> {
        crate::distgeom::dg_bounds_matrix(self)
    }

    pub fn avalon_fingerprint(
        &self,
        params: &crate::avalon_fingerprint::AvalonFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::avalon_fingerprint::avalon_fingerprint(self, params)
    }

    pub fn morgan_fingerprint(
        &self,
        params: &crate::MorganFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::morgan_fingerprint(self, params)
    }

    pub fn morgan_fingerprint_with_output(
        &self,
        params: &crate::MorganFingerprintParams,
    ) -> Result<crate::MorganFingerprintOutput, crate::FingerprintError> {
        crate::fingerprint::morgan_fingerprint_with_output(self, params)
    }

    pub fn atom_pair_fingerprint(
        &self,
        params: &crate::AtomPairFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::atom_pair_fingerprint(self, params)
    }

    pub fn atom_pair_fingerprint_with_output(
        &self,
        params: &crate::AtomPairFingerprintParams,
    ) -> Result<crate::AtomPairFingerprintOutput, crate::FingerprintError> {
        crate::fingerprint::atom_pair_fingerprint_with_output(self, params)
    }

    pub fn topological_fingerprint(
        &self,
        params: &crate::fingerprint::TopologicalFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::topological_fingerprint(self, params)
    }

    pub fn topological_fingerprint_with_output(
        &self,
        params: &crate::fingerprint::TopologicalFingerprintParams,
        request: crate::fingerprint::TopologicalFingerprintOutputRequest,
    ) -> Result<crate::fingerprint::TopologicalFingerprintResult, crate::FingerprintError> {
        crate::fingerprint::topological_fingerprint_with_output(self, params, request)
    }

    pub fn maccs_fingerprint(
        &self,
        params: &crate::fingerprint::MaccsFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        crate::fingerprint::maccs_fingerprint(self, params)
    }

    pub fn hash(&self) -> Result<u64, crate::mol_hash::HashError> {
        crate::mol_hash::mol_hash(self)
    }

    pub fn hash_with_ranks(&self, ranks: &[u32]) -> Result<u64, crate::mol_hash::HashError> {
        crate::mol_hash::mol_hash_with_ranks(self, ranks)
    }

    pub fn fragments(&self) -> Result<Vec<Molecule>, crate::fragment::FragmentError> {
        crate::fragment::get_mol_frags(self)
    }

    pub fn largest_fragment(&self) -> Result<Molecule, crate::fragment::FragmentError> {
        crate::fragment::get_largest_fragment(self)
    }

    pub fn murcko_scaffold(&self) -> Result<Molecule, crate::mol_hash::HashError> {
        crate::mol_hash::mol_murcko_scaffold(self)
    }

    pub fn net_scaffold(&self) -> Result<Molecule, crate::mol_hash::HashError> {
        crate::mol_hash::mol_net_scaffold(self)
    }

    pub fn to_svg(&self, width: u32, height: u32) -> Result<String, crate::SvgDrawError> {
        crate::draw::mol_to_svg(self, width, height)
    }

    pub fn to_png(&self, width: u32, height: u32) -> Result<Vec<u8>, crate::SvgDrawError> {
        crate::draw::mol_to_png(self, width, height)
    }

    pub fn to_pdb_block(&self, conf_id: i32, flavor: u32) -> String {
        crate::pdb_writer::mol_to_pdb_block(self, conf_id, flavor)
    }

    pub(crate) fn prepared_for_drawing_parity(
        &self,
    ) -> Result<crate::draw::PreparedDrawMolecule, crate::SvgDrawError> {
        crate::draw::prepare_mol_for_drawing_parity(self)
    }

    pub fn tetrahedral_stereo(&self) -> Result<Vec<crate::TetrahedralStereo>, crate::StereoError> {
        crate::stereo::tetrahedral_stereo(self)
    }

    pub fn perceive_stereochemistry(&self) -> Result<(), crate::StereoError> {
        crate::stereo::perceive_stereochemistry(self)
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

    pub(crate) fn topological_distance_matrix_cache_or_init(
        &self,
        initialize: impl FnOnce() -> Vec<f64>,
    ) -> Arc<[f64]> {
        self.computed_properties
            .topological_distance_matrix_or_init(initialize)
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
        let block = std::mem::replace(
            &mut self.properties,
            Arc::new(MoleculeProperties::default()),
        );
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
    use super::Molecule;
    use crate::avalon_fingerprint::{AvalonFingerprintParams, avalon_fingerprint};
    use crate::fingerprint::{MaccsFingerprintParams, TopologicalFingerprintParams};
    use crate::{
        AtomSpec, BondOrder, BondSpec, Conformer3D, Element, MoleculeBuilder,
        assign_stereochemistry, fragment, mol_hash, pdb_writer, perceive_stereochemistry,
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
        let murcko_module =
            crate::mol_hash::mol_murcko_scaffold(&mol).expect("module murcko scaffold");
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
            mol.maccs_fingerprint(&maccs_params)
                .expect("method MACCS fingerprint"),
            crate::fingerprint::maccs_fingerprint(&mol, &maccs_params)
                .expect("module MACCS fingerprint")
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

        assert_eq!(
            mol.to_pdb_block(-1, 0),
            pdb_writer::mol_to_pdb_block(&mol, -1, 0)
        );
    }

    #[test]
    fn molecule_perceive_stereochemistry_matches_module_function() {
        let mol = Molecule::from_smiles("Cl[C@H](Br)I").expect("failed to parse stereo molecule");

        assert_eq!(
            mol.perceive_stereochemistry(),
            perceive_stereochemistry(&mol)
        );
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
        let _ = mol
            .perceive_stereochemistry()
            .expect("perceive stereochemistry");
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
        let mol = builder
            .build()
            .expect("failed to build fingerprint molecule");
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
            mol.topological_fingerprint(&params)
                .expect("Molecule helper result"),
            crate::fingerprint::topological_fingerprint(&mol, &params)
                .expect("module function result")
        );
    }
}
