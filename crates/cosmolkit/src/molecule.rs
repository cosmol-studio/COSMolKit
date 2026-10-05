//! Runtime-owned live molecule state.

use std::fmt;
use std::sync::Arc;

use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension,
    MoleculeProperties, SdfPropertyListTarget, TopologyBlock, TopologyValidationError,
};

use crate::MoleculeBuilder;
use crate::ops::{DerivedState, OperationError};

/// Private runtime authority for derived-state cache storage.
///
/// RUN-effects extends the contents and transition rules. Keeping the cache in
/// the live state now ensures detached model values cannot acquire cache
/// authority and that state replacement also separates cache ownership.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct DerivedCacheBlock {
    valid: DerivedState,
    #[cfg(any(
        feature = "cap-valence",
        feature = "cap-fingerprints",
        feature = "cap-hydrogens",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-descriptors",
        feature = "cap-forcefields"
    ))]
    valence: Option<cosmolkit_core::ValenceAssignment>,
    #[cfg(any(
        feature = "cap-rings",
        feature = "cap-descriptors",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-hydrogens",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-fingerprints"
    ))]
    rings: Option<cosmolkit_core::RingInfo>,
    #[cfg(any(feature = "cap-rings", feature = "cap-fingerprints"))]
    ring_families: Option<cosmolkit_core::RingInfo>,
}

impl DerivedCacheBlock {
    fn is_empty(&self) -> bool {
        self.valid == DerivedState::NONE
            && {
                #[cfg(any(
                    feature = "cap-valence",
                    feature = "cap-fingerprints",
                    feature = "cap-hydrogens",
                    feature = "cap-smiles",
                    feature = "cap-sanitize",
                    feature = "cap-descriptors",
                    feature = "cap-forcefields"
                ))]
                {
                    self.valence.is_none()
                }
                #[cfg(not(any(
                    feature = "cap-valence",
                    feature = "cap-fingerprints",
                    feature = "cap-hydrogens",
                    feature = "cap-smiles",
                    feature = "cap-sanitize",
                    feature = "cap-descriptors",
                    feature = "cap-forcefields"
                )))]
                {
                    true
                }
            }
            && {
                #[cfg(any(
                    feature = "cap-rings",
                    feature = "cap-descriptors",
                    feature = "cap-smiles",
                    feature = "cap-sanitize",
                    feature = "cap-hydrogens",
                    feature = "cap-kekulize",
                    feature = "cap-aromaticity",
                    feature = "cap-fingerprints"
                ))]
                {
                    #[cfg(any(feature = "cap-rings", feature = "cap-fingerprints"))]
                    {
                        self.rings.is_none() && self.ring_families.is_none()
                    }
                    #[cfg(not(any(feature = "cap-rings", feature = "cap-fingerprints")))]
                    {
                        self.rings.is_none()
                    }
                }
                #[cfg(not(any(
                    feature = "cap-rings",
                    feature = "cap-descriptors",
                    feature = "cap-smiles",
                    feature = "cap-sanitize",
                    feature = "cap-hydrogens",
                    feature = "cap-kekulize",
                    feature = "cap-aromaticity",
                    feature = "cap-fingerprints"
                )))]
                {
                    true
                }
            }
    }

    pub(crate) fn valid_states(&self) -> DerivedState {
        self.valid
    }

    pub(crate) fn mark_valid(&mut self, states: DerivedState) {
        self.valid = self.valid.union(states);
    }

    pub(crate) fn clear(&mut self, states: DerivedState) {
        self.valid = self.valid.difference(states);
        #[cfg(any(
            feature = "cap-valence",
            feature = "cap-fingerprints",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors",
            feature = "cap-forcefields"
        ))]
        if states.intersects(DerivedState::VALENCE) {
            self.valence = None;
        }
        #[cfg(any(
            feature = "cap-rings",
            feature = "cap-descriptors",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-hydrogens",
            feature = "cap-kekulize",
            feature = "cap-aromaticity",
            feature = "cap-fingerprints"
        ))]
        if states.intersects(DerivedState::RINGS) {
            self.rings = None;
        }
        #[cfg(any(feature = "cap-rings", feature = "cap-fingerprints"))]
        if states.intersects(DerivedState::RING_FAMILIES) {
            self.ring_families = None;
        }
    }

    #[cfg(any(
        feature = "cap-valence",
        feature = "cap-fingerprints",
        feature = "cap-hydrogens",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-descriptors",
        feature = "cap-forcefields"
    ))]
    pub(crate) fn install_valence_assignment(
        &mut self,
        assignment: cosmolkit_core::ValenceAssignment,
    ) {
        self.valence = Some(assignment);
    }

    #[cfg(any(
        feature = "cap-valence",
        feature = "cap-fingerprints",
        feature = "cap-hydrogens",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-descriptors",
        feature = "cap-forcefields"
    ))]
    pub(crate) fn valence_assignment(&self) -> Option<&cosmolkit_core::ValenceAssignment> {
        if self.valid.contains(DerivedState::VALENCE) {
            self.valence.as_ref()
        } else {
            None
        }
    }

    #[cfg(any(
        feature = "cap-rings",
        feature = "cap-descriptors",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-hydrogens",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-fingerprints"
    ))]
    pub(crate) fn install_ring_info(&mut self, rings: cosmolkit_core::RingInfo) {
        self.rings = Some(rings);
    }

    #[cfg(any(
        feature = "cap-rings",
        feature = "cap-descriptors",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-hydrogens",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-fingerprints"
    ))]
    pub(crate) fn ring_info(&self) -> Option<&cosmolkit_core::RingInfo> {
        self.rings.as_ref()
    }

    /// Valid-gated ordinary-ring read: an installed payload without the
    /// RINGS validity bit is invisible to consumers. Malformed pairs are
    /// rejected by construction/commit validation instead of being served.
    #[cfg(any(
        feature = "cap-rings",
        feature = "cap-descriptors",
        feature = "cap-smiles",
        feature = "cap-sanitize",
        feature = "cap-hydrogens",
        feature = "cap-kekulize",
        feature = "cap-aromaticity",
        feature = "cap-fingerprints"
    ))]
    pub(crate) fn valid_ring_info(&self) -> Option<&cosmolkit_core::RingInfo> {
        if self.valid.contains(DerivedState::RINGS) {
            self.rings.as_ref()
        } else {
            None
        }
    }

    #[cfg(any(feature = "cap-rings", feature = "cap-fingerprints"))]
    pub(crate) fn install_ring_family_info(&mut self, families: cosmolkit_core::RingInfo) {
        self.ring_families = Some(families);
    }

    #[cfg(any(feature = "cap-rings", feature = "cap-fingerprints"))]
    pub(crate) fn ring_family_info(&self) -> Option<&cosmolkit_core::RingInfo> {
        self.ring_families.as_ref()
    }

    pub(crate) fn validate_for_atom_count(&self, atom_count: usize) -> Result<(), OperationError> {
        #[cfg(any(
            feature = "cap-valence",
            feature = "cap-fingerprints",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors",
            feature = "cap-forcefields"
        ))]
        {
            let valid = self.valid.contains(DerivedState::VALENCE);
            match (valid, self.valence.as_ref()) {
                (false, None) => {}
                (true, Some(assignment)) => {
                    for (field, actual) in [
                        ("explicit_valence", assignment.explicit_valence.len()),
                        ("implicit_hydrogens", assignment.implicit_hydrogens.len()),
                    ] {
                        if actual != atom_count {
                            return Err(OperationError::InvalidDerivedCache {
                                state: "valence",
                                field,
                                actual,
                                expected: atom_count,
                            });
                        }
                    }
                }
                (true, None) => {
                    return Err(OperationError::InvalidDerivedCache {
                        state: "valence",
                        field: "assignment",
                        actual: 0,
                        expected: 1,
                    });
                }
                (false, Some(_)) => {
                    return Err(OperationError::InvalidDerivedCache {
                        state: "valence",
                        field: "validity_bit",
                        actual: 0,
                        expected: 1,
                    });
                }
            }
        }
        let _ = atom_count;
        Ok(())
    }

    pub(crate) fn validate_for_topology(
        &self,
        topology: &TopologyBlock,
    ) -> Result<(), OperationError> {
        self.validate_for_atom_count(topology.atoms.len())?;
        #[cfg(any(
            feature = "cap-rings",
            feature = "cap-descriptors",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-hydrogens",
            feature = "cap-kekulize",
            feature = "cap-aromaticity",
            feature = "cap-fingerprints"
        ))]
        {
            let valid = self.valid.contains(DerivedState::RINGS);
            match (valid, self.rings.as_ref()) {
                (false, None) => {}
                (true, Some(rings)) => {
                    if !rings.is_initialized() {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "rings",
                            field: "initialized",
                            actual: 0,
                            expected: 1,
                        });
                    }
                    if rings.atom_rings().len() != rings.bond_rings().len() {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "rings",
                            field: "ring_rows",
                            actual: rings.bond_rings().len(),
                            expected: rings.atom_rings().len(),
                        });
                    }
                    if rings.are_ring_families_initialized()
                        || !rings.atom_ring_families().is_empty()
                        || !rings.bond_ring_families().is_empty()
                    {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "rings",
                            field: "ring_families",
                            actual: 1,
                            expected: 0,
                        });
                    }
                    for (atoms, bonds) in rings.atom_rings().iter().zip(rings.bond_rings()) {
                        if atoms.len() != bonds.len() {
                            return Err(OperationError::InvalidDerivedCache {
                                state: "rings",
                                field: "ring_size",
                                actual: bonds.len(),
                                expected: atoms.len(),
                            });
                        }
                    }
                    for atom in rings.atom_rings().iter().flatten() {
                        if atom.index() >= topology.atoms.len() {
                            return Err(OperationError::InvalidDerivedCache {
                                state: "rings",
                                field: "atom_id",
                                actual: atom.index(),
                                expected: topology.atoms.len(),
                            });
                        }
                    }
                    for bond in rings.bond_rings().iter().flatten() {
                        if bond.index() >= topology.bonds.len() {
                            return Err(OperationError::InvalidDerivedCache {
                                state: "rings",
                                field: "bond_id",
                                actual: bond.index(),
                                expected: topology.bonds.len(),
                            });
                        }
                    }
                }
                (true, None) => {
                    return Err(OperationError::InvalidDerivedCache {
                        state: "rings",
                        field: "assignment",
                        actual: 0,
                        expected: 1,
                    });
                }
                (false, Some(_)) => {
                    return Err(OperationError::InvalidDerivedCache {
                        state: "rings",
                        field: "validity_bit",
                        actual: 0,
                        expected: 1,
                    });
                }
            }
        }
        #[cfg(feature = "cap-rings")]
        {
            let valid = self.valid.contains(DerivedState::RING_FAMILIES);
            match (valid, self.ring_families.as_ref()) {
                (false, None) => {}
                (true, Some(families)) => {
                    if !families.is_initialized() || !families.are_ring_families_initialized() {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "ring_families",
                            field: "initialized",
                            actual: 0,
                            expected: 1,
                        });
                    }
                    if families.find_type() != cosmolkit_core::RingFindType::OtherOrUnknown {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "ring_families",
                            field: "find_type",
                            actual: families.find_type() as usize,
                            expected: cosmolkit_core::RingFindType::OtherOrUnknown as usize,
                        });
                    }
                    if !families.atom_rings().is_empty() || !families.bond_rings().is_empty() {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "ring_families",
                            field: "ordinary_rings",
                            actual: families.atom_rings().len() + families.bond_rings().len(),
                            expected: 0,
                        });
                    }
                    if families.atom_ring_families().len() != families.bond_ring_families().len() {
                        return Err(OperationError::InvalidDerivedCache {
                            state: "ring_families",
                            field: "family_rows",
                            actual: families.bond_ring_families().len(),
                            expected: families.atom_ring_families().len(),
                        });
                    }
                    for atom in families.atom_ring_families().iter().flatten() {
                        if atom.index() >= topology.atoms.len() {
                            return Err(OperationError::InvalidDerivedCache {
                                state: "ring_families",
                                field: "atom_id",
                                actual: atom.index(),
                                expected: topology.atoms.len(),
                            });
                        }
                    }
                    for bond in families.bond_ring_families().iter().flatten() {
                        if bond.index() >= topology.bonds.len() {
                            return Err(OperationError::InvalidDerivedCache {
                                state: "ring_families",
                                field: "bond_id",
                                actual: bond.index(),
                                expected: topology.bonds.len(),
                            });
                        }
                    }
                    families
                        .num_relevant_cycles()
                        .map_err(OperationError::Rings)?;
                }
                (true, None) => {
                    return Err(OperationError::InvalidDerivedCache {
                        state: "ring_families",
                        field: "assignment",
                        actual: 0,
                        expected: 1,
                    });
                }
                (false, Some(_)) => {
                    return Err(OperationError::InvalidDerivedCache {
                        state: "ring_families",
                        field: "validity_bit",
                        actual: 0,
                        expected: 1,
                    });
                }
            }
        }
        Ok(())
    }
}

#[derive(Clone)]
struct MoleculeState {
    topology: Arc<TopologyBlock>,
    coordinates: Arc<CoordinateBlock>,
    properties: Arc<MoleculeProperties>,
    derived_cache: Arc<DerivedCacheBlock>,
}

impl MoleculeState {
    fn validate_parts(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Result<(), OperationError> {
        topology
            .validate()
            .map_err(|error: TopologyValidationError| OperationError::InvalidTopology(error))?;
        coordinates
            .validate_for_atom_count(topology.atoms.len())
            .map_err(OperationError::InvalidCoordinates)?;
        for property_list in properties.sdf_property_lists() {
            let (target, expected) = match property_list.target() {
                SdfPropertyListTarget::Atom => ("atom", topology.atoms.len()),
                SdfPropertyListTarget::Bond => ("bond", topology.bonds.len()),
            };
            if property_list.values().len() != expected {
                return Err(OperationError::InvalidPropertyList {
                    target,
                    name: property_list.name().to_owned(),
                    values: property_list.values().len(),
                    expected,
                });
            }
        }
        Ok(())
    }

    fn try_new(
        topology: TopologyBlock,
        mut coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        coordinates.source_coordinate_dim = if coordinates.conformers_3d.is_empty() {
            (!coordinates.conformers_2d.is_empty()).then_some(CoordinateDimension::TwoD)
        } else {
            Some(CoordinateDimension::ThreeD)
        };
        Self::validate_parts(&topology, &coordinates, &properties)?;

        Ok(Self {
            topology: Arc::new(topology),
            coordinates: Arc::new(coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::new(DerivedCacheBlock::default()),
        })
    }

    fn try_new_with_cache(
        topology: Arc<TopologyBlock>,
        coordinates: Arc<CoordinateBlock>,
        properties: Arc<MoleculeProperties>,
        derived_cache: Arc<DerivedCacheBlock>,
    ) -> Result<Self, OperationError> {
        Self::validate_parts(&topology, &coordinates, &properties)?;
        derived_cache.validate_for_topology(&topology)?;
        Ok(Self {
            topology,
            coordinates,
            properties,
            derived_cache,
        })
    }
}

/// The authoritative molecule value exposed by the top-level crate.
///
/// The model blocks are deliberately private. Algorithms can only receive
/// detached blocks through an operation context, and only this runtime can
/// install a validated result back into a live molecule.
pub struct Molecule {
    state: Arc<MoleculeState>,
    // Private computed descriptor properties; never live topology, coordinate,
    // property or derived-validity mutation authority for algorithm owners.
    #[cfg(feature = "cap-descriptors")]
    descriptor_queries: std::sync::Mutex<cosmolkit_descriptors::DescriptorComputedState>,
    #[cfg(feature = "cap-descriptors")]
    descriptor_queries_poisoned: bool,
}

impl Clone for Molecule {
    fn clone(&self) -> Self {
        // RDKit source (ROMol.cpp, initFromOther non-quick copy):
        // RDKit✔️✔️: d_props = other.d_props;
        // Copy only detached computed rows; the four runtime blocks share.
        // Infallible Clone retains poisoned rows AND the poison condition.
        // It never presents damaged state as a successful query or defaults it.
        #[cfg(feature = "cap-descriptors")]
        let (memo, poisoned) = match self.descriptor_queries.lock() {
            Ok(rows) => (rows.clone(), self.descriptor_queries_poisoned),
            Err(rows) => (rows.into_inner().clone(), true),
        };
        Self {
            state: Arc::clone(&self.state),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries: std::sync::Mutex::new(memo),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries_poisoned: poisoned,
        }
    }
}

impl fmt::Debug for Molecule {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Molecule")
            .field("topology", self.state.topology.as_ref())
            .field("coordinates", self.state.coordinates.as_ref())
            .field("properties", self.state.properties.as_ref())
            .field(
                "derived_cache_is_empty",
                &self.state.derived_cache.is_empty(),
            )
            .finish()
    }
}

impl PartialEq for Molecule {
    fn eq(&self, other: &Self) -> bool {
        self.state.topology == other.state.topology
            && self.state.coordinates == other.state.coordinates
            && self.state.properties == other.state.properties
    }
}

impl Molecule {
    /// Query-specific detached memo borrow. This grants no authority over the
    /// four operation-managed runtime blocks. Algorithms see only typed state.
    #[cfg(feature = "cap-descriptors")]
    pub(crate) fn descriptor_queries_runtime(
        &self,
    ) -> Result<
        std::sync::MutexGuard<'_, cosmolkit_descriptors::DescriptorComputedState>,
        crate::DescriptorReadError,
    > {
        if self.descriptor_queries_poisoned {
            return Err(crate::DescriptorReadError::CachePoisoned);
        }
        self.descriptor_queries
            .lock()
            .map_err(|_| crate::DescriptorReadError::CachePoisoned)
    }

    /// Private operation snapshot: copy computed rows only when the existing
    /// source policy preserves them. Clear policies allocate no cached vectors.
    /// A poisoned query state is a structural operation error before mutation.
    pub(crate) fn operation_snapshot_runtime(
        &self,
        preserve_queries: bool,
        operation: &'static str,
    ) -> Result<Self, OperationError> {
        #[cfg(feature = "cap-descriptors")]
        let memo = {
            let rows = self.descriptor_queries_runtime().map_err(|_| {
                OperationError::OperationContract {
                    operation,
                    field: "descriptor_query_cache",
                    issue: "descriptor query cache is poisoned",
                    expected: 0,
                    actual: 1,
                }
            })?;
            // RDKit✔️✔️: d_props = other.d_props;
            // The declared clearComputedProps policy discards these properties;
            // do not copy vectors just to clear them in the same transaction.
            if preserve_queries {
                rows.clone()
            } else {
                Default::default()
            }
        };
        #[cfg(not(feature = "cap-descriptors"))]
        let _ = (preserve_queries, operation);
        Ok(Self {
            state: Arc::clone(&self.state),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries: std::sync::Mutex::new(memo),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries_poisoned: false,
        })
    }

    /// Move the already copied transaction memo into its validated result.
    /// The operation snapshot is private, unpoisoned and never queried by bodies.
    #[cfg(feature = "cap-descriptors")]
    pub(crate) fn take_descriptor_queries_runtime(
        &mut self,
    ) -> std::sync::Mutex<cosmolkit_descriptors::DescriptorComputedState> {
        std::mem::take(&mut self.descriptor_queries)
    }

    #[cfg(feature = "cap-descriptors")]
    pub(crate) fn install_descriptor_queries_runtime(
        &mut self,
        memo: std::sync::Mutex<cosmolkit_descriptors::DescriptorComputedState>,
    ) {
        self.descriptor_queries = memo;
        self.descriptor_queries_poisoned = false;
    }

    /// Creates an empty, structurally valid molecule.
    #[must_use]
    pub fn new() -> Self {
        Self::from_validated_parts(
            TopologyBlock::default(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .expect("the statically empty molecule blocks are valid")
    }

    /// Constructs a live molecule after validating all local block invariants.
    pub fn from_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        MoleculeBuilder::from_parts(topology, coordinates, properties).build()
    }

    pub(crate) fn from_validated_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        Ok(Self {
            state: Arc::new(MoleculeState::try_new(topology, coordinates, properties)?),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries: std::sync::Mutex::default(),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries_poisoned: false,
        })
    }

    pub(crate) fn from_runtime_parts(
        topology: Arc<TopologyBlock>,
        coordinates: Arc<CoordinateBlock>,
        properties: Arc<MoleculeProperties>,
        derived_cache: Arc<DerivedCacheBlock>,
    ) -> Result<Self, OperationError> {
        Ok(Self {
            state: Arc::new(MoleculeState::try_new_with_cache(
                topology,
                coordinates,
                properties,
                derived_cache,
            )?),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries: std::sync::Mutex::default(),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries_poisoned: false,
        })
    }

    /// Constructor-only transport of final detached chemistry state. Cache
    /// authority remains here, not in the parser or algorithm owner.
    #[cfg(feature = "cap-smiles")]
    pub(super) fn from_smiles_parts_with_derived_state(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        valence: Option<cosmolkit_core::ValenceAssignment>,
        rings: Option<cosmolkit_core::RingInfo>,
    ) -> Result<Self, OperationError> {
        let mut state = MoleculeState::try_new(topology, coordinates, properties)?;
        let mut cache = DerivedCacheBlock::default();
        if let Some(assignment) = valence {
            cache.install_valence_assignment(assignment);
            cache.mark_valid(DerivedState::VALENCE);
        }
        // Moved final ring carrier: an initialized Some (including
        // initialized-empty and Other quality) is stored and marked RINGS;
        // None or an uninitialized reset clears ordinary storage and never
        // marks valid or finds rings. No finder runs at this seam.
        match rings {
            Some(rings) if rings.is_initialized() => {
                // Actual-site cfg(test) observation of the moved row
                // buffers immediately before installation.
                #[cfg(test)]
                ring_install_probe::record(&rings);
                cache.install_ring_info(rings);
                cache.mark_valid(DerivedState::RINGS);
            }
            _ => {
                cache.clear(DerivedState::RINGS);
            }
        }
        cache.validate_for_topology(&state.topology)?;
        state.derived_cache = Arc::new(cache);
        Ok(Self {
            state: Arc::new(state),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries: std::sync::Mutex::default(),
            #[cfg(feature = "cap-descriptors")]
            descriptor_queries_poisoned: false,
        })
    }

    /// Returns detached semantic blocks for checked construction of a new value.
    #[must_use]
    pub fn to_builder(&self) -> MoleculeBuilder {
        MoleculeBuilder::from_parts(
            self.state.topology.as_ref().clone(),
            self.state.coordinates.as_ref().clone(),
            self.state.properties.as_ref().clone(),
        )
    }

    /// Returns the immutable topology value.
    #[must_use]
    pub fn topology(&self) -> &TopologyBlock {
        self.state.topology.as_ref()
    }

    /// Returns the coordinates of the first stored 2D conformer, if present.
    ///
    /// The rows follow atom order. This borrows existing coordinates without
    /// generating a layout or falling back to a 3D conformer.
    ///
    /// Coordinate access on a molecule must specify the dimension:
    ///
    /// ```compile_fail,E0599
    /// let molecule = cosmolkit::Molecule::new();
    /// molecule.coordinates();
    /// ```
    ///
    /// Complete coordinate blocks are private runtime state:
    ///
    /// ```compile_fail,E0624
    /// let molecule = cosmolkit::Molecule::new();
    /// molecule.coordinate_block_runtime();
    /// ```
    #[must_use]
    pub fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.state
            .coordinates
            .conformers_2d
            .first()
            .map(Conformer2D::coordinates)
    }

    /// Complete block access for runtime validation and detached owner calls.
    pub(crate) fn coordinate_block_runtime(&self) -> &CoordinateBlock {
        self.state.coordinates.as_ref()
    }

    /// Returns the immutable molecule properties.
    #[must_use]
    pub fn properties(&self) -> &MoleculeProperties {
        self.state.properties.as_ref()
    }

    /// Returns the number of atoms in the authoritative topology.
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.state.topology.atoms.len()
    }

    /// Returns the number of bonds in the authoritative topology.
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.state.topology.bonds.len()
    }

    /// Returns all atoms in canonical row order.
    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.state.topology.atoms
    }

    /// Returns all bonds in canonical row order.
    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.state.topology.bonds
    }

    /// Returns the atom with the requested stable row identifier.
    #[must_use]
    pub fn atom(&self, atom_id: AtomId) -> Option<&Atom> {
        self.state.topology.atoms.get(atom_id.index())
    }

    /// Returns the bond with the requested stable row identifier.
    #[must_use]
    pub fn bond(&self, bond_id: BondId) -> Option<&Bond> {
        self.state.topology.bonds.get(bond_id.index())
    }

    /// Returns all stored 3D conformers in their original order, without cloning.
    ///
    /// The public API has no dimension-ambiguous conformer tuple:
    ///
    /// ```compile_fail,E0599
    /// let molecule = cosmolkit::Molecule::new();
    /// molecule.conformers();
    /// ```
    #[must_use]
    pub fn conformers_3d(&self) -> &[Conformer3D] {
        &self.state.coordinates.conformers_3d
    }

    /// Returns an ordinary molecule property by key.
    #[must_use]
    pub fn property(&self, key: &str) -> Option<&str> {
        self.state.properties.prop(key)
    }

    pub(crate) fn derived_cache_runtime(&self) -> &DerivedCacheBlock {
        self.state.derived_cache.as_ref()
    }

    pub(crate) fn topology_arc_runtime(&self) -> Arc<TopologyBlock> {
        Arc::clone(&self.state.topology)
    }

    pub(crate) fn coordinates_arc_runtime(&self) -> Arc<CoordinateBlock> {
        Arc::clone(&self.state.coordinates)
    }

    pub(crate) fn properties_arc_runtime(&self) -> Arc<MoleculeProperties> {
        Arc::clone(&self.state.properties)
    }

    pub(crate) fn derived_cache_arc_runtime(&self) -> Arc<DerivedCacheBlock> {
        Arc::clone(&self.state.derived_cache)
    }

    #[cfg(test)]
    pub(crate) const fn runtime_constructions(&self) -> usize {
        1
    }
}

/// Actual-site cfg(test) history for constructor ring installs: append-only,
/// never reset. Each entry records BOTH row-buffer addresses of the final
/// carrier immediately before it is moved into the derived cache. The seam
/// itself runs NO finder.
#[cfg(all(test, feature = "cap-smiles"))]
pub(crate) mod ring_install_probe {
    use std::cell::RefCell;

    thread_local! {
        static INSTALLS: RefCell<Vec<(usize, usize)>> = const { RefCell::new(Vec::new()) };
    }

    pub(crate) fn record(rings: &cosmolkit_core::RingInfo) {
        INSTALLS.with(|installs| {
            installs.borrow_mut().push((
                rings.atom_rings().as_ptr() as usize,
                rings.bond_rings().as_ptr() as usize,
            ));
        });
    }

    pub(crate) fn len() -> usize {
        INSTALLS.with(|installs| installs.borrow().len())
    }

    pub(crate) fn after(baseline: usize) -> Vec<(usize, usize)> {
        INSTALLS.with(|installs| {
            let installs = installs.borrow();
            installs[baseline.min(installs.len())..].to_vec()
        })
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

/// C1: twenty real constructor calls installing the moved final ring
/// carrier. RINGS validity is proven independent of VALENCE; the frozen
/// final-quality table is checked against the LIVE cache; actual install
/// row-buffer addresses are observed at the seam and preserved.
#[cfg(all(test, feature = "cap-smiles"))]
mod ring_live_constructor_tests {
    use super::*;

    #[test]
    fn ring_live_constructor_state_move_and_valence_independence() {
        let mut calls = 0usize;
        for input in ["", "CC", "c1ccccc1", "[H]C1CCCCC1", "[2H]C1CCCCC1"] {
            for (profile, sanitize, remove_hydrogens) in [
                ("bothfalse", false, false),
                ("remove-only", false, true),
                ("sanitize-only", true, false),
                ("bothtrue", true, true),
            ] {
                let label = format!("{profile}/{input:?}");
                let params = cosmolkit_smiles::SmilesParseParams {
                    sanitize,
                    remove_hydrogens,
                    ..cosmolkit_smiles::SmilesParseParams::default()
                };
                let probe_before = ring_install_probe::len();
                let molecule = Molecule::from_smiles_with_params(input, &params)
                    .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                calls += 1;
                let installs = ring_install_probe::after(probe_before);

                // Literal final topology counts and element/isotope identity.
                let deuterium = input == "[2H]C1CCCCC1";
                let (want_atoms, want_bonds, hydrogen_kept) = match input {
                    "" => (0usize, 0usize, false),
                    "CC" => (2, 1, false),
                    "c1ccccc1" => (6, 6, false),
                    _ if deuterium => (7, 7, true),
                    _ => {
                        if remove_hydrogens {
                            (6, 6, false)
                        } else {
                            (7, 7, true)
                        }
                    }
                };
                assert_eq!(molecule.num_atoms(), want_atoms, "{label}: atoms");
                assert_eq!(molecule.num_bonds(), want_bonds, "{label}: bonds");
                match input {
                    "" | "CC" => {}
                    "c1ccccc1" => {
                        for atom in molecule.atoms() {
                            assert_eq!(
                                atom.element(),
                                cosmolkit_model::Element::C,
                                "{label}: element"
                            );
                            assert!(atom.isotope().is_none(), "{label}: isotope");
                        }
                    }
                    _ => {
                        if hydrogen_kept {
                            assert_eq!(
                                molecule.atoms()[0].element(),
                                cosmolkit_model::Element::H,
                                "{label}: H element"
                            );
                            assert_eq!(
                                molecule.atoms()[0].isotope(),
                                deuterium.then_some(2),
                                "{label}: H isotope"
                            );
                        }
                        let carbon_start = usize::from(hydrogen_kept);
                        for index in carbon_start..want_atoms {
                            assert_eq!(
                                molecule.atoms()[index].element(),
                                cosmolkit_model::Element::C,
                                "{label}: C element {index}"
                            );
                            assert!(
                                molecule.atoms()[index].isotope().is_none(),
                                "{label}: C isotope {index}"
                            );
                        }
                    }
                }

                // Frozen final ring-quality table for the LIVE cache.
                let cache = molecule.derived_cache_runtime();
                let expected_quality = match profile {
                    "bothfalse" => None,
                    "remove-only" => Some(cosmolkit_core::RingFindType::Fast),
                    "sanitize-only" => Some(cosmolkit_core::RingFindType::SymmSssr),
                    _ if input.is_empty() => Some(cosmolkit_core::RingFindType::Fast),
                    _ => Some(cosmolkit_core::RingFindType::SymmSssr),
                };
                assert_eq!(
                    cache.valid_states().contains(DerivedState::RINGS),
                    expected_quality.is_some(),
                    "{label}: RINGS validity"
                );
                // RINGS is independent of the CK-VALENCE-001 sanitize gate.
                assert_eq!(
                    cache.valid_states().contains(DerivedState::VALENCE),
                    sanitize,
                    "{label}: VALENCE validity"
                );
                match expected_quality {
                    None => {
                        assert!(cache.valid_ring_info().is_none(), "{label}: absent");
                        assert_eq!(installs.len(), 0, "{label}: no seam install");
                        assert!(cache.ring_info().is_none(), "{label}: no storage");
                    }
                    Some(quality) => {
                        let rings = cache.valid_ring_info().expect("{label}: installed");
                        assert!(rings.is_initialized(), "{label}: initialized");
                        assert_eq!(rings.find_type(), quality, "{label}: find type");
                        assert_eq!(rings.atom_row_count(), want_atoms, "{label}: atom dims");
                        assert_eq!(rings.bond_row_count(), want_bonds, "{label}: bond dims");
                        // Exactly one moved install; no finder at the seam.
                        assert_eq!(installs.len(), 1, "{label}: one seam install");
                        if !rings.atom_rings().is_empty() {
                            assert_eq!(
                                rings.atom_rings().as_ptr() as usize,
                                installs[0].0,
                                "{label}: atom buffers moved"
                            );
                            assert_eq!(
                                rings.bond_rings().as_ptr() as usize,
                                installs[0].1,
                                "{label}: bond buffers moved"
                            );
                        }
                        // Ordered paired rows and every membership entry.
                        let cycle = matches!(input, "c1ccccc1" | "[H]C1CCCCC1" | "[2H]C1CCCCC1");
                        if cycle {
                            assert_eq!(rings.atom_rings().len(), 1, "{label}: ring count");
                            assert_eq!(rings.bond_rings().len(), 1, "{label}: bond rows");
                            let start = usize::from(hydrogen_kept);
                            let expected_ids: Vec<usize> = (start..start + 6).collect();
                            let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                                .iter()
                                .map(|atom| atom.index())
                                .collect();
                            atoms_row.sort_unstable();
                            let mut bonds_row: Vec<usize> = rings.bond_rings()[0]
                                .iter()
                                .map(|bond| bond.index())
                                .collect();
                            bonds_row.sort_unstable();
                            assert_eq!(atoms_row, expected_ids, "{label}: cycle atoms");
                            assert_eq!(bonds_row, expected_ids, "{label}: cycle bonds");
                            for index in 0..want_atoms {
                                let expected: &[usize] = if hydrogen_kept && index < start {
                                    &[]
                                } else {
                                    &[0]
                                };
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    expected,
                                    "{label}: atom member {index}"
                                );
                            }
                            for index in 0..want_bonds {
                                let expected: &[usize] = if hydrogen_kept && index < start {
                                    &[]
                                } else {
                                    &[0]
                                };
                                assert_eq!(
                                    rings.bond_members(BondId::new(index)),
                                    expected,
                                    "{label}: bond member {index}"
                                );
                            }
                        } else {
                            assert!(rings.atom_rings().is_empty(), "{label}: rows");
                            assert!(rings.bond_rings().is_empty(), "{label}: bond rows");
                            for index in 0..want_atoms {
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    &[] as &[usize],
                                    "{label}: empty member {index}"
                                );
                            }
                            for index in 0..want_bonds {
                                assert_eq!(
                                    rings.bond_members(BondId::new(index)),
                                    &[] as &[usize],
                                    "{label}: empty bond member {index}"
                                );
                            }
                        }
                    }
                }

                // Ordinary computed metadata and coordinates are preserved
                // by the transport seam.
                if sanitize || remove_hydrogens {
                    assert_eq!(
                        molecule.properties().prop("_StereochemDone"),
                        Some("1"),
                        "{label}: done marker"
                    );
                }
                assert!(molecule.coordinates_2d().is_none(), "{label}: 2d");
                assert!(molecule.conformers_3d().is_empty(), "{label}: 3d");
            }
        }
        assert_eq!(calls, 20, "exact census");
    }
}

/// C6: one real registered-operation lifecycle over the live ring state,
/// builder isolation, malformed-cache rejection at construction, and
/// value/in-place failure semantics.
#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-descriptors",
    feature = "cap-rings",
    feature = "cap-hydrogens",
    feature = "cap-depict",
    feature = "cap-kekulize"
))]
mod ring_live_lifecycle_tests {
    use super::*;

    fn ring_rows(rings: &cosmolkit_core::RingInfo) -> Vec<usize> {
        let mut row: Vec<usize> = rings.atom_rings()[0].iter().map(|a| a.index()).collect();
        row.sort_unstable();
        row
    }

    #[test]
    fn ring_live_lifecycle_operation_sequence() {
        // Default constructor: Symm, count 1.
        let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
        let peer = molecule.clone();
        let peer_cache = peer.derived_cache_runtime().clone();
        assert_eq!(molecule.num_rings().unwrap(), 1, "initial count");
        {
            let rings = molecule.derived_cache_runtime().valid_ring_info().unwrap();
            assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::SymmSssr);
        }
        assert_eq!(molecule.num_rings().unwrap(), 1, "repeat query");

        // Coordinate value operation preserves rows/type and old IDs.
        let with_coords = molecule.with_2d_coordinates().unwrap();
        {
            let rings = with_coords
                .derived_cache_runtime()
                .valid_ring_info()
                .unwrap();
            assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::SymmSssr);
            assert_eq!(ring_rows(rings), vec![0, 1, 2, 3, 4, 5], "old IDs");
        }
        assert!(with_coords.coordinates_2d().is_some(), "2d installed");

        // AddHs preserves ring state through LeafAtomAppend.
        let with_hs = with_coords.with_hydrogens().unwrap();
        {
            let rings = with_hs.derived_cache_runtime().valid_ring_info().unwrap();
            assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::SymmSssr);
            assert_eq!(ring_rows(rings), vec![0, 1, 2, 3, 4, 5], "old IDs");
            assert!(with_hs.num_atoms() > 6, "hydrogens appended");
        }

        // RemoveHs with sanitize=false leaves rings ABSENT: num_rings is
        // the typed error, every other classifier the source-defined 0.
        let removed = with_hs
            .without_hydrogens_with_params(&cosmolkit_core::RemoveHsParams {
                sanitize: false,
                ..cosmolkit_core::RemoveHsParams::default()
            })
            .unwrap();
        assert!(
            matches!(
                removed.num_rings(),
                Err(crate::DescriptorReadError::MissingInitializedRings)
            ),
            "absent after sanitize=false removal"
        );
        assert_eq!(removed.num_heterocycles().unwrap(), 0);
        assert_eq!(removed.num_aromatic_rings().unwrap(), 0);
        assert_eq!(removed.num_saturated_rings().unwrap(), 0);
        assert_eq!(removed.num_aliphatic_rings().unwrap(), 0);
        assert_eq!(removed.num_aromatic_heterocycles().unwrap(), 0);
        assert_eq!(removed.num_aromatic_carbocycles().unwrap(), 0);
        assert_eq!(removed.num_aliphatic_heterocycles().unwrap(), 0);
        assert_eq!(removed.num_aliphatic_carbocycles().unwrap(), 0);
        assert_eq!(removed.num_saturated_heterocycles().unwrap(), 0);
        assert_eq!(removed.num_saturated_carbocycles().unwrap(), 0);

        // Explicit assignment installs Fast, count 1.
        let assigned = removed.with_assigned_rings().unwrap();
        {
            let rings = assigned.derived_cache_runtime().valid_ring_info().unwrap();
            assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::Fast);
            assert_eq!(assigned.num_rings().unwrap(), 1, "explicit Fast count");
        }

        // The peer never changed and never gained a writable cache view.
        assert_eq!(peer.derived_cache_runtime(), &peer_cache);
        assert_eq!(peer.num_rings().unwrap(), 1, "peer intact");
        assert!(std::sync::Arc::ptr_eq(
            &peer.derived_cache_arc_runtime(),
            &molecule.derived_cache_arc_runtime()
        ));
    }

    #[test]
    fn ring_live_lifecycle_builder_starts_absent() {
        let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
        let mut builder = molecule.to_builder();
        let _ = builder.add_atom(cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::C));
        let rebuilt = builder.build().unwrap();
        // The builder path starts from a fresh cache: no stale rows are
        // ever transferred from the old molecule.
        assert!(
            matches!(
                rebuilt.num_rings(),
                Err(crate::DescriptorReadError::MissingInitializedRings)
            ),
            "builder starts absent"
        );
        assert_eq!(rebuilt.num_heterocycles().unwrap(), 0);
        assert_eq!(molecule.num_rings().unwrap(), 1, "source untouched");
    }

    #[test]
    fn ring_live_lifecycle_malformed_cache_rejected_at_construction() {
        // Malformed paired rows / out-of-range IDs are UNREACHABLE through
        // public or seam construction (RingInfo exposes no row-level API);
        // the reachable malformed pair is the validity/storage mismatch,
        // rejected by the EXISTING construction validation.
        let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();

        // Validity bit without storage.
        let mut bad = DerivedCacheBlock::default();
        bad.mark_valid(DerivedState::RINGS);
        assert!(
            matches!(
                Molecule::from_runtime_parts(
                    topology.clone(),
                    coordinates.clone(),
                    properties.clone(),
                    Arc::new(bad)
                ),
                Err(OperationError::InvalidDerivedCache {
                    state: "rings",
                    field: "assignment",
                    ..
                })
            ),
            "validity without storage rejected"
        );

        // Storage without the validity bit.
        let mut stale = DerivedCacheBlock::default();
        stale.install_ring_info(cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::SymmSssr,
            6,
            6,
        ));
        assert!(
            matches!(
                Molecule::from_runtime_parts(topology, coordinates, properties, Arc::new(stale)),
                Err(OperationError::InvalidDerivedCache {
                    state: "rings",
                    field: "validity_bit",
                    ..
                })
            ),
            "storage without validity rejected"
        );
    }

    #[test]
    fn ring_live_lifecycle_failure_semantics() {
        // Value failure preserves the source molecule exactly.
        let base = Molecule::from_smiles("c1ccccc1").unwrap();
        let molecule = Molecule::from_smiles_parts_with_derived_state(
            base.topology().clone(),
            base.coordinate_block_runtime().clone(),
            base.properties().clone(),
            None,
            Some(cosmolkit_core::RingInfo::new(
                cosmolkit_core::RingFindType::Fast,
                6,
                6,
            )),
        )
        .unwrap();
        let observer = molecule.clone();
        let error = molecule
            .with_kekulized_bonds_with_params(&crate::KekulizeParams {
                mark_atoms_bonds: true,
                canonical: false,
                max_backtracks: crate::KekulizeParams::default().max_backtracks,
            })
            .unwrap_err();
        assert!(
            matches!(
                &error,
                crate::OperationError::Kekulize(
                    cosmolkit_core::KekulizeError::AromaticAtomOutsideRing { .. }
                )
            ),
            "got {error:?}"
        );
        assert_eq!(molecule, observer, "value failure preserves source");

        // In-place failure is BASIC safety: the receiver stays whole and
        // queryable, not rolled back.
        let mut target = observer.clone();
        let in_place = target.kekulize_bonds_with_params_(&crate::KekulizeParams {
            mark_atoms_bonds: true,
            canonical: false,
            max_backtracks: crate::KekulizeParams::default().max_backtracks,
        });
        assert!(in_place.is_err(), "same typed failure in place");
        assert_eq!(target.num_atoms(), 6, "storage complete");
        assert_eq!(target.num_heterocycles().unwrap(), 0, "queryable");
    }
}

#[cfg(all(test, feature = "cap-valence"))]
mod valence_cache_tests {
    use super::*;

    #[cfg(all(feature = "cap-smiles", feature = "cap-descriptors"))]
    #[test]
    fn valence_transport_smiles_flag_product_and_read_only_descriptors() {
        // Fixed 2026.03.1 rows, including a falsely tagged radical center whose
        // stereo cleanup removes its explicit H and refreshes its valence.
        for (smiles, formula, mass) in [
            ("CCO", "C2H6O", 46.069),
            ("[C@H](C)C", "C3H7", 43.089),
            ("C[C@H](F)C", "C3H7F", 62.087),
            ("[13CH3][NH3+]", "CH6N+", 33.05835484),
            ("c1cc[nH]c1", "C4H5N", 67.091),
            ("[H]OC", "CH4O", 32.042),
            ("", "", 0.0),
        ] {
            for sanitize in [false, true] {
                for remove_hydrogens in [false, true] {
                    let molecule = Molecule::from_smiles_with_params(
                        smiles,
                        &cosmolkit_smiles::SmilesParseParams {
                            sanitize,
                            remove_hydrogens,
                            ..Default::default()
                        },
                    )
                    .unwrap();
                    let cache = molecule.derived_cache_arc_runtime();
                    assert_eq!(
                        cache.valid_states().contains(DerivedState::VALENCE),
                        sanitize,
                        "{smiles} sanitize={sanitize} remove_hydrogens={remove_hydrogens}"
                    );
                    if sanitize {
                        let expected = cosmolkit_core::assign_valence_with_options_for_topology(
                            molecule.topology(),
                            cosmolkit_core::ValenceModel::RdkitLike,
                            false,
                        )
                        .unwrap();
                        assert_eq!(cache.valence_assignment(), Some(&expected), "{smiles}");
                        assert_eq!(
                            molecule.molecular_formula().unwrap(),
                            formula,
                            "{smiles} sanitize={sanitize} remove_hydrogens={remove_hydrogens} atoms={:?}",
                            molecule.atoms()
                        );
                        assert!((molecule.molecular_weight().unwrap() - mass).abs() < 1e-9);
                        let peer = molecule.clone();
                        let topology = molecule.topology_arc_runtime();
                        let coordinates = molecule.coordinates_arc_runtime();
                        let properties = molecule.properties_arc_runtime();
                        for only_heavy in [false, true] {
                            molecule.molecular_weight_with_params(only_heavy).unwrap();
                            molecule
                                .exact_molecular_weight_with_params(only_heavy)
                                .unwrap();
                        }
                        for separate in [false, true] {
                            for abbreviate in [false, true] {
                                molecule
                                    .molecular_formula_with_params(separate, abbreviate)
                                    .unwrap();
                            }
                        }
                        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
                        assert!(Arc::ptr_eq(&cache, &peer.derived_cache_arc_runtime()));
                        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
                        assert!(Arc::ptr_eq(
                            &coordinates,
                            &molecule.coordinates_arc_runtime()
                        ));
                        assert!(Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()));
                    } else {
                        assert!(cache.valence_assignment().is_none());
                    }
                }
            }
        }
    }

    #[cfg(feature = "cap-hydrogens")]
    #[test]
    fn removal_final_valence_obeys_sanitize_and_preserves_source_cache() {
        let topology = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(
                    AtomId::new(0),
                    cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::C),
                ),
                Atom::from_spec(
                    AtomId::new(1),
                    cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::H),
                ),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                cosmolkit_model::BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(1),
                    cosmolkit_model::BondOrder::Single,
                ),
            )],
            vec![],
            vec![],
        )
        .unwrap();
        let source = Molecule::from_parts(
            topology,
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap()
        .with_assigned_valence()
        .unwrap();
        let original_cache = source.derived_cache_arc_runtime();
        assert_eq!(
            original_cache
                .valence_assignment()
                .unwrap()
                .explicit_valence,
            [1, 1]
        );
        for sanitize in [false, true] {
            let params = cosmolkit_core::RemoveHsParams {
                sanitize,
                ..Default::default()
            };
            let output = source.without_hydrogens_with_params(&params).unwrap();
            let cache = output.derived_cache_runtime();
            assert_eq!(
                cache.valid_states().contains(DerivedState::VALENCE),
                sanitize
            );
            if sanitize {
                let expected =
                    cosmolkit_core::assign_valence(output.topology(), &Default::default()).unwrap();
                assert_eq!(cache.valence_assignment(), Some(&expected));
            } else {
                assert!(cache.valence_assignment().is_none());
                // Repeating a no-op deletion must not revive an invalid cache.
                let repeated = output.without_hydrogens_with_params(&params).unwrap();
                assert!(
                    repeated
                        .derived_cache_runtime()
                        .valence_assignment()
                        .is_none()
                );
                let recomputed = output.with_assigned_valence().unwrap();
                assert!(
                    recomputed
                        .derived_cache_runtime()
                        .valid_states()
                        .contains(DerivedState::VALENCE)
                );
                assert!(
                    output
                        .derived_cache_runtime()
                        .valence_assignment()
                        .is_none()
                );
            }
            cache.validate_for_atom_count(output.num_atoms()).unwrap();
            let mut in_place = source.clone();
            in_place.remove_hydrogens_with_params_(&params).unwrap();
            assert_eq!(in_place.derived_cache_runtime(), cache);
            assert!(Arc::ptr_eq(
                &source.derived_cache_arc_runtime(),
                &original_cache
            ));
        }
    }

    #[test]
    fn valence_payload_and_validity_bit_are_one_validated_cache_state() {
        let mut cache = DerivedCacheBlock::default();
        cache.install_valence_assignment(cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![1, 2],
            implicit_hydrogens: vec![3, 2],
        });
        assert!(matches!(
            cache.validate_for_atom_count(2),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "validity_bit",
                ..
            })
        ));
        cache.mark_valid(DerivedState::VALENCE);
        assert_eq!(cache.validate_for_atom_count(2), Ok(()));
        assert_eq!(
            cache.valence_assignment().unwrap().explicit_valence,
            vec![1, 2]
        );
        assert!(matches!(
            cache.validate_for_atom_count(3),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "explicit_valence",
                actual: 2,
                expected: 3,
            })
        ));
        cache.clear(DerivedState::VALENCE);
        assert_eq!(cache.valence_assignment(), None);
        assert_eq!(cache.validate_for_atom_count(2), Ok(()));
    }
}

#[cfg(all(test, feature = "cap-forcefields"))]
mod forcefields_valence_cache_tests {
    use super::*;

    #[test]
    fn uff_param_p04_forcefields_stores_reads_validates_and_clears_assignment() {
        let mut cache = DerivedCacheBlock::default();
        assert_eq!(cache.validate_for_atom_count(2), Ok(()));
        assert!(!cache.valid_states().contains(DerivedState::VALENCE));
        assert_eq!(cache.valence_assignment(), None);

        cache.install_valence_assignment(cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![1, 2],
            implicit_hydrogens: vec![3, 2],
        });
        cache.mark_valid(DerivedState::VALENCE);

        assert!(cache.valid_states().contains(DerivedState::VALENCE));
        let stored = cache
            .valence_assignment()
            .expect("forcefields-only cache keeps the detached assignment");
        assert_eq!(stored.explicit_valence, [1, 2]);
        assert_eq!(stored.implicit_hydrogens, [3, 2]);
        assert_eq!(cache.validate_for_atom_count(2), Ok(()));

        cache.clear(DerivedState::VALENCE);
        assert!(!cache.valid_states().contains(DerivedState::VALENCE));
        assert_eq!(cache.valence_assignment(), None);
        assert_eq!(cache.validate_for_atom_count(2), Ok(()));
    }

    #[test]
    fn uff_param_p04_forcefields_rejects_inconsistent_cache_pairs() {
        let mut missing_validity = DerivedCacheBlock::default();
        missing_validity.install_valence_assignment(cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![1],
            implicit_hydrogens: vec![0],
        });
        assert!(matches!(
            missing_validity.validate_for_atom_count(1),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "validity_bit",
                actual: 0,
                expected: 1,
            })
        ));

        let mut missing_assignment = DerivedCacheBlock::default();
        missing_assignment.mark_valid(DerivedState::VALENCE);
        assert!(matches!(
            missing_assignment.validate_for_atom_count(1),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "assignment",
                actual: 0,
                expected: 1,
            })
        ));

        let mut explicit_length = DerivedCacheBlock::default();
        explicit_length.install_valence_assignment(cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![1, 2],
            implicit_hydrogens: vec![0],
        });
        explicit_length.mark_valid(DerivedState::VALENCE);
        assert!(matches!(
            explicit_length.validate_for_atom_count(1),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "explicit_valence",
                actual: 2,
                expected: 1,
            })
        ));

        let mut implicit_length = DerivedCacheBlock::default();
        implicit_length.install_valence_assignment(cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![1],
            implicit_hydrogens: vec![0, 2],
        });
        implicit_length.mark_valid(DerivedState::VALENCE);
        assert!(matches!(
            implicit_length.validate_for_atom_count(1),
            Err(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "implicit_hydrogens",
                actual: 2,
                expected: 1,
            })
        ));
    }
}
