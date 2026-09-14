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
    #[cfg(any(feature = "valence", feature = "hydrogens"))]
    valence: Option<cosmolkit_core::ValenceAssignment>,
    #[cfg(feature = "rings")]
    rings: Option<cosmolkit_core::RingInfo>,
    #[cfg(feature = "rings")]
    ring_families: Option<cosmolkit_core::RingInfo>,
}

impl DerivedCacheBlock {
    fn is_empty(&self) -> bool {
        self.valid == DerivedState::NONE
            && {
                #[cfg(any(feature = "valence", feature = "hydrogens"))]
                {
                    self.valence.is_none()
                }
                #[cfg(not(any(feature = "valence", feature = "hydrogens")))]
                {
                    true
                }
            }
            && {
                #[cfg(feature = "rings")]
                {
                    self.rings.is_none() && self.ring_families.is_none()
                }
                #[cfg(not(feature = "rings"))]
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
        #[cfg(any(feature = "valence", feature = "hydrogens"))]
        if states.intersects(DerivedState::VALENCE) {
            self.valence = None;
        }
        #[cfg(feature = "rings")]
        if states.intersects(DerivedState::RINGS) {
            self.rings = None;
        }
        #[cfg(feature = "rings")]
        if states.intersects(DerivedState::RING_FAMILIES) {
            self.ring_families = None;
        }
    }

    #[cfg(any(feature = "valence", feature = "hydrogens"))]
    pub(crate) fn install_valence_assignment(
        &mut self,
        assignment: cosmolkit_core::ValenceAssignment,
    ) {
        self.valence = Some(assignment);
    }

    #[cfg(any(feature = "valence", feature = "hydrogens"))]
    pub(crate) fn valence_assignment(&self) -> Option<&cosmolkit_core::ValenceAssignment> {
        self.valence.as_ref()
    }

    #[cfg(feature = "rings")]
    pub(crate) fn install_ring_info(&mut self, rings: cosmolkit_core::RingInfo) {
        self.rings = Some(rings);
    }

    #[cfg(feature = "rings")]
    pub(crate) fn ring_info(&self) -> Option<&cosmolkit_core::RingInfo> {
        self.rings.as_ref()
    }

    #[cfg(feature = "rings")]
    pub(crate) fn install_ring_family_info(&mut self, families: cosmolkit_core::RingInfo) {
        self.ring_families = Some(families);
    }

    #[cfg(feature = "rings")]
    pub(crate) fn ring_family_info(&self) -> Option<&cosmolkit_core::RingInfo> {
        self.ring_families.as_ref()
    }

    pub(crate) fn validate_for_atom_count(&self, atom_count: usize) -> Result<(), OperationError> {
        #[cfg(any(feature = "valence", feature = "hydrogens"))]
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
        #[cfg(feature = "rings")]
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
#[derive(Clone)]
pub struct Molecule {
    state: Arc<MoleculeState>,
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

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(all(test, feature = "valence"))]
mod valence_cache_tests {
    use super::*;

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
