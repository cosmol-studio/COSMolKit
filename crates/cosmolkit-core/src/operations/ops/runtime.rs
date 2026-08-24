//! Runtime capability object for registered molecule operations.
//!
//! This module owns the working molecule state. Operation bodies live in sibling
//! modules, so Rust privacy prevents them from accessing `OpParts` fields.

use super::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperationTrace {
    touched_blocks: BlockSet,
    claimed_write_blocks: BlockSet,
    recorded_topology_edit: TopologyEditKind,
    remapped_blocks: BlockSet,
    preserved_cache: DerivedState,
    cleared_cache: DerivedState,
    updated_cache: DerivedState,
}

impl OperationTrace {
    #[must_use]
    pub const fn touched_blocks(&self) -> BlockSet {
        self.touched_blocks
    }

    #[must_use]
    pub const fn remapped_blocks(&self) -> BlockSet {
        self.remapped_blocks
    }

    #[must_use]
    pub const fn preserved_cache(&self) -> DerivedState {
        self.preserved_cache
    }

    #[must_use]
    pub const fn cleared_cache(&self) -> DerivedState {
        self.cleared_cache
    }

    #[must_use]
    pub const fn updated_cache(&self) -> DerivedState {
        self.updated_cache
    }
}

#[cfg(feature = "op-contracts")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BlockLifecycle {
    Available,
    Begun,
    Committed,
}

#[cfg(feature = "op-contracts")]
#[derive(Debug, Clone, PartialEq, Eq)]
struct ContractSourceSnapshot {
    atom_count: usize,
    bond_endpoints: Vec<(AtomId, AtomId)>,
    cip_state: CipStateSnapshot,
}

#[cfg(feature = "op-contracts")]
#[derive(Debug, Clone, PartialEq, Eq)]
struct CipStateSnapshot {
    computed: CipPropertySnapshot,
    atom_codes: Vec<CipPropertySnapshot>,
    atom_neighbor_orders: Vec<CipPropertySnapshot>,
    atom_ranks: Vec<CipPropertySnapshot>,
    bond_codes: Vec<CipPropertySnapshot>,
    bond_neighbor_orders: Vec<CipPropertySnapshot>,
}

#[cfg(feature = "op-contracts")]
#[derive(Debug, Clone, PartialEq, Eq)]
struct CipPropertySnapshot {
    value: Option<String>,
    computed: bool,
}

#[cfg(feature = "op-contracts")]
impl CipPropertySnapshot {
    fn new(value: Option<&str>, computed: bool) -> Self {
        Self {
            value: value.map(str::to_owned),
            computed,
        }
    }

    fn after_clear_computed(&self) -> Self {
        if self.computed {
            Self {
                value: None,
                computed: false,
            }
        } else {
            self.clone()
        }
    }
}

#[cfg(feature = "op-contracts")]
impl CipStateSnapshot {
    fn from_molecule(molecule: &Molecule) -> Self {
        Self {
            computed: CipPropertySnapshot::new(
                molecule.prop("_CIPComputed"),
                molecule.is_prop_computed("_CIPComputed"),
            ),
            atom_codes: molecule
                .atoms()
                .iter()
                .map(|atom| {
                    CipPropertySnapshot::new(
                        atom.prop("_CIPCode"),
                        atom.is_prop_computed("_CIPCode"),
                    )
                })
                .collect(),
            atom_neighbor_orders: molecule
                .atoms()
                .iter()
                .map(|atom| {
                    CipPropertySnapshot::new(
                        atom.prop("_CIPNeighborOrder"),
                        atom.is_prop_computed("_CIPNeighborOrder"),
                    )
                })
                .collect(),
            atom_ranks: molecule
                .atoms()
                .iter()
                .map(|atom| {
                    CipPropertySnapshot::new(
                        atom.prop("_CIPRank"),
                        atom.is_prop_computed("_CIPRank"),
                    )
                })
                .collect(),
            bond_codes: molecule
                .bonds()
                .iter()
                .map(|bond| {
                    CipPropertySnapshot::new(
                        bond.prop("_CIPCode"),
                        bond.is_prop_computed("_CIPCode"),
                    )
                })
                .collect(),
            bond_neighbor_orders: molecule
                .bonds()
                .iter()
                .map(|bond| {
                    CipPropertySnapshot::new(
                        bond.prop("_CIPNeighborOrder"),
                        bond.is_prop_computed("_CIPNeighborOrder"),
                    )
                })
                .collect(),
        }
    }
}

#[cfg(feature = "op-contracts")]
impl ContractSourceSnapshot {
    fn from_molecule(molecule: &Molecule) -> Self {
        Self {
            atom_count: molecule.num_atoms(),
            bond_endpoints: molecule
                .bonds()
                .iter()
                .map(|bond| (bond.begin(), bond.end()))
                .collect(),
            cip_state: CipStateSnapshot::from_molecule(molecule),
        }
    }

    const fn num_atoms(&self) -> usize {
        self.atom_count
    }

    fn num_bonds(&self) -> usize {
        self.bond_endpoints.len()
    }

    fn bond_endpoint(&self, index: usize) -> Option<(AtomId, AtomId)> {
        self.bond_endpoints.get(index).copied()
    }

    fn cip_state(&self) -> CipStateSnapshot {
        self.cip_state.clone()
    }
}

#[cfg(feature = "op-contracts")]
enum ContractSource<'a> {
    Borrowed(&'a Molecule),
    Snapshot(ContractSourceSnapshot),
}

#[cfg(feature = "op-contracts")]
impl ContractSource<'_> {
    fn num_atoms(&self) -> usize {
        match self {
            Self::Borrowed(molecule) => molecule.num_atoms(),
            Self::Snapshot(snapshot) => snapshot.num_atoms(),
        }
    }

    fn num_bonds(&self) -> usize {
        match self {
            Self::Borrowed(molecule) => molecule.num_bonds(),
            Self::Snapshot(snapshot) => snapshot.num_bonds(),
        }
    }

    fn bond_endpoint(&self, index: usize) -> Option<(AtomId, AtomId)> {
        match self {
            Self::Borrowed(molecule) => molecule
                .bonds()
                .get(index)
                .map(|bond| (bond.begin(), bond.end())),
            Self::Snapshot(snapshot) => snapshot.bond_endpoint(index),
        }
    }

    fn cip_state(&self) -> CipStateSnapshot {
        match self {
            Self::Borrowed(molecule) => CipStateSnapshot::from_molecule(molecule),
            Self::Snapshot(snapshot) => snapshot.cip_state(),
        }
    }
}

pub struct OpParts<'a> {
    spec: &'static MoleculeOpSpec,
    #[cfg(feature = "op-contracts")]
    contract_source: ContractSource<'a>,
    working: Molecule,
    in_place_target: Option<&'a mut Molecule>,
    topology_mapping: Option<TopologyMapping>,
    #[cfg(feature = "op-contracts")]
    topology_lifecycle: BlockLifecycle,
    #[cfg(feature = "op-contracts")]
    coordinates_lifecycle: BlockLifecycle,
    #[cfg(feature = "op-contracts")]
    properties_lifecycle: BlockLifecycle,

    #[cfg(feature = "op-contracts")]
    trace: OperationTrace,
}

impl<'a> OpParts<'a> {
    pub(crate) fn new(
        source: &'a Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        validate_semantic_preconditions(source, spec)?;
        Ok(Self {
            spec,
            #[cfg(feature = "op-contracts")]
            contract_source: ContractSource::Borrowed(source),
            working: source.clone(),
            in_place_target: None,
            topology_mapping: None,
            #[cfg(feature = "op-contracts")]
            topology_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            coordinates_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            properties_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            trace: OperationTrace {
                touched_blocks: BlockSet::NONE,
                claimed_write_blocks: BlockSet::NONE,
                recorded_topology_edit: TopologyEditKind::None,
                remapped_blocks: BlockSet::NONE,
                preserved_cache: DerivedState::NONE,
                cleared_cache: DerivedState::NONE,
                updated_cache: DerivedState::NONE,
            },
        })
    }

    pub(crate) fn new_in_place(
        target: &'a mut Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        validate_semantic_preconditions(target, spec)?;
        #[cfg(feature = "op-contracts")]
        let contract_source =
            ContractSource::Snapshot(ContractSourceSnapshot::from_molecule(target));
        let working = std::mem::take(target);
        Ok(Self {
            spec,
            #[cfg(feature = "op-contracts")]
            contract_source,
            working,
            in_place_target: Some(target),
            topology_mapping: None,
            #[cfg(feature = "op-contracts")]
            topology_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            coordinates_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            properties_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            trace: OperationTrace {
                touched_blocks: BlockSet::NONE,
                claimed_write_blocks: BlockSet::NONE,
                recorded_topology_edit: TopologyEditKind::None,
                remapped_blocks: BlockSet::NONE,
                preserved_cache: DerivedState::NONE,
                cleared_cache: DerivedState::NONE,
                updated_cache: DerivedState::NONE,
            },
        })
    }

    // Agent guardrail:
    // OpParts is wrapper-owned state migration and contract-recording machinery.
    // Do not add chemistry, perception, sanitization, or operation-specific
    // behavior here. Operation impl_fn bodies compute behavior themselves or
    // call domain modules, then use OpParts only to read the working molecule,
    // apply state changes, and record contract effects.
    pub(crate) fn begin_topology_read(&self) -> Result<MoleculeReadParts<'_>, OperationError> {
        self.validate_access_spec()?;
        #[cfg(feature = "op-contracts")]
        {
            if !self.spec.access.read().contains(BlockSet::TOPOLOGY)
                || self.spec.access.write().contains(BlockSet::TOPOLOGY)
            {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation attempted to read topology outside its registry read access",
                });
            }
        }
        Ok(MoleculeReadParts::from_molecule_with_access(
            &self.working,
            crate::read_parts::MoleculeReadAccess::TOPOLOGY,
        ))
    }

    pub(crate) fn begin_operation_read(&self) -> Result<MoleculeReadParts<'_>, OperationError> {
        self.validate_access_spec()?;
        #[cfg(feature = "op-contracts")]
        {
            let required = BlockSet::TOPOLOGY
                .union(BlockSet::COORDINATES)
                .union(BlockSet::PROPERTIES)
                .union(BlockSet::DERIVED_CACHE);
            if !self.spec.access.can_read(required) {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation attempted a full read outside its registry access",
                });
            }
        }
        Ok(MoleculeReadParts::from_molecule_with_access(
            &self.working,
            crate::read_parts::MoleculeReadAccess::ALL,
        ))
    }

    pub(crate) fn with_topology_read_parts<R>(
        &self,
        topology: TopologyBlock,
        f: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let view = self.read_parts_for_topology(topology)?;
        f(MoleculeReadParts::from_molecule_with_access(
            &view,
            self.read_access(),
        ))
    }

    pub(crate) fn with_block_read_parts<R>(
        &self,
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        f: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let view = self.read_parts_for_blocks(topology, coordinates, properties)?;
        f(MoleculeReadParts::from_molecule_with_access(
            &view,
            self.read_access(),
        ))
    }

    pub(crate) fn with_optional_block_read_parts<R>(
        &self,
        topology: TopologyBlock,
        coordinates: Option<&CoordinateBlock>,
        properties: Option<&MoleculeProperties>,
        f: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let view = self.read_parts_for_optional_blocks(topology, coordinates, properties)?;
        f(MoleculeReadParts::from_molecule_with_access(
            &view,
            self.read_access(),
        ))
    }

    pub(crate) fn with_borrowed_optional_block_read_parts<R>(
        &self,
        topology: &TopologyBlock,
        coordinates: Option<&CoordinateBlock>,
        properties: Option<&MoleculeProperties>,
        f: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        self.validate_access_spec()?;
        let read = MoleculeReadParts::from_blocks_with_access(
            topology,
            coordinates.unwrap_or_else(|| self.working.coordinate_block()),
            properties.unwrap_or_else(|| self.working.properties()),
            self.working.derived_cache(),
            self.working.capabilities(),
            self.read_access(),
        );
        f(read)
    }

    pub(crate) fn with_coordinate_update_read_parts<R>(
        &self,
        f: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        self.validate_access_spec()?;
        #[cfg(feature = "op-contracts")]
        if !self.spec.access.can_write(BlockSet::COORDINATES) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation attempted to read a local coordinate block without coordinate write access",
            });
        }
        let read = MoleculeReadParts::from_blocks_with_access(
            self.working.topology_block(),
            self.working.coordinate_block(),
            self.working.properties(),
            self.working.derived_cache(),
            self.working.capabilities(),
            self.declared_read_access(),
        );
        f(read)
    }

    fn read_access(&self) -> crate::read_parts::MoleculeReadAccess {
        let blocks = self.spec.access.read().union(self.spec.access.write());
        Self::read_access_for_blocks(blocks)
    }

    fn declared_read_access(&self) -> crate::read_parts::MoleculeReadAccess {
        Self::read_access_for_blocks(self.spec.access.read())
    }

    fn read_access_for_blocks(blocks: BlockSet) -> crate::read_parts::MoleculeReadAccess {
        let mut access = crate::read_parts::MoleculeReadAccess::NONE;
        if blocks.contains(BlockSet::TOPOLOGY) {
            access = access.union(crate::read_parts::MoleculeReadAccess::TOPOLOGY);
        }
        if blocks.contains(BlockSet::COORDINATES) {
            access = access.union(crate::read_parts::MoleculeReadAccess::COORDINATES);
        }
        if blocks.contains(BlockSet::PROPERTIES) {
            access = access.union(crate::read_parts::MoleculeReadAccess::PROPERTIES);
        }
        if blocks.contains(BlockSet::DERIVED_CACHE) {
            access = access.union(crate::read_parts::MoleculeReadAccess::DERIVED_CACHE);
        }
        access
    }

    fn read_parts_for_topology(&self, topology: TopologyBlock) -> Result<Molecule, OperationError> {
        self.read_parts_for_blocks(
            topology,
            self.working.coordinate_block().clone(),
            self.working.properties().clone(),
        )
    }

    fn read_parts_for_blocks(
        &self,
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Molecule, OperationError> {
        Molecule::from_operation_blocks(
            topology,
            coordinates,
            properties,
            self.working.derived_cache().clone(),
            self.working.capabilities_block(),
        )
        .map_err(|failure| OperationError::InvariantViolation {
            operation: self.spec,
            failure,
        })
    }

    fn read_parts_for_optional_blocks(
        &self,
        topology: TopologyBlock,
        coordinates: Option<&CoordinateBlock>,
        properties: Option<&MoleculeProperties>,
    ) -> Result<Molecule, OperationError> {
        self.read_parts_for_blocks(
            topology,
            coordinates
                .cloned()
                .unwrap_or_else(|| self.working.coordinate_block().clone()),
            properties
                .cloned()
                .unwrap_or_else(|| self.working.properties().clone()),
        )
    }

    pub(crate) fn begin_topology_mut(&mut self) -> Result<TopologyBlock, OperationError> {
        self.begin_block_mut(BlockSet::TOPOLOGY)?;
        #[cfg(feature = "op-contracts")]
        {
            self.topology_lifecycle = BlockLifecycle::Begun;
        }
        Ok(if self.in_place_target.is_some() {
            self.working.take_topology_block_or_clone_from_operation()
        } else {
            self.working.topology_block().clone()
        })
    }

    pub(crate) fn commit_topology(
        &mut self,
        mut topology: TopologyBlock,
    ) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            if self.topology_lifecycle != BlockLifecycle::Begun {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "topology block was not begun before commit",
                });
            }
        }
        topology.adjacency =
            crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
        self.record_mutation(BlockSet::TOPOLOGY);
        self.working.replace_topology_block_from_operation(topology);
        #[cfg(feature = "op-contracts")]
        {
            self.topology_lifecycle = BlockLifecycle::Committed;
        }
        Ok(())
    }

    pub(crate) fn begin_coordinates_mut(&mut self) -> Result<CoordinateBlock, OperationError> {
        self.begin_block_mut(BlockSet::COORDINATES)?;
        #[cfg(feature = "op-contracts")]
        {
            self.coordinates_lifecycle = BlockLifecycle::Begun;
        }
        Ok(if self.in_place_target.is_some() {
            self.working.take_coordinate_block_or_clone()
        } else {
            self.working.coordinate_block().clone()
        })
    }

    pub(crate) fn commit_coordinates(
        &mut self,
        coordinates: CoordinateBlock,
    ) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            if self.coordinates_lifecycle != BlockLifecycle::Begun {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "coordinate block was not begun before commit",
                });
            }
        }
        self.record_mutation(BlockSet::COORDINATES);
        self.working.replace_coordinate_block(coordinates);
        #[cfg(feature = "op-contracts")]
        {
            self.coordinates_lifecycle = BlockLifecycle::Committed;
        }
        Ok(())
    }

    pub(crate) fn begin_properties_mut(&mut self) -> Result<MoleculeProperties, OperationError> {
        self.begin_block_mut(BlockSet::PROPERTIES)?;
        #[cfg(feature = "op-contracts")]
        {
            self.properties_lifecycle = BlockLifecycle::Begun;
        }
        Ok(if self.in_place_target.is_some() {
            self.working.take_properties_or_clone()
        } else {
            self.working.properties().clone()
        })
    }

    pub(crate) fn commit_properties(
        &mut self,
        properties: MoleculeProperties,
    ) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            if self.properties_lifecycle != BlockLifecycle::Begun {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "properties block was not begun before commit",
                });
            }
        }
        self.record_mutation(BlockSet::PROPERTIES);
        self.working.replace_properties(properties);
        #[cfg(feature = "op-contracts")]
        {
            self.properties_lifecycle = BlockLifecycle::Committed;
        }
        Ok(())
    }

    pub(crate) fn with_topology_mut<R>(
        &mut self,
        mutate: impl FnOnce(&mut Self, &mut TopologyBlock) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let mut topology = self.begin_topology_mut()?;
        let result = mutate(self, &mut topology);
        self.commit_topology(topology)?;
        result
    }

    pub(crate) fn with_coordinates_mut<R>(
        &mut self,
        mutate: impl FnOnce(&mut Self, &mut CoordinateBlock) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let mut coordinates = self.begin_coordinates_mut()?;
        let result = mutate(self, &mut coordinates);
        self.commit_coordinates(coordinates)?;
        result
    }

    pub(crate) fn with_topology_and_properties_mut<R>(
        &mut self,
        mutate: impl FnOnce(
            &mut Self,
            &mut TopologyBlock,
            &mut MoleculeProperties,
        ) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let mut topology = self.begin_topology_mut()?;
        let mut properties = match self.begin_properties_mut() {
            Ok(properties) => properties,
            Err(error) => {
                self.commit_topology(topology)?;
                return Err(error);
            }
        };
        let result = mutate(self, &mut topology, &mut properties);
        self.commit_topology(topology)?;
        self.commit_properties(properties)?;
        result
    }

    pub(crate) fn with_topology_coordinates_properties_mut<R>(
        &mut self,
        mutate: impl FnOnce(
            &mut Self,
            &mut TopologyBlock,
            &mut CoordinateBlock,
            &mut MoleculeProperties,
        ) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        let mut topology = self.begin_topology_mut()?;
        let mut coordinates = match self.begin_coordinates_mut() {
            Ok(coordinates) => coordinates,
            Err(error) => {
                self.commit_topology(topology)?;
                return Err(error);
            }
        };
        let mut properties = match self.begin_properties_mut() {
            Ok(properties) => properties,
            Err(error) => {
                self.commit_topology(topology)?;
                self.commit_coordinates(coordinates)?;
                return Err(error);
            }
        };
        let result = mutate(self, &mut topology, &mut coordinates, &mut properties);
        self.commit_topology(topology)?;
        self.commit_coordinates(coordinates)?;
        self.commit_properties(properties)?;
        result
    }

    pub(crate) fn record_topology_edit(
        &mut self,
        kind: TopologyEditKind,
    ) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            if kind == TopologyEditKind::Local
                && matches!(
                    self.spec.topology_edit,
                    TopologyEditKind::Appending
                        | TopologyEditKind::Compacting
                        | TopologyEditKind::Renumbering
                        | TopologyEditKind::Merge
                )
            {
                return Ok(());
            }
            if self.spec.topology_edit != kind {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "recorded topology edit does not match registry declaration",
                });
            }
            self.trace.recorded_topology_edit = kind;
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = kind;
        }
        Ok(())
    }

    pub(crate) fn record_topology_mapping(&mut self, mapping: TopologyMapping) {
        self.topology_mapping = Some(mapping);
        self.record_remapped(self.spec.auto_remap);
    }

    /// Check that `state` is in `recompute` or the guarded
    /// `operation_defined` category (write permission).
    /// Panics with a clear message on violation — this is a programming error.
    #[cfg(feature = "op-contracts")]
    fn check_cache_write_permission(&self, state: DerivedState) {
        let effects = self.spec.derived_effects;
        let allowed = effects.recompute().union(effects.operation_defined());
        if !allowed.contains(state) {
            panic!(
                "cache write permission violation: operation `{}` attempted to write \
                 derived state `{:?}` but only has `recompute({:?})` and \
                 `operation_defined({:?})` \
                 permissions",
                self.spec.method,
                state,
                effects.recompute(),
                effects.operation_defined(),
            );
        }
    }

    /// Check that every bit set in `states` is in `invalidate | recompute |
    /// operation_defined` (clear permission).
    /// Panics with a clear message on violation.
    #[cfg(feature = "op-contracts")]
    fn check_cache_clear_permission(&self, states: DerivedState) {
        let effects = self.spec.derived_effects;
        let allowed = effects
            .invalidate()
            .union(effects.recompute())
            .union(effects.operation_defined());
        let forbidden = states.bits() & !allowed.bits();
        if forbidden != 0 {
            panic!(
                "cache clear permission violation: operation `{}` attempted to clear \
                 derived state bits `{:#010b}` but only has `invalidate({:?})`, \
                 `recompute({:?})`, and `operation_defined({:?})` \
                 permissions",
                self.spec.method,
                forbidden,
                effects.invalidate(),
                effects.recompute(),
                effects.operation_defined(),
            );
        }
    }

    pub(crate) fn set_rings_cache(&mut self, rings: crate::RingInfo) {
        #[cfg(feature = "op-contracts")]
        self.check_cache_write_permission(DerivedState::RINGS);
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().rings = Some(rings);
        self.record_updated_cache(DerivedState::RINGS);
    }

    pub(crate) fn set_ring_families_cache(&mut self, ring_families: crate::RingInfo) {
        #[cfg(feature = "op-contracts")]
        self.check_cache_write_permission(DerivedState::RING_FAMILIES);
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().ring_families = Some(ring_families);
        self.record_updated_cache(DerivedState::RING_FAMILIES);
    }

    pub(crate) fn set_valence_cache(&mut self, valence: crate::ValenceAssignment) {
        #[cfg(feature = "op-contracts")]
        self.check_cache_write_permission(DerivedState::VALENCE);
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().valence = Some(valence);
        self.record_updated_cache(DerivedState::VALENCE);
    }

    pub(crate) fn mark_aromaticity_valid(&mut self) {
        #[cfg(feature = "op-contracts")]
        self.check_cache_write_permission(DerivedState::AROMATICITY);
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().aromaticity_valid = true;
        self.record_updated_cache(DerivedState::AROMATICITY);
    }

    #[allow(dead_code)]
    pub(crate) fn mark_stereo_handled(&mut self) {
        #[cfg(feature = "op-contracts")]
        self.check_cache_write_permission(DerivedState::STEREO);
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().stereo_valid = true;
        self.record_updated_cache(DerivedState::STEREO);
    }

    pub(crate) fn clear_cache(&mut self, states: DerivedState) {
        #[cfg(feature = "op-contracts")]
        {
            self.check_cache_clear_permission(states);
            self.trace.cleared_cache |= states;
        }
        if states.touches_cache() {
            self.record_mutation(BlockSet::DERIVED_CACHE);
            self.working.derived_cache_mut().invalidate(states);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = states;
        }
    }

    pub(crate) fn clear_computed_properties(&mut self) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.clear_computed_property_cache();
    }

    pub(crate) fn prove_preserved(
        &mut self,
        states: DerivedState,
        proof: PreservationProof,
    ) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            if !self.spec.derived_effects.preserve().contains(states) {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation attempted to prove preservation for undeclared derived states",
                });
            }
            match proof {
                PreservationProof::LeafAtomAppend => {
                    self.validate_leaf_atom_append_preservation()?
                }
            }
            self.trace.preserved_cache |= states;
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = states;
            let _ = proof;
        }
        Ok(())
    }

    pub(crate) fn finish(self) -> Result<Molecule, OperationError> {
        debug_assert!(self.in_place_target.is_none());
        #[cfg(feature = "op-contracts")]
        {
            let this = self;
            this.validate_contract()?;
            enforce_molecule_invariants(&this.working).map_err(|failure| {
                OperationError::InvariantViolation {
                    operation: this.spec,
                    failure,
                }
            })?;
            Ok(this.working)
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            enforce_molecule_invariants(&self.working).map_err(|failure| {
                OperationError::InvariantViolation {
                    operation: self.spec,
                    failure,
                }
            })?;
            Ok(self.working)
        }
    }

    pub(crate) fn abort_in_place(mut self) {
        let Some(target) = self.in_place_target.take() else {
            return;
        };
        *target = self.working;
    }

    pub(crate) fn finish_in_place(self) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            let mut this = self;
            let validation = this.validate_contract().and_then(|()| {
                enforce_molecule_invariants(&this.working).map_err(|failure| {
                    OperationError::InvariantViolation {
                        operation: this.spec,
                        failure,
                    }
                })
            });
            let target = this
                .in_place_target
                .take()
                .ok_or(OperationError::InvalidInput {
                    operation: this.spec,
                    message: "in-place operation was finished without an in-place target",
                })?;
            *target = this.working;
            validation
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let mut this = self;
            let validation = enforce_molecule_invariants(&this.working).map_err(|failure| {
                OperationError::InvariantViolation {
                    operation: this.spec,
                    failure,
                }
            });
            let target = this
                .in_place_target
                .take()
                .ok_or(OperationError::InvalidInput {
                    operation: this.spec,
                    message: "in-place operation was finished without an in-place target",
                })?;
            *target = this.working;
            validation
        }
    }

    fn record_mutation(&mut self, block: BlockSet) {
        #[cfg(feature = "op-contracts")]
        {
            assert!(
                self.spec.access.can_write(block) && self.spec.may_mutate.contains(block),
                "operation `{}` attempted to mutate a block outside its registry permissions",
                self.spec.method
            );
            self.trace.touched_blocks = self.trace.touched_blocks.union(block);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = block;
        }
    }

    #[cfg(feature = "op-contracts")]
    fn validate_access_spec(&self) -> Result<(), OperationError> {
        if self.spec.access.has_overlapping_read_write() {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation access declares the same block as both read and write",
            });
        }
        if self.spec.access.write() != self.spec.may_mutate {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation access write set must match may_mutate",
            });
        }
        let operation_defined = self.spec.derived_effects.operation_defined();
        if operation_defined != DerivedState::NONE {
            let approved_operation = matches!(
                self.spec.method,
                "without_hydrogens_with_sanitize" | "without_hydrogens_with_params"
            );
            if !approved_operation || operation_defined != DerivedState::VALENCE {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation-defined derived effects are currently permitted only for valence in the hydrogen-removal operation family",
                });
            }
            if !self.spec.access.can_write(BlockSet::DERIVED_CACHE) {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation-defined derived effects require derived-cache write access",
                });
            }
        }
        match self.spec.cip_state {
            CipStatePolicy::Assign if self.spec.method != "with_cip_labels_with_options" => {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "modern CIP assignment state is owned only by with_cip_labels_with_options",
                });
            }
            CipStatePolicy::ClearComputed
                if !self.spec.access.can_write(BlockSet::TOPOLOGY)
                    || !self.spec.access.can_write(BlockSet::PROPERTIES) =>
            {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "clearing computed CIP state requires topology and properties write access",
                });
            }
            _ => {}
        }
        Ok(())
    }

    #[cfg(not(feature = "op-contracts"))]
    fn validate_access_spec(&self) -> Result<(), OperationError> {
        Ok(())
    }

    fn begin_block_mut(&mut self, block: BlockSet) -> Result<(), OperationError> {
        self.validate_access_spec()?;
        #[cfg(feature = "op-contracts")]
        {
            if !self.spec.access.can_write(block) {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation attempted to write a block outside its registry access",
                });
            }
            let lifecycle = if block == BlockSet::TOPOLOGY {
                self.topology_lifecycle
            } else if block == BlockSet::COORDINATES {
                self.coordinates_lifecycle
            } else if block == BlockSet::PROPERTIES {
                self.properties_lifecycle
            } else {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation attempted to begin an unknown block",
                });
            };
            if lifecycle != BlockLifecycle::Available {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "operation attempted to begin the same writable block twice",
                });
            }
            self.trace.claimed_write_blocks = self.trace.claimed_write_blocks.union(block);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = block;
        }
        Ok(())
    }

    fn record_updated_cache(&mut self, state: DerivedState) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.updated_cache = self.trace.updated_cache.union(state);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = state;
        }
    }

    #[cfg(feature = "op-contracts")]
    fn validate_leaf_atom_append_preservation(&self) -> Result<(), OperationError> {
        if self.spec.topology_edit != TopologyEditKind::Appending {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "leaf-append preservation proof requires an appending topology operation",
            });
        }
        let Some(mapping) = &self.topology_mapping else {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "leaf-append preservation proof requires a topology mapping",
            });
        };
        let old_atom_count = self.contract_source.num_atoms();
        let old_bond_count = self.contract_source.num_bonds();
        if mapping.atoms().old_to_new().len() != old_atom_count
            || mapping.bonds().old_to_new().len() != old_bond_count
            || mapping.atoms().new_to_old().len() != self.working.num_atoms()
            || mapping.bonds().new_to_old().len() != self.working.num_bonds()
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "leaf-append preservation proof has inconsistent mapping dimensions",
            });
        }
        for (old_idx, mapped) in mapping.atoms().old_to_new().iter().enumerate() {
            if *mapped != Some(AtomId::new(old_idx)) {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "leaf-append preservation proof requires identity mapping for old atoms",
                });
            }
        }
        for (old_idx, mapped) in mapping.bonds().old_to_new().iter().enumerate() {
            if *mapped != Some(crate::BondId::new(old_idx)) {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "leaf-append preservation proof requires identity mapping for old bonds",
                });
            }
        }
        for old_idx in 0..old_bond_count {
            let Some((before_begin, before_end)) = self.contract_source.bond_endpoint(old_idx)
            else {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "leaf-append preservation proof found missing source bond endpoint",
                });
            };
            let after = &self.working.bonds()[old_idx];
            if before_begin != after.begin() || before_end != after.end() {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "leaf-append preservation proof detected changed old bond endpoints",
                });
            }
        }

        let mut appended_degrees =
            vec![0usize; self.working.num_atoms().saturating_sub(old_atom_count)];
        for bond in &self.working.bonds()[old_bond_count..] {
            let begin_old = bond.begin().index() < old_atom_count;
            let end_old = bond.end().index() < old_atom_count;
            if begin_old == end_old {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "leaf-append preservation proof requires every appended bond to connect one old atom and one appended atom",
                });
            }
            let appended_idx = if begin_old {
                bond.end().index() - old_atom_count
            } else {
                bond.begin().index() - old_atom_count
            };
            if appended_idx >= appended_degrees.len() {
                return Err(OperationError::InvalidInput {
                    operation: self.spec,
                    message: "leaf-append preservation proof found appended bond referencing an out-of-range atom",
                });
            }
            appended_degrees[appended_idx] += 1;
        }
        if appended_degrees.iter().any(|degree| *degree != 1) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "leaf-append preservation proof requires every appended atom to be a degree-one leaf",
            });
        }
        Ok(())
    }

    #[cfg(feature = "op-contracts")]
    fn validate_contract(&self) -> Result<(), OperationError> {
        self.validate_access_spec()?;
        if self.topology_lifecycle == BlockLifecycle::Begun
            || self.coordinates_lifecycle == BlockLifecycle::Begun
            || self.properties_lifecycle == BlockLifecycle::Begun
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation finished while a writable block was still checked out",
            });
        }
        if !self
            .trace
            .touched_blocks
            .contains(self.trace.claimed_write_blocks)
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation did not commit every claimed writable block",
            });
        }
        let effects = self.spec.derived_effects;
        let recompute_ds = effects.recompute();
        let operation_defined_ds = effects.operation_defined();
        if recompute_ds.intersects(effects.preserve())
            || recompute_ds.intersects(effects.invalidate())
            || recompute_ds.intersects(operation_defined_ds)
            || effects.preserve().intersects(effects.invalidate())
            || effects.preserve().intersects(operation_defined_ds)
            || effects.invalidate().intersects(operation_defined_ds)
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation derived_effects contains overlapping effect categories",
            });
        }

        let requires_effect_handling = self.trace.touched_blocks != BlockSet::NONE;
        let updated_or_cleared = self.trace.cleared_cache | self.trace.updated_cache;
        if requires_effect_handling && !updated_or_cleared.contains(self.spec.needs_update()) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not clear or update every required cache state",
            });
        }

        if requires_effect_handling
            && !self
                .trace
                .preserved_cache
                .contains(self.spec.derived_effects.preserve())
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not prove every declared preserved derived state",
            });
        }

        if self.spec.requires_mapping == MappingRequirement::Required
            && self.topology_mapping.is_none()
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "strong topology operation did not record a topology mapping",
            });
        }

        if !self.trace.remapped_blocks.contains(self.spec.auto_remap) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation did not remap every registry-required block",
            });
        }

        if self.trace.touched_blocks.contains(BlockSet::TOPOLOGY)
            && self.spec.topology_edit != TopologyEditKind::None
            && self.trace.recorded_topology_edit == TopologyEditKind::None
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation touched topology without recording the registry topology edit",
            });
        }

        self.validate_cip_state_lifecycle()?;

        Ok(())
    }

    #[cfg(feature = "op-contracts")]
    fn validate_cip_state_lifecycle(&self) -> Result<(), OperationError> {
        let before = self.contract_source.cip_state();
        let after = CipStateSnapshot::from_molecule(&self.working);
        match self.spec.cip_state {
            CipStatePolicy::Assign => Ok(()),
            CipStatePolicy::Preserve if before == after => Ok(()),
            CipStatePolicy::Preserve => Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation declared CIP state preservation but changed observable CIP state",
            }),
            CipStatePolicy::ClearComputed => {
                if after.computed.computed
                    || after.atom_codes.iter().any(|state| state.computed)
                    || after
                        .atom_neighbor_orders
                        .iter()
                        .any(|state| state.computed)
                    || after.atom_ranks.iter().any(|state| state.computed)
                    || after.bond_codes.iter().any(|state| state.computed)
                    || after
                        .bond_neighbor_orders
                        .iter()
                        .any(|state| state.computed)
                {
                    return Err(OperationError::InvalidInput {
                        operation: self.spec,
                        message: "operation declared computed CIP clearing but retained computed CIP state",
                    });
                }
                if after.computed != before.computed.after_clear_computed() {
                    return Err(OperationError::InvalidInput {
                        operation: self.spec,
                        message: "operation did not reproduce molecule computed-property clearing for _CIPComputed",
                    });
                }
                self.validate_mapped_cip_properties(&before, &after)
            }
        }
    }

    #[cfg(feature = "op-contracts")]
    fn validate_mapped_cip_properties(
        &self,
        before: &CipStateSnapshot,
        after: &CipStateSnapshot,
    ) -> Result<(), OperationError> {
        for old_index in 0..before.atom_codes.len() {
            let new_index = match &self.topology_mapping {
                Some(mapping) => mapping
                    .atoms()
                    .old_to_new()
                    .get(old_index)
                    .copied()
                    .flatten()
                    .map(AtomId::index),
                None => (old_index < after.atom_codes.len()).then_some(old_index),
            };
            if let Some(new_index) = new_index {
                for (before_values, after_values, property) in [
                    (&before.atom_codes, &after.atom_codes, "_CIPCode"),
                    (
                        &before.atom_neighbor_orders,
                        &after.atom_neighbor_orders,
                        "_CIPNeighborOrder",
                    ),
                    (&before.atom_ranks, &after.atom_ranks, "_CIPRank"),
                ] {
                    if after_values.get(new_index)
                        != before_values
                            .get(old_index)
                            .map(|state| state.after_clear_computed())
                            .as_ref()
                    {
                        return Err(OperationError::InvalidInput {
                            operation: self.spec,
                            message: match property {
                                "_CIPCode" => {
                                    "operation did not reproduce atom computed-property clearing for _CIPCode"
                                }
                                "_CIPNeighborOrder" => {
                                    "operation did not reproduce atom computed-property clearing for _CIPNeighborOrder"
                                }
                                _ => {
                                    "operation did not reproduce atom computed-property clearing for _CIPRank"
                                }
                            },
                        });
                    }
                }
            }
        }
        for old_index in 0..before.bond_codes.len() {
            let new_index = match &self.topology_mapping {
                Some(mapping) => mapping
                    .bonds()
                    .old_to_new()
                    .get(old_index)
                    .copied()
                    .flatten()
                    .map(crate::BondId::index),
                None => (old_index < after.bond_codes.len()).then_some(old_index),
            };
            if let Some(new_index) = new_index {
                for (before_values, after_values, property) in [
                    (&before.bond_codes, &after.bond_codes, "_CIPCode"),
                    (
                        &before.bond_neighbor_orders,
                        &after.bond_neighbor_orders,
                        "_CIPNeighborOrder",
                    ),
                ] {
                    if after_values.get(new_index)
                        != before_values
                            .get(old_index)
                            .map(|state| state.after_clear_computed())
                            .as_ref()
                    {
                        return Err(OperationError::InvalidInput {
                            operation: self.spec,
                            message: match property {
                                "_CIPCode" => {
                                    "operation did not reproduce bond computed-property clearing for _CIPCode"
                                }
                                _ => {
                                    "operation did not reproduce bond computed-property clearing for _CIPNeighborOrder"
                                }
                            },
                        });
                    }
                }
            }
        }
        Ok(())
    }

    fn record_remapped(&mut self, block: BlockSet) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.remapped_blocks = self.trace.remapped_blocks.union(block);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = block;
        }
    }
}

fn validate_semantic_preconditions(
    molecule: &Molecule,
    spec: &'static MoleculeOpSpec,
) -> Result<(), OperationError> {
    let preconditions = spec.semantic_preconditions;
    if preconditions.contains(SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY)
        && molecule.num_atoms() != 0
        && molecule.topology_trust() != TopologyTrust::TrustedGraph
    {
        return Err(OperationError::Precondition {
            operation: spec,
            requirement: SemanticPrecondition::TrustedBondTopology,
            message: "operation requires a molecule with trusted bond topology",
        });
    }
    if preconditions.contains(SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED)
        && !hydrogen_ownership_is_represented(molecule)
    {
        return Err(OperationError::Precondition {
            operation: spec,
            requirement: SemanticPrecondition::HydrogenOwnershipRepresented,
            message: "explicit hydrogen atoms must be connected to their owning heavy atoms",
        });
    }
    Ok(())
}

fn hydrogen_ownership_is_represented(molecule: &Molecule) -> bool {
    if !molecule
        .atoms()
        .iter()
        .any(|atom| atom.atomic_number() != 1)
    {
        return true;
    }
    molecule
        .atoms()
        .iter()
        .filter(|atom| atom.atomic_number() == 1)
        .all(|atom| {
            molecule
                .topology_block()
                .adjacency
                .neighbors_of(atom.id().index())
                .iter()
                .any(|neighbor| {
                    molecule
                        .atom(AtomId::new(neighbor.atom_index))
                        .is_some_and(|neighbor_atom| neighbor_atom.atomic_number() != 1)
                })
        })
}
