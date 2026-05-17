//! Macro-controlled molecule operations.
//!
//! This module is the canonical operation-contract surface for molecular state
//! migration. Operation bodies must be written as no-argument functions annotated
//! with `#[mol_op_body(op_name, parts)]`; the attribute macro injects the only
//! allowed mutable capability object.
//!
//! Agent guardrail: operation bodies must not take `&mut Molecule`, must not
//! access `Molecule` storage directly, and must not update derived-state caches
//! or invalidation flags by hand. If a future change appears to require bypassing
//! those rules, the agent is not allowed to continue silently. It must stop and
//! confirm the design exception with the human author before editing code.
//!
//! RDKit marker convention is defined in `dev/source_reproduction_protocol.md`.

use std::fmt;

use cosmolkit_macros::{mol_op_body, molecule_ops};

use crate::{
    Atom, AtomId, Bond, DerivedState, InvariantError, Molecule, MoleculeProperties, SupportStatus,
    invariants::enforce_molecule_invariants,
    molecule::{CoordinateBlock, DerivedCacheBlock, TopologyBlock, TopologyMapping},
};

pub(crate) use crate::read_parts::MoleculeReadParts;

trait MoleculeReadAccess<'a>: Copy {
    fn atoms(self) -> &'a [Atom];
    fn bonds(self) -> &'a [Bond];
    fn atom(self, atom: AtomId) -> Option<&'a Atom>;
    fn num_atoms(self) -> usize;
    fn derived_cache(self) -> &'a DerivedCacheBlock;
    fn assign_valence_with_options(
        self,
        model: crate::ValenceModel,
        strict: bool,
    ) -> Result<crate::ValenceAssignment, crate::ValenceError>;
    fn rank_mol_atoms(self) -> Result<Vec<usize>, crate::KekulizeError>;
}

impl<'a> MoleculeReadAccess<'a> for MoleculeReadParts<'a> {
    fn atoms(self) -> &'a [Atom] {
        MoleculeReadParts::atoms(self)
    }

    fn bonds(self) -> &'a [Bond] {
        MoleculeReadParts::bonds(self)
    }

    fn atom(self, atom: AtomId) -> Option<&'a Atom> {
        MoleculeReadParts::atom(self, atom)
    }

    fn num_atoms(self) -> usize {
        MoleculeReadParts::num_atoms(self)
    }

    fn derived_cache(self) -> &'a DerivedCacheBlock {
        MoleculeReadParts::derived_cache(self)
    }

    fn assign_valence_with_options(
        self,
        model: crate::ValenceModel,
        strict: bool,
    ) -> Result<crate::ValenceAssignment, crate::ValenceError> {
        MoleculeReadParts::assign_valence_with_options(self, model, strict)
    }

    fn rank_mol_atoms(self) -> Result<Vec<usize>, crate::KekulizeError> {
        MoleculeReadParts::rank_mol_atoms(self)
    }
}

impl<'a> MoleculeReadAccess<'a> for &'a Molecule {
    fn atoms(self) -> &'a [Atom] {
        self.atoms()
    }

    fn bonds(self) -> &'a [Bond] {
        self.bonds()
    }

    fn atom(self, atom: AtomId) -> Option<&'a Atom> {
        self.atom(atom)
    }

    fn num_atoms(self) -> usize {
        self.num_atoms()
    }

    fn derived_cache(self) -> &'a DerivedCacheBlock {
        self.derived_cache()
    }

    fn assign_valence_with_options(
        self,
        model: crate::ValenceModel,
        strict: bool,
    ) -> Result<crate::ValenceAssignment, crate::ValenceError> {
        crate::assign_valence_with_options(self, model, strict)
    }

    fn rank_mol_atoms(self) -> Result<Vec<usize>, crate::KekulizeError> {
        crate::canon_rank::rank_mol_atoms(self)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoleculeOpKind {
    Strong,
    Weak,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopologyEditKind {
    None,
    Local,
    Compacting,
    Appending,
    Renumbering,
    Merge,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperationDomain {
    Topology,
    Coordinate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParityPolicy {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockSet(u8);

impl BlockSet {
    pub const NONE: Self = Self(0);
    pub const TOPOLOGY: Self = Self(1 << 0);
    pub const COORDINATES: Self = Self(1 << 1);
    pub const PROPERTIES: Self = Self(1 << 2);
    pub const DERIVED_CACHE: Self = Self(1 << 3);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }

    #[must_use]
    pub const fn intersects(self, other: Self) -> bool {
        (self.0 & other.0) != 0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockAccess {
    read: BlockSet,
    write: BlockSet,
}

impl BlockAccess {
    pub const NONE: Self = Self {
        read: BlockSet::NONE,
        write: BlockSet::NONE,
    };

    #[must_use]
    pub const fn new(read: BlockSet, write: BlockSet) -> Self {
        Self { read, write }
    }

    #[must_use]
    pub const fn read(self) -> BlockSet {
        self.read
    }

    #[must_use]
    pub const fn write(self) -> BlockSet {
        self.write
    }

    #[must_use]
    pub const fn can_read(self, block: BlockSet) -> bool {
        self.read.contains(block) || self.write.contains(block)
    }

    #[must_use]
    pub const fn can_write(self, block: BlockSet) -> bool {
        self.write.contains(block)
    }

    #[must_use]
    pub const fn has_overlapping_read_write(self) -> bool {
        self.read.intersects(self.write)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreservationProof {
    LeafAtomAppend,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedEffects {
    recompute: DerivedState,
    preserve: DerivedState,
    invalidate: DerivedState,
}

impl DerivedEffects {
    pub const NONE: Self = Self {
        recompute: DerivedState::NONE,
        preserve: DerivedState::NONE,
        invalidate: DerivedState::NONE,
    };

    #[must_use]
    pub const fn new(
        recompute: DerivedState,
        preserve: DerivedState,
        invalidate: DerivedState,
    ) -> Self {
        Self {
            recompute,
            preserve,
            invalidate,
        }
    }

    #[must_use]
    pub const fn recompute(self) -> DerivedState {
        self.recompute
    }

    #[must_use]
    pub const fn preserve(self) -> DerivedState {
        self.preserve
    }

    #[must_use]
    pub const fn invalidate(self) -> DerivedState {
        self.invalidate
    }

    #[must_use]
    pub const fn needs_update(self) -> DerivedState {
        self.invalidate.union(self.recompute)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MappingRequirement {
    None,
    Identity,
    Required,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MoleculeOpSpec {
    pub method: &'static str,
    pub impl_fn: &'static str,
    pub domain: OperationDomain,
    pub kind: MoleculeOpKind,
    pub topology_edit: TopologyEditKind,
    pub access: BlockAccess,
    pub may_mutate: BlockSet,
    pub auto_remap: BlockSet,
    pub derived_effects: DerivedEffects,
    pub requires_mapping: MappingRequirement,
    pub allows_noop: bool,
    pub support: SupportStatus,
    pub parity: ParityPolicy,
    pub io_roundtrip: bool,
}

impl fmt::Display for MoleculeOpSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.method)
    }
}

impl MoleculeOpSpec {
    #[must_use]
    pub const fn needs_update(self: &Self) -> DerivedState {
        self.derived_effects.needs_update()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OpOutcome {
    Changed,
    NoOp { reason: &'static str },
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum OperationError {
    #[error("{operation}: unsupported operation: {reason}")]
    Unsupported {
        operation: &'static MoleculeOpSpec,
        reason: &'static str,
    },
    #[error("{operation}: invalid input: {message}")]
    InvalidInput {
        operation: &'static MoleculeOpSpec,
        message: &'static str,
    },
    #[error("{operation}: chemistry error: {message}")]
    Chemistry {
        operation: &'static MoleculeOpSpec,
        message: &'static str,
    },
    #[error("{operation}: {source}")]
    UnsupportedFeature {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::UnsupportedFeatureError,
    },
    #[error("{operation}: valence error: {source}")]
    Valence {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::ValenceError,
    },
    #[error("{operation}: ring finding error: {source}")]
    RingFinding {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::RingFindingError,
    },
    #[error("{operation}: aromaticity error: {source}")]
    Aromaticity {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::AromaticityError,
    },
    #[error("{operation}: {source}")]
    Sanitize {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::SanitizeError,
    },
    #[error("{operation}: kekulize error: {source}")]
    Kekulize {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::KekulizeError,
    },
    #[error("{operation}: hydrogen removal error: {source}")]
    RemoveHydrogens {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::RemoveHydrogensError,
    },
    #[error("{operation}: hydrogen addition error: {source}")]
    AddHydrogens {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::AddHydrogensError,
    },
    #[error("{operation}: invariant violation: {failure}")]
    InvariantViolation {
        operation: &'static MoleculeOpSpec,
        failure: InvariantError,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvariantCheckSet(u32);

impl InvariantCheckSet {
    pub const NONE: Self = Self(0);
    pub const GRAPH_INDEX_VALIDITY: Self = Self(1 << 0);
    pub const COORDINATE_ROW_ALIGNMENT: Self = Self(1 << 1);
    pub const MAPPING_RECORDED_OR_EXPLICITLY_UNSUPPORTED: Self = Self(1 << 2);
    pub const STEREO_VALIDITY: Self = Self(1 << 3);
    pub const CACHE_INVALIDATION: Self = Self(1 << 4);
    pub const PROPERTY_POLICY: Self = Self(1 << 5);
    pub const SOURCE_UNCHANGED: Self = Self(1 << 6);
    pub const UNSUPPORTED_IS_EXPLICIT: Self = Self(1 << 7);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SupportMatrixEntry {
    pub feature: &'static crate::FeatureSpec,
    pub operation: Option<&'static MoleculeOpSpec>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OperationInvariantEntry {
    pub operation: &'static MoleculeOpSpec,
    pub profile: &'static str,
    pub required_checks: InvariantCheckSet,
}

impl OperationInvariantEntry {
    #[must_use]
    pub const fn for_operation(operation: &'static MoleculeOpSpec, profile: &'static str) -> Self {
        Self {
            operation,
            profile,
            required_checks: InvariantCheckSet::GRAPH_INDEX_VALIDITY
                .union(InvariantCheckSet::COORDINATE_ROW_ALIGNMENT)
                .union(InvariantCheckSet::MAPPING_RECORDED_OR_EXPLICITLY_UNSUPPORTED)
                .union(InvariantCheckSet::STEREO_VALIDITY)
                .union(InvariantCheckSet::CACHE_INVALIDATION)
                .union(InvariantCheckSet::PROPERTY_POLICY)
                .union(InvariantCheckSet::SOURCE_UNCHANGED)
                .union(InvariantCheckSet::UNSUPPORTED_IS_EXPLICIT),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParityMatrixEntry {
    pub operation: &'static MoleculeOpSpec,
    pub feature: &'static crate::FeatureSpec,
    pub profile: &'static str,
    pub rdkit_version: Option<&'static str>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperationTrace {
    touched_blocks: BlockSet,
    claimed_write_blocks: BlockSet,
    recorded_topology_edit: TopologyEditKind,
    remapped_blocks: BlockSet,
    preserved_cache: DerivedState,
    read_cache: DerivedState,
    cleared_cache: DerivedState,
    updated_cache: DerivedState,
    outcome: Option<OpOutcome>,
}

impl OperationTrace {
    #[must_use]
    pub const fn touched_blocks(&self) -> BlockSet {
        self.touched_blocks
    }

    #[must_use]
    pub const fn read_cache(&self) -> DerivedState {
        self.read_cache
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

    #[must_use]
    pub fn outcome(&self) -> Option<&OpOutcome> {
        self.outcome.as_ref()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BlockLifecycle {
    Available,
    Begun,
    Committed,
}

pub struct OpParts<'a> {
    spec: &'static MoleculeOpSpec,
    source: &'a Molecule,
    working: Molecule,
    topology_mapping: Option<TopologyMapping>,
    topology_lifecycle: BlockLifecycle,
    coordinates_lifecycle: BlockLifecycle,
    properties_lifecycle: BlockLifecycle,

    #[cfg(feature = "op-contracts")]
    trace: OperationTrace,
}

impl<'a> OpParts<'a> {
    pub(crate) fn new(source: &'a Molecule, spec: &'static MoleculeOpSpec) -> Self {
        Self {
            spec,
            source,
            working: source.clone(),
            topology_mapping: None,
            topology_lifecycle: BlockLifecycle::Available,
            coordinates_lifecycle: BlockLifecycle::Available,
            properties_lifecycle: BlockLifecycle::Available,
            #[cfg(feature = "op-contracts")]
            trace: OperationTrace {
                touched_blocks: BlockSet::NONE,
                claimed_write_blocks: BlockSet::NONE,
                recorded_topology_edit: TopologyEditKind::None,
                remapped_blocks: BlockSet::NONE,
                preserved_cache: DerivedState::NONE,
                read_cache: DerivedState::NONE,
                cleared_cache: DerivedState::NONE,
                updated_cache: DerivedState::NONE,
                outcome: None,
            },
        }
    }

    // Agent guardrail:
    // OpParts is wrapper-owned state migration and contract-recording machinery.
    // Do not add chemistry, perception, sanitization, or operation-specific
    // behavior here. Operation impl_fn bodies compute behavior themselves or
    // call domain modules, then use OpParts only to read the working molecule,
    // apply state changes, and record contract effects.
    pub(crate) fn begin_topology_read(&self) -> Result<MoleculeReadParts<'_>, OperationError> {
        self.validate_access_spec()?;
        if !self.spec.access.read().contains(BlockSet::TOPOLOGY)
            || self.spec.access.write().contains(BlockSet::TOPOLOGY)
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation attempted to read topology outside its registry read access",
            });
        }
        Ok(MoleculeReadParts::from_molecule(&self.working))
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
        self.topology_lifecycle = BlockLifecycle::Begun;
        Ok(self.working.topology_block().clone())
    }

    pub(crate) fn commit_topology(
        &mut self,
        mut topology: TopologyBlock,
    ) -> Result<(), OperationError> {
        if self.topology_lifecycle != BlockLifecycle::Begun {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "topology block was not begun before commit",
            });
        }
        topology.adjacency =
            crate::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
        self.record_mutation(BlockSet::TOPOLOGY);
        self.working.replace_topology_block(topology);
        self.topology_lifecycle = BlockLifecycle::Committed;
        Ok(())
    }

    pub(crate) fn begin_coordinates_mut(&mut self) -> Result<CoordinateBlock, OperationError> {
        self.begin_block_mut(BlockSet::COORDINATES)?;
        self.coordinates_lifecycle = BlockLifecycle::Begun;
        Ok(self.working.coordinate_block().clone())
    }

    pub(crate) fn commit_coordinates(
        &mut self,
        coordinates: CoordinateBlock,
    ) -> Result<(), OperationError> {
        if self.coordinates_lifecycle != BlockLifecycle::Begun {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "coordinate block was not begun before commit",
            });
        }
        self.record_mutation(BlockSet::COORDINATES);
        self.working.replace_coordinate_block(coordinates);
        self.coordinates_lifecycle = BlockLifecycle::Committed;
        Ok(())
    }

    pub(crate) fn begin_properties_mut(&mut self) -> Result<MoleculeProperties, OperationError> {
        self.begin_block_mut(BlockSet::PROPERTIES)?;
        self.properties_lifecycle = BlockLifecycle::Begun;
        Ok(self.working.properties().clone())
    }

    pub(crate) fn commit_properties(
        &mut self,
        properties: MoleculeProperties,
    ) -> Result<(), OperationError> {
        if self.properties_lifecycle != BlockLifecycle::Begun {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "properties block was not begun before commit",
            });
        }
        self.record_mutation(BlockSet::PROPERTIES);
        self.working.replace_properties(properties);
        self.properties_lifecycle = BlockLifecycle::Committed;
        Ok(())
    }

    pub(crate) fn record_topology_edit(
        &mut self,
        kind: TopologyEditKind,
    ) -> Result<(), OperationError> {
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
        #[cfg(feature = "op-contracts")]
        {
            self.trace.recorded_topology_edit = kind;
        }
        Ok(())
    }

    pub(crate) fn record_topology_mapping(&mut self, mapping: TopologyMapping) {
        self.topology_mapping = Some(mapping);
        self.record_remapped(self.spec.auto_remap);
    }

    /// Check that `state` is in `requires | recompute` (write permission).
    /// Check that `state` is in `recompute` (write permission).
    /// Panics with a clear message on violation — this is a programming error.
    #[cfg(feature = "op-contracts")]
    fn check_cache_write_permission(&self, state: DerivedState) {
        let effects = self.spec.derived_effects;
        let allowed = effects.recompute();
        if !allowed.contains(state) {
            panic!(
                "cache write permission violation: operation `{}` attempted to write \
                 derived state `{:?}` but only has `recompute({:?})` \
                 permissions",
                self.spec.method,
                state,
                effects.recompute(),
            );
        }
    }

    /// Check that every bit set in `states` is in `invalidate | recompute` (clear permission).
    /// Panics with a clear message on violation.
    #[cfg(feature = "op-contracts")]
    fn check_cache_clear_permission(&self, states: DerivedState) {
        let effects = self.spec.derived_effects;
        let allowed = effects.invalidate().union(effects.recompute());
        let forbidden = states.bits() & !allowed.bits();
        if forbidden != 0 {
            panic!(
                "cache clear permission violation: operation `{}` attempted to clear \
                 derived state bits `{:#010b}` but only has `invalidate({:?})` and `recompute({:?})` \
                 permissions",
                self.spec.method,
                forbidden,
                effects.invalidate(),
                effects.recompute(),
            );
        }
    }

    /// Check that `state` is in `preserve` (read permission).
    /// Panics with a clear message on violation.
    #[cfg(feature = "op-contracts")]
    fn check_cache_read_permission(&self, state: DerivedState) {
        let effects = self.spec.derived_effects;
        let allowed = effects.preserve();
        if !allowed.contains(state) {
            panic!(
                "cache read permission violation: operation `{}` attempted to read \
                 derived state `{:?}` but only has `preserve({:?})` \
                 permissions",
                self.spec.method,
                state,
                effects.preserve(),
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

    pub(crate) fn prove_preserved(
        &mut self,
        states: DerivedState,
        proof: PreservationProof,
    ) -> Result<(), OperationError> {
        if !self.spec.derived_effects.preserve().contains(states) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation attempted to prove preservation for undeclared derived states",
            });
        }
        match proof {
            PreservationProof::LeafAtomAppend => self.validate_leaf_atom_append_preservation()?,
        }
        #[cfg(feature = "op-contracts")]
        {
            self.trace.preserved_cache |= states;
        }
        Ok(())
    }

    pub(crate) fn finish(self, outcome: OpOutcome) -> Result<Molecule, OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            let mut this = self;
            this.trace.outcome = Some(outcome);
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
            let _ = outcome;
            enforce_molecule_invariants(&self.working).map_err(|failure| {
                OperationError::InvariantViolation {
                    operation: self.spec,
                    failure,
                }
            })?;
            Ok(self.working)
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
        Ok(())
    }

    fn begin_block_mut(&mut self, block: BlockSet) -> Result<(), OperationError> {
        self.validate_access_spec()?;
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
        #[cfg(feature = "op-contracts")]
        {
            self.trace.claimed_write_blocks = self.trace.claimed_write_blocks.union(block);
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
        let old_atom_count = self.source.num_atoms();
        let old_bond_count = self.source.num_bonds();
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
            let before = &self.source.bonds()[old_idx];
            let after = &self.working.bonds()[old_idx];
            if before.begin() != after.begin() || before.end() != after.end() {
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
        let effects = self.spec.derived_effects;
        let recompute_ds = effects.recompute();
        if recompute_ds.intersects(effects.preserve())
            || recompute_ds.intersects(effects.invalidate())
            || effects.preserve().intersects(effects.invalidate())
        {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation derived_effects contains overlapping effect categories",
            });
        }

        let updated_or_cleared = self.trace.cleared_cache | self.trace.updated_cache;
        if !updated_or_cleared.contains(self.spec.needs_update()) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not clear or update every required cache state",
            });
        }

        if !self
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

molecule_ops! {
    op with_hydrogens(params: crate::hydrogens::AddHsParams) {
        method: with_hydrogens_with_params,
        impl_fn: with_hydrogens_impl,
        default_method: with_hydrogens,
        default_args: [crate::hydrogens::AddHsParams::default()],
        domain: topology,
        kind: strong,
        topology_edit: appending,
        access: { read: [], write: [topology, coordinates, properties, derived_cache] },
        may_mutate: [topology, coordinates, properties, derived_cache],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [],
            preserve: [rings, ring_families],
            invalidate: [valence, aromaticity, stereo, drawing, fingerprint],
        },
        requires_mapping: required,
        allows_noop: false,
        feature: HYDROGENS_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "strong_topology_with_coordinates",
        parity_profile: "add_hs_default_rdkit",
    }

    op without_hydrogens(sanitize: bool) {
        method: without_hydrogens_with_sanitize,
        impl_fn: without_hydrogens_impl,
        default_method: without_hydrogens,
        default_args: [true],
        domain: topology,
        kind: strong,
        topology_edit: compacting,
        access: { read: [], write: [topology, coordinates, properties, derived_cache] },
        may_mutate: [topology, coordinates, properties, derived_cache],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [valence, aromaticity, rings],
            preserve: [],
            invalidate: [ring_families, stereo, drawing, fingerprint],
        },
        requires_mapping: required,
        allows_noop: true,
        feature: HYDROGENS_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "strong_topology_with_coordinates",
        parity_profile: "remove_hs_default_rdkit",
    }

    op without_hydrogens_with_params(params: crate::hydrogens::RemoveHsParams, sanitize: bool) {
        method: without_hydrogens_with_params,
        impl_fn: without_hydrogens_with_params_impl,
        domain: topology,
        kind: strong,
        topology_edit: compacting,
        access: { read: [], write: [topology, coordinates, properties, derived_cache] },
        may_mutate: [topology, coordinates, properties, derived_cache],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [valence, aromaticity, rings],
            preserve: [],
            invalidate: [ring_families, stereo, drawing, fingerprint],
        },
        requires_mapping: required,
        allows_noop: true,
        feature: HYDROGENS_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "strong_topology_with_coordinates",
        parity_profile: "remove_hs_parameterized_rdkit",
    }

    op with_kekulized_bonds(clear_aromatic_flags: bool) {
        method: with_kekulized_bonds,
        impl_fn: with_kekulized_bonds_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [], write: [topology, derived_cache] },
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings, valence],
            preserve: [],
            invalidate: [aromaticity, drawing, fingerprint],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: KEKULIZE_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "kekulize_clear_aromatic_flags",
    }

    op sanitized(ops: crate::SanitizeOps) {
        method: sanitized_with_ops,
        impl_fn: sanitized_impl,
        default_method: sanitized,
        default_args: [crate::SanitizeOps::ALL],
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [], write: [topology, derived_cache] },
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings, valence, aromaticity],
            preserve: [],
            invalidate: [ring_families, stereo, drawing, fingerprint],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: SANITIZE_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "sanitize_default_rdkit",
    }

    op assigned_valence() {
        method: with_assigned_valence,
        impl_fn: assigned_valence_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [derived_cache] },
        may_mutate: [derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [valence],
            preserve: [],
            invalidate: [],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: VALENCE_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_derived_cache_update",
        parity_profile: "update_property_cache_rdkit",
    }

    op assigned_rings() {
        method: with_assigned_rings,
        impl_fn: assigned_rings_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [derived_cache] },
        may_mutate: [derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings],
            preserve: [],
            invalidate: [],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: RINGS_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_derived_cache_update",
        parity_profile: "symmetrize_sssr_rdkit",
    }

    op assigned_ring_families() {
        method: with_assigned_ring_families,
        impl_fn: assigned_ring_families_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [derived_cache] },
        may_mutate: [derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [ring_families],
            preserve: [],
            invalidate: [],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: RINGS_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_derived_cache_update",
        parity_profile: "ring_families_rdkit_urf",
    }

    op assigned_aromaticity() {
        method: with_assigned_aromaticity,
        impl_fn: assigned_aromaticity_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [], write: [topology, derived_cache] },
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings, valence, aromaticity],
            preserve: [],
            invalidate: [drawing, fingerprint],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: AROMATICITY_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "set_aromaticity_default_rdkit",
    }

    op assigned_radicals() {
        method: with_assigned_radicals,
        impl_fn: assigned_radicals_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [], write: [topology, derived_cache] },
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [valence],
            preserve: [],
            invalidate: [],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: VALENCE_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "assign_radicals_rdkit",
    }

    op with_2d_coordinates() {
        method: with_2d_coordinates,
        impl_fn: with_2d_coordinates_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: COORDINATE_2D_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_generation",
        parity_profile: "compute_2d_coords_default_rdkit",
    }
}

#[mol_op_body(with_hydrogens, parts)]
fn with_hydrogens_impl(params: crate::hydrogens::AddHsParams) -> Result<OpOutcome, OperationError> {
    let mut topology = parts.begin_topology_mut()?;
    let mut coordinates = parts.begin_coordinates_mut()?;
    let mut properties = parts.begin_properties_mut()?;
    let view =
        parts.read_parts_for_blocks(topology.clone(), coordinates.clone(), properties.clone())?;
    let assignment = MoleculeReadParts::from_molecule(&view)
        .add_hs_assignment(&params)
        .map_err(|source| OperationError::AddHydrogens {
            operation: &WITH_HYDROGENS_SPEC,
            source,
        })?;

    let changed = apply_add_hs_assignment(
        parts,
        &mut topology,
        &mut coordinates,
        &mut properties,
        &assignment,
    )?;
    parts.commit_topology(topology)?;
    parts.commit_coordinates(coordinates)?;
    parts.commit_properties(properties)?;
    parts.prove_preserved(
        DerivedState::RINGS | DerivedState::RING_FAMILIES,
        PreservationProof::LeafAtomAppend,
    )?;

    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "no hydrogens were added by AddHsParameters",
        }
    })
}

fn apply_add_hs_assignment(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    coordinates: &mut CoordinateBlock,
    properties: &mut MoleculeProperties,
    assignment: &crate::hydrogens::AddHsAssignment,
) -> Result<bool, OperationError> {
    let old_atom_count = topology.atoms.len();
    let old_bond_count = topology.bonds.len();
    let added_atom_count = assignment.hydrogens_to_add.len();
    let added_bond_count = assignment.hydrogens_to_add.len();
    let changed = added_atom_count != 0
        || !assignment.atom_explicit_hydrogen_updates.is_empty()
        || !assignment.atom_pdb_residue_info_updates.is_empty()
        || !assignment.clear_isotopic_hydrogen_properties.is_empty();

    let view =
        parts.read_parts_for_blocks(topology.clone(), coordinates.clone(), properties.clone())?;
    let read = MoleculeReadParts::from_molecule(&view);
    let coords_2d_to_append = read
        .coords_2d()
        .map(|coords| add_hs_terminal_coords_2d(read, assignment, coords));
    let coords_2d_to_append = match coords_2d_to_append {
        Some(coords) => Some(coords?),
        None => None,
    };
    let conformer_coords_to_append = read
        .conformers_3d()
        .iter()
        .map(|conformer| add_hs_terminal_coords_3d(read, assignment, conformer.coords()))
        .collect::<Result<Vec<_>, _>>()?;

    for update in &assignment.atom_explicit_hydrogen_updates {
        let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs explicit hydrogen update references an out-of-range atom",
            });
        };
        atom.set_explicit_hydrogens(update.explicit_hydrogens);
    }
    for atom_id in &assignment.clear_isotopic_hydrogen_properties {
        let Some(atom) = topology.atoms.get_mut(atom_id.index()) else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs isotope property cleanup references an out-of-range atom",
            });
        };
        atom.set_tracked_isotopic_hydrogens(Vec::new());
    }
    for update in &assignment.atom_pdb_residue_info_updates {
        let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs PDB residue info update references an out-of-range atom",
            });
        };
        atom.set_pdb_residue_info(Some(update.pdb_residue_info.clone()));
    }

    for hydrogen in &assignment.hydrogens_to_add {
        if topology.atoms.get(hydrogen.heavy_atom.index()).is_none() {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs hydrogen references an out-of-range heavy atom",
            });
        }
        let atom_id = crate::AtomId::new(topology.atoms.len());
        let mut spec = crate::AtomSpec::new(crate::Element::H);
        if let Some(isotope) = hydrogen.isotope {
            spec = spec.with_isotope(isotope);
        }
        if hydrogen.is_implicit {
            spec = spec.with_implicit_hydrogen(true);
        }
        for (key, value) in &hydrogen.props {
            spec = spec.with_prop(key, value);
        }
        if let Some(info) = hydrogen.pdb_residue_info.clone() {
            spec = spec.with_pdb_residue_info(info);
        }
        topology.atoms.push(crate::Atom::from_spec(atom_id, spec));

        let bond_id = crate::BondId::new(topology.bonds.len());
        topology.bonds.push(crate::Bond::from_spec(
            bond_id,
            crate::BondSpec::new(hydrogen.heavy_atom, atom_id, crate::BondOrder::Single),
        ));
    }

    if let Some(coords) = coords_2d_to_append {
        let Some(existing) = &mut coordinates.coords_2d else {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs lost 2D coordinate block while appending hydrogen coordinates",
            });
        };
        existing.extend(coords);
    }
    if !conformer_coords_to_append.is_empty() {
        for (conformer, coords) in coordinates
            .conformers_3d
            .iter_mut()
            .zip(conformer_coords_to_append)
        {
            for coord in coords {
                conformer.push_coord(coord);
            }
        }
    }

    let mapping = TopologyMapping::with_appended(
        old_atom_count,
        old_bond_count,
        added_atom_count,
        added_bond_count,
    );
    properties.remap_topology(&mapping);

    parts.record_topology_edit(TopologyEditKind::Appending)?;
    parts.record_topology_mapping(mapping);
    parts.clear_cache(WITH_HYDROGENS_SPEC.needs_update());
    Ok(changed)
}

fn add_hs_terminal_coords_2d(
    read_parts: MoleculeReadParts<'_>,
    assignment: &crate::hydrogens::AddHsAssignment,
    coords: &[[f64; 2]],
) -> Result<Vec<[f64; 2]>, OperationError> {
    let molecule = read_parts;
    let coords_3d = coords
        .iter()
        .map(|coord| [coord[0], coord[1], 0.0])
        .collect::<Vec<_>>();
    Ok(
        add_hs_terminal_coords(molecule, assignment, &coords_3d, false)?
            .into_iter()
            .map(|coord| [coord[0], coord[1]])
            .collect(),
    )
}

fn add_hs_terminal_coords_3d(
    read_parts: MoleculeReadParts<'_>,
    assignment: &crate::hydrogens::AddHsAssignment,
    coords: &[[f64; 3]],
) -> Result<Vec<[f64; 3]>, OperationError> {
    let molecule = read_parts;
    add_hs_terminal_coords(molecule, assignment, coords, true)
}

fn add_hs_terminal_coords<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    assignment: &crate::hydrogens::AddHsAssignment,
    coords: &[[f64; 3]],
    is_3d: bool,
) -> Result<Vec<[f64; 3]>, OperationError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: void setTerminalAtomCoords(ROMol &mol, unsigned int idx, unsigned int otherIdx)
    // RDKit✔️✔️: void setTerminalAtomCoords(ROMol &mol, unsigned int idx,
    // RDKit✔️✔️:                            unsigned int otherIdx) {
    // RDKit✔️✔️:   PRECONDITION(otherIdx != idx, "degenerate atoms");
    // RDKit✔️✔️:   Atom *atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:   PRECONDITION(mol.getAtomDegree(atom) == 1, "bad atom degree");
    // RDKit✔️✔️:   const Bond *bond = mol.getBondBetweenAtoms(otherIdx, idx);
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond between atoms");
    // RDKit✔️✔️:   const Atom *otherAtom = mol.getAtomWithIdx(otherIdx);
    // RDKit✔️✔️:   double bondLength =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getRb0(1) +
    // RDKit✔️✔️:       PeriodicTable::getTable()->getRb0(otherAtom->getAtomicNum());
    // RDKit✔️✔️:   RDGeom::Point3D dirVect(0, 0, 0);
    // RDKit✔️✔️:   switch (otherAtom->getDegree()) {
    // COSMolKit computes the same sequential coordinate assignment here, then
    // writes the finished coordinates through OpParts after topology mutation.
    let old_atom_count = molecule.num_atoms();
    let mut working_coords = coords.to_vec();
    let mut adjacency = add_hs_virtual_adjacency(molecule, assignment.hydrogens_to_add.len());
    let mut appended = Vec::with_capacity(assignment.hydrogens_to_add.len());

    for (offset, hydrogen) in assignment.hydrogens_to_add.iter().enumerate() {
        let heavy_index = hydrogen.heavy_atom.index();
        let hydrogen_index = old_atom_count + offset;
        if coords.get(heavy_index).is_none() {
            return Err(OperationError::InvalidInput {
                operation: &WITH_HYDROGENS_SPEC,
                message: "AddHs coordinate append references an out-of-range heavy atom",
            });
        }
        adjacency[heavy_index].push((hydrogen_index, None));
        adjacency[hydrogen_index].push((heavy_index, None));

        let coord = if assignment.add_terminal_coordinates {
            add_hs_set_terminal_atom_coord(
                molecule,
                &adjacency,
                &working_coords,
                hydrogen_index,
                heavy_index,
                is_3d,
            )?
        } else {
            working_coords[heavy_index]
        };
        working_coords.push(coord);
        appended.push(coord);
    }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: void setTerminalAtomCoords(ROMol &mol, unsigned int idx, unsigned int otherIdx)
    Ok(appended)
}

fn add_hs_virtual_adjacency<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    hydrogens_to_add: usize,
) -> Vec<Vec<(usize, Option<crate::BondId>)>> {
    let mut adjacency = vec![Vec::new(); molecule.num_atoms() + hydrogens_to_add];
    for bond in molecule.bonds() {
        adjacency[bond.begin().index()].push((bond.end().index(), Some(bond.id())));
        adjacency[bond.end().index()].push((bond.begin().index(), Some(bond.id())));
    }
    adjacency
}

fn add_hs_set_terminal_atom_coord<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &[Vec<(usize, Option<crate::BondId>)>],
    coords: &[[f64; 3]],
    atom_index: usize,
    other_index: usize,
    is_3d: bool,
) -> Result<[f64; 3], OperationError> {
    if atom_index == other_index {
        return Err(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment received degenerate atoms",
        });
    }
    let atom_neighbors = adjacency
        .get(atom_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment references an out-of-range atom",
        })?;
    if atom_neighbors.len() != 1 {
        return Err(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment requires terminal degree one",
        });
    }
    if atom_neighbors[0].0 != other_index {
        return Err(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs terminal coordinate assignment requires a bond to the heavy atom",
        });
    }
    let other_atom = molecule
        .atoms()
        .get(other_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs coordinate append references an out-of-range heavy atom",
        })?;
    let other_pos = *coords
        .get(other_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs coordinate append references an out-of-range heavy atom",
        })?;
    let bond_length = if is_3d {
        add_hs_rdkit_rb0(1) + add_hs_rdkit_rb0(other_atom.atomic_number())
    } else {
        1.0
    };
    let neighbors = adjacency
        .get(other_index)
        .ok_or(OperationError::InvalidInput {
            operation: &WITH_HYDROGENS_SPEC,
            message: "AddHs coordinate append references an out-of-range heavy atom",
        })?
        .iter()
        .filter_map(|(neighbor, _)| (*neighbor != atom_index).then_some(*neighbor))
        .collect::<Vec<_>>();

    let mut dir: [f64; 3];
    match neighbors.len() + 1 {
        1 => {
            // RDKit✔️✔️:     case 1:
            // RDKit✔️✔️:       // No other atoms present.
            // RDKit✔️✔️: if ((*cfi)->is3D()) { dirVect.z = 1; } else { dirVect.x = 1; }
            dir = if is_3d {
                [0.0, 0.0, 1.0]
            } else {
                [1.0, 0.0, 0.0]
            };
        }
        2 => {
            // RDKit✔️✔️:     case 2:
            // RDKit✔️✔️:       // One other neighbor.
            // RDKit✔️✔️: nbr1Vect = nbr1Pos - otherPos; ... nbr1Vect *= -1;
            let nbr1 = neighbors[0];
            let Some(mut nbr1_vec) = add_hs_direction_away(coords, other_pos, nbr1) else {
                return Ok(other_pos);
            };
            match other_atom.hybridization() {
                crate::Hybridization::Sp3 => {
                    // RDKit✔️✔️: tform.SetRotation((180 - 109.471) * M_PI / 180., perpVect);
                    let perp = if is_3d {
                        add_hs_perpendicular(nbr1_vec)
                    } else {
                        [0.0, 0.0, 1.0]
                    };
                    dir = add_hs_rotate(nbr1_vec, perp, (180.0 - 109.471_f64).to_radians());
                }
                crate::Hybridization::Sp2 => {
                    // RDKit✔️✔️: default 3D position is an arbitrary perpendicular; for 2D z-normal.
                    let mut perp = if is_3d {
                        add_hs_perpendicular(nbr1_vec)
                    } else {
                        [0.0, 0.0, 1.0]
                    };
                    if adjacency[nbr1].len() > 1
                        && add_hs_bond_is_pi_like(molecule, adjacency, other_index, nbr1)
                        && let Some(nbr2) = adjacency[nbr1].iter().find_map(|(neighbor, _)| {
                            (*neighbor != other_index).then_some(*neighbor)
                        })
                        && let (Some(nbr1_pos), Some(nbr2_pos)) =
                            (coords.get(nbr1).copied(), coords.get(nbr2).copied())
                        && let Some(nbr2_vec) =
                            add_hs_normalized(add_hs_vec_sub(nbr2_pos, nbr1_pos))
                    {
                        let cross = add_hs_cross(nbr2_vec, nbr1_vec);
                        if add_hs_len_sq(cross) >= ADD_HS_SQ_DIST_ZERO_TOL {
                            perp = cross;
                        }
                    }
                    perp = add_hs_normalized(perp).unwrap_or(perp);
                    dir = add_hs_rotate(nbr1_vec, perp, 60.0_f64.to_radians());
                }
                crate::Hybridization::Sp => {
                    // RDKit✔️✔️: just lay the H along the vector.
                    dir = nbr1_vec;
                }
                _ => {
                    // RDKit✔️✔️: default branch lays the H along the vector.
                    dir = nbr1_vec;
                }
            }
            nbr1_vec = dir;
            dir = nbr1_vec;
        }
        3 => {
            // RDKit✔️✔️:     case 3:
            // RDKit✔️✔️:       // Two other neighbors.
            // RDKit✔️✔️: start along the average of the two vectors.
            let Some(nbr1_vec) = add_hs_direction_away(coords, other_pos, neighbors[0]) else {
                return Ok(other_pos);
            };
            let Some(nbr2_vec) = add_hs_direction_away(coords, other_pos, neighbors[1]) else {
                return Ok(other_pos);
            };
            let Some(mut avg_dir) = add_hs_normalized(add_hs_vec_add(nbr1_vec, nbr2_vec)) else {
                return Ok(other_pos);
            };
            if is_3d {
                match other_atom.hybridization() {
                    crate::Hybridization::Sp3 => {
                        // RDKit✔️✔️: rotate about the perpendicular-to-neighbors axis.
                        let nbr_perp = add_hs_cross(nbr1_vec, nbr2_vec);
                        let rotn_axis = add_hs_cross(nbr_perp, avg_dir);
                        if let Some(axis) = add_hs_normalized(rotn_axis) {
                            avg_dir =
                                add_hs_rotate(avg_dir, axis, (109.471_f64 / 2.0).to_radians());
                        }
                    }
                    crate::Hybridization::Sp2 => {}
                    _ => {}
                }
            }
            dir = avg_dir;
        }
        4 => {
            // RDKit✔️✔️:     case 4:
            // RDKit✔️✔️:       // Three other neighbors.
            // RDKit✔️✔️: three other neighbors.
            let Some(nbr1_vec) = add_hs_direction_away(coords, other_pos, neighbors[0]) else {
                return Ok(other_pos);
            };
            let Some(nbr2_vec) = add_hs_direction_away(coords, other_pos, neighbors[1]) else {
                return Ok(other_pos);
            };
            let Some(nbr3_vec) = add_hs_direction_away(coords, other_pos, neighbors[2]) else {
                return Ok(other_pos);
            };
            if is_3d {
                if add_hs_dot(nbr3_vec, add_hs_cross(nbr1_vec, nbr2_vec)).abs() < 0.1 {
                    dir = add_hs_cross(nbr1_vec, nbr2_vec);
                    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
                        dir = add_hs_vec_scale(add_hs_cross(nbr1_vec, nbr3_vec), -1.0);
                    }
                    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
                        dir = add_hs_cross(nbr2_vec, nbr3_vec);
                    }
                    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
                        return Ok(other_pos);
                    }
                    if other_atom.prop("_CIPCode").is_some() {
                        let v1 = add_hs_vec_sub(dir, nbr3_vec);
                        let v2 = add_hs_vec_sub(nbr1_vec, nbr3_vec);
                        let v3 = add_hs_vec_sub(nbr2_vec, nbr3_vec);
                        let vol = add_hs_dot(v1, add_hs_cross(v2, v3));
                        if (other_atom.chiral_tag() == crate::ChiralTag::TetrahedralCcw
                            && vol < 0.0)
                            || (other_atom.chiral_tag() == crate::ChiralTag::TetrahedralCw
                                && vol > 0.0)
                        {
                            dir = add_hs_vec_scale(dir, -1.0);
                        }
                    }
                } else {
                    dir = add_hs_vec_add(add_hs_vec_add(nbr1_vec, nbr2_vec), nbr3_vec);
                }
            } else {
                // RDKit✔️✔️: pick the bisector opposite the smallest accumulated outer angle.
                let angle12 = add_hs_angle(nbr1_vec, nbr2_vec);
                let angle13 = add_hs_angle(nbr1_vec, nbr3_vec);
                let angle23 = add_hs_angle(nbr2_vec, nbr3_vec);
                let accum1 = angle12 + angle13;
                let accum2 = angle12 + angle23;
                let accum3 = angle13 + angle23;
                dir = if accum1 <= accum2 && accum1 <= accum3 {
                    add_hs_pick_bisector(nbr1_vec, nbr2_vec, nbr3_vec)
                } else if accum2 <= accum1 && accum2 <= accum3 {
                    add_hs_pick_bisector(nbr2_vec, nbr1_vec, nbr3_vec)
                } else {
                    add_hs_pick_bisector(nbr3_vec, nbr1_vec, nbr2_vec)
                };
            }
        }
        _ => {
            // RDKit✔️✔️:     default:
            // RDKit✔️✔️:       atomPos = otherPos + dirVect * bondLength;
            // RDKit✔️✔️:       for (...) { (*cfi)->setAtomPos(idx, atomPos); }
            return Ok([0.0, 0.0, 0.0]);
        }
    }

    let dir = add_hs_normalized(dir).unwrap_or(dir);
    Ok(add_hs_vec_add(
        other_pos,
        add_hs_vec_scale(dir, bond_length),
    ))
}

const ADD_HS_SQ_DIST_ZERO_TOL: f64 = 1.0e-8;

fn add_hs_rdkit_rb0(atomic_number: u8) -> f64 {
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
    // RDKit✔️✔️: double Rb0() const { return rB0; }
    // RDKit✔️✔️: double getRb0(UINT atomicNumber) const { ... return byanum[atomicNumber].Rb0(); }
    const RB0: [f64; 113] = [
        0.0, 0.33, 0.7, 1.23, 0.9, 0.82, 0.77, 0.7, 0.66, 0.611, 0.7, 1.54, 1.36, 1.18, 0.937,
        0.89, 1.04, 0.997, 1.74, 2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17,
        1.25, 1.26, 1.188, 1.2, 1.17, 1.167, 1.91, 2.16, 1.91, 1.62, 1.45, 1.34, 1.3, 1.27, 1.25,
        1.25, 1.28, 1.34, 1.48, 1.44, 1.385, 1.4, 1.378, 1.387, 1.98, 2.35, 1.98, 1.69, 1.83, 1.82,
        1.81, 1.8, 1.8, 1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72, 1.44, 1.34, 1.3,
        1.28, 1.26, 1.27, 1.3, 1.34, 1.49, 1.48, 1.48, 1.45, 1.46, 1.45, 2.4, 2.0, 1.9, 1.88, 1.79,
        1.61, 1.58, 1.55, 1.53, 1.07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    RB0.get(usize::from(atomic_number)).copied().unwrap_or(0.0)
    // END RDKIT CPP FUNCTION PeriodicTable::getRb0 / atomicData::Rb0
}

fn add_hs_bond_is_pi_like<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &[Vec<(usize, Option<crate::BondId>)>],
    begin: usize,
    end: usize,
) -> bool {
    adjacency[begin]
        .iter()
        .find(|(neighbor, _)| *neighbor == end)
        .and_then(|(_, bond)| *bond)
        .and_then(|bond| molecule.bonds().get(bond.index()))
        .is_some_and(|bond| {
            bond.is_aromatic() || bond.order() == crate::BondOrder::Double || bond.is_conjugated()
        })
}

fn add_hs_direction_away(
    coords: &[[f64; 3]],
    other_pos: [f64; 3],
    neighbor: usize,
) -> Option<[f64; 3]> {
    let neighbor_pos = *coords.get(neighbor)?;
    add_hs_normalized(add_hs_vec_scale(
        add_hs_vec_sub(neighbor_pos, other_pos),
        -1.0,
    ))
}

fn add_hs_pick_bisector(nbr1_vec: [f64; 3], nbr2_vec: [f64; 3], nbr3_vec: [f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: RDGeom::Point3D pickBisector(const RDGeom::Point3D&, const RDGeom::Point3D&, const RDGeom::Point3D&)
    // RDKit✔️✔️: auto dirVect = nbr2Vect + nbr3Vect;
    let mut dir = add_hs_vec_add(nbr2_vec, nbr3_vec);
    if add_hs_len_sq(dir) < ADD_HS_SQ_DIST_ZERO_TOL {
        // RDKit✔️✔️: std::swap(dirVect.x, dirVect.y); dirVect.x *= -1;
        dir = [-nbr2_vec[1], nbr2_vec[0], nbr2_vec[2]];
    }
    if add_hs_dot(dir, nbr1_vec) < 0.0 {
        dir = add_hs_vec_scale(dir, -1.0);
    }
    dir
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/AddHs.cpp :: RDGeom::Point3D pickBisector(const RDGeom::Point3D&, const RDGeom::Point3D&, const RDGeom::Point3D&)
}

fn add_hs_vec_add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn add_hs_vec_sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn add_hs_vec_scale(vec: [f64; 3], scale: f64) -> [f64; 3] {
    [vec[0] * scale, vec[1] * scale, vec[2] * scale]
}

fn add_hs_dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn add_hs_cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn add_hs_len_sq(vec: [f64; 3]) -> f64 {
    add_hs_dot(vec, vec)
}

fn add_hs_normalized(vec: [f64; 3]) -> Option<[f64; 3]> {
    let len_sq = add_hs_len_sq(vec);
    if len_sq < ADD_HS_SQ_DIST_ZERO_TOL {
        return None;
    }
    Some(add_hs_vec_scale(vec, 1.0 / len_sq.sqrt()))
}

fn add_hs_perpendicular(vec: [f64; 3]) -> [f64; 3] {
    let axis = if vec[0].abs() < vec[1].abs() {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    add_hs_normalized(add_hs_cross(vec, axis)).unwrap_or([0.0, 0.0, 1.0])
}

fn add_hs_rotate(vec: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let axis = add_hs_normalized(axis).unwrap_or([0.0, 0.0, 1.0]);
    let cos = angle.cos();
    let sin = angle.sin();
    let term1 = add_hs_vec_scale(vec, cos);
    let term2 = add_hs_vec_scale(add_hs_cross(axis, vec), sin);
    let term3 = add_hs_vec_scale(axis, add_hs_dot(axis, vec) * (1.0 - cos));
    add_hs_vec_add(add_hs_vec_add(term1, term2), term3)
}

fn add_hs_angle(left: [f64; 3], right: [f64; 3]) -> f64 {
    let denom = (add_hs_len_sq(left) * add_hs_len_sq(right)).sqrt();
    if denom < ADD_HS_SQ_DIST_ZERO_TOL {
        return 0.0;
    }
    (add_hs_dot(left, right) / denom).clamp(-1.0, 1.0).acos()
}

#[mol_op_body(without_hydrogens, parts)]
fn without_hydrogens_impl(sanitize: bool) -> Result<OpOutcome, OperationError> {
    let params = crate::hydrogens::RemoveHsParams::default();
    without_hydrogens_apply(parts, &params, sanitize)
}

#[mol_op_body(without_hydrogens_with_params, parts)]
fn without_hydrogens_with_params_impl(
    params: crate::hydrogens::RemoveHsParams,
    sanitize: bool,
) -> Result<OpOutcome, OperationError> {
    without_hydrogens_apply(parts, &params, sanitize)
}

fn without_hydrogens_apply(
    parts: &mut OpParts<'_>,
    params: &crate::hydrogens::RemoveHsParams,
    sanitize: bool,
) -> Result<OpOutcome, OperationError> {
    let operation = parts.spec;
    let mut topology = parts.begin_topology_mut()?;
    let mut coordinates = parts.begin_coordinates_mut()?;
    let mut properties = parts.begin_properties_mut()?;
    let view =
        parts.read_parts_for_blocks(topology.clone(), coordinates.clone(), properties.clone())?;
    let assignment = MoleculeReadParts::from_molecule(&view)
        .remove_hs_assignment(params, sanitize)
        .map_err(|source| OperationError::RemoveHydrogens { operation, source })?;

    let mut changed = apply_remove_hs_assignment(parts, &mut topology, &assignment)?;
    let atoms_to_remove = assignment.atoms_to_remove.clone();
    if atoms_to_remove.is_empty() {
        let atom_count = topology.atoms.len();
        let bond_count = topology.bonds.len();
        parts.record_topology_edit(TopologyEditKind::Compacting)?;
        parts.record_topology_mapping(TopologyMapping::identity(atom_count, bond_count));
        parts.clear_cache(WITHOUT_HYDROGENS_SPEC.needs_update());
    } else {
        let mapping = topology.remove_atoms_with_mapping(&atoms_to_remove);
        coordinates.remap_topology(&mapping);
        properties.remap_topology(&mapping);
        parts.record_topology_mapping(mapping);
        parts.record_topology_edit(TopologyEditKind::Compacting)?;
        parts.clear_cache(WITHOUT_HYDROGENS_SPEC.needs_update());
        changed = true;
    }
    if assignment.sanitize_after_removal {
        changed |=
            sanitize_after_remove_hs_removal(parts, &mut topology, &coordinates, &properties)?;
    }
    parts.commit_topology(topology)?;
    parts.commit_coordinates(coordinates)?;
    parts.commit_properties(properties)?;

    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "no removable explicit hydrogens matched RemoveHsParameters",
        }
    })
}

fn apply_remove_hs_assignment(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    assignment: &crate::hydrogens::RemoveHsAssignment,
) -> Result<bool, OperationError> {
    let operation = parts.spec;
    let has_local_updates = !assignment.atom_explicit_hydrogen_updates.is_empty()
        || !assignment.atom_chirality_inversions.is_empty()
        || !assignment.atom_property_updates.is_empty()
        || !assignment.bond_direction_updates.is_empty()
        || !assignment.bond_stereo_updates.is_empty()
        || !assignment.bond_stereo_atom_replacements.is_empty()
        || !assignment.sgroup_updates.is_empty();

    if has_local_updates {
        for update in &assignment.atom_explicit_hydrogen_updates {
            let Some(atom) = topology.atoms.get_mut(update.atom.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H explicit hydrogen update references an out-of-range atom",
                });
            };
            atom.set_explicit_hydrogens(update.explicit_hydrogens);
        }
        for atom_id in &assignment.atom_chirality_inversions {
            let Some(atom) = topology.atoms.get_mut(atom_id.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H chirality update references an out-of-range atom",
                });
            };
            atom.set_chiral_tag(match atom.chiral_tag() {
                crate::ChiralTag::TetrahedralCw => crate::ChiralTag::TetrahedralCcw,
                crate::ChiralTag::TetrahedralCcw => crate::ChiralTag::TetrahedralCw,
                other => other,
            });
        }
        for update in &assignment.atom_property_updates {
            match update {
                crate::hydrogens::AtomPropertyUpdate::SetUnknownStereo { atom } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H atom property update references an out-of-range atom",
                        });
                    };
                    atom.set_unknown_stereo(true);
                }
                crate::hydrogens::AtomPropertyUpdate::SetIsotopicHydrogens { atom, isotopes } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H isotope property update references an out-of-range atom",
                        });
                    };
                    atom.set_tracked_isotopic_hydrogens(isotopes.clone());
                }
                crate::hydrogens::AtomPropertyUpdate::ClearIsotopicHydrogens { atom } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H isotope property cleanup references an out-of-range atom",
                        });
                    };
                    atom.set_tracked_isotopic_hydrogens(Vec::new());
                }
                crate::hydrogens::AtomPropertyUpdate::ClearExcessChiralExplicitHydrogens {
                    atom,
                } => {
                    let Some(atom) = topology.atoms.get_mut(atom.index()) else {
                        return Err(OperationError::InvalidInput {
                            operation,
                            message: "remove-H chiral explicit-H cleanup references an out-of-range atom",
                        });
                    };
                    atom.set_explicit_hydrogens(0);
                }
            }
        }
        for update in &assignment.bond_direction_updates {
            let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H bond direction update references an out-of-range bond",
                });
            };
            bond.set_direction(update.direction);
        }
        for update in &assignment.bond_stereo_updates {
            let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H bond stereo update references an out-of-range bond",
                });
            };
            bond.set_stereo_atoms(update.stereo_atoms);
            bond.set_stereo(update.stereo);
        }
        for update in &assignment.bond_stereo_atom_replacements {
            let Some(bond) = topology.bonds.get_mut(update.bond.index()) else {
                return Err(OperationError::InvalidInput {
                    operation,
                    message: "remove-H stereo atom replacement references an out-of-range bond",
                });
            };
            let Some([left, right]) = bond.stereo_atoms() else {
                continue;
            };
            let stereo_atoms = [
                if left == update.old_atom {
                    update.new_atom
                } else {
                    left
                },
                if right == update.old_atom {
                    update.new_atom
                } else {
                    right
                },
            ];
            bond.set_stereo_atoms(Some(stereo_atoms));
        }
        for update in &assignment.sgroup_updates {
            for substance_group in &mut topology.substance_groups {
                match update {
                    crate::hydrogens::SGroupRemoveHsUpdate::RemoveBond { bond } => {
                        substance_group.remove_bond(*bond);
                    }
                    crate::hydrogens::SGroupRemoveHsUpdate::RemoveAtom { atom } => {
                        substance_group.remove_atom(*atom);
                    }
                    crate::hydrogens::SGroupRemoveHsUpdate::RemoveParentAtom { atom } => {
                        substance_group.remove_parent_atom(*atom);
                    }
                    crate::hydrogens::SGroupRemoveHsUpdate::ClearAttachPointLeavingAtom {
                        atom,
                    } => {
                        substance_group.clear_attach_point_leaving_atom(*atom);
                    }
                }
            }
        }
    }
    Ok(has_local_updates)
}

fn sanitize_after_remove_hs_removal(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<bool, OperationError> {
    // BEGIN RDKIT CPP BODY MolOps::removeHs sanitize-after-removal
    // RDKit✔️✔️:   if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    // RDKit✔️✔️:     sanitizeMol(mol);
    // RDKit✔️✔️:   }
    let changed = run_sanitize_pipeline_on_topology(
        parts,
        topology,
        Some(coordinates),
        Some(properties),
        crate::SanitizeOps::ALL,
        &WITHOUT_HYDROGENS_SPEC,
    )?;
    // END RDKIT CPP BODY MolOps::removeHs sanitize-after-removal
    Ok(changed)
}

#[mol_op_body(with_kekulized_bonds, parts)]
fn with_kekulized_bonds_impl(clear_aromatic_flags: bool) -> Result<OpOutcome, OperationError> {
    // RDKit✔️✔️: void kekulizeMol(ROMol &mol, bool clearAromaticFlags = false,
    // RDKit✔️✔️:                  bool canonical = true) {
    // RDKit✔️✔️:   auto &wmol = static_cast<RWMol &>(mol);
    // RDKit✔️✔️:   MolOps::Kekulize(wmol, clearAromaticFlags, canonical);
    // RDKit✔️✔️: }
    let mut topology = parts.begin_topology_mut()?;
    let view = parts.read_parts_for_topology(topology.clone())?;
    let rings = MoleculeReadParts::from_molecule(&view)
        .symmetrize_sssr()
        .map_err(|source| OperationError::RingFinding {
            operation: &WITH_KEKULIZED_BONDS_SPEC,
            source,
        })?;
    parts.set_rings_cache(rings);
    let view = parts.read_parts_for_topology(topology.clone())?;
    let ring_info = MoleculeReadParts::from_molecule(&view)
        .derived_cache()
        .rings
        .as_ref()
        .expect("rings were recomputed immediately above")
        .clone();
    let assignment = MoleculeReadParts::from_molecule(&view)
        .kekulize_assignment(Some(&ring_info), clear_aromatic_flags, true, 100)
        .map_err(|source| OperationError::Kekulize {
            operation: &WITH_KEKULIZED_BONDS_SPEC,
            source,
        })?;

    let changed = crate::kekulize::apply_kekulize_assignment(&mut topology, &assignment);
    let view = parts.read_parts_for_topology(topology.clone())?;
    let valence = MoleculeReadParts::from_molecule(&view)
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
        .map_err(|source| OperationError::Valence {
            operation: &WITH_KEKULIZED_BONDS_SPEC,
            source,
        })?;
    parts.commit_topology(topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(DerivedState::AROMATICITY);
    parts.set_valence_cache(valence);
    parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "kekulization assignment produced no effective topology-state change",
        }
    })
}

#[mol_op_body(sanitized, parts)]
fn sanitized_impl(ops: crate::SanitizeOps) -> Result<OpOutcome, OperationError> {
    let mut parts = parts;
    let changed = run_sanitize_pipeline(&mut parts, ops, &SANITIZED_SPEC)?;
    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "requested sanitize steps only updated derived state",
        }
    })
}

fn sanitize_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    message: String,
    unsupported: Option<crate::UnsupportedFeatureError>,
) -> OperationError {
    OperationError::Sanitize {
        operation,
        source: crate::SanitizeError {
            step,
            message,
            unsupported,
        },
    }
}

fn sanitize_invalid_input_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    message: &'static str,
) -> OperationError {
    sanitize_error(operation, step, message.to_string(), None)
}

fn sanitize_valence_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::ValenceError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_ring_finding_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::RingFindingError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_kekulize_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::KekulizeError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_aromaticity_error(
    operation: &'static MoleculeOpSpec,
    step: crate::SanitizeStep,
    source: crate::AromaticityError,
) -> OperationError {
    sanitize_error(operation, step, source.to_string(), None)
}

fn sanitize_stage<T, E>(
    step: crate::SanitizeStep,
    stage: impl FnOnce() -> Result<T, E>,
    map_error: impl FnOnce(crate::SanitizeStep, E) -> OperationError,
) -> Result<T, OperationError> {
    stage().map_err(|source| map_error(step, source))
}

fn sanitize_properties_stage(
    parts: &mut OpParts<'_>,
    topology: &TopologyBlock,
    coordinates: Option<&CoordinateBlock>,
    properties: Option<&MoleculeProperties>,
    strict: bool,
    operation: &'static MoleculeOpSpec,
) -> Result<crate::ValenceAssignment, OperationError> {
    sanitize_stage(
        crate::SanitizeStep::Properties,
        || {
            sanitize_recompute_property_cache(
                parts,
                topology,
                coordinates,
                properties,
                strict,
                crate::SanitizeStep::Properties,
                operation,
            )
        },
        |_step, error| error,
    )
}

fn sanitize_recompute_property_cache(
    parts: &mut OpParts<'_>,
    topology: &TopologyBlock,
    coordinates: Option<&CoordinateBlock>,
    properties: Option<&MoleculeProperties>,
    strict: bool,
    step: crate::SanitizeStep,
    operation: &'static MoleculeOpSpec,
) -> Result<crate::ValenceAssignment, OperationError> {
    let view = parts.read_parts_for_optional_blocks(topology.clone(), coordinates, properties)?;
    let valence = MoleculeReadParts::from_molecule(&view)
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, strict)
        .map_err(|source| sanitize_valence_error(operation, step, source))?;
    parts.set_valence_cache(valence.clone());
    Ok(valence)
}

fn run_sanitize_pipeline(
    parts: &mut OpParts<'_>,
    ops: crate::SanitizeOps,
    operation: &'static MoleculeOpSpec,
) -> Result<bool, OperationError> {
    let mut topology = parts.begin_topology_mut()?;
    let changed =
        run_sanitize_pipeline_on_topology(parts, &mut topology, None, None, ops, operation)?;
    parts.commit_topology(topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    Ok(changed)
}

fn run_sanitize_pipeline_on_topology(
    parts: &mut OpParts<'_>,
    topology: &mut TopologyBlock,
    coordinates: Option<&CoordinateBlock>,
    properties: Option<&MoleculeProperties>,
    ops: crate::SanitizeOps,
    operation: &'static MoleculeOpSpec,
) -> Result<bool, OperationError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::sanitizeMol
    // RDKit✔️✔️: void sanitizeMol(RWMol &mol) {
    // RDKit✔️✔️:   unsigned int failedOp = 0;
    // RDKit✔️✔️:   sanitizeMol(mol, failedOp, SANITIZE_ALL);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
    // RDKit✔️✔️:                  unsigned int sanitizeOps) {
    // RDKit✔️✔️:   // clear out any cached properties
    // RDKit✔️✔️:   mol.clearComputedProps();
    parts.clear_cache(SANITIZED_SPEC.needs_update());
    macro_rules! sanitize_read {
        ($read:ident => $body:expr) => {{
            let view = parts
                .read_parts_for_optional_blocks((*topology).clone(), coordinates, properties)
                .expect("sanitize working topology should satisfy molecule invariants");
            let $read = MoleculeReadParts::from_molecule(&view);
            $body
        }};
    }

    let mut changed = false;
    let property_cache_strict = ops.contains(crate::SanitizeOps::PROPERTIES);

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUP;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     // clean up things like nitro groups
    // RDKit✔️✔️:     cleanUp(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP) {
        let cleanup = sanitize_stage(
            crate::SanitizeStep::Cleanup,
            || sanitize_read!(read => sanitize_cleanup_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let cleanup_changed = sanitize_read!(read => cleanup_changes_molecule(read, &cleanup));
        if cleanup_changed {
            apply_sanitize_cleanup_assignment(topology, cleanup);
        }
        changed |= cleanup_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     cleanUpOrganometallics(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS) {
        let organometallics = sanitize_stage(
            crate::SanitizeStep::CleanupOrganometallics,
            || sanitize_read!(read => sanitize_organometallic_cleanup_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let organometallics_changed =
            sanitize_read!(read => organometallic_cleanup_changes_molecule(read, &organometallics));
        if organometallics_changed {
            apply_sanitize_organometallic_cleanup_assignment(topology, organometallics);
        }
        changed |= organometallics_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     mol.updatePropertyCache(true);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    let valence = sanitize_properties_stage(
        parts,
        topology,
        coordinates,
        properties,
        property_cache_strict,
        operation,
    )?;
    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SYMMRINGS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     VECT_INT_VECT arings;
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol, arings);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SYMMRINGS) {
        let rings = sanitize_stage(
            crate::SanitizeStep::SymmRings,
            || sanitize_read!(read => read.symmetrize_sssr()),
            |step, source| sanitize_ring_finding_error(operation, step, source),
        )?;
        parts.set_rings_cache(rings);
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_KEKULIZE;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     Kekulize(mol, true, false);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::KEKULIZE) {
        if sanitize_read!(read => read.derived_cache().rings.is_none()) {
            let rings = sanitize_stage(
                crate::SanitizeStep::Kekulize,
                || sanitize_read!(read => read.symmetrize_sssr()),
                |step, source| sanitize_ring_finding_error(operation, step, source),
            )?;
            parts.set_rings_cache(rings);
        }
        let ring_info = sanitize_read!(read => {
            read.derived_cache()
                .rings
                .as_ref()
                .expect("rings were recomputed immediately above")
                .clone()
        });
        let assignment = sanitize_stage(
            crate::SanitizeStep::Kekulize,
            || sanitize_read!(read => read.kekulize_assignment(Some(&ring_info), true, false, 100)),
            |step, source| sanitize_kekulize_error(operation, step, source),
        )?;
        let kekulize_changed = crate::kekulize::apply_kekulize_assignment(topology, &assignment);
        if kekulize_changed {
            sanitize_stage(
                crate::SanitizeStep::Kekulize,
                || {
                    sanitize_recompute_property_cache(
                        parts,
                        topology,
                        coordinates,
                        properties,
                        property_cache_strict,
                        crate::SanitizeStep::Kekulize,
                        operation,
                    )
                },
                |_step, error| error,
            )?;
            parts.clear_cache(
                DerivedState::AROMATICITY | DerivedState::DRAWING | DerivedState::FINGERPRINT,
            );
        }
        changed |= kekulize_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_FINDRADICALS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     assignRadicals(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::FIND_RADICALS) {
        let radicals = sanitize_stage(
            crate::SanitizeStep::FindRadicals,
            || sanitize_read!(read => read.assign_radicals()),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let radicals_changed = sanitize_read!(read => read
            .atoms()
            .iter()
            .zip(radicals.iter().copied())
            .any(|(atom, radical)| atom.radical_electrons() != radical));
        if radicals_changed {
            for (atom, radical) in topology.atoms.iter_mut().zip(radicals) {
                atom.set_radical_electrons(radical);
            }
        }
        changed |= radicals_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SETAROMATICITY;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     setAromaticity(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SET_AROMATICITY) {
        if sanitize_read!(read => read.derived_cache().rings.is_none()) {
            let rings = sanitize_stage(
                crate::SanitizeStep::SetAromaticity,
                || sanitize_read!(read => read.symmetrize_sssr()),
                |step, source| sanitize_ring_finding_error(operation, step, source),
            )?;
            parts.set_rings_cache(rings);
        }
        let assignment = sanitize_stage(
            crate::SanitizeStep::SetAromaticity,
            || sanitize_read!(read => read.set_aromaticity(crate::AromaticityModel::Default)),
            |step, source| sanitize_aromaticity_error(operation, step, source),
        )?;
        for (atom, is_aromatic) in topology
            .atoms
            .iter_mut()
            .zip(assignment.atom_aromatic.iter().copied())
        {
            atom.set_aromatic(is_aromatic);
        }
        for (bond, is_aromatic) in topology
            .bonds
            .iter_mut()
            .zip(assignment.bond_aromatic.iter().copied())
        {
            bond.set_aromatic(is_aromatic);
            if is_aromatic
                && matches!(
                    bond.order(),
                    crate::BondOrder::Single | crate::BondOrder::Double
                )
            {
                bond.set_order(crate::BondOrder::Aromatic);
            }
        }
        parts.mark_aromaticity_valid();
        parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        changed = true;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SETCONJUGATION;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     setConjugation(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SET_CONJUGATION) {
        let conjugation = sanitize_stage(
            crate::SanitizeStep::SetConjugation,
            || sanitize_read!(read => sanitize_conjugation_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let conjugation_changed = sanitize_read!(read => read
            .bonds()
            .iter()
            .zip(conjugation.iter().copied())
            .any(|(bond, is_conjugated)| bond.is_conjugated() != is_conjugated));
        if conjugation_changed {
            for (bond, is_conjugated) in topology.bonds.iter_mut().zip(conjugation) {
                bond.set_conjugated(is_conjugated);
            }
            parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        }
        changed |= conjugation_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_SETHYBRIDIZATION;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     setHybridization(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::SET_HYBRIDIZATION) {
        let hybridization = sanitize_stage(
            crate::SanitizeStep::SetHybridization,
            || sanitize_read!(read => sanitize_hybridization_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let hybridization_changed = sanitize_read!(read => read
            .atoms()
            .iter()
            .zip(hybridization.iter().copied())
            .any(|(atom, hybridization)| atom.hybridization() != hybridization));
        if hybridization_changed {
            for (atom, hybridization) in topology.atoms.iter_mut().zip(hybridization) {
                atom.set_hybridization(hybridization);
            }
            parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        }
        changed |= hybridization_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUPATROPISOMERS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     cleanupAtropisomers(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_ATROPISOMERS) {
        if sanitize_read!(read => read.derived_cache().rings.is_none()) {
            let rings = sanitize_stage(
                crate::SanitizeStep::CleanupAtropisomers,
                || sanitize_read!(read => read.symmetrize_sssr()),
                |step, source| sanitize_ring_finding_error(operation, step, source),
            )?;
            parts.set_rings_cache(rings);
        }
        let ring_info = sanitize_read!(read => {
            read.derived_cache()
                .rings
                .as_ref()
                .expect("rings were recomputed immediately above")
                .clone()
        });
        let atropisomer_cleanup =
            sanitize_read!(read => sanitize_cleanup_atropisomers_assignment(read, &ring_info));
        let atropisomer_cleanup_changed = sanitize_read!(read => {
            atropisomer_cleanup
                .iter()
                .any(|(bond, _)| read.bonds()[bond.index()].stereo() != crate::BondStereo::None)
        });
        if atropisomer_cleanup_changed {
            apply_sanitize_cleanup_atropisomers_assignment(topology, &atropisomer_cleanup);
            parts.clear_cache(
                DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT,
            );
        }
        changed |= atropisomer_cleanup_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     cleanupChirality(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_CHIRALITY) {
        let chirality_cleanup = sanitize_stage(
            crate::SanitizeStep::CleanupChirality,
            || sanitize_read!(read => sanitize_cleanup_chirality_assignment(read)),
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let chirality_cleanup_changed =
            sanitize_read!(read => cleanup_chirality_changes_molecule(read, &chirality_cleanup));
        if chirality_cleanup_changed {
            apply_sanitize_cleanup_chirality_assignment(topology, &chirality_cleanup);
            parts.clear_cache(
                DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT,
            );
        }
        changed |= chirality_cleanup_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_ADJUSTHS;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     adjustHs(mol);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::ADJUST_HYDROGENS) {
        let hydrogen_adjustment = sanitize_stage(
            crate::SanitizeStep::AdjustHydrogens,
            || {
                sanitize_read!(read => {
                    sanitize_adjust_hydrogens_assignment(read)
                })
            },
            |step, source| sanitize_valence_error(operation, step, source),
        )?;
        let hydrogen_adjustment_changed =
            sanitize_read!(read => adjust_hydrogens_changes_molecule(read, &hydrogen_adjustment));
        if hydrogen_adjustment_changed {
            apply_sanitize_adjust_hydrogens_assignment(topology, &hydrogen_adjustment);
            parts.clear_cache(DerivedState::VALENCE);
            parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        }
        changed |= hydrogen_adjustment_changed;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     mol.updatePropertyCache(true);
    // RDKit✔️✔️:   }
    if ops.contains(crate::SanitizeOps::PROPERTIES) {
        sanitize_properties_stage(parts, topology, coordinates, properties, true, operation)?;
    }
    // RDKit✔️✔️:   operationThatFailed = 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::sanitizeMol
    Ok(changed)
}

fn sanitize_conjugation_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<bool>, crate::ValenceError> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION MolOps::setConjugation
    // RDKit✔️✔️: void setConjugation(ROMol &mol) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     bond->setIsConjugated(bond->getIsAromatic());
    // RDKit✔️✔️:   }
    let adjacency = sanitize_adjacency(molecule)?;
    let valence = sanitize_valence_facts(molecule)?;
    let mut conjugated = molecule
        .bonds()
        .iter()
        .map(crate::Bond::is_aromatic)
        .collect::<Vec<_>>();

    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     markConjAtomBonds(atom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for atom in molecule.atoms() {
        sanitize_mark_conjugated_atom_bonds(
            molecule,
            &adjacency,
            &valence,
            atom.id(),
            &mut conjugated,
        )?;
    }
    Ok(conjugated)
    // END RDKIT CPP FUNCTION MolOps::setConjugation
}

fn sanitize_mark_conjugated_atom_bonds<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
    conjugated: &mut [bool],
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION markConjAtomBonds
    // RDKit✔️✔️: void markConjAtomBonds(Atom *at) {
    // RDKit✔️✔️:   if (!isAtomConjugCand(at)) { return; }
    if !sanitize_is_atom_conjugation_candidate(molecule, adjacency, valence, atom)? {
        return Ok(());
    }

    // RDKit✔️✔️:   int sbo = at->getDegree() + at->getTotalNumHs();
    // RDKit✔️✔️:   if ((sbo < 2) || (sbo > 3)) { return; }
    let sbo = sanitize_total_degree(molecule, adjacency, valence, atom)?;
    if !(2..=3).contains(&sbo) {
        return Ok(());
    }

    // RDKit✔️✔️:   for (const auto bnd1 : mol.atomBonds(at)) {
    for neighbor_1 in adjacency.neighbors_of(atom.index()) {
        let bond_1 = molecule.bonds().get(neighbor_1.bond.index()).ok_or(
            crate::ValenceError::UnsupportedBranch {
                reason: "conjugation adjacency bond index out of range",
            },
        )?;
        let other_1 = AtomId::new(neighbor_1.atom_index);
        // RDKit✔️✔️:     if (bnd1->getValenceContrib(at) < 1.5 ||
        // RDKit✔️✔️:         !isAtomConjugCand(bnd1->getOtherAtom(at))) { continue; }
        if crate::valence::bond_valence_contrib(bond_1, atom)? < 1.5
            || !sanitize_is_atom_conjugation_candidate(molecule, adjacency, valence, other_1)?
        {
            continue;
        }

        // RDKit✔️✔️:     for (const auto bnd2 : mol.atomBonds(at)) {
        for neighbor_2 in adjacency.neighbors_of(atom.index()) {
            if neighbor_1.bond == neighbor_2.bond {
                continue;
            }
            let atom_2 = AtomId::new(neighbor_2.atom_index);
            // RDKit✔️✔️:       sbo = at2->getDegree() + at2->getTotalNumHs();
            // RDKit✔️✔️:       if (sbo > 3) { continue; }
            if sanitize_total_degree(molecule, adjacency, valence, atom_2)? > 3 {
                continue;
            }
            // RDKit✔️✔️:       if (isAtomConjugCand(at2)) {
            // RDKit✔️✔️:         bnd1->setIsConjugated(true);
            // RDKit✔️✔️:         bnd2->setIsConjugated(true);
            // RDKit✔️✔️:       }
            if sanitize_is_atom_conjugation_candidate(molecule, adjacency, valence, atom_2)? {
                conjugated[neighbor_1.bond.index()] = true;
                conjugated[neighbor_2.bond.index()] = true;
            }
        }
    }
    Ok(())
    // END RDKIT CPP FUNCTION markConjAtomBonds
}

fn sanitize_is_atom_conjugation_candidate<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<bool, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION isAtomConjugCand
    // RDKit✔️✔️: bool isAtomConjugCand(const Atom *at) {
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "conjugation atom index out of range",
        })?;
    // RDKit✔️✔️:   const auto &vals =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
    let valence_list = crate::rdkit_valence_list(atom_data.atomic_number())?.ok_or(
        crate::ValenceError::UnsupportedBranch {
            reason: "conjugation valence-list atomic number out of range",
        },
    )?;
    let total_valence = sanitize_total_valence(valence, atom);
    // RDKit✔️✔️:   if (!at->getFormalCharge() && vals.front() >= 0 &&
    // RDKit✔️✔️:       at->getTotalValence() > static_cast<unsigned int>(vals.front())) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if atom_data.formal_charge() == 0
        && valence_list.first().is_some_and(|value| *value >= 0)
        && total_valence > valence_list[0]
    {
        return Ok(false);
    }
    // RDKit✔️✔️:   int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:   auto res = ((at->getAtomicNum() <= 10) || (nouter != 5 && nouter != 6) ||
    // RDKit✔️✔️:               (nouter == 6 && at->getTotalDegree() < 2u)) &&
    // RDKit✔️✔️:              MolOps::countAtomElec(at) > 0;
    let n_outer = crate::valence::periodic_table_outer_electrons(atom_data.atomic_number())?;
    Ok(((atom_data.atomic_number() <= 10)
        || (n_outer != 5 && n_outer != 6)
        || (n_outer == 6 && sanitize_total_degree(molecule, adjacency, valence, atom)? < 2))
        && sanitize_count_atom_electrons(molecule, adjacency, valence, atom)? > 0)
    // END RDKIT CPP FUNCTION isAtomConjugCand
}

fn sanitize_hybridization_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<crate::Hybridization>, crate::ValenceError> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION MolOps::setHybridization
    // RDKit✔️✔️: void setHybridization(ROMol &mol) {
    let adjacency = sanitize_adjacency(molecule)?;
    let valence = sanitize_valence_facts(molecule)?;
    molecule
        .atoms()
        .iter()
        .map(|atom| {
            // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
            // RDKit✔️✔️:     if (atom->getAtomicNum() == 0) {
            // RDKit✔️✔️:       atom->setHybridization(Atom::UNSPECIFIED);
            if atom.atomic_number() == 0 {
                return Ok(crate::Hybridization::Unspecified);
            }

            let total_degree = sanitize_total_degree(molecule, &adjacency, &valence, atom.id())?;
            // RDKit✔️✔️:       switch (atom->getChiralTag()) { ... stereo spec matches coordination number ... }
            match atom.chiral_tag() {
                crate::ChiralTag::Tetrahedral
                | crate::ChiralTag::TetrahedralCw
                | crate::ChiralTag::TetrahedralCcw
                    if total_degree == 4 =>
                {
                    return Ok(crate::Hybridization::Sp3);
                }
                crate::ChiralTag::SquarePlanar if (2..=4).contains(&total_degree) => {
                    return Ok(crate::Hybridization::Sp2d);
                }
                crate::ChiralTag::TrigonalBipyramidal if (2..=5).contains(&total_degree) => {
                    return Ok(crate::Hybridization::Sp3d);
                }
                crate::ChiralTag::Octahedral if (2..=6).contains(&total_degree) => {
                    return Ok(crate::Hybridization::Sp3d2);
                }
                _ => {}
            }

            // RDKit✔️✔️:       if (atom->getAtomicNum() < 89) {
            // RDKit✔️✔️:         norbs = numBondsPlusLonePairs(atom);
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         norbs = atom->getTotalDegree();
            // RDKit✔️✔️:       }
            let norbs = if atom.atomic_number() < 89 {
                sanitize_num_bonds_plus_lone_pairs(molecule, &adjacency, &valence, atom.id())?
            } else {
                total_degree
            };
            // RDKit✔️✔️:       switch (norbs) { ... }
            Ok(match norbs {
                0 | 1 => crate::Hybridization::S,
                2 => crate::Hybridization::Sp,
                3 => crate::Hybridization::Sp2,
                4 => {
                    if total_degree > 3 || !sanitize_atom_has_conjugated_bond(molecule, atom.id()) {
                        crate::Hybridization::Sp3
                    } else {
                        crate::Hybridization::Sp2
                    }
                }
                5 => crate::Hybridization::Sp3d,
                6 => crate::Hybridization::Sp3d2,
                _ => crate::Hybridization::Unspecified,
            })
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
        })
        .collect()
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::setHybridization
}

fn sanitize_num_bonds_plus_lone_pairs<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION numBondsPlusLonePairs
    // RDKit✔️✔️: int numBondsPlusLonePairs(Atom *at) {
    // RDKit✔️✔️:   int deg = at->getTotalDegree();
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "hybridization atom index out of range",
        })?;
    let mut degree = sanitize_total_degree(molecule, adjacency, valence, atom)?;
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::ZERO ||
    // RDKit✔️✔️:         (isDative(*bond) && at->getIdx() != bond->getEndAtomIdx())) { --deg; }
    // RDKit✔️✔️:   }
    for neighbor in adjacency.neighbors_of(atom.index()) {
        let bond = molecule.bonds().get(neighbor.bond.index()).ok_or(
            crate::ValenceError::UnsupportedBranch {
                reason: "hybridization adjacency bond index out of range",
            },
        )?;
        if bond.order() == crate::BondOrder::Zero
            || (sanitize_is_dative_bond(bond.order()) && atom != bond.end())
        {
            degree -= 1;
        }
    }
    // RDKit✔️✔️:   if (at->getAtomicNum() <= 1) { return deg; }
    if atom_data.atomic_number() <= 1 {
        return Ok(degree);
    }
    // RDKit✔️✔️:   int nouter = PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:   int totalValence = at->getTotalValence();
    // RDKit✔️✔️:   int chg = at->getFormalCharge();
    let n_outer = crate::valence::periodic_table_outer_electrons(atom_data.atomic_number())?;
    let total_valence = sanitize_total_valence(valence, atom);
    let charge = i32::from(atom_data.formal_charge());
    // RDKit✔️✔️:   int numFreeElectrons = nouter - (totalValence + chg);
    let free_electrons = n_outer - (total_valence + charge);
    // RDKit✔️✔️:   if (totalValence + nouter - chg < 8) { ... } else { ... }
    if total_valence + n_outer - charge < 8 {
        let radicals = i32::from(atom_data.radical_electrons());
        let lone_pairs = (free_electrons - radicals) / 2;
        Ok(degree + lone_pairs + radicals)
    } else {
        Ok(degree + free_electrons / 2)
    }
    // END RDKIT CPP FUNCTION numBondsPlusLonePairs
}

fn sanitize_count_atom_electrons<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "countAtomElec atom index out of range",
        })?;
    let default_valence = crate::valence::rdkit_default_valence(atom_data.atomic_number())?;
    if default_valence <= 1 {
        return Ok(-1);
    }
    let mut degree = sanitize_total_degree(molecule, adjacency, valence, atom)?;
    for neighbor in adjacency.neighbors_of(atom.index()) {
        let bond = molecule.bonds().get(neighbor.bond.index()).ok_or(
            crate::ValenceError::UnsupportedBranch {
                reason: "countAtomElec adjacency bond index out of range",
            },
        )?;
        if crate::valence::bond_valence_contrib(bond, atom)? == 0.0
            && !sanitize_is_bond_order_query(bond.query())
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return Ok(-1);
    }
    let n_outer = crate::valence::periodic_table_outer_electrons(atom_data.atomic_number())?;
    let lone_pairs = (n_outer - default_valence - i32::from(atom_data.formal_charge())).max(0);
    let radicals = i32::from(atom_data.radical_electrons());
    let mut result = (default_valence - degree) + lone_pairs - radicals;
    if result > 1 {
        let graph_degree =
            i32::try_from(adjacency.neighbors_of(atom.index()).len()).map_err(|_| {
                crate::ValenceError::UnsupportedBranch {
                    reason: "countAtomElec atom degree out of range",
                }
            })?;
        if valence.explicit_valence[atom.index()] - graph_degree > 1 {
            result = 1;
        }
    }
    Ok(result)
}

fn sanitize_valence_facts<'a>(
    molecule: impl MoleculeReadAccess<'a>,
) -> Result<crate::ValenceAssignment, crate::ValenceError> {
    if let Some(valence) = &molecule.derived_cache().valence {
        Ok(valence.clone())
    } else {
        molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
    }
}

fn sanitize_adjacency<'a>(
    molecule: impl MoleculeReadAccess<'a>,
) -> Result<crate::AdjacencyList, crate::ValenceError> {
    crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds()).map_err(|_| {
        crate::ValenceError::UnsupportedBranch {
            reason: "sanitize topology bond atom index out of range",
        }
    })
}

fn sanitize_total_degree<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    let atom_data = molecule
        .atom(atom)
        .ok_or(crate::ValenceError::UnsupportedBranch {
            reason: "sanitize total-degree atom index out of range",
        })?;
    let graph_degree = i32::try_from(adjacency.neighbors_of(atom.index()).len()).map_err(|_| {
        crate::ValenceError::UnsupportedBranch {
            reason: "sanitize atom degree out of range",
        }
    })?;
    Ok(graph_degree
        + i32::from(atom_data.explicit_hydrogens())
        + valence.implicit_hydrogens[atom.index()].max(0))
}

fn sanitize_total_valence(valence: &crate::ValenceAssignment, atom: AtomId) -> i32 {
    valence.explicit_valence[atom.index()] + valence.implicit_hydrogens[atom.index()].max(0)
}

fn sanitize_atom_has_conjugated_bond<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    atom: AtomId,
) -> bool {
    molecule
        .bonds()
        .iter()
        .any(|bond| bond.is_conjugated() && (bond.begin() == atom || bond.end() == atom))
}

fn sanitize_is_dative_bond(order: crate::BondOrder) -> bool {
    matches!(
        order,
        crate::BondOrder::Dative
            | crate::BondOrder::DativeOne
            | crate::BondOrder::DativeLeft
            | crate::BondOrder::DativeRight
    )
}

fn sanitize_is_bond_order_query(
    query: Option<&crate::QueryNode<crate::BondQueryPredicate>>,
) -> bool {
    match query {
        Some(query) => {
            crate::valence::has_bond_type_query(query)
                && !crate::valence::has_complex_bond_type_query(query)
        }
        None => false,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SanitizeCleanupAssignment {
    atom_formal_charges: Vec<i8>,
    bond_orders: Vec<crate::BondOrder>,
}

fn sanitize_cleanup_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<SanitizeCleanupAssignment, crate::ValenceError> {
    let molecule = read_parts;
    let adjacency = sanitize_adjacency(molecule)?;
    let mut assignment = SanitizeCleanupAssignment {
        atom_formal_charges: molecule
            .atoms()
            .iter()
            .map(crate::Atom::formal_charge)
            .collect(),
        bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
    };

    sanitize_nitrogens_cleanup_assignment(molecule, &adjacency, &mut assignment)?;
    for atom in molecule.atoms() {
        match atom.atomic_number() {
            15 => sanitize_phosphorus_cleanup_assignment(
                molecule,
                &adjacency,
                atom.id(),
                &mut assignment,
            )?,
            17 | 35 | 53 => sanitize_halogen_cleanup_assignment(
                molecule,
                &adjacency,
                atom.id(),
                &mut assignment,
            )?,
            _ => {}
        }
    }

    Ok(assignment)
}

fn sanitize_nitrogens_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    assignment: &mut SanitizeCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION nitrogensCleanup
    // RDKit✔️✔️: void nitrogensCleanup(RWMol &mol) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> nitrogensToConsider(mol.getNumAtoms());
    let mut nitrogens_to_consider = Vec::new();

    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->getAtomicNum() != 7) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->getFormalCharge()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->calcExplicitValence(false) != 5) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     nitrogensToConsider.set(atom->getIdx());
    for atom in molecule.atoms() {
        if atom.atomic_number() != 7 || assignment.atom_formal_charges[atom.id().index()] != 0 {
            continue;
        }
        if sanitize_cleanup_explicit_valence(molecule, adjacency, assignment, atom.id())? != 5 {
            continue;
        }
        nitrogens_to_consider.push(atom.id());

        // RDKit✔️✔️:     bool updateNeeded = false;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
        // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
        // RDKit✔️✔️:            Bond::DOUBLE)) {
        // RDKit✔️✔️:         auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        // RDKit✔️✔️:         b->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:         atom->setFormalCharge(1);
        // RDKit✔️✔️:         nbr->setFormalCharge(-1);
        // RDKit✔️✔️:         updateNeeded = true;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom.id()) {
            let bond = &molecule.bonds()[bond_index];
            let neighbor = sanitize_cleanup_other_atom(bond, atom.id());
            let neighbor_atom = &molecule.atoms()[neighbor.index()];
            if neighbor_atom.atomic_number() == 8
                && assignment.atom_formal_charges[neighbor.index()] == 0
                && assignment.bond_orders[bond_index] == crate::BondOrder::Double
            {
                assignment.bond_orders[bond_index] = crate::BondOrder::Single;
                assignment.atom_formal_charges[atom.id().index()] = 1;
                assignment.atom_formal_charges[neighbor.index()] = -1;
                break;
            }
        }
    }

    // RDKit✔️✔️:   for (auto aid = nitrogensToConsider.find_first();
    // RDKit✔️✔️:        aid != boost::dynamic_bitset<>::npos;
    // RDKit✔️✔️:        aid = nitrogensToConsider.find_next(aid)) {
    for atom in nitrogens_to_consider {
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 7) && (nbr->getFormalCharge() == 0) &&
        // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
        // RDKit✔️✔️:            Bond::TRIPLE)) {
        // RDKit✔️✔️:         auto b = mol.getBondBetweenAtoms(aid, nbr->getIdx());
        // RDKit✔️✔️:         b->setBondType(Bond::DOUBLE);
        // RDKit✔️✔️:         atom->setFormalCharge(1);
        // RDKit✔️✔️:         nbr->setFormalCharge(-1);
        // RDKit✔️✔️:         updateNeeded = true;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
            let bond = &molecule.bonds()[bond_index];
            let neighbor = sanitize_cleanup_other_atom(bond, atom);
            let neighbor_atom = &molecule.atoms()[neighbor.index()];
            if neighbor_atom.atomic_number() == 7
                && assignment.atom_formal_charges[neighbor.index()] == 0
                && assignment.bond_orders[bond_index] == crate::BondOrder::Triple
            {
                assignment.bond_orders[bond_index] = crate::BondOrder::Double;
                assignment.atom_formal_charges[atom.index()] = 1;
                assignment.atom_formal_charges[neighbor.index()] = -1;
                break;
            }
        }
    }
    // END RDKIT CPP FUNCTION nitrogensCleanup
    Ok(())
}

fn sanitize_phosphorus_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    atom: AtomId,
    assignment: &mut SanitizeCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION phosphorusCleanup
    // RDKit✔️✔️: void phosphorusCleanup(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   if (atom->getFormalCharge()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if assignment.atom_formal_charges[atom.index()] != 0 {
        return Ok(());
    }
    // RDKit✔️✔️:   if (atom->calcExplicitValence(false) == 5 && atom->getDegree() == 3) {
    if sanitize_cleanup_explicit_valence(molecule, adjacency, assignment, atom)? != 5
        || sanitize_cleanup_degree(adjacency, atom) != 3
    {
        return Ok(());
    }

    // RDKit✔️✔️:     Bond *dbl_to_O = nullptr;
    // RDKit✔️✔️:     Atom *O_atom = nullptr;
    // RDKit✔️✔️:     bool hasDoubleToCorN = false;
    let mut double_to_oxygen = None;
    let mut oxygen_atom = None;
    let mut has_double_to_carbon_or_nitrogen = false;
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
    for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
        let bond = &molecule.bonds()[bond_index];
        let neighbor = sanitize_cleanup_other_atom(bond, atom);
        let neighbor_atom = &molecule.atoms()[neighbor.index()];
        if neighbor_atom.atomic_number() == 8
            && assignment.atom_formal_charges[neighbor.index()] == 0
            && assignment.bond_orders[bond_index] == crate::BondOrder::Double
        {
            // RDKit✔️✔️:       if ((nbr->getAtomicNum() == 8) && (nbr->getFormalCharge() == 0) &&
            // RDKit✔️✔️:           (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:            Bond::DOUBLE)) {
            // RDKit✔️✔️:         dbl_to_O = mol.getBondBetweenAtoms(aid, nbr->getIdx());
            // RDKit✔️✔️:         O_atom = nbr;
            double_to_oxygen = Some(bond_index);
            oxygen_atom = Some(neighbor);
        } else if matches!(neighbor_atom.atomic_number(), 6 | 7)
            && sanitize_cleanup_degree(adjacency, neighbor) >= 2
            && assignment.bond_orders[bond_index] == crate::BondOrder::Double
        {
            // RDKit✔️✔️:       } else if ((nbr->getAtomicNum() == 6 || nbr->getAtomicNum() == 7) &&
            // RDKit✔️✔️:                  (nbr->getDegree() >= 2) &&
            // RDKit✔️✔️:                  (mol.getBondBetweenAtoms(aid, nbr->getIdx())->getBondType() ==
            // RDKit✔️✔️:                   Bond::DOUBLE)) {
            // RDKit✔️✔️:         hasDoubleToCorN = true;
            has_double_to_carbon_or_nitrogen = true;
        }
    }

    // RDKit✔️✔️:     if (hasDoubleToCorN && dbl_to_O != nullptr) {
    // RDKit✔️✔️:       O_atom->setFormalCharge(-1);
    // RDKit✔️✔️:       dbl_to_O->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:       atom->setFormalCharge(1);
    // RDKit✔️✔️:     }
    if has_double_to_carbon_or_nitrogen {
        if let (Some(bond_index), Some(oxygen)) = (double_to_oxygen, oxygen_atom) {
            assignment.atom_formal_charges[oxygen.index()] = -1;
            assignment.bond_orders[bond_index] = crate::BondOrder::Single;
            assignment.atom_formal_charges[atom.index()] = 1;
        }
    }
    // END RDKIT CPP FUNCTION phosphorusCleanup
    Ok(())
}

fn sanitize_halogen_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    atom: AtomId,
    assignment: &mut SanitizeCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION halogenCleanup
    // RDKit✔️✔️: void halogenCleanup(RWMol &mol, Atom *atom) {
    // RDKit✔️✔️:   int ev = atom->calcExplicitValence(false);
    // RDKit✔️✔️:   if (atom->getFormalCharge() == 0 && (ev == 7 || ev == 5 || ev == 3)) {
    let explicit_valence =
        sanitize_cleanup_explicit_valence(molecule, adjacency, assignment, atom)?;
    if assignment.atom_formal_charges[atom.index()] != 0 || !matches!(explicit_valence, 3 | 5 | 7) {
        return Ok(());
    }

    // RDKit✔️✔️:     bool neighborsAllO = true;
    // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
    // RDKit✔️✔️:       if (nbr->getAtomicNum() != 8) {
    // RDKit✔️✔️:         neighborsAllO = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    if !sanitize_cleanup_incident_bonds(adjacency, atom)
        .into_iter()
        .all(|bond_index| {
            let neighbor = sanitize_cleanup_other_atom(&molecule.bonds()[bond_index], atom);
            molecule.atoms()[neighbor.index()].atomic_number() == 8
        })
    {
        return Ok(());
    }

    // RDKit✔️✔️:     if (neighborsAllO) {
    // RDKit✔️✔️:       int formalCharge = 0;
    // RDKit✔️✔️:       for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:         if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:           auto otherAtom = bond->getOtherAtom(atom);
    // RDKit✔️✔️:           formalCharge++;
    // RDKit✔️✔️:           otherAtom->setFormalCharge(-1);
    // RDKit✔️✔️:           otherAtom->calcExplicitValence(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->setFormalCharge(formalCharge);
    let mut formal_charge = 0;
    for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
        if assignment.bond_orders[bond_index] == crate::BondOrder::Double {
            assignment.bond_orders[bond_index] = crate::BondOrder::Single;
            let other_atom = sanitize_cleanup_other_atom(&molecule.bonds()[bond_index], atom);
            formal_charge += 1;
            assignment.atom_formal_charges[other_atom.index()] = -1;
        }
    }
    assignment.atom_formal_charges[atom.index()] = formal_charge;
    // END RDKIT CPP FUNCTION halogenCleanup
    Ok(())
}

fn sanitize_cleanup_explicit_valence<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    assignment: &SanitizeCleanupAssignment,
    atom: AtomId,
) -> Result<i32, crate::ValenceError> {
    let mut accum = f64::from(molecule.atoms()[atom.index()].explicit_hydrogens());
    for bond_index in sanitize_cleanup_incident_bonds(adjacency, atom) {
        let mut bond = molecule.bonds()[bond_index].clone();
        bond.set_order(assignment.bond_orders[bond_index]);
        accum += crate::valence::bond_valence_contrib(&bond, atom)?;
    }
    Ok((accum + 0.1) as i32)
}

fn sanitize_cleanup_degree(adjacency: &crate::AdjacencyList, atom: AtomId) -> usize {
    adjacency.neighbors_of(atom.index()).len()
}

fn sanitize_cleanup_incident_bonds(adjacency: &crate::AdjacencyList, atom: AtomId) -> Vec<usize> {
    // RDKit✔️✔️: cleanUp() neighbors/bonds helper equivalent
    adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| neighbor.bond.index())
        .collect()
}

fn sanitize_cleanup_other_atom(bond: &crate::Bond, atom: AtomId) -> AtomId {
    if bond.begin() == atom {
        bond.end()
    } else {
        bond.begin()
    }
}

fn cleanup_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &SanitizeCleanupAssignment,
) -> bool {
    let molecule = read_parts;
    molecule
        .atoms()
        .iter()
        .zip(assignment.atom_formal_charges.iter().copied())
        .any(|(atom, formal_charge)| atom.formal_charge() != formal_charge)
        || molecule
            .bonds()
            .iter()
            .zip(assignment.bond_orders.iter().copied())
            .any(|(bond, order)| bond.order() != order)
}

fn apply_sanitize_cleanup_assignment(
    topology: &mut TopologyBlock,
    assignment: SanitizeCleanupAssignment,
) {
    for (atom, formal_charge) in topology
        .atoms
        .iter_mut()
        .zip(assignment.atom_formal_charges)
    {
        atom.set_formal_charge(formal_charge);
    }
    for (bond, order) in topology.bonds.iter_mut().zip(assignment.bond_orders) {
        bond.set_order(order);
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SanitizeOrganometallicCleanupAssignment {
    bond_orders: Vec<crate::BondOrder>,
    bond_endpoints: Vec<(AtomId, AtomId)>,
}

fn sanitize_organometallic_cleanup_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<SanitizeOrganometallicCleanupAssignment, crate::ValenceError> {
    let molecule = read_parts;
    let mut assignment = SanitizeOrganometallicCleanupAssignment {
        bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
        bond_endpoints: molecule
            .bonds()
            .iter()
            .map(|bond| (bond.begin(), bond.end()))
            .collect(),
    };

    // BEGIN RDKIT CPP FUNCTION MolOps::cleanUpOrganometallics
    // RDKit✔️✔️: void cleanUpOrganometallics(RWMol &mol) {
    // RDKit✔️✔️:   bool needsFixing = false;
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (isHypervalentNonMetal(atom) && !noDative(atom)) {
    let valence = molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?;
    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| crate::ValenceError::UnsupportedBranch {
            reason: "organometallic cleanup adjacency could not be built",
        })?;

    let mut needs_fixing = false;
    for atom in molecule.atoms() {
        if sanitize_is_hypervalent_nonmetal(molecule, &adjacency, &valence, atom.id())?
            && !sanitize_no_dative(atom)
        {
            // RDKit✔️✔️:       for (auto bond : mol.atomBonds(atom)) {
            // RDKit✔️✔️:         if (bond->getBondType() == Bond::BondType::SINGLE &&
            // RDKit✔️✔️:             QueryOps::isMetal(*bond->getOtherAtom(atom))) {
            // RDKit✔️✔️:           needsFixing = true;
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            needs_fixing = !sanitize_organometallic_single_bonded_metals(
                molecule,
                &adjacency,
                &assignment,
                atom.id(),
            )
            .is_empty();
        }
        if needs_fixing {
            break;
        }
    }
    // RDKit✔️✔️:   if (!needsFixing) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if !needs_fixing {
        return Ok(assignment);
    }

    // RDKit✔️✔️:   mol.updatePropertyCache(false);
    // RDKit✔️✔️:   std::vector<unsigned int> ranks(mol.getNumAtoms());
    // RDKit✔️✔️:   RDKit::Canon::rankMolAtoms(mol, ranks);
    let ranks = MoleculeReadAccess::rank_mol_atoms(molecule).map_err(|_| {
        crate::ValenceError::UnsupportedBranch {
            reason: "organometallic cleanup canonical ranking failed",
        }
    })?;
    // RDKit✔️✔️:   std::sort(atom_ranks.begin(), atom_ranks.end(),
    // RDKit✔️✔️:             [](const std::pair<int, int> &p1, std::pair<int, int> &p2) -> bool {
    // RDKit✔️✔️:               return p1.second < p2.second;
    // RDKit✔️✔️:             });
    let mut atom_order: Vec<usize> = (0..molecule.num_atoms()).collect();
    atom_order.sort_by_key(|atom| ranks[*atom]);
    // RDKit✔️✔️:   for (auto ar : atom_ranks) {
    // RDKit✔️✔️:     auto atom = mol.getAtomWithIdx(ar.first);
    // RDKit✔️✔️:     metalBondCleanup(mol, atom, ranks);
    // RDKit✔️✔️:   }
    for atom_index in atom_order {
        sanitize_metal_bond_cleanup_assignment(
            molecule,
            &adjacency,
            &valence,
            &ranks,
            AtomId::new(atom_index),
            &mut assignment,
        )?;
    }
    // END RDKIT CPP FUNCTION MolOps::cleanUpOrganometallics

    Ok(assignment)
}

fn sanitize_metal_bond_cleanup_assignment<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    ranks: &[usize],
    atom: AtomId,
    assignment: &mut SanitizeOrganometallicCleanupAssignment,
) -> Result<(), crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION metalBondCleanup
    // RDKit✔️✔️: void metalBondCleanup(RWMol &mol, Atom *atom,
    // RDKit✔️✔️:                       const std::vector<unsigned int> &ranks) {
    // RDKit✔️✔️:   if (isHypervalentNonMetal(atom) && !noDative(atom)) {
    if !sanitize_is_hypervalent_nonmetal(molecule, adjacency, valence, atom)?
        || sanitize_no_dative(&molecule.atoms()[atom.index()])
    {
        return Ok(());
    }

    // RDKit✔️✔️:     std::vector<Atom *> metals;
    // RDKit✔️✔️:     for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:           QueryOps::isMetal(*bond->getOtherAtom(atom))) {
    // RDKit✔️✔️:         metals.push_back(bond->getOtherAtom(atom));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    let mut metals =
        sanitize_organometallic_single_bonded_metals(molecule, adjacency, assignment, atom);
    if metals.is_empty() {
        return Ok(());
    }

    // RDKit✔️✔️:       std::sort(metals.begin(), metals.end(),
    // RDKit✔️✔️:                 [&](const Atom *a1, const Atom *a2) -> bool {
    // RDKit✔️✔️:                   int nda1 = numDativeBonds(a1);
    // RDKit✔️✔️:                   int nda2 = numDativeBonds(a2);
    // RDKit✔️✔️:                   if (nda1 == nda2) {
    // RDKit✔️✔️:                     return ranks[a1->getIdx()] > ranks[a2->getIdx()];
    // RDKit✔️✔️:                   } else {
    // RDKit✔️✔️:                     return nda1 < nda2;
    // RDKit✔️✔️:                   }
    // RDKit✔️✔️:                 });
    metals.sort_by(|left, right| {
        let left_dative = sanitize_num_dative_bonds(adjacency, assignment, *left);
        let right_dative = sanitize_num_dative_bonds(adjacency, assignment, *right);
        left_dative
            .cmp(&right_dative)
            .then_with(|| ranks[right.index()].cmp(&ranks[left.index()]))
    });
    let metal = metals[0];

    // RDKit✔️✔️:       auto bond =
    // RDKit✔️✔️:           mol.getBondBetweenAtoms(atom->getIdx(), metals.front()->getIdx());
    // RDKit✔️✔️:       if (bond) {
    // RDKit✔️✔️:         bond->setBondType(RDKit::Bond::BondType::DATIVE);
    // RDKit✔️✔️:         bond->setBeginAtom(atom);
    // RDKit✔️✔️:         bond->setEndAtom(metals.front());
    // RDKit✔️✔️:       }
    if let Some(bond_index) = adjacency
        .neighbors_of(atom.index())
        .iter()
        .find(|neighbor| neighbor.atom_index == metal.index())
        .map(|neighbor| neighbor.bond.index())
    {
        assignment.bond_orders[bond_index] = crate::BondOrder::Dative;
        assignment.bond_endpoints[bond_index] = (atom, metal);
    }
    // END RDKIT CPP FUNCTION metalBondCleanup
    Ok(())
}

fn sanitize_is_hypervalent_nonmetal<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<bool, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION isHypervalentNonMetal
    // RDKit✔️✔️: bool isHypervalentNonMetal(Atom *atom) {
    // RDKit✔️✔️:   if (QueryOps::isMetal(*atom)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    let atom_data = &molecule.atoms()[atom.index()];
    if is_metal_atomic_number(atom_data.atomic_number()) {
        return Ok(false);
    }
    // RDKit✔️✔️:   int ev = atom->getValence(Atom::ValenceType::EXPLICIT);
    let explicit_valence = valence.explicit_valence[atom.index()];
    // RDKit✔️✔️:   int effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
    let effective_atomic_number =
        i16::from(atom_data.atomic_number()) - i16::from(atom_data.formal_charge());
    // RDKit✔️✔️:   if (effAtomicNum <= 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if effective_atomic_number <= 0 {
        return Ok(false);
    }
    // RDKit✔️✔️:   const auto &otherValens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(effAtomicNum);
    let Some(valences) = crate::rdkit_valence_list(effective_atomic_number as u8)? else {
        return Ok(false);
    };
    let Some(&max_valence) = valences.last() else {
        return Ok(false);
    };
    // RDKit✔️✔️:   if (maxV > 0 && (ev > maxV || (ev == maxV && atom->getIsAromatic() &&
    // RDKit✔️✔️:                                  atom->getTotalDegree() == 4))) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if max_valence > 0
        && (explicit_valence > max_valence
            || (explicit_valence == max_valence
                && atom_data.is_aromatic()
                && sanitize_total_degree(molecule, adjacency, valence, atom)? == 4))
    {
        return Ok(true);
    }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isHypervalentNonMetal
    Ok(false)
}

fn sanitize_num_dative_bonds(
    adjacency: &crate::AdjacencyList,
    assignment: &SanitizeOrganometallicCleanupAssignment,
    atom: AtomId,
) -> usize {
    adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| sanitize_is_dative_bond(assignment.bond_orders[neighbor.bond.index()]))
        .count()
}

fn sanitize_no_dative(atom: &crate::Atom) -> bool {
    matches!(atom.atomic_number(), 1 | 2 | 9 | 10)
}

fn sanitize_organometallic_single_bonded_metals<'a>(
    molecule: impl MoleculeReadAccess<'a>,
    adjacency: &crate::AdjacencyList,
    assignment: &SanitizeOrganometallicCleanupAssignment,
    atom: AtomId,
) -> Vec<AtomId> {
    // RDKit✔️✔️: cleanUpOrganometallics()/metalBondCleanup() single-bonded metal helper
    adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter_map(|neighbor| {
            if assignment.bond_orders[neighbor.bond.index()] != crate::BondOrder::Single {
                return None;
            }
            let other = AtomId::new(neighbor.atom_index);
            is_metal_atomic_number(molecule.atoms()[other.index()].atomic_number()).then_some(other)
        })
        .collect()
}

fn organometallic_cleanup_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &SanitizeOrganometallicCleanupAssignment,
) -> bool {
    let molecule = read_parts;
    molecule
        .bonds()
        .iter()
        .zip(assignment.bond_orders.iter().copied())
        .any(|(bond, order)| bond.order() != order)
        || molecule
            .bonds()
            .iter()
            .zip(assignment.bond_endpoints.iter().copied())
            .any(|(bond, (begin, end))| bond.begin() != begin || bond.end() != end)
}

fn apply_sanitize_organometallic_cleanup_assignment(
    topology: &mut TopologyBlock,
    assignment: SanitizeOrganometallicCleanupAssignment,
) {
    for (bond, (order, (begin, end))) in topology.bonds.iter_mut().zip(
        assignment
            .bond_orders
            .into_iter()
            .zip(assignment.bond_endpoints),
    ) {
        bond.set_order(order);
        bond.set_endpoints(begin, end);
    }
}

fn sanitize_cleanup_atropisomers_assignment(
    read_parts: MoleculeReadParts<'_>,
    rings: &crate::RingInfo,
) -> Vec<(crate::BondId, crate::BondStereo)> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION MolOps::cleanupAtropisomers
    // RDKit✔️✔️: void cleanupAtropisomers(RWMol &mol) {
    // RDKit✔️✔️:   auto hybs = MolOps::Hybridizations(mol);
    // RDKit✔️✔️:   MolOps::cleanupAtropisomers(mol, hybs);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void cleanupAtropisomers(RWMol &mol, MolOps::Hybridizations &hybs) {
    let hybridization = molecule
        .atoms()
        .iter()
        .map(crate::Atom::hybridization)
        .collect::<Vec<_>>();
    let mut assignment = Vec::new();

    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     switch (bond->getStereo()) {
    // RDKit✔️✔️:       case Bond::BondStereo::STEREOATROPCW:
    // RDKit✔️✔️:       case Bond::BondStereo::STEREOATROPCCW:
    for bond in molecule.bonds() {
        if !matches!(
            bond.stereo(),
            crate::BondStereo::AtropCw | crate::BondStereo::AtropCcw
        ) {
            continue;
        }
        // BEGIN RDKIT CPP FUNCTION checkBond
        // RDKit✔️✔️:   if (hybs[bond->getBeginAtomIdx()] != Atom::SP2 ||
        // RDKit✔️✔️:       hybs[bond->getEndAtomIdx()] != Atom::SP2 ||
        // RDKit✔️✔️:       (ri->numBondRings(bond->getIdx()) > 0 &&
        // RDKit✔️✔️:        ri->minBondRingSize(bond->getIdx()) < 8)) {
        // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREONONE);
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION checkBond
        if hybridization[bond.begin().index()] != crate::Hybridization::Sp2
            || hybridization[bond.end().index()] != crate::Hybridization::Sp2
            || (rings.num_bond_rings(bond.id()) > 0 && rings.min_bond_ring_size(bond.id()) < 8)
        {
            assignment.push((bond.id(), crate::BondStereo::None));
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::cleanupAtropisomers
    assignment
}

fn apply_sanitize_cleanup_atropisomers_assignment(
    topology: &mut TopologyBlock,
    assignment: &[(crate::BondId, crate::BondStereo)],
) {
    for &(bond_id, stereo) in assignment {
        let bond = &mut topology.bonds[bond_id.index()];
        bond.set_stereo(stereo);
        if stereo == crate::BondStereo::None {
            bond.set_stereo_atoms(None);
            for stereo_group in &mut topology.stereo_groups {
                stereo_group.remove_bond(bond_id);
            }
        }
    }
    topology.stereo_groups.retain(|group| !group.is_empty());
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SanitizeCleanupChiralityAssignment {
    atom_chiral_tags: Vec<crate::ChiralTag>,
    chiral_permutations: Vec<Option<u32>>,
    cleanup_stereo_groups: bool,
}

fn sanitize_cleanup_chirality_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<SanitizeCleanupChiralityAssignment, crate::ValenceError> {
    let molecule = read_parts;
    let valence = molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?;
    let adjacency = crate::AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())
        .map_err(|_| crate::ValenceError::UnsupportedBranch {
            reason: "chirality cleanup adjacency could not be built",
        })?;
    let mut assignment = SanitizeCleanupChiralityAssignment {
        atom_chiral_tags: molecule
            .atoms()
            .iter()
            .map(crate::Atom::chiral_tag)
            .collect(),
        chiral_permutations: molecule
            .atoms()
            .iter()
            .map(crate::Atom::chiral_permutation)
            .collect(),
        cleanup_stereo_groups: false,
    };

    // BEGIN RDKIT CPP FUNCTION cleanupChirality
    // RDKit✔️✔️: void cleanupChirality(RWMol &mol) {
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     switch (atom->getChiralTag()) {
    for atom in molecule.atoms() {
        let atom_id = atom.id();
        match atom.chiral_tag() {
            // RDKit✔️✔️:       case Atom::CHI_TETRAHEDRAL_CW:
            // RDKit✔️✔️:       case Atom::CHI_TETRAHEDRAL_CCW:
            // RDKit✔️✔️:         if (atom->getHybridization() != Atom::SP3) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::TetrahedralCw | crate::ChiralTag::TetrahedralCcw => {
                if atom.hybridization() != crate::Hybridization::Sp3 {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                    assignment.cleanup_stereo_groups = true;
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_TETRAHEDRAL:
            // RDKit✔️✔️:         if (atom->getHybridization() != Atom::SP3) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           perm = 0;
            // RDKit✔️✔️:           atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
            // RDKit✔️✔️:           if (perm > 2) {
            // RDKit✔️✔️:             perm = 0;
            // RDKit✔️✔️:             atom->setProp(common_properties::_chiralPermutation, perm);
            crate::ChiralTag::Tetrahedral => {
                if atom.hybridization() != crate::Hybridization::Sp3 {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                    assignment.cleanup_stereo_groups = true;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 2);
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_SQUAREPLANAR:
            // RDKit✔️✔️:         degree = atom->getTotalDegree();
            // RDKit✔️✔️:         if (degree < 2 || degree > 4) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::SquarePlanar => {
                let degree = sanitize_total_degree(molecule, &adjacency, &valence, atom_id)?;
                if !(2..=4).contains(&degree) {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 3);
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_TRIGONALBIPYRAMIDAL:
            // RDKit✔️✔️:         degree = atom->getTotalDegree();
            // RDKit✔️✔️:         if (degree < 2 || degree > 5) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::TrigonalBipyramidal => {
                let degree = sanitize_total_degree(molecule, &adjacency, &valence, atom_id)?;
                if !(2..=5).contains(&degree) {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 20);
                }
            }
            // RDKit✔️✔️:       case Atom::CHI_OCTAHEDRAL:
            // RDKit✔️✔️:         degree = atom->getTotalDegree();
            // RDKit✔️✔️:         if (degree < 2 || degree > 6) {
            // RDKit✔️✔️:           atom->setChiralTag(Atom::CHI_UNSPECIFIED);
            crate::ChiralTag::Octahedral => {
                let degree = sanitize_total_degree(molecule, &adjacency, &valence, atom_id)?;
                if !(2..=6).contains(&degree) {
                    assignment.atom_chiral_tags[atom_id.index()] = crate::ChiralTag::Unspecified;
                } else {
                    sanitize_reset_chiral_permutation_above(&mut assignment, atom_id, 30);
                }
            }
            // RDKit✔️✔️:       default:
            // RDKit✔️✔️:         break;
            crate::ChiralTag::Unspecified | crate::ChiralTag::Other | crate::ChiralTag::Allene => {}
        }
    }
    // RDKit✔️✔️:   if (needCleanupStereoGroups) {
    // RDKit✔️✔️:     Chirality::cleanupStereoGroups(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION cleanupChirality
    Ok(assignment)
}

fn sanitize_reset_chiral_permutation_above(
    assignment: &mut SanitizeCleanupChiralityAssignment,
    atom: AtomId,
    limit: u32,
) {
    if assignment.chiral_permutations[atom.index()].is_some_and(|value| value > limit) {
        assignment.chiral_permutations[atom.index()] = Some(0);
    }
}

fn cleanup_chirality_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    assignment: &SanitizeCleanupChiralityAssignment,
) -> bool {
    let molecule = read_parts;
    molecule
        .atoms()
        .iter()
        .zip(assignment.atom_chiral_tags.iter().copied())
        .any(|(atom, tag)| atom.chiral_tag() != tag)
        || molecule
            .atoms()
            .iter()
            .zip(assignment.chiral_permutations.iter().copied())
            .any(|(atom, permutation)| atom.chiral_permutation() != permutation)
}

fn apply_sanitize_cleanup_chirality_assignment(
    topology: &mut TopologyBlock,
    assignment: &SanitizeCleanupChiralityAssignment,
) {
    let mut cleared_atoms = Vec::new();
    for (idx, (atom, tag)) in topology
        .atoms
        .iter_mut()
        .zip(assignment.atom_chiral_tags.iter().copied())
        .enumerate()
    {
        if atom.chiral_tag() != tag && tag == crate::ChiralTag::Unspecified {
            cleared_atoms.push(AtomId::new(idx));
        }
        atom.set_chiral_tag(tag);
    }
    for (atom, permutation) in topology
        .atoms
        .iter_mut()
        .zip(assignment.chiral_permutations.iter().copied())
    {
        atom.set_chiral_permutation(permutation);
    }
    if assignment.cleanup_stereo_groups {
        for atom in cleared_atoms {
            for stereo_group in &mut topology.stereo_groups {
                stereo_group.remove_atom(atom);
            }
        }
        topology.stereo_groups.retain(|group| !group.is_empty());
    }
}

fn sanitize_adjust_hydrogens_assignment(
    read_parts: MoleculeReadParts<'_>,
) -> Result<Vec<u8>, crate::ValenceError> {
    let molecule = read_parts;
    // BEGIN RDKIT CPP FUNCTION adjustHs
    // RDKit✔️✔️: void adjustHs(RWMol &mol) {
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     int origImplicitV = atom->getValence(Atom::ValenceType::IMPLICIT);
    // RDKit✔️✔️:     atom->calcExplicitValence(false);
    // RDKit✔️✔️:     int origExplicitV = atom->getNumExplicitHs();
    // RDKit✔️✔️:     int newImplicitV = atom->calcImplicitValence(false);
    let original_valence = sanitize_valence_facts(molecule)?;
    let current_valence =
        molecule.assign_valence_with_options(crate::ValenceModel::RdkitLike, false)?;
    let mut explicit_hydrogens = molecule
        .atoms()
        .iter()
        .map(crate::Atom::explicit_hydrogens)
        .collect::<Vec<_>>();

    for atom in molecule.atoms() {
        let original_implicit = original_valence.implicit_hydrogens[atom.id().index()];
        let new_implicit = current_valence.implicit_hydrogens[atom.id().index()];
        // RDKit✔️✔️:     if (newImplicitV < origImplicitV) {
        // RDKit✔️✔️:       atom->setNumExplicitHs(origExplicitV + (origImplicitV - newImplicitV));
        // RDKit✔️✔️:       atom->calcExplicitValence(false);
        // RDKit✔️✔️:     }
        if new_implicit < original_implicit {
            let delta = u8::try_from(original_implicit - new_implicit).map_err(|_| {
                crate::ValenceError::UnsupportedBranch {
                    reason: "adjustHs implicit hydrogen delta out of range",
                }
            })?;
            explicit_hydrogens[atom.id().index()] =
                explicit_hydrogens[atom.id().index()].saturating_add(delta);
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION adjustHs
    Ok(explicit_hydrogens)
}

fn adjust_hydrogens_changes_molecule(
    read_parts: MoleculeReadParts<'_>,
    explicit_hydrogens: &[u8],
) -> bool {
    let molecule = read_parts;
    molecule
        .atoms()
        .iter()
        .zip(explicit_hydrogens.iter().copied())
        .any(|(atom, explicit_hydrogens)| atom.explicit_hydrogens() != explicit_hydrogens)
}

fn apply_sanitize_adjust_hydrogens_assignment(
    topology: &mut TopologyBlock,
    explicit_hydrogens: &[u8],
) {
    for (atom, explicit_hydrogens) in topology.atoms.iter_mut().zip(explicit_hydrogens.iter()) {
        atom.set_explicit_hydrogens(*explicit_hydrogens);
    }
}

fn is_metal_atomic_number(atomic_number: u8) -> bool {
    matches!(
        atomic_number,
        3 | 4
            | 11
            | 12
            | 13
            | 19..=31
            | 37..=50
            | 55..=83
            | 87..=118
    )
}

#[mol_op_body(assigned_valence, parts)]
fn assigned_valence_impl() -> Result<OpOutcome, OperationError> {
    let read = parts.begin_topology_read()?;
    let valence = read
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
        .map_err(|source| OperationError::Valence {
            operation: &ASSIGNED_VALENCE_SPEC,
            source,
        })?;
    parts.set_valence_cache(valence);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_rings, parts)]
fn assigned_rings_impl() -> Result<OpOutcome, OperationError> {
    let read = parts.begin_topology_read()?;
    let rings = read
        .symmetrize_sssr()
        .map_err(|source| OperationError::RingFinding {
            operation: &ASSIGNED_RINGS_SPEC,
            source,
        })?;
    parts.set_rings_cache(rings);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_ring_families, parts)]
fn assigned_ring_families_impl() -> Result<OpOutcome, OperationError> {
    let read = parts.begin_topology_read()?;
    let ring_families =
        read.find_ring_families(false, false)
            .map_err(|source| OperationError::RingFinding {
                operation: &ASSIGNED_RING_FAMILIES_SPEC,
                source,
            })?;
    parts.set_ring_families_cache(ring_families);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_aromaticity, parts)]
fn assigned_aromaticity_impl() -> Result<OpOutcome, OperationError> {
    let mut topology = parts.begin_topology_mut()?;
    let view = parts.read_parts_for_topology(topology.clone())?;
    let read = MoleculeReadParts::from_molecule(&view);
    let rings = read
        .symmetrize_sssr()
        .map_err(|source| OperationError::RingFinding {
            operation: &ASSIGNED_AROMATICITY_SPEC,
            source,
        })?;
    parts.set_rings_cache(rings);
    let view = parts.read_parts_for_topology(topology.clone())?;
    let assignment = MoleculeReadParts::from_molecule(&view)
        .set_aromaticity(crate::AromaticityModel::Default)
        .map_err(|source| OperationError::Aromaticity {
            operation: &ASSIGNED_AROMATICITY_SPEC,
            source,
        })?;
    for (atom, is_aromatic) in topology
        .atoms
        .iter_mut()
        .zip(assignment.atom_aromatic.iter().copied())
    {
        atom.set_aromatic(is_aromatic);
    }
    for (bond, is_aromatic) in topology
        .bonds
        .iter_mut()
        .zip(assignment.bond_aromatic.iter().copied())
    {
        bond.set_aromatic(is_aromatic);
        if is_aromatic
            && matches!(
                bond.order(),
                crate::BondOrder::Single | crate::BondOrder::Double
            )
        {
            bond.set_order(crate::BondOrder::Aromatic);
        }
    }
    let view = parts.read_parts_for_topology(topology.clone())?;
    let valence = MoleculeReadParts::from_molecule(&view)
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
        .map_err(|source| OperationError::Valence {
            operation: &ASSIGNED_AROMATICITY_SPEC,
            source,
        })?;
    parts.commit_topology(topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.set_valence_cache(valence);
    parts.mark_aromaticity_valid();
    parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_radicals, parts)]
fn assigned_radicals_impl() -> Result<OpOutcome, OperationError> {
    let mut topology = parts.begin_topology_mut()?;
    let view = parts.read_parts_for_topology(topology.clone())?;
    let read = MoleculeReadParts::from_molecule(&view);

    let radicals = read
        .assign_radicals()
        .map_err(|source| OperationError::Valence {
            operation: &ASSIGNED_RADICALS_SPEC,
            source,
        })?;

    let changed = read
        .atoms()
        .iter()
        .zip(radicals.iter().copied())
        .any(|(atom, radical)| atom.radical_electrons() != radical);

    if changed {
        for (atom, radical) in topology.atoms.iter_mut().zip(radicals) {
            atom.set_radical_electrons(radical);
        }
    }

    let view = parts.read_parts_for_topology(topology.clone())?;
    let valence = MoleculeReadParts::from_molecule(&view)
        .assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
        .map_err(|source| OperationError::Valence {
            operation: &ASSIGNED_RADICALS_SPEC,
            source,
        })?;
    parts.commit_topology(topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.set_valence_cache(valence);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(with_2d_coordinates, parts)]
fn with_2d_coordinates_impl() -> Result<OpOutcome, OperationError> {
    let (atoms, bonds) = {
        let read = parts.begin_topology_read()?;
        (read.atoms(), read.bonds())
    };
    let coords =
        crate::coordinates::compute_2d_coords(atoms, bonds).map_err(|source| match source {
            crate::coordinates::Coordinate2DError::InvalidInput(message) => {
                OperationError::InvalidInput {
                    operation: &WITH_2D_COORDINATES_SPEC,
                    message,
                }
            }
            crate::coordinates::Coordinate2DError::UnsupportedFeature(_) => {
                OperationError::UnsupportedFeature {
                    operation: &WITH_2D_COORDINATES_SPEC,
                    source: crate::UnsupportedFeatureError::from_spec(
                        &crate::COORDINATE_2D_FEATURE,
                    ),
                }
            }
        })?;
    let mut coord_block = parts.begin_coordinates_mut()?;
    coord_block.coords_2d = Some(coords);
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(OpOutcome::Changed)
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::{BondOrder, BondQueryPredicate, QueryNode};

    const TEST_NEEDS_VALENCE_UPDATE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_needs_valence_update",
        impl_fn: "test_needs_valence_update_impl",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        access: BlockAccess::new(BlockSet::NONE, BlockSet::DERIVED_CACHE),
        may_mutate: BlockSet::DERIVED_CACHE,
        auto_remap: BlockSet::NONE,
        derived_effects: DerivedEffects::new(
            DerivedState::NONE,                               // recompute
            DerivedState::NONE,                               // preserve
            DerivedState::VALENCE.union(DerivedState::RINGS), // invalidate: needs_update target + clear-permitted
        ),
        requires_mapping: MappingRequirement::None,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_RECOMPUTE_VALENCE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_recompute_valence",
        impl_fn: "test_recompute_valence_impl",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        access: BlockAccess::new(BlockSet::NONE, BlockSet::DERIVED_CACHE),
        may_mutate: BlockSet::DERIVED_CACHE,
        auto_remap: BlockSet::NONE,
        derived_effects: DerivedEffects::new(
            DerivedState::VALENCE, // recompute
            DerivedState::NONE,    // preserve
            DerivedState::NONE,    // invalidate
        ),
        requires_mapping: MappingRequirement::None,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_OVERLAPPING_ACCESS_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_overlapping_access",
        impl_fn: "test_overlapping_access_impl",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::Local,
        access: BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::TOPOLOGY),
        may_mutate: BlockSet::TOPOLOGY,
        auto_remap: BlockSet::NONE,
        derived_effects: DerivedEffects::NONE,
        requires_mapping: MappingRequirement::None,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    #[test]
    fn molecule_read_parts_does_not_expose_raw_molecule_escape() {
        let ops_source = include_str!("ops.rs");
        let read_parts_source = include_str!("../model/read_parts.rs");
        assert!(!ops_source.contains(concat!("read_parts", ".", "molecule")));
        assert!(!ops_source.contains(concat!(".", "molecule", "()")));
        assert!(!read_parts_source.contains(concat!("pub(crate) fn ", "molecule")));
    }

    #[test]
    fn begin_topology_mut_rejects_second_begin() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);
        let _topology = parts
            .begin_topology_mut()
            .expect("first topology begin should succeed");
        let err = match parts.begin_topology_mut() {
            Ok(_) => panic!("second topology begin must be rejected"),
            Err(err) => err,
        };
        assert!(
            matches!(err, OperationError::InvalidInput { message, .. } if message.contains("same writable block twice"))
        );
    }

    #[test]
    fn begin_topology_mut_rejects_second_begin_before_commit() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);
        let _topology = parts
            .begin_topology_mut()
            .expect("first topology begin should succeed");
        let err = match parts.begin_topology_mut() {
            Ok(_) => panic!("second topology begin must be rejected"),
            Err(err) => err,
        };
        assert!(
            matches!(err, OperationError::InvalidInput { message, .. } if message.contains("same writable block twice"))
        );
    }

    #[test]
    fn begin_access_rejects_overlapping_read_and_write_blocks() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_OVERLAPPING_ACCESS_SPEC);
        let err = match parts.begin_topology_mut() {
            Ok(_) => panic!("overlapping read/write access must be rejected"),
            Err(err) => err,
        };
        assert!(
            matches!(err, OperationError::InvalidInput { message, .. } if message.contains("both read and write"))
        );
    }

    #[test]
    fn sanitize_is_bond_order_query_matches_rdkit_simple_vs_complex_split() {
        assert!(sanitize_is_bond_order_query(Some(&QueryNode::predicate(
            BondQueryPredicate::Order(BondOrder::Single),
        ))));
        assert!(sanitize_is_bond_order_query(Some(&QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        ]))));
        assert!(!sanitize_is_bond_order_query(Some(&QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Double]),
        ))));
        assert!(!sanitize_is_bond_order_query(Some(&QueryNode::not(
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        ))));
        assert!(!sanitize_is_bond_order_query(Some(&QueryNode::predicate(
            BondQueryPredicate::Any,
        ))));
    }

    fn same_operation(left: &'static MoleculeOpSpec, right: &'static MoleculeOpSpec) -> bool {
        left.method == right.method && left.impl_fn == right.impl_fn
    }

    fn support_matrix_contains(operation: &'static MoleculeOpSpec) -> bool {
        SUPPORT_MATRIX.iter().any(|entry| {
            entry
                .operation
                .is_some_and(|candidate| same_operation(candidate, operation))
        })
    }

    fn invariant_matrix_contains(operation: &'static MoleculeOpSpec) -> bool {
        OPERATION_INVARIANT_MATRIX
            .iter()
            .any(|entry| same_operation(entry.operation, operation))
    }

    fn parity_matrix_contains(operation: &'static MoleculeOpSpec) -> bool {
        PARITY_MATRIX
            .iter()
            .any(|entry| same_operation(entry.operation, operation))
    }

    #[allow(dead_code)]
    fn assert_unsupported_feature(
        result: Result<crate::Molecule, OperationError>,
        operation: &'static MoleculeOpSpec,
        feature: &'static crate::FeatureSpec,
    ) {
        match result {
            Err(OperationError::UnsupportedFeature {
                operation: actual_operation,
                source,
            }) => {
                assert!(same_operation(actual_operation, operation));
                assert_eq!(source.feature, feature.name);
            }
            other => panic!(
                "expected UnsupportedFeature for {}, got {other:?}",
                operation.method
            ),
        }
    }

    #[test]
    fn registered_ops_have_unique_methods() {
        let mut methods = HashSet::new();
        for operation in MOLECULE_OPS.iter().copied() {
            assert!(
                methods.insert(operation.method),
                "duplicate registered operation method `{}`",
                operation.method
            );
        }
    }

    #[test]
    fn registered_ops_have_support_and_invariant_entries() {
        assert_eq!(SUPPORT_MATRIX.len(), MOLECULE_OPS.len());
        assert_eq!(OPERATION_INVARIANT_MATRIX.len(), MOLECULE_OPS.len());
        for operation in MOLECULE_OPS.iter().copied() {
            assert!(
                support_matrix_contains(operation),
                "missing support matrix entry for `{}`",
                operation.method
            );
            assert!(
                invariant_matrix_contains(operation),
                "missing invariant matrix entry for `{}`",
                operation.method
            );
        }
    }

    #[test]
    fn parity_registered_ops_have_parity_entries() {
        for operation in MOLECULE_OPS.iter().copied() {
            if operation.parity != ParityPolicy::NotApplicable {
                assert!(
                    parity_matrix_contains(operation),
                    "missing parity matrix entry for `{}`",
                    operation.method
                );
            }
        }
    }

    #[test]
    fn sanitized_all_runs_through_operation_pipeline_without_changing_source() {
        let molecule = crate::Molecule::new();
        let original = molecule.clone();

        let sanitized = molecule.sanitized().unwrap();
        let sanitized_with_all = molecule
            .sanitized_with_ops(crate::SanitizeOps::ALL)
            .unwrap();

        assert_eq!(sanitized.num_atoms(), 0);
        assert_eq!(sanitized_with_all.num_atoms(), 0);
        assert!(matches!(
            molecule.with_2d_coordinates(),
            Err(OperationError::InvalidInput {
                operation: &WITH_2D_COORDINATES_SPEC,
                ..
            })
        ));

        assert_eq!(molecule, original);
    }

    #[test]
    fn with_hydrogens_adds_implicit_hydrogens_through_operation_pipeline() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder.set_2d_coordinates(vec![[0.5, 1.0]]).unwrap();
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let result = molecule.with_hydrogens().unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.num_atoms(), 5);
        assert_eq!(result.num_bonds(), 4);
        assert_eq!(result.atoms()[carbon.index()].explicit_hydrogens(), 0);
        assert!(
            result.atoms()[1..]
                .iter()
                .all(|atom| atom.atomic_number() == 1 && atom.implicit_hydrogen())
        );
        assert_eq!(
            result.coords_2d(),
            Some(&[[0.5, 1.0], [0.5, 1.0], [0.5, 1.0], [0.5, 1.0], [0.5, 1.0]][..])
        );
    }

    #[test]
    fn add_hs_terminal_coords_follow_rdkit_sequential_append() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
        let molecule = builder.build().unwrap();
        let params = crate::hydrogens::AddHsParams {
            add_coords: true,
            ..crate::hydrogens::AddHsParams::default()
        };
        let read_parts = MoleculeReadParts::from_molecule(&molecule);
        let assignment = crate::hydrogens::add_hs_assignment(read_parts, &params).unwrap();

        let coords = add_hs_terminal_coords_2d(
            MoleculeReadParts::from_molecule(&molecule),
            &assignment,
            molecule.coords_2d().unwrap(),
        )
        .unwrap();

        assert_eq!(coords.len(), 4);
        assert_eq!(coords[0], [1.0, 0.0]);
        assert_ne!(coords[1], [0.0, 0.0]);
    }

    #[test]
    fn add_hs_terminal_coords_3d_uses_rdkit_rb0_bond_length() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder.add_3d_conformer(vec![[0.0, 0.0, 0.0]]).unwrap();
        let molecule = builder.build().unwrap();
        let assignment = crate::hydrogens::AddHsAssignment {
            hydrogens_to_add: vec![crate::hydrogens::AddHydrogen {
                heavy_atom: carbon,
                isotope: None,
                is_implicit: true,
                props: Default::default(),
                pdb_residue_info: None,
            }],
            add_terminal_coordinates: true,
            ..crate::hydrogens::AddHsAssignment::default()
        };

        let coords = add_hs_terminal_coords_3d(
            MoleculeReadParts::from_molecule(&molecule),
            &assignment,
            molecule.conformers_3d()[0].coords(),
        )
        .unwrap();

        assert_eq!(coords.len(), 1);
        assert!((coords[0][0] - 0.0).abs() < 1.0e-12);
        assert!((coords[0][1] - 0.0).abs() < 1.0e-12);
        assert!((coords[0][2] - 1.10).abs() < 1.0e-12);
    }

    #[test]
    fn add_hs_terminal_coords_default_degree_branch_matches_rdkit_zero_vector() {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let mut coords = vec![[5.0, 5.0, 5.0]];
        for index in 0..5 {
            let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(
                    center,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
            coords.push([index as f64, 0.0, 0.0]);
        }
        builder.add_3d_conformer(coords).unwrap();
        let molecule = builder.build().unwrap();
        let assignment = crate::hydrogens::AddHsAssignment {
            hydrogens_to_add: vec![crate::hydrogens::AddHydrogen {
                heavy_atom: center,
                isotope: None,
                is_implicit: true,
                props: Default::default(),
                pdb_residue_info: None,
            }],
            add_terminal_coordinates: true,
            ..crate::hydrogens::AddHsAssignment::default()
        };

        let coords = add_hs_terminal_coords_3d(
            MoleculeReadParts::from_molecule(&molecule),
            &assignment,
            molecule.conformers_3d()[0].coords(),
        )
        .unwrap();

        assert_eq!(coords, vec![[0.0, 0.0, 0.0]]);
    }

    #[test]
    fn add_hs_terminal_coords_rejects_degenerate_or_nonterminal_virtual_atom_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_3d_conformer(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            .unwrap();
        let molecule = builder.build().unwrap();
        let adjacency = vec![vec![(1, None), (2, None)], vec![(0, None)], vec![(0, None)]];

        let degenerate = add_hs_set_terminal_atom_coord(
            MoleculeReadParts::from_molecule(&molecule),
            &adjacency,
            molecule.conformers_3d()[0].coords(),
            0,
            0,
            true,
        )
        .unwrap_err();
        assert!(degenerate.to_string().contains("degenerate atoms"));

        let nonterminal = add_hs_set_terminal_atom_coord(
            MoleculeReadParts::from_molecule(&molecule),
            &adjacency,
            molecule.conformers_3d()[0].coords(),
            0,
            1,
            true,
        )
        .unwrap_err();
        assert!(nonterminal.to_string().contains("degree one"));
    }

    #[test]
    fn with_hydrogens_with_params_materializes_add_coords_branch() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
        let molecule = builder.build().unwrap();
        let params = crate::AddHsParams {
            add_coords: true,
            ..crate::AddHsParams::default()
        };

        let result = molecule.with_hydrogens_with_params(params).unwrap();

        assert_eq!(result.num_atoms(), 5);
        let coords = result.coords_2d().unwrap();
        assert_eq!(coords.len(), 5);
        assert_eq!(coords[0], [0.0, 0.0]);
        assert_eq!(coords[1], [1.0, 0.0]);
        assert_eq!(coords[2], [-1.0, 0.0]);
    }

    #[test]
    fn with_hydrogens_materializes_explicit_h_count_and_clears_heavy_atom_count() {
        let mut builder = crate::MoleculeBuilder::new();
        let nitrogen =
            builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_explicit_hydrogens(2));
        let molecule = builder.build().unwrap();

        let result = molecule.with_hydrogens().unwrap();

        assert_eq!(result.num_atoms(), 4);
        assert_eq!(result.atoms()[nitrogen.index()].explicit_hydrogens(), 0);
        assert_eq!(
            result.atoms()[1..]
                .iter()
                .filter(|atom| !atom.implicit_hydrogen())
                .count(),
            2
        );
        assert_eq!(
            result.atoms()[1..]
                .iter()
                .filter(|atom| atom.implicit_hydrogen())
                .count(),
            1
        );
    }

    #[test]
    fn with_hydrogens_commits_topology_with_rebuilt_adjacency_for_valence_followups() {
        let molecule = crate::Molecule::from_smiles("C=C").unwrap();

        let result = molecule.with_hydrogens().unwrap();
        let assignment = crate::assign_valence(&result, crate::ValenceModel::RdkitLike).unwrap();

        assert_eq!(assignment.explicit_valence, vec![4, 4, 1, 1, 1, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![0, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn with_hydrogens_replays_tracked_isotopes_and_clears_tracking_property() {
        let mut builder = crate::MoleculeBuilder::new();
        let nitrogen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_explicit_hydrogens(2)
                .with_tracked_isotopic_hydrogens(vec![2, 3]),
        );
        let molecule = builder.build().unwrap();

        let result = molecule.with_hydrogens().unwrap();

        assert_eq!(result.atoms()[nitrogen.index()].prop("_isotopicHs"), None);
        assert_eq!(
            result.atoms()[nitrogen.index()].tracked_isotopic_hydrogens(),
            &[] as &[u16]
        );
        assert_eq!(result.atoms()[nitrogen.index()].explicit_hydrogens(), 0);
        assert_eq!(
            result.atoms()[1..]
                .iter()
                .map(crate::Atom::isotope)
                .collect::<Vec<_>>(),
            vec![Some(2), Some(3), None]
        );
    }

    #[test]
    fn add_hs_operation_materializes_typed_pdb_residue_info() {
        let mut builder = crate::MoleculeBuilder::new();
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let molecule = builder.build().unwrap();
        let assignment = crate::hydrogens::AddHsAssignment {
            hydrogens_to_add: vec![crate::hydrogens::AddHydrogen {
                heavy_atom: nitrogen,
                isotope: None,
                is_implicit: false,
                props: Default::default(),
                pdb_residue_info: Some(crate::AtomPdbResidueInfo::new(
                    " H1 ", 12, "GLY", 3, "A", false,
                )),
            }],
            ..crate::hydrogens::AddHsAssignment::default()
        };
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC);
        let mut topology = parts.begin_topology_mut().unwrap();
        let mut coordinates = parts.begin_coordinates_mut().unwrap();
        let mut properties = parts.begin_properties_mut().unwrap();

        let changed = apply_add_hs_assignment(
            &mut parts,
            &mut topology,
            &mut coordinates,
            &mut properties,
            &assignment,
        )
        .unwrap();
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();
        parts
            .prove_preserved(
                DerivedState::RINGS | DerivedState::RING_FAMILIES,
                PreservationProof::LeafAtomAppend,
            )
            .unwrap();
        let result = parts.finish(OpOutcome::Changed).unwrap();

        assert!(changed);
        let info = result.atoms()[1].pdb_residue_info().unwrap();
        assert_eq!(info.atom_name(), " H1 ");
        assert_eq!(info.serial_number(), 12);
        assert_eq!(info.residue_name(), "GLY");
        assert_eq!(info.residue_number(), 3);
        assert_eq!(info.chain_id(), "A");
    }

    #[test]
    fn with_hydrogens_with_params_materializes_add_residue_info_branch() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(
            crate::AtomSpec::new(crate::Element::N).with_pdb_residue_info(
                crate::AtomPdbResidueInfo::new(" N  ", 10, "GLY", 3, "A", false),
            ),
        );
        let molecule = builder.build().unwrap();
        let params = crate::AddHsParams {
            add_residue_info: true,
            ..crate::AddHsParams::default()
        };

        let result = molecule.with_hydrogens_with_params(params).unwrap();

        assert_eq!(result.num_atoms(), 4);
        let first_h_info = result.atoms()[1].pdb_residue_info().unwrap();
        assert_eq!(first_h_info.atom_name(), " H1 ");
        assert_eq!(first_h_info.residue_name(), "GLY");
        assert_eq!(first_h_info.residue_number(), 3);
        assert_eq!(first_h_info.chain_id(), "A");
    }

    #[test]
    fn add_hs_operation_materializes_existing_atom_pdb_residue_info_updates() {
        let mut builder = crate::MoleculeBuilder::new();
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let molecule = builder.build().unwrap();
        let assignment = crate::hydrogens::AddHsAssignment {
            atom_pdb_residue_info_updates: vec![crate::hydrogens::AtomPdbResidueInfoUpdate {
                atom: hydrogen,
                pdb_residue_info: crate::AtomPdbResidueInfo::new(" H1 ", 12, "GLY", 3, "A", false),
            }],
            ..crate::hydrogens::AddHsAssignment::default()
        };
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC);
        let mut topology = parts.begin_topology_mut().unwrap();
        let mut coordinates = parts.begin_coordinates_mut().unwrap();
        let mut properties = parts.begin_properties_mut().unwrap();

        let changed = apply_add_hs_assignment(
            &mut parts,
            &mut topology,
            &mut coordinates,
            &mut properties,
            &assignment,
        )
        .unwrap();
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();
        parts
            .prove_preserved(
                DerivedState::RINGS | DerivedState::RING_FAMILIES,
                PreservationProof::LeafAtomAppend,
            )
            .unwrap();
        let result = parts.finish(OpOutcome::Changed).unwrap();

        assert!(changed);
        let info = result.atoms()[hydrogen.index()].pdb_residue_info().unwrap();
        assert_eq!(info.atom_name(), " H1 ");
        assert_eq!(info.serial_number(), 12);
        assert_eq!(info.residue_name(), "GLY");
        assert_eq!(info.residue_number(), 3);
        assert_eq!(info.chain_id(), "A");
    }

    #[test]
    fn without_hydrogens_removes_basic_explicit_hydrogen_through_operation_pipeline() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_pdb_residue_info(
                crate::AtomPdbResidueInfo::new(" C  ", 7, "GLY", 3, "A", false),
            ),
        );
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        builder
            .add_bond(crate::BondSpec::new(
                carbon,
                hydrogen,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let result = molecule.without_hydrogens_with_sanitize(false).unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.num_bonds(), 0);
        assert_eq!(result.coords_2d(), Some(&[[0.0, 0.0]][..]));
        assert_eq!(
            result.atoms()[0]
                .pdb_residue_info()
                .unwrap()
                .serial_number(),
            7
        );
    }

    #[test]
    fn without_hydrogens_materializes_unknown_stereo_as_typed_atom_state() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        builder
            .add_bond(
                crate::BondSpec::new(carbon, hydrogen, crate::BondOrder::Single)
                    .with_direction(crate::BondDirection::Unknown),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.without_hydrogens_with_sanitize(false).unwrap();

        assert_eq!(result.num_atoms(), 1);
        assert!(result.atoms()[0].unknown_stereo());
        assert_eq!(result.atoms()[0].prop("_UnknownStereo"), None);
    }

    #[test]
    fn without_hydrogens_with_params_materializes_remove_and_track_isotopes_branch() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_no_implicit(true));
        let protium = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let deuterium = builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        builder
            .add_bond(crate::BondSpec::new(
                carbon,
                protium,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                carbon,
                deuterium,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = crate::RemoveHsParams {
            remove_and_track_isotopes: true,
            ..crate::RemoveHsParams::default()
        };

        let result = molecule
            .without_hydrogens_with_params(params, false)
            .unwrap();

        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.atoms()[0].tracked_isotopic_hydrogens(), &[2]);
        assert_eq!(result.atoms()[0].prop("_isotopicHs"), None);
    }

    #[test]
    fn without_hydrogens_updates_sgroup_membership_before_topology_compaction() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        builder
            .add_bond(crate::BondSpec::new(
                carbon,
                hydrogen,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::Superatom,
                )
                .with_atoms(vec![carbon, hydrogen]),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.without_hydrogens_with_sanitize(false).unwrap();

        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.substance_groups().len(), 1);
        assert_eq!(
            result.substance_groups()[0].atoms(),
            &[crate::AtomId::new(0)]
        );
    }

    #[test]
    fn without_hydrogens_with_sanitize_runs_full_sanitize_on_aromatic_result() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false)
            .unwrap()
            .with_hydrogens()
            .unwrap();
        let original = molecule.clone();

        let result = molecule.without_hydrogens_with_sanitize(true).unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.num_atoms(), 6);
        assert!(result.atoms().iter().all(crate::Atom::is_aromatic));
        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| bond.is_aromatic() && bond.order() == crate::BondOrder::Aromatic)
        );
        assert!(result.derived_cache().valence.is_some());
        assert!(result.derived_cache().rings.is_some());
    }

    #[test]
    fn without_hydrogens_without_sanitize_skips_full_sanitize_pipeline() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false)
            .unwrap()
            .with_hydrogens()
            .unwrap();

        let result = molecule.without_hydrogens_with_sanitize(false).unwrap();

        assert_eq!(result.num_atoms(), 6);
        assert!(result.derived_cache().valence.is_none());
        assert!(result.derived_cache().rings.is_none());
    }

    #[test]
    fn with_hydrogens_preserves_existing_ring_caches_by_leaf_append_proof() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("C1CCCCC1", false)
            .unwrap()
            .with_assigned_rings()
            .unwrap()
            .with_assigned_ring_families()
            .unwrap();
        let rings_before = molecule.derived_cache().rings.clone();
        let ring_families_before = molecule.derived_cache().ring_families.clone();

        let result = molecule.with_hydrogens().unwrap();

        assert_eq!(result.derived_cache().rings, rings_before);
        assert_eq!(result.derived_cache().ring_families, ring_families_before);
        assert!(result.derived_cache().valence.is_none());
    }

    #[test]
    fn sanitized_supported_subset_updates_cache_through_operation_pipeline() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let original = molecule.clone();
        let ops = crate::SanitizeOps::CLEANUP
            | crate::SanitizeOps::PROPERTIES
            | crate::SanitizeOps::SYMMRINGS
            | crate::SanitizeOps::FIND_RADICALS
            | crate::SanitizeOps::SET_AROMATICITY
            | crate::SanitizeOps::SET_CONJUGATION;

        let result = molecule.sanitized_with_ops(ops).unwrap();

        assert_eq!(molecule, original);
        let cache = result.derived_cache();
        assert!(cache.valence.is_some());
        assert!(cache.rings.is_some());
        assert!(cache.ring_families.is_none());
        assert!(cache.aromaticity_valid);
    }

    #[test]
    fn sanitized_set_aromaticity_recomputes_valence_after_aromatic_bond_updates() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("C1=CC=CC=C1", false).unwrap();
        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES
                    | crate::SanitizeOps::SYMMRINGS
                    | crate::SanitizeOps::SET_AROMATICITY,
            )
            .unwrap();

        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| bond.order() == crate::BondOrder::Aromatic)
        );
        let expected_valence =
            crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, true)
                .unwrap();
        assert_eq!(result.derived_cache().valence, Some(expected_valence));
    }

    #[test]
    fn sanitized_kekulize_runs_kekulize_assignment_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::SYMMRINGS | crate::SanitizeOps::KEKULIZE)
            .unwrap();

        assert!(result.bonds().iter().all(|bond| !bond.is_aromatic()));
        assert!(result.bonds().iter().all(|bond| matches!(
            bond.order(),
            crate::BondOrder::Single | crate::BondOrder::Double
        )));
        assert_eq!(
            result
                .bonds()
                .iter()
                .filter(|bond| bond.order() == crate::BondOrder::Double)
                .count(),
            3
        );
        assert!(result.derived_cache().rings.is_some());
    }

    #[test]
    fn sanitized_kekulize_materializes_ring_cache_without_explicit_symmrings_step_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::KEKULIZE)
            .unwrap();

        assert!(result.derived_cache().rings.is_some());
        assert!(result.bonds().iter().all(|bond| matches!(
            bond.order(),
            crate::BondOrder::Single | crate::BondOrder::Double
        )));
    }

    #[test]
    fn sanitized_reports_kekulize_failure_step_like_rdkit_operation_that_failed() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c", false).unwrap();

        let err = molecule
            .sanitized_with_ops(crate::SanitizeOps::KEKULIZE)
            .unwrap_err();

        match err {
            OperationError::Sanitize { source, .. } => {
                assert_eq!(source.step, crate::SanitizeStep::Kekulize);
                assert!(source.message.contains("aromatic"));
            }
            other => panic!("expected sanitize error, got {other:?}"),
        }
    }

    #[test]
    fn sanitized_reports_properties_before_later_requested_steps() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();

        let err = molecule
            .sanitized_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::KEKULIZE)
            .unwrap_err();

        match err {
            OperationError::Sanitize { source, .. } => {
                assert_eq!(source.step, crate::SanitizeStep::Properties);
                assert!(source.message.contains("greater than permitted"));
            }
            other => panic!("expected sanitize error, got {other:?}"),
        }
    }

    #[test]
    fn sanitized_reports_properties_before_multiple_later_requested_steps() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();

        let err = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES
                    | crate::SanitizeOps::KEKULIZE
                    | crate::SanitizeOps::SET_AROMATICITY
                    | crate::SanitizeOps::SET_CONJUGATION,
            )
            .unwrap_err();

        match err {
            OperationError::Sanitize { source, .. } => {
                assert_eq!(source.step, crate::SanitizeStep::Properties);
                assert!(source.message.contains("greater than permitted"));
            }
            other => panic!("expected sanitize error, got {other:?}"),
        }
    }

    #[test]
    fn sanitize_stage_maps_requested_step_through_shared_helper() {
        let err = sanitize_stage(
            crate::SanitizeStep::Cleanup,
            || -> Result<(), crate::ValenceError> {
                Err(crate::ValenceError::UnsupportedBranch {
                    reason: "helper step mapping regression",
                })
            },
            |step, source| sanitize_valence_error(&SANITIZED_SPEC, step, source),
        )
        .unwrap_err();

        match err {
            OperationError::Sanitize { source, .. } => {
                assert_eq!(source.step, crate::SanitizeStep::Cleanup);
                assert!(source.message.contains("helper step mapping regression"));
            }
            other => panic!("expected sanitize error, got {other:?}"),
        }
    }

    #[test]
    fn sanitized_without_properties_uses_non_strict_property_cache_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("C(=O)(=O)(=O)", false).unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::NONE)
            .unwrap();
        let expected =
            crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, false)
                .unwrap();

        assert_eq!(result.derived_cache().valence, Some(expected));
    }

    #[test]
    fn sanitized_cleanup_converts_neutral_nitro_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let oxygen_single = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let oxygen_double = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                carbon,
                nitrogen,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                oxygen_single,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                oxygen_double,
                crate::BondOrder::Double,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap();

        assert_eq!(result.atoms()[nitrogen.index()].formal_charge(), 1);
        assert_eq!(
            result
                .atoms()
                .iter()
                .filter(|atom| atom.atomic_number() == 8 && atom.formal_charge() == -1)
                .count(),
            1
        );
        assert_eq!(
            result
                .bonds()
                .iter()
                .filter(|bond| bond.order() == crate::BondOrder::Double)
                .count(),
            1
        );
    }

    #[test]
    fn sanitized_nitrogens_cleanup_rewrites_neutral_nitrogen_triple_bond_branch_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let carbon_one = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let nitrogen_center = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let carbon_two = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let nitrogen_terminal = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_one,
                nitrogen_center,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen_center,
                carbon_two,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let triple_bond = builder
            .add_bond(crate::BondSpec::new(
                nitrogen_center,
                nitrogen_terminal,
                crate::BondOrder::Triple,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap();

        assert_eq!(result.atoms()[nitrogen_center.index()].formal_charge(), 1);
        assert_eq!(
            result.atoms()[nitrogen_terminal.index()].formal_charge(),
            -1
        );
        assert_eq!(
            result.bonds()[triple_bond.index()].order(),
            crate::BondOrder::Double
        );
    }

    #[test]
    fn sanitized_cleanup_converts_phosphorus_oxo_like_rdkit() {
        let phosphorus_element =
            crate::Element::from_atomic_number(15).expect("phosphorus atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let carbon_single = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let phosphorus = builder.add_atom(crate::AtomSpec::new(phosphorus_element));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let carbon_double = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_single,
                phosphorus,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                phosphorus,
                oxygen,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                phosphorus,
                nitrogen,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                carbon_double,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap();

        assert_eq!(result.atoms()[phosphorus.index()].formal_charge(), 1);
        assert_eq!(result.atoms()[oxygen.index()].formal_charge(), -1);
        assert_eq!(result.atoms()[nitrogen.index()].formal_charge(), 0);
        assert_eq!(result.bonds()[1].order(), crate::BondOrder::Single);
        assert_eq!(result.bonds()[2].order(), crate::BondOrder::Double);
    }

    #[test]
    fn sanitized_phosphorus_cleanup_leaves_double_oxo_without_double_cn_branch_unchanged_like_rdkit()
     {
        let phosphorus_element =
            crate::Element::from_atomic_number(15).expect("phosphorus atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let phosphorus = builder.add_atom(crate::AtomSpec::new(phosphorus_element));
        let oxygen_one = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let oxygen_two = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let bond_one = builder
            .add_bond(crate::BondSpec::new(
                phosphorus,
                oxygen_one,
                crate::BondOrder::Double,
            ))
            .unwrap();
        let bond_two = builder
            .add_bond(crate::BondSpec::new(
                phosphorus,
                oxygen_two,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                phosphorus,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap();

        assert_eq!(result.atoms()[phosphorus.index()].formal_charge(), 0);
        assert_eq!(result.atoms()[oxygen_one.index()].formal_charge(), 0);
        assert_eq!(result.atoms()[oxygen_two.index()].formal_charge(), 0);
        assert_eq!(
            result.bonds()[bond_one.index()].order(),
            crate::BondOrder::Double
        );
        assert_eq!(
            result.bonds()[bond_two.index()].order(),
            crate::BondOrder::Double
        );
    }

    #[test]
    fn sanitized_cleanup_converts_hypervalent_halogen_oxo_like_rdkit() {
        let chlorine =
            crate::Element::from_atomic_number(17).expect("chlorine atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(chlorine));
        let oxygen_one = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let oxygen_two = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                center,
                oxygen_one,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                center,
                oxygen_two,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap();

        assert_eq!(result.atoms()[center.index()].formal_charge(), 1);
        assert_eq!(result.atoms()[oxygen_one.index()].formal_charge(), -1);
        assert_eq!(result.atoms()[oxygen_two.index()].formal_charge(), 0);
        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| bond.order() == crate::BondOrder::Single)
        );
    }

    #[test]
    fn sanitized_halogen_cleanup_skips_non_oxo_neighbor_branch_like_rdkit() {
        let chlorine =
            crate::Element::from_atomic_number(17).expect("chlorine atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(chlorine));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let double_bond = builder
            .add_bond(crate::BondSpec::new(
                center,
                oxygen,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                center,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap();

        assert_eq!(result.atoms()[center.index()].formal_charge(), 0);
        assert_eq!(result.atoms()[oxygen.index()].formal_charge(), 0);
        assert_eq!(
            result.bonds()[double_bond.index()].order(),
            crate::BondOrder::Double
        );
    }

    #[test]
    fn cleanup_incident_bonds_returns_only_local_bond_indices() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        let center_left = builder
            .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Double))
            .unwrap();
        let center_right = builder
            .add_bond(crate::BondSpec::new(a1, a3, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let adjacency = sanitize_adjacency(&molecule).unwrap();

        let incident = sanitize_cleanup_incident_bonds(&adjacency, a1);

        assert_eq!(incident, vec![0, center_left.index(), center_right.index()]);
    }

    #[test]
    fn cleanup_incident_bonds_explicit_valence_uses_assignment_bond_orders() {
        let mut builder = crate::MoleculeBuilder::new();
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen_bond = builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                oxygen,
                crate::BondOrder::Double,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        let adjacency = sanitize_adjacency(&molecule).unwrap();
        let mut assignment = SanitizeCleanupAssignment {
            atom_formal_charges: molecule
                .atoms()
                .iter()
                .map(crate::Atom::formal_charge)
                .collect(),
            bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
        };
        assignment.bond_orders[oxygen_bond.index()] = crate::BondOrder::Single;

        let valence =
            sanitize_cleanup_explicit_valence(&molecule, &adjacency, &assignment, nitrogen)
                .unwrap();

        assert_eq!(valence, 2);
    }

    #[test]
    fn sanitized_organometallic_cleanup_converts_single_metal_bond_to_dative_like_rdkit() {
        let iron = crate::Element::from_atomic_number(26).expect("iron atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let carbon_one = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_two = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_three = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal = builder.add_atom(crate::AtomSpec::new(iron));
        for neighbor in [carbon_one, carbon_two, carbon_three, metal] {
            builder
                .add_bond(crate::BondSpec::new(
                    nitrogen,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
            .unwrap();
        let metal_bond = result
            .bonds()
            .iter()
            .find(|bond| {
                (bond.begin() == nitrogen && bond.end() == metal)
                    || (bond.begin() == metal && bond.end() == nitrogen)
            })
            .unwrap();

        assert_eq!(metal_bond.order(), crate::BondOrder::Dative);
        assert_eq!(metal_bond.begin(), nitrogen);
        assert_eq!(metal_bond.end(), metal);
    }

    #[test]
    fn sanitized_organometallic_cleanup_prefers_metal_with_fewer_existing_dative_bonds() {
        let iron = crate::Element::from_atomic_number(26).expect("iron atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();

        let donor = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let donor_c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let donor_c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let donor_c3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal_busy = builder.add_atom(crate::AtomSpec::new(iron));
        let metal_open = builder.add_atom(crate::AtomSpec::new(iron));

        for neighbor in [donor_c1, donor_c2, donor_c3, metal_busy, metal_open] {
            builder
                .add_bond(crate::BondSpec::new(
                    donor,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }

        let donor_busy = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let busy_c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let busy_c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let busy_c3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for neighbor in [busy_c1, busy_c2, busy_c3] {
            builder
                .add_bond(crate::BondSpec::new(
                    donor_busy,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .add_bond(crate::BondSpec::new(
                donor_busy,
                metal_busy,
                crate::BondOrder::Dative,
            ))
            .unwrap();

        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
            .unwrap();

        let donor_to_busy = result
            .bonds()
            .iter()
            .find(|bond| {
                (bond.begin() == donor && bond.end() == metal_busy)
                    || (bond.begin() == metal_busy && bond.end() == donor)
            })
            .unwrap();
        let donor_to_open = result
            .bonds()
            .iter()
            .find(|bond| {
                (bond.begin() == donor && bond.end() == metal_open)
                    || (bond.begin() == metal_open && bond.end() == donor)
            })
            .unwrap();

        assert_eq!(donor_to_busy.order(), crate::BondOrder::Single);
        assert_eq!(donor_to_open.order(), crate::BondOrder::Dative);
        assert_eq!(donor_to_open.begin(), donor);
        assert_eq!(donor_to_open.end(), metal_open);
    }

    #[test]
    fn sanitized_organometallic_cleanup_skips_non_hypervalent_donor_like_rdkit() {
        let iron = crate::Element::from_atomic_number(26).expect("iron atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let oxygen =
            builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_no_implicit(true));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal = builder.add_atom(crate::AtomSpec::new(iron));
        builder
            .add_bond(crate::BondSpec::new(
                oxygen,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let metal_bond = builder
            .add_bond(crate::BondSpec::new(
                oxygen,
                metal,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
            .unwrap();

        assert_eq!(
            result.bonds()[metal_bond.index()].order(),
            crate::BondOrder::Single
        );
    }

    #[test]
    fn metal_bond_cleanup_prefers_higher_rank_when_dative_counts_tie() {
        let iron = crate::Element::from_atomic_number(26).expect("iron atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let donor = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal_plain = builder.add_atom(crate::AtomSpec::new(iron));
        let metal_substituted = builder.add_atom(crate::AtomSpec::new(iron));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        for neighbor in [c1, c2, c3, metal_plain, metal_substituted] {
            builder
                .add_bond(crate::BondSpec::new(
                    donor,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .add_bond(crate::BondSpec::new(
                metal_substituted,
                hydrogen,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        let adjacency = sanitize_adjacency(&molecule).unwrap();
        let valence = molecule
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
            .unwrap();
        let ranks = molecule.rank_mol_atoms().unwrap();
        let mut assignment = SanitizeOrganometallicCleanupAssignment {
            bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
            bond_endpoints: molecule
                .bonds()
                .iter()
                .map(|bond| (bond.begin(), bond.end()))
                .collect(),
        };

        sanitize_metal_bond_cleanup_assignment(
            &molecule,
            &adjacency,
            &valence,
            &ranks,
            donor,
            &mut assignment,
        )
        .unwrap();

        let chosen_metal = assignment
            .bond_endpoints
            .iter()
            .zip(assignment.bond_orders.iter())
            .find_map(|(&(begin, end), &order)| {
                (order == crate::BondOrder::Dative && begin == donor).then_some(end)
            })
            .unwrap();
        let expected = [metal_plain, metal_substituted]
            .into_iter()
            .max_by_key(|atom| ranks[atom.index()])
            .unwrap();
        assert_eq!(chosen_metal, expected);
    }

    #[test]
    fn hypervalent_nonmetal_predicate_matches_metal_and_aromatic_degree_four_branches() {
        let carbon = crate::Element::C;
        let iron = crate::Element::from_atomic_number(26).unwrap();

        let mut aromatic_builder = crate::MoleculeBuilder::new();
        let sulfur_atom = aromatic_builder.add_atom(
            crate::AtomSpec::new(carbon)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        let mut aromatic_neighbors = Vec::new();
        for _ in 0..4 {
            let carbon = aromatic_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            aromatic_builder
                .add_bond(crate::BondSpec::new(
                    sulfur_atom,
                    carbon,
                    crate::BondOrder::Single,
                ))
                .unwrap();
            aromatic_neighbors.push(carbon);
        }
        let aromatic = aromatic_builder.build().unwrap();
        let aromatic_adj = sanitize_adjacency(&aromatic).unwrap();
        let aromatic_valence = aromatic
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
            .unwrap();

        assert!(
            sanitize_is_hypervalent_nonmetal(
                &aromatic,
                &aromatic_adj,
                &aromatic_valence,
                sulfur_atom
            )
            .unwrap()
        );

        let mut metal_builder = crate::MoleculeBuilder::new();
        let metal_atom = metal_builder.add_atom(crate::AtomSpec::new(iron));
        let ligand = metal_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        metal_builder
            .add_bond(crate::BondSpec::new(
                metal_atom,
                ligand,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let metal_molecule = metal_builder.build().unwrap();
        let metal_adj = sanitize_adjacency(&metal_molecule).unwrap();
        let metal_valence = metal_molecule
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
            .unwrap();

        assert!(
            !sanitize_is_hypervalent_nonmetal(
                &metal_molecule,
                &metal_adj,
                &metal_valence,
                metal_atom
            )
            .unwrap()
        );
    }

    #[test]
    fn single_bonded_metals_filters_non_single_and_non_metal_neighbors() {
        let iron = crate::Element::from_atomic_number(26).unwrap();
        let mut builder = crate::MoleculeBuilder::new();
        let donor = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let metal_single = builder.add_atom(crate::AtomSpec::new(iron));
        let metal_rewritten = builder.add_atom(crate::AtomSpec::new(iron));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let keep_bond = builder
            .add_bond(crate::BondSpec::new(
                donor,
                metal_single,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let rewritten_bond = builder
            .add_bond(crate::BondSpec::new(
                donor,
                metal_rewritten,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                donor,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        let adjacency = sanitize_adjacency(&molecule).unwrap();
        let mut assignment = SanitizeOrganometallicCleanupAssignment {
            bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
            bond_endpoints: molecule
                .bonds()
                .iter()
                .map(|bond| (bond.begin(), bond.end()))
                .collect(),
        };
        assignment.bond_orders[rewritten_bond.index()] = crate::BondOrder::Dative;

        let metals =
            sanitize_organometallic_single_bonded_metals(&molecule, &adjacency, &assignment, donor);

        assert_eq!(metals, vec![metal_single]);
        assert_eq!(keep_bond.index(), 0);
    }

    #[test]
    fn sanitized_cleanup_atropisomers_clears_non_sp2_atrop_bond_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let left = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let right = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(left, right, crate::BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCw),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_ATROPISOMERS)
            .unwrap();

        assert_eq!(result.bonds()[0].stereo(), crate::BondStereo::None);
        assert_eq!(result.bonds()[0].stereo_atoms(), None);
    }

    #[test]
    fn sanitized_cleanup_atropisomers_clears_small_ring_atrop_stereo_and_group_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let atoms = (0..6)
            .map(|_| {
                builder.add_atom(
                    crate::AtomSpec::new(crate::Element::C)
                        .with_hybridization(crate::Hybridization::Sp2),
                )
            })
            .collect::<Vec<_>>();
        let atrop_bond = builder
            .add_bond(
                crate::BondSpec::new(atoms[0], atoms[1], crate::BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCcw),
            )
            .unwrap();
        for idx in 1..6 {
            builder
                .add_bond(crate::BondSpec::new(
                    atoms[idx],
                    atoms[(idx + 1) % 6],
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .add_stereo_group(crate::StereoGroup::new(
                crate::StereoGroupKind::Or,
                Vec::new(),
                vec![atrop_bond],
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_ATROPISOMERS)
            .unwrap();

        assert_eq!(
            result.bonds()[atrop_bond.index()].stereo(),
            crate::BondStereo::None
        );
        assert_eq!(result.bonds()[atrop_bond.index()].stereo_atoms(), None);
        assert!(result.stereo_groups().is_empty());
    }

    #[test]
    fn sanitized_cleanup_chirality_clears_non_sp3_tetrahedral_tag_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_chiral_tag(crate::ChiralTag::TetrahedralCw),
        );
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
            .unwrap();

        assert_eq!(
            result.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
    }

    #[test]
    fn sanitized_cleanup_chirality_cleans_stereo_groups_for_non_sp3_tetrahedral_tags_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let atom = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_chiral_tag(crate::ChiralTag::TetrahedralCw),
        );
        builder
            .add_stereo_group(crate::StereoGroup::new(
                crate::StereoGroupKind::Absolute,
                vec![atom],
                Vec::new(),
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
            .unwrap();

        assert_eq!(
            result.atoms()[atom.index()].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert!(result.stereo_groups().is_empty());
    }

    #[test]
    fn sanitized_cleanup_chirality_resets_tetrahedral_permutation_above_limit_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_chiral_tag(crate::ChiralTag::Tetrahedral)
                .with_hybridization(crate::Hybridization::Sp3)
                .with_chiral_permutation(7),
        );
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
            .unwrap();

        assert_eq!(
            result.atoms()[0].chiral_tag(),
            crate::ChiralTag::Tetrahedral
        );
        assert_eq!(result.atoms()[0].chiral_permutation(), Some(0));
    }

    #[test]
    fn sanitized_cleanup_chirality_resets_square_planar_permutation_above_limit_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_no_implicit(true)
                .with_chiral_tag(crate::ChiralTag::SquarePlanar)
                .with_chiral_permutation(7),
        );
        let left = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let right = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(center, left, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                center,
                right,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
            .unwrap();

        assert_eq!(
            result.atoms()[center.index()].chiral_tag(),
            crate::ChiralTag::SquarePlanar
        );
        assert_eq!(result.atoms()[center.index()].chiral_permutation(), Some(0));
    }

    #[test]
    fn sanitized_cleanup_chirality_leaves_invalid_square_planar_stereo_group_untouched_without_tetrahedral_cleanup_flag_like_rdkit()
     {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_no_implicit(true)
                .with_chiral_tag(crate::ChiralTag::SquarePlanar)
                .with_chiral_permutation(3),
        );
        let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(
                center,
                neighbor,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_stereo_group(crate::StereoGroup::new(
                crate::StereoGroupKind::Absolute,
                vec![center],
                Vec::new(),
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
            .unwrap();

        assert_eq!(
            result.atoms()[center.index()].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(result.stereo_groups().len(), 1);
        assert_eq!(result.stereo_groups()[0].atoms(), &[center]);
    }

    #[test]
    fn sanitized_sets_conjugation_for_butadiene_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("C=CC=C", false).unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_CONJUGATION,
            )
            .unwrap();

        assert!(result.bonds().iter().all(crate::Bond::is_conjugated));
    }

    #[test]
    fn sanitized_set_conjugation_keeps_aromatic_bonds_conjugated_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::SET_CONJUGATION)
            .unwrap();

        assert!(result.bonds().iter().all(crate::Bond::is_conjugated));
    }

    #[test]
    fn sanitized_set_conjugation_uses_heteroatom_lone_pair_candidate_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("NC=O", false).unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_CONJUGATION,
            )
            .unwrap();

        assert_eq!(result.num_bonds(), 2);
        assert!(result.bonds().iter().all(crate::Bond::is_conjugated));
    }

    #[test]
    fn sanitized_sets_hybridization_after_conjugation_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES
                    | crate::SanitizeOps::SET_CONJUGATION
                    | crate::SanitizeOps::SET_HYBRIDIZATION,
            )
            .unwrap();

        assert_eq!(result.atoms()[0].hybridization(), crate::Hybridization::Sp3);
        assert_eq!(result.atoms()[1].hybridization(), crate::Hybridization::Sp3);
        assert_eq!(result.atoms()[2].hybridization(), crate::Hybridization::Sp3);
    }

    #[test]
    fn sanitized_set_hybridization_uses_chiral_tag_coordination_override() {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_chiral_tag(crate::ChiralTag::TetrahedralCw),
        );
        for _ in 0..4 {
            let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(
                    center,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION,
            )
            .unwrap();

        assert_eq!(
            result.atoms()[center.index()].hybridization(),
            crate::Hybridization::Sp3
        );
    }

    #[test]
    fn sanitized_set_hybridization_excludes_dative_bonds_from_num_bonds_plus_lone_pairs() {
        let iron = crate::Element::from_atomic_number(26).unwrap();
        let mut builder = crate::MoleculeBuilder::new();
        let nitrogen =
            builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_no_implicit(true));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal = builder.add_atom(crate::AtomSpec::new(iron));
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                nitrogen,
                metal,
                crate::BondOrder::Dative,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION,
            )
            .unwrap();

        assert_eq!(
            result.atoms()[nitrogen.index()].hybridization(),
            crate::Hybridization::Sp2
        );
    }

    #[test]
    fn sanitized_set_hybridization_excludes_zero_bonds_from_num_bonds_plus_lone_pairs() {
        let mut builder = crate::MoleculeBuilder::new();
        let oxygen =
            builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_no_implicit(true));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let dummy = builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        builder
            .add_bond(crate::BondSpec::new(
                oxygen,
                carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(oxygen, dummy, crate::BondOrder::Zero))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION,
            )
            .unwrap();

        assert_eq!(
            result.atoms()[oxygen.index()].hybridization(),
            crate::Hybridization::Sp2
        );
    }

    #[test]
    fn sanitized_set_hybridization_uses_conjugated_bond_sp2_branch_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("NC=O", false).unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES
                    | crate::SanitizeOps::SET_CONJUGATION
                    | crate::SanitizeOps::SET_HYBRIDIZATION,
            )
            .unwrap();

        assert_eq!(result.atoms()[0].hybridization(), crate::Hybridization::Sp2);
    }

    #[test]
    fn sanitized_set_hybridization_uses_atomic_number_cutoff_for_actinides() {
        let actinium = crate::Element::from_atomic_number(89).unwrap();
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(actinium).with_no_implicit(true));
        for _ in 0..2 {
            let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(
                    center,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION,
            )
            .unwrap();

        assert_eq!(
            result.atoms()[center.index()].hybridization(),
            crate::Hybridization::Sp
        );
    }

    #[test]
    fn sanitized_find_radicals_recomputes_property_cache_after_topology_state_update() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(3),
        );
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::FIND_RADICALS)
            .unwrap();

        assert_eq!(result.atoms()[0].radical_electrons(), 1);
        let expected =
            crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, false)
                .unwrap();
        assert_eq!(result.derived_cache().valence, Some(expected));
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_requested_step_does_not_mutate_source() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let original = molecule.clone();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
            .unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.num_atoms(), molecule.num_atoms());
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_materializes_disappearing_pyrrolic_hydrogen() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("N1C=CC=C1", false).unwrap();

        let result = molecule
            .sanitized_with_ops(
                crate::SanitizeOps::PROPERTIES
                    | crate::SanitizeOps::SYMMRINGS
                    | crate::SanitizeOps::KEKULIZE
                    | crate::SanitizeOps::SET_AROMATICITY
                    | crate::SanitizeOps::ADJUST_HYDROGENS,
            )
            .unwrap();

        assert!(result.atoms()[0].is_aromatic());
        assert_eq!(result.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(
            result
                .derived_cache()
                .valence
                .as_ref()
                .unwrap()
                .implicit_hydrogens[0],
            0
        );
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_preserves_existing_explicit_hydrogen_when_delta_is_zero()
     {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1),
        );
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
            .unwrap();

        assert_eq!(result.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_leaves_stable_explicit_hydrogens_unchanged_like_rdkit()
    {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

        let result = molecule
            .sanitized_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
            .unwrap();

        let explicit_hs = result
            .atoms()
            .iter()
            .map(crate::Atom::explicit_hydrogens)
            .collect::<Vec<_>>();
        assert_eq!(explicit_hs, vec![0, 0, 0]);
        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(crate::Atom::explicit_hydrogens)
                .collect::<Vec<_>>(),
            explicit_hs
        );
    }

    #[test]
    fn experimental_kekulize_runs_through_operation_pipeline_without_changing_source() {
        assert_eq!(
            WITH_KEKULIZED_BONDS_SPEC.support,
            SupportStatus::Experimental
        );

        let molecule = crate::Molecule::new();
        let original = molecule.clone();
        let result = molecule
            .with_kekulized_bonds(true)
            .expect("experimental kekulize skeleton should satisfy op contract");

        assert_eq!(molecule, original);
        assert_eq!(result.atoms(), original.atoms());
        assert_eq!(result.bonds(), original.bonds());
        assert_eq!(result.coords_2d(), original.coords_2d());
        assert_eq!(result.conformers_3d(), original.conformers_3d());
        assert_eq!(
            result.source_coordinate_dim(),
            original.source_coordinate_dim()
        );
        assert_eq!(result.properties(), original.properties());
        assert_eq!(
            result.derived_cache().valence,
            Some(crate::ValenceAssignment {
                explicit_valence: Vec::new(),
                implicit_hydrogens: Vec::new(),
            })
        );
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn op_parts_rejects_permission_violation_under_strict_checks() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);
        let err = parts
            .begin_coordinates_mut()
            .expect_err("coordinate begin should be rejected");
        assert!(
            matches!(err, OperationError::InvalidInput { message, .. } if message.contains("outside its registry access"))
        );
    }

    #[test]
    fn needs_update_clears_matching_derived_cache_entries() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &SANITIZED_SPEC);
        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
        parts.mark_aromaticity_valid();

        parts.clear_cache(SANITIZED_SPEC.needs_update());
        let result = parts
            .finish(OpOutcome::Changed)
            .expect("cache invalidation should satisfy operation contract");

        let cache = result.derived_cache();
        assert!(cache.valence.is_none());
        assert!(cache.rings.is_none());
        assert!(cache.ring_families.is_none());
        assert!(!cache.aromaticity_valid);
        assert!(!cache.stereo_valid);
    }

    #[test]
    fn needs_update_accepts_cache_updates_without_prior_clear() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &SANITIZED_SPEC);

        parts.set_rings_cache(crate::RingInfo::new(crate::RingFindType::SymmSssr, 0, 0));
        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
        parts.mark_aromaticity_valid();
        parts.clear_cache(
            DerivedState::RING_FAMILIES
                | DerivedState::STEREO
                | DerivedState::DRAWING
                | DerivedState::FINGERPRINT,
        );
        let result = parts
            .finish(OpOutcome::Changed)
            .expect("updated cache entries should satisfy needs_update without clear first");

        let cache = result.derived_cache();
        assert!(cache.valence.is_some());
        assert!(cache.rings.is_some());
        assert!(cache.ring_families.is_none());
        assert!(cache.aromaticity_valid);
        assert!(!cache.stereo_valid);
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn finish_rejects_missing_needs_update_handling() {
        let molecule = crate::Molecule::new();
        let parts = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC);

        let err = parts
            .finish(OpOutcome::NoOp {
                reason: "intentionally missed needs_update",
            })
            .expect_err("needs_update must be cleared or updated before finish");

        assert!(matches!(
            err,
            OperationError::InvalidInput {
                message: "operation body did not clear or update every required cache state",
                ..
            }
        ));
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn finish_rejects_unrelated_cache_clear_for_needs_update() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC);

        parts.clear_cache(DerivedState::RINGS);
        let err = parts
            .finish(OpOutcome::Changed)
            .expect_err("clearing rings must not satisfy a valence needs_update contract");

        assert!(matches!(
            err,
            OperationError::InvalidInput {
                message: "operation body did not clear or update every required cache state",
                ..
            }
        ));
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn finish_rejects_declared_preservation_without_proof() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC);
        let topology = parts.begin_topology_mut().unwrap();
        let coordinates = parts.begin_coordinates_mut().unwrap();
        let properties = parts.begin_properties_mut().unwrap();

        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();
        parts
            .record_topology_edit(TopologyEditKind::Appending)
            .unwrap();
        parts.record_topology_mapping(TopologyMapping::with_appended(0, 0, 0, 0));
        parts.clear_cache(WITH_HYDROGENS_SPEC.needs_update());

        let err = parts
            .finish(OpOutcome::Changed)
            .expect_err("declared preserve states require an explicit preservation proof");

        assert!(matches!(
            err,
            OperationError::InvalidInput {
                message: "operation body did not prove every declared preserved derived state",
                ..
            }
        ));
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn leaf_append_preservation_proof_rejects_non_leaf_appended_atom() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().unwrap();
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC);
        let mut topology = parts.begin_topology_mut().unwrap();
        let coordinates = parts.begin_coordinates_mut().unwrap();
        let properties = parts.begin_properties_mut().unwrap();

        let appended = AtomId::new(topology.atoms.len());
        topology.atoms.push(crate::Atom::from_spec(
            appended,
            crate::AtomSpec::new(crate::Element::H),
        ));
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();
        parts
            .record_topology_edit(TopologyEditKind::Appending)
            .unwrap();
        parts.record_topology_mapping(TopologyMapping::with_appended(1, 0, 1, 0));

        let err = parts
            .prove_preserved(
                DerivedState::RINGS | DerivedState::RING_FAMILIES,
                PreservationProof::LeafAtomAppend,
            )
            .expect_err("an appended atom with no appended leaf bond must not preserve rings");

        assert!(matches!(
            err,
            OperationError::InvalidInput {
                message: "leaf-append preservation proof requires every appended atom to be a degree-one leaf",
                ..
            }
        ));
    }

    #[test]
    fn finish_accepts_update_path_for_recompute() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_RECOMPUTE_VALENCE_SPEC);

        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
        let result = parts
            .finish(OpOutcome::Changed)
            .expect("setting valence cache should satisfy recompute requirement");

        assert_eq!(
            result.derived_cache().valence,
            Some(crate::ValenceAssignment {
                explicit_valence: Vec::new(),
                implicit_hydrogens: Vec::new(),
            })
        );
        assert_eq!(molecule.derived_cache().valence, None);
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "cache clear permission violation")]
    fn clear_cache_panics_without_derived_cache_permission_when_cache_is_touched() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC);
        parts.clear_cache(DerivedState::VALENCE);
    }

    #[cfg(not(feature = "op-contracts"))]
    #[test]
    fn op_contract_checks_are_disabled_without_feature() {
        let molecule = crate::Molecule::new();

        let mut unauthorized = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC);
        unauthorized.clear_cache(DerivedState::VALENCE);
        unauthorized
            .finish(OpOutcome::Changed)
            .expect("without op-contracts, cache permission checks are disabled");

        let missing_update = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC);
        missing_update
            .finish(OpOutcome::Changed)
            .expect("without op-contracts, needs_update checks are disabled");
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "cache write permission violation")]
    fn set_valence_cache_panics_without_requires_or_recompute() {
        // WITH_2D_COORDINATES_SPEC has no requires/recompute — set_valence_cache must panic
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC);
        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "cache clear permission violation")]
    fn clear_cache_panics_without_invalidate_or_recompute() {
        // TEST_RECOMPUTE_VALENCE_SPEC has VALENCE in recompute (not invalidate)
        // RINGS is in neither invalidate nor recompute — clear_cache(DerivedState::RINGS) must panic
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_RECOMPUTE_VALENCE_SPEC);
        parts.clear_cache(DerivedState::RINGS);
    }

    #[test]
    fn op_parts_cow_mutation_changes_result_without_changing_source() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);

        let mut topology = parts.begin_topology_mut().unwrap();
        topology.atoms.push(crate::Atom::from_spec(
            crate::AtomId::new(0),
            crate::AtomSpec::new(crate::Element::C),
        ));
        let view = parts.read_parts_for_topology(topology.clone()).unwrap();
        let valence = MoleculeReadParts::from_molecule(&view)
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
            .unwrap();
        parts.commit_topology(topology).unwrap();
        parts.record_topology_edit(TopologyEditKind::Local).unwrap();
        parts.set_rings_cache(crate::RingInfo::new(crate::RingFindType::SymmSssr, 1, 0));
        parts.set_valence_cache(valence);
        parts.clear_cache(DerivedState::AROMATICITY);
        parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);

        let result = parts
            .finish(OpOutcome::Changed)
            .expect("COW topology edit should produce a valid molecule");

        assert_eq!(molecule.num_atoms(), 0);
        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.atomic_numbers(), vec![6]);
    }

    #[test]
    fn compacting_edit_uses_begin_commit_blocks_and_records_mapping() {
        let mut builder = crate::Molecule::builder();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o1 = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let n2 = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        builder
            .add_bond(crate::BondSpec::new(c0, o1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(o1, n2, crate::BondOrder::Double))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
            .unwrap();
        let properties = crate::MoleculeProperties::default()
            .with_sdf_property_list(crate::SdfPropertyList::new(
                crate::SdfPropertyListTarget::Atom,
                "atom_tag",
                vec![
                    Some("c0".to_string()),
                    Some("o1".to_string()),
                    Some("n2".to_string()),
                ],
            ))
            .with_sdf_property_list(crate::SdfPropertyList::new(
                crate::SdfPropertyListTarget::Bond,
                "bond_tag",
                vec![Some("c-o".to_string()), Some("o-n".to_string())],
            ));
        builder = builder.with_properties(properties);
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let mut parts = OpParts::new(&molecule, &WITHOUT_HYDROGENS_SPEC);
        let mut topology = parts.begin_topology_mut().unwrap();
        let mut coordinates = parts.begin_coordinates_mut().unwrap();
        let mut properties = parts.begin_properties_mut().unwrap();
        let mapping = topology.remove_atoms_with_mapping(&[o1]);
        coordinates.remap_topology(&mapping);
        properties.remap_topology(&mapping);
        parts
            .record_topology_edit(TopologyEditKind::Compacting)
            .unwrap();
        parts.record_topology_mapping(mapping.clone());
        assert_eq!(
            mapping.atoms().old_to_new(),
            &[
                Some(crate::AtomId::new(0)),
                None,
                Some(crate::AtomId::new(1))
            ]
        );
        assert_eq!(mapping.bonds().old_to_new(), &[None, None]);
        parts.clear_cache(WITHOUT_HYDROGENS_SPEC.needs_update());
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();

        let result = parts
            .finish(OpOutcome::Changed)
            .expect("strong compacting edit should satisfy registry contract");

        assert_eq!(molecule, original);
        assert_eq!(result.atomic_numbers(), vec![6, 7]);
        assert_eq!(result.num_bonds(), 0);
        assert_eq!(result.coords_2d().unwrap(), &[[0.0, 0.0], [2.0, 0.0]]);
        assert_eq!(
            result.properties().sdf_property_lists()[0].values(),
            &[Some("c0".to_string()), Some("n2".to_string())]
        );
        assert_eq!(result.properties().sdf_property_lists()[1].values(), &[]);
    }

    #[test]
    fn strong_remove_atoms_remaps_surviving_sgroup_parent_relationships() {
        let mut builder = crate::Molecule::builder();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o1 = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let n2 = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::Superatom,
                )
                .with_atoms(vec![c0]),
            )
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(1),
                    crate::SubstanceGroupKind::Data,
                )
                .with_atoms(vec![o1])
                .with_parent(crate::SubstanceGroupId::new(0)),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let mut parts = OpParts::new(&molecule, &WITHOUT_HYDROGENS_SPEC);
        let mut topology = parts.begin_topology_mut().unwrap();
        let mut coordinates = parts.begin_coordinates_mut().unwrap();
        let mut properties = parts.begin_properties_mut().unwrap();
        let mapping = topology.remove_atoms_with_mapping(&[n2]);
        coordinates.remap_topology(&mapping);
        properties.remap_topology(&mapping);
        parts
            .record_topology_edit(TopologyEditKind::Compacting)
            .unwrap();
        parts.record_topology_mapping(mapping);
        parts.clear_cache(WITHOUT_HYDROGENS_SPEC.needs_update());
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();

        let result = parts
            .finish(OpOutcome::Changed)
            .expect("strong compacting edit should preserve surviving SGroup parent links");

        assert_eq!(result.substance_groups().len(), 2);
        assert_eq!(
            result.substance_groups()[0].atoms(),
            &[crate::AtomId::new(0)]
        );
        assert_eq!(
            result.substance_groups()[1].atoms(),
            &[crate::AtomId::new(1)]
        );
        assert_eq!(
            result.substance_groups()[1].parent(),
            Some(crate::SubstanceGroupId::new(0))
        );
    }

    #[test]
    fn with_hydrogens_extends_sdf_property_lists_for_appended_atoms_and_bonds() {
        let mut builder = crate::Molecule::builder();
        let carbon = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_explicit_hydrogens(1)
                .with_no_implicit(true),
        );
        let properties = crate::MoleculeProperties::default()
            .with_sdf_property_list(crate::SdfPropertyList::new(
                crate::SdfPropertyListTarget::Atom,
                "atom_tag",
                vec![Some("c0".to_string())],
            ))
            .with_sdf_property_list(crate::SdfPropertyList::new(
                crate::SdfPropertyListTarget::Bond,
                "bond_tag",
                Vec::new(),
            ));
        builder = builder.with_properties(properties);
        let molecule = builder.build().unwrap();

        let result = molecule.with_hydrogens().unwrap();

        assert_eq!(result.num_atoms(), 2);
        assert_eq!(
            result.properties().sdf_property_lists()[0].values(),
            &[Some("c0".to_string()), None]
        );
        assert_eq!(
            result.properties().sdf_property_lists()[1].values(),
            &[None]
        );
        assert_eq!(result.atoms()[carbon.index()].explicit_hydrogens(), 0);
    }

    #[test]
    fn compacting_topology_edit_record_is_rejected_for_weak_operations() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);
        let err = parts
            .record_topology_edit(TopologyEditKind::Compacting)
            .expect_err("weak operations must not record compacting topology edits");
        assert!(matches!(err, OperationError::InvalidInput { .. }));
    }

    #[test]
    fn with_kekulized_bonds_uses_rdkit_wrapper_canonical_default_for_benzene() {
        let molecule = crate::Molecule::from_smiles("C1=CC=CC=C1").unwrap();

        let kekulized = molecule.with_kekulized_bonds(false).unwrap();
        let bond_orders = kekulized
            .bonds()
            .iter()
            .map(|bond| bond.order())
            .collect::<Vec<_>>();

        assert_eq!(
            bond_orders,
            vec![
                crate::BondOrder::Single,
                crate::BondOrder::Double,
                crate::BondOrder::Single,
                crate::BondOrder::Double,
                crate::BondOrder::Single,
                crate::BondOrder::Double
            ]
        );
        assert!(kekulized.bonds().iter().all(|bond| bond.is_aromatic()));
        assert!(kekulized.atoms().iter().all(|atom| atom.is_aromatic()));
    }
}
