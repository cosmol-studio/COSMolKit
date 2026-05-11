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

use std::{fmt, marker::PhantomData};

use cosmolkit_macros::{mol_op_body, molecule_ops};

use crate::{
    AtomId, DerivedState, InvariantError, Molecule, MoleculeProperties, SupportStatus,
    invariants::enforce_molecule_invariants,
    molecule::{CoordinateBlock, DerivedCacheBlock, TopologyBlock, TopologyMapping},
};

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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OperationReportSet(u8);

impl OperationReportSet {
    pub const NONE: Self = Self(0);
    pub const ATOM_MAPPING: Self = Self(1 << 0);
    pub const BOND_MAPPING: Self = Self(1 << 1);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedStateSet(u32);

impl DerivedStateSet {
    pub const NONE: Self = Self(0);
    pub const ADJACENCY: Self = Self(1 << 0);
    pub const RINGS: Self = Self(1 << 1);
    pub const RING_FAMILIES: Self = Self(1 << 2);
    pub const VALENCE: Self = Self(1 << 3);
    pub const AROMATICITY: Self = Self(1 << 4);
    pub const STEREO: Self = Self(1 << 5);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }

    #[must_use]
    pub const fn as_derived_state(self) -> DerivedState {
        let mut states = DerivedState::NONE;
        if self.contains(Self::ADJACENCY) {
            states = states.union(DerivedState::ADJACENCY);
        }
        if self.contains(Self::RINGS) {
            states = states.union(DerivedState::RINGS);
        }
        if self.contains(Self::RING_FAMILIES) {
            states = states.union(DerivedState::RING_FAMILIES);
        }
        if self.contains(Self::VALENCE) {
            states = states.union(DerivedState::VALENCE);
        }
        if self.contains(Self::AROMATICITY) {
            states = states.union(DerivedState::AROMATICITY);
        }
        if self.contains(Self::STEREO) {
            states = states.union(DerivedState::STEREO);
        }
        states
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
    pub may_mutate: BlockSet,
    pub auto_remap: BlockSet,
    pub must_handle: DerivedStateSet,
    pub needs_update: DerivedState,
    pub requires_mapping: MappingRequirement,
    pub report: OperationReportSet,
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
    #[error("{operation}: kekulize error: {source}")]
    Kekulize {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::KekulizeError,
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
    remapped_blocks: BlockSet,
    handled: DerivedStateSet,
    cleared_cache: DerivedState,
    updated_cache: DerivedStateSet,
    reported: OperationReportSet,
    outcome: Option<OpOutcome>,
}

impl OperationTrace {
    #[must_use]
    pub const fn touched_blocks(&self) -> BlockSet {
        self.touched_blocks
    }

    #[must_use]
    pub const fn handled(&self) -> DerivedStateSet {
        self.handled
    }

    #[must_use]
    pub const fn remapped_blocks(&self) -> BlockSet {
        self.remapped_blocks
    }

    #[must_use]
    pub const fn cleared_cache(&self) -> DerivedState {
        self.cleared_cache
    }

    #[must_use]
    pub const fn updated_cache(&self) -> DerivedStateSet {
        self.updated_cache
    }

    #[must_use]
    pub fn outcome(&self) -> Option<&OpOutcome> {
        self.outcome.as_ref()
    }
}

pub struct OpParts<'a> {
    spec: &'static MoleculeOpSpec,
    working: Molecule,
    topology_mapping: Option<TopologyMapping>,
    _source: PhantomData<&'a Molecule>,

    #[cfg(feature = "op-contracts")]
    trace: OperationTrace,
}

impl<'a> OpParts<'a> {
    pub(crate) fn new(source: &'a Molecule, spec: &'static MoleculeOpSpec) -> Self {
        Self {
            spec,
            working: source.clone(),
            topology_mapping: None,
            _source: PhantomData,
            #[cfg(feature = "op-contracts")]
            trace: OperationTrace {
                touched_blocks: BlockSet::NONE,
                remapped_blocks: BlockSet::NONE,
                handled: DerivedStateSet::NONE,
                cleared_cache: DerivedState::NONE,
                updated_cache: DerivedStateSet::NONE,
                reported: OperationReportSet::NONE,
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
    #[must_use]
    pub(crate) fn molecule(&self) -> &Molecule {
        &self.working
    }

    #[allow(dead_code)]
    pub(crate) fn topology(&self) -> &TopologyBlock {
        self.working.topology_block()
    }

    pub(crate) fn topology_mut(&mut self) -> &mut TopologyBlock {
        self.record_mutation(BlockSet::TOPOLOGY);
        self.working.topology_block_mut()
    }

    // These accessors are part of the operation body API. They may be unused
    // until the corresponding operation family is ported.
    #[allow(dead_code)]
    pub(crate) fn coordinates_mut(&mut self) -> &mut CoordinateBlock {
        self.record_mutation(BlockSet::COORDINATES);
        self.working.coordinate_block_mut()
    }

    #[allow(dead_code)]
    pub(crate) fn properties_mut(&mut self) -> &mut MoleculeProperties {
        self.record_mutation(BlockSet::PROPERTIES);
        self.working.properties_mut()
    }

    #[allow(dead_code)]
    pub(crate) fn derived_cache_mut(&mut self) -> &mut DerivedCacheBlock {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut()
    }

    #[allow(dead_code)]
    pub(crate) fn remove_atoms(
        &mut self,
        atoms_to_remove: &[AtomId],
    ) -> Result<&TopologyMapping, OperationError> {
        self.assert_compacting_topology_edit_allowed()?;

        let mapping = {
            self.record_mutation(BlockSet::TOPOLOGY);
            self.working
                .topology_block_mut()
                .remove_atoms_with_mapping(atoms_to_remove)
        };

        self.remap_blocks_required_by_spec(&mapping)?;
        self.clear_cache(self.spec.needs_update);
        self.set_topology_mapping(mapping);

        Ok(self
            .topology_mapping
            .as_ref()
            .expect("topology mapping was just recorded"))
    }

    pub(crate) fn set_adjacency_cache(&mut self, adjacency: crate::AdjacencyList) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().adjacency = Some(adjacency);
        self.record_updated_cache(DerivedStateSet::ADJACENCY);
    }

    pub(crate) fn set_rings_cache(&mut self, rings: crate::RingInfo) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().rings = Some(rings);
        self.record_updated_cache(DerivedStateSet::RINGS);
    }

    pub(crate) fn set_ring_families_cache(&mut self, ring_families: crate::RingInfo) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().ring_families = Some(ring_families);
        self.record_updated_cache(DerivedStateSet::RING_FAMILIES);
    }

    pub(crate) fn set_valence_cache(&mut self, valence: crate::ValenceAssignment) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().valence = Some(valence);
        self.record_updated_cache(DerivedStateSet::VALENCE);
    }

    pub(crate) fn mark_aromaticity_valid(&mut self) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().aromaticity_valid = true;
        self.record_updated_cache(DerivedStateSet::AROMATICITY);
    }

    #[allow(dead_code)]
    pub(crate) fn mark_stereo_handled(&mut self) {
        self.record_mutation(BlockSet::DERIVED_CACHE);
        self.working.derived_cache_mut().stereo_valid = true;
        self.record_updated_cache(DerivedStateSet::STEREO);
    }

    pub(crate) fn clear_cache(&mut self, states: DerivedState) {
        #[cfg(feature = "op-contracts")]
        {
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

    fn remap_blocks_required_by_spec(
        &mut self,
        mapping: &TopologyMapping,
    ) -> Result<(), OperationError> {
        if self.spec.auto_remap.contains(BlockSet::COORDINATES) {
            self.record_mutation(BlockSet::COORDINATES);
            self.working.coordinate_block_mut().remap_topology(mapping);
            self.record_remapped(BlockSet::COORDINATES);
        }

        if self.spec.auto_remap.contains(BlockSet::PROPERTIES) {
            self.record_mutation(BlockSet::PROPERTIES);
            self.record_remapped(BlockSet::PROPERTIES);
        }

        Ok(())
    }

    fn set_topology_mapping(&mut self, mapping: TopologyMapping) {
        self.topology_mapping = Some(mapping);
        self.record_reported(
            OperationReportSet::ATOM_MAPPING.union(OperationReportSet::BOND_MAPPING),
        );
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
                self.spec.may_mutate.contains(block),
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

    fn record_handled(&mut self, state: DerivedStateSet) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.handled = self.trace.handled.union(state);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = state;
        }
    }

    fn record_updated_cache(&mut self, state: DerivedStateSet) {
        self.record_handled(state);
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
    fn validate_contract(&self) -> Result<(), OperationError> {
        if !self.trace.handled.contains(self.spec.must_handle) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not handle every required derived state",
            });
        }

        let updated_or_cleared =
            self.trace.cleared_cache | self.trace.updated_cache.as_derived_state();
        if !updated_or_cleared.contains(self.spec.needs_update) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not clear or update every required cache state",
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

        if !self.trace.reported.contains(self.spec.report) {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation did not record every registry-required report",
            });
        }

        Ok(())
    }

    fn assert_compacting_topology_edit_allowed(&self) -> Result<(), OperationError> {
        if self.spec.kind != MoleculeOpKind::Strong {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "compacting topology edits require a strong operation",
            });
        }
        if self.spec.topology_edit != TopologyEditKind::Compacting {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "operation registry does not allow compacting topology edits",
            });
        }
        if self.spec.requires_mapping != MappingRequirement::Required {
            return Err(OperationError::InvalidInput {
                operation: self.spec,
                message: "compacting topology edits must require a topology mapping",
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

    fn record_reported(&mut self, report: OperationReportSet) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.reported = self.trace.reported.union(report);
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = report;
        }
    }
}

molecule_ops! {
    op with_hydrogens() {
        method: with_hydrogens,
        impl_fn: with_hydrogens_impl,
        domain: topology,
        kind: strong,
        topology_edit: compacting,
        may_mutate: [topology, coordinates, derived_cache],
        auto_remap: [coordinates],
        must_handle: [valence, aromaticity, stereo],
        needs_update: [adjacency, rings, ring_families, valence, aromaticity, stereo, drawing, fingerprint],
        requires_mapping: required,
        report: [atom_mapping, bond_mapping],
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
        may_mutate: [topology, coordinates, derived_cache],
        auto_remap: [coordinates],
        must_handle: [valence, aromaticity, stereo],
        needs_update: [adjacency, rings, ring_families, valence, aromaticity, stereo, drawing, fingerprint],
        requires_mapping: required,
        report: [atom_mapping, bond_mapping],
        allows_noop: true,
        feature: HYDROGENS_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "strong_topology_with_coordinates",
        parity_profile: "remove_hs_default_rdkit",
    }

    op with_kekulized_bonds(clear_aromatic_flags: bool) {
        method: with_kekulized_bonds,
        impl_fn: with_kekulized_bonds_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        must_handle: [valence],
        needs_update: [aromaticity, valence, drawing, fingerprint],
        requires_mapping: none,
        report: [],
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
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        must_handle: [],
        needs_update: [rings, ring_families, valence, aromaticity, stereo, drawing, fingerprint],
        requires_mapping: none,
        report: [],
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
        may_mutate: [derived_cache],
        auto_remap: [],
        must_handle: [adjacency, valence],
        needs_update: [adjacency, valence],
        requires_mapping: none,
        report: [],
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
        may_mutate: [derived_cache],
        auto_remap: [],
        must_handle: [adjacency, rings],
        needs_update: [adjacency, rings],
        requires_mapping: none,
        report: [],
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
        may_mutate: [derived_cache],
        auto_remap: [],
        must_handle: [adjacency, ring_families],
        needs_update: [adjacency, ring_families],
        requires_mapping: none,
        report: [],
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
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        must_handle: [adjacency, rings, valence, aromaticity],
        needs_update: [adjacency, rings, valence, aromaticity, drawing, fingerprint],
        requires_mapping: none,
        report: [],
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
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        must_handle: [adjacency, valence],
        needs_update: [adjacency, valence],
        requires_mapping: none,
        report: [],
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
        may_mutate: [coordinates],
        auto_remap: [],
        must_handle: [],
        needs_update: [drawing],
        requires_mapping: none,
        report: [],
        allows_noop: true,
        feature: COORDINATE_2D_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_generation",
        parity_profile: "compute_2d_coords_default_rdkit",
    }
}

#[mol_op_body(with_hydrogens, parts)]
fn with_hydrogens_impl() -> Result<OpOutcome, OperationError> {
    let _ = parts;
    Err(OperationError::UnsupportedFeature {
        operation: &WITH_HYDROGENS_SPEC,
        source: crate::UnsupportedFeatureError::from_spec(&crate::HYDROGENS_FEATURE),
    })
}

#[mol_op_body(without_hydrogens, parts)]
fn without_hydrogens_impl(sanitize: bool) -> Result<OpOutcome, OperationError> {
    let _ = (parts, sanitize);
    Err(OperationError::UnsupportedFeature {
        operation: &WITHOUT_HYDROGENS_SPEC,
        source: crate::UnsupportedFeatureError::from_spec(&crate::HYDROGENS_FEATURE),
    })
}

#[mol_op_body(with_kekulized_bonds, parts)]
fn with_kekulized_bonds_impl(clear_aromatic_flags: bool) -> Result<OpOutcome, OperationError> {
    if parts.molecule().derived_cache().rings.is_none() {
        let rings = crate::symmetrize_sssr(parts.molecule()).map_err(|source| {
            OperationError::RingFinding {
                operation: &WITH_KEKULIZED_BONDS_SPEC,
                source,
            }
        })?;
        parts.set_rings_cache(rings);
    }
    let ring_info = parts
        .molecule()
        .derived_cache()
        .rings
        .as_ref()
        .expect("rings were recomputed immediately above");
    let assignment = crate::kekulize::kekulize_assignment(
        parts.molecule(),
        Some(ring_info),
        clear_aromatic_flags,
        false,
        100,
    )
    .map_err(|source| OperationError::Kekulize {
        operation: &WITH_KEKULIZED_BONDS_SPEC,
        source,
    })?;

    let changed = crate::kekulize::apply_kekulize_assignment(parts.topology_mut(), &assignment);
    parts.clear_cache(DerivedState::AROMATICITY);
    let valence =
        crate::assign_valence_with_options(parts.molecule(), crate::ValenceModel::RdkitLike, true)
            .map_err(|source| OperationError::Valence {
                operation: &WITH_KEKULIZED_BONDS_SPEC,
                source,
            })?;
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
    // BEGIN RDKIT CPP FUNCTION MolOps::sanitizeMol
    // RDKit❗✔️: void sanitizeMol(RWMol &mol) {
    // RDKit❗✔️:   unsigned int failedOp = 0;
    // RDKit❗✔️:   sanitizeMol(mol, failedOp, SANITIZE_ALL);
    // RDKit❗✔️: }
    // RDKit❗✔️: void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
    // RDKit❗✔️:                  unsigned int sanitizeOps) {
    // RDKit❗✔️:   // clear out any cached properties
    // RDKit❗✔️:   mol.clearComputedProps();
    parts.clear_cache(SANITIZED_SPEC.needs_update);

    let mut changed = false;

    // RDKit❌❌:   operationThatFailed = SANITIZE_CLEANUP;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     // clean up things like nitro groups
    // RDKit❌❌:     cleanUp(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::CLEANUP)
        && sanitize_cleanup_may_need_real_work(parts.molecule())
    {
        return Err(unsupported_sanitize_step("SANITIZE_CLEANUP is not ported"));
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     cleanUpOrganometallics(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
        && sanitize_organometallic_cleanup_may_need_real_work(parts.molecule())
    {
        return Err(unsupported_sanitize_step(
            "SANITIZE_CLEANUP_ORGANOMETALLICS is not ported",
        ));
    }

    // RDKit❗✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit❗✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit❗✔️:     mol.updatePropertyCache(true);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     mol.updatePropertyCache(false);
    // RDKit❗✔️:   }
    let adjacency = crate::AdjacencyList::try_from_topology(
        parts.molecule().num_atoms(),
        parts.molecule().bonds(),
    )
    .map_err(|source| OperationError::InvalidInput {
        operation: &SANITIZED_SPEC,
        message: match source {
            crate::AdjacencyError::BondAtomOutOfRange { .. } => {
                "topology bond atom index out of range while recomputing adjacency"
            }
        },
    })?;
    parts.set_adjacency_cache(adjacency);
    let valence = crate::assign_valence_with_options(
        parts.molecule(),
        crate::ValenceModel::RdkitLike,
        ops.contains(crate::SanitizeOps::PROPERTIES),
    )
    .map_err(|source| OperationError::Valence {
        operation: &SANITIZED_SPEC,
        source,
    })?;
    parts.set_valence_cache(valence);

    // RDKit❗✔️:   operationThatFailed = SANITIZE_SYMMRINGS;
    // RDKit❗✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit❗✔️:     VECT_INT_VECT arings;
    // RDKit❗✔️:     MolOps::symmetrizeSSSR(mol, arings);
    // RDKit❗✔️:   }
    if ops.contains(crate::SanitizeOps::SYMMRINGS) {
        let rings = crate::symmetrize_sssr(parts.molecule()).map_err(|source| {
            OperationError::RingFinding {
                operation: &SANITIZED_SPEC,
                source,
            }
        })?;
        parts.set_rings_cache(rings);
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_KEKULIZE;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     Kekulize(mol, true, false);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::KEKULIZE) {
        return Err(unsupported_sanitize_step("SANITIZE_KEKULIZE is not ported"));
    }

    // RDKit❗✔️:   operationThatFailed = SANITIZE_FINDRADICALS;
    // RDKit❗✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit❗✔️:     assignRadicals(mol);
    // RDKit❗✔️:   }
    if ops.contains(crate::SanitizeOps::FIND_RADICALS) {
        let radicals =
            crate::assign_radicals(parts.molecule()).map_err(|source| OperationError::Valence {
                operation: &SANITIZED_SPEC,
                source,
            })?;
        let radicals_changed = parts
            .molecule()
            .atoms()
            .iter()
            .zip(radicals.iter().copied())
            .any(|(atom, radical)| atom.radical_electrons() != radical);
        if radicals_changed {
            let topology = parts.topology_mut();
            for (atom, radical) in topology.atoms.iter_mut().zip(radicals) {
                atom.set_radical_electrons(radical);
            }
            let valence = crate::assign_valence_with_options(
                parts.molecule(),
                crate::ValenceModel::RdkitLike,
                true,
            )
            .map_err(|source| OperationError::Valence {
                operation: &SANITIZED_SPEC,
                source,
            })?;
            parts.set_valence_cache(valence);
        }
        changed |= radicals_changed;
    }

    // RDKit❗✔️:   operationThatFailed = SANITIZE_SETAROMATICITY;
    // RDKit❗✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit❗✔️:     setAromaticity(mol);
    // RDKit❗✔️:   }
    if ops.contains(crate::SanitizeOps::SET_AROMATICITY) {
        if parts.molecule().derived_cache().rings.is_none() {
            let rings = crate::symmetrize_sssr(parts.molecule()).map_err(|source| {
                OperationError::RingFinding {
                    operation: &SANITIZED_SPEC,
                    source,
                }
            })?;
            parts.set_rings_cache(rings);
        }
        let assignment = crate::set_aromaticity(parts.molecule(), crate::AromaticityModel::Default)
            .map_err(|source| OperationError::Aromaticity {
                operation: &SANITIZED_SPEC,
                source,
            })?;
        {
            let topology = parts.topology_mut();
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
        }
        let valence = crate::assign_valence_with_options(
            parts.molecule(),
            crate::ValenceModel::RdkitLike,
            ops.contains(crate::SanitizeOps::PROPERTIES),
        )
        .map_err(|source| OperationError::Valence {
            operation: &SANITIZED_SPEC,
            source,
        })?;
        parts.set_valence_cache(valence);
        parts.mark_aromaticity_valid();
        parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
        changed = true;
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_SETCONJUGATION;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     setConjugation(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::SET_CONJUGATION)
        && sanitize_conjugation_may_need_real_work(parts.molecule())
    {
        return Err(unsupported_sanitize_step(
            "SANITIZE_SETCONJUGATION is not ported",
        ));
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_SETHYBRIDIZATION;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     setHybridization(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::SET_HYBRIDIZATION) {
        return Err(unsupported_sanitize_step(
            "SANITIZE_SETHYBRIDIZATION is not ported",
        ));
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_CLEANUPATROPISOMERS;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     cleanupAtropisomers(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_ATROPISOMERS) {
        return Err(unsupported_sanitize_step(
            "SANITIZE_CLEANUPATROPISOMERS is not ported",
        ));
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     cleanupChirality(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::CLEANUP_CHIRALITY) {
        return Err(unsupported_sanitize_step(
            "SANITIZE_CLEANUPCHIRALITY is not ported",
        ));
    }

    // RDKit❌❌:   operationThatFailed = SANITIZE_ADJUSTHS;
    // RDKit❌❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit❌❌:     adjustHs(mol);
    // RDKit❌❌:   }
    if ops.contains(crate::SanitizeOps::ADJUST_HYDROGENS) {
        return Err(unsupported_sanitize_step("SANITIZE_ADJUSTHS is not ported"));
    }

    // RDKit❗✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit❗✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit❗✔️:     mol.updatePropertyCache(true);
    // RDKit❗✔️:   }
    if ops.contains(crate::SanitizeOps::PROPERTIES) {
        let valence = crate::assign_valence_with_options(
            parts.molecule(),
            crate::ValenceModel::RdkitLike,
            true,
        )
        .map_err(|source| OperationError::Valence {
            operation: &SANITIZED_SPEC,
            source,
        })?;
        parts.set_valence_cache(valence);
    }
    // RDKit❗✔️:   operationThatFailed = 0;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::sanitizeMol
    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "requested sanitize steps only updated derived state",
        }
    })
}

fn unsupported_sanitize_step(reason: &'static str) -> OperationError {
    OperationError::Unsupported {
        operation: &SANITIZED_SPEC,
        reason,
    }
}

fn sanitize_cleanup_may_need_real_work(molecule: &Molecule) -> bool {
    molecule
        .atoms()
        .iter()
        .any(|atom| atom.formal_charge() != 0 || atom.radical_electrons() != 0)
}

fn sanitize_organometallic_cleanup_may_need_real_work(molecule: &Molecule) -> bool {
    molecule
        .atoms()
        .iter()
        .any(|atom| is_metal_atomic_number(atom.atomic_number()))
}

fn sanitize_conjugation_may_need_real_work(molecule: &Molecule) -> bool {
    molecule.bonds().iter().any(|bond| {
        matches!(
            bond.order(),
            crate::BondOrder::Double
                | crate::BondOrder::Triple
                | crate::BondOrder::Aromatic
                | crate::BondOrder::OneAndHalf
                | crate::BondOrder::TwoAndHalf
                | crate::BondOrder::ThreeAndHalf
                | crate::BondOrder::FourAndHalf
                | crate::BondOrder::FiveAndHalf
        ) || bond.is_aromatic()
    })
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
    let adjacency = crate::AdjacencyList::try_from_topology(
        parts.molecule().num_atoms(),
        parts.molecule().bonds(),
    )
    .map_err(|source| OperationError::InvalidInput {
        operation: &ASSIGNED_VALENCE_SPEC,
        message: match source {
            crate::AdjacencyError::BondAtomOutOfRange { .. } => {
                "topology bond atom index out of range while recomputing adjacency"
            }
        },
    })?;
    parts.set_adjacency_cache(adjacency);
    let valence =
        crate::assign_valence_with_options(parts.molecule(), crate::ValenceModel::RdkitLike, true)
            .map_err(|source| OperationError::Valence {
                operation: &ASSIGNED_VALENCE_SPEC,
                source,
            })?;
    parts.set_valence_cache(valence);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_rings, parts)]
fn assigned_rings_impl() -> Result<OpOutcome, OperationError> {
    let adjacency = crate::AdjacencyList::try_from_topology(
        parts.molecule().num_atoms(),
        parts.molecule().bonds(),
    )
    .map_err(|source| OperationError::InvalidInput {
        operation: &ASSIGNED_RINGS_SPEC,
        message: match source {
            crate::AdjacencyError::BondAtomOutOfRange { .. } => {
                "topology bond atom index out of range while recomputing adjacency"
            }
        },
    })?;
    parts.set_adjacency_cache(adjacency);
    let rings =
        crate::symmetrize_sssr(parts.molecule()).map_err(|source| OperationError::RingFinding {
            operation: &ASSIGNED_RINGS_SPEC,
            source,
        })?;
    parts.set_rings_cache(rings);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_ring_families, parts)]
fn assigned_ring_families_impl() -> Result<OpOutcome, OperationError> {
    let adjacency = crate::AdjacencyList::try_from_topology(
        parts.molecule().num_atoms(),
        parts.molecule().bonds(),
    )
    .map_err(|source| OperationError::InvalidInput {
        operation: &ASSIGNED_RING_FAMILIES_SPEC,
        message: match source {
            crate::AdjacencyError::BondAtomOutOfRange { .. } => {
                "topology bond atom index out of range while recomputing adjacency"
            }
        },
    })?;
    parts.set_adjacency_cache(adjacency);
    let ring_families =
        crate::find_ring_families(parts.molecule(), false, false).map_err(|source| {
            OperationError::RingFinding {
                operation: &ASSIGNED_RING_FAMILIES_SPEC,
                source,
            }
        })?;
    parts.set_ring_families_cache(ring_families);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_aromaticity, parts)]
fn assigned_aromaticity_impl() -> Result<OpOutcome, OperationError> {
    let adjacency = crate::AdjacencyList::try_from_topology(
        parts.molecule().num_atoms(),
        parts.molecule().bonds(),
    )
    .map_err(|source| OperationError::InvalidInput {
        operation: &ASSIGNED_AROMATICITY_SPEC,
        message: match source {
            crate::AdjacencyError::BondAtomOutOfRange { .. } => {
                "topology bond atom index out of range while recomputing adjacency"
            }
        },
    })?;
    parts.set_adjacency_cache(adjacency);
    let rings =
        crate::symmetrize_sssr(parts.molecule()).map_err(|source| OperationError::RingFinding {
            operation: &ASSIGNED_AROMATICITY_SPEC,
            source,
        })?;
    parts.set_rings_cache(rings);
    let assignment = crate::set_aromaticity(parts.molecule(), crate::AromaticityModel::Default)
        .map_err(|source| OperationError::Aromaticity {
            operation: &ASSIGNED_AROMATICITY_SPEC,
            source,
        })?;
    {
        let topology = parts.topology_mut();
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
    }
    let valence =
        crate::assign_valence_with_options(parts.molecule(), crate::ValenceModel::RdkitLike, true)
            .map_err(|source| OperationError::Valence {
                operation: &ASSIGNED_AROMATICITY_SPEC,
                source,
            })?;
    parts.set_valence_cache(valence);
    parts.mark_aromaticity_valid();
    parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(assigned_radicals, parts)]
fn assigned_radicals_impl() -> Result<OpOutcome, OperationError> {
    let adjacency = crate::AdjacencyList::try_from_topology(
        parts.molecule().num_atoms(),
        parts.molecule().bonds(),
    )
    .map_err(|source| OperationError::InvalidInput {
        operation: &ASSIGNED_RADICALS_SPEC,
        message: match source {
            crate::AdjacencyError::BondAtomOutOfRange { .. } => {
                "topology bond atom index out of range while recomputing adjacency"
            }
        },
    })?;
    parts.set_adjacency_cache(adjacency);
    let radicals =
        crate::assign_radicals(parts.molecule()).map_err(|source| OperationError::Valence {
            operation: &ASSIGNED_RADICALS_SPEC,
            source,
        })?;

    let changed = parts
        .molecule()
        .atoms()
        .iter()
        .zip(radicals.iter().copied())
        .any(|(atom, radical)| atom.radical_electrons() != radical);

    if changed {
        let topology = parts.topology_mut();
        for (atom, radical) in topology.atoms.iter_mut().zip(radicals) {
            atom.set_radical_electrons(radical);
        }
    }

    let valence =
        crate::assign_valence_with_options(parts.molecule(), crate::ValenceModel::RdkitLike, true)
            .map_err(|source| OperationError::Valence {
                operation: &ASSIGNED_RADICALS_SPEC,
                source,
            })?;
    parts.set_valence_cache(valence);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(with_2d_coordinates, parts)]
fn with_2d_coordinates_impl() -> Result<OpOutcome, OperationError> {
    let _ = parts;
    Err(OperationError::UnsupportedFeature {
        operation: &WITH_2D_COORDINATES_SPEC,
        source: crate::UnsupportedFeatureError::from_spec(&crate::COORDINATE_2D_FEATURE),
    })
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    const TEST_NEEDS_VALENCE_UPDATE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_needs_valence_update",
        impl_fn: "test_needs_valence_update_impl",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        may_mutate: BlockSet::DERIVED_CACHE,
        auto_remap: BlockSet::NONE,
        must_handle: DerivedStateSet::NONE,
        needs_update: DerivedState::VALENCE,
        requires_mapping: MappingRequirement::None,
        report: OperationReportSet::NONE,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_MUST_HANDLE_VALENCE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_must_handle_valence",
        impl_fn: "test_must_handle_valence_impl",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::None,
        may_mutate: BlockSet::DERIVED_CACHE,
        auto_remap: BlockSet::NONE,
        must_handle: DerivedStateSet::VALENCE,
        needs_update: DerivedState::VALENCE,
        requires_mapping: MappingRequirement::None,
        report: OperationReportSet::NONE,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

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
        assert_eq!(MOLECULE_OPS.len(), 10);
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
    fn unsupported_operations_return_structured_errors_without_changing_source() {
        let molecule = crate::Molecule::new();
        let original = molecule.clone();

        assert_unsupported_feature(
            molecule.with_hydrogens(),
            &WITH_HYDROGENS_SPEC,
            &crate::HYDROGENS_FEATURE,
        );
        assert_unsupported_feature(
            molecule.without_hydrogens(),
            &WITHOUT_HYDROGENS_SPEC,
            &crate::HYDROGENS_FEATURE,
        );
        assert_unsupported_feature(
            molecule.without_hydrogens_with_sanitize(false),
            &WITHOUT_HYDROGENS_SPEC,
            &crate::HYDROGENS_FEATURE,
        );
        assert!(matches!(
            molecule.sanitized(),
            Err(OperationError::Unsupported {
                operation,
                reason
            }) if same_operation(operation, &SANITIZED_SPEC)
                && reason == "SANITIZE_KEKULIZE is not ported"
        ));
        assert!(matches!(
            molecule.sanitized_with_ops(crate::SanitizeOps::ALL),
            Err(OperationError::Unsupported {
                operation,
                reason
            }) if same_operation(operation, &SANITIZED_SPEC)
                && reason == "SANITIZE_KEKULIZE is not ported"
        ));
        assert_unsupported_feature(
            molecule.with_2d_coordinates(),
            &WITH_2D_COORDINATES_SPEC,
            &crate::COORDINATE_2D_FEATURE,
        );

        assert_eq!(molecule, original);
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
        assert!(cache.adjacency.is_some());
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
    fn sanitized_unported_requested_step_fails_without_committing_cache_changes() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let original = molecule.clone();

        let err = molecule
            .sanitized_with_ops(crate::SanitizeOps::SET_HYBRIDIZATION)
            .unwrap_err();

        assert!(matches!(
            err,
            OperationError::Unsupported {
                operation,
                reason
            } if same_operation(operation, &SANITIZED_SPEC)
                && reason == "SANITIZE_SETHYBRIDIZATION is not ported"
        ));
        assert_eq!(molecule, original);
        assert!(molecule.derived_cache().valence.is_none());
        assert!(molecule.derived_cache().rings.is_none());
        assert!(molecule.derived_cache().ring_families.is_none());
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
    #[should_panic(expected = "attempted to mutate a block outside its registry permissions")]
    fn op_parts_panics_on_permission_violation_under_strict_checks() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);
        let _ = parts.coordinates_mut();
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "attempted to mutate a block outside its registry permissions")]
    fn op_parts_panics_when_derived_cache_mutation_is_not_registered() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC);
        let _ = parts.derived_cache_mut();
    }

    #[test]
    fn needs_update_clears_matching_derived_cache_entries() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &SANITIZED_SPEC);
        {
            let cache = parts.derived_cache_mut();
            cache.adjacency = Some(crate::AdjacencyList::from_topology(0, &[]));
            cache.valence = Some(crate::ValenceAssignment {
                explicit_valence: Vec::new(),
                implicit_hydrogens: Vec::new(),
            });
            cache.rings = Some(crate::RingInfo::new(crate::RingFindType::SymmSssr, 0, 0));
            cache.ring_families = Some(crate::RingInfo::new(crate::RingFindType::SymmSssr, 0, 0));
            cache.aromaticity_valid = true;
            cache.stereo_valid = true;
        }

        parts.clear_cache(SANITIZED_SPEC.needs_update | DerivedState::ADJACENCY);
        let result = parts
            .finish(OpOutcome::Changed)
            .expect("cache invalidation should satisfy operation contract");

        let cache = result.derived_cache();
        assert!(cache.adjacency.is_none());
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
        {
            let cache = parts.derived_cache_mut();
            cache.valence = Some(crate::ValenceAssignment {
                explicit_valence: Vec::new(),
                implicit_hydrogens: Vec::new(),
            });
            cache.rings = Some(crate::RingInfo::new(crate::RingFindType::SymmSssr, 0, 0));
            cache.ring_families = Some(crate::RingInfo::new(crate::RingFindType::SymmSssr, 0, 0));
            cache.aromaticity_valid = true;
            cache.stereo_valid = true;
        }

        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
        parts.mark_aromaticity_valid();
        parts.mark_stereo_handled();
        parts.clear_cache(
            DerivedState::RINGS
                | DerivedState::RING_FAMILIES
                | DerivedState::DRAWING
                | DerivedState::FINGERPRINT,
        );
        let result = parts
            .finish(OpOutcome::Changed)
            .expect("updated cache entries should satisfy needs_update without clear first");

        let cache = result.derived_cache();
        assert!(cache.valence.is_some());
        assert!(cache.rings.is_none());
        assert!(cache.ring_families.is_none());
        assert!(cache.aromaticity_valid);
        assert!(cache.stereo_valid);
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
    fn finish_rejects_cleared_cache_when_must_handle_requires_update() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_MUST_HANDLE_VALENCE_SPEC);

        parts.clear_cache(DerivedState::VALENCE);
        let err = parts
            .finish(OpOutcome::Changed)
            .expect_err("must_handle valence requires an update path, not just clear");

        assert!(matches!(
            err,
            OperationError::InvalidInput {
                message: "operation body did not handle every required derived state",
                ..
            }
        ));
    }

    #[test]
    fn finish_accepts_update_path_for_needs_update_and_must_handle() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_MUST_HANDLE_VALENCE_SPEC);

        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
        let result = parts
            .finish(OpOutcome::Changed)
            .expect("setting valence cache should satisfy both needs_update and must_handle");

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
    #[should_panic(expected = "attempted to mutate a block outside its registry permissions")]
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

    #[test]
    fn op_parts_cow_mutation_changes_result_without_changing_source() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);

        {
            let topology = parts.topology_mut();
            topology.atoms.push(crate::Atom::from_spec(
                crate::AtomId::new(0),
                crate::AtomSpec::new(crate::Element::C),
            ));
        }
        let valence = crate::assign_valence_with_options(
            parts.molecule(),
            crate::ValenceModel::RdkitLike,
            true,
        )
        .unwrap();
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
    fn strong_remove_atoms_atomically_remaps_coordinates_and_records_mapping() {
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
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let mut parts = OpParts::new(&molecule, &WITHOUT_HYDROGENS_SPEC);
        assert_eq!(parts.topology().atoms.len(), 3);
        let mapping = parts.remove_atoms(&[o1]).unwrap();
        assert_eq!(
            mapping.atoms().old_to_new(),
            &[
                Some(crate::AtomId::new(0)),
                None,
                Some(crate::AtomId::new(1))
            ]
        );
        assert_eq!(mapping.bonds().old_to_new(), &[None, None]);
        let valence = crate::assign_valence_with_options(
            parts.molecule(),
            crate::ValenceModel::RdkitLike,
            true,
        )
        .unwrap();
        parts.set_valence_cache(valence);
        parts.mark_aromaticity_valid();
        parts.mark_stereo_handled();

        let result = parts
            .finish(OpOutcome::Changed)
            .expect("strong compacting edit should satisfy registry contract");

        assert_eq!(molecule, original);
        assert_eq!(result.atomic_numbers(), vec![6, 7]);
        assert_eq!(result.num_bonds(), 0);
        assert_eq!(result.coords_2d().unwrap(), &[[0.0, 0.0], [2.0, 0.0]]);
    }

    #[test]
    fn compacting_topology_edit_is_rejected_for_weak_operations() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC);
        let err = parts
            .remove_atoms(&[crate::AtomId::new(0)])
            .expect_err("weak operations must not compact topology");
        assert!(matches!(err, OperationError::InvalidInput { .. }));
    }
}
