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
    TopologyTrust,
    invariants::enforce_molecule_invariants,
    molecule::{CoordinateBlock, DerivedCacheBlock, TopologyBlock, TopologyMapping},
};

pub(crate) use crate::read_parts::MoleculeReadParts;

pub(crate) trait MoleculeReadAccess<'a>: Copy {
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
pub enum SemanticPrecondition {
    TrustedBondTopology,
    HydrogenOwnershipRepresented,
    ValenceComputable,
}

impl fmt::Display for SemanticPrecondition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TrustedBondTopology => f.write_str("trusted_bond_topology"),
            Self::HydrogenOwnershipRepresented => f.write_str("hydrogen_ownership_represented"),
            Self::ValenceComputable => f.write_str("valence_computable"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SemanticPreconditionSet(u8);

impl SemanticPreconditionSet {
    pub const NONE: Self = Self(0);
    pub const TRUSTED_BOND_TOPOLOGY: Self = Self(1 << 0);
    pub const HYDROGEN_OWNERSHIP_REPRESENTED: Self = Self(1 << 1);
    pub const VALENCE_COMPUTABLE: Self = Self(1 << 2);

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
    pub semantic_preconditions: SemanticPreconditionSet,
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
    #[error("{operation}: semantic precondition `{requirement}` failed: {message}")]
    Precondition {
        operation: &'static MoleculeOpSpec,
        requirement: SemanticPrecondition,
        message: &'static str,
    },
    #[error("{operation}: chemistry error: {message}")]
    Chemistry {
        operation: &'static MoleculeOpSpec,
        message: &'static str,
    },
    #[error("{operation}: distance-geometry error: {source}")]
    DistanceGeometry {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::DgBoundsError,
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
                read_cache: DerivedState::NONE,
                cleared_cache: DerivedState::NONE,
                updated_cache: DerivedState::NONE,
                outcome: None,
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
                read_cache: DerivedState::NONE,
                cleared_cache: DerivedState::NONE,
                updated_cache: DerivedState::NONE,
                outcome: None,
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
            self.working.take_topology_block_or_clone()
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
        self.working.replace_topology_block(topology);
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

    pub(crate) fn finish(self, outcome: OpOutcome) -> Result<Molecule, OperationError> {
        debug_assert!(self.in_place_target.is_none());
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

    pub(crate) fn abort_in_place(mut self) {
        let Some(target) = self.in_place_target.take() else {
            return;
        };
        *target = self.working;
    }

    pub(crate) fn finish_in_place(self, outcome: OpOutcome) -> Result<(), OperationError> {
        #[cfg(feature = "op-contracts")]
        {
            let mut this = self;
            this.trace.outcome = Some(outcome);
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
            let _ = outcome;
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

molecule_ops! {
    op with_hydrogens(params: crate::hydrogens::AddHsParams) {
        method: with_hydrogens_with_params,
        impl_fn: with_hydrogens_impl,
        default_method: with_hydrogens,
        default_args: [crate::hydrogens::AddHsParams::default()],
        inplace: true,
        inplace_method: add_hydrogens_with_params_,
        default_inplace_method: add_hydrogens_,
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
        semantic_preconditions: [
            trusted_bond_topology,
            hydrogen_ownership_represented,
            valence_computable,
        ],
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
        inplace: true,
        inplace_method: remove_hydrogens_with_sanitize_,
        default_inplace_method: remove_hydrogens_,
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
        inplace: true,
        inplace_method: remove_hydrogens_with_params_,
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
        inplace: true,
        inplace_method: kekulize_,
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

    op sanitize(ops: crate::SanitizeOps) {
        method: sanitize_with_ops,
        impl_fn: sanitize_impl,
        default_method: sanitize,
        default_args: [crate::SanitizeOps::ALL],
        inplace: true,
        inplace_method: sanitize_with_ops_,
        default_inplace_method: sanitize_,
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
        inplace: true,
        inplace_method: assign_valence_,
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
        inplace: true,
        inplace_method: assign_rings_,
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
        inplace: true,
        inplace_method: assign_ring_families_,
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
        inplace: true,
        inplace_method: assign_aromaticity_,
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
        inplace: true,
        inplace_method: assign_radicals_,
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

    op with_2d_coordinates(params: crate::With2DCoordinatesParams) {
        method: with_2d_coordinates_with_params,
        impl_fn: with_2d_coordinates_impl,
        default_method: with_2d_coordinates,
        default_args: [crate::With2DCoordinatesParams::default()],
        inplace: true,
        inplace_method: compute_2d_coordinates_with_params_,
        default_inplace_method: compute_2d_coordinates_,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        semantic_preconditions: [trusted_bond_topology],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_generation",
        parity_profile: "compute_2d_coords_default_rdkit",
    }

    op with_3d_conformer(params: crate::EmbedParameters) {
        method: with_3d_conformer_with_params,
        impl_fn: with_3d_conformer_impl,
        default_method: with_3d_conformer,
        default_args: [crate::EmbedParameters::etkdg_v3()],
        inplace: true,
        inplace_method: embed_3d_conformer_with_params_,
        default_inplace_method: embed_3d_conformer_,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        semantic_preconditions: [trusted_bond_topology],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: CONFORMER_GENERATION_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_generation_3d",
        parity_profile: "embed_molecule_etkdgv3_rdkit",
    }

    op with_3d_conformers(num_confs: usize, params: crate::EmbedParameters) {
        method: with_3d_conformers_with_params,
        impl_fn: with_3d_conformers_impl,
        inplace: true,
        inplace_method: embed_3d_conformers_with_params_,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        semantic_preconditions: [trusted_bond_topology],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: CONFORMER_GENERATION_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_generation_3d_multi",
        parity_profile: "embed_multiple_confs_rdkit",
    }

    op with_2d_coordinate_block(coords: Vec<[f64; 2]>) {
        method: with_2d_coordinate_block,
        impl_fn: with_2d_coordinate_block_impl,
        inplace: true,
        inplace_method: set_2d_coordinate_block_,
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
        allows_noop: false,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_set_2d",
        parity_profile: "set_2d_coordinates_manual",
    }

    op with_3d_coordinates(coords: Vec<[f64; 3]>, conformer_index: usize) {
        method: with_3d_coordinates,
        impl_fn: with_3d_coordinates_impl,
        inplace: true,
        inplace_method: set_3d_coordinates_,
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
        allows_noop: false,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_set_3d",
        parity_profile: "set_3d_coordinates_manual",
    }

    op with_cleared_3d_conformers() {
        method: with_cleared_3d_conformers,
        impl_fn: with_cleared_3d_conformers_impl,
        inplace: true,
        inplace_method: clear_3d_conformers_,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        requires_mapping: none,
        allows_noop: true,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_clear_3d_conformers",
        parity_profile: "remove_all_conformers_manual",
    }

    op with_added_3d_conformer(coords: Vec<[f64; 3]>, is_3d: bool) {
        method: with_added_3d_conformer,
        impl_fn: with_added_3d_conformer_impl,
        inplace: true,
        inplace_method: add_3d_conformer_,
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
        allows_noop: false,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_add_3d_conformer",
        parity_profile: "add_3d_conformer_manual",
    }

    op with_only_3d_conformer(coords: Vec<[f64; 3]>, is_3d: bool) {
        method: with_only_3d_conformer,
        impl_fn: with_only_3d_conformer_impl,
        inplace: true,
        inplace_method: set_only_3d_conformer_,
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
        allows_noop: false,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_set_single_3d_conformer",
        parity_profile: "remove_all_then_add_conformer_manual",
    }
}

pub(crate) mod hydrogens;
mod sanitize_pipeline;

use self::{hydrogens::*, sanitize_pipeline::*};

pub(crate) use self::sanitize_pipeline::{
    sanitize_conjugation_assignment, sanitize_hybridization_assignment,
};

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
fn with_2d_coordinates_impl(
    params: crate::With2DCoordinatesParams,
) -> Result<OpOutcome, OperationError> {
    let (atoms, bonds) = {
        let read = parts.begin_topology_read()?;
        (read.atoms(), read.bonds())
    };
    let coords = crate::coordinates::compute_2d_coords_with_params(
        atoms,
        bonds,
        &params.as_compute_params(),
    )
    .map_err(|source| match source {
        crate::coordinates::Coordinate2DError::InvalidInput(message) => {
            OperationError::InvalidInput {
                operation: &WITH_2D_COORDINATES_SPEC,
                message,
            }
        }
        crate::coordinates::Coordinate2DError::UnsupportedFeature(_) => {
            OperationError::UnsupportedFeature {
                operation: &WITH_2D_COORDINATES_SPEC,
                source: crate::UnsupportedFeatureError::from_spec(&crate::COORDINATE_2D_FEATURE),
            }
        }
    })?;
    let mut coord_block = parts.begin_coordinates_mut()?;
    if params.clear_confs {
        coord_block.conformers_2d.clear();
    }
    coord_block.conformers_2d.push(crate::Conformer2D::new(
        coord_block.conformers_2d.len(),
        coords,
    ));
    coord_block.source_coordinate_dim = Some(crate::CoordinateDimension::TwoD);
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(with_2d_coordinate_block, parts)]
fn with_2d_coordinate_block_impl(coords: Vec<[f64; 2]>) -> Result<OpOutcome, OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_2D_COORDINATE_BLOCK_SPEC,
            message: "2D coordinate row count mismatch",
        });
    }

    let mut coord_block = parts.begin_coordinates_mut()?;
    coord_block.conformers_2d.clear();
    coord_block
        .conformers_2d
        .push(crate::Conformer2D::new(0, coords));
    coord_block.source_coordinate_dim = Some(crate::CoordinateDimension::TwoD);
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(OpOutcome::Changed)
}

fn source_coordinate_dim_for_block(
    coord_block: &CoordinateBlock,
) -> Option<crate::CoordinateDimension> {
    if coord_block
        .conformers_3d
        .iter()
        .any(crate::Conformer3D::is_3d)
    {
        Some(crate::CoordinateDimension::ThreeD)
    } else if !coord_block.conformers_2d.is_empty() || !coord_block.conformers_3d.is_empty() {
        Some(crate::CoordinateDimension::TwoD)
    } else {
        None
    }
}

#[mol_op_body(with_3d_conformer, parts)]
fn with_3d_conformer_impl(mut params: crate::EmbedParameters) -> Result<OpOutcome, OperationError> {
    let source = parts.working.clone();
    let (embedded, id) =
        crate::distgeom::embed_molecule(&source, &mut params).map_err(|source| {
            OperationError::DistanceGeometry {
                operation: &WITH_3D_CONFORMER_SPEC,
                source,
            }
        })?;
    let _coord_block = parts.begin_coordinates_mut()?;
    let coord_block = embedded.coordinate_block().clone();
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    if id < 0 {
        Ok(OpOutcome::NoOp {
            reason: "RDKit EmbedMolecule returned -1",
        })
    } else {
        Ok(OpOutcome::Changed)
    }
}

#[mol_op_body(with_3d_coordinates, parts)]
fn with_3d_coordinates_impl(
    coords: Vec<[f64; 3]>,
    conformer_index: usize,
) -> Result<OpOutcome, OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_3D_COORDINATES_SPEC,
            message: "3D conformer row count mismatch",
        });
    }

    let mut coord_block = parts.begin_coordinates_mut()?;
    if conformer_index >= coord_block.conformers_3d.len() {
        return Err(OperationError::InvalidInput {
            operation: &WITH_3D_COORDINATES_SPEC,
            message: "3D conformer index out of range",
        });
    }
    let existing = &coord_block.conformers_3d[conformer_index];
    coord_block.conformers_3d[conformer_index] =
        crate::Conformer3D::new(existing.id(), coords, existing.is_3d());
    coord_block.source_coordinate_dim = source_coordinate_dim_for_block(&coord_block);
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(with_cleared_3d_conformers, parts)]
fn with_cleared_3d_conformers_impl() -> Result<OpOutcome, OperationError> {
    let mut coord_block = parts.begin_coordinates_mut()?;
    let changed = !coord_block.conformers_3d.is_empty()
        || matches!(
            coord_block.source_coordinate_dim,
            Some(crate::CoordinateDimension::ThreeD)
        );
    coord_block.conformers_3d.clear();
    coord_block.source_coordinate_dim = source_coordinate_dim_for_block(&coord_block);
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(if changed {
        OpOutcome::Changed
    } else {
        OpOutcome::NoOp {
            reason: "molecule has no 3D conformers to clear",
        }
    })
}

#[mol_op_body(with_3d_conformers, parts)]
fn with_3d_conformers_impl(
    num_confs: usize,
    mut params: crate::EmbedParameters,
) -> Result<OpOutcome, OperationError> {
    let source = parts.working.clone();
    let num_confs = u32::try_from(num_confs).map_err(|_| OperationError::InvalidInput {
        operation: &WITH_3D_CONFORMERS_SPEC,
        message: "num_confs does not fit in RDKit unsigned int parameter",
    })?;
    let (embedded, ids) = crate::distgeom::embed_multiple_confs(&source, num_confs, &mut params)
        .map_err(|source| OperationError::DistanceGeometry {
            operation: &WITH_3D_CONFORMERS_SPEC,
            source,
        })?;
    let _coord_block = parts.begin_coordinates_mut()?;
    let coord_block = embedded.coordinate_block().clone();
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    if ids.is_empty() {
        Ok(OpOutcome::NoOp {
            reason: "RDKit EmbedMultipleConfs returned no conformer IDs",
        })
    } else {
        Ok(OpOutcome::Changed)
    }
}

#[mol_op_body(with_added_3d_conformer, parts)]
fn with_added_3d_conformer_impl(
    coords: Vec<[f64; 3]>,
    is_3d: bool,
) -> Result<OpOutcome, OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_ADDED_3D_CONFORMER_SPEC,
            message: "3D conformer row count mismatch",
        });
    }

    let mut coord_block = parts.begin_coordinates_mut()?;
    let next_id = coord_block
        .conformers_3d
        .iter()
        .map(crate::Conformer3D::id)
        .max()
        .map_or(0, |max_id| max_id + 1);
    coord_block
        .conformers_3d
        .push(crate::Conformer3D::new(next_id, coords, is_3d));
    coord_block.source_coordinate_dim = source_coordinate_dim_for_block(&coord_block);
    parts.commit_coordinates(coord_block)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(OpOutcome::Changed)
}

#[mol_op_body(with_only_3d_conformer, parts)]
fn with_only_3d_conformer_impl(
    coords: Vec<[f64; 3]>,
    is_3d: bool,
) -> Result<OpOutcome, OperationError> {
    let atom_count = {
        let read = parts.begin_topology_read()?;
        read.num_atoms()
    };
    if coords.len() != atom_count {
        return Err(OperationError::InvalidInput {
            operation: &WITH_ONLY_3D_CONFORMER_SPEC,
            message: "3D conformer row count mismatch",
        });
    }

    let mut coord_block = parts.begin_coordinates_mut()?;
    coord_block.conformers_3d.clear();
    coord_block
        .conformers_3d
        .push(crate::Conformer3D::new(0, coords, is_3d));
    coord_block.source_coordinate_dim = source_coordinate_dim_for_block(&coord_block);
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
        semantic_preconditions: SemanticPreconditionSet::NONE,
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
        semantic_preconditions: SemanticPreconditionSet::NONE,
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
        semantic_preconditions: SemanticPreconditionSet::NONE,
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

    #[cfg(feature = "op-contracts")]
    #[test]
    fn begin_topology_mut_rejects_second_begin() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();
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

    #[cfg(feature = "op-contracts")]
    #[test]
    fn begin_topology_mut_rejects_second_begin_before_commit() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();
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

    #[cfg(feature = "op-contracts")]
    #[test]
    fn begin_access_rejects_overlapping_read_and_write_blocks() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &TEST_OVERLAPPING_ACCESS_SPEC).unwrap();
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

    fn support_feature_for(
        operation: &'static MoleculeOpSpec,
    ) -> Option<&'static crate::FeatureSpec> {
        SUPPORT_MATRIX.iter().find_map(|entry| {
            entry
                .operation
                .is_some_and(|candidate| same_operation(candidate, operation))
                .then_some(entry.feature)
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
    fn conformer_generation_operation_registry_exposes_3d_generation_entries() {
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_3D_CONFORMER_SPEC))
        );
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_3D_CONFORMERS_SPEC))
        );
        assert_eq!(
            support_feature_for(&WITH_3D_CONFORMER_SPEC).map(|feature| feature.name),
            Some(crate::CONFORMER_GENERATION_FEATURE.name)
        );
        assert_eq!(
            support_feature_for(&WITH_3D_CONFORMERS_SPEC).map(|feature| feature.name),
            Some(crate::CONFORMER_GENERATION_FEATURE.name)
        );
    }

    #[test]
    fn coordinate_edit_operations_are_registered_as_operation_contracts() {
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_2D_COORDINATE_BLOCK_SPEC))
        );
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_3D_COORDINATES_SPEC))
        );
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_CLEARED_3D_CONFORMERS_SPEC))
        );
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_ADDED_3D_CONFORMER_SPEC))
        );
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_ONLY_3D_CONFORMER_SPEC))
        );
        assert_eq!(
            support_feature_for(&WITH_2D_COORDINATE_BLOCK_SPEC).map(|feature| feature.name),
            Some(crate::COORDINATE_EDIT_FEATURE.name)
        );
        assert_eq!(
            support_feature_for(&WITH_3D_COORDINATES_SPEC).map(|feature| feature.name),
            Some(crate::COORDINATE_EDIT_FEATURE.name)
        );
        assert_eq!(
            support_feature_for(&WITH_CLEARED_3D_CONFORMERS_SPEC).map(|feature| feature.name),
            Some(crate::COORDINATE_EDIT_FEATURE.name)
        );
        assert_eq!(
            support_feature_for(&WITH_ADDED_3D_CONFORMER_SPEC).map(|feature| feature.name),
            Some(crate::COORDINATE_EDIT_FEATURE.name)
        );
        assert_eq!(
            support_feature_for(&WITH_ONLY_3D_CONFORMER_SPEC).map(|feature| feature.name),
            Some(crate::COORDINATE_EDIT_FEATURE.name)
        );
        assert_eq!(
            WITH_2D_COORDINATE_BLOCK_SPEC.domain,
            OperationDomain::Coordinate
        );
        assert_eq!(WITH_3D_COORDINATES_SPEC.domain, OperationDomain::Coordinate);
        assert_eq!(
            WITH_CLEARED_3D_CONFORMERS_SPEC.domain,
            OperationDomain::Coordinate
        );
        assert_eq!(
            WITH_ADDED_3D_CONFORMER_SPEC.domain,
            OperationDomain::Coordinate
        );
        assert_eq!(
            WITH_ONLY_3D_CONFORMER_SPEC.domain,
            OperationDomain::Coordinate
        );
        assert!(support_matrix_contains(&WITH_2D_COORDINATE_BLOCK_SPEC));
        assert!(support_matrix_contains(&WITH_3D_COORDINATES_SPEC));
        assert!(support_matrix_contains(&WITH_CLEARED_3D_CONFORMERS_SPEC));
        assert!(support_matrix_contains(&WITH_ADDED_3D_CONFORMER_SPEC));
        assert!(support_matrix_contains(&WITH_ONLY_3D_CONFORMER_SPEC));
        assert!(invariant_matrix_contains(&WITH_2D_COORDINATE_BLOCK_SPEC));
        assert!(invariant_matrix_contains(&WITH_3D_COORDINATES_SPEC));
        assert!(invariant_matrix_contains(&WITH_CLEARED_3D_CONFORMERS_SPEC));
        assert!(invariant_matrix_contains(&WITH_ADDED_3D_CONFORMER_SPEC));
        assert!(invariant_matrix_contains(&WITH_ONLY_3D_CONFORMER_SPEC));
        assert!(parity_matrix_contains(&WITH_2D_COORDINATE_BLOCK_SPEC));
        assert!(parity_matrix_contains(&WITH_3D_COORDINATES_SPEC));
        assert!(parity_matrix_contains(&WITH_CLEARED_3D_CONFORMERS_SPEC));
        assert!(parity_matrix_contains(&WITH_ADDED_3D_CONFORMER_SPEC));
        assert!(parity_matrix_contains(&WITH_ONLY_3D_CONFORMER_SPEC));
    }

    #[test]
    fn conformer_generation_operation_registry_exposes_3d_coordinate_generation() {
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_3D_CONFORMER_SPEC))
        );
        assert!(
            MOLECULE_OPS
                .iter()
                .any(|operation| same_operation(operation, &WITH_3D_CONFORMERS_SPEC))
        );
        assert_eq!(
            support_feature_for(&WITH_3D_CONFORMER_SPEC).map(|feature| feature.name),
            Some(crate::CONFORMER_GENERATION_FEATURE.name)
        );
        assert_eq!(
            support_feature_for(&WITH_3D_CONFORMERS_SPEC).map(|feature| feature.name),
            Some(crate::CONFORMER_GENERATION_FEATURE.name)
        );
        assert_eq!(WITH_3D_CONFORMER_SPEC.domain, OperationDomain::Coordinate);
        assert_eq!(WITH_3D_CONFORMERS_SPEC.domain, OperationDomain::Coordinate);
        assert!(support_matrix_contains(&WITH_3D_CONFORMER_SPEC));
        assert!(support_matrix_contains(&WITH_3D_CONFORMERS_SPEC));
        assert!(invariant_matrix_contains(&WITH_3D_CONFORMER_SPEC));
        assert!(invariant_matrix_contains(&WITH_3D_CONFORMERS_SPEC));
        assert!(parity_matrix_contains(&WITH_3D_CONFORMER_SPEC));
        assert!(parity_matrix_contains(&WITH_3D_CONFORMERS_SPEC));

        let molecule = crate::Molecule::from_smiles("CC").expect("ethane");
        let generated = molecule
            .with_3d_conformer()
            .expect("default ETKDGv3 conformer");
        assert!(molecule.conformers_3d().is_empty());
        assert_eq!(generated.conformers_3d().len(), 1);

        let mut params = crate::EmbedParameters::etkdg();
        params.random_seed = 0xf00d;
        params.num_threads = 1;
        let generated_multi = molecule
            .with_3d_conformers_with_params(2, params)
            .expect("multi conformer generation");
        assert_eq!(generated_multi.conformers_3d().len(), 2);
    }

    #[test]
    fn sanitized_all_runs_through_operation_pipeline_without_changing_source() {
        let molecule = crate::Molecule::new();
        let original = molecule.clone();

        let sanitized = molecule.sanitize().unwrap();
        let sanitized_with_all = molecule.sanitize_with_ops(crate::SanitizeOps::ALL).unwrap();

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
    fn with_2d_coordinates_with_params_preserves_source_and_uses_parameterized_surface() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let result = molecule
            .with_2d_coordinates_with_params(crate::With2DCoordinatesParams {
                force_rdkit: true,
                use_ring_templates: true,
                ..crate::With2DCoordinatesParams::default()
            })
            .unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.conformers_2d().len(), 1);
        assert_eq!(
            result.source_coordinate_dim(),
            Some(crate::CoordinateDimension::TwoD)
        );
    }

    #[test]
    fn three_d_conformer_clear_and_single_assignment_run_through_coordinate_ops() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let first = vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]];
        let second = vec![[0.0, 0.0, 1.0], [1.2, 0.0, 1.0]];
        let replacement = vec![[0.5, 0.5, 0.5], [1.7, 0.5, 0.5]];

        let with_two = molecule
            .with_added_3d_conformer(first.clone(), true)
            .unwrap()
            .with_added_3d_conformer(second.clone(), true)
            .unwrap();
        let cleared = with_two.with_cleared_3d_conformers().unwrap();
        let only = with_two
            .with_only_3d_conformer(replacement.clone(), true)
            .unwrap();

        assert_eq!(with_two.conformers_3d().len(), 2);
        assert!(cleared.conformers_3d().is_empty());
        assert_eq!(only.conformers_3d().len(), 1);
        assert_eq!(only.conformers_3d()[0].coordinates(), replacement);
        assert_eq!(
            only.source_coordinate_dim(),
            Some(crate::CoordinateDimension::ThreeD)
        );

        let mut in_place = with_two.clone();
        in_place.clear_3d_conformers_().unwrap();
        assert!(in_place.conformers_3d().is_empty());
        in_place
            .set_only_3d_conformer_(second.clone(), true)
            .unwrap();
        assert_eq!(in_place.conformers_3d().len(), 1);
        assert_eq!(in_place.conformers_3d()[0].coordinates(), second);

        let two_d_flagged = with_two
            .with_only_3d_conformer(replacement.clone(), false)
            .unwrap();
        assert_eq!(
            two_d_flagged.source_coordinate_dim(),
            Some(crate::CoordinateDimension::TwoD)
        );
        assert!(!two_d_flagged.conformers_3d()[0].is_3d());
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
            result.coordinates_2d(),
            Some(&[[0.5, 1.0], [0.5, 1.0], [0.5, 1.0], [0.5, 1.0], [0.5, 1.0]][..])
        );
    }

    #[test]
    fn add_hydrogens_in_place_matches_value_operation() {
        let molecule = crate::Molecule::from_smiles("CCO").unwrap();
        let expected = molecule.with_hydrogens().unwrap();

        let mut in_place = molecule.clone();
        in_place.add_hydrogens_().unwrap();

        assert_eq!(in_place, expected);
    }

    #[test]
    fn add_hydrogens_in_place_preserves_shared_source_value() {
        let mut molecule = crate::Molecule::from_smiles("CCO").unwrap();
        let shared = molecule.clone();

        molecule.add_hydrogens_().unwrap();

        assert_eq!(shared.num_atoms(), 3);
        assert_eq!(shared.num_bonds(), 2);
        assert_eq!(molecule.num_atoms(), 9);
        assert_eq!(molecule.num_bonds(), 8);
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
            molecule.coordinates_2d().unwrap(),
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
            molecule.conformers_3d()[0].coordinates(),
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
            molecule.conformers_3d()[0].coordinates(),
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
            molecule.conformers_3d()[0].coordinates(),
            0,
            0,
            true,
        )
        .unwrap_err();
        assert!(degenerate.to_string().contains("degenerate atoms"));

        let nonterminal = add_hs_set_terminal_atom_coord(
            MoleculeReadParts::from_molecule(&molecule),
            &adjacency,
            molecule.conformers_3d()[0].coordinates(),
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
        let coords = result.coordinates_2d().unwrap();
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
    fn with_hydrogens_clears_atom_cip_ranks_like_rdkit_addhs() {
        let smiles =
            "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12";
        let molecule = crate::Molecule::from_smiles(smiles).unwrap();
        assert!(
            molecule
                .atoms()
                .iter()
                .any(|atom| atom.prop("_CIPRank").is_some()),
            "SMILES sanitize path should assign legacy _CIPRank before AddHs"
        );

        let result = molecule
            .with_hydrogens_with_params(crate::AddHsParams {
                only_on_atoms: Some(
                    [4usize, 6, 8, 10, 14, 16, 23]
                        .into_iter()
                        .map(crate::AtomId::new)
                        .collect(),
                ),
                ..crate::AddHsParams::default()
            })
            .unwrap();

        assert!(
            result
                .atoms()
                .iter()
                .all(|atom| atom.prop("_CIPRank").is_none()),
            "RDKit AddHs clears atom _CIPRank computed props before depiction"
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
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC).unwrap();
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
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC).unwrap();
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
        assert_eq!(result.coordinates_2d(), Some(&[[0.0, 0.0]][..]));
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

        let result = molecule.sanitize_with_ops(ops).unwrap();

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
            .sanitize_with_ops(
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
            .sanitize_with_ops(crate::SanitizeOps::SYMMRINGS | crate::SanitizeOps::KEKULIZE)
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
            .sanitize_with_ops(crate::SanitizeOps::KEKULIZE)
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
            .sanitize_with_ops(crate::SanitizeOps::KEKULIZE)
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
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::KEKULIZE)
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
            .sanitize_with_ops(
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
            |step, source| sanitize_valence_error(&SANITIZE_SPEC, step, source),
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
            .sanitize_with_ops(crate::SanitizeOps::NONE)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_ATROPISOMERS)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_ATROPISOMERS)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
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
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
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
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_CONJUGATION)
            .unwrap();

        assert!(result.bonds().iter().all(crate::Bond::is_conjugated));
    }

    #[test]
    fn sanitized_set_conjugation_keeps_aromatic_bonds_conjugated_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::SET_CONJUGATION)
            .unwrap();

        assert!(result.bonds().iter().all(crate::Bond::is_conjugated));
    }

    #[test]
    fn sanitized_set_conjugation_uses_heteroatom_lone_pair_candidate_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("NC=O", false).unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_CONJUGATION)
            .unwrap();

        assert_eq!(result.num_bonds(), 2);
        assert!(result.bonds().iter().all(crate::Bond::is_conjugated));
    }

    #[test]
    fn sanitized_sets_hybridization_after_conjugation_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

        let result = molecule
            .sanitize_with_ops(
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
            .sanitize_with_ops(
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
            .sanitize_with_ops(
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
            .sanitize_with_ops(
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
            .sanitize_with_ops(
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
            .sanitize_with_ops(
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
            .sanitize_with_ops(crate::SanitizeOps::FIND_RADICALS)
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
            .sanitize_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
            .unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.num_atoms(), molecule.num_atoms());
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_materializes_disappearing_pyrrolic_hydrogen() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("N1C=CC=C1", false).unwrap();

        let result = molecule
            .sanitize_with_ops(
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
            .sanitize_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
            .unwrap();

        assert_eq!(result.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_leaves_stable_explicit_hydrogens_unchanged_like_rdkit()
    {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
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
        assert_eq!(result.coordinates_2d(), original.coordinates_2d());
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
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();
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
        let mut parts = OpParts::new(&molecule, &SANITIZE_SPEC).unwrap();
        parts.set_valence_cache(crate::ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        });
        parts.mark_aromaticity_valid();

        parts.clear_cache(SANITIZE_SPEC.needs_update());
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
        let mut parts = OpParts::new(&molecule, &SANITIZE_SPEC).unwrap();

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
        let parts = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC).unwrap();

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
        let mut parts = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC).unwrap();

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
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC).unwrap();
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
        let mut parts = OpParts::new(&molecule, &WITH_HYDROGENS_SPEC).unwrap();
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
        let mut parts = OpParts::new(&molecule, &TEST_RECOMPUTE_VALENCE_SPEC).unwrap();

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
        let mut parts = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC).unwrap();
        parts.clear_cache(DerivedState::VALENCE);
    }

    #[cfg(not(feature = "op-contracts"))]
    #[test]
    fn op_contract_checks_are_disabled_without_feature() {
        let molecule = crate::Molecule::new();

        let mut unauthorized = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC).unwrap();
        unauthorized.clear_cache(DerivedState::VALENCE);
        unauthorized
            .finish(OpOutcome::Changed)
            .expect("without op-contracts, cache permission checks are disabled");

        let missing_update = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC).unwrap();
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
        let mut parts = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC).unwrap();
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
        let mut parts = OpParts::new(&molecule, &TEST_RECOMPUTE_VALENCE_SPEC).unwrap();
        parts.clear_cache(DerivedState::RINGS);
    }

    #[test]
    fn op_parts_cow_mutation_changes_result_without_changing_source() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();

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

        let mut parts = OpParts::new(&molecule, &WITHOUT_HYDROGENS_SPEC).unwrap();
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
        assert_eq!(result.coordinates_2d().unwrap(), &[[0.0, 0.0], [2.0, 0.0]]);
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

        let mut parts = OpParts::new(&molecule, &WITHOUT_HYDROGENS_SPEC).unwrap();
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

    #[cfg(feature = "op-contracts")]
    #[test]
    fn compacting_topology_edit_record_is_rejected_for_weak_operations() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();
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
