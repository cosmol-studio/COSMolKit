//! Macro-controlled molecule operations.
//!
//! This module is the canonical operation-contract surface for molecular state
//! migration. Single-output operation bodies use
//! `#[mol_op_body(op_name, parts)]`; multiple-output bodies use
//! `#[mol_multi_op_body(op_name, parts)]`. Both attribute macros inject the
//! only allowed capability object while preserving declared operation
//! parameters.
//!
//! Agent guardrail: operation bodies must not take `&mut Molecule`, must not
//! access `Molecule` storage directly, and must not update derived-state caches
//! or invalidation flags by hand. If a future change appears to require bypassing
//! those rules, the agent is not allowed to continue silently. It must stop and
//! confirm the design exception with the human author before editing code.
//!
//! RDKit marker convention is defined in `dev/source_reproduction_protocol.md`.

use std::fmt;

#[allow(unused_imports)]
use cosmolkit_macros::mol_multi_op_body;
use cosmolkit_macros::{mol_op_body, molecule_ops};

use crate::{
    AtomId, DerivedState, InvariantError, Molecule, MoleculeProperties, SupportStatus, TopologyTrust,
    invariants::enforce_molecule_invariants,
    molecule::{CoordinateBlock, TopologyBlock, TopologyMapping},
};

pub(crate) use crate::read_parts::MoleculeReadParts;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoleculeOpKind {
    Strong,
    Weak,
}

/// Number of molecule values produced by a registered operation wrapper.
///
/// Multiple-output operations still validate every emitted molecule through
/// an independent `OpParts` lifecycle; this field does not relax mutation or
/// derived-state authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoleculeOpOutput {
    Single,
    Multiple,
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

/// Source-defined lifecycle of modern CIP observable state for an operation.
///
/// This is separate from `DerivedState::STEREO`: `_CIPCode` is deliberately a
/// non-computed RDKit property, while `_CIPNeighborOrder` and `_CIPComputed`
/// are computed properties with different clearing behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipStatePolicy {
    /// Preserve all modern and legacy CIP properties byte-for-byte.
    Preserve,
    /// Reproduce `ROMol::clearComputedProps()`: preserve `_CIPCode`, clear
    /// `_CIPNeighborOrder`, `_CIPComputed`, and legacy `_CIPRank`.
    ClearComputed,
    /// The operation is the sole modern CIP assignment producer.
    Assign,
    /// Reproduce the selective source-defined CIP/stereo lifecycle of
    /// tautomer enumeration. This policy is allow-listed to that operation;
    /// it grants no block access or mutation capability by itself.
    TautomerSourceTransition,
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
        self.read.union(self.write).contains(block)
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
    operation_defined: DerivedState,
}

impl DerivedEffects {
    pub const NONE: Self = Self {
        recompute: DerivedState::NONE,
        preserve: DerivedState::NONE,
        invalidate: DerivedState::NONE,
        operation_defined: DerivedState::NONE,
    };

    #[must_use]
    pub const fn new(
        recompute: DerivedState,
        preserve: DerivedState,
        invalidate: DerivedState,
        operation_defined: DerivedState,
    ) -> Self {
        Self {
            recompute,
            preserve,
            invalidate,
            operation_defined,
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

    /// States whose source-required transition is implemented by the
    /// operation body instead of one of the three standard mechanisms.
    ///
    /// This category grants no read access. The registry guardrail currently
    /// permits it only for valence in the hydrogen-removal operation family.
    #[must_use]
    pub const fn operation_defined(self) -> DerivedState {
        self.operation_defined
    }

    #[must_use]
    pub const fn needs_update(self) -> DerivedState {
        self.invalidate.union(self.recompute).union(self.operation_defined)
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
}

impl fmt::Display for SemanticPrecondition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TrustedBondTopology => f.write_str("trusted_bond_topology"),
            Self::HydrogenOwnershipRepresented => f.write_str("hydrogen_ownership_represented"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SemanticPreconditionSet(u8);

impl SemanticPreconditionSet {
    pub const NONE: Self = Self(0);
    pub const TRUSTED_BOND_TOPOLOGY: Self = Self(1 << 0);
    pub const HYDROGEN_OWNERSHIP_REPRESENTED: Self = Self(1 << 1);

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
    pub output: MoleculeOpOutput,
    pub result_type: &'static str,
    pub domain: OperationDomain,
    pub kind: MoleculeOpKind,
    pub topology_edit: TopologyEditKind,
    pub access: BlockAccess,
    pub may_mutate: BlockSet,
    pub auto_remap: BlockSet,
    pub derived_effects: DerivedEffects,
    pub cip_state: CipStatePolicy,
    pub semantic_preconditions: SemanticPreconditionSet,
    pub requires_mapping: MappingRequirement,
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
    #[error("{operation}: molecule transform failed: {source}")]
    MolTransform {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::mol_transforms::MolTransformError,
    },
    #[error("{operation}: molecule alignment failed: {source}")]
    Alignment {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::AlignmentError,
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
    #[error("{operation}: stereo error: {source}")]
    Stereo {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::StereoError,
    },
    #[error("{operation}: modern CIP labeling error: {source}")]
    CipLabeler {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: crate::CipLabelerError,
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
    #[error("{operation}: tautomer enumeration failed: {source}")]
    Tautomer {
        operation: &'static MoleculeOpSpec,
        #[source]
        source: Box<crate::chemistry::tautomer::TautomerRunError>,
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
            recompute: [valence],
            preserve: [rings, ring_families],
            invalidate: [aromaticity, stereo, drawing, fingerprint],
        },
        cip_state: clear_computed,
        semantic_preconditions: [
            trusted_bond_topology,
        ],
        requires_mapping: required,
        feature: HYDROGENS_FEATURE,
        parity: required_now,
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
            recompute: [aromaticity, rings],
            preserve: [],
            invalidate: [ring_families, stereo, drawing, fingerprint],
            operation_defined: [valence],
        },
        cip_state: clear_computed,
        requires_mapping: required,
        feature: HYDROGENS_FEATURE,
        parity: required_now,
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
            recompute: [aromaticity, rings],
            preserve: [],
            invalidate: [ring_families, stereo, drawing, fingerprint],
            operation_defined: [valence],
        },
        cip_state: clear_computed,
        requires_mapping: required,
        feature: HYDROGENS_FEATURE,
        parity: required_now,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: KEKULIZE_FEATURE,
        parity: required_now,
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
        access: { read: [], write: [topology, properties, derived_cache] },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings, valence, aromaticity],
            preserve: [],
            invalidate: [ring_families, stereo, drawing, fingerprint],
        },
        cip_state: clear_computed,
        requires_mapping: none,
        feature: SANITIZE_FEATURE,
        parity: required_now,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "sanitize_default_rdkit",
    }

    op assigned_valence(strict: bool) {
        method: with_assigned_valence_strict,
        impl_fn: assigned_valence_impl,
        default_method: with_assigned_valence,
        default_args: [true],
        inplace: true,
        inplace_method: assign_valence_strict_,
        default_inplace_method: assign_valence_,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: VALENCE_FEATURE,
        parity: required_now,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: RINGS_FEATURE,
        parity: required_now,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: RINGS_FEATURE,
        parity: required_now,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: AROMATICITY_FEATURE,
        parity: required_now,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: VALENCE_FEATURE,
        parity: required_now,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "assign_radicals_rdkit",
    }

    op with_chiral_tags_from_structure(conformer_id: i32, replace_existing_tags: bool) {
        method: with_chiral_tags_from_structure,
        docs: "Return a new molecule with atom chiral tags assigned from a 3D conformer. A negative conformer id selects the first conformer. Existing tags are replaced only when `replace_existing_tags` is true. Atom and bond ordering and all conformer coordinates are preserved. This stable operation has exact full-state parity with pinned RDKit 2026.03.1 `assignChiralTypesFrom3D` across all 77 fixed oracle records for tetrahedral C/S/Se and environment-enabled square-planar, trigonal-bipyramidal, and octahedral centers, including properties, no-ops, and defined errors. It does not perform `assignStereochemistryFrom3D`, 3D double-bond direction/E-Z assignment, CIP orchestration, or distinct-substituent validation.",
        impl_fn: with_chiral_tags_from_structure_impl,
        inplace: true,
        inplace_method: assign_chiral_tags_from_structure_,
        inplace_docs: "Assign atom chiral tags from a selected 3D conformer in place. This stable mutating form of `Molecule::with_chiral_tags_from_structure` has the same pinned-RDKit parity scope; failures are transactional and leave the molecule unchanged.",
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [coordinates], write: [topology, properties, derived_cache] },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [stereo, drawing, fingerprint],
        },
        cip_state: preserve,
        semantic_preconditions: [trusted_bond_topology],
        requires_mapping: none,
        feature: STEREO_FEATURE,
        parity: required_now,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "assign_atom_chiral_tags_from_structure_rdkit",
    }

    op assigned_cip_labels(options: crate::CipLabelOptions) {
        method: with_cip_labels_with_options,
        docs: "Return a new molecule with modern molecular-context CIP labels assigned according to `options`. When both selections are absent or empty, all eligible atoms and bonds are labeled, matching the pinned wrapper's truth-value dispatch. Once either selection is non-empty, an absent or empty category selects no entries. The source molecule, topology indexing, and coordinates are preserved.",
        impl_fn: assigned_cip_labels_impl,
        default_method: with_cip_labels,
        default_args: [crate::CipLabelOptions::default()],
        inplace: true,
        inplace_method: assign_cip_labels_with_options_,
        default_inplace_method: assign_cip_labels_,
        inplace_docs: "Assign modern molecular-context CIP labels in place. Errors leave internally complete molecule storage; callers must not assume rollback of source-ordered CIP state.",
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [], write: [topology, properties, derived_cache] },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [stereo, drawing, fingerprint],
        },
        cip_state: assign,
        semantic_preconditions: [trusted_bond_topology],
        requires_mapping: none,
        feature: CIP_LABELER_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "weak_topology_state",
        parity_profile: "modern_ciplabeler_rdkit",
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
        access: { read: [topology, properties], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        semantic_preconditions: [trusted_bond_topology],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: COORDINATE_2D_FEATURE,
        parity: required_now,
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
        access: { read: [topology, properties, derived_cache], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        semantic_preconditions: [trusted_bond_topology],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: CONFORMER_GENERATION_FEATURE,
        parity: required_now,
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
        access: { read: [topology, properties, derived_cache], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        semantic_preconditions: [trusted_bond_topology],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: CONFORMER_GENERATION_FEATURE,
        parity: required_now,
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
        cip_state: preserve,
        requires_mapping: none,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_set_3d",
        parity_profile: "set_3d_coordinates_manual",
    }

    op with_alignment_to(reference: &crate::Molecule, params: &crate::AlignmentParameters) {
        method: with_alignment_to,
        impl_fn: with_alignment_to_impl,
        result_type: crate::AlignmentResult,
        inplace: true,
        inplace_method: align_to_,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [topology, properties, derived_cache], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [drawing],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: MOLALIGN_FEATURE,
        parity: required_now,
        io_roundtrip: true,
        invariant_profile: "coordinate_align_molecule",
        parity_profile: "align_mol_rdkit",
    }

    op with_aligned_conformers(params: crate::ConformerAlignmentParameters) {
        method: with_aligned_conformers_with_params,
        impl_fn: with_aligned_conformers_impl,
        result_type: crate::ConformerAlignmentReport,
        default_method: with_aligned_conformers,
        default_args: [crate::ConformerAlignmentParameters::default()],
        inplace: true,
        inplace_method: align_conformers_with_params_,
        default_inplace_method: align_conformers_,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: MOLALIGN_FEATURE,
        parity: required_now,
        io_roundtrip: true,
        invariant_profile: "coordinate_align_conformers",
        parity_profile: "align_mol_conformers_rdkit",
    }

    op with_atom_position(atom: usize, position: [f64; 3], conformer_index: usize) {
        method: with_atom_position,
        impl_fn: with_atom_position_impl,
        inplace: true,
        inplace_method: set_atom_position_,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_set_atom_position",
        parity_profile: "set_atom_position_rdkit",
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
        cip_state: preserve,
        requires_mapping: none,
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
        cip_state: preserve,
        requires_mapping: none,
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
        cip_state: preserve,
        requires_mapping: none,
        feature: COORDINATE_EDIT_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "coordinate_set_single_3d_conformer",
        parity_profile: "remove_all_then_add_conformer_manual",
    }

    op enumerate_tautomers_with_options(enumerator: &crate::chemistry::tautomer::TautomerEnumerator<'_>) {
        method: enumerate_tautomers_with_options,
        impl_fn: enumerate_tautomers_with_options_impl,
        default_method: enumerate_tautomers,
        default_args: [&crate::chemistry::tautomer::TautomerEnumerator::new()],
        output: multiple,
        result_type: crate::chemistry::tautomer::TautomerEnumeration,
        assemble_fn: assemble_tautomer_enumeration,
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: { read: [coordinates], write: [topology, properties, derived_cache] },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings, ring_families, valence, aromaticity, stereo],
            preserve: [],
            invalidate: [drawing, fingerprint],
        },
        cip_state: tautomer_source_transition,
        requires_mapping: none,
        feature: TAUTOMER_ENUMERATION_FEATURE,
        parity: required_when_supported,
        io_roundtrip: true,
        invariant_profile: "tautomer_enumeration_multi_output",
        parity_profile: "rdkit_tautomer_enumeration_full_plan",
        docs: "Enumerate ordered tautomers with the configured source-aligned enumerator.",
    }
}

mod bodies;
pub(crate) mod hydrogens;
mod runtime;
mod sanitize_pipeline;
mod tautomer;

// These names are resolved by macro-expanded multiple-output wrappers and
// bodies once the registry contains its first multiple-output operation.
pub(crate) use self::runtime::MultiMoleculeOpParts;
use self::runtime::OpParts;
pub use self::runtime::OperationTrace;
use self::{bodies::*, hydrogens::*, sanitize_pipeline::*, tautomer::*};

pub(crate) use self::sanitize_pipeline::{sanitize_conjugation_assignment, sanitize_hybridization_assignment};

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::{BondOrder, BondQueryPredicate, QueryNode};

    mod generated_multiple_output_wrapper {
        use super::*;

        #[cosmolkit_macros::mol_multi_op_body(test_generated_multiple_output, parts)]
        fn test_generated_multiple_output_impl() -> Result<usize, OperationError> {
            let first = parts.derive_from_source(|_branch| Ok(()))?;
            parts.emit(first)?;
            let second = parts.derive_from_branch(first, |branch| {
                branch.with_topology_mut(|_branch, topology| {
                    topology.atoms[0].set_formal_charge(1);
                    Ok(())
                })?;
                branch.record_topology_edit(TopologyEditKind::Local)?;
                branch.clear_cache(TEST_GENERATED_MULTIPLE_OUTPUT_SPEC.needs_update());
                Ok(())
            })?;
            parts.emit(second)?;
            Ok(2)
        }

        fn assemble_test_generated_multiple_output(
            molecules: Vec<Molecule>,
            count: usize,
        ) -> Result<(Vec<Molecule>, usize), OperationError> {
            Ok((molecules, count))
        }

        cosmolkit_macros::molecule_ops! {
            op test_generated_multiple_output() {
                method: test_generated_multiple_output,
                impl_fn: test_generated_multiple_output_impl,
                output: multiple,
                result_type: (Vec<crate::Molecule>, usize),
                assemble_fn: assemble_test_generated_multiple_output,
                domain: topology,
                kind: weak,
                topology_edit: local,
                access: { read: [], write: [topology, derived_cache] },
                may_mutate: [topology, derived_cache],
                auto_remap: [],
                derived_effects: {
                    recompute: [],
                    preserve: [],
                    invalidate: [valence, aromaticity, stereo, drawing, fingerprint],
                },
                cip_state: preserve,
                requires_mapping: none,
                feature: HYDROGENS_FEATURE,
                parity: not_applicable,
                io_roundtrip: false,
                invariant_profile: "test_multiple_output",
            }
        }

        #[test]
        fn generated_multiple_output_wrapper_validates_and_assembles_typed_results() {
            let source = crate::Molecule::from_smiles("C").unwrap();

            let (molecules, count) = source.test_generated_multiple_output().unwrap();

            assert_eq!(count, 2);
            assert_eq!(molecules.len(), 2);
            assert_eq!(molecules[0].atoms()[0].formal_charge(), 0);
            assert_eq!(molecules[1].atoms()[0].formal_charge(), 1);
            assert_eq!(source.atoms()[0].formal_charge(), 0);
            assert_eq!(TEST_GENERATED_MULTIPLE_OUTPUT_SPEC.output, MoleculeOpOutput::Multiple);
            assert_eq!(
                TEST_GENERATED_MULTIPLE_OUTPUT_SPEC.result_type,
                "(Vec < crate :: Molecule > , usize)"
            );
        }
    }

    const TEST_NEEDS_VALENCE_UPDATE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_needs_valence_update",
        impl_fn: "test_needs_valence_update_impl",
        output: MoleculeOpOutput::Single,
        result_type: "Molecule",
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
            DerivedState::NONE,                               // operation_defined
        ),
        cip_state: CipStatePolicy::Preserve,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_RECOMPUTE_VALENCE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_recompute_valence",
        impl_fn: "test_recompute_valence_impl",
        output: MoleculeOpOutput::Single,
        result_type: "Molecule",
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
            DerivedState::NONE,    // operation_defined
        ),
        cip_state: CipStatePolicy::Preserve,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_OVERLAPPING_ACCESS_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_overlapping_access",
        impl_fn: "test_overlapping_access_impl",
        output: MoleculeOpOutput::Single,
        result_type: "Molecule",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::Local,
        access: BlockAccess::new(BlockSet::TOPOLOGY, BlockSet::TOPOLOGY),
        may_mutate: BlockSet::TOPOLOGY,
        auto_remap: BlockSet::NONE,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        derived_effects: DerivedEffects::NONE,
        cip_state: CipStatePolicy::Preserve,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_MULTIPLE_OUTPUT_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "test_multiple_output",
        impl_fn: "test_multiple_output_impl",
        output: MoleculeOpOutput::Multiple,
        result_type: "Vec<Molecule>",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::Local,
        access: BlockAccess::new(BlockSet::NONE, BlockSet::TOPOLOGY.union(BlockSet::DERIVED_CACHE)),
        may_mutate: BlockSet::TOPOLOGY.union(BlockSet::DERIVED_CACHE),
        auto_remap: BlockSet::NONE,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        derived_effects: DerivedEffects::new(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::VALENCE
                .union(DerivedState::AROMATICITY)
                .union(DerivedState::STEREO)
                .union(DerivedState::DRAWING)
                .union(DerivedState::FINGERPRINT),
            DerivedState::NONE,
        ),
        cip_state: CipStatePolicy::Preserve,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    const TEST_TAUTOMER_CIP_SPEC: MoleculeOpSpec = MoleculeOpSpec {
        method: "enumerate_tautomers_with_options",
        impl_fn: "test_tautomer_cip_impl",
        output: MoleculeOpOutput::Multiple,
        result_type: "Vec<Molecule>",
        domain: OperationDomain::Topology,
        kind: MoleculeOpKind::Weak,
        topology_edit: TopologyEditKind::Local,
        access: BlockAccess::new(
            BlockSet::NONE,
            BlockSet::TOPOLOGY
                .union(BlockSet::PROPERTIES)
                .union(BlockSet::DERIVED_CACHE),
        ),
        may_mutate: BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE),
        auto_remap: BlockSet::NONE,
        semantic_preconditions: SemanticPreconditionSet::NONE,
        derived_effects: DerivedEffects::new(
            DerivedState::NONE,
            DerivedState::NONE,
            DerivedState::VALENCE
                .union(DerivedState::RINGS)
                .union(DerivedState::RING_FAMILIES)
                .union(DerivedState::AROMATICITY)
                .union(DerivedState::STEREO)
                .union(DerivedState::DRAWING)
                .union(DerivedState::FINGERPRINT),
            DerivedState::NONE,
        ),
        cip_state: CipStatePolicy::TautomerSourceTransition,
        requires_mapping: MappingRequirement::None,
        support: SupportStatus::Experimental,
        parity: ParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    #[test]
    fn molecule_read_parts_does_not_expose_raw_molecule_escape() {
        let ops_source = include_str!("ops.rs");
        let read_parts_source = include_str!("../model/read_parts.rs");
        assert!(!ops_source.contains(concat!("Molecule", "ReadAccess")));
        assert!(!ops_source.contains(concat!("read_parts", ".", "molecule")));
        assert!(!ops_source.contains(concat!(".", "molecule", "()")));
        assert!(!read_parts_source.contains(concat!("with_", "molecule", "_read")));
        assert!(!read_parts_source.contains(concat!("pub(crate) fn ", "molecule")));
    }

    fn op_body_sources() -> [(&'static str, &'static str); 4] {
        [
            ("bodies.rs", include_str!("ops/bodies.rs")),
            ("hydrogens.rs", include_str!("ops/hydrogens.rs")),
            ("sanitize_pipeline.rs", include_str!("ops/sanitize_pipeline.rs")),
            ("tautomer.rs", include_str!("ops/tautomer.rs")),
        ]
    }

    fn assert_sources_do_not_contain(sources: &[(&str, &str)], forbidden: &[&str]) {
        for (name, source) in sources {
            for pattern in forbidden {
                assert!(
                    !source.contains(pattern),
                    "{name} must not contain forbidden operation-body pattern `{pattern}`"
                );
            }
        }
    }

    fn parts_method_names(source: &str) -> Vec<&str> {
        let mut names = Vec::new();
        let mut remaining = source;
        while let Some(index) = remaining.find("parts.") {
            let after = &remaining[index + "parts.".len()..];
            let len = after
                .find(|ch: char| !(ch.is_ascii_alphanumeric() || ch == '_'))
                .unwrap_or(after.len());
            if len > 0 {
                names.push(&after[..len]);
            }
            remaining = &after[len..];
        }
        names
    }

    #[test]
    fn op_parts_working_state_is_runtime_private() {
        let ops_source = include_str!("ops.rs");
        let runtime_source = include_str!("ops/runtime.rs");

        assert!(ops_source.contains("mod runtime;"));
        assert!(runtime_source.contains("pub struct OpParts"));
        assert!(runtime_source.contains("working: Molecule"));
        assert!(ops_source.contains(concat!("use self::runtime::", "OpParts;")));
        assert!(!ops_source.contains(concat!("pub(crate) use self::runtime::", "OpParts")));
        assert!(!ops_source.contains(concat!("pub use self::runtime::", "OpParts")));
    }

    #[test]
    fn operation_body_modules_do_not_access_runtime_or_raw_molecule_escape_hatches() {
        let sources = op_body_sources();
        assert_sources_do_not_contain(
            &sources,
            &[
                concat!("parts", ".", "working"),
                concat!(".", "working"),
                concat!("parts", ".", "spec"),
                concat!("&", "Molecule", ","),
                concat!("&", "Molecule", ")"),
                concat!("&", "Molecule", "\n"),
                concat!("Molecule", "::"),
                concat!("MoleculeReadParts", "::", "from_molecule"),
                concat!("from_", "operation", "_blocks"),
                concat!("read_parts", "_for_", "topology"),
                concat!("read_parts", "_for_", "blocks"),
                concat!("read_parts", "_for_", "optional", "_blocks"),
                concat!("OpParts", "::", "new"),
                concat!("new_", "in_place"),
            ],
        );
    }

    #[test]
    fn operation_body_parts_usage_stays_within_capability_methods() {
        const ALLOWED: &[&str] = &[
            "begin_operation_read",
            "begin_topology_read",
            "begin_topology_mut",
            "commit_topology",
            "begin_coordinates_mut",
            "commit_coordinates",
            "begin_properties_mut",
            "commit_properties",
            "with_topology_mut",
            "with_coordinates_mut",
            "with_topology_and_properties_mut",
            "with_topology_coordinates_properties_mut",
            "record_topology_edit",
            "record_topology_mapping",
            "clear_cache",
            "clear_computed_properties",
            "set_rings_cache",
            "set_ring_families_cache",
            "set_valence_cache",
            "mark_aromaticity_valid",
            "mark_stereo_handled",
            "prove_preserved",
            "with_topology_read_parts",
            "with_block_read_parts",
            "with_optional_block_read_parts",
            "with_borrowed_optional_block_read_parts",
            "with_coordinate_update_read_parts",
            "derive_from_source",
            "derive_from_branch",
            "with_source_read_parts",
            "with_branch_read_parts",
            "with_source_and_branch_read_parts",
            "with_source_and_branches_read_parts",
            "emit",
        ];

        for (source_name, source) in op_body_sources() {
            for method in parts_method_names(source) {
                assert!(
                    ALLOWED.contains(&method),
                    "{source_name} uses non-capability OpParts method `{method}`"
                );
            }
        }
    }

    #[test]
    fn conformer_alignment_operation_calls_the_alignment_core_once() {
        let molecule = crate::Molecule::from_smiles("CCC").expect("molecule graph");
        let mut builder = molecule.to_builder();
        builder
            .add_conformer(crate::Conformer3D::new(
                7,
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
                true,
            ))
            .expect("reference conformer");
        builder
            .add_conformer(crate::Conformer3D::new(
                17,
                vec![[3.0, 4.0, 5.0], [4.0, 4.0, 5.0], [3.0, 6.0, 5.0]],
                true,
            ))
            .expect("probe conformer");
        let molecule = builder.build().expect("multi-conformer molecule");

        crate::mol_align::reset_align_conformers_call_count();
        let (_aligned, report) = molecule.with_aligned_conformers().expect("value conformer alignment");
        assert_eq!(report.rmsds.len(), 1);
        assert_eq!(crate::mol_align::align_conformers_call_count(), 1);

        let mut inplace = molecule;
        crate::mol_align::reset_align_conformers_call_count();
        let report = inplace.align_conformers_().expect("in-place conformer alignment");
        assert_eq!(report.rmsds.len(), 1);
        assert_eq!(crate::mol_align::align_conformers_call_count(), 1);
    }

    #[test]
    fn fallible_operation_bodies_do_not_leave_raw_blocks_checked_out() {
        let sources = op_body_sources();
        assert_sources_do_not_contain(
            &sources,
            &[
                "parts.begin_topology_mut",
                "parts.commit_topology",
                "parts.begin_coordinates_mut",
                "parts.commit_coordinates",
                "parts.begin_properties_mut",
                "parts.commit_properties",
            ],
        );
    }

    #[test]
    fn scoped_in_place_error_returns_current_topology_block() {
        let mut molecule = crate::Molecule::from_smiles("C").unwrap();
        let mut parts = OpParts::new_in_place(&mut molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();

        let error = parts
            .with_topology_mut(|_parts, topology| {
                topology.atoms[0].set_formal_charge(1);
                Err::<(), _>(OperationError::InvalidInput {
                    operation: &WITH_KEKULIZED_BONDS_SPEC,
                    message: "synthetic failure after a partial topology update",
                })
            })
            .expect_err("the scoped mutation must return the operation error");
        assert!(matches!(error, OperationError::InvalidInput { .. }));

        parts.abort_in_place();
        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].formal_charge(), 1);
        enforce_molecule_invariants(&molecule)
            .expect("the current topology block must be structurally complete after error");
    }

    #[test]
    fn failed_in_place_sanitize_does_not_leak_checkout_placeholder() {
        let mut molecule = crate::Molecule::from_smiles_with_sanitize("c", false).unwrap();

        molecule
            .sanitize_()
            .expect_err("an isolated aromatic atom must fail kekulization");

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].element(), crate::Element::C);
        enforce_molecule_invariants(&molecule)
            .expect("failed sanitization must return a complete current topology block");
    }

    #[test]
    fn multiple_output_branches_preserve_source_and_validate_each_emitted_molecule() {
        let source = crate::Molecule::from_smiles("C").unwrap();
        let mut parts = MultiMoleculeOpParts::new(&source, &TEST_MULTIPLE_OUTPUT_SPEC).unwrap();

        let unchanged = parts
            .derive_from_source(|_branch| Ok(()))
            .expect("the unchanged source is a valid branch");
        parts.emit(unchanged).unwrap();

        let charged = parts
            .derive_from_branch(unchanged, |branch| {
                branch.with_topology_mut(|_branch, topology| {
                    topology.atoms[0].set_formal_charge(1);
                    Ok(())
                })?;
                branch.record_topology_edit(TopologyEditKind::Local)?;
                branch.clear_cache(TEST_MULTIPLE_OUTPUT_SPEC.needs_update());
                Ok(())
            })
            .expect("a child branch must complete its own contract lifecycle");

        let observed_charge = parts
            .with_branch_read_parts(charged, |read| Ok(read.atoms()[0].formal_charge()))
            .unwrap();
        assert_eq!(observed_charge, 1);
        parts.emit(charged).unwrap();

        let outputs = parts.finish().unwrap();
        assert_eq!(outputs.len(), 2);
        assert_eq!(outputs[0].atoms()[0].formal_charge(), 0);
        assert_eq!(outputs[1].atoms()[0].formal_charge(), 1);
        assert_eq!(source.atoms()[0].formal_charge(), 0);
        for output in &outputs {
            enforce_molecule_invariants(output).expect("every emitted branch must satisfy molecule invariants");
        }
    }

    #[test]
    fn multiple_output_executor_rejects_cardinality_and_unknown_branch_misuse() {
        let source = crate::Molecule::from_smiles("C").unwrap();

        assert!(matches!(
            MultiMoleculeOpParts::new(&source, &TEST_RECOMPUTE_VALENCE_SPEC),
            Err(OperationError::InvalidInput {
                message: "multiple-output capability used with a single-output operation",
                ..
            })
        ));
        assert!(matches!(
            OpParts::new(&source, &TEST_MULTIPLE_OUTPUT_SPEC),
            Err(OperationError::InvalidInput {
                message: "single-output OpParts used with a multiple-output operation",
                ..
            })
        ));

        let mut other = MultiMoleculeOpParts::new(&source, &TEST_MULTIPLE_OUTPUT_SPEC).unwrap();
        let foreign_branch = other.derive_from_source(|_branch| Ok(())).unwrap();
        let mut parts = MultiMoleculeOpParts::new(&source, &TEST_MULTIPLE_OUTPUT_SPEC).unwrap();
        parts
            .derive_from_source(|_branch| Ok(()))
            .expect("the local branch deliberately occupies the same index");
        assert!(matches!(
            parts.emit(foreign_branch),
            Err(OperationError::InvalidInput {
                message: "multiple-output operation used a branch handle from another execution",
                ..
            })
        ));
    }

    #[test]
    fn tautomer_multiple_output_spec_rejects_foreign_handles_after_source_and_child_derivation() {
        let source = crate::Molecule::from_smiles("CC(C)=O").unwrap();
        let mut foreign = MultiMoleculeOpParts::new(&source, &ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC).unwrap();
        let foreign_branch = foreign.derive_from_source(|_branch| Ok(())).unwrap();

        let mut parts = MultiMoleculeOpParts::new(&source, &ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC).unwrap();
        let source_branch = parts.derive_from_source(|_branch| Ok(())).unwrap();
        let child_branch = parts.derive_from_branch(source_branch, |_branch| Ok(())).unwrap();
        parts.emit(source_branch).unwrap();
        parts.emit(child_branch).unwrap();
        assert!(matches!(
            parts.emit(foreign_branch),
            Err(OperationError::InvalidInput {
                message: "multiple-output operation used a branch handle from another execution",
                ..
            })
        ));

        let outputs = parts.finish().unwrap();
        assert_eq!(outputs, [source.clone(), source]);
    }

    #[test]
    fn multiple_output_duplicate_emissions_preserve_requested_order() {
        let source = crate::Molecule::from_smiles("C").unwrap();
        let mut parts = MultiMoleculeOpParts::new(&source, &TEST_MULTIPLE_OUTPUT_SPEC).unwrap();
        let branch = parts.derive_from_source(|_branch| Ok(())).unwrap();

        parts.emit(branch).unwrap();
        parts.emit(branch).unwrap();

        let outputs = parts.finish().unwrap();
        assert_eq!(outputs.len(), 2);
        assert_eq!(outputs[0].to_smiles(true).unwrap(), "C");
        assert_eq!(outputs[1].to_smiles(true).unwrap(), "C");
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn multiple_output_branch_cannot_bypass_per_branch_contract_validation() {
        let source = crate::Molecule::from_smiles("C").unwrap();
        let mut parts = MultiMoleculeOpParts::new(&source, &TEST_MULTIPLE_OUTPUT_SPEC).unwrap();

        let error = parts
            .derive_from_source(|branch| {
                branch.with_topology_mut(|_branch, topology| {
                    topology.atoms[0].set_formal_charge(1);
                    Ok(())
                })?;
                Ok(())
            })
            .expect_err("an unrecorded topology edit must fail branch finalization");

        assert!(matches!(
            error,
            OperationError::InvalidInput {
                message: "operation body did not clear or update every required cache state",
                ..
            }
        ));
        assert!(parts.finish().unwrap().is_empty());
        assert_eq!(source.atoms()[0].formal_charge(), 0);
    }

    #[test]
    fn distgeom_operation_bodies_use_coordinate_update_boundary() {
        let bodies_source = include_str!("ops/bodies.rs");
        assert!(bodies_source.contains("embed_molecule_coordinate_update"));
        assert!(bodies_source.contains("embed_multiple_confs_coordinate_update"));
        assert!(!bodies_source.contains(concat!("embed_", "molecule", "(")));
        assert!(!bodies_source.contains(concat!("embed_", "multiple", "_confs", "(")));
        assert!(!bodies_source.contains(concat!("working", ".", "clone", "()")));
    }

    #[test]
    fn operation_body_functions_live_outside_runtime_module() {
        let runtime_source = include_str!("ops/runtime.rs");
        let bodies_source = include_str!("ops/bodies.rs");
        let hydrogens_source = include_str!("ops/hydrogens.rs");
        let sanitize_source = include_str!("ops/sanitize_pipeline.rs");
        let tautomer_source = include_str!("ops/tautomer.rs");

        assert!(!runtime_source.contains("#[mol_op_body"));
        assert!(!runtime_source.contains("#[mol_multi_op_body"));
        assert!(bodies_source.contains("#[mol_op_body"));
        assert!(hydrogens_source.contains("#[mol_op_body"));
        assert!(sanitize_source.contains("#[mol_op_body"));
        assert!(tautomer_source.contains("#[cosmolkit_macros::mol_multi_op_body"));
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn begin_topology_mut_rejects_second_begin() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();
        let _topology = parts.begin_topology_mut().expect("first topology begin should succeed");
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
        let _topology = parts.begin_topology_mut().expect("first topology begin should succeed");
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
    fn finish_rejects_a_writable_block_left_checked_out() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();
        let _topology = parts.begin_topology_mut().unwrap();

        let error = parts
            .finish()
            .expect_err("finish must reject an unclosed begin/commit lifecycle");

        assert!(matches!(
            error,
            OperationError::InvalidInput {
                message: "operation finished while a writable block was still checked out",
                ..
            }
        ));
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
        assert!(matches!(err, OperationError::InvalidInput { message, .. } if message.contains("both read and write")));
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

    fn support_feature_for(operation: &'static MoleculeOpSpec) -> Option<&'static crate::FeatureSpec> {
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
            other => panic!("expected UnsupportedFeature for {}, got {other:?}", operation.method),
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
    fn operation_defined_effect_is_limited_to_hydrogen_removal_valence() {
        let operation_defined = MOLECULE_OPS
            .iter()
            .copied()
            .filter(|operation| operation.derived_effects.operation_defined() != DerivedState::NONE)
            .collect::<Vec<_>>();

        assert_eq!(
            operation_defined
                .iter()
                .map(|operation| operation.method)
                .collect::<Vec<_>>(),
            vec!["without_hydrogens_with_sanitize", "without_hydrogens_with_params"],
            "OperationDefined is approved only for the two generated specs of the single hydrogen-removal operation family"
        );
        for operation in operation_defined {
            assert_eq!(operation.derived_effects.operation_defined(), DerivedState::VALENCE);
            assert!(operation.access.can_write(BlockSet::DERIVED_CACHE));
            assert_eq!(operation.parity, ParityPolicy::RequiredNow);
        }
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn strict_runtime_rejects_unapproved_operation_defined_spec() {
        let mut unapproved = WITHOUT_HYDROGENS_SPEC;
        unapproved.method = "synthetic_unapproved_operation";
        let unapproved = Box::leak(Box::new(unapproved));
        let molecule = crate::Molecule::new();
        let parts = OpParts::new(&molecule, unapproved).unwrap();

        let error = parts
            .finish()
            .expect_err("strict runtime must enforce the OperationDefined allow-list");
        assert!(matches!(
            error,
            OperationError::InvalidInput {
                message: "operation-defined derived effects are currently permitted only for valence in the hydrogen-removal operation family",
                ..
            }
        ));
    }

    #[test]
    fn tautomer_cip_registry_policy_is_absent_or_owned_by_the_single_named_operation() {
        let owners = MOLECULE_OPS
            .iter()
            .copied()
            .filter(|operation| operation.cip_state == CipStatePolicy::TautomerSourceTransition)
            .collect::<Vec<_>>();

        assert!(owners.len() <= 1);
        assert!(
            owners
                .iter()
                .all(|operation| operation.method == "enumerate_tautomers_with_options")
        );
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn tautomer_cip_strict_runtime_accepts_unchanged_and_source_defined_changed_branches() {
        let source = crate::Molecule::from_smiles("C").unwrap();
        let mut parts = MultiMoleculeOpParts::new(&source, &TEST_TAUTOMER_CIP_SPEC).unwrap();

        let unchanged = parts.derive_from_source(|_branch| Ok(())).unwrap();
        parts.emit(unchanged).unwrap();
        let changed = parts
            .derive_from_source(|branch| {
                branch.with_topology_and_properties_mut(|_branch, topology, properties| {
                    topology.atoms[0].set_prop("_CIPCode", "R");
                    properties.set_prop("_StereochemDone", "1");
                    Ok(())
                })?;
                branch.record_topology_edit(TopologyEditKind::Local)?;
                branch.clear_cache(TEST_TAUTOMER_CIP_SPEC.needs_update());
                Ok(())
            })
            .unwrap();
        parts.emit(changed).unwrap();

        let outputs = parts.finish().unwrap();
        assert_eq!(outputs.len(), 2);
        assert!(outputs[0].atoms()[0].prop("_CIPCode").is_none());
        assert_eq!(outputs[1].atoms()[0].prop("_CIPCode"), Some("R"));
        assert_eq!(outputs[1].prop("_StereochemDone"), Some("1"));
        assert!(source.atoms()[0].prop("_CIPCode").is_none());
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn tautomer_cip_strict_runtime_rejects_every_registered_non_tautomer_operation() {
        let molecule = crate::Molecule::from_smiles("C").unwrap();
        for registered in MOLECULE_OPS.iter().copied() {
            if registered.method == "enumerate_tautomers_with_options" {
                continue;
            }
            let mut unapproved = *registered;
            unapproved.cip_state = CipStatePolicy::TautomerSourceTransition;
            let unapproved = Box::leak(Box::new(unapproved));
            let error = OpParts::new(&molecule, unapproved)
                .and_then(OpParts::finish)
                .expect_err("non-tautomer operation must reject the tautomer CIP policy");
            assert!(matches!(
                error,
                OperationError::InvalidInput {
                    message: "source-defined tautomer CIP transitions are owned only by enumerate_tautomers_with_options",
                    ..
                }
            ));
        }
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn tautomer_cip_strict_runtime_requires_both_write_capabilities() {
        let molecule = crate::Molecule::from_smiles("C").unwrap();
        for write in [
            BlockSet::PROPERTIES.union(BlockSet::DERIVED_CACHE),
            BlockSet::TOPOLOGY.union(BlockSet::DERIVED_CACHE),
        ] {
            let mut invalid = TEST_TAUTOMER_CIP_SPEC;
            invalid.access = BlockAccess::new(BlockSet::NONE, write);
            invalid.may_mutate = invalid.access.write();
            let invalid = Box::leak(Box::new(invalid));
            let mut parts = MultiMoleculeOpParts::new(&molecule, invalid).unwrap();
            let error = parts
                .derive_from_source(|_branch| Ok(()))
                .expect_err("missing CIP write capability must fail");
            assert!(matches!(
                error,
                OperationError::InvalidInput {
                    message: "source-defined tautomer CIP transitions require topology and properties write access",
                    ..
                }
            ));
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
        assert_eq!(WITH_2D_COORDINATE_BLOCK_SPEC.domain, OperationDomain::Coordinate);
        assert_eq!(WITH_3D_COORDINATES_SPEC.domain, OperationDomain::Coordinate);
        assert_eq!(WITH_CLEARED_3D_CONFORMERS_SPEC.domain, OperationDomain::Coordinate);
        assert_eq!(WITH_ADDED_3D_CONFORMER_SPEC.domain, OperationDomain::Coordinate);
        assert_eq!(WITH_ONLY_3D_CONFORMER_SPEC.domain, OperationDomain::Coordinate);
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
    fn molalign_operations_have_their_own_feature_boundary() {
        let molalign_operations = MOLECULE_OPS
            .iter()
            .copied()
            .filter(|operation| {
                support_feature_for(operation).map(|feature| feature.name) == Some(crate::MOLALIGN_FEATURE.name)
            })
            .map(|operation| operation.method)
            .collect::<Vec<_>>();
        assert_eq!(
            molalign_operations,
            vec![WITH_ALIGNMENT_TO_SPEC.method, WITH_ALIGNED_CONFORMERS_SPEC.method,]
        );
        for operation in [&WITH_ALIGNMENT_TO_SPEC, &WITH_ALIGNED_CONFORMERS_SPEC] {
            assert_eq!(operation.domain, OperationDomain::Coordinate);
            assert_eq!(operation.parity, ParityPolicy::RequiredNow);
            assert!(support_matrix_contains(operation));
            assert!(invariant_matrix_contains(operation));
            assert!(parity_matrix_contains(operation));
        }
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
        let generated = molecule.with_3d_conformer().expect("default ETKDGv3 conformer");
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
        assert_eq!(result.source_coordinate_dim(), Some(crate::CoordinateDimension::TwoD));
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
        let only = with_two.with_only_3d_conformer(replacement.clone(), true).unwrap();

        assert_eq!(with_two.conformers_3d().len(), 2);
        assert!(cleared.conformers_3d().is_empty());
        assert_eq!(only.conformers_3d().len(), 1);
        assert_eq!(only.conformers_3d()[0].coordinates(), replacement);
        assert_eq!(only.source_coordinate_dim(), Some(crate::CoordinateDimension::ThreeD));

        let mut in_place = with_two.clone();
        in_place.clear_3d_conformers_().unwrap();
        assert!(in_place.conformers_3d().is_empty());
        in_place.set_only_3d_conformer_(second.clone(), true).unwrap();
        assert_eq!(in_place.conformers_3d().len(), 1);
        assert_eq!(in_place.conformers_3d()[0].coordinates(), second);

        let two_d_flagged = with_two.with_only_3d_conformer(replacement.clone(), false).unwrap();
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
        assert_eq!(coords[0], [0.0, 0.0, 1.10]);
    }

    #[test]
    fn add_hs_terminal_coords_default_degree_branch_matches_rdkit_zero_vector() {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let mut coords = vec![[5.0, 5.0, 5.0]];
        for index in 0..5 {
            let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(center, neighbor, crate::BondOrder::Single))
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
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_explicit_hydrogens(2));
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
    fn with_hydrogens_materializes_bracket_hydrogen_count_like_rdkit() {
        let molecule = crate::Molecule::from_smiles("[HH]").unwrap();
        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);

        let result = molecule.with_hydrogens().unwrap();

        assert_eq!(result.num_atoms(), 2);
        assert_eq!(result.num_bonds(), 1);
        assert!(result.atoms().iter().all(|atom| atom.atomic_number() == 1));
        assert!(result.atoms().iter().all(|atom| atom.explicit_hydrogens() == 0));
        assert_eq!(result.bonds()[0].order(), crate::BondOrder::Single);
    }

    #[test]
    fn with_hydrogens_commits_topology_with_rebuilt_adjacency_for_valence_followups() {
        let molecule = crate::Molecule::from_smiles("C=C").unwrap();

        let result = molecule.with_hydrogens().unwrap();
        let assignment = crate::assign_valence(&result, crate::ValenceModel::RdkitLike).unwrap();

        assert_eq!(assignment.explicit_valence, vec![4, 4, 1, 1, 1, 1]);
        assert_eq!(assignment.implicit_hydrogens, vec![0, 0, 0, 0, 0, 0]);
        assert_eq!(result.derived_cache().valence, Some(assignment));
    }

    #[test]
    fn with_hydrogens_immediately_supports_modern_and_legacy_topological_torsion() {
        let molecule = crate::Molecule::from_smiles("CCCC").unwrap();
        let result = molecule.with_hydrogens().unwrap();
        let expected = std::collections::BTreeMap::from([(4_706_550_819_u64, 1)]);

        let modern = crate::topological_torsion_sparse_count_fingerprint(
            &result,
            &crate::TopologicalTorsionFingerprintParams::default(),
        )
        .unwrap();
        assert_eq!(modern.size(), 1_u64 << 36);
        assert_eq!(modern.nonzero_elements(), &expected);

        let legacy =
            crate::topological_torsion_legacy_fingerprint(&result, &crate::TopologicalTorsionLegacyParams::default())
                .unwrap();
        let crate::TopologicalTorsionLegacyResult::SparseCount(legacy) = legacy else {
            panic!("default legacy Topological Torsion must return unfolded counts");
        };
        assert_eq!(legacy.size(), (1_u64 << 36) - 1);
        assert_eq!(legacy.nonzero_elements(), &expected);
    }

    #[test]
    fn with_hydrogens_partial_selection_updates_only_source_on_atoms_cache_entries() {
        let molecule = crate::Molecule::from_smiles("CCO").unwrap();
        let before = molecule.derived_cache().valence.clone().unwrap();
        assert_eq!(before.explicit_valence, vec![1, 2, 1]);
        assert_eq!(before.implicit_hydrogens, vec![3, 2, 1]);

        let result = molecule
            .with_hydrogens_with_params(crate::AddHsParams {
                only_on_atoms: Some(vec![crate::AtomId::new(0)]),
                ..crate::AddHsParams::default()
            })
            .unwrap();
        let after = result.derived_cache().valence.as_ref().unwrap();

        assert_eq!(after.explicit_valence, vec![4, 2, 1, 1, 1, 1]);
        assert_eq!(after.implicit_hydrogens, vec![0, 2, 1, 0, 0, 0]);
    }

    #[test]
    fn with_hydrogens_skip_queries_preserves_skipped_atom_cache_entry() {
        let mut builder = crate::MoleculeBuilder::new();
        let query_carbon = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_query(crate::QueryNode::predicate(crate::AtomQueryPredicate::Any)),
        );
        let plain_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(
                query_carbon,
                plain_carbon,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap().with_assigned_valence_strict(false).unwrap();
        let before = molecule.derived_cache().valence.clone().unwrap();

        let result = molecule
            .with_hydrogens_with_params(crate::AddHsParams {
                skip_queries: true,
                ..crate::AddHsParams::default()
            })
            .unwrap();
        let after = result.derived_cache().valence.as_ref().unwrap();

        assert_eq!(after.explicit_valence[query_carbon.index()], before.explicit_valence[0]);
        assert_eq!(
            after.implicit_hydrogens[query_carbon.index()],
            before.implicit_hydrogens[0]
        );
        assert_eq!(after.explicit_valence[plain_carbon.index()], 4);
        assert_eq!(after.implicit_hydrogens[plain_carbon.index()], 0);
        assert_eq!(after.explicit_valence.len(), result.num_atoms());
        assert_eq!(after.implicit_hydrogens.len(), result.num_atoms());
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
            result.atoms()[1..].iter().map(crate::Atom::isotope).collect::<Vec<_>>(),
            vec![Some(2), Some(3), None]
        );
    }

    #[test]
    fn with_hydrogens_clears_atom_cip_ranks_like_rdkit_addhs() {
        let smiles = "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12";
        let molecule = crate::Molecule::from_smiles(smiles).unwrap();
        assert!(
            molecule.atoms().iter().any(|atom| atom.prop("_CIPRank").is_some()),
            "SMILES sanitize path should assign legacy _CIPRank before AddHs"
        );
        assert_eq!(molecule.prop("_StereochemDone"), Some("1"));

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
            result.atoms().iter().all(|atom| atom.prop("_CIPRank").is_none()),
            "RDKit AddHs clears atom _CIPRank computed props before depiction"
        );
        assert_eq!(
            result.prop("_StereochemDone"),
            None,
            "RDKit AddHs clears the molecule-level stereochemistry computed marker"
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
                pdb_residue_info: Some(crate::AtomPdbResidueInfo::new(" H1 ", 12, "GLY", 3, "A", false)),
            }],
            atoms_to_update_property_cache: vec![nitrogen],
            valence_before_add_hs: Some(
                crate::assign_valence_with_options(&molecule, crate::ValenceModel::RdkitLike, false).unwrap(),
            ),
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
        let result = parts.finish().unwrap();

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
            crate::AtomSpec::new(crate::Element::N)
                .with_pdb_residue_info(crate::AtomPdbResidueInfo::new(" N  ", 10, "GLY", 3, "A", false)),
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
            atoms_to_update_property_cache: vec![hydrogen],
            valence_before_add_hs: Some(
                crate::assign_valence_with_options(&molecule, crate::ValenceModel::RdkitLike, false).unwrap(),
            ),
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
        let result = parts.finish().unwrap();

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
            crate::AtomSpec::new(crate::Element::C)
                .with_pdb_residue_info(crate::AtomPdbResidueInfo::new(" C  ", 7, "GLY", 3, "A", false)),
        );
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        builder
            .add_bond(crate::BondSpec::new(carbon, hydrogen, crate::BondOrder::Single))
            .unwrap();
        builder.set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]]).unwrap();
        let molecule = builder.build().unwrap();
        let original = molecule.clone();

        let result = molecule.without_hydrogens_with_sanitize(false).unwrap();

        assert_eq!(molecule, original);
        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.num_bonds(), 0);
        assert_eq!(result.coordinates_2d(), Some(&[[0.0, 0.0]][..]));
        assert_eq!(result.atoms()[0].pdb_residue_info().unwrap().serial_number(), 7);
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
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_no_implicit(true));
        let protium = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let deuterium = builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        builder
            .add_bond(crate::BondSpec::new(carbon, protium, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(carbon, deuterium, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = crate::RemoveHsParams {
            remove_and_track_isotopes: true,
            ..crate::RemoveHsParams::default()
        };

        let result = molecule.without_hydrogens_with_params(params, false).unwrap();

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
            .add_bond(crate::BondSpec::new(carbon, hydrogen, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(crate::SubstanceGroupId::new(0), crate::SubstanceGroupKind::Superatom)
                    .with_atoms(vec![carbon, hydrogen]),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.without_hydrogens_with_sanitize(false).unwrap();

        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.substance_groups().len(), 1);
        assert_eq!(result.substance_groups()[0].atoms(), &[crate::AtomId::new(0)]);
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
        assert!(result.derived_cache().valence.is_some());
        assert!(result.derived_cache().rings.is_none());
    }

    #[test]
    fn without_hydrogens_without_sanitize_retains_rdkit_pre_removal_valence_cache() {
        let molecule = crate::Molecule::from_smiles("CCO").unwrap().with_hydrogens().unwrap();
        let original = molecule.clone();

        let value = molecule.without_hydrogens_with_sanitize(false).unwrap();
        let mut in_place = molecule.clone();
        in_place.remove_hydrogens_with_sanitize_(false).unwrap();

        assert_eq!(molecule, original);
        for result in [&value, &in_place] {
            assert_eq!(
                result.derived_cache().valence,
                Some(crate::ValenceAssignment {
                    explicit_valence: vec![4, 4, 2],
                    implicit_hydrogens: vec![0, 0, 0],
                })
            );
            assert_eq!(
                crate::assign_valence_with_options(result, crate::ValenceModel::RdkitLike, false,).unwrap(),
                crate::ValenceAssignment {
                    explicit_valence: vec![1, 2, 1],
                    implicit_hydrogens: vec![3, 2, 1],
                },
                "the retained source cache must remain distinguishable from a post-removal recomputation"
            );
        }
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
        let expected = crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, false).unwrap();
        assert_eq!(result.derived_cache().valence, Some(expected));
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
                crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SYMMRINGS | crate::SanitizeOps::SET_AROMATICITY,
            )
            .unwrap();

        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| bond.order() == crate::BondOrder::Aromatic)
        );
        let expected_valence =
            crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, true).unwrap();
        assert_eq!(result.derived_cache().valence, Some(expected_valence));
    }

    #[test]
    fn sanitized_kekulize_runs_kekulize_assignment_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::SYMMRINGS | crate::SanitizeOps::KEKULIZE)
            .unwrap();

        assert!(result.bonds().iter().all(|bond| !bond.is_aromatic()));
        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| matches!(bond.order(), crate::BondOrder::Single | crate::BondOrder::Double))
        );
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
        let molecule = crate::Molecule::from_smiles("Cc1nc2c(nc1C)C(=O)C1=C(C2=O)C2C=CC1CC2").unwrap();
        assert_eq!(
            molecule
                .derived_cache()
                .rings
                .as_ref()
                .expect("default sanitize should initialize SymmSSSR")
                .find_type(),
            crate::RingFindType::SymmSssr
        );

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::KEKULIZE).unwrap();

        let rings = result
            .derived_cache()
            .rings
            .as_ref()
            .expect("Kekulize should initialize ordinary SSSR ring information");
        assert_eq!(rings.find_type(), crate::RingFindType::Sssr);
        assert_eq!(rings.num_rings(), 4);
        assert!(
            result
                .bonds()
                .iter()
                .all(|bond| matches!(bond.order(), crate::BondOrder::Single | crate::BondOrder::Double))
        );

        let symmetrized = molecule
            .sanitize_with_ops(crate::SanitizeOps::SYMMRINGS | crate::SanitizeOps::KEKULIZE)
            .unwrap();
        let rings = symmetrized
            .derived_cache()
            .rings
            .as_ref()
            .expect("explicit SymmSSSR should initialize symmetric ring information");
        assert_eq!(rings.find_type(), crate::RingFindType::SymmSssr);
        assert_eq!(rings.num_rings(), 5);
    }

    #[test]
    fn sanitized_reports_kekulize_failure_step_like_rdkit_operation_that_failed() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("c", false).unwrap();

        let err = molecule.sanitize_with_ops(crate::SanitizeOps::KEKULIZE).unwrap_err();

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

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::NONE).unwrap();
        let expected = crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, false).unwrap();

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
            .add_bond(crate::BondSpec::new(carbon, nitrogen, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(nitrogen, oxygen_single, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(nitrogen, oxygen_double, crate::BondOrder::Double))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::CLEANUP).unwrap();

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

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::CLEANUP).unwrap();

        assert_eq!(result.atoms()[nitrogen_center.index()].formal_charge(), 1);
        assert_eq!(result.atoms()[nitrogen_terminal.index()].formal_charge(), -1);
        assert_eq!(result.bonds()[triple_bond.index()].order(), crate::BondOrder::Double);
    }

    #[test]
    fn sanitized_cleanup_converts_phosphorus_oxo_like_rdkit() {
        let phosphorus_element = crate::Element::from_atomic_number(15).expect("phosphorus atomic number is valid");
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
            .add_bond(crate::BondSpec::new(phosphorus, oxygen, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(phosphorus, nitrogen, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(nitrogen, carbon_double, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::CLEANUP).unwrap();

        assert_eq!(result.atoms()[phosphorus.index()].formal_charge(), 1);
        assert_eq!(result.atoms()[oxygen.index()].formal_charge(), -1);
        assert_eq!(result.atoms()[nitrogen.index()].formal_charge(), 0);
        assert_eq!(result.bonds()[1].order(), crate::BondOrder::Single);
        assert_eq!(result.bonds()[2].order(), crate::BondOrder::Double);
    }

    #[test]
    fn sanitized_phosphorus_cleanup_leaves_double_oxo_without_double_cn_branch_unchanged_like_rdkit() {
        let phosphorus_element = crate::Element::from_atomic_number(15).expect("phosphorus atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let phosphorus = builder.add_atom(crate::AtomSpec::new(phosphorus_element));
        let oxygen_one = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let oxygen_two = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let bond_one = builder
            .add_bond(crate::BondSpec::new(phosphorus, oxygen_one, crate::BondOrder::Double))
            .unwrap();
        let bond_two = builder
            .add_bond(crate::BondSpec::new(phosphorus, oxygen_two, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(phosphorus, carbon, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::CLEANUP).unwrap();

        assert_eq!(result.atoms()[phosphorus.index()].formal_charge(), 0);
        assert_eq!(result.atoms()[oxygen_one.index()].formal_charge(), 0);
        assert_eq!(result.atoms()[oxygen_two.index()].formal_charge(), 0);
        assert_eq!(result.bonds()[bond_one.index()].order(), crate::BondOrder::Double);
        assert_eq!(result.bonds()[bond_two.index()].order(), crate::BondOrder::Double);
    }

    #[test]
    fn sanitized_cleanup_converts_hypervalent_halogen_oxo_like_rdkit() {
        let chlorine = crate::Element::from_atomic_number(17).expect("chlorine atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(chlorine));
        let oxygen_one = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let oxygen_two = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(center, oxygen_one, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(center, oxygen_two, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::CLEANUP).unwrap();

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
        let chlorine = crate::Element::from_atomic_number(17).expect("chlorine atomic number is valid");
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(chlorine));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let double_bond = builder
            .add_bond(crate::BondSpec::new(center, oxygen, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(center, carbon, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::CLEANUP).unwrap();

        assert_eq!(result.atoms()[center.index()].formal_charge(), 0);
        assert_eq!(result.atoms()[oxygen.index()].formal_charge(), 0);
        assert_eq!(result.bonds()[double_bond.index()].order(), crate::BondOrder::Double);
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
        let read = MoleculeReadParts::from_molecule(&molecule);
        let adjacency = sanitize_adjacency(read).unwrap();

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
            .add_bond(crate::BondSpec::new(nitrogen, oxygen, crate::BondOrder::Double))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(nitrogen, carbon, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let read = MoleculeReadParts::from_molecule(&molecule);
        let adjacency = sanitize_adjacency(read).unwrap();
        let mut assignment = SanitizeCleanupAssignment {
            atom_formal_charges: molecule.atoms().iter().map(crate::Atom::formal_charge).collect(),
            bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
        };
        assignment.bond_orders[oxygen_bond.index()] = crate::BondOrder::Single;

        let valence = sanitize_cleanup_explicit_valence(read, &adjacency, &assignment, nitrogen).unwrap();

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
                .add_bond(crate::BondSpec::new(nitrogen, neighbor, crate::BondOrder::Single))
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
                (bond.begin() == nitrogen && bond.end() == metal) || (bond.begin() == metal && bond.end() == nitrogen)
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
                .add_bond(crate::BondSpec::new(donor, neighbor, crate::BondOrder::Single))
                .unwrap();
        }

        let donor_busy = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let busy_c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let busy_c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let busy_c3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for neighbor in [busy_c1, busy_c2, busy_c3] {
            builder
                .add_bond(crate::BondSpec::new(donor_busy, neighbor, crate::BondOrder::Single))
                .unwrap();
        }
        builder
            .add_bond(crate::BondSpec::new(donor_busy, metal_busy, crate::BondOrder::Dative))
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
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_no_implicit(true));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal = builder.add_atom(crate::AtomSpec::new(iron));
        builder
            .add_bond(crate::BondSpec::new(oxygen, carbon, crate::BondOrder::Single))
            .unwrap();
        let metal_bond = builder
            .add_bond(crate::BondSpec::new(oxygen, metal, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS)
            .unwrap();

        assert_eq!(result.bonds()[metal_bond.index()].order(), crate::BondOrder::Single);
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
                .add_bond(crate::BondSpec::new(donor, neighbor, crate::BondOrder::Single))
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
        let read = MoleculeReadParts::from_molecule(&molecule);
        let adjacency = sanitize_adjacency(read).unwrap();
        let valence = read
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
            .unwrap();
        let ranks = read.rank_mol_atoms().unwrap();
        let mut assignment = SanitizeOrganometallicCleanupAssignment {
            bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
            bond_endpoints: molecule.bonds().iter().map(|bond| (bond.begin(), bond.end())).collect(),
        };

        sanitize_metal_bond_cleanup_assignment(read, &adjacency, &valence, &ranks, donor, &mut assignment).unwrap();

        let chosen_metal = assignment
            .bond_endpoints
            .iter()
            .zip(assignment.bond_orders.iter())
            .find_map(|(&(begin, end), &order)| (order == crate::BondOrder::Dative && begin == donor).then_some(end))
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
        let sulfur_atom =
            aromatic_builder.add_atom(crate::AtomSpec::new(carbon).with_aromatic(true).with_no_implicit(true));
        let mut aromatic_neighbors = Vec::new();
        for _ in 0..4 {
            let carbon = aromatic_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            aromatic_builder
                .add_bond(crate::BondSpec::new(sulfur_atom, carbon, crate::BondOrder::Single))
                .unwrap();
            aromatic_neighbors.push(carbon);
        }
        let aromatic = aromatic_builder.build().unwrap();
        let aromatic_read = MoleculeReadParts::from_molecule(&aromatic);
        let aromatic_adj = sanitize_adjacency(aromatic_read).unwrap();
        let aromatic_valence = aromatic_read
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
            .unwrap();

        assert!(
            sanitize_is_hypervalent_nonmetal(aromatic_read, &aromatic_adj, &aromatic_valence, sulfur_atom).unwrap()
        );

        let mut metal_builder = crate::MoleculeBuilder::new();
        let metal_atom = metal_builder.add_atom(crate::AtomSpec::new(iron));
        let ligand = metal_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        metal_builder
            .add_bond(crate::BondSpec::new(metal_atom, ligand, crate::BondOrder::Single))
            .unwrap();
        let metal_molecule = metal_builder.build().unwrap();
        let metal_read = MoleculeReadParts::from_molecule(&metal_molecule);
        let metal_adj = sanitize_adjacency(metal_read).unwrap();
        let metal_valence = metal_read
            .assign_valence_with_options(crate::ValenceModel::RdkitLike, false)
            .unwrap();

        assert!(!sanitize_is_hypervalent_nonmetal(metal_read, &metal_adj, &metal_valence, metal_atom).unwrap());
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
            .add_bond(crate::BondSpec::new(donor, metal_single, crate::BondOrder::Single))
            .unwrap();
        let rewritten_bond = builder
            .add_bond(crate::BondSpec::new(donor, metal_rewritten, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(donor, carbon, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let read = MoleculeReadParts::from_molecule(&molecule);
        let adjacency = sanitize_adjacency(read).unwrap();
        let mut assignment = SanitizeOrganometallicCleanupAssignment {
            bond_orders: molecule.bonds().iter().map(crate::Bond::order).collect(),
            bond_endpoints: molecule.bonds().iter().map(|bond| (bond.begin(), bond.end())).collect(),
        };
        assignment.bond_orders[rewritten_bond.index()] = crate::BondOrder::Dative;

        let metals = sanitize_organometallic_single_bonded_metals(read, &adjacency, &assignment, donor);

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
                crate::BondSpec::new(left, right, crate::BondOrder::Single).with_stereo(crate::BondStereo::AtropCw),
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
                builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_hybridization(crate::Hybridization::Sp2))
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

        assert_eq!(result.bonds()[atrop_bond.index()].stereo(), crate::BondStereo::None);
        assert_eq!(result.bonds()[atrop_bond.index()].stereo_atoms(), None);
        assert!(result.stereo_groups().is_empty());
    }

    #[test]
    fn sanitized_cleanup_chirality_clears_non_sp3_tetrahedral_tag_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_chiral_tag(crate::ChiralTag::TetrahedralCw));
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP_CHIRALITY)
            .unwrap();

        assert_eq!(result.atoms()[0].chiral_tag(), crate::ChiralTag::Unspecified);
    }

    #[test]
    fn sanitized_cleanup_chirality_cleans_stereo_groups_for_non_sp3_tetrahedral_tags_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let atom =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_chiral_tag(crate::ChiralTag::TetrahedralCw));
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

        assert_eq!(result.atoms()[atom.index()].chiral_tag(), crate::ChiralTag::Unspecified);
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

        assert_eq!(result.atoms()[0].chiral_tag(), crate::ChiralTag::Tetrahedral);
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
            .add_bond(crate::BondSpec::new(center, right, crate::BondOrder::Single))
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
            .add_bond(crate::BondSpec::new(center, neighbor, crate::BondOrder::Single))
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

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::SET_CONJUGATION).unwrap();

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
        let center =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_chiral_tag(crate::ChiralTag::TetrahedralCw));
        for _ in 0..4 {
            let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(center, neighbor, crate::BondOrder::Single))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION)
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
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_no_implicit(true));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let metal = builder.add_atom(crate::AtomSpec::new(iron));
        builder
            .add_bond(crate::BondSpec::new(nitrogen, carbon, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(nitrogen, metal, crate::BondOrder::Dative))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION)
            .unwrap();

        assert_eq!(
            result.atoms()[nitrogen.index()].hybridization(),
            crate::Hybridization::Sp2
        );
    }

    #[test]
    fn sanitized_set_hybridization_excludes_zero_bonds_from_num_bonds_plus_lone_pairs() {
        let mut builder = crate::MoleculeBuilder::new();
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_no_implicit(true));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let dummy = builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        builder
            .add_bond(crate::BondSpec::new(oxygen, carbon, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(oxygen, dummy, crate::BondOrder::Zero))
            .unwrap();
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION)
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
                .add_bond(crate::BondSpec::new(center, neighbor, crate::BondOrder::Single))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::PROPERTIES | crate::SanitizeOps::SET_HYBRIDIZATION)
            .unwrap();

        assert_eq!(result.atoms()[center.index()].hybridization(), crate::Hybridization::Sp);
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

        let result = molecule.sanitize_with_ops(crate::SanitizeOps::FIND_RADICALS).unwrap();

        assert_eq!(result.atoms()[0].radical_electrons(), 1);
        let expected = crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, false).unwrap();
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
            result.derived_cache().valence.as_ref().unwrap().implicit_hydrogens[0],
            0
        );
    }

    #[test]
    fn sanitize_adjust_hydrogens_assignment_preserves_existing_explicit_hydrogen_when_delta_is_zero() {
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
    fn sanitize_adjust_hydrogens_assignment_leaves_stable_explicit_hydrogens_unchanged_like_rdkit() {
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
    fn sanitize_adjust_hydrogens_repopulates_property_cache_without_properties_flag_like_rdkit() {
        let molecule = crate::Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

        let result = molecule
            .sanitize_with_ops(crate::SanitizeOps::ADJUST_HYDROGENS)
            .unwrap();
        let expected = crate::assign_valence_with_options(&result, crate::ValenceModel::RdkitLike, false).unwrap();

        assert_eq!(result.derived_cache().valence, Some(expected));
        assert!(!crate::valence::needs_update_property_cache(
            crate::read_parts::MoleculeReadParts::from_molecule(&result)
        ));
    }

    #[test]
    fn supported_kekulize_runs_through_operation_pipeline_without_changing_source() {
        assert_eq!(
            WITH_KEKULIZED_BONDS_SPEC.support,
            SupportStatus::SupportedWithRdkitParity {
                rdkit_version: "2026.03.1"
            }
        );

        let molecule = crate::Molecule::new();
        let original = molecule.clone();
        let result = molecule
            .with_kekulized_bonds(true)
            .expect("supported kekulize should satisfy its operation contract");

        assert_eq!(molecule, original);
        assert_eq!(result.atoms(), original.atoms());
        assert_eq!(result.bonds(), original.bonds());
        assert_eq!(result.coordinates_2d(), original.coordinates_2d());
        assert_eq!(result.conformers_3d(), original.conformers_3d());
        assert_eq!(result.source_coordinate_dim(), original.source_coordinate_dim());
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

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "operation read capability does not permit access to the coordinates block")]
    fn topology_read_capability_does_not_expose_coordinates() {
        let molecule = crate::Molecule::new();
        let parts = OpParts::new(&molecule, &ASSIGNED_VALENCE_SPEC).unwrap();
        let read = parts.begin_topology_read().unwrap();

        let _ = read.coordinates();
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
            .finish()
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
            DerivedState::RING_FAMILIES | DerivedState::STEREO | DerivedState::DRAWING | DerivedState::FINGERPRINT,
        );
        let result = parts
            .finish()
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
    fn finish_accepts_untouched_operation_without_derived_state_updates() {
        let molecule = crate::Molecule::new();
        let parts = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC).unwrap();

        let result = parts
            .finish()
            .expect("an untouched operation has no derived-state effects to fulfill");
        assert_eq!(result, molecule);
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    fn committed_write_cannot_bypass_cache_obligations() {
        let molecule = crate::Molecule::new();
        let mut parts = OpParts::new(&molecule, &WITH_2D_COORDINATES_SPEC).unwrap();
        parts.with_coordinates_mut(|_parts, _coordinates| Ok(())).unwrap();

        let error = parts
            .finish()
            .expect_err("a committed write must satisfy derived-state obligations");

        assert!(matches!(
            error,
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
            .finish()
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
        parts.record_topology_edit(TopologyEditKind::Appending).unwrap();
        parts.record_topology_mapping(TopologyMapping::with_appended(0, 0, 0, 0));
        parts.clear_cache(WITH_HYDROGENS_SPEC.needs_update());

        let err = parts
            .finish()
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
        parts.record_topology_edit(TopologyEditKind::Appending).unwrap();
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
            .finish()
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
            .finish()
            .expect("without op-contracts, cache permission checks are disabled");

        let missing_update = OpParts::new(&molecule, &TEST_NEEDS_VALENCE_UPDATE_SPEC).unwrap();
        missing_update
            .finish()
            .expect("without op-contracts, needs_update checks are disabled");
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "cache write permission violation")]
    fn set_valence_cache_panics_without_recompute() {
        // WITH_2D_COORDINATES_SPEC has no recompute permission, so
        // set_valence_cache must panic.
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
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("COW source molecule must be valid");
        let mut parts = OpParts::new(&molecule, &WITH_KEKULIZED_BONDS_SPEC).unwrap();

        let mut topology = parts.begin_topology_mut().unwrap();
        topology.atoms[0].set_prop("cow", "changed");
        let valence = parts
            .with_topology_read_parts(topology.clone(), |read| {
                read.assign_valence_with_options(crate::ValenceModel::RdkitLike, true)
                    .map_err(|source| OperationError::Valence {
                        operation: &WITH_KEKULIZED_BONDS_SPEC,
                        source,
                    })
            })
            .unwrap();
        parts.commit_topology(topology).unwrap();
        parts.record_topology_edit(TopologyEditKind::Local).unwrap();
        parts.set_rings_cache(crate::RingInfo::new(crate::RingFindType::SymmSssr, 1, 0));
        parts.set_valence_cache(valence);
        parts.clear_cache(DerivedState::AROMATICITY);
        parts.clear_cache(DerivedState::DRAWING | DerivedState::FINGERPRINT);

        let result = parts
            .finish()
            .expect("COW topology edit should produce a valid molecule");

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.atomic_numbers(), vec![6]);
        assert_eq!(molecule.atoms()[0].prop("cow"), None);
        assert_eq!(result.atoms()[0].prop("cow"), Some("changed"));
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
                vec![Some("c0".to_string()), Some("o1".to_string()), Some("n2".to_string())],
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
        coordinates.remap_topology(&mapping.retained_atom_indices());
        properties.remap_topology(&mapping.atoms.new_to_old, &mapping.bonds.new_to_old);
        parts.record_topology_edit(TopologyEditKind::Compacting).unwrap();
        parts.record_topology_mapping(mapping.clone());
        assert_eq!(
            mapping.atoms().old_to_new(),
            &[Some(crate::AtomId::new(0)), None, Some(crate::AtomId::new(1))]
        );
        assert_eq!(mapping.bonds().old_to_new(), &[None, None]);
        parts.clear_cache(WITHOUT_HYDROGENS_SPEC.needs_update());
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();

        let result = parts
            .finish()
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
                crate::SubstanceGroup::new(crate::SubstanceGroupId::new(0), crate::SubstanceGroupKind::Superatom)
                    .with_atoms(vec![c0]),
            )
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(crate::SubstanceGroupId::new(1), crate::SubstanceGroupKind::Data)
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
        coordinates.remap_topology(&mapping.retained_atom_indices());
        properties.remap_topology(&mapping.atoms.new_to_old, &mapping.bonds.new_to_old);
        parts.record_topology_edit(TopologyEditKind::Compacting).unwrap();
        parts.record_topology_mapping(mapping);
        parts.clear_cache(WITHOUT_HYDROGENS_SPEC.needs_update());
        parts.commit_topology(topology).unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        parts.commit_properties(properties).unwrap();

        let result = parts
            .finish()
            .expect("strong compacting edit should preserve surviving SGroup parent links");

        assert_eq!(result.substance_groups().len(), 2);
        assert_eq!(result.substance_groups()[0].atoms(), &[crate::AtomId::new(0)]);
        assert_eq!(result.substance_groups()[1].atoms(), &[crate::AtomId::new(1)]);
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
        assert_eq!(result.properties().sdf_property_lists()[1].values(), &[None]);
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
        let bond_orders = kekulized.bonds().iter().map(|bond| bond.order()).collect::<Vec<_>>();

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
