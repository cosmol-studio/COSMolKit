//! Stable operation metadata values generated from `molecule_ops!`.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SupportStatus {
    Supported,
    SupportedWithRdkitParity,
    PreservedOnly,
    Experimental,
    Unsupported { reason: &'static str },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FeatureSpec {
    pub name: &'static str,
    pub category: &'static str,
    pub status: SupportStatus,
    pub rdkit_parity_sensitive: bool,
    pub docs: &'static str,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct UnsupportedFeatureError {
    pub feature: &'static str,
    pub reason: &'static str,
}

impl UnsupportedFeatureError {
    #[must_use]
    pub const fn from_spec(spec: &'static FeatureSpec) -> Self {
        let reason = match spec.status {
            SupportStatus::Unsupported { reason } => reason,
            _ => "feature is not marked unsupported",
        };
        Self {
            feature: spec.name,
            reason,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoleculeOpOutput {
    Single,
    Multiple,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OperationDomain {
    Topology,
    Coordinate,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MoleculeOpKind {
    Weak,
    Strong,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TopologyEditKind {
    None,
    Local,
    Compacting,
    Appending,
    Renumbering,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
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

    #[must_use]
    pub const fn intersection(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }

    #[must_use]
    pub const fn difference(self, other: Self) -> Self {
        Self(self.0 & !other.0)
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    #[must_use]
    pub const fn bits(self) -> u8 {
        self.0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BlockAccess {
    read: BlockSet,
    write: BlockSet,
}

impl BlockAccess {
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

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct DerivedState(u16);

impl DerivedState {
    pub const NONE: Self = Self(0);
    pub const RINGS: Self = Self(1 << 0);
    pub const RING_FAMILIES: Self = Self(1 << 1);
    pub const VALENCE: Self = Self(1 << 2);
    pub const AROMATICITY: Self = Self(1 << 3);
    pub const STEREO: Self = Self(1 << 4);
    pub const COORDINATES: Self = Self(1 << 5);
    pub const DRAWING: Self = Self(1 << 6);
    pub const FINGERPRINT: Self = Self(1 << 7);

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

    #[must_use]
    pub const fn intersection(self, other: Self) -> Self {
        Self(self.0 & other.0)
    }

    #[must_use]
    pub const fn difference(self, other: Self) -> Self {
        Self(self.0 & !other.0)
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    #[must_use]
    pub const fn bits(self) -> u16 {
        self.0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DerivedEffects {
    pub recompute: DerivedState,
    pub preserve: DerivedState,
    pub invalidate: DerivedState,
    pub operation_defined: DerivedState,
}

impl DerivedEffects {
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
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CipStatePolicy {
    Preserve,
    ClearComputed,
    Assign,
    TautomerSourceTransition,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
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

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    #[must_use]
    pub const fn bits(self) -> u8 {
        self.0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MappingRequirement {
    None,
    Identity,
    Required,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParityPolicy {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SupportMatrixEntry {
    pub feature: &'static FeatureSpec,
    pub operation: Option<&'static MoleculeOpSpec>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OperationInvariantEntry {
    pub operation: &'static MoleculeOpSpec,
    pub profile: &'static str,
}

impl OperationInvariantEntry {
    #[must_use]
    pub const fn for_operation(operation: &'static MoleculeOpSpec, profile: &'static str) -> Self {
        Self { operation, profile }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ParityMatrixEntry {
    pub operation: &'static MoleculeOpSpec,
    pub feature: &'static FeatureSpec,
    pub profile: &'static str,
    pub rdkit_version: Option<&'static str>,
}

/// Allocation-free, declaration-ordered projection of unique feature specs.
#[derive(Clone, Debug)]
pub struct FeatureSpecIter {
    entries: &'static [SupportMatrixEntry],
    index: usize,
}

impl Iterator for FeatureSpecIter {
    type Item = &'static FeatureSpec;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(entry) = self.entries.get(self.index) {
            let current_index = self.index;
            self.index += 1;

            let already_seen = self.entries[..current_index]
                .iter()
                .any(|previous| core::ptr::eq(previous.feature, entry.feature));
            if !already_seen {
                return Some(entry.feature);
            }
        }
        None
    }
}

/// Returns each feature referenced by the generated support matrix once.
#[must_use]
pub fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter {
        entries: super::registry::SUPPORT_MATRIX,
        index: 0,
    }
}

/// Looks up a generated feature by its exact canonical name.
#[must_use]
pub fn feature_spec(name: &str) -> Option<&'static FeatureSpec> {
    feature_specs().find(|feature| feature.name == name)
}

/// Returns the generated operation declarations in declaration order.
#[must_use]
pub fn operation_specs() -> &'static [&'static MoleculeOpSpec] {
    super::registry::MOLECULE_OPS
}

/// Looks up a generated operation by its exact canonical method name.
#[must_use]
pub fn operation_spec(method: &str) -> Option<&'static MoleculeOpSpec> {
    operation_specs()
        .iter()
        .copied()
        .find(|operation| operation.method == method)
}

/// Returns the generated support matrix without copying its rows.
#[must_use]
pub fn support_matrix() -> &'static [SupportMatrixEntry] {
    super::registry::SUPPORT_MATRIX
}

/// Returns the generated operation-invariant matrix without copying its rows.
#[must_use]
pub fn operation_invariant_matrix() -> &'static [OperationInvariantEntry] {
    super::registry::OPERATION_INVARIANT_MATRIX
}

/// Returns the generated parity matrix without copying its rows.
#[must_use]
pub fn parity_matrix() -> &'static [ParityMatrixEntry] {
    super::registry::PARITY_MATRIX
}

/// Looks up an invariant row by the exact generated operation method.
#[must_use]
pub fn operation_invariant(method: &str) -> Option<&'static OperationInvariantEntry> {
    operation_invariant_matrix()
        .iter()
        .find(|entry| entry.operation.method == method)
}

/// Looks up a parity row by the exact generated operation method.
#[must_use]
pub fn operation_parity(method: &str) -> Option<&'static ParityMatrixEntry> {
    parity_matrix()
        .iter()
        .find(|entry| entry.operation.method == method)
}
