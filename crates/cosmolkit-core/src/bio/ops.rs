//! Operation contract system for BioStructure.
//!
//! Mirrors the `MoleculeOpSpec` / `OpParts` pattern in `ops.rs`, adapted for
//! the BioStructure hierarchy. Operation bodies must be registered via the
//! `bio_structure_ops!` macro. `BioOpParts` is wrapper-owned migration and
//! contract-recording machinery; operation-specific behavior belongs in the
//! operation body or in bio domain modules, not in `BioOpParts`.

use std::marker::PhantomData;

use crate::{
    SupportStatus,
    bio::{AssemblyId, AtomId, BioStructure, BondId, ChainId, EntityId, ModelId, ResidueId},
};

// ---------------------------------------------------------------------------
// BioBlockSet — bitmask of mutable blocks
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioBlockSet(u32);

impl BioBlockSet {
    pub const NONE: Self = Self(0);
    pub const ATOMS: Self = Self(1 << 0);
    pub const RESIDUES: Self = Self(1 << 1);
    pub const CHAINS: Self = Self(1 << 2);
    pub const ENTITIES: Self = Self(1 << 3);
    pub const MODELS: Self = Self(1 << 4);
    pub const COORDINATES: Self = Self(1 << 5);
    pub const BONDS: Self = Self(1 << 6);
    pub const ASSEMBLIES: Self = Self(1 << 7);
    pub const ANNOTATIONS: Self = Self(1 << 8);
    pub const DERIVED_CACHE: Self = Self(1 << 9);
    pub const PROPERTIES: Self = Self(1 << 10);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

// ---------------------------------------------------------------------------
// BioStateSet — bitmask of structural state that must be handled
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioStateSet(u32);

impl BioStateSet {
    pub const NONE: Self = Self(0);
    pub const HIERARCHY: Self = Self(1 << 0);
    pub const RESIDUE_SPANS: Self = Self(1 << 1);
    pub const CHAIN_SPANS: Self = Self(1 << 2);
    pub const MODEL_SPANS: Self = Self(1 << 3);
    pub const COORDINATE_ALIGNMENT: Self = Self(1 << 4);
    pub const ENTITY_MAPPING: Self = Self(1 << 5);
    pub const ALTLOC_GROUPS: Self = Self(1 << 6);
    pub const ASSEMBLY_REFERENCES: Self = Self(1 << 7);
    pub const BOND_REFERENCES: Self = Self(1 << 8);
    pub const SELECTION_PROVENANCE: Self = Self(1 << 9);
    pub const POLYMER_ANNOTATION: Self = Self(1 << 10);
    pub const SECONDARY_STRUCTURE: Self = Self(1 << 11);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

// ---------------------------------------------------------------------------
// BioDerivedState — bitmask of derived cache entries
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioDerivedState(u64);

impl BioDerivedState {
    pub const NONE: Self = Self(0);
    pub const ATOM_INDEX: Self = Self(1 << 0);
    pub const RESIDUE_INDEX: Self = Self(1 << 1);
    pub const CHAIN_INDEX: Self = Self(1 << 2);
    pub const ENTITY_INDEX: Self = Self(1 << 3);
    pub const SEQUENCE_CACHE: Self = Self(1 << 4);
    pub const POLYMER_CACHE: Self = Self(1 << 5);
    pub const ALTLOC_CACHE: Self = Self(1 << 6);
    pub const ASSEMBLY_CACHE: Self = Self(1 << 7);
    pub const BOND_CACHE: Self = Self(1 << 8);
    pub const BACKBONE_GEOMETRY: Self = Self(1 << 9);
    pub const SIDECHAIN_GEOMETRY: Self = Self(1 << 10);
    pub const NUCLEIC_GEOMETRY: Self = Self(1 << 11);
    pub const SECONDARY_STRUCTURE: Self = Self(1 << 12);
    pub const CONTACT_MAP: Self = Self(1 << 13);
    pub const GRAPH_CACHE: Self = Self(1 << 14);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

impl std::ops::BitOr for BioDerivedState {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

// ---------------------------------------------------------------------------
// Operation classification enums
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioOpKind {
    /// Does not change row identity (coordinate transforms, annotations, cache).
    Weak,
    /// Changes row topology or hierarchy identity (remove residues, assembly
    /// expansion, altloc resolution, merge). Requires a mapping.
    Strong,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioEditKind {
    None,
    Local,
    Compacting,
    Expanding,
    Renumbering,
    Splitting,
    Merging,
    Transforming,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioOpDomain {
    Selection,
    Hierarchy,
    Coordinate,
    Assembly,
    Annotation,
    Bonding,
    Polymer,
    ChemistryBridge,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioParityPolicy {
    NotApplicable,
    GemmiWhenApplicable,
    BiopythonWhenApplicable,
    PdbSpecRequired,
    RequiredNow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MappingRequirement {
    None,
    Identity,
    Required,
}

// ---------------------------------------------------------------------------
// BioStructureOpSpec
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioStructureOpSpec {
    pub method: &'static str,
    pub impl_fn: &'static str,
    pub domain: BioOpDomain,
    pub kind: BioOpKind,
    pub edit_kind: BioEditKind,
    pub may_mutate: BioBlockSet,
    pub auto_remap: BioBlockSet,
    pub must_handle: BioStateSet,
    pub needs_update: BioDerivedState,
    pub requires_mapping: MappingRequirement,
    pub allows_noop: bool,
    pub support: SupportStatus,
    pub parity: BioParityPolicy,
    pub io_roundtrip: bool,
}

impl std::fmt::Display for BioStructureOpSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.method)
    }
}

// ---------------------------------------------------------------------------
// Operation outcome and errors
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BioOpOutcome {
    Changed,
    NoOp { reason: &'static str },
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BioOperationError {
    #[error("{operation}: unsupported operation: {reason}")]
    Unsupported {
        operation: &'static BioStructureOpSpec,
        reason: &'static str,
    },
    #[error("{operation}: invalid input: {message}")]
    InvalidInput {
        operation: &'static BioStructureOpSpec,
        message: &'static str,
    },
    #[error("{operation}: invariant violation: {message}")]
    InvariantViolation {
        operation: &'static BioStructureOpSpec,
        message: &'static str,
    },
}

// ---------------------------------------------------------------------------
// Mapping system
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioRowMapping<T: Copy> {
    pub old_to_new: Vec<Option<T>>,
    pub new_to_old: Vec<T>,
}

impl<T: Copy> BioRowMapping<T> {
    #[must_use]
    pub fn identity(len: usize, make: impl Fn(u32) -> T) -> Self {
        let ids: Vec<T> = (0..len as u32).map(&make).collect();
        Self {
            old_to_new: ids.iter().copied().map(Some).collect(),
            new_to_old: ids,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioStructureMapping {
    pub atoms: BioRowMapping<AtomId>,
    pub residues: BioRowMapping<ResidueId>,
    pub chains: BioRowMapping<ChainId>,
    pub entities: BioRowMapping<EntityId>,
    pub models: BioRowMapping<ModelId>,
    pub bonds: BioRowMapping<BondId>,
    pub assemblies: BioRowMapping<AssemblyId>,
}

// ---------------------------------------------------------------------------
// BioOpParts — the only mutable capability object for operation bodies
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioOperationTrace {
    touched_blocks: BioBlockSet,
    remapped_blocks: BioBlockSet,
    handled: BioStateSet,
    cleared_cache: BioDerivedState,
    updated_cache: BioDerivedState,
    outcome_recorded: bool,
}

pub struct BioOpParts<'a> {
    spec: &'static BioStructureOpSpec,
    working: BioStructure,
    mapping: Option<BioStructureMapping>,
    _source: PhantomData<&'a BioStructure>,

    #[cfg(feature = "op-contracts")]
    trace: BioOperationTrace,
}

impl<'a> BioOpParts<'a> {
    pub(crate) fn new(source: &'a BioStructure, spec: &'static BioStructureOpSpec) -> Self {
        Self {
            spec,
            working: source.clone(),
            mapping: None,
            _source: PhantomData,
            #[cfg(feature = "op-contracts")]
            trace: BioOperationTrace {
                touched_blocks: BioBlockSet::NONE,
                remapped_blocks: BioBlockSet::NONE,
                handled: BioStateSet::NONE,
                cleared_cache: BioDerivedState::NONE,
                updated_cache: BioDerivedState::NONE,
                outcome_recorded: false,
            },
        }
    }

    #[must_use]
    pub(crate) fn structure(&self) -> &BioStructure {
        &self.working
    }

    pub(crate) fn clear_cache(&mut self, states: BioDerivedState) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.cleared_cache = self.trace.cleared_cache | states;
        }
        let _ = states;
    }

    fn mark_handled(&mut self, states: BioStateSet) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.handled = self.trace.handled.union(states);
        }
        let _ = states;
    }

    fn record_remapped(&mut self, blocks: BioBlockSet) {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.remapped_blocks = self.trace.remapped_blocks.union(blocks);
        }
        let _ = blocks;
    }

    pub(crate) fn record_identity_mapping(&mut self) {
        self.mapping = Some(BioStructureMapping {
            atoms: BioRowMapping::identity(self.working.atoms.len(), AtomId::new),
            residues: BioRowMapping::identity(self.working.residues.len(), ResidueId::new),
            chains: BioRowMapping::identity(self.working.chains.len(), ChainId::new),
            entities: BioRowMapping {
                old_to_new: vec![],
                new_to_old: vec![],
            },
            models: BioRowMapping::identity(self.working.models.len(), ModelId::new),
            bonds: BioRowMapping {
                old_to_new: vec![],
                new_to_old: vec![],
            },
            assemblies: BioRowMapping {
                old_to_new: vec![],
                new_to_old: vec![],
            },
        });
    }

    pub(crate) fn mark_hierarchy_contract_handled(&mut self) {
        self.mark_handled(
            BioStateSet::HIERARCHY
                .union(BioStateSet::RESIDUE_SPANS)
                .union(BioStateSet::CHAIN_SPANS)
                .union(BioStateSet::MODEL_SPANS)
                .union(BioStateSet::COORDINATE_ALIGNMENT),
        );
        self.record_remapped(BioBlockSet::COORDINATES);
    }

    pub(crate) fn remove_residues(
        &mut self,
        residues_to_remove: &[ResidueId],
    ) -> Result<&BioStructureMapping, BioOperationError> {
        self.assert_compacting_hierarchy_edit_allowed()?;

        let mut remove_residue = vec![false; self.working.residues.len()];
        for residue in residues_to_remove {
            if let Some(slot) = remove_residue.get_mut(residue.index() as usize) {
                *slot = true;
            }
        }

        let keep_residue: Vec<bool> = remove_residue.iter().map(|remove| !remove).collect();
        let keep_atom: Vec<bool> = self
            .working
            .atoms
            .iter()
            .map(|atom| keep_residue[atom.residue_id.index() as usize])
            .collect();

        let mut atom_old_to_new = vec![None; keep_atom.len()];
        let mut atom_new_to_old = Vec::new();
        for (old, keep) in keep_atom.iter().copied().enumerate() {
            if keep {
                let new_id = AtomId::new(atom_new_to_old.len() as u32);
                atom_old_to_new[old] = Some(new_id);
                atom_new_to_old.push(AtomId::new(old as u32));
            }
        }

        let mut residue_old_to_new = vec![None; keep_residue.len()];
        let mut residue_new_to_old = Vec::new();
        for (old, keep) in keep_residue.iter().copied().enumerate() {
            if keep {
                let new_id = ResidueId::new(residue_new_to_old.len() as u32);
                residue_old_to_new[old] = Some(new_id);
                residue_new_to_old.push(ResidueId::new(old as u32));
            }
        }

        let new_atoms: Vec<_> = atom_new_to_old
            .iter()
            .map(|old_id| {
                let mut row = self.working.atoms[old_id.index() as usize].clone();
                row.residue_id = residue_old_to_new[row.residue_id.index() as usize]
                    .expect("kept atom must belong to a kept residue");
                row
            })
            .collect();

        let new_residues: Vec<_> = residue_new_to_old
            .iter()
            .map(|old_id| {
                let old_row = &self.working.residues[old_id.index() as usize];
                let new_start = (old_row.atom_span.start..old_row.atom_span.end())
                    .find_map(|idx| atom_old_to_new[idx as usize].map(AtomId::index))
                    .unwrap_or(new_atoms.len() as u32);
                let new_len = (old_row.atom_span.start..old_row.atom_span.end())
                    .filter(|idx| atom_old_to_new[*idx as usize].is_some())
                    .count() as u32;
                let mut row = old_row.clone();
                row.atom_span = crate::bio::RowSpan::new(new_start, new_len);
                row
            })
            .collect();

        let new_chains: Vec<_> = self
            .working
            .chains
            .iter()
            .map(|chain| {
                let new_start = (chain.residue_span.start..chain.residue_span.end())
                    .find_map(|idx| residue_old_to_new[idx as usize].map(ResidueId::index))
                    .unwrap_or(new_residues.len() as u32);
                let new_len = (chain.residue_span.start..chain.residue_span.end())
                    .filter(|idx| residue_old_to_new[*idx as usize].is_some())
                    .count() as u32;
                let mut row = chain.clone();
                row.residue_span = crate::bio::RowSpan::new(new_start, new_len);
                row
            })
            .collect();

        let new_positions: Vec<_> = atom_new_to_old
            .iter()
            .map(|old_id| self.working.coordinates.positions[old_id.index() as usize])
            .collect();

        self.record_mutation(BioBlockSet::ATOMS);
        self.working.atoms = new_atoms;
        self.record_mutation(BioBlockSet::RESIDUES);
        self.working.residues = new_residues;
        self.record_mutation(BioBlockSet::CHAINS);
        self.working.chains = new_chains;
        self.record_mutation(BioBlockSet::COORDINATES);
        self.working.coordinates.positions = new_positions;

        self.mark_hierarchy_contract_handled();
        self.clear_cache(self.spec.needs_update);
        self.mapping = Some(BioStructureMapping {
            atoms: BioRowMapping {
                old_to_new: atom_old_to_new,
                new_to_old: atom_new_to_old,
            },
            residues: BioRowMapping {
                old_to_new: residue_old_to_new,
                new_to_old: residue_new_to_old,
            },
            chains: BioRowMapping::identity(self.working.chains.len(), ChainId::new),
            entities: BioRowMapping {
                old_to_new: vec![],
                new_to_old: vec![],
            },
            models: BioRowMapping::identity(self.working.models.len(), ModelId::new),
            bonds: BioRowMapping {
                old_to_new: vec![],
                new_to_old: vec![],
            },
            assemblies: BioRowMapping {
                old_to_new: vec![],
                new_to_old: vec![],
            },
        });

        Ok(self.mapping.as_ref().expect("mapping was just recorded"))
    }

    pub(crate) fn finish(
        #[cfg_attr(not(feature = "op-contracts"), allow(unused_mut))] mut self,
        outcome: BioOpOutcome,
    ) -> Result<BioStructure, BioOperationError> {
        #[cfg(feature = "op-contracts")]
        {
            self.trace.outcome_recorded = true;
            let _ = outcome;
            self.validate_contract()?;
        }
        #[cfg(not(feature = "op-contracts"))]
        {
            let _ = outcome;
        }
        if self.spec.requires_mapping == MappingRequirement::Required && self.mapping.is_none() {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "strong operation did not record a BioStructureMapping",
            });
        }
        crate::bio_invariants::enforce_bio_structure_invariants(&self.working).map_err(
            |message| BioOperationError::InvariantViolation {
                operation: self.spec,
                message,
            },
        )?;
        Ok(self.working)
    }

    fn record_mutation(&mut self, block: BioBlockSet) {
        #[cfg(feature = "op-contracts")]
        {
            assert!(
                self.spec.may_mutate.contains(block),
                "bio operation `{}` attempted to mutate a block outside its registry permissions",
                self.spec.method
            );
            self.trace.touched_blocks = self.trace.touched_blocks.union(block);
        }
        let _ = block;
    }

    fn assert_compacting_hierarchy_edit_allowed(&self) -> Result<(), BioOperationError> {
        if self.spec.kind != BioOpKind::Strong {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "compacting hierarchy edits require a strong operation",
            });
        }
        if self.spec.edit_kind != BioEditKind::Compacting {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "operation registry does not allow compacting hierarchy edits",
            });
        }
        if self.spec.requires_mapping != MappingRequirement::Required {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "compacting hierarchy edits must require a mapping",
            });
        }
        Ok(())
    }

    #[cfg(feature = "op-contracts")]
    fn validate_contract(&self) -> Result<(), BioOperationError> {
        if !self.trace.handled.contains(self.spec.must_handle) {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not handle every required BioStructure state",
            });
        }
        let updated_or_cleared = self.trace.cleared_cache | self.trace.updated_cache;
        if !updated_or_cleared.contains(self.spec.needs_update) {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "operation body did not clear or update every required BioStructure cache state",
            });
        }
        if !self.trace.remapped_blocks.contains(self.spec.auto_remap) {
            return Err(BioOperationError::InvalidInput {
                operation: self.spec,
                message: "operation did not remap every registry-required BioStructure block",
            });
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Registry tables (populated by bio_structure_ops! macro)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioSupportMatrixEntry {
    pub feature: &'static crate::FeatureSpec,
    pub operation: &'static BioStructureOpSpec,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioOperationInvariantEntry {
    pub operation: &'static BioStructureOpSpec,
    pub profile: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioParityMatrixEntry {
    pub operation: &'static BioStructureOpSpec,
    pub profile: &'static str,
}

// ---------------------------------------------------------------------------
// Operation registry
// ---------------------------------------------------------------------------

use cosmolkit_macros::bio_op_body;
use cosmolkit_macros::bio_structure_ops;

bio_structure_ops! {
    op remove_waters() {
        method: without_waters,
        impl_fn: remove_waters_impl,
        domain: selection,
        kind: strong,
        edit_kind: compacting,
        may_mutate: [atoms, residues, chains, models, coordinates],
        auto_remap: [coordinates],
        must_handle: [hierarchy, residue_spans, chain_spans, model_spans, coordinate_alignment],
        needs_update: [atom_index, residue_index, chain_index],
        requires_mapping: required,
        allows_noop: true,
        feature: BIO_SELECTION_FEATURE,
        parity: not_applicable,
        io_roundtrip: false,
        invariant_profile: "strong_bio_hierarchy",
    }
}

#[bio_op_body(remove_waters, parts)]
fn remove_waters_impl() -> Result<BioOpOutcome, BioOperationError> {
    use crate::bio::ResidueKind;

    let water_residue_ids: Vec<ResidueId> = parts
        .structure()
        .residues
        .iter()
        .enumerate()
        .filter(|(_, r)| r.kind == ResidueKind::Water)
        .map(|(index, _)| ResidueId::new(index as u32))
        .collect();

    if water_residue_ids.is_empty() {
        parts.record_identity_mapping();
        parts.mark_hierarchy_contract_handled();
        parts.clear_cache(BIO_REMOVE_WATERS_SPEC.needs_update);
        return Ok(BioOpOutcome::NoOp {
            reason: "no water residues found",
        });
    }

    parts.remove_residues(&water_residue_ids)?;
    Ok(BioOpOutcome::Changed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio::*;

    const TEST_WEAK_COMPACT_SPEC: BioStructureOpSpec = BioStructureOpSpec {
        method: "test_weak_compact",
        impl_fn: "test_weak_compact_impl",
        domain: BioOpDomain::Selection,
        kind: BioOpKind::Weak,
        edit_kind: BioEditKind::Compacting,
        may_mutate: BioBlockSet::ATOMS,
        auto_remap: BioBlockSet::NONE,
        must_handle: BioStateSet::NONE,
        needs_update: BioDerivedState::NONE,
        requires_mapping: MappingRequirement::Required,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: BioParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    #[cfg(feature = "op-contracts")]
    const TEST_UNAUTHORIZED_REMOVE_RESIDUES_SPEC: BioStructureOpSpec = BioStructureOpSpec {
        method: "test_unauthorized_remove_residues",
        impl_fn: "test_unauthorized_remove_residues_impl",
        domain: BioOpDomain::Selection,
        kind: BioOpKind::Strong,
        edit_kind: BioEditKind::Compacting,
        may_mutate: BioBlockSet::NONE,
        auto_remap: BioBlockSet::NONE,
        must_handle: BioStateSet::NONE,
        needs_update: BioDerivedState::NONE,
        requires_mapping: MappingRequirement::Required,
        allows_noop: true,
        support: SupportStatus::Experimental,
        parity: BioParityPolicy::NotApplicable,
        io_roundtrip: false,
    };

    fn make_structure_with_waters() -> BioStructure {
        // Model 0 → Chain 0 → Residue 0 (ALA, atom 0), Residue 1 (HOH, atom 1)
        let mut s = BioStructure::new();
        s.models.push(ModelRow {
            chain_span: RowSpan::new(0, 1),
            source_model_number: Some(1),
        });
        s.chains.push(ChainRow {
            model_id: ModelId::new(0),
            entity_id: None,
            residue_span: RowSpan::new(0, 2),
            kind: ChainKind::Mixed,
            source: ChainSourceIds {
                auth_chain_id: None,
                label_asym_id: None,
            },
        });
        s.residues.push(ResidueRow {
            chain_id: ChainId::new(0),
            atom_span: RowSpan::new(0, 1),
            name: ResidueName([b'A', b'L', b'A', 0], 3),
            kind: ResidueKind::AminoAcid,
            source: ResidueSourceIds { seq_id: None },
        });
        s.residues.push(ResidueRow {
            chain_id: ChainId::new(0),
            atom_span: RowSpan::new(1, 1),
            name: ResidueName([b'H', b'O', b'H', 0], 3),
            kind: ResidueKind::Water,
            source: ResidueSourceIds { seq_id: None },
        });
        s.atoms.push(AtomRow {
            residue_id: ResidueId::new(0),
            name: AtomName([b' ', b'C', b'A', b' ']),
            element: crate::Element::C,
            altloc: None,
            occupancy: None,
            b_iso: None,
            formal_charge: None,
            anisou: None,
            source: AtomSourceIds { serial: None },
        });
        s.atoms.push(AtomRow {
            residue_id: ResidueId::new(1),
            name: AtomName([b' ', b'O', b' ', b' ']),
            element: crate::Element::O,
            altloc: None,
            occupancy: None,
            b_iso: None,
            formal_charge: None,
            anisou: None,
            source: AtomSourceIds { serial: None },
        });
        s.coordinates.positions = vec![[1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        s
    }

    #[test]
    fn registered_bio_ops_have_matrix_entries() {
        assert_eq!(BIO_STRUCTURE_OPS.len(), 1);
        for operation in BIO_STRUCTURE_OPS {
            assert!(
                BIO_SUPPORT_MATRIX
                    .iter()
                    .any(|entry| std::ptr::eq(entry.operation, *operation)),
                "missing bio support matrix entry for {}",
                operation.method
            );
            assert!(
                BIO_OPERATION_INVARIANT_MATRIX
                    .iter()
                    .any(|entry| std::ptr::eq(entry.operation, *operation)),
                "missing bio invariant matrix entry for {}",
                operation.method
            );
        }
    }

    #[test]
    fn remove_waters_removes_water_residue_and_atom() {
        let s = make_structure_with_waters();
        let result = s.without_waters().expect("remove_waters should succeed");

        assert_eq!(result.num_atoms(), 1);
        assert_eq!(result.num_residues(), 1);
        assert_eq!(result.residues[0].kind, ResidueKind::AminoAcid);
        assert_eq!(result.coordinates.positions, vec![[1.0, 0.0, 0.0]]);
    }

    #[test]
    fn remove_waters_is_noop_on_structure_without_waters() {
        let mut s = BioStructure::new();
        s.models.push(ModelRow {
            chain_span: RowSpan::new(0, 0),
            source_model_number: None,
        });
        let result = s.without_waters().expect("noop should succeed");
        assert_eq!(result.num_atoms(), 0);
    }

    #[test]
    fn remove_waters_preserves_source_invariants() {
        let s = make_structure_with_waters();
        let result = s.without_waters().unwrap();
        crate::bio_invariants::enforce_bio_structure_invariants(&result)
            .expect("result must satisfy invariants");
    }

    #[test]
    fn remove_residues_rejects_weak_operation_specs() {
        let s = BioStructure::new();
        let mut parts = BioOpParts::new(&s, &TEST_WEAK_COMPACT_SPEC);
        let err = parts
            .remove_residues(&[])
            .expect_err("weak operation must not compact hierarchy");
        assert!(matches!(err, BioOperationError::InvalidInput { .. }));
    }

    #[cfg(feature = "op-contracts")]
    #[test]
    #[should_panic(expected = "attempted to mutate a block outside its registry permissions")]
    fn remove_residues_panics_when_registry_does_not_allow_mutation() {
        let s = make_structure_with_waters();
        let mut parts = BioOpParts::new(&s, &TEST_UNAUTHORIZED_REMOVE_RESIDUES_SPEC);
        let water = ResidueId::new(1);
        let _ = parts.remove_residues(&[water]);
    }
}
