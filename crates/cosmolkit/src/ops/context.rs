//! Lazy detached-block access for generated operation capabilities.

use std::{marker::PhantomData, sync::Arc};

use cosmolkit_model::{
    CoordinateBlock, CoordinateDimension, MoleculeProperties, SdfPropertyListTarget, TopologyBlock,
    TopologyMapping,
};

use super::{
    BlockSet, CipStatePolicy, DerivedState, MappingRequirement, MoleculeOpOutput, MoleculeOpSpec,
    OperationError, TopologyEditKind,
};
use crate::Molecule;
use crate::molecule::DerivedCacheBlock;
use crate::strict::{OPERATION_CONTRACTS_ENABLED, RUNTIME_INVARIANTS_ENABLED};

enum WorkingBlock<T> {
    Shared,
    CheckedOut,
    Installed(T),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum PreservationProof {
    UnchangedInput,
    LeafAtomAppend,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct EffectTrace {
    updated: DerivedState,
    cleared: DerivedState,
    preserved: DerivedState,
    cip_applied: bool,
}

impl EffectTrace {
    fn handled(self) -> DerivedState {
        self.updated.union(self.cleared).union(self.preserved)
    }
}

/// One private transaction whose access type is emitted by `molecule_ops!`.
///
/// Construction clones only the Arc-backed `Molecule`. Individual semantic
/// blocks detach lazily when their generated capability checks them out.
pub(crate) struct OpParts<'a, Access> {
    spec: &'static MoleculeOpSpec,
    source: Molecule,
    topology: WorkingBlock<TopologyBlock>,
    coordinates: WorkingBlock<CoordinateBlock>,
    properties: WorkingBlock<MoleculeProperties>,
    derived_cache: WorkingBlock<DerivedCacheBlock>,
    topology_edit: Option<TopologyEditKind>,
    topology_mapping: Option<TopologyMapping>,
    remapped_blocks: BlockSet,
    effect_trace: EffectTrace,
    in_place_target: Option<&'a mut Molecule>,
    access: PhantomData<Access>,
}

impl<'a, Access> OpParts<'a, Access> {
    pub(crate) fn new(
        source: &'a Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        Self::from_source(source.clone(), None, spec)
    }

    pub(crate) fn new_in_place(
        target: &'a mut Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        let source = target.clone();
        Self::from_source(source, Some(target), spec)
    }

    fn from_source(
        source: Molecule,
        in_place_target: Option<&'a mut Molecule>,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if spec.output != MoleculeOpOutput::Single {
            return Err(OperationError::OutputMismatch {
                operation: spec.method,
                expected: MoleculeOpOutput::Single,
                actual: spec.output,
            });
        }
        Self::validate_semantic_preconditions(spec)?;
        Ok(Self {
            spec,
            source,
            topology: WorkingBlock::Shared,
            coordinates: WorkingBlock::Shared,
            properties: WorkingBlock::Shared,
            derived_cache: WorkingBlock::Shared,
            topology_edit: None,
            topology_mapping: None,
            remapped_blocks: BlockSet::NONE,
            effect_trace: EffectTrace::default(),
            in_place_target,
            access: PhantomData,
        })
    }

    pub(super) fn validate_semantic_preconditions(
        spec: &'static MoleculeOpSpec,
    ) -> Result<(), OperationError> {
        if spec.semantic_preconditions.is_empty() {
            return Ok(());
        }
        Err(OperationError::SemanticPreconditionContract {
            operation: spec.method,
            missing: spec.semantic_preconditions,
            issue: "the live runtime has no authoritative provenance evidence for this precondition",
        })
    }

    fn ensure_read_access(
        &self,
        block: BlockSet,
        name: &'static str,
    ) -> Result<(), OperationError> {
        if self.spec.access.can_read(block) {
            Ok(())
        } else {
            Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: name,
            })
        }
    }

    fn ensure_write_access(
        &self,
        block: BlockSet,
        name: &'static str,
    ) -> Result<(), OperationError> {
        if self.spec.access.can_write(block) {
            Ok(())
        } else {
            Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: name,
            })
        }
    }

    pub(super) fn read_topology_runtime(&self) -> Result<&TopologyBlock, OperationError> {
        self.ensure_read_access(BlockSet::TOPOLOGY, "topology")?;
        match &self.topology {
            WorkingBlock::Shared => Ok(self.source.topology()),
            WorkingBlock::Installed(topology) => Ok(topology),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "topology",
            }),
        }
    }

    pub(super) fn read_coordinates_runtime(&self) -> Result<&CoordinateBlock, OperationError> {
        self.ensure_read_access(BlockSet::COORDINATES, "coordinates")?;
        match &self.coordinates {
            WorkingBlock::Shared => Ok(self.source.coordinates()),
            WorkingBlock::Installed(coordinates) => Ok(coordinates),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "coordinates",
            }),
        }
    }

    pub(super) fn read_properties_runtime(&self) -> Result<&MoleculeProperties, OperationError> {
        self.ensure_read_access(BlockSet::PROPERTIES, "properties")?;
        match &self.properties {
            WorkingBlock::Shared => Ok(self.source.properties()),
            WorkingBlock::Installed(properties) => Ok(properties),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "properties",
            }),
        }
    }

    pub(super) fn read_derived_cache_runtime(&self) -> Result<&DerivedCacheBlock, OperationError> {
        self.ensure_read_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        match &self.derived_cache {
            WorkingBlock::Shared => Ok(self.source.derived_cache_runtime()),
            WorkingBlock::Installed(cache) => Ok(cache),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "derived_cache",
            }),
        }
    }

    pub(super) fn checkout_topology_runtime(&mut self) -> Result<TopologyBlock, OperationError> {
        self.ensure_write_access(BlockSet::TOPOLOGY, "topology")?;
        match std::mem::replace(&mut self.topology, WorkingBlock::CheckedOut) {
            WorkingBlock::Shared => Ok(self.source.topology().clone()),
            WorkingBlock::Installed(topology) => Ok(topology),
            WorkingBlock::CheckedOut => {
                self.topology = WorkingBlock::CheckedOut;
                Err(OperationError::BlockCheckedOut {
                    operation: self.spec.method,
                    block: "topology",
                })
            }
        }
    }

    pub(super) fn checkout_coordinates_runtime(
        &mut self,
    ) -> Result<CoordinateBlock, OperationError> {
        self.ensure_write_access(BlockSet::COORDINATES, "coordinates")?;
        match std::mem::replace(&mut self.coordinates, WorkingBlock::CheckedOut) {
            WorkingBlock::Shared => Ok(self.source.coordinates().clone()),
            WorkingBlock::Installed(coordinates) => Ok(coordinates),
            WorkingBlock::CheckedOut => {
                self.coordinates = WorkingBlock::CheckedOut;
                Err(OperationError::BlockCheckedOut {
                    operation: self.spec.method,
                    block: "coordinates",
                })
            }
        }
    }

    pub(super) fn checkout_properties_runtime(
        &mut self,
    ) -> Result<MoleculeProperties, OperationError> {
        self.ensure_write_access(BlockSet::PROPERTIES, "properties")?;
        match std::mem::replace(&mut self.properties, WorkingBlock::CheckedOut) {
            WorkingBlock::Shared => Ok(self.source.properties().clone()),
            WorkingBlock::Installed(properties) => Ok(properties),
            WorkingBlock::CheckedOut => {
                self.properties = WorkingBlock::CheckedOut;
                Err(OperationError::BlockCheckedOut {
                    operation: self.spec.method,
                    block: "properties",
                })
            }
        }
    }

    pub(super) fn checkout_derived_cache_runtime(
        &mut self,
    ) -> Result<DerivedCacheBlock, OperationError> {
        self.ensure_write_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        match std::mem::replace(&mut self.derived_cache, WorkingBlock::CheckedOut) {
            WorkingBlock::Shared => Ok(self.source.derived_cache_runtime().clone()),
            WorkingBlock::Installed(cache) => Ok(cache),
            WorkingBlock::CheckedOut => {
                self.derived_cache = WorkingBlock::CheckedOut;
                Err(OperationError::BlockCheckedOut {
                    operation: self.spec.method,
                    block: "derived_cache",
                })
            }
        }
    }

    fn require_checkout<T>(
        spec: &'static MoleculeOpSpec,
        block: &'static str,
        slot: &mut WorkingBlock<T>,
    ) -> Result<(), OperationError> {
        if matches!(slot, WorkingBlock::CheckedOut) {
            Ok(())
        } else {
            Err(OperationError::BlockNotCheckedOut {
                operation: spec.method,
                block,
            })
        }
    }

    pub(super) fn install_topology_runtime(
        &mut self,
        topology: TopologyBlock,
    ) -> Result<(), OperationError> {
        self.ensure_write_access(BlockSet::TOPOLOGY, "topology")?;
        Self::require_checkout(self.spec, "topology", &mut self.topology)?;
        topology
            .validate()
            .map_err(OperationError::InvalidTopology)?;
        self.topology = WorkingBlock::Installed(topology);
        Ok(())
    }

    pub(super) fn install_coordinates_runtime(
        &mut self,
        mut coordinates: CoordinateBlock,
    ) -> Result<(), OperationError> {
        self.ensure_write_access(BlockSet::COORDINATES, "coordinates")?;
        Self::require_checkout(self.spec, "coordinates", &mut self.coordinates)?;
        let local_row_count = coordinates
            .conformers_2d
            .first()
            .map(|conformer| conformer.coordinates().len())
            .or_else(|| {
                coordinates
                    .conformers_3d
                    .first()
                    .map(|conformer| conformer.coordinates().len())
            })
            .unwrap_or(0);
        coordinates
            .validate_for_atom_count(local_row_count)
            .map_err(OperationError::InvalidCoordinates)?;
        coordinates.source_coordinate_dim = if coordinates.conformers_3d.is_empty() {
            (!coordinates.conformers_2d.is_empty()).then_some(CoordinateDimension::TwoD)
        } else {
            Some(CoordinateDimension::ThreeD)
        };
        self.coordinates = WorkingBlock::Installed(coordinates);
        Ok(())
    }

    pub(super) fn install_properties_runtime(
        &mut self,
        properties: MoleculeProperties,
    ) -> Result<(), OperationError> {
        self.ensure_write_access(BlockSet::PROPERTIES, "properties")?;
        Self::require_checkout(self.spec, "properties", &mut self.properties)?;
        self.properties = WorkingBlock::Installed(properties);
        Ok(())
    }

    pub(super) fn install_derived_cache_runtime(
        &mut self,
        cache: DerivedCacheBlock,
    ) -> Result<(), OperationError> {
        self.ensure_write_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        Self::require_checkout(self.spec, "derived_cache", &mut self.derived_cache)?;
        self.derived_cache = WorkingBlock::Installed(cache);
        Ok(())
    }

    pub(super) fn record_topology_edit_runtime(
        &mut self,
        edit: TopologyEditKind,
    ) -> Result<(), OperationError> {
        if self.topology_edit.is_some() {
            return Err(OperationError::TopologyEditContract {
                operation: self.spec.method,
                issue: "duplicate-record",
                expected: self.spec.topology_edit,
                actual: Some(edit),
            });
        }
        if edit != self.spec.topology_edit || edit == TopologyEditKind::None {
            return Err(OperationError::TopologyEditContract {
                operation: self.spec.method,
                issue: "declaration-mismatch",
                expected: self.spec.topology_edit,
                actual: Some(edit),
            });
        }
        self.topology_edit = Some(edit);
        Ok(())
    }

    pub(super) fn record_topology_mapping_runtime(
        &mut self,
        mapping: TopologyMapping,
    ) -> Result<(), OperationError> {
        if self.topology_mapping.is_some() {
            return Err(OperationError::MappingContract {
                operation: self.spec.method,
                issue: "duplicate mapping record",
                requirement: self.spec.requires_mapping,
            });
        }
        if self.spec.requires_mapping == MappingRequirement::None {
            return Err(OperationError::MappingContract {
                operation: self.spec.method,
                issue: "mapping was not expected",
                requirement: self.spec.requires_mapping,
            });
        }
        self.topology_mapping = Some(mapping);
        Ok(())
    }

    fn validate_property_lists(
        properties: &MoleculeProperties,
        atom_count: usize,
        bond_count: usize,
    ) -> Result<(), OperationError> {
        for property_list in properties.sdf_property_lists() {
            let (target, expected) = match property_list.target() {
                SdfPropertyListTarget::Atom => ("atom", atom_count),
                SdfPropertyListTarget::Bond => ("bond", bond_count),
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

    fn validate_mapping_obligation(
        spec: &'static MoleculeOpSpec,
        before: &TopologyBlock,
        after: &TopologyBlock,
        actual_edit: Option<TopologyEditKind>,
        mapping: Option<&TopologyMapping>,
    ) -> Result<(), OperationError> {
        let expected_edit = spec.topology_edit;
        let expected_actual = (expected_edit != TopologyEditKind::None).then_some(expected_edit);
        if actual_edit != expected_actual {
            return Err(OperationError::TopologyEditContract {
                operation: spec.method,
                issue: if actual_edit.is_none() {
                    "missing-record"
                } else {
                    "declaration-mismatch"
                },
                expected: expected_edit,
                actual: actual_edit,
            });
        }

        let row_identity_unchanged = before.atoms.len() == after.atoms.len()
            && before.bonds.len() == after.bonds.len()
            && before
                .atoms
                .iter()
                .zip(&after.atoms)
                .all(|(source, candidate)| source.id() == candidate.id())
            && before
                .bonds
                .iter()
                .zip(&after.bonds)
                .all(|(source, candidate)| source.id() == candidate.id());
        if expected_edit == TopologyEditKind::None && !row_identity_unchanged {
            return Err(OperationError::TopologyEditContract {
                operation: spec.method,
                issue: "undeclared-topology-change",
                expected: expected_edit,
                actual: actual_edit,
            });
        }

        let old_atom_count = before.atoms.len();
        let new_atom_count = after.atoms.len();
        let old_bond_count = before.bonds.len();
        let new_bond_count = after.bonds.len();
        match spec.requires_mapping {
            MappingRequirement::None => {
                if mapping.is_some() {
                    return Err(OperationError::MappingContract {
                        operation: spec.method,
                        issue: "mapping was not expected",
                        requirement: spec.requires_mapping,
                    });
                }
                if old_atom_count != new_atom_count || old_bond_count != new_bond_count {
                    return Err(OperationError::MappingContract {
                        operation: spec.method,
                        issue: "row counts changed without a mapping",
                        requirement: spec.requires_mapping,
                    });
                }
            }
            MappingRequirement::Identity | MappingRequirement::Required => {
                let mapping = mapping.ok_or(OperationError::MappingContract {
                    operation: spec.method,
                    issue: "required mapping was not recorded",
                    requirement: spec.requires_mapping,
                })?;
                mapping
                    .validate_for_counts(
                        old_atom_count,
                        new_atom_count,
                        old_bond_count,
                        new_bond_count,
                    )
                    .map_err(|source| OperationError::InvalidTopologyMapping {
                        operation: spec.method,
                        source,
                    })?;
                if spec.requires_mapping == MappingRequirement::Identity
                    && *mapping != TopologyMapping::identity(old_atom_count, old_bond_count)
                {
                    return Err(OperationError::MappingContract {
                        operation: spec.method,
                        issue: "mapping is not identity",
                        requirement: spec.requires_mapping,
                    });
                }
            }
        }
        Ok(())
    }

    pub(super) fn apply_runtime_remap_runtime(&mut self) -> Result<(), OperationError> {
        let source_topology = self.source.topology();
        source_topology
            .validate()
            .map_err(OperationError::InvalidTopology)?;
        self.source
            .coordinates()
            .validate_for_atom_count(source_topology.atoms.len())
            .map_err(OperationError::InvalidCoordinates)?;
        Self::validate_property_lists(
            self.source.properties(),
            source_topology.atoms.len(),
            source_topology.bonds.len(),
        )?;

        let candidate_topology = match &self.topology {
            WorkingBlock::Shared => source_topology,
            WorkingBlock::Installed(topology) => topology,
            WorkingBlock::CheckedOut => {
                return Err(OperationError::BlockCheckedOut {
                    operation: self.spec.method,
                    block: "topology",
                });
            }
        };
        candidate_topology
            .validate()
            .map_err(OperationError::InvalidTopology)?;
        Self::validate_mapping_obligation(
            self.spec,
            source_topology,
            candidate_topology,
            self.topology_edit,
            self.topology_mapping.as_ref(),
        )?;

        if self.spec.auto_remap.intersects(BlockSet::TOPOLOGY)
            || self.spec.auto_remap.intersects(BlockSet::DERIVED_CACHE)
        {
            return Err(OperationError::AutoRemapContract {
                operation: self.spec.method,
                block: "non-dependent-state",
                issue: "only coordinates and properties support generic topology remap",
            });
        }

        let mapping = self.topology_mapping.as_ref();
        let new_atom_count = candidate_topology.atoms.len();
        let new_bond_count = candidate_topology.bonds.len();

        let coordinate_candidate = if self.spec.auto_remap.contains(BlockSet::COORDINATES) {
            self.ensure_write_access(BlockSet::COORDINATES, "coordinates")?;
            if self.remapped_blocks.contains(BlockSet::COORDINATES) {
                return Err(OperationError::AutoRemapContract {
                    operation: self.spec.method,
                    block: "coordinates",
                    issue: "remap was already applied",
                });
            }
            let mut coordinates = match &self.coordinates {
                WorkingBlock::Shared => self.source.coordinates().clone(),
                WorkingBlock::Installed(coordinates) => {
                    coordinates
                        .validate_for_atom_count(new_atom_count)
                        .map_err(OperationError::InvalidCoordinates)?;
                    coordinates.clone()
                }
                WorkingBlock::CheckedOut => {
                    return Err(OperationError::BlockCheckedOut {
                        operation: self.spec.method,
                        block: "coordinates",
                    });
                }
            };
            if matches!(&self.coordinates, WorkingBlock::Shared) {
                let mapping = mapping.ok_or(OperationError::MappingContract {
                    operation: self.spec.method,
                    issue: "auto-remap requires a recorded mapping",
                    requirement: self.spec.requires_mapping,
                })?;
                let has_appended_atoms = mapping.atoms().new_to_old().iter().any(Option::is_none);
                let has_conformers =
                    !coordinates.conformers_2d.is_empty() || !coordinates.conformers_3d.is_empty();
                if has_appended_atoms && has_conformers {
                    return Err(OperationError::CoordinateAppendRequiresValues {
                        operation: self.spec.method,
                    });
                }
                coordinates.remap_topology(&mapping.retained_atom_indices());
            }
            coordinates
                .validate_for_atom_count(new_atom_count)
                .map_err(OperationError::InvalidCoordinates)?;
            Some(coordinates)
        } else {
            None
        };

        let property_candidate = if self.spec.auto_remap.contains(BlockSet::PROPERTIES) {
            self.ensure_write_access(BlockSet::PROPERTIES, "properties")?;
            if self.remapped_blocks.contains(BlockSet::PROPERTIES) {
                return Err(OperationError::AutoRemapContract {
                    operation: self.spec.method,
                    block: "properties",
                    issue: "remap was already applied",
                });
            }
            let mut properties = match &self.properties {
                WorkingBlock::Shared => self.source.properties().clone(),
                WorkingBlock::Installed(properties) => properties.clone(),
                WorkingBlock::CheckedOut => {
                    return Err(OperationError::BlockCheckedOut {
                        operation: self.spec.method,
                        block: "properties",
                    });
                }
            };
            if matches!(&self.properties, WorkingBlock::Shared) {
                let mapping = mapping.ok_or(OperationError::MappingContract {
                    operation: self.spec.method,
                    issue: "auto-remap requires a recorded mapping",
                    requirement: self.spec.requires_mapping,
                })?;
                properties
                    .remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
            }
            Self::validate_property_lists(&properties, new_atom_count, new_bond_count)?;
            Some(properties)
        } else {
            None
        };

        if let Some(coordinates) = coordinate_candidate {
            self.coordinates = WorkingBlock::Installed(coordinates);
            self.remapped_blocks = self.remapped_blocks.union(BlockSet::COORDINATES);
        }
        if let Some(properties) = property_candidate {
            self.properties = WorkingBlock::Installed(properties);
            self.remapped_blocks = self.remapped_blocks.union(BlockSet::PROPERTIES);
        }
        Ok(())
    }

    fn effect_error(
        spec: &'static MoleculeOpSpec,
        action: &'static str,
        states: DerivedState,
        issue: &'static str,
    ) -> OperationError {
        OperationError::DerivedEffectContract {
            operation: spec.method,
            action,
            states,
            issue,
        }
    }

    pub(super) fn validate_effect_contract(
        spec: &'static MoleculeOpSpec,
    ) -> Result<(), OperationError> {
        let effects = spec.derived_effects;
        let overlap = effects
            .recompute
            .intersection(effects.preserve)
            .union(effects.recompute.intersection(effects.invalidate))
            .union(effects.recompute.intersection(effects.operation_defined))
            .union(effects.preserve.intersection(effects.invalidate))
            .union(effects.preserve.intersection(effects.operation_defined))
            .union(effects.invalidate.intersection(effects.operation_defined));
        if !overlap.is_empty() {
            return Err(Self::effect_error(
                spec,
                "declaration",
                overlap,
                "effect categories overlap",
            ));
        }

        let all = effects
            .recompute
            .union(effects.preserve)
            .union(effects.invalidate)
            .union(effects.operation_defined);
        if !all.is_empty()
            && (!spec.access.can_write(BlockSet::DERIVED_CACHE)
                || !spec.may_mutate.contains(BlockSet::DERIVED_CACHE))
        {
            return Err(Self::effect_error(
                spec,
                "declaration",
                all,
                "derived effects require derived_cache write and may_mutate authority",
            ));
        }

        if !effects.operation_defined.is_empty()
            && (!(matches!(
                spec.method,
                "without_hydrogens" | "without_hydrogens_with_params"
            )) || effects.operation_defined != DerivedState::VALENCE)
        {
            return Err(Self::effect_error(
                spec,
                "operation_defined",
                effects.operation_defined,
                "only valence in the hydrogen-removal family is allow-listed",
            ));
        }

        let writes_cip_blocks = spec.access.can_write(BlockSet::TOPOLOGY)
            && spec.access.can_write(BlockSet::PROPERTIES)
            && spec.may_mutate.contains(BlockSet::TOPOLOGY)
            && spec.may_mutate.contains(BlockSet::PROPERTIES);
        match spec.cip_state {
            CipStatePolicy::Preserve => {}
            CipStatePolicy::ClearComputed if !writes_cip_blocks => {
                return Err(OperationError::CipStateContract {
                    operation: spec.method,
                    policy: spec.cip_state,
                    issue: "clear requires topology and properties write and may_mutate authority",
                });
            }
            CipStatePolicy::Assign
                if spec.method != "with_cip_labels_with_options" || !writes_cip_blocks =>
            {
                return Err(OperationError::CipStateContract {
                    operation: spec.method,
                    policy: spec.cip_state,
                    issue: "assign is reserved for with_cip_labels_with_options with topology and properties write authority",
                });
            }
            CipStatePolicy::TautomerSourceTransition
                if spec.method != "enumerate_tautomers_with_options"
                    || spec.output != MoleculeOpOutput::Multiple
                    || !writes_cip_blocks =>
            {
                return Err(OperationError::CipStateContract {
                    operation: spec.method,
                    policy: spec.cip_state,
                    issue: "tautomer transition requires the exact multiple-output operation and topology/properties write authority",
                });
            }
            _ => {}
        }
        Ok(())
    }

    fn validate_effect_action(
        &self,
        action: &'static str,
        states: DerivedState,
        allowed: DerivedState,
    ) -> Result<(), OperationError> {
        Self::validate_effect_contract(self.spec)?;
        if states.is_empty() {
            return Err(Self::effect_error(
                self.spec,
                action,
                states,
                "empty action masks cannot satisfy an obligation",
            ));
        }
        let undeclared = states.difference(allowed);
        if !undeclared.is_empty() {
            return Err(Self::effect_error(
                self.spec,
                action,
                undeclared,
                "state is not declared for this action",
            ));
        }
        let duplicate = states.intersection(self.effect_trace.handled());
        if !duplicate.is_empty() {
            return Err(Self::effect_error(
                self.spec,
                action,
                duplicate,
                "state was already handled",
            ));
        }
        Ok(())
    }

    fn current_cache_candidate(&self) -> Result<DerivedCacheBlock, OperationError> {
        match &self.derived_cache {
            WorkingBlock::Shared => Ok(self.source.derived_cache_runtime().clone()),
            WorkingBlock::Installed(cache) => Ok(cache.clone()),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "derived_cache",
            }),
        }
    }

    pub(super) fn mark_cache_updated_runtime(
        &mut self,
        states: DerivedState,
    ) -> Result<(), OperationError> {
        let allowed = self
            .spec
            .derived_effects
            .recompute
            .union(self.spec.derived_effects.operation_defined);
        self.validate_effect_action("update", states, allowed)?;
        self.ensure_write_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        let mut candidate = self.current_cache_candidate()?;
        candidate.mark_valid(states);
        self.derived_cache = WorkingBlock::Installed(candidate);
        self.effect_trace.updated = self.effect_trace.updated.union(states);
        Ok(())
    }

    pub(super) fn clear_cache_runtime(
        &mut self,
        states: DerivedState,
    ) -> Result<(), OperationError> {
        let allowed = self
            .spec
            .derived_effects
            .recompute
            .union(self.spec.derived_effects.invalidate)
            .union(self.spec.derived_effects.operation_defined);
        self.validate_effect_action("clear", states, allowed)?;
        self.ensure_write_access(BlockSet::DERIVED_CACHE, "derived_cache")?;
        let mut candidate = self.current_cache_candidate()?;
        candidate.clear(states);
        self.derived_cache = WorkingBlock::Installed(candidate);
        self.effect_trace.cleared = self.effect_trace.cleared.union(states);
        Ok(())
    }

    fn current_topology_candidate(&self) -> Result<&TopologyBlock, OperationError> {
        match &self.topology {
            WorkingBlock::Shared => Ok(self.source.topology()),
            WorkingBlock::Installed(topology) => Ok(topology),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "topology",
            }),
        }
    }

    fn current_coordinates_candidate(&self) -> Result<&CoordinateBlock, OperationError> {
        match &self.coordinates {
            WorkingBlock::Shared => Ok(self.source.coordinates()),
            WorkingBlock::Installed(coordinates) => Ok(coordinates),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "coordinates",
            }),
        }
    }

    fn current_properties_candidate(&self) -> Result<&MoleculeProperties, OperationError> {
        match &self.properties {
            WorkingBlock::Shared => Ok(self.source.properties()),
            WorkingBlock::Installed(properties) => Ok(properties),
            WorkingBlock::CheckedOut => Err(OperationError::BlockCheckedOut {
                operation: self.spec.method,
                block: "properties",
            }),
        }
    }

    pub(super) fn prove_preserved_runtime(
        &mut self,
        states: DerivedState,
        proof: PreservationProof,
    ) -> Result<(), OperationError> {
        self.validate_effect_action("preserve", states, self.spec.derived_effects.preserve)?;
        match proof {
            PreservationProof::UnchangedInput => {
                let unchanged = self.current_topology_candidate()? == self.source.topology()
                    && self.current_coordinates_candidate()? == self.source.coordinates()
                    && self.current_properties_candidate()? == self.source.properties()
                    && self
                        .current_cache_candidate()?
                        .valid_states()
                        .intersection(states)
                        == self
                            .source
                            .derived_cache_runtime()
                            .valid_states()
                            .intersection(states);
                if !unchanged {
                    return Err(Self::effect_error(
                        self.spec,
                        "preserve",
                        states,
                        "unchanged-input proof failed",
                    ));
                }
            }
            PreservationProof::LeafAtomAppend => {
                return Err(Self::effect_error(
                    self.spec,
                    "preserve",
                    states,
                    "leaf-atom-append proof is not implemented by this runtime unit",
                ));
            }
        }
        self.effect_trace.preserved = self.effect_trace.preserved.union(states);
        Ok(())
    }

    fn cip_properties_equal(&self) -> Result<bool, OperationError> {
        const ATOM_KEYS: [&str; 3] = ["_CIPCode", "_CIPNeighborOrder", "_CIPRank"];
        const BOND_KEYS: [&str; 2] = ["_CIPCode", "_CIPNeighborOrder"];
        let candidate_topology = self.current_topology_candidate()?;
        let candidate_properties = self.current_properties_candidate()?;
        if candidate_topology.atoms.len() != self.source.topology().atoms.len()
            || candidate_topology.bonds.len() != self.source.topology().bonds.len()
        {
            return Ok(false);
        }
        for (candidate, source) in candidate_topology
            .atoms
            .iter()
            .zip(&self.source.topology().atoms)
        {
            for key in ATOM_KEYS {
                if candidate.prop(key) != source.prop(key)
                    || candidate.is_prop_computed(key) != source.is_prop_computed(key)
                {
                    return Ok(false);
                }
            }
        }
        for (candidate, source) in candidate_topology
            .bonds
            .iter()
            .zip(&self.source.topology().bonds)
        {
            for key in BOND_KEYS {
                if candidate.prop(key) != source.prop(key)
                    || candidate.is_prop_computed(key) != source.is_prop_computed(key)
                {
                    return Ok(false);
                }
            }
        }
        Ok(candidate_properties.prop("_CIPComputed")
            == self.source.properties().prop("_CIPComputed")
            && candidate_properties.is_prop_computed("_CIPComputed")
                == self.source.properties().is_prop_computed("_CIPComputed"))
    }

    pub(super) fn apply_cip_policy_runtime(&mut self) -> Result<(), OperationError> {
        Self::validate_effect_contract(self.spec)?;
        if self.effect_trace.cip_applied {
            return Err(OperationError::CipStateContract {
                operation: self.spec.method,
                policy: self.spec.cip_state,
                issue: "CIP transition was already applied",
            });
        }
        match self.spec.cip_state {
            CipStatePolicy::Preserve => {
                if !self.cip_properties_equal()? {
                    return Err(OperationError::CipStateContract {
                        operation: self.spec.method,
                        policy: self.spec.cip_state,
                        issue: "candidate changed CIP-observable state",
                    });
                }
            }
            CipStatePolicy::ClearComputed => {
                self.ensure_write_access(BlockSet::TOPOLOGY, "topology")?;
                self.ensure_write_access(BlockSet::PROPERTIES, "properties")?;
                let mut topology = self.current_topology_candidate()?.clone();
                let mut properties = self.current_properties_candidate()?.clone();
                for atom in &mut topology.atoms {
                    atom.clear_computed_props();
                }
                for bond in &mut topology.bonds {
                    bond.clear_computed_props();
                }
                properties.clear_computed_props();
                self.topology = WorkingBlock::Installed(topology);
                self.properties = WorkingBlock::Installed(properties);
            }
            CipStatePolicy::Assign => {
                let properties = self.current_properties_candidate()?;
                if properties.prop("_CIPComputed").is_none()
                    || !properties.is_prop_computed("_CIPComputed")
                {
                    return Err(OperationError::CipStateContract {
                        operation: self.spec.method,
                        policy: self.spec.cip_state,
                        issue: "assignment did not install computed _CIPComputed evidence",
                    });
                }
            }
            CipStatePolicy::TautomerSourceTransition => {
                // The source-backed tautomer owner produces the complete
                // candidate topology and properties. The generic runtime does
                // not reproduce that chemistry; it only enforces the exact
                // declaration allow-list above and validates the resulting
                // detached candidate before construction.
            }
        }
        self.effect_trace.cip_applied = true;
        Ok(())
    }

    fn has_staged_state(&self) -> bool {
        !matches!(self.topology, WorkingBlock::Shared)
            || !matches!(self.coordinates, WorkingBlock::Shared)
            || !matches!(self.properties, WorkingBlock::Shared)
            || !matches!(self.derived_cache, WorkingBlock::Shared)
            || !self.effect_trace.handled().is_empty()
    }

    fn validate_effect_completion(&self) -> Result<(), OperationError> {
        Self::validate_effect_contract(self.spec)?;
        if !self.has_staged_state() {
            return Ok(());
        }
        let effects = self.spec.derived_effects;
        let updated_or_cleared = self.effect_trace.updated.union(self.effect_trace.cleared);
        for (action, missing) in [
            (
                "recompute",
                effects.recompute.difference(updated_or_cleared),
            ),
            (
                "preserve",
                effects.preserve.difference(self.effect_trace.preserved),
            ),
            (
                "invalidate",
                effects.invalidate.difference(self.effect_trace.cleared),
            ),
            (
                "operation_defined",
                effects.operation_defined.difference(updated_or_cleared),
            ),
        ] {
            if !missing.is_empty() {
                return Err(Self::effect_error(
                    self.spec,
                    action,
                    missing,
                    "declared effect was not completed",
                ));
            }
        }
        if !self.effect_trace.cip_applied {
            return Err(OperationError::CipStateContract {
                operation: self.spec.method,
                policy: self.spec.cip_state,
                issue: "CIP policy was not applied",
            });
        }
        Ok(())
    }

    fn validate_operation_contract(&self) -> Result<(), OperationError> {
        Self::validate_operation_spec(self.spec)
    }

    pub(super) fn validate_operation_spec(
        spec: &'static MoleculeOpSpec,
    ) -> Result<(), OperationError> {
        let read_write_overlap = spec.access.read().intersection(spec.access.write());
        if !read_write_overlap.is_empty() {
            return Err(OperationError::OperationContract {
                operation: spec.method,
                field: "access",
                issue: "read and write sets overlap",
                expected: 0,
                actual: read_write_overlap.bits(),
            });
        }
        if spec.may_mutate != spec.access.write() {
            return Err(OperationError::OperationContract {
                operation: spec.method,
                field: "may_mutate",
                issue: "must equal access.write",
                expected: spec.access.write().bits(),
                actual: spec.may_mutate.bits(),
            });
        }
        let unauthorized_remap = spec.auto_remap.difference(spec.access.write());
        if !unauthorized_remap.is_empty() {
            return Err(OperationError::OperationContract {
                operation: spec.method,
                field: "auto_remap",
                issue: "contains a block without write authority",
                expected: spec.access.write().bits(),
                actual: spec.auto_remap.bits(),
            });
        }
        Ok(())
    }

    fn ensure_complete_blocks(&self) -> Result<(), OperationError> {
        for (block, checked_out) in [
            (
                "topology",
                matches!(self.topology, WorkingBlock::CheckedOut),
            ),
            (
                "coordinates",
                matches!(self.coordinates, WorkingBlock::CheckedOut),
            ),
            (
                "properties",
                matches!(self.properties, WorkingBlock::CheckedOut),
            ),
            (
                "derived_cache",
                matches!(self.derived_cache, WorkingBlock::CheckedOut),
            ),
        ] {
            if checked_out {
                return Err(OperationError::IncompleteCommit {
                    operation: self.spec.method,
                    block,
                });
            }
        }
        Ok(())
    }

    fn validate_mapping_completion(&self) -> Result<(), OperationError> {
        let candidate_topology = self.current_topology_candidate()?;
        Self::validate_mapping_obligation(
            self.spec,
            self.source.topology(),
            candidate_topology,
            self.topology_edit,
            self.topology_mapping.as_ref(),
        )?;

        let missing = self.spec.auto_remap.difference(self.remapped_blocks);
        if !missing.is_empty() {
            return Err(OperationError::OperationContract {
                operation: self.spec.method,
                field: "auto_remap",
                issue: "declared remap was not completed",
                expected: self.spec.auto_remap.bits(),
                actual: self.remapped_blocks.bits(),
            });
        }
        let unexpected = self.remapped_blocks.difference(self.spec.auto_remap);
        if !unexpected.is_empty() {
            return Err(OperationError::OperationContract {
                operation: self.spec.method,
                field: "auto_remap",
                issue: "an undeclared remap was recorded",
                expected: self.spec.auto_remap.bits(),
                actual: self.remapped_blocks.bits(),
            });
        }
        Ok(())
    }

    fn validate_detached_candidate_invariants(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Result<(), OperationError> {
        topology
            .validate()
            .map_err(OperationError::InvalidTopology)?;
        coordinates
            .validate_for_atom_count(topology.atoms.len())
            .map_err(OperationError::InvalidCoordinates)?;
        Self::validate_property_lists(properties, topology.atoms.len(), topology.bonds.len())
    }

    fn validate_candidate_invariants(&self) -> Result<(), OperationError> {
        Self::validate_detached_candidate_invariants(
            self.current_topology_candidate()?,
            self.current_coordinates_candidate()?,
            self.current_properties_candidate()?,
        )
    }

    fn validate_candidate(&self) -> Result<(), OperationError> {
        self.ensure_complete_blocks()?;
        if OPERATION_CONTRACTS_ENABLED {
            self.validate_operation_contract()?;
            self.validate_mapping_completion()?;
            self.validate_effect_completion()?;
        }
        if RUNTIME_INVARIANTS_ENABLED {
            self.validate_candidate_invariants()?;
        }
        Ok(())
    }

    fn materialize<T>(
        slot: WorkingBlock<T>,
        source: Arc<T>,
        operation: &'static str,
        block: &'static str,
    ) -> Result<Arc<T>, OperationError> {
        match slot {
            WorkingBlock::Shared => Ok(source),
            WorkingBlock::Installed(value) => Ok(Arc::new(value)),
            WorkingBlock::CheckedOut => Err(OperationError::IncompleteCommit { operation, block }),
        }
    }

    pub(crate) fn finish(self) -> Result<Molecule, OperationError> {
        let (topology, coordinates, properties, derived_cache) = self.finish_parts()?;
        Molecule::from_runtime_parts(topology, coordinates, properties, derived_cache)
    }

    fn finish_parts(
        self,
    ) -> Result<
        (
            Arc<TopologyBlock>,
            Arc<CoordinateBlock>,
            Arc<MoleculeProperties>,
            Arc<DerivedCacheBlock>,
        ),
        OperationError,
    > {
        self.validate_candidate()?;
        let operation = self.spec.method;
        let topology = Self::materialize(
            self.topology,
            self.source.topology_arc_runtime(),
            operation,
            "topology",
        )?;
        let coordinates = Self::materialize(
            self.coordinates,
            self.source.coordinates_arc_runtime(),
            operation,
            "coordinates",
        )?;
        let properties = Self::materialize(
            self.properties,
            self.source.properties_arc_runtime(),
            operation,
            "properties",
        )?;
        let derived_cache = Self::materialize(
            self.derived_cache,
            self.source.derived_cache_arc_runtime(),
            operation,
            "derived_cache",
        )?;
        Ok((topology, coordinates, properties, derived_cache))
    }

    pub(crate) fn abort_in_place(&mut self) {
        if self.in_place_target.is_none() {
            return;
        }

        // The authoritative target is not replaced until `finish_in_place`.
        // Explicitly discard every staged value and its contract evidence so
        // the body-error path cannot retain a checked-out block, mapping, or
        // derived-state transition and cannot accidentally be finalized by a
        // later internal call. The Arc-backed source is the same cheap COW
        // input used by the operation body; no second rollback runtime or
        // mutable live-state handle is created here.
        self.topology = WorkingBlock::Shared;
        self.coordinates = WorkingBlock::Shared;
        self.properties = WorkingBlock::Shared;
        self.derived_cache = WorkingBlock::Shared;
        self.topology_edit = None;
        self.topology_mapping = None;
        self.remapped_blocks = BlockSet::NONE;
        self.effect_trace = EffectTrace::default();
    }

    pub(crate) fn finish_in_place(mut self) -> Result<(), OperationError> {
        let target = self
            .in_place_target
            .take()
            .ok_or(OperationError::IncompleteCommit {
                operation: self.spec.method,
                block: "in_place_target",
            })?;
        let replacement = self.finish()?;
        *target = replacement;
        Ok(())
    }
}

pub(super) fn validate_multiple_candidate(
    source: &Molecule,
    spec: &'static MoleculeOpSpec,
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<
    (
        Arc<TopologyBlock>,
        Arc<CoordinateBlock>,
        Arc<MoleculeProperties>,
        Arc<DerivedCacheBlock>,
    ),
    OperationError,
> {
    OpParts::<()>::validate_operation_spec(spec)?;

    let row_identity_unchanged = source.topology().atoms.len() == topology.atoms.len()
        && source.topology().bonds.len() == topology.bonds.len()
        && source
            .topology()
            .atoms
            .iter()
            .zip(&topology.atoms)
            .all(|(before, after)| before.id() == after.id())
        && source
            .topology()
            .bonds
            .iter()
            .zip(&topology.bonds)
            .all(|(before, after)| before.id() == after.id());
    let mapping = match spec.requires_mapping {
        MappingRequirement::None => None,
        MappingRequirement::Identity | MappingRequirement::Required if row_identity_unchanged => {
            Some(TopologyMapping::identity(
                topology.atoms.len(),
                topology.bonds.len(),
            ))
        }
        MappingRequirement::Identity | MappingRequirement::Required => {
            return Err(OperationError::MappingContract {
                operation: spec.method,
                issue: "multiple-output tuple has no non-identity mapping evidence",
                requirement: spec.requires_mapping,
            });
        }
    };

    // A declared identity/required mapping conflict is classified before
    // inspecting the candidate's local row ids because the tuple carries no
    // non-identity mapping evidence. All other malformed detached tuples keep
    // the shared runtime's exact invariant error ahead of write projection.
    if RUNTIME_INVARIANTS_ENABLED {
        OpParts::<()>::validate_detached_candidate_invariants(
            &topology,
            &coordinates,
            &properties,
        )?;
    }

    for (changed, block, name) in [
        (
            topology != *source.topology(),
            BlockSet::TOPOLOGY,
            "topology",
        ),
        (
            coordinates != *source.coordinates(),
            BlockSet::COORDINATES,
            "coordinates",
        ),
        (
            properties != *source.properties(),
            BlockSet::PROPERTIES,
            "properties",
        ),
    ] {
        if changed && !spec.access.can_write(block) {
            return Err(OperationError::AccessDenied {
                operation: spec.method,
                block: name,
            });
        }
    }

    let mut candidate = OpParts::<()> {
        spec,
        source: source.clone(),
        topology: WorkingBlock::Installed(topology),
        coordinates: WorkingBlock::Installed(coordinates),
        properties: WorkingBlock::Installed(properties),
        derived_cache: WorkingBlock::Shared,
        topology_edit: (spec.topology_edit != TopologyEditKind::None).then_some(spec.topology_edit),
        topology_mapping: mapping,
        remapped_blocks: spec.auto_remap,
        effect_trace: EffectTrace::default(),
        in_place_target: None,
        access: PhantomData,
    };

    OpParts::<()>::validate_effect_contract(spec)?;
    let clear = spec
        .derived_effects
        .recompute
        .union(spec.derived_effects.invalidate)
        .union(spec.derived_effects.operation_defined);
    if !clear.is_empty() {
        candidate.clear_cache_runtime(clear)?;
    }
    candidate.apply_cip_policy_runtime()?;
    if !spec.derived_effects.preserve.is_empty() {
        candidate.prove_preserved_runtime(
            spec.derived_effects.preserve,
            PreservationProof::UnchangedInput,
        )?;
    }
    candidate.finish_parts()
}

#[cfg(test)]
#[path = "../../tests/support/run_access_internal.rs"]
mod tests;

#[cfg(test)]
#[path = "../../tests/support/run_mapping_internal.rs"]
mod mapping_tests;

#[cfg(test)]
#[path = "../../tests/support/run_effects_internal.rs"]
mod effects_tests;

#[cfg(all(test, feature = "op-contracts-strict"))]
#[path = "../../tests/support/run_commit_internal.rs"]
mod commit_tests;

#[cfg(all(test, feature = "op-contracts-strict"))]
#[path = "../../tests/support/run_failure_internal.rs"]
mod failure_tests;
