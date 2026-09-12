#![cfg(not(feature = "hydrogens"))]

use std::sync::Arc;

use cosmolkit_macros::{mol_op_body, molecule_ops};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

pub use cosmolkit::ops::{
    CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, OperationInvariantEntry,
    ParityMatrixEntry, ParityPolicy, SemanticPreconditionSet, SupportMatrixEntry,
};
pub use cosmolkit::{
    BlockAccess, BlockSet, FeatureSpec, MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec,
    OperationDomain, OperationError, SupportStatus, TopologyEditKind, UnsupportedFeatureError,
};

mod ops {
    pub use crate::{
        BlockAccess, BlockSet, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement,
        MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec, OperationDomain, OperationError,
        OperationInvariantEntry, ParityMatrixEntry, ParityPolicy, SemanticPreconditionSet,
        SupportMatrixEntry, SupportStatus, TopologyEditKind, UnsupportedFeatureError,
    };
}

mod strict {
    pub(crate) const RUNTIME_INVARIANTS_ENABLED: bool = cfg!(feature = "runtime-invariants");
    pub(crate) const OPERATION_CONTRACTS_ENABLED: bool = cfg!(feature = "op-contracts");
}

mod molecule {
    use super::*;

    #[derive(Clone, Debug, Default, Eq, PartialEq)]
    pub(crate) struct DerivedCacheBlock {
        valid: DerivedState,
    }

    impl DerivedCacheBlock {
        pub(crate) fn valid_states(&self) -> DerivedState {
            self.valid
        }

        pub(crate) fn mark_valid(&mut self, states: DerivedState) {
            self.valid = self.valid.union(states);
        }

        pub(crate) fn clear(&mut self, states: DerivedState) {
            self.valid = self.valid.difference(states);
        }
    }

    #[derive(Clone, Debug, PartialEq)]
    struct MoleculeState {
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        derived_cache: DerivedCacheBlock,
        runtime_constructions: usize,
    }

    #[derive(Clone, Debug, PartialEq)]
    pub struct Molecule {
        state: Arc<MoleculeState>,
    }

    impl Molecule {
        pub fn from_parts(
            topology: TopologyBlock,
            coordinates: CoordinateBlock,
            properties: MoleculeProperties,
        ) -> Result<Self, OperationError> {
            topology
                .validate()
                .map_err(OperationError::InvalidTopology)?;
            coordinates
                .validate_for_atom_count(topology.atoms.len())
                .map_err(OperationError::InvalidCoordinates)?;
            Ok(Self {
                state: Arc::new(MoleculeState {
                    topology,
                    coordinates,
                    properties,
                    derived_cache: DerivedCacheBlock::default(),
                    runtime_constructions: 0,
                }),
            })
        }

        pub(crate) fn from_runtime_parts(
            topology: Arc<TopologyBlock>,
            coordinates: Arc<CoordinateBlock>,
            properties: Arc<MoleculeProperties>,
            derived_cache: Arc<DerivedCacheBlock>,
        ) -> Result<Self, OperationError> {
            let mut molecule = Self::from_parts(
                topology.as_ref().clone(),
                coordinates.as_ref().clone(),
                properties.as_ref().clone(),
            )?;
            Arc::get_mut(&mut molecule.state)
                .expect("new runtime state is uniquely owned")
                .derived_cache = derived_cache.as_ref().clone();
            Arc::get_mut(&mut molecule.state)
                .expect("new runtime state is uniquely owned")
                .runtime_constructions = 1;
            Ok(molecule)
        }

        pub fn topology(&self) -> &TopologyBlock {
            &self.state.topology
        }

        pub fn coordinates(&self) -> &CoordinateBlock {
            &self.state.coordinates
        }

        pub fn properties(&self) -> &MoleculeProperties {
            &self.state.properties
        }

        pub fn num_atoms(&self) -> usize {
            self.state.topology.atoms.len()
        }

        pub(crate) fn derived_cache_runtime(&self) -> &DerivedCacheBlock {
            &self.state.derived_cache
        }

        pub(crate) fn topology_arc_runtime(&self) -> Arc<TopologyBlock> {
            Arc::new(self.state.topology.clone())
        }

        pub(crate) fn coordinates_arc_runtime(&self) -> Arc<CoordinateBlock> {
            Arc::new(self.state.coordinates.clone())
        }

        pub(crate) fn properties_arc_runtime(&self) -> Arc<MoleculeProperties> {
            Arc::new(self.state.properties.clone())
        }

        pub(crate) fn derived_cache_arc_runtime(&self) -> Arc<DerivedCacheBlock> {
            Arc::new(self.state.derived_cache.clone())
        }

        pub(crate) fn runtime_constructions(&self) -> usize {
            self.state.runtime_constructions
        }
    }
}

pub use molecule::Molecule;

#[path = "../src/ops/context.rs"]
mod context;
pub(crate) use context::{OpParts, PreservationProof};

pub const SYNTHETIC_FEATURE: FeatureSpec = FeatureSpec {
    name: "synthetic-runtime-access",
    category: "internal-test",
    status: SupportStatus::Experimental,
    rdkit_parity_sensitive: false,
    docs: "Compile-only operation projection used by RUN-access tests.",
};

#[mol_op_body(synthetic_access, parts)]
fn synthetic_access_impl() -> Result<(), OperationError> {
    let _ = parts.topology()?;
    let coordinates = parts.checkout_coordinates()?;
    parts.install_coordinates(coordinates)?;
    parts.apply_cip_policy()
}

molecule_ops! {
    op synthetic_access {
        method: synthetic_access,
        impl_fn: crate::synthetic_access_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [topology], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [], preserve: [], invalidate: [], operation_defined: [],
        },
        cip_state: preserve,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::SYNTHETIC_FEATURE,
        parity: not_applicable,
        io_roundtrip: false,
        invariant_profile: "run-access-synthetic",
    }
}

#[test]
fn generated_marker_spec_matrices_and_wrapper_compile_at_the_runtime_callsite() {
    assert_eq!(SyntheticAccessAccess::__COSMOLKIT_BODY_CLASS, 0);
    assert_eq!(SyntheticAccessAccess::__COSMOLKIT_ACCESS_READ, 1);
    assert_eq!(SyntheticAccessAccess::__COSMOLKIT_ACCESS_WRITE, 2);
    assert_eq!(SYNTHETIC_ACCESS_SPEC.access.read().bits(), 1);
    assert_eq!(SYNTHETIC_ACCESS_SPEC.access.write().bits(), 2);
    assert_eq!(MOLECULE_OPS, &[&SYNTHETIC_ACCESS_SPEC]);
    assert_eq!(SUPPORT_MATRIX.len(), 1);
    assert_eq!(SUPPORT_MATRIX[0].operation, Some(&SYNTHETIC_ACCESS_SPEC));
    assert_eq!(OPERATION_INVARIANT_MATRIX.len(), 1);
    assert_eq!(
        OPERATION_INVARIANT_MATRIX[0],
        OperationInvariantEntry::for_operation(&SYNTHETIC_ACCESS_SPEC, "run-access-synthetic")
    );
    assert!(PARITY_MATRIX.is_empty());

    let source = Molecule::from_parts(
        TopologyBlock::default(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();
    let result = source.synthetic_access().unwrap();
    assert_eq!(result.num_atoms(), 0);
    assert_eq!(source.num_atoms(), 0);
}

#[test]
fn source_guards_enforce_private_fields_and_one_generated_capability_source() {
    let context = include_str!("../src/ops/context.rs");
    let registry = include_str!("../src/ops/registry.rs");
    let hydrogens = include_str!("../src/ops/hydrogens.rs");

    assert!(context.contains("pub(crate) struct OpParts<'a, Access>"));
    for forbidden in [
        "pub struct OpParts",
        "pub(crate) spec:",
        "pub(crate) source:",
        "pub(crate) topology:",
        "pub(crate) coordinates:",
        "pub(crate) properties:",
        "pub(crate) derived_cache:",
        "pub(crate) in_place_target:",
        "pub(crate) fn read_topology_runtime",
        "pub(crate) fn checkout_topology_runtime",
        "pub(crate) fn install_topology_runtime",
        "pub fn working",
        "pub(crate) fn working",
    ] {
        assert!(
            !context.contains(forbidden),
            "found forbidden source: {forbidden}"
        );
    }

    assert_eq!(registry.matches("molecule_ops!").count(), 1);
    assert!(!registry.contains("pub(crate) struct WithHydrogensAccess"));
    assert!(!registry.contains("pub(crate) struct WithoutHydrogensAccess"));
    assert!(!registry.contains("const WITH_HYDROGENS_SPEC: MoleculeOpSpec ="));
    assert!(!registry.contains("const WITHOUT_HYDROGENS_SPEC: MoleculeOpSpec ="));
    assert!(!hydrogens.contains("impl Molecule"));
    assert!(!context.contains("cosmolkit_core"));

    // The declaration macro is the only capability table. Runtime source has
    // no marker-specific grant list that can drift from registry access.
    assert!(!context.contains("impl_all_block_write_access"));
    assert!(!context.contains("WithHydrogensAccess"));
    assert!(!context.contains("WithoutHydrogensAccess"));
}
