//! Public Morgan configuration values for the read-only Molecule facade.

use std::fmt;

use cosmolkit_core::{RingInfo, RingSearchParams, ValenceAssignment, symmetrized_sssr};
use cosmolkit_fingerprints::{
    Fingerprint, FingerprintAdditionalOutput, MorganAtomInvariants, MorganCall, MorganError,
    MorganParams, MorganPreparedInput, SparseBitFingerprint, SparseCountFingerprint,
    SparseCountFingerprint32, morgan_bits, morgan_count, morgan_sparse_bits, morgan_sparse_count,
};

use crate::{DerivedState, FingerprintPreparationError, Molecule, QueryGraph};

/// Typed failure from the public Morgan fingerprint preparation and generation
/// boundary.
#[derive(Debug)]
pub enum MorganReadError {
    /// Shared molecule-state preparation failed before Morgan generation.
    Preparation(FingerprintPreparationError),
    /// Canonical Morgan generation failed with its original owner error.
    Generator(MorganError),
}

impl fmt::Display for MorganReadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Preparation(error) => write!(formatter, "Morgan preparation failed: {error}"),
            Self::Generator(error) => write!(formatter, "Morgan generation failed: {error}"),
        }
    }
}

impl std::error::Error for MorganReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Preparation(error) => Some(error),
            Self::Generator(error) => Some(error),
        }
    }
}

enum MorganReadRings<'a> {
    Cached(&'a RingInfo),
    Temporary(RingInfo),
}

pub(super) struct MorganReadPreparation<'a> {
    molecule: &'a Molecule,
    valence: &'a ValenceAssignment,
    rings: MorganReadRings<'a>,
}

impl MorganReadPreparation<'_> {
    pub(super) fn owner_input(&self) -> MorganPreparedInput<'_> {
        let rings = match &self.rings {
            MorganReadRings::Cached(rings) => *rings,
            MorganReadRings::Temporary(rings) => rings,
        };
        MorganPreparedInput {
            topology: self.molecule.topology(),
            coordinates: self.molecule.coordinate_block_runtime(),
            properties: self.molecule.properties(),
            valence: self.valence,
            rings,
        }
    }
}

pub(super) fn prepare_morgan_read_input(
    molecule: &Molecule,
) -> Result<MorganReadPreparation<'_>, FingerprintPreparationError> {
    let cache = molecule.derived_cache_runtime();
    let valence = cache
        .valence_assignment()
        .ok_or(FingerprintPreparationError::MissingPreparedValence)?;

    let cached_rings = cache
        .valid_states()
        .contains(DerivedState::RINGS)
        .then(|| cache.ring_info())
        .flatten()
        .filter(|rings| rings.is_initialized());
    let rings = match cached_rings {
        Some(rings) => MorganReadRings::Cached(rings),
        None => MorganReadRings::Temporary(
            symmetrized_sssr(molecule.topology(), &RingSearchParams::default())
                .map_err(FingerprintPreparationError::RingPreparation)?,
        ),
    };

    Ok(MorganReadPreparation {
        molecule,
        valence,
        rings,
    })
}

/// Select the atom invariant source for a Morgan fingerprint call.
///
/// `FeaturePatterns` owns the already-built query graphs. The facade passes
/// their borrowed slice to the canonical fingerprint owner without parsing or
/// replacing those values.
#[derive(Debug, Clone, PartialEq)]
pub enum MorganInvariants {
    /// Use the source connectivity invariant generator.
    Connectivity,
    /// Use the source feature invariant generator.
    Features,
    /// Match the provided query graphs as feature patterns.
    FeaturePatterns(Vec<QueryGraph>),
}

/// Configuration shared by the four read-only Morgan fingerprint methods.
///
/// Optional vectors preserve `None` separately from `Some(Vec::new())` and
/// are borrowed when projected into the fingerprint owner's call arguments.
#[derive(Debug, Clone, PartialEq)]
pub struct MorganFingerprintParams {
    /// Canonical Morgan generator options and source defaults.
    pub generator: MorganParams,
    /// Source atom indices at which to start fingerprint expansion.
    pub from_atoms: Option<Vec<u32>>,
    /// Source atom indices to ignore.
    pub ignore_atoms: Option<Vec<u32>>,
    /// Optional caller-supplied atom invariants.
    pub custom_atom_invariants: Option<Vec<u32>>,
    /// Optional caller-supplied bond invariants.
    pub custom_bond_invariants: Option<Vec<u32>>,
    /// Source conformer identifier; the source generator currently ignores it.
    pub conformer_id: i32,
    /// Atom invariant generation policy.
    pub invariants: MorganInvariants,
}

impl Default for MorganFingerprintParams {
    fn default() -> Self {
        // Reuse the canonical source-shaped generator defaults. The per-call
        // values match MorganCall::default while retaining owned vectors here.
        Self {
            generator: MorganParams::default(),
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            conformer_id: -1,
            invariants: MorganInvariants::Connectivity,
        }
    }
}

impl Molecule {
    /// Computes a sparse Morgan count fingerprint with source defaults.
    pub fn morgan_sparse_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint, MorganReadError> {
        self.morgan_sparse_count_fingerprint_with_params(&MorganFingerprintParams::default(), None)
    }

    /// Computes a sparse Morgan count fingerprint using explicit parameters.
    pub fn morgan_sparse_count_fingerprint_with_params(
        &self,
        params: &MorganFingerprintParams,
        additional_output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint, MorganReadError> {
        let prepared = prepare_morgan_read_input(self).map_err(MorganReadError::Preparation)?;
        let input = prepared.owner_input();
        let call = MorganCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
        };
        let invariants = match &params.invariants {
            MorganInvariants::Connectivity => MorganAtomInvariants::Connectivity,
            MorganInvariants::Features => MorganAtomInvariants::Features,
            MorganInvariants::FeaturePatterns(patterns) => {
                MorganAtomInvariants::FeaturePatterns(patterns)
            }
        };

        morgan_sparse_count(
            &input,
            &params.generator,
            &call,
            invariants,
            additional_output,
        )
        .map_err(MorganReadError::Generator)
    }

    /// Computes a sparse Morgan presence fingerprint with source defaults.
    pub fn morgan_sparse_fingerprint(&self) -> Result<SparseBitFingerprint, MorganReadError> {
        self.morgan_sparse_fingerprint_with_params(&MorganFingerprintParams::default(), None)
    }

    /// Computes a sparse Morgan presence fingerprint using explicit parameters.
    pub fn morgan_sparse_fingerprint_with_params(
        &self,
        params: &MorganFingerprintParams,
        additional_output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseBitFingerprint, MorganReadError> {
        let prepared = prepare_morgan_read_input(self).map_err(MorganReadError::Preparation)?;
        let input = prepared.owner_input();
        let call = MorganCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
        };
        let invariants = match &params.invariants {
            MorganInvariants::Connectivity => MorganAtomInvariants::Connectivity,
            MorganInvariants::Features => MorganAtomInvariants::Features,
            MorganInvariants::FeaturePatterns(patterns) => {
                MorganAtomInvariants::FeaturePatterns(patterns)
            }
        };

        morgan_sparse_bits(
            &input,
            &params.generator,
            &call,
            invariants,
            additional_output,
        )
        .map_err(MorganReadError::Generator)
    }

    /// Computes a hashed Morgan count fingerprint with source defaults.
    pub fn morgan_count_fingerprint(&self) -> Result<SparseCountFingerprint32, MorganReadError> {
        self.morgan_count_fingerprint_with_params(&MorganFingerprintParams::default(), None)
    }

    /// Computes a hashed Morgan count fingerprint using explicit parameters.
    pub fn morgan_count_fingerprint_with_params(
        &self,
        params: &MorganFingerprintParams,
        additional_output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, MorganReadError> {
        let prepared = prepare_morgan_read_input(self).map_err(MorganReadError::Preparation)?;
        let input = prepared.owner_input();
        let call = MorganCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
        };
        let invariants = match &params.invariants {
            MorganInvariants::Connectivity => MorganAtomInvariants::Connectivity,
            MorganInvariants::Features => MorganAtomInvariants::Features,
            MorganInvariants::FeaturePatterns(patterns) => {
                MorganAtomInvariants::FeaturePatterns(patterns)
            }
        };

        morgan_count(
            &input,
            &params.generator,
            &call,
            invariants,
            additional_output,
        )
        .map_err(MorganReadError::Generator)
    }

    /// Computes a dense Morgan bit fingerprint with source defaults.
    pub fn morgan_fingerprint(&self) -> Result<Fingerprint, MorganReadError> {
        self.morgan_fingerprint_with_params(&MorganFingerprintParams::default(), None)
    }

    /// Computes a dense Morgan bit fingerprint using explicit parameters.
    pub fn morgan_fingerprint_with_params(
        &self,
        params: &MorganFingerprintParams,
        additional_output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<Fingerprint, MorganReadError> {
        let prepared = prepare_morgan_read_input(self).map_err(MorganReadError::Preparation)?;
        let input = prepared.owner_input();
        let call = MorganCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
        };
        let invariants = match &params.invariants {
            MorganInvariants::Connectivity => MorganAtomInvariants::Connectivity,
            MorganInvariants::Features => MorganAtomInvariants::Features,
            MorganInvariants::FeaturePatterns(patterns) => {
                MorganAtomInvariants::FeaturePatterns(patterns)
            }
        };

        morgan_bits(
            &input,
            &params.generator,
            &call,
            invariants,
            additional_output,
        )
        .map_err(MorganReadError::Generator)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use std::sync::Arc;

    use cosmolkit_core::{RingFindType, RingFindingError, RingInfo, ValenceAssignment};
    use cosmolkit_fingerprints::{FingerprintError, MorganError};

    use crate::molecule::DerivedCacheBlock;
    use crate::ops::OperationError;
    use crate::{
        BINDING_CONTRACT, BindingItem, BindingKind, BindingOwner, BindingReceiver, BindingTypeRole,
        DerivedState, Fingerprint, FingerprintAdditionalOutput, FingerprintPreparationError,
        FunctionStatus, Molecule, MorganFingerprintParams, MorganInvariants, MorganParams,
        MorganReadError, QueryGraph, SparseBitFingerprint, StateModel,
    };

    use super::{MorganReadRings, prepare_morgan_read_input};

    fn empty_molecule_with_cache(cache: DerivedCacheBlock) -> Result<Molecule, OperationError> {
        let empty = Molecule::new();
        Molecule::from_runtime_parts(
            empty.topology_arc_runtime(),
            empty.coordinates_arc_runtime(),
            empty.properties_arc_runtime(),
            Arc::new(cache),
        )
    }

    fn valid_empty_valence() -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        }
    }

    fn public_binding_entry(semantic_id: &str) -> &'static crate::BindingContractEntry {
        BINDING_CONTRACT
            .iter()
            .find(|entry| entry.semantic_id == semantic_id)
            .unwrap_or_else(|| panic!("missing public binding entry {semantic_id}"))
    }

    #[test]
    fn morgan_public_preparation_uses_temporary_symmetrized_rings_without_cache_write() {
        let mut cache = DerivedCacheBlock::default();
        cache.install_valence_assignment(valid_empty_valence());
        cache.mark_valid(DerivedState::VALENCE);
        let molecule = empty_molecule_with_cache(cache)
            .expect("the fixed empty topology accepts empty valid valence rows");
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();
        let before_validity = cache.valid_states();

        let prepared = prepare_morgan_read_input(&molecule)
            .expect("valid valence permits temporary canonical ring preparation");
        let MorganReadRings::Temporary(rings) = &prepared.rings else {
            panic!("absent valid ring cache selects one temporary assignment");
        };
        assert!(rings.is_initialized());
        assert_eq!(rings.find_type(), RingFindType::SymmSssr);
        assert!(rings.atom_rings().is_empty());
        assert!(rings.bond_rings().is_empty());
        assert_eq!(rings.num_rings(), 0);

        let owner_input = prepared.owner_input();
        assert!(std::ptr::eq(owner_input.topology, molecule.topology()));
        assert!(std::ptr::eq(
            owner_input.coordinates,
            molecule.coordinate_block_runtime()
        ));
        assert!(std::ptr::eq(owner_input.properties, molecule.properties()));
        assert!(std::ptr::eq(
            owner_input.valence,
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .expect("the runtime-valid empty valence rows remain borrowed")
        ));
        assert_eq!(owner_input.valence, &valid_empty_valence());
        assert_eq!(owner_input.rings.find_type(), RingFindType::SymmSssr);

        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert!(Arc::ptr_eq(
            &coordinates,
            &molecule.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
        assert_eq!(
            molecule.derived_cache_runtime().valid_states(),
            before_validity
        );
        assert_eq!(before_validity, DerivedState::VALENCE);
        assert_eq!(molecule.derived_cache_runtime().ring_info(), None);
    }

    #[test]
    fn morgan_public_preparation_borrows_valid_cached_ring_rows_unchanged() {
        let cached_rings = RingInfo::new(RingFindType::Fast, 0, 0);
        let mut cache = DerivedCacheBlock::default();
        cache.install_valence_assignment(valid_empty_valence());
        cache.mark_valid(DerivedState::VALENCE);
        cache.install_ring_info(cached_rings);
        cache.mark_valid(DerivedState::RINGS);
        let molecule = empty_molecule_with_cache(cache)
            .expect("the fixed empty topology accepts initialized empty cached rows");
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();
        let before_validity = cache.valid_states();
        let cached_rings = cache.ring_info().expect("the ring payload is present");

        let prepared = prepare_morgan_read_input(&molecule)
            .expect("valid valence and initialized valid ring rows are borrowable");
        let MorganReadRings::Cached(rings) = &prepared.rings else {
            panic!("valid ring cache rows must not be recomputed");
        };
        assert!(std::ptr::eq(*rings, cached_rings));
        assert_eq!(rings.find_type(), RingFindType::Fast);
        assert!(rings.atom_rings().is_empty());
        assert!(rings.bond_rings().is_empty());
        assert_eq!(rings.num_rings(), 0);

        let owner_input = prepared.owner_input();
        assert!(std::ptr::eq(owner_input.rings, cached_rings));
        assert!(std::ptr::eq(owner_input.topology, molecule.topology()));
        assert!(std::ptr::eq(
            owner_input.coordinates,
            molecule.coordinate_block_runtime()
        ));
        assert!(std::ptr::eq(owner_input.properties, molecule.properties()));
        assert!(std::ptr::eq(
            owner_input.valence,
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .expect("valid valence remains borrowed")
        ));
        assert_eq!(owner_input.rings.find_type(), RingFindType::Fast);

        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert!(Arc::ptr_eq(
            &coordinates,
            &molecule.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
        assert_eq!(
            molecule.derived_cache_runtime().valid_states(),
            before_validity
        );
        assert_eq!(
            before_validity,
            DerivedState::VALENCE.union(DerivedState::RINGS)
        );
        assert!(std::ptr::eq(
            molecule.derived_cache_runtime().ring_info().unwrap(),
            cached_rings
        ));
    }

    #[test]
    fn morgan_public_preparation_reports_missing_valid_valence_without_mutation() {
        let molecule = Molecule::new();
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();

        assert!(matches!(
            prepare_morgan_read_input(&molecule),
            Err(FingerprintPreparationError::MissingPreparedValence)
        ));
        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert!(Arc::ptr_eq(
            &coordinates,
            &molecule.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
        assert_eq!(
            molecule.derived_cache_runtime().valid_states(),
            DerivedState::NONE
        );
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_none()
        );
        assert!(molecule.derived_cache_runtime().ring_info().is_none());
    }

    #[test]
    fn morgan_public_preparation_rejects_each_malformed_valence_cache_shape() {
        let cases = [
            (
                "validity_bit",
                Some(ValenceAssignment {
                    explicit_valence: Vec::new(),
                    implicit_hydrogens: Vec::new(),
                }),
                false,
                "validity_bit",
                0,
                1,
            ),
            (
                "explicit_valence",
                Some(ValenceAssignment {
                    explicit_valence: vec![0],
                    implicit_hydrogens: Vec::new(),
                }),
                true,
                "explicit_valence",
                1,
                0,
            ),
            (
                "implicit_hydrogens",
                Some(ValenceAssignment {
                    explicit_valence: Vec::new(),
                    implicit_hydrogens: vec![0],
                }),
                true,
                "implicit_hydrogens",
                1,
                0,
            ),
            ("assignment", None, true, "assignment", 0, 1),
        ];

        for (case, assignment, mark_valid, expected_field, expected_actual, expected_length) in
            cases
        {
            let mut cache = DerivedCacheBlock::default();
            if let Some(assignment) = assignment {
                cache.install_valence_assignment(assignment);
            }
            if mark_valid {
                cache.mark_valid(DerivedState::VALENCE);
            }

            let error = empty_molecule_with_cache(cache)
                .expect_err("runtime construction rejects malformed valence cache state");
            match error {
                OperationError::InvalidDerivedCache {
                    state,
                    field,
                    actual,
                    expected,
                } => {
                    assert_eq!(state, "valence", "case={case}");
                    assert_eq!(field, expected_field, "case={case}");
                    assert_eq!(actual, expected_actual, "case={case}");
                    assert_eq!(expected, expected_length, "case={case}");
                }
                other => panic!("case={case}: unexpected cache construction error {other:?}"),
            }
        }
    }

    #[test]
    fn morgan_public_configuration_defaults_and_empty_option_vectors() {
        let params = MorganFingerprintParams::default();

        assert_eq!(params.generator.radius, 3);
        assert!(!params.generator.include_chirality);
        assert!(params.generator.use_bond_types);
        assert!(params.generator.include_ring_membership);
        assert!(!params.generator.only_nonzero_invariants);
        assert!(!params.generator.include_redundant_environments);
        assert_eq!(params.generator.fp_size, 2048);
        assert!(!params.generator.count_simulation);
        assert_eq!(params.generator.count_bounds, vec![1, 2, 4, 8]);
        assert_eq!(params.generator.bits_per_feature, 1);
        assert!(params.from_atoms.is_none());
        assert!(params.ignore_atoms.is_none());
        assert!(params.custom_atom_invariants.is_none());
        assert!(params.custom_bond_invariants.is_none());
        assert_eq!(params.conformer_id, -1);
        assert!(matches!(params.invariants, MorganInvariants::Connectivity));
        assert_eq!(params.generator, MorganParams::default());

        let empty_vectors = MorganFingerprintParams {
            from_atoms: Some(Vec::new()),
            ignore_atoms: Some(Vec::new()),
            custom_atom_invariants: Some(Vec::new()),
            custom_bond_invariants: Some(Vec::new()),
            ..MorganFingerprintParams::default()
        };
        assert_eq!(empty_vectors.from_atoms.as_deref(), Some(&[][..]));
        assert_eq!(empty_vectors.ignore_atoms.as_deref(), Some(&[][..]));
        assert_eq!(
            empty_vectors.custom_atom_invariants.as_deref(),
            Some(&[][..])
        );
        assert_eq!(
            empty_vectors.custom_bond_invariants.as_deref(),
            Some(&[][..])
        );

        let features = MorganInvariants::Features;
        let empty_patterns = MorganInvariants::FeaturePatterns(Vec::new());
        assert_ne!(features, empty_patterns);
        assert!(matches!(features, MorganInvariants::Features));
        assert!(matches!(
            empty_patterns,
            MorganInvariants::FeaturePatterns(ref patterns) if patterns.is_empty()
        ));
    }

    #[test]
    fn morgan_public_configuration_preserves_existing_query_graph_identity() {
        let graph = QueryGraph::from_parts(
            Vec::new(),
            Vec::new(),
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("the empty query graph is structurally valid");
        let patterns = vec![graph];
        let original_graphs = patterns.as_ptr();

        let params = MorganFingerprintParams {
            invariants: MorganInvariants::FeaturePatterns(patterns),
            ..MorganFingerprintParams::default()
        };

        let MorganInvariants::FeaturePatterns(retained) = &params.invariants else {
            panic!("feature-pattern selection remains explicit");
        };
        assert_eq!(retained.len(), 1);
        assert_eq!(retained.as_ptr(), original_graphs);
    }

    #[test]
    fn morgan_public_error_preserves_typed_causes_and_missing_leaf() {
        let ring_source = RingFindingError::Value {
            message: "fixed ring preparation failure",
        };
        let ring_errors: [Box<dyn std::error::Error>; 3] = [
            Box::new(MorganReadError::Preparation(
                FingerprintPreparationError::RingPreparation(ring_source.clone()),
            )),
            Box::new(crate::AtomPairReadError::Preparation(
                FingerprintPreparationError::RingPreparation(ring_source.clone()),
            )),
            Box::new(crate::TopologicalTorsionReadError::Preparation(
                FingerprintPreparationError::RingPreparation(ring_source.clone()),
            )),
        ];
        for error in &ring_errors {
            let preparation = error
                .source()
                .unwrap()
                .downcast_ref::<FingerprintPreparationError>()
                .expect("each family retains the shared preparation error");
            let ring_cause = std::error::Error::source(preparation)
                .expect("ring preparation retains its concrete cause");
            assert_eq!(
                ring_cause.downcast_ref::<RingFindingError>(),
                Some(&ring_source)
            );
        }

        let generator_error = MorganReadError::Generator(MorganError::Fingerprint(
            FingerprintError::InvalidArguments {
                reason: "fixed Morgan generation failure",
            },
        ));
        let generator_cause = std::error::Error::source(&generator_error)
            .expect("Morgan generation retains the canonical owner error");
        let owner_error = generator_cause
            .downcast_ref::<MorganError>()
            .expect("the first source remains MorganError");
        assert!(matches!(
            owner_error,
            MorganError::Fingerprint(FingerprintError::InvalidArguments {
                reason: "fixed Morgan generation failure"
            })
        ));
        let fingerprint_cause = std::error::Error::source(owner_error)
            .expect("MorganError retains its concrete fingerprint cause");
        assert!(matches!(
            fingerprint_cause.downcast_ref::<FingerprintError>(),
            Some(FingerprintError::InvalidArguments {
                reason: "fixed Morgan generation failure"
            })
        ));

        let molecule = Molecule::new();
        let errors: [Box<dyn std::error::Error>; 12] = [
            Box::new(molecule.morgan_fingerprint().unwrap_err()),
            Box::new(molecule.morgan_sparse_fingerprint().unwrap_err()),
            Box::new(molecule.morgan_count_fingerprint().unwrap_err()),
            Box::new(molecule.morgan_sparse_count_fingerprint().unwrap_err()),
            Box::new(molecule.atom_pair_fingerprint().unwrap_err()),
            Box::new(molecule.atom_pair_sparse_fingerprint().unwrap_err()),
            Box::new(molecule.atom_pair_count_fingerprint().unwrap_err()),
            Box::new(molecule.atom_pair_sparse_count_fingerprint().unwrap_err()),
            Box::new(molecule.topological_torsion_fingerprint().unwrap_err()),
            Box::new(
                molecule
                    .topological_torsion_sparse_fingerprint()
                    .unwrap_err(),
            ),
            Box::new(
                molecule
                    .topological_torsion_count_fingerprint()
                    .unwrap_err(),
            ),
            Box::new(
                molecule
                    .topological_torsion_sparse_count_fingerprint()
                    .unwrap_err(),
            ),
        ];
        for error in &errors {
            let preparation = error
                .source()
                .unwrap()
                .downcast_ref::<FingerprintPreparationError>()
                .expect("public read errors retain a shared preparation cause");
            assert!(matches!(
                preparation,
                FingerprintPreparationError::MissingPreparedValence
            ));
            assert!(std::error::Error::source(preparation).is_none());
        }
    }

    #[test]
    fn morgan_public_metadata_contract_additional_output_bindings() {
        for (semantic_id, rust_path, python_name, javascript_name, role) in [
            (
                "types.FingerprintPreparationError",
                "crate::FingerprintPreparationError",
                "FingerprintPreparationError",
                "FingerprintPreparationError",
                BindingTypeRole::Error,
            ),
            (
                "types.MorganParams",
                "crate::MorganParams",
                "MorganParams",
                "MorganParams",
                BindingTypeRole::Parameter,
            ),
            (
                "types.MorganInvariants",
                "crate::MorganInvariants",
                "MorganInvariants",
                "MorganInvariants",
                BindingTypeRole::Parameter,
            ),
            (
                "types.MorganFingerprintParams",
                "crate::MorganFingerprintParams",
                "MorganFingerprintParams",
                "MorganFingerprintParams",
                BindingTypeRole::Parameter,
            ),
            (
                "types.MorganReadError",
                "crate::MorganReadError",
                "MorganReadError",
                "MorganReadError",
                BindingTypeRole::Error,
            ),
            (
                "types.FingerprintAdditionalOutput",
                "crate::FingerprintAdditionalOutput",
                "FingerprintAdditionalOutput",
                "FingerprintAdditionalOutput",
                BindingTypeRole::Value,
            ),
        ] {
            let entry = public_binding_entry(semantic_id);
            assert_eq!(entry.item, BindingItem::Type, "{semantic_id}");
            assert_eq!(entry.owner, BindingOwner::Type, "{semantic_id}");
            assert_eq!(entry.rust_path.replace(' ', ""), rust_path, "{semantic_id}");
            assert_eq!(entry.python_name, python_name, "{semantic_id}");
            assert_eq!(entry.javascript_name, javascript_name, "{semantic_id}");
            assert_eq!(entry.feature, "cap-fingerprints", "{semantic_id}");
            assert_eq!(entry.status, FunctionStatus::Experimental, "{semantic_id}");
            assert_eq!(entry.type_role, Some(role), "{semantic_id}");
            assert!(entry.callable.is_none(), "{semantic_id}");
        }

        let default_entry = public_binding_entry("FingerprintAdditionalOutput.default");
        assert_eq!(default_entry.item, BindingItem::Callable);
        assert_eq!(default_entry.owner, BindingOwner::Type);
        assert_eq!(
            default_entry.rust_path.replace(' ', ""),
            "crate::FingerprintAdditionalOutput::default"
        );
        assert_eq!(default_entry.python_name, "default");
        assert_eq!(default_entry.javascript_name, "default");
        assert_eq!(default_entry.feature, "cap-fingerprints");
        assert_eq!(default_entry.status, FunctionStatus::Experimental);
        let default = default_entry.callable.unwrap();
        assert_eq!(default.kind, BindingKind::Static);
        assert_eq!(default.receiver, None);
        assert!(default.parameters.is_empty());
        assert_eq!(
            default.output_type.replace(' ', ""),
            "crate::FingerprintAdditionalOutput"
        );
        let constructor_entry = public_binding_entry("FingerprintAdditionalOutput.new");
        assert_eq!(
            constructor_entry.rust_path.replace(' ', ""),
            "crate::FingerprintAdditionalOutput::new"
        );
        let constructor = constructor_entry.callable.unwrap();
        assert_eq!(
            constructor.output_type.replace(' ', ""),
            "crate::FingerprintAdditionalOutput"
        );
        assert_eq!(constructor.error_type, None);
        assert!(constructor.parameters.is_empty());
        assert_eq!(default.error_type, None);
        assert_eq!(default.state_model, StateModel::ValueReturning);
        assert_eq!(default.operation_semantic_id, None);

        let allocator_and_getter_contracts = [
            (
                "FingerprintAdditionalOutput.allocate_atom_counts",
                "crate::FingerprintAdditionalOutput::allocate_atom_counts",
                "allocate_atom_counts",
                "allocateAtomCounts",
                BindingReceiver::Mutable,
                "()",
                StateModel::InPlace,
            ),
            (
                "FingerprintAdditionalOutput.allocate_atom_to_bits",
                "crate::FingerprintAdditionalOutput::allocate_atom_to_bits",
                "allocate_atom_to_bits",
                "allocateAtomToBits",
                BindingReceiver::Mutable,
                "()",
                StateModel::InPlace,
            ),
            (
                "FingerprintAdditionalOutput.allocate_bit_info_map",
                "crate::FingerprintAdditionalOutput::allocate_bit_info_map",
                "allocate_bit_info_map",
                "allocateBitInfoMap",
                BindingReceiver::Mutable,
                "()",
                StateModel::InPlace,
            ),
            (
                "FingerprintAdditionalOutput.allocate_bit_paths",
                "crate::FingerprintAdditionalOutput::allocate_bit_paths",
                "allocate_bit_paths",
                "allocateBitPaths",
                BindingReceiver::Mutable,
                "()",
                StateModel::InPlace,
            ),
            (
                "FingerprintAdditionalOutput.allocate_atoms_per_bit",
                "crate::FingerprintAdditionalOutput::allocate_atoms_per_bit",
                "allocate_atoms_per_bit",
                "allocateAtomsPerBit",
                BindingReceiver::Mutable,
                "()",
                StateModel::InPlace,
            ),
            (
                "FingerprintAdditionalOutput.atom_counts",
                "crate::FingerprintAdditionalOutput::atom_counts",
                "atom_counts",
                "atomCounts",
                BindingReceiver::Shared,
                "Option<&[u32]>",
                StateModel::ReadOnly,
            ),
            (
                "FingerprintAdditionalOutput.atom_to_bits",
                "crate::FingerprintAdditionalOutput::atom_to_bits",
                "atom_to_bits",
                "atomToBits",
                BindingReceiver::Shared,
                "Option<&[Vec<u64>]>",
                StateModel::ReadOnly,
            ),
            (
                "FingerprintAdditionalOutput.bit_info_map",
                "crate::FingerprintAdditionalOutput::bit_info_map",
                "bit_info_map",
                "bitInfoMap",
                BindingReceiver::Shared,
                "Option<&std::collections::BTreeMap<u64,Vec<(u32,u32)>>>",
                StateModel::ReadOnly,
            ),
            (
                "FingerprintAdditionalOutput.bit_paths",
                "crate::FingerprintAdditionalOutput::bit_paths",
                "bit_paths",
                "bitPaths",
                BindingReceiver::Shared,
                "Option<&std::collections::BTreeMap<u64,Vec<Vec<i32>>>>",
                StateModel::ReadOnly,
            ),
            (
                "FingerprintAdditionalOutput.atoms_per_bit",
                "crate::FingerprintAdditionalOutput::atoms_per_bit",
                "atoms_per_bit",
                "atomsPerBit",
                BindingReceiver::Shared,
                "Option<&std::collections::BTreeMap<u64,Vec<Vec<i32>>>>",
                StateModel::ReadOnly,
            ),
        ];

        for (semantic_id, rust_path, python_name, javascript_name, receiver, output, state) in
            allocator_and_getter_contracts
        {
            let entry = public_binding_entry(semantic_id);
            assert_eq!(entry.item, BindingItem::Callable, "{semantic_id}");
            assert_eq!(entry.owner, BindingOwner::Type, "{semantic_id}");
            assert_eq!(entry.rust_path.replace(' ', ""), rust_path, "{semantic_id}");
            assert_eq!(entry.python_name, python_name, "{semantic_id}");
            assert_eq!(entry.javascript_name, javascript_name, "{semantic_id}");
            assert_eq!(entry.feature, "cap-fingerprints", "{semantic_id}");
            assert_eq!(entry.status, FunctionStatus::Experimental, "{semantic_id}");
            assert!(entry.type_role.is_none(), "{semantic_id}");
            let callable = entry.callable.unwrap();
            assert_eq!(callable.kind, BindingKind::Instance, "{semantic_id}");
            assert_eq!(callable.receiver, Some(receiver), "{semantic_id}");
            assert!(callable.parameters.is_empty(), "{semantic_id}");
            assert_eq!(
                callable.output_type.replace(' ', ""),
                output,
                "{semantic_id}"
            );
            assert_eq!(callable.error_type, None, "{semantic_id}");
            assert_eq!(callable.state_model, state, "{semantic_id}");
            assert_eq!(callable.operation_semantic_id, None, "{semantic_id}");
        }

        let _: fn() -> FingerprintAdditionalOutput = FingerprintAdditionalOutput::default;
        let _: fn(&mut FingerprintAdditionalOutput) =
            FingerprintAdditionalOutput::allocate_atom_counts;
        let _: fn(&mut FingerprintAdditionalOutput) =
            FingerprintAdditionalOutput::allocate_atom_to_bits;
        let _: fn(&mut FingerprintAdditionalOutput) =
            FingerprintAdditionalOutput::allocate_bit_info_map;
        let _: fn(&mut FingerprintAdditionalOutput) =
            FingerprintAdditionalOutput::allocate_bit_paths;
        let _: fn(&mut FingerprintAdditionalOutput) =
            FingerprintAdditionalOutput::allocate_atoms_per_bit;
        let _: for<'a> fn(&'a FingerprintAdditionalOutput) -> Option<&'a [u32]> =
            FingerprintAdditionalOutput::atom_counts;
        let _: for<'a> fn(&'a FingerprintAdditionalOutput) -> Option<&'a [Vec<u64>]> =
            FingerprintAdditionalOutput::atom_to_bits;
        let _: for<'a> fn(
            &'a FingerprintAdditionalOutput,
        ) -> Option<&'a BTreeMap<u64, Vec<(u32, u32)>>> = FingerprintAdditionalOutput::bit_info_map;
        let _: for<'a> fn(
            &'a FingerprintAdditionalOutput,
        ) -> Option<&'a BTreeMap<u64, Vec<Vec<i32>>>> = FingerprintAdditionalOutput::bit_paths;
        let _: for<'a> fn(
            &'a FingerprintAdditionalOutput,
        ) -> Option<&'a BTreeMap<u64, Vec<Vec<i32>>>> = FingerprintAdditionalOutput::atoms_per_bit;

        let default = FingerprintAdditionalOutput::default();
        assert!(default.atom_counts().is_none());
        assert!(default.atom_to_bits().is_none());
        assert!(default.bit_info_map().is_none());
        assert!(default.bit_paths().is_none());
        assert!(default.atoms_per_bit().is_none());
    }

    #[test]
    fn morgan_public_metadata_contract_all_additional_output_allocation_masks() {
        let mut masks_checked = 0;
        for mask in 0_u8..32 {
            let allocated = [
                mask & 0b00001 != 0,
                mask & 0b00010 != 0,
                mask & 0b00100 != 0,
                mask & 0b01000 != 0,
                mask & 0b10000 != 0,
            ];
            let mut output = FingerprintAdditionalOutput::default();
            assert!(
                output.atom_counts().is_none(),
                "default counts, mask {mask:#07b}"
            );
            assert!(
                output.atom_to_bits().is_none(),
                "default atomToBits, mask {mask:#07b}"
            );
            assert!(
                output.bit_info_map().is_none(),
                "default bitInfoMap, mask {mask:#07b}"
            );
            assert!(
                output.bit_paths().is_none(),
                "default bitPaths, mask {mask:#07b}"
            );
            assert!(
                output.atoms_per_bit().is_none(),
                "default atomsPerBit, mask {mask:#07b}"
            );

            if allocated[0] {
                output.allocate_atom_counts();
            }
            if allocated[1] {
                output.allocate_atom_to_bits();
            }
            if allocated[2] {
                output.allocate_bit_info_map();
            }
            if allocated[3] {
                output.allocate_bit_paths();
            }
            if allocated[4] {
                output.allocate_atoms_per_bit();
            }

            assert_eq!(
                output.atom_counts().map(|values| values.len()),
                allocated[0].then_some(0),
                "atom_counts None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.atom_to_bits().map(|rows| rows.len()),
                allocated[1].then_some(0),
                "atom_to_bits None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.bit_info_map().map(BTreeMap::len),
                allocated[2].then_some(0),
                "bit_info_map None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.bit_paths().map(BTreeMap::len),
                allocated[3].then_some(0),
                "bit_paths None/Some(empty), mask {mask:#07b}"
            );
            assert_eq!(
                output.atoms_per_bit().map(BTreeMap::len),
                allocated[4].then_some(0),
                "atoms_per_bit None/Some(empty), mask {mask:#07b}"
            );
            masks_checked += 1;
        }
        assert_eq!(masks_checked, 32, "all five-bit public allocation masks");
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    fn morgan_public_sparse_count_cases() -> [(&'static str, &'static [(u64, i32)]); 5] {
        [
            ("", &[]),
            (
                "CCO",
                &[
                    (864_662_311, 1),
                    (1_535_166_686, 1),
                    (2_245_384_272, 1),
                    (2_246_728_737, 1),
                    (3_542_456_614, 1),
                    (4_018_048_386, 1),
                ],
            ),
            (
                "c1cc[nH]c1",
                &[
                    (98_513_984, 2),
                    (116_898_731, 2),
                    (1_482_649_460, 2),
                    (2_132_511_834, 1),
                    (2_266_426_494, 1),
                    (2_293_755_984, 1),
                    (2_654_043_257, 1),
                    (2_753_863_138, 2),
                    (3_218_693_969, 4),
                ],
            ),
            (
                "c1ccccc1",
                &[
                    (98_513_984, 6),
                    (2_763_854_213, 6),
                    (3_218_693_969, 6),
                    (3_741_631_696, 1),
                ],
            ),
            ("[Na+].[Cl-]", &[(3_737_048_253, 1), (3_855_292_234, 1)]),
        ]
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    fn morgan_public_sparse_count_molecule(smiles: &str) -> Molecule {
        Molecule::from_smiles(smiles)
            .expect("the fixed source-literal SMILES parses")
            .with_assigned_valence()
            .expect("the fixed source-literal molecule has prepared valence")
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    fn morgan_public_sparse_count_output_mask(mask: u8) -> FingerprintAdditionalOutput {
        let mut output = FingerprintAdditionalOutput::default();
        if mask & 0b00001 != 0 {
            output.allocate_atom_counts();
        }
        if mask & 0b00010 != 0 {
            output.allocate_atom_to_bits();
        }
        if mask & 0b00100 != 0 {
            output.allocate_bit_info_map();
        }
        if mask & 0b01000 != 0 {
            output.allocate_bit_paths();
        }
        if mask & 0b10000 != 0 {
            output.allocate_atoms_per_bit();
        }
        output
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    fn assert_morgan_public_sparse_count_output_mask(
        actual: &FingerprintAdditionalOutput,
        complete: &FingerprintAdditionalOutput,
        mask: u8,
        smiles: &str,
    ) {
        assert_eq!(
            actual.atom_counts(),
            (mask & 0b00001 != 0).then(|| complete.atom_counts().unwrap()),
            "atom counts, mask {mask:#07b}, {smiles}"
        );
        assert_eq!(
            actual.atom_to_bits(),
            (mask & 0b00010 != 0).then(|| complete.atom_to_bits().unwrap()),
            "atom-to-bits, mask {mask:#07b}, {smiles}"
        );
        assert_eq!(
            actual.bit_info_map(),
            (mask & 0b00100 != 0).then(|| complete.bit_info_map().unwrap()),
            "bit-info map, mask {mask:#07b}, {smiles}"
        );
        assert_eq!(
            actual.bit_paths(),
            (mask & 0b01000 != 0).then(|| complete.bit_paths().unwrap()),
            "bit paths, mask {mask:#07b}, {smiles}"
        );
        assert_eq!(
            actual.atoms_per_bit(),
            (mask & 0b10000 != 0).then(|| complete.atoms_per_bit().unwrap()),
            "atoms per bit, mask {mask:#07b}, {smiles}"
        );
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_sparse_count_fixed_literals_and_short_defaults() {
        for (smiles, expected_entries) in morgan_public_sparse_count_cases() {
            let molecule = morgan_public_sparse_count_molecule(smiles);
            let expected: BTreeMap<u64, i32> = expected_entries.iter().copied().collect();
            let params = MorganFingerprintParams::default();
            let configured = molecule
                .morgan_sparse_count_fingerprint_with_params(&params, None)
                .expect("the fixed source-literal sparse-count call succeeds");
            let short = molecule
                .morgan_sparse_count_fingerprint()
                .expect("the short source-default sparse-count call succeeds");

            assert_eq!(configured.nonzero_elements(), &expected, "{smiles}");
            assert_eq!(short.nonzero_elements(), &expected, "{smiles}");
            assert_eq!(
                short.nonzero_elements(),
                configured.nonzero_elements(),
                "{smiles}"
            );
            assert_eq!(short.length(), configured.length(), "{smiles}");
        }
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_sparse_count_all_output_masks_none_reuse_and_preservation() {
        for (smiles, expected_entries) in morgan_public_sparse_count_cases() {
            let molecule = morgan_public_sparse_count_molecule(smiles);
            let topology = molecule.topology_arc_runtime();
            let coordinates = molecule.coordinates_arc_runtime();
            let properties = molecule.properties_arc_runtime();
            let cache = molecule.derived_cache_arc_runtime();
            let valid_states = molecule.derived_cache_runtime().valid_states();
            let params = MorganFingerprintParams::default();
            let expected: BTreeMap<u64, i32> = expected_entries.iter().copied().collect();

            let mut complete_output = morgan_public_sparse_count_output_mask(0b1_1111);
            let complete = molecule
                .morgan_sparse_count_fingerprint_with_params(&params, Some(&mut complete_output))
                .expect("the all-channel source-literal call succeeds");
            assert_eq!(complete.nonzero_elements(), &expected, "{smiles}");

            let without_output = molecule
                .morgan_sparse_count_fingerprint_with_params(&params, None)
                .expect("a null FingerprintAdditionalOutput remains valid");
            assert_eq!(without_output.nonzero_elements(), &expected, "{smiles}");

            let mut masks_checked = 0;
            for mask in 0_u8..32 {
                let mut output = morgan_public_sparse_count_output_mask(mask);
                let actual = molecule
                    .morgan_sparse_count_fingerprint_with_params(&params, Some(&mut output))
                    .expect("each independent output-allocation mask succeeds");
                assert_eq!(
                    actual.nonzero_elements(),
                    &expected,
                    "mask={mask:#07b}, {smiles}"
                );
                assert_morgan_public_sparse_count_output_mask(
                    &output,
                    &complete_output,
                    mask,
                    smiles,
                );
                masks_checked += 1;
            }
            assert_eq!(masks_checked, 32, "all five-bit masks, {smiles}");

            assert!(
                Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()),
                "{smiles}"
            );
            assert!(
                Arc::ptr_eq(&coordinates, &molecule.coordinates_arc_runtime()),
                "{smiles}"
            );
            assert!(
                Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()),
                "{smiles}"
            );
            assert!(
                Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()),
                "{smiles}"
            );
            assert_eq!(
                molecule.derived_cache_runtime().valid_states(),
                valid_states,
                "{smiles}"
            );
        }

        let mut reused = morgan_public_sparse_count_output_mask(0b1_1111);
        let params = MorganFingerprintParams::default();
        let cco = morgan_public_sparse_count_molecule("CCO");
        cco.morgan_sparse_count_fingerprint_with_params(&params, Some(&mut reused))
            .expect("the first reused-output call populates Morgan channels");
        let retained_atoms_per_bit = reused
            .atoms_per_bit()
            .expect("the allocated atoms-per-bit channel remains present")
            .clone();
        assert!(!retained_atoms_per_bit.is_empty());

        let empty = morgan_public_sparse_count_molecule("");
        let reset = empty
            .morgan_sparse_count_fingerprint_with_params(&params, Some(&mut reused))
            .expect("the second reused-output call succeeds for the empty molecule");
        assert!(reset.nonzero_elements().is_empty());
        assert_eq!(reused.atom_counts(), Some(&[][..]));
        assert_eq!(reused.atom_to_bits(), Some(&[][..]));
        assert!(reused.bit_info_map().is_some_and(BTreeMap::is_empty));
        assert!(reused.bit_paths().is_some_and(BTreeMap::is_empty));
        assert_eq!(reused.atoms_per_bit(), Some(&retained_atoms_per_bit));
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_count_fixed_literals_and_short_defaults() {
        let molecule = morgan_public_sparse_count_molecule("CCO");
        let params = MorganFingerprintParams::default();
        let expected_default = BTreeMap::from([
            (80_u32, 1_i32),
            (222, 1),
            (294, 1),
            (807, 1),
            (1_057, 1),
            (1_410, 1),
        ]);

        let configured = molecule
            .morgan_count_fingerprint_with_params(&params, None)
            .expect("the fixed source-default hashed count call succeeds");
        let short = molecule
            .morgan_count_fingerprint()
            .expect("the short source-default hashed count call succeeds");
        assert_eq!(configured.length(), 2048);
        assert_eq!(configured.nonzero_elements(), &expected_default);
        assert_eq!(short.nonzero_elements(), &expected_default);
        assert_eq!(short.nonzero_elements(), configured.nonzero_elements());
        assert_eq!(short.length(), configured.length());

        // W02 fixes all six default-radius CCO environments folding to index
        // zero at size one; O03 separately fixes three custom radius-zero IDs
        // colliding at that same index.
        let size_one_params = MorganFingerprintParams {
            generator: MorganParams {
                fp_size: 1,
                ..MorganParams::default()
            },
            ..MorganFingerprintParams::default()
        };
        let size_one = molecule
            .morgan_count_fingerprint_with_params(&size_one_params, None)
            .expect("the fixed source size-one collision call succeeds");
        assert_eq!(size_one.length(), 1);
        assert_eq!(size_one.nonzero_elements(), &BTreeMap::from([(0, 6)]));

        let custom_collision_params = MorganFingerprintParams {
            generator: MorganParams {
                radius: 0,
                fp_size: 1,
                ..MorganParams::default()
            },
            custom_atom_invariants: Some(vec![1, 2, 3]),
            ..MorganFingerprintParams::default()
        };
        let custom_collision = molecule
            .morgan_count_fingerprint_with_params(&custom_collision_params, None)
            .expect("the fixed source custom-invariant collision call succeeds");
        assert_eq!(custom_collision.length(), 1);
        assert_eq!(
            custom_collision.nonzero_elements(),
            &BTreeMap::from([(0, 3)])
        );
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_count_all_output_masks_none_and_source_order() {
        let molecule = morgan_public_sparse_count_molecule("CCO");
        let params = MorganFingerprintParams::default();
        let expected = BTreeMap::from([
            (80_u32, 1_i32),
            (222, 1),
            (294, 1),
            (807, 1),
            (1_057, 1),
            (1_410, 1),
        ]);

        let mut complete_output = morgan_public_sparse_count_output_mask(0b1_1111);
        let complete = molecule
            .morgan_count_fingerprint_with_params(&params, Some(&mut complete_output))
            .expect("the all-channel fixed source count call succeeds");
        assert_eq!(complete.length(), 2048);
        assert_eq!(complete.nonzero_elements(), &expected);

        let without_output = molecule
            .morgan_count_fingerprint_with_params(&params, None)
            .expect("a null FingerprintAdditionalOutput remains valid for hashed counts");
        assert_eq!(without_output.length(), 2048);
        assert_eq!(without_output.nonzero_elements(), &expected);

        let mut masks_checked = 0;
        for mask in 0_u8..32 {
            let mut output = morgan_public_sparse_count_output_mask(mask);
            let actual = molecule
                .morgan_count_fingerprint_with_params(&params, Some(&mut output))
                .expect("each independent count-output allocation mask succeeds");
            assert_eq!(actual.length(), 2048, "mask={mask:#07b}");
            assert_eq!(actual.nonzero_elements(), &expected, "mask={mask:#07b}");
            // Compare borrowed nested channels directly: their vector order is
            // source traversal order and must not be normalized or sorted.
            assert_morgan_public_sparse_count_output_mask(&output, &complete_output, mask, "CCO");
            masks_checked += 1;
        }
        assert_eq!(masks_checked, 32, "all five-bit count-output masks");
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    fn morgan_public_sparse_bits_cases() -> [(&'static str, &'static [i32]); 5] {
        [
            ("", &[]),
            (
                "CCO",
                &[
                    -2_049_583_024,
                    -2_048_238_559,
                    -752_510_682,
                    -276_918_910,
                    864_662_311,
                    1_535_166_686,
                ],
            ),
            (
                "c1cc[nH]c1",
                &[
                    -2_028_540_802,
                    -2_001_211_312,
                    -1_640_924_039,
                    -1_541_104_158,
                    -1_076_273_327,
                    98_513_984,
                    116_898_731,
                    1_482_649_460,
                    2_132_511_834,
                ],
            ),
            (
                "c1ccccc1",
                &[-1_531_113_083, -1_076_273_327, -553_335_600, 98_513_984],
            ),
            ("[Na+].[Cl-]", &[-557_919_043, -439_675_062]),
        ]
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_sparse_bits_fixed_literals_and_short_defaults() {
        for (smiles, expected_bits) in morgan_public_sparse_bits_cases() {
            let molecule = morgan_public_sparse_count_molecule(smiles);
            let params = MorganFingerprintParams::default();
            let configured = molecule
                .morgan_sparse_fingerprint_with_params(&params, None)
                .expect("the fixed source-literal sparse-bit call succeeds");
            let short = molecule
                .morgan_sparse_fingerprint()
                .expect("the short source-default sparse-bit call succeeds");

            for (label, actual) in [("configured", &configured), ("short", &short)] {
                assert_eq!(actual.n_bits(), u32::MAX, "{label}, {smiles}");
                assert_eq!(
                    actual.num_on_bits(),
                    expected_bits.len() as u32,
                    "{label}, {smiles}"
                );
                assert_eq!(
                    actual.on_bits().as_slice(),
                    expected_bits,
                    "{label}, {smiles}"
                );
            }
            assert_eq!(short.on_bits(), configured.on_bits(), "{smiles}");
            assert_eq!(short.n_bits(), configured.n_bits(), "{smiles}");
        }
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_sparse_bits_all_output_masks_none_and_preservation() {
        for (smiles, expected_bits) in morgan_public_sparse_bits_cases() {
            let molecule = morgan_public_sparse_count_molecule(smiles);
            let topology = molecule.topology_arc_runtime();
            let coordinates = molecule.coordinates_arc_runtime();
            let properties = molecule.properties_arc_runtime();
            let cache = molecule.derived_cache_arc_runtime();
            let valid_states = molecule.derived_cache_runtime().valid_states();
            let params = MorganFingerprintParams::default();

            let mut complete_output = morgan_public_sparse_count_output_mask(0b1_1111);
            let complete = molecule
                .morgan_sparse_fingerprint_with_params(&params, Some(&mut complete_output))
                .expect("the all-channel source-literal sparse-bit call succeeds");
            assert_eq!(complete.n_bits(), u32::MAX, "{smiles}");
            assert_eq!(complete.on_bits().as_slice(), expected_bits, "{smiles}");

            let without_output = molecule
                .morgan_sparse_fingerprint_with_params(&params, None)
                .expect("a null FingerprintAdditionalOutput remains valid for sparse bits");
            assert_eq!(without_output.n_bits(), u32::MAX, "{smiles}");
            assert_eq!(
                without_output.on_bits().as_slice(),
                expected_bits,
                "{smiles}"
            );

            let mut masks_checked = 0;
            for mask in 0_u8..32 {
                let mut output = morgan_public_sparse_count_output_mask(mask);
                let actual = molecule
                    .morgan_sparse_fingerprint_with_params(&params, Some(&mut output))
                    .expect("each independent sparse-bit output-allocation mask succeeds");
                assert_eq!(actual.n_bits(), u32::MAX, "mask={mask:#07b}, {smiles}");
                assert_eq!(
                    actual.on_bits().as_slice(),
                    expected_bits,
                    "mask={mask:#07b}, {smiles}"
                );
                assert_morgan_public_sparse_count_output_mask(
                    &output,
                    &complete_output,
                    mask,
                    smiles,
                );
                masks_checked += 1;
            }
            assert_eq!(masks_checked, 32, "all five-bit masks, {smiles}");

            assert!(
                Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()),
                "{smiles}"
            );
            assert!(
                Arc::ptr_eq(&coordinates, &molecule.coordinates_arc_runtime()),
                "{smiles}"
            );
            assert!(
                Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()),
                "{smiles}"
            );
            assert!(
                Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()),
                "{smiles}"
            );
            assert_eq!(
                molecule.derived_cache_runtime().valid_states(),
                valid_states,
                "{smiles}"
            );
        }
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_bits_fixed_literals_and_short_defaults() {
        let molecule = morgan_public_sparse_count_molecule("CCO");
        let params = MorganFingerprintParams::default();
        let expected_default = [80_u32, 222, 294, 807, 1_057, 1_410];

        let configured = molecule
            .morgan_fingerprint_with_params(&params, None)
            .expect("the fixed source-default dense-bit call succeeds");
        let short = molecule
            .morgan_fingerprint()
            .expect("the short source-default dense-bit call succeeds");
        assert_eq!(configured.n_bits(), 2048);
        assert_eq!(configured.on_bits().as_slice(), &expected_default);
        assert_eq!(short.n_bits(), 2048);
        assert_eq!(short.on_bits().as_slice(), &expected_default);
        assert_eq!(short, configured);

        let collision_params = MorganFingerprintParams {
            generator: MorganParams {
                fp_size: 1,
                ..MorganParams::default()
            },
            ..MorganFingerprintParams::default()
        };
        let collision = molecule
            .morgan_fingerprint_with_params(&collision_params, None)
            .expect("the fixed source size-one dense collision call succeeds");
        assert_eq!(collision.n_bits(), 1);
        assert_eq!(collision.on_bits(), vec![0]);

        let count_simulation_params = MorganFingerprintParams {
            generator: MorganParams {
                radius: 0,
                fp_size: 128,
                count_simulation: true,
                count_bounds: vec![1, 2, 4, 8],
                ..MorganParams::default()
            },
            custom_atom_invariants: Some(vec![1, 2, 1]),
            ..MorganFingerprintParams::default()
        };
        let simulated = molecule
            .morgan_fingerprint_with_params(&count_simulation_params, None)
            .expect("the fixed O04 dense count-simulation call succeeds");
        assert_eq!(simulated.n_bits(), 128);
        assert_eq!(simulated.on_bits(), vec![4, 5, 8]);
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_bits_all_output_masks_none_and_exact_lengths() {
        let molecule = morgan_public_sparse_count_molecule("CCO");
        let params = MorganFingerprintParams {
            generator: MorganParams {
                radius: 0,
                fp_size: 128,
                count_simulation: true,
                count_bounds: vec![1, 2, 4, 8],
                ..MorganParams::default()
            },
            custom_atom_invariants: Some(vec![1, 2, 1]),
            ..MorganFingerprintParams::default()
        };
        let expected_bits = [4_u32, 5, 8];

        let mut complete_output = morgan_public_sparse_count_output_mask(0b1_1111);
        let complete = molecule
            .morgan_fingerprint_with_params(&params, Some(&mut complete_output))
            .expect("the all-channel fixed O04 dense call succeeds");
        assert_eq!(complete.n_bits(), 128);
        assert_eq!(complete.on_bits().as_slice(), &expected_bits);

        let expected_atom_counts = [1_u32, 1, 1];
        let expected_atom_to_bits = vec![vec![4_u64, 5], vec![8], vec![4, 5]];
        let expected_bit_info: BTreeMap<u64, Vec<(u32, u32)>> = BTreeMap::from([
            (4, vec![(0, 0), (2, 0)]),
            (5, vec![(0, 0), (2, 0)]),
            (8, vec![(1, 0)]),
        ]);
        let expected_bit_paths: BTreeMap<u64, Vec<Vec<i32>>> = BTreeMap::new();
        let expected_atoms_per_bit: BTreeMap<u64, Vec<Vec<i32>>> = BTreeMap::from([
            (4, vec![vec![0], vec![2]]),
            (5, vec![vec![0], vec![2]]),
            (8, vec![vec![1]]),
        ]);
        assert_eq!(
            complete_output.atom_counts(),
            Some(&expected_atom_counts[..])
        );
        assert_eq!(
            complete_output.atom_to_bits(),
            Some(&expected_atom_to_bits[..])
        );
        assert_eq!(complete_output.bit_info_map(), Some(&expected_bit_info));
        assert_eq!(complete_output.bit_paths(), Some(&expected_bit_paths));
        assert_eq!(
            complete_output.atoms_per_bit(),
            Some(&expected_atoms_per_bit)
        );

        let without_output = molecule
            .morgan_fingerprint_with_params(&params, None)
            .expect("a null FingerprintAdditionalOutput remains valid for dense bits");
        assert_eq!(without_output.n_bits(), 128);
        assert_eq!(without_output.on_bits().as_slice(), &expected_bits);

        let mut masks_checked = 0;
        for mask in 0_u8..32 {
            let mut output = morgan_public_sparse_count_output_mask(mask);
            let actual = molecule
                .morgan_fingerprint_with_params(&params, Some(&mut output))
                .expect("each independent dense output-allocation mask succeeds");
            assert_eq!(actual.n_bits(), 128, "mask={mask:#07b}");
            assert_eq!(
                actual.on_bits().as_slice(),
                &expected_bits,
                "mask={mask:#07b}"
            );
            assert_morgan_public_sparse_count_output_mask(&output, &complete_output, mask, "CCO");
            masks_checked += 1;
        }
        assert_eq!(masks_checked, 32, "all five-bit dense-output masks");
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn morgan_public_call_product_all_configured_families() {
        #[derive(Clone, Copy)]
        struct Expected {
            sparse_counts: &'static [(u64, i32)],
            sparse_bits: &'static [i32],
            counts: &'static [(u32, i32)],
            bits: &'static [u32],
        }

        const DEFAULT_RADIUS_ZERO: Expected = Expected {
            sparse_counts: &[(864_662_311, 1), (2_245_384_272, 1), (2_246_728_737, 1)],
            sparse_bits: &[-2_049_583_024, -2_048_238_559, 864_662_311],
            counts: &[(80, 1), (807, 1), (1_057, 1)],
            bits: &[80, 807, 1_057],
        };
        const EMPTY: Expected = Expected {
            sparse_counts: &[],
            sparse_bits: &[],
            counts: &[],
            bits: &[],
        };
        const CUSTOM_ATOMS: Expected = Expected {
            sparse_counts: &[(1, 2), (2, 1)],
            sparse_bits: &[1, 2],
            counts: &[(1, 2), (2, 1)],
            bits: &[1, 2],
        };
        const FEATURES: Expected = Expected {
            sparse_counts: &[(0, 2), (3, 1)],
            sparse_bits: &[0, 3],
            counts: &[(0, 2), (3, 1)],
            bits: &[0, 3],
        };
        const EMPTY_FEATURE_PATTERNS: Expected = Expected {
            sparse_counts: &[(0, 3)],
            sparse_bits: &[0],
            counts: &[(0, 3)],
            bits: &[0],
        };

        let molecule = morgan_public_sparse_count_molecule("CCO");
        let topology = molecule.topology_arc_runtime();
        let coordinates = molecule.coordinates_arc_runtime();
        let properties = molecule.properties_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();
        let valid_states = molecule.derived_cache_runtime().valid_states();
        let assert_source_unchanged = || {
            assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
            assert!(Arc::ptr_eq(
                &coordinates,
                &molecule.coordinates_arc_runtime()
            ));
            assert!(Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()));
            assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
            assert_eq!(
                molecule.derived_cache_runtime().valid_states(),
                valid_states
            );
        };

        let assert_outputs = |params: &MorganFingerprintParams, expected: Expected, label: &str| {
            let sparse_counts = molecule
                .morgan_sparse_count_fingerprint_with_params(params, None)
                .unwrap_or_else(|error| panic!("{label}: sparse-count call failed: {error:?}"));
            assert_eq!(
                sparse_counts.length(),
                u64::MAX,
                "{label}, sparse count length"
            );
            assert_eq!(
                sparse_counts.nonzero_elements(),
                &expected
                    .sparse_counts
                    .iter()
                    .copied()
                    .collect::<BTreeMap<_, _>>(),
                "{label}, source raw sparse counts"
            );
            assert_source_unchanged();

            let sparse_bits = molecule
                .morgan_sparse_fingerprint_with_params(params, None)
                .unwrap_or_else(|error| panic!("{label}: sparse-bit call failed: {error:?}"));
            assert_eq!(sparse_bits.n_bits(), u32::MAX, "{label}, sparse bit length");
            assert_eq!(
                sparse_bits.on_bits().as_slice(),
                expected.sparse_bits,
                "{label}, source signed sparse bits"
            );
            assert_source_unchanged();

            let counts = molecule
                .morgan_count_fingerprint_with_params(params, None)
                .unwrap_or_else(|error| panic!("{label}: count call failed: {error:?}"));
            assert_eq!(counts.length(), 2048, "{label}, hashed count length");
            assert_eq!(
                counts.nonzero_elements(),
                &expected.counts.iter().copied().collect::<BTreeMap<_, _>>(),
                "{label}, source hashed counts"
            );
            assert_source_unchanged();

            let bits = molecule
                .morgan_fingerprint_with_params(params, None)
                .unwrap_or_else(|error| panic!("{label}: dense-bit call failed: {error:?}"));
            assert_eq!(bits.n_bits(), 2048, "{label}, dense bit length");
            assert_eq!(
                bits.on_bits().as_slice(),
                expected.bits,
                "{label}, dense bits"
            );
            assert_source_unchanged();
        };

        macro_rules! assert_all_four_preconditions {
            ($params:expr, $reason:expr, $label:expr) => {{
                let error = match molecule
                    .morgan_sparse_count_fingerprint_with_params($params, None)
                {
                    Err(error) => error,
                    Ok(_) => panic!("{}: sparse-count call unexpectedly succeeded", $label),
                };
                assert!(matches!(
                    error,
                    MorganReadError::Generator(MorganError::Fingerprint(
                        FingerprintError::PreconditionViolation { what }
                    )) if what == $reason
                ), "{}: unexpected sparse-count error: {error:?}", $label);
                assert_source_unchanged();

                let error = match molecule.morgan_sparse_fingerprint_with_params($params, None) {
                    Err(error) => error,
                    Ok(_) => panic!("{}: sparse-bit call unexpectedly succeeded", $label),
                };
                assert!(matches!(
                    error,
                    MorganReadError::Generator(MorganError::Fingerprint(
                        FingerprintError::PreconditionViolation { what }
                    )) if what == $reason
                ), "{}: unexpected sparse-bit error: {error:?}", $label);
                assert_source_unchanged();

                let error = match molecule.morgan_count_fingerprint_with_params($params, None) {
                    Err(error) => error,
                    Ok(_) => panic!("{}: count call unexpectedly succeeded", $label),
                };
                assert!(matches!(
                    error,
                    MorganReadError::Generator(MorganError::Fingerprint(
                        FingerprintError::PreconditionViolation { what }
                    )) if what == $reason
                ), "{}: unexpected count error: {error:?}", $label);
                assert_source_unchanged();

                let error = match molecule.morgan_fingerprint_with_params($params, None) {
                    Err(error) => error,
                    Ok(_) => panic!("{}: dense-bit call unexpectedly succeeded", $label),
                };
                assert!(matches!(
                    error,
                    MorganReadError::Generator(MorganError::Fingerprint(
                        FingerprintError::PreconditionViolation { what }
                    )) if what == $reason
                ), "{}: unexpected dense-bit error: {error:?}", $label);
                assert_source_unchanged();
            }};
        }

        let radius_zero = MorganParams {
            radius: 0,
            ..MorganParams::default()
        };

        // Literal argument-shape tables freeze None, present-empty, and
        // duplicate fromAtoms, plus absent/empty/short/exact/long custom
        // atom and bond invariant vectors before any public call is made.
        const FROM_ATOM_SHAPES: [(&str, Option<&[u32]>, Expected); 3] = [
            ("fromAtoms-none", None, DEFAULT_RADIUS_ZERO),
            ("fromAtoms-empty", Some(&[]), EMPTY),
            (
                "fromAtoms-duplicates",
                Some(&[2, 0, 1, 2, 0]),
                DEFAULT_RADIUS_ZERO,
            ),
        ];
        for (label, from_atoms, expected) in FROM_ATOM_SHAPES {
            let params = MorganFingerprintParams {
                generator: radius_zero.clone(),
                from_atoms: from_atoms.map(|indices| indices.to_vec()),
                ..MorganFingerprintParams::default()
            };
            assert_outputs(&params, expected, label);
        }

        let ignored_arguments = MorganFingerprintParams {
            generator: radius_zero.clone(),
            ignore_atoms: Some(vec![0, 1, 2]),
            conformer_id: 12345,
            ..MorganFingerprintParams::default()
        };
        assert_outputs(
            &ignored_arguments,
            DEFAULT_RADIUS_ZERO,
            "ignored ignoreAtoms and absent confId",
        );

        const ATOM_INVARIANT_SHAPES: [(&str, Option<&[u32]>, Option<&str>, Expected); 5] = [
            ("atom-invariants-none", None, None, DEFAULT_RADIUS_ZERO),
            (
                "atom-invariants-empty",
                Some(&[]),
                Some("bad atom invariants size"),
                EMPTY,
            ),
            (
                "atom-invariants-short",
                Some(&[11, 12]),
                Some("bad atom invariants size"),
                EMPTY,
            ),
            (
                "atom-invariants-exact",
                Some(&[1, 2, 1]),
                None,
                CUSTOM_ATOMS,
            ),
            (
                "atom-invariants-long",
                Some(&[1, 2, 1, 99]),
                None,
                CUSTOM_ATOMS,
            ),
        ];
        for (label, values, expected_error, expected) in ATOM_INVARIANT_SHAPES {
            let params = MorganFingerprintParams {
                generator: radius_zero.clone(),
                custom_atom_invariants: values.map(|items| items.to_vec()),
                ..MorganFingerprintParams::default()
            };
            if let Some(reason) = expected_error {
                assert_all_four_preconditions!(&params, reason, label);
            } else {
                assert_outputs(&params, expected, label);
            }
        }

        const BOND_INVARIANT_SHAPES: [(&str, Option<&[u32]>, Option<&str>, Expected); 5] = [
            ("bond-invariants-none", None, None, DEFAULT_RADIUS_ZERO),
            (
                "bond-invariants-empty",
                Some(&[]),
                Some("bad bond invariants size"),
                EMPTY,
            ),
            (
                "bond-invariants-short",
                Some(&[21]),
                Some("bad bond invariants size"),
                EMPTY,
            ),
            (
                "bond-invariants-exact",
                Some(&[21, 22]),
                None,
                DEFAULT_RADIUS_ZERO,
            ),
            (
                "bond-invariants-long",
                Some(&[21, 22, 99]),
                None,
                DEFAULT_RADIUS_ZERO,
            ),
        ];
        for (label, values, expected_error, expected) in BOND_INVARIANT_SHAPES {
            let params = MorganFingerprintParams {
                generator: radius_zero.clone(),
                custom_bond_invariants: values.map(|items| items.to_vec()),
                ..MorganFingerprintParams::default()
            };
            if let Some(reason) = expected_error {
                assert_all_four_preconditions!(&params, reason, label);
            } else {
                assert_outputs(&params, expected, label);
            }
        }

        let both_short = MorganFingerprintParams {
            generator: radius_zero.clone(),
            custom_atom_invariants: Some(vec![11, 12]),
            custom_bond_invariants: Some(vec![21]),
            ..MorganFingerprintParams::default()
        };
        assert_all_four_preconditions!(
            &both_short,
            "bad atom invariants size",
            "atom invariant error precedes bond invariant error"
        );

        let invariant_shapes = [
            (
                "connectivity",
                MorganInvariants::Connectivity,
                DEFAULT_RADIUS_ZERO,
            ),
            ("features", MorganInvariants::Features, FEATURES),
            (
                "feature-patterns-empty",
                MorganInvariants::FeaturePatterns(Vec::new()),
                EMPTY_FEATURE_PATTERNS,
            ),
        ];
        for (label, invariants, expected) in invariant_shapes {
            let params = MorganFingerprintParams {
                generator: radius_zero.clone(),
                invariants,
                ..MorganFingerprintParams::default()
            };
            assert_outputs(&params, expected, label);
        }
    }

    #[test]
    fn morgan_public_vector_contract_root_signatures() {
        let _: fn(&Fingerprint) -> u32 = Fingerprint::n_bits;
        let _: fn(&Fingerprint) -> Vec<u32> = Fingerprint::on_bits;
        let _: fn(&SparseBitFingerprint) -> u32 = SparseBitFingerprint::n_bits;
        let _: fn(&SparseBitFingerprint) -> Vec<i32> = SparseBitFingerprint::on_bits;
    }

    #[test]
    fn morgan_public_vector_contract_empty_and_dense_high_bits() {
        let empty = Fingerprint::empty();
        assert_eq!(empty.n_bits(), 0);
        assert_eq!(empty.on_bits(), Vec::<u32>::new());

        let mut dense = Fingerprint::new(65);
        for bit in [64, 0, 63, 31] {
            assert!(!dense.set_bit(bit).unwrap());
        }
        assert_eq!(dense.n_bits(), 65);
        assert_eq!(dense.on_bits(), vec![0, 31, 63, 64]);
    }

    #[test]
    fn morgan_public_vector_contract_sparse_high_and_source_maximum() {
        let mut high = SparseBitFingerprint::new(u32::MAX);
        for bit in [5, (1_u32 << 31) + 1] {
            assert!(!high.set_bit(bit).unwrap());
        }
        assert_eq!(high.n_bits(), u32::MAX);
        assert_eq!(high.on_bits(), vec![-2_147_483_647, 5]);

        let mut maximum = SparseBitFingerprint::new(u32::MAX);
        for bit in [7, u32::MAX] {
            assert!(!maximum.set_bit(bit).unwrap());
        }
        assert_eq!(maximum.n_bits(), u32::MAX);
        assert_eq!(maximum.on_bits(), vec![-1, 7]);

        let empty = SparseBitFingerprint::new(0);
        assert_eq!(empty.n_bits(), 0);
        assert_eq!(empty.on_bits(), Vec::<i32>::new());
    }
}
