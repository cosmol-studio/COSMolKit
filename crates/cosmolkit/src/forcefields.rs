use std::error::Error;
use std::fmt;

pub use cosmolkit_forcefields::{UffParameterError, UffParameterErrorKind};

use crate::ops::OperationError;

#[derive(Debug)]
pub enum UffParameterQueryError {
    Cache(OperationError),
    Parameters(cosmolkit_forcefields::UffParameterError),
}

impl fmt::Display for UffParameterQueryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cache(source) => fmt::Display::fmt(source, formatter),
            Self::Parameters(source) => fmt::Display::fmt(source, formatter),
        }
    }
}

impl Error for UffParameterQueryError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Cache(source) => Some(source),
            Self::Parameters(source) => Some(source),
        }
    }
}

impl crate::Molecule {
    /// Reports whether all atom types have UFF parameters using prepared state.
    pub fn uff_has_all_molecule_params(&self) -> Result<bool, UffParameterQueryError> {
        let topology = self.topology();
        let atom_count = topology.atoms.len();

        if atom_count == 0 {
            let empty_assignment = cosmolkit_core::ValenceAssignment {
                explicit_valence: Vec::new(),
                implicit_hydrogens: Vec::new(),
            };
            return cosmolkit_forcefields::uff_has_all_molecule_params(topology, &empty_assignment)
                .map_err(UffParameterQueryError::Parameters);
        }

        let cache = self.derived_cache_runtime();
        cache
            .validate_for_atom_count(atom_count)
            .map_err(UffParameterQueryError::Cache)?;
        let assignment = cache.valence_assignment().ok_or_else(|| {
            UffParameterQueryError::Cache(OperationError::InvalidDerivedCache {
                state: "valence",
                field: "assignment",
                actual: 0,
                expected: 1,
            })
        })?;

        cosmolkit_forcefields::uff_has_all_molecule_params(topology, assignment)
            .map_err(UffParameterQueryError::Parameters)
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error;
    use std::sync::Arc;

    use cosmolkit_core::ValenceAssignment;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Conformer2D, Conformer3D, CoordinateDimension, Element,
        Hybridization, PropertyValue, SdfPropertyList, SdfPropertyListTarget, TopologyBlock,
    };

    use super::UffParameterQueryError;
    use crate::molecule::DerivedCacheBlock;
    use crate::ops::{DerivedState, OperationError};
    use crate::{CoordinateBlock, Molecule, MoleculeProperties};

    fn run_query(molecule: &Molecule) -> Result<bool, UffParameterQueryError> {
        let query: fn(&Molecule) -> Result<bool, UffParameterQueryError> =
            Molecule::uff_has_all_molecule_params;
        query(molecule)
    }

    fn one_carbon_topology() -> TopologyBlock {
        let carbon = Element::from_atomic_number(6).expect("carbon is in the element table");
        let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(carbon));
        TopologyBlock::try_from_parts(vec![atom], Vec::new(), Vec::new(), Vec::new())
            .expect("one isolated carbon is structurally valid")
    }

    fn one_parameter_query_topology(
        atomic_number: u8,
        no_implicit: bool,
        dummy_label: Option<&str>,
    ) -> TopologyBlock {
        let element = Element::from_atomic_number(atomic_number)
            .expect("fixed parameter-query element is in the model table");
        let mut spec = AtomSpec::new(element)
            .with_hybridization(Hybridization::Sp3)
            .with_no_implicit(no_implicit);
        if let Some(label) = dummy_label {
            spec = spec
                .with_prop("dummyLabel", label)
                .expect("fixed dummyLabel is a valid atom property");
        }
        let atom = Atom::from_spec(AtomId::new(0), spec);
        TopologyBlock::try_from_parts(vec![atom], Vec::new(), Vec::new(), Vec::new())
            .expect("one source-aligned parameter-query atom is structurally valid")
    }

    fn run_counted_query(
        molecule: &Molecule,
        actual_calls: &mut usize,
    ) -> Result<bool, UffParameterQueryError> {
        *actual_calls += 1;
        run_query(molecule)
    }

    fn assert_exact_uncached_query_error(error: &UffParameterQueryError) {
        let UffParameterQueryError::Cache(OperationError::InvalidDerivedCache {
            state,
            field,
            actual,
            expected,
        }) = error
        else {
            panic!("uncached nonempty query did not retain its exact cache error: {error:?}");
        };
        assert_eq!(
            (*state, *field, *actual, *expected),
            ("valence", "assignment", 0, 1)
        );
        let stored = match error {
            UffParameterQueryError::Cache(stored) => stored,
            UffParameterQueryError::Parameters(_) => unreachable!(),
        };
        let reported = Error::source(error).expect("cache error remains the borrowed cause");
        assert!(std::ptr::eq(
            reported
                .downcast_ref::<OperationError>()
                .expect("concrete cache cause"),
            stored,
        ));
        assert!(Error::source(stored).is_none());
    }

    #[derive(Clone, Copy)]
    enum QueryCoordinateFixture {
        None,
        TwoD,
        ThreeD,
    }

    fn p09_coordinate_block(fixture: QueryCoordinateFixture, atom_count: usize) -> CoordinateBlock {
        match fixture {
            QueryCoordinateFixture::None => CoordinateBlock::default(),
            QueryCoordinateFixture::TwoD => CoordinateBlock {
                conformers_2d: vec![
                    Conformer2D::new(
                        17,
                        (0..atom_count)
                            .map(|row| [-0.0, row as f64 + 0.5])
                            .collect(),
                    )
                    .with_prop("uff-param-p09-dimension", "2d"),
                ],
                source_coordinate_dim: Some(CoordinateDimension::TwoD),
                ..CoordinateBlock::default()
            },
            QueryCoordinateFixture::ThreeD => CoordinateBlock {
                conformers_3d: vec![
                    Conformer3D::new(
                        29,
                        (0..atom_count)
                            .map(|row| [-0.0, row as f64 + 0.25, -(row as f64 + 0.75)])
                            .collect(),
                        true,
                    )
                    .with_prop("uff-param-p09-dimension", "3d"),
                ],
                source_coordinate_dim: Some(CoordinateDimension::ThreeD),
                ..CoordinateBlock::default()
            },
        }
    }

    fn p09_coordinate_bits(coordinates: &CoordinateBlock) -> Vec<u64> {
        let mut bits = Vec::new();
        for conformer in &coordinates.conformers_2d {
            for row in conformer.coordinates() {
                bits.extend(row.iter().map(|value| value.to_bits()));
            }
        }
        for conformer in &coordinates.conformers_3d {
            for row in conformer.coordinates() {
                bits.extend(row.iter().map(|value| value.to_bits()));
            }
        }
        bits
    }

    fn p09_properties(
        fixture_id: usize,
        atom_count: usize,
        bond_count: usize,
    ) -> MoleculeProperties {
        let atom_values = (0..atom_count)
            .map(|row| {
                (row % 2 == 0)
                    .then(|| PropertyValue::String(format!("atom-{fixture_id}-{row}").into()))
            })
            .collect();
        let bond_values = (0..bond_count)
            .map(|row| {
                (row % 2 == 0)
                    .then(|| PropertyValue::String(format!("bond-{fixture_id}-{row}").into()))
            })
            .collect();

        MoleculeProperties::default()
            .with_name(format!("UFF-P09 source fixture {fixture_id}"))
            .with_prop("raw-source-metadata", format!("row-{fixture_id}"))
            .expect("fixed raw source metadata key is valid")
            .with_computed_prop("computed-source-metadata", format!("computed-{fixture_id}"))
            .expect("fixed computed metadata key is valid")
            .with_sdf_data_field("P09_DUPLICATE", format!("first-{fixture_id}"))
            .with_sdf_data_field("P09_DUPLICATE", format!("second-{fixture_id}"))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "P09_ATOM_ROWS",
                atom_values,
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "P09_BOND_ROWS",
                bond_values,
            ))
    }

    struct QueryStorageSnapshot {
        topology_arc: Arc<TopologyBlock>,
        coordinate_arc: Arc<CoordinateBlock>,
        properties_arc: Arc<MoleculeProperties>,
        cache_arc: Arc<DerivedCacheBlock>,
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        cache: DerivedCacheBlock,
        coordinate_bits: Vec<u64>,
        valid_states: DerivedState,
        assignment: Option<ValenceAssignment>,
    }

    impl QueryStorageSnapshot {
        fn capture(molecule: &Molecule) -> Self {
            let coordinates = molecule.coordinate_block_runtime().clone();
            let cache = molecule.derived_cache_runtime();
            Self {
                topology_arc: molecule.topology_arc_runtime(),
                coordinate_arc: molecule.coordinates_arc_runtime(),
                properties_arc: molecule.properties_arc_runtime(),
                cache_arc: molecule.derived_cache_arc_runtime(),
                topology: molecule.topology().clone(),
                coordinate_bits: p09_coordinate_bits(&coordinates),
                coordinates,
                properties: molecule.properties().clone(),
                cache: cache.clone(),
                valid_states: cache.valid_states(),
                assignment: cache.valence_assignment().cloned(),
            }
        }

        fn assert_unchanged(&self, molecule: &Molecule) {
            assert!(Arc::ptr_eq(
                &self.topology_arc,
                &molecule.topology_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &self.coordinate_arc,
                &molecule.coordinates_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &self.properties_arc,
                &molecule.properties_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &self.cache_arc,
                &molecule.derived_cache_arc_runtime()
            ));
            assert_eq!(molecule.topology(), &self.topology);
            assert_eq!(molecule.coordinate_block_runtime(), &self.coordinates);
            assert_eq!(
                p09_coordinate_bits(molecule.coordinate_block_runtime()),
                self.coordinate_bits
            );
            assert_eq!(molecule.properties(), &self.properties);
            let cache = molecule.derived_cache_runtime();
            assert_eq!(cache, &self.cache);
            assert_eq!(cache.valid_states(), self.valid_states);
            assert_eq!(cache.valence_assignment(), self.assignment.as_ref());
        }

        fn assert_shared_blocks(&self, peer: &Self) {
            assert!(Arc::ptr_eq(&self.topology_arc, &peer.topology_arc));
            assert!(Arc::ptr_eq(&self.coordinate_arc, &peer.coordinate_arc));
            assert!(Arc::ptr_eq(&self.properties_arc, &peer.properties_arc));
            assert!(Arc::ptr_eq(&self.cache_arc, &peer.cache_arc));
        }
    }

    fn assert_query_does_not_change_prepared_state(molecule: &Molecule, expected: bool) {
        let cache = molecule.derived_cache_runtime();
        let cache_arc = molecule.derived_cache_arc_runtime();
        let valid_states = cache.valid_states();
        let assignment = cache.valence_assignment().cloned();
        let atoms = molecule.atoms().to_vec();

        assert_eq!(
            run_query(molecule).expect("prepared query succeeds"),
            expected
        );

        assert!(Arc::ptr_eq(
            &cache_arc,
            &molecule.derived_cache_arc_runtime()
        ));
        assert_eq!(
            molecule.derived_cache_runtime().valid_states(),
            valid_states
        );
        assert_eq!(
            molecule.derived_cache_runtime().valence_assignment(),
            assignment.as_ref()
        );
        assert_eq!(molecule.atoms(), atoms);
    }

    fn assert_uncached_query_fails_without_preparation(molecule: &Molecule) {
        let cache_arc = molecule.derived_cache_arc_runtime();
        let atoms = molecule.atoms().to_vec();
        assert!(
            !molecule
                .derived_cache_runtime()
                .valid_states()
                .contains(DerivedState::VALENCE)
        );
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_none()
        );

        let error = run_query(molecule).expect_err("nonempty query requires stored valence");
        assert_eq!(
            error.to_string(),
            "derived cache state `valence` field `assignment` has 0 entries, expected 1"
        );
        let stored = match &error {
            UffParameterQueryError::Cache(OperationError::InvalidDerivedCache {
                state,
                field,
                actual,
                expected,
            }) => {
                assert_eq!(
                    (*state, *field, *actual, *expected),
                    ("valence", "assignment", 0, 1)
                );
                match &error {
                    UffParameterQueryError::Cache(stored) => stored,
                    UffParameterQueryError::Parameters(_) => unreachable!(),
                }
            }
            UffParameterQueryError::Cache(other) => {
                panic!("unexpected cache error: {other:?}")
            }
            UffParameterQueryError::Parameters(other) => {
                panic!("uncached state was mapped to domain parameters: {other}")
            }
        };
        let reported = Error::source(&error).expect("cache error remains the borrowed cause");
        // Compare concrete borrowed storage, not a dyn Error vtable address:
        // distinct codegen units may emit different vtables for the same type.
        assert!(std::ptr::eq(
            reported
                .downcast_ref::<OperationError>()
                .expect("concrete cache cause"),
            stored,
        ));
        assert!(Error::source(stored).is_none());

        assert!(Arc::ptr_eq(
            &cache_arc,
            &molecule.derived_cache_arc_runtime()
        ));
        assert!(
            !molecule
                .derived_cache_runtime()
                .valid_states()
                .contains(DerivedState::VALENCE)
        );
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_none()
        );
        assert_eq!(molecule.atoms(), atoms);
    }

    fn molecule_with_cached_assignment(
        topology: TopologyBlock,
        assignment: ValenceAssignment,
    ) -> Molecule {
        let mut cache = DerivedCacheBlock::default();
        cache.install_valence_assignment(assignment);
        cache.mark_valid(DerivedState::VALENCE);
        Molecule::from_runtime_parts(
            Arc::new(topology),
            Arc::new(CoordinateBlock::default()),
            Arc::new(MoleculeProperties::default()),
            Arc::new(cache),
        )
        .expect("the synthetic cache fixture has valid row counts")
    }

    // This manually constructed child proves wrapper trait dispatch only; it
    // is not evidence for a live Molecule cache failure.
    #[test]
    fn uff_param_p05_manual_cache_error_trait_dispatch_preserves_stored_child() {
        let child = OperationError::InvalidDerivedCache {
            state: "valence",
            field: "explicit_valence",
            actual: 2,
            expected: 5,
        };
        let expected_display =
            "derived cache state `valence` field `explicit_valence` has 2 entries, expected 5";
        let error = UffParameterQueryError::Cache(child);

        assert_eq!(error.to_string(), expected_display);
        let stored = match &error {
            UffParameterQueryError::Cache(stored) => stored,
            UffParameterQueryError::Parameters(_) => unreachable!(),
        };
        let reported = Error::source(&error).expect("cache cause remains available");
        assert!(std::ptr::eq(
            reported
                .downcast_ref::<OperationError>()
                .expect("concrete cache cause"),
            stored,
        ));
        assert_eq!(
            reported.downcast_ref::<OperationError>(),
            Some(stored),
            "the wrapper retains the original OperationError value"
        );
        assert!(Error::source(stored).is_none());
    }

    #[test]
    fn uff_param_p05_actual_detached_preparation_error_keeps_nested_source() {
        let topology = TopologyBlock::default();
        let assignment = ValenceAssignment {
            explicit_valence: vec![1],
            implicit_hydrogens: Vec::new(),
        };
        let child = cosmolkit_forcefields::uff_has_all_molecule_params(&topology, &assignment)
            .expect_err("the detached helper validates the supplied cache row count");
        let expected_display =
            "ValenceAssignmentLengthMismatch { field: Explicit, expected: 0, actual: 1 }";
        assert_eq!(child.to_string(), expected_display);

        let error = UffParameterQueryError::Parameters(child);
        assert_eq!(error.to_string(), expected_display);
        let stored = match &error {
            UffParameterQueryError::Parameters(stored) => stored,
            UffParameterQueryError::Cache(_) => unreachable!(),
        };
        assert_eq!(
            stored.kind(),
            cosmolkit_forcefields::UffParameterErrorKind::Preparation
        );
        let reported = Error::source(&error).expect("domain preparation error remains available");
        let downcast = reported
            .downcast_ref::<cosmolkit_forcefields::UffParameterError>()
            .expect("the root wrapper retains the concrete domain error type");
        assert!(
            std::ptr::eq(downcast, stored),
            "the root wrapper retains the original domain error value"
        );

        let preparation =
            Error::source(stored).expect("domain error retains its concrete preparation cause");
        assert_eq!(preparation.to_string(), expected_display);
        assert!(
            Error::source(preparation).is_none(),
            "the source preparation variant has no nested cause"
        );
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn uff_param_p06_cached_known_type_returns_true_without_state_change() {
        let molecule = Molecule::from_smiles("C")
            .expect("fixed carbon SMILES parses")
            .with_assigned_valence()
            .expect("the existing valence operation prepares carbon");
        assert_query_does_not_change_prepared_state(&molecule, true);
    }

    #[cfg(all(feature = "cap-smiles", feature = "cap-valence"))]
    #[test]
    fn uff_param_p06_cached_missing_type_returns_false_without_state_change() {
        let molecule = Molecule::from_smiles("[Cu+]")
            .expect("fixed copper SMILES parses")
            .with_assigned_valence()
            .expect("the existing valence operation prepares copper");
        assert_query_does_not_change_prepared_state(&molecule, false);
    }

    #[cfg(feature = "cap-smiles")]
    #[test]
    fn uff_param_p06_from_parts_and_smiles_require_cached_valence() {
        let from_parts = Molecule::from_parts(
            one_carbon_topology(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .expect("valid detached parts construct a molecule");
        assert_uncached_query_fails_without_preparation(&from_parts);

        // Sanitized construction retains its final assignment. No second
        // with_assigned_valence operation is needed for this unchanged graph.
        let from_smiles = Molecule::from_smiles("C").expect("fixed carbon SMILES parses");
        assert_query_does_not_change_prepared_state(&from_smiles, true);

        // Constructor policy, not a UFF heuristic, controls cache validity.
        // Keep the prior missing-cache/error/sharing proof for unsanitized
        // inputs and cover the entire two-flag product without dropping rows.
        let mut calls = 0;
        for sanitize in [false, true] {
            for remove_hs in [false, true] {
                let molecule = Molecule::from_smiles_with_params(
                    "C",
                    &crate::SmilesParseParams {
                        sanitize,
                        remove_hs,
                        ..Default::default()
                    },
                )
                .expect("fixed carbon parses under both constructor policies");
                if sanitize || remove_hs {
                    assert_eq!(
                        molecule.derived_cache_runtime().valence_assignment(),
                        Some(&ValenceAssignment {
                            explicit_valence: vec![0],
                            implicit_hydrogens: vec![4]
                        })
                    );
                    // Unsanitized carbon retains Unspecified hybridization: the
                    // source UFF label has no matching parameter row.
                    assert_query_does_not_change_prepared_state(&molecule, sanitize);
                } else {
                    assert_uncached_query_fails_without_preparation(&molecule);
                }
                calls += 1;
            }
        }
        assert_eq!(calls, 4);
    }

    #[test]
    fn uff_param_p06_empty_uncached_query_is_true_without_state_change() {
        let molecule = Molecule::new();
        assert!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .is_none()
        );
        assert!(
            !molecule
                .derived_cache_runtime()
                .valid_states()
                .contains(DerivedState::VALENCE)
        );
        assert_query_does_not_change_prepared_state(&molecule, true);
    }

    #[test]
    fn uff_param_p06_real_domain_preparation_error_remains_parameters() {
        // A private row-count-valid cache fixture exercises forwarding of the
        // real detached signed-byte-domain error; it is not a public cache
        // construction path or a claim that with_assigned_valence emits 128.
        let molecule = molecule_with_cached_assignment(
            one_carbon_topology(),
            ValenceAssignment {
                explicit_valence: vec![128],
                implicit_hydrogens: vec![0],
            },
        );
        let error = run_query(&molecule).expect_err("the detached source cache domain is enforced");
        let UffParameterQueryError::Parameters(parameters) = &error else {
            panic!("source preparation failure was not retained as Parameters: {error}");
        };
        assert_eq!(
            parameters.kind(),
            cosmolkit_forcefields::UffParameterErrorKind::Preparation
        );
        let forwarded = Error::source(&error).expect("domain error is the borrowed root cause");
        assert!(std::ptr::eq(
            forwarded
                .downcast_ref::<cosmolkit_forcefields::UffParameterError>()
                .expect("concrete domain cause"),
            parameters,
        ));
        let preparation =
            Error::source(parameters).expect("typed preparation cause remains available");
        assert_eq!(
            preparation.to_string(),
            "SourceValenceOutOfRange { atom_id: AtomId(0), field: Explicit, value: 128 }"
        );
        assert!(Error::source(preparation).is_none());
    }

    #[cfg(feature = "cap-valence")]
    #[test]
    fn uff_param_p08_public_query_cache_and_atom_type_matrix() {
        // Fixed source-aligned cases: C_3 is in the default table; the dummy
        // symbol "*" has no row. `no_implicit` is an independent Atom state.
        let rows = [
            (6, false, None, true),
            (6, true, None, true),
            (0, false, Some("*"), false),
            (0, true, Some("*"), false),
        ];
        let mut actual_calls = 0;

        for (atomic_number, no_implicit, dummy_label, expected) in rows {
            let molecule = Molecule::from_parts(
                one_parameter_query_topology(atomic_number, no_implicit, dummy_label),
                CoordinateBlock::default(),
                MoleculeProperties::default(),
            )
            .expect("fixed source-aligned AtomSpec builds a valid live molecule");
            let uncached_arc = molecule.derived_cache_arc_runtime();
            let uncached_states = molecule.derived_cache_runtime().valid_states();
            let uncached_atoms = molecule.atoms().to_vec();
            assert!(!uncached_states.contains(DerivedState::VALENCE));
            assert!(
                molecule
                    .derived_cache_runtime()
                    .valence_assignment()
                    .is_none()
            );

            let error = run_counted_query(&molecule, &mut actual_calls)
                .expect_err("nonempty public query requires an existing valid assignment");
            assert_exact_uncached_query_error(&error);
            assert!(Arc::ptr_eq(
                &uncached_arc,
                &molecule.derived_cache_arc_runtime()
            ));
            assert_eq!(
                molecule.derived_cache_runtime().valid_states(),
                uncached_states
            );
            assert!(
                molecule
                    .derived_cache_runtime()
                    .valence_assignment()
                    .is_none()
            );
            assert_eq!(molecule.atoms(), uncached_atoms);

            let prepared = molecule
                .with_assigned_valence()
                .expect("the existing public operation installs the valid valence state");
            assert!(
                !molecule
                    .derived_cache_runtime()
                    .valid_states()
                    .contains(DerivedState::VALENCE)
            );
            assert!(
                molecule
                    .derived_cache_runtime()
                    .valence_assignment()
                    .is_none()
            );
            let prepared_arc = prepared.derived_cache_arc_runtime();
            let prepared_states = prepared.derived_cache_runtime().valid_states();
            let prepared_assignment = prepared
                .derived_cache_runtime()
                .valence_assignment()
                .cloned();
            let prepared_atoms = prepared.atoms().to_vec();
            assert!(prepared_states.contains(DerivedState::VALENCE));
            assert!(prepared_assignment.is_some());

            assert_eq!(
                run_counted_query(&prepared, &mut actual_calls)
                    .expect("prepared known/unknown type is a bool result"),
                expected
            );
            assert!(Arc::ptr_eq(
                &prepared_arc,
                &prepared.derived_cache_arc_runtime()
            ));
            assert_eq!(
                prepared.derived_cache_runtime().valid_states(),
                prepared_states
            );
            assert_eq!(
                prepared.derived_cache_runtime().valence_assignment(),
                prepared_assignment.as_ref()
            );
            assert_eq!(prepared.atoms(), prepared_atoms);
        }

        let empty = Molecule::new();
        let prepared_empty = empty
            .with_assigned_valence()
            .expect("the existing valence operation accepts an empty molecule");
        for (molecule, expect_cached) in [(&empty, false), (&prepared_empty, true)] {
            let cache_arc = molecule.derived_cache_arc_runtime();
            let valid_states = molecule.derived_cache_runtime().valid_states();
            let assignment = molecule
                .derived_cache_runtime()
                .valence_assignment()
                .cloned();
            assert_eq!(valid_states.contains(DerivedState::VALENCE), expect_cached);
            if expect_cached {
                assert_eq!(
                    assignment,
                    Some(ValenceAssignment {
                        explicit_valence: Vec::new(),
                        implicit_hydrogens: Vec::new(),
                    })
                );
            } else {
                assert!(assignment.is_none());
            }

            assert!(
                run_counted_query(molecule, &mut actual_calls)
                    .expect("empty input reads no atom cache and returns true")
            );
            assert!(Arc::ptr_eq(
                &cache_arc,
                &molecule.derived_cache_arc_runtime()
            ));
            assert_eq!(
                molecule.derived_cache_runtime().valid_states(),
                valid_states
            );
            assert_eq!(
                molecule.derived_cache_runtime().valence_assignment(),
                assignment.as_ref()
            );
            assert!(molecule.atoms().is_empty());
        }

        assert_eq!(actual_calls, 10, "8 nonempty and 2 empty public calls");
    }

    #[test]
    fn uff_param_p08_detached_atom_cache_component_boundaries() {
        #[derive(Clone, Copy)]
        enum Expected {
            Value(bool),
            SourcePrecondition { field: &'static str, value: i32 },
        }

        let rows = [
            (false, 0, 4, Expected::Value(true)),
            (
                false,
                -1,
                4,
                Expected::SourcePrecondition {
                    field: "Explicit",
                    value: -1,
                },
            ),
            (
                false,
                0,
                -1,
                Expected::SourcePrecondition {
                    field: "ImplicitHydrogen",
                    value: -1,
                },
            ),
            (
                false,
                -1,
                -1,
                Expected::SourcePrecondition {
                    field: "Explicit",
                    value: -1,
                },
            ),
            (true, 0, 0, Expected::Value(true)),
            (
                true,
                -1,
                0,
                Expected::SourcePrecondition {
                    field: "Explicit",
                    value: -1,
                },
            ),
            (true, 0, -1, Expected::Value(true)),
            (
                true,
                -1,
                -1,
                Expected::SourcePrecondition {
                    field: "Explicit",
                    value: -1,
                },
            ),
        ];
        let mut actual_calls = 0;

        for (no_implicit, explicit, implicit, expected) in rows {
            let topology = one_parameter_query_topology(6, no_implicit, None);
            let assignment = ValenceAssignment {
                explicit_valence: vec![explicit],
                implicit_hydrogens: vec![implicit],
            };
            actual_calls += 1;
            let result = cosmolkit_forcefields::uff_has_all_molecule_params(&topology, &assignment);

            match expected {
                Expected::Value(value) => assert_eq!(result.expect("valid source state"), value),
                Expected::SourcePrecondition { field, value } => {
                    let error = result.expect_err("source getValence precondition is preserved");
                    assert_eq!(error.kind(), super::UffParameterErrorKind::Preparation);
                    let expected_display = format!(
                        "SourceValencePrecondition {{ atom_id: AtomId(0), field: {field}, value: {value} }}"
                    );
                    assert_eq!(error.to_string(), expected_display);
                    let first_source = Error::source(&error)
                        .expect("preparation errors retain their concrete source");
                    let second_source = Error::source(&error)
                        .expect("the same stored preparation source remains available");
                    assert!(std::ptr::eq(first_source, second_source));
                    assert_eq!(first_source.to_string(), expected_display);
                    assert!(Error::source(first_source).is_none());
                }
            }
        }

        assert_eq!(
            actual_calls, 8,
            "two no-implicit x four component-state calls"
        );
    }

    #[cfg(all(
        feature = "cap-smiles",
        feature = "cap-hydrogens",
        feature = "cap-valence"
    ))]
    #[test]
    fn uff_param_p09_query_preserves_shared_storage_across_coordinate_states() {
        let cases: [(&str, [bool; 2]); 10] = [
            ("", [true, true]),
            ("CC", [true, true]),
            ("CCO", [true, true]),
            ("c1ccccc1", [true, true]),
            ("CC(=O)O", [true, true]),
            ("N", [true, true]),
            ("[Na+].[Cl-]", [true, true]),
            ("*", [false, false]),
            ("[Cu+]", [false, false]),
            ("[SiH4]", [true, true]),
        ];
        let mut primary_calls = 0;
        let mut repeat_calls = 0;

        for (case_index, (smiles, expected_by_hydrogen_option)) in cases.into_iter().enumerate() {
            for (add_hydrogens, expected) in [
                (false, expected_by_hydrogen_option[0]),
                (true, expected_by_hydrogen_option[1]),
            ] {
                let parsed = Molecule::from_smiles(smiles)
                    .unwrap_or_else(|error| panic!("fixed P07 SMILES {smiles:?}: {error}"));
                let chemistry = if add_hydrogens {
                    parsed
                        .with_hydrogens()
                        .unwrap_or_else(|error| panic!("P07 AddHs for {smiles:?}: {error}"))
                } else {
                    parsed
                };
                let topology = chemistry.topology().clone();
                let properties = p09_properties(
                    case_index * 2 + usize::from(add_hydrogens),
                    topology.atoms.len(),
                    topology.bonds.len(),
                );
                let with_source_metadata =
                    Molecule::from_parts(topology.clone(), CoordinateBlock::default(), properties)
                        .expect("P09 metadata fixture validates against the P07 topology");
                let prepared = with_source_metadata
                    .with_assigned_valence()
                    .unwrap_or_else(|error| panic!("P09 valence for {smiles:?}: {error}"));
                assert_eq!(prepared.topology(), &topology);

                let trusted_topology = prepared.topology_arc_runtime();
                let trusted_properties = prepared.properties_arc_runtime();
                let trusted_cache = prepared.derived_cache_arc_runtime();
                let trusted_assignment = prepared
                    .derived_cache_runtime()
                    .valence_assignment()
                    .cloned()
                    .expect("the existing valence operation installed a trusted assignment");
                assert!(
                    prepared
                        .derived_cache_runtime()
                        .valid_states()
                        .contains(DerivedState::VALENCE)
                );

                let mut variants = vec![prepared];
                for fixture in [QueryCoordinateFixture::TwoD, QueryCoordinateFixture::ThreeD] {
                    let coordinates = p09_coordinate_block(fixture, topology.atoms.len());
                    let coordinate_bits = p09_coordinate_bits(&coordinates);
                    if topology.atoms.is_empty() {
                        assert!(
                            coordinate_bits.is_empty(),
                            "empty topology has valid zero-row coordinate conformers"
                        );
                    } else {
                        assert!(coordinate_bits.contains(&(-0.0_f64).to_bits()));
                    }
                    let coordinate_arc = Arc::new(coordinates);
                    let variant = Molecule::from_runtime_parts(
                        Arc::clone(&trusted_topology),
                        coordinate_arc,
                        Arc::clone(&trusted_properties),
                        Arc::clone(&trusted_cache),
                    )
                    .expect("existing validated runtime constructor accepts aligned coordinates");
                    assert!(Arc::ptr_eq(
                        &trusted_topology,
                        &variant.topology_arc_runtime()
                    ));
                    assert!(Arc::ptr_eq(
                        &trusted_properties,
                        &variant.properties_arc_runtime()
                    ));
                    assert!(Arc::ptr_eq(
                        &trusted_cache,
                        &variant.derived_cache_arc_runtime()
                    ));
                    assert_eq!(
                        variant.derived_cache_runtime().valence_assignment(),
                        Some(&trusted_assignment)
                    );
                    variants.push(variant);
                }

                assert_eq!(variants.len(), 3, "no-coordinate, 2D, and 3D variants");
                for original in &variants {
                    let peer = original.clone();
                    let original_snapshot = QueryStorageSnapshot::capture(original);
                    let peer_snapshot = QueryStorageSnapshot::capture(&peer);
                    original_snapshot.assert_shared_blocks(&peer_snapshot);

                    original_snapshot.assert_unchanged(original);
                    peer_snapshot.assert_unchanged(&peer);
                    assert_eq!(
                        run_counted_query(original, &mut primary_calls).unwrap_or_else(
                            |error| panic!("P09 primary query {smiles:?}: {error}")
                        ),
                        expected,
                        "input={smiles:?}, add_h={add_hydrogens}"
                    );
                    original_snapshot.assert_unchanged(original);
                    peer_snapshot.assert_unchanged(&peer);

                    original_snapshot.assert_unchanged(original);
                    peer_snapshot.assert_unchanged(&peer);
                    assert_eq!(
                        run_counted_query(&peer, &mut repeat_calls)
                            .unwrap_or_else(|error| panic!("P09 peer repeat {smiles:?}: {error}")),
                        expected,
                        "peer repeat input={smiles:?}, add_h={add_hydrogens}"
                    );
                    original_snapshot.assert_unchanged(original);
                    peer_snapshot.assert_unchanged(&peer);
                }
            }
        }

        assert_eq!(primary_calls, 60, "20 fixed P07 rows x 3 coordinate states");
        assert_eq!(
            repeat_calls, 60,
            "one separately counted peer repeat per fixture"
        );

        let uncached_topology = Molecule::from_smiles("C")
            .expect("fixed uncached P09 source row parses")
            .topology()
            .clone();
        let mut uncached_primary_calls = 0;
        let mut uncached_repeat_calls = 0;
        for (fixture_index, fixture) in [
            QueryCoordinateFixture::None,
            QueryCoordinateFixture::TwoD,
            QueryCoordinateFixture::ThreeD,
        ]
        .into_iter()
        .enumerate()
        {
            let uncached = Molecule::from_parts(
                uncached_topology.clone(),
                p09_coordinate_block(fixture, uncached_topology.atoms.len()),
                p09_properties(
                    100 + fixture_index,
                    uncached_topology.atoms.len(),
                    uncached_topology.bonds.len(),
                ),
            )
            .expect("uncached error fixture is a validated live molecule");
            assert!(
                !uncached
                    .derived_cache_runtime()
                    .valid_states()
                    .contains(DerivedState::VALENCE)
            );
            assert!(
                uncached
                    .derived_cache_runtime()
                    .valence_assignment()
                    .is_none()
            );
            let peer = uncached.clone();
            let original_snapshot = QueryStorageSnapshot::capture(&uncached);
            let peer_snapshot = QueryStorageSnapshot::capture(&peer);
            original_snapshot.assert_shared_blocks(&peer_snapshot);

            original_snapshot.assert_unchanged(&uncached);
            peer_snapshot.assert_unchanged(&peer);
            let error = run_counted_query(&uncached, &mut uncached_primary_calls)
                .expect_err("real nonempty uncached query preserves its cache error");
            assert_exact_uncached_query_error(&error);
            original_snapshot.assert_unchanged(&uncached);
            peer_snapshot.assert_unchanged(&peer);

            original_snapshot.assert_unchanged(&uncached);
            peer_snapshot.assert_unchanged(&peer);
            let error = run_counted_query(&peer, &mut uncached_repeat_calls)
                .expect_err("peer repeat preserves the same real cache error");
            assert_exact_uncached_query_error(&error);
            original_snapshot.assert_unchanged(&uncached);
            peer_snapshot.assert_unchanged(&peer);
        }

        assert_eq!(
            uncached_primary_calls, 3,
            "one real uncached failure per dimension"
        );
        assert_eq!(
            uncached_repeat_calls, 3,
            "one separately counted uncached peer repeat per dimension"
        );
        assert_eq!(
            primary_calls + repeat_calls + uncached_primary_calls + uncached_repeat_calls,
            126,
            "60 matrix calls, 60 prepared repeats, and 6 uncached error calls"
        );
    }
}

pub use cosmolkit_forcefields::{
    MmffAtomProperties, MmffMolPropertiesError, MmffProperties, MmffPropertiesParams, MmffVariant,
};

impl crate::Molecule {
    /// Whether the default MMFF94 parameter set covers all input atoms.
    pub fn mmff_has_all_molecule_params(&self) -> Result<bool, MmffMolPropertiesError> {
        cosmolkit_forcefields::mmff_has_all_molecule_params(
            self.topology(),
            self.properties().prop("_MMFFSanitized").is_some(),
            self.derived_cache_runtime().valid_ring_info(),
        )
    }

    /// MMFF94 atom types, formal charges and partial charges in atom order.
    pub fn mmff_properties(&self) -> Result<MmffProperties, MmffMolPropertiesError> {
        self.mmff_properties_with_params(&MmffPropertiesParams::default())
    }

    /// Atom types and charges using the source-defined variant selection.
    pub fn mmff_properties_with_params(
        &self,
        params: &MmffPropertiesParams,
    ) -> Result<MmffProperties, MmffMolPropertiesError> {
        cosmolkit_forcefields::mmff_properties(
            self.topology(),
            self.properties().prop("_MMFFSanitized").is_some(),
            params,
            self.derived_cache_runtime().valid_ring_info(),
        )
    }
}

pub use cosmolkit_forcefields::{MmffEnergyGradient, MmffEvaluationParams};
impl crate::Molecule {
    pub fn mmff_energy_gradient(&self) -> Result<Option<MmffEnergyGradient>, OperationError> {
        self.mmff_energy_gradient_with_params(&MmffEvaluationParams::default())
    }
    pub fn mmff_energy_gradient_with_params(
        &self,
        params: &MmffEvaluationParams,
    ) -> Result<Option<MmffEnergyGradient>, OperationError> {
        cosmolkit_forcefields::evaluate_mmff(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            params,
            self.derived_cache_runtime().valid_ring_info(),
        )
        .map_err(crate::ops::mmff_optimization::owner_error)
    }
}
