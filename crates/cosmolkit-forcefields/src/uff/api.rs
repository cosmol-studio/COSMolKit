use std::error::Error;
use std::fmt;

use super::atom_typer::{UffAtomStateRef, UffTypingError};
use super::builder::UffBuilderError;
use super::params::UffParamError;

#[cfg(test)]
std::thread_local! {
    static PREPARE_PARAMETER_QUERY_CALLS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
}

#[cfg(test)]
pub(super) fn reset_prepare_parameter_query_calls() {
    PREPARE_PARAMETER_QUERY_CALLS.with(|calls| calls.set(0));
}

#[cfg(test)]
pub(super) fn prepare_parameter_query_calls() -> usize {
    PREPARE_PARAMETER_QUERY_CALLS.with(std::cell::Cell::get)
}

#[derive(Debug)]
pub struct UffParameterError {
    cause: UffParameterCause,
}

#[derive(Debug)]
enum UffParameterCause {
    Preparation(UffBuilderError),
    ParameterTable(UffParamError),
    Typing(UffTypingError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UffParameterErrorKind {
    Preparation,
    ParameterTable,
    Typing,
}

impl UffParameterError {
    pub fn kind(&self) -> UffParameterErrorKind {
        match &self.cause {
            UffParameterCause::Preparation(_) => UffParameterErrorKind::Preparation,
            UffParameterCause::ParameterTable(_) => UffParameterErrorKind::ParameterTable,
            UffParameterCause::Typing(_) => UffParameterErrorKind::Typing,
        }
    }
}

impl fmt::Display for UffParameterError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.cause {
            UffParameterCause::Preparation(source) => fmt::Display::fmt(source, formatter),
            UffParameterCause::ParameterTable(source) => fmt::Display::fmt(source, formatter),
            UffParameterCause::Typing(source) => fmt::Display::fmt(source, formatter),
        }
    }
}

impl Error for UffParameterError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.cause {
            UffParameterCause::Preparation(source) => Some(source),
            UffParameterCause::ParameterTable(source) => Some(source),
            UffParameterCause::Typing(source) => Some(source),
        }
    }
}

pub(super) struct PreparedParameterQuery<'a> {
    pub(super) typing_state: UffAtomStateRef<'a>,
}

pub(super) fn prepare_parameter_query<'a>(
    topology: &'a cosmolkit_model::TopologyBlock,
    assignment: &'a cosmolkit_core::ValenceAssignment,
) -> Result<PreparedParameterQuery<'a>, UffParameterError> {
    #[cfg(test)]
    PREPARE_PARAMETER_QUERY_CALLS.with(|calls| calls.set(calls.get() + 1));

    // RDKit❗✔️: int totalValence = atom->getTotalValence();
    // RDKit❗✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit❗✔️: for (const auto bnd : mol.atomBonds(at)) {
    // RDKit❗✔️:   if (bnd->getIsConjugated()) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: return false;
    // The borrowed cached state retains validated row lengths, value ranges,
    // no_implicit handling, and topology-validation order. The shared typer
    // reads total valence by row and follows incident bonds only when the
    // source label branch asks for conjugation; neither projection is stored.
    let typing_state =
        UffAtomStateRef::cached(topology, assignment).map_err(|source| UffParameterError {
            cause: UffParameterCause::Preparation(source),
        })?;

    Ok(PreparedParameterQuery { typing_state })
}

pub fn uff_has_all_molecule_params(
    topology: &cosmolkit_model::TopologyBlock,
    assignment: &cosmolkit_core::ValenceAssignment,
) -> Result<bool, UffParameterError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFFHasAllMoleculeParams (rdForceFields.cpp:107-112)
    // RDKit❗❌: bool UFFHasAllMoleculeParams(const ROMol &mol) {
    // RDKit❗❌:   UFF::AtomicParamVect types;
    // RDKit❗❌:   bool foundAll;
    // RDKit❗❌:   boost::tie(types, foundAll) = UFF::getAtomTypes(mol);
    // RDKit❗❌:   return foundAll;
    // RDKit❗❌: }
    // The source typer obtains its cached default table before traversing atom
    // rows. Preserve that boundary order before validating the supplied cache.
    let params =
        super::params::ParamCollection::get_params("").map_err(|source| UffParameterError {
            cause: UffParameterCause::ParameterTable(source),
        })?;

    let prepared = prepare_parameter_query(topology, assignment)?;
    let mut diagnostics = Vec::new();
    let (_, found_all) = super::atom_typer::get_atom_types_from_state(
        topology,
        prepared.typing_state,
        params.as_ref(),
        &mut diagnostics,
    )
    .map_err(|source| UffParameterError {
        cause: UffParameterCause::Typing(source),
    })?;
    // The shared typer materializes the source N nullable parameter slots and
    // appends private diagnostics in atom order. Cached preparation no longer
    // allocates total-valence or conjugation projections; topology validation
    // still builds temporary adjacency, and the typer/diagnostics still own
    // their respective vectors. No zero-allocation claim is made.
    Ok(found_all)
    // END RDKIT CPP FUNCTION RDKit::UFFHasAllMoleculeParams
}

/// Error from the narrow UFF rest-length projection used by distance geometry.
#[derive(Debug)]
pub struct UffBoundsError {
    cause: UffBoundsCause,
}
#[derive(Debug)]
enum UffBoundsCause {
    Parameter(UffParameterError),
    BondOrder(cosmolkit_core::ValenceError),
    BondMath(super::bond::BondMathError),
    AssignmentLength {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
}
impl fmt::Display for UffBoundsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.cause {
            UffBoundsCause::Parameter(e) => fmt::Display::fmt(e, f),
            UffBoundsCause::BondOrder(e) => fmt::Display::fmt(e, f),
            UffBoundsCause::BondMath(e) => fmt::Display::fmt(e, f),
            UffBoundsCause::AssignmentLength {
                field,
                actual,
                expected,
            } => write!(
                f,
                "{field} assignment has {actual} rows; expected {expected}"
            ),
        }
    }
}
impl Error for UffBoundsError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.cause {
            UffBoundsCause::Parameter(e) => Some(e),
            UffBoundsCause::BondOrder(e) => Some(e),
            UffBoundsCause::BondMath(e) => Some(e),
            UffBoundsCause::AssignmentLength { .. } => None,
        }
    }
}
/// Ordered source UFF bond rest lengths, borrowing independently computed atom assignments.
/// `None` is exactly the set12Bounds branch for missing endpoint parameters or nonpositive order.
/// The caller owns the source-defined van der Waals fallback and distance-bound policy.
pub fn uff_bond_rest_lengths(
    topology: &cosmolkit_model::TopologyBlock,
    assignment: &cosmolkit_core::ValenceAssignment,
    hybridizations: &[cosmolkit_model::Hybridization],
    conjugated_bonds: &[bool],
) -> Result<Vec<Option<f64>>, UffBoundsError> {
    // RDKit❗✔️:   auto [atomParams, foundAll] = UFF::getAtomTypes(mol);
    // RDKit❗✔️:   CHECK_INVARIANT(atomParams.size() == mol.getNumAtoms(),
    // RDKit❗✔️:                   "parameter vector size mismatch");
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     auto begId = bond->getBeginAtomIdx();
    // RDKit❗✔️:     auto endId = bond->getEndAtomIdx();
    // RDKit❗✔️:     auto bOrder = bond->getBondTypeAsDouble();
    // RDKit❗✔️:     if (atomParams[begId] && atomParams[endId] && bOrder > 0) {
    // RDKit❗✔️:       auto bl = ForceFields::UFF::Utils::calcBondRestLength(
    // RDKit❗✔️:           bOrder, atomParams[begId], atomParams[endId]);
    // Source: pinned BoundsMatrixBuilder.cpp:244-267. The caller retains all
    // squish/bound/fallback writes; this projection borrows one existing typer
    // and one existing formula. No topology, parameter table, or assignment clone.
    // Additional ordered E optional values are the cross-owner output cost.
    let params =
        super::params::ParamCollection::get_params("").map_err(|source| UffBoundsError {
            cause: UffBoundsCause::Parameter(UffParameterError {
                cause: UffParameterCause::ParameterTable(source),
            }),
        })?;
    let prepared =
        prepare_parameter_query(topology, assignment).map_err(|source| UffBoundsError {
            cause: UffBoundsCause::Parameter(source),
        })?;
    for (field, actual, expected) in [
        ("hybridization", hybridizations.len(), topology.atoms.len()),
        ("conjugation", conjugated_bonds.len(), topology.bonds.len()),
    ] {
        if actual != expected {
            return Err(UffBoundsError {
                cause: UffBoundsCause::AssignmentLength {
                    field,
                    actual,
                    expected,
                },
            });
        }
    }
    let mut diagnostics = Vec::new();
    let (atom_params, _found_all) = super::atom_typer::get_atom_types_with_assignments(
        topology,
        prepared.typing_state,
        &params,
        &mut diagnostics,
        Some((hybridizations, conjugated_bonds)),
    )
    .map_err(|source| UffBoundsError {
        cause: UffBoundsCause::Parameter(UffParameterError {
            cause: UffParameterCause::Typing(source),
        }),
    })?;
    let mut lengths = Vec::with_capacity(topology.bonds.len());
    for bond in &topology.bonds {
        let beg = bond.begin().index();
        let end = bond.end().index();
        let order =
            cosmolkit_core::bond_type_as_double(bond.order()).map_err(|source| UffBoundsError {
                cause: UffBoundsCause::BondOrder(source),
            })?;
        let length = match (atom_params[beg], atom_params[end]) {
            (Some(p1), Some(p2)) if order > 0.0 => Some(
                super::bond::calc_bond_rest_length(order, p1, p2).map_err(|source| {
                    UffBoundsError {
                        cause: UffBoundsCause::BondMath(source),
                    }
                })?,
            ),
            _ => None,
        };
        lengths.push(length);
    }
    Ok(lengths)
}

#[cfg(test)]
mod tests {
    use std::error::Error;
    use std::sync::Arc;

    use cosmolkit_core::ValenceAssignment;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element, Hybridization,
        TopologyBlock,
    };

    use super::{
        PreparedParameterQuery, UffParameterCause, UffParameterError, UffParameterErrorKind,
        prepare_parameter_query_calls, reset_prepare_parameter_query_calls,
    };
    use crate::uff::atom_typer::{self, UffAtomStateRef, UffTypingError, UffTypingInput};
    use crate::uff::builder::{
        self, PreparedValenceField, UffBuilderError, conjugation_projection_vec_constructions,
        reset_conjugation_projection_vec_constructions,
        reset_typing_valence_projection_vec_constructions,
        typing_valence_projection_vec_constructions,
    };
    use crate::uff::params::{ParamCollection, UffParamError};

    fn assert_source_is_stored(error: &UffParameterError) -> &(dyn Error + 'static) {
        let reported = Error::source(error).expect("parameter errors retain their typed cause");
        let stored = match &error.cause {
            UffParameterCause::Preparation(source) => source as &(dyn Error + 'static),
            UffParameterCause::ParameterTable(source) => source as &(dyn Error + 'static),
            UffParameterCause::Typing(source) => source as &(dyn Error + 'static),
        };
        assert!(std::ptr::eq(reported, stored));
        reported
    }

    fn one_carbon_topology() -> TopologyBlock {
        let carbon = Element::from_atomic_number(6).expect("carbon is in the element table");
        let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(carbon));
        TopologyBlock::try_from_parts(vec![atom], Vec::new(), Vec::new(), Vec::new())
            .expect("one isolated carbon is a valid topology")
    }

    fn empty_topology() -> TopologyBlock {
        TopologyBlock::try_from_parts(Vec::new(), Vec::new(), Vec::new(), Vec::new())
            .expect("empty topology is valid")
    }

    fn topology_with_no_implicit(no_implicit: &[bool]) -> TopologyBlock {
        let carbon = Element::from_atomic_number(6).expect("carbon is in the element table");
        let atoms = no_implicit
            .iter()
            .enumerate()
            .map(|(row, &no_implicit)| {
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(carbon).with_no_implicit(no_implicit),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed atom-only topology is valid")
    }

    fn topology_with_bonds(atom_count: usize, edges: &[(usize, usize, bool)]) -> TopologyBlock {
        let carbon = Element::from_atomic_number(6).expect("carbon is in the element table");
        let atoms = (0..atom_count)
            .map(|row| Atom::from_spec(AtomId::new(row), AtomSpec::new(carbon)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end, is_conjugated))| {
                let mut bond = Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                );
                bond.set_conjugated(is_conjugated);
                bond
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed bond topology is structurally valid")
    }

    fn label_atom(
        atomic_number: u8,
        id: usize,
        hybridization: Hybridization,
        dummy_label: Option<&str>,
    ) -> Atom {
        let element = Element::from_atomic_number(atomic_number)
            .expect("fixed label test atomic number is in the model range");
        let spec = AtomSpec::new(element).with_hybridization(hybridization);
        let spec = if let Some(dummy_label) = dummy_label {
            spec.with_prop("dummyLabel", dummy_label)
                .expect("dummy label property has a nonempty key")
        } else {
            spec
        };
        Atom::from_spec(AtomId::new(id), spec)
    }

    fn topology_from_atoms(atoms: Vec<Atom>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed atom-only topology is structurally valid")
    }

    fn assignment(explicit_valence: &[i32], implicit_hydrogens: &[i32]) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit_valence.to_vec(),
            implicit_hydrogens: implicit_hydrogens.to_vec(),
        }
    }

    fn assert_cached_rows(
        prepared: &PreparedParameterQuery<'_>,
        topology: &TopologyBlock,
        assignment: &ValenceAssignment,
        expected_total_valences: &[i32],
        expected_conjugated_presence: &[bool],
    ) {
        let UffAtomStateRef::Cached {
            topology: cached_topology,
            assignment: cached_assignment,
        } = prepared.typing_state
        else {
            panic!("prepared query retains the original cached inputs");
        };
        assert!(std::ptr::eq(cached_topology, topology));
        assert!(std::ptr::eq(cached_assignment, assignment));
        assert_eq!(expected_total_valences.len(), topology.atoms.len());
        assert_eq!(expected_conjugated_presence.len(), topology.atoms.len());

        for (atom_index, &expected) in expected_total_valences.iter().enumerate() {
            assert_eq!(prepared.typing_state.total_valence_at(atom_index), expected);
        }
        for (atom_index, &expected) in expected_conjugated_presence.iter().enumerate() {
            assert_eq!(
                prepared.typing_state.conjugated_presence_at(atom_index),
                expected
            );
        }
    }

    fn assert_preparation_source(error: &UffParameterError, expected: UffBuilderError) {
        assert_eq!(error.kind(), UffParameterErrorKind::Preparation);
        let source = Error::source(error).expect("preparation errors retain their typed cause");
        assert_eq!(
            source.downcast_ref::<UffBuilderError>(),
            Some(&expected),
            "preparation error must preserve its exact cause payload"
        );
    }

    #[test]
    fn uff_param_p01_preparation_cause_borrows_actual_length_failure() {
        let topology = one_carbon_topology();
        let assignment = ValenceAssignment {
            explicit_valence: Vec::new(),
            implicit_hydrogens: Vec::new(),
        };
        let actual = builder::prepare_typing_valence(&topology, &assignment)
            .expect_err("the existing preparation helper reports its explicit row mismatch");
        let display = actual.to_string();
        let error = UffParameterError {
            cause: UffParameterCause::Preparation(actual),
        };

        assert_eq!(error.kind(), UffParameterErrorKind::Preparation);
        assert_eq!(error.to_string(), display);
        let source = assert_source_is_stored(&error);
        let downcast = source
            .downcast_ref::<UffBuilderError>()
            .expect("the preparation cause keeps its concrete type");
        assert_eq!(
            downcast,
            &UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::Explicit,
                expected: 1,
                actual: 0,
            }
        );
    }

    #[test]
    fn uff_param_p01_parameter_table_cause_wraps_custom_parse_failure() {
        // This real custom-data parse failure exercises the immutable table
        // wrapper branch; the built-in default table itself is valid and is
        // not represented as failing in the public query.
        let actual = ParamCollection::get_params("\n")
            .expect_err("an empty terminated record is rejected by the parser");
        let display = actual.to_string();
        let error = UffParameterError {
            cause: UffParameterCause::ParameterTable(actual),
        };

        assert_eq!(error.kind(), UffParameterErrorKind::ParameterTable);
        assert_eq!(error.to_string(), display);
        let source = assert_source_is_stored(&error);
        let downcast = source
            .downcast_ref::<UffParamError>()
            .expect("the parameter-table cause keeps its concrete type");
        assert_eq!(downcast, &UffParamError::EmptyLine { line_number: 1 });
    }

    #[test]
    fn uff_param_p01_typing_cause_borrows_actual_default_table_typing_failure() {
        let topology = empty_topology();
        let params = ParamCollection::get_params("").expect("the default table is valid");
        let actual = atom_typer::uff_has_all_molecule_parameters(
            &topology,
            &[0],
            &[],
            params.as_ref(),
            &mut Vec::new(),
        )
        .expect_err("the existing typer rejects the mismatched prepared-state length");
        let display = actual.to_string();
        let error = UffParameterError {
            cause: UffParameterCause::Typing(actual),
        };

        assert_eq!(error.kind(), UffParameterErrorKind::Typing);
        assert_eq!(error.to_string(), display);
        let source = assert_source_is_stored(&error);
        let downcast = source
            .downcast_ref::<UffTypingError>()
            .expect("the typing cause keeps its concrete type");
        assert_eq!(
            downcast,
            &UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 0,
                actual: 1,
            }
        );
    }

    #[test]
    fn uff_param_p02_empty_and_isolated_b01_rows_preserve_borrowed_inputs() {
        let empty = empty_topology();
        let empty_assignment = assignment(&[], &[]);
        let prepared = super::prepare_parameter_query(&empty, &empty_assignment)
            .expect("empty B01 rows prepare without atom-cache reads");
        assert_cached_rows(&prepared, &empty, &empty_assignment, &[], &[]);

        // Retain B01's ordered no_implicit fixture, including its stored -1
        // implicit sentinel on the row whose source flag suppresses that read.
        let isolated = topology_with_no_implicit(&[false, true, false]);
        let cached = assignment(&[4, 3, 1], &[1, -1, 2]);
        let prepared = super::prepare_parameter_query(&isolated, &cached)
            .expect("B01 no_implicit rows ignore only the suppressed sentinel");
        assert_cached_rows(
            &prepared,
            &isolated,
            &cached,
            &[5, 3, 3],
            &[false, false, false],
        );
    }

    #[test]
    fn uff_param_p02_disconnected_b02_conjugated_rows_keep_atom_order() {
        // Retain B02's disconnected edges and its one conjugated endpoint pair.
        let topology = topology_with_bonds(4, &[(0, 2, true), (3, 1, false)]);
        let cached = assignment(&[1, 1, 1, 1], &[3, 3, 3, 3]);
        let prepared = super::prepare_parameter_query(&topology, &cached)
            .expect("source-aligned valence and conjugation rows prepare");

        assert_cached_rows(
            &prepared,
            &topology,
            &cached,
            &[4, 4, 4, 4],
            &[true, false, true, false],
        );
    }

    #[test]
    fn uff_integrate_i01_prepared_identity_and_cached_facts() {
        let empty = empty_topology();
        let empty_assignment = assignment(&[], &[]);
        let prepared = super::prepare_parameter_query(&empty, &empty_assignment)
            .expect("empty source state keeps both original borrows");
        assert_cached_rows(&prepared, &empty, &empty_assignment, &[], &[]);

        let disconnected = topology_with_bonds(4, &[(0, 2, true), (3, 1, false)]);
        let cached = assignment(&[1, 1, 1, 1], &[3, 3, 3, 3]);
        let prepared = super::prepare_parameter_query(&disconnected, &cached)
            .expect("disconnected source state is projected in atom-row order");
        assert_cached_rows(
            &prepared,
            &disconnected,
            &cached,
            &[4, 4, 4, 4],
            &[true, false, true, false],
        );
    }

    #[test]
    fn uff_prepare_p06_actual_query_keeps_borrowed_state_and_nullable_rows() {
        let topology = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp3, None),
            label_atom(0, 1, Hybridization::Unspecified, Some("unknown-middle")),
            label_atom(6, 2, Hybridization::Sp3, None),
        ]);
        let cached = assignment(&[4, 0, 4], &[0, 0, 0]);
        let prepared = super::prepare_parameter_query(&topology, &cached)
            .expect("valid cached values are borrowed directly");
        assert_cached_rows(
            &prepared,
            &topology,
            &cached,
            &[4, 0, 4],
            &[false, false, false],
        );

        let params = ParamCollection::get_params("").expect("the pinned default table is valid");
        let (slots, found_all) = atom_typer::get_atom_types_from_state(
            &topology,
            prepared.typing_state,
            params.as_ref(),
            &mut Vec::new(),
        )
        .expect("all three source rows are visited");
        assert_eq!(slots.len(), 3);
        assert!(slots[0].is_some());
        assert!(slots[1].is_none());
        assert!(slots[2].is_some());
        assert!(!found_all);

        reset_prepare_parameter_query_calls();
        reset_typing_valence_projection_vec_constructions();
        reset_conjugation_projection_vec_constructions();
        assert!(
            !super::uff_has_all_molecule_params(&topology, &cached)
                .expect("an unrecognized middle row remains a false query result")
        );
        assert_eq!(prepare_parameter_query_calls(), 1);
        assert_eq!(typing_valence_projection_vec_constructions(), 0);
        assert_eq!(conjugation_projection_vec_constructions(), 0);
    }

    #[test]
    fn uff_prepare_p06_actual_query_does_not_recompute_missing_cache_rows() {
        let topology = one_carbon_topology();
        let cases = [
            (
                assignment(&[], &[0]),
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::Explicit,
                    expected: 1,
                    actual: 0,
                },
            ),
            (
                assignment(&[0], &[]),
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::ImplicitHydrogen,
                    expected: 1,
                    actual: 0,
                },
            ),
        ];

        for (cached, expected) in cases {
            reset_prepare_parameter_query_calls();
            let error = super::uff_has_all_molecule_params(&topology, &cached)
                .expect_err("missing stored cache rows remain typed errors");
            assert_preparation_source(&error, expected);
            assert!(
                assert_source_is_stored(&error)
                    .downcast_ref::<UffBuilderError>()
                    .is_some()
            );
            assert_eq!(prepare_parameter_query_calls(), 1);
        }
    }

    #[test]
    fn uff_param_p02_b01_length_mismatches_preserve_field_and_direction() {
        let topology = one_carbon_topology();
        let cases = [
            (
                assignment(&[], &[0]),
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::Explicit,
                    expected: 1,
                    actual: 0,
                },
            ),
            (
                assignment(&[0, 1], &[0]),
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::Explicit,
                    expected: 1,
                    actual: 2,
                },
            ),
            (
                assignment(&[0], &[]),
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::ImplicitHydrogen,
                    expected: 1,
                    actual: 0,
                },
            ),
            (
                assignment(&[0], &[0, 1]),
                UffBuilderError::ValenceAssignmentLengthMismatch {
                    field: PreparedValenceField::ImplicitHydrogen,
                    expected: 1,
                    actual: 2,
                },
            ),
        ];

        for (cached, expected) in cases {
            let error = super::prepare_parameter_query(&topology, &cached)
                .err()
                .expect("both source cache lengths must match the atom rows");
            assert_preparation_source(&error, expected);
        }
    }

    #[test]
    fn uff_param_p02_b01_negative_and_out_of_byte_cache_rows_remain_typed() {
        let topology = one_carbon_topology();
        let cases = [
            (
                assignment(&[-1], &[0]),
                UffBuilderError::SourceValencePrecondition {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::Explicit,
                    value: -1,
                },
            ),
            (
                assignment(&[128], &[0]),
                UffBuilderError::SourceValenceOutOfRange {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::Explicit,
                    value: 128,
                },
            ),
            (
                assignment(&[0], &[-1]),
                UffBuilderError::SourceValencePrecondition {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::ImplicitHydrogen,
                    value: -1,
                },
            ),
            (
                assignment(&[0], &[128]),
                UffBuilderError::SourceValenceOutOfRange {
                    atom_id: AtomId::new(0),
                    field: PreparedValenceField::ImplicitHydrogen,
                    value: 128,
                },
            ),
        ];

        for (cached, expected) in cases {
            let error = super::prepare_parameter_query(&topology, &cached)
                .err()
                .expect("invalid stored source-cache components must not be normalized");
            assert_preparation_source(&error, expected);
        }
    }

    #[test]
    fn uff_param_p03_empty_and_known_rows_reuse_cached_default_table() {
        let cached_params = ParamCollection::get_params("").expect("pinned default UFF table");
        let empty = empty_topology();
        assert!(
            super::uff_has_all_molecule_params(&empty, &assignment(&[], &[]))
                .expect("empty topology has no atom reads")
        );

        let known = topology_from_atoms(vec![
            label_atom(6, 0, Hybridization::Sp3, None),
            label_atom(6, 1, Hybridization::Sp3, None),
        ]);
        let known_cache = assignment(&[4, 4], &[0, 0]);
        assert!(
            super::uff_has_all_molecule_params(&known, &known_cache)
                .expect("both carbon labels are present in the default table")
        );
        assert!(Arc::ptr_eq(
            &cached_params,
            &ParamCollection::get_params("").expect("cached default UFF table remains available"),
        ));

        assert!(
            super::uff_has_all_molecule_params(&known, &known_cache)
                .expect("repeated query reuses the same valid table")
        );
        assert!(Arc::ptr_eq(
            &cached_params,
            &ParamCollection::get_params("").expect("repeated query keeps the cached table"),
        ));
    }

    #[test]
    fn uff_param_p03_unknown_first_middle_last_rows_return_false() {
        // W09 remains the independent assertion of nullable slots and ordered
        // diagnostics while this test exercises the detached bool entry.
        let topology = topology_from_atoms(vec![
            label_atom(0, 0, Hybridization::Unspecified, Some("unknown-first")),
            label_atom(6, 1, Hybridization::Sp3, None),
            label_atom(0, 2, Hybridization::Unspecified, Some("unknown-middle")),
            label_atom(6, 3, Hybridization::Sp3, None),
            label_atom(0, 4, Hybridization::Unspecified, Some("unknown-last")),
        ]);
        let cached = assignment(&[0, 4, 0, 4, 0], &[0; 5]);

        assert!(
            !super::uff_has_all_molecule_params(&topology, &cached)
                .expect("missing rows produce false after source traversal")
        );
    }

    #[test]
    fn uff_param_p03_actual_preparation_failure_keeps_typed_source_chain() {
        let topology = one_carbon_topology();
        let invalid_cache = assignment(&[], &[0]);
        let error = super::uff_has_all_molecule_params(&topology, &invalid_cache)
            .expect_err("a missing explicit cache row is a typed precondition failure");
        assert_preparation_source(
            &error,
            UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::Explicit,
                expected: 1,
                actual: 0,
            },
        );
        let source = assert_source_is_stored(&error);
        assert_eq!(
            source
                .downcast_ref::<UffBuilderError>()
                .expect("actual preparation source remains concrete"),
            &UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::Explicit,
                expected: 1,
                actual: 0,
            },
        );
    }
}

#[cfg(test)]
mod conformer_bond_length_tests {
    use super::*;
    use cosmolkit_core::{
        ValenceModel, assign_conjugation_flags, assign_hybridization_with_conjugation,
        assign_valence_for_topology,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element, TopologyBlock,
    };
    fn lengths(graph: &TopologyBlock) -> Vec<Option<f64>> {
        let valence = assign_valence_for_topology(graph, ValenceModel::RdkitLike).unwrap();
        let conjugated = assign_conjugation_flags(graph, &valence).unwrap();
        let hybridization =
            assign_hybridization_with_conjugation(graph, &valence, &conjugated).unwrap();
        uff_bond_rest_lengths(graph, &valence, &hybridization.values, &conjugated).unwrap()
    }
    #[test]
    fn original_explicit_sulfur_hydrogen_rest_lengths() {
        let record = cosmolkit_smiles::parse_smiles("C[SH+][O-]", &Default::default()).unwrap();
        let result = lengths(&record.topology);
        assert!((result[0].unwrap() - 1.776_838_813_449_356_2).abs() < 1.0e-14);
        assert!((result[1].unwrap() - 1.679_472_654_030_108).abs() < 1.0e-14);
    }
    #[test]
    fn original_unsanitized_ethane_computed_assignment_is_used() {
        let record = cosmolkit_smiles::parse_smiles(
            "CC",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let stale = record
            .topology
            .atoms
            .iter()
            .map(|a| a.hybridization())
            .collect::<Vec<_>>();
        let result = lengths(&record.topology);
        assert!(result[0].unwrap() > 1.4);
        assert!(result[0].unwrap() < 1.6);
        assert_eq!(
            stale,
            record
                .topology
                .atoms
                .iter()
                .map(|a| a.hybridization())
                .collect::<Vec<_>>()
        );
    }
    #[test]
    fn original_dummy_endpoint_has_no_uff_parameter() {
        let atoms = [Element::C, Element::DUMMY]
            .into_iter()
            .enumerate()
            .map(|(i, e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
            .collect();
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        assert_eq!(lengths(&graph), vec![None]);
    }
    #[test]
    fn detached_assignment_shape_error_keeps_typed_preparation_order() {
        let record = cosmolkit_smiles::parse_smiles("CC", &Default::default()).unwrap();
        let graph = &record.topology;
        let valid = assign_valence_for_topology(graph, ValenceModel::RdkitLike).unwrap();
        let error = uff_bond_rest_lengths(graph, &valid, &[], &[false]).unwrap_err();
        assert!(matches!(
            error.cause,
            UffBoundsCause::AssignmentLength {
                field: "hybridization",
                actual: 0,
                expected: 2
            }
        ));
        let invalid = cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![],
            implicit_hydrogens: vec![],
        };
        let error = uff_bond_rest_lengths(graph, &invalid, &[], &[]).unwrap_err();
        assert!(
            Error::source(&error)
                .unwrap()
                .downcast_ref::<UffParameterError>()
                .is_some()
        );
    }
}

/// Typed source hydrogen-getter failure retained behind the existing FF owner.
#[derive(Debug)]
pub struct MissingExplicitHydrogensError {
    source: UffBuilderError,
}
impl fmt::Display for MissingExplicitHydrogensError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.source, f)
    }
}
impl Error for MissingExplicitHydrogensError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.source)
    }
}
/// Reuse the existing source needsHs owner without emitting a UFF warning.
#[doc(hidden)]
pub fn needs_explicit_hydrogens(
    topology: &cosmolkit_model::TopologyBlock,
    implicit_hydrogens: &[i32],
) -> Result<bool, MissingExplicitHydrogensError> {
    super::builder::needs_hydrogens(topology, implicit_hydrogens)
        .map_err(|source| MissingExplicitHydrogensError { source })
}
