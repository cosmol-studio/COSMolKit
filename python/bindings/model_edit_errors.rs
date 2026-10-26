//! Concrete model error dependencies; public roots expose them only through
//! their actual typed source/payload relationships.
use crate::error_framework::error_projection;
use crate::property_values::PropertyValueError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(AdjacencyError, AdjacencyErrorKind, ck::AdjacencyError, "model", {
    BondAtomOutOfRange { bond, endpoint: &'static str, atom: usize => |value: &ck::AtomId| value.index(), atom_count: usize },
    DuplicateBondId { bond, first_position: usize, second_position: usize },
    DuplicateEdge { first_bond: usize => |value: &ck::BondId| value.index(), second_bond: usize => |value: &ck::BondId| value.index(), begin: usize => |value: &ck::AtomId| value.index(), end: usize => |value: &ck::AtomId| value.index() },
}, cause = |_: &ck::AdjacencyError| None::<PyErr>, methods = {
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to AdjacencyError.BondAtomOutOfRange, AdjacencyError.DuplicateBondId. None when this payload is absent."] fn bond(&self) -> Option<usize> { match self.inner {
        ck::AdjacencyError::BondAtomOutOfRange { bond, .. }
        | ck::AdjacencyError::DuplicateBondId { bond, .. } => Some(bond.index()), _ => None,
    } }
});

error_projection!(QueryAtomConversionError, QueryAtomConversionErrorKind, ck::QueryAtomConversionError, "model", {
    NonElementAtomicNumber { atom: usize => |value: &ck::AtomId| value.index(), atomic_number: u8 },
}, cause = |_: &ck::QueryAtomConversionError| None::<PyErr>);

fn property_payload(
    py: Python<'_>,
    value: &ck::PropertyValueError,
) -> PyResult<crate::error_framework::ErrorObject<PropertyValueError>> {
    PropertyValueError::from_rust(value.clone())
        .value(py)
        .extract()
        .map_err(Into::into)
}

error_projection!(AtomPropertyError, AtomPropertyErrorKind, ck::AtomPropertyError, "property", {
    EmptyKey, ComputedListKind(_),
}, cause = |_: &ck::AtomPropertyError| None::<PyErr>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured PropertyValueError payload of AtomPropertyError.ComputedListKind; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn computed_list_kind(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyValueError>>> {
        match &self.inner { ck::AtomPropertyError::ComputedListKind(value) => property_payload(py, value).map(Some), _ => Ok(None) }
    }
});
error_projection!(BondValueError, BondValueErrorKind, ck::BondValueError, "model", {
    StereoAtomsRequired, EmptyPropertyKey, ComputedListKind(_),
}, cause = |_: &ck::BondValueError| None::<PyErr>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured PropertyValueError payload of BondValueError.ComputedListKind; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn computed_list_kind(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyValueError>>> {
        match &self.inner { ck::BondValueError::ComputedListKind(value) => property_payload(py, value).map(Some), _ => Ok(None) }
    }
});
error_projection!(MoleculePropertyError, MoleculePropertyErrorKind, ck::MoleculePropertyError, "property", {
    EmptyKey, ComputedListKind(_),
}, cause = |_: &ck::MoleculePropertyError| None::<PyErr>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured PropertyValueError payload of MoleculePropertyError.ComputedListKind; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn computed_list_kind(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyValueError>>> {
        match &self.inner { ck::MoleculePropertyError::ComputedListKind(value) => property_payload(py, value).map(Some), _ => Ok(None) }
    }
});

error_projection!(MappingValidationError, MappingValidationErrorKind, ck::MappingValidationError, "model", {
    Length { entity, direction, actual: usize, expected: usize },
    OutOfRange { entity, direction, row, mapped, target_count: usize },
    InverseMismatch { entity, direction, row, mapped },
}, cause = |_: &ck::MappingValidationError| None::<PyErr>, methods = {
    #[getter]
    #[doc = "Structural entity involved in the failure. Applies to MappingValidationError.Length, MappingValidationError.OutOfRange, MappingValidationError.InverseMismatch."] fn entity(&self) -> &'static str { match self.inner {
        ck::MappingValidationError::Length { entity, .. }
        | ck::MappingValidationError::OutOfRange { entity, .. }
        | ck::MappingValidationError::InverseMismatch { entity, .. } => entity,
    } }
    #[getter]
    #[doc = "Direction of the atom or bond mapping that failed validation."] fn direction(&self) -> &'static str { match self.inner {
        ck::MappingValidationError::Length { direction, .. }
        | ck::MappingValidationError::OutOfRange { direction, .. }
        | ck::MappingValidationError::InverseMismatch { direction, .. } => direction,
    } }
    #[getter]
    #[doc = "Row reported by the native failure. Applies to MappingValidationError.OutOfRange, MappingValidationError.InverseMismatch. None when this payload is absent."] fn row(&self) -> Option<usize> { match self.inner {
        ck::MappingValidationError::OutOfRange { row, .. }
        | ck::MappingValidationError::InverseMismatch { row, .. } => Some(row), _ => None,
    } }
    #[getter]
    #[doc = "Mapped identifier rejected during validation. Present for `MappingValidationError.OutOfRange`, `MappingValidationError.InverseMismatch`; None when unavailable."] fn mapped(&self) -> Option<usize> { match self.inner {
        ck::MappingValidationError::OutOfRange { mapped, .. }
        | ck::MappingValidationError::InverseMismatch { mapped, .. } => Some(mapped), _ => None,
    } }
});

#[derive(IntoPyObject)]
enum IndexPayload {
    Index(usize),
    Endpoints((usize, usize)),
}

error_projection!(QueryStateError, QueryStateErrorKind, ck::QueryStateError, "model", {
    NonElementAtomIdentity { position, atom: usize => |value: &ck::AtomId| value.index(), atomic_number: u8 },
    AtomCount { actual, expected }, BondCount { actual, expected },
    AtomId { position, actual, expected }, BondId { position, actual, expected },
    BondEndpoints { position, actual, expected }, Mapping(_),
    AppendedRow { entity, position }, AppendedRowCount { entity, actual, expected },
}, cause = |inner: &ck::QueryStateError| match inner {
    ck::QueryStateError::Mapping(source) => Some(MappingValidationError::from_rust(source.clone())),
    ck::QueryStateError::NonElementAtomIdentity { .. } | ck::QueryStateError::AtomCount { .. }
    | ck::QueryStateError::BondCount { .. } | ck::QueryStateError::AtomId { .. }
    | ck::QueryStateError::BondId { .. } | ck::QueryStateError::BondEndpoints { .. }
    | ck::QueryStateError::AppendedRow { .. } | ck::QueryStateError::AppendedRowCount { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<MappingValidationError>>, methods = {
    #[getter]
    #[doc = "Position reported by the native failure. Applies to QueryStateError.NonElementAtomIdentity, QueryStateError.AtomId, QueryStateError.BondId, QueryStateError.BondEndpoints, QueryStateError.AppendedRow. None when this payload is absent."] fn position(&self) -> Option<usize> { match self.inner {
        ck::QueryStateError::NonElementAtomIdentity { position, .. }
        | ck::QueryStateError::AtomId { position, .. }
        | ck::QueryStateError::BondId { position, .. }
        | ck::QueryStateError::BondEndpoints { position, .. }
        | ck::QueryStateError::AppendedRow { position, .. } => Some(position), _ => None,
    } }
    #[getter]
    #[doc = "Structural entity involved in the failure. Applies to QueryStateError.AppendedRow, QueryStateError.AppendedRowCount. None when this payload is absent."] fn entity(&self) -> Option<&'static str> { match self.inner {
        ck::QueryStateError::AppendedRow { entity, .. }
        | ck::QueryStateError::AppendedRowCount { entity, .. } => Some(entity), _ => None,
    } }
    #[getter]
    #[doc = "Observed value reported by the native check. Applies to QueryStateError.AtomCount, QueryStateError.BondCount, QueryStateError.AppendedRowCount, QueryStateError.AtomId, QueryStateError.BondId, QueryStateError.BondEndpoints. None when this payload is absent."] fn actual(&self) -> Option<IndexPayload> { match self.inner {
        ck::QueryStateError::AtomCount { actual, .. }
        | ck::QueryStateError::BondCount { actual, .. }
        | ck::QueryStateError::AppendedRowCount { actual, .. } => Some(IndexPayload::Index(actual)),
        ck::QueryStateError::AtomId { actual, .. } => Some(IndexPayload::Index(actual.index())),
        ck::QueryStateError::BondId { actual, .. } => Some(IndexPayload::Index(actual.index())),
        ck::QueryStateError::BondEndpoints { actual, .. } => Some(IndexPayload::Endpoints((actual.0.index(),actual.1.index()))),
        _ => None,
    } }
    #[getter]
    #[doc = "Required value reported by the native check. Applies to QueryStateError.AtomCount, QueryStateError.BondCount, QueryStateError.AppendedRowCount, QueryStateError.AtomId, QueryStateError.BondId, QueryStateError.BondEndpoints. None when this payload is absent."] fn expected(&self) -> Option<IndexPayload> { match self.inner {
        ck::QueryStateError::AtomCount { expected, .. }
        | ck::QueryStateError::BondCount { expected, .. }
        | ck::QueryStateError::AppendedRowCount { expected, .. } => Some(IndexPayload::Index(expected)),
        ck::QueryStateError::AtomId { expected, .. } => Some(IndexPayload::Index(expected.index())),
        ck::QueryStateError::BondId { expected, .. } => Some(IndexPayload::Index(expected.index())),
        ck::QueryStateError::BondEndpoints { expected, .. } => Some(IndexPayload::Endpoints((expected.0.index(),expected.1.index()))),
        _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_adjacency_and_non_element_failures_keep_exact_model_fields() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let bond = ck::Bond::from_spec(
                ck::BondId::new(9),
                ck::BondSpec::new(
                    ck::AtomId::new(0),
                    ck::AtomId::new(5),
                    ck::BondOrder::Single,
                ),
            );
            let native = ck::AdjacencyList::try_from_topology(2, &[bond]).unwrap_err();
            let display = native.to_string();
            assert!(native.source().is_none());
            let error = AdjacencyError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(error.value(py).getattr("bond")?.extract::<usize>()?, 9);
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 5);
            assert_eq!(
                error.value(py).getattr("endpoint")?.extract::<String>()?,
                "end"
            );
            assert_eq!(
                error.value(py).getattr("atom_count")?.extract::<usize>()?,
                2
            );
            assert!(error.cause(py).is_none());
            assert!(
                error
                    .value(py)
                    .setattr("atom", 0)
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );

            let query = ck::QueryAtom::from_identity_parts(
                ck::AtomId::new(17),
                ck::QueryAtomIdentity::from_atomic_number(255),
                ck::QueryNode::predicate(ck::AtomQueryPredicate::AtomicNumber(255)),
            );
            let native = query.try_to_atom().unwrap_err();
            let display = native.to_string();
            assert!(native.source().is_none());
            let error = QueryAtomConversionError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 17);
            assert_eq!(
                error.value(py).getattr("atomic_number")?.extract::<u8>()?,
                255
            );
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn model_error_dependencies_keep_real_payloads_and_source_shapes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let mut atom =
                ck::Atom::from_spec(ck::AtomId::new(0), ck::AtomSpec::new(ck::Element::C));
            let native = atom.set_prop("", ck::PropertyValue::Int(1)).unwrap_err();
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = AtomPropertyError::from_rust(native);
            assert!(error.cause(py).is_none());
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let native = ck::MoleculeProperties::default()
                .with_prop("", 1_i32)
                .unwrap_err();
            let error = MoleculePropertyError::from_rust(native);
            assert!(error.value(py).getattr("computed_list_kind")?.is_none());
            let value_error = ck::PropertyValue::Int(7).as_bool().unwrap_err();
            for error in [
                AtomPropertyError::from_rust(ck::AtomPropertyError::ComputedListKind(
                    value_error.clone(),
                )),
                BondValueError::from_rust(ck::BondValueError::ComputedListKind(
                    value_error.clone(),
                )),
                MoleculePropertyError::from_rust(ck::MoleculePropertyError::ComputedListKind(
                    value_error.clone(),
                )),
            ] {
                assert!(error.cause(py).is_none());
                let payload = error.value(py).getattr("computed_list_kind")?;
                assert!(payload.is_instance(<PropertyValueError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                assert_eq!(payload.str()?.to_str()?, value_error.to_string());
                assert_eq!(
                    error
                        .value(py)
                        .getattr("source_message")?
                        .extract::<Option<String>>()?,
                    None
                );
            }
            let native = ck::MappingValidationError::Length {
                entity: "atom",
                direction: "old-to-new",
                actual: 2,
                expected: 3,
            };
            let display = native.to_string();
            let native = ck::QueryStateError::Mapping(native);
            assert!(native.source().unwrap().is::<ck::MappingValidationError>());
            let error = QueryStateError::from_rust(native);
            let cause = error.cause(py).unwrap();
            assert!(cause.is_instance(py, <MappingValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(cause.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            assert_eq!(cause.value(py).getattr("actual")?.extract::<usize>()?, 2);
            assert_eq!(cause.value(py).getattr("expected")?.extract::<usize>()?, 3);
            let error = QueryStateError::from_rust(ck::QueryStateError::BondEndpoints {
                position: 4,
                actual: (ck::AtomId::new(2), ck::AtomId::new(3)),
                expected: (ck::AtomId::new(1), ck::AtomId::new(0)),
            });
            assert_eq!(
                error
                    .value(py)
                    .getattr("actual")?
                    .extract::<(usize, usize)>()?,
                (2, 3)
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("expected")?
                    .extract::<(usize, usize)>()?,
                (1, 0)
            );
            let error = QueryStateError::from_rust(ck::QueryStateError::AtomId {
                position: 1,
                actual: ck::AtomId::new(7),
                expected: ck::AtomId::new(1),
            });
            assert_eq!(error.value(py).getattr("actual")?.extract::<usize>()?, 7);
            assert_eq!(error.value(py).getattr("expected")?.extract::<usize>()?, 1);
            Ok(())
        })
        .unwrap();
    }
}
