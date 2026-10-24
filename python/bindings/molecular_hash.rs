//! Native molecular hashing and its closed, concrete error graph.
use crate::chemical_values::BondOrder;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{MappingValidationError, QueryStateError};
use crate::model_errors::TopologyValidationError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(CipRankError, CipRankErrorKind, ck::CipRankError, "stereo", {
    InvalidTopology(_), InvalidQueryState(_),
    ValenceRowCount { field: &'static str, actual, atom_count },
    NegativeImplicitHydrogen { atom, value },
    AtomMapOutOfRange { atom, map_number: u32 },
    InvariantCount { actual, atom_count },
    InvariantOutOfRange { atom, value },
    TooManyNeighbors { atom, degree: usize, maximum_supported: usize },
    UnsupportedBondOrder { bond: usize => |value: &ck::BondId| value.index(), order: BondOrder => |value: &ck::BondOrder| (*value).into() },
}, cause = |inner: &ck::CipRankError| Python::attach(|py| match inner {
    // The query-state wrapper is transparent: retain its *source*, rather
    // than adding the wrapped QueryStateError as an invented cause.
    ck::CipRankError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).cause(py),
    ck::CipRankError::InvalidTopology(_) | ck::CipRankError::ValenceRowCount { .. }
    | ck::CipRankError::NegativeImplicitHydrogen { .. } | ck::CipRankError::AtomMapOutOfRange { .. }
    | ck::CipRankError::InvariantCount { .. } | ck::CipRankError::InvariantOutOfRange { .. }
    | ck::CipRankError::TooManyNeighbors { .. } | ck::CipRankError::UnsupportedBondOrder { .. } => None,
}), source = Option<crate::error_framework::ErrorObject<MappingValidationError>>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::CipRankError::NegativeImplicitHydrogen { atom, .. }
        | ck::CipRankError::AtomMapOutOfRange { atom, .. }
        | ck::CipRankError::InvariantOutOfRange { atom, .. }
        | ck::CipRankError::TooManyNeighbors { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn value(&self) -> Option<i64> { match self.inner {
        ck::CipRankError::NegativeImplicitHydrogen { value, .. } => Some(i64::from(value)),
        ck::CipRankError::InvariantOutOfRange { value, .. } => Some(value), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::CipRankError::ValenceRowCount { actual, .. }
        | ck::CipRankError::InvariantCount { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::CipRankError::ValenceRowCount { atom_count, .. }
        | ck::CipRankError::InvariantCount { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> {
        match &self.inner {
            ck::CipRankError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_state(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<QueryStateError>>> {
        match &self.inner {
            ck::CipRankError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});

error_projection!(MoleculeHashError, MoleculeHashErrorKind, ck::MoleculeHashError, "hashing", {
    EmptyMolecule, MissingPreparedValence, CipRanks(_), InvalidTopology(_),
    RankCount { actual: usize, atom_count: usize },
}, cause = |inner: &ck::MoleculeHashError| match inner {
    ck::MoleculeHashError::CipRanks(value) => Some(CipRankError::from_rust(value.clone())),
    ck::MoleculeHashError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::MoleculeHashError::EmptyMolecule | ck::MoleculeHashError::MissingPreparedValence
    | ck::MoleculeHashError::RankCount { .. } => None,
}, source = Option<Source<crate::error_framework::ErrorObject<CipRankError>, crate::error_framework::ErrorObject<TopologyValidationError>>>);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;

    #[test]
    fn real_hash_queries_and_errors_preserve_native_results_and_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CCO").unwrap();
            let before = native.clone();
            let molecule = Py::new(
                py,
                Molecule {
                    inner: native.clone(),
                },
            )?;
            assert_eq!(
                molecule
                    .bind(py)
                    .call_method0("molecular_hash")?
                    .extract::<u64>()?,
                native.molecular_hash().unwrap()
            );
            let ranks = vec![2, 0, 1];
            assert_eq!(
                molecule
                    .bind(py)
                    .call_method1("molecular_hash_with_ranks", (ranks.clone(),))?
                    .extract::<u64>()?,
                native.molecular_hash_with_ranks(&ranks).unwrap()
            );
            assert_eq!(molecule.borrow(py).inner, before);
            assert!(std::ptr::eq(
                molecule.borrow(py).inner.topology(),
                before.topology()
            ));

            let mut builder = ck::MoleculeBuilder::new();
            builder.add_atom(ck::AtomSpec::new(ck::Element::C));
            for (value, expected) in [
                (ck::Molecule::new(), MoleculeHashErrorKind::EmptyMolecule),
                (
                    builder.build().unwrap(),
                    MoleculeHashErrorKind::MissingPreparedValence,
                ),
            ] {
                let display = value.molecular_hash().unwrap_err().to_string();
                let value = Py::new(py, Molecule { inner: value })?;
                let error = value.bind(py).call_method0("molecular_hash").unwrap_err();
                assert!(error.is_instance(py, <MoleculeHashError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                assert!(
                    error
                        .value(py)
                        .getattr("variant")?
                        .extract::<MoleculeHashErrorKind>()?
                        == expected
                );
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(error.cause(py).is_none());
            }
            let error = molecule
                .bind(py)
                .call_method1("molecular_hash_with_ranks", (vec![1_u32],))
                .unwrap_err();
            assert!(error.is_instance(py, <MoleculeHashError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).getattr("actual")?.extract::<usize>()?, 1);
            assert_eq!(
                error.value(py).getattr("atom_count")?.extract::<usize>()?,
                3
            );

            let mapping = ck::MappingValidationError::Length {
                entity: "atom",
                direction: "old-to-new",
                actual: 2,
                expected: 3,
            };
            let rank = ck::CipRankError::InvalidQueryState(ck::QueryStateError::Mapping(mapping));
            assert!(rank.source().unwrap().is::<ck::MappingValidationError>());
            let error = MoleculeHashError::from_rust(ck::MoleculeHashError::CipRanks(rank));
            let rank = error.cause(py).unwrap();
            assert!(rank.is_instance(py, <CipRankError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            let source = rank.cause(py).unwrap();
            assert!(source.is_instance(py, <MappingValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(rank.value(py).getattr("source")?.is(source.value(py)));
            assert!(
                rank.value(py)
                    .getattr("query_state")?
                    .is_instance(<QueryStateError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert!(error.value(py).getattr("source")?.is(rank.value(py)));
            assert!(source.cause(py).is_none());
            assert_eq!(source.value(py).getattr("actual")?.extract::<usize>()?, 2);
            Ok(())
        })
        .unwrap();
    }
}
