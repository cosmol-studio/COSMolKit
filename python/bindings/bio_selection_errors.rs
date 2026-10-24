use crate::bio_errors::BioStructureError;
use crate::error_framework::{Source, error_projection};
use crate::protein_values::ProteinProjectionError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(BioRowModelError, BioRowModelErrorKind, ck::BioRowModelError, "bio", {
    MissingModelNumber,
}, cause = |inner: &ck::BioRowModelError| match inner {
    ck::BioRowModelError::MissingModelNumber => None::<PyErr>,
});

error_projection!(BioRowChainError, BioRowChainErrorKind, ck::BioRowChainError, "bio", {
    MissingCanonicalChainName,
}, cause = |inner: &ck::BioRowChainError| match inner {
    ck::BioRowChainError::MissingCanonicalChainName => None::<PyErr>,
});

error_projection!(BioRowTraverseError, BioRowTraverseErrorKind, ck::BioRowTraverseError, "bio", {
    Model(_), Chain(_),
}, cause = |inner: &ck::BioRowTraverseError| match inner {
    ck::BioRowTraverseError::Model(source) => Some(BioRowModelError::from_rust(*source)),
    ck::BioRowTraverseError::Chain(source) => Some(BioRowChainError::from_rust(*source)),
}, source = Source<crate::error_framework::ErrorObject<BioRowModelError>, crate::error_framework::ErrorObject<BioRowChainError>>);

error_projection!(BioSelectionCopyCause, BioSelectionCopyCauseKind, ck::BioSelectionCopyCause, "bio", {
    Structure(_), Traverse(_),
}, cause = |inner: &ck::BioSelectionCopyCause| match inner {
    ck::BioSelectionCopyCause::Structure(source) => Some(BioStructureError::from_rust(source.clone())),
    ck::BioSelectionCopyCause::Traverse(source) => Some(BioRowTraverseError::from_rust(*source)),
}, source = Source<crate::error_framework::ErrorObject<BioStructureError>, crate::error_framework::ErrorObject<BioRowTraverseError>>);

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct BioSelectionMatchError {
    inner: ck::BioSelectionMatchError,
}
}

crate::error_framework::error_projection! {methods impl BioSelectionMatchError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<crate::error_framework::ErrorObject<BioRowTraverseError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

impl BioSelectionMatchError {
    pub(crate) fn from_rust(inner: ck::BioSelectionMatchError) -> PyErr {
        Python::attach(|py| {
            let Some(source) = inner
                .source()
                .and_then(|error| error.downcast_ref::<ck::BioRowTraverseError>())
            else {
                return pyo3::exceptions::PyTypeError::new_err(
                    "BioSelectionMatchError source is not BioRowTraverseError",
                );
            };
            let cause = BioRowTraverseError::from_rust(*source);
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, Some(cause));
                    error
                }
                Err(error) => error,
            }
        })
    }
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct BioSelectionCopyError {
    inner: ck::BioSelectionCopyError,
}
}

crate::error_framework::error_projection! {methods impl BioSelectionCopyError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(
        slf: Bound<'_, pyo3::exceptions::PyBaseException>,
    ) -> PyResult<Source<crate::error_framework::ErrorObject<BioStructureError>, crate::error_framework::ErrorObject<BioRowTraverseError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

crate::enums::binding_enum!(BioOperationErrorKind, {
    Structure,
    Protein,
    Selection,
});

crate::error_framework::error_projection! {record
pub struct BioOperationError {
    kind: BioOperationErrorKind,
    message: String,
}
}

crate::error_framework::error_projection! {methods impl BioOperationError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn variant(&self) -> BioOperationErrorKind {
        self.kind
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(
        slf: Bound<'_, pyo3::exceptions::PyBaseException>,
    ) -> PyResult<
        Source<
            crate::error_framework::ErrorObject<BioStructureError>,
            Source<crate::error_framework::ErrorObject<ProteinProjectionError>, crate::error_framework::ErrorObject<BioSelectionCopyError>>,
        >,
    > {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> &str {
        &self.message
    }
}
}

impl BioOperationError {
    pub(crate) fn from_rust(inner: ck::BioOperationError) -> PyErr {
        Python::attach(|py| {
            let message = inner.to_string();
            // Consume each native cause once. SelectionCopyError is not Clone;
            // retaining a second root error would force duplicate ownership.
            let (kind, cause) = match inner {
                ck::BioOperationError::Structure(source) => (
                    BioOperationErrorKind::Structure,
                    BioStructureError::from_rust(source),
                ),
                ck::BioOperationError::Protein(source) => (
                    BioOperationErrorKind::Protein,
                    ProteinProjectionError::from_rust(source),
                ),
                ck::BioOperationError::Selection(source) => (
                    BioOperationErrorKind::Selection,
                    BioSelectionCopyError::from_rust(source),
                ),
            };
            match Py::new(py, Self { kind, message }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, Some(cause));
                    error
                }
                Err(error) => error,
            }
        })
    }
}

impl BioSelectionCopyError {
    pub(crate) fn from_rust(inner: ck::BioSelectionCopyError) -> PyErr {
        Python::attach(|py| {
            // Rust Error::source skips SelectionCopyCause itself. Preserve that
            // exact graph rather than inserting an extra enum exception.
            let cause = match inner.cause() {
                ck::BioSelectionCopyCause::Structure(source) => {
                    BioStructureError::from_rust(source.clone())
                }
                ck::BioSelectionCopyCause::Traverse(source) => {
                    BioRowTraverseError::from_rust(*source)
                }
            };
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, Some(cause));
                    error
                }
                Err(error) => error,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registered_bio_source_graph_preserves_all_concrete_branches() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let model = BioRowModelError::from_rust(ck::BioRowModelError::MissingModelNumber);
            assert!(model.cause(py).is_none());
            assert!(
                model
                    .value(py)
                    .getattr("variant")?
                    .extract::<BioRowModelErrorKind>()?
                    == BioRowModelErrorKind::MissingModelNumber
            );
            let chain =
                BioRowChainError::from_rust(ck::BioRowChainError::MissingCanonicalChainName);
            assert!(chain.cause(py).is_none());
            assert!(
                chain
                    .value(py)
                    .getattr("variant")?
                    .extract::<BioRowChainErrorKind>()?
                    == BioRowChainErrorKind::MissingCanonicalChainName
            );
            for (native, expected) in [
                (
                    ck::BioRowTraverseError::Model(ck::BioRowModelError::MissingModelNumber),
                    BioRowTraverseErrorKind::Model,
                ),
                (
                    ck::BioRowTraverseError::Chain(ck::BioRowChainError::MissingCanonicalChainName),
                    BioRowTraverseErrorKind::Chain,
                ),
            ] {
                let error =
                    BioSelectionCopyCause::from_rust(ck::BioSelectionCopyCause::Traverse(native));
                assert!(
                    error
                        .value(py)
                        .getattr("variant")?
                        .extract::<BioSelectionCopyCauseKind>()?
                        == BioSelectionCopyCauseKind::Traverse
                );
                let cause = error.cause(py).unwrap();
                assert!(cause.is_instance(py, <BioRowTraverseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                assert!(
                    cause
                        .value(py)
                        .getattr("variant")?
                        .extract::<BioRowTraverseErrorKind>()?
                        == expected
                );
                let leaf = cause.cause(py).unwrap();
                match native {
                    ck::BioRowTraverseError::Model(_) => {
                        assert!(leaf.is_instance(py, <BioRowModelError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()))
                    }
                    ck::BioRowTraverseError::Chain(_) => {
                        assert!(leaf.is_instance(py, <BioRowChainError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()))
                    }
                }
                assert!(leaf.cause(py).is_none());
            }
            let error = BioSelectionCopyCause::from_rust(ck::BioSelectionCopyCause::Structure(
                ck::BioStructureError::RowIndexTooLarge { value: usize::MAX },
            ));
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<BioSelectionCopyCauseKind>()?
                    == BioSelectionCopyCauseKind::Structure
            );
            let cause = error.cause(py).unwrap();
            assert!(cause.is_instance(py, <BioStructureError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(
                cause.value(py).getattr("value")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(cause.cause(py).is_none());
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn source_only_records_keep_native_cause_and_identity() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let text = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n";
            let selection = ck::BioSelection::from_cid("/").unwrap();
            for missing_model in [true, false] {
                let mut parts = ck::BioStructure::from_pdb(text).unwrap().into_parts();
                if missing_model {
                    parts.models[0] = ck::BioModelRow::new(parts.models[0].chain_span(), None);
                } else {
                    let chain = &parts.chains[0];
                    parts.chains[0] = ck::BioChainRow::new(chain.model_id(), chain.entity_id(),
                        chain.residue_span(), chain.kind(), ck::ChainSourceIds::new(None, None));
                }
                let structure = ck::BioStructure::from_parts(parts).unwrap();
                let native = structure.selected_atom_ids(&selection).unwrap_err();
                let display = native.to_string();
                let error = BioSelectionMatchError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let cause = error.cause(py).unwrap();
                assert!(cause.is_instance(py, <BioRowTraverseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                assert!(error.value(py).getattr("source")?.is(cause.value(py)));
                assert!(cause.value(py).getattr("source")?.is(cause.cause(py).unwrap().value(py)));
                assert!(error.value(py).getattr("variant").is_err());
                assert!(error.value(py).setattr("source", py.None()).is_err());
                let native = structure.with_selection(&selection).unwrap_err();
                let ck::BioOperationError::Selection(copy) = native else { panic!("{native}") };
                let display = copy.to_string();
                let error = BioSelectionCopyError::from_rust(copy);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let cause = error.cause(py).unwrap();
                assert!(cause.is_instance(py, <BioRowTraverseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                assert!(error.value(py).getattr("source")?.is(cause.value(py)));
                assert!(error.value(py).getattr("variant").is_err());
            }
            Ok(())
        }).unwrap();
    }
}
