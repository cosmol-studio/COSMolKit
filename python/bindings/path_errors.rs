//! Source-shaped path errors; transparent tuple values are not invented causes.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::error_framework::{Source, error_projection};
use crate::matrix_errors::MatrixError;
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
};
use crate::model_errors::{StereoGroupError, TopologyValidationError};
use crate::query_model_errors::QueryGraphError;
use cosmolkit as ck;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use std::error::Error;

crate::enums::binding_enum!(PeriodicTableErrorKind, {
    UnknownIsotope,
});

crate::error_framework::error_projection! {record
pub struct PeriodicTableError {
    inner: ck::PeriodicTableError,
}
}

impl PeriodicTableError {
    fn fields(&self) -> PyResult<(u8, u16)> {
        match self.inner {
            ck::PeriodicTableError::UnknownIsotope { element, isotope } => {
                Ok((element.atomic_number(), isotope))
            }
            // The native public enum is non_exhaustive; never classify a
            // future Rust variant as the one known by this binding version.
            _ => Err(PyTypeError::new_err(
                "unrecognized native PeriodicTableError variant",
            )),
        }
    }
    pub(crate) fn from_rust(inner: ck::PeriodicTableError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}
crate::error_framework::error_projection! {methods impl PeriodicTableError {
    #[getter]
    fn domain(&self) -> &'static str {
        "periodic_table"
    }
    #[getter]
    fn variant(&self) -> PyResult<PeriodicTableErrorKind> {
        self.fields()?;
        Ok(PeriodicTableErrorKind::UnknownIsotope)
    }
    #[getter]
    fn element(&self) -> PyResult<u8> {
        self.fields().map(|fields| fields.0)
    }
    #[getter]
    fn isotope(&self) -> PyResult<u16> {
        self.fields().map(|fields| fields.1)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

pub(crate) type PathSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<MoleculePropertyError>,
        Source<
            crate::error_framework::ErrorObject<BondValueError>,
            Source<
                crate::error_framework::ErrorObject<AtomPropertyError>,
                crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
            >,
        >,
    >,
>;

error_projection!(PathError, PathErrorKind, ck::PathError, "paths", {
    StereoGroup(_), MoleculeProperty(_), BondProperty(_), AtomProperty(_), InvalidTopology(_),
    AtomOutOfRange { role: &'static str, atom, atom_count },
    EqualShortestPathEndpoints { atom }, InvalidLengthRange { lower: usize, upper: usize },
    LengthOverflow { length: usize }, AtomPathOutOfRange { position: usize, atom, atom_count },
    IgnoredAtomMaskLength { actual, expected }, ExtraInvariantLength { actual, expected },
    PeriodicTable(_), InvalidSubsetTopology(_), QueryGraph(_), Mapping(_),
    TemplateAttachmentRemap { carrier: usize => |value: &ck::AtomId| value.index(), source }, Matrix(_),
}, cause = |inner: &ck::PathError| Python::attach(|py| match inner {
    ck::PathError::StereoGroup(source) => Some(StereoGroupError::from_rust(*source)),
    ck::PathError::MoleculeProperty(source) => Some(MoleculePropertyError::from_rust(source.clone())),
    ck::PathError::BondProperty(source) => Some(BondValueError::from_rust(source.clone())),
    ck::PathError::AtomProperty(source) => Some(AtomPropertyError::from_rust(source.clone())),
    ck::PathError::InvalidTopology(source) => TopologyValidationError::from_rust(source.clone()).cause(py),
    ck::PathError::QueryGraph(source) => QueryGraphError::from_rust(source.clone()).cause(py),
    ck::PathError::Matrix(source) => MatrixError::from_rust(source.clone()).cause(py),
    ck::PathError::TemplateAttachmentRemap { source, .. } => Some(TemplateAttachmentOrderError::from_rust(source.clone())),
    ck::PathError::AtomOutOfRange { .. } | ck::PathError::EqualShortestPathEndpoints { .. }
    | ck::PathError::InvalidLengthRange { .. } | ck::PathError::LengthOverflow { .. }
    | ck::PathError::AtomPathOutOfRange { .. } | ck::PathError::IgnoredAtomMaskLength { .. }
    | ck::PathError::ExtraInvariantLength { .. } | ck::PathError::PeriodicTable(_)
    | ck::PathError::InvalidSubsetTopology(_) | ck::PathError::Mapping(_) => None,
}), source = Option<PathSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::PathError::AtomOutOfRange { atom, .. } | ck::PathError::EqualShortestPathEndpoints { atom }
        | ck::PathError::AtomPathOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::PathError::AtomOutOfRange { atom_count, .. }
        | ck::PathError::AtomPathOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::PathError::IgnoredAtomMaskLength { actual, .. } | ck::PathError::ExtraInvariantLength { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::PathError::IgnoredAtomMaskLength { expected, .. } | ck::PathError::ExtraInvariantLength { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> {
        match &self.inner {
            ck::PathError::InvalidTopology(value) | ck::PathError::InvalidSubsetTopology(value) =>
                TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn periodic_table(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PeriodicTableError>>> {
        match self.inner { ck::PathError::PeriodicTable(value) =>
            PeriodicTableError::from_rust(value).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_graph(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<QueryGraphError>>> {
        match &self.inner { ck::PathError::QueryGraph(value) =>
            QueryGraphError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn mapping(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MappingValidationError>>> {
        match &self.inner { ck::PathError::Mapping(value) =>
            MappingValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn matrix(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MatrixError>>> {
        match &self.inner { ck::PathError::Matrix(value) =>
            MatrixError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_attachment_rejection_and_transparent_path_payload_are_not_confused() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::TemplateAttachmentOrder::new(vec![
                ck::TemplateAttachment::new(ck::AtomId::new(2), "first"),
                ck::TemplateAttachment::new(ck::AtomId::new(2), "second"),
            ])
            .unwrap_err();
            let native = ck::PathError::TemplateAttachmentRemap {
                carrier: ck::AtomId::new(7),
                source,
            };
            let display = native.to_string();
            let error = PathError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(error.value(py).getattr("carrier")?.extract::<usize>()?, 7);
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(
                source
                    .extract::<crate::error_framework::ErrorObject<TemplateAttachmentOrderError>>()
                    .is_ok()
            );
            assert_eq!(source.getattr("target")?.extract::<usize>()?, 2);
            let native = ck::PathError::PeriodicTable(ck::PeriodicTableError::UnknownIsotope {
                element: ck::Element::from_atomic_number(6).unwrap(),
                isotope: 999,
            });
            assert!(native.source().is_none());
            let error = PathError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let payload = error.value(py).getattr("periodic_table")?;
            assert_eq!(payload.getattr("element")?.extract::<u8>()?, 6);
            assert_eq!(payload.getattr("isotope")?.extract::<u16>()?, 999);
            assert!(payload.getattr("__cause__")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
