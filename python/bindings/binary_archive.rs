use crate::error_framework::error_projection;
use crate::model_errors::StereoGroupError;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(PickleError, PickleErrorKind, ck::PickleError, "serialization", {
    StereoGroup(_), UnexpectedEof, UnsupportedVersion(_),
    UnsupportedArchiveVersion { major: u16, minor: u16 },
    UnsupportedSectionVersion { section, version },
    InvalidArchive(_), MissingRequiredSection(_), DuplicateSection(_), UnknownRequiredSection(_),
    DataLengthMismatch { expected: usize, actual: usize },
    InvalidEnumValue { value: u8, type_name: &'static str },
    InvalidMolecule(_), TooManyAtoms(_), TooManyBonds(_), StringTooLong(_),
}, cause = |inner: &ck::PickleError| match inner {
    ck::PickleError::StereoGroup(source) => Some(StereoGroupError::from_rust(*source)),
    ck::PickleError::UnexpectedEof | ck::PickleError::UnsupportedVersion(_)
    | ck::PickleError::UnsupportedArchiveVersion { .. }
    | ck::PickleError::UnsupportedSectionVersion { .. }
    | ck::PickleError::InvalidArchive(_) | ck::PickleError::MissingRequiredSection(_)
    | ck::PickleError::DuplicateSection(_) | ck::PickleError::UnknownRequiredSection(_)
    | ck::PickleError::DataLengthMismatch { .. } | ck::PickleError::InvalidEnumValue { .. }
    | ck::PickleError::InvalidMolecule(_) | ck::PickleError::TooManyAtoms(_)
    | ck::PickleError::TooManyBonds(_) | ck::PickleError::StringTooLong(_) => None,
}, source = Option<crate::error_framework::ErrorObject<StereoGroupError>>, methods = {
    #[getter]
    fn version(&self) -> Option<u16> {
        match self.inner {
            ck::PickleError::UnsupportedVersion(value) => Some(u16::from(value)),
            ck::PickleError::UnsupportedSectionVersion { version, .. } => Some(version),
            _ => None,
        }
    }
    #[getter]
    fn section(&self) -> Option<u16> {
        match self.inner {
            ck::PickleError::UnsupportedSectionVersion { section, .. }
            | ck::PickleError::MissingRequiredSection(section)
            | ck::PickleError::DuplicateSection(section)
            | ck::PickleError::UnknownRequiredSection(section) => Some(section),
            _ => None,
        }
    }
    #[getter]
    fn message(&self) -> Option<&str> {
        match &self.inner {
            ck::PickleError::InvalidArchive(message)
            | ck::PickleError::InvalidMolecule(message) => Some(message),
            _ => None,
        }
    }
    #[getter]
    fn length(&self) -> Option<usize> {
        match self.inner {
            ck::PickleError::TooManyAtoms(length) | ck::PickleError::TooManyBonds(length)
            | ck::PickleError::StringTooLong(length) => Some(length),
            _ => None,
        }
    }
});

#[pyfunction]
pub fn mol_from_binary(data: &[u8]) -> PyResult<Molecule> {
    Molecule::from_binary(data)
}

// Stable reconstruction entry used by existing cosmolkit pickle data.
#[pyfunction]
#[doc = r#"{"role":"private_export"}"#]
pub fn _molecule_from_binary(data: &[u8]) -> PyResult<Molecule> {
    Molecule::from_binary(data)
}

crate::callbacks::typed_callback!(
    MoleculeRebuild,
    [pyo3::type_hint_identifier!("builtins", "bytes")],
    <Py<Molecule> as IntoPyObject<'_>>::OUTPUT_TYPE
);

/// Python's index protocol without a fixed-width integer restriction.
pub struct PickleProtocol;
impl FromPyObject<'_, '_> for PickleProtocol {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr =
        pyo3::type_hint_identifier!("typing", "SupportsIndex");
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self> {
        value
            .py()
            .import("operator")?
            .getattr("index")?
            .call1((value,))?;
        Ok(Self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn archive_projection_retains_native_payloads_and_concrete_source() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_binary(&[0]).unwrap_err();
            let display = native.to_string();
            let error = PickleError::from_rust(native);
            assert_eq!(error.value(py).getattr("version")?.extract::<u16>()?, 0);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.cause(py).is_none());
            let native = ck::PickleError::StereoGroup(ck::StereoGroupError::DuplicateAtom);
            assert!(native.source().unwrap().is::<ck::StereoGroupError>());
            let error = PickleError::from_rust(native);
            let cause = error.cause(py).unwrap();
            assert!(
                cause.is_instance(
                    py,
                    <StereoGroupError as crate::error_framework::ExceptionPayload>::exception_type(
                        py
                    )
                    .expect("error backend")
                    .as_any()
                )
            );
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            let error = PickleError::from_rust(ck::PickleError::InvalidArchive(
                "exact\0payload".to_owned(),
            ));
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                "exact\0payload"
            );
            let error = PickleError::from_rust(ck::PickleError::TooManyAtoms(usize::MAX));
            assert_eq!(
                error.value(py).getattr("length")?.extract::<usize>()?,
                usize::MAX
            );
            Ok(())
        })
        .unwrap();
    }
}
