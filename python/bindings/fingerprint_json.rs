//! Immutable native fingerprint configuration JSON, including exact parse causes.
use crate::io_source_errors::{ErrorRef, IoErrorKind, StandardSource, standard_source};
use crate::retained_error::RetainedError;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::enums::binding_enum!(SerdeJsonCategory, {
    Io = serde_json::error::Category::Io as isize,
    Syntax = serde_json::error::Category::Syntax as isize,
    Data = serde_json::error::Category::Data as isize,
    Eof = serde_json::error::Category::Eof as isize,
});
impl From<serde_json::error::Category> for SerdeJsonCategory {
    fn from(value: serde_json::error::Category) -> Self {
        match value {
            serde_json::error::Category::Io => Self::Io,
            serde_json::error::Category::Syntax => Self::Syntax,
            serde_json::error::Category::Data => Self::Data,
            serde_json::error::Category::Eof => Self::Eof,
        }
    }
}
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct SerdeJsonError {
    node: ErrorRef,
}
}
crate::error_framework::error_projection! {methods impl SerdeJsonError {
    #[getter]
    fn domain(&self) -> &'static str {
        "serde_json"
    }
    #[getter]
    fn line(&self) -> usize {
        self.native().line()
    }
    #[getter]
    fn column(&self) -> usize {
        self.native().column()
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn classify(&self) -> SerdeJsonCategory {
        self.native().classify().into()
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn io_error_kind(&self) -> Option<IoErrorKind> {
        self.native().io_error_kind().map(IoErrorKind::from_native)
    }
    #[getter]
    fn is_io(&self) -> bool {
        self.native().is_io()
    }
    #[getter]
    fn is_syntax(&self) -> bool {
        self.native().is_syntax()
    }
    #[getter]
    fn is_data(&self) -> bool {
        self.native().is_data()
    }
    #[getter]
    fn is_eof(&self) -> bool {
        self.native().is_eof()
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<StandardSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.native().to_string()
    }
}
}
impl SerdeJsonError {
    fn native(&self) -> &serde_json::Error {
        self.node
            .error()
            .downcast_ref()
            .expect("actual retained JSON parse error")
    }
    fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            let cause = node.source().map(standard_source);
            match Py::new(py, Self { node }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}

crate::enums::binding_enum!(FingerprintJsonErrorKind, {
    Parse,
    Invalid,
    UnsupportedComponent,
});
crate::error_framework::error_projection! {record
pub struct FingerprintJsonError {
    inner: RetainedError<ck::FingerprintJsonError>,
}
}
crate::error_framework::error_projection! {methods impl FingerprintJsonError {
    #[getter]
    fn domain(&self) -> &'static str {
        "fingerprint"
    }
    #[getter]
    fn variant(&self) -> FingerprintJsonErrorKind {
        match self.inner.as_ref() {
            ck::FingerprintJsonError::Parse(_) => FingerprintJsonErrorKind::Parse,
            ck::FingerprintJsonError::Invalid(_) => FingerprintJsonErrorKind::Invalid,
            ck::FingerprintJsonError::UnsupportedComponent {
                component: _,
                source_type: _,
            } => FingerprintJsonErrorKind::UnsupportedComponent,
        }
    }
    #[getter]
    fn invalid(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::FingerprintJsonError::Invalid(value) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn component(&self) -> Option<&'static str> {
        match self.inner.as_ref() {
            ck::FingerprintJsonError::UnsupportedComponent { component, .. } => Some(component),
            _ => None,
        }
    }
    #[getter]
    fn source_type(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::FingerprintJsonError::UnsupportedComponent { source_type, .. } => Some(source_type),
            _ => None,
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<crate::error_framework::ErrorObject<SerdeJsonError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl FingerprintJsonError {
    pub(crate) fn from_rust(inner: ck::FingerprintJsonError) -> PyErr {
        Self::from_ref(RetainedError::new(inner).node().clone())
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let inner = RetainedError::<ck::FingerprintJsonError>::from_ref(node);
        Python::attach(|py| {
            let cause = match inner.as_ref() {
                ck::FingerprintJsonError::Parse(_) => Some(SerdeJsonError::from_ref(
                    inner.node().source().expect("native JSON parse source"),
                )),
                ck::FingerprintJsonError::Invalid(_)
                | ck::FingerprintJsonError::UnsupportedComponent { .. } => None,
            };
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
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
    use crate::fingerprint_configuration::{
        AtomPairGeneratorParams, MorganGeneratorParams, TopologicalTorsionGeneratorParams,
    };
    #[test]
    fn three_native_json_configuration_families_preserve_updates_defaults_errors_and_owners() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            macro_rules! verify {
                ($name:ident) => {{
                    let native = ck::$name::default();
                    let original = native.to_json();
                    let params = Py::new(py, $name::from_rust(py, native.clone())?)?;
                    for text in [
                        "",
                        "   ",
                        "{}",
                        "{\"fpSize\":\"128\",\"countBounds\":[1,4,9]}",
                        "{\"countBounds\":[]}",
                        "{\"fpSize\":\"not an integer\"}",
                    ] {
                        let expected = native.with_json(text).unwrap();
                        let actual = params
                            .bind(py)
                            .call_method1("with_json", (text,))?
                            .extract::<Py<$name>>()?;
                        assert_eq!(actual.borrow(py).to_rust(py)?, expected);
                        assert_eq!(params.borrow(py).to_rust(py)?.to_json(), original);
                    }
                    let parse = "{\n  bad";
                    let native_error = native.with_json(parse).unwrap_err();
                    let ck::FingerprintJsonError::Parse(native_parse) = &native_error else {
                        panic!("actual parse variant")
                    };
                    let error = params
                        .bind(py)
                        .call_method1("with_json", (parse,))
                        .unwrap_err();
                    assert!(error.is_instance(py, <FingerprintJsonError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                    assert_eq!(error.value(py).str()?.to_str()?, native_error.to_string());
                    let source = error
                        .value(py)
                        .getattr("source")?
                        .extract::<crate::error_framework::ErrorObject<SerdeJsonError>>()?;
                    assert!(source.bind(py).is(&error.value(py).getattr("__cause__")?));
                    assert_eq!(
                        source.bind(py).getattr("line")?.extract::<usize>()?,
                        native_parse.line()
                    );
                    assert_eq!(
                        source.bind(py).getattr("column")?.extract::<usize>()?,
                        native_parse.column()
                    );
                    assert!(
                        source
                            .bind(py)
                            .getattr("classify")?
                            .extract::<SerdeJsonCategory>()?
                            == SerdeJsonCategory::from(native_parse.classify())
                    );
                    assert!(source.bind(py).getattr("io_error_kind")?.is_none());
                    assert!(source.bind(py).getattr("source")?.is_none());
                    assert!(source.bind(py).setattr("line", 0).is_err());
                    drop(error);
                    drop(params);
                    assert!(!source.bind(py).str()?.is_empty()?);
                    let error = native.with_json("{\"countBounds\":[\"bad\"]}").unwrap_err();
                    let error = FingerprintJsonError::from_rust(error);
                    assert!(
                        error
                            .value(py)
                            .getattr("variant")?
                            .extract::<FingerprintJsonErrorKind>()?
                            == FingerprintJsonErrorKind::Invalid
                    );
                    assert!(error.value(py).getattr("source")?.is_none());
                    assert!(
                        !error
                            .value(py)
                            .getattr("invalid")?
                            .extract::<String>()?
                            .is_empty()
                    );
                }};
            }
            verify!(MorganGeneratorParams);
            verify!(AtomPairGeneratorParams);
            verify!(TopologicalTorsionGeneratorParams);
            Ok(())
        })
        .unwrap();
    }
}
