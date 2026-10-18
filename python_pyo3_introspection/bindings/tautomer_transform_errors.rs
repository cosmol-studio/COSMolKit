//! Complete internal catalog-transform error dependency, not a catalog fallback.
use crate::error_framework::error_projection;
use crate::smarts_errors::SmartsParseError;
use crate::string_errors::QueryCompileError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(TautomerTransformError, TautomerTransformErrorKind, ck::TautomerTransformError, "tautomer", {
    ChargeSymbolNotRecognised, CannotParseSmarts { smarts: String, source }, QueryCompile(_),
}, cause = |inner: &ck::TautomerTransformError| match inner {
    ck::TautomerTransformError::CannotParseSmarts { source, .. } => Some(SmartsParseError::from_rust(source.clone())),
    ck::TautomerTransformError::ChargeSymbolNotRecognised
    | ck::TautomerTransformError::QueryCompile(_) => None,
}, source = Option<Py<SmartsParseError>>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_compile_error(&self, py: Python<'_>) -> PyResult<Option<Py<QueryCompileError>>> {
        match &self.inner {
            ck::TautomerTransformError::QueryCompile(value) => QueryCompileError::from_rust(value.clone())
                .value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_catalog_transform_failures_keep_text_and_precise_transparent_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let error =
                ck::TautomerParams::from_transform_data(&[("bad", "[", "", "")]).unwrap_err();
            let ck::TautomerCatalogError::Transform(native) = error else {
                panic!("expected actual native transform error")
            };
            assert!(matches!(
                native,
                ck::TautomerTransformError::CannotParseSmarts { .. }
            ));
            assert!(native.source().unwrap().is::<ck::SmartsParseError>());
            let display = native.to_string();
            let error = TautomerTransformError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(error.value(py).getattr("smarts")?.extract::<String>()?, "[");
            let source = error.value(py).getattr("source")?;
            assert!(source.is_instance_of::<SmartsParseError>());
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(error.value(py).setattr("smarts", "C").is_err());
            assert!(error.value(py).getattr("query_compile_error")?.is_none());

            let error =
                ck::TautomerParams::from_transform_data(&[("charge", "C", "", "x")]).unwrap_err();
            let ck::TautomerCatalogError::Transform(native) = error else {
                panic!("expected actual native charge transform error")
            };
            assert!(matches!(
                native,
                ck::TautomerTransformError::ChargeSymbolNotRecognised
            ));
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = TautomerTransformError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).getattr("smarts")?.is_none());

            let native = ck::TautomerTransformError::QueryCompile(
                ck::QueryCompileError::InvalidGraph("bad\0graph".into()),
            );
            assert!(native.source().is_none());
            let error = TautomerTransformError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let value = error.value(py).getattr("query_compile_error")?;
            assert!(value.is_instance_of::<QueryCompileError>());
            assert_eq!(value.getattr("message")?.extract::<String>()?, "bad\0graph");
            assert!(
                error
                    .value(py)
                    .setattr("query_compile_error", py.None())
                    .is_err()
            );
            Ok(())
        })
        .unwrap();
    }
}
