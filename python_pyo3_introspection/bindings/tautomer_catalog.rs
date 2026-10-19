//! Catalog factories directly call the public Rust facade with exact failures.
use crate::error_framework::Source;
use crate::io_source_errors::{ErrorRef, IoError};
use crate::smarts_errors::SmartsParseError;
use crate::tautomer_configuration::TautomerParams;
use crate::tautomer_transform_errors::TautomerTransformError;
use crate::utf8_errors::FromUtf8Error;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::Arc;

#[pyclass(module = "cosmolkit", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum TautomerCatalogErrorKind {
    BadInputFile,
    BadStreamContents,
    InvalidUtf8,
    Transform,
    TransformIndexOutOfRange,
    DeserializationUnderConstruction,
}

#[pyclass(module = "cosmolkit", extends=PyValueError)]
pub struct TautomerCatalogError {
    inner: Arc<ck::TautomerCatalogError>,
}
#[pymethods]
impl TautomerCatalogError {
    #[getter]
    fn domain(&self) -> &'static str {
        "tautomer"
    }
    #[getter]
    fn variant(&self) -> TautomerCatalogErrorKind {
        match self.inner.as_ref() {
            ck::TautomerCatalogError::BadInputFile { path: _, source: _ } => {
                TautomerCatalogErrorKind::BadInputFile
            }
            ck::TautomerCatalogError::BadStreamContents(_) => {
                TautomerCatalogErrorKind::BadStreamContents
            }
            ck::TautomerCatalogError::InvalidUtf8(_) => TautomerCatalogErrorKind::InvalidUtf8,
            ck::TautomerCatalogError::Transform(_) => TautomerCatalogErrorKind::Transform,
            ck::TautomerCatalogError::TransformIndexOutOfRange { index: _, len: _ } => {
                TautomerCatalogErrorKind::TransformIndexOutOfRange
            }
            ck::TautomerCatalogError::DeserializationUnderConstruction => {
                TautomerCatalogErrorKind::DeserializationUnderConstruction
            }
        }
    }
    #[getter]
    fn path(&self) -> Option<PathBuf> {
        match self.inner.as_ref() {
            ck::TautomerCatalogError::BadInputFile { path, .. } => Some(path.clone()),
            _ => None,
        }
    }
    #[getter]
    fn index(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::TautomerCatalogError::TransformIndexOutOfRange { index, .. } => Some(*index),
            _ => None,
        }
    }
    #[getter]
    fn len(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::TautomerCatalogError::TransformIndexOutOfRange { len, .. } => Some(*len),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn transform(&self, py: Python<'_>) -> PyResult<Option<Py<TautomerTransformError>>> {
        match self.inner.as_ref() {
            ck::TautomerCatalogError::Transform(value) => {
                TautomerTransformError::from_rust(value.clone())
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(
        slf: Bound<'_, Self>,
    ) -> PyResult<Option<Source<Py<IoError>, Source<Py<FromUtf8Error>, Py<SmartsParseError>>>>>
    {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl TautomerCatalogError {
    pub(crate) fn from_rust(inner: ck::TautomerCatalogError) -> PyErr {
        let inner = Arc::new(inner);
        Python::attach(|py| {
            let cause = match inner.as_ref() {
                ck::TautomerCatalogError::BadInputFile { .. }
                | ck::TautomerCatalogError::BadStreamContents(_) => {
                    let root: Arc<dyn std::error::Error + Send + Sync> = inner.clone();
                    Some(IoError::from_ref(
                        ErrorRef::new(root)
                            .source()
                            .expect("native catalog IO source"),
                    ))
                }
                ck::TautomerCatalogError::InvalidUtf8(value) => {
                    Some(FromUtf8Error::from_rust(value.clone()))
                }
                // Transparent native errors forward the child's actual source,
                // while the wrapped transform remains a separate typed payload.
                ck::TautomerCatalogError::Transform(value) => {
                    TautomerTransformError::from_rust(value.clone()).cause(py)
                }
                ck::TautomerCatalogError::TransformIndexOutOfRange { .. }
                | ck::TautomerCatalogError::DeserializationUnderConstruction => None,
            };
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}

pub(crate) fn v1(py: Python<'_>) -> PyResult<TautomerParams> {
    let native = ck::TautomerParams::v1().map_err(TautomerCatalogError::from_rust)?;
    TautomerParams::from_native(py, native)
}
pub(crate) fn from_transform_data(
    py: Python<'_>,
    data: Vec<(String, String, String, String)>,
) -> PyResult<TautomerParams> {
    let rows: Vec<_> = data
        .iter()
        .map(|(name, smarts, bonds, charges)| {
            (
                name.as_str(),
                smarts.as_str(),
                bonds.as_str(),
                charges.as_str(),
            )
        })
        .collect();
    let native =
        ck::TautomerParams::from_transform_data(&rows).map_err(TautomerCatalogError::from_rust)?;
    TautomerParams::from_native(py, native)
}
pub(crate) fn from_transform_file(py: Python<'_>, path: PathBuf) -> PyResult<TautomerParams> {
    let native =
        ck::TautomerParams::from_transform_file(path).map_err(TautomerCatalogError::from_rust)?;
    TautomerParams::from_native(py, native)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::ffi::c_str;
    use pyo3::types::PyModule;

    #[test]
    fn real_catalog_factories_preserve_defaults_paths_typed_causes_and_owned_lifetimes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "tautomer_catalog_test")?;
            module.add_class::<TautomerParams>()?;
            module.add_class::<TautomerCatalogError>()?;
            module.add_class::<TautomerCatalogErrorKind>()?;
            module.add_class::<IoError>()?;
            module.add_class::<crate::io_source_errors::IoErrorKind>()?;
            module.add_class::<FromUtf8Error>()?;
            py.run(
                c_str!(
                    r#"
from pathlib import Path
from tempfile import TemporaryDirectory
import os
old = TautomerParams.v1()
assert old.transform_count() > 0
custom = TautomerParams.from_transform_data([('identity', 'C', '', '')])
assert custom.transform_count() == 1
assert custom.max_tautomers == TautomerParams().max_tautomers
assert TautomerParams.from_transform_data([]).transform_count() == 0
with TemporaryDirectory() as folder:
    path = Path(folder) / 'transforms.tsv'
    path.write_text('identity\tC\n', encoding='utf-8')
    assert TautomerParams.from_transform_file(path).transform_count() == 1
    assert TautomerParams.from_transform_file(os.fspath(path)).transform_count() == 1
    missing = Path(folder) / 'missing.tsv'
    try:
        TautomerParams.from_transform_file(missing)
    except TautomerCatalogError as error:
        assert error.variant == TautomerCatalogErrorKind.BadInputFile
        assert os.fspath(error.path) == os.fspath(missing)
        saved = error.source
        assert saved is error.__cause__ and isinstance(saved, IoError)
        assert saved.kind == IoErrorKind.NotFound
        assert isinstance(saved.raw_os_error, int)
        assert saved.source is None
        try: saved.raw_os_error = 0
        except AttributeError: pass
        else: raise AssertionError('mutable IO payload')
    else: raise AssertionError('missing catalog accepted')
    path.write_bytes(b'\xff\n')
    try: TautomerParams.from_transform_file(path)
    except TautomerCatalogError as error:
        assert error.variant == TautomerCatalogErrorKind.InvalidUtf8
        assert isinstance(error.source, FromUtf8Error)
        assert error.source is error.__cause__
        assert error.source.as_bytes() == b'\xff'
    else: raise AssertionError('invalid UTF-8 accepted')
assert saved.kind == IoErrorKind.NotFound and str(saved)
try: TautomerParams.from_transform_data([('bad', '[', '', '')])
except TautomerCatalogError as error:
    assert error.variant == TautomerCatalogErrorKind.Transform
    assert error.transform is not None
    assert error.source is error.__cause__
    assert error.source is not error.transform
    assert error.source is not None
else: raise AssertionError('invalid transform accepted')
"#
                ),
                Some(&module.dict()),
                None,
            )?;
            let old = module.getattr("old")?.extract::<PyRef<TautomerParams>>()?;
            assert_eq!(
                old.to_rust(py)?.transform_count(),
                ck::TautomerParams::v1().unwrap().transform_count()
            );
            Ok(())
        })
        .unwrap();
    }
}
