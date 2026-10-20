//! Full mmCIF writer projection; serialization and defaults remain in Rust IO.
use crate::bio_errors::BioStructureError;
use crate::chemical_values::enum_value;
use crate::error_framework::Source;
use crate::io_source_errors::{ErrorRef, IoError};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::Arc;

enum_value!(CifReadErrorKind, { Syntax, MissingValue, DuplicateName, InvalidLoop, InvalidValue, OutOfRange });

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit", extends=PyValueError)]
pub struct CifReadError {
    inner: ck::CifReadError,
}
#[pymethods]
impl CifReadError {
    #[getter]
    fn domain(&self) -> &'static str {
        "cif"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> CifReadErrorKind {
        self.inner.kind().into()
    }
    #[getter]
    fn source(&self) -> &str {
        self.inner.source()
    }
    #[getter]
    fn line(&self) -> usize {
        self.inner.line()
    }
    #[getter]
    fn column(&self) -> usize {
        self.inner.column()
    }
    #[getter]
    fn message(&self) -> &str {
        self.inner.message()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl CifReadError {
    pub(crate) fn from_rust(inner: ck::CifReadError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

crate::enums::binding_enum!(BioMmcifWriteErrorKind, {
    Structure,
    Cif,
    InvalidText,
    Io,
    FileWrite,
});
#[pyclass(module="cosmolkit", extends=PyValueError)]
pub struct BioMmcifWriteError {
    inner: Arc<ck::BioMmcifWriteError>,
}
#[pymethods]
impl BioMmcifWriteError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn variant(&self) -> BioMmcifWriteErrorKind {
        match self.inner.as_ref() {
            ck::BioMmcifWriteError::Structure(_) => BioMmcifWriteErrorKind::Structure,
            ck::BioMmcifWriteError::Cif(_) => BioMmcifWriteErrorKind::Cif,
            ck::BioMmcifWriteError::InvalidText { field: _ } => BioMmcifWriteErrorKind::InvalidText,
            ck::BioMmcifWriteError::Io(_) => BioMmcifWriteErrorKind::Io,
            ck::BioMmcifWriteError::FileWrite { path: _, source: _ } => {
                BioMmcifWriteErrorKind::FileWrite
            }
        }
    }
    #[getter]
    fn field(&self) -> Option<&'static str> {
        match self.inner.as_ref() {
            ck::BioMmcifWriteError::InvalidText { field } => Some(field),
            _ => None,
        }
    }
    #[getter]
    fn path(&self) -> Option<PathBuf> {
        match self.inner.as_ref() {
            ck::BioMmcifWriteError::FileWrite { path, .. } => Some(path.clone()),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(
        slf: Bound<'_, Self>,
    ) -> PyResult<Option<Source<Py<BioStructureError>, Source<Py<CifReadError>, Py<IoError>>>>>
    {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio_structure::BioStructure;
    use crate::bio_write_configuration::BioMmcifWriteParams;
    use pyo3::types::PyDict;

    #[test]
    fn actual_mmcif_text_files_options_and_typed_sources_preserve_the_structure() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::BioStructure::from_pdb("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n").unwrap();
            let before = native.to_pdb().unwrap();
            let baseline = native.to_mmcif().unwrap();
            let structure = Py::new(py, BioStructure { inner: native })?;
            assert_eq!(structure.bind(py).call_method0("to_mmcif")?.extract::<String>()?, baseline);
            let directory = std::env::temp_dir().join(format!("ck-pyo3-mmcif-{}", std::process::id()));
            std::fs::create_dir_all(&directory).unwrap();
            let path = directory.join("out.cif");
            let as_path = py.import("pathlib")?.getattr("Path")?.call1((path.clone(),))?;
            for (field, value) in [
                ("atoms", false.into_pyobject(py)?.to_owned().into_any()),
                ("compact", true.into_pyobject(py)?.to_owned().into_any()),
                ("auth_all", true.into_pyobject(py)?.to_owned().into_any()),
                ("misuse_hash", true.into_pyobject(py)?.to_owned().into_any()),
                ("align_pairs", 12_u16.into_pyobject(py)?.into_any()),
                ("align_loops", 6_u16.into_pyobject(py)?.into_any()),
            ] {
                let params = Py::new(py, BioMmcifWriteParams { inner: crate::parameters::NativeField::Owned(ck::BioMmcifWriteParams::default() )})?;
                params.bind(py).setattr(field, value.clone())?;
                let expected = structure.borrow(py).inner.to_mmcif_with_params(&params.borrow(py).to_rust(py)?).unwrap();
                let kwargs = PyDict::new(py);
                kwargs.set_item(field, value)?;
                for actual in [
                    structure.bind(py).call_method1("to_mmcif", (params.clone_ref(py),))?,
                    structure.bind(py).call_method1("to_mmcif_with_params", (params.clone_ref(py),))?,
                    structure.bind(py).call_method("to_mmcif", (), Some(&kwargs))?,
                ] { assert_eq!(actual.extract::<String>()?, expected); }
                structure.bind(py).call_method1("write_mmcif", (&as_path, params.clone_ref(py)))?;
                assert_eq!(std::fs::read_to_string(&path).unwrap(), expected);
                structure.bind(py).call_method1("write_mmcif_with_params", (path.clone(), params.clone_ref(py)))?;
                assert_eq!(std::fs::read_to_string(&path).unwrap(), expected);
                structure.bind(py).call_method("write_mmcif", (&as_path,), Some(&kwargs))?;
                assert_eq!(std::fs::read_to_string(&path).unwrap(), expected);
                assert!(structure.bind(py).call_method("to_mmcif", (params.clone_ref(py),), Some(&kwargs)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
                assert!(structure.bind(py).call_method("write_mmcif", (&as_path, params), Some(&kwargs)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            structure.bind(py).call_method1("write_mmcif", (&as_path,))?;
            assert_eq!(std::fs::read_to_string(&path).unwrap(), baseline);
            let missing = directory.join("missing").join("out.cif");
            let error = structure.bind(py).call_method1("write_mmcif", (missing.clone(),)).unwrap_err();
            assert!(error.is_instance_of::<BioMmcifWriteError>(py));
            assert!(error.value(py).getattr("variant")?.extract::<BioMmcifWriteErrorKind>()? == BioMmcifWriteErrorKind::FileWrite);
            assert_eq!(error.value(py).getattr("path")?.extract::<PathBuf>()?, missing);
            let retained = error.value(py).getattr("source")?.unbind();
            assert!(retained.bind(py).is_instance_of::<IoError>());
            assert!(retained.bind(py).is(&error.value(py).getattr("__cause__")?));
            assert!(error.value(py).setattr("path", "different").is_err());
            drop(error);
            std::fs::remove_file(&path).unwrap();
            std::fs::remove_dir(&directory).unwrap();
            assert!(!retained.bind(py).str()?.is_empty()?);
            assert_eq!(structure.borrow(py).inner.to_pdb().unwrap(), before);
            assert_eq!(structure.borrow(py).inner.to_mmcif().unwrap(), baseline);
            let error = BioMmcifWriteError::from_rust(ck::BioMmcifWriteError::InvalidText { field: "name\0field" });
            assert_eq!(error.value(py).getattr("field")?.extract::<String>()?, "name\0field");
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("field", "changed").is_err());
            let error = BioMmcifWriteError::from_rust(ck::BioMmcifWriteError::Structure(ck::BioStructureError::CoordinateCountMismatch { atom_count: 2, coordinate_count: 1 }));
            assert!(error.value(py).getattr("source")?.is_instance_of::<BioStructureError>());
            assert!(error.value(py).getattr("source")?.is(&error.value(py).getattr("__cause__")?));
            Ok(())
        }).unwrap();
    }
}
impl BioMmcifWriteError {
    pub(crate) fn from_rust(inner: ck::BioMmcifWriteError) -> PyErr {
        let inner = Arc::new(inner);
        Python::attach(|py| {
            let cause = match inner.as_ref() {
                ck::BioMmcifWriteError::Structure(value) => {
                    Some(BioStructureError::from_rust(value.clone()))
                }
                ck::BioMmcifWriteError::Cif(value) => Some(CifReadError::from_rust(value.clone())),
                ck::BioMmcifWriteError::Io(_) | ck::BioMmcifWriteError::FileWrite { .. } => {
                    let root: Arc<dyn std::error::Error + Send + Sync> = inner.clone();
                    Some(IoError::from_ref(
                        ErrorRef::new(root)
                            .source()
                            .expect("native mmCIF IO source"),
                    ))
                }
                ck::BioMmcifWriteError::InvalidText { .. } => None,
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
