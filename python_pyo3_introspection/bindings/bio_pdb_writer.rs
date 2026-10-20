//! BIO PDB serialization forwards the public owner and preserves native failures.
use crate::bio_errors::BioStructureError;
use crate::error_framework::Source;
use crate::io_source_errors::{ErrorRef, IoError};
use crate::params::BioPdbWriteParams;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::Arc;

crate::enums::binding_enum!(BioPdbWriteErrorKind, {
    InvalidStructure,
    ChainNameTooLong,
    NegativeSerial,
    SerialIncrementOverflow,
    SerialOffsetOverflow,
    NegativeBase36Value,
    UnrepresentableField,
    Io,
});

#[pyclass(module = "cosmolkit", extends=PyValueError)]
pub struct BioPdbWriteError {
    inner: Arc<ck::BioPdbWriteError>,
}
#[pymethods]
impl BioPdbWriteError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn variant(&self) -> BioPdbWriteErrorKind {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::InvalidStructure(_) => BioPdbWriteErrorKind::InvalidStructure,
            ck::BioPdbWriteError::ChainNameTooLong {
                chain: _,
                length: _,
            } => BioPdbWriteErrorKind::ChainNameTooLong,
            ck::BioPdbWriteError::NegativeSerial { serial: _ } => {
                BioPdbWriteErrorKind::NegativeSerial
            }
            ck::BioPdbWriteError::SerialIncrementOverflow { from: _ } => {
                BioPdbWriteErrorKind::SerialIncrementOverflow
            }
            ck::BioPdbWriteError::SerialOffsetOverflow { serial: _ } => {
                BioPdbWriteErrorKind::SerialOffsetOverflow
            }
            ck::BioPdbWriteError::NegativeBase36Value { value: _ } => {
                BioPdbWriteErrorKind::NegativeBase36Value
            }
            ck::BioPdbWriteError::UnrepresentableField {
                field: _,
                detail: _,
            } => BioPdbWriteErrorKind::UnrepresentableField,
            ck::BioPdbWriteError::Io { path: _, source: _ } => BioPdbWriteErrorKind::Io,
        }
    }
    #[getter]
    fn chain(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::ChainNameTooLong { chain, .. } => Some(chain),
            _ => None,
        }
    }
    #[getter]
    fn length(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::ChainNameTooLong { length, .. } => Some(*length),
            _ => None,
        }
    }
    #[getter]
    fn serial(&self) -> Option<i32> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::NegativeSerial { serial }
            | ck::BioPdbWriteError::SerialOffsetOverflow { serial } => Some(*serial),
            _ => None,
        }
    }
    #[getter]
    fn from_(&self) -> Option<i32> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::SerialIncrementOverflow { from } => Some(*from),
            _ => None,
        }
    }
    #[getter]
    fn value(&self) -> Option<i32> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::NegativeBase36Value { value } => Some(*value),
            _ => None,
        }
    }
    #[getter]
    fn field(&self) -> Option<&'static str> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::UnrepresentableField { field, .. } => Some(field),
            _ => None,
        }
    }
    #[getter]
    fn detail(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::UnrepresentableField { detail, .. } => Some(detail),
            _ => None,
        }
    }
    #[getter]
    fn path(&self) -> Option<PathBuf> {
        match self.inner.as_ref() {
            ck::BioPdbWriteError::Io { path, .. } => Some(path.clone()),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(
        slf: Bound<'_, Self>,
    ) -> PyResult<Option<Source<Py<BioStructureError>, Py<IoError>>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl BioPdbWriteError {
    pub(crate) fn from_rust(inner: ck::BioPdbWriteError) -> PyErr {
        let inner = Arc::new(inner);
        Python::attach(|py| {
            let cause = match inner.as_ref() {
                ck::BioPdbWriteError::InvalidStructure(value) => {
                    Some(BioStructureError::from_rust(value.clone()))
                }
                ck::BioPdbWriteError::Io { .. } => {
                    let root: Arc<dyn std::error::Error + Send + Sync> = inner.clone();
                    Some(IoError::from_ref(
                        ErrorRef::new(root)
                            .source()
                            .expect("native BIO PDB IO source"),
                    ))
                }
                ck::BioPdbWriteError::ChainNameTooLong { .. }
                | ck::BioPdbWriteError::NegativeSerial { .. }
                | ck::BioPdbWriteError::SerialIncrementOverflow { .. }
                | ck::BioPdbWriteError::SerialOffsetOverflow { .. }
                | ck::BioPdbWriteError::NegativeBase36Value { .. }
                | ck::BioPdbWriteError::UnrepresentableField { .. } => None,
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

pub(crate) fn keywords(
    ter_records: Option<bool>,
    numbered_ter: Option<bool>,
    ter_ignores_type: Option<bool>,
    preserve_serial: Option<bool>,
    end_record: Option<bool>,
) -> BioPdbWriteParams {
    let defaults = ck::BioPdbWriteParams::default();
    BioPdbWriteParams {
        inner: crate::parameters::NativeField::Owned(ck::BioPdbWriteParams {
            ter_records: ter_records.unwrap_or(defaults.ter_records),
            numbered_ter: numbered_ter.unwrap_or(defaults.numbered_ter),
            ter_ignores_type: ter_ignores_type.unwrap_or(defaults.ter_ignores_type),
            preserve_serial: preserve_serial.unwrap_or(defaults.preserve_serial),
            end_record: end_record.unwrap_or(defaults.end_record),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio_structure::BioStructure;
    use pyo3::types::PyDict;

    #[test]
    fn complete_pdb_writer_family_forwards_text_files_options_exact_errors_and_keeps_source() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::BioStructure::from_pdb("ATOM     17  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n").unwrap();
            let before = native.to_mmcif().unwrap();
            let baseline = native.to_pdb().unwrap();
            let structure = Py::new(py, BioStructure { inner: native })?;
            assert_eq!(structure.bind(py).call_method0("to_pdb")?.extract::<String>()?, baseline);
            let folder = py.import("tempfile")?.getattr("TemporaryDirectory")?.call0()?;
            let directory = PathBuf::from(folder.getattr("name")?.extract::<String>()?);
            let path = directory.join("output.pdb");
            let as_path = py.import("pathlib")?.getattr("Path")?.call1((path.clone(),))?;
            for bits in 0..32 {
                let native_params = ck::BioPdbWriteParams {
                    ter_records: bits & 1 != 0, numbered_ter: bits & 2 != 0,
                    ter_ignores_type: bits & 4 != 0, preserve_serial: bits & 8 != 0,
                    end_record: bits & 16 != 0,
                };
                let expected = structure.borrow(py).inner.to_pdb_with_params(&native_params).unwrap();
                let params = Py::new(py, BioPdbWriteParams { inner: crate::parameters::NativeField::Owned(native_params )})?;
                let kwargs = PyDict::new(py);
                kwargs.set_item("ter_records", native_params.ter_records)?;
                kwargs.set_item("numbered_ter", native_params.numbered_ter)?;
                kwargs.set_item("ter_ignores_type", native_params.ter_ignores_type)?;
                kwargs.set_item("preserve_serial", native_params.preserve_serial)?;
                kwargs.set_item("end_record", native_params.end_record)?;
                for actual in [
                    structure.bind(py).call_method1("to_pdb", (params.clone_ref(py),))?,
                    structure.bind(py).call_method1("to_pdb_with_params", (params.clone_ref(py),))?,
                    structure.bind(py).call_method("to_pdb", (), Some(&kwargs))?,
                ] { assert_eq!(actual.extract::<String>()?, expected); }
                structure.bind(py).call_method1("write_pdb_with_params", (&as_path, params.clone_ref(py)))?;
                assert_eq!(std::fs::read_to_string(&path).unwrap(), expected);
                structure.bind(py).call_method1("write_pdb", (path.clone(), params.clone_ref(py)))?;
                assert_eq!(std::fs::read_to_string(&path).unwrap(), expected);
                structure.bind(py).call_method("write_pdb", (&as_path,), Some(&kwargs))?;
                assert_eq!(std::fs::read_to_string(&path).unwrap(), expected);
                assert!(structure.bind(py).call_method("to_pdb", (params.clone_ref(py),), Some(&kwargs)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
                assert!(structure.bind(py).call_method("write_pdb", (&as_path, params), Some(&kwargs)).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            structure.bind(py).call_method1("write_pdb", (&as_path,))?;
            assert_eq!(std::fs::read_to_string(&path).unwrap(), baseline);
            let missing = directory.join("missing").join("out.pdb");
            let error = structure.bind(py).call_method1("write_pdb", (missing.clone(),)).unwrap_err();
            assert!(error.is_instance_of::<BioPdbWriteError>(py));
            assert_eq!(error.value(py).getattr("path")?.extract::<PathBuf>()?, missing);
            let retained = error.value(py).getattr("source")?.unbind();
            assert!(retained.bind(py).is_instance_of::<IoError>());
            assert!(retained.bind(py).is(&error.value(py).getattr("__cause__")?));
            assert!(retained.bind(py).getattr("source")?.is_none());
            drop(error);
            assert!(!retained.bind(py).str()?.is_empty()?);
            assert_eq!(structure.borrow(py).inner.to_mmcif().unwrap(), before);
            assert_eq!(structure.borrow(py).inner.to_pdb().unwrap(), baseline);

            let cases = [
                (ck::BioPdbWriteError::ChainNameTooLong { chain: "AB\0C".into(), length: 4 }, BioPdbWriteErrorKind::ChainNameTooLong, "chain", "AB\0C".into_pyobject(py)?.into_any().unbind()),
                (ck::BioPdbWriteError::NegativeSerial { serial: -1 }, BioPdbWriteErrorKind::NegativeSerial, "serial", (-1_i32).into_pyobject(py)?.into_any().unbind()),
                (ck::BioPdbWriteError::SerialIncrementOverflow { from: i32::MAX }, BioPdbWriteErrorKind::SerialIncrementOverflow, "from_", i32::MAX.into_pyobject(py)?.into_any().unbind()),
                (ck::BioPdbWriteError::SerialOffsetOverflow { serial: i32::MIN }, BioPdbWriteErrorKind::SerialOffsetOverflow, "serial", i32::MIN.into_pyobject(py)?.into_any().unbind()),
                (ck::BioPdbWriteError::NegativeBase36Value { value: -2 }, BioPdbWriteErrorKind::NegativeBase36Value, "value", (-2_i32).into_pyobject(py)?.into_any().unbind()),
                (ck::BioPdbWriteError::UnrepresentableField { field: "atom_name", detail: "counted\0text".into() }, BioPdbWriteErrorKind::UnrepresentableField, "detail", "counted\0text".into_pyobject(py)?.into_any().unbind()),
            ];
            for (native, variant, field, expected) in cases {
                let display = native.to_string();
                let error = BioPdbWriteError::from_rust(native);
                assert!(error.value(py).getattr("variant")?.extract::<BioPdbWriteErrorKind>()? == variant);
                assert!(error.value(py).getattr(field)?.eq(expected.bind(py))?);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(error.value(py).getattr("source")?.is_none());
                assert!(error.value(py).setattr(field, py.None()).is_err());
            }
            let error = BioPdbWriteError::from_rust(ck::BioPdbWriteError::InvalidStructure(ck::BioStructureError::CoordinateCountMismatch { atom_count: 2, coordinate_count: 1 }));
            assert!(error.value(py).getattr("source")?.is_instance_of::<BioStructureError>());
            assert!(error.value(py).getattr("source")?.is(&error.value(py).getattr("__cause__")?));
            Ok(())
        }).unwrap();
    }
}
