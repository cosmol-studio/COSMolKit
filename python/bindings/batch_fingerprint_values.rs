use crate::error_framework::error_projection;
use crate::native::Fingerprint;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::error::Error;
use std::sync::Arc;

ck::project_configuration!(ck::TopologicalFingerprintParams);

error_projection!(BatchFingerprintOutputError, BatchFingerprintOutputErrorKind, ck::BatchFingerprintOutputError, "fingerprint", {
    MissingAdditionalOutput { fingerprint_kind: &'static str },
}, cause = |_: &ck::BatchFingerprintOutputError| None::<PyErr>);

#[pyclass(module = "cosmolkit", frozen)]
/// One batch record's fingerprint and optional native provenance output.
pub struct BatchFingerprintOutput {
    inner: Arc<ck::BatchFingerprintOutput>,
}

impl BatchFingerprintOutput {
    pub(crate) fn from_rust(inner: ck::BatchFingerprintOutput) -> Self {
        Self {
            inner: Arc::new(inner),
        }
    }
}

#[pymethods]
impl BatchFingerprintOutput {
    #[getter]
    /// Bit fingerprint for this batch record.
    fn fingerprint(&self) -> Fingerprint {
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    #[getter]
    /// Owner-retaining view of collected provenance; raises BatchFingerprintOutputError when provenance was not requested.
    fn additional_output(&self) -> PyResult<BatchFingerprintAdditionalOutput> {
        self.inner
            .additional_output()
            .map_err(BatchFingerprintOutputError::from_rust)?;
        Ok(BatchFingerprintAdditionalOutput {
            owner: Arc::clone(&self.inner),
        })
    }
}

/// Immutable view keeps the actual result alive; constructing it never copies
/// provenance. Individual reads only allocate their requested Python value.
#[pyclass(module = "cosmolkit", frozen)]
pub struct BatchFingerprintAdditionalOutput {
    owner: Arc<ck::BatchFingerprintOutput>,
}

impl BatchFingerprintAdditionalOutput {
    fn inner(&self) -> PyResult<&ck::BatchFingerprintAdditionalOutput> {
        self.owner
            .additional_output()
            .map_err(BatchFingerprintOutputError::from_rust)
    }
}

#[pymethods]
impl BatchFingerprintAdditionalOutput {
    #[getter]
    /// Per-atom feature contribution counts, or None when not collected.
    fn atom_counts(&self) -> PyResult<Option<Vec<u32>>> {
        Ok(self.inner()?.atom_counts().map(<[u32]>::to_vec))
    }
    #[getter]
    /// Contributing fingerprint bit identifiers for each atom, or None when not collected.
    fn atom_to_bits(&self) -> PyResult<Option<Vec<Vec<u64>>>> {
        Ok(self.inner()?.atom_to_bits().map(<[Vec<u64>]>::to_vec))
    }
    #[getter]
    /// Bit identifiers mapped to native feature descriptors (center/radius or atom-pair endpoints, depending on the generator), or None when not collected.
    fn bit_info_map(&self) -> PyResult<Option<BTreeMap<u64, Vec<(u32, u32)>>>> {
        Ok(self.inner()?.bit_info_map().cloned())
    }
    #[getter]
    /// Bit identifiers mapped to contributing native paths, or None when not collected.
    fn bit_paths(&self) -> PyResult<Option<BTreeMap<u64, Vec<Vec<i32>>>>> {
        Ok(self.inner()?.bit_paths().cloned())
    }
    #[getter]
    /// Bit identifiers mapped to participating atom sets, or None when not collected.
    fn atoms_per_bit(&self) -> PyResult<Option<BTreeMap<u64, Vec<Vec<i32>>>>> {
        Ok(self.inner()?.atoms_per_bit().cloned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_batch_provenance_keeps_owner_and_exact_native_results() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let batch = ck::MoleculeBatch::from_smiles_list(&["CCO".to_owned()]).unwrap();
            for collect in [false, true] {
                let mut records = batch
                    .fingerprint_morgan_with_output_with_params(
                        &ck::MorganFingerprintParams::default(),
                        collect,
                        &ck::BatchQueryParams {
                            n_jobs: Some(2),
                            ..Default::default()
                        },
                    )
                    .unwrap();
                let native = records.remove(0).unwrap();
                let expected_bits = native.fingerprint().on_bits();
                let wrapped = BatchFingerprintOutput::from_rust(native);
                assert_eq!(wrapped.fingerprint().inner.on_bits(), expected_bits);
                if collect {
                    let view = wrapped.additional_output()?;
                    assert!(Arc::ptr_eq(&view.owner, &wrapped.inner));
                    let native = wrapped.inner.additional_output().unwrap();
                    assert_eq!(
                        view.atom_counts()?,
                        native.atom_counts().map(<[u32]>::to_vec)
                    );
                    assert_eq!(
                        view.atom_to_bits()?,
                        native.atom_to_bits().map(<[Vec<u64>]>::to_vec)
                    );
                    assert_eq!(view.bit_info_map()?, native.bit_info_map().cloned());
                    assert_eq!(view.bit_paths()?, native.bit_paths().cloned());
                    assert_eq!(view.atoms_per_bit()?, native.atoms_per_bit().cloned());
                    let expected_counts = view.atom_counts()?;
                    drop(wrapped);
                    let view = Py::new(py, view)?;
                    let observed = view
                        .bind(py)
                        .getattr("atom_counts")?
                        .extract::<Option<Vec<u32>>>()?;
                    assert_eq!(observed, expected_counts);
                } else {
                    let error = match wrapped.additional_output() {
                        Ok(_) => panic!("native absent output must not become an empty collector"),
                        Err(error) => error,
                    };
                    assert!(error.is_instance(py, <BatchFingerprintOutputError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                    assert_eq!(
                        error
                            .value(py)
                            .getattr("fingerprint_kind")?
                            .extract::<String>()?,
                        "Morgan"
                    );
                    assert!(error.cause(py).is_none());
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
