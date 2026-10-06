//! Frozen batch fingerprint values project the canonical read-only owner snapshots.
use crate::canonical_values::Fingerprint;
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;
pyo3::create_exception!(cosmolkit, BatchFingerprintOutputError, PyValueError);
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct BatchFingerprintAdditionalOutput {
    pub(crate) inner: ck::BatchFingerprintAdditionalOutput,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchFingerprintAdditionalOutput {
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts().map(<[u32]>::to_vec)
    }
    fn atom_to_bits(&self) -> Option<Vec<Vec<u64>>> {
        self.inner.atom_to_bits().map(<[Vec<u64>]>::to_vec)
    }
    fn bit_info_map(&self) -> Option<BTreeMap<u64, Vec<(u32, u32)>>> {
        self.inner.bit_info_map().cloned()
    }
    fn bit_paths(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.bit_paths().cloned()
    }
    fn atoms_per_bit(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.atoms_per_bit().cloned()
    }
    fn __repr__(&self) -> String {
        format!(
            "BatchFingerprintAdditionalOutput(atom_counts={},atom_to_bits={},bit_info_map={},bit_paths={},atoms_per_bit={})",
            self.inner.atom_counts().is_some(),
            self.inner.atom_to_bits().is_some(),
            self.inner.bit_info_map().is_some(),
            self.inner.bit_paths().is_some(),
            self.inner.atoms_per_bit().is_some()
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct BatchFingerprintOutput {
    pub(crate) inner: ck::BatchFingerprintOutput,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchFingerprintOutput {
    fn fingerprint(&self) -> Fingerprint {
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    fn additional_output(&self) -> PyResult<BatchFingerprintAdditionalOutput> {
        self.inner
            .additional_output()
            .map(|inner| BatchFingerprintAdditionalOutput {
                inner: inner.clone(),
            })
            .map_err(|source| BatchFingerprintOutputError::new_err(source.to_string()))
    }
    fn __repr__(&self) -> String {
        format!(
            "BatchFingerprintOutput(n_bits={},has_additional_output={})",
            self.inner.fingerprint().n_bits(),
            self.inner.additional_output.is_some()
        )
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<BatchFingerprintOutput>()?;
    module.add_class::<BatchFingerprintAdditionalOutput>()?;
    module.add(
        "BatchFingerprintOutputError",
        module.py().get_type::<BatchFingerprintOutputError>(),
    )?;
    Ok(())
}
