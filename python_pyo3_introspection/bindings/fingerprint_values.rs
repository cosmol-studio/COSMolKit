//! Typed sparse fingerprint values; every operation delegates to the facade.

/// Namespace for the facade's source-defined atom-pair encoding vocabulary.
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct AtomPairsParameters;

#[pymethods]
impl AtomPairsParameters {
    #[staticmethod]
    fn version() -> &'static str {
        ck::AtomPairsParameters::version()
    }
    #[staticmethod]
    fn num_type_bits() -> u32 {
        ck::AtomPairsParameters::num_type_bits()
    }
    #[staticmethod]
    fn num_pi_bits() -> u32 {
        ck::AtomPairsParameters::num_pi_bits()
    }
    #[staticmethod]
    fn num_branch_bits() -> u32 {
        ck::AtomPairsParameters::num_branch_bits()
    }
    #[staticmethod]
    fn num_chiral_bits() -> u32 {
        ck::AtomPairsParameters::num_chiral_bits()
    }
    #[staticmethod]
    fn code_size() -> u32 {
        ck::AtomPairsParameters::code_size()
    }
    #[staticmethod]
    fn num_path_bits() -> u32 {
        ck::AtomPairsParameters::num_path_bits()
    }
    #[staticmethod]
    fn max_path_length() -> u32 {
        ck::AtomPairsParameters::max_path_length()
    }
    #[staticmethod]
    fn num_atom_pair_fingerprint_bits() -> u32 {
        ck::AtomPairsParameters::num_atom_pair_fingerprint_bits()
    }
    #[staticmethod]
    fn atom_types() -> Vec<u32> {
        ck::AtomPairsParameters::atom_types()
    }
}

use crate::errors::FingerprintError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct SparseBitFingerprint {
    pub(crate) inner: ck::SparseBitFingerprint,
}

#[pymethods]
impl SparseBitFingerprint {
    #[new]
    fn new(n_bits: u32) -> Self {
        Self {
            inner: ck::SparseBitFingerprint::new(n_bits),
        }
    }
    fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    fn on_bits(&self) -> Vec<i32> {
        self.inner.on_bits()
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct SparseCountFingerprint {
    pub(crate) inner: ck::SparseCountFingerprint,
}

#[pymethods]
impl SparseCountFingerprint {
    #[staticmethod]
    fn new(length: u64) -> Self {
        Self {
            inner: ck::SparseCountFingerprint::new(length),
        }
    }
    fn length(&self) -> u64 {
        self.inner.length()
    }
    fn value(&self, index: u64) -> PyResult<i32> {
        self.inner.value(index).map_err(FingerprintError::from_rust)
    }
    fn set_value(&mut self, index: u64, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(FingerprintError::from_rust)
    }
    fn nonzero_elements(&self) -> BTreeMap<u64, i32> {
        self.inner.nonzero_elements().clone()
    }
    #[pyo3(signature = (use_abs=false))]
    fn total_value(&self, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(FingerprintError::from_rust)
    }
    fn fuzzy_and(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn fuzzy_or(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_added(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_subtracted(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_added_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_subtracted_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_multiplied_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_divided_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_divided_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct SparseCountFingerprint32 {
    pub(crate) inner: ck::SparseCountFingerprint32,
}

#[pymethods]
impl SparseCountFingerprint32 {
    #[staticmethod]
    fn new(length: u32) -> Self {
        Self {
            inner: ck::SparseCountFingerprint32::new(length),
        }
    }
    fn length(&self) -> u32 {
        self.inner.length()
    }
    fn value(&self, index: u32) -> PyResult<i32> {
        self.inner.value(index).map_err(FingerprintError::from_rust)
    }
    fn set_value(&mut self, index: u32, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(FingerprintError::from_rust)
    }
    fn nonzero_elements(&self) -> BTreeMap<u32, i32> {
        self.inner.nonzero_elements().clone()
    }
    #[pyo3(signature = (use_abs=false))]
    fn total_value(&self, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(FingerprintError::from_rust)
    }
    fn fuzzy_and(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn fuzzy_or(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_added(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_subtracted(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_added_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_subtracted_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_multiplied_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    fn with_divided_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_divided_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
}
