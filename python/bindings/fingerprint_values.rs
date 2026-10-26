//! Typed sparse fingerprint values; every operation delegates to the facade.

/// Namespace for the facade's source-defined atom-pair encoding vocabulary.
#[pyclass(module = "cosmolkit", frozen)]
pub struct AtomPairsParameters;

#[pymethods]
impl AtomPairsParameters {
    #[staticmethod]
    /// Version of the native atom-pair encoding scheme.
    fn version() -> &'static str {
        ck::AtomPairsParameters::version()
    }
    #[staticmethod]
    /// Number of bits reserved for the atom type in an atom code.
    fn num_type_bits() -> u32 {
        ck::AtomPairsParameters::num_type_bits()
    }
    #[staticmethod]
    /// Number of bits reserved for pi-electron information in an atom code.
    fn num_pi_bits() -> u32 {
        ck::AtomPairsParameters::num_pi_bits()
    }
    #[staticmethod]
    /// Number of bits reserved for atom branching information.
    fn num_branch_bits() -> u32 {
        ck::AtomPairsParameters::num_branch_bits()
    }
    #[staticmethod]
    /// Number of bits reserved for chirality information.
    fn num_chiral_bits() -> u32 {
        ck::AtomPairsParameters::num_chiral_bits()
    }
    #[staticmethod]
    /// Number of bits in the native atom code.
    fn code_size() -> u32 {
        ck::AtomPairsParameters::code_size()
    }
    #[staticmethod]
    /// Number of bits reserved for pair path lengths.
    fn num_path_bits() -> u32 {
        ck::AtomPairsParameters::num_path_bits()
    }
    #[staticmethod]
    /// Largest path length representable by the atom-pair encoding.
    fn max_path_length() -> u32 {
        ck::AtomPairsParameters::max_path_length()
    }
    #[staticmethod]
    /// Size of the native atom-pair fingerprint identifier space.
    fn num_atom_pair_fingerprint_bits() -> u32 {
        ck::AtomPairsParameters::num_atom_pair_fingerprint_bits()
    }
    #[staticmethod]
    /// Atomic numbers represented by the atom-pair type encoding, in encoding order.
    fn atom_types() -> Vec<u32> {
        ck::AtomPairsParameters::atom_types()
    }
}

use crate::errors::FingerprintError;
use cosmolkit as ck;

crate::configuration::type_projection! {
    #[pyo3::pyclass(module = "cosmolkit", frozen)]
    /// Immutable bit fingerprint with an explicit bit-index capacity.
    /// on_bits() returns the set positions; to_numpy() expands bits to an
    /// independent uint8 array. Similarity and bit operations use native values.
    pub struct Fingerprint {
        pub(crate) inner: ck::Fingerprint,
    }
}
use pyo3::prelude::*;
use std::collections::BTreeMap;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Immutable sparse bit vector storing set positions without a dense bit buffer.
    pub struct SparseBitFingerprint {
        pub(crate) inner: ck::SparseBitFingerprint,
    }
}

#[pymethods]
impl SparseBitFingerprint {
    #[new]
    /// Create an empty sparse bit vector with the given index-space capacity.
    fn new(n_bits: u32) -> Self {
        Self {
            inner: ck::SparseBitFingerprint::new(n_bits),
        }
    }
    #[getter]
    /// Declared capacity of the sparse bit fingerprint.
    fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    #[getter]
    /// Set bit indices returned as an independent list.
    fn on_bits(&self) -> Vec<i32> {
        self.inner.on_bits()
    }
}

#[pyclass(module = "cosmolkit")]
/// Mutable sparse integer counts with 64-bit indices and explicit index-space length.
pub struct SparseCountFingerprint {
    pub(crate) inner: ck::SparseCountFingerprint,
}

#[pymethods]
impl SparseCountFingerprint {
    #[new]
    fn py_new(length: u64) -> Self {
        Self::new(length)
    }
    #[staticmethod]
    /// Create a sparse count vector with no nonzero entries and the given index-space length.
    fn new(length: u64) -> Self {
        Self {
            inner: ck::SparseCountFingerprint::new(length),
        }
    }
    #[getter]
    /// Declared index-space length, not the number of stored nonzero entries.
    fn length(&self) -> u64 {
        self.inner.length()
    }
    /// Return the count at an index; absent entries are zero and out-of-range indices raise FingerprintError.
    fn value(&self, index: u64) -> PyResult<i32> {
        self.inner.value(index).map_err(FingerprintError::from_rust)
    }
    /// Assign a count at an index; zero removes the stored entry. Out-of-range indices raise FingerprintError.
    fn set_value(&mut self, index: u64, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(FingerprintError::from_rust)
    }
    #[getter]
    /// Independent dictionary mapping nonzero entry indices to their counts.
    fn nonzero_elements(&self) -> BTreeMap<u64, i32> {
        self.inner.nonzero_elements().clone()
    }
    #[pyo3(signature = (use_abs=false))]
    /// Sum stored counts, optionally using their absolute values; native arithmetic errors propagate.
    fn total_value(&self, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint using the native componentwise fuzzy intersection.
    fn fuzzy_and(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint using the native componentwise fuzzy union.
    fn fuzzy_or(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint with corresponding counts added; both inputs remain unchanged.
    fn with_added(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint with corresponding counts subtracted; both inputs remain unchanged.
    fn with_subtracted(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after adding the scalar using native sparse-vector semantics.
    fn with_added_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after subtracting the scalar using native sparse-vector semantics.
    fn with_subtracted_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after multiplying counts by the scalar.
    fn with_multiplied_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after dividing counts by the scalar; native arithmetic errors propagate.
    fn with_divided_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_divided_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
}

#[pyclass(module = "cosmolkit")]
/// Mutable sparse integer counts with 32-bit indices and explicit index-space length.
pub struct SparseCountFingerprint32 {
    pub(crate) inner: ck::SparseCountFingerprint32,
}

#[pymethods]
impl SparseCountFingerprint32 {
    #[new]
    fn py_new(length: u32) -> Self {
        Self::new(length)
    }
    #[staticmethod]
    /// Create a sparse count vector with no nonzero entries and the given index-space length.
    fn new(length: u32) -> Self {
        Self {
            inner: ck::SparseCountFingerprint32::new(length),
        }
    }
    #[getter]
    /// Declared index-space length, not the number of stored nonzero entries.
    fn length(&self) -> u32 {
        self.inner.length()
    }
    /// Return the count at an index; absent entries are zero and out-of-range indices raise FingerprintError.
    fn value(&self, index: u32) -> PyResult<i32> {
        self.inner.value(index).map_err(FingerprintError::from_rust)
    }
    /// Assign a count at an index; zero removes the stored entry. Out-of-range indices raise FingerprintError.
    fn set_value(&mut self, index: u32, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(FingerprintError::from_rust)
    }
    #[getter]
    /// Independent dictionary mapping nonzero entry indices to their counts.
    fn nonzero_elements(&self) -> BTreeMap<u32, i32> {
        self.inner.nonzero_elements().clone()
    }
    #[pyo3(signature = (use_abs=false))]
    /// Sum stored counts, optionally using their absolute values; native arithmetic errors propagate.
    fn total_value(&self, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint using the native componentwise fuzzy intersection.
    fn fuzzy_and(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint using the native componentwise fuzzy union.
    fn fuzzy_or(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint with corresponding counts added; both inputs remain unchanged.
    fn with_added(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint with corresponding counts subtracted; both inputs remain unchanged.
    fn with_subtracted(&self, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after adding the scalar using native sparse-vector semantics.
    fn with_added_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after subtracting the scalar using native sparse-vector semantics.
    fn with_subtracted_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after multiplying counts by the scalar.
    fn with_multiplied_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
    /// Return a new fingerprint after dividing counts by the scalar; native arithmetic errors propagate.
    fn with_divided_scalar(&self, value: i32) -> PyResult<Self> {
        self.inner
            .with_divided_scalar(value)
            .map(|inner| Self { inner })
            .map_err(FingerprintError::from_rust)
    }
}
