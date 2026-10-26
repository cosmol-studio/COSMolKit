//! Static declaration metadata stays owned and ordered by the Rust facade.
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "cosmolkit", frozen)]
/// Read-only public feature declaration from the compiled Rust facade.
pub struct FeatureSpec {
    pub(crate) inner: &'static ck::FeatureSpec,
}
#[pymethods]
impl FeatureSpec {
    #[getter]
    /// Canonical public feature name.
    fn name(&self) -> &'static str {
        self.inner.name
    }
    #[getter]
    /// Category assigned to this feature declaration.
    fn category(&self) -> &'static str {
        self.inner.category
    }
    #[getter]
    /// Documentation text retained by the canonical feature declaration.
    fn docs(&self) -> &'static str {
        self.inner.docs
    }
}

#[pyclass(module = "cosmolkit")]
/// Iterator over the compiled feature declarations in native declaration order.
pub struct FeatureSpecIter {
    inner: ck::FeatureSpecIter,
}
#[pymethods]
impl FeatureSpecIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self) -> Option<FeatureSpec> {
        self.inner.next().map(|inner| FeatureSpec { inner })
    }
}

#[pyfunction]
/// Iterate the canonical public feature declarations compiled into this library.
pub fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter {
        inner: ck::feature_specs(),
    }
}
#[pyfunction]
/// Look up a public feature declaration by its exact name; return None for an unknown feature.
pub fn feature_spec(name: &str) -> Option<FeatureSpec> {
    ck::feature_spec(name).map(|inner| FeatureSpec { inner })
}
