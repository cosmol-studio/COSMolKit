//! Static declaration metadata stays owned and ordered by the Rust facade.
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "cosmolkit", frozen)]
pub struct FeatureSpec {
    pub(crate) inner: &'static ck::FeatureSpec,
}
#[pymethods]
impl FeatureSpec {
    #[getter]
    fn name(&self) -> &'static str {
        self.inner.name
    }
    #[getter]
    fn category(&self) -> &'static str {
        self.inner.category
    }
    #[getter]
    fn docs(&self) -> &'static str {
        self.inner.docs
    }
}

#[pyclass(module = "cosmolkit")]
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
pub fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter {
        inner: ck::feature_specs(),
    }
}
#[pyfunction]
pub fn feature_spec(name: &str) -> Option<FeatureSpec> {
    ck::feature_spec(name).map(|inner| FeatureSpec { inner })
}
