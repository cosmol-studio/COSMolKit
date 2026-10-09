//! Thin projections of the facade's valence configuration.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ValenceModel {
    RdkitLike = 0,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct ValenceParams {
    pub(crate) inner: ck::ValenceParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ValenceParams {
    #[new]
    #[pyo3(signature = (model=ValenceModel::RdkitLike, strict=true))]
    fn new(model: ValenceModel, strict: bool) -> Self {
        Self {
            inner: ck::ValenceParams {
                model: match model {
                    ValenceModel::RdkitLike => ck::ValenceModel::RdkitLike,
                },
                strict,
            },
        }
    }

    #[getter]
    fn model(&self) -> ValenceModel {
        match self.inner.model {
            ck::ValenceModel::RdkitLike => ValenceModel::RdkitLike,
        }
    }

    #[getter]
    fn strict(&self) -> bool {
        self.inner.strict
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<ValenceModel>()?;
    module.add_class::<ValenceParams>()?;
    Ok(())
}
