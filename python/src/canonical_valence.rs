//! Thin projections of the facade's valence configuration.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

/// Valence calculation model used when assigning atom valence.
///
/// Declared values: ``RdkitLike``.
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ValenceModel {
    RdkitLike = 0,
}

/// Writable configuration for valence assignment.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ValenceParams {
    pub(crate) inner: ck::ValenceParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ValenceParams {
    /// Configure valence assignment; omitted fields use the defaults shown in the signature.
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

    /// ValenceModel selecting the valence assignment rules.
    #[getter]
    fn model(&self) -> ValenceModel {
        match self.inner.model {
            ck::ValenceModel::RdkitLike => ValenceModel::RdkitLike,
        }
    }

    /// Whether invalid valence assignments are rejected rather than retained.
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
