//! Checked Python projections for Gasteiger options and detached atom rows.

use cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyclass(module = "cosmolkit", from_py_object)]
#[derive(Clone)]
/// Checked settings for the Gasteiger partial-charge iteration.
pub struct GasteigerChargeParams {
    pub(crate) inner: ck::GasteigerChargeParams,
}

#[pymethods]
impl GasteigerChargeParams {
    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        crate::representation::configuration(slf.as_any(), "GasteigerChargeParams", &[])
    }

    #[new]
    #[pyo3(signature = (iterations=12, throw_on_param_failure=false))]
    pub(crate) fn new(iterations: i64, throw_on_param_failure: bool) -> PyResult<Self> {
        Ok(Self {
            inner: ck::GasteigerChargeParams::new(iterations, throw_on_param_failure)
                .map_err(|error| PyValueError::new_err(error.to_string()))?,
        })
    }

    #[getter]
    /// Number of charge-equalization iterations, from 0 through 2147483647.
    fn iterations(&self) -> u32 {
        self.inner.iterations()
    }

    #[setter]
    fn set_iterations(&mut self, value: i64) -> PyResult<()> {
        self.inner
            .set_iterations(value)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[getter]
    /// Whether to raise an error when an atom has no Gasteiger parameters.
    fn throw_on_param_failure(&self) -> bool {
        self.inner.throw_on_param_failure()
    }

    #[setter]
    fn set_throw_on_param_failure(&mut self, value: bool) {
        self.inner.set_throw_on_param_failure(value);
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached Gasteiger charge values for one graph atom.
    pub struct GasteigerCharge {
        pub(crate) inner: ck::GasteigerCharge,
    }
}

#[pymethods]
impl GasteigerCharge {
    #[getter]
    /// Charge on this graph atom, in elementary-charge units.
    fn atom_charge(&self) -> f64 {
        self.inner.atom_charge()
    }

    #[getter]
    /// Sum of charges on this atom's attached hydrogens that are not graph vertices.
    fn attached_hydrogen_charge(&self) -> f64 {
        self.inner.attached_hydrogen_charge()
    }
}
