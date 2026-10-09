//! Thin UFF projections; all chemistry stays in the public Rust facade.
use crate::drawing_binding::{Molecule, operation_pyerr};
use ::cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, UffOptimizationError, PyValueError);
pyo3::create_exception!(cosmolkit, UffParameterQueryError, PyValueError);
pyo3::create_exception!(cosmolkit, UffParameterError, PyValueError);
fn cause_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    if let Some(source) = source.downcast_ref::<ck::OperationError>() {
        return operation_pyerr(py, source.clone());
    }
    let error = PyValueError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|next| cause_pyerr(py, next)));
    error
}
pub(crate) fn optimization_pyerr(py: Python<'_>, source: &ck::UffOptimizationError) -> PyErr {
    let error = UffOptimizationError::new_err(source.to_string());
    let (_, requested) = match source.kind() {
        ck::UffOptimizationErrorKind::MissingConformer { requested } => {
            ("MissingConformer", requested)
        }
        ck::UffOptimizationErrorKind::Evaluation => ("Evaluation", None),
        ck::UffOptimizationErrorKind::Rings => ("Rings", None),
        ck::UffOptimizationErrorKind::Optimization => ("Optimization", None),
        ck::UffOptimizationErrorKind::ConformerOptimization => ("ConformerOptimization", None),
    };
    let kind = match Py::new(
        py,
        crate::canonical_error_values::UffOptimizationErrorKind {
            inner: source.kind(),
        },
    ) {
        Ok(kind) => kind,
        Err(e) => return e,
    };
    if let Err(e) = error
        .value(py)
        .setattr("domain", "uff_optimization")
        .and_then(|()| error.value(py).setattr("_kind", kind))
        .and_then(|()| error.value(py).setattr("requested", requested))
    {
        return e;
    }
    error.set_cause(
        py,
        std::error::Error::source(source).map(|next| cause_pyerr(py, next)),
    );
    error
}
pub(crate) fn parameter_query_pyerr(py: Python<'_>, source: ck::UffParameterQueryError) -> PyErr {
    let error = UffParameterQueryError::new_err(source.to_string());
    let (kind, cause) = match &source {
        ck::UffParameterQueryError::Cache(cause) => ("Cache", operation_pyerr(py, cause.clone())),
        ck::UffParameterQueryError::Parameters(cause) => {
            let error = UffParameterError::new_err(cause.to_string());
            let kind = match Py::new(
                py,
                crate::canonical_error_values::UffParameterErrorKind::from(cause.kind()),
            ) {
                Ok(kind) => kind,
                Err(e) => return e,
            };
            if let Err(e) = error.value(py).setattr("_kind", kind) {
                return e;
            }
            error.set_cause(
                py,
                std::error::Error::source(cause).map(|next| cause_pyerr(py, next)),
            );
            ("Parameters", error)
        }
    };
    if let Err(e) = error
        .value(py)
        .setattr("domain", "uff_parameters")
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return e;
    }
    error.set_cause(py, Some(cause));
    error
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct UffOptimizationParams {
    pub(crate) inner: ck::UffOptimizationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffOptimizationParams {
    #[new]
    #[pyo3(signature=(max_iterations=1000,vdw_threshold=10.0,conformer_id=None,ignore_interfragment_interactions=true))]
    fn new(
        max_iterations: i32,
        vdw_threshold: f64,
        conformer_id: Option<i32>,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        // RDKit✔️✔️:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
        // RDKit✔️✔️:     bool ignoreInterfragInteractions = true) {
        // Source defaults match the registered Rust parameter value. Scalar
        // projection has constant cost and performs no chemistry work.
        let conformer_id = conformer_id.and_then(crate::mmff_binding::source_conformer_id);
        Self {
            inner: ck::UffOptimizationParams {
                max_iterations,
                vdw_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            },
        }
    }
    #[getter]
    fn max_iterations(&self) -> i32 {
        self.inner.max_iterations
    }
    #[getter]
    fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct UffConformerOptimizationParams {
    pub(crate) inner: ck::UffConformerOptimizationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffConformerOptimizationParams {
    #[new]
    #[pyo3(signature=(num_threads=1,max_iterations=1000,vdw_threshold=10.0,ignore_interfragment_interactions=true))]
    fn new(
        num_threads: i32,
        max_iterations: i32,
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        // RDKit✔️✔️:     ROMol &mol, int maxIters = 1000, double vdwThresh = 10.0, int confId = -1,
        // RDKit✔️✔️:     bool ignoreInterfragInteractions = true) {
        // Source defaults match the registered Rust parameter value. Scalar
        // projection has constant cost and performs no chemistry work.
        Self {
            inner: ck::UffConformerOptimizationParams {
                num_threads,
                max_iterations,
                vdw_threshold,
                ignore_interfragment_interactions,
            },
        }
    }
    #[getter]
    fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }
    #[getter]
    fn max_iterations(&self) -> i32 {
        self.inner.max_iterations
    }
    #[getter]
    fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct UffOptimizationResult {
    pub(crate) inner: ck::UffOptimizationResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffOptimizationResult {
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    fn __repr__(&self) -> String {
        format!(
            "UffOptimizationResult(needs_more={}, energy={})",
            self.inner.status_code(),
            self.inner.energy()
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct UffConformerResult {
    pub(crate) inner: ck::UffConformerResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffConformerResult {
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    fn conformer_id(&self) -> usize {
        self.inner.conformer_id()
    }
    fn __repr__(&self) -> String {
        format!(
            "UffConformerResult(needs_more={}, energy={})",
            self.inner.status_code(),
            self.inner.energy()
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct UffConformerOptimizationResult {
    pub(crate) inner: ck::UffConformerOptimizationResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffConformerOptimizationResult {
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conformer_results(&self) -> Vec<UffConformerResult> {
        self.inner
            .conformer_results()
            .iter()
            .copied()
            .map(|inner| UffConformerResult { inner })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "UffConformerOptimizationResult(conformers={})",
            self.inner.conformers.len()
        )
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<crate::canonical_error_values::UffParameterErrorKind>()?;
    module.add_class::<crate::canonical_error_values::UffOptimizationErrorKind>()?;
    crate::canonical_error_accessors::attach(
        module.py().get_type::<UffParameterError>().as_any(),
        &[("kind", "_kind")],
    )?;
    crate::canonical_error_accessors::attach(
        module.py().get_type::<UffOptimizationError>().as_any(),
        &[("kind", "_kind")],
    )?;
    module.add_class::<UffEvaluationParams>()?;
    module.add_class::<UffEnergyGradient>()?;
    module.add_class::<UffOptimizationParams>()?;
    module.add_class::<UffConformerOptimizationParams>()?;
    module.add_class::<UffOptimizationResult>()?;
    module.add_class::<UffConformerResult>()?;
    module.add_class::<UffConformerOptimizationResult>()?;
    module.add(
        "UffOptimizationError",
        module.py().get_type::<UffOptimizationError>(),
    )?;
    module.add(
        "UffParameterQueryError",
        module.py().get_type::<UffParameterQueryError>(),
    )?;
    module.add(
        "UffParameterError",
        module.py().get_type::<UffParameterError>(),
    )?;
    Ok(())
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct UffEvaluationParams {
    pub(crate) inner: ck::UffEvaluationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffEvaluationParams {
    #[new]
    #[pyo3(signature=(vdw_threshold=10.0,conformer_id=None,ignore_interfragment_interactions=true))]
    fn new(
        vdw_threshold: f64,
        conformer_id: Option<i32>,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            inner: ck::UffEvaluationParams {
                vdw_threshold,
                conformer_id: conformer_id.and_then(crate::mmff_binding::source_conformer_id),
                ignore_interfragment_interactions,
            },
        }
    }
    #[getter]
    fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct UffEnergyGradient {
    pub(crate) inner: ck::UffEnergyGradient,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffEnergyGradient {
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    fn gradient(&self) -> Vec<f64> {
        self.inner.gradient().to_vec()
    }
}
