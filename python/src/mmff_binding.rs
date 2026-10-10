//! Thin MMFF projections of canonical public COSMolKit types.
use crate::drawing_binding::{Molecule, operation_pyerr};
use ::cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
fn cause_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    if let Some(source) = source.downcast_ref::<ck::MmffMolPropertiesError>() {
        return properties_pyerr_ref(py, source);
    }
    let error = PyValueError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|next| cause_pyerr(py, next)));
    error
}
pyo3::create_exception!(cosmolkit, MmffMolPropertiesError, PyValueError);
pub(crate) fn properties_pyerr(py: Python<'_>, source: ck::MmffMolPropertiesError) -> PyErr {
    properties_pyerr_ref(py, &source)
}
pub(crate) fn properties_pyerr_ref(py: Python<'_>, source: &ck::MmffMolPropertiesError) -> PyErr {
    let kind = match source {
        ck::MmffMolPropertiesError::Params(_) => "Params",
        ck::MmffMolPropertiesError::Kekulize(_) => "Kekulize",
        ck::MmffMolPropertiesError::Aromaticity(_) => "Aromaticity",
        ck::MmffMolPropertiesError::RingFinding(_) => "RingFinding",
        ck::MmffMolPropertiesError::Valence(_) => "Valence",
        ck::MmffMolPropertiesError::AtomIndexOutOfRange { .. } => "AtomIndexOutOfRange",
        ck::MmffMolPropertiesError::AtomTypePropertiesMissing { .. } => "AtomTypePropertiesMissing",
        ck::MmffMolPropertiesError::AtomTypePbciMissing { .. } => "AtomTypePbciMissing",
    };
    let error = MmffMolPropertiesError::new_err(source.to_string());
    if let Err(attribute_error) = error
        .value(py)
        .setattr("domain", "mmff_properties")
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(source).map(|cause| cause_pyerr(py, cause)),
    );
    error
}
pyo3::create_exception!(cosmolkit, MmffOptimizationError, PyValueError);
pub(crate) fn optimization_pyerr(py: Python<'_>, source: &ck::MmffOptimizationError) -> PyErr {
    let error = MmffOptimizationError::new_err(source.to_string());
    if let Err(attribute_error) = error
        .value(py)
        .setattr("domain", "mmff_optimization")
        .and_then(|()| error.value(py).setattr("kind", "MmffOptimization"))
    {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(source).map(|cause| cause_pyerr(py, cause)),
    );
    error
}
pub(crate) fn source_conformer_id(id: i32) -> Option<usize> {
    // RDKit✔️✔️:   if (id < 0) {
    // RDKit✔️✔️:     return *(d_confs.front());
    // RDKit✔️✔️:   }
    // ROMol.cpp, pinned 351f8f378f8ad6bbd517980c38896e66bf907af8.
    // All signed source selectors below zero designate the first conformer;
    // None is the canonical Rust selection. Constant-time scalar conversion.
    if id < 0 { None } else { Some(id as usize) }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct MmffOptimizationParams {
    pub(crate) inner: ck::MmffOptimizationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffOptimizationParams {
    #[new]
    #[pyo3(signature=(mmff_variant="MMFF94",max_iterations=200,non_bonded_threshold=100.0,conformer_id=None,ignore_interfragment_interactions=true))]
    fn new(
        mmff_variant: &str,
        max_iterations: i32,
        non_bonded_threshold: f64,
        conformer_id: Option<i32>,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        let conformer_id = conformer_id.and_then(source_conformer_id);
        Self {
            inner: ck::MmffOptimizationParams {
                mmff_variant: mmff_variant.into(),
                max_iterations,
                non_bonded_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            },
        }
    }
    #[getter]
    fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
    #[getter]
    fn max_iterations(&self) -> i32 {
        self.inner.max_iterations
    }
    #[getter]
    fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold
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
pub(crate) struct MmffConformerOptimizationParams {
    pub(crate) inner: ck::MmffConformerOptimizationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffConformerOptimizationParams {
    #[new]
    #[pyo3(signature=(num_threads=1,max_iterations=1000,mmff_variant="MMFF94",non_bonded_threshold=10.0,ignore_interfragment_interactions=true))]
    fn new(
        num_threads: i32,
        max_iterations: i32,
        mmff_variant: &str,
        non_bonded_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            inner: ck::MmffConformerOptimizationParams {
                num_threads,
                max_iterations,
                mmff_variant: mmff_variant.into(),
                non_bonded_threshold,
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
    fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
    #[getter]
    fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct MmffPropertiesParams {
    pub(crate) inner: ck::MmffPropertiesParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffPropertiesParams {
    #[new]
    #[pyo3(signature=(mmff_variant="MMFF94"))]
    fn new(mmff_variant: &str) -> Self {
        Self {
            inner: ck::MmffPropertiesParams {
                mmff_variant: mmff_variant.into(),
            },
        }
    }
    #[getter]
    fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MmffOptimizeMoleculeResult {
    pub(crate) inner: ck::MmffOptimizeMoleculeResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffOptimizeMoleculeResult {
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule.clone(),
        }
    }
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn __repr__(&self) -> String {
        format!(
            "MmffOptimizeMoleculeResult(needs_more={})",
            self.inner.needs_more
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MmffOptimizeMoleculeConfResult {
    pub(crate) inner: ck::MmffOptimizeMoleculeConfResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffOptimizeMoleculeConfResult {
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn energy(&self) -> f64 {
        self.inner.energy
    }
    fn __repr__(&self) -> String {
        format!(
            "MmffOptimizeMoleculeConfResult(needs_more={}, energy={})",
            self.inner.needs_more, self.inner.energy
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MmffOptimizeMoleculeConfsResult {
    pub(crate) inner: ck::MmffOptimizeMoleculeConfsResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffOptimizeMoleculeConfsResult {
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule.clone(),
        }
    }
    fn conformer_results(&self) -> Vec<MmffOptimizeMoleculeConfResult> {
        self.inner
            .conformer_results
            .iter()
            .copied()
            .map(|inner| MmffOptimizeMoleculeConfResult { inner })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "MmffOptimizeMoleculeConfsResult(conformers={})",
            self.inner.conformer_results.len()
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MmffProperties {
    pub(crate) inner: ck::MmffProperties,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffProperties {
    fn is_valid(&self) -> bool {
        self.inner.is_valid()
    }
    fn variant(&self) -> MmffVariant {
        match self.inner.variant() {
            ck::MmffVariant::Mmff94 => MmffVariant::Mmff94,
            ck::MmffVariant::Mmff94s => MmffVariant::Mmff94s,
        }
    }
    fn atoms(&self) -> Vec<MmffAtomProperties> {
        self.inner
            .atoms()
            .iter()
            .copied()
            .map(|inner| MmffAtomProperties { inner })
            .collect()
    }
    fn atom_type(&self, py: Python<'_>, atom_index: usize) -> PyResult<u8> {
        self.inner
            .atom_type(atom_index)
            .map_err(|e| properties_pyerr(py, e))
    }
    fn formal_charge(&self, py: Python<'_>, atom_index: usize) -> PyResult<f64> {
        self.inner
            .formal_charge(atom_index)
            .map_err(|e| properties_pyerr(py, e))
    }
    fn partial_charge(&self, py: Python<'_>, atom_index: usize) -> PyResult<f64> {
        self.inner
            .partial_charge(atom_index)
            .map_err(|e| properties_pyerr(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct Conformer3D {
    pub(crate) inner: ck::Conformer3D,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Conformer3D {
    fn id(&self) -> usize {
        self.inner.id()
    }
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
    /// Return an independent float64 NumPy array (N, 3).
    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]", imports = ("numpy", "typing")))]
    fn coordinates<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        crate::canonical_coordinate_input::coordinate_array(py, self.inner.coordinates())
    }
    fn props(&self, py: Python<'_>) -> PyResult<std::collections::BTreeMap<String, String>> {
        self.inner
            .props()
            .iter()
            .map(|(key, value)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, key)?,
                    crate::canonical_sdf::decode_source_text(py, value)?,
                ))
            })
            .collect()
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<MmffVariant>()?;
    module.add_class::<MmffAtomProperties>()?;
    module.add_class::<MmffEvaluationParams>()?;
    module.add_class::<MmffEnergyGradient>()?;
    module.add(
        "MmffMolPropertiesError",
        module.py().get_type::<MmffMolPropertiesError>(),
    )?;
    module.add(
        "MmffOptimizationError",
        module.py().get_type::<MmffOptimizationError>(),
    )?;
    module.add_class::<MmffOptimizationParams>()?;
    module.add_class::<MmffConformerOptimizationParams>()?;
    module.add_class::<MmffPropertiesParams>()?;
    module.add_class::<MmffOptimizeMoleculeResult>()?;
    module.add_class::<MmffOptimizeMoleculeConfResult>()?;
    module.add_class::<MmffOptimizeMoleculeConfsResult>()?;
    module.add_class::<MmffProperties>()?;
    module.add_class::<Conformer3D>()?;
    Ok(())
}
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) enum MmffVariant {
    Mmff94,
    Mmff94s,
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MmffAtomProperties {
    inner: ck::MmffAtomProperties,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffAtomProperties {
    fn atom_type(&self) -> u8 {
        self.inner.atom_type
    }
    fn formal_charge(&self) -> f64 {
        self.inner.formal_charge
    }
    fn partial_charge(&self) -> f64 {
        self.inner.partial_charge
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct MmffEvaluationParams {
    pub(crate) inner: ck::MmffEvaluationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffEvaluationParams {
    #[new]
    #[pyo3(signature=(mmff_variant="MMFF94",non_bonded_threshold=100.0,conformer_id=None,ignore_interfragment_interactions=true))]
    fn new(
        mmff_variant: &str,
        non_bonded_threshold: f64,
        conformer_id: Option<i32>,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        let conformer_id = conformer_id.and_then(source_conformer_id);
        Self {
            inner: ck::MmffEvaluationParams {
                mmff_variant: mmff_variant.into(),
                non_bonded_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            },
        }
    }
    #[getter]
    fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
    #[getter]
    fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold
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
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct MmffEnergyGradient {
    pub(crate) inner: ck::MmffEnergyGradient,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffEnergyGradient {
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    fn gradient(&self) -> Vec<f64> {
        self.inner.gradient().to_vec()
    }
}
