//! Thin projections of owned public Rust force fields. No chemistry owner here.
use crate::drawing_binding::operation_pyerr;
use ::cosmolkit as ck;
use numpy::{
    AllowTypeChange, IntoPyArray, PyArray2, PyArrayLike,
    ndarray::{Array2, IxDyn},
};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyTuple};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MolecularForceFieldErrorKind {
    MissingConformer,
    Parameterization,
    InvalidParameterization,
    InvalidAtomIndex,
    InvalidFixedAtom,
    CoordinateCount,
    CoordinateShape,
    NonFiniteCoordinate,
    InvalidTolerance,
    Preparation,
    Construction,
    Initialization,
    Rings,
    Kernel,
}
impl From<ck::MolecularForceFieldErrorKind> for MolecularForceFieldErrorKind {
    fn from(kind: ck::MolecularForceFieldErrorKind) -> Self {
        match kind {
            ck::MolecularForceFieldErrorKind::MissingConformer => Self::MissingConformer,
            ck::MolecularForceFieldErrorKind::Parameterization => Self::Parameterization,
            ck::MolecularForceFieldErrorKind::InvalidParameterization => {
                Self::InvalidParameterization
            }
            ck::MolecularForceFieldErrorKind::InvalidAtomIndex => Self::InvalidAtomIndex,
            ck::MolecularForceFieldErrorKind::InvalidFixedAtom => Self::InvalidFixedAtom,
            ck::MolecularForceFieldErrorKind::CoordinateCount => Self::CoordinateCount,
            ck::MolecularForceFieldErrorKind::CoordinateShape => Self::CoordinateShape,
            ck::MolecularForceFieldErrorKind::NonFiniteCoordinate => Self::NonFiniteCoordinate,
            ck::MolecularForceFieldErrorKind::InvalidTolerance => Self::InvalidTolerance,
            ck::MolecularForceFieldErrorKind::Preparation => Self::Preparation,
            ck::MolecularForceFieldErrorKind::Construction => Self::Construction,
            ck::MolecularForceFieldErrorKind::Initialization => Self::Initialization,
            ck::MolecularForceFieldErrorKind::Rings => Self::Rings,
            ck::MolecularForceFieldErrorKind::Kernel => Self::Kernel,
        }
    }
}

pyo3::create_exception!(cosmolkit, ForceFieldError, PyValueError);
pyo3::create_exception!(cosmolkit, MmffForceFieldError, PyValueError);
pyo3::create_exception!(cosmolkit, UffForceFieldError, PyValueError);
fn cause_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    if let Some(source) = source.downcast_ref::<ck::ForceFieldError>() {
        return force_pyerr(py, source);
    }
    if let Some(source) = source.downcast_ref::<ck::OperationError>() {
        return operation_pyerr(py, source.clone());
    }
    if let Some(source) = source.downcast_ref::<ck::MmffMolPropertiesError>() {
        return crate::mmff_binding::properties_pyerr_ref(py, source);
    }
    let error = PyValueError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|cause| cause_pyerr(py, cause)));
    error
}
fn annotate(
    py: Python<'_>,
    error: PyErr,
    kind: ck::MolecularForceFieldErrorKind,
    atom_index: Option<usize>,
    component: Option<usize>,
    actual: Option<usize>,
    expected: Option<usize>,
) -> PyErr {
    let value = error.value(py);
    if let Err(error) = value
        .setattr("domain", "persistent_forcefield")
        .and_then(|()| value.setattr("_requested", py.None()))
        .and_then(|()| value.setattr("_kind", MolecularForceFieldErrorKind::from(kind)))
        .and_then(|()| value.setattr("_atom_index", atom_index))
        .and_then(|()| value.setattr("_component", component))
        .and_then(|()| value.setattr("_actual", actual))
        .and_then(|()| value.setattr("_expected", expected))
    {
        return error;
    }
    error
}
fn force_pyerr(py: Python<'_>, source: &ck::ForceFieldError) -> PyErr {
    let error = annotate(
        py,
        ForceFieldError::new_err(source.to_string()),
        source.kind(),
        source.atom_index(),
        source.component(),
        source.actual(),
        source.expected(),
    );
    if let Err(attribute_error) = error.value(py).setattr("_requested", source.requested()) {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(source).map(|cause| cause_pyerr(py, cause)),
    );
    error
}
pub(crate) fn mmff_pyerr(py: Python<'_>, source: ck::MmffForceFieldError) -> PyErr {
    let error = annotate(
        py,
        MmffForceFieldError::new_err(source.to_string()),
        source.kind(),
        source.atom_index(),
        source.component(),
        source.actual(),
        source.expected(),
    );
    if let Err(attribute_error) = error.value(py).setattr("_requested", source.requested()) {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source).map(|cause| cause_pyerr(py, cause)),
    );
    error
}
pub(crate) fn uff_pyerr(py: Python<'_>, source: ck::UffForceFieldError) -> PyErr {
    let error = UffForceFieldError::new_err(source.to_string());
    let value = error.value(py);
    if let Err(error) = value
        .setattr("domain", "persistent_forcefield")
        .and_then(|()| value.setattr("_kind", MolecularForceFieldErrorKind::from(source.kind())))
        .and_then(|()| value.setattr("_requested", source.requested()))
    {
        return error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source).map(|cause| cause_pyerr(py, cause)),
    );
    error
}
fn input_pyerr(
    py: Python<'_>,
    kind: ck::MolecularForceFieldErrorKind,
    cause: Option<PyErr>,
    actual: Option<usize>,
    expected: Option<usize>,
) -> PyErr {
    let error = annotate(
        py,
        ForceFieldError::new_err(format!("{kind:?}: invalid persistent force field input")),
        kind,
        None,
        None,
        actual,
        expected,
    );
    error.set_cause(py, cause);
    error
}
fn parse_atom_id(
    py: Python<'_>,
    value: &Bound<'_, PyAny>,
    kind: ck::MolecularForceFieldErrorKind,
) -> PyResult<ck::AtomId> {
    let index = value.extract::<usize>().map_err(|cause| {
        let error = input_pyerr(py, kind, Some(cause), None, None);
        // The registered atom_index property is int | None. Preserve negative
        // and oversized Python integers, while a noninteger input stays in the
        // conversion cause instead of violating that property's type.
        if !value.is_instance_of::<pyo3::types::PyInt>() {
            return error;
        }
        match error.value(py).setattr("_atom_index", value) {
            Ok(()) => error,
            Err(attribute_error) => attribute_error,
        }
    })?;
    Ok(ck::AtomId::new(index))
}
fn position_row(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<[f64; 3]> {
    let values = value.extract::<Vec<f64>>().map_err(|cause| {
        input_pyerr(
            py,
            ck::MolecularForceFieldErrorKind::CoordinateShape,
            Some(cause),
            None,
            Some(3),
        )
    })?;
    let actual = values.len();
    values.try_into().map_err(|_| {
        input_pyerr(
            py,
            ck::MolecularForceFieldErrorKind::CoordinateShape,
            None,
            Some(actual),
            Some(3),
        )
    })
}
fn position_rows(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<Vec<[f64; 3]>> {
    let converted = value
        .extract::<PyArrayLike<'_, f64, IxDyn, AllowTypeChange>>()
        .map_err(|cause| {
            input_pyerr(
                py,
                ck::MolecularForceFieldErrorKind::CoordinateShape,
                Some(cause),
                None,
                Some(3),
            )
        })?;
    let array = converted.as_array();
    let shape = array.shape();
    if shape.len() != 2 {
        return Err(input_pyerr(
            py,
            ck::MolecularForceFieldErrorKind::CoordinateShape,
            None,
            Some(shape.len()),
            Some(2),
        ));
    }
    if shape[1] != 3 {
        return Err(input_pyerr(
            py,
            ck::MolecularForceFieldErrorKind::CoordinateShape,
            None,
            Some(shape[1]),
            Some(3),
        ));
    }
    Ok(array
        .outer_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect())
}
fn snapshot<'py>(py: Python<'py>, rows: &[[f64; 3]]) -> Bound<'py, PyArray2<f64>> {
    let mut output = Array2::zeros((rows.len(), 3));
    for (i, row) in rows.iter().enumerate() {
        for j in 0..3 {
            output[[i, j]] = row[j];
        }
    }
    output.into_pyarray(py)
}
// Every omitted keyword obtains its value from the same registered Rust params.
pub(crate) fn default_conformer() -> Option<usize> {
    ck::MmffForceFieldParams::default().conformer_id()
}
pub(crate) fn default_variant() -> String {
    ck::MmffForceFieldParams::default().mmff_variant().into()
}
pub(crate) fn default_non_bonded() -> f64 {
    ck::MmffForceFieldParams::default().non_bonded_threshold()
}
pub(crate) fn default_mmff_ignore() -> bool {
    ck::MmffForceFieldParams::default().ignore_interfragment_interactions()
}
pub(crate) fn default_uff_conformer() -> Option<usize> {
    ck::UffForceFieldParams::default().conformer_id()
}
pub(crate) fn default_vdw() -> f64 {
    ck::UffForceFieldParams::default().vdw_threshold()
}
pub(crate) fn default_uff_ignore() -> bool {
    ck::UffForceFieldParams::default().ignore_interfragment_interactions()
}
pub(crate) fn default_iterations() -> u32 {
    ck::ForceFieldMinimizeParams::default().max_iterations()
}
pub(crate) fn default_energy_tolerance() -> f64 {
    ck::ForceFieldMinimizeParams::default().energy_tolerance()
}
fn default_force_tolerance() -> f64 {
    ck::ForceFieldMinimizeParams::default().force_tolerance()
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct MmffForceFieldParams {
    pub(crate) inner: ck::MmffForceFieldParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MmffForceFieldParams {
    #[new]
    #[pyo3(signature=(*,conformer_id=default_conformer(),mmff_variant=default_variant(),non_bonded_threshold=default_non_bonded(),ignore_interfragment_interactions=default_mmff_ignore()))]
    fn new(
        conformer_id: Option<usize>,
        mmff_variant: String,
        non_bonded_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            inner: ck::MmffForceFieldParams::new(
                conformer_id,
                mmff_variant,
                non_bonded_threshold,
                ignore_interfragment_interactions,
            ),
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id()
    }
    #[getter]
    fn mmff_variant(&self) -> String {
        self.inner.mmff_variant().into()
    }
    #[getter]
    fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold()
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct UffForceFieldParams {
    pub(crate) inner: ck::UffForceFieldParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl UffForceFieldParams {
    #[new]
    #[pyo3(signature=(*,conformer_id=default_uff_conformer(),vdw_threshold=default_vdw(),ignore_interfragment_interactions=default_uff_ignore()))]
    fn new(
        conformer_id: Option<usize>,
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            inner: ck::UffForceFieldParams::new(
                conformer_id,
                vdw_threshold,
                ignore_interfragment_interactions,
            ),
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id()
    }
    #[getter]
    fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold()
    }
    #[getter]
    fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct ForceFieldMinimizeParams {
    pub(crate) inner: ck::ForceFieldMinimizeParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ForceFieldMinimizeParams {
    #[new]
    #[pyo3(signature=(*,max_iterations=default_iterations(),force_tolerance=default_force_tolerance(),energy_tolerance=default_energy_tolerance()))]
    fn new(max_iterations: u32, force_tolerance: f64, energy_tolerance: f64) -> Self {
        Self {
            inner: ck::ForceFieldMinimizeParams::new(
                max_iterations,
                force_tolerance,
                energy_tolerance,
            ),
        }
    }
    #[getter]
    fn energy_tolerance(&self) -> f64 {
        self.inner.energy_tolerance()
    }
    #[getter]
    fn max_iterations(&self) -> u32 {
        self.inner.max_iterations()
    }
    #[getter]
    fn force_tolerance(&self) -> f64 {
        self.inner.force_tolerance()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct ForceFieldEnergyGradient {
    inner: ck::ForceFieldEnergyGradient,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ForceFieldEnergyGradient {
    #[getter]
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    #[gen_stub(override_return_type(type_repr="numpy.typing.NDArray[numpy.float64]",imports=("numpy","numpy.typing")))]
    #[getter]
    fn gradient<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        snapshot(py, self.inner.gradient())
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct ForceFieldMinimizeOutcome {
    inner: ck::ForceFieldMinimizeOutcome,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ForceFieldMinimizeOutcome {
    #[getter]
    fn converged(&self) -> bool {
        self.inner.converged()
    }
    #[getter]
    fn iterations(&self) -> u32 {
        self.inner.iterations()
    }
    #[getter]
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MolecularForceField {
    pub(crate) inner: ck::MolecularForceField,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MolecularForceField {
    fn position(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "builtins.int"))] atom_id: &Bound<'_, PyAny>,
    ) -> PyResult<(f64, f64, f64)> {
        let row = self
            .inner
            .position(parse_atom_id(
                py,
                atom_id,
                ck::MolecularForceFieldErrorKind::InvalidAtomIndex,
            )?)
            .map_err(|e| force_pyerr(py, &e))?;
        Ok((row[0], row[1], row[2]))
    }
    fn set_position_(
        &mut self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "builtins.int"))] atom_id: &Bound<'_, PyAny>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.float] | numpy.typing.NDArray[numpy.float64]",imports=("typing","numpy","numpy.typing")))]
        position: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let id = parse_atom_id(
            py,
            atom_id,
            ck::MolecularForceFieldErrorKind::InvalidAtomIndex,
        )?;
        let row = position_row(py, position)?;
        self.inner
            .set_position_(id, row)
            .map_err(|e| force_pyerr(py, &e))
    }
    #[gen_stub(override_return_type(type_repr="numpy.typing.NDArray[numpy.float64]",imports=("numpy","numpy.typing")))]
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        snapshot(py, &self.inner.positions())
    }
    fn set_positions_(
        &mut self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[typing.Sequence[builtins.float]] | numpy.typing.NDArray[numpy.float64]",imports=("typing","numpy","numpy.typing")))]
        positions: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let rows = position_rows(py, positions)?;
        self.inner
            .set_positions_(&rows)
            .map_err(|e| force_pyerr(py, &e))
    }
    #[gen_stub(override_return_type(type_repr = "builtins.tuple[builtins.int, ...]"))]
    fn fixed_atoms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.inner.fixed_atoms().iter().map(|id| id.index()))
    }
    fn set_fixed_atoms_(
        &mut self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Iterable[builtins.int]",imports=("typing")))]
        atom_ids: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let items = atom_ids.try_iter().map_err(|e| {
            input_pyerr(
                py,
                ck::MolecularForceFieldErrorKind::InvalidFixedAtom,
                Some(e),
                None,
                None,
            )
        })?;
        let ids = items
            .map(|item| {
                item.map_err(|e| {
                    input_pyerr(
                        py,
                        ck::MolecularForceFieldErrorKind::InvalidFixedAtom,
                        Some(e),
                        None,
                        None,
                    )
                })
                .and_then(|item| {
                    parse_atom_id(
                        py,
                        &item,
                        ck::MolecularForceFieldErrorKind::InvalidFixedAtom,
                    )
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        self.inner
            .set_fixed_atoms_(&ids)
            .map_err(|e| force_pyerr(py, &e))
    }
    fn energy(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner.energy().map_err(|e| force_pyerr(py, &e))
    }
    #[gen_stub(override_return_type(type_repr="numpy.typing.NDArray[numpy.float64]",imports=("numpy","numpy.typing")))]
    fn gradient<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.inner
            .gradient()
            .map(|rows| snapshot(py, &rows))
            .map_err(|e| force_pyerr(py, &e))
    }
    /// Full energy derivative as an independent float64 (N, 3) array.
    /// Fixed-atom rows are not zeroed; the fixed set and positions are unchanged.
    /// Physical force is the negative of this gradient.
    #[gen_stub(override_return_type(type_repr="numpy.typing.NDArray[numpy.float64]",imports=("numpy","numpy.typing")))]
    fn gradient_unconstrained<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.inner
            .gradient_unconstrained()
            .map(|rows| snapshot(py, &rows))
            .map_err(|e| force_pyerr(py, &e))
    }
    fn energy_gradient(&self, py: Python<'_>) -> PyResult<ForceFieldEnergyGradient> {
        self.inner
            .energy_gradient()
            .map(|inner| ForceFieldEnergyGradient { inner })
            .map_err(|e| force_pyerr(py, &e))
    }
    #[pyo3(signature=(*,max_iterations=default_iterations(),force_tolerance=default_force_tolerance(),energy_tolerance=default_energy_tolerance()))]
    fn minimize_(
        &mut self,
        py: Python<'_>,
        max_iterations: u32,
        force_tolerance: f64,
        energy_tolerance: f64,
    ) -> PyResult<ForceFieldMinimizeOutcome> {
        let params =
            ck::ForceFieldMinimizeParams::new(max_iterations, force_tolerance, energy_tolerance);
        self.inner
            .minimize_with_params_(&params)
            .map(|inner| ForceFieldMinimizeOutcome { inner })
            .map_err(|e| force_pyerr(py, &e))
    }
    fn minimize_with_params_(
        &mut self,
        py: Python<'_>,
        params: &ForceFieldMinimizeParams,
    ) -> PyResult<ForceFieldMinimizeOutcome> {
        self.inner
            .minimize_with_params_(&params.inner)
            .map(|inner| ForceFieldMinimizeOutcome { inner })
            .map_err(|e| force_pyerr(py, &e))
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<MmffForceFieldParams>()?;
    module.add_class::<UffForceFieldParams>()?;
    module.add_class::<ForceFieldMinimizeParams>()?;
    module.add_class::<MolecularForceFieldErrorKind>()?;
    module.add_class::<ForceFieldEnergyGradient>()?;
    module.add_class::<ForceFieldMinimizeOutcome>()?;
    module.add_class::<MolecularForceField>()?;
    module.add("ForceFieldError", module.py().get_type::<ForceFieldError>())?;
    module.add(
        "MmffForceFieldError",
        module.py().get_type::<MmffForceFieldError>(),
    )?;
    module.add(
        "UffForceFieldError",
        module.py().get_type::<UffForceFieldError>(),
    )?;
    // Publish the registered payload getters as read-only Python properties.
    let descriptors = PyModule::from_code(
        module.py(),
        c"def getter(field):\n    return property(lambda self: getattr(self, '_' + field))\n",
        c"_persistent_error_properties",
        c"_persistent_error_properties",
    )?;
    for name in [
        "ForceFieldError",
        "MmffForceFieldError",
        "UffForceFieldError",
    ] {
        let class = module.getattr(name)?;
        let entry = ck::BINDING_CONTRACT
            .iter()
            .find(|entry| entry.python_name == name && entry.item == ck::BindingItem::Type)
            .expect("registered persistent exception");
        for field in ck::BINDING_CONTRACT_PROPERTIES
            .iter()
            .filter(|field| field.type_semantic_id == entry.semantic_id)
        {
            class.setattr(
                field.name,
                descriptors.getattr("getter")?.call1((field.name,))?,
            )?;
        }
    }
    Ok(())
}
