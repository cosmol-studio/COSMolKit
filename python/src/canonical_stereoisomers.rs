//! Source-defined stereoisomer options and lazy protocol projected through the facade.
use crate::drawing_binding::{Molecule, operation_pyerr};
use ::cosmolkit as ck;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::sync::{Arc, Mutex};

pyo3::create_exception!(
    cosmolkit,
    EnumerationError,
    PyValueError,
    "Stereoisomer enumeration failed for the graph or enumeration options."
);
pyo3::create_exception!(
    cosmolkit,
    EnumerationRunError,
    PyValueError,
    "Stereoisomer enumeration could not prepare or process the supplied molecule."
);

// Private binding failures retain their category, reason and boundary context.
// They are not the chemistry owner's shared-provider RandomSourcePoisoned.
#[derive(Debug)]
struct BindingLockFailure {
    context: &'static str,
    reason: String,
}
impl std::fmt::Display for BindingLockFailure {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "stereoisomer binding lock poisoned during {}: {}",
            self.context, self.reason
        )
    }
}
impl std::error::Error for BindingLockFailure {}

fn binding_lock_pyerr(
    py: Python<'_>,
    context: &'static str,
    reason: String,
    original: Option<PyErr>,
) -> PyErr {
    let source = BindingLockFailure { context, reason };
    let error = crate::canonical_values::annotate(
        py,
        PyRuntimeError::new_err(source.to_string()),
        "stereoisomers",
        "BindingLockPoisoned",
        &source,
    );
    error.set_cause(py, original);
    let attributes = error
        .value(py)
        .setattr("context", source.context)
        .and_then(|()| error.value(py).setattr("reason", source.reason));
    match attributes {
        Ok(()) => error,
        Err(attribute_error) => {
            attribute_error.set_cause(py, Some(error));
            attribute_error
        }
    }
}

type CallbackError = Arc<Mutex<Option<PyErr>>>;
fn take_callback_error(py: Python<'_>, slot: &CallbackError) -> PyResult<Option<PyErr>> {
    match slot.lock() {
        Ok(mut error) => Ok(error.take()),
        Err(cause) => {
            let reason = cause.to_string();
            // Read poisoned storage solely to retain a captured exception as
            // the failure's cause; never clear poison or resume enumeration.
            let original = cause.into_inner().take();
            Err(binding_lock_pyerr(
                py,
                "callback exception retrieval",
                reason,
                original,
            ))
        }
    }
}
fn save_callback_error(py: Python<'_>, slot: &CallbackError, error: PyErr) -> String {
    let message = error.to_string();
    match slot.lock() {
        Ok(mut saved) => {
            *saved = Some(error);
            message
        }
        Err(cause) => {
            let reason = cause.to_string();
            let failure = binding_lock_pyerr(py, "callback exception capture", reason, Some(error));
            let message = failure.to_string();
            // Retain the original Python exception identity through this typed
            // binding failure. The provider receives Err, never a random value.
            *cause.into_inner() = Some(failure);
            message
        }
    }
}

type BitsCallback =
    Box<dyn FnMut(usize) -> Result<num_bigint::BigUint, String> + Send + Sync + 'static>;
fn python_bits_callback(random: Py<PyAny>, call_method: bool) -> (BitsCallback, CallbackError) {
    let errors = Arc::new(Mutex::new(None));
    let slot = errors.clone();
    let callback: BitsCallback = Box::new(move |width| {
        Python::attach(|py| {
            let result = if call_method {
                random.bind(py).call_method1("getrandbits", (width,))
            } else {
                random.bind(py).call1((width,))
            };
            result
                .and_then(|value| value.extract::<num_bigint::BigUint>())
                .map_err(|error| save_callback_error(py, &slot, error))
        })
    });
    (callback, errors)
}
fn python_bits_source(
    random: Py<PyAny>,
    call_method: bool,
) -> (ck::StereoisomerRandomSource, CallbackError) {
    let (callback, errors) = python_bits_callback(random, call_method);
    (
        ck::StereoisomerRandomSource::from_random_bits(callback),
        errors,
    )
}

pub(crate) fn error_pyerr(py: Python<'_>, source: &ck::EnumerationError) -> PyErr {
    use ck::EnumerationError as E;
    let kind = match source {
        E::Valence(_) => "Valence",
        E::Rings(_) => "Rings",
        E::Stereo(_) => "Stereo",
        E::DoubleBond(_) => "DoubleBond",
        E::SmilesWrite(_) => "SmilesWrite",
        E::PotentialStereo(_) => "PotentialStereo",
        E::Hydrogens(_) => "Hydrogens",
        E::Embedding(_) => "Embedding",
        E::AtomProperty(_) => "AtomProperty",
        E::BondProperty(_) => "BondProperty",
        E::MoleculeProperty(_) => "MoleculeProperty",
        E::Topology(_) => "Topology",
        E::Coordinates(_) => "Coordinates",
        E::InvalidFlipperAtom { .. } => "InvalidFlipperAtom",
        E::InvalidFlipperBond { .. } => "InvalidFlipperBond",
        E::RandomBitsSource { .. } => "RandomBitsSource",
        E::RandomSourcePoisoned(_) => "RandomSourcePoisoned",
        E::EmbeddingCoordinates => "EmbeddingCoordinates",
    };
    let error = crate::canonical_values::annotate(
        py,
        EnumerationError::new_err(source.to_string()),
        "stereoisomers",
        kind,
        source,
    );
    let attributes = || -> PyResult<()> {
        match source {
            E::RandomBitsSource { bit_count, message } => {
                error.value(py).setattr("bit_count", *bit_count)?;
                error.value(py).setattr("message", message)?;
            }
            E::InvalidFlipperAtom { atom } => error.value(py).setattr("atom", atom.index())?,
            E::InvalidFlipperBond { bond } => error.value(py).setattr("bond", bond.index())?,
            _ => {}
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(error) => error,
    }
}
pub(crate) fn run_pyerr(py: Python<'_>, source: &ck::EnumerationRunError) -> PyErr {
    crate::canonical_values::annotate(
        py,
        EnumerationRunError::new_err(source.to_string()),
        "stereoisomers",
        "Enumeration",
        source,
    )
}
fn iterator_error(
    py: Python<'_>,
    source: ck::OperationError,
    slot: Option<&CallbackError>,
) -> PyErr {
    let original = match slot.map(|slot| take_callback_error(py, slot)).transpose() {
        Ok(error) => error.flatten(),
        Err(error) => Some(error),
    };
    let error = operation_pyerr(py, source);
    if let Some(original) = original {
        // The callback's typed Python exception is retained as the cause of the
        // source-specific enumeration error, below the runtime operation error.
        if let Some(run) = error.cause(py) {
            if let Some(enumeration) = run.cause(py) {
                enumeration.set_cause(py, Some(original));
            } else {
                run.set_cause(py, Some(original));
            }
        } else {
            error.set_cause(py, Some(original));
        }
    }
    error
}

/// Source-defined seed or shared random-bits provider.
///
/// Cloning a seed does not share an advanced generator. Cloning a provider
/// retains its identity and state without locking or drawing any bits.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "StereoisomerRandomSource", frozen, skip_from_py_object)]
pub(crate) struct StereoisomerRandomSource {
    inner: ck::StereoisomerRandomSource,
    errors: Option<CallbackError>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl StereoisomerRandomSource {
    /// Construct a reproducible stereoisomer random source from an integer seed.
    #[staticmethod]
    fn from_integer_seed(
        #[gen_stub(override_type(type_repr = "builtins.int", imports = ("builtins")))]
        value: num_bigint::BigInt,
    ) -> Self {
        Self {
            inner: ck::StereoisomerRandomSource::from_integer_seed(value),
            errors: None,
        }
    }
    /// Construct a stereoisomer random source from a callable accepting a bit count and returning a nonnegative integer.
    #[staticmethod]
    fn from_random_bits(value: Py<PyAny>) -> Self {
        let (inner, errors) = python_bits_source(value, false);
        Self {
            inner,
            errors: Some(errors),
        }
    }
    fn __repr__(&self) -> String {
        format!("StereoisomerRandomSource({:?})", self.inner)
    }
}

/// Writable configuration for stereoisomer enumeration and sampling.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "StereoisomerOptions", skip_from_py_object, dict, weakref)]
pub(crate) struct StereoisomerOptions {
    inner: ck::StereoisomerOptions,
    original_random_source: Option<Py<PyAny>>,
}
impl Default for StereoisomerOptions {
    fn default() -> Self {
        Self {
            inner: ck::StereoisomerOptions::default(),
            original_random_source: None,
        }
    }
}
impl StereoisomerOptions {
    fn enumeration_options(
        &self,
        py: Python<'_>,
    ) -> PyResult<(ck::StereoisomerOptions, Option<CallbackError>)> {
        let mut options = self.inner.clone();
        let Some(original) = &self.original_random_source else {
            return Ok((options, None));
        };
        if original
            .bind(py)
            .is_instance_of::<StereoisomerRandomSource>()
        {
            let source = original
                .bind(py)
                .extract::<PyRef<'_, StereoisomerRandomSource>>()?;
            options.set_random_source(Some(source.inner.clone()));
            return Ok((options, source.errors.clone()));
        }
        // Original d892 PyStereoisomerOptions::random_source (source projection):
        // let random_class = py.import("random")?.getattr("Random")?;
        // if rand.bind(py).is_instance(&random_class)? {
        //     Ok(Some(rand.clone_ref(py)))
        // } else {
        //     Ok(Some(random_class.call1((rand.bind(py),))?.unbind()))
        // Normalize at enumeration invocation, including exhaustive/no-center
        // branches. Construction, setters, getters and counting never normalize.
        let random_class = py.import("random")?.getattr("Random")?;
        let random = if original.bind(py).is_instance(&random_class)? {
            original.clone_ref(py)
        } else {
            random_class.call1((original.bind(py),))?.unbind()
        };
        let (source, errors) = python_bits_source(random, true);
        options.set_random_source(Some(source));
        Ok((options, Some(errors)))
    }
}
#[cosmolkit_macros::python_configuration(existing_setters)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl StereoisomerOptions {
    /// Configure stereoisomer enumeration and sampling; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature = (try_embedding=false, only_unassigned=true, max_isomers=1024, random_source=None, unique=true, only_stereo_groups=false))]
    fn new(
        try_embedding: bool,
        only_unassigned: bool,
        max_isomers: usize,
        random_source: Option<Py<PyAny>>,
        unique: bool,
        only_stereo_groups: bool,
    ) -> Self {
        Self {
            inner: ck::StereoisomerOptions::new(
                try_embedding,
                only_unassigned,
                max_isomers,
                None,
                unique,
                only_stereo_groups,
            ),
            original_random_source: random_source,
        }
    }
    /// Return a new StereoisomerOptions value with the standard enumeration defaults.
    #[staticmethod]
    #[pyo3(name = "default")]
    fn default_options() -> Self {
        <Self as Default>::default()
    }
    /// Whether candidate stereoisomers are embedded to reject infeasible configurations.
    #[getter]
    fn try_embedding(&self) -> bool {
        self.inner.try_embedding()
    }
    #[setter]
    fn set_try_embedding(&mut self, value: bool) {
        self.inner.set_try_embedding(value);
    }
    /// Whether enumeration changes only currently unspecified stereocenters.
    #[getter]
    fn only_unassigned(&self) -> bool {
        self.inner.only_unassigned()
    }
    #[setter]
    fn set_only_unassigned(&mut self, value: bool) {
        self.inner.set_only_unassigned(value);
    }
    /// Whether enumeration is restricted to enhanced stereo groups.
    #[getter]
    fn only_stereo_groups(&self) -> bool {
        self.inner.only_stereo_groups()
    }
    #[setter]
    fn set_only_stereo_groups(&mut self, value: bool) {
        self.inner.set_only_stereo_groups(value);
    }
    /// Maximum number of stereoisomers to produce; zero requests exhaustive enumeration.
    #[getter]
    fn max_isomers(&self) -> usize {
        self.inner.max_isomers()
    }
    #[setter]
    fn set_max_isomers(&mut self, value: usize) {
        self.inner.set_max_isomers(value);
    }
    /// Seeded or callback-based random source used when stereoisomer sampling is required.
    #[getter]
    fn random_source(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.original_random_source
            .as_ref()
            .map(|value| value.clone_ref(py))
    }
    #[setter]
    fn set_random_source(&mut self, value: Option<Py<PyAny>>) {
        self.original_random_source = value;
    }
    /// Whether equivalent stereoisomers are deduplicated.
    #[getter]
    fn unique(&self) -> bool {
        self.inner.unique()
    }
    #[setter]
    fn set_unique(&mut self, value: bool) {
        self.inner.set_unique(value);
    }
    fn __repr__(&self) -> String {
        format!(
            "StereoisomerOptions(try_embedding={}, only_unassigned={}, max_isomers={}, random_source={}, unique={}, only_stereo_groups={})",
            self.inner.try_embedding(),
            self.inner.only_unassigned(),
            self.inner.max_isomers(),
            if self.original_random_source.is_some() {
                "..."
            } else {
                "None"
            },
            self.inner.unique(),
            self.inner.only_stereo_groups()
        )
    }
}

/// Lazy iterator over source-ordered stereoisomers.
///
/// Construction performs the source-defined preprocessing and candidate
/// discovery. Configuration application, uniqueness, optional embedding, and
/// their errors are deferred until ``next()`` requests an output.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(name = "StereoisomerIterator", skip_from_py_object)]
pub(crate) struct StereoisomerIterator {
    inner: Mutex<ck::StereoisomerIterator>,
    errors: Option<CallbackError>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl StereoisomerIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    #[gen_stub(override_return_type(type_repr = "Molecule"))]
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<Molecule>> {
        self.next(py)
    }
    /// Return the next stereoisomer, or None when enumeration is exhausted.
    fn next(&mut self, py: Python<'_>) -> PyResult<Option<Molecule>> {
        self.inner
            .get_mut()
            .map_err(|cause| binding_lock_pyerr(py, "iterator.next", cause.to_string(), None))?
            .next()
            .transpose()
            .map(|row| row.map(Molecule::from_inner))
            .map_err(|error| iterator_error(py, error, self.errors.as_ref()))
    }
    /// Number of stereoisomers already yielded by this iterator.
    #[getter]
    fn yielded_count(&self, py: Python<'_>) -> PyResult<usize> {
        self.inner
            .lock()
            .map(|inner| inner.yielded_count())
            .map_err(|cause| {
                binding_lock_pyerr(py, "iterator.yielded_count", cause.to_string(), None)
            })
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "StereoisomerIterator(yielded_count={})",
            self.yielded_count(py)?
        ))
    }
}

pub(crate) fn enumerate_stereoisomers(
    molecule: &Molecule,
    py: Python<'_>,
) -> PyResult<StereoisomerIterator> {
    let inner = molecule
        .inner
        .enumerate_stereoisomers()
        .map_err(|error| operation_pyerr(py, error))?;
    Ok(StereoisomerIterator {
        inner: Mutex::new(inner),
        errors: None,
    })
}
pub(crate) fn enumerate_stereoisomers_with_options(
    molecule: &Molecule,
    py: Python<'_>,
    options: &StereoisomerOptions,
) -> PyResult<StereoisomerIterator> {
    let (options, errors) = options.enumeration_options(py)?;
    let inner = molecule
        .inner
        .enumerate_stereoisomers_with_options(&options)
        .map_err(|error| iterator_error(py, error, errors.as_ref()))?;
    Ok(StereoisomerIterator {
        inner: Mutex::new(inner),
        errors,
    })
}
pub(crate) fn enumerate_stereoisomers_with_random_bits(
    molecule: &Molecule,
    py: Python<'_>,
    options: &StereoisomerOptions,
    callback: Py<PyAny>,
) -> PyResult<StereoisomerIterator> {
    let (callback, errors) = python_bits_callback(callback, false);
    let inner = molecule
        .inner
        .enumerate_stereoisomers_with_random_bits(&options.inner, callback)
        .map_err(|error| iterator_error(py, error, Some(&errors)))?;
    Ok(StereoisomerIterator {
        inner: Mutex::new(inner),
        errors: Some(errors),
    })
}
pub(crate) fn stereoisomer_count(
    molecule: &Molecule,
    py: Python<'_>,
) -> PyResult<num_bigint::BigUint> {
    molecule
        .inner
        .stereoisomer_count()
        .map_err(|error| error_pyerr(py, &error))
}
pub(crate) fn stereoisomer_count_with_options(
    molecule: &Molecule,
    py: Python<'_>,
    options: &StereoisomerOptions,
) -> PyResult<num_bigint::BigUint> {
    // Counting does not normalize Python seed objects or initialize generators.
    molecule
        .inner
        .stereoisomer_count_with_options(&options.inner)
        .map_err(|error| error_pyerr(py, &error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<StereoisomerOptions>()?;
    module.add_class::<StereoisomerRandomSource>()?;
    module.add_class::<StereoisomerIterator>()?;
    module.add(
        "EnumerationError",
        module.py().get_type::<EnumerationError>(),
    )?;
    module.add(
        "EnumerationRunError",
        module.py().get_type::<EnumerationRunError>(),
    )?;
    Ok(())
}
