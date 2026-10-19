//! Exact native enumeration options and random-source representation.
use crate::callbacks::{CallbackErrors, typed_callback};
use cosmolkit as ck;
use num_bigint::{BigInt, BigUint};
use pyo3::prelude::*;
use pyo3::types::{PyWeakrefMethods, PyWeakrefReference};
use std::cell::RefCell;
use std::sync::Arc;

typed_callback!(
    RandomBitsCallback,
    [pyo3::type_hint_identifier!("builtins", "int")],
    pyo3::type_hint_identifier!("builtins", "int")
);

thread_local! {
    // The native random source invokes its provider synchronously on next().
    // A stack isolates nested calls and threads without changing source identity.
    static RANDOM_ERRORS: RefCell<Vec<Arc<CallbackErrors>>> = const { RefCell::new(Vec::new()) };
}

struct RandomBitsCall;
impl Drop for RandomBitsCall {
    fn drop(&mut self) {
        RANDOM_ERRORS.with(|stack| {
            stack.borrow_mut().pop();
        });
    }
}

pub(crate) fn with_random_bits_errors<T>(run: impl FnOnce() -> T) -> PyResult<T> {
    let errors = Arc::new(CallbackErrors::new());
    RANDOM_ERRORS.with(|stack| stack.borrow_mut().push(Arc::clone(&errors)));
    let scope = RandomBitsCall;
    let result = run();
    drop(scope);
    match errors.take() {
        Some(error) => Err(error),
        None => Ok(result),
    }
}

fn draw_random_bits(run: impl FnOnce(Python<'_>) -> PyResult<BigUint>) -> Result<BigUint, String> {
    let errors = RANDOM_ERRORS
        .with(|stack| stack.borrow().last().cloned())
        .ok_or_else(|| "Python random-bits callback requires its binding invocation".to_owned())?;
    Python::attach(run).map_err(|error| {
        errors.record(error);
        // Rust's callback contract carries a String; the active binding
        // invocation restores the original Python error, not this sentinel.
        "Python random-bits callback failed".to_owned()
    })
}

#[pyclass(module = "cosmolkit", weakref, eq)]
pub struct StereoisomerRandomSource {
    pub(crate) inner: ck::StereoisomerRandomSource,
    // This is the only strong Python callback owner. Native source clones
    // share their provider identity but hold only a weak reference back here.
    callback: Option<RandomBitsCallback>,
}

impl PartialEq for StereoisomerRandomSource {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

#[pymethods]
impl StereoisomerRandomSource {
    #[staticmethod]
    fn from_integer_seed(value: BigInt) -> Self {
        Self {
            inner: ck::StereoisomerRandomSource::from_integer_seed(value),
            callback: None,
        }
    }
    #[staticmethod]
    pub(crate) fn from_random_bits(
        py: Python<'_>,
        value: RandomBitsCallback,
    ) -> PyResult<Py<Self>> {
        // Allocate the Python owner before creating its weak native provider.
        // The temporary seed cannot be observed: the source is returned only
        // after the actual provider has been installed, without drawing bits.
        let owner = Bound::new(
            py,
            Self {
                inner: ck::StereoisomerRandomSource::from_integer_seed(BigInt::from(0)),
                callback: Some(value),
            },
        )?;
        let provider = Self::provider(&owner)?;
        owner.borrow_mut().inner = ck::StereoisomerRandomSource::from_random_bits(provider);
        Ok(owner.unbind())
    }
    fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        if let Some(callback) = &self.callback {
            callback.traverse(&visit)?;
        }
        Ok(())
    }
    fn __clear__(&mut self) {
        self.callback = None;
    }
    fn __repr__(&self) -> String {
        format!("StereoisomerRandomSource({:?})", self.inner)
    }
}

impl StereoisomerRandomSource {
    // Reuse the sole checked callback/weak-owner path for the native explicit
    // random-bits overload as well as the native source selector factory.
    pub(crate) fn provider(
        owner: &Bound<'_, Self>,
    ) -> PyResult<Box<dyn FnMut(usize) -> Result<BigUint, String> + Send + Sync + 'static>> {
        let weak = PyWeakrefReference::new(owner.as_any())?.unbind();
        Ok(Box::new(move |bits| {
            draw_random_bits(|py| {
                let owner = weak.bind(py).upgrade_as::<Self>()?.ok_or_else(|| {
                    pyo3::exceptions::PyReferenceError::new_err("random source owner has expired")
                })?;
                let callback = owner
                    .try_borrow()?
                    .callback
                    .as_ref()
                    .ok_or_else(|| {
                        pyo3::exceptions::PyReferenceError::new_err(
                            "random source callback has been cleared",
                        )
                    })?
                    .clone_ref(py);
                callback.bind(py).call1((bits,))?.extract::<BigUint>()
            })
        }))
    }
}

#[pyclass(module = "cosmolkit")]
pub struct StereoisomerOptions {
    pub(crate) inner: ck::StereoisomerOptions,
    // Keep the exact source alive rather than an untraversable native clone.
    source_owner: Option<Py<StereoisomerRandomSource>>,
}

macro_rules! stereoisomer_configuration {
    ($($field:ident: $ty:ty),+ $(,)?) => {
        paste::paste! {
            #[pymethods]
            impl StereoisomerOptions {
                #[new]
                #[pyo3(signature = (*, $($field=ck::StereoisomerOptions::default().$field,)+ random_source=None))]
                pub(crate) fn new(py: Python<'_>, $($field: $ty,)+ random_source: Option<Py<StereoisomerRandomSource>>) -> Self {
                    let mut inner = ck::StereoisomerOptions::default();
                    $(inner.$field = $field;)+
                    inner.set_random_source(random_source.as_ref().map(|source| source.borrow(py).inner.clone()));
                    Self { inner, source_owner: random_source }
                }
                #[staticmethod]
                fn default() -> Self {
                    Self { inner: ck::StereoisomerOptions::default(), source_owner: None }
                }
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.[<set_ $field>](value); })+
    #[getter]
    pub(crate) fn random_source(&self, py: Python<'_>) -> Option<Py<StereoisomerRandomSource>> {
        self.source_owner
            .as_ref()
            .map(|source| source.clone_ref(py))
    }
    #[setter]
    fn set_random_source(&mut self, py: Python<'_>, value: Option<Py<StereoisomerRandomSource>>) {
        self.inner
            .set_random_source(value.as_ref().map(|source| source.borrow(py).inner.clone()));
        self.source_owner = value;
    }
    fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        visit.call(&self.source_owner)
    }
    fn __clear__(&mut self) {
        self.inner.set_random_source(None);
        self.source_owner = None;
    }

                fn __repr__(&self) -> String { format!("{:?}", self.inner) }
            }
        }
    };
}
ck::project_configuration!(ck::StereoisomerOptions, using stereoisomer_configuration);

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::PyDict;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    #[test]
    fn random_provider_owners_are_retained_and_unreachable_cycles_are_collected() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            use pyo3::ffi::c_str;
            let module = pyo3::types::PyModule::new(py, "random_gc_test")?;
            module.add_class::<StereoisomerRandomSource>()?;
            module.add_class::<StereoisomerOptions>()?;
            let scope = PyDict::new(py);
            scope.set_item("ck", module)?;
            py.run(
                c_str!(
                    r#"
import gc
import weakref
class Provider:
    def __init__(self):
        self.calls = 0
    def __call__(self, bits):
        self.calls += 1
        return 0

def source_cycle():
    provider = Provider()
    source = ck.StereoisomerRandomSource.from_random_bits(provider)
    provider.source = source
    return weakref.ref(provider)

def options_cycle():
    provider = Provider()
    source = ck.StereoisomerRandomSource.from_random_bits(provider)
    options = ck.StereoisomerOptions(random_source=source)
    assert options.random_source is source
    provider.options = options
    return weakref.ref(provider)

for create_cycle in (source_cycle, options_cycle):
    reference = create_cycle()
    gc.collect()
    assert reference() is None

provider = Provider()
reference = weakref.ref(provider)
source = ck.StereoisomerRandomSource.from_random_bits(provider)
assert source != object()
options = ck.StereoisomerOptions(max_isomers=1, random_source=source)
assert options.random_source is source
del provider, source
gc.collect()
assert reference() is not None
"#
                ),
                Some(&scope),
                None,
            )?;
            // Only Options retains the Python provider owner. Cloning the
            // native configuration must not clone/reset its shared provider.
            let options = scope
                .get_item("options")?
                .unwrap()
                .extract::<Py<StereoisomerOptions>>()?;
            let native = options.borrow(py).inner.clone();
            let molecule = ck::Molecule::from_smiles("FC(Cl)C(Br)I").unwrap();
            assert!(
                with_random_bits_errors(|| molecule
                    .enumerate_stereoisomers_with_options(&native)
                    .unwrap()
                    .next()
                    .unwrap())?
                .is_ok()
            );
            py.run(
                c_str!(
                    r#"
assert reference().calls > 0
retained = options.random_source
options.random_source = None
gc.collect()
assert reference() is not None
options.random_source = retained
assert options.random_source is retained
del retained
options.random_source = None
gc.collect()
assert reference() is None
"#
                ),
                Some(&scope),
                None,
            )?;
            // A bare native provider clone is not a Python lifetime owner.
            // Expiry becomes the original structured language error, never
            // a numerical fallback or a dangling callback access.
            let error = with_random_bits_errors(|| {
                molecule
                    .enumerate_stereoisomers_with_options(&native)
                    .unwrap()
                    .next()
                    .unwrap()
            })
            .unwrap_err();
            assert!(error.is_instance_of::<pyo3::exceptions::PyReferenceError>(py));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn arbitrary_integer_sources_and_provider_errors_preserve_actual_values_and_identity() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            use pyo3::ffi::c_str;
            use pyo3::types::PyModule;
            let probe = PyModule::from_code(
                py,
                c_str!(
                    r#"
huge = (1 << 20000) + 17
seen = []
fail = False
negative = False
marker = ValueError('original random provider failure')
def wide(bits):
    return huge
def provider(bits):
    seen.append(bits)
    if fail:
        raise marker
    return -1 if negative else 0
"#
                ),
                c_str!("random_source_probe.py"),
                c_str!("random_source_probe"),
            )?;
            let class = py.get_type::<StereoisomerRandomSource>();
            let native_positive: BigInt = (BigInt::from(1_u8) << 20000_usize) + 17_u8;
            for expected in [&native_positive, &-native_positive.clone()] {
                let input = expected.into_pyobject(py)?;
                let source = class
                    .call_method1("from_integer_seed", (input,))?
                    .extract::<Py<StereoisomerRandomSource>>()?;
                assert_eq!(
                    source.borrow(py).inner,
                    ck::StereoisomerRandomSource::from_integer_seed(expected.clone())
                );
            }
            let wide: RandomBitsCallback = probe.getattr("wide")?.extract()?;
            let expected_unsigned = (BigUint::from(1_u8) << 20000_usize) + 17_u8;
            assert_eq!(
                with_random_bits_errors(|| draw_random_bits(|py| wide
                    .bind(py)
                    .call1((20001,))?
                    .extract::<BigUint>()))?
                .unwrap(),
                expected_unsigned
            );
            assert_eq!(probe.getattr("huge")?.extract::<BigInt>()?, native_positive);
            let source = class
                .call_method1("from_random_bits", (probe.getattr("provider")?,))?
                .extract::<Py<StereoisomerRandomSource>>()?;
            let same_source = source.borrow(py).inner.clone();
            let options = ck::StereoisomerOptions {
                max_isomers: 1,
                random_source: Some(same_source.clone()),
                ..ck::StereoisomerOptions::default()
            };
            let molecule = ck::Molecule::from_smiles("FC(Cl)C(Br)I").unwrap();
            let run = || {
                let mut iterator = molecule
                    .enumerate_stereoisomers_with_options(&options)
                    .unwrap();
                iterator.next().unwrap()
            };
            assert!(with_random_bits_errors(run)?.is_ok());
            assert!(probe.getattr("seen")?.len()? > 0);
            probe.setattr("fail", true)?;
            let error = with_random_bits_errors(run).unwrap_err();
            assert!(error.value(py).is(&probe.getattr("marker")?));
            probe.setattr("fail", false)?;
            probe.setattr("negative", true)?;
            assert!(
                with_random_bits_errors(run)
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyOverflowError>(py)
            );
            probe.setattr("negative", false)?;
            assert!(with_random_bits_errors(run)?.is_ok());
            assert_eq!(source.borrow(py).inner, same_source);
            // A nested invocation's failure belongs only to its own scope.
            probe.setattr("fail", true)?;
            let nested = with_random_bits_errors(|| with_random_bits_errors(run))?;
            assert!(nested.unwrap_err().value(py).is(&probe.getattr("marker")?));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn full_options_preserve_native_configuration_and_existing_source_state() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let class = py.get_type::<StereoisomerOptions>();
            let options = class.call0()?.extract::<Py<StereoisomerOptions>>()?;
            assert_eq!(options.borrow(py).inner, ck::StereoisomerOptions::default());
            let default = class
                .call_method0("default")?
                .extract::<Py<StereoisomerOptions>>()?;
            assert_eq!(default.borrow(py).inner, options.borrow(py).inner);
            let molecule = ck::Molecule::from_smiles("F[C@](Cl)(Br)I").unwrap();
            let original_count = molecule
                .stereoisomer_count_with_options(&options.borrow(py).inner)
                .unwrap();
            options.bind(py).setattr("only_unassigned", false)?;
            let supplied = molecule
                .stereoisomer_count_with_options(&options.borrow(py).inner)
                .unwrap();
            assert_ne!(original_count, supplied);
            let mut native = ck::StereoisomerOptions::default();
            native.set_only_unassigned(false);
            assert_eq!(
                supplied,
                molecule.stereoisomer_count_with_options(&native).unwrap()
            );

            let calls = Arc::new(AtomicUsize::new(0));
            let observed = Arc::clone(&calls);
            let source = ck::StereoisomerRandomSource::from_random_bits(Box::new(move |_| {
                observed.fetch_add(1, Ordering::SeqCst);
                Ok("0".parse().unwrap())
            }));
            let expected_source = source.clone();
            let source = Py::new(
                py,
                StereoisomerRandomSource {
                    inner: source,
                    callback: None,
                },
            )?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("try_embedding", false)?;
            kwargs.set_item("only_unassigned", false)?;
            kwargs.set_item("max_isomers", 1)?;
            kwargs.set_item("random_source", &source)?;
            kwargs.set_item("unique", false)?;
            kwargs.set_item("only_stereo_groups", false)?;
            let options = class
                .call((), Some(&kwargs))?
                .extract::<Py<StereoisomerOptions>>()?;
            assert_eq!(
                options.borrow(py).inner,
                ck::StereoisomerOptions::new(
                    false,
                    false,
                    1,
                    Some(expected_source.clone()),
                    false,
                    false
                )
            );
            let retained = options
                .bind(py)
                .getattr("random_source")?
                .extract::<Py<StereoisomerRandomSource>>()?;
            assert_eq!(retained.borrow(py).inner, expected_source);
            options.bind(py).setattr("random_source", py.None())?;
            assert!(options.borrow(py).inner.random_source().is_none());
            options.bind(py).setattr("random_source", &retained)?;
            assert_eq!(
                options.borrow(py).inner.random_source(),
                Some(&expected_source)
            );
            drop(expected_source);
            drop(source);
            drop(retained);
            drop(kwargs);
            let molecule = ck::Molecule::from_smiles("FC(Cl)C(Br)I").unwrap();
            let mut iterator = molecule
                .enumerate_stereoisomers_with_options(&options.borrow(py).inner)
                .unwrap();
            assert!(iterator.next().unwrap().is_ok());
            assert!(calls.load(Ordering::SeqCst) > 0);
            Ok(())
        })
        .unwrap();
    }
}
