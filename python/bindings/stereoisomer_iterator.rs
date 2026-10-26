//! Lazy native iteration retains the exact Python random-provider owner.
use crate::native::Molecule;
use crate::operation_errors::OperationError;
use crate::stereoisomer_options::{
    RandomBitsCallback, StereoisomerOptions, StereoisomerRandomSource, with_random_bits_errors,
};
use pyo3::prelude::*;
use std::sync::{Mutex, MutexGuard};

#[pyclass(module = "cosmolkit")]
/// Lazy iterator producing stereoisomer molecules while retaining its random source.
pub struct StereoisomerIterator {
    inner: Mutex<cosmolkit::StereoisomerIterator>,
    source_owner: Option<Py<StereoisomerRandomSource>>,
    callback_owner: Option<Py<StereoisomerRandomSource>>,
}
impl StereoisomerIterator {
    fn borrow(&self) -> PyResult<MutexGuard<'_, cosmolkit::StereoisomerIterator>> {
        // A callback can reenter its iterator on the same Python thread. Never
        // wait on that active native loan; this is a binding borrow error.
        self.inner.try_lock().map_err(|error| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "stereoisomer iterator is unavailable: {error}"
            ))
        })
    }
}
#[pymethods]
impl StereoisomerIterator {
    /// Return the next stereoisomer, or None when enumeration is exhausted.
    fn next(&self, py: Python<'_>) -> PyResult<Option<Molecule>> {
        let mut inner = self.borrow()?;
        let native = &mut *inner;
        with_random_bits_errors(|| py.detach(|| native.next()))?
            .transpose()
            .map(|value| value.map(|inner| Molecule { inner }))
            .map_err(OperationError::from_rust)
    }
    fn __next__(&self, py: Python<'_>) -> PyResult<Option<Molecule>> {
        self.next(py)
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    #[getter]
    /// Number of stereoisomers yielded so far.
    fn yielded_count(&self) -> PyResult<usize> {
        Ok(self.borrow()?.yielded_count())
    }
    fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        if let Some(value) = &self.source_owner {
            visit.call(value)?;
        }
        if let Some(value) = &self.callback_owner {
            visit.call(value)?;
        }
        Ok(())
    }
    fn __clear__(&mut self) {
        self.source_owner = None;
        self.callback_owner = None;
    }
}

pub(crate) fn enumerate(
    py: Python<'_>,
    molecule: &Molecule,
    options: Option<Py<StereoisomerOptions>>,
) -> PyResult<StereoisomerIterator> {
    let (inner, source_owner) = match options {
        Some(options) => {
            let (options, owner) = {
                let options = options.borrow(py);
                (options.inner.clone(), options.random_source(py))
            };
            (
                with_random_bits_errors(|| {
                    py.detach(|| {
                        molecule
                            .inner
                            .enumerate_stereoisomers_with_options(&options)
                    })
                })?
                .map_err(OperationError::from_rust)?,
                owner,
            )
        }
        None => (
            with_random_bits_errors(|| py.detach(|| molecule.inner.enumerate_stereoisomers()))?
                .map_err(OperationError::from_rust)?,
            None,
        ),
    };
    Ok(StereoisomerIterator {
        inner: Mutex::new(inner),
        source_owner,
        callback_owner: None,
    })
}
pub(crate) fn with_random_bits(
    py: Python<'_>,
    molecule: &Molecule,
    options: Py<StereoisomerOptions>,
    callback: RandomBitsCallback,
) -> PyResult<StereoisomerIterator> {
    let (options, source_owner) = {
        let options = options.borrow(py);
        (options.inner.clone(), options.random_source(py))
    };
    let callback_owner = StereoisomerRandomSource::from_random_bits(py, callback)?;
    let provider = StereoisomerRandomSource::provider(callback_owner.bind(py))?;
    let inner = with_random_bits_errors(|| {
        py.detach(|| {
            molecule
                .inner
                .enumerate_stereoisomers_with_random_bits(&options, provider)
        })
    })?
    .map_err(OperationError::from_rust)?;
    Ok(StereoisomerIterator {
        inner: Mutex::new(inner),
        source_owner,
        callback_owner: Some(callback_owner),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit as ck;
    use pyo3::ffi::c_str;
    use pyo3::types::{PyDict, PyModule};

    #[test]
    fn lazy_stereoisomers_preserve_native_order_provider_owner_gc_and_original_python_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native=ck::Molecule::from_smiles("FC(Cl)C(Br)I").unwrap();
            let before=native.to_binary().unwrap();
            let expected=native.enumerate_stereoisomers().unwrap().map(|value|crate::property_values::decode_source_text(py,&value.unwrap().to_smiles().unwrap())).collect::<PyResult<Vec<_>>>()?;
            let module=PyModule::new(py,"stereoisomer_iterator_test")?;
            module.add_class::<Molecule>()?;
            module.add_class::<StereoisomerOptions>()?;
            module.add_class::<StereoisomerRandomSource>()?;
            let locals=PyDict::new(py);
            locals.set_item("api",module)?;
            locals.set_item("molecule",Py::new(py,Molecule{inner:native})?)?;
            locals.set_item("before",before)?;
            locals.set_item("expected",expected)?;
            py.run(c_str!(r#"
import gc
import weakref
options = api.StereoisomerOptions()
for stream in [molecule.enumerate_stereoisomers(), molecule.enumerate_stereoisomers(options),
               molecule.enumerate_stereoisomers(max_isomers=1024),
               molecule.enumerate_stereoisomers_with_options(options)]:
    assert iter(stream) is stream and stream.yielded_count == 0
    first = stream.next()
    assert first is not None and stream.yielded_count == 1
    actual = [first.to_smiles()] + [value.to_smiles() for value in stream]
    assert actual == expected
    assert stream.next() is None and stream.next() is None
    assert stream.yielded_count == len(expected)
assert bytes(molecule.to_binary()) == bytes(before)
try: molecule.enumerate_stereoisomers(options, max_isomers=1)
except TypeError: pass
else: raise AssertionError('mixed iterator configuration accepted')
class Provider:
    def __init__(self): self.calls = []; self.iterator = None
    def __call__(self, bits): self.calls.append(bits); return 0
provider = Provider()
reference = weakref.ref(provider)
source = api.StereoisomerRandomSource.from_random_bits(provider)
options = api.StereoisomerOptions(max_isomers=1, random_source=source)
stream = molecule.enumerate_stereoisomers(options)
del source, options, provider, molecule
gc.collect()
assert reference() is not None
result = stream.next()
assert result is not None and reference().calls
assert stream.next() is None and stream.yielded_count == 1
del stream
gc.collect()
assert reference() is None
molecule = api.Molecule.from_smiles('FC(Cl)C(Br)I')
provider = Provider()
reference = weakref.ref(provider)
source = api.StereoisomerRandomSource.from_random_bits(provider)
options = api.StereoisomerOptions(max_isomers=1, random_source=source)
stream = molecule.enumerate_stereoisomers(options)
provider.iterator = stream
del stream, source, options, provider
gc.collect()
assert reference() is None  # native closure must not hide a strong Python edge
provider = Provider()
reference = weakref.ref(provider)
options = api.StereoisomerOptions(max_isomers=1)
stream = molecule.enumerate_stereoisomers_with_random_bits(options, provider)
del provider, options, molecule
gc.collect()
assert reference() is not None
assert stream.next() is not None and reference().calls
assert stream.next() is None
del stream
gc.collect()
assert reference() is None
molecule = api.Molecule.from_smiles('FC(Cl)C(Br)I')
marker = LookupError('original random-bits exception')
def bad(bits): raise marker
stream = None
try:
    stream = molecule.enumerate_stereoisomers_with_random_bits(api.StereoisomerOptions(max_isomers=1), bad)
    stream.next()
except LookupError as error: assert error is marker
else: raise AssertionError('random-bits exception swallowed')
if stream is not None: assert stream.next() is None
assert molecule.enumerate_stereoisomers().next() is not None
"#),Some(&locals),None)?;
            Ok(())
        }).unwrap();
    }
}
