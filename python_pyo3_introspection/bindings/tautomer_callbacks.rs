//! Synchronous native callback loans. Only explicit snapshots retain native data.
use crate::atom_bond_values::{Atom, Bond};
use crate::attachment_values::AtomMetadata;
use crate::callbacks::{CallbackErrors, typed_callback};
use crate::leaf_values::TautomerEnumerationStatus;
use crate::metadata_values::MoleculeProperties;
use crate::property_values::decode_source_text;
use crate::tautomer_errors::TautomerRunError;
use crate::tautomer_values::TautomerScore;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::{PyReferenceError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use std::collections::BTreeSet;
use std::ffi::c_void;
use std::ptr::NonNull;
use std::sync::{Arc, Mutex};
use std::thread::{self, ThreadId};

struct NativeLoan<'source, 'progress> {
    source: ck::TautomerMoleculeView<'source>,
    progress: Option<ck::TautomerProgress<'progress>>,
}

struct NativePointer(NonNull<c_void>);
// SAFETY: moving this private pointer handle never moves or accesses the
// pointee. Lease::with_native checks the immutable origin ThreadId before any
// dereference; no other function dereferences it. Mutex serializes access and
// invalidation. The origin scope clears the pointer before dropping its native
// stack loan. Neither Send for the pointee nor Sync for this wrapper is claimed.
unsafe impl Send for NativePointer {}

/// All aliases of a single callback use the same dynamic exclusive borrow.
#[derive(Clone)]
struct Lease {
    origin: ThreadId,
    slot: Arc<Mutex<Option<NativePointer>>>,
}

impl Lease {
    fn with_native<R>(
        &self,
        run: impl for<'source, 'progress> FnOnce(&mut NativeLoan<'source, 'progress>) -> R,
    ) -> PyResult<R> {
        let mut slot = self.slot.try_lock().map_err(|_| {
            PyRuntimeError::new_err("tautomer callback view is borrowed or unavailable")
        })?;
        let pointer = slot.as_mut().ok_or_else(|| {
            PyReferenceError::new_err(
                "tautomer callback view has expired; retain snapshot() instead",
            )
        })?;
        // A foreign thread may inspect the guarded expired/active flag, never
        // the pointee. Expiry takes precedence over affinity on every thread.
        if thread::current().id() != self.origin {
            return Err(PyRuntimeError::new_err(
                "live tautomer callback views can only be used on their callback thread; retain snapshot() instead",
            ));
        }
        // SAFETY: scope() is the only constructor. The pointer addresses its
        // stack-local NativeLoan, which is not moved until the guard invalidates
        // this slot. The origin-thread check precedes this sole dereference;
        // a foreign thread can retain/drop a handle but never access the loan.
        // MutexGuard excludes every other native access, including reentry
        // through a candidate alias. The HRTB closure cannot return a
        // reference tied to either erased native lifetime. No Python call or
        // conversion runs inside these closures. No 'static reference is made.
        let loan = unsafe { pointer.0.cast::<NativeLoan<'_, '_>>().as_mut() };
        Ok(run(loan))
    }
}

struct Expire(Lease);
impl Drop for Expire {
    fn drop(&mut self) {
        // with_native cannot export its MutexGuard. During unwinding it is dropped
        // before this outer guard. Invalidate before NativeLoan is destroyed.
        self.0
            .slot
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .take();
    }
}

fn scope<R>(
    source: ck::TautomerMoleculeView<'_>,
    progress: Option<ck::TautomerProgress<'_>>,
    run: impl FnOnce(Lease) -> PyResult<R>,
) -> PyResult<R> {
    let mut native = NativeLoan { source, progress };
    let lease = Lease {
        origin: thread::current().id(),
        slot: Arc::new(Mutex::new(Some(NativePointer(
            NonNull::from(&mut native).cast(),
        )))),
    };
    let guard = Expire(lease.clone());
    let result = run(lease);
    drop(guard);
    result
}

#[derive(Clone)]
enum ProgressStorage {
    Live(Lease),
    Owned(Arc<Mutex<ck::TautomerProgress<'static>>>),
}
impl ProgressStorage {
    fn with_native<R>(
        &self,
        run: impl for<'a> FnOnce(&mut ck::TautomerProgress<'a>) -> R,
    ) -> PyResult<R> {
        match self {
            Self::Live(lease) => lease.with_native(|loan| {
                run(loan.progress.as_mut().expect("progress loan constructor"))
            }),
            Self::Owned(owner) => {
                let mut native = owner.try_lock().map_err(|_| {
                    PyRuntimeError::new_err("tautomer progress is borrowed or unavailable")
                })?;
                Ok(run(&mut native))
            }
        }
    }
}

enum ViewStorage {
    Source(Lease),
    Entry {
        progress: ProgressStorage,
        index: usize,
    },
    Owned(Mutex<ck::TautomerMoleculeView<'static>>),
}

#[pyclass(module = "fixed_projection._native")]
pub struct TautomerMoleculeView {
    storage: ViewStorage,
}
impl TautomerMoleculeView {
    fn with_native<R>(
        &self,
        run: impl for<'a> FnOnce(&mut ck::TautomerMoleculeView<'a>) -> R,
    ) -> PyResult<R> {
        match &self.storage {
            ViewStorage::Source(lease) => lease.with_native(|loan| run(&mut loan.source)),
            ViewStorage::Entry { progress, index } => {
                let result = progress.with_native(|native| {
                    native.entry(*index).map(|(_, mut value)| run(&mut value))
                })?;
                result
                    .ok_or_else(|| PyReferenceError::new_err("tautomer candidate no longer exists"))
            }
            ViewStorage::Owned(owner) => {
                let mut value = owner.try_lock().map_err(|_| {
                    PyRuntimeError::new_err("tautomer snapshot is borrowed or unavailable")
                })?;
                Ok(run(&mut value))
            }
        }
    }
}

#[pymethods]
impl TautomerMoleculeView {
    /// Explicitly own source blocks, or share the actual candidate snapshot.
    fn snapshot(&self) -> PyResult<Self> {
        let value = self.with_native(|native| native.snapshot())?;
        Ok(Self {
            storage: ViewStorage::Owned(Mutex::new(value)),
        })
    }
    fn to_owned(&self) -> PyResult<Self> {
        let value = self.with_native(|native| native.to_owned())?;
        Ok(Self {
            storage: ViewStorage::Owned(Mutex::new(value)),
        })
    }
    fn num_atoms(&self) -> PyResult<usize> {
        self.with_native(|native| native.num_atoms())
    }
    fn num_bonds(&self) -> PyResult<usize> {
        self.with_native(|native| native.num_bonds())
    }
    fn to_smiles(&self, py: Python<'_>) -> PyResult<String> {
        let text = self
            .with_native(|native| native.to_smiles())?
            .map_err(TautomerRunError::from_rust)?;
        decode_source_text(py, &text)
    }
    fn tautomer_score(&self) -> PyResult<TautomerScore> {
        // Only the native short loan can update its score cache. Release it
        // before converting values/errors; no Python callback runs under it.
        self.with_native(|native| native.tautomer_score())?
            .map(|inner| TautomerScore { inner })
            .map_err(TautomerRunError::from_rust)
    }
    fn atom_degree(&self, id: usize) -> PyResult<Option<usize>> {
        self.with_native(|native| native.atom_degree(ck::AtomId::new(id)))
    }
    fn atom_metadata(&self) -> PyResult<Vec<AtomMetadata>> {
        self.with_native(|native| native.atom_metadata())?
            .map(|values| {
                values
                    .into_iter()
                    .map(|inner| AtomMetadata { inner })
                    .collect()
            })
            .map_err(ValenceError::from_rust)
    }
    /// Return detached atom values, not borrowed Python views.
    fn atoms(&self) -> PyResult<Vec<Atom>> {
        self.with_native(|native| {
            native
                .atoms()
                .iter()
                .cloned()
                .map(|inner| Atom { inner })
                .collect()
        })
    }
    fn bonds(&self) -> PyResult<Vec<Bond>> {
        self.with_native(|native| {
            native
                .bonds()
                .iter()
                .cloned()
                .map(|inner| Bond { inner })
                .collect()
        })
    }
    fn atom(&self, id: usize) -> PyResult<Option<Atom>> {
        self.with_native(|native| {
            native
                .atom(ck::AtomId::new(id))
                .cloned()
                .map(|inner| Atom { inner })
        })
    }
    fn bond(&self, id: usize) -> PyResult<Option<Bond>> {
        self.with_native(|native| {
            native
                .bond(ck::BondId::new(id))
                .cloned()
                .map(|inner| Bond { inner })
        })
    }
    /// Copy the detached property value; this is not a zero-copy query.
    fn properties(&self) -> PyResult<MoleculeProperties> {
        self.with_native(|native| MoleculeProperties {
            inner: native.properties().clone(),
        })
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct TautomerProgress {
    storage: ProgressStorage,
}
#[pymethods]
impl TautomerProgress {
    fn snapshot(&self) -> PyResult<Self> {
        let native = self.storage.with_native(|value| value.snapshot())?;
        Ok(Self {
            storage: ProgressStorage::Owned(Arc::new(Mutex::new(native))),
        })
    }
    fn to_owned(&self) -> PyResult<Self> {
        let native = self.storage.with_native(|value| value.to_owned())?;
        Ok(Self {
            storage: ProgressStorage::Owned(Arc::new(Mutex::new(native))),
        })
    }
    fn len(&self) -> PyResult<usize> {
        self.storage.with_native(|native| native.len())
    }
    fn is_empty(&self) -> PyResult<bool> {
        self.storage.with_native(|native| native.is_empty())
    }
    fn status(&self) -> PyResult<TautomerEnumerationStatus> {
        self.storage.with_native(|native| native.status().into())
    }
    fn num_transforms(&self) -> PyResult<u32> {
        self.storage.with_native(|native| native.num_transforms())
    }
    fn modified_atoms(&self) -> PyResult<BTreeSet<usize>> {
        self.storage.with_native(|native| {
            native
                .modified_atoms()
                .iter()
                .map(|id| id.index())
                .collect()
        })
    }
    fn modified_bonds(&self) -> PyResult<BTreeSet<usize>> {
        self.storage.with_native(|native| {
            native
                .modified_bonds()
                .iter()
                .map(|id| id.index())
                .collect()
        })
    }
    fn entry(
        &self,
        py: Python<'_>,
        index: usize,
    ) -> PyResult<Option<(String, TautomerMoleculeView)>> {
        let key = self
            .storage
            .with_native(|native| native.entry(index).map(|(key, _)| key.clone()))?;
        key.map(|key| {
            Ok((
                decode_source_text(py, &key)?,
                TautomerMoleculeView {
                    storage: ViewStorage::Entry {
                        progress: self.storage.clone(),
                        index,
                    },
                },
            ))
        })
        .transpose()
    }
    #[doc = r#"{"role":"return_type","imports":["collections.abc"],"type":"collections.abc.Iterator[tuple[str, TautomerMoleculeView]]"}"#]
    fn entries<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        // Validate expiry now, including creation of an otherwise unused iterator.
        self.storage.with_native(|_| ())?;
        Py::new(
            py,
            Entries {
                progress: self.storage.clone(),
                index: 0,
            },
        )?
        .into_bound(py)
        .into_any()
        .cast_into::<PyIterator>()
        .map_err(Into::into)
    }
}

#[pyclass(module = "fixed_projection._native")]
struct Entries {
    progress: ProgressStorage,
    index: usize,
}
#[pymethods]
impl Entries {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<(String, TautomerMoleculeView)>> {
        let value = TautomerProgress {
            storage: self.progress.clone(),
        }
        .entry(py, self.index)?;
        if value.is_some() {
            self.index += 1;
        }
        Ok(value)
    }
}

typed_callback!(
    EnumerationCallback,
    [
        <Py<TautomerMoleculeView> as FromPyObject>::INPUT_TYPE,
        <Py<TautomerProgress> as FromPyObject>::INPUT_TYPE,
    ]
);
typed_callback!(
    Scorer,
    [<Py<TautomerMoleculeView> as FromPyObject>::INPUT_TYPE],
    pyo3::type_hint_identifier!("builtins", "int")
);

pub(crate) struct EnumerationAdapter {
    pub(crate) callback: EnumerationCallback,
    pub(crate) errors: Arc<CallbackErrors>,
}
impl ck::TautomerEnumerationCallback for EnumerationAdapter {
    fn should_continue(
        &self,
        source: ck::TautomerMoleculeView<'_>,
        progress: ck::TautomerProgress<'_>,
    ) -> Result<bool, ck::TautomerRunError> {
        let result = Python::attach(|py| {
            scope(source, Some(progress), |lease| {
                let source = Py::new(
                    py,
                    TautomerMoleculeView {
                        storage: ViewStorage::Source(lease.clone()),
                    },
                )?;
                let progress = Py::new(
                    py,
                    TautomerProgress {
                        storage: ProgressStorage::Live(lease),
                    },
                )?;
                self.callback
                    .bind(py)
                    .call1((source, progress))?
                    .extract::<bool>()
            })
        });
        result.map_err(|error| {
            self.errors.record(error);
            ck::TautomerRunError::Callback("Python callback failed".to_owned())
        })
    }
}

pub(crate) struct ScorerAdapter {
    pub(crate) callback: Scorer,
    pub(crate) errors: Arc<CallbackErrors>,
}
impl ck::TautomerScorer for ScorerAdapter {
    fn score(&self, source: ck::TautomerMoleculeView<'_>) -> Result<i32, ck::TautomerRunError> {
        let result = Python::attach(|py| {
            scope(source, None, |lease| {
                let source = Py::new(
                    py,
                    TautomerMoleculeView {
                        storage: ViewStorage::Source(lease),
                    },
                )?;
                self.callback.bind(py).call1((source,))?.extract::<i32>()
            })
        });
        result.map_err(|error| {
            self.errors.record(error);
            ck::TautomerRunError::Callback("Python scorer failed".to_owned())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::PyModule;

    #[test]
    fn callback_and_scorer_write_and_score_actual_loans_and_retained_snapshots() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CC(=O)C").unwrap();
            let expected = source.enumerate_tautomers().unwrap();
            let mut scores = std::collections::BTreeMap::new();
            for index in 0..expected.len() {
                let candidate = expected.get(index).unwrap();
                let key = decode_source_text(py, &candidate.to_smiles().unwrap())?;
                let score = candidate.tautomer_score().unwrap();
                scores.insert(
                    key,
                    (
                        score.ring(),
                        score.substructure(),
                        score.hetero_hydrogen(),
                        score.total(),
                    ),
                );
            }
            let module = PyModule::from_code(
                py,
                c"
live = []
retained = []
seen = []
def score_parts(value):
    score = value.tautomer_score()
    return (score.ring(), score.substructure(), score.hetero_hydrogen(), score.total())
def check(value):
    key = value.to_smiles()
    assert score_parts(value) == expected[key]
    snapshot = value.snapshot()
    assert snapshot.to_smiles() == key
    assert score_parts(snapshot) == expected[key]
    retained.append((key, snapshot))
    live.append(value)
    seen.append(key)
    return expected[key][-1]
def callback(source, progress):
    check(source)
    for key, value in progress.entries():
        assert value.to_smiles() == key
        check(value)
    return True
def scorer(value):
    return check(value)
def verify_after():
    for key, snapshot in retained:
        assert snapshot.to_smiles() == key
        assert score_parts(snapshot) == expected[key]
    for value in live:
        for method in ('to_smiles', 'tautomer_score'):
            try:
                getattr(value, method)()
            except ReferenceError:
                pass
            else:
                raise AssertionError('expired native write/score loan remained usable')
    assert set(seen) == set(expected)
",
                c"<tautomer-native-methods>",
                c"tautomer_native_methods",
            )?;
            module.setattr("expected", scores)?;
            let errors = Arc::new(CallbackErrors::new());
            let mut params = ck::TautomerParams::default();
            params.set_callback(Some(Arc::new(EnumerationAdapter {
                callback: module
                    .getattr("callback")?
                    .extract::<EnumerationCallback>()?,
                errors: errors.clone(),
            })));
            let enumeration = source.enumerate_tautomers_with_params(&params);
            if let Some(error) = errors.take() {
                return Err(error);
            }
            assert_eq!(enumeration.unwrap().len(), expected.len());
            params.set_callback(None);
            params.set_scorer(Some(Arc::new(ScorerAdapter {
                callback: module.getattr("scorer")?.extract::<Scorer>()?,
                errors: errors.clone(),
            })));
            let canonical = source.canonical_tautomer_with_params(&params);
            if let Some(error) = errors.take() {
                return Err(error);
            }
            assert_eq!(canonical.unwrap().num_atoms(), 4);
            drop(params);
            drop(source);
            drop(expected);
            module.getattr("verify_after")?.call0()?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn actual_scorer_calls_expire_individually_and_keep_snapshots_and_original_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::from_code(
                py,
                c"
live = []
snapshots = []
failure = None
marker = LookupError('scorer marker')
def verify_expired():
    for previous in live:
        try:
            previous.num_atoms()
        except ReferenceError:
            pass
        else:
            raise AssertionError('each scoring loan must expire before the next call')
def scorer(candidate):
    verify_expired()
    live.append(candidate)
    snapshots.append(candidate.snapshot())
    if failure is not None:
        raise failure
    return candidate.num_atoms()
def verify_after():
    verify_expired()
    assert len(snapshots) >= 2
    assert all(snapshot.num_atoms() == 4 for snapshot in snapshots)
",
                c"<tautomer-scorer>",
                c"scorer_lease_test",
            )?;
            for should_fail in [false, true] {
                if should_fail {
                    module.setattr("failure", module.getattr("marker")?)?;
                }
                // A new invocation owns its own error slot, even when it uses
                // the same Python callable and retains earlier snapshots.
                let errors = Arc::new(CallbackErrors::new());
                let mut params = ck::TautomerParams::default();
                params.set_scorer(Some(Arc::new(ScorerAdapter {
                    callback: module.getattr("scorer")?.extract::<Scorer>()?,
                    errors: errors.clone(),
                })));
                let source = ck::Molecule::from_smiles("CC(=O)C").unwrap();
                let result = source.canonical_tautomer_with_params(&params);
                if should_fail {
                    assert!(matches!(
                        result,
                        Err(ck::OperationError::Tautomer(
                            ck::TautomerRunError::Callback(_)
                        ))
                    ));
                    let error = errors.take().expect("retain original scoring failure");
                    assert!(error.value(py).is(&module.getattr("marker")?));
                } else {
                    if let Some(error) = errors.take() {
                        return Err(error);
                    }
                    assert_eq!(result.unwrap().num_atoms(), 4);
                }
                drop(params);
                drop(source);
                module.getattr("verify_after")?.call0()?;
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn live_loans_reject_foreign_threads_but_owned_snapshots_keep_native_thread_behavior() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::from_code(
                py,
                c"
from threading import Thread
retained = []
failures = []
results = []
def during(live, snapshot, progress):
    try:
        try:
            live.num_atoms()
        except RuntimeError:
            pass
        else:
            raise AssertionError('foreign thread accessed a borrowed native view')
        results.extend([snapshot.num_atoms(), progress.entry(0)[1].num_atoms(),
                        next(progress.entries())[1].num_atoms()])
    except BaseException as error:
        failures.append(error)
def callback(source, progress):
    snapshot = source.snapshot()
    owned_progress = progress.snapshot()
    retained.extend([snapshot, owned_progress, owned_progress.entry(0)[1], owned_progress.entries(),
                     source, progress, progress.entry(0)[1], progress.entries()])
    worker = Thread(target=during, args=(source, snapshot, owned_progress))
    worker.start()
    worker.join()
    assert not failures, failures
    return False
def after():
    try:
        for live, method in [(retained[4], 'num_atoms'), (retained[5], 'len'),
                             (retained[6], 'num_atoms'), (retained[7], '__next__')]:
            try:
                getattr(live, method)()
            except ReferenceError:
                pass
            else:
                raise AssertionError('expired loan must expire on every thread')
        results.extend([retained[0].num_atoms(), retained[1].len(),
                        retained[2].num_atoms(), next(retained[3])[1].num_atoms()])
    except BaseException as error:
        failures.append(error)
def verify_after():
    worker = Thread(target=after)
    worker.start()
    worker.join()
    assert not failures, failures
    assert results[:3] == [4, 4, 4]
    assert results[3] == results[5] == results[6] == 4
    assert results[4] > 0
",
                c"<tautomer-thread>",
                c"lease_thread_test",
            )?;
            let errors = Arc::new(CallbackErrors::new());
            let callback = module
                .getattr("callback")?
                .extract::<EnumerationCallback>()?;
            let mut params = ck::TautomerParams::default();
            params.set_callback(Some(Arc::new(EnumerationAdapter {
                callback,
                errors: errors.clone(),
            })));
            let result = ck::Molecule::from_smiles("CC(=O)C")
                .unwrap()
                .enumerate_tautomers_with_params(&params);
            if let Some(error) = errors.take() {
                return Err(error);
            }
            assert_eq!(
                result.unwrap().status(),
                ck::TautomerEnumerationStatus::Canceled
            );
            module.getattr("verify_after")?.call0()?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn actual_callback_loans_expire_but_explicit_snapshots_and_detached_values_survive() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::from_code(
                py,
                c"
saved = []
def callback(source, progress):
    iterator = progress.entries()
    key, candidate = next(iterator)
    assert source.num_atoms() == candidate.num_atoms()
    assert progress.len() > 0  # The iterator must not keep progress borrowed.
    assert progress.entry(0)[0] == key
    saved.extend([source, progress, iterator, candidate,
                  source.snapshot(), progress.snapshot(), candidate.snapshot(),
                  source.atom(0), source.properties()])
    return False
",
                c"<tautomer-lease>",
                c"lease_test",
            )?;
            let callback = module
                .getattr("callback")?
                .extract::<EnumerationCallback>()?;
            let errors = Arc::new(CallbackErrors::new());
            let mut params = ck::TautomerParams::default();
            params.set_callback(Some(Arc::new(EnumerationAdapter {
                callback,
                errors: errors.clone(),
            })));
            let source = ck::Molecule::from_smiles("CC(=O)C").unwrap();
            let result = source.enumerate_tautomers_with_params(&params).unwrap();
            assert_eq!(result.status(), ck::TautomerEnumerationStatus::Canceled);
            assert!(errors.take().is_none());
            let saved = module.getattr("saved")?;
            for (index, method) in [
                (0, "num_atoms"),
                (1, "len"),
                (2, "__next__"),
                (3, "num_atoms"),
            ] {
                let error = saved.get_item(index)?.call_method0(method).unwrap_err();
                assert!(error.is_instance_of::<PyReferenceError>(py));
            }
            assert_eq!(
                saved
                    .get_item(4)?
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                4
            );
            let owned = saved.get_item(5)?;
            assert!(owned.call_method0("len")?.extract::<usize>()? > 0);
            let child = owned.call_method1("entry", (0,))?.get_item(1)?;
            let owned_iterator = owned.call_method0("entries")?;
            saved.del_item(5)?;
            drop(owned);
            assert_eq!(child.call_method0("num_atoms")?.extract::<usize>()?, 4);
            assert_eq!(
                owned_iterator
                    .call_method0("__next__")?
                    .get_item(1)?
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                4
            );
            // Removing index 5 shifts the retained candidate/atom/property values.
            assert_eq!(
                saved
                    .get_item(5)?
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                4
            );
            assert_eq!(
                saved
                    .get_item(6)?
                    .call_method0("atomic_number")?
                    .extract::<u8>()?,
                6
            );
            saved.get_item(7)?.call_method0("name")?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn actual_callback_error_keeps_original_exception_and_expires_aliases() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::from_code(
                py,
                c"
saved = []
marker = ValueError('callback marker')
def callback(source, progress):
    saved.extend([source, progress, progress.entries(), progress.entry(0)[1]])
    raise marker
",
                c"<tautomer-error>",
                c"lease_error_test",
            )?;
            let errors = Arc::new(CallbackErrors::new());
            let callback = module
                .getattr("callback")?
                .extract::<EnumerationCallback>()?;
            let mut params = ck::TautomerParams::default();
            params.set_callback(Some(Arc::new(EnumerationAdapter {
                callback,
                errors: errors.clone(),
            })));
            let native = ck::Molecule::from_smiles("CC(=O)C")
                .unwrap()
                .enumerate_tautomers_with_params(&params)
                .unwrap_err();
            assert!(matches!(
                native,
                ck::OperationError::Tautomer(ck::TautomerRunError::Callback(_))
            ));
            let original = errors.take().expect("retain actual Python failure");
            assert!(original.value(py).is(&module.getattr("marker")?));
            for (index, method) in [
                (0, "num_atoms"),
                (1, "len"),
                (2, "__next__"),
                (3, "num_atoms"),
            ] {
                let error = module
                    .getattr("saved")?
                    .get_item(index)?
                    .call_method0(method)
                    .unwrap_err();
                assert!(error.is_instance_of::<PyReferenceError>(py));
            }
            Ok(())
        })
        .unwrap();
    }
}
