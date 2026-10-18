//! Tautomer configuration forwards native policy and reuses scoped callbacks.
use crate::callbacks::CallbackErrors;
use crate::tautomer_callbacks::{EnumerationAdapter, EnumerationCallback, Scorer, ScorerAdapter};
use crate::tautomer_values::TautomerScoreParams;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyclass(module = "fixed_projection._native")]
pub struct TautomerParams {
    pub(crate) inner: ck::TautomerParams,
    callback: Option<EnumerationCallback>,
    scorer: Option<Scorer>,
}

// These methods only forward native policy access. Callback ownership and
// invocation errors use the same mechanism as SubstructMatchParams.
macro_rules! policy_fields {
    ($($field:ident: $ty:ty),+ $(,)?) => {
        paste::paste! {
            #[pymethods]
            impl TautomerParams {
                #[new]
                #[pyo3(signature = (*, $($field=ck::TautomerParams::default().$field(),)+ callback=None, scorer=None, score_params=None))]
                pub(crate) fn new($($field: $ty,)+ callback: Option<EnumerationCallback>, scorer: Option<Scorer>, score_params: Option<&TautomerScoreParams>) -> Self {
                    let mut inner = ck::TautomerParams::default();
                    $(inner.[<set_ $field>]($field);)+
                    if let Some(params) = score_params { inner.score_params = params.inner.clone(); }
                    Self { inner, callback, scorer }
                }
                $(
                    #[getter]
                    fn $field(&self) -> $ty { self.inner.$field() }
                    #[setter($field)]
                    fn [<assign_ $field>](&mut self, value: $ty) { self.inner.[<set_ $field>](value); }
                    fn [<set_ $field>](&mut self, value: $ty) { self.inner.[<set_ $field>](value); }
                    fn [<with_ $field>](&self, py: Python<'_>, value: $ty) -> Self {
                        Self {
                            inner: self.inner.clone().[<with_ $field>](value),
                            callback: self.callback.as_ref().map(|callback| callback.clone_ref(py)),
                            scorer: self.scorer.as_ref().map(|scorer| scorer.clone_ref(py)),
                        }
                    }
                )+
                fn transform_count(&self) -> usize { self.inner.transform_count() }
                #[staticmethod]
                fn v1() -> PyResult<Self> { crate::tautomer_catalog::v1() }
                #[staticmethod]
                fn from_transform_data(data: Vec<(String, String, String, String)>) -> PyResult<Self> {
                    crate::tautomer_catalog::from_transform_data(data)
                }
                #[staticmethod]
                fn from_transform_file(path: std::path::PathBuf) -> PyResult<Self> {
                    crate::tautomer_catalog::from_transform_file(path)
                }
                #[getter]
                fn score_params(&self) -> TautomerScoreParams { TautomerScoreParams { inner: self.inner.score_params.clone() } }
                #[setter]
                fn set_score_params(&mut self, value: Option<&TautomerScoreParams>) {
                    self.inner.score_params = value.map(|params| params.inner.clone()).unwrap_or_else(|| ck::TautomerParams::default().score_params);
                }
                #[getter]
                fn callback(&self, py: Python<'_>) -> Option<EnumerationCallback> {
                    self.callback.as_ref().map(|callback| callback.clone_ref(py))
                }
                #[setter(callback)]
                fn assign_callback(&mut self, value: Option<EnumerationCallback>) { self.set_callback(value); }
                fn set_callback(&mut self, callback: Option<EnumerationCallback>) { self.callback = callback; }
                #[getter]
                fn scorer(&self, py: Python<'_>) -> Option<Scorer> {
                    self.scorer.as_ref().map(|scorer| scorer.clone_ref(py))
                }
                #[setter(scorer)]
                fn assign_scorer(&mut self, value: Option<Scorer>) { self.set_scorer(value); }
                fn set_scorer(&mut self, scorer: Option<Scorer>) { self.scorer = scorer; }
                fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
                    if let Some(callback) = &self.callback { callback.traverse(&visit)?; }
                    if let Some(scorer) = &self.scorer { scorer.traverse(&visit)?; }
                    Ok(())
                }
                fn __clear__(&mut self) { self.callback = None; self.scorer = None; }
            }
        }
    };
}
policy_fields! {
    max_tautomers: u32,
    max_transforms: u32,
    remove_sp3_stereo: bool,
    remove_bond_stereo: bool,
    remove_isotopic_hydrogens: bool,
    reassign_stereo: bool,
}

impl TautomerParams {
    pub(crate) fn clone_for_call(&self, py: Python<'_>) -> Self {
        Self {
            inner: self.inner.clone(),
            callback: self.callback.as_ref().map(|value| value.clone_ref(py)),
            scorer: self.scorer.as_ref().map(|value| value.clone_ref(py)),
        }
    }
    pub(crate) fn from_native(inner: ck::TautomerParams) -> Self {
        Self {
            inner,
            callback: None,
            scorer: None,
        }
    }
    pub(crate) fn with_callbacks<T, E>(
        &self,
        py: Python<'_>,
        run: impl FnOnce(&ck::TautomerParams) -> Result<T, E>,
    ) -> PyResult<Result<T, E>> {
        let mut params = self.inner.clone();
        let errors = Arc::new(CallbackErrors::new());
        if let Some(callback) = &self.callback {
            params.set_callback(Some(Arc::new(EnumerationAdapter {
                callback: callback.clone_ref(py),
                errors: Arc::clone(&errors),
            })));
        }
        if let Some(scorer) = &self.scorer {
            params.set_scorer(Some(Arc::new(ScorerAdapter {
                callback: scorer.clone_ref(py),
                errors: Arc::clone(&errors),
            })));
        }
        let result = run(&params);
        match errors.take() {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::ffi::c_str;
    use pyo3::types::PyModule;

    #[test]
    fn complete_policy_configuration_preserves_native_values_callbacks_and_owned_copies() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "tautomer_configuration_test")?;
            module.add_class::<TautomerParams>()?;
            module.add_class::<TautomerScoreParams>()?;
            py.run(
                c_str!(
                    r#"
import gc
import weakref
defaults = TautomerParams()
params = TautomerParams(max_tautomers=23, max_transforms=29, remove_sp3_stereo=False,
    remove_bond_stereo=False, remove_isotopic_hydrogens=False, reassign_stereo=False,
    score_params=TautomerScoreParams(terms=[]))
for field, value in [('max_tautomers', 31), ('max_transforms', 37),
    ('remove_sp3_stereo', True), ('remove_bond_stereo', True),
    ('remove_isotopic_hydrogens', True), ('reassign_stereo', True)]:
    original = getattr(params, field)
    copied = getattr(params, 'with_' + field)(value)
    assert getattr(copied, field) == value and getattr(params, field) == original
    getattr(params, 'set_' + field)(value)
    assert getattr(params, field) == value
    setattr(params, field, original)
    assert getattr(params, field) == original
assert params.score_params.terms == []
params.score_params = None
assert params.score_params.terms is None
try:
    params.max_tautomers = -1
except OverflowError:
    pass
else:
    raise AssertionError('negative policy count accepted')
assert params.max_tautomers == 23
events = []
def callback(source, progress):
    events.append((source.num_atoms(), progress.len()))
    return True
params.callback = callback
assert params.callback is callback
assert params.with_max_tautomers(7).callback is callback
params.set_callback(None)
assert params.callback is None
params.set_callback(callback)
def scorer(source):
    return source.num_atoms()
params.scorer = scorer
assert params.scorer is scorer
assert params.with_reassign_stereo(True).scorer is scorer
params.set_scorer(None)
assert params.scorer is None
class Callback:
    def __init__(self, params): self.params = params
    def __call__(self, source, progress): return True
cycle = TautomerParams()
callback_cycle = Callback(cycle)
observer = weakref.ref(callback_cycle)
cycle.callback = callback_cycle
del callback_cycle, cycle
gc.collect()
assert observer() is None
class ExpectedError(Exception): pass
sentinel = ExpectedError('actual callback exception')
def failing(source, progress): raise sentinel
"#
                ),
                Some(&module.dict()),
                None,
            )?;
            let default = module
                .getattr("defaults")?
                .extract::<PyRef<TautomerParams>>()?;
            let native_default = ck::TautomerParams::default();
            assert_eq!(
                default.inner.max_tautomers(),
                native_default.max_tautomers()
            );
            assert_eq!(
                default.inner.max_transforms(),
                native_default.max_transforms()
            );
            assert_eq!(
                default.inner.remove_sp3_stereo(),
                native_default.remove_sp3_stereo()
            );
            assert_eq!(
                default.inner.remove_bond_stereo(),
                native_default.remove_bond_stereo()
            );
            assert_eq!(
                default.inner.remove_isotopic_hydrogens(),
                native_default.remove_isotopic_hydrogens()
            );
            assert_eq!(
                default.inner.reassign_stereo(),
                native_default.reassign_stereo()
            );
            assert_eq!(default.transform_count(), native_default.transform_count());
            let params_object = module.getattr("params")?;
            let params = params_object.extract::<PyRef<TautomerParams>>()?;
            let molecule = ck::Molecule::from_smiles("CC(=O)C").unwrap();
            params
                .with_callbacks(py, |native| {
                    py.detach(|| molecule.enumerate_tautomers_with_params(native))
                })?
                .unwrap();
            assert!(!module.getattr("events")?.is_empty()?);
            drop(params);
            params_object.setattr("callback", module.getattr("failing")?)?;
            let params = params_object.extract::<PyRef<TautomerParams>>()?;
            let error = params
                .with_callbacks(py, |native| {
                    py.detach(|| molecule.enumerate_tautomers_with_params(native))
                })
                .unwrap_err();
            assert!(error.value(py).is(&module.getattr("sentinel")?));
            Ok(())
        })
        .unwrap();
    }
}
