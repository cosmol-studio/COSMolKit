use crate::atom_bond_values::{Atom, Bond};
use crate::callbacks::{CallbackErrors, typed_callback};
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::sync::Arc;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct QueryAtom {
        pub(crate) inner: ck::QueryAtom,
    }
}
#[pymethods]
impl QueryAtom {
    #[getter]
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[getter]
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    #[getter]
    fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    #[getter]
    fn isotope(&self) -> Option<u16> {
        self.inner.isotope()
    }
    #[getter]
    fn atom_map(&self) -> Option<u32> {
        self.inner.atom_map()
    }
    #[getter]
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    #[getter]
    fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    #[getter]
    fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
}

typed_callback!(
    FinalMatch,
    [
        <Py<Molecule> as FromPyObject>::INPUT_TYPE,
        pyo3::type_hint_subscript!(
            pyo3::type_hint_identifier!("collections.abc", "Sequence"),
            pyo3::type_hint_identifier!("builtins", "int")
        )
    ]
);
typed_callback!(
    AtomMatch,
    [
        <Py<QueryAtom> as FromPyObject>::INPUT_TYPE,
        <Py<Atom> as FromPyObject>::INPUT_TYPE
    ]
);
typed_callback!(
    BondMatch,
    [
        <Py<Bond> as FromPyObject>::INPUT_TYPE,
        <Py<Bond> as FromPyObject>::INPUT_TYPE
    ]
);

#[pyclass(module = "cosmolkit")]
pub struct SubstructMatchParams {
    // Actual native inputs plus collector-visible Python callback owners.
    inner: ck::SubstructMatchParams,
    final_match: Option<FinalMatch>,
    atom_match: Option<AtomMatch>,
    bond_match: Option<BondMatch>,
}

/// One invocation owns its extracted native inputs and callback handles.
pub(crate) struct MatchInvocation {
    inner: ck::SubstructMatchParams,
    final_match: Option<FinalMatch>,
    atom_match: Option<AtomMatch>,
    bond_match: Option<BondMatch>,
}
impl MatchInvocation {
    fn new(
        inner: ck::SubstructMatchParams,
        final_match: Option<FinalMatch>,
        atom_match: Option<AtomMatch>,
        bond_match: Option<BondMatch>,
        extra_final_check: Option<FinalMatch>,
        extra_atom_check: Option<AtomMatch>,
        extra_bond_check: Option<BondMatch>,
    ) -> PyResult<Self> {
        if final_match.is_some() && extra_final_check.is_some()
            || atom_match.is_some() && extra_atom_check.is_some()
            || bond_match.is_some() && extra_bond_check.is_some()
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "callback and its alias are mutually exclusive",
            ));
        }
        Ok(Self {
            inner,
            final_match: final_match.or(extra_final_check),
            atom_match: atom_match.or(extra_atom_check),
            bond_match: bond_match.or(extra_bond_check),
        })
    }
    pub(crate) fn with_callbacks<T, E>(
        self,
        py: Python<'_>,
        target: &Py<Molecule>,
        run: impl FnOnce(&ck::SubstructMatchParams) -> Result<T, E>,
    ) -> PyResult<Result<T, E>> {
        let target = target.try_borrow(py)?;
        self.with_callback_values(py, &[&target.inner], run)
    }
    pub(crate) fn with_callback_values<T, E>(
        self,
        py: Python<'_>,
        targets: &[&ck::Molecule],
        run: impl FnOnce(&ck::SubstructMatchParams) -> Result<T, E>,
    ) -> PyResult<Result<T, E>> {
        let mut params = self.inner;
        let errors = Arc::new(CallbackErrors::new());
        if let Some(callback) = &self.final_match {
            let callback = callback.clone_ref(py);
            // Retain native COW values, never aliases of the Python receiver.
            // Each callback then owns its own value, including when retained
            // or mutated after the matcher and its input owners have returned.
            let targets: Vec<_> = targets.iter().map(|target| (*target).clone()).collect();
            let errors = Arc::clone(&errors);
            params.extra_final_check = Some(Arc::new(move |target, indices| {
                errors.call(|py| {
                    let index = if targets.len() == 1 {
                        0
                    } else {
                        target.input_index().ok_or_else(|| {
                            pyo3::exceptions::PyRuntimeError::new_err(
                                "native multi-input matcher omitted its input index",
                            )
                        })?
                    };
                    let target = targets.get(index).ok_or_else(|| {
                        pyo3::exceptions::PyRuntimeError::new_err(
                            "native matcher input index is outside this invocation",
                        )
                    })?;
                    callback
                        .bind(py)
                        .call1((
                            Py::new(
                                py,
                                Molecule {
                                    inner: target.clone(),
                                },
                            )?,
                            indices.to_vec(),
                        ))?
                        .extract::<bool>()
                })
            }));
        }
        if let Some(callback) = &self.atom_match {
            let callback = callback.clone_ref(py);
            let errors = Arc::clone(&errors);
            params.extra_atom_check = Some(Arc::new(move |_, query, _, atom| {
                errors.call(|py| {
                    callback
                        .bind(py)
                        .call1((
                            Py::new(
                                py,
                                QueryAtom {
                                    inner: query.clone(),
                                },
                            )?,
                            Py::new(
                                py,
                                Atom {
                                    inner: atom.clone(),
                                },
                            )?,
                        ))?
                        .extract::<bool>()
                })
            }));
        }
        if let Some(callback) = &self.bond_match {
            let callback = callback.clone_ref(py);
            let errors = Arc::clone(&errors);
            params.extra_bond_check = Some(Arc::new(move |query, target| {
                errors.call(|py| {
                    callback
                        .bind(py)
                        .call1((
                            Py::new(
                                py,
                                Bond {
                                    inner: query.clone(),
                                },
                            )?,
                            Py::new(
                                py,
                                Bond {
                                    inner: target.clone(),
                                },
                            )?,
                        ))?
                        .extract::<bool>()
                })
            }));
        }
        let result = run(&params);
        match errors.take() {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }
}
macro_rules! match_parameters {
    ($($field:ident: $ty:ty),+ $(,)?) => {
        impl SubstructMatchParams {
            pub(crate) fn from_keywords($($field: Option<$ty>,)+ atom_properties: Option<Vec<String>>, bond_properties: Option<Vec<String>>, final_match: Option<FinalMatch>, atom_match: Option<AtomMatch>, bond_match: Option<BondMatch>, extra_final_check: Option<FinalMatch>, extra_atom_check: Option<AtomMatch>, extra_bond_check: Option<BondMatch>) -> PyResult<MatchInvocation> {
                let defaults = ck::SubstructMatchParams::default();
                MatchInvocation::new(ck::SubstructMatchParams {
                    $($field: $field.unwrap_or(defaults.$field),)+
                    atom_properties: atom_properties.unwrap_or(defaults.atom_properties),
                    bond_properties: bond_properties.unwrap_or(defaults.bond_properties),
                    extra_atom_check: None, extra_bond_check: None, extra_final_check: None,
                }, final_match, atom_match, bond_match, extra_final_check, extra_atom_check, extra_bond_check)
            }
        }
        paste::paste! {
            #[pymethods]
            impl SubstructMatchParams {
                #[new]
                #[pyo3(signature=(*, $($field=ck::SubstructMatchParams::default().$field,)+ atom_properties=None, bond_properties=None, final_match=None, atom_match=None, bond_match=None, extra_final_check=None, extra_atom_check=None, extra_bond_check=None))]
                pub(crate) fn new(py: Python<'_>, $($field: $ty,)+ atom_properties: Option<Vec<String>>, bond_properties: Option<Vec<String>>, final_match: Option<FinalMatch>, atom_match: Option<AtomMatch>, bond_match: Option<BondMatch>, extra_final_check: Option<FinalMatch>, extra_atom_check: Option<AtomMatch>, extra_bond_check: Option<BondMatch>) -> PyResult<Self> {
                    let invocation = Self::from_keywords($(Some($field),)+ atom_properties, bond_properties, final_match, atom_match, bond_match, extra_final_check, extra_atom_check, extra_bond_check)?;
                    Self::from_invocation(py, invocation)
                }
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
                #[getter] fn atom_properties(&self) -> crate::parameters::ContainerOutput<'_, Vec<String>> { crate::parameters::ContainerOutput(&self.inner.atom_properties) }
                #[setter] fn set_atom_properties(&mut self, value: Option<Vec<String>>) -> PyResult<()> {
                    self.inner.atom_properties = value.unwrap_or_else(|| ck::SubstructMatchParams::default().atom_properties);
                    Ok(())
                }
                #[getter] fn bond_properties(&self) -> crate::parameters::ContainerOutput<'_, Vec<String>> { crate::parameters::ContainerOutput(&self.inner.bond_properties) }
                #[setter] fn set_bond_properties(&mut self, value: Option<Vec<String>>) -> PyResult<()> {
                    self.inner.bond_properties = value.unwrap_or_else(|| ck::SubstructMatchParams::default().bond_properties);
                    Ok(())
                }
                #[getter] fn final_match(&self, py: Python<'_>) -> Option<FinalMatch> { self.final_match.as_ref().map(|value| value.clone_ref(py)) }
                #[setter] fn set_final_match(&mut self, value: Option<FinalMatch>) { self.final_match = value; }
                #[getter] fn atom_match(&self, py: Python<'_>) -> Option<AtomMatch> { self.atom_match.as_ref().map(|value| value.clone_ref(py)) }
                #[setter] fn set_atom_match(&mut self, value: Option<AtomMatch>) { self.atom_match = value; }
                #[getter] fn bond_match(&self, py: Python<'_>) -> Option<BondMatch> { self.bond_match.as_ref().map(|value| value.clone_ref(py)) }
                #[setter] fn set_bond_match(&mut self, value: Option<BondMatch>) { self.bond_match = value; }
                #[getter] fn extra_final_check(&self, py: Python<'_>) -> Option<FinalMatch> { self.final_match(py) }
                #[setter] fn set_extra_final_check(&mut self, value: Option<FinalMatch>) { self.set_final_match(value); }
                #[getter] fn extra_atom_check(&self, py: Python<'_>) -> Option<AtomMatch> { self.atom_match(py) }
                #[setter] fn set_extra_atom_check(&mut self, value: Option<AtomMatch>) { self.set_atom_match(value); }
                #[getter] fn extra_bond_check(&self, py: Python<'_>) -> Option<BondMatch> { self.bond_match(py) }
                #[setter] fn set_extra_bond_check(&mut self, value: Option<BondMatch>) { self.set_bond_match(value); }
                fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
                    if let Some(value) = &self.final_match { value.traverse(&visit)?; }
                    if let Some(value) = &self.atom_match { value.traverse(&visit)?; }
                    if let Some(value) = &self.bond_match { value.traverse(&visit)?; }
                    Ok(())
                }
                fn __clear__(&mut self) {
                    self.final_match = None; self.atom_match = None; self.bond_match = None;
                }
            }
        }
    };
}
ck::project_configuration!(ck::SubstructMatchParams, using match_parameters);

impl SubstructMatchParams {
    fn from_invocation(py: Python<'_>, value: MatchInvocation) -> PyResult<Self> {
        let mut result = Self::from_native(py, value.inner)?;
        result.final_match = value.final_match;
        result.atom_match = value.atom_match;
        result.bond_match = value.bond_match;
        Ok(result)
    }
    pub(crate) fn from_native(_py: Python<'_>, inner: ck::SubstructMatchParams) -> PyResult<Self> {
        Ok(Self {
            inner,
            final_match: None,
            atom_match: None,
            bond_match: None,
        })
    }
    pub(crate) fn to_rust(&self, _py: Python<'_>) -> PyResult<ck::SubstructMatchParams> {
        Ok(self.inner.clone())
    }
    pub(crate) fn capture(&self, py: Python<'_>) -> PyResult<MatchInvocation> {
        Ok(MatchInvocation {
            inner: self.to_rust(py)?,
            final_match: self.final_match.as_ref().map(|value| value.clone_ref(py)),
            atom_match: self.atom_match.as_ref().map(|value| value.clone_ref(py)),
            bond_match: self.bond_match.as_ref().map(|value| value.clone_ref(py)),
        })
    }
    pub(crate) fn with_native(
        &self,
        py: Python<'_>,
        inner: ck::SubstructMatchParams,
    ) -> PyResult<Self> {
        Self::from_invocation(
            py,
            MatchInvocation {
                inner,
                final_match: self.final_match.as_ref().map(|value| value.clone_ref(py)),
                atom_match: self.atom_match.as_ref().map(|value| value.clone_ref(py)),
                bond_match: self.bond_match.as_ref().map(|value| value.clone_ref(py)),
            },
        )
    }
    pub(crate) fn has_final_match(&self) -> bool {
        self.final_match.is_some()
    }
    pub(crate) fn with_callbacks<T, E>(
        &self,
        py: Python<'_>,
        target: &Py<Molecule>,
        run: impl FnOnce(&ck::SubstructMatchParams) -> Result<T, E>,
    ) -> PyResult<Result<T, E>> {
        self.capture(py)?.with_callbacks(py, target, run)
    }
    pub(crate) fn with_callback_values<T, E>(
        &self,
        py: Python<'_>,
        targets: &[&ck::Molecule],
        run: impl FnOnce(&ck::SubstructMatchParams) -> Result<T, E>,
    ) -> PyResult<Result<T, E>> {
        self.capture(py)?.with_callback_values(py, targets, run)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::ffi::c_str;
    use pyo3::types::PyModule;

    #[test]
    fn molecule_copy_and_callback_values_share_blocks_and_isolate_writes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let original = ck::Molecule::from_smiles("CCC").unwrap();
            let target = Py::new(
                py,
                Molecule {
                    inner: original.clone(),
                },
            )?;
            let copied = target
                .bind(py)
                .call_method0("copy")?
                .extract::<Py<Molecule>>()?;
            assert!(!copied.is(&target));
            let source_blocks = original.to_builder();
            let copied_blocks = copied.borrow(py).inner.to_builder();
            assert!(std::ptr::eq(
                original.topology(),
                copied.borrow(py).inner.topology()
            ));
            assert!(std::ptr::eq(
                original.properties(),
                copied.borrow(py).inner.properties()
            ));
            assert!(std::ptr::eq(
                source_blocks.coordinates(),
                copied_blocks.coordinates()
            ));
            let module = PyModule::new(py, "cow_callback_probe")?;
            module.add_class::<Molecule>()?;
            module.add_class::<crate::query_compiler::QueryGraph>()?;
            module.add_class::<SubstructMatchParams>()?;
            let scope = pyo3::types::PyDict::new(py);
            scope.set_item("ck", module)?;
            scope.set_item("mol", &target)?;
            scope.set_item("copied", &copied)?;
            py.run(
                c_str!(
                    r#"
import copy
other = copy.copy(mol)
assert other is not mol and other.to_binary() == mol.to_binary()
copied.add_hydrogens_()
assert copied.num_atoms() == 11 and mol.num_atoms() == 3
seen = []
def accept(target, indices):
    assert target is not mol and target.num_atoms() == 3
    target.add_hydrogens_()
    assert target.num_atoms() == 11
    seen.append(target)
    return True
matches = mol.substruct_matches(ck.QueryGraph.from_smarts("C"), final_match=accept)
assert [value.atom_mapping for value in matches] == [[0], [1], [2]]
assert len(seen) == 3 and len({id(value) for value in seen}) == 3
assert mol.num_atoms() == 3
seen[0].remove_hydrogens_()
assert seen[0].num_atoms() == 3 and seen[1].num_atoms() == 11
mol.add_hydrogens_()
assert other.num_atoms() == 3
del mol, copied, other
assert seen[2].num_atoms() == 11
"#
                ),
                Some(&scope),
                None,
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn public_matchers_forward_configurations_callbacks_and_collect_cycles() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "matcher_producer_test")?;
            module.add_class::<Molecule>()?;
            module.add_class::<crate::query_compiler::QueryGraph>()?;
            module.add_class::<SubstructMatchParams>()?;
            module.add_function(wrap_pyfunction!(crate::query_compiler::compile_query, &module)?)?;
            let scope = pyo3::types::PyDict::new(py);
            scope.set_item("ck", &module)?;
            py.run(c_str!(r#"
import gc
mol = ck.Molecule.from_smiles("CCCC")
query = ck.QueryGraph.from_smarts("C")
assert [match.atom_mapping for match in mol.substruct_matches(query)] == [[0], [1], [2], [3]]
assert mol.substruct_match(query).atom_mapping == [0]
assert mol.has_substruct_match(query) is True
params = ck.SubstructMatchParams(max_matches=2)
assert [match.atom_mapping for match in mol.substruct_matches(query, params)] == [[0], [1]]
assert [match.atom_mapping for match in mol.substruct_matches(query, max_matches=2)] == [[0], [1]]
assert [match.atom_mapping for match in mol.substruct_matches_with_params(query, params)] == [[0], [1]]
assert mol.substruct_match(query, params).atom_mapping == [0]
assert mol.substruct_match(query, max_matches=2).atom_mapping == [0]
assert mol.substruct_match_with_params(query, params).atom_mapping == [0]
assert mol.has_substruct_match(query, params) is True
assert mol.has_substruct_match(query, max_matches=2) is True
assert mol.has_substruct_match_with_params(query, params) is True
assert [match.atom_mapping for match in mol.substruct_matches_compiled(ck.compile_query(query))] == [[0], [1], [2], [3]]
assert [match.atom_mapping for match in mol.substruct_matches(query, max_matches=0)] == [[0], [1], [2], [3]]
seen = []
def final_match(target, indices):
    assert target is not mol
    assert target.num_atoms() == mol.num_atoms()
    seen.append(tuple(indices))
    return True
params.final_match = final_match
for method in (mol.substruct_matches, mol.substruct_match, mol.has_substruct_match):
    seen.clear()
    method(query, params)
    assert seen == [(0,), (1,)]
    seen.clear()
    method(query, max_matches=2, final_match=final_match)
    assert seen == [(0,), (1,)]
    try:
        method(query, params, use_chirality=False)
    except TypeError:
        pass
    else:
        raise AssertionError("mixed configuration accepted")
marker = LookupError("original producer callback failure")
def fail(target, indices):
    raise marker
for method in (mol.substruct_matches, mol.substruct_match, mol.has_substruct_match):
    try:
        method(query, final_match=fail)
    except LookupError as error:
        assert error is marker
    else:
        raise AssertionError("callback exception lost")
collected = []
class Callback:
    def __init__(self, params):
        self.params = params
    def __call__(self, target, indices):
        return True
    def __del__(self):
        collected.append(True)
def create_cycle():
    params = ck.SubstructMatchParams()
    params.final_match = Callback(params)
create_cycle()
gc.collect()
assert collected == [True]
"#), Some(&scope), None)?;
            let error = crate::match_errors::MatchError::from_rust(
                ck::MatchError::Substruct(ck::SubstructMatchError::PropertyInteger {
                    property: "count", source: ck::PropertyIntReadError::InvalidKind { kind: ck::PropertyValueKind::Bool },
                })
            );
            assert!(error.value(py).getattr("substruct")?.is_instance_of::<crate::match_errors::SubstructMatchError>());
            assert!(error.value(py).getattr("source")?.is(&error.value(py).getattr("__cause__")?));
            assert!(error.cause(py).unwrap().is_instance_of::<crate::property_numeric_errors::PropertyIntReadError>(py));
            Ok(())
        }).unwrap();
    }

    #[test]
    fn rust_matcher_invokes_typed_callbacks_and_preserves_python_exception() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let callbacks = PyModule::from_code(py, c_str!(r#"
seen = []
atoms = []
bonds = []
marker = ValueError("original callback failure")
def final_filter(target, indices):
    seen.append((target, tuple(indices)))
    return True
def atom_filter(query, target):
    atoms.append((query.atomic_number(), target.atomic_number()))
    return True
def bond_filter(query, target):
    bonds.append((query.id(), target.id()))
    return True
def fail(target, indices):
    raise marker
def wrong_type(target, indices):
    return "not a bool"
"#), c_str!("callback_probe.py"), c_str!("callback_probe"))?;
            let target = Py::new(py, Molecule { inner: ck::Molecule::from_smiles("CCO").unwrap() })?;
            let query = ck::parse_smarts("CO").unwrap();
            callbacks.setattr("target", target.clone_ref(py))?;
            let mut params = SubstructMatchParams::from_native(py, ck::SubstructMatchParams::default())?;
            params.final_match = Some(callbacks.getattr("final_filter")?.extract()?);
            params.atom_match = Some(callbacks.getattr("atom_filter")?.extract()?);
            params.bond_match = Some(callbacks.getattr("bond_filter")?.extract()?);
            let run = |native: &ck::SubstructMatchParams| target.borrow(py).inner.substruct_matches_with_params(&query, native);
            let matches = params.with_callbacks(py, &target, run)?.unwrap();
            assert_eq!(matches.len(), 1);
            assert_eq!(matches[0].atom_mapping, vec![1, 2]);
            py.run(c_str!("assert seen[0][0] is not target\nassert seen[0][0].to_binary() == target.to_binary()\nassert seen[0][1] == (1, 2)\nassert atoms and bonds"), Some(&callbacks.dict()), None)?;

            params.final_match = Some(callbacks.getattr("fail")?.extract()?);
            let error = params.with_callbacks(py, &target, run).unwrap_err();
            assert!(error.value(py).is(&callbacks.getattr("marker")?));
            params.final_match = Some(callbacks.getattr("final_filter")?.extract()?);
            assert_eq!(params.with_callbacks(py, &target, run)?.unwrap().len(), 1);
            params.final_match = Some(callbacks.getattr("wrong_type")?.extract()?);
            assert!(params.with_callbacks(py, &target, run).unwrap_err().is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            Ok(())
        }).unwrap();
    }
}
