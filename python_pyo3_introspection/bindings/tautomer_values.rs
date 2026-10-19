use crate::leaf_values::TautomerEnumerationStatus;
use crate::native::Molecule;
use crate::property_values::decode_source_text;
use cosmolkit as ck;
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use std::collections::BTreeSet;
use std::sync::Arc;

#[pyclass(module = "cosmolkit", frozen)]
pub struct TautomerEnumeration {
    pub(crate) inner: Arc<ck::TautomerEnumeration>,
}
#[pymethods]
impl TautomerEnumeration {
    #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_tautomers=crate::arguments::Argument::Missing,max_transforms=crate::arguments::Argument::Missing,remove_sp3_stereo=crate::arguments::Argument::Missing,remove_bond_stereo=crate::arguments::Argument::Missing,remove_isotopic_hydrogens=crate::arguments::Argument::Missing,reassign_stereo=crate::arguments::Argument::Missing,callback=crate::arguments::Argument::Missing,scorer=crate::arguments::Argument::Missing,score_params=crate::arguments::Argument::Missing))]
    fn canonical_tautomer(
        &self,
        py: Python<'_>,
        params: crate::arguments::Argument<
            Option<Py<crate::tautomer_configuration::TautomerParams>>,
        >,
        max_tautomers: crate::arguments::Argument<Option<u32>>,
        max_transforms: crate::arguments::Argument<Option<u32>>,
        remove_sp3_stereo: crate::arguments::Argument<Option<bool>>,
        remove_bond_stereo: crate::arguments::Argument<Option<bool>>,
        remove_isotopic_hydrogens: crate::arguments::Argument<Option<bool>>,
        reassign_stereo: crate::arguments::Argument<Option<bool>>,
        callback: crate::arguments::Argument<
            Option<crate::tautomer_callbacks::EnumerationCallback>,
        >,
        scorer: crate::arguments::Argument<Option<crate::tautomer_callbacks::Scorer>>,
        score_params: crate::arguments::Argument<
            Option<Py<crate::tautomer_values::TautomerScoreParams>>,
        >,
    ) -> PyResult<Molecule> {
        let keyword_form = max_tautomers.is_supplied()
            || max_transforms.is_supplied()
            || remove_sp3_stereo.is_supplied()
            || remove_bond_stereo.is_supplied()
            || remove_isotopic_hydrogens.is_supplied()
            || reassign_stereo.is_supplied()
            || callback.is_supplied()
            || scorer.is_supplied()
            || score_params.is_supplied();
        crate::arguments::configuration_arguments!(params; max_tautomers, max_transforms, remove_sp3_stereo, remove_bond_stereo, remove_isotopic_hydrogens, reassign_stereo, callback, scorer, score_params);
        if keyword_form {
            self._canonical_tautomer_keywords(
                py,
                max_tautomers,
                max_transforms,
                remove_sp3_stereo,
                remove_bond_stereo,
                remove_isotopic_hydrogens,
                reassign_stereo,
                callback,
                scorer,
                score_params,
            )
        } else {
            self._canonical_tautomer_params(py, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"TautomerEnumeration.canonical_tautomer"}"#]
    #[pyo3(signature=(params=None))]
    fn _canonical_tautomer_params(
        &self,
        py: Python<'_>,
        params: Option<Py<crate::tautomer_configuration::TautomerParams>>,
    ) -> PyResult<Molecule> {
        let params = params
            .as_ref()
            .map(|value| value.try_borrow(py)?.capture(py))
            .transpose()?;
        crate::tautomer_producers::canonical_enumeration(py, self, params)
    }
    #[doc = r#"{"role":"overload","target":"TautomerEnumeration.canonical_tautomer"}"#]
    #[pyo3(signature=(*,max_tautomers=None,max_transforms=None,remove_sp3_stereo=None,remove_bond_stereo=None,remove_isotopic_hydrogens=None,reassign_stereo=None,callback=None,scorer=None,score_params=None))]
    fn _canonical_tautomer_keywords(
        &self,
        py: Python<'_>,
        max_tautomers: Option<u32>,
        max_transforms: Option<u32>,
        remove_sp3_stereo: Option<bool>,
        remove_bond_stereo: Option<bool>,
        remove_isotopic_hydrogens: Option<bool>,
        reassign_stereo: Option<bool>,
        callback: Option<crate::tautomer_callbacks::EnumerationCallback>,
        scorer: Option<crate::tautomer_callbacks::Scorer>,
        score_params: Option<Py<crate::tautomer_values::TautomerScoreParams>>,
    ) -> PyResult<Molecule> {
        if max_tautomers.is_none()
            && max_transforms.is_none()
            && remove_sp3_stereo.is_none()
            && remove_bond_stereo.is_none()
            && remove_isotopic_hydrogens.is_none()
            && reassign_stereo.is_none()
            && callback.is_none()
            && scorer.is_none()
            && score_params.is_none()
        {
            return self._canonical_tautomer_params(py, None);
        }
        let params = crate::tautomer_producers::keywords(
            py,
            max_tautomers,
            max_transforms,
            remove_sp3_stereo,
            remove_bond_stereo,
            remove_isotopic_hydrogens,
            reassign_stereo,
            callback,
            scorer,
            score_params,
        )?;
        crate::tautomer_producers::canonical_enumeration(py, self, Some(params))
    }
    fn canonical_tautomer_with_params(
        &self,
        py: Python<'_>,
        params: Py<crate::tautomer_configuration::TautomerParams>,
    ) -> PyResult<Molecule> {
        self._canonical_tautomer_params(py, Some(params))
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn status(&self) -> TautomerEnumerationStatus {
        self.inner.status().into()
    }
    fn modified_atoms(&self) -> BTreeSet<usize> {
        self.inner
            .modified_atoms()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn modified_bonds(&self) -> BTreeSet<usize> {
        self.inner
            .modified_bonds()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn canonical_smiles(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        self.inner
            .canonical_smiles()
            .into_iter()
            .map(|text| decode_source_text(py, text))
            .collect()
    }
    fn get(&self, index: usize) -> Option<Molecule> {
        self.inner.get(index).map(|inner| Molecule {
            inner: inner.clone(),
        })
    }
    #[doc = r#"{"role":"return_type","imports":["collections.abc"],"type":"collections.abc.Iterator[Molecule]"}"#]
    fn iter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        Py::new(
            py,
            TautomerMoleculeIterator {
                owner: Arc::clone(&self.inner),
                index: 0,
            },
        )?
        .into_bound(py)
        .into_any()
        .cast_into::<PyIterator>()
        .map_err(Into::into)
    }
    #[doc = r#"{"role":"return_type","imports":["collections.abc"],"type":"collections.abc.Iterator[tuple[str, Molecule]]"}"#]
    fn entries<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        Py::new(
            py,
            TautomerEntryIterator {
                owner: Arc::clone(&self.inner),
                index: 0,
            },
        )?
        .into_bound(py)
        .into_any()
        .cast_into::<PyIterator>()
        .map_err(Into::into)
    }
}

#[pyclass(module = "cosmolkit")]
struct TautomerMoleculeIterator {
    owner: Arc<ck::TautomerEnumeration>,
    index: usize,
}

#[pymethods]
impl TautomerMoleculeIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self) -> Option<Molecule> {
        let inner = self.owner.get(self.index)?;
        self.index += 1;
        Some(Molecule {
            inner: inner.clone(),
        })
    }
}

#[pyclass(module = "cosmolkit")]
struct TautomerEntryIterator {
    owner: Arc<ck::TautomerEnumeration>,
    index: usize,
}

#[pymethods]
impl TautomerEntryIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<(String, Molecule)>> {
        // The native iterator is Map<slice::Iter<_>, _>: nth indexes in O(1),
        // without copying the remaining keys or candidates into Python values.
        let Some((key, molecule)) = self.owner.entries().nth(self.index) else {
            return Ok(None);
        };
        self.index += 1;
        Ok(Some((
            decode_source_text(py, key)?,
            Molecule {
                inner: molecule.clone(),
            },
        )))
    }
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct TautomerScore {
    pub(crate) inner: ck::TautomerScore,
}

#[pymethods]
impl TautomerScore {
    fn ring(&self) -> i32 {
        self.inner.ring()
    }
    fn substructure(&self) -> i32 {
        self.inner.substructure()
    }
    fn hetero_hydrogen(&self) -> i32 {
        self.inner.hetero_hydrogen()
    }
    fn total(&self) -> i32 {
        self.inner.total()
    }
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct TautomerScoreTerm {
    pub(crate) inner: ck::TautomerScoreTerm,
}

#[pymethods]
impl TautomerScoreTerm {
    #[new]
    fn construct(name: String, smarts: String, score: i32) -> Self {
        Self {
            inner: ck::TautomerScoreTerm::new(name, smarts, score),
        }
    }
    #[staticmethod]
    fn new(name: String, smarts: String, score: i32) -> Self {
        Self::construct(name, smarts, score)
    }
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn smarts(&self) -> &str {
        self.inner.smarts()
    }
    fn score(&self) -> i32 {
        self.inner.score()
    }
}

#[pyfunction]
pub fn default_tautomer_score_terms() -> Vec<TautomerScoreTerm> {
    ck::default_tautomer_score_terms()
        .iter()
        .cloned()
        .map(|inner| TautomerScoreTerm { inner })
        .collect()
}

impl crate::parameters::Configuration for TautomerScoreTerm {
    type Native = ck::TautomerScoreTerm;
    fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
        Ok(Self { inner })
    }
    fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
        Ok(self.inner.clone())
    }
}
ck::project_configuration!(ck::TautomerScoreParams);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entries_keep_native_order_shared_owner_and_independent_cow_values() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CC(=O)C").unwrap();
            let enumeration = Arc::new(source.enumerate_tautomers().unwrap());
            let expected = enumeration
                .entries()
                .map(|(key, molecule)| {
                    Ok((decode_source_text(py, key)?, molecule.to_binary().unwrap()))
                })
                .collect::<PyResult<Vec<_>>>()?;
            assert!(!expected.is_empty());
            let projected = Py::new(
                py,
                TautomerEnumeration {
                    inner: Arc::clone(&enumeration),
                },
            )?;
            let iterator = projected.bind(py).call_method0("entries")?;
            assert!(iterator.call_method0("__iter__")?.is(&iterator));
            drop(projected);
            drop(source);
            for ((key, binary), (_, native)) in expected.iter().zip(enumeration.entries()) {
                let (actual_key, actual) =
                    iterator
                        .call_method0("__next__")?
                        .extract::<(String, Py<Molecule>)>()?;
                assert_eq!(&actual_key, key);
                assert!(std::ptr::eq(
                    actual.borrow(py).inner.topology(),
                    native.topology()
                ));
                assert_eq!(actual.borrow(py).inner.to_binary().unwrap(), *binary);
                actual.borrow_mut(py).inner.add_hydrogens_().unwrap();
                assert!(actual.borrow(py).inner.num_atoms() > native.num_atoms());
                assert_eq!(native.to_binary().unwrap(), *binary);
            }
            drop(enumeration);
            for _ in 0..2 {
                assert!(
                    iterator
                        .call_method0("__next__")
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyStopIteration>(py)
                );
            }
            let empty = Py::new(
                py,
                TautomerEnumeration {
                    inner: Arc::new(ck::TautomerEnumeration::default()),
                },
            )?;
            assert!(
                empty
                    .bind(py)
                    .call_method0("entries")?
                    .call_method0("__next__")
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyStopIteration>(py)
            );
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn actual_enumeration_queries_preserve_keys_sets_and_cow_values() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CC(=O)C")
                .unwrap()
                .enumerate_tautomers()
                .unwrap();
            let expected_len = source.len();
            let expected_status = source.status();
            let expected_keys = source
                .canonical_smiles()
                .into_iter()
                .map(|key| decode_source_text(py, key))
                .collect::<PyResult<Vec<_>>>()?;
            let expected_atoms = source
                .modified_atoms()
                .iter()
                .map(|id| id.index())
                .collect::<BTreeSet<_>>();
            let expected_bonds = source
                .modified_bonds()
                .iter()
                .map(|id| id.index())
                .collect::<BTreeSet<_>>();
            let projected = TautomerEnumeration {
                inner: Arc::new(source),
            };
            assert_eq!(projected.len(), expected_len);
            assert!(!projected.is_empty());
            assert!(projected.status() == expected_status.into());
            assert_eq!(projected.canonical_smiles(py)?, expected_keys);
            assert_eq!(projected.modified_atoms(), expected_atoms);
            assert_eq!(projected.modified_bonds(), expected_bonds);
            assert!(projected.get(expected_len).is_none());
            let candidate = projected.get(0).unwrap();
            assert!(std::ptr::eq(
                candidate.inner.topology(),
                projected.inner.get(0).unwrap().topology()
            ));
            let projected = Py::new(py, projected)?;
            assert_eq!(
                projected
                    .bind(py)
                    .call_method0("canonical_smiles")?
                    .extract::<Vec<String>>()?,
                expected_keys
            );
            assert_eq!(
                projected
                    .bind(py)
                    .call_method0("modified_atoms")?
                    .extract::<BTreeSet<usize>>()?,
                expected_atoms
            );
            assert_eq!(
                projected
                    .bind(py)
                    .call_method0("modified_bonds")?
                    .extract::<BTreeSet<usize>>()?,
                expected_bonds
            );
            let value = projected.bind(py).call_method1("get", (0,))?;
            assert_eq!(
                value.call_method0("num_atoms")?.extract::<usize>()?,
                candidate.inner.num_atoms()
            );
            assert!(
                projected
                    .bind(py)
                    .call_method1("get", (expected_len,))?
                    .is_none()
            );
            assert!(
                projected
                    .bind(py)
                    .call_method1("get", (-1,))
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyOverflowError>(py)
            );
            let iterator = projected.bind(py).call_method0("iter")?;
            assert!(iterator.call_method0("__iter__")?.is(&iterator));
            let native_owner = Arc::clone(&projected.borrow(py).inner);
            drop(projected);
            assert_eq!(candidate.inner.num_atoms(), 4);
            for expected in native_owner.iter() {
                let next = iterator.call_method0("__next__")?;
                let next = next.extract::<PyRef<'_, Molecule>>()?;
                assert!(std::ptr::eq(next.inner.topology(), expected.topology()));
            }
            drop(native_owner);
            for _ in 0..2 {
                assert!(
                    iterator
                        .call_method0("__next__")
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyStopIteration>(py)
                );
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn actual_native_tautomer_score_getters_remain_exact() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CC(=O)C").unwrap();
            let source = molecule.tautomer_score().unwrap();
            let expected = [
                source.ring(),
                source.substructure(),
                source.hetero_hydrogen(),
                source.total(),
            ];
            let result = Py::new(py, TautomerScore { inner: source })?;
            for (method, value) in ["ring", "substructure", "hetero_hydrogen", "total"]
                .into_iter()
                .zip(expected)
            {
                assert_eq!(
                    result.bind(py).call_method0(method)?.extract::<i32>()?,
                    value
                );
            }
            Ok(())
        })
        .unwrap();
    }
}
