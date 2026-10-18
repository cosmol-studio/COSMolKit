//! Source-owned Morgan selectors and configuration, without fingerprint logic.
use crate::fingerprint_configuration::MorganParams;
use crate::fingerprint_providers::AtomPairAtomInvariantsGenerator;
use crate::query_compiler::QueryGraph;
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native", frozen, from_py_object, eq)]
#[derive(Clone, PartialEq)]
pub struct MorganInvariants {
    pub(crate) inner: ck::MorganInvariants,
}

#[pymethods]
impl MorganInvariants {
    #[classattr]
    #[pyo3(name = "Connectivity")]
    const CONNECTIVITY: Self = Self {
        inner: ck::MorganInvariants::Connectivity,
    };
    #[classattr]
    #[pyo3(name = "Features")]
    const FEATURES: Self = Self {
        inner: ck::MorganInvariants::Features,
    };
    #[staticmethod]
    #[pyo3(name = "FeaturePatterns")]
    fn feature_patterns(py: Python<'_>, patterns: Vec<Py<QueryGraph>>) -> Self {
        Self {
            inner: ck::MorganInvariants::FeaturePatterns(
                patterns
                    .iter()
                    .map(|value| value.borrow(py).inner().clone())
                    .collect(),
            ),
        }
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct MorganAtomInvariantsGenerator {
    pub(crate) inner: ck::MorganAtomInvariantsGenerator,
}

#[pymethods]
impl MorganAtomInvariantsGenerator {
    #[staticmethod]
    #[pyo3(signature=(include_ring_membership=true))]
    fn connectivity(include_ring_membership: bool) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::connectivity(include_ring_membership),
        }
    }
    #[staticmethod]
    #[pyo3(signature=(patterns=None))]
    fn features(py: Python<'_>, patterns: Option<Vec<Py<QueryGraph>>>) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::features(patterns.map(|values| {
                values
                    .iter()
                    .map(|value| value.borrow(py).inner().clone())
                    .collect()
            })),
        }
    }
    #[staticmethod]
    fn atom_pair(generator: &AtomPairAtomInvariantsGenerator) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::atom_pair(generator.inner.clone()),
        }
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct MorganFingerprintParams {
    pub(crate) inner: ck::MorganFingerprintParams,
}

macro_rules! morgan_fingerprint_configuration {
    ($($field:ident),+ $(,)?) => {
        paste::paste! {
            #[pymethods]
            impl MorganFingerprintParams {
                #[new]
                #[pyo3(signature=(*, generator=None, $($field=None,)+ conformer_id=ck::MorganFingerprintParams::default().conformer_id, invariants=None))]
                fn new(generator: Option<&MorganParams>, $($field: Option<Vec<u32>>,)+ conformer_id: i32, invariants: Option<&MorganInvariants>) -> Self {
                    Self { inner: ck::MorganFingerprintParams {
                        generator: generator.map(|value| value.inner.clone()).unwrap_or_else(|| ck::MorganFingerprintParams::default().generator),
                        $($field,)+ conformer_id,
                        invariants: invariants.map(|value| value.inner.clone()).unwrap_or_else(|| ck::MorganFingerprintParams::default().invariants),
                    } }
                }
                #[getter]
                fn generator(&self) -> MorganParams { MorganParams { inner: self.inner.generator.clone() } }
                #[setter]
                fn set_generator(&mut self, value: Option<&MorganParams>) {
                    self.inner.generator = value.map(|value| value.inner.clone()).unwrap_or_else(|| ck::MorganFingerprintParams::default().generator);
                }
                $(#[getter]
                fn $field(&self) -> Option<Vec<u32>> { self.inner.$field.clone() }
                #[setter]
                fn [<set_ $field>](&mut self, value: Option<Vec<u32>>) { self.inner.$field = value; })+
                #[getter]
                fn conformer_id(&self) -> i32 { self.inner.conformer_id }
                #[setter]
                fn set_conformer_id(&mut self, value: i32) { self.inner.conformer_id = value; }
                #[getter]
                fn invariants(&self) -> MorganInvariants { MorganInvariants { inner: self.inner.invariants.clone() } }
                #[setter]
                fn set_invariants(&mut self, value: Option<&MorganInvariants>) {
                    self.inner.invariants = value.map(|value| value.inner.clone()).unwrap_or_else(|| ck::MorganFingerprintParams::default().invariants);
                }
            }
        }
    };
}
morgan_fingerprint_configuration!(
    from_atoms,
    ignore_atoms,
    custom_atom_invariants,
    custom_bond_invariants
);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_providers_capture_exact_graphs_and_configuration() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for patterns in [
                None,
                Some(vec![]),
                Some(vec![ck::parse_smarts("[#6]").unwrap()]),
            ] {
                let expected = ck::MorganAtomInvariantsGenerator::features(patterns.clone());
                let values = patterns
                    .clone()
                    .map(|values| {
                        values
                            .into_iter()
                            .map(|value| Py::new(py, QueryGraph::from_rust(value)))
                            .collect::<PyResult<Vec<_>>>()
                    })
                    .transpose()?;
                let actual = MorganAtomInvariantsGenerator::features(py, values);
                let expected =
                    ck::MorganFingerprintGenerator::new(None, Some(&expected), None).unwrap();
                let actual =
                    ck::MorganFingerprintGenerator::new(None, Some(&actual.inner), None).unwrap();
                assert_eq!(actual.to_json().unwrap(), expected.to_json().unwrap());
            }
            for include_rings in [false, true] {
                let actual = MorganAtomInvariantsGenerator::connectivity(include_rings);
                let expected = ck::MorganAtomInvariantsGenerator::connectivity(include_rings);
                let actual =
                    ck::MorganFingerprintGenerator::new(None, Some(&actual.inner), None).unwrap();
                let expected =
                    ck::MorganFingerprintGenerator::new(None, Some(&expected), None).unwrap();
                assert_eq!(actual.to_json().unwrap(), expected.to_json().unwrap());
            }
            let invariants = MorganInvariants::feature_patterns(
                py,
                vec![Py::new(
                    py,
                    QueryGraph::from_rust(ck::parse_smarts("O").unwrap()),
                )?],
            );
            let params = MorganFingerprintParams::new(
                None,
                Some(vec![]),
                None,
                Some(vec![u32::MAX]),
                None,
                -1,
                Some(&invariants),
            );
            assert_eq!(params.inner.from_atoms, Some(vec![]));
            assert_eq!(params.inner.ignore_atoms, None);
            assert_eq!(params.inner.custom_atom_invariants, Some(vec![u32::MAX]));
            assert_eq!(params.inner.invariants, invariants.inner);
            Ok(())
        })
        .unwrap();
    }
}
