//! Source-owned Morgan selectors and configuration, without fingerprint logic.
use crate::fingerprint_providers::AtomPairAtomInvariantsGenerator;
use crate::query_compiler::QueryGraph;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
    #[derive(Clone, PartialEq)]
/// Morgan atom-invariant selection: connectivity, standard features, or supplied feature query patterns.
pub struct MorganInvariants {
        pub(crate) inner: ck::MorganInvariants,
    }
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
    /// Use the supplied query patterns to define atom-feature invariants.
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

#[pyclass(module = "cosmolkit", frozen)]
/// Reusable atom-invariant provider for Morgan fingerprints.
pub struct MorganAtomInvariantsGenerator {
    pub(crate) inner: ck::MorganAtomInvariantsGenerator,
}

#[pymethods]
impl MorganAtomInvariantsGenerator {
    #[staticmethod]
    #[pyo3(signature=(include_ring_membership=true))]
    /// Create a connectivity-invariant provider, optionally including ring membership.
    fn connectivity(include_ring_membership: bool) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::connectivity(include_ring_membership),
        }
    }
    #[staticmethod]
    #[pyo3(signature=(patterns=None))]
    /// Create a feature-invariant provider; None selects the native feature patterns.
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
    /// Create a Morgan atom provider from an independent copy of the atom-pair invariant configuration.
    fn atom_pair(py: Python<'_>, generator: &AtomPairAtomInvariantsGenerator) -> PyResult<Self> {
        Ok(Self {
            inner: ck::MorganAtomInvariantsGenerator::atom_pair(generator.to_rust(py)?),
        })
    }
}

ck::project_configuration!(ck::MorganFingerprintParams);

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
            let mut native = ck::MorganFingerprintParams::default();
            native.from_atoms = Some(vec![]);
            native.custom_atom_invariants = Some(vec![u32::MAX]);
            native.invariants = invariants.inner.clone();
            let params = MorganFingerprintParams::from_rust(py, native)?;
            let params = params.to_rust(py)?;
            assert_eq!(params.from_atoms, Some(vec![]));
            assert_eq!(params.ignore_atoms, None);
            assert_eq!(params.custom_atom_invariants, Some(vec![u32::MAX]));
            assert_eq!(params.invariants, invariants.inner);
            Ok(())
        })
        .unwrap();
    }
}
