//! Actual native MCS result transport with shared query ownership.
use crate::property_values::decode_source_text;
use crate::query_compiler::QueryGraph;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::sync::Arc;

#[pyclass(module = "cosmolkit", frozen)]
/// Native maximum-common-substructure result, completion state and retained query graphs.
pub struct McsResult {
    pub(crate) inner: Arc<ck::McsResult>,
}

#[pymethods]
impl McsResult {
    #[getter]
    /// Maximum-common-substructure query graph, or None if no query was produced.
    fn query(&self) -> Option<QueryGraph> {
        QueryGraph::from_mcs_query(Arc::clone(&self.inner))
    }
    #[getter]
    /// Number of atoms in the reported common substructure.
    fn atom_count(&self) -> usize {
        self.inner.atom_count
    }
    #[getter]
    /// Number of bonds in the reported common substructure.
    fn bond_count(&self) -> usize {
        self.inner.bond_count
    }
    #[getter]
    /// Whether the native search completed rather than stopping early.
    fn completed(&self) -> bool {
        self.inner.completed
    }
    #[getter]
    /// SMARTS representation of the reported common substructure.
    fn smarts(&self, py: Python<'_>) -> PyResult<String> {
        decode_source_text(py, &self.inner.smarts)
    }
    #[getter]
    /// Alternative maximum-common-substructure query graphs keyed by their SMARTS.
    fn degenerate(&self, py: Python<'_>) -> PyResult<BTreeMap<String, QueryGraph>> {
        self.inner
            .degenerate
            .keys()
            .map(|key| {
                Ok((
                    decode_source_text(py, key)?,
                    QueryGraph::from_mcs_degenerate(Arc::clone(&self.inner), key.clone()),
                ))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_mcs_result_fields_and_query_views_preserve_native_owner() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let first = ck::Molecule::from_smiles("CCO").unwrap();
            let second = ck::Molecule::from_smiles("CCN").unwrap();
            let result = ck::maximum_common_substructure_with_params(
                &[&first, &second],
                &ck::McsParameters {
                    store_all: true,
                    ..Default::default()
                },
            )
            .unwrap();
            let expected_atoms = result.atom_count;
            let expected_bonds = result.bond_count;
            let expected_completed = result.completed;
            let expected_smarts = decode_source_text(py, &result.smarts)?;
            let expected_degenerate = result.degenerate.len();
            let projected = McsResult {
                inner: Arc::new(result),
            };
            let query = projected.query();
            if let Some(query) = &query {
                assert!(std::ptr::eq(
                    query.inner(),
                    projected.inner.query.as_ref().unwrap()
                ));
            }
            let degenerate = projected.degenerate(py)?;
            assert_eq!(degenerate.len(), expected_degenerate);
            for (key, value) in &degenerate {
                assert!(std::ptr::eq(
                    value.inner(),
                    &projected.inner.degenerate[key.as_bytes()]
                ));
            }
            let projected = Py::new(py, projected)?;
            assert_eq!(
                projected
                    .bind(py)
                    .getattr("atom_count")?
                    .extract::<usize>()?,
                expected_atoms
            );
            assert_eq!(
                projected
                    .bind(py)
                    .getattr("bond_count")?
                    .extract::<usize>()?,
                expected_bonds
            );
            assert_eq!(
                projected.bind(py).getattr("completed")?.extract::<bool>()?,
                expected_completed
            );
            assert_eq!(
                projected.bind(py).getattr("smarts")?.extract::<String>()?,
                expected_smarts
            );
            assert!(projected.bind(py).setattr("completed", false).is_err());
            drop(projected);
            if let Some(query) = query {
                let query = Py::new(py, query)?;
                assert_eq!(
                    query.bind(py).getattr("num_atoms")?.extract::<usize>()?,
                    expected_atoms
                );
            }
            for query in degenerate.into_values() {
                let query = Py::new(py, query)?;
                assert_eq!(
                    query.bind(py).getattr("num_bonds")?.extract::<usize>()?,
                    expected_bonds
                );
            }
            Ok(())
        })
        .unwrap();
    }
}
