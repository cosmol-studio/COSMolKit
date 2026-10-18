//! Native SDF graph alternatives, with owned molecule values and query owners.
use crate::native::Molecule;
use crate::query_compiler::QueryGraph;
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native", frozen)]
pub enum SdfGraph {
    Molecule { value: Molecule },
    Query { value: QueryGraph },
}

impl SdfGraph {
    pub(crate) fn from_rust(inner: ck::SdfGraph) -> Self {
        match inner {
            ck::SdfGraph::Molecule(inner) => Self::Molecule {
                value: Molecule { inner },
            },
            ck::SdfGraph::Query(inner) => Self::Query {
                value: QueryGraph::from_rust(inner),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_graph_variants_keep_value_cow_and_query_owner_lifetime() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CCO").unwrap();
            let expected = native.clone();
            let original = Py::new(py, Molecule { inner: native })?;
            let class = py.get_type::<SdfGraph>();
            let graph = class.getattr("Molecule")?.call1((&original,))?;
            assert!(graph.is_instance(&class.getattr("Molecule")?)?);
            let value = graph.getattr("value")?.extract::<Py<Molecule>>()?;
            assert_eq!(value.borrow(py).inner, expected);
            assert!(std::ptr::eq(
                value.borrow(py).inner.topology(),
                expected.topology()
            ));
            original.borrow_mut(py).inner.add_hydrogens_().unwrap();
            assert_eq!(original.borrow(py).inner.num_atoms(), 9);
            assert_eq!(value.borrow(py).inner.num_atoms(), 3);
            let value_again = graph.getattr("value")?.extract::<Py<Molecule>>()?;
            assert_eq!(value_again.borrow(py).inner, expected);
            drop(original);
            drop(graph);
            assert_eq!(value.borrow(py).inner.num_atoms(), 3);

            let query = QueryGraph::from_rust(ck::parse_smarts("C=O").unwrap());
            let query = Py::new(py, query)?;
            let graph = class.getattr("Query")?.call1((&query,))?;
            assert!(graph.is_instance(&class.getattr("Query")?)?);
            let value = graph.getattr("value")?.extract::<Py<QueryGraph>>()?;
            assert!(std::ptr::eq(
                value.borrow(py).inner(),
                query.borrow(py).inner()
            ));
            drop(query);
            drop(graph);
            assert_eq!(
                value
                    .bind(py)
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                2
            );

            let native = ck::SdfGraph::Query(ck::parse_smarts("N").unwrap());
            let graph = SdfGraph::from_rust(native).into_pyobject(py)?.into_any();
            assert!(graph.is_instance(&class.getattr("Query")?)?);
            assert_eq!(
                graph
                    .getattr("value")?
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                1
            );
            Ok(())
        })
        .unwrap();
    }
}
