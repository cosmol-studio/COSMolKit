//! Immutable topology view retains the native molecule's shared blocks.
use crate::atom_bond_values::{Atom, Bond};
use crate::group_values::SubstanceGroup;
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "cosmolkit", frozen)]
/// Read-only topology view retaining the molecule's COW storage. Requested row lists are detached copies.
pub struct TopologyBlock {
    owner: ck::Molecule,
}

impl TopologyBlock {
    pub(crate) fn from_molecule(owner: &ck::Molecule) -> Self {
        // Molecule::clone retains its COW blocks. No topology rows are copied
        // until a caller explicitly requests owned Python row values.
        Self {
            owner: owner.clone(),
        }
    }
    pub(crate) fn inner(&self) -> &ck::TopologyBlock {
        self.owner.topology()
    }
}

#[pymethods]
impl TopologyBlock {
    #[getter]
    /// Independent detached atom rows in stored order.
    fn atoms(&self) -> Vec<Atom> {
        self.inner()
            .atoms
            .iter()
            .cloned()
            .map(|inner| Atom { inner })
            .collect()
    }
    #[getter]
    /// Independent detached bond rows in stored order.
    fn bonds(&self) -> Vec<Bond> {
        self.inner()
            .bonds
            .iter()
            .cloned()
            .map(|inner| Bond { inner })
            .collect()
    }
    #[getter]
    /// Independent copies of the stored substance groups.
    fn substance_groups(&self) -> Vec<SubstanceGroup> {
        self.inner()
            .substance_groups
            .iter()
            .cloned()
            .map(|inner| SubstanceGroup { inner })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn topology_result_retains_native_shared_block_after_source_is_dropped() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let result = TopologyBlock::from_molecule(&molecule);
            assert!(std::ptr::eq(result.inner(), molecule.topology()));
            assert_eq!(result.inner(), molecule.topology());
            let expected_ids = molecule
                .topology()
                .atoms
                .iter()
                .map(|atom| atom.id().index())
                .collect::<Vec<_>>();
            drop(molecule);
            let result = Py::new(py, result)?;
            let atoms = result
                .bind(py)
                .getattr("atoms")?
                .extract::<Vec<Py<Atom>>>()?;
            let actual_ids = atoms
                .iter()
                .map(|atom| {
                    atom.bind(py)
                        .getattr("id")
                        .and_then(|value| value.extract::<usize>())
                })
                .collect::<PyResult<Vec<_>>>()?;
            assert_eq!(actual_ids, expected_ids);
            assert!(result.bind(py).setattr("atoms", py.None()).is_err());
            Ok(())
        })
        .unwrap();
    }
}
