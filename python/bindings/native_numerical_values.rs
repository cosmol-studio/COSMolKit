//! Readonly transport of actual facade numerical results.
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Readonly dense square matrix with native row/column indexing.
    pub struct DenseMatrix {
        pub(crate) inner: ck::DenseMatrix,
    }
}

#[pymethods]
impl DenseMatrix {
    #[getter]
    /// Number of rows and columns in the square matrix.
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }
    #[getter]
    /// Independent flat list of matrix entries in native storage order.
    fn values(&self) -> Vec<f64> {
        self.inner.values().to_vec()
    }
    /// Return the matrix entry, or None when either index is out of range.
    fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.inner.get(row, column)
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Energy and per-atom Cartesian gradient evaluated by a persistent native force field.
    pub struct ForceFieldEnergyGradient {
        pub(crate) inner: ck::ForceFieldEnergyGradient,
    }
}

#[pymethods]
impl ForceFieldEnergyGradient {
    #[getter]
    /// Evaluated force-field energy.
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    #[getter]
    /// Independent list of Cartesian gradient triples in atom order.
    fn gradient(&self) -> Vec<[f64; 3]> {
        self.inner.gradient().to_vec()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Persistent force-field minimization outcome with convergence, iteration count and final energy.
    pub struct ForceFieldMinimizeOutcome {
        pub(crate) inner: ck::ForceFieldMinimizeOutcome,
    }
}

#[pymethods]
impl ForceFieldMinimizeOutcome {
    #[getter]
    /// Whether minimization reached the requested convergence criteria.
    fn converged(&self) -> bool {
        self.inner.converged()
    }
    #[getter]
    /// Number of minimization iterations performed.
    fn iterations(&self) -> u32 {
        self.inner.iterations()
    }
    #[getter]
    /// Evaluated force-field energy.
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_facade_numerical_results_preserve_bits_and_readonly_fields() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let matrix = molecule.distance_matrix().unwrap();
            let expected = matrix
                .values()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>();
            let dimension = matrix.dimension();
            let projected = Py::new(py, DenseMatrix { inner: matrix })?;
            let values = projected
                .bind(py)
                .getattr("values")?
                .extract::<Vec<f64>>()?;
            assert_eq!(
                values.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                expected
            );
            assert_eq!(
                projected
                    .bind(py)
                    .getattr("dimension")?
                    .extract::<usize>()?,
                dimension
            );
            assert_eq!(
                projected
                    .bind(py)
                    .call_method1("get", (dimension, 0))?
                    .extract::<Option<f64>>()?,
                None
            );

            let mut builder = molecule.to_builder();
            builder
                .add_3d_conformer(vec![[0., 0., 0.], [1.9, 0.2, 0.], [2.9, 0.4, 0.1]])
                .unwrap();
            let molecule = builder.build().unwrap().with_assigned_valence().unwrap();
            let mut field = molecule.uff_force_field().unwrap();
            let result = field.energy_gradient().unwrap();
            let energy = result.energy().to_bits();
            let gradient = result
                .gradient()
                .iter()
                .flatten()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>();
            let projected = Py::new(py, ForceFieldEnergyGradient { inner: result })?;
            assert_eq!(
                projected
                    .bind(py)
                    .getattr("energy")?
                    .extract::<f64>()?
                    .to_bits(),
                energy
            );
            let actual = projected
                .bind(py)
                .getattr("gradient")?
                .extract::<Vec<[f64; 3]>>()?;
            assert_eq!(
                actual
                    .iter()
                    .flatten()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                gradient
            );
            assert!(projected.bind(py).setattr("energy", 0.).is_err());

            let result = field
                .minimize_with_params_(&ck::ForceFieldMinimizeParams::new(2, 1e-4, 1e-6))
                .unwrap();
            let expected = (
                result.converged(),
                result.iterations(),
                result.energy().to_bits(),
            );
            let projected = Py::new(py, ForceFieldMinimizeOutcome { inner: result })?;
            assert_eq!(
                projected.bind(py).getattr("converged")?.extract::<bool>()?,
                expected.0
            );
            assert_eq!(
                projected.bind(py).getattr("iterations")?.extract::<u32>()?,
                expected.1
            );
            assert_eq!(
                projected
                    .bind(py)
                    .getattr("energy")?
                    .extract::<f64>()?
                    .to_bits(),
                expected.2
            );
            assert!(projected.bind(py).setattr("iterations", 4).is_err());
            Ok(())
        })
        .unwrap();
    }
}
