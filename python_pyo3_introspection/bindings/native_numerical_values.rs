//! Readonly transport of actual facade numerical results.
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "cosmolkit", frozen)]
pub struct DenseMatrix {
    pub(crate) inner: ck::DenseMatrix,
}

#[pymethods]
impl DenseMatrix {
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }
    fn values(&self) -> Vec<f64> {
        self.inner.values().to_vec()
    }
    fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.inner.get(row, column)
    }
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct ForceFieldEnergyGradient {
    pub(crate) inner: ck::ForceFieldEnergyGradient,
}

#[pymethods]
impl ForceFieldEnergyGradient {
    #[getter]
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
    #[getter]
    fn gradient(&self) -> Vec<[f64; 3]> {
        self.inner.gradient().to_vec()
    }
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct ForceFieldMinimizeOutcome {
    pub(crate) inner: ck::ForceFieldMinimizeOutcome,
}

#[pymethods]
impl ForceFieldMinimizeOutcome {
    #[getter]
    fn converged(&self) -> bool {
        self.inner.converged()
    }
    #[getter]
    fn iterations(&self) -> u32 {
        self.inner.iterations()
    }
    #[getter]
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
                .call_method0("values")?
                .extract::<Vec<f64>>()?;
            assert_eq!(
                values.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                expected
            );
            assert_eq!(
                projected
                    .bind(py)
                    .call_method0("dimension")?
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
