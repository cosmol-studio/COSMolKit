//! Detached source result transport. No energy or descriptor calculation lives here.
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    pub struct OptimizationOutcome { pub(crate) inner: ck::OptimizationOutcome }
}
#[pymethods]
impl OptimizationOutcome {
    #[staticmethod]
    fn converged_with_energy(energy: f64) -> Self {
        Self {
            inner: ck::OptimizationOutcome::converged_with_energy(energy),
        }
    }
    #[staticmethod]
    fn not_converged_with_energy(energy: f64) -> Self {
        Self {
            inner: ck::OptimizationOutcome::not_converged_with_energy(energy),
        }
    }
    #[staticmethod]
    fn missing_parameters() -> Self {
        Self {
            inner: ck::OptimizationOutcome::missing_parameters(),
        }
    }
    #[getter]
    fn status(&self) -> crate::leaf_values::OptimizationStatus {
        self.inner.status().into()
    }
    #[getter]
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
    #[getter]
    fn converged(&self) -> bool {
        self.inner.converged()
    }
    #[getter]
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
}

macro_rules! gradient_result {
    ($name:ident) => {
        crate::configuration::type_projection! {
            #[pyclass(module = "cosmolkit", frozen, from_py_object)]
            #[derive(Clone)]
            pub struct $name {
                pub(crate) inner: ck::$name,
            }
        }
        #[pymethods]
        impl $name {
            #[new]
            fn new(energy: f64, gradient: Vec<f64>) -> Self {
                Self {
                    inner: ck::$name { energy, gradient },
                }
            }
            #[getter]
            fn energy(&self) -> f64 {
                self.inner.energy()
            }
            #[getter]
            fn gradient(&self) -> Vec<f64> {
                self.inner.gradient().to_vec()
            }
        }
    };
}
gradient_result!(UffEnergyGradient);
gradient_result!(MmffEnergyGradient);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct UffOptimizationResult {
        pub(crate) inner: ck::UffOptimizationResult,
    }
}
#[pymethods]
impl UffOptimizationResult {
    #[getter]
    fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome {
            inner: self.inner.outcome(),
        }
    }
    #[new]
    fn new(molecule: &Molecule, outcome: &OptimizationOutcome) -> Self {
        Self {
            inner: ck::UffOptimizationResult {
                molecule: molecule.inner.clone(),
                outcome: outcome.inner,
            },
        }
    }
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct UffConformerResult {
        pub(crate) inner: ck::UffConformerResult,
    }
}
#[pymethods]
impl UffConformerResult {
    #[getter]
    fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome {
            inner: self.inner.outcome(),
        }
    }
    #[new]
    fn new(conformer_id: usize, outcome: &OptimizationOutcome) -> Self {
        Self {
            inner: ck::UffConformerResult {
                conformer_id,
                outcome: outcome.inner,
            },
        }
    }
    #[getter]
    fn conformer_id(&self) -> usize {
        self.inner.conformer_id()
    }
    #[getter]
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct UffConformerOptimizationResult {
        pub(crate) inner: ck::UffConformerOptimizationResult,
    }
}
#[pymethods]
impl UffConformerOptimizationResult {
    #[new]
    fn new(molecule: &Molecule, conformers: Vec<UffConformerResult>) -> Self {
        Self {
            inner: ck::UffConformerOptimizationResult {
                molecule: molecule.inner.clone(),
                conformers: conformers.into_iter().map(|value| value.inner).collect(),
            },
        }
    }
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    fn conformer_results(&self) -> Vec<UffConformerResult> {
        self.inner
            .conformer_results()
            .iter()
            .cloned()
            .map(|inner| UffConformerResult { inner })
            .collect()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct MmffOptimizeMoleculeResult {
        pub(crate) inner: ck::MmffOptimizeMoleculeResult,
    }
}
#[pymethods]
impl MmffOptimizeMoleculeResult {
    #[getter]
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
    #[getter]
    fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome {
            inner: self.inner.outcome(),
        }
    }
    #[new]
    fn new(molecule: &Molecule, outcome: &OptimizationOutcome) -> Self {
        Self {
            inner: ck::MmffOptimizeMoleculeResult {
                molecule: molecule.inner.clone(),
                outcome: outcome.inner,
            },
        }
    }
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct MmffOptimizeMoleculeConfResult {
        pub(crate) inner: ck::MmffOptimizeMoleculeConfResult,
    }
}
#[pymethods]
impl MmffOptimizeMoleculeConfResult {
    #[getter]
    fn outcome(&self) -> OptimizationOutcome {
        OptimizationOutcome {
            inner: self.inner.outcome(),
        }
    }
    #[new]
    fn new(outcome: &OptimizationOutcome) -> Self {
        Self {
            inner: ck::MmffOptimizeMoleculeConfResult {
                outcome: outcome.inner,
            },
        }
    }
    #[getter]
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct MmffOptimizeMoleculeConfsResult {
        pub(crate) inner: ck::MmffOptimizeMoleculeConfsResult,
    }
}
#[pymethods]
impl MmffOptimizeMoleculeConfsResult {
    #[new]
    fn new(molecule: &Molecule, conformer_results: Vec<MmffOptimizeMoleculeConfResult>) -> Self {
        Self {
            inner: ck::MmffOptimizeMoleculeConfsResult {
                molecule: molecule.inner.clone(),
                conformer_results: conformer_results
                    .into_iter()
                    .map(|value| value.inner)
                    .collect(),
            },
        }
    }
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    fn conformer_results(&self) -> Vec<MmffOptimizeMoleculeConfResult> {
        self.inner
            .conformer_results()
            .iter()
            .cloned()
            .map(|inner| MmffOptimizeMoleculeConfResult { inner })
            .collect()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct CrippenTotals {
        pub(crate) inner: ck::CrippenTotals,
    }
}
#[pymethods]
impl CrippenTotals {
    #[new]
    fn new(logp: f64, molar_refractivity: f64) -> Self {
        Self {
            inner: ck::CrippenTotals {
                logp,
                molar_refractivity,
            },
        }
    }
    #[getter]
    fn logp(&self) -> f64 {
        self.inner.logp
    }
    #[getter]
    fn molar_refractivity(&self) -> f64 {
        self.inner.molar_refractivity
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct LabuteAsaContributions {
        pub(crate) inner: ck::LabuteAsaContributions,
    }
}
#[pymethods]
impl LabuteAsaContributions {
    #[new]
    fn new(asa: f64, atom_contributions: Vec<f64>, hydrogen_contribution: f64) -> Self {
        Self {
            inner: ck::LabuteAsaContributions {
                asa,
                atom_contributions,
                hydrogen_contribution,
            },
        }
    }
    #[getter]
    fn asa(&self) -> f64 {
        self.inner.asa
    }
    #[getter]
    fn atom_contributions(&self) -> Vec<f64> {
        self.inner.atom_contributions.clone()
    }
    #[getter]
    fn hydrogen_contribution(&self) -> f64 {
        self.inner.hydrogen_contribution
    }
}
