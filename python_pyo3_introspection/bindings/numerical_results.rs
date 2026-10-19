//! Detached source result transport. No energy or descriptor calculation lives here.
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! gradient_result {
    ($name:ident) => {
        #[pyclass(module = "cosmolkit", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $name {
            pub(crate) inner: ck::$name,
        }
        #[pymethods]
        impl $name {
            #[new]
            fn new(energy: f64, gradient: Vec<f64>) -> Self {
                Self {
                    inner: ck::$name { energy, gradient },
                }
            }
            fn energy(&self) -> f64 {
                self.inner.energy()
            }
            fn gradient(&self) -> Vec<f64> {
                self.inner.gradient().to_vec()
            }
        }
    };
}
gradient_result!(UffEnergyGradient);
gradient_result!(MmffEnergyGradient);

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct UffOptimizationResult {
    pub(crate) inner: ck::UffOptimizationResult,
}
#[pymethods]
impl UffOptimizationResult {
    #[new]
    fn new(molecule: &Molecule, status: i32, energy: f64) -> Self {
        Self {
            inner: ck::UffOptimizationResult {
                molecule: molecule.inner.clone(),
                status,
                energy,
            },
        }
    }
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct UffConformerResult {
    pub(crate) inner: ck::UffConformerResult,
}
#[pymethods]
impl UffConformerResult {
    #[new]
    fn new(conformer_id: usize, status: i32, energy: f64) -> Self {
        Self {
            inner: ck::UffConformerResult {
                conformer_id,
                status,
                energy,
            },
        }
    }
    fn conformer_id(&self) -> usize {
        self.inner.conformer_id()
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct UffConformerOptimizationResult {
    pub(crate) inner: ck::UffConformerOptimizationResult,
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
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conformer_results(&self) -> Vec<UffConformerResult> {
        self.inner
            .conformer_results()
            .iter()
            .cloned()
            .map(|inner| UffConformerResult { inner })
            .collect()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct MmffOptimizeMoleculeResult {
    pub(crate) inner: ck::MmffOptimizeMoleculeResult,
}
#[pymethods]
impl MmffOptimizeMoleculeResult {
    #[new]
    fn new(molecule: &Molecule, needs_more: i32) -> Self {
        Self {
            inner: ck::MmffOptimizeMoleculeResult {
                molecule: molecule.inner.clone(),
                needs_more,
            },
        }
    }
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct MmffOptimizeMoleculeConfResult {
    pub(crate) inner: ck::MmffOptimizeMoleculeConfResult,
}
#[pymethods]
impl MmffOptimizeMoleculeConfResult {
    #[new]
    fn new(needs_more: i32, energy: f64) -> Self {
        Self {
            inner: ck::MmffOptimizeMoleculeConfResult { needs_more, energy },
        }
    }
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    fn energy(&self) -> f64 {
        self.inner.energy()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct MmffOptimizeMoleculeConfsResult {
    pub(crate) inner: ck::MmffOptimizeMoleculeConfsResult,
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
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conformer_results(&self) -> Vec<MmffOptimizeMoleculeConfResult> {
        self.inner
            .conformer_results()
            .iter()
            .cloned()
            .map(|inner| MmffOptimizeMoleculeConfResult { inner })
            .collect()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct CrippenTotals {
    pub(crate) inner: ck::CrippenTotals,
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

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct LabuteAsaContributions {
    pub(crate) inner: ck::LabuteAsaContributions,
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
