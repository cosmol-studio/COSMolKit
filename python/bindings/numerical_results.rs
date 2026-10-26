//! Detached source result transport. No energy or descriptor calculation lives here.
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Force-field optimization status with an energy only when parameters and an evaluated outcome are available.
    pub struct OptimizationOutcome { pub(crate) inner: ck::OptimizationOutcome }
}
#[pymethods]
impl OptimizationOutcome {
    #[staticmethod]
    /// Construct a converged outcome with its evaluated energy.
    fn converged_with_energy(energy: f64) -> Self {
        Self {
            inner: ck::OptimizationOutcome::converged_with_energy(energy),
        }
    }
    #[staticmethod]
    /// Construct a nonconverged outcome retaining its last evaluated energy.
    fn not_converged_with_energy(energy: f64) -> Self {
        Self {
            inner: ck::OptimizationOutcome::not_converged_with_energy(energy),
        }
    }
    #[staticmethod]
    /// Construct an outcome indicating that required force-field parameters were unavailable.
    fn missing_parameters() -> Self {
        Self {
            inner: ck::OptimizationOutcome::missing_parameters(),
        }
    }
    #[getter]
    /// Typed optimization status: converged, not converged or missing parameters.
    fn status(&self) -> crate::leaf_values::OptimizationStatus {
        self.inner.status().into()
    }
    #[getter]
    /// Evaluated force-field energy; optimization outcomes without parameters return None.
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
    #[getter]
    /// Whether the optimizer reached its convergence criteria.
    fn converged(&self) -> bool {
        self.inner.converged()
    }
    #[getter]
    /// Whether the native outcome reports nonconvergence requiring further optimization.
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    /// Native optimization status code corresponding to the typed outcome.
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
}

macro_rules! gradient_result {
    ($name:ident) => {
        crate::configuration::type_projection! {
            #[pyclass(module = "cosmolkit", frozen, from_py_object)]
            #[derive(Clone)]
            /// Force-field energy and Cartesian gradient at the evaluated coordinates.
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
            /// Evaluated force-field energy; optimization outcomes without parameters return None.
            fn energy(&self) -> f64 {
                self.inner.energy()
            }
            #[getter]
            /// Flat Cartesian energy-gradient components in atom order (x, y, z per atom).
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
    /// UFF optimization result containing the output molecule and optimization outcome.
    pub struct UffOptimizationResult {
        pub(crate) inner: ck::UffOptimizationResult,
    }
}
#[pymethods]
impl UffOptimizationResult {
    #[getter]
    /// Typed optimization outcome retaining both convergence status and optional energy.
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
    /// Output molecule as a copy-on-write value; editing this value does not alter the result.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    /// Native optimization status code corresponding to the typed outcome.
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    /// Whether the native outcome reports nonconvergence requiring further optimization.
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    /// Evaluated force-field energy; optimization outcomes without parameters return None.
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// UFF optimization outcome associated with one conformer identifier.
    pub struct UffConformerResult {
        pub(crate) inner: ck::UffConformerResult,
    }
}
#[pymethods]
impl UffConformerResult {
    #[getter]
    /// Typed optimization outcome retaining both convergence status and optional energy.
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
    /// Identifier of the conformer associated with this outcome.
    fn conformer_id(&self) -> usize {
        self.inner.conformer_id()
    }
    #[getter]
    /// Native optimization status code corresponding to the typed outcome.
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    /// Whether the native outcome reports nonconvergence requiring further optimization.
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    /// Evaluated force-field energy; optimization outcomes without parameters return None.
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// UFF multi-conformer optimization result with the output molecule and ordered per-conformer outcomes.
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
    /// Output molecule as a copy-on-write value; editing this value does not alter the result.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    /// Independent list of per-conformer outcomes in native result order.
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
    /// MMFF optimization result containing the output molecule and optimization outcome.
    pub struct MmffOptimizeMoleculeResult {
        pub(crate) inner: ck::MmffOptimizeMoleculeResult,
    }
}
#[pymethods]
impl MmffOptimizeMoleculeResult {
    #[getter]
    /// Evaluated force-field energy; optimization outcomes without parameters return None.
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
    #[getter]
    /// Typed optimization outcome retaining both convergence status and optional energy.
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
    /// Output molecule as a copy-on-write value; editing this value does not alter the result.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    /// Native optimization status code corresponding to the typed outcome.
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    /// Whether the native outcome reports nonconvergence requiring further optimization.
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// MMFF outcome for one conformer optimization.
    pub struct MmffOptimizeMoleculeConfResult {
        pub(crate) inner: ck::MmffOptimizeMoleculeConfResult,
    }
}
#[pymethods]
impl MmffOptimizeMoleculeConfResult {
    #[getter]
    /// Typed optimization outcome retaining both convergence status and optional energy.
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
    /// Native optimization status code corresponding to the typed outcome.
    fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[getter]
    /// Whether the native outcome reports nonconvergence requiring further optimization.
    fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    #[getter]
    /// Evaluated force-field energy; optimization outcomes without parameters return None.
    fn energy(&self) -> Option<f64> {
        self.inner.energy()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// MMFF multi-conformer optimization result with the output molecule and ordered per-conformer outcomes.
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
    /// Output molecule as a copy-on-write value; editing this value does not alter the result.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    /// Independent list of per-conformer outcomes in native result order.
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
    /// Molecular Crippen logP and molar refractivity totals.
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
    /// Crippen estimate of the octanol/water partition coefficient logarithm.
    fn logp(&self) -> f64 {
        self.inner.logp
    }
    #[getter]
    /// Crippen molecular refractivity estimate.
    fn molar_refractivity(&self) -> f64 {
        self.inner.molar_refractivity
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Labute approximate surface area with per-atom and aggregate hydrogen contributions.
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
    /// Total Labute approximate molecular surface area.
    fn asa(&self) -> f64 {
        self.inner.asa
    }
    #[getter]
    /// Independent list of Labute surface-area contributions in atom order.
    fn atom_contributions(&self) -> Vec<f64> {
        self.inner.atom_contributions.clone()
    }
    #[getter]
    /// Aggregate hydrogen contribution to the Labute surface-area result.
    fn hydrogen_contribution(&self) -> f64 {
        self.inner.hydrogen_contribution
    }
}
