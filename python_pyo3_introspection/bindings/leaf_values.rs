use crate::chemical_values::enum_value;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

enum_value!(TautomerEnumerationStatus, { Completed, MaxTautomersReached, MaxTransformsReached, Canceled });
enum_value!(UffParameterErrorKind, { Preparation, ParameterTable, Typing });
enum_value!(MmffVariant, { Mmff94, Mmff94s });

#[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
pub struct UffOptimizationErrorKind {
    pub(crate) inner: ck::UffOptimizationErrorKind,
}
#[pymethods]
impl UffOptimizationErrorKind {
    #[classattr]
    #[pyo3(name = "Rings")]
    const RINGS: Self = Self {
        inner: ck::UffOptimizationErrorKind::Rings,
    };
    #[classattr]
    #[pyo3(name = "Optimization")]
    const OPTIMIZATION: Self = Self {
        inner: ck::UffOptimizationErrorKind::Optimization,
    };
    #[classattr]
    #[pyo3(name = "ConformerOptimization")]
    const CONFORMER_OPTIMIZATION: Self = Self {
        inner: ck::UffOptimizationErrorKind::ConformerOptimization,
    };
    #[classattr]
    #[pyo3(name = "Evaluation")]
    const EVALUATION: Self = Self {
        inner: ck::UffOptimizationErrorKind::Evaluation,
    };
    #[staticmethod]
    #[pyo3(name="MissingConformer", signature=(requested))]
    fn missing_conformer(requested: Option<usize>) -> Self {
        Self {
            inner: ck::UffOptimizationErrorKind::MissingConformer { requested },
        }
    }
    #[getter]
    fn requested(&self) -> Option<usize> {
        match self.inner {
            ck::UffOptimizationErrorKind::MissingConformer { requested } => requested,
            ck::UffOptimizationErrorKind::Rings
            | ck::UffOptimizationErrorKind::Optimization
            | ck::UffOptimizationErrorKind::ConformerOptimization
            | ck::UffOptimizationErrorKind::Evaluation => None,
        }
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct TopologicalFingerprintOutput {
    pub(crate) inner: ck::TopologicalFingerprintOutput,
}
#[pymethods]
impl TopologicalFingerprintOutput {
    #[new]
    #[pyo3(signature=(atom_bits, bit_info))]
    fn new(
        atom_bits: Option<Vec<Vec<u32>>>,
        bit_info: Option<BTreeMap<u32, Vec<Vec<i32>>>>,
    ) -> Self {
        Self {
            inner: ck::TopologicalFingerprintOutput {
                atom_bits,
                bit_info,
            },
        }
    }
    #[getter]
    fn atom_bits(&self) -> Option<Vec<Vec<u32>>> {
        self.inner.atom_bits.clone()
    }
    #[getter]
    fn bit_info(&self) -> Option<BTreeMap<u32, Vec<Vec<i32>>>> {
        self.inner.bit_info.clone()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct MmffAtomProperties {
    pub(crate) inner: ck::MmffAtomProperties,
}
#[pymethods]
impl MmffAtomProperties {
    #[new]
    fn new(atom_type: u8, formal_charge: f64, partial_charge: f64) -> Self {
        Self {
            inner: ck::MmffAtomProperties {
                atom_type,
                formal_charge,
                partial_charge,
            },
        }
    }
    fn atom_type(&self) -> u8 {
        self.inner.atom_type()
    }
    fn formal_charge(&self) -> f64 {
        self.inner.formal_charge()
    }
    fn partial_charge(&self) -> f64 {
        self.inner.partial_charge()
    }
}
