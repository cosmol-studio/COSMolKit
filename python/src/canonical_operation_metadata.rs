//! Readonly projections of the facade's generated operation metadata.
//!
//! Enum-valued declaration fields use their Rust variant names. Block and
//! derived-state sets retain their exact bit masks; no permissions are granted.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};
use std::collections::BTreeMap;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FeatureSpec {
    inner: &'static ck::FeatureSpec,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FeatureSpec {
    #[getter]
    fn name(&self) -> &'static str {
        self.inner.name
    }
    #[getter]
    fn category(&self) -> &'static str {
        self.inner.category
    }
    #[getter]
    fn docs(&self) -> &'static str {
        self.inner.docs
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FunctionStatus {
    inner: ck::FunctionStatus,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FunctionStatus {
    #[getter]
    fn kind(&self) -> &'static str {
        match self.inner {
            ck::FunctionStatus::Parity { .. } => "Parity",
            ck::FunctionStatus::ParityWithDifferences { .. } => "ParityWithDifferences",
            ck::FunctionStatus::Native => "Native",
            ck::FunctionStatus::Experimental => "Experimental",
        }
    }
    #[getter]
    fn reference(&self) -> Option<&'static str> {
        match self.inner {
            ck::FunctionStatus::Parity { reference }
            | ck::FunctionStatus::ParityWithDifferences { reference, .. } => Some(reference),
            ck::FunctionStatus::Native | ck::FunctionStatus::Experimental => None,
        }
    }
    #[getter]
    fn explanation(&self) -> Option<&'static str> {
        match self.inner {
            ck::FunctionStatus::ParityWithDifferences { explanation, .. } => Some(explanation),
            _ => None,
        }
    }
}

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub(crate) enum ParityPolicy {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}

impl From<ck::ParityPolicy> for ParityPolicy {
    fn from(value: ck::ParityPolicy) -> Self {
        match value {
            ck::ParityPolicy::NotApplicable => Self::NotApplicable,
            ck::ParityPolicy::RequiredWhenSupported => Self::RequiredWhenSupported,
            ck::ParityPolicy::RequiredNow => Self::RequiredNow,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MoleculeOpSpec {
    inner: &'static ck::MoleculeOpSpec,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeOpSpec {
    #[getter]
    fn method(&self) -> &'static str {
        self.inner.method
    }
    #[getter]
    fn impl_fn(&self) -> &'static str {
        self.inner.impl_fn
    }
    #[getter]
    fn output(&self) -> String {
        format!("{:?}", self.inner.output)
    }
    #[getter]
    fn result_type(&self) -> &'static str {
        self.inner.result_type
    }
    #[getter]
    fn domain(&self) -> String {
        format!("{:?}", self.inner.domain)
    }
    #[getter]
    fn kind(&self) -> String {
        format!("{:?}", self.inner.kind)
    }
    #[getter]
    fn topology_edit(&self) -> String {
        format!("{:?}", self.inner.topology_edit)
    }
    /// Exact read/write block bit masks from the generated declaration.
    #[getter]
    fn access(&self) -> BTreeMap<String, u8> {
        BTreeMap::from([
            ("read".into(), self.inner.access.read().bits()),
            ("write".into(), self.inner.access.write().bits()),
        ])
    }
    #[getter]
    fn may_mutate(&self) -> u8 {
        self.inner.may_mutate.bits()
    }
    #[getter]
    fn auto_remap(&self) -> u8 {
        self.inner.auto_remap.bits()
    }
    /// Exact derived-state masks, including operation-defined transitions.
    #[getter]
    fn derived_effects(&self) -> BTreeMap<String, u16> {
        let effects = self.inner.derived_effects;
        BTreeMap::from([
            ("recompute".into(), effects.recompute.bits()),
            ("preserve".into(), effects.preserve.bits()),
            ("invalidate".into(), effects.invalidate.bits()),
            ("operation_defined".into(), effects.operation_defined.bits()),
        ])
    }
    #[getter]
    fn cip_state(&self) -> String {
        format!("{:?}", self.inner.cip_state)
    }
    #[getter]
    fn semantic_preconditions(&self) -> u8 {
        self.inner.semantic_preconditions.bits()
    }
    #[getter]
    fn requires_mapping(&self) -> String {
        format!("{:?}", self.inner.requires_mapping)
    }
    #[getter]
    fn status(&self) -> FunctionStatus {
        FunctionStatus {
            inner: self.inner.status,
        }
    }
    #[getter]
    fn parity(&self) -> ParityPolicy {
        self.inner.parity.into()
    }
    #[getter]
    fn io_roundtrip(&self) -> bool {
        self.inner.io_roundtrip
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SupportMatrixEntry {
    inner: &'static ck::SupportMatrixEntry,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SupportMatrixEntry {
    #[getter]
    fn feature(&self) -> FeatureSpec {
        FeatureSpec {
            inner: self.inner.feature,
        }
    }
    #[getter]
    fn operation(&self) -> Option<MoleculeOpSpec> {
        self.inner.operation.map(|inner| MoleculeOpSpec { inner })
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct OperationInvariantEntry {
    inner: &'static ck::OperationInvariantEntry,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl OperationInvariantEntry {
    #[getter]
    fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec {
            inner: self.inner.operation,
        }
    }
    #[getter]
    fn profile(&self) -> &'static str {
        self.inner.profile
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ParityMatrixEntry {
    inner: &'static ck::ParityMatrixEntry,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ParityMatrixEntry {
    #[getter]
    fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec {
            inner: self.inner.operation,
        }
    }
    #[getter]
    fn feature(&self) -> FeatureSpec {
        FeatureSpec {
            inner: self.inner.feature,
        }
    }
    #[getter]
    fn profile(&self) -> &'static str {
        self.inner.profile
    }
    #[getter]
    fn rdkit_version(&self) -> Option<&'static str> {
        self.inner.rdkit_version
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct FeatureSpecIter {
    inner: ck::FeatureSpecIter,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FeatureSpecIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self) -> Option<FeatureSpec> {
        self.inner.next().map(|inner| FeatureSpec { inner })
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter {
        inner: ck::feature_specs(),
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn feature_spec(name: &str) -> Option<FeatureSpec> {
    ck::feature_spec(name).map(|inner| FeatureSpec { inner })
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_specs() -> Vec<MoleculeOpSpec> {
    ck::operation_specs()
        .iter()
        .map(|&inner| MoleculeOpSpec { inner })
        .collect()
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_spec(method: &str) -> Option<MoleculeOpSpec> {
    ck::operation_spec(method).map(|inner| MoleculeOpSpec { inner })
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn support_matrix() -> Vec<SupportMatrixEntry> {
    ck::support_matrix()
        .iter()
        .map(|inner| SupportMatrixEntry { inner })
        .collect()
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_invariant_matrix() -> Vec<OperationInvariantEntry> {
    ck::operation_invariant_matrix()
        .iter()
        .map(|inner| OperationInvariantEntry { inner })
        .collect()
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parity_matrix() -> Vec<ParityMatrixEntry> {
    ck::parity_matrix()
        .iter()
        .map(|inner| ParityMatrixEntry { inner })
        .collect()
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_invariant(method: &str) -> Option<OperationInvariantEntry> {
    ck::operation_invariant(method).map(|inner| OperationInvariantEntry { inner })
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_parity(method: &str) -> Option<ParityMatrixEntry> {
    ck::operation_parity(method).map(|inner| ParityMatrixEntry { inner })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<FeatureSpec>()?;
    module.add_class::<FeatureSpecIter>()?;
    module.add_class::<MoleculeOpSpec>()?;
    module.add_class::<SupportMatrixEntry>()?;
    module.add_class::<OperationInvariantEntry>()?;
    module.add_class::<ParityMatrixEntry>()?;
    module.add_class::<FunctionStatus>()?;
    module.add_class::<ParityPolicy>()?;
    module.add_function(wrap_pyfunction!(feature_specs, module)?)?;
    module.add_function(wrap_pyfunction!(feature_spec, module)?)?;
    module.add_function(wrap_pyfunction!(operation_specs, module)?)?;
    module.add_function(wrap_pyfunction!(operation_spec, module)?)?;
    module.add_function(wrap_pyfunction!(support_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(operation_invariant_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(parity_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(operation_invariant, module)?)?;
    module.add_function(wrap_pyfunction!(operation_parity, module)?)?;
    Ok(())
}
