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

/// Read-only description of a compiled public feature selection.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FeatureSpec {
    inner: &'static ck::FeatureSpec,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FeatureSpec {
    /// Stored name of this value.
    #[getter]
    fn name(&self) -> &'static str {
        self.inner.name
    }
    /// Category of this compiled public feature.
    #[getter]
    fn category(&self) -> &'static str {
        self.inner.category
    }
    /// Description of the feature from its canonical declaration.
    #[getter]
    fn docs(&self) -> &'static str {
        self.inner.docs
    }
}

/// Declared behavior commitment and reference library, not a test-pass result.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FunctionStatus {
    inner: ck::FunctionStatus,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FunctionStatus {
    /// Return the declared status: Parity, ParityWithDifferences, Native or Experimental; this is not a test-pass claim.
    #[getter]
    fn kind(&self) -> &'static str {
        match self.inner {
            ck::FunctionStatus::Parity { .. } => "Parity",
            ck::FunctionStatus::ParityWithDifferences { .. } => "ParityWithDifferences",
            ck::FunctionStatus::Native => "Native",
            ck::FunctionStatus::Experimental => "Experimental",
        }
    }
    /// Reference library named by the function behavior commitment, or None for native/experimental behavior.
    #[getter]
    fn reference(&self) -> Option<&'static str> {
        match self.inner {
            ck::FunctionStatus::Parity { reference }
            | ck::FunctionStatus::ParityWithDifferences { reference, .. } => Some(reference),
            ck::FunctionStatus::Native | ck::FunctionStatus::Experimental => None,
        }
    }
    /// Approved behavior difference explanation, or None when no difference is declared.
    #[getter]
    fn explanation(&self) -> Option<&'static str> {
        match self.inner {
            ck::FunctionStatus::ParityWithDifferences { explanation, .. } => Some(explanation),
            _ => None,
        }
    }
}

/// Declared reference-validation policy for an operation; does not indicate that any tests have run.
///
/// Declared values: ``NotApplicable``, ``RequiredWhenSupported``, ``RequiredNow``.
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

/// Read-only operation declaration metadata; inspecting it does not grant mutation or runtime access.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MoleculeOpSpec {
    inner: &'static ck::MoleculeOpSpec,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeOpSpec {
    /// Experimental method, or canonical operation name in operation metadata.
    #[getter]
    fn method(&self) -> &'static str {
        self.inner.method
    }
    /// Implementation function name recorded by the operation declaration; not a callable Python entry point.
    #[getter]
    fn impl_fn(&self) -> &'static str {
        self.inner.impl_fn
    }
    /// Output multiplicity/category recorded by the operation declaration.
    #[getter]
    fn output(&self) -> String {
        format!("{:?}", self.inner.output)
    }
    /// Declared operation result type name.
    #[getter]
    fn result_type(&self) -> &'static str {
        self.inner.result_type
    }
    /// Chemical/API domain associated with the operation or error.
    #[getter]
    fn domain(&self) -> String {
        format!("{:?}", self.inner.domain)
    }
    /// Classification/discriminant of this value as defined by its owning type.
    #[getter]
    fn kind(&self) -> String {
        format!("{:?}", self.inner.kind)
    }
    /// Declared topology-edit category of this operation.
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
    /// Declared graph/coordinate/property blocks that the operation may change.
    #[getter]
    fn may_mutate(&self) -> u8 {
        self.inner.may_mutate.bits()
    }
    /// Whether the runtime automatically remaps referenced data for this operation.
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
    /// Declared effect on cached CIP assignment state.
    #[getter]
    fn cip_state(&self) -> String {
        format!("{:?}", self.inner.cip_state)
    }
    /// Semantic preconditions declared for the operation.
    #[getter]
    fn semantic_preconditions(&self) -> u8 {
        self.inner.semantic_preconditions.bits()
    }
    /// Whether the operation must produce an atom/bond topology mapping.
    #[getter]
    fn requires_mapping(&self) -> String {
        format!("{:?}", self.inner.requires_mapping)
    }
    /// Declared status of the operation/result; inspect the typed value rather than inferring success from a message.
    #[getter]
    fn status(&self) -> FunctionStatus {
        FunctionStatus {
            inner: self.inner.status,
        }
    }
    /// Declared reference-comparison policy; not a test execution result.
    #[getter]
    fn parity(&self) -> ParityPolicy {
        self.inner.parity.into()
    }
    /// Declared input/output preservation commitment.
    #[getter]
    fn io_roundtrip(&self) -> bool {
        self.inner.io_roundtrip
    }
}

/// Read-only feature/support declaration associated with a public operation.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SupportMatrixEntry {
    inner: &'static ck::SupportMatrixEntry,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SupportMatrixEntry {
    /// Owning compiled capability recorded for this declaration.
    #[getter]
    fn feature(&self) -> FeatureSpec {
        FeatureSpec {
            inner: self.inner.feature,
        }
    }
    /// Return the operation declaration associated with this support entry, when available.
    #[getter]
    fn operation(&self) -> Option<MoleculeOpSpec> {
        self.inner.operation.map(|inner| MoleculeOpSpec { inner })
    }
}

/// Read-only invariant-validation profile associated with a public operation.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct OperationInvariantEntry {
    inner: &'static ck::OperationInvariantEntry,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl OperationInvariantEntry {
    /// Name of the operation that produced this error.
    #[getter]
    fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec {
            inner: self.inner.operation,
        }
    }
    /// Validation/reference profile recorded by this matrix entry.
    #[getter]
    fn profile(&self) -> &'static str {
        self.inner.profile
    }
}

/// Read-only reference-comparison declaration for a public operation; not a test result.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ParityMatrixEntry {
    inner: &'static ck::ParityMatrixEntry,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ParityMatrixEntry {
    /// Name of the operation that produced this error.
    #[getter]
    fn operation(&self) -> MoleculeOpSpec {
        MoleculeOpSpec {
            inner: self.inner.operation,
        }
    }
    /// Owning compiled capability recorded for this declaration.
    #[getter]
    fn feature(&self) -> FeatureSpec {
        FeatureSpec {
            inner: self.inner.feature,
        }
    }
    /// Validation/reference profile recorded by this matrix entry.
    #[getter]
    fn profile(&self) -> &'static str {
        self.inner.profile
    }
    /// Reference RDKit version recorded by this parity matrix entry.
    #[getter]
    fn rdkit_version(&self) -> Option<&'static str> {
        self.inner.rdkit_version
    }
}

/// Allocation-free, declaration-ordered projection of unique feature specs.
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

/// Returns each feature referenced by the generated support matrix once.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn feature_specs() -> FeatureSpecIter {
    FeatureSpecIter {
        inner: ck::feature_specs(),
    }
}
/// Looks up a generated feature by its exact canonical name.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn feature_spec(name: &str) -> Option<FeatureSpec> {
    ck::feature_spec(name).map(|inner| FeatureSpec { inner })
}
/// Returns the generated operation declarations in declaration order.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_specs() -> Vec<MoleculeOpSpec> {
    ck::operation_specs()
        .iter()
        .map(|&inner| MoleculeOpSpec { inner })
        .collect()
}
/// Looks up a generated operation by its exact canonical method name.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_spec(method: &str) -> Option<MoleculeOpSpec> {
    ck::operation_spec(method).map(|inner| MoleculeOpSpec { inner })
}
/// Returns the generated support matrix without copying its rows.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn support_matrix() -> Vec<SupportMatrixEntry> {
    ck::support_matrix()
        .iter()
        .map(|inner| SupportMatrixEntry { inner })
        .collect()
}
/// Returns the generated operation-invariant matrix without copying its rows.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_invariant_matrix() -> Vec<OperationInvariantEntry> {
    ck::operation_invariant_matrix()
        .iter()
        .map(|inner| OperationInvariantEntry { inner })
        .collect()
}
/// Returns the generated parity matrix without copying its rows.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parity_matrix() -> Vec<ParityMatrixEntry> {
    ck::parity_matrix()
        .iter()
        .map(|inner| ParityMatrixEntry { inner })
        .collect()
}
/// Looks up an invariant row by the exact generated operation method.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn operation_invariant(method: &str) -> Option<OperationInvariantEntry> {
    ck::operation_invariant(method).map(|inner| OperationInvariantEntry { inner })
}
/// Looks up a parity row by the exact generated operation method.
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
