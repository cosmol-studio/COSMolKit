//! Readonly projections of the facade's generated declaration tables.
use crate::feature_metadata::FeatureSpec;
use crate::public_selectors::ParityPolicy;
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
pub struct FunctionStatus {
    inner: ck::FunctionStatus,
}

#[pymethods]
impl FunctionStatus {
    #[classattr]
    #[pyo3(name = "Native")]
    const NATIVE: Self = Self {
        inner: ck::FunctionStatus::Native,
    };
    #[classattr]
    #[pyo3(name = "Experimental")]
    const EXPERIMENTAL: Self = Self {
        inner: ck::FunctionStatus::Experimental,
    };
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
            ck::FunctionStatus::Parity { .. }
            | ck::FunctionStatus::Native
            | ck::FunctionStatus::Experimental => None,
        }
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

#[pyclass(module = "cosmolkit", frozen)]
pub struct MoleculeOpSpec {
    inner: &'static ck::MoleculeOpSpec,
}

impl MoleculeOpSpec {
    pub(crate) fn from_native(inner: &'static ck::MoleculeOpSpec) -> Self {
        Self { inner }
    }
}

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
    fn result_type(&self) -> &'static str {
        self.inner.result_type
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

#[pyclass(module = "cosmolkit", frozen)]
pub struct SupportMatrixEntry {
    inner: &'static ck::SupportMatrixEntry,
}

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

#[pyclass(module = "cosmolkit", frozen)]
pub struct OperationInvariantEntry {
    inner: &'static ck::OperationInvariantEntry,
}

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

#[pyclass(module = "cosmolkit", frozen)]
pub struct ParityMatrixEntry {
    inner: &'static ck::ParityMatrixEntry,
}

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

#[pyfunction]
pub fn operation_specs() -> Vec<MoleculeOpSpec> {
    ck::operation_specs()
        .iter()
        .map(|&inner| MoleculeOpSpec { inner })
        .collect()
}
#[pyfunction]
pub fn operation_spec(method: &str) -> Option<MoleculeOpSpec> {
    ck::operation_spec(method).map(|inner| MoleculeOpSpec { inner })
}
#[pyfunction]
pub fn support_matrix() -> Vec<SupportMatrixEntry> {
    ck::support_matrix()
        .iter()
        .map(|inner| SupportMatrixEntry { inner })
        .collect()
}
#[pyfunction]
pub fn operation_invariant_matrix() -> Vec<OperationInvariantEntry> {
    ck::operation_invariant_matrix()
        .iter()
        .map(|inner| OperationInvariantEntry { inner })
        .collect()
}
#[pyfunction]
pub fn parity_matrix() -> Vec<ParityMatrixEntry> {
    ck::parity_matrix()
        .iter()
        .map(|inner| ParityMatrixEntry { inner })
        .collect()
}
#[pyfunction]
pub fn operation_invariant(method: &str) -> Option<OperationInvariantEntry> {
    ck::operation_invariant(method).map(|inner| OperationInvariantEntry { inner })
}
#[pyfunction]
pub fn operation_parity(method: &str) -> Option<ParityMatrixEntry> {
    ck::operation_parity(method).map(|inner| ParityMatrixEntry { inner })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn operation_metadata_keeps_actual_table_order_and_static_links() {
        assert_eq!(operation_specs().len(), ck::operation_specs().len());
        for (actual, &native) in operation_specs().iter().zip(ck::operation_specs()) {
            assert!(std::ptr::eq(actual.inner, native));
            assert_eq!(actual.method(), native.method);
            assert_eq!(actual.status().inner, native.status);
            assert!(std::ptr::eq(
                operation_spec(native.method).unwrap().inner,
                native
            ));
        }
        for (actual, native) in support_matrix().iter().zip(ck::support_matrix()) {
            assert!(std::ptr::eq(actual.inner, native));
            assert!(std::ptr::eq(actual.feature().inner, native.feature));
            assert_eq!(
                actual.operation().map(|value| value.inner),
                native.operation
            );
        }
        assert_eq!(support_matrix().len(), ck::support_matrix().len());
        assert_eq!(
            operation_invariant_matrix().len(),
            ck::operation_invariant_matrix().len()
        );
        for (actual, native) in operation_invariant_matrix()
            .iter()
            .zip(ck::operation_invariant_matrix())
        {
            assert!(std::ptr::eq(actual.inner, native));
            assert!(std::ptr::eq(
                operation_invariant(native.operation.method).unwrap().inner,
                native
            ));
        }
        assert_eq!(parity_matrix().len(), ck::parity_matrix().len());
        for (actual, native) in parity_matrix().iter().zip(ck::parity_matrix()) {
            assert!(std::ptr::eq(actual.inner, native));
            assert!(std::ptr::eq(
                operation_parity(native.operation.method).unwrap().inner,
                native
            ));
        }
        assert!(operation_spec("not-a-native-operation").is_none());
        assert!(operation_invariant("not-a-native-operation").is_none());
        assert!(operation_parity("not-a-native-operation").is_none());
        for (inner, reference, explanation) in [
            (ck::FunctionStatus::Native, None, None),
            (ck::FunctionStatus::Experimental, None, None),
            (
                ck::FunctionStatus::Parity {
                    reference: "source",
                },
                Some("source"),
                None,
            ),
            (
                ck::FunctionStatus::ParityWithDifferences {
                    reference: "source",
                    explanation: "exact scope",
                },
                Some("source"),
                Some("exact scope"),
            ),
        ] {
            let actual = FunctionStatus { inner };
            assert_eq!(actual.reference(), reference);
            assert_eq!(actual.explanation(), explanation);
        }
    }
}
