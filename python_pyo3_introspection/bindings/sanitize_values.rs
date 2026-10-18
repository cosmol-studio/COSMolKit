use crate::chemical_values::enum_value;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(SanitizeStage, {
    #[pyo3(name = "None_")] None,
    Cleanup, Properties, SymmRings, Kekulize, FindRadicals, SetAromaticity,
    SetConjugation, SetHybridization, CleanupChirality, AdjustHs,
    CleanupOrganometallics, CleanupAtropisomers
});

macro_rules! sanitize_flags {
    ($($flag:ident),+ $(,)?) => {
        #[pyclass(module="fixed_projection._native", frozen, from_py_object, eq)]
        #[derive(Clone, Copy, PartialEq)]
        pub struct SanitizeOperations { pub(crate) inner: ck::SanitizeOperations }
        #[pymethods]
        impl SanitizeOperations {
            $(#[classattr] #[allow(non_snake_case)] fn $flag() -> Self { Self { inner: ck::SanitizeOperations::$flag } })+
            fn bits(&self) -> u32 { self.inner.bits() }
            fn contains(&self, operation: &Self) -> bool { self.inner.contains(operation.inner) }
            fn is_empty(&self) -> bool { self.inner.is_empty() }
            fn __index__(&self) -> u32 { self.inner.bits() }
            fn __or__(&self, other: &Self) -> Self { Self { inner: self.inner | other.inner } }
            fn __and__(&self, other: &Self) -> Self { Self { inner: self.inner & other.inner } }
        }
        impl SanitizeOperationsInput {
            pub(crate) fn into_rust(self) -> PyResult<ck::SanitizeOperations> {
                match self {
                    Self::Operations(value) => Ok(value.inner),
                    Self::Name(value) => {
                        $(if value == stringify!($flag).to_ascii_lowercase() { return Ok(ck::SanitizeOperations::$flag); })+
                        Err(pyo3::exceptions::PyValueError::new_err(format!("unknown SanitizeOperations: {value}")))
                    }
                }
            }
        }
    };
}
#[derive(FromPyObject, IntoPyObject)]
pub enum SanitizeOperationsInput {
    Operations(SanitizeOperations),
    Name(String),
}
sanitize_flags! { NONE, CLEANUP, PROPERTIES, SYMM_RINGS, KEKULIZE, FIND_RADICALS,
SET_AROMATICITY, SET_CONJUGATION, SET_HYBRIDIZATION, CLEANUP_CHIRALITY,
ADJUST_HS, CLEANUP_ORGANOMETALLICS, CLEANUP_ATROPISOMERS, ALL }

#[pyclass(module = "fixed_projection._native")]
pub struct SanitizeParams {
    pub(crate) inner: ck::SanitizeParams,
}
#[pymethods]
impl SanitizeParams {
    #[new]
    #[pyo3(signature=(*, operations=None))]
    fn new(operations: Option<SanitizeOperationsInput>) -> PyResult<Self> {
        Ok(Self {
            inner: ck::SanitizeParams {
                operations: operations
                    .map(SanitizeOperationsInput::into_rust)
                    .transpose()?
                    .unwrap_or_else(|| ck::SanitizeParams::default().operations),
            },
        })
    }
    #[getter]
    fn operations(&self) -> SanitizeOperations {
        SanitizeOperations {
            inner: self.inner.operations,
        }
    }
    #[setter]
    fn set_operations(&mut self, value: Option<SanitizeOperationsInput>) -> PyResult<()> {
        self.inner.operations = value
            .map(SanitizeOperationsInput::into_rust)
            .transpose()?
            .unwrap_or_else(|| ck::SanitizeParams::default().operations);
        Ok(())
    }
}
