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
        #[pyclass(module="cosmolkit", frozen, from_py_object, eq)]
        #[derive(Clone, Copy, PartialEq)]
        /// Bit mask selecting sanitization stages; combine constants with the bitwise OR operator.
        pub struct SanitizeOperations { pub(crate) inner: ck::SanitizeOperations }
        #[pymethods]
        impl SanitizeOperations {
            $(#[classattr] #[allow(non_snake_case)] fn $flag() -> Self { Self { inner: ck::SanitizeOperations::$flag } })+
            #[staticmethod]
            /// Construct a sanitization mask, rejecting bits not recognized by the native API.
            fn from_bits(bits: u32) -> PyResult<Self> {
                ck::SanitizeOperations::from_bits(bits)
                    .map(|inner| Self { inner })
                    .map_err(crate::sanitize_errors::SanitizeError::from_rust)
            }
            #[getter]
            /// Integer representation of this option mask.
            fn bits(&self) -> u32 { self.inner.bits() }
            /// Whether every bit in the supplied mask is present.
            fn contains(&self, operation: &Self) -> bool { self.inner.contains(operation.inner) }
            #[getter]
            /// Whether the mask contains no enabled options.
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

crate::configuration::type_projection! {
    impl Configuration for SanitizeOperations {
        type Native = ck::SanitizeOperations;
        const INPUT: pyo3::inspect::PyStaticExpr =
            <SanitizeOperationsInput as FromPyObject<'_, '_>>::INPUT_TYPE;
        fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
            value.extract::<SanitizeOperationsInput>()?.into_rust()
        }
        fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
            Ok(Self { inner })
        }
        fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
            Ok(self.inner)
        }
    }
}

ck::project_configuration!(ck::SanitizeParams);
