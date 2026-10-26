use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! cx_fields {
    ($($flag:ident),+ $(,)?) => {
        #[pyclass(module="cosmolkit", frozen, from_py_object, eq)]
        #[derive(Clone, Copy, PartialEq)]
        /// Bit mask selecting CXSMILES extension fields; combine constants with the bitwise OR operator.
        pub struct CxSmilesFields { pub(crate) inner: ck::CxSmilesFields }
        #[pymethods]
        impl CxSmilesFields {
            $(#[classattr] #[allow(non_snake_case)] fn $flag() -> Self { Self { inner: ck::CxSmilesFields::$flag } })+
            #[getter]
            /// Integer representation of this option mask.
            fn bits(&self) -> u32 { self.inner.bits() }
            /// Whether every bit in the supplied mask is present.
            fn contains(&self, other: &Self) -> bool { self.inner.contains(other.inner) }
            fn __or__(&self, other: &Self) -> Self { Self { inner: self.inner | other.inner } }
            fn __index__(&self) -> u32 { self.inner.bits() }
        }
        impl CxSmilesFieldsInput {
            pub(crate) fn into_rust(self) -> PyResult<ck::CxSmilesFields> {
                match self {
                    Self::Fields(value) => Ok(value.inner),
                    Self::Name(value) => {
                        $(if value == stringify!($flag).to_ascii_lowercase() { return Ok(ck::CxSmilesFields::$flag); })+
                        Err(pyo3::exceptions::PyValueError::new_err(format!("unknown CxSmilesFields: {value}")))
                    }
                }
            }
        }
    };
}
#[derive(FromPyObject, IntoPyObject)]
pub enum CxSmilesFieldsInput {
    Fields(CxSmilesFields),
    Name(String),
}
cx_fields! { NONE, ATOM_LABELS, MOLFILE_VALUES, COORDS, RADICALS, ATOM_PROPS,
LINKNODES, ENHANCED_STEREO, SGROUPS, POLYMER, BOND_CFG, BOND_ATROPISOMER,
COORDINATE_BONDS, HYDROGEN_BONDS, ZERO_BONDS, ALL, ALL_BUT_COORDS }

#[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
/// Coordinate selection for CXSMILES: Auto, TwoD(id), or ThreeD(id).
pub struct CxCoordinateSelection {
    pub(crate) inner: ck::CxCoordinateSelection,
}
#[pymethods]
impl CxCoordinateSelection {
    #[classattr]
    #[pyo3(name = "Auto")]
    fn auto() -> Self {
        Self {
            inner: ck::CxCoordinateSelection::Auto,
        }
    }
    #[staticmethod]
    #[pyo3(name = "TwoD")]
    /// Select a stored two-dimensional layout by its dimension-scoped identifier.
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::TwoD { id },
        }
    }
    #[staticmethod]
    #[pyo3(name = "ThreeD")]
    /// Select a stored three-dimensional conformer by its dimension-scoped identifier.
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::ThreeD { id },
        }
    }
    #[getter]
    /// Dimension-scoped coordinate identifier, or None for Auto selection.
    fn id(&self) -> Option<usize> {
        match self.inner {
            ck::CxCoordinateSelection::Auto => None,
            ck::CxCoordinateSelection::TwoD { id } | ck::CxCoordinateSelection::ThreeD { id } => {
                Some(id)
            }
        }
    }
}
#[derive(FromPyObject, IntoPyObject)]
pub enum CxCoordinateSelectionInput {
    Selection(CxCoordinateSelection),
    Name(String),
}
impl CxCoordinateSelectionInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::CxCoordinateSelection> {
        match self {
            Self::Selection(value) => Ok(value.inner),
            Self::Name(value) if value == "auto" => Ok(ck::CxCoordinateSelection::Auto),
            Self::Name(value) => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown coordinate selection {value}; use TwoD(id) or ThreeD(id) for explicit conformer IDs"
            ))),
        }
    }
}

crate::configuration::type_projection! {
    impl Configuration for CxSmilesFields {
        type Native = ck::CxSmilesFields;
        const INPUT: pyo3::inspect::PyStaticExpr =
            <CxSmilesFieldsInput as FromPyObject<'_, '_>>::INPUT_TYPE;
        fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
            value.extract::<CxSmilesFieldsInput>()?.into_rust()
        }
        fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
            Ok(Self { inner })
        }
        fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
            Ok(self.inner)
        }
    }
}
crate::configuration::type_projection! {
    impl Configuration for CxCoordinateSelection {
        type Native = ck::CxCoordinateSelection;
        const INPUT: pyo3::inspect::PyStaticExpr =
            <CxCoordinateSelectionInput as FromPyObject<'_, '_>>::INPUT_TYPE;
        fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
            value.extract::<CxCoordinateSelectionInput>()?.into_rust()
        }
        fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
            Ok(Self { inner })
        }
        fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
            Ok(self.inner)
        }
    }
}

ck::project_configuration!(ck::CxSmilesWriteParams);

ck::project_configuration!(ck::FragmentSmilesWriteParams);
ck::project_configuration!(ck::FragmentCxSmilesWriteParams);

ck::project_configuration!(ck::SmartsParseParams);
