use crate::chemical_values::{EnumValue, enum_value};
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(SdfFormat, { V2000, V3000 });
enum_value!(Mol2Type, { Corina });
enum_value!(BatchErrorMode, { Strict, KeepErrors });

ck::project_configuration!(ck::BatchImageParams);
#[derive(FromPyObject, IntoPyObject)]
pub enum SdfFormatInput {
    Enum(SdfFormat),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum Mol2TypeInput {
    Enum(Mol2Type),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum BatchErrorModeInput {
    Enum(BatchErrorMode),
    Name(String),
}
macro_rules! input {
    ($input:ident, $value:ident) => {
        impl $input {
            pub(crate) fn into_rust(self) -> PyResult<ck::$value> {
                match self {
                    Self::Enum(value) => Ok(value.into()),
                    Self::Name(value) => $value::from_name(&value).map(Into::into),
                }
            }
        }
    };
}
input!(SdfFormatInput, SdfFormat);
input!(Mol2TypeInput, Mol2Type);
input!(BatchErrorModeInput, BatchErrorMode);

#[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
pub struct MolCoordinateSelection {
    pub(crate) inner: ck::MolCoordinateSelection,
}
impl From<ck::MolCoordinateSelection> for MolCoordinateSelection {
    fn from(inner: ck::MolCoordinateSelection) -> Self {
        Self { inner }
    }
}
#[pymethods]
impl MolCoordinateSelection {
    #[classattr]
    #[pyo3(name = "Auto")]
    fn auto() -> Self {
        ck::MolCoordinateSelection::Auto.into()
    }
    #[staticmethod]
    #[pyo3(name = "TwoD")]
    fn two_d(id: usize) -> Self {
        ck::MolCoordinateSelection::TwoD { id }.into()
    }
    #[staticmethod]
    #[pyo3(name = "ThreeD")]
    fn three_d(id: usize) -> Self {
        ck::MolCoordinateSelection::ThreeD { id }.into()
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        match self.inner {
            ck::MolCoordinateSelection::Auto => None,
            ck::MolCoordinateSelection::TwoD { id } | ck::MolCoordinateSelection::ThreeD { id } => {
                Some(id)
            }
        }
    }
}
#[derive(FromPyObject, IntoPyObject)]
pub enum MolCoordinateSelectionInput {
    Selection(MolCoordinateSelection),
    Name(String),
}
impl MolCoordinateSelectionInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::MolCoordinateSelection> {
        match self {
            Self::Selection(value) => Ok(value.inner),
            Self::Name(value) if value == "auto" => Ok(ck::MolCoordinateSelection::Auto),
            Self::Name(value) => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown coordinate selection {value}; use TwoD(id) or ThreeD(id) for explicit IDs"
            ))),
        }
    }
}
impl crate::parameters::Configuration for MolCoordinateSelection {
    type Native = ck::MolCoordinateSelection;
    const INPUT: pyo3::inspect::PyStaticExpr =
        <MolCoordinateSelectionInput as FromPyObject<'_, '_>>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        value.extract::<MolCoordinateSelectionInput>()?.into_rust()
    }
    fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
        Ok(Self { inner })
    }
    fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
        Ok(self.inner)
    }
}
ck::project_configuration!(ck::MolBlockWriteParams);
ck::project_configuration!(ck::PdbWriteParams);
ck::project_configuration!(ck::Mol2ReadParams);

ck::project_configuration!(ck::BatchParams);
ck::project_configuration!(ck::BatchExportParams);
