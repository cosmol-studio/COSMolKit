use crate::chemical_values::{EnumValue, enum_value};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::path::PathBuf;

enum_value!(SdfFormat, { V2000, V3000 });
enum_value!(Mol2Type, { Corina });
enum_value!(BatchErrorMode, { Strict, KeepErrors });

#[pyclass(module = "fixed_projection._native")]
pub struct BatchImageParams {
    pub(crate) inner: ck::BatchImageParams,
}
#[pymethods]
impl BatchImageParams {
    #[new]
    #[pyo3(signature=(*, format=ck::BatchImageParams::default().format, width=ck::BatchImageParams::default().width, height=ck::BatchImageParams::default().height, execution=None, filenames=None, report_path=None))]
    fn new(
        format: String,
        width: u32,
        height: u32,
        execution: Option<&BatchParams>,
        filenames: Option<Vec<Option<String>>>,
        report_path: Option<PathBuf>,
    ) -> Self {
        Self {
            inner: ck::BatchImageParams {
                format,
                width,
                height,
                execution: execution
                    .map(|value| value.inner.clone())
                    .unwrap_or_else(|| ck::BatchImageParams::default().execution),
                filenames,
                report_path,
            },
        }
    }
    #[getter]
    fn format(&self) -> &str {
        &self.inner.format
    }
    #[setter]
    fn set_format(&mut self, value: String) {
        self.inner.format = value;
    }
    #[getter]
    fn width(&self) -> u32 {
        self.inner.width
    }
    #[setter]
    fn set_width(&mut self, value: u32) {
        self.inner.width = value;
    }
    #[getter]
    fn height(&self) -> u32 {
        self.inner.height
    }
    #[setter]
    fn set_height(&mut self, value: u32) {
        self.inner.height = value;
    }
    #[getter]
    fn execution(&self) -> BatchParams {
        BatchParams {
            inner: self.inner.execution.clone(),
        }
    }
    #[setter]
    fn set_execution(&mut self, value: Option<&BatchParams>) {
        self.inner.execution = value
            .map(|value| value.inner.clone())
            .unwrap_or_else(|| ck::BatchImageParams::default().execution);
    }
    #[getter]
    fn filenames(&self) -> Option<Vec<Option<String>>> {
        self.inner.filenames.clone()
    }
    #[setter]
    fn set_filenames(&mut self, value: Option<Vec<Option<String>>>) {
        self.inner.filenames = value;
    }
    #[getter]
    fn report_path(&self) -> Option<PathBuf> {
        self.inner.report_path.clone()
    }
    #[setter]
    fn set_report_path(&mut self, value: Option<PathBuf>) {
        self.inner.report_path = value;
    }
}
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

#[pyclass(module = "fixed_projection._native", frozen, from_py_object, eq)]
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
macro_rules! io_parameters {
    ($name:ident, scalars {$($field:ident: $ty:ty),+ $(,)?}, selectors {$($selector:ident: $input:ident => $output:ident),+ $(,)?}) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct $name { pub(crate) inner: ck::$name }
        paste::paste! {
            #[pymethods]
            impl $name {
                #[new]
                #[pyo3(signature=(*, $($field=ck::$name::default().$field,)+ $($selector=None),+))]
                fn new($($field: $ty,)+ $($selector: Option<$input>),+) -> PyResult<Self> {
                    Ok(Self { inner: ck::$name { $($field,)+ $($selector: $selector.map($input::into_rust).transpose()?.unwrap_or_else(|| ck::$name::default().$selector)),+ } })
                }
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
                $(#[getter] fn $selector(&self) -> $output { self.inner.$selector.into() }
                #[setter] fn [<set_ $selector>](&mut self, value: Option<$input>) -> PyResult<()> { self.inner.$selector = value.map($input::into_rust).transpose()?.unwrap_or_else(|| ck::$name::default().$selector); Ok(()) })+
            }
        }
    };
}
io_parameters!(MolBlockWriteParams,
    scalars { force_2d: bool, include_stereo: bool, kekulize: bool, precision: usize, include_coordinates: bool },
    selectors { format: SdfFormatInput => SdfFormat, coordinate_selection: MolCoordinateSelectionInput => MolCoordinateSelection });
io_parameters!(PdbWriteParams, scalars { all_conformers: bool, flavor: u32 },
    selectors { coordinate_selection: MolCoordinateSelectionInput => MolCoordinateSelection });
io_parameters!(Mol2ReadParams, scalars { sanitize: bool, remove_hs: bool, cleanup_substructures: bool },
    selectors { variant: Mol2TypeInput => Mol2Type });

#[pyclass(module = "fixed_projection._native")]
pub struct BatchParams {
    pub(crate) inner: ck::BatchParams,
}
#[pymethods]
impl BatchParams {
    #[new]
    #[pyo3(signature=(*, errors=None, n_jobs=ck::BatchParams::default().n_jobs, progress_bar=ck::BatchParams::default().progress_bar))]
    fn new(
        errors: Option<BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BatchParams {
                errors: errors.map(BatchErrorModeInput::into_rust).transpose()?,
                n_jobs,
                progress_bar,
            },
        })
    }
    #[getter]
    fn errors(&self) -> Option<BatchErrorMode> {
        self.inner.errors.map(Into::into)
    }
    #[setter]
    fn set_errors(&mut self, value: Option<BatchErrorModeInput>) -> PyResult<()> {
        self.inner.errors = value.map(BatchErrorModeInput::into_rust).transpose()?;
        Ok(())
    }
    #[getter]
    fn n_jobs(&self) -> Option<usize> {
        self.inner.n_jobs
    }
    #[setter]
    fn set_n_jobs(&mut self, value: Option<usize>) {
        self.inner.n_jobs = value;
    }
    #[getter]
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar
    }
    #[setter]
    fn set_progress_bar(&mut self, value: Option<bool>) {
        self.inner.progress_bar = value;
    }
}
#[pyclass(module = "fixed_projection._native")]
pub struct BatchExportParams {
    pub(crate) inner: ck::BatchExportParams,
}
#[pymethods]
impl BatchExportParams {
    #[new]
    #[pyo3(signature=(*, format=None, errors=None, n_jobs=ck::BatchExportParams::default().n_jobs, progress_bar=ck::BatchExportParams::default().progress_bar))]
    fn new(
        format: Option<SdfFormatInput>,
        errors: Option<BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BatchExportParams {
                format: format
                    .map(SdfFormatInput::into_rust)
                    .transpose()?
                    .unwrap_or_else(|| ck::BatchExportParams::default().format),
                errors: errors.map(BatchErrorModeInput::into_rust).transpose()?,
                n_jobs,
                progress_bar,
            },
        })
    }
    #[getter]
    fn format(&self) -> SdfFormat {
        self.inner.format.into()
    }
    #[setter]
    fn set_format(&mut self, value: Option<SdfFormatInput>) -> PyResult<()> {
        self.inner.format = value
            .map(SdfFormatInput::into_rust)
            .transpose()?
            .unwrap_or_else(|| ck::BatchExportParams::default().format);
        Ok(())
    }
    #[getter]
    fn errors(&self) -> Option<BatchErrorMode> {
        self.inner.errors.map(Into::into)
    }
    #[setter]
    fn set_errors(&mut self, value: Option<BatchErrorModeInput>) -> PyResult<()> {
        self.inner.errors = value.map(BatchErrorModeInput::into_rust).transpose()?;
        Ok(())
    }
    #[getter]
    fn n_jobs(&self) -> Option<usize> {
        self.inner.n_jobs
    }
    #[setter]
    fn set_n_jobs(&mut self, value: Option<usize>) {
        self.inner.n_jobs = value;
    }
    #[getter]
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar
    }
    #[setter]
    fn set_progress_bar(&mut self, value: Option<bool>) {
        self.inner.progress_bar = value;
    }
}
