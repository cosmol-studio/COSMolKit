//! Canonical IO values and errors; chemistry lives exclusively in cosmolkit.
use ::cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, MolecularIoError, PyValueError);
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct XyzWriteParams {
    pub(crate) inner: ck::XyzWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl XyzWriteParams {
    #[new]
    #[pyo3(signature = (conformer_id=None, precision=6))]
    fn new(conformer_id: Option<usize>, precision: u32) -> Self {
        Self {
            inner: ck::XyzWriteParams {
                conformer_id,
                precision,
            },
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id
    }
    #[getter]
    fn precision(&self) -> u32 {
        self.inner.precision
    }
}
pub(crate) fn error_pyerr(py: Python<'_>, source: ck::MolecularIoError) -> PyErr {
    use ck::MolecularIoError as E;
    if let E::Sdf(error) = source {
        return crate::canonical_sdf::sdf_pyerr(py, error);
    }
    let kind = match &source {
        E::Io { .. } => "Io",
        E::XyzRead(_) => "XyzRead",
        E::XyzWrite(_) => "XyzWrite",
        E::Construction(_) => "Construction",
        E::Mol2Read(_) => "Mol2Read",
        E::Mol2Post(_) => "Mol2Post",
        E::MolWrite(_) => "MolWrite",
        E::Sdf(_) => "Sdf",
        E::NoRecord { .. } => "NoRecord",
        E::OutputUtf8 { .. } => "OutputUtf8",
        E::Parameter { .. } => "Parameter",
    };
    let error = crate::canonical_values::annotate(
        py,
        MolecularIoError::new_err(source.to_string()),
        "io",
        kind,
        &source,
    );
    let fields = || -> PyResult<()> {
        let value = error.value(py);
        match &source {
            E::Io { path, source } => {
                // Transport filesystem path identity with the interpreter's
                // filesystem decoder, including non-UTF8 paths on Unix.
                value.setattr("path", path)?;
                value.setattr("filename", path)?;
                value.setattr("os_error_code", source.raw_os_error())?;
                value.setattr("errno", source.raw_os_error())?;
            }
            E::OutputUtf8 { format, source } => {
                let utf8 = source.utf8_error();
                value.setattr("format", *format)?;
                value.setattr("bytes", PyBytes::new(py, source.as_bytes()))?;
                value.setattr("valid_up_to", utf8.valid_up_to())?;
                value.setattr("error_len", utf8.error_len())?;
            }
            E::Parameter { name, detail } => {
                value.setattr("parameter", *name)?;
                value.setattr("detail", *detail)?;
            }
            E::XyzRead(e) => {
                use ck::XyzReadError as R;
                match e {
                    R::AtomCount { value: text } => {
                        value.setattr("value", text)?;
                    }
                    R::MissingCoordinates { line } => {
                        value.setattr("line", *line)?;
                    }
                    R::Coordinate {
                        value: text, line, ..
                    } => {
                        value.setattr("value", text)?;
                        value.setattr("line", *line)?;
                    }
                    R::AtomSymbol { message } => {
                        value.setattr("message", message)?;
                    }
                    R::EmptyBlock
                    | R::UnexpectedEof
                    | R::Topology(_)
                    | R::Coordinates(_)
                    | R::MoleculeProperty(_) => {}
                }
            }
            E::XyzWrite(ck::XyzWriteError::ConformerNotFound { id }) => {
                value.setattr("conformer_id", *id)?;
            }
            E::NoRecord { format } => {
                value.setattr("format", *format)?;
            }
            E::Mol2Read(_)
            | E::Mol2Post(_)
            | E::MolWrite(_)
            | E::Sdf(_)
            | E::XyzWrite(_)
            | E::Construction(_) => {}
        }
        Ok(())
    };
    if let Err(attribute_error) = fields() {
        return attribute_error;
    }
    if let ck::MolecularIoError::MolWrite(cause) = &source {
        let attrs = || -> PyResult<()> {
            match cause {
                ck::MolWriteError::AmbiguousCoordinates { two_d, three_d } => {
                    error
                        .value(py)
                        .setattr("coordinate_kind", "AmbiguousCoordinates")?;
                    error.value(py).setattr("two_d", *two_d)?;
                    error.value(py).setattr("three_d", *three_d)?;
                }
                ck::MolWriteError::MissingCoordinate { dimension, id } => {
                    error
                        .value(py)
                        .setattr("coordinate_kind", "MissingCoordinate")?;
                    error.value(py).setattr(
                        "dimension",
                        match dimension {
                            ck::CoordinateDimension::TwoD => "2d",
                            ck::CoordinateDimension::ThreeD => "3d",
                        },
                    )?;
                    error.value(py).setattr("coordinate_id", *id)?;
                }
                ck::MolWriteError::CoordinateDimensionMismatch => {
                    error
                        .value(py)
                        .setattr("coordinate_kind", "CoordinateDimensionMismatch")?;
                }
                _ => {}
            }
            Ok(())
        };
        if let Err(e) = attrs() {
            return e;
        }
    }
    if let ck::MolecularIoError::Mol2Post(source) = &source {
        if let Err(e) = error.value(py).setattr("stage", source.stage) {
            return e;
        }
    }
    if let ck::MolecularIoError::Construction(cause) = source {
        error.set_cause(py, Some(crate::drawing_binding::operation_pyerr(py, cause)));
    }
    error
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "MolecularIoError",
        module.py().get_type::<MolecularIoError>(),
    )?;
    module.add_class::<XyzWriteParams>()?;
    module.add_class::<Mol2Type>()?;
    module.add_class::<Mol2ReadParams>()?;
    module.add_class::<SdfFormat>()?;
    module.add_class::<MolCoordinateSelection>()?;
    module.add_class::<MolBlockWriteParams>()?;
    Ok(())
}

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Mol2Type {
    Corina,
}
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum SdfFormat {
    V2000,
    V3000,
}
impl From<SdfFormat> for ck::SdfFormat {
    fn from(value: SdfFormat) -> Self {
        match value {
            SdfFormat::V2000 => Self::V2000,
            SdfFormat::V3000 => Self::V3000,
        }
    }
}
impl From<ck::SdfFormat> for SdfFormat {
    fn from(value: ck::SdfFormat) -> Self {
        match value {
            ck::SdfFormat::V2000 => Self::V2000,
            ck::SdfFormat::V3000 => Self::V3000,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct Mol2ReadParams {
    pub(crate) inner: ck::Mol2ReadParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Mol2ReadParams {
    #[new]
    #[pyo3(signature = (*, sanitize=true, remove_hs=true, variant=Mol2Type::Corina, cleanup_substructures=true))]
    fn new(
        sanitize: bool,
        remove_hs: bool,
        variant: Mol2Type,
        cleanup_substructures: bool,
    ) -> Self {
        let variant = match variant {
            Mol2Type::Corina => ck::Mol2Type::Corina,
        };
        Self {
            inner: ck::Mol2ReadParams {
                sanitize,
                remove_hs,
                variant,
                cleanup_substructures,
            },
        }
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[getter]
    fn variant(&self) -> Mol2Type {
        match self.inner.variant {
            ck::Mol2Type::Corina => Mol2Type::Corina,
        }
    }
    #[getter]
    fn cleanup_substructures(&self) -> bool {
        self.inner.cleanup_substructures
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct MolCoordinateSelection {
    pub(crate) inner: ck::MolCoordinateSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MolCoordinateSelection {
    #[new]
    fn new() -> Self {
        Self {
            inner: ck::MolCoordinateSelection::Auto,
        }
    }
    #[staticmethod]
    fn auto() -> Self {
        Self::new()
    }
    #[staticmethod]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::MolCoordinateSelection::TwoD { id },
        }
    }
    #[staticmethod]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::MolCoordinateSelection::ThreeD { id },
        }
    }
    #[getter]
    fn dimension(&self) -> Option<crate::canonical_sdf::CoordinateDimension> {
        match self.inner {
            ck::MolCoordinateSelection::Auto => None,
            ck::MolCoordinateSelection::TwoD { .. } => {
                Some(crate::canonical_sdf::CoordinateDimension::TwoD)
            }
            ck::MolCoordinateSelection::ThreeD { .. } => {
                Some(crate::canonical_sdf::CoordinateDimension::ThreeD)
            }
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        match self.inner {
            ck::MolCoordinateSelection::Auto => None,
            ck::MolCoordinateSelection::TwoD { id } | ck::MolCoordinateSelection::ThreeD { id } => {
                Some(id)
            }
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct MolBlockWriteParams {
    pub(crate) inner: ck::MolBlockWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MolBlockWriteParams {
    #[new]
    #[pyo3(signature = (*, format=SdfFormat::V2000, force_2d=false, include_stereo=true, kekulize=true, precision=6, coordinate_selection=None, include_coordinates=true))]
    fn new(
        format: SdfFormat,
        force_2d: bool,
        include_stereo: bool,
        kekulize: bool,
        precision: usize,
        coordinate_selection: Option<&MolCoordinateSelection>,
        include_coordinates: bool,
    ) -> Self {
        Self {
            inner: ck::MolBlockWriteParams {
                format: format.into(),
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection: coordinate_selection
                    .map_or(ck::MolCoordinateSelection::Auto, |value| value.inner),
                include_coordinates,
            },
        }
    }
    #[getter]
    fn format(&self) -> SdfFormat {
        self.inner.format.into()
    }
    #[getter]
    fn force_2d(&self) -> bool {
        self.inner.force_2d
    }
    #[getter]
    fn include_stereo(&self) -> bool {
        self.inner.include_stereo
    }
    #[getter]
    fn kekulize(&self) -> bool {
        self.inner.kekulize
    }
    #[getter]
    fn precision(&self) -> usize {
        self.inner.precision
    }
    #[getter]
    fn coordinate_selection(&self) -> MolCoordinateSelection {
        MolCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
    #[getter]
    fn include_coordinates(&self) -> bool {
        self.inner.include_coordinates
    }
}
