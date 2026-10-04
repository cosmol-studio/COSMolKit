//! Experimental drawing projections; all chemistry stays in the public facade.

use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, DrawingError, PyValueError);

pyo3::create_exception!(cosmolkit, OperationError, PyValueError);

// Transport actual source messages only; Python cannot retain Rust downcast identity.
fn source_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    let error = PyValueError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|cause| source_pyerr(py, cause)));
    error
}

fn drawing_pyerr(py: Python<'_>, source: ck::DrawingError) -> PyErr {
    use ck::DrawingError as E;
    let error = DrawingError::new_err(source.to_string());
    let kind = match &source {
        E::Property(..) => "Property",
        E::Topology(..) => "Topology",
        E::Coordinates(..) => "Coordinates",
        E::Mapping(..) => "Mapping",
        E::Kekulize(..) => "Kekulize",
        E::Hydrogen(..) => "Hydrogen",
        E::Wedge(..) => "Wedge",
        E::Valence(..) => "Valence",
        E::CoordinateGeneration(..) => "CoordinateGeneration",
        E::SvgParse(..) => "SvgParse",
        E::PngEncode(..) => "PngEncode",
        E::StateRows { .. } => "StateRows",
        E::HydrogenAppend { .. } => "HydrogenAppend",
        E::InvalidDimensions { .. } => "InvalidDimensions",
        E::PixmapAllocation { .. } => "PixmapAllocation",
    };
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        value.setattr("domain", "drawing")?;
        value.setattr("kind", kind)?;
        match &source {
            E::InvalidDimensions { width, height } | E::PixmapAllocation { width, height } => {
                value.setattr("width", *width)?;
                value.setattr("height", *height)?;
            }
            E::StateRows {
                field,
                actual,
                expected,
            } => {
                value.setattr("field", *field)?;
                value.setattr("actual", *actual)?;
                value.setattr("expected", *expected)?;
            }
            E::HydrogenAppend { row, reason } => {
                value.setattr("row", *row)?;
                value.setattr("reason", *reason)?;
            }
            _ => {}
        }
        Ok(())
    };
    if let Err(attribute_error) = attributes() {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source).map(|cause| source_pyerr(py, cause)),
    );
    error
}

fn operation_pyerr(py: Python<'_>, source: ck::OperationError) -> PyErr {
    use ck::OperationError as E;
    let kind = match &source {
        E::UnsupportedFeature { .. } => "UnsupportedFeature",
        E::Unsupported { .. } => "Unsupported",
        E::OutputMismatch { .. } => "OutputMismatch",
        E::AccessDenied { .. } => "AccessDenied",
        E::BlockCheckedOut { .. } => "BlockCheckedOut",
        E::BlockNotCheckedOut { .. } => "BlockNotCheckedOut",
        E::IncompleteCommit { .. } => "IncompleteCommit",
        E::TopologyEditContract { .. } => "TopologyEditContract",
        E::MappingContract { .. } => "MappingContract",
        E::InvalidTopologyMapping { .. } => "InvalidTopologyMapping",
        E::AutoRemapContract { .. } => "AutoRemapContract",
        E::OperationContract { .. } => "OperationContract",
        E::SemanticPreconditionContract { .. } => "SemanticPreconditionContract",
        E::CoordinateAppendRequiresValues { .. } => "CoordinateAppendRequiresValues",
        E::DerivedEffectContract { .. } => "DerivedEffectContract",
        E::CipStateContract { .. } => "CipStateContract",
        E::InvalidPropertyList { .. } => "InvalidPropertyList",
        E::InvalidDerivedCache { .. } => "InvalidDerivedCache",
        E::InvalidAlgorithmResult { .. } => "InvalidAlgorithmResult",
        E::Algorithm { .. } => "Algorithm",
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidTopologyEdit(..) => "InvalidTopologyEdit",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::InvalidProperty(..) => "InvalidProperty",
        E::Valence(..) => "Valence",
        E::Radical(..) => "Radical",
        E::Rings(..) => "Rings",
        E::PotentialStereo(..) => "PotentialStereo",
        E::Stereo(..) => "Stereo",
        E::CipLabeler(..) => "CipLabeler",
        E::Transform(..) => "Transform",
        E::Coordinate2D(..) => "Coordinate2D",
        E::Kekulize(..) => "Kekulize",
        E::Aromaticity(..) => "Aromaticity",
        E::Sanitize(..) => "Sanitize",
        E::Hydrogen(..) => "Hydrogen",
    };
    let error = OperationError::new_err(source.to_string());
    if let Err(attribute_error) = error
        .value(py)
        .setattr("domain", "operation")
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source).map(|cause| source_pyerr(py, cause)),
    );
    error
}

fn expand_user_path(path: &str) -> PyResult<std::path::PathBuf> {
    if path == "~" || path.starts_with("~/") {
        let home = std::env::var_os("HOME")
            .ok_or_else(|| PyValueError::new_err("cannot expand '~': HOME is not set"))?;
        let mut expanded = std::path::PathBuf::from(home);
        if let Some(rest) = path.strip_prefix("~/") {
            expanded.push(rest);
        }
        Ok(expanded)
    } else {
        Ok(std::path::PathBuf::from(path))
    }
}

fn write_drawing_file(path: &str, bytes: &[u8]) -> PyResult<()> {
    let expanded = expand_user_path(path)?;
    std::fs::write(&expanded, bytes).map_err(|error| {
        pyo3::exceptions::PyOSError::new_err((
            error.raw_os_error(),
            error.to_string(),
            expanded.to_string_lossy().into_owned(),
        ))
    })
}

/// Immutable detached parameters projected from the public facade.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
struct Coordinate2DParams {
    inner: ck::Coordinate2DParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Coordinate2DParams {
    #[new]
    #[pyo3(signature = (coordinate_map=None, *, canonical_orientation=false, clear_existing_2d=true, flips_per_sample=0, samples=0, sample_seed=0, permute_degree_four=false, force_rdkit=false, use_ring_templates=false))]
    fn new(
        coordinate_map: Option<std::collections::BTreeMap<usize, [f64; 2]>>,
        canonical_orientation: bool,
        clear_existing_2d: bool,
        flips_per_sample: u32,
        samples: u32,
        sample_seed: i32,
        permute_degree_four: bool,
        force_rdkit: bool,
        use_ring_templates: bool,
    ) -> Self {
        Self {
            inner: ck::Coordinate2DParams {
                coordinate_map: coordinate_map.unwrap_or_default(),
                canonical_orientation,
                clear_existing_2d,
                flips_per_sample,
                samples,
                sample_seed,
                permute_degree_four,
                force_rdkit,
                use_ring_templates,
            },
        }
    }

    #[getter]
    fn coordinate_map(&self) -> std::collections::BTreeMap<usize, [f64; 2]> {
        self.inner.coordinate_map.clone()
    }

    #[getter]
    fn canonical_orientation(&self) -> bool {
        self.inner.canonical_orientation
    }

    #[getter]
    fn clear_existing_2d(&self) -> bool {
        self.inner.clear_existing_2d
    }

    #[getter]
    fn flips_per_sample(&self) -> u32 {
        self.inner.flips_per_sample
    }

    #[getter]
    fn samples(&self) -> u32 {
        self.inner.samples
    }

    #[getter]
    fn sample_seed(&self) -> i32 {
        self.inner.sample_seed
    }

    #[getter]
    fn permute_degree_four(&self) -> bool {
        self.inner.permute_degree_four
    }

    #[getter]
    fn force_rdkit(&self) -> bool {
        self.inner.force_rdkit
    }

    #[getter]
    fn use_ring_templates(&self) -> bool {
        self.inner.use_ring_templates
    }
}

/// Python ownership wraps the ONE live runtime value, not detached chemistry.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
struct Molecule {
    inner: ck::Molecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Molecule {
    #[staticmethod]
    fn from_smiles(smiles: &str) -> PyResult<Self> {
        ck::Molecule::from_smiles(smiles)
            .map(|inner| Self { inner })
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }

    fn to_smiles(&self) -> PyResult<String> {
        self.inner
            .to_smiles()
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn coordinates_2d(&self) -> Option<Vec<[f64; 2]>> {
        // Python receives an owned copy and cannot mutate the runtime block.
        self.inner.coordinates_2d().map(<[_]>::to_vec)
    }

    fn with_2d_coordinates(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_2d_coordinates()
            .map(|inner| Self { inner })
            .map_err(|error| operation_pyerr(py, error))
    }

    fn with_2d_coordinates_with_params(
        &self,
        py: Python<'_>,
        params: &Coordinate2DParams,
    ) -> PyResult<Self> {
        self.inner
            .with_2d_coordinates_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|error| operation_pyerr(py, error))
    }

    /// Experimental SVG; required dimensions follow the current registry.
    #[pyo3(signature = (width, height))]
    fn to_svg(&self, py: Python<'_>, width: u32, height: u32) -> PyResult<String> {
        self.inner
            .to_svg(width, height)
            .map_err(|error| drawing_pyerr(py, error))
    }

    /// Experimental PNG bytes from the same public drawing owner.
    #[gen_stub(override_return_type(type_repr = "builtins.bytes", imports = ("builtins")))]
    #[pyo3(signature = (width, height))]
    fn to_png<'py>(
        &self,
        py: Python<'py>,
        width: u32,
        height: u32,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let png = self
            .inner
            .to_png(width, height)
            .map_err(|error| drawing_pyerr(py, error))?;
        Ok(PyBytes::new(py, &png))
    }

    #[pyo3(signature = (path, width, height))]
    fn write_svg(&self, py: Python<'_>, path: &str, width: u32, height: u32) -> PyResult<()> {
        let svg = self.to_svg(py, width, height)?;
        write_drawing_file(path, svg.as_bytes())
    }

    #[pyo3(signature = (path, width, height))]
    fn write_png(&self, py: Python<'_>, path: &str, width: u32, height: u32) -> PyResult<()> {
        let png = self.to_png(py, width, height)?;
        write_drawing_file(path, png.as_bytes())
    }
}

#[pymodule]
fn cosmolkit(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", ck::version())?;
    module.add("_binding_profile", "drawing-bindings")?;
    module.add("DrawingError", module.py().get_type::<DrawingError>())?;
    module.add("OperationError", module.py().get_type::<OperationError>())?;
    module.add_class::<Molecule>()?;
    module.add_class::<Coordinate2DParams>()?;
    Ok(())
}
