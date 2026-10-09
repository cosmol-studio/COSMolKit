//! Concrete registered source errors, preserved through native error chains.
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

macro_rules! source_errors {
    ($($name:ident),* $(,)?) => {
        $(pyo3::create_exception!(cosmolkit, $name, PyValueError);)*
        pub(crate) fn register(module: &Bound<'_, pyo3::types::PyModule>) -> PyResult<()> {
            $(module.add(stringify!($name), module.py().get_type::<$name>())?;)*
            Ok(())
        }
        pub(crate) fn convert(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> Option<PyErr> {
            $(if let Some(source) = source.downcast_ref::<ck::$name>() {
                let error = $name::new_err(source.to_string());
                error.set_cause(py, std::error::Error::source(source).map(|cause| crate::canonical_values::source_pyerr(py, cause)));
                return Some(error);
            })*
            None
        }
    }
}

source_errors!(
    BioSelectionMatchError,
    BioSelectionCopyError,
    BioSelectionCopyCause,
    BioRowTraverseError,
    BioRowModelError,
    BioRowChainError,
    TemplateAttachmentOrderError,
    Coordinate2DError,
    Coordinate2DTemplateError,
    Coordinate2DLayoutError,
    TransformError,
    AromaticityError,
    StereoError,
    CipLabelerError,
    HydrogenError,
    ResidueSequenceError,
    ProteinProjectionError,
    SdfReadError,
    SmilesStereoError,
    Mol2ReadError,
    Mol2PostError,
    XyzReadError,
    XyzWriteError,
    MolWriteError,
);
