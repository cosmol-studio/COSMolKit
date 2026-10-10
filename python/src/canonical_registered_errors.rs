//! Concrete registered source errors, preserved through native error chains.
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

macro_rules! source_errors {
    ($($name:ident => $doc:literal),* $(,)?) => {
        $(pyo3::create_exception!(cosmolkit, $name, PyValueError, $doc);)*
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
    BioSelectionMatchError => "A structural selection could not be evaluated against the biological hierarchy.",
    BioSelectionCopyError => "Selected biological rows could not be copied into a consistent structure.",
    BioSelectionCopyCause => "Underlying hierarchy or reference error encountered while copying a structural selection.",
    BioRowTraverseError => "The biological hierarchy could not be traversed with the supplied row references.",
    BioRowModelError => "A biological model row reference is invalid.",
    BioRowChainError => "A biological chain row reference is invalid.",
    TemplateAttachmentOrderError => "Reaction-template attachment order metadata is invalid.",
    Coordinate2DError => "2D coordinate generation or assignment failed.",
    Coordinate2DTemplateError => "A 2D depiction template could not be used for the supplied molecule.",
    Coordinate2DLayoutError => "The molecular 2D layout could not be computed.",
    TransformError => "A coordinate transform could not be applied to the requested coordinates.",
    AromaticityError => "Aromaticity assignment failed for the graph or selected aromaticity model.",
    StereoError => "Stereochemical assignment or manipulation failed.",
    CipLabelerError => "CIP stereochemical labeling failed or exceeded the requested iteration limit.",
    HydrogenError => "Explicit hydrogen addition or removal failed for the supplied graph and options.",
    ResidueSequenceError => "Residue sequence conversion failed for the selected residue classification.",
    ProteinProjectionError => "The biological structure could not be projected as the requested protein view.",
    SdfReadError => "An SDF record or property list could not be read consistently.",
    SmilesStereoError => "SMILES stereochemical annotations could not be interpreted consistently.",
    Mol2ReadError => "MOL2 molecular text could not be parsed.",
    Mol2PostError => "A parsed MOL2 graph could not be chemically prepared as requested.",
    XyzReadError => "XYZ coordinate text could not be parsed into a molecular coordinate graph.",
    XyzWriteError => "The selected molecular coordinates could not be serialized as XYZ.",
    MolWriteError => "The molecular graph could not be serialized as a MOL/SDF block with the requested coordinates and options.",
);
