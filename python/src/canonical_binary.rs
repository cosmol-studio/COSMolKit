//! Python delegates exclusively to canonical public native serialization.
use ::cosmolkit as ck;
use pyo3::prelude::*;

pyo3::create_exception!(cosmolkit, PickleError, pyo3::exceptions::PyValueError);

pub(crate) fn error_pyerr(py: Python<'_>, source: &ck::PickleError) -> PyErr {
    use ck::PickleError as E;
    let kind = match source {
        E::StereoGroup(_) => "StereoGroup",
        E::UnexpectedEof => "UnexpectedEof",
        E::UnsupportedVersion(_) => "UnsupportedVersion",
        E::UnsupportedArchiveVersion { .. } => "UnsupportedArchiveVersion",
        E::UnsupportedSectionVersion { .. } => "UnsupportedSectionVersion",
        E::InvalidArchive(_) => "InvalidArchive",
        E::MissingRequiredSection(_) => "MissingRequiredSection",
        E::DuplicateSection(_) => "DuplicateSection",
        E::UnknownRequiredSection(_) => "UnknownRequiredSection",
        E::DataLengthMismatch { .. } => "DataLengthMismatch",
        E::InvalidEnumValue { .. } => "InvalidEnumValue",
        E::InvalidMolecule(_) => "InvalidMolecule",
        E::TooManyAtoms(_) => "TooManyAtoms",
        E::TooManyBonds(_) => "TooManyBonds",
        E::StringTooLong(_) => "StringTooLong",
    };
    let error = crate::canonical_values::annotate(
        py,
        PickleError::new_err(source.to_string()),
        "serialization",
        kind,
        source,
    );
    let fields = || -> PyResult<()> {
        let value = error.value(py);
        match source {
            E::UnsupportedVersion(version) => value.setattr("version", *version)?,
            E::UnsupportedArchiveVersion { major, minor } => {
                value.setattr("major", *major)?;
                value.setattr("minor", *minor)?;
            }
            E::UnsupportedSectionVersion { section, version } => {
                value.setattr("section", *section)?;
                value.setattr("version", *version)?;
            }
            E::MissingRequiredSection(section)
            | E::DuplicateSection(section)
            | E::UnknownRequiredSection(section) => value.setattr("section", *section)?,
            E::DataLengthMismatch { expected, actual } => {
                value.setattr("expected", *expected)?;
                value.setattr("actual", *actual)?;
            }
            E::InvalidEnumValue {
                value: invalid,
                type_name,
            } => {
                value.setattr("value", *invalid)?;
                value.setattr("type_name", *type_name)?;
            }
            E::TooManyAtoms(count) | E::TooManyBonds(count) | E::StringTooLong(count) => {
                value.setattr("count", *count)?
            }
            E::InvalidArchive(message) | E::InvalidMolecule(message) => {
                value.setattr("message", message)?
            }
            E::StereoGroup(_) | E::UnexpectedEof => {}
        }
        Ok(())
    };
    match fields() {
        Ok(()) => error,
        Err(error) => error,
    }
}

/// Private pickle reconstruction hook; no chemistry or alternative decoder.
#[pyfunction]
fn _molecule_from_binary(
    py: Python<'_>,
    data: &[u8],
) -> PyResult<crate::drawing_binding::Molecule> {
    ck::Molecule::from_binary(data)
        .map(|inner| crate::drawing_binding::Molecule { inner })
        .map_err(|error| error_pyerr(py, &error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("PickleError", module.py().get_type::<PickleError>())?;
    module.add_function(wrap_pyfunction!(_molecule_from_binary, module)?)?;
    Ok(())
}
