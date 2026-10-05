//! Thin canonical native hash and legacy CIP error projections.
use ::cosmolkit as ck;
use pyo3::prelude::*;
pyo3::create_exception!(cosmolkit, MoleculeHashError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(cosmolkit, CipRankError, pyo3::exceptions::PyValueError);
fn cip_error_pyerr(py: Python<'_>, source: &ck::CipRankError) -> PyErr {
    use ck::CipRankError as E;
    let kind = match source {
        E::InvalidTopology(_) => "InvalidTopology",
        E::InvalidQueryState(_) => "InvalidQueryState",
        E::ValenceRowCount { .. } => "ValenceRowCount",
        E::NegativeImplicitHydrogen { .. } => "NegativeImplicitHydrogen",
        E::AtomMapOutOfRange { .. } => "AtomMapOutOfRange",
        E::InvariantCount { .. } => "InvariantCount",
        E::InvariantOutOfRange { .. } => "InvariantOutOfRange",
        E::TooManyNeighbors { .. } => "TooManyNeighbors",
        E::UnsupportedBondOrder { .. } => "UnsupportedBondOrder",
    };
    let error = crate::canonical_values::annotate(
        py,
        CipRankError::new_err(source.to_string()),
        "cip_ranking",
        kind,
        source,
    );
    let fields = || -> PyResult<()> {
        let value = error.value(py);
        match source {
            E::ValenceRowCount {
                field,
                actual,
                atom_count,
            } => {
                value.setattr("field", *field)?;
                value.setattr("actual", *actual)?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::NegativeImplicitHydrogen {
                atom,
                value: invalid,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("value", *invalid)?;
            }
            E::AtomMapOutOfRange { atom, map_number } => {
                value.setattr("atom", atom.index())?;
                value.setattr("map_number", *map_number)?;
            }
            E::InvariantCount { actual, atom_count } => {
                value.setattr("actual", *actual)?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::InvariantOutOfRange {
                atom,
                value: invalid,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("value", *invalid)?;
            }
            E::TooManyNeighbors {
                atom,
                degree,
                maximum_supported,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("degree", *degree)?;
                value.setattr("maximum_supported", *maximum_supported)?;
            }
            E::UnsupportedBondOrder { bond, order } => {
                value.setattr("bond", bond.index())?;
                value.setattr("order", order.rdkit_code())?;
            }
            E::InvalidTopology(cause) => {
                error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
            }
            E::InvalidQueryState(cause) => {
                error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
            }
        }
        Ok(())
    };
    match fields() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pub(crate) fn error_pyerr(py: Python<'_>, source: &ck::MoleculeHashError) -> PyErr {
    use ck::MoleculeHashError as E;
    let kind = match source {
        E::EmptyMolecule => "EmptyMolecule",
        E::MissingPreparedValence => "MissingPreparedValence",
        E::CipRanks(_) => "CipRanks",
        E::InvalidTopology(_) => "InvalidTopology",
        E::RankCount { .. } => "RankCount",
    };
    let error = crate::canonical_values::annotate(
        py,
        MoleculeHashError::new_err(source.to_string()),
        "molecular_hash",
        kind,
        source,
    );
    let fields = || -> PyResult<()> {
        match source {
            E::RankCount { actual, atom_count } => {
                error.value(py).setattr("actual", *actual)?;
                error.value(py).setattr("atom_count", *atom_count)?;
            }
            E::CipRanks(cause) => error.set_cause(py, Some(cip_error_pyerr(py, cause))),
            E::EmptyMolecule | E::MissingPreparedValence | E::InvalidTopology(_) => {}
        }
        Ok(())
    };
    match fields() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "MoleculeHashError",
        module.py().get_type::<MoleculeHashError>(),
    )?;
    module.add("CipRankError", module.py().get_type::<CipRankError>())?;
    Ok(())
}
