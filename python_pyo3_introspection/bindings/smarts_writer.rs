//! Complete facade SMARTS serialization and concrete native failure transport.
use crate::canonical_errors::{PropertyStringError, RequiredPropertyStringError};
use crate::chemical_values::BondDirection;
use crate::error_framework::Source;
use crate::io_source_errors::ErrorRef;
use crate::model_edit_errors::{AtomPropertyError, MoleculePropertyError};
use crate::model_errors::{CoordinateValidationError, StereoGroupError};
use crate::params::SmartsWriteParams;
use crate::property_numeric_errors::{PropertyIntReadError, PropertyUIntReadError};
use crate::property_values::{PropertyValueError, PropertyValueKind};
use crate::query_compiler::QueryGraph;
use crate::query_model_errors::QueryGraphError;
use crate::query_payloads::{AtomQueryPredicate, BondQueryPredicate, value_kind};
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::smarts_writer_support_errors::{SmartsTraversalError, WedgeError};
use crate::smiles_parse_errors::SmilesParseError;
use crate::stereo_errors::AtropisomerError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use std::error::Error;

value_kind!(SmartsWriteErrorKind, SmartsWriteError, {
    StereoGroup(_), Traversal(_), CxCoordinates(_), CxRingInfo(_), CxRingAtomOrderIndex { index, count },
    CxRingStereoReferenceMissing { bond, side }, CxWedge(_), CxBondConfigAtropMissingCarriers { bond },
    CxAtomPropertyOutput(_), CxMissingConformer, CxCoordinateSource(_), CxCoordinateOutput(_),
    CxSourceBondOutOfRange { bond, bond_count }, CxMoleculePropertyUInt { property, source },
    CxSgroupVectorCast { group, property, actual, expected }, CxSgroupPropertyWrite { group, property, source },
    CxStereoGroup(_), CxSourceAtomOutOfRange { atom, atom_count }, MoleculePropertyWrite(_),
    SourceAtomCount { count }, AtomPropertyWrite { atom, source }, CanonicalTraversal(_),
    UnwritableBondQuery { description }, SourceAtomToLeftIndex { index }, SourceBondBeginIndex { index },
    AtomMapInt { atom, source }, AtomTypeAtomicNumber { atomic_number }, ChargeMagnitudeOverflow { value },
    PropertyValue(_), CxRequiredProperty { atom, source }, CxPropertyList { atom, source },
    CxCoordinateSelectionArity { expected, actual }, CxMissingOutputOrder { template },
    CxOutputOrderPropertyType { template, property }, CxAtomPropertyKind { template, atom, property, kind },
    CxAtomPropertyUInt { template, atom, property, source }, CxAtomPropertyWrite { template, atom, property, source },
    CxCoordinateStorage { template, source }, CxCoordinateSelection { template, source }, CxComposition(_),
    CxOutputOrder { template, kind, row, count }, CxRowCount { kind, count },
    CxBondPropertyUInt { bond, property, source }, InvalidPropertyKind { bond, property, kind },
    CxAtomPropertyInt { atom, property, source }, CxSgroupPropertyUInt { group, property, source },
    Property(_), Valence(_), InvalidGraph(_), QueryGraphTraversalUnsupported { detail },
    OrAboveAndBelowAnd, UnknownCombination { description }, MissingRecursiveQueryMolecule,
    SourceBondDirection { direction }, UnsupportedBondQuery { predicate }, UnsupportedAtomQuery { predicate },
    CompositeChildCount { kind }, XorComposite, QueryGraphCxExtensionsUnsupported { detail },
    RootedAtomOutOfRange { atom }, EmptyAtomSelection, EmptyBondSelection,
    FragmentAtomOutOfRange { atom }, FragmentBondOutOfRange { bond }, BondAtomNotEndpoint { bond, atom },
});
pub(crate) type WriterSource = Source<
    Source<
        Py<StereoGroupError>,
        Source<
            Py<SmartsTraversalError>,
            Source<Py<CoordinateValidationError>, Py<RingFindingError>>,
        >,
    >,
    Source<
        Source<
            Py<WedgeError>,
            Source<
                Py<SmilesParseError>,
                Source<Py<PropertyUIntReadError>, Py<MoleculePropertyError>>,
            >,
        >,
        Source<
            Source<
                Py<AtropisomerError>,
                Source<
                    Py<AtomPropertyError>,
                    Source<Py<PropertyIntReadError>, Py<PropertyValueError>>,
                >,
            >,
            Source<
                Py<RequiredPropertyStringError>,
                Source<Py<QueryGraphError>, Source<Py<PropertyStringError>, Py<ValenceError>>>,
            >,
        >,
    >,
>;

#[pyclass(module="cosmolkit", extends=PyValueError)]
pub struct SmartsWriteError {
    inner: RetainedError<ck::SmartsWriteError>,
}
impl SmartsWriteError {
    pub(crate) fn from_rust(inner: ck::SmartsWriteError) -> PyErr {
        Self::from_retained(RetainedError::new(inner))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Self::from_retained(RetainedError::from_ref(node))
    }
    fn from_retained(inner: RetainedError<ck::SmartsWriteError>) -> PyErr {
        use ck::SmartsWriteError as E;
        let cause = match inner.as_ref() {
            E::StereoGroup(v) => Some(StereoGroupError::from_rust(*v)),
            E::Traversal(v) => Some(SmartsTraversalError::from_rust(v.clone())),
            E::CxCoordinates(v)
            | E::CxCoordinateSource(v)
            | E::CxCoordinateStorage { source: v, .. } => {
                Some(CoordinateValidationError::from_rust(v.clone()))
            }
            E::CxRingInfo(v) => Some(RingFindingError::from_rust(v.clone())),
            E::CxWedge(v) => Some(WedgeError::from_rust(v.clone())),
            E::CxAtomPropertyOutput(v)
            | E::CxCoordinateOutput(v)
            | E::CanonicalTraversal(v)
            | E::CxCoordinateSelection { source: v, .. } => {
                Some(SmilesParseError::from_rust(v.clone()))
            }
            E::CxMoleculePropertyUInt { source: v, .. }
            | E::CxAtomPropertyUInt { source: v, .. }
            | E::CxBondPropertyUInt { source: v, .. }
            | E::CxSgroupPropertyUInt { source: v, .. } => {
                Some(PropertyUIntReadError::from_rust(v.clone()))
            }
            E::CxSgroupPropertyWrite { source: v, .. } | E::MoleculePropertyWrite(v) => {
                Some(MoleculePropertyError::from_rust(v.clone()))
            }
            E::CxStereoGroup(v) => Some(AtropisomerError::from_rust(v.clone())),
            E::AtomPropertyWrite { source: v, .. }
            | E::CxPropertyList { source: v, .. }
            | E::CxAtomPropertyWrite { source: v, .. } => {
                Some(AtomPropertyError::from_rust(v.clone()))
            }
            E::AtomMapInt { source: v, .. } | E::CxAtomPropertyInt { source: v, .. } => {
                Some(PropertyIntReadError::from_rust(v.clone()))
            }
            E::PropertyValue(v) => Some(PropertyValueError::from_rust(v.clone())),
            E::CxRequiredProperty { source: v, .. } => {
                Some(RequiredPropertyStringError::from_rust(v.clone()))
            }
            E::CxComposition(v) => Some(QueryGraphError::from_rust(v.clone())),
            E::Property(v) => Some(PropertyStringError::from_rust(v.clone())),
            E::Valence(v) => Some(ValenceError::from_rust(v.clone())),
            _ => {
                debug_assert!(inner.source().is_none());
                None
            }
        };
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => {
                let error = PyErr::from_value(value.into_bound(py).into_any());
                error.set_cause(py, cause);
                error
            }
            Err(error) => error,
        })
    }
}
#[pymethods]
impl SmartsWriteError {
    #[getter]
    fn domain(&self) -> &'static str {
        "search"
    }
    #[getter]
    fn variant(&self) -> SmartsWriteErrorKind {
        self.inner.as_ref().into()
    }
    #[getter]
    fn source_message(&self) -> Option<String> {
        self.inner.source().map(ToString::to_string)
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<WriterSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn predicate(&self) -> Option<Source<AtomQueryPredicate, BondQueryPredicate>> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::UnsupportedAtomQuery { .. } => {
                Some(Source::First(AtomQueryPredicate {
                    owner: self.inner.clone(),
                }))
            }
            ck::SmartsWriteError::UnsupportedBondQuery { .. } => {
                Some(Source::Second(BondQueryPredicate {
                    owner: self.inner.clone(),
                }))
            }
            _ => None,
        }
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        use ck::SmartsWriteError as E;
        match self.inner.as_ref() {
            E::CxSourceAtomOutOfRange { atom, .. }
            | E::AtomPropertyWrite { atom, .. }
            | E::AtomMapInt { atom, .. }
            | E::CxRequiredProperty { atom, .. }
            | E::CxPropertyList { atom, .. }
            | E::CxAtomPropertyKind { atom, .. }
            | E::CxAtomPropertyUInt { atom, .. }
            | E::CxAtomPropertyWrite { atom, .. }
            | E::CxAtomPropertyInt { atom, .. } => Some(atom.index()),
            E::RootedAtomOutOfRange { atom }
            | E::FragmentAtomOutOfRange { atom }
            | E::BondAtomNotEndpoint { atom, .. } => Some(*atom),
            _ => None,
        }
    }
    #[getter]
    fn bond(&self) -> Option<usize> {
        use ck::SmartsWriteError as E;
        match self.inner.as_ref() {
            E::CxRingStereoReferenceMissing { bond, .. }
            | E::CxBondConfigAtropMissingCarriers { bond }
            | E::CxSourceBondOutOfRange { bond, .. }
            | E::CxBondPropertyUInt { bond, .. }
            | E::InvalidPropertyKind { bond, .. } => Some(bond.index()),
            E::FragmentBondOutOfRange { bond } | E::BondAtomNotEndpoint { bond, .. } => Some(*bond),
            _ => None,
        }
    }
    #[getter]
    fn group(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxSgroupVectorCast { group, .. }
            | ck::SmartsWriteError::CxSgroupPropertyWrite { group, .. }
            | ck::SmartsWriteError::CxSgroupPropertyUInt { group, .. } => Some(group.index()),
            _ => None,
        }
    }
    #[getter]
    fn property(&self) -> Option<&'static str> {
        use ck::SmartsWriteError as E;
        match self.inner.as_ref() {
            E::CxMoleculePropertyUInt { property, .. }
            | E::CxSgroupVectorCast { property, .. }
            | E::CxSgroupPropertyWrite { property, .. }
            | E::CxOutputOrderPropertyType { property, .. }
            | E::CxAtomPropertyKind { property, .. }
            | E::CxAtomPropertyUInt { property, .. }
            | E::CxAtomPropertyWrite { property, .. }
            | E::CxBondPropertyUInt { property, .. }
            | E::InvalidPropertyKind { property, .. }
            | E::CxAtomPropertyInt { property, .. }
            | E::CxSgroupPropertyUInt { property, .. } => Some(property),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<Source<&'static str, PropertyValueKind>> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxOutputOrder { kind, .. }
            | ck::SmartsWriteError::CxRowCount { kind, .. }
            | ck::SmartsWriteError::CompositeChildCount { kind } => Some(Source::First(kind)),
            ck::SmartsWriteError::CxAtomPropertyKind { kind, .. }
            | ck::SmartsWriteError::InvalidPropertyKind { kind, .. } => {
                Some(Source::Second((*kind).into()))
            }
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn actual(&self) -> Option<Source<usize, PropertyValueKind>> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxCoordinateSelectionArity { actual, .. } => {
                Some(Source::First(*actual))
            }
            ck::SmartsWriteError::CxSgroupVectorCast { actual, .. } => {
                Some(Source::Second((*actual).into()))
            }
            _ => None,
        }
    }
    #[getter]
    fn expected(&self) -> Option<Source<usize, &'static str>> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxCoordinateSelectionArity { expected, .. } => {
                Some(Source::First(*expected))
            }
            ck::SmartsWriteError::CxSgroupVectorCast { expected, .. } => {
                Some(Source::Second(expected))
            }
            _ => None,
        }
    }
    #[getter]
    fn template(&self) -> Option<usize> {
        use ck::SmartsWriteError as E;
        match self.inner.as_ref() {
            E::CxMissingOutputOrder { template }
            | E::CxOutputOrderPropertyType { template, .. }
            | E::CxAtomPropertyKind { template, .. }
            | E::CxAtomPropertyUInt { template, .. }
            | E::CxAtomPropertyWrite { template, .. }
            | E::CxCoordinateStorage { template, .. }
            | E::CxCoordinateSelection { template, .. }
            | E::CxOutputOrder { template, .. } => Some(*template),
            _ => None,
        }
    }
    #[getter]
    fn index(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxRingAtomOrderIndex { index, .. }
            | ck::SmartsWriteError::SourceAtomToLeftIndex { index }
            | ck::SmartsWriteError::SourceBondBeginIndex { index } => Some(*index),
            _ => None,
        }
    }
    #[getter]
    fn count(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxRingAtomOrderIndex { count, .. }
            | ck::SmartsWriteError::SourceAtomCount { count }
            | ck::SmartsWriteError::CxOutputOrder { count, .. }
            | ck::SmartsWriteError::CxRowCount { count, .. } => Some(*count),
            _ => None,
        }
    }
    #[getter]
    fn atom_count(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxSourceAtomOutOfRange { atom_count, .. } => Some(*atom_count),
            _ => None,
        }
    }
    #[getter]
    fn bond_count(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxSourceBondOutOfRange { bond_count, .. } => Some(*bond_count),
            _ => None,
        }
    }
    #[getter]
    fn side(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxRingStereoReferenceMissing { side, .. } => Some(*side),
            _ => None,
        }
    }
    #[getter]
    fn row(&self) -> Option<usize> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::CxOutputOrder { row, .. } => Some(*row),
            _ => None,
        }
    }
    #[getter]
    fn atomic_number(&self) -> Option<i32> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::AtomTypeAtomicNumber { atomic_number } => Some(*atomic_number),
            _ => None,
        }
    }
    #[getter]
    fn value(&self) -> Option<Source<i32, String>> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::ChargeMagnitudeOverflow { value } => Some(Source::First(*value)),
            ck::SmartsWriteError::InvalidGraph(value) => Some(Source::Second(value.clone())),
            _ => None,
        }
    }
    #[getter]
    fn description(&self) -> Option<String> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::UnwritableBondQuery { description } => {
                Some((*description).into())
            }
            ck::SmartsWriteError::UnknownCombination { description } => Some(description.clone()),
            _ => None,
        }
    }
    #[getter]
    fn detail(&self) -> Option<&'static str> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::QueryGraphTraversalUnsupported { detail }
            | ck::SmartsWriteError::QueryGraphCxExtensionsUnsupported { detail } => Some(detail),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn direction(&self) -> Option<BondDirection> {
        match self.inner.as_ref() {
            ck::SmartsWriteError::SourceBondDirection { direction } => Some((*direction).into()),
            _ => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

pub(crate) fn keywords(
    include_atom_maps: Option<bool>,
    isomeric_smiles: Option<bool>,
    include_dative_bonds: Option<bool>,
    rooted_at_atom: Option<usize>,
) -> SmartsWriteParams {
    let defaults = ck::SmartsWriteParams::default();
    SmartsWriteParams {
        inner: crate::parameters::NativeField::Owned(ck::SmartsWriteParams {
            include_atom_maps: include_atom_maps.unwrap_or(defaults.include_atom_maps),
            isomeric_smiles: isomeric_smiles.unwrap_or(defaults.isomeric_smiles),
            include_dative_bonds: include_dative_bonds.unwrap_or(defaults.include_dative_bonds),
            rooted_at_atom,
        }),
    }
}
pub(crate) fn text(value: ck::PropertyText) -> PyResult<String> {
    Python::attach(|py| crate::property_values::decode_source_text(py, &value))
}
fn mixed(params: Option<&SmartsWriteParams>, values: [bool; 4]) -> PyResult<()> {
    if params.is_some() && values.into_iter().any(|value| value) {
        return Err(PyTypeError::new_err(
            "Params and keyword configuration are mutually exclusive",
        ));
    }
    Ok(())
}
#[pyfunction]
#[pyo3(signature=(query, params=crate::arguments::Argument::Missing, *, include_atom_maps=crate::arguments::Argument::Missing, isomeric_smiles=crate::arguments::Argument::Missing, include_dative_bonds=crate::arguments::Argument::Missing, rooted_at_atom=crate::arguments::Argument::Missing))]
pub fn write_smarts(
    query: &QueryGraph,
    params: crate::arguments::Argument<Option<PyRef<'_, SmartsWriteParams>>>,
    include_atom_maps: crate::arguments::Argument<Option<bool>>,
    isomeric_smiles: crate::arguments::Argument<Option<bool>>,
    include_dative_bonds: crate::arguments::Argument<Option<bool>>,
    rooted_at_atom: crate::arguments::Argument<Option<usize>>,
) -> PyResult<String> {
    let keyword_form = include_atom_maps.is_supplied()
        || isomeric_smiles.is_supplied()
        || include_dative_bonds.is_supplied()
        || rooted_at_atom.is_supplied();
    crate::arguments::configuration_arguments!(params; include_atom_maps, isomeric_smiles, include_dative_bonds, rooted_at_atom);
    let params = params.as_deref();
    if keyword_form {
        _write_smarts_keywords(
            query,
            include_atom_maps,
            isomeric_smiles,
            include_dative_bonds,
            rooted_at_atom,
        )
    } else {
        _write_smarts_params(query, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"write_smarts"}"#]
#[pyo3(signature=(query, params=None))]
pub fn _write_smarts_params(
    query: &QueryGraph,
    params: Option<&SmartsWriteParams>,
) -> PyResult<String> {
    let params = Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;
    let defaults = ck::SmartsWriteParams::default();
    ck::write_smarts(query.inner(), params.as_ref().unwrap_or(&defaults))
        .map_err(SmartsWriteError::from_rust)
        .and_then(text)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"write_smarts"}"#]
#[pyo3(signature=(query, *, include_atom_maps=None, isomeric_smiles=None, include_dative_bonds=None, rooted_at_atom=None))]
pub fn _write_smarts_keywords(
    query: &QueryGraph,
    include_atom_maps: Option<bool>,
    isomeric_smiles: Option<bool>,
    include_dative_bonds: Option<bool>,
    rooted_at_atom: Option<usize>,
) -> PyResult<String> {
    _write_smarts_params(
        query,
        Some(&keywords(
            include_atom_maps,
            isomeric_smiles,
            include_dative_bonds,
            rooted_at_atom,
        )),
    )
}
#[pyfunction]
#[pyo3(signature=(query, params=crate::arguments::Argument::Missing, *, include_atom_maps=crate::arguments::Argument::Missing, isomeric_smiles=crate::arguments::Argument::Missing, include_dative_bonds=crate::arguments::Argument::Missing, rooted_at_atom=crate::arguments::Argument::Missing))]
pub fn write_cx_smarts(
    query: &QueryGraph,
    params: crate::arguments::Argument<Option<PyRef<'_, SmartsWriteParams>>>,
    include_atom_maps: crate::arguments::Argument<Option<bool>>,
    isomeric_smiles: crate::arguments::Argument<Option<bool>>,
    include_dative_bonds: crate::arguments::Argument<Option<bool>>,
    rooted_at_atom: crate::arguments::Argument<Option<usize>>,
) -> PyResult<String> {
    let keyword_form = include_atom_maps.is_supplied()
        || isomeric_smiles.is_supplied()
        || include_dative_bonds.is_supplied()
        || rooted_at_atom.is_supplied();
    crate::arguments::configuration_arguments!(params; include_atom_maps, isomeric_smiles, include_dative_bonds, rooted_at_atom);
    let params = params.as_deref();
    if keyword_form {
        _write_cx_smarts_keywords(
            query,
            include_atom_maps,
            isomeric_smiles,
            include_dative_bonds,
            rooted_at_atom,
        )
    } else {
        _write_cx_smarts_params(query, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"write_cx_smarts"}"#]
#[pyo3(signature=(query, params=None))]
pub fn _write_cx_smarts_params(
    query: &QueryGraph,
    params: Option<&SmartsWriteParams>,
) -> PyResult<String> {
    let params = Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;
    let defaults = ck::SmartsWriteParams::default();
    ck::write_cx_smarts(query.inner(), params.as_ref().unwrap_or(&defaults))
        .map_err(SmartsWriteError::from_rust)
        .and_then(text)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"write_cx_smarts"}"#]
#[pyo3(signature=(query, *, include_atom_maps=None, isomeric_smiles=None, include_dative_bonds=None, rooted_at_atom=None))]
pub fn _write_cx_smarts_keywords(
    query: &QueryGraph,
    include_atom_maps: Option<bool>,
    isomeric_smiles: Option<bool>,
    include_dative_bonds: Option<bool>,
    rooted_at_atom: Option<usize>,
) -> PyResult<String> {
    _write_cx_smarts_params(
        query,
        Some(&keywords(
            include_atom_maps,
            isomeric_smiles,
            include_dative_bonds,
            rooted_at_atom,
        )),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arguments::Argument::{Missing, Value};
    use pyo3::types::PyDict;

    #[test]
    fn smarts_writers_forward_all_options_and_preserve_rich_native_failure_owners() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native_query = ck::parse_smarts("[C:1]-[O:2]").unwrap();
            let query = QueryGraph::from_rust(native_query.clone());
            let native_mol = ck::Molecule::from_smiles("CO").unwrap();
            let before = native_mol.to_binary().unwrap();
            let mol = Py::new(py, crate::native::Molecule { inner: native_mol })?;
            for flags in 0..8 {
                for root in [None, Some(1)] {
                    let params = keywords(
                        Some(flags & 1 != 0),
                        Some(flags & 2 != 0),
                        Some(flags & 4 != 0),
                        root,
                    );
                    let expected = ck::write_smarts(&native_query, &params.to_rust(py)?).unwrap();
                    let expected = String::from_utf8(expected.into_bytes()).unwrap();
                    assert_eq!(_write_smarts_params(&query, Some(&params))?, expected);
                    assert_eq!(
                        _write_smarts_keywords(
                            &query,
                            Some(flags & 1 != 0),
                            Some(flags & 2 != 0),
                            Some(flags & 4 != 0),
                            root
                        )?,
                        expected
                    );
                    let expected =
                        ck::write_cx_smarts(&native_query, &params.to_rust(py)?).unwrap();
                    let expected = String::from_utf8(expected.into_bytes()).unwrap();
                    assert_eq!(_write_cx_smarts_params(&query, Some(&params))?, expected);
                    assert_eq!(
                        _write_cx_smarts_keywords(
                            &query,
                            Some(flags & 1 != 0),
                            Some(flags & 2 != 0),
                            Some(flags & 4 != 0),
                            root
                        )?,
                        expected
                    );
                    let kwargs = PyDict::new(py);
                    kwargs.set_item("include_atom_maps", flags & 1 != 0)?;
                    kwargs.set_item("isomeric_smiles", flags & 2 != 0)?;
                    kwargs.set_item("include_dative_bonds", flags & 4 != 0)?;
                    kwargs.set_item("rooted_at_atom", root)?;
                    let owned_params = Py::new(py, params)?;
                    for (name, alias, expected) in [
                        (
                            "to_smarts",
                            "to_smarts_with_params",
                            String::from_utf8(
                                mol.borrow(py)
                                    .inner
                                    .to_smarts_with_params(&owned_params.borrow(py).to_rust(py)?)
                                    .unwrap()
                                    .into_bytes(),
                            )
                            .unwrap(),
                        ),
                        (
                            "to_cx_smarts",
                            "to_cx_smarts_with_params",
                            String::from_utf8(
                                mol.borrow(py)
                                    .inner
                                    .to_cx_smarts_with_params(&owned_params.borrow(py).to_rust(py)?)
                                    .unwrap()
                                    .into_bytes(),
                            )
                            .unwrap(),
                        ),
                    ] {
                        assert_eq!(
                            mol.bind(py)
                                .call_method1(name, (&owned_params,))?
                                .extract::<String>()?,
                            expected
                        );
                        assert_eq!(
                            mol.bind(py)
                                .call_method1(alias, (&owned_params,))?
                                .extract::<String>()?,
                            expected
                        );
                        assert_eq!(
                            mol.bind(py)
                                .call_method(name, (), Some(&kwargs))?
                                .extract::<String>()?,
                            expected
                        );
                        assert!(
                            mol.bind(py)
                                .call_method(name, (&owned_params,), Some(&kwargs))
                                .unwrap_err()
                                .is_instance_of::<PyTypeError>(py)
                        );
                    }
                }
            }
            assert_eq!(mol.borrow(py).inner.to_binary().unwrap(), before);
            let params = keywords(None, None, None, Some(999));
            let native = ck::write_smarts(&native_query, &params.to_rust(py)?).unwrap_err();
            let error = _write_smarts_params(&query, Some(&params)).unwrap_err();
            assert!(error.is_instance_of::<SmartsWriteError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, native.to_string());
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 999);
            let params = Py::new(py, params)?;
            assert!(
                write_smarts(
                    &query,
                    Value(Some(params.borrow(py))),
                    Value(Some(true)),
                    Missing,
                    Missing,
                    Missing
                )
                .unwrap_err()
                .is_instance_of::<PyTypeError>(py)
            );

            let mut recursive = ck::RecursiveStructureQuery::from_query_graph(native_query, 17);
            recursive.set_initialized(true);
            let error = SmartsWriteError::from_rust(ck::SmartsWriteError::UnsupportedAtomQuery {
                predicate: ck::AtomQueryPredicate::RecursiveSmarts(recursive),
            });
            let predicate = error.value(py).getattr("predicate")?;
            let recursive = predicate.getattr("value")?;
            let graph = recursive.getattr("query_graph")?;
            assert!(recursive.getattr("initialized")?.extract::<bool>()?);
            assert_eq!(recursive.getattr("serial_number")?.extract::<u32>()?, 17);
            assert!(predicate.setattr("value", 0).is_err());
            assert!(recursive.setattr("serial_number", 0).is_err());
            drop(error);
            drop(predicate);
            drop(recursive);
            assert_eq!(graph.call_method0("num_atoms")?.extract::<usize>()?, 2);

            let bounds = ck::AtomRangeBounds::Inclusive {
                lower: -3,
                upper: 7,
                lower_open: true,
                upper_open: false,
            };
            let data = ck::AtomRangeDataFunction::AtomRingSize {
                lower: 4,
                upper: 9,
                lower_open: false,
                upper_open: true,
            };
            let error = SmartsWriteError::from_rust(ck::SmartsWriteError::UnsupportedAtomQuery {
                predicate: ck::AtomQueryPredicate::Range(ck::AtomRangeQuery::new(bounds, data)),
            });
            let range = error.value(py).getattr("predicate")?.getattr("value")?;
            assert_eq!(
                range
                    .getattr("bounds")?
                    .getattr("lower")?
                    .extract::<i32>()?,
                -3
            );
            assert!(
                range
                    .getattr("bounds")?
                    .getattr("lower_open")?
                    .extract::<bool>()?
            );
            assert_eq!(
                range
                    .getattr("data_function")?
                    .getattr("upper")?
                    .extract::<i32>()?,
                9
            );
            assert!(
                range
                    .getattr("data_function")?
                    .getattr("upper_open")?
                    .extract::<bool>()?
            );
            let error = SmartsWriteError::from_rust(ck::SmartsWriteError::UnsupportedBondQuery {
                predicate: ck::BondQueryPredicate::OrderIn(vec![
                    ck::BondOrder::Single,
                    ck::BondOrder::Double,
                ]),
            });
            assert_eq!(
                error
                    .value(py)
                    .getattr("predicate")?
                    .getattr("value")?
                    .len()?,
                2
            );
            assert!(error.value(py).getattr("__cause__")?.is_none());
            let native =
                ck::SmartsWriteError::CxCoordinates(ck::CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: 2,
                    rows: 1,
                    atom_count: 3,
                });
            let display = native.to_string();
            let error = SmartsWriteError::from_rust(native);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance_of::<CoordinateValidationError>());
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(cause.getattr("atom_count")?.extract::<usize>()?, 3);
            Ok(())
        })
        .unwrap();
    }
}
