//! Concrete fragment copy failures retain the native transparent source graph.
use crate::canonical_errors::PropertyStringError;
use crate::chemical_values::enum_value;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::path_errors::{PathError, PathSource};
use crate::property_numeric_errors::{PropertyUIntReadError, UIntLexicalReadError};
use crate::sanitize_errors::{SanitizeError, SanitizeSource};
use crate::topology_edit_errors::{EditSource, TopologyEditError};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

enum_value!(SubstanceGroupAtomList, { Atoms, ParentAtoms });

error_projection!(SelectedBondCopyError, SelectedBondCopyErrorKind, ck::SelectedBondCopyError, "fragments", {
    Stereo(_), MissingEndpointMapping {
        bond: usize => |v: &ck::BondId| v.index(), atom: usize => |v: &ck::AtomId| v.index()
    },
}, cause = |inner: &ck::SelectedBondCopyError| match inner {
    ck::SelectedBondCopyError::Stereo(value) => Some(BondValueError::from_rust(value.clone())),
    ck::SelectedBondCopyError::MissingEndpointMapping { .. } => None,
}, source = Option<Py<BondValueError>>);

error_projection!(SelectedSubstanceGroupCopyError, SelectedSubstanceGroupCopyErrorKind, ck::SelectedSubstanceGroupCopyError, "fragments", {
    MissingAtomMapping { group, list, atom: usize => |v: &ck::AtomId| v.index() },
    AtomOutOfRange { group, list,
        source_atom: usize => |v: &ck::AtomId| v.index(), mapped_atom: usize => |v: &ck::AtomId| v.index(), atom_count: usize },
    MissingBondMapping { group, bond: usize => |v: &ck::BondId| v.index() },
    BondOutOfRange { group, source_bond: usize => |v: &ck::BondId| v.index(),
        mapped_bond: usize => |v: &ck::BondId| v.index(), bond_count: usize },
}, cause = |_: &ck::SelectedSubstanceGroupCopyError| None::<PyErr>, methods = {
    #[getter]
    fn group(&self) -> usize { match self.inner {
        ck::SelectedSubstanceGroupCopyError::MissingAtomMapping { group, .. }
        | ck::SelectedSubstanceGroupCopyError::AtomOutOfRange { group, .. }
        | ck::SelectedSubstanceGroupCopyError::MissingBondMapping { group, .. }
        | ck::SelectedSubstanceGroupCopyError::BondOutOfRange { group, .. } => group.index(),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn list(&self) -> Option<SubstanceGroupAtomList> { match self.inner {
        ck::SelectedSubstanceGroupCopyError::MissingAtomMapping { list, .. }
        | ck::SelectedSubstanceGroupCopyError::AtomOutOfRange { list, .. } => Some(list.into()), _ => None,
    } }
});

type FullSource = Source<EditSource, Py<UIntLexicalReadError>>;
error_projection!(FullCopyComponentError, FullCopyComponentErrorKind, ck::FullCopyComponentError, "fragments", {
    SourceUInt(_), SourceText(_), TopologyEdit(_), MappingValidation(_), CoordinateValidation(_),
    AtomMaskSize { expected: usize, actual: usize }, AddedAtomRow { row: usize },
    MissingCoordinateRow { conformer_id: usize, source_row: usize },
}, cause = |inner: &ck::FullCopyComponentError| Python::attach(|py| match inner {
    ck::FullCopyComponentError::SourceUInt(value) => PropertyUIntReadError::from_rust(value.clone()).cause(py),
    ck::FullCopyComponentError::SourceText(value) => PropertyStringError::from_rust(*value).cause(py),
    ck::FullCopyComponentError::TopologyEdit(value) => TopologyEditError::from_rust(value.clone()).cause(py),
    ck::FullCopyComponentError::MappingValidation(value) => MappingValidationError::from_rust(value.clone()).cause(py),
    ck::FullCopyComponentError::CoordinateValidation(value) => CoordinateValidationError::from_rust(value.clone()).cause(py),
    ck::FullCopyComponentError::AtomMaskSize { .. } | ck::FullCopyComponentError::AddedAtomRow { .. }
    | ck::FullCopyComponentError::MissingCoordinateRow { .. } => None,
}), source = Option<FullSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn source_uint(&self, py: Python<'_>) -> PyResult<Option<Py<PropertyUIntReadError>>> { match &self.inner {
        ck::FullCopyComponentError::SourceUInt(value) => PropertyUIntReadError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn source_text(&self, py: Python<'_>) -> PyResult<Option<Py<PropertyStringError>>> { match &self.inner {
        ck::FullCopyComponentError::SourceText(value) => PropertyStringError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology_edit(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyEditError>>> { match &self.inner {
        ck::FullCopyComponentError::TopologyEdit(value) => TopologyEditError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn mapping_validation(&self, py: Python<'_>) -> PyResult<Option<Py<MappingValidationError>>> { match &self.inner {
        ck::FullCopyComponentError::MappingValidation(value) => MappingValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinate_validation(&self, py: Python<'_>) -> PyResult<Option<Py<CoordinateValidationError>>> { match &self.inner {
        ck::FullCopyComponentError::CoordinateValidation(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

type SubsetSource = Source<
    SanitizeSource,
    Source<EditSource, Source<Py<BondValueError>, Py<MoleculePropertyError>>>,
>;
error_projection!(AtomPathSubsetCopyError, AtomPathSubsetCopyErrorKind, ck::AtomPathSubsetCopyError, "fragments", {
    StereoGroup(_), ComputedProperties(_), MoleculeProperty(_), BondProperty(_), AtomProperty(_),
    BondCopy(_), SubstanceGroupCopy(_), TopologyValidation(_), MappingValidation(_), Sanitize(_),
}, cause = |inner: &ck::AtomPathSubsetCopyError| Python::attach(|py| match inner {
    ck::AtomPathSubsetCopyError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::AtomPathSubsetCopyError::ComputedProperties(value) => TopologyEditError::from_rust(value.clone()).cause(py),
    ck::AtomPathSubsetCopyError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::AtomPathSubsetCopyError::BondProperty(value) => Some(BondValueError::from_rust(value.clone())),
    ck::AtomPathSubsetCopyError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::AtomPathSubsetCopyError::BondCopy(value) => SelectedBondCopyError::from_rust(value.clone()).cause(py),
    ck::AtomPathSubsetCopyError::SubstanceGroupCopy(value) => SelectedSubstanceGroupCopyError::from_rust(value.clone()).cause(py),
    ck::AtomPathSubsetCopyError::TopologyValidation(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::AtomPathSubsetCopyError::MappingValidation(value) => MappingValidationError::from_rust(value.clone()).cause(py),
    ck::AtomPathSubsetCopyError::Sanitize(value) => SanitizeError::from_rust(value.clone()).cause(py),
}), source = Option<SubsetSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn computed_properties(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyEditError>>> { match &self.inner {
        ck::AtomPathSubsetCopyError::ComputedProperties(value) => TopologyEditError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn bond_copy(&self, py: Python<'_>) -> PyResult<Option<Py<SelectedBondCopyError>>> { match &self.inner {
        ck::AtomPathSubsetCopyError::BondCopy(value) => SelectedBondCopyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn substance_group_copy(&self, py: Python<'_>) -> PyResult<Option<Py<SelectedSubstanceGroupCopyError>>> { match &self.inner {
        ck::AtomPathSubsetCopyError::SubstanceGroupCopy(value) => SelectedSubstanceGroupCopyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology_validation(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::AtomPathSubsetCopyError::TopologyValidation(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn mapping_validation(&self, py: Python<'_>) -> PyResult<Option<Py<MappingValidationError>>> { match &self.inner {
        ck::AtomPathSubsetCopyError::MappingValidation(value) => MappingValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn sanitize(&self, py: Python<'_>) -> PyResult<Option<Py<SanitizeError>>> { match &self.inner {
        ck::AtomPathSubsetCopyError::Sanitize(value) => SanitizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

type BuildSource =
    Source<PathSource, Source<Py<AtomPathSubsetCopyError>, Py<FullCopyComponentError>>>;
error_projection!(OrderedFragmentBuildError, OrderedFragmentBuildErrorKind, ck::OrderedFragmentBuildError, "fragments", {
    ComponentLabels(_), FastSubset { component_index, source }, SlowFullCopy { component_index, source },
}, cause = |inner: &ck::OrderedFragmentBuildError| Python::attach(|py| match inner {
    ck::OrderedFragmentBuildError::ComponentLabels(value) => PathError::from_rust(value.clone()).cause(py),
    ck::OrderedFragmentBuildError::FastSubset { source, .. } => Some(AtomPathSubsetCopyError::from_rust(source.clone())),
    ck::OrderedFragmentBuildError::SlowFullCopy { source, .. } => Some(FullCopyComponentError::from_rust(source.clone())),
}), source = Option<BuildSource>, methods = {
    #[getter]
    fn component_index(&self) -> Option<usize> { match self.inner {
        ck::OrderedFragmentBuildError::FastSubset { component_index, .. }
        | ck::OrderedFragmentBuildError::SlowFullCopy { component_index, .. } => Some(component_index), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn component_labels(&self, py: Python<'_>) -> PyResult<Option<Py<PathError>>> { match &self.inner {
        ck::OrderedFragmentBuildError::ComponentLabels(value) => PathError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

type FailureSource = Source<BuildSource, Source<Py<MoleculePropertyError>, Py<SanitizeError>>>;
error_projection!(MoleculeFragmentsFailure, MoleculeFragmentsFailureKind, ck::MoleculeFragmentsFailure, "fragments", {
    FinalComputedProperties { component_index, source }, Build(_), FinalSanitize { component_index, source },
}, cause = |inner: &ck::MoleculeFragmentsFailure| Python::attach(|py| match inner {
    ck::MoleculeFragmentsFailure::FinalComputedProperties { source, .. } => Some(MoleculePropertyError::from_rust(source.clone())),
    ck::MoleculeFragmentsFailure::Build(value) => OrderedFragmentBuildError::from_rust(value.clone()).cause(py),
    ck::MoleculeFragmentsFailure::FinalSanitize { source, .. } => Some(SanitizeError::from_rust(source.clone())),
}), source = Option<FailureSource>, methods = {
    #[getter]
    fn component_index(&self) -> Option<usize> { match self.inner {
        ck::MoleculeFragmentsFailure::FinalComputedProperties { component_index, .. }
        | ck::MoleculeFragmentsFailure::FinalSanitize { component_index, .. } => Some(component_index), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn build(&self, py: Python<'_>) -> PyResult<Option<Py<OrderedFragmentBuildError>>> { match &self.inner {
        ck::MoleculeFragmentsFailure::Build(value) => OrderedFragmentBuildError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[pyclass(module="fixed_projection._native", extends=PyValueError)]
#[doc = r#"{"role":"error_record"}"#]
pub struct MoleculeFragmentsError {
    inner: ck::MoleculeFragmentsError,
}
impl MoleculeFragmentsError {
    pub(crate) fn from_rust(inner: ck::MoleculeFragmentsError) -> PyErr {
        Python::attach(|py| {
            let cause = MoleculeFragmentsFailure::from_rust(inner.failure().clone());
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, Some(cause));
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl MoleculeFragmentsError {
    #[getter]
    fn domain(&self) -> &'static str {
        "fragments"
    }
    #[getter]
    fn component_index(&self) -> Option<usize> {
        self.inner.component_index()
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Py<MoleculeFragmentsFailure>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn failure(slf: Bound<'_, Self>) -> PyResult<Py<MoleculeFragmentsFailure>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragment_failures_preserve_native_copy_layers_fields_and_actual_facade_failure() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let inner = ck::OrderedFragmentBuildError::SlowFullCopy {
                component_index: usize::MAX,
                source: ck::FullCopyComponentError::SourceUInt(
                    ck::PropertyUIntReadError::Lexical {
                        value: "wrong".into(),
                        source: ck::UIntLexicalReadError::Character {
                            position: 3,
                            byte: b'x',
                        },
                    },
                ),
            };
            let native = ck::MoleculeFragmentsError::from(inner);
            assert!(
                native
                    .source()
                    .unwrap()
                    .is::<ck::MoleculeFragmentsFailure>()
            );
            let message = native.to_string();
            let error = MoleculeFragmentsError::from_rust(native);
            let value = error.value(py);
            assert_eq!(value.str()?.to_str()?, message);
            assert_eq!(
                value.getattr("component_index")?.extract::<usize>()?,
                usize::MAX
            );
            let failure = value.getattr("source")?;
            assert!(failure.is(&value.getattr("__cause__")?));
            assert!(failure.is(&value.getattr("failure")?));
            assert!(failure.is_instance_of::<MoleculeFragmentsFailure>());
            let full = failure.getattr("source")?;
            assert!(full.is_instance_of::<FullCopyComponentError>());
            assert!(
                failure
                    .getattr("build")?
                    .is_instance_of::<OrderedFragmentBuildError>()
            );
            assert!(
                full.getattr("source_uint")?
                    .is_instance_of::<PropertyUIntReadError>()
            );
            assert!(
                full.getattr("source")?
                    .is_instance_of::<UIntLexicalReadError>()
            );
            assert_eq!(
                full.getattr("source")?
                    .getattr("position")?
                    .extract::<usize>()?,
                3
            );
            assert!(value.setattr("component_index", 0).is_err());

            let native = ck::AtomPathSubsetCopyError::SubstanceGroupCopy(
                ck::SelectedSubstanceGroupCopyError::AtomOutOfRange {
                    group: ck::SubstanceGroupId::new(7),
                    list: ck::SubstanceGroupAtomList::ParentAtoms,
                    source_atom: ck::AtomId::new(4),
                    mapped_atom: ck::AtomId::new(8),
                    atom_count: 6,
                },
            );
            assert!(native.source().is_none());
            let error = AtomPathSubsetCopyError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let payload = error.value(py).getattr("substance_group_copy")?;
            assert_eq!(payload.getattr("group")?.extract::<usize>()?, 7);
            assert_eq!(payload.getattr("source_atom")?.extract::<usize>()?, 4);
            assert_eq!(payload.getattr("mapped_atom")?.extract::<usize>()?, 8);
            assert!(
                payload
                    .getattr("list")?
                    .extract::<SubstanceGroupAtomList>()?
                    == SubstanceGroupAtomList::ParentAtoms
            );
            assert!(
                payload
                    .setattr("list", SubstanceGroupAtomList::Atoms)
                    .is_err()
            );

            let molecule = ck::Molecule::from_smiles_with_params(
                "c1cccc1",
                &ck::SmilesParseParams {
                    sanitize: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let before = molecule.to_binary().unwrap();
            let native = molecule.fragments().unwrap_err();
            let ck::OperationError::Fragments(native) = native else {
                panic!("expected actual fragments sanitation failure")
            };
            let expected = native.to_string();
            let error = MoleculeFragmentsError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, expected);
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .getattr("source")?
                    .is_instance_of::<SanitizeError>()
            );
            assert_eq!(molecule.to_binary().unwrap(), before);
            Ok(())
        })
        .unwrap();
    }
}
