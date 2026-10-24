//! Exact product/run/apply failures; transparent wrappers retain typed payloads.
use crate::canonical_errors::PropertyStringError;
use crate::double_stereo_errors::{DoubleBondStereoError, DoubleStereoSource};
use crate::error_framework::{Source, error_projection};
use crate::match_errors::SubstructMatchError;
use crate::model_edit_errors::{
    AdjacencyError, AtomPropertyError, BondValueError, MoleculePropertyError,
    QueryAtomConversionError,
};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::property_numeric_errors::{
    PropertyIntReadError, PropertyUIntReadError, UIntLexicalReadError,
};
use crate::reaction_validation::{ReactionInitializationError, ReactionValidationError};
use crate::ring_errors::{RingFindingError, RingSource};
use crate::smiles_parse_errors::{SmilesParseError, SmilesSource};
use crate::stereo_order_errors::StereoOrderError;
use crate::topology_edit_errors::{EditSource, TopologyEditError};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type ProductSource = Source<
    Source<EditSource, RingSource>,
    Source<
        Source<DoubleStereoSource, SmilesSource>,
        Source<
            Source<
                crate::error_framework::ErrorObject<StereoGroupError>,
                crate::error_framework::ErrorObject<UIntLexicalReadError>,
            >,
            Source<
                crate::error_framework::ErrorObject<TopologyValidationError>,
                Source<
                    crate::error_framework::ErrorObject<SubstructMatchError>,
                    Source<
                        crate::error_framework::ErrorObject<PropertyIntReadError>,
                        Source<
                            crate::error_framework::ErrorObject<PropertyUIntReadError>,
                            crate::error_framework::ErrorObject<AtomPropertyError>,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(ReactionProductError, ReactionProductErrorKind, ck::ReactionProductError, "reaction", {
    StereoGroup(_), Coordinate(_), TopologyEdit(_), Adjacency(_), RingFinding(_), SourceUInt(_),
    Valence(_), StereoOrder(_), DoubleBondStereo(_), StereoGetter { atom, source },
    CoordinateSelection(_), CarrierIdentity(_), PropertyInt { atom, key, source },
    PropertyUInt { atom, key, source }, MissingProperty { atom, key }, MoleculeProperty(_),
    AtomProperty(_), BondValue(_), TemplateProperty(_),
    Invariant { stage: &'static str, detail: &'static str, reactant_atom: Option<usize>, product_atom: Option<usize>, bond },
    RowOverflow { kind: &'static str, index: usize }, QueryReactantBond { bond },
}, cause = |inner: &ck::ReactionProductError| Python::attach(|py| match inner {
    ck::ReactionProductError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::ReactionProductError::Coordinate(value) => CoordinateValidationError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::TopologyEdit(value) => TopologyEditError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::Adjacency(value) => AdjacencyError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::SourceUInt(value) => PropertyUIntReadError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::DoubleBondStereo(value) => DoubleBondStereoError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::StereoGetter { source, .. } => Some(SubstructMatchError::from_rust(source.clone())),
    ck::ReactionProductError::CoordinateSelection(value) => SmilesParseError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::CarrierIdentity(value) => QueryAtomConversionError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::PropertyInt { source, .. } => Some(PropertyIntReadError::from_rust(source.clone())),
    ck::ReactionProductError::PropertyUInt { source, .. } => Some(PropertyUIntReadError::from_rust(source.clone())),
    ck::ReactionProductError::MoleculeProperty(value) => MoleculePropertyError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::BondValue(value) => BondValueError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::TemplateProperty(value) => ReactionValidationError::from_rust(value.clone()).cause(py),
    ck::ReactionProductError::MissingProperty { .. } | ck::ReactionProductError::Invariant { .. }
    | ck::ReactionProductError::RowOverflow { .. } | ck::ReactionProductError::QueryReactantBond { .. } => None,
}), source = Option<ProductSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::ReactionProductError::StereoGetter { atom, .. } | ck::ReactionProductError::PropertyInt { atom, .. }
        | ck::ReactionProductError::PropertyUInt { atom, .. } | ck::ReactionProductError::MissingProperty { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn key(&self) -> Option<&'static str> { match self.inner {
        ck::ReactionProductError::PropertyInt { key, .. } | ck::ReactionProductError::PropertyUInt { key, .. }
        | ck::ReactionProductError::MissingProperty { key, .. } => Some(key), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::ReactionProductError::Invariant { bond, .. } => bond.map(ck::BondId::index),
        ck::ReactionProductError::QueryReactantBond { bond } => Some(bond.index()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinate(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>> { match &self.inner {
        ck::ReactionProductError::Coordinate(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology_edit(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyEditError>>> { match &self.inner {
        ck::ReactionProductError::TopologyEdit(value) => TopologyEditError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn adjacency(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AdjacencyError>>> { match &self.inner {
        ck::ReactionProductError::Adjacency(value) => AdjacencyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RingFindingError>>> { match &self.inner {
        ck::ReactionProductError::RingFinding(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn source_uint(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyUIntReadError>>> { match &self.inner {
        ck::ReactionProductError::SourceUInt(value) => PropertyUIntReadError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::ReactionProductError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn stereo_order(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<StereoOrderError>>> { match &self.inner {
        ck::ReactionProductError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn double_bond_stereo(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<DoubleBondStereoError>>> { match &self.inner {
        ck::ReactionProductError::DoubleBondStereo(value) => DoubleBondStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinate_selection(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<SmilesParseError>>> { match &self.inner {
        ck::ReactionProductError::CoordinateSelection(value) => SmilesParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn carrier_identity(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<QueryAtomConversionError>>> { match &self.inner {
        ck::ReactionProductError::CarrierIdentity(value) => QueryAtomConversionError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn molecule_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MoleculePropertyError>>> { match &self.inner {
        ck::ReactionProductError::MoleculeProperty(value) => MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atom_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtomPropertyError>>> { match &self.inner {
        ck::ReactionProductError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn bond_value(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<BondValueError>>> { match &self.inner {
        ck::ReactionProductError::BondValue(value) => BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn template_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ReactionValidationError>>> { match &self.inner {
        ck::ReactionProductError::TemplateProperty(value) => ReactionValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

pub(crate) type RunSource = Source<
    crate::error_framework::ErrorObject<ReactionProductError>,
    Source<
        crate::error_framework::ErrorObject<ReactionValidationError>,
        crate::error_framework::ErrorObject<SubstructMatchError>,
    >,
>;

error_projection!(ReactionRunError, ReactionRunErrorKind, ck::ReactionRunError, "reaction", {
    Product { set: usize, template, source }, CoordinateSelectionArity { expected, actual },
    Initialization(_), NeedsInitialization, ReactantArity { expected, actual },
    ReactantTemplateIndex { index: usize, count }, Matching { reactant, template, source },
    MatchAtomIndex { reactant, template, atom: usize, atom_count: usize },
    CombinationLevel { level: u32, count }, CombinationSize { expected, actual }, EmptyCombinationLevels,
}, cause = |inner: &ck::ReactionRunError| Python::attach(|py| match inner {
    ck::ReactionRunError::Product { source, .. } => Some(ReactionProductError::from_rust(source.clone())),
    ck::ReactionRunError::Initialization(value) => ReactionInitializationError::from_rust(value.clone()).cause(py),
    ck::ReactionRunError::Matching { source, .. } => Some(SubstructMatchError::from_rust(source.clone())),
    ck::ReactionRunError::CoordinateSelectionArity { .. } | ck::ReactionRunError::NeedsInitialization
    | ck::ReactionRunError::ReactantArity { .. } | ck::ReactionRunError::ReactantTemplateIndex { .. }
    | ck::ReactionRunError::MatchAtomIndex { .. } | ck::ReactionRunError::CombinationLevel { .. }
    | ck::ReactionRunError::CombinationSize { .. } | ck::ReactionRunError::EmptyCombinationLevels => None,
}), source = Option<RunSource>, methods = {
    #[getter]
    fn template(&self) -> Option<usize> { match self.inner {
        ck::ReactionRunError::Product { template, .. } | ck::ReactionRunError::Matching { template, .. }
        | ck::ReactionRunError::MatchAtomIndex { template, .. } => Some(template), _ => None,
    } }
    #[getter]
    fn reactant(&self) -> Option<usize> { match self.inner {
        ck::ReactionRunError::Matching { reactant, .. } | ck::ReactionRunError::MatchAtomIndex { reactant, .. } => Some(reactant), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::ReactionRunError::CoordinateSelectionArity { expected, .. } | ck::ReactionRunError::ReactantArity { expected, .. } => Some(expected),
        ck::ReactionRunError::CombinationSize { expected, .. } => Some(expected as usize), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::ReactionRunError::CoordinateSelectionArity { actual, .. } | ck::ReactionRunError::ReactantArity { actual, .. }
        | ck::ReactionRunError::CombinationSize { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn count(&self) -> Option<usize> { match self.inner {
        ck::ReactionRunError::ReactantTemplateIndex { count, .. } => Some(count),
        ck::ReactionRunError::CombinationLevel { count, .. } => Some(count as usize), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn initialization(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ReactionInitializationError>>> { match &self.inner {
        ck::ReactionRunError::Initialization(value) => ReactionInitializationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_reaction_failures_preserve_transparent_causes_and_wrapped_payloads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CC").unwrap();
            let before = molecule.to_binary().unwrap();
            let mut reaction = ck::parse_smirks("[C:1]>>[N:1]").unwrap();
            let failure = molecule
                .reaction_products(&mut reaction, usize::MAX)
                .unwrap_err();
            let ck::OperationError::ReactionRun(native) = failure else {
                panic!("actual run error expected")
            };
            assert!(matches!(
                native,
                ck::ReactionRunError::ReactantTemplateIndex {
                    index: usize::MAX,
                    count: 1
                }
            ));
            let display = native.to_string();
            let error = ReactionRunError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("index")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(error.value(py).getattr("count")?.extract::<usize>()?, 1);
            assert!(error.value(py).getattr("source")?.is_none());
            assert_eq!(molecule.to_binary().unwrap(), before);

            let native = ck::ReactionApplyError::Matching(ck::ReactionRunError::Product {
                set: 4,
                template: 2,
                source: ck::ReactionProductError::PropertyInt {
                    atom: ck::AtomId::new(usize::MAX),
                    key: "tracking",
                    source: ck::PropertyIntReadError::Lexical {
                        value: "bad".into(),
                        source: ck::UIntLexicalReadError::Character {
                            position: 0,
                            byte: b'b',
                        },
                    },
                },
            });
            // Apply's transparent Matching skips the Run wrapper as a cause,
            // but the readonly matching payload still owns that full value.
            assert!(native.source().unwrap().is::<ck::ReactionProductError>());
            let display = native.to_string();
            let error = ReactionApplyError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let product = error.value(py).getattr("source")?;
            assert!(product.is(&error.value(py).getattr("__cause__")?));
            assert!(product.is_instance(<ReactionProductError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(product.getattr("atom")?.extract::<usize>()?, usize::MAX);
            assert_eq!(product.getattr("key")?.extract::<String>()?, "tracking");
            let number = product.getattr("source")?;
            assert!(number.is_instance(<PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(number.getattr("value")?.extract::<String>()?, "bad");
            let lexical = number.getattr("source")?;
            assert!(lexical.is_instance(<UIntLexicalReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(lexical.getattr("byte")?.extract::<u8>()?, b'b');
            let matching = error.value(py).getattr("matching")?;
            assert!(matching.is_instance(<ReactionRunError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(matching.getattr("set")?.extract::<usize>()?, 4);
            assert_eq!(matching.getattr("template")?.extract::<usize>()?, 2);
            assert!(error.value(py).setattr("matching", py.None()).is_err());

            let report = ck::Reaction::new().validate().unwrap();
            assert!(!report.is_valid());
            let expected_errors = report.num_errors();
            let native =
                ck::ReactionRunError::Initialization(ck::ReactionInitializationError::Invalid {
                    report,
                });
            assert!(native.source().is_none());
            let error = ReactionRunError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let initialization = error.value(py).getattr("initialization")?;
            assert!(initialization.is_instance(<ReactionInitializationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(
                initialization
                    .getattr("report")?
                    .getattr("num_errors")?
                    .extract::<usize>()?,
                expected_errors
            );

            let native =
                ck::ReactionProductError::TemplateProperty(ck::ReactionValidationError::Property {
                    role: ck::ReactionRole::Reactant,
                    template: 3,
                    atom: ck::AtomId::new(2),
                    property: "molAtomMapNumber",
                    source: ck::PropertyIntReadError::InvalidKind {
                        kind: ck::PropertyValueKind::Bool,
                    },
                });
            assert!(native.source().unwrap().is::<ck::PropertyIntReadError>());
            let error = ReactionProductError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance(<PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            let payload = error.value(py).getattr("template_property")?;
            assert!(payload.is_instance(<ReactionValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(payload.getattr("template")?.extract::<usize>()?, 3);

            let native =
                ck::ReactionApplyError::SourceText(ck::PropertyStringError::UnsupportedKind {
                    kind: ck::PropertyValueKind::Bool,
                });
            assert!(native.source().is_none());
            let error = ReactionApplyError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let text = error.value(py).getattr("source_text")?;
            assert!(text.is_instance(<PropertyStringError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                text.getattr("kind")?
                    .extract::<crate::property_values::PropertyValueKind>()?
                    == crate::property_values::PropertyValueKind::Bool
            );
            let error = ReactionProductError::from_rust(ck::ReactionProductError::RowOverflow {
                kind: "bond",
                index: usize::MAX,
            });
            assert_eq!(
                error.value(py).getattr("kind")?.extract::<String>()?,
                "bond"
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<ReactionProductErrorKind>()?
                    == ReactionProductErrorKind::RowOverflow
            );
            assert_eq!(
                error.value(py).getattr("index")?.extract::<usize>()?,
                usize::MAX
            );
            Ok(())
        })
        .unwrap();
    }
}

type ApplySource = Source<RunSource, Source<ProductSource, EditSource>>;
error_projection!(ReactionApplyError, ReactionApplyErrorKind, ck::ReactionApplyError, "reaction", {
    SourceText(_), Initialization(_), ApplicabilityArity { reactants: usize, products: usize },
    AddsProductAtom { atom: usize => |value: &ck::AtomId| value.index() }, Matching(_), Product(_), Edit(_),
}, cause = |inner: &ck::ReactionApplyError| Python::attach(|py| match inner {
    ck::ReactionApplyError::SourceText(value) => PropertyStringError::from_rust(*value).cause(py),
    ck::ReactionApplyError::Initialization(value) => ReactionInitializationError::from_rust(value.clone()).cause(py),
    ck::ReactionApplyError::Matching(value) => ReactionRunError::from_rust(value.clone()).cause(py),
    ck::ReactionApplyError::Product(value) => ReactionProductError::from_rust(value.clone()).cause(py),
    ck::ReactionApplyError::Edit(value) => TopologyEditError::from_rust(value.clone()).cause(py),
    ck::ReactionApplyError::ApplicabilityArity { .. } | ck::ReactionApplyError::AddsProductAtom { .. } => None,
}), source = Option<ApplySource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn source_text(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyStringError>>> { match &self.inner {
        ck::ReactionApplyError::SourceText(value) => PropertyStringError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn initialization(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ReactionInitializationError>>> { match &self.inner {
        ck::ReactionApplyError::Initialization(value) => ReactionInitializationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn matching(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ReactionRunError>>> { match &self.inner {
        ck::ReactionApplyError::Matching(value) => ReactionRunError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn product(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ReactionProductError>>> { match &self.inner {
        ck::ReactionApplyError::Product(value) => ReactionProductError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn edit(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyEditError>>> { match &self.inner {
        ck::ReactionApplyError::Edit(value) => TopologyEditError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});
