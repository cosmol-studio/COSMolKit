//! Concrete CIP and atom-code failures retain native payloads and cause levels.
use crate::canonical_errors::PropertyStringError;
use crate::chemical_values::{BondOrder, ChiralTag};
use crate::double_stereo_errors::DoubleBondStereoError;
use crate::error_framework::{Source, error_projection};
use crate::kekulize_errors::KekulizeError;
use crate::model_edit_errors::{AtomPropertyError, BondValueError, MoleculePropertyError};
use crate::model_errors::TopologyValidationError;
use crate::path_errors::PeriodicTableError;
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::property_values::PropertyValueKind;
use crate::ring_errors::RingFindingError;
use crate::stereo_errors::{AtropisomerError, LegacySource};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

#[derive(IntoPyObject)]
enum CipKind {
    Property(PropertyValueKind),
    Index(&'static str),
}

pub(crate) type CipSource = Source<
    LegacySource,
    Source<
        crate::error_framework::ErrorObject<PropertyUIntReadError>,
        crate::error_framework::ErrorObject<MoleculePropertyError>,
    >,
>;

error_projection!(CipLabelerError, CipLabelerErrorKind, ck::CipLabelerError, "cip", {
    UnsignedPropertyRead { atom, property, source }, InvalidPropertyKind { atom, property, kind },
    AtomIndexOutOfRange { index, atom_count: usize }, BondIndexOutOfRange { index, bond_count: usize },
    SourceIndexWidthExceeded { kind, index }, BondNotIncident { bond, atom },
    NonIntegerBondOrder { order: BondOrder => |value: &ck::BondOrder| (*value).into() },
    EdgeEndpointMismatch { edge: usize, node: usize }, TooManyNodes { limit: usize },
    MaxIterationsExceeded, ControlCCaught, UnexpectedUpEdgeOrdering, NoSequenceRuleProvided, EmptyPairListReference, DescriptorListLengthMismatch, InvalidStereoDescriptor, SubstituentsShouldBeTopologicallyEquivalent, SomethingUnexpected, Rule4bInstanceNotInRuleSet, Rule5NewInstanceNotInRuleSet, ParityVectorsMustHaveSize4, EmptyConfigurationFoci, BadTetrahedralConfig, TetrahedralConfigurationMustHave4Carriers, DescriptorNotSupportedForAtoms, InvalidAtomDescriptor, CarrierMismatch, BadSp2BondFoci, BadSp2BondConfig, IncorrectNumberOfStereoAtoms, DescriptorNotSupportedForDoubleBonds, InvalidBondDescriptor, BadAtropisomerBondFoci, BadAtropisomerBondConfig, DescriptorNotSupportedForAtropisomerBonds,
    InvalidInternalState { detail: &'static str },
    UnsupportedConfiguration { atom, tag: ChiralTag => |value: &ck::ChiralTag| (*value).into() },
    InvalidTopology(_), AtomProperty { atom, source }, BondValue { bond, source }, MoleculeProperty(_),
    Atropisomer(_), DoubleStereo(_), RingFinding(_), Valence(_), Kekulize(_), PeriodicTable(_),
}, cause = |inner: &ck::CipLabelerError| Python::attach(|py| match inner {
    ck::CipLabelerError::UnsignedPropertyRead { source, .. } => Some(PropertyUIntReadError::from_rust(source.clone())),
    ck::CipLabelerError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::CipLabelerError::AtomProperty { source, .. } => Some(AtomPropertyError::from_rust(source.clone())),
    ck::CipLabelerError::BondValue { source, .. } => Some(BondValueError::from_rust(source.clone())),
    ck::CipLabelerError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::CipLabelerError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).cause(py),
    ck::CipLabelerError::DoubleStereo(value) => DoubleBondStereoError::from_rust(value.clone()).cause(py),
    ck::CipLabelerError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::CipLabelerError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::CipLabelerError::Kekulize(value) => KekulizeError::from_rust(value.clone()).cause(py),
    ck::CipLabelerError::PeriodicTable(value) => PeriodicTableError::from_rust(value.clone()).cause(py),
    ck::CipLabelerError::InvalidPropertyKind { .. }
    | ck::CipLabelerError::AtomIndexOutOfRange { .. }
    | ck::CipLabelerError::BondIndexOutOfRange { .. }
    | ck::CipLabelerError::SourceIndexWidthExceeded { .. }
    | ck::CipLabelerError::BondNotIncident { .. }
    | ck::CipLabelerError::NonIntegerBondOrder { .. }
    | ck::CipLabelerError::EdgeEndpointMismatch { .. }
    | ck::CipLabelerError::TooManyNodes { .. }
    | ck::CipLabelerError::InvalidInternalState { .. }
    | ck::CipLabelerError::UnsupportedConfiguration { .. }
    | ck::CipLabelerError::MaxIterationsExceeded
    | ck::CipLabelerError::ControlCCaught
    | ck::CipLabelerError::UnexpectedUpEdgeOrdering
    | ck::CipLabelerError::NoSequenceRuleProvided
    | ck::CipLabelerError::EmptyPairListReference
    | ck::CipLabelerError::DescriptorListLengthMismatch
    | ck::CipLabelerError::InvalidStereoDescriptor
    | ck::CipLabelerError::SubstituentsShouldBeTopologicallyEquivalent
    | ck::CipLabelerError::SomethingUnexpected
    | ck::CipLabelerError::Rule4bInstanceNotInRuleSet
    | ck::CipLabelerError::Rule5NewInstanceNotInRuleSet
    | ck::CipLabelerError::ParityVectorsMustHaveSize4
    | ck::CipLabelerError::EmptyConfigurationFoci
    | ck::CipLabelerError::BadTetrahedralConfig
    | ck::CipLabelerError::TetrahedralConfigurationMustHave4Carriers
    | ck::CipLabelerError::DescriptorNotSupportedForAtoms
    | ck::CipLabelerError::InvalidAtomDescriptor
    | ck::CipLabelerError::CarrierMismatch
    | ck::CipLabelerError::BadSp2BondFoci
    | ck::CipLabelerError::BadSp2BondConfig
    | ck::CipLabelerError::IncorrectNumberOfStereoAtoms
    | ck::CipLabelerError::DescriptorNotSupportedForDoubleBonds
    | ck::CipLabelerError::InvalidBondDescriptor
    | ck::CipLabelerError::BadAtropisomerBondFoci
    | ck::CipLabelerError::BadAtropisomerBondConfig
    | ck::CipLabelerError::DescriptorNotSupportedForAtropisomerBonds => None,
}), source = Option<CipSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::CipLabelerError::UnsignedPropertyRead { atom, .. }
        | ck::CipLabelerError::InvalidPropertyKind { atom, .. }
        | ck::CipLabelerError::BondNotIncident { atom, .. }
        | ck::CipLabelerError::UnsupportedConfiguration { atom, .. }
        | ck::CipLabelerError::AtomProperty { atom, .. } => Some(atom), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::CipLabelerError::BondNotIncident { bond, .. }
        | ck::CipLabelerError::BondValue { bond, .. } => Some(bond), _ => None,
    } }
    #[getter]
    fn index(&self) -> Option<usize> { match self.inner {
        ck::CipLabelerError::AtomIndexOutOfRange { index, .. }
        | ck::CipLabelerError::BondIndexOutOfRange { index, .. }
        | ck::CipLabelerError::SourceIndexWidthExceeded { index, .. } => Some(index), _ => None,
    } }
    #[getter]
    fn property(&self) -> Option<&'static str> { match self.inner {
        ck::CipLabelerError::UnsignedPropertyRead { property, .. }
        | ck::CipLabelerError::InvalidPropertyKind { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<CipKind> { match self.inner {
        ck::CipLabelerError::InvalidPropertyKind { kind, .. } => Some(CipKind::Property(kind.into())),
        ck::CipLabelerError::SourceIndexWidthExceeded { kind, .. } => Some(CipKind::Index(kind)),
        _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::CipLabelerError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn molecule_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MoleculePropertyError>>> { match &self.inner {
        ck::CipLabelerError::MoleculeProperty(value) => MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atropisomer(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtropisomerError>>> { match &self.inner {
        ck::CipLabelerError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn double_stereo(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<DoubleBondStereoError>>> { match &self.inner {
        ck::CipLabelerError::DoubleStereo(value) => DoubleBondStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RingFindingError>>> { match &self.inner {
        ck::CipLabelerError::RingFinding(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::CipLabelerError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kekulize(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<KekulizeError>>> { match &self.inner {
        ck::CipLabelerError::Kekulize(value) => KekulizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn periodic_table(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PeriodicTableError>>> { match &self.inner {
        ck::CipLabelerError::PeriodicTable(value) => PeriodicTableError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(AtomCodeError, AtomCodeErrorKind, ck::AtomCodeError, "fingerprint", {
    Precondition { what }, Postcondition { what, code: u32 }, Valence(_), Cip(_), PropertyString(_),
}, cause = |inner: &ck::AtomCodeError| match inner {
    ck::AtomCodeError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::AtomCodeError::Cip(value) => Some(CipLabelerError::from_rust(value.clone())),
    ck::AtomCodeError::PropertyString(value) => Some(PropertyStringError::from_rust(*value)),
    ck::AtomCodeError::Precondition { .. } | ck::AtomCodeError::Postcondition { .. } => None,
}, source = Option<Source<crate::error_framework::ErrorObject<ValenceError>, Source<crate::error_framework::ErrorObject<CipLabelerError>, crate::error_framework::ErrorObject<PropertyStringError>>>>, methods = {
    #[getter]
    fn what(&self) -> Option<&'static str> { match self.inner {
        ck::AtomCodeError::Precondition { what } | ck::AtomCodeError::Postcondition { what, .. } => Some(what),
        _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::property_numeric_errors::UIntLexicalReadError;

    #[test]
    fn atom_code_and_cip_keep_exact_cause_depth_and_native_kind_values() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::AtomCodeError::Cip(ck::CipLabelerError::UnsignedPropertyRead {
                atom: 3,
                property: "_CIPRank",
                source: ck::PropertyUIntReadError::Lexical {
                    value: "bad".into(),
                    source: ck::UIntLexicalReadError::Character {
                        position: 1,
                        byte: b'a',
                    },
                },
            });
            assert!(native.source().unwrap().is::<ck::CipLabelerError>());
            let message = native.to_string();
            let error = AtomCodeError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let cip = error.value(py).getattr("source")?;
            assert!(cip.is_instance(<CipLabelerError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(cip.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(cip.getattr("atom")?.extract::<usize>()?, 3);
            assert_eq!(cip.getattr("property")?.extract::<String>()?, "_CIPRank");
            let numeric = cip.getattr("source")?;
            assert!(numeric.is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                numeric
                    .getattr("source")?
                    .is_instance(<UIntLexicalReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert!(cip.setattr("atom", 7).is_err());

            let error = CipLabelerError::from_rust(ck::CipLabelerError::InvalidPropertyKind {
                atom: 2,
                property: "_CIPCode",
                kind: ck::PropertyValueKind::Bool,
            });
            assert!(
                error
                    .value(py)
                    .getattr("kind")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            let error = CipLabelerError::from_rust(ck::CipLabelerError::SourceIndexWidthExceeded {
                kind: "atom",
                index: usize::MAX,
            });
            assert_eq!(
                error.value(py).getattr("kind")?.extract::<String>()?,
                "atom"
            );
            assert_eq!(
                error.value(py).getattr("index")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<CipLabelerErrorKind>()?
                    == CipLabelerErrorKind::SourceIndexWidthExceeded
            );

            let native = ck::CipLabelerError::Valence(ck::ValenceError::InvalidTopology {
                source: ck::TopologyValidationError::AdjacencyMismatch,
            });
            assert!(native.source().unwrap().is::<ck::TopologyValidationError>());
            let error = CipLabelerError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert!(
                error
                    .value(py)
                    .getattr("valence")?
                    .is_instance(<ValenceError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            Ok(())
        })
        .unwrap();
    }
}
