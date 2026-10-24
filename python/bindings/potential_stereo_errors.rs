//! Potential-stereo error values preserve native fields and exact cause edges.
use crate::canonical_errors::CanonicalRankError;
use crate::chemical_values::BondOrder;
use crate::double_stereo_errors::DoubleBondStereoError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::BondValueError;
use crate::model_errors::TopologyValidationError;
use crate::property_values::PropertyValueKind;
use crate::ring_errors::RingFindingError;
use crate::stereo_order_errors::StereoOrderError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use pyo3::types::PyInt;
use std::error::Error;

pub(crate) type PotentialSource = Source<
    crate::error_framework::ErrorObject<RingFindingError>,
    Source<
        crate::error_framework::ErrorObject<CanonicalRankError>,
        Source<
            crate::error_framework::ErrorObject<TopologyValidationError>,
            Source<
                crate::error_framework::ErrorObject<ValenceError>,
                Source<
                    crate::error_framework::ErrorObject<StereoOrderError>,
                    Source<
                        crate::error_framework::ErrorObject<DoubleBondStereoError>,
                        crate::error_framework::ErrorObject<BondValueError>,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(PotentialStereoError, PotentialStereoErrorKind, ck::PotentialStereoError, "stereo", {
    InvalidPropertyKind { atom, property: &'static str, kind },
    InvalidRingStereoReference { atom, value }, EmptyRingStereoReferences { atom },
    Ring(_), Rank(_), InvalidTopology(_),
    InvalidValence { field, actual: usize, atom_count }, InvalidValenceValue { field, atom, value },
    AtomOutOfRange { atom, atom_count },
    InvalidRingInfo { reason, row: usize, value, limit: usize }, InvalidAtomDegree { atom, degree },
    InvalidBondDegree { bond, endpoint: &'static str, degree }, InvalidStereoReferences { bond, reason },
    UnsupportedBondOrder { bond, order: BondOrder => |value: &ck::BondOrder| (*value).into() },
    InvalidChiralPermutation { atom, permutation: u32 }, AtropisomerDependencyUnavailable { bond },
    RefinementDidNotConverge { iterations: usize }, Valence(_), StereoOrder(_), DoubleStereo(_), BondValue(_),
}, cause = |inner: &ck::PotentialStereoError| match inner {
    ck::PotentialStereoError::Ring(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::PotentialStereoError::Rank(value) => Some(CanonicalRankError::from_rust(value.clone())),
    ck::PotentialStereoError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::PotentialStereoError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::PotentialStereoError::StereoOrder(value) => Some(StereoOrderError::from_rust(value.clone())),
    ck::PotentialStereoError::DoubleStereo(value) => Some(DoubleBondStereoError::from_rust(value.clone())),
    ck::PotentialStereoError::BondValue(value) => Some(BondValueError::from_rust(value.clone())),
    ck::PotentialStereoError::InvalidPropertyKind { .. } | ck::PotentialStereoError::InvalidRingStereoReference { .. }
    | ck::PotentialStereoError::EmptyRingStereoReferences { .. } | ck::PotentialStereoError::InvalidValence { .. }
    | ck::PotentialStereoError::InvalidValenceValue { .. } | ck::PotentialStereoError::AtomOutOfRange { .. }
    | ck::PotentialStereoError::InvalidRingInfo { .. } | ck::PotentialStereoError::InvalidAtomDegree { .. }
    | ck::PotentialStereoError::InvalidBondDegree { .. } | ck::PotentialStereoError::InvalidStereoReferences { .. }
    | ck::PotentialStereoError::UnsupportedBondOrder { .. } | ck::PotentialStereoError::InvalidChiralPermutation { .. }
    | ck::PotentialStereoError::AtropisomerDependencyUnavailable { .. }
    | ck::PotentialStereoError::RefinementDidNotConverge { .. } => None,
}, source = Option<PotentialSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::PotentialStereoError::InvalidPropertyKind { atom, .. }
        | ck::PotentialStereoError::InvalidRingStereoReference { atom, .. }
        | ck::PotentialStereoError::EmptyRingStereoReferences { atom }
        | ck::PotentialStereoError::InvalidValenceValue { atom, .. }
        | ck::PotentialStereoError::AtomOutOfRange { atom, .. }
        | ck::PotentialStereoError::InvalidAtomDegree { atom, .. }
        | ck::PotentialStereoError::InvalidChiralPermutation { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::PotentialStereoError::InvalidBondDegree { bond, .. }
        | ck::PotentialStereoError::InvalidStereoReferences { bond, .. }
        | ck::PotentialStereoError::UnsupportedBondOrder { bond, .. }
        | ck::PotentialStereoError::AtropisomerDependencyUnavailable { bond } => Some(bond.index()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::PotentialStereoError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    fn field(&self) -> Option<&'static str> { match self.inner {
        ck::PotentialStereoError::InvalidValence { field, .. }
        | ck::PotentialStereoError::InvalidValenceValue { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::PotentialStereoError::InvalidValence { atom_count, .. }
        | ck::PotentialStereoError::AtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn value<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyInt>>> { match self.inner {
        ck::PotentialStereoError::InvalidRingStereoReference { value, .. }
        | ck::PotentialStereoError::InvalidValenceValue { value, .. } => value.into_pyobject(py).map(Some).map_err(Into::into),
        ck::PotentialStereoError::InvalidRingInfo { value, .. } => value.into_pyobject(py).map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    fn reason(&self) -> Option<&'static str> { match self.inner {
        ck::PotentialStereoError::InvalidRingInfo { reason, .. }
        | ck::PotentialStereoError::InvalidStereoReferences { reason, .. } => Some(reason), _ => None,
    } }
    #[getter]
    fn degree(&self) -> Option<usize> { match self.inner {
        ck::PotentialStereoError::InvalidAtomDegree { degree, .. }
        | ck::PotentialStereoError::InvalidBondDegree { degree, .. } => Some(degree), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::double_stereo_errors::DoubleBondStereoErrorKind;
    use crate::property_numeric_errors::PropertyIntReadError;

    #[test]
    fn stereo_error_chain_preserves_native_context_and_integer_width() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::PotentialStereoError::DoubleStereo(
                ck::DoubleBondStereoError::NumericPropertyRead {
                    atom: Some(ck::AtomId::new(2)),
                    bond: None,
                    property: "_UnknownStereo",
                    source: ck::PropertyIntReadError::UnsignedOverflow { value: u32::MAX },
                },
            );
            assert!(native.source().is_some());
            let display = native.to_string();
            let projected = PotentialStereoError::from_rust(native);
            let value = projected.value(py);
            assert_eq!(value.str()?.to_str()?, display);
            let double = value.getattr("source")?;
            assert!(double.is(&value.getattr("__cause__")?));
            assert!(double.is_instance(<DoubleBondStereoError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                double
                    .getattr("variant")?
                    .extract::<DoubleBondStereoErrorKind>()?
                    == DoubleBondStereoErrorKind::NumericPropertyRead
            );
            assert_eq!(double.getattr("atom")?.extract::<usize>()?, 2);
            assert!(double.getattr("bond")?.is_none());
            assert_eq!(
                double.getattr("property")?.extract::<String>()?,
                "_UnknownStereo"
            );
            assert!(double.setattr("property", "other").is_err());
            let number = double.getattr("source")?;
            assert!(number.is(&double.getattr("__cause__")?));
            assert!(number.is_instance(<PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(number.getattr("value")?.extract::<u32>()?, u32::MAX);
            assert!(number.getattr("__cause__")?.is_none());
            for native in [
                ck::PotentialStereoError::InvalidRingStereoReference {
                    atom: ck::AtomId::new(0),
                    value: -7,
                },
                ck::PotentialStereoError::InvalidRingInfo {
                    reason: "index",
                    row: 1,
                    value: usize::MAX,
                    limit: 5,
                },
            ] {
                let signed = matches!(
                    native,
                    ck::PotentialStereoError::InvalidRingStereoReference { .. }
                );
                assert!(native.source().is_none());
                let error = PotentialStereoError::from_rust(native);
                if signed {
                    assert_eq!(error.value(py).getattr("value")?.extract::<i32>()?, -7);
                } else {
                    assert_eq!(
                        error.value(py).getattr("value")?.extract::<usize>()?,
                        usize::MAX
                    );
                }
                assert!(error.cause(py).is_none());
            }
            Ok(())
        })
        .unwrap();
    }
}
