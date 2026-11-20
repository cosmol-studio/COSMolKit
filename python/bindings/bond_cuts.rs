//! Python projection of a cut instruction and its structured errors.
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{ErrorObject, Source, error_projection};
use crate::model_edit_errors::{AtomPropertyError, BondValueError};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::property_numeric_errors::{PropertyUIntReadError, UIntLexicalReadError};
use crate::property_values::PropertyValueKind;
use crate::stereo_order_errors::StereoOrderError;
use crate::topology_edit_errors::{EditSource, TopologyEditError};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

macro_rules! project_brics_match {
    ($name:ident, $doc:literal) => {
        crate::configuration::type_projection! {
            value(eq)
            #[doc = $doc]
            pub struct $name { pub(crate) inner: ck::$name }
        }

        #[pymethods]
        impl $name {
            /// Zero-based index of the bond in the molecule used for discovery.
            #[getter]
            fn bond(&self) -> usize {
                self.inner.bond().index()
            }

            /// Two (source atom index, BRICS environment label) pairs in rule order.
            /// The first atom need not be the stored bond's begin atom.
            #[getter]
            fn ends(&self) -> ((usize, u8), (usize, u8)) {
                let [left, right] = self.inner.ends().map(|(a, e)| (a.index(), e));
                (left, right)
            }

            /// Select this source bond with each dummy label on its matched atom's side.
            fn to_bond_cut(&self) -> BondCut {
                BondCut {
                    inner: self.inner.to_bond_cut(),
                }
            }
        }
    };
}
project_brics_match!(
    BricsBondMatch,
    "A site with oriented standard BRICS endpoint labels."
);
project_brics_match!(
    RBricsBondMatch,
    "A site with oriented regular r-BRICS endpoint labels."
);

crate::configuration::type_projection! {
    value(eq)
    /// A bond selected in the input molecule, with optional attachment-side isotope labels.
    pub struct BondCut { pub(crate) inner: ck::BondCut }
}

#[pymethods]
impl BondCut {
    #[staticmethod]
    #[pyo3(signature = (bond, begin_label=None, end_label=None))]
    /// Select a bond by its zero-based index. Supply both labels or omit both.
    /// An explicit label stays on its named endpoint's side. By default each dummy
    /// takes the index of the atom it replaces. Both cut methods check that the
    /// effective labels fit in the isotope field (0 through 65535).
    fn at_bond(
        bond: usize,
        begin_label: Option<usize>,
        end_label: Option<usize>,
    ) -> PyResult<Self> {
        ck::BondCut::at_bond(ck::BondId::new(bond), begin_label, end_label)
            .map(|inner| Self { inner })
            .map_err(BondCutError::from_rust)
    }
}

type BricsSource = Source<
    ErrorObject<crate::reaction_execution_errors::ReactionRunError>,
    Source<
        ErrorObject<crate::sanitize_errors::SanitizeError>,
        Source<
            ErrorObject<crate::smiles_parse_errors::SmilesParseError>,
            ErrorObject<crate::utf8_errors::FromUtf8Error>,
        >,
    >,
>;

error_projection!(BricsDecompositionError, BricsDecompositionErrorKind, ck::BricsDecompositionError, "fragmentation", {
    Reaction(_), UnkeyedProduct(_), Smiles(_), Encoding(_),
}, cause = |error: &ck::BricsDecompositionError| Some(match error {
    ck::BricsDecompositionError::Reaction(e) => crate::reaction_execution_errors::ReactionRunError::from_rust(e.clone()),
    ck::BricsDecompositionError::UnkeyedProduct(e) => crate::sanitize_errors::SanitizeError::from_rust(e.clone()),
    ck::BricsDecompositionError::Smiles(e) => crate::smiles_parse_errors::SmilesParseError::from_rust(e.clone()),
    ck::BricsDecompositionError::Encoding(e) => crate::utf8_errors::FromUtf8Error::from_rust(e.clone()),
}), source = BricsSource);

type CutSource = Source<
    EditSource,
    Source<
        ErrorObject<PropertyUIntReadError>,
        Source<
            ErrorObject<AtomPropertyError>,
            Source<ErrorObject<UIntLexicalReadError>, ErrorObject<TopologyValidationError>>,
        >,
    >,
>;

error_projection!(BondCutError, BondCutErrorKind, ck::BondCutError, "fragmentation", {
    IncompleteLabels { bond }, EmptySelection, Duplicate { bond }, BondOutOfRange { bond, count: usize },
    EndpointsChanged { bond, expected: (usize,usize) => |v: &[ck::AtomId;2]| (v[0].index(),v[1].index()), actual: (usize,usize) => |v: &[ck::AtomId;2]| (v[0].index(),v[1].index()) },
    IsotopeOutOfRange { bond, attachment: usize => |v: &ck::AtomId| v.index(), label: usize },
    IncidentOrderProperty { atom: usize => |v: &ck::AtomId| v.index(), kind: PropertyValueKind => |v: &ck::PropertyValueKind| (*v).into() },
    Edit(_), Coordinates(_), Stereo(_), BondProperty(_), AtomProperty(_), NumericProperty(_), TextProperty(_), Valence(_),
}, cause = |inner: &ck::BondCutError| Python::attach(|py| match inner {
    ck::BondCutError::Edit(e) => TopologyEditError::from_rust(e.clone()).cause(py),
    ck::BondCutError::Coordinates(e) => CoordinateValidationError::from_rust(e.clone()).cause(py),
    ck::BondCutError::Stereo(e) => StereoOrderError::from_rust(e.clone()).cause(py),
    ck::BondCutError::BondProperty(e) => BondValueError::from_rust(e.clone()).cause(py),
    ck::BondCutError::AtomProperty(e) => AtomPropertyError::from_rust(e.clone()).cause(py),
    ck::BondCutError::NumericProperty(e) => PropertyUIntReadError::from_rust(e.clone()).cause(py),
    ck::BondCutError::TextProperty(e) => PropertyStringError::from_rust(e.clone()).cause(py),
    ck::BondCutError::Valence(e) => ValenceError::from_rust(e.clone()).cause(py),
    ck::BondCutError::IncompleteLabels{..} | ck::BondCutError::EmptySelection | ck::BondCutError::Duplicate{..}
    | ck::BondCutError::BondOutOfRange{..} | ck::BondCutError::EndpointsChanged{..} | ck::BondCutError::IsotopeOutOfRange{..} | ck::BondCutError::IncidentOrderProperty{..} => None,
}), source = Option<CutSource>, methods = {
    #[getter]
    /// Source bond selected by the failed instruction, when available.
    fn bond(&self) -> Option<usize> { match &self.inner {
        ck::BondCutError::IncompleteLabels{bond} | ck::BondCutError::Duplicate{bond} | ck::BondCutError::BondOutOfRange{bond,..}
        | ck::BondCutError::EndpointsChanged{bond,..} | ck::BondCutError::IsotopeOutOfRange{bond,..} => Some(bond.index()), _=>None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native TopologyEditError for the Edit variant; otherwise None."]
    fn edit(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<TopologyEditError>>> {
        match &self.inner { ck::BondCutError::Edit(e) => TopologyEditError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native CoordinateValidationError for the Coordinates variant; otherwise None."]
    fn coordinates(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<CoordinateValidationError>>> {
        match &self.inner { ck::BondCutError::Coordinates(e) => CoordinateValidationError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native StereoOrderError for the Stereo variant; otherwise None."]
    fn stereo(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<StereoOrderError>>> {
        match &self.inner { ck::BondCutError::Stereo(e) => StereoOrderError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native BondValueError for the BondProperty variant; otherwise None."]
    fn bond_property(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<BondValueError>>> {
        match &self.inner { ck::BondCutError::BondProperty(e) => BondValueError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native AtomPropertyError for the AtomProperty variant; otherwise None."]
    fn atom_property(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<AtomPropertyError>>> {
        match &self.inner { ck::BondCutError::AtomProperty(e) => AtomPropertyError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native PropertyUIntReadError for the NumericProperty variant; otherwise None."]
    fn numeric_property(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<PropertyUIntReadError>>> {
        match &self.inner { ck::BondCutError::NumericProperty(e) => PropertyUIntReadError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native PropertyStringError for the TextProperty variant; otherwise None."]
    fn text_property(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<PropertyStringError>>> {
        match &self.inner { ck::BondCutError::TextProperty(e) => PropertyStringError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native ValenceError for the Valence variant; otherwise None."]
    fn valence(&self, py:Python<'_>) -> PyResult<Option<ErrorObject<ValenceError>>> {
        match &self.inner { ck::BondCutError::Valence(e) => ValenceError::from_rust(e.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None) }
    }
});
