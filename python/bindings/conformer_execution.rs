//! Native bounds/count queries and their complete generation error owners.
use crate::alignment_errors::ThreadCountError;
use crate::chemical_values::enum_value;
use crate::embedding_values::ConformerError;
use crate::error_framework::{Source, error_projection};
use crate::fragment_errors::MoleculeFragmentsFailure;
use crate::io_source_errors::{ErrorRef, IoError};
use crate::kekulize_errors::KekulizeSource;
use crate::match_errors::MatchSource;
use crate::matrix_errors::MatrixError;
use crate::retained_error::RetainedError;
use crate::smarts_errors::SmartsSource;
use crate::smiles_parse_errors::SmilesSource;
use crate::stereo_errors::{AtropSource, LegacySource};
use crate::uff_parameter_errors::{
    ForceFieldKernelError, UffBondMathError, UffBuilderError, UffInversionContributionError,
    UffParameterError,
};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

fn project(node: ErrorRef) -> PyErr {
    macro_rules! retained {
        ($($name:ident),+ $(,)?)=>{$(
            if node.error().is::<ck::$name>() {return $name::from_ref(node);}
        )+};
    }
    retained!(
        GraphBoundsError,
        UffBoundsError,
        GenerationFailure,
        GenerationError,
        PruningError,
        BoundsQueryError,
        ConformerRunError,
        EnumerationError,
        EnumerationRunError,
        UffParameterError,
        MissingExplicitHydrogensError,
        ThreadCountError
    );
    if node.error().is::<std::io::Error>() {
        return IoError::from_ref(node);
    }
    if node.error().is::<std::sync::PoisonError<()>>() {
        return PoisonError::from_ref(node);
    }
    macro_rules! copied {
        ($($module:ident::$name:ident),+ $(,)?)=>{$(
            if let Some(value)=node.error().downcast_ref::<ck::$name>(){
                return crate::$module::$name::from_rust(value.clone());
            }
        )+};
    }
    copied!(
        conformer_execution::BoundsMatrixError,
        conformer_execution::EmbedSeedError,
        conformer_execution::CrystalffTorsionPreferencesError,
        conformer_execution::ConformerOptimizerError,
        model_errors::TopologyValidationError,
        model_errors::StereoGroupError,
        model_errors::CoordinateValidationError,
        model_edit_errors::AtomPropertyError,
        model_edit_errors::BondValueError,
        model_edit_errors::MoleculePropertyError,
        model_edit_errors::MappingValidationError,
        attachment_values::TemplateAttachmentOrderError,
        valence_errors::ValenceError,
        ring_errors::RingFindingError,
        canonical_errors::CanonicalRankError,
        canonical_errors::PropertyStringError,
        canonical_errors::RequiredPropertyStringError,
        sanitize_support_errors::ConjugationError,
        sanitize_leaf_errors::HybridizationError,
        stereo_order_errors::StereoOrderError,
        stereo_errors::AtropisomerError,
        stereo_errors::LegacyStereoError,
        double_stereo_errors::DoubleBondStereoError,
        potential_stereo_errors::PotentialStereoError,
        topology_edit_errors::TopologyEditError,
        property_numeric_errors::PropertyIntReadError,
        property_numeric_errors::PropertyUIntReadError,
        property_numeric_errors::PropertyULongReadError,
        property_values::PropertyValueError,
        smiles_stereo_errors::SmilesStereoError,
        parser_errors::ParserCarrierError,
        parser_errors::ParserStereoOrderError,
        hydrogen_errors::HydrogenError,
        kekulize_errors::KekulizeError,
        sanitize_errors::SanitizeError,
        matrix_errors::MatrixError,
        embedding_values::ConformerError,
        uff_parameter_errors::UffBuilderError,
        uff_parameter_errors::ForceFieldKernelError,
        uff_parameter_errors::UffInversionContributionError,
        uff_parameter_errors::UffBondMathError,
        match_errors::SubstructMatchError,
        match_errors::QueryMatchContextError,
        string_errors::QueryCompileError,
        smarts_errors::SmartsParseError,
        query_model_errors::QueryGraphError,
        fragment_errors::MoleculeFragmentsError,
        fragment_errors::MoleculeFragmentsFailure
    );
    unreachable!("unmodeled concrete source in native conformer/enumeration graph")
}

enum_value!(MatrixAxis,{Row,Column});
enum_value!(BoundKind,{Upper,Lower});
error_projection!(BoundsMatrixError,BoundsMatrixErrorKind,ck::BoundsMatrixError,"conformer",{
    DimensionOverflow {dimension},AllocationFailed {elements:usize},DataLengthTooShort {expected:usize,actual:usize},
    IndexOutOfBounds {axis,index:usize,dimension},InvalidBound {kind,row:usize,column:usize},
},cause=|_:&ck::BoundsMatrixError|None::<PyErr>,methods={
    #[getter]
    #[doc = "Coordinate or matrix dimension involved in the failure. Applies to BoundsMatrixError.DimensionOverflow, BoundsMatrixError.IndexOutOfBounds. None when this payload is absent."] fn dimension(&self)->Option<usize>{match self.inner{ck::BoundsMatrixError::DimensionOverflow {dimension}|ck::BoundsMatrixError::IndexOutOfBounds {dimension,..}=>Some(dimension),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Coordinate axis involved in the failure. Applies to BoundsMatrixError.IndexOutOfBounds. None when this payload is absent."] fn axis(&self)->Option<MatrixAxis>{match self.inner{ck::BoundsMatrixError::IndexOutOfBounds {axis,..}=>Some(axis.into()),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to BoundsMatrixError.InvalidBound. None when this payload is absent."] fn kind(&self)->Option<BoundKind>{match self.inner{ck::BoundsMatrixError::InvalidBound {kind,..}=>Some(kind.into()),_=>None}}
});
error_projection!(EmbedSeedError,EmbedSeedErrorKind,ck::EmbedSeedError,"conformer",{
    InvalidRandomSeed {seed},ResultInvariant {seed,index:usize,value:i32},
},cause=|_:&ck::EmbedSeedError|None::<PyErr>,methods={
    #[getter]
    #[doc = "Random seed associated with the failed embedding seed validation."] fn seed(&self)->i32{match self.inner{ck::EmbedSeedError::InvalidRandomSeed {seed}|ck::EmbedSeedError::ResultInvariant {seed,..}=>seed}}
});

type BoundsSource = Source<
    crate::error_framework::ErrorObject<BoundsMatrixError>,
    Source<
        crate::error_framework::ErrorObject<UffBoundsError>,
        Source<
            crate::error_framework::ErrorObject<MatrixError>,
            crate::error_framework::ErrorObject<ValenceError>,
        >,
    >,
>;
error_projection!(retained GraphBoundsError,GraphBoundsErrorKind,ck::GraphBoundsError,"conformer",{
    Bounds(_),Uff(_),Matrix(_),Valence(_),InvalidBounds(_),Detail(_),Input(_),AtomicNumber(_),
},cause=|inner:&RetainedError<ck::GraphBoundsError>|inner.node().source().map(project),source=Option<BoundsSource>,methods={
    #[getter]
    #[doc = "Detailed native diagnostic message. Applies to GraphBoundsError.InvalidBounds, GraphBoundsError.Detail, GraphBoundsError.Input. None when this payload is absent."] fn message(&self)->Option<&str>{match self.inner.as_ref(){ck::GraphBoundsError::InvalidBounds(value)|ck::GraphBoundsError::Detail(value)=>Some(value),ck::GraphBoundsError::Input(value)=>Some(value),_=>None}}
    #[getter]
    #[doc = "Atomic number involved in the failed chemistry operation. Applies to GraphBoundsError.AtomicNumber. None when this payload is absent."] fn atomic_number(&self)->Option<u8>{match self.inner.as_ref(){ck::GraphBoundsError::AtomicNumber(value)=>Some(*value),_=>None}}
});
type PreparationSource = Source<
    crate::error_framework::ErrorObject<UffParameterError>,
    Source<
        crate::error_framework::ErrorObject<ValenceError>,
        crate::error_framework::ErrorObject<UffBondMathError>,
    >,
>;
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct UffBoundsError {
    inner: RetainedError<ck::UffBoundsError>,
}
}
impl UffBoundsError {
    pub(crate) fn from_rust(inner: ck::UffBoundsError) -> PyErr {
        Self::from_ref(RetainedError::new(inner).node().clone())
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let cause = node.source().map(project);
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl UffBoundsError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured UFF bounds-construction failure details."] fn cause(&self) -> UffBoundsCause {
        UffBoundsCause {
            inner: self.inner.clone(),
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<PreparationSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
crate::enums::binding_enum!(UffBoundsCauseKind, {
    Parameter,
    BondOrder,
    BondMath,
    AssignmentLength,
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
#[doc = "Structured cause of UFF bounds construction failure."]
pub struct UffBoundsCause {
    inner: RetainedError<ck::UffBoundsError>,
}

crate::enums::binding_enum!(DistGeomForceFieldErrorKind, {
    Kernel,
    Inversion,
    IndexOverflow,
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
#[doc = "Distance-geometry force-field construction diagnostic."]
pub struct DistGeomForceFieldError {
    inner: ck::DistGeomForceFieldError,
}
#[pymethods]
impl DistGeomForceFieldError {
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
    fn variant(&self) -> DistGeomForceFieldErrorKind {
        match self.inner {
            ck::DistGeomForceFieldError::Kernel(_) => DistGeomForceFieldErrorKind::Kernel,
            ck::DistGeomForceFieldError::Inversion(_) => DistGeomForceFieldErrorKind::Inversion,
            ck::DistGeomForceFieldError::IndexOverflow(_) => {
                DistGeomForceFieldErrorKind::IndexOverflow
            }
        }
    }
    #[getter]
    #[doc = "Index reported by the failing native check. Applies to DistGeomForceFieldError.IndexOverflow. None when this payload is absent."]
    fn index(&self) -> Option<usize> {
        match self.inner {
            ck::DistGeomForceFieldError::IndexOverflow(value) => Some(value),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Structured native error payload. Present for `DistGeomForceFieldError.Kernel`; None when unavailable."]
    fn kernel(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<ForceFieldKernelError>>> {
        match self.inner {
            ck::DistGeomForceFieldError::Kernel(value) => ForceFieldKernelError::from_rust(value)
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = "Structured native error payload. Present for `DistGeomForceFieldError.Inversion`; None when unavailable."]
    fn inversion(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<UffInversionContributionError>>> {
        match self.inner {
            ck::DistGeomForceFieldError::Inversion(value) => {
                UffInversionContributionError::from_rust(value)
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
}
crate::enums::binding_enum!(ConformerOptimizerCauseKind, {
    Input,
    Construction,
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
#[doc = "Structured cause of conformer optimization failure, preserving its native variant."]
pub struct ConformerOptimizerCause {
    inner: ck::ConformerOptimizerCause,
}
#[pymethods]
impl ConformerOptimizerCause {
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
    fn variant(&self) -> ConformerOptimizerCauseKind {
        match self.inner {
            ck::ConformerOptimizerCause::Input(_) => ConformerOptimizerCauseKind::Input,
            ck::ConformerOptimizerCause::Construction(_) => {
                ConformerOptimizerCauseKind::Construction
            }
        }
    }
    #[getter]
    #[doc = "Native diagnostic message. Present for `ConformerOptimizerCause.Input`; None when unavailable."]
    fn message(&self) -> Option<&'static str> {
        match self.inner {
            ck::ConformerOptimizerCause::Input(value) => Some(value),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Force-field construction details for Construction failures; None otherwise."]
    fn construction(&self) -> Option<DistGeomForceFieldError> {
        match self.inner {
            ck::ConformerOptimizerCause::Construction(inner) => {
                Some(DistGeomForceFieldError { inner })
            }
            _ => None,
        }
    }
}
type OptimizerSource = Source<
    crate::error_framework::ErrorObject<ForceFieldKernelError>,
    crate::error_framework::ErrorObject<UffInversionContributionError>,
>;
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct ConformerOptimizerError {
    inner: ck::ConformerOptimizerError,
}
}
impl ConformerOptimizerError {
    pub(crate) fn from_rust(inner: ck::ConformerOptimizerError) -> PyErr {
        let cause = inner.source().map(|source| {
            if let Some(value) = source.downcast_ref::<ck::ForceFieldKernelError>() {
                ForceFieldKernelError::from_rust(*value)
            } else if let Some(value) = source.downcast_ref::<ck::UffInversionContributionError>() {
                UffInversionContributionError::from_rust(*value)
            } else {
                unreachable!("native optimizer source")
            }
        });
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => {
                let error = crate::error_framework::from_payload(value);
                error.set_cause(py, cause);
                error
            }
            Err(error) => error,
        })
    }
}
crate::error_framework::error_projection! {methods impl ConformerOptimizerError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured native cause of conformer optimization failure."] fn cause(&self) -> ConformerOptimizerCause {
        ConformerOptimizerCause {
            inner: *self.inner.cause(),
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<OptimizerSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct MissingExplicitHydrogensError {
    inner: RetainedError<ck::MissingExplicitHydrogensError>,
}
}
impl MissingExplicitHydrogensError {
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let cause = node.source().map(project);
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl MissingExplicitHydrogensError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<crate::error_framework::ErrorObject<UffBuilderError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct ConformerRunError {
    inner: RetainedError<ck::ConformerRunError>,
}
}
impl ConformerRunError {
    pub(crate) fn from_rust(inner: ck::ConformerRunError) -> PyErr {
        Self::from_ref(RetainedError::new(inner).node().clone())
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let cause = node.source().map(project);
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl ConformerRunError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Source<crate::error_framework::ErrorObject<BoundsQueryError>, crate::error_framework::ErrorObject<GenerationError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct EnumerationRunError {
    inner: RetainedError<ck::EnumerationRunError>,
}
}
impl EnumerationRunError {
    pub(crate) fn from_rust(inner: ck::EnumerationRunError) -> PyErr {
        Self::from_ref(RetainedError::new(inner).node().clone())
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let cause = node.source().map(project);
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl EnumerationRunError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "stereo"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<crate::error_framework::ErrorObject<EnumerationError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
#[pymethods]
impl UffBoundsCause {
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
    fn variant(&self) -> UffBoundsCauseKind {
        match self.inner.cause() {
            ck::UffBoundsCause::Parameter(_) => UffBoundsCauseKind::Parameter,
            ck::UffBoundsCause::BondOrder(_) => UffBoundsCauseKind::BondOrder,
            ck::UffBoundsCause::BondMath(_) => UffBoundsCauseKind::BondMath,
            ck::UffBoundsCause::AssignmentLength {
                field: _,
                actual: _,
                expected: _,
            } => UffBoundsCauseKind::AssignmentLength,
        }
    }
    #[getter]
    #[doc = "Input field associated with the failure. Present for `UffBoundsCause.AssignmentLength`; None when unavailable."]
    fn field(&self) -> Option<&'static str> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { field, .. } => Some(*field),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Observed count or value. Present for `UffBoundsCause.AssignmentLength`; None when unavailable."]
    fn actual(&self) -> Option<usize> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { actual, .. } => Some(*actual),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Expected count or value. Present for `UffBoundsCause.AssignmentLength`; None when unavailable."]
    fn expected(&self) -> Option<usize> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { expected, .. } => Some(*expected),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Structured native error payload. Present for `UffBoundsCause.Parameter`; None when unavailable."]
    fn parameter(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<UffParameterError>>> {
        if matches!(self.inner.cause(), ck::UffBoundsCause::Parameter(_)) {
            UffParameterError::from_ref(
                self.inner
                    .node()
                    .source()
                    .expect("actual native UFF bounds parameter source"),
            )
            .value(py)
            .extract()
            .map(Some)
            .map_err(Into::into)
        } else {
            Ok(None)
        }
    }
    #[getter]
    #[doc = "Structured native error payload. Present for `UffBoundsCause.BondOrder`; None when unavailable."]
    fn bond_order(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> {
        match self.inner.cause() {
            ck::UffBoundsCause::BondOrder(value) => ValenceError::from_rust(value.clone())
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = "Structured native error payload. Present for `UffBoundsCause.BondMath`; None when unavailable."]
    fn bond_math(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<UffBondMathError>>> {
        match self.inner.cause() {
            ck::UffBoundsCause::BondMath(value) => UffBondMathError::from_rust(*value)
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
}

type PruningSource = Source<
    LegacySource,
    Source<
        MatchSource,
        Source<SmartsSource, crate::error_framework::ErrorObject<ThreadCountError>>,
    >,
>;
error_projection!(retained PruningError,PruningErrorKind,ck::PruningError,"conformer",{
    Hydrogens(_),Projection(_),Terminal(_),Match(_),Context(_),Input(_),
},cause=|inner:&RetainedError<ck::PruningError>|inner.node().source().map(project),source=Option<PruningSource>,methods={
    #[getter]
    #[doc = "Detailed native diagnostic message. Applies to PruningError.Input. None when this payload is absent."] fn message(&self)->Option<&'static str>{match self.inner.as_ref(){ck::PruningError::Input(value)=>Some(*value),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Hydrogen count or hydrogen-related payload reported by the native check. Applies to PruningError.Hydrogens. None when this payload is absent."] fn hydrogens(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::hydrogen_errors::HydrogenError>>>{match self.inner.as_ref(){ck::PruningError::Hydrogens(value)=>crate::hydrogen_errors::HydrogenError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `PruningError.Projection`; None when unavailable."] fn projection(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::query_model_errors::QueryGraphError>>>{match self.inner.as_ref(){ck::PruningError::Projection(value)=>crate::query_model_errors::QueryGraphError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `PruningError.Terminal`; None when unavailable."] fn terminal(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::alignment_errors::AlignmentError>>>{match self.inner.as_ref(){ck::PruningError::Terminal(value)=>crate::alignment_errors::AlignmentError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `PruningError.Match`; None when unavailable."] fn matching(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::match_errors::SubstructMatchError>>>{match self.inner.as_ref(){ck::PruningError::Match(value)=>crate::match_errors::SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `PruningError.Context`; None when unavailable."] fn context(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::match_errors::QueryMatchContextError>>>{match self.inner.as_ref(){ck::PruningError::Context(value)=>crate::match_errors::QueryMatchContextError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
});
error_projection!(CrystalffTorsionPreferencesError,CrystalffTorsionPreferencesErrorKind,ck::CrystalffTorsionPreferencesError,"conformer",{
    InvalidExperimentalTorsionVersion {version:u32},MissingSmartsPattern {line},IncompleteParameterLine {line},
    InvalidIntegerToken {token,line},InvalidFloatToken {token,line},InvalidParameterSourceFormat {constant_name},
    UnsupportedEscapeSequence {constant_name,escape:String=>|value:&char|value.to_string()},
    UnterminatedQuotedLiteral {constant_name},TrailingEscape {constant_name},EmptyMolecule,PreparedTarget(_),Match(_),
    QueryBuildFailed {smarts,reason},RingInfoUnavailable {reason},MissingCentralBond {smarts,aid2:usize,aid3:usize},
},cause=|inner:&ck::CrystalffTorsionPreferencesError|RetainedError::new(inner.clone()).node().source().map(project),source=Option<MatchSource>,methods={
    #[getter]
    #[doc = "Source line involved in the failure. Applies to CrystalffTorsionPreferencesError.MissingSmartsPattern, CrystalffTorsionPreferencesError.IncompleteParameterLine, CrystalffTorsionPreferencesError.InvalidIntegerToken, CrystalffTorsionPreferencesError.InvalidFloatToken. None when this payload is absent."] fn line(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::MissingSmartsPattern {line}|ck::CrystalffTorsionPreferencesError::IncompleteParameterLine {line}|ck::CrystalffTorsionPreferencesError::InvalidIntegerToken {line,..}|ck::CrystalffTorsionPreferencesError::InvalidFloatToken {line,..}=>Some(line),_=>None}}
    #[getter]
    #[doc = "Input token rejected by the parser. Applies to CrystalffTorsionPreferencesError.InvalidIntegerToken, CrystalffTorsionPreferencesError.InvalidFloatToken. None when this payload is absent."] fn token(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::InvalidIntegerToken {token,..}|ck::CrystalffTorsionPreferencesError::InvalidFloatToken {token,..}=>Some(token),_=>None}}
    #[getter]
    #[doc = "Named parameter-table constant involved in the failure. Applies to CrystalffTorsionPreferencesError.InvalidParameterSourceFormat, CrystalffTorsionPreferencesError.UnsupportedEscapeSequence, CrystalffTorsionPreferencesError.UnterminatedQuotedLiteral, CrystalffTorsionPreferencesError.TrailingEscape. None when this payload is absent."] fn constant_name(&self)->Option<&'static str>{match self.inner{ck::CrystalffTorsionPreferencesError::InvalidParameterSourceFormat {constant_name}|ck::CrystalffTorsionPreferencesError::UnsupportedEscapeSequence {constant_name,..}|ck::CrystalffTorsionPreferencesError::UnterminatedQuotedLiteral {constant_name}|ck::CrystalffTorsionPreferencesError::TrailingEscape {constant_name}=>Some(constant_name),_=>None}}
    #[getter]
    #[doc = "SMARTS expression implicated in the parsing or query failure. Applies to CrystalffTorsionPreferencesError.QueryBuildFailed, CrystalffTorsionPreferencesError.MissingCentralBond. None when this payload is absent."] fn smarts(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::QueryBuildFailed {smarts,..}|ck::CrystalffTorsionPreferencesError::MissingCentralBond {smarts,..}=>Some(smarts),_=>None}}
    #[getter]
    #[doc = "Native explanation for rejecting the operation. Applies to CrystalffTorsionPreferencesError.QueryBuildFailed, CrystalffTorsionPreferencesError.RingInfoUnavailable. None when this payload is absent."] fn reason(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::QueryBuildFailed {reason,..}|ck::CrystalffTorsionPreferencesError::RingInfoUnavailable {reason}=>Some(reason),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `CrystalffTorsionPreferencesError.PreparedTarget`; None when unavailable."] fn prepared_target(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::match_errors::QueryMatchContextError>>>{match &self.inner{ck::CrystalffTorsionPreferencesError::PreparedTarget(value)=>crate::match_errors::QueryMatchContextError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `CrystalffTorsionPreferencesError.Match`; None when unavailable."] fn matching(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::match_errors::SubstructMatchError>>>{match &self.inner{ck::CrystalffTorsionPreferencesError::Match(value)=>crate::match_errors::SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
});

type GenerationNumericSource = Source<
    BoundsSource,
    Source<
        crate::error_framework::ErrorObject<ForceFieldKernelError>,
        crate::error_framework::ErrorObject<UffInversionContributionError>,
    >,
>;
type GenerationPreparationSource = Source<
    PruningSource,
    Source<
        crate::error_framework::ErrorObject<MoleculeFragmentsFailure>,
        Source<crate::error_framework::ErrorObject<UffBuilderError>, AtropSource>,
    >,
>;
pub(crate) type GenerationSource = Source<
    Source<SmilesSource, LegacySource>,
    Source<
        GenerationNumericSource,
        Source<
            GenerationPreparationSource,
            Source<MatchSource, crate::error_framework::ErrorObject<IoError>>,
        >,
    >,
>;
error_projection!(retained GenerationFailure,GenerationFailureKind,ck::GenerationFailure,"conformer",{
    TopologyBounds(_),TorsionPreferences(_),
},cause=|inner:&RetainedError<ck::GenerationFailure>|inner.node().source().map(project),source=Option<Source<BoundsSource,MatchSource>>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured GraphBoundsError payload of GenerationFailure.TopologyBounds; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn topology_bounds(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<GraphBoundsError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationFailure::TopologyBounds(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationFailure::TopologyBounds(value)=parent.downcast_ref::<ck::GenerationFailure>().expect("actual generation failure owner") else{unreachable!()};value});
            GraphBoundsError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured CrystalffTorsionPreferencesError payload of GenerationFailure.TorsionPreferences; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn torsion_preferences(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<CrystalffTorsionPreferencesError>>>{match self.inner.as_ref(){ck::GenerationFailure::TorsionPreferences(value)=>CrystalffTorsionPreferencesError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
});
error_projection!(retained BoundsQueryError,BoundsQueryErrorKind,ck::BoundsQueryError,"conformer",{
    Topology(_),Valence(_),Rings(_),Conjugation(_),Hybridization(_),Bounds(_),
},cause=|inner:&RetainedError<ck::BoundsQueryError>|inner.node().source().map(project),source=Option<Source<KekulizeSource,BoundsSource>>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `BoundsQueryError.Topology`; None when unavailable."] fn topology(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_errors::TopologyValidationError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Topology(value)=>crate::model_errors::TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of BoundsQueryError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `BoundsQueryError.Rings`; None when unavailable."] fn rings(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::ring_errors::RingFindingError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Rings(value)=>crate::ring_errors::RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `BoundsQueryError.Conjugation`; None when unavailable."] fn conjugation(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::sanitize_support_errors::ConjugationError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Conjugation(value)=>crate::sanitize_support_errors::ConjugationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `BoundsQueryError.Hybridization`; None when unavailable."] fn hybridization(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::sanitize_leaf_errors::HybridizationError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Hybridization(value)=>crate::sanitize_leaf_errors::HybridizationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured GraphBoundsError payload of BoundsQueryError.Bounds; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn bounds(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<GraphBoundsError>>>{
        if matches!(self.inner.as_ref(),ck::BoundsQueryError::Bounds(_)){
            let node=self.inner.node().project(|parent|{let ck::BoundsQueryError::Bounds(value)=parent.downcast_ref::<ck::BoundsQueryError>().expect("actual bounds query owner") else{unreachable!()};value});
            GraphBoundsError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
});
error_projection!(retained GenerationError,GenerationErrorKind,ck::GenerationError,"conformer",{
    Topology(_),Coordinates(_),Valence(_),Rings(_),Conjugation(_),Hybridization(_),Stereo(_),Fragments(_),
    Hydrogens(_),Pruning(_),Seed(_),GenerationFailed(_),Atropisomer(_),Bounds(_),ForceField(_),Numeric(_),
    ThreadCount(_),WorkerLaunch(_),Input(_),CoordinateMapIndex {atom_index,atom_count:usize},
    CoordinateMapNonFinite {atom_index},UnsupportedImplementation {implementation:&'static str},
},cause=|inner:&RetainedError<ck::GenerationError>|inner.node().source().map(project),source=Option<GenerationSource>,methods={
    #[getter]
    #[doc = "Atom index involved in the native failure. Applies to GenerationError.CoordinateMapIndex, GenerationError.CoordinateMapNonFinite. None when this payload is absent."] fn atom_index(&self)->Option<i32>{match self.inner.as_ref(){ck::GenerationError::CoordinateMapIndex {atom_index,..}|ck::GenerationError::CoordinateMapNonFinite {atom_index}=>Some(*atom_index),_=>None}}
    #[getter]
    #[doc = "Detailed native diagnostic message. Applies to GenerationError.Input. None when this payload is absent."] fn message(&self)->Option<&'static str>{match self.inner.as_ref(){ck::GenerationError::Input(value)=>Some(*value),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Topology`; None when unavailable."] fn topology(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_errors::TopologyValidationError>>>{match self.inner.as_ref(){ck::GenerationError::Topology(value)=>crate::model_errors::TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Coordinates`; None when unavailable."] fn coordinates(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_errors::CoordinateValidationError>>>{match self.inner.as_ref(){ck::GenerationError::Coordinates(value)=>crate::model_errors::CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of GenerationError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>>{match self.inner.as_ref(){ck::GenerationError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Rings`; None when unavailable."] fn rings(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::ring_errors::RingFindingError>>>{match self.inner.as_ref(){ck::GenerationError::Rings(value)=>crate::ring_errors::RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Conjugation`; None when unavailable."] fn conjugation(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::sanitize_support_errors::ConjugationError>>>{match self.inner.as_ref(){ck::GenerationError::Conjugation(value)=>crate::sanitize_support_errors::ConjugationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Hybridization`; None when unavailable."] fn hybridization(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::sanitize_leaf_errors::HybridizationError>>>{match self.inner.as_ref(){ck::GenerationError::Hybridization(value)=>crate::sanitize_leaf_errors::HybridizationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Stereo`; None when unavailable."] fn stereo(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::stereo_errors::LegacyStereoError>>>{match self.inner.as_ref(){ck::GenerationError::Stereo(value)=>crate::stereo_errors::LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Fragments`; None when unavailable."] fn fragments(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::fragment_errors::MoleculeFragmentsError>>>{match self.inner.as_ref(){ck::GenerationError::Fragments(value)=>crate::fragment_errors::MoleculeFragmentsError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured EmbedSeedError payload of GenerationError.Seed; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn seed(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<EmbedSeedError>>>{match self.inner.as_ref(){ck::GenerationError::Seed(value)=>EmbedSeedError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `GenerationError.Atropisomer`; None when unavailable."] fn atropisomer(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::stereo_errors::AtropisomerError>>>{match self.inner.as_ref(){ck::GenerationError::Atropisomer(value)=>crate::stereo_errors::AtropisomerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured BoundsMatrixError payload of GenerationError.Bounds; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn bounds(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<BoundsMatrixError>>>{match self.inner.as_ref(){ck::GenerationError::Bounds(value)=>BoundsMatrixError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured ConformerOptimizerError payload of GenerationError.ForceField; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn force_field(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ConformerOptimizerError>>>{match self.inner.as_ref(){ck::GenerationError::ForceField(value)=>ConformerOptimizerError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured ConformerError payload of GenerationError.Numeric; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn numeric(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ConformerError>>>{match self.inner.as_ref(){ck::GenerationError::Numeric(value)=>ConformerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured MissingExplicitHydrogensError payload of GenerationError.Hydrogens; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn hydrogens(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<MissingExplicitHydrogensError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::Hydrogens(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::Hydrogens(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            MissingExplicitHydrogensError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured PruningError payload of GenerationError.Pruning; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn pruning(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PruningError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::Pruning(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::Pruning(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            PruningError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured GenerationFailure payload of GenerationError.GenerationFailed; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn generation_failed(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<GenerationFailure>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::GenerationFailed(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::GenerationFailed(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            GenerationFailure::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured ThreadCountError payload of GenerationError.ThreadCount; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn thread_count(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ThreadCountError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::ThreadCount(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::ThreadCount(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            ThreadCountError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
});

// The concrete standard poison value carries unit data and no source/kind.
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_leaf"}"#]
pub struct PoisonError {
    inner: RetainedError<std::sync::PoisonError<()>>,
}
}
impl PoisonError {
    fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => crate::error_framework::from_payload(value),
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl PoisonError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "standard"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

pub(crate) type EnumerationSource = Source<
    GenerationSource,
    Source<AtropSource, Source<crate::error_framework::ErrorObject<PoisonError>, SmilesSource>>,
>;
error_projection!(retained EnumerationError,EnumerationErrorKind,ck::EnumerationError,"stereo",{
    Valence(_),Rings(_),Stereo(_),DoubleBond(_),SmilesWrite(_),PotentialStereo(_),Hydrogens(_),Embedding(_),
    AtomProperty(_),BondProperty(_),MoleculeProperty(_),Topology(_),Coordinates(_),
    InvalidFlipperAtom {atom:usize=>|value:&ck::AtomId|value.index()},InvalidFlipperBond {bond:usize=>|value:&ck::BondId|value.index()},
    RandomBitsSource {bit_count:usize,message:String},RandomSourcePoisoned(_),EmbeddingCoordinates,
},cause=|inner:&RetainedError<ck::EnumerationError>|inner.node().source().map(project),source=Option<EnumerationSource>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of EnumerationError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>>{match self.inner.as_ref(){ck::EnumerationError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.Rings`; None when unavailable."] fn rings(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::ring_errors::RingFindingError>>>{match self.inner.as_ref(){ck::EnumerationError::Rings(value)=>crate::ring_errors::RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.Stereo`; None when unavailable."] fn stereo(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::stereo_errors::LegacyStereoError>>>{match self.inner.as_ref(){ck::EnumerationError::Stereo(value)=>crate::stereo_errors::LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.DoubleBond`; None when unavailable."] fn double_bond(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::double_stereo_errors::DoubleBondStereoError>>>{match self.inner.as_ref(){ck::EnumerationError::DoubleBond(value)=>crate::double_stereo_errors::DoubleBondStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.SmilesWrite`; None when unavailable."] fn smiles_write(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::smiles_parse_errors::SmilesParseError>>>{match self.inner.as_ref(){ck::EnumerationError::SmilesWrite(value)=>crate::smiles_parse_errors::SmilesParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.PotentialStereo`; None when unavailable."] fn potential_stereo(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::potential_stereo_errors::PotentialStereoError>>>{match self.inner.as_ref(){ck::EnumerationError::PotentialStereo(value)=>crate::potential_stereo_errors::PotentialStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Hydrogen count or hydrogen-related payload reported by the native check. Applies to EnumerationError.Hydrogens. None when this payload is absent."] fn hydrogens(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::hydrogen_errors::HydrogenError>>>{match self.inner.as_ref(){ck::EnumerationError::Hydrogens(value)=>crate::hydrogen_errors::HydrogenError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.AtomProperty`; None when unavailable."] fn atom_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_edit_errors::AtomPropertyError>>>{match self.inner.as_ref(){ck::EnumerationError::AtomProperty(value)=>crate::model_edit_errors::AtomPropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.BondProperty`; None when unavailable."] fn bond_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_edit_errors::BondValueError>>>{match self.inner.as_ref(){ck::EnumerationError::BondProperty(value)=>crate::model_edit_errors::BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.MoleculeProperty`; None when unavailable."] fn molecule_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_edit_errors::MoleculePropertyError>>>{match self.inner.as_ref(){ck::EnumerationError::MoleculeProperty(value)=>crate::model_edit_errors::MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.Topology`; None when unavailable."] fn topology(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_errors::TopologyValidationError>>>{match self.inner.as_ref(){ck::EnumerationError::Topology(value)=>crate::model_errors::TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured native error payload. Present for `EnumerationError.Coordinates`; None when unavailable."] fn coordinates(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<crate::model_errors::CoordinateValidationError>>>{match self.inner.as_ref(){ck::EnumerationError::Coordinates(value)=>crate::model_errors::CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    #[doc = "Structured GenerationError payload of EnumerationError.Embedding; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn embedding(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<GenerationError>>>{
        if matches!(self.inner.as_ref(),ck::EnumerationError::Embedding(_)){
            let node=self.inner.node().project(|parent|{let ck::EnumerationError::Embedding(value)=parent.downcast_ref::<ck::EnumerationError>().expect("actual enumeration owner")else{unreachable!()};value});
            GenerationError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
});
