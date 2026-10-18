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
use pyo3::exceptions::PyValueError;
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
    fn dimension(&self)->Option<usize>{match self.inner{ck::BoundsMatrixError::DimensionOverflow {dimension}|ck::BoundsMatrixError::IndexOutOfBounds {dimension,..}=>Some(dimension),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn axis(&self)->Option<MatrixAxis>{match self.inner{ck::BoundsMatrixError::IndexOutOfBounds {axis,..}=>Some(axis.into()),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn kind(&self)->Option<BoundKind>{match self.inner{ck::BoundsMatrixError::InvalidBound {kind,..}=>Some(kind.into()),_=>None}}
});
error_projection!(EmbedSeedError,EmbedSeedErrorKind,ck::EmbedSeedError,"conformer",{
    InvalidRandomSeed {seed},ResultInvariant {seed,index:usize,value:i32},
},cause=|_:&ck::EmbedSeedError|None::<PyErr>,methods={
    #[getter]
    fn seed(&self)->i32{match self.inner{ck::EmbedSeedError::InvalidRandomSeed {seed}|ck::EmbedSeedError::ResultInvariant {seed,..}=>seed}}
});

type BoundsSource = Source<
    Py<BoundsMatrixError>,
    Source<Py<UffBoundsError>, Source<Py<MatrixError>, Py<ValenceError>>>,
>;
error_projection!(retained GraphBoundsError,GraphBoundsErrorKind,ck::GraphBoundsError,"conformer",{
    Bounds(_),Uff(_),Matrix(_),Valence(_),InvalidBounds(_),Detail(_),Input(_),AtomicNumber(_),
},cause=|inner:&RetainedError<ck::GraphBoundsError>|inner.node().source().map(project),source=Option<BoundsSource>,methods={
    #[getter]
    fn message(&self)->Option<&str>{match self.inner.as_ref(){ck::GraphBoundsError::InvalidBounds(value)|ck::GraphBoundsError::Detail(value)=>Some(value),ck::GraphBoundsError::Input(value)=>Some(value),_=>None}}
    #[getter]
    fn atomic_number(&self)->Option<u8>{match self.inner.as_ref(){ck::GraphBoundsError::AtomicNumber(value)=>Some(*value),_=>None}}
});
type PreparationSource =
    Source<Py<UffParameterError>, Source<Py<ValenceError>, Py<UffBondMathError>>>;
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct UffBoundsError {
    inner: RetainedError<ck::UffBoundsError>,
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
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl UffBoundsError {
    #[getter]
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn cause(&self) -> UffBoundsCause {
        UffBoundsCause {
            inner: self.inner.clone(),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<PreparationSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum UffBoundsCauseKind {
    Parameter,
    BondOrder,
    BondMath,
    AssignmentLength,
}
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct UffBoundsCause {
    inner: RetainedError<ck::UffBoundsError>,
}

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum DistGeomForceFieldErrorKind {
    Kernel,
    Inversion,
    IndexOverflow,
}
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct DistGeomForceFieldError {
    inner: ck::DistGeomForceFieldError,
}
#[pymethods]
impl DistGeomForceFieldError {
    #[getter]
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
    fn index(&self) -> Option<usize> {
        match self.inner {
            ck::DistGeomForceFieldError::IndexOverflow(value) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn kernel(&self, py: Python<'_>) -> PyResult<Option<Py<ForceFieldKernelError>>> {
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
    fn inversion(&self, py: Python<'_>) -> PyResult<Option<Py<UffInversionContributionError>>> {
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
#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum ConformerOptimizerCauseKind {
    Input,
    Construction,
}
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct ConformerOptimizerCause {
    inner: ck::ConformerOptimizerCause,
}
#[pymethods]
impl ConformerOptimizerCause {
    #[getter]
    fn variant(&self) -> ConformerOptimizerCauseKind {
        match self.inner {
            ck::ConformerOptimizerCause::Input(_) => ConformerOptimizerCauseKind::Input,
            ck::ConformerOptimizerCause::Construction(_) => {
                ConformerOptimizerCauseKind::Construction
            }
        }
    }
    #[getter]
    fn message(&self) -> Option<&'static str> {
        match self.inner {
            ck::ConformerOptimizerCause::Input(value) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn construction(&self) -> Option<DistGeomForceFieldError> {
        match self.inner {
            ck::ConformerOptimizerCause::Construction(inner) => {
                Some(DistGeomForceFieldError { inner })
            }
            _ => None,
        }
    }
}
type OptimizerSource = Source<Py<ForceFieldKernelError>, Py<UffInversionContributionError>>;
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct ConformerOptimizerError {
    inner: ck::ConformerOptimizerError,
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
                let error = PyErr::from_value(value.into_bound(py).into_any());
                error.set_cause(py, cause);
                error
            }
            Err(error) => error,
        })
    }
}
#[pymethods]
impl ConformerOptimizerError {
    #[getter]
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn cause(&self) -> ConformerOptimizerCause {
        ConformerOptimizerCause {
            inner: *self.inner.cause(),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<OptimizerSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct MissingExplicitHydrogensError {
    inner: RetainedError<ck::MissingExplicitHydrogensError>,
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
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl MissingExplicitHydrogensError {
    #[getter]
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Py<UffBuilderError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct ConformerRunError {
    inner: RetainedError<ck::ConformerRunError>,
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
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl ConformerRunError {
    #[getter]
    fn domain(&self) -> &'static str {
        "conformer"
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Source<Py<BoundsQueryError>, Py<GenerationError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct EnumerationRunError {
    inner: RetainedError<ck::EnumerationRunError>,
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
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl EnumerationRunError {
    #[getter]
    fn domain(&self) -> &'static str {
        "stereo"
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Py<EnumerationError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
#[pymethods]
impl UffBoundsCause {
    #[getter]
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
    fn field(&self) -> Option<&'static str> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { field, .. } => Some(*field),
            _ => None,
        }
    }
    #[getter]
    fn actual(&self) -> Option<usize> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { actual, .. } => Some(*actual),
            _ => None,
        }
    }
    #[getter]
    fn expected(&self) -> Option<usize> {
        match self.inner.cause() {
            ck::UffBoundsCause::AssignmentLength { expected, .. } => Some(*expected),
            _ => None,
        }
    }
    #[getter]
    fn parameter(&self, py: Python<'_>) -> PyResult<Option<Py<UffParameterError>>> {
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
    fn bond_order(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> {
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
    fn bond_math(&self, py: Python<'_>) -> PyResult<Option<Py<UffBondMathError>>> {
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

type PruningSource =
    Source<LegacySource, Source<MatchSource, Source<SmartsSource, Py<ThreadCountError>>>>;
error_projection!(retained PruningError,PruningErrorKind,ck::PruningError,"conformer",{
    Hydrogens(_),Projection(_),Terminal(_),Match(_),Context(_),Input(_),
},cause=|inner:&RetainedError<ck::PruningError>|inner.node().source().map(project),source=Option<PruningSource>,methods={
    #[getter]
    fn message(&self)->Option<&'static str>{match self.inner.as_ref(){ck::PruningError::Input(value)=>Some(*value),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn hydrogens(&self,py:Python<'_>)->PyResult<Option<Py<crate::hydrogen_errors::HydrogenError>>>{match self.inner.as_ref(){ck::PruningError::Hydrogens(value)=>crate::hydrogen_errors::HydrogenError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn projection(&self,py:Python<'_>)->PyResult<Option<Py<crate::query_model_errors::QueryGraphError>>>{match self.inner.as_ref(){ck::PruningError::Projection(value)=>crate::query_model_errors::QueryGraphError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn terminal(&self,py:Python<'_>)->PyResult<Option<Py<crate::alignment_errors::AlignmentError>>>{match self.inner.as_ref(){ck::PruningError::Terminal(value)=>crate::alignment_errors::AlignmentError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn matching(&self,py:Python<'_>)->PyResult<Option<Py<crate::match_errors::SubstructMatchError>>>{match self.inner.as_ref(){ck::PruningError::Match(value)=>crate::match_errors::SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn context(&self,py:Python<'_>)->PyResult<Option<Py<crate::match_errors::QueryMatchContextError>>>{match self.inner.as_ref(){ck::PruningError::Context(value)=>crate::match_errors::QueryMatchContextError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
});
error_projection!(CrystalffTorsionPreferencesError,CrystalffTorsionPreferencesErrorKind,ck::CrystalffTorsionPreferencesError,"conformer",{
    InvalidExperimentalTorsionVersion {version:u32},MissingSmartsPattern {line},IncompleteParameterLine {line},
    InvalidIntegerToken {token,line},InvalidFloatToken {token,line},InvalidParameterSourceFormat {constant_name},
    UnsupportedEscapeSequence {constant_name,escape:String=>|value:&char|value.to_string()},
    UnterminatedQuotedLiteral {constant_name},TrailingEscape {constant_name},EmptyMolecule,PreparedTarget(_),Match(_),
    QueryBuildFailed {smarts,reason},RingInfoUnavailable {reason},MissingCentralBond {smarts,aid2:usize,aid3:usize},
},cause=|inner:&ck::CrystalffTorsionPreferencesError|RetainedError::new(inner.clone()).node().source().map(project),source=Option<MatchSource>,methods={
    #[getter]
    fn line(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::MissingSmartsPattern {line}|ck::CrystalffTorsionPreferencesError::IncompleteParameterLine {line}|ck::CrystalffTorsionPreferencesError::InvalidIntegerToken {line,..}|ck::CrystalffTorsionPreferencesError::InvalidFloatToken {line,..}=>Some(line),_=>None}}
    #[getter]
    fn token(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::InvalidIntegerToken {token,..}|ck::CrystalffTorsionPreferencesError::InvalidFloatToken {token,..}=>Some(token),_=>None}}
    #[getter]
    fn constant_name(&self)->Option<&'static str>{match self.inner{ck::CrystalffTorsionPreferencesError::InvalidParameterSourceFormat {constant_name}|ck::CrystalffTorsionPreferencesError::UnsupportedEscapeSequence {constant_name,..}|ck::CrystalffTorsionPreferencesError::UnterminatedQuotedLiteral {constant_name}|ck::CrystalffTorsionPreferencesError::TrailingEscape {constant_name}=>Some(constant_name),_=>None}}
    #[getter]
    fn smarts(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::QueryBuildFailed {smarts,..}|ck::CrystalffTorsionPreferencesError::MissingCentralBond {smarts,..}=>Some(smarts),_=>None}}
    #[getter]
    fn reason(&self)->Option<&str>{match &self.inner{ck::CrystalffTorsionPreferencesError::QueryBuildFailed {reason,..}|ck::CrystalffTorsionPreferencesError::RingInfoUnavailable {reason}=>Some(reason),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn prepared_target(&self,py:Python<'_>)->PyResult<Option<Py<crate::match_errors::QueryMatchContextError>>>{match &self.inner{ck::CrystalffTorsionPreferencesError::PreparedTarget(value)=>crate::match_errors::QueryMatchContextError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn matching(&self,py:Python<'_>)->PyResult<Option<Py<crate::match_errors::SubstructMatchError>>>{match &self.inner{ck::CrystalffTorsionPreferencesError::Match(value)=>crate::match_errors::SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
});

type GenerationNumericSource =
    Source<BoundsSource, Source<Py<ForceFieldKernelError>, Py<UffInversionContributionError>>>;
type GenerationPreparationSource = Source<
    PruningSource,
    Source<Py<MoleculeFragmentsFailure>, Source<Py<UffBuilderError>, AtropSource>>,
>;
pub(crate) type GenerationSource = Source<
    Source<SmilesSource, LegacySource>,
    Source<
        GenerationNumericSource,
        Source<GenerationPreparationSource, Source<MatchSource, Py<IoError>>>,
    >,
>;
error_projection!(retained GenerationFailure,GenerationFailureKind,ck::GenerationFailure,"conformer",{
    TopologyBounds(_),TorsionPreferences(_),
},cause=|inner:&RetainedError<ck::GenerationFailure>|inner.node().source().map(project),source=Option<Source<BoundsSource,MatchSource>>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn topology_bounds(&self,py:Python<'_>)->PyResult<Option<Py<GraphBoundsError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationFailure::TopologyBounds(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationFailure::TopologyBounds(value)=parent.downcast_ref::<ck::GenerationFailure>().expect("actual generation failure owner") else{unreachable!()};value});
            GraphBoundsError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn torsion_preferences(&self,py:Python<'_>)->PyResult<Option<Py<CrystalffTorsionPreferencesError>>>{match self.inner.as_ref(){ck::GenerationFailure::TorsionPreferences(value)=>CrystalffTorsionPreferencesError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
});
error_projection!(retained BoundsQueryError,BoundsQueryErrorKind,ck::BoundsQueryError,"conformer",{
    Topology(_),Valence(_),Rings(_),Conjugation(_),Hybridization(_),Bounds(_),
},cause=|inner:&RetainedError<ck::BoundsQueryError>|inner.node().source().map(project),source=Option<Source<KekulizeSource,BoundsSource>>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn topology(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_errors::TopologyValidationError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Topology(value)=>crate::model_errors::TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn valence(&self,py:Python<'_>)->PyResult<Option<Py<ValenceError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn rings(&self,py:Python<'_>)->PyResult<Option<Py<crate::ring_errors::RingFindingError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Rings(value)=>crate::ring_errors::RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn conjugation(&self,py:Python<'_>)->PyResult<Option<Py<crate::sanitize_support_errors::ConjugationError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Conjugation(value)=>crate::sanitize_support_errors::ConjugationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn hybridization(&self,py:Python<'_>)->PyResult<Option<Py<crate::sanitize_leaf_errors::HybridizationError>>>{match self.inner.as_ref(){ck::BoundsQueryError::Hybridization(value)=>crate::sanitize_leaf_errors::HybridizationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn bounds(&self,py:Python<'_>)->PyResult<Option<Py<GraphBoundsError>>>{
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
    fn atom_index(&self)->Option<i32>{match self.inner.as_ref(){ck::GenerationError::CoordinateMapIndex {atom_index,..}|ck::GenerationError::CoordinateMapNonFinite {atom_index}=>Some(*atom_index),_=>None}}
    #[getter]
    fn message(&self)->Option<&'static str>{match self.inner.as_ref(){ck::GenerationError::Input(value)=>Some(*value),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn topology(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_errors::TopologyValidationError>>>{match self.inner.as_ref(){ck::GenerationError::Topology(value)=>crate::model_errors::TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn coordinates(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_errors::CoordinateValidationError>>>{match self.inner.as_ref(){ck::GenerationError::Coordinates(value)=>crate::model_errors::CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn valence(&self,py:Python<'_>)->PyResult<Option<Py<ValenceError>>>{match self.inner.as_ref(){ck::GenerationError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn rings(&self,py:Python<'_>)->PyResult<Option<Py<crate::ring_errors::RingFindingError>>>{match self.inner.as_ref(){ck::GenerationError::Rings(value)=>crate::ring_errors::RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn conjugation(&self,py:Python<'_>)->PyResult<Option<Py<crate::sanitize_support_errors::ConjugationError>>>{match self.inner.as_ref(){ck::GenerationError::Conjugation(value)=>crate::sanitize_support_errors::ConjugationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn hybridization(&self,py:Python<'_>)->PyResult<Option<Py<crate::sanitize_leaf_errors::HybridizationError>>>{match self.inner.as_ref(){ck::GenerationError::Hybridization(value)=>crate::sanitize_leaf_errors::HybridizationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn stereo(&self,py:Python<'_>)->PyResult<Option<Py<crate::stereo_errors::LegacyStereoError>>>{match self.inner.as_ref(){ck::GenerationError::Stereo(value)=>crate::stereo_errors::LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn fragments(&self,py:Python<'_>)->PyResult<Option<Py<crate::fragment_errors::MoleculeFragmentsError>>>{match self.inner.as_ref(){ck::GenerationError::Fragments(value)=>crate::fragment_errors::MoleculeFragmentsError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn seed(&self,py:Python<'_>)->PyResult<Option<Py<EmbedSeedError>>>{match self.inner.as_ref(){ck::GenerationError::Seed(value)=>EmbedSeedError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn atropisomer(&self,py:Python<'_>)->PyResult<Option<Py<crate::stereo_errors::AtropisomerError>>>{match self.inner.as_ref(){ck::GenerationError::Atropisomer(value)=>crate::stereo_errors::AtropisomerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn bounds(&self,py:Python<'_>)->PyResult<Option<Py<BoundsMatrixError>>>{match self.inner.as_ref(){ck::GenerationError::Bounds(value)=>BoundsMatrixError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn force_field(&self,py:Python<'_>)->PyResult<Option<Py<ConformerOptimizerError>>>{match self.inner.as_ref(){ck::GenerationError::ForceField(value)=>ConformerOptimizerError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn numeric(&self,py:Python<'_>)->PyResult<Option<Py<ConformerError>>>{match self.inner.as_ref(){ck::GenerationError::Numeric(value)=>ConformerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn hydrogens(&self,py:Python<'_>)->PyResult<Option<Py<MissingExplicitHydrogensError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::Hydrogens(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::Hydrogens(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            MissingExplicitHydrogensError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn pruning(&self,py:Python<'_>)->PyResult<Option<Py<PruningError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::Pruning(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::Pruning(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            PruningError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn generation_failed(&self,py:Python<'_>)->PyResult<Option<Py<GenerationFailure>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::GenerationFailed(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::GenerationFailed(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            GenerationFailure::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn thread_count(&self,py:Python<'_>)->PyResult<Option<Py<ThreadCountError>>>{
        if matches!(self.inner.as_ref(),ck::GenerationError::ThreadCount(_)){
            let node=self.inner.node().project(|parent|{let ck::GenerationError::ThreadCount(value)=parent.downcast_ref::<ck::GenerationError>().expect("actual generation owner")else{unreachable!()};value});
            ThreadCountError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
});

// The concrete standard poison value carries unit data and no source/kind.
#[doc = r#"{"role":"error_leaf"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct PoisonError {
    inner: RetainedError<std::sync::PoisonError<()>>,
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
                Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl PoisonError {
    #[getter]
    fn domain(&self) -> &'static str {
        "standard"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

pub(crate) type EnumerationSource =
    Source<GenerationSource, Source<AtropSource, Source<Py<PoisonError>, SmilesSource>>>;
error_projection!(retained EnumerationError,EnumerationErrorKind,ck::EnumerationError,"stereo",{
    Valence(_),Rings(_),Stereo(_),DoubleBond(_),SmilesWrite(_),PotentialStereo(_),Hydrogens(_),Embedding(_),
    AtomProperty(_),BondProperty(_),MoleculeProperty(_),Topology(_),Coordinates(_),
    InvalidFlipperAtom {atom:usize=>|value:&ck::AtomId|value.index()},InvalidFlipperBond {bond:usize=>|value:&ck::BondId|value.index()},
    RandomBitsSource {bit_count:usize,message:String},RandomSourcePoisoned(_),EmbeddingCoordinates,
},cause=|inner:&RetainedError<ck::EnumerationError>|inner.node().source().map(project),source=Option<EnumerationSource>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn valence(&self,py:Python<'_>)->PyResult<Option<Py<ValenceError>>>{match self.inner.as_ref(){ck::EnumerationError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn rings(&self,py:Python<'_>)->PyResult<Option<Py<crate::ring_errors::RingFindingError>>>{match self.inner.as_ref(){ck::EnumerationError::Rings(value)=>crate::ring_errors::RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn stereo(&self,py:Python<'_>)->PyResult<Option<Py<crate::stereo_errors::LegacyStereoError>>>{match self.inner.as_ref(){ck::EnumerationError::Stereo(value)=>crate::stereo_errors::LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn double_bond(&self,py:Python<'_>)->PyResult<Option<Py<crate::double_stereo_errors::DoubleBondStereoError>>>{match self.inner.as_ref(){ck::EnumerationError::DoubleBond(value)=>crate::double_stereo_errors::DoubleBondStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn smiles_write(&self,py:Python<'_>)->PyResult<Option<Py<crate::smiles_parse_errors::SmilesParseError>>>{match self.inner.as_ref(){ck::EnumerationError::SmilesWrite(value)=>crate::smiles_parse_errors::SmilesParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn potential_stereo(&self,py:Python<'_>)->PyResult<Option<Py<crate::potential_stereo_errors::PotentialStereoError>>>{match self.inner.as_ref(){ck::EnumerationError::PotentialStereo(value)=>crate::potential_stereo_errors::PotentialStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn hydrogens(&self,py:Python<'_>)->PyResult<Option<Py<crate::hydrogen_errors::HydrogenError>>>{match self.inner.as_ref(){ck::EnumerationError::Hydrogens(value)=>crate::hydrogen_errors::HydrogenError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn atom_property(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_edit_errors::AtomPropertyError>>>{match self.inner.as_ref(){ck::EnumerationError::AtomProperty(value)=>crate::model_edit_errors::AtomPropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn bond_property(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_edit_errors::BondValueError>>>{match self.inner.as_ref(){ck::EnumerationError::BondProperty(value)=>crate::model_edit_errors::BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn molecule_property(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_edit_errors::MoleculePropertyError>>>{match self.inner.as_ref(){ck::EnumerationError::MoleculeProperty(value)=>crate::model_edit_errors::MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn topology(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_errors::TopologyValidationError>>>{match self.inner.as_ref(){ck::EnumerationError::Topology(value)=>crate::model_errors::TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn coordinates(&self,py:Python<'_>)->PyResult<Option<Py<crate::model_errors::CoordinateValidationError>>>{match self.inner.as_ref(){ck::EnumerationError::Coordinates(value)=>crate::model_errors::CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn embedding(&self,py:Python<'_>)->PyResult<Option<Py<GenerationError>>>{
        if matches!(self.inner.as_ref(),ck::EnumerationError::Embedding(_)){
            let node=self.inner.node().project(|parent|{let ck::EnumerationError::Embedding(value)=parent.downcast_ref::<ck::EnumerationError>().expect("actual enumeration owner")else{unreachable!()};value});
            GenerationError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
});
