//! Exact construction/optimization graph; immutable views retain actual owners.
use crate::alignment_errors::ThreadCountError;
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::error_framework::{Source, error_projection};
use crate::fragment_errors::MoleculeFragmentsError;
use crate::io_source_errors::ErrorRef;
use crate::kekulize_errors::KekulizeSource;
use crate::leaf_values::UffOptimizationErrorKind;
use crate::match_errors::MatchError;
use crate::mmff_properties::{MmffMolPropertiesError, MmffParamParseError};
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
};
use crate::model_errors::{StereoGroupError, TopologyValidationError};
use crate::path_errors::{PathError, PathSource};
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::smarts_errors::SmartsParseError;
use crate::string_errors::QueryCompileError;
use crate::uff_parameter_errors::{
    ForceFieldKernelError, UffBuilderError, UffDefaultTorsionQueryError, UffParamError,
    UffParameterError, UffTorsionBondQueryError, UffTypingError,
};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

// Only actual concrete source nodes from the enums below enter this function.
// Source-free variants do not reach it; unknown errors are never downgraded.
fn project(node: ErrorRef) -> PyErr {
    macro_rules! retained {
        ($($name:ident),+ $(,)?) => { $(
            if node.error().is::<ck::$name>() { return $name::from_ref(node); }
        )+ };
    }
    retained!(
        NonbondedFragmentMappingError,
        NonbondedAssemblyError,
        ForceFieldConstructionError,
        AutomaticForceFieldConstructionError,
        SingleConformerOptimizationError,
        SerialUffOptimizationError,
        DispatchedUffOptimizationError,
        UffPreparedOptimizationError,
        UffParameterError,
        MmffBuilderError,
        MmffOptimizationFailure,
        MmffOptimizationSourceError,
        UffOptimizationError,
        UffSingleError,
        UffConformerError,
        UffEvaluationError,
        ThreadCountError
    );
    if let Some(value) = node.error().downcast_ref::<ck::MoleculeFragmentsError>() {
        return MoleculeFragmentsError::from_rust(value.clone());
    }
    macro_rules! copied {
        ($($name:ident),+ $(,)?) => { $(
            if let Some(value) = node.error().downcast_ref::<ck::$name>() {
                return $name::from_rust(value.clone());
            }
        )+ };
    }
    copied!(
        FragmentCoordinateViewError,
        UffBuilderError,
        UffParamError,
        UffTypingError,
        ForceFieldKernelError,
        OptimizationStageError,
        SerialResultCapacityError,
        SerialCoordinateCountError,
        SerialConformerOptimizationError,
        UffThreadCountError,
        MmffMolPropertiesError,
        MmffParamParseError,
        SmartsParseError,
        QueryCompileError,
        MatchError,
        UffTorsionBondQueryError,
        UffDefaultTorsionQueryError,
        PathError,
        TopologyValidationError,
        MappingValidationError,
        MoleculePropertyError,
        StereoGroupError,
        TemplateAttachmentOrderError,
        AtomPropertyError,
        BondValueError,
        RingFindingError
    );
    unreachable!("unmodeled concrete native force-field execution source")
}

error_projection!(FragmentCoordinateViewError,FragmentCoordinateViewErrorKind,ck::FragmentCoordinateViewError,"forcefield",{
    SelectedConformerMissing {id},SelectedConformerDuplicated {id},
    KernelCoordinateRowCount {id,expected:usize,actual},
    KernelCoordinateRowWidth {id,row:usize,actual},
},cause=|_:&ck::FragmentCoordinateViewError|None::<PyErr>,methods={
    #[getter]
    fn id(&self)->usize { match self.inner {
        ck::FragmentCoordinateViewError::SelectedConformerMissing {id}
        |ck::FragmentCoordinateViewError::SelectedConformerDuplicated {id}
        |ck::FragmentCoordinateViewError::KernelCoordinateRowCount {id,..}
        |ck::FragmentCoordinateViewError::KernelCoordinateRowWidth {id,..}=>id,
    }}
    #[getter]
    fn actual(&self)->Option<usize>{match self.inner {
        ck::FragmentCoordinateViewError::KernelCoordinateRowCount {actual,..}
        |ck::FragmentCoordinateViewError::KernelCoordinateRowWidth {actual,..}=>Some(actual),_=>None,
    }}
});

type FragmentSource = Source<Py<MoleculeFragmentsError>, Py<FragmentCoordinateViewError>>;
error_projection!(retained NonbondedFragmentMappingError,NonbondedFragmentMappingErrorKind,ck::NonbondedFragmentMappingError,"forcefield",{
    FragmentCopy(_),CoordinateView(_),RequestedMappingMissing,
    AtomOutOfRange {component_index:usize,atom_index,atom_count:usize},
    DuplicateAtom {atom_index,first_component:usize,second_component:usize},
    MissingAtom {atom_index},
},cause=|inner:&RetainedError<ck::NonbondedFragmentMappingError>|inner.node().source().map(project),source=Option<FragmentSource>,methods={
    #[getter]
    fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){
        ck::NonbondedFragmentMappingError::AtomOutOfRange {atom_index,..}
        |ck::NonbondedFragmentMappingError::DuplicateAtom {atom_index,..}
        |ck::NonbondedFragmentMappingError::MissingAtom {atom_index}=>Some(*atom_index),_=>None,
    }}
});
error_projection!(retained NonbondedAssemblyError,NonbondedAssemblyErrorKind,ck::NonbondedAssemblyError,"forcefield",{
    FragmentMapping(_),PairBuilder(_),
},cause=|inner:&RetainedError<ck::NonbondedAssemblyError>|inner.node().source().map(project),
source=Option<Source<Py<NonbondedFragmentMappingError>,Py<UffBuilderError>>>);
error_projection!(retained ForceFieldConstructionError,ForceFieldConstructionErrorKind,ck::ForceFieldConstructionError,"forcefield",{
    Builder(_),Nonbonded(_),
},cause=|inner:&RetainedError<ck::ForceFieldConstructionError>|inner.node().source().map(project),
source=Option<Source<Py<UffBuilderError>,Py<NonbondedAssemblyError>>>);
type AutomaticSource =
    Source<Py<UffParamError>, Source<Py<UffTypingError>, Py<ForceFieldConstructionError>>>;
error_projection!(retained AutomaticForceFieldConstructionError,AutomaticForceFieldConstructionErrorKind,ck::AutomaticForceFieldConstructionError,"forcefield",{
    ParameterTable(_),Typing(_),Construction(_),
},cause=|inner:&RetainedError<ck::AutomaticForceFieldConstructionError>|inner.node().source().map(project),source=Option<AutomaticSource>);

error_projection!(OptimizationStageError,OptimizationStageErrorKind,ck::OptimizationStageError,"forcefield",{
    Initialize(_),Minimize(_),FinalEnergy(_),
},cause=|inner:&ck::OptimizationStageError| match inner {
    ck::OptimizationStageError::Initialize(v)|ck::OptimizationStageError::Minimize(v)
    |ck::OptimizationStageError::FinalEnergy(v)=>Some(ForceFieldKernelError::from_rust(*v)),
},source=Option<Py<ForceFieldKernelError>>);

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct SerialResultCapacityError {
    inner: ck::SerialResultCapacityError,
}
impl SerialResultCapacityError {
    pub(crate) fn from_rust(inner: ck::SerialResultCapacityError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}
#[pymethods]
impl SerialResultCapacityError {
    #[getter]
    fn domain(&self) -> &'static str {
        "forcefield"
    }
    #[getter]
    fn conformers(&self) -> usize {
        self.inner.conformers()
    }
    #[getter]
    fn result_slots(&self) -> usize {
        self.inner.result_slots()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct SerialCoordinateCountError {
    inner: ck::SerialCoordinateCountError,
}
impl SerialCoordinateCountError {
    pub(crate) fn from_rust(inner: ck::SerialCoordinateCountError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}
#[pymethods]
impl SerialCoordinateCountError {
    #[getter]
    fn domain(&self) -> &'static str {
        "forcefield"
    }
    #[getter]
    fn conformer_id(&self) -> usize {
        self.inner.conformer_id()
    }
    #[getter]
    fn atoms(&self) -> usize {
        self.inner.atoms()
    }
    #[getter]
    fn coordinates(&self) -> usize {
        self.inner.coordinates()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
type SerialSource = Source<
    Py<SerialResultCapacityError>,
    Source<Py<SerialCoordinateCountError>, Py<OptimizationStageError>>,
>;
error_projection!(SerialConformerOptimizationError,SerialConformerOptimizationErrorKind,ck::SerialConformerOptimizationError,"forcefield",{
    ResultCapacity(_),CoordinateCount {input_index,source},
    Optimization {input_index,conformer_id:usize,source},
},cause=|inner:&ck::SerialConformerOptimizationError|match inner {
    ck::SerialConformerOptimizationError::ResultCapacity(v)=>Some(SerialResultCapacityError::from_rust(*v)),
    ck::SerialConformerOptimizationError::CoordinateCount {source,..}=>Some(SerialCoordinateCountError::from_rust(*source)),
    ck::SerialConformerOptimizationError::Optimization {source,..}=>Some(OptimizationStageError::from_rust(*source)),
},source=Option<SerialSource>,methods={
    #[getter]
    fn input_index(&self)->Option<usize>{match self.inner {
        ck::SerialConformerOptimizationError::CoordinateCount {input_index,..}
        |ck::SerialConformerOptimizationError::Optimization {input_index,..}=>Some(input_index),_=>None,
    }}
});
error_projection!(retained SingleConformerOptimizationError,SingleConformerOptimizationErrorKind,ck::SingleConformerOptimizationError,"uff",{
    Construction(_),Optimization(_),
},cause=|inner:&RetainedError<ck::SingleConformerOptimizationError>|inner.node().source().map(project),
source=Option<Source<Py<AutomaticForceFieldConstructionError>,Py<OptimizationStageError>>>);
error_projection!(retained SerialUffOptimizationError,SerialUffOptimizationErrorKind,ck::SerialUffOptimizationError,"uff",{
    Construction(_),Optimization(_),
},cause=|inner:&RetainedError<ck::SerialUffOptimizationError>|inner.node().source().map(project),
source=Option<Source<Py<AutomaticForceFieldConstructionError>,Py<SerialConformerOptimizationError>>>);
error_projection!(UffThreadCountError,UffThreadCountErrorKind,ck::UffThreadCountError,"uff",{
    UndefinedSignedNegation,
},cause=|_:&ck::UffThreadCountError|None::<PyErr>);
type DispatchSource = Source<
    Py<AutomaticForceFieldConstructionError>,
    Source<Py<UffThreadCountError>, Py<ThreadCountError>>,
>;
error_projection!(retained DispatchedUffOptimizationError,DispatchedUffOptimizationErrorKind,ck::DispatchedUffOptimizationError,"uff",{
    Construction(_),ThreadCount(_),ThreadObservation(_),
},cause=|inner:&RetainedError<ck::DispatchedUffOptimizationError>|inner.node().source().map(project),source=Option<DispatchSource>);
type PreparedSource = Source<
    Py<UffParameterError>,
    Source<
        Py<SingleConformerOptimizationError>,
        Source<
            Py<SerialUffOptimizationError>,
            Source<Py<SerialConformerOptimizationError>, Py<DispatchedUffOptimizationError>>,
        >,
    >,
>;
error_projection!(retained UffPreparedOptimizationError,UffPreparedOptimizationErrorKind,ck::UffPreparedOptimizationError,"uff",{
    Preparation(_),Single(_),Serial(_),ConformerStage(_),Dispatch(_),
},cause=|inner:&RetainedError<ck::UffPreparedOptimizationError>|inner.node().source().map(project),source=Option<PreparedSource>);

// Native records without a public discriminator retain their real concrete
// source and Display. No artificial singleton variant is introduced.
macro_rules! source_record {
    ($name:ident,$native:ident,$source:ty) => {
        #[doc = r#"{"role":"error_record"}"#]
        #[pyclass(module="cosmolkit",extends=PyValueError)]
        pub struct $name {
            inner: RetainedError<ck::$native>,
        }
        impl $name {
            pub(crate) fn from_rust(inner: ck::$native) -> PyErr {
                Self::from_retained(RetainedError::new(inner))
            }
            pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
                Self::from_retained(RetainedError::from_ref(node))
            }
            fn from_retained(inner: RetainedError<ck::$native>) -> PyErr {
                let cause = inner.node().source().map(project);
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
        impl $name {
            #[getter]
            fn domain(&self) -> &'static str {
                "forcefield"
            }
            #[getter]
            #[doc = r#"{"role":"error_source"}"#]
            fn source(slf: Bound<'_, Self>) -> PyResult<$source> {
                Ok(slf.getattr("__cause__")?.extract()?)
            }
            fn __str__(&self) -> String {
                self.inner.to_string()
            }
        }
    };
}
source_record!(
    UffSingleError,
    UffSingleError,
    Py<UffPreparedOptimizationError>
);
source_record!(
    UffConformerError,
    UffConformerError,
    Py<UffPreparedOptimizationError>
);
source_record!(
    MmffOptimizationSourceError,
    MmffOptimizationSourceError,
    Py<MmffOptimizationFailure>
);
source_record!(
    MmffOptimizationError,
    MmffOptimizationError,
    Py<MmffOptimizationSourceError>
);

type TorsionSource = Source<
    Py<UffDefaultTorsionQueryError>,
    Source<Py<SmartsParseError>, Source<Py<QueryCompileError>, Py<MatchError>>>,
>;
type MmffBuilderSource = Source<
    Source<KekulizeSource, Py<MmffParamParseError>>,
    Source<TorsionSource, Source<FragmentSource, PathSource>>,
>;
error_projection!(retained MmffBuilderError,MmffBuilderErrorKind,ck::MmffBuilderError,"mmff",{
    MolProperties(_),ForceFieldKernel(_),TorsionBondQuery(_),Fragment(_),
    Components(_),Valence(_),Rings(_),Missing3dConformer {conf_id},InvalidConformerId {conf_id},
},cause=|inner:&RetainedError<ck::MmffBuilderError>|inner.node().source().map(project),source=Option<MmffBuilderSource>,methods={
    #[getter]
    fn conf_id(&self)->Option<isize>{match self.inner.as_ref(){
        ck::MmffBuilderError::Missing3dConformer {conf_id}|ck::MmffBuilderError::InvalidConformerId {conf_id}=>Some(*conf_id),_=>None,
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn mol_properties(&self,py:Python<'_>)->PyResult<Option<Py<MmffMolPropertiesError>>>{match self.inner.as_ref(){
        ck::MmffBuilderError::MolProperties(value)=>MmffMolPropertiesError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn force_field_kernel(&self,py:Python<'_>)->PyResult<Option<Py<ForceFieldKernelError>>>{match self.inner.as_ref(){
        ck::MmffBuilderError::ForceFieldKernel(value)=>ForceFieldKernelError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn torsion_bond_query(&self,py:Python<'_>)->PyResult<Option<Py<UffTorsionBondQueryError>>>{match self.inner.as_ref(){
        ck::MmffBuilderError::TorsionBondQuery(value)=>UffTorsionBondQueryError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn fragment(&self,py:Python<'_>)->PyResult<Option<Py<NonbondedFragmentMappingError>>>{
        if matches!(self.inner.as_ref(),ck::MmffBuilderError::Fragment(_)) {
            let node=self.inner.node().project(|parent|{
                let ck::MmffBuilderError::Fragment(value)=parent.downcast_ref::<ck::MmffBuilderError>().expect("native MMFF builder payload owner") else {unreachable!("native fragment payload branch")};
                value
            });
            NonbondedFragmentMappingError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn components(&self,py:Python<'_>)->PyResult<Option<Py<PathError>>>{match self.inner.as_ref(){
        ck::MmffBuilderError::Components(value)=>PathError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn valence(&self,py:Python<'_>)->PyResult<Option<Py<ValenceError>>>{match self.inner.as_ref(){
        ck::MmffBuilderError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn rings(&self,py:Python<'_>)->PyResult<Option<Py<RingFindingError>>>{match self.inner.as_ref(){
        ck::MmffBuilderError::Rings(value)=>RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
});
type MmffFailureSource = Source<MmffBuilderSource, Source<Py<ForceFieldKernelError>, SerialSource>>;
error_projection!(retained MmffOptimizationFailure,MmffOptimizationFailureKind,ck::MmffOptimizationFailure,"mmff",{
    Properties(_),Builder(_),Kernel(_),Stage(_),Serial(_),Threads(_),
},cause=|inner:&RetainedError<ck::MmffOptimizationFailure>|inner.node().source().map(project),source=Option<MmffFailureSource>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn properties(&self,py:Python<'_>)->PyResult<Option<Py<MmffMolPropertiesError>>>{match self.inner.as_ref(){
        ck::MmffOptimizationFailure::Properties(value)=>MmffMolPropertiesError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn builder(&self,py:Python<'_>)->PyResult<Option<Py<MmffBuilderError>>>{
        if matches!(self.inner.as_ref(),ck::MmffOptimizationFailure::Builder(_)) {
            let node=self.inner.node().project(|parent|{
                let ck::MmffOptimizationFailure::Builder(value)=parent.downcast_ref::<ck::MmffOptimizationFailure>().expect("native MMFF failure payload owner") else {unreachable!("native builder payload branch")};
                value
            });
            MmffBuilderError::from_ref(node).value(py).extract().map(Some).map_err(Into::into)
        }else{Ok(None)}
    }
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn kernel(&self,py:Python<'_>)->PyResult<Option<Py<ForceFieldKernelError>>>{match self.inner.as_ref(){
        ck::MmffOptimizationFailure::Kernel(value)=>ForceFieldKernelError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn stage(&self,py:Python<'_>)->PyResult<Option<Py<OptimizationStageError>>>{match self.inner.as_ref(){
        ck::MmffOptimizationFailure::Stage(value)=>OptimizationStageError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn serial(&self,py:Python<'_>)->PyResult<Option<Py<SerialConformerOptimizationError>>>{match self.inner.as_ref(){
        ck::MmffOptimizationFailure::Serial(value)=>SerialConformerOptimizationError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn threads(&self,py:Python<'_>)->PyResult<Option<Py<UffThreadCountError>>>{match self.inner.as_ref(){
        ck::MmffOptimizationFailure::Threads(value)=>UffThreadCountError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None),
    }}
});

crate::enums::binding_enum!(UffEvaluationFailureKind, {
    MissingConformer,
    Preparation,
    Construction,
    Kernel,
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct UffEvaluationFailure {
    inner: RetainedError<ck::UffEvaluationError>,
}
#[pymethods]
impl UffEvaluationFailure {
    #[getter]
    fn variant(&self) -> UffEvaluationFailureKind {
        match self.inner.failure() {
            ck::UffEvaluationFailure::MissingConformer { requested: _ } => {
                UffEvaluationFailureKind::MissingConformer
            }
            ck::UffEvaluationFailure::Preparation(_) => UffEvaluationFailureKind::Preparation,
            ck::UffEvaluationFailure::Construction(_) => UffEvaluationFailureKind::Construction,
            ck::UffEvaluationFailure::Kernel(_) => UffEvaluationFailureKind::Kernel,
        }
    }
    #[getter]
    fn requested(&self) -> Option<usize> {
        match self.inner.failure() {
            ck::UffEvaluationFailure::MissingConformer { requested } => *requested,
            _ => None,
        }
    }
    #[getter]
    fn preparation(&self, py: Python<'_>) -> PyResult<Option<Py<UffParameterError>>> {
        if matches!(
            self.inner.failure(),
            ck::UffEvaluationFailure::Preparation(_)
        ) {
            project(
                self.inner
                    .node()
                    .source()
                    .expect("actual native preparation source"),
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
    fn construction(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<Py<AutomaticForceFieldConstructionError>>> {
        if matches!(
            self.inner.failure(),
            ck::UffEvaluationFailure::Construction(_)
        ) {
            project(
                self.inner
                    .node()
                    .source()
                    .expect("actual native construction source"),
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
    fn kernel(&self, py: Python<'_>) -> PyResult<Option<Py<ForceFieldKernelError>>> {
        if matches!(self.inner.failure(), ck::UffEvaluationFailure::Kernel(_)) {
            project(
                self.inner
                    .node()
                    .source()
                    .expect("actual native kernel source"),
            )
            .value(py)
            .extract()
            .map(Some)
            .map_err(Into::into)
        } else {
            Ok(None)
        }
    }
}
type EvaluationSource = Source<
    Py<UffParameterError>,
    Source<Py<AutomaticForceFieldConstructionError>, Py<ForceFieldKernelError>>,
>;
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct UffEvaluationError {
    inner: RetainedError<ck::UffEvaluationError>,
}
impl UffEvaluationError {
    pub(crate) fn from_rust(inner: ck::UffEvaluationError) -> PyErr {
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
impl UffEvaluationError {
    #[getter]
    fn domain(&self) -> &'static str {
        "uff"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn failure(&self) -> UffEvaluationFailure {
        UffEvaluationFailure {
            inner: self.inner.clone(),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<EvaluationSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
type UffSource = Source<
    Py<RingFindingError>,
    Source<Py<UffSingleError>, Source<Py<UffConformerError>, Py<UffEvaluationError>>>,
>;
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct UffOptimizationError {
    inner: RetainedError<ck::UffOptimizationError>,
}
impl UffOptimizationError {
    pub(crate) fn from_rust(inner: ck::UffOptimizationError) -> PyErr {
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
impl UffOptimizationError {
    #[getter]
    fn domain(&self) -> &'static str {
        "uff"
    }
    #[getter]
    fn kind(&self) -> UffOptimizationErrorKind {
        UffOptimizationErrorKind {
            inner: self.inner.kind(),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<UffSource>> {
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
    fn execution_graph_preserves_actual_native_failures_and_transparent_payload_owners() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let before = molecule.to_binary().unwrap();
            let native = match molecule.uff_energy_gradient().unwrap_err() {
                ck::OperationError::UffOptimization(error) => error,
                error => panic!("unexpected native evaluation error: {error:?}"),
            };
            let message = native.to_string();
            let error = UffOptimizationError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert!(error.value(py).getattr("source")?.is_none());
            assert_eq!(
                error
                    .value(py)
                    .call_method0("kind")?
                    .extract::<PyRef<'_, UffOptimizationErrorKind>>()?
                    .inner,
                ck::UffOptimizationErrorKind::MissingConformer { requested: None }
            );
            assert_eq!(molecule.to_binary().unwrap(), before);

            let native = ck::MmffOptimizationFailure::Builder(ck::MmffBuilderError::Fragment(
                ck::NonbondedFragmentMappingError::MissingAtom {
                    atom_index: usize::MAX,
                },
            ));
            assert!(native.source().is_none());
            let error = MmffOptimizationFailure::from_rust(native);
            let parent = error.value(py);
            assert!(parent.getattr("__cause__")?.is_none());
            let builder = parent
                .getattr("builder")?
                .extract::<Py<MmffBuilderError>>()?;
            let payload = builder
                .bind(py)
                .getattr("fragment")?
                .extract::<Py<NonbondedFragmentMappingError>>()?;
            {
                let parent = parent.extract::<PyRef<'_, MmffOptimizationFailure>>()?;
                let ck::MmffOptimizationFailure::Builder(ck::MmffBuilderError::Fragment(original)) =
                    parent.inner.as_ref()
                else {
                    unreachable!()
                };
                let child = payload.borrow(py);
                assert!(std::ptr::eq(original, child.inner.as_ref()));
            }
            drop(error);
            drop(builder);
            assert_eq!(
                payload.bind(py).getattr("atom_index")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(payload.bind(py).getattr("source")?.is_none());
            assert!(payload.bind(py).setattr("atom_index", 0).is_err());

            let io = std::fs::File::open("/proc/self/this-is-not-a-cosmolkit-file").unwrap_err();
            let io_message = io.to_string();
            let raw = io.raw_os_error();
            let native = ck::UffPreparedOptimizationError::Dispatch(
                ck::DispatchedUffOptimizationError::ThreadObservation(
                    ck::ThreadCountError::ObservationIo(io),
                ),
            );
            let error = UffPreparedOptimizationError::from_rust(native);
            let dispatch = error.value(py).getattr("source")?;
            assert!(dispatch.is(&error.value(py).getattr("__cause__")?));
            let observation = dispatch.getattr("source")?;
            assert!(observation.is(&dispatch.getattr("__cause__")?));
            let leaf = observation
                .getattr("source")?
                .extract::<Py<crate::io_source_errors::IoError>>()?;
            drop(error);
            assert_eq!(leaf.bind(py).str()?.to_str()?, io_message);
            assert_eq!(
                leaf.bind(py)
                    .getattr("raw_os_error")?
                    .extract::<Option<i32>>()?,
                raw
            );

            let native = ck::MmffOptimizationFailure::Builder(ck::MmffBuilderError::Fragment(
                ck::NonbondedFragmentMappingError::CoordinateView(
                    ck::FragmentCoordinateViewError::KernelCoordinateRowWidth {
                        id: usize::MAX,
                        row: 7,
                        actual: 2,
                    },
                ),
            ));
            let error = MmffOptimizationFailure::from_rust(native);
            let source = error.value(py).getattr("source")?;
            assert!(source.is_instance_of::<FragmentCoordinateViewError>());
            let builder = error.value(py).getattr("builder")?;
            let fragment = builder.getattr("fragment")?;
            let separate_source = fragment.getattr("source")?;
            assert!(separate_source.is(&fragment.getattr("__cause__")?));
            assert_eq!(source.getattr("id")?.extract::<usize>()?, usize::MAX);
            assert_eq!(separate_source.getattr("row")?.extract::<usize>()?, 7);
            assert_eq!(separate_source.getattr("actual")?.extract::<usize>()?, 2);
            assert!(separate_source.setattr("actual", 3).is_err());

            let native = ck::SerialConformerOptimizationError::Optimization {
                input_index: usize::MAX,
                conformer_id: 71,
                source: ck::OptimizationStageError::Minimize(ck::ForceFieldKernelError::NoPoints),
            };
            let error = SerialConformerOptimizationError::from_rust(native);
            assert_eq!(
                error.value(py).getattr("input_index")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("conformer_id")?
                    .extract::<usize>()?,
                71
            );
            let stage = error.value(py).getattr("source")?;
            assert!(stage.is(&error.value(py).getattr("__cause__")?));
            assert!(
                stage
                    .getattr("source")?
                    .is_instance_of::<ForceFieldKernelError>()
            );
            Ok(())
        })
        .unwrap();
    }
}
