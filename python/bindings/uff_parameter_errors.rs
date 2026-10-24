//! Full native UFF parameter/builder failure values; no force-field state is exposed.
use crate::canonical_errors::PropertyStringError;
use crate::chemical_values::{Hybridization, enum_value};
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::ErrorRef;
use crate::leaf_values::UffParameterErrorKind;
use crate::match_errors::MatchError;
use crate::model_errors::TopologyValidationError;
use crate::retained_error::RetainedError;
use crate::smarts_errors::SmartsParseError;
use crate::string_errors::QueryCompileError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

enum_value!(ForceFieldIndexArgument,{I,J});
enum_value!(ForceFieldBondIndexArgument,{First,Second});
enum_value!(ForceFieldAngleIndexArgument,{First,Second,Third});
enum_value!(ForceFieldTorsionIndexArgument,{First,Second,Third,Fourth});
enum_value!(ForceFieldAngleRangeBound,{Minimum,Maximum});
enum_value!(UffInversionIndexArgument,{First,Second,Third,Fourth});
enum_value!(UffPreparedValenceField,{Explicit,ImplicitHydrogen});
enum_value!(UffTypingInput,{TotalValence,ConjugatedBondPresence});

type KernelArgument = Source<
    ForceFieldIndexArgument,
    Source<
        ForceFieldBondIndexArgument,
        Source<ForceFieldAngleIndexArgument, ForceFieldTorsionIndexArgument>,
    >,
>;

error_projection!(ForceFieldKernelError,ForceFieldKernelErrorKind,ck::ForceFieldKernelError,"forcefield",{
    NoPoints,NoDistanceMatrix,MatrixSizeMismatch,NotInitialized,BadBounds,BadBondOrder,
    PositionCountMismatch { num_points:u32,position_count:usize }, IndexOutOfRange {argument,index,upper_bound},
    BondIndexOutOfRange {argument,index,upper_bound},AngleDegeneratePoints,AngleIndexOutOfRange {argument,index,upper_bound},
    AngleBadOrder {order},TorsionIndexOutOfRange {argument,index,upper_bound},TorsionDegeneratePoints,TorsionBadHybridizations,
    TorsionBadOrder {order},AngleOutOfRange {bound},AngleBoundsOrder,PackedAngleBoundsOrder,TorsionBoundsOrder,
    BadIndex,OptimizerBadTolerance,OptimizerBadDirection,BadFixedPoint {index,upper_bound},TransferPostcondition,
    OopParametersMissing,TorsionParametersOutsideDefinedSource,TorsionAddressOutsideDefinedSource {index,buffer_len},
    OopAddressOutsideDefinedSource {index,buffer_len},
},cause=|_:&ck::ForceFieldKernelError|None::<PyErr>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn argument(&self)->Option<KernelArgument>{match self.inner {
        ck::ForceFieldKernelError::IndexOutOfRange {argument,..}=>Some(Source::First(argument.into())),
        ck::ForceFieldKernelError::BondIndexOutOfRange {argument,..}=>Some(Source::Second(Source::First(argument.into()))),
        ck::ForceFieldKernelError::AngleIndexOutOfRange {argument,..}=>Some(Source::Second(Source::Second(Source::First(argument.into())))),
        ck::ForceFieldKernelError::TorsionIndexOutOfRange {argument,..}=>Some(Source::Second(Source::Second(Source::Second(argument.into())))),_=>None,
    }}
    #[getter]
    fn index(&self)->Option<i64>{match self.inner {
        ck::ForceFieldKernelError::IndexOutOfRange {index,..}|ck::ForceFieldKernelError::BondIndexOutOfRange {index,..}
        |ck::ForceFieldKernelError::AngleIndexOutOfRange {index,..}|ck::ForceFieldKernelError::TorsionIndexOutOfRange {index,..}=>Some(i64::from(index)),
        ck::ForceFieldKernelError::BadFixedPoint {index,..}|ck::ForceFieldKernelError::OopAddressOutsideDefinedSource {index,..}=>Some(i64::from(index)),
        ck::ForceFieldKernelError::TorsionAddressOutsideDefinedSource {index,..}=>Some(i64::from(index)),_=>None,
    }}
    #[getter]
    fn upper_bound(&self)->Option<usize>{match self.inner {
        ck::ForceFieldKernelError::IndexOutOfRange {upper_bound,..}|ck::ForceFieldKernelError::BadFixedPoint {upper_bound,..}=>Some(upper_bound as usize),
        ck::ForceFieldKernelError::BondIndexOutOfRange {upper_bound,..}|ck::ForceFieldKernelError::AngleIndexOutOfRange {upper_bound,..}
        |ck::ForceFieldKernelError::TorsionIndexOutOfRange {upper_bound,..}=>Some(upper_bound),_=>None,
    }}
    #[getter]
    fn order(&self)->Option<u32>{match self.inner {ck::ForceFieldKernelError::AngleBadOrder {order}|ck::ForceFieldKernelError::TorsionBadOrder {order}=>Some(order),_=>None}}
    #[getter]
    fn buffer_len(&self)->Option<usize>{match self.inner {ck::ForceFieldKernelError::TorsionAddressOutsideDefinedSource {buffer_len,..}|ck::ForceFieldKernelError::OopAddressOutsideDefinedSource {buffer_len,..}=>Some(buffer_len),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn bound(&self)->Option<ForceFieldAngleRangeBound>{match self.inner {ck::ForceFieldKernelError::AngleOutOfRange {bound}=>Some(bound.into()),_=>None}}
});

error_projection!(UffInversionContributionError,UffInversionContributionErrorKind,ck::UffInversionContributionError,"uff",{
    IndexOutOfRange { argument,index:u32,upper_bound:usize },
},cause=|_:&ck::UffInversionContributionError|None::<PyErr>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn argument(&self)->UffInversionIndexArgument { let ck::UffInversionContributionError::IndexOutOfRange {argument,..}=self.inner;argument.into() }
});
error_projection!(UffBondMathError,UffBondMathErrorKind,ck::UffBondMathError,"uff",{
    InvalidBondOrder {bond_order:f64},
},cause=|_:&ck::UffBondMathError|None::<PyErr>);
error_projection!(UffParamError,UffParamErrorKind,ck::UffParamError,"uff",{
    EmptyLine {line_number},MalformedLine {line_number,column_count:usize},ParseFloat {line_number,column_name:&'static str,value:String},
},cause=|_:&ck::UffParamError|None::<PyErr>,methods={
    #[getter]
    fn line_number(&self)->usize{match self.inner {ck::UffParamError::EmptyLine {line_number}|ck::UffParamError::MalformedLine {line_number,..}|ck::UffParamError::ParseFloat {line_number,..}=>line_number}}
});

type TypingSource = Source<
    crate::error_framework::ErrorObject<ValenceError>,
    Source<
        crate::error_framework::ErrorObject<PropertyStringError>,
        crate::error_framework::ErrorObject<UffBondMathError>,
    >,
>;
error_projection!(UffTypingError,UffTypingErrorKind,ck::UffTypingError,"uff",{
    CoreValence(_),CorePropertyString(_),BondMath(_),PreparedStateLength {input,expected:usize,actual:usize},
    AtomIdOutOfBounds {atom_id:usize=>|v:&ck::AtomId|v.index(),atom_count},AtomIndexOutOfBounds {index:usize,atom_count},
    BondIdOutOfBounds {bond_id:usize=>|v:&ck::BondId|v.index(),bond_count:usize},CenterBondMissingAfterSourceSuccess {idx2:usize,idx3:usize},
    CentralParameterSlotMissingAfterSourceSuccess {atom_index:usize},
},cause=|inner:&ck::UffTypingError|match inner {
    ck::UffTypingError::CoreValence(v)=>Some(ValenceError::from_rust(v.clone())),
    ck::UffTypingError::CorePropertyString(v)=>Some(PropertyStringError::from_rust(*v)),
    ck::UffTypingError::BondMath(v)=>Some(UffBondMathError::from_rust(*v)),
    ck::UffTypingError::PreparedStateLength {..}|ck::UffTypingError::AtomIdOutOfBounds {..}|ck::UffTypingError::AtomIndexOutOfBounds {..}
    |ck::UffTypingError::BondIdOutOfBounds {..}|ck::UffTypingError::CenterBondMissingAfterSourceSuccess {..}|ck::UffTypingError::CentralParameterSlotMissingAfterSourceSuccess {..}=>None,
},source=Option<TypingSource>,methods={
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn input(&self)->Option<UffTypingInput>{match self.inner {ck::UffTypingError::PreparedStateLength {input,..}=>Some(input.into()),_=>None}}
    #[getter]
    fn atom_count(&self)->Option<usize>{match self.inner {ck::UffTypingError::AtomIdOutOfBounds {atom_count,..}|ck::UffTypingError::AtomIndexOutOfBounds {atom_count,..}=>Some(atom_count),_=>None}}
});

error_projection!(UffDefaultTorsionQueryError,UffDefaultTorsionQueryErrorKind,ck::UffDefaultTorsionQueryError,"uff",{
    Parse(_),Compile(_),
},cause=|inner:&ck::UffDefaultTorsionQueryError|match inner {
    ck::UffDefaultTorsionQueryError::Parse(v)=>Some(SmartsParseError::from_rust(v.clone())),
    ck::UffDefaultTorsionQueryError::Compile(v)=>Some(QueryCompileError::from_rust(v.clone())),
},source=Option<Source<crate::error_framework::ErrorObject<SmartsParseError>,crate::error_framework::ErrorObject<QueryCompileError>>>);

type TorsionQuerySource = Source<
    crate::error_framework::ErrorObject<UffDefaultTorsionQueryError>,
    Source<
        crate::error_framework::ErrorObject<SmartsParseError>,
        Source<
            crate::error_framework::ErrorObject<QueryCompileError>,
            crate::error_framework::ErrorObject<MatchError>,
        >,
    >,
>;
error_projection!(UffTorsionBondQueryError,UffTorsionBondQueryErrorKind,ck::UffTorsionBondQueryError,"uff",{
    DefaultQuery(_),Parse(_),Compile(_),Match(_),MatchArity {smarts:String,actual:usize},MissingMatchedBond {begin_atom_index:usize,end_atom_index:usize},
},cause=|inner:&ck::UffTorsionBondQueryError|match inner {
    ck::UffTorsionBondQueryError::DefaultQuery(v)=>Some(UffDefaultTorsionQueryError::from_rust(v.clone())),
    ck::UffTorsionBondQueryError::Parse(v)=>Some(SmartsParseError::from_rust(v.clone())),
    ck::UffTorsionBondQueryError::Compile(v)=>Some(QueryCompileError::from_rust(v.clone())),
    ck::UffTorsionBondQueryError::Match(v)=>Some(MatchError::from_rust(v.clone())),
    ck::UffTorsionBondQueryError::MatchArity {..}|ck::UffTorsionBondQueryError::MissingMatchedBond {..}=>None,
},source=Option<TorsionQuerySource>);

type BuilderSource = Source<
    crate::error_framework::ErrorObject<UffTorsionBondQueryError>,
    Source<
        crate::error_framework::ErrorObject<UffInversionContributionError>,
        Source<
            crate::error_framework::ErrorObject<ValenceError>,
            Source<
                crate::error_framework::ErrorObject<ForceFieldKernelError>,
                crate::error_framework::ErrorObject<TopologyValidationError>,
            >,
        >,
    >,
>;
error_projection!(UffBuilderError,UffBuilderErrorKind,ck::UffBuilderError,"uff",{
    NeighborMatrixIndexOutOfRange {n_atoms,i:usize,j:usize},NeighborMatrixIndexOverflow {n_atoms},
    NeighborMatrixAllocationFailed {byte_count:usize},NeighborMatrixCellPositionOverflow {position},NeighborMatrixStorageOutOfRange {position,byte_index:usize,storage_len:usize},
    ValenceAssignmentLengthMismatch {field,expected:usize,actual},SourceValencePrecondition {atom_id,field,value},SourceValenceOutOfRange {atom_id,field,value},
    ParamsLengthMismatch {atoms,params:usize},SelectedThreeDimensionalConformerNotFound {conformer_id},SelectedConformerCoordinateCountMismatch {conformer_id,atoms,coordinates:usize},
    SourceBondIndexOverflow {atom_index},SourceTorsionIndexOverflow {atom_index},SourceAngleIndexOverflow {atom_index},SourceInversionIndexOverflow {atom_index},
    SourceDirectionVectorBelowTolerance {center_atom_index,neighbor_atom_index:usize},SourceTbpAtomPrecondition {center_atom_index},SourceTbpHybridizationPrecondition {center_atom_index,actual},
    SourceTbpDegreePrecondition {center_atom_index,actual_degree:usize},SourceTbpAxialBondNotFound {center_atom_index},SourceTbpEquatorialBondNotFound {center_atom_index,role:u8},
    SourceTbpCenterParamsMissing {center_atom_index},TorsionBondQuery(_),InversionContribution(_),Valence(_),ForceFieldKernel(_),TopologyValidation(_),
},cause=|inner:&ck::UffBuilderError|match inner {
    ck::UffBuilderError::TorsionBondQuery(v)=>Some(UffTorsionBondQueryError::from_rust(v.clone())),
    ck::UffBuilderError::InversionContribution(v)=>Some(UffInversionContributionError::from_rust(*v)),
    ck::UffBuilderError::Valence(v)=>Some(ValenceError::from_rust(v.clone())),
    ck::UffBuilderError::ForceFieldKernel(v)=>Some(ForceFieldKernelError::from_rust(*v)),
    ck::UffBuilderError::TopologyValidation(v)=>Some(TopologyValidationError::from_rust(v.clone())),
    ck::UffBuilderError::NeighborMatrixIndexOutOfRange {..}|ck::UffBuilderError::NeighborMatrixIndexOverflow {..}|ck::UffBuilderError::NeighborMatrixAllocationFailed {..}
    |ck::UffBuilderError::NeighborMatrixCellPositionOverflow {..}|ck::UffBuilderError::NeighborMatrixStorageOutOfRange {..}|ck::UffBuilderError::ValenceAssignmentLengthMismatch {..}
    |ck::UffBuilderError::SourceValencePrecondition {..}|ck::UffBuilderError::SourceValenceOutOfRange {..}|ck::UffBuilderError::ParamsLengthMismatch {..}
    |ck::UffBuilderError::SelectedThreeDimensionalConformerNotFound {..}|ck::UffBuilderError::SelectedConformerCoordinateCountMismatch {..}
    |ck::UffBuilderError::SourceBondIndexOverflow {..}|ck::UffBuilderError::SourceTorsionIndexOverflow {..}|ck::UffBuilderError::SourceAngleIndexOverflow {..}
    |ck::UffBuilderError::SourceInversionIndexOverflow {..}|ck::UffBuilderError::SourceDirectionVectorBelowTolerance {..}|ck::UffBuilderError::SourceTbpAtomPrecondition {..}
    |ck::UffBuilderError::SourceTbpHybridizationPrecondition {..}|ck::UffBuilderError::SourceTbpDegreePrecondition {..}|ck::UffBuilderError::SourceTbpAxialBondNotFound {..}
    |ck::UffBuilderError::SourceTbpEquatorialBondNotFound {..}|ck::UffBuilderError::SourceTbpCenterParamsMissing {..}=>None,
},source=Option<BuilderSource>,methods={
    #[getter]
    fn n_atoms(&self)->Option<usize>{match self.inner {ck::UffBuilderError::NeighborMatrixIndexOutOfRange {n_atoms,..}|ck::UffBuilderError::NeighborMatrixIndexOverflow {n_atoms}=>Some(n_atoms),_=>None}}
    #[getter]
    fn position(&self)->Option<usize>{match self.inner {ck::UffBuilderError::NeighborMatrixCellPositionOverflow {position}|ck::UffBuilderError::NeighborMatrixStorageOutOfRange {position,..}=>Some(position),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn field(&self)->Option<UffPreparedValenceField>{match self.inner {ck::UffBuilderError::ValenceAssignmentLengthMismatch {field,..}|ck::UffBuilderError::SourceValencePrecondition {field,..}|ck::UffBuilderError::SourceValenceOutOfRange {field,..}=>Some(field.into()),_=>None}}
    #[getter]
    fn atom_id(&self)->Option<usize>{match self.inner {ck::UffBuilderError::SourceValencePrecondition {atom_id,..}|ck::UffBuilderError::SourceValenceOutOfRange {atom_id,..}=>Some(atom_id.index()),_=>None}}
    #[getter]
    fn value(&self)->Option<i32>{match self.inner {ck::UffBuilderError::SourceValencePrecondition {value,..}|ck::UffBuilderError::SourceValenceOutOfRange {value,..}=>Some(value),_=>None}}
    #[getter]
    #[doc=r#"{"role":"error_payload"}"#]
    fn actual(&self)->Option<Source<usize,Hybridization>>{match self.inner {ck::UffBuilderError::ValenceAssignmentLengthMismatch {actual,..}=>Some(Source::First(actual)),ck::UffBuilderError::SourceTbpHybridizationPrecondition {actual,..}=>Some(Source::Second(actual.into())),_=>None}}
    #[getter]
    fn atoms(&self)->Option<usize>{match self.inner {ck::UffBuilderError::ParamsLengthMismatch {atoms,..}|ck::UffBuilderError::SelectedConformerCoordinateCountMismatch {atoms,..}=>Some(atoms),_=>None}}
    #[getter]
    fn conformer_id(&self)->Option<usize>{match self.inner {ck::UffBuilderError::SelectedThreeDimensionalConformerNotFound {conformer_id}|ck::UffBuilderError::SelectedConformerCoordinateCountMismatch {conformer_id,..}=>Some(conformer_id),_=>None}}
    #[getter]
    fn atom_index(&self)->Option<usize>{match self.inner {ck::UffBuilderError::SourceBondIndexOverflow {atom_index}|ck::UffBuilderError::SourceTorsionIndexOverflow {atom_index}|ck::UffBuilderError::SourceAngleIndexOverflow {atom_index}|ck::UffBuilderError::SourceInversionIndexOverflow {atom_index}=>Some(atom_index),_=>None}}
    #[getter]
    fn center_atom_index(&self)->Option<usize>{match self.inner {ck::UffBuilderError::SourceDirectionVectorBelowTolerance {center_atom_index,..}|ck::UffBuilderError::SourceTbpAtomPrecondition {center_atom_index}|ck::UffBuilderError::SourceTbpHybridizationPrecondition {center_atom_index,..}|ck::UffBuilderError::SourceTbpDegreePrecondition {center_atom_index,..}|ck::UffBuilderError::SourceTbpAxialBondNotFound {center_atom_index}|ck::UffBuilderError::SourceTbpEquatorialBondNotFound {center_atom_index,..}|ck::UffBuilderError::SourceTbpCenterParamsMissing {center_atom_index}=>Some(center_atom_index),_=>None}}
});

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct UffParameterError {
    inner: RetainedError<ck::UffParameterError>,
}
}
impl UffParameterError {
    pub(crate) fn from_rust(inner: ck::UffParameterError) -> PyErr {
        Self::from_retained(RetainedError::new(inner))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Self::from_retained(RetainedError::from_ref(node))
    }
    fn from_retained(inner: RetainedError<ck::UffParameterError>) -> PyErr {
        // The native public record always returns precisely one of these three
        // concrete source types. Retain the actual non-Clone root and select
        // by its native public kind, not its Display/Debug text.
        let source = inner
            .source()
            .expect("UffParameterError always has its concrete source");
        let cause = match inner.kind() {
            ck::UffParameterErrorKind::Preparation => UffBuilderError::from_rust(
                source
                    .downcast_ref::<ck::UffBuilderError>()
                    .expect("native preparation source")
                    .clone(),
            ),
            ck::UffParameterErrorKind::ParameterTable => UffParamError::from_rust(
                source
                    .downcast_ref::<ck::UffParamError>()
                    .expect("native parameter-table source")
                    .clone(),
            ),
            ck::UffParameterErrorKind::Typing => UffTypingError::from_rust(
                source
                    .downcast_ref::<ck::UffTypingError>()
                    .expect("native typing source")
                    .clone(),
            ),
        };
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => {
                let error = crate::error_framework::from_payload(value);
                error.set_cause(py, Some(cause));
                error
            }
            Err(error) => error,
        })
    }
}
crate::error_framework::error_projection! {methods impl UffParameterError {
    #[getter]
    fn domain(&self) -> &'static str {
        "uff"
    }
    #[getter]
    fn kind(&self) -> UffParameterErrorKind {
        self.inner.kind().into()
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(
        slf: Bound<'_, pyo3::exceptions::PyBaseException>,
    ) -> PyResult<Source<crate::error_framework::ErrorObject<UffBuilderError>, Source<crate::error_framework::ErrorObject<UffParamError>, crate::error_framework::ErrorObject<UffTypingError>>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_uff_parameter_graph_keeps_native_sources_arguments_and_integer_widths() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            assert!(
                ck::Molecule::from_smiles("CCO")
                    .unwrap()
                    .uff_has_all_molecule_params()
                    .unwrap()
            );
            let parse = ck::parse_smarts("[").unwrap_err();
            let native =
                ck::UffBuilderError::TorsionBondQuery(ck::UffTorsionBondQueryError::DefaultQuery(
                    ck::UffDefaultTorsionQueryError::Parse(parse),
                ));
            assert!(
                native
                    .source()
                    .unwrap()
                    .is::<ck::UffTorsionBondQueryError>()
            );
            let message = native.to_string();
            let error = UffBuilderError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let torsion = error.value(py).getattr("source")?;
            assert!(torsion.is(&error.value(py).getattr("__cause__")?));
            assert!(torsion.is_instance(<UffTorsionBondQueryError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let default = torsion.getattr("source")?;
            assert!(default.is_instance(<UffDefaultTorsionQueryError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                default
                    .getattr("source")?
                    .is_instance(<SmartsParseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );

            let error =
                UffBuilderError::from_rust(ck::UffBuilderError::SourceValencePrecondition {
                    atom_id: ck::AtomId::new(usize::MAX),
                    field: ck::UffPreparedValenceField::ImplicitHydrogen,
                    value: i32::MIN,
                });
            assert_eq!(
                error.value(py).getattr("atom_id")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(
                error.value(py).getattr("value")?.extract::<i32>()?,
                i32::MIN
            );
            assert!(
                error
                    .value(py)
                    .getattr("field")?
                    .extract::<UffPreparedValenceField>()?
                    == UffPreparedValenceField::ImplicitHydrogen
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .setattr("field", UffPreparedValenceField::Explicit)
                    .is_err()
            );

            for (native, index) in [
                (
                    ck::ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                        index: i16::MIN,
                        buffer_len: usize::MAX,
                    },
                    i64::from(i16::MIN),
                ),
                (
                    ck::ForceFieldKernelError::OopAddressOutsideDefinedSource {
                        index: i32::MIN,
                        buffer_len: usize::MAX,
                    },
                    i64::from(i32::MIN),
                ),
            ] {
                let error = ForceFieldKernelError::from_rust(native);
                assert_eq!(error.value(py).getattr("index")?.extract::<i64>()?, index);
                assert_eq!(
                    error.value(py).getattr("buffer_len")?.extract::<usize>()?,
                    usize::MAX
                );
                assert!(error.value(py).setattr("index", 0).is_err());
            }
            let error =
                ForceFieldKernelError::from_rust(ck::ForceFieldKernelError::IndexOutOfRange {
                    argument: ck::ForceFieldIndexArgument::J,
                    index: u32::MAX,
                    upper_bound: u32::MAX,
                });
            assert_eq!(
                error.value(py).getattr("index")?.extract::<u32>()?,
                u32::MAX
            );
            assert!(
                error
                    .value(py)
                    .getattr("argument")?
                    .extract::<ForceFieldIndexArgument>()?
                    == ForceFieldIndexArgument::J
            );
            let error = UffTypingError::from_rust(ck::UffTypingError::PreparedStateLength {
                input: ck::UffTypingInput::ConjugatedBondPresence,
                expected: usize::MAX,
                actual: 0,
            });
            assert_eq!(
                error.value(py).getattr("expected")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(
                error
                    .value(py)
                    .getattr("input")?
                    .extract::<UffTypingInput>()?
                    == UffTypingInput::ConjugatedBondPresence
            );
            let error = UffParamError::from_rust(ck::UffParamError::ParseFloat {
                line_number: usize::MAX,
                column_name: "r1",
                value: "x\0y".into(),
            });
            assert_eq!(
                error.value(py).getattr("line_number")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(
                error.value(py).getattr("value")?.extract::<String>()?,
                "x\0y"
            );
            Ok(())
        })
        .unwrap();
    }
}
