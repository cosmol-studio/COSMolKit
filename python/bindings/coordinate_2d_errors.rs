//! Exact 2D layout failures; borrowed native causes retain their real owner.
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::{ErrorRef, IoError};
use crate::match_errors::SubstructMatchError;
use crate::matrix_errors::MatrixError;
use crate::model_edit_errors::QueryAtomConversionError;
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::path_errors::PathError;
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::smarts_errors::SmartsParseError;
use crate::stereo_errors::LegacyStereoError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

fn project(node: ErrorRef) -> PyErr {
    if node.error().is::<ck::Coordinate2DError>() {
        return Coordinate2DError::from_ref(node);
    }
    if node.error().is::<ck::Coordinate2DTemplateError>() {
        return Coordinate2DTemplateError::from_ref(node);
    }
    if node.error().is::<ck::Coordinate2DLayoutError>() {
        return Coordinate2DLayoutError::from_ref(node);
    }
    if node.error().is::<std::io::Error>() {
        return IoError::from_ref(node);
    }
    macro_rules! leaf {
        ($($name:ident),+ $(,)?) => { $(
            if let Some(value) = node.error().downcast_ref::<ck::$name>() {
                return $name::from_rust(value.clone());
            }
        )+ };
    }
    leaf!(
        TopologyValidationError,
        CoordinateValidationError,
        QueryAtomConversionError,
        ValenceError,
        RingFindingError,
        LegacyStereoError,
        SmartsParseError,
        PathError,
        MatrixError,
        SubstructMatchError,
        PropertyUIntReadError
    );
    unreachable!("unmodeled concrete native 2D error source")
}

type TemplateSource = Source<
    Source<
        crate::error_framework::ErrorObject<SmartsParseError>,
        crate::error_framework::ErrorObject<TopologyValidationError>,
    >,
    Source<
        crate::error_framework::ErrorObject<QueryAtomConversionError>,
        Source<
            crate::error_framework::ErrorObject<RingFindingError>,
            Source<
                crate::error_framework::ErrorObject<PathError>,
                crate::error_framework::ErrorObject<IoError>,
            >,
        >,
    >,
>;
error_projection!(retained Coordinate2DTemplateError, Coordinate2DTemplateErrorKind, ck::Coordinate2DTemplateError, "depict", {
    InvalidDefaultRow { index, source }, InvalidTopology { index, source },
    NonElementIdentity { index, source }, RingInitialization { index, source },
    ConnectedComponents { index, source }, ExternalOpen { path, source },
    ExternalRead { path, source }, ExternalInvalidSmarts { path, line:usize, source },
    MissingCoordinates { row }, ThreeDimensionalCoordinates { row },
    MultipleFragments { row }, NotRingSystem { row },
}, cause=|inner:&RetainedError<ck::Coordinate2DTemplateError>|inner.node().source().map(project),
source=Option<TemplateSource>, methods={
#[getter]
#[doc = "Index associated with the failure. Present for `Coordinate2DTemplateError.InvalidDefaultRow`, `Coordinate2DTemplateError.InvalidTopology`, `Coordinate2DTemplateError.NonElementIdentity`, `Coordinate2DTemplateError.RingInitialization`, `Coordinate2DTemplateError.ConnectedComponents`; None when unavailable."] fn index(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DTemplateError::InvalidDefaultRow { index, .. }|ck::Coordinate2DTemplateError::InvalidTopology { index, .. }|ck::Coordinate2DTemplateError::NonElementIdentity { index, .. }|ck::Coordinate2DTemplateError::RingInitialization { index, .. }|ck::Coordinate2DTemplateError::ConnectedComponents { index, .. } => Some(*index), _=>None,
}}
#[getter]
#[doc = "Source path associated with the failure. Present for `Coordinate2DTemplateError.ExternalOpen`, `Coordinate2DTemplateError.ExternalRead`, `Coordinate2DTemplateError.ExternalInvalidSmarts`; None when unavailable."] fn path(&self)->Option<String>{ match self.inner.as_ref(){
ck::Coordinate2DTemplateError::ExternalOpen { path, .. }|ck::Coordinate2DTemplateError::ExternalRead { path, .. }|ck::Coordinate2DTemplateError::ExternalInvalidSmarts { path, .. } => Some(path.clone()), _=>None,
}}
#[getter]
#[doc = "Source row associated with the failure. Present for `Coordinate2DTemplateError.MissingCoordinates`, `Coordinate2DTemplateError.ThreeDimensionalCoordinates`, `Coordinate2DTemplateError.MultipleFragments`, `Coordinate2DTemplateError.NotRingSystem`; None when unavailable."] fn row(&self)->Option<String>{ match self.inner.as_ref(){
ck::Coordinate2DTemplateError::MissingCoordinates { row, .. }|ck::Coordinate2DTemplateError::ThreeDimensionalCoordinates { row, .. }|ck::Coordinate2DTemplateError::MultipleFragments { row, .. }|ck::Coordinate2DTemplateError::NotRingSystem { row, .. } => Some(row.clone()), _=>None,
}}
});
type LayoutSource = Source<
    crate::error_framework::ErrorObject<PropertyUIntReadError>,
    Source<
        crate::error_framework::ErrorObject<SubstructMatchError>,
        Source<
            crate::error_framework::ErrorObject<PathError>,
            crate::error_framework::ErrorObject<MatrixError>,
        >,
    >,
>;
error_projection!(retained Coordinate2DLayoutError,Coordinate2DLayoutErrorKind,ck::Coordinate2DLayoutError,"depict",{
    AtomIndexOutOfRange { atom, atom_count }, AtomAlreadyEmbedded { atom }, AtomNotEmbedded { atom },
    NotEnoughEmbeddedNeighbors { atom, count }, CoincidentPoints, InvalidAngle,
    NoCommonAtoms, MismatchedTopology, EmptyAttachment { atom },
    NonTetrahedralNoLigand { centre }, NonTetrahedralLigandOverflow { centre },
    CisTransBondInvalid { bond }, CollisionBondInvalid { bond },
    UndefinedSamplingDistance { first:usize, second:usize }, AtomCountTooLarge { atom_count },
    InvalidRankProperty { atom, key:&'static str, source }, GeometryNotEnoughNeighbors { atom, count },
    TemplateMatch(_), GraphPath(_), GraphDistance(_),
},cause=|inner:&RetainedError<ck::Coordinate2DLayoutError>|inner.node().source().map(project),
source=Option<LayoutSource>,methods={
#[getter]
#[doc = "Atom index involved in the failure. Present for `Coordinate2DLayoutError.AtomIndexOutOfRange`, `Coordinate2DLayoutError.AtomAlreadyEmbedded`, `Coordinate2DLayoutError.AtomNotEmbedded`, `Coordinate2DLayoutError.NotEnoughEmbeddedNeighbors`, `Coordinate2DLayoutError.EmptyAttachment`, `Coordinate2DLayoutError.InvalidRankProperty`, `Coordinate2DLayoutError.GeometryNotEnoughNeighbors`; None when unavailable."] fn atom(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::AtomIndexOutOfRange { atom, .. }|ck::Coordinate2DLayoutError::AtomAlreadyEmbedded { atom, .. }|ck::Coordinate2DLayoutError::AtomNotEmbedded { atom, .. }|ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors { atom, .. }|ck::Coordinate2DLayoutError::EmptyAttachment { atom, .. }|ck::Coordinate2DLayoutError::InvalidRankProperty { atom, .. }|ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors { atom, .. } => Some(*atom), _=>None,
}}
#[getter]
#[doc = "Number of atoms available to the operation. Present for `Coordinate2DLayoutError.AtomIndexOutOfRange`, `Coordinate2DLayoutError.AtomCountTooLarge`; None when unavailable."] fn atom_count(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::AtomIndexOutOfRange { atom_count, .. }|ck::Coordinate2DLayoutError::AtomCountTooLarge { atom_count, .. } => Some(*atom_count), _=>None,
}}
#[getter]
#[doc = "Count observed when the failure occurred. Present for `Coordinate2DLayoutError.NotEnoughEmbeddedNeighbors`, `Coordinate2DLayoutError.GeometryNotEnoughNeighbors`; None when unavailable."] fn count(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors { count, .. }|ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors { count, .. } => Some(*count), _=>None,
}}
#[getter]
#[doc = "Center atom involved in the coordinate-layout failure. Present for `Coordinate2DLayoutError.NonTetrahedralNoLigand`, `Coordinate2DLayoutError.NonTetrahedralLigandOverflow`; None when unavailable."] fn centre(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::NonTetrahedralNoLigand { centre, .. }|ck::Coordinate2DLayoutError::NonTetrahedralLigandOverflow { centre, .. } => Some(*centre), _=>None,
}}
#[getter]
#[doc = "Bond index involved in the failure. Present for `Coordinate2DLayoutError.CisTransBondInvalid`, `Coordinate2DLayoutError.CollisionBondInvalid`; None when unavailable."] fn bond(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::CisTransBondInvalid { bond, .. }|ck::Coordinate2DLayoutError::CollisionBondInvalid { bond, .. } => Some(*bond), _=>None,
}}
});
pub(crate) type CoordinateSource = Source<
    Source<
        crate::error_framework::ErrorObject<TopologyValidationError>,
        Source<
            crate::error_framework::ErrorObject<ValenceError>,
            crate::error_framework::ErrorObject<RingFindingError>,
        >,
    >,
    Source<
        crate::error_framework::ErrorObject<LegacyStereoError>,
        Source<
            crate::error_framework::ErrorObject<Coordinate2DTemplateError>,
            Source<
                crate::error_framework::ErrorObject<Coordinate2DLayoutError>,
                crate::error_framework::ErrorObject<CoordinateValidationError>,
            >,
        >,
    >,
>;
error_projection!(retained Coordinate2DError,Coordinate2DErrorKind,ck::Coordinate2DError,"depict",{
    InvalidTopology(_),PropertyCache(_),RingFinding(_),StereoAssignment(_),
    TemplateLoading(_),Fragment(_),CoordinateValidation(_),
    ConformerIdOverflow {max_id:usize},CoordGenUnavailable,
},cause=|inner:&RetainedError<ck::Coordinate2DError>|inner.node().source().map(project),source=Option<CoordinateSource>);

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn coordinate_error_graph_preserves_actual_layout_failure_and_retained_template_causes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let before = molecule.to_binary().unwrap();
            let options = ck::Coordinate2DParams {
                coordinate_map: BTreeMap::from([(usize::MAX, [1.0, 2.0])]),
                ..Default::default()
            };
            let native = match molecule
                .with_2d_coordinates_with_params(&options)
                .unwrap_err()
            {
                ck::OperationError::Coordinate2D(error) => error,
                error => panic!("unexpected native coordinate failure: {error:?}"),
            };
            let message = native.to_string();
            let error = Coordinate2DError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let layout = error.value(py).getattr("source")?;
            assert!(layout.is_instance(<Coordinate2DLayoutError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(layout.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(layout.getattr("atom")?.extract::<usize>()?, usize::MAX);
            assert!(layout.setattr("atom", 0).is_err());
            assert_eq!(molecule.to_binary().unwrap(), before);

            let native = ck::Coordinate2DError::TemplateLoading(
                ck::Coordinate2DTemplateError::ExternalOpen {
                    path: "template\\0原始".to_string(),
                    source: std::sync::Arc::new(std::io::Error::from_raw_os_error(13)),
                },
            );
            let message = native.to_string();
            let error = Coordinate2DError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let template = error
                .value(py)
                .getattr("source")?
                .extract::<crate::error_framework::ErrorObject<Coordinate2DTemplateError>>()?;
            let io = template
                .bind(py)
                .getattr("source")?
                .extract::<crate::error_framework::ErrorObject<IoError>>()?;
            assert!(io.bind(py).is(&template.bind(py).getattr("__cause__")?));
            drop(error);
            assert_eq!(
                template.bind(py).getattr("path")?.extract::<String>()?,
                "template\\0原始"
            );
            assert_eq!(io.bind(py).getattr("raw_os_error")?.extract::<i32>()?, 13);
            assert!(template.bind(py).setattr("path", "changed").is_err());

            let native =
                ck::Coordinate2DError::Fragment(ck::Coordinate2DLayoutError::InvalidRankProperty {
                    atom: usize::MAX,
                    key: "_CIPRank",
                    source: ck::PropertyUIntReadError::Negative { value: -7 },
                });
            let error = Coordinate2DError::from_rust(native);
            let layout = error.value(py).getattr("source")?;
            assert_eq!(layout.getattr("key")?.extract::<String>()?, "_CIPRank");
            let numeric = layout.getattr("source")?;
            assert!(numeric.is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(numeric.getattr("value")?.extract::<i64>()?, -7);
            assert!(numeric.is(&layout.getattr("__cause__")?));

            let error = Coordinate2DError::from_rust(ck::Coordinate2DError::ConformerIdOverflow {
                max_id: usize::MAX,
            });
            assert_eq!(
                error.value(py).getattr("max_id")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(error.value(py).getattr("source")?.is_none());
            let error = Coordinate2DTemplateError::from_rust(
                ck::Coordinate2DTemplateError::NotRingSystem {
                    row: "CC".to_string(),
                },
            );
            assert_eq!(error.value(py).getattr("row")?.extract::<String>()?, "CC");
            assert!(error.value(py).getattr("source")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
