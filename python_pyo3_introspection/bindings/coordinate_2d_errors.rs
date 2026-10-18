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
use pyo3::exceptions::PyValueError;
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
    Source<Py<SmartsParseError>, Py<TopologyValidationError>>,
    Source<
        Py<QueryAtomConversionError>,
        Source<Py<RingFindingError>, Source<Py<PathError>, Py<IoError>>>,
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
fn index(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DTemplateError::InvalidDefaultRow { index, .. }|ck::Coordinate2DTemplateError::InvalidTopology { index, .. }|ck::Coordinate2DTemplateError::NonElementIdentity { index, .. }|ck::Coordinate2DTemplateError::RingInitialization { index, .. }|ck::Coordinate2DTemplateError::ConnectedComponents { index, .. } => Some(*index), _=>None,
}}
#[getter]
fn path(&self)->Option<String>{ match self.inner.as_ref(){
ck::Coordinate2DTemplateError::ExternalOpen { path, .. }|ck::Coordinate2DTemplateError::ExternalRead { path, .. }|ck::Coordinate2DTemplateError::ExternalInvalidSmarts { path, .. } => Some(path.clone()), _=>None,
}}
#[getter]
fn row(&self)->Option<String>{ match self.inner.as_ref(){
ck::Coordinate2DTemplateError::MissingCoordinates { row, .. }|ck::Coordinate2DTemplateError::ThreeDimensionalCoordinates { row, .. }|ck::Coordinate2DTemplateError::MultipleFragments { row, .. }|ck::Coordinate2DTemplateError::NotRingSystem { row, .. } => Some(row.clone()), _=>None,
}}
});
type LayoutSource = Source<
    Py<PropertyUIntReadError>,
    Source<Py<SubstructMatchError>, Source<Py<PathError>, Py<MatrixError>>>,
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
fn atom(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::AtomIndexOutOfRange { atom, .. }|ck::Coordinate2DLayoutError::AtomAlreadyEmbedded { atom, .. }|ck::Coordinate2DLayoutError::AtomNotEmbedded { atom, .. }|ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors { atom, .. }|ck::Coordinate2DLayoutError::EmptyAttachment { atom, .. }|ck::Coordinate2DLayoutError::InvalidRankProperty { atom, .. }|ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors { atom, .. } => Some(*atom), _=>None,
}}
#[getter]
fn atom_count(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::AtomIndexOutOfRange { atom_count, .. }|ck::Coordinate2DLayoutError::AtomCountTooLarge { atom_count, .. } => Some(*atom_count), _=>None,
}}
#[getter]
fn count(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors { count, .. }|ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors { count, .. } => Some(*count), _=>None,
}}
#[getter]
fn centre(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::NonTetrahedralNoLigand { centre, .. }|ck::Coordinate2DLayoutError::NonTetrahedralLigandOverflow { centre, .. } => Some(*centre), _=>None,
}}
#[getter]
fn bond(&self)->Option<usize>{ match self.inner.as_ref(){
ck::Coordinate2DLayoutError::CisTransBondInvalid { bond, .. }|ck::Coordinate2DLayoutError::CollisionBondInvalid { bond, .. } => Some(*bond), _=>None,
}}
});
pub(crate) type CoordinateSource = Source<
    Source<Py<TopologyValidationError>, Source<Py<ValenceError>, Py<RingFindingError>>>,
    Source<
        Py<LegacyStereoError>,
        Source<
            Py<Coordinate2DTemplateError>,
            Source<Py<Coordinate2DLayoutError>, Py<CoordinateValidationError>>,
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
            assert!(layout.is_instance_of::<Coordinate2DLayoutError>());
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
                .extract::<Py<Coordinate2DTemplateError>>()?;
            let io = template
                .bind(py)
                .getattr("source")?
                .extract::<Py<IoError>>()?;
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
            assert!(numeric.is_instance_of::<PropertyUIntReadError>());
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
