//! Exact PDB/XYZ writer failures; parent IO entrypoints remain separate work.
use crate::canonical_errors::PropertyStringError;
use crate::embedding_values::CoordinateDimension;
use crate::error_framework::{Source, error_projection};
use crate::kekulize_errors::KekulizeError;
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type PdbSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<CoordinateValidationError>,
        crate::error_framework::ErrorObject<KekulizeError>,
    >,
>;

error_projection!(PdbWriteError, PdbWriteErrorKind, ck::PdbWriteError, "io", {
    Property(_), ConformerNotFound { id }, Unsupported { feature: &'static str },
    Topology(_), Coordinates(_), AmbiguousCoordinates { two_d: usize, three_d: usize },
    MissingCoordinate { dimension: CoordinateDimension => |value: &ck::CoordinateDimension| (*value).into(), id },
    ConflictingCoordinateSelection, Kekulize(_),
}, cause = |inner: &ck::PdbWriteError| match inner {
    ck::PdbWriteError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::PdbWriteError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::PdbWriteError::Kekulize(value) => Some(KekulizeError::from_rust(value.clone())),
    // Native Property is transparent and its concrete child has no source.
    ck::PdbWriteError::Property(_) | ck::PdbWriteError::ConformerNotFound { .. }
    | ck::PdbWriteError::Unsupported { .. } | ck::PdbWriteError::AmbiguousCoordinates { .. }
    | ck::PdbWriteError::MissingCoordinate { .. } | ck::PdbWriteError::ConflictingCoordinateSelection => None,
}, source = Option<PdbSource>, methods = {
    #[getter]
    fn id(&self) -> Option<usize> { match self.inner {
        ck::PdbWriteError::ConformerNotFound { id }
        | ck::PdbWriteError::MissingCoordinate { id, .. } => Some(id),
        _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn property_error(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyStringError>>> {
        match self.inner {
            ck::PdbWriteError::Property(value) => PropertyStringError::from_rust(value)
                .value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});

type XyzSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<CoordinateValidationError>,
        crate::error_framework::ErrorObject<PropertyStringError>,
    >,
>;

error_projection!(XyzWriteError, XyzWriteErrorKind, ck::XyzWriteError, "io", {
    Topology(_), Coordinates(_), ConformerNotFound { id: usize }, Property(_),
}, cause = |inner: &ck::XyzWriteError| match inner {
    ck::XyzWriteError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::XyzWriteError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::XyzWriteError::Property(value) => Some(PropertyStringError::from_rust(*value)),
    ck::XyzWriteError::ConformerNotFound { .. } => None,
}, source = Option<XyzSource>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn real_writer_selection_failures_and_exact_causal_levels_remain_structured() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CC").unwrap();
            let before = molecule.to_binary().unwrap();
            let error = molecule
                .to_pdb_with_params(&ck::PdbWriteParams {
                    coordinate_selection: ck::MolCoordinateSelection::ThreeD { id: usize::MAX },
                    all_conformers: true,
                    ..Default::default()
                })
                .unwrap_err();
            let ck::MolecularIoError::PdbWrite(native) = error else {
                panic!("expected native PDB writer failure")
            };
            assert!(matches!(
                native,
                ck::PdbWriteError::ConflictingCoordinateSelection
            ));
            let display = native.to_string();
            let error = PdbWriteError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());

            // Native XYZ returns an empty block before checking an explicit
            // selector when no conformer exists (xyz.rs:376–377).
            assert_eq!(
                molecule
                    .to_xyz_with_params(&ck::XyzWriteParams {
                        conformer_id: Some(usize::MAX),
                        ..Default::default()
                    })
                    .unwrap(),
                ""
            );
            let xyz_source = molecule
                .with_only_3d_conformer(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
                .unwrap();
            let xyz_before = xyz_source.to_binary().unwrap();
            let error = xyz_source
                .to_xyz_with_params(&ck::XyzWriteParams {
                    conformer_id: Some(usize::MAX),
                    ..Default::default()
                })
                .unwrap_err();
            let ck::MolecularIoError::XyzWrite(native) = error else {
                panic!("expected native XYZ writer failure")
            };
            assert!(
                matches!(native, ck::XyzWriteError::ConformerNotFound { id } if id == usize::MAX)
            );
            let display = native.to_string();
            let error = XyzWriteError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("id")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(error.value(py).setattr("id", 0).is_err());
            assert!(error.cause(py).is_none());
            assert_eq!(molecule.to_binary().unwrap(), before);
            assert_eq!(xyz_source.to_binary().unwrap(), xyz_before);

            let property = ck::PropertyStringError::UnsupportedKind {
                kind: ck::PropertyValueKind::StringVector,
            };
            let transparent = ck::PdbWriteError::Property(property);
            assert!(transparent.source().is_none());
            let error = PdbWriteError::from_rust(transparent);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .getattr("property_error")?
                    .is_instance(<PropertyStringError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert!(
                error
                    .value(py)
                    .setattr("property_error", py.None())
                    .is_err()
            );
            let explicit = ck::XyzWriteError::Property(property);
            assert!(explicit.source().unwrap().is::<ck::PropertyStringError>());
            let error = XyzWriteError::from_rust(explicit);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance(<PropertyStringError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.getattr("__cause__")?.is_none());

            let error = PdbWriteError::from_rust(ck::PdbWriteError::MissingCoordinate {
                dimension: ck::CoordinateDimension::TwoD,
                id: usize::MAX,
            });
            assert_eq!(
                error.value(py).getattr("id")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(
                error
                    .value(py)
                    .getattr("dimension")?
                    .extract::<CoordinateDimension>()?
                    == ck::CoordinateDimension::TwoD.into()
            );
            assert!(error.value(py).setattr("dimension", py.None()).is_err());
            for error in [
                PdbWriteError::from_rust(ck::PdbWriteError::Coordinates(
                    ck::CoordinateValidationError::Missing3DConformer { id: usize::MAX },
                )),
                XyzWriteError::from_rust(ck::XyzWriteError::Coordinates(
                    ck::CoordinateValidationError::Missing3DConformer { id: usize::MAX },
                )),
            ] {
                let source = error.value(py).getattr("source")?;
                assert!(source.is_instance(<CoordinateValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                assert_eq!(source.getattr("id")?.extract::<usize>()?, usize::MAX);
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert!(error.value(py).setattr("source", py.None()).is_err());
            }
            Ok(())
        })
        .unwrap();
    }
}
