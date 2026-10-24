//! Concrete detached transform errors and their actual model sources.
use crate::error_framework::{Source, error_projection};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(TransformError, TransformErrorKind, ck::TransformError, "coordinates", {
    InvalidTopology(_), InvalidCoordinates(_),
    AtomOutOfRange { role, atom, atom_count: usize },
    WeightCount { actual: usize, required: usize },
    NonFiniteWeight { atom, value },
    InvalidWeightSum { sum: f64 }, NoSelectedAtoms,
    NonFinitePoint { role, axis: &'static str },
    NonFiniteTransform { row: usize, column: usize },
    EigenDecompositionDidNotConverge, No3dConformer,
    ConformerNotFound { conformer_id: usize },
    AtomsNotBonded { first_role, first, second_role, second },
    RingBondNotMovable { bond: usize => |id: &ck::BondId| id.index() },
    BothAngleBondsInRing { first_bond: usize => |id: &ck::BondId| id.index(), second_bond: usize => |id: &ck::BondId| id.index() },
    CoincidentCoordinates { first_role, first, second_role, second },
    UndefinedRotationAxis,
    InvalidTargetValue { quantity: &'static str, value },
}, cause = |inner: &ck::TransformError| match inner {
    ck::TransformError::InvalidTopology(source) => Some(TopologyValidationError::from_rust(source.clone())),
    ck::TransformError::InvalidCoordinates(source) => Some(CoordinateValidationError::from_rust(source.clone())),
    ck::TransformError::AtomOutOfRange { .. } | ck::TransformError::WeightCount { .. }
    | ck::TransformError::NonFiniteWeight { .. } | ck::TransformError::InvalidWeightSum { .. }
    | ck::TransformError::NoSelectedAtoms | ck::TransformError::NonFinitePoint { .. }
    | ck::TransformError::NonFiniteTransform { .. } | ck::TransformError::EigenDecompositionDidNotConverge
    | ck::TransformError::No3dConformer | ck::TransformError::ConformerNotFound { .. }
    | ck::TransformError::AtomsNotBonded { .. } | ck::TransformError::RingBondNotMovable { .. }
    | ck::TransformError::BothAngleBondsInRing { .. } | ck::TransformError::CoincidentCoordinates { .. }
    | ck::TransformError::UndefinedRotationAxis | ck::TransformError::InvalidTargetValue { .. } => None,
}, source = Option<Source<crate::error_framework::ErrorObject<TopologyValidationError>, crate::error_framework::ErrorObject<CoordinateValidationError>>>, methods = {
    #[getter]
    fn role(&self) -> Option<&'static str> {
        match self.inner {
            ck::TransformError::AtomOutOfRange { role, .. } | ck::TransformError::NonFinitePoint { role, .. } => Some(role),
            _ => None,
        }
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        match self.inner {
            ck::TransformError::AtomOutOfRange { atom, .. } | ck::TransformError::NonFiniteWeight { atom, .. } => Some(atom.index()),
            _ => None,
        }
    }
    #[getter]
    fn value(&self) -> Option<f64> {
        match self.inner {
            ck::TransformError::NonFiniteWeight { value, .. } | ck::TransformError::InvalidTargetValue { value, .. } => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn first_role(&self) -> Option<&'static str> {
        match self.inner {
            ck::TransformError::AtomsNotBonded { first_role, .. } | ck::TransformError::CoincidentCoordinates { first_role, .. } => Some(first_role),
            _ => None,
        }
    }
    #[getter]
    fn first(&self) -> Option<usize> {
        match self.inner {
            ck::TransformError::AtomsNotBonded { first, .. } | ck::TransformError::CoincidentCoordinates { first, .. } => Some(first.index()),
            _ => None,
        }
    }
    #[getter]
    fn second_role(&self) -> Option<&'static str> {
        match self.inner {
            ck::TransformError::AtomsNotBonded { second_role, .. } | ck::TransformError::CoincidentCoordinates { second_role, .. } => Some(second_role),
            _ => None,
        }
    }
    #[getter]
    fn second(&self) -> Option<usize> {
        match self.inner {
            ck::TransformError::AtomsNotBonded { second, .. } | ck::TransformError::CoincidentCoordinates { second, .. } => Some(second.index()),
            _ => None,
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_transform_failure_keeps_typed_model_sources_and_payload_bits() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CC").unwrap();
            let failure = molecule
                .with_atom_position(ck::AtomId::new(0), [0.; 3])
                .unwrap_err();
            let ck::OperationError::Transform(source) = failure else {
                panic!("expected native transform cause")
            };
            let expected = source.to_string();
            let error = TransformError::from_rust(source);
            assert!(error.is_instance(py, <TransformError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).str()?.to_str()?, expected);
            assert!(error.cause(py).is_none());
            for source in [
                ck::TransformError::InvalidTopology(ck::TopologyValidationError::AdjacencyMismatch),
                ck::TransformError::InvalidCoordinates(ck::CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: 3,
                    rows: 1,
                    atom_count: 2,
                }),
            ] {
                let error = TransformError::from_rust(source);
                let cause = error.cause(py).unwrap();
                assert!(error.value(py).getattr("source")?.is(cause.value(py)));
                assert!(
                    cause.is_instance(py, <TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())
                        || cause.is_instance(py, <CoordinateValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())
                );
                assert!(error.value(py).setattr("source", py.None()).is_err());
            }
            let bits = 0x7ff8_0000_0000_0042;
            let error = TransformError::from_rust(ck::TransformError::NonFiniteWeight {
                atom: ck::AtomId::new(usize::MAX),
                value: f64::from_bits(bits),
            });
            assert_eq!(
                error
                    .value(py)
                    .getattr("atom")?
                    .extract::<Option<usize>>()?,
                Some(usize::MAX)
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("value")?
                    .extract::<f64>()?
                    .to_bits(),
                bits
            );
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
