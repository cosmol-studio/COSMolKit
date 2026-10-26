use crate::attachment_values::TemplateAttachmentOrderError;
use crate::chemical_values::BondOrder;
use crate::error_framework::{Source, error_projection};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(MatrixError, MatrixErrorKind, ck::MatrixError, "matrix", {
    InvalidTopology(_), InvalidCoordinates(_), ActiveBondsWithoutAtoms,
    ActiveAtomOutOfRange { position, atom, atom_count },
    DuplicateActiveAtom { atom, first_position, second_position },
    ActiveBondOutOfRange { position, bond, bond_count },
    DuplicateActiveBond { bond, first_position, second_position },
    ActiveBondEndpointMissing { bond, endpoint: &'static str, atom },
    UnsupportedBondOrder { bond, order: BondOrder => |value: &ck::BondOrder| (*value).into() },
    No3dConformer, ConformerNotFound { conformer_id: usize },
    MatrixDimensionOverflow { dimension: usize },
}, cause = |inner: &ck::MatrixError| Python::attach(|py| match inner {
    // Native transparent wrappers expose the wrapped error's source, not the
    // wrapped error itself. Reuse the existing exhaustive model conversion.
    ck::MatrixError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::MatrixError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).cause(py),
    ck::MatrixError::ActiveBondsWithoutAtoms | ck::MatrixError::ActiveAtomOutOfRange { .. }
    | ck::MatrixError::DuplicateActiveAtom { .. } | ck::MatrixError::ActiveBondOutOfRange { .. }
    | ck::MatrixError::DuplicateActiveBond { .. } | ck::MatrixError::ActiveBondEndpointMissing { .. }
    | ck::MatrixError::UnsupportedBondOrder { .. } | ck::MatrixError::No3dConformer
    | ck::MatrixError::ConformerNotFound { .. } | ck::MatrixError::MatrixDimensionOverflow { .. } => None,
}), source = Option<Source<crate::error_framework::ErrorObject<StereoGroupError>, crate::error_framework::ErrorObject<TemplateAttachmentOrderError>>>, methods = {
    #[getter]
    #[doc = "Position reported by the native failure. Applies to MatrixError.ActiveAtomOutOfRange, MatrixError.ActiveBondOutOfRange. None when this payload is absent."] fn position(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::ActiveAtomOutOfRange { position, .. } | ck::MatrixError::ActiveBondOutOfRange { position, .. } => Some(position), _ => None }
    }
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to MatrixError.ActiveAtomOutOfRange, MatrixError.DuplicateActiveAtom, MatrixError.ActiveBondEndpointMissing. None when this payload is absent."] fn atom(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::ActiveAtomOutOfRange { atom, .. } | ck::MatrixError::DuplicateActiveAtom { atom, .. } | ck::MatrixError::ActiveBondEndpointMissing { atom, .. } => Some(atom.index()), _ => None }
    }
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to MatrixError.ActiveBondOutOfRange, MatrixError.DuplicateActiveBond, MatrixError.ActiveBondEndpointMissing, MatrixError.UnsupportedBondOrder. None when this payload is absent."] fn bond(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::ActiveBondOutOfRange { bond, .. } | ck::MatrixError::DuplicateActiveBond { bond, .. } | ck::MatrixError::ActiveBondEndpointMissing { bond, .. } | ck::MatrixError::UnsupportedBondOrder { bond, .. } => Some(bond.index()), _ => None }
    }
    #[getter]
    #[doc = "Atom count used when validating the reported atom reference. Applies to MatrixError.ActiveAtomOutOfRange. None when this payload is absent."] fn atom_count(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::ActiveAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None }
    }
    #[getter]
    #[doc = "Bond count used when validating the reported bond reference. Applies to MatrixError.ActiveBondOutOfRange. None when this payload is absent."] fn bond_count(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::ActiveBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None }
    }
    #[getter]
    #[doc = "First position involved in the reported conflict. Applies to MatrixError.DuplicateActiveAtom, MatrixError.DuplicateActiveBond. None when this payload is absent."] fn first_position(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::DuplicateActiveAtom { first_position, .. } | ck::MatrixError::DuplicateActiveBond { first_position, .. } => Some(first_position), _ => None }
    }
    #[getter]
    #[doc = "Second position involved in the reported conflict. Applies to MatrixError.DuplicateActiveAtom, MatrixError.DuplicateActiveBond. None when this payload is absent."] fn second_position(&self) -> Option<usize> {
        match self.inner { ck::MatrixError::DuplicateActiveAtom { second_position, .. } | ck::MatrixError::DuplicateActiveBond { second_position, .. } => Some(second_position), _ => None }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured TopologyValidationError payload of MatrixError.InvalidTopology; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> {
        match &self.inner { ck::MatrixError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract::<crate::error_framework::ErrorObject<TopologyValidationError>>().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured CoordinateValidationError payload of MatrixError.InvalidCoordinates; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn coordinates(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>> {
        match &self.inner { ck::MatrixError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract::<crate::error_framework::ErrorObject<CoordinateValidationError>>().map(Some).map_err(Into::into), _ => Ok(None) }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use crate::params::{DistanceMatrix3dParams, DistanceMatrixParams};
    use pyo3::types::PyDict;

    #[test]
    fn native_matrix_calls_keep_bits_and_transparent_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CC=O").unwrap();
            let mut builder = native.to_builder();
            builder
                .add_3d_conformer(vec![[0., 0., 0.], [1.9, 0.2, 0.], [2.9, 0.4, 0.1]])
                .unwrap();
            let native = builder.build().unwrap();
            let molecule = Py::new(
                py,
                Molecule {
                    inner: native.clone(),
                },
            )?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("use_atom_weights", true)?;
            kwargs.set_item("use_bond_order", true)?;
            let params = ck::DistanceMatrixParams {
                use_atom_weights: true,
                use_bond_order: true,
            };
            let expected = native.distance_matrix_with_params(&params).unwrap();
            let options = Py::new(
                py,
                DistanceMatrixParams {
                    inner: crate::parameters::NativeField::Owned(params),
                },
            )?;
            for result in [
                molecule
                    .bind(py)
                    .call_method("distance_matrix", (), Some(&kwargs))?,
                molecule
                    .bind(py)
                    .call_method1("distance_matrix", (options.clone_ref(py),))?,
                molecule
                    .bind(py)
                    .call_method1("distance_matrix_with_params", (options.clone_ref(py),))?,
            ] {
                assert_eq!(
                    result
                        .getattr("values")?
                        .extract::<Vec<f64>>()?
                        .iter()
                        .map(|value| value.to_bits())
                        .collect::<Vec<_>>(),
                    expected
                        .values()
                        .iter()
                        .map(|value| value.to_bits())
                        .collect::<Vec<_>>()
                );
            }
            assert!(
                molecule
                    .bind(py)
                    .call_method("distance_matrix", (options,), Some(&kwargs))
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
            );
            let params = ck::DistanceMatrix3dParams {
                conformer_id: None,
                use_atom_weights: true,
            };
            let expected = native.distance_matrix_3d_with_params(&params).unwrap();
            let options = Py::new(
                py,
                DistanceMatrix3dParams {
                    inner: crate::parameters::NativeField::Owned(params),
                },
            )?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("use_atom_weights", true)?;
            for result in [
                molecule
                    .bind(py)
                    .call_method("distance_matrix_3d", (), Some(&kwargs))?,
                molecule
                    .bind(py)
                    .call_method1("distance_matrix_3d", (options.clone_ref(py),))?,
                molecule
                    .bind(py)
                    .call_method1("distance_matrix_3d_with_params", (options,))?,
            ] {
                assert_eq!(
                    result
                        .getattr("values")?
                        .extract::<Vec<f64>>()?
                        .iter()
                        .map(|value| value.to_bits())
                        .collect::<Vec<_>>(),
                    expected
                        .values()
                        .iter()
                        .map(|value| value.to_bits())
                        .collect::<Vec<_>>()
                );
            }
            let native = ck::MatrixError::InvalidTopology(
                ck::TopologyValidationError::StereoGroup(ck::StereoGroupError::DuplicateAtom),
            );
            assert!(native.source().unwrap().is::<ck::StereoGroupError>());
            let display = native.to_string();
            let error = MatrixError::from_rust(native);
            assert!(
                error
                    .cause(py)
                    .unwrap()
                    .is_instance(py, <StereoGroupError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())
            );
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is(error.cause(py).unwrap().value(py))
            );
            assert!(
                error
                    .value(py)
                    .getattr("topology")?
                    .is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let native = ck::MatrixError::InvalidCoordinates(
                ck::CoordinateValidationError::Missing3DConformer { id: 7 },
            );
            assert!(native.source().is_none());
            let error = MatrixError::from_rust(native);
            assert!(error.cause(py).is_none());
            assert_eq!(
                error
                    .value(py)
                    .getattr("coordinates")?
                    .getattr("id")?
                    .extract::<usize>()?,
                7
            );
            Ok(())
        })
        .unwrap();
    }
}
