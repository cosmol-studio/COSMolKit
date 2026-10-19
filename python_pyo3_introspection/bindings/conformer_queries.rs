//! Representation-only bounds and stereoisomer-count facade forwarding.
use crate::conformer_execution::{ConformerRunError, EnumerationError};
use crate::native::Molecule;
use crate::stereoisomer_options::{
    StereoisomerOptions, StereoisomerRandomSource, with_random_bits_errors,
};
use cosmolkit as ck;
use num_bigint::BigUint;
use pyo3::prelude::*;

pub(crate) fn bounds<'py>(
    molecule: &Molecule,
    py: Python<'py>,
) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
    let values = py
        .detach(|| molecule.inner.dg_bounds_matrix())
        .map_err(ConformerRunError::from_rust)?;
    numpy::PyArray2::from_vec2(py, &values).map_err(Into::into)
}
pub(crate) fn count(
    py: Python<'_>,
    molecule: &Molecule,
    options: Option<Py<StereoisomerOptions>>,
) -> PyResult<BigUint> {
    match options {
        Some(options) => {
            let (inner, owner) = {
                let options = options.borrow(py);
                (options.inner.clone(), options.random_source(py))
            };
            let result = with_random_bits_errors(|| {
                py.detach(|| molecule.inner.stereoisomer_count_with_options(&inner))
            })?;
            drop(owner);
            result.map_err(EnumerationError::from_rust)
        }
        None => with_random_bits_errors(|| py.detach(|| molecule.inner.stereoisomer_count()))?
            .map_err(EnumerationError::from_rust),
    }
}
pub(crate) fn count_keywords(
    py: Python<'_>,
    molecule: &Molecule,
    try_embedding: Option<bool>,
    only_unassigned: Option<bool>,
    max_isomers: Option<usize>,
    random_source: Option<Py<StereoisomerRandomSource>>,
    unique: Option<bool>,
    only_stereo_groups: Option<bool>,
) -> PyResult<BigUint> {
    if try_embedding.is_none()
        && only_unassigned.is_none()
        && max_isomers.is_none()
        && random_source.is_none()
        && unique.is_none()
        && only_stereo_groups.is_none()
    {
        return count(py, molecule, None);
    }
    let defaults = ck::StereoisomerOptions::default();
    let options = StereoisomerOptions::new(
        py,
        try_embedding.unwrap_or(defaults.try_embedding),
        only_unassigned.unwrap_or(defaults.only_unassigned),
        only_stereo_groups.unwrap_or(defaults.only_stereo_groups),
        max_isomers.unwrap_or(defaults.max_isomers),
        unique.unwrap_or(defaults.unique),
        random_source,
    );
    count(py, molecule, Some(Py::new(py, options)?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer_execution::{
        BoundsMatrixError, ConformerOptimizerError, EnumerationErrorKind, GenerationError,
        GenerationFailure, PoisonError,
    };
    use numpy::PyArrayMethods;
    use pyo3::types::PyDict;
    use std::error::Error;

    #[test]
    fn bounds_and_stereoisomer_queries_forward_options_arrays_and_retained_concrete_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for text in ["CCO", "FC(Cl)(Br)C", "C[C@H](F)Cl"] {
                let native = ck::Molecule::from_smiles(text).unwrap();
                let before = native.to_binary().unwrap();
                let expected = native.dg_bounds_matrix().unwrap();
                let expected_count = native.stereoisomer_count().unwrap();
                let molecule = Py::new(py, Molecule { inner: native })?;
                let matrix = molecule
                    .bind(py)
                    .call_method0("dg_bounds_matrix")?
                    .cast_into::<numpy::PyArray2<f64>>()?;
                assert_eq!(
                    matrix.readonly().as_array().shape(),
                    &[expected.len(), expected.len()]
                );
                let readonly = matrix.readonly();
                let array = readonly.as_array();
                for (actual, expected) in array.iter().zip(expected.iter().flatten()) {
                    assert_eq!(actual.to_bits(), expected.to_bits());
                }
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method0("stereoisomer_count")?
                        .extract::<BigUint>()?,
                    expected_count
                );
                let source = py
                    .get_type::<StereoisomerRandomSource>()
                    .call_method1("from_integer_seed", (num_bigint::BigInt::from(7),))?
                    .extract::<Py<StereoisomerRandomSource>>()?;
                let options = Py::new(
                    py,
                    StereoisomerOptions::new(
                        py,
                        false,
                        false,
                        false,
                        17,
                        false,
                        Some(source.clone_ref(py)),
                    ),
                )?;
                let expected = molecule
                    .borrow(py)
                    .inner
                    .stereoisomer_count_with_options(&options.borrow(py).inner)
                    .unwrap();
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method1("stereoisomer_count", (options.clone_ref(py),))?
                        .extract::<BigUint>()?,
                    expected
                );
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method1("stereoisomer_count_with_options", (options.clone_ref(py),))?
                        .extract::<BigUint>()?,
                    expected
                );
                let kwargs = PyDict::new(py);
                kwargs.set_item("only_unassigned", false)?;
                kwargs.set_item("max_isomers", 17)?;
                kwargs.set_item("random_source", source.clone_ref(py))?;
                kwargs.set_item("unique", false)?;
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method("stereoisomer_count", (), Some(&kwargs))?
                        .extract::<BigUint>()?,
                    expected
                );
                assert!(
                    molecule
                        .bind(py)
                        .call_method(
                            "stereoisomer_count",
                            (options.clone_ref(py),),
                            Some(&kwargs)
                        )
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            }

            let invalid = ck::Molecule::from_smiles_with_params(
                "[CH5]",
                &ck::SmilesParseParams {
                    sanitize: false,
                    remove_hs: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let message = invalid.dg_bounds_matrix().unwrap_err().to_string();
            let object = Py::new(py, Molecule { inner: invalid })?;
            let error = object
                .bind(py)
                .call_method0("dg_bounds_matrix")
                .unwrap_err();
            assert!(error.is_instance_of::<ConformerRunError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let query = error.value(py).getattr("source")?;
            assert!(query.is(&error.value(py).getattr("__cause__")?));
            assert!(
                query
                    .getattr("valence")?
                    .is_instance_of::<crate::valence_errors::ValenceError>()
            );
            drop(error);
            drop(object);
            assert!(query.getattr("source")?.is_none());

            let native = ck::EnumerationError::Embedding(ck::GenerationError::GenerationFailed(
                ck::GenerationFailure::TopologyBounds(ck::GraphBoundsError::Bounds(
                    ck::BoundsMatrixError::IndexOutOfBounds {
                        axis: ck::MatrixAxis::Column,
                        index: usize::MAX,
                        dimension: 9,
                    },
                )),
            ));
            let message = native.to_string();
            let error = EnumerationError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<EnumerationErrorKind>()?
                    == EnumerationErrorKind::Embedding
            );
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance_of::<BoundsMatrixError>());
            let embedding = error
                .value(py)
                .getattr("embedding")?
                .extract::<Py<GenerationError>>()?;
            let failure = embedding
                .bind(py)
                .getattr("generation_failed")?
                .extract::<Py<GenerationFailure>>()?;
            let bounds = failure.bind(py).getattr("topology_bounds")?;
            let retained = bounds
                .getattr("source")?
                .extract::<Py<BoundsMatrixError>>()?;
            drop(error);
            drop(embedding);
            drop(failure);
            drop(bounds);
            assert_eq!(
                retained.bind(py).getattr("index")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(
                retained
                    .bind(py)
                    .getattr("dimension")?
                    .extract::<Option<usize>>()?,
                Some(9)
            );
            assert!(retained.bind(py).setattr("index", 0).is_err());

            let kernel = ck::ForceFieldKernelError::NoPoints;
            let native = ck::ConformerOptimizerError::from(kernel);
            let error = ConformerOptimizerError::from_rust(native);
            let cause = error.value(py).getattr("cause")?;
            let value = cause.getattr("construction")?.getattr("kernel")?;
            assert!(value.is_instance_of::<crate::uff_parameter_errors::ForceFieldKernelError>());
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is(&error.value(py).getattr("__cause__")?)
            );
            assert!(cause.setattr("construction", py.None()).is_err());

            let native =
                ck::EnumerationError::RandomSourcePoisoned(std::sync::PoisonError::new(()));
            let message = native.source().unwrap().to_string();
            let error = EnumerationError::from_rust(native);
            let cause = error
                .value(py)
                .getattr("source")?
                .extract::<Py<PoisonError>>()?;
            drop(error);
            assert_eq!(cause.bind(py).str()?.to_str()?, message);
            assert!(cause.bind(py).getattr("__cause__")?.is_none());
            assert!(cause.bind(py).setattr("domain", "changed").is_err());
            Ok(())
        })
        .unwrap();
    }
}
