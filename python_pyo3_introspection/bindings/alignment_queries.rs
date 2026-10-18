//! Readonly alignment: native configuration transport, never alignment math.
use crate::alignment_values::{
    AllConformerRmsdParameters, BestAlignmentParameters, CoordinateRmsdParameters,
};
use crate::collection_params::{AlignmentAtomMap, AlignmentParameters};
use cosmolkit as ck;

pub(crate) fn transform_keywords(
    probe_conformer_id: Option<i32>,
    reference_conformer_id: Option<i32>,
    atom_map: Option<Vec<AlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    reflect: Option<bool>,
    max_iterations: Option<u32>,
) -> AlignmentParameters {
    let defaults = ck::AlignmentParameters::default();
    AlignmentParameters {
        inner: ck::AlignmentParameters {
            probe_conformer_id: probe_conformer_id.unwrap_or(defaults.probe_conformer_id),
            reference_conformer_id: reference_conformer_id
                .unwrap_or(defaults.reference_conformer_id),
            atom_map: atom_map.map(|items| items.into_iter().map(|value| value.inner).collect()),
            weights,
            reflect: reflect.unwrap_or(defaults.reflect),
            max_iterations: max_iterations.unwrap_or(defaults.max_iterations),
        },
    }
}
pub(crate) fn best_keywords(
    probe_conformer_id: Option<i32>,
    reference_conformer_id: Option<i32>,
    atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
    weights: Option<Vec<f64>>,
    reflect: Option<bool>,
    max_iterations: Option<u32>,
    max_matches: Option<i32>,
    symmetrize_conjugated_terminal_groups: Option<bool>,
    ignore_hydrogens: Option<bool>,
    num_threads: Option<i32>,
) -> BestAlignmentParameters {
    let defaults = ck::BestAlignmentParameters::default();
    BestAlignmentParameters {
        inner: ck::BestAlignmentParameters {
            probe_conformer_id: probe_conformer_id.unwrap_or(defaults.probe_conformer_id),
            reference_conformer_id: reference_conformer_id
                .unwrap_or(defaults.reference_conformer_id),
            atom_maps: atom_maps
                .map(|maps| {
                    maps.into_iter()
                        .map(|row| row.into_iter().map(|value| value.inner).collect())
                        .collect()
                })
                .unwrap_or(defaults.atom_maps),
            weights,
            reflect: reflect.unwrap_or(defaults.reflect),
            max_iterations: max_iterations.unwrap_or(defaults.max_iterations),
            max_matches: max_matches.unwrap_or(defaults.max_matches),
            symmetrize_conjugated_terminal_groups: symmetrize_conjugated_terminal_groups
                .unwrap_or(defaults.symmetrize_conjugated_terminal_groups),
            ignore_hydrogens: ignore_hydrogens.unwrap_or(defaults.ignore_hydrogens),
            num_threads: num_threads.unwrap_or(defaults.num_threads),
        },
    }
}
pub(crate) fn coordinate_keywords(
    probe_conformer_id: Option<i32>,
    reference_conformer_id: Option<i32>,
    atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
    weights: Option<Vec<f64>>,
    max_matches: Option<i32>,
    symmetrize_conjugated_terminal_groups: Option<bool>,
) -> CoordinateRmsdParameters {
    let defaults = ck::CoordinateRmsdParameters::default();
    CoordinateRmsdParameters {
        inner: ck::CoordinateRmsdParameters {
            probe_conformer_id: probe_conformer_id.unwrap_or(defaults.probe_conformer_id),
            reference_conformer_id: reference_conformer_id
                .unwrap_or(defaults.reference_conformer_id),
            atom_maps: atom_maps
                .map(|maps| {
                    maps.into_iter()
                        .map(|row| row.into_iter().map(|value| value.inner).collect())
                        .collect()
                })
                .unwrap_or(defaults.atom_maps),
            weights,
            max_matches: max_matches.unwrap_or(defaults.max_matches),
            symmetrize_conjugated_terminal_groups: symmetrize_conjugated_terminal_groups
                .unwrap_or(defaults.symmetrize_conjugated_terminal_groups),
        },
    }
}
pub(crate) fn all_keywords(
    atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
    weights: Option<Vec<f64>>,
    max_matches: Option<i32>,
    symmetrize_conjugated_terminal_groups: Option<bool>,
    ignore_hydrogens: Option<bool>,
    num_threads: Option<i32>,
) -> AllConformerRmsdParameters {
    let defaults = ck::AllConformerRmsdParameters::default();
    AllConformerRmsdParameters {
        inner: ck::AllConformerRmsdParameters {
            atom_maps: atom_maps
                .map(|maps| {
                    maps.into_iter()
                        .map(|row| row.into_iter().map(|value| value.inner).collect())
                        .collect()
                })
                .unwrap_or(defaults.atom_maps),
            weights,
            max_matches: max_matches.unwrap_or(defaults.max_matches),
            symmetrize_conjugated_terminal_groups: symmetrize_conjugated_terminal_groups
                .unwrap_or(defaults.symmetrize_conjugated_terminal_groups),
            ignore_hydrogens: ignore_hydrogens.unwrap_or(defaults.ignore_hydrogens),
            num_threads: num_threads.unwrap_or(defaults.num_threads),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment_errors::{AlignmentError, AlignmentErrorKind, ThreadCountError};
    use crate::alignment_values::{AlignmentResult, ConformerRmsd};
    use crate::native::Molecule;
    use pyo3::prelude::*;
    use pyo3::types::PyDict;

    fn molecule(shift: f64) -> ck::Molecule {
        let molecule = ck::Molecule::from_smiles("CCO").unwrap();
        let mut builder = molecule.to_builder();
        builder
            .add_3d_conformer(vec![
                [shift, 0., 0.],
                [shift + 1., 1., 0.],
                [shift + 2., 0., 1.],
            ])
            .unwrap();
        builder
            .add_3d_conformer(vec![
                [shift, 1., 0.],
                [shift + 1., 2., 0.],
                [shift + 2., 1., 1.],
            ])
            .unwrap();
        builder.build().unwrap()
    }
    fn exact_result(actual: &ck::AlignmentResult, expected: &ck::AlignmentResult) {
        assert_eq!(actual.rmsd().to_bits(), expected.rmsd().to_bits());
        assert_eq!(actual.atom_map(), expected.atom_map());
        for (a, b) in actual
            .transform()
            .matrix()
            .iter()
            .flatten()
            .zip(expected.transform().matrix().iter().flatten())
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    #[test]
    fn complete_readonly_alignment_family_keeps_exact_native_results_options_errors_and_source_blocks()
     {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let probe_native = molecule(3.);
            let reference_native = molecule(0.);
            let before_probe = probe_native.to_binary().unwrap();
            let before_reference = reference_native.to_binary().unwrap();
            let probe = Py::new(
                py,
                Molecule {
                    inner: probe_native,
                },
            )?;
            let reference = Py::new(
                py,
                Molecule {
                    inner: reference_native,
                },
            )?;

            macro_rules! result_query {
                ($name:ident, $alias:ident, $params:ident) => {{
                    let expected = probe
                        .borrow(py)
                        .inner
                        .$name(&reference.borrow(py).inner)
                        .unwrap();
                    let actual = probe
                        .bind(py)
                        .call_method1(stringify!($name), (reference.clone_ref(py),))?
                        .extract::<Py<AlignmentResult>>()?;
                    exact_result(&actual.borrow(py).inner, &expected);
                    let native = ck::$params {
                        reflect: true,
                        max_iterations: 20,
                        ..Default::default()
                    };
                    let expected = probe
                        .borrow(py)
                        .inner
                        .$alias(&reference.borrow(py).inner, &native)
                        .unwrap();
                    let params = Py::new(py, $params { inner: native })?;
                    let kwargs = PyDict::new(py);
                    kwargs.set_item("reflect", true)?;
                    kwargs.set_item("max_iterations", 20)?;
                    for actual in [
                        probe.bind(py).call_method1(
                            stringify!($name),
                            (reference.clone_ref(py), params.clone_ref(py)),
                        )?,
                        probe.bind(py).call_method1(
                            stringify!($alias),
                            (reference.clone_ref(py), params.clone_ref(py)),
                        )?,
                        probe.bind(py).call_method(
                            stringify!($name),
                            (reference.clone_ref(py),),
                            Some(&kwargs),
                        )?,
                    ] {
                        exact_result(
                            &actual.extract::<Py<AlignmentResult>>()?.borrow(py).inner,
                            &expected,
                        );
                    }
                    assert!(
                        probe
                            .bind(py)
                            .call_method(
                                stringify!($name),
                                (reference.clone_ref(py), params),
                                Some(&kwargs)
                            )
                            .unwrap_err()
                            .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                    );
                }};
            }
            result_query!(
                alignment_transform_to,
                alignment_transform_to_with_params,
                AlignmentParameters
            );
            result_query!(
                best_alignment_to,
                best_alignment_to_with_params,
                BestAlignmentParameters
            );

            macro_rules! scalar_query {
                ($name:ident, $alias:ident, $params:ident) => {{
                    let expected = probe
                        .borrow(py)
                        .inner
                        .$name(&reference.borrow(py).inner)
                        .unwrap();
                    assert_eq!(
                        probe
                            .bind(py)
                            .call_method1(stringify!($name), (reference.clone_ref(py),))?
                            .extract::<f64>()?
                            .to_bits(),
                        expected.to_bits()
                    );
                    let native = ck::$params {
                        weights: Some(vec![1., 2., 3.]),
                        ..Default::default()
                    };
                    let expected = probe
                        .borrow(py)
                        .inner
                        .$alias(&reference.borrow(py).inner, &native)
                        .unwrap();
                    let params = Py::new(py, $params { inner: native })?;
                    let kwargs = PyDict::new(py);
                    kwargs.set_item("weights", vec![1., 2., 3.])?;
                    for actual in [
                        probe.bind(py).call_method1(
                            stringify!($name),
                            (reference.clone_ref(py), params.clone_ref(py)),
                        )?,
                        probe.bind(py).call_method1(
                            stringify!($alias),
                            (reference.clone_ref(py), params.clone_ref(py)),
                        )?,
                        probe.bind(py).call_method(
                            stringify!($name),
                            (reference.clone_ref(py),),
                            Some(&kwargs),
                        )?,
                    ] {
                        assert_eq!(actual.extract::<f64>()?.to_bits(), expected.to_bits());
                    }
                    assert!(
                        probe
                            .bind(py)
                            .call_method(
                                stringify!($name),
                                (reference.clone_ref(py), params),
                                Some(&kwargs)
                            )
                            .unwrap_err()
                            .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                    );
                }};
            }
            scalar_query!(
                best_rmsd_to,
                best_rmsd_to_with_params,
                BestAlignmentParameters
            );
            scalar_query!(
                coordinate_rmsd_to,
                coordinate_rmsd_to_with_params,
                CoordinateRmsdParameters
            );
            for num_threads in [1, 2] {
                let native = ck::AllConformerRmsdParameters {
                    num_threads,
                    ..Default::default()
                };
                let expected = probe
                    .borrow(py)
                    .inner
                    .all_conformer_best_rmsds_with_params(&native)
                    .unwrap();
                let params = Py::new(py, AllConformerRmsdParameters { inner: native })?;
                let kwargs = PyDict::new(py);
                kwargs.set_item("num_threads", num_threads)?;
                let outputs = [
                    probe
                        .bind(py)
                        .call_method1("all_conformer_best_rmsds", (params.clone_ref(py),))?,
                    probe.bind(py).call_method1(
                        "all_conformer_best_rmsds_with_params",
                        (params.clone_ref(py),),
                    )?,
                    probe
                        .bind(py)
                        .call_method("all_conformer_best_rmsds", (), Some(&kwargs))?,
                    probe.bind(py).call_method0("all_conformer_best_rmsds")?,
                ];
                for actual in outputs {
                    let actual = actual.extract::<Vec<ConformerRmsd>>()?;
                    assert_eq!(actual.len(), expected.len());
                    for (actual, expected) in actual.iter().zip(&expected) {
                        assert_eq!(
                            actual.inner.probe_conformer_id(),
                            expected.probe_conformer_id()
                        );
                        assert_eq!(
                            actual.inner.reference_conformer_id(),
                            expected.reference_conformer_id()
                        );
                        assert_eq!(actual.inner.rmsd().to_bits(), expected.rmsd().to_bits());
                    }
                }
                assert!(
                    probe
                        .bind(py)
                        .call_method("all_conformer_best_rmsds", (params,), Some(&kwargs))
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );
            }
            let kwargs = PyDict::new(py);
            kwargs.set_item("atom_map", Vec::<AlignmentAtomMap>::new())?;
            let error = probe
                .bind(py)
                .call_method(
                    "alignment_transform_to",
                    (reference.clone_ref(py),),
                    Some(&kwargs),
                )
                .unwrap_err();
            assert!(error.is_instance_of::<AlignmentError>(py));
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<AlignmentErrorKind>()?
                    == AlignmentErrorKind::EmptyAtomMap
            );
            assert!(error.value(py).getattr("source")?.is_none());
            let kwargs = PyDict::new(py);
            kwargs.set_item("probe_conformer_id", i32::MAX)?;
            let error = probe
                .bind(py)
                .call_method(
                    "coordinate_rmsd_to",
                    (reference.clone_ref(py),),
                    Some(&kwargs),
                )
                .unwrap_err();
            assert!(error.is_instance_of::<AlignmentError>(py));
            assert_eq!(error.value(py).getattr("id")?.extract::<i32>()?, i32::MAX);
            assert!(error.value(py).setattr("id", 0).is_err());
            let kwargs = PyDict::new(py);
            kwargs.set_item("num_threads", i32::MIN)?;
            let error = probe
                .bind(py)
                .call_method("all_conformer_best_rmsds", (), Some(&kwargs))
                .unwrap_err();
            assert!(error.is_instance_of::<AlignmentError>(py));
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance_of::<ThreadCountError>()
            );
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is(&error.value(py).getattr("__cause__")?)
            );
            assert_eq!(probe.borrow(py).inner.to_binary().unwrap(), before_probe);
            assert_eq!(
                reference.borrow(py).inner.to_binary().unwrap(),
                before_reference
            );
            drop(probe);
            drop(reference);
            assert!(!error.value(py).str()?.is_empty()?);
            Ok(())
        })
        .unwrap();
    }
}
