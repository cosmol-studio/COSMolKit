//! Borrow actual molecule inputs; the facade owns the MCS search and result.
use crate::mcs::McsResult;
use crate::mcs_errors::McsError;
use crate::native::Molecule;
use crate::params::{McsAtomCompareParameters, McsBondCompareParameters};
use crate::search_reaction_params::{
    McsAtomComparatorInput, McsBondComparatorInput, McsParameters,
};
use cosmolkit as ck;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use std::sync::Arc;

#[pyfunction]
#[pyo3(signature=(inputs, params=crate::arguments::Argument::Missing, *, store_all=crate::arguments::Argument::Missing, maximize_bonds=crate::arguments::Argument::Missing, threshold=crate::arguments::Argument::Missing, timeout=crate::arguments::Argument::Missing, verbose=crate::arguments::Argument::Missing, atom_compare_parameters=crate::arguments::Argument::Missing, bond_compare_parameters=crate::arguments::Argument::Missing, atom_comparator=crate::arguments::Argument::Missing, bond_comparator=crate::arguments::Argument::Missing, initial_seed=crate::arguments::Argument::Missing))]
pub fn maximum_common_substructure(
    py: Python<'_>,
    inputs: Vec<Py<Molecule>>,
    params: crate::arguments::Argument<Option<PyRef<'_, McsParameters>>>,
    store_all: crate::arguments::Argument<Option<bool>>,
    maximize_bonds: crate::arguments::Argument<Option<bool>>,
    threshold: crate::arguments::Argument<Option<f64>>,
    timeout: crate::arguments::Argument<Option<u32>>,
    verbose: crate::arguments::Argument<Option<bool>>,
    atom_compare_parameters: crate::arguments::Argument<
        Option<PyRef<'_, McsAtomCompareParameters>>,
    >,
    bond_compare_parameters: crate::arguments::Argument<
        Option<PyRef<'_, McsBondCompareParameters>>,
    >,
    atom_comparator: crate::arguments::Argument<Option<McsAtomComparatorInput>>,
    bond_comparator: crate::arguments::Argument<Option<McsBondComparatorInput>>,
    initial_seed: crate::arguments::Argument<Option<String>>,
) -> PyResult<McsResult> {
    let keyword_form = store_all.is_supplied()
        || maximize_bonds.is_supplied()
        || threshold.is_supplied()
        || timeout.is_supplied()
        || verbose.is_supplied()
        || atom_compare_parameters.is_supplied()
        || bond_compare_parameters.is_supplied()
        || atom_comparator.is_supplied()
        || bond_comparator.is_supplied()
        || initial_seed.is_supplied();
    crate::arguments::configuration_arguments!(params; store_all, maximize_bonds, threshold, timeout, verbose, atom_compare_parameters, bond_compare_parameters, atom_comparator, bond_comparator, initial_seed);
    let params = params.as_deref();
    let atom_compare_parameters = atom_compare_parameters.as_deref();
    let bond_compare_parameters = bond_compare_parameters.as_deref();
    if keyword_form {
        _maximum_common_substructure_keywords(
            py,
            inputs,
            store_all,
            maximize_bonds,
            threshold,
            timeout,
            verbose,
            atom_compare_parameters,
            bond_compare_parameters,
            atom_comparator,
            bond_comparator,
            initial_seed,
        )
    } else {
        _maximum_common_substructure_params(py, inputs, params)
    }
}

#[pyfunction]
#[doc = r#"{"role":"overload","target":"maximum_common_substructure"}"#]
#[pyo3(signature=(inputs, params=None))]
pub fn _maximum_common_substructure_params(
    py: Python<'_>,
    inputs: Vec<Py<Molecule>>,
    params: Option<&McsParameters>,
) -> PyResult<McsResult> {
    let params = params.map(|value| value.to_rust(py)).transpose()?;
    run_mcs(py, inputs, params.as_ref())
}

fn run_mcs(
    py: Python<'_>,
    inputs: Vec<Py<Molecule>>,
    params: Option<&ck::McsParameters>,
) -> PyResult<McsResult> {
    // Input handles and their immutable Python borrow guards live until the
    // synchronous Rust call has finished. Only Rust references cross detach;
    // no Python object access, graph cloning or target reconstruction occurs.
    let guards = inputs
        .iter()
        .map(|value| value.try_borrow(py))
        .collect::<Result<Vec<_>, _>>()?;
    let native: Vec<_> = guards.iter().map(|value| &value.inner).collect();
    let result = match params {
        Some(params) => py.detach(|| ck::maximum_common_substructure_with_params(&native, params)),
        None => py.detach(|| ck::maximum_common_substructure(&native)),
    };
    result
        .map(|inner| McsResult {
            inner: Arc::new(inner),
        })
        .map_err(McsError::from_rust)
}

#[pyfunction]
#[doc = r#"{"role":"overload","target":"maximum_common_substructure"}"#]
#[pyo3(signature=(inputs, *, store_all=None, maximize_bonds=None, threshold=None, timeout=None, verbose=None, atom_compare_parameters=None, bond_compare_parameters=None, atom_comparator=None, bond_comparator=None, initial_seed=None))]
pub fn _maximum_common_substructure_keywords(
    py: Python<'_>,
    inputs: Vec<Py<Molecule>>,
    store_all: Option<bool>,
    maximize_bonds: Option<bool>,
    threshold: Option<f64>,
    timeout: Option<u32>,
    verbose: Option<bool>,
    atom_compare_parameters: Option<&McsAtomCompareParameters>,
    bond_compare_parameters: Option<&McsBondCompareParameters>,
    atom_comparator: Option<McsAtomComparatorInput>,
    bond_comparator: Option<McsBondComparatorInput>,
    initial_seed: Option<String>,
) -> PyResult<McsResult> {
    let defaults = ck::McsParameters::default();
    let params = ck::McsParameters {
        store_all: store_all.unwrap_or(defaults.store_all),
        maximize_bonds: maximize_bonds.unwrap_or(defaults.maximize_bonds),
        threshold: threshold.unwrap_or(defaults.threshold),
        timeout: timeout.unwrap_or(defaults.timeout),
        verbose: verbose.unwrap_or(defaults.verbose),
        initial_seed: initial_seed.unwrap_or(defaults.initial_seed),
        atom_compare_parameters: atom_compare_parameters
            .map(|value| value.to_rust(py))
            .transpose()?
            .unwrap_or(defaults.atom_compare_parameters),
        bond_compare_parameters: bond_compare_parameters
            .map(|value| value.to_rust(py))
            .transpose()?
            .unwrap_or(defaults.bond_compare_parameters),
        atom_comparator: atom_comparator
            .map(McsAtomComparatorInput::into_rust)
            .transpose()?
            .unwrap_or(defaults.atom_comparator),
        bond_comparator: bond_comparator
            .map(McsBondComparatorInput::into_rust)
            .transpose()?
            .unwrap_or(defaults.bond_comparator),
    };
    run_mcs(py, inputs, Some(&params))
}

#[pyfunction]
pub fn maximum_common_substructure_with_params(
    py: Python<'_>,
    inputs: Vec<Py<Molecule>>,
    params: &McsParameters,
) -> PyResult<McsResult> {
    _maximum_common_substructure_params(py, inputs, Some(params))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query_compiler::QueryGraph;
    use pyo3::types::{PyDict, PyList, PyModule};
    use std::collections::BTreeMap;

    #[test]
    fn actual_mcs_producers_keep_options_partial_result_fields_sources_and_query_owners() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "mcs_probe")?;
            module.add_function(wrap_pyfunction!(maximum_common_substructure, &module)?)?;
            module.add_function(wrap_pyfunction!(
                maximum_common_substructure_with_params,
                &module
            )?)?;
            let first = Py::new(
                py,
                Molecule {
                    inner: ck::Molecule::from_smiles("CCO").unwrap(),
                },
            )?;
            let second = Py::new(
                py,
                Molecule {
                    inner: ck::Molecule::from_smiles("CCN").unwrap(),
                },
            )?;
            let before = [
                first.borrow(py).inner.to_binary().unwrap(),
                second.borrow(py).inner.to_binary().unwrap(),
            ];
            let inputs = PyList::new(py, [first.clone_ref(py), second.clone_ref(py)])?;
            let expected = ck::maximum_common_substructure(&[
                &first.borrow(py).inner,
                &second.borrow(py).inner,
            ])
            .unwrap();
            let baseline = module
                .getattr("maximum_common_substructure")?
                .call1((&inputs,))?
                .extract::<Py<McsResult>>()?;
            assert_eq!(baseline.borrow(py).inner.as_ref(), &expected);
            let mut queries = Vec::new();
            for native_params in [
                ck::McsParameters::default(),
                ck::McsParameters {
                    store_all: true,
                    maximize_bonds: false,
                    threshold: 0.5,
                    initial_seed: "C".into(),
                    ..Default::default()
                },
                ck::McsParameters {
                    atom_comparator: ck::McsAtomComparator::AtomCompareAny,
                    bond_comparator: ck::McsBondComparator::BondCompareAny,
                    atom_compare_parameters: ck::McsAtomCompareParameters {
                        match_formal_charge: true,
                        match_valences: true,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                ck::McsParameters {
                    bond_compare_parameters: ck::McsBondCompareParameters {
                        ring_matches_ring_only: true,
                        complete_rings_only: true,
                        ..Default::default()
                    },
                    ..Default::default()
                },
            ] {
                let expected = ck::maximum_common_substructure_with_params(
                    &[&first.borrow(py).inner, &second.borrow(py).inner],
                    &native_params,
                )
                .unwrap();
                let original_params = native_params.clone();
                let kwargs = PyDict::new(py);
                kwargs.set_item("store_all", native_params.store_all)?;
                kwargs.set_item("maximize_bonds", native_params.maximize_bonds)?;
                kwargs.set_item("threshold", native_params.threshold)?;
                kwargs.set_item("timeout", native_params.timeout)?;
                kwargs.set_item("verbose", native_params.verbose)?;
                kwargs.set_item("initial_seed", native_params.initial_seed.clone())?;
                kwargs.set_item(
                    "atom_compare_parameters",
                    Py::new(
                        py,
                        McsAtomCompareParameters {
                            inner: crate::parameters::NativeField::Owned(
                                native_params.atom_compare_parameters.clone(),
                            ),
                        },
                    )?,
                )?;
                kwargs.set_item(
                    "bond_compare_parameters",
                    Py::new(
                        py,
                        McsBondCompareParameters {
                            inner: crate::parameters::NativeField::Owned(
                                native_params.bond_compare_parameters.clone(),
                            ),
                        },
                    )?,
                )?;
                kwargs.set_item(
                    "atom_comparator",
                    crate::chemical_values::McsAtomComparator::from(native_params.atom_comparator),
                )?;
                kwargs.set_item(
                    "bond_comparator",
                    crate::chemical_values::McsBondComparator::from(native_params.bond_comparator),
                )?;
                let params = Py::new(py, McsParameters::from_rust(py, native_params)?)?;
                for actual in [
                    module
                        .getattr("maximum_common_substructure")?
                        .call1((&inputs, params.clone_ref(py)))?,
                    module
                        .getattr("maximum_common_substructure")?
                        .call((&inputs,), Some(&kwargs))?,
                    module
                        .getattr("maximum_common_substructure_with_params")?
                        .call1((&inputs, params.clone_ref(py)))?,
                ] {
                    let actual = actual.extract::<Py<McsResult>>()?;
                    assert_eq!(actual.borrow(py).inner.as_ref(), &expected);
                    if let Some(native_query) = expected.query.as_ref() {
                        let query = actual
                            .bind(py)
                            .getattr("query")?
                            .extract::<Py<QueryGraph>>()?;
                        assert_eq!(query.borrow(py).inner(), native_query);
                        queries.push((query, native_query.num_atoms()));
                    } else {
                        assert!(actual.bind(py).getattr("query")?.is_none());
                    }
                    let degenerate = actual.bind(py).getattr("degenerate")?.extract::<BTreeMap<
                        String,
                        Py<QueryGraph>,
                    >>(
                    )?;
                    assert_eq!(degenerate.len(), expected.degenerate.len());
                    for (key, query) in degenerate {
                        let native = &expected.degenerate[key.as_bytes()];
                        assert_eq!(query.borrow(py).inner(), native);
                        queries.push((query, native.num_atoms()));
                    }
                }
                assert!(
                    module
                        .getattr("maximum_common_substructure")?
                        .call((&inputs, params.clone_ref(py)), Some(&kwargs))
                        .unwrap_err()
                        .is_instance_of::<PyTypeError>(py)
                );
                assert_eq!(params.borrow(py).to_rust(py)?, original_params);
            }
            for values in [vec![], vec![first.clone_ref(py)]] {
                let count = values.len();
                let error = module
                    .getattr("maximum_common_substructure")?
                    .call1((values,))
                    .unwrap_err();
                assert!(
                    error.is_instance(
                        py,
                        <McsError as crate::error_framework::ExceptionPayload>::exception_type(py)
                            .expect("error backend")
                            .as_any()
                    )
                );
                assert_eq!(
                    error
                        .value(py)
                        .getattr("state")?
                        .getattr("count")?
                        .extract::<usize>()?,
                    count
                );
                assert!(error.cause(py).is_none());
            }
            let invalid = PyDict::new(py);
            invalid.set_item("threshold", 1.1)?;
            let expected_error = ck::maximum_common_substructure_with_params(
                &[&first.borrow(py).inner, &second.borrow(py).inner],
                &ck::McsParameters {
                    threshold: 1.1,
                    ..Default::default()
                },
            )
            .unwrap_err();
            let error = module
                .getattr("maximum_common_substructure")?
                .call((&inputs,), Some(&invalid))
                .unwrap_err();
            assert!(
                error.is_instance(
                    py,
                    <McsError as crate::error_framework::ExceptionPayload>::exception_type(py)
                        .expect("error backend")
                        .as_any()
                )
            );
            assert_eq!(error.value(py).str()?.to_str()?, expected_error.to_string());
            assert_eq!(first.borrow(py).inner.to_binary().unwrap(), before[0]);
            assert_eq!(second.borrow(py).inner.to_binary().unwrap(), before[1]);
            drop(inputs);
            drop(first);
            drop(second);
            drop(baseline);
            for (query, count) in queries {
                assert_eq!(
                    query
                        .bind(py)
                        .call_method0("num_atoms")?
                        .extract::<usize>()?,
                    count
                );
            }
            Ok(())
        })
        .unwrap();
    }
}
