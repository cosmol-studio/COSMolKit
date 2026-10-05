//! Canonical descriptor errors projected only through the public facade.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pyo3::create_exception!(cosmolkit, DescriptorReadError, PyValueError);
pyo3::create_exception!(cosmolkit, DescriptorError, PyValueError);

fn domain_pyerr(py: Python<'_>, source: &ck::DescriptorError) -> PyErr {
    use ck::DescriptorError as E;
    let error = DescriptorError::new_err(source.to_string());
    let attrs = || -> PyResult<()> {
        let value = error.value(py);
        value.setattr("domain", "descriptors")?;
        let kind = match source {
            E::InvalidConnectivityPath {
                function,
                expected_rows,
                actual_rows,
            } => {
                value.setattr("function", *function)?;
                value.setattr("expected_rows", *expected_rows)?;
                value.setattr("actual_rows", *actual_rows)?;
                "InvalidConnectivityPath"
            }
            E::InvalidTopology { function, .. } => {
                value.setattr("function", *function)?;
                "InvalidTopology"
            }
            E::Path { function, .. } => {
                value.setattr("function", *function)?;
                "Path"
            }
            E::InvalidHallKierContributionRows { actual, minimum } => {
                value.setattr("actual", *actual)?;
                value.setattr("minimum", *minimum)?;
                "InvalidHallKierContributionRows"
            }
            E::MissingComputedScalar {
                function,
                include_sulfur_phosphorus,
            } => {
                value.setattr("function", *function)?;
                value.setattr("include_sulfur_phosphorus", *include_sulfur_phosphorus)?;
                "MissingComputedScalar"
            }
            E::MissingLabuteHydrogens { function } => {
                value.setattr("function", *function)?;
                "MissingLabuteHydrogens"
            }
            E::MissingLabuteAsa { function } => {
                value.setattr("function", *function)?;
                "MissingLabuteAsa"
            }
            E::MismatchedBinArrays {
                contribs_len,
                bin_prop_len,
                bins_len,
            } => {
                value.setattr("contribs_len", *contribs_len)?;
                value.setattr("bin_prop_len", *bin_prop_len)?;
                value.setattr("bins_len", *bins_len)?;
                "MismatchedBinArrays"
            }
            E::CrippenParamNumeric { field, cell } => {
                value.setattr("field", *field)?;
                value.setattr("cell", cell.as_str())?;
                "CrippenParamNumeric"
            }
            E::MissingCrippenMrContributions { function } => {
                value.setattr("function", *function)?;
                "MissingCrippenMrContributions"
            }
            E::MissingCrippenMr { function } => {
                value.setattr("function", *function)?;
                "MissingCrippenMr"
            }
            E::InvalidCrippenOptionalRows {
                field,
                actual,
                expected,
            } => {
                value.setattr("field", *field)?;
                value.setattr("actual", *actual)?;
                value.setattr("expected", *expected)?;
                "InvalidCrippenOptionalRows"
            }
            E::MissingCrippenDefaultPattern { row } => {
                value.setattr("row", *row)?;
                "MissingCrippenDefaultPattern"
            }
            E::TopologyEdit { function, .. } => {
                value.setattr("function", *function)?;
                "TopologyEdit"
            }
            E::MissingFinalHydrogenState { field } => {
                value.setattr("field", *field)?;
                "MissingFinalHydrogenState"
            }
            E::Hydrogens { function, .. } => {
                value.setattr("function", *function)?;
                "Hydrogens"
            }
            E::CountOverflow { function, field } => {
                value.setattr("function", *function)?;
                value.setattr("field", *field)?;
                "CountOverflow"
            }
            E::Valence { function, .. } => {
                value.setattr("function", *function)?;
                "Valence"
            }
            E::InvalidValenceRows {
                function,
                field,
                actual,
                expected,
            } => {
                value.setattr("function", *function)?;
                value.setattr("field", *field)?;
                value.setattr("actual", *actual)?;
                value.setattr("expected", *expected)?;
                "InvalidValenceRows"
            }
            E::Unsupported { function, detail } => {
                value.setattr("function", *function)?;
                value.setattr("detail", detail.as_str())?;
                "Unsupported"
            }
            E::Ring { function, .. } => {
                value.setattr("function", *function)?;
                "Ring"
            }
            E::Search { function, .. } => {
                value.setattr("function", *function)?;
                "Search"
            }
            E::Stereo { function, .. } => {
                value.setattr("function", *function)?;
                "Stereo"
            }
        };
        value.setattr("kind", kind)?;
        Ok(())
    };
    if let Err(attribute_error) = attrs() {
        return attribute_error;
    }
    // Preserve the real borrowed source chain; deeper Rust downcast identity
    // has no Python equivalent and is not invented from message parsing.
    error.set_cause(
        py,
        std::error::Error::source(source)
            .map(|cause| crate::canonical_values::source_pyerr(py, cause)),
    );
    error
}

pub(crate) fn descriptor_pyerr(py: Python<'_>, source: ck::DescriptorReadError) -> PyErr {
    use ck::DescriptorReadError as E;
    let error = DescriptorReadError::new_err(source.to_string());
    let kind = match &source {
        E::MissingPreparedValence => "MissingPreparedValence",
        E::CachePoisoned => "CachePoisoned",
        E::MissingInitializedRings => "MissingInitializedRings",
        E::Algorithm { .. } => "Algorithm",
    };
    if let Err(attribute_error) = error
        .value(py)
        .setattr("domain", "descriptors")
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return attribute_error;
    }
    if let E::Algorithm { source } = &source {
        error.set_cause(py, Some(domain_pyerr(py, source)));
    }
    error
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "DescriptorReadError",
        module.py().get_type::<DescriptorReadError>(),
    )?;
    module.add("DescriptorError", module.py().get_type::<DescriptorError>())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d02_missing_state_conversion_has_no_fabricated_child() {
        Python::attach(|py| {
            for source in [
                ck::DescriptorReadError::MissingPreparedValence,
                ck::DescriptorReadError::MissingInitializedRings,
            ] {
                let expected_kind = match &source {
                    ck::DescriptorReadError::MissingPreparedValence => "MissingPreparedValence",
                    _ => "MissingInitializedRings",
                };
                let error = descriptor_pyerr(py, source);
                assert!(error.is_instance_of::<DescriptorReadError>(py));
                assert_eq!(
                    error
                        .value(py)
                        .getattr("kind")
                        .unwrap()
                        .extract::<String>()
                        .unwrap(),
                    expected_kind
                );
                assert!(error.cause(py).is_none());
            }
        });
    }

    #[test]
    fn d02_algorithm_conversion_preserves_typed_child_context() {
        Python::attach(|py| {
            let error = descriptor_pyerr(
                py,
                ck::DescriptorReadError::Algorithm {
                    source: ck::DescriptorError::InvalidHallKierContributionRows {
                        actual: 1,
                        minimum: 2,
                    },
                },
            );
            assert!(error.is_instance_of::<DescriptorReadError>(py));
            assert_eq!(
                error
                    .value(py)
                    .getattr("kind")
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "Algorithm"
            );
            let child = error.cause(py).unwrap();
            assert!(child.is_instance_of::<DescriptorError>(py));
            assert_eq!(
                child
                    .value(py)
                    .getattr("kind")
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "InvalidHallKierContributionRows"
            );
            assert_eq!(
                child
                    .value(py)
                    .getattr("actual")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                1
            );
            assert_eq!(
                child
                    .value(py)
                    .getattr("minimum")
                    .unwrap()
                    .extract::<usize>()
                    .unwrap(),
                2
            );
            assert!(child.cause(py).is_none());
        });
    }

    #[test]
    fn d02_source_unsupported_kind_is_not_inferred_from_other_errors() {
        Python::attach(|py| {
            let error = domain_pyerr(
                py,
                &ck::DescriptorError::Unsupported {
                    function: "fixture",
                    detail: "independent source capability".to_owned(),
                },
            );
            assert!(error.is_instance_of::<DescriptorError>(py));
            assert_eq!(
                error
                    .value(py)
                    .getattr("kind")
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "Unsupported"
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("detail")
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "independent source capability"
            );
        });
    }
}
