//! Public producers use the existing operation runtime and scoped callbacks.
use crate::native::Molecule;
use crate::operation_errors::OperationError;
use crate::tautomer_configuration::TautomerParams;
use crate::tautomer_values::TautomerEnumeration;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::sync::Arc;

pub(crate) fn keywords(
    py: Python<'_>,
    max_tautomers: Option<u32>,
    max_transforms: Option<u32>,
    remove_sp3_stereo: Option<bool>,
    remove_bond_stereo: Option<bool>,
    remove_isotopic_hydrogens: Option<bool>,
    reassign_stereo: Option<bool>,
    callback: Option<crate::tautomer_callbacks::EnumerationCallback>,
    scorer: Option<crate::tautomer_callbacks::Scorer>,
    score_params: Option<Py<crate::tautomer_values::TautomerScoreParams>>,
) -> TautomerParams {
    let defaults = ck::TautomerParams::default();
    let score_params = score_params.as_ref().map(|value| value.borrow(py));
    TautomerParams::new(
        max_tautomers.unwrap_or_else(|| defaults.max_tautomers()),
        max_transforms.unwrap_or_else(|| defaults.max_transforms()),
        remove_sp3_stereo.unwrap_or_else(|| defaults.remove_sp3_stereo()),
        remove_bond_stereo.unwrap_or_else(|| defaults.remove_bond_stereo()),
        remove_isotopic_hydrogens.unwrap_or_else(|| defaults.remove_isotopic_hydrogens()),
        reassign_stereo.unwrap_or_else(|| defaults.reassign_stereo()),
        callback,
        scorer,
        score_params.as_deref(),
    )
}

pub(crate) fn enumerate(
    py: Python<'_>,
    molecule: &Molecule,
    params: Option<&TautomerParams>,
) -> PyResult<TautomerEnumeration> {
    let inner = match params {
        None => py.detach(|| molecule.inner.enumerate_tautomers()),
        Some(params) => params.with_callbacks(py, |params| {
            py.detach(|| molecule.inner.enumerate_tautomers_with_params(params))
        })?,
    }
    .map_err(OperationError::from_rust)?;
    Ok(TautomerEnumeration {
        inner: Arc::new(inner),
    })
}
pub(crate) fn canonical(
    py: Python<'_>,
    molecule: &Molecule,
    params: Option<&TautomerParams>,
) -> PyResult<Molecule> {
    let inner = match params {
        None => py.detach(|| molecule.inner.canonical_tautomer()),
        Some(params) => params.with_callbacks(py, |params| {
            py.detach(|| molecule.inner.canonical_tautomer_with_params(params))
        })?,
    }
    .map_err(OperationError::from_rust)?;
    Ok(Molecule { inner })
}
pub(crate) fn canonical_enumeration(
    py: Python<'_>,
    result: &TautomerEnumeration,
    params: Option<&TautomerParams>,
) -> PyResult<Molecule> {
    let inner = match params {
        None => py.detach(|| result.inner.canonical_tautomer()),
        Some(params) => params.with_callbacks(py, |params| {
            py.detach(|| result.inner.canonical_tautomer_with_params(params))
        })?,
    }
    .map_err(OperationError::from_rust)?;
    Ok(Molecule { inner })
}

#[pyfunction]
#[pyo3(signature=(molecules, params=None))]
pub fn canonical_tautomer_from_molecules(
    py: Python<'_>,
    molecules: Vec<Py<Molecule>>,
    params: Option<Py<TautomerParams>>,
) -> PyResult<Molecule> {
    let params = params.map(|value| value.borrow(py).clone_for_call(py));
    let guards = molecules
        .iter()
        .map(|value| value.try_borrow(py))
        .collect::<Result<Vec<_>, _>>()
        .map_err(PyErr::from)?;
    let native = guards.iter().map(|value| &value.inner).collect::<Vec<_>>();
    let invoke = |params: &ck::TautomerParams| {
        py.detach(|| {
            ck::canonical_tautomer_from_molecule_hosts_with_params(
                native.len(),
                |index, run| run(native[index]),
                params,
            )
        })
    };
    let inner = match params.as_ref() {
        Some(params) => params.with_callbacks(py, invoke)?,
        None => invoke(&ck::TautomerParams::default()),
    }
    .map_err(OperationError::from_rust)?;
    Ok(Molecule { inner })
}
#[pyfunction]
pub fn canonical_tautomer_from_molecules_with_params(
    py: Python<'_>,
    molecules: Vec<Py<Molecule>>,
    params: Py<TautomerParams>,
) -> PyResult<Molecule> {
    canonical_tautomer_from_molecules(py, molecules, Some(params))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::operation_errors::{OperationErrorKind, SmilesError};
    use pyo3::ffi::c_str;
    use pyo3::types::{PyDict, PyModule};

    #[test]
    fn public_tautomer_producers_preserve_native_results_scoped_callbacks_and_operation_payloads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CC(=O)C").unwrap();
            let before = native.to_binary().unwrap();
            let expected = native.enumerate_tautomers().unwrap();
            let expected_text = crate::property_values::decode_source_text(py, &expected.canonical_tautomer().unwrap().to_smiles().unwrap())?;
            let molecule = Py::new(py, Molecule { inner: native })?;
            let locals = PyDict::new(py);
            locals.set_item("molecule", molecule.bind(py))?;
            locals.set_item("expected_count", expected.len())?;
            locals.set_item("expected_text", expected_text)?;
            locals.set_item("before", before.clone())?;
            let module = PyModule::new(py, "tautomer_producer_test")?;
            module.add_class::<TautomerParams>()?;
            module.add_function(wrap_pyfunction!(canonical_tautomer_from_molecules, &module)?)?;
            module.add_function(wrap_pyfunction!(canonical_tautomer_from_molecules_with_params, &module)?)?;
            locals.set_item("api", module)?;
            py.run(c_str!(r#"
import gc
params = api.TautomerParams(max_tautomers=23)
default = molecule.enumerate_tautomers()
assert default.len() == expected_count
configured = molecule.enumerate_tautomers(params)
assert configured.canonical_smiles() == default.canonical_smiles()
assert molecule.enumerate_tautomers(max_tautomers=23).canonical_smiles() == configured.canonical_smiles()
assert molecule.enumerate_tautomers_with_params(params).canonical_smiles() == configured.canonical_smiles()
for result in [molecule.canonical_tautomer(), molecule.canonical_tautomer(params),
    molecule.canonical_tautomer(max_tautomers=23), molecule.canonical_tautomer_with_params(params),
    default.canonical_tautomer(), default.canonical_tautomer(params),
    default.canonical_tautomer(max_tautomers=23), default.canonical_tautomer_with_params(params),
    api.canonical_tautomer_from_molecules(list(default.iter())),
    api.canonical_tautomer_from_molecules(list(default.iter()), params),
    api.canonical_tautomer_from_molecules_with_params(list(default.iter()), params)]:
    assert result.to_smiles() == expected_text
for call in [lambda: molecule.enumerate_tautomers(params, max_tautomers=23),
             lambda: molecule.canonical_tautomer(params, max_tautomers=23),
             lambda: default.canonical_tautomer(params, max_tautomers=23)]:
    try: call()
    except TypeError: pass
    else: raise AssertionError('mixed configuration accepted')
live = []
snapshots = []
scores = []
def callback(source, progress):
    params.max_tautomers = 23  # the producer must not hold a PyClass borrow
    live.append(source)
    snapshots.append(progress.snapshot())
    return True
def scorer(source):
    scores.append(source)
    return source.tautomer_score().total()
params.callback = callback
params.scorer = scorer
retained = molecule.enumerate_tautomers(params)
assert live and snapshots
chosen = retained.canonical_tautomer(params)
assert scores and chosen.to_smiles() == expected_text
assert molecule.canonical_tautomer(params).to_smiles() == expected_text
for source in live + scores:
    try: source.num_atoms()
    except ReferenceError: pass
    else: raise AssertionError('callback loan escaped')
assert any(snapshot.len() for snapshot in snapshots)
marker = LookupError('same invocation exception')
def bad(source, progress): raise marker
params.callback = bad
try: molecule.enumerate_tautomers(params)
except LookupError as error: assert error is marker
else: raise AssertionError('Python exception lost')
params.callback = None
params.scorer = None
assert molecule.enumerate_tautomers(params).len() == expected_count
assert bytes(molecule.to_binary()) == bytes(before)
del molecule, params, configured, default
gc.collect()
assert retained.canonical_tautomer().to_smiles() == expected_text
assert chosen.to_smiles() == expected_text
"#), Some(&locals), None)?;
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);

            let spec = ck::operation_specs()[0];
            let unsupported = OperationError::from_rust(ck::OperationError::UnsupportedFeature {
                operation: spec,
                source: ck::UnsupportedFeatureError { feature: "exact\0feature", reason: "reason" },
            });
            let value = unsupported.value(py);
            assert!(unsupported.cause(py).is_none());
            assert!(value.getattr("cause")?.is_none());
            assert_eq!(value.getattr("source")?.getattr("feature")?.extract::<String>()?, "exact\0feature");
            assert_eq!(value.getattr("operation")?.getattr("method")?.extract::<String>()?, spec.method);
            assert!(value.getattr("source")?.setattr("reason", "changed").is_err());
            let mapping = ck::MappingValidationError::Length { entity: "atom", direction: "source", actual: 2, expected: 3 };
            let error = OperationError::from_rust(ck::OperationError::InvalidTopologyMapping { operation: "actual", source: mapping });
            assert!(error.cause(py).is_none());
            assert_eq!(error.value(py).getattr("source")?.getattr("actual")?.extract::<usize>()?, 2);
            let derived = OperationError::from_rust(ck::OperationError::DerivedEffectContract { operation: spec.method, action: "preserve", states: spec.derived_effects.preserve, issue: "exact" });
            assert_eq!(derived.value(py).getattr("states")?.getattr("bits")?.extract::<u16>()?, spec.derived_effects.preserve.bits());
            let missing = OperationError::from_rust(ck::OperationError::SemanticPreconditionContract { operation: spec.method, missing: spec.semantic_preconditions, issue: "exact" });
            assert_eq!(missing.value(py).getattr("missing")?.getattr("bits")?.extract::<u16>()?, u16::from(spec.semantic_preconditions.bits()));
            let policy = OperationError::from_rust(ck::OperationError::CipStateContract { operation: spec.method, policy: spec.cip_state, issue: "exact" });
            assert_eq!(policy.value(py).getattr("policy")?.getattr("discriminant")?.extract::<isize>()?, spec.cip_state as isize);
            let requirement = OperationError::from_rust(ck::OperationError::MappingContract { operation: spec.method, requirement: spec.requires_mapping, issue: "exact" });
            assert_eq!(requirement.value(py).getattr("requirement")?.getattr("discriminant")?.extract::<isize>()?, spec.requires_mapping as isize);
            let empty = ck::TautomerEnumeration::default().canonical_tautomer().unwrap_err();
            let text = empty.to_string();
            let error = OperationError::from_rust(empty);
            assert_eq!(error.value(py).str()?.to_str()?, text);
            assert!(error.value(py).getattr("variant")?.extract::<OperationErrorKind>()? == OperationErrorKind::Tautomer);
            assert!(error.value(py).getattr("cause")?.is(&error.value(py).getattr("__cause__")?));
            let parse = SmilesError::from_rust(ck::Molecule::from_smiles("CC(").unwrap_err());
            assert!(parse.cause(py).unwrap().is_instance_of::<crate::smiles_parse_errors::SmilesParseError>(py));
            let operation = ck::OperationError::Algorithm { operation: "exact", detail: "detail".into() };
            let construction = SmilesError::from_rust(ck::SmilesError::Construction(operation));
            assert!(construction.cause(py).unwrap().is_instance_of::<OperationError>(py));
            Ok(())
        }).unwrap();
    }
}
