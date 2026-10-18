use crate::chemical_values::{BondOrder, enum_value};
use crate::error_framework::error_projection;
use crate::model_errors::TopologyValidationError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

enum_value!(ValencePhase, { EffectiveAtomicNumber, Explicit, Implicit });

error_projection!(ValenceError, ValenceErrorKind, ck::ValenceError, "valence", {
    PiElectronExplicitValenceCacheNotInitialized { atom },
    PiElectronInvariant { atom, explicit_valence: u32, physical_bonds: u32 },
    InvalidValence { atom, atomic_number, formal_charge: i8, phase, calculated: Option<i32>, reason: &'static str, message: String },
    InvalidTopology { source }, AtomOutOfRange { atom, atom_count },
    AdjacencyAtomOutOfRange { atom, neighbor_atom, atom_count },
    AdjacencyBondOutOfRange { atom, bond, bond_count: usize },
    AdjacencyEndpointMismatch { atom, neighbor_atom, bond, begin: usize => |value: &ck::AtomId| value.index(), end: usize => |value: &ck::AtomId| value.index() },
    InvalidExplicitValenceInput { atom, value: i32 },
    PeriodicTableLookup { atomic_number, field: &'static str },
    ExplicitValenceCacheNotInitialized { atom }, ImplicitValenceCacheNotInitialized { atom },
    HydrogenCountOverflow { atom, explicit: u32, implicit: u32, neighbor_hydrogens: usize },
    BadBondType { bond, order: BondOrder => |value: &ck::BondOrder| (*value).into() },
}, cause = |inner: &ck::ValenceError| match inner {
    ck::ValenceError::InvalidTopology { source } => Some(TopologyValidationError::from_rust(source.clone())),
    ck::ValenceError::PiElectronExplicitValenceCacheNotInitialized { .. }
    | ck::ValenceError::PiElectronInvariant { .. } | ck::ValenceError::InvalidValence { .. }
    | ck::ValenceError::AtomOutOfRange { .. } | ck::ValenceError::AdjacencyAtomOutOfRange { .. }
    | ck::ValenceError::AdjacencyBondOutOfRange { .. } | ck::ValenceError::AdjacencyEndpointMismatch { .. }
    | ck::ValenceError::InvalidExplicitValenceInput { .. } | ck::ValenceError::PeriodicTableLookup { .. }
    | ck::ValenceError::ExplicitValenceCacheNotInitialized { .. } | ck::ValenceError::ImplicitValenceCacheNotInitialized { .. }
    | ck::ValenceError::HydrogenCountOverflow { .. } | ck::ValenceError::BadBondType { .. } => None,
}, source = Option<Py<TopologyValidationError>>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::ValenceError::PiElectronExplicitValenceCacheNotInitialized { atom }
        | ck::ValenceError::PiElectronInvariant { atom, .. } | ck::ValenceError::InvalidValence { atom, .. }
        | ck::ValenceError::AtomOutOfRange { atom, .. } | ck::ValenceError::AdjacencyAtomOutOfRange { atom, .. }
        | ck::ValenceError::AdjacencyBondOutOfRange { atom, .. } | ck::ValenceError::AdjacencyEndpointMismatch { atom, .. }
        | ck::ValenceError::InvalidExplicitValenceInput { atom, .. } | ck::ValenceError::ExplicitValenceCacheNotInitialized { atom }
        | ck::ValenceError::ImplicitValenceCacheNotInitialized { atom } | ck::ValenceError::HydrogenCountOverflow { atom, .. } => Some(atom.index()),
        _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn phase(&self) -> Option<ValencePhase> {
        match self.inner { ck::ValenceError::InvalidValence { phase, .. } => Some(phase.into()), _ => None }
    }
    #[getter]
    fn atomic_number(&self) -> Option<u8> { match self.inner {
        ck::ValenceError::InvalidValence { atomic_number, .. }
        | ck::ValenceError::PeriodicTableLookup { atomic_number, .. } => Some(atomic_number), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::ValenceError::AtomOutOfRange { atom_count, .. }
        | ck::ValenceError::AdjacencyAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn neighbor_atom(&self) -> Option<usize> { match self.inner {
        ck::ValenceError::AdjacencyAtomOutOfRange { neighbor_atom, .. }
        | ck::ValenceError::AdjacencyEndpointMismatch { neighbor_atom, .. } => Some(neighbor_atom), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::ValenceError::AdjacencyBondOutOfRange { bond, .. }
        | ck::ValenceError::AdjacencyEndpointMismatch { bond, .. } => Some(bond.index()),
        ck::ValenceError::BadBondType { bond, .. } => bond.map(ck::BondId::index), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;

    #[test]
    fn real_valence_queries_keep_precise_phase_payload_and_source() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let invalid = ck::Molecule::from_smiles_with_params(
                "[CH5]",
                &ck::SmilesParseParams {
                    sanitize: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let native = invalid.atom_metadata(true).unwrap_err();
            let display = native.to_string();
            let ck::ValenceError::InvalidValence {
                phase,
                calculated,
                atom,
                ..
            } = &native
            else {
                panic!("expected native invalid valence")
            };
            let expected_phase: ValencePhase = (*phase).into();
            let expected_calculated = *calculated;
            let expected_atom = atom.index();
            assert!(native.source().is_none());
            let molecule = Py::new(py, Molecule { inner: invalid })?;
            let error = molecule.bind(py).call_method0("atom_metadata").unwrap_err();
            assert!(error.is_instance_of::<ValenceError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(
                error
                    .value(py)
                    .getattr("phase")?
                    .extract::<ValencePhase>()?
                    == expected_phase
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("calculated")?
                    .extract::<Option<i32>>()?,
                expected_calculated
            );
            assert_eq!(
                error.value(py).getattr("atom")?.extract::<usize>()?,
                expected_atom
            );
            assert!(error.cause(py).is_none());
            assert!(
                error
                    .value(py)
                    .setattr("phase", expected_phase)
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );
            assert!(
                molecule
                    .bind(py)
                    .call_method1("has_valence_violation", (0,))?
                    .extract::<bool>()?
            );
            let native = ck::Molecule::from_smiles("CCO").unwrap();
            let expected = native.atom_metadata(true).unwrap();
            let projected = Py::new(py, Molecule { inner: native })?;
            let actual = projected.bind(py).call_method0("atom_metadata")?;
            for (index, row) in expected.iter().enumerate() {
                let item = actual.get_item(index)?;
                for (field, value) in [
                    ("explicit_valence", row.explicit_valence),
                    ("implicit_hydrogens", row.implicit_hydrogens),
                    ("total_hydrogens", row.total_hydrogens),
                    ("total_valence", row.total_valence),
                ] {
                    assert_eq!(item.getattr(field)?.extract::<i32>()?, value);
                }
                assert_eq!(item.getattr("degree")?.extract::<usize>()?, row.degree);
            }
            let error = ValenceError::from_rust(ck::ValenceError::InvalidTopology {
                source: ck::TopologyValidationError::AdjacencyMismatch,
            });
            let cause = error.cause(py).unwrap();
            assert!(cause.is_instance_of::<TopologyValidationError>(py));
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            assert!(error.value(py).getattr("phase")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
