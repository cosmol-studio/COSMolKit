//! Path-score failure representation; packing and input checks stay native.
use crate::cip_atom_code_errors::AtomCodeError;
use crate::error_framework::{Source, error_projection};
use crate::errors::FingerprintError;
use crate::model_errors::TopologyValidationError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(TopologicalTorsionPathScoreError, TopologicalTorsionPathScoreErrorKind, ck::TopologicalTorsionPathScoreError, "fingerprint", {
    ZeroSize, ShortPath { actual, required },
    ShortAtomCodes { actual, required },
    AtomIndexOutOfRange { index, atom_count: usize },
    AtomCodeUnderflow { index, code: u32, subtract: u32 },
    InvalidTopology(_), AtomCode(_), PackedCode(_),
}, cause = |inner: &ck::TopologicalTorsionPathScoreError| match inner {
    ck::TopologicalTorsionPathScoreError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::TopologicalTorsionPathScoreError::AtomCode(value) => Some(AtomCodeError::from_rust(value.clone())),
    ck::TopologicalTorsionPathScoreError::PackedCode(value) => Some(FingerprintError::from_rust(value.clone())),
    ck::TopologicalTorsionPathScoreError::ZeroSize
    | ck::TopologicalTorsionPathScoreError::ShortPath { .. }
    | ck::TopologicalTorsionPathScoreError::ShortAtomCodes { .. }
    | ck::TopologicalTorsionPathScoreError::AtomIndexOutOfRange { .. }
    | ck::TopologicalTorsionPathScoreError::AtomCodeUnderflow { .. } => None,
}, source = Option<Source<Py<TopologyValidationError>, Source<Py<AtomCodeError>, Py<FingerprintError>>>>, methods = {
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::TopologicalTorsionPathScoreError::ShortPath { actual, .. }
        | ck::TopologicalTorsionPathScoreError::ShortAtomCodes { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn required(&self) -> Option<usize> { match self.inner {
        ck::TopologicalTorsionPathScoreError::ShortPath { required, .. }
        | ck::TopologicalTorsionPathScoreError::ShortAtomCodes { required, .. } => Some(required), _ => None,
    } }
    #[getter]
    fn index(&self) -> Option<usize> { match self.inner {
        ck::TopologicalTorsionPathScoreError::AtomIndexOutOfRange { index, .. }
        | ck::TopologicalTorsionPathScoreError::AtomCodeUnderflow { index, .. } => Some(index), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_path_scores_keep_native_order_optional_codes_and_exact_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CCCC").unwrap();
            let before = native.to_binary().unwrap();
            let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
            for path in [
                vec![0, 1, 2, 3],
                vec![3, 2, 1, 0],
                vec![0, 1, 2, 3, usize::MAX],
            ] {
                for codes in [
                    None,
                    Some(vec![]),
                    Some(vec![20, 40, 60, 80]),
                    Some(vec![u32::MAX; 4]),
                ] {
                    let expected = molecule
                        .borrow(py)
                        .inner
                        .topological_torsion_path_score(&path, 4, codes.as_deref())
                        .unwrap();
                    let actual = molecule
                        .bind(py)
                        .call_method1("topological_torsion_path_score", (path.clone(), 4, codes))?
                        .extract::<u64>()?;
                    assert_eq!(actual, expected);
                }
            }
            let expected = molecule
                .borrow(py)
                .inner
                .topological_torsion_path_score(&[0, 1, 2, 3], 4, None)
                .unwrap();
            assert_eq!(
                molecule
                    .bind(py)
                    .call_method1("topological_torsion_path_score", (vec![0, 1, 2, 3], 4))?
                    .extract::<u64>()?,
                expected
            );
            for (path, size, codes) in [
                (vec![], 0, None),
                (vec![0], 4, None),
                (vec![0, 1, 2, 3], 4, Some(vec![5])),
                (vec![usize::MAX, 1, 2, 3], 4, Some(vec![5; 4])),
                (vec![0, 1, 2, 3], 4, Some(vec![0; 4])),
            ] {
                let native = molecule
                    .borrow(py)
                    .inner
                    .topological_torsion_path_score(&path, size, codes.as_deref())
                    .unwrap_err();
                let expected = native.to_string();
                let error = molecule
                    .bind(py)
                    .call_method1("topological_torsion_path_score", (path, size, codes))
                    .unwrap_err();
                assert!(error.is_instance_of::<TopologicalTorsionPathScoreError>(py));
                assert_eq!(error.value(py).str()?.to_str()?, expected);
                assert!(error.value(py).getattr("source")?.is_none());
                if let ck::TopologicalTorsionPathScoreError::AtomIndexOutOfRange { index, .. } =
                    native
                {
                    assert_eq!(error.value(py).getattr("index")?.extract::<usize>()?, index);
                    assert!(error.value(py).setattr("index", 0).is_err());
                }
            }
            let error = TopologicalTorsionPathScoreError::from_rust(
                ck::TopologicalTorsionPathScoreError::AtomCode(ck::AtomCodeError::Precondition {
                    what: "real native category",
                }),
            );
            let source = error.value(py).getattr("source")?;
            assert!(source.is_instance_of::<AtomCodeError>());
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(
                source.getattr("what")?.extract::<String>()?,
                "real native category"
            );
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            Ok(())
        })
        .unwrap();
    }
}
