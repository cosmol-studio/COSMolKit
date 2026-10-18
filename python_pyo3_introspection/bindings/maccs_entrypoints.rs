//! The registered raw vector and explicit Params entrypoints use the facade.
use crate::native::{Fingerprint, Molecule};
use crate::params::MaccsFingerprintParams;
use cosmolkit as ck;
use pyo3::prelude::*;

#[test]
fn raw_and_explicit_params_preserve_complete_native_vectors_and_errors() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        for text in ["CCO", "c1ccccc1", "[Na+].[Cl-]"] {
            let native = ck::Molecule::from_smiles(text).unwrap();
            let before = native.to_binary().unwrap();
            let raw = native.fingerprint_maccs_raw().unwrap();
            let expected = native.fingerprint_maccs().unwrap();
            let molecule = Py::new(py, Molecule { inner: native })?;
            let params = Py::new(
                py,
                MaccsFingerprintParams {
                    inner: ck::MaccsFingerprintParams::default(),
                },
            )?;
            let projected_raw = molecule
                .bind(py)
                .call_method0("fingerprint_maccs_raw")?
                .extract::<Py<Fingerprint>>()?;
            assert_eq!(projected_raw.borrow(py).inner, raw);
            assert_eq!(projected_raw.borrow(py).inner.n_bits(), 167);
            assert!(!projected_raw.borrow(py).inner.get_bit(0).unwrap());
            for result in [
                molecule.bind(py).call_method0("fingerprint_maccs")?,
                molecule
                    .bind(py)
                    .call_method1("fingerprint_maccs", (params.clone_ref(py),))?,
                molecule
                    .bind(py)
                    .call_method1("fingerprint_maccs_with_params", (params.clone_ref(py),))?,
            ] {
                assert_eq!(result.extract::<PyRef<'_, Fingerprint>>()?.inner, expected);
            }
            params.bind(py).setattr("n_bits", 200)?;
            let expected = molecule
                .borrow(py)
                .inner
                .fingerprint_maccs_with_params(&params.borrow(py).inner)
                .unwrap_err();
            let error = molecule
                .bind(py)
                .call_method1("fingerprint_maccs_with_params", (params,))
                .unwrap_err();
            assert!(error.is_instance_of::<crate::maccs_errors::MaccsFingerprintError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
        }
        Ok(())
    })
    .unwrap();
}
