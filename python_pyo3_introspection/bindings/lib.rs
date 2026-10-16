#[path = "../framework/errors.rs"]
mod error_framework;
mod errors;
#[path = "../framework/indexing.rs"]
mod indexing;
#[path = "../framework/parameters.rs"]
mod parameters;
mod params;

#[pyo3::pymodule(name = "_native")]
pub mod native {
    use super::{errors, indexing};
    use cosmolkit as ck;
    use numpy::IntoPyArray;
    use pyo3::prelude::*;

    #[pymodule_export]
    use super::errors::{
        FingerprintError, FingerprintErrorKind, MaccsFingerprintError, MaccsFingerprintErrorKind,
        SmilesError, SmilesErrorKind,
    };

    #[pymodule_export]
    use super::params::MaccsFingerprintParams;

    #[pyfunction]
    fn mol_from_smiles(text: &str) -> PyResult<Molecule> {
        Molecule::from_smiles(text)
    }

    #[pyfunction]
    fn version() -> &'static str {
        ck::version()
    }

    #[pyclass(module = "fixed_projection._native", frozen)]
    pub struct Molecule {
        pub(crate) inner: ck::Molecule,
    }

    #[pymethods]
    impl Molecule {
        #[staticmethod]
        fn from_smiles(text: &str) -> PyResult<Self> {
            ck::Molecule::from_smiles(text)
                .map(|inner| Self { inner })
                .map_err(errors::SmilesError::from_rust)
        }
        fn num_atoms(&self) -> usize {
            self.inner.num_atoms()
        }
        #[pyo3(signature = (params=None, *, n_bits=None))]
        fn fingerprint_maccs(
            &self,
            params: Option<&MaccsFingerprintParams>,
            n_bits: Option<usize>,
        ) -> PyResult<Fingerprint> {
            if let Some(params) = params {
                if n_bits.is_some() {
                    return Err(pyo3::exceptions::PyTypeError::new_err(
                        "parameter object and configuration keywords are mutually exclusive",
                    ));
                }
                self._fingerprint_maccs_params(Some(params))
            } else {
                self._fingerprint_maccs_keywords(n_bits)
            }
        }

        #[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
        #[pyo3(signature = (params=None))]
        fn _fingerprint_maccs_params(
            &self,
            params: Option<&MaccsFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            let defaults;
            let effective = if let Some(params) = params {
                &params.inner
            } else {
                defaults = ck::MaccsFingerprintParams::default();
                &defaults
            };
            self.inner
                .fingerprint_maccs_with_params(effective)
                .map(|inner| Fingerprint { inner })
                .map_err(errors::MaccsFingerprintError::from_rust)
        }

        #[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
        #[pyo3(signature = (*, n_bits=None))]
        fn _fingerprint_maccs_keywords(&self, n_bits: Option<usize>) -> PyResult<Fingerprint> {
            let params = MaccsFingerprintParams {
                inner: ck::MaccsFingerprintParams {
                    n_bits: n_bits.unwrap_or_else(|| ck::MaccsFingerprintParams::default().n_bits),
                },
            };
            self._fingerprint_maccs_params(Some(&params))
        }
    }

    #[pyclass(module = "fixed_projection._native", frozen)]
    pub struct Fingerprint {
        inner: ck::Fingerprint,
    }

    #[pymethods]
    impl Fingerprint {
        #[staticmethod]
        fn from_on_bits(n_bits: u32, on_bits: Vec<u32>) -> PyResult<Self> {
            ck::Fingerprint::from_on_bits(n_bits, on_bits)
                .map(|inner| Self { inner })
                .map_err(errors::FingerprintError::from_rust)
        }
        fn n_bits(&self) -> u32 {
            self.inner.n_bits()
        }
        fn on_bits(&self) -> Vec<u32> {
            self.inner.on_bits()
        }
        fn get_bit(&self, bit: u32) -> PyResult<bool> {
            self.inner
                .get_bit(bit)
                .map_err(errors::FingerprintError::from_rust)
        }
        #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int], numpy.dtype[numpy.uint8]]"}"#]
        fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray1<u8>>> {
            let values = (0..self.inner.n_bits())
                .map(|bit| self.inner.get_bit(bit).map(u8::from))
                .collect::<Result<Vec<_>, _>>()
                .map_err(errors::FingerprintError::from_rust)?;
            Ok(values.into_pyarray(py))
        }
    }

    #[pyclass(module = "fixed_projection._native", frozen)]
    pub struct FingerprintBatch {
        values: Vec<Option<Py<Fingerprint>>>,
    }

    #[pymethods]
    impl FingerprintBatch {
        #[new]
        fn new(values: Vec<Option<Py<Fingerprint>>>) -> Self {
            Self { values }
        }
        fn __len__(&self) -> usize {
            self.values.len()
        }
        fn __getitem__(
            &self,
            py: Python<'_>,
            index: indexing::Index<'_>,
        ) -> PyResult<indexing::Indexed<Option<Py<Fingerprint>>, Py<Self>>> {
            indexing::get(py, index, &self.values, |values| Self { values })
        }
        #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.uint8]]"}"#]
        fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray2<u8>>> {
            let width = self
                .values
                .first()
                .and_then(Option::as_ref)
                .map_or(0, |fp| fp.borrow(py).inner.n_bits() as usize);
            let mut data = numpy::ndarray::Array2::zeros((self.values.len(), width));
            for (row_index, (item, mut row)) in self.values.iter().zip(data.rows_mut()).enumerate()
            {
                let fingerprint = item
                    .as_ref()
                    .ok_or_else(|| {
                        pyo3::exceptions::PyValueError::new_err(format!(
                            "fingerprint row {row_index} is None"
                        ))
                    })?
                    .borrow(py);
                if fingerprint.inner.n_bits() as usize != width {
                    return Err(pyo3::exceptions::PyValueError::new_err(format!(
                        "fingerprint row {row_index} has a different width"
                    )));
                }
                for (bit, value) in row.iter_mut().enumerate() {
                    *value = u8::from(
                        fingerprint
                            .inner
                            .get_bit(bit as u32)
                            .map_err(errors::FingerprintError::from_rust)?,
                    );
                }
            }
            Ok(data.into_pyarray(py))
        }
    }

    // Register SMILES configuration calls and the molecule collection.
    #[cfg(feature = "extension")]
    #[pymodule_export]
    use super::extension::{MoleculeBatch, SmilesParseParams, mol_from_smiles_configured};

    #[cfg(feature = "extension")]
    #[pymodule_export]
    use super::extension::{
        _mol_from_smiles_configured_keywords, _mol_from_smiles_configured_params,
    };
}

#[cfg(feature = "extension")]
mod extension;
