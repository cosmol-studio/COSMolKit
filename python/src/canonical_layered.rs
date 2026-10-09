//! Layered projections call the single canonical Rust facade.
use crate::canonical_values::Fingerprint;
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::error::Error;

pyo3::create_exception!(cosmolkit, LayeredFingerprintError, PyValueError);
pub(crate) fn layered_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::LayeredFingerprintError>,
) -> PyErr {
    let source = source.borrow();
    use ck::LayeredFingerprintError as E;
    let err = LayeredFingerprintError::new_err(source.to_string());
    let attrs = || -> PyResult<()> {
        let object = err.value(py);
        object.setattr("domain", "Fingerprint")?;
        object.setattr(
            "kind",
            match source {
                E::InvalidArguments { .. } => "InvalidArguments",
                E::Topology(_) => "Topology",
                E::Query(_) => "Query",
                E::Rings(_) => "Rings",
                E::Paths(_) => "Paths",
                E::Value(_) => "Value",
            },
        )?;
        if let E::InvalidArguments { reason } = source {
            object.setattr("reason", reason)?;
        }
        Ok(())
    };
    if let Err(e) = attrs() {
        return e;
    }
    err.set_cause(
        py,
        source
            .source()
            .map(|e| crate::canonical_values::source_pyerr(py, e)),
    );
    err
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct LayeredFingerprintLayers {
    inner: ck::LayeredFingerprintLayers,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LayeredFingerprintLayers {
    #[staticmethod]
    fn from_bits_retain(bits: u32) -> Self {
        Self {
            inner: ck::LayeredFingerprintLayers::from_bits_retain(bits),
        }
    }
    fn bits(&self) -> u32 {
        self.inner.bits()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct LayeredFingerprintParams {
    pub(crate) inner: ck::LayeredFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LayeredFingerprintParams {
    #[new]
    #[pyo3(signature=(*,layers=0xffff_ffff,min_path=1,max_path=7,fp_size=2048,atom_counts=None,set_only_bits=None,branched_paths=true,from_atoms=None))]
    fn new(
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        atom_counts: Option<Vec<u32>>,
        set_only_bits: Option<&Fingerprint>,
        branched_paths: bool,
        from_atoms: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::LayeredFingerprintParams {
                layers: ck::LayeredFingerprintLayers::from_bits_retain(layers),
                min_path,
                max_path,
                fp_size,
                atom_counts,
                set_only_bits: set_only_bits.map(|value| value.inner.clone()),
                branched_paths,
                from_atoms,
            },
        }
    }
    #[getter]
    fn layers(&self) -> u32 {
        self.inner.layers.bits()
    }
    #[getter]
    fn min_path(&self) -> u32 {
        self.inner.min_path
    }
    #[getter]
    fn max_path(&self) -> u32 {
        self.inner.max_path
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[getter]
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts.clone()
    }
    #[getter]
    fn set_only_bits(&self) -> Option<Fingerprint> {
        self.inner
            .set_only_bits
            .clone()
            .map(|inner| Fingerprint { inner })
    }
    #[getter]
    fn branched_paths(&self) -> bool {
        self.inner.branched_paths
    }
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct LayeredFingerprintResult {
    pub(crate) inner: ck::LayeredFingerprintResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LayeredFingerprintResult {
    fn fingerprint(&self) -> Fingerprint {
        // COSMolKit✔️✔️:     fn fingerprint(&self) -> Fingerprint {
        // COSMolKit✔️✔️:         self.fingerprint.clone()
        // COSMolKit✔️✔️:     }
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    fn atom_counts(&self) -> Option<Vec<u32>> {
        // COSMolKit✔️✔️:     fn atom_counts(&self) -> Option<Vec<u32>> {
        // COSMolKit✔️✔️:         self.atom_counts.clone()
        // COSMolKit✔️✔️:     }
        self.inner.atom_counts().map(<[u32]>::to_vec)
    }
    fn __repr__(&self) -> String {
        // COSMolKit✔️✔️:     fn __repr__(&self) -> String {
        // COSMolKit✔️✔️:         format!(
        // COSMolKit✔️✔️:             "LayeredFingerprintResult(n_bits={}, has_atom_counts={})",
        // COSMolKit✔️✔️:             self.fingerprint.inner.n_bits(),
        // COSMolKit✔️✔️:             self.atom_counts.is_some()
        // COSMolKit✔️✔️:         )
        // COSMolKit✔️✔️:     }
        format!(
            "LayeredFingerprintResult(n_bits={}, has_atom_counts={})",
            self.inner.fingerprint().n_bits(),
            self.inner.atom_counts().is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn fingerprint_layered_query_with_params(
    py: Python<'_>,
    query: &crate::canonical_search::QueryGraph,
    params: &LayeredFingerprintParams,
) -> PyResult<Fingerprint> {
    ck::fingerprint_layered_query_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| layered_pyerr(py, e))
}
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn fingerprint_layered_query_with_output_with_params(
    py: Python<'_>,
    query: &crate::canonical_search::QueryGraph,
    params: &LayeredFingerprintParams,
) -> PyResult<LayeredFingerprintResult> {
    ck::fingerprint_layered_query_with_output_with_params(&query.inner, &params.inner)
        .map(|inner| LayeredFingerprintResult { inner })
        .map_err(|e| layered_pyerr(py, e))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<LayeredFingerprintLayers>()?;
    module.add_class::<LayeredFingerprintParams>()?;
    module.add_class::<LayeredFingerprintResult>()?;
    module.add(
        "LayeredFingerprintError",
        module.py().get_type::<LayeredFingerprintError>(),
    )?;
    module.add_function(wrap_pyfunction!(
        fingerprint_layered_query_with_params,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        fingerprint_layered_query_with_output_with_params,
        module
    )?)?;
    Ok(())
}
