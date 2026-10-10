//! RDKFingerprint projections of canonical Rust values; no chemistry algorithm.
use crate::canonical_values::Fingerprint;
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::{collections::BTreeMap, error::Error};
pyo3::create_exception!(
    cosmolkit,
    TopologicalFingerprintError,
    PyValueError,
    "Topological fingerprint generation failed for the supplied graph or options."
);
pub(crate) fn topological_pyerr(py: Python<'_>, source: ck::TopologicalFingerprintError) -> PyErr {
    use ck::TopologicalFingerprintError as E;
    let err = TopologicalFingerprintError::new_err(source.to_string());
    let attributes = || -> PyResult<()> {
        let object = err.value(py);
        object.setattr("domain", "Fingerprint")?;
        object.setattr(
            "kind",
            match &source {
                E::InvalidArguments { .. } => "InvalidArguments",
                E::OutputNotRequested { .. } => "OutputNotRequested",
                E::Topology(_) => "Topology",
                E::Query(_) => "Query",
                E::Paths(_) => "Paths",
                E::Value(_) => "Value",
            },
        )?;
        match &source {
            E::InvalidArguments { reason } => object.setattr("reason", reason)?,
            E::OutputNotRequested { field } => object.setattr("field", field)?,
            _ => {}
        }
        Ok(())
    };
    if let Err(error) = attributes() {
        return error;
    }
    err.set_cause(
        py,
        source
            .source()
            .map(|cause| crate::canonical_values::source_pyerr(py, cause)),
    );
    err
}
/// Writable configuration for RDKit topological path fingerprints.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct TopologicalFingerprintParams {
    pub(crate) inner: ck::TopologicalFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalFingerprintParams {
    /// Configure RDKit topological path fingerprints; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,min_path=1,max_path=7,fp_size=2048,num_bits_per_feature=2,use_hs=true,target_density=0.0,min_size=128,branched_paths=true,use_bond_order=true,atom_invariants=None,from_atoms=None))]
    fn new(
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        num_bits_per_feature: u32,
        use_hs: bool,
        target_density: f64,
        min_size: u32,
        branched_paths: bool,
        use_bond_order: bool,
        atom_invariants: Option<Vec<u32>>,
        from_atoms: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::TopologicalFingerprintParams {
                ignore_atoms: None,
                min_path,
                max_path,
                fp_size,
                num_bits_per_feature,
                use_hs,
                target_density,
                min_size,
                branched_paths,
                use_bond_order,
                atom_invariants,
                from_atoms,
            },
        }
    }
    /// Minimum number of bonds in an enumerated fingerprint path.
    #[getter]
    fn min_path(&self) -> u32 {
        self.inner.min_path
    }
    /// Maximum number of bonds in an enumerated fingerprint path.
    #[getter]
    fn max_path(&self) -> u32 {
        self.inner.max_path
    }
    /// Number of bins/bits in the folded fingerprint.
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    /// Number of hashed bit positions generated per path feature.
    #[getter]
    fn num_bits_per_feature(&self) -> u32 {
        self.inner.num_bits_per_feature
    }
    /// Whether explicit graph hydrogens participate in fingerprint paths.
    #[getter]
    fn use_hs(&self) -> bool {
        self.inner.use_hs
    }
    /// Target on-bit density for optional iterative fingerprint folding.
    #[getter]
    fn target_density(&self) -> f64 {
        self.inner.target_density
    }
    /// Minimum folded fingerprint width.
    #[getter]
    fn min_size(&self) -> u32 {
        self.inner.min_size
    }
    /// Whether branched subgraphs, not just linear paths, are included.
    #[getter]
    fn branched_paths(&self) -> bool {
        self.inner.branched_paths
    }
    /// Whether bond-order weights are included in the calculation.
    #[getter]
    fn use_bond_order(&self) -> bool {
        self.inner.use_bond_order
    }
    /// Atom-indexed invariants used for path fingerprint hashing.
    #[getter]
    fn atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.atom_invariants.clone()
    }
    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
}
/// Flags selecting atom-to-bit participation and bit-path metadata for a topological fingerprint.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct TopologicalFingerprintOutputRequest {
    pub(crate) inner: ck::TopologicalFingerprintOutputRequest,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalFingerprintOutputRequest {
    /// Construct a TopologicalFingerprintOutputRequest value from the supplied inputs.
    #[new]
    #[pyo3(signature=(*,atom_bits=false,bit_info=false))]
    fn new(atom_bits: bool, bit_info: bool) -> Self {
        Self {
            inner: ck::TopologicalFingerprintOutputRequest {
                atom_bits,
                bit_info,
            },
        }
    }
    /// Atom-indexed fingerprint bit participation, or request flag enabling that sink.
    #[getter]
    fn atom_bits(&self) -> bool {
        self.inner.atom_bits
    }
    /// Map from fingerprint bit indices to contributing paths, or request flag enabling that sink.
    #[getter]
    fn bit_info(&self) -> bool {
        self.inner.bit_info
    }
}
/// Typed provenance returned by the source ``atomBits`` and ``bitInfo`` outputs.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologicalFingerprintOutput {
    inner: ck::TopologicalFingerprintOutput,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalFingerprintOutput {
    /// Atom-indexed fingerprint bit participation, or request flag enabling that sink.
    #[getter]
    fn atom_bits(&self) -> Option<Vec<Vec<u32>>> {
        self.inner.atom_bits.clone()
    }
    /// Map from fingerprint bit indices to contributing paths, or request flag enabling that sink.
    #[getter]
    fn bit_info(&self) -> Option<BTreeMap<u32, Vec<Vec<i32>>>> {
        self.inner.bit_info.clone()
    }
}
/// Topological bit fingerprint plus the requested atom/bit-path metadata.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologicalFingerprintResult {
    pub(crate) inner: ck::TopologicalFingerprintResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalFingerprintResult {
    /// Return the bit fingerprint stored in this result.
    fn fingerprint(&self) -> Fingerprint {
        // COSMolKit❗✔️:     fn fingerprint(&self) -> Fingerprint {
        // COSMolKit❗✔️:         self.fingerprint.clone()
        // COSMolKit❗✔️:     }
        // Complexity: native transport clones the requested detached value once; no regeneration.
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    /// Atom-indexed fingerprint bit participation, or request flag enabling that sink.
    fn atom_bits(&self, py: Python<'_>) -> PyResult<Vec<Vec<u32>>> {
        // COSMolKit❗✔️:     fn atom_bits(&self) -> PyResult<Vec<Vec<u32>>> {
        // COSMolKit❗✔️:         self.atom_bits.clone().ok_or_else(|| {
        // COSMolKit❗✔️:             PyValueError::new_err(
        // COSMolKit❗✔️:                 "topological atom_bits output was not requested for this fingerprint result",
        // COSMolKit❗✔️:             )
        // COSMolKit❗✔️:         })
        // COSMolKit❗✔️:     }
        // Complexity: native transport clones the requested detached value once; no regeneration.
        self.inner
            .atom_bits()
            .map(<[Vec<u32>]>::to_vec)
            .map_err(|error| topological_pyerr(py, error))
    }
    /// Map from fingerprint bit indices to contributing paths, or request flag enabling that sink.
    fn bit_info(&self, py: Python<'_>) -> PyResult<BTreeMap<u32, Vec<Vec<i32>>>> {
        // COSMolKit❗✔️:     fn bit_info(&self) -> PyResult<BTreeMap<u32, Vec<Vec<i32>>>> {
        // COSMolKit❗✔️:         self.bit_info.clone().ok_or_else(|| {
        // COSMolKit❗✔️:             PyValueError::new_err(
        // COSMolKit❗✔️:                 "topological bit_info output was not requested for this fingerprint result",
        // COSMolKit❗✔️:             )
        // COSMolKit❗✔️:         })
        // COSMolKit❗✔️:     }
        // Complexity: native transport clones the requested detached value once; no regeneration.
        self.inner
            .bit_info()
            .cloned()
            .map_err(|error| topological_pyerr(py, error))
    }
    fn __repr__(&self) -> String {
        // COSMolKit❗✔️:     fn __repr__(&self) -> String {
        // COSMolKit❗✔️:         format!(
        // COSMolKit❗✔️:             "TopologicalFingerprintResult(n_bits={}, has_atom_bits={}, has_bit_info={})",
        // COSMolKit❗✔️:             self.fingerprint.inner.n_bits(),
        // COSMolKit❗✔️:             self.atom_bits.is_some(),
        // COSMolKit❗✔️:             self.bit_info.is_some()
        // COSMolKit❗✔️:         )
        // COSMolKit❗✔️:     }
        // Complexity: native transport clones the requested detached value once; no regeneration.
        self.inner.to_string()
    }
}
/// Compute topological fixed-width bit fingerprints for a QueryGraph without discarding its predicates.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn fingerprint_topological_query_with_params(
    py: Python<'_>,
    query: &crate::canonical_search::QueryGraph,
    params: &TopologicalFingerprintParams,
) -> PyResult<Fingerprint> {
    ck::fingerprint_topological_query_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| topological_pyerr(py, e))
}
/// Compute topological fixed-width bit fingerprints for a QueryGraph without discarding its predicates; include the requested atom/bit environment metadata.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn fingerprint_topological_query_with_output_with_params(
    py: Python<'_>,
    query: &crate::canonical_search::QueryGraph,
    params: &TopologicalFingerprintParams,
    request: &TopologicalFingerprintOutputRequest,
) -> PyResult<TopologicalFingerprintResult> {
    ck::fingerprint_topological_query_with_output_with_params(
        &query.inner,
        &params.inner,
        request.inner,
    )
    .map(|inner| TopologicalFingerprintResult { inner })
    .map_err(|e| topological_pyerr(py, e))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<TopologicalFingerprintParams>()?;
    module.add_class::<TopologicalFingerprintOutputRequest>()?;
    module.add_class::<TopologicalFingerprintOutput>()?;
    module.add_class::<TopologicalFingerprintResult>()?;
    module.add(
        "TopologicalFingerprintError",
        module.py().get_type::<TopologicalFingerprintError>(),
    )?;
    module.add_function(wrap_pyfunction!(
        fingerprint_topological_query_with_params,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        fingerprint_topological_query_with_output_with_params,
        module
    )?)?;
    Ok(())
}
