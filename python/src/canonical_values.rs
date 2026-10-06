//! Detached value and option projections. All operations delegate to cosmolkit.

use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;

pyo3::create_exception!(cosmolkit, SmilesError, PyValueError);
pyo3::create_exception!(cosmolkit, SmilesWriteError, PyValueError);
pyo3::create_exception!(cosmolkit, MorganReadError, PyValueError);
pyo3::create_exception!(cosmolkit, FingerprintPreparationError, PyValueError);
pyo3::create_exception!(cosmolkit, AtomPairReadError, PyValueError);
pyo3::create_exception!(cosmolkit, TopologicalTorsionReadError, PyValueError);
pyo3::create_exception!(cosmolkit, FingerprintError, PyValueError);
pyo3::create_exception!(cosmolkit, FingerprintJsonError, PyValueError);

pub(crate) fn source_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    if let Some(error) = source.downcast_ref::<ck::BioStructureError>() {
        return crate::canonical_bio_binding::structure_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ValenceError>() {
        return crate::canonical_atom_bond::valence_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::KekulizeError>() {
        return crate::canonical_chemistry_values::kekulize_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::SanitizeError>() {
        return crate::canonical_chemistry_values::sanitize_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::MatrixError>() {
        return crate::canonical_chemistry_values::matrix_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::CoordinateInputError>() {
        return crate::canonical_coordinate_input::error_pyerr(py, error);
    }
    if let Some(preparation) = source.downcast_ref::<ck::FingerprintPreparationError>() {
        let kind = match preparation {
            ck::FingerprintPreparationError::MissingPreparedValence => "MissingPreparedValence",
            ck::FingerprintPreparationError::RingPreparation(_) => "RingPreparation",
        };
        return annotate(
            py,
            FingerprintPreparationError::new_err(preparation.to_string()),
            "fingerprints",
            kind,
            preparation,
        );
    }
    let error = PyValueError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|cause| source_pyerr(py, cause)));
    error
}

pub(crate) fn annotate(
    py: Python<'_>,
    error: PyErr,
    domain: &str,
    kind: &str,
    source: &(dyn std::error::Error + 'static),
) -> PyErr {
    if let Err(e) = error
        .value(py)
        .setattr("domain", domain)
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return e;
    }
    error.set_cause(py, source.source().map(|cause| source_pyerr(py, cause)));
    error
}

pub(crate) fn smiles_pyerr(py: Python<'_>, source: ck::SmilesError) -> PyErr {
    let kind = match &source {
        ck::SmilesError::Parse(_) => "Parse",
        ck::SmilesError::Hydrogen(_) => "Hydrogen",
        ck::SmilesError::Sanitize(_) => "Sanitize",
        ck::SmilesError::Stereo(_) => "Stereo",
        ck::SmilesError::Construction(_) => "Construction",
    };
    annotate(
        py,
        SmilesError::new_err(source.to_string()),
        "smiles",
        kind,
        &source,
    )
}

pub(crate) fn smiles_write_pyerr(py: Python<'_>, source: ck::SmilesWriteError) -> PyErr {
    let kind = match &source {
        ck::SmilesWriteError::Write(_) => "Write",
        ck::SmilesWriteError::Fragment(_) => "Fragment",
    };
    annotate(
        py,
        SmilesWriteError::new_err(source.to_string()),
        "smiles",
        kind,
        &source,
    )
}

pub(crate) fn fingerprint_json_pyerr(py: Python<'_>, source: ck::FingerprintJsonError) -> PyErr {
    let kind = match &source {
        ck::FingerprintJsonError::Parse(_) => "Parse",
        ck::FingerprintJsonError::Invalid(_) => "Invalid",
        ck::FingerprintJsonError::UnsupportedComponent { .. } => "UnsupportedComponent",
    };
    annotate(
        py,
        FingerprintJsonError::new_err(source.to_string()),
        "fingerprints",
        kind,
        &source,
    )
}

pub(crate) fn morgan_pyerr(py: Python<'_>, source: ck::MorganReadError) -> PyErr {
    let kind = match &source {
        ck::MorganReadError::Preparation(_) => "Preparation",
        ck::MorganReadError::Generator(_) => "Generator",
    };
    annotate(
        py,
        MorganReadError::new_err(source.to_string()),
        "fingerprints",
        kind,
        &source,
    )
}

pub(crate) fn atom_pair_pyerr(py: Python<'_>, source: ck::AtomPairReadError) -> PyErr {
    let kind = match &source {
        ck::AtomPairReadError::Preparation(_) => "Preparation",
        ck::AtomPairReadError::Generator(_) => "Generator",
    };
    annotate(
        py,
        AtomPairReadError::new_err(source.to_string()),
        "fingerprints",
        kind,
        &source,
    )
}

pub(crate) fn topological_torsion_pyerr(
    py: Python<'_>,
    source: ck::TopologicalTorsionReadError,
) -> PyErr {
    let kind = match &source {
        ck::TopologicalTorsionReadError::Preparation(_) => "Preparation",
        ck::TopologicalTorsionReadError::Generator(_) => "Generator",
    };
    annotate(
        py,
        TopologicalTorsionReadError::new_err(source.to_string()),
        "fingerprints",
        kind,
        &source,
    )
}

fn fingerprint_pyerr(py: Python<'_>, source: ck::FingerprintError) -> PyErr {
    use ck::FingerprintError as E;
    let kind = match source {
        E::Unsupported => "Unsupported",
        E::SparseIndexOutOfRange { .. } => "SparseIndexOutOfRange",
        E::BitLengthMismatch { .. } => "BitLengthMismatch",
        E::InvalidFoldFactor { .. } => "InvalidFoldFactor",
        E::RangeError { .. } => "RangeError",
        E::UndefinedArithmetic { .. } => "UndefinedArithmetic",
        E::PreconditionViolation { .. } => "PreconditionViolation",
        E::InvalidArguments { .. } => "InvalidArguments",
    };
    let error = annotate(
        py,
        FingerprintError::new_err(source.to_string()),
        "fingerprints",
        kind,
        &source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match source {
            E::SparseIndexOutOfRange { index, size } => {
                value.setattr("index", index)?;
                value.setattr("size", size)?;
            }
            E::BitLengthMismatch { left, right } => {
                value.setattr("left", left)?;
                value.setattr("right", right)?;
            }
            E::InvalidFoldFactor { factor, n_bits } => {
                value.setattr("factor", factor)?;
                value.setattr("n_bits", n_bits)?;
            }
            E::RangeError { value: parameter } => value.setattr("value", parameter)?,
            E::UndefinedArithmetic { site } => value.setattr("site", site)?,
            E::PreconditionViolation { what } => value.setattr("what", what)?,
            E::InvalidArguments { reason } => value.setattr("reason", reason)?,
            E::Unsupported => (),
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(e) => e,
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SmilesParseParams {
    pub(crate) inner: ck::SmilesParseParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmilesParseParams {
    #[new]
    #[pyo3(signature = (*, sanitize=true, allow_cxsmiles=true, strict_cxsmiles=true, parse_name=true, remove_hydrogens=true, skip_cleanup=false, debug_parse=false, replacements=None))]
    fn new(
        sanitize: bool,
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
        parse_name: bool,
        remove_hydrogens: bool,
        skip_cleanup: bool,
        debug_parse: bool,
        replacements: Option<BTreeMap<String, String>>,
    ) -> Self {
        // Project every registered canonical field without a parsing policy here.
        Self {
            inner: ck::SmilesParseParams {
                sanitize,
                allow_cxsmiles,
                strict_cxsmiles,
                parse_name,
                remove_hydrogens,
                skip_cleanup,
                debug_parse,
                replacements: replacements.unwrap_or_default(),
            },
        }
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[getter]
    fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles
    }
    #[getter]
    fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles
    }
    #[getter]
    fn parse_name(&self) -> bool {
        self.inner.parse_name
    }
    #[getter]
    fn remove_hydrogens(&self) -> bool {
        self.inner.remove_hydrogens
    }
    #[getter]
    fn skip_cleanup(&self) -> bool {
        self.inner.skip_cleanup
    }
    #[getter]
    fn debug_parse(&self) -> bool {
        self.inner.debug_parse
    }
    #[getter]
    fn replacements(&self) -> BTreeMap<String, String> {
        self.inner.replacements.clone()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SmilesWriteParams {
    pub(crate) inner: ck::SmilesWriteParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmilesWriteParams {
    #[new]
    #[pyo3(signature = (*, do_isomeric_smiles=true, do_kekule=false, canonical=true, clean_stereo=true, rooted_at_atom=None, all_bonds_explicit=false, all_hydrogens_explicit=false, include_dative_bonds=true, ignore_atom_map_numbers=false))]
    fn new(
        do_isomeric_smiles: bool,
        do_kekule: bool,
        canonical: bool,
        clean_stereo: bool,
        rooted_at_atom: Option<usize>,
        all_bonds_explicit: bool,
        all_hydrogens_explicit: bool,
        include_dative_bonds: bool,
        ignore_atom_map_numbers: bool,
    ) -> Self {
        Self {
            inner: ck::SmilesWriteParams {
                do_isomeric_smiles,
                do_kekule,
                canonical,
                clean_stereo,
                rooted_at_atom: rooted_at_atom.map(ck::AtomId::new),
                all_bonds_explicit,
                all_hydrogens_explicit,
                include_dative_bonds,
                ignore_atom_map_numbers,
            },
        }
    }
    #[getter]
    fn do_isomeric_smiles(&self) -> bool {
        self.inner.do_isomeric_smiles
    }
    #[getter]
    fn do_kekule(&self) -> bool {
        self.inner.do_kekule
    }
    #[getter]
    fn canonical(&self) -> bool {
        self.inner.canonical
    }
    #[getter]
    fn clean_stereo(&self) -> bool {
        self.inner.clean_stereo
    }
    #[getter]
    fn all_bonds_explicit(&self) -> bool {
        self.inner.all_bonds_explicit
    }
    #[getter]
    fn all_hydrogens_explicit(&self) -> bool {
        self.inner.all_hydrogens_explicit
    }
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[getter]
    fn ignore_atom_map_numbers(&self) -> bool {
        self.inner.ignore_atom_map_numbers
    }
    #[getter]
    fn rooted_at_atom(&self) -> Option<usize> {
        self.inner.rooted_at_atom.map(|id| id.index())
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct Fingerprint {
    pub(crate) inner: ck::Fingerprint,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Fingerprint {
    #[staticmethod]
    fn from_on_bits(py: Python<'_>, n_bits: u32, on_bits: Vec<u32>) -> PyResult<Self> {
        ck::Fingerprint::from_on_bits(n_bits, on_bits)
            .map(|inner| Self { inner })
            .map_err(|error| fingerprint_pyerr(py, error))
    }
    fn tanimoto(&self, py: Python<'_>, other: &Self) -> PyResult<f64> {
        self.inner
            .tanimoto(&other.inner)
            .map_err(|error| fingerprint_pyerr(py, error))
    }
    fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    fn on_bits(&self) -> Vec<u32> {
        self.inner.on_bits()
    }
    fn __len__(&self) -> usize {
        self.inner.n_bits() as usize
    }
    fn __repr__(&self) -> String {
        format!("Fingerprint(n_bits={})", self.inner.n_bits())
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SparseBitFingerprint {
    pub(crate) inner: ck::SparseBitFingerprint,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SparseBitFingerprint {
    fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    fn on_bits(&self) -> Vec<i32> {
        self.inner.on_bits()
    }
    fn __len__(&self) -> usize {
        self.inner.n_bits() as usize
    }
    fn __repr__(&self) -> String {
        format!("SparseBitFingerprint(n_bits={})", self.inner.n_bits())
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SparseCountFingerprint {
    pub(crate) inner: ck::SparseCountFingerprint,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SparseCountFingerprint {
    fn __len__(&self) -> PyResult<usize> {
        // Original d892 binding uses usize::try_from; the PyO3 len slot then
        // enforces Python Py_ssize_t bounds on the returned usize itself.
        usize::try_from(self.inner.length()).map_err(|_| {
            pyo3::exceptions::PyOverflowError::new_err(
                "fingerprint size exceeds Python platform size",
            )
        })
    }
    #[staticmethod]
    fn new(length: u64) -> Self {
        Self {
            inner: ck::SparseCountFingerprint::new(length),
        }
    }
    fn length(&self) -> u64 {
        self.inner.length()
    }
    fn value(&self, py: Python<'_>, index: u64) -> PyResult<i32> {
        self.inner
            .value(index)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn set_value(&mut self, py: Python<'_>, index: u64, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn nonzero_elements(&self) -> BTreeMap<u64, i32> {
        self.inner.nonzero_elements().clone()
    }
    #[pyo3(signature = (use_abs=false))]
    fn total_value(&self, py: Python<'_>, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn fuzzy_and(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn fuzzy_or(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_added(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_subtracted(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_added_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_subtracted_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_multiplied_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_divided_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_divided_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn __repr__(&self) -> String {
        format!("SparseCountFingerprint(length={})", self.inner.length())
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SparseCountFingerprint32 {
    pub(crate) inner: ck::SparseCountFingerprint32,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SparseCountFingerprint32 {
    #[staticmethod]
    fn new(length: u32) -> Self {
        Self {
            inner: ck::SparseCountFingerprint32::new(length),
        }
    }
    fn length(&self) -> u32 {
        self.inner.length()
    }
    fn value(&self, py: Python<'_>, index: u32) -> PyResult<i32> {
        self.inner
            .value(index)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn set_value(&mut self, py: Python<'_>, index: u32, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn nonzero_elements(&self) -> BTreeMap<u32, i32> {
        self.inner.nonzero_elements().clone()
    }
    #[pyo3(signature = (use_abs=false))]
    fn total_value(&self, py: Python<'_>, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn fuzzy_and(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn fuzzy_or(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_added(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_subtracted(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_added_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_subtracted_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_multiplied_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn with_divided_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_divided_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    fn __repr__(&self) -> String {
        format!("SparseCountFingerprint32(length={})", self.inner.length())
    }
}

#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn version() -> &'static str {
    ck::version()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::canonical_fingerprint_values::register(module)?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_class::<SmilesParseParams>()?;
    module.add_class::<SmilesWriteParams>()?;
    module.add_class::<Fingerprint>()?;
    module.add_class::<SparseBitFingerprint>()?;
    module.add_class::<SparseCountFingerprint>()?;
    module.add_class::<SparseCountFingerprint32>()?;
    module.add("SmilesError", module.py().get_type::<SmilesError>())?;
    module.add(
        "SmilesWriteError",
        module.py().get_type::<SmilesWriteError>(),
    )?;
    module.add(
        "AtomPairReadError",
        module.py().get_type::<AtomPairReadError>(),
    )?;
    module.add(
        "TopologicalTorsionReadError",
        module.py().get_type::<TopologicalTorsionReadError>(),
    )?;
    module.add(
        "FingerprintJsonError",
        module.py().get_type::<FingerprintJsonError>(),
    )?;
    module.add("MorganReadError", module.py().get_type::<MorganReadError>())?;
    module.add(
        "FingerprintPreparationError",
        module.py().get_type::<FingerprintPreparationError>(),
    )?;
    module.add(
        "FingerprintError",
        module.py().get_type::<FingerprintError>(),
    )?;
    Ok(())
}
