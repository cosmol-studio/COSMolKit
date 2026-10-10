//! Detached value and option projections. All operations delegate to cosmolkit.

use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;

pyo3::create_exception!(
    cosmolkit,
    SmilesError,
    PyValueError,
    "SMILES text could not be parsed or chemically prepared as requested."
);
pyo3::create_exception!(
    cosmolkit,
    SmilesWriteError,
    PyValueError,
    "The molecule could not be serialized as SMILES with the supplied options."
);
pyo3::create_exception!(
    cosmolkit,
    MorganReadError,
    PyValueError,
    "Morgan fingerprint preparation could not read the required molecular state."
);
pyo3::create_exception!(
    cosmolkit,
    FingerprintPreparationError,
    PyValueError,
    "The molecular state required by fingerprint generation could not be prepared."
);
pyo3::create_exception!(
    cosmolkit,
    AtomPairReadError,
    PyValueError,
    "Atom-pair fingerprint preparation could not read the required molecular state."
);
pyo3::create_exception!(
    cosmolkit,
    TopologicalTorsionReadError,
    PyValueError,
    "Topological-torsion fingerprint preparation could not read the required molecular state."
);
pyo3::create_exception!(
    cosmolkit,
    FingerprintError,
    PyValueError,
    "Fingerprint construction, indexing or conversion failed."
);
pyo3::create_exception!(
    cosmolkit,
    FingerprintJsonError,
    PyValueError,
    "A fingerprint JSON representation could not be parsed or serialized."
);

pub(crate) fn source_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    if let Some(error) = source.downcast_ref::<ck::SmartsParseError>() {
        return crate::canonical_search::parse_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::SmartsWriteError>() {
        return crate::canonical_search::write_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::SubstructMatchError>() {
        return crate::canonical_search::substruct_pyerr(py, error.clone());
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionModelError>() {
        return crate::canonical_reaction::model_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionParseError>() {
        return crate::canonical_reaction::parse_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionRunError>() {
        return crate::canonical_reaction::run_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionApplyError>() {
        return crate::canonical_reaction::apply_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionProductError>() {
        return crate::canonical_reaction::product_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionWriteError>() {
        return crate::canonical_reaction::write_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionValidationError>() {
        return crate::canonical_reaction::validation_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::ReactionInitializationError>() {
        return crate::canonical_reaction::initialization_error(py, error);
    }
    // Preserve each batch row's typed scalar cause through the same canonical
    // mapper as a direct scalar call, without copying or reconstructing errors.
    if let Some(error) = source.downcast_ref::<ck::SmilesError>() {
        return smiles_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::SmilesWriteError>() {
        return smiles_write_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::MorganReadError>() {
        return morgan_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::AtomPairReadError>() {
        return atom_pair_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::LayeredFingerprintError>() {
        return crate::canonical_layered::layered_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::PatternFingerprintError>() {
        return crate::canonical_pattern::pattern_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::OperationError>() {
        return crate::drawing_binding::operation_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::DrawingError>() {
        return crate::drawing_binding::drawing_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::DrawingWriteError>() {
        return crate::drawing_binding::drawing_write_pyerr(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::FingerprintError>() {
        return fingerprint_pyerr(py, *error);
    }
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
    if let Some(error) = source.downcast_ref::<ck::BatchValidationError>() {
        return crate::canonical_batch::batch_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::BatchImageError>() {
        return crate::canonical_batch::batch_image_error(py, error);
    }
    if let Some(error) = source.downcast_ref::<ck::TopologicalTorsionReadError>() {
        return topological_torsion_pyerr(py, error);
    }
    if let Some(source) = source.downcast_ref::<ck::BatchError>() {
        let error = annotate(
            py,
            PyValueError::new_err(source.to_string()),
            "batch",
            "Record",
            source,
        );
        let attributes = || -> PyResult<()> {
            error.value(py).setattr("index", source.index)?;
            error.value(py).setattr("operation", source.operation)?;
            error.value(py).setattr("message", &source.message)?;
            Ok(())
        };
        return match attributes() {
            Ok(()) => error,
            Err(error) => error,
        };
    }

    if let Some(error) = source.downcast_ref::<ck::EnumerationError>() {
        return crate::canonical_stereoisomers::error_pyerr(py, error);
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
    if let Some(error) = crate::canonical_registered_errors::convert(py, source) {
        return error;
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

pub(crate) fn smiles_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::SmilesError>,
) -> PyErr {
    let source = source.borrow();
    let kind = match source {
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
        source,
    )
}

pub(crate) fn smiles_write_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::SmilesWriteError>,
) -> PyErr {
    let source = source.borrow();
    let kind = match source {
        ck::SmilesWriteError::Write(_) => "Write",
        ck::SmilesWriteError::Fragment(_) => "Fragment",
    };
    annotate(
        py,
        SmilesWriteError::new_err(source.to_string()),
        "smiles",
        kind,
        source,
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

pub(crate) fn morgan_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::MorganReadError>,
) -> PyErr {
    let source = source.borrow();
    let kind = match source {
        ck::MorganReadError::Preparation(_) => "Preparation",
        ck::MorganReadError::Generator(_) => "Generator",
    };
    annotate(
        py,
        MorganReadError::new_err(source.to_string()),
        "fingerprints",
        kind,
        source,
    )
}

pub(crate) fn atom_pair_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::AtomPairReadError>,
) -> PyErr {
    let source = source.borrow();
    let kind = match source {
        ck::AtomPairReadError::Preparation(_) => "Preparation",
        ck::AtomPairReadError::Generator(_) => "Generator",
    };
    annotate(
        py,
        AtomPairReadError::new_err(source.to_string()),
        "fingerprints",
        kind,
        source,
    )
}

pub(crate) fn topological_torsion_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::TopologicalTorsionReadError>,
) -> PyErr {
    let source = source.borrow();
    let kind = match &source {
        ck::TopologicalTorsionReadError::Preparation(_) => "Preparation",
        ck::TopologicalTorsionReadError::Generator(_) => "Generator",
    };
    annotate(
        py,
        TopologicalTorsionReadError::new_err(source.to_string()),
        "fingerprints",
        kind,
        source,
    )
}

fn fingerprint_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::FingerprintError>,
) -> PyErr {
    let source = source.borrow();
    use ck::FingerprintError as E;
    let kind = match source {
        E::EmptyFingerprint => "EmptyFingerprint",
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
        source,
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
            E::EmptyFingerprint | E::Unsupported => (),
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(e) => e,
    }
}

/// Writable configuration for SMILES parsing.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct SmilesParseParams {
    pub(crate) inner: ck::SmilesParseParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmilesParseParams {
    /// Configure SMILES parsing; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature = (*, sanitize=true, allow_cxsmiles=true, strict_cxsmiles=true, parse_name=true, remove_hs=true, skip_cleanup=false, debug_parse=false, replacements=None))]
    fn new(
        sanitize: bool,
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
        parse_name: bool,
        remove_hs: bool,
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
                remove_hs,
                skip_cleanup,
                debug_parse,
                replacements: replacements.unwrap_or_default(),
            },
        }
    }
    /// Apply to a new molecule and return the result: perform the selected chemical sanitization stages. The source molecule is unchanged.
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    /// Whether a CX extension following the graph notation is parsed.
    #[getter]
    fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles
    }
    /// Whether malformed CX extension data is rejected.
    #[getter]
    fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles
    }
    /// Whether trailing text is interpreted as the molecule/query name.
    #[getter]
    fn parse_name(&self) -> bool {
        self.inner.parse_name
    }
    /// Whether removable explicit hydrogens are removed during input conversion.
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    /// Whether parser post-processing is skipped.
    #[getter]
    fn skip_cleanup(&self) -> bool {
        self.inner.skip_cleanup
    }
    /// Parser debug-output level.
    #[getter]
    fn debug_parse(&self) -> bool {
        self.inner.debug_parse
    }
    /// Text substitutions applied before parsing.
    #[getter]
    fn replacements(&self) -> BTreeMap<String, String> {
        self.inner.replacements.clone()
    }
}

/// Writable configuration for SMILES serialization.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct SmilesWriteParams {
    pub(crate) inner: ck::SmilesWriteParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmilesWriteParams {
    /// Configure SMILES serialization; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature = (*, isomeric_smiles=true, kekule=false, canonical=true, clean_stereo=true, rooted_at_atom=None, all_bonds_explicit=false, all_hydrogens_explicit=false, include_dative_bonds=true, ignore_atom_map_numbers=false))]
    fn new(
        isomeric_smiles: bool,
        kekule: bool,
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
                isomeric_smiles,
                kekule,
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
    /// Whether isotope and stereochemical information is included in the output notation.
    #[getter]
    fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles
    }
    /// Whether aromatic systems are written with explicit single/double bonds.
    #[getter]
    fn kekule(&self) -> bool {
        self.inner.kekule
    }
    /// Whether canonical atom traversal is used for output.
    #[getter]
    fn canonical(&self) -> bool {
        self.inner.canonical
    }
    /// Whether invalid stereochemical annotations are cleaned before SMILES output.
    #[getter]
    fn clean_stereo(&self) -> bool {
        self.inner.clean_stereo
    }
    /// Whether all bonds, including single bonds, have explicit output symbols.
    #[getter]
    fn all_bonds_explicit(&self) -> bool {
        self.inner.all_bonds_explicit
    }
    /// Whether hydrogen counts are written explicitly on every atom.
    #[getter]
    fn all_hydrogens_explicit(&self) -> bool {
        self.inner.all_hydrogens_explicit
    }
    /// Whether dative bonds are included in the requested graph operation/output.
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    /// Whether atom-map numbers are excluded from canonical ranking.
    #[getter]
    fn ignore_atom_map_numbers(&self) -> bool {
        self.inner.ignore_atom_map_numbers
    }
    /// Atom index at which output traversal starts, or None for the default traversal.
    #[getter]
    fn rooted_at_atom(&self) -> Option<usize> {
        self.inner.rooted_at_atom.map(|id| id.index())
    }
}

/// Fixed-width binary fingerprint backed by Rust.
///
/// Use on_bits() for set-bit indices and to_numpy() for a uint8 vector of shape
/// (n_bits,). Logical bit width is distinct from population count.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct Fingerprint {
    pub(crate) inner: ck::Fingerprint,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Fingerprint {
    /// Construct a bit fingerprint of the given width from zero-based set-bit indices; out-of-range indices raise FingerprintError.
    #[staticmethod]
    fn from_on_bits(py: Python<'_>, n_bits: u32, on_bits: Vec<u32>) -> PyResult<Self> {
        ck::Fingerprint::from_on_bits(n_bits, on_bits)
            .map(|inner| Self { inner })
            .map_err(|error| fingerprint_pyerr(py, error))
    }
    /// Return intersection/union Tanimoto similarity with an equal-width bit fingerprint; unequal widths raise FingerprintError.
    fn tanimoto(&self, py: Python<'_>, other: &Self) -> PyResult<f64> {
        self.inner
            .tanimoto(&other.inner)
            .map_err(|error| fingerprint_pyerr(py, error))
    }
    /// Return the logical fingerprint width in bits, not the number of set bits.
    fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    /// Return set-bit indices in increasing logical bit order.
    fn on_bits(&self) -> Vec<u32> {
        self.inner.on_bits()
    }
    /// Return an independent C-contiguous uint8 NumPy array of shape (n_bits,), containing 0 or 1 in logical bit-index order.
    #[gen_stub(override_return_type(type_repr = "numpy.typing.NDArray[numpy.uint8]", imports = ("numpy", "numpy.typing")))]
    fn to_numpy<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<u8>> {
        crate::fingerprint_numpy::to_numpy(py, &self.inner)
    }
    fn __len__(&self) -> usize {
        self.inner.n_bits() as usize
    }
    fn __repr__(&self) -> String {
        format!("Fingerprint(n_bits={})", self.inner.n_bits())
    }
}

/// Sparse bit vector (``SparseBitVect`` equivalent).
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SparseBitFingerprint {
    pub(crate) inner: ck::SparseBitFingerprint,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SparseBitFingerprint {
    /// Logical width of the fingerprint in bits.
    fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    /// On bits as raw signed storage values, in the set's signed ascending
    /// order (source ``IntVect`` is ``std::vector<int>``).
    ///
    /// Source: ``void SparseBitVect::getOnBits(IntVect &v) const``
    /// (SparseBitVect.cpp:252-263): the ``set<int>`` iterates in signed
    /// order, so indices >= 2^31 appear first as wrapped negative values
    /// (oracle: {5, 2^31+1} -> [-2147483647, 5]; {7, u32::MAX} -> [-1, 7]).
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

/// Sparse count fingerprint with 64-bit feature indices and integer counts for nonzero entries.
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
    /// Construct a SparseCountFingerprint value from the supplied inputs.
    #[staticmethod]
    fn new(length: u64) -> Self {
        Self {
            inner: ck::SparseCountFingerprint::new(length),
        }
    }
    /// Logical size of the sparse fingerprint index space.
    fn length(&self) -> u64 {
        self.inner.length()
    }
    /// Return the count at the supplied sparse fingerprint index; absent entries have count zero.
    fn value(&self, py: Python<'_>, index: u64) -> PyResult<i32> {
        self.inner
            .value(index)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Set the count at one sparse fingerprint index; a zero count removes the nonzero entry.
    fn set_value(&mut self, py: Python<'_>, index: u64, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a dictionary mapping sparse fingerprint feature indices to nonzero counts.
    fn nonzero_elements(&self) -> BTreeMap<u64, i32> {
        self.inner.nonzero_elements().clone()
    }
    /// Return the sum of the stored sparse fingerprint counts.
    #[pyo3(signature = (use_abs=false))]
    fn total_value(&self, py: Python<'_>, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a new sparse count fingerprint containing the element-wise minimum counts.
    fn fuzzy_and(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a new sparse count fingerprint containing the element-wise maximum counts.
    fn fuzzy_or(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return the element-wise sum of two compatible sparse count fingerprints.
    fn with_added(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return the element-wise difference of two compatible sparse count fingerprints.
    fn with_subtracted(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with the scalar added to stored nonzero counts.
    fn with_added_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with the scalar subtracted from stored nonzero counts.
    fn with_subtracted_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with stored counts multiplied by the scalar.
    fn with_multiplied_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with stored counts divided by the scalar using integer arithmetic.
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

/// Sparse count fingerprint with 32-bit feature indices and integer counts for nonzero entries.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SparseCountFingerprint32 {
    pub(crate) inner: ck::SparseCountFingerprint32,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SparseCountFingerprint32 {
    /// Construct a SparseCountFingerprint32 value from the supplied inputs.
    #[staticmethod]
    fn new(length: u32) -> Self {
        Self {
            inner: ck::SparseCountFingerprint32::new(length),
        }
    }
    /// Logical size of the sparse fingerprint index space.
    fn length(&self) -> u32 {
        self.inner.length()
    }
    /// Return the count at the supplied sparse fingerprint index; absent entries have count zero.
    fn value(&self, py: Python<'_>, index: u32) -> PyResult<i32> {
        self.inner
            .value(index)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Set the count at one sparse fingerprint index; a zero count removes the nonzero entry.
    fn set_value(&mut self, py: Python<'_>, index: u32, value: i32) -> PyResult<()> {
        self.inner
            .set_value(index, value)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a dictionary mapping sparse fingerprint feature indices to nonzero counts.
    fn nonzero_elements(&self) -> BTreeMap<u32, i32> {
        self.inner.nonzero_elements().clone()
    }
    /// Return the sum of the stored sparse fingerprint counts.
    #[pyo3(signature = (use_abs=false))]
    fn total_value(&self, py: Python<'_>, use_abs: bool) -> PyResult<i32> {
        self.inner
            .total_value(use_abs)
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a new sparse count fingerprint containing the element-wise minimum counts.
    fn fuzzy_and(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a new sparse count fingerprint containing the element-wise maximum counts.
    fn fuzzy_or(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return the element-wise sum of two compatible sparse count fingerprints.
    fn with_added(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return the element-wise difference of two compatible sparse count fingerprints.
    fn with_subtracted(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with the scalar added to stored nonzero counts.
    fn with_added_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_added_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with the scalar subtracted from stored nonzero counts.
    fn with_subtracted_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_subtracted_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with stored counts multiplied by the scalar.
    fn with_multiplied_scalar(&self, py: Python<'_>, value: i32) -> PyResult<Self> {
        self.inner
            .with_multiplied_scalar(value)
            .map(|inner| Self { inner })
            .map_err(|e| fingerprint_pyerr(py, e))
    }
    /// Return a sparse count fingerprint with stored counts divided by the scalar using integer arithmetic.
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

/// Return the version identifier recorded by this value or implementation.
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
    module.add_class::<crate::fingerprint_numpy::FingerprintBatch>()?;
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
