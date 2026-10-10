//! Residue vocabulary and immutable result projections of the public BIO owner.

pyo3::create_exception!(
    cosmolkit,
    ResidueCodeParseError,
    pyo3::exceptions::PyValueError,
    "A residue dictionary code name is not recognized; input() returns the supplied name."
);
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyMapping, PyMappingProxy};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};

pub(crate) fn enum_member<'py>(
    py: Python<'py>,
    name: &str,
    code: i64,
) -> PyResult<Bound<'py, PyAny>> {
    py.import("cosmolkit")?.getattr(name)?.call1((code,))
}

fn kind_from_code(code: i64) -> PyResult<ck::ResidueInfoKind> {
    use ck::ResidueInfoKind as K;
    match code {
        0 => Ok(K::Unknown),
        1 => Ok(K::Aa),
        2 => Ok(K::Aad),
        3 => Ok(K::Paa),
        4 => Ok(K::Maa),
        5 => Ok(K::Rna),
        6 => Ok(K::Dna),
        7 => Ok(K::Buf),
        8 => Ok(K::Hoh),
        9 => Ok(K::Pyr),
        10 => Ok(K::Ket),
        11 => Ok(K::Els),
        _ => Err(PyValueError::new_err(format!(
            "unsupported ResidueInfoKind code {code}"
        ))),
    }
}

/// Read-only residue dictionary entry containing identity, class, sequence codes, composition and standard/modified status.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ResidueInfo {
    pub(crate) inner: ck::ResidueInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ResidueInfo {
    /// Stored residue/atom code; its interpretation is defined by the owning value type.
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueCode", i64::from(self.inner.code.as_u16()))
    }
    /// Stored name of this value.
    fn name(&self) -> &'static str {
        self.inner.name
    }
    /// Return the residue dictionary chemical classification.
    #[gen_stub(override_return_type(type_repr = "ResidueInfoKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueInfoKind", self.inner.kind as i64)
    }
    /// Return the text name of the residue dictionary classification.
    fn kind_name(&self) -> &'static str {
        self.inner.kind.name()
    }
    /// Residue polymer-linking classification.
    fn linking_type(&self) -> u8 {
        self.inner.linking_type
    }
    /// One-letter residue code recorded by the residue dictionary.
    fn one_letter_code(&self) -> String {
        self.inner.one_letter_code.to_string()
    }
    /// Hydrogen count recorded in the residue dictionary.
    fn hydrogen_count(&self) -> u8 {
        self.inner.hydrogen_count
    }
    /// Residue molecular weight from the residue dictionary.
    fn weight(&self) -> f32 {
        self.inner.weight
    }
    /// Whether the residue name was found in the built-in dictionary.
    fn found(&self) -> bool {
        self.inner.found()
    }
    /// Whether the residue dictionary identifies water.
    fn is_water(&self) -> bool {
        self.inner.is_water()
    }
    /// Return whether the dictionary classifies this residue as DNA.
    fn is_dna(&self) -> bool {
        self.inner.is_dna()
    }
    /// Return whether the dictionary classifies this residue as RNA.
    fn is_rna(&self) -> bool {
        self.inner.is_rna()
    }
    /// Whether the residue dictionary identifies a nucleic-acid residue.
    fn is_nucleic_acid(&self) -> bool {
        self.inner.is_nucleic_acid()
    }
    /// Whether the residue dictionary identifies an amino acid.
    fn is_amino_acid(&self) -> bool {
        self.inner.is_amino_acid()
    }
    /// Return whether the residue is classified as buffer or water.
    fn is_buffer_or_water(&self) -> bool {
        self.inner.is_buffer_or_water()
    }
    /// Whether the residue belongs to the standard residue set.
    fn is_standard(&self) -> bool {
        self.inner.is_standard()
    }
    /// Residue code suitable for FASTA sequence output.
    fn fasta_code(&self) -> String {
        self.inner.fasta_code().to_string()
    }
    /// One-letter code of the corresponding standard parent residue.
    fn canonical_one_letter_code(&self) -> Option<String> {
        self.inner
            .canonical_one_letter_code()
            .map(|c| c.to_string())
    }
    /// Standard parent residue code for a modified residue, when known.
    #[gen_stub(override_return_type(type_repr = "typing.Optional[ResidueCode]"))]
    fn parent_standard_code<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .parent_standard_code()
            .map(|c| enum_member(py, "ResidueCode", i64::from(c.as_u16())))
            .transpose()
    }
    /// Whether the residue dictionary identifies a modified amino acid.
    fn is_modified_amino_acid(&self) -> bool {
        self.inner.is_modified_amino_acid()
    }
    /// Return whether the residue participates in peptide polymer linkage.
    fn is_peptide_linking(&self) -> bool {
        self.inner.is_peptide_linking()
    }
    /// Return whether the residue participates in nucleic-acid polymer linkage.
    fn is_na_linking(&self) -> bool {
        self.inner.is_na_linking()
    }
    fn __repr__(&self) -> String {
        format!(
            "ResidueInfo(name='{}', code='{:?}', kind='{}')",
            self.inner.name,
            self.inner.code,
            self.inner.kind.name()
        )
    }
}

/// Look up a residue code and return its dictionary information, preserving the not-found state.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn find_residue_info(name: &str) -> ResidueInfo {
    ResidueInfo {
        inner: ck::find_residue_info(name),
    }
}

/// Return the built-in dictionary index for a residue code when found.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn find_residue_info_index(name: &str) -> usize {
    ck::find_residue_info_index(name)
}

/// Return the residue dictionary entry at the supplied index.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn residue_info(index: usize) -> PyResult<ResidueInfo> {
    ck::residue_info_checked(index)
        .map(|inner| ResidueInfo { inner })
        .ok_or_else(|| PyIndexError::new_err(format!("residue info index {index} out of range")))
}

/// Return the residue dictionary entry at the supplied index; reject an invalid index.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn residue_info_checked(index: usize) -> Option<ResidueInfo> {
    ck::residue_info_checked(index).map(|inner| ResidueInfo { inner })
}

/// Residue name and its dictionary classification, preserving unrecognized names.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ResidueIdentity {
    inner: ck::ResidueIdentity,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ResidueIdentity {
    /// Construct a ResidueIdentity value from the supplied inputs.
    #[staticmethod]
    fn new(name: String) -> Self {
        Self {
            inner: ck::ResidueIdentity::new(name),
        }
    }
    /// Stored name of this value.
    fn name(&self) -> &str {
        self.inner.name()
    }
    /// Stored residue/atom code; its interpretation is defined by the owning value type.
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueCode", i64::from(self.inner.code().as_u16()))
    }
    /// Return the associated residue dictionary information, including its classification and sequence codes.
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.inner.info(),
        }
    }
    /// Return whether the residue name appears in the built-in residue dictionary.
    fn is_tabulated(&self) -> bool {
        self.inner.is_tabulated()
    }
}

/// Return the canonical residue code from a residue identity.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(type_repr = "ResidueCode"))]
fn residue_code<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    enum_member(
        py,
        "ResidueCode",
        i64::from(ck::residue_code(name).as_u16()),
    )
}

/// Expand a one-letter residue code to a three-letter residue code using the chosen polymer convention.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn expand_one_letter(
    code: &str,
    #[gen_stub(override_type(type_repr = "ResidueInfoKind | builtins.str | builtins.int"))]
    kind: &Bound<'_, PyAny>,
) -> PyResult<Option<String>> {
    let kind = crate::canonical_atom_bond::enum_code(kind, "ResidueInfoKind")?;
    let mut chars = code.chars();
    let c = chars
        .next()
        .filter(|_| chars.next().is_none())
        .ok_or_else(|| PyValueError::new_err("code must contain exactly one character"))?;
    Ok(ck::expand_one_letter(c, kind_from_code(kind)?).map(str::to_owned))
}

/// Expand a sequence of one-letter residue codes to residue names in sequence order.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn expand_one_letter_sequence(
    seq: &str,
    #[gen_stub(override_type(type_repr = "ResidueInfoKind | builtins.str | builtins.int"))]
    kind: &Bound<'_, PyAny>,
) -> PyResult<Vec<String>> {
    let kind = crate::canonical_atom_bond::enum_code(kind, "ResidueInfoKind")?;
    ck::expand_one_letter_sequence(seq, kind_from_code(kind)?)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    let int_enum = py.import("enum")?.getattr("IntEnum")?;
    let members = PyDict::new(py);
    let names = PyDict::new(py);
    // The existing checked table owner supplies every original enum member;
    // this projection adds no chemistry lookup table or alias algorithm.
    let mut index = 0;
    while let Some(info) = ck::residue_info_checked(index) {
        let name = format!("{:?}", info.code);
        members.set_item(&name, info.code.as_u16())?;
        names.set_item(info.name, name)?;
        index += 1;
    }
    let codes = int_enum.call1(("ResidueCode", members))?;
    codes.setattr("__module__", "cosmolkit")?;
    codes.setattr("__doc__", "Residue identities from the built-in chemical dictionary, including amino acids, nucleotides, water and non-polymer entries.")?;
    let code_map = PyDict::new(py);
    for (key, name) in names.iter() {
        code_map.set_item(key, codes.getattr(name.extract::<String>()?.as_str())?)?;
    }
    // These are the original source vocabulary MAP keys (not callable aliases).
    for (key, name) in [
        ("TRY", "TRP"),
        ("WAT", "HOH"),
        ("H2O", "HOH"),
        ("+A", "DA"),
        ("+C", "DC"),
        ("+G", "DG"),
        ("+I", "DI"),
        ("+T", "DT"),
        ("+U", "DU"),
        ("+N", "DN"),
    ] {
        code_map.set_item(key, codes.getattr(name)?)?;
    }
    module.add("ResidueCode", codes)?;
    module.add(
        "RESIDUE_CODE_MAP",
        PyMappingProxy::new(py, code_map.cast::<PyMapping>()?),
    )?;
    let members = PyDict::new(py);
    for code in 0..12 {
        let kind = kind_from_code(code)?;
        members.set_item(kind.name(), code)?;
    }
    let kinds = int_enum.call1(("ResidueInfoKind", members))?;
    kinds.setattr("__module__", "cosmolkit")?;
    kinds.setattr("__doc__", "Chemical classification of a residue dictionary entry, including peptide, nucleic-acid and non-polymer classes.")?;
    // Preserve IntEnum.name as a string property, including mapping/pickle use.
    let error = py.get_type::<ResidueCodeParseError>();
    crate::canonical_error_accessors::residue_error(error.as_any())?;
    module.add("ResidueCodeParseError", error)?;
    let kind_map = PyDict::new(py);
    for code in 0..12 {
        let name = kind_from_code(code)?.name();
        kind_map.set_item(name, kinds.getattr(name)?)?;
    }
    module.add("ResidueInfoKind", kinds)?;
    module.add(
        "RESIDUE_INFO_KIND_MAP",
        PyMappingProxy::new(py, kind_map.cast::<PyMapping>()?),
    )?;
    module.add_class::<ResidueInfo>()?;
    module.add_class::<ResidueIdentity>()?;
    module.add_function(wrap_pyfunction!(find_residue_info, module)?)?;
    module.add_function(wrap_pyfunction!(find_residue_info_index, module)?)?;
    module.add_function(wrap_pyfunction!(residue_info, module)?)?;
    module.add_function(wrap_pyfunction!(residue_info_checked, module)?)?;
    module.add_function(wrap_pyfunction!(residue_code, module)?)?;
    module.add_function(wrap_pyfunction!(expand_one_letter, module)?)?;
    module.add_function(wrap_pyfunction!(expand_one_letter_sequence, module)?)?;
    Ok(())
}
