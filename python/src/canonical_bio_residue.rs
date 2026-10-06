//! Residue vocabulary and immutable result projections of the public BIO owner.

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

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ResidueInfo {
    pub(crate) inner: ck::ResidueInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ResidueInfo {
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueCode", i64::from(self.inner.code.as_u16()))
    }
    fn name(&self) -> &'static str {
        self.inner.name
    }
    #[gen_stub(override_return_type(type_repr = "ResidueInfoKind"))]
    fn kind<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueInfoKind", self.inner.kind as i64)
    }
    fn kind_name(&self) -> &'static str {
        self.inner.kind.name()
    }
    fn linking_type(&self) -> u8 {
        self.inner.linking_type
    }
    fn one_letter_code(&self) -> String {
        self.inner.one_letter_code.to_string()
    }
    fn hydrogen_count(&self) -> u8 {
        self.inner.hydrogen_count
    }
    fn weight(&self) -> f32 {
        self.inner.weight
    }
    fn found(&self) -> bool {
        self.inner.found()
    }
    fn is_water(&self) -> bool {
        self.inner.is_water()
    }
    fn is_dna(&self) -> bool {
        self.inner.is_dna()
    }
    fn is_rna(&self) -> bool {
        self.inner.is_rna()
    }
    fn is_nucleic_acid(&self) -> bool {
        self.inner.is_nucleic_acid()
    }
    fn is_amino_acid(&self) -> bool {
        self.inner.is_amino_acid()
    }
    fn is_buffer_or_water(&self) -> bool {
        self.inner.is_buffer_or_water()
    }
    fn is_standard(&self) -> bool {
        self.inner.is_standard()
    }
    fn fasta_code(&self) -> String {
        self.inner.fasta_code().to_string()
    }
    fn canonical_one_letter_code(&self) -> Option<String> {
        self.inner
            .canonical_one_letter_code()
            .map(|c| c.to_string())
    }
    #[gen_stub(override_return_type(type_repr = "typing.Optional[ResidueCode]"))]
    fn parent_standard_code<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.inner
            .parent_standard_code()
            .map(|c| enum_member(py, "ResidueCode", i64::from(c.as_u16())))
            .transpose()
    }
    fn is_modified_amino_acid(&self) -> bool {
        self.inner.is_modified_amino_acid()
    }
    fn is_peptide_linking(&self) -> bool {
        self.inner.is_peptide_linking()
    }
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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn find_residue_info(name: &str) -> ResidueInfo {
    ResidueInfo {
        inner: ck::find_residue_info(name),
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn find_residue_info_index(name: &str) -> usize {
    ck::find_residue_info_index(name)
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn residue_info(index: usize) -> PyResult<ResidueInfo> {
    ck::residue_info_checked(index)
        .map(|inner| ResidueInfo { inner })
        .ok_or_else(|| PyIndexError::new_err(format!("residue info index {index} out of range")))
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn residue_info_checked(index: usize) -> Option<ResidueInfo> {
    ck::residue_info_checked(index).map(|inner| ResidueInfo { inner })
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ResidueIdentity {
    inner: ck::ResidueIdentity,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ResidueIdentity {
    #[staticmethod]
    fn new(name: String) -> Self {
        Self {
            inner: ck::ResidueIdentity::new(name),
        }
    }
    fn name(&self) -> &str {
        self.inner.name()
    }
    #[gen_stub(override_return_type(type_repr = "ResidueCode"))]
    fn code<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ResidueCode", i64::from(self.inner.code().as_u16()))
    }
    fn info(&self) -> ResidueInfo {
        ResidueInfo {
            inner: self.inner.info(),
        }
    }
    fn is_tabulated(&self) -> bool {
        self.inner.is_tabulated()
    }
}

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

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn expand_one_letter(
    code: &str,
    #[gen_stub(override_type(type_repr = "ResidueInfoKind"))] kind: i64,
) -> PyResult<Option<String>> {
    let mut chars = code.chars();
    let c = chars
        .next()
        .filter(|_| chars.next().is_none())
        .ok_or_else(|| PyValueError::new_err("code must contain exactly one character"))?;
    Ok(ck::expand_one_letter(c, kind_from_code(kind)?).map(str::to_owned))
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn expand_one_letter_sequence(
    seq: &str,
    #[gen_stub(override_type(type_repr = "ResidueInfoKind"))] kind: i64,
) -> PyResult<Vec<String>> {
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
