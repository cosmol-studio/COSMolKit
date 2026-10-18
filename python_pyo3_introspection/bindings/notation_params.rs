use crate::chemical_params::SmilesWriteParams;
use crate::parameters::StringMap;
use crate::property_values::decode_source_text;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

macro_rules! cx_fields {
    ($($flag:ident),+ $(,)?) => {
        #[pyclass(module="fixed_projection._native", frozen, from_py_object, eq)]
        #[derive(Clone, Copy, PartialEq)]
        pub struct CxSmilesFields { pub(crate) inner: ck::CxSmilesFields }
        #[pymethods]
        impl CxSmilesFields {
            $(#[classattr] #[allow(non_snake_case)] fn $flag() -> Self { Self { inner: ck::CxSmilesFields::$flag } })+
            fn bits(&self) -> u32 { self.inner.bits() }
            fn contains(&self, other: &Self) -> bool { self.inner.contains(other.inner) }
            fn __or__(&self, other: &Self) -> Self { Self { inner: self.inner | other.inner } }
            fn __index__(&self) -> u32 { self.inner.bits() }
        }
        impl CxSmilesFieldsInput {
            pub(crate) fn into_rust(self) -> PyResult<ck::CxSmilesFields> {
                match self {
                    Self::Fields(value) => Ok(value.inner),
                    Self::Name(value) => {
                        $(if value == stringify!($flag).to_ascii_lowercase() { return Ok(ck::CxSmilesFields::$flag); })+
                        Err(pyo3::exceptions::PyValueError::new_err(format!("unknown CxSmilesFields: {value}")))
                    }
                }
            }
        }
    };
}
#[derive(FromPyObject, IntoPyObject)]
pub enum CxSmilesFieldsInput {
    Fields(CxSmilesFields),
    Name(String),
}
cx_fields! { NONE, ATOM_LABELS, MOLFILE_VALUES, COORDS, RADICALS, ATOM_PROPS,
LINKNODES, ENHANCED_STEREO, SGROUPS, POLYMER, BOND_CFG, BOND_ATROPISOMER,
COORDINATE_BONDS, HYDROGEN_BONDS, ZERO_BONDS, ALL, ALL_BUT_COORDS }

#[pyclass(module = "fixed_projection._native", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
pub struct CxCoordinateSelection {
    pub(crate) inner: ck::CxCoordinateSelection,
}
#[pymethods]
impl CxCoordinateSelection {
    #[classattr]
    #[pyo3(name = "Auto")]
    fn auto() -> Self {
        Self {
            inner: ck::CxCoordinateSelection::Auto,
        }
    }
    #[staticmethod]
    #[pyo3(name = "TwoD")]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::TwoD { id },
        }
    }
    #[staticmethod]
    #[pyo3(name = "ThreeD")]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::ThreeD { id },
        }
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        match self.inner {
            ck::CxCoordinateSelection::Auto => None,
            ck::CxCoordinateSelection::TwoD { id } | ck::CxCoordinateSelection::ThreeD { id } => {
                Some(id)
            }
        }
    }
}
#[derive(FromPyObject, IntoPyObject)]
pub enum CxCoordinateSelectionInput {
    Selection(CxCoordinateSelection),
    Name(String),
}
impl CxCoordinateSelectionInput {
    pub(crate) fn into_rust(self) -> PyResult<ck::CxCoordinateSelection> {
        match self {
            Self::Selection(value) => Ok(value.inner),
            Self::Name(value) if value == "auto" => Ok(ck::CxCoordinateSelection::Auto),
            Self::Name(value) => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown coordinate selection {value}; use TwoD(id) or ThreeD(id) for explicit conformer IDs"
            ))),
        }
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct CxSmilesWriteParams {
    pub(crate) inner: ck::CxSmilesWriteParams,
}
#[pymethods]
impl CxSmilesWriteParams {
    #[new]
    #[pyo3(signature=(*, smiles=None, fields=None, coordinate_selection=None))]
    fn new(
        smiles: Option<&SmilesWriteParams>,
        fields: Option<CxSmilesFieldsInput>,
        coordinate_selection: Option<CxCoordinateSelectionInput>,
    ) -> PyResult<Self> {
        let defaults = ck::CxSmilesWriteParams::default();
        Ok(Self {
            inner: ck::CxSmilesWriteParams {
                smiles: smiles.map(|value| value.inner).unwrap_or(defaults.smiles),
                fields: fields
                    .map(CxSmilesFieldsInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.fields),
                coordinate_selection: coordinate_selection
                    .map(CxCoordinateSelectionInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.coordinate_selection),
            },
        })
    }
    #[getter]
    fn smiles(&self) -> SmilesWriteParams {
        SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
    #[setter]
    fn set_smiles(&mut self, value: Option<&SmilesWriteParams>) {
        self.inner.smiles = value
            .map(|value| value.inner)
            .unwrap_or_else(|| ck::CxSmilesWriteParams::default().smiles);
    }
    #[getter]
    fn fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.fields,
        }
    }
    #[setter]
    fn set_fields(&mut self, value: Option<CxSmilesFieldsInput>) -> PyResult<()> {
        self.inner.fields = value
            .map(CxSmilesFieldsInput::into_rust)
            .transpose()?
            .unwrap_or_else(|| ck::CxSmilesWriteParams::default().fields);
        Ok(())
    }
    #[getter]
    fn coordinate_selection(&self) -> CxCoordinateSelection {
        CxCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
    #[setter]
    fn set_coordinate_selection(
        &mut self,
        value: Option<CxCoordinateSelectionInput>,
    ) -> PyResult<()> {
        self.inner.coordinate_selection = value
            .map(CxCoordinateSelectionInput::into_rust)
            .transpose()?
            .unwrap_or_else(|| ck::CxSmilesWriteParams::default().coordinate_selection);
        Ok(())
    }
}

macro_rules! fragment_parameters {
    ($name:ident, $field:ident: $nested:ident) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct $name { pub(crate) inner: ck::$name }
        paste::paste! {
            #[pymethods]
            impl $name {
                #[new]
                #[pyo3(signature=(*, atoms, $field=None, bonds=None, atom_symbols=None, bond_symbols=None))]
                fn new(atoms: Vec<usize>, $field: Option<&$nested>, bonds: Option<Vec<usize>>, atom_symbols: Option<Vec<String>>, bond_symbols: Option<Vec<String>>) -> Self {
                    Self { inner: ck::$name {
                        atoms: atoms.into_iter().map(ck::AtomId::new).collect(),
                        $field: $field.map(|value| value.inner.clone()).unwrap_or_else(|| ck::$name::default().$field),
                        bonds: bonds.map(|values| values.into_iter().map(ck::BondId::new).collect()), atom_symbols, bond_symbols,
                    } }
                }
                #[getter] fn atoms(&self) -> Vec<usize> { self.inner.atoms.iter().copied().map(ck::AtomId::index).collect() }
                #[setter] fn set_atoms(&mut self, value: Vec<usize>) { self.inner.atoms = value.into_iter().map(ck::AtomId::new).collect(); }
                #[getter] fn $field(&self) -> $nested { $nested { inner: self.inner.$field.clone() } }
                #[setter] fn [<set_ $field>](&mut self, value: Option<&$nested>) { self.inner.$field = value.map(|value| value.inner.clone()).unwrap_or_else(|| ck::$name::default().$field); }
                #[getter] fn bonds(&self) -> Option<Vec<usize>> { self.inner.bonds.as_ref().map(|values| values.iter().copied().map(ck::BondId::index).collect()) }
                #[setter] fn set_bonds(&mut self, value: Option<Vec<usize>>) { self.inner.bonds = value.map(|values| values.into_iter().map(ck::BondId::new).collect()); }
                #[getter] fn atom_symbols(&self) -> Option<Vec<String>> { self.inner.atom_symbols.clone() }
                #[setter] fn set_atom_symbols(&mut self, value: Option<Vec<String>>) { self.inner.atom_symbols = value; }
                #[getter] fn bond_symbols(&self) -> Option<Vec<String>> { self.inner.bond_symbols.clone() }
                #[setter] fn set_bond_symbols(&mut self, value: Option<Vec<String>>) { self.inner.bond_symbols = value; }
            }
        }
    };
}
fragment_parameters!(FragmentSmilesWriteParams, smiles: SmilesWriteParams);
fragment_parameters!(FragmentCxSmilesWriteParams, cx: CxSmilesWriteParams);

#[pyclass(module = "fixed_projection._native")]
pub struct SmartsParseParams {
    pub(crate) inner: ck::SmartsParseParams,
}
macro_rules! smarts_parameters {
    ($($field:ident),+ $(,)?) => {
        paste::paste! {
            #[pymethods]
            impl SmartsParseParams {
                #[new]
                #[pyo3(signature=(*, $($field=ck::SmartsParseParams::default().$field,)+ replacements=None))]
                fn new($($field: bool,)+ replacements: Option<StringMap>) -> Self {
                    Self { inner: ck::SmartsParseParams { $($field,)+ replacements: replacements.map(|values| values.0.into_iter().map(|(key,value)| (key.into(),value.into())).collect()).unwrap_or_else(|| ck::SmartsParseParams::default().replacements) } }
                }
                $(#[getter] fn $field(&self) -> bool { self.inner.$field }
                #[setter] fn [<set_ $field>](&mut self, value: bool) { self.inner.$field = value; })+
                #[getter] fn replacements(&self, py: Python<'_>) -> PyResult<StringMap> {
                    self.inner.replacements.iter().map(|(key, value)| Ok((decode_source_text(py,key)?,decode_source_text(py,value)?))).collect::<PyResult<BTreeMap<String,String>>>().map(StringMap)
                }
                #[setter] fn set_replacements(&mut self, values: Option<StringMap>) {
                    self.inner.replacements = values.map(|values| values.0.into_iter().map(|(key,value)| (key.into(),value.into())).collect()).unwrap_or_else(|| ck::SmartsParseParams::default().replacements);
                }
            }
        }
    };
}
smarts_parameters! { allow_cxsmiles, strict_cxsmiles, parse_name, merge_hs, skip_cleanup, debug_parse }
