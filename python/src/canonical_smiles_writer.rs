//! SMILES writer options; serialization delegates to the public facade.
use crate::canonical_values::SmilesWriteParams;
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct CxSmilesFields {
    pub(crate) inner: ck::CxSmilesFields,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CxSmilesFields {
    fn bits(&self) -> u32 {
        self.inner.bits()
    }
    fn contains(&self, other: &Self) -> bool {
        self.inner.contains(other.inner)
    }
    fn __or__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }
    #[classattr]
    #[pyo3(name = "NONE")]
    fn flag_none() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::NONE,
        }
    }
    #[classattr]
    #[pyo3(name = "ATOM_LABELS")]
    fn flag_atom_labels() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ATOM_LABELS,
        }
    }
    #[classattr]
    #[pyo3(name = "MOLFILE_VALUES")]
    fn flag_molfile_values() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::MOLFILE_VALUES,
        }
    }
    #[classattr]
    #[pyo3(name = "COORDS")]
    fn flag_coords() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::COORDS,
        }
    }
    #[classattr]
    #[pyo3(name = "RADICALS")]
    fn flag_radicals() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::RADICALS,
        }
    }
    #[classattr]
    #[pyo3(name = "ATOM_PROPS")]
    fn flag_atom_props() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ATOM_PROPS,
        }
    }
    #[classattr]
    #[pyo3(name = "LINKNODES")]
    fn flag_linknodes() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::LINKNODES,
        }
    }
    #[classattr]
    #[pyo3(name = "ENHANCED_STEREO")]
    fn flag_enhanced_stereo() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ENHANCED_STEREO,
        }
    }
    #[classattr]
    #[pyo3(name = "SGROUPS")]
    fn flag_sgroups() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::SGROUPS,
        }
    }
    #[classattr]
    #[pyo3(name = "POLYMER")]
    fn flag_polymer() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::POLYMER,
        }
    }
    #[classattr]
    #[pyo3(name = "BOND_CFG")]
    fn flag_bond_cfg() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::BOND_CFG,
        }
    }
    #[classattr]
    #[pyo3(name = "BOND_ATROPISOMER")]
    fn flag_bond_atropisomer() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::BOND_ATROPISOMER,
        }
    }
    #[classattr]
    #[pyo3(name = "COORDINATE_BONDS")]
    fn flag_coordinate_bonds() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::COORDINATE_BONDS,
        }
    }
    #[classattr]
    #[pyo3(name = "HYDROGEN_BONDS")]
    fn flag_hydrogen_bonds() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::HYDROGEN_BONDS,
        }
    }
    #[classattr]
    #[pyo3(name = "ZERO_BONDS")]
    fn flag_zero_bonds() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ZERO_BONDS,
        }
    }
    #[classattr]
    #[pyo3(name = "ALL")]
    fn flag_all() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ALL,
        }
    }
    #[classattr]
    #[pyo3(name = "ALL_BUT_COORDS")]
    fn flag_all_but_coords() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ALL_BUT_COORDS,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct CxCoordinateSelection {
    inner: ck::CxCoordinateSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CxCoordinateSelection {
    #[staticmethod]
    fn auto() -> Self {
        Self {
            inner: ck::CxCoordinateSelection::Auto,
        }
    }
    #[staticmethod]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::TwoD { id },
        }
    }
    #[staticmethod]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::ThreeD { id },
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct CxSmilesWriteParams {
    pub(crate) inner: ck::CxSmilesWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CxSmilesWriteParams {
    #[new]
    #[pyo3(signature = (*, smiles=None, fields=None, coordinate_selection=None))]
    fn new(
        smiles: Option<&SmilesWriteParams>,
        fields: Option<&CxSmilesFields>,
        coordinate_selection: Option<&CxCoordinateSelection>,
    ) -> Self {
        let defaults = ck::CxSmilesWriteParams::default();
        Self {
            inner: ck::CxSmilesWriteParams {
                smiles: smiles.map_or(defaults.smiles, |value| value.inner),
                fields: fields.map_or(defaults.fields, |value| value.inner),
                coordinate_selection: coordinate_selection
                    .map_or(defaults.coordinate_selection, |value| value.inner),
            },
        }
    }
    #[getter]
    fn smiles(&self) -> SmilesWriteParams {
        SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
    #[getter]
    fn fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.fields,
        }
    }
    #[getter]
    fn coordinate_selection(&self) -> CxCoordinateSelection {
        CxCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct RandomSmilesWriteParams {
    pub(crate) inner: ck::RandomSmilesWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl RandomSmilesWriteParams {
    #[new]
    #[pyo3(signature = (*, do_isomeric_smiles=true, do_kekule=false, all_bonds_explicit=false, all_hydrogens_explicit=false))]
    fn new(
        do_isomeric_smiles: bool,
        do_kekule: bool,
        all_bonds_explicit: bool,
        all_hydrogens_explicit: bool,
    ) -> Self {
        Self {
            inner: ck::RandomSmilesWriteParams {
                do_isomeric_smiles,
                do_kekule,
                all_bonds_explicit,
                all_hydrogens_explicit,
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
    fn all_bonds_explicit(&self) -> bool {
        self.inner.all_bonds_explicit
    }
    #[getter]
    fn all_hydrogens_explicit(&self) -> bool {
        self.inner.all_hydrogens_explicit
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FragmentSmilesWriteParams {
    pub(crate) inner: ck::FragmentSmilesWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FragmentSmilesWriteParams {
    #[new]
    #[pyo3(signature = (atoms, *, smiles=None, bonds=None, atom_symbols=None, bond_symbols=None))]
    fn new(
        atoms: Vec<usize>,
        smiles: Option<&SmilesWriteParams>,
        bonds: Option<Vec<usize>>,
        atom_symbols: Option<Vec<String>>,
        bond_symbols: Option<Vec<String>>,
    ) -> Self {
        Self {
            inner: ck::FragmentSmilesWriteParams {
                smiles: smiles.map_or_else(ck::SmilesWriteParams::default, |value| value.inner),
                atoms: atoms.into_iter().map(ck::AtomId::new).collect(),
                bonds: bonds.map(|rows| rows.into_iter().map(ck::BondId::new).collect()),
                atom_symbols,
                bond_symbols,
            },
        }
    }
    #[getter]
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms.iter().map(|id| id.index()).collect()
    }
    #[getter]
    fn bonds(&self) -> Option<Vec<usize>> {
        self.inner
            .bonds
            .as_ref()
            .map(|rows| rows.iter().map(|id| id.index()).collect())
    }
    #[getter]
    fn atom_symbols(&self) -> Option<Vec<String>> {
        self.inner.atom_symbols.clone()
    }
    #[getter]
    fn bond_symbols(&self) -> Option<Vec<String>> {
        self.inner.bond_symbols.clone()
    }
    #[getter]
    fn smiles(&self) -> SmilesWriteParams {
        SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FragmentCxSmilesWriteParams {
    pub(crate) inner: ck::FragmentCxSmilesWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FragmentCxSmilesWriteParams {
    #[new]
    #[pyo3(signature = (atoms, *, cx=None, bonds=None, atom_symbols=None, bond_symbols=None))]
    fn new(
        atoms: Vec<usize>,
        cx: Option<&CxSmilesWriteParams>,
        bonds: Option<Vec<usize>>,
        atom_symbols: Option<Vec<String>>,
        bond_symbols: Option<Vec<String>>,
    ) -> Self {
        Self {
            inner: ck::FragmentCxSmilesWriteParams {
                cx: cx.map_or_else(ck::CxSmilesWriteParams::default, |value| value.inner),
                atoms: atoms.into_iter().map(ck::AtomId::new).collect(),
                bonds: bonds.map(|rows| rows.into_iter().map(ck::BondId::new).collect()),
                atom_symbols,
                bond_symbols,
            },
        }
    }
    #[getter]
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms.iter().map(|id| id.index()).collect()
    }
    #[getter]
    fn bonds(&self) -> Option<Vec<usize>> {
        self.inner
            .bonds
            .as_ref()
            .map(|rows| rows.iter().map(|id| id.index()).collect())
    }
    #[getter]
    fn atom_symbols(&self) -> Option<Vec<String>> {
        self.inner.atom_symbols.clone()
    }
    #[getter]
    fn bond_symbols(&self) -> Option<Vec<String>> {
        self.inner.bond_symbols.clone()
    }
    #[getter]
    fn cx(&self) -> CxSmilesWriteParams {
        CxSmilesWriteParams {
            inner: self.inner.cx,
        }
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<CxSmilesFields>()?;
    module.add_class::<CxCoordinateSelection>()?;
    module.add_class::<CxSmilesWriteParams>()?;
    module.add_class::<RandomSmilesWriteParams>()?;
    module.add_class::<FragmentSmilesWriteParams>()?;
    module.add_class::<FragmentCxSmilesWriteParams>()?;
    Ok(())
}
