//! SMILES writer options; serialization delegates to the public facade.
use crate::canonical_values::SmilesWriteParams;
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
/// Bit mask selecting CXSMILES/CXSMARTS annotation fields. Combine flags with bitwise OR; coordinate fields use explicit coordinate selection.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct CxSmilesFields {
    pub(crate) inner: ck::CxSmilesFields,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CxSmilesFields {
    /// Return the underlying unsigned integer bit mask.
    fn bits(&self) -> u32 {
        self.inner.bits()
    }
    /// Return whether every bit in the supplied flag value is included in this mask.
    fn contains(&self, other: &Self) -> bool {
        self.inner.contains(other.inner)
    }
    fn __or__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }
    /// CxSmilesFields value selecting none.
    #[classattr]
    #[pyo3(name = "NONE")]
    fn flag_none() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::NONE,
        }
    }
    /// CxSmilesFields value selecting atom labels.
    #[classattr]
    #[pyo3(name = "ATOM_LABELS")]
    fn flag_atom_labels() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ATOM_LABELS,
        }
    }
    /// CxSmilesFields value selecting molfile values.
    #[classattr]
    #[pyo3(name = "MOLFILE_VALUES")]
    fn flag_molfile_values() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::MOLFILE_VALUES,
        }
    }
    /// CxSmilesFields value selecting coords.
    #[classattr]
    #[pyo3(name = "COORDS")]
    fn flag_coords() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::COORDS,
        }
    }
    /// CxSmilesFields value selecting radicals.
    #[classattr]
    #[pyo3(name = "RADICALS")]
    fn flag_radicals() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::RADICALS,
        }
    }
    /// CxSmilesFields value selecting atom props.
    #[classattr]
    #[pyo3(name = "ATOM_PROPS")]
    fn flag_atom_props() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ATOM_PROPS,
        }
    }
    /// CxSmilesFields value selecting linknodes.
    #[classattr]
    #[pyo3(name = "LINKNODES")]
    fn flag_linknodes() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::LINKNODES,
        }
    }
    /// CxSmilesFields value selecting enhanced stereo.
    #[classattr]
    #[pyo3(name = "ENHANCED_STEREO")]
    fn flag_enhanced_stereo() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ENHANCED_STEREO,
        }
    }
    /// CxSmilesFields value selecting sgroups.
    #[classattr]
    #[pyo3(name = "SGROUPS")]
    fn flag_sgroups() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::SGROUPS,
        }
    }
    /// CxSmilesFields value selecting polymer.
    #[classattr]
    #[pyo3(name = "POLYMER")]
    fn flag_polymer() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::POLYMER,
        }
    }
    /// CxSmilesFields value selecting bond cfg.
    #[classattr]
    #[pyo3(name = "BOND_CFG")]
    fn flag_bond_cfg() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::BOND_CFG,
        }
    }
    /// CxSmilesFields value selecting bond atropisomer.
    #[classattr]
    #[pyo3(name = "BOND_ATROPISOMER")]
    fn flag_bond_atropisomer() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::BOND_ATROPISOMER,
        }
    }
    /// CxSmilesFields value selecting coordinate bonds.
    #[classattr]
    #[pyo3(name = "COORDINATE_BONDS")]
    fn flag_coordinate_bonds() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::COORDINATE_BONDS,
        }
    }
    /// CxSmilesFields value selecting hydrogen bonds.
    #[classattr]
    #[pyo3(name = "HYDROGEN_BONDS")]
    fn flag_hydrogen_bonds() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::HYDROGEN_BONDS,
        }
    }
    /// CxSmilesFields value selecting zero bonds.
    #[classattr]
    #[pyo3(name = "ZERO_BONDS")]
    fn flag_zero_bonds() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ZERO_BONDS,
        }
    }
    /// CxSmilesFields value selecting all.
    #[classattr]
    #[pyo3(name = "ALL")]
    fn flag_all() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ALL,
        }
    }
    /// CxSmilesFields value selecting all but coords.
    #[classattr]
    #[pyo3(name = "ALL_BUT_COORDS")]
    fn flag_all_but_coords() -> CxSmilesFields {
        Self {
            inner: ck::CxSmilesFields::ALL_BUT_COORDS,
        }
    }
}
/// Selects which stored coordinate set a CXSMILES export uses.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct CxCoordinateSelection {
    inner: ck::CxCoordinateSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CxCoordinateSelection {
    /// Construct a coordinate selection for automatic coordinate resolution, rejecting an ambiguous selection.
    #[staticmethod]
    fn auto() -> Self {
        Self {
            inner: ck::CxCoordinateSelection::Auto,
        }
    }
    /// Construct a coordinate selection for the separate stored 2D conformer.
    #[staticmethod]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::TwoD { id },
        }
    }
    /// Construct a coordinate selection for the stored 3D conformer with the supplied ID.
    #[staticmethod]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::CxCoordinateSelection::ThreeD { id },
        }
    }
}
/// Writable configuration for CXSMILES serialization.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct CxSmilesWriteParams {
    pub(crate) inner: ck::CxSmilesWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CxSmilesWriteParams {
    /// Configure CXSMILES serialization; omitted fields use the defaults shown in the signature.
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
    /// Nested SMILES writer configuration.
    #[getter]
    fn smiles(&self) -> SmilesWriteParams {
        SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
    /// CX annotation field-selection bit mask.
    #[getter]
    fn fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.fields,
        }
    }
    /// Explicit selection of stored 2D coordinates, a 3D conformer, or automatic resolution.
    #[getter]
    fn coordinate_selection(&self) -> CxCoordinateSelection {
        CxCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
}
/// Source options accepted by RDKit's random-SMILES vector writer.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct RandomSmilesWriteParams {
    pub(crate) inner: ck::RandomSmilesWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl RandomSmilesWriteParams {
    /// Construct a RandomSmilesWriteParams value from the supplied inputs.
    #[new]
    #[pyo3(signature = (*, isomeric_smiles=true, kekule=false, all_bonds_explicit=false, all_hydrogens_explicit=false))]
    fn new(
        isomeric_smiles: bool,
        kekule: bool,
        all_bonds_explicit: bool,
        all_hydrogens_explicit: bool,
    ) -> Self {
        Self {
            inner: ck::RandomSmilesWriteParams {
                isomeric_smiles,
                kekule,
                all_bonds_explicit,
                all_hydrogens_explicit,
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
}
/// Original-index fragment selection; symbol arrays are indexed by the full molecule.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct FragmentSmilesWriteParams {
    pub(crate) inner: ck::FragmentSmilesWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FragmentSmilesWriteParams {
    /// Construct a FragmentSmilesWriteParams value from the supplied inputs.
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
    /// Zero-based atom indices to include in fragment SMILES output.
    #[getter]
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms.iter().map(|id| id.index()).collect()
    }
    /// Zero-based bond indices to include, or None to include the bonds between selected atoms.
    #[getter]
    fn bonds(&self) -> Option<Vec<usize>> {
        self.inner
            .bonds
            .as_ref()
            .map(|rows| rows.iter().map(|id| id.index()).collect())
    }
    /// Optional atom-indexed output symbol overrides for fragment serialization.
    #[getter]
    fn atom_symbols(&self) -> Option<Vec<String>> {
        self.inner.atom_symbols.clone()
    }
    /// Optional bond-indexed output symbol overrides for fragment serialization.
    #[getter]
    fn bond_symbols(&self) -> Option<Vec<String>> {
        self.inner.bond_symbols.clone()
    }
    /// Nested SMILES writer configuration.
    #[getter]
    fn smiles(&self) -> SmilesWriteParams {
        SmilesWriteParams {
            inner: self.inner.smiles,
        }
    }
}
/// Fragment selection with explicit CX fields and dimension-scoped coordinates.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct FragmentCxSmilesWriteParams {
    pub(crate) inner: ck::FragmentCxSmilesWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FragmentCxSmilesWriteParams {
    /// Construct a FragmentCxSmilesWriteParams value from the supplied inputs.
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
    /// Zero-based atom indices to include in fragment CXSMILES output.
    #[getter]
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms.iter().map(|id| id.index()).collect()
    }
    /// Zero-based bond indices to include, or None to include the bonds between selected atoms.
    #[getter]
    fn bonds(&self) -> Option<Vec<usize>> {
        self.inner
            .bonds
            .as_ref()
            .map(|rows| rows.iter().map(|id| id.index()).collect())
    }
    /// Optional atom-indexed output symbol overrides for fragment serialization.
    #[getter]
    fn atom_symbols(&self) -> Option<Vec<String>> {
        self.inner.atom_symbols.clone()
    }
    /// Optional bond-indexed output symbol overrides for fragment serialization.
    #[getter]
    fn bond_symbols(&self) -> Option<Vec<String>> {
        self.inner.bond_symbols.clone()
    }
    /// Nested CXSMILES writer configuration.
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
