use crate::chemical_values::{EnumValue, McsAtomComparator, McsBondComparator};
use crate::construction_values::ReactionCoordinateSelection;
use crate::notation_params::{CxSmilesFields, CxSmilesFieldsInput};
use crate::params::{McsAtomCompareParameters, McsBondCompareParameters};
use cosmolkit as ck;
use pyo3::prelude::*;

#[derive(FromPyObject, IntoPyObject)]
pub enum McsAtomComparatorInput {
    Enum(McsAtomComparator),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum McsBondComparatorInput {
    Enum(McsBondComparator),
    Name(String),
}
impl McsAtomComparatorInput {
    fn into_rust(self) -> PyResult<ck::McsAtomComparator> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(value) => McsAtomComparator::from_name(&value).map(Into::into),
        }
    }
}
impl McsBondComparatorInput {
    fn into_rust(self) -> PyResult<ck::McsBondComparator> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(value) => McsBondComparator::from_name(&value).map(Into::into),
        }
    }
}

macro_rules! mcs_parameters {
    ($($field:ident: $ty:ty),+ $(,)?) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct McsParameters { pub(crate) inner: ck::McsParameters }
        paste::paste! {
            #[pymethods]
            impl McsParameters {
                #[new]
                #[pyo3(signature=(*, $($field=ck::McsParameters::default().$field,)+ atom_compare_parameters=None, bond_compare_parameters=None, atom_comparator=None, bond_comparator=None))]
                pub(crate) fn new($($field: $ty,)+ atom_compare_parameters: Option<&McsAtomCompareParameters>, bond_compare_parameters: Option<&McsBondCompareParameters>, atom_comparator: Option<McsAtomComparatorInput>, bond_comparator: Option<McsBondComparatorInput>) -> PyResult<Self> {
                    let defaults = ck::McsParameters::default();
                    Ok(Self { inner: ck::McsParameters { $($field,)+
                        atom_compare_parameters: atom_compare_parameters.map(|value| value.inner.clone()).unwrap_or(defaults.atom_compare_parameters),
                        bond_compare_parameters: bond_compare_parameters.map(|value| value.inner.clone()).unwrap_or(defaults.bond_compare_parameters),
                        atom_comparator: atom_comparator.map(McsAtomComparatorInput::into_rust).transpose()?.unwrap_or(defaults.atom_comparator),
                        bond_comparator: bond_comparator.map(McsBondComparatorInput::into_rust).transpose()?.unwrap_or(defaults.bond_comparator),
                    } })
                }
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
                #[getter] fn atom_compare_parameters(&self) -> McsAtomCompareParameters { McsAtomCompareParameters { inner: self.inner.atom_compare_parameters.clone() } }
                #[setter] fn set_atom_compare_parameters(&mut self, value: Option<&McsAtomCompareParameters>) { self.inner.atom_compare_parameters = value.map(|value| value.inner.clone()).unwrap_or_else(|| ck::McsParameters::default().atom_compare_parameters); }
                #[getter] fn bond_compare_parameters(&self) -> McsBondCompareParameters { McsBondCompareParameters { inner: self.inner.bond_compare_parameters.clone() } }
                #[setter] fn set_bond_compare_parameters(&mut self, value: Option<&McsBondCompareParameters>) { self.inner.bond_compare_parameters = value.map(|value| value.inner.clone()).unwrap_or_else(|| ck::McsParameters::default().bond_compare_parameters); }
                #[getter] fn atom_comparator(&self) -> McsAtomComparator { self.inner.atom_comparator.into() }
                #[setter] fn set_atom_comparator(&mut self, value: Option<McsAtomComparatorInput>) -> PyResult<()> { self.inner.atom_comparator = value.map(McsAtomComparatorInput::into_rust).transpose()?.unwrap_or_else(|| ck::McsParameters::default().atom_comparator); Ok(()) }
                #[getter] fn bond_comparator(&self) -> McsBondComparator { self.inner.bond_comparator.into() }
                #[setter] fn set_bond_comparator(&mut self, value: Option<McsBondComparatorInput>) -> PyResult<()> { self.inner.bond_comparator = value.map(McsBondComparatorInput::into_rust).transpose()?.unwrap_or_else(|| ck::McsParameters::default().bond_comparator); Ok(()) }
            }
        }
    };
}
mcs_parameters! { store_all: bool, maximize_bonds: bool, threshold: f64,
timeout: u32, verbose: bool, initial_seed: String }

#[pyclass(module = "fixed_projection._native")]
pub struct ReactionSingleRunParams {
    pub(crate) inner: ck::ReactionSingleRunParams,
}
#[pymethods]
impl ReactionSingleRunParams {
    #[new]
    #[pyo3(signature=(*, coordinate_selection=None))]
    fn new(coordinate_selection: Option<&ReactionCoordinateSelection>) -> Self {
        Self {
            inner: ck::ReactionSingleRunParams {
                coordinate_selection: coordinate_selection
                    .map(|value| value.inner)
                    .unwrap_or_else(|| ck::ReactionSingleRunParams::default().coordinate_selection),
            },
        }
    }
    #[getter]
    fn coordinate_selection(&self) -> ReactionCoordinateSelection {
        ReactionCoordinateSelection {
            inner: self.inner.coordinate_selection,
        }
    }
    #[setter]
    fn set_coordinate_selection(&mut self, value: Option<&ReactionCoordinateSelection>) {
        self.inner.coordinate_selection = value
            .map(|value| value.inner)
            .unwrap_or_else(|| ck::ReactionSingleRunParams::default().coordinate_selection);
    }
}

macro_rules! reaction_write_parameters {
    ($($field:ident: $ty:ty),+ $(,)?) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct ReactionWriteParams { pub(crate) inner: ck::ReactionWriteParams }
        paste::paste! {
            #[pymethods]
            impl ReactionWriteParams {
                #[new]
                #[pyo3(signature=(*, $($field=ck::ReactionWriteParams::default().$field,)+ cx_fields=None, coordinate_selections=None))]
                pub(crate) fn new($($field: $ty,)+ cx_fields: Option<CxSmilesFieldsInput>, coordinate_selections: Option<Vec<ReactionCoordinateSelection>>) -> PyResult<Self> {
                    Ok(Self { inner: ck::ReactionWriteParams { $($field,)+
                        cx_fields: cx_fields.map(CxSmilesFieldsInput::into_rust).transpose()?.unwrap_or_else(|| ck::ReactionWriteParams::default().cx_fields),
                        coordinate_selections: coordinate_selections.map(|values| values.into_iter().map(|value| value.inner).collect()).unwrap_or_else(|| ck::ReactionWriteParams::default().coordinate_selections),
                    } })
                }
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
                #[getter] fn cx_fields(&self) -> CxSmilesFields { CxSmilesFields { inner: self.inner.cx_fields } }
                #[setter] fn set_cx_fields(&mut self, value: Option<CxSmilesFieldsInput>) -> PyResult<()> { self.inner.cx_fields = value.map(CxSmilesFieldsInput::into_rust).transpose()?.unwrap_or_else(|| ck::ReactionWriteParams::default().cx_fields); Ok(()) }
                #[getter] fn coordinate_selections(&self) -> Vec<ReactionCoordinateSelection> { self.inner.coordinate_selections.iter().copied().map(|inner| ReactionCoordinateSelection { inner }).collect() }
                #[setter] fn set_coordinate_selections(&mut self, value: Option<Vec<ReactionCoordinateSelection>>) { self.inner.coordinate_selections = value.map(|values| values.into_iter().map(|value| value.inner).collect()).unwrap_or_else(|| ck::ReactionWriteParams::default().coordinate_selections); }
            }
        }
    };
}
reaction_write_parameters! { canonical: bool, isomeric_smiles: bool, rooted_at_atom: Option<usize>,
include_dative_bonds: bool, include_cx: bool }
