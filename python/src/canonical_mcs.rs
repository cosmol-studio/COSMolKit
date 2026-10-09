//! Thin MCS projections. The search owner implements all chemistry.
use crate::canonical_search::QueryGraph;
use crate::drawing_binding::Molecule;
use ::cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};
use std::collections::BTreeMap;

#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum McsAtomComparator {
    Any,
    Elements,
    Isotopes,
    AnyHeavyAtom,
}
impl McsAtomComparator {
    fn core(self) -> ck::McsAtomComparator {
        match self {
            Self::Any => ck::McsAtomComparator::AtomCompareAny,
            Self::Elements => ck::McsAtomComparator::AtomCompareElements,
            Self::Isotopes => ck::McsAtomComparator::AtomCompareIsotopes,
            Self::AnyHeavyAtom => ck::McsAtomComparator::AtomCompareAnyHeavyAtom,
        }
    }
    fn from_core(value: ck::McsAtomComparator) -> Self {
        match value {
            ck::McsAtomComparator::AtomCompareAny => Self::Any,
            ck::McsAtomComparator::AtomCompareElements => Self::Elements,
            ck::McsAtomComparator::AtomCompareIsotopes => Self::Isotopes,
            ck::McsAtomComparator::AtomCompareAnyHeavyAtom => Self::AnyHeavyAtom,
        }
    }
}
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum McsBondComparator {
    Any,
    Order,
    OrderExact,
}
impl McsBondComparator {
    fn core(self) -> ck::McsBondComparator {
        match self {
            Self::Any => ck::McsBondComparator::BondCompareAny,
            Self::Order => ck::McsBondComparator::BondCompareOrder,
            Self::OrderExact => ck::McsBondComparator::BondCompareOrderExact,
        }
    }
    fn from_core(value: ck::McsBondComparator) -> Self {
        match value {
            ck::McsBondComparator::BondCompareAny => Self::Any,
            ck::McsBondComparator::BondCompareOrder => Self::Order,
            ck::McsBondComparator::BondCompareOrderExact => Self::OrderExact,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct McsAtomCompareParameters {
    inner: ck::McsAtomCompareParameters,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsAtomCompareParameters {
    #[new]
    #[pyo3(signature = (*, match_valences=false, match_chiral_tag=false, match_formal_charge=false, ring_matches_ring_only=false, complete_rings_only=false, match_isotope=false, max_distance=-1.0))]
    fn new(
        match_valences: bool,
        match_chiral_tag: bool,
        match_formal_charge: bool,
        ring_matches_ring_only: bool,
        complete_rings_only: bool,
        match_isotope: bool,
        max_distance: f64,
    ) -> Self {
        Self {
            inner: ck::McsAtomCompareParameters {
                match_valences,
                match_chiral_tag,
                match_formal_charge,
                ring_matches_ring_only,
                complete_rings_only,
                match_isotope,
                max_distance,
            },
        }
    }
    #[getter]
    fn match_valences(&self) -> bool {
        self.inner.match_valences
    }
    #[getter]
    fn match_chiral_tag(&self) -> bool {
        self.inner.match_chiral_tag
    }
    #[getter]
    fn match_formal_charge(&self) -> bool {
        self.inner.match_formal_charge
    }
    #[getter]
    fn ring_matches_ring_only(&self) -> bool {
        self.inner.ring_matches_ring_only
    }
    #[getter]
    fn complete_rings_only(&self) -> bool {
        self.inner.complete_rings_only
    }
    #[getter]
    fn match_isotope(&self) -> bool {
        self.inner.match_isotope
    }
    #[getter]
    fn max_distance(&self) -> f64 {
        self.inner.max_distance
    }
    fn __repr__(&self) -> String {
        format!("McsAtomCompareParameters({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct McsBondCompareParameters {
    inner: ck::McsBondCompareParameters,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsBondCompareParameters {
    #[new]
    #[pyo3(signature = (*, ring_matches_ring_only=false, complete_rings_only=false, match_fused_rings=false, match_fused_rings_strict=false, match_stereo=false))]
    fn new(
        ring_matches_ring_only: bool,
        complete_rings_only: bool,
        match_fused_rings: bool,
        match_fused_rings_strict: bool,
        match_stereo: bool,
    ) -> Self {
        Self {
            inner: ck::McsBondCompareParameters {
                ring_matches_ring_only,
                complete_rings_only,
                match_fused_rings,
                match_fused_rings_strict,
                match_stereo,
            },
        }
    }
    #[getter]
    fn ring_matches_ring_only(&self) -> bool {
        self.inner.ring_matches_ring_only
    }
    #[getter]
    fn complete_rings_only(&self) -> bool {
        self.inner.complete_rings_only
    }
    #[getter]
    fn match_fused_rings(&self) -> bool {
        self.inner.match_fused_rings
    }
    #[getter]
    fn match_fused_rings_strict(&self) -> bool {
        self.inner.match_fused_rings_strict
    }
    #[getter]
    fn match_stereo(&self) -> bool {
        self.inner.match_stereo
    }
    fn __repr__(&self) -> String {
        format!("McsBondCompareParameters({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct McsParameters {
    inner: ck::McsParameters,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsParameters {
    #[new]
    #[pyo3(signature = (*, store_all=false, maximize_bonds=true, threshold=1.0, timeout=0, verbose=false, atom_compare_parameters=None, bond_compare_parameters=None, atom_comparator=McsAtomComparator::Elements, bond_comparator=McsBondComparator::Order, initial_seed=String::new()))]
    fn new(
        store_all: bool,
        maximize_bonds: bool,
        threshold: f64,
        timeout: u32,
        verbose: bool,
        atom_compare_parameters: Option<&McsAtomCompareParameters>,
        bond_compare_parameters: Option<&McsBondCompareParameters>,
        atom_comparator: McsAtomComparator,
        bond_comparator: McsBondComparator,
        initial_seed: String,
    ) -> Self {
        Self {
            inner: ck::McsParameters {
                store_all,
                maximize_bonds,
                threshold,
                timeout,
                verbose,
                atom_compare_parameters: atom_compare_parameters
                    .map(|v| v.inner.clone())
                    .unwrap_or_default(),
                bond_compare_parameters: bond_compare_parameters
                    .map(|v| v.inner.clone())
                    .unwrap_or_default(),
                atom_comparator: atom_comparator.core(),
                bond_comparator: bond_comparator.core(),
                initial_seed,
            },
        }
    }
    #[getter]
    fn store_all(&self) -> bool {
        self.inner.store_all
    }
    #[getter]
    fn maximize_bonds(&self) -> bool {
        self.inner.maximize_bonds
    }
    #[getter]
    fn threshold(&self) -> f64 {
        self.inner.threshold
    }
    #[getter]
    fn timeout(&self) -> u32 {
        self.inner.timeout
    }
    #[getter]
    fn verbose(&self) -> bool {
        self.inner.verbose
    }
    #[getter]
    fn atom_compare_parameters(&self) -> McsAtomCompareParameters {
        McsAtomCompareParameters {
            inner: self.inner.atom_compare_parameters.clone(),
        }
    }
    #[getter]
    fn bond_compare_parameters(&self) -> McsBondCompareParameters {
        McsBondCompareParameters {
            inner: self.inner.bond_compare_parameters.clone(),
        }
    }
    #[getter]
    fn atom_comparator(&self) -> McsAtomComparator {
        McsAtomComparator::from_core(self.inner.atom_comparator)
    }
    #[getter]
    fn bond_comparator(&self) -> McsBondComparator {
        McsBondComparator::from_core(self.inner.bond_comparator)
    }
    #[getter]
    fn initial_seed(&self) -> String {
        self.inner.initial_seed.clone()
    }
    fn __repr__(&self) -> String {
        format!("McsParameters({:?})", self.inner)
    }
}

pyo3::create_exception!(cosmolkit, McsError, PyValueError);

fn mcs_pyerr(py: Python<'_>, source: ck::McsError) -> PyErr {
    let kind = match &source {
        ck::McsError::State(_) => "State",
        ck::McsError::Progress(_) => "Progress",
        ck::McsError::StereoOrder(_) => "StereoOrder",
        ck::McsError::TargetTableCount { .. } => "TargetTableCount",
        ck::McsError::ThresholdCountOutOfRange { .. } => "ThresholdCountOutOfRange",
        ck::McsError::MatchTableOutOfRange { .. } => "MatchTableOutOfRange",
        ck::McsError::RingMembershipMissing { .. } => "RingMembershipMissing",
        ck::McsError::RingOutOfRange { .. } => "RingOutOfRange",
        ck::McsError::MappedBondMissing { .. } => "MappedBondMissing",
        ck::McsError::TooManyNewBonds { .. } => "TooManyNewBonds",
        ck::McsError::InitialSeedParse { .. } => "InitialSeedParse",
        ck::McsError::InitialSeedMatch { .. } => "InitialSeedMatch",
        ck::McsError::InitialSeedBondMissing { .. } => "InitialSeedBondMissing",
        ck::McsError::ResultValueOutOfRange { .. } => "ResultValueOutOfRange",
        ck::McsError::ResultQueryGraph { .. } => "ResultQueryGraph",
        ck::McsError::ResultSmarts { .. } => "ResultSmarts",
        ck::McsError::SeedReconstructionBondMissing { .. } => "SeedReconstructionBondMissing",
        ck::McsError::ResultContextMissing => "ResultContextMissing",
    };
    crate::canonical_values::annotate(
        py,
        McsError::new_err(source.to_string()),
        "mcs",
        kind,
        &source,
    )
}

/// A detached query result. Query and retained alternatives are independent snapshots.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct McsResult {
    inner: ck::McsResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsResult {
    #[getter]
    fn query(&self) -> Option<QueryGraph> {
        self.inner.query.clone().map(|inner| QueryGraph { inner })
    }
    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count
    }
    #[getter]
    fn bond_count(&self) -> usize {
        self.inner.bond_count
    }
    #[getter]
    fn completed(&self) -> bool {
        self.inner.completed
    }
    #[getter]
    fn smarts(&self, py: Python<'_>) -> PyResult<String> {
        crate::canonical_sdf::decode_source_text(py, &self.inner.smarts)
    }
    #[getter]
    fn degenerate(&self, py: Python<'_>) -> PyResult<BTreeMap<String, QueryGraph>> {
        self.inner
            .degenerate
            .iter()
            .map(|(text, query)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, text)?,
                    QueryGraph {
                        inner: query.clone(),
                    },
                ))
            })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "McsResult(atom_count={}, bond_count={}, completed={})",
            self.inner.atom_count, self.inner.bond_count, self.inner.completed
        )
    }
}

/// Experimental maximum common query search; inputs are not modified.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn maximum_common_substructure(
    py: Python<'_>,
    inputs: Vec<PyRef<'_, Molecule>>,
) -> PyResult<McsResult> {
    let inputs: Vec<_> = inputs.iter().map(|mol| &mol.inner).collect();
    ck::maximum_common_substructure(&inputs)
        .map(|inner| McsResult { inner })
        .map_err(|error| mcs_pyerr(py, error))
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn maximum_common_substructure_with_params(
    py: Python<'_>,
    inputs: Vec<PyRef<'_, Molecule>>,
    params: &McsParameters,
) -> PyResult<McsResult> {
    let inputs: Vec<_> = inputs.iter().map(|mol| &mol.inner).collect();
    ck::maximum_common_substructure_with_params(&inputs, &params.inner)
        .map(|inner| McsResult { inner })
        .map_err(|error| mcs_pyerr(py, error))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<McsAtomComparator>()?;
    module.add_class::<McsBondComparator>()?;
    module.add_class::<McsAtomCompareParameters>()?;
    module.add_class::<McsBondCompareParameters>()?;
    module.add_class::<McsParameters>()?;
    module.add_class::<McsResult>()?;
    module.add("McsError", module.py().get_type::<McsError>())?;
    module.add_function(wrap_pyfunction!(maximum_common_substructure, module)?)?;
    module.add_function(wrap_pyfunction!(
        maximum_common_substructure_with_params,
        module
    )?)?;
    Ok(())
}
