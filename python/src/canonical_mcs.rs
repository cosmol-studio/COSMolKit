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

/// MCS atom comparison rule: Any, element identity, isotope identity, or any non-hydrogen atom.
///
/// Declared values: ``Any``, ``Elements``, ``Isotopes``, ``AnyHeavyAtom``.
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
/// MCS bond comparison rule: Any, bond order, or exact bond-order matching.
///
/// Declared values: ``Any``, ``Order``, ``OrderExact``.
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

/// Writable configuration for MCS atom comparison.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct McsAtomCompareParameters {
    inner: ck::McsAtomCompareParameters,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsAtomCompareParameters {
    /// Configure MCS atom comparison; omitted fields use the defaults shown in the signature.
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
    /// Whether MCS atom matching requires equal valence.
    #[getter]
    fn match_valences(&self) -> bool {
        self.inner.match_valences
    }
    /// Whether MCS atom matching considers chiral tags.
    #[getter]
    fn match_chiral_tag(&self) -> bool {
        self.inner.match_chiral_tag
    }
    /// Whether MCS atom matching requires equal formal charge.
    #[getter]
    fn match_formal_charge(&self) -> bool {
        self.inner.match_formal_charge
    }
    /// Whether ring atoms/bonds may match only other ring atoms/bonds.
    #[getter]
    fn ring_matches_ring_only(&self) -> bool {
        self.inner.ring_matches_ring_only
    }
    /// Whether incomplete ring fragments are excluded from the MCS.
    #[getter]
    fn complete_rings_only(&self) -> bool {
        self.inner.complete_rings_only
    }
    /// Whether MCS atom matching requires equal isotope labels.
    #[getter]
    fn match_isotope(&self) -> bool {
        self.inner.match_isotope
    }
    /// Maximum permitted distance between matched atoms in angstroms; a negative value disables coordinate-distance matching.
    #[getter]
    fn max_distance(&self) -> f64 {
        self.inner.max_distance
    }
    fn __repr__(&self) -> String {
        format!("McsAtomCompareParameters({:?})", self.inner)
    }
}

/// Writable configuration for MCS bond comparison.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct McsBondCompareParameters {
    inner: ck::McsBondCompareParameters,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsBondCompareParameters {
    /// Configure MCS bond comparison; omitted fields use the defaults shown in the signature.
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
    /// Whether ring atoms/bonds may match only other ring atoms/bonds.
    #[getter]
    fn ring_matches_ring_only(&self) -> bool {
        self.inner.ring_matches_ring_only
    }
    /// Whether incomplete ring fragments are excluded from the MCS.
    #[getter]
    fn complete_rings_only(&self) -> bool {
        self.inner.complete_rings_only
    }
    /// Whether MCS matching considers ring-fusion relationships.
    #[getter]
    fn match_fused_rings(&self) -> bool {
        self.inner.match_fused_rings
    }
    /// Whether MCS matching requires identical ring-fusion relationships.
    #[getter]
    fn match_fused_rings_strict(&self) -> bool {
        self.inner.match_fused_rings_strict
    }
    /// Whether MCS bond matching considers bond stereochemistry.
    #[getter]
    fn match_stereo(&self) -> bool {
        self.inner.match_stereo
    }
    fn __repr__(&self) -> String {
        format!("McsBondCompareParameters({:?})", self.inner)
    }
}

/// Writable configuration for maximum-common-substructure search.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct McsParameters {
    inner: ck::McsParameters,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsParameters {
    /// Configure maximum-common-substructure search; omitted fields use the defaults shown in the signature.
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
    /// Whether all degenerate maximum common substructures are retained.
    #[getter]
    fn store_all(&self) -> bool {
        self.inner.store_all
    }
    /// Whether MCS search maximizes bond count rather than atom count.
    #[getter]
    fn maximize_bonds(&self) -> bool {
        self.inner.maximize_bonds
    }
    /// Minimum fraction of input molecules that must contain the common substructure.
    #[getter]
    fn threshold(&self) -> f64 {
        self.inner.threshold
    }
    /// Time limit in seconds for the corresponding search or embedding operation.
    #[getter]
    fn timeout(&self) -> u32 {
        self.inner.timeout
    }
    /// Whether diagnostic output is enabled.
    #[getter]
    fn verbose(&self) -> bool {
        self.inner.verbose
    }
    /// Nested MCS atom-comparison configuration; edits update this parameter object.
    #[getter]
    fn atom_compare_parameters(&self) -> McsAtomCompareParameters {
        McsAtomCompareParameters {
            inner: self.inner.atom_compare_parameters.clone(),
        }
    }
    /// Nested MCS bond-comparison configuration; edits update this parameter object.
    #[getter]
    fn bond_compare_parameters(&self) -> McsBondCompareParameters {
        McsBondCompareParameters {
            inner: self.inner.bond_compare_parameters.clone(),
        }
    }
    /// Atom comparison rule used by maximum-common-substructure search.
    #[getter]
    fn atom_comparator(&self) -> McsAtomComparator {
        McsAtomComparator::from_core(self.inner.atom_comparator)
    }
    /// Bond comparison rule used by maximum-common-substructure search.
    #[getter]
    fn bond_comparator(&self) -> McsBondComparator {
        McsBondComparator::from_core(self.inner.bond_comparator)
    }
    /// Optional SMARTS pattern used to seed MCS search.
    #[getter]
    fn initial_seed(&self) -> String {
        self.inner.initial_seed.clone()
    }
    fn __repr__(&self) -> String {
        format!("McsParameters({:?})", self.inner)
    }
}

pyo3::create_exception!(
    cosmolkit,
    McsError,
    PyValueError,
    "Maximum-common-substructure search failed for the supplied molecules or parameters."
);

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

/// Maximum-common-substructure result with its query and retained alternative queries as independent snapshots.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct McsResult {
    inner: ck::McsResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl McsResult {
    /// Return the QueryGraph represented by this compiled query.
    #[getter]
    fn query(&self) -> Option<QueryGraph> {
        self.inner.query.clone().map(|inner| QueryGraph { inner })
    }
    /// Number of atoms in the maximum common substructure.
    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count
    }
    /// Number of bonds in the maximum common substructure.
    #[getter]
    fn bond_count(&self) -> usize {
        self.inner.bond_count
    }
    /// Whether MCS search finished without reaching its timeout.
    #[getter]
    fn completed(&self) -> bool {
        self.inner.completed
    }
    /// SMARTS text defining the query or scoring pattern.
    #[getter]
    fn smarts(&self, py: Python<'_>) -> PyResult<String> {
        crate::canonical_sdf::decode_source_text(py, &self.inner.smarts)
    }
    /// All retained degenerate MCS queries when store_all is enabled.
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
/// Find an MCS with explicit options. ``timeout`` is in seconds; an interrupted
/// search returns its best partial result with ``completed == false``.
///
/// Options requiring ring or valence state use only valid existing assignments;
/// absent state produces the owner's typed error, never an implicit write-back.
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
