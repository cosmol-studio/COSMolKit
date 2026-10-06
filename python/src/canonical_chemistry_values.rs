//! Canonical configuration and result values for chemistry facade projections.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct KekulizeParams {
    pub(crate) inner: ck::KekulizeParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl KekulizeParams {
    #[new]
    #[pyo3(signature = (*, mark_atoms_bonds=true, canonical=true, max_backtracks=100))]
    fn new(mark_atoms_bonds: bool, canonical: bool, max_backtracks: u32) -> Self {
        Self {
            inner: ck::KekulizeParams {
                mark_atoms_bonds,
                canonical,
                max_backtracks,
            },
        }
    }
    #[getter]
    fn mark_atoms_bonds(&self) -> bool {
        self.inner.mark_atoms_bonds
    }
    #[getter]
    fn canonical(&self) -> bool {
        self.inner.canonical
    }
    #[getter]
    fn max_backtracks(&self) -> u32 {
        self.inner.max_backtracks
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct RingSearchParams {
    pub(crate) inner: ck::RingSearchParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl RingSearchParams {
    #[new]
    #[pyo3(signature = (*, include_dative_bonds=false, include_hydrogen_bonds=false))]
    fn new(include_dative_bonds: bool, include_hydrogen_bonds: bool) -> Self {
        Self {
            inner: ck::RingSearchParams {
                include_dative_bonds,
                include_hydrogen_bonds,
            },
        }
    }
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[getter]
    fn include_hydrogen_bonds(&self) -> bool {
        self.inner.include_hydrogen_bonds
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct StructureTagParams {
    pub(crate) inner: ck::StructureTagParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl StructureTagParams {
    #[new]
    #[pyo3(signature = (*, conformer_id=-1, replace_existing_tags=true))]
    fn new(conformer_id: i32, replace_existing_tags: bool) -> Self {
        Self {
            inner: ck::StructureTagParams {
                conformer_id,
                replace_existing_tags,
            },
        }
    }
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[getter]
    fn replace_existing_tags(&self) -> bool {
        self.inner.replace_existing_tags
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct DistanceMatrixParams {
    pub(crate) inner: ck::DistanceMatrixParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl DistanceMatrixParams {
    #[new]
    #[pyo3(signature = (*, use_bond_order=false, use_atom_weights=false))]
    fn new(use_bond_order: bool, use_atom_weights: bool) -> Self {
        Self {
            inner: ck::DistanceMatrixParams {
                use_bond_order,
                use_atom_weights,
            },
        }
    }
    #[getter]
    fn use_bond_order(&self) -> bool {
        self.inner.use_bond_order
    }
    #[getter]
    fn use_atom_weights(&self) -> bool {
        self.inner.use_atom_weights
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct DistanceMatrix3dParams {
    pub(crate) inner: ck::DistanceMatrix3dParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl DistanceMatrix3dParams {
    #[new]
    #[pyo3(signature = (*, conformer_id=None, use_atom_weights=false))]
    fn new(conformer_id: Option<usize>, use_atom_weights: bool) -> Self {
        Self {
            inner: ck::DistanceMatrix3dParams {
                conformer_id,
                use_atom_weights,
            },
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id
    }
    #[getter]
    fn use_atom_weights(&self) -> bool {
        self.inner.use_atom_weights
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPositionParams {
    pub(crate) inner: ck::AtomPositionParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPositionParams {
    #[new]
    #[pyo3(signature = (*, conformer_id=None))]
    fn new(conformer_id: Option<usize>) -> Self {
        Self {
            inner: ck::AtomPositionParams { conformer_id },
        }
    }
    #[getter]
    fn conformer_id(&self) -> Option<usize> {
        self.inner.conformer_id
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct RemoveHsParams {
    pub(crate) inner: ck::RemoveHsParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl RemoveHsParams {
    #[new]
    #[pyo3(signature = (*, remove_degree_zero=false, remove_higher_degrees=false, remove_only_h_neighbors=false, remove_isotopes=false, remove_and_track_isotopes=false, remove_dummy_neighbors=false, remove_defining_bond_stereo=false, remove_with_wedged_bond=true, remove_with_query=false, remove_mapped=true, remove_in_sgroups=true, show_warnings=true, remove_nonimplicit=true, update_explicit_count=false, remove_hydrides=false, remove_nontetrahedral_neighbors=false, sanitize=true))]
    fn new(
        remove_degree_zero: bool,
        remove_higher_degrees: bool,
        remove_only_h_neighbors: bool,
        remove_isotopes: bool,
        remove_and_track_isotopes: bool,
        remove_dummy_neighbors: bool,
        remove_defining_bond_stereo: bool,
        remove_with_wedged_bond: bool,
        remove_with_query: bool,
        remove_mapped: bool,
        remove_in_sgroups: bool,
        show_warnings: bool,
        remove_nonimplicit: bool,
        update_explicit_count: bool,
        remove_hydrides: bool,
        remove_nontetrahedral_neighbors: bool,
        sanitize: bool,
    ) -> Self {
        Self {
            inner: ck::RemoveHsParams {
                remove_degree_zero,
                remove_higher_degrees,
                remove_only_h_neighbors,
                remove_isotopes,
                remove_and_track_isotopes,
                remove_dummy_neighbors,
                remove_defining_bond_stereo,
                remove_with_wedged_bond,
                remove_with_query,
                remove_mapped,
                remove_in_sgroups,
                show_warnings,
                remove_nonimplicit,
                update_explicit_count,
                remove_hydrides,
                remove_nontetrahedral_neighbors,
                sanitize,
            },
        }
    }
    #[getter]
    fn remove_degree_zero(&self) -> bool {
        self.inner.remove_degree_zero
    }
    #[getter]
    fn remove_higher_degrees(&self) -> bool {
        self.inner.remove_higher_degrees
    }
    #[getter]
    fn remove_only_h_neighbors(&self) -> bool {
        self.inner.remove_only_h_neighbors
    }
    #[getter]
    fn remove_isotopes(&self) -> bool {
        self.inner.remove_isotopes
    }
    #[getter]
    fn remove_and_track_isotopes(&self) -> bool {
        self.inner.remove_and_track_isotopes
    }
    #[getter]
    fn remove_dummy_neighbors(&self) -> bool {
        self.inner.remove_dummy_neighbors
    }
    #[getter]
    fn remove_defining_bond_stereo(&self) -> bool {
        self.inner.remove_defining_bond_stereo
    }
    #[getter]
    fn remove_with_wedged_bond(&self) -> bool {
        self.inner.remove_with_wedged_bond
    }
    #[getter]
    fn remove_with_query(&self) -> bool {
        self.inner.remove_with_query
    }
    #[getter]
    fn remove_mapped(&self) -> bool {
        self.inner.remove_mapped
    }
    #[getter]
    fn remove_in_sgroups(&self) -> bool {
        self.inner.remove_in_sgroups
    }
    #[getter]
    fn show_warnings(&self) -> bool {
        self.inner.show_warnings
    }
    #[getter]
    fn remove_nonimplicit(&self) -> bool {
        self.inner.remove_nonimplicit
    }
    #[getter]
    fn update_explicit_count(&self) -> bool {
        self.inner.update_explicit_count
    }
    #[getter]
    fn remove_hydrides(&self) -> bool {
        self.inner.remove_hydrides
    }
    #[getter]
    fn remove_nontetrahedral_neighbors(&self) -> bool {
        self.inner.remove_nontetrahedral_neighbors
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum AromaticityModel {
    Rdkit,
    Simple,
    Mdl,
    Mmff94,
    Custom,
}
impl AromaticityModel {
    fn core(self) -> ck::AromaticityModel {
        match self {
            Self::Rdkit => ck::AromaticityModel::Rdkit,
            Self::Simple => ck::AromaticityModel::Simple,
            Self::Mdl => ck::AromaticityModel::Mdl,
            Self::Mmff94 => ck::AromaticityModel::Mmff94,
            Self::Custom => ck::AromaticityModel::Custom,
        }
    }
    fn from_core(value: ck::AromaticityModel) -> Self {
        match value {
            ck::AromaticityModel::Rdkit => Self::Rdkit,
            ck::AromaticityModel::Simple => Self::Simple,
            ck::AromaticityModel::Mdl => Self::Mdl,
            ck::AromaticityModel::Mmff94 => Self::Mmff94,
            ck::AromaticityModel::Custom => Self::Custom,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AromaticityParams {
    pub(crate) inner: ck::AromaticityParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AromaticityParams {
    #[new]
    #[pyo3(signature = (model=AromaticityModel::Rdkit))]
    fn new(model: AromaticityModel) -> Self {
        Self {
            inner: ck::AromaticityParams {
                model: model.core(),
            },
        }
    }
    #[getter]
    fn model(&self) -> AromaticityModel {
        AromaticityModel::from_core(self.inner.model)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AddHsParams {
    pub(crate) inner: ck::AddHsParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AddHsParams {
    #[new]
    #[pyo3(signature = (*, explicit_only=false, add_coords=false, add_residue_info=false, skip_queries=false, only_on_atoms=None))]
    fn new(
        explicit_only: bool,
        add_coords: bool,
        add_residue_info: bool,
        skip_queries: bool,
        only_on_atoms: Option<Vec<usize>>,
    ) -> Self {
        Self {
            inner: ck::AddHsParams {
                explicit_only,
                add_coords,
                add_residue_info,
                skip_queries,
                only_on_atoms: only_on_atoms
                    .map(|rows| rows.into_iter().map(ck::AtomId::new).collect()),
            },
        }
    }
    #[getter]
    fn explicit_only(&self) -> bool {
        self.inner.explicit_only
    }
    #[getter]
    fn add_coords(&self) -> bool {
        self.inner.add_coords
    }
    #[getter]
    fn add_residue_info(&self) -> bool {
        self.inner.add_residue_info
    }
    #[getter]
    fn skip_queries(&self) -> bool {
        self.inner.skip_queries
    }
    #[getter]
    fn only_on_atoms(&self) -> Option<Vec<usize>> {
        self.inner
            .only_on_atoms
            .as_ref()
            .map(|rows| rows.iter().map(|id| id.index()).collect())
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct DenseMatrix {
    pub(crate) inner: ck::DenseMatrix,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl DenseMatrix {
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }
    fn values(&self) -> Vec<f64> {
        self.inner.values().to_vec()
    }
    fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.inner.get(row, column)
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "ChemistryProblemError",
        module.py().get_type::<ChemistryProblemError>(),
    )?;
    module.add("SanitizeError", module.py().get_type::<SanitizeError>())?;
    module.add("MatrixError", module.py().get_type::<MatrixError>())?;
    module.add("KekulizeError", module.py().get_type::<KekulizeError>())?;
    module.add_class::<SanitizeStage>()?;
    module.add_class::<SanitizeOperations>()?;
    module.add_class::<SanitizeParams>()?;
    module.add_class::<ChemistryProblem>()?;
    module.add_class::<ChemistryProblemReport>()?;
    module.add_class::<KekulizeParams>()?;
    module.add_class::<RingSearchParams>()?;
    module.add_class::<StructureTagParams>()?;
    module.add_class::<DistanceMatrixParams>()?;
    module.add_class::<DistanceMatrix3dParams>()?;
    module.add_class::<AtomPositionParams>()?;
    module.add_class::<RemoveHsParams>()?;
    module.add_class::<AromaticityModel>()?;
    module.add_class::<AromaticityParams>()?;
    module.add_class::<AddHsParams>()?;
    module.add_class::<DenseMatrix>()?;
    Ok(())
}

pyo3::create_exception!(cosmolkit, SanitizeError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(cosmolkit, MatrixError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(cosmolkit, KekulizeError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(
    cosmolkit,
    ChemistryProblemError,
    pyo3::exceptions::PyValueError
);

pub(crate) fn sanitize_pyerr(py: Python<'_>, source: ck::SanitizeError) -> PyErr {
    use ck::SanitizeError as E;
    let kind = match &source {
        E::InvalidOperations { .. } => "InvalidOperations",
        E::InvalidTopology { .. } => "InvalidTopology",
        E::InvalidQueryState { .. } => "InvalidQueryState",
        E::Cleanup { .. } => "Cleanup",
        E::Properties { .. } => "Properties",
        E::Rings { .. } => "Rings",
        E::Kekulize { .. } => "Kekulize",
        E::Radicals { .. } => "Radicals",
        E::Aromaticity { .. } => "Aromaticity",
        E::Conjugation { .. } => "Conjugation",
        E::Hybridization { .. } => "Hybridization",
        E::Atropisomers { .. } => "Atropisomers",
        E::Chirality { .. } => "Chirality",
        E::AdjustHs { .. } => "AdjustHs",
    };
    let error = crate::canonical_values::annotate(
        py,
        SanitizeError::new_err(source.to_string()),
        "sanitize",
        kind,
        &source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match &source {
            E::InvalidOperations { bits, unknown_bits } => {
                value.setattr("bits", *bits)?;
                value.setattr("unknown_bits", *unknown_bits)?;
            }
            E::InvalidTopology { stage, .. }
            | E::InvalidQueryState { stage, .. }
            | E::Cleanup { stage, .. }
            | E::Properties { stage, .. }
            | E::Rings { stage, .. }
            | E::Kekulize { stage, .. }
            | E::Radicals { stage, .. }
            | E::Aromaticity { stage, .. }
            | E::Conjugation { stage, .. }
            | E::Hybridization { stage, .. }
            | E::Atropisomers { stage, .. }
            | E::Chirality { stage, .. }
            | E::AdjustHs { stage, .. } => {
                value.setattr("stage", SanitizeStage::from(*stage))?;
            }
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(error) => error,
    }
}
pub(crate) fn matrix_pyerr(py: Python<'_>, source: ck::MatrixError) -> PyErr {
    use ck::MatrixError as E;
    let kind = match &source {
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::ActiveBondsWithoutAtoms => "ActiveBondsWithoutAtoms",
        E::ActiveAtomOutOfRange { .. } => "ActiveAtomOutOfRange",
        E::DuplicateActiveAtom { .. } => "DuplicateActiveAtom",
        E::ActiveBondOutOfRange { .. } => "ActiveBondOutOfRange",
        E::DuplicateActiveBond { .. } => "DuplicateActiveBond",
        E::ActiveBondEndpointMissing { .. } => "ActiveBondEndpointMissing",
        E::UnsupportedBondOrder { .. } => "UnsupportedBondOrder",
        E::No3dConformer => "No3dConformer",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::MatrixDimensionOverflow { .. } => "MatrixDimensionOverflow",
    };
    let error = crate::canonical_values::annotate(
        py,
        MatrixError::new_err(source.to_string()),
        "matrices",
        kind,
        &source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match &source {
            E::ActiveAtomOutOfRange {
                position,
                atom,
                atom_count,
            } => {
                value.setattr("position", *position)?;
                value.setattr("atom", atom.index())?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::DuplicateActiveAtom {
                atom,
                first_position,
                second_position,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("first_position", *first_position)?;
                value.setattr("second_position", *second_position)?;
            }
            E::ActiveBondOutOfRange {
                position,
                bond,
                bond_count,
            } => {
                value.setattr("position", *position)?;
                value.setattr("bond", bond.index())?;
                value.setattr("bond_count", *bond_count)?;
            }
            E::DuplicateActiveBond {
                bond,
                first_position,
                second_position,
            } => {
                value.setattr("bond", bond.index())?;
                value.setattr("first_position", *first_position)?;
                value.setattr("second_position", *second_position)?;
            }
            E::ActiveBondEndpointMissing {
                bond,
                endpoint,
                atom,
            } => {
                value.setattr("bond", bond.index())?;
                value.setattr("endpoint", *endpoint)?;
                value.setattr("atom", atom.index())?;
            }
            E::UnsupportedBondOrder { bond, order } => {
                value.setattr("bond", bond.index())?;
                value.setattr(
                    "order",
                    crate::canonical_atom_bond::enum_member(py, "BondOrder", *order as i64)?,
                )?;
            }
            E::ConformerNotFound { conformer_id } => {
                value.setattr("conformer_id", *conformer_id)?
            }
            E::MatrixDimensionOverflow { dimension } => value.setattr("dimension", *dimension)?,
            E::InvalidTopology(_)
            | E::InvalidCoordinates(_)
            | E::ActiveBondsWithoutAtoms
            | E::No3dConformer => {}
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(error) => error,
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum SanitizeStage {
    #[pyo3(name = "None_")]
    None = 0x000,
    Cleanup = 0x001,
    Properties = 0x002,
    SymmRings = 0x004,
    Kekulize = 0x008,
    FindRadicals = 0x010,
    SetAromaticity = 0x020,
    SetConjugation = 0x040,
    SetHybridization = 0x080,
    CleanupChirality = 0x100,
    AdjustHs = 0x200,
    CleanupOrganometallics = 0x400,
    CleanupAtropisomers = 0x800,
}

pub(crate) fn kekulize_pyerr(py: Python<'_>, source: ck::KekulizeError) -> PyErr {
    use ck::KekulizeError as E;
    let kind = match &source {
        E::InvalidTopology(_) => "InvalidTopology",
        E::AtomSelectionLength { .. } => "AtomSelectionLength",
        E::BondSelectionLength { .. } => "BondSelectionLength",
        E::CandidateAtomOutOfRange { .. } => "CandidateAtomOutOfRange",
        E::DuplicateCandidateAtom { .. } => "DuplicateCandidateAtom",
        E::MatchingStateLength { .. } => "MatchingStateLength",
        E::DuplicateDoneAtom { .. } => "DuplicateDoneAtom",
        E::MissingCandidateBond { .. } => "MissingCandidateBond",
        E::MissingBacktrackAnchor { .. } => "MissingBacktrackAnchor",
        E::QuestionSubsetOverflow { .. } => "QuestionSubsetOverflow",
        E::AromaticAtomOutsideRing { .. } => "AromaticAtomOutsideRing",
        E::NotKekulizable { .. } => "NotKekulizable",
        E::PostconditionValenceMismatch { .. } => "PostconditionValenceMismatch",
        E::UnsupportedQueryState { .. } => "UnsupportedQueryState",
        E::InvalidQueryState(_) => "InvalidQueryState",
        E::IntegerOverflow { .. } => "IntegerOverflow",
        E::RingFinding(_) => "RingFinding",
        E::Valence(_) => "Valence",
        E::CanonicalRank(_) => "CanonicalRank",
    };
    let error = crate::canonical_values::annotate(
        py,
        KekulizeError::new_err(source.to_string()),
        "kekulize",
        kind,
        &source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match &source {
            E::AtomSelectionLength { expected, actual }
            | E::BondSelectionLength { expected, actual } => {
                value.setattr("expected", *expected)?;
                value.setattr("actual", *actual)?;
            }
            E::CandidateAtomOutOfRange { atom, atom_count } => {
                value.setattr("atom", atom.index())?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::DuplicateCandidateAtom { atom }
            | E::DuplicateDoneAtom { atom }
            | E::MissingBacktrackAnchor { atom }
            | E::AromaticAtomOutsideRing { atom } => value.setattr("atom", atom.index())?,
            E::MatchingStateLength {
                field,
                expected,
                actual,
            } => {
                value.setattr("field", *field)?;
                value.setattr("expected", *expected)?;
                value.setattr("actual", *actual)?;
            }
            E::MissingCandidateBond { begin, end } => {
                value.setattr("begin", begin.index())?;
                value.setattr("end", end.index())?;
            }
            E::QuestionSubsetOverflow {
                questions,
                bit_width,
            } => {
                value.setattr("questions", *questions)?;
                value.setattr("bit_width", *bit_width)?;
            }
            E::NotKekulizable { problem_atoms } => value.setattr(
                "problem_atoms",
                problem_atoms
                    .iter()
                    .map(|id| id.index())
                    .collect::<Vec<_>>(),
            )?,
            E::PostconditionValenceMismatch {
                atom,
                before,
                after,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("before", *before)?;
                value.setattr("after", *after)?;
            }
            E::UnsupportedQueryState { bond, detail } => {
                value.setattr("bond", bond.index())?;
                value.setattr("detail", *detail)?;
            }
            E::IntegerOverflow { atom, field } => {
                value.setattr("atom", atom.index())?;
                value.setattr("field", *field)?;
            }
            E::Valence(source) => error.set_cause(
                py,
                Some(crate::canonical_atom_bond::valence_pyerr(
                    py,
                    source.clone(),
                )),
            ),
            E::InvalidTopology(_)
            | E::InvalidQueryState(_)
            | E::RingFinding(_)
            | E::CanonicalRank(_) => {}
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(error) => error,
    }
}

impl From<ck::SanitizeStage> for SanitizeStage {
    fn from(stage: ck::SanitizeStage) -> Self {
        match stage {
            ck::SanitizeStage::None => Self::None,
            ck::SanitizeStage::Cleanup => Self::Cleanup,
            ck::SanitizeStage::Properties => Self::Properties,
            ck::SanitizeStage::SymmRings => Self::SymmRings,
            ck::SanitizeStage::Kekulize => Self::Kekulize,
            ck::SanitizeStage::FindRadicals => Self::FindRadicals,
            ck::SanitizeStage::SetAromaticity => Self::SetAromaticity,
            ck::SanitizeStage::SetConjugation => Self::SetConjugation,
            ck::SanitizeStage::SetHybridization => Self::SetHybridization,
            ck::SanitizeStage::CleanupChirality => Self::CleanupChirality,
            ck::SanitizeStage::AdjustHs => Self::AdjustHs,
            ck::SanitizeStage::CleanupOrganometallics => Self::CleanupOrganometallics,
            ck::SanitizeStage::CleanupAtropisomers => Self::CleanupAtropisomers,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SanitizeOperations {
    inner: ck::SanitizeOperations,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SanitizeOperations {
    #[staticmethod]
    fn from_bits(py: Python<'_>, bits: u32) -> PyResult<Self> {
        ck::SanitizeOperations::from_bits(bits)
            .map(|inner| Self { inner })
            .map_err(|e| sanitize_pyerr(py, e))
    }
    fn bits(&self) -> u32 {
        self.inner.bits()
    }
    fn contains(&self, operation: &Self) -> bool {
        self.inner.contains(operation.inner)
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn __or__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }
    fn __and__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner & other.inner,
        }
    }
    #[classattr]
    #[pyo3(name = "NONE")]
    fn flag_none() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::NONE,
        }
    }
    #[classattr]
    #[pyo3(name = "CLEANUP")]
    fn flag_cleanup() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::CLEANUP,
        }
    }
    #[classattr]
    #[pyo3(name = "PROPERTIES")]
    fn flag_properties() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::PROPERTIES,
        }
    }
    #[classattr]
    #[pyo3(name = "SYMM_RINGS")]
    fn flag_symm_rings() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::SYMM_RINGS,
        }
    }
    #[classattr]
    #[pyo3(name = "KEKULIZE")]
    fn flag_kekulize() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::KEKULIZE,
        }
    }
    #[classattr]
    #[pyo3(name = "FIND_RADICALS")]
    fn flag_find_radicals() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::FIND_RADICALS,
        }
    }
    #[classattr]
    #[pyo3(name = "SET_AROMATICITY")]
    fn flag_set_aromaticity() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::SET_AROMATICITY,
        }
    }
    #[classattr]
    #[pyo3(name = "SET_CONJUGATION")]
    fn flag_set_conjugation() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::SET_CONJUGATION,
        }
    }
    #[classattr]
    #[pyo3(name = "SET_HYBRIDIZATION")]
    fn flag_set_hybridization() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::SET_HYBRIDIZATION,
        }
    }
    #[classattr]
    #[pyo3(name = "CLEANUP_CHIRALITY")]
    fn flag_cleanup_chirality() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::CLEANUP_CHIRALITY,
        }
    }
    #[classattr]
    #[pyo3(name = "ADJUST_HS")]
    fn flag_adjust_hs() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::ADJUST_HS,
        }
    }
    #[classattr]
    #[pyo3(name = "CLEANUP_ORGANOMETALLICS")]
    fn flag_cleanup_organometallics() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::CLEANUP_ORGANOMETALLICS,
        }
    }
    #[classattr]
    #[pyo3(name = "CLEANUP_ATROPISOMERS")]
    fn flag_cleanup_atropisomers() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::CLEANUP_ATROPISOMERS,
        }
    }
    #[classattr]
    #[pyo3(name = "ALL")]
    fn flag_all() -> SanitizeOperations {
        Self {
            inner: ck::SanitizeOperations::ALL,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SanitizeParams {
    pub(crate) inner: ck::SanitizeParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SanitizeParams {
    #[new]
    #[pyo3(signature = (operations=None))]
    fn new(operations: Option<&SanitizeOperations>) -> Self {
        Self {
            inner: ck::SanitizeParams {
                operations: operations.map_or(ck::SanitizeOperations::ALL, |value| value.inner),
            },
        }
    }
    #[getter]
    fn operations(&self) -> SanitizeOperations {
        SanitizeOperations {
            inner: self.inner.operations,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ChemistryProblem {
    inner: ck::ChemistryProblem,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl ChemistryProblem {
    #[getter]
    fn operation(&self) -> SanitizeStage {
        self.inner.operation.into()
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr = "ChemistryProblemError"))]
    fn error(&self, py: Python<'_>) -> Py<PyAny> {
        let kind = match &self.inner.error {
            ck::ChemistryProblemError::Valence(..) => "Valence",
            ck::ChemistryProblemError::Kekulize(..) => "Kekulize",
        };
        let error = crate::canonical_values::annotate(
            py,
            ChemistryProblemError::new_err(self.inner.error.to_string()),
            "sanitize",
            kind,
            &self.inner.error,
        );
        // Transparent Rust error wrappers may skip their inner error in
        // Error::source(). Preserve that typed payload at the language boundary.
        let cause = match &self.inner.error {
            ck::ChemistryProblemError::Valence(source) => {
                crate::canonical_atom_bond::valence_pyerr(py, source.clone())
            }
            ck::ChemistryProblemError::Kekulize(source) => kekulize_pyerr(py, source.clone()),
        };
        error.set_cause(py, Some(cause));
        error.value(py).clone().into_any().unbind()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ChemistryProblemReport {
    pub(crate) inner: ck::ChemistryProblemReport,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ChemistryProblemReport {
    #[getter]
    fn problems(&self) -> Vec<ChemistryProblem> {
        self.inner
            .problems
            .iter()
            .cloned()
            .map(|inner| ChemistryProblem { inner })
            .collect()
    }
}
