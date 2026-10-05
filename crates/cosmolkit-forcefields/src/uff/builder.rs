use std::collections::BTreeMap;
use std::sync::OnceLock;

use super::angle::AngleBendContrib;
use super::atom_typer::{
    NEEDS_EXPLICIT_HYDROGENS_WARNING_MESSAGE, UffAtomStateRef, UffTypingDiagnostic,
    UffTypingDiagnosticKind, UffTypingError, get_atom_types_from_state,
};
use super::bond::BondStretchContrib;
use super::inversion::{InversionContrib, InversionContributionError};
use super::nonbonded::{VdwContrib, calc_nonbonded_minimum};
use super::params::{AtomicParams, ParamCollection, UffParamError};
use super::torsion::TorsionAngleContrib;
use crate::geometry::{Point3, direction_vector};
use crate::kernel::{BondIndexArgument, ForceField, ForceFieldKernelError};
use cosmolkit_core::{
    FragmentCoordinateView, FragmentCoordinateViewError, MoleculeFragmentsError, RingInfo,
    ValenceAssignment, ValenceError, bond_type_as_double,
    get_molecule_fragments_with_coordinate_view,
};
use cosmolkit_model::{
    AtomId, BondId, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, Hybridization,
    MoleculeProperties, NeighborRef, TopologyBlock, TopologyValidationError,
};
use cosmolkit_search::{
    CompiledQuery, MatchError, QueryCompileError, SearchTarget, SmartsParseError,
    SmartsParseParams, SubstructMatchParams, build_query_match_context, parse_smarts,
    try_get_substruct_atom_matches_with_compiled_query_and_context,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum PreparedValenceField {
    Explicit,
    ImplicitHydrogen,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum UffBuilderError {
    NeighborMatrixIndexOutOfRange {
        n_atoms: usize,
        i: usize,
        j: usize,
    },
    NeighborMatrixIndexOverflow {
        n_atoms: usize,
    },
    NeighborMatrixAllocationFailed {
        byte_count: usize,
    },
    NeighborMatrixCellPositionOverflow {
        position: usize,
    },
    NeighborMatrixStorageOutOfRange {
        position: usize,
        byte_index: usize,
        storage_len: usize,
    },
    ValenceAssignmentLengthMismatch {
        field: PreparedValenceField,
        expected: usize,
        actual: usize,
    },
    SourceValencePrecondition {
        atom_id: AtomId,
        field: PreparedValenceField,
        value: i32,
    },
    SourceValenceOutOfRange {
        atom_id: AtomId,
        field: PreparedValenceField,
        value: i32,
    },
    ParamsLengthMismatch {
        atoms: usize,
        params: usize,
    },
    SelectedThreeDimensionalConformerNotFound {
        conformer_id: usize,
    },
    SelectedConformerCoordinateCountMismatch {
        conformer_id: usize,
        atoms: usize,
        coordinates: usize,
    },
    SourceBondIndexOverflow {
        atom_index: usize,
    },
    SourceTorsionIndexOverflow {
        atom_index: usize,
    },
    SourceAngleIndexOverflow {
        atom_index: usize,
    },
    SourceInversionIndexOverflow {
        atom_index: usize,
    },
    SourceDirectionVectorBelowTolerance {
        center_atom_index: usize,
        neighbor_atom_index: usize,
    },
    SourceTbpAtomPrecondition {
        center_atom_index: usize,
    },
    SourceTbpHybridizationPrecondition {
        center_atom_index: usize,
        actual: Hybridization,
    },
    SourceTbpDegreePrecondition {
        center_atom_index: usize,
        actual_degree: usize,
    },
    SourceTbpAxialBondNotFound {
        center_atom_index: usize,
    },
    SourceTbpEquatorialBondNotFound {
        center_atom_index: usize,
        role: u8,
    },
    SourceTbpCenterParamsMissing {
        center_atom_index: usize,
    },
    TorsionBondQuery(TorsionBondQueryError),
    InversionContribution(InversionContributionError),
    Valence(ValenceError),
    ForceFieldKernel(ForceFieldKernelError),
    TopologyValidation(TopologyValidationError),
}

impl std::fmt::Display for UffBuilderError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for UffBuilderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::NeighborMatrixIndexOutOfRange { .. } => None,
            Self::NeighborMatrixIndexOverflow { .. } => None,
            Self::NeighborMatrixAllocationFailed { .. } => None,
            Self::NeighborMatrixCellPositionOverflow { .. } => None,
            Self::NeighborMatrixStorageOutOfRange { .. } => None,
            Self::ValenceAssignmentLengthMismatch { .. } => None,
            Self::SourceValencePrecondition { .. } => None,
            Self::SourceValenceOutOfRange { .. } => None,
            Self::ParamsLengthMismatch { .. } => None,
            Self::SelectedThreeDimensionalConformerNotFound { .. } => None,
            Self::SelectedConformerCoordinateCountMismatch { .. } => None,
            Self::SourceBondIndexOverflow { .. } => None,
            Self::SourceTorsionIndexOverflow { .. } => None,
            Self::SourceAngleIndexOverflow { .. } => None,
            Self::SourceInversionIndexOverflow { .. } => None,
            Self::SourceDirectionVectorBelowTolerance { .. } => None,
            Self::SourceTbpAtomPrecondition { .. } => None,
            Self::SourceTbpHybridizationPrecondition { .. } => None,
            Self::SourceTbpDegreePrecondition { .. } => None,
            Self::SourceTbpAxialBondNotFound { .. } => None,
            Self::SourceTbpEquatorialBondNotFound { .. } => None,
            Self::SourceTbpCenterParamsMissing { .. } => None,
            Self::TorsionBondQuery(source) => Some(source),
            Self::InversionContribution(source) => Some(source),
            Self::Valence(source) => Some(source),
            Self::ForceFieldKernel(source) => Some(source),
            Self::TopologyValidation(source) => Some(source),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum DefaultTorsionQueryError {
    Parse(SmartsParseError),
    Compile(QueryCompileError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum TorsionBondQueryError {
    DefaultQuery(DefaultTorsionQueryError),
    Parse(SmartsParseError),
    Compile(QueryCompileError),
    Match(MatchError),
    MatchArity {
        smarts: String,
        actual: usize,
    },
    MissingMatchedBond {
        begin_atom_index: usize,
        end_atom_index: usize,
    },
}

impl std::fmt::Display for DefaultTorsionQueryError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for DefaultTorsionQueryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse(source) => Some(source),
            Self::Compile(source) => Some(source),
        }
    }
}

impl std::fmt::Display for TorsionBondQueryError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for TorsionBondQueryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::DefaultQuery(source) => Some(source),
            Self::Parse(source) => Some(source),
            Self::Compile(source) => Some(source),
            Self::Match(source) => Some(source),
            Self::MatchArity { .. } => None,
            Self::MissingMatchedBond { .. } => None,
        }
    }
}

pub(crate) const DEFAULT_TORSION_BOND_SMARTS: &str = "[!$(*#*)&!D1]~[!$(*#*)&!D1]";

/// Source `AtomicParamVect` rows borrow the canonical cached parameter table.
type UffParamsByAtom<'a> = [Option<&'a AtomicParams>];

#[derive(Debug)]
pub(crate) enum NonbondedFragmentMappingError {
    FragmentCopy(MoleculeFragmentsError),
    CoordinateView(FragmentCoordinateViewError),
    RequestedMappingMissing,
    AtomOutOfRange {
        component_index: usize,
        atom_index: usize,
        atom_count: usize,
    },
    DuplicateAtom {
        atom_index: usize,
        first_component: usize,
        second_component: usize,
    },
    MissingAtom {
        atom_index: usize,
    },
}

impl From<MoleculeFragmentsError> for NonbondedFragmentMappingError {
    fn from(error: MoleculeFragmentsError) -> Self {
        Self::FragmentCopy(error)
    }
}

#[derive(Debug)]
pub(super) enum NonbondedAssemblyError {
    FragmentMapping(NonbondedFragmentMappingError),
    PairBuilder(UffBuilderError),
}

impl std::fmt::Display for NonbondedFragmentMappingError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for NonbondedFragmentMappingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::FragmentCopy(source) => Some(source),
            Self::CoordinateView(source) => Some(source),
            Self::RequestedMappingMissing => None,
            Self::AtomOutOfRange { .. } => None,
            Self::DuplicateAtom { .. } => None,
            Self::MissingAtom { .. } => None,
        }
    }
}

impl std::fmt::Display for NonbondedAssemblyError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for NonbondedAssemblyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::FragmentMapping(source) => Some(source),
            Self::PairBuilder(source) => Some(source),
        }
    }
}

#[derive(Debug)]
pub(super) enum ForceFieldConstructionError {
    Builder(UffBuilderError),
    Nonbonded(NonbondedAssemblyError),
}

#[derive(Debug)]
pub(super) enum AutomaticForceFieldConstructionError {
    ParameterTable(UffParamError),
    Typing(UffTypingError),
    Construction(ForceFieldConstructionError),
}

impl std::fmt::Display for ForceFieldConstructionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for ForceFieldConstructionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Builder(source) => Some(source),
            Self::Nonbonded(source) => Some(source),
        }
    }
}

impl std::fmt::Display for AutomaticForceFieldConstructionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for AutomaticForceFieldConstructionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ParameterTable(source) => Some(source),
            Self::Typing(source) => Some(source),
            Self::Construction(source) => Some(source),
        }
    }
}

impl From<UffBuilderError> for ForceFieldConstructionError {
    fn from(error: UffBuilderError) -> Self {
        Self::Builder(error)
    }
}

impl From<NonbondedAssemblyError> for ForceFieldConstructionError {
    fn from(error: NonbondedAssemblyError) -> Self {
        Self::Nonbonded(error)
    }
}

pub(crate) fn prepare_nonbonded_fragment_mapping(
    topology: &TopologyBlock,
    coordinates: &FragmentCoordinateView<'_>,
    molecule_properties: &MoleculeProperties,
    ignore_interfragment_interactions: bool,
) -> Result<Option<Vec<usize>>, NonbondedFragmentMappingError> {
    // BEGIN RDKIT CPP FUNCTION UFF::addNonbonded fragment preparation
    // RDKit❗❌:   INT_VECT fragMapping;
    // RDKit❗❌:   if (ignoreInterfragInteractions) {
    // RDKit❗❌:     std::vector<ROMOL_SPTR> molFrags =
    // RDKit❗❌:         MolOps::getMolFrags(mol, true, &fragMapping);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION UFF::addNonbonded fragment preparation
    // Behavior: false returns before component construction. True completes
    // the source-ordered sanitized F18 copy and final-sanitize pipeline before
    // translating its component rows to one label per source atom. Fragment
    // values are dropped only after extraction and map construction succeed.
    // The map checks return typed internal consistency errors instead of
    // leaving a default or sentinel label in the result.
    // Complexity: core extraction retains its O(FV) component projection,
    // detached-copy and sanitation work. Translation adds O(V+F) time and O(V)
    // label storage; the coordinate view is borrowed and avoids a second input
    // coordinate-block clone.
    if !ignore_interfragment_interactions {
        return Ok(None);
    }

    let fragments = get_molecule_fragments_with_coordinate_view(
        topology,
        coordinates,
        molecule_properties,
        true,
        true,
    )?;
    let mut labels = vec![None; topology.atoms.len()];
    for (component_index, fragment) in fragments.iter().enumerate() {
        for atom in fragment.component_atoms() {
            let atom_index = atom.index();
            let Some(slot) = labels.get_mut(atom_index) else {
                return Err(NonbondedFragmentMappingError::AtomOutOfRange {
                    component_index,
                    atom_index,
                    atom_count: topology.atoms.len(),
                });
            };
            if let Some(first_component) = *slot {
                return Err(NonbondedFragmentMappingError::DuplicateAtom {
                    atom_index,
                    first_component,
                    second_component: component_index,
                });
            }
            *slot = Some(component_index);
        }
    }
    let labels = labels
        .into_iter()
        .enumerate()
        .map(|(atom_index, component)| {
            component.ok_or(NonbondedFragmentMappingError::MissingAtom { atom_index })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(labels))
}

fn nonbonded_pair_is_eligible(
    n_atoms: usize,
    atom_i: usize,
    atom_j: usize,
    params: &UffParamsByAtom<'_>,
    fragment_mapping: &[usize],
    ignore_interfragment_interactions: bool,
    neighbor_matrix: &[u8],
) -> Result<bool, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::addNonbonded outer parameter guard
    // (Builder.cpp:447-450; the outer index is supplied by the caller)
    // RDKit❗✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit❗✔️:     if (!params[i]) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // END RDKIT CPP FUNCTION UFF::addNonbonded outer parameter guard
    if params[atom_i].is_none() {
        return Ok(false);
    }
    nonbonded_pair_is_eligible_after_i_parameter(
        n_atoms,
        atom_i,
        atom_j,
        params[atom_j],
        fragment_mapping,
        ignore_interfragment_interactions,
        neighbor_matrix,
    )
}

fn nonbonded_pair_is_eligible_after_i_parameter(
    n_atoms: usize,
    atom_i: usize,
    atom_j: usize,
    params_j: Option<&AtomicParams>,
    fragment_mapping: &[usize],
    ignore_interfragment_interactions: bool,
    neighbor_matrix: &[u8],
) -> Result<bool, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::addNonbonded pair eligibility
    // (Builder.cpp:447-459; loop variables are supplied by the caller)
    // RDKit❗✔️:     for (unsigned int j = i + 1; j < nAtoms; j++) {
    // RDKit❗✔️:       if (!params[j] ||
    // RDKit❗✔️:           (ignoreInterfragInteractions && fragMapping[i] != fragMapping[j])) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j)) >=
    // RDKit❗✔️:           RELATION_1_4) {
    // END RDKIT CPP FUNCTION UFF::addNonbonded pair eligibility
    // Source-call invariants: addNonbonded checks atom/parameter length before
    // entering these loops, i and j are in-range with i < j, and F19 supplies
    // a complete fragment mapping only when the flag is true. Preserve each
    // early continue before reading later pair state.
    if params_j.is_none() {
        return Ok(false);
    }
    if ignore_interfragment_interactions && fragment_mapping[atom_i] != fragment_mapping[atom_j] {
        return Ok(false);
    }

    const RELATION_1_4: u8 = 2;
    let relation_position = two_bit_cell_pos(n_atoms, atom_i, atom_j)?;
    let relation = get_two_bit_cell(neighbor_matrix, relation_position)?;
    Ok(relation >= RELATION_1_4)
}

fn append_nonbonded_pair_if_within_threshold(
    field: &mut ForceField<'_>,
    atom_i: u32,
    atom_j: u32,
    params_i: &AtomicParams,
    params_j: &AtomicParams,
    vdw_threshold: f64,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::addNonbonded distance and append
    // (Builder.cpp:459-465)
    // RDKit❗✔️:         double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
    // RDKit❗✔️:         if (dist < vdwThresh *
    // RDKit❗✔️:                        UFF::Utils::calcNonbondedMinimum(params[i], params[j])) {
    // RDKit❗✔️:           vdWContrib *contrib;
    // RDKit❗✔️:           contrib = new vdWContrib(field, i, j, params[i], params[j]);
    // RDKit❗✔️:           field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:         }
    // END RDKIT CPP FUNCTION UFF::addNonbonded distance and append
    // BEGIN RDKIT CPP HELPER UFF::Utils::calcNonbondedMinimum
    // (ForceField/UFF/Nonbonded.cpp:21-24)
    // RDKit❗✔️: double calcNonbondedMinimum(const AtomicParams *at1Params,
    // RDKit❗✔️:                             const AtomicParams *at2Params) {
    // RDKit❗✔️:   return sqrt(at1Params->x1 * at2Params->x1);
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER UFF::Utils::calcNonbondedMinimum

    let (point_i, point_j) = {
        let positions = field.positions();
        let position_i = positions.get(atom_i as usize).ok_or_else(|| {
            UffBuilderError::ForceFieldKernel(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: atom_i,
                upper_bound: positions.len(),
            })
        })?;
        let position_j = positions.get(atom_j as usize).ok_or_else(|| {
            UffBuilderError::ForceFieldKernel(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: atom_j,
                upper_bound: positions.len(),
            })
        })?;
        (
            Point3 {
                x: position_i[0],
                y: position_i[1],
                z: position_i[2],
            },
            Point3 {
                x: position_j[0],
                y: position_j[1],
                z: position_j[2],
            },
        )
    };
    let distance = Point3::difference(&point_i, &point_j).length();
    if distance < vdw_threshold * calc_nonbonded_minimum(params_i, params_j) {
        let contribution = VdwContrib::new(field.positions(), atom_i, atom_j, params_i, params_j)
            .map_err(UffBuilderError::ForceFieldKernel)?;
        field.add_contribution(Box::new(contribution));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn add_nonbonded(
    topology: &TopologyBlock,
    conformers_2d: &[Conformer2D],
    conformers_3d_before: &[Conformer3D],
    selected_conformer_id: usize,
    selected_conformer_is_3d: bool,
    selected_conformer_props: &BTreeMap<String, String>,
    conformers_3d_after: &[Conformer3D],
    source_coordinate_dim: Option<CoordinateDimension>,
    molecule_properties: &MoleculeProperties,
    params: &UffParamsByAtom<'_>,
    field: &mut ForceField<'_>,
    neighbor_matrix: &[u8],
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
) -> Result<(), NonbondedAssemblyError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addNonbonded (Builder.cpp:433-465)
    // RDKit❗✔️: void addNonbonded(const ROMol &mol, int confId, const AtomicParamVect &params,
    // RDKit❗✔️:                   ForceFields::ForceField *field,
    // RDKit❗✔️:                   boost::shared_array<std::uint8_t> neighborMatrix,
    // RDKit❗✔️:                   double vdwThresh, bool ignoreInterfragInteractions) {
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗✔️:   PRECONDITION(field, "bad forcefield");
    // RDKit❗✔️:
    // RDKit❗✔️:   INT_VECT fragMapping;
    // RDKit❗✔️:   if (ignoreInterfragInteractions) {
    // RDKit❗✔️:     std::vector<ROMOL_SPTR> molFrags =
    // RDKit❗✔️:         MolOps::getMolFrags(mol, true, &fragMapping);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit❗✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit❗✔️:     if (!params[i]) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (unsigned int j = i + 1; j < nAtoms; j++) {
    // RDKit❗✔️:       if (!params[j] ||
    // RDKit❗✔️:           (ignoreInterfragInteractions && fragMapping[i] != fragMapping[j])) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j)) >=
    // RDKit❗✔️:           RELATION_1_4) {
    // RDKit❗✔️:         double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
    // RDKit❗✔️:         if (dist < vdwThresh *
    // RDKit❗✔️:                        UFF::Utils::calcNonbondedMinimum(params[i], params[j])) {
    // RDKit❗✔️:           vdWContrib *contrib;
    // RDKit❗✔️:           contrib = new vdWContrib(field, i, j, params[i], params[j]);
    // RDKit❗✔️:           field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addNonbonded
    // Rust borrows encode the non-null ForceField precondition. Source
    // getConformer(confId) has already selected these mutable kernel rows; the
    // split conformer metadata below lets F19 borrow them for fragment copying.
    // The narrow borrow block ends before any contribution is appended.
    if topology.atoms.len() != params.len() {
        return Err(NonbondedAssemblyError::PairBuilder(
            UffBuilderError::ParamsLengthMismatch {
                atoms: topology.atoms.len(),
                params: params.len(),
            },
        ));
    }

    let prepared_fragment_mapping = if ignore_interfragment_interactions {
        let selected_kernel_rows = field
            .positions()
            .iter()
            .map(|coordinates| &coordinates[..])
            .collect::<Vec<_>>();
        let coordinate_view = FragmentCoordinateView::from_split_conformers(
            conformers_2d,
            conformers_3d_before,
            selected_conformer_id,
            selected_conformer_is_3d,
            selected_conformer_props,
            &selected_kernel_rows,
            conformers_3d_after,
            source_coordinate_dim,
        )
        .map_err(|error| {
            NonbondedAssemblyError::FragmentMapping(NonbondedFragmentMappingError::CoordinateView(
                error,
            ))
        })?;
        prepare_nonbonded_fragment_mapping(topology, &coordinate_view, molecule_properties, true)
            .map_err(NonbondedAssemblyError::FragmentMapping)?
    } else {
        None
    };
    let fragment_mapping = match prepared_fragment_mapping.as_deref() {
        Some(mapping) => mapping,
        None if !ignore_interfragment_interactions => &[],
        None => {
            return Err(NonbondedAssemblyError::FragmentMapping(
                NonbondedFragmentMappingError::RequestedMappingMissing,
            ));
        }
    };

    let n_atoms = topology.atoms.len();
    for atom_i in 0..n_atoms {
        let Some(params_i) = params[atom_i] else {
            continue;
        };
        for atom_j in atom_i + 1..n_atoms {
            let params_j = params[atom_j];
            let eligible = nonbonded_pair_is_eligible_after_i_parameter(
                n_atoms,
                atom_i,
                atom_j,
                params_j,
                fragment_mapping,
                ignore_interfragment_interactions,
                neighbor_matrix,
            )
            .map_err(NonbondedAssemblyError::PairBuilder)?;
            if !eligible {
                continue;
            }
            let Some(params_j) = params_j else {
                continue;
            };
            let atom_i = u32::try_from(atom_i).map_err(|_| {
                NonbondedAssemblyError::PairBuilder(UffBuilderError::NeighborMatrixIndexOverflow {
                    n_atoms,
                })
            })?;
            let atom_j = u32::try_from(atom_j).map_err(|_| {
                NonbondedAssemblyError::PairBuilder(UffBuilderError::NeighborMatrixIndexOverflow {
                    n_atoms,
                })
            })?;
            append_nonbonded_pair_if_within_threshold(
                field,
                atom_i,
                atom_j,
                params_i,
                params_j,
                vdw_threshold,
            )
            .map_err(NonbondedAssemblyError::PairBuilder)?;
        }
    }

    // Behavior marker — RDKit❗✔️: complete the optional sanitized-fragment
    // stage before pairs, preserve nullable parameter/label/matrix/distance
    // short circuits, and append each qualifying pair in ascending order.
    // Complexity marker — RDKit❗✔️: one O(FV) source fragment pipeline and
    // O(A) labels only when requested, followed by O(A^2) candidate visits;
    // each accepted candidate performs fixed scalar work and at most one box.
    Ok(())
}

fn default_torsion_query() -> Result<&'static CompiledQuery, DefaultTorsionQueryError> {
    // BEGIN RDKIT CPP HELPER Builder.cpp::DefaultTorsionBondSmarts::create/query
    // RDKit❗✔️: const std::string DefaultTorsionBondSmarts::ds_string =
    // RDKit❗✔️:     "[!$(*#*)&!D1]~[!$(*#*)&!D1]";
    // RDKit❗✔️: void DefaultTorsionBondSmarts::create() {
    // RDKit❗✔️:   ds_instance.reset(SmartsToMol(ds_string));
    // RDKit❗✔️: }
    // RDKit❗✔️: const ROMol *DefaultTorsionBondSmarts::query() {
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:   std::call_once(ds_flag, create);
    // RDKit❗✔️: #else
    // RDKit❗✔️:   static bool created = false;
    // RDKit❗✔️:   if (!created) {
    // RDKit❗✔️:     created = true;
    // RDKit❗✔️:     create();
    // RDKit❗✔️:   }
    // RDKit❗✔️: #endif
    // RDKit❗✔️:   return ds_instance.get();
    // RDKit❗✔️: }
    // RDKit❗✔️: static const std::string &string() { return ds_string; }
    // END RDKIT CPP HELPER Builder.cpp::DefaultTorsionBondSmarts::create/query

    // Keep the exact source SMARTS as one immutable search-owned execution plan.
    // OnceLock retains a typed construction failure rather than replacing it
    // with an empty query or a hand-written bond predicate.
    static QUERY: OnceLock<Result<CompiledQuery, DefaultTorsionQueryError>> = OnceLock::new();
    QUERY
        .get_or_init(|| {
            let query = parse_smarts(DEFAULT_TORSION_BOND_SMARTS, &SmartsParseParams::default())
                .map_err(DefaultTorsionQueryError::Parse)?;
            CompiledQuery::compile(query).map_err(DefaultTorsionQueryError::Compile)
        })
        .as_ref()
        .map_err(Clone::clone)
}

pub(crate) fn torsion_bond_matches(
    topology: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    torsion_bond_smarts: &str,
) -> Result<Vec<Vec<usize>>, TorsionBondQueryError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTorsions query selection (Builder.cpp:505-515)
    // RDKit❗❌:   // find all of the torsion bonds:
    // RDKit❗❌:   std::vector<MatchVectType> matchVect;
    // RDKit❗❌:   const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
    // RDKit❗❌:   const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
    // RDKit❗❌:                            ? defaultQuery
    // RDKit❗❌:                            : SmartsToMol(torsionBondSmarts);
    // RDKit❗❌:   TEST_ASSERT(query);
    // RDKit❗❌:   unsigned int nHits = SubstructMatch(mol, *query, matchVect);
    // RDKit❗❌:   if (query != defaultQuery) {
    // RDKit❗❌:     delete query;
    // RDKit❗❌:   }
    // RDKit❗❌:   for (unsigned int i = 0; i < nHits; i++) {
    // RDKit❗❌:     MatchVectType match = matchVect[i];
    // RDKit❗❌:     TEST_ASSERT(match.size() == 2);
    // RDKit❗❌:     int idx1 = match[0].second;
    // RDKit❗❌:     int idx2 = match[1].second;
    // END RDKIT CPP FUNCTION UFF::Tools::addTorsions query selection

    let default_query = default_torsion_query().map_err(TorsionBondQueryError::DefaultQuery)?;
    let custom_query = if torsion_bond_smarts == DEFAULT_TORSION_BOND_SMARTS {
        None
    } else {
        let query = parse_smarts(torsion_bond_smarts, &SmartsParseParams::default())
            .map_err(TorsionBondQueryError::Parse)?;
        Some(CompiledQuery::compile(query).map_err(TorsionBondQueryError::Compile)?)
    };

    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        topology,
        &coordinates,
        &topology.stereo_groups,
        Some(rings),
        Some(valence),
    );
    let query = custom_query.as_ref().unwrap_or(default_query);
    let match_params = SubstructMatchParams::default();
    let query_context = build_query_match_context(&target);
    let matches = try_get_substruct_atom_matches_with_compiled_query_and_context(
        &target,
        query,
        &match_params,
        &query_context,
    )
    .map_err(MatchError::from)
    .map_err(TorsionBondQueryError::Match)?;
    // The custom compiled query is dropped after matching and before ordered
    // match validation, mirroring RDKit's delete immediately after SubstructMatch.
    drop(custom_query);

    for matched in &matches {
        if matched.len() != 2 {
            return Err(TorsionBondQueryError::MatchArity {
                smarts: torsion_bond_smarts.to_owned(),
                actual: matched.len(),
            });
        }
    }

    // Reuse the canonical atom rows in source order. The atom-only projection
    // avoids materializing Search's separate per-hit bond mapping.
    Ok(matches)
}

/// Resolve the source `getBondBetweenAtoms` invariant after the source
/// per-center parameter guard has accepted both matched atoms.
pub(crate) fn source_torsion_bond_index(
    topology: &TopologyBlock,
    begin_atom_index: usize,
    end_atom_index: usize,
) -> Result<usize, TorsionBondQueryError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTorsions target bond lookup (Builder.cpp:525-527)
    // RDKit❗❌:     const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
    // RDKit❗❌:     TEST_ASSERT(bond);
    // END RDKIT CPP FUNCTION UFF::Tools::addTorsions target bond lookup

    topology
        .adjacency
        .neighbors_of(begin_atom_index)
        .iter()
        .find(|neighbor| neighbor.atom_index == end_atom_index)
        .map(|neighbor| neighbor.bond.index())
        .ok_or(TorsionBondQueryError::MissingMatchedBond {
            begin_atom_index,
            end_atom_index,
        })
}

fn torsions_for_bond(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    begin_atom_index: usize,
    end_atom_index: usize,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTorsions per-match expansion (Builder.cpp:522-584)
    // RDKit❗✔️:     if (!params[idx1] || !params[idx2]) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
    // RDKit❗✔️:     std::vector<TorsionAngleContrib *> contribsHere;
    // RDKit❗✔️:     TEST_ASSERT(bond);
    // RDKit❗✔️:     const Atom *atom1 = mol.getAtomWithIdx(idx1);
    // RDKit❗✔️:     const Atom *atom2 = mol.getAtomWithIdx(idx2);
    // RDKit❗✔️:     if ((atom1->getHybridization() == Atom::SP2 ||
    // RDKit❗✔️:          atom1->getHybridization() == Atom::SP3) &&
    // RDKit❗✔️:         (atom2->getHybridization() == Atom::SP2 ||
    // RDKit❗✔️:          atom2->getHybridization() == Atom::SP3)) {
    // RDKit❗✔️:       ROMol::OEDGE_ITER beg1, end1;
    // RDKit❗✔️:       boost::tie(beg1, end1) = mol.getAtomBonds(atom1);
    // RDKit❗✔️:       while (beg1 != end1) {
    // RDKit❗✔️:         const Bond *tBond1 = mol[*beg1];
    // RDKit❗✔️:         if (tBond1 != bond) {
    // RDKit❗✔️:           int bIdx = tBond1->getOtherAtomIdx(idx1);
    // RDKit❗✔️:           ROMol::OEDGE_ITER beg2, end2;
    // RDKit❗✔️:           boost::tie(beg2, end2) = mol.getAtomBonds(atom2);
    // RDKit❗✔️:           while (beg2 != end2) {
    // RDKit❗✔️:             const Bond *tBond2 = mol[*beg2];
    // RDKit❗✔️:             if (tBond2 != bond && tBond2 != tBond1) {
    // RDKit❗✔️:               int eIdx = tBond2->getOtherAtomIdx(idx2);
    // RDKit❗✔️:               // make sure this isn't a three-membered ring:
    // RDKit❗✔️:               if (eIdx != bIdx) {
    // RDKit❗✔️:                 // we now have a torsion involving atoms (bonds):
    // RDKit❗✔️:                 //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
    // RDKit❗✔️:                 TorsionAngleContrib *contrib;
    // RDKit❗✔️:                 // if either of the end atoms is SP2 hybridized, set a flag
    // RDKit❗✔️:                 // here.
    // RDKit❗✔️:                 bool hasSP2 = false;
    // RDKit❗✔️:                 if (mol.getAtomWithIdx(bIdx)->getHybridization() == Atom::SP2 ||
    // RDKit❗✔️:                     mol.getAtomWithIdx(eIdx)->getHybridization() == Atom::SP2) {
    // RDKit❗✔️:                   hasSP2 = true;
    // RDKit❗✔️:                 }
    // RDKit❗✔️:                 // std::cout << "Torsion: " << bIdx << "-" << idx1 << "-" <<
    // RDKit❗✔️:                 // idx2 << "-" << eIdx << std::endl;
    // RDKit❗✔️:                 // if(okToIncludeTorsion(mol,bond,bIdx,idx1,idx2,eIdx)){
    // RDKit❗✔️:                 // std::cout << "  INCLUDED" << std::endl;
    // RDKit❗✔️:                 contrib = new TorsionAngleContrib(
    // RDKit❗✔️:                     field, bIdx, idx1, idx2, eIdx, bond->getBondTypeAsDouble(),
    // RDKit❗✔️:                     atom1->getAtomicNum(), atom2->getAtomicNum(),
    // RDKit❗✔️:                     atom1->getHybridization(), atom2->getHybridization(),
    // RDKit❗✔️:                     params[idx1], params[idx2], hasSP2);
    // RDKit❗✔️:                 field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:                 contribsHere.push_back(contrib);
    // RDKit❗✔️:                 //}
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:             beg2++;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         beg1++;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // now divide the force constant for each contribution to the torsion energy
    // RDKit❗✔️:     // about this bond by the number of contribs about this bond:
    // RDKit❗✔️:     for (auto chI = contribsHere.begin(); chI != contribsHere.end(); ++chI) {
    // RDKit❗✔️:       (*chI)->scaleForceConstant(contribsHere.size());
    // RDKit❗✔️:     }
    // END RDKIT CPP FUNCTION UFF::Tools::addTorsions per-match expansion

    let Some(begin_params) = params[begin_atom_index] else {
        return Ok(());
    };
    let Some(end_params) = params[end_atom_index] else {
        return Ok(());
    };

    let central_bond_index = source_torsion_bond_index(topology, begin_atom_index, end_atom_index)
        .map_err(UffBuilderError::TorsionBondQuery)?;
    let central_bond = &topology.bonds[central_bond_index];
    let atom1 = &topology.atoms[begin_atom_index];
    let atom2 = &topology.atoms[end_atom_index];
    let mut contributions_here = Vec::new();
    let mut construction_error = None;

    if (atom1.hybridization() == Hybridization::Sp2 || atom1.hybridization() == Hybridization::Sp3)
        && (atom2.hybridization() == Hybridization::Sp2
            || atom2.hybridization() == Hybridization::Sp3)
    {
        'begin_neighbors: for first_neighbor in topology.adjacency.neighbors_of(begin_atom_index) {
            let first_bond_index = first_neighbor.bond.index();
            if first_bond_index != central_bond_index {
                let first_terminal_index = first_neighbor.atom_index;
                for second_neighbor in topology.adjacency.neighbors_of(end_atom_index) {
                    let second_bond_index = second_neighbor.bond.index();
                    if second_bond_index != central_bond_index
                        && second_bond_index != first_bond_index
                    {
                        let second_terminal_index = second_neighbor.atom_index;
                        if second_terminal_index != first_terminal_index {
                            let has_sp2 = topology.atoms[first_terminal_index].hybridization()
                                == Hybridization::Sp2
                                || topology.atoms[second_terminal_index].hybridization()
                                    == Hybridization::Sp2;

                            let contribution = (|| {
                                let first_terminal_index = u32::try_from(first_terminal_index)
                                    .map_err(|_| UffBuilderError::SourceTorsionIndexOverflow {
                                        atom_index: first_terminal_index,
                                    })?;
                                let begin_source_index =
                                    u32::try_from(begin_atom_index).map_err(|_| {
                                        UffBuilderError::SourceTorsionIndexOverflow {
                                            atom_index: begin_atom_index,
                                        }
                                    })?;
                                let end_source_index =
                                    u32::try_from(end_atom_index).map_err(|_| {
                                        UffBuilderError::SourceTorsionIndexOverflow {
                                            atom_index: end_atom_index,
                                        }
                                    })?;
                                let second_terminal_index = u32::try_from(second_terminal_index)
                                    .map_err(|_| UffBuilderError::SourceTorsionIndexOverflow {
                                        atom_index: second_terminal_index,
                                    })?;
                                let bond_order = bond_type_as_double(central_bond.order())
                                    .map_err(UffBuilderError::Valence)?;
                                let atomic_number1 = i32::from(atom1.atomic_number());
                                let atomic_number2 = i32::from(atom2.atomic_number());
                                let hybridization1 = atom1.hybridization();
                                let hybridization2 = atom2.hybridization();

                                TorsionAngleContrib::new(
                                    field.positions(),
                                    first_terminal_index,
                                    begin_source_index,
                                    end_source_index,
                                    second_terminal_index,
                                    bond_order,
                                    atomic_number1,
                                    atomic_number2,
                                    hybridization1,
                                    hybridization2,
                                    begin_params,
                                    end_params,
                                    has_sp2,
                                )
                                .map_err(|error| {
                                    UffBuilderError::ForceFieldKernel(ForceFieldKernelError::from(
                                        error,
                                    ))
                                })
                            })();

                            match contribution {
                                Ok(contribution) => contributions_here.push(contribution),
                                Err(error) => {
                                    construction_error = Some(error);
                                    break 'begin_neighbors;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if let Some(error) = construction_error {
        for contribution in contributions_here {
            field.add_contribution(Box::new(contribution));
        }
        return Err(error);
    }

    let contribution_count = contributions_here.len() as u32;
    for contribution in &mut contributions_here {
        contribution.scale_force_constant(contribution_count);
    }
    for contribution in contributions_here {
        field.add_contribution(Box::new(contribution));
    }

    // Behavior marker — RDKit❗✔️: preserve query-side endpoint orientation,
    // nested source adjacency order, center-only parameter skips, exclusions,
    // partial-prefix errors, and this match's independent final scaling.
    // Complexity marker — RDKit✔️✔️: O(degree(begin) * degree(end)) traversal,
    // one local per-match term vector, and one boxed kernel term per survivor.
    Ok(())
}

pub(super) fn add_torsions(
    topology: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    params: &UffParamsByAtom<'_>,
    torsion_bond_smarts: &str,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTorsions (Builder.cpp:499-585)
    // RDKit❗❌: void addTorsions(const ROMol &mol, const AtomicParamVect &params,
    // RDKit❗❌:                  ForceFields::ForceField *field,
    // RDKit❗❌:                  const std::string &torsionBondSmarts) {
    // RDKit❗❌:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗❌:   PRECONDITION(field, "bad forcefield");
    // RDKit❗❌:
    // RDKit❗❌:   // find all of the torsion bonds:
    // RDKit❗❌:   std::vector<MatchVectType> matchVect;
    // RDKit❗❌:   const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
    // RDKit❗❌:   const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
    // RDKit❗❌:                            ? defaultQuery
    // RDKit❗❌:                            : SmartsToMol(torsionBondSmarts);
    // RDKit❗❌:   TEST_ASSERT(query);
    // RDKit❗❌:   unsigned int nHits = SubstructMatch(mol, *query, matchVect);
    // RDKit❗❌:   if (query != defaultQuery) {
    // RDKit❗❌:     delete query;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int i = 0; i < nHits; i++) {
    // RDKit❗❌:     MatchVectType match = matchVect[i];
    // RDKit❗❌:     TEST_ASSERT(match.size() == 2);
    // RDKit❗❌:     int idx1 = match[0].second;
    // RDKit❗❌:     int idx2 = match[1].second;
    // RDKit❗❌:     if (!params[idx1] || !params[idx2]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
    // RDKit❗❌:     std::vector<TorsionAngleContrib *> contribsHere;
    // RDKit❗❌:     TEST_ASSERT(bond);
    // RDKit❗❌:     const Atom *atom1 = mol.getAtomWithIdx(idx1);
    // RDKit❗❌:     const Atom *atom2 = mol.getAtomWithIdx(idx2);
    // RDKit❗❌:
    // RDKit❗❌:     if ((atom1->getHybridization() == Atom::SP2 ||
    // RDKit❗❌:          atom1->getHybridization() == Atom::SP3) &&
    // RDKit❗❌:         (atom2->getHybridization() == Atom::SP2 ||
    // RDKit❗❌:          atom2->getHybridization() == Atom::SP3)) {
    // RDKit❗❌:       ROMol::OEDGE_ITER beg1, end1;
    // RDKit❗❌:       boost::tie(beg1, end1) = mol.getAtomBonds(atom1);
    // RDKit❗❌:       while (beg1 != end1) {
    // RDKit❗❌:         const Bond *tBond1 = mol[*beg1];
    // RDKit❗❌:         if (tBond1 != bond) {
    // RDKit❗❌:           int bIdx = tBond1->getOtherAtomIdx(idx1);
    // RDKit❗❌:           ROMol::OEDGE_ITER beg2, end2;
    // RDKit❗❌:           boost::tie(beg2, end2) = mol.getAtomBonds(atom2);
    // RDKit❗❌:           while (beg2 != end2) {
    // RDKit❗❌:             const Bond *tBond2 = mol[*beg2];
    // RDKit❗❌:             if (tBond2 != bond && tBond2 != tBond1) {
    // RDKit❗❌:               int eIdx = tBond2->getOtherAtomIdx(idx2);
    // RDKit❗❌:               // make sure this isn't a three-membered ring:
    // RDKit❗❌:               if (eIdx != bIdx) {
    // RDKit❗❌:                 // we now have a torsion involving atoms (bonds):
    // RDKit❗❌:                 //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
    // RDKit❗❌:                 TorsionAngleContrib *contrib;
    // RDKit❗❌:
    // RDKit❗❌:                 // if either of the end atoms is SP2 hybridized, set a flag
    // RDKit❗❌:                 // here.
    // RDKit❗❌:                 bool hasSP2 = false;
    // RDKit❗❌:                 if (mol.getAtomWithIdx(bIdx)->getHybridization() == Atom::SP2 ||
    // RDKit❗❌:                     mol.getAtomWithIdx(eIdx)->getHybridization() == Atom::SP2) {
    // RDKit❗❌:                   hasSP2 = true;
    // RDKit❗❌:                 }
    // RDKit❗❌:                 // std::cout << "Torsion: " << bIdx << "-" << idx1 << "-" <<
    // RDKit❗❌:                 // idx2 << "-" << eIdx << std::endl;
    // RDKit❗❌:                 // if(okToIncludeTorsion(mol,bond,bIdx,idx1,idx2,eIdx)){
    // RDKit❗❌:                 // std::cout << "  INCLUDED" << std::endl;
    // RDKit❗❌:                 contrib = new TorsionAngleContrib(
    // RDKit❗❌:                     field, bIdx, idx1, idx2, eIdx, bond->getBondTypeAsDouble(),
    // RDKit❗❌:                     atom1->getAtomicNum(), atom2->getAtomicNum(),
    // RDKit❗❌:                     atom1->getHybridization(), atom2->getHybridization(),
    // RDKit❗❌:                     params[idx1], params[idx2], hasSP2);
    // RDKit❗❌:                 field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗❌:                 contribsHere.push_back(contrib);
    // RDKit❗❌:                 //}
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             beg2++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         beg1++;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // now divide the force constant for each contribution to the torsion energy
    // RDKit❗❌:     // about this bond by the number of contribs about this bond:
    // RDKit❗❌:     for (auto chI = contribsHere.begin(); chI != contribsHere.end(); ++chI) {
    // RDKit❗❌:       (*chI)->scaleForceConstant(contribsHere.size());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION UFF::Tools::addTorsions

    // PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    // PRECONDITION(field, "bad forcefield");
    // The mutable reference makes the source non-null field precondition
    // unrepresentable here.

    let matches = torsion_bond_matches(topology, rings, valence, torsion_bond_smarts)
        .map_err(UffBuilderError::TorsionBondQuery)?;
    for matched in matches {
        let [begin_atom_index, end_atom_index] = matched.as_slice() else {
            return Err(UffBuilderError::TorsionBondQuery(
                TorsionBondQueryError::MatchArity {
                    smarts: torsion_bond_smarts.to_owned(),
                    actual: matched.len(),
                },
            ));
        };
        torsions_for_bond(topology, params, *begin_atom_index, *end_atom_index, field)?;
    }

    // Behavior marker — RDKit❗❌: preserves parameter/query/error ordering,
    // match endpoint order, and the immediate per-match append/error boundary.
    // Complexity marker — RDKit❗❌: atom-only projection avoids one per-hit
    // bond map; Search still retains separate raw query/target map rows and the
    // source-shaped nested adjacency cost recorded by B25.
    Ok(())
}

fn angle_order(
    hybridization: Hybridization,
    atom_i: AtomId,
    atom_j: AtomId,
    atom_k: AtomId,
    rings: &RingInfo,
) -> u32 {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addAngles order switch (Builder.cpp:181-230)
    // RDKit✔️✔️:           int order = 0;
    // RDKit✔️✔️:           switch (atomJ->getHybridization()) {
    let order = match hybridization {
        // RDKit✔️✔️:             case Atom::SP:
        // RDKit✔️✔️:               order = 1;
        // RDKit✔️✔️:               break;
        Hybridization::Sp => 1,

        // RDKit✔️✔️:             case Atom::SP2:
        // RDKit✔️✔️:               order = 3;
        Hybridization::Sp2 => {
            let mut order = 3;

            // RDKit✔️✔️:               // the following is a hack to get decent geometries
            // RDKit✔️✔️:               // with 3- and 4-membered rings incorporating sp2 atoms
            // RDKit✔️✔️:               // if the central atom is in a ring of size 3
            // RDKit✔️✔️:               if (rings->isAtomInRingOfSize(j, 3)) {
            if rings.is_atom_in_ring_of_size(atom_j, 3) {
                // RDKit✔️✔️:                 // if the central atom and one of the bonded atoms, but not the
                // RDKit✔️✔️:                 //  other one are inside a ring, then this angle is between a
                // RDKit✔️✔️:                 // ring substituent and a ring edge
                // RDKit✔️✔️:                 if (rings->isAtomInRingOfSize(i, 3) !=
                // RDKit✔️✔️:                     rings->isAtomInRingOfSize(k, 3)) {
                if rings.is_atom_in_ring_of_size(atom_i, 3)
                    != rings.is_atom_in_ring_of_size(atom_k, 3)
                {
                    // RDKit✔️✔️:                   order = 30;
                    order = 30;
                    // RDKit✔️✔️:                 }
                    // RDKit✔️✔️:                 // if all atoms are inside the ring, then this is one of ring
                    // RDKit✔️✔️:                 // angles
                    // RDKit✔️✔️:                 else if (rings->isAtomInRingOfSize(i, 3) &&
                    // RDKit✔️✔️:                          rings->isAtomInRingOfSize(k, 3)) {
                } else if rings.is_atom_in_ring_of_size(atom_i, 3)
                    && rings.is_atom_in_ring_of_size(atom_k, 3)
                {
                    // RDKit✔️✔️:                   order = 35;
                    order = 35;
                    // RDKit✔️✔️:                 }
                }
            // RDKit✔️✔️:               }
            // RDKit✔️✔️:               // if the central atom is in a ring of size 4
            // RDKit✔️✔️:               else if (rings->isAtomInRingOfSize(j, 4)) {
            } else if rings.is_atom_in_ring_of_size(atom_j, 4) {
                // RDKit✔️✔️:                 // if the central atom and one of the bonded atoms, but not the
                // RDKit✔️✔️:                 //  other one are inside a ring, then this angle is between a
                // RDKit✔️✔️:                 // ring substituent and a ring edge
                // RDKit✔️✔️:                 if (rings->isAtomInRingOfSize(i, 4) !=
                // RDKit✔️✔️:                     rings->isAtomInRingOfSize(k, 4)) {
                if rings.is_atom_in_ring_of_size(atom_i, 4)
                    != rings.is_atom_in_ring_of_size(atom_k, 4)
                {
                    // RDKit✔️✔️:                   order = 40;
                    order = 40;
                    // RDKit✔️✔️:                 }
                    // RDKit✔️✔️:                 // if all atoms are inside the ring, then this is one of ring
                    // RDKit✔️✔️:                 // angles
                    // RDKit✔️✔️:                 else if (rings->isAtomInRingOfSize(i, 4) &&
                    // RDKit✔️✔️:                          rings->isAtomInRingOfSize(k, 4)) {
                } else if rings.is_atom_in_ring_of_size(atom_i, 4)
                    && rings.is_atom_in_ring_of_size(atom_k, 4)
                {
                    // RDKit✔️✔️:                   order = 45;
                    order = 45;
                    // RDKit✔️✔️:                 }
                }
            }
            // RDKit✔️✔️:               }
            // RDKit✔️✔️:               // end of the hack
            // RDKit✔️✔️:               break;
            order
        }

        // RDKit✔️✔️:             case Atom::SP3D2:
        // RDKit✔️✔️:               order = 4;
        // RDKit✔️✔️:               break;
        Hybridization::Sp3d2 => 4,

        // RDKit✔️✔️:             default:
        // RDKit✔️✔️:               order = 0;
        // RDKit✔️✔️:               break;
        _ => 0,
    };
    // RDKit✔️✔️:           }
    // END RDKIT CPP FUNCTION UFF::Tools::addAngles order switch
    // Complexity: only SP2 performs ring membership lookups (at most five);
    // each lookup scans one atom's ring-membership row and allocates nothing.
    order
}

fn select_tbp_axial<P: AsRef<[f64]>>(
    center_atom_index: usize,
    neighbors: &[NeighborRef],
    positions: &[P],
) -> Result<(BondId, BondId), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles axial selection (Builder.cpp:263-293)
    // RDKit❗✔️:   double mostNeg = 100.0;
    let mut most_neg = 100.0;
    // RDKit❗✔️:   ROMol::OEDGE_ITER beg1, end1;
    // RDKit❗✔️:   boost::tie(beg1, end1) = mol.getAtomBonds(atom);
    // RDKit❗✔️:   unsigned int aid = atom->getIdx();
    let center_row = positions[center_atom_index].as_ref();
    let center_point = Point3 {
        x: center_row[0],
        y: center_row[1],
        z: center_row[2],
    };
    let mut axial_pair = None;
    // RDKit❗✔️:   while (beg1 != end1) {
    // RDKit❗✔️:     const Bond *bond1 = mol[*beg1];
    // RDKit❗✔️:     unsigned int oaid = bond1->getOtherAtomIdx(aid);
    // RDKit❗✔️:     RDGeom::Point3D v1 =
    // RDKit❗✔️:         conf.getAtomPos(aid).directionVector(conf.getAtomPos(oaid));
    for neighbor1 in neighbors {
        let neighbor1_row = positions[neighbor1.atom_index].as_ref();
        let neighbor1_point = Point3 {
            x: neighbor1_row[0],
            y: neighbor1_row[1],
            z: neighbor1_row[2],
        };
        let v1 = direction_vector(&center_point, &neighbor1_point).map_err(|_| {
            UffBuilderError::SourceDirectionVectorBelowTolerance {
                center_atom_index,
                neighbor_atom_index: neighbor1.atom_index,
            }
        })?;

        // RDKit❗✔️:     ROMol::OEDGE_ITER beg2, end2;
        // RDKit❗✔️:     boost::tie(beg2, end2) = mol.getAtomBonds(atom);
        // RDKit❗✔️:     while (beg2 != end2) {
        // RDKit❗✔️:       const Bond *bond2 = mol[*beg2];
        for neighbor2 in neighbors {
            // RDKit❗✔️:       if (bond2->getIdx() > bond1->getIdx()) {
            if neighbor2.bond > neighbor1.bond {
                // RDKit❗✔️:         unsigned int oaid2 = bond2->getOtherAtomIdx(aid);
                // RDKit❗✔️:         RDGeom::Point3D v2 =
                // RDKit❗✔️:             conf.getAtomPos(aid).directionVector(conf.getAtomPos(oaid2));
                let neighbor2_row = positions[neighbor2.atom_index].as_ref();
                let neighbor2_point = Point3 {
                    x: neighbor2_row[0],
                    y: neighbor2_row[1],
                    z: neighbor2_row[2],
                };
                let v2 = direction_vector(&center_point, &neighbor2_point).map_err(|_| {
                    UffBuilderError::SourceDirectionVectorBelowTolerance {
                        center_atom_index,
                        neighbor_atom_index: neighbor2.atom_index,
                    }
                })?;
                // RDKit❗✔️:         double dot = v1.dotProduct(v2);
                let dot = v1.dot_product(&v2);
                // RDKit❗✔️:         if (dot < mostNeg) {
                if dot < most_neg {
                    // RDKit❗✔️:           mostNeg = dot;
                    most_neg = dot;
                    // RDKit❗✔️:           ax1 = bond1;
                    // RDKit❗✔️:           ax2 = bond2;
                    axial_pair = Some((neighbor1.bond, neighbor2.bond));
                    // RDKit❗✔️:         }
                }
            }
            // RDKit❗✔️:       ++beg2;
        }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     ++beg1;
    }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   CHECK_INVARIANT(ax1, "axial bond not found");
    // RDKit❗✔️:   CHECK_INVARIANT(ax2, "axial bond not found");
    let (axial1, axial2) =
        axial_pair.ok_or(UffBuilderError::SourceTbpAxialBondNotFound { center_atom_index })?;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles axial selection
    // Behavior marker — RDKit❗✔️: B16 fixed axial-pair regressions follow.
    // Complexity marker — RDKit✔️✔️: nested source-order scans and direct
    // BondId comparison retain O(d^2) work; coordinates and vectors stay borrowed/stack-local.
    Ok((axial1, axial2))
}

fn select_tbp_equatorial(
    center_atom_index: usize,
    neighbors: &[NeighborRef],
    axial_bonds: (BondId, BondId),
) -> Result<[BondId; 3], UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles equatorial selection (Builder.cpp:295-313)
    // RDKit✔️✔️:   boost::tie(beg1, end1) = mol.getAtomBonds(atom);
    // RDKit✔️✔️:   while (beg1 != end1) {
    // RDKit✔️✔️:     const Bond *bond = mol[*beg1];
    // RDKit✔️✔️:     ++beg1;
    // RDKit✔️✔️:     if (bond == ax1 || bond == ax2) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!eq1) {
    // RDKit✔️✔️:       eq1 = bond;
    // RDKit✔️✔️:     } else if (!eq2) {
    // RDKit✔️✔️:       eq2 = bond;
    // RDKit✔️✔️:     } else if (!eq3) {
    // RDKit✔️✔️:       eq3 = bond;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut equatorial = [None; 3];
    for neighbor in neighbors {
        if neighbor.bond == axial_bonds.0 || neighbor.bond == axial_bonds.1 {
            continue;
        }
        if equatorial[0].is_none() {
            equatorial[0] = Some(neighbor.bond);
        } else if equatorial[1].is_none() {
            equatorial[1] = Some(neighbor.bond);
        } else if equatorial[2].is_none() {
            equatorial[2] = Some(neighbor.bond);
        }
    }

    // RDKit✔️✔️:   CHECK_INVARIANT(eq1, "equatorial bond not found");
    let eq1 = equatorial[0].ok_or(UffBuilderError::SourceTbpEquatorialBondNotFound {
        center_atom_index,
        role: 1,
    })?;
    // RDKit✔️✔️:   CHECK_INVARIANT(eq2, "equatorial bond not found");
    let eq2 = equatorial[1].ok_or(UffBuilderError::SourceTbpEquatorialBondNotFound {
        center_atom_index,
        role: 2,
    })?;
    // RDKit✔️✔️:   CHECK_INVARIANT(eq3, "equatorial bond not found");
    let eq3 = equatorial[2].ok_or(UffBuilderError::SourceTbpEquatorialBondNotFound {
        center_atom_index,
        role: 3,
    })?;
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles equatorial selection
    // Behavior marker — RDKit✔️✔️: fixed regressions cover all source roles,
    // pair exclusions, source order, and malformed direct-slice branches.
    // Complexity marker — RDKit✔️✔️: one source-order pass and three fixed
    // stack slots; no collection allocation or sorting.
    Ok([eq1, eq2, eq3])
}

fn append_tbp_axial_angle(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    center_atom_index: usize,
    axial_neighbors: [NeighborRef; 2],
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles Axial-Axial (Builder.cpp:321-329)
    // RDKit❗✔️:   // Axial-Axial
    // RDKit❗✔️:   i = ax1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = ax2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax1->getBondTypeAsDouble(),
    // RDKit❗✔️:         ax2->getBondTypeAsDouble(), params[i], params[atomIdx], params[j], 2);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles Axial-Axial

    // BEGIN RDKIT CPP HELPER RDKit::Bond::getOtherAtomIdx (Bond.cpp:80-88)
    // RDKit❗✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit❗✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_endAtomIdx;
    // RDKit❗✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_beginAtomIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit❗✔️:   // above, but no need to be redundant, so just throw.
    // RDKit❗✔️:   POSTCONDITION(false, "bad index");
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDKit::Bond::getOtherAtomIdx
    // The validated topology adjacency already resolves getOtherAtomIdx into
    // each NeighborRef; preserve the selected axial pair order without rescanning.
    let axial1 = axial_neighbors[0];
    let axial2 = axial_neighbors[1];
    // RDKit❗✔️:   if (params[i] && params[j]) {
    let Some(params_i) = params[axial1.atom_index] else {
        return Ok(());
    };
    let Some(params_j) = params[axial2.atom_index] else {
        return Ok(());
    };
    let center_params = params[center_atom_index]
        .as_ref()
        .ok_or(UffBuilderError::SourceTbpCenterParamsMissing { center_atom_index })?;

    let bond_order_ij = bond_type_as_double(topology.bonds[axial1.bond.index()].order())
        .map_err(UffBuilderError::Valence)?;
    let bond_order_jk = bond_type_as_double(topology.bonds[axial2.bond.index()].order())
        .map_err(UffBuilderError::Valence)?;
    let source_i = u32::try_from(axial1.atom_index).map_err(|_| {
        UffBuilderError::SourceAngleIndexOverflow {
            atom_index: axial1.atom_index,
        }
    })?;
    let source_center = u32::try_from(center_atom_index).map_err(|_| {
        UffBuilderError::SourceAngleIndexOverflow {
            atom_index: center_atom_index,
        }
    })?;
    let source_j = u32::try_from(axial2.atom_index).map_err(|_| {
        UffBuilderError::SourceAngleIndexOverflow {
            atom_index: axial2.atom_index,
        }
    })?;
    let contrib = AngleBendContrib::new(
        field.positions(),
        source_i,
        source_center,
        source_j,
        bond_order_ij,
        bond_order_jk,
        params_i,
        center_params,
        params_j,
        2,
    )
    .map_err(|error| UffBuilderError::ForceFieldKernel(ForceFieldKernelError::from(error)))?;
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    field.add_contribution(Box::new(contrib));
    // Behavior marker — RDKit❗✔️: fixed endpoint-mask, center-precondition,
    // kernel-value, and append-order regressions follow in B18.
    // Complexity marker — RDKit✔️✔️: two direct endpoint/bond lookups, borrowed
    // parameters, one canonical constructor, one Box allocation and Vec append.
    Ok(())
}

fn append_tbp_equatorial_angles(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    center_atom_index: usize,
    equatorial_neighbors: [NeighborRef; 3],
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles Equatorial-Equatorial (Builder.cpp:330-354)
    // RDKit❗✔️:   // Equatorial-Equatorial
    // RDKit❗✔️:   i = eq1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, eq1->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq2->getBondTypeAsDouble(), params[i], params[atomIdx], params[j], 3);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = eq1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq3->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, eq1->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq3->getBondTypeAsDouble(), params[i], params[atomIdx], params[j], 3);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = eq2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq3->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, eq2->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq3->getBondTypeAsDouble(), params[i], params[atomIdx], params[j], 3);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles Equatorial-Equatorial

    // BEGIN RDKIT CPP HELPER RDKit::Bond::getOtherAtomIdx (Bond.cpp:80-88)
    // RDKit❗✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit❗✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_endAtomIdx;
    // RDKit❗✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_beginAtomIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit❗✔️:   // above, but no need to be redundant, so just throw.
    // RDKit❗✔️:   POSTCONDITION(false, "bad index");
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDKit::Bond::getOtherAtomIdx
    // Each NeighborRef is the validated adjacency result of getOtherAtomIdx;
    // preserve source role order without rescanning bonds.
    let [eq1, eq2, eq3] = equatorial_neighbors;

    let mut append_pair = |first: NeighborRef, second: NeighborRef| {
        let Some(first_params) = params[first.atom_index] else {
            return Ok(());
        };
        let Some(second_params) = params[second.atom_index] else {
            return Ok(());
        };

        let first_bond_order = bond_type_as_double(topology.bonds[first.bond.index()].order())
            .map_err(UffBuilderError::Valence)?;
        let second_bond_order = bond_type_as_double(topology.bonds[second.bond.index()].order())
            .map_err(UffBuilderError::Valence)?;
        let center_params = params[center_atom_index]
            .as_ref()
            .ok_or(UffBuilderError::SourceTbpCenterParamsMissing { center_atom_index })?;
        let first_atom_index = u32::try_from(first.atom_index).map_err(|_| {
            UffBuilderError::SourceAngleIndexOverflow {
                atom_index: first.atom_index,
            }
        })?;
        let center_index = u32::try_from(center_atom_index).map_err(|_| {
            UffBuilderError::SourceAngleIndexOverflow {
                atom_index: center_atom_index,
            }
        })?;
        let second_atom_index = u32::try_from(second.atom_index).map_err(|_| {
            UffBuilderError::SourceAngleIndexOverflow {
                atom_index: second.atom_index,
            }
        })?;
        let contribution = AngleBendContrib::new(
            field.positions(),
            first_atom_index,
            center_index,
            second_atom_index,
            first_bond_order,
            second_bond_order,
            first_params,
            center_params,
            second_params,
            3,
        )
        .map_err(|error| UffBuilderError::ForceFieldKernel(ForceFieldKernelError::from(error)))?;
        field.add_contribution(Box::new(contribution));
        Ok(())
    };

    append_pair(eq1, eq2)?;
    append_pair(eq1, eq3)?;
    append_pair(eq2, eq3)?;
    // Behavior marker — RDKit✔️✔️: all eight endpoint masks, the first active
    // pair failure, a later typed failure, and fixed real-kernel energy/gradient
    // plus source-order-sensitive summation pass in B19 regressions.
    // Complexity marker — RDKit✔️✔️: three fixed source-order checks, direct
    // adjacency/bond lookups, one allocation per surviving source term.
    Ok(())
}

fn append_tbp_mixed_angles(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    center_atom_index: usize,
    axial_neighbors: [NeighborRef; 2],
    equatorial_neighbors: [NeighborRef; 3],
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles Axial-Equatorial (Builder.cpp:356-405)
    // RDKit❗✔️:   // Axial-Equatorial
    // RDKit❗✔️:   i = ax1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax1->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq1->getBondTypeAsDouble(), params[i], params[atomIdx], params[j]);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = ax1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax1->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq2->getBondTypeAsDouble(), params[i], params[atomIdx], params[j]);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = ax1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq3->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax1->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq3->getBondTypeAsDouble(), params[i], params[atomIdx], params[j]);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = ax2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq1->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax2->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq1->getBondTypeAsDouble(), params[i], params[atomIdx], params[j]);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = ax2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax2->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq2->getBondTypeAsDouble(), params[i], params[atomIdx], params[j]);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   i = ax2->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   j = eq3->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:   if (params[i] && params[j]) {
    // RDKit❗✔️:     contrib = new AngleBendContrib(
    // RDKit❗✔️:         field, i, atomIdx, j, ax2->getBondTypeAsDouble(),
    // RDKit❗✔️:         eq3->getBondTypeAsDouble(), params[i], params[atomIdx], params[j]);
    // RDKit❗✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles Axial-Equatorial

    // BEGIN RDKIT CPP HELPER RDKit::Bond::getOtherAtomIdx (Bond.cpp:80-88)
    // RDKit❗✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit❗✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_endAtomIdx;
    // RDKit❗✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_beginAtomIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit❗✔️:   // above, but no need to be redundant, so just throw.
    // RDKit❗✔️:   POSTCONDITION(false, "bad index");
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDKit::Bond::getOtherAtomIdx

    // BEGIN RDKIT CPP DEFAULT ARGUMENT ForceFields::UFF::AngleBendContrib::AngleBendContrib (AngleBend.h:46-49)
    // RDKit❗✔️:   AngleBendContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
    // RDKit❗✔️:                    unsigned int idx3, double bondOrder12, double bondOrder23,
    // RDKit❗✔️:                    const AtomicParams *at1Params, const AtomicParams *at2Params,
    // RDKit❗✔️:                    const AtomicParams *at3Params, unsigned int order = 0);
    // END RDKIT CPP DEFAULT ARGUMENT ForceFields::UFF::AngleBendContrib::AngleBendContrib

    // The validated NeighborRef values are the source Bond::getOtherAtomIdx
    // results; keep axial and equatorial role order without rescanning.
    let [axial1, axial2] = axial_neighbors;
    let [equatorial1, equatorial2, equatorial3] = equatorial_neighbors;

    let mut append_pair = |first: NeighborRef, second: NeighborRef| {
        let Some(first_params) = params[first.atom_index] else {
            return Ok(());
        };
        let Some(second_params) = params[second.atom_index] else {
            return Ok(());
        };

        let first_bond_order = bond_type_as_double(topology.bonds[first.bond.index()].order())
            .map_err(UffBuilderError::Valence)?;
        let second_bond_order = bond_type_as_double(topology.bonds[second.bond.index()].order())
            .map_err(UffBuilderError::Valence)?;
        let center_params = params[center_atom_index]
            .as_ref()
            .ok_or(UffBuilderError::SourceTbpCenterParamsMissing { center_atom_index })?;
        let first_atom_index = u32::try_from(first.atom_index).map_err(|_| {
            UffBuilderError::SourceAngleIndexOverflow {
                atom_index: first.atom_index,
            }
        })?;
        let center_index = u32::try_from(center_atom_index).map_err(|_| {
            UffBuilderError::SourceAngleIndexOverflow {
                atom_index: center_atom_index,
            }
        })?;
        let second_atom_index = u32::try_from(second.atom_index).map_err(|_| {
            UffBuilderError::SourceAngleIndexOverflow {
                atom_index: second.atom_index,
            }
        })?;
        let contribution = AngleBendContrib::new(
            field.positions(),
            first_atom_index,
            center_index,
            second_atom_index,
            first_bond_order,
            second_bond_order,
            first_params,
            center_params,
            second_params,
            0,
        )
        .map_err(|error| UffBuilderError::ForceFieldKernel(ForceFieldKernelError::from(error)))?;
        field.add_contribution(Box::new(contribution));
        Ok(())
    };

    append_pair(axial1, equatorial1)?;
    append_pair(axial1, equatorial2)?;
    append_pair(axial1, equatorial3)?;
    append_pair(axial2, equatorial1)?;
    append_pair(axial2, equatorial2)?;
    append_pair(axial2, equatorial3)?;
    // Behavior marker — RDKit✔️✔️: all 32 terminal masks, typed first and later
    // failures, source append order, and real kernel energy/gradient outputs pass.
    // Complexity marker — RDKit✔️✔️: six fixed pair calls use borrowed role
    // arrays and direct bond lookups; each active pair allocates one Box.
    Ok(())
}

fn add_trigonal_bipyramid_angles(
    topology: &TopologyBlock,
    center_atom_index: usize,
    params: &UffParamsByAtom<'_>,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles (Builder.cpp:248-405)
    // RDKit✔️✔️: void addTrigonalBipyramidAngles(const Atom *atom, const ROMol &mol, int confId,
    // RDKit✔️✔️:                                 const AtomicParamVect &params,
    // RDKit✔️✔️:                                 ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    let atom = topology
        .atoms
        .get(center_atom_index)
        .ok_or(UffBuilderError::SourceTbpAtomPrecondition { center_atom_index })?;
    // RDKit✔️✔️:   PRECONDITION(atom->getHybridization() == Atom::SP3D, "bad hybridization");
    if atom.hybridization() != Hybridization::Sp3d {
        return Err(UffBuilderError::SourceTbpHybridizationPrecondition {
            center_atom_index,
            actual: atom.hybridization(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(atom->getDegree() == 5, "bad degree");
    let neighbors = topology.adjacency.neighbors_of(center_atom_index);
    if neighbors.len() != 5 {
        return Err(UffBuilderError::SourceTbpDegreePrecondition {
            center_atom_index,
            actual_degree: neighbors.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(field, "bad forcefield");
    // A mutable reference is the non-null field precondition.

    // RDKit✔️✔️:   const Bond *ax1 = nullptr, *ax2 = nullptr;
    // RDKit✔️✔️:   const Bond *eq1 = nullptr, *eq2 = nullptr, *eq3 = nullptr;
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   double mostNeg = 100.0;
    // Axial selection borrows the already selected 3D rows from the same
    // field later used by the angle contributions; no coordinate block copy.
    let axial_bonds = select_tbp_axial(center_atom_index, neighbors, field.positions())?;
    // RDKit✔️✔️:   CHECK_INVARIANT(ax1, "axial bond not found");
    // RDKit✔️✔️:   CHECK_INVARIANT(ax2, "axial bond not found");

    // RDKit✔️✔️:   if (bond == ax1 || bond == ax2) { continue; }
    // RDKit✔️✔️:   if (!eq1) { eq1 = bond; }
    // RDKit✔️✔️:   else if (!eq2) { eq2 = bond; }
    // RDKit✔️✔️:   else if (!eq3) { eq3 = bond; }
    let equatorial_bonds = select_tbp_equatorial(center_atom_index, neighbors, axial_bonds)?;
    // RDKit✔️✔️:   CHECK_INVARIANT(eq1, "equatorial bond not found");
    // RDKit✔️✔️:   CHECK_INVARIANT(eq2, "equatorial bond not found");
    // RDKit✔️✔️:   CHECK_INVARIANT(eq3, "equatorial bond not found");

    // The validated incident rows resolve source Bond pointers by identity;
    // retaining the returned order does not assume adjacency order == BondId.
    let neighbor_for_bond = |bond: BondId| {
        neighbors
            .iter()
            .find(|neighbor| neighbor.bond == bond)
            .copied()
    };
    let axial_neighbors = [
        neighbor_for_bond(axial_bonds.0)
            .ok_or(UffBuilderError::SourceTbpAxialBondNotFound { center_atom_index })?,
        neighbor_for_bond(axial_bonds.1)
            .ok_or(UffBuilderError::SourceTbpAxialBondNotFound { center_atom_index })?,
    ];
    let equatorial_neighbors = [
        neighbor_for_bond(equatorial_bonds[0]).ok_or(
            UffBuilderError::SourceTbpEquatorialBondNotFound {
                center_atom_index,
                role: 1,
            },
        )?,
        neighbor_for_bond(equatorial_bonds[1]).ok_or(
            UffBuilderError::SourceTbpEquatorialBondNotFound {
                center_atom_index,
                role: 2,
            },
        )?,
        neighbor_for_bond(equatorial_bonds[2]).ok_or(
            UffBuilderError::SourceTbpEquatorialBondNotFound {
                center_atom_index,
                role: 3,
            },
        )?,
    ];

    // RDKit✔️✔️:   // Axial-Axial
    append_tbp_axial_angle(topology, params, center_atom_index, axial_neighbors, field)?;
    // RDKit✔️✔️:   // Equatorial-Equatorial
    append_tbp_equatorial_angles(
        topology,
        params,
        center_atom_index,
        equatorial_neighbors,
        field,
    )?;
    // RDKit✔️✔️:   // Axial-Equatorial
    append_tbp_mixed_angles(
        topology,
        params,
        center_atom_index,
        axial_neighbors,
        equatorial_neighbors,
        field,
    )?;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addTrigonalBipyramidAngles
    // Behavior marker — RDKit✔️✔️: B21's composed fixed-regression gate passes.
    // Complexity marker — RDKit✔️✔️: source O(d^2) pair selection, O(d) role
    // resolution and ten fixed branches; borrowed rows and no term-list copy.
    Ok(())
}

fn inversion_center(
    topology: &TopologyBlock,
    center_atom_index: usize,
) -> Option<(u8, [usize; 4], bool)> {
    // BEGIN RDKIT CPP HELPER UFF::Tools::addInversions center selection (Builder.cpp:606-638)
    // RDKit❗✔️:     int at2AtomicNum = atom[1]->getAtomicNum();
    // RDKit❗✔️:     // if the central atom is not carbon, nitrogen, oxygen,
    // RDKit❗✔️:     // phosphorous, arsenic, antimonium or bismuth, skip it
    // RDKit❗✔️:     if (((at2AtomicNum != 6) && (at2AtomicNum != 7) && (at2AtomicNum != 8) &&
    // RDKit❗✔️:          (at2AtomicNum != 15) && (at2AtomicNum != 33) && (at2AtomicNum != 51) &&
    // RDKit❗✔️:          (at2AtomicNum != 83)) ||
    // RDKit❗✔️:         (atom[1]->getDegree() != 3)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // if the central atom is carbon, nitrogen or oxygen
    // RDKit❗✔️:     // but hybridization is not sp2, skip it
    // RDKit❗✔️:     if (((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) &&
    // RDKit❗✔️:         (atom[1]->getHybridization() != Atom::SP2)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
    // RDKit❗✔️:     unsigned int i = 0;
    // RDKit❗✔️:     bool isBoundToSP2O = false;
    // RDKit❗✔️:     for (; nbrIdx != endNbrs; ++nbrIdx) {
    // RDKit❗✔️:       atom[i] = mol[*nbrIdx];
    // RDKit❗✔️:       idx[i] = atom[i]->getIdx();
    // RDKit❗✔️:       // if the central atom is sp2 carbon and is
    // RDKit❗✔️:       // bound to sp2 oxygen, set a flag
    // RDKit❗✔️:       if (!isBoundToSP2O) {
    // RDKit❗✔️:         isBoundToSP2O =
    // RDKit❗✔️:             ((at2AtomicNum == 6) && (atom[i]->getAtomicNum() == 8) &&
    // RDKit❗✔️:              (atom[i]->getHybridization() == Atom::SP2));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (!i) {
    // RDKit❗✔️:         ++i;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       ++i;
    // RDKit❗✔️:     }

    let center_atom = &topology.atoms[center_atom_index];
    let at2_atomic_num = center_atom.atomic_number();
    if !matches!(at2_atomic_num, 6 | 7 | 8 | 15 | 33 | 51 | 83) {
        return None;
    }

    let neighbors = topology.adjacency.neighbors_of(center_atom_index);
    if neighbors.len() != 3 {
        return None;
    }

    if matches!(at2_atomic_num, 6 | 7 | 8) && center_atom.hybridization() != Hybridization::Sp2 {
        return None;
    }

    let mut idx = [0_usize; 4];
    idx[1] = center_atom_index;
    let mut i = 0;
    let mut is_bound_to_sp2_o = false;
    for neighbor in neighbors {
        let neighbor_atom = &topology.atoms[neighbor.atom_index];
        idx[i] = neighbor.atom_index;
        if !is_bound_to_sp2_o {
            is_bound_to_sp2_o = at2_atomic_num == 6
                && neighbor_atom.atomic_number() == 8
                && neighbor_atom.hybridization() == Hybridization::Sp2;
        }
        if i == 0 {
            i += 1;
        }
        i += 1;
    }

    // Behavior marker — RDKit❗✔️: preserve the exact allowed elements,
    // degree/hybridization filters, neighbor-slot order, and carbon-only
    // SP2-oxygen flag; B27's fixed branch matrix is the next gate.
    // Complexity marker — RDKit✔️✔️: one constant-time filter sequence and
    // one borrowed scan of exactly three neighbors; fixed stack storage and
    // no allocation, cloning, or chemistry recomputation.
    // END RDKIT CPP HELPER UFF::Tools::addInversions center selection
    Some((at2_atomic_num, idx, is_bound_to_sp2_o))
}

pub(super) fn add_inversions(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addInversions (Builder.cpp:593-666)
    // RDKit❗✔️: void addInversions(const ROMol &mol, const AtomicParamVect &params,
    // RDKit❗✔️:                    ForceFields::ForceField *field) {
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    // RDKit❗✔️:   PRECONDITION(field, "bad forcefield");
    // The mutable field reference preserves the source non-null precondition.

    // RDKit❗✔️:   unsigned int idx[4];
    // RDKit❗✔️:   unsigned int n[4];
    // RDKit❗✔️:   const Atom *atom[4];
    // RDKit❗✔️:   ROMol::ADJ_ITER nbrIdx;
    // RDKit❗✔️:   ROMol::ADJ_ITER endNbrs;
    // RDKit❗✔️:   for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
    for center_atom_index in 0..topology.atoms.len() {
        let Some((at2_atomic_num, idx, is_bound_to_sp2_o)) =
            inversion_center(topology, center_atom_index)
        else {
            continue;
        };

        // RDKit❗✔️:     for (unsigned int i = 0; i < 3; ++i) {
        // RDKit❗✔️:       n[1] = 1;
        // RDKit❗✔️:       switch (i) {
        // RDKit❗✔️:         case 0:
        // RDKit❗✔️:           n[0] = 0;
        // RDKit❗✔️:           n[2] = 2;
        // RDKit❗✔️:           n[3] = 3;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:         case 1:
        // RDKit❗✔️:           n[0] = 0;
        // RDKit❗✔️:           n[2] = 3;
        // RDKit❗✔️:           n[3] = 2;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:         case 2:
        // RDKit❗✔️:           n[0] = 2;
        // RDKit❗✔️:           n[2] = 3;
        // RDKit❗✔️:           n[3] = 0;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:       }
        for [idx1, idx2, idx3, idx4] in [
            [idx[0], idx[1], idx[2], idx[3]],
            [idx[0], idx[1], idx[3], idx[2]],
            [idx[2], idx[1], idx[3], idx[0]],
        ] {
            let idx1 =
                u32::try_from(idx1).map_err(|_| UffBuilderError::SourceInversionIndexOverflow {
                    atom_index: idx1 as usize,
                })?;
            let idx2 =
                u32::try_from(idx2).map_err(|_| UffBuilderError::SourceInversionIndexOverflow {
                    atom_index: idx2 as usize,
                })?;
            let idx3 =
                u32::try_from(idx3).map_err(|_| UffBuilderError::SourceInversionIndexOverflow {
                    atom_index: idx3 as usize,
                })?;
            let idx4 =
                u32::try_from(idx4).map_err(|_| UffBuilderError::SourceInversionIndexOverflow {
                    atom_index: idx4 as usize,
                })?;

            // RDKit❗✔️:       InversionContrib *contrib;
            // RDKit❗✔️:       contrib = new InversionContrib(field, idx[n[0]], idx[n[1]], idx[n[2]],
            // RDKit❗✔️:                                      idx[n[3]], at2AtomicNum, isBoundToSP2O);
            let contribution = InversionContrib::new(
                field.positions(),
                idx1,
                idx2,
                idx3,
                idx4,
                i32::from(at2_atomic_num),
                is_bound_to_sp2_o,
            )
            .map_err(UffBuilderError::InversionContribution)?;

            // RDKit❗✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
            field.add_contribution(Box::new(contribution));
        }
        // RDKit❗✔️:     }
    }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Behavior marker — RDKit❗✔️: preserve the length-first precondition,
    // source center order, no missing-row skip, three permutation order,
    // constructor error timing, and one append after each successful term.
    // Complexity marker — RDKit✔️✔️: O(V) center scan, three fixed neighbor
    // inspections per eligible center, and exactly three scalar allocations
    // and appends per eligible center; no temporary term vector or row clone.
    // END RDKIT CPP FUNCTION UFF::Tools::addInversions
    Ok(())
}

pub(super) fn add_angle_special_cases(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addAngleSpecialCases (Builder.cpp:412-430)
    // RDKit✔️✔️: void addAngleSpecialCases(const ROMol &mol, int confId,
    // RDKit✔️✔️:                           const AtomicParamVect &params,
    // RDKit✔️✔️:                           ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(field, "bad forcefield");
    // The existing mutable borrow is the non-null field precondition; selected
    // coordinate rows already live in this ForceField, so no confId is needed.

    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    for (atom_idx, atom) in topology.atoms.iter().enumerate() {
        // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(i);
        // RDKit✔️✔️:     // trigonal bipyramidal:
        // RDKit✔️✔️:     if ((atom->getHybridization() == Atom::SP3D && atom->getDegree() == 5)) {
        if atom.hybridization() == Hybridization::Sp3d
            && topology.adjacency.neighbors_of(atom_idx).len() == 5
        {
            // RDKit✔️✔️:       addTrigonalBipyramidAngles(atom, mol, confId, params, field);
            add_trigonal_bipyramid_angles(topology, atom_idx, params, field)?;
            // RDKit✔️✔️:     }
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addAngleSpecialCases
    // Behavior marker — RDKit✔️✔️: the exact dispatch matrix, error order,
    // successful two-center terms, and prefix retention pass fixed regressions.
    // Complexity marker — RDKit✔️✔️: one atom pass, O(1) stored-state/degree
    // checks, no adjacency rebuild or temporary allocation.
    Ok(())
}

pub(super) fn add_angles(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    rings: &RingInfo,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addAngles (Builder.cpp:138-241)
    // RDKit❗✔️: void addAngles(const ROMol &mol, const AtomicParamVect &params,
    // RDKit❗✔️:                ForceFields::ForceField *field) {
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    // RDKit❗✔️:   PRECONDITION(field, "bad forcefield");
    // A mutable reference is the non-null field precondition. RingInfo is
    // borrowed from preparation instead of looked up or recomputed here.
    // RDKit❗✔️:   ROMol::ADJ_ITER nbr1Idx;
    // RDKit❗✔️:   ROMol::ADJ_ITER end1Nbrs;
    // RDKit❗✔️:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit❗✔️:   ROMol::ADJ_ITER end2Nbrs;
    // RDKit❗✔️:   RingInfo *rings = mol.getRingInfo();
    // RDKit❗✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = topology.atoms.len();
    // RDKit❗✔️:   for (unsigned int j = 0; j < nAtoms; j++) {
    for j in 0..n_atoms {
        // RDKit❗✔️:     if (!params[j]) {
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
        let Some(center_params) = params[j] else {
            continue;
        };
        // RDKit❗✔️:     const Atom *atomJ = mol.getAtomWithIdx(j);
        let atom_j = &topology.atoms[j];
        // RDKit❗✔️:     if (atomJ->getDegree() == 1) {
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
        let neighbors = topology.adjacency.neighbors_of(j);
        if neighbors.len() == 1 {
            continue;
        }
        // RDKit❗✔️:     boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(atomJ);
        // RDKit❗✔️:     for (; nbr1Idx != end1Nbrs; nbr1Idx++) {
        for (nbr1_pos, nbr1) in neighbors.iter().enumerate() {
            // RDKit❗✔️:       const Atom *atomI = mol[*nbr1Idx];
            // RDKit❗✔️:       unsigned int i = atomI->getIdx();
            let i = nbr1.atom_index;
            // RDKit❗✔️:       if (!params[i]) {
            // RDKit❗✔️:         continue;
            // RDKit❗✔️:       }
            let Some(atom_i_params) = params[i] else {
                continue;
            };
            // RDKit❗✔️:       boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(atomJ);
            // RDKit❗✔️:       for (; nbr2Idx != end2Nbrs; nbr2Idx++) {
            for (nbr2_pos, nbr2) in neighbors.iter().enumerate() {
                // RDKit❗✔️:         if (nbr2Idx < (nbr1Idx + 1)) {
                // RDKit❗✔️:           continue;
                // RDKit❗✔️:         }
                if nbr2_pos < nbr1_pos + 1 {
                    continue;
                }
                // RDKit❗✔️:         const Atom *atomK = mol[*nbr2Idx];
                // RDKit❗✔️:         unsigned int k = atomK->getIdx();
                let k = nbr2.atom_index;
                // RDKit❗✔️:         if (!params[k]) {
                // RDKit❗✔️:           continue;
                // RDKit❗✔️:         }
                let Some(atom_k_params) = params[k] else {
                    continue;
                };
                // RDKit❗✔️:         // skip special cases:
                // RDKit❗✔️:         if (!(atomJ->getHybridization() == Atom::SP3D &&
                // RDKit❗✔️:               atomJ->getDegree() == 5)) {
                if !(atom_j.hybridization() == Hybridization::Sp3d && neighbors.len() == 5) {
                    // RDKit❗✔️:           const Bond *b1 = mol.getBondBetweenAtoms(i, j);
                    // RDKit❗✔️:           const Bond *b2 = mol.getBondBetweenAtoms(k, j);
                    let bond_ij = &topology.bonds[nbr1.bond.index()];
                    let bond_kj = &topology.bonds[nbr2.bond.index()];
                    // RDKit❗✔️:           // FIX: recognize amide bonds here.
                    // RDKit❗✔️:           AngleBendContrib *contrib;
                    // RDKit❗✔️:           int order = 0;
                    // RDKit❗✔️:           switch (atomJ->getHybridization()) {
                    let order = angle_order(
                        atom_j.hybridization(),
                        topology.atoms[i].id(),
                        atom_j.id(),
                        topology.atoms[k].id(),
                        rings,
                    );
                    // RDKit❗✔️:             case Atom::SP:
                    // RDKit❗✔️:               order = 1;
                    // RDKit❗✔️:               break;
                    // RDKit❗✔️:             case Atom::SP2:
                    // RDKit❗✔️:               order = 3;
                    // RDKit❗✔️:               // the following is a hack to get decent geometries
                    // RDKit❗✔️:               // with 3- and 4-membered rings incorporating sp2 atoms
                    // RDKit❗✔️:               // if the central atom is in a ring of size 3
                    // RDKit❗✔️:               if (rings->isAtomInRingOfSize(j, 3)) {
                    // RDKit❗✔️:                 // if the central atom and one of the bonded atoms, but not the
                    // RDKit❗✔️:                 //  other one are inside the ring, then this angle is between a
                    // RDKit❗✔️:                 // ring substituent and a ring edge
                    // RDKit❗✔️:                 if (rings->isAtomInRingOfSize(i, 3) !=
                    // RDKit❗✔️:                     rings->isAtomInRingOfSize(k, 3)) {
                    // RDKit❗✔️:                   order = 30;
                    // RDKit❗✔️:                 }
                    // RDKit❗✔️:                 // if all atoms are inside the ring, then this is one of ring
                    // RDKit❗✔️:                 // angles
                    // RDKit❗✔️:                 else if (rings->isAtomInRingOfSize(i, 3) &&
                    // RDKit❗✔️:                          rings->isAtomInRingOfSize(k, 3)) {
                    // RDKit❗✔️:                   order = 35;
                    // RDKit❗✔️:                 }
                    // RDKit❗✔️:               }
                    // RDKit❗✔️:               // if the central atom is in a ring of size 4
                    // RDKit❗✔️:               else if (rings->isAtomInRingOfSize(j, 4)) {
                    // RDKit❗✔️:                 // if the central atom and one of the bonded atoms, but not the
                    // RDKit❗✔️:                 //  other one are inside the ring, then this angle is between a
                    // RDKit❗✔️:                 // ring substituent and a ring edge
                    // RDKit❗✔️:                 if (rings->isAtomInRingOfSize(i, 4) !=
                    // RDKit❗✔️:                     rings->isAtomInRingOfSize(k, 4)) {
                    // RDKit❗✔️:                   order = 40;
                    // RDKit❗✔️:                 }
                    // RDKit❗✔️:                 // if all atoms are inside the ring, then this is one of ring
                    // RDKit❗✔️:                 // angles
                    // RDKit❗✔️:                 else if (rings->isAtomInRingOfSize(i, 4) &&
                    // RDKit❗✔️:                          rings->isAtomInRingOfSize(k, 4)) {
                    // RDKit❗✔️:                   order = 45;
                    // RDKit❗✔️:                 }
                    // RDKit❗✔️:               }
                    // RDKit❗✔️:               // end of the hack
                    // RDKit❗✔️:               break;
                    // RDKit❗✔️:             case Atom::SP3D2:
                    // RDKit❗✔️:               order = 4;
                    // RDKit❗✔️:               break;
                    // RDKit❗✔️:             default:
                    // RDKit❗✔️:               order = 0;
                    // RDKit❗✔️:               break;
                    // RDKit❗✔️:           }

                    // RDKit❗✔️:           contrib =
                    // RDKit❗✔️:               new AngleBendContrib(field, i, j, k, b1->getBondTypeAsDouble(),
                    // RDKit❗✔️:                                    b2->getBondTypeAsDouble(), params[i],
                    // RDKit❗✔️:                                    params[j], params[k], order);
                    let bond_order_ij =
                        bond_type_as_double(bond_ij.order()).map_err(UffBuilderError::Valence)?;
                    let bond_order_kj =
                        bond_type_as_double(bond_kj.order()).map_err(UffBuilderError::Valence)?;
                    let source_i = u32::try_from(i)
                        .map_err(|_| UffBuilderError::SourceAngleIndexOverflow { atom_index: i })?;
                    let source_j = u32::try_from(j)
                        .map_err(|_| UffBuilderError::SourceAngleIndexOverflow { atom_index: j })?;
                    let source_k = u32::try_from(k)
                        .map_err(|_| UffBuilderError::SourceAngleIndexOverflow { atom_index: k })?;
                    let contrib = AngleBendContrib::new(
                        field.positions(),
                        source_i,
                        source_j,
                        source_k,
                        bond_order_ij,
                        bond_order_kj,
                        atom_i_params,
                        center_params,
                        atom_k_params,
                        order,
                    )
                    .map_err(|error| {
                        UffBuilderError::ForceFieldKernel(ForceFieldKernelError::from(error))
                    })?;
                    // RDKit❗✔️:           field->contribs().push_back(ForceFields::ContribPtr(contrib));
                    field.add_contribution(Box::new(contrib));
                    // RDKit❗✔️:         }
                }
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:     }
        }
        // RDKit❗✔️:   }
    }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addAngles
    // Behavior marker — RDKit❗✔️: fixed B15 tests pass; C++ leaves relative
    // evaluation order unspecified when both bond-order arguments fail.
    // Complexity marker — RDKit✔️✔️: nested scans preserve source degree-pair work;
    // borrowed adjacency, parameters, ring state, and coordinates allocate no copies.
    Ok(())
}

pub(super) fn add_bonds(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    field: &mut ForceField<'_>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::addBonds (Builder.cpp:32-52)
    // RDKit✔️✔️: void addBonds(const ROMol &mol, const AtomicParamVect &params,
    // RDKit✔️✔️:               ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    // RDKit✔️✔️:   PRECONDITION(field, "bad forcefield");
    // The mutable reference makes RDKit's non-null field precondition explicit.

    // RDKit✔️✔️:   for (ROMol::ConstBondIterator bi = mol.beginBonds(); bi != mol.endBonds();
    // RDKit✔️✔️:        bi++) {
    for bond in &topology.bonds {
        // RDKit✔️✔️:     int idx1 = (*bi)->getBeginAtomIdx();
        let idx1 = bond.begin().index();
        // RDKit✔️✔️:     int idx2 = (*bi)->getEndAtomIdx();
        let idx2 = bond.end().index();

        // RDKit✔️✔️:     // FIX: recognize amide bonds here.

        // RDKit✔️✔️:     if (params[idx1] && params[idx2]) {
        let Some(end1_params) = params[idx1] else {
            continue;
        };
        let Some(end2_params) = params[idx2] else {
            continue;
        };

        // RDKit✔️✔️:       BondStretchContrib *contrib;
        // RDKit✔️✔️:       contrib = new BondStretchContrib(field, idx1, idx2,
        // RDKit✔️✔️:                                        (*bi)->getBondTypeAsDouble(),
        // RDKit✔️✔️:                                        params[idx1], params[idx2]);
        let bond_order = bond_type_as_double(bond.order()).map_err(UffBuilderError::Valence)?;
        let source_idx1 = u32::try_from(idx1)
            .map_err(|_| UffBuilderError::SourceBondIndexOverflow { atom_index: idx1 })?;
        let source_idx2 = u32::try_from(idx2)
            .map_err(|_| UffBuilderError::SourceBondIndexOverflow { atom_index: idx2 })?;
        let contrib = BondStretchContrib::new(
            field.positions(),
            source_idx1,
            source_idx2,
            bond_order,
            end1_params,
            end2_params,
        )
        .map_err(UffBuilderError::ForceFieldKernel)?;

        // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
        field.add_contribution(Box::new(contrib));
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::addBonds
    // Complexity: one source-order pass over bonds, constant-time endpoint
    // parameter/index handling, and one allocation plus append per emitted
    // contribution; parameter rows and force-field positions are borrowed.
    Ok(())
}

pub(super) fn two_bit_cell_pos(
    n_atoms: usize,
    i: usize,
    j: usize,
) -> Result<usize, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::twoBitCellPos (Builder.cpp:54-60)
    // RDKit❗✔️: unsigned int twoBitCellPos(unsigned int nAtoms, int i, int j) {
    // RDKit❗✔️:   if (j < i) {
    if n_atoms == 0 || i >= n_atoms || j >= n_atoms {
        return Err(UffBuilderError::NeighborMatrixIndexOutOfRange { n_atoms, i, j });
    }

    // Source inputs use unsigned int nAtoms and return an unsigned int cell
    // position. Reject unrepresentable atom counts before doing host arithmetic.
    let source_atom_count = u32::try_from(n_atoms)
        .map_err(|_| UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let source_atom_count = u64::from(source_atom_count);
    let source_atom_count_u32 = u32::try_from(source_atom_count)
        .map_err(|_| UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let next_atom_count = source_atom_count_u32
        .checked_add(1)
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let allocation_numerator = source_atom_count_u32
        .checked_mul(next_atom_count)
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let source_allocation_bytes = allocation_numerator
        .checked_sub(1)
        .map(|cells_eight| cells_eight / 8)
        .and_then(|bytes_before_last| bytes_before_last.checked_add(1))
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let _host_allocation_bytes = usize::try_from(source_allocation_bytes)
        .map_err(|_| UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;

    let source_cell_count = source_atom_count
        .checked_add(1)
        .and_then(|next| source_atom_count.checked_mul(next))
        .map(|product| product / 2)
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let _host_cell_count = usize::try_from(source_cell_count)
        .map_err(|_| UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    if source_cell_count > u64::from(u32::MAX) + 1 {
        return Err(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms });
    }

    let (mut row, mut column) = (i, j);
    if j < i {
        // RDKit❗✔️:     std::swap(i, j);
        std::mem::swap(&mut row, &mut column);
        // RDKit❗✔️:   }
    }

    // RDKit❗✔️:   return i * (nAtoms - 1) + i * (1 - i) / 2 + j;
    // Algebraically factor the source row offset as i*(2*nAtoms-i-1)/2.
    // Divide the even factor first so intermediates remain host-size checked.
    let row_span = n_atoms
        .checked_mul(2)
        .and_then(|twice_n| {
            row.checked_add(1)
                .and_then(|row_end| twice_n.checked_sub(row_end))
        })
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let row_start = if row % 2 == 0 {
        (row / 2).checked_mul(row_span)
    } else {
        row.checked_mul(row_span / 2)
    }
    .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let cell_pos = row_start
        .checked_add(column)
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;

    if u64::try_from(cell_pos).map_or(true, |position| position > u64::from(u32::MAX)) {
        return Err(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms });
    }

    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::twoBitCellPos
    // Complexity: constant-time validation and arithmetic, with no allocation
    // or per-pair scan; invalid host/source domains fail as typed input errors.
    Ok(cell_pos)
}

pub(super) fn set_two_bit_cell(
    res: &mut [u8],
    pos: usize,
    value: u8,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::setTwoBitCell (Builder.cpp:62-68)
    // RDKit❗✔️: void setTwoBitCell(boost::shared_array<std::uint8_t> &res, unsigned int pos,
    // RDKit❗✔️:                    std::uint8_t value) {
    // RDKit❗✔️:   unsigned int twoBitPos = pos / 4;
    let source_pos = u32::try_from(pos)
        .map_err(|_| UffBuilderError::NeighborMatrixCellPositionOverflow { position: pos })?;
    let two_bit_pos = pos / 4;
    let storage_len = res.len();
    let byte =
        res.get_mut(two_bit_pos)
            .ok_or(UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: pos,
                byte_index: two_bit_pos,
                storage_len,
            })?;

    // RDKit❗✔️:   unsigned int shift = 2 * (pos % 4);
    let shift = 2 * (source_pos % 4);
    // RDKit❗✔️:   std::uint8_t twoBitMask = 3 << shift;
    let two_bit_mask = 3_u8 << shift;
    // RDKit❗✔️:   res[twoBitPos] = ((res[twoBitPos] & (~twoBitMask)) | (value << shift));
    *byte = (*byte & !two_bit_mask) | (value << shift);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::setTwoBitCell
    // Complexity: O(1) source-width and slice-boundary checks followed by one
    // indexed byte read/modify/write; no allocation or scan is added.
    Ok(())
}

pub(super) fn get_two_bit_cell(res: &[u8], pos: usize) -> Result<u8, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::getTwoBitCell (Builder.cpp:70-78)
    // RDKit❗✔️: std::uint8_t getTwoBitCell(boost::shared_array<std::uint8_t> &res,
    // RDKit❗✔️:                            unsigned int pos) {
    // RDKit❗✔️:   unsigned int twoBitPos = pos / 4;
    let source_pos = u32::try_from(pos)
        .map_err(|_| UffBuilderError::NeighborMatrixCellPositionOverflow { position: pos })?;
    let two_bit_pos = pos / 4;
    let storage_len = res.len();
    let byte = *res
        .get(two_bit_pos)
        .ok_or(UffBuilderError::NeighborMatrixStorageOutOfRange {
            position: pos,
            byte_index: two_bit_pos,
            storage_len,
        })?;

    // RDKit❗✔️:   unsigned int shift = 2 * (pos % 4);
    let shift = 2 * (source_pos % 4);
    // RDKit❗✔️:   std::uint8_t twoBitMask = 3 << shift;
    let two_bit_mask = 3_u8 << shift;
    // RDKit❗✔️:   return ((res[twoBitPos] & twoBitMask) >> shift);
    let relation = (byte & two_bit_mask) >> shift;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::getTwoBitCell
    // Complexity: O(1) source-width and slice-boundary checks followed by one
    // indexed byte read; no allocation, scan or mutation is added.
    Ok(relation)
}

pub(super) fn build_neighbor_matrix(topology: &TopologyBlock) -> Result<Vec<u8>, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::Tools::buildNeighborMatrix (Builder.cpp:91-130)
    // RDKit❗✔️: boost::shared_array<std::uint8_t> buildNeighborMatrix(const ROMol &mol) {
    // RDKit❗✔️: enum { RELATION_1_2 = 0, RELATION_1_3 = 1, RELATION_1_4 = 2, RELATION_1_X = 3 };
    // RDKit❗✔️:   const std::uint8_t RELATION_1_X_INIT = RELATION_1_X | (RELATION_1_X << 2) |
    // RDKit❗✔️:                                          (RELATION_1_X << 4) |
    // RDKit❗✔️:                                          (RELATION_1_X << 6);
    const RELATION_1_2: u8 = 0;
    const RELATION_1_3: u8 = 1;
    const RELATION_1_X: u8 = 3;
    let relation_1_x_init =
        RELATION_1_X | (RELATION_1_X << 2) | (RELATION_1_X << 4) | (RELATION_1_X << 6);

    // RDKit❗✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = topology.atoms.len();
    if n_atoms != 0 {
        // Validate source-width allocation and cell-position bounds before
        // allocating; buildNeighborMatrix itself assumes a valid ROMol.
        let _first_cell = two_bit_cell_pos(n_atoms, 0, 0)?;
    }

    // RDKit❗✔️:   unsigned nTwoBitCells = (nAtoms * (nAtoms + 1) - 1) / 8 + 1;
    let source_n_atoms = u32::try_from(n_atoms)
        .map_err(|_| UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let source_next_n_atoms = source_n_atoms
        .checked_add(1)
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    let source_triangular_numerator = source_n_atoms
        .checked_mul(source_next_n_atoms)
        .ok_or(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;
    // The pinned unsigned subtraction wraps for an empty ROMol, yielding
    // 536,870,912 initialized bytes. Keep that edge explicit; positive counts
    // were validated through the checked source-index helper above.
    let source_byte_numerator = source_triangular_numerator.wrapping_sub(1);
    let source_n_two_bit_cells = source_byte_numerator / 8 + 1;
    let n_two_bit_cells = usize::try_from(source_n_two_bit_cells)
        .map_err(|_| UffBuilderError::NeighborMatrixIndexOverflow { n_atoms })?;

    // RDKit❗✔️:   boost::shared_array<std::uint8_t> res(new std::uint8_t[nTwoBitCells]);
    // RDKit❗✔️:   std::memset(res.get(), RELATION_1_X_INIT, nTwoBitCells);
    let mut res = Vec::new();
    res.try_reserve_exact(n_two_bit_cells).map_err(|_| {
        UffBuilderError::NeighborMatrixAllocationFailed {
            byte_count: n_two_bit_cells,
        }
    })?;
    res.resize(n_two_bit_cells, relation_1_x_init);

    // RDKit❗✔️:   for (ROMol::ConstBondIterator bondi = mol.beginBonds();
    // RDKit❗✔️:        bondi != mol.endBonds(); ++bondi) {
    for (bond_i_index, bond_i) in topology.bonds.iter().enumerate() {
        // RDKit❗✔️:     setTwoBitCell(res,
        // RDKit❗✔️:                   twoBitCellPos(nAtoms, (*bondi)->getBeginAtomIdx(),
        // RDKit❗✔️:                                 (*bondi)->getEndAtomIdx()),
        // RDKit❗✔️:                   RELATION_1_2);
        set_two_bit_cell(
            &mut res,
            two_bit_cell_pos(n_atoms, bond_i.begin().index(), bond_i.end().index())?,
            RELATION_1_2,
        )?;
        // RDKit❗✔️:     unsigned int bondiBeginAtomIdx = (*bondi)->getBeginAtomIdx();
        let bond_i_begin_atom_idx = bond_i.begin().index();
        // RDKit❗✔️:     unsigned int bondiEndAtomIdx = (*bondi)->getEndAtomIdx();
        let bond_i_end_atom_idx = bond_i.end().index();
        // RDKit❗✔️:     for (ROMol::ConstBondIterator bondj = bondi; ++bondj != mol.endBonds();) {
        for bond_j in topology.bonds.iter().skip(bond_i_index + 1) {
            // RDKit❗✔️:       int idx1 = -1;
            // RDKit❗✔️:       int idx3 = -1;
            let mut idx1 = None;
            let mut idx3 = None;
            // RDKit❗✔️:       unsigned int bondjBeginAtomIdx = (*bondj)->getBeginAtomIdx();
            let bond_j_begin_atom_idx = bond_j.begin().index();
            // RDKit❗✔️:       unsigned int bondjEndAtomIdx = (*bondj)->getEndAtomIdx();
            let bond_j_end_atom_idx = bond_j.end().index();
            // RDKit❗✔️:       if (bondiBeginAtomIdx == bondjBeginAtomIdx) {
            if bond_i_begin_atom_idx == bond_j_begin_atom_idx {
                // RDKit❗✔️:         idx1 = bondiEndAtomIdx;
                idx1 = Some(bond_i_end_atom_idx);
                // RDKit❗✔️:         idx3 = bondjEndAtomIdx;
                idx3 = Some(bond_j_end_atom_idx);
                // RDKit❗✔️:       } else if (bondiBeginAtomIdx == bondjEndAtomIdx) {
            } else if bond_i_begin_atom_idx == bond_j_end_atom_idx {
                // RDKit❗✔️:         idx1 = bondiEndAtomIdx;
                idx1 = Some(bond_i_end_atom_idx);
                // RDKit❗✔️:         idx3 = bondjBeginAtomIdx;
                idx3 = Some(bond_j_begin_atom_idx);
                // RDKit❗✔️:       } else if (bondiEndAtomIdx == bondjBeginAtomIdx) {
            } else if bond_i_end_atom_idx == bond_j_begin_atom_idx {
                // RDKit❗✔️:         idx1 = bondiBeginAtomIdx;
                idx1 = Some(bond_i_begin_atom_idx);
                // RDKit❗✔️:         idx3 = bondjEndAtomIdx;
                idx3 = Some(bond_j_end_atom_idx);
                // RDKit❗✔️:       } else if (bondiEndAtomIdx == bondjEndAtomIdx) {
            } else if bond_i_end_atom_idx == bond_j_end_atom_idx {
                // RDKit❗✔️:         idx1 = bondiBeginAtomIdx;
                idx1 = Some(bond_i_begin_atom_idx);
                // RDKit❗✔️:         idx3 = bondjBeginAtomIdx;
                idx3 = Some(bond_j_begin_atom_idx);
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:       if (idx1 > -1) {
            if let (Some(idx1), Some(idx3)) = (idx1, idx3) {
                // RDKit❗✔️:         setTwoBitCell(res, twoBitCellPos(nAtoms, idx1, idx3), RELATION_1_3);
                set_two_bit_cell(
                    &mut res,
                    two_bit_cell_pos(n_atoms, idx1, idx3)?,
                    RELATION_1_3,
                )?;
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:     }
        }
        // RDKit❗✔️:   }
    }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::Tools::buildNeighborMatrix
    // Complexity: O(A^2) packed initialization plus O(B^2) source bond-pair
    // comparisons; one output allocation and no allocation inside either loop.
    Ok(res)
}

#[cfg(test)]
#[derive(Debug)]
pub(super) struct PreparedTypingValence<'a> {
    pub(super) total_valences: Vec<i32>,
    topology: &'a TopologyBlock,
    assignment: &'a ValenceAssignment,
}

#[cfg(test)]
std::thread_local! {
    static TYPING_VALENCE_PROJECTION_VEC_CONSTRUCTIONS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
    static CONJUGATION_PROJECTION_VEC_CONSTRUCTIONS: std::cell::Cell<usize> =
        const { std::cell::Cell::new(0) };
}

#[cfg(test)]
pub(super) fn reset_typing_valence_projection_vec_constructions() {
    TYPING_VALENCE_PROJECTION_VEC_CONSTRUCTIONS.with(|count| count.set(0));
}

#[cfg(test)]
pub(super) fn typing_valence_projection_vec_constructions() -> usize {
    TYPING_VALENCE_PROJECTION_VEC_CONSTRUCTIONS.with(std::cell::Cell::get)
}

#[cfg(test)]
pub(super) fn reset_conjugation_projection_vec_constructions() {
    CONJUGATION_PROJECTION_VEC_CONSTRUCTIONS.with(|count| count.set(0));
}

#[cfg(test)]
pub(super) fn conjugation_projection_vec_constructions() -> usize {
    CONJUGATION_PROJECTION_VEC_CONSTRUCTIONS.with(std::cell::Cell::get)
}

#[cfg(test)]
impl PreparedTypingValence<'_> {
    pub(super) fn implicit_hydrogens(&self) -> impl ExactSizeIterator<Item = i32> + '_ {
        // BEGIN RDKIT CPP FUNCTION Atom::getNumImplicitHs (Atom.cpp:297-305)
        // RDKit❗✔️: unsigned int Atom::getNumImplicitHs() const {
        // RDKit❗✔️:   if (df_noImplicit) {
        // RDKit❗✔️:     return 0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   PRECONDITION(d_implicitValence > -1,
        // RDKit❗✔️:                "getNumImplicitHs() called without preceding call to "
        // RDKit❗✔️:                "calcImplicitValence()");
        // RDKit❗✔️:   return getValence(ValenceType::IMPLICIT);
        // RDKit❗✔️: }
        self.assignment
            .implicit_hydrogens
            .iter()
            .enumerate()
            .map(|(row, value)| {
                if self.topology.atoms[row].no_implicit() {
                    0
                } else {
                    *value
                }
            })
        // END RDKIT CPP FUNCTION Atom::getNumImplicitHs
    }
}

#[cfg(test)]
pub(super) fn prepare_typing_valence<'a>(
    topology: &'a TopologyBlock,
    assignment: &'a ValenceAssignment,
) -> Result<PreparedTypingValence<'a>, UffBuilderError> {
    let atom_count = topology.atoms.len();
    validate_typing_valence_cache(topology, assignment)?;

    #[cfg(test)]
    TYPING_VALENCE_PROJECTION_VEC_CONSTRUCTIONS.with(|count| count.set(count.get() + 1));
    let mut total_valences = Vec::with_capacity(atom_count);
    for (row, atom) in topology.atoms.iter().enumerate() {
        let explicit = assignment.explicit_valence[row];
        let implicit = if atom.no_implicit() {
            0
        } else {
            assignment.implicit_hydrogens[row]
        };

        // Valid signed-byte cache components sum to at most 254.
        total_valences.push(explicit + implicit);
    }

    // RDKit❗❌: unsigned int Atom::getTotalValence() const {
    // RDKit❗❌:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit❗❌: }
    // Complexity: validation and projection are two ordered O(V) passes plus
    // one output allocation. This test-only helper materializes the retained
    // regression projection; production reads cached values through borrowed
    // rows.
    Ok(PreparedTypingValence {
        total_valences,
        topology,
        assignment,
    })
}

#[cfg(test)]
pub(super) fn prepare_conjugated_presence(
    topology: &TopologyBlock,
) -> Result<Vec<bool>, UffBuilderError> {
    topology
        .validate()
        .map_err(UffBuilderError::TopologyValidation)?;

    #[cfg(test)]
    CONJUGATION_PROJECTION_VEC_CONSTRUCTIONS.with(|count| count.set(count.get() + 1));
    let mut atom_has_conjugated_bond = Vec::with_capacity(topology.atoms.len());
    for atom_index in 0..topology.atoms.len() {
        let has_conjugated_bond =
            atom_has_conjugated_bond_from_validated_topology(topology, atom_index);
        atom_has_conjugated_bond.push(has_conjugated_bond);
    }

    // Complexity: model validation rebuilds expected adjacency once, then
    // this ordered projection visits V rows and at most 2E incident entries.
    // It allocates V result values plus the validator's temporary adjacency.
    // This test-only helper materializes every row for retained regression
    // cases; production evaluates the same predicate at its source branch.
    // The production topology validation and temporary adjacency remain costs.
    Ok(atom_has_conjugated_bond)
}

pub(super) fn atom_has_conjugated_bond_from_validated_topology(
    topology: &TopologyBlock,
    atom_index: usize,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION MolOps::atomHasConjugatedBond (ConjugHybrid.cpp:116-127)
    // RDKit✔️✔️: bool atomHasConjugatedBond(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   for (const auto bnd : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (bnd->getIsConjugated()) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // The caller validated adjacency once at the borrowed cached-state
    // boundary; this shared leaf retains the source incident-bond scan.
    topology
        .adjacency
        .neighbors_of(atom_index)
        .iter()
        .any(|neighbor| topology.bonds[neighbor.bond.index()].is_conjugated())
    // END RDKIT CPP FUNCTION MolOps::atomHasConjugatedBond
}

pub(super) fn validate_typing_valence_cache(
    topology: &TopologyBlock,
    assignment: &ValenceAssignment,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getValence/getTotalValence (Atom.cpp:316-336)
    // RDKit❗❌: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit❗❌:   if (!dp_mol) {
    // RDKit❗❌:     return 0;
    // RDKit❗❌:   }
    // RDKit❗❌:   PRECONDITION(
    // RDKit❗❌:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit❗❌:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit❗❌:   PRECONDITION(
    // RDKit❗❌:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit❗❌:        d_implicitValence > -1),
    // RDKit❗❌:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit❗❌:   if (which == ValenceType::EXPLICIT) {
    // RDKit❗❌:     return d_explicitValence;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // RDKit❗❌: unsigned int Atom::getTotalValence() const {
    // RDKit❗❌:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit❗❌: }
    let atom_count = topology.atoms.len();
    if assignment.explicit_valence.len() != atom_count {
        return Err(UffBuilderError::ValenceAssignmentLengthMismatch {
            field: PreparedValenceField::Explicit,
            expected: atom_count,
            actual: assignment.explicit_valence.len(),
        });
    }
    if assignment.implicit_hydrogens.len() != atom_count {
        return Err(UffBuilderError::ValenceAssignmentLengthMismatch {
            field: PreparedValenceField::ImplicitHydrogen,
            expected: atom_count,
            actual: assignment.implicit_hydrogens.len(),
        });
    }

    // Source Atom stores computed cache components in signed int8_t fields,
    // initialized to -1. Keep explicit-before-implicit and noImplicit order.
    for (row, atom) in topology.atoms.iter().enumerate() {
        validate_source_cache_value(
            atom.id(),
            PreparedValenceField::Explicit,
            assignment.explicit_valence[row],
        )?;
        if !atom.no_implicit() {
            validate_source_cache_value(
                atom.id(),
                PreparedValenceField::ImplicitHydrogen,
                assignment.implicit_hydrogens[row],
            )?;
        }
    }

    // Complexity: one ordered O(V) validation pass with no allocation.
    Ok(())
    // END RDKIT CPP FUNCTION Atom::getValence/getTotalValence
}

fn validate_source_cache_value(
    atom_id: AtomId,
    field: PreparedValenceField,
    value: i32,
) -> Result<(), UffBuilderError> {
    if value < 0 {
        return Err(UffBuilderError::SourceValencePrecondition {
            atom_id,
            field,
            value,
        });
    }
    if value > i32::from(i8::MAX) {
        return Err(UffBuilderError::SourceValenceOutOfRange {
            atom_id,
            field,
            value,
        });
    }
    Ok(())
}

pub(super) fn needs_hydrogens_warning(
    topology: &TopologyBlock,
    implicit_hydrogens: &[i32],
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<bool, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::needsHs (AddHs.cpp:1340-1348)
    // RDKit❗✔️: bool needsHs(const ROMol &mol) {
    // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:     bool includeNeighbors = false;
    // RDKit❗✔️:     if (atom->getTotalNumHs(includeNeighbors)) {
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::needsHs

    // BEGIN RDKIT CPP FUNCTION Atom::getNumExplicitHs (Atom.h:221-224)
    // RDKit❗✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // END RDKIT CPP FUNCTION Atom::getNumExplicitHs

    // BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs (Atom.cpp:286-294)
    // RDKit❗✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit❗✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit❗✔️:   if (includeNeighbors && dp_mol) {
    // RDKit❗✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit❗✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit❗✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit❗✔️:     });
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalNumHs

    // BEGIN RDKIT CPP FUNCTION Atom::getNumImplicitHs (Atom.cpp:297-305)
    // RDKit❗✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit❗✔️:   if (df_noImplicit) {
    // RDKit❗✔️:     return 0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit❗✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit❗✔️:                "calcImplicitValence()");
    // RDKit❗✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Atom::getNumImplicitHs

    // BEGIN RDKIT CPP FUNCTION constructForceField missing-explicit-H warning (Builder.cpp:678-684)
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗✔️:
    // RDKit❗✔️:   if (MolOps::needsHs(mol)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION constructForceField missing-explicit-H warning

    let atom_count = topology.atoms.len();
    if implicit_hydrogens.len() != atom_count {
        return Err(UffBuilderError::ValenceAssignmentLengthMismatch {
            field: PreparedValenceField::ImplicitHydrogen,
            expected: atom_count,
            actual: implicit_hydrogens.len(),
        });
    }

    for (row, atom) in topology.atoms.iter().enumerate() {
        let explicit_hydrogens = i32::from(atom.explicit_hydrogens());
        let implicit_hydrogens = if atom.no_implicit() {
            0
        } else {
            let value = implicit_hydrogens[row];
            validate_source_cache_value(atom.id(), PreparedValenceField::ImplicitHydrogen, value)?;
            value
        };

        if explicit_hydrogens + implicit_hydrogens > 0 {
            diagnostics.push(UffTypingDiagnostic {
                atom_id: None,
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: NEEDS_EXPLICIT_HYDROGENS_WARNING_MESSAGE,
            });
            return Ok(true);
        }
    }

    // Complexity: the borrowed scan visits rows in source order and exits at
    // the first positive total; it allocates nothing unless the one warning is
    // appended to the existing ordered diagnostic sink.
    Ok(false)
}

fn validate_force_field_preamble(
    topology: &TopologyBlock,
    params: &UffParamsByAtom<'_>,
    implicit_hydrogens: &[i32],
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::constructForceField(params) preconditions
    // (Builder.cpp:678-684)
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗✔️:
    // RDKit❗✔️:   if (MolOps::needsHs(mol)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION UFF::constructForceField(params) preconditions
    if topology.atoms.len() != params.len() {
        return Err(UffBuilderError::ParamsLengthMismatch {
            atoms: topology.atoms.len(),
            params: params.len(),
        });
    }
    needs_hydrogens_warning(topology, implicit_hydrogens, diagnostics)?;
    Ok(())
}

fn append_selected_atom_positions<'a>(
    field: &mut ForceField<'a>,
    selected: &'a mut Conformer3D,
    atom_count: usize,
) -> Result<(), UffBuilderError> {
    // BEGIN RDKIT CPP HELPER Conformer::getAtomPos(unsigned int) (Conformer.cpp:39-45)
    // RDKit❗✔️: RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) {
    // RDKit❗✔️:   PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    // RDKit❗✔️:   URANGE_CHECK(atomId, d_positions.size());
    // RDKit❗✔️:   return d_positions[atomId];
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER Conformer::getAtomPos(unsigned int)
    // The source checks the attached molecule's complete atom-row count on
    // each requested row. The zero-atom loop makes no getAtomPos call, so it
    // does not validate an otherwise unused conformer row count.
    let coordinate_count = selected.coordinates().len();
    if atom_count != 0 && coordinate_count != atom_count {
        return Err(UffBuilderError::SelectedConformerCoordinateCountMismatch {
            conformer_id: selected.id(),
            atoms: atom_count,
            coordinates: coordinate_count,
        });
    }
    field.positions_mut().extend(
        selected
            .coordinates_mut()
            .iter_mut()
            .take(atom_count)
            .map(|position| &mut position[..]),
    );
    Ok(())
}

fn prepare_force_field_preamble<'a>(
    topology: &TopologyBlock,
    coordinates: &'a mut CoordinateBlock,
    selected_3d_conformer_id: usize,
    params: &UffParamsByAtom<'_>,
    implicit_hydrogens: &[i32],
    diagnostics: &mut Vec<UffTypingDiagnostic>,
) -> Result<ForceField<'a>, UffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION UFF::constructForceField(params) preamble
    // RDKit❗✔️: ForceFields::ForceField *constructForceField(ROMol &mol,
    // RDKit❗✔️:                                              const AtomicParamVect &params,
    // RDKit❗✔️:                                              double vdwThresh, int confId,
    // RDKit❗✔️:                                              bool ignoreInterfragInteractions) {
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗✔️:
    // RDKit❗✔️:   if (MolOps::needsHs(mol)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> res(new ForceFields::ForceField());
    // RDKit❗✔️:
    // RDKit❗✔️:   // add the atomic positions:
    // RDKit❗✔️:   Conformer &conf = mol.getConformer(confId);
    // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❗✔️:     res->positions().push_back(&conf.getAtomPos(i));
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION UFF::constructForceField(params) preamble
    // Reuse the source-ordered precondition/warning owner. It checks parameter
    // count before appending any needsHs diagnostic.
    validate_force_field_preamble(topology, params, implicit_hydrogens, diagnostics)?;

    // BEGIN RDKIT CPP HELPER ROMol::getConformer(int) (ROMol.cpp:630-655)
    // RDKit❗✔️: const Conformer &ROMol::getConformer(int id) const {
    // RDKit❗✔️:   // make sure we have more than one conformation
    // RDKit❗✔️:   if (d_confs.size() == 0) {
    // RDKit❗✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (id < 0) {
    // RDKit❗✔️:     return *(d_confs.front());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto cid = (unsigned int)id;
    // RDKit❗✔️:   for (auto conf : d_confs) {
    // RDKit❗✔️:     if (conf->getId() == cid) {
    // RDKit❗✔️:       return *conf;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // we did not find a conformation with the specified ID
    // RDKit❗✔️:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit❗✔️:   mesg += id;
    // RDKit❗✔️:   throw ConformerException(mesg);
    // RDKit❗✔️: }
    // RDKit❗✔️: Conformer &ROMol::getConformer(int id) {
    // RDKit❗✔️:   return const_cast<Conformer &>(
    // RDKit❗✔️:       static_cast<const ROMol *>(this)->getConformer(id));
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER ROMol::getConformer(int)
    // The detached forcefield entry receives an already explicit 3D ID, so
    // the source's negative-ID first-conformer fallback is not part of this
    // helper contract. For a represented nonnegative ID, retain the source
    // linear ID search, not vector-position indexing. The separate 3D model
    // collection is the dimension boundary; the stored is_3d metadata flag
    // does not alter RDKit's getConformer lookup.
    let mut field = ForceField::new(3);
    let Some(conformer) = coordinates
        .conformers_3d
        .iter_mut()
        .find(|conformer| conformer.id() == selected_3d_conformer_id)
    else {
        return Err(UffBuilderError::SelectedThreeDimensionalConformerNotFound {
            conformer_id: selected_3d_conformer_id,
        });
    };
    append_selected_atom_positions(&mut field, conformer, topology.atoms.len())?;

    // Behavior marker — RDKit❗✔️: preserve parameter failure before the
    // existing needsHs warning, warning before ID lookup, selected-3D-ID
    // lookup, and atom-row position order. Rust models explicit nonnegative
    // 3D selection; it does not use a dimension-mixing or negative-ID default.
    // Complexity marker — RDKit✔️✔️: one warning scan, one linear ID scan,
    // and one borrowed row insertion pass; no coordinate values or block are
    // copied, and the kernel remains uninitialized as in the source preamble.
    Ok(field)
}

#[allow(clippy::too_many_arguments)]
fn construct_force_field_with_params<'a>(
    topology: &TopologyBlock,
    coordinates: &'a mut CoordinateBlock,
    selected_3d_conformer_id: usize,
    params: &UffParamsByAtom<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    torsion_bond_smarts: &str,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
) -> Result<ForceField<'a>, ForceFieldConstructionError> {
    // BEGIN RDKIT CPP FUNCTION UFF::constructForceField(params)
    // (Builder.cpp:674-702)
    // RDKit❗✔️: ForceFields::ForceField *constructForceField(ROMol &mol,
    // RDKit❗✔️:                                              const AtomicParamVect &params,
    // RDKit❗✔️:                                              double vdwThresh, int confId,
    // RDKit❗✔️:                                              bool ignoreInterfragInteractions) {
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗✔️:
    // RDKit❗✔️:   if (MolOps::needsHs(mol)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> res(new ForceFields::ForceField());
    // RDKit❗✔️:
    // RDKit❗✔️:   // add the atomic positions:
    // RDKit❗✔️:   Conformer &conf = mol.getConformer(confId);
    // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❗✔️:     res->positions().push_back(&conf.getAtomPos(i));
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   Tools::addBonds(mol, params, res.get());
    // RDKit❗✔️:   Tools::addAngles(mol, params, res.get());
    // RDKit❗✔️:   Tools::addAngleSpecialCases(mol, confId, params, res.get());
    // RDKit❗✔️:   boost::shared_array<std::uint8_t> neighborMat =
    // RDKit❗✔️:       Tools::buildNeighborMatrix(mol);
    // RDKit❗✔️:   Tools::addNonbonded(mol, confId, params, res.get(), neighborMat,
    // RDKit❗✔️:                       vdwThresh, ignoreInterfragInteractions);
    // RDKit❗✔️:   Tools::addTorsions(mol, params, res.get());
    // RDKit❗✔️:   Tools::addInversions(mol, params, res.get());
    // RDKit❗✔️:
    // RDKit❗✔️:   return res.release();
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::constructForceField(params)

    validate_force_field_preamble(topology, params, &valence.implicit_hydrogens, diagnostics)?;

    // The source's default ForceField dimension is three. Construct exactly
    // one field before its conformer lookup and keep it uninitialized here.
    let field = ForceField::new(3);
    let CoordinateBlock {
        conformers_2d,
        conformers_3d,
        source_coordinate_dim,
    } = coordinates;
    let source_coordinate_dim = *source_coordinate_dim;
    let conformers_2d = conformers_2d.as_slice();
    // ROMol::getConformer searches by ID in insertion order; the explicit
    // detached constructor input is already a selected 3D ID.
    let Some(selected_index) = conformers_3d
        .iter()
        .position(|conformer| conformer.id() == selected_3d_conformer_id)
    else {
        return Err(UffBuilderError::SelectedThreeDimensionalConformerNotFound {
            conformer_id: selected_3d_conformer_id,
        }
        .into());
    };
    let (conformers_3d_before, selected_and_after) = conformers_3d.split_at_mut(selected_index);
    let (selected_slice, conformers_3d_after) = selected_and_after.split_at_mut(1);
    let selected = &mut selected_slice[0];
    let selected_conformer_is_3d = selected.is_3d();

    // F19 needs the selected conformer's properties alongside mutable
    // selected coordinate rows. Model accessors intentionally keep those
    // fields separate; copy only this small metadata map when that branch is
    // requested, while every other conformer and all coordinate rows remain
    // borrowed. This bounded O(P) metadata clone does not copy coordinates.
    let selected_conformer_props = if ignore_interfragment_interactions {
        selected.props().clone()
    } else {
        BTreeMap::new()
    };
    construct_force_field_with_selected_params(
        field,
        topology,
        selected,
        &UffConformerContext {
            two_d: conformers_2d,
            before: conformers_3d_before,
            selected_id: selected_3d_conformer_id,
            selected_is_3d: selected_conformer_is_3d,
            selected_props: &selected_conformer_props,
            after: conformers_3d_after,
            source_dimension: source_coordinate_dim,
        },
        params,
        rings,
        valence,
        molecule_properties,
        torsion_bond_smarts,
        vdw_threshold,
        ignore_interfragment_interactions,
    )
}

pub(super) struct UffConformerContext<'a> {
    pub(super) two_d: &'a [Conformer2D],
    pub(super) before: &'a [Conformer3D],
    pub(super) selected_id: usize,
    pub(super) selected_is_3d: bool,
    pub(super) selected_props: &'a BTreeMap<String, String>,
    pub(super) after: &'a [Conformer3D],
    pub(super) source_dimension: Option<CoordinateDimension>,
}

#[allow(clippy::too_many_arguments)]
fn construct_force_field_with_selected_params<'a>(
    mut field: ForceField<'a>,
    topology: &TopologyBlock,
    selected: &'a mut Conformer3D,
    context: &UffConformerContext<'_>,
    params: &UffParamsByAtom<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    torsion_bond_smarts: &str,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
) -> Result<ForceField<'a>, ForceFieldConstructionError> {
    // BEGIN RDKIT CPP FUNCTION UFF::constructForceField(params)
    // (Builder.cpp:674-702)
    // RDKit❗✔️: ForceFields::ForceField *constructForceField(ROMol &mol,
    // RDKit❗✔️:                                              const AtomicParamVect &params,
    // RDKit❗✔️:                                              double vdwThresh, int confId,
    // RDKit❗✔️:                                              bool ignoreInterfragInteractions) {
    // RDKit❗✔️:   PRECONDITION(mol.getNumAtoms() == params.size(), "bad parameters");
    // RDKit❗✔️:
    // RDKit❗✔️:   if (MolOps::needsHs(mol)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> res(new ForceFields::ForceField());
    // RDKit❗✔️:
    // RDKit❗✔️:   // add the atomic positions:
    // RDKit❗✔️:   Conformer &conf = mol.getConformer(confId);
    // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❗✔️:     res->positions().push_back(&conf.getAtomPos(i));
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   Tools::addBonds(mol, params, res.get());
    // RDKit❗✔️:   Tools::addAngles(mol, params, res.get());
    // RDKit❗✔️:   Tools::addAngleSpecialCases(mol, confId, params, res.get());
    // RDKit❗✔️:   boost::shared_array<std::uint8_t> neighborMat =
    // RDKit❗✔️:       Tools::buildNeighborMatrix(mol);
    // RDKit❗✔️:   Tools::addNonbonded(mol, confId, params, res.get(), neighborMat,
    // RDKit❗✔️:                       vdwThresh, ignoreInterfragInteractions);
    // RDKit❗✔️:   Tools::addTorsions(mol, params, res.get());
    // RDKit❗✔️:   Tools::addInversions(mol, params, res.get());
    // RDKit❗✔️:
    // RDKit❗✔️:   return res.release();
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::constructForceField(params)

    // Behavior: the same source contribution owners and order serve both
    // writable optimization and read-only evaluation; all context remains
    // borrowed. Complexity: O(A) position handles, no unrelated coordinate,
    // topology or runtime-state cloning.
    append_selected_atom_positions(&mut field, selected, topology.atoms.len())?;

    // Keep the exact source contribution and failure order. Each stage is the
    // existing private owner; errors drop this local field and its term prefix.
    add_bonds(topology, params, &mut field)?;
    add_angles(topology, params, rings, &mut field)?;
    add_angle_special_cases(topology, params, &mut field)?;
    let neighbor_matrix = build_neighbor_matrix(topology)?;
    add_nonbonded(
        topology,
        context.two_d,
        context.before,
        context.selected_id,
        context.selected_is_3d,
        context.selected_props,
        context.after,
        context.source_dimension,
        molecule_properties,
        params,
        &mut field,
        &neighbor_matrix,
        vdw_threshold,
        ignore_interfragment_interactions,
    )?;
    add_torsions(
        topology,
        rings,
        valence,
        params,
        torsion_bond_smarts,
        &mut field,
    )?;
    add_inversions(topology, params, &mut field)?;

    Ok(field)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn construct_force_field_with_automatic_typing<'a>(
    topology: &TopologyBlock,
    coordinates: &'a mut CoordinateBlock,
    selected_3d_conformer_id: usize,
    typing_state: UffAtomStateRef<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    torsion_bond_smarts: &str,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
) -> Result<ForceField<'a>, AutomaticForceFieldConstructionError> {
    // BEGIN RDKIT CPP FUNCTION UFF::constructForceField(automatic typing)
    // (Builder.cpp:712-719)
    // RDKit❗✔️: ForceFields::ForceField *constructForceField(ROMol &mol, double vdwThresh,
    // RDKit❗✔️:                                              int confId,
    // RDKit❗✔️:                                              bool ignoreInterfragInteractions) {
    // RDKit❗✔️:   bool foundAll;
    // RDKit❗✔️:   AtomicParamVect params;
    // RDKit❗✔️:   boost::tie(params, foundAll) = getAtomTypes(mol);
    // RDKit❗✔️:   return constructForceField(mol, params, vdwThresh, confId,
    // RDKit❗✔️:                              ignoreInterfragInteractions);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::constructForceField(automatic typing)

    with_automatic_atom_types(
        topology,
        typing_state,
        diagnostics,
        |params, diagnostics| {
            construct_force_field_with_params(
                topology,
                coordinates,
                selected_3d_conformer_id,
                params,
                rings,
                valence,
                molecule_properties,
                diagnostics,
                torsion_bond_smarts,
                vdw_threshold,
                ignore_interfragment_interactions,
            )
        },
    )
}

fn with_automatic_atom_types<'a>(
    topology: &TopologyBlock,
    typing_state: UffAtomStateRef<'_>,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    build: impl FnOnce(
        &UffParamsByAtom<'_>,
        &mut Vec<UffTypingDiagnostic>,
    ) -> Result<ForceField<'a>, ForceFieldConstructionError>,
) -> Result<ForceField<'a>, AutomaticForceFieldConstructionError> {
    // BEGIN RDKIT CPP HELPER UFF::getAtomTypes (AtomTyper.cpp:507-533)
    // RDKit❗✔️: std::pair<AtomicParamVect, bool> getAtomTypes(const ROMol &mol,
    // RDKit❗✔️:                                               const std::string &) {
    // RDKit❗✔️:   bool foundAll = true;
    // RDKit❗✔️:   auto params = ParamCollection::getParams();
    // RDKit❗✔️:
    // RDKit❗✔️:   AtomicParamVect paramVect;
    // RDKit❗✔️:   paramVect.resize(mol.getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❗✔️:     const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit❗✔️:
    // RDKit❗✔️:     // construct the atom key:
    // RDKit❗✔️:     std::string atomKey = Tools::getAtomLabel(atom);
    // RDKit❗✔️:
    // RDKit❗✔️:     // ok, we've got the atom key, now get the parameters:
    // RDKit❗✔️:     const AtomicParams *theParams = (*params)(atomKey);
    // RDKit❗✔️:     if (!theParams) {
    // RDKit❗✔️:       foundAll = false;
    // RDKit❗✔️:       BOOST_LOG(rdErrorLog) << "UFFTYPER: Unrecognized atom type: " << atomKey
    // RDKit❗✔️:                             << " (" << i << ")" << std::endl;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     paramVect[i] = theParams;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return std::make_pair(paramVect, foundAll);
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER UFF::getAtomTypes

    // This detached constructor consumes the already-prepared valence and
    // conjugation rows once, through the existing source-shaped atom typer.
    // It does not sanitize the topology or recompute either chemistry state.
    // RDKit ignores `foundAll` in the overload above and always delegates the
    // nullable vector, so missing rows remain nonfatal diagnostic slots here.
    let params = ParamCollection::get_params("")
        .map_err(AutomaticForceFieldConstructionError::ParameterTable)?;
    // Supplied rows remain a thin adapter, but their shape validation stays
    // after default-table acquisition, matching source getAtomTypes order.
    let typing_state = match typing_state {
        UffAtomStateRef::Cached { .. } => Ok(typing_state),
        UffAtomStateRef::SuppliedRows {
            total_valences,
            conjugated_presence,
        } => UffAtomStateRef::supplied_rows(topology, total_valences, conjugated_presence),
    }
    .map_err(AutomaticForceFieldConstructionError::Typing)?;
    let (params_by_atom, _found_all) =
        get_atom_types_from_state(topology, typing_state, params.as_ref(), diagnostics)
            .map_err(AutomaticForceFieldConstructionError::Typing)?;

    // `params` remains alive through delegation, while all generated kernel
    // terms own their source scalar parameters and retain only the selected
    // coordinate-row borrow. Typing diagnostics therefore precede the
    // delegated parameter preamble, needs-H warning, and construction errors.
    build(&params_by_atom, diagnostics).map_err(AutomaticForceFieldConstructionError::Construction)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn construct_force_field_with_automatic_typing_from_selected<'a>(
    topology: &TopologyBlock,
    selected: &'a mut Conformer3D,
    context: &UffConformerContext<'_>,
    typing_state: UffAtomStateRef<'_>,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
) -> Result<ForceField<'a>, AutomaticForceFieldConstructionError> {
    // BEGIN RDKIT CPP FUNCTION UFF::constructForceField(automatic typing)
    // (Builder.cpp:712-719)
    // RDKit❗✔️: ForceFields::ForceField *constructForceField(ROMol &mol, double vdwThresh,
    // RDKit❗✔️:                                              int confId,
    // RDKit❗✔️:                                              bool ignoreInterfragInteractions) {
    // RDKit❗✔️:   bool foundAll;
    // RDKit❗✔️:   AtomicParamVect params;
    // RDKit❗✔️:   boost::tie(params, foundAll) = getAtomTypes(mol);
    // RDKit❗✔️:   return constructForceField(mol, params, vdwThresh, confId,
    // RDKit❗✔️:                              ignoreInterfragInteractions);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION UFF::constructForceField(automatic typing)

    // Scalar parameters and all chemistry stages reuse their unique owners.
    // Only the selected conformer's cloned row is mutable; peers are borrowed.
    with_automatic_atom_types(
        topology,
        typing_state,
        diagnostics,
        |params, diagnostics| {
            validate_force_field_preamble(
                topology,
                params,
                &valence.implicit_hydrogens,
                diagnostics,
            )?;
            construct_force_field_with_selected_params(
                ForceField::new(3),
                topology,
                selected,
                context,
                params,
                rings,
                valence,
                molecule_properties,
                DEFAULT_TORSION_BOND_SMARTS,
                vdw_threshold,
                ignore_interfragment_interactions,
            )
        },
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn construct_force_field_with_automatic_typing_from_rows<'a>(
    topology: &TopologyBlock,
    coordinates: &'a mut CoordinateBlock,
    selected_3d_conformer_id: usize,
    total_valences: &[i32],
    atom_has_conjugated_bond: &[bool],
    rings: &RingInfo,
    valence: &ValenceAssignment,
    molecule_properties: &MoleculeProperties,
    diagnostics: &mut Vec<UffTypingDiagnostic>,
    torsion_bond_smarts: &str,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
) -> Result<ForceField<'a>, AutomaticForceFieldConstructionError> {
    construct_force_field_with_automatic_typing(
        topology,
        coordinates,
        selected_3d_conformer_id,
        UffAtomStateRef::SuppliedRows {
            total_valences,
            conjugated_presence: atom_has_conjugated_bond,
        },
        rings,
        valence,
        molecule_properties,
        diagnostics,
        torsion_bond_smarts,
        vdw_threshold,
        ignore_interfragment_interactions,
    )
}

#[cfg(test)]
mod tests {
    use super::super::atom_typer::{UffAtomStateRef, UffTypingInput, get_atom_types};
    use super::super::convenience::{
        OptimizationOutcome, OptimizationStageError, SerialConformer, SerialUffOptimizationError,
        SingleConformerOptimizationError, SingleConformerOptions, optimize_serial_uff,
        optimize_single_conformer,
    };
    use super::super::inversion::InversionIndexArgument;
    use super::construct_force_field_with_automatic_typing_from_rows as construct_force_field_with_automatic_typing;
    use super::*;
    use crate::kernel::{
        AngleIndexArgument, BondIndexArgument, Cf3dFragAcceptContributionIdentity,
        EvaluationContext, ForceFieldContribution, ForceFieldKernelError, cf3d_bld_b05_calc_energy,
        cf3d_bld_b05_calc_grad, cf3d_frag_accept_contribution_identities,
        cf3d_frag_f24_contribution_energies,
    };
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomQueryPredicate, AtomSpec, Bond, BondId, BondOrder,
        BondQueryPredicate, BondSpec, Conformer2D, Conformer3D, CoordinateDimension, Element,
        MoleculeProperties, QueryNode, SGroupAttachPoint, SGroupBondRole, SGroupBracket,
        SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData, SGroupDisplay,
        SdfPropertyList, SdfPropertyListTarget, StereoGroup, StereoGroupKind, SubstanceGroup,
        SubstanceGroupId, SubstanceGroupKind,
    };

    #[test]
    fn uff_error_e07_builder_causes_borrow_stored_children() {
        // This literal table checks trait dispatch for every existing variant.
        // B02/B15/B24/B26/B28 continue to exercise actual caller failures.
        let cases = [
            UffBuilderError::NeighborMatrixIndexOutOfRange {
                n_atoms: 19,
                i: 3,
                j: 17,
            },
            UffBuilderError::NeighborMatrixIndexOverflow { n_atoms: 23 },
            UffBuilderError::NeighborMatrixAllocationFailed { byte_count: 4_073 },
            UffBuilderError::NeighborMatrixCellPositionOverflow { position: 577 },
            UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: 59,
                byte_index: 7,
                storage_len: 53,
            },
            UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::ImplicitHydrogen,
                expected: 19,
                actual: 23,
            },
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(29),
                field: PreparedValenceField::Explicit,
                value: -11,
            },
            UffBuilderError::SourceValenceOutOfRange {
                atom_id: AtomId::new(31),
                field: PreparedValenceField::ImplicitHydrogen,
                value: 127,
            },
            UffBuilderError::ParamsLengthMismatch {
                atoms: 31,
                params: 29,
            },
            UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 47 },
            UffBuilderError::SelectedConformerCoordinateCountMismatch {
                conformer_id: 53,
                atoms: 59,
                coordinates: 61,
            },
            UffBuilderError::SourceBondIndexOverflow { atom_index: 67 },
            UffBuilderError::SourceTorsionIndexOverflow { atom_index: 71 },
            UffBuilderError::SourceAngleIndexOverflow { atom_index: 73 },
            UffBuilderError::SourceInversionIndexOverflow { atom_index: 79 },
            UffBuilderError::SourceDirectionVectorBelowTolerance {
                center_atom_index: 83,
                neighbor_atom_index: 89,
            },
            UffBuilderError::SourceTbpAtomPrecondition {
                center_atom_index: 97,
            },
            UffBuilderError::SourceTbpHybridizationPrecondition {
                center_atom_index: 101,
                actual: Hybridization::Sp2,
            },
            UffBuilderError::SourceTbpDegreePrecondition {
                center_atom_index: 103,
                actual_degree: 5,
            },
            UffBuilderError::SourceTbpAxialBondNotFound {
                center_atom_index: 107,
            },
            UffBuilderError::SourceTbpEquatorialBondNotFound {
                center_atom_index: 109,
                role: 3,
            },
            UffBuilderError::SourceTbpCenterParamsMissing {
                center_atom_index: 113,
            },
            UffBuilderError::TorsionBondQuery(TorsionBondQueryError::MissingMatchedBond {
                begin_atom_index: 17,
                end_atom_index: 29,
            }),
            UffBuilderError::InversionContribution(InversionContributionError::IndexOutOfRange {
                argument: InversionIndexArgument::Fourth,
                index: 31,
                upper_bound: 47,
            }),
            UffBuilderError::Valence(ValenceError::BadBondType {
                bond: Some(BondId::new(19)),
                order: BondOrder::DativeLeft,
            }),
            UffBuilderError::ForceFieldKernel(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Third,
                index: 13,
                upper_bound: 29,
            }),
            UffBuilderError::TopologyValidation(TopologyValidationError::AdjacencyMismatch),
        ];

        assert_eq!(cases.len(), 27);
        for error in cases {
            assert_eq!(error.to_string(), format!("{error:?}"));
            let downcast = (&error as &dyn std::error::Error)
                .downcast_ref::<UffBuilderError>()
                .expect("the outer error retains its concrete type");
            assert_eq!(downcast, &error);

            match &error {
                UffBuilderError::NeighborMatrixIndexOutOfRange { .. }
                | UffBuilderError::NeighborMatrixIndexOverflow { .. }
                | UffBuilderError::NeighborMatrixAllocationFailed { .. }
                | UffBuilderError::NeighborMatrixCellPositionOverflow { .. }
                | UffBuilderError::NeighborMatrixStorageOutOfRange { .. }
                | UffBuilderError::ValenceAssignmentLengthMismatch { .. }
                | UffBuilderError::SourceValencePrecondition { .. }
                | UffBuilderError::SourceValenceOutOfRange { .. }
                | UffBuilderError::ParamsLengthMismatch { .. }
                | UffBuilderError::SelectedThreeDimensionalConformerNotFound { .. }
                | UffBuilderError::SelectedConformerCoordinateCountMismatch { .. }
                | UffBuilderError::SourceBondIndexOverflow { .. }
                | UffBuilderError::SourceTorsionIndexOverflow { .. }
                | UffBuilderError::SourceAngleIndexOverflow { .. }
                | UffBuilderError::SourceInversionIndexOverflow { .. }
                | UffBuilderError::SourceDirectionVectorBelowTolerance { .. }
                | UffBuilderError::SourceTbpAtomPrecondition { .. }
                | UffBuilderError::SourceTbpHybridizationPrecondition { .. }
                | UffBuilderError::SourceTbpDegreePrecondition { .. }
                | UffBuilderError::SourceTbpAxialBondNotFound { .. }
                | UffBuilderError::SourceTbpEquatorialBondNotFound { .. }
                | UffBuilderError::SourceTbpCenterParamsMissing { .. } => {
                    assert!(std::error::Error::source(&error).is_none());
                }
                UffBuilderError::TorsionBondQuery(stored) => {
                    let exposed = std::error::Error::source(&error)
                        .expect("TorsionBondQuery exposes its stored child")
                        .downcast_ref::<TorsionBondQueryError>()
                        .expect("the torsion child keeps its concrete type");
                    assert!(std::ptr::eq(stored, exposed));
                }
                UffBuilderError::InversionContribution(stored) => {
                    let exposed = std::error::Error::source(&error)
                        .expect("InversionContribution exposes its stored child")
                        .downcast_ref::<InversionContributionError>()
                        .expect("the inversion child keeps its concrete type");
                    assert!(std::ptr::eq(stored, exposed));
                }
                UffBuilderError::Valence(stored) => {
                    let exposed = std::error::Error::source(&error)
                        .expect("Valence exposes its stored child")
                        .downcast_ref::<ValenceError>()
                        .expect("the valence child keeps its concrete type");
                    assert!(std::ptr::eq(stored, exposed));
                }
                UffBuilderError::ForceFieldKernel(stored) => {
                    let exposed = std::error::Error::source(&error)
                        .expect("ForceFieldKernel exposes its stored child")
                        .downcast_ref::<ForceFieldKernelError>()
                        .expect("the kernel child keeps its concrete type");
                    assert!(std::ptr::eq(stored, exposed));
                }
                UffBuilderError::TopologyValidation(stored) => {
                    let exposed = std::error::Error::source(&error)
                        .expect("TopologyValidation exposes its stored child")
                        .downcast_ref::<TopologyValidationError>()
                        .expect("the topology child keeps its concrete type");
                    assert!(std::ptr::eq(stored, exposed));
                }
            }
        }
    }

    #[test]
    fn uff_error_e06_mapping_and_assembly_causes_preserve_identity() {
        // The manually constructed mapping rows below check trait dispatch
        // only. The actual F19/F22 caller failures remain separate source
        // evidence for opaque fragment-copy and PairBuilder errors.
        let topology = f19_connected_overvalent_topology();
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let coordinates = cf3d_frag_f22_coordinates(&points);
        let properties = MoleculeProperties::default().with_name("E06-fragment-copy");
        let all_missing = vec![None; topology.atoms.len()];
        let mut rows = coordinates.conformers_3d[1].coordinates().to_vec();
        let actual_fragment_assembly = match cf3d_frag_f22_append(
            &topology,
            &coordinates,
            &properties,
            &all_missing,
            &mut rows,
            &[],
            2.0,
            true,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("the pinned F19 sanitize fixture must fail before pair skips")
            }
        };
        assert_eq!(
            actual_fragment_assembly.to_string(),
            format!("{actual_fragment_assembly:?}")
        );
        let stored_mapping = match &actual_fragment_assembly {
            NonbondedAssemblyError::FragmentMapping(mapping) => mapping,
            NonbondedAssemblyError::PairBuilder(error) => {
                panic!("expected the actual F19 mapping error, got {error:?}")
            }
        };
        let exposed_mapping = std::error::Error::source(&actual_fragment_assembly)
            .expect("FragmentMapping exposes its stored child")
            .downcast_ref::<NonbondedFragmentMappingError>()
            .expect("the assembly source retains the mapping error type");
        assert!(std::ptr::eq(stored_mapping, exposed_mapping));
        let stored_fragment_copy = match stored_mapping {
            NonbondedFragmentMappingError::FragmentCopy(error) => error,
            _ => panic!("the F19 fixture must retain its opaque copy error"),
        };
        assert_eq!(stored_fragment_copy.component_index(), Some(0));
        let exposed_fragment_copy = std::error::Error::source(stored_mapping)
            .expect("FragmentCopy exposes its stored core error")
            .downcast_ref::<MoleculeFragmentsError>()
            .expect("the mapping source retains MoleculeFragmentsError");
        assert!(std::ptr::eq(stored_fragment_copy, exposed_fragment_copy));
        assert!(
            std::error::Error::source(exposed_fragment_copy).is_some(),
            "MoleculeFragmentsError keeps its existing downstream failure source"
        );

        let mapping_cases = [
            NonbondedFragmentMappingError::CoordinateView(
                FragmentCoordinateViewError::SelectedConformerDuplicated { id: 71 },
            ),
            NonbondedFragmentMappingError::RequestedMappingMissing,
            NonbondedFragmentMappingError::AtomOutOfRange {
                component_index: 73,
                atom_index: 79,
                atom_count: 83,
            },
            NonbondedFragmentMappingError::DuplicateAtom {
                atom_index: 89,
                first_component: 97,
                second_component: 101,
            },
            NonbondedFragmentMappingError::MissingAtom { atom_index: 103 },
        ];
        assert_eq!(mapping_cases.len() + 1, 6);
        for mapping_error in mapping_cases {
            assert_eq!(mapping_error.to_string(), format!("{mapping_error:?}"));
            let downcast = (&mapping_error as &dyn std::error::Error)
                .downcast_ref::<NonbondedFragmentMappingError>()
                .expect("the mapping error retains its concrete type");
            assert!(std::ptr::eq(&mapping_error, downcast));
            match &mapping_error {
                NonbondedFragmentMappingError::FragmentCopy(_) => {
                    panic!("opaque copy errors are exercised through the F19 caller")
                }
                NonbondedFragmentMappingError::CoordinateView(stored) => {
                    let exposed = std::error::Error::source(&mapping_error)
                        .expect("CoordinateView exposes its stored child")
                        .downcast_ref::<FragmentCoordinateViewError>()
                        .expect("the coordinate child retains its concrete type");
                    assert!(std::ptr::eq(stored, exposed));
                    assert!(matches!(
                        exposed,
                        FragmentCoordinateViewError::SelectedConformerDuplicated { id: 71 }
                    ));
                }
                NonbondedFragmentMappingError::RequestedMappingMissing => {
                    assert!(std::error::Error::source(&mapping_error).is_none());
                }
                NonbondedFragmentMappingError::AtomOutOfRange {
                    component_index,
                    atom_index,
                    atom_count,
                } => {
                    assert_eq!((*component_index, *atom_index, *atom_count), (73, 79, 83));
                    assert!(std::error::Error::source(&mapping_error).is_none());
                }
                NonbondedFragmentMappingError::DuplicateAtom {
                    atom_index,
                    first_component,
                    second_component,
                } => {
                    assert_eq!(
                        (*atom_index, *first_component, *second_component),
                        (89, 97, 101)
                    );
                    assert!(std::error::Error::source(&mapping_error).is_none());
                }
                NonbondedFragmentMappingError::MissingAtom { atom_index } => {
                    assert_eq!(*atom_index, 103);
                    assert!(std::error::Error::source(&mapping_error).is_none());
                }
            }
        }

        let short_params = vec![None; topology.atoms.len() - 1];
        let mut rows = coordinates.conformers_3d[1].coordinates().to_vec();
        let actual_pair_assembly = match cf3d_frag_f22_append(
            &topology,
            &coordinates,
            &properties,
            &short_params,
            &mut rows,
            &[],
            2.0,
            true,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("the source parameter-count check must precede fragment copying")
            }
        };
        assert_eq!(
            actual_pair_assembly.to_string(),
            format!("{actual_pair_assembly:?}")
        );
        let stored_pair_builder = match &actual_pair_assembly {
            NonbondedAssemblyError::PairBuilder(error) => error,
            NonbondedAssemblyError::FragmentMapping(error) => {
                panic!("expected the actual PairBuilder error, got {error:?}")
            }
        };
        assert!(matches!(
            stored_pair_builder,
            UffBuilderError::ParamsLengthMismatch {
                atoms: 4,
                params: 3,
            }
        ));
        let exposed_pair_builder = std::error::Error::source(&actual_pair_assembly)
            .expect("PairBuilder exposes its stored child")
            .downcast_ref::<UffBuilderError>()
            .expect("the assembly source retains UffBuilderError");
        assert!(std::ptr::eq(stored_pair_builder, exposed_pair_builder));
        assert!(std::error::Error::source(exposed_pair_builder).is_none());
    }

    #[test]
    fn uff_error_e08_construction_causes_borrow_stored_children() {
        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let exposed = std::error::Error::source(parent)
                .expect("the source-bearing variant exposes its stored child")
                .downcast_ref::<Child>()
                .expect("the source keeps the child's concrete type");
            assert!(std::ptr::eq(stored, exposed));
        }

        // These hand-built values test trait dispatch only. Existing U06
        // prepared-state and fragment failures and U09 selected-ID failure
        // exercise real callers. The immutable default parameter table does
        // not provide a source-backed ParameterTable failure fixture.
        let parameter_table =
            AutomaticForceFieldConstructionError::ParameterTable(UffParamError::EmptyLine {
                line_number: 149,
            });
        assert_eq!(
            parameter_table.to_string(),
            "ParameterTable(EmptyLine { line_number: 149 })"
        );
        let parameter_error = match &parameter_table {
            AutomaticForceFieldConstructionError::ParameterTable(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&parameter_table, parameter_error);
        assert!(std::error::Error::source(parameter_error).is_none());

        let typing =
            AutomaticForceFieldConstructionError::Typing(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 151,
                actual: 157,
            });
        assert_eq!(
            typing.to_string(),
            "Typing(PreparedStateLength { input: TotalValence, expected: 151, actual: 157 })"
        );
        let typing_error = match &typing {
            AutomaticForceFieldConstructionError::Typing(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&typing, typing_error);
        assert!(std::error::Error::source(typing_error).is_none());

        let builder_construction = AutomaticForceFieldConstructionError::Construction(
            ForceFieldConstructionError::Builder(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::OptimizerBadDirection,
            )),
        );
        assert_eq!(
            builder_construction.to_string(),
            "Construction(Builder(ForceFieldKernel(OptimizerBadDirection)))"
        );
        let construction = match &builder_construction {
            AutomaticForceFieldConstructionError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&builder_construction, construction);
        assert_eq!(
            construction.to_string(),
            "Builder(ForceFieldKernel(OptimizerBadDirection))"
        );
        let builder_error = match construction {
            ForceFieldConstructionError::Builder(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(construction, builder_error);
        let kernel_error = match builder_error {
            UffBuilderError::ForceFieldKernel(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(builder_error, kernel_error);
        assert!(std::error::Error::source(kernel_error).is_none());

        let nonbonded_construction = AutomaticForceFieldConstructionError::Construction(
            ForceFieldConstructionError::Nonbonded(NonbondedAssemblyError::PairBuilder(
                UffBuilderError::ForceFieldKernel(ForceFieldKernelError::OptimizerBadDirection),
            )),
        );
        assert_eq!(
            nonbonded_construction.to_string(),
            "Construction(Nonbonded(PairBuilder(ForceFieldKernel(OptimizerBadDirection))))"
        );
        let construction = match &nonbonded_construction {
            AutomaticForceFieldConstructionError::Construction(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(&nonbonded_construction, construction);
        assert_eq!(
            construction.to_string(),
            "Nonbonded(PairBuilder(ForceFieldKernel(OptimizerBadDirection)))"
        );
        let assembly_error = match construction {
            ForceFieldConstructionError::Nonbonded(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(construction, assembly_error);
        let builder_error = match assembly_error {
            NonbondedAssemblyError::PairBuilder(stored) => stored,
            NonbondedAssemblyError::FragmentMapping(_) => unreachable!(),
        };
        assert_stored_child(assembly_error, builder_error);
        let kernel_error = match builder_error {
            UffBuilderError::ForceFieldKernel(stored) => stored,
            _ => unreachable!(),
        };
        assert_stored_child(builder_error, kernel_error);
        assert!(std::error::Error::source(kernel_error).is_none());
    }

    #[test]
    fn uff_error_e05_query_errors_preserve_borrowed_search_causes() {
        // These manually built wrappers test trait dispatch only. The existing
        // B24/B26 tests below continue to exercise actual parse, match-arity,
        // and missing-bond callers.
        let default_parse_child = SmartsParseError::TemplateAttachmentRemap {
            carrier: 29,
            source: cosmolkit_model::TemplateAttachmentOrderError::Empty,
        };
        let default_parse = DefaultTorsionQueryError::Parse(default_parse_child);
        assert_eq!(default_parse.to_string(), format!("{default_parse:?}"));
        let stored_parse = match &default_parse {
            DefaultTorsionQueryError::Parse(child) => child,
            _ => unreachable!(),
        };
        let exposed_parse = std::error::Error::source(&default_parse)
            .expect("default Parse must expose its stored Search child")
            .downcast_ref::<SmartsParseError>()
            .expect("default Parse source keeps its original type");
        assert!(std::ptr::eq(stored_parse, exposed_parse));
        assert_eq!(
            exposed_parse.to_string(),
            "template attachment remap failed for carrier atom 29: template attachment order must contain at least one entry"
        );
        let stored_template_source = match stored_parse {
            SmartsParseError::TemplateAttachmentRemap { source, .. } => source,
            _ => unreachable!(),
        };
        let exposed_template_source = std::error::Error::source(stored_parse)
            .expect("Search parse error keeps its existing nested source")
            .downcast_ref::<cosmolkit_model::TemplateAttachmentOrderError>()
            .expect("Search parse source keeps its original type");
        assert!(std::ptr::eq(
            stored_template_source,
            exposed_template_source
        ));
        assert_eq!(
            exposed_template_source.to_string(),
            "template attachment order must contain at least one entry"
        );

        let default_compile_child =
            QueryCompileError::InvalidGraph("default compile detail".to_owned());
        let default_compile = DefaultTorsionQueryError::Compile(default_compile_child);
        assert_eq!(default_compile.to_string(), format!("{default_compile:?}"));
        let stored_default_compile = match &default_compile {
            DefaultTorsionQueryError::Compile(child) => child,
            _ => unreachable!(),
        };
        let exposed_default_compile = std::error::Error::source(&default_compile)
            .expect("default Compile must expose its stored Search child")
            .downcast_ref::<QueryCompileError>()
            .expect("default Compile source keeps its original type");
        assert!(std::ptr::eq(
            stored_default_compile,
            exposed_default_compile
        ));
        assert_eq!(
            exposed_default_compile.to_string(),
            "query graph is invalid: default compile detail"
        );
        assert!(std::error::Error::source(exposed_default_compile).is_none());

        let nested_default = DefaultTorsionQueryError::Compile(QueryCompileError::InvalidGraph(
            "wrapped default detail".to_owned(),
        ));
        let default_query = TorsionBondQueryError::DefaultQuery(nested_default);
        assert_eq!(default_query.to_string(), format!("{default_query:?}"));
        let stored_default_query = match &default_query {
            TorsionBondQueryError::DefaultQuery(child) => child,
            _ => unreachable!(),
        };
        let exposed_default_query = std::error::Error::source(&default_query)
            .expect("DefaultQuery must expose its stored error")
            .downcast_ref::<DefaultTorsionQueryError>()
            .expect("DefaultQuery source keeps its original type");
        assert!(std::ptr::eq(stored_default_query, exposed_default_query));
        let stored_nested_compile = match stored_default_query {
            DefaultTorsionQueryError::Compile(child) => child,
            _ => unreachable!(),
        };
        let exposed_nested_compile = std::error::Error::source(stored_default_query)
            .expect("wrapped DefaultQuery keeps its existing child")
            .downcast_ref::<QueryCompileError>()
            .expect("wrapped compile source keeps its original type");
        assert!(std::ptr::eq(stored_nested_compile, exposed_nested_compile));

        let torsion_parse_child = SmartsParseError::UnclosedBracket(47);
        let torsion_parse = TorsionBondQueryError::Parse(torsion_parse_child);
        assert_eq!(torsion_parse.to_string(), format!("{torsion_parse:?}"));
        let stored_torsion_parse = match &torsion_parse {
            TorsionBondQueryError::Parse(child) => child,
            _ => unreachable!(),
        };
        let exposed_torsion_parse = std::error::Error::source(&torsion_parse)
            .expect("Parse must expose its stored Search child")
            .downcast_ref::<SmartsParseError>()
            .expect("Parse source keeps its original type");
        assert!(std::ptr::eq(stored_torsion_parse, exposed_torsion_parse));
        assert_eq!(
            exposed_torsion_parse.to_string(),
            "unclosed bracket at position 47"
        );
        assert!(std::error::Error::source(exposed_torsion_parse).is_none());

        let torsion_compile_child =
            QueryCompileError::InvalidGraph("torsion compile detail".to_owned());
        let torsion_compile = TorsionBondQueryError::Compile(torsion_compile_child);
        assert_eq!(torsion_compile.to_string(), format!("{torsion_compile:?}"));
        let stored_torsion_compile = match &torsion_compile {
            TorsionBondQueryError::Compile(child) => child,
            _ => unreachable!(),
        };
        let exposed_torsion_compile = std::error::Error::source(&torsion_compile)
            .expect("Compile must expose its stored Search child")
            .downcast_ref::<QueryCompileError>()
            .expect("Compile source keeps its original type");
        assert!(std::ptr::eq(
            stored_torsion_compile,
            exposed_torsion_compile
        ));
        assert_eq!(
            exposed_torsion_compile.to_string(),
            "query graph is invalid: torsion compile detail"
        );
        assert!(std::error::Error::source(exposed_torsion_compile).is_none());

        let match_child =
            MatchError::Substruct(cosmolkit_search::SubstructMatchError::Unsupported {
                branch: "E05 fixed trait case",
                rdkit_function: "SubstructMatch(E05)",
            });
        let match_error = TorsionBondQueryError::Match(match_child);
        assert_eq!(match_error.to_string(), format!("{match_error:?}"));
        let stored_match = match &match_error {
            TorsionBondQueryError::Match(child) => child,
            _ => unreachable!(),
        };
        let exposed_match = std::error::Error::source(&match_error)
            .expect("Match must expose its stored Search error")
            .downcast_ref::<MatchError>()
            .expect("Match source keeps its original type");
        assert!(std::ptr::eq(stored_match, exposed_match));
        assert_eq!(
            exposed_match.to_string(),
            "RDKit substructure matching branch E05 fixed trait case is unsupported until SubstructMatch(E05) is source-ported"
        );
        let stored_substruct = match exposed_match {
            MatchError::Substruct(child) => child,
            _ => unreachable!(),
        };
        assert_eq!(
            stored_substruct,
            &cosmolkit_search::SubstructMatchError::Unsupported {
                branch: "E05 fixed trait case",
                rdkit_function: "SubstructMatch(E05)",
            }
        );
        assert!(std::error::Error::source(exposed_match).is_none());
        assert!(std::error::Error::source(stored_substruct).is_none());

        let arity_smarts = "[C;D7]~[N;D5]".to_owned();
        let arity_error = TorsionBondQueryError::MatchArity {
            smarts: arity_smarts,
            actual: 7,
        };
        assert_eq!(arity_error.to_string(), format!("{arity_error:?}"));
        assert!(std::error::Error::source(&arity_error).is_none());
        match &arity_error {
            TorsionBondQueryError::MatchArity { smarts, actual } => {
                assert_eq!(smarts, "[C;D7]~[N;D5]");
                assert_eq!(*actual, 7);
            }
            _ => unreachable!(),
        }

        let missing_bond = TorsionBondQueryError::MissingMatchedBond {
            begin_atom_index: 71,
            end_atom_index: 73,
        };
        assert_eq!(missing_bond.to_string(), format!("{missing_bond:?}"));
        assert!(std::error::Error::source(&missing_bond).is_none());
        match &missing_bond {
            TorsionBondQueryError::MissingMatchedBond {
                begin_atom_index,
                end_atom_index,
            } => assert_eq!((*begin_atom_index, *end_atom_index), (71, 73)),
            _ => unreachable!(),
        }
    }

    fn topology(no_implicit: &[bool]) -> TopologyBlock {
        let atoms = no_implicit
            .iter()
            .enumerate()
            .map(|(row, &no_implicit)| {
                let element =
                    Element::from_atomic_number(6).expect("carbon is present in the element table");
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(element).with_no_implicit(no_implicit),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed atom-only topology is valid")
    }

    fn topology_with_hydrogen_counts(
        explicit_hydrogens: &[u8],
        no_implicit: &[bool],
    ) -> TopologyBlock {
        let atoms = explicit_hydrogens
            .iter()
            .zip(no_implicit)
            .enumerate()
            .map(|(row, (&explicit_hydrogens, &no_implicit))| {
                let element =
                    Element::from_atomic_number(6).expect("carbon is present in the element table");
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(element)
                        .with_explicit_hydrogens(explicit_hydrogens)
                        .with_no_implicit(no_implicit),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed hydrogen-count topology is valid")
    }

    fn expected_needs_hydrogen_warning() -> UffTypingDiagnostic {
        UffTypingDiagnostic {
            atom_id: None,
            kind: UffTypingDiagnosticKind::Warning,
            message_prefix: "Molecule does not have explicit Hs. Consider calling AddHs()",
        }
    }

    fn assignment(explicit: &[i32], implicit: &[i32]) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: explicit.to_vec(),
            implicit_hydrogens: implicit.to_vec(),
        }
    }

    fn topology_with_bonds(atom_count: usize, edges: &[(usize, usize, bool)]) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|row| {
                let element =
                    Element::from_atomic_number(6).expect("carbon is present in the element table");
                Atom::from_spec(AtomId::new(row), AtomSpec::new(element))
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end, is_conjugated))| {
                let mut bond = Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                );
                bond.set_conjugated(is_conjugated);
                bond
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed bond topology is structurally valid")
    }

    fn topology_with_bond_orders(
        atom_count: usize,
        edges: &[(usize, usize, BondOrder)],
    ) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|row| {
                let element =
                    Element::from_atomic_number(6).expect("carbon is present in the element table");
                Atom::from_spec(AtomId::new(row), AtomSpec::new(element))
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end, order))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed bond-order topology is structurally valid")
    }

    fn atomic_params(r1: f64, z1: f64) -> AtomicParams {
        AtomicParams {
            r1,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1,
            v1: 0.0,
            u1: 0.0,
            gmp_xi: 1.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn cf3d_bld_borrowed_params(params: &[Option<AtomicParams>]) -> Vec<Option<&AtomicParams>> {
        params.iter().map(Option::as_ref).collect()
    }

    fn cf3d_frag_f21_params(x1: f64, d1: f64) -> AtomicParams {
        let mut params = atomic_params(1.0, 6.0);
        params.x1 = x1;
        params.d1 = d1;
        params
    }

    fn attach_positions<'a>(field: &mut ForceField<'a>, rows: &'a mut [[f64; 3]]) {
        field
            .positions_mut()
            .extend(rows.iter_mut().map(|row| &mut row[..]));
    }

    fn cf3d_bld_b15_topology(
        hybridizations: &[Hybridization],
        edges: &[(usize, usize, BondOrder)],
    ) -> TopologyBlock {
        let atoms = hybridizations
            .iter()
            .enumerate()
            .map(|(row, &hybridization)| {
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(Element::C).with_hybridization(hybridization),
                )
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end, order))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed B15 topology is structurally valid")
    }

    fn cf3d_bld_b15_ring_info(topology: &TopologyBlock) -> RingInfo {
        cosmolkit_core::fast_find_rings(topology).expect("fixed B15 ring input is valid")
    }

    fn cf3d_bld_b15_atomic_params() -> AtomicParams {
        AtomicParams {
            r1: 0.5,
            theta0: std::f64::consts::FRAC_PI_2,
            x1: 1.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 1.0,
            v1: 0.0,
            u1: 0.0,
            gmp_xi: 1.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn cf3d_bld_b15_full_params(atom_count: usize) -> Vec<Option<AtomicParams>> {
        (0..atom_count)
            .map(|_| Some(cf3d_bld_b15_atomic_params()))
            .collect()
    }

    fn cf3d_bld_b15_points(atom_count: usize) -> Vec<[f64; 3]> {
        let fixed_points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 1.0, 2.0],
            [2.0, -1.0, 3.0],
            [-2.0, -3.0, 1.0],
        ];
        fixed_points[..atom_count].to_vec()
    }

    fn cf3d_bld_b15_append_expected_terms(
        field: &mut ForceField<'_>,
        params: &[Option<AtomicParams>],
        terms: &[(usize, usize, usize, u32)],
    ) {
        for &(i, j, k, order) in terms {
            let contribution = {
                let positions = field.positions();
                AngleBendContrib::new(
                    positions,
                    u32::try_from(i).expect("fixed expected index fits source width"),
                    u32::try_from(j).expect("fixed expected index fits source width"),
                    u32::try_from(k).expect("fixed expected index fits source width"),
                    1.0,
                    1.0,
                    params[i].as_ref().expect("fixed source parameter at i"),
                    params[j].as_ref().expect("fixed source parameter at j"),
                    params[k].as_ref().expect("fixed source parameter at k"),
                    order,
                )
            }
            .expect("source-fixed expected angle is valid");
            field.add_contribution(Box::new(contribution));
        }
    }

    fn cf3d_bld_b15_assert_kernel_matches_terms(
        topology: &TopologyBlock,
        params: &[Option<AtomicParams>],
        rings: &RingInfo,
        points: &[[f64; 3]],
        expected_terms: &[(usize, usize, usize, u32)],
    ) {
        let coordinates = points
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points.to_vec();
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);
        let borrowed_params = cf3d_bld_borrowed_params(params);
        add_angles(topology, &borrowed_params, rings, &mut actual)
            .expect("source-fixed B15 topology appends its expected terms");
        actual.initialize().expect("actual B15 field initializes");

        let mut expected_points = points.to_vec();
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        cf3d_bld_b15_append_expected_terms(&mut expected, params, expected_terms);
        expected
            .initialize()
            .expect("source-fixed field initializes");

        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &coordinates)
            .expect("actual B15 field energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &coordinates)
            .expect("source-fixed expected energy evaluates");
        assert_eq!(
            actual_energy, expected_energy,
            "fixed source term list {expected_terms:?}"
        );

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &coordinates, &mut actual_gradient)
            .expect("actual B15 field gradient evaluates");
        cf3d_bld_b05_calc_grad(&mut expected, &coordinates, &mut expected_gradient)
            .expect("source-fixed expected gradient evaluates");
        assert_eq!(actual_gradient, expected_gradient);
    }

    #[test]
    fn cf3d_bld_b15_parameter_masks_and_error_order_are_source_defined() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 3],
            &[(0, 1, BondOrder::DativeLeft), (0, 2, BondOrder::Single)],
        );
        let rings = cf3d_bld_b15_ring_info(&topology);
        let points = cf3d_bld_b15_points(3);
        let coordinates = points
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let masks = [
            (false, false, false),
            (false, false, true),
            (false, true, false),
            (false, true, true),
            (true, false, false),
            (true, false, true),
            (true, true, false),
            (true, true, true),
        ];

        for (has_i, has_j, has_k) in masks {
            let mut params = vec![None; 3];
            params[0] = has_j.then(cf3d_bld_b15_atomic_params);
            params[1] = has_i.then(cf3d_bld_b15_atomic_params);
            params[2] = has_k.then(cf3d_bld_b15_atomic_params);

            let mut field = ForceField::new(3);
            let mut actual_points = points.clone();
            attach_positions(&mut field, &mut actual_points);
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            let result = add_angles(&topology, &borrowed_params, &rings, &mut field);

            if has_i && has_j && has_k {
                assert_eq!(
                    result,
                    Err(UffBuilderError::Valence(ValenceError::BadBondType {
                        bond: None,
                        order: BondOrder::DativeLeft,
                    }))
                );
            } else {
                assert_eq!(result, Ok(()));
                field.initialize().expect("empty source mask initializes");
                assert_eq!(
                    cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                        .expect("empty source mask energy evaluates"),
                    0.0
                );
                let mut gradient = vec![0.0; coordinates.len()];
                cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
                    .expect("empty source mask gradient evaluates");
                assert_eq!(gradient, vec![0.0; coordinates.len()]);
            }
        }

        let short_params = vec![Some(cf3d_bld_b15_atomic_params()); 2];
        let mut field = ForceField::new(3);
        let short_borrowed_params = cf3d_bld_borrowed_params(&short_params);
        assert_eq!(
            add_angles(&topology, &short_borrowed_params, &rings, &mut field),
            Err(UffBuilderError::ParamsLengthMismatch {
                atoms: 3,
                params: 2,
            })
        );
    }

    #[test]
    fn cf3d_bld_b15_constructor_bounds_remain_typed() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 3],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        let rings = cf3d_bld_b15_ring_info(&topology);
        let params = cf3d_bld_b15_full_params(3);
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut field = ForceField::new(3);
        let mut short_points = cf3d_bld_b15_points(2);
        attach_positions(&mut field, &mut short_points);

        assert_eq!(
            add_angles(&topology, &borrowed_params, &rings, &mut field),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::Third,
                    index: 2,
                    upper_bound: 2,
                }
            ))
        );
    }

    #[test]
    fn cf3d_bld_b15_degree_matrix_and_sp3d_exclusion_match_source() {
        let empty = cf3d_bld_b15_topology(&[Hybridization::Sp3], &[]);
        cf3d_bld_b15_assert_kernel_matches_terms(
            &empty,
            &cf3d_bld_b15_full_params(1),
            &cf3d_bld_b15_ring_info(&empty),
            &cf3d_bld_b15_points(1),
            &[],
        );

        let degree_one =
            cf3d_bld_b15_topology(&[Hybridization::Sp3; 2], &[(0, 1, BondOrder::Single)]);
        cf3d_bld_b15_assert_kernel_matches_terms(
            &degree_one,
            &cf3d_bld_b15_full_params(2),
            &cf3d_bld_b15_ring_info(&degree_one),
            &cf3d_bld_b15_points(2),
            &[],
        );

        let degree_two = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 3],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &degree_two,
            &cf3d_bld_b15_full_params(3),
            &cf3d_bld_b15_ring_info(&degree_two),
            &cf3d_bld_b15_points(3),
            &[(1, 0, 2, 0)],
        );

        let degree_five = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 6],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &degree_five,
            &cf3d_bld_b15_full_params(6),
            &cf3d_bld_b15_ring_info(&degree_five),
            &cf3d_bld_b15_points(6),
            &[
                (1, 0, 2, 0),
                (1, 0, 3, 0),
                (1, 0, 4, 0),
                (1, 0, 5, 0),
                (2, 0, 3, 0),
                (2, 0, 4, 0),
                (2, 0, 5, 0),
                (3, 0, 4, 0),
                (3, 0, 5, 0),
                (4, 0, 5, 0),
            ],
        );

        let degree_six = cf3d_bld_b15_topology(
            &[Hybridization::Sp3d; 7],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 4, BondOrder::Single),
                (0, 5, BondOrder::Single),
                (0, 6, BondOrder::Single),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &degree_six,
            &cf3d_bld_b15_full_params(7),
            &cf3d_bld_b15_ring_info(&degree_six),
            &cf3d_bld_b15_points(7),
            &[
                (1, 0, 2, 0),
                (1, 0, 3, 0),
                (1, 0, 4, 0),
                (1, 0, 5, 0),
                (1, 0, 6, 0),
                (2, 0, 3, 0),
                (2, 0, 4, 0),
                (2, 0, 5, 0),
                (2, 0, 6, 0),
                (3, 0, 4, 0),
                (3, 0, 5, 0),
                (3, 0, 6, 0),
                (4, 0, 5, 0),
                (4, 0, 6, 0),
                (5, 0, 6, 0),
            ],
        );

        let sp3d_degree_five = cf3d_bld_b15_topology(
            &[Hybridization::Sp3d; 6],
            &[
                (0, 1, BondOrder::DativeLeft),
                (0, 2, BondOrder::DativeLeft),
                (0, 3, BondOrder::DativeLeft),
                (0, 4, BondOrder::DativeLeft),
                (0, 5, BondOrder::DativeLeft),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &sp3d_degree_five,
            &cf3d_bld_b15_full_params(6),
            &cf3d_bld_b15_ring_info(&sp3d_degree_five),
            &cf3d_bld_b15_points(6),
            &[],
        );
    }

    #[test]
    fn cf3d_bld_b15_all_b14_orders_reach_real_kernel_terms() {
        for (hybridization, order) in [
            (Hybridization::Sp3, 0),
            (Hybridization::Sp, 1),
            (Hybridization::Sp2, 3),
            (Hybridization::Sp3d2, 4),
        ] {
            let topology = cf3d_bld_b15_topology(
                &[hybridization; 3],
                &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
            );
            cf3d_bld_b15_assert_kernel_matches_terms(
                &topology,
                &cf3d_bld_b15_full_params(3),
                &cf3d_bld_b15_ring_info(&topology),
                &cf3d_bld_b15_points(3),
                &[(1, 0, 2, order)],
            );
        }

        let triangle = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 3],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &triangle,
            &cf3d_bld_b15_full_params(3),
            &cf3d_bld_b15_ring_info(&triangle),
            &cf3d_bld_b15_points(3),
            &[(1, 0, 2, 35), (0, 1, 2, 35), (1, 2, 0, 35)],
        );

        let triangle_with_terminal = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 4],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &triangle_with_terminal,
            &cf3d_bld_b15_full_params(4),
            &cf3d_bld_b15_ring_info(&triangle_with_terminal),
            &cf3d_bld_b15_points(4),
            &[
                (1, 0, 2, 35),
                (1, 0, 3, 30),
                (2, 0, 3, 30),
                (0, 1, 2, 35),
                (1, 2, 0, 35),
            ],
        );

        let square = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 4],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 0, BondOrder::Single),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &square,
            &cf3d_bld_b15_full_params(4),
            &cf3d_bld_b15_ring_info(&square),
            &cf3d_bld_b15_points(4),
            &[(1, 0, 3, 45), (0, 1, 2, 45), (1, 2, 3, 45), (2, 3, 0, 45)],
        );

        let square_with_terminal = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 5],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 0, BondOrder::Single),
                (0, 4, BondOrder::Single),
            ],
        );
        cf3d_bld_b15_assert_kernel_matches_terms(
            &square_with_terminal,
            &cf3d_bld_b15_full_params(5),
            &cf3d_bld_b15_ring_info(&square_with_terminal),
            &cf3d_bld_b15_points(5),
            &[
                (1, 0, 3, 45),
                (1, 0, 4, 40),
                (3, 0, 4, 40),
                (0, 1, 2, 45),
                (1, 2, 3, 45),
                (2, 3, 0, 45),
            ],
        );
    }

    #[test]
    fn cf3d_bld_b15_neighbor_slot_order_preserves_partial_append_prefix() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (0, 1, BondOrder::DativeLeft),
            ],
        );
        let rings = cf3d_bld_b15_ring_info(&topology);
        let params = cf3d_bld_b15_full_params(4);
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let coordinates = [
            [0.0, 0.0, 0.0],
            [-1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        ];
        let flat_coordinates = coordinates
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = coordinates;
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);

        assert_eq!(
            add_angles(&topology, &borrowed_params, &rings, &mut actual),
            Err(UffBuilderError::Valence(ValenceError::BadBondType {
                bond: None,
                order: BondOrder::DativeLeft,
            }))
        );
        actual
            .initialize()
            .expect("successful source-prefix contribution initializes");
        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &flat_coordinates)
            .expect("source-prefix energy evaluates");
        let mut actual_gradient = vec![0.0; flat_coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &flat_coordinates, &mut actual_gradient)
            .expect("source-prefix gradient evaluates");

        let mut expected_points = coordinates;
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        cf3d_bld_b15_append_expected_terms(&mut expected, &params, &[(2, 0, 3, 0)]);
        expected
            .initialize()
            .expect("fixed one-term source prefix initializes");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &flat_coordinates)
            .expect("fixed one-term energy evaluates");
        let mut expected_gradient = vec![0.0; flat_coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut expected, &flat_coordinates, &mut expected_gradient)
            .expect("fixed one-term gradient evaluates");
        assert_eq!(actual_energy, expected_energy);
        assert_eq!(actual_gradient, expected_gradient);
    }

    #[test]
    fn cf3d_bld_b01_empty_rows_remain_empty() {
        let topology = topology(&[]);
        let assignment = assignment(&[], &[]);
        let before = assignment.clone();

        let prepared = prepare_typing_valence(&topology, &assignment)
            .expect("empty source-aligned rows are valid");

        assert!(prepared.total_valences.is_empty());
        assert_eq!(
            prepared.implicit_hydrogens().collect::<Vec<_>>(),
            Vec::<i32>::new()
        );
        assert_eq!(assignment, before);
        assert!(std::ptr::eq(prepared.assignment, &assignment));
        assert!(std::ptr::eq(prepared.topology, &topology));
    }

    #[test]
    fn cf3d_bld_b01_single_row_sums_explicit_then_implicit() {
        let topology = topology(&[false]);
        let assignment = assignment(&[6], &[2]);
        let before = assignment.clone();

        let prepared = prepare_typing_valence(&topology, &assignment)
            .expect("source-computed single row is valid");

        assert_eq!(prepared.total_valences, [8]);
        assert_eq!(prepared.implicit_hydrogens().collect::<Vec<_>>(), [2]);
        assert_eq!(assignment, before);
        assert!(std::ptr::eq(prepared.assignment, &assignment));
    }

    #[test]
    fn cf3d_bld_b01_multiple_rows_keep_atom_order_and_no_implicit_projection() {
        let topology = topology(&[false, true, false]);
        let assignment = assignment(&[4, 3, 1], &[1, -1, 2]);

        let prepared = prepare_typing_valence(&topology, &assignment)
            .expect("noImplicit suppresses the second row's implicit cache");

        assert_eq!(prepared.total_valences, [5, 3, 3]);
        assert_eq!(prepared.implicit_hydrogens().collect::<Vec<_>>(), [1, 0, 2]);
    }

    #[test]
    fn cf3d_bld_b01_explicit_length_mismatch_is_typed() {
        let topology = topology(&[false]);
        let assignment = assignment(&[], &[0]);

        assert_eq!(
            prepare_typing_valence(&topology, &assignment).unwrap_err(),
            UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::Explicit,
                expected: 1,
                actual: 0,
            }
        );
        assert!(assignment.explicit_valence.is_empty());
    }

    #[test]
    fn cf3d_bld_b01_implicit_length_mismatch_is_typed() {
        let topology = topology(&[false]);
        let assignment = assignment(&[0], &[]);

        assert_eq!(
            prepare_typing_valence(&topology, &assignment).unwrap_err(),
            UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::ImplicitHydrogen,
                expected: 1,
                actual: 0,
            }
        );
        assert!(assignment.implicit_hydrogens.is_empty());
    }

    #[test]
    fn cf3d_bld_b01_missing_explicit_cache_is_a_source_precondition_error() {
        let topology = topology(&[false]);
        let assignment = assignment(&[-1], &[0]);
        let before = assignment.clone();

        assert_eq!(
            prepare_typing_valence(&topology, &assignment).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: -1,
            }
        );
        assert_eq!(assignment, before);
    }

    #[test]
    fn cf3d_bld_b01_missing_implicit_cache_is_a_source_precondition_error() {
        let topology = topology(&[false]);
        let assignment = assignment(&[0], &[-1]);

        assert_eq!(
            prepare_typing_valence(&topology, &assignment).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: -1,
            }
        );
    }

    #[test]
    fn cf3d_bld_b01_no_implicit_skips_implicit_cache_but_not_explicit_getter() {
        let topology = topology(&[true]);
        let ignored_implicit = assignment(&[7], &[-1]);
        let prepared = prepare_typing_valence(&topology, &ignored_implicit)
            .expect("getNumImplicitHs returns zero before inspecting its cache");

        assert_eq!(prepared.total_valences, [7]);
        assert_eq!(prepared.implicit_hydrogens().collect::<Vec<_>>(), [0]);

        let missing_explicit = assignment(&[-1], &[-1]);
        assert_eq!(
            prepare_typing_valence(&topology, &missing_explicit).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: -1,
            }
        );
    }

    #[test]
    fn cf3d_bld_b01_signed_byte_maxima_sum_without_narrowing() {
        let topology = topology(&[false]);
        let assignment = assignment(&[127], &[127]);

        let prepared = prepare_typing_valence(&topology, &assignment)
            .expect("127 is the largest nonnegative signed-byte source cache value");

        assert_eq!(prepared.total_valences, [254]);
        assert_eq!(prepared.implicit_hydrogens().collect::<Vec<_>>(), [127]);
    }

    #[test]
    fn cf3d_bld_b01_explicit_cache_above_signed_byte_maximum_is_typed() {
        let topology = topology(&[false]);
        let assignment = assignment(&[128], &[0]);

        assert_eq!(
            prepare_typing_valence(&topology, &assignment).unwrap_err(),
            UffBuilderError::SourceValenceOutOfRange {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: 128,
            }
        );
    }

    #[test]
    fn cf3d_bld_b01_implicit_cache_above_signed_byte_maximum_is_typed() {
        let topology = topology(&[false]);
        let assignment = assignment(&[0], &[128]);

        assert_eq!(
            prepare_typing_valence(&topology, &assignment).unwrap_err(),
            UffBuilderError::SourceValenceOutOfRange {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: 128,
            }
        );
    }

    #[test]
    fn uff_prepare_p01_empty_rows_and_shape_errors_keep_source_order() {
        let empty_topology = topology(&[]);
        let empty_assignment = assignment(&[], &[]);
        let empty = prepare_typing_valence(&empty_topology, &empty_assignment)
            .expect("empty aligned source rows are valid");
        assert!(empty.total_valences.is_empty());

        let one_atom = topology(&[false]);
        let both_short = assignment(&[], &[]);
        assert_eq!(
            prepare_typing_valence(&one_atom, &both_short).unwrap_err(),
            UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::Explicit,
                expected: 1,
                actual: 0,
            }
        );

        let explicit_aligned = assignment(&[0], &[]);
        assert_eq!(
            prepare_typing_valence(&one_atom, &explicit_aligned).unwrap_err(),
            UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::ImplicitHydrogen,
                expected: 1,
                actual: 0,
            }
        );

        let both_invalid = assignment(&[-1], &[-1]);
        assert_eq!(
            prepare_typing_valence(&one_atom, &both_invalid).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: -1,
            }
        );
    }

    #[test]
    fn uff_prepare_p01_explicit_cache_keeps_signed_byte_bounds() {
        let implicit_allowed = topology(&[false]);
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[0], &[0]))
                .expect("zero explicit and implicit caches are valid")
                .total_valences,
            [0]
        );
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[127], &[0]))
                .expect("127 is the largest stored signed-byte cache value")
                .total_valences,
            [127]
        );
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[-1], &[0])).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: -1,
            }
        );
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[128], &[0])).unwrap_err(),
            UffBuilderError::SourceValenceOutOfRange {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: 128,
            }
        );

        let no_implicit = topology(&[true]);
        assert_eq!(
            prepare_typing_valence(&no_implicit, &assignment(&[0], &[-1]))
                .expect("explicit cache is checked while noImplicit skips implicit cache")
                .total_valences,
            [0]
        );
        assert_eq!(
            prepare_typing_valence(&no_implicit, &assignment(&[-1], &[-1])).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::Explicit,
                value: -1,
            }
        );
    }

    #[test]
    fn uff_prepare_p01_implicit_cache_keeps_no_implicit_short_circuit() {
        let implicit_allowed = topology(&[false]);
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[0], &[0]))
                .expect("zero implicit cache is valid")
                .total_valences,
            [0]
        );
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[0], &[127]))
                .expect("127 is the largest stored signed-byte cache value")
                .total_valences,
            [127]
        );
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[0], &[-1])).unwrap_err(),
            UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: -1,
            }
        );
        assert_eq!(
            prepare_typing_valence(&implicit_allowed, &assignment(&[0], &[128])).unwrap_err(),
            UffBuilderError::SourceValenceOutOfRange {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: 128,
            }
        );

        let no_implicit = topology(&[true]);
        assert_eq!(
            prepare_typing_valence(&no_implicit, &assignment(&[0], &[128]))
                .expect("noImplicit reads zero without validating the implicit cache")
                .total_valences,
            [0]
        );
    }

    #[test]
    fn uff_prepare_p01_cached_values_full32_product() {
        #[derive(Clone, Copy)]
        enum ExpectedFailure {
            SourcePrecondition {
                field: PreparedValenceField,
                value: i32,
            },
            SourceOutOfRange {
                field: PreparedValenceField,
                value: i32,
            },
        }

        use ExpectedFailure::{SourceOutOfRange, SourcePrecondition};
        type Expected = Result<i32, ExpectedFailure>;

        let cases: [(i32, bool, [Expected; 4]); 8] = [
            (
                -1,
                false,
                [Err(SourcePrecondition {
                    field: PreparedValenceField::Explicit,
                    value: -1,
                }); 4],
            ),
            (
                -1,
                true,
                [Err(SourcePrecondition {
                    field: PreparedValenceField::Explicit,
                    value: -1,
                }); 4],
            ),
            (
                0,
                false,
                [
                    Err(SourcePrecondition {
                        field: PreparedValenceField::ImplicitHydrogen,
                        value: -1,
                    }),
                    Ok(0),
                    Ok(127),
                    Err(SourceOutOfRange {
                        field: PreparedValenceField::ImplicitHydrogen,
                        value: 128,
                    }),
                ],
            ),
            (0, true, [Ok(0); 4]),
            (
                127,
                false,
                [
                    Err(SourcePrecondition {
                        field: PreparedValenceField::ImplicitHydrogen,
                        value: -1,
                    }),
                    Ok(127),
                    Ok(254),
                    Err(SourceOutOfRange {
                        field: PreparedValenceField::ImplicitHydrogen,
                        value: 128,
                    }),
                ],
            ),
            (127, true, [Ok(127); 4]),
            (
                128,
                false,
                [Err(SourceOutOfRange {
                    field: PreparedValenceField::Explicit,
                    value: 128,
                }); 4],
            ),
            (
                128,
                true,
                [Err(SourceOutOfRange {
                    field: PreparedValenceField::Explicit,
                    value: 128,
                }); 4],
            ),
        ];
        let implicit_values = [-1, 0, 127, 128];
        let mut real_calls = 0;

        for (explicit, no_implicit, expected_by_implicit) in cases {
            for (implicit, expected) in implicit_values.into_iter().zip(expected_by_implicit) {
                real_calls += 1;
                let input_topology = topology(&[no_implicit]);
                let input_assignment = assignment(&[explicit], &[implicit]);
                let topology_before = input_topology.clone();
                let assignment_before = input_assignment.clone();

                reset_typing_valence_projection_vec_constructions();
                reset_conjugation_projection_vec_constructions();
                let actual = UffAtomStateRef::cached(&input_topology, &input_assignment);

                match expected {
                    Ok(expected_total) => {
                        let state = actual.expect("literal table marks this cached row valid");
                        match state {
                            UffAtomStateRef::Cached {
                                topology: borrowed_topology,
                                assignment: borrowed_assignment,
                            } => {
                                assert!(std::ptr::eq(borrowed_topology, &input_topology));
                                assert!(std::ptr::eq(borrowed_assignment, &input_assignment));
                            }
                            UffAtomStateRef::SuppliedRows { .. } => {
                                panic!("cached factory returned supplied rows")
                            }
                        }
                        assert_eq!(
                            state.total_valence_at(0),
                            expected_total,
                            "explicit={explicit}, implicit={implicit}, noImplicit={no_implicit}"
                        );
                    }
                    Err(failure) => {
                        let expected_error = match failure {
                            SourcePrecondition { field, value } => {
                                UffBuilderError::SourceValencePrecondition {
                                    atom_id: AtomId::new(0),
                                    field,
                                    value,
                                }
                            }
                            SourceOutOfRange { field, value } => {
                                UffBuilderError::SourceValenceOutOfRange {
                                    atom_id: AtomId::new(0),
                                    field,
                                    value,
                                }
                            }
                        };
                        assert_eq!(
                            actual.unwrap_err(),
                            expected_error,
                            "explicit={explicit}, implicit={implicit}, noImplicit={no_implicit}"
                        );
                    }
                }

                assert_eq!(input_topology, topology_before);
                assert_eq!(input_assignment, assignment_before);
                assert_eq!(typing_valence_projection_vec_constructions(), 0);
                assert_eq!(conjugation_projection_vec_constructions(), 0);
            }
        }

        assert_eq!(real_calls, 32);
    }

    #[test]
    fn uff_prepare_p01_cached_lengths_full8_product() {
        #[derive(Clone, Copy)]
        enum ExpectedShape {
            ExplicitLengthMismatch,
            ImplicitLengthMismatch,
            TotalZero,
        }

        let cases = [
            (false, 0, 0, ExpectedShape::ExplicitLengthMismatch),
            (false, 0, 1, ExpectedShape::ExplicitLengthMismatch),
            (false, 1, 0, ExpectedShape::ImplicitLengthMismatch),
            (false, 1, 1, ExpectedShape::TotalZero),
            (true, 0, 0, ExpectedShape::ExplicitLengthMismatch),
            (true, 0, 1, ExpectedShape::ExplicitLengthMismatch),
            (true, 1, 0, ExpectedShape::ImplicitLengthMismatch),
            (true, 1, 1, ExpectedShape::TotalZero),
        ];
        let mut real_calls = 0;

        for (no_implicit, explicit_len, implicit_len, expected) in cases {
            real_calls += 1;
            let input_topology = topology(&[no_implicit]);
            let input_assignment = assignment(&vec![0; explicit_len], &vec![0; implicit_len]);
            let topology_before = input_topology.clone();
            let assignment_before = input_assignment.clone();

            reset_typing_valence_projection_vec_constructions();
            reset_conjugation_projection_vec_constructions();
            let actual = UffAtomStateRef::cached(&input_topology, &input_assignment);

            match expected {
                ExpectedShape::ExplicitLengthMismatch => {
                    assert_eq!(
                        actual.unwrap_err(),
                        UffBuilderError::ValenceAssignmentLengthMismatch {
                            field: PreparedValenceField::Explicit,
                            expected: 1,
                            actual: 0,
                        },
                        "explicit_len={explicit_len}, implicit_len={implicit_len}, noImplicit={no_implicit}"
                    );
                }
                ExpectedShape::ImplicitLengthMismatch => {
                    assert_eq!(
                        actual.unwrap_err(),
                        UffBuilderError::ValenceAssignmentLengthMismatch {
                            field: PreparedValenceField::ImplicitHydrogen,
                            expected: 1,
                            actual: 0,
                        },
                        "explicit_len={explicit_len}, implicit_len={implicit_len}, noImplicit={no_implicit}"
                    );
                }
                ExpectedShape::TotalZero => {
                    let state = actual.expect("literal table marks aligned zero rows valid");
                    match state {
                        UffAtomStateRef::Cached {
                            topology: borrowed_topology,
                            assignment: borrowed_assignment,
                        } => {
                            assert!(std::ptr::eq(borrowed_topology, &input_topology));
                            assert!(std::ptr::eq(borrowed_assignment, &input_assignment));
                        }
                        UffAtomStateRef::SuppliedRows { .. } => {
                            panic!("cached factory returned supplied rows")
                        }
                    }
                    assert_eq!(state.total_valence_at(0), 0);
                }
            }

            assert_eq!(input_topology, topology_before);
            assert_eq!(input_assignment, assignment_before);
            assert_eq!(typing_valence_projection_vec_constructions(), 0);
            assert_eq!(conjugation_projection_vec_constructions(), 0);
        }

        assert_eq!(real_calls, 8);

        let empty_topology = topology(&[]);
        let empty_assignment = assignment(&[], &[]);
        let empty_topology_before = empty_topology.clone();
        let empty_assignment_before = empty_assignment.clone();
        reset_typing_valence_projection_vec_constructions();
        reset_conjugation_projection_vec_constructions();
        let empty_state = UffAtomStateRef::cached(&empty_topology, &empty_assignment)
            .expect("separate aligned-empty cached state is valid");
        match empty_state {
            UffAtomStateRef::Cached {
                topology: borrowed_topology,
                assignment: borrowed_assignment,
            } => {
                assert!(std::ptr::eq(borrowed_topology, &empty_topology));
                assert!(std::ptr::eq(borrowed_assignment, &empty_assignment));
            }
            UffAtomStateRef::SuppliedRows { .. } => {
                panic!("cached factory returned supplied rows")
            }
        }
        assert_eq!(empty_topology, empty_topology_before);
        assert_eq!(empty_assignment, empty_assignment_before);
        assert_eq!(typing_valence_projection_vec_constructions(), 0);
        assert_eq!(conjugation_projection_vec_constructions(), 0);
    }

    #[test]
    fn cf3d_bld_b02_empty_and_isolated_atoms_have_no_conjugated_bonds() {
        let empty = topology(&[]);
        assert_eq!(
            prepare_conjugated_presence(&empty).expect("empty topology is valid"),
            Vec::<bool>::new()
        );

        let isolates = topology(&[false, false, false]);
        let before = isolates.clone();
        assert_eq!(
            prepare_conjugated_presence(&isolates).expect("isolated atoms are valid"),
            [false, false, false]
        );
        assert_eq!(isolates, before);
    }

    #[test]
    fn cf3d_bld_b02_nonconjugated_bonds_leave_all_rows_false() {
        let topology = topology_with_bonds(4, &[(0, 1, false), (3, 2, false)]);

        assert_eq!(
            prepare_conjugated_presence(&topology).expect("nonconjugated bonds are valid"),
            [false, false, false, false]
        );
    }

    #[test]
    fn cf3d_bld_b02_one_conjugated_bond_marks_only_its_endpoints() {
        let topology = topology_with_bonds(4, &[(0, 2, true), (3, 1, false)]);

        assert_eq!(
            prepare_conjugated_presence(&topology).expect("one conjugated bond is valid"),
            [true, false, true, false]
        );
    }

    #[test]
    fn cf3d_bld_b02_all_conjugated_bonds_cover_both_endpoint_orientations() {
        let topology = topology_with_bonds(4, &[(0, 1, true), (3, 2, true)]);

        assert_eq!(
            prepare_conjugated_presence(&topology).expect("all bonds are conjugated"),
            [true, true, true, true]
        );
    }

    #[test]
    fn cf3d_bld_b02_adjacency_mismatch_is_typed_and_input_is_unchanged() {
        let mut topology = topology_with_bonds(2, &[(0, 1, true)]);
        topology.adjacency = AdjacencyList::default();
        let before = topology.clone();

        assert_eq!(
            prepare_conjugated_presence(&topology).unwrap_err(),
            UffBuilderError::TopologyValidation(TopologyValidationError::AdjacencyMismatch)
        );
        assert_eq!(topology, before);
    }

    #[test]
    fn cf3d_bld_b03_explicit_implicit_hydrogen_truth_matrix() {
        let cases = [(0, 0, false), (0, 1, true), (1, 0, true), (1, 1, true)];

        for (explicit_hydrogens, implicit_hydrogens, expected) in cases {
            let topology = topology_with_hydrogen_counts(&[explicit_hydrogens], &[false]);
            let topology_before = topology.clone();
            let implicit_hydrogens = [implicit_hydrogens];
            let implicit_before = implicit_hydrogens;
            let mut diagnostics = Vec::new();

            assert_eq!(
                needs_hydrogens_warning(&topology, &implicit_hydrogens, &mut diagnostics),
                Ok(expected),
                "explicit={explicit_hydrogens}, implicit={implicit_hydrogens:?}"
            );
            assert_eq!(
                diagnostics,
                if expected {
                    vec![expected_needs_hydrogen_warning()]
                } else {
                    Vec::new()
                }
            );
            assert_eq!(topology, topology_before);
            assert_eq!(implicit_hydrogens, implicit_before);
        }
    }

    #[test]
    fn cf3d_bld_b03_empty_and_explicit_neighbor_hydrogen_do_not_warn() {
        let empty = topology(&[]);
        let mut diagnostics = Vec::new();
        assert_eq!(
            needs_hydrogens_warning(&empty, &[], &mut diagnostics),
            Ok(false)
        );
        assert!(diagnostics.is_empty());

        let carbon =
            Element::from_atomic_number(6).expect("carbon is present in the element table");
        let hydrogen =
            Element::from_atomic_number(1).expect("hydrogen is present in the element table");
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(carbon).with_no_implicit(true)),
            Atom::from_spec(
                AtomId::new(1),
                AtomSpec::new(hydrogen).with_no_implicit(true),
            ),
        ];
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        );
        let neighbor_hydrogen =
            TopologyBlock::try_from_parts(atoms, vec![bond], Vec::new(), Vec::new())
                .expect("fixed carbon-hydrogen topology is valid");
        let before = neighbor_hydrogen.clone();

        assert_eq!(
            needs_hydrogens_warning(&neighbor_hydrogen, &[-1, -1], &mut diagnostics),
            Ok(false),
            "includeNeighbors=false excludes a separately bonded H atom"
        );
        assert!(diagnostics.is_empty());
        assert_eq!(neighbor_hydrogen, before);
    }

    #[test]
    fn cf3d_bld_b03_no_implicit_skips_cache_before_explicit_total() {
        let no_hydrogens = topology_with_hydrogen_counts(&[0], &[true]);
        let mut diagnostics = Vec::new();
        assert_eq!(
            needs_hydrogens_warning(&no_hydrogens, &[-1], &mut diagnostics),
            Ok(false),
            "noImplicit returns zero without inspecting its cache"
        );
        assert!(diagnostics.is_empty());

        let explicit_hydrogen_count = topology_with_hydrogen_counts(&[1], &[true]);
        assert_eq!(
            needs_hydrogens_warning(&explicit_hydrogen_count, &[-1], &mut diagnostics),
            Ok(true),
            "the explicit H total is tested after noImplicit supplies zero"
        );
        assert_eq!(diagnostics, vec![expected_needs_hydrogen_warning()]);
    }

    #[test]
    fn cf3d_bld_b03_first_and_last_positive_rows_keep_source_scan_order() {
        let first_positive = topology_with_hydrogen_counts(&[1, 0, 0], &[false; 3]);
        let mut diagnostics = Vec::new();
        assert_eq!(
            needs_hydrogens_warning(&first_positive, &[0, 0, -1], &mut diagnostics),
            Ok(true),
            "a later missing getter is not reached after the first positive total"
        );
        assert_eq!(diagnostics, vec![expected_needs_hydrogen_warning()]);

        diagnostics.clear();
        let last_positive = topology_with_hydrogen_counts(&[0, 0, 0], &[false; 3]);
        assert_eq!(
            needs_hydrogens_warning(&last_positive, &[0, 0, 1], &mut diagnostics),
            Ok(true),
            "the last row is checked when all earlier totals are zero"
        );
        assert_eq!(diagnostics, vec![expected_needs_hydrogen_warning()]);
    }

    #[test]
    fn cf3d_bld_b03_getter_failures_precede_positive_total_and_warning() {
        let explicit_before_failed_implicit = topology_with_hydrogen_counts(&[1], &[false]);
        let mut diagnostics = Vec::new();
        assert_eq!(
            needs_hydrogens_warning(&explicit_before_failed_implicit, &[-1], &mut diagnostics,),
            Err(UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: -1,
            }),
            "getTotalNumHs obtains implicit H before testing the explicit-plus-implicit sum"
        );
        assert!(diagnostics.is_empty());

        let failed_first_row = topology_with_hydrogen_counts(&[0, 0], &[false; 2]);
        assert_eq!(
            needs_hydrogens_warning(&failed_first_row, &[-1, 1], &mut diagnostics),
            Err(UffBuilderError::SourceValencePrecondition {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: -1,
            }),
            "a failed earlier getter exits before a later positive row"
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_bld_b03_alignment_and_source_cache_range_errors_are_typed() {
        let topology = topology_with_hydrogen_counts(&[0], &[false]);
        let mut diagnostics = Vec::new();

        assert_eq!(
            needs_hydrogens_warning(&topology, &[], &mut diagnostics),
            Err(UffBuilderError::ValenceAssignmentLengthMismatch {
                field: PreparedValenceField::ImplicitHydrogen,
                expected: 1,
                actual: 0,
            })
        );
        assert!(diagnostics.is_empty());

        assert_eq!(
            needs_hydrogens_warning(&topology, &[128], &mut diagnostics),
            Err(UffBuilderError::SourceValenceOutOfRange {
                atom_id: AtomId::new(0),
                field: PreparedValenceField::ImplicitHydrogen,
                value: 128,
            })
        );
        assert!(diagnostics.is_empty());
    }

    #[test]
    fn cf3d_bld_b03_warning_is_one_global_event_after_existing_diagnostics() {
        let topology = topology_with_hydrogen_counts(&[1, 1], &[false; 2]);
        let existing = UffTypingDiagnostic {
            atom_id: Some(AtomId::new(99)),
            kind: UffTypingDiagnosticKind::Error,
            message_prefix: "existing source diagnostic",
        };
        let mut diagnostics = vec![existing];

        assert_eq!(
            needs_hydrogens_warning(&topology, &[0, 0], &mut diagnostics),
            Ok(true)
        );
        assert_eq!(
            diagnostics,
            vec![existing, expected_needs_hydrogen_warning()]
        );
    }

    #[test]
    fn cf3d_bld_b09_symmetric_triangular_indices_cover_diagonal_and_order() {
        // Fixed cells from the pinned four-atom upper triangle:
        // row 0 => 0..3, row 1 => 4..6, row 2 => 7..8, row 3 => 9.
        for (i, j, expected) in [
            (0, 0, 0),
            (0, 1, 1),
            (1, 2, 5),
            (2, 1, 5),
            (2, 2, 7),
            (3, 0, 3),
            (3, 3, 9),
        ] {
            assert_eq!(two_bit_cell_pos(4, i, j), Ok(expected));
        }
    }

    #[test]
    fn cf3d_bld_b09_validates_empty_single_and_source_allocation_boundaries() {
        assert_eq!(
            two_bit_cell_pos(0, 0, 0),
            Err(UffBuilderError::NeighborMatrixIndexOutOfRange {
                n_atoms: 0,
                i: 0,
                j: 0,
            })
        );
        assert_eq!(two_bit_cell_pos(1, 0, 0), Ok(0));
        assert_eq!(
            two_bit_cell_pos(4, 0, 4),
            Err(UffBuilderError::NeighborMatrixIndexOutOfRange {
                n_atoms: 4,
                i: 0,
                j: 4,
            })
        );

        // The last cell at n=65,535 remains representable by the pinned
        // unsigned allocation formula and the helper's unsigned result.
        assert_eq!(two_bit_cell_pos(65_535, 65_534, 65_534), Ok(2_147_450_879));
        // The source allocation numerator n*(n+1) first overflows at 65,536.
        assert_eq!(
            two_bit_cell_pos(65_536, 0, 0),
            Err(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms: 65_536 })
        );
        assert_eq!(
            two_bit_cell_pos(usize::MAX, 0, 0),
            Err(UffBuilderError::NeighborMatrixIndexOverflow {
                n_atoms: usize::MAX
            })
        );
    }

    #[test]
    fn cf3d_bld_b10_replaces_every_cell_value_and_preserves_neighbors() {
        const SHIFTS: [u32; 4] = [0, 2, 4, 6];
        const OTHER_BITS: [u8; 4] = [0xe4, 0x93, 0x4e, 0x21];
        const EXPECTED: [[u8; 4]; 4] = [
            [0xe4, 0xe5, 0xe6, 0xe7],
            [0x93, 0x97, 0x9b, 0x9f],
            [0x4e, 0x5e, 0x6e, 0x7e],
            [0x21, 0x61, 0xa1, 0xe1],
        ];

        for slot in 0..4 {
            for prior in 0..4 {
                let shift = SHIFTS[slot];
                let before = OTHER_BITS[slot] | ((prior as u8) << shift);
                assert_eq!((before >> shift) & 0b11, prior as u8);

                for relation in 0..4 {
                    let mut storage = [before];
                    assert_eq!(set_two_bit_cell(&mut storage, slot, relation as u8), Ok(()));
                    assert_eq!(storage, [EXPECTED[slot][relation]]);
                    assert_eq!(storage[0] & !(0b11 << shift), OTHER_BITS[slot]);
                }
            }
        }
    }

    #[test]
    fn cf3d_bld_b10_rejects_source_and_storage_positions_without_mutation() {
        let mut empty = [];
        assert_eq!(
            set_two_bit_cell(&mut empty, 0, 1),
            Err(UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: 0,
                byte_index: 0,
                storage_len: 0,
            })
        );

        let mut one_byte = [0x96];
        assert_eq!(
            set_two_bit_cell(&mut one_byte, 4, 2),
            Err(UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: 4,
                byte_index: 1,
                storage_len: 1,
            })
        );
        assert_eq!(one_byte, [0x96]);

        if usize::MAX > u32::MAX as usize {
            assert_eq!(
                set_two_bit_cell(&mut one_byte, usize::MAX, 3),
                Err(UffBuilderError::NeighborMatrixCellPositionOverflow {
                    position: usize::MAX,
                })
            );
            assert_eq!(one_byte, [0x96]);
        }
    }

    #[test]
    fn cf3d_bld_b11_reads_all_relations_at_each_shift_and_storage_edge() {
        const PACKED: [u8; 4] = [0xe4, 0x1b, 0xb1, 0x4e];
        const EXPECTED: [[u8; 4]; 4] = [[0, 1, 2, 3], [3, 2, 1, 0], [1, 0, 3, 2], [2, 3, 0, 1]];

        for byte_index in 0..4 {
            for slot in 0..4 {
                let position = byte_index * 4 + slot;
                assert_eq!(
                    get_two_bit_cell(&PACKED, position),
                    Ok(EXPECTED[byte_index][slot])
                );
            }
        }
        assert_eq!(get_two_bit_cell(&PACKED, 0), Ok(0));
        assert_eq!(get_two_bit_cell(&PACKED, 15), Ok(1));
    }

    #[test]
    fn cf3d_bld_b11_returns_typed_errors_for_invalid_storage_positions() {
        let empty = [];
        assert_eq!(
            get_two_bit_cell(&empty, 0),
            Err(UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: 0,
                byte_index: 0,
                storage_len: 0,
            })
        );

        let one_byte = [0xe4];
        assert_eq!(
            get_two_bit_cell(&one_byte, 4),
            Err(UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: 4,
                byte_index: 1,
                storage_len: 1,
            })
        );
        if usize::MAX > u32::MAX as usize {
            assert_eq!(
                get_two_bit_cell(&one_byte, usize::MAX),
                Err(UffBuilderError::NeighborMatrixCellPositionOverflow {
                    position: usize::MAX,
                })
            );
        }
    }

    #[test]
    fn cf3d_bld_b12_preserves_empty_unsigned_wrap_and_single_diagonal() {
        let empty = build_neighbor_matrix(&topology_with_bonds(0, &[]))
            .expect("the pinned zero-atom unsigned allocation path is reproducible");
        assert_eq!(empty.len(), 536_870_912);
        assert!(empty.iter().all(|&byte| byte == 0xff));

        let single = build_neighbor_matrix(&topology_with_bonds(1, &[]))
            .expect("one atom allocates one initialized byte");
        assert_eq!(single.as_slice(), &[0xff]);
        assert_eq!(
            get_two_bit_cell(&single, two_bit_cell_pos(1, 0, 0).unwrap()),
            Ok(3)
        );
    }

    #[test]
    fn cf3d_bld_b12_matches_fixed_chain_disconnected_and_square_bytes() {
        let chain = topology_with_bonds(4, &[(0, 1, false), (1, 2, false), (2, 3, false)]);
        assert_eq!(
            build_neighbor_matrix(&chain).unwrap().as_slice(),
            &[0xd3, 0xd3, 0xfc]
        );

        let disconnected = topology_with_bonds(5, &[(0, 1, false), (1, 2, false), (3, 4, false)]);
        assert_eq!(
            build_neighbor_matrix(&disconnected).unwrap().as_slice(),
            &[0xd3, 0xcf, 0xff, 0xf3]
        );

        let square = topology_with_bonds(
            4,
            &[(0, 1, false), (1, 2, false), (2, 3, false), (3, 0, false)],
        );
        assert_eq!(
            build_neighbor_matrix(&square).unwrap().as_slice(),
            &[0x13, 0xd3, 0xfc]
        );
    }

    #[test]
    fn cf3d_bld_b12_cycle_overwrites_depend_on_source_bond_order_and_orientation() {
        let first_order = topology_with_bonds(3, &[(0, 1, false), (1, 2, false), (0, 2, false)]);
        assert_eq!(
            build_neighbor_matrix(&first_order).unwrap().as_slice(),
            &[0xc7, 0xfc]
        );

        let reordered_reversed =
            topology_with_bonds(3, &[(2, 0, false), (2, 1, false), (1, 0, false)]);
        assert_eq!(
            build_neighbor_matrix(&reordered_reversed)
                .unwrap()
                .as_slice(),
            &[0xd3, 0xfc]
        );
    }

    #[test]
    fn cf3d_bld_b12_covers_all_four_shared_endpoint_branches_and_padding() {
        let topology = topology_with_bonds(
            6,
            &[
                (0, 1, false),
                (0, 2, false),
                (3, 0, false),
                (1, 4, false),
                (5, 1, false),
            ],
        );
        let matrix = build_neighbor_matrix(&topology).unwrap();
        assert_eq!(matrix.as_slice(), &[0x03, 0x75, 0xc1, 0xfd, 0x7f, 0xff]);

        for atom in 0..6 {
            let diagonal = two_bit_cell_pos(6, atom, atom).unwrap();
            assert_eq!(get_two_bit_cell(&matrix, diagonal), Ok(3));
        }
        for padding_cell in 21..24 {
            assert_eq!(get_two_bit_cell(&matrix, padding_cell), Ok(3));
        }
    }

    #[test]
    fn cf3d_bld_b12_rejects_positive_source_allocation_overflow_before_reserving() {
        let too_many_atoms = topology_with_bonds(65_536, &[]);
        assert_eq!(
            build_neighbor_matrix(&too_many_atoms),
            Err(UffBuilderError::NeighborMatrixIndexOverflow { n_atoms: 65_536 })
        );
    }

    #[test]
    fn cf3d_bld_b13_crosses_endpoint_parameter_masks_with_source_positive_orders() {
        // Fixed values from Builder.cpp + BondStretch.cpp for r1=.75, equal
        // GMP_Xi=1, Z1=1 and a fixed 2.5-unit endpoint separation.
        const SOURCE_ORDERS: [(BondOrder, f64, f64); 14] = [
            (BondOrder::Single, 98.38814814814815, 196.7762962962963),
            (BondOrder::Double, 170.53514119351411, 299.58105975653194),
            (BondOrder::Triple, 235.20478661078965, 385.73884224635083),
            (BondOrder::Quadruple, 295.99671633644476, 463.5880630411399),
            (BondOrder::Quintuple, 354.3871345071998, 536.3140642552067),
            (BondOrder::Hextuple, 411.155679024958, 605.5340711473896),
            (BondOrder::OneAndHalf, 135.812987407545, 251.270099458445),
            (
                BondOrder::TwoAndHalf,
                203.49690045132726,
                344.01355734300944,
            ),
            (
                BondOrder::ThreeAndHalf,
                265.96914818914445,
                425.4478349151615,
            ),
            (BondOrder::FourAndHalf, 325.4334882526594, 500.4687717749772),
            (
                BondOrder::FiveAndHalf,
                382.93976274503427,
                571.2923862851911,
            ),
            (BondOrder::Aromatic, 135.812987407545, 251.270099458445),
            (BondOrder::Dative, 98.38814814814815, 196.7762962962963),
            (BondOrder::DativeOne, 98.38814814814815, 196.7762962962963),
        ];
        let masks = [(false, false), (false, true), (true, false), (true, true)];

        for (order, expected_energy, expected_gradient) in SOURCE_ORDERS {
            for (has_begin, has_end) in masks {
                let topology = topology_with_bond_orders(2, &[(0, 1, order)]);
                let params = [
                    has_begin.then_some(atomic_params(0.75, 1.0)),
                    has_end.then_some(atomic_params(0.75, 1.0)),
                ];
                let borrowed_params = cf3d_bld_borrowed_params(&params);
                let mut stored_rows = [[0.0; 3]; 2];
                let mut field = ForceField::new(3);
                attach_positions(&mut field, &mut stored_rows);

                assert_eq!(add_bonds(&topology, &borrowed_params, &mut field), Ok(()));
                field.initialize().expect("two fixed points initialize");
                let coordinates = [0.0, 0.0, 0.0, 2.5, 0.0, 0.0];
                let energy = cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                    .expect("fixed bond geometry evaluates");
                let mut gradient = [0.0; 6];
                cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
                    .expect("fixed bond gradient evaluates");

                if has_begin && has_end {
                    assert!(
                        (energy - expected_energy).abs() < 1.0e-10,
                        "wrong energy for {order:?}: {energy}"
                    );
                    let expected = [-expected_gradient, 0.0, 0.0, expected_gradient, 0.0, 0.0];
                    for (actual, expected) in gradient.into_iter().zip(expected) {
                        assert!(
                            (actual - expected).abs() < 1.0e-10,
                            "wrong gradient for {order:?}: {gradient:?}"
                        );
                    }
                } else {
                    assert_eq!(energy, 0.0, "skipped {order:?} term emitted energy");
                    assert_eq!(
                        gradient, [0.0; 6],
                        "skipped {order:?} term emitted gradient"
                    );
                }
            }
        }
    }

    #[test]
    fn cf3d_bld_b13_preserves_bond_row_contribution_order_and_real_evaluation() {
        let topology = topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        assert_eq!(
            topology
                .bonds
                .iter()
                .map(|bond| bond.id().index())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        let large_z = 2.0_f64.powi(54) / (2.0 * 332.06);
        let small_z = 1.0 / (2.0 * 332.06);
        let params = [
            Some(atomic_params(0.5, 1.0)),
            Some(atomic_params(0.5, large_z)),
            Some(atomic_params(0.5, small_z)),
            Some(atomic_params(0.5, small_z)),
        ];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut stored_rows = [[0.0; 3]; 4];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut stored_rows);

        add_bonds(&topology, &borrowed_params, &mut field)
            .expect("all three parameterized rows append");
        field.initialize().expect("four fixed points initialize");
        let coordinates = [0.0, 0.0, 0.0, -2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates),
            Ok(2.0_f64.powi(53))
        );
        let mut gradient = [0.0; 12];
        cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
            .expect("fixed source-order gradient evaluates");
        assert_eq!(
            gradient,
            [
                2.0_f64.powi(54),
                0.0,
                0.0,
                -2.0_f64.powi(54),
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
            ]
        );
    }

    #[test]
    fn cf3d_bld_b13_keeps_typed_skip_constructor_and_first_error_behavior() {
        let missing_parameter_topology =
            topology_with_bond_orders(2, &[(0, 1, BondOrder::DativeLeft)]);
        for params in [
            [None, None],
            [None, Some(atomic_params(0.75, 1.0))],
            [Some(atomic_params(0.75, 1.0)), None],
        ] {
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            let mut stored_rows = [[0.0; 3]; 2];
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut stored_rows);
            assert_eq!(
                add_bonds(&missing_parameter_topology, &borrowed_params, &mut field),
                Ok(())
            );
            field.initialize().expect("two fixed points initialize");
            assert_eq!(cf3d_bld_b05_calc_energy(&mut field, &[0.0; 6]), Ok(0.0));
        }

        for order in [
            BondOrder::Unspecified,
            BondOrder::Ionic,
            BondOrder::Zero,
            BondOrder::Hydrogen,
        ] {
            let topology = topology_with_bond_orders(2, &[(0, 1, order)]);
            let params = [
                Some(atomic_params(0.75, 1.0)),
                Some(atomic_params(0.75, 1.0)),
            ];
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            let mut stored_rows = [[0.0; 3]; 2];
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut stored_rows);
            assert_eq!(
                add_bonds(&topology, &borrowed_params, &mut field),
                Err(UffBuilderError::ForceFieldKernel(
                    ForceFieldKernelError::BadBondOrder
                )),
                "zero-valued source mapping must reach the constructor for {order:?}"
            );
        }

        for order in [
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::ThreeCenter,
            BondOrder::Other,
        ] {
            let topology = topology_with_bond_orders(2, &[(0, 1, order)]);
            let params = [
                Some(atomic_params(0.75, 1.0)),
                Some(atomic_params(0.75, 1.0)),
            ];
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            let mut stored_rows = [[0.0; 3]; 2];
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut stored_rows);
            assert_eq!(
                add_bonds(&topology, &borrowed_params, &mut field),
                Err(UffBuilderError::Valence(ValenceError::BadBondType {
                    bond: None,
                    order,
                }))
            );
        }

        let length_topology = topology_with_bond_orders(2, &[(0, 1, BondOrder::Single)]);
        let mut no_rows = [];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut no_rows);
        assert_eq!(
            add_bonds(&length_topology, &[], &mut field),
            Err(UffBuilderError::ParamsLengthMismatch {
                atoms: 2,
                params: 0,
            })
        );

        for (row_count, argument, index, upper_bound) in [
            (0, BondIndexArgument::First, 0, 0),
            (1, BondIndexArgument::Second, 1, 1),
        ] {
            let mut stored_rows = vec![[0.0; 3]; row_count];
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut stored_rows);
            let params = [
                Some(atomic_params(0.75, 1.0)),
                Some(atomic_params(0.75, 1.0)),
            ];
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            assert_eq!(
                add_bonds(&length_topology, &borrowed_params, &mut field),
                Err(UffBuilderError::ForceFieldKernel(
                    ForceFieldKernelError::BondIndexOutOfRange {
                        argument,
                        index,
                        upper_bound,
                    }
                ))
            );
        }

        let valid_then_invalid =
            topology_with_bond_orders(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Zero)]);
        let params = [
            Some(atomic_params(0.75, 1.0)),
            Some(atomic_params(0.75, 1.0)),
            Some(atomic_params(0.75, 1.0)),
        ];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut stored_rows = [[0.0; 3]; 3];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut stored_rows);
        assert_eq!(
            add_bonds(&valid_then_invalid, &borrowed_params, &mut field),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::BadBondOrder
            ))
        );
        field
            .initialize()
            .expect("first valid term remains evaluable");
        let coordinates = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        assert!(
            (cf3d_bld_b05_calc_energy(&mut field, &coordinates).unwrap() - 24.597037037037037)
                .abs()
                < 1.0e-10
        );
        let mut gradient = [0.0; 9];
        cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient).unwrap();
        assert_eq!(
            gradient,
            [
                -98.38814814814815,
                0.0,
                0.0,
                98.38814814814815,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0
            ]
        );

        let invalid_then_valid =
            topology_with_bond_orders(3, &[(1, 2, BondOrder::Zero), (0, 1, BondOrder::Single)]);
        let mut stored_rows = [[0.0; 3]; 3];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut stored_rows);
        assert_eq!(
            add_bonds(&invalid_then_valid, &borrowed_params, &mut field),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::BadBondOrder
            ))
        );
        field
            .initialize()
            .expect("empty failed prefix remains evaluable");
        assert_eq!(cf3d_bld_b05_calc_energy(&mut field, &coordinates), Ok(0.0));
    }

    fn cf3d_bld_b14_ring_info(atom_count: usize, edges: &[(usize, usize, bool)]) -> RingInfo {
        let topology = topology_with_bonds(atom_count, edges);
        cosmolkit_core::fast_find_rings(&topology).expect("fixed ring topology is valid")
    }

    #[test]
    fn cf3d_bld_b14_all_hybridizations_keep_source_defaults() {
        let rings = cf3d_bld_b14_ring_info(3, &[]);
        let cases = [
            (Hybridization::Unspecified, 0),
            (Hybridization::S, 0),
            (Hybridization::Sp, 1),
            (Hybridization::Sp2, 3),
            (Hybridization::Sp3, 0),
            (Hybridization::Sp2d, 0),
            (Hybridization::Sp3d, 0),
            (Hybridization::Sp3d2, 4),
            (Hybridization::Other, 0),
        ];

        for (hybridization, expected) in cases {
            assert_eq!(
                angle_order(
                    hybridization,
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    &rings,
                ),
                expected,
                "hybridization {hybridization:?}"
            );
        }
    }

    #[test]
    fn cf3d_bld_b14_sp2_uses_three_ring_terminal_membership_matrix() {
        let rings = cf3d_bld_b14_ring_info(
            5,
            &[
                (0, 1, false),
                (1, 2, false),
                (2, 0, false),
                (0, 3, false),
                (0, 4, false),
            ],
        );
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(0), 3));
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(1), 3));
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(2), 3));
        assert!(!rings.is_atom_in_ring_of_size(AtomId::new(3), 3));
        assert!(!rings.is_atom_in_ring_of_size(AtomId::new(4), 3));

        let cases = [
            (AtomId::new(3), AtomId::new(4), 3),
            (AtomId::new(1), AtomId::new(3), 30),
            (AtomId::new(3), AtomId::new(1), 30),
            (AtomId::new(1), AtomId::new(2), 35),
        ];
        for (atom_i, atom_k, expected) in cases {
            assert_eq!(
                angle_order(Hybridization::Sp2, atom_i, AtomId::new(0), atom_k, &rings),
                expected
            );
        }
    }

    #[test]
    fn cf3d_bld_b14_sp2_uses_four_ring_terminal_membership_matrix() {
        let rings = cf3d_bld_b14_ring_info(
            6,
            &[
                (0, 1, false),
                (1, 2, false),
                (2, 3, false),
                (3, 0, false),
                (0, 4, false),
                (0, 5, false),
            ],
        );
        assert!(!rings.is_atom_in_ring_of_size(AtomId::new(0), 3));
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(0), 4));
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(1), 4));
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(3), 4));
        assert!(!rings.is_atom_in_ring_of_size(AtomId::new(4), 4));
        assert!(!rings.is_atom_in_ring_of_size(AtomId::new(5), 4));

        let cases = [
            (AtomId::new(4), AtomId::new(5), 3),
            (AtomId::new(1), AtomId::new(4), 40),
            (AtomId::new(4), AtomId::new(1), 40),
            (AtomId::new(1), AtomId::new(3), 45),
        ];
        for (atom_i, atom_k, expected) in cases {
            assert_eq!(
                angle_order(Hybridization::Sp2, atom_i, AtomId::new(0), atom_k, &rings),
                expected
            );
        }
    }

    #[test]
    fn cf3d_bld_b14_sp2_uses_individual_ring_membership_and_three_ring_priority() {
        let two_triangles = cf3d_bld_b14_ring_info(
            5,
            &[
                (0, 1, false),
                (1, 2, false),
                (2, 0, false),
                (0, 3, false),
                (3, 4, false),
                (4, 0, false),
            ],
        );
        assert!(two_triangles.is_atom_in_ring_of_size(AtomId::new(1), 3));
        assert!(two_triangles.is_atom_in_ring_of_size(AtomId::new(3), 3));
        assert!(!two_triangles.are_atoms_in_same_ring_of_size(AtomId::new(1), AtomId::new(3), 3));
        assert_eq!(
            angle_order(
                Hybridization::Sp2,
                AtomId::new(1),
                AtomId::new(0),
                AtomId::new(3),
                &two_triangles,
            ),
            35
        );

        let both_sizes = cf3d_bld_b14_ring_info(
            6,
            &[
                (0, 1, false),
                (1, 2, false),
                (2, 0, false),
                (0, 3, false),
                (3, 4, false),
                (4, 5, false),
                (5, 0, false),
            ],
        );
        assert!(both_sizes.is_atom_in_ring_of_size(AtomId::new(0), 3));
        assert!(both_sizes.is_atom_in_ring_of_size(AtomId::new(0), 4));
        assert!(both_sizes.is_atom_in_ring_of_size(AtomId::new(3), 4));
        assert!(both_sizes.is_atom_in_ring_of_size(AtomId::new(5), 4));
        assert!(!both_sizes.is_atom_in_ring_of_size(AtomId::new(3), 3));
        assert!(!both_sizes.is_atom_in_ring_of_size(AtomId::new(5), 3));
        assert_eq!(
            angle_order(
                Hybridization::Sp2,
                AtomId::new(3),
                AtomId::new(0),
                AtomId::new(5),
                &both_sizes,
            ),
            3
        );

        for hybridization in [
            Hybridization::Unspecified,
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Other,
        ] {
            let expected = match hybridization {
                Hybridization::Sp => 1,
                Hybridization::Sp3d2 => 4,
                _ => 0,
            };
            assert_eq!(
                angle_order(
                    hybridization,
                    AtomId::new(3),
                    AtomId::new(0),
                    AtomId::new(5),
                    &both_sizes,
                ),
                expected,
                "hybridization {hybridization:?}"
            );
        }
    }

    fn cf3d_bld_b16_neighbors(order: [usize; 5]) -> Vec<NeighborRef> {
        order
            .into_iter()
            .map(|bond| NeighborRef {
                atom_index: bond + 1,
                bond: BondId::new(bond),
            })
            .collect()
    }

    fn cf3d_bld_b16_positions(
        center: [f64; 3],
        directions_by_bond: [[f64; 3]; 5],
    ) -> Vec<[f64; 3]> {
        let mut positions = vec![[0.0; 3]; 6];
        positions[0] = center;
        for (bond, direction) in directions_by_bond.into_iter().enumerate() {
            positions[bond + 1] = [
                center[0] + direction[0],
                center[1] + direction[1],
                center[2] + direction[2],
            ];
        }
        positions
    }

    #[test]
    fn cf3d_bld_b16_selects_each_of_the_ten_source_bond_pairs() {
        let pairs = [
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 3),
            (2, 4),
            (3, 4),
        ];

        for (first, second) in pairs {
            let mut directions = [[0.0, 0.0, 1.0]; 5];
            directions[first] = [1.0, 0.0, 0.0];
            directions[second] = [-1.0, 0.0, 0.0];
            let positions = cf3d_bld_b16_positions([0.0; 3], directions);
            assert_eq!(
                select_tbp_axial(0, &cf3d_bld_b16_neighbors([0, 1, 2, 3, 4]), &positions,),
                Ok((BondId::new(first), BondId::new(second))),
                "fixed antipodal pair ({first}, {second})"
            );
        }
    }

    #[test]
    fn cf3d_bld_b16_equal_minima_retain_first_nested_source_pair() {
        let directions = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ];
        let positions = cf3d_bld_b16_positions([0.0; 3], directions);
        assert_eq!(
            select_tbp_axial(0, &cf3d_bld_b16_neighbors([2, 0, 1, 3, 4]), &positions,),
            Ok((BondId::new(2), BondId::new(4)))
        );
    }

    #[test]
    fn cf3d_bld_b16_uses_bond_ids_across_adjacency_order() {
        let mut directions = [[0.0, 0.0, 1.0]; 5];
        directions[0] = [1.0, 0.0, 0.0];
        directions[2] = [-1.0, 0.0, 0.0];
        let positions = cf3d_bld_b16_positions([0.0; 3], directions);

        assert_eq!(
            select_tbp_axial(0, &cf3d_bld_b16_neighbors([4, 2, 3, 0, 1]), &positions,),
            Ok((BondId::new(0), BondId::new(2)))
        );
    }

    #[test]
    fn cf3d_bld_b16_is_invariant_to_fixed_translation_and_positive_scale() {
        let center = [17.0, -3.5, 91.0];
        let directions = [
            [33.0, 0.0, 0.0],
            [-33.0, 0.0, 0.0],
            [0.0, 22.0, 0.0],
            [0.0, 0.0, 44.0],
            [0.0, 11.0, 11.0],
        ];
        let positions = cf3d_bld_b16_positions(center, directions);

        assert_eq!(
            select_tbp_axial(0, &cf3d_bld_b16_neighbors([0, 1, 2, 3, 4]), &positions,),
            Ok((BondId::new(0), BondId::new(1)))
        );
    }

    #[test]
    fn cf3d_bld_b16_zero_and_sub_tolerance_vectors_return_typed_source_error() {
        for direction in [[0.0, 0.0, 0.0], [1.0e-17, 0.0, 0.0]] {
            let mut directions = [[0.0, 0.0, 1.0]; 5];
            directions[0] = direction;
            let positions = cf3d_bld_b16_positions([0.0; 3], directions);

            assert_eq!(
                select_tbp_axial(0, &cf3d_bld_b16_neighbors([0, 1, 2, 3, 4]), &positions,),
                Err(UffBuilderError::SourceDirectionVectorBelowTolerance {
                    center_atom_index: 0,
                    neighbor_atom_index: 1,
                })
            );
        }
    }

    #[test]
    fn cf3d_bld_b16_nan_vectors_reach_the_source_missing_axis_error() {
        let directions = [[f64::NAN, 0.0, 0.0]; 5];
        let positions = cf3d_bld_b16_positions([0.0; 3], directions);

        assert_eq!(
            select_tbp_axial(0, &cf3d_bld_b16_neighbors([0, 1, 2, 3, 4]), &positions,),
            Err(UffBuilderError::SourceTbpAxialBondNotFound {
                center_atom_index: 0,
            })
        );
    }

    fn cf3d_bld_b17_neighbors(order: &[usize]) -> Vec<NeighborRef> {
        order
            .iter()
            .map(|&bond| NeighborRef {
                atom_index: bond + 1,
                bond: BondId::new(bond),
            })
            .collect()
    }

    #[test]
    fn cf3d_bld_b17_all_axial_pairs_keep_the_three_source_roles() {
        let cases = [
            ((0, 1), [2, 3, 4]),
            ((0, 2), [1, 3, 4]),
            ((0, 3), [1, 2, 4]),
            ((0, 4), [1, 2, 3]),
            ((1, 2), [0, 3, 4]),
            ((1, 3), [0, 2, 4]),
            ((1, 4), [0, 2, 3]),
            ((2, 3), [0, 1, 4]),
            ((2, 4), [0, 1, 3]),
            ((3, 4), [0, 1, 2]),
        ];
        let neighbors = cf3d_bld_b17_neighbors(&[0, 1, 2, 3, 4]);

        for ((axial1, axial2), expected) in cases {
            assert_eq!(
                select_tbp_equatorial(0, &neighbors, (BondId::new(axial1), BondId::new(axial2)),),
                Ok(expected.map(BondId::new)),
                "axial pair ({axial1}, {axial2})"
            );
        }
    }

    #[test]
    fn cf3d_bld_b17_preserves_incident_order_across_bond_id_permutations() {
        assert_eq!(
            select_tbp_equatorial(
                0,
                &cf3d_bld_b17_neighbors(&[4, 2, 0, 3, 1]),
                (BondId::new(3), BondId::new(0)),
            ),
            Ok([BondId::new(4), BondId::new(2), BondId::new(1)])
        );
        assert_eq!(
            select_tbp_equatorial(
                0,
                &cf3d_bld_b17_neighbors(&[3, 1, 4, 0, 2]),
                (BondId::new(4), BondId::new(1)),
            ),
            Ok([BondId::new(3), BondId::new(0), BondId::new(2)])
        );
    }

    #[test]
    fn cf3d_bld_b17_reports_each_missing_source_role_in_order() {
        let cases: [(&[usize], u8); 3] = [(&[0, 1], 1), (&[0, 1, 2], 2), (&[0, 1, 2, 3], 3)];
        for (order, role) in cases {
            assert_eq!(
                select_tbp_equatorial(
                    7,
                    &cf3d_bld_b17_neighbors(order),
                    (BondId::new(0), BondId::new(1)),
                ),
                Err(UffBuilderError::SourceTbpEquatorialBondNotFound {
                    center_atom_index: 7,
                    role,
                })
            );
        }
    }

    #[test]
    fn cf3d_bld_b17_does_not_deduplicate_malformed_incident_roles() {
        assert_eq!(
            select_tbp_equatorial(
                0,
                &cf3d_bld_b17_neighbors(&[0, 1, 2, 2, 3]),
                (BondId::new(0), BondId::new(1)),
            ),
            Ok([BondId::new(2), BondId::new(2), BondId::new(3)])
        );
    }

    struct Cf3dBldB18FixedEnergy(f64);

    impl ForceFieldContribution for Cf3dBldB18FixedEnergy {
        fn get_energy(
            &self,
            _context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            Ok(self.0)
        }

        fn get_grad(
            &self,
            _context: &mut EvaluationContext<'_>,
            _gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(Self(self.0))
        }
    }

    fn cf3d_bld_b18_topology() -> TopologyBlock {
        topology_with_bond_orders(3, &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)])
    }

    fn cf3d_bld_b18_axial_neighbors() -> [NeighborRef; 2] {
        [
            NeighborRef {
                atom_index: 1,
                bond: BondId::new(0),
            },
            NeighborRef {
                atom_index: 2,
                bond: BondId::new(1),
            },
        ]
    }

    fn cf3d_bld_b18_field<'a>(rows: &'a mut [[f64; 3]]) -> ForceField<'a> {
        let mut field = ForceField::new(3);
        attach_positions(&mut field, rows);
        field.initialize().expect("fixed B18 field initializes");
        field
    }

    #[test]
    fn cf3d_bld_b18_all_endpoint_masks_keep_the_source_guard_order() {
        let topology = cf3d_bld_b18_topology();
        let neighbors = cf3d_bld_b18_axial_neighbors();
        let coordinates = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 0.8, 0.0];

        for (first_present, second_present) in [(false, false), (false, true), (true, false)] {
            let mut rows = [[0.0; 3]; 3];
            let mut field = cf3d_bld_b18_field(&mut rows);
            field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
            let endpoint = atomic_params(0.5, 1.0);
            let params = vec![
                None,
                first_present.then_some(endpoint),
                second_present.then_some(endpoint),
            ];
            let borrowed_params = cf3d_bld_borrowed_params(&params);

            assert_eq!(
                append_tbp_axial_angle(&topology, &borrowed_params, 0, neighbors, &mut field),
                Ok(()),
                "endpoint presence ({first_present}, {second_present})"
            );
            assert_eq!(
                cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                    .expect("guarded B18 field evaluates"),
                7.25,
                "missing endpoint must append no angle"
            );
            let mut gradient = [3.0; 9];
            cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
                .expect("guarded B18 gradient evaluates");
            assert_eq!(gradient, [3.0; 9]);
        }
    }

    #[test]
    fn cf3d_bld_b18_missing_center_fails_only_after_both_endpoints_exist() {
        let topology = cf3d_bld_b18_topology();
        let neighbors = cf3d_bld_b18_axial_neighbors();
        let coordinates = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 0.8, 0.0];
        let mut rows = [[0.0; 3]; 3];
        let mut field = cf3d_bld_b18_field(&mut rows);
        field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let endpoint = atomic_params(0.5, 1.0);
        let params = [None, Some(endpoint), Some(endpoint)];
        let borrowed_params = cf3d_bld_borrowed_params(&params);

        assert_eq!(
            append_tbp_axial_angle(&topology, &borrowed_params, 0, neighbors, &mut field),
            Err(UffBuilderError::SourceTbpCenterParamsMissing {
                center_atom_index: 0,
            })
        );
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("failed constructor leaves the existing field evaluable"),
            7.25,
            "constructor failure must append no angle"
        );
    }

    #[test]
    fn cf3d_bld_b18_appends_order_two_angle_with_fixed_kernel_values() {
        let topology = cf3d_bld_b18_topology();
        let neighbors = cf3d_bld_b18_axial_neighbors();
        let coordinates = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 0.8, 0.0];
        let mut rows = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.6, 0.8, 0.0]];
        let mut field = cf3d_bld_b18_field(&mut rows);
        // Source-order summation makes the B18 term observable only when it is
        // appended after both existing terms: 1e20 + -1e20 + E == E.
        field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(1.0e20)));
        field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(-1.0e20)));
        let mut center = atomic_params(0.5, 1.0);
        center.theta0 = std::f64::consts::FRAC_PI_2;
        let endpoint = atomic_params(0.5, 1.0);
        let params = [Some(center), Some(endpoint), Some(endpoint)];
        let borrowed_params = cf3d_bld_borrowed_params(&params);

        append_tbp_axial_angle(&topology, &borrowed_params, 0, neighbors, &mut field)
            .expect("complete B18 endpoint and center parameters append the angle");

        let energy = cf3d_bld_b05_calc_energy(&mut field, &coordinates)
            .expect("real B18 angle energy evaluates");
        assert!(
            (energy - 112.70490132518644).abs() < 1.0e-9,
            "energy {energy}"
        );

        let mut gradient = [0.0; 9];
        cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
            .expect("real B18 angle gradient evaluates");
        let expected_gradient = [
            135.24588159022372,
            67.62294079511186,
            0.0,
            0.0,
            -169.05735198777964,
            0.0,
            -135.24588159022372,
            101.43441119266778,
            0.0,
        ];
        for (axis, (actual, expected)) in gradient.into_iter().zip(expected_gradient).enumerate() {
            assert!(
                (actual - expected).abs() < 1.0e-9,
                "gradient component {axis}: actual {actual}, expected {expected}"
            );
        }
    }

    fn cf3d_bld_b19_topology(third_bond_order: BondOrder) -> TopologyBlock {
        topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, third_bond_order),
            ],
        )
    }

    fn cf3d_bld_b19_neighbors() -> [NeighborRef; 3] {
        [
            NeighborRef {
                atom_index: 1,
                bond: BondId::new(0),
            },
            NeighborRef {
                atom_index: 2,
                bond: BondId::new(1),
            },
            NeighborRef {
                atom_index: 3,
                bond: BondId::new(2),
            },
        ]
    }

    fn cf3d_bld_b19_field<'a>(rows: &'a mut [[f64; 3]]) -> ForceField<'a> {
        let mut field = ForceField::new(3);
        attach_positions(&mut field, rows);
        field.initialize().expect("fixed B19 field initializes");
        field
    }

    fn cf3d_bld_b19_center_params() -> AtomicParams {
        let mut center = atomic_params(0.5, 1.0);
        center.theta0 = std::f64::consts::FRAC_PI_2;
        center
    }

    const CF3D_BLD_B19_COORDINATES: [f64; 12] =
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 0.8, 0.0, -0.8, 0.6, 0.0];

    #[test]
    fn cf3d_bld_b19_all_eight_endpoint_masks_keep_fixed_kernel_values() {
        let topology = cf3d_bld_b19_topology(BondOrder::Single);
        let neighbors = cf3d_bld_b19_neighbors();
        let center = cf3d_bld_b19_center_params();
        let endpoint = atomic_params(0.5, 1.0);
        let expected_energy = [
            7.25,
            7.25,
            7.25,
            83.01273922415311,
            7.25,
            32.60860279816696,
            46.383646293467514,
            147.50498831578759,
        ];
        const PAIR_MASKS: [u8; 3] = [0b011, 0b101, 0b110];
        const PAIR_GRADIENTS: [[f64; 12]; 3] = [
            [
                33.06010438872132,
                16.530052194360664,
                0.0,
                0.0,
                -41.325130485901646,
                0.0,
                -33.06010438872132,
                24.79507829154099,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            [
                65.93236727523407,
                197.79710182570233,
                0.0,
                0.0,
                -109.88727879205683,
                0.0,
                0.0,
                0.0,
                0.0,
                -65.93236727523407,
                -87.90982303364545,
                0.0,
            ],
            [
                23.480187776080516,
                -164.36131443256355,
                0.0,
                0.0,
                0.0,
                0.0,
                -93.92075110432204,
                70.44056332824152,
                0.0,
                70.44056332824152,
                93.92075110432204,
                0.0,
            ],
        ];

        for mask in 0..8_u8 {
            let mut rows = [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.6, 0.8, 0.0],
                [-0.8, 0.6, 0.0],
            ];
            let mut field = cf3d_bld_b19_field(&mut rows);
            field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
            let params = [
                Some(center),
                (mask & 0b001 != 0).then_some(endpoint),
                (mask & 0b010 != 0).then_some(endpoint),
                (mask & 0b100 != 0).then_some(endpoint),
            ];
            let borrowed_params = cf3d_bld_borrowed_params(&params);

            append_tbp_equatorial_angles(&topology, &borrowed_params, 0, neighbors, &mut field)
                .expect("valid B19 endpoint mask follows the source guards");

            let energy = cf3d_bld_b05_calc_energy(&mut field, &CF3D_BLD_B19_COORDINATES)
                .expect("fixed B19 energy evaluates");
            assert!(
                (energy - expected_energy[usize::from(mask)]).abs() < 1.0e-9,
                "mask {mask:03b}: energy {energy}, expected {}",
                expected_energy[usize::from(mask)]
            );

            let mut expected_gradient = [3.0; 12];
            for (pair_mask, pair_gradient) in PAIR_MASKS.into_iter().zip(PAIR_GRADIENTS) {
                if mask & pair_mask == pair_mask {
                    for (expected, delta) in expected_gradient.iter_mut().zip(pair_gradient) {
                        *expected += delta;
                    }
                }
            }
            let mut gradient = [3.0; 12];
            cf3d_bld_b05_calc_grad(&mut field, &CF3D_BLD_B19_COORDINATES, &mut gradient)
                .expect("fixed B19 gradient evaluates");
            for (axis, (actual, expected)) in
                gradient.into_iter().zip(expected_gradient).enumerate()
            {
                assert!(
                    (actual - expected).abs() < 1.0e-9,
                    "mask {mask:03b}, gradient component {axis}: actual {actual}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn cf3d_bld_b19_missing_center_is_skipped_without_pairs_and_fails_on_first_pair() {
        let topology = cf3d_bld_b19_topology(BondOrder::Single);
        let neighbors = cf3d_bld_b19_neighbors();
        let endpoint = atomic_params(0.5, 1.0);

        let mut no_pair_rows = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.6, 0.8, 0.0],
            [-0.8, 0.6, 0.0],
        ];
        let mut no_pair_field = cf3d_bld_b19_field(&mut no_pair_rows);
        no_pair_field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let no_pair_params = [None, Some(endpoint), None, None];
        let no_pair_borrowed_params = cf3d_bld_borrowed_params(&no_pair_params);
        assert_eq!(
            append_tbp_equatorial_angles(
                &topology,
                &no_pair_borrowed_params,
                0,
                neighbors,
                &mut no_pair_field,
            ),
            Ok(()),
            "the absent center is not read while every source endpoint pair is guarded out"
        );
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut no_pair_field, &CF3D_BLD_B19_COORDINATES)
                .expect("no-pair field remains evaluable"),
            7.25
        );

        let mut first_pair_rows = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.6, 0.8, 0.0],
            [-0.8, 0.6, 0.0],
        ];
        let mut first_pair_field = cf3d_bld_b19_field(&mut first_pair_rows);
        first_pair_field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let first_pair_params = [None, Some(endpoint), Some(endpoint), Some(endpoint)];
        let first_pair_borrowed_params = cf3d_bld_borrowed_params(&first_pair_params);
        assert_eq!(
            append_tbp_equatorial_angles(
                &topology,
                &first_pair_borrowed_params,
                0,
                neighbors,
                &mut first_pair_field,
            ),
            Err(UffBuilderError::SourceTbpCenterParamsMissing {
                center_atom_index: 0,
            }),
            "the first active source pair observes the absent center parameter"
        );
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut first_pair_field, &CF3D_BLD_B19_COORDINATES)
                .expect("first-pair failure leaves the prior field evaluable"),
            7.25,
            "the first constructor failure appends no contribution"
        );
    }

    #[test]
    fn cf3d_bld_b19_later_pair_failure_preserves_the_source_order_prefix() {
        let topology = cf3d_bld_b19_topology(BondOrder::Other);
        let neighbors = cf3d_bld_b19_neighbors();
        let params = [
            Some(cf3d_bld_b19_center_params()),
            Some(atomic_params(0.5, 1.0)),
            Some(atomic_params(0.5, 1.0)),
            Some(atomic_params(0.5, 1.0)),
        ];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut rows = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.6, 0.8, 0.0],
            [-0.8, 0.6, 0.0],
        ];
        let mut field = cf3d_bld_b19_field(&mut rows);
        field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));

        assert_eq!(
            append_tbp_equatorial_angles(&topology, &borrowed_params, 0, neighbors, &mut field),
            Err(UffBuilderError::Valence(ValenceError::BadBondType {
                bond: None,
                order: BondOrder::Other,
            })),
            "eq1-eq2 appends before the eq1-eq3 bond-type failure"
        );
        let energy = cf3d_bld_b05_calc_energy(&mut field, &CF3D_BLD_B19_COORDINATES)
            .expect("the source-order prefix remains evaluable after the later failure");
        assert!(
            (energy - 83.01273922415311).abs() < 1.0e-9,
            "retained eq1-eq2 prefix energy is {energy}"
        );

        let mut gradient = [3.0; 12];
        cf3d_bld_b05_calc_grad(&mut field, &CF3D_BLD_B19_COORDINATES, &mut gradient)
            .expect("the retained source-order prefix gradient evaluates");
        let expected_gradient = [
            36.06010438872132,
            19.530052194360664,
            3.0,
            3.0,
            -38.325130485901646,
            3.0,
            -30.06010438872132,
            27.79507829154099,
            3.0,
            3.0,
            3.0,
            3.0,
        ];
        for (axis, (actual, expected)) in gradient.into_iter().zip(expected_gradient).enumerate() {
            assert!(
                (actual - expected).abs() < 1.0e-9,
                "retained prefix gradient component {axis}: actual {actual}, expected {expected}"
            );
        }
    }

    #[test]
    fn cf3d_bld_b19_kernel_sum_distinguishes_eq1_eq2_eq3_append_order() {
        let topology = cf3d_bld_b19_topology(BondOrder::Single);
        let neighbors = cf3d_bld_b19_neighbors();
        let params = [
            Some(cf3d_bld_b19_center_params()),
            Some(atomic_params(0.5, 1_427_199_427_767_656.5)),
            Some(atomic_params(0.5, 0.099_961_021_381_882_05)),
            Some(atomic_params(0.5, 2_763_058_092_158_183.0)),
        ];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut rows = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.6, 0.8, 0.0],
            [-0.8, 0.6, 0.0],
        ];
        let mut field = cf3d_bld_b19_field(&mut rows);

        append_tbp_equatorial_angles(&topology, &borrowed_params, 0, neighbors, &mut field)
            .expect("complete B19 parameters append all three source pairs");
        let energy = cf3d_bld_b05_calc_energy(&mut field, &CF3D_BLD_B19_COORDINATES)
            .expect("source-order B19 sum evaluates");
        let large_term = 1.0e32_f64;
        let expected_source_order = f64::from_bits(large_term.to_bits() + 2);
        assert_eq!(
            energy.to_bits(),
            expected_source_order.to_bits(),
            "eq1-eq2, eq1-eq3, eq2-eq3 summation must retain source order; got {energy:?}"
        );
    }

    const CF3D_BLD_B20_PAIR_ENDPOINT_MASKS: [u8; 6] =
        [0b0_0101, 0b0_1001, 0b1_0001, 0b0_0110, 0b0_1010, 0b1_0010];

    const CF3D_BLD_B20_PAIR_ENERGIES: [f64; 6] = [
        56.246819856718076,
        65.621289832837746,
        771.68707211099547,
        7.4107041844419257,
        627.66103006271692,
        9.8809389125892366,
    ];

    const CF3D_BLD_B20_PAIR_GRADIENTS: [[f64; 18]; 6] = [
        [
            310.44910264415972,
            186.26946158649582,
            248.3592821153278,
            0.0,
            -186.26946158649582,
            -248.3592821153278,
            0.0,
            0.0,
            0.0,
            -310.44910264415972,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            362.19061975151971,
            -289.75249580121579,
            217.31437185091181,
            0.0,
            289.75249580121579,
            -217.31437185091181,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -362.19061975151971,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            478.00661180263717,
            0.0,
            -1434.0198354079121,
            0.0,
            0.0,
            796.6776863377288,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -478.00661180263717,
            0.0,
            637.34214907018304,
        ],
        [
            39.30798320382997,
            -91.718627475603256,
            -52.4106442717733,
            0.0,
            0.0,
            0.0,
            36.284292188150737,
            27.213219141113051,
            100.78970052264096,
            -75.592275391980706,
            64.505408334490212,
            -48.379056250867649,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        [
            -1232.3313723535146,
            0.0,
            1232.3313723535146,
            0.0,
            0.0,
            0.0,
            480.9098038452741,
            360.68235288395539,
            -751.42156850824063,
            0.0,
            0.0,
            0.0,
            751.42156850824063,
            -360.68235288395539,
            -480.9098038452741,
            0.0,
            0.0,
            0.0,
        ],
        [
            122.29150330080437,
            -69.880859029031086,
            52.410644271773307,
            0.0,
            0.0,
            0.0,
            -86.007211112653636,
            -64.505408334490227,
            -100.78970052264097,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -36.284292188150744,
            134.38626736352131,
            48.379056250867663,
        ],
    ];

    const CF3D_BLD_B20_COORDINATES: [f64; 18] = [
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.6, 0.8, 0.0, 0.0, 0.6, 0.8, 0.0, -0.8, 0.6, -0.8, 0.0,
        -0.6,
    ];

    fn cf3d_bld_b20_topology(bad_bond: Option<usize>) -> TopologyBlock {
        let mut orders = [BondOrder::Single; 5];
        if let Some(index) = bad_bond {
            orders[index] = BondOrder::Other;
        }
        topology_with_bond_orders(
            6,
            &[
                (0, 1, orders[0]),
                (0, 2, orders[1]),
                (0, 3, orders[2]),
                (0, 4, orders[3]),
                (0, 5, orders[4]),
            ],
        )
    }

    fn cf3d_bld_b20_neighbors() -> ([NeighborRef; 2], [NeighborRef; 3]) {
        (
            [
                NeighborRef {
                    atom_index: 1,
                    bond: BondId::new(0),
                },
                NeighborRef {
                    atom_index: 2,
                    bond: BondId::new(1),
                },
            ],
            [
                NeighborRef {
                    atom_index: 3,
                    bond: BondId::new(2),
                },
                NeighborRef {
                    atom_index: 4,
                    bond: BondId::new(3),
                },
                NeighborRef {
                    atom_index: 5,
                    bond: BondId::new(4),
                },
            ],
        )
    }

    fn cf3d_bld_b20_params(
        terminal_mask: u8,
        center: Option<AtomicParams>,
    ) -> [Option<AtomicParams>; 6] {
        let endpoints = [
            atomic_params(0.5, 0.8),
            atomic_params(0.5, 1.0),
            atomic_params(0.5, 1.2),
            atomic_params(0.5, 1.4),
            atomic_params(0.5, 1.6),
        ];
        [
            center,
            (terminal_mask & 0b0_0001 != 0).then_some(endpoints[0]),
            (terminal_mask & 0b0_0010 != 0).then_some(endpoints[1]),
            (terminal_mask & 0b0_0100 != 0).then_some(endpoints[2]),
            (terminal_mask & 0b0_1000 != 0).then_some(endpoints[3]),
            (terminal_mask & 0b1_0000 != 0).then_some(endpoints[4]),
        ]
    }

    fn cf3d_bld_b20_rows() -> [[f64; 3]; 6] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-0.6, 0.8, 0.0],
            [0.0, 0.6, 0.8],
            [0.0, -0.8, 0.6],
            [-0.8, 0.0, -0.6],
        ]
    }

    fn cf3d_bld_b20_field<'a>(rows: &'a mut [[f64; 3]]) -> ForceField<'a> {
        let mut field = ForceField::new(3);
        attach_positions(&mut field, rows);
        field.initialize().expect("fixed B20 field initializes");
        field
    }

    fn cf3d_bld_b20_center_params() -> AtomicParams {
        let mut center = atomic_params(0.5, 1.0);
        center.theta0 = 1.2;
        center
    }

    fn cf3d_bld_b20_assert_kernel_outputs(
        field: &mut ForceField<'_>,
        expected_pair_mask: u8,
        context: &str,
    ) {
        let mut expected_energy = 7.25;
        let mut expected_gradient = [3.0; 18];
        for pair_index in 0..6 {
            if expected_pair_mask & (1 << pair_index) == 0 {
                continue;
            }
            expected_energy += CF3D_BLD_B20_PAIR_ENERGIES[pair_index];
            for (expected, delta) in expected_gradient
                .iter_mut()
                .zip(CF3D_BLD_B20_PAIR_GRADIENTS[pair_index])
            {
                *expected += delta;
            }
        }

        let energy = cf3d_bld_b05_calc_energy(field, &CF3D_BLD_B20_COORDINATES)
            .expect("fixed B20 energy evaluates");
        assert!(
            (energy - expected_energy).abs() < 1.0e-8,
            "{context}: energy {energy}, expected {expected_energy}"
        );

        let mut gradient = [3.0; 18];
        cf3d_bld_b05_calc_grad(field, &CF3D_BLD_B20_COORDINATES, &mut gradient)
            .expect("fixed B20 gradient evaluates");
        for (axis, (actual, expected)) in gradient.into_iter().zip(expected_gradient).enumerate() {
            assert!(
                (actual - expected).abs() < 1.0e-8,
                "{context}: gradient component {axis}: actual {actual}, expected {expected}"
            );
        }
    }

    #[test]
    fn cf3d_bld_b20_all_32_terminal_masks_match_fixed_kernel_outputs() {
        let topology = cf3d_bld_b20_topology(None);
        let (axial_neighbors, equatorial_neighbors) = cf3d_bld_b20_neighbors();
        let center = cf3d_bld_b20_center_params();

        for terminal_mask in 0..32_u8 {
            let mut rows = cf3d_bld_b20_rows();
            let mut field = cf3d_bld_b20_field(&mut rows);
            field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
            let params = cf3d_bld_b20_params(terminal_mask, Some(center));
            let borrowed_params = cf3d_bld_borrowed_params(&params);

            append_tbp_mixed_angles(
                &topology,
                &borrowed_params,
                0,
                axial_neighbors,
                equatorial_neighbors,
                &mut field,
            )
            .expect("every B20 terminal presence mask follows its independent pair guards");

            let expected_pair_mask = CF3D_BLD_B20_PAIR_ENDPOINT_MASKS.iter().enumerate().fold(
                0_u8,
                |active_pairs, (pair_index, endpoints)| {
                    if terminal_mask & endpoints == *endpoints {
                        active_pairs | (1 << pair_index)
                    } else {
                        active_pairs
                    }
                },
            );
            cf3d_bld_b20_assert_kernel_outputs(
                &mut field,
                expected_pair_mask,
                &format!("terminal mask {terminal_mask:05b}"),
            );
        }
    }

    #[test]
    fn cf3d_bld_b20_missing_center_is_read_only_for_active_pairs() {
        let topology = cf3d_bld_b20_topology(None);
        let (axial_neighbors, equatorial_neighbors) = cf3d_bld_b20_neighbors();

        let mut skipped_rows = cf3d_bld_b20_rows();
        let mut skipped_field = cf3d_bld_b20_field(&mut skipped_rows);
        skipped_field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let skipped_params = cf3d_bld_b20_params(0b0_0001, None);
        let skipped_borrowed_params = cf3d_bld_borrowed_params(&skipped_params);
        assert_eq!(
            append_tbp_mixed_angles(
                &topology,
                &skipped_borrowed_params,
                0,
                axial_neighbors,
                equatorial_neighbors,
                &mut skipped_field,
            ),
            Ok(()),
            "a missing center is not read when all active source pairs are guarded out"
        );
        cf3d_bld_b20_assert_kernel_outputs(&mut skipped_field, 0, "inactive missing center");

        let mut active_rows = cf3d_bld_b20_rows();
        let mut active_field = cf3d_bld_b20_field(&mut active_rows);
        active_field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let active_params = cf3d_bld_b20_params(0b0_0101, None);
        let active_borrowed_params = cf3d_bld_borrowed_params(&active_params);
        assert_eq!(
            append_tbp_mixed_angles(
                &topology,
                &active_borrowed_params,
                0,
                axial_neighbors,
                equatorial_neighbors,
                &mut active_field,
            ),
            Err(UffBuilderError::SourceTbpCenterParamsMissing {
                center_atom_index: 0,
            }),
            "the first active pair observes its missing center parameter"
        );
        cf3d_bld_b20_assert_kernel_outputs(
            &mut active_field,
            0,
            "first active pair center failure",
        );
    }

    #[test]
    fn cf3d_bld_b20_later_bond_failures_keep_the_exact_source_prefix() {
        let (axial_neighbors, equatorial_neighbors) = cf3d_bld_b20_neighbors();
        let center = cf3d_bld_b20_center_params();
        let cases = [
            (0, 0b0_0101, 0b00_0000, "ax1-eq1"),
            (3, 0b0_1101, 0b00_0001, "ax1-eq2"),
            (4, 0b1_1101, 0b00_0011, "ax1-eq3"),
            (1, 0b0_0110, 0b00_0000, "ax2-eq1"),
            (3, 0b0_1110, 0b00_1000, "ax2-eq2"),
            (4, 0b1_1110, 0b01_1000, "ax2-eq3"),
        ];

        for (bad_bond, terminal_mask, expected_pair_mask, pair_name) in cases {
            let topology = cf3d_bld_b20_topology(Some(bad_bond));
            let params = cf3d_bld_b20_params(terminal_mask, Some(center));
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            let mut rows = cf3d_bld_b20_rows();
            let mut field = cf3d_bld_b20_field(&mut rows);
            field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));

            assert_eq!(
                append_tbp_mixed_angles(
                    &topology,
                    &borrowed_params,
                    0,
                    axial_neighbors,
                    equatorial_neighbors,
                    &mut field,
                ),
                Err(UffBuilderError::Valence(ValenceError::BadBondType {
                    bond: None,
                    order: BondOrder::Other,
                })),
                "the first invalid source bond for {pair_name} stops later pairs"
            );
            cf3d_bld_b20_assert_kernel_outputs(
                &mut field,
                expected_pair_mask,
                &format!("failure at {pair_name}"),
            );
        }
    }

    #[test]
    fn cf3d_bld_b20_first_constructor_failure_appends_nothing() {
        let topology = cf3d_bld_b20_topology(None);
        let (axial_neighbors, mut equatorial_neighbors) = cf3d_bld_b20_neighbors();
        equatorial_neighbors[0] = axial_neighbors[0];
        let params = cf3d_bld_b20_params(0b1_1111, Some(cf3d_bld_b20_center_params()));
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut rows = cf3d_bld_b20_rows();
        let mut field = cf3d_bld_b20_field(&mut rows);
        field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));

        assert_eq!(
            append_tbp_mixed_angles(
                &topology,
                &borrowed_params,
                0,
                axial_neighbors,
                equatorial_neighbors,
                &mut field,
            ),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::AngleDegeneratePoints
            )),
            "the first source constructor rejects a repeated endpoint before append"
        );
        cf3d_bld_b20_assert_kernel_outputs(&mut field, 0, "first constructor failure");
    }

    #[test]
    fn cf3d_bld_b20_axial_major_gradient_sum_distinguishes_pair_order() {
        let topology = cf3d_bld_b20_topology(None);
        let (axial_neighbors, equatorial_neighbors) = cf3d_bld_b20_neighbors();
        let mut params = cf3d_bld_b20_params(0b1_1111, Some(cf3d_bld_b20_center_params()));
        params[1] = Some(atomic_params(0.5, 2.912094357152274e-7));
        params[2] = Some(atomic_params(0.5, 1.0695364405790272e-7));
        params[3] = Some(atomic_params(0.5, 5.058278107424415e-7));
        params[4] = Some(atomic_params(0.5, 0.008774016718120612));
        params[5] = Some(atomic_params(0.5, 1.0332019698466278e-5));
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut rows = cf3d_bld_b20_rows();
        let mut field = cf3d_bld_b20_field(&mut rows);

        append_tbp_mixed_angles(
            &topology,
            &borrowed_params,
            0,
            axial_neighbors,
            equatorial_neighbors,
            &mut field,
        )
        .expect("all six B20 source terms append");
        let mut gradient = [0.0; 18];
        cf3d_bld_b05_calc_grad(&mut field, &CF3D_BLD_B20_COORDINATES, &mut gradient)
            .expect("all six B20 gradients evaluate");

        // Pinned AngleBend.cpp source arithmetic gives this center-x sum in
        // axial-major order; reverse pair order gives bits ...6acc instead.
        assert_eq!(gradient[0].to_bits(), 0x3e19_d667_74af_6a62);
    }

    const CF3D_BLD_B21_SOURCE_TERMS: [(usize, usize, usize, u32); 10] = [
        (1, 0, 5, 2),
        (2, 0, 3, 3),
        (2, 0, 4, 3),
        (3, 0, 4, 3),
        (1, 0, 2, 0),
        (1, 0, 3, 0),
        (1, 0, 4, 0),
        (5, 0, 2, 0),
        (5, 0, 3, 0),
        (5, 0, 4, 0),
    ];

    fn cf3d_bld_b21_topology(hybridization: Hybridization, degree: usize) -> TopologyBlock {
        let hybridizations = [
            hybridization,
            Hybridization::Sp3,
            Hybridization::Sp3,
            Hybridization::Sp3,
            Hybridization::Sp3,
            Hybridization::Sp3,
        ];
        let edges = [
            (0, 1, BondOrder::Single),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (0, 5, BondOrder::Single),
        ];
        cf3d_bld_b15_topology(&hybridizations, &edges[..degree])
    }

    fn cf3d_bld_b21_assert_same_kernel_results(
        actual: &mut ForceField<'_>,
        expected: &mut ForceField<'_>,
        context: &str,
    ) {
        let actual_energy = cf3d_bld_b05_calc_energy(actual, &CF3D_BLD_B20_COORDINATES)
            .expect("fixed B21 actual energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(expected, &CF3D_BLD_B20_COORDINATES)
            .expect("fixed B21 source-listed energy evaluates");
        assert_eq!(
            actual_energy.to_bits(),
            expected_energy.to_bits(),
            "{context}: energy actual {actual_energy:?}, expected {expected_energy:?}"
        );

        let mut actual_gradient = [0.0; 18];
        let mut expected_gradient = [0.0; 18];
        cf3d_bld_b05_calc_grad(actual, &CF3D_BLD_B20_COORDINATES, &mut actual_gradient)
            .expect("fixed B21 actual gradient evaluates");
        cf3d_bld_b05_calc_grad(expected, &CF3D_BLD_B20_COORDINATES, &mut expected_gradient)
            .expect("fixed B21 source-listed gradient evaluates");
        assert_eq!(
            actual_gradient.map(f64::to_bits),
            expected_gradient.map(f64::to_bits),
            "{context}: source-listed term order must match component-wise"
        );
    }

    #[test]
    fn cf3d_bld_b21_all_32_endpoint_masks_match_the_fixed_ten_term_source_list() {
        let topology = cf3d_bld_b21_topology(Hybridization::Sp3d, 5);
        let center = cf3d_bld_b20_center_params();

        for terminal_mask in 0..32_u8 {
            let params = cf3d_bld_b20_params(terminal_mask, Some(center));
            let borrowed_params = cf3d_bld_borrowed_params(&params);
            let mut actual_rows = cf3d_bld_b20_rows();
            let original_rows = actual_rows;
            let mut actual = cf3d_bld_b20_field(&mut actual_rows);
            actual.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
            add_trigonal_bipyramid_angles(&topology, 0, &borrowed_params, &mut actual)
                .expect("valid B21 inputs append the source-guarded terms");

            let expected_terms: Vec<_> = CF3D_BLD_B21_SOURCE_TERMS
                .into_iter()
                .filter(|(first, _, third, _)| params[*first].is_some() && params[*third].is_some())
                .collect();
            if terminal_mask == 0b1_1111 {
                assert_eq!(expected_terms, CF3D_BLD_B21_SOURCE_TERMS);
                assert_eq!(expected_terms.len(), 10);
            }
            let mut expected_rows = cf3d_bld_b20_rows();
            let mut expected = cf3d_bld_b20_field(&mut expected_rows);
            expected.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
            cf3d_bld_b15_append_expected_terms(&mut expected, &params, &expected_terms);

            cf3d_bld_b21_assert_same_kernel_results(
                &mut actual,
                &mut expected,
                &format!("endpoint mask {terminal_mask:05b}"),
            );
            for (row, original) in original_rows.iter().enumerate() {
                assert_eq!(
                    &actual.positions()[row][..],
                    &original[..],
                    "endpoint mask {terminal_mask:05b} must not change coordinate row {row}"
                );
            }
        }
    }

    #[test]
    fn cf3d_bld_b21_equal_axial_minima_keep_the_first_source_pair_and_roles() {
        let topology = cf3d_bld_b21_topology(Hybridization::Sp3d, 5);
        let params = cf3d_bld_b20_params(0b1_1111, Some(cf3d_bld_b20_center_params()));
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let tie_rows = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ];
        let original_rows = tie_rows;
        let tie_coordinates = [
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0,
            0.0,
        ];
        let mut actual_rows = tie_rows;
        let mut actual = cf3d_bld_b20_field(&mut actual_rows);
        let mut expected_rows = tie_rows;
        let mut expected = cf3d_bld_b20_field(&mut expected_rows);
        let expected_tie_terms = [
            (1, 0, 2, 2),
            (3, 0, 4, 3),
            (3, 0, 5, 3),
            (4, 0, 5, 3),
            (1, 0, 3, 0),
            (1, 0, 4, 0),
            (1, 0, 5, 0),
            (2, 0, 3, 0),
            (2, 0, 4, 0),
            (2, 0, 5, 0),
        ];

        add_trigonal_bipyramid_angles(&topology, 0, &borrowed_params, &mut actual)
            .expect("fixed equal-minimum B21 geometry is valid");
        cf3d_bld_b15_append_expected_terms(&mut expected, &params, &expected_tie_terms);
        let mut actual_gradient = [0.0; 18];
        let mut expected_gradient = [0.0; 18];
        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &tie_coordinates)
            .expect("fixed tie actual energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &tie_coordinates)
            .expect("fixed tie source-listed energy evaluates");
        cf3d_bld_b05_calc_grad(&mut actual, &tie_coordinates, &mut actual_gradient)
            .expect("fixed tie actual gradient evaluates");
        cf3d_bld_b05_calc_grad(&mut expected, &tie_coordinates, &mut expected_gradient)
            .expect("fixed tie source-listed gradient evaluates");
        assert_eq!(actual_energy.to_bits(), expected_energy.to_bits());
        assert_eq!(
            actual_gradient.map(f64::to_bits),
            expected_gradient.map(f64::to_bits)
        );
        for (row, original) in original_rows.iter().enumerate() {
            assert_eq!(
                &actual.positions()[row][..],
                &original[..],
                "tie geometry must not change coordinate row {row}"
            );
        }
    }

    #[test]
    fn cf3d_bld_b21_missing_center_fails_only_for_an_active_source_term() {
        let topology = cf3d_bld_b21_topology(Hybridization::Sp3d, 5);
        let mut no_endpoint_rows = cf3d_bld_b20_rows();
        let mut no_endpoint_field = cf3d_bld_b20_field(&mut no_endpoint_rows);
        no_endpoint_field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let no_endpoint_params = cf3d_bld_b20_params(0, None);
        let no_endpoint_borrowed_params = cf3d_bld_borrowed_params(&no_endpoint_params);
        assert_eq!(
            add_trigonal_bipyramid_angles(
                &topology,
                0,
                &no_endpoint_borrowed_params,
                &mut no_endpoint_field
            ),
            Ok(())
        );
        cf3d_bld_b20_assert_kernel_outputs(
            &mut no_endpoint_field,
            0,
            "inactive source terms do not read the missing center",
        );

        let mut active_rows = cf3d_bld_b20_rows();
        let mut active_field = cf3d_bld_b20_field(&mut active_rows);
        active_field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        let active_params = cf3d_bld_b20_params(0b1_0001, None);
        let active_borrowed_params = cf3d_bld_borrowed_params(&active_params);
        assert_eq!(
            add_trigonal_bipyramid_angles(&topology, 0, &active_borrowed_params, &mut active_field),
            Err(UffBuilderError::SourceTbpCenterParamsMissing {
                center_atom_index: 0,
            })
        );
        cf3d_bld_b20_assert_kernel_outputs(
            &mut active_field,
            0,
            "first active axial pair fails before appending without center params",
        );
    }

    #[test]
    fn cf3d_bld_b21_source_preconditions_keep_source_order_and_categories() {
        let valid_topology = cf3d_bld_b21_topology(Hybridization::Sp3d, 5);
        let valid_params = cf3d_bld_b20_params(0b1_1111, Some(cf3d_bld_b20_center_params()));
        let valid_borrowed_params = cf3d_bld_borrowed_params(&valid_params);
        let wrong_hybridization = cf3d_bld_b21_topology(Hybridization::Sp2, 5);
        let wrong_degree = cf3d_bld_b21_topology(Hybridization::Sp3d, 4);
        let short_params = [None; 5];
        let mut rows = cf3d_bld_b20_rows();
        let mut field = cf3d_bld_b20_field(&mut rows);

        assert_eq!(
            add_trigonal_bipyramid_angles(&valid_topology, 6, &valid_borrowed_params, &mut field),
            Err(UffBuilderError::SourceTbpAtomPrecondition {
                center_atom_index: 6,
            })
        );
        assert_eq!(
            add_trigonal_bipyramid_angles(&wrong_hybridization, 0, &short_params, &mut field),
            Err(UffBuilderError::SourceTbpHybridizationPrecondition {
                center_atom_index: 0,
                actual: Hybridization::Sp2,
            })
        );
        assert_eq!(
            add_trigonal_bipyramid_angles(&wrong_degree, 0, &short_params, &mut field),
            Err(UffBuilderError::SourceTbpDegreePrecondition {
                center_atom_index: 0,
                actual_degree: 4,
            })
        );
        assert_eq!(
            add_trigonal_bipyramid_angles(&valid_topology, 0, &short_params, &mut field),
            Err(UffBuilderError::ParamsLengthMismatch {
                atoms: 6,
                params: 5,
            })
        );
    }

    fn cf3d_bld_b22_topology(hybridization: Hybridization, degree: usize) -> TopologyBlock {
        let mut hybridizations = vec![Hybridization::Sp3; degree + 1];
        hybridizations[0] = hybridization;
        let edges = (1..=degree)
            .map(|neighbor| (0, neighbor, BondOrder::Single))
            .collect::<Vec<_>>();
        cf3d_bld_b15_topology(&hybridizations, &edges)
    }

    fn cf3d_bld_b22_rows(atom_count: usize) -> Vec<[f64; 3]> {
        let fixed_rows = cf3d_bld_b20_rows();
        (0..atom_count)
            .map(|atom| fixed_rows.get(atom).copied().unwrap_or([2.0, 2.0, 2.0]))
            .collect()
    }

    fn cf3d_bld_b22_coordinates(rows: &[[f64; 3]]) -> Vec<f64> {
        rows.iter().flat_map(|row| row.iter().copied()).collect()
    }

    fn cf3d_bld_b22_two_center_topology() -> TopologyBlock {
        let mut hybridizations = vec![Hybridization::Sp3; 12];
        hybridizations[0] = Hybridization::Sp3d;
        hybridizations[6] = Hybridization::Sp3d;
        let mut edges = Vec::with_capacity(10);
        for (center, first_neighbor) in [(0, 1), (6, 7)] {
            for neighbor in first_neighbor..first_neighbor + 5 {
                edges.push((center, neighbor, BondOrder::Single));
            }
        }
        cf3d_bld_b15_topology(&hybridizations, &edges)
    }

    fn cf3d_bld_b22_two_center_params() -> Vec<Option<AtomicParams>> {
        let mut params = (0..12)
            .map(|atom| Some(atomic_params(0.5, 0.8 + f64::from(atom as u32) * 0.1)))
            .collect::<Vec<_>>();
        params[0] = Some(cf3d_bld_b20_center_params());
        params[6] = Some(cf3d_bld_b20_center_params());
        params
    }

    fn cf3d_bld_b22_two_center_rows() -> Vec<[f64; 3]> {
        let first_center_rows = cf3d_bld_b20_rows();
        let mut rows = first_center_rows.to_vec();
        rows.extend(
            first_center_rows
                .iter()
                .map(|[x, y, z]| [x + 10.0, y - 3.0, z + 2.0]),
        );
        rows
    }

    fn cf3d_bld_b22_assert_same_kernel_results(
        actual: &mut ForceField<'_>,
        expected: &mut ForceField<'_>,
        coordinates: &[f64],
        context: &str,
    ) {
        let actual_energy = cf3d_bld_b05_calc_energy(actual, coordinates)
            .expect("fixed B22 actual energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(expected, coordinates)
            .expect("fixed B22 source-listed energy evaluates");
        assert_eq!(
            actual_energy.to_bits(),
            expected_energy.to_bits(),
            "{context}: actual energy {actual_energy:?}, expected {expected_energy:?}"
        );

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(actual, coordinates, &mut actual_gradient)
            .expect("fixed B22 actual gradient evaluates");
        cf3d_bld_b05_calc_grad(expected, coordinates, &mut expected_gradient)
            .expect("fixed B22 source-listed gradient evaluates");
        assert_eq!(
            actual_gradient
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            expected_gradient
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            "{context}: source-listed term order must match component-wise"
        );
    }

    #[test]
    fn cf3d_bld_b22_dispatches_exactly_the_sp3d_degree_five_case() {
        let hybridizations = [
            Hybridization::Unspecified,
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Other,
        ];

        for hybridization in hybridizations {
            for degree in 0..=6 {
                let topology = cf3d_bld_b22_topology(hybridization, degree);
                let endpoint = atomic_params(0.5, 1.0);
                let mut params = vec![Some(endpoint); degree + 1];
                params[0] = None;
                let borrowed_params = cf3d_bld_borrowed_params(&params);
                let mut rows = cf3d_bld_b22_rows(degree + 1);
                let original_rows = rows.clone();
                let coordinates = cf3d_bld_b22_coordinates(&rows);
                let mut field = cf3d_bld_b20_field(&mut rows);
                field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));

                let result = add_angle_special_cases(&topology, &borrowed_params, &mut field);
                if hybridization == Hybridization::Sp3d && degree == 5 {
                    assert_eq!(
                        result,
                        Err(UffBuilderError::SourceTbpCenterParamsMissing {
                            center_atom_index: 0,
                        }),
                        "the source-dispatched TBP must reach B21 without a center-parameter skip"
                    );
                } else {
                    assert_eq!(
                        result,
                        Ok(()),
                        "{hybridization:?} degree {degree} does not dispatch"
                    );
                }

                assert_eq!(
                    cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                        .expect("fixed B22 matrix field evaluates")
                        .to_bits(),
                    7.25_f64.to_bits(),
                    "{hybridization:?} degree {degree}: no angle is appended before this result"
                );
                for (row, original) in field.positions().iter().zip(&original_rows) {
                    assert_eq!(
                        &row[..],
                        &original[..],
                        "{hybridization:?} degree {degree}: dispatch preserves coordinates"
                    );
                }
            }
        }
    }

    #[test]
    fn cf3d_bld_b22_checks_parameter_length_before_dispatch() {
        let topology = cf3d_bld_b22_topology(Hybridization::Sp3d, 5);
        let mut rows = cf3d_bld_b20_rows();
        let coordinates = cf3d_bld_b22_coordinates(&rows);
        let mut field = cf3d_bld_b20_field(&mut rows);
        field.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));

        assert_eq!(
            add_angle_special_cases(&topology, &[None; 5], &mut field),
            Err(UffBuilderError::ParamsLengthMismatch {
                atoms: 6,
                params: 5,
            })
        );
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("fixed B22 precondition field evaluates")
                .to_bits(),
            7.25_f64.to_bits(),
            "the failed global precondition appends no special term"
        );
    }

    #[test]
    fn cf3d_bld_b22_orders_multiple_centers_and_retains_the_prior_prefix_on_error() {
        let topology = cf3d_bld_b22_two_center_topology();
        let params = cf3d_bld_b22_two_center_params();
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let rows = cf3d_bld_b22_two_center_rows();
        let coordinates = cf3d_bld_b22_coordinates(&rows);

        let mut actual_rows = rows.clone();
        let mut actual = cf3d_bld_b20_field(&mut actual_rows);
        actual.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        add_angle_special_cases(&topology, &borrowed_params, &mut actual)
            .expect("both source-ordered special centers append their complete terms");

        let mut expected_terms = CF3D_BLD_B21_SOURCE_TERMS.to_vec();
        expected_terms.extend(
            CF3D_BLD_B21_SOURCE_TERMS
                .into_iter()
                .map(|(first, center, third, order)| (first + 6, center + 6, third + 6, order)),
        );
        let mut expected_rows = rows.clone();
        let mut expected = cf3d_bld_b20_field(&mut expected_rows);
        expected.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        cf3d_bld_b15_append_expected_terms(&mut expected, &params, &expected_terms);
        cf3d_bld_b22_assert_same_kernel_results(
            &mut actual,
            &mut expected,
            &coordinates,
            "two special centers append in ascending atom order",
        );

        let mut failing_params = params.clone();
        failing_params[6] = None;
        let failing_borrowed_params = cf3d_bld_borrowed_params(&failing_params);
        let mut actual_prefix_rows = rows.clone();
        let mut actual_prefix = cf3d_bld_b20_field(&mut actual_prefix_rows);
        actual_prefix.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        assert_eq!(
            add_angle_special_cases(&topology, &failing_borrowed_params, &mut actual_prefix),
            Err(UffBuilderError::SourceTbpCenterParamsMissing {
                center_atom_index: 6,
            }),
            "the second center reports its active missing parameter after the first center"
        );

        let mut expected_prefix_rows = rows;
        let mut expected_prefix = cf3d_bld_b20_field(&mut expected_prefix_rows);
        expected_prefix.add_contribution(Box::new(Cf3dBldB18FixedEnergy(7.25)));
        cf3d_bld_b15_append_expected_terms(
            &mut expected_prefix,
            &failing_params,
            &CF3D_BLD_B21_SOURCE_TERMS,
        );
        cf3d_bld_b22_assert_same_kernel_results(
            &mut actual_prefix,
            &mut expected_prefix,
            &coordinates,
            "first center's source prefix remains appended before second-center error",
        );
    }

    fn cf3d_bld_b23_atom_predicate_has_source_shape(
        predicate: &QueryNode<AtomQueryPredicate>,
    ) -> bool {
        let QueryNode::And(terms) = predicate else {
            return false;
        };
        let [recursive, degree] = terms.as_slice() else {
            return false;
        };
        matches!(
            recursive,
            QueryNode::Not(child)
                if matches!(
                    child.as_ref(),
                    QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(_))
                )
        ) && matches!(
            degree,
            QueryNode::Not(child)
                if matches!(
                    child.as_ref(),
                    QueryNode::Predicate(AtomQueryPredicate::ExplicitDegree(1))
                )
        )
    }

    fn cf3d_bld_b23_matched_pairs(topology: &TopologyBlock) -> Vec<(usize, usize)> {
        let compiled = default_torsion_query().expect("the pinned default SMARTS compiles");
        let mut pairs = compiled
            .matches(topology)
            .expect("the compiled default SMARTS matches a valid detached target")
            .into_iter()
            .map(|matched| {
                let [first, second] = matched.atom_mapping.as_slice() else {
                    panic!("the pinned query has exactly two atom mappings")
                };
                (usize::min(*first, *second), usize::max(*first, *second))
            })
            .collect::<Vec<_>>();
        pairs.sort_unstable();
        pairs
    }

    #[test]
    fn cf3d_bld_b23_repeated_calls_share_one_cached_compiled_query() {
        let first = default_torsion_query().expect("the pinned default SMARTS compiles");
        let second = default_torsion_query().expect("the pinned default SMARTS remains cached");

        assert!(std::ptr::eq(first, second));
    }

    #[test]
    fn cf3d_bld_b23_compiles_the_exact_source_query_shape() {
        let compiled = default_torsion_query().expect("the pinned default SMARTS compiles");

        assert_eq!(DEFAULT_TORSION_BOND_SMARTS, "[!$(*#*)&!D1]~[!$(*#*)&!D1]");
        assert_eq!(compiled.num_atoms(), 2);
        assert_eq!(compiled.num_bonds(), 1);
        assert!(
            compiled
                .query()
                .atoms()
                .iter()
                .all(|atom| cf3d_bld_b23_atom_predicate_has_source_shape(atom.predicate()))
        );
        assert_eq!(
            compiled.query().bonds()[0].predicate(),
            &QueryNode::predicate(BondQueryPredicate::Any)
        );
    }

    #[test]
    fn cf3d_bld_b23_matches_fixed_pairs_for_terminals_triples_branches_cycles_and_fragments() {
        let chain = topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        assert_eq!(cf3d_bld_b23_matched_pairs(&chain), vec![(1, 2)]);

        let triple_neighbors = topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Triple),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        assert!(cf3d_bld_b23_matched_pairs(&triple_neighbors).is_empty());

        let branch = topology_with_bond_orders(
            6,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (4, 5, BondOrder::Single),
            ],
        );
        assert_eq!(cf3d_bld_b23_matched_pairs(&branch), vec![(1, 4)]);

        let cycle = topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 0, BondOrder::Single),
            ],
        );
        assert_eq!(
            cf3d_bld_b23_matched_pairs(&cycle),
            vec![(0, 1), (0, 3), (1, 2), (2, 3)]
        );

        let disconnected = topology_with_bond_orders(
            8,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (4, 5, BondOrder::Single),
                (5, 6, BondOrder::Single),
                (6, 7, BondOrder::Single),
            ],
        );
        assert_eq!(
            cf3d_bld_b23_matched_pairs(&disconnected),
            vec![(1, 2), (5, 6)]
        );
    }

    fn cf3d_bld_b24_matches(
        topology: &TopologyBlock,
        smarts: &str,
    ) -> Result<Vec<Vec<usize>>, TorsionBondQueryError> {
        let rings = RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            topology.atoms.len(),
            topology.bonds.len(),
        );
        let zero_rows = vec![0; topology.atoms.len()];
        let valence = assignment(&zero_rows, &zero_rows);
        torsion_bond_matches(topology, &rings, &valence, smarts)
    }

    fn cf3d_bld_b24_atom_mappings(matches: &[Vec<usize>]) -> Vec<Vec<usize>> {
        matches.to_vec()
    }

    #[test]
    fn cf3d_bld_b24_default_and_equivalent_custom_query_keep_identical_results() {
        let chain = topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );

        let default = cf3d_bld_b24_matches(&chain, DEFAULT_TORSION_BOND_SMARTS)
            .expect("the pinned default query matches the fixed chain");
        let custom = cf3d_bld_b24_matches(&chain, "[!D1&!$(*#*)]~[!D1&!$(*#*)]")
            .expect("the equivalent custom query compiles and matches");

        assert_eq!(default, custom);
        assert_eq!(cf3d_bld_b24_atom_mappings(&default), vec![vec![1, 2]]);
    }

    #[test]
    fn cf3d_bld_b24_two_atom_custom_pattern_preserves_query_mapping_order() {
        let chain =
            topology_with_bond_orders(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);

        let matches = cf3d_bld_b24_matches(&chain, "*-*")
            .expect("the two-atom custom query matches both single bonds");

        assert_eq!(
            cf3d_bld_b24_atom_mappings(&matches),
            vec![vec![0, 1], vec![1, 2]]
        );
    }

    #[test]
    fn cf3d_bld_b24_invalid_custom_syntax_retains_typed_parse_error() {
        let chain = topology_with_bond_orders(2, &[(0, 1, BondOrder::Single)]);

        assert!(matches!(
            cf3d_bld_b24_matches(&chain, "["),
            Err(TorsionBondQueryError::Parse(_))
        ));
    }

    #[test]
    fn cf3d_bld_b24_wrong_match_arity_is_the_source_ordered_typed_failure() {
        let chain = topology_with_bond_orders(2, &[(0, 1, BondOrder::Single)]);

        assert_eq!(
            cf3d_bld_b24_matches(&chain, "*").unwrap_err(),
            TorsionBondQueryError::MatchArity {
                smarts: "*".to_owned(),
                actual: 1,
            }
        );
    }

    #[test]
    fn cf3d_bld_b24_disconnected_match_retains_the_source_bond_invariant() {
        let disconnected = topology_with_bond_orders(2, &[]);
        let matches = cf3d_bld_b24_matches(&disconnected, "*.*")
            .expect("the disconnected two-atom query has a valid atom mapping");

        assert_eq!(matches.len(), 1);
        let [begin_atom_index, end_atom_index] = matches[0].as_slice() else {
            panic!("the source match is checked to contain two query atoms")
        };
        assert_eq!(
            source_torsion_bond_index(&disconnected, *begin_atom_index, *end_atom_index,),
            Err(TorsionBondQueryError::MissingMatchedBond {
                begin_atom_index: *begin_atom_index,
                end_atom_index: *end_atom_index,
            })
        );
    }

    #[test]
    fn cf3d_bld_b24_unique_symmetric_cycle_matches_keep_source_orientation_and_order() {
        let cycle = topology_with_bond_orders(
            4,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 0, BondOrder::Single),
            ],
        );

        let matches = cf3d_bld_b24_matches(&cycle, "*~*")
            .expect("the symmetric two-atom query matches the four-cycle");

        assert_eq!(
            cf3d_bld_b24_atom_mappings(&matches),
            vec![vec![0, 1], vec![0, 3], vec![1, 2], vec![2, 3]]
        );
    }

    #[test]
    fn cf3d_bld_b24_custom_recursive_predicate_uses_full_query_matching() {
        let triple_neighbors =
            topology_with_bond_orders(3, &[(0, 1, BondOrder::Triple), (1, 2, BondOrder::Single)]);

        let matches = cf3d_bld_b24_matches(&triple_neighbors, "[$(*#*)]~[*]")
            .expect("the custom recursive predicate matches the triple-bond atom");

        assert_eq!(
            cf3d_bld_b24_atom_mappings(&matches),
            vec![vec![0, 1], vec![1, 2]]
        );
    }

    #[test]
    fn uff_one_fix_custom_query_with_no_hits_returns_empty_atom_rows() {
        let chain =
            topology_with_bond_orders(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);

        assert_eq!(cf3d_bld_b24_matches(&chain, "[O]~[O]"), Ok(Vec::new()));
    }

    #[test]
    fn cf3d_bld_b24_match_results_stop_at_the_source_default_limit() {
        let atom_count = 1_002;
        let edges = (0..atom_count)
            .map(|atom| (atom, (atom + 1) % atom_count, BondOrder::Single))
            .collect::<Vec<_>>();
        let cycle = topology_with_bond_orders(atom_count, &edges);

        let matches = cf3d_bld_b24_matches(&cycle, "*~*")
            .expect("the fixed high-symmetry cycle matches within the source cap");

        assert_eq!(matches.len(), 1_000);
    }

    #[test]
    fn uff_one_fix_default_custom_atom_rows_and_constructor_order() {
        cf3d_bld_b24_default_and_equivalent_custom_query_keep_identical_results();
        cf3d_bld_b26_default_and_custom_matches_append_overlapping_terms_in_source_order();
    }

    #[test]
    fn uff_one_fix_reversed_and_symmetric_atom_rows() {
        cf3d_bld_b24_unique_symmetric_cycle_matches_keep_source_orientation_and_order();
        cf3d_bld_b26_custom_query_preserves_reversed_endpoint_and_nested_order();
    }

    #[test]
    fn uff_one_fix_nonbonded_pairs_and_typed_error_order() {
        cf3d_bld_b24_disconnected_match_retains_the_source_bond_invariant();
        cf3d_bld_b26_preserves_precondition_query_and_bond_error_order();
    }

    #[test]
    fn uff_one_fix_recursive_rows_and_default_match_limit() {
        cf3d_bld_b24_custom_recursive_predicate_uses_full_query_matching();
        cf3d_bld_b24_match_results_stop_at_the_source_default_limit();
    }

    fn cf3d_bld_b25_params(atom_count: usize) -> Vec<Option<AtomicParams>> {
        (0..atom_count)
            .map(|_| {
                let mut params = cf3d_bld_b15_atomic_params();
                params.v1 = 1.0;
                params.u1 = 1.0;
                Some(params)
            })
            .collect()
    }

    fn cf3d_bld_b25_four_points() -> [[f64; 3]; 4] {
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    }

    fn cf3d_bld_b25_append_expected_terms(
        topology: &TopologyBlock,
        params: &[Option<AtomicParams>],
        field: &mut ForceField<'_>,
        terms: &[(usize, usize, usize, usize, bool)],
    ) {
        let scale = u32::try_from(terms.len()).expect("the fixed source match count fits u32");
        for &(first, begin, end, last, terminal_sp2) in terms {
            let mut contribution = {
                TorsionAngleContrib::new(
                    field.positions(),
                    u32::try_from(first).expect("fixed atom index fits u32"),
                    u32::try_from(begin).expect("fixed atom index fits u32"),
                    u32::try_from(end).expect("fixed atom index fits u32"),
                    u32::try_from(last).expect("fixed atom index fits u32"),
                    1.0,
                    i32::from(topology.atoms[begin].atomic_number()),
                    i32::from(topology.atoms[end].atomic_number()),
                    topology.atoms[begin].hybridization(),
                    topology.atoms[end].hybridization(),
                    params[begin]
                        .as_ref()
                        .expect("source-fixed begin-center parameters are present"),
                    params[end]
                        .as_ref()
                        .expect("source-fixed end-center parameters are present"),
                    terminal_sp2,
                )
            }
            .expect("source-fixed expected torsion is valid");
            contribution.scale_force_constant(scale);
            field.add_contribution(Box::new(contribution));
        }
    }

    fn cf3d_bld_b25_assert_kernel_terms(
        topology: &TopologyBlock,
        params: &[Option<AtomicParams>],
        begin: usize,
        end: usize,
        points: &[[f64; 3]],
        expected_terms: &[(usize, usize, usize, usize, bool)],
    ) -> f64 {
        let coordinates = points
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();

        let mut actual_points = points.to_vec();
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);
        let borrowed_params = cf3d_bld_borrowed_params(params);
        torsions_for_bond(topology, &borrowed_params, begin, end, &mut actual)
            .expect("the source-fixed matched bond expands successfully");
        actual
            .initialize()
            .expect("actual source field initializes");

        let mut expected_points = points.to_vec();
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        cf3d_bld_b25_append_expected_terms(topology, params, &mut expected, expected_terms);
        expected
            .initialize()
            .expect("source-fixed expected field initializes");

        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &coordinates)
            .expect("actual source field energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &coordinates)
            .expect("source-fixed expected energy evaluates");
        assert_eq!(actual_energy, expected_energy);

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &coordinates, &mut actual_gradient)
            .expect("actual source field gradient evaluates");
        cf3d_bld_b05_calc_grad(&mut expected, &coordinates, &mut expected_gradient)
            .expect("source-fixed expected gradient evaluates");
        assert_eq!(actual_gradient, expected_gradient);
        actual_energy
    }

    #[test]
    fn cf3d_bld_b25_all_center_hybridization_pairs_follow_source_gate() {
        let hybridizations = [
            Hybridization::Unspecified,
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Other,
        ];
        let points = cf3d_bld_b25_four_points();

        for &begin_hybridization in &hybridizations {
            for &end_hybridization in &hybridizations {
                let topology = cf3d_bld_b15_topology(
                    &[
                        Hybridization::Sp2,
                        begin_hybridization,
                        end_hybridization,
                        Hybridization::Sp3,
                    ],
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Single),
                        (2, 3, BondOrder::Single),
                    ],
                );
                let mut params = cf3d_bld_b25_params(4);
                params[0] = None;
                params[3] = None;
                let source_accepts_centers =
                    matches!(begin_hybridization, Hybridization::Sp2 | Hybridization::Sp3)
                        && matches!(end_hybridization, Hybridization::Sp2 | Hybridization::Sp3);
                let expected_terms = if source_accepts_centers {
                    vec![(0, 1, 2, 3, true)]
                } else {
                    Vec::new()
                };

                cf3d_bld_b25_assert_kernel_terms(
                    &topology,
                    &params,
                    1,
                    2,
                    &points,
                    &expected_terms,
                );
            }
        }
    }

    #[test]
    fn cf3d_bld_b25_only_both_central_parameters_guard_the_match() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let points = cf3d_bld_b25_four_points();

        for (begin_present, end_present) in
            [(false, false), (false, true), (true, false), (true, true)]
        {
            let mut params = cf3d_bld_b25_params(4);
            params[0] = None;
            params[3] = None;
            if !begin_present {
                params[1] = None;
            }
            if !end_present {
                params[2] = None;
            }
            let expected_terms = if begin_present && end_present {
                vec![(0, 1, 2, 3, false)]
            } else {
                Vec::new()
            };

            cf3d_bld_b25_assert_kernel_terms(&topology, &params, 1, 2, &points, &expected_terms);
        }

        let disconnected = topology_with_bond_orders(2, &[]);
        let mut points = [[0.0; 3]; 2];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut points);
        let mut params = cf3d_bld_b25_params(2);
        params[0] = None;
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        assert_eq!(
            torsions_for_bond(&disconnected, &borrowed_params, 0, 1, &mut field),
            Ok(())
        );
        let present_params = cf3d_bld_b25_params(2);
        let present_borrowed_params = cf3d_bld_borrowed_params(&present_params);
        assert_eq!(
            torsions_for_bond(&disconnected, &present_borrowed_params, 0, 1, &mut field),
            Err(UffBuilderError::TorsionBondQuery(
                TorsionBondQueryError::MissingMatchedBond {
                    begin_atom_index: 0,
                    end_atom_index: 1,
                }
            ))
        );
    }

    #[test]
    fn cf3d_bld_b25_terminal_sp2_flag_checks_each_end_and_both_center_orders() {
        let points = cf3d_bld_b25_four_points();
        for &(begin_center, end_center) in &[
            (Hybridization::Sp3, Hybridization::Sp2),
            (Hybridization::Sp2, Hybridization::Sp3),
        ] {
            for (first_terminal_sp2, second_terminal_sp2) in
                [(false, false), (true, false), (false, true), (true, true)]
            {
                let topology = cf3d_bld_b15_topology(
                    &[
                        if first_terminal_sp2 {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                        begin_center,
                        end_center,
                        if second_terminal_sp2 {
                            Hybridization::Sp2
                        } else {
                            Hybridization::Sp3
                        },
                    ],
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Single),
                        (2, 3, BondOrder::Single),
                    ],
                );
                let mut params = cf3d_bld_b25_params(4);
                params[0] = None;
                params[3] = None;
                let terminal_sp2 = first_terminal_sp2 || second_terminal_sp2;
                let energy = cf3d_bld_b25_assert_kernel_terms(
                    &topology,
                    &params,
                    1,
                    2,
                    &points,
                    &[(0, 1, 2, 3, terminal_sp2)],
                );
                let expected_energy = if terminal_sp2 {
                    0.2928932188134522
                } else {
                    0.4999999999999992
                };
                assert_eq!(energy, expected_energy);
            }
        }
    }

    #[test]
    fn cf3d_bld_b25_triangle_exclusion_and_terminal_branches_keep_adjacency_order() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 6],
            &[
                (1, 0, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (2, 5, BondOrder::Single),
            ],
        );
        let params = cf3d_bld_b25_params(6);
        let points = [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, -2.0, 0.0],
        ];

        cf3d_bld_b25_assert_kernel_terms(
            &topology,
            &params,
            1,
            2,
            &points,
            &[
                (0, 1, 2, 3, false),
                (0, 1, 2, 5, false),
                (4, 1, 2, 0, false),
                (4, 1, 2, 3, false),
                (4, 1, 2, 5, false),
            ],
        );

        let triangle = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 3],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        let triangle_points = [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        cf3d_bld_b25_assert_kernel_terms(
            &triangle,
            &cf3d_bld_b25_params(3),
            1,
            2,
            &triangle_points,
            &[],
        );

        let bond_only =
            cf3d_bld_b15_topology(&[Hybridization::Sp3; 2], &[(0, 1, BondOrder::Single)]);
        let bond_only_points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        cf3d_bld_b25_assert_kernel_terms(
            &bond_only,
            &cf3d_bld_b25_params(2),
            0,
            1,
            &bond_only_points,
            &[],
        );
    }

    #[test]
    fn cf3d_bld_b25_per_match_scaling_has_fixed_source_energy_and_prefix_error() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 6],
            &[
                (1, 0, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (2, 5, BondOrder::Single),
            ],
        );
        let mut params = cf3d_bld_b25_params(6);
        params[1].as_mut().expect("begin parameters").v1 = 4.0;
        params[2].as_mut().expect("end parameters").v1 = 9.0;
        params[0] = None;
        params[3] = None;
        params[4] = None;
        params[5] = None;
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let points = [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 2.0, 0.0],
        ];
        let energy = cf3d_bld_b25_assert_kernel_terms(
            &topology,
            &params,
            1,
            2,
            &points,
            &[
                (0, 1, 2, 3, false),
                (0, 1, 2, 5, false),
                (4, 1, 2, 3, false),
                (4, 1, 2, 5, false),
            ],
        );
        assert_eq!(energy, 6.0);

        let mut prefix_points = points[..4].to_vec();
        let mut prefix_field = ForceField::new(3);
        attach_positions(&mut prefix_field, &mut prefix_points);
        assert_eq!(
            torsions_for_bond(&topology, &borrowed_params, 1, 2, &mut prefix_field),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument: crate::kernel::TorsionIndexArgument::Fourth,
                    index: 5,
                    upper_bound: 4,
                }
            ))
        );
        prefix_field
            .initialize()
            .expect("source prefix remains a valid partial kernel field");
        let prefix_coordinates = points[..4]
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut prefix_field, &prefix_coordinates),
            Ok(6.0)
        );
    }

    fn cf3d_bld_b26_state(topology: &TopologyBlock) -> (RingInfo, ValenceAssignment) {
        let rings = cf3d_bld_b15_ring_info(topology);
        let rows = vec![0; topology.atoms.len()];
        (rings, assignment(&rows, &rows))
    }

    fn cf3d_bld_b26_distinct_params(atom_count: usize) -> Vec<Option<AtomicParams>> {
        let mut params = cf3d_bld_b25_params(atom_count);
        for (index, row) in params.iter_mut().enumerate() {
            let params = row.as_mut().expect("source-fixed parameter row exists");
            params.v1 = (index + 1) as f64;
            params.u1 = (index + 2) as f64;
        }
        params
    }

    fn cf3d_bld_b26_assert_source_terms(
        topology: &TopologyBlock,
        rings: &RingInfo,
        valence: &ValenceAssignment,
        params: &[Option<AtomicParams>],
        smarts: &str,
        points: &[[f64; 3]],
        expected_match_terms: &[&[(usize, usize, usize, usize, bool)]],
    ) -> (f64, Vec<f64>) {
        let coordinates = points
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();

        let mut actual_points = points.to_vec();
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);
        let borrowed_params = cf3d_bld_borrowed_params(params);
        add_torsions(
            topology,
            rings,
            valence,
            &borrowed_params,
            smarts,
            &mut actual,
        )
        .expect("source-fixed torsion dispatch succeeds");
        actual
            .initialize()
            .expect("actual source field initializes");

        let mut expected_points = points.to_vec();
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        for terms in expected_match_terms {
            cf3d_bld_b25_append_expected_terms(topology, params, &mut expected, terms);
        }
        expected
            .initialize()
            .expect("source-fixed expected field initializes");

        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &coordinates)
            .expect("actual source field energy evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &coordinates)
            .expect("source-fixed expected energy evaluates");
        assert_eq!(actual_energy, expected_energy);

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &coordinates, &mut actual_gradient)
            .expect("actual source field gradient evaluates");
        cf3d_bld_b05_calc_grad(&mut expected, &coordinates, &mut expected_gradient)
            .expect("source-fixed expected gradient evaluates");
        assert_eq!(actual_gradient, expected_gradient);

        (actual_energy, actual_gradient)
    }

    #[test]
    fn cf3d_bld_b26_default_and_custom_matches_append_overlapping_terms_in_source_order() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 6],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Single),
            ],
        );
        let (rings, valence) = cf3d_bld_b26_state(&topology);
        let params = cf3d_bld_b26_distinct_params(6);
        let points = [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 1.0, 1.0],
            [2.0, 2.0, 0.0],
        ];
        let default_matches =
            torsion_bond_matches(&topology, &rings, &valence, DEFAULT_TORSION_BOND_SMARTS)
                .expect("the pinned default query matches the fixed chain");
        assert_eq!(default_matches, vec![vec![1, 2], vec![2, 3], vec![3, 4]]);

        let first = [(0, 1, 2, 3, false)];
        let second = [(1, 2, 3, 4, false)];
        let third = [(2, 3, 4, 5, false)];
        let expected = [&first[..], &second[..], &third[..]];
        let default_result = cf3d_bld_b26_assert_source_terms(
            &topology,
            &rings,
            &valence,
            &params,
            DEFAULT_TORSION_BOND_SMARTS,
            &points,
            &expected,
        );
        let custom_result = cf3d_bld_b26_assert_source_terms(
            &topology, &rings, &valence, &params, "*~*", &points, &expected,
        );
        assert_eq!(default_result, custom_result);
    }

    #[test]
    fn cf3d_bld_b26_custom_query_preserves_reversed_endpoint_and_nested_order() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 5],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (2, 4, BondOrder::Single),
            ],
        );
        let (rings, valence) = cf3d_bld_b26_state(&topology);
        let matches = torsion_bond_matches(&topology, &rings, &valence, "[D2]~[D3]")
            .expect("the degree-directed custom pattern matches");
        assert_eq!(matches, vec![vec![2, 1]]);

        let params = cf3d_bld_b26_distinct_params(5);
        let points = [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [2.0, 1.0, 0.0],
        ];
        let reversed_terms = [(4, 2, 1, 0, false), (4, 2, 1, 3, false)];
        let expected = [&reversed_terms[..]];
        cf3d_bld_b26_assert_source_terms(
            &topology,
            &rings,
            &valence,
            &params,
            "[D2]~[D3]",
            &points,
            &expected,
        );
    }

    #[test]
    fn cf3d_bld_b26_preserves_precondition_query_and_bond_error_order() {
        let chain = cf3d_bld_b15_topology(&[Hybridization::Sp3; 2], &[(0, 1, BondOrder::Single)]);
        let (rings, valence) = cf3d_bld_b26_state(&chain);
        let params = cf3d_bld_b25_params(2);
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut points);

        assert_eq!(
            add_torsions(
                &chain,
                &rings,
                &valence,
                &borrowed_params[..1],
                "[",
                &mut field
            ),
            Err(UffBuilderError::ParamsLengthMismatch {
                atoms: 2,
                params: 1,
            })
        );
        assert!(matches!(
            add_torsions(&chain, &rings, &valence, &borrowed_params, "[", &mut field),
            Err(UffBuilderError::TorsionBondQuery(
                TorsionBondQueryError::Parse(_)
            ))
        ));
        assert_eq!(
            add_torsions(&chain, &rings, &valence, &borrowed_params, "*", &mut field),
            Err(UffBuilderError::TorsionBondQuery(
                TorsionBondQueryError::MatchArity {
                    smarts: "*".to_owned(),
                    actual: 1,
                }
            ))
        );

        let disconnected = cf3d_bld_b15_topology(&[Hybridization::Sp3; 2], &[]);
        let (rings, valence) = cf3d_bld_b26_state(&disconnected);
        let matched = torsion_bond_matches(&disconnected, &rings, &valence, "*.*")
            .expect("the fixed disconnected query produces one source match");
        assert_eq!(matched[0], vec![0, 1]);
        assert_eq!(
            add_torsions(
                &disconnected,
                &rings,
                &valence,
                &borrowed_params,
                "*.*",
                &mut field
            ),
            Err(UffBuilderError::TorsionBondQuery(
                TorsionBondQueryError::MissingMatchedBond {
                    begin_atom_index: 0,
                    end_atom_index: 1,
                }
            ))
        );

        let mut missing_center = params.clone();
        missing_center[0] = None;
        let missing_center_borrowed_params = cf3d_bld_borrowed_params(&missing_center);
        assert_eq!(
            add_torsions(
                &disconnected,
                &rings,
                &valence,
                &missing_center_borrowed_params,
                "*.*",
                &mut field,
            ),
            Ok(())
        );
    }

    #[test]
    fn cf3d_bld_b26_later_constructor_error_keeps_prior_match_terms() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 5],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
            ],
        );
        let (rings, valence) = cf3d_bld_b26_state(&topology);
        let params = cf3d_bld_b25_params(5);
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let points = cf3d_bld_b25_four_points();
        let coordinates = points
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points;
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);

        assert_eq!(
            add_torsions(
                &topology,
                &rings,
                &valence,
                &borrowed_params,
                DEFAULT_TORSION_BOND_SMARTS,
                &mut actual,
            ),
            Err(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument: crate::kernel::TorsionIndexArgument::Fourth,
                    index: 4,
                    upper_bound: 4,
                }
            ))
        );

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        cf3d_bld_b25_append_expected_terms(
            &topology,
            &params,
            &mut expected,
            &[(0, 1, 2, 3, false)],
        );
        actual.initialize().expect("the retained term initializes");
        expected
            .initialize()
            .expect("the source-fixed retained term initializes");

        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &coordinates)
            .expect("the retained source term evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &coordinates)
            .expect("the source-fixed term evaluates");
        assert_eq!(actual_energy, expected_energy);
        // Pinned TorsionAngle.cpp selects V=1, n=3, cosTerm=-1 for this
        // single C(sp3)-C(sp3) bond with both source V1 values equal to 1.
        // The fixed coordinates yield cosPhi=0.7071067811865475,
        // sinPhiSq=0.5000000000000001 and cosNPhi=-0.7071067811865478;
        // the source energy equation therefore gives 0.1464466094067261.
        assert_eq!(actual_energy, 0.1464466094067261);

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &coordinates, &mut actual_gradient)
            .expect("the retained source gradient evaluates");
        cf3d_bld_b05_calc_grad(&mut expected, &coordinates, &mut expected_gradient)
            .expect("the source-fixed gradient evaluates");
        assert_eq!(actual_gradient, expected_gradient);
    }

    fn cf3d_bld_b27_topology(
        center_atomic_number: u8,
        center_hybridization: Hybridization,
        neighbor_states: &[(u8, Hybridization)],
        neighbor_order: &[usize],
    ) -> TopologyBlock {
        assert_eq!(neighbor_states.len(), neighbor_order.len());

        let center_element = Element::from_atomic_number(center_atomic_number)
            .expect("the fixed atomic number is in the source periodic table");
        let mut atoms = Vec::with_capacity(neighbor_states.len() + 1);
        atoms.push(Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(center_element).with_hybridization(center_hybridization),
        ));
        for (neighbor_offset, &(atomic_number, hybridization)) in neighbor_states.iter().enumerate()
        {
            let element = Element::from_atomic_number(atomic_number)
                .expect("the fixed neighbor atomic number is in the source periodic table");
            atoms.push(Atom::from_spec(
                AtomId::new(neighbor_offset + 1),
                AtomSpec::new(element).with_hybridization(hybridization),
            ));
        }

        let bonds = neighbor_order
            .iter()
            .enumerate()
            .map(|(bond_row, &neighbor_offset)| {
                Bond::from_spec(
                    BondId::new(bond_row),
                    BondSpec::new(
                        AtomId::new(0),
                        AtomId::new(neighbor_offset + 1),
                        BondOrder::Single,
                    ),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed B27 topology is structurally valid")
    }

    #[test]
    fn cf3d_bld_b27_element_degree_hybridization_matrix() {
        for center_atomic_number in 0_u8..=118 {
            for degree in 2_usize..=4 {
                let neighbor_order: Vec<_> = (0..degree).rev().collect();
                let neighbors = vec![(1, Hybridization::Sp3); degree];

                for center_hybridization in [Hybridization::Sp2, Hybridization::Sp3] {
                    let topology = cf3d_bld_b27_topology(
                        center_atomic_number,
                        center_hybridization,
                        &neighbors,
                        &neighbor_order,
                    );
                    let source_element_is_allowed =
                        matches!(center_atomic_number, 6 | 7 | 8 | 15 | 33 | 51 | 83);
                    let source_hybridization_is_allowed =
                        !matches!(center_atomic_number, 6 | 7 | 8)
                            || center_hybridization == Hybridization::Sp2;
                    let expected = if degree == 3
                        && source_element_is_allowed
                        && source_hybridization_is_allowed
                    {
                        Some((center_atomic_number, [3, 0, 2, 1], false))
                    } else {
                        None
                    };

                    assert_eq!(
                        inversion_center(&topology, 0),
                        expected,
                        "Z={center_atomic_number}, degree={degree}, hybridization={center_hybridization:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn cf3d_bld_b27_neighbor_order_and_carbon_sp2_oxygen_flag() {
        let neighbor_order = [2, 0, 1];
        for oxygen_position in 0..3 {
            let mut neighbors = [(1, Hybridization::Sp3); 3];
            neighbors[neighbor_order[oxygen_position]] = (8, Hybridization::Sp2);
            let carbon = cf3d_bld_b27_topology(6, Hybridization::Sp2, &neighbors, &neighbor_order);
            assert_eq!(
                inversion_center(&carbon, 0),
                Some((6, [3, 0, 1, 2], true)),
                "SP2 oxygen in source-neighbor position {oxygen_position}"
            );

            neighbors[neighbor_order[oxygen_position]] = (8, Hybridization::Sp3);
            let carbon_without_sp2_oxygen =
                cf3d_bld_b27_topology(6, Hybridization::Sp2, &neighbors, &neighbor_order);
            assert_eq!(
                inversion_center(&carbon_without_sp2_oxygen, 0),
                Some((6, [3, 0, 1, 2], false)),
                "SP3 oxygen in source-neighbor position {oxygen_position}"
            );

            for (center_atomic_number, center_hybridization) in [
                (7, Hybridization::Sp2),
                (8, Hybridization::Sp2),
                (15, Hybridization::Sp3),
                (33, Hybridization::Sp3),
                (51, Hybridization::Sp3),
                (83, Hybridization::Sp3),
            ] {
                neighbors[neighbor_order[oxygen_position]] = (8, Hybridization::Sp2);
                let noncarbon = cf3d_bld_b27_topology(
                    center_atomic_number,
                    center_hybridization,
                    &neighbors,
                    &neighbor_order,
                );
                assert_eq!(
                    inversion_center(&noncarbon, 0),
                    Some((center_atomic_number, [3, 0, 1, 2], false)),
                    "central Z={center_atomic_number} must not acquire carbon-only oxygen flag"
                );
            }
        }
    }

    fn cf3d_bld_b28_topology(center_atomic_numbers: &[u8]) -> TopologyBlock {
        let mut atoms = Vec::with_capacity(center_atomic_numbers.len() * 4);
        let mut bonds = Vec::with_capacity(center_atomic_numbers.len() * 3);

        for (group, &center_atomic_number) in center_atomic_numbers.iter().enumerate() {
            let base = group * 4;
            let center_element = Element::from_atomic_number(center_atomic_number)
                .expect("fixed source inversion elements are present");
            let center_hybridization = if matches!(center_atomic_number, 6 | 7 | 8) {
                Hybridization::Sp2
            } else {
                Hybridization::Sp3
            };
            atoms.push(Atom::from_spec(
                AtomId::new(base),
                AtomSpec::new(center_element).with_hybridization(center_hybridization),
            ));

            for neighbor_offset in 0..3 {
                let (atomic_number, hybridization) =
                    if center_atomic_number == 6 && neighbor_offset == 0 {
                        (8, Hybridization::Sp2)
                    } else {
                        (1, Hybridization::Sp3)
                    };
                let element = Element::from_atomic_number(atomic_number)
                    .expect("fixed source inversion neighbor is present");
                atoms.push(Atom::from_spec(
                    AtomId::new(base + neighbor_offset + 1),
                    AtomSpec::new(element).with_hybridization(hybridization),
                ));
            }

            // Source adjacency order is neighbor 2, neighbor 0, neighbor 1.
            for neighbor_offset in [2, 0, 1] {
                let bond_row = bonds.len();
                bonds.push(Bond::from_spec(
                    BondId::new(bond_row),
                    BondSpec::new(
                        AtomId::new(base),
                        AtomId::new(base + neighbor_offset + 1),
                        BondOrder::Single,
                    ),
                ));
            }
        }

        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed B28 disconnected centers are structurally valid")
    }

    fn cf3d_bld_b28_one_center_topology(
        center_atom_index: usize,
        source_neighbor_indices: [usize; 3],
    ) -> TopologyBlock {
        let atoms = (0..4)
            .map(|row| {
                if row == center_atom_index {
                    Atom::from_spec(
                        AtomId::new(row),
                        AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
                    )
                } else {
                    Atom::from_spec(
                        AtomId::new(row),
                        AtomSpec::new(Element::H).with_hybridization(Hybridization::Sp3),
                    )
                }
            })
            .collect();
        let bonds = source_neighbor_indices
            .into_iter()
            .enumerate()
            .map(|(bond_row, neighbor_index)| {
                Bond::from_spec(
                    BondId::new(bond_row),
                    BondSpec::new(
                        AtomId::new(center_atom_index),
                        AtomId::new(neighbor_index),
                        BondOrder::Single,
                    ),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed B28 range-error topology is structurally valid")
    }

    fn cf3d_bld_b28_points(center_count: usize) -> Vec<[f64; 3]> {
        let one_center = cf3d_bld_b25_four_points();
        let mut points = Vec::with_capacity(center_count * 4);
        for group in 0..center_count {
            let x_offset = group as f64 * 4.0;
            for point in &one_center {
                points.push([point[0] + x_offset, point[1], point[2]]);
            }
        }
        points
    }

    fn cf3d_bld_b28_append_expected_terms(
        field: &mut ForceField<'_>,
        terms: [[usize; 4]; 3],
        atomic_number: u8,
        is_c_bound_to_o: bool,
    ) {
        for term in terms {
            let [idx1, idx2, idx3, idx4] = term
                .map(|index| u32::try_from(index).expect("fixed B28 source atom index fits u32"));
            let contribution = InversionContrib::new(
                field.positions(),
                idx1,
                idx2,
                idx3,
                idx4,
                i32::from(atomic_number),
                is_c_bound_to_o,
            )
            .expect("fixed source term indices are present in the field");
            field.add_contribution(Box::new(contribution));
        }
    }

    #[test]
    fn cf3d_bld_b28_appends_source_permutations_for_all_centers_without_parameter_rows() {
        let center_atomic_numbers = [6, 7, 8, 15, 33, 51, 83];
        let topology = cf3d_bld_b28_topology(&center_atomic_numbers);
        let params = vec![None; topology.atoms.len()];
        let points = cf3d_bld_b28_points(center_atomic_numbers.len());
        let coordinates = points
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();

        let mut actual_points = points.clone();
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);
        assert_eq!(add_inversions(&topology, &params, &mut actual), Ok(()));

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        for (group, &atomic_number) in center_atomic_numbers.iter().enumerate() {
            let base = group * 4;
            cf3d_bld_b28_append_expected_terms(
                &mut expected,
                [
                    [base + 3, base, base + 1, base + 2],
                    [base + 3, base, base + 2, base + 1],
                    [base + 1, base, base + 2, base + 3],
                ],
                atomic_number,
                atomic_number == 6,
            );
        }

        actual
            .initialize()
            .expect("all fixed source terms initialize");
        expected
            .initialize()
            .expect("the source-fixed terms initialize");
        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &coordinates)
            .expect("all source inversion terms evaluate");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &coordinates)
            .expect("all fixed source terms evaluate");
        assert_ne!(actual_energy, 0.0);
        assert_eq!(actual_energy, expected_energy);

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &coordinates, &mut actual_gradient)
            .expect("all source inversion gradients evaluate");
        cf3d_bld_b05_calc_grad(&mut expected, &coordinates, &mut expected_gradient)
            .expect("all source-fixed gradients evaluate");
        assert_eq!(actual_gradient, expected_gradient);
    }

    #[test]
    fn cf3d_bld_b28_preserves_error_order_center_order_and_prior_terms() {
        let argument_cases = [
            (0, [3, 1, 2], InversionIndexArgument::First),
            (3, [0, 1, 2], InversionIndexArgument::Second),
            (0, [1, 3, 2], InversionIndexArgument::Third),
            (0, [1, 2, 3], InversionIndexArgument::Fourth),
        ];
        for (center_index, neighbor_indices, argument) in argument_cases {
            let topology = cf3d_bld_b28_one_center_topology(center_index, neighbor_indices);
            let params = vec![None; topology.atoms.len()];
            let mut points = [[0.0; 3]; 3];
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut points);

            assert_eq!(
                add_inversions(&topology, &params, &mut field),
                Err(UffBuilderError::InversionContribution(
                    InversionContributionError::IndexOutOfRange {
                        argument,
                        index: 3,
                        upper_bound: 3,
                    }
                )),
                "the constructor keeps its source argument check order"
            );
        }

        let topology = cf3d_bld_b28_topology(&[6, 7]);
        let params = vec![None; topology.atoms.len()];
        let points = cf3d_bld_b28_points(1);
        let coordinates = points
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();
        let mut actual_points = points.clone();
        let mut actual = ForceField::new(3);
        attach_positions(&mut actual, &mut actual_points);

        assert_eq!(
            add_inversions(&topology, &params, &mut actual),
            Err(UffBuilderError::InversionContribution(
                InversionContributionError::IndexOutOfRange {
                    argument: InversionIndexArgument::First,
                    index: 7,
                    upper_bound: 4,
                }
            )),
            "the later center fails after the first center's terms were appended"
        );

        let mut expected_points = points;
        let mut expected = ForceField::new(3);
        attach_positions(&mut expected, &mut expected_points);
        cf3d_bld_b28_append_expected_terms(
            &mut expected,
            [[3, 0, 1, 2], [3, 0, 2, 1], [1, 0, 2, 3]],
            6,
            true,
        );
        actual
            .initialize()
            .expect("the earlier terms remain installed");
        expected
            .initialize()
            .expect("the source-fixed earlier terms initialize");
        let actual_energy = cf3d_bld_b05_calc_energy(&mut actual, &coordinates)
            .expect("the retained source prefix evaluates");
        let expected_energy = cf3d_bld_b05_calc_energy(&mut expected, &coordinates)
            .expect("the fixed source prefix evaluates");
        assert_ne!(actual_energy, 0.0);
        assert_eq!(actual_energy, expected_energy);

        let mut actual_gradient = vec![0.0; coordinates.len()];
        let mut expected_gradient = vec![0.0; coordinates.len()];
        cf3d_bld_b05_calc_grad(&mut actual, &coordinates, &mut actual_gradient)
            .expect("the retained source gradients evaluate");
        cf3d_bld_b05_calc_grad(&mut expected, &coordinates, &mut expected_gradient)
            .expect("the source-fixed gradients evaluate");
        assert_eq!(actual_gradient, expected_gradient);
    }

    fn f19_topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = elements
            .iter()
            .copied()
            .enumerate()
            .map(|(row, element)| Atom::from_spec(AtomId::new(row), AtomSpec::new(element)))
            .collect();
        let bonds = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(row, (begin, end))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed F19 topology is structurally valid")
    }

    fn f19_connected_overvalent_topology() -> TopologyBlock {
        f19_topology(
            &[Element::O, Element::C, Element::C, Element::C],
            &[(0, 1), (0, 2), (0, 3)],
        )
    }

    fn f19_coordinates(atom_count: usize) -> CoordinateBlock {
        let conformers_2d = vec![
            Conformer2D::new(
                41,
                (0..atom_count)
                    .map(|row| [2.0 * row as f64, 2.0 * row as f64 + 1.0])
                    .collect(),
            )
            .with_prop("origin", "2d-41"),
        ];
        let conformers_3d = [(73, false), (11, true), (5, false)]
            .into_iter()
            .map(|(id, is_3d)| {
                let offset = match id {
                    73 => [0.0, 10.0, 20.0],
                    11 => [1000.0, 1010.0, 1020.0],
                    5 => [0.0, 20.0, 40.0],
                    _ => unreachable!("fixed F19 conformer ID"),
                };
                Conformer3D::new(
                    id,
                    (0..atom_count)
                        .map(|row| {
                            [
                                row as f64 + offset[0],
                                row as f64 + offset[1],
                                row as f64 + offset[2],
                            ]
                        })
                        .collect(),
                    is_3d,
                )
                .with_prop("origin", format!("3d-{id}"))
            })
            .collect();
        CoordinateBlock {
            conformers_2d,
            conformers_3d,
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
        }
    }

    fn cf3d_frag_f22_coordinates(points: &[[f64; 3]]) -> CoordinateBlock {
        let mut coordinates = f19_coordinates(points.len());
        coordinates.conformers_3d[1]
            .coordinates_mut()
            .copy_from_slice(points);
        coordinates
    }

    fn cf3d_frag_f22_append<'a>(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        molecule_properties: &MoleculeProperties,
        params: &[Option<AtomicParams>],
        rows: &'a mut [[f64; 3]],
        neighbor_matrix: &[u8],
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Result<ForceField<'a>, NonbondedAssemblyError> {
        let borrowed_params = cf3d_bld_borrowed_params(params);
        let selected_conformer = &coordinates.conformers_3d[1];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, rows);
        add_nonbonded(
            topology,
            &coordinates.conformers_2d,
            &coordinates.conformers_3d[..1],
            selected_conformer.id(),
            selected_conformer.is_3d(),
            selected_conformer.props(),
            &coordinates.conformers_3d[2..],
            coordinates.source_coordinate_dim,
            molecule_properties,
            &borrowed_params,
            &mut field,
            neighbor_matrix,
            vdw_threshold,
            ignore_interfragment_interactions,
        )?;
        Ok(field)
    }

    #[test]
    fn cf3d_frag_f19_false_skips_fragment_copy_for_valid_and_invalid_inputs() {
        let coordinates = f19_coordinates(4);
        let properties = MoleculeProperties::default();
        for topology in [
            topology_with_bonds(4, &[(0, 1, false), (2, 3, false)]),
            f19_connected_overvalent_topology(),
        ] {
            let original = (topology.clone(), coordinates.clone(), properties.clone());
            let coordinate_view = FragmentCoordinateView::from_coordinate_block(&coordinates);
            assert!(prepare_nonbonded_fragment_mapping(
                &topology,
                &coordinate_view,
                &properties,
                false,
            )
            .unwrap()
            .is_none());
            assert_eq!(
                (topology, coordinates.clone(), properties.clone()),
                original
            );
        }
    }

    #[test]
    fn cf3d_frag_f19_true_propagates_connected_sanitize_failure() {
        let topology = f19_connected_overvalent_topology();
        let coordinates = f19_coordinates(4);
        let properties = MoleculeProperties::default();
        let original = (topology.clone(), coordinates.clone(), properties.clone());
        let coordinate_view = FragmentCoordinateView::from_coordinate_block(&coordinates);

        let error =
            prepare_nonbonded_fragment_mapping(&topology, &coordinate_view, &properties, true)
                .unwrap_err();

        assert!(matches!(
            error,
            NonbondedFragmentMappingError::FragmentCopy(source)
                if source.component_index() == Some(0)
        ));
        assert_eq!((topology, coordinates, properties), original);
    }

    #[test]
    fn cf3d_frag_f19_true_propagates_disconnected_first_failing_fragment() {
        let topology = f19_topology(
            &[
                Element::C,
                Element::C,
                Element::O,
                Element::C,
                Element::C,
                Element::C,
            ],
            &[(0, 1), (2, 3), (2, 4), (2, 5)],
        );
        let coordinates = f19_coordinates(6);
        let properties = MoleculeProperties::default();
        let original = (topology.clone(), coordinates.clone(), properties.clone());
        let coordinate_view = FragmentCoordinateView::from_coordinate_block(&coordinates);

        let error =
            prepare_nonbonded_fragment_mapping(&topology, &coordinate_view, &properties, true)
                .unwrap_err();

        assert!(matches!(
            error,
            NonbondedFragmentMappingError::FragmentCopy(source)
                if source.component_index() == Some(1)
        ));
        assert_eq!((topology, coordinates, properties), original);
    }

    #[test]
    fn cf3d_frag_f19_true_returns_source_ordered_labels_and_copies_kernel_view() {
        let topology = topology_with_bonds(4, &[(0, 1, false), (2, 3, false)]);
        let coordinates = f19_coordinates(4);
        let original_coordinates = coordinates.clone();
        let properties = MoleculeProperties::default();
        let selected_id = 11;
        let selected_is_3d = coordinates.conformers_3d[1].is_3d();
        let selected_props = coordinates.conformers_3d[1].props().clone();
        let source_dimension = coordinates.source_coordinate_dim;
        let mut kernel_rows = [
            [100.0, 101.0, 102.0],
            [110.0, 111.0, 112.0],
            [120.0, 121.0, 122.0],
            [130.0, 131.0, 132.0],
        ];
        let original_kernel_rows = kernel_rows;
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut kernel_rows);
        let selected_kernel_rows = field
            .positions()
            .iter()
            .map(|row| &row[..])
            .collect::<Vec<_>>();
        let coordinate_view = FragmentCoordinateView::from_split_conformers(
            &coordinates.conformers_2d,
            &coordinates.conformers_3d[..1],
            selected_id,
            selected_is_3d,
            &selected_props,
            &selected_kernel_rows,
            &coordinates.conformers_3d[2..],
            source_dimension,
        )
        .expect("fixed kernel rows have three values each");

        let labels =
            prepare_nonbonded_fragment_mapping(&topology, &coordinate_view, &properties, true)
                .unwrap()
                .expect("the true source branch creates a mapping");
        assert_eq!(labels, vec![0, 0, 1, 1]);

        let fragments = get_molecule_fragments_with_coordinate_view(
            &topology,
            &coordinate_view,
            &properties,
            true,
            true,
        )
        .unwrap();
        assert_eq!(fragments.len(), 2);
        assert_eq!(
            fragments[0].component_atoms(),
            &[AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(
            fragments[1].component_atoms(),
            &[AtomId::new(2), AtomId::new(3)]
        );
        for fragment in &fragments {
            let copied = fragment.coordinates();
            assert_eq!(copied.source_coordinate_dim, source_dimension);
            assert_eq!(copied.conformers_2d.len(), 1);
            assert_eq!(copied.conformers_2d[0].id(), 41);
            assert_eq!(
                copied.conformers_2d[0].props().get("origin").unwrap(),
                "2d-41"
            );
            assert_eq!(
                copied
                    .conformers_3d
                    .iter()
                    .map(Conformer3D::id)
                    .collect::<Vec<_>>(),
                vec![73, 11, 5]
            );
            assert_eq!(
                copied
                    .conformers_3d
                    .iter()
                    .map(Conformer3D::is_3d)
                    .collect::<Vec<_>>(),
                vec![false, true, false]
            );
            assert_eq!(
                copied
                    .conformers_3d
                    .iter()
                    .map(|conformer| conformer.props().get("origin").unwrap().as_str())
                    .collect::<Vec<_>>(),
                vec!["3d-73", "3d-11", "3d-5"]
            );
        }
        assert_eq!(
            fragments[0].coordinates().conformers_2d[0].coordinates(),
            &[[0.0, 1.0], [2.0, 3.0]]
        );
        assert_eq!(
            fragments[0].coordinates().conformers_3d[0].coordinates(),
            &[[0.0, 10.0, 20.0], [1.0, 11.0, 21.0]]
        );
        assert_eq!(
            fragments[0].coordinates().conformers_3d[1].coordinates(),
            &[[100.0, 101.0, 102.0], [110.0, 111.0, 112.0]]
        );
        assert_eq!(
            fragments[0].coordinates().conformers_3d[2].coordinates(),
            &[[0.0, 20.0, 40.0], [1.0, 21.0, 41.0]]
        );
        assert_eq!(coordinates, original_coordinates);
        drop(coordinate_view);
        drop(selected_kernel_rows);
        assert_eq!(
            field
                .positions()
                .iter()
                .map(|row| [row[0], row[1], row[2]])
                .collect::<Vec<_>>(),
            original_kernel_rows
        );
    }

    #[test]
    fn cf3d_frag_f20_pair_guard_covers_parameter_flag_component_and_relation_matrix() {
        let topology = topology_with_bonds(2, &[]);
        let cell = two_bit_cell_pos(2, 0, 1).expect("fixed two-atom pair cell");
        assert_eq!(cell, 1);
        let initialized_matrix = build_neighbor_matrix(&topology).expect("fixed matrix");

        for relation in 0_u8..=3 {
            let mut neighbor_matrix = initialized_matrix.clone();
            set_two_bit_cell(&mut neighbor_matrix, cell, relation)
                .expect("fixed relation cell fits");
            for (has_i_params, has_j_params) in
                [(false, false), (false, true), (true, false), (true, true)]
            {
                let parameter_rows = [
                    has_i_params.then(|| atomic_params(1.0, 6.0)),
                    has_j_params.then(|| atomic_params(1.5, 8.0)),
                ];
                let params = cf3d_bld_borrowed_params(&parameter_rows);
                for ignore_interfragment_interactions in [false, true] {
                    for different_components in [false, true] {
                        let fragment_mapping = [0, if different_components { 1 } else { 0 }];
                        let expected = has_i_params
                            && has_j_params
                            && (!ignore_interfragment_interactions || !different_components)
                            && relation >= 2;
                        assert_eq!(
                            nonbonded_pair_is_eligible(
                                2,
                                0,
                                1,
                                &params,
                                &fragment_mapping,
                                ignore_interfragment_interactions,
                                &neighbor_matrix,
                            ),
                            Ok(expected),
                            "i params={has_i_params}, j params={has_j_params}, ignore fragments={ignore_interfragment_interactions}, different components={different_components}, relation={relation}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn cf3d_frag_f20_pair_guard_preserves_short_circuit_and_matrix_error_order() {
        let first = atomic_params(1.0, 6.0);
        let second = atomic_params(1.5, 8.0);

        let missing_i_params = [None, Some(second.clone())];
        let missing_i_params = cf3d_bld_borrowed_params(&missing_i_params);
        assert_eq!(
            nonbonded_pair_is_eligible(2, 0, 1, &missing_i_params, &[], true, &[]),
            Ok(false),
            "missing i params skips j, fragment-map and matrix reads"
        );

        let missing_j_params = [Some(first.clone()), None];
        let missing_j_params = cf3d_bld_borrowed_params(&missing_j_params);
        assert_eq!(
            nonbonded_pair_is_eligible(2, 0, 1, &missing_j_params, &[], true, &[]),
            Ok(false),
            "missing j params skips fragment-map and matrix reads"
        );

        let both_params = [Some(first), Some(second)];
        let both_params = cf3d_bld_borrowed_params(&both_params);
        assert_eq!(
            nonbonded_pair_is_eligible(2, 0, 1, &both_params, &[0, 1], true, &[]),
            Ok(false),
            "different component labels skip the matrix read"
        );

        let cell = two_bit_cell_pos(2, 0, 1).expect("fixed pair cell");
        assert_eq!(cell, 1);
        assert_eq!(
            nonbonded_pair_is_eligible(2, 0, 1, &both_params, &[], false, &[]),
            Err(UffBuilderError::NeighborMatrixStorageOutOfRange {
                position: 1,
                byte_index: 0,
                storage_len: 0,
            }),
            "with fragment checks disabled, valid params reach the relation read"
        );
    }

    #[test]
    fn cf3d_frag_f21_uses_strict_distance_cutoff_and_canonical_vdw_params() {
        let params_i = cf3d_frag_f21_params(2.0, 3.0);
        let params_j = cf3d_frag_f21_params(8.0, 12.0);

        // The source geometric means are xij=sqrt(2*8)=4 and
        // D=sqrt(3*12)=6. At distance 2, r=2 and D*(r^12-2*r^6)=23808.
        // The addNonbonded pair cutoff is strictly 2*xij=8.
        for (distance, expected_energy) in [(2.0, 23808.0), (8.0, 0.0), (9.0, 0.0)] {
            let mut rows = [[0.0, 0.0, 0.0], [distance, 0.0, 0.0]];
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut rows);

            append_nonbonded_pair_if_within_threshold(&mut field, 0, 1, &params_i, &params_j, 2.0)
                .expect("source-valid pair positions append or skip without error");
            field.initialize().expect("fixed F21 field initializes");

            let coordinates = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
            let energy = cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("fixed F21 field energy evaluates");
            assert_eq!(energy, expected_energy, "pair distance {distance}");
        }
    }

    #[test]
    fn cf3d_frag_f21_preserves_coincident_zero_and_negative_threshold_branches() {
        let cases: [([[f64; 3]; 2], [f64; 6], f64, [f64; 6]); 3] = [
            (
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                [0.0; 6],
                2.0,
                [100.0, 100.0, 100.0, -100.0, -100.0, -100.0],
            ),
            ([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [0.0; 6], 0.0, [0.0; 6]),
            (
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                -1.0,
                [0.0; 6],
            ),
        ];
        let params_i = cf3d_frag_f21_params(1.0, 1.0);
        let params_j = cf3d_frag_f21_params(1.0, 1.0);

        for (rows, coordinates, threshold, expected_gradient) in cases {
            let mut rows = rows;
            let mut field = ForceField::new(3);
            attach_positions(&mut field, &mut rows);
            append_nonbonded_pair_if_within_threshold(
                &mut field, 0, 1, &params_i, &params_j, threshold,
            )
            .expect("source-valid pair positions append or skip without error");
            field.initialize().expect("fixed F21 field initializes");

            assert_eq!(
                cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                    .expect("fixed F21 field energy evaluates"),
                0.0
            );
            let mut gradient = [0.0; 6];
            cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
                .expect("fixed F21 field gradient evaluates");
            assert_eq!(gradient, expected_gradient);
        }
    }

    #[test]
    fn cf3d_frag_f21_leaves_nonfinite_parameter_flow_to_source_arithmetic() {
        let coordinates = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let mut rows = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut rows);
        let nan_minimum = cf3d_frag_f21_params(f64::NAN, 1.0);
        let finite = cf3d_frag_f21_params(1.0, 1.0);

        append_nonbonded_pair_if_within_threshold(&mut field, 0, 1, &nan_minimum, &finite, 2.0)
            .expect("NaN source comparison is a threshold miss");
        field.initialize().expect("fixed F21 field initializes");
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates),
            Ok(0.0),
            "distance < NaN is false, so no contribution is appended"
        );

        let mut rows = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut rows);
        let nan_depth = cf3d_frag_f21_params(1.0, f64::NAN);
        append_nonbonded_pair_if_within_threshold(&mut field, 0, 1, &nan_depth, &finite, 2.0)
            .expect("finite source minimum permits the pair");
        field.initialize().expect("fixed F21 field initializes");
        assert!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("source nonfinite energy propagates")
                .is_nan(),
            "AtomicParams has no finite-value guard; the canonical depth arithmetic remains NaN"
        );
    }

    #[test]
    fn cf3d_frag_f21_appends_terms_in_source_pair_order() {
        let cancelling_depth = 31.0 * 2.0_f64.powi(59);
        let params_0 = cf3d_frag_f21_params(1.0, 1.0);
        let params_1 = cf3d_frag_f21_params(1.0, 2.0_f64.powi(104));
        let params_2 = cf3d_frag_f21_params(1.0, cancelling_depth * cancelling_depth);
        let params_3 = cf3d_frag_f21_params(1.0, 1.0);
        let mut rows = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        let mut field = ForceField::new(3);
        attach_positions(&mut field, &mut rows);

        // These disjoint candidate pairs are called in the source's ascending
        // nested order: (0,1), (0,2), (0,3). Their fixed energies are
        // +31*2^59, -31*2^59, -1; ordered accumulation preserves the final -1.
        append_nonbonded_pair_if_within_threshold(&mut field, 0, 1, &params_0, &params_1, 2.0)
            .expect("first source-ordered pair appends");
        append_nonbonded_pair_if_within_threshold(&mut field, 0, 2, &params_0, &params_2, 2.0)
            .expect("second source-ordered pair appends");
        append_nonbonded_pair_if_within_threshold(&mut field, 0, 3, &params_0, &params_3, 2.0)
            .expect("third source-ordered pair appends");
        field.initialize().expect("fixed F21 field initializes");

        let coordinates = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("source-ordered F21 terms evaluate"),
            -1.0
        );
    }

    #[test]
    fn cf3d_frag_f22_zero_and_one_atom_inputs_skip_pair_state() {
        let one_atom_params = cf3d_frag_f21_params(1.0, 1.0);
        let cases = [
            (Vec::<Element>::new(), Vec::<Option<AtomicParams>>::new()),
            (vec![Element::C], vec![None]),
            (vec![Element::C], vec![Some(one_atom_params)]),
        ];

        for (elements, params) in cases {
            let topology = f19_topology(&elements, &[]);
            let points = (0..elements.len())
                .map(|row| [row as f64, 0.0, 0.0])
                .collect::<Vec<_>>();
            let coordinates = cf3d_frag_f22_coordinates(&points);
            let molecule_properties = MoleculeProperties::default().with_name("F22-small");
            let original = (
                topology.clone(),
                coordinates.clone(),
                molecule_properties.clone(),
            );
            let mut rows = coordinates.conformers_3d[1].coordinates().to_vec();
            let original_rows = rows.clone();

            let field = cf3d_frag_f22_append(
                &topology,
                &coordinates,
                &molecule_properties,
                &params,
                &mut rows,
                &[],
                2.0,
                false,
            )
            .expect("zero/one atom source loops never read a pair cell");
            drop(field);

            assert_eq!(
                (topology, coordinates, molecule_properties),
                original,
                "the detached inputs remain unchanged"
            );
            assert_eq!(rows, original_rows, "borrowed kernel rows remain unchanged");
        }
    }

    #[test]
    fn cf3d_frag_f22_pair_matrix_preserves_masks_relations_fragments_and_cutoff() {
        // Builder.cpp visits each (i,j) only after both source parameter
        // optionals, the requested fragment labels, relation >= 1_4, and the
        // strict distance threshold have passed. Fixed energy/gradient values
        // below come from the existing UFF vdW formula for x1=d1=1 at r=0.5.
        const INCLUDED_ENERGY: f64 = 3968.0;
        const INCLUDED_GRADIENT: [f64; 6] = [96_768.0, 0.0, 0.0, -96_768.0, 0.0, 0.0];
        let molecule_properties = MoleculeProperties::default().with_name("F22-matrix");

        for parameter_mask in 0_u8..4 {
            for ignore_interfragment_interactions in [false, true] {
                for same_component in [false, true] {
                    let topology = if same_component {
                        topology_with_bonds(2, &[(0, 1, false)])
                    } else {
                        topology_with_bonds(2, &[])
                    };
                    let original_topology = topology.clone();
                    for relation in 0_u8..4 {
                        let mut neighbor_matrix = vec![0_u8; 1];
                        let pair_cell = two_bit_cell_pos(2, 0, 1).expect("fixed pair cell");
                        set_two_bit_cell(&mut neighbor_matrix, pair_cell, relation)
                            .expect("one byte holds all three two-atom cells");

                        for (distance, expected_pair) in [(0.5, true), (2.0, false), (3.0, false)] {
                            let points = [[0.0, 0.0, 0.0], [distance, 0.0, 0.0]];
                            let coordinates = cf3d_frag_f22_coordinates(&points);
                            let original_coordinates = coordinates.clone();
                            let params = [
                                (parameter_mask & 1 != 0).then(|| cf3d_frag_f21_params(1.0, 1.0)),
                                (parameter_mask & 2 != 0).then(|| cf3d_frag_f21_params(1.0, 1.0)),
                            ];
                            let eligible = expected_pair
                                && parameter_mask == 3
                                && relation >= 2
                                && (!ignore_interfragment_interactions || same_component);
                            let mut rows = coordinates.conformers_3d[1].coordinates().to_vec();
                            let original_rows = rows.clone();
                            let mut field = cf3d_frag_f22_append(
                                &topology,
                                &coordinates,
                                &molecule_properties,
                                &params,
                                &mut rows,
                                &neighbor_matrix,
                                2.0,
                                ignore_interfragment_interactions,
                            )
                            .expect("valid fixed F22 source state assembles");
                            field
                                .initialize()
                                .expect("fixed F22 kernel field initializes");

                            let coordinates_flat = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
                            let energy = cf3d_bld_b05_calc_energy(&mut field, &coordinates_flat)
                                .expect("real F22 kernel energy evaluates");
                            let mut gradient = [0.0; 6];
                            cf3d_bld_b05_calc_grad(&mut field, &coordinates_flat, &mut gradient)
                                .expect("real F22 kernel gradient evaluates");
                            assert_eq!(
                                energy,
                                if eligible { INCLUDED_ENERGY } else { 0.0 },
                                "mask={parameter_mask}, ignore={ignore_interfragment_interactions}, same={same_component}, relation={relation}, distance={distance}"
                            );
                            assert_eq!(
                                gradient,
                                if eligible {
                                    INCLUDED_GRADIENT
                                } else {
                                    [0.0; 6]
                                },
                                "gradient mask={parameter_mask}, ignore={ignore_interfragment_interactions}, same={same_component}, relation={relation}, distance={distance}"
                            );
                            drop(field);
                            assert_eq!(topology, original_topology);
                            assert_eq!(coordinates, original_coordinates);
                            assert_eq!(rows, original_rows);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn cf3d_frag_f22_preconditions_and_requested_sanitize_precede_pair_skips() {
        let topology = f19_connected_overvalent_topology();
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let coordinates = cf3d_frag_f22_coordinates(&points);
        let molecule_properties = MoleculeProperties::default().with_name("F22-errors");
        let original = (
            topology.clone(),
            coordinates.clone(),
            molecule_properties.clone(),
        );
        let all_missing = vec![None; 4];
        let mut rows = coordinates.conformers_3d[1].coordinates().to_vec();
        let original_rows = rows.clone();

        let error = match cf3d_frag_f22_append(
            &topology,
            &coordinates,
            &molecule_properties,
            &all_missing,
            &mut rows,
            &[],
            2.0,
            true,
        ) {
            Err(error) => error,
            Ok(_) => panic!("source fragment copying sanitizes before nullable-pair skips"),
        };
        assert!(matches!(
            error,
            NonbondedAssemblyError::FragmentMapping(
                NonbondedFragmentMappingError::FragmentCopy(source)
            ) if source.component_index() == Some(0)
        ));
        assert_eq!(
            (
                topology.clone(),
                coordinates.clone(),
                molecule_properties.clone()
            ),
            original,
            "failed fragment preparation does not alter model inputs"
        );
        assert_eq!(
            rows, original_rows,
            "failed preparation leaves kernel rows intact"
        );

        // The source parameter-count precondition comes before its fragment
        // pipeline, so the shorter row must win over the same sanitize error.
        let short_params = vec![None; 3];
        let error = match cf3d_frag_f22_append(
            &topology,
            &coordinates,
            &molecule_properties,
            &short_params,
            &mut rows,
            &[],
            2.0,
            true,
        ) {
            Err(error) => error,
            Ok(_) => panic!("the source parameter-count precondition must fail"),
        };
        assert!(matches!(
            error,
            NonbondedAssemblyError::PairBuilder(UffBuilderError::ParamsLengthMismatch {
                atoms: 4,
                params: 3,
            })
        ));

        // With the flag off, source code never requests sanitized fragments;
        // missing i parameters then skip every matrix read as well.
        let field = cf3d_frag_f22_append(
            &topology,
            &coordinates,
            &molecule_properties,
            &all_missing,
            &mut rows,
            &[],
            2.0,
            false,
        )
        .expect("disabled interfragment filtering performs no F19 or matrix work");
        drop(field);
        assert_eq!(rows, original_rows);
        assert_eq!(
            (topology, coordinates, molecule_properties),
            original,
            "both success and error paths preserve detached input values"
        );
    }

    #[test]
    fn cf3d_frag_f23_parameter_count_precedes_warning_and_conformer_lookup() {
        let topology = topology_with_hydrogen_counts(&[0], &[false]);
        let mut coordinates = f19_coordinates(1);
        let original_coordinates = coordinates.clone();
        let params: [Option<&AtomicParams>; 0] = [];
        let mut diagnostics = Vec::new();

        let error = match prepare_force_field_preamble(
            &topology,
            &mut coordinates,
            999,
            &params,
            &[1],
            &mut diagnostics,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("the source parameter precondition must fail");
            }
        };
        assert_eq!(
            error,
            UffBuilderError::ParamsLengthMismatch {
                atoms: 1,
                params: 0,
            },
            "the source parameter precondition wins over a needed-H warning and missing ID"
        );
        assert!(diagnostics.is_empty());
        assert_eq!(coordinates, original_coordinates);
    }

    #[test]
    fn cf3d_frag_f23_hydrogen_warning_precedes_missing_conformer_error() {
        let params: [Option<&AtomicParams>; 1] = [None];
        let cases = [
            (0, 1, true, "prepared implicit H"),
            (1, 0, true, "explicit H"),
            (0, 0, false, "no source H total"),
        ];

        for (explicit_hydrogens, implicit_hydrogens, warns, label) in cases {
            let topology = topology_with_hydrogen_counts(&[explicit_hydrogens], &[false]);
            let mut coordinates = f19_coordinates(1);
            let original_coordinates = coordinates.clone();
            let mut diagnostics = vec![UffTypingDiagnostic {
                atom_id: Some(AtomId::new(71)),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: "earlier source diagnostic",
            }];

            let error = match prepare_force_field_preamble(
                &topology,
                &mut coordinates,
                999,
                &params,
                &[implicit_hydrogens],
                &mut diagnostics,
            ) {
                Err(error) => error,
                Ok(field) => {
                    drop(field);
                    panic!("the requested 3D conformer ID is absent");
                }
            };
            assert_eq!(
                error,
                UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 999 },
                "{label} retains the source lookup error"
            );
            let mut expected_diagnostics = vec![UffTypingDiagnostic {
                atom_id: Some(AtomId::new(71)),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: "earlier source diagnostic",
            }];
            if warns {
                expected_diagnostics.push(expected_needs_hydrogen_warning());
            }
            assert_eq!(diagnostics, expected_diagnostics, "{label}");
            assert_eq!(coordinates, original_coordinates, "{label}");
        }
    }

    #[test]
    fn cf3d_frag_f23_selects_non_index_3d_id_and_borrows_rows_in_order() {
        let topology = f19_topology(&[Element::C, Element::O], &[]);
        let mut coordinates = f19_coordinates(2);
        coordinates.conformers_2d.push(
            Conformer2D::new(5, vec![[80.0, 81.0], [82.0, 83.0]])
                .with_prop("origin", "2d-id-collision"),
        );
        let original_coordinates = coordinates.clone();
        let selected = &coordinates.conformers_3d[2];
        assert_eq!(selected.id(), 5);
        assert!(
            !selected.is_3d(),
            "the storage collection defines dimension"
        );
        let expected_row_ptrs = selected
            .coordinates()
            .iter()
            .map(|row| row.as_ptr())
            .collect::<Vec<_>>();
        let expected_rows = [[0.0, 20.0, 40.0], [1.0, 21.0, 41.0]];
        let params: [Option<&AtomicParams>; 2] = [None, None];
        let mut diagnostics = Vec::new();

        let field = prepare_force_field_preamble(
            &topology,
            &mut coordinates,
            5,
            &params,
            &[0, 0],
            &mut diagnostics,
        )
        .expect("select ID 5 from the 3D collection, not ID/vector position in 2D");
        assert!(diagnostics.is_empty());
        assert_eq!(field.positions().len(), topology.atoms.len());
        for row in 0..topology.atoms.len() {
            assert_eq!(field.positions()[row].as_ptr(), expected_row_ptrs[row]);
            assert_eq!(&field.positions()[row][..], &expected_rows[row][..]);
        }
        drop(field);
        assert_eq!(coordinates, original_coordinates);
    }

    fn cf3d_frag_f24_coordinates(points: &[[f64; 3]]) -> CoordinateBlock {
        let mut coordinates = CoordinateBlock::default();
        coordinates.conformers_2d.push(
            Conformer2D::new(
                31,
                points
                    .iter()
                    .map(|point| [point[0] + 100.0, point[1]])
                    .collect(),
            )
            .with_prop("origin", "2d-id-collision"),
        );
        coordinates.conformers_3d.push(
            Conformer3D::new(
                12,
                points
                    .iter()
                    .map(|point| [point[0] + 20.0, point[1], point[2]])
                    .collect(),
                false,
            )
            .with_prop("origin", "before-selected"),
        );
        coordinates
            .conformers_3d
            .push(Conformer3D::new(31, points.to_vec(), true).with_prop("origin", "selected"));
        coordinates.conformers_3d.push(
            Conformer3D::new(
                50,
                points
                    .iter()
                    .map(|point| [point[0], point[1] + 20.0, point[2]])
                    .collect(),
                false,
            )
            .with_prop("origin", "after-selected"),
        );
        coordinates.source_coordinate_dim = Some(CoordinateDimension::ThreeD);
        coordinates
    }

    fn cf3d_frag_f24_params(atom_count: usize) -> Vec<Option<AtomicParams>> {
        cf3d_bld_b25_params(atom_count)
            .into_iter()
            .map(|row| {
                row.map(|mut params| {
                    params.x1 = 1.0;
                    params.d1 = 1.0;
                    params.v1 = 1.0;
                    params.u1 = 1.0;
                    params
                })
            })
            .collect()
    }

    macro_rules! cf3d_frag_accept_f24_bond {
        ($end1_idx:literal, $end2_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::BondStretch {
                end1_idx: $end1_idx,
                end2_idx: $end2_idx,
                rest_len: 1.0,
                force_constant: 664.12,
            }
        };
    }

    macro_rules! cf3d_frag_accept_f24_angle {
        ($at1_idx:literal, $at2_idx:literal, $at3_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::AngleBend {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                at3_idx: $at3_idx,
                order: 0,
                force_constant: 352.20281664120762,
                c0: 0.25,
                c1: -6.123233995736766e-17,
                c2: 0.25,
                theta0: 1.5707963267948966,
            }
        };
    }

    macro_rules! cf3d_frag_accept_f24_ordered_angle {
        ($at1_idx:literal, $at2_idx:literal, $at3_idx:literal, $order:literal) => {
            Cf3dFragAcceptContributionIdentity::AngleBend {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                at3_idx: $at3_idx,
                order: $order,
                force_constant: 352.20281664120762,
                c0: 0.0,
                c1: 0.0,
                c2: 0.0,
                theta0: 1.5707963267948966,
            }
        };
    }

    macro_rules! cf3d_frag_accept_f24_ring_angle {
        ($at1_idx:literal, $at2_idx:literal, $at3_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::AngleBend {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                at3_idx: $at3_idx,
                order: 0,
                force_constant: 1162.2100000000003,
                c0: 0.50000000000000011,
                c1: -0.66666666666666685,
                c2: 0.33333333333333337,
                theta0: 1.0471975511965976,
            }
        };
    }

    macro_rules! cf3d_frag_accept_f24_vdw {
        ($at1_idx:literal, $at2_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::Vdw {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                x_ij: 1.0,
                well_depth: 1.0,
                threshold: 10.0,
            }
        };
    }

    macro_rules! cf3d_frag_accept_f24_torsion {
        ($at1_idx:literal, $at2_idx:literal, $at3_idx:literal, $at4_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::TorsionAngle {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                at3_idx: $at3_idx,
                at4_idx: $at4_idx,
                order: 3,
                force_constant: 1.0,
                cos_term: -1.0,
            }
        };
    }

    macro_rules! cf3d_frag_accept_f24_inversion {
        ($at1_idx:literal, $at2_idx:literal, $at3_idx:literal, $at4_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::Inversion {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                at3_idx: $at3_idx,
                at4_idx: $at4_idx,
                force_constant: 2.0,
                c0: 1.0,
                c1: -1.0,
                c2: 0.0,
            }
        };
    }

    macro_rules! cf3d_frag_accept_c3_bond {
        ($end1_idx:literal, $end2_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::BondStretch {
                end1_idx: $end1_idx,
                end2_idx: $end2_idx,
                rest_len: 1.514,
                force_constant: 699.591798712679,
            }
        };
    }

    macro_rules! cf3d_frag_accept_c3_vdw {
        ($at1_idx:literal, $at2_idx:literal) => {
            Cf3dFragAcceptContributionIdentity::Vdw {
                at1_idx: $at1_idx,
                at2_idx: $at2_idx,
                x_ij: 3.851,
                well_depth: 0.105,
                threshold: 38.509999999999998,
            }
        };
    }

    fn cf3d_frag_accept_identities_match(
        observed: &[Cf3dFragAcceptContributionIdentity],
        expected: &[Cf3dFragAcceptContributionIdentity],
    ) -> bool {
        observed == expected
    }

    fn cf3d_frag_accept_assert_fixed_force_field_reference(
        force_field: &mut ForceField<'_>,
        expected_identities: &[Cf3dFragAcceptContributionIdentity],
        expected_term_energies: &[f64],
        expected_total_energy: f64,
        expected_gradient: &[f64],
        label: &str,
    ) {
        assert_eq!(
            crate::kernel::cf3d_frag_f24_contribution_energies(force_field),
            Err(ForceFieldKernelError::NotInitialized),
            "constructor returns the source-uninitialized field: {label}"
        );

        let observed_identities =
            crate::kernel::cf3d_frag_accept_contribution_identities(force_field);
        assert!(
            cf3d_frag_accept_identities_match(&observed_identities, expected_identities),
            "{label}: pinned source contribution identities remain in exact append order; expected {expected_identities:?}, observed {observed_identities:?}"
        );

        force_field
            .initialize()
            .expect("fixed acceptance force field initializes");
        let actual_term_energies = crate::kernel::cf3d_frag_f24_contribution_energies(force_field)
            .expect("fixed acceptance contribution energies evaluate");
        assert_eq!(
            actual_term_energies.len(),
            expected_term_energies.len(),
            "{label}: contribution count"
        );
        for (term_index, (actual, expected)) in actual_term_energies
            .iter()
            .zip(expected_term_energies)
            .enumerate()
        {
            assert!(
                (*actual - *expected).abs() <= 1.0e-10,
                "{label}: source term {term_index} expected {expected:.17e}, got {actual:.17e}"
            );
        }

        let (actual_total_energy, actual_gradient) =
            crate::kernel::cf3d_frag_accept_current_energy_and_gradient(force_field)
                .expect("fixed acceptance total energy and gradient evaluate");
        assert!(
            (actual_total_energy - expected_total_energy).abs() <= 1.0e-10,
            "{label}: total expected {expected_total_energy:.17e}, got {actual_total_energy:.17e}"
        );
        assert_eq!(
            actual_gradient.len(),
            expected_gradient.len(),
            "{label}: complete atom-major gradient length"
        );
        for (component, (actual, expected)) in
            actual_gradient.iter().zip(expected_gradient).enumerate()
        {
            assert!(
                (*actual - *expected).abs() <= 1.0e-10,
                "{label}: gradient component {component} expected {expected:.17e}, got {actual:.17e}"
            );
        }
    }

    fn cf3d_frag_accept_assert_fixed_source_reference(
        topology: TopologyBlock,
        points: &[[f64; 3]],
        params: Vec<Option<AtomicParams>>,
        valence: ValenceAssignment,
        ignore_interfragment_interactions: bool,
        expected_identities: &[Cf3dFragAcceptContributionIdentity],
        expected_term_energies: &[f64],
        expected_total_energy: f64,
        expected_gradient: &[f64],
        label: &str,
    ) {
        let original_topology = topology.clone();
        let mut coordinates = cf3d_frag_f24_coordinates(points);
        let original_coordinates = coordinates.clone();
        let molecule_properties = MoleculeProperties::default().with_name("F24-fixed");
        let original_properties = molecule_properties.clone();
        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed acceptance topology yields source ring information");
        let selected_index = coordinates
            .conformers_3d
            .iter()
            .position(|conformer| conformer.id() == 31)
            .expect("fixed acceptance coordinate block has selected ID 31");
        let selected_row_ptrs = coordinates.conformers_3d[selected_index]
            .coordinates()
            .iter()
            .map(|row| row.as_ptr())
            .collect::<Vec<_>>();
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut diagnostics = Vec::new();

        let mut force_field = construct_force_field_with_params(
            &topology,
            &mut coordinates,
            31,
            &borrowed_params,
            &rings,
            &valence,
            &molecule_properties,
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            ignore_interfragment_interactions,
        )
        .unwrap_or_else(|error| panic!("{label}: fixed source state failed: {error:?}"));
        assert!(
            diagnostics.is_empty(),
            "{label}: no source H warning is due"
        );
        assert_eq!(
            force_field.positions().len(),
            topology.atoms.len(),
            "{label}"
        );
        for (row, expected_ptr) in selected_row_ptrs.iter().enumerate() {
            assert_eq!(
                force_field.positions()[row].as_ptr(),
                *expected_ptr,
                "{label}: selected ID31 row {row} remains borrowed"
            );
        }
        cf3d_frag_accept_assert_fixed_force_field_reference(
            &mut force_field,
            expected_identities,
            expected_term_energies,
            expected_total_energy,
            expected_gradient,
            label,
        );

        drop(force_field);
        assert_eq!(topology, original_topology, "{label}: topology is borrowed");
        assert_eq!(
            coordinates, original_coordinates,
            "{label}: coordinates are unchanged"
        );
        assert_eq!(
            molecule_properties, original_properties,
            "{label}: properties are unchanged"
        );
    }

    #[test]
    fn cf3d_frag_accept_rejects_swapped_equal_energy_bond_identities() {
        let params = cf3d_bld_b15_atomic_params();
        let mut rows = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.5, 0.0, 0.0],
        ];
        let positions = rows
            .iter_mut()
            .map(|point| &mut point[..])
            .collect::<Vec<_>>();
        let first = BondStretchContrib::new(&positions, 0, 1, 1.0, &params, &params)
            .expect("fixed source first single bond contribution");
        let second = BondStretchContrib::new(&positions, 2, 3, 1.0, &params, &params)
            .expect("fixed source second single bond contribution");

        let mut force_field = ForceField::new(3);
        force_field.positions_mut().extend(positions);
        force_field.add_contribution(Box::new(first));
        force_field.add_contribution(Box::new(second));
        force_field
            .initialize()
            .expect("two source bonds initialize");
        let energies = crate::kernel::cf3d_frag_f24_contribution_energies(&mut force_field)
            .expect("fixed source bond energies evaluate");
        assert_eq!(energies.len(), 2);
        assert_eq!(
            energies[0], energies[1],
            "the two source bond terms have equal energy"
        );

        const PINNED_SOURCE_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 2] = [
            cf3d_frag_accept_f24_bond!(0, 1),
            cf3d_frag_accept_f24_bond!(2, 3),
        ];
        let observed = crate::kernel::cf3d_frag_accept_contribution_identities(&force_field);
        assert!(cf3d_frag_accept_identities_match(
            &observed,
            &PINNED_SOURCE_IDENTITIES
        ));

        let mut swapped_identities = observed;
        swapped_identities.swap(0, 1);
        assert!(
            !cf3d_frag_accept_identities_match(&swapped_identities, &PINNED_SOURCE_IDENTITIES),
            "same-category equal-energy bond swaps must fail the fixed identity comparison"
        );
    }

    #[test]
    fn cf3d_frag_accept_acyclic_matches_pinned_source_terms_total_and_gradient() {
        const EXPECTED_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 7] = [
            cf3d_frag_accept_f24_bond!(0, 1),
            cf3d_frag_accept_f24_bond!(1, 2),
            cf3d_frag_accept_f24_bond!(2, 3),
            cf3d_frag_accept_f24_angle!(0, 1, 2),
            cf3d_frag_accept_f24_angle!(1, 2, 3),
            cf3d_frag_accept_f24_vdw!(0, 3),
            cf3d_frag_accept_f24_torsion!(0, 1, 2, 3),
        ];
        const EXPECTED_TERM_ENERGIES: [f64; 7] = [
            0.0,
            0.0,
            56.97248895678014,
            0.0,
            0.0,
            -0.23437499999999989,
            0.1464466094067261,
        ];
        const EXPECTED_TOTAL_ENERGY: f64 = 56.884560566186863;
        const EXPECTED_GRADIENT: [f64; 12] = [
            -0.65625000000002109,
            0.0,
            0.40441017177982119,
            2.1316402421976193e-14,
            6.3165944978342198e-15,
            -1.0606601717798361,
            -1.4999807924141973e-14,
            -195.04657456428001,
            -193.98591439250015,
            0.65625000000001477,
            195.04657456428001,
            194.64216439250018,
        ];

        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let points = cf3d_bld_b25_four_points();
        for ignore_interfragment_interactions in [false, true] {
            cf3d_frag_accept_assert_fixed_source_reference(
                topology.clone(),
                &points,
                cf3d_frag_f24_params(4),
                assignment(&[1, 2, 2, 1], &[0; 4]),
                ignore_interfragment_interactions,
                &EXPECTED_IDENTITIES,
                &EXPECTED_TERM_ENERGIES,
                EXPECTED_TOTAL_ENERGY,
                &EXPECTED_GRADIENT,
                if ignore_interfragment_interactions {
                    "acyclic chain, ignore interfragment interactions"
                } else {
                    "acyclic chain, include interfragment interactions"
                },
            );
        }
    }

    #[test]
    fn cf3d_frag_accept_cyclic_matches_pinned_source_terms_total_and_gradient() {
        const EXPECTED_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 6] = [
            cf3d_frag_accept_f24_bond!(0, 1),
            cf3d_frag_accept_f24_bond!(1, 2),
            cf3d_frag_accept_f24_bond!(2, 0),
            cf3d_frag_accept_f24_ring_angle!(1, 0, 2),
            cf3d_frag_accept_f24_ring_angle!(0, 1, 2),
            cf3d_frag_accept_f24_ring_angle!(1, 2, 0),
        ];
        const EXPECTED_TERM_ENERGIES: [f64; 6] = [
            0.0,
            1.0638450578568244,
            1.0638450578568244,
            0.69727672331441859,
            0.69727672331441859,
            2.9589574127846334,
        ];
        const EXPECTED_TOTAL_ENERGY: f64 = 6.4812009751271198;
        const EXPECTED_GRADIENT: [f64; 9] = [
            -92.885600635316962,
            102.38194407475203,
            0.0,
            92.885600635316962,
            102.38194407475203,
            0.0,
            0.0,
            -204.763888149504,
            0.0,
        ];

        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 3],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]];
        for ignore_interfragment_interactions in [false, true] {
            cf3d_frag_accept_assert_fixed_source_reference(
                topology.clone(),
                &points,
                cf3d_frag_f24_params(3),
                assignment(&[2, 2, 2], &[0; 3]),
                ignore_interfragment_interactions,
                &EXPECTED_IDENTITIES,
                &EXPECTED_TERM_ENERGIES,
                EXPECTED_TOTAL_ENERGY,
                &EXPECTED_GRADIENT,
                if ignore_interfragment_interactions {
                    "cyclic triangle, ignore interfragment interactions"
                } else {
                    "cyclic triangle, include interfragment interactions"
                },
            );
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn cf3d_frag_accept_TBP_matches_pinned_source_terms_total_and_gradient() {
        const EXPECTED_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 15] = [
            cf3d_frag_accept_f24_bond!(0, 1),
            cf3d_frag_accept_f24_bond!(0, 2),
            cf3d_frag_accept_f24_bond!(0, 3),
            cf3d_frag_accept_f24_bond!(0, 4),
            cf3d_frag_accept_f24_bond!(0, 5),
            cf3d_frag_accept_f24_ordered_angle!(1, 0, 5, 2),
            cf3d_frag_accept_f24_ordered_angle!(2, 0, 3, 3),
            cf3d_frag_accept_f24_ordered_angle!(2, 0, 4, 3),
            cf3d_frag_accept_f24_ordered_angle!(3, 0, 4, 3),
            cf3d_frag_accept_f24_angle!(1, 0, 2),
            cf3d_frag_accept_f24_angle!(1, 0, 3),
            cf3d_frag_accept_f24_angle!(1, 0, 4),
            cf3d_frag_accept_f24_angle!(5, 0, 2),
            cf3d_frag_accept_f24_angle!(5, 0, 3),
            cf3d_frag_accept_f24_angle!(5, 0, 4),
        ];
        const EXPECTED_TERM_ENERGIES: [f64; 15] = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            63.396506995417347,
            78.174624112512092,
            5.0316477058288935,
            39.133646293467514,
            63.39650699541739,
            0.0,
            0.0,
            40.573764477067094,
            40.57376447706713,
            22.822742518350243,
        ];
        const EXPECTED_TOTAL_ENERGY: f64 = 353.10320357512768;
        const EXPECTED_GRADIENT: [f64; 18] = [
            -252.1073396762832,
            229.11553980594496,
            312.64940720413813,
            0.0,
            -169.05735198777964,
            -169.0573519877797,
            -190.37210292639884,
            -142.77907719479913,
            -139.04027993483837,
            129.72334142528956,
            -137.89670086939211,
            103.42252565204407,
            146.40366682141726,
            85.371708655802195,
            113.82894487440292,
            166.35243435597519,
            135.24588159022372,
            -221.80324580796696,
        ];

        let topology = cf3d_bld_b21_topology(Hybridization::Sp3d, 5);
        let points = CF3D_BLD_B20_COORDINATES
            .chunks_exact(3)
            .map(|point| [point[0], point[1], point[2]])
            .collect::<Vec<_>>();
        for ignore_interfragment_interactions in [false, true] {
            if ignore_interfragment_interactions {
                // Pinned Builder.cpp calls getMolFrags(mol, true) before the
                // nonbonded pair loop for this flag. The frozen carbon TBP
                // center has explicit valence five, so source sanitation
                // fails here and no field terms or numerical result exist.
                let mut coordinates = cf3d_frag_f24_coordinates(&points);
                let original_coordinates = coordinates.clone();
                let original_topology = topology.clone();
                let molecule_properties = MoleculeProperties::default().with_name("F24-fixed");
                let original_properties = molecule_properties.clone();
                let rings = cosmolkit_core::fast_find_rings(&topology)
                    .expect("fixed TBP topology yields source ring information");
                let params = cf3d_frag_f24_params(6);
                let borrowed_params = cf3d_bld_borrowed_params(&params);
                let valence = assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]);
                let mut diagnostics = Vec::new();
                let result = construct_force_field_with_params(
                    &topology,
                    &mut coordinates,
                    31,
                    &borrowed_params,
                    &rings,
                    &valence,
                    &molecule_properties,
                    &mut diagnostics,
                    DEFAULT_TORSION_BOND_SMARTS,
                    1000.0,
                    true,
                );
                let source_error = match result {
                    Err(ForceFieldConstructionError::Nonbonded(
                        NonbondedAssemblyError::FragmentMapping(
                            NonbondedFragmentMappingError::FragmentCopy(source_error),
                        ),
                    )) => source_error,
                    Err(error) => panic!(
                        "trigonal bipyramid, ignored interfragment interactions: wrong source error: {error:?}"
                    ),
                    Ok(_) => panic!(
                        "trigonal bipyramid, ignored interfragment interactions: source fragment sanitation must fail"
                    ),
                };
                assert_eq!(source_error.component_index(), Some(0));
                let source_error_debug = format!("{source_error:?}");
                for expected_fragment in [
                    "FinalSanitize",
                    "atom: AtomId(0)",
                    "atomic_number: 6",
                    "phase: Explicit",
                    "calculated: Some(5)",
                    "greater than permitted",
                ] {
                    assert!(
                        source_error_debug.contains(expected_fragment),
                        "missing source failure detail {expected_fragment:?} in {source_error_debug}"
                    );
                }
                assert!(diagnostics.is_empty());
                assert_eq!(topology, original_topology);
                assert_eq!(coordinates, original_coordinates);
                assert_eq!(molecule_properties, original_properties);
            } else {
                cf3d_frag_accept_assert_fixed_source_reference(
                    topology.clone(),
                    &points,
                    cf3d_frag_f24_params(6),
                    assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]),
                    false,
                    &EXPECTED_IDENTITIES,
                    &EXPECTED_TERM_ENERGIES,
                    EXPECTED_TOTAL_ENERGY,
                    &EXPECTED_GRADIENT,
                    "trigonal bipyramid, include interfragment interactions",
                );
            }
        }
    }

    #[test]
    fn cf3d_frag_accept_inversion_matches_pinned_source_terms_total_and_gradient() {
        const EXPECTED_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 9] = [
            cf3d_frag_accept_f24_bond!(0, 1),
            cf3d_frag_accept_f24_bond!(0, 2),
            cf3d_frag_accept_f24_bond!(0, 3),
            cf3d_frag_accept_f24_ordered_angle!(1, 0, 2, 3),
            cf3d_frag_accept_f24_ordered_angle!(1, 0, 3, 3),
            cf3d_frag_accept_f24_ordered_angle!(2, 0, 3, 3),
            cf3d_frag_accept_f24_inversion!(1, 0, 2, 3),
            cf3d_frag_accept_f24_inversion!(1, 0, 3, 2),
            cf3d_frag_accept_f24_inversion!(2, 0, 3, 1),
        ];
        const EXPECTED_TERM_ENERGIES: [f64; 9] = [
            0.0,
            0.0,
            0.0,
            39.133646293467514,
            39.133646293467514,
            39.133646293467514,
            2.0,
            2.0,
            2.0,
        ];
        const EXPECTED_TOTAL_ENERGY: f64 = 123.40093888040255;
        const EXPECTED_GRADIENT: [f64; 12] = [
            -234.80187776080507,
            -234.80187776080507,
            -234.80187776080507,
            0.0,
            117.40093888040253,
            117.40093888040253,
            117.40093888040253,
            0.0,
            117.40093888040253,
            117.40093888040253,
            117.40093888040253,
            0.0,
        ];

        let topology = cf3d_bld_b15_topology(
            &[
                Hybridization::Sp2,
                Hybridization::Sp3,
                Hybridization::Sp3,
                Hybridization::Sp3,
            ],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        for ignore_interfragment_interactions in [false, true] {
            cf3d_frag_accept_assert_fixed_source_reference(
                topology.clone(),
                &points,
                cf3d_frag_f24_params(4),
                assignment(&[3, 1, 1, 1], &[0; 4]),
                ignore_interfragment_interactions,
                &EXPECTED_IDENTITIES,
                &EXPECTED_TERM_ENERGIES,
                EXPECTED_TOTAL_ENERGY,
                &EXPECTED_GRADIENT,
                if ignore_interfragment_interactions {
                    "inversion center, ignore interfragment interactions"
                } else {
                    "inversion center, include interfragment interactions"
                },
            );
        }
    }

    #[test]
    fn cf3d_frag_accept_disconnected_matches_pinned_source_terms_total_and_gradient() {
        const INCLUDE_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 6] = [
            cf3d_frag_accept_c3_bond!(0, 1),
            cf3d_frag_accept_c3_bond!(2, 3),
            cf3d_frag_accept_c3_vdw!(0, 2),
            cf3d_frag_accept_c3_vdw!(0, 3),
            cf3d_frag_accept_c3_vdw!(1, 2),
            cf3d_frag_accept_c3_vdw!(1, 3),
        ];
        const INCLUDE_TERM_ENERGIES: [f64; 6] = [
            0.068559996273842655,
            0.068559996273842655,
            -0.039261529774872787,
            -0.0088855879334485757,
            -0.042044400055051313,
            -0.039261529774872787,
        ];
        const INCLUDE_TOTAL_ENERGY: f64 = 0.0076669450094398392;
        const INCLUDE_GRADIENT: [f64; 12] = [
            9.744641055138084,
            0.0,
            0.0,
            -9.341304215726641,
            0.0,
            0.0,
            9.341304215726641,
            0.0,
            0.0,
            -9.744641055138084,
            0.0,
            0.0,
        ];
        const IGNORE_IDENTITIES: [Cf3dFragAcceptContributionIdentity; 2] = [
            cf3d_frag_accept_c3_bond!(0, 1),
            cf3d_frag_accept_c3_bond!(2, 3),
        ];
        const IGNORE_TERM_ENERGIES: [f64; 2] = [0.068559996273842655, 0.068559996273842655];
        const IGNORE_TOTAL_ENERGY: f64 = 0.13711999254768531;
        const IGNORE_GRADIENT: [f64; 12] = [
            9.7942851819775143,
            0.0,
            0.0,
            -9.7942851819775143,
            0.0,
            0.0,
            9.7942851819775143,
            0.0,
            0.0,
            -9.7942851819775143,
            0.0,
            0.0,
        ];

        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)],
        );
        let points = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.5, 0.0, 0.0],
        ];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1; 4], &[0; 4]);
        let total_valences = [4; 4];
        let conjugated = [false; 4];
        let molecule_properties = MoleculeProperties::default().with_name("F24-fixed");
        let default_params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut typing_diagnostics = Vec::new();
        let (typed_params, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugated,
            default_params.as_ref(),
            &mut typing_diagnostics,
        )
        .expect("frozen disconnected fixture uses source UFF typing");
        assert!(found_all);
        assert!(typing_diagnostics.is_empty());

        for ignore_interfragment_interactions in [false, true] {
            let (
                expected_identities,
                expected_term_energies,
                expected_total_energy,
                expected_gradient,
            ): (&[Cf3dFragAcceptContributionIdentity], &[f64], f64, &[f64]) =
                if ignore_interfragment_interactions {
                    (
                        &IGNORE_IDENTITIES,
                        &IGNORE_TERM_ENERGIES,
                        IGNORE_TOTAL_ENERGY,
                        &IGNORE_GRADIENT,
                    )
                } else {
                    (
                        &INCLUDE_IDENTITIES,
                        &INCLUDE_TERM_ENERGIES,
                        INCLUDE_TOTAL_ENERGY,
                        &INCLUDE_GRADIENT,
                    )
                };
            for use_automatic_typing in [false, true] {
                let label = format!(
                    "disconnected dimers, ignore_interfragment_interactions={ignore_interfragment_interactions}, automatic_typing={use_automatic_typing}"
                );
                let mut coordinates = cf3d_frag_f24_coordinates(&points);
                let original_coordinates = coordinates.clone();
                let selected_index = coordinates
                    .conformers_3d
                    .iter()
                    .position(|conformer| conformer.id() == 31)
                    .expect("frozen acceptance coordinate block has selected ID 31");
                let selected_row_ptrs = coordinates.conformers_3d[selected_index]
                    .coordinates()
                    .iter()
                    .map(|row| row.as_ptr())
                    .collect::<Vec<_>>();
                let mut diagnostics = Vec::new();
                let mut force_field = if use_automatic_typing {
                    construct_force_field_with_automatic_typing(
                        &topology,
                        &mut coordinates,
                        31,
                        &total_valences,
                        &conjugated,
                        &rings,
                        &valence,
                        &molecule_properties,
                        &mut diagnostics,
                        DEFAULT_TORSION_BOND_SMARTS,
                        1000.0,
                        ignore_interfragment_interactions,
                    )
                    .unwrap_or_else(|error| panic!("{label}: automatic construction: {error:?}"))
                } else {
                    construct_force_field_with_params(
                        &topology,
                        &mut coordinates,
                        31,
                        &typed_params,
                        &rings,
                        &valence,
                        &molecule_properties,
                        &mut diagnostics,
                        DEFAULT_TORSION_BOND_SMARTS,
                        1000.0,
                        ignore_interfragment_interactions,
                    )
                    .unwrap_or_else(|error| panic!("{label}: supplied construction: {error:?}"))
                };
                assert_eq!(
                    force_field.positions().len(),
                    topology.atoms.len(),
                    "{label}"
                );
                for (row, expected_ptr) in selected_row_ptrs.iter().enumerate() {
                    assert_eq!(
                        force_field.positions()[row].as_ptr(),
                        *expected_ptr,
                        "{label}: selected ID31 row {row} remains borrowed"
                    );
                }
                cf3d_frag_accept_assert_fixed_force_field_reference(
                    &mut force_field,
                    expected_identities,
                    expected_term_energies,
                    expected_total_energy,
                    expected_gradient,
                    &label,
                );
                drop(force_field);
                assert_eq!(
                    coordinates, original_coordinates,
                    "{label}: coordinates unchanged"
                );
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn cf3d_frag_f24_shared_helper_assembly_order_reference<'a>(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        selected_conformer_id: usize,
        params: &UffParamsByAtom<'_>,
        rings: &RingInfo,
        valence: &ValenceAssignment,
        molecule_properties: &MoleculeProperties,
        torsion_bond_smarts: &str,
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
        rows: &'a mut [[f64; 3]],
    ) -> Result<ForceField<'a>, ForceFieldConstructionError> {
        // This test assembly calls the same production stage helpers as the
        // constructor, in the pinned Builder.cpp stage order. It checks
        // assembly wiring and order only; it is not an independent chemistry
        // or numerical reference.
        let mut field = ForceField::new(3);
        attach_positions(&mut field, rows);
        add_bonds(topology, params, &mut field)?;
        add_angles(topology, params, rings, &mut field)?;
        add_angle_special_cases(topology, params, &mut field)?;
        let neighbor_matrix = build_neighbor_matrix(topology)?;
        let selected_index = coordinates
            .conformers_3d
            .iter()
            .position(|conformer| conformer.id() == selected_conformer_id)
            .ok_or(UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                conformer_id: selected_conformer_id,
            })?;
        let selected = &coordinates.conformers_3d[selected_index];
        add_nonbonded(
            topology,
            &coordinates.conformers_2d,
            &coordinates.conformers_3d[..selected_index],
            selected_conformer_id,
            selected.is_3d(),
            selected.props(),
            &coordinates.conformers_3d[selected_index + 1..],
            coordinates.source_coordinate_dim,
            molecule_properties,
            params,
            &mut field,
            &neighbor_matrix,
            vdw_threshold,
            ignore_interfragment_interactions,
        )?;
        add_torsions(
            topology,
            rings,
            valence,
            params,
            torsion_bond_smarts,
            &mut field,
        )?;
        add_inversions(topology, params, &mut field)?;
        Ok(field)
    }

    #[allow(clippy::too_many_arguments)]
    fn cf3d_frag_f24_assert_shared_helper_assembly_order(
        topology: TopologyBlock,
        mut coordinates: CoordinateBlock,
        params: Vec<Option<AtomicParams>>,
        valence: ValenceAssignment,
        ignore_interfragment_interactions: bool,
        expected_contribution_count: usize,
        label: &str,
    ) {
        let original_topology = topology.clone();
        let original_coordinates = coordinates.clone();
        let molecule_properties = MoleculeProperties::default().with_name("F24-fixed");
        let original_properties = molecule_properties.clone();
        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed F24 topology yields source ring information");
        let selected_index = coordinates
            .conformers_3d
            .iter()
            .position(|conformer| conformer.id() == 31)
            .expect("fixed F24 coordinate block has selected ID 31");
        let selected_row_ptrs = coordinates.conformers_3d[selected_index]
            .coordinates()
            .iter()
            .map(|row| row.as_ptr())
            .collect::<Vec<_>>();
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let mut diagnostics = Vec::new();

        let mut actual = construct_force_field_with_params(
            &topology,
            &mut coordinates,
            31,
            &borrowed_params,
            &rings,
            &valence,
            &molecule_properties,
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            ignore_interfragment_interactions,
        )
        .unwrap_or_else(|error| panic!("{label}: fixed source state failed: {error:?}"));
        assert!(
            diagnostics.is_empty(),
            "{label}: no source H warning is due"
        );
        assert_eq!(actual.positions().len(), topology.atoms.len(), "{label}");
        for (row, expected_ptr) in selected_row_ptrs.iter().enumerate() {
            assert_eq!(
                actual.positions()[row].as_ptr(),
                *expected_ptr,
                "{label} row {row}"
            );
        }
        assert_eq!(
            crate::kernel::cf3d_frag_f24_contribution_energies(&mut actual),
            Err(ForceFieldKernelError::NotInitialized),
            "F24 returns the source-uninitialized field: {label}"
        );

        let reference_coordinates = original_coordinates.clone();
        let mut reference_rows = reference_coordinates.conformers_3d[selected_index]
            .coordinates()
            .to_vec();
        let mut shared_helper_reference = cf3d_frag_f24_shared_helper_assembly_order_reference(
            &topology,
            &reference_coordinates,
            31,
            &borrowed_params,
            &rings,
            &valence,
            &molecule_properties,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            ignore_interfragment_interactions,
            &mut reference_rows,
        )
        .unwrap_or_else(|error| {
            panic!("{label}: shared-helper assembly-order reference failed: {error:?}")
        });

        actual
            .initialize()
            .expect("actual F24 kernel initializes in test");
        shared_helper_reference
            .initialize()
            .expect("shared-helper assembly-order F24 kernel initializes in test");
        let actual_energies = crate::kernel::cf3d_frag_f24_contribution_energies(&mut actual)
            .expect("actual F24 term energies evaluate");
        let reference_energies =
            crate::kernel::cf3d_frag_f24_contribution_energies(&mut shared_helper_reference)
                .expect("shared-helper assembly-order term energies evaluate");
        assert_eq!(
            actual_energies.len(),
            expected_contribution_count,
            "{label}"
        );
        assert_eq!(
            reference_energies.len(),
            expected_contribution_count,
            "{label}"
        );
        assert_eq!(
            actual_energies
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            reference_energies
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            "F24 and its shared-helper assembly-order reference must preserve the same category and append order: {label}"
        );

        drop(actual);
        drop(shared_helper_reference);
        assert_eq!(topology, original_topology, "{label}: topology is borrowed");
        assert_eq!(
            coordinates, original_coordinates,
            "{label}: coordinates are unchanged"
        );
        assert_eq!(
            reference_rows,
            original_coordinates.conformers_3d[selected_index].coordinates()
        );
        assert_eq!(
            molecule_properties, original_properties,
            "{label}: properties are unchanged"
        );
    }

    #[test]
    fn cf3d_frag_f24_acyclic_cyclic_tbp_and_inversion_share_assembly_order() {
        let acyclic = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let acyclic_points = cf3d_bld_b25_four_points();
        cf3d_frag_f24_assert_shared_helper_assembly_order(
            acyclic,
            cf3d_frag_f24_coordinates(&acyclic_points),
            cf3d_frag_f24_params(4),
            assignment(&[1, 2, 2, 1], &[0; 4]),
            false,
            7,
            "acyclic chain: three bonds, two angles, one 1-4 vdW pair, one torsion",
        );

        let cyclic = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 3],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        let cyclic_points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]];
        cf3d_frag_f24_assert_shared_helper_assembly_order(
            cyclic,
            cf3d_frag_f24_coordinates(&cyclic_points),
            cf3d_frag_f24_params(3),
            assignment(&[2, 2, 2], &[0; 3]),
            false,
            6,
            "three-membered ring: three bonds and three source ring angles",
        );

        let tbp = cf3d_bld_b21_topology(Hybridization::Sp3d, 5);
        let tbp_points = CF3D_BLD_B20_COORDINATES
            .chunks_exact(3)
            .map(|point| [point[0], point[1], point[2]])
            .collect::<Vec<_>>();
        cf3d_frag_f24_assert_shared_helper_assembly_order(
            tbp,
            cf3d_frag_f24_coordinates(&tbp_points),
            cf3d_frag_f24_params(6),
            assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]),
            false,
            15,
            "trigonal bipyramid: five bonds and ten special angles",
        );

        let inversion = cf3d_bld_b15_topology(
            &[
                Hybridization::Sp2,
                Hybridization::Sp3,
                Hybridization::Sp3,
                Hybridization::Sp3,
            ],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let inversion_points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        cf3d_frag_f24_assert_shared_helper_assembly_order(
            inversion,
            cf3d_frag_f24_coordinates(&inversion_points),
            cf3d_frag_f24_params(4),
            assignment(&[3, 1, 1, 1], &[0; 4]),
            false,
            9,
            "sp2 carbon inversion center: three bonds, three angles, three permutations",
        );
    }

    #[test]
    fn cf3d_frag_f24_interfragment_flags_and_nullable_rows_keep_source_guards() {
        let topology = topology_with_bonds(4, &[(0, 1, false), (2, 3, false)]);
        let points = [
            [0.0, 0.0, 0.0],
            [0.4, 0.0, 0.0],
            [0.8, 0.0, 0.0],
            [1.2, 0.0, 0.0],
        ];
        for (ignore_interfragment_interactions, expected_count) in [(false, 6), (true, 2)] {
            cf3d_frag_f24_assert_shared_helper_assembly_order(
                topology.clone(),
                cf3d_frag_f24_coordinates(&points),
                cf3d_frag_f24_params(4),
                assignment(&[1, 1, 1, 1], &[0; 4]),
                ignore_interfragment_interactions,
                expected_count,
                if ignore_interfragment_interactions {
                    "fragment filtering retains only two bond terms"
                } else {
                    "unfiltered two dimers append four cross-fragment vdW terms"
                },
            );
        }

        let chain = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let chain_points = cf3d_bld_b25_four_points();
        let mut partial = cf3d_frag_f24_params(4);
        partial[1] = None;
        partial[2] = None;
        cf3d_frag_f24_assert_shared_helper_assembly_order(
            chain.clone(),
            cf3d_frag_f24_coordinates(&chain_points),
            partial,
            assignment(&[1, 2, 2, 1], &[0; 4]),
            false,
            1,
            "only the source-eligible atom 0-3 vdW pair survives missing rows",
        );
        cf3d_frag_f24_assert_shared_helper_assembly_order(
            chain,
            cf3d_frag_f24_coordinates(&chain_points),
            vec![None; 4],
            assignment(&[1, 2, 2, 1], &[0; 4]),
            false,
            0,
            "all-null parameters skip every contribution without skipping construction",
        );
    }

    #[test]
    fn cf3d_frag_f24_preserves_first_stage_errors_and_coordinate_preconditions() {
        let one_h_topology = topology_with_hydrogen_counts(&[1], &[false]);
        let one_param = [Some(atomic_params(0.5, 1.0))];
        let borrowed_one_param = cf3d_bld_borrowed_params(&one_param);
        let mut missing_coordinates = CoordinateBlock::default();
        let mut diagnostics = Vec::new();
        let missing_id = match construct_force_field_with_params(
            &one_h_topology,
            &mut missing_coordinates,
            700,
            &borrowed_one_param,
            &RingInfo::new(cosmolkit_core::RingFindType::Fast, 1, 0),
            &assignment(&[0], &[0]),
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            2.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("missing explicit conformer ID must fail after the source H warning");
            }
        };
        assert!(matches!(
            missing_id,
            ForceFieldConstructionError::Builder(
                UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 700 }
            )
        ));
        assert_eq!(diagnostics, vec![expected_needs_hydrogen_warning()]);

        let two_atom_topology = topology_with_bond_orders(2, &[(0, 1, BondOrder::Zero)]);
        let two_params = cf3d_frag_f24_params(2);
        let borrowed_two_params = cf3d_bld_borrowed_params(&two_params);
        let mut malformed_bond_coordinates =
            cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let malformed_bond_error = match construct_force_field_with_params(
            &two_atom_topology,
            &mut malformed_bond_coordinates,
            31,
            &borrowed_two_params,
            &cf3d_bld_b15_ring_info(&two_atom_topology),
            &assignment(&[1, 1], &[0, 0]),
            &MoleculeProperties::default(),
            &mut Vec::new(),
            "[",
            2.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("the first source bond constructor must reject zero order");
            }
        };
        assert!(matches!(
            malformed_bond_error,
            ForceFieldConstructionError::Builder(UffBuilderError::ForceFieldKernel(
                ForceFieldKernelError::BadBondOrder
            ))
        ));

        let two_atom_topology = topology_with_bonds(2, &[(0, 1, false)]);
        let mut short_coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0]]);
        let coordinate_error = match construct_force_field_with_params(
            &two_atom_topology,
            &mut short_coordinates,
            31,
            &borrowed_two_params,
            &cf3d_bld_b15_ring_info(&two_atom_topology),
            &assignment(&[1, 1], &[0, 0]),
            &MoleculeProperties::default(),
            &mut Vec::new(),
            "[",
            2.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("source getAtomPos checks the selected row count before addBonds");
            }
        };
        assert!(matches!(
            coordinate_error,
            ForceFieldConstructionError::Builder(
                UffBuilderError::SelectedConformerCoordinateCountMismatch {
                    conformer_id: 31,
                    atoms: 2,
                    coordinates: 1,
                }
            )
        ));

        let invalid_fragment_topology = f19_connected_overvalent_topology();
        let mut invalid_fragment_coordinates = f19_coordinates(4);
        let invalid_fragment_params: [Option<&AtomicParams>; 4] = [None; 4];
        let invalid_fragment_valence = assignment(&[0; 4], &[0; 4]);
        let invalid_fragment_error = match construct_force_field_with_params(
            &invalid_fragment_topology,
            &mut invalid_fragment_coordinates,
            11,
            &invalid_fragment_params,
            &cf3d_bld_b15_ring_info(&invalid_fragment_topology),
            &invalid_fragment_valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            "[",
            2.0,
            true,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("source fragment sanitization precedes torsion query parsing");
            }
        };
        assert!(matches!(
            invalid_fragment_error,
            ForceFieldConstructionError::Nonbonded(
                NonbondedAssemblyError::FragmentMapping(
                    NonbondedFragmentMappingError::FragmentCopy(source)
                )
            ) if source.component_index() == Some(0)
        ));
    }

    #[test]
    fn cf3d_frag_f25_automatic_matches_supplied_stage_results_and_selected_rows() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)],
        );
        let points = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.5, 0.0, 0.0],
        ];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1, 1, 1], &[0; 4]);
        let total_valences = [4; 4];
        let conjugated = [false; 4];
        let molecule_properties = MoleculeProperties::default();
        let default_params = ParamCollection::get_params("").expect("pinned default UFF table");

        for (ignore_interfragment_interactions, expected_term_count) in [(false, 6), (true, 2)] {
            let mut automatic_coordinates = cf3d_frag_f24_coordinates(&points);
            let original_automatic_coordinates = automatic_coordinates.clone();
            let automatic_row_ptrs = automatic_coordinates.conformers_3d[1]
                .coordinates()
                .iter()
                .map(|row| row.as_ptr())
                .collect::<Vec<_>>();
            let mut automatic_diagnostics = Vec::new();
            let mut automatic = construct_force_field_with_automatic_typing(
                &topology,
                &mut automatic_coordinates,
                31,
                &total_valences,
                &conjugated,
                &rings,
                &valence,
                &molecule_properties,
                &mut automatic_diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                1000.0,
                ignore_interfragment_interactions,
            )
            .expect("all four source C_3 rows are available");
            assert!(automatic_diagnostics.is_empty());
            assert_eq!(automatic.positions().len(), topology.atoms.len());
            for (row, expected_ptr) in automatic_row_ptrs.iter().enumerate() {
                assert_eq!(automatic.positions()[row].as_ptr(), *expected_ptr);
            }

            let mut expected_typing_diagnostics = Vec::new();
            let (typed_params, found_all) = get_atom_types(
                &topology,
                &total_valences,
                &conjugated,
                default_params.as_ref(),
                &mut expected_typing_diagnostics,
            )
            .expect("the same source-prepared state yields parameter slots");
            assert!(found_all);
            assert!(expected_typing_diagnostics.is_empty());

            let mut supplied_coordinates = cf3d_frag_f24_coordinates(&points);
            let original_supplied_coordinates = supplied_coordinates.clone();
            let supplied_row_ptrs = supplied_coordinates.conformers_3d[1]
                .coordinates()
                .iter()
                .map(|row| row.as_ptr())
                .collect::<Vec<_>>();
            let mut supplied_diagnostics = Vec::new();
            let mut supplied = construct_force_field_with_params(
                &topology,
                &mut supplied_coordinates,
                31,
                &typed_params,
                &rings,
                &valence,
                &molecule_properties,
                &mut supplied_diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                1000.0,
                ignore_interfragment_interactions,
            )
            .expect("the source-typed rows construct through F24");
            assert!(supplied_diagnostics.is_empty());
            for (row, expected_ptr) in supplied_row_ptrs.iter().enumerate() {
                assert_eq!(supplied.positions()[row].as_ptr(), *expected_ptr);
            }

            automatic.initialize().expect("automatic field initializes");
            supplied.initialize().expect("supplied field initializes");
            let automatic_terms =
                crate::kernel::cf3d_frag_f24_contribution_energies(&mut automatic)
                    .expect("automatic term energies evaluate");
            let supplied_terms = crate::kernel::cf3d_frag_f24_contribution_energies(&mut supplied)
                .expect("supplied term energies evaluate");
            assert_eq!(automatic_terms.len(), expected_term_count);
            assert_eq!(supplied_terms.len(), expected_term_count);
            assert_eq!(
                automatic_terms
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                supplied_terms
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                "automatic typing must delegate the same source-ordered terms"
            );

            let flat_coordinates = points
                .iter()
                .flat_map(|point| point.iter().copied())
                .collect::<Vec<_>>();
            let automatic_energy = cf3d_bld_b05_calc_energy(&mut automatic, &flat_coordinates)
                .expect("automatic total energy evaluates");
            let supplied_energy = cf3d_bld_b05_calc_energy(&mut supplied, &flat_coordinates)
                .expect("supplied total energy evaluates");
            assert_eq!(automatic_energy.to_bits(), supplied_energy.to_bits());
            let mut automatic_gradient = vec![0.0; flat_coordinates.len()];
            let mut supplied_gradient = vec![0.0; flat_coordinates.len()];
            cf3d_bld_b05_calc_grad(&mut automatic, &flat_coordinates, &mut automatic_gradient)
                .expect("automatic gradient evaluates");
            cf3d_bld_b05_calc_grad(&mut supplied, &flat_coordinates, &mut supplied_gradient)
                .expect("supplied gradient evaluates");
            assert_eq!(
                automatic_gradient
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                supplied_gradient
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>()
            );

            drop(automatic);
            drop(supplied);
            assert_eq!(automatic_coordinates, original_automatic_coordinates);
            assert_eq!(supplied_coordinates, original_supplied_coordinates);
        }
    }

    #[test]
    fn uff_prepare_p05_cached_and_supplied_builder_keep_source_energy_and_error_order() {
        let ion_atom = |row, atomic_number, charge| {
            let element = Element::from_atomic_number(atomic_number)
                .expect("fixed W06 ion element exists in the source table");
            Atom::from_spec(
                AtomId::new(row),
                AtomSpec::new(element)
                    .with_hybridization(Hybridization::Unspecified)
                    .with_formal_charge(charge)
                    .with_no_implicit(true),
            )
        };
        let w06_topology = TopologyBlock::try_from_parts(
            vec![ion_atom(0, 11, 1), ion_atom(1, 17, -1)],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed W06 disconnected-ion topology is structurally valid");
        let w06_valence = assignment(&[0, 0], &[0, 0]);
        let w06_total_valences = [0, 0];
        let w06_conjugation = [false; 2];
        let w06_states = [
            UffAtomStateRef::cached(&w06_topology, &w06_valence)
                .expect("fixed W06 cached rows validate"),
            UffAtomStateRef::supplied_rows(&w06_topology, &w06_total_valences, &w06_conjugation)
                .expect("fixed W06 supplied rows validate"),
        ];
        let rings = cf3d_bld_b15_ring_info(&w06_topology);
        let properties = MoleculeProperties::default();
        let minimum = (2.983_f64 * 3.947).sqrt();
        let well_depth = (0.03_f64 * 0.227).sqrt();
        let ratio = minimum / 4.0;
        let ratio3 = (ratio * ratio) * ratio;
        let ratio6 = ratio3 * ratio3;
        let ratio12 = ratio6 * ratio6;
        let expected_source_energy = well_depth * (ratio12 - 2.0 * ratio6);
        let points = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let flat_points = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        let mut observed_source_energies = [0.0; 2];

        for (state_index, typing_state) in w06_states.into_iter().enumerate() {
            let mut coordinates = cf3d_frag_f24_coordinates(&points);
            let original_coordinates = coordinates.clone();
            let mut diagnostics = Vec::new();
            let mut field = super::construct_force_field_with_automatic_typing(
                &w06_topology,
                &mut coordinates,
                31,
                typing_state,
                &rings,
                &w06_valence,
                &properties,
                &mut diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                10.0,
                false,
            )
            .expect("both source-shaped typing inputs resolve Na+1 and Cl-1");
            assert!(diagnostics.is_empty());
            assert_eq!(field.positions().len(), 2);
            field
                .initialize()
                .expect("fixed W06 source field initializes");
            let term_energies =
                cf3d_frag_f24_contribution_energies(&mut field).expect("one W06 pair evaluates");
            assert_eq!(term_energies.len(), 1);
            assert!(
                (term_energies[0] - expected_source_energy).abs() <= 1.0e-12,
                "fixed W06 source VDW term for cached/supplied input {state_index}"
            );
            let total_energy = cf3d_bld_b05_calc_energy(&mut field, &flat_points)
                .expect("fixed W06 source total energy evaluates");
            assert!(
                (total_energy - expected_source_energy).abs() <= 1.0e-12,
                "fixed W06 source total energy for cached/supplied input {state_index}"
            );
            observed_source_energies[state_index] = total_energy;
            drop(field);
            assert_eq!(coordinates, original_coordinates);
        }
        assert_eq!(
            observed_source_energies[0].to_bits(),
            observed_source_energies[1].to_bits(),
            "cached and supplied input preserve identical source parameters and energy"
        );

        let missing_element =
            Element::from_atomic_number(0).expect("fixed source dummy element is modeled");
        let missing_atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(missing_element)
                .with_hybridization(Hybridization::Unspecified)
                .with_prop("dummyLabel", "f25-unrecognized")
                .expect("fixed source dummy label is valid"),
        );
        let carbon = Element::C;
        let first_carbon = Atom::from_spec(
            AtomId::new(1),
            AtomSpec::new(carbon).with_hybridization(Hybridization::Sp3),
        );
        let second_carbon = Atom::from_spec(
            AtomId::new(2),
            AtomSpec::new(carbon).with_hybridization(Hybridization::Sp3),
        );
        let invalid_bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Zero),
        );
        let error_topology = TopologyBlock::try_from_parts(
            vec![missing_atom, first_carbon, second_carbon],
            vec![invalid_bond],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed F25 source error topology is structurally valid");
        let error_valence = assignment(&[0, 1, 1], &[1, 0, 0]);
        let supplied_error_totals = [0, 1, 1];
        let supplied_error_conjugation = [false; 3];
        let error_states = [
            UffAtomStateRef::cached(&error_topology, &error_valence)
                .expect("fixed F25 cached error rows validate"),
            UffAtomStateRef::supplied_rows(
                &error_topology,
                &supplied_error_totals,
                &supplied_error_conjugation,
            )
            .expect("fixed F25 supplied error rows validate"),
        ];
        let error_rings = cf3d_bld_b15_ring_info(&error_topology);
        let error_points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.5, 0.0, 0.0]];

        for typing_state in error_states {
            let mut coordinates = cf3d_frag_f24_coordinates(&error_points);
            let original_coordinates = coordinates.clone();
            let mut diagnostics = Vec::new();
            let error = match super::construct_force_field_with_automatic_typing(
                &error_topology,
                &mut coordinates,
                31,
                typing_state,
                &error_rings,
                &error_valence,
                &properties,
                &mut diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                10.0,
                false,
            ) {
                Err(error) => error,
                Ok(field) => {
                    drop(field);
                    panic!("source zero-order bond must fail after atom typing");
                }
            };
            assert!(matches!(
                error,
                AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(UffBuilderError::ForceFieldKernel(
                        ForceFieldKernelError::BadBondOrder
                    ))
                )
            ));
            assert_eq!(
                diagnostics,
                [
                    UffTypingDiagnostic {
                        atom_id: Some(AtomId::new(0)),
                        kind: UffTypingDiagnosticKind::Error,
                        message_prefix: "UFFTYPER: Unrecognized atom type: ",
                    },
                    expected_needs_hydrogen_warning(),
                ],
                "cached and supplied routes retain typing diagnostics before the builder error"
            );
            assert_eq!(coordinates, original_coordinates);
        }
    }

    #[test]
    fn cf3d_frag_f25_missing_type_is_nonfatal_and_diagnostics_keep_source_order() {
        let missing_element =
            Element::from_atomic_number(0).expect("the source dummy element is modeled");
        let missing_atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(missing_element)
                .with_hybridization(Hybridization::Unspecified)
                .with_prop("dummyLabel", "f25-unrecognized")
                .expect("fixed dummy label is valid"),
        );
        let carbon = Element::C;
        let atom1 = Atom::from_spec(
            AtomId::new(1),
            AtomSpec::new(carbon).with_hybridization(Hybridization::Sp3),
        );
        let atom2 = Atom::from_spec(
            AtomId::new(2),
            AtomSpec::new(carbon).with_hybridization(Hybridization::Sp3),
        );
        let invalid_bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Zero),
        );
        let topology = TopologyBlock::try_from_parts(
            vec![missing_atom, atom1, atom2],
            vec![invalid_bond],
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed F25 typed-error topology is structurally valid");
        let total_valences = [0, 4, 4];
        let conjugated = [false; 3];
        let params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut source_typing_diagnostics = Vec::new();
        let (source_slots, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugated,
            params.as_ref(),
            &mut source_typing_diagnostics,
        )
        .expect("a missing table key is a nullable row, not a typing exception");
        assert!(!found_all);
        assert_eq!(source_slots.len(), 3);
        assert!(source_slots[0].is_none());
        assert!(source_slots[1].is_some());
        assert!(source_slots[2].is_some());
        assert_eq!(
            source_typing_diagnostics,
            [UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: "UFFTYPER: Unrecognized atom type: ",
            }],
        );

        let mut coordinates =
            cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.5, 0.0, 0.0]]);
        let original_coordinates = coordinates.clone();
        let valence = assignment(&[0, 1, 1], &[1, 0, 0]);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let mut diagnostics = Vec::new();
        let error = match construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("source addBonds must report the supported zero-order bond error");
            }
        };
        assert!(matches!(
            error,
            AutomaticForceFieldConstructionError::Construction(
                ForceFieldConstructionError::Builder(UffBuilderError::ForceFieldKernel(
                    ForceFieldKernelError::BadBondOrder
                ))
            )
        ));
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(0)),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: "UFFTYPER: Unrecognized atom type: ",
                },
                expected_needs_hydrogen_warning(),
            ],
            "typing error logs precede the delegated needs-H warning and builder error"
        );
        assert_eq!(coordinates, original_coordinates);
    }

    #[test]
    fn cf3d_frag_f25_uses_prepared_total_valence_without_recomputation() {
        let magnesium = Element::from_atomic_number(12).expect("magnesium is modeled");
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(magnesium)
                .with_hybridization(Hybridization::Sp3)
                .with_explicit_hydrogens(1),
        );
        let topology =
            TopologyBlock::try_from_parts(vec![atom], Vec::new(), Vec::new(), Vec::new())
                .expect("fixed magnesium atom topology is structurally valid");
        let total_valences = [0];
        let conjugated = [false];
        let mut coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0]]);
        let original_coordinates = coordinates.clone();
        let valence = assignment(&[2], &[0]);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let mut diagnostics = Vec::new();
        let mut field = construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        )
        .expect("the prepared Mg3+2 row exists in the fixed table");
        assert_eq!(
            diagnostics,
            [
                UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(0)),
                    kind: UffTypingDiagnosticKind::Error,
                    message_prefix: "UFFTYPER: Unrecognized charge state for atom: ",
                },
                expected_needs_hydrogen_warning(),
            ],
            "the supplied total-valence row controls typing before needs-H diagnostics"
        );
        assert_eq!(field.positions().len(), 1);
        assert_eq!(
            crate::kernel::cf3d_frag_f24_contribution_energies(&mut field),
            Err(ForceFieldKernelError::NotInitialized),
            "the automatic overload delegates without initializing the field"
        );
        drop(field);
        assert_eq!(coordinates, original_coordinates);
    }

    fn cf3d_frag_integration_group(component: usize, fragment_local: bool) -> SubstanceGroup {
        let source_atom = component * 2;
        let atom_start = if fragment_local { 0 } else { source_atom };
        let bond_row = if fragment_local { 0 } else { component };
        let group_id = if fragment_local { 0 } else { component };
        let bond_id = BondId::new(bond_row);
        SubstanceGroup::new(
            SubstanceGroupId::new(group_id),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_rdkit_sequence_id(70 + component as u32)
        .with_external_id(90 + component as u32)
        .with_atoms(vec![AtomId::new(atom_start), AtomId::new(atom_start + 1)])
        .with_bonds(vec![bond_id])
        .with_bond_role(bond_id, SGroupBondRole::Contained)
        .with_head_crossing_bonds(vec![bond_id])
        .with_crossing_bond_correspondence(vec![bond_id])
        .with_parent_atoms(vec![AtomId::new(atom_start + 1)])
        .with_label(format!("UFF-FRAG integration group {component}"))
        .with_connection(SGroupConnection::HeadToTail)
        .with_subtype("SUP")
        .with_bracket_style(SGroupBracketStyle::Bracket)
        .with_expansion_state("expanded")
        .with_class(format!("class-{component}"))
        .with_component_number(100 + component as u32)
        .with_display(SGroupDisplay {
            brackets: vec![SGroupBracket::new([
                [component as f64, 0.0, 0.0],
                [component as f64 + 1.0, 0.0, 0.0],
                [component as f64, 1.0, 0.0],
            ])],
            field_position: Some([component as f64 + 0.25, 0.75]),
            display_tag: Some(format!("display-{component}")),
        })
        .with_data(SGroupData {
            field_name: Some(format!("field-{component}")),
            field_type: Some("STRING".to_owned()),
            field_info: Some(format!("info-{component}")),
            field_display: Some("plain".to_owned()),
            units: Some("unit".to_owned()),
            query_type: Some("query".to_owned()),
            query_op: Some("equal".to_owned()),
            values: vec![format!("value-{component}")],
        })
        .with_attach_points(vec![SGroupAttachPoint {
            atom: AtomId::new(atom_start),
            leaving_atom: Some(AtomId::new(atom_start + 1)),
            label: Some(format!("attach-{component}")),
            order: Some(1),
        }])
        .with_cstates(vec![SGroupCState::new(
            bond_id,
            [component as f64 + 0.5, 1.25, -2.0],
        )])
        .with_prop("sgroup-property", format!("sgroup-{component}"))
        .with_data_field(format!("data-field-{component}"))
    }

    fn cf3d_frag_integration_stereo_group(component: usize, fragment_local: bool) -> StereoGroup {
        let atom_index = if fragment_local { 0 } else { component * 2 };
        let bond_index = if fragment_local { 0 } else { component };
        StereoGroup::new(
            if component == 0 {
                StereoGroupKind::Or
            } else {
                StereoGroupKind::And
            },
            vec![AtomId::new(atom_index)],
            vec![BondId::new(bond_index)],
        )
        .with_id(110 + component as u32)
        .with_write_id(210 + component as u32)
    }

    fn cf3d_frag_integration_metadata_source()
    -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
        let mut atoms = Vec::with_capacity(4);
        for row in 0..4 {
            let spec = AtomSpec::new(Element::C)
                .with_hybridization(Hybridization::Sp3)
                .with_prop("source-atom", format!("atom-{row}"))
                .expect("fixed atom metadata key is valid")
                .with_computed_prop("computed-atom", format!("atom-cache-{row}"))
                .expect("fixed computed atom metadata key is valid");
            let mut atom = Atom::from_spec(AtomId::new(row), spec);
            atom.set_temporary_flags(0xF0F0_0000_0000_0000 | row as u64);
            atoms.push(atom);
        }

        let mut bonds = Vec::with_capacity(2);
        for component in 0..2 {
            let spec = BondSpec::new(
                AtomId::new(component * 2),
                AtomId::new(component * 2 + 1),
                BondOrder::Single,
            )
            .with_prop("source-bond", format!("bond-{component}"))
            .expect("fixed bond metadata key is valid")
            .with_computed_prop("computed-bond", format!("bond-cache-{component}"))
            .expect("fixed computed bond metadata key is valid");
            let mut bond = Bond::from_spec(BondId::new(component), spec);
            bond.set_temporary_flags(0x0F0F_0000_0000_0000 | component as u64);
            bonds.push(bond);
        }

        let groups = (0..2)
            .map(|component| cf3d_frag_integration_group(component, false))
            .collect();
        let stereo_groups = (0..2)
            .map(|component| cf3d_frag_integration_stereo_group(component, false))
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, groups, stereo_groups)
            .expect("fixed metadata-rich two-dimer topology is valid");
        let points = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.5, 0.0, 0.0],
        ];
        let coordinates = cf3d_frag_f24_coordinates(&points);
        let properties = MoleculeProperties::default()
            .with_name("UFF-FRAG integration source")
            .with_prop("ordinary-molecule", "retained")
            .expect("fixed molecule metadata key is valid")
            .with_computed_prop("computed-molecule", "clear-after-sanitize")
            .expect("fixed computed molecule metadata key is valid")
            .with_sdf_data_field("SOURCE", "integration")
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "atom-values",
                [
                    Some(cosmolkit_model::PropertyValue::String("atom-0".to_owned())),
                    Some(cosmolkit_model::PropertyValue::String("atom-1".to_owned())),
                    Some(cosmolkit_model::PropertyValue::String("atom-2".to_owned())),
                    Some(cosmolkit_model::PropertyValue::String("atom-3".to_owned())),
                ]
                .to_vec(),
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "bond-values",
                [
                    Some(cosmolkit_model::PropertyValue::String("bond-0".to_owned())),
                    Some(cosmolkit_model::PropertyValue::String("bond-1".to_owned())),
                ]
                .to_vec(),
            ));
        (topology, coordinates, properties)
    }

    #[test]
    fn cf3d_frag_integration_both_constructors_keep_source_order_and_kernel_results() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp3; 4],
            &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)],
        );
        let points = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [6.5, 0.0, 0.0],
        ];
        let molecule_properties = MoleculeProperties::default().with_name("integration-fixed");
        let valence = assignment(&[1; 4], &[0; 4]);
        let total_valences = [4; 4];
        let conjugated = [false; 4];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let default_params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut typing_diagnostics = Vec::new();
        let (typed_params, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugated,
            default_params.as_ref(),
            &mut typing_diagnostics,
        )
        .expect("the fixed prepared state types in source atom order");
        assert!(found_all);
        assert!(typing_diagnostics.is_empty());

        let flat_coordinates = points
            .iter()
            .flat_map(|point| point.iter().copied())
            .collect::<Vec<_>>();
        for (ignore_interfragment_interactions, expected_count) in [(false, 6), (true, 2)] {
            let mut automatic_coordinates = cf3d_frag_f24_coordinates(&points);
            let original_automatic_coordinates = automatic_coordinates.clone();
            let selected_row_ptrs = automatic_coordinates.conformers_3d[1]
                .coordinates()
                .iter()
                .map(|row| row.as_ptr())
                .collect::<Vec<_>>();
            let mut automatic_diagnostics = Vec::new();
            let mut automatic = construct_force_field_with_automatic_typing(
                &topology,
                &mut automatic_coordinates,
                31,
                &total_valences,
                &conjugated,
                &rings,
                &valence,
                &molecule_properties,
                &mut automatic_diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                1000.0,
                ignore_interfragment_interactions,
            )
            .expect("automatic constructor completes for shared-helper assembly-order comparison");
            assert!(automatic_diagnostics.is_empty());
            for (row, expected_ptr) in selected_row_ptrs.iter().enumerate() {
                assert_eq!(automatic.positions()[row].as_ptr(), *expected_ptr);
            }

            let mut supplied_coordinates = cf3d_frag_f24_coordinates(&points);
            let original_supplied_coordinates = supplied_coordinates.clone();
            let mut supplied_diagnostics = Vec::new();
            let mut supplied = construct_force_field_with_params(
                &topology,
                &mut supplied_coordinates,
                31,
                &typed_params,
                &rings,
                &valence,
                &molecule_properties,
                &mut supplied_diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                1000.0,
                ignore_interfragment_interactions,
            )
            .expect("supplied constructor completes for shared-helper assembly-order comparison");
            assert!(supplied_diagnostics.is_empty());

            let reference_coordinates = cf3d_frag_f24_coordinates(&points);
            let mut reference_rows = points;
            let mut shared_helper_reference = cf3d_frag_f24_shared_helper_assembly_order_reference(
                &topology,
                &reference_coordinates,
                31,
                &typed_params,
                &rings,
                &valence,
                &molecule_properties,
                DEFAULT_TORSION_BOND_SMARTS,
                1000.0,
                ignore_interfragment_interactions,
                &mut reference_rows,
            )
            .expect("existing shared-helper assembly-order reference constructs");

            automatic
                .initialize()
                .expect("automatic field initializes with the real kernel");
            supplied
                .initialize()
                .expect("supplied field initializes with the real kernel");
            shared_helper_reference.initialize().expect(
                "fixed shared-helper assembly-order stages initialize with the real kernel",
            );
            let automatic_terms =
                crate::kernel::cf3d_frag_f24_contribution_energies(&mut automatic)
                    .expect("automatic constructor terms evaluate");
            let supplied_terms = crate::kernel::cf3d_frag_f24_contribution_energies(&mut supplied)
                .expect("supplied constructor terms evaluate");
            let expected_terms =
                crate::kernel::cf3d_frag_f24_contribution_energies(&mut shared_helper_reference)
                    .expect("shared-helper assembly-order reference terms evaluate");
            assert_eq!(automatic_terms.len(), expected_count);
            assert_eq!(supplied_terms.len(), expected_count);
            assert_eq!(expected_terms.len(), expected_count);
            let expected_term_bits = expected_terms
                .iter()
                .map(|energy| energy.to_bits())
                .collect::<Vec<_>>();
            assert_eq!(
                automatic_terms
                    .iter()
                    .map(|energy| energy.to_bits())
                    .collect::<Vec<_>>(),
                expected_term_bits,
                "automatic constructor must match the shared-helper assembly-order reference"
            );
            assert_eq!(
                supplied_terms
                    .iter()
                    .map(|energy| energy.to_bits())
                    .collect::<Vec<_>>(),
                expected_term_bits,
                "supplied constructor must match the shared-helper assembly-order reference"
            );

            let automatic_energy = cf3d_bld_b05_calc_energy(&mut automatic, &flat_coordinates)
                .expect("automatic full-kernel energy evaluates");
            let supplied_energy = cf3d_bld_b05_calc_energy(&mut supplied, &flat_coordinates)
                .expect("supplied full-kernel energy evaluates");
            let expected_energy =
                cf3d_bld_b05_calc_energy(&mut shared_helper_reference, &flat_coordinates)
                    .expect("shared-helper assembly-order full-kernel energy evaluates");
            assert_eq!(
                automatic_energy.to_bits(),
                expected_energy.to_bits(),
                "automatic total energy must match the shared-helper assembly-order reference"
            );
            assert_eq!(
                supplied_energy.to_bits(),
                expected_energy.to_bits(),
                "supplied total energy must match the shared-helper assembly-order reference"
            );

            let mut automatic_gradient = vec![0.0; flat_coordinates.len()];
            let mut supplied_gradient = vec![0.0; flat_coordinates.len()];
            let mut expected_gradient = vec![0.0; flat_coordinates.len()];
            cf3d_bld_b05_calc_grad(&mut automatic, &flat_coordinates, &mut automatic_gradient)
                .expect("automatic full-kernel gradient evaluates");
            cf3d_bld_b05_calc_grad(&mut supplied, &flat_coordinates, &mut supplied_gradient)
                .expect("supplied full-kernel gradient evaluates");
            cf3d_bld_b05_calc_grad(
                &mut shared_helper_reference,
                &flat_coordinates,
                &mut expected_gradient,
            )
            .expect("shared-helper assembly-order full-kernel gradient evaluates");
            let expected_gradient_bits = expected_gradient
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>();
            assert_eq!(
                automatic_gradient
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                expected_gradient_bits,
                "automatic gradient must match the shared-helper assembly-order reference"
            );
            assert_eq!(
                supplied_gradient
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                expected_gradient_bits,
                "supplied gradient must match the shared-helper assembly-order reference"
            );

            drop(automatic);
            drop(supplied);
            drop(shared_helper_reference);
            assert_eq!(automatic_coordinates, original_automatic_coordinates);
            assert_eq!(supplied_coordinates, original_supplied_coordinates);
        }
    }

    #[test]
    fn cf3d_frag_integration_requested_fragments_keep_complete_metadata() {
        let (topology, coordinates, properties) = cf3d_frag_integration_metadata_source();
        let original = (topology.clone(), coordinates.clone(), properties.clone());
        let coordinate_view = FragmentCoordinateView::from_coordinate_block(&coordinates);
        let fragments = get_molecule_fragments_with_coordinate_view(
            &topology,
            &coordinate_view,
            &properties,
            true,
            true,
        )
        .expect("the full-copy path retains complete source component metadata");
        assert_eq!(fragments.len(), 2);

        for component in 0..2 {
            let first_source_atom = component * 2;
            let fragment = &fragments[component];
            assert_eq!(
                fragment.component_atoms(),
                &[
                    AtomId::new(first_source_atom),
                    AtomId::new(first_source_atom + 1)
                ]
            );
            assert_eq!(fragment.topology().atoms.len(), 2);
            assert_eq!(fragment.topology().bonds.len(), 1);
            assert_eq!(fragment.topology().atoms[0].id(), AtomId::new(0));
            assert_eq!(fragment.topology().atoms[1].id(), AtomId::new(1));
            for (local_row, source_row) in [first_source_atom, first_source_atom + 1]
                .into_iter()
                .enumerate()
            {
                let atom = &fragment.topology().atoms[local_row];
                assert_eq!(
                    atom.prop("source-atom"),
                    Some(&cosmolkit_model::PropertyValue::String(format!(
                        "atom-{source_row}"
                    )))
                );
                assert!(atom.prop("computed-atom").is_none());
                assert!(atom.computed_prop_names().is_empty());
                assert_eq!(
                    atom.temporary_flags(),
                    0xF0F0_0000_0000_0000 | source_row as u64
                );
            }
            let bond = &fragment.topology().bonds[0];
            assert_eq!(bond.id(), BondId::new(0));
            assert_eq!(bond.begin(), AtomId::new(0));
            assert_eq!(bond.end(), AtomId::new(1));
            assert_eq!(
                bond.prop("source-bond"),
                Some(&cosmolkit_model::PropertyValue::String(format!(
                    "bond-{component}"
                )))
            );
            assert!(bond.prop("computed-bond").is_none());
            assert!(bond.computed_prop_names().is_empty());
            assert_eq!(
                bond.temporary_flags(),
                0x0F0F_0000_0000_0000 | component as u64
            );

            assert_eq!(
                fragment.topology().substance_groups,
                vec![cf3d_frag_integration_group(component, true)]
            );
            assert_eq!(
                fragment.topology().stereo_groups,
                vec![cf3d_frag_integration_stereo_group(component, true)]
            );
            assert!(fragment.topology().validate().is_ok());
            assert_eq!(
                fragment.topology_mapping().atoms().new_to_old(),
                &[
                    Some(AtomId::new(first_source_atom)),
                    Some(AtomId::new(first_source_atom + 1)),
                ]
            );
            assert_eq!(
                fragment.topology_mapping().bonds().new_to_old(),
                &[Some(BondId::new(component))]
            );

            let copied_coordinates = fragment.coordinates();
            assert_eq!(
                copied_coordinates.source_coordinate_dim,
                Some(CoordinateDimension::ThreeD)
            );
            assert_eq!(
                copied_coordinates
                    .conformers_2d
                    .iter()
                    .map(Conformer2D::id)
                    .collect::<Vec<_>>(),
                vec![31]
            );
            assert_eq!(
                copied_coordinates
                    .conformers_3d
                    .iter()
                    .map(Conformer3D::id)
                    .collect::<Vec<_>>(),
                vec![12, 31, 50]
            );
            assert_eq!(
                copied_coordinates
                    .conformers_3d
                    .iter()
                    .map(Conformer3D::is_3d)
                    .collect::<Vec<_>>(),
                vec![false, true, false]
            );
            for (source_conf, copied_conf) in coordinates
                .conformers_3d
                .iter()
                .zip(&copied_coordinates.conformers_3d)
            {
                assert_eq!(copied_conf.props(), source_conf.props());
                assert_eq!(
                    copied_conf.coordinates(),
                    &[
                        source_conf.coordinates()[first_source_atom],
                        source_conf.coordinates()[first_source_atom + 1],
                    ]
                );
            }
            assert_eq!(
                copied_coordinates.conformers_2d[0].props(),
                coordinates.conformers_2d[0].props()
            );
            assert_eq!(
                copied_coordinates.conformers_2d[0].coordinates(),
                &[
                    coordinates.conformers_2d[0].coordinates()[first_source_atom],
                    coordinates.conformers_2d[0].coordinates()[first_source_atom + 1],
                ]
            );

            let copied_properties = fragment.molecule_properties();
            assert_eq!(
                copied_properties.name(),
                Some("UFF-FRAG integration source")
            );
            assert_eq!(
                copied_properties.prop("ordinary-molecule"),
                Some("retained")
            );
            assert!(copied_properties.prop("computed-molecule").is_none());
            assert!(!copied_properties.is_prop_computed("computed-molecule"));
            assert_eq!(
                copied_properties.sdf_data_fields(),
                &[("SOURCE".to_owned(), "integration".to_owned())]
            );
            assert_eq!(copied_properties.sdf_property_lists().len(), 2);
            assert_eq!(
                copied_properties.sdf_property_lists()[0].values(),
                &[
                    Some(cosmolkit_model::PropertyValue::String(format!(
                        "atom-{first_source_atom}"
                    ))),
                    Some(cosmolkit_model::PropertyValue::String(format!(
                        "atom-{}",
                        first_source_atom + 1
                    ))),
                ]
            );
            assert_eq!(
                copied_properties.sdf_property_lists()[1].values(),
                &[Some(cosmolkit_model::PropertyValue::String(format!(
                    "bond-{component}"
                )))]
            );
        }

        let default_params = ParamCollection::get_params("").expect("pinned default UFF table");
        let total_valences = [4; 4];
        let conjugated = [false; 4];
        let valence = assignment(&[1; 4], &[0; 4]);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let mut typing_diagnostics = Vec::new();
        let (typed_params, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugated,
            default_params.as_ref(),
            &mut typing_diagnostics,
        )
        .expect("the metadata fixture's prepared state types");
        assert!(found_all);
        assert!(typing_diagnostics.is_empty());

        let mut automatic_coordinates = coordinates.clone();
        let mut automatic_diagnostics = Vec::new();
        let automatic = construct_force_field_with_automatic_typing(
            &topology,
            &mut automatic_coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut automatic_diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            true,
        )
        .expect("automatic full constructor completes the metadata fragment path");
        assert!(automatic_diagnostics.is_empty());
        drop(automatic);

        let mut supplied_coordinates = coordinates.clone();
        let mut supplied_diagnostics = Vec::new();
        let supplied = construct_force_field_with_params(
            &topology,
            &mut supplied_coordinates,
            31,
            &typed_params,
            &rings,
            &valence,
            &properties,
            &mut supplied_diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            true,
        )
        .expect("supplied full constructor completes the metadata fragment path");
        assert!(supplied_diagnostics.is_empty());
        drop(supplied);
        assert_eq!((topology, coordinates, properties), original);
        assert_eq!(automatic_coordinates, original.1);
        assert_eq!(supplied_coordinates, original.1);
    }

    #[test]
    fn cf3d_frag_integration_fragment_sanitize_errors_keep_inputs_unchanged() {
        let atoms = (0..6)
            .map(|row| {
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
                )
            })
            .collect();
        let bonds = (1..6)
            .map(|row| {
                Bond::from_spec(
                    BondId::new(row - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(row), BondOrder::Single),
                )
            })
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed overvalent carbon graph is structurally valid");
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let coordinates = cf3d_frag_f24_coordinates(&points);
        let original_coordinates = coordinates.clone();
        let properties = MoleculeProperties::default().with_name("fragment-error-atomicity");
        let original = (topology.clone(), coordinates.clone(), properties.clone());
        let total_valences = [4; 6];
        let conjugated = [false; 6];
        let valence = assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let default_params = ParamCollection::get_params("").expect("pinned default UFF table");
        let mut typing_diagnostics = Vec::new();
        let (typed_params, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugated,
            default_params.as_ref(),
            &mut typing_diagnostics,
        )
        .expect("the source-prepared C_3 rows are valid");
        assert!(found_all);
        assert!(typing_diagnostics.is_empty());

        let mut supplied_coordinates = coordinates.clone();
        let mut supplied_diagnostics = Vec::new();
        let supplied_error = match construct_force_field_with_params(
            &topology,
            &mut supplied_coordinates,
            31,
            &typed_params,
            &rings,
            &valence,
            &properties,
            &mut supplied_diagnostics,
            "[",
            1000.0,
            true,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("the source fragment sanitize error precedes torsion query parsing");
            }
        };
        assert!(matches!(
            supplied_error,
            ForceFieldConstructionError::Nonbonded(
                NonbondedAssemblyError::FragmentMapping(
                    NonbondedFragmentMappingError::FragmentCopy(source)
                )
            ) if source.component_index() == Some(0)
        ));
        assert!(supplied_diagnostics.is_empty());
        assert_eq!(supplied_coordinates, original_coordinates);

        let mut automatic_coordinates = coordinates.clone();
        let mut automatic_diagnostics = Vec::new();
        let automatic_error = match construct_force_field_with_automatic_typing(
            &topology,
            &mut automatic_coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut automatic_diagnostics,
            "[",
            1000.0,
            true,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("automatic typing must reach the same source fragment failure");
            }
        };
        assert!(matches!(
            automatic_error,
            AutomaticForceFieldConstructionError::Construction(
                ForceFieldConstructionError::Nonbonded(
                    NonbondedAssemblyError::FragmentMapping(
                        NonbondedFragmentMappingError::FragmentCopy(source)
                    )
                )
            ) if source.component_index() == Some(0)
        ));
        assert!(automatic_diagnostics.is_empty());
        assert_eq!(automatic_coordinates, original_coordinates);
        assert_eq!((topology, coordinates, properties), original);
    }

    #[test]
    fn uff_one_u06_misaligned_prepared_atom_rows_keep_typed_causes() {
        let topology = cf3d_bld_b15_topology(&[Hybridization::Sp2], &[]);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0], &[0]);
        let original_topology = topology.clone();

        let mut total_valence_coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0]]);
        let total_valence_error = match construct_force_field_with_automatic_typing(
            &topology,
            &mut total_valence_coordinates,
            31,
            &[],
            &[false],
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("a short prepared total-valence vector must be rejected");
            }
        };
        assert!(matches!(
            total_valence_error,
            AutomaticForceFieldConstructionError::Typing(UffTypingError::PreparedStateLength {
                input: UffTypingInput::TotalValence,
                expected: 1,
                actual: 0,
            })
        ));

        let mut conjugation_coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0]]);
        let conjugation_error = match construct_force_field_with_automatic_typing(
            &topology,
            &mut conjugation_coordinates,
            31,
            &[4],
            &[],
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        ) {
            Err(error) => error,
            Ok(field) => {
                drop(field);
                panic!("a short prepared conjugation vector must be rejected");
            }
        };
        assert!(matches!(
            conjugation_error,
            AutomaticForceFieldConstructionError::Typing(UffTypingError::PreparedStateLength {
                input: UffTypingInput::ConjugatedBondPresence,
                expected: 1,
                actual: 0,
            })
        ));

        assert_eq!(topology, original_topology);
    }

    #[test]
    fn uff_one_u06_aromatic_and_prepared_conjugation_select_source_rows() {
        let atoms = (0..5)
            .map(|row| {
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(Element::C)
                        .with_hybridization(Hybridization::Sp2)
                        .with_aromatic(row < 3),
                )
            })
            .collect();
        let mut bonds = [(0, 1), (1, 2), (2, 0)]
            .into_iter()
            .enumerate()
            .map(|(row, (begin, end))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic),
                )
            })
            .collect::<Vec<_>>();
        let mut conjugated_bond = Bond::from_spec(
            BondId::new(3),
            BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Double),
        );
        conjugated_bond.set_conjugated(true);
        bonds.push(conjugated_bond);
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed aromatic ring and conjugated bond topology is valid");
        let original_topology = topology.clone();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let total_valences = [4; 5];
        let no_conjugation = [false; 5];
        let conjugation = [false, false, false, true, true];
        let table = ParamCollection::get_params("").expect("pinned default UFF table");

        let mut base_diagnostics = Vec::new();
        let (base_params, base_found_all) = get_atom_types(
            &topology,
            &total_valences,
            &no_conjugation,
            table.as_ref(),
            &mut base_diagnostics,
        )
        .expect("fixed aromatic/conjugation rows align");
        assert!(base_found_all);
        assert!(base_diagnostics.is_empty());
        let c_r = table.get("C_R").expect("pinned source C_R row");
        let c_2 = table.get("C_2").expect("pinned source C_2 row");
        for atom_index in 0..3 {
            assert!(std::ptr::eq(
                base_params[atom_index].expect("aromatic C_R"),
                c_r
            ));
        }
        for atom_index in 3..5 {
            assert!(std::ptr::eq(
                base_params[atom_index].expect("nonconjugated C_2"),
                c_2
            ));
        }

        let mut typing_diagnostics = Vec::new();
        let (typed_params, found_all) = get_atom_types(
            &topology,
            &total_valences,
            &conjugation,
            table.as_ref(),
            &mut typing_diagnostics,
        )
        .expect("fixed prepared conjugation rows align");
        assert!(found_all);
        assert!(typing_diagnostics.is_empty());
        for (atom_index, params) in typed_params.iter().enumerate() {
            assert!(
                std::ptr::eq(params.expect("aromatic/conjugated C_R"), c_r),
                "source aromatic/conjugation predicate at atom {atom_index} selects C_R"
            );
        }

        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.866_025_403_784_438_6, 0.0],
            [5.0, 0.0, 0.0],
            [6.3, 0.0, 0.0],
        ];
        let valence = assignment(&[3, 3, 3, 2, 2], &[1, 1, 1, 2, 2]);
        let properties = MoleculeProperties::default();
        let mut automatic_coordinates = cf3d_frag_f24_coordinates(&points);
        let automatic_before = automatic_coordinates.clone();
        let mut automatic_diagnostics = Vec::new();
        let mut automatic = construct_force_field_with_automatic_typing(
            &topology,
            &mut automatic_coordinates,
            31,
            &total_valences,
            &conjugation,
            &rings,
            &valence,
            &properties,
            &mut automatic_diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        )
        .unwrap_or_else(|error| panic!("source automatic construction: {error:?}"));
        let mut supplied_coordinates = cf3d_frag_f24_coordinates(&points);
        let supplied_before = supplied_coordinates.clone();
        let mut supplied_diagnostics = Vec::new();
        let mut supplied = construct_force_field_with_params(
            &topology,
            &mut supplied_coordinates,
            31,
            &typed_params,
            &rings,
            &valence,
            &properties,
            &mut supplied_diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        )
        .unwrap_or_else(|error| panic!("source parameter constructor: {error:?}"));
        assert_eq!(automatic_diagnostics, [expected_needs_hydrogen_warning()]);
        assert_eq!(supplied_diagnostics, automatic_diagnostics);
        automatic.initialize().expect("automatic field initializes");
        supplied.initialize().expect("supplied field initializes");
        assert_eq!(
            cf3d_frag_accept_contribution_identities(&automatic),
            cf3d_frag_accept_contribution_identities(&supplied),
            "automatic overload forwards the source-typed nullable parameter rows"
        );
        let automatic_energies = cf3d_frag_f24_contribution_energies(&mut automatic)
            .expect("automatic contribution energies evaluate");
        let supplied_energies = cf3d_frag_f24_contribution_energies(&mut supplied)
            .expect("supplied contribution energies evaluate");
        assert_eq!(
            automatic_energies
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            supplied_energies
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        );
        drop(automatic);
        drop(supplied);
        assert_eq!(automatic_coordinates, automatic_before);
        assert_eq!(supplied_coordinates, supplied_before);
        assert_eq!(topology, original_topology);
        assert_eq!(valence, assignment(&[3, 3, 3, 2, 2], &[1, 1, 1, 2, 2]));
    }

    #[test]
    fn uff_one_u06_prepared_ring_and_valence_reach_source_constructor_stages() {
        let topology = cf3d_bld_b15_topology(
            &[Hybridization::Sp2; 3],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        let original_topology = topology.clone();
        let rings = cf3d_bld_b15_ring_info(&topology);
        assert_eq!(
            angle_order(
                Hybridization::Sp2,
                AtomId::new(1),
                AtomId::new(0),
                AtomId::new(2),
                &rings,
            ),
            35,
            "source Builder.cpp selects ring order 35 when all three atoms are in a 3-ring"
        );
        let total_valences = [4; 3];
        let conjugation = [false; 3];
        let valence = assignment(&[2; 3], &[2; 3]);
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.866_025_403_784_438_6, 0.0],
        ];
        let mut coordinates = cf3d_frag_f24_coordinates(&points);
        let original_coordinates = coordinates.clone();
        let mut diagnostics = Vec::new();
        let force_field = construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugation,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        )
        .unwrap_or_else(|error| panic!("source prepared-state construction: {error:?}"));
        assert_eq!(diagnostics, [expected_needs_hydrogen_warning()]);

        let ring_angles = cf3d_frag_accept_contribution_identities(&force_field)
            .into_iter()
            .filter_map(|identity| match identity {
                Cf3dFragAcceptContributionIdentity::AngleBend {
                    order,
                    theta0,
                    c0,
                    c1,
                    c2,
                    ..
                } => Some((order, theta0, c0, c1, c2)),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(ring_angles.len(), 3);
        let source_theta0 = 60.0 / 180.0 * std::f64::consts::PI;
        let source_sin_theta0 = source_theta0.sin();
        let source_cos_theta0 = source_theta0.cos();
        let source_c2 = 1.0 / (4.0 * (source_sin_theta0 * source_sin_theta0).max(1e-8));
        let source_c1 = -4.0 * source_c2 * source_cos_theta0;
        let source_c0 = source_c2 * (2.0 * source_cos_theta0 * source_cos_theta0 + 1.0);
        for (order, theta0, c0, c1, c2) in ring_angles {
            assert_eq!(order, 0, "source AngleBendContrib normalizes ring order 35");
            assert!((theta0 - source_theta0).abs() <= 1e-15);
            assert!((c0 - source_c0).abs() <= 1e-15);
            assert!((c1 - source_c1).abs() <= 1e-15);
            assert!((c2 - source_c2).abs() <= 1e-15);
        }
        drop(force_field);
        assert_eq!(coordinates, original_coordinates);
        assert_eq!(topology, original_topology);
    }

    #[test]
    fn uff_one_u07_selects_noncontiguous_three_d_id_and_borrows_its_rows() {
        let topology = cf3d_bld_b15_topology(&[Hybridization::Sp3], &[]);
        let params = [Some(cf3d_bld_b15_atomic_params())];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0], &[0]);
        let mut coordinates = CoordinateBlock::default();
        coordinates
            .conformers_2d
            .push(Conformer2D::new(17, vec![[70.0, 71.0]]));
        coordinates
            .conformers_3d
            .push(Conformer3D::new(2, vec![[2.0, 2.1, 2.2]], false));
        coordinates
            .conformers_3d
            .push(Conformer3D::new(17, vec![[17.0, 17.1, 17.2]], true));
        coordinates.conformers_3d.push(Conformer3D::new(
            1000,
            vec![[1000.0, 1000.1, 1000.2]],
            false,
        ));
        coordinates.source_coordinate_dim = Some(CoordinateDimension::ThreeD);
        let original_coordinates = coordinates.clone();
        let selected_row_ptr = coordinates.conformers_3d[1].coordinates()[0]
            .as_ptr()
            .cast::<f64>();

        let force_field = construct_force_field_with_params(
            &topology,
            &mut coordinates,
            17,
            &borrowed_params,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        )
        .unwrap_or_else(|error| panic!("source ID17 constructor: {error:?}"));
        assert_eq!(force_field.positions().len(), 1);
        assert_eq!(force_field.positions()[0].as_ptr(), selected_row_ptr);
        assert_eq!(force_field.positions()[0], [17.0, 17.1, 17.2]);
        drop(force_field);
        assert_eq!(coordinates, original_coordinates);
    }

    #[test]
    fn uff_one_u07_selects_the_single_explicit_three_d_row() {
        let topology = cf3d_bld_b15_topology(&[Hybridization::Sp3], &[]);
        let params = [Some(cf3d_bld_b15_atomic_params())];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0], &[0]);
        let mut coordinates = CoordinateBlock::default();
        coordinates
            .conformers_3d
            .push(Conformer3D::new(5, vec![[5.0, 6.0, 7.0]], false));
        coordinates.source_coordinate_dim = Some(CoordinateDimension::ThreeD);
        let selected_row_ptr = coordinates.conformers_3d[0].coordinates()[0]
            .as_ptr()
            .cast::<f64>();

        let force_field = construct_force_field_with_params(
            &topology,
            &mut coordinates,
            5,
            &borrowed_params,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut Vec::new(),
            DEFAULT_TORSION_BOND_SMARTS,
            1000.0,
            false,
        )
        .unwrap_or_else(|error| panic!("source ID5 constructor: {error:?}"));
        assert_eq!(force_field.positions().len(), 1);
        assert_eq!(force_field.positions()[0].as_ptr(), selected_row_ptr);
        assert_eq!(force_field.positions()[0], [5.0, 6.0, 7.0]);
    }

    #[test]
    fn uff_one_u07_missing_three_d_id_ignores_a_same_id_two_d_row() {
        let topology = cf3d_bld_b15_topology(&[Hybridization::Sp3], &[]);
        let params = [Some(cf3d_bld_b15_atomic_params())];
        let borrowed_params = cf3d_bld_borrowed_params(&params);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0], &[0]);

        for conformers_3d in [
            Vec::new(),
            vec![
                Conformer3D::new(12, vec![[12.0, 12.1, 12.2]], false),
                Conformer3D::new(50, vec![[50.0, 50.1, 50.2]], false),
            ],
        ] {
            let mut coordinates = CoordinateBlock::default();
            coordinates
                .conformers_2d
                .push(Conformer2D::new(31, vec![[31.0, 31.1]]));
            coordinates.conformers_3d = conformers_3d;
            let original_coordinates = coordinates.clone();
            let error = match construct_force_field_with_params(
                &topology,
                &mut coordinates,
                31,
                &borrowed_params,
                &rings,
                &valence,
                &MoleculeProperties::default(),
                &mut Vec::new(),
                DEFAULT_TORSION_BOND_SMARTS,
                1000.0,
                false,
            ) {
                Err(error) => error,
                Ok(field) => {
                    drop(field);
                    panic!("2D ID31 and unrelated 3D rows cannot satisfy the source 3D lookup");
                }
            };
            assert!(matches!(
                error,
                ForceFieldConstructionError::Builder(
                    UffBuilderError::SelectedThreeDimensionalConformerNotFound { conformer_id: 31 }
                )
            ));
            assert_eq!(coordinates, original_coordinates);
        }
    }

    fn uff_one_u08_atom(
        row: usize,
        atomic_number: u8,
        hybridization: Hybridization,
        aromatic: bool,
        formal_charge: i8,
        no_implicit: bool,
    ) -> Atom {
        let element = Element::from_atomic_number(atomic_number)
            .expect("fixed U08 atomic number exists in the element table");
        Atom::from_spec(
            AtomId::new(row),
            AtomSpec::new(element)
                .with_hybridization(hybridization)
                .with_aromatic(aromatic)
                .with_formal_charge(formal_charge)
                .with_no_implicit(no_implicit),
        )
    }

    fn uff_one_u08_topology(
        atoms: Vec<Atom>,
        edges: &[(usize, usize, BondOrder)],
    ) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end, order))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed U08 topology is structurally valid")
    }

    fn uff_one_u08_source_bond_parameters(
        ri: f64,
        rj: f64,
        xi: f64,
        xj: f64,
        zi: f64,
        zj: f64,
        bond_order: f64,
    ) -> (f64, f64) {
        // Independent pinned UFF equations from BondStretch.cpp:22-46,
        // with literal parameter rows supplied by Params.cpp at each call.
        let bond_order_correction = -0.1332 * (ri + rj) * bond_order.ln();
        let sqrt_xi = xi.sqrt();
        let sqrt_xj = xj.sqrt();
        let electronegativity_correction =
            ri * rj * (sqrt_xi - sqrt_xj).powi(2) / (xi * ri + xj * rj);
        let rest_length = ri + rj + bond_order_correction - electronegativity_correction;
        let force_constant = 2.0 * 332.06 * zi * zj / rest_length.powi(3);
        (rest_length, force_constant)
    }

    fn uff_one_u09_two_carbon_topology() -> TopologyBlock {
        uff_one_u08_topology(
            vec![
                uff_one_u08_atom(0, 6, Hybridization::Sp3, false, 0, true),
                uff_one_u08_atom(1, 6, Hybridization::Sp3, false, 0, true),
            ],
            &[(0, 1, BondOrder::Single)],
        )
    }

    fn uff_one_u09_optimize(
        topology: &TopologyBlock,
        coordinates: &mut CoordinateBlock,
        molecule_properties: &MoleculeProperties,
        diagnostics: &mut Vec<UffTypingDiagnostic>,
        options: SingleConformerOptions,
    ) -> Result<super::super::convenience::OptimizationOutcome, SingleConformerOptimizationError>
    {
        let total_valences = [4, 4];
        let conjugated_bonds = [false; 2];
        let rings = cf3d_bld_b15_ring_info(topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        optimize_single_conformer(
            topology,
            coordinates,
            &total_valences,
            &conjugated_bonds,
            &rings,
            &valence,
            molecule_properties,
            diagnostics,
            options,
        )
    }

    #[test]
    fn uff_one_u08_missing_type_remains_nullable_during_constructor_delegation() {
        // AtomTyper.cpp:507-533 leaves an unrecognized row null and returns
        // foundAll=false; Builder.cpp:712-719 delegates that vector anyway.
        let missing_element =
            Element::from_atomic_number(0).expect("the source dummy element is modeled");
        let missing_atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(missing_element)
                .with_hybridization(Hybridization::Unspecified)
                .with_no_implicit(true)
                .with_prop("dummyLabel", "u08-unrecognized")
                .expect("fixed U08 dummy label is valid"),
        );
        let topology = uff_one_u08_topology(
            vec![
                missing_atom,
                uff_one_u08_atom(1, 6, Hybridization::Sp3, false, 0, true),
                uff_one_u08_atom(2, 6, Hybridization::Sp3, false, 0, true),
            ],
            &[(1, 2, BondOrder::Single)],
        );
        let total_valences = [0, 4, 4];
        let conjugated = [false; 3];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0, 1, 1], &[0; 3]);
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.5, 0.0, 0.0]];
        let mut coordinates = cf3d_frag_f24_coordinates(&points);
        let mut diagnostics = Vec::new();
        let force_field = construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            0.0,
            false,
        )
        .expect("source automatic constructor delegates a nullable parameter row");
        assert_eq!(
            diagnostics,
            [UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: "UFFTYPER: Unrecognized atom type: ",
            }]
        );
        let identities = cf3d_frag_accept_contribution_identities(&force_field);
        assert_eq!(identities.len(), 1);
        let Cf3dFragAcceptContributionIdentity::BondStretch {
            end1_idx,
            end2_idx,
            rest_len,
            force_constant,
        } = identities[0]
        else {
            panic!("source addBonds retains the supported C_3-C_3 term");
        };
        assert_eq!((end1_idx, end2_idx), (1, 2));
        let source_c3_bond =
            uff_one_u08_source_bond_parameters(0.757, 0.757, 5.343, 5.343, 1.912, 1.912, 1.0);
        assert!((rest_len - source_c3_bond.0).abs() <= 1.0e-12);
        assert!((force_constant - source_c3_bond.1).abs() <= 1.0e-12);
    }

    #[test]
    fn uff_one_u08_aromatic_constructor_uses_the_pinned_c_r_row() {
        // RDKit testUFFHelpers.cpp:35-51 expects aromatic carbons to type C_R.
        let atoms = (0..6)
            .map(|row| uff_one_u08_atom(row, 6, Hybridization::Sp2, true, 0, true))
            .collect();
        let topology = uff_one_u08_topology(
            atoms,
            &[
                (0, 1, BondOrder::Aromatic),
                (1, 2, BondOrder::Aromatic),
                (2, 3, BondOrder::Aromatic),
                (3, 4, BondOrder::Aromatic),
                (4, 5, BondOrder::Aromatic),
                (5, 0, BondOrder::Aromatic),
            ],
        );
        let total_valences = [4; 6];
        let conjugated = [false; 6];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[3; 6], &[0; 6]);
        let sqrt3_over2 = 3.0_f64.sqrt() / 2.0;
        let points = [
            [1.0, 0.0, 0.0],
            [0.5, sqrt3_over2, 0.0],
            [-0.5, sqrt3_over2, 0.0],
            [-1.0, 0.0, 0.0],
            [-0.5, -sqrt3_over2, 0.0],
            [0.5, -sqrt3_over2, 0.0],
        ];
        let mut coordinates = cf3d_frag_f24_coordinates(&points);
        let mut diagnostics = Vec::new();
        let force_field = construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            0.0,
            false,
        )
        .expect("pinned aromatic C_R rows construct the source field");
        assert!(diagnostics.is_empty());

        let identities = cf3d_frag_accept_contribution_identities(&force_field);
        assert_eq!(identities.len(), 18);
        let source_bond_parameters =
            uff_one_u08_source_bond_parameters(0.729, 0.729, 5.343, 5.343, 1.912, 1.912, 1.5);
        for (identity, (expected_begin, expected_end)) in
            identities[..6]
                .iter()
                .zip([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)])
        {
            let Cf3dFragAcceptContributionIdentity::BondStretch {
                end1_idx,
                end2_idx,
                rest_len,
                force_constant,
            } = identity
            else {
                panic!("source addBonds emits the six aromatic bonds first");
            };
            assert_eq!((*end1_idx, *end2_idx), (expected_begin, expected_end));
            assert!((*rest_len - source_bond_parameters.0).abs() <= 1.0e-12);
            assert!((*force_constant - source_bond_parameters.1).abs() <= 1.0e-12);
        }
        assert!(identities[6..12].iter().all(|identity| matches!(
            identity,
            Cf3dFragAcceptContributionIdentity::AngleBend { .. }
        )));
        assert!(identities[12..].iter().all(|identity| matches!(
            identity,
            Cf3dFragAcceptContributionIdentity::TorsionAngle { .. }
        )));
    }

    #[test]
    fn uff_one_u08_charged_mg_constructor_uses_the_source_valence_suffix() {
        // RDKit testUFFHelpers.cpp:83-101 expects Mg3+2 for [Mg]=C.
        let topology = uff_one_u08_topology(
            vec![
                uff_one_u08_atom(0, 12, Hybridization::Sp3, false, 0, true),
                uff_one_u08_atom(1, 6, Hybridization::Sp2, false, 0, true),
            ],
            &[(0, 1, BondOrder::Double)],
        );
        let total_valences = [2, 2];
        let conjugated = [false; 2];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[2, 2], &[0, 0]);
        let points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut coordinates = cf3d_frag_f24_coordinates(&points);
        let mut diagnostics = Vec::new();
        let force_field = construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            0.0,
            false,
        )
        .expect("pinned Mg3+2 and C_2 rows form the source bond");
        assert!(diagnostics.is_empty());
        let identities = cf3d_frag_accept_contribution_identities(&force_field);
        assert_eq!(identities.len(), 1);
        let Cf3dFragAcceptContributionIdentity::BondStretch {
            end1_idx,
            end2_idx,
            rest_len,
            force_constant,
        } = identities[0]
        else {
            panic!("source addBonds emits the recognized Mg=C bond");
        };
        assert_eq!((end1_idx, end2_idx), (0, 1));
        let source_bond_parameters =
            uff_one_u08_source_bond_parameters(1.421, 0.732, 3.951, 5.343, 1.787, 1.912, 2.0);
        assert!((rest_len - source_bond_parameters.0).abs() <= 1.0e-12);
        assert!((force_constant - source_bond_parameters.1).abs() <= 1.0e-12);
    }

    #[test]
    fn uff_one_u08_salt_constructor_keeps_source_fragment_interactions() {
        // AtomTyper.cpp leaves alkali-metal and halogen labels unhybridized;
        // the pinned Na and Cl rows are in ForceField/UFF/Params.cpp.
        let topology = uff_one_u08_topology(
            vec![
                uff_one_u08_atom(0, 11, Hybridization::Unspecified, false, 1, true),
                uff_one_u08_atom(1, 17, Hybridization::Unspecified, false, -1, true),
            ],
            &[],
        );
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0, 0], &[0, 0]);
        let points = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let expected_minimum = (2.983_f64 * 3.947).sqrt();
        let expected_depth = (0.03_f64 * 0.227).sqrt();

        for ignore_interfragment_interactions in [true, false] {
            let mut coordinates = cf3d_frag_f24_coordinates(&points);
            let mut diagnostics = Vec::new();
            let force_field = construct_force_field_with_automatic_typing(
                &topology,
                &mut coordinates,
                31,
                &total_valences,
                &conjugated,
                &rings,
                &valence,
                &MoleculeProperties::default(),
                &mut diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                10.0,
                ignore_interfragment_interactions,
            )
            .unwrap_or_else(|error| panic!("source salt construction failed: {error:?}"));
            assert!(diagnostics.is_empty());
            let identities = cf3d_frag_accept_contribution_identities(&force_field);
            if ignore_interfragment_interactions {
                assert!(identities.is_empty());
            } else {
                assert_eq!(identities.len(), 1);
                let Cf3dFragAcceptContributionIdentity::Vdw {
                    at1_idx,
                    at2_idx,
                    x_ij,
                    well_depth,
                    threshold,
                } = identities[0]
                else {
                    panic!("source addNonbonded emits the close interfragment Na-Cl pair");
                };
                assert_eq!((at1_idx, at2_idx), (0, 1));
                assert_eq!(x_ij.to_bits(), expected_minimum.to_bits());
                assert_eq!(well_depth.to_bits(), expected_depth.to_bits());
                assert_eq!(threshold.to_bits(), (10.0 * expected_minimum).to_bits());
            }
        }
    }

    #[test]
    fn uff_one_u08_tbp_constructor_uses_source_p_type_and_angle_order() {
        // testUFFHelpers.cpp:95-101 pins P_3+5; Builder.cpp:412-430 and
        // 248-405 pin the special TBP dispatch and ten ordered angles.
        let mut atoms = vec![uff_one_u08_atom(0, 15, Hybridization::Sp3d, false, 0, true)];
        atoms
            .extend((1..6).map(|row| uff_one_u08_atom(row, 6, Hybridization::Sp3, false, 0, true)));
        let edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
            .map(|(begin, end)| (begin, end, BondOrder::Single));
        let topology = uff_one_u08_topology(atoms, &edges);
        let total_valences = [5, 1, 1, 1, 1, 1];
        let conjugated = [false; 6];
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]);
        let points = CF3D_BLD_B20_COORDINATES
            .chunks_exact(3)
            .map(|point| [point[0], point[1], point[2]])
            .collect::<Vec<_>>();
        let mut coordinates = cf3d_frag_f24_coordinates(&points);
        let mut diagnostics = Vec::new();
        let force_field = construct_force_field_with_automatic_typing(
            &topology,
            &mut coordinates,
            31,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            DEFAULT_TORSION_BOND_SMARTS,
            0.0,
            false,
        )
        .expect("pinned P_3+5 and C_3 rows construct the TBP field");
        assert_eq!(
            diagnostics,
            [UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Warning,
                message_prefix: "UFFTYPER: Warning: hybridization set to SP3 for atom ",
            }]
        );
        let identities = cf3d_frag_accept_contribution_identities(&force_field);
        assert_eq!(identities.len(), 15);
        let source_bond_parameters =
            uff_one_u08_source_bond_parameters(1.056, 0.757, 5.463, 5.343, 2.863, 1.912, 1.0);
        for (row, identity) in identities[..5].iter().enumerate() {
            let Cf3dFragAcceptContributionIdentity::BondStretch {
                end1_idx,
                end2_idx,
                rest_len,
                force_constant,
            } = identity
            else {
                panic!("source addBonds emits five P-C bonds before TBP angles");
            };
            assert_eq!((*end1_idx, *end2_idx), (0, row as u32 + 1));
            assert!((*rest_len - source_bond_parameters.0).abs() <= 1.0e-12);
            assert!((*force_constant - source_bond_parameters.1).abs() <= 1.0e-12);
        }
        for (identity, (expected_at1, expected_at2, expected_at3, expected_order)) in
            identities[5..].iter().zip(CF3D_BLD_B21_SOURCE_TERMS)
        {
            let Cf3dFragAcceptContributionIdentity::AngleBend {
                at1_idx,
                at2_idx,
                at3_idx,
                order,
                ..
            } = identity
            else {
                panic!("source TBP stage emits its ten angle terms after bonds");
            };
            assert_eq!(
                (
                    *at1_idx as usize,
                    *at2_idx as usize,
                    *at3_idx as usize,
                    *order
                ),
                (expected_at1, expected_at2, expected_at3, expected_order)
            );
        }
    }

    #[test]
    fn uff_one_u09_converged_path_mutates_only_the_selected_three_d_rows() {
        // UFF.h:40-47 delegates selected conformer construction to the
        // ordinary OptimizeMolecule path; testUFFHelpers.cpp:502-519 expects
        // status zero and mutation of that selected conformer's positions.
        let topology = uff_one_u09_two_carbon_topology();
        let molecule_properties = MoleculeProperties::default()
            .with_name("fixed-u09")
            .with_prop("source", "preserve")
            .expect("fixed property key is valid");
        let original_molecule_properties = molecule_properties.clone();
        let mut coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let original_two_d = coordinates.conformers_2d[0].clone();
        let original_before = coordinates.conformers_3d[0].clone();
        let original_selected = coordinates.conformers_3d[1].clone();
        let original_after = coordinates.conformers_3d[2].clone();
        let original_dimension = coordinates.source_coordinate_dim;
        let mut diagnostics = Vec::new();
        let outcome = uff_one_u09_optimize(
            &topology,
            &mut coordinates,
            &molecule_properties,
            &mut diagnostics,
            SingleConformerOptions::for_conformer(31),
        )
        .expect("the source C_3-C_3 convenience path converges");

        assert_eq!(outcome.status, 0);
        assert!(outcome.energy.is_finite());
        assert_ne!(
            coordinates.conformers_3d[1].coordinates(),
            original_selected.coordinates()
        );
        assert_eq!(coordinates.conformers_2d[0], original_two_d);
        assert_eq!(coordinates.conformers_3d[0], original_before);
        assert_eq!(coordinates.conformers_3d[2], original_after);
        assert_eq!(coordinates.conformers_3d[1].id(), original_selected.id());
        assert_eq!(
            coordinates.conformers_3d[1].is_3d(),
            original_selected.is_3d()
        );
        assert_eq!(
            coordinates.conformers_3d[1].props(),
            original_selected.props()
        );
        assert_eq!(coordinates.source_coordinate_dim, original_dimension);
        assert_eq!(molecule_properties, original_molecule_properties);
    }

    #[test]
    fn uff_serial_s08_reference_id_may_be_in_or_outside_processed_rows() {
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let mut reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let properties = MoleculeProperties::default();
        let mut diagnostics = Vec::new();

        let mut reference_id_rows = [[0.0, 0.0, 0.0], [1.514, 0.0, 0.0]];
        let mut includes_reference = [SerialConformer {
            id: 31,
            positions: &mut reference_id_rows,
        }];
        let mut reference_result = [OptimizationOutcome {
            status: -9,
            energy: 123.0,
        }];
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = 0;
        options.ignore_interfragment_interactions = false;
        optimize_serial_uff(
            &topology,
            &mut reference,
            &mut includes_reference,
            &mut reference_result,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect("reference ID may also occur in the explicitly processed set");
        assert_eq!(reference_result[0].energy, 0.0);

        let mut other_id_rows = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut excludes_reference = [SerialConformer {
            id: 777,
            positions: &mut other_id_rows,
        }];
        let mut other_result = [OptimizationOutcome {
            status: -9,
            energy: 123.0,
        }];
        optimize_serial_uff(
            &topology,
            &mut reference,
            &mut excludes_reference,
            &mut other_result,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect("processing an ID outside the reference block remains explicit");
        let rest_length: f64 = 0.757 + 0.757;
        let force_constant = 2.0 * 332.06 * 1.912 * 1.912 / rest_length.powi(3);
        let displacement = 2.0 - rest_length;
        let expected_energy = 0.5 * force_constant * displacement * displacement;
        assert!((other_result[0].energy - expected_energy).abs() <= 1.0e-12);
        assert_eq!(reference.conformers_3d[1].id(), 31);
        assert_eq!(other_id_rows, [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    }

    #[test]
    fn uff_serial_s08_missing_reference_fails_before_empty_traversal() {
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let mut reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let original_reference = reference.clone();
        let properties = MoleculeProperties::default();
        let mut diagnostics = Vec::new();
        let mut results = [];
        let mut options = SingleConformerOptions::for_conformer(99);
        options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let error = optimize_serial_uff(
            &topology,
            &mut reference,
            &mut [],
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect_err("construction must resolve the explicit reference before empty traversal");
        assert!(matches!(
            error,
            SerialUffOptimizationError::Construction(
                AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 99
                        }
                    )
                )
            )
        ));
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts().0, 1);
        assert!(results.is_empty());
        assert_eq!(reference, original_reference);
    }

    #[test]
    fn uff_serial_s08_empty_processed_set_constructs_once_then_returns() {
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let mut reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let properties = MoleculeProperties::default();
        let mut diagnostics = Vec::new();
        let mut results = [];
        let mut options = SingleConformerOptions::for_conformer(31);
        options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        optimize_serial_uff(
            &topology,
            &mut reference,
            &mut [],
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect("source constructs its field even when no conformers are processed");
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
        assert!(results.is_empty());
        assert_eq!(reference.conformers_3d[1].id(), 31);
    }

    #[test]
    fn uff_serial_s09_reference_geometry_selects_one_vdw_set_for_every_row() {
        // Builder.cpp:433-469 selects vdW terms from the construction conformer
        // before serial traversal. Nonbonded.cpp:21-27 and 54-69 pin the
        // geometric-mean parameters and energy; Params.cpp pins Na and Cl.
        let topology = uff_one_u08_topology(
            vec![
                uff_one_u08_atom(0, 11, Hybridization::Unspecified, false, 1, true),
                uff_one_u08_atom(1, 17, Hybridization::Unspecified, false, -1, true),
            ],
            &[],
        );
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0, 0], &[0, 0]);
        let total_valences = [0, 0];
        let conjugated = [false; 2];
        let molecule_properties = MoleculeProperties::default();
        let expected_minimum = (2.983_f64 * 3.947).sqrt();
        let expected_depth = (0.03_f64 * 0.227).sqrt();
        let source_vdw_energy = |distance: f64| {
            let ratio = expected_minimum / distance;
            let ratio3 = ratio * ratio * ratio;
            let ratio6 = ratio3 * ratio3;
            expected_depth * (ratio6 * ratio6 - 2.0 * ratio6)
        };
        let expected_energies = [source_vdw_energy(2.0), source_vdw_energy(20.0)];
        assert!(20.0 < 10.0 * expected_minimum);

        for (reference_distance, reference_includes_pair) in [(2.0, true), (40.0, false)] {
            for ignore_interfragment_interactions in [true, false] {
                let pair_is_included =
                    reference_includes_pair && !ignore_interfragment_interactions;
                let expected = if pair_is_included {
                    expected_energies
                } else {
                    [0.0, 0.0]
                };
                let expected_terms = usize::from(pair_is_included);
                let mut reference =
                    cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [reference_distance, 0.0, 0.0]]);
                let original_reference = reference.clone();
                let mut near_row = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
                let mut far_row = [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]];
                let mut conformers = [
                    SerialConformer {
                        id: 80,
                        positions: &mut near_row,
                    },
                    SerialConformer {
                        id: 2,
                        positions: &mut far_row,
                    },
                ];
                let mut results = [OptimizationOutcome {
                    status: -9,
                    energy: 123.0,
                }; 2];
                let mut diagnostics = Vec::new();
                let mut options = SingleConformerOptions::for_conformer(31);
                options.max_iterations = 0;
                options.ignore_interfragment_interactions = ignore_interfragment_interactions;

                crate::kernel::cf3d_uff_one_kernel_counts_reset();
                optimize_serial_uff(
                    &topology,
                    &mut reference,
                    &mut conformers,
                    &mut results,
                    &total_valences,
                    &conjugated,
                    &rings,
                    &valence,
                    &molecule_properties,
                    &mut diagnostics,
                    options,
                )
                .expect("fixed Na-Cl serial inputs construct and evaluate");

                assert!(diagnostics.is_empty());
                assert!((results[0].energy - expected[0]).abs() <= 1.0e-12);
                assert!((results[1].energy - expected[1]).abs() <= 1.0e-12);
                assert_eq!(
                    crate::kernel::cf3d_uff_one_kernel_counts(),
                    (1, expected_terms, 0)
                );
                assert_eq!(near_row, [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
                assert_eq!(far_row, [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]);
                assert_eq!(reference, original_reference);
            }
        }
    }

    #[test]
    fn uff_serial_s10_reference_tbp_angles_survive_rebinding_and_input_order() {
        // Builder.cpp:248-405 selects axial/equatorial bonds from confId once;
        // Builder.cpp:412-426 dispatches that special-case at construction.
        let mut atoms = vec![uff_one_u08_atom(0, 15, Hybridization::Sp3d, false, 0, true)];
        atoms
            .extend((1..6).map(|row| uff_one_u08_atom(row, 6, Hybridization::Sp3, false, 0, true)));
        let edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
            .map(|(begin, end)| (begin, end, BondOrder::Single));
        let topology = uff_one_u08_topology(atoms, &edges);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]);
        let total_valences = [5, 1, 1, 1, 1, 1];
        let conjugated = [false; 6];
        let properties = MoleculeProperties::default();
        let original_points: [[f64; 3]; 6] = CF3D_BLD_B20_COORDINATES
            .chunks_exact(3)
            .map(|point| [point[0], point[1], point[2]])
            .collect::<Vec<_>>()
            .try_into()
            .expect("the fixed B20 reference has six points");
        let diagonal = 1.0 / 3.0_f64.sqrt();
        let alternate_points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
            [diagonal, diagonal, diagonal],
        ];
        // Source bond-ID order selects axial atoms 2 and 4 for this geometry:
        // their unit-vector dot is -1, and all remaining neighbors stay in source order.
        const ALTERNATE_SOURCE_TERMS: [(usize, usize, usize, u32); 10] = [
            (2, 0, 4, 2),
            (1, 0, 3, 3),
            (1, 0, 5, 3),
            (3, 0, 5, 3),
            (2, 0, 1, 0),
            (2, 0, 3, 0),
            (2, 0, 5, 0),
            (4, 0, 1, 0),
            (4, 0, 3, 0),
            (4, 0, 5, 0),
        ];

        for (reference_label, reference_points, source_terms) in [
            ("B20", original_points, CF3D_BLD_B21_SOURCE_TERMS),
            ("alternate", alternate_points, ALTERNATE_SOURCE_TERMS),
        ] {
            let mut expected_reference = cf3d_frag_f24_coordinates(&reference_points);
            let mut expected_diagnostics = Vec::new();
            let mut expected_field = construct_force_field_with_automatic_typing(
                &topology,
                &mut expected_reference,
                31,
                &total_valences,
                &conjugated,
                &rings,
                &valence,
                &properties,
                &mut expected_diagnostics,
                DEFAULT_TORSION_BOND_SMARTS,
                0.0,
                false,
            )
            .expect("pinned P_3+5 reference field constructs");
            assert_eq!(
                expected_diagnostics,
                [UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(0)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: "UFFTYPER: Warning: hybridization set to SP3 for atom ",
                }],
                "{reference_label}: source typing diagnostic is preserved"
            );
            let expected_identities = cf3d_frag_accept_contribution_identities(&expected_field);
            assert_eq!(expected_identities.len(), 15, "{reference_label}");
            for (identity, (expected_at1, expected_at2, expected_at3, expected_order)) in
                expected_identities[5..].iter().zip(source_terms)
            {
                let Cf3dFragAcceptContributionIdentity::AngleBend {
                    at1_idx,
                    at2_idx,
                    at3_idx,
                    order,
                    ..
                } = identity
                else {
                    panic!("the pinned constructor appends five bonds before ten TBP angles");
                };
                assert_eq!(
                    (
                        *at1_idx as usize,
                        *at2_idx as usize,
                        *at3_idx as usize,
                        *order,
                    ),
                    (expected_at1, expected_at2, expected_at3, expected_order),
                    "{reference_label}: source term identity and append order"
                );
            }
            expected_field
                .initialize()
                .expect("source-listed expected TBP field initializes");

            let expected_original = cf3d_bld_b05_calc_energy(
                &mut expected_field,
                &original_points
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<Vec<_>>(),
            )
            .expect("source-listed angles evaluate at the fixed B20 geometry");
            let expected_alternate = cf3d_bld_b05_calc_energy(
                &mut expected_field,
                &alternate_points
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<Vec<_>>(),
            )
            .expect("source-listed angles evaluate at the alternate TBP geometry");
            drop(expected_field);

            for original_first in [false, true] {
                let (
                    first_points,
                    second_points,
                    first_expected,
                    second_expected,
                    first_id,
                    second_id,
                ) = if original_first {
                    (
                        original_points,
                        alternate_points,
                        expected_original,
                        expected_alternate,
                        2,
                        80,
                    )
                } else {
                    (
                        alternate_points,
                        original_points,
                        expected_alternate,
                        expected_original,
                        80,
                        2,
                    )
                };
                let mut reference = cf3d_frag_f24_coordinates(&reference_points);
                let mut first_row = first_points;
                let mut second_row = second_points;
                let original_first_row = first_row;
                let original_second_row = second_row;
                let mut conformers = [
                    SerialConformer {
                        id: first_id,
                        positions: &mut first_row,
                    },
                    SerialConformer {
                        id: second_id,
                        positions: &mut second_row,
                    },
                ];
                let mut results = [OptimizationOutcome {
                    status: -9,
                    energy: 123.0,
                }; 2];
                let mut diagnostics = Vec::new();
                let mut options = SingleConformerOptions::for_conformer(31);
                options.max_iterations = 0;
                options.vdw_threshold = 0.0;
                options.ignore_interfragment_interactions = false;

                crate::kernel::cf3d_uff_one_kernel_counts_reset();
                optimize_serial_uff(
                    &topology,
                    &mut reference,
                    &mut conformers,
                    &mut results,
                    &total_valences,
                    &conjugated,
                    &rings,
                    &valence,
                    &properties,
                    &mut diagnostics,
                    options,
                )
                .expect("serial TBP optimization retains the fixed reference term set");

                assert_eq!(diagnostics, expected_diagnostics, "{reference_label}");
                assert_eq!(results[0].status, 1, "{reference_label}, first slot");
                assert_eq!(results[1].status, 1, "{reference_label}, second slot");
                assert_eq!(results[0].energy.to_bits(), first_expected.to_bits());
                assert_eq!(results[1].energy.to_bits(), second_expected.to_bits());
                assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 15, 0));
                assert_eq!(first_row, original_first_row);
                assert_eq!(second_row, original_second_row);
            }
        }
    }

    #[test]
    fn uff_serial_s12_complete_entry_preserves_reference_metadata_on_success() {
        // UFF.h:69-78 constructs from the explicit reference before the
        // FFConvenience.h:66-77 serial traversal. Detached input rows are the
        // only mutable coordinates in this entry; other coordinate sets and
        // their conformer metadata remain caller-owned.
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let properties = MoleculeProperties::default()
            .with_name("serial-s12")
            .with_prop("owner", "fixed-regression")
            .expect("fixed S12 molecule property key is valid");
        let mut reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let reference_before = reference.clone();
        let mut first = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut second = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let first_before = first;
        let second_before = second;
        let mut conformers = [
            SerialConformer {
                id: 80,
                positions: &mut first,
            },
            SerialConformer {
                id: 2,
                positions: &mut second,
            },
        ];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut results = [sentinel; 2];
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = 0;
        options.ignore_interfragment_interactions = false;

        optimize_serial_uff(
            &topology,
            &mut reference,
            &mut conformers,
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect("the complete source-ordered serial entry returns normal status 1");

        assert!(diagnostics.is_empty());
        assert!(results.iter().all(|result| result.status == 1));
        assert!(results.iter().all(|result| result.energy.is_finite()));
        assert_eq!(
            conformers
                .iter()
                .map(|conformer| conformer.id)
                .collect::<Vec<_>>(),
            [80, 2]
        );
        drop(conformers);
        assert_eq!(first, first_before);
        assert_eq!(second, second_before);
        assert_eq!(reference, reference_before);
    }

    #[test]
    fn uff_serial_s12_complete_entry_preserves_metadata_on_typed_failures() {
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let properties = MoleculeProperties::default();
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };

        let mut missing_reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let missing_reference_before = missing_reference.clone();
        let mut missing_reference_row = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let missing_reference_row_before = missing_reference_row;
        let mut missing_reference_conformers = [SerialConformer {
            id: 80,
            positions: &mut missing_reference_row,
        }];
        let mut missing_reference_results = [sentinel];
        let mut diagnostics = Vec::new();
        let mut missing_reference_options = SingleConformerOptions::for_conformer(99);
        missing_reference_options.max_iterations = 0;
        missing_reference_options.ignore_interfragment_interactions = false;
        let construction_error = optimize_serial_uff(
            &topology,
            &mut missing_reference,
            &mut missing_reference_conformers,
            &mut missing_reference_results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            missing_reference_options,
        )
        .expect_err("the explicit source construction reference must resolve before traversal");
        assert!(matches!(
            construction_error,
            SerialUffOptimizationError::Construction(
                AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 99
                        }
                    )
                )
            )
        ));
        assert_eq!(missing_reference, missing_reference_before);
        assert_eq!(missing_reference_row, missing_reference_row_before);
        assert_eq!(missing_reference_results, [sentinel]);

        let equilibrium_distance = 0.757_f64 + 0.757;
        let mut equilibrium_reference =
            cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [equilibrium_distance, 0.0, 0.0]]);
        let equilibrium_reference_before = equilibrium_reference.clone();
        let mut equilibrium_row = [[0.0, 0.0, 0.0], [equilibrium_distance, 0.0, 0.0]];
        let equilibrium_row_before = equilibrium_row;
        let mut equilibrium_conformers = [SerialConformer {
            id: 41,
            positions: &mut equilibrium_row,
        }];
        let mut equilibrium_results = [sentinel];
        let mut equilibrium_options = SingleConformerOptions::for_conformer(31);
        equilibrium_options.ignore_interfragment_interactions = false;
        let optimization_error = optimize_serial_uff(
            &topology,
            &mut equilibrium_reference,
            &mut equilibrium_conformers,
            &mut equilibrium_results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            equilibrium_options,
        )
        .expect_err("the source zero-gradient minimization error remains typed");
        assert!(matches!(
            optimization_error,
            SerialUffOptimizationError::Optimization(
                super::super::convenience::SerialConformerOptimizationError::Optimization {
                    input_index: 0,
                    conformer_id: 41,
                    source: OptimizationStageError::Minimize(
                        ForceFieldKernelError::OptimizerBadDirection
                    ),
                }
            )
        ));
        assert_eq!(equilibrium_reference, equilibrium_reference_before);
        assert_eq!(equilibrium_row, equilibrium_row_before);
        assert_eq!(equilibrium_results, [sentinel]);
    }

    #[test]
    fn uff_serial_s12_missing_parameters_keep_source_diagnostic_and_partial_terms() {
        // AtomTyper.cpp:507-533 leaves an unknown row null and logs it;
        // Builder.cpp:712-719 delegates that parameter vector regardless of
        // foundAll. The supported C3-C3 bond remains the one stored term.
        let missing_element =
            Element::from_atomic_number(0).expect("the source dummy element is modeled");
        let missing_atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(missing_element)
                .with_hybridization(Hybridization::Unspecified)
                .with_no_implicit(true)
                .with_prop("dummyLabel", "s12-unrecognized")
                .expect("fixed S12 dummy label is valid"),
        );
        let topology = uff_one_u08_topology(
            vec![
                missing_atom,
                uff_one_u08_atom(1, 6, Hybridization::Sp3, false, 0, true),
                uff_one_u08_atom(2, 6, Hybridization::Sp3, false, 0, true),
            ],
            &[(1, 2, BondOrder::Single)],
        );
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[0, 1, 1], &[0; 3]);
        let total_valences = [0, 4, 4];
        let conjugated = [false; 3];
        let properties = MoleculeProperties::default();
        let points = [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut reference = cf3d_frag_f24_coordinates(&points);
        let reference_before = reference.clone();
        let mut first = points;
        let mut second = [[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.5, 0.0, 0.0]];
        let first_before = first;
        let second_before = second;
        let mut conformers = [
            SerialConformer {
                id: 80,
                positions: &mut first,
            },
            SerialConformer {
                id: 2,
                positions: &mut second,
            },
        ];
        let mut results = [OptimizationOutcome {
            status: -9,
            energy: 123.0,
        }; 2];
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = 0;
        options.ignore_interfragment_interactions = true;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();

        optimize_serial_uff(
            &topology,
            &mut reference,
            &mut conformers,
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect("the nullable source parameter row is nonfatal");

        drop(conformers);
        assert_eq!(
            diagnostics,
            [UffTypingDiagnostic {
                atom_id: Some(AtomId::new(0)),
                kind: UffTypingDiagnosticKind::Error,
                message_prefix: "UFFTYPER: Unrecognized atom type: ",
            }]
        );
        let (rest_length, force_constant) =
            uff_one_u08_source_bond_parameters(0.757, 0.757, 5.343, 5.343, 1.912, 1.912, 1.0);
        for (row, result) in [(&first, results[0]), (&second, results[1])] {
            let distance = (row[1][0] - row[2][0]).abs();
            let independently_computed_energy =
                0.5 * force_constant * (distance - rest_length).powi(2);
            assert_eq!(result.status, 1);
            assert!((result.energy - independently_computed_energy).abs() <= 1.0e-12);
        }
        assert_eq!(first, first_before);
        assert_eq!(second, second_before);
        assert_eq!(reference, reference_before);
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
    }

    #[test]
    fn uff_serial_s12_late_optimizer_failure_preserves_prior_and_later_slots() {
        // FFConvenience.h:66-77 writes each result after initialize/minimize/
        // calcEnergy. The source BFGS invariant failure on a later row leaves
        // prior slots written and later rows/results untouched.
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let properties = MoleculeProperties::default();
        let (rest_length, force_constant) =
            uff_one_u08_source_bond_parameters(0.757, 0.757, 5.343, 5.343, 1.912, 1.912, 1.0);
        let mut reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [rest_length, 0.0, 0.0]]);
        let reference_before = reference.clone();
        let mut first = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut equilibrium = [[0.0, 0.0, 0.0], [rest_length, 0.0, 0.0]];
        let mut later = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let first_before = first;
        let equilibrium_before = equilibrium;
        let later_before = later;
        let mut conformers = [
            SerialConformer {
                id: 80,
                positions: &mut first,
            },
            SerialConformer {
                id: 2,
                positions: &mut equilibrium,
            },
            SerialConformer {
                id: 900,
                positions: &mut later,
            },
        ];
        let sentinel = OptimizationOutcome {
            status: -9,
            energy: 123.0,
        };
        let mut results = [sentinel; 3];
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();

        let error = optimize_serial_uff(
            &topology,
            &mut reference,
            &mut conformers,
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect_err("the second source-ordered row has the zero-gradient invariant failure");
        assert!(matches!(
            error,
            SerialUffOptimizationError::Optimization(
                super::super::convenience::SerialConformerOptimizationError::Optimization {
                    input_index: 1,
                    conformer_id: 2,
                    source: OptimizationStageError::Minimize(
                        ForceFieldKernelError::OptimizerBadDirection
                    ),
                }
            )
        ));
        assert!(diagnostics.is_empty());
        assert_eq!(results[0].status, 0);
        assert!(results[0].energy.is_finite());
        drop(conformers);
        let first_distance = (first[1][0] - first[0][0]).abs();
        let independent_first_energy =
            0.5 * force_constant * (first_distance - rest_length).powi(2);
        assert!((results[0].energy - independent_first_energy).abs() <= 1.0e-12);
        assert_ne!(first[1][0].to_bits(), first_before[1][0].to_bits());
        assert_eq!(
            equilibrium.map(|row| row.map(f64::to_bits)),
            equilibrium_before.map(|row| row.map(f64::to_bits)),
        );
        assert_eq!(
            later.map(|row| row.map(f64::to_bits)),
            later_before.map(|row| row.map(f64::to_bits)),
        );
        assert_eq!(&results[1..], &[sentinel; 2]);
        assert_eq!(reference, reference_before);
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
    }

    #[test]
    fn uff_serial_s13_reuses_source_position_buffer_and_constructs_once() {
        // FFConvenience.h:71-77 overwrites ff.positions() for each row, then
        // initializes and writes results in source traversal order. The
        // source-sized pointer vector and one UFF field survive both rows.
        let topology = uff_one_u09_two_carbon_topology();
        let rings = cf3d_bld_b15_ring_info(&topology);
        let valence = assignment(&[1, 1], &[0, 0]);
        let total_valences = [4, 4];
        let conjugated = [false; 2];
        let properties = MoleculeProperties::default();
        let mut reference = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let reference_before = reference.clone();
        let mut first = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let mut second = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let first_before = first;
        let second_before = second;
        let mut conformers = [
            SerialConformer {
                id: 80,
                positions: &mut first,
            },
            SerialConformer {
                id: 2,
                positions: &mut second,
            },
        ];
        let mut results = [OptimizationOutcome {
            status: -9,
            energy: 123.0,
        }; 2];
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = 0;
        options.ignore_interfragment_interactions = false;
        crate::kernel::cf3d_uff_one_kernel_counts_reset();

        optimize_serial_uff(
            &topology,
            &mut reference,
            &mut conformers,
            &mut results,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &properties,
            &mut diagnostics,
            options,
        )
        .expect("two source-ordered rows reuse one source-built field");

        assert!(diagnostics.is_empty());
        assert_eq!(
            conformers.iter().map(|row| row.id).collect::<Vec<_>>(),
            [80, 2]
        );
        assert_eq!(results.map(|result| result.status), [1, 1]);
        drop(conformers);
        let (rest_length, force_constant) =
            uff_one_u08_source_bond_parameters(0.757, 0.757, 5.343, 5.343, 1.912, 1.912, 1.0);
        for (row, result) in [(&first, results[0]), (&second, results[1])] {
            let distance = (row[1][0] - row[0][0]).abs();
            let expected_energy = 0.5 * force_constant * (distance - rest_length).powi(2);
            assert!((result.energy - expected_energy).abs() <= 1.0e-12);
        }
        assert_eq!(first, first_before);
        assert_eq!(second, second_before);
        assert_eq!(reference, reference_before);
        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));

        let work = crate::kernel::cf3d_uff_one_serial_work_counts();
        assert_eq!(work.position_buffer_growths, 0);
        assert_eq!(work.initialize_calls, 2);
        assert_ne!(work.initial_position_buffer_address, 0);
        assert_eq!(
            work.initial_position_buffer_address,
            work.final_position_buffer_address,
        );
        assert!(work.initial_position_buffer_capacity >= 2);
        assert_eq!(
            work.initial_position_buffer_capacity,
            work.final_position_buffer_capacity,
        );
    }

    #[test]
    fn uff_one_u09_zero_iterations_preserve_source_status_energy_and_coordinates() {
        let topology = uff_one_u09_two_carbon_topology();
        let molecule_properties = MoleculeProperties::default();
        let mut coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let original_coordinates = coordinates.clone();
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = 0;
        options.ignore_interfragment_interactions = false;

        let outcome = uff_one_u09_optimize(
            &topology,
            &mut coordinates,
            &molecule_properties,
            &mut diagnostics,
            options,
        )
        .expect("the source zero-iteration path returns a normal status");

        // Params.cpp C_3 row and BondStretch.cpp:22-46,69-76 give equal
        // endpoint parameters, rest length .757 + .757, and the fixed bond
        // energy at distance 2.0. The bonded pair is excluded from nonbonded.
        let rest_length: f64 = 0.757 + 0.757;
        let force_constant = 2.0 * 332.06 * 1.912 * 1.912 / rest_length.powi(3);
        let displacement = 2.0 - rest_length;
        let expected_energy = 0.5 * force_constant * displacement * displacement;
        assert_eq!(outcome.status, 1);
        assert!((outcome.energy - expected_energy).abs() <= 1.0e-12);
        assert_eq!(coordinates, original_coordinates);
    }

    #[test]
    fn uff_one_u09_preserves_the_selected_builder_error_cause() {
        let topology = uff_one_u09_two_carbon_topology();
        let molecule_properties = MoleculeProperties::default();
        let mut coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let original_coordinates = coordinates.clone();
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(99);
        options.ignore_interfragment_interactions = false;

        let result = uff_one_u09_optimize(
            &topology,
            &mut coordinates,
            &molecule_properties,
            &mut diagnostics,
            options,
        );

        assert!(matches!(
            result,
            Err(SingleConformerOptimizationError::Construction(
                AutomaticForceFieldConstructionError::Construction(
                    ForceFieldConstructionError::Builder(
                        UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                            conformer_id: 99
                        }
                    )
                )
            ))
        ));
        assert_eq!(coordinates, original_coordinates);
    }

    #[test]
    fn uff_one_u09_preserves_the_minimizer_kernel_error_cause() {
        let topology = uff_one_u09_two_carbon_topology();
        let molecule_properties = MoleculeProperties::default();
        let mut coordinates = cf3d_frag_f24_coordinates(&[[0.0, 0.0, 0.0], [1.514, 0.0, 0.0]]);
        let original_coordinates = coordinates.clone();
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.ignore_interfragment_interactions = false;

        let result = uff_one_u09_optimize(
            &topology,
            &mut coordinates,
            &molecule_properties,
            &mut diagnostics,
            options,
        );

        assert!(matches!(
            result,
            Err(SingleConformerOptimizationError::Optimization(
                OptimizationStageError::Minimize(ForceFieldKernelError::OptimizerBadDirection)
            ))
        ));
        assert_eq!(coordinates, original_coordinates);
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum UffOneU10Family {
        Linear,
        BranchedTorsion,
        Salt,
        Aromatic,
        Tbp,
    }

    struct UffOneU10Reference<'a> {
        max_iterations: i32,
        ignore_interfragment_interactions: bool,
        status: i32,
        energy: f64,
        coordinates: &'a [[f64; 3]],
    }

    fn uff_one_u10_fixture(
        family: UffOneU10Family,
    ) -> (
        TopologyBlock,
        Vec<[f64; 3]>,
        Vec<i32>,
        Vec<bool>,
        cosmolkit_core::ValenceAssignment,
    ) {
        match family {
            UffOneU10Family::Linear => (
                uff_one_u09_two_carbon_topology(),
                vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
                vec![4, 4],
                vec![false; 2],
                assignment(&[1, 1], &[0, 0]),
            ),
            UffOneU10Family::BranchedTorsion => {
                let atoms = (0..5)
                    .map(|row| uff_one_u08_atom(row, 6, Hybridization::Sp3, false, 0, true))
                    .collect();
                let topology = uff_one_u08_topology(
                    atoms,
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Single),
                        (1, 3, BondOrder::Single),
                        (3, 4, BondOrder::Single),
                    ],
                );
                (
                    topology,
                    vec![
                        [-1.4, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 1.4, 0.0],
                        [1.4, 0.0, 0.0],
                        [2.1, 1.2, 0.9],
                    ],
                    vec![4; 5],
                    vec![false; 5],
                    assignment(&[1, 3, 1, 2, 1], &[0; 5]),
                )
            }
            UffOneU10Family::Salt => {
                let topology = uff_one_u08_topology(
                    vec![
                        uff_one_u08_atom(0, 11, Hybridization::Unspecified, false, 1, true),
                        uff_one_u08_atom(1, 17, Hybridization::Unspecified, false, -1, true),
                    ],
                    &[],
                );
                (
                    topology,
                    vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
                    vec![0, 0],
                    vec![false; 2],
                    assignment(&[0, 0], &[0, 0]),
                )
            }
            UffOneU10Family::Aromatic => {
                let atoms = (0..6)
                    .map(|row| uff_one_u08_atom(row, 6, Hybridization::Sp2, true, 0, true))
                    .collect();
                let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
                    .into_iter()
                    .enumerate()
                    .map(|(row, (begin, end))| {
                        Bond::from_spec(
                            BondId::new(row),
                            BondSpec::new(
                                AtomId::new(begin),
                                AtomId::new(end),
                                BondOrder::Aromatic,
                            )
                            .with_aromatic(true),
                        )
                    })
                    .collect();
                let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
                    .expect("source aromatic atoms and bonds carry both aromatic flags");
                let h = 0.866_025_403_784_438_6;
                (
                    topology,
                    vec![
                        [1.0, 0.0, 0.0],
                        [0.5, h, 0.0],
                        [-0.5, h, 0.0],
                        [-1.0, 0.0, 0.0],
                        [-0.5, -h, 0.0],
                        [0.5, -h, 0.0],
                    ],
                    vec![4; 6],
                    vec![false; 6],
                    assignment(&[3; 6], &[0; 6]),
                )
            }
            UffOneU10Family::Tbp => {
                let mut atoms = vec![uff_one_u08_atom(0, 15, Hybridization::Sp3d, false, 0, true)];
                atoms.extend(
                    (1..6).map(|row| uff_one_u08_atom(row, 6, Hybridization::Sp3, false, 0, true)),
                );
                let topology = uff_one_u08_topology(
                    atoms,
                    &[
                        (0, 1, BondOrder::Single),
                        (0, 2, BondOrder::Single),
                        (0, 3, BondOrder::Single),
                        (0, 4, BondOrder::Single),
                        (0, 5, BondOrder::Single),
                    ],
                );
                (
                    topology,
                    vec![
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [-0.6, 0.8, 0.0],
                        [0.0, 0.6, 0.8],
                        [0.0, -0.8, 0.6],
                        [-0.8, 0.0, -0.6],
                    ],
                    vec![5, 1, 1, 1, 1, 1],
                    vec![false; 6],
                    assignment(&[5, 1, 1, 1, 1, 1], &[0; 6]),
                )
            }
        }
    }

    fn uff_one_u10_run(
        family: UffOneU10Family,
        max_iterations: i32,
        ignore_interfragment_interactions: bool,
    ) -> (
        super::super::convenience::OptimizationOutcome,
        Vec<[f64; 3]>,
        Vec<UffTypingDiagnostic>,
    ) {
        let (topology, points, total_valences, conjugated, valence) = uff_one_u10_fixture(family);
        let rings = cf3d_bld_b15_ring_info(&topology);
        let mut coordinates = cf3d_frag_f24_coordinates(&points);
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = max_iterations;
        options.ignore_interfragment_interactions = ignore_interfragment_interactions;
        let outcome = optimize_single_conformer(
            &topology,
            &mut coordinates,
            &total_valences,
            &conjugated,
            &rings,
            &valence,
            &MoleculeProperties::default(),
            &mut diagnostics,
            options,
        )
        .unwrap_or_else(|error| {
            panic!("fixed source {family:?} UFF optimization failed: {error:?}")
        });
        let selected = coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == 31)
            .expect("the fixed fixture retains selected source conformer ID31");
        (outcome, selected.coordinates().to_vec(), diagnostics)
    }

    fn assert_uff_one_u10_matrix(
        family: UffOneU10Family,
        references: &[UffOneU10Reference<'_>; 4],
    ) {
        // RDKit testUFFHelpers.cpp:515-517 uses feq's default absolute
        // tolerance; RDGeneral/utils.cpp:38 defines that tolerance comparison.
        const SOURCE_FEQ_TOLERANCE: f64 = 1.0e-4;
        for reference in references {
            let (outcome, actual_coordinates, diagnostics) = uff_one_u10_run(
                family,
                reference.max_iterations,
                reference.ignore_interfragment_interactions,
            );
            assert_eq!(outcome.status, reference.status);
            assert!(
                (outcome.energy - reference.energy).abs() <= SOURCE_FEQ_TOLERANCE,
                "{family:?}, maxIters={}, ignoreInterfrag={}: expected energy {}, got {}",
                reference.max_iterations,
                reference.ignore_interfragment_interactions,
                reference.energy,
                outcome.energy
            );
            assert_eq!(actual_coordinates.len(), reference.coordinates.len());
            for (row, (actual, expected)) in actual_coordinates
                .iter()
                .zip(reference.coordinates)
                .enumerate()
            {
                for axis in 0..3 {
                    assert!(
                        (actual[axis] - expected[axis]).abs() <= SOURCE_FEQ_TOLERANCE,
                        "{family:?}, maxIters={}, ignoreInterfrag={}, row={row}, axis={axis}: expected {}, got {}",
                        reference.max_iterations,
                        reference.ignore_interfragment_interactions,
                        expected[axis],
                        actual[axis]
                    );
                }
            }
            let expected_diagnostics = if family == UffOneU10Family::Tbp {
                vec![UffTypingDiagnostic {
                    atom_id: Some(AtomId::new(0)),
                    kind: UffTypingDiagnosticKind::Warning,
                    message_prefix: "UFFTYPER: Warning: hybridization set to SP3 for atom ",
                }]
            } else {
                Vec::new()
            };
            assert_eq!(diagnostics, expected_diagnostics);
        }
    }

    #[test]
    fn uff_one_u10_linear_source_optimization_matrix() {
        let initial = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let optimized = [
            [0.242_999_999_999_999_77, 0.0, 0.0],
            [1.757_000_000_000_000_1, 0.0, 0.0],
        ];
        assert_uff_one_u10_matrix(
            UffOneU10Family::Linear,
            &[
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: false,
                    status: 1,
                    energy: 82.620_392_244_369_953,
                    coordinates: &initial,
                },
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: true,
                    status: 1,
                    energy: 82.620_392_244_369_953,
                    coordinates: &initial,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: false,
                    status: 0,
                    energy: 1.724_626_936_305_249_4e-29,
                    coordinates: &optimized,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: true,
                    status: 0,
                    energy: 1.724_626_936_305_249_4e-29,
                    coordinates: &optimized,
                },
            ],
        );
    }

    #[test]
    fn uff_one_u10_branched_torsion_source_optimization_matrix() {
        let source_coordinates = [
            [-1.4, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.4, 0.0],
            [1.4, 0.0, 0.0],
            [2.1, 1.2, 0.9],
        ];
        assert_uff_one_u10_matrix(
            UffOneU10Family::BranchedTorsion,
            &[
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: false,
                    status: 1,
                    energy: 150.720_215_359_750_51,
                    coordinates: &source_coordinates,
                },
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: true,
                    status: 1,
                    energy: 150.720_215_359_750_51,
                    coordinates: &source_coordinates,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: false,
                    status: 0,
                    energy: 150.720_122_816_400_25,
                    coordinates: &source_coordinates,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: true,
                    status: 0,
                    energy: 150.720_122_816_400_25,
                    coordinates: &source_coordinates,
                },
            ],
        );
    }

    #[test]
    fn uff_one_u10_salt_source_optimization_matrix() {
        let source_coordinates = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let optimized_with_interactions = [
            [-3.080_850_145_936_969_5, 0.0, 0.0],
            [5.080_850_145_936_974_4, 0.0, 0.0],
        ];
        assert_uff_one_u10_matrix(
            UffOneU10Family::Salt,
            &[
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: false,
                    status: 1,
                    energy: 49.461_474_187_943_104,
                    coordinates: &source_coordinates,
                },
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: true,
                    status: 0,
                    energy: 0.0,
                    coordinates: &source_coordinates,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: false,
                    status: 1,
                    energy: -0.000_908_824_740_464_385_04,
                    coordinates: &optimized_with_interactions,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: true,
                    status: 0,
                    energy: 0.0,
                    coordinates: &source_coordinates,
                },
            ],
        );
    }

    #[test]
    fn uff_one_u10_aromatic_source_optimization_matrix() {
        let h = 0.866_025_403_784_438_6;
        let initial = [
            [1.0, 0.0, 0.0],
            [0.5, h, 0.0],
            [-0.5, h, 0.0],
            [-1.0, 0.0, 0.0],
            [-0.5, -h, 0.0],
            [0.5, -h, 0.0],
        ];
        let optimized = [
            [1.398_526_890_296_198_3, -1.568_066_543_708_494_3e-9, 0.0],
            [0.699_263_300_991_334_16, 1.211_159_771_947_467_8, 0.0],
            [-0.699_263_302_310_720_42, 1.211_159_770_386_722_7, 0.0],
            [-1.398_526_887_700_742_4, -1.430_662_258_580_250_6e-9, 0.0],
            [-0.699_263_306_136_564_4, -1.211_159_768_987_247_3, 0.0],
            [0.699_263_304_860_495_59, -1.211_159_770_348_214_2, 0.0],
        ];
        assert_uff_one_u10_matrix(
            UffOneU10Family::Aromatic,
            &[
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: false,
                    status: 1,
                    energy: 1185.319_763_568_795_9,
                    coordinates: &initial,
                },
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: true,
                    status: 1,
                    energy: 1185.319_763_568_795_9,
                    coordinates: &initial,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: false,
                    status: 0,
                    energy: 11.354_132_646_913_435,
                    coordinates: &optimized,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: true,
                    status: 0,
                    energy: 11.354_132_646_913_435,
                    coordinates: &optimized,
                },
            ],
        );
    }

    #[test]
    fn uff_one_u10_tbp_source_optimization_matrix() {
        let initial = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-0.6, 0.8, 0.0],
            [0.0, 0.6, 0.8],
            [0.0, -0.8, 0.6],
            [-0.8, 0.0, -0.6],
        ];
        let optimized = [
            [
                0.061_731_843_648_203_187,
                -0.161_020_032_109_450_51,
                -0.222_038_857_274_712_59,
            ],
            [
                1.878_168_286_799_183_1,
                -0.334_869_680_420_552_71,
                -0.021_858_341_307_527_486,
            ],
            [
                -0.643_576_809_023_930_33,
                1.438_009_642_745_417_6,
                0.310_574_265_889_973_39,
            ],
            [
                0.188_918_662_711_506_98,
                1.354_535_155_844_982_1,
                0.735_992_240_341_444_84,
            ],
            [
                -0.207_137_701_374_027_48,
                -1.569_881_871_849_908_4,
                0.989_689_955_013_540_29,
            ],
            [
                -1.678_104_282_761_317_7,
                -0.126_773_214_208_941_58,
                -0.992_359_262_661_794_64,
            ],
        ];
        assert_uff_one_u10_matrix(
            UffOneU10Family::Tbp,
            &[
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: false,
                    status: 1,
                    energy: 1142.363_631_778_066_3,
                    coordinates: &initial,
                },
                UffOneU10Reference {
                    max_iterations: 0,
                    ignore_interfragment_interactions: true,
                    status: 1,
                    energy: 1142.363_631_778_066_3,
                    coordinates: &initial,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: false,
                    status: 0,
                    energy: 72.914_702_361_769_429,
                    coordinates: &optimized,
                },
                UffOneU10Reference {
                    max_iterations: 200,
                    ignore_interfragment_interactions: true,
                    status: 0,
                    energy: 72.914_702_361_769_429,
                    coordinates: &optimized,
                },
            ],
        );
    }

    #[test]
    fn uff_one_u11_selected_path_builds_one_field_and_reuses_its_term() {
        const UFF_ONE_SOURCE_FEQ_TOLERANCE: f64 = 1.0e-4;
        // The pinned linear C_3-C_3 reference has one bond term, no accepted
        // nonbonded pair, no torsion, and no inversion (UFF.h plus Builder.cpp
        // source construction order). Disabling interfragment interactions
        // takes the source early branch and avoids the required F18 component
        // copies; the selected API still receives borrowed detached blocks.
        let topology = uff_one_u09_two_carbon_topology();
        let molecule_properties = MoleculeProperties::default();
        let initial = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let expected = [
            [0.242_999_999_999_999_77, 0.0, 0.0],
            [1.757_000_000_000_000_1, 0.0, 0.0],
        ];
        let mut coordinates = cf3d_frag_f24_coordinates(&initial);
        let mut diagnostics = Vec::new();
        let mut options = SingleConformerOptions::for_conformer(31);
        options.max_iterations = 200;
        options.ignore_interfragment_interactions = false;

        crate::kernel::cf3d_uff_one_kernel_counts_reset();
        let outcome = uff_one_u09_optimize(
            &topology,
            &mut coordinates,
            &molecule_properties,
            &mut diagnostics,
            options,
        )
        .expect("the pinned fixed linear input has a UFF field");

        assert_eq!(crate::kernel::cf3d_uff_one_kernel_counts(), (1, 1, 0));
        assert_eq!(outcome.status, 0);
        assert!(
            (outcome.energy - 1.724_626_936_305_249_4e-29).abs() <= UFF_ONE_SOURCE_FEQ_TOLERANCE
        );
        let unselected = coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == 12)
            .expect("source ID12 remains present");
        assert_eq!(unselected.coordinates()[0], [20.0, 0.0, 0.0]);
        let actual = coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == 31)
            .expect("the optimizer borrows selected source ID31")
            .coordinates();
        for (actual, expected) in actual.iter().zip(expected) {
            for axis in 0..3 {
                assert!(
                    (actual[axis] - expected[axis]).abs() <= UFF_ONE_SOURCE_FEQ_TOLERANCE,
                    "row {axis}: expected {expected:?}, got {actual:?}"
                );
            }
        }
        assert!(diagnostics.is_empty());
    }
}
