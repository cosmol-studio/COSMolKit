//! MolBlock writer implementation.

use cosmolkit_model::{PropertyText, PropertyValue};
use std::borrow::Cow;

use cosmolkit_core::{RingInfo, ValenceAssignment, WedgeAssignments};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, Bond, BondDirection, BondId, BondOrder, BondQueryPredicate,
    BondStereo, ChiralTag, CoordinateDimension, QueryNode, SGroupBondRole, SGroupBracketStyle,
    SGroupConnection, SGroupData, StereoGroupKind, SubstanceGroup, SubstanceGroupKind,
};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

/// Borrowed detached format input; no cache installation or live authority.
#[derive(Clone, Copy)]
pub struct MolWriteInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    pub rings: Option<&'a RingInfo>,
}

#[derive(Clone)]
struct MolWriteContext<'a> {
    topology: Cow<'a, TopologyBlock>,
    coordinates: Cow<'a, CoordinateBlock>,
    properties: &'a MoleculeProperties,
    valence: Option<Cow<'a, ValenceAssignment>>,
    rings: Option<&'a RingInfo>,
    query: Option<&'a cosmolkit_model::QueryGraph>,
}
impl MolWriteContext<'_> {
    fn query_atom(&self, atom: AtomId) -> Option<&cosmolkit_model::QueryAtom> {
        self.query
            .and_then(|q| q.atoms().get(atom.index()))
            .filter(|row| !row.predicate_is_carrier_derived())
    }
    fn query_bond(&self, bond: BondId) -> Option<&cosmolkit_model::QueryBond> {
        self.query
            .and_then(|q| q.bonds().get(bond.index()))
            .filter(|row| !row.predicate_is_carrier_derived())
    }
    fn atoms(&self) -> &[Atom] {
        &self.topology.atoms
    }
    fn bonds(&self) -> &[Bond] {
        &self.topology.bonds
    }
    fn num_atoms(&self) -> usize {
        self.topology.atoms.len()
    }
    fn num_bonds(&self) -> usize {
        self.topology.bonds.len()
    }
    fn topology_block(&self) -> &TopologyBlock {
        &self.topology
    }
    fn substance_groups(&self) -> &[SubstanceGroup] {
        &self.topology.substance_groups
    }
    fn stereo_groups(&self) -> &[cosmolkit_model::StereoGroup] {
        &self.topology.stereo_groups
    }
    fn properties(&self) -> &MoleculeProperties {
        self.properties
    }
    fn prop(&self, key: &str) -> Option<&PropertyValue> {
        self.properties.prop(key)
    }
    fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.coordinates
            .conformers_2d
            .first()
            .map(|c| c.coordinates())
    }
    fn conformers_3d(&self) -> &[cosmolkit_model::Conformer3D] {
        &self.coordinates.conformers_3d
    }
    fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.coordinates.source_coordinate_dim
    }
}
pub fn write_mol_block_with_params(
    data: MolWriteInput<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let MolWriteInput {
        topology,
        coordinates,
        properties,
        rings,
    } = data;
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let input = MolWriteContext {
        topology: Cow::Borrowed(topology),
        coordinates: Cow::Borrowed(coordinates),
        properties,
        valence: None,
        rings,
        query: None,
    };
    mol_to_mol_block_with_params(&input, params)
}
pub fn write_sdf_with_params(
    data: MolWriteInput<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let MolWriteInput {
        topology,
        coordinates,
        properties,
        rings,
    } = data;
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let input = MolWriteContext {
        topology: Cow::Borrowed(topology),
        coordinates: Cow::Borrowed(coordinates),
        properties,
        valence: None,
        rings,
        query: None,
    };
    mol_to_sdf_record_with_params(&input, params)
}
pub fn write_sdf_3d_with_params(
    data: MolWriteInput<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let MolWriteInput {
        topology,
        coordinates,
        properties,
        rings,
    } = data;
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let input = MolWriteContext {
        topology: Cow::Borrowed(topology),
        coordinates: Cow::Borrowed(coordinates),
        properties,
        valence: None,
        rings,
        query: None,
    };
    let selection = export_selection(params, Some(CoordinateDimension::ThreeD))?;
    let block = match params.format {
        SdfFormat::V2000 => mol_to_v2000_block_with_params(&input, selection, params)?,
        SdfFormat::V3000 => mol_to_v3000_block_with_params(&input, selection, params)?,
    };
    Ok(append_sdf_record_fields(block, &input))
}

/// Explicit detached query record input for the existing MOL/SDF formatter.
#[derive(Clone, Copy)]
pub struct QueryMolWriteInput<'a> {
    pub query: &'a cosmolkit_model::QueryGraph,
    pub properties: &'a MoleculeProperties,
    pub source_coordinate_dim: Option<CoordinateDimension>,
    pub rings: Option<&'a RingInfo>,
}

fn query_write_context(data: QueryMolWriteInput<'_>) -> Result<MolWriteContext<'_>, MolWriteError> {
    // RDKit❗❌:   RWMol trwmol(mol);
    // The source writer uses one mutable scratch copy. Here the existing
    // detached carrier conversion retains the separately borrowed query AST;
    // it never constructs a concrete Molecule or installs live derived state.
    // One O(V+E) scratch allocation plus coordinate cloning matches that copy.
    data.query.validate()?;
    let topology = TopologyBlock::try_from_parts(
        data.query
            .atoms()
            .iter()
            .map(cosmolkit_model::QueryAtom::try_to_atom)
            .collect::<Result<Vec<_>, _>>()?,
        data.query
            .bonds()
            .iter()
            .map(|row| row.bond().clone())
            .collect(),
        cosmolkit_model::query_substance_groups(data.query).to_vec(),
        data.query.stereo_groups().to_vec(),
    )?;
    cosmolkit_model::QueryStateRef::try_for_topology(
        data.query.atoms(),
        data.query.bonds(),
        &topology,
    )?;
    let coordinates = data.query.coordinate_block(data.source_coordinate_dim);
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    Ok(MolWriteContext {
        topology: Cow::Owned(topology),
        coordinates: Cow::Owned(coordinates),
        properties: data.properties,
        valence: None,
        rings: data.rings,
        query: Some(data.query),
    })
}

pub fn write_query_mol_block_with_params(
    data: QueryMolWriteInput<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let input = query_write_context(data)?;
    mol_to_mol_block_with_params(&input, params)
}

pub fn write_query_sdf_with_params(
    data: QueryMolWriteInput<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let input = query_write_context(data)?;
    mol_to_sdf_record_with_params(&input, params)
}

const MIN_V2000_COORD: f64 = -10_000.0;
const MAX_V2000_COORD: f64 = 100_000.0;

#[derive(Debug, thiserror::Error)]
pub enum MolWriteError {
    #[error("MolBlock writing subset is not supported: {0}")]
    UnsupportedSubset(&'static str),
    #[error("MolBlock writing failed: {0}")]
    Value(String),
    #[error(
        "MOL/SDF coordinate selection is ambiguous: {two_d} 2D layouts and {three_d} 3D conformers"
    )]
    AmbiguousCoordinates { two_d: usize, three_d: usize },
    #[error("MOL/SDF selected {dimension:?} coordinate ID {id} is missing")]
    MissingCoordinate {
        dimension: CoordinateDimension,
        id: usize,
    },
    #[error("MOL/SDF coordinate selection conflicts with the requested output dimension")]
    CoordinateDimensionMismatch,
    #[error(transparent)]
    QueryGraph(#[from] cosmolkit_model::QueryGraphError),
    #[error(transparent)]
    QueryAtom(#[from] cosmolkit_model::QueryAtomConversionError),
    #[error(transparent)]
    QueryState(#[from] cosmolkit_model::QueryStateError),
    #[error(transparent)]
    QuerySmarts(#[from] cosmolkit_search::SmartsWriteError),
    #[error("MolBlock valence assignment failed: {0}")]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error(transparent)]
    Kekulize(#[from] cosmolkit_core::KekulizeError),
    #[error(transparent)]
    Atropisomer(#[from] cosmolkit_core::AtropisomerError),
    #[error(transparent)]
    Wedge(#[from] cosmolkit_core::WedgeError),
    #[error(transparent)]
    Depict(#[from] cosmolkit_depict::DepictError),
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error(transparent)]
    Property(#[from] cosmolkit_model::PropertyValueError),
    #[error(transparent)]
    UnsignedProperty(#[from] cosmolkit_core::PropertyUIntReadError),
    #[error(transparent)]
    Detached(#[from] crate::SdfWriteError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdfFormat {
    V2000,
    V3000,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MolBlockWriteParams {
    pub format: SdfFormat,
    pub force_2d: bool,
    pub include_stereo: bool,
    pub kekulize: bool,
    pub precision: usize,
    pub coordinate_selection: MolCoordinateSelection,
    pub include_coordinates: bool,
}

impl Default for MolBlockWriteParams {
    fn default() -> Self {
        Self {
            format: SdfFormat::V2000,
            force_2d: false,
            include_stereo: true,
            kekulize: true,
            precision: 6,
            coordinate_selection: MolCoordinateSelection::Auto,
            include_coordinates: true,
        }
    }
}

/// Canonical unique-or-explicit geometry selection for MOL and SDF export.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum MolCoordinateSelection {
    #[default]
    Auto,
    TwoD {
        id: usize,
    },
    ThreeD {
        id: usize,
    },
}
#[derive(Clone, Copy)]
enum CoordinateSelection {
    Auto,
    TwoD(Option<usize>),
    ThreeD(Option<usize>),
    Disabled,
}
fn export_selection(
    params: &MolBlockWriteParams,
    dimension: Option<CoordinateDimension>,
) -> Result<CoordinateSelection, MolWriteError> {
    if !params.include_coordinates {
        return Ok(CoordinateSelection::Disabled);
    }
    match (dimension, params.coordinate_selection) {
        (Some(CoordinateDimension::TwoD), MolCoordinateSelection::ThreeD { .. })
        | (Some(CoordinateDimension::ThreeD), MolCoordinateSelection::TwoD { .. }) => {
            Err(MolWriteError::CoordinateDimensionMismatch)
        }
        (_, MolCoordinateSelection::TwoD { id }) => Ok(CoordinateSelection::TwoD(Some(id))),
        (_, MolCoordinateSelection::ThreeD { id }) => Ok(CoordinateSelection::ThreeD(Some(id))),
        (Some(CoordinateDimension::TwoD), MolCoordinateSelection::Auto) => {
            Ok(CoordinateSelection::TwoD(None))
        }
        (Some(CoordinateDimension::ThreeD), MolCoordinateSelection::Auto) => {
            Ok(CoordinateSelection::ThreeD(None))
        }
        (None, MolCoordinateSelection::Auto) => Ok(if params.force_2d {
            CoordinateSelection::TwoD(None)
        } else {
            CoordinateSelection::Auto
        }),
    }
}

struct SelectedCoordinates {
    coords: Option<Vec<[f64; 3]>>,
    is_3d: bool,
    selection: Option<MolCoordinateSelection>,
}

fn mol_to_v2000_2d_block(molecule: &MolWriteContext<'_>) -> Result<PropertyText, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        force_2d: true,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::TwoD(None), &params)
}

fn mol_to_v2000_3d_block(molecule: &MolWriteContext<'_>) -> Result<PropertyText, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::ThreeD(None), &params)
}

fn mol_to_v2000_block(molecule: &MolWriteContext<'_>) -> Result<PropertyText, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::Auto, &params)
}

fn mol_to_v3000_block(molecule: &MolWriteContext<'_>) -> Result<PropertyText, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::Auto, &params)
}

fn mol_to_v3000_2d_block(molecule: &MolWriteContext<'_>) -> Result<PropertyText, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        force_2d: true,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::TwoD(None), &params)
}

fn mol_to_v3000_3d_block(molecule: &MolWriteContext<'_>) -> Result<PropertyText, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::ThreeD(None), &params)
}

fn mol_to_mol_block_with_params(
    molecule: &MolWriteContext<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let selection = export_selection(params, None)?;
    match params.format {
        SdfFormat::V2000 => mol_to_v2000_block_with_params(molecule, selection, params),
        SdfFormat::V3000 => mol_to_v3000_block_with_params(molecule, selection, params),
    }
}

fn mol_to_2d_sdf_record(
    molecule: &MolWriteContext<'_>,
    format: SdfFormat,
) -> Result<PropertyText, MolWriteError> {
    let block = match format {
        SdfFormat::V2000 => mol_to_v2000_2d_block(molecule)?,
        SdfFormat::V3000 => mol_to_v3000_2d_block(molecule)?,
    };
    Ok(append_sdf_record_fields(block, molecule))
}

fn mol_to_sdf_record_with_params(
    molecule: &MolWriteContext<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let block = mol_to_mol_block_with_params(molecule, params)?;
    Ok(append_sdf_record_fields(block, molecule))
}

fn mol_to_3d_sdf_record(
    molecule: &MolWriteContext<'_>,
    format: SdfFormat,
) -> Result<PropertyText, MolWriteError> {
    let block = match format {
        SdfFormat::V2000 => mol_to_v2000_3d_block(molecule)?,
        SdfFormat::V3000 => mol_to_v3000_3d_block(molecule)?,
    };
    Ok(append_sdf_record_fields(block, molecule))
}

fn should_auto_upgrade_to_v3000(
    molecule: &MolWriteContext<'_>,
    selected: &SelectedCoordinates,
) -> bool {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: std::string outputMolToMolBlock(const RWMol &tmol, int confId,
    // RDKit❗❌:                                 MolFileFormat whichFormat,
    // RDKit❗❌:                                 unsigned int precision,
    // RDKit❗❌:                                 const boost::dynamic_bitset<> &aromaticBonds) {
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   unsigned int nAtoms, nBonds, nLists, chiralFlag, nsText, nRxnComponents;
    // RDKit❗❌:   unsigned int nReactants, nProducts, nIntermediates;
    // RDKit❗❌:   nAtoms = tmol.getNumAtoms();
    // RDKit❗❌:   nBonds = tmol.getNumBonds();
    // RDKit❗❌:   nLists = 0;
    // RDKit❗❌:
    // RDKit❗❌:   const auto &sgroups = getSubstanceGroups(tmol);
    // RDKit❗❌:   unsigned int nSGroups = sgroups.size();
    // RDKit❗❌:
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V2000 &&
    // RDKit❗❌:       (nAtoms > 999 || nBonds > 999 || nSGroups > 999)) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "V2000 format does not support more than 999 atoms, bonds or SGroups.");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   chiralFlag = 0;
    // RDKit❗❌:   nsText = 0;
    // RDKit❗❌:   nRxnComponents = 0;
    // RDKit❗❌:   nReactants = 0;
    // RDKit❗❌:   nProducts = 0;
    // RDKit❗❌:   nIntermediates = 0;
    // RDKit❗❌:
    // RDKit❗❌:   tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit❗❌:
    // RDKit❗❌:   const Conformer *conf;
    // RDKit❗❌:   if (confId < 0 && tmol.getNumConformers() == 0) {
    // RDKit❗❌:     conf = nullptr;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     conf = &(tmol.getConformer(confId));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool coordMagnitudeTooLargeForV2K = false;
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     for (auto &pos : conf->getPositions()) {
    // RDKit❗❌:       if ((pos.x >= MAX_V2000_COORD || pos.x <= MIN_V2000_COORD) ||
    // RDKit❗❌:           (pos.y >= MAX_V2000_COORD || pos.y <= MIN_V2000_COORD) ||
    // RDKit❗❌:           (pos.z >= MAX_V2000_COORD || pos.z <= MIN_V2000_COORD)) {
    // RDKit❗❌:         coordMagnitudeTooLargeForV2K = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V2000 && coordMagnitudeTooLargeForV2K) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "V2000 format does not support atom positions <= " +
    // RDKit❗❌:         std::to_string((int)MIN_V2000_COORD) +
    // RDKit❗❌:         " or >= " + std::to_string((int)MAX_V2000_COORD));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string text;
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:
    // RDKit❗❌:   // info
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     std::stringstream ss;
    // RDKit❗❌:     ss << "  " << std::setw(8) << "RDKit";
    // RDKit❗❌:     ss << std::setw(10) << "";
    // RDKit❗❌:     if (conf) {
    // RDKit❗❌:       if (conf->is3D()) {
    // RDKit❗❌:         ss << "3D";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         ss << common_properties::TWOD;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     res += ss.str();
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:   // comments
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:
    // RDKit❗❌:   bool hasDative = false;
    // RDKit❗❌:   for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:     if (bond->getBondType() == Bond::DATIVE) {
    // RDKit❗❌:       hasDative = true;
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool isV3000 = false;
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V3000) {
    // RDKit❗❌:     isV3000 = true;
    // RDKit❗❌:   } else if (whichFormat == MolFileFormat::unspecified &&
    // RDKit❗❌:              (coordMagnitudeTooLargeForV2K || hasDative || nAtoms > 999 ||
    // RDKit❗❌:               nBonds > 999 || nSGroups > 999 ||
    // RDKit❗❌:               !tmol.getStereoGroups().empty())) {
    // RDKit❗❌:     isV3000 = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // the counts line:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   if (isV3000) {
    // RDKit❗❌:     // All counts in the V3000 info line should be 0
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << "999 V3000\n";
    // RDKit❗❌:   } else {
    // RDKit❗❌:     ss << std::setw(3) << nAtoms;
    // RDKit❗❌:     ss << std::setw(3) << nBonds;
    // RDKit❗❌:     ss << std::setw(3) << nLists;
    // RDKit❗❌:     ss << std::setw(3) << nSGroups;
    // RDKit❗❌:     ss << std::setw(3) << chiralFlag;
    // RDKit❗❌:     ss << std::setw(3) << nsText;
    // RDKit❗❌:     ss << std::setw(3) << nRxnComponents;
    // RDKit❗❌:     ss << std::setw(3) << nReactants;
    // RDKit❗❌:     ss << std::setw(3) << nProducts;
    // RDKit❗❌:     ss << std::setw(3) << nIntermediates;
    // RDKit❗❌:     ss << "999 V2000\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   res += ss.str();
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> queryListAtoms(tmol.getNumAtoms());
    // RDKit❗❌:   if (!isV3000) {
    // RDKit❗❌:     // V2000 output.
    // RDKit❗❌:     for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
    // RDKit❗❌:          atomIt != tmol.endAtoms(); ++atomIt) {
    // RDKit❗❌:       res += GetMolFileAtomLine(*atomIt, conf, queryListAtoms);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // RDKit❗❌:
    // RDKit❗❌:     for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:       res += GetMolFileBondLine(bond, wedgeBonds, conf,
    // RDKit❗❌:                                 aromaticBonds[bond->getIdx()]);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     res += GetMolFileChargeInfo(tmol);
    // RDKit❗❌:     res += GetMolFileRGroupInfo(tmol);
    // RDKit❗❌:     res += GetMolFileQueryInfo(tmol, queryListAtoms);
    // RDKit❗❌:     res += GetMolFileAliasInfo(tmol);
    // RDKit❗❌:     res += GetMolFileZBOInfo(tmol);
    // RDKit❗❌:
    // RDKit❗❌:     res += GetMolFilePXAInfo(tmol);
    // RDKit❗❌:     res += GetMolFileSGroupInfo(tmol);
    // RDKit❗❌:
    // RDKit❗❌:     // FIX: R-group logic, SGroups and 3D features etc.
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // V3000 output.
    // RDKit❗❌:     res +=
    // RDKit❗❌:         FileParserUtils::getV3000CTAB(tmol, aromaticBonds, confId, precision);
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "M  END\n";
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    molecule
        .bonds()
        .iter()
        .any(|bond| bond.order() == BondOrder::Dative)
        || molecule.num_atoms() > 999
        || molecule.num_bonds() > 999
        || molecule.substance_groups().len() > 999
        || !molecule.stereo_groups().is_empty()
        || selected.coords.as_deref().is_some_and(|rows| {
            rows.iter()
                .flatten()
                .any(|value| *value <= MIN_V2000_COORD || *value >= MAX_V2000_COORD)
        })
}

/// Prepared molecule with aromatic-bond bookkeeping.
/// Before kekulization, the set of bonds that were aromatic is recorded so
/// crossed-bond stereo output can be suppressed for those bonds after
/// kekulization (matching RDKit's `prepareMol` → `aromaticBonds`
/// bookkeeping).
struct PreparedMol<'a> {
    molecule: MolWriteContext<'a>,
    /// Indices (in the bond table) of bonds that were aromatic before
    /// kekulization. After kekulization, these should still be written
    /// as bond type 4 (aromatic) in the molfile output.
    aromatic_bonds: Vec<usize>,
    wedge_bonds: WedgeAssignments,
    selected: SelectedCoordinates,
    ring_info: RingInfo,
}

#[derive(Clone, Copy)]
struct MolfileStereoContext<'a> {
    crossed: cosmolkit_core::CrossedBondContext<'a>,
    conformer: Option<cosmolkit_core::AtropisomerConformer<'a>>,
}

fn prepare_mol_for_writing<'a>(
    molecule: &'a MolWriteContext<'a>,
    params: &MolBlockWriteParams,
    selection: CoordinateSelection,
) -> Result<PreparedMol<'a>, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: void prepareMol(RWMol &trwmol, const MolWriterParams &params,
    // RDKit❗❌:                 boost::dynamic_bitset<> &aromaticBonds) {
    // RDKit❗❌:   // NOTE: kekulize the molecule before writing it out
    // RDKit❗❌:   // because of the way mol files handle aromaticity
    // RDKit❗❌:   if (trwmol.needsUpdatePropertyCache()) {
    // RDKit❗❌:     trwmol.updatePropertyCache(false);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (params.kekulize && trwmol.getNumBonds()) {
    // RDKit❗❌:     for (const auto bond : trwmol.bonds()) {
    // RDKit❗❌:       if (bond->getIsAromatic()) {
    // RDKit❗❌:         aromaticBonds.set(bond->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     MolOps::Kekulize(trwmol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (params.includeStereo && !trwmol.getNumConformers()) {
    // RDKit❗❌:     // generate coordinates so that the stereo we generate makes sense
    // RDKit❗❌:     RDDepict::compute2DCoords(trwmol);
    // RDKit❗❌:   }
    // RDKit❗❌:   FileParserUtils::moveAdditionalPropertiesToSGroups(trwmol);
    // RDKit❗❌: }

    let aromatic_bonds = molecule
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, b)| b.is_aromatic() || b.order() == BondOrder::Aromatic)
        .map(|(i, _)| i)
        .collect();
    let mut mol = molecule.clone();
    if params.kekulize && molecule.num_bonds() != 0 {
        mol.topology = Cow::Owned(
            cosmolkit_core::kekulize_with_query_state(
                &molecule.topology,
                &cosmolkit_core::KekulizeParams::default(),
                molecule
                    .query
                    .map(|q| {
                        cosmolkit_model::QueryStateRef::try_for_topology(
                            q.atoms(),
                            q.bonds(),
                            &molecule.topology,
                        )
                    })
                    .transpose()?,
            )?
            .topology,
        );
    }
    // RDKit✔️✔️:   if (params.includeStereo && !trwmol.getNumConformers()) {
    // RDKit✔️✔️:     // generate coordinates so that the stereo we generate makes sense
    // RDKit✔️✔️:     RDDepict::compute2DCoords(trwmol);
    // RDKit✔️✔️:   }
    // Restore source prepareMol behavior on the temporary writer value. This
    // does not install geometry on the input or change unique/explicit
    // selection for existing conformers. CK's explicit missing-ID and
    // coordinate-disabled options retain their declared behavior.
    if params.include_stereo
        && mol.coordinates.conformers_2d.is_empty()
        && mol.coordinates.conformers_3d.is_empty()
        && matches!(
            selection,
            CoordinateSelection::Auto | CoordinateSelection::TwoD(None)
        )
    {
        let conformer = cosmolkit_depict::compute_2d_coordinates(
            &mol.topology,
            mol.properties,
            &Default::default(),
        )?;
        let coordinates = mol.coordinates.to_mut();
        coordinates.record_source_conformer_append(CoordinateDimension::TwoD)?;
        coordinates.conformers_2d.push(conformer);
    }
    let selected = select_coordinates(&mol, selection)?;
    let conformer = writer_conformer(&mol, &selected);
    let (wedge_bonds, ring_info) = cosmolkit_core::pick_bonds_to_wedge_with_existing_ring_info(
        &mol.topology,
        conformer,
        mol.rings.cloned(),
    )?;
    let valence = molblock_valence_assignment(&mol)?.into_owned();
    mol.valence = Some(Cow::Owned(valence));
    Ok(PreparedMol {
        molecule: mol,
        aromatic_bonds,
        wedge_bonds,
        selected,
        ring_info,
    })
}

fn writer_conformer<'a>(
    molecule: &'a MolWriteContext<'_>,
    selected: &SelectedCoordinates,
) -> Option<cosmolkit_core::AtropisomerConformer<'a>> {
    // The canonical CORE wedge APIs borrow the selected storage carrier itself,
    // preserving Conformer3D::is_3d() even when it is false. Do not construct a
    // flag-forced conformer or turn its XYZ rows into a 2D carrier: atrop wedging
    // branches on the stored bit while tetrahedral wedging uses source XY.
    // Existing dimension-local ID lookup is O(C), with no row allocation.
    match selected.selection {
        Some(MolCoordinateSelection::TwoD { id }) => molecule
            .coordinates
            .conformers_2d
            .iter()
            .find(|c| c.id() == id)
            .map(cosmolkit_core::AtropisomerConformer::TwoD),
        Some(MolCoordinateSelection::ThreeD { id }) => molecule
            .coordinates
            .conformers_3d
            .iter()
            .find(|c| c.id() == id)
            .map(cosmolkit_core::AtropisomerConformer::ThreeD),
        _ => None,
    }
}

fn mol_to_v3000_block_with_params(
    molecule: &MolWriteContext<'_>,
    selection: CoordinateSelection,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: std::string getV3000CTAB(const ROMol &tmol,
    // RDKit❗❌:                          const boost::dynamic_bitset<> &wasAromatic, int confId,
    // RDKit❗❌:                          unsigned int precision) {
    // RDKit❗❌:   auto nAtoms = tmol.getNumAtoms();
    // RDKit❗❌:   auto nBonds = tmol.getNumBonds();
    // RDKit❗❌:   const auto &sgroups = getSubstanceGroups(tmol);
    // RDKit❗❌:   auto nSGroups = sgroups.size();
    // RDKit❗❌:
    // RDKit❗❌:   unsigned chiralFlag = 0;
    // RDKit❗❌:   tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit❗❌:
    // RDKit❗❌:   const Conformer *conf = nullptr;
    // RDKit❗❌:   if (confId >= 0 || tmol.getNumConformers()) {
    // RDKit❗❌:     conf = &(tmol.getConformer(confId));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string res = "M  V30 BEGIN CTAB\n";
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   int num3DConstraints = 0;  //< not implemented
    // RDKit❗❌:   ss << "M  V30 COUNTS " << nAtoms << " " << nBonds << " " << nSGroups << " "
    // RDKit❗❌:      << num3DConstraints << " " << chiralFlag << "\n";
    // RDKit❗❌:
    // RDKit❗❌:   res += ss.str();
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> queryListAtoms(tmol.getNumAtoms());
    // RDKit❗❌:   res += "M  V30 BEGIN ATOM\n";
    // RDKit❗❌:   for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
    // RDKit❗❌:        atomIt != tmol.endAtoms(); ++atomIt) {
    // RDKit❗❌:     res += GetV3000MolFileAtomLine(*atomIt, conf, queryListAtoms, precision);
    // RDKit❗❌:     res += "\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "M  V30 END ATOM\n";
    // RDKit❗❌:
    // RDKit❗❌:   auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // RDKit❗❌:   if (tmol.getNumBonds()) {
    // RDKit❗❌:     res += "M  V30 BEGIN BOND\n";
    // RDKit❗❌:
    // RDKit❗❌:     for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:       res += GetV3000MolFileBondLine(bond, wedgeBonds, conf,
    // RDKit❗❌:                                      wasAromatic[bond->getIdx()]);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:     res += "M  V30 END BOND\n";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (nSGroups > 0) {
    // RDKit❗❌:     res += "M  V30 BEGIN SGROUP\n";
    // RDKit❗❌:     unsigned int idx = 0;
    // RDKit❗❌:     for (const auto &sgroup : sgroups) {
    // RDKit❗❌:       res += GetV3000MolFileSGroupLines(++idx, sgroup);
    // RDKit❗❌:     }
    // RDKit❗❌:     res += "M  V30 END SGROUP\n";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (tmol.hasProp(common_properties::molFileLinkNodes)) {
    // RDKit❗❌:     auto pval = tmol.getProp<std::string>(common_properties::molFileLinkNodes);
    // RDKit❗❌:
    // RDKit❗❌:     std::vector<std::string> linknodes;
    // RDKit❗❌:     boost::split(linknodes, pval, boost::is_any_of("|"));
    // RDKit❗❌:     for (const auto &linknode : linknodes) {
    // RDKit❗❌:       res += "M  V30 LINKNODE " + linknode + "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   appendEnhancedStereoGroups(res, tmol, wedgeBonds);
    // RDKit❗❌:
    // RDKit❗❌:   res += "M  V30 END CTAB\n";
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    let prepared = prepare_mol_for_writing(molecule, params, selection)?;
    render_v3000(&prepared, params)
}

fn render_v3000(
    prepared: &PreparedMol<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    // RDKit❗❌: void GetMolFileAtomProperties(const Atom *atom, const Conformer *conf,
    // RDKit❗❌:                               int &totValence, int &atomMapNumber,
    // RDKit❗❌:                               unsigned int &parityFlag, double &x, double &y,
    // RDKit❗❌:                               double &z) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   totValence = 0;
    // RDKit❗❌:   atomMapNumber = 0;
    // RDKit❗❌:   parityFlag = 0;
    // RDKit❗❌:   x = y = z = 0.0;
    // RDKit❗❌:
    // RDKit❗❌:   if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                               atomMapNumber)) {
    // RDKit❗❌:     // XXX FIX ME->should we fail here? previously we would not assign
    // RDKit❗❌:     // the atomMapNumber if it didn't exist which could result in garbage
    // RDKit❗❌:     //  values.
    // RDKit❗❌:     atomMapNumber = 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:     x = pos.x;
    // RDKit❗❌:     y = pos.y;
    // RDKit❗❌:     z = pos.z;
    // RDKit❗❌:     if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:         atom->getChiralTag() != Atom::CHI_OTHER && atom->getDegree() >= 3 &&
    // RDKit❗❌:         atom->getTotalDegree() == 4) {
    // RDKit❗❌:       parityFlag = getAtomParityFlag(atom, conf);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (hasNonDefaultValence(atom)) {
    // RDKit❗❌:     if (atom->getTotalDegree() == 0) {
    // RDKit❗❌:       // Specify zero valence for elements/metals without neighbors
    // RDKit❗❌:       // or hydrogens (degree 0) instead of writing them as radicals.
    // RDKit❗❌:       totValence = 15;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // write the total valence for other atoms
    // RDKit❗❌:       totValence = atom->getTotalValence() % 15;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Selection supplies both XYZ rows and the independent stored is3D bit.
    // The source always copies x/y/z; only parity is gated by is3D. The chiral
    // tag and degree guards live in the existing parity helper. Wedge selection
    // and rendering borrow this same source conformer through writer_conformer,
    // retaining its flag and full rows for the canonical CORE helper branches.
    // These guards are O(1); no new scan, row clone or geometry conversion.
    let molecule = &prepared.molecule;
    let aromatic_bonds = &prepared.aromatic_bonds;
    let wedge_bonds = &prepared.wedge_bonds;
    let selected = &prepared.selected;
    let stereo_context = MolfileStereoContext {
        crossed: cosmolkit_core::CrossedBondContext::new(
            &molecule.topology,
            molecule.valence.as_deref().expect("prepared valence"),
            &prepared.ring_info,
            true,
        )?,
        conformer: writer_conformer(molecule, selected),
    };
    validate_v3000_writer_subset(molecule, params.include_stereo, stereo_context)?;
    let chiral_flag = molfile_chiral_flag(molecule)?;
    let generated_sgroups = v3000_generated_zbo_sgroups(
        molecule,
        molecule.valence.as_deref().expect("prepared valence"),
    );

    let mut out = PropertyText::new();
    out.extend_bytes(
        (molecule
            .properties()
            .name()
            .map_or(&b""[..], PropertyText::as_bytes))
        .as_ref(),
    );
    out.push_byte(b'\n');
    out.extend_bytes((&molfile_info_line(molecule, selected)?).as_ref());
    out.push_byte(b'\n');
    out.extend_bytes(
        (molecule
            .prop("_MolFileComments")
            .map(crate::sdf::model_string_property)
            .transpose()?
            .unwrap_or_default())
        .as_ref(),
    );
    out.push_byte(b'\n');
    out.extend_bytes(("  0  0  0  0  0  0  0  0  0  0999 V3000\n").as_ref());
    out.extend_bytes(("M  V30 BEGIN CTAB\n").as_ref());
    out.extend_bytes(
        (&format!(
            "M  V30 COUNTS {} {} {} 0 {}\n",
            molecule.num_atoms(),
            molecule.num_bonds(),
            molecule.substance_groups().len() + generated_sgroups.len(),
            chiral_flag
        ))
            .as_ref(),
    );
    let v3k_parity_flags: Vec<u32> = if selected.is_3d {
        if let Some(ref coords_3d) = selected.coords {
            let valence = molblock_valence_assignment(molecule)?;
            molecule
                .atoms()
                .iter()
                .map(|atom| get_atom_parity_flag(molecule, atom, coords_3d, &valence))
                .collect()
        } else {
            vec![0u32; molecule.num_atoms()]
        }
    } else {
        vec![0u32; molecule.num_atoms()]
    };
    out.extend_bytes(("M  V30 BEGIN ATOM\n").as_ref());
    for atom in molecule.atoms() {
        let coord = selected
            .coords
            .as_ref()
            .and_then(|coords| coords.get(atom.id().index()).copied())
            .unwrap_or([0.0, 0.0, 0.0]);
        out.extend_bytes(
            (&v3000_atom_line(molecule, atom, coord, params.precision, &v3k_parity_flags)?)
                .as_ref(),
        );
        out.push_byte(b'\n');
    }
    out.extend_bytes(("M  V30 END ATOM\n").as_ref());
    if molecule.num_bonds() != 0 {
        out.extend_bytes(("M  V30 BEGIN BOND\n").as_ref());
        for bond in molecule.bonds() {
            out.extend_bytes(
                (&v3000_bond_line(
                    molecule,
                    bond,
                    params.include_stereo,
                    aromatic_bonds,
                    wedge_bonds,
                    selected.coords.as_deref(),
                    stereo_context,
                )?)
                    .as_ref(),
            );
            out.push_byte(b'\n');
        }
        out.extend_bytes(("M  V30 END BOND\n").as_ref());
    }
    append_v3000_sgroup_lines(&mut out, molecule, &generated_sgroups)?;
    // RDKit❗✔️:   appendEnhancedStereoGroups(res, tmol, wedgeBonds);
    // getV3000CTAB emits existing groups independently of includeStereo.
    // This calls the same group formatter once, with no additional scan/copy.
    // Behavioral validation is proposed for independent p1 review.
    append_v3000_collection_lines(&mut out, molecule, &prepared.wedge_bonds)?;
    out.extend_bytes(("M  V30 END CTAB\n").as_ref());
    out.extend_bytes(("M  END\n").as_ref());
    Ok(out)
}

fn mol_to_v2000_block_with_params(
    molecule: &MolWriteContext<'_>,
    selection: CoordinateSelection,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: std::string outputMolToMolBlock(const RWMol &tmol, int confId,
    // RDKit❗❌:                                 MolFileFormat whichFormat,
    // RDKit❗❌:                                 unsigned int precision,
    // RDKit❗❌:                                 const boost::dynamic_bitset<> &aromaticBonds) {
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   unsigned int nAtoms, nBonds, nLists, chiralFlag, nsText, nRxnComponents;
    // RDKit❗❌:   unsigned int nReactants, nProducts, nIntermediates;
    // RDKit❗❌:   nAtoms = tmol.getNumAtoms();
    // RDKit❗❌:   nBonds = tmol.getNumBonds();
    // RDKit❗❌:   nLists = 0;
    // RDKit❗❌:
    // RDKit❗❌:   const auto &sgroups = getSubstanceGroups(tmol);
    // RDKit❗❌:   unsigned int nSGroups = sgroups.size();
    // RDKit❗❌:
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V2000 &&
    // RDKit❗❌:       (nAtoms > 999 || nBonds > 999 || nSGroups > 999)) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "V2000 format does not support more than 999 atoms, bonds or SGroups.");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   chiralFlag = 0;
    // RDKit❗❌:   nsText = 0;
    // RDKit❗❌:   nRxnComponents = 0;
    // RDKit❗❌:   nReactants = 0;
    // RDKit❗❌:   nProducts = 0;
    // RDKit❗❌:   nIntermediates = 0;
    // RDKit❗❌:
    // RDKit❗❌:   tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit❗❌:
    // RDKit❗❌:   const Conformer *conf;
    // RDKit❗❌:   if (confId < 0 && tmol.getNumConformers() == 0) {
    // RDKit❗❌:     conf = nullptr;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     conf = &(tmol.getConformer(confId));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool coordMagnitudeTooLargeForV2K = false;
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     for (auto &pos : conf->getPositions()) {
    // RDKit❗❌:       if ((pos.x >= MAX_V2000_COORD || pos.x <= MIN_V2000_COORD) ||
    // RDKit❗❌:           (pos.y >= MAX_V2000_COORD || pos.y <= MIN_V2000_COORD) ||
    // RDKit❗❌:           (pos.z >= MAX_V2000_COORD || pos.z <= MIN_V2000_COORD)) {
    // RDKit❗❌:         coordMagnitudeTooLargeForV2K = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V2000 && coordMagnitudeTooLargeForV2K) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "V2000 format does not support atom positions <= " +
    // RDKit❗❌:         std::to_string((int)MIN_V2000_COORD) +
    // RDKit❗❌:         " or >= " + std::to_string((int)MAX_V2000_COORD));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string text;
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:
    // RDKit❗❌:   // info
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     std::stringstream ss;
    // RDKit❗❌:     ss << "  " << std::setw(8) << "RDKit";
    // RDKit❗❌:     ss << std::setw(10) << "";
    // RDKit❗❌:     if (conf) {
    // RDKit❗❌:       if (conf->is3D()) {
    // RDKit❗❌:         ss << "3D";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         ss << common_properties::TWOD;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     res += ss.str();
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:   // comments
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:
    // RDKit❗❌:   bool hasDative = false;
    // RDKit❗❌:   for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:     if (bond->getBondType() == Bond::DATIVE) {
    // RDKit❗❌:       hasDative = true;
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool isV3000 = false;
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V3000) {
    // RDKit❗❌:     isV3000 = true;
    // RDKit❗❌:   } else if (whichFormat == MolFileFormat::unspecified &&
    // RDKit❗❌:              (coordMagnitudeTooLargeForV2K || hasDative || nAtoms > 999 ||
    // RDKit❗❌:               nBonds > 999 || nSGroups > 999 ||
    // RDKit❗❌:               !tmol.getStereoGroups().empty())) {
    // RDKit❗❌:     isV3000 = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // the counts line:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   if (isV3000) {
    // RDKit❗❌:     // All counts in the V3000 info line should be 0
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << "999 V3000\n";
    // RDKit❗❌:   } else {
    // RDKit❗❌:     ss << std::setw(3) << nAtoms;
    // RDKit❗❌:     ss << std::setw(3) << nBonds;
    // RDKit❗❌:     ss << std::setw(3) << nLists;
    // RDKit❗❌:     ss << std::setw(3) << nSGroups;
    // RDKit❗❌:     ss << std::setw(3) << chiralFlag;
    // RDKit❗❌:     ss << std::setw(3) << nsText;
    // RDKit❗❌:     ss << std::setw(3) << nRxnComponents;
    // RDKit❗❌:     ss << std::setw(3) << nReactants;
    // RDKit❗❌:     ss << std::setw(3) << nProducts;
    // RDKit❗❌:     ss << std::setw(3) << nIntermediates;
    // RDKit❗❌:     ss << "999 V2000\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   res += ss.str();
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> queryListAtoms(tmol.getNumAtoms());
    // RDKit❗❌:   if (!isV3000) {
    // RDKit❗❌:     // V2000 output.
    // RDKit❗❌:     for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
    // RDKit❗❌:          atomIt != tmol.endAtoms(); ++atomIt) {
    // RDKit❗❌:       res += GetMolFileAtomLine(*atomIt, conf, queryListAtoms);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // RDKit❗❌:
    // RDKit❗❌:     for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:       res += GetMolFileBondLine(bond, wedgeBonds, conf,
    // RDKit❗❌:                                 aromaticBonds[bond->getIdx()]);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     res += GetMolFileChargeInfo(tmol);
    // RDKit❗❌:     res += GetMolFileRGroupInfo(tmol);
    // RDKit❗❌:     res += GetMolFileQueryInfo(tmol, queryListAtoms);
    // RDKit❗❌:     res += GetMolFileAliasInfo(tmol);
    // RDKit❗❌:     res += GetMolFileZBOInfo(tmol);
    // RDKit❗❌:
    // RDKit❗❌:     res += GetMolFilePXAInfo(tmol);
    // RDKit❗❌:     res += GetMolFileSGroupInfo(tmol);
    // RDKit❗❌:
    // RDKit❗❌:     // FIX: R-group logic, SGroups and 3D features etc.
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // V3000 output.
    // RDKit❗❌:     res +=
    // RDKit❗❌:         FileParserUtils::getV3000CTAB(tmol, aromaticBonds, confId, precision);
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "M  END\n";
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    let prepared = prepare_mol_for_writing(molecule, params, selection)?;
    if should_auto_upgrade_to_v3000(&prepared.molecule, &prepared.selected) {
        return render_v3000(&prepared, params);
    }
    render_v2000(&prepared, params)
}

fn render_v2000(
    prepared: &PreparedMol<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    // RDKit❗❌: void GetMolFileAtomProperties(const Atom *atom, const Conformer *conf,
    // RDKit❗❌:                               int &totValence, int &atomMapNumber,
    // RDKit❗❌:                               unsigned int &parityFlag, double &x, double &y,
    // RDKit❗❌:                               double &z) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   totValence = 0;
    // RDKit❗❌:   atomMapNumber = 0;
    // RDKit❗❌:   parityFlag = 0;
    // RDKit❗❌:   x = y = z = 0.0;
    // RDKit❗❌:
    // RDKit❗❌:   if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                               atomMapNumber)) {
    // RDKit❗❌:     // XXX FIX ME->should we fail here? previously we would not assign
    // RDKit❗❌:     // the atomMapNumber if it didn't exist which could result in garbage
    // RDKit❗❌:     //  values.
    // RDKit❗❌:     atomMapNumber = 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:     x = pos.x;
    // RDKit❗❌:     y = pos.y;
    // RDKit❗❌:     z = pos.z;
    // RDKit❗❌:     if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:         atom->getChiralTag() != Atom::CHI_OTHER && atom->getDegree() >= 3 &&
    // RDKit❗❌:         atom->getTotalDegree() == 4) {
    // RDKit❗❌:       parityFlag = getAtomParityFlag(atom, conf);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (hasNonDefaultValence(atom)) {
    // RDKit❗❌:     if (atom->getTotalDegree() == 0) {
    // RDKit❗❌:       // Specify zero valence for elements/metals without neighbors
    // RDKit❗❌:       // or hydrogens (degree 0) instead of writing them as radicals.
    // RDKit❗❌:       totValence = 15;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // write the total valence for other atoms
    // RDKit❗❌:       totValence = atom->getTotalValence() % 15;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Selection supplies both XYZ rows and the independent stored is3D bit.
    // The source always copies x/y/z; only parity is gated by is3D. The chiral
    // tag and degree guards live in the existing parity helper. Wedge selection
    // and rendering borrow this same source conformer through writer_conformer,
    // retaining its flag and full rows for the canonical CORE helper branches.
    // These guards are O(1); no new scan, row clone or geometry conversion.
    let molecule = &prepared.molecule;
    let aromatic_bonds = &prepared.aromatic_bonds;
    let wedge_bonds = &prepared.wedge_bonds;
    let selected = &prepared.selected;
    let stereo_context = MolfileStereoContext {
        crossed: cosmolkit_core::CrossedBondContext::new(
            &molecule.topology,
            molecule.valence.as_deref().expect("prepared valence"),
            &prepared.ring_info,
            true,
        )?,
        conformer: writer_conformer(molecule, selected),
    };
    validate_v2000_writer_subset(molecule, params.include_stereo, stereo_context)?;
    validate_v2000_coordinate_range(selected.coords.as_deref())?;
    let chiral_flag = molfile_chiral_flag(molecule)?;

    let mut out = PropertyText::new();
    out.extend_bytes(
        (molecule
            .properties()
            .name()
            .map_or(&b""[..], PropertyText::as_bytes))
        .as_ref(),
    );
    out.push_byte(b'\n');
    out.extend_bytes((&molfile_info_line(molecule, selected)?).as_ref());
    out.push_byte(b'\n');
    out.extend_bytes(
        (molecule
            .prop("_MolFileComments")
            .map(crate::sdf::model_string_property)
            .transpose()?
            .unwrap_or_default())
        .as_ref(),
    );
    out.push_byte(b'\n');
    out.extend_bytes(
        (&format!(
            "{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}999 V2000\n",
            molecule.num_atoms(),
            molecule.num_bonds(),
            0,
            molecule.substance_groups().len(),
            chiral_flag,
            0,
            0,
            0,
            0,
            0
        ))
            .as_ref(),
    );

    let parity_flags: Vec<u32> = if selected.is_3d {
        if let Some(ref coords_3d) = selected.coords {
            let valence = molblock_valence_assignment(molecule)?;
            molecule
                .atoms()
                .iter()
                .map(|atom| get_atom_parity_flag(molecule, atom, coords_3d, &valence))
                .collect()
        } else {
            vec![0u32; molecule.num_atoms()]
        }
    } else {
        vec![0u32; molecule.num_atoms()]
    };
    for atom in molecule.atoms() {
        let coord = selected
            .coords
            .as_ref()
            .and_then(|coords| coords.get(atom.id().index()).copied())
            .unwrap_or([0.0, 0.0, 0.0]);
        out.extend_bytes((&v2000_atom_line(atom, coord, molecule, &parity_flags)?).as_ref());
        out.push_byte(b'\n');
    }
    for bond in molecule.bonds() {
        out.extend_bytes(
            (&v2000_bond_line(
                molecule,
                bond,
                params.include_stereo,
                aromatic_bonds,
                wedge_bonds,
                selected.coords.as_deref(),
                stereo_context,
            )?)
                .as_ref(),
        );
        out.push_byte(b'\n');
    }
    append_v2000_property_lines(&mut out, molecule)?;
    append_v2000_rgroup_lines(&mut out, molecule)?;
    append_v2000_value_lines(&mut out, molecule)?;
    append_v2000_alias_lines(&mut out, molecule)?;
    append_v2000_zbo_lines(
        &mut out,
        molecule,
        molecule.valence.as_deref().expect("prepared valence"),
    );
    append_v2000_pxa_lines(&mut out, molecule)?;
    append_v2000_sgroup_lines(&mut out, molecule)?;
    out.extend_bytes(("M  END\n").as_ref());
    Ok(out)
}

/// Validate that the molecule's query atoms/bonds can be written in
/// the V2000/V3000 format. Rejects only unsupported recursive SMARTS,
/// RGroupLabel, and MolFileAlias queries.
fn molfile_chiral_flag(molecule: &MolWriteContext<'_>) -> Result<u32, MolWriteError> {
    Ok(molecule
        .prop("_MolFileChiralFlag")
        .map(cosmolkit_core::property_value_to_uint)
        .transpose()?
        .unwrap_or(0))
}

fn select_coordinates(
    molecule: &MolWriteContext<'_>,
    selection: CoordinateSelection,
) -> Result<SelectedCoordinates, MolWriteError> {
    // RDKit❗❌: const Conformer &ROMol::getConformer(int id) const {
    // RDKit❗❌:   // make sure we have more than one conformation
    // RDKit❗❌:   if (d_confs.size() == 0) {
    // RDKit❗❌:     throw ConformerException("No conformations available on the molecule");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (id < 0) {
    // RDKit❗❌:     return *(d_confs.front());
    // RDKit❗❌:   }
    // RDKit❗❌:   auto cid = (unsigned int)id;
    // RDKit❗❌:   for (auto conf : d_confs) {
    // RDKit❗❌:     if (conf->getId() == cid) {
    // RDKit❗❌:       return *conf;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // we did not find a conformation with the specified ID
    // RDKit❗❌:   std::string mesg = "Can't find conformation with ID: ";
    // RDKit❗❌:   mesg += id;
    // RDKit❗❌:   throw ConformerException(mesg);
    // RDKit❗❌: }
    // Choosing XYZ storage and ID preserves the selected conformer's
    // independent is_3d flag, as in the source writer's conf->is3D() checks.
    // Full XYZ rows remain available even when the generated label is 2D.
    // Canonical difference authorized by the approved Selection contract and
    // ROOT IO44 decision: do not use source dimension, position, smallest ID,
    // or insertion order. Extra materialization of 2D as 3D remains a known
    // allocation cost, so the performance axis stays ❌.
    let none = || SelectedCoordinates {
        coords: None,
        is_3d: false,
        selection: None,
    };
    let two_d = molecule.coordinates.conformers_2d.len();
    let three_d = molecule.coordinates.conformers_3d.len();
    let selection = match selection {
        CoordinateSelection::Disabled => return Ok(none()),
        CoordinateSelection::Auto => match (two_d, three_d) {
            (0, 0) => return Ok(none()),
            (1, 0) => CoordinateSelection::TwoD(Some(molecule.coordinates.conformers_2d[0].id())),
            (0, 1) => CoordinateSelection::ThreeD(Some(molecule.coordinates.conformers_3d[0].id())),
            _ => return Err(MolWriteError::AmbiguousCoordinates { two_d, three_d }),
        },
        CoordinateSelection::TwoD(None) => {
            if two_d == 1 {
                CoordinateSelection::TwoD(Some(molecule.coordinates.conformers_2d[0].id()))
            } else if two_d == 0 {
                return Err(MolWriteError::MissingCoordinate {
                    dimension: CoordinateDimension::TwoD,
                    id: 0,
                });
            } else {
                return Err(MolWriteError::AmbiguousCoordinates { two_d, three_d: 0 });
            }
        }
        CoordinateSelection::ThreeD(None) => {
            if three_d == 1 {
                CoordinateSelection::ThreeD(Some(molecule.coordinates.conformers_3d[0].id()))
            } else if three_d == 0 {
                return Err(MolWriteError::MissingCoordinate {
                    dimension: CoordinateDimension::ThreeD,
                    id: 0,
                });
            } else {
                return Err(MolWriteError::AmbiguousCoordinates { two_d: 0, three_d });
            }
        }
        explicit => explicit,
    };
    match selection {
        CoordinateSelection::TwoD(Some(id)) => {
            let conformer = molecule
                .coordinates
                .conformers_2d
                .iter()
                .find(|c| c.id() == id)
                .ok_or(MolWriteError::MissingCoordinate {
                    dimension: CoordinateDimension::TwoD,
                    id,
                })?;
            Ok(SelectedCoordinates {
                coords: Some(
                    conformer
                        .coordinates()
                        .iter()
                        .map(|p| [p[0], p[1], 0.0])
                        .collect(),
                ),
                is_3d: false,
                selection: Some(MolCoordinateSelection::TwoD { id }),
            })
        }
        CoordinateSelection::ThreeD(Some(id)) => {
            let conformer = molecule
                .coordinates
                .conformers_3d
                .iter()
                .find(|c| c.id() == id)
                .ok_or(MolWriteError::MissingCoordinate {
                    dimension: CoordinateDimension::ThreeD,
                    id,
                })?;
            Ok(SelectedCoordinates {
                coords: Some(conformer.coordinates().to_vec()),
                is_3d: conformer.is_3d(),
                selection: Some(MolCoordinateSelection::ThreeD { id }),
            })
        }
        _ => unreachable!("all implicit selections resolved"),
    }
}

fn validate_v2000_writer_subset(
    molecule: &MolWriteContext<'_>,
    include_stereo: bool,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<(), MolWriteError> {
    if molecule.num_atoms() > 999
        || molecule.num_bonds() > 999
        || molecule.substance_groups().len() > 999
    {
        return Err(MolWriteError::Value(
            "V2000 format does not support more than 999 atoms, bonds or SGroups".to_string(),
        ));
    }
    if include_stereo && !molecule.stereo_groups().is_empty() {
        return Err(MolWriteError::UnsupportedSubset(
            "MolBlock enhanced stereo writing is not ported",
        ));
    }
    for bond in molecule.bonds() {
        if include_stereo {
            let empty_wedge_bonds = WedgeAssignments::default();
            v2000_bond_stereo_code(molecule, bond, &empty_wedge_bonds, None, stereo_context)?;
        }
        v2000_bond_type_code(bond)?;
    }
    Ok(())
}

fn validate_v3000_writer_subset(
    molecule: &MolWriteContext<'_>,
    include_stereo: bool,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<(), MolWriteError> {
    if include_stereo {
        validate_v3000_stereo_groups(molecule)?;
    }
    for bond in molecule.bonds() {
        if include_stereo {
            let empty_wedge_bonds = WedgeAssignments::default();
            v3000_bond_cfg_code(molecule, bond, &empty_wedge_bonds, None, stereo_context)?;
        }
        v3000_bond_type_code(bond)?;
    }
    Ok(())
}

fn validate_v3000_stereo_groups(molecule: &MolWriteContext<'_>) -> Result<(), MolWriteError> {
    // Atropisomer (bond-based) stereo groups are now supported via atom
    // collection from the bond endpoints. No separate validation needed.
    // display; COSMolKit does not model wedge bonds for stereo groups yet.
    let _ = molecule;
    Ok(())
}

fn validate_v2000_coordinate_range(coords: Option<&[[f64; 3]]>) -> Result<(), MolWriteError> {
    let Some(coords) = coords else {
        return Ok(());
    };
    for coord in coords {
        for value in coord {
            if *value >= MAX_V2000_COORD || *value <= MIN_V2000_COORD {
                return Err(MolWriteError::Value(
                    "V2000 atom positions must be > -100000 and < 1000000".to_string(),
                ));
            }
        }
    }
    Ok(())
}

fn molfile_info_line(
    molecule: &MolWriteContext<'_>,
    selected: &SelectedCoordinates,
) -> Result<PropertyText, MolWriteError> {
    // RDKit❗❌: std::string outputMolToMolBlock(const RWMol &tmol, int confId,
    // RDKit❗❌:                                 MolFileFormat whichFormat,
    // RDKit❗❌:                                 unsigned int precision,
    // RDKit❗❌:                                 const boost::dynamic_bitset<> &aromaticBonds) {
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   unsigned int nAtoms, nBonds, nLists, chiralFlag, nsText, nRxnComponents;
    // RDKit❗❌:   unsigned int nReactants, nProducts, nIntermediates;
    // RDKit❗❌:   nAtoms = tmol.getNumAtoms();
    // RDKit❗❌:   nBonds = tmol.getNumBonds();
    // RDKit❗❌:   nLists = 0;
    // RDKit❗❌:
    // RDKit❗❌:   const auto &sgroups = getSubstanceGroups(tmol);
    // RDKit❗❌:   unsigned int nSGroups = sgroups.size();
    // RDKit❗❌:
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V2000 &&
    // RDKit❗❌:       (nAtoms > 999 || nBonds > 999 || nSGroups > 999)) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "V2000 format does not support more than 999 atoms, bonds or SGroups.");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   chiralFlag = 0;
    // RDKit❗❌:   nsText = 0;
    // RDKit❗❌:   nRxnComponents = 0;
    // RDKit❗❌:   nReactants = 0;
    // RDKit❗❌:   nProducts = 0;
    // RDKit❗❌:   nIntermediates = 0;
    // RDKit❗❌:
    // RDKit❗❌:   tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit❗❌:
    // RDKit❗❌:   const Conformer *conf;
    // RDKit❗❌:   if (confId < 0 && tmol.getNumConformers() == 0) {
    // RDKit❗❌:     conf = nullptr;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     conf = &(tmol.getConformer(confId));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool coordMagnitudeTooLargeForV2K = false;
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     for (auto &pos : conf->getPositions()) {
    // RDKit❗❌:       if ((pos.x >= MAX_V2000_COORD || pos.x <= MIN_V2000_COORD) ||
    // RDKit❗❌:           (pos.y >= MAX_V2000_COORD || pos.y <= MIN_V2000_COORD) ||
    // RDKit❗❌:           (pos.z >= MAX_V2000_COORD || pos.z <= MIN_V2000_COORD)) {
    // RDKit❗❌:         coordMagnitudeTooLargeForV2K = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V2000 && coordMagnitudeTooLargeForV2K) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "V2000 format does not support atom positions <= " +
    // RDKit❗❌:         std::to_string((int)MIN_V2000_COORD) +
    // RDKit❗❌:         " or >= " + std::to_string((int)MAX_V2000_COORD));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string text;
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:
    // RDKit❗❌:   // info
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     std::stringstream ss;
    // RDKit❗❌:     ss << "  " << std::setw(8) << "RDKit";
    // RDKit❗❌:     ss << std::setw(10) << "";
    // RDKit❗❌:     if (conf) {
    // RDKit❗❌:       if (conf->is3D()) {
    // RDKit❗❌:         ss << "3D";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         ss << common_properties::TWOD;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     res += ss.str();
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:   // comments
    // RDKit❗❌:   if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) {
    // RDKit❗❌:     res += text;
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "\n";
    // RDKit❗❌:
    // RDKit❗❌:   bool hasDative = false;
    // RDKit❗❌:   for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:     if (bond->getBondType() == Bond::DATIVE) {
    // RDKit❗❌:       hasDative = true;
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool isV3000 = false;
    // RDKit❗❌:   if (whichFormat == MolFileFormat::V3000) {
    // RDKit❗❌:     isV3000 = true;
    // RDKit❗❌:   } else if (whichFormat == MolFileFormat::unspecified &&
    // RDKit❗❌:              (coordMagnitudeTooLargeForV2K || hasDative || nAtoms > 999 ||
    // RDKit❗❌:               nBonds > 999 || nSGroups > 999 ||
    // RDKit❗❌:               !tmol.getStereoGroups().empty())) {
    // RDKit❗❌:     isV3000 = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // the counts line:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   if (isV3000) {
    // RDKit❗❌:     // All counts in the V3000 info line should be 0
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << std::setw(3) << 0;
    // RDKit❗❌:     ss << "999 V3000\n";
    // RDKit❗❌:   } else {
    // RDKit❗❌:     ss << std::setw(3) << nAtoms;
    // RDKit❗❌:     ss << std::setw(3) << nBonds;
    // RDKit❗❌:     ss << std::setw(3) << nLists;
    // RDKit❗❌:     ss << std::setw(3) << nSGroups;
    // RDKit❗❌:     ss << std::setw(3) << chiralFlag;
    // RDKit❗❌:     ss << std::setw(3) << nsText;
    // RDKit❗❌:     ss << std::setw(3) << nRxnComponents;
    // RDKit❗❌:     ss << std::setw(3) << nReactants;
    // RDKit❗❌:     ss << std::setw(3) << nProducts;
    // RDKit❗❌:     ss << std::setw(3) << nIntermediates;
    // RDKit❗❌:     ss << "999 V2000\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   res += ss.str();
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> queryListAtoms(tmol.getNumAtoms());
    // RDKit❗❌:   if (!isV3000) {
    // RDKit❗❌:     // V2000 output.
    // RDKit❗❌:     for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
    // RDKit❗❌:          atomIt != tmol.endAtoms(); ++atomIt) {
    // RDKit❗❌:       res += GetMolFileAtomLine(*atomIt, conf, queryListAtoms);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // RDKit❗❌:
    // RDKit❗❌:     for (const auto bond : tmol.bonds()) {
    // RDKit❗❌:       res += GetMolFileBondLine(bond, wedgeBonds, conf,
    // RDKit❗❌:                                 aromaticBonds[bond->getIdx()]);
    // RDKit❗❌:       res += "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     res += GetMolFileChargeInfo(tmol);
    // RDKit❗❌:     res += GetMolFileRGroupInfo(tmol);
    // RDKit❗❌:     res += GetMolFileQueryInfo(tmol, queryListAtoms);
    // RDKit❗❌:     res += GetMolFileAliasInfo(tmol);
    // RDKit❗❌:     res += GetMolFileZBOInfo(tmol);
    // RDKit❗❌:
    // RDKit❗❌:     res += GetMolFilePXAInfo(tmol);
    // RDKit❗❌:     res += GetMolFileSGroupInfo(tmol);
    // RDKit❗❌:
    // RDKit❗❌:     // FIX: R-group logic, SGroups and 3D features etc.
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // V3000 output.
    // RDKit❗❌:     res +=
    // RDKit❗❌:         FileParserUtils::getV3000CTAB(tmol, aromaticBonds, confId, precision);
    // RDKit❗❌:   }
    // RDKit❗❌:   res += "M  END\n";
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // The generated header reads the actual flag of the selected carrier.
    // XYZ storage dimension is independent of this bit, including nonzero Z.
    // The existing CK explicit-selection metadata policy and program spelling
    // remain distinct from source implicit selection and literal RDKit text.
    // Cost: one O(1) flag branch and one header string allocation; selected XYZ
    // copying elsewhere remains an O(A) cost over source borrowed storage.
    let label = selected
        .coords
        .as_ref()
        .map(|_| if selected.is_3d { "3D" } else { "2D" });
    if let Some(label) = label {
        let mut line = format!("  {:>8}{:>10}", "COSMolKit", "");
        line.push_str(label);
        return Ok(line.into());
    }
    if let Some(info) = molecule.prop("_MolFileInfo") {
        return Ok(crate::sdf::model_string_property(info)?);
    }
    if let Some(info) = molecule.prop("_MolFileInfoLine") {
        return Ok(crate::sdf::model_string_property(info)?);
    }
    let mut line = format!("  {:>8}{:>10}", "COSMolKit", "");
    Ok(line.into())
}

fn v2000_atom_line(
    atom: &Atom,
    coord: [f64; 3],
    molecule: &MolWriteContext<'_>,
    parity_flags: &[u32],
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFileAtomLine(const Atom *atom, const Conformer *conf,
    // RDKit❗❌:                                      boost::dynamic_bitset<> &queryListAtoms) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   int totValence, atomMapNumber;
    // RDKit❗❌:   unsigned int parityFlag;
    // RDKit❗❌:   double x, y, z;
    // RDKit❗❌:   GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
    // RDKit❗❌:                            y, z);
    // RDKit❗❌:
    // RDKit❗❌:   if ((x >= MAX_V2000_COORD || x <= MIN_V2000_COORD) ||
    // RDKit❗❌:       (y >= MAX_V2000_COORD || y <= MIN_V2000_COORD) ||
    // RDKit❗❌:       (z >= MAX_V2000_COORD || z <= MIN_V2000_COORD)) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "MolFile coordinates must be in (-100000, 1000000)");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   int massDiff, chg, stereoCare, hCount, rxnComponentType, rxnComponentNumber,
    // RDKit❗❌:       inversionFlag, exactChangeFlag;
    // RDKit❗❌:   massDiff = 0;
    // RDKit❗❌:   chg = 0;
    // RDKit❗❌:   stereoCare = 0;
    // RDKit❗❌:   hCount = 0;
    // RDKit❗❌:   rxnComponentType = 0;
    // RDKit❗❌:   rxnComponentNumber = 0;
    // RDKit❗❌:   inversionFlag = 0;
    // RDKit❗❌:   exactChangeFlag = 0;
    // RDKit❗❌:
    // RDKit❗❌:   atom->getPropIfPresent(common_properties::molRxnRole, rxnComponentType);
    // RDKit❗❌:   atom->getPropIfPresent(common_properties::molRxnComponent,
    // RDKit❗❌:                          rxnComponentNumber);
    // RDKit❗❌:
    // RDKit❗❌:   std::string symbol = AtomGetMolFileSymbol(atom, true, queryListAtoms);
    // RDKit❗❌:   // it feels ugly to use snprintf instead of boost::format, but at least of the
    // RDKit❗❌:   // time of this writing (with boost 1.55), the snprintf version runs in 20% of
    // RDKit❗❌:   // the time.
    // RDKit❗❌:   char dest[128];
    // RDKit❗❌: #ifndef _MSC_VER
    // RDKit❗❌:   snprintf(dest, 128,
    // RDKit❗❌:            "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
    // RDKit❗❌:            z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
    // RDKit❗❌:            totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
    // RDKit❗❌:            inversionFlag, exactChangeFlag);
    // RDKit❗❌: #else
    // RDKit❗❌:   // ok, technically we should be being more careful about this, but given that
    // RDKit❗❌:   // the format string makes it impossible for this to overflow, I think we're
    // RDKit❗❌:   // safe. I just used the snprintf above to prevent linters from complaining
    // RDKit❗❌:   // about use of sprintf
    // RDKit❗❌:   sprintf_s(dest, 128,
    // RDKit❗❌:             "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
    // RDKit❗❌:             z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
    // RDKit❗❌:             totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
    // RDKit❗❌:             inversionFlag, exactChangeFlag);
    // RDKit❗❌:
    // RDKit❗❌: #endif
    // RDKit❗❌:   res += dest;
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    let symbol = v2000_atom_symbol(atom, true, molecule.query_atom(atom.id()))?;
    let atom_idx = atom.id().index();
    let parity_flag = parity_flags.get(atom_idx).copied().unwrap_or(0);
    let tot_valence = molfile_total_valence_field(molecule, atom)?;
    crate::sdf::format_v2000_atom_line(
        atom,
        coord,
        symbol.as_bytes(),
        parity_flag as i32,
        0,
        0,
        tot_valence as i32,
        0,
        0,
    )
    .map_err(MolWriteError::Detached)
}

fn get_atom_parity_flag(
    molecule: &MolWriteContext<'_>,
    atom: &Atom,
    coords_3d: &[[f64; 3]],
    valence: &ValenceAssignment,
) -> u32 {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf) {
    // RDKit❗❌:   PRECONDITION(atom, "bad atom");
    // RDKit❗❌:   PRECONDITION(conf, "bad conformer");
    // RDKit❗❌:   if (!conf->is3D() ||
    // RDKit❗❌:       !(atom->getDegree() >= 3 && atom->getTotalDegree() == 4)) {
    // RDKit❗❌:     return 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   const ROMol &mol = atom->getOwningMol();
    // RDKit❗❌:   RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:   std::vector<std::pair<unsigned int, RDGeom::Point3D>> vs;
    // RDKit❗❌:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit❗❌:   boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit❗❌:   while (nbrIdx != endNbrs) {
    // RDKit❗❌:     const Atom *at = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:     unsigned int idx = at->getIdx();
    // RDKit❗❌:     RDGeom::Point3D v = conf->getAtomPos(idx);
    // RDKit❗❌:     v -= pos;
    // RDKit❗❌:     if (at->getAtomicNum() == 1) {
    // RDKit❗❌:       idx += mol.getNumAtoms();
    // RDKit❗❌:     }
    // RDKit❗❌:     vs.emplace_back(idx, v);
    // RDKit❗❌:     ++nbrIdx;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::sort(vs.begin(), vs.end(), Rankers::pairLess);
    // RDKit❗❌:   double vol;
    // RDKit❗❌:   if (vs.size() == 4) {
    // RDKit❗❌:     vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     vol = -vs[0].second.crossProduct(vs[1].second).dotProduct(vs[2].second);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (vol < 0) {
    // RDKit❗❌:     return 2;
    // RDKit❗❌:   } else if (vol > 0) {
    // RDKit❗❌:     return 1;
    // RDKit❗❌:   }
    // RDKit❗❌:   return 0;
    // RDKit❗❌: }

    let atom_idx = atom.id().index();
    if matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) {
        return 0;
    }
    let adjacency = &molecule.topology_block().adjacency;
    let neighbors = adjacency.neighbors_of(atom_idx);
    if !(neighbors.len() >= 3 && molfile_total_degree(molecule, atom, valence) == 4) {
        return 0;
    }
    // Compute vectors from the center atom to each neighbor.
    let pos = coords_3d[atom_idx];
    let mut vs: Vec<(usize, [f64; 3])> = neighbors
        .iter()
        .map(|nbr| {
            let npos = coords_3d[nbr.atom_index];
            let idx = nbr.atom_index;
            let v = [npos[0] - pos[0], npos[1] - pos[1], npos[2] - pos[2]];
            let sort_key = if molecule.atoms()[idx].atomic_number() == 1 {
                idx + molecule.num_atoms()
            } else {
                idx
            };
            (sort_key, v)
        })
        .collect();
    vs.sort_by(|a, b| a.0.cmp(&b.0));
    let cross = [
        vs[0].1[1] * vs[1].1[2] - vs[0].1[2] * vs[1].1[1],
        vs[0].1[2] * vs[1].1[0] - vs[0].1[0] * vs[1].1[2],
        vs[0].1[0] * vs[1].1[1] - vs[0].1[1] * vs[1].1[0],
    ];
    let vol = if vs.len() == 4 {
        cross[0] * vs[3].1[0] + cross[1] * vs[3].1[1] + cross[2] * vs[3].1[2]
    } else {
        -(cross[0] * vs[2].1[0] + cross[1] * vs[2].1[1] + cross[2] * vs[2].1[2])
    };
    // The source tests the exact sign, including signed zero and NaN falling
    // through to zero; no geometry-wide tolerance is involved.
    if vol < 0.0 {
        2
    } else if vol > 0.0 {
        1
    } else {
        0
    }
}

fn has_non_default_valence(
    molecule: &MolWriteContext<'_>,
    atom: &Atom,
    valence: &ValenceAssignment,
) -> Result<bool, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: bool hasNonDefaultValence(const Atom *atom) {
    // RDKit❗❌:   if (atom->getNumRadicalElectrons() != 0) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   // for queries and atoms which don't have computed properties, the answer is
    // RDKit❗❌:   // always no:
    // RDKit❗❌:   if (atom->hasQuery() || atom->needsUpdatePropertyCache()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (atom->getAtomicNum() == 1 ||
    // RDKit❗❌:       SmilesWrite ::inOrganicSubset(atom->getAtomicNum())) {
    // RDKit❗❌:     // for the ones we "know", we may have to specify the valence if it's
    // RDKit❗❌:     // not the default value
    // RDKit❗❌:     auto effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
    // RDKit❗❌:     return atom->getNoImplicit() &&
    // RDKit❗❌:            (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) !=
    // RDKit❗❌:             PeriodicTable::getTable()->getDefaultValence(effAtomicNum));
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }

    if atom.radical_electrons() != 0 {
        return Ok(true);
    }

    if molecule.query_atom(atom.id()).is_some() {
        return Ok(false);
    }

    if atom.atomic_number() == 1 || cosmolkit_core::is_rdkit_organic_subset(atom.atomic_number()) {
        if !atom.no_implicit() {
            return Ok(false);
        }
        let effective_atomic_num =
            u32::from(atom.atomic_number()).wrapping_sub(i32::from(atom.formal_charge()) as u32);
        let effective_atomic_num = u8::try_from(effective_atomic_num)
            .map_err(|_| MolWriteError::Value("Atomic number not found".to_string()))?;
        let default_valence = cosmolkit_core::rdkit_default_valence(effective_atomic_num)
            .map_err(|_| MolWriteError::Value("Atomic number not found".to_string()))?;
        let explicit_valence = valence.explicit_valence[atom.id().index()];
        return Ok(explicit_valence != default_valence);
    }
    let _ = molecule;
    Ok(true)
}

fn molblock_valence_assignment<'a>(
    molecule: &'a MolWriteContext<'_>,
) -> Result<Cow<'a, ValenceAssignment>, MolWriteError> {
    if let Some(valence) = &molecule.valence {
        return Ok(Cow::Borrowed(valence));
    }
    cosmolkit_core::assign_valence(
        &molecule.topology,
        &cosmolkit_core::ValenceParams {
            strict: false,
            ..Default::default()
        },
    )
    .map(Cow::Owned)
    .map_err(MolWriteError::Valence)
}

fn molfile_total_valence_field(
    molecule: &MolWriteContext<'_>,
    atom: &Atom,
) -> Result<u32, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: void GetMolFileAtomProperties(const Atom *atom, const Conformer *conf,
    // RDKit❗❌:                               int &totValence, int &atomMapNumber,
    // RDKit❗❌:                               unsigned int &parityFlag, double &x, double &y,
    // RDKit❗❌:                               double &z) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   totValence = 0;
    // RDKit❗❌:   atomMapNumber = 0;
    // RDKit❗❌:   parityFlag = 0;
    // RDKit❗❌:   x = y = z = 0.0;
    // RDKit❗❌:
    // RDKit❗❌:   if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                               atomMapNumber)) {
    // RDKit❗❌:     // XXX FIX ME->should we fail here? previously we would not assign
    // RDKit❗❌:     // the atomMapNumber if it didn't exist which could result in garbage
    // RDKit❗❌:     //  values.
    // RDKit❗❌:     atomMapNumber = 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:     x = pos.x;
    // RDKit❗❌:     y = pos.y;
    // RDKit❗❌:     z = pos.z;
    // RDKit❗❌:     if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:         atom->getChiralTag() != Atom::CHI_OTHER && atom->getDegree() >= 3 &&
    // RDKit❗❌:         atom->getTotalDegree() == 4) {
    // RDKit❗❌:       parityFlag = getAtomParityFlag(atom, conf);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (hasNonDefaultValence(atom)) {
    // RDKit❗❌:     if (atom->getTotalDegree() == 0) {
    // RDKit❗❌:       // Specify zero valence for elements/metals without neighbors
    // RDKit❗❌:       // or hydrogens (degree 0) instead of writing them as radicals.
    // RDKit❗❌:       totValence = 15;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // write the total valence for other atoms
    // RDKit❗❌:       totValence = atom->getTotalValence() % 15;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }

    let valence = molblock_valence_assignment(molecule)?;
    if !has_non_default_valence(molecule, atom, &valence)? {
        return Ok(0);
    }
    let total_degree = molfile_total_degree(molecule, atom, &valence);
    if total_degree == 0 {
        Ok(15)
    } else {
        let total_valence = valence.explicit_valence[atom.id().index()]
            + valence.implicit_hydrogens[atom.id().index()].max(0);
        Ok((total_valence as u32) % 15)
    }
}

fn molfile_total_hydrogens(atom: &Atom, valence: &ValenceAssignment) -> i32 {
    // BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs
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
    // END RDKIT CPP FUNCTION
    // BEGIN RDKIT CPP FUNCTION Atom::getNumImplicitHs
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
    // RDKit❗✔️:
    // END RDKIT CPP FUNCTION
    // This private writer specialization uses source includeNeighbors=false.
    // PreparedMol already computed canonical non-strict valence and propagated
    // its errors. Reuse that cache, respecting noImplicit, with O(1) indexed
    // lookup and no new hydrogen perception, tree clone, or topology scan.
    i32::from(atom.explicit_hydrogens())
        + if atom.no_implicit() {
            0
        } else {
            valence.implicit_hydrogens[atom.id().index()]
        }
}

fn molfile_total_degree(
    molecule: &MolWriteContext<'_>,
    atom: &Atom,
    valence: &ValenceAssignment,
) -> i32 {
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.id().index())
        .len() as i32
        + i32::from(atom.explicit_hydrogens())
        + valence.implicit_hydrogens[atom.id().index()].max(0)
}

fn v2000_bond_line(
    molecule: &MolWriteContext<'_>,
    bond: &Bond,
    include_stereo: bool,
    aromatic_bonds: &[usize],
    wedge_bonds: &WedgeAssignments,
    coords: Option<&[[f64; 3]]>,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFileBondLine(
    // RDKit❗❌:     const Bond *bond,
    // RDKit❗❌:     const std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds,
    // RDKit❗❌:     const Conformer *conf, bool wasAromatic) {
    // RDKit❗❌:   PRECONDITION(bond, "");
    // RDKit❗❌:
    // RDKit❗❌:   int dirCode = 0;
    // RDKit❗❌:   bool reverse = false;
    // RDKit❗❌:   RDKit::Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode,
    // RDKit❗❌:                                              reverse);
    // RDKit❗❌:   // do not cross bonds which were aromatic before kekulization
    // RDKit❗❌:   if (wasAromatic && dirCode == 3) {
    // RDKit❗❌:     dirCode = 0;
    // RDKit❗❌:   }
    // RDKit❗❌:   int symbol = BondGetMolFileSymbol(bond);
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   if (reverse) {
    // RDKit❗❌:     // switch the begin and end atoms on the bond line
    // RDKit❗❌:     ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    // RDKit❗❌:     ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    // RDKit❗❌:     ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    // RDKit❗❌:   }
    // RDKit❗❌:   ss << std::setw(3) << symbol;
    // RDKit❗❌:   ss << " " << std::setw(2) << dirCode;
    // RDKit❗❌:
    // RDKit❗❌:   if (bond->hasQuery()) {
    // RDKit❗❌:     int topol = getQueryBondTopology(bond);
    // RDKit❗❌:     if (topol) {
    // RDKit❗❌:       ss << " " << std::setw(2) << 0 << " " << std::setw(2) << topol;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return ss.str();
    // RDKit❗❌: }

    let (mut stereo_code, reverse) =
        v2000_bond_stereo_code(molecule, bond, wedge_bonds, coords, stereo_context)?;
    if aromatic_bonds.contains(&bond.id().index()) && stereo_code == 3 {
        stereo_code = 0;
    }
    let type_code = query_bond_type_code(molecule.query_bond(bond.id()))
        .filter(|code| *code != 0)
        .map_or_else(|| v2000_bond_type_code(bond), Ok)?;
    let (begin_idx, end_idx) = if reverse {
        (bond.end().index(), bond.begin().index())
    } else {
        (bond.begin().index(), bond.end().index())
    };
    let mut line = format!(
        "{:>3}{:>3}{:>3} {:>2}",
        begin_idx + 1,
        end_idx + 1,
        type_code,
        stereo_code
    );

    let topology = molecule
        .query_bond(bond.id())
        .map_or(0, |row| query_bond_topology(row.predicate()));
    if topology != 0 {
        line.push_str(&format!(" {:>2} {:>2}", 0, topology));
    }
    Ok(PropertyText::from(line))
}

fn v2000_bond_stereo_code(
    molecule: &MolWriteContext<'_>,
    bond: &Bond,
    wedge_bonds: &WedgeAssignments,
    coords: Option<&[[f64; 3]]>,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<(u32, bool), MolWriteError> {
    //
    //
    let (dir, reverse) =
        molfile_bond_stereo_info(molecule, bond, wedge_bonds, coords, stereo_context)?;
    Ok((molfile_bond_dir_code(dir), reverse))
}

fn molfile_bond_stereo_info(
    molecule: &MolWriteContext<'_>,
    bond: &Bond,
    wedge_bonds: &WedgeAssignments,
    coords: Option<&[[f64; 3]]>,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<(BondDirection, bool), MolWriteError> {
    let info = cosmolkit_core::get_molfile_bond_stereo_info(
        &stereo_context.crossed,
        wedge_bonds,
        bond.id(),
        stereo_context.conformer,
    )?;
    Ok((info.direction, info.reverse))
}

fn can_have_direction_for_molfile(bond: &Bond) -> bool {
    matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
}

fn molfile_bond_dir_code(dir: BondDirection) -> u32 {
    match dir {
        BondDirection::None => 0,
        BondDirection::BeginWedge => 1,
        BondDirection::BeginDash => 6,
        BondDirection::Unknown => 4,
        BondDirection::EitherDouble => 3,
        BondDirection::EndUpRight | BondDirection::EndDownRight => 0,
    }
}

fn v2000_bond_type_code(bond: &Bond) -> Result<u32, MolWriteError> {
    match bond.order() {
        BondOrder::Single => Ok(if bond.is_aromatic() { 4 } else { 1 }),
        BondOrder::Double => Ok(if bond.is_aromatic() { 4 } else { 2 }),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Zero => Ok(1),
        BondOrder::Dative => Ok(9),
        _ => Ok(0),
    }
}

fn v3000_atom_line(
    molecule: &MolWriteContext<'_>,
    atom: &Atom,
    coord: [f64; 3],
    precision: usize,
    parity_flags: &[u32],
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetV3000MolFileAtomLine(
    // RDKit❗❌:     const Atom *atom, const Conformer *conf,
    // RDKit❗❌:     boost::dynamic_bitset<> &queryListAtoms, unsigned int precision) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   int totValence, atomMapNumber;
    // RDKit❗❌:   unsigned int parityFlag;
    // RDKit❗❌:   double x, y, z;
    // RDKit❗❌:   GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
    // RDKit❗❌:                            y, z);
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   ss << "M  V30 " << atom->getIdx() + 1;
    // RDKit❗❌:
    // RDKit❗❌:   std::string symbol = AtomGetMolFileSymbol(atom, false, queryListAtoms);
    // RDKit❗❌:   if (!isAtomListQuery(atom) || queryListAtoms[atom->getIdx()]) {
    // RDKit❗❌:     ss << " " << symbol;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     INT_VECT vals;
    // RDKit❗❌:     getAtomListQueryVals(atom->getQuery(), vals);
    // RDKit❗❌:     if (atom->getQuery()->getNegation()) {
    // RDKit❗❌:       ss << " "
    // RDKit❗❌:          << "\"NOT";
    // RDKit❗❌:     }
    // RDKit❗❌:     ss << " [";
    // RDKit❗❌:     for (unsigned int i = 0; i < vals.size(); ++i) {
    // RDKit❗❌:       if (i != 0) {
    // RDKit❗❌:         ss << ",";
    // RDKit❗❌:       }
    // RDKit❗❌:       ss << PeriodicTable::getTable()->getElementSymbol(vals[i]);
    // RDKit❗❌:     }
    // RDKit❗❌:     ss << "]";
    // RDKit❗❌:     if (atom->getQuery()->getNegation()) {
    // RDKit❗❌:       ss << "\"";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::streamsize currentPrecision = ss.precision();
    // RDKit❗❌:   ss << std::fixed;
    // RDKit❗❌:   ss << std::setprecision(precision);
    // RDKit❗❌:   ss << " " << x << " " << y << " " << z;
    // RDKit❗❌:   ss << std::setprecision(currentPrecision);
    // RDKit❗❌:   ss << std::defaultfloat;
    // RDKit❗❌:   ss << " " << atomMapNumber;
    // RDKit❗❌:
    // RDKit❗❌:   // Extra atom properties.
    // RDKit❗❌:   int chg = atom->getFormalCharge();
    // RDKit❗❌:   int isotope = atom->getIsotope();
    // RDKit❗❌:   if (parityFlag != 0) {
    // RDKit❗❌:     ss << " CFG=" << parityFlag;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (chg != 0) {
    // RDKit❗❌:     ss << " CHG=" << chg;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (isotope != 0 && !isAtomRGroup(*atom)) {
    // RDKit❗❌:     // the documentation for V3000 CTABs says that this should contain the
    // RDKit❗❌:     // "absolute atomic weight" (whatever that means).
    // RDKit❗❌:     // Online examples seem to have integer (isotope) values and Marvin won't
    // RDKit❗❌:     // even read something that has a float.
    // RDKit❗❌:     // We'll go with the int.
    // RDKit❗❌:     int mass = static_cast<int>(std::round(atom->getMass()));
    // RDKit❗❌:     // dummies may have an isotope set but they always have a mass of zero:
    // RDKit❗❌:     if (!mass) {
    // RDKit❗❌:       mass = isotope;
    // RDKit❗❌:     }
    // RDKit❗❌:     ss << " MASS=" << mass;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int nRadEs = atom->getNumRadicalElectrons();
    // RDKit❗❌:   if (nRadEs != 0 && atom->getTotalDegree() != 0) {
    // RDKit❗❌:     if (nRadEs % 2) {
    // RDKit❗❌:       nRadEs = 2;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       nRadEs = 3;  // we use triplets, not singlets:
    // RDKit❗❌:     }
    // RDKit❗❌:     ss << " RAD=" << nRadEs;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (totValence != 0) {
    // RDKit❗❌:     if (totValence == 15) {
    // RDKit❗❌:       ss << " VAL=-1";
    // RDKit❗❌:     } else {
    // RDKit❗❌:       ss << " VAL=" << totValence;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (symbol == "R#") {
    // RDKit❗❌:     unsigned int rLabel = 1;
    // RDKit❗❌:     atom->getPropIfPresent(common_properties::_MolFileRLabel, rLabel);
    // RDKit❗❌:     ss << " RGROUPS=(1 " << rLabel << ")";
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   {
    // RDKit❗❌:     int iprop;
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molAttachOrder, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " ATTCHORD=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molAttachPoint, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " ATTCHPT=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molAtomSeqId, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " SEQID=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     {
    // RDKit❗❌:       std::string sprop;
    // RDKit❗❌:       if (atom->getPropIfPresent(common_properties::molAtomSeqName, sprop)) {
    // RDKit❗❌:         ss << " SEQNAME=" << sprop;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molRxnExactChange, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " EXACHG=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molInversionFlag, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       if (iprop == 1 || iprop == 2) {
    // RDKit❗❌:         ss << " INVRET=" << iprop;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molStereoCare, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " STBOX=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molSubstCount, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " SUBST=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molRingBondCount, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " RBCNT=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   {
    // RDKit❗❌:     std::string sprop;
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molAtomClass, sprop)) {
    // RDKit❗❌:       ss << " CLASS=" << sprop;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // HCOUNT - *query* hydrogen count. Not written by this writer.
    // RDKit❗❌:
    // RDKit❗❌:   return ss.str();
    // RDKit❗❌: }

    let symbol = v3000_atom_symbol(atom, molecule.query_atom(atom.id()))?;
    let parity_flag = parity_flags[atom.id().index()];
    let tot_valence = molfile_total_valence_field(molecule, atom)?;
    let mut out: PropertyText = format!("M  V30 {} ", atom.id().index() + 1).into();
    out.extend_bytes(symbol.as_bytes());
    out.extend_bytes(
        (&format!(
            " {0:.1$} {2:.1$} {3:.1$} {4}",
            coord[0],
            precision,
            coord[1],
            coord[2],
            atom.atom_map().unwrap_or(0)
        ))
            .as_ref(),
    );
    if parity_flag != 0 {
        out.extend_bytes((&format!(" CFG={parity_flag}")).as_ref());
    }
    if atom.formal_charge() != 0 {
        out.extend_bytes((&format!(" CHG={}", atom.formal_charge())).as_ref());
    }
    if let Some(isotope) = atom.isotope()
        && !is_atom_rgroup(atom)
    {
        out.extend_bytes((&format!(" MASS={isotope}")).as_ref());
    }
    let electrons = atom.radical_electrons();
    let valence = molblock_valence_assignment(molecule)?;
    if electrons != 0 && molfile_total_degree(molecule, atom, &valence) != 0 {
        let code = if electrons % 2 == 1 { 2 } else { 3 };
        out.extend_bytes((&format!(" RAD={code}")).as_ref());
    }
    if tot_valence != 0 {
        if tot_valence == 15 {
            out.extend_bytes((" VAL=-1").as_ref());
        } else {
            out.extend_bytes((&format!(" VAL={tot_valence}")).as_ref());
        }
    }
    if atom.prop("_MolFileRLabel").is_some() {
        let label = molfile_rlabel(atom)?.expect("present label");
        out.extend_bytes((&format!(" RGROUPS=(1 {label})")).as_ref());
    }
    crate::sdf::append_v3000_atom_properties(&mut out, atom, true)?;
    Ok(out)
}

fn is_atom_rgroup(atom: &Atom) -> bool {
    // BEGIN RDKIT CPP FUNCTION isAtomRGroup
    // RDKit✔️✔️: bool isAtomRGroup(const Atom &atom) {
    // RDKit✔️✔️:   return atom.getAtomicNum() == 0 &&
    // RDKit✔️✔️:          atom.hasProp(common_properties::_MolFileRLabel);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    // Exact dummy identity and property presence; typed label conversion remains
    // in molfile_rlabel. O(1) scalar access plus existing property lookup.
    atom.atomic_number() == 0 && atom.prop("_MolFileRLabel").is_some()
}

fn molfile_rlabel(atom: &Atom) -> Result<Option<u32>, MolWriteError> {
    // RDKit❗❌: template <>
    // RDKit❗❌: inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
    // RDKit❗❌:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit❗❌:     return v.value.u;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (rdvalue_is<int>(v)) {
    // RDKit❗❌:     return boost::numeric_cast<unsigned int>(v.value.i);
    // RDKit❗❌:   }
    // RDKit❗❌:   throw std::bad_any_cast();
    // RDKit❗❌: }
    // RDKit❗❌: template <class T>
    // RDKit❗❌: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit❗❌:     RDValue_cast_t arg) {
    // RDKit❗❌:   T res;
    // RDKit❗❌:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit❗❌:     Utils::LocaleSwitcher ls;
    // RDKit❗❌:     try {
    // RDKit❗❌:       res = rdvalue_cast<T>(arg);
    // RDKit❗❌:     } catch (const std::bad_any_cast &exc) {
    // RDKit❗❌:       try {
    // RDKit❗❌: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit❗❌: 	// trim only the right characters, this mimics how SD values
    // RDKit❗❌: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit❗❌: 	boost::trim_right(val);
    // RDKit❗❌:         res = boost::lexical_cast<T>(val);
    // RDKit❗❌:       } catch (...) {
    // RDKit❗❌:         throw exc;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res = rdvalue_cast<T>(arg);
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    let Some(value) = atom.prop("_MolFileRLabel") else {
        return Ok(None);
    };
    let label = cosmolkit_core::property_value_to_uint(value)?;
    Ok(Some(label))
}

fn v3000_atom_symbol(
    atom: &Atom,
    query: Option<&cosmolkit_model::QueryAtom>,
) -> Result<PropertyText, MolWriteError> {
    // RDKit❗❌: const std::string GetV3000MolFileAtomLine(
    // RDKit❗❌:     const Atom *atom, const Conformer *conf,
    // RDKit❗❌:     boost::dynamic_bitset<> &queryListAtoms, unsigned int precision) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   int totValence, atomMapNumber;
    // RDKit❗❌:   unsigned int parityFlag;
    // RDKit❗❌:   double x, y, z;
    // RDKit❗❌:   GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
    // RDKit❗❌:                            y, z);
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   ss << "M  V30 " << atom->getIdx() + 1;
    // RDKit❗❌:
    // RDKit❗❌:   std::string symbol = AtomGetMolFileSymbol(atom, false, queryListAtoms);
    // RDKit❗❌:   if (!isAtomListQuery(atom) || queryListAtoms[atom->getIdx()]) {
    // RDKit❗❌:     ss << " " << symbol;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     INT_VECT vals;
    // RDKit❗❌:     getAtomListQueryVals(atom->getQuery(), vals);
    // RDKit❗❌:     if (atom->getQuery()->getNegation()) {
    // RDKit❗❌:       ss << " "
    // RDKit❗❌:          << "\"NOT";
    // RDKit❗❌:     }
    // RDKit❗❌:     ss << " [";
    // RDKit❗❌:     for (unsigned int i = 0; i < vals.size(); ++i) {
    // RDKit❗❌:       if (i != 0) {
    // RDKit❗❌:         ss << ",";
    // RDKit❗❌:       }
    // RDKit❗❌:       ss << PeriodicTable::getTable()->getElementSymbol(vals[i]);
    // RDKit❗❌:     }
    // RDKit❗❌:     ss << "]";
    // RDKit❗❌:     if (atom->getQuery()->getNegation()) {
    // RDKit❗❌:       ss << "\"";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    if let Some(row) = query
        && query_atom_special_symbol(row).is_none()
        && let Some((numbers, negated)) = atom_list_query(row)
    {
        let values = numbers
            .into_iter()
            .map(molfile_atom_symbol)
            .collect::<Result<Vec<_>, _>>()?
            .join(",");
        return Ok(if negated {
            format!("\"NOT [{values}]\"").into()
        } else {
            format!("[{values}]").into()
        });
    }
    v2000_atom_symbol(atom, false, query)
}

fn v3000_bond_line(
    molecule: &MolWriteContext<'_>,
    bond: &Bond,
    include_stereo: bool,
    aromatic_bonds: &[usize],
    wedge_bonds: &WedgeAssignments,
    coords: Option<&[[f64; 3]]>,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetV3000MolFileBondLine(
    // RDKit❗❌:     const Bond *bond,
    // RDKit❗❌:     const std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds,
    // RDKit❗❌:     const Conformer *conf, bool wasAromatic) {
    // RDKit❗❌:   PRECONDITION(bond, "");
    // RDKit❗❌:
    // RDKit❗❌:   int dirCode = 0;
    // RDKit❗❌:   bool reverse = false;
    // RDKit❗❌:   RDKit::Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode,
    // RDKit❗❌:                                              reverse);
    // RDKit❗❌:   // do not cross bonds which were aromatic before kekulization
    // RDKit❗❌:   if (wasAromatic && dirCode == 3) {
    // RDKit❗❌:     dirCode = 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   ss << "M  V30 " << bond->getIdx() + 1;
    // RDKit❗❌:   ss << " " << GetV3000BondCode(bond);
    // RDKit❗❌:   if (reverse) {
    // RDKit❗❌:     // switch the begin and end atoms on the bond line
    // RDKit❗❌:     ss << " " << bond->getEndAtomIdx() + 1;
    // RDKit❗❌:     ss << " " << bond->getBeginAtomIdx() + 1;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     ss << " " << bond->getBeginAtomIdx() + 1;
    // RDKit❗❌:     ss << " " << bond->getEndAtomIdx() + 1;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (dirCode != 0) {
    // RDKit❗❌:     ss << " CFG=" << BondStereoCodeV2000ToV3000(dirCode);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (bond->hasQuery()) {
    // RDKit❗❌:     int topol = getQueryBondTopology(bond);
    // RDKit❗❌:     if (topol) {
    // RDKit❗❌:       ss << " TOPO=" << topol;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   {
    // RDKit❗❌:     int iprop;
    // RDKit❗❌:     if (bond->getPropIfPresent(common_properties::molReactStatus, iprop) &&
    // RDKit❗❌:         iprop) {
    // RDKit❗❌:       ss << " RXCTR=" << iprop;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   {
    // RDKit❗❌:     std::string sprop;
    // RDKit❗❌:     if (bond->getPropIfPresent(common_properties::molStereoCare, sprop) &&
    // RDKit❗❌:         sprop != "0") {
    // RDKit❗❌:       ss << " STBOX=" << sprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, sprop) &&
    // RDKit❗❌:         sprop != "0") {
    // RDKit❗❌:       ss << " ENDPTS=" << sprop;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (bond->getPropIfPresent(common_properties::_MolFileBondAttach, sprop) &&
    // RDKit❗❌:         sprop != "0") {
    // RDKit❗❌:       ss << " ATTACH=" << sprop;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return ss.str();
    // RDKit❗❌: }

    let (mut cfg, reverse) =
        v3000_bond_cfg_code(molecule, bond, wedge_bonds, coords, stereo_context)?;
    let (begin_idx, end_idx) = if reverse {
        (bond.end().index(), bond.begin().index())
    } else {
        (bond.begin().index(), bond.end().index())
    };
    let type_code = query_bond_type_code(molecule.query_bond(bond.id()))
        .filter(|code| *code != 0)
        .map_or_else(|| v3000_bond_type_code(bond), Ok)?;
    let mut out: PropertyText = format!(
        "M  V30 {} {} {} {}",
        bond.id().index() + 1,
        type_code,
        begin_idx + 1,
        end_idx + 1
    )
    .into();
    if aromatic_bonds.contains(&bond.id().index())
        && cfg == Some(2)
        && matches!(
            molfile_bond_stereo_info(molecule, bond, wedge_bonds, coords, stereo_context)?.0,
            BondDirection::EitherDouble
        )
    {
        cfg = None;
    }
    if let Some(cfg) = cfg {
        out.extend_bytes((&format!(" CFG={cfg}")).as_ref());
    }

    let topology = molecule
        .query_bond(bond.id())
        .map_or(0, |row| query_bond_topology(row.predicate()));
    if topology != 0 {
        out.extend_bytes((&format!(" TOPO={topology}")).as_ref());
    }
    crate::sdf::append_v3000_bond_properties(&mut out, bond)?;
    Ok(out)
}

fn v3000_bond_cfg_code(
    molecule: &MolWriteContext<'_>,
    bond: &Bond,
    wedge_bonds: &WedgeAssignments,
    coords: Option<&[[f64; 3]]>,
    stereo_context: MolfileStereoContext<'_>,
) -> Result<(Option<u32>, bool), MolWriteError> {
    //
    //
    let (dir, reverse) =
        molfile_bond_stereo_info(molecule, bond, wedge_bonds, coords, stereo_context)?;
    let cfg = match molfile_bond_dir_code(dir) {
        0 => None,
        1 => Some(1),
        3 | 4 => Some(2),
        6 => Some(3),
        _ => None,
    };
    Ok((cfg, reverse))
}

fn v3000_bond_type_code(bond: &Bond) -> Result<u32, MolWriteError> {
    match bond.order() {
        BondOrder::Single => Ok(if bond.is_aromatic() { 4 } else { 1 }),
        BondOrder::Double => Ok(if bond.is_aromatic() { 4 } else { 2 }),
        BondOrder::Triple => Ok(3),
        BondOrder::Aromatic => Ok(4),
        BondOrder::Zero => Ok(1),
        BondOrder::Dative => Ok(9),
        BondOrder::Hydrogen => Ok(10),
        _ => Ok(0),
    }
}

fn append_v3000_sgroup_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
    generated_sgroups: &[SubstanceGroup],
) -> Result<(), MolWriteError> {
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() && generated_sgroups.is_empty() {
        return Ok(());
    }
    out.extend_bytes(("M  V30 BEGIN SGROUP\n").as_ref());
    for (idx, sgroup) in sgroups.iter().enumerate() {
        out.extend_bytes(
            (&crate::sdf_sgroups::write_v3000_sgroup(idx + 1, sgroup, molecule.bonds())?).as_ref(),
        );
    }
    let offset = sgroups.len();
    for (idx, sgroup) in generated_sgroups.iter().enumerate() {
        out.extend_bytes(
            (&crate::sdf_sgroups::write_v3000_sgroup(offset + idx + 1, sgroup, molecule.bonds())?)
                .as_ref(),
        );
    }
    out.extend_bytes(("M  V30 END SGROUP\n").as_ref());
    Ok(())
}

fn v3000_generated_zbo_sgroups(
    molecule: &MolWriteContext<'_>,
    valence: &ValenceAssignment,
) -> Vec<SubstanceGroup> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: void createZBOSubstanceGroups(ROMol &mol) {
    // RDKit❗❌:   SubstanceGroup bsg(&mol, "DAT");
    // RDKit❗❌:   bsg.setProp("FIELDNAME", "ZBO");
    // RDKit❗❌:   boost::dynamic_bitset<> atomsAffected(mol.getNumAtoms(), 0);
    // RDKit❗❌:   for (const auto bond : mol.bonds()) {
    // RDKit❗❌:     if (bond->getBondType() == Bond::ZERO) {
    // RDKit❗❌:       bsg.addBondWithIdx(bond->getIdx());
    // RDKit❗❌:       atomsAffected[bond->getBeginAtomIdx()] = 1;
    // RDKit❗❌:       atomsAffected[bond->getEndAtomIdx()] = 1;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (atomsAffected.any()) {
    // RDKit❗❌:     for (auto i = 0u; i < atomsAffected.size(); ++i) {
    // RDKit❗❌:       if (atomsAffected[i]) {
    // RDKit❗❌:         bsg.addAtomWithIdx(i);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     SubstanceGroup asg(&mol, "DAT");
    // RDKit❗❌:     asg.setProp("FIELDNAME", "HYD");
    // RDKit❗❌:     SubstanceGroup zsg(&mol, "DAT");
    // RDKit❗❌:     zsg.setProp("FIELDNAME", "ZCH");
    // RDKit❗❌:     std::string asgText;
    // RDKit❗❌:     std::string zsgText;
    // RDKit❗❌:     for (auto i = 0u; i < atomsAffected.size(); ++i) {
    // RDKit❗❌:       if (atomsAffected[i]) {
    // RDKit❗❌:         const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit❗❌:         asg.addAtomWithIdx(i);
    // RDKit❗❌:         if (!asgText.empty()) {
    // RDKit❗❌:           asgText += ";";
    // RDKit❗❌:         }
    // RDKit❗❌:         asgText += std::to_string(atom->getTotalNumHs());
    // RDKit❗❌:         zsg.addAtomWithIdx(i);
    // RDKit❗❌:         if (!zsgText.empty()) {
    // RDKit❗❌:           zsgText += ";";
    // RDKit❗❌:         }
    // RDKit❗❌:         zsgText += std::to_string(atom->getFormalCharge());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     addSubstanceGroup(mol, bsg);
    // RDKit❗❌:
    // RDKit❗❌:     std::vector<std::string> aDataFields{asgText};
    // RDKit❗❌:
    // RDKit❗❌:     asg.setProp("DATAFIELDS", aDataFields);
    // RDKit❗❌:     addSubstanceGroup(mol, asg);
    // RDKit❗❌:     std::vector<std::string> zDataFields{zsgText};
    // RDKit❗❌:     zsg.setProp("DATAFIELDS", zDataFields);
    // RDKit❗❌:     addSubstanceGroup(mol, zsg);
    // RDKit❗❌:   }
    // RDKit❗❌: }

    let zero_bonds = molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero)
        .map(Bond::id)
        .collect::<Vec<_>>();
    if zero_bonds.is_empty() {
        return Vec::new();
    }
    let mut affected = vec![false; molecule.num_atoms()];
    for bond in molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero)
    {
        affected[bond.begin().index()] = true;
        affected[bond.end().index()] = true;
    }
    let atoms = affected
        .iter()
        .enumerate()
        .filter_map(|(idx, is_affected)| is_affected.then_some(AtomId::new(idx)))
        .collect::<Vec<_>>();
    let hydrogens = atoms
        .iter()
        .map(|atom| molfile_total_hydrogens(&molecule.atoms()[atom.index()], valence).to_string())
        .collect::<Vec<_>>()
        .join(";");
    let charges = atoms
        .iter()
        .map(|atom| molecule.atoms()[atom.index()].formal_charge().to_string())
        .collect::<Vec<_>>()
        .join(";");
    vec![
        SubstanceGroup::new(
            cosmolkit_model::SubstanceGroupId::new(0),
            SubstanceGroupKind::Data,
        )
        .with_atoms(atoms.clone())
        .with_bonds(zero_bonds)
        .with_data(SGroupData {
            field_name: Some(PropertyText::from("ZBO")),
            ..SGroupData::default()
        }),
        SubstanceGroup::new(
            cosmolkit_model::SubstanceGroupId::new(1),
            SubstanceGroupKind::Data,
        )
        .with_atoms(atoms.clone())
        .with_data(SGroupData {
            field_name: Some(PropertyText::from("HYD")),
            values: vec![PropertyText::from(hydrogens)],
            ..SGroupData::default()
        }),
        SubstanceGroup::new(
            cosmolkit_model::SubstanceGroupId::new(2),
            SubstanceGroupKind::Data,
        )
        .with_atoms(atoms)
        .with_data(SGroupData {
            field_name: Some(PropertyText::from("ZCH")),
            values: vec![PropertyText::from(charges)],
            ..SGroupData::default()
        }),
    ]
}

fn append_v3000_collection_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
    wedge_bonds: &WedgeAssignments,
) -> Result<(), MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: void appendEnhancedStereoGroups(
    // RDKit❗❌:     std::string &res, const RWMol &tmol,
    // RDKit❗❌:     std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds) {
    // RDKit❗❌:   if (!tmol.getStereoGroups().empty()) {
    // RDKit❗❌:     auto stereo_groups = tmol.getStereoGroups();
    // RDKit❗❌:     assignStereoGroupIds(stereo_groups);
    // RDKit❗❌:     res += "M  V30 BEGIN COLLECTION\n";
    // RDKit❗❌:     std::string tmp;
    // RDKit❗❌:     tmp.reserve(80);
    // RDKit❗❌:     for (auto &&group : stereo_groups) {
    // RDKit❗❌:       tmp += "M  V30 MDLV30/";
    // RDKit❗❌:       switch (group.getGroupType()) {
    // RDKit❗❌:         case RDKit::StereoGroupType::STEREO_ABSOLUTE:
    // RDKit❗❌:           tmp += "STEABS";
    // RDKit❗❌:           break;
    // RDKit❗❌:         case RDKit::StereoGroupType::STEREO_OR:
    // RDKit❗❌:           tmp += "STEREL";
    // RDKit❗❌:           tmp += std::to_string(group.getWriteId());
    // RDKit❗❌:           break;
    // RDKit❗❌:         case RDKit::StereoGroupType::STEREO_AND:
    // RDKit❗❌:           tmp += "STERAC";
    // RDKit❗❌:           tmp += std::to_string(group.getWriteId());
    // RDKit❗❌:           break;
    // RDKit❗❌:       }
    // RDKit❗❌:       tmp += " ATOMS=(";
    // RDKit❗❌:
    // RDKit❗❌:       std::vector<unsigned int> atomIds;
    // RDKit❗❌:       Atropisomers::getAllAtomIdsForStereoGroup(tmol, group, atomIds,
    // RDKit❗❌:                                                 wedgeBonds);
    // RDKit❗❌:
    // RDKit❗❌:       tmp += std::to_string(atomIds.size());
    // RDKit❗❌:       for (auto &&atom : atomIds) {
    // RDKit❗❌:         tmp += ' ';
    // RDKit❗❌:         // atoms are 1 indexed in molfiles
    // RDKit❗❌:         auto idxStr = std::to_string(atom + 1);
    // RDKit❗❌:         if (tmp.size() + idxStr.size() >= 78) {
    // RDKit❗❌:           res += tmp + "-\n";
    // RDKit❗❌:           tmp = "M  V30 ";
    // RDKit❗❌:         }
    // RDKit❗❌:         tmp += idxStr;
    // RDKit❗❌:       }
    // RDKit❗❌:       res += tmp + ")\n";
    // RDKit❗❌:       tmp.clear();
    // RDKit❗❌:     }
    // RDKit❗❌:     res += tmp + "M  V30 END COLLECTION\n";
    // RDKit❗❌:   }
    // RDKit❗❌: }

    let atoms = cosmolkit_core::get_all_atom_ids_for_stereo_groups(
        &molecule.topology,
        molecule.stereo_groups(),
        wedge_bonds,
    )?;
    out.extend_bytes(
        (&crate::sdf_sgroups::write_v3000_collection_rows(molecule.stereo_groups(), &atoms))
            .as_ref(),
    );
    Ok(())
}

fn append_v2000_property_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
) -> Result<(), MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFileChargeInfo(const RWMol &mol) {
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   std::stringstream chgss;
    // RDKit❗❌:   std::stringstream radss;
    // RDKit❗❌:   std::stringstream massdiffss;
    // RDKit❗❌:   unsigned int nChgs = 0;
    // RDKit❗❌:   unsigned int nRads = 0;
    // RDKit❗❌:   unsigned int nMassDiffs = 0;
    // RDKit❗❌:   for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
    // RDKit❗❌:        atomIt != mol.endAtoms(); ++atomIt) {
    // RDKit❗❌:     const Atom *atom = *atomIt;
    // RDKit❗❌:     if (atom->getFormalCharge() != 0) {
    // RDKit❗❌:       ++nChgs;
    // RDKit❗❌:       chgss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
    // RDKit❗❌:                    atom->getFormalCharge();
    // RDKit❗❌:       if (nChgs == 8) {
    // RDKit❗❌:         res << boost::format("M  CHG%3d") % nChgs << chgss.str() << "\n";
    // RDKit❗❌:         chgss.str("");
    // RDKit❗❌:         nChgs = 0;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     unsigned int nRadEs = atom->getNumRadicalElectrons();
    // RDKit❗❌:     if (nRadEs != 0 && atom->getTotalDegree() != 0) {
    // RDKit❗❌:       ++nRads;
    // RDKit❗❌:       if (nRadEs % 2) {
    // RDKit❗❌:         nRadEs = 2;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         nRadEs = 3;  // we use triplets, not singlets:
    // RDKit❗❌:       }
    // RDKit❗❌:       radss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % nRadEs;
    // RDKit❗❌:       if (nRads == 8) {
    // RDKit❗❌:         res << boost::format("M  RAD%3d") % nRads << radss.str() << "\n";
    // RDKit❗❌:         radss.str("");
    // RDKit❗❌:         nRads = 0;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!isAtomRGroup(*atom)) {
    // RDKit❗❌:       int isotope = atom->getIsotope();
    // RDKit❗❌:       if (isotope != 0) {
    // RDKit❗❌:         ++nMassDiffs;
    // RDKit❗❌:         massdiffss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
    // RDKit❗❌:                           isotope;
    // RDKit❗❌:         if (nMassDiffs == 8) {
    // RDKit❗❌:           res << boost::format("M  ISO%3d") % nMassDiffs << massdiffss.str()
    // RDKit❗❌:               << "\n";
    // RDKit❗❌:           massdiffss.str("");
    // RDKit❗❌:           nMassDiffs = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (nChgs) {
    // RDKit❗❌:     res << boost::format("M  CHG%3d") % nChgs << chgss.str() << "\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   if (nRads) {
    // RDKit❗❌:     res << boost::format("M  RAD%3d") % nRads << radss.str() << "\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   if (nMassDiffs) {
    // RDKit❗❌:     res << boost::format("M  ISO%3d") % nMassDiffs << massdiffss.str() << "\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   return res.str();
    // RDKit❗❌: }

    let charges = molecule
        .atoms()
        .iter()
        .filter_map(|atom| {
            (atom.formal_charge() != 0)
                .then_some((atom.id().index() + 1, i32::from(atom.formal_charge())))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "CHG", &charges);

    let valence = molblock_valence_assignment(molecule)?;
    let radicals = molecule
        .atoms()
        .iter()
        .filter_map(|atom| {
            let electrons = atom.radical_electrons();
            if electrons == 0 || molfile_total_degree(molecule, atom, &valence) == 0 {
                return None;
            }
            let code = if electrons % 2 == 1 { 2 } else { 3 };
            Some((atom.id().index() + 1, code))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "RAD", &radicals);

    let isotopes = molecule
        .atoms()
        .iter()
        .filter(|atom| !is_atom_rgroup(atom))
        .filter_map(|atom| {
            atom.isotope()
                .map(|isotope| (atom.id().index() + 1, i32::from(isotope)))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ISO", &isotopes);
    Ok(())
}

fn append_v2000_rgroup_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
) -> Result<(), MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFileRGroupInfo(const RWMol &mol) {
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   unsigned int nEntries = 0;
    // RDKit❗❌:   for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
    // RDKit❗❌:        atomIt != mol.endAtoms(); ++atomIt) {
    // RDKit❗❌:     unsigned int lbl;
    // RDKit❗❌:     if ((*atomIt)->getPropIfPresent(common_properties::_MolFileRLabel, lbl)) {
    // RDKit❗❌:       ss << " " << std::setw(3) << (*atomIt)->getIdx() + 1 << " "
    // RDKit❗❌:          << std::setw(3) << lbl;
    // RDKit❗❌:       ++nEntries;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::stringstream ss2;
    // RDKit❗❌:   if (nEntries) {
    // RDKit❗❌:     ss2 << "M  RGP" << std::setw(3) << nEntries << ss.str() << "\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   return ss2.str();
    // RDKit❗❌: }

    let entries = molecule
        .atoms()
        .iter()
        .filter(|atom| atom.prop("_MolFileRLabel").is_some())
        .map(|atom| {
            Ok((
                atom.id().index() + 1,
                molfile_rlabel(atom)?.expect("present label"),
            ))
        })
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    if !entries.is_empty() {
        out.extend_bytes((&format!("M  RGP{:>3}", entries.len())).as_ref());
        for (idx, label) in entries {
            out.extend_bytes((&format!(" {:>3} {:>3}", idx, label)).as_ref());
        }
        out.push_byte(b'\n');
    }
    Ok(())
}

fn append_v2000_value_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
) -> Result<(), MolWriteError> {
    // RDKit❗❌: const std::string GetMolFileQueryInfo(
    // RDKit❗❌:     const RWMol &mol, const boost::dynamic_bitset<> &queryListAtoms) {
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   boost::dynamic_bitset<> listQs(mol.getNumAtoms());
    // RDKit❗❌:   for (const auto atom : mol.atoms()) {
    // RDKit❗❌:     if (isAtomListQuery(atom) && !queryListAtoms[atom->getIdx()]) {
    // RDKit❗❌:       listQs.set(atom->getIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto atom : mol.atoms()) {
    // RDKit❗❌:     bool wrote_query = false;
    // RDKit❗❌:     if (!listQs[atom->getIdx()] && !queryListAtoms[atom->getIdx()] &&
    // RDKit❗❌:         hasComplexQuery(atom)) {
    // RDKit❗❌:       std::string sma =
    // RDKit❗❌:           SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(atom));
    // RDKit❗❌:       ss << "V  " << std::setw(3) << atom->getIdx() + 1 << " " << sma << "\n";
    // RDKit❗❌:       wrote_query = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string molFileValue;
    // RDKit❗❌:     if (!wrote_query &&
    // RDKit❗❌:         atom->getPropIfPresent(common_properties::molFileValue, molFileValue)) {
    // RDKit❗❌:       ss << "V  " << std::setw(3) << atom->getIdx() + 1 << " " << molFileValue
    // RDKit❗❌:          << "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto atom : mol.atoms()) {
    // RDKit❗❌:     if (listQs[atom->getIdx()]) {
    // RDKit❗❌:       INT_VECT vals;
    // RDKit❗❌:       getAtomListQueryVals(atom->getQuery(), vals);
    // RDKit❗❌:       ss << "M  ALS " << std::setw(3) << atom->getIdx() + 1 << " ";
    // RDKit❗❌:       ss << std::setw(2) << vals.size();
    // RDKit❗❌:       if (atom->getQuery()->getNegation()) {
    // RDKit❗❌:         ss << " T ";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         ss << " F ";
    // RDKit❗❌:       }
    // RDKit❗❌:       for (auto val : vals) {
    // RDKit❗❌:         ss << std::setw(4) << std::left
    // RDKit❗❌:            << (PeriodicTable::getTable()->getElementSymbol(val));
    // RDKit❗❌:       }
    // RDKit❗❌:       ss << "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return ss.str();
    // RDKit❗❌: }
    for atom in molecule.atoms() {
        let query = molecule.query_atom(atom.id());
        let complex = query.is_some_and(|row| {
            query_atom_special_symbol(row).is_none()
                && atom_list_query(row).is_none()
                && has_complex_atom_query(row)
        });
        if complex {
            let sma = cosmolkit_search::query_atom_to_smarts(
                query.expect("checked query"),
                &Default::default(),
            )?;
            out.extend_bytes(format!("V  {:>3} ", atom.id().index() + 1).as_bytes());
            out.extend_bytes(sma.as_bytes());
            out.push_byte(b'\n');
        } else if let Some(value) = atom.prop("molFileValue") {
            let value = crate::sdf::model_string_property(value)?;
            out.extend_bytes(format!("V  {:>3} ", atom.id().index() + 1).as_bytes());
            out.extend_bytes(value.as_bytes());
            out.push_byte(b'\n');
        }
    }
    for atom in molecule.atoms() {
        if let Some(query) = molecule.query_atom(atom.id())
            && query_atom_special_symbol(query).is_none()
            && let Some((numbers, negated)) = atom_list_query(query)
        {
            out.extend_bytes(
                (&format!(
                    "M  ALS {:>3} {:>2} {} ",
                    atom.id().index() + 1,
                    numbers.len(),
                    if negated { "T" } else { "F" }
                ))
                    .as_ref(),
            );
            for n in numbers {
                out.extend_bytes((&format!("{:<4}", molfile_atom_symbol(n)?)).as_ref());
            }
            out.push_byte(b'\n');
        }
    }
    Ok(())
}

fn append_v2000_zbo_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
    valence: &ValenceAssignment,
) {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFileZBOInfo(const RWMol &mol) {
    // RDKit❗❌:   std::stringstream res;
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   unsigned int nEntries = 0;
    // RDKit❗❌:   boost::dynamic_bitset<> atomsAffected(mol.getNumAtoms(), 0);
    // RDKit❗❌:   for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
    // RDKit❗❌:        bondIt != mol.endBonds(); ++bondIt) {
    // RDKit❗❌:     if ((*bondIt)->getBondType() == Bond::ZERO) {
    // RDKit❗❌:       ++nEntries;
    // RDKit❗❌:       ss << " " << std::setw(3) << (*bondIt)->getIdx() + 1 << " "
    // RDKit❗❌:          << std::setw(3) << 0;
    // RDKit❗❌:       if (nEntries == 8) {
    // RDKit❗❌:         res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << "\n";
    // RDKit❗❌:         nEntries = 0;
    // RDKit❗❌:         ss.str("");
    // RDKit❗❌:       }
    // RDKit❗❌:       atomsAffected[(*bondIt)->getBeginAtomIdx()] = 1;
    // RDKit❗❌:       atomsAffected[(*bondIt)->getEndAtomIdx()] = 1;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (nEntries) {
    // RDKit❗❌:     res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << "\n";
    // RDKit❗❌:   }
    // RDKit❗❌:   if (atomsAffected.count()) {
    // RDKit❗❌:     std::stringstream hydss;
    // RDKit❗❌:     unsigned int nhyd = 0;
    // RDKit❗❌:     std::stringstream zchss;
    // RDKit❗❌:     unsigned int nzch = 0;
    // RDKit❗❌:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:       if (!atomsAffected[i]) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit❗❌:       nhyd++;
    // RDKit❗❌:       hydss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
    // RDKit❗❌:                    atom->getTotalNumHs();
    // RDKit❗❌:       if (nhyd == 8) {
    // RDKit❗❌:         res << boost::format("M  HYD%3d") % nhyd << hydss.str() << "\n";
    // RDKit❗❌:         hydss.str("");
    // RDKit❗❌:         nhyd = 0;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (atom->getFormalCharge()) {
    // RDKit❗❌:         nzch++;
    // RDKit❗❌:         zchss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
    // RDKit❗❌:                      atom->getFormalCharge();
    // RDKit❗❌:         if (nzch == 8) {
    // RDKit❗❌:           res << boost::format("M  ZCH%3d") % nzch << zchss.str() << "\n";
    // RDKit❗❌:           zchss.str("");
    // RDKit❗❌:           nzch = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (nhyd) {
    // RDKit❗❌:       res << boost::format("M  HYD%3d") % nhyd << hydss.str() << "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:     if (nzch) {
    // RDKit❗❌:       res << boost::format("M  ZCH%3d") % nzch << zchss.str() << "\n";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res.str();
    // RDKit❗❌: }

    let zbo_entries = molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero || bond.prop("_ZBO").is_some())
        .map(|bond| (bond.id().index() + 1, 0))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ZBO", &zbo_entries);
    if zbo_entries.is_empty() {
        return;
    }
    let mut affected = vec![false; molecule.num_atoms()];
    for bond in molecule
        .bonds()
        .iter()
        .filter(|bond| bond.order() == BondOrder::Zero || bond.prop("_ZBO").is_some())
    {
        affected[bond.begin().index()] = true;
        affected[bond.end().index()] = true;
    }
    let hydrogens = molecule
        .atoms()
        .iter()
        .filter(|atom| affected[atom.id().index()])
        .map(|atom| {
            (
                atom.id().index() + 1,
                molfile_total_hydrogens(atom, valence),
            )
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "HYD", &hydrogens);
    let zcharges = molecule
        .atoms()
        .iter()
        .filter(|atom| affected[atom.id().index()] && atom.formal_charge() != 0)
        .map(|atom| (atom.id().index() + 1, i32::from(atom.formal_charge())))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ZCH", &zcharges);
}

fn append_v2000_pxa_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
) -> Result<(), MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFilePXAInfo(const RWMol &mol) {
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   for (const auto atom : mol.atoms()) {
    // RDKit❗❌:     if (atom->hasProp("_MolFile_PXA")) {
    // RDKit❗❌:       res +=
    // RDKit❗❌:           boost::str(boost::format("M  PXA % 3d%s\n") % (atom->getIdx() + 1) %
    // RDKit❗❌:                      atom->getProp<std::string>("_MolFile_PXA"));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    for atom in molecule.atoms() {
        if let Some(pxa) = atom.prop("_MolFile_PXA") {
            let pxa = crate::sdf::model_string_property(pxa)?;
            out.extend_bytes(format!("M  PXA {:>3}", atom.id().index() + 1).as_bytes());
            out.extend_bytes(pxa.as_bytes());
            out.push_byte(b'\n');
        }
    }
    Ok(())
}

fn append_v2000_sgroup_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
) -> Result<(), MolWriteError> {
    out.extend_bytes((&crate::sdf_sgroups::write_v2000_sgroups(&molecule.topology)?).as_ref());
    Ok(())
}

fn v2000_int_field(value: usize) -> String {
    format!(" {value:>3}")
}

fn v2000_double_field(value: f64) -> String {
    format!("{value:>10.4}")
}

fn v2000_string_field(
    value: &[u8],
    field_size: usize,
    pad: bool,
    add_separator: bool,
) -> PropertyText {
    let mut out = PropertyText::new();
    if add_separator {
        out.push_byte(b' ');
    }
    out.extend_bytes(&value[..value.len().min(field_size)]);
    if pad {
        for _ in value.len()..field_size {
            out.push_byte(b' ');
        }
    }
    out
}

fn append_v2000_alias_lines(
    out: &mut PropertyText,
    molecule: &MolWriteContext<'_>,
) -> Result<(), MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string GetMolFileAliasInfo(const RWMol &mol) {
    // RDKit❗❌:   std::stringstream ss;
    // RDKit❗❌:   for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
    // RDKit❗❌:        atomIt != mol.endAtoms(); ++atomIt) {
    // RDKit❗❌:     std::string lbl;
    // RDKit❗❌:     if ((*atomIt)->getPropIfPresent(common_properties::molFileAlias, lbl)) {
    // RDKit❗❌:       if (!lbl.empty()) {
    // RDKit❗❌:         ss << "A  " << std::setw(3) << (*atomIt)->getIdx() + 1 << "\n"
    // RDKit❗❌:            << lbl << "\n";
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return ss.str();
    // RDKit❗❌: }

    for atom in molecule.atoms() {
        if let Some(alias) = atom.prop("molFileAlias") {
            let alias = crate::sdf::model_string_property(alias)?;
            if alias.is_empty() {
                continue;
            }
            out.extend_bytes(format!("A  {:>3}\n", atom.id().index() + 1).as_bytes());
            out.extend_bytes(alias.as_bytes());
            out.push_byte(b'\n');
        }
    }
    Ok(())
}

fn append_v2000_counted_property(out: &mut PropertyText, label: &str, entries: &[(usize, i32)]) {
    for chunk in entries.chunks(8) {
        out.extend_bytes((&format!("M  {label}{:>3}", chunk.len())).as_ref());
        for (idx, value) in chunk {
            out.extend_bytes((&format!(" {:>3} {:>3}", idx, value)).as_ref());
        }
        out.push_byte(b'\n');
    }
}

fn atom_degree(molecule: &MolWriteContext<'_>, atom_index: usize) -> usize {
    molecule
        .bonds()
        .iter()
        .filter(|bond| bond.begin().index() == atom_index || bond.end().index() == atom_index)
        .count()
}

fn append_sdf_record_fields(
    mut block: PropertyText,
    molecule: &MolWriteContext<'_>,
) -> PropertyText {
    // RDKit❗✔️: void _writePropToStream(std::ostream *dp_ostream, const ROMol &mol,
    // RDKit❗✔️:                         const std::string &name, int d_molid) {
    // RDKit❗✔️:   PRECONDITION(dp_ostream, "no output stream");
    // RDKit❗✔️:
    // RDKit❗✔️:   // write the property value
    // RDKit❗✔️:   // FIX: we will assume for now that the desired property value is
    // RDKit❗✔️:   // catable to a string
    // RDKit❗✔️:   std::string pval;
    // RDKit❗✔️:   try {
    // RDKit❗✔️:     mol.getProp(name, pval);
    // RDKit❗✔️:   } catch (std::bad_any_cast &) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // warn and skip if we include a new line
    // RDKit❗✔️:   if (name.find("\n") != std::string::npos) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "WARNING: Skipping property " << name
    // RDKit❗✔️:         << " because the name includes a newline" << std::endl;
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (pval.find("\r\n\r\n") != std::string::npos ||
    // RDKit❗✔️:       pval.find("\n\n") != std::string::npos) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "WARNING: Skipping property " << name
    // RDKit❗✔️:         << " because the value includes an illegal blank line" << std::endl;
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // write the property header line
    // RDKit❗✔️:   (*dp_ostream) << ">  <" << name << ">  ";
    // RDKit❗✔️:   if (d_molid >= 0) {
    // RDKit❗✔️:     (*dp_ostream) << "(" << d_molid + 1 << ") ";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   (*dp_ostream) << "\n";
    // RDKit❗✔️:
    // RDKit❗✔️:   (*dp_ostream) << pval << "\n";
    // RDKit❗✔️:
    // RDKit❗✔️:   // empty line after the property
    // RDKit❗✔️:   (*dp_ostream) << "\n";
    // RDKit❗✔️: }
    // RDKit❗✔️: void _MolToSDStream(std::ostream *dp_ostream, const ROMol &mol, int confId,
    // RDKit❗✔️:                     bool df_kekulize, bool df_forceV3000, int d_molid,
    // RDKit❗✔️:                     STR_VECT *props) {
    // RDKit❗✔️:   PRECONDITION(dp_ostream, "no output stream");
    // RDKit❗✔️:
    // RDKit❗✔️:   // write the molecule
    // RDKit❗✔️:   (*dp_ostream) << MolToMolBlock(mol, true, confId, df_kekulize, df_forceV3000);
    // RDKit❗✔️:
    // RDKit❗✔️:   // now write the properties
    // RDKit❗✔️:   STR_VECT_CI pi;
    // RDKit❗✔️:   if (props && props->size() > 0) {
    // RDKit❗✔️:     // check if we have any properties the user specified to write out
    // RDKit❗✔️:     // in which loop over them and write them out
    // RDKit❗✔️:     for (pi = props->begin(); pi != props->end(); pi++) {
    // RDKit❗✔️:       if (mol.hasProp(*pi)) {
    // RDKit❗✔️:         _writePropToStream(dp_ostream, mol, (*pi), d_molid);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // if use did not specify any properties, write all non computed properties
    // RDKit❗✔️:     // out to the file
    // RDKit❗✔️:     STR_VECT properties = mol.getPropList();
    // RDKit❗✔️:     STR_VECT compLst;
    // RDKit❗✔️:     mol.getPropIfPresent(RDKit::detail::computedPropName, compLst);
    // RDKit❗✔️:
    // RDKit❗✔️:     STR_VECT_CI pi;
    // RDKit❗✔️:     for (pi = properties.begin(); pi != properties.end(); pi++) {
    // RDKit❗✔️:       // ignore any of the following properties
    // RDKit❗✔️:       if (((*pi) == RDKit::detail::computedPropName) ||
    // RDKit❗✔️:           ((*pi) == common_properties::_Name) || ((*pi) == "_MolFileInfo") ||
    // RDKit❗✔️:           ((*pi) == "_MolFileComments") ||
    // RDKit❗✔️:           ((*pi) == common_properties::_MolFileChiralFlag)) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       // check if this property is not computed
    // RDKit❗✔️:       if (std::find(compLst.begin(), compLst.end(), (*pi)) == compLst.end()) {
    // RDKit❗✔️:         _writePropToStream(dp_ostream, mol, (*pi), d_molid);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // add the $$$$ that marks the end of a molecule
    // RDKit❗✔️:   (*dp_ostream) << "$$$$\n";
    // RDKit❗✔️: }
    // Modeled detached ordered string fields use the source newline/blank-line
    // guards and framing with d_molid < 0. Full source bodies above also retain
    // property enumeration/computed-property branches for independent review.
    // Local cost review: one field traversal, linear name/value substring
    // checks, amortized String appends; no per-field molecule clone or nested
    // field scan. Existing behavior and unsupported boundaries are unchanged.
    for (name, value) in molecule.properties().sdf_data_fields() {
        if name.as_bytes().contains(&b'\n')
            || value.as_bytes().windows(4).any(|w| w == b"\r\n\r\n")
            || value.as_bytes().windows(2).any(|w| w == b"\n\n")
        {
            continue;
        }
        block.extend_bytes((">  <").as_ref());
        block.extend_bytes((name).as_ref());
        block.extend_bytes((">  \n").as_ref());
        block.extend_bytes((value).as_ref());
        block.extend_bytes(("\n\n").as_ref());
    }
    block.extend_bytes(("$$$$\n").as_ref());
    block
}

fn v2000_atom_symbol(
    atom: &Atom,
    pad_with_spaces: bool,
    query: Option<&cosmolkit_model::QueryAtom>,
) -> Result<PropertyText, MolWriteError> {
    // Full pinned source; detached adaptation and acceptance remain under review.
    // RDKit❗❌: const std::string AtomGetMolFileSymbol(
    // RDKit❗❌:     const Atom *atom, bool padWithSpaces,
    // RDKit❗❌:     boost::dynamic_bitset<> &queryListAtoms) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   if (atom->hasProp(common_properties::_MolFileRLabel)) {
    // RDKit❗❌:     res = "R#";
    // RDKit❗❌:     //    } else if(!atom->hasQuery() && atom->getAtomicNum()){
    // RDKit❗❌:   } else if (atom->getAtomicNum()) {
    // RDKit❗❌:     res = atom->getSymbol();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     if (!atom->hasProp(common_properties::dummyLabel)) {
    // RDKit❗❌:       if (atom->hasQuery() &&
    // RDKit❗❌:           (atom->getQuery()->getTypeLabel() == "A" ||
    // RDKit❗❌:            (atom->getQuery()->getNegation() &&
    // RDKit❗❌:             atom->getQuery()->getDescription() == "AtomAtomicNum" &&
    // RDKit❗❌:             static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() ==
    // RDKit❗❌:                 1))) {
    // RDKit❗❌:         res = "A";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() &&
    // RDKit❗❌:                  (atom->getQuery()->getTypeLabel() == "Q" ||
    // RDKit❗❌:                   (atom->getQuery()->getNegation() &&
    // RDKit❗❌:                    atom->getQuery()->getDescription() == "AtomOr" &&
    // RDKit❗❌:                    atom->getQuery()->endChildren() -
    // RDKit❗❌:                            atom->getQuery()->beginChildren() ==
    // RDKit❗❌:                        2 &&
    // RDKit❗❌:                    (*atom->getQuery()->beginChildren())->getDescription() ==
    // RDKit❗❌:                        "AtomAtomicNum" &&
    // RDKit❗❌:                    static_cast<ATOM_EQUALS_QUERY *>(
    // RDKit❗❌:                        (*atom->getQuery()->beginChildren()).get())
    // RDKit❗❌:                            ->getVal() == 6 &&
    // RDKit❗❌:                    (*++(atom->getQuery()->beginChildren()))->getDescription() ==
    // RDKit❗❌:                        "AtomAtomicNum" &&
    // RDKit❗❌:                    static_cast<ATOM_EQUALS_QUERY *>(
    // RDKit❗❌:                        (*++(atom->getQuery()->beginChildren())).get())
    // RDKit❗❌:                            ->getVal() == 1))) {
    // RDKit❗❌:         res = "Q";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "X") {
    // RDKit❗❌:         res = "X";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "M") {
    // RDKit❗❌:         res = "M";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "AH") {
    // RDKit❗❌:         res = "AH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "QH") {
    // RDKit❗❌:         res = "QH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "XH") {
    // RDKit❗❌:         res = "XH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "MH") {
    // RDKit❗❌:         res = "MH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (hasComplexQuery(atom)) {
    // RDKit❗❌:         if (isAtomListQuery(atom)) {
    // RDKit❗❌:           res = "L";
    // RDKit❗❌:         } else {
    // RDKit❗❌:           res = "*";
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res = "R";
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       std::string symb;
    // RDKit❗❌:       atom->getProp(common_properties::dummyLabel, symb);
    // RDKit❗❌:       if (symb == "*") {
    // RDKit❗❌:         res = "R";
    // RDKit❗❌:       } else if (symb == "X") {
    // RDKit❗❌:         res = "R";
    // RDKit❗❌:       } else if (symb == "Xa") {
    // RDKit❗❌:         res = "R1";
    // RDKit❗❌:       } else if (symb == "Xb") {
    // RDKit❗❌:         res = "R2";
    // RDKit❗❌:       } else if (symb == "Xc") {
    // RDKit❗❌:         res = "R3";
    // RDKit❗❌:       } else if (symb == "Xd") {
    // RDKit❗❌:         res = "R4";
    // RDKit❗❌:       } else if (symb == "Xf") {
    // RDKit❗❌:         res = "R5";
    // RDKit❗❌:       } else if (symb == "Xg") {
    // RDKit❗❌:         res = "R6";
    // RDKit❗❌:       } else if (symb == "Xh") {
    // RDKit❗❌:         res = "R7";
    // RDKit❗❌:       } else if (symb == "Xi") {
    // RDKit❗❌:         res = "R8";
    // RDKit❗❌:       } else if (symb == "Xj") {
    // RDKit❗❌:         res = "R9";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res = symb;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // pad the end with spaces
    // RDKit❗❌:   if (padWithSpaces) {
    // RDKit❗❌:     while (res.size() < 3) {
    // RDKit❗❌:       res += " ";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }

    let mut symbol: PropertyText = if atom.prop("_MolFileRLabel").is_some() {
        "R#".into()
    } else if atom.atomic_number() != 0 {
        molfile_atom_symbol(atom.atomic_number())?.into()
    } else if let Some(dummy_label) = atom.prop("dummyLabel") {
        let label = crate::sdf::model_string_property(dummy_label)?;
        match label.as_bytes() {
            b"*" | b"X" => "R".into(),
            b"Xa" => "R1".into(),
            b"Xb" => "R2".into(),
            b"Xc" => "R3".into(),
            b"Xd" => "R4".into(),
            b"Xf" => "R5".into(),
            b"Xg" => "R6".into(),
            b"Xh" => "R7".into(),
            b"Xi" => "R8".into(),
            b"Xj" => "R9".into(),
            _ => label,
        }
    } else if let Some(row) = query {
        if let Some(symbol) = query_atom_special_symbol(row) {
            symbol.into()
        } else if has_complex_atom_query(row) {
            if atom_list_query(row).is_some() {
                "L".into()
            } else {
                "*".into()
            }
        } else {
            "R".into()
        }
    } else {
        "R".into()
    };
    if pad_with_spaces {
        while symbol.len() < 3 {
            symbol.push_byte(b' ');
        }
    }
    Ok(symbol)
}

fn query_atom_special_symbol(row: &cosmolkit_model::QueryAtom) -> Option<&'static str> {
    // RDKit❗❌: const std::string AtomGetMolFileSymbol(
    // RDKit❗❌:     const Atom *atom, bool padWithSpaces,
    // RDKit❗❌:     boost::dynamic_bitset<> &queryListAtoms) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:   if (atom->hasProp(common_properties::_MolFileRLabel)) {
    // RDKit❗❌:     res = "R#";
    // RDKit❗❌:     //    } else if(!atom->hasQuery() && atom->getAtomicNum()){
    // RDKit❗❌:   } else if (atom->getAtomicNum()) {
    // RDKit❗❌:     res = atom->getSymbol();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     if (!atom->hasProp(common_properties::dummyLabel)) {
    // RDKit❗❌:       if (atom->hasQuery() &&
    // RDKit❗❌:           (atom->getQuery()->getTypeLabel() == "A" ||
    // RDKit❗❌:            (atom->getQuery()->getNegation() &&
    // RDKit❗❌:             atom->getQuery()->getDescription() == "AtomAtomicNum" &&
    // RDKit❗❌:             static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() ==
    // RDKit❗❌:                 1))) {
    // RDKit❗❌:         res = "A";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() &&
    // RDKit❗❌:                  (atom->getQuery()->getTypeLabel() == "Q" ||
    // RDKit❗❌:                   (atom->getQuery()->getNegation() &&
    // RDKit❗❌:                    atom->getQuery()->getDescription() == "AtomOr" &&
    // RDKit❗❌:                    atom->getQuery()->endChildren() -
    // RDKit❗❌:                            atom->getQuery()->beginChildren() ==
    // RDKit❗❌:                        2 &&
    // RDKit❗❌:                    (*atom->getQuery()->beginChildren())->getDescription() ==
    // RDKit❗❌:                        "AtomAtomicNum" &&
    // RDKit❗❌:                    static_cast<ATOM_EQUALS_QUERY *>(
    // RDKit❗❌:                        (*atom->getQuery()->beginChildren()).get())
    // RDKit❗❌:                            ->getVal() == 6 &&
    // RDKit❗❌:                    (*++(atom->getQuery()->beginChildren()))->getDescription() ==
    // RDKit❗❌:                        "AtomAtomicNum" &&
    // RDKit❗❌:                    static_cast<ATOM_EQUALS_QUERY *>(
    // RDKit❗❌:                        (*++(atom->getQuery()->beginChildren())).get())
    // RDKit❗❌:                            ->getVal() == 1))) {
    // RDKit❗❌:         res = "Q";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "X") {
    // RDKit❗❌:         res = "X";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "M") {
    // RDKit❗❌:         res = "M";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "AH") {
    // RDKit❗❌:         res = "AH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "QH") {
    // RDKit❗❌:         res = "QH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "XH") {
    // RDKit❗❌:         res = "XH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (atom->hasQuery() && atom->getQuery()->getTypeLabel() == "MH") {
    // RDKit❗❌:         res = "MH";
    // RDKit❗❌:         queryListAtoms.set(atom->getIdx());
    // RDKit❗❌:       } else if (hasComplexQuery(atom)) {
    // RDKit❗❌:         if (isAtomListQuery(atom)) {
    // RDKit❗❌:           res = "L";
    // RDKit❗❌:         } else {
    // RDKit❗❌:           res = "*";
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res = "R";
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       std::string symb;
    // RDKit❗❌:       atom->getProp(common_properties::dummyLabel, symb);
    // RDKit❗❌:       if (symb == "*") {
    // RDKit❗❌:         res = "R";
    // RDKit❗❌:       } else if (symb == "X") {
    // RDKit❗❌:         res = "R";
    // RDKit❗❌:       } else if (symb == "Xa") {
    // RDKit❗❌:         res = "R1";
    // RDKit❗❌:       } else if (symb == "Xb") {
    // RDKit❗❌:         res = "R2";
    // RDKit❗❌:       } else if (symb == "Xc") {
    // RDKit❗❌:         res = "R3";
    // RDKit❗❌:       } else if (symb == "Xd") {
    // RDKit❗❌:         res = "R4";
    // RDKit❗❌:       } else if (symb == "Xf") {
    // RDKit❗❌:         res = "R5";
    // RDKit❗❌:       } else if (symb == "Xg") {
    // RDKit❗❌:         res = "R6";
    // RDKit❗❌:       } else if (symb == "Xh") {
    // RDKit❗❌:         res = "R7";
    // RDKit❗❌:       } else if (symb == "Xi") {
    // RDKit❗❌:         res = "R8";
    // RDKit❗❌:       } else if (symb == "Xj") {
    // RDKit❗❌:         res = "R9";
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res = symb;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // pad the end with spaces
    // RDKit❗❌:   if (padWithSpaces) {
    // RDKit❗❌:     while (res.size() < 3) {
    // RDKit❗❌:       res += " ";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: namespace {
    // RDKit❗❌: unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf) {
    // RDKit❗❌:   PRECONDITION(atom, "bad atom");
    // RDKit❗❌:   PRECONDITION(conf, "bad conformer");
    // RDKit❗❌:   if (!conf->is3D() ||
    // RDKit❗❌:       !(atom->getDegree() >= 3 && atom->getTotalDegree() == 4)) {
    // RDKit❗❌:     return 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   const ROMol &mol = atom->getOwningMol();
    // RDKit❗❌:   RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:   std::vector<std::pair<unsigned int, RDGeom::Point3D>> vs;
    // RDKit❗❌:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit❗❌:   boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit❗❌:   while (nbrIdx != endNbrs) {
    // RDKit❗❌:     const Atom *at = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:     unsigned int idx = at->getIdx();
    // RDKit❗❌:     RDGeom::Point3D v = conf->getAtomPos(idx);
    // RDKit❗❌:     v -= pos;
    // RDKit❗❌:     if (at->getAtomicNum() == 1) {
    // RDKit❗❌:       idx += mol.getNumAtoms();
    // RDKit❗❌:     }
    // RDKit❗❌:     vs.emplace_back(idx, v);
    // RDKit❗❌:     ++nbrIdx;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::sort(vs.begin(), vs.end(), Rankers::pairLess);
    // RDKit❗❌:   double vol;
    // RDKit❗❌:   if (vs.size() == 4) {
    // RDKit❗❌:     vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     vol = -vs[0].second.crossProduct(vs[1].second).dotProduct(vs[2].second);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (vol < 0) {
    // RDKit❗❌:     return 2;
    // RDKit❗❌:   } else if (vol > 0) {
    // RDKit❗❌:     return 1;
    // RDKit❗❌:   }
    // RDKit❗❌:   return 0;
    // RDKit❗❌: }
    // RDKit❗❌: }  // namespace
    // RDKit❗❌:
    // RDKit❗❌: bool hasNonDefaultValence(const Atom *atom) {
    // RDKit❗❌:   if (atom->getNumRadicalElectrons() != 0) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   // for queries and atoms which don't have computed properties, the answer is
    // RDKit❗❌:   // always no:
    // RDKit❗❌:   if (atom->hasQuery() || atom->needsUpdatePropertyCache()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (atom->getAtomicNum() == 1 ||
    // RDKit❗❌:       SmilesWrite ::inOrganicSubset(atom->getAtomicNum())) {
    // RDKit❗❌:     // for the ones we "know", we may have to specify the valence if it's
    // RDKit❗❌:     // not the default value
    // RDKit❗❌:     auto effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
    // RDKit❗❌:     return atom->getNoImplicit() &&
    // RDKit❗❌:            (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) !=
    // RDKit❗❌:             PeriodicTable::getTable()->getDefaultValence(effAtomicNum));
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: void GetMolFileAtomProperties(const Atom *atom, const Conformer *conf,
    // RDKit❗❌:                               int &totValence, int &atomMapNumber,
    // RDKit❗❌:                               unsigned int &parityFlag, double &x, double &y,
    // RDKit❗❌:                               double &z) {
    // RDKit❗❌:   PRECONDITION(atom, "");
    // RDKit❗❌:   totValence = 0;
    // RDKit❗❌:   atomMapNumber = 0;
    // RDKit❗❌:   parityFlag = 0;
    // RDKit❗❌:   x = y = z = 0.0;
    // RDKit❗❌:
    // RDKit❗❌:   if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                               atomMapNumber)) {
    // RDKit❗❌:     // XXX FIX ME->should we fail here? previously we would not assign
    // RDKit❗❌:     // the atomMapNumber if it didn't exist which could result in garbage
    // RDKit❗❌:     //  values.
    // RDKit❗❌:     atomMapNumber = 0;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    // RDKit❗❌:     x = pos.x;
    // RDKit❗❌:     y = pos.y;
    // RDKit❗❌:     z = pos.z;
    // RDKit❗❌:     if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:         atom->getChiralTag() != Atom::CHI_OTHER && atom->getDegree() >= 3 &&
    // RDKit❗❌:         atom->getTotalDegree() == 4) {
    // RDKit❗❌:       parityFlag = getAtomParityFlag(atom, conf);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (hasNonDefaultValence(atom)) {
    // RDKit❗❌:     if (atom->getTotalDegree() == 0) {
    // RDKit❗❌:       // Specify zero valence for elements/metals without neighbors
    // RDKit❗❌:       // or hydrogens (degree 0) instead of writing them as radicals.
    // RDKit❗❌:       totValence = 15;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // write the total valence for other atoms
    // RDKit❗❌:       totValence = atom->getTotalValence() % 15;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    match row.predicate() {
        QueryNode::Not(child)
            if matches!(
                child.as_ref(),
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(1))
            ) =>
        {
            Some("A")
        }
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberNotIn(numbers))
            if numbers == &[6, 1] =>
        {
            Some("Q")
        }
        _ => None,
    }
}

fn has_complex_atom_query(row: &cosmolkit_model::QueryAtom) -> bool {
    // RDKit❗❌: bool hasComplexQuery(const Atom *atom) {
    // RDKit❗❌:   PRECONDITION(atom, "bad atom");
    // RDKit❗❌:   bool res = false;
    // RDKit❗❌:   if (atom->hasQuery()) {
    // RDKit❗❌:     res = true;
    // RDKit❗❌:     // counter examples:
    // RDKit❗❌:     //  1) atomic number
    // RDKit❗❌:     //  2) the smarts parser inserts AtomAnd queries
    // RDKit❗❌:     //     for "C" or "c":
    // RDKit❗❌:     //
    // RDKit❗❌:     std::string descr = atom->getQuery()->getDescription();
    // RDKit❗❌:     if (descr == "AtomAtomicNum" &&
    // RDKit❗❌:         static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() ==
    // RDKit❗❌:             atom->getAtomicNum()) {
    // RDKit❗❌:       res = false;
    // RDKit❗❌:     } else if (descr == "AtomAnd") {
    // RDKit❗❌:       if ((*atom->getQuery()->beginChildren())->getDescription() ==
    // RDKit❗❌:           "AtomAtomicNum") {
    // RDKit❗❌:         res = false;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    match row.predicate() {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(n)) => *n != row.atomic_number(),
        QueryNode::And(children)
            if children.first().is_some_and(|first| {
                matches!(
                    first,
                    QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(_))
                )
            }) =>
        {
            false
        }
        _ => true,
    }
}

fn atom_list_query(row: &cosmolkit_model::QueryAtom) -> Option<(Vec<u8>, bool)> {
    // RDKit❗❌: template <typename T>
    // RDKit❗❌: bool _atomListQueryHelper(const T query, bool ignoreNegation) {
    // RDKit❗❌:   PRECONDITION(query, "no query");
    // RDKit❗❌:   if (!ignoreNegation && query->getNegation()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (query->getDescription() == "AtomAtomicNum" ||
    // RDKit❗❌:       query->getDescription() == "AtomType") {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (query->getDescription() == "AtomOr") {
    // RDKit❗❌:     for (const auto &child : boost::make_iterator_range(query->beginChildren(),
    // RDKit❗❌:                                                         query->endChildren())) {
    // RDKit❗❌:       if (!_atomListQueryHelper(child, ignoreNegation)) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   return false;
    // RDKit❗❌: }
    // RDKit❗❌: }  // namespace
    // RDKit❗❌: bool isAtomListQuery(const Atom *a) {
    // RDKit❗❌:   PRECONDITION(a, "bad atom");
    // RDKit❗❌:   if (!a->hasQuery()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (a->getQuery()->getDescription() == "AtomOr") {
    // RDKit❗❌:     for (const auto &child : boost::make_iterator_range(
    // RDKit❗❌:              a->getQuery()->beginChildren(), a->getQuery()->endChildren())) {
    // RDKit❗❌:       if (!_atomListQueryHelper(child, false)) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     return true;
    // RDKit❗❌:   } else if (a->getQuery()->getNegation() &&
    // RDKit❗❌:              _atomListQueryHelper(a->getQuery(), true)) {
    // RDKit❗❌:     // this was github #5930: negated list queries containing a single atom were
    // RDKit❗❌:     // being lost on output
    // RDKit❗❌:     return true;
    // RDKit❗❌:   } else if (a->getQuery()->getDescription() == "AtomAtomicNum" &&
    // RDKit❗❌:              static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal() !=
    // RDKit❗❌:                  a->getAtomicNum()) {
    // RDKit❗❌:     // when reading single-member atom lists from CTABs we end up with simple
    // RDKit❗❌:     // AtomAtomicNum queries where the atomic number of the atom itself is zero.
    // RDKit❗❌:     // Recognize this case.
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   return false;
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
    // RDKit❗❌:                           std::vector<int> &vals) {
    // RDKit❗❌:   // list queries are series of nested ors of AtomAtomicNum queries
    // RDKit❗❌:   PRECONDITION(q, "bad query");
    // RDKit❗❌:   auto descr = q->getDescription();
    // RDKit❗❌:   if (descr == "AtomOr") {
    // RDKit❗❌:     for (const auto &child :
    // RDKit❗❌:          boost::make_iterator_range(q->beginChildren(), q->endChildren())) {
    // RDKit❗❌:       auto descr = child->getDescription();
    // RDKit❗❌:       if (child->getNegation() ||
    // RDKit❗❌:           (descr != "AtomOr" && descr != "AtomAtomicNum" &&
    // RDKit❗❌:            descr != "AtomType")) {
    // RDKit❗❌:         throw ValueErrorException("bad query type1");
    // RDKit❗❌:       }
    // RDKit❗❌:       // we don't allow negation of any children of the query:
    // RDKit❗❌:       if (descr == "AtomOr") {
    // RDKit❗❌:         getAtomListQueryVals(child.get(), vals);
    // RDKit❗❌:       } else if (descr == "AtomAtomicNum") {
    // RDKit❗❌:         vals.push_back(static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal());
    // RDKit❗❌:       } else if (descr == "AtomType") {
    // RDKit❗❌:         auto v = static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal();
    // RDKit❗❌:         // aromatic AtomType queries add 1000 to the atomic number;
    // RDKit❗❌:         // correct for that:
    // RDKit❗❌:         if (v >= 1000) {
    // RDKit❗❌:           v -= 1000;
    // RDKit❗❌:         }
    // RDKit❗❌:         vals.push_back(v);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else if (descr == "AtomAtomicNum") {
    // RDKit❗❌:     vals.push_back(static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal());
    // RDKit❗❌:   } else if (descr == "AtomType") {
    // RDKit❗❌:     auto v = static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal();
    // RDKit❗❌:     // aromatic AtomType queries add 1000 to the atomic number;
    // RDKit❗❌:     // correct for that:
    // RDKit❗❌:     if (v >= 1000) {
    // RDKit❗❌:       v -= 1000;
    // RDKit❗❌:     }
    // RDKit❗❌:     vals.push_back(v);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     CHECK_INVARIANT(0, "bad query type");
    // RDKit❗❌:   }
    // RDKit❗❌: }
    fn values(node: &QueryNode<AtomQueryPredicate>, output: &mut Vec<u8>) -> bool {
        // RDKit❗❌: void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
        // RDKit❗❌:                           std::vector<int> &vals) {
        // RDKit❗❌:   // list queries are series of nested ors of AtomAtomicNum queries
        // RDKit❗❌:   PRECONDITION(q, "bad query");
        // RDKit❗❌:   auto descr = q->getDescription();
        // RDKit❗❌:   if (descr == "AtomOr") {
        // RDKit❗❌:     for (const auto &child :
        // RDKit❗❌:          boost::make_iterator_range(q->beginChildren(), q->endChildren())) {
        // RDKit❗❌:       auto descr = child->getDescription();
        // RDKit❗❌:       if (child->getNegation() ||
        // RDKit❗❌:           (descr != "AtomOr" && descr != "AtomAtomicNum" &&
        // RDKit❗❌:            descr != "AtomType")) {
        // RDKit❗❌:         throw ValueErrorException("bad query type1");
        // RDKit❗❌:       }
        // RDKit❗❌:       // we don't allow negation of any children of the query:
        // RDKit❗❌:       if (descr == "AtomOr") {
        // RDKit❗❌:         getAtomListQueryVals(child.get(), vals);
        // RDKit❗❌:       } else if (descr == "AtomAtomicNum") {
        // RDKit❗❌:         vals.push_back(static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal());
        // RDKit❗❌:       } else if (descr == "AtomType") {
        // RDKit❗❌:         auto v = static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal();
        // RDKit❗❌:         // aromatic AtomType queries add 1000 to the atomic number;
        // RDKit❗❌:         // correct for that:
        // RDKit❗❌:         if (v >= 1000) {
        // RDKit❗❌:           v -= 1000;
        // RDKit❗❌:         }
        // RDKit❗❌:         vals.push_back(v);
        // RDKit❗❌:       }
        // RDKit❗❌:     }
        // RDKit❗❌:   } else if (descr == "AtomAtomicNum") {
        // RDKit❗❌:     vals.push_back(static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal());
        // RDKit❗❌:   } else if (descr == "AtomType") {
        // RDKit❗❌:     auto v = static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal();
        // RDKit❗❌:     // aromatic AtomType queries add 1000 to the atomic number;
        // RDKit❗❌:     // correct for that:
        // RDKit❗❌:     if (v >= 1000) {
        // RDKit❗❌:       v -= 1000;
        // RDKit❗❌:     }
        // RDKit❗❌:     vals.push_back(v);
        // RDKit❗❌:   } else {
        // RDKit❗❌:     CHECK_INVARIANT(0, "bad query type");
        // RDKit❗❌:   }
        // RDKit❗❌: }
        match node {
            QueryNode::Predicate(
                AtomQueryPredicate::AtomicNumber(n)
                | AtomQueryPredicate::AtomType {
                    atomic_number: n, ..
                },
            ) => {
                output.push(*n);
                true
            }
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(numbers)) => {
                output.extend(numbers);
                true
            }
            QueryNode::Or(children) => children.iter().all(|c| values(c, output)),
            _ => false,
        }
    }
    let (node, negated) = match row.predicate() {
        QueryNode::Not(child) => (child.as_ref(), true),
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberNotIn(numbers)) => {
            return Some((numbers.clone(), true));
        }
        node => (node, false),
    };
    let is_list = negated
        || matches!(
            node,
            QueryNode::Or(_) | QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(_))
        )
        || matches!(node,QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(n)) if *n!=row.atomic_number());
    let mut result = Vec::new();
    (is_list && values(node, &mut result)).then_some((result, negated))
}

fn query_bond_topology(query: &QueryNode<BondQueryPredicate>) -> u32 {
    // RDKit❗❌: int getQueryBondTopology(const Bond *bond) {
    // RDKit❗❌:   PRECONDITION(bond, "no bond");
    // RDKit❗❌:   PRECONDITION(bond->hasQuery(), "no query");
    // RDKit❗❌:   int res = 0;
    // RDKit❗❌:   Bond::QUERYBOND_QUERY *qry = bond->getQuery();
    // RDKit❗❌:   // start by catching combined bond order + bond topology queries
    // RDKit❗❌:
    // RDKit❗❌:   if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
    // RDKit❗❌:       qry->endChildren() - qry->beginChildren() == 2) {
    // RDKit❗❌:     auto child1 = qry->beginChildren();
    // RDKit❗❌:     auto child2 = child1 + 1;
    // RDKit❗❌:     if (((*child1)->getDescription() == "BondInRing") !=
    // RDKit❗❌:         ((*child2)->getDescription() == "BondInRing")) {
    // RDKit❗❌:       if ((*child1)->getDescription() != "BondInRing") {
    // RDKit❗❌:         std::swap(child1, child2);
    // RDKit❗❌:       }
    // RDKit❗❌:       qry = child1->get();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (qry->getDescription() == "BondInRing") {
    // RDKit❗❌:     if (qry->getNegation()) {
    // RDKit❗❌:       res = 2;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res = 1;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    let query = if let QueryNode::And(children) = query {
        if children.len() == 2 {
            // Source tests the child's BondInRing description independently of negation.
            let is_ring = |node: &QueryNode<BondQueryPredicate>| match node {
                QueryNode::Predicate(BondQueryPredicate::IsInRing(_)) => true,
                QueryNode::Not(child) => matches!(
                    child.as_ref(),
                    QueryNode::Predicate(BondQueryPredicate::IsInRing(_))
                ),
                _ => false,
            };
            let a = is_ring(&children[0]);
            let b = is_ring(&children[1]);
            if a != b {
                if a { &children[0] } else { &children[1] }
            } else {
                query
            }
        } else {
            query
        }
    } else {
        query
    };
    match query {
        QueryNode::Predicate(BondQueryPredicate::IsInRing(true)) => 1,
        QueryNode::Predicate(BondQueryPredicate::IsInRing(false)) => 2,
        QueryNode::Not(child) => match child.as_ref() {
            QueryNode::Predicate(BondQueryPredicate::IsInRing(true)) => 2,
            QueryNode::Predicate(BondQueryPredicate::IsInRing(false)) => 1,
            _ => 0,
        },
        _ => 0,
    }
}

fn query_bond_type_code(row: Option<&cosmolkit_model::QueryBond>) -> Option<u32> {
    // RDKit❗❌: int getQueryBondSymbol(const Bond *bond) {
    // RDKit❗❌:   PRECONDITION(bond, "no bond");
    // RDKit❗❌:   PRECONDITION(bond->hasQuery(), "no query");
    // RDKit❗❌:   int res = 8;
    // RDKit❗❌:
    // RDKit❗❌:   Bond::QUERYBOND_QUERY *qry = bond->getQuery();
    // RDKit❗❌:   if (qry->getDescription() == "BondOrder" || getQueryBondTopology(bond)) {
    // RDKit❗❌:     // trap the simple bond-order query
    // RDKit❗❌:     res = 0;
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // start by catching combined bond order + bond topology queries
    // RDKit❗❌:     if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
    // RDKit❗❌:         qry->endChildren() - qry->beginChildren() == 2) {
    // RDKit❗❌:       auto child1 = qry->beginChildren();
    // RDKit❗❌:       auto child2 = child1 + 1;
    // RDKit❗❌:       if ((*child2)->getDescription() == "BondInRing") {
    // RDKit❗❌:         qry = child1->get();
    // RDKit❗❌:       } else if ((*child1)->getDescription() == "BondInRing") {
    // RDKit❗❌:         qry = child2->get();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (qry->getDescription() == "BondOr" && !qry->getNegation()) {
    // RDKit❗❌:       if (qry->endChildren() - qry->beginChildren() == 2) {
    // RDKit❗❌:         auto child1 = qry->beginChildren();
    // RDKit❗❌:         auto child2 = child1 + 1;
    // RDKit❗❌:         if ((*child1)->getDescription() == "BondOrder" &&
    // RDKit❗❌:             !(*child1)->getNegation() &&
    // RDKit❗❌:             (*child2)->getDescription() == "BondOrder" &&
    // RDKit❗❌:             !(*child2)->getNegation()) {
    // RDKit❗❌:           // ok, it's a bond query we have a chance of dealing with
    // RDKit❗❌:           int t1 = static_cast<BOND_EQUALS_QUERY *>(child1->get())->getVal();
    // RDKit❗❌:           int t2 = static_cast<BOND_EQUALS_QUERY *>(child2->get())->getVal();
    // RDKit❗❌:           if (t1 > t2) {
    // RDKit❗❌:             std::swap(t1, t2);
    // RDKit❗❌:           }
    // RDKit❗❌:           if (t1 == Bond::SINGLE && t2 == Bond::DOUBLE) {
    // RDKit❗❌:             res = 5;
    // RDKit❗❌:           } else if (t1 == Bond::SINGLE && t2 == Bond::AROMATIC) {
    // RDKit❗❌:             res = 6;
    // RDKit❗❌:           } else if (t1 == Bond::DOUBLE && t2 == Bond::AROMATIC) {
    // RDKit❗❌:             res = 7;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     } else if (qry->getDescription() == "SingleOrAromaticBond" &&
    // RDKit❗❌:                !qry->getNegation()) {
    // RDKit❗❌:       res = 6;
    // RDKit❗❌:     } else if (qry->getDescription() == "SingleOrDoubleBond" &&
    // RDKit❗❌:                !qry->getNegation()) {
    // RDKit❗❌:       res = 5;
    // RDKit❗❌:     } else if (qry->getDescription() == "DoubleOrAromaticBond" &&
    // RDKit❗❌:                !qry->getNegation()) {
    // RDKit❗❌:       res = 7;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    let row = row?;
    let query = row.predicate();
    if matches!(query, QueryNode::Predicate(BondQueryPredicate::Order(_)))
        || query_bond_topology(query) != 0
    {
        return Some(0);
    }
    let orders = match query {
        QueryNode::Predicate(BondQueryPredicate::OrderIn(orders)) => Some(orders.clone()),
        QueryNode::Or(children) if children.len() == 2 => children
            .iter()
            .map(|c| {
                if let QueryNode::Predicate(BondQueryPredicate::Order(order)) = c {
                    Some(*order)
                } else {
                    None
                }
            })
            .collect::<Option<Vec<_>>>(),
        _ => None,
    };
    Some(match orders.as_deref() {
        Some(orders)
            if orders.len() == 2
                && orders.contains(&BondOrder::Single)
                && orders.contains(&BondOrder::Double) =>
        {
            5
        }
        Some(orders)
            if orders.len() == 2
                && orders.contains(&BondOrder::Single)
                && orders.contains(&BondOrder::Aromatic) =>
        {
            6
        }
        Some(orders)
            if orders.len() == 2
                && orders.contains(&BondOrder::Double)
                && orders.contains(&BondOrder::Aromatic) =>
        {
            7
        }
        _ => 8,
    })
}

fn molfile_atom_symbol(atomic_number: u8) -> Result<&'static str, MolWriteError> {
    cosmolkit_core::rdkit_element_symbol(atomic_number)
        .map_err(|_| MolWriteError::Value("Atomic number not found".to_owned()))
}

// ── Wedge bond helpers (RDKit source-level port) ────────────────────────────

/// The source Python 2D export prepares an ephemeral 2D conformer even when
/// includeStereo is false; the original value is borrowed and never changed.
pub fn write_sdf_2d_with_params(
    data: MolWriteInput<'_>,
    params: &MolBlockWriteParams,
) -> Result<PropertyText, MolWriteError> {
    let MolWriteInput {
        topology,
        coordinates,
        properties,
        rings,
    } = data;
    let mut params = *params;
    params.force_2d = true;
    let selection = export_selection(&params, Some(CoordinateDimension::TwoD))?;
    if !params.include_coordinates
        || !matches!(params.coordinate_selection, MolCoordinateSelection::Auto)
        || !coordinates.conformers_2d.is_empty()
    {
        topology.validate()?;
        coordinates.validate_for_atom_count(topology.atoms.len())?;
        let input = MolWriteContext {
            topology: Cow::Borrowed(topology),
            coordinates: Cow::Borrowed(coordinates),
            properties,
            valence: None,
            rings,
            query: None,
        };
        let block = match params.format {
            SdfFormat::V2000 => mol_to_v2000_block_with_params(&input, selection, &params)?,
            SdfFormat::V3000 => mol_to_v3000_block_with_params(&input, selection, &params)?,
        };
        return Ok(append_sdf_record_fields(block, &input));
    }
    let conformer =
        cosmolkit_depict::compute_2d_coordinates(topology, properties, &Default::default())?;
    let mut prepared = coordinates.clone();
    prepared.record_source_conformer_append(CoordinateDimension::TwoD)?;
    prepared.conformers_2d.push(conformer);
    write_sdf_2d_with_params(
        MolWriteInput {
            topology,
            coordinates: &prepared,
            properties,
            rings,
        },
        &params,
    )
}

#[cfg(test)]
mod original_private_writer_conditions {
    use super::*;
    use cosmolkit_model::{AtomSpec, Element};
    #[test]
    fn non_default_valence_skips_periodic_table_lookup_when_implicit_hs_are_allowed() {
        let atom_id = AtomId::new(0);
        let topology = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(
                atom_id,
                AtomSpec::new(Element::C).with_formal_charge(9),
            )],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let molecule = MolWriteContext {
            topology: Cow::Borrowed(&topology),
            coordinates: Cow::Borrowed(&coordinates),
            properties: &properties,
            valence: None,
            rings: None,
            query: None,
        };
        let atom = &molecule.atoms()[atom_id.index()];
        assert!(!atom.no_implicit());
        let valence = molblock_valence_assignment(&molecule).unwrap();
        assert!(!has_non_default_valence(&molecule, atom, &valence).unwrap());
    }
    #[test]
    fn molfile_total_valence_field_tracks_rdkit_r_dummy_after_v3000_read() {
        let input = concat!(
            "\n",
            "     RDKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 R -0.750000 0.000000 0.000000 1 VAL=1\n",
            "M  V30 2 C 0.750000 -0.000000 0.000000 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );

        let crate::MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } = crate::read_mol_graph_record_detached_with_params(
            input,
            crate::SdfDataReadParams::default(),
        )
        .unwrap()
        .finish_mol_post(crate::MolPostParams::default())
        .unwrap()
        .mol_block
        else {
            panic!("original concrete R-dummy record must remain concrete");
        };
        let context = MolWriteContext {
            topology: Cow::Borrowed(&topology),
            coordinates: Cow::Borrowed(&coordinates),
            properties: &properties,
            valence: None,
            rings: None,
            query: None,
        };
        let molecule = &context;
        let atom = &molecule.atoms()[0];

        assert_eq!(atom.atomic_number(), 0);
        assert_eq!(
            atom.prop("dummyLabel"),
            Some(&cosmolkit_model::PropertyValue::String("R".into()))
        );
        assert_eq!(atom.prop("_MolFileRLabel"), None);
        assert_eq!(atom.atom_map(), Some(1));
        // The tagged MolBlockRecord::Molecule match above retains the original no-query condition.
        assert!(atom.no_implicit());
        assert_eq!(atom.explicit_hydrogens(), 0);
        assert_eq!(molecule.topology_block().adjacency.neighbors_of(0).len(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
        let valence = molblock_valence_assignment(molecule).unwrap();
        assert!(has_non_default_valence(molecule, atom, &valence).unwrap());
        assert_eq!(molfile_total_valence_field(molecule, atom).unwrap(), 1);
    }
}
