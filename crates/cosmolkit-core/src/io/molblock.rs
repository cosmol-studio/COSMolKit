//! MolBlock writer implementation.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::f64::consts::PI;

use crate::{
    Atom, AtomId, AtomQueryPredicate, Bond, BondDirection, BondId, BondOrder, BondQueryPredicate,
    BondStereo, ChiralTag, CoordinateDimension, QueryNode, SGroupBondRole, SGroupBracketStyle,
    SGroupConnection, SGroupData, StereoGroupKind, SubstanceGroup, SubstanceGroupKind,
};
use crate::{Molecule, RingInfo, UnsupportedFeatureError, find_sssr};

const MIN_V2000_COORD: f64 = -100_000.0;
const MAX_V2000_COORD: f64 = 1_000_000.0;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MolWriteError {
    #[error("MolBlock writing is not implemented")]
    NotImplemented,
    #[error("MolBlock writing subset is not supported: {0}")]
    UnsupportedSubset(&'static str),
    #[error("MolBlock writing failed: {0}")]
    Value(String),
    #[error(transparent)]
    Operation(#[from] crate::OperationError),
    #[error(transparent)]
    UnsupportedFeature(#[from] UnsupportedFeatureError),
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
}

impl Default for MolBlockWriteParams {
    fn default() -> Self {
        Self {
            format: SdfFormat::V2000,
            force_2d: false,
            include_stereo: true,
            kekulize: true,
            precision: 6,
        }
    }
}

enum CoordinateSelection {
    Auto,
    TwoD,
    ThreeD,
}

struct SelectedCoordinates {
    coords: Option<Vec<[f64; 3]>>,
    is_3d: bool,
    label: Option<&'static str>,
}

pub fn mol_to_v2000_2d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        force_2d: true,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::TwoD, &params)
}

pub fn mol_to_v2000_3d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::ThreeD, &params)
}

pub fn mol_to_v2000_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        ..Default::default()
    };
    mol_to_v2000_block_with_params(molecule, CoordinateSelection::Auto, &params)
}

pub fn mol_to_v3000_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::Auto, &params)
}

pub fn mol_to_v3000_2d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        force_2d: true,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::TwoD, &params)
}

pub fn mol_to_v3000_3d_block(molecule: &Molecule) -> Result<String, MolWriteError> {
    let params = MolBlockWriteParams {
        format: SdfFormat::V3000,
        ..Default::default()
    };
    mol_to_v3000_block_with_params(molecule, CoordinateSelection::ThreeD, &params)
}

pub fn mol_to_mol_block_with_params(
    molecule: &Molecule,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/FileWriters.h :: MolToMolBlock inline overloads
    // RDKit❗✔️: MolWriterParams params{includeStereo, kekulize, forceV3000};
    // RDKit❗✔️: return MolToMolBlock(mol, params, confId);
    // RDKit❗✔️: MolWriterParams v3KParams{params};
    // RDKit❗✔️: v3KParams.forceV3000 = true;
    // RDKit❗✔️: return MolToMolBlock(mol, v3KParams, confId);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/FileWriters.h :: MolToMolBlock inline overloads
    let selection = if params.force_2d {
        CoordinateSelection::TwoD
    } else {
        CoordinateSelection::Auto
    };
    match params.format {
        SdfFormat::V2000 => {
            if should_auto_upgrade_to_v3000(molecule) {
                mol_to_v3000_block_with_params(molecule, selection, params)
            } else {
                mol_to_v2000_block_with_params(molecule, selection, params)
            }
        }
        SdfFormat::V3000 => mol_to_v3000_block_with_params(molecule, selection, params),
    }
}

pub fn mol_to_2d_sdf_record(
    molecule: &Molecule,
    format: SdfFormat,
) -> Result<String, MolWriteError> {
    let block = match format {
        SdfFormat::V2000 => mol_to_v2000_2d_block(molecule)?,
        SdfFormat::V3000 => mol_to_v3000_2d_block(molecule)?,
    };
    Ok(append_sdf_record_fields(block, molecule))
}

pub fn mol_to_sdf_record_with_params(
    molecule: &Molecule,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    let block = mol_to_mol_block_with_params(molecule, params)?;
    Ok(append_sdf_record_fields(block, molecule))
}

pub fn mol_to_3d_sdf_record(
    molecule: &Molecule,
    format: SdfFormat,
) -> Result<String, MolWriteError> {
    let block = match format {
        SdfFormat::V2000 => mol_to_v2000_3d_block(molecule)?,
        SdfFormat::V3000 => mol_to_v3000_3d_block(molecule)?,
    };
    Ok(append_sdf_record_fields(block, molecule))
}

fn should_auto_upgrade_to_v3000(molecule: &Molecule) -> bool {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock auto-V3000 selection
    // RDKit✔️✔️: bool hasDative = false;
    // RDKit✔️✔️: for (const auto bond : tmol.bonds()) { if (bond->getBondType() == Bond::DATIVE) { hasDative = true; break; } }
    // RDKit✔️✔️: if (whichFormat == MolFileFormat::unspecified &&
    // RDKit✔️✔️:     (coordMagnitudeTooLargeForV2K || hasDative || nAtoms > 999 ||
    // RDKit✔️✔️:      nBonds > 999 || nSGroups > 999 || !tmol.getStereoGroups().empty())) {
    // RDKit✔️✔️:   isV3000 = true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock auto-V3000 selection
    molecule
        .bonds()
        .iter()
        .any(|bond| bond.order() == BondOrder::Dative)
        || molecule.num_atoms() > 999
        || molecule.num_bonds() > 999
        || molecule.substance_groups().len() > 999
        || !molecule.stereo_groups().is_empty()
}

/// Prepared molecule with aromatic-bond bookkeeping.
/// Before kekulization, the set of bonds that were aromatic is recorded so
/// crossed-bond stereo output can be suppressed for those bonds after
/// kekulization (matching RDKit's `prepareMol` → `aromaticBonds`
/// bookkeeping).
struct PreparedMol<'a> {
    molecule: Cow<'a, Molecule>,
    /// Indices (in the bond table) of bonds that were aromatic before
    /// kekulization. After kekulization, these should still be written
    /// as bond type 4 (aromatic) in the molfile output.
    aromatic_bonds: Vec<usize>,
    wedge_bonds: BTreeMap<usize, usize>,
    selected: SelectedCoordinates,
}

fn prepare_mol_for_writing<'a>(
    molecule: &'a Molecule,
    params: &MolBlockWriteParams,
    selection: CoordinateSelection,
) -> Result<PreparedMol<'a>, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: prepareMol
    // RDKit❗✔️: RWMol trwmol(mol);
    // RDKit❗✔️: if (params.kekulize && trwmol.getNumBonds()) { ... MolOps::Kekulize(trwmol); }
    // RDKit✔️✔️: if (params.includeStereo && !trwmol.getNumConformers()) {
    // RDKit✔️✔️:   RDDepict::compute2DCoords(trwmol);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: prepareMol
    // Record aromatic bonds before kekulization.
    let aromatic_bonds: Vec<usize> = molecule
        .bonds()
        .iter()
        .enumerate()
        .filter(|(_, bond)| bond.is_aromatic() || bond.order() == BondOrder::Aromatic)
        .map(|(idx, _)| idx)
        .collect();

    let mut mol = if params.kekulize {
        Cow::Owned(molecule.with_kekulized_bonds(true)?)
    } else {
        Cow::Borrowed(molecule)
    };

    let mut selected = select_coordinates(mol.as_ref(), selection)?;
    if params.include_stereo && selected.coords.is_none() {
        let generated = mol.as_ref().with_2d_coordinates().map_err(|source| {
            MolWriteError::Value(format!(
                "compute2DCoords failed during MolBlock write: {source}"
            ))
        })?;
        mol = Cow::Owned(generated);
        selected = select_coordinates(mol.as_ref(), CoordinateSelection::TwoD)?;
    }

    // outputMolToMolBlock() always calls pickBondsToWedge() before serializing
    // bond lines; this is not gated on includeStereo. The writer uses this map
    // transiently in GetMolFileBondLine instead of mutating bond endpoints.
    let ring_info = match mol.as_ref().derived_cache().rings.as_ref() {
        Some(ri) if ri.is_sssr_or_better() => ri.clone(),
        _ => find_sssr(mol.as_ref()).map_err(|e| MolWriteError::Value(e.to_string()))?,
    };
    let wedge_bonds = pick_bonds_to_wedge(mol.as_ref(), &ring_info);

    Ok(PreparedMol {
        molecule: mol,
        aromatic_bonds,
        wedge_bonds,
        selected,
    })
}

pub(crate) fn wedge_molecule_bonds_like_rdkit(
    molecule: &mut Molecule,
    coords: Option<&[[f64; 3]]>,
) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP CALL third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock
    // RDKit❗✔️: auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    let ring_info = match molecule.derived_cache().rings.as_ref() {
        Some(ri) if ri.is_sssr_or_better() => ri.clone(),
        _ => find_sssr(molecule).map_err(|e| MolWriteError::Value(e.to_string()))?,
    };
    let wedge_bonds = pick_bonds_to_wedge(molecule, &ring_info);

    if wedge_bonds.is_empty() {
        return Ok(());
    }

    let owned_coords;
    let coords = if let Some(coords) = coords {
        Some(coords)
    } else {
        owned_coords = select_coordinates(molecule, CoordinateSelection::Auto)?;
        owned_coords.coords.as_deref()
    };

    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: wedgeMolBonds
    // RDKit❗✔️: auto wedgeBonds = Chirality::pickBondsToWedge(mol, params, conf);
    // RDKit❗✔️: for (const auto &[wbi, wedgeInfo] : wedgeBonds) {
    // RDKit❗✔️:   auto bond = mol.getBondWithIdx(wbi);
    // RDKit❗✔️:   auto dir =
    // RDKit❗✔️:       detail::determineBondWedgeState(bond, wedgeInfo->getIdx(), conf);
    // RDKit❗✔️:   if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    // RDKit❗✔️:     bond->setBondDir(dir);
    // RDKit❗✔️:     if (static_cast<unsigned int>(wedgeInfo->getIdx()) !=
    // RDKit❗✔️:         bond->getBeginAtomIdx()) {
    // RDKit❗✔️:       auto tmp = bond->getBeginAtomIdx();
    // RDKit❗✔️:       bond->setBeginAtomIdx(bond->getEndAtomIdx());
    // RDKit❗✔️:       bond->setEndAtomIdx(tmp);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: wedgeMolBonds
    let mut dirs: Vec<(usize, BondDirection, usize)> = Vec::new();
    for (&bond_idx, &chiral_atom_idx) in &wedge_bonds {
        if let Some(bond) = molecule.bonds().get(bond_idx) {
            let dir = determine_bond_wedge_state(molecule, bond, chiral_atom_idx, coords)?;
            if matches!(dir, BondDirection::BeginWedge | BondDirection::BeginDash) {
                dirs.push((bond_idx, dir, chiral_atom_idx));
            }
        }
    }
    for (bond_idx, dir, chiral_atom_idx) in dirs {
        if let Some(bond) = molecule.topology_block_mut().bonds.get_mut(bond_idx) {
            bond.set_direction(dir);
            if bond.begin().index() != chiral_atom_idx {
                bond.set_endpoints(bond.end(), bond.begin());
            }
        }
    }
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: determineBondWedgeState
// RDKit❗✔️: Bond::BondDir determineBondWedgeState(const Bond *bond,
// RDKit❗✔️:                                       unsigned int fromAtomIdx,
// RDKit❗✔️:                                       const Conformer *conf) {
fn determine_bond_wedge_state(
    molecule: &Molecule,
    bond: &Bond,
    from_atom_idx: usize,
    coords: Option<&[[f64; 3]]>,
) -> Result<BondDirection, MolWriteError> {
    // RDKit❗✔️: auto res = bond->getBondDir();
    let mut res = bond.direction();
    // RDKit❗✔️: if (!conf) { return res; }
    let Some(coords) = coords else {
        return Ok(res);
    };
    if bond.order() != BondOrder::Single {
        return Ok(res);
    }

    // RDKit❗✔️: if (bond->getBeginAtom()->getIdx() == fromAtomIdx) {
    // RDKit❗✔️:   atom = bond->getBeginAtom();
    // RDKit❗✔️:   bondAtom = bond->getEndAtom();
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   atom = bond->getEndAtom();
    // RDKit❗✔️:   bondAtom = bond->getBeginAtom();
    // RDKit❗✔️: }
    let (atom_idx, bond_atom_idx) = if bond.begin().index() == from_atom_idx {
        (bond.begin().index(), bond.end().index())
    } else if bond.end().index() == from_atom_idx {
        (bond.end().index(), bond.begin().index())
    } else {
        return Ok(res);
    };
    let atom = &molecule.atoms()[atom_idx];

    // RDKit❗✔️: auto chiralType = atom->getChiralTag();
    // RDKit❗✔️: TEST_ASSERT(chiralType == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:             chiralType == Atom::CHI_TETRAHEDRAL_CCW);
    let chiral_type = atom.chiral_tag();
    if !matches!(
        chiral_type,
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    ) {
        return Ok(res);
    }

    let center_loc = *coords.get(atom_idx).ok_or_else(|| {
        MolWriteError::Value("wedge-state coordinate count does not match atom count".to_string())
    })?;
    let bond_atom_loc = *coords.get(bond_atom_idx).ok_or_else(|| {
        MolWriteError::Value("wedge-state coordinate count does not match atom count".to_string())
    })?;
    // RDKit❗✔️: centerLoc.z = 0.0;
    // RDKit❗✔️: tmpPt.z = 0.0;
    let center_loc = [center_loc[0], center_loc[1], 0.0];
    let bond_atom_loc = [bond_atom_loc[0], bond_atom_loc[1], 0.0];

    // RDKit❗✔️: refVect = centerLoc.directionVector(tmpPt);
    let Some(ref_vect) = normalized_direction(center_loc, bond_atom_loc) else {
        return Ok(res);
    };

    // RDKit❗✔️: neighborBondIndices.push_back(bond->getIdx());
    // RDKit❗✔️: neighborBondAngles.push_back(0.0);
    let mut neighbor_bond_indices = vec![bond.id().index()];
    let mut neighbor_bond_angles = vec![0.0];
    // RDKit❗✔️: for (const auto nbrBond : mol->atomBonds(atom)) {
    for (nbr_bond_idx, other_atom_idx) in atom_bonds(molecule, atom_idx) {
        if nbr_bond_idx == bond.id().index() {
            continue;
        }
        let other_loc = *coords.get(other_atom_idx).ok_or_else(|| {
            MolWriteError::Value(
                "wedge-state coordinate count does not match atom count".to_string(),
            )
        })?;
        let other_loc = [other_loc[0], other_loc[1], 0.0];
        // RDKit❗✔️: tmpVect = centerLoc.directionVector(tmpPt);
        let Some(tmp_vect) = normalized_direction(center_loc, other_loc) else {
            return Ok(res);
        };
        // RDKit✔️✔️: auto angle = refVect.signedAngleTo(tmpVect);
        // RDKit❗✔️: if (angle < 0.0) { angle += 2. * M_PI; }
        let mut angle = rdkit_point3d_signed_angle_to(ref_vect, tmp_vect);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        // RDKit❗✔️: while (angleIt != neighborBondAngles.end() && angle > (*angleIt)) {
        let insert_idx = neighbor_bond_angles
            .iter()
            .position(|existing| angle <= *existing)
            .unwrap_or(neighbor_bond_angles.len());
        neighbor_bond_angles.insert(insert_idx, angle);
        neighbor_bond_indices.insert(insert_idx, nbr_bond_idx);
    }

    // RDKit❗✔️: int nSwaps = atom->getPerturbationOrder(neighborBondIndices);
    let mut n_swaps =
        molfile_perturbation_order_from_bond_indices(molecule, atom_idx, &neighbor_bond_indices)?;

    // RDKit❗✔️: if (neighborBondAngles.size() == 3) {
    // RDKit❗✔️:   double angle1 = (*angleIt);
    // RDKit❗✔️:   double angle2 = (*angleIt);
    // RDKit❗✔️:   constexpr double angleTol = M_PI * 1.9 / 180.;
    // RDKit❗✔️:   if (angle2 - angle1 >= (M_PI - angleTol)) { nSwaps++; }
    // RDKit❗✔️: }
    if neighbor_bond_angles.len() == 3 {
        let angle1 = neighbor_bond_angles[1];
        let angle2 = neighbor_bond_angles[2];
        let angle_tol = PI * 1.9 / 180.0;
        if angle2 - angle1 >= (PI - angle_tol) {
            n_swaps = n_swaps.saturating_add(1);
        }
    }
    // RDKit❗✔️: if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:   if (nSwaps % 2 == 1) { res = Bond::BEGINDASH; }
    // RDKit❗✔️:   else { res = Bond::BEGINWEDGE; }
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   if (nSwaps % 2 == 1) { res = Bond::BEGINWEDGE; }
    // RDKit❗✔️:   else { res = Bond::BEGINDASH; }
    // RDKit❗✔️: }
    match chiral_type {
        ChiralTag::TetrahedralCcw => {
            res = if n_swaps % 2 == 1 {
                BondDirection::BeginDash
            } else {
                BondDirection::BeginWedge
            };
        }
        ChiralTag::TetrahedralCw => {
            res = if n_swaps % 2 == 1 {
                BondDirection::BeginWedge
            } else {
                BondDirection::BeginDash
            };
        }
        _ => {}
    }

    Ok(res)
}
// END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: determineBondWedgeState

fn normalized_direction(from: [f64; 3], to: [f64; 3]) -> Option<[f64; 3]> {
    let dx = to[0] - from[0];
    let dy = to[1] - from[1];
    let dz = to[2] - from[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len == 0.0 {
        return None;
    }
    Some([dx / len, dy / len, dz / len])
}

fn rdkit_point3d_signed_angle_to(reference: [f64; 3], other: [f64; 3]) -> f64 {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/Geometry/point.h :: Point3D::angleTo / signedAngleTo
    // RDKit✔️✔️: double angleTo(const Point3D &other) const {
    // RDKit✔️✔️:   double lsq = lengthSq() * other.lengthSq();
    // RDKit✔️✔️:   double dotProd = dotProduct(other);
    // RDKit✔️✔️:   dotProd /= sqrt(lsq);
    // RDKit✔️✔️:   if (dotProd <= -1.0) { return M_PI; }
    // RDKit✔️✔️:   if (dotProd >= 1.0) { return 0.0; }
    // RDKit✔️✔️:   return acos(dotProd);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double signedAngleTo(const Point3D &other) const {
    // RDKit✔️✔️:   double res = this->angleTo(other);
    // RDKit✔️✔️:   if ((this->x * other.y - this->y * other.x) < -zero_tolerance) {
    // RDKit✔️✔️:     res = 2.0 * M_PI - res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/Geometry/point.h :: Point3D::angleTo / signedAngleTo
    let lsq =
        (reference[0] * reference[0] + reference[1] * reference[1] + reference[2] * reference[2])
            * (other[0] * other[0] + other[1] * other[1] + other[2] * other[2]);
    let mut dot = reference[0] * other[0] + reference[1] * other[1] + reference[2] * other[2];
    dot /= lsq.sqrt();
    let mut res = if dot <= -1.0 {
        PI
    } else if dot >= 1.0 {
        0.0
    } else {
        dot.acos()
    };
    let cross_z = reference[0] * other[1] - reference[1] * other[0];
    if cross_z < -1.0e-16 {
        res = 2.0 * PI - res;
    }
    res
}

fn molfile_perturbation_order_from_bond_indices(
    molecule: &Molecule,
    atom_idx: usize,
    probe: &[usize],
) -> Result<u32, MolWriteError> {
    let reference = molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom_idx)
        .iter()
        .map(|neighbor| neighbor.bond.index())
        .collect::<Vec<_>>();
    if probe.len() != reference.len() {
        return Err(MolWriteError::Value(
            "Atom::getPerturbationOrder probe/reference length mismatch".to_string(),
        ));
    }
    let mut work = probe.to_vec();
    let mut swaps = 0_u32;
    for (idx, expected) in reference.iter().copied().enumerate() {
        if work[idx] == expected {
            continue;
        }
        let Some(found_idx) = work[idx..]
            .iter()
            .position(|bond_idx| *bond_idx == expected)
            .map(|offset| idx + offset)
        else {
            return Err(MolWriteError::Value(
                "Atom::getPerturbationOrder expected bond missing from probe order".to_string(),
            ));
        };
        work.swap(idx, found_idx);
        swaps = swaps.saturating_add(1);
    }
    Ok(swaps)
}

fn mol_to_v3000_block_with_params(
    molecule: &Molecule,
    selection: CoordinateSelection,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock / FileParserUtils::getV3000CTAB
    // RDKit❗✔️: if (isV3000) { ss << std::setw(3) << 0; ... ss << "999 V3000\n"; }
    // RDKit❗✔️: std::string res = "M  V30 BEGIN CTAB\n";
    // RDKit❗✔️: ss << "M  V30 COUNTS " << nAtoms << " " << nBonds << " " << nSGroups << " " << num3DConstraints << " " << chiralFlag << "\n";
    // RDKit❗✔️: res += "M  V30 BEGIN ATOM\n";
    // RDKit❗✔️: res += GetV3000MolFileAtomLine(*atomIt, conf, queryListAtoms, precision);
    // RDKit❗✔️: if (tmol.getNumBonds()) { res += "M  V30 BEGIN BOND\n"; ... res += GetV3000MolFileBondLine(...); }
    // RDKit❗✔️: if (nSGroups > 0) { res += "M  V30 BEGIN SGROUP\n"; ... }
    // RDKit❗✔️: appendEnhancedStereoGroups(res, tmol, wedgeBonds);
    // RDKit❗✔️: res += "M  V30 END CTAB\n";
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: outputMolToMolBlock / FileParserUtils::getV3000CTAB
    let prepared = prepare_mol_for_writing(molecule, params, selection)?;
    let molecule = prepared.molecule.as_ref();
    let aromatic_bonds = &prepared.aromatic_bonds;
    let wedge_bonds = &prepared.wedge_bonds;
    let selected = &prepared.selected;
    validate_v3000_writer_subset(molecule, params.include_stereo)?;
    let chiral_flag = molfile_chiral_flag(molecule)?;
    let generated_sgroups = v3000_generated_zbo_sgroups(molecule);

    let mut out = String::new();
    out.push_str(molecule.properties().name().unwrap_or_default());
    out.push('\n');
    out.push_str(&molfile_info_line(molecule, selected.label));
    out.push('\n');
    out.push_str(molecule.prop("_MolFileComments").unwrap_or_default());
    out.push('\n');
    out.push_str("  0  0  0  0  0  0  0  0  0  0999 V3000\n");
    out.push_str("M  V30 BEGIN CTAB\n");
    out.push_str(&format!(
        "M  V30 COUNTS {} {} {} 0 {}\n",
        molecule.num_atoms(),
        molecule.num_bonds(),
        molecule.substance_groups().len() + generated_sgroups.len(),
        chiral_flag
    ));
    let v3k_parity_flags: Vec<u32> = if selected.is_3d {
        if let Some(ref coords_3d) = selected.coords {
            let valence = molblock_valence_assignment(molecule).map_err(|source| {
                MolWriteError::Value(format!("V3000 atom parity assignment failed: {source}"))
            })?;
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
    out.push_str("M  V30 BEGIN ATOM\n");
    for atom in molecule.atoms() {
        let coord = selected
            .coords
            .as_ref()
            .and_then(|coords| coords.get(atom.id().index()).copied())
            .unwrap_or([0.0, 0.0, 0.0]);
        out.push_str(&v3000_atom_line(
            molecule,
            atom,
            coord,
            params.precision,
            &v3k_parity_flags,
        )?);
        out.push('\n');
    }
    out.push_str("M  V30 END ATOM\n");
    if molecule.num_bonds() != 0 {
        out.push_str("M  V30 BEGIN BOND\n");
        for bond in molecule.bonds() {
            out.push_str(&v3000_bond_line(
                molecule,
                bond,
                params.include_stereo,
                aromatic_bonds,
                wedge_bonds,
                selected.coords.as_deref(),
            )?);
            out.push('\n');
        }
        out.push_str("M  V30 END BOND\n");
    }
    append_v3000_sgroup_lines(&mut out, molecule, &generated_sgroups)?;
    if params.include_stereo {
        append_v3000_collection_lines(&mut out, molecule)?;
    }
    out.push_str("M  V30 END CTAB\n");
    out.push_str("M  END\n");
    Ok(out)
}

fn mol_to_v2000_block_with_params(
    molecule: &Molecule,
    selection: CoordinateSelection,
    params: &MolBlockWriteParams,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: MolToV2KMolBlock / outputMolToMolBlock
    // RDKit❗✔️: RWMol trwmol(mol);
    // RDKit❗✔️: prepareMol(trwmol, params, aromaticBonds);
    // RDKit❗✔️: return outputMolToMolBlock(trwmol, confId, MolFileFormat::V2000,
    // RDKit❗✔️:                              params.precision, aromaticBonds);
    // RDKit❗✔️: nAtoms = tmol.getNumAtoms();
    // RDKit❗✔️: nBonds = tmol.getNumBonds();
    // RDKit❗✔️: if (whichFormat == MolFileFormat::V2000 &&
    // RDKit❗✔️:     (nAtoms > 999 || nBonds > 999 || nSGroups > 999)) {
    // RDKit❗✔️:   throw ValueErrorException(
    // RDKit❗✔️:       "V2000 format does not support more than 999 atoms, bonds or SGroups.");
    // RDKit❗✔️: }
    // RDKit❗✔️: tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);
    // RDKit❗✔️: if (tmol.getPropIfPresent(common_properties::_Name, text)) { res += text; }
    // RDKit❗✔️: res += "\n";
    // RDKit❗✔️: // info
    // RDKit❗✔️: if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    // RDKit❗✔️:   res += text;
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   std::stringstream ss;
    // RDKit❗✔️:   ss << "  " << std::setw(8) << "RDKit";
    // RDKit❗✔️:   ss << std::setw(10) << "";
    // RDKit❗✔️:   if (conf) { if (conf->is3D()) { ss << "3D"; } else { ss << common_properties::TWOD; } }
    // RDKit❗✔️:   res += ss.str();
    // RDKit❗✔️: }
    // RDKit❗✔️: res += "\n";
    // RDKit❗✔️: if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) { res += text; }
    // RDKit❗✔️: res += "\n";
    // RDKit❗✔️: ss << std::setw(3) << nAtoms;
    // RDKit❗✔️: ss << std::setw(3) << nBonds;
    // RDKit❗✔️: ss << std::setw(3) << nLists;
    // RDKit❗✔️: ss << std::setw(3) << nSGroups;
    // RDKit❗✔️: ss << std::setw(3) << chiralFlag;
    // RDKit❗✔️: ss << std::setw(3) << nsText;
    // RDKit❗✔️: ss << std::setw(3) << nRxnComponents;
    // RDKit❗✔️: ss << std::setw(3) << nReactants;
    // RDKit❗✔️: ss << std::setw(3) << nProducts;
    // RDKit❗✔️: ss << std::setw(3) << nIntermediates;
    // RDKit❗✔️: ss << "999 V2000\n";
    // RDKit❗✔️: for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms(); atomIt != tmol.endAtoms(); ++atomIt) {
    // RDKit❗✔️:   res += GetMolFileAtomLine(*atomIt, conf, queryListAtoms);
    // RDKit❗✔️:   res += "\n";
    // RDKit❗✔️: }
    // RDKit❗✔️: auto wedgeBonds = Chirality::pickBondsToWedge(tmol, nullptr, conf);
    // RDKit❗✔️: for (const auto bond : tmol.bonds()) {
    // RDKit❗✔️:   res += GetMolFileBondLine(bond, wedgeBonds, conf, aromaticBonds[bond->getIdx()]);
    // RDKit❗✔️:   res += "\n";
    // RDKit❗✔️: }
    // RDKit❗✔️: res += GetMolFileChargeInfo(tmol);
    // RDKit❗✔️: res += GetMolFileRGroupInfo(tmol);
    // RDKit❗✔️: res += GetMolFileRGroupInfo(tmol);
    // RDKit❗✔️: res += GetMolFileQueryInfo(tmol, queryListAtoms);
    // RDKit❗✔️: res += GetMolFileAliasInfo(tmol);
    // RDKit❗✔️: res += GetMolFileZBOInfo(tmol);
    // RDKit❗✔️: res += GetMolFilePXAInfo(tmol);
    // RDKit❗✔️: res += GetMolFileSGroupInfo(tmol);
    // RDKit❗✔️: res += "M  END\n";
    // RDKit❗✔️: return res;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: MolToV2KMolBlock / outputMolToMolBlock
    let prepared = prepare_mol_for_writing(molecule, params, selection)?;
    let molecule = prepared.molecule.as_ref();
    let aromatic_bonds = &prepared.aromatic_bonds;
    let wedge_bonds = &prepared.wedge_bonds;
    let selected = &prepared.selected;
    validate_v2000_writer_subset(molecule, params.include_stereo)?;
    validate_v2000_coordinate_range(selected.coords.as_deref())?;
    let chiral_flag = molfile_chiral_flag(molecule)?;

    let mut out = String::new();
    out.push_str(molecule.properties().name().unwrap_or_default());
    out.push('\n');
    out.push_str(&molfile_info_line(molecule, selected.label));
    out.push('\n');
    out.push_str(molecule.prop("_MolFileComments").unwrap_or_default());
    out.push('\n');
    out.push_str(&format!(
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
    ));

    let parity_flags: Vec<u32> = if selected.is_3d {
        if let Some(ref coords_3d) = selected.coords {
            let valence = molblock_valence_assignment(molecule).map_err(|source| {
                MolWriteError::Value(format!("MolBlock atom parity assignment failed: {source}"))
            })?;
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
        out.push_str(&v2000_atom_line(atom, coord, molecule, &parity_flags)?);
        out.push('\n');
    }
    for bond in molecule.bonds() {
        out.push_str(&v2000_bond_line(
            molecule,
            bond,
            params.include_stereo,
            aromatic_bonds,
            wedge_bonds,
            selected.coords.as_deref(),
        )?);
        out.push('\n');
    }
    append_v2000_property_lines(&mut out, molecule)?;
    append_v2000_rgroup_lines(&mut out, molecule)?;
    append_v2000_value_lines(&mut out, molecule);
    append_v2000_alias_lines(&mut out, molecule);
    append_v2000_zbo_lines(&mut out, molecule);
    append_v2000_pxa_lines(&mut out, molecule);
    append_v2000_sgroup_lines(&mut out, molecule)?;
    // Write V2000 SSS query records for atoms and bonds with complex queries.
    // A single counter is shared across all SSS records in the CTfile.
    let mut sss_counter = 0u32;
    for atom in molecule.atoms() {
        out.push_str(&v2000_query_atom_sss_lines(atom, &mut sss_counter)?);
    }
    for bond in molecule.bonds() {
        out.push_str(&v2000_query_bond_sss_lines(bond, &mut sss_counter)?);
    }
    out.push_str("M  END\n");
    Ok(out)
}

/// Walk a query node and collect all leaf predicates, handling And/Or/Not.
fn flatten_query_predicates<T: Clone>(node: &QueryNode<T>) -> Vec<QueryFlattenNode<T>> {
    let mut result = Vec::new();
    flatten_query_predicates_inner(node, &mut result, false);
    result
}

fn flatten_query_predicates_inner<T: Clone>(
    node: &QueryNode<T>,
    result: &mut Vec<QueryFlattenNode<T>>,
    negated: bool,
) {
    match node {
        QueryNode::Predicate(p) => {
            result.push(QueryFlattenNode {
                predicate: p.clone(),
                negated,
            });
        }
        QueryNode::And(children) => {
            for child in children {
                flatten_query_predicates_inner(child, result, negated);
            }
        }
        QueryNode::Or(children) => {
            for child in children {
                flatten_query_predicates_inner(child, result, negated);
            }
        }
        QueryNode::Not(child) => {
            flatten_query_predicates_inner(child, result, !negated);
        }
    }
}

#[derive(Clone)]
struct QueryFlattenNode<T: Clone> {
    predicate: T,
    negated: bool,
}

/// Write V2000 SSS SAP records for a query atom. Returns the SAP lines.
fn v2000_query_atom_sss_lines(atom: &Atom, sss_counter: &mut u32) -> Result<String, MolWriteError> {
    let query = match atom.query() {
        Some(q) => q,
        None => return Ok(String::new()),
    };
    // If it's a simple atom-list query, the "L" symbol is already
    // written in the atom line; no SSS records needed.
    if v2000_atom_list_query(atom).is_some() {
        return Ok(String::new());
    }

    let predicates = flatten_query_predicates(query);
    if predicates.is_empty() {
        return Ok(String::new());
    }

    let mut out = String::new();
    let atom_idx = atom.id().index() + 1; // 1-based in MolBlock

    for flat in &predicates {
        match &flat.predicate {
            AtomQueryPredicate::RGroupLabel(id) => {
                // RDKit writes: M  ALS <nAtoms> <atomIdx> R<id>
                out.push_str(&format!("M  ALS{:>4}{:>4} R{}\\n", 1, atom_idx, id));
                continue;
            }
            AtomQueryPredicate::MolFileAlias(text) => {
                // RDKit writes: M  AAL <nAtoms> <atomIdx1> ... <aliasStr>
                out.push_str(&format!("M  AAL{:>4}{:>4} {}\\n", 1, atom_idx, text));
                continue;
            }
            AtomQueryPredicate::RecursiveSmarts(smarts) => {
                // RDKit writes: M  SMS <atomIdx> <SMARTS>
                out.push_str(&format!("M  SMS{:>4} {}\\n", atom_idx, smarts));
                continue;
            }
            _ => {}
        }
        let code = atom_query_predicate_to_sap_code(&flat.predicate);
        if code == 0 {
            continue;
        }
        *sss_counter += 1;
        let value = atom_query_predicate_to_sap_value(&flat.predicate);

        if flat.negated {
            out.push_str(&format!(
                "M  SAP{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter,
                atom_idx,
                code,
                -(value as i32)
            ));
        } else {
            out.push_str(&format!(
                "M  SAP{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter, atom_idx, code, value
            ));
        }
    }

    Ok(out)
}

/// Convert an AtomQueryPredicate to its V2000 SAP property code.
fn atom_query_predicate_to_sap_code(pred: &AtomQueryPredicate) -> u32 {
    match pred {
        AtomQueryPredicate::AtomicNumber(_) => 1,
        AtomQueryPredicate::AtomicNumberIn(_) => 1,
        AtomQueryPredicate::AtomicNumberNotIn(_) => 1,
        AtomQueryPredicate::FormalCharge(_) => 2,
        AtomQueryPredicate::Isotope(_) => 3,
        AtomQueryPredicate::ExplicitDegree(_) => 10,
        AtomQueryPredicate::ExplicitDegreeLessEqual(_) => 10,
        AtomQueryPredicate::NonHydrogenDegree(_) => 5,
        AtomQueryPredicate::ImplicitHydrogenCount(_) => 6,
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(_) => 6,
        AtomQueryPredicate::RingBondCount(_) => 9,
        AtomQueryPredicate::RingBondCountLessEqual(_) => 9,
        AtomQueryPredicate::RingBondCountNeedsScan => 9,
        AtomQueryPredicate::IsAromatic(_) => 12,
        AtomQueryPredicate::IsUnsaturated => 11,
        AtomQueryPredicate::RGroupLabel(_) => 20, // ALS (atom list symbol)
        AtomQueryPredicate::MolFileAlias(_) => 21, // AAL (atom alias)
        AtomQueryPredicate::RecursiveSmarts(_) => 22, // SMS (substructure SMARTS)
        AtomQueryPredicate::Any => 0,
        AtomQueryPredicate::UnsupportedFeature(_) => 0,
        // Phase A7 additions
        _ => 0,
    }
}

/// Convert an AtomQueryPredicate to its V2000 SAP property value.
fn atom_query_predicate_to_sap_value(pred: &AtomQueryPredicate) -> u32 {
    match pred {
        AtomQueryPredicate::AtomicNumber(z) => *z as u32,
        AtomQueryPredicate::AtomicNumberIn(_) => 0, // Handled by ALS
        AtomQueryPredicate::AtomicNumberNotIn(_) => 0,
        AtomQueryPredicate::FormalCharge(c) => (*c as i8) as u32,
        AtomQueryPredicate::Isotope(i) => *i as u32,
        AtomQueryPredicate::ExplicitDegree(d) => *d as u32,
        AtomQueryPredicate::ExplicitDegreeLessEqual(d) => *d as u32,
        AtomQueryPredicate::NonHydrogenDegree(d) => *d as u32,
        AtomQueryPredicate::ImplicitHydrogenCount(h) => *h as u32,
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(h) => *h as u32,
        AtomQueryPredicate::RingBondCount(c) => *c as u32,
        AtomQueryPredicate::RingBondCountLessEqual(c) => *c as u32,
        AtomQueryPredicate::RingBondCountNeedsScan => 1,
        AtomQueryPredicate::IsAromatic(true) => 1,
        AtomQueryPredicate::IsAromatic(false) => 0,
        AtomQueryPredicate::IsUnsaturated => 1,
        AtomQueryPredicate::RGroupLabel(id) => *id,
        AtomQueryPredicate::MolFileAlias(_) => 0, // String value, handled inline
        AtomQueryPredicate::RecursiveSmarts(_) => 0, // String value, handled inline
        AtomQueryPredicate::Any => 0,
        AtomQueryPredicate::UnsupportedFeature(_) => 0,
        // Phase A7 additions
        _ => 0,
    }
}

/// Write V2000 SSS SBT records for a query bond. Returns the SBT lines.
fn v2000_query_bond_sss_lines(bond: &Bond, sss_counter: &mut u32) -> Result<String, MolWriteError> {
    let query = match bond.query() {
        Some(q) => q,
        None => return Ok(String::new()),
    };
    // Simple bond query symbols and topology are already handled
    // by v2000_bond_query_symbol and v2000_bond_topology_code.
    // Only complex queries need SSS records.
    if v2000_bond_query_symbol(bond).is_some() {
        return Ok(String::new());
    }

    let predicates = flatten_query_predicates(query);
    if predicates.is_empty() {
        return Ok(String::new());
    }

    let mut out = String::new();
    let bond_idx = bond.id().index() + 1;

    for flat in &predicates {
        let (code, value) = bond_query_predicate_to_sbt_code(&flat.predicate);
        if code == 0 {
            continue;
        }
        *sss_counter += 1;
        if flat.negated {
            out.push_str(&format!(
                "M  SBT{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter,
                bond_idx,
                code,
                -(value as i32)
            ));
        } else {
            out.push_str(&format!(
                "M  SBT{:>4}{:>4}{:>4}{:>4}\n",
                *sss_counter, bond_idx, code, value
            ));
        }
    }

    Ok(out)
}

fn bond_query_predicate_to_sbt_code(pred: &BondQueryPredicate) -> (u32, u32) {
    match pred {
        BondQueryPredicate::Order(o) => match o {
            crate::BondOrder::Single => (1, 1),
            crate::BondOrder::Double => (1, 2),
            crate::BondOrder::Triple => (1, 3),
            crate::BondOrder::Quadruple => (1, 4),
            crate::BondOrder::Aromatic => (1, 5),
            _ => (0, 0),
        },
        BondQueryPredicate::IsAromatic(true) => (2, 1),
        BondQueryPredicate::IsAromatic(false) => (2, 0),
        BondQueryPredicate::IsInRing(true) => (3, 1),
        BondQueryPredicate::IsInRing(false) => (3, 0),
        BondQueryPredicate::Any => (0, 0), // Handled by symbol=8
        BondQueryPredicate::MolFileQueryCode(c) => (1, *c),
        _ => (0, 0),
    }
}

/// Validate that the molecule's query atoms/bonds can be written in
/// the V2000/V3000 format. Rejects only unsupported recursive SMARTS,
/// RGroupLabel, and MolFileAlias queries.
fn molfile_chiral_flag(molecule: &Molecule) -> Result<i32, MolWriteError> {
    molecule
        .prop("_MolFileChiralFlag")
        .map(|value| {
            value.trim().parse::<i32>().map_err(|_| {
                MolWriteError::Value(format!("invalid _MolFileChiralFlag value '{value}'"))
            })
        })
        .transpose()
        .map(|flag| flag.unwrap_or(0))
}

fn select_coordinates(
    molecule: &Molecule,
    selection: CoordinateSelection,
) -> Result<SelectedCoordinates, MolWriteError> {
    match selection {
        CoordinateSelection::TwoD => {
            let coords = molecule
                .coordinates_2d()
                .ok_or(MolWriteError::UnsupportedSubset(
                    "2D coordinates are required for V2000 2D output",
                ))?;
            if coords.len() != molecule.num_atoms() {
                return Err(MolWriteError::Value(
                    "2D coordinate count does not match atom count".to_string(),
                ));
            }
            Ok(SelectedCoordinates {
                coords: Some(
                    coords
                        .iter()
                        .map(|coord| [coord[0], coord[1], 0.0])
                        .collect(),
                ),
                is_3d: false,
                label: Some("2D"),
            })
        }
        CoordinateSelection::ThreeD => {
            let conformer =
                molecule
                    .conformers_3d()
                    .first()
                    .ok_or(MolWriteError::UnsupportedSubset(
                        "3D conformer coordinates are required for V2000 3D output",
                    ))?;
            if conformer.coordinates().len() != molecule.num_atoms() {
                return Err(MolWriteError::Value(
                    "3D coordinate count does not match atom count".to_string(),
                ));
            }
            Ok(SelectedCoordinates {
                coords: Some(conformer.coordinates().to_vec()),
                is_3d: true,
                label: Some("3D"),
            })
        }
        CoordinateSelection::Auto => {
            if matches!(
                molecule.source_coordinate_dim(),
                Some(CoordinateDimension::ThreeD)
            ) && !molecule.conformers_3d().is_empty()
            {
                return select_coordinates(molecule, CoordinateSelection::ThreeD);
            }
            if molecule.coordinates_2d().is_some() {
                return select_coordinates(molecule, CoordinateSelection::TwoD);
            }
            if !molecule.conformers_3d().is_empty() {
                return select_coordinates(molecule, CoordinateSelection::ThreeD);
            }
            Ok(SelectedCoordinates {
                coords: None,
                is_3d: false,
                label: None,
            })
        }
    }
}

fn validate_v2000_writer_subset(
    molecule: &Molecule,
    include_stereo: bool,
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
    for atom in molecule.atoms() {
        if atom.query().is_some() && v2000_atom_list_query(atom).is_none() {
            // Complex query atoms are written via SSS SAP records.
            // Only reject truly unsupported query types.
            if let Some(query) = atom.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let code = atom_query_predicate_to_sap_code(&flat.predicate);
                    // Accept: standard SAP codes (1-17) + special codes (20=RGroup, 21=Alias, 22=SMARTS)
                    if code == 0 || code > 22 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query atom predicate",
                        ));
                    }
                }
            }
        }
        if !atom.tracked_isotopic_hydrogens().is_empty() {
            return Err(MolWriteError::UnsupportedSubset(
                "collapsed hydrogen MolBlock writing is not ported",
            ));
        }
    }
    for bond in molecule.bonds() {
        if bond.query().is_some()
            && v2000_bond_query_symbol(bond).is_none()
            && v2000_bond_topology_code(bond).is_none()
        {
            // Complex query bonds are written via SSS SBT records.
            if let Some(query) = bond.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let (code, _) = bond_query_predicate_to_sbt_code(&flat.predicate);
                    if code == 0 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query bond predicate",
                        ));
                    }
                }
            }
        }
        if include_stereo {
            let empty_wedge_bonds = BTreeMap::new();
            v2000_bond_stereo_code(molecule, bond, &empty_wedge_bonds, None)?;
            if bond.unknown_stereo() {
                return Err(MolWriteError::UnsupportedSubset(
                    "bond stereochemistry MolBlock writing is not ported",
                ));
            }
        }
        v2000_bond_type_code(bond)?;
    }
    Ok(())
}

fn validate_v3000_writer_subset(
    molecule: &Molecule,
    include_stereo: bool,
) -> Result<(), MolWriteError> {
    if include_stereo {
        validate_v3000_stereo_groups(molecule)?;
    }
    for atom in molecule.atoms() {
        if atom.query().is_some() && v2000_atom_list_query(atom).is_none() {
            // Complex query atoms in V3000 are written via query
            // property annotations. Only reject truly unsupported types.
            if let Some(query) = atom.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let code = atom_query_predicate_to_sap_code(&flat.predicate);
                    if code == 0 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query atom V3000 predicate",
                        ));
                    }
                }
            }
        }
        if include_stereo && atom.unknown_stereo() {
            return Err(MolWriteError::UnsupportedSubset(
                "atom stereochemistry V3000 writing is not ported",
            ));
        }
        if !atom.tracked_isotopic_hydrogens().is_empty() {
            return Err(MolWriteError::UnsupportedSubset(
                "collapsed hydrogen V3000 writing is not ported",
            ));
        }
    }
    for bond in molecule.bonds() {
        if bond.query().is_some()
            && v2000_bond_query_symbol(bond).is_none()
            && v2000_bond_topology_code(bond).is_none()
        {
            // Complex query bonds in V3000 are supported via SBT equivalents.
            if let Some(query) = bond.query() {
                let predicates = flatten_query_predicates(query);
                for flat in &predicates {
                    let (code, _) = bond_query_predicate_to_sbt_code(&flat.predicate);
                    if code == 0 {
                        return Err(MolWriteError::UnsupportedSubset(
                            "unsupported query bond V3000 predicate",
                        ));
                    }
                }
            }
        }
        if include_stereo {
            let empty_wedge_bonds = BTreeMap::new();
            v3000_bond_cfg_code(molecule, bond, &empty_wedge_bonds, None)?;
            if bond.unknown_stereo() {
                return Err(MolWriteError::UnsupportedSubset(
                    "bond stereochemistry V3000 writing is not ported",
                ));
            }
        }
        v3000_bond_type_code(bond)?;
    }
    Ok(())
}

fn validate_v3000_stereo_groups(molecule: &Molecule) -> Result<(), MolWriteError> {
    // Atropisomer (bond-based) stereo groups are now supported via atom
    // collection from the bond endpoints. No separate validation needed.
    // RDKit additionally validates wedge bond consistency for atropisomer
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

fn molfile_info_line(molecule: &Molecule, label: Option<&'static str>) -> String {
    if let Some(label) = label {
        let mut line = format!("  {:>8}{:>10}", "COSMolKit", "");
        line.push_str(label);
        return line;
    }
    if let Some(info) = molecule.prop("_MolFileInfo") {
        return info.to_string();
    }
    if let Some(info) = molecule.prop("_MolFileInfoLine") {
        return info.to_string();
    }
    let mut line = format!("  {:>8}{:>10}", "COSMolKit", "");
    line
}

// Standard integer atomic weights (mass numbers of the most common isotope)
// indexed by atomic number (Z). Used to compute massDiff for molfile output:
//   massDiff = isotope - atomicWeight
const MOLFILE_ATOMIC_WEIGHT: [i32; 119] = [
    0, 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51, 52,
    55, 56, 59, 59, 64, 65, 70, 73, 75, 79, 80, 79, 85, 88, 89, 91, 93, 96, 98, 101, 103, 106, 108,
    112, 115, 119, 122, 128, 127, 132, 133, 138, 139, 140, 141, 144, 145, 152, 153, 157, 159, 162,
    165, 167, 169, 173, 175, 178, 181, 184, 187, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222,
    223, 226, 227, 232, 231, 238, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262, 263, 268,
    269, 275, 278, 281, 282, 285, 286, 289, 290, 293, 294, 294, 294, 294,
];

fn v2000_atom_line(
    atom: &Atom,
    coord: [f64; 3],
    molecule: &Molecule,
    parity_flags: &[u32],
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAtomProperties / GetMolFileAtomLine
    // RDKit❗✔️: totValence = 0;
    // RDKit❗✔️: atomMapNumber = 0;
    // RDKit❗✔️: parityFlag = 0;
    // RDKit❗✔️: x = y = z = 0.0;
    // RDKit❗✔️: if (conf) { const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx()); x = pos.x; y = pos.y; z = pos.z; ... }
    // RDKit✔️✔️: if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED && ...) { parityFlag = getAtomParityFlag(atom, conf); }
    // RDKit✔️✔️: if (hasNonDefaultValence(atom)) { ... }
    // RDKit❗✔️: snprintf(dest, 128,
    // RDKit❗✔️:          "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
    // RDKit❗✔️:          z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
    // RDKit❗✔️:          totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
    // RDKit❗✔️:          inversionFlag, exactChangeFlag);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAtomProperties / GetMolFileAtomLine
    let symbol = v2000_atom_symbol(atom, true)?;
    let atom_idx = atom.id().index();
    let parity_flag = parity_flags.get(atom_idx).copied().unwrap_or(0);
    let chg = 0;
    let tot_valence = molfile_total_valence_field(molecule, atom)?;
    let mass_diff = 0i32;
    let h_count = atom
        .prop("_MolFileHCount")
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    let stereo_care = atom
        .prop("_MolFileStereoCare")
        .and_then(|v| v.parse::<i32>().ok())
        .unwrap_or(0);
    Ok(format!(
        "{:>10.4}{:>10.4}{:>10.4} {}{:>2}{:>3}{:>3}{:>3}{:>3}{:>3}  0{:>3}{:>3}{:>3}{:>3}{:>3}",
        coord[0],
        coord[1],
        coord[2],
        symbol,
        mass_diff,
        chg,
        parity_flag,
        h_count,
        stereo_care,
        tot_valence,
        0, // rxnComponentType
        0, // rxnComponentNumber
        atom.atom_map().unwrap_or(0),
        0, // inversionFlag
        0, // exactChangeFlag
    ))
}

// RDKit✔️✔️: unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf) {
// RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
// RDKit✔️✔️:   PRECONDITION(conf, "bad conformer");
// RDKit✔️✔️:   if (!conf->is3D() ||
// RDKit✔️✔️:       !(atom->getDegree() >= 3 && atom->getTotalDegree() == 4)) {
// RDKit✔️✔️:     return 0;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   const ROMol &mol = atom->getOwningMol();
// RDKit✔️✔️:   RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
// RDKit✔️✔️:   std::vector<std::pair<unsigned int, RDGeom::Point3D>> vs;
// RDKit✔️✔️:   ... (neighbor loop with vector subtraction)
// RDKit✔️✔️:   std::sort(vs.begin(), vs.end(), Rankers::pairLess);
// RDKit✔️✔️:   double vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
// RDKit✔️✔️:   return (vol < 0) ? 2 : (vol > 0) ? 1 : 0;
// RDKit✔️✔️: }
// END RDKIT FUNCTION getAtomParityFlag
fn get_atom_parity_flag(
    molecule: &Molecule,
    atom: &Atom,
    coords_3d: &[[f64; 3]],
    valence: &crate::ValenceAssignment,
) -> u32 {
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
            // RDKit shifts H-atom indices by numAtoms so they sort after heavy atoms.
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
    if vol.abs() < 1e-10 {
        0
    } else if vol < 0.0 {
        2
    } else {
        1
    }
}

// RDKit✔️✔️: bool hasNonDefaultValence(const Atom *atom) {
// RDKit✔️✔️:   if (atom->getNumRadicalElectrons() != 0) {
// RDKit✔️✔️:     return true;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // for queries and atoms which don't have computed properties, the answer is
// RDKit✔️✔️:   // always no:
// RDKit✔️✔️:   if (atom->hasQuery() || atom->needsUpdatePropertyCache()) {
// RDKit✔️✔️:     return false;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (atom->getAtomicNum() == 1 ||
// RDKit✔️✔️:       SmilesWrite ::inOrganicSubset(atom->getAtomicNum())) {
// RDKit✔️✔️:     // for the ones we "know", we may have to specify the valence if it's
// RDKit✔️✔️:     // not the default value
// RDKit✔️✔️:     auto effAtomicNum = atom->getAtomicNum() - atom->getFormalCharge();
// RDKit✔️✔️:     return atom->getNoImplicit() &&
// RDKit✔️✔️:            (static_cast<int>(atom->getValence(Atom::ValenceType::EXPLICIT)) !=
// RDKit✔️✔️:             PeriodicTable::getTable()->getDefaultValence(effAtomicNum));
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return true;
// RDKit✔️✔️: }
// END RDKIT FUNCTION hasNonDefaultValence
fn has_non_default_valence(
    molecule: &Molecule,
    atom: &Atom,
    valence: &crate::ValenceAssignment,
) -> Result<bool, MolWriteError> {
    if atom.radical_electrons() != 0 {
        return Ok(true);
    }
    if atom.query().is_some() {
        return Ok(false);
    }
    if atom.atomic_number() == 1
        || crate::notation::smiles_write::in_organic_subset(atom.atomic_number()).unwrap_or(false)
    {
        let effective_atomic_num =
            i32::from(atom.atomic_number()) - i32::from(atom.formal_charge());
        let default_valence = u8::try_from(effective_atomic_num)
            .ok()
            .and_then(|atomic_number| {
                crate::chemistry::valence::rdkit_default_valence(atomic_number).ok()
            })
            .unwrap_or(-1);
        let explicit_valence = valence.explicit_valence[atom.id().index()];
        return Ok(atom.no_implicit() && explicit_valence != default_valence);
    }
    let _ = molecule;
    Ok(true)
}

fn molblock_valence_assignment(
    molecule: &Molecule,
) -> Result<crate::ValenceAssignment, MolWriteError> {
    match molecule.derived_cache().valence.as_ref() {
        Some(valence) => Ok(valence.clone()),
        None => crate::assign_valence_with_options(molecule, crate::ValenceModel::RdkitLike, false)
            .map_err(|source| {
                MolWriteError::Value(format!("MolBlock valence assignment failed: {source}"))
            }),
    }
}

fn molfile_total_valence_field(molecule: &Molecule, atom: &Atom) -> Result<u32, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAtomProperties valence branch
    // RDKit✔️✔️: totValence = 0;
    // RDKit✔️✔️: if (hasNonDefaultValence(atom)) {
    // RDKit✔️✔️:   if (atom->getTotalDegree() == 0) {
    // RDKit✔️✔️:     totValence = 15;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     totValence = atom->getTotalValence() % 15;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAtomProperties valence branch
    let valence = molblock_valence_assignment(molecule).map_err(|source| {
        MolWriteError::Value(format!(
            "MolBlock VAL field assignment failed for atom {}: {source}",
            atom.id().index()
        ))
    })?;
    if let Some(value) = atom
        .prop("molTotValence")
        .and_then(|value| value.parse::<i32>().ok())
    {
        return Ok(match value {
            -1 => 15,
            value if value > 0 => value as u32,
            _ => 0,
        });
    }
    if !has_non_default_valence(molecule, atom, &valence)? {
        return Ok(0);
    }
    // RDKit✔️✔️:   if (atom->getTotalDegree() == 0) {
    // RDKit✔️✔️:     totValence = 15;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     totValence = atom->getTotalValence() % 15;
    // RDKit✔️✔️:   }
    let total_degree = molfile_total_degree(molecule, atom, &valence);
    if total_degree == 0 {
        Ok(15)
    } else {
        let total_valence = valence.explicit_valence[atom.id().index()]
            + valence.implicit_hydrogens[atom.id().index()].max(0);
        Ok((total_valence as u32) % 15)
    }
}

fn molfile_total_degree(
    molecule: &Molecule,
    atom: &Atom,
    valence: &crate::ValenceAssignment,
) -> i32 {
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.id().index())
        .len() as i32
        + i32::from(atom.explicit_hydrogens())
        + valence.implicit_hydrogens[atom.id().index()].max(0)
}

fn should_be_crossed_bond_for_writer(
    molecule: &Molecule,
    bond: &Bond,
) -> Result<bool, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: shouldBeACrossedBond
    // RDKit✔️✔️: if (bond->getStereo() == Bond::STEREOANY) { ... return true; }
    // RDKit✔️✔️: if (bond->getStereo() != Bond::BondStereo::STEREONONE) { return false; }
    // RDKit✔️✔️: if (!Chirality::detail::isBondPotentialStereoBond(bond)) { return false; }
    // RDKit✔️✔️: if (bond->getBondDir() == Bond::EITHERDOUBLE) { return true; }
    // RDKit✔️✔️: if (beginAtom->getDegree() > 1 && endAtom->getDegree() > 1 &&
    // RDKit✔️✔️:     (beginAtom->getTotalValence() - beginAtom->getTotalDegree()) == 1 &&
    // RDKit✔️✔️:     (endAtom->getTotalValence() - endAtom->getTotalDegree()) == 1) {
    // RDKit❗✔️:   if (canBeStereoBond(bond)) { return true; }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return false;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: shouldBeACrossedBond
    if bond.order() != BondOrder::Double {
        return Ok(false);
    }

    let begin_idx = bond.begin().index();
    let end_idx = bond.end().index();
    let adjacency = &molecule.topology_block().adjacency;

    if bond.stereo() == BondStereo::Any {
        for neighbor in adjacency.neighbors_of(begin_idx) {
            let nbr_bond = &molecule.bonds()[neighbor.bond.index()];
            if nbr_bond.direction() == BondDirection::Unknown
                && nbr_bond.begin().index() == begin_idx
            {
                return Ok(false);
            }
        }
        for neighbor in adjacency.neighbors_of(end_idx) {
            let nbr_bond = &molecule.bonds()[neighbor.bond.index()];
            if nbr_bond.direction() == BondDirection::Unknown && nbr_bond.begin().index() == end_idx
            {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    if bond.stereo() != BondStereo::None {
        return Ok(false);
    }
    if !crate::stereo::is_bond_candidate_for_stereo(molecule, bond.id().index()) {
        return Ok(false);
    }
    if bond.direction() == BondDirection::EitherDouble {
        return Ok(true);
    }

    // RDKit✔️✔️: const auto beginAtom = bond->getBeginAtom();
    // RDKit✔️✔️: const auto endAtom = bond->getEndAtom();
    // RDKit✔️✔️: if (beginAtom->getDegree() > 1 && endAtom->getDegree() > 1 &&
    // RDKit✔️✔️:     (beginAtom->getTotalValence() - beginAtom->getTotalDegree()) == 1 &&
    // RDKit✔️✔️:     (endAtom->getTotalValence() - endAtom->getTotalDegree()) == 1) {
    let valence = molblock_valence_assignment(molecule)?;
    let has_one_unsaturation = |atom_idx: usize| {
        let atom = &molecule.atoms()[atom_idx];
        let degree = adjacency.neighbors_of(atom_idx).len();
        let total_valence =
            valence.explicit_valence[atom_idx] + valence.implicit_hydrogens[atom_idx].max(0);
        let total_degree = molfile_total_degree(molecule, atom, &valence);
        degree > 1 && total_valence - total_degree == 1
    };
    if has_one_unsaturation(begin_idx) && has_one_unsaturation(end_idx) {
        // RDKit✔️✔️:   // we only do this if each atom only has one unsaturation
        // RDKit✔️✔️:   // FIX: this is the fix for github #2649, but we will need to
        // RDKit✔️✔️:   // change it once we start handling allenes properly
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (canBeStereoBond(bond)) {
        // RDKit✔️✔️:     return true;  // crossed double bond
        // RDKit✔️✔️:   }
        if can_be_stereo_bond_for_writer(molecule, bond)? {
            return Ok(true);
        }
        // RDKit✔️✔️: }
    }

    // RDKit✔️✔️: return false;  // NOT crossed double bond
    Ok(false)
}

fn can_be_stereo_bond_for_writer(molecule: &Molecule, bond: &Bond) -> Result<bool, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: canBeStereoBond
    // RDKit✔️✔️: if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:     bond->getBondType() != Bond::BondType::AROMATIC) { return false; }
    // RDKit✔️✔️: for (const auto atom : {beginAtom, endAtom}) {
    // RDKit✔️✔️:   std::vector<int> nbrRanks;
    // RDKit✔️✔️:   for (auto nbrBond : bond->getOwningMol().atomBonds(atom)) {
    // RDKit✔️✔️:     if (nbrBond == bond) { continue; }
    // RDKit✔️✔️:     if (nbrBond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:       if (nbrBond->getBondDir() == Bond::ENDUPRIGHT ||
    // RDKit✔️✔️:           nbrBond->getBondDir() == Bond::ENDDOWNRIGHT) { return false; }
    // RDKit✔️✔️:       if (nbrBond->getBondDir() == Bond::BondDir::UNKNOWN &&
    // RDKit✔️✔️:           nbrBond->getBeginAtom() == atom) { return false; }
    // RDKit❗✔️:       if (rank >= 0 && rank already seen) { return false; } else { nbrRanks.push_back(rank); }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return true;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: canBeStereoBond
    if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic) {
        return Ok(false);
    }
    let adjacency = &molecule.topology_block().adjacency;
    let mut computed_cip_ranks: Option<Vec<u32>> = None;
    for atom_idx in [bond.begin().index(), bond.end().index()] {
        let mut nbr_ranks = Vec::new();
        for neighbor in adjacency.neighbors_of(atom_idx) {
            if neighbor.bond == bond.id() {
                continue;
            }
            let nbr_bond = &molecule.bonds()[neighbor.bond.index()];
            if nbr_bond.order() != BondOrder::Single {
                continue;
            }
            if matches!(
                nbr_bond.direction(),
                BondDirection::EndUpRight | BondDirection::EndDownRight
            ) {
                return Ok(false);
            }
            if nbr_bond.direction() == BondDirection::Unknown
                && nbr_bond.begin().index() == atom_idx
            {
                return Ok(false);
            }
            let other_atom = &molecule.atoms()[neighbor.atom_index];
            let rank = other_atom
                .prop("_ChiralAtomRank")
                .or_else(|| other_atom.prop("_CIPRank"))
                .and_then(|value| value.parse::<i32>().ok())
                .or_else(|| {
                    if computed_cip_ranks.is_none() {
                        computed_cip_ranks = Some(
                            crate::stereo::assign_atom_cip_ranks(molecule)
                                .map_err(|source| {
                                    MolWriteError::Value(format!(
                                        "MolBlock CIP rank assignment failed: {source}"
                                    ))
                                })
                                .ok()?,
                        );
                    }
                    computed_cip_ranks
                        .as_ref()
                        .and_then(|ranks| ranks.get(neighbor.atom_index))
                        .copied()
                        .map(|rank| rank as i32)
                })
                .unwrap_or(-1);
            if rank >= 0 {
                if nbr_ranks.contains(&rank) {
                    return Ok(false);
                }
                nbr_ranks.push(rank);
            }
        }
    }
    Ok(true)
}

fn v2000_bond_line(
    molecule: &Molecule,
    bond: &Bond,
    include_stereo: bool,
    aromatic_bonds: &[usize],
    wedge_bonds: &BTreeMap<usize, usize>,
    coords: Option<&[[f64; 3]]>,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: BondGetMolFileSymbol / GetMolFileBondLine
    // RDKit❗✔️: if (bond->hasQuery()) { res = getQueryBondSymbol(bond); }
    // RDKit❗✔️: switch (bond->getBondType()) { case Bond::SINGLE: ... case Bond::DATIVE: res = 9; ... }
    // RDKit❗✔️: RDKit::Chirality::GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode, reverse);
    // RDKit❗✔️: ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    // RDKit❗✔️: ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    // RDKit❗✔️: ss << std::setw(3) << symbol;
    // RDKit❗✔️: ss << " " << std::setw(2) << dirCode;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: BondGetMolFileSymbol / GetMolFileBondLine
    let (mut stereo_code, reverse) = v2000_bond_stereo_code(molecule, bond, wedge_bonds, coords)?;
    // RDKit✔️✔️: do not cross bonds which were aromatic before kekulization.
    if aromatic_bonds.contains(&bond.id().index()) && stereo_code == 3 {
        stereo_code = 0;
    }
    let type_code = v2000_bond_type_code(bond)?;
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
    if let Some(topology) = v2000_bond_topology_code(bond) {
        line.push_str(&format!(" {:>2} {:>2}", 0, topology));
    }
    Ok(line)
}

fn v2000_bond_stereo_code(
    molecule: &Molecule,
    bond: &Bond,
    wedge_bonds: &BTreeMap<usize, usize>,
    coords: Option<&[[f64; 3]]>,
) -> Result<(u32, bool), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: GetMolFileBondStereoInfo
    // RDKit❗✔️: reverse = false;
    // RDKit❗✔️: dir = Bond::NONE;
    // RDKit❗✔️: if (canHaveDirection(*bond)) { dir = determineBondWedgeState(...); }
    // RDKit❗✔️: else if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit❗✔️:   if (Chirality::shouldBeACrossedBond(bond)) { dir = Bond::BondDir::EITHERDOUBLE; }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: GetMolFileBondStereoInfo
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: BondGetDirCode
    // RDKit✔️✔️: switch (dir) {
    // RDKit✔️✔️:   case Bond::NONE: res = 0; break;
    // RDKit✔️✔️:   case Bond::BEGINWEDGE: res = 1; break;
    // RDKit✔️✔️:   case Bond::BEGINDASH: res = 6; break;
    // RDKit✔️✔️:   case Bond::UNKNOWN: res = 4; break;
    // RDKit✔️✔️:   case Bond::BondDir::EITHERDOUBLE: res = 3; break;
    // RDKit✔️✔️:   default: break;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: BondGetDirCode
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileBondLine stereo dirCode branch
    // RDKit❗✔️: BEGINWEDGE -> 1
    // RDKit❗✔️: EITHERDOUBLE double bond -> 3
    // RDKit❗✔️: UNKNOWN single bond -> 4
    // RDKit❗✔️: BEGINDASH -> 6
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileBondLine stereo dirCode branch
    let (dir, reverse) = molfile_bond_stereo_info(molecule, bond, wedge_bonds, coords)?;
    Ok((molfile_bond_dir_code(dir), reverse))
}

fn molfile_bond_stereo_info(
    molecule: &Molecule,
    bond: &Bond,
    wedge_bonds: &BTreeMap<usize, usize>,
    coords: Option<&[[f64; 3]]>,
) -> Result<(BondDirection, bool), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: GetMolFileBondStereoInfo
    // RDKit✔️✔️: reverse = false;
    // RDKit✔️✔️: dir = Bond::NONE;
    // RDKit✔️✔️: if (canHaveDirection(*bond)) {
    // RDKit✔️✔️:   dir = Chirality::detail::determineBondWedgeState(bond, wedgeBonds, conf);
    // RDKit✔️✔️:   if ((dir == Bond::BEGINDASH) ||
    // RDKit✔️✔️:       (dir == Bond::BEGINWEDGE || dir == Bond::UNKNOWN)) {
    // RDKit✔️✔️:     auto wbi = wedgeBonds.find(bond->getIdx());
    // RDKit✔️✔️:     if (wbi != wedgeBonds.end() && ... &&
    // RDKit✔️✔️:         static_cast<unsigned int>(wbi->second->getIdx()) !=
    // RDKit✔️✔️:             bond->getBeginAtomIdx()) {
    // RDKit✔️✔️:       reverse = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:   if (Chirality::shouldBeACrossedBond(bond)) {
    // RDKit✔️✔️:     dir = Bond::BondDir::EITHERDOUBLE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: GetMolFileBondStereoInfo
    let mut dir = BondDirection::None;
    let mut reverse = false;
    if can_have_direction_for_molfile(bond) {
        dir = if let Some(&from_atom_idx) = wedge_bonds.get(&bond.id().index()) {
            determine_bond_wedge_state(molecule, bond, from_atom_idx, coords)?
        } else {
            bond.direction()
        };
        if matches!(
            dir,
            BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown
        ) && wedge_bonds
            .get(&bond.id().index())
            .is_some_and(|&from_atom_idx| from_atom_idx != bond.begin().index())
        {
            reverse = true;
        }
    } else if bond.order() == BondOrder::Double
        && should_be_crossed_bond_for_writer(molecule, bond)?
    {
        dir = BondDirection::EitherDouble;
    }
    Ok((dir, reverse))
}

fn can_have_direction_for_molfile(bond: &Bond) -> bool {
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
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
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: BondGetMolFileSymbol
    // RDKit✔️✔️: int res = 0;
    // RDKit❗✔️: if (bond->hasQuery()) { res = getQueryBondSymbol(bond); }
    // RDKit✔️✔️: if (!res) {
    // RDKit✔️✔️:   switch (bond->getBondType()) {
    // RDKit✔️✔️:     case Bond::SINGLE: ... res = bond->getIsAromatic() ? 4 : 1; break;
    // RDKit✔️✔️:     case Bond::DOUBLE: ... res = bond->getIsAromatic() ? 4 : 2; break;
    // RDKit✔️✔️:     case Bond::TRIPLE: res = 3; break;
    // RDKit✔️✔️:     case Bond::AROMATIC: res = 4; break;
    // RDKit✔️✔️:     case Bond::ZERO: res = 1; break;
    // RDKit✔️✔️:     case Bond::DATIVE: res = 9; break;
    // RDKit✔️✔️:     default: break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return res;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: BondGetMolFileSymbol
    if let Some(code) = v2000_bond_query_symbol(bond) {
        return Ok(code);
    }
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
    molecule: &Molecule,
    atom: &Atom,
    coord: [f64; 3],
    precision: usize,
    parity_flags: &[u32],
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileAtomLine
    // RDKit❗✔️: ss << "M  V30 " << atom->getIdx() + 1;
    // RDKit❗✔️: std::string symbol = AtomGetMolFileSymbol(atom, false, queryListAtoms);
    // RDKit❗✔️: if (!isAtomListQuery(atom) || queryListAtoms[atom->getIdx()]) { ss << " " << symbol; } else { ... ss << " [" << symbols << "]"; }
    // RDKit❗✔️: ss << std::fixed; ss << std::setprecision(precision); ss << " " << x << " " << y << " " << z;
    // RDKit❗✔️: ss << " " << atomMapNumber;
    // RDKit✔️✔️: if (parityFlag != 0) { ss << " CFG=" << parityFlag; }
    // RDKit❗✔️: if (chg != 0) { ss << " CHG=" << chg; }
    // RDKit❗✔️: if (isotope != 0 && !isAtomRGroup(*atom)) { ss << " MASS=" << mass; }
    // RDKit❗✔️: if (nRadEs != 0 && atom->getTotalDegree() != 0) { ... ss << " RAD=" << nRadEs; }
    // RDKit❗✔️: if (totValence != 0) { if (totValence == 15) { ss << " VAL=-1"; } else { ss << " VAL=" << totValence; } }
    // RDKit❗✔️: if (symbol == "R#") { ... ss << " RGROUPS=(1 " << rLabel << ")"; }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileAtomLine
    let symbol = v3000_atom_symbol(atom)?;
    let parity_flag = parity_flags.get(atom.id().index()).copied().unwrap_or(0);
    let tot_valence = molfile_total_valence_field(molecule, atom)?;
    let mut out = format!("M  V30 {} {}", atom.id().index() + 1, symbol);
    out.push_str(&format!(
        " {0:.1$} {2:.1$} {3:.1$} {4}",
        coord[0],
        precision,
        coord[1],
        coord[2],
        atom.atom_map().unwrap_or(0)
    ));
    if parity_flag != 0 {
        out.push_str(&format!(" CFG={parity_flag}"));
    }
    if atom.formal_charge() != 0 {
        out.push_str(&format!(" CHG={}", atom.formal_charge()));
    }
    if let Some(isotope) = atom.isotope()
        && atom.prop("_MolFileRLabel").is_none()
    {
        out.push_str(&format!(" MASS={isotope}"));
    }
    let electrons = atom.radical_electrons();
    let valence = molblock_valence_assignment(molecule).map_err(|source| {
        MolWriteError::Value(format!(
            "V3000 MolBlock RAD field assignment failed for atom {}: {source}",
            atom.id().index()
        ))
    })?;
    if electrons != 0 && molfile_total_degree(molecule, atom, &valence) != 0 {
        let code = if electrons % 2 == 1 { 2 } else { 3 };
        out.push_str(&format!(" RAD={code}"));
    }
    if tot_valence != 0 {
        if tot_valence == 15 {
            out.push_str(" VAL=-1");
        } else {
            out.push_str(&format!(" VAL={tot_valence}"));
        }
    }
    if atom.prop("_MolFileRLabel").is_some() {
        let label = atom
            .prop("_MolFileRLabel")
            .and_then(|value| value.parse::<u32>().ok())
            .unwrap_or(1);
        out.push_str(&format!(" RGROUPS=(1 {label})"));
    }
    append_v3000_atom_int_prop(&mut out, atom, "molAttachOrder", "ATTCHORD");
    append_v3000_atom_int_prop(&mut out, atom, "molAttachPoint", "ATTCHPT");
    append_v3000_atom_int_prop(&mut out, atom, "molAtomSeqId", "SEQID");
    if let Some(value) = atom.prop("molAtomSeqName") {
        out.push_str(&format!(" SEQNAME={value}"));
    }
    append_v3000_atom_int_prop(&mut out, atom, "molRxnExactChange", "EXACHG");
    if let Some(value) = atom.mol_inversion_flag() {
        out.push_str(&format!(" INVRET={value}"));
    }
    append_v3000_atom_int_prop(&mut out, atom, "molStereoCare", "STBOX");
    append_v3000_atom_int_prop(&mut out, atom, "molSubstCount", "SUBST");
    append_v3000_atom_int_prop(&mut out, atom, "molRingBondCount", "RBCNT");
    if let Some(value) = atom.prop("molAtomClass") {
        out.push_str(&format!(" CLASS={value}"));
    }
    Ok(out)
}

fn append_v3000_atom_int_prop(out: &mut String, atom: &Atom, prop: &str, label: &str) {
    if let Some(value) = atom.prop(prop)
        && value != "0"
    {
        out.push_str(&format!(" {label}={value}"));
    }
}

fn v3000_atom_symbol(atom: &Atom) -> Result<String, MolWriteError> {
    if let Some((atomic_numbers, negated)) = v2000_atom_list_query(atom) {
        let mut value = String::new();
        if negated {
            value.push_str("\"NOT ");
        }
        value.push('[');
        for (idx, atomic_number) in atomic_numbers.iter().enumerate() {
            if idx != 0 {
                value.push(',');
            }
            value.push_str(molfile_atom_symbol(*atomic_number)?);
        }
        value.push(']');
        if negated {
            value.push('"');
        }
        return Ok(value);
    }
    Ok(v2000_atom_symbol(atom, false)?)
}

fn v3000_bond_line(
    molecule: &Molecule,
    bond: &Bond,
    include_stereo: bool,
    aromatic_bonds: &[usize],
    wedge_bonds: &BTreeMap<usize, usize>,
    coords: Option<&[[f64; 3]]>,
) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine
    // RDKit❗✔️: ss << "M  V30 " << bond->getIdx() + 1;
    // RDKit❗✔️: ss << " " << GetV3000BondCode(bond);
    // RDKit❗✔️: ss << " " << bond->getBeginAtomIdx() + 1 << " " << bond->getEndAtomIdx() + 1;
    // RDKit❗✔️: if (includeStereo) { ... ss << " CFG=" << cfg; }
    // RDKit❗✔️: if (bond->hasQuery()) { int topol = getQueryBondTopology(bond); if (topol) { ss << " TOPO=" << topol; } }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::molReactStatus, iprop) && iprop) { ss << " RXCTR=" << iprop; }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::molStereoCare, sprop) && sprop != "0") { ss << " STBOX=" << sprop; }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, sprop) && sprop != "0") { ss << " ENDPTS=" << sprop; }
    // RDKit❗✔️: if (bond->getPropIfPresent(common_properties::_MolFileBondAttach, sprop) && sprop != "0") { ss << " ATTACH=" << sprop; }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine
    let (mut cfg, reverse) = v3000_bond_cfg_code(molecule, bond, wedge_bonds, coords)?;
    let (begin_idx, end_idx) = if reverse {
        (bond.end().index(), bond.begin().index())
    } else {
        (bond.begin().index(), bond.end().index())
    };
    let type_code = v3000_bond_type_code(bond)?;
    let mut out = format!(
        "M  V30 {} {} {} {}",
        bond.id().index() + 1,
        type_code,
        begin_idx + 1,
        end_idx + 1
    );
    if aromatic_bonds.contains(&bond.id().index())
        && cfg == Some(2)
        && matches!(
            molfile_bond_stereo_info(molecule, bond, wedge_bonds, coords)?.0,
            BondDirection::EitherDouble
        )
    {
        cfg = None;
    }
    if let Some(cfg) = cfg {
        out.push_str(&format!(" CFG={cfg}"));
    }
    if let Some(topology) = v2000_bond_topology_code(bond) {
        out.push_str(&format!(" TOPO={topology}"));
    }
    append_v3000_bond_prop(&mut out, bond, "molReactStatus", "RXCTR");
    append_v3000_bond_prop(&mut out, bond, "molStereoCare", "STBOX");
    append_v3000_bond_prop(&mut out, bond, "_MolFileBondEndPts", "ENDPTS");
    append_v3000_bond_prop(&mut out, bond, "_MolFileBondAttach", "ATTACH");
    Ok(out)
}

fn v3000_bond_cfg_code(
    molecule: &Molecule,
    bond: &Bond,
    wedge_bonds: &BTreeMap<usize, usize>,
    coords: Option<&[[f64; 3]]>,
) -> Result<(Option<u32>, bool), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: GetMolFileBondStereoInfo
    // RDKit❗✔️: reverse = false;
    // RDKit❗✔️: dir = Bond::NONE;
    // RDKit❗✔️: if (canHaveDirection(*bond)) { dir = determineBondWedgeState(...); }
    // RDKit❗✔️: else if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit❗✔️:   if (Chirality::shouldBeACrossedBond(bond)) { dir = Bond::BondDir::EITHERDOUBLE; }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: GetMolFileBondStereoInfo
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: BondGetDirCode
    // RDKit✔️✔️: switch (dir) {
    // RDKit✔️✔️:   case Bond::NONE: res = 0; break;
    // RDKit✔️✔️:   case Bond::BEGINWEDGE: res = 1; break;
    // RDKit✔️✔️:   case Bond::BEGINDASH: res = 6; break;
    // RDKit✔️✔️:   case Bond::UNKNOWN: res = 4; break;
    // RDKit✔️✔️:   case Bond::BondDir::EITHERDOUBLE: res = 3; break;
    // RDKit✔️✔️:   default: break;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Chirality.cpp :: BondGetDirCode
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine stereo CFG branch
    // RDKit❗✔️: BEGINWEDGE -> CFG=1
    // RDKit❗✔️: UNKNOWN single bond -> CFG=2
    // RDKit❗✔️: EITHERDOUBLE double bond -> CFG=2
    // RDKit❗✔️: BEGINDASH -> CFG=3
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000MolFileBondLine stereo CFG branch
    let (dir, reverse) = molfile_bond_stereo_info(molecule, bond, wedge_bonds, coords)?;
    let cfg = match molfile_bond_dir_code(dir) {
        0 => None,
        1 => Some(1),
        3 | 4 => Some(2),
        6 => Some(3),
        _ => None,
    };
    Ok((cfg, reverse))
}

fn append_v3000_bond_prop(out: &mut String, bond: &Bond, prop: &str, label: &str) {
    if let Some(value) = bond.prop(prop)
        && value != "0"
    {
        out.push_str(&format!(" {label}={value}"));
    }
}

fn v3000_bond_type_code(bond: &Bond) -> Result<u32, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000BondCode
    // RDKit✔️✔️: int res = 0;
    // RDKit❗✔️: if (bond->hasQuery()) { res = getQueryBondSymbol(bond); }
    // RDKit✔️✔️: if (!res) {
    // RDKit✔️✔️:   switch (bond->getBondType()) {
    // RDKit✔️✔️:     case Bond::SINGLE: ... res = bond->getIsAromatic() ? 4 : 1; break;
    // RDKit✔️✔️:     case Bond::DOUBLE: ... res = bond->getIsAromatic() ? 4 : 2; break;
    // RDKit✔️✔️:     case Bond::TRIPLE: res = 3; break;
    // RDKit✔️✔️:     case Bond::AROMATIC: res = 4; break;
    // RDKit✔️✔️:     case Bond::DATIVE: res = 9; break;
    // RDKit✔️✔️:     case Bond::HYDROGEN: res = 10; break;
    // RDKit✔️✔️:     case Bond::ZERO: res = 1; break;
    // RDKit✔️✔️:     default: res = 0; break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return res;
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetV3000BondCode
    if let Some(code) = v2000_bond_query_symbol(bond) {
        return Ok(code);
    }
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
    out: &mut String,
    molecule: &Molecule,
    generated_sgroups: &[SubstanceGroup],
) -> Result<(), MolWriteError> {
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() && generated_sgroups.is_empty() {
        return Ok(());
    }
    out.push_str("M  V30 BEGIN SGROUP\n");
    for (idx, sgroup) in sgroups.iter().enumerate() {
        out.push_str(&v3000_sgroup_lines(idx + 1, sgroup)?);
    }
    let offset = sgroups.len();
    for (idx, sgroup) in generated_sgroups.iter().enumerate() {
        out.push_str(&v3000_sgroup_lines(offset + idx + 1, sgroup)?);
    }
    out.push_str("M  V30 END SGROUP\n");
    Ok(())
}

fn v3000_generated_zbo_sgroups(molecule: &Molecule) -> Vec<SubstanceGroup> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: createZBOSubstanceGroups
    // RDKit❗✔️: SubstanceGroup bsg(&mol, "DAT"); bsg.setProp("FIELDNAME", "ZBO");
    // RDKit❗✔️: for (const auto bond : mol.bonds()) { if (bond->getBondType() == Bond::ZERO) { bsg.addBondWithIdx(...); atomsAffected[...] = 1; } }
    // RDKit❗✔️: if (atomsAffected.any()) { for affected atoms { bsg.addAtomWithIdx(i); ... asgText += getTotalNumHs(); zsgText += getFormalCharge(); } }
    // RDKit❗✔️: addSubstanceGroup(mol, bsg); asg.setProp("DATAFIELDS", aDataFields); addSubstanceGroup(mol, asg); zsg.setProp("DATAFIELDS", zDataFields); addSubstanceGroup(mol, zsg);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: createZBOSubstanceGroups
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
        .filter_map(|(idx, is_affected)| is_affected.then_some(crate::AtomId::new(idx)))
        .collect::<Vec<_>>();
    let hydrogens = atoms
        .iter()
        .map(|atom| {
            molecule.atoms()[atom.index()]
                .explicit_hydrogens()
                .to_string()
        })
        .collect::<Vec<_>>()
        .join(";");
    let charges = atoms
        .iter()
        .map(|atom| molecule.atoms()[atom.index()].formal_charge().to_string())
        .collect::<Vec<_>>()
        .join(";");
    vec![
        SubstanceGroup::new(crate::SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(atoms.clone())
            .with_bonds(zero_bonds)
            .with_data(SGroupData {
                field_name: Some("ZBO".to_string()),
                ..SGroupData::default()
            }),
        SubstanceGroup::new(crate::SubstanceGroupId::new(1), SubstanceGroupKind::Data)
            .with_atoms(atoms.clone())
            .with_data(SGroupData {
                field_name: Some("HYD".to_string()),
                values: vec![hydrogens],
                ..SGroupData::default()
            }),
        SubstanceGroup::new(crate::SubstanceGroupId::new(2), SubstanceGroupKind::Data)
            .with_atoms(atoms)
            .with_data(SGroupData {
                field_name: Some("ZCH".to_string()),
                values: vec![charges],
                ..SGroupData::default()
            }),
    ]
}

fn v3000_sgroup_lines(idx: usize, sgroup: &SubstanceGroup) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetV3000MolFileSGroupLines
    // RDKit❗✔️: std::string currLine = (boost::format("M  V30 %d %s %d") % idx % sgroup.getProp<std::string>("TYPE") % id).str();
    // RDKit❗✔️: addBlockToSGroupString(BuildV3000IdxVectorDataBlock("ATOMS", sgroup.getAtoms()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(BuildV3000BondsBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(BuildV3000IdxVectorDataBlock("PATOMS", sgroup.getParentAtoms()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("SUBTYPE", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("MULT", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("CONNECT", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000ParentBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CompNoBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("LABEL", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000BracketBlock(sgroup.getBrackets()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("ESTATE", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000CStateBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDNAME", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDINFO", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("FIELDDISP", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("QUERYTYPE", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("QUERYOP", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000FieldDataBlock(sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("CLASS", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000AttachPointBlock(sgroup.getAttachPoints()), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("BRKTYP", sgroup), currLine, os);
    // RDKit❗✔️: addBlockToSGroupString(FormatV3000StringPropertyBlock("SEQID", sgroup), currLine, os);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetV3000MolFileSGroupLines
    let mut os = String::new();
    let mut current = format!(
        "M  V30 {} {} {}",
        idx,
        sgroup_kind_code(sgroup)?,
        sgroup.external_id().unwrap_or(0)
    );
    add_v3000_sgroup_block(
        v3000_index_vector_block("ATOMS", sgroup.atoms().iter().map(|atom| atom.index() + 1)),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(v3000_sgroup_bonds_block(sgroup), &mut current, &mut os);
    add_v3000_sgroup_block(
        v3000_index_vector_block(
            "PATOMS",
            sgroup.parent_atoms().iter().map(|atom| atom.index() + 1),
        ),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("SUBTYPE", sgroup.subtype()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("MULT", sgroup.props().get("MULT").map(String::as_str)),
        &mut current,
        &mut os,
    );
    let connect = sgroup
        .connection()
        .map(sgroup_connection_code)
        .transpose()?;
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("CONNECT", connect),
        &mut current,
        &mut os,
    );
    if let Some(parent) = sgroup.parent() {
        add_v3000_sgroup_block(
            format!(" PARENT={}", parent.index() + 1),
            &mut current,
            &mut os,
        );
    }
    if let Some(compno) = sgroup.component_number() {
        add_v3000_sgroup_block(format!(" COMPNO={compno}"), &mut current, &mut os);
    }
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("LABEL", sgroup.label()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(v3000_sgroup_bracket_block(sgroup), &mut current, &mut os);
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("ESTATE", sgroup.expansion_state()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(v3000_sgroup_cstate_block(sgroup), &mut current, &mut os);
    if let Some(data) = sgroup.data() {
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("FIELDNAME", data.field_name.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("FIELDINFO", data.field_info.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("FIELDDISP", data.field_display.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("QUERYTYPE", data.query_type.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_string_property_block("QUERYOP", data.query_op.as_deref()),
            &mut current,
            &mut os,
        );
        add_v3000_sgroup_block(
            v3000_sgroup_field_data_block(data.values.iter().map(String::as_str)),
            &mut current,
            &mut os,
        );
    }
    add_v3000_sgroup_block(
        v3000_sgroup_field_data_block(sgroup.data_fields().iter().map(String::as_str)),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("CLASS", sgroup.class()),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_attach_point_block(sgroup),
        &mut current,
        &mut os,
    );
    let bracket_style = v3000_sgroup_bracket_style_value(sgroup.bracket_style())?;
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block("BRKTYP", bracket_style),
        &mut current,
        &mut os,
    );
    add_v3000_sgroup_block(
        v3000_sgroup_string_property_block(
            "SEQID",
            sgroup.props().get("SEQID").map(String::as_str),
        ),
        &mut current,
        &mut os,
    );
    if !current.is_empty() && current != "M  V30" {
        os.push_str(&current);
        os.push('\n');
    }
    Ok(os)
}

fn v3000_index_vector_block(label: &str, values: impl IntoIterator<Item = usize>) -> String {
    let values = values.into_iter().collect::<Vec<_>>();
    if values.is_empty() {
        return String::new();
    }
    let mut out = format!(" {label}=({}", values.len());
    for value in values {
        out.push_str(&format!(" {value}"));
    }
    out.push(')');
    out
}

fn v3000_sgroup_bonds_block(sgroup: &SubstanceGroup) -> String {
    let crossing = sgroup
        .bonds()
        .iter()
        .filter(|bond| sgroup.bond_role(**bond) == SGroupBondRole::Crossing)
        .map(|bond| bond.index() + 1);
    let contained = sgroup
        .bonds()
        .iter()
        .filter(|bond| sgroup.bond_role(**bond) == SGroupBondRole::Contained)
        .map(|bond| bond.index() + 1);
    let mut out = v3000_index_vector_block("XBONDS", crossing);
    out.push_str(&v3000_index_vector_block("CBONDS", contained));
    out
}

fn v3000_sgroup_string_property_block(label: &str, value: Option<&str>) -> String {
    let Some(value) = value else {
        return String::new();
    };
    if value.is_empty() {
        return String::new();
    }
    let mut out = format!(" {label}=");
    out.push_str(&v3000_quote_string_property(value));
    out
}

fn v3000_quote_string_property(value: &str) -> String {
    let needs_quotes = value.contains(' ') || value.contains('"') || value.contains('(');
    let mut out = String::new();
    if needs_quotes {
        out.push('"');
    }
    for ch in value.chars() {
        out.push(ch);
        if ch == '"' {
            out.push(ch);
        }
    }
    if needs_quotes {
        out.push('"');
    }
    out
}

fn v3000_sgroup_bracket_block(sgroup: &SubstanceGroup) -> String {
    let Some(display) = sgroup.display() else {
        return String::new();
    };
    let mut out = String::new();
    for bracket in &display.brackets {
        out.push_str(" BRKXYZ=(9");
        out.push_str(&format!(
            " {} {} 0 {} {} 0 0 0 0)",
            v3000_double_field(bracket.p1[0]),
            v3000_double_field(bracket.p1[1]),
            v3000_double_field(bracket.p2[0]),
            v3000_double_field(bracket.p2[1])
        ));
    }
    out
}

fn v3000_sgroup_cstate_block(sgroup: &SubstanceGroup) -> String {
    let mut out = String::new();
    for cstate in sgroup.cstates() {
        out.push_str(" CSTATE=(");
        if matches!(sgroup.kind(), SubstanceGroupKind::Superatom) {
            out.push_str(&format!(
                "4 {} {} {} 0",
                cstate.bond.index() + 1,
                v3000_double_field(cstate.vector[0]),
                v3000_double_field(cstate.vector[1])
            ));
        } else {
            out.push_str(&format!("1 {}", cstate.bond.index() + 1));
        }
        out.push(')');
    }
    out
}

fn v3000_sgroup_field_data_block<'a>(values: impl IntoIterator<Item = &'a str>) -> String {
    let mut out = String::new();
    for value in values {
        out.push_str(" FIELDDATA=\"");
        out.push_str(value);
        out.push('"');
    }
    out
}

fn v3000_sgroup_attach_point_block(sgroup: &SubstanceGroup) -> String {
    let mut out = String::new();
    for attach_point in sgroup.attach_points() {
        out.push_str(&format!(" SAP=(3 {}", attach_point.atom.index() + 1));
        if attach_point.leaving_atom == Some(attach_point.atom) {
            out.push_str(" aidx");
        } else {
            let leaving = attach_point.leaving_atom.map_or(0, |atom| atom.index() + 1);
            out.push_str(&format!(" {leaving}"));
        }
        out.push(' ');
        out.push_str(attach_point.label.as_deref().unwrap_or(""));
        out.push(')');
    }
    out
}

fn v3000_sgroup_bracket_style_value(
    style: Option<&SGroupBracketStyle>,
) -> Result<Option<&str>, MolWriteError> {
    Ok(match style {
        Some(SGroupBracketStyle::Bracket) => Some("BRACKET"),
        Some(SGroupBracketStyle::Parenthesis) => Some("PAREN"),
        Some(SGroupBracketStyle::None) | None => None,
        Some(SGroupBracketStyle::Unknown(value)) => {
            if value.is_empty() {
                None
            } else {
                Some(value.as_str())
            }
        }
    })
}

fn add_v3000_sgroup_block(block: String, current_line: &mut String, out: &mut String) {
    if block.is_empty() {
        return;
    }
    if current_line.is_empty() {
        current_line.push_str("M  V30");
    }
    if current_line.len() + block.len() < 78 {
        current_line.push_str(&block);
        return;
    }
    out.push_str(current_line);
    out.push_str(" -\n");
    let mut start = 0_usize;
    while block.len().saturating_sub(start) >= 73 {
        out.push_str("M  V30");
        if start != 0 {
            out.push(' ');
        }
        let end = safe_char_boundary(&block, start + 72);
        out.push_str(&block[start..end]);
        start = end;
        if start < block.len() {
            out.push_str("-\n");
        }
    }
    if start < block.len() {
        *current_line = "M  V30".to_string();
        if start != 0 {
            current_line.push(' ');
        }
        current_line.push_str(&block[start..]);
    } else {
        current_line.clear();
    }
}

fn safe_char_boundary(value: &str, mut index: usize) -> usize {
    index = index.min(value.len());
    while index > 0 && !value.is_char_boundary(index) {
        index -= 1;
    }
    index
}

fn v3000_double_field(value: f64) -> String {
    v2000_double_field(value).trim().to_string()
}

fn append_v3000_collection_lines(
    out: &mut String,
    molecule: &Molecule,
) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: appendEnhancedStereoGroups
    // RDKit❗✔️: if (!tmol.getStereoGroups().empty()) {
    // RDKit❗✔️:   auto stereo_groups = tmol.getStereoGroups();
    // RDKit❗✔️:   assignStereoGroupIds(stereo_groups);
    // RDKit❗✔️:   res += "M  V30 BEGIN COLLECTION\n";
    // RDKit❗✔️:   for (auto &&group : stereo_groups) {
    // RDKit❗✔️:     tmp += "M  V30 MDLV30/";
    // RDKit❗✔️:     switch (group.getGroupType()) { ... STEABS / STEREL[id] / STERAC[id] ... }
    // RDKit❗✔️:     Atropisomers::getAllAtomIdsForStereoGroup(tmol, group, atomIds, wedgeBonds);
    // RDKit❗✔️:     tmp += std::to_string(atomIds.size()); ... wrap at 78 chars ...
    // RDKit❗✔️:   }
    // RDKit❗✔️:   res += tmp + "M  V30 END COLLECTION\n";
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: appendEnhancedStereoGroups
    let groups = molecule.stereo_groups();
    if groups.is_empty() {
        return Ok(());
    }
    let write_ids = assigned_v3000_stereo_group_ids(groups);
    out.push_str("M  V30 BEGIN COLLECTION\n");
    for (group, write_id) in groups.iter().zip(write_ids) {
        let label = match group.kind() {
            StereoGroupKind::Absolute => "STEABS".to_string(),
            StereoGroupKind::Or => format!("STEREL{}", write_id.unwrap_or(1)),
            StereoGroupKind::And => format!("STERAC{}", write_id.unwrap_or(1)),
        };
        // Collect atoms: if the group has explicit atoms use them;
        // otherwise collect atoms from all bonds (atropisomer/bond-based groups).
        let atom_ids: Vec<crate::AtomId> = if !group.atoms().is_empty() {
            group.atoms().to_vec()
        } else {
            group
                .bonds()
                .iter()
                .flat_map(|bond_id| {
                    let bond = &molecule.bonds()[bond_id.index()];
                    [bond.begin(), bond.end()]
                })
                .collect()
        };
        append_v3000_collection_group_line(out, &label, atom_ids.iter().copied());
    }
    out.push_str("M  V30 END COLLECTION\n");
    Ok(())
}

fn assigned_v3000_stereo_group_ids(groups: &[crate::StereoGroup]) -> Vec<Option<u32>> {
    let mut or_ids = Vec::<u32>::new();
    let mut and_ids = Vec::<u32>::new();
    let mut assigned = groups
        .iter()
        .map(crate::StereoGroup::id)
        .collect::<Vec<_>>();
    for (idx, group) in groups.iter().enumerate() {
        let Some(id) = assigned[idx] else {
            continue;
        };
        let ids = match group.kind() {
            StereoGroupKind::Or => &mut or_ids,
            StereoGroupKind::And => &mut and_ids,
            StereoGroupKind::Absolute => continue,
        };
        if id != 0 && ids.contains(&id) {
            assigned[idx] = None;
        } else if id != 0 {
            ids.push(id);
        }
    }
    let mut next_or = 0_u32;
    let mut next_and = 0_u32;
    for (idx, group) in groups.iter().enumerate() {
        if group.kind() == StereoGroupKind::Absolute || assigned[idx].is_some() {
            continue;
        }
        let (next, ids) = match group.kind() {
            StereoGroupKind::Or => (&mut next_or, &mut or_ids),
            StereoGroupKind::And => (&mut next_and, &mut and_ids),
            StereoGroupKind::Absolute => unreachable!(),
        };
        *next += 1;
        while ids.contains(&*next) {
            *next += 1;
        }
        ids.push(*next);
        assigned[idx] = Some(*next);
    }
    assigned
}

fn append_v3000_collection_group_line(
    out: &mut String,
    label: &str,
    atoms: impl IntoIterator<Item = crate::AtomId>,
) {
    let atoms = atoms.into_iter().collect::<Vec<_>>();
    let mut line = format!("M  V30 MDLV30/{label} ATOMS=({}", atoms.len());
    for atom in atoms {
        let idx = (atom.index() + 1).to_string();
        if line.len() + idx.len() >= 78 {
            out.push_str(&line);
            out.push_str("-\n");
            line = "M  V30 ".to_string();
        } else {
            line.push(' ');
        }
        line.push_str(&idx);
    }
    out.push_str(&line);
    out.push_str(")\n");
}

fn v2000_bond_query_symbol(bond: &Bond) -> Option<u32> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondSymbol / BondGetMolFileSymbol
    // RDKit❗✔️: int res = 8;
    // RDKit❗✔️: if (qry->getDescription() == "BondOrder" || getQueryBondTopology(bond)) { res = 0; }
    // RDKit❗✔️: if (qry->getDescription() == "BondOr" && !qry->getNegation()) { ... SINGLE+DOUBLE => 5; SINGLE+AROMATIC => 6; DOUBLE+AROMATIC => 7; }
    // RDKit❗✔️: else if (qry->getDescription() == "SingleOrAromaticBond" && !qry->getNegation()) { res = 6; }
    // RDKit❗✔️: else if (qry->getDescription() == "SingleOrDoubleBond" && !qry->getNegation()) { res = 5; }
    // RDKit❗✔️: else if (qry->getDescription() == "DoubleOrAromaticBond" && !qry->getNegation()) { res = 7; }
    // RDKit❗✔️: if (bond->hasQuery()) { res = getQueryBondSymbol(bond); }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondSymbol / BondGetMolFileSymbol
    let query = bond.query()?;
    v2000_bond_query_symbol_from_query(query)
}

fn v2000_bond_query_symbol_from_query(query: &QueryNode<BondQueryPredicate>) -> Option<u32> {
    match query {
        QueryNode::Predicate(BondQueryPredicate::Any) => Some(8),
        QueryNode::Predicate(BondQueryPredicate::MolFileQueryCode(code)) => Some(*code),
        QueryNode::Predicate(BondQueryPredicate::OrderIn(orders)) => {
            v2000_bond_order_query_symbol(orders)
        }
        QueryNode::And(children) => {
            let non_topology = children
                .iter()
                .filter(|child| v2000_bond_topology_code_from_query(child).is_none())
                .collect::<Vec<_>>();
            if non_topology.len() == 1 {
                v2000_bond_query_symbol_from_query(non_topology[0])
            } else {
                None
            }
        }
        _ => None,
    }
}

fn v2000_bond_topology_code(bond: &Bond) -> Option<u32> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondTopology / GetMolFileBondLine
    // RDKit❗✔️: if (qry->getDescription() == "BondAnd" && !qry->getNegation() && qry->endChildren() - qry->beginChildren() == 2) { ... qry = BondInRing child; }
    // RDKit❗✔️: if (qry->getDescription() == "BondInRing") { if (qry->getNegation()) { res = 2; } else { res = 1; } }
    // RDKit❗✔️: if (bond->hasQuery()) { int topol = getQueryBondTopology(bond); if (topol) { ss << " " << std::setw(2) << 0 << " " << std::setw(2) << topol; } }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: getQueryBondTopology / GetMolFileBondLine
    v2000_bond_topology_code_from_query(bond.query()?)
}

fn v2000_bond_topology_code_from_query(query: &QueryNode<BondQueryPredicate>) -> Option<u32> {
    match query {
        QueryNode::Predicate(BondQueryPredicate::IsInRing(true)) => Some(1),
        QueryNode::Predicate(BondQueryPredicate::IsInRing(false)) => Some(2),
        QueryNode::And(children) if children.len() == 2 => children
            .iter()
            .find_map(v2000_bond_topology_code_from_query),
        _ => None,
    }
}

fn v2000_bond_order_query_symbol(orders: &[BondOrder]) -> Option<u32> {
    let has_single = orders.contains(&BondOrder::Single);
    let has_double = orders.contains(&BondOrder::Double);
    let has_aromatic = orders.contains(&BondOrder::Aromatic);
    match (orders.len(), has_single, has_double, has_aromatic) {
        (2, true, true, false) => Some(5),
        (2, true, false, true) => Some(6),
        (2, false, true, true) => Some(7),
        _ => None,
    }
}

fn append_v2000_property_lines(out: &mut String, molecule: &Molecule) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileChargeInfo
    // RDKit❗✔️: if (atom->getFormalCharge() != 0) { ++nChgs; chgss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % atom->getFormalCharge(); ... }
    // RDKit❗✔️: unsigned int nRadEs = atom->getNumRadicalElectrons();
    // RDKit❗✔️: if (nRadEs != 0 && atom->getTotalDegree() != 0) { ... if (nRadEs % 2) { nRadEs = 2; } else { nRadEs = 3; } ... }
    // RDKit❗✔️: int isotope = atom->getIsotope();
    // RDKit❗✔️: if (isotope != 0) { ++nMassDiffs; massdiffss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % isotope; ... }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileChargeInfo
    let charges = molecule
        .atoms()
        .iter()
        .filter_map(|atom| {
            (atom.formal_charge() != 0)
                .then_some((atom.id().index() + 1, i32::from(atom.formal_charge())))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "CHG", &charges);

    let valence = molblock_valence_assignment(molecule).map_err(|source| {
        MolWriteError::Value(format!("MolBlock RAD field assignment failed: {source}"))
    })?;
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
        .filter_map(|atom| {
            atom.isotope()
                .map(|isotope| (atom.id().index() + 1, i32::from(isotope)))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "ISO", &isotopes);
    Ok(())
}

fn append_v2000_rgroup_lines(out: &mut String, molecule: &Molecule) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileRGroupInfo
    // RDKit❗✔️: unsigned int nEntries = 0;
    // RDKit❗✔️: if ((*atomIt)->getPropIfPresent(common_properties::_MolFileRLabel, lbl)) {
    // RDKit❗✔️:   ss << " " << std::setw(3) << (*atomIt)->getIdx() + 1 << " "
    // RDKit❗✔️:      << std::setw(3) << lbl;
    // RDKit❗✔️:   ++nEntries;
    // RDKit❗✔️: }
    // RDKit❗✔️: if (nEntries) { ss2 << "M  RGP" << std::setw(3) << nEntries << ss.str() << "\n"; }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileRGroupInfo
    let entries = molecule
        .atoms()
        .iter()
        .filter_map(|atom| atom.prop("_MolFileRLabel").map(|label| (atom, label)))
        .map(|(atom, label)| {
            let label = label.parse::<i32>().map_err(|_| {
                MolWriteError::Value(format!(
                    "invalid _MolFileRLabel value '{}' on atom {}",
                    label,
                    atom.id().index()
                ))
            })?;
            Ok((atom.id().index() + 1, label))
        })
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    if !entries.is_empty() {
        out.push_str(&format!("M  RGP{:>3}", entries.len()));
        for (idx, label) in entries {
            out.push_str(&format!(" {:>3} {:>3}", idx, label));
        }
        out.push('\n');
    }
    Ok(())
}

fn append_v2000_value_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    // RDKit❗✔️: std::string molFileValue;
    // RDKit❗✔️: if (!wrote_query &&
    // RDKit❗✔️:     atom->getPropIfPresent(common_properties::molFileValue, molFileValue)) {
    // RDKit❗✔️:   ss << "V  " << std::setw(3) << atom->getIdx() + 1 << " " << molFileValue
    // RDKit❗✔️:      << "\n";
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    for atom in molecule.atoms() {
        if let Some(value) = atom.prop("molFileValue") {
            out.push_str(&format!("V  {:>3} {value}\n", atom.id().index() + 1));
        }
    }
    append_v2000_atom_list_lines(out, molecule);
}

fn append_v2000_atom_list_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    // RDKit❗✔️: for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:   if (listQs[atom->getIdx()]) {
    // RDKit❗✔️:     INT_VECT vals;
    // RDKit❗✔️:     getAtomListQueryVals(atom->getQuery(), vals);
    // RDKit❗✔️:     ss << "M  ALS " << std::setw(3) << atom->getIdx() + 1 << " ";
    // RDKit❗✔️:     ss << std::setw(2) << vals.size();
    // RDKit❗✔️:     if (atom->getQuery()->getNegation()) { ss << " T "; } else { ss << " F "; }
    // RDKit❗✔️:     for (auto val : vals) { ss << std::setw(4) << std::left << (PeriodicTable::getTable()->getElementSymbol(val)); }
    // RDKit❗✔️:     ss << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileQueryInfo
    for atom in molecule.atoms() {
        if let Some((atomic_numbers, negated)) = v2000_atom_list_query(atom) {
            out.push_str(&format!(
                "M  ALS {:>3} {:>2} {} ",
                atom.id().index() + 1,
                atomic_numbers.len(),
                if negated { "T" } else { "F" }
            ));
            for atomic_number in atomic_numbers {
                let symbol = molfile_atom_symbol(*atomic_number).unwrap_or("*");
                out.push_str(&format!("{symbol:<4}"));
            }
            out.push('\n');
        }
    }
}

fn append_v2000_zbo_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileZBOInfo
    // RDKit❗✔️: if ((*bondIt)->getBondType() == Bond::ZERO) {
    // RDKit❗✔️:   ss << " " << std::setw(3) << (*bondIt)->getIdx() + 1 << " " << std::setw(3) << 0;
    // RDKit❗✔️:   atomsAffected[(*bondIt)->getBeginAtomIdx()] = 1;
    // RDKit❗✔️:   atomsAffected[(*bondIt)->getEndAtomIdx()] = 1;
    // RDKit❗✔️: }
    // RDKit❗✔️: if (nEntries) { res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << "\n"; }
    // RDKit❗✔️: hydss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % atom->getTotalNumHs();
    // RDKit❗✔️: if (atom->getFormalCharge()) { zchss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % atom->getFormalCharge(); }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileZBOInfo
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
        .map(|atom| (atom.id().index() + 1, i32::from(atom.explicit_hydrogens())))
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

fn append_v2000_pxa_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFilePXAInfo
    // RDKit❗✔️: if (atom->hasProp("_MolFile_PXA")) {
    // RDKit❗✔️:   res += boost::str(boost::format("M  PXA % 3d%s\n") % (atom->getIdx() + 1) %
    // RDKit❗✔️:                    atom->getProp<std::string>("_MolFile_PXA"));
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFilePXAInfo
    for atom in molecule.atoms() {
        if let Some(pxa) = atom.prop("_MolFile_PXA") {
            out.push_str(&format!("M  PXA {:>3}{pxa}\n", atom.id().index() + 1));
        }
    }
}

fn append_v2000_sgroup_lines(out: &mut String, molecule: &Molecule) -> Result<(), MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetMolFileSGroupInfo
    // RDKit❗✔️: ret << BuildV2000STYLines(mol);
    // RDKit❗✔️: ret << BuildV2000SLBLines(mol);
    // RDKit❗✔️: ret << BuildV2000StringPropLines(8, mol, "SUBTYPE", "SST", 3);
    // RDKit❗✔️: ret << BuildV2000StringPropLines(8, mol, "CONNECT", "SCN", 3);
    // RDKit❗✔️: ret << BuildV2000SDSLines(mol);
    // RDKit❗✔️: ret << BuildV2000SPLLines(mol);
    // RDKit❗✔️: ret << BuildV2000SNCLines(mol);
    // RDKit❗✔️: ret << BuildV2000SBTLines(mol);
    // RDKit❗✔️: ret << BuildV2000IdxVectorDataLines(15, idx, "SAL", sgroup.getAtoms());
    // RDKit❗✔️: ret << BuildV2000IdxVectorDataLines(15, idx, "SPA", sgroup.getParentAtoms());
    // RDKit❗✔️: ret << BuildV2000IdxVectorDataLines(15, idx, "SBL", sgroup.getBonds());
    // RDKit❗✔️: ret << BuildV2000SDILine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SMTLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SBVLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SDTLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SDDLine(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SCDSEDLines(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SAPLines(idx, sgroup);
    // RDKit❗✔️: ret << BuildV2000SCLLine(idx, sgroup);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolSGroupWriting.cpp :: GetMolFileSGroupInfo
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() {
        return Ok(());
    }

    let sty_entries = sgroups
        .iter()
        .enumerate()
        .map(|(idx, sgroup)| Ok((idx + 1, sgroup_kind_code(sgroup)?)))
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    append_v2000_sgroup_string_entries(out, "STY", 8, 3, &sty_entries);

    let slb_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.external_id().map(|id| (idx + 1, id as i32)))
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "SLB", &slb_entries);

    let sst_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.subtype().map(|value| (idx + 1, value)))
        .collect::<Vec<_>>();
    append_v2000_sgroup_string_entries(out, "SST", 8, 3, &sst_entries);

    let scn_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.connection().map(|value| (idx + 1, value)))
        .map(|(idx, connection)| sgroup_connection_code(connection).map(|code| (idx, code)))
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    append_v2000_sgroup_string_entries(out, "SCN", 8, 3, &scn_entries);

    let sds_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| (sgroup.expansion_state() == Some("E")).then_some(idx + 1))
        .collect::<Vec<_>>();
    append_v2000_sgroup_index_entries(out, "SDS EXP", 15, &sds_entries);

    let spl_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| {
            sgroup
                .parent()
                .map(|parent| (idx + 1, parent.index() as i32 + 1))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "SPL", &spl_entries);

    let snc_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| {
            sgroup
                .component_number()
                .map(|value| (idx + 1, value as i32))
        })
        .collect::<Vec<_>>();
    append_v2000_counted_property(out, "SNC", &snc_entries);

    let sbt_entries = sgroups
        .iter()
        .enumerate()
        .filter_map(|(idx, sgroup)| sgroup.bracket_style().map(|style| (idx + 1, style)))
        .map(|(idx, style)| sgroup_bracket_style_code(style).map(|code| (idx, code)))
        .collect::<Result<Vec<_>, MolWriteError>>()?;
    append_v2000_counted_property(out, "SBT", &sbt_entries);

    for (idx, sgroup) in sgroups.iter().enumerate() {
        let idx = idx + 1;
        append_v2000_sgroup_member_lines(
            out,
            idx,
            "SAL",
            sgroup.atoms().iter().map(|atom| atom.index() + 1),
        );
        append_v2000_sgroup_member_lines(
            out,
            idx,
            "SPA",
            sgroup.parent_atoms().iter().map(|atom| atom.index() + 1),
        );
        append_v2000_sgroup_member_lines(
            out,
            idx,
            "SBL",
            sgroup.bonds().iter().map(|bond| bond.index() + 1),
        );
        append_v2000_sgroup_sdi_lines(out, idx, sgroup);
        append_v2000_sgroup_smt_line(out, idx, sgroup);
        append_v2000_sgroup_sbv_lines(out, idx, sgroup);
        append_v2000_sgroup_sdt_line(out, idx, sgroup);
        append_v2000_sgroup_sdd_line(out, idx, sgroup);
        append_v2000_sgroup_scd_sed_lines(out, idx, sgroup)?;
        append_v2000_sgroup_sap_lines(out, idx, sgroup);
        append_v2000_sgroup_scl_line(out, idx, sgroup);
    }
    Ok(())
}

fn sgroup_kind_code(sgroup: &SubstanceGroup) -> Result<&str, MolWriteError> {
    Ok(match sgroup.kind() {
        SubstanceGroupKind::Data => "DAT",
        SubstanceGroupKind::Superatom => "SUP",
        SubstanceGroupKind::MultipleGroup => "MUL",
        SubstanceGroupKind::StructuralRepeatUnit => "SRU",
        SubstanceGroupKind::Monomer => "MON",
        SubstanceGroupKind::Copolymer => "COP",
        SubstanceGroupKind::Crosslink => "CRO",
        SubstanceGroupKind::Graft => "GRA",
        SubstanceGroupKind::Modification => "MOD",
        SubstanceGroupKind::Mer => "MER",
        SubstanceGroupKind::AnyPolymer => "ANY",
        SubstanceGroupKind::MixtureComponent => "COM",
        SubstanceGroupKind::Mixture => "MIX",
        SubstanceGroupKind::Formulation => "FOR",
        SubstanceGroupKind::Generic(value) => value.as_str(),
    })
}

fn sgroup_connection_code(connection: &SGroupConnection) -> Result<&str, MolWriteError> {
    Ok(match connection {
        SGroupConnection::HeadToHead => "HH",
        SGroupConnection::HeadToTail => "HT",
        SGroupConnection::Either => "EU",
        SGroupConnection::Unknown(value) => value.as_str(),
    })
}

fn sgroup_bracket_style_code(style: &SGroupBracketStyle) -> Result<i32, MolWriteError> {
    match style {
        SGroupBracketStyle::Bracket | SGroupBracketStyle::None => Ok(0),
        SGroupBracketStyle::Parenthesis => Ok(1),
        SGroupBracketStyle::Unknown(value) => Err(MolWriteError::Value(format!(
            "invalid SGroup BRKTYP value '{value}'"
        ))),
    }
}

fn append_v2000_sgroup_string_entries<T: AsRef<str>>(
    out: &mut String,
    code: &str,
    entries_per_line: usize,
    field_width: usize,
    entries: &[(usize, T)],
) {
    for chunk in entries.chunks(entries_per_line) {
        out.push_str(&format!("M  {code}{:>3}", chunk.len()));
        for (idx, value) in chunk {
            out.push_str(&v2000_int_field(*idx));
            out.push_str(&v2000_string_field(value.as_ref(), field_width, true, true));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_index_entries(
    out: &mut String,
    code: &str,
    entries_per_line: usize,
    entries: &[usize],
) {
    for chunk in entries.chunks(entries_per_line) {
        out.push_str(&format!("M  {code}{:>3}", chunk.len()));
        for idx in chunk {
            out.push_str(&v2000_int_field(*idx));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_member_lines<I>(out: &mut String, sgroup_idx: usize, code: &str, values: I)
where
    I: IntoIterator<Item = usize>,
{
    let values = values.into_iter().collect::<Vec<_>>();
    for chunk in values.chunks(15) {
        out.push_str(&format!(
            "M  {code}{}{:>3}",
            v2000_int_field(sgroup_idx),
            chunk.len()
        ));
        for value in chunk {
            out.push_str(&v2000_int_field(*value));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_sdi_lines(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    if let Some(display) = sgroup.display() {
        for bracket in &display.brackets {
            out.push_str(&format!(
                "M  SDI{}{:>3}{}{}{}{}\n",
                v2000_int_field(idx),
                4,
                v2000_double_field(bracket.p1[0]),
                v2000_double_field(bracket.p1[1]),
                v2000_double_field(bracket.p2[0]),
                v2000_double_field(bracket.p2[1])
            ));
        }
    }
}

fn append_v2000_sgroup_smt_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    let value = if matches!(sgroup.kind(), SubstanceGroupKind::MultipleGroup) {
        sgroup.props().get("MULT").map(String::as_str)
    } else {
        sgroup.label()
    };
    if let Some(value) = value {
        out.push_str(&format!(
            "M  SMT{}{}\n",
            v2000_int_field(idx),
            v2000_string_field(value, 69, false, true)
        ));
    }
}

fn append_v2000_sgroup_sbv_lines(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    for cstate in sgroup.cstates() {
        out.push_str(&format!(
            "M  SBV{}{}",
            v2000_int_field(idx),
            v2000_int_field(cstate.bond.index() + 1)
        ));
        if matches!(sgroup.kind(), SubstanceGroupKind::Superatom) {
            out.push_str(&v2000_double_field(cstate.vector[0]));
            out.push_str(&v2000_double_field(cstate.vector[1]));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_sdt_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    let Some(data) = sgroup.data() else {
        return;
    };
    let Some(field_name) = data.field_name.as_deref() else {
        return;
    };
    out.push_str(&format!("M  SDT{}", v2000_int_field(idx)));
    out.push_str(&v2000_string_field(field_name, 30, true, true));
    out.push_str(&v2000_string_field(
        data.field_type.as_deref().unwrap_or("T"),
        2,
        true,
        false,
    ));
    if let Some(value) = data.field_info.as_deref() {
        out.push_str(&v2000_string_field(value, 20, true, false));
    }
    if let Some(value) = data.query_type.as_deref() {
        out.push_str(&v2000_string_field(value, 2, true, false));
    }
    if let Some(value) = data.query_op.as_deref() {
        out.push_str(&v2000_string_field(value, 15, true, false));
    }
    out.push('\n');
}

fn append_v2000_sgroup_sdd_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    if let Some(value) = sgroup.data().and_then(|data| data.field_display.as_deref()) {
        out.push_str(&format!(
            "M  SDD{}{}\n",
            v2000_int_field(idx),
            v2000_string_field(value, 69, false, true)
        ));
    }
}

fn append_v2000_sgroup_scd_sed_lines(
    out: &mut String,
    idx: usize,
    sgroup: &SubstanceGroup,
) -> Result<(), MolWriteError> {
    let values = sgroup
        .data()
        .map(|data| data.values.as_slice())
        .unwrap_or(&[]);
    for value in values.iter().chain(sgroup.data_fields()) {
        if value.len() > 200 {
            return Err(MolWriteError::Value(format!(
                "Data field '{value}' in SGroup {idx} is longer than limit of 200 characters."
            )));
        }
        if value.is_empty() {
            out.push_str(&format!("M  SED{} \n", v2000_int_field(idx)));
            continue;
        }
        let mut chunks = value.as_bytes().chunks(69).peekable();
        while let Some(chunk) = chunks.next() {
            let text = String::from_utf8_lossy(chunk);
            let has_more = chunks.peek().is_some();
            let code = if has_more { "SCD" } else { "SED" };
            out.push_str(&format!(
                "M  {code}{}{}\n",
                v2000_int_field(idx),
                v2000_string_field(&text, 69, has_more, true)
            ));
        }
    }
    Ok(())
}

fn append_v2000_sgroup_sap_lines(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    for chunk in sgroup.attach_points().chunks(6) {
        out.push_str(&format!("M  SAP{}{:>3}", v2000_int_field(idx), chunk.len()));
        for attach_point in chunk {
            out.push_str(&v2000_int_field(attach_point.atom.index() + 1));
            let leaving = attach_point.leaving_atom.map_or(0, |atom| atom.index() + 1);
            out.push_str(&v2000_int_field(leaving));
            out.push_str(&v2000_string_field(
                attach_point.label.as_deref().unwrap_or(""),
                2,
                false,
                true,
            ));
        }
        out.push('\n');
    }
}

fn append_v2000_sgroup_scl_line(out: &mut String, idx: usize, sgroup: &SubstanceGroup) {
    if let Some(value) = sgroup.class() {
        out.push_str(&format!(
            "M  SCL{}{}\n",
            v2000_int_field(idx),
            v2000_string_field(value, 69, false, true)
        ));
    }
}

fn v2000_int_field(value: usize) -> String {
    format!(" {value:>3}")
}

fn v2000_double_field(value: f64) -> String {
    format!("{value:>10.4}")
}

fn v2000_string_field(value: &str, field_size: usize, pad: bool, add_separator: bool) -> String {
    let mut out = String::new();
    if add_separator {
        out.push(' ');
    }
    if value.len() >= field_size {
        out.push_str(&value[..field_size]);
    } else if pad {
        out.push_str(&format!("{value:<field_size$}"));
    } else {
        out.push_str(value);
    }
    out
}

fn append_v2000_alias_lines(out: &mut String, molecule: &Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAliasInfo
    // RDKit❗✔️: if ((*atomIt)->getPropIfPresent(common_properties::molFileAlias, lbl)) {
    // RDKit❗✔️:   if (!lbl.empty()) {
    // RDKit❗✔️:     ss << "A  " << std::setw(3) << (*atomIt)->getIdx() + 1 << "\n"
    // RDKit❗✔️:        << lbl << "\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: GetMolFileAliasInfo
    for atom in molecule.atoms() {
        if let Some(alias) = atom.prop("molFileAlias")
            && !alias.is_empty()
        {
            out.push_str(&format!("A  {:>3}\n{alias}\n", atom.id().index() + 1));
        }
    }
}

fn append_v2000_counted_property(out: &mut String, label: &str, entries: &[(usize, i32)]) {
    for chunk in entries.chunks(8) {
        out.push_str(&format!("M  {label}{:>3}", chunk.len()));
        for (idx, value) in chunk {
            out.push_str(&format!(" {:>3} {:>3}", idx, value));
        }
        out.push('\n');
    }
}

fn atom_degree(molecule: &Molecule, atom_index: usize) -> usize {
    molecule
        .bonds()
        .iter()
        .filter(|bond| bond.begin().index() == atom_index || bond.end().index() == atom_index)
        .count()
}

fn append_sdf_record_fields(mut block: String, molecule: &Molecule) -> String {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/SDWriter.cpp :: _MolToSDStream / _writePropToStream
    // RDKit❗✔️: (*dp_ostream) << MolToMolBlock(mol, true, confId, df_kekulize, df_forceV3000);
    // RDKit❗✔️: (*dp_ostream) << ">  <" << name << ">  ";
    // RDKit❗✔️: (*dp_ostream) << "\n";
    // RDKit❗✔️: (*dp_ostream) << pval << "\n";
    // RDKit❗✔️: (*dp_ostream) << "\n";
    // RDKit❗✔️: (*dp_ostream) << "$$$$\n";
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/SDWriter.cpp :: _MolToSDStream / _writePropToStream
    for (name, value) in molecule.properties().sdf_data_fields() {
        if name.contains('\n') || value.contains("\r\n\r\n") || value.contains("\n\n") {
            continue;
        }
        block.push_str(">  <");
        block.push_str(name);
        block.push_str(">  \n");
        block.push_str(value);
        block.push_str("\n\n");
    }
    block.push_str("$$$$\n");
    block
}

fn v2000_atom_symbol(atom: &Atom, pad_with_spaces: bool) -> Result<String, MolWriteError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: AtomGetMolFileSymbol
    // RDKit❗✔️: if (atom->hasProp(common_properties::_MolFileRLabel)) { res = "R#"; }
    // RDKit❗✔️: else if (atom->getAtomicNum()) { res = atom->getSymbol(); }
    // RDKit❗✔️: else if (atom->hasProp(common_properties::dummyLabel)) { ... if (symb == "*") { res = "R"; } else if (symb == "X") { res = "R"; } else if (symb == "Xa") { res = "R1"; } ... else { res = symb; } }
    // RDKit❗✔️: else { if (hasComplexQuery(atom)) { if (isAtomListQuery(atom)) { res = "L"; } else { res = "*"; } } else { res = "R"; } }
    // RDKit❗✔️: if (padWithSpaces) { while (res.size() < 3) { res += " "; } }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp :: AtomGetMolFileSymbol
    let mut symbol = if atom.prop("_MolFileRLabel").is_some() {
        "R#".to_string()
    } else if atom.atomic_number() != 0 {
        molfile_atom_symbol(atom.atomic_number())?.to_string()
    } else if let Some(dummy_label) = atom.prop("dummyLabel") {
        match dummy_label {
            "*" | "X" => "R".to_string(),
            "Xa" => "R1".to_string(),
            "Xb" => "R2".to_string(),
            "Xc" => "R3".to_string(),
            "Xd" => "R4".to_string(),
            "Xf" => "R5".to_string(),
            "Xg" => "R6".to_string(),
            "Xh" => "R7".to_string(),
            "Xi" => "R8".to_string(),
            "Xj" => "R9".to_string(),
            other => other.to_string(),
        }
    } else if v2000_atom_list_query(atom).is_some() {
        "L".to_string()
    } else {
        molfile_atom_symbol(atom.atomic_number())?.to_string()
    };
    if pad_with_spaces {
        while symbol.len() < 3 {
            symbol.push(' ');
        }
    }
    Ok(symbol)
}

fn v2000_atom_list_query(atom: &Atom) -> Option<(&[u8], bool)> {
    match atom.query()? {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(atomic_numbers)) => {
            Some((atomic_numbers.as_slice(), false))
        }
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberNotIn(atomic_numbers)) => {
            Some((atomic_numbers.as_slice(), true))
        }
        _ => None,
    }
}

fn molfile_atom_symbol(atomic_number: u8) -> Result<&'static str, MolWriteError> {
    match atomic_number {
        0 => Ok("*"),
        1 => Ok("H"),
        2 => Ok("He"),
        3 => Ok("Li"),
        4 => Ok("Be"),
        5 => Ok("B"),
        6 => Ok("C"),
        7 => Ok("N"),
        8 => Ok("O"),
        9 => Ok("F"),
        10 => Ok("Ne"),
        11 => Ok("Na"),
        12 => Ok("Mg"),
        13 => Ok("Al"),
        14 => Ok("Si"),
        15 => Ok("P"),
        16 => Ok("S"),
        17 => Ok("Cl"),
        18 => Ok("Ar"),
        19 => Ok("K"),
        20 => Ok("Ca"),
        21 => Ok("Sc"),
        22 => Ok("Ti"),
        23 => Ok("V"),
        24 => Ok("Cr"),
        25 => Ok("Mn"),
        26 => Ok("Fe"),
        27 => Ok("Co"),
        28 => Ok("Ni"),
        29 => Ok("Cu"),
        30 => Ok("Zn"),
        31 => Ok("Ga"),
        32 => Ok("Ge"),
        33 => Ok("As"),
        34 => Ok("Se"),
        35 => Ok("Br"),
        36 => Ok("Kr"),
        37 => Ok("Rb"),
        38 => Ok("Sr"),
        39 => Ok("Y"),
        40 => Ok("Zr"),
        41 => Ok("Nb"),
        42 => Ok("Mo"),
        43 => Ok("Tc"),
        44 => Ok("Ru"),
        45 => Ok("Rh"),
        46 => Ok("Pd"),
        47 => Ok("Ag"),
        48 => Ok("Cd"),
        49 => Ok("In"),
        50 => Ok("Sn"),
        51 => Ok("Sb"),
        52 => Ok("Te"),
        53 => Ok("I"),
        54 => Ok("Xe"),
        55 => Ok("Cs"),
        56 => Ok("Ba"),
        57 => Ok("La"),
        58 => Ok("Ce"),
        59 => Ok("Pr"),
        60 => Ok("Nd"),
        61 => Ok("Pm"),
        62 => Ok("Sm"),
        63 => Ok("Eu"),
        64 => Ok("Gd"),
        65 => Ok("Tb"),
        66 => Ok("Dy"),
        67 => Ok("Ho"),
        68 => Ok("Er"),
        69 => Ok("Tm"),
        70 => Ok("Yb"),
        71 => Ok("Lu"),
        72 => Ok("Hf"),
        73 => Ok("Ta"),
        74 => Ok("W"),
        75 => Ok("Re"),
        76 => Ok("Os"),
        77 => Ok("Ir"),
        78 => Ok("Pt"),
        79 => Ok("Au"),
        80 => Ok("Hg"),
        81 => Ok("Tl"),
        82 => Ok("Pb"),
        83 => Ok("Bi"),
        84 => Ok("Po"),
        85 => Ok("At"),
        86 => Ok("Rn"),
        87 => Ok("Fr"),
        88 => Ok("Ra"),
        89 => Ok("Ac"),
        90 => Ok("Th"),
        91 => Ok("Pa"),
        92 => Ok("U"),
        93 => Ok("Np"),
        94 => Ok("Pu"),
        95 => Ok("Am"),
        96 => Ok("Cm"),
        97 => Ok("Bk"),
        98 => Ok("Cf"),
        99 => Ok("Es"),
        100 => Ok("Fm"),
        101 => Ok("Md"),
        102 => Ok("No"),
        103 => Ok("Lr"),
        104 => Ok("Rf"),
        105 => Ok("Db"),
        106 => Ok("Sg"),
        107 => Ok("Bh"),
        108 => Ok("Hs"),
        109 => Ok("Mt"),
        110 => Ok("Ds"),
        111 => Ok("Rg"),
        112 => Ok("Cn"),
        113 => Ok("Nh"),
        114 => Ok("Fl"),
        115 => Ok("Mc"),
        116 => Ok("Lv"),
        117 => Ok("Ts"),
        118 => Ok("Og"),
        _ => Err(MolWriteError::UnsupportedSubset(
            "unknown element MolBlock writing is not ported",
        )),
    }
}

// ── Wedge bond helpers (RDKit source-level port) ────────────────────────────

/// ADAPTED FROM draw.rs: neighbours of an atom
fn atom_neighbors(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    let mut ns = Vec::new();
    for b in mol.bonds() {
        if b.begin().index() == atom_idx {
            ns.push(b.end().index());
        } else if b.end().index() == atom_idx {
            ns.push(b.begin().index());
        }
    }
    ns
}

/// Bonds incident to an atom, returns Vec of (bond_idx, other_atom_idx).
fn atom_bonds(mol: &Molecule, atom_idx: usize) -> Vec<(usize, usize)> {
    let mut res = Vec::new();
    for (i, b) in mol.bonds().iter().enumerate() {
        if b.begin().index() == atom_idx {
            res.push((i, b.end().index()));
        } else if b.end().index() == atom_idx {
            res.push((i, b.begin().index()));
        }
    }
    res
}

/// Counts of stereo-spectified double bonds around an atom.
///
/// BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: getDoubleBondPresence
/// RDKit❌❌: std::tuple<unsigned int, unsigned int, unsigned int> getDoubleBondPresence(
/// RDKit❌❌:     const ROMol &mol, const Atom &atom) {
fn get_double_bond_presence(mol: &Molecule, atom_idx: usize) -> (u32, u32, u32) {
    let mut has_double = 0u32;
    let mut has_known_double = 0u32;
    let mut has_any_double = 0u32;
    for b in mol.bonds() {
        if b.begin().index() != atom_idx && b.end().index() != atom_idx {
            continue;
        }
        // RDKit❌❌: if (bond->getBondType() == Bond::BondType::DOUBLE) {
        if b.order() == BondOrder::Double {
            has_double += 1;
            // RDKit❌❌: if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
            if b.stereo() == BondStereo::Any {
                has_any_double += 1;
            // RDKit❌❌: } else if (bond->getStereo() > Bond::BondStereo::STEREOANY) {
            } else if matches!(
                b.stereo(),
                BondStereo::Z | BondStereo::E | BondStereo::Cis | BondStereo::Trans
            ) {
                has_known_double += 1;
            }
        }
    }
    // RDKit❌❌: return std::make_tuple(hasDouble, hasKnownDouble, hasAnyDouble);
    (has_double, has_known_double, has_any_double)
}
/// END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: getDoubleBondPresence

/// Rank chiral atoms by the number of chiral neighbors / Hs they have.
///
/// BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: countChiralNbrs
/// RDKit❌❌: std::pair<bool, INT_VECT> countChiralNbrs(const ROMol &mol, int noNbrs) {
fn count_chiral_nbrs(mol: &Molecule, no_nbrs: i32) -> (bool, Vec<i32>) {
    // RDKit❌❌: INT_VECT nChiralNbrs(mol.getNumAtoms(), noNbrs);
    let mut n_chiral_nbrs = vec![no_nbrs; mol.num_atoms()];

    // RDKit❌❌: // start by looking for bonds that are already wedged
    for bond in mol.bonds() {
        // RDKit❌❌: if (bond->getBondDir() == Bond::BEGINWEDGE ||
        // RDKit❌❌:     bond->getBondDir() == Bond::BEGINDASH ||
        // RDKit❌❌:     bond->getBondDir() == Bond::UNKNOWN) {
        let already_wedged = matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash | BondDirection::Unknown
        );
        if !already_wedged {
            continue;
        }
        // RDKit❌❌: if (bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit❌❌:     bond->getBeginAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
        if mol.atom(bond.begin()).is_some_and(|a| {
            matches!(
                a.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            )
        }) {
            // RDKit❌❌: nChiralNbrs[bond->getBeginAtomIdx()] = noNbrs + 1;
            n_chiral_nbrs[bond.begin().index()] = no_nbrs + 1;
        // RDKit❌❌: } else if (bond->getEndAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit❌❌:                bond->getEndAtom()->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
        } else if mol.atom(bond.end()).is_some_and(|a| {
            matches!(
                a.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            )
        }) {
            // RDKit❌❌: nChiralNbrs[bond->getEndAtomIdx()] = noNbrs + 1;
            n_chiral_nbrs[bond.end().index()] = no_nbrs + 1;
        }
    }

    // RDKit❌❌: // now rank atoms by the number of chiral neighbors or Hs they have:
    // RDKit❌❌: bool chiNbrs = false;
    let mut chi_nbrs = false;
    for (idx, at) in mol.atoms().iter().enumerate() {
        // RDKit❌❌: if (nChiralNbrs[at->getIdx()] > noNbrs) { continue; }
        if n_chiral_nbrs[idx] > no_nbrs {
            continue;
        }
        // RDKit❌❌: auto type = at->getChiralTag();
        // RDKit❌❌: if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
        // RDKit❌❌:   continue;
        // RDKit❌❌: }
        if !matches!(
            at.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            continue;
        }
        // RDKit❌❌: nChiralNbrs[at->getIdx()] = 0;
        n_chiral_nbrs[idx] = 0;
        // RDKit❌❌: chiNbrs = true;
        chi_nbrs = true;
        // RDKit❌❌: for (const auto nat : mol.atomNeighbors(at)) {
        for nbr_idx in atom_neighbors(mol, idx) {
            let Some(nat) = mol.atom(AtomId::new(nbr_idx)) else {
                continue;
            };
            // RDKit❌❌: if (nat->getAtomicNum() == 1) {
            // RDKit❌❌:   nChiralNbrs[at->getIdx()] -= 10;
            // RDKit❌❌:   continue;
            // RDKit❌❌: }
            if nat.atomic_number() == 1 {
                n_chiral_nbrs[idx] -= 10;
                continue;
            }
            // RDKit❌❌: type = nat->getChiralTag();
            // RDKit❌❌: if (type != Atom::CHI_TETRAHEDRAL_CW &&
            // RDKit❌❌:     type != Atom::CHI_TETRAHEDRAL_CCW) {
            // RDKit❌❌:   continue;
            // RDKit❌❌: }
            if !matches!(
                nat.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                continue;
            }
            // RDKit❌❌: nChiralNbrs[at->getIdx()] -= 1;
            n_chiral_nbrs[idx] -= 1;
        }
    }
    // RDKit❌❌: return std::make_pair(chiNbrs, nChiralNbrs);
    (chi_nbrs, n_chiral_nbrs)
}
/// END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: countChiralNbrs

/// Pick a single bond to wedge at a chiral center.
///
/// BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: pickBondToWedge
/// RDKit❌❌: int pickBondToWedge(
/// RDKit❌❌:     const Atom *atom, const ROMol &mol, const INT_VECT &nChiralNbrs,
/// RDKit❌❌:     const std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> &wedgeBonds,
/// RDKit❌❌:     int noNbrs) {
fn pick_bond_to_wedge(
    mol: &Molecule,
    atom_idx: usize,
    n_chiral_nbrs: &[i32],
    wedge_bonds: &BTreeMap<usize, usize>,
    no_nbrs: i32,
    ring_info: &RingInfo,
) -> Option<usize> {
    // RDKit❌❌: std::vector<std::pair<int, int>> nbrScores;
    let mut nbr_scores: Vec<(i32, usize)> = Vec::new();

    for (bid, oaidx) in atom_bonds(mol, atom_idx) {
        let Some(bond) = mol.bonds().get(bid) else {
            continue;
        };
        // RDKit❌❌: if (bond->getBondType() != Bond::SINGLE) { continue; }
        if bond.order() != BondOrder::Single {
            continue;
        }

        // RDKit❌❌: if (wedgeBonds.find(bid) == wedgeBonds.end()) {
        if wedge_bonds.contains_key(&bid) {
            continue;
        }

        // RDKit❌❌: auto *oatom = bond->getOtherAtom(atom);
        let Some(oatom) = mol.atom(AtomId::new(oaidx)) else {
            continue;
        };

        // RDKit❌❌: if (oatom->getAtomicNum() == 1) {
        if oatom.atomic_number() == 1 {
            // RDKit❌❌: nbrScores.emplace_back(-1000000, bid);
            nbr_scores.push((-1_000_000, bid));
            continue;
        }

        // RDKit❌❌: int nbrScore = oatom->getAtomicNum() + 100 * oatom->getDegree() +
        // RDKit❌❌:                1000 * ((oatom->getChiralTag() != Atom::CHI_UNSPECIFIED));
        let degree = atom_degree(mol, oaidx);
        let chiral_penalty = if oatom.chiral_tag() != ChiralTag::Unspecified {
            1
        } else {
            0
        };
        let mut nbr_score = i32::from(oatom.atomic_number())
            + 100 * i32::try_from(degree).unwrap_or(0)
            + 1000 * chiral_penalty;

        // RDKit❌❌: int oIdx = oatom->getIdx();
        // RDKit❌❌: if (nChiralNbrs[oIdx] < noNbrs) {
        if n_chiral_nbrs[oaidx] < no_nbrs {
            // RDKit❌❌: nbrScore -= 100000 * nChiralNbrs[oIdx];
            nbr_score -= 100_000 * n_chiral_nbrs[oaidx];
        }

        // RDKit❌❌: nbrScore += 10000 * mol.getRingInfo()->numAtomRings(oIdx);
        nbr_score +=
            10_000 * i32::try_from(ring_info.num_atom_rings(AtomId::new(oaidx))).unwrap_or(0);
        // RDKit❌❌: nbrScore += 20000 * mol.getRingInfo()->numBondRings(bid);
        nbr_score +=
            20_000 * i32::try_from(ring_info.num_bond_rings(BondId::new(bid))).unwrap_or(0);

        // RDKit❌❌: auto [hasDoubleBond, hasKnownDoubleBond, hasAnyDoubleBond] =
        // RDKit❌❌:     getDoubleBondPresence(mol, *oatom);
        let (has_db, has_known_db, has_any_db) = get_double_bond_presence(mol, oaidx);
        // RDKit❌❌: nbrScore += 11000 * hasDoubleBond;
        nbr_score += 11_000 * i32::try_from(has_db).unwrap_or(0);
        // RDKit❌❌: nbrScore += 12000 * hasKnownDoubleBond;
        nbr_score += 12_000 * i32::try_from(has_known_db).unwrap_or(0);
        // RDKit❌❌: nbrScore += 23000 * hasAnyDoubleBond;
        nbr_score += 23_000 * i32::try_from(has_any_db).unwrap_or(0);

        // RDKit❌❌: if (oatom->hasProp(common_properties::_fromAttachPoint)) {
        if oatom.prop("_fromAttachPoint").is_some() {
            // RDKit❌❌: nbrScore += 500000;
            nbr_score += 500_000;
        }

        // RDKit❌❌: nbrScores.emplace_back(nbrScore, bid);
        nbr_scores.push((nbr_score, bid));
    }

    // RDKit❌❌: if (nbrScores.empty()) { return -1; }
    // RDKit✔️✔️: auto minPr = std::min_element(nbrScores.begin(), nbrScores.end());
    // RDKit✔️✔️: return minPr->second;
    nbr_scores
        .iter()
        .min_by_key(|(score, bid)| (*score, *bid))
        .map(|(_, bid)| *bid)
}
/// END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: pickBondToWedge

/// Orchestrator: return map of bond_idx → chiral_center_atom_idx for wedging.
///
/// BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: pickBondsToWedge
/// RDKit❌❌: std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> pickBondsToWedge(
/// RDKit❌❌:     const ROMol &mol, const BondWedgingParameters *params,
/// RDKit❌❌:     const Conformer *conf) {
fn pick_bonds_to_wedge(mol: &Molecule, ring_info: &RingInfo) -> BTreeMap<usize, usize> {
    // RDKit❌❌: std::vector<unsigned int> indices(mol.getNumAtoms());
    // RDKit❌❌: std::iota(indices.begin(), indices.end(), 0);
    let mut indices: Vec<usize> = (0..mol.num_atoms()).collect();
    // RDKit❌❌: static int noNbrs = 100;
    const NO_NBRS: i32 = 100;

    // RDKit❌❌: auto [chiNbrs, nChiralNbrs] = detail::countChiralNbrs(mol, noNbrs);
    let (chi_nbrs, n_chiral_nbrs) = count_chiral_nbrs(mol, NO_NBRS);

    // RDKit❌❌: if (chiNbrs) {
    // RDKit❌❌:   std::sort(indices.begin(), indices.end(),
    // RDKit❌❌:             [&nChiralNbrs = nChiralNbrs](auto i1, auto i2) {
    // RDKit❌❌:                 return nChiralNbrs[i1] < nChiralNbrs[i2];
    // RDKit❌❌:             });
    // RDKit❌❌: }
    if chi_nbrs {
        rdkit_std_sort_indices_by_chiral_neighbor_count(&mut indices, &n_chiral_nbrs);
    }

    // RDKit❌❌: std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> wedgeInfo;
    let mut wedge_info: BTreeMap<usize, usize> = BTreeMap::new();

    // RDKit❌❌: for (auto idx : indices) {
    for &idx in &indices {
        // RDKit❌❌: if (nChiralNbrs[idx] > noNbrs) { continue; }
        if n_chiral_nbrs[idx] > NO_NBRS {
            continue;
        }

        let Some(atom) = mol.atom(AtomId::new(idx)) else {
            continue;
        };
        // RDKit❌❌: auto type = atom->getChiralTag();
        // RDKit❌❌: if (type != Atom::CHI_TETRAHEDRAL_CW && type != Atom::CHI_TETRAHEDRAL_CCW) {
        // RDKit❌❌:   break;
        // RDKit❌❌: }
        if !matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            break;
        }

        // RDKit❌❌: auto bnd1 = detail::pickBondToWedge(atom, mol, nChiralNbrs, wedgeInfo, noNbrs);
        let bnd1 = pick_bond_to_wedge(mol, idx, &n_chiral_nbrs, &wedge_info, NO_NBRS, ring_info);

        // RDKit❌❌: if (bnd1 >= 0) {
        if let Some(bond_idx) = bnd1 {
            // RDKit❌❌: wedgeInfo[bnd1] = std::move(wi);
            wedge_info.insert(bond_idx, idx);
        }
    }

    // RDKit❌❌: return wedgeInfo;
    wedge_info
}
/// END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/WedgeBonds.cpp :: pickBondsToWedge

fn rdkit_std_sort_indices_by_chiral_neighbor_count(indices: &mut [usize], n_chiral_nbrs: &[i32]) {
    // BEGIN LIBSTDC++ CPP FUNCTIONS /usr/include/c++/13/bits/stl_algo.h :: std::__sort helpers
    // libstdc++✔️✔️: std::__introsort_loop(__first, __last,
    // libstdc++✔️✔️:                       std::__lg(__last - __first) * 2, __comp);
    // libstdc++✔️✔️: std::__final_insertion_sort(__first, __last, __comp);
    // libstdc++✔️✔️: while (__last - __first > int(_S_threshold)) {
    // libstdc++✔️✔️:   --__depth_limit;
    // libstdc++✔️✔️:   _RandomAccessIterator __cut =
    // libstdc++✔️✔️:     std::__unguarded_partition_pivot(__first, __last, __comp);
    // libstdc++✔️✔️:   std::__introsort_loop(__cut, __last, __depth_limit, __comp);
    // libstdc++✔️✔️:   __last = __cut;
    // libstdc++✔️✔️: }
    // libstdc++✔️✔️: std::__move_median_to_first(__first, __first + 1, __mid, __last - 1,
    // libstdc++✔️✔️:                             __comp);
    // libstdc++✔️✔️: return std::__unguarded_partition(__first + 1, __last, __first, __comp);
    // libstdc++✔️✔️: enum { _S_threshold = 16 };
    // END LIBSTDC++ CPP FUNCTIONS /usr/include/c++/13/bits/stl_algo.h :: std::__sort helpers
    if indices.is_empty() {
        return;
    }
    let depth_limit = (usize::BITS - indices.len().leading_zeros() - 1) as usize * 2;
    libstdcxx_introsort_loop(indices, 0, indices.len(), depth_limit, n_chiral_nbrs);
    libstdcxx_final_insertion_sort(indices, 0, indices.len(), n_chiral_nbrs);
}

fn libstdcxx_less(indices: &[usize], left: usize, right: usize, n_chiral_nbrs: &[i32]) -> bool {
    n_chiral_nbrs[indices[left]] < n_chiral_nbrs[indices[right]]
}

fn libstdcxx_value_less(
    value: usize,
    right: usize,
    indices: &[usize],
    n_chiral_nbrs: &[i32],
) -> bool {
    n_chiral_nbrs[value] < n_chiral_nbrs[indices[right]]
}

fn libstdcxx_move_median_to_first(
    indices: &mut [usize],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    n_chiral_nbrs: &[i32],
) {
    if libstdcxx_less(indices, a, b, n_chiral_nbrs) {
        if libstdcxx_less(indices, b, c, n_chiral_nbrs) {
            indices.swap(result, b);
        } else if libstdcxx_less(indices, a, c, n_chiral_nbrs) {
            indices.swap(result, c);
        } else {
            indices.swap(result, a);
        }
    } else if libstdcxx_less(indices, a, c, n_chiral_nbrs) {
        indices.swap(result, a);
    } else if libstdcxx_less(indices, b, c, n_chiral_nbrs) {
        indices.swap(result, c);
    } else {
        indices.swap(result, b);
    }
}

fn libstdcxx_unguarded_partition(
    indices: &mut [usize],
    mut first: usize,
    last: usize,
    pivot: usize,
    n_chiral_nbrs: &[i32],
) -> usize {
    let mut last = last;
    loop {
        while libstdcxx_less(indices, first, pivot, n_chiral_nbrs) {
            first += 1;
        }
        last -= 1;
        while libstdcxx_less(indices, pivot, last, n_chiral_nbrs) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        indices.swap(first, last);
        first += 1;
    }
}

fn libstdcxx_unguarded_partition_pivot(
    indices: &mut [usize],
    first: usize,
    last: usize,
    n_chiral_nbrs: &[i32],
) -> usize {
    let mid = first + (last - first) / 2;
    libstdcxx_move_median_to_first(indices, first, first + 1, mid, last - 1, n_chiral_nbrs);
    libstdcxx_unguarded_partition(indices, first + 1, last, first, n_chiral_nbrs)
}

fn libstdcxx_introsort_loop(
    indices: &mut [usize],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    n_chiral_nbrs: &[i32],
) {
    const S_THRESHOLD: usize = 16;
    while last - first > S_THRESHOLD {
        if depth_limit == 0 {
            indices[first..last].sort_unstable_by_key(|&idx| n_chiral_nbrs[idx]);
            return;
        }
        depth_limit -= 1;
        let cut = libstdcxx_unguarded_partition_pivot(indices, first, last, n_chiral_nbrs);
        libstdcxx_introsort_loop(indices, cut, last, depth_limit, n_chiral_nbrs);
        last = cut;
    }
}

fn libstdcxx_insertion_sort(
    indices: &mut [usize],
    first: usize,
    last: usize,
    n_chiral_nbrs: &[i32],
) {
    if first == last {
        return;
    }
    for i in (first + 1)..last {
        if libstdcxx_less(indices, i, first, n_chiral_nbrs) {
            let value = indices[i];
            for pos in (first..i).rev() {
                indices[pos + 1] = indices[pos];
            }
            indices[first] = value;
        } else {
            libstdcxx_unguarded_linear_insert(indices, i, n_chiral_nbrs);
        }
    }
}

fn libstdcxx_unguarded_linear_insert(
    indices: &mut [usize],
    mut last: usize,
    n_chiral_nbrs: &[i32],
) {
    let value = indices[last];
    let mut next = last - 1;
    while libstdcxx_value_less(value, next, indices, n_chiral_nbrs) {
        indices[last] = indices[next];
        last = next;
        next -= 1;
    }
    indices[last] = value;
}

fn libstdcxx_unguarded_insertion_sort(
    indices: &mut [usize],
    first: usize,
    last: usize,
    n_chiral_nbrs: &[i32],
) {
    for i in first..last {
        libstdcxx_unguarded_linear_insert(indices, i, n_chiral_nbrs);
    }
}

fn libstdcxx_final_insertion_sort(
    indices: &mut [usize],
    first: usize,
    last: usize,
    n_chiral_nbrs: &[i32],
) {
    const S_THRESHOLD: usize = 16;
    if last - first > S_THRESHOLD {
        libstdcxx_insertion_sort(indices, first, first + S_THRESHOLD, n_chiral_nbrs);
        libstdcxx_unguarded_insertion_sort(indices, first + S_THRESHOLD, last, n_chiral_nbrs);
    } else {
        libstdcxx_insertion_sort(indices, first, last, n_chiral_nbrs);
    }
}

#[cfg(test)]
mod tests;
