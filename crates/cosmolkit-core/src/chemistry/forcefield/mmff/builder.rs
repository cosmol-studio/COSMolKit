//! Source-backed RDKit MMFF force-field builder helpers.

use std::collections::VecDeque;

use crate::chemistry::forcefield::core::ForceField;
use crate::chemistry::forcefield::core::ForceFieldVec3;
use crate::chemistry::forcefield::torsion_query::{
    DEFAULT_TORSION_BOND_SMARTS, TorsionBondQueryError, match_torsion_bonds,
};
use crate::chemistry::mol_transforms::{MolTransformError, get_atom_positions};
use crate::notation::fragment::{get_fragment_atom_mapping, get_mol_frags};
use crate::{AdjacencyList, Hybridization, Molecule};

use super::angle_bend::AngleBendContrib;
use super::bond_stretch::BondStretchContrib;
use super::mol_properties::{MmffMolProperties, MmffMolPropertiesError};
use super::nonbonded::NonbondedContrib;
use super::oop_bend::OopBendContrib;
use super::stretch_bend::StretchBendContrib;
use super::torsion_angle::TorsionAngleContrib;

const RELATION_1_2: u8 = 0;
const RELATION_1_3: u8 = 1;
const RELATION_1_4: u8 = 2;
const RELATION_1_X: u8 = 3;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MmffBuilderError {
    #[error(transparent)]
    MolProperties(#[from] MmffMolPropertiesError),
    #[error(transparent)]
    MolTransform(#[from] MolTransformError),
    #[error(transparent)]
    Fragment(#[from] crate::fragment::FragmentError),
    #[error(transparent)]
    Operation(#[from] crate::OperationError),
    #[error("MMFF constructForceField requires a 3D conformer when conf_id={conf_id}")]
    Missing3dConformer { conf_id: isize },
    #[error("MMFF constructForceField conf_id must be >= -1, got {conf_id}")]
    InvalidConformerId { conf_id: isize },
    #[error("RDKit MMFF buildNeighborMatrix is unsupported for empty molecules")]
    EmptyMolecule,
    #[error(transparent)]
    TorsionBondQuery(#[from] TorsionBondQueryError),
}

pub(crate) fn select_mmff_conformer_index(
    mol: &Molecule,
    conf_id: isize,
) -> Result<usize, MmffBuilderError> {
    let conformers = mol.conformers_3d();
    if conformers.is_empty() {
        return Err(MmffBuilderError::Missing3dConformer { conf_id });
    }
    if conf_id == -1 {
        return Ok(0);
    }
    let requested =
        usize::try_from(conf_id).map_err(|_| MmffBuilderError::InvalidConformerId { conf_id })?;
    conformers
        .iter()
        .position(|conformer| conformer.id() == requested)
        .ok_or(MmffBuilderError::Missing3dConformer { conf_id })
}

pub(crate) fn construct_force_field(
    mol: &Molecule,
    non_bonded_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> Result<ForceField, MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::constructForceField(ROMol&, double, int, bool) (Builder.cpp:1093-1103)
    // RDKit✔️❌: ForceFields::ForceField *constructForceField(ROMol &mol, double nonBondedThresh,
    // RDKit✔️❌:                                              int confId,
    // RDKit✔️❌:                                              bool ignoreInterfragInteractions) {
    // RDKit✔️❌:   MMFFMolProperties mmffMolProperties(mol);
    let mmff_mol_properties = MmffMolProperties::new(mol, "MMFF94", 0)?;
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties.isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );
    // RDKit✔️❌:   ForceFields::ForceField *res =
    // RDKit✔️❌:       constructForceField(mol, &mmffMolProperties, nonBondedThresh, confId,
    // RDKit✔️❌:                           ignoreInterfragInteractions);
    let res = construct_force_field_with_props(
        &mmff_mol_properties.molecule,
        &mmff_mol_properties,
        non_bonded_thresh,
        conf_id,
        ignore_interfrag_interactions,
    )?;
    // RDKit✔️❌:
    // RDKit✔️❌:   return res;
    Ok(res)
    // RDKit✔️❌: }
}

pub(crate) fn construct_force_field_with_props(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    non_bonded_thresh: f64,
    conf_id: isize,
    ignore_interfrag_interactions: bool,
) -> Result<ForceField, MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::constructForceField(ROMol&, MMFFMolProperties*, double, int, bool) (Builder.cpp:1111-1142)
    // RDKit✔️❌: ForceFields::ForceField *constructForceField(
    // RDKit✔️❌:     ROMol &mol, MMFFMolProperties *mmffMolProperties, double nonBondedThresh,
    // RDKit✔️❌:     int confId, bool ignoreInterfragInteractions) {
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust reference models RDKit's non-null MMFFMolProperties precondition.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit✔️❌:   std::unique_ptr<ForceFields::ForceField> res(new ForceFields::ForceField());
    let mut res = ForceField::new(3);
    // RDKit✔️❌:   // add the atomic positions:
    // RDKit✔️❌:   Conformer &conf = mol.getConformer(confId);
    let conf_index = select_mmff_conformer_index(mol, conf_id)?;
    let positions = get_atom_positions(mol, conf_index)?;
    // RDKit✔️❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    for coords in positions {
        // RDKit✔️❌:     res->positions().push_back(&(conf.getAtomPos(i)));
        res.positions_mut()
            .push(ForceFieldVec3::new(coords[0], coords[1], coords[2]));
    }

    // RDKit✔️❌:   res->initialize();
    res.initialize();
    // RDKit✔️❌:   if (mmffMolProperties->getMMFFBondTerm()) {
    if mmff_mol_properties.bond_term {
        // RDKit✔️❌:     Tools::addBonds(mol, mmffMolProperties, res.get());
        add_bonds(mol, mmff_mol_properties, &mut res)?;
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   if (mmffMolProperties->getMMFFAngleTerm()) {
    if mmff_mol_properties.angle_term {
        // RDKit✔️❌:     Tools::addAngles(mol, mmffMolProperties, res.get());
        add_angles(mol, mmff_mol_properties, &mut res)?;
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   if (mmffMolProperties->getMMFFStretchBendTerm()) {
    if mmff_mol_properties.stretch_bend_term {
        // RDKit✔️❌:     Tools::addStretchBend(mol, mmffMolProperties, res.get());
        add_stretch_bend(mol, mmff_mol_properties, &mut res)?;
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   if (mmffMolProperties->getMMFFOopTerm()) {
    if mmff_mol_properties.oop_term {
        // RDKit✔️❌:     Tools::addOop(mol, mmffMolProperties, res.get());
        add_oop(mol, mmff_mol_properties, &mut res)?;
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   if (mmffMolProperties->getMMFFTorsionTerm()) {
    if mmff_mol_properties.torsion_term {
        // RDKit✔️❌:     Tools::addTorsions(mol, mmffMolProperties, res.get());
        add_torsions(
            mol,
            mmff_mol_properties,
            &mut res,
            DEFAULT_TORSION_BOND_SMARTS,
        )?;
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   if (mmffMolProperties->getMMFFVdWTerm() ||
    // RDKit✔️❌:       mmffMolProperties->getMMFFEleTerm()) {
    if mmff_mol_properties.vdw_term || mmff_mol_properties.ele_term {
        // RDKit✔️❌:     boost::shared_array<std::uint8_t> neighborMat =
        // RDKit✔️❌:         Tools::buildNeighborMatrix(mol);
        let neighbor_matrix = build_neighbor_matrix(mol);
        // RDKit✔️❌:     Tools::addNonbonded(mol, confId, mmffMolProperties, res.get(), neighborMat,
        // RDKit✔️❌:                         nonBondedThresh, ignoreInterfragInteractions);
        add_nonbonded(
            mol,
            mmff_mol_properties,
            &mut res,
            &neighbor_matrix,
            non_bonded_thresh,
            ignore_interfrag_interactions,
        )?;
        // RDKit✔️❌:   }
    }

    // RDKit✔️❌:   return res.release();
    Ok(res)
    // RDKit✔️❌: }
}

pub(crate) fn add_bonds(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    field: &mut ForceField,
) -> Result<(), MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::addBonds (Builder.cpp:26-97)
    // RDKit✔️✔️: void addBonds(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
    // RDKit✔️✔️:               ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(field, "bad ForceField");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust references reproduce RDKit's non-null force-field and
    // MMFFMolProperties preconditions.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit❌❌:   std::ostream &oStream = mmffMolProperties->getMMFFOStream();
    // RDKit❌❌:   double totalBondStretchEnergy = 0.0;
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << "\n"
    // RDKit❌❌:                  "B O N D   S T R E T C H I N G\n\n"
    // RDKit❌❌:                  "------ATOMS------   ATOM TYPES   FF     BOND     IDEAL       "
    // RDKit❌❌:                  "                 FORCE\n"
    // RDKit❌❌:                  "  I        J          I    J   CLASS   LENGTH   LENGTH    "
    // RDKit❌❌:                  "DIFF.    ENERGY   CONSTANT\n"
    // RDKit❌❌:                  "-------------------------------------------------------------"
    // RDKit❌❌:                  "------------------------"
    // RDKit❌❌:               << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF verbosity stream output.

    // RDKit✔️✔️:   auto contrib = std::make_unique<BondStretchContrib>(field);
    let mut contrib = BondStretchContrib::new(field);
    // RDKit✔️✔️:   bool hasContrib = false;
    let mut has_contrib = false;
    // RDKit✔️✔️:   for (ROMol::ConstBondIterator bi = mol.beginBonds(); bi != mol.endBonds();
    // RDKit✔️✔️:        ++bi) {
    for bond in mol.bonds() {
        // RDKit✔️✔️:     unsigned int idx1 = (*bi)->getBeginAtomIdx();
        let idx1 = bond.begin().index();
        // RDKit✔️✔️:     unsigned int idx2 = (*bi)->getEndAtomIdx();
        let idx2 = bond.end().index();
        // RDKit✔️✔️:     unsigned int bondType;
        // RDKit✔️✔️:     MMFFBond mmffBondParams;
        // RDKit✔️✔️:     if (mmffMolProperties->getMMFFBondStretchParams(mol, idx1, idx2, bondType,
        // RDKit✔️✔️:                                                     mmffBondParams)) {
        if let Some((_bond_type, mmff_bond_params)) =
            mmff_mol_properties.get_mmff_bond_stretch_params(idx1, idx2)?
        {
            // RDKit✔️✔️:       contrib->addTerm(idx1, idx2, &mmffBondParams);
            contrib.add_term(idx1, idx2, &mmff_bond_params);
            // RDKit✔️✔️:       hasContrib = true;
            has_contrib = true;
            // RDKit❌❌:       if (mmffMolProperties->getMMFFVerbosity()) {
            // RDKit❌❌:         unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx1);
            // RDKit❌❌:         unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx2);
            // RDKit❌❌:         const Atom *iAtom = (*bi)->getBeginAtom();
            // RDKit❌❌:         const Atom *jAtom = (*bi)->getEndAtom();
            // RDKit❌❌:         const double dist = field->distance(idx1, idx2);
            // RDKit❌❌:         const double bondStretchEnergy = MMFF::Utils::calcBondStretchEnergy(
            // RDKit❌❌:             mmffBondParams.r0, mmffBondParams.kb, dist);
            // RDKit❌❌:         if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            // RDKit❌❌:           oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
            // RDKit❌❌:                   << std::setw(5) << idx1 + 1 << std::setw(2)
            // RDKit❌❌:                   << jAtom->getSymbol() << " #" << std::setw(5) << idx2 + 1
            // RDKit❌❌:                   << std::right << std::setw(5) << iAtomType << std::setw(5)
            // RDKit❌❌:                   << jAtomType << std::setw(6) << bondType << "  " << std::fixed
            // RDKit❌❌:                   << std::setprecision(3) << std::setw(9) << dist
            // RDKit❌❌:                   << std::setw(9) << mmffBondParams.r0 << std::setw(9)
            // RDKit❌❌:                   << dist - mmffBondParams.r0 << std::setw(10)
            // RDKit❌❌:                   << bondStretchEnergy << std::setw(10) << mmffBondParams.kb
            // RDKit❌❌:                   << std::endl;
            // RDKit❌❌:         }
            // RDKit❌❌:         totalBondStretchEnergy += bondStretchEnergy;
            // RDKit❌❌:       }
            // Verbosity-dependent reporting is not modeled yet.
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     oStream << "TOTAL BOND STRETCH ENERGY      =" << std::right << std::setw(16)
    // RDKit❌❌:             << std::fixed << std::setprecision(4) << totalBondStretchEnergy
    // RDKit❌❌:             << std::endl;
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF summary stream output.
    // RDKit✔️✔️:   if (hasContrib) {
    if has_contrib {
        // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
        field.add_contrib(Box::new(contrib));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    Ok(())
}

pub(crate) fn add_angles(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    field: &mut ForceField,
) -> Result<(), MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::addAngles (Builder.cpp:141-257)
    // RDKit✔️✔️: void addAngles(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
    // RDKit✔️✔️:                ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(field, "bad ForceField");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust references reproduce RDKit's non-null force-field and
    // MMFFMolProperties preconditions.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit❌❌:   std::ostream &oStream = mmffMolProperties->getMMFFOStream();
    // RDKit✔️✔️:   unsigned int idx[3];
    // RDKit❌❌:   const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
    // RDKit❌❌:   ROMol::ADJ_ITER nbr1Idx;
    // RDKit❌❌:   ROMol::ADJ_ITER end1Nbrs;
    // RDKit❌❌:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit❌❌:   ROMol::ADJ_ITER end2Nbrs;
    // COSMolKit does not model RDKit's ostream or iterator types directly.

    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.num_atoms();
    let adjacency = AdjacencyList::from_topology(n_atoms, mol.bonds());
    // RDKit❌❌:   double totalAngleBendEnergy = 0.0;
    // RDKit❌❌:   RDGeom::PointPtrVect points;
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << "\n"
    // RDKit❌❌:                  "A N G L E   B E N D I N G\n\n"
    // RDKit❌❌:                  "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE "
    // RDKit❌❌:                  "   IDEAL                          FORCE\n"
    // RDKit❌❌:                  "  I        J        K          I    J    K   CLASS    ANGLE  "
    // RDKit❌❌:                  "   ANGLE      DIFF.    ENERGY   CONSTANT\n"
    // RDKit❌❌:                  "-------------------------------------------------------------"
    // RDKit❌❌:                  "-----------------------------------------"
    // RDKit❌❌:               << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     points = field->positions();
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF verbosity stream output.
    // RDKit✔️✔️:   auto contrib = std::make_unique<AngleBendContrib>(field);
    let mut contrib = AngleBendContrib::new(field);
    // RDKit✔️✔️:   bool hasContrib = false;
    let mut has_contrib = false;
    // RDKit✔️✔️:   for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
    for idx1 in 0..n_atoms {
        // RDKit✔️✔️:     const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
        let neighbors = adjacency.neighbors_of(idx1);
        // RDKit✔️✔️:     if (jAtom->getDegree() == 1) {
        if neighbors.len() == 1 {
            // RDKit✔️✔️:       continue;
            continue;
        }
        // RDKit✔️✔️:     unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
        let _j_atom_type = mmff_mol_properties.get_mmff_atom_type(idx1)?;
        // RDKit❌❌:     const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
        let mmff_prop_params_central_atom = mmff_mol_properties.get_mmff_central_atom_prop(idx1)?;
        // RDKit✔️✔️:     boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
        // RDKit✔️✔️:     for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
        for nbr1_pos in 0..neighbors.len() {
            // RDKit✔️✔️:       const Atom *iAtom = mol[*nbr1Idx];
            // RDKit✔️✔️:       idx[0] = iAtom->getIdx();
            let idx0 = neighbors[nbr1_pos].atom_index;
            // RDKit✔️✔️:       boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
            // RDKit✔️✔️:       for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
            for nbr2_pos in 0..neighbors.len() {
                // RDKit✔️✔️:         if (nbr2Idx < (nbr1Idx + 1)) {
                if nbr2_pos < (nbr1_pos + 1) {
                    // RDKit✔️✔️:           continue;
                    continue;
                }
                // RDKit✔️✔️:         const Atom *kAtom = mol[*nbr2Idx];
                // RDKit✔️✔️:         idx[2] = kAtom->getIdx();
                let idx2 = neighbors[nbr2_pos].atom_index;
                // RDKit✔️✔️:         unsigned int angleType;
                // RDKit✔️✔️:         MMFFAngle mmffAngleParams;
                // RDKit✔️✔️:         if (mmffMolProperties->getMMFFAngleBendParams(
                // RDKit✔️✔️:                 mol, idx[0], idx[1], idx[2], angleType, mmffAngleParams)) {
                if let Some((_angle_type, mmff_angle_params)) =
                    mmff_mol_properties.get_mmff_angle_bend_params(idx0, idx1, idx2)?
                {
                    // RDKit✔️✔️:           hasContrib = true;
                    has_contrib = true;
                    // RDKit✔️✔️:           contrib->addTerm(idx[0], idx[1], idx[2], &mmffAngleParams,
                    // RDKit✔️✔️:                            mmffPropParamsCentralAtom);
                    contrib.add_term(
                        idx0,
                        idx1,
                        idx2,
                        &mmff_angle_params,
                        &mmff_prop_params_central_atom,
                    );
                    // RDKit❌❌:           if (mmffMolProperties->getMMFFVerbosity()) {
                    // RDKit❌❌:             unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
                    // RDKit❌❌:             unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
                    // RDKit❌❌:             const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                    // RDKit❌❌:                                      (*(points[idx[0]]))[1],
                    // RDKit❌❌:                                      (*(points[idx[0]]))[2]);
                    // RDKit❌❌:             const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                    // RDKit❌❌:                                      (*(points[idx[1]]))[1],
                    // RDKit❌❌:                                      (*(points[idx[1]]))[2]);
                    // RDKit❌❌:             const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                    // RDKit❌❌:                                      (*(points[idx[2]]))[1],
                    // RDKit❌❌:                                      (*(points[idx[2]]))[2]);
                    // RDKit❌❌:             const double cosTheta = MMFF::Utils::calcCosTheta(
                    // RDKit❌❌:                 p1, p2, p3, field->distance(idx[0], idx[1]),
                    // RDKit❌❌:                 field->distance(idx[1], idx[2]));
                    // RDKit❌❌:             const double theta = RAD2DEG * acos(cosTheta);
                    // RDKit❌❌:             const double angleBendEnergy = MMFF::Utils::calcAngleBendEnergy(
                    // RDKit❌❌:                 mmffAngleParams.theta0, mmffAngleParams.ka,
                    // RDKit❌❌:                 mmffPropParamsCentralAtom->linh, cosTheta);
                    // RDKit❌❌:             if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
                    // RDKit❌❌:               oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                    // RDKit❌❌:                       << std::setw(5) << idx[0] + 1 << std::setw(2)
                    // RDKit❌❌:                       << jAtom->getSymbol() << " #" << std::setw(5)
                    // RDKit❌❌:                       << idx[1] + 1 << std::setw(2) << kAtom->getSymbol()
                    // RDKit❌❌:                       << " #" << std::setw(5) << idx[2] + 1 << std::right
                    // RDKit❌❌:                       << std::setw(5) << iAtomType << std::setw(5) << jAtomType
                    // RDKit❌❌:                       << std::setw(5) << kAtomType << std::setw(6) << angleType
                    // RDKit❌❌:                       << "  " << std::fixed << std::setprecision(3)
                    // RDKit❌❌:                       << std::setw(10) << theta << std::setw(10)
                    // RDKit❌❌:                       << mmffAngleParams.theta0 << std::setw(10)
                    // RDKit❌❌:                       << theta - mmffAngleParams.theta0 << std::setw(10)
                    // RDKit❌❌:                       << angleBendEnergy << std::setw(10) << mmffAngleParams.ka
                    // RDKit❌❌:                       << std::endl;
                    // RDKit❌❌:             }
                    // RDKit❌❌:             totalAngleBendEnergy += angleBendEnergy;
                    // RDKit❌❌:           }
                    // Verbosity-dependent reporting is not modeled yet.
                    // RDKit✔️✔️:         }
                }
            }
        }
    }
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     oStream << "TOTAL ANGLE BEND ENERGY        =" << std::right << std::setw(16)
    // RDKit❌❌:             << std::fixed << std::setprecision(4) << totalAngleBendEnergy
    // RDKit❌❌:             << std::endl;
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF summary stream output.
    // RDKit✔️✔️:   if (hasContrib) {
    if has_contrib {
        // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
        field.add_contrib(Box::new(contrib));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    Ok(())
}

pub(crate) fn add_stretch_bend(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    field: &mut ForceField,
) -> Result<(), MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::addStretchBend (Builder.cpp:259-407)
    // RDKit✔️✔️: void addStretchBend(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
    // RDKit✔️✔️:                     ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(field, "bad ForceField");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust references reproduce RDKit's non-null force-field and
    // MMFFMolProperties preconditions.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit❌❌:   std::ostream &oStream = mmffMolProperties->getMMFFOStream();
    // RDKit✔️✔️:   unsigned int idx[3];
    // RDKit❌❌:   const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
    // RDKit❌❌:   ROMol::ADJ_ITER nbr1Idx;
    // RDKit❌❌:   ROMol::ADJ_ITER end1Nbrs;
    // RDKit❌❌:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit❌❌:   ROMol::ADJ_ITER end2Nbrs;
    // COSMolKit does not model RDKit's ostream or iterator types directly.

    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.num_atoms();
    let adjacency = AdjacencyList::from_topology(n_atoms, mol.bonds());
    // RDKit❌❌:   double totalStretchBendEnergy = 0.0;
    // RDKit❌❌:   RDGeom::PointPtrVect points;
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << "\n"
    // RDKit❌❌:                  "S T R E T C H   B E N D I N G\n\n"
    // RDKit❌❌:                  "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE "
    // RDKit❌❌:                  "    DELTA     DELTA     DELTA               F CON\n"
    // RDKit❌❌:                  "  I        J        K          I    J    K   CLASS    ANGLE  "
    // RDKit❌❌:                  "    ANGLE     R(I,J)    R(J,K)   ENERGY    I-J (J-K)\n"
    // RDKit❌❌:                  "-------------------------------------------------------------"
    // RDKit❌❌:                  "-----------------------------------------------------"
    // RDKit❌❌:               << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     points = field->positions();
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF verbosity stream output.

    // RDKit✔️✔️:   auto contrib = std::make_unique<StretchBendContrib>(field);
    let mut contrib = StretchBendContrib::new(field);
    // RDKit✔️✔️:   bool contribAdded = false;
    let mut contrib_added = false;

    // RDKit✔️✔️:   for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
    for idx1 in 0..n_atoms {
        let neighbors = adjacency.neighbors_of(idx1);
        // RDKit✔️✔️:     const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
        // RDKit✔️✔️:     if (jAtom->getDegree() == 1) {
        if neighbors.len() == 1 {
            // RDKit✔️✔️:       continue;
            continue;
        }
        // RDKit✔️✔️:     unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
        let _j_atom_type = mmff_mol_properties.get_mmff_atom_type(idx1)?;
        // RDKit❌❌:     const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
        let mmff_prop_params_central_atom = mmff_mol_properties.get_mmff_central_atom_prop(idx1)?;
        // RDKit✔️✔️:     if (mmffPropParamsCentralAtom->linh) {
        if mmff_prop_params_central_atom.linh != 0 {
            // RDKit✔️✔️:       continue;
            continue;
        }
        // RDKit✔️✔️:     boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
        // RDKit✔️✔️:     unsigned int i = 0;
        for nbr1_pos in 0..neighbors.len() {
            let idx0 = neighbors[nbr1_pos].atom_index;
            // RDKit✔️✔️:       const Atom *iAtom = mol[*nbr1Idx];
            // RDKit✔️✔️:       boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
            // RDKit✔️✔️:       unsigned int j = 0;
            for nbr2_pos in 0..neighbors.len() {
                let idx2 = neighbors[nbr2_pos].atom_index;
                // RDKit✔️✔️:         const Atom *kAtom = mol[*nbr2Idx];
                // RDKit✔️✔️:         if (j < (i + 1)) {
                if nbr2_pos < (nbr1_pos + 1) {
                    // RDKit✔️✔️:           ++j;
                    // RDKit✔️✔️:           continue;
                    continue;
                }
                // RDKit✔️✔️:         idx[0] = iAtom->getIdx();
                // RDKit✔️✔️:         idx[2] = kAtom->getIdx();
                // RDKit✔️✔️:         unsigned int stretchBendType;
                // RDKit✔️✔️:         MMFFStbn mmffStbnParams;
                // RDKit✔️✔️:         MMFFBond mmffBondParams[2];
                // RDKit✔️✔️:         MMFFAngle mmffAngleParams;
                // RDKit✔️✔️:         if (mmffMolProperties->getMMFFStretchBendParams(
                // RDKit✔️✔️:                 mol, idx[0], idx[1], idx[2], stretchBendType, mmffStbnParams,
                // RDKit✔️✔️:                 mmffBondParams, mmffAngleParams)) {
                if let Some((
                    _stretch_bend_type,
                    mmff_stbn_params,
                    mmff_bond_params,
                    mmff_angle_params,
                )) = mmff_mol_properties.get_mmff_stretch_bend_params(idx0, idx1, idx2)?
                {
                    // RDKit✔️✔️:           contribAdded = true;
                    contrib_added = true;
                    // RDKit✔️✔️:           contrib->addTerm(idx[0], idx[1], idx[2], &mmffStbnParams,
                    // RDKit✔️✔️:                            &mmffAngleParams, &mmffBondParams[0],
                    // RDKit✔️✔️:                            &mmffBondParams[1]);
                    contrib.add_term(
                        idx0,
                        idx1,
                        idx2,
                        &mmff_stbn_params,
                        &mmff_angle_params,
                        &mmff_bond_params[0],
                        &mmff_bond_params[1],
                    );
                    // RDKit❌❌:           if (mmffMolProperties->getMMFFVerbosity()) {
                    // RDKit❌❌:             ...
                    // RDKit❌❌:             totalStretchBendEnergy +=
                    // RDKit❌❌:                 (stretchBendEnergies.first + stretchBendEnergies.second);
                    // RDKit❌❌:           }
                    // Verbosity-dependent reporting is not modeled yet.
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         ++j;
            }
            // RDKit✔️✔️:       ++i;
        }
    }
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     oStream << "TOTAL STRETCH-BEND ENERGY      =" << std::right << std::setw(16)
    // RDKit❌❌:             << std::fixed << std::setprecision(4) << totalStretchBendEnergy
    // RDKit❌❌:             << std::endl;
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF summary stream output.
    // RDKit✔️✔️:   if (contribAdded) {
    if contrib_added {
        // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
        field.add_contrib(Box::new(contrib));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    Ok(())
}

pub(crate) fn add_oop(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    field: &mut ForceField,
) -> Result<(), MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::addOop (Builder.cpp:433-536)
    // RDKit✔️✔️: void addOop(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
    // RDKit✔️✔️:             ForceFields::ForceField *field) {
    // RDKit✔️✔️:   PRECONDITION(field, "bad ForceField");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust references reproduce RDKit's non-null force-field and
    // MMFFMolProperties preconditions.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit❌❌:   std::ostream &oStream = mmffMolProperties->getMMFFOStream();
    // RDKit✔️✔️:   unsigned int idx[4];
    // RDKit✔️✔️:   unsigned int atomType[4];
    // RDKit✔️✔️:   unsigned int n[4];
    // RDKit❌❌:   const Atom *atom[4];
    // RDKit❌❌:   ROMol::ADJ_ITER nbrIdx;
    // RDKit❌❌:   ROMol::ADJ_ITER endNbrs;
    // COSMolKit does not model RDKit's ostream, Atom pointer, or iterator types directly.
    let adjacency = AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());

    // RDKit❌❌:   double totalOopBendEnergy = 0.0;
    // RDKit❌❌:   RDGeom::PointPtrVect points;
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << "\n"
    // RDKit❌❌:                  "O U T - O F - P L A N E   B E N D I N G\n\n"
    // RDKit❌❌:                  "--------------ATOMS---------------         ATOM TYPES        "
    // RDKit❌❌:                  " OOP                FORCE\n"
    // RDKit❌❌:                  "  I        J        K        L          I    J    K    L     "
    // RDKit❌❌:                  "ANGLE    ENERGY   CONSTANT\n"
    // RDKit❌❌:                  "-------------------------------------------------------------"
    // RDKit❌❌:                  "-----------------------------"
    // RDKit❌❌:               << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     points = field->positions();
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF verbosity stream output.

    // RDKit✔️✔️:   bool hasContrib = false;
    let mut has_contrib = false;
    // RDKit✔️✔️:   auto contrib = std::make_unique<OopBendContrib>(field);
    let mut contrib = OopBendContrib::new(field);
    // RDKit✔️✔️:   for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
    for idx1 in 0..mol.num_atoms() {
        let neighbors = adjacency.neighbors_of(idx1);
        // RDKit✔️✔️:     atom[1] = mol.getAtomWithIdx(idx[1]);
        // RDKit✔️✔️:     if (atom[1]->getDegree() != 3) {
        if neighbors.len() != 3 {
            // RDKit✔️✔️:       continue;
            continue;
        }

        let mut idx = [0usize; 4];
        idx[1] = idx1;
        let mut atom_type = [0u8; 4];
        // RDKit✔️✔️:     atomType[1] = mmffMolProperties->getMMFFAtomType(idx[1]);
        atom_type[1] = mmff_mol_properties.get_mmff_atom_type(idx1)?;
        // RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
        // RDKit✔️✔️:     unsigned int i = 0;
        let mut i = 0usize;
        // RDKit✔️✔️:     for (; nbrIdx != endNbrs; ++nbrIdx) {
        for neighbor in neighbors {
            // RDKit✔️✔️:       atom[i] = mol[*nbrIdx];
            // RDKit✔️✔️:       idx[i] = atom[i]->getIdx();
            idx[i] = neighbor.atom_index;
            // RDKit✔️✔️:       atomType[i] = mmffMolProperties->getMMFFAtomType(idx[i]);
            atom_type[i] = mmff_mol_properties.get_mmff_atom_type(idx[i])?;
            // RDKit✔️✔️:       if (!i) {
            if i == 0 {
                // RDKit✔️✔️:         ++i;
                i += 1;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       ++i;
            i += 1;
            // RDKit✔️✔️:     }
        }

        // RDKit✔️✔️:     MMFFOop mmffOopParams;
        // RDKit✔️✔️:     if (!(mmffMolProperties->getMMFFOopBendParams(mol, idx[0], idx[1], idx[2],
        // RDKit✔️✔️:                                                   idx[3], mmffOopParams))) {
        let Some(mmff_oop_params) =
            mmff_mol_properties.get_mmff_oop_bend_params(idx[0], idx[1], idx[2], idx[3])?
        else {
            // RDKit✔️✔️:       continue;
            continue;
        };

        // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
        for i in 0..3 {
            // RDKit✔️✔️:       n[1] = 1;
            let mut n = [0usize; 4];
            n[1] = 1;
            // RDKit✔️✔️:       switch (i) {
            match i {
                // RDKit✔️✔️:         case 0:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 2;
                // RDKit✔️✔️:           n[3] = 3;
                // RDKit✔️✔️:           break;
                0 => {
                    n[0] = 0;
                    n[2] = 2;
                    n[3] = 3;
                }
                // RDKit✔️✔️:         case 1:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 2;
                // RDKit✔️✔️:           break;
                1 => {
                    n[0] = 0;
                    n[2] = 3;
                    n[3] = 2;
                }
                // RDKit✔️✔️:         case 2:
                // RDKit✔️✔️:           n[0] = 2;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 0;
                // RDKit✔️✔️:           break;
                2 => {
                    n[0] = 2;
                    n[2] = 3;
                    n[3] = 0;
                }
                _ => unreachable!(),
            }
            // RDKit✔️✔️:       contrib->addTerm(idx[n[0]], idx[n[1]], idx[n[2]], idx[n[3]],
            // RDKit✔️✔️:                        &mmffOopParams);
            contrib.add_term(idx[n[0]], idx[n[1]], idx[n[2]], idx[n[3]], &mmff_oop_params);
            // RDKit✔️✔️:       hasContrib = true;
            has_contrib = true;
            // RDKit❌❌:       if (mmffMolProperties->getMMFFVerbosity()) {
            // RDKit❌❌:         const RDGeom::Point3D p1((*(points[idx[n[0]]]))[0],
            // RDKit❌❌:                                  (*(points[idx[n[0]]]))[1],
            // RDKit❌❌:                                  (*(points[idx[n[0]]]))[2]);
            // RDKit❌❌:         const RDGeom::Point3D p2((*(points[idx[n[1]]]))[0],
            // RDKit❌❌:                                  (*(points[idx[n[1]]]))[1],
            // RDKit❌❌:                                  (*(points[idx[n[1]]]))[2]);
            // RDKit❌❌:         const RDGeom::Point3D p3((*(points[idx[n[2]]]))[0],
            // RDKit❌❌:                                  (*(points[idx[n[2]]]))[1],
            // RDKit❌❌:                                  (*(points[idx[n[2]]]))[2]);
            // RDKit❌❌:         const RDGeom::Point3D p4((*(points[idx[n[3]]]))[0],
            // RDKit❌❌:                                  (*(points[idx[n[3]]]))[1],
            // RDKit❌❌:                                  (*(points[idx[n[3]]]))[2]);
            // RDKit❌❌:         const double chi = MMFF::Utils::calcOopChi(p1, p2, p3, p4);
            // RDKit❌❌:         const double oopBendEnergy =
            // RDKit❌❌:             MMFF::Utils::calcOopBendEnergy(chi, mmffOopParams.koop);
            // RDKit❌❌:         if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            // RDKit❌❌:           oStream << std::left << std::setw(2) << atom[0]->getSymbol() << " #"
            // RDKit❌❌:                   << std::setw(5) << idx[n[0]] + 1 << std::setw(2)
            // RDKit❌❌:                   << atom[1]->getSymbol() << " #" << std::setw(5)
            // RDKit❌❌:                   << idx[n[1]] + 1 << std::setw(2) << atom[2]->getSymbol()
            // RDKit❌❌:                   << " #" << std::setw(5) << idx[n[2]] + 1 << std::setw(2)
            // RDKit❌❌:                   << atom[3]->getSymbol() << " #" << std::setw(5)
            // RDKit❌❌:                   << idx[n[3]] + 1 << std::right << std::setw(5)
            // RDKit❌❌:                   << atomType[n[0]] << std::setw(5) << atomType[n[1]]
            // RDKit❌❌:                   << std::setw(5) << atomType[n[2]] << std::setw(5)
            // RDKit❌❌:                   << atomType[n[3]] << std::fixed << std::setprecision(3)
            // RDKit❌❌:                   << std::setw(10) << chi << std::setw(10) << oopBendEnergy
            // RDKit❌❌:                   << std::setw(10) << mmffOopParams.koop << std::endl;
            // RDKit❌❌:         }
            // RDKit❌❌:         totalOopBendEnergy += oopBendEnergy;
            // RDKit❌❌:       }
            // Verbosity-dependent reporting is not modeled yet.
        }
        let _ = atom_type;
        // RDKit✔️✔️:   }
    }
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❌❌:       oStream << std::endl;
    // RDKit❌❌:     }
    // RDKit❌❌:     oStream << "TOTAL OUT-OF-PLANE BEND ENERGY =" << std::right << std::setw(16)
    // RDKit❌❌:             << std::fixed << std::setprecision(4) << totalOopBendEnergy
    // RDKit❌❌:             << std::endl;
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF summary stream output.
    // RDKit✔️✔️:   if (hasContrib) {
    if has_contrib {
        // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
        field.add_contrib(Box::new(contrib));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    Ok(())
}

fn is_mmff_double_zero(value: f64) -> bool {
    value < 1.0e-10 && value > -1.0e-10
}

fn two_bit_cell_pos(n_atoms: usize, mut i: usize, mut j: usize) -> usize {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::twoBitCellPos (Builder.cpp:102-108)
    // RDKit✔️✔️: unsigned int twoBitCellPos(unsigned int N, int i, int j) {
    // RDKit✔️✔️:   if (j < i) {
    if j < i {
        // RDKit✔️✔️:     std::swap(i, j);
        std::mem::swap(&mut i, &mut j);
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return i * (N - 1) - (i - 1) * i / 2 + j;
    // The rearranged triangular-number form avoids unsigned `i - 1` at i=0.
    i * n_atoms - i * i.saturating_sub(1) / 2 + (j - i)
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MMFF::Tools::twoBitCellPos
}

fn set_two_bit_cell(res: &mut [u8], pos: usize, value: u8) {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::setTwoBitCell (Builder.cpp:110-113)
    // RDKit✔️✔️: void setTwoBitCell(boost::shared_array<std::uint8_t> &res, unsigned int pos,
    // RDKit✔️✔️:                    std::uint8_t value) {
    // RDKit✔️✔️:   res[pos] = value;
    res[pos] = value;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MMFF::Tools::setTwoBitCell
}

fn get_two_bit_cell(res: &[u8], pos: usize) -> u8 {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::getTwoBitCell (Builder.cpp:115-119)
    // RDKit✔️✔️: std::uint8_t getTwoBitCell(boost::shared_array<std::uint8_t> &res,
    // RDKit✔️✔️:                            unsigned int pos) {
    // RDKit✔️✔️:   return res[pos];
    res[pos]
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MMFF::Tools::getTwoBitCell
}

fn build_neighbor_matrix(mol: &Molecule) -> Vec<u8> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::buildNeighborMatrix (Builder.cpp:129-158)
    // RDKit✔️✔️: // ------------------------------------------------------------------------
    // RDKit✔️✔️: //
    // RDKit✔️✔️: // the two-bit matrix returned by this contains:
    // RDKit✔️✔️: //   0: if atoms i and j are directly connected
    // RDKit✔️✔️: //   1: if atoms i and j are connected via an atom
    // RDKit✔️✔️: //   2: if atoms i and j are in a 1,4 relationship
    // RDKit✔️✔️: //   3: otherwise
    // RDKit✔️✔️: //
    // RDKit✔️✔️: // ------------------------------------------------------------------------
    // RDKit✔️✔️: boost::shared_array<std::uint8_t> buildNeighborMatrix(const ROMol &mol) {
    // RDKit✔️✔️:   const std::uint8_t RELATION_1_X_INIT = RELATION_1_X | (RELATION_1_X << 2) |
    // RDKit✔️✔️:                                          (RELATION_1_X << 4) |
    // RDKit✔️✔️:                                          (RELATION_1_X << 6);
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.num_atoms();
    // RDKit✔️✔️:   unsigned nTwoBitCells = nAtoms * (nAtoms + 1) / 2;
    let n_two_bit_cells = n_atoms * (n_atoms + 1) / 2;
    // RDKit✔️✔️:   boost::shared_array<std::uint8_t> res(new std::uint8_t[nTwoBitCells]);
    // RDKit✔️✔️:   std::memset(res.get(), RELATION_1_X_INIT, nTwoBitCells);
    let mut res = vec![RELATION_1_X; n_two_bit_cells];

    // RDKit✔️✔️:
    // RDKit✔️✔️:   constexpr bool useBO = false;
    // RDKit✔️✔️:   constexpr bool useAtomWts = false;
    // RDKit✔️✔️:   auto dmat = MolOps::getDistanceMat(mol, useBO, useAtomWts);
    let mut dmat = vec![usize::MAX; n_atoms * n_atoms];
    let mut queue = VecDeque::with_capacity(n_atoms);
    for source in 0..n_atoms {
        dmat[source * n_atoms + source] = 0;
        queue.push_back(source);
        while let Some(current) = queue.pop_front() {
            let distance = dmat[source * n_atoms + current];
            for neighbor in mol.topology_block().adjacency.neighbors_of(current) {
                let entry = &mut dmat[source * n_atoms + neighbor.atom_index];
                if *entry == usize::MAX {
                    *entry = distance + 1;
                    queue.push_back(neighbor.atom_index);
                }
            }
        }
    }
    // RDKit✔️✔️:   for (unsigned i = 0; i < nAtoms; ++i) {
    for i in 0..n_atoms {
        // RDKit✔️✔️:     setTwoBitCell(res, twoBitCellPos(nAtoms, i, i), RELATION_1_X);
        set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, i, i), RELATION_1_X);
        // RDKit✔️✔️:     for (unsigned j = i + 1; j < nAtoms; ++j) {
        for j in (i + 1)..n_atoms {
            let distance = dmat[i * n_atoms + j];
            // RDKit✔️✔️:       if (dmat[i * nAtoms + j] == 1.0) {
            if distance == 1 {
                // RDKit✔️✔️:         setTwoBitCell(res, twoBitCellPos(nAtoms, i, j), RELATION_1_2);
                set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, i, j), RELATION_1_2);
                // RDKit✔️✔️:       } else if (dmat[i * nAtoms + j] == 2.0) {
            } else if distance == 2 {
                // RDKit✔️✔️:         setTwoBitCell(res, twoBitCellPos(nAtoms, i, j), RELATION_1_3);
                set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, i, j), RELATION_1_3);
                // RDKit✔️✔️:       } else if (dmat[i * nAtoms + j] == 3.0) {
            } else if distance == 3 {
                // RDKit✔️✔️:         setTwoBitCell(res, twoBitCellPos(nAtoms, i, j), RELATION_1_4);
                set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, i, j), RELATION_1_4);
                // RDKit✔️✔️:       } else if (dmat[i * nAtoms + j] < 1e7) {
            } else if distance != usize::MAX {
                // RDKit✔️✔️:         // the distance matrix sets the distance to 1e8 for atoms that have no
                // RDKit✔️✔️:         // connecting path, so here we know that there's at least a connection
                // RDKit✔️✔️:         setTwoBitCell(res, twoBitCellPos(nAtoms, i, j), RELATION_1_X);
                set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, i, j), RELATION_1_X);
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    res
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MMFF::Tools::buildNeighborMatrix
}

pub(crate) fn add_nonbonded(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    field: &mut ForceField,
    neighbor_matrix: &[u8],
    non_bonded_thresh: f64,
    ignore_interfrag_interactions: bool,
) -> Result<(), MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::addNonbonded (Builder.cpp:923-1087)
    // RDKit✔️✔️: void addNonbonded(const ROMol &mol, int confId,
    // RDKit✔️✔️:                   MMFFMolProperties *mmffMolProperties,
    // RDKit✔️✔️:                   ForceFields::ForceField *field,
    // RDKit✔️✔️:                   boost::shared_array<std::uint8_t> neighborMatrix,
    // RDKit✔️✔️:                   double nonBondedThresh, bool ignoreInterfragInteractions) {
    // RDKit✔️✔️:   PRECONDITION(field, "bad ForceField");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust references reproduce RDKit's non-null force-field and
    // MMFFMolProperties preconditions.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit✔️✔️:   INT_VECT fragMapping;
    // RDKit✔️✔️:   if (ignoreInterfragInteractions) {
    // RDKit✔️✔️:     std::vector<ROMOL_SPTR> molFrags =
    // RDKit✔️✔️:         MolOps::getMolFrags(mol, true, &fragMapping);
    let frag_mapping = if ignore_interfrag_interactions {
        let mapping = get_fragment_atom_mapping(mol);
        for fragment in get_mol_frags(mol)? {
            let _sanitized_fragment = fragment.sanitize()?;
        }
        Some(mapping)
    } else {
        None
    };
    // RDKit✔️✔️:   }

    // RDKit❌❌:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = mol.num_atoms();
    // RDKit❌❌:   // FIX: need a solution for the verbosity here
    // RDKit❌❌:   // std::ostream &oStream = mmffMolProperties->getMMFFOStream();
    // RDKit❌❌:   std::stringstream vdwStream;
    // RDKit❌❌:   std::stringstream eleStream;
    // RDKit❌❌:   double totalVdWEnergy = 0.0;
    // RDKit❌❌:   double totalEleEnergy = 0.0;
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     ...
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF verbosity stream output.

    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // COSMolKit's builder path reads coordinates from the target force field
    // positions rather than from a Molecule conformer handle.
    // RDKit✔️✔️:   auto contrib = std::make_unique<NonbondedContrib>(field);
    let mut contrib = NonbondedContrib::new(field);
    // RDKit✔️✔️:   bool hasContrib = false;
    let mut has_contrib = false;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    for i in 0..n_atoms {
        // RDKit✔️✔️:     for (unsigned int j = i + 1; j < nAtoms; ++j) {
        for j in (i + 1)..n_atoms {
            // RDKit✔️✔️:       std::uint8_t cell =
            // RDKit✔️✔️:           getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j));
            let cell = get_two_bit_cell(neighbor_matrix, two_bit_cell_pos(n_atoms, i, j));
            // RDKit✔️✔️:       if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
            if ignore_interfrag_interactions
                && frag_mapping
                    .as_ref()
                    .is_some_and(|mapping| mapping[i] != mapping[j])
            {
                // RDKit✔️✔️:         continue;
                continue;
            }
            // RDKit✔️✔️:       if (cell >= RELATION_1_4) {
            if cell >= RELATION_1_4 {
                // RDKit✔️✔️:         bool is1_4 = (cell == RELATION_1_4);
                let is_1_4 = cell == RELATION_1_4;
                let pos_i = field.positions()[i];
                let pos_j = field.positions()[j];
                // RDKit✔️✔️:         double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
                let dist = (pos_i - pos_j).length();
                // RDKit✔️✔️:         if (dist > nonBondedThresh) {
                if dist > non_bonded_thresh {
                    // RDKit✔️✔️:           continue;
                    continue;
                }
                // RDKit✔️✔️:         MMFFVdWRijstarEps *vdwConstants = nullptr;
                let mut vdw_constants = None;
                // RDKit✔️✔️:         MMFFVdWRijstarEps mmffVdWConstants;
                // RDKit✔️✔️:         if (mmffMolProperties->getMMFFVdWTerm() &&
                // RDKit✔️✔️:             mmffMolProperties->getMMFFVdWParams(i, j, mmffVdWConstants)) {
                if mmff_mol_properties.vdw_term {
                    vdw_constants = mmff_mol_properties.get_mmff_vdw_params(i, j)?;
                    if vdw_constants.is_some() {
                        // RDKit✔️✔️:           vdwConstants = &mmffVdWConstants;
                        // RDKit✔️✔️:           hasContrib = true;
                        has_contrib = true;
                        // RDKit❌❌:           if (mmffMolProperties->getMMFFVerbosity()) {
                        // RDKit❌❌:             ...
                        // RDKit❌❌:             totalVdWEnergy += vdWEnergy;
                        // RDKit❌❌:           }
                        // Verbosity-dependent reporting is not modeled yet.
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         bool includeCharge =
                // RDKit✔️✔️:             mmffMolProperties->getMMFFEleTerm() &&
                // RDKit✔️✔️:             !(isDoubleZero(mmffMolProperties->getMMFFPartialCharge(i)) ||
                // RDKit✔️✔️:               isDoubleZero(mmffMolProperties->getMMFFPartialCharge(j)));
                let partial_charge_i = mmff_mol_properties.get_mmff_partial_charge(i)?;
                let partial_charge_j = mmff_mol_properties.get_mmff_partial_charge(j)?;
                let include_charge = mmff_mol_properties.ele_term
                    && !(is_mmff_double_zero(partial_charge_i)
                        || is_mmff_double_zero(partial_charge_j));
                // RDKit✔️✔️:         double chargeTerm = 0.0;
                let mut charge_term = 0.0;
                // RDKit✔️✔️:         if (includeCharge) {
                if include_charge {
                    // RDKit✔️✔️:           chargeTerm = mmffMolProperties->getMMFFPartialCharge(i) *
                    // RDKit✔️✔️:                        mmffMolProperties->getMMFFPartialCharge(j) /
                    // RDKit✔️✔️:                        mmffMolProperties->getMMFFDielectricConstant();
                    charge_term = partial_charge_i * partial_charge_j
                        / mmff_mol_properties.dielectric_constant;
                    // RDKit✔️✔️:           hasContrib = true;
                    has_contrib = true;
                    // RDKit❌❌:           if (mmffMolProperties->getMMFFVerbosity()) {
                    // RDKit❌❌:             ...
                    // RDKit❌❌:             totalEleEnergy += eleEnergy;
                    // RDKit❌❌:           }
                    // Verbosity-dependent reporting is not modeled yet.
                }
                // RDKit✔️✔️:         contrib->addTerm(i, j, vdwConstants, includeCharge, chargeTerm,
                // RDKit✔️✔️:                          mmffMolProperties->getMMFFDielectricModel(), is1_4);
                contrib.add_term(
                    i,
                    j,
                    vdw_constants.as_ref(),
                    include_charge,
                    charge_term,
                    mmff_mol_properties.dielectric_model,
                    is_1_4,
                );
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit❌❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❌❌:     ...
    // RDKit❌❌:   }
    // COSMolKit does not currently model RDKit MMFF summary stream output.
    // RDKit✔️✔️:   if (hasContrib) {
    if has_contrib {
        // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
        field.add_contrib(Box::new(contrib));
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    Ok(())
}

pub(crate) fn add_torsions(
    mol: &Molecule,
    mmff_mol_properties: &MmffMolProperties,
    field: &mut ForceField,
    torsion_bond_smarts: &str,
) -> Result<(), MmffBuilderError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::Tools::addTorsions (Builder.cpp:568-704)
    // RDKit❗❌: void addTorsions(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
    // RDKit❗❌:                  ForceFields::ForceField *field,
    // RDKit❗❌:                  const std::string &torsionBondSmarts) {
    // RDKit✔️✔️:   PRECONDITION(field, "bad ForceField");
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
    // Rust references reproduce RDKit's non-null force-field and
    // MMFFMolProperties preconditions.
    // RDKit✔️✔️:   PRECONDITION(mmffMolProperties->isValid(),
    // RDKit✔️✔️:                "missing atom types - invalid force-field");
    assert!(
        mmff_mol_properties.is_valid(),
        "missing atom types - invalid force-field"
    );

    // RDKit❗❌:   std::ostream &oStream = mmffMolProperties->getMMFFOStream();
    // RDKit❗❌:   ROMol::ADJ_ITER nbr1Idx;
    // RDKit❗❌:   ROMol::ADJ_ITER end1Nbrs;
    // RDKit❗❌:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit❗❌:   ROMol::ADJ_ITER end2Nbrs;
    // RDKit❗❌:   double totalTorsionEnergy = 0.0;
    // RDKit❗❌:   RDGeom::PointPtrVect points;
    // RDKit❗❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❗❌:     ...
    // RDKit❗❌:     points = field->positions();
    // RDKit❗❌:   }
    // COSMolKit does not currently model RDKit MMFF verbosity stream output.

    let torsion_bond_matches = match_torsion_bonds(mol, torsion_bond_smarts)?;

    // RDKit✔️✔️:   auto contrib = std::make_unique<TorsionAngleContrib>(field);
    let mut contrib = TorsionAngleContrib::new(field);
    // RDKit✔️✔️:   bool contribAdded = false;
    let mut contrib_added = false;

    // RDKit✔️✔️:   for (unsigned int i = 0; i < nHits; ++i) {
    for matched in &torsion_bond_matches {
        let idx2 = matched.begin_atom_index;
        let idx3 = matched.end_atom_index;
        let bond_idx = matched.bond_index;
        // RDKit✔️✔️:     MatchVectType match = matchVect[i];
        // RDKit✔️✔️:     TEST_ASSERT(match.size() == 2);
        // RDKit✔️✔️:     int idx2 = match[0].second;
        // RDKit✔️✔️:     int idx3 = match[1].second;
        // RDKit✔️✔️:     const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
        let bond = &mol.bonds()[bond_idx];
        let j_atom = &mol.atoms()[idx2];
        let k_atom = &mol.atoms()[idx3];
        // RDKit✔️✔️:     if (((jAtom->getHybridization() == Atom::SP2) ||
        // RDKit✔️✔️:          (jAtom->getHybridization() == Atom::SP3)) &&
        // RDKit✔️✔️:         ((kAtom->getHybridization() == Atom::SP2) ||
        // RDKit✔️✔️:          (kAtom->getHybridization() == Atom::SP3))) {
        if (j_atom.hybridization() == Hybridization::Sp2
            || j_atom.hybridization() == Hybridization::Sp3)
            && (k_atom.hybridization() == Hybridization::Sp2
                || k_atom.hybridization() == Hybridization::Sp3)
        {
            // RDKit✔️✔️:       ROMol::OEDGE_ITER beg1, end1;
            // RDKit✔️✔️:       boost::tie(beg1, end1) = mol.getAtomBonds(jAtom);
            // RDKit✔️✔️:       while (beg1 != end1) {
            for first_neighbor in mol.topology_block().adjacency.neighbors_of(idx2) {
                let t_bond1 = &mol.bonds()[first_neighbor.bond.index()];
                // RDKit✔️✔️:         const Bond *tBond1 = mol[*beg1];
                // RDKit✔️✔️:         if (tBond1 != bond) {
                if t_bond1.id() != bond.id() {
                    // RDKit✔️✔️:           int idx1 = tBond1->getOtherAtomIdx(idx2);
                    let idx1 = first_neighbor.atom_index;
                    // RDKit✔️✔️:           ROMol::OEDGE_ITER beg2, end2;
                    // RDKit✔️✔️:           boost::tie(beg2, end2) = mol.getAtomBonds(kAtom);
                    // RDKit✔️✔️:           while (beg2 != end2) {
                    for second_neighbor in mol.topology_block().adjacency.neighbors_of(idx3) {
                        let t_bond2 = &mol.bonds()[second_neighbor.bond.index()];
                        // RDKit✔️✔️:             const Bond *tBond2 = mol[*beg2];
                        // RDKit✔️✔️:             if ((tBond2 != bond) && (tBond2 != tBond1)) {
                        if t_bond2.id() != bond.id() && t_bond2.id() != t_bond1.id() {
                            // RDKit✔️✔️:               int idx4 = tBond2->getOtherAtomIdx(idx3);
                            let idx4 = second_neighbor.atom_index;
                            // RDKit✔️✔️:               // make sure this isn't a three-membered ring:
                            // RDKit✔️✔️:               if (idx4 != idx1) {
                            if idx4 != idx1 {
                                // RDKit✔️✔️:                 unsigned int torType;
                                // RDKit✔️✔️:                 MMFFTor mmffTorParams;
                                // RDKit✔️✔️:                 if (mmffMolProperties->getMMFFTorsionParams(
                                // RDKit✔️✔️:                         mol, idx1, idx2, idx3, idx4, torType, mmffTorParams)) {
                                if let Some((_tor_type, mmff_tor_params)) = mmff_mol_properties
                                    .get_mmff_torsion_params(idx1, idx2, idx3, idx4)?
                                {
                                    // RDKit✔️✔️:                   contrib->addTerm(idx1, idx2, idx3, idx4, &mmffTorParams);
                                    contrib.add_term(idx1, idx2, idx3, idx4, &mmff_tor_params);
                                    // RDKit✔️✔️:                   contribAdded = true;
                                    contrib_added = true;
                                    // RDKit❗❌:                   if (mmffMolProperties->getMMFFVerbosity()) {
                                    // RDKit❗❌:                     ...
                                    // RDKit❗❌:                     totalTorsionEnergy += torsionEnergy;
                                    // RDKit❗❌:                   }
                                    // Verbosity-dependent reporting is not modeled yet.
                                }
                            }
                        }
                    }
                }
            }
        }
        // RDKit✔️✔️:     }
    }
    // RDKit❗❌:   if (mmffMolProperties->getMMFFVerbosity()) {
    // RDKit❗❌:     if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    // RDKit❗❌:       oStream << std::endl;
    // RDKit❗❌:     }
    // RDKit❗❌:     oStream << "TOTAL TORSIONAL ENERGY         =" << std::right << std::setw(16)
    // RDKit❗❌:             << std::fixed << std::setprecision(4) << totalTorsionEnergy
    // RDKit❗❌:             << std::endl;
    // RDKit❗❌:   }
    // COSMolKit does not currently model RDKit MMFF summary stream output.
    // RDKit✔️✔️:   if (contribAdded) {
    if contrib_added {
        // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib.release()));
        field.add_contrib(Box::new(contrib));
        // RDKit✔️✔️:   }
    }
    // RDKit❗❌: }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        MmffBuilderError, RELATION_1_2, RELATION_1_3, RELATION_1_4, RELATION_1_X, add_angles,
        add_bonds, add_nonbonded, add_oop, add_stretch_bend, add_torsions, build_neighbor_matrix,
        construct_force_field, construct_force_field_with_props, get_two_bit_cell,
        set_two_bit_cell, two_bit_cell_pos,
    };
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
    use crate::chemistry::forcefield::mmff::mol_properties::{
        MMFF_DIELECTRIC_CONSTANT, MMFF_DIELECTRIC_DISTANCE, MmffAtomProperties, MmffMolProperties,
        MmffVariant,
    };
    use crate::{
        AromaticityAssignment, AtomSpec, BondOrder, BondSpec, Conformer3D, Element, Hybridization,
        Molecule, MoleculeBuilder,
    };

    fn force_field_with_positions(positions: &[[f64; 3]]) -> ForceField {
        let mut force_field = ForceField::new(3);
        force_field.positions_mut().extend(
            positions
                .iter()
                .map(|coords| ForceFieldVec3::new(coords[0], coords[1], coords[2])),
        );
        force_field.initialize();
        force_field
    }

    fn mmff_props_for_molecule_and_atom_types(
        molecule: Molecule,
        atom_types: &[u8],
    ) -> MmffMolProperties {
        let num_bonds = molecule.num_bonds();
        MmffMolProperties {
            molecule,
            valid: true,
            variant: MmffVariant::Mmff94,
            bond_term: true,
            angle_term: true,
            stretch_bend_term: true,
            oop_term: true,
            torsion_term: true,
            vdw_term: true,
            ele_term: true,
            dielectric_constant: 1.0,
            dielectric_model: MMFF_DIELECTRIC_CONSTANT,
            verbosity: 0,
            atom_properties: atom_types
                .iter()
                .copied()
                .map(|atom_type| MmffAtomProperties {
                    atom_type,
                    formal_charge: 0.0,
                    partial_charge: 0.0,
                })
                .collect(),
            aromaticity: AromaticityAssignment {
                atom_aromatic: Vec::new(),
                bond_aromatic: vec![false; num_bonds],
                aromatic_ring_count: 0,
            },
        }
    }

    fn two_atom_molecule(first: Element, second: Element, order: BondOrder) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let first = builder.add_atom(AtomSpec::new(first));
        let second = builder.add_atom(AtomSpec::new(second));
        builder
            .add_bond(BondSpec::new(first, second, order))
            .expect("test molecule bond endpoints are valid");
        builder.build().expect("test molecule is valid")
    }

    fn linear_three_atom_molecule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::H));
        let a2 = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder.build().expect("test molecule is valid")
    }

    fn three_atom_angle_molecule(first: Element, center: Element, third: Element) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(first));
        let a1 = builder.add_atom(AtomSpec::new(center));
        let a2 = builder.add_atom(AtomSpec::new(third));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder.build().expect("test molecule is valid")
    }

    fn three_neighbor_center_molecule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::H));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::H));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder.build().expect("test molecule is valid")
    }

    fn source_bond_stretch_energy(r0: f64, kb: f64, distance: f64) -> f64 {
        let dist_term = distance - r0;
        let dist_term2 = dist_term * dist_term;
        let c1 = 143.9325;
        let cs = -2.0;
        let c3 = 7.0 / 12.0;
        0.5 * c1 * kb * dist_term2 * (1.0 + cs * dist_term + c3 * cs * cs * dist_term2)
    }

    fn source_angle_bend_energy(theta0: f64, ka: f64, is_linear: bool, cos_theta: f64) -> f64 {
        let angle = 180.0 / std::f64::consts::PI * cos_theta.acos() - theta0;
        let cb = -0.006981317;
        let c2 = 143.9325 * (std::f64::consts::PI / 180.0) * (std::f64::consts::PI / 180.0);
        if is_linear {
            143.9325 * ka * (1.0 + cos_theta)
        } else {
            0.5 * c2 * ka * angle * angle * (1.0 + cb * angle)
        }
    }

    fn source_stretch_bend_energy(
        delta_dist1: f64,
        delta_dist2: f64,
        delta_theta: f64,
        force_constants: (f64, f64),
    ) -> f64 {
        let factor = 143.9325 * (std::f64::consts::PI / 180.0) * delta_theta;
        factor * force_constants.0 * delta_dist1 + factor * force_constants.1 * delta_dist2
    }

    fn source_oop_chi(
        i_point: ForceFieldVec3,
        j_point: ForceFieldVec3,
        k_point: ForceFieldVec3,
        l_point: ForceFieldVec3,
    ) -> f64 {
        let mut r_ji = i_point - j_point;
        let mut r_jk = k_point - j_point;
        let mut r_jl = l_point - j_point;
        r_ji /= r_ji.length();
        r_jk /= r_jk.length();
        r_jl /= r_jl.length();

        let mut n = r_ji.cross_product(r_jk);
        n /= n.length();
        let mut sin_chi = n.dot_product(r_jl);
        sin_chi = sin_chi.clamp(-1.0, 1.0);

        (180.0 / std::f64::consts::PI) * sin_chi.asin()
    }

    fn source_oop_bend_energy(chi: f64, koop: f64) -> f64 {
        let c2 = 143.9325 * (std::f64::consts::PI / 180.0) * (std::f64::consts::PI / 180.0);
        0.5 * c2 * koop * chi * chi
    }

    fn source_torsion_energy(v1: f64, v2: f64, v3: f64, cos_phi: f64) -> f64 {
        let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
        let cos3_phi = cos_phi * (2.0 * cos2_phi - 1.0);
        0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2_phi) + v3 * (1.0 + cos3_phi))
    }

    fn source_vdw_energy(dist: f64, r_star_ij: f64, well_depth: f64) -> f64 {
        let vdw1 = 1.07;
        let vdw1m1 = vdw1 - 1.0;
        let vdw2 = 1.12;
        let vdw2m1 = vdw2 - 1.0;
        let dist2 = dist * dist;
        let dist7 = dist2 * dist2 * dist2 * dist;
        let a_term = vdw1 * r_star_ij / (dist + vdw1m1 * r_star_ij);
        let a_term2 = a_term * a_term;
        let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
        let r_star_ij2 = r_star_ij * r_star_ij;
        let r_star_ij7 = r_star_ij2 * r_star_ij2 * r_star_ij2 * r_star_ij;
        let b_term = vdw2 * r_star_ij7 / (dist7 + vdw2m1 * r_star_ij7) - 2.0;
        well_depth * a_term7 * b_term
    }

    fn source_ele_energy(dist: f64, charge_term: f64, dielectric_model: u8, is_1_4: bool) -> f64 {
        let mut corr_dist = dist + 0.05;
        let diel = 332.0716;
        let sc1_4 = 0.75;
        if dielectric_model == MMFF_DIELECTRIC_DISTANCE {
            corr_dist *= corr_dist;
        }
        diel * charge_term / corr_dist * if is_1_4 { sc1_4 } else { 1.0 }
    }

    fn source_torsion_cos_phi(
        i_point: ForceFieldVec3,
        j_point: ForceFieldVec3,
        k_point: ForceFieldVec3,
        l_point: ForceFieldVec3,
    ) -> f64 {
        let r1 = i_point - j_point;
        let r2 = k_point - j_point;
        let r3 = j_point - k_point;
        let r4 = l_point - k_point;
        let t1 = r1.cross_product(r2);
        let t2 = r3.cross_product(r4);
        let t1_len = t1.length();
        let t2_len = t2.length();
        if t1_len < 1.0e-10 || t2_len < 1.0e-10 {
            return 0.0;
        }
        (t1.dot_product(t2) / (t1_len * t2_len)).clamp(-1.0, 1.0)
    }

    fn torsion_chain_molecule(
        atoms: [(Element, Hybridization); 4],
        bonds: [BondOrder; 3],
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(atoms[0].0).with_hybridization(atoms[0].1));
        let a1 = builder.add_atom(AtomSpec::new(atoms[1].0).with_hybridization(atoms[1].1));
        let a2 = builder.add_atom(AtomSpec::new(atoms[2].0).with_hybridization(atoms[2].1));
        let a3 = builder.add_atom(AtomSpec::new(atoms[3].0).with_hybridization(atoms[3].1));
        for (begin, end, order) in [(a0, a1, bonds[0]), (a1, a2, bonds[1]), (a2, a3, bonds[2])] {
            builder
                .add_bond(BondSpec::new(begin, end, order))
                .expect("test molecule bond endpoints are valid");
        }
        builder.build().expect("test molecule is valid")
    }

    fn triangle_ring_molecule() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        let a1 = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        let a2 = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a0)] {
            builder
                .add_bond(BondSpec::new(begin, end, BondOrder::Single))
                .expect("test molecule bond endpoints are valid");
        }
        builder.build().expect("test molecule is valid")
    }

    fn linear_four_atom_molecule(elements: [Element; 4]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(elements[0]));
        let a1 = builder.add_atom(AtomSpec::new(elements[1]));
        let a2 = builder.add_atom(AtomSpec::new(elements[2]));
        let a3 = builder.add_atom(AtomSpec::new(elements[3]));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a3)] {
            builder
                .add_bond(BondSpec::new(begin, end, BondOrder::Single))
                .expect("test molecule bond endpoints are valid");
        }
        builder.build().expect("test molecule is valid")
    }

    fn disconnected_two_fragment_chain() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::H));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder.build().expect("test molecule is valid")
    }

    fn disconnected_two_atom_fragments_with_conformer() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::H));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Single))
            .expect("test molecule bond endpoints are valid");
        builder
            .add_3d_conformer(vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [5.2, 0.0, 0.0],
            ])
            .expect("3D conformer row count should match");
        builder.build().expect("test molecule is valid")
    }

    fn linear_four_atom_molecule_with_two_conformers() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::H));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::H));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a3)] {
            builder
                .add_bond(BondSpec::new(begin, end, BondOrder::Single))
                .expect("test molecule bond endpoints are valid");
        }
        builder
            .add_3d_conformer(vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.2, 0.0, 0.0],
            ])
            .expect("first 3D conformer row count should match");
        builder
            .add_3d_conformer(vec![
                [10.0, 0.0, 0.0],
                [11.0, 0.0, 0.0],
                [12.0, 0.0, 0.0],
                [13.2, 0.0, 0.0],
            ])
            .expect("second 3D conformer row count should match");
        builder.build().expect("test molecule is valid")
    }

    fn empty_molecule_with_named_conformers() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder
            .add_conformer(Conformer3D::new(7, Vec::new(), true))
            .expect("empty 3D conformer row count should match");
        builder
            .add_conformer(Conformer3D::new(9, Vec::new(), true))
            .expect("empty 3D conformer row count should match");
        builder.build().expect("empty test molecule is valid")
    }

    fn manual_neighbor_matrix(n_atoms: usize, pairs: &[(usize, usize, u8)]) -> Vec<u8> {
        let mut res = vec![RELATION_1_X; n_atoms * (n_atoms + 1) / 2];
        for atom_idx in 0..n_atoms {
            set_two_bit_cell(
                &mut res,
                two_bit_cell_pos(n_atoms, atom_idx, atom_idx),
                RELATION_1_X,
            );
        }
        for &(i, j, relation) in pairs {
            set_two_bit_cell(&mut res, two_bit_cell_pos(n_atoms, i, j), relation);
        }
        res
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "actual={actual}, expected={expected}"
        );
    }

    #[test]
    fn mmff_builder_neighbor_matrix_distinguishes_all_topological_relations() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let matrix = build_neighbor_matrix(&mol);

        assert_eq!(matrix.len(), mol.num_atoms() * (mol.num_atoms() + 1) / 2);
        assert_eq!(
            get_two_bit_cell(&matrix, two_bit_cell_pos(mol.num_atoms(), 0, 1)),
            RELATION_1_2
        );
        assert_eq!(
            get_two_bit_cell(&matrix, two_bit_cell_pos(mol.num_atoms(), 0, 2)),
            RELATION_1_3
        );
        assert_eq!(
            get_two_bit_cell(&matrix, two_bit_cell_pos(mol.num_atoms(), 0, 3)),
            RELATION_1_4
        );
        assert_eq!(
            get_two_bit_cell(&matrix, two_bit_cell_pos(mol.num_atoms(), 0, 0)),
            RELATION_1_X
        );
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_builds_from_default_conformer() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);

        let ff = construct_force_field_with_props(&mol, &props, 100.0, -1, false)
            .expect("default conformer should build MMFF force field");

        assert_eq!(ff.positions().len(), 4);
        assert_eq!(ff.positions()[0], ForceFieldVec3::new(0.0, 0.0, 0.0));
        assert_eq!(ff.positions()[3], ForceFieldVec3::new(3.2, 0.0, 0.0));
        assert_eq!(ff.contribs().len(), 4);
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_selects_requested_conformer_id() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);

        let ff = construct_force_field_with_props(&mol, &props, 100.0, 1, false)
            .expect("explicit conformer id should select matching conformer");

        assert_eq!(ff.positions().len(), 4);
        assert_eq!(ff.positions()[0], ForceFieldVec3::new(10.0, 0.0, 0.0));
        assert_eq!(ff.positions()[3], ForceFieldVec3::new(13.2, 0.0, 0.0));
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_respects_term_toggles() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.bond_term = false;
        props.angle_term = false;
        props.stretch_bend_term = false;
        props.oop_term = false;
        props.torsion_term = false;
        props.vdw_term = false;
        props.ele_term = false;

        let ff = construct_force_field_with_props(&mol, &props, 100.0, -1, false)
            .expect("disabling all terms should still construct a force field");

        assert_eq!(ff.positions().len(), 4);
        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_respects_nonbonded_threshold() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.bond_term = false;
        props.angle_term = false;
        props.stretch_bend_term = false;
        props.oop_term = false;
        props.torsion_term = false;
        props.atom_properties[0].partial_charge = 0.25;
        props.atom_properties[3].partial_charge = -0.5;

        let ff = construct_force_field_with_props(&mol, &props, 3.0, -1, false)
            .expect("threshold miss should still build force field");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_respects_fragment_filter() {
        let mol = disconnected_two_atom_fragments_with_conformer();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5, 1, 5]);
        props.bond_term = false;
        props.angle_term = false;
        props.stretch_bend_term = false;
        props.oop_term = false;
        props.torsion_term = false;
        props.atom_properties[0].partial_charge = 0.2;
        props.atom_properties[3].partial_charge = -0.4;

        let ff = construct_force_field_with_props(&mol, &props, 100.0, -1, true)
            .expect("interfragment filter should skip disconnected nonbonded pair");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_rejects_missing_3d_conformer() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);

        let err = match construct_force_field_with_props(&mol, &props, 100.0, -1, false) {
            Ok(_) => panic!("missing 3D conformer should fail explicitly"),
            Err(err) => err,
        };

        assert_eq!(err, MmffBuilderError::Missing3dConformer { conf_id: -1 });
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_rejects_invalid_negative_conformer_id() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);

        let err = match construct_force_field_with_props(&mol, &props, 100.0, -2, false) {
            Ok(_) => panic!("conf_id below -1 should fail"),
            Err(err) => err,
        };

        assert_eq!(err, MmffBuilderError::InvalidConformerId { conf_id: -2 });
    }

    #[test]
    fn mmff_builder_construct_force_field_with_props_rejects_unknown_conformer_id() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);

        let err = match construct_force_field_with_props(&mol, &props, 100.0, 7, false) {
            Ok(_) => panic!("unknown conformer id should fail explicitly"),
            Err(err) => err,
        };

        assert_eq!(err, MmffBuilderError::Missing3dConformer { conf_id: 7 });
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_construct_force_field_with_props_rejects_invalid_mmff_properties() {
        let mol = linear_four_atom_molecule_with_two_conformers();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.valid = false;

        let _ = construct_force_field_with_props(&mol, &props, 100.0, -1, false);
    }

    #[test]
    fn mmff_builder_construct_force_field_supports_ported_empirical_rules_after_atom_typing() {
        let mol = linear_four_atom_molecule_with_two_conformers();

        let ff = construct_force_field(&mol, 100.0, -1, false)
            .expect("ported empirical rules should construct the MMFF force field");

        assert_eq!(ff.positions().len(), mol.num_atoms());
        assert!(!ff.contribs().is_empty());
    }

    #[test]
    fn mmff_builder_construct_force_field_rejects_missing_3d_conformer_for_empty_molecule() {
        let mol = Molecule::new();

        let err = match construct_force_field(&mol, 100.0, -1, false) {
            Ok(_) => panic!("empty molecule without 3D conformer should fail explicitly"),
            Err(err) => err,
        };

        assert_eq!(err, MmffBuilderError::Missing3dConformer { conf_id: -1 });
    }

    #[test]
    fn mmff_builder_construct_force_field_rejects_invalid_negative_conformer_id() {
        let mol = empty_molecule_with_named_conformers();

        let err = match construct_force_field(&mol, 100.0, -2, false) {
            Ok(_) => panic!("conf_id below -1 should fail explicitly"),
            Err(err) => err,
        };

        assert_eq!(err, MmffBuilderError::InvalidConformerId { conf_id: -2 });
    }

    #[test]
    fn mmff_builder_construct_force_field_rejects_unknown_conformer_id() {
        let mol = empty_molecule_with_named_conformers();

        let err = match construct_force_field(&mol, 100.0, 8, false) {
            Ok(_) => panic!("unknown conformer id should fail explicitly"),
            Err(err) => err,
        };

        assert_eq!(err, MmffBuilderError::Missing3dConformer { conf_id: 8 });
    }

    #[test]
    fn mmff_builder_add_bonds_adds_single_bond_stretch_contrib() {
        let mol = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5]);
        let mut ff = force_field_with_positions(&[[0.0, 0.0, 0.0], [1.343, 0.0, 0.0]]);

        add_bonds(&mol, &props, &mut ff).expect("tabulated MMFF bond should add contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.343, 0.0, 0.0];
        let expected = source_bond_stretch_energy(1.093, 4.766, 1.343);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_bonds_keeps_force_field_empty_for_molecule_without_bonds() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.add_atom(AtomSpec::new(Element::H));
        let mol = builder.build().expect("test molecule is valid");
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5]);
        let mut ff = force_field_with_positions(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);

        add_bonds(&mol, &props, &mut ff).expect("molecule without bonds should be a no-op");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_bonds_accumulates_multiple_bonds_into_one_contrib() {
        let mol = linear_three_atom_molecule();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5, 5]);
        let mut ff =
            force_field_with_positions(&[[0.0, 0.0, 0.0], [1.343, 0.0, 0.0], [-1.293, 0.0, 0.0]]);

        add_bonds(&mol, &props, &mut ff).expect("multiple tabulated bonds should add one contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.343, 0.0, 0.0, -1.293, 0.0, 0.0];
        let expected = source_bond_stretch_energy(1.093, 4.766, 1.343)
            + source_bond_stretch_energy(1.093, 4.766, 1.293);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_add_bonds_rejects_invalid_mmff_properties() {
        let mol = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5]);
        props.valid = false;
        let mut ff = force_field_with_positions(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);

        let _ = add_bonds(&mol, &props, &mut ff);
    }

    #[test]
    fn mmff_builder_add_bonds_adds_empirical_bond_stretch_contrib() {
        let mol = two_atom_molecule(Element::C, Element::O, BondOrder::Single);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[60, 6]);
        let mut ff = force_field_with_positions(&[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]);

        add_bonds(&mol, &props, &mut ff).expect("empirical bond parameters should add a term");

        assert_eq!(ff.contribs().len(), 1);
        let expected = source_bond_stretch_energy(1.405, 5.129115902527102, 1.2);
        assert_close(
            ff.contribs()[0].get_energy(&[0.0, 0.0, 0.0, 1.2, 0.0, 0.0]),
            expected,
        );
    }

    #[test]
    fn mmff_builder_add_angles_adds_single_angle_bend_contrib() {
        let mol = three_atom_angle_molecule(Element::H, Element::C, Element::C);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1]);
        let mut ff =
            force_field_with_positions(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        add_angles(&mol, &props, &mut ff).expect("tabulated MMFF angle should add contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let expected = source_angle_bend_energy(110.549, 0.636, false, 0.0);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_angles_keeps_force_field_empty_without_angle_terms() {
        let mol = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5]);
        let mut ff = force_field_with_positions(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);

        add_angles(&mol, &props, &mut ff).expect("molecule without angles should be a no-op");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_angles_accumulates_multiple_angles_into_one_contrib() {
        let mol = three_neighbor_center_molecule();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 5, 1]);
        let mut ff = force_field_with_positions(&[
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        add_angles(&mol, &props, &mut ff)
            .expect("multiple tabulated angles should add one contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let hch = props
            .get_mmff_angle_bend_params(0, 1, 2)
            .expect("H-C-H lookup should parse")
            .expect("H-C-H angle should be tabulated")
            .1;
        let hcc = props
            .get_mmff_angle_bend_params(0, 1, 3)
            .expect("H-C-C lookup should parse")
            .expect("H-C-C angle should be tabulated")
            .1;
        let expected = source_angle_bend_energy(hch.theta0, hch.ka, false, 0.0)
            + source_angle_bend_energy(hcc.theta0, hcc.ka, false, 0.0)
            + source_angle_bend_energy(hcc.theta0, hcc.ka, false, 0.0);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_add_angles_rejects_invalid_mmff_properties() {
        let mol = three_atom_angle_molecule(Element::H, Element::C, Element::C);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1]);
        props.valid = false;
        let mut ff =
            force_field_with_positions(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        let _ = add_angles(&mol, &props, &mut ff);
    }

    #[test]
    fn mmff_builder_add_angles_adds_empirical_angle_contrib() {
        let mol = three_atom_angle_molecule(Element::F, Element::C, Element::CL);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[11, 1, 12]);
        let mut ff =
            force_field_with_positions(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        add_angles(&mol, &props, &mut ff).expect("empirical angle should add a contribution");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let expected = source_angle_bend_energy(108.9, 1.2566039721725888, false, 0.0);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_stretch_bend_adds_single_stretch_bend_contrib() {
        let mol = three_atom_angle_molecule(Element::H, Element::C, Element::C);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1]);
        let mut ff =
            force_field_with_positions(&[[1.1, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.2, 0.0]]);

        add_stretch_bend(&mol, &props, &mut ff)
            .expect("tabulated MMFF stretch-bend should add contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0];
        let theta = 90.0;
        let expected =
            source_stretch_bend_energy(1.1 - 1.093, 1.2 - 1.508, theta - 110.549, (0.070, 0.227));
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_stretch_bend_keeps_force_field_empty_without_terms() {
        let mol = two_atom_molecule(Element::C, Element::H, BondOrder::Single);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5]);
        let mut ff = force_field_with_positions(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);

        add_stretch_bend(&mol, &props, &mut ff)
            .expect("molecule without stretch-bend terms should be a no-op");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_stretch_bend_accumulates_multiple_terms_into_one_contrib() {
        let mol = three_neighbor_center_molecule();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 5, 1]);
        let mut ff = force_field_with_positions(&[
            [1.1, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.2, 0.0],
            [0.0, 0.0, 1.6],
        ]);

        add_stretch_bend(&mol, &props, &mut ff)
            .expect("multiple tabulated stretch-bends should add one contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, 1.6];
        let hch = props
            .get_mmff_stretch_bend_params(0, 1, 2)
            .expect("H-C-H stretch-bend lookup should parse")
            .expect("H-C-H stretch-bend should be tabulated");
        let hcc = props
            .get_mmff_stretch_bend_params(0, 1, 3)
            .expect("H-C-C stretch-bend lookup should parse")
            .expect("H-C-C stretch-bend should be tabulated");
        let expected = source_stretch_bend_energy(
            1.1 - hch.2[0].r0,
            1.2 - hch.2[1].r0,
            90.0 - hch.3.theta0,
            (hch.1.kba_ijk, hch.1.kba_kji),
        ) + source_stretch_bend_energy(
            1.1 - hcc.2[0].r0,
            1.6 - hcc.2[1].r0,
            90.0 - hcc.3.theta0,
            (hcc.1.kba_ijk, hcc.1.kba_kji),
        ) + source_stretch_bend_energy(
            1.2 - hcc.2[0].r0,
            1.6 - hcc.2[1].r0,
            90.0 - hcc.3.theta0,
            (hcc.1.kba_ijk, hcc.1.kba_kji),
        );
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_stretch_bend_skips_linear_central_atom() {
        let mol = three_atom_angle_molecule(Element::C, Element::C, Element::C);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 4, 1]);
        let mut ff =
            force_field_with_positions(&[[1.2, 0.0, 0.0], [0.0, 0.0, 0.0], [-1.2, 0.0, 0.0]]);

        add_stretch_bend(&mol, &props, &mut ff)
            .expect("linear central atom should skip stretch-bend addition");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_add_stretch_bend_rejects_invalid_mmff_properties() {
        let mol = three_atom_angle_molecule(Element::H, Element::C, Element::C);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1]);
        props.valid = false;
        let mut ff =
            force_field_with_positions(&[[1.1, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.2, 0.0]]);

        let _ = add_stretch_bend(&mol, &props, &mut ff);
    }

    #[test]
    fn mmff_builder_add_stretch_bend_adds_empirical_angle_contrib() {
        let mol = three_atom_angle_molecule(Element::F, Element::C, Element::CL);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[11, 1, 12]);
        let mut ff =
            force_field_with_positions(&[[1.36, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.773, 0.0]]);

        add_stretch_bend(&mol, &props, &mut ff)
            .expect("stretch-bend should consume empirical angle parameters");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.36, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.773, 0.0];
        assert_close(ff.contribs()[0].get_energy(&pos), 0.0);
    }

    #[test]
    fn mmff_builder_add_oop_adds_single_oop_bend_contrib() {
        let mol = three_neighbor_center_molecule();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 2, 1, 2]);
        let positions = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_oop(&mol, &props, &mut ff).expect("tabulated MMFF OOP term should add contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let params = props
            .get_mmff_oop_bend_params(0, 1, 2, 3)
            .expect("OOP lookup should parse")
            .expect("OOP params should be tabulated");
        let expected = source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 1.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
            ),
            params.koop,
        ) + source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
                ForceFieldVec3::new(0.0, 1.0, 0.0),
            ),
            params.koop,
        ) + source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(0.0, 1.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
                ForceFieldVec3::new(1.0, 0.0, 0.0),
            ),
            params.koop,
        );
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_oop_keeps_force_field_empty_without_terms() {
        let mol = three_neighbor_center_molecule();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 1]);
        let mut ff = force_field_with_positions(&[
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        add_oop(&mol, &props, &mut ff).expect("missing OOP parameters should be a no-op");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_oop_skips_non_trigonal_center() {
        let mol = three_atom_angle_molecule(Element::H, Element::C, Element::C);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1]);
        let mut ff =
            force_field_with_positions(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);

        add_oop(&mol, &props, &mut ff).expect("non-degree-3 center should skip OOP addition");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_add_oop_rejects_invalid_mmff_properties() {
        let mol = three_neighbor_center_molecule();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 2, 1, 2]);
        props.valid = false;
        let mut ff = force_field_with_positions(&[
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        let _ = add_oop(&mol, &props, &mut ff);
    }

    #[test]
    fn mmff_builder_add_oop_uses_mmff94s_parameters() {
        let mol = three_neighbor_center_molecule();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 10, 1, 1]);
        props.variant = MmffVariant::Mmff94s;
        let positions = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_oop(&mol, &props, &mut ff).expect("MMFF94s OOP row should add contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let params = props
            .get_mmff_oop_bend_params(0, 1, 2, 3)
            .expect("MMFF94s OOP lookup should parse")
            .expect("MMFF94s OOP params should be tabulated");
        let expected = source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 1.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
            ),
            params.koop,
        ) + source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
                ForceFieldVec3::new(0.0, 1.0, 0.0),
            ),
            params.koop,
        ) + source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(0.0, 1.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
                ForceFieldVec3::new(1.0, 0.0, 0.0),
            ),
            params.koop,
        );
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_torsions_accepts_custom_bond_smarts() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 1]);
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);

        add_torsions(&mol, &props, &mut ff, "[*:1]-[*:2]")
            .expect("custom torsion SMARTS should select the three single bonds");

        assert_eq!(ff.contribs().len(), 1);
    }

    #[test]
    fn mmff_builder_add_nonbonded_adds_single_combined_term_for_1_4_pair() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.atom_properties[0].partial_charge = 0.25;
        props.atom_properties[3].partial_charge = -0.5;
        props.dielectric_constant = 2.0;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_2),
                (1, 3, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_4),
            ],
        );
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.2, 0.0, 0.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false)
            .expect("1-4 pair with vdW and electrostatics should add one contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.2, 0.0, 0.0];
        let vdw = props
            .get_mmff_vdw_params(0, 3)
            .expect("vdW lookup should parse")
            .expect("H-H vdW constants should exist");
        let charge_term = 0.25 * -0.5 / 2.0;
        let expected = source_vdw_energy(3.2, vdw.r_ij_star, vdw.epsilon)
            + source_ele_energy(3.2, charge_term, MMFF_DIELECTRIC_CONSTANT, true);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_nonbonded_skips_pairs_below_1_4_relation() {
        let mol = linear_three_atom_molecule();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5, 5]);
        props.atom_properties[0].partial_charge = 0.2;
        props.atom_properties[1].partial_charge = -0.3;
        props.atom_properties[2].partial_charge = 0.4;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_2),
                (1, 2, RELATION_1_3),
            ],
        );
        let mut ff =
            force_field_with_positions(&[[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [0.0, 1.2, 0.0]]);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false)
            .expect("1-2 and 1-3 relations should not add nonbonded terms");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_nonbonded_skips_pairs_beyond_threshold() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.atom_properties[0].partial_charge = 0.25;
        props.atom_properties[3].partial_charge = -0.5;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_2),
                (1, 3, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_4),
            ],
        );
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.2, 0.0, 0.0],
        ]);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 3.0, false)
            .expect("distance threshold miss should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_nonbonded_respects_fragment_interaction_filter() {
        let mol = disconnected_two_fragment_chain();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5, 1, 5]);
        props.atom_properties[0].partial_charge = 0.2;
        props.atom_properties[3].partial_charge = -0.4;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_X),
            ],
        );
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.2, 0.0, 0.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, true)
            .expect("interfragment filter should skip disconnected pairs");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_nonbonded_allows_interfragment_pairs_when_not_ignored() {
        let mol = disconnected_two_fragment_chain();
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 5, 1, 5]);
        props.atom_properties[0].partial_charge = 0.2;
        props.atom_properties[3].partial_charge = -0.4;
        props.dielectric_model = MMFF_DIELECTRIC_DISTANCE;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (2, 3, RELATION_1_2),
                (0, 2, RELATION_1_2),
                (1, 2, RELATION_1_3),
                (1, 3, RELATION_1_3),
                (0, 3, RELATION_1_X),
            ],
        );
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [5.2, 0.0, 0.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false)
            .expect("interfragment pairs should be included when not ignored");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 0.0, 5.2, 0.0, 0.0];
        let vdw = props
            .get_mmff_vdw_params(0, 3)
            .expect("vdW lookup should parse")
            .expect("C-H vdW constants should exist");
        let charge_term = 0.2 * -0.4 / props.dielectric_constant;
        let expected = source_vdw_energy(5.2, vdw.r_ij_star, vdw.epsilon)
            + source_ele_energy(5.2, charge_term, MMFF_DIELECTRIC_DISTANCE, false);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_nonbonded_adds_electrostatics_without_vdw_when_vdw_term_disabled() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.vdw_term = false;
        props.atom_properties[0].partial_charge = 0.25;
        props.atom_properties[3].partial_charge = -0.5;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_2),
                (1, 3, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_4),
            ],
        );
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.2, 0.0, 0.0],
        ]);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false)
            .expect("electrostatics-only pair should add a term");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.2, 0.0, 0.0];
        let expected = source_ele_energy(3.2, 0.25 * -0.5, MMFF_DIELECTRIC_CONSTANT, true);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_nonbonded_adds_vdw_without_electrostatics_for_zero_charge_pair() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_2),
                (1, 3, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_4),
            ],
        );
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.2, 0.0, 0.0],
        ]);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false)
            .expect("zero partial charges should still allow vdW-only term");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.2, 0.0, 0.0];
        let vdw = props
            .get_mmff_vdw_params(0, 3)
            .expect("vdW lookup should parse")
            .expect("H-H vdW constants should exist");
        let expected = source_vdw_energy(3.2, vdw.r_ij_star, vdw.epsilon);
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_nonbonded_keeps_force_field_empty_when_both_terms_disabled() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.vdw_term = false;
        props.ele_term = false;
        props.atom_properties[0].partial_charge = 0.25;
        props.atom_properties[3].partial_charge = -0.5;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_2),
                (1, 3, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_4),
            ],
        );
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.2, 0.0, 0.0],
        ]);

        add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false)
            .expect("disabled vdW and electrostatics should be a no-op");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_add_nonbonded_rejects_invalid_mmff_properties() {
        let mol = linear_four_atom_molecule([Element::H, Element::C, Element::C, Element::H]);
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 5]);
        props.valid = false;
        let neighbor_matrix = manual_neighbor_matrix(
            mol.num_atoms(),
            &[
                (0, 1, RELATION_1_2),
                (0, 2, RELATION_1_3),
                (1, 2, RELATION_1_2),
                (1, 3, RELATION_1_3),
                (2, 3, RELATION_1_2),
                (0, 3, RELATION_1_4),
            ],
        );
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.2, 0.0, 0.0],
        ]);

        let _ = add_nonbonded(&mol, &props, &mut ff, &neighbor_matrix, 100.0, false);
    }

    #[test]
    fn mmff_builder_add_torsions_adds_single_torsion_contrib_for_default_query_hit() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 1]);
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("default torsion query hit should add a contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 1.0];
        let params = props
            .get_mmff_torsion_params(0, 1, 2, 3)
            .expect("torsion lookup should parse")
            .expect("torsion params should be tabulated");
        let expected = source_torsion_energy(
            params.1.v1,
            params.1.v2,
            params.1.v3,
            source_torsion_cos_phi(
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(1.0, 1.0, 0.0),
                ForceFieldVec3::new(2.0, 1.0, 1.0),
            ),
        );
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    fn mmff_builder_add_torsions_skips_default_query_hits_with_triple_adjacent_atom() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp2),
                (Element::C, Hybridization::Sp2),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Triple, BondOrder::Single, BondOrder::Single],
        );
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 1]);
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("triple-adjacent default-query miss should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_torsions_skips_non_sp2_sp3_central_atoms() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3d2),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 1]);
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("unsupported central hybridization should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_torsions_skips_three_membered_ring_torsions() {
        let mol = triangle_ring_molecule();
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1]);
        let mut ff =
            force_field_with_positions(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.8, 0.0]]);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("three-membered ring torsions should skip");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_torsions_keeps_force_field_empty_without_terms() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 2, 3]);
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("zero or missing torsion parameters should be a no-op");

        assert_eq!(ff.contribs().len(), 0);
    }

    #[test]
    fn mmff_builder_add_torsions_uses_mmff94s_parameters() {
        let mol = torsion_chain_molecule(
            [
                (Element::H, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::F, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[5, 1, 1, 10]);
        props.variant = MmffVariant::Mmff94s;
        let positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ];
        let mut ff = force_field_with_positions(&positions);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("MMFF94s torsion row should add a contrib");

        assert_eq!(ff.contribs().len(), 1);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2.0, 1.0, 1.0];
        let params = props
            .get_mmff_torsion_params(0, 1, 2, 3)
            .expect("MMFF94s torsion lookup should parse")
            .expect("MMFF94s torsion params should be tabulated");
        let expected = source_torsion_energy(
            params.1.v1,
            params.1.v2,
            params.1.v3,
            source_torsion_cos_phi(
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(1.0, 1.0, 0.0),
                ForceFieldVec3::new(2.0, 1.0, 1.0),
            ),
        );
        assert_close(ff.contribs()[0].get_energy(&pos), expected);
    }

    #[test]
    #[should_panic(expected = "missing atom types - invalid force-field")]
    fn mmff_builder_add_torsions_rejects_invalid_mmff_properties() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let mut props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 1]);
        props.valid = false;
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);

        let _ = add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]");
    }

    #[test]
    fn mmff_builder_add_torsions_accepts_ported_empirical_fallback() {
        let mol = torsion_chain_molecule(
            [
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::C, Hybridization::Sp3),
                (Element::MG, Hybridization::Unspecified),
            ],
            [BondOrder::Single, BondOrder::Single, BondOrder::Single],
        );
        let props = mmff_props_for_molecule_and_atom_types(mol.clone(), &[1, 1, 1, 99]);
        let mut ff = force_field_with_positions(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 1.0],
        ]);

        add_torsions(&mol, &props, &mut ff, "[!$(*#*)&!D1]~[!$(*#*)&!D1]")
            .expect("ported empirical torsion fallback should be supported");

        assert_eq!(ff.contribs().len(), 1);
    }
}
