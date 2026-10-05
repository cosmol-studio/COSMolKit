//! RDKit MMFF properties over detached topology values.
//! Original modern Rust source: 114677c74afca511a02fdbf9fcd9de394fb5ad75c15c35b4c2e325a8fbcb3d56.
//! Upstream pin: 351f8f378f8ad6bbd517980c38896e66bf907af8.

use super::nonbonded::{calc_unscaled_vdw_minimum, calc_unscaled_vdw_well_depth, scale_vdw_params};
use super::params::*;
use super::params_text::MmffParamParseError;
use cosmolkit_core::{
    AromaticityError, AromaticityModel, AromaticityParams, KekulizeError, KekulizeParams,
    RingFindingError, RingInfo, RingSearchParams, ValenceError, ValenceModel,
    assign_mmff_aromaticity_prepared, assign_valence_for_topology, kekulize, symmetrized_sssr,
};
use cosmolkit_model::{Atom, AtomId, Bond, BondOrder, TopologyBlock};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MmffParamError {
    #[error(transparent)]
    Table(#[from] MmffParamParseError),
    #[error(transparent)]
    Vdw(#[from] MmffVdwParseError),
    #[error(transparent)]
    Angle(#[from] MmffAngleLookupError),
    #[error(transparent)]
    Torsion(#[from] MmffTorLookupError),
    #[error(transparent)]
    Oop(#[from] MmffOopLookupError),
    #[error("missing MMFF default data: {symbol}")]
    MissingDefaultData { symbol: &'static str },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MmffVdwRijstarEps {
    pub r_ij_star_unscaled: f64,
    pub epsilon_unscaled: f64,
    pub r_ij_star: f64,
    pub epsilon: f64,
}

pub const MMFF_SANITIZED_PROP: &str = "_MMFFSanitized";
pub const MMFF_DIELECTRIC_CONSTANT: u8 = 1;
pub const MMFF_DIELECTRIC_DISTANCE: u8 = 2;
pub const MMFF_VERBOSITY_NONE: u8 = 0;
pub const MMFF_VERBOSITY_LOW: u8 = 1;
pub const MMFF_VERBOSITY_HIGH: u8 = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MmffVariant {
    Mmff94,
    Mmff94s,
}

impl MmffVariant {
    #[must_use]
    pub fn from_rdkit_constructor_arg(mmff_variant: &str) -> Self {
        if mmff_variant == "MMFF94s" {
            Self::Mmff94s
        } else {
            Self::Mmff94
        }
    }

    #[must_use]
    pub const fn as_rdkit_str(self) -> &'static str {
        match self {
            Self::Mmff94 => "MMFF94",
            Self::Mmff94s => "MMFF94s",
        }
    }

    #[must_use]
    pub const fn is_mmff94s(self) -> bool {
        matches!(self, Self::Mmff94s)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MmffAtomProperties {
    pub atom_type: u8,
    pub formal_charge: f64,
    pub partial_charge: f64,
}

impl Default for MmffAtomProperties {
    fn default() -> Self {
        Self {
            // RDKit✔️✔️: std::uint8_t mmffAtomType{0};
            atom_type: 0,
            // RDKit✔️✔️: double mmffFormalCharge{0.0};
            formal_charge: 0.0,
            // RDKit✔️✔️: double mmffPartialCharge{0.0};
            partial_charge: 0.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MmffMolProperties {
    pub topology: TopologyBlock,
    pub valid: bool,
    pub variant: MmffVariant,
    pub bond_term: bool,
    pub angle_term: bool,
    pub stretch_bend_term: bool,
    pub oop_term: bool,
    pub torsion_term: bool,
    pub vdw_term: bool,
    pub ele_term: bool,
    pub dielectric_constant: f64,
    pub dielectric_model: u8,
    pub verbosity: u8,
    pub atom_properties: Vec<MmffAtomProperties>,
    pub aromatic_ring_count: usize,
    pub acquired_rings: Option<RingInfo>,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MmffMolPropertiesError {
    #[error(transparent)]
    Params(#[from] MmffParamError),
    #[error(transparent)]
    Kekulize(#[from] KekulizeError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error("MMFF atom property index {atom_index} is out of range for {atoms} atoms")]
    AtomIndexOutOfRange { atom_index: usize, atoms: usize },
    #[error("MMFF atom type {atom_type} has no property row")]
    AtomTypePropertiesMissing { atom_type: u8 },
    #[error("MMFF atom type {atom_type} has no PBCI row")]
    AtomTypePbciMissing { atom_type: u8 },
}

fn is_double_zero(value: f64) -> bool {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::isDoubleZero (Params.h:41-43)
    // RDKit✔️✔️: inline bool isDoubleZero(const double x) {
    // RDKit✔️✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit✔️✔️: }
    value < 1.0e-10 && value > -1.0e-10
}

fn is_mmff_atom_n_oxide(
    mol: &TopologyBlock,
    implicit_hydrogens: &[i32],
    atom_index: usize,
) -> Result<bool, MmffMolPropertiesError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::isAtomNOxide (AtomTyper.cpp:335-352)
    // RDKit✔️✔️: // returns true if the atom is an N-oxide
    // RDKit✔️✔️: bool isAtomNOxide(const Atom *atom) {
    // RDKit✔️✔️:   bool isNOxide = false;
    let mut is_n_oxide = false;
    // RDKit✔️✔️:   const ROMol &mol = atom->getOwningMol();
    let atom = &mol.atoms.as_slice()[atom_index];
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx;
    // RDKit✔️✔️:   ROMol::ADJ_ITER endNbrs;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((atom->getAtomicNum() == 7) && (atom->getTotalDegree() >= 3)) {
    if atom.atomic_number() == 7 && mmff_total_degree(mol, implicit_hydrogens, atom_index)? >= 3 {
        // RDKit✔️✔️:     // loop over neighbors
        // RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
        // RDKit✔️✔️:     for (; (!isNOxide) && (nbrIdx != endNbrs); ++nbrIdx) {
        for neighbor in mol.adjacency.neighbors_of(atom_index) {
            if is_n_oxide {
                break;
            }
            // RDKit✔️✔️:       const Atom *nbrAtom = mol[*nbrIdx];
            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
            // RDKit✔️✔️:       isNOxide =
            // RDKit✔️✔️:           ((nbrAtom->getAtomicNum() == 8) && (nbrAtom->getTotalDegree() == 1));
            is_n_oxide = nbr_atom.atomic_number() == 8
                && mmff_total_degree(mol, implicit_hydrogens, neighbor.atom_index)? == 1;
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return isNOxide;
    Ok(is_n_oxide)
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MMFF::isAtomNOxide
}

impl MmffMolProperties {
    pub fn new(
        mol: &TopologyBlock,
        mmff_sanitized: bool,
        mmff_variant: &str,
        verbosity: u8,
    ) -> Result<Self, MmffMolPropertiesError> {
        Self::new_prepared(mol, mmff_sanitized, mmff_variant, verbosity, None)
    }

    pub fn new_prepared(
        mol: &TopologyBlock,
        mmff_sanitized: bool,
        mmff_variant: &str,
        verbosity: u8,
        supplied_rings: Option<&RingInfo>,
    ) -> Result<Self, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::MMFFMolProperties::MMFFMolProperties (AtomTyper.cpp:2328-2408)
        // RDKit❗✔️: // constructs a MMFFMolProperties object for ROMol mol filled
        // RDKit❗✔️: // with MMFF atom types, formal and partial charges
        // RDKit❗✔️: // in case atom types are missing, d_valid is set to false,
        // RDKit❗✔️: // charges are set to 0.0 and the force-field is unusable
        // RDKit❗✔️: MMFFMolProperties::MMFFMolProperties(ROMol &mol, const std::string &mmffVariant,
        // RDKit❗✔️:                                      std::uint8_t verbosity,
        // RDKit❗✔️:                                      std::ostream &oStream)
        // RDKit❗✔️:     : d_valid(true),
        let mut valid = true;
        // RDKit✔️✔️:       d_mmffs(mmffVariant == "MMFF94s"),
        let variant = MmffVariant::from_rdkit_constructor_arg(mmff_variant);
        // RDKit✔️✔️:       d_bondTerm(true),
        let bond_term = true;
        // RDKit✔️✔️:       d_angleTerm(true),
        let angle_term = true;
        // RDKit✔️✔️:       d_stretchBendTerm(true),
        let stretch_bend_term = true;
        // RDKit✔️✔️:       d_oopTerm(true),
        let oop_term = true;
        // RDKit✔️✔️:       d_torsionTerm(true),
        let torsion_term = true;
        // RDKit✔️✔️:       d_vdWTerm(true),
        let vdw_term = true;
        // RDKit✔️✔️:       d_eleTerm(true),
        let ele_term = true;
        // RDKit✔️✔️:       d_dielConst(1.0),
        let dielectric_constant = 1.0;
        // RDKit✔️✔️:       d_dielModel(CONSTANT),
        let dielectric_model = MMFF_DIELECTRIC_CONSTANT;
        // RDKit✔️✔️:       d_verbosity(verbosity),
        // RDKit✔️✔️:       d_oStream(&oStream),
        // RDKit✔️✔️:       d_MMFFAtomPropertiesPtrVect(mol.getNumAtoms()) {
        let mut atom_properties = vec![MmffAtomProperties::default(); mol.atoms.len()];
        // RDKit❗✔️:   if (MolOps::needsHs(mol)) {
        // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
        // RDKit❗✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
        // RDKit❗✔️:         << std::endl;
        // RDKit❗✔️:   }
        // COSMolKit does not currently model RDKit warning-log side effects.
        // RDKit✔️✔️:   if (!mol.hasProp(common_properties::_MMFFSanitized)) {
        // RDKit✔️✔️:     bool isAromaticSet = false;
        // RDKit✔️✔️:     for (const auto atom : mol.atoms()) {
        // RDKit✔️✔️:       if (atom->getIsAromatic()) {
        // RDKit✔️✔️:         isAromaticSet = true;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut prepared_molecule = mol.clone();
        if !mmff_sanitized {
            let is_aromatic_set = prepared_molecule
                .atoms
                .as_slice()
                .iter()
                .any(Atom::is_aromatic);
            // RDKit✔️✔️:     if (isAromaticSet) {
            // RDKit✔️✔️:       MolOps::Kekulize((RWMol &)mol, true);
            // RDKit✔️✔️:     }
            if is_aromatic_set {
                prepared_molecule = kekulize(
                    &prepared_molecule,
                    &KekulizeParams {
                        mark_atoms_bonds: true,
                        ..KekulizeParams::default()
                    },
                )?
                .topology;
            }
            // RDKit✔️✔️:     mol.setProp(common_properties::_MMFFSanitized, 1, true);
            // RDKit✔️✔️:   }
            // The caller owns _MMFFSanitized property installation/lifecycle.
        }
        // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        // RDKit✔️✔️:     d_MMFFAtomPropertiesPtrVect[i] =
        // RDKit✔️✔️:         MMFFAtomPropertiesPtr(new MMFFAtomProperties());
        // RDKit✔️✔️:   }
        // RDKit❗✔️:   MolOps::setMMFFAromaticity((RWMol &)mol);
        // Behavior: the existing core MMFF owner includes the noncarbon
        // explicit-H adjustment using pre-aromaticity valence. Move its final
        // detached topology; never duplicate aromaticity or implicit-H logic.
        // Complexity: retain the existing core traversal and one detached
        // topology; no second prepared topology/aromaticity snapshot is kept.
        // RDKit❗✔️: RingMembershipSize::RingMembershipSize(const ROMol &mol) {
        // RDKit❗✔️:   static const unsigned int MAX_NUM_RINGS = (0xFFFFFFFF >> 1);
        // RDKit❗✔️:   const RingInfo *ringInfo = mol.getRingInfo();
        // RDKit❗✔️:   const VECT_INT_VECT &atomRings = ringInfo->atomRings();
        // RDKit❗✔️:   PRECONDITION(atomRings.size() < MAX_NUM_RINGS, "Too many rings");
        // RDKit❗✔️:   for (std::uint32_t ringIdx = 0; ringIdx < atomRings.size(); ++ringIdx) {
        // RDKit❗✔️:     unsigned int ringSize = atomRings[ringIdx].size();
        // RDKit❗✔️:     std::uint32_t ringIdxWithAromaticFlag = ringIdx;
        // RDKit❗✔️:     bool ringIsAromatic = isRingAromatic(mol, atomRings[ringIdx]);
        // RDKit❗✔️:     if (ringIsAromatic) {
        // RDKit❗✔️:       ringIdxWithAromaticFlag |= IS_AROMATIC_BIT;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     auto it = d_ringSizeMembershipMap.find(ringSize);
        // RDKit❗✔️:     if (it == d_ringSizeMembershipMap.end()) {
        // RDKit❗✔️:       it = d_ringSizeMembershipMap
        // RDKit❗✔️:                .insert(std::make_pair(ringSize, RingMembershipMap()))
        // RDKit❗✔️:                .first;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     for (int atomIdxIt : atomRings[ringIdx]) {
        // RDKit❗✔️:       auto it2 = it->second.find(atomIdxIt);
        // RDKit❗✔️:       if (it2 == it->second.end()) {
        // RDKit❗✔️:         it2 = it->second.insert(std::make_pair(atomIdxIt, RingMembership()))
        // RDKit❗✔️:                   .first;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       it2->second.getRingIdxSet().insert(ringIdxWithAromaticFlag);
        // RDKit❗✔️:       if (ringIsAromatic) {
        // RDKit❗✔️:         it2->second.setIsInAromaticRing(true);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Source setMMFFAromaticity and RingMembershipSize both read the same
        // stored ring rows. Preserve supplied initialized rows, including empty
        // rows and Fast state. The detached boundary acquires missing rows once.
        // Source initialized RingInfo membership queries return empty beyond
        // the stored dimensions. AddHs preserves those existing ring rows.
        // Transport that state to the core exact-dimension carrier by growing
        // only empty membership rows, without finding or replacing rings.
        // RDKit❗❌: void RingInfo::preallocate(unsigned int numAtoms, unsigned int numBonds) {
        // RDKit❗❌:   d_atomMembers.resize(numAtoms);
        // RDKit❗❌:   d_bondMembers.resize(numBonds);
        // RDKit❗❌: }
        // Behavior: ring/family rows and find quality are unchanged.
        // Complexity: undersized carriers are cloned at the detached transport;
        // the source grows in place. Full-size carriers remain borrowed.
        let acquired_rings = match supplied_rings {
            None => Some(symmetrized_sssr(
                &prepared_molecule,
                &RingSearchParams::default(),
            )?),
            Some(rings)
                if rings.is_initialized()
                    && rings.atom_row_count() <= prepared_molecule.atoms.len()
                    && rings.bond_row_count() <= prepared_molecule.bonds.len()
                    && (rings.atom_row_count() < prepared_molecule.atoms.len()
                        || rings.bond_row_count() < prepared_molecule.bonds.len()) =>
            {
                let mut transported = rings.clone();
                transported
                    .preallocate(prepared_molecule.atoms.len(), prepared_molecule.bonds.len());
                Some(transported)
            }
            Some(_) => None,
        };
        let ring_info = acquired_rings
            .as_ref()
            .or(supplied_rings)
            .expect("supplied or acquired rings");
        let aromaticity = assign_mmff_aromaticity_prepared(&prepared_molecule, ring_info)?;
        let aromatic_ring_count = aromaticity.aromatic_ring_count;
        prepared_molecule = aromaticity.topology;
        // RDKit❗✔️:   RingMembershipSize rmSize(mol);

        let valence = assign_valence_for_topology(&prepared_molecule, ValenceModel::RdkitLike)?;
        // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
        // RDKit❗✔️:     if (atom->getAtomicNum() != 1) {
        // RDKit❗✔️:       this->setMMFFHeavyAtomType(rmSize, atom);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        for atom in prepared_molecule.atoms.as_slice() {
            if atom.atomic_number() != 1 {
                atom_properties[atom.id().index()] = set_mmff_heavy_atom_type(
                    &prepared_molecule,
                    &ring_info,
                    &valence.explicit_valence,
                    &valence.implicit_hydrogens,
                    atom,
                )?;
                if atom_properties[atom.id().index()].atom_type == 0 {
                    valid = false;
                }
            }
        }
        // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
        // RDKit✔️✔️:     if (atom->getAtomicNum() == 1) {
        // RDKit✔️✔️:       this->setMMFFHydrogenType(atom);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        for atom in prepared_molecule.atoms.as_slice() {
            if atom.atomic_number() == 1 {
                let properties =
                    set_mmff_hydrogen_type(&prepared_molecule, &atom_properties, atom)?;
                if properties.atom_type == 0 {
                    valid = false;
                }
                atom_properties[atom.id().index()] = properties;
            }
        }
        // RDKit✔️✔️:   if (this->isValid()) {
        // RDKit✔️✔️:     this->computeMMFFCharges(mol);
        // RDKit✔️✔️:   }
        // RDKit❌❌:   if (verbosity == MMFF_VERBOSITY_HIGH) {
        // RDKit❌❌:     oStream << "\n"
        // RDKit❌❌:                "A T O M   T Y P E S   A N D   C H A R G E S\n\n"
        // RDKit❌❌:                "          ATOM    FORMAL   PARTIAL\n"
        // RDKit❌❌:                " ATOM     TYPE    CHARGE    CHARGE\n"
        // RDKit❌❌:                "-----------------------------------"
        // RDKit❌❌:             << std::endl;
        // RDKit❌❌:     for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
        // RDKit❌❌:       oStream << std::left << std::setw(2)
        // RDKit❌❌:               << mol.getAtomWithIdx(idx)->getSymbol() << std::left << " #"
        // RDKit❌❌:               << std::setw(5) << idx + 1 << std::right << std::setw(5)
        // RDKit❌❌:               << (unsigned int)(this->getMMFFAtomType(idx)) << std::right
        // RDKit❌❌:               << std::setw(10) << std::fixed << std::setprecision(3)
        // RDKit❌❌:               << this->getMMFFFormalCharge(idx) << std::right << std::setw(10)
        // RDKit❌❌:               << this->getMMFFPartialCharge(idx) << std::endl;
        // RDKit❌❌:     }
        // RDKit❌❌:     if (!(this->isValid())) {
        // RDKit❌❌:       oStream << "\nMissing atom types - charges were not computed"
        // RDKit❌❌:               << std::endl;
        // RDKit❌❌:     }
        // RDKit❌❌:   }
        // RDKit❗✔️: }
        let mut properties = Self {
            topology: prepared_molecule,
            valid,
            variant,
            bond_term,
            angle_term,
            stretch_bend_term,
            oop_term,
            torsion_term,
            vdw_term,
            ele_term,
            dielectric_constant,
            dielectric_model,
            verbosity,
            atom_properties,
            aromatic_ring_count,
            acquired_rings: None,
        };
        if properties.is_valid() {
            properties.compute_mmff_charges(&ring_info)?;
        }
        properties.acquired_rings = acquired_rings;
        Ok(properties)
    }

    fn compute_mmff_charges(&mut self, ring_info: &RingInfo) -> Result<(), MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::computeMMFFCharges (AtomTyper.cpp:3071-3487)
        // RDKit✔️✔️: void MMFFMolProperties::computeMMFFCharges(const ROMol &mol) {
        // RDKit✔️✔️:   PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
        assert!(self.is_valid(), "missing atom types - invalid force-field");

        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int idx;
        // RDKit✔️✔️:   unsigned int i;
        // RDKit✔️✔️:   unsigned int j;
        // RDKit✔️✔️:   unsigned int atomType;
        // RDKit✔️✔️:   unsigned int nbrAtomType;
        // RDKit✔️✔️:   unsigned int nConj = 0;
        // RDKit✔️✔️:   unsigned int old_nConj = 0;
        // RDKit✔️✔️:   double pChg = 0.0;
        // RDKit✔️✔️:   double fChg = 0.0;
        // RDKit✔️✔️:   boost::dynamic_bitset<> conjNBitVect(mol.getNumAtoms());
        // RDKit✔️✔️:   VECT_INT_VECT atomRings = mol.getRingInfo()->atomRings();
        // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx;
        // RDKit✔️✔️:   ROMol::ADJ_ITER endNbrs;
        // RDKit✔️✔️:   ROMol::ADJ_ITER nbr2Idx;
        // RDKit✔️✔️:   ROMol::ADJ_ITER end2Nbrs;
        // RDKit✔️✔️:   const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
        let mmff_prop = default_mmff_prop().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:   const MMFFPBCICollection *mmffPBCI = DefaultParameters::getMMFFPBCI();
        let mmff_pbci = default_mmff_pbci().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:   const MMFFChgCollection *mmffChg = DefaultParameters::getMMFFChg();
        let mmff_chg = default_mmff_chg().map_err(|e| MmffParamError::from(e.clone()))?;
        let mol = &self.topology;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   // We need to set formal charges upfront
        // RDKit✔️✔️:   for (idx = 0; idx < mol.getNumAtoms(); ++idx) {
        for idx in 0..mol.atoms.len() {
            // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(idx);
            let atom = &mol.atoms.as_slice()[idx];
            // RDKit✔️✔️:     atomType = this->getMMFFAtomType(idx);
            let atom_type = self.atom_properties[idx].atom_type;
            // RDKit✔️✔️:     fChg = 0.0;
            let mut formal_charge = 0.0;
            // RDKit✔️✔️:     switch (atomType) {
            match atom_type {
                // RDKit✔️✔️:       // special cases
                // RDKit✔️✔️:       case 32:
                // RDKit✔️✔️:       // O2CM
                // RDKit✔️✔️:       // Oxygen in carboxylate group
                // RDKit✔️✔️:       case 72:
                // RDKit✔️✔️:         // SM
                // RDKit✔️✔️:         // Anionic terminal sulfur
                32 | 72 => {
                    // RDKit✔️✔️:         // loop over neighbors
                    // RDKit✔️✔️:         boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit✔️✔️:         for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(idx) {
                        // RDKit✔️✔️:           const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        // RDKit✔️✔️:           nbrAtomType = this->getMMFFAtomType(nbrAtom->getIdx());
                        let nbr_atom_type = self.atom_properties[neighbor.atom_index].atom_type;
                        // RDKit✔️✔️:           // loop over neighbors of the neighbor
                        // RDKit✔️✔️:           // count how many terminal oxygen/sulfur atoms
                        // RDKit✔️✔️:           // or secondary nitrogens
                        // RDKit✔️✔️:           // are bonded to the neighbor of ipso
                        // RDKit✔️✔️:           int nSecNbondedToNbr = 0;
                        let mut n_sec_n_bonded_to_nbr = 0_i32;
                        // RDKit✔️✔️:           int nTermOSbondedToNbr = 0;
                        let mut n_term_os_bonded_to_nbr = 0_i32;
                        // RDKit✔️✔️:           boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                        // RDKit✔️✔️:           for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                        for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                            // RDKit✔️✔️:             const Atom *nbr2Atom = mol[*nbr2Idx];
                            let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                            // RDKit✔️✔️:             // if it's nitrogen with 2 neighbors and it is not aromatic,
                            // RDKit✔️✔️:             // increment the counter of secondary nitrogens
                            // RDKit✔️✔️:             if ((nbr2Atom->getAtomicNum() == 7) &&
                            // RDKit✔️✔️:                 (nbr2Atom->getDegree() == 2) &&
                            // RDKit✔️✔️:                 (!(nbr2Atom->getIsAromatic()))) {
                            if nbr2_atom.atomic_number() == 7
                                && mol.adjacency.neighbors_of(neighbor2.atom_index).len() == 2
                                && !nbr2_atom.is_aromatic()
                            {
                                // RDKit✔️✔️:               ++nSecNbondedToNbr;
                                n_sec_n_bonded_to_nbr += 1;
                                // RDKit✔️✔️:             }
                            }
                            // RDKit✔️✔️:             // if it's terminal oxygen/sulfur,
                            // RDKit✔️✔️:             // increment the terminal oxygen/sulfur counter
                            // RDKit✔️✔️:             if (((nbr2Atom->getAtomicNum() == 8) ||
                            // RDKit✔️✔️:                  (nbr2Atom->getAtomicNum() == 16)) &&
                            // RDKit✔️✔️:                 (nbr2Atom->getDegree() == 1)) {
                            if matches!(nbr2_atom.atomic_number(), 8 | 16)
                                && mol.adjacency.neighbors_of(neighbor2.atom_index).len() == 1
                            {
                                // RDKit✔️✔️:               ++nTermOSbondedToNbr;
                                n_term_os_bonded_to_nbr += 1;
                                // RDKit✔️✔️:             }
                            }
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           // in case its sulfur with two terminal oxygen/sulfur atoms and one
                        // RDKit✔️✔️:           // secondary
                        // RDKit✔️✔️:           // nitrogen, this is a deprotonated sulfonamide, so we should not
                        // RDKit✔️✔️:           // consider
                        // RDKit✔️✔️:           // nitrogen as a replacement for oxygen/sulfur in a sulfone
                        // RDKit✔️✔️:           if ((nbrAtom->getAtomicNum() == 16) && (nTermOSbondedToNbr == 2) &&
                        // RDKit✔️✔️:               (nSecNbondedToNbr == 1)) {
                        if nbr_atom.atomic_number() == 16
                            && n_term_os_bonded_to_nbr == 2
                            && n_sec_n_bonded_to_nbr == 1
                        {
                            // RDKit✔️✔️:             nSecNbondedToNbr = 0;
                            n_sec_n_bonded_to_nbr = 0;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           // if the neighbor is carbon
                        // RDKit✔️✔️:           if ((nbrAtom->getAtomicNum() == 6) && nTermOSbondedToNbr) {
                        if nbr_atom.atomic_number() == 6 && n_term_os_bonded_to_nbr != 0 {
                            // RDKit✔️✔️:             // O2CM
                            // RDKit✔️✔️:             // Oxygen in (thio)carboxylate group: charge is shared
                            // RDKit✔️✔️:             // across 2 oxygens/sulfur atoms in (thio)carboxylate,
                            // RDKit✔️✔️:             // 3 oxygen/sulfur atoms in (thio)carbonate
                            // RDKit✔️✔️:             // SM
                            // RDKit✔️✔️:             // Anionic terminal sulfur: charge is localized
                            // RDKit✔️✔️:             fChg = ((nTermOSbondedToNbr == 1)
                            // RDKit✔️✔️:                         ? -1.0
                            // RDKit✔️✔️:                         : -((double)(nTermOSbondedToNbr - 1) /
                            // RDKit✔️✔️:                             (double)nTermOSbondedToNbr));
                            formal_charge = if n_term_os_bonded_to_nbr == 1 {
                                -1.0
                            } else {
                                -f64::from(n_term_os_bonded_to_nbr - 1)
                                    / f64::from(n_term_os_bonded_to_nbr)
                            };
                            // RDKit✔️✔️:             break;
                            break;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           // if the neighbor is NO2 or NO3
                        // RDKit✔️✔️:           if ((nbrAtomType == 45) && (nTermOSbondedToNbr == 3)) {
                        if nbr_atom_type == 45 && n_term_os_bonded_to_nbr == 3 {
                            // RDKit✔️✔️:             // O3N
                            // RDKit✔️✔️:             // Nitrate anion oxygen
                            // RDKit✔️✔️:             fChg = -1.0 / 3.0;
                            formal_charge = -1.0 / 3.0;
                            // RDKit✔️✔️:             break;
                            break;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           // if the neighbor is PO2, PO3, PO4
                        // RDKit✔️✔️:           if ((nbrAtomType == 25) && nTermOSbondedToNbr) {
                        if nbr_atom_type == 25 && n_term_os_bonded_to_nbr != 0 {
                            // RDKit✔️✔️:             // OP
                            // RDKit✔️✔️:             // Oxygen in phosphine oxide
                            // RDKit✔️✔️:             // O2P
                            // RDKit✔️✔️:             // One of 2 terminal O's on P
                            // RDKit✔️✔️:             // O3P
                            // RDKit✔️✔️:             // One of 3 terminal O's on P
                            // RDKit✔️✔️:             // O4P
                            // RDKit✔️✔️:             // One of 4 terminal O's on P
                            // RDKit✔️✔️:             fChg = ((nTermOSbondedToNbr == 1)
                            // RDKit✔️✔️:                         ? 0.0
                            // RDKit✔️✔️:                         : -((double)(nTermOSbondedToNbr - 1) /
                            // RDKit✔️✔️:                             (double)nTermOSbondedToNbr));
                            formal_charge = if n_term_os_bonded_to_nbr == 1 {
                                0.0
                            } else {
                                -f64::from(n_term_os_bonded_to_nbr - 1)
                                    / f64::from(n_term_os_bonded_to_nbr)
                            };
                            // RDKit✔️✔️:             break;
                            break;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           // if the neighbor is SO2, SO2N, SO3, SO4, SO2M, SSOM
                        // RDKit✔️✔️:           if ((nbrAtomType == 18) && nTermOSbondedToNbr) {
                        if nbr_atom_type == 18 && n_term_os_bonded_to_nbr != 0 {
                            // RDKit✔️✔️:             // SO2
                            // RDKit✔️✔️:             // Sulfone sulfur
                            // RDKit✔️✔️:             // SO2N
                            // RDKit✔️✔️:             // Sulfonamide sulfur
                            // RDKit✔️✔️:             // SO3
                            // RDKit✔️✔️:             // Sulfonate group sulfur
                            // RDKit✔️✔️:             // SO4
                            // RDKit✔️✔️:             // Sulfate group sulfur
                            // RDKit✔️✔️:             // SNO
                            // RDKit✔️✔️:             // Sulfur in nitrogen analog of a sulfone
                            // RDKit✔️✔️:             fChg =
                            // RDKit✔️✔️:                 (((nSecNbondedToNbr + nTermOSbondedToNbr) == 2)
                            // RDKit✔️✔️:                      ? 0.0
                            // RDKit✔️✔️:                      : -((double)((nSecNbondedToNbr + nTermOSbondedToNbr) - 2) /
                            // RDKit✔️✔️:                          (double)nTermOSbondedToNbr));
                            formal_charge = if n_sec_n_bonded_to_nbr + n_term_os_bonded_to_nbr == 2
                            {
                                0.0
                            } else {
                                -f64::from(n_sec_n_bonded_to_nbr + n_term_os_bonded_to_nbr - 2)
                                    / f64::from(n_term_os_bonded_to_nbr)
                            };
                            // RDKit✔️✔️:             break;
                            break;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           if ((nbrAtomType == 73) && nTermOSbondedToNbr) {
                        if nbr_atom_type == 73 && n_term_os_bonded_to_nbr != 0 {
                            // RDKit✔️✔️:             // SO2M
                            // RDKit✔️✔️:             // Sulfur in anionic sulfinate group
                            // RDKit✔️✔️:             // SSOM
                            // RDKit✔️✔️:             // Tricoordinate sulfur in anionic thiosulfinate group
                            // RDKit✔️✔️:             fChg = ((nTermOSbondedToNbr == 1)
                            // RDKit✔️✔️:                         ? 0.0
                            // RDKit✔️✔️:                         : -((double)(nTermOSbondedToNbr - 1) /
                            // RDKit✔️✔️:                             (double)nTermOSbondedToNbr));
                            formal_charge = if n_term_os_bonded_to_nbr == 1 {
                                0.0
                            } else {
                                -f64::from(n_term_os_bonded_to_nbr - 1)
                                    / f64::from(n_term_os_bonded_to_nbr)
                            };
                            // RDKit✔️✔️:             break;
                            break;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:           if ((nbrAtomType == 77) && nTermOSbondedToNbr) {
                        if nbr_atom_type == 77 && n_term_os_bonded_to_nbr != 0 {
                            // RDKit✔️✔️:             // O4Cl
                            // RDKit✔️✔️:             // Oxygen in perchlorate anion
                            // RDKit✔️✔️:             fChg = -(1.0 / (double)nTermOSbondedToNbr);
                            formal_charge = -1.0 / f64::from(n_term_os_bonded_to_nbr);
                            // RDKit✔️✔️:             break;
                            break;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       case 76:
                76 => {
                    // RDKit✔️✔️:         // N5M
                    // RDKit✔️✔️:         // Nitrogen in 5-ring aromatic anion
                    // RDKit✔️✔️:         // we don't need to bother about the neighbors with N5M
                    // RDKit✔️✔️:         for (i = 0; i < atomRings.size(); ++i) {
                    // RDKit✔️✔️:           if ((std::find(atomRings[i].begin(), atomRings[i].end(), idx) !=
                    // RDKit✔️✔️:                atomRings[i].end())) {
                    // RDKit✔️✔️:             break;
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:         }
                    let atom_ring = ring_info
                        .atom_rings()
                        .iter()
                        .find(|ring| ring.contains(&AtomId::new(idx)));
                    // RDKit✔️✔️:         // find how many nitrogens with atom type 76 we have
                    // RDKit✔️✔️:         // and share the formal charge accordingly
                    // RDKit✔️✔️:         if (i < atomRings.size()) {
                    if let Some(atom_ring) = atom_ring {
                        // RDKit✔️✔️:           unsigned int nNitrogensIn5Ring = 0;
                        // RDKit✔️✔️:           for (j = 0; j < atomRings[i].size(); ++j) {
                        // RDKit✔️✔️:             if (this->getMMFFAtomType(atomRings[i][j]) == 76) {
                        // RDKit✔️✔️:               ++nNitrogensIn5Ring;
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:           }
                        let n_nitrogens_in_5_ring = atom_ring
                            .iter()
                            .filter(|atom_id| self.atom_properties[atom_id.index()].atom_type == 76)
                            .count();
                        // RDKit✔️✔️:           if (nNitrogensIn5Ring) {
                        if n_nitrogens_in_5_ring != 0 {
                            // RDKit✔️✔️:             fChg = -(1.0 / (double)nNitrogensIn5Ring);
                            formal_charge = -(1.0 / n_nitrogens_in_5_ring as f64);
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       case 55:
                // RDKit✔️✔️:       case 56:
                // RDKit✔️✔️:       case 81:
                55 | 56 | 81 => {
                    // RDKit✔️✔️:         // NIM+
                    // RDKit✔️✔️:         // Aromatic nitrogen in imidazolium
                    // RDKit✔️✔️:         // N5A+
                    // RDKit✔️✔️:         // Positive nitrogen in 5-ring alpha position
                    // RDKit✔️✔️:         // N5B+
                    // RDKit✔️✔️:         // Positive nitrogen in 5-ring beta position
                    // RDKit✔️✔️:         // N5+
                    // RDKit✔️✔️:         // Positive nitrogen in other 5-ring position
                    // RDKit✔️✔️:         // we need to loop over all molecule atoms
                    // RDKit✔️✔️:         // and find all those nitrogens with atom type
                    // RDKit✔️✔️:         // 81, 55 or 56, check whether they are conjugated
                    // RDKit✔️✔️:         // with ipso and keep on looping until no more
                    // RDKit✔️✔️:         // conjugated atoms can be found. Finally, we divide
                    // RDKit✔️✔️:         // the total formal charge that was found on the
                    // RDKit✔️✔️:         // conjugated system by the number of conjugated nitrogens
                    // RDKit✔️✔️:         // of types 81, 55 or 56 that were found.
                    // RDKit✔️✔️:         // This is not strictly what is described
                    // RDKit✔️✔️:         // in the MMFF papers, but it is the only way to get an
                    // RDKit✔️✔️:         // integer total formal charge, which makes sense to me
                    // RDKit✔️✔️:         // probably such conjugated systems are anyway out of the
                    // RDKit✔️✔️:         // scope of MMFF, but this is an attempt to correctly
                    // RDKit✔️✔️:         // deal with them somehow
                    // RDKit✔️✔️:         fChg = (double)(atom->getFormalCharge());
                    formal_charge = f64::from(atom.formal_charge());
                    // RDKit✔️✔️:         nConj = 1;
                    let mut n_conj = 1_usize;
                    // RDKit✔️✔️:         old_nConj = 0;
                    let mut old_n_conj = 0_usize;
                    // RDKit✔️✔️:         conjNBitVect.reset();
                    let mut conj_n_bit_vec = vec![false; mol.atoms.len()];
                    // RDKit✔️✔️:         conjNBitVect[idx] = 1;
                    conj_n_bit_vec[idx] = true;
                    // RDKit✔️✔️:         while (nConj > old_nConj) {
                    while n_conj > old_n_conj {
                        // RDKit✔️✔️:           old_nConj = nConj;
                        old_n_conj = n_conj;
                        // RDKit✔️✔️:           for (i = 0; i < mol.getNumAtoms(); ++i) {
                        for i in 0..mol.atoms.len() {
                            // RDKit✔️✔️:             // if this atom is not marked as conj, move on
                            // RDKit✔️✔️:             if (!conjNBitVect[i]) {
                            if !conj_n_bit_vec[i] {
                                // RDKit✔️✔️:               continue;
                                continue;
                                // RDKit✔️✔️:             }
                            }
                            // RDKit✔️✔️:             // loop over neighbors
                            // RDKit✔️✔️:             boost::tie(nbrIdx, endNbrs) =
                            // RDKit✔️✔️:                 mol.getAtomNeighbors(mol.getAtomWithIdx(i));
                            // RDKit✔️✔️:             for (; nbrIdx != endNbrs; ++nbrIdx) {
                            for neighbor in mol.adjacency.neighbors_of(i) {
                                // RDKit✔️✔️:               const Atom *nbrAtom = mol[*nbrIdx];
                                // RDKit✔️✔️:               nbrAtomType = this->getMMFFAtomType(nbrAtom->getIdx());
                                let nbr_atom_type =
                                    self.atom_properties[neighbor.atom_index].atom_type;
                                // RDKit✔️✔️:               // if atom type is not 80 or 57, move on
                                // RDKit✔️✔️:               if ((nbrAtomType != 57) && (nbrAtomType != 80)) {
                                if nbr_atom_type != 57 && nbr_atom_type != 80 {
                                    // RDKit✔️✔️:                 continue;
                                    continue;
                                    // RDKit✔️✔️:               }
                                }
                                // RDKit✔️✔️:               // loop over neighbors of the neighbor
                                // RDKit✔️✔️:               // if they are nitrogens of type 81, 55 or 56 and
                                // RDKit✔️✔️:               // they are not not marked as conjugated yet, do it
                                // RDKit✔️✔️:               // and increment the nConj counter by 1
                                // RDKit✔️✔️:               boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                                // RDKit✔️✔️:               for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                                for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                                    // RDKit✔️✔️:                 const Atom *nbr2Atom = mol[*nbr2Idx];
                                    let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                                    // RDKit✔️✔️:                 // if atom type is not 81, 55 or 56, move on
                                    // RDKit✔️✔️:                 nbrAtomType = this->getMMFFAtomType(nbr2Atom->getIdx());
                                    let nbr2_atom_type =
                                        self.atom_properties[neighbor2.atom_index].atom_type;
                                    // RDKit✔️✔️:                 if ((nbrAtomType != 55) && (nbrAtomType != 56) &&
                                    // RDKit✔️✔️:                     (nbrAtomType != 81)) {
                                    if !matches!(nbr2_atom_type, 55 | 56 | 81) {
                                        // RDKit✔️✔️:                   continue;
                                        continue;
                                        // RDKit✔️✔️:                 }
                                    }
                                    // RDKit✔️✔️:                 j = nbr2Atom->getIdx();
                                    let j = neighbor2.atom_index;
                                    // RDKit✔️✔️:                 // if this nitrogen is not yet marked as conjugated,
                                    // RDKit✔️✔️:                 // mark it and increment the counter and eventually
                                    // RDKit✔️✔️:                 // adjust the total formal charge of the conjugated system
                                    // RDKit✔️✔️:                 if (!conjNBitVect[j]) {
                                    if !conj_n_bit_vec[j] {
                                        // RDKit✔️✔️:                   conjNBitVect[j] = 1;
                                        conj_n_bit_vec[j] = true;
                                        // RDKit✔️✔️:                   fChg += (double)(nbr2Atom->getFormalCharge());
                                        formal_charge += f64::from(nbr2_atom.formal_charge());
                                        // RDKit✔️✔️:                   ++nConj;
                                        n_conj += 1;
                                        // RDKit✔️✔️:                 }
                                    }
                                    // RDKit✔️✔️:               }
                                }
                                // RDKit✔️✔️:             }
                            }
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:         if (nConj) {
                    if n_conj != 0 {
                        // RDKit✔️✔️:           fChg /= (double)nConj;
                        formal_charge /= n_conj as f64;
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       case 61:
                61 => {
                    // RDKit✔️✔️:         // loop over neighbors
                    // RDKit✔️✔️:         boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit✔️✔️:         for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(idx) {
                        // RDKit✔️✔️:           const Atom *nbrAtom = mol[*nbrIdx];
                        // RDKit✔️✔️:           // if it is diazonium, set a +1 formal charge on
                        // RDKit✔️✔️:           // the secondary nitrogen
                        // RDKit✔️✔️:           if (this->getMMFFAtomType(nbrAtom->getIdx()) == 42) {
                        if self.atom_properties[neighbor.atom_index].atom_type == 42 {
                            // RDKit✔️✔️:             fChg = 1.0;
                            formal_charge = 1.0;
                            // RDKit✔️✔️:           }
                        }
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       // non-complicated +1 atom types
                // RDKit✔️✔️:       case 34:
                // RDKit✔️✔️:       // NR+
                // RDKit✔️✔️:       // Quaternary nitrogen
                // RDKit✔️✔️:       case 49:
                // RDKit✔️✔️:       // O+
                // RDKit✔️✔️:       // Oxonium oxygen
                // RDKit✔️✔️:       case 51:
                // RDKit✔️✔️:       // O=+
                // RDKit✔️✔️:       // Oxenium oxygen
                // RDKit✔️✔️:       case 54:
                // RDKit✔️✔️:       // N+=C
                // RDKit✔️✔️:       // Iminium nitrogen
                // RDKit✔️✔️:       // N+=N
                // RDKit✔️✔️:       // Positively charged nitrogen doubly bonded to N
                // RDKit✔️✔️:       case 58:
                // RDKit✔️✔️:       // NPD+
                // RDKit✔️✔️:       // Aromatic nitrogen in pyridinium
                // RDKit✔️✔️:       case 92:
                // RDKit✔️✔️:       // LI+
                // RDKit✔️✔️:       // Lithium cation
                // RDKit✔️✔️:       case 93:
                // RDKit✔️✔️:       // NA+
                // RDKit✔️✔️:       // Sodium cation
                // RDKit✔️✔️:       case 94:
                // RDKit✔️✔️:       // K+
                // RDKit✔️✔️:       // Potassium cation
                // RDKit✔️✔️:       case 97:
                // RDKit✔️✔️:         // CU+1
                // RDKit✔️✔️:         // Monopositive copper cation
                34 | 49 | 51 | 54 | 58 | 92 | 93 | 94 | 97 => {
                    // RDKit✔️✔️:         fChg = 1.0;
                    formal_charge = 1.0;
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       // non-complicated +2 atom types
                // RDKit✔️✔️:       case 87:
                // RDKit✔️✔️:       // FE+2
                // RDKit✔️✔️:       // Dipositive iron cation
                // RDKit✔️✔️:       case 95:
                // RDKit✔️✔️:       // ZN+2
                // RDKit✔️✔️:       // Dipositive zinc cation
                // RDKit✔️✔️:       case 96:
                // RDKit✔️✔️:       // CA+2
                // RDKit✔️✔️:       // Dipositive calcium cation
                // RDKit✔️✔️:       case 98:
                // RDKit✔️✔️:       // CU+2
                // RDKit✔️✔️:       // Dipositive copper cation
                // RDKit✔️✔️:       case 99:
                // RDKit✔️✔️:         // MG+2
                // RDKit✔️✔️:         // Dipositive magnesium cation
                87 | 95 | 96 | 98 | 99 => {
                    // RDKit✔️✔️:         fChg = 2.0;
                    formal_charge = 2.0;
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       // non-complicated +3 atom types
                // RDKit✔️✔️:       case 88:
                // RDKit✔️✔️:         // FE+3
                // RDKit✔️✔️:         // Tripositive iron cation
                88 => {
                    // RDKit✔️✔️:         fChg = 3.0;
                    formal_charge = 3.0;
                    // RDKit✔️✔️:         break;
                }

                // RDKit✔️✔️:       // non-complicated -1 atom types
                // RDKit✔️✔️:       case 35:
                // RDKit✔️✔️:       // OM
                // RDKit✔️✔️:       // Oxide oxygen on sp3 carbon
                // RDKit✔️✔️:       // OM2
                // RDKit✔️✔️:       // Oxide oxygen on sp2 carbon
                // RDKit✔️✔️:       // OM
                // RDKit✔️✔️:       // Oxide oxygen on sp3 nitrogen (not in original MMFF.I Table III)
                // RDKit✔️✔️:       // OM2
                // RDKit✔️✔️:       // Oxide oxygen on sp2 nitrogen (not in original MMFF.I Table III)
                // RDKit✔️✔️:       case 62:
                // RDKit✔️✔️:       // NM
                // RDKit✔️✔️:       // Anionic divalent nitrogen
                // RDKit✔️✔️:       case 89:
                // RDKit✔️✔️:       // F-
                // RDKit✔️✔️:       // Fluoride anion
                // RDKit✔️✔️:       case 90:
                // RDKit✔️✔️:       // Cl-
                // RDKit✔️✔️:       // Chloride anion
                // RDKit✔️✔️:       case 91:
                // RDKit✔️✔️:         // BR-
                // RDKit✔️✔️:         // Bromide anion
                35 | 62 | 89 | 90 | 91 => {
                    // RDKit✔️✔️:         fChg = -1.0;
                    formal_charge = -1.0;
                    // RDKit✔️✔️:         break;
                }
                _ => {} // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     this->setMMFFFormalCharge(idx, fChg);
            self.atom_properties[idx].formal_charge = formal_charge;
            // RDKit✔️✔️:   }
        }

        // RDKit✔️✔️:   // now we compute partial charges
        // RDKit✔️✔️:   // See Halgren, T. MMFF.V, J. Comput. Chem. 1996, 17, 616-641
        // RDKit✔️✔️:   // https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<616::AID-JCC5>3.0.CO;2-X
        // RDKit✔️✔️:   for (idx = 0; idx < mol.getNumAtoms(); ++idx) {
        for idx in 0..mol.atoms.len() {
            // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(idx);
            // RDKit✔️✔️:     atomType = this->getMMFFAtomType(idx);
            let atom_type = self.atom_properties[idx].atom_type;
            // RDKit✔️✔️:     double q0 = this->getMMFFFormalCharge(idx);
            let mut q0 = self.atom_properties[idx].formal_charge;
            // RDKit✔️✔️:     auto M = (double)((*mmffProp)(atomType)->crd);
            let atom_prop = mmff_prop
                .get(u32::from(atom_type))
                .ok_or(MmffMolPropertiesError::AtomTypePropertiesMissing { atom_type })?;
            let m = f64::from(atom_prop.crd);
            // RDKit✔️✔️:     double v = (*mmffPBCI)(atomType)->fcadj;
            let atom_pbci = mmff_pbci
                .get(u32::from(atom_type))
                .ok_or(MmffMolPropertiesError::AtomTypePbciMissing { atom_type })?;
            let v = atom_pbci.fcadj;
            // RDKit✔️✔️:     double sumFormalCharge = 0.0;
            let mut sum_formal_charge = 0.0;
            // RDKit✔️✔️:     double sumPartialCharge = 0.0;
            let mut sum_partial_charge = 0.0;
            // RDKit✔️✔️:     double nbrFormalCharge;
            // RDKit✔️✔️:     std::pair<int, const MMFFChg *> mmffChgParams;

            // RDKit✔️✔️:
            // RDKit✔️✔️:     if (isDoubleZero(v)) {
            if is_double_zero(v) {
                // RDKit✔️✔️:       // loop over neighbors
                // RDKit✔️✔️:       boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                // RDKit✔️✔️:       for (; nbrIdx != endNbrs; ++nbrIdx) {
                for neighbor in mol.adjacency.neighbors_of(idx) {
                    // RDKit✔️✔️:         const Atom *nbrAtom = mol[*nbrIdx];
                    // RDKit✔️✔️:         nbrFormalCharge = this->getMMFFFormalCharge(nbrAtom->getIdx());
                    let nbr_formal_charge = self.atom_properties[neighbor.atom_index].formal_charge;
                    // RDKit✔️✔️:         // if neighbors have a negative formal charge, the latter
                    // RDKit✔️✔️:         // influences the charge on ipso
                    // RDKit✔️✔️:         if (nbrFormalCharge < 0.0) {
                    if nbr_formal_charge < 0.0 {
                        // RDKit✔️✔️:           q0 += (nbrFormalCharge / (2.0 * (double)(nbrAtom->getDegree())));
                        let nbr_degree = mol.adjacency.neighbors_of(neighbor.atom_index).len();
                        q0 += nbr_formal_charge / (2.0 * nbr_degree as f64);
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:       }
                }
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     // there is a special case for anionic divalent nitrogen
            // RDKit✔️✔️:     // with positively charged neighbor
            // RDKit✔️✔️:     if (atomType == 62) {
            if atom_type == 62 {
                // RDKit✔️✔️:       boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                // RDKit✔️✔️:       for (; nbrIdx != endNbrs; ++nbrIdx) {
                for neighbor in mol.adjacency.neighbors_of(idx) {
                    // RDKit✔️✔️:         const Atom *nbrAtom = mol[*nbrIdx];
                    // RDKit✔️✔️:         nbrFormalCharge = this->getMMFFFormalCharge(nbrAtom->getIdx());
                    let nbr_formal_charge = self.atom_properties[neighbor.atom_index].formal_charge;
                    // RDKit✔️✔️:         if (nbrFormalCharge > 0.0) {
                    if nbr_formal_charge > 0.0 {
                        // RDKit✔️✔️:           q0 -= (nbrFormalCharge / 2.0);
                        q0 -= nbr_formal_charge / 2.0;
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:       }
                }
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     // loop over neighbors
            // RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
            // RDKit✔️✔️:     for (; nbrIdx != endNbrs; ++nbrIdx) {
            for neighbor in mol.adjacency.neighbors_of(idx) {
                // RDKit✔️✔️:       const Atom *nbrAtom = mol[*nbrIdx];
                // RDKit✔️✔️:       const Bond *bond =
                // RDKit✔️✔️:           mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx());
                let bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                // RDKit✔️✔️:       // we need to determine the sign of bond charge
                // RDKit✔️✔️:       // increments depending on the bonding relationship
                // RDKit✔️✔️:       // i.e. we have parameters for [a,b] bonds
                // RDKit✔️✔️:       // but it depends whether ipso is a or b
                // RDKit✔️✔️:       unsigned int nbrAtomType = this->getMMFFAtomType(nbrAtom->getIdx());
                let nbr_atom_type = self.atom_properties[neighbor.atom_index].atom_type;
                // RDKit✔️✔️:       unsigned int bondType = this->getMMFFBondType(bond);
                let bond_type = self.get_mmff_bond_type(bond, mmff_prop)?;
                // RDKit✔️✔️:       mmffChgParams =
                // RDKit✔️✔️:           mmffChg->getMMFFChgParams(bondType, atomType, nbrAtomType);
                let mmff_chg_params =
                    mmff_chg.get(bond_type, u32::from(atom_type), u32::from(nbr_atom_type));
                // RDKit✔️✔️:       sumPartialCharge +=
                // RDKit✔️✔️:           (mmffChgParams.second
                // RDKit✔️✔️:                ? (double)(mmffChgParams.first) * ((mmffChgParams.second)->bci)
                // RDKit✔️✔️:                : ((*mmffPBCI)(atomType)->pbci -
                // RDKit✔️✔️:                   (*mmffPBCI)(nbrAtomType)->pbci));
                sum_partial_charge += if let Some(chg_params) = mmff_chg_params.1 {
                    f64::from(mmff_chg_params.0) * chg_params.bci
                } else {
                    let nbr_pbci = mmff_pbci.get(u32::from(nbr_atom_type)).ok_or(
                        MmffMolPropertiesError::AtomTypePbciMissing {
                            atom_type: nbr_atom_type,
                        },
                    )?;
                    atom_pbci.pbci - nbr_pbci.pbci
                };
                // RDKit✔️✔️:       nbrFormalCharge = this->getMMFFFormalCharge(nbrAtom->getIdx());
                let nbr_formal_charge = self.atom_properties[neighbor.atom_index].formal_charge;
                // RDKit✔️✔️:       sumFormalCharge += nbrFormalCharge;
                sum_formal_charge += nbr_formal_charge;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     // we compute ipso partial charge according to
            // RDKit✔️✔️:     // equation 15, page 622 MMFF.V paper
            // RDKit✔️✔️:     pChg = (1.0 - M * v) * q0 + v * sumFormalCharge + sumPartialCharge;
            let partial_charge = (1.0 - m * v) * q0 + v * sum_formal_charge + sum_partial_charge;
            // RDKit✔️✔️:     this->setMMFFPartialCharge(atom->getIdx(), pChg);
            self.atom_properties[idx].partial_charge = partial_charge;
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::computeMMFFCharges
        Ok(())
    }

    #[must_use]
    pub const fn is_valid(&self) -> bool {
        self.valid
    }

    #[must_use]
    pub const fn mmff_variant(&self) -> MmffVariant {
        self.variant
    }

    pub fn get_mmff_atom_type(&self, idx: usize) -> Result<u8, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFAtomType (AtomTyper.h:94-98)
        // RDKit✔️✔️: std::uint8_t getMMFFAtomType(const unsigned int idx) {
        // RDKit✔️✔️:   URANGE_CHECK(idx, this->d_MMFFAtomPropertiesPtrVect.size());
        let atom_properties =
            self.atom_properties
                .get(idx)
                .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                    atom_index: idx,
                    atoms: self.atom_properties.len(),
                })?;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   return this->d_MMFFAtomPropertiesPtrVect[idx]->mmffAtomType;
        // RDKit✔️✔️: }
        Ok(atom_properties.atom_type)
    }

    pub fn get_mmff_central_atom_prop(
        &self,
        idx: usize,
    ) -> Result<MmffProp, MmffMolPropertiesError> {
        let atom_type = self.get_mmff_atom_type(idx)?;
        let mmff_prop = default_mmff_prop().map_err(|e| MmffParamError::from(e.clone()))?;
        mmff_prop
            .get(u32::from(atom_type))
            .copied()
            .ok_or(MmffMolPropertiesError::AtomTypePropertiesMissing { atom_type })
    }

    pub fn get_mmff_formal_charge(&self, idx: usize) -> Result<f64, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFFormalCharge (AtomTyper.h:99-103)
        // RDKit✔️✔️: double getMMFFFormalCharge(const unsigned int idx) {
        // RDKit✔️✔️:   URANGE_CHECK(idx, this->d_MMFFAtomPropertiesPtrVect.size());
        let atom_properties =
            self.atom_properties
                .get(idx)
                .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                    atom_index: idx,
                    atoms: self.atom_properties.len(),
                })?;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   return this->d_MMFFAtomPropertiesPtrVect[idx]->mmffFormalCharge;
        // RDKit✔️✔️: }
        Ok(atom_properties.formal_charge)
    }

    pub fn get_mmff_partial_charge(&self, idx: usize) -> Result<f64, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFPartialCharge (AtomTyper.h:104-108)
        // RDKit✔️✔️: double getMMFFPartialCharge(const unsigned int idx) {
        // RDKit✔️✔️:   URANGE_CHECK(idx, this->d_MMFFAtomPropertiesPtrVect.size());
        let atom_properties =
            self.atom_properties
                .get(idx)
                .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                    atom_index: idx,
                    atoms: self.atom_properties.len(),
                })?;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   return this->d_MMFFAtomPropertiesPtrVect[idx]->mmffPartialCharge;
        // RDKit✔️✔️: }
        Ok(atom_properties.partial_charge)
    }

    fn get_mmff_bond_type(
        &self,
        bond: &Bond,
        mmff_prop: &MmffPropCollection,
    ) -> Result<u32, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFBondType (AtomTyper.cpp:2456-2477)
        // RDKit✔️✔️: // returns the MMFF bond type of the bond
        // RDKit✔️✔️: unsigned int MMFFMolProperties::getMMFFBondType(const Bond *bond) {
        // RDKit✔️✔️:   PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
        // RDKit✔️✔️:   PRECONDITION(bond, "invalid bond");
        debug_assert!(self.is_valid());

        // RDKit✔️✔️:
        // RDKit✔️✔️:   const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
        // RDKit✔️✔️:   const ForceFields::MMFF::MMFFProp *mmffPropAtom1 =
        // RDKit✔️✔️:       (*mmffProp)(this->getMMFFAtomType(bond->getBeginAtomIdx()));
        let atom_type_1 = self.get_mmff_atom_type(bond.begin().index())?;
        let mmff_prop_atom_1 = mmff_prop.get(u32::from(atom_type_1)).ok_or(
            MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: atom_type_1,
            },
        )?;

        // RDKit✔️✔️:   const ForceFields::MMFF::MMFFProp *mmffPropAtom2 =
        // RDKit✔️✔️:       (*mmffProp)(this->getMMFFAtomType(bond->getEndAtomIdx()));
        let atom_type_2 = self.get_mmff_atom_type(bond.end().index())?;
        let mmff_prop_atom_2 = mmff_prop.get(u32::from(atom_type_2)).ok_or(
            MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: atom_type_2,
            },
        )?;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   // return 1 if the bond is single and the properties for this
        // RDKit✔️✔️:   // single bond match either those of sbmb or aromatic bonds
        // RDKit✔️✔️:   // for this atom pair, 0 if they don't
        // RDKit✔️✔️:   return (unsigned int)(((bond->getBondType() == Bond::SINGLE) &&
        // RDKit✔️✔️:                          ((mmffPropAtom1->sbmb && mmffPropAtom2->sbmb) ||
        // RDKit✔️✔️:                           (mmffPropAtom1->arom && mmffPropAtom2->arom)))
        // RDKit✔️✔️:                             ? 1
        // RDKit✔️✔️:                             : 0);
        // RDKit✔️✔️: }
        Ok(
            if bond.order() == BondOrder::Single
                && ((mmff_prop_atom_1.sbmb != 0 && mmff_prop_atom_2.sbmb != 0)
                    || (mmff_prop_atom_1.arom != 0 && mmff_prop_atom_2.arom != 0))
            {
                1
            } else {
                0
            },
        )
    }

    fn get_mmff_angle_type(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        mmff_prop: &MmffPropCollection,
    ) -> Result<Option<u32>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFAngleType (AtomTyper.cpp:2412-2454)
        // RDKit✔️✔️: // returns the MMFF angle type of the angle formed
        // RDKit✔️✔️: // by atoms with indexes idx1, idx2, idx3
        // RDKit✔️✔️: unsigned int MMFFMolProperties::getMMFFAngleType(const ROMol &mol,
        // RDKit✔️✔️:                                                  const unsigned int idx1,
        // RDKit✔️✔️:                                                  const unsigned int idx2,
        // RDKit✔️✔️:                                                  const unsigned int idx3) {
        // RDKit✔️✔️:   PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
        debug_assert!(self.is_valid());

        // RDKit✔️✔️:
        // RDKit✔️✔️:   // ftp://ftp.wiley.com/public/journals/jcc/suppmat/17/553/MMFF-III_AppendixA.html
        // RDKit✔️✔️:   //
        // RDKit✔️✔️:   // AT[IJK]    Structural significance
        // RDKit✔️✔️:   //--------------------------------------------------------------------------
        // RDKit✔️✔️:   //  0		      The angle i-j-k is a "normal" bond angle
        // RDKit✔️✔️:   //  1 		    Either bond i-j or bond j-k has a bond type of 1
        // RDKit✔️✔️:   //  2		      Bonds i-j and j-k each have bond types of 1; the
        // RDKit✔️✔️:   //  sum is 2. 3		      The angle occurs in a three-membered ring
        // RDKit✔️✔️:   //  4		      The angle occurs in a four-membered ring
        // RDKit✔️✔️:   //  5		      Is in a three-membered ring and the sum of the
        // RDKit✔️✔️:   //  bond types is
        // RDKit✔️✔️:   //  1
        // RDKit✔️✔️:   //  6		      Is in a three-membered ring and the sum of the
        // RDKit✔️✔️:   //  bond types is
        // RDKit✔️✔️:   //  2
        // RDKit✔️✔️:   //  7		      Is in a four-membered ring and the sum of the bond
        // RDKit✔️✔️:   //  types is
        // RDKit✔️✔️:   //  1
        // RDKit✔️✔️:   //  8		      Is in a four-membered ring and the sum of the bond
        // RDKit✔️✔️:   //  types is
        // RDKit✔️✔️:   //  2
        // RDKit✔️✔️:
        let Some(bond_12) = self.bond_between_atom_indices(idx1, idx2) else {
            return Ok(None);
        };
        let Some(bond_23) = self.bond_between_atom_indices(idx2, idx3) else {
            return Ok(None);
        };

        // RDKit✔️✔️:   unsigned int bondTypeSum =
        // RDKit✔️✔️:       this->getMMFFBondType(mol.getBondBetweenAtoms(idx1, idx2)) +
        // RDKit✔️✔️:       this->getMMFFBondType(mol.getBondBetweenAtoms(idx2, idx3));
        let bond_type_sum = self.get_mmff_bond_type(bond_12, mmff_prop)?
            + self.get_mmff_bond_type(bond_23, mmff_prop)?;
        // RDKit✔️✔️:   unsigned int angleType = bondTypeSum;
        let mut angle_type = bond_type_sum;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int size = isAngleInRingOfSize3or4(mol, idx1, idx2, idx3);
        // RDKit✔️✔️:   if (size) {
        if let Some(size) = self.angle_ring_size_3_or_4(idx1, idx2, idx3) {
            // RDKit✔️✔️:     angleType = size;
            angle_type = size;
            // RDKit✔️✔️:     if (bondTypeSum) {
            if bond_type_sum != 0 {
                // RDKit✔️✔️:       angleType += (bondTypeSum + size - 2);
                angle_type += bond_type_sum + size - 2;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return angleType;
        // RDKit✔️✔️: }
        Ok(Some(angle_type))
    }

    pub fn get_mmff_bond_stretch_params(
        &self,
        idx1: usize,
        idx2: usize,
    ) -> Result<Option<(u32, MmffBond)>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFBondStretchParams (AtomTyper.cpp:3490-3517)
        // RDKit✔️✔️: bool MMFFMolProperties::getMMFFBondStretchParams(
        // RDKit✔️✔️:     const ROMol &mol, const unsigned int idx1, const unsigned int idx2,
        // RDKit✔️✔️:     unsigned int &bondType, MMFFBond &mmffBondStretchParams) {
        // RDKit✔️✔️:   const MMFFBondCollection *mmffBond = DefaultParameters::getMMFFBond();
        let mmff_bond = default_mmff_bond().map_err(|e| MmffParamError::from(e.clone()))?;
        let mmff_prop = default_mmff_prop().map_err(|e| MmffParamError::from(e.clone()))?;

        // RDKit✔️✔️:   bool res = false;
        // RDKit✔️✔️:   if (isValid()) {
        if !self.is_valid() {
            // RDKit returns false without touching output parameters.
            return Ok(None);
        }

        // RDKit✔️✔️:     unsigned int iAtomType = getMMFFAtomType(idx1);
        let i_atom_type = u32::from(self.get_mmff_atom_type(idx1)?);
        // RDKit✔️✔️:     unsigned int jAtomType = getMMFFAtomType(idx2);
        let j_atom_type = u32::from(self.get_mmff_atom_type(idx2)?);
        // RDKit✔️✔️:     const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
        // RDKit✔️✔️:     if (bond) {
        let Some(bond) = self.bond_between_atom_indices(idx1, idx2) else {
            return Ok(None);
        };

        // RDKit✔️✔️:       bondType = getMMFFBondType(bond);
        let bond_type = self.get_mmff_bond_type(bond, mmff_prop)?;
        // RDKit✔️✔️:       bool areMMFFBondParamsEmpirical = false;
        // RDKit✔️✔️:       const MMFFBond *mmffBondParams =
        // RDKit✔️✔️:           (*mmffBond)(bondType, iAtomType, jAtomType);
        if let Some(mmff_bond_params) = mmff_bond.get(bond_type, i_atom_type, j_atom_type) {
            // RDKit✔️✔️:       if (mmffBondParams) {
            // RDKit✔️✔️:         mmffBondStretchParams = *mmffBondParams;
            // RDKit✔️✔️:         if (areMMFFBondParamsEmpirical) {
            // RDKit✔️✔️:           delete mmffBondParams;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         res = true;
            return Ok(Some((bond_type, *mmff_bond_params)));
        }

        // RDKit❗✔️:       if (!mmffBondParams) {
        // RDKit❗✔️:         mmffBondParams = getMMFFBondStretchEmpiricalRuleParams(mol, bond);
        // RDKit❗✔️:         areMMFFBondParamsEmpirical = true;
        // RDKit❗✔️:       }
        let mmff_bond_params =
            self.get_mmff_bond_stretch_empirical_rule_params(idx1, idx2, mmff_prop)?;
        Ok(Some((bond_type, mmff_bond_params)))
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
    }

    fn get_mmff_bond_stretch_empirical_rule_params(
        &self,
        idx1: usize,
        idx2: usize,
        mmff_prop: &MmffPropCollection,
    ) -> Result<MmffBond, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFBondStretchEmpiricalRuleParams (AtomTyper.cpp:2573-2729)
        // RDKit❗✔️: // empirical rule to compute bond stretching parameters if
        // RDKit❗✔️: // tabulated parameters could not be found. The returned
        // RDKit❗✔️: // pointer to a MMFFBond object must be freed by the caller
        // RDKit❗✔️: const ForceFields::MMFF::MMFFBond *
        // RDKit❗✔️: MMFFMolProperties::getMMFFBondStretchEmpiricalRuleParams(const ROMol &,
        // RDKit❗✔️:                                                          const Bond *bond) {
        // RDKit❗✔️:   PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
        debug_assert!(self.is_valid());

        // RDKit❗✔️:   const MMFFBond *mmffBndkParams;
        // RDKit❗✔️:   const MMFFHerschbachLaurie *mmffHerschbachLaurieParams;
        // RDKit❗✔️:   const MMFFProp *mmffAtomPropParams[2];
        // RDKit❗✔️:   const MMFFCovRadPauEle *mmffAtomCovRadPauEleParams[2];
        // RDKit❗✔️:   const MMFFBndkCollection *mmffBndk = DefaultParameters::getMMFFBndk();
        // RDKit❗✔️:   const MMFFHerschbachLaurieCollection *mmffHerschbachLaurie =
        // RDKit❗✔️:       DefaultParameters::getMMFFHerschbachLaurie();
        // RDKit❗✔️:   const MMFFCovRadPauEleCollection *mmffCovRadPauEle =
        // RDKit❗✔️:       DefaultParameters::getMMFFCovRadPauEle();
        // RDKit❗✔️:   const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
        let atomic_num_1 = u32::from(self.topology.atoms.as_slice()[idx1].atomic_number());
        let atomic_num_2 = u32::from(self.topology.atoms.as_slice()[idx2].atomic_number());
        let atom_type_1 = u32::from(self.get_mmff_atom_type(idx1)?);
        let atom_type_2 = u32::from(self.get_mmff_atom_type(idx2)?);

        // RDKit❗✔️:   unsigned int atomicNum1 = bond->getBeginAtom()->getAtomicNum();
        // RDKit❗✔️:   unsigned int atomicNum2 = bond->getEndAtom()->getAtomicNum();
        // RDKit❗✔️:   mmffBndkParams = (*mmffBndk)(atomicNum1, atomicNum2);
        let mmff_bndk_params = default_mmff_bndk()
            .map_err(|e| MmffParamError::from(e.clone()))?
            .get(atomic_num_1 as i32, atomic_num_2 as i32);
        // RDKit❗✔️:   mmffAtomCovRadPauEleParams[0] = (*mmffCovRadPauEle)(atomicNum1);
        // RDKit❗✔️:   mmffAtomCovRadPauEleParams[1] = (*mmffCovRadPauEle)(atomicNum2);
        let cov_rad_pau_ele_1 = default_mmff_cov_rad_pau_ele()
            .map_err(|e| MmffParamError::from(e.clone()))?
            .get(atomic_num_1)
            .ok_or(MmffParamError::MissingDefaultData {
                symbol: "defaultMMFFCovRadPauEle atom 1",
            })?;
        let cov_rad_pau_ele_2 = default_mmff_cov_rad_pau_ele()
            .map_err(|e| MmffParamError::from(e.clone()))?
            .get(atomic_num_2)
            .ok_or(MmffParamError::MissingDefaultData {
                symbol: "defaultMMFFCovRadPauEle atom 2",
            })?;
        // RDKit❗✔️:   mmffAtomPropParams[0] =
        // RDKit❗✔️:       (*mmffProp)(this->getMMFFAtomType(bond->getBeginAtomIdx()));
        // RDKit❗✔️:   mmffAtomPropParams[1] =
        // RDKit❗✔️:       (*mmffProp)(this->getMMFFAtomType(bond->getEndAtomIdx()));
        mmff_prop
            .get(atom_type_1)
            .ok_or(MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: atom_type_1 as u8,
            })?;
        mmff_prop
            .get(atom_type_2)
            .ok_or(MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: atom_type_2 as u8,
            })?;

        // RDKit❗✔️:   PRECONDITION(mmffAtomCovRadPauEleParams[0],
        // RDKit❗✔️:                "covalent radius/Pauling electronegativity parameters for atom "
        // RDKit❗✔️:                "1 not found");
        // RDKit❗✔️:   PRECONDITION(mmffAtomCovRadPauEleParams[1],
        // RDKit❗✔️:                "covalent radius/Pauling electronegativity parameters for atom "
        // RDKit❗✔️:                "2 not found");
        // RDKit❗✔️:   PRECONDITION(mmffAtomPropParams[0],
        // RDKit❗✔️:                "property parameters for atom 1 not found");
        // RDKit❗✔️:   PRECONDITION(mmffAtomPropParams[1],
        // RDKit❗✔️:                "property parameters for atom 2 not found");
        // RDKit❗✔️:   auto *mmffBondParams = new ForceFields::MMFF::MMFFBond();
        // RDKit❗✔️:   const double c = (((atomicNum1 == 1) || (atomicNum2 == 1)) ? 0.050 : 0.085);
        let c = if atomic_num_1 == 1 || atomic_num_2 == 1 {
            0.050
        } else {
            0.085
        };
        // RDKit❗✔️:   const double n = 1.4;
        let n = 1.4_f64;
        // RDKit❗✔️: #if 0
        // RDKit❗✔️:       const double delta = 0.008;
        // RDKit❗✔️: #endif
        // RDKit❗✔️: #if 1
        // RDKit❗✔️:   const double delta = 0.0;
        // RDKit❗✔️: #endif
        let delta = 0.0_f64;
        // RDKit❗✔️:   double r0_i[2];
        // RDKit❗✔️:   // MMFF.V, page 625
        // RDKit❗✔️:   for (unsigned int i = 0; i < 2; ++i) {
        // RDKit❗✔️:     r0_i[i] = mmffAtomCovRadPauEleParams[i]->r0;
        // RDKit❗✔️:   }
        let r0_i = [cov_rad_pau_ele_1.r0, cov_rad_pau_ele_2.r0];

        // RDKit❗✔️: // the part of the empirical rule concerning H
        // RDKit❗✔️: // parameters appears not to be used - tests are
        // RDKit❗✔️: // passed only in its absence, hence it is
        // RDKit❗✔️: // currently excluded
        // RDKit❗✔️: #if 0
        // RDKit❗✔️:         switch (mmffAtomPropParams[i]->mltb) {
        // RDKit❗✔️:           case 1:
        // RDKit❗✔️:           case 2:
        // RDKit❗✔️:             H_i[i] = 2;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:
        // RDKit❗✔️:           case 3:
        // RDKit❗✔️:             H_i[i] = 1;
        // RDKit❗✔️:
        // RDKit❗✔️:           default:
        // RDKit❗✔️:             H_i[i] = 3;
        // RDKit❗✔️:         }
        // RDKit❗✔️: #endif
        // RDKit❗✔️:   }
        // RDKit❗✔️: // also the part of the empirical rule concerning BO
        // RDKit❗✔️: // parameters appears not to be used - tests are
        // RDKit❗✔️: // passed only in its absence, hence it is
        // RDKit❗✔️: // currently excluded
        // RDKit❗✔️: #if 0
        // RDKit❗✔️:       unsigned int BO_ij = (unsigned int)(bond->getBondTypeAsDouble());
        // RDKit❗✔️:       if ((mmffAtomPropParams[0]->mltb == 1)
        // RDKit❗✔️:         && (mmffAtomPropParams[1]->mltb == 1)) {
        // RDKit❗✔️:         BO_ij = 4;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (((mmffAtomPropParams[0]->mltb == 1)
        // RDKit❗✔️:         && (mmffAtomPropParams[1]->mltb == 2))
        // RDKit❗✔️:         || ((mmffAtomPropParams[0]->mltb == 2)
        // RDKit❗✔️:         && (mmffAtomPropParams[1]->mltb == 1))) {
        // RDKit❗✔️:         BO_ij = 5;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (areAtomsInSameAromaticRing(mol,
        // RDKit❗✔️:         bond->getBeginAtomIdx(), bond->getEndAtomIdx())) {
        // RDKit❗✔️:         BO_ij = (((mmffAtomPropParams[0]->pilp == 0)
        // RDKit❗✔️:           && (mmffAtomPropParams[1]->pilp == 0)) ? 4 : 5);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (BO_ij == 1) {
        // RDKit❗✔️:         for (unsigned int i = 0; i < 2; ++i) {
        // RDKit❗✔️:           std::cout << "H" << i << "=" << H_i[i] << std::endl;
        // RDKit❗✔️:           switch (H_i[i]) {
        // RDKit❗✔️:             case 1:
        // RDKit❗✔️:               r0_i[i] -= 0.08;
        // RDKit❗✔️:             break;
        // RDKit❗✔️:
        // RDKit❗✔️:             case 2:
        // RDKit❗✔️:               r0_i[i] -= 0.03;
        // RDKit❗✔️:             break;
        // RDKit❗✔️:           }
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:       else {
        // RDKit❗✔️:         double dec = 0.0;
        // RDKit❗✔️:         switch (BO_ij) {
        // RDKit❗✔️:           case 5:
        // RDKit❗✔️:             dec = 0.04;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:
        // RDKit❗✔️:           case 4:
        // RDKit❗✔️:             dec = 0.075;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:
        // RDKit❗✔️:           case 3:
        // RDKit❗✔️:             dec = 0.17;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:
        // RDKit❗✔️:           case 2:
        // RDKit❗✔️:             dec = 0.10;
        // RDKit❗✔️:           break;
        // RDKit❗✔️:         }
        // RDKit❗✔️:         r0_i[0] -= dec;
        // RDKit❗✔️:         r0_i[1] -= dec;
        // RDKit❗✔️:       }
        // RDKit❗✔️: #endif
        // Both adjustment blocks are disabled in the pinned source, so the
        // executable Rust path intentionally uses the unadjusted radii.
        // RDKit❗✔️:   // equation (18) - MMFF.V, page 625
        // RDKit❗✔️:   mmffBondParams->r0 = (r0_i[0] + r0_i[1] -
        // RDKit❗✔️:                         c * pow(fabs(mmffAtomCovRadPauEleParams[0]->chi -
        // RDKit❗✔️:                                      mmffAtomCovRadPauEleParams[1]->chi),
        // RDKit❗✔️:                                 n) -
        // RDKit❗✔️:                         delta);
        let r0 = r0_i[0] + r0_i[1]
            - c * (cov_rad_pau_ele_1.chi - cov_rad_pau_ele_2.chi)
                .abs()
                .powf(n)
            - delta;
        // RDKit❗✔️:   if (mmffBndkParams) {
        let kb = if let Some(mmff_bndk_params) = mmff_bndk_params {
            // RDKit❗✔️:     // equation (19) - MMFF.V, page 625
            // RDKit❗✔️:     double coeff = mmffBndkParams->r0 / mmffBondParams->r0;
            // RDKit❗✔️:     double coeff2 = coeff * coeff;
            // RDKit❗✔️:     double coeff6 = coeff2 * coeff2 * coeff2;
            // RDKit❗✔️:     mmffBondParams->kb = mmffBndkParams->kb * coeff6;
            let coeff = mmff_bndk_params.r0 / r0;
            let coeff_2 = coeff * coeff;
            let coeff_6 = coeff_2 * coeff_2 * coeff_2;
            mmff_bndk_params.kb * coeff_6
        // RDKit❗✔️:   } else {
        } else {
            // RDKit❗✔️:     // MMFF.V, page 627
            // RDKit❗✔️:     // Herschbach-Laurie version of Badger's rule
            // RDKit❗✔️:     // J. Chem. Phys. 35, 458 (1961); https://doi.org/10.1063/1.1731952
            // RDKit❗✔️:     // equation (8), page 5
            // RDKit❗✔️:     mmffHerschbachLaurieParams = (*mmffHerschbachLaurie)(
            // RDKit❗✔️:         getPeriodicTableRowHL(atomicNum1), getPeriodicTableRowHL(atomicNum2));
            let herschbach_laurie = default_mmff_herschbach_laurie()
                .map_err(|e| MmffParamError::from(e.clone()))?
                .get(
                    mmff_periodic_table_row_hl(atomic_num_1 as u8) as i32,
                    mmff_periodic_table_row_hl(atomic_num_2 as u8) as i32,
                )
                .ok_or(MmffParamError::MissingDefaultData {
                    symbol: "defaultMMFFHerschbachLaurie row pair",
                })?;
            // RDKit❗✔️:     mmffBondParams->kb =
            // RDKit❗✔️:         pow(10.0, -(mmffBondParams->r0 - mmffHerschbachLaurieParams->a_ij) /
            // RDKit❗✔️:                       mmffHerschbachLaurieParams->d_ij);
            10.0_f64.powf(-(r0 - herschbach_laurie.a_ij) / herschbach_laurie.d_ij)
            // RDKit❗✔️:   }
        };

        // RDKit❗✔️:   return (const ForceFields::MMFF::MMFFBond *)mmffBondParams;
        // RDKit❗✔️: }
        Ok(MmffBond { kb, r0 })
        // END RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFBondStretchEmpiricalRuleParams
    }

    fn get_mmff_angle_bend_empirical_rule_params(
        &self,
        old_mmff_angle_params: Option<&MmffAngle>,
        mmff_prop_params_central_atom: &MmffProp,
        mmff_bond_params_1: &MmffBond,
        mmff_bond_params_2: &MmffBond,
        idx1: usize,
        idx2: usize,
        idx3: usize,
    ) -> MmffAngle {
        // BEGIN RDKIT CPP FUNCTION RDKit::MMFF::getMMFFAngleBendEmpiricalRuleParams (AtomTyper.cpp:2731-2872)
        // RDKit❗✔️: // empirical rule to compute angle bending parameters if
        // RDKit❗✔️: // tabulated parameters could not be found. The returned
        // RDKit❗✔️: // pointer to a MMFFAngle object must be freed by the caller
        // RDKit❗✔️: const ForceFields::MMFF::MMFFAngle *getMMFFAngleBendEmpiricalRuleParams(
        // RDKit❗✔️:     const ROMol &mol, const ForceFields::MMFF::MMFFAngle *oldMMFFAngleParams,
        // RDKit❗✔️:     const ForceFields::MMFF::MMFFProp *mmffPropParamsCentralAtom,
        // RDKit❗✔️:     const ForceFields::MMFF::MMFFBond *mmffBondParams1,
        // RDKit❗✔️:     const ForceFields::MMFF::MMFFBond *mmffBondParams2, unsigned int idx1,
        // RDKit❗✔️:     unsigned int idx2, unsigned int idx3) {
        // RDKit❗✔️:   int atomicNum[3];
        // RDKit❗✔️:   atomicNum[0] = mol.getAtomWithIdx(idx1)->getAtomicNum();
        // RDKit❗✔️:   atomicNum[1] = mol.getAtomWithIdx(idx2)->getAtomicNum();
        // RDKit❗✔️:   atomicNum[2] = mol.getAtomWithIdx(idx3)->getAtomicNum();
        let atomic_num = [
            self.topology.atoms.as_slice()[idx1].atomic_number(),
            self.topology.atoms.as_slice()[idx2].atomic_number(),
            self.topology.atoms.as_slice()[idx3].atomic_number(),
        ];
        // RDKit❗✔️:   auto *mmffAngleParams = new ForceFields::MMFF::MMFFAngle();
        // RDKit❗✔️:   unsigned int ringSize = isAngleInRingOfSize3or4(mol, idx1, idx2, idx3);
        let ring_size = self.angle_ring_size_3_or_4(idx1, idx2, idx3).unwrap_or(0);
        // RDKit❗✔️:   if (!oldMMFFAngleParams) {
        let theta0 = if old_mmff_angle_params.is_none() {
            // RDKit❗✔️:     // angle rest value empirical rule
            // RDKit❗✔️:     mmffAngleParams->theta0 = 120.0;
            let mut theta0 = 120.0;
            // RDKit❗✔️:     switch (mmffPropParamsCentralAtom->crd) {
            match mmff_prop_params_central_atom.crd {
                // RDKit❗✔️:       case 4:
                4 => {
                    // RDKit❗✔️:         // if the central atom has crd = 4
                    // RDKit❗✔️:         mmffAngleParams->theta0 = 109.45;
                    theta0 = 109.45;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       case 2:
                2 => {
                    // RDKit❗✔️:         // if the central atom is oxygen
                    // RDKit❗✔️:         if (atomicNum[1] == 8) {
                    if atomic_num[1] == 8 {
                        // RDKit❗✔️:           mmffAngleParams->theta0 = 105.0;
                        theta0 = 105.0;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         // if the central atom is linear
                    // RDKit❗✔️:         else if (mmffPropParamsCentralAtom->linh == 1) {
                    } else if mmff_prop_params_central_atom.linh == 1 {
                        // RDKit❗✔️:           mmffAngleParams->theta0 = 180.0;
                        theta0 = 180.0;
                        // RDKit❗✔️:         }
                    }
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       case 3:
                3 => {
                    // RDKit❗✔️:         if ((mmffPropParamsCentralAtom->val == 3) &&
                    // RDKit❗✔️:             (mmffPropParamsCentralAtom->mltb == 0)) {
                    if mmff_prop_params_central_atom.val == 3
                        && mmff_prop_params_central_atom.mltb == 0
                    {
                        // RDKit❗✔️:           // if the central atom is nitrogen
                        // RDKit❗✔️:           if (atomicNum[1] == 7) {
                        if atomic_num[1] == 7 {
                            // RDKit❗✔️:             mmffAngleParams->theta0 = 107.0;
                            theta0 = 107.0;
                        // RDKit❗✔️:           } else {
                        } else {
                            // RDKit❗✔️:             mmffAngleParams->theta0 = 92.0;
                            theta0 = 92.0;
                            // RDKit❗✔️:           }
                        }
                        // RDKit❗✔️:         }
                    }
                    // RDKit❗✔️:         break;
                }
                _ => {}
            }
            // RDKit❗✔️:     }
            // RDKit❗✔️:     if (ringSize == 3) {
            if ring_size == 3 {
                // RDKit❗✔️:       mmffAngleParams->theta0 = 60.0;
                theta0 = 60.0;
            // RDKit❗✔️:     } else if (ringSize == 4) {
            } else if ring_size == 4 {
                // RDKit❗✔️:       mmffAngleParams->theta0 = 90.0;
                theta0 = 90.0;
                // RDKit❗✔️:     }
            }
            theta0
        // RDKit❗✔️:   } else {
        } else {
            // RDKit❗✔️:     mmffAngleParams->theta0 = oldMMFFAngleParams->theta0;
            let theta0 = old_mmff_angle_params
                .expect("source else branch requires old MMFF angle parameters")
                .theta0;
            // RDKit❗✔️:   }
            theta0
        };
        // RDKit❗✔️:   // angle force constant empirical rule
        // RDKit❗✔️:   double Z[3] = {0.0, 0.0, 0.0};
        let mut z = [0.0; 3];
        // RDKit❗✔️:   double C[3] = {0.0, 0.0, 0.0};
        let mut c = [0.0; 3];
        // RDKit❗✔️:   double beta = 1.75;
        let mut beta = 1.75;
        // RDKit❗✔️:   for (unsigned int i = 0; i < 3; ++i) {
        for i in 0..3 {
            // RDKit❗✔️:     // Table VI - MMFF.V, page 628
            // RDKit❗✔️:     switch (atomicNum[i]) {
            match atomic_num[i] {
                // RDKit❗✔️:       // Hydrogen
                // RDKit❗✔️:       case 1:
                1 => {
                    // RDKit❗✔️:         Z[i] = 1.395;
                    z[i] = 1.395;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Carbon
                // RDKit❗✔️:       case 6:
                6 => {
                    // RDKit❗✔️:         Z[i] = 2.494;
                    z[i] = 2.494;
                    // RDKit❗✔️:         C[i] = 1.016;
                    c[i] = 1.016;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Nitrogen
                // RDKit❗✔️:       case 7:
                7 => {
                    // RDKit❗✔️:         Z[i] = 2.711;
                    z[i] = 2.711;
                    // RDKit❗✔️:         C[i] = 1.113;
                    c[i] = 1.113;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Oxygen
                // RDKit❗✔️:       case 8:
                8 => {
                    // RDKit❗✔️:         Z[i] = 3.045;
                    z[i] = 3.045;
                    // RDKit❗✔️:         C[i] = 1.337;
                    c[i] = 1.337;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Fluorine
                // RDKit❗✔️:       case 9:
                9 => {
                    // RDKit❗✔️:         Z[i] = 2.847;
                    z[i] = 2.847;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Silicon
                // RDKit❗✔️:       case 14:
                14 => {
                    // RDKit❗✔️:         Z[i] = 2.350;
                    z[i] = 2.350;
                    // RDKit❗✔️:         C[i] = 0.811;
                    c[i] = 0.811;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Phosphorus
                // RDKit❗✔️:       case 15:
                15 => {
                    // RDKit❗✔️:         Z[i] = 2.350;
                    z[i] = 2.350;
                    // RDKit❗✔️:         C[i] = 1.068;
                    c[i] = 1.068;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Sulfur
                // RDKit❗✔️:       case 16:
                16 => {
                    // RDKit❗✔️:         Z[i] = 2.980;
                    z[i] = 2.980;
                    // RDKit❗✔️:         C[i] = 1.249;
                    c[i] = 1.249;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Chlorine
                // RDKit❗✔️:       case 17:
                17 => {
                    // RDKit❗✔️:         Z[i] = 2.909;
                    z[i] = 2.909;
                    // RDKit❗✔️:         C[i] = 1.078;
                    c[i] = 1.078;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Bromine
                // RDKit❗✔️:       case 35:
                35 => {
                    // RDKit❗✔️:         Z[i] = 3.017;
                    z[i] = 3.017;
                    // RDKit❗✔️:         break;
                }
                // RDKit❗✔️:       // Iodine
                // RDKit❗✔️:       case 53:
                53 => {
                    // RDKit❗✔️:         Z[i] = 3.086;
                    z[i] = 3.086;
                    // RDKit❗✔️:         break;
                }
                _ => {}
            }
            // RDKit❗✔️:     }
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   double r0_ij = mmffBondParams1->r0;
        let r0_ij = mmff_bond_params_1.r0;
        // RDKit❗✔️:   double r0_jk = mmffBondParams2->r0;
        let r0_jk = mmff_bond_params_2.r0;
        // RDKit❗✔️:   double D =
        // RDKit❗✔️:       (r0_ij - r0_jk) * (r0_ij - r0_jk) / ((r0_ij + r0_jk) * (r0_ij + r0_jk));
        let d = (r0_ij - r0_jk) * (r0_ij - r0_jk) / ((r0_ij + r0_jk) * (r0_ij + r0_jk));
        // RDKit❗✔️:   double theta0_rad = DEG2RAD * mmffAngleParams->theta0;
        let theta0_rad = std::f64::consts::PI / 180.0 * theta0;
        // RDKit❗✔️:   if (ringSize == 4) {
        if ring_size == 4 {
            // RDKit❗✔️:     beta *= 0.85;
            beta *= 0.85;
        // RDKit❗✔️:   } else if (ringSize == 3) {
        } else if ring_size == 3 {
            // RDKit❗✔️:     beta *= 0.05;
            beta *= 0.05;
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   // equation (20) - MMFF.V, page 628
        // RDKit❗✔️:   mmffAngleParams->ka =
        // RDKit❗✔️:       beta * Z[0] * C[1] * Z[2] /
        // RDKit❗✔️:       ((r0_ij + r0_jk) * theta0_rad * theta0_rad * exp(2.0 * D));
        let ka = beta * z[0] * c[1] * z[2]
            / ((r0_ij + r0_jk) * theta0_rad * theta0_rad * (2.0 * d).exp());

        // RDKit❗✔️:   return (const ForceFields::MMFF::MMFFAngle *)mmffAngleParams;
        // RDKit❗✔️: }
        MmffAngle { ka, theta0 }
        // END RDKIT CPP FUNCTION RDKit::MMFF::getMMFFAngleBendEmpiricalRuleParams
    }

    pub fn get_mmff_angle_bend_params(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
    ) -> Result<Option<(u32, MmffAngle)>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFAngleBendParams (AtomTyper.cpp:3520-3562)
        // RDKit✔️✔️: bool MMFFMolProperties::getMMFFAngleBendParams(const ROMol &mol,
        // RDKit✔️✔️:                                                const unsigned int idx1,
        // RDKit✔️✔️:                                                const unsigned int idx2,
        // RDKit✔️✔️:                                                const unsigned int idx3,
        // RDKit✔️✔️:                                                unsigned int &angleType,
        // RDKit✔️✔️:                                                MMFFAngle &mmffAngleBendParams) {
        // RDKit✔️✔️:   bool res = false;
        // RDKit✔️✔️:   if (isValid() && mol.getBondBetweenAtoms(idx1, idx2) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx2, idx3)) {
        if !self.is_valid() {
            return Ok(None);
        }
        if self.bond_between_atom_indices(idx1, idx2).is_none()
            || self.bond_between_atom_indices(idx2, idx3).is_none()
        {
            return Ok(None);
        }

        // RDKit✔️✔️:     const MMFFAngleCollection *mmffAngle = DefaultParameters::getMMFFAngle();
        let mmff_angle = default_mmff_angle().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
        let mmff_prop = default_mmff_prop().map_err(|e| MmffParamError::from(e.clone()))?;
        let mmff_def = default_mmff_def().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     unsigned int idx[3] = {idx1, idx2, idx3};
        // RDKit❌❌:     MMFFBond mmffBondParams[2];
        // RDKit✔️✔️:     unsigned int atomType[3];
        // RDKit✔️✔️:     unsigned int i;
        // RDKit✔️✔️:     angleType = getMMFFAngleType(mol, idx1, idx2, idx3);
        let Some(angle_type) = self.get_mmff_angle_type(idx1, idx2, idx3, mmff_prop)? else {
            return Ok(None);
        };
        // RDKit✔️✔️:     bool areMMFFAngleParamsEmpirical = false;
        // RDKit✔️✔️:     for (i = 0; i < 3; ++i) {
        // RDKit✔️✔️:       atomType[i] = getMMFFAtomType(idx[i]);
        // RDKit✔️✔️:     }
        let atom_type_1 = u32::from(self.get_mmff_atom_type(idx1)?);
        let atom_type_2 = u32::from(self.get_mmff_atom_type(idx2)?);
        let atom_type_3 = u32::from(self.get_mmff_atom_type(idx3)?);

        // RDKit✔️✔️:     const MMFFAngle *mmffAngleParams =
        // RDKit✔️✔️:         (*mmffAngle)(DefaultParameters::getMMFFDef(), angleType, atomType[0],
        // RDKit✔️✔️:                      atomType[1], atomType[2]);
        let mmff_angle_params = mmff_angle
            .get(mmff_def, angle_type, atom_type_1, atom_type_2, atom_type_3)
            .map_err(MmffParamError::from)?;
        // RDKit✔️✔️:     const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(atomType[1]);
        let mmff_prop_params_central_atom = mmff_prop.get(atom_type_2).ok_or(
            MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: self.get_mmff_atom_type(idx2)?,
            },
        )?;

        if let Some(mmff_angle_params) = mmff_angle_params {
            if !is_double_zero(mmff_angle_params.ka) {
                // RDKit✔️✔️:     if (mmffAngleParams) {
                // RDKit✔️✔️:       mmffAngleBendParams = *mmffAngleParams;
                // RDKit✔️✔️:       res = true;
                // RDKit✔️✔️:       if (areMMFFAngleParamsEmpirical) {
                // RDKit✔️✔️:         delete mmffAngleParams;
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:     }
                return Ok(Some((angle_type, *mmff_angle_params)));
            }
        }

        // RDKit❗✔️:     if ((!mmffAngleParams) || (isDoubleZero(mmffAngleParams->ka))) {
        // RDKit❗✔️:       areMMFFAngleParamsEmpirical = true;
        // RDKit❗✔️:       for (i = 0; areMMFFAngleParamsEmpirical && (i < 2); ++i) {
        // RDKit❗✔️:         unsigned int bondType;
        // RDKit❗✔️:         areMMFFAngleParamsEmpirical = getMMFFBondStretchParams(
        // RDKit❗✔️:             mol, idx[i], idx[i + 1], bondType, mmffBondParams[i]);
        // RDKit❗✔️:       }
        let Some((_bond_type_1, mmff_bond_params_1)) =
            self.get_mmff_bond_stretch_params(idx1, idx2)?
        else {
            return Ok(None);
        };
        let Some((_bond_type_2, mmff_bond_params_2)) =
            self.get_mmff_bond_stretch_params(idx2, idx3)?
        else {
            return Ok(None);
        };
        // RDKit❗✔️:       if (areMMFFAngleParamsEmpirical) {
        // RDKit❗✔️:         mmffAngleParams = getMMFFAngleBendEmpiricalRuleParams(
        // RDKit❗✔️:             mol, mmffAngleParams, mmffPropParamsCentralAtom, &mmffBondParams[0],
        // RDKit❗✔️:             &mmffBondParams[1], idx[0], idx[1], idx[2]);
        // RDKit❗✔️:       }
        let mmff_angle_params = self.get_mmff_angle_bend_empirical_rule_params(
            mmff_angle_params,
            mmff_prop_params_central_atom,
            &mmff_bond_params_1,
            &mmff_bond_params_2,
            idx1,
            idx2,
            idx3,
        );
        // RDKit❗✔️:     }
        Ok(Some((angle_type, mmff_angle_params)))
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
    }

    pub fn get_mmff_stretch_bend_params(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
    ) -> Result<Option<(u32, MmffStbn, [MmffBond; 2], MmffAngle)>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFStretchBendParams (AtomTyper.cpp:3568-3625)
        // RDKit✔️✔️: bool MMFFMolProperties::getMMFFStretchBendParams(
        // RDKit✔️✔️:     const ROMol &mol, const unsigned int idx1, const unsigned int idx2,
        // RDKit✔️✔️:     const unsigned int idx3, unsigned int &stretchBendType,
        // RDKit✔️✔️:     MMFFStbn &mmffStretchBendParams, MMFFBond mmffBondStretchParams[2],
        // RDKit✔️✔️:     MMFFAngle &mmffAngleBendParams) {
        // RDKit✔️✔️:   bool res = false;
        // RDKit✔️✔️:   if (isValid()) {
        if !self.is_valid() {
            return Ok(None);
        }

        // RDKit✔️✔️:     const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
        let mmff_prop = default_mmff_prop().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     const MMFFStbnCollection *mmffStbn = DefaultParameters::getMMFFStbn();
        let mmff_stbn = default_mmff_stbn().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     const MMFFDfsbCollection *mmffDfsb = DefaultParameters::getMMFFDfsb();
        let mmff_dfsb = default_mmff_dfsb().map_err(|e| MmffParamError::from(e.clone()))?;

        // RDKit✔️✔️:     unsigned int idx[3] = {idx1, idx2, idx3};
        // RDKit✔️✔️:     unsigned int atomType[3];
        // RDKit✔️✔️:     unsigned int bondType[2];
        // RDKit✔️✔️:     unsigned int angleType;
        let idx = [idx1, idx2, idx3];

        // RDKit✔️✔️:     const MMFFProp *mmffPropParamsCentralAtom =
        // RDKit✔️✔️:         (*mmffProp)(getMMFFAtomType(idx[1]));
        let central_atom_type = self.get_mmff_atom_type(idx[1])?;
        let mmff_prop_params_central_atom = mmff_prop.get(u32::from(central_atom_type)).ok_or(
            MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: central_atom_type,
            },
        )?;

        // RDKit✔️✔️:     if (!(mmffPropParamsCentralAtom->linh)) {
        if mmff_prop_params_central_atom.linh != 0 {
            // RDKit returns false for linear-H central atom rows without
            // touching the output stretch-bend parameters.
            return Ok(None);
        }

        // RDKit✔️✔️:       res = true;
        // RDKit✔️✔️:       unsigned int i = 0;
        // RDKit✔️✔️:       for (i = 0; i < 3; ++i) {
        // RDKit✔️✔️:         atomType[i] = getMMFFAtomType(idx[i]);
        // RDKit✔️✔️:       }
        let atom_type = [
            u32::from(self.get_mmff_atom_type(idx[0])?),
            u32::from(central_atom_type),
            u32::from(self.get_mmff_atom_type(idx[2])?),
        ];

        // RDKit✔️✔️:       for (i = 0; res && (i < 2); ++i) {
        // RDKit✔️✔️:         res = getMMFFBondStretchParams(mol, idx[i], idx[i + 1], bondType[i],
        // RDKit✔️✔️:                                        mmffBondStretchParams[i]);
        // RDKit✔️✔️:       }
        let Some((bond_type_1, bond_params_1)) =
            self.get_mmff_bond_stretch_params(idx[0], idx[1])?
        else {
            return Ok(None);
        };
        let Some((bond_type_2, bond_params_2)) =
            self.get_mmff_bond_stretch_params(idx[1], idx[2])?
        else {
            return Ok(None);
        };
        let bond_type = [bond_type_1, bond_type_2];

        // RDKit✔️✔️:       if (res) {
        // RDKit✔️✔️:         res = getMMFFAngleBendParams(mol, idx1, idx2, idx3, angleType,
        // RDKit✔️✔️:                                      mmffAngleBendParams);
        // RDKit✔️✔️:       }
        let Some((angle_type, angle_params)) = self.get_mmff_angle_bend_params(idx1, idx2, idx3)?
        else {
            return Ok(None);
        };

        // RDKit✔️✔️:       std::pair<bool, const MMFFStbn *> mmffStbnParams;
        // RDKit✔️✔️:       if (res) {
        // RDKit✔️✔️:         stretchBendType = getMMFFStretchBendType(
        // RDKit✔️✔️:             angleType, (atomType[0] <= atomType[2]) ? bondType[0] : bondType[1],
        // RDKit✔️✔️:             (atomType[0] < atomType[2]) ? bondType[1] : bondType[0]);
        let stretch_bend_type = mmff_stretch_bend_type(
            angle_type,
            if atom_type[0] <= atom_type[2] {
                bond_type[0]
            } else {
                bond_type[1]
            },
            if atom_type[0] < atom_type[2] {
                bond_type[1]
            } else {
                bond_type[0]
            },
        );

        // RDKit✔️✔️:         mmffStbnParams = mmffStbn->getMMFFStbnParams(
        // RDKit✔️✔️:             stretchBendType, bondType[0], bondType[1], atomType[0], atomType[1],
        // RDKit✔️✔️:             atomType[2]);
        let mut stbn_params = mmff_stbn.get(
            stretch_bend_type,
            bond_type[0],
            bond_type[1],
            atom_type[0],
            atom_type[1],
            atom_type[2],
        );

        // RDKit✔️✔️:         if (!(mmffStbnParams.second)) {
        // RDKit✔️✔️:           mmffStbnParams = mmffDfsb->getMMFFDfsbParams(
        // RDKit✔️✔️:               getPeriodicTableRow(mol.getAtomWithIdx(idx1)->getAtomicNum()),
        // RDKit✔️✔️:               getPeriodicTableRow(mol.getAtomWithIdx(idx2)->getAtomicNum()),
        // RDKit✔️✔️:               getPeriodicTableRow(mol.getAtomWithIdx(idx3)->getAtomicNum()));
        // RDKit✔️✔️:         }
        if stbn_params.1.is_none() {
            let row_1 =
                mmff_periodic_table_row(self.topology.atoms.as_slice()[idx1].atomic_number());
            let row_2 =
                mmff_periodic_table_row(self.topology.atoms.as_slice()[idx2].atomic_number());
            let row_3 =
                mmff_periodic_table_row(self.topology.atoms.as_slice()[idx3].atomic_number());
            stbn_params = mmff_dfsb.get(row_1, row_2, row_3);
        }

        // RDKit✔️✔️:         res = (mmffStbnParams.second &&
        // RDKit✔️✔️:                !(isDoubleZero((mmffStbnParams.second)->kbaIJK) &&
        // RDKit✔️✔️:                  isDoubleZero((mmffStbnParams.second)->kbaKJI)));
        // RDKit✔️✔️:       }
        let Some(stbn_param) = stbn_params.1 else {
            return Ok(None);
        };
        if is_double_zero(stbn_param.kba_ijk) && is_double_zero(stbn_param.kba_kji) {
            return Ok(None);
        }

        // RDKit✔️✔️:       if (res) {
        // RDKit✔️✔️:         if (mmffStbnParams.first) {
        // RDKit✔️✔️:           mmffStretchBendParams.kbaIJK = (mmffStbnParams.second)->kbaKJI;
        // RDKit✔️✔️:           mmffStretchBendParams.kbaKJI = (mmffStbnParams.second)->kbaIJK;
        // RDKit✔️✔️:         } else {
        // RDKit✔️✔️:           mmffStretchBendParams = *(mmffStbnParams.second);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        let stretch_bend_params = if stbn_params.0 {
            MmffStbn {
                kba_ijk: stbn_param.kba_kji,
                kba_kji: stbn_param.kba_ijk,
            }
        } else {
            *stbn_param
        };

        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Ok(Some((
            stretch_bend_type,
            stretch_bend_params,
            [bond_params_1, bond_params_2],
            angle_params,
        )))
    }

    fn get_mmff_torsion_type(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        mmff_prop: &MmffPropCollection,
    ) -> Result<(u32, u32), MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFTorsionType (AtomTyper.cpp:2529-2575)
        // RDKit✔️✔️: // given a dihedral angle formed by 4 atoms with indexes
        // RDKit✔️✔️: // idx1, idx2, idx3, idx4, it returns a std::pair whose first element
        // RDKit✔️✔️: // is the principal torsion type, and the second is the secondary
        // RDKit✔️✔️: // torsion type, to be used only if parameters could not be found
        // RDKit✔️✔️: // (empirically found - this is not mentioned either in MMFF.IV
        // RDKit✔️✔️: // nor in MMFF.V)
        // RDKit✔️✔️: const std::pair<unsigned int, unsigned int>
        // RDKit✔️✔️: MMFFMolProperties::getMMFFTorsionType(const ROMol &mol, const unsigned int idx1,
        // RDKit✔️✔️:                                       const unsigned int idx2,
        // RDKit✔️✔️:                                       const unsigned int idx3,
        // RDKit✔️✔️:                                       const unsigned int idx4) {
        // RDKit✔️✔️:   PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
        debug_assert!(self.is_valid());

        // RDKit✔️✔️:
        // RDKit✔️✔️:   const Bond *bondJK = mol.getBondBetweenAtoms(idx2, idx3);
        let Some(bond_jk) = self.bond_between_atom_indices(idx2, idx3) else {
            return Ok((0, 0));
        };
        // RDKit✔️✔️:   unsigned int bondTypeIJ =
        // RDKit✔️✔️:       this->getMMFFBondType(mol.getBondBetweenAtoms(idx1, idx2));
        let Some(bond_ij) = self.bond_between_atom_indices(idx1, idx2) else {
            return Ok((0, 0));
        };
        let bond_type_ij = self.get_mmff_bond_type(bond_ij, mmff_prop)?;
        // RDKit✔️✔️:   unsigned int bondTypeJK = this->getMMFFBondType(bondJK);
        let bond_type_jk = self.get_mmff_bond_type(bond_jk, mmff_prop)?;
        // RDKit✔️✔️:   unsigned int bondTypeKL =
        // RDKit✔️✔️:       this->getMMFFBondType(mol.getBondBetweenAtoms(idx3, idx4));
        let Some(bond_kl) = self.bond_between_atom_indices(idx3, idx4) else {
            return Ok((0, 0));
        };
        let bond_type_kl = self.get_mmff_bond_type(bond_kl, mmff_prop)?;
        // RDKit✔️✔️:   unsigned int torsionType = bondTypeJK;
        let mut torsion_type = bond_type_jk;
        // RDKit✔️✔️:   unsigned int secondTorsionType = 0;
        let mut second_torsion_type = 0;

        // RDKit✔️✔️:
        // RDKit✔️✔️:   // according to MMFF.IV page 609 the condition should be as simple as
        // RDKit✔️✔️:   // if ((bondTypeJK == 0) && ((bondTypeIJ == 1) || (bondTypeKL == 1))) {
        // RDKit✔️✔️:   // but CYGUAN01 fails the test, so the following condition was
        // RDKit✔️✔️:   // empirically determined to be the correct one
        // RDKit✔️✔️:   if ((bondTypeJK == 0) && (bondJK->getBondType() == Bond::SINGLE) &&
        // RDKit✔️✔️:       ((bondTypeIJ == 1) || (bondTypeKL == 1))) {
        if bond_type_jk == 0
            && bond_jk.order() == BondOrder::Single
            && (bond_type_ij == 1 || bond_type_kl == 1)
        {
            // RDKit✔️✔️:     torsionType = 2;
            torsion_type = 2;
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   unsigned int size = isTorsionInRingOfSize4or5(mol, idx1, idx2, idx3, idx4);
        let size = self.torsion_ring_size_4_or_5(idx1, idx2, idx3, idx4);
        // RDKit✔️✔️:   // the additional check on the existence of a bond between I and K or J and
        // RDKit✔️✔️:   // L is to avoid assigning torsionType 4 to those torsions in a 4-membered
        // RDKit✔️✔️:   // ring constituted by the fusion of two 3-membered rings, even though it
        // RDKit✔️✔️:   // would be harmless for the energy calculation since parameters for
        // RDKit✔️✔️:   // 4,22,22,22,22 and 0,22,22,22,22 are identical
        // RDKit✔️✔️:   if ((size == 4) && (!(mol.getBondBetweenAtoms(idx1, idx3) ||
        // RDKit✔️✔️:                         mol.getBondBetweenAtoms(idx2, idx4)))) {
        if size == Some(4)
            && self.bond_between_atom_indices(idx1, idx3).is_none()
            && self.bond_between_atom_indices(idx2, idx4).is_none()
        {
            // RDKit✔️✔️:     secondTorsionType = torsionType;
            // RDKit✔️✔️:     torsionType = 4;
            second_torsion_type = torsion_type;
            torsion_type = 4;
        // RDKit✔️✔️:   } else if ((size == 5) && ((this->getMMFFAtomType(idx1) == 1) ||
        // RDKit✔️✔️:                              (this->getMMFFAtomType(idx2) == 1) ||
        // RDKit✔️✔️:                              (this->getMMFFAtomType(idx3) == 1) ||
        // RDKit✔️✔️:                              (this->getMMFFAtomType(idx4) == 1))) {
        } else if size == Some(5)
            && (self.get_mmff_atom_type(idx1)? == 1
                || self.get_mmff_atom_type(idx2)? == 1
                || self.get_mmff_atom_type(idx3)? == 1
                || self.get_mmff_atom_type(idx4)? == 1)
        {
            // RDKit✔️✔️:     secondTorsionType = torsionType;
            // RDKit✔️✔️:     torsionType = 5;
            second_torsion_type = torsion_type;
            torsion_type = 5;
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return std::make_pair(torsionType, secondTorsionType);
        // RDKit✔️✔️: }
        Ok((torsion_type, second_torsion_type))
    }

    fn get_mmff_torsion_empirical_rule_params(
        &self,
        idx2: usize,
        idx3: usize,
        mmff_prop: &MmffPropCollection,
    ) -> Result<MmffTor, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFTorsionEmpiricalRuleParams (AtomTyper.cpp:2874-3067)
        // RDKit❗✔️: // empirical rule to compute torsional parameters if
        // RDKit❗✔️: // tabulated parameters could not be found
        // RDKit❗✔️: // the indexes of the two central atoms J and K
        // RDKit❗✔️: // idx2 and idx3 must be supplied. The returned pointer
        // RDKit❗✔️: // to a MMFFTor object must be freed by the caller
        // RDKit❗✔️: const ForceFields::MMFF::MMFFTor *
        // RDKit❗✔️: MMFFMolProperties::getMMFFTorsionEmpiricalRuleParams(const ROMol &mol,
        // RDKit❗✔️:                                                      unsigned int idx2,
        // RDKit❗✔️:                                                      unsigned int idx3) {
        // RDKit❗✔️:   PRECONDITION(this->isValid(), "missing atom types - invalid force-field");
        debug_assert!(self.is_valid());

        // RDKit❗✔️:
        // RDKit❗✔️:   const MMFFPropCollection *mmffProp = DefaultParameters::getMMFFProp();
        // RDKit❗✔️:   const MMFFAromCollection *mmffArom = DefaultParameters::getMMFFArom();
        // RDKit❗✔️:   auto *mmffTorParams = new ForceFields::MMFF::MMFFTor();
        let mut mmff_tor_params = MmffTor {
            v1: 0.0,
            v2: 0.0,
            v3: 0.0,
        };
        // RDKit❗✔️:   unsigned int jAtomType = this->getMMFFAtomType(idx2);
        let j_atom_type = self.get_mmff_atom_type(idx2)?;
        // RDKit❗✔️:   unsigned int kAtomType = this->getMMFFAtomType(idx3);
        let k_atom_type = self.get_mmff_atom_type(idx3)?;
        // RDKit❗✔️:   const MMFFProp *jMMFFProp = (*mmffProp)(jAtomType);
        let j_mmff_prop = mmff_prop.get(u32::from(j_atom_type)).ok_or(
            MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: j_atom_type,
            },
        )?;
        // RDKit❗✔️:   const MMFFProp *kMMFFProp = (*mmffProp)(kAtomType);
        let k_mmff_prop = mmff_prop.get(u32::from(k_atom_type)).ok_or(
            MmffMolPropertiesError::AtomTypePropertiesMissing {
                atom_type: k_atom_type,
            },
        )?;
        // RDKit❗✔️:   const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
        let bond = self
            .bond_between_atom_indices(idx2, idx3)
            .expect("torsion empirical rules require the validated central bond");
        // RDKit❗✔️:   double U[2] = {0.0, 0.0};
        let mut u = [0.0_f64; 2];
        // RDKit❗✔️:   double V[2] = {0.0, 0.0};
        let mut v = [0.0_f64; 2];
        // RDKit❗✔️:   double W[2] = {0.0, 0.0};
        let mut w = [0.0_f64; 2];
        // RDKit❗✔️:   double beta = 0.0;
        // RDKit❗✔️:   double pi_jk = 0.0;
        // RDKit❗✔️:   const auto N_jk = (double)((jMMFFProp->crd - 1) * (kMMFFProp->crd - 1));
        let n_jk = (f64::from(j_mmff_prop.crd) - 1.0) * (f64::from(k_mmff_prop.crd) - 1.0);
        // RDKit❗✔️:   int atomicNum[2] = {mol.getAtomWithIdx(idx2)->getAtomicNum(),
        // RDKit❗✔️:                       mol.getAtomWithIdx(idx3)->getAtomicNum()};
        let atomic_num = [
            self.topology.atoms.as_slice()[idx2].atomic_number(),
            self.topology.atoms.as_slice()[idx3].atomic_number(),
        ];

        // RDKit❗✔️:
        // RDKit❗✔️:   for (unsigned int i = 0; i < 2; ++i) {
        // RDKit❗✔️:     switch (atomicNum[i]) {
        for i in 0..2 {
            match atomic_num[i] {
                // RDKit❗✔️:       // carbon
                // RDKit❗✔️:       case 6:
                // RDKit❗✔️:         U[i] = 2.0;
                // RDKit❗✔️:         V[i] = 2.12;
                // RDKit❗✔️:         break;
                6 => {
                    u[i] = 2.0;
                    v[i] = 2.12;
                }

                // RDKit❗✔️:
                // RDKit❗✔️:       // nitrogen
                // RDKit❗✔️:       case 7:
                // RDKit❗✔️:         U[i] = 2.0;
                // RDKit❗✔️:         V[i] = 1.5;
                // RDKit❗✔️:         break;
                7 => {
                    u[i] = 2.0;
                    v[i] = 1.5;
                }

                // RDKit❗✔️:
                // RDKit❗✔️:       // oxygen
                // RDKit❗✔️:       case 8:
                // RDKit❗✔️:         U[i] = 2.0;
                // RDKit❗✔️:         V[i] = 0.2;
                // RDKit❗✔️:         W[i] = 2.0;
                // RDKit❗✔️:         break;
                8 => {
                    u[i] = 2.0;
                    v[i] = 0.2;
                    w[i] = 2.0;
                }

                // RDKit❗✔️:
                // RDKit❗✔️:       // silicon
                // RDKit❗✔️:       case 14:
                // RDKit❗✔️:         U[i] = 1.25;
                // RDKit❗✔️:         V[i] = 1.22;
                // RDKit❗✔️:         break;
                14 => {
                    u[i] = 1.25;
                    v[i] = 1.22;
                }

                // RDKit❗✔️:
                // RDKit❗✔️:       // phosphorus
                // RDKit❗✔️:       case 15:
                // RDKit❗✔️:         U[i] = 1.25;
                // RDKit❗✔️:         V[i] = 2.40;
                // RDKit❗✔️:         break;
                15 => {
                    u[i] = 1.25;
                    v[i] = 2.40;
                }

                // RDKit❗✔️:
                // RDKit❗✔️:       // sulfur
                // RDKit❗✔️:       case 16:
                // RDKit❗✔️:         U[i] = 1.25;
                // RDKit❗✔️:         V[i] = 0.49;
                // RDKit❗✔️:         W[i] = 8.0;
                // RDKit❗✔️:         break;
                16 => {
                    u[i] = 1.25;
                    v[i] = 0.49;
                    w[i] = 8.0;
                }
                _ => {}
            }
            // RDKit❗✔️:     }
            // RDKit❗✔️:   }
        }

        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (a)
        // RDKit❗✔️:   if (jMMFFProp->linh || kMMFFProp->linh) {
        if j_mmff_prop.linh != 0 || k_mmff_prop.linh != 0 {
            // RDKit❗✔️:     mmffTorParams->V1 = 0.0;
            // RDKit❗✔️:     mmffTorParams->V2 = 0.0;
            // RDKit❗✔️:     mmffTorParams->V3 = 0.0;
            // RDKit❗✔️:   }
            // RDKit❗✔️:
            // RDKit❗✔️:   // rule (b)
            // RDKit❗✔️:   else if (mmffArom->isMMFFAromatic(jAtomType) &&
            // RDKit❗✔️:            mmffArom->isMMFFAromatic(kAtomType) && bond->getIsAromatic()) {
        } else if mmff_is_aromatic_atom_type(j_atom_type)
            && mmff_is_aromatic_atom_type(k_atom_type)
            && bond.is_aromatic()
        {
            // RDKit❗✔️:     beta = ((((jMMFFProp->val == 3) && (kMMFFProp->val == 4)) ||
            // RDKit❗✔️:              ((jMMFFProp->val == 4) && (kMMFFProp->val == 3)))
            // RDKit❗✔️:                 ? 3.0
            // RDKit❗✔️:                 : 6.0);
            let beta = if (j_mmff_prop.val == 3 && k_mmff_prop.val == 4)
                || (j_mmff_prop.val == 4 && k_mmff_prop.val == 3)
            {
                3.0
            } else {
                6.0
            };
            // RDKit❗✔️:     pi_jk = (((jMMFFProp->pilp == 0) && (kMMFFProp->pilp == 0)) ? 0.5 : 0.3);
            let pi_jk = if j_mmff_prop.pilp == 0 && k_mmff_prop.pilp == 0 {
                0.5
            } else {
                0.3
            };
            // RDKit❗✔️:     mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
            mmff_tor_params.v2 = beta * pi_jk * (u[0] * u[1]).sqrt();
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (c)
        // RDKit❗✔️:   else if (bond->getBondType() == Bond::DOUBLE) {
        } else if bond.order() == BondOrder::Double {
            // RDKit❗✔️:     beta = 6.0;
            let beta = 6.0;
            // RDKit❗✔️:     pi_jk = (((jMMFFProp->mltb == 2) && (kMMFFProp->mltb == 2)) ? 1.0 : 0.4);
            let pi_jk = if j_mmff_prop.mltb == 2 && k_mmff_prop.mltb == 2 {
                1.0
            } else {
                0.4
            };
            // RDKit❗✔️:     mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
            mmff_tor_params.v2 = beta * pi_jk * (u[0] * u[1]).sqrt();
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (d)
        // RDKit❗✔️:   else if ((jMMFFProp->crd == 4) && (kMMFFProp->crd == 4)) {
        } else if j_mmff_prop.crd == 4 && k_mmff_prop.crd == 4 {
            // RDKit❗✔️:     mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
            mmff_tor_params.v3 = (v[0] * v[1]).sqrt() / n_jk;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (e)
        // RDKit❗✔️:   else if ((jMMFFProp->crd == 4) && (kMMFFProp->crd != 4)) {
        } else if j_mmff_prop.crd == 4 && k_mmff_prop.crd != 4 {
            // RDKit❗✔️:     if (((kMMFFProp->crd == 3) &&
            // RDKit❗✔️:          (((kMMFFProp->val == 4) || (kMMFFProp->val == 34)) ||
            // RDKit❗✔️:           kMMFFProp->mltb)) ||
            // RDKit❗✔️:         ((kMMFFProp->crd == 2) && ((kMMFFProp->val == 3) || kMMFFProp->mltb))) {
            if (k_mmff_prop.crd == 3
                && (k_mmff_prop.val == 4 || k_mmff_prop.val == 34 || k_mmff_prop.mltb != 0))
                || (k_mmff_prop.crd == 2 && (k_mmff_prop.val == 3 || k_mmff_prop.mltb != 0))
            {
                // RDKit❗✔️:       mmffTorParams->V1 = 0.0;
                // RDKit❗✔️:       mmffTorParams->V2 = 0.0;
                // RDKit❗✔️:       mmffTorParams->V3 = 0.0;
                // RDKit❗✔️:     } else {
            } else {
                // RDKit❗✔️:       mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
                mmff_tor_params.v3 = (v[0] * v[1]).sqrt() / n_jk;
            }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (f)
        // RDKit❗✔️:   else if ((kMMFFProp->crd == 4) && (jMMFFProp->crd != 4)) {
        } else if k_mmff_prop.crd == 4 && j_mmff_prop.crd != 4 {
            // RDKit❗✔️:     if (((jMMFFProp->crd == 3) &&
            // RDKit❗✔️:          (((jMMFFProp->val == 4) || (jMMFFProp->val == 34)) ||
            // RDKit❗✔️:           jMMFFProp->mltb)) ||
            // RDKit❗✔️:         ((jMMFFProp->crd == 2) && ((jMMFFProp->val == 3) || jMMFFProp->mltb))) {
            if (j_mmff_prop.crd == 3
                && (j_mmff_prop.val == 4 || j_mmff_prop.val == 34 || j_mmff_prop.mltb != 0))
                || (j_mmff_prop.crd == 2 && (j_mmff_prop.val == 3 || j_mmff_prop.mltb != 0))
            {
                // RDKit❗✔️:       mmffTorParams->V1 = 0.0;
                // RDKit❗✔️:       mmffTorParams->V2 = 0.0;
                // RDKit❗✔️:       mmffTorParams->V3 = 0.0;
                // RDKit❗✔️:     } else {
            } else {
                // RDKit❗✔️:       mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
                mmff_tor_params.v3 = (v[0] * v[1]).sqrt() / n_jk;
            }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (g)
        // RDKit❗✔️:   else if (((bond->getBondType() == Bond::SINGLE) && jMMFFProp->mltb &&
        // RDKit❗✔️:             kMMFFProp->mltb) ||
        // RDKit❗✔️:            (jMMFFProp->mltb && kMMFFProp->pilp) ||
        // RDKit❗✔️:            (jMMFFProp->pilp && kMMFFProp->mltb)) {
        } else if (bond.order() == BondOrder::Single
            && j_mmff_prop.mltb != 0
            && k_mmff_prop.mltb != 0)
            || (j_mmff_prop.mltb != 0 && k_mmff_prop.pilp != 0)
            || (j_mmff_prop.pilp != 0 && k_mmff_prop.mltb != 0)
        {
            // RDKit❗✔️:     // case (1)
            // RDKit❗✔️:     if (jMMFFProp->pilp && kMMFFProp->pilp) {
            if j_mmff_prop.pilp != 0 && k_mmff_prop.pilp != 0 {
                // RDKit❗✔️:       mmffTorParams->V1 = 0.0;
                // RDKit❗✔️:       mmffTorParams->V2 = 0.0;
                // RDKit❗✔️:       mmffTorParams->V3 = 0.0;
                // RDKit❗✔️:     }
                // RDKit❗✔️:     // case (2)
                // RDKit❗✔️:     else if (jMMFFProp->pilp && kMMFFProp->mltb) {
            } else if j_mmff_prop.pilp != 0 && k_mmff_prop.mltb != 0 {
                // RDKit❗✔️:       beta = 6.0;
                let beta = 6.0;
                // RDKit❗✔️:       if (jMMFFProp->mltb == 1) {
                // RDKit❗✔️:         pi_jk = 0.5;
                let pi_jk = if j_mmff_prop.mltb == 1 {
                    0.5
                // RDKit❗✔️:       } else if ((getPeriodicTableRow(atomicNum[0]) == 2) &&
                // RDKit❗✔️:                  (getPeriodicTableRow(atomicNum[1]) == 2)) {
                // RDKit❗✔️:         pi_jk = 0.3;
                } else if mmff_periodic_table_row(atomic_num[0]) == 2
                    && mmff_periodic_table_row(atomic_num[1]) == 2
                {
                    0.3
                // RDKit❗✔️:       } else if ((getPeriodicTableRow(atomicNum[0]) != 2) ||
                // RDKit❗✔️:                  (getPeriodicTableRow(atomicNum[1]) != 2)) {
                // RDKit❗✔️:         pi_jk = 0.15;
                } else {
                    0.15
                };
                // RDKit❗✔️:       }
                // RDKit❗✔️:       mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
                mmff_tor_params.v2 = beta * pi_jk * (u[0] * u[1]).sqrt();
            // RDKit❗✔️:     }
            // RDKit❗✔️:     // case (3)
            // RDKit❗✔️:     else if (kMMFFProp->pilp && jMMFFProp->mltb) {
            } else if k_mmff_prop.pilp != 0 && j_mmff_prop.mltb != 0 {
                // RDKit❗✔️:       beta = 6.0;
                let beta = 6.0;
                // RDKit❗✔️:       if (kMMFFProp->mltb == 1) {
                // RDKit❗✔️:         pi_jk = 0.5;
                let pi_jk = if k_mmff_prop.mltb == 1 {
                    0.5
                // RDKit❗✔️:       } else if ((getPeriodicTableRow(atomicNum[0]) == 2) &&
                // RDKit❗✔️:                  (getPeriodicTableRow(atomicNum[1]) == 2)) {
                // RDKit❗✔️:         pi_jk = 0.3;
                } else if mmff_periodic_table_row(atomic_num[0]) == 2
                    && mmff_periodic_table_row(atomic_num[1]) == 2
                {
                    0.3
                // RDKit❗✔️:       } else if ((getPeriodicTableRow(atomicNum[0]) != 2) ||
                // RDKit❗✔️:                  (getPeriodicTableRow(atomicNum[1]) != 2)) {
                // RDKit❗✔️:         pi_jk = 0.15;
                } else {
                    0.15
                };
                // RDKit❗✔️:       }
                // RDKit❗✔️:       mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
                mmff_tor_params.v2 = beta * pi_jk * (u[0] * u[1]).sqrt();
            // RDKit❗✔️:     }
            // RDKit❗✔️:     // case (4)
            // RDKit❗✔️:     else if (((jMMFFProp->mltb == 1) || (kMMFFProp->mltb == 1)) &&
            // RDKit❗✔️:              ((atomicNum[0] != 6) || (atomicNum[1] != 6))) {
            } else if (j_mmff_prop.mltb == 1 || k_mmff_prop.mltb == 1)
                && (atomic_num[0] != 6 || atomic_num[1] != 6)
            {
                // RDKit❗✔️:       beta = 6.0;
                // RDKit❗✔️:       pi_jk = 0.4;
                // RDKit❗✔️:       mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
                mmff_tor_params.v2 = 6.0 * 0.4 * (u[0] * u[1]).sqrt();
            // RDKit❗✔️:     }
            // RDKit❗✔️:     // case (5)
            // RDKit❗✔️:     else {
            } else {
                // RDKit❗✔️:       beta = 6.0;
                // RDKit❗✔️:       pi_jk = 0.15;
                // RDKit❗✔️:       mmffTorParams->V2 = beta * pi_jk * sqrt(U[0] * U[1]);
                mmff_tor_params.v2 = 6.0 * 0.15 * (u[0] * u[1]).sqrt();
            }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // rule (h)
        // RDKit❗✔️:   else {
        } else {
            // RDKit❗✔️:     if (((atomicNum[0] == 8) || (atomicNum[0] == 16)) &&
            // RDKit❗✔️:         ((atomicNum[1] == 8) || (atomicNum[1] == 16))) {
            if (atomic_num[0] == 8 || atomic_num[0] == 16)
                && (atomic_num[1] == 8 || atomic_num[1] == 16)
            {
                // RDKit❗✔️:       mmffTorParams->V2 = -sqrt(W[0] * W[1]);
                mmff_tor_params.v2 = -(w[0] * w[1]).sqrt();
            // RDKit❗✔️:     } else {
            } else {
                // RDKit❗✔️:       mmffTorParams->V3 = sqrt(V[0] * V[1]) / N_jk;
                mmff_tor_params.v3 = (v[0] * v[1]).sqrt() / n_jk;
            }
            // RDKit❗✔️:     }
            // RDKit❗✔️:   }
        }

        // RDKit❗✔️:
        // RDKit❗✔️:   return (const MMFFTor *)mmffTorParams;
        // RDKit❗✔️: }
        Ok(mmff_tor_params)
    }

    pub fn get_mmff_torsion_params(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
    ) -> Result<Option<(u32, MmffTor)>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFTorsionParams (AtomTyper.cpp:3629-3670)
        // RDKit✔️✔️: bool MMFFMolProperties::getMMFFTorsionParams(
        // RDKit✔️✔️:     const ROMol &mol, const unsigned int idx1, const unsigned int idx2,
        // RDKit✔️✔️:     const unsigned int idx3, const unsigned int idx4, unsigned int &torsionType,
        // RDKit✔️✔️:     MMFFTor &mmffTorsionParams) {
        // RDKit✔️✔️:   bool res = false;
        // RDKit✔️✔️:   if (isValid() && mol.getBondBetweenAtoms(idx1, idx2) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx2, idx3) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx3, idx4)) {
        if !self.is_valid() {
            return Ok(None);
        }
        if self.bond_between_atom_indices(idx1, idx2).is_none()
            || self.bond_between_atom_indices(idx2, idx3).is_none()
            || self.bond_between_atom_indices(idx3, idx4).is_none()
        {
            return Ok(None);
        }

        // RDKit✔️✔️:     unsigned int i;
        // RDKit✔️✔️:     unsigned int idx[4] = {idx1, idx2, idx3, idx4};
        let idx = [idx1, idx2, idx3, idx4];
        // RDKit✔️✔️:     unsigned int atomType[4];
        // RDKit✔️✔️:     const MMFFTorCollection *mmffTor =
        // RDKit✔️✔️:         DefaultParameters::getMMFFTor(getMMFFVariant() == "MMFF94s");
        let mmff_tor = default_mmff_tor(self.variant.is_mmff94s())
            .map_err(|e| MmffParamError::from(e.clone()))?;
        let mmff_def = default_mmff_def().map_err(|e| MmffParamError::from(e.clone()))?;
        let mmff_prop = default_mmff_prop().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     for (i = 0; i < 4; ++i) {
        // RDKit✔️✔️:       atomType[i] = getMMFFAtomType(idx[i]);
        // RDKit✔️✔️:     }
        let atom_type = [
            u32::from(self.get_mmff_atom_type(idx[0])?),
            u32::from(self.get_mmff_atom_type(idx[1])?),
            u32::from(self.get_mmff_atom_type(idx[2])?),
            u32::from(self.get_mmff_atom_type(idx[3])?),
        ];
        // RDKit✔️✔️:     const std::pair<unsigned int, unsigned int> torTypePair =
        // RDKit✔️✔️:         getMMFFTorsionType(mol, idx1, idx2, idx3, idx4);
        let tor_type_pair = self.get_mmff_torsion_type(idx1, idx2, idx3, idx4, mmff_prop)?;
        // RDKit✔️✔️:     bool areMMFFTorParamsEmpirical = false;
        // RDKit✔️✔️:     const std::pair<const unsigned int, const MMFFTor *> mmffTorPair =
        // RDKit✔️✔️:         mmffTor->getMMFFTorParams(DefaultParameters::getMMFFDef(), torTypePair,
        // RDKit✔️✔️:                                   atomType[0], atomType[1], atomType[2],
        // RDKit✔️✔️:                                   atomType[3]);
        let mmff_tor_pair = mmff_tor
            .get(
                mmff_def,
                tor_type_pair,
                atom_type[0],
                atom_type[1],
                atom_type[2],
                atom_type[3],
            )
            .map_err(MmffParamError::from)?;
        // RDKit✔️✔️:     torsionType = (mmffTorPair.first ? mmffTorPair.first : torTypePair.first);
        let mut torsion_type = if mmff_tor_pair.0 != 0 {
            mmff_tor_pair.0
        } else {
            tor_type_pair.0
        };
        // RDKit✔️✔️:     const MMFFTor *mmffTorParams = mmffTorPair.second;
        let mmff_tor_params = if let Some(mmff_tor_params) = mmff_tor_pair.1 {
            *mmff_tor_params
        } else {
            // RDKit❗✔️:     if (!mmffTorParams) {
            // RDKit❗✔️:       torsionType = torTypePair.first;
            torsion_type = tor_type_pair.0;
            // RDKit❗✔️:       mmffTorParams = getMMFFTorsionEmpiricalRuleParams(mol, idx2, idx3);
            let mmff_tor_params =
                self.get_mmff_torsion_empirical_rule_params(idx2, idx3, mmff_prop)?;
            // RDKit❗✔️:       areMMFFTorParamsEmpirical = true;
            // RDKit❗✔️:     }
            mmff_tor_params
        };
        // RDKit✔️✔️:     res =
        // RDKit✔️✔️:         (!(isDoubleZero(mmffTorParams->V1) && isDoubleZero(mmffTorParams->V2) &&
        // RDKit✔️✔️:            isDoubleZero(mmffTorParams->V3)));
        if is_double_zero(mmff_tor_params.v1)
            && is_double_zero(mmff_tor_params.v2)
            && is_double_zero(mmff_tor_params.v3)
        {
            return Ok(None);
        }
        // RDKit✔️✔️:     if (res) {
        // RDKit❗✔️:       mmffTorsionParams = *mmffTorParams;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (areMMFFTorParamsEmpirical) {
        // RDKit✔️✔️:       delete mmffTorParams;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Ok(Some((torsion_type, mmff_tor_params)))
    }

    pub fn get_mmff_oop_bend_params(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
    ) -> Result<Option<MmffOop>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFOopBendParams (AtomTyper.cpp:3672-3703)
        // RDKit✔️✔️: bool MMFFMolProperties::getMMFFOopBendParams(const ROMol &mol,
        // RDKit✔️✔️:                                              const unsigned int idx1,
        // RDKit✔️✔️:                                              const unsigned int idx2,
        // RDKit✔️✔️:                                              const unsigned int idx3,
        // RDKit✔️✔️:                                              const unsigned int idx4,
        // RDKit✔️✔️:                                              MMFFOop &mmffOopBendParams) {
        // RDKit✔️✔️:   bool res = false;
        // RDKit✔️✔️:   if (isValid() && mol.getBondBetweenAtoms(idx1, idx2) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx2, idx3) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx2, idx4)) {
        if !self.is_valid() {
            return Ok(None);
        }
        if self.bond_between_atom_indices(idx1, idx2).is_none()
            || self.bond_between_atom_indices(idx2, idx3).is_none()
            || self.bond_between_atom_indices(idx2, idx4).is_none()
        {
            return Ok(None);
        }

        // RDKit✔️✔️:     unsigned int i;
        // RDKit✔️✔️:     unsigned int idx[4] = {idx1, idx2, idx3, idx4};
        let idx = [idx1, idx2, idx3, idx4];
        // RDKit✔️✔️:     unsigned int atomType[4];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     const MMFFOopCollection *mmffOop =
        // RDKit✔️✔️:         DefaultParameters::getMMFFOop(getMMFFVariant() == "MMFF94s");
        let mmff_oop = default_mmff_oop(self.variant.is_mmff94s())
            .map_err(|e| MmffParamError::from(e.clone()))?;
        let mmff_def = default_mmff_def().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     for (i = 0; i < 4; ++i) {
        // RDKit✔️✔️:       atomType[i] = getMMFFAtomType(idx[i]);
        // RDKit✔️✔️:     }
        let atom_type = [
            u32::from(self.get_mmff_atom_type(idx[0])?),
            u32::from(self.get_mmff_atom_type(idx[1])?),
            u32::from(self.get_mmff_atom_type(idx[2])?),
            u32::from(self.get_mmff_atom_type(idx[3])?),
        ];
        // RDKit✔️✔️:     const MMFFOop *mmffOopParams =
        // RDKit✔️✔️:         (*mmffOop)(DefaultParameters::getMMFFDef(), atomType[0], atomType[1],
        // RDKit✔️✔️:                    atomType[2], atomType[3]);
        let mmff_oop_params = mmff_oop
            .get(
                mmff_def,
                atom_type[0],
                atom_type[1],
                atom_type[2],
                atom_type[3],
            )
            .map_err(MmffParamError::from)?;
        // RDKit✔️✔️:     // if no parameters could be found, we exclude this term (SURDOX02)
        // RDKit✔️✔️:     if (mmffOopParams) {
        let Some(mmff_oop_params) = mmff_oop_params else {
            return Ok(None);
        };
        // RDKit✔️✔️:       mmffOopBendParams = *mmffOopParams;
        // RDKit✔️✔️:       res = true;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Ok(Some(*mmff_oop_params))
    }

    pub fn get_mmff_vdw_params(
        &self,
        idx1: usize,
        idx2: usize,
    ) -> Result<Option<MmffVdwRijstarEps>, MmffMolPropertiesError> {
        // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::getMMFFVdWParams (AtomTyper.cpp:3703-3729)
        // RDKit✔️✔️: bool MMFFMolProperties::getMMFFVdWParams(const unsigned int idx1,
        // RDKit✔️✔️:                                          const unsigned int idx2,
        // RDKit✔️✔️:                                          MMFFVdWRijstarEps &mmffVdWParams) {
        // RDKit✔️✔️:   bool res = false;
        // RDKit✔️✔️:   if (isValid()) {
        if !self.is_valid() {
            return Ok(None);
        }

        // RDKit✔️✔️:     const MMFFVdWCollection *mmffVdW = DefaultParameters::getMMFFVdW();
        let mmff_vdw = default_mmff_vdw().map_err(|e| MmffParamError::from(e.clone()))?;
        // RDKit✔️✔️:     const unsigned int iAtomType = getMMFFAtomType(idx1);
        // RDKit✔️✔️:     const unsigned int jAtomType = getMMFFAtomType(idx2);
        let i_atom_type = u32::from(self.get_mmff_atom_type(idx1)?);
        let j_atom_type = u32::from(self.get_mmff_atom_type(idx2)?);
        // RDKit✔️✔️:     const MMFFVdW *mmffVdWParamsIAtom = (*mmffVdW)(iAtomType);
        // RDKit✔️✔️:     const MMFFVdW *mmffVdWParamsJAtom = (*mmffVdW)(jAtomType);
        let mmff_vdw_params_i_atom = mmff_vdw.get(i_atom_type);
        let mmff_vdw_params_j_atom = mmff_vdw.get(j_atom_type);
        // RDKit✔️✔️:     if (mmffVdWParamsIAtom && mmffVdWParamsJAtom) {
        let (Some(mmff_vdw_params_i_atom), Some(mmff_vdw_params_j_atom)) =
            (mmff_vdw_params_i_atom, mmff_vdw_params_j_atom)
        else {
            return Ok(None);
        };

        // RDKit✔️✔️:       mmffVdWParams.R_ij_starUnscaled = MMFF::Utils::calcUnscaledVdWMinimum(
        // RDKit✔️✔️:           mmffVdW, mmffVdWParamsIAtom, mmffVdWParamsJAtom);
        let r_ij_star_unscaled =
            calc_unscaled_vdw_minimum(mmff_vdw, mmff_vdw_params_i_atom, mmff_vdw_params_j_atom);
        // RDKit✔️✔️:       mmffVdWParams.epsilonUnscaled = MMFF::Utils::calcUnscaledVdWWellDepth(
        // RDKit✔️✔️:           mmffVdWParams.R_ij_starUnscaled, mmffVdWParamsIAtom,
        // RDKit✔️✔️:           mmffVdWParamsJAtom);
        let epsilon_unscaled = calc_unscaled_vdw_well_depth(
            r_ij_star_unscaled,
            mmff_vdw_params_i_atom,
            mmff_vdw_params_j_atom,
        );
        // RDKit✔️✔️:       mmffVdWParams.R_ij_star = mmffVdWParams.R_ij_starUnscaled;
        // RDKit✔️✔️:       mmffVdWParams.epsilon = mmffVdWParams.epsilonUnscaled;
        let mut params = MmffVdwRijstarEps {
            r_ij_star_unscaled,
            epsilon_unscaled,
            r_ij_star: r_ij_star_unscaled,
            epsilon: epsilon_unscaled,
        };
        // RDKit✔️✔️:       MMFF::Utils::scaleVdWParams(mmffVdWParams.R_ij_star,
        // RDKit✔️✔️:                                   mmffVdWParams.epsilon, mmffVdW,
        // RDKit✔️✔️:                                   mmffVdWParamsIAtom, mmffVdWParamsJAtom);
        scale_vdw_params(
            &mut params.r_ij_star,
            &mut params.epsilon,
            mmff_vdw,
            mmff_vdw_params_i_atom,
            mmff_vdw_params_j_atom,
        );
        // RDKit✔️✔️:       res = true;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Ok(Some(params))
    }

    fn angle_ring_size_3_or_4(&self, idx1: usize, idx2: usize, idx3: usize) -> Option<u32> {
        // BEGIN RDKIT CPP HELPER RDKit::MMFF::isAngleInRingOfSize3or4 (AtomTyper.cpp:355-397)
        // RDKit✔️✔️: // if the angle formed by atoms with indexes idx1, idx2, idx3
        // RDKit✔️✔️: // is in a ring of {3,4} atoms returns 3 or 4, respectively;
        // RDKit✔️✔️: // otherwise it returns 0
        // RDKit✔️✔️: unsigned int isAngleInRingOfSize3or4(const ROMol &mol, const unsigned int idx1,
        // RDKit✔️✔️:                                      const unsigned int idx2,
        // RDKit✔️✔️:                                      const unsigned int idx3) {
        // RDKit✔️✔️:   unsigned int ringSize = 0;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (mol.getBondBetweenAtoms(idx1, idx2) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx2, idx3)) {
        if self.bond_between_atom_indices(idx1, idx2).is_some()
            && self.bond_between_atom_indices(idx2, idx3).is_some()
        {
            // RDKit✔️✔️:     if (mol.getBondBetweenAtoms(idx3, idx1)) {
            // RDKit✔️✔️:       ringSize = 3;
            if self.bond_between_atom_indices(idx3, idx1).is_some() {
                return Some(3);
            }
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       std::set<unsigned int> s1;
            // RDKit✔️✔️:       std::set<unsigned int> s2;
            // RDKit✔️✔️:       std::vector<int> intersect;
            // RDKit✔️✔️:       ROMol::ADJ_ITER nbrIdx;
            // RDKit✔️✔️:       ROMol::ADJ_ITER endNbrs;
            // RDKit✔️✔️:       unsigned int newIdx;
            // RDKit✔️✔️:       boost::tie(nbrIdx, endNbrs) =
            // RDKit✔️✔️:           mol.getAtomNeighbors(mol.getAtomWithIdx(idx1));
            // RDKit✔️✔️:       for (; nbrIdx != endNbrs; ++nbrIdx) {
            // RDKit✔️✔️:         newIdx = mol[*nbrIdx]->getIdx();
            // RDKit✔️✔️:         if (newIdx != idx2) {
            // RDKit✔️✔️:           s1.insert(newIdx);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            let neighbors_1 = self.atom_neighbor_indices_except(idx1, idx2);
            // RDKit✔️✔️:       boost::tie(nbrIdx, endNbrs) =
            // RDKit✔️✔️:           mol.getAtomNeighbors(mol.getAtomWithIdx(idx3));
            // RDKit✔️✔️:       for (; nbrIdx != endNbrs; ++nbrIdx) {
            // RDKit✔️✔️:         newIdx = mol[*nbrIdx]->getIdx();
            // RDKit✔️✔️:         if (newIdx != idx2) {
            // RDKit✔️✔️:           s2.insert(newIdx);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            let neighbors_3 = self.atom_neighbor_indices_except(idx3, idx2);
            // RDKit✔️✔️:       std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
            // RDKit✔️✔️:                             std::back_inserter(intersect));
            // RDKit✔️✔️:       if (intersect.size()) {
            // RDKit✔️✔️:         ringSize = 4;
            if neighbors_1
                .iter()
                .any(|neighbor| neighbors_3.contains(neighbor))
            {
                return Some(4);
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return ringSize;
        // RDKit✔️✔️: }
        None
    }

    fn torsion_ring_size_4_or_5(
        &self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
    ) -> Option<u32> {
        // BEGIN RDKIT CPP HELPER RDKit::MMFF::isTorsionInRingOfSize4or5 (AtomTyper.cpp:403-445)
        // RDKit✔️✔️: unsigned int isTorsionInRingOfSize4or5(const ROMol &mol,
        // RDKit✔️✔️:                                        const unsigned int idx1,
        // RDKit✔️✔️:                                        const unsigned int idx2,
        // RDKit✔️✔️:                                        const unsigned int idx3,
        // RDKit✔️✔️:                                        const unsigned int idx4) {
        // RDKit✔️✔️:   unsigned int ringSize = 0;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (mol.getBondBetweenAtoms(idx1, idx2) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx2, idx3) &&
        // RDKit✔️✔️:       mol.getBondBetweenAtoms(idx3, idx4)) {
        if self.bond_between_atom_indices(idx1, idx2).is_some()
            && self.bond_between_atom_indices(idx2, idx3).is_some()
            && self.bond_between_atom_indices(idx3, idx4).is_some()
        {
            // RDKit✔️✔️:     if (mol.getBondBetweenAtoms(idx4, idx1)) {
            // RDKit✔️✔️:       ringSize = 4;
            if self.bond_between_atom_indices(idx4, idx1).is_some() {
                return Some(4);
            }
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       std::set<unsigned int> s1;
            // RDKit✔️✔️:       std::set<unsigned int> s2;
            // RDKit✔️✔️:       std::vector<int> intersect;
            // RDKit✔️✔️:       ROMol::ADJ_ITER nbrIdx;
            // RDKit✔️✔️:       ROMol::ADJ_ITER endNbrs;
            // RDKit✔️✔️:       unsigned int newIdx;
            // RDKit✔️✔️:       boost::tie(nbrIdx, endNbrs) =
            // RDKit✔️✔️:           mol.getAtomNeighbors(mol.getAtomWithIdx(idx1));
            // RDKit✔️✔️:       for (; nbrIdx != endNbrs; ++nbrIdx) {
            // RDKit✔️✔️:         newIdx = mol[*nbrIdx]->getIdx();
            // RDKit✔️✔️:         if (newIdx != idx2) {
            // RDKit✔️✔️:           s1.insert(newIdx);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            let neighbors_1 = self.atom_neighbor_indices_except(idx1, idx2);
            // RDKit✔️✔️:       boost::tie(nbrIdx, endNbrs) =
            // RDKit✔️✔️:           mol.getAtomNeighbors(mol.getAtomWithIdx(idx4));
            // RDKit✔️✔️:       for (; nbrIdx != endNbrs; ++nbrIdx) {
            // RDKit✔️✔️:         newIdx = mol[*nbrIdx]->getIdx();
            // RDKit✔️✔️:         if (newIdx != idx3) {
            // RDKit✔️✔️:           s2.insert(newIdx);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            let neighbors_4 = self.atom_neighbor_indices_except(idx4, idx3);
            // RDKit✔️✔️:       std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
            // RDKit✔️✔️:                             std::back_inserter(intersect));
            // RDKit✔️✔️:       if (intersect.size()) {
            // RDKit✔️✔️:         ringSize = 5;
            if neighbors_1
                .iter()
                .any(|neighbor| neighbors_4.contains(neighbor))
            {
                return Some(5);
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return ringSize;
        // RDKit✔️✔️: }
        None
    }

    fn atom_neighbor_indices_except(&self, atom_index: usize, excluded: usize) -> Vec<usize> {
        self.topology
            .bonds
            .as_slice()
            .iter()
            .filter_map(|bond| {
                if bond.begin().index() == atom_index {
                    Some(bond.end().index())
                } else if bond.end().index() == atom_index {
                    Some(bond.begin().index())
                } else {
                    None
                }
            })
            .filter(|neighbor| *neighbor != excluded)
            .collect()
    }

    fn bond_between_atom_indices(&self, idx1: usize, idx2: usize) -> Option<&Bond> {
        let atom_1 = AtomId::new(idx1);
        let atom_2 = AtomId::new(idx2);
        self.topology.bonds.as_slice().iter().find(|bond| {
            (bond.begin() == atom_1 && bond.end() == atom_2)
                || (bond.begin() == atom_2 && bond.end() == atom_1)
        })
    }
}

fn mmff_is_atom_in_aromatic_ring_of_size(
    mol: &TopologyBlock,
    ring_info: &RingInfo,
    atom_id: AtomId,
    ring_size: usize,
) -> bool {
    // BEGIN RDKIT CPP METHOD RDKit::MMFF::RingMembershipSize::isAtomInAromaticRingOfSize (AtomTyper.cpp:174-188)
    // RDKit✔️✔️: bool RingMembershipSize::isAtomInAromaticRingOfSize(
    // RDKit✔️✔️:     const Atom *atom, const unsigned int ringSize) const {
    // RDKit✔️✔️:   auto it = d_ringSizeMembershipMap.find(ringSize);
    // RDKit✔️✔️:   bool isAromatic = (it != d_ringSizeMembershipMap.end());
    // RDKit✔️✔️:   if (isAromatic) {
    // RDKit✔️✔️:     auto it2 = it->second.find(atom->getIdx());
    // RDKit✔️✔️:     isAromatic = (it2 != it->second.end());
    // RDKit✔️✔️:     if (isAromatic) {
    // RDKit✔️✔️:       isAromatic = it2->second.getIsInAromaticRing();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return isAromatic;
    // RDKit✔️✔️: }
    // END RDKIT CPP METHOD RDKit::MMFF::RingMembershipSize::isAtomInAromaticRingOfSize
    ring_info
        .atom_rings()
        .iter()
        .filter(|ring| ring.len() == ring_size && ring.contains(&atom_id))
        .any(|ring| {
            ring.iter().enumerate().all(|(idx, begin)| {
                let end = ring[(idx + 1) % ring.len()];
                mol.adjacency
                    .neighbors_of(begin.index())
                    .iter()
                    .find(|neighbor| neighbor.atom_index == end.index())
                    .is_some_and(|neighbor| {
                        mol.bonds.as_slice()[neighbor.bond.index()].is_aromatic()
                    })
            })
        })
}

fn set_mmff_heavy_atom_type(
    mol: &TopologyBlock,
    ring_info: &RingInfo,
    explicit_valence: &[i32],
    implicit_hydrogens: &[i32],
    atom: &Atom,
) -> Result<MmffAtomProperties, MmffMolPropertiesError> {
    // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::setMMFFHeavyAtomType (AtomTyper.cpp:504-2080)
    // RDKit❌❌: void MMFFMolProperties::setMMFFHeavyAtomType(const RingMembershipSize &rmSize,
    // RDKit❌❌:                                              const Atom *atom) {
    // RDKit❌❌:   unsigned int atomType = 0;
    let mut atom_type = 0_u8;
    // RDKit❌❌:   unsigned int i;
    // RDKit❌❌:   unsigned int j;
    // RDKit❌❌:   unsigned int nTermObondedToN = 0;
    // RDKit❌❌:   bool alphaOrBetaInSameRing = false;
    // RDKit❌❌:   bool isAlphaOS = false;
    // RDKit❌❌:   bool isBetaOS = false;
    // RDKit❗✔️:   bool isNSO2orNSO3orNCN = false;
    let mut is_nso2_or_nso3_or_ncn = false;
    // RDKit✔️✔️:   const ROMol &mol = atom->getOwningMol();
    // RDKit❗✔️:   RingInfo *ringInfo = mol.getRingInfo();
    // RDKit❗✔️:   ROMol::ADJ_ITER nbrIdx;
    // RDKit❗✔️:   ROMol::ADJ_ITER endNbrs;
    // RDKit❗✔️:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit❗✔️:   ROMol::ADJ_ITER end2Nbrs;
    // RDKit❗✔️:   ROMol::ADJ_ITER nbr3Idx;
    // RDKit❗✔️:   ROMol::ADJ_ITER end3Nbrs;
    // RDKit❌❌:   std::vector<const Atom *> alphaHet;
    // RDKit❌❌:   std::vector<const Atom *> betaHet;
    //
    // RDKit❗✔️:   if (atom->getIsAromatic()) {
    if atom.is_aromatic() {
        // RDKit❗✔️:     if (rmSize.isAtomInAromaticRingOfSize(atom, 5)) {
        // RDKit❗✔️:       // 5-membered aromatic rings
        // RDKit❌❌:       // if ipso is carbon or nitrogen, find eventual alpha and beta heteroatoms
        // RDKit❌❌:       if ((atom->getAtomicNum() == 6) || (atom->getAtomicNum() == 7)) {
        // RDKit❌❌:         // loop over alpha neighbors
        // RDKit❌❌:         boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
        // RDKit❌❌:         for (; nbrIdx != endNbrs; ++nbrIdx) {
        // RDKit❌❌:           const Atom *nbrAtom = mol[*nbrIdx];
        // RDKit❌❌:           // if the alpha neighbor is not in a 5-membered aromatic
        // RDKit❌❌:           // ring, skip to the next neighbor
        // RDKit❌❌:           if (!rmSize.isAtomInAromaticRingOfSize(nbrAtom, 5)) {
        // RDKit❌❌:             continue;
        // RDKit❌❌:           }
        // RDKit❌❌:           // if the alpha neighbor belongs to the same ring of ipso atom
        // RDKit❌❌:           // and it is either oxygen, sulfur, or non-N-oxide trivalent nitrogen,
        // RDKit❌❌:           // add it to the alpha atom vector
        // RDKit❌❌:           if (rmSize.areAtomsInSameRingOfSize(5, 2, atom, nbrAtom) &&
        // RDKit❌❌:               ((nbrAtom->getAtomicNum() == 8) ||
        // RDKit❌❌:                (nbrAtom->getAtomicNum() == 16) ||
        // RDKit❌❌:                ((nbrAtom->getAtomicNum() == 7) &&
        // RDKit❌❌:                 (nbrAtom->getTotalDegree() == 3) &&
        // RDKit❌❌:                 (!isAtomNOxide(nbrAtom))))) {
        // RDKit❌❌:             alphaHet.push_back(nbrAtom);
        // RDKit❌❌:           }
        // RDKit❌❌:           // loop over beta neighbors
        // RDKit❌❌:           boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
        // RDKit❌❌:           for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
        // RDKit❌❌:             const Atom *nbr2Atom = mol[*nbr2Idx];
        // RDKit❌❌:             // if we have gone back to the ipso atom, move on
        // RDKit❌❌:             if (nbr2Atom->getIdx() == atom->getIdx()) {
        // RDKit❌❌:               continue;
        // RDKit❌❌:             }
        // RDKit❌❌:             // if the beta neighbor is not in a 5-membered aromatic
        // RDKit❌❌:             // ring, skip to the next neighbor
        // RDKit❌❌:             if (!rmSize.isAtomInAromaticRingOfSize(nbr2Atom, 5)) {
        // RDKit❌❌:               continue;
        // RDKit❌❌:             }
        // RDKit❌❌:             // if the beta neighbor belongs to the same ring of ipso atom
        // RDKit❌❌:             // and it is either oxygen, sulfur, or non-N-oxide trivalent
        // RDKit❌❌:             // nitrogen,
        // RDKit❌❌:             // add it to the beta atom vector
        // RDKit❌❌:             if (rmSize.areAtomsInSameRingOfSize(5, 2, atom, nbr2Atom) &&
        // RDKit❌❌:                 ((nbr2Atom->getAtomicNum() == 8) ||
        // RDKit❌❌:                  (nbr2Atom->getAtomicNum() == 16) ||
        // RDKit❌❌:                  ((nbr2Atom->getAtomicNum() == 7) &&
        // RDKit❌❌:                   (nbr2Atom->getTotalDegree() == 3) &&
        // RDKit❌❌:                   (!isAtomNOxide(nbr2Atom))))) {
        // RDKit❌❌:               betaHet.push_back(nbr2Atom);
        // RDKit❌❌:             }
        // RDKit❌❌:           }
        // RDKit❌❌:         }
        // RDKit❌❌:       }
        // RDKit❗✔️:       switch (atom->getAtomicNum()) {
        // RDKit❌❌:         // Carbon
        // RDKit❌❌:         case 6:
        // RDKit❌❌:           // CB
        // RDKit❌❌:           // Aromatic carbon, e.g., in benzene
        // RDKit❌❌:           atomType = 37;
        // RDKit❌❌:           break;
        // RDKit❌❌:         // Nitrogen
        // RDKit❌❌:         case 7:
        // RDKit❌❌:           if (isAtomNOxide(atom)) {
        // RDKit❌❌:             // N5AX
        // RDKit❌❌:             // N-oxide nitrogen in 5-ring alpha position
        // RDKit❌❌:             // N5BX
        // RDKit❌❌:             // N-oxide nitrogen in 5-ring beta position
        // RDKit❌❌:             // N5OX
        // RDKit❌❌:             // N-oxide nitrogen in other 5-ring position
        // RDKit❌❌:             atomType = 82;
        // RDKit❌❌:             break;
        // RDKit❗✔️:           }
        // RDKit❗✔️:           // if there are neither alpha nor beta heteroatoms
        // RDKit❗✔️:           if ((!(alphaHet.size())) && (!(betaHet.size()))) {
        // RDKit❗✔️:             // if it is nitrogen
        // RDKit❗✔️:             // if valence is 3, it's pyrrole nitrogen
        // RDKit❗✔️:             if (atom->getTotalDegree() == 3) {
        // RDKit❗✔️:               // NPYL
        // RDKit❗✔️:               // Aromatic 5-ring nitrogen with pi lone pair
        // RDKit❗✔️:               atomType = 39;
        // RDKit❗✔️:               break;
        // RDKit❗✔️:             }
        // RDKit❌❌:             // otherwise it is anionic
        // RDKit❌❌:             // N5M
        // RDKit❌❌:             // Nitrogen in 5-ring aromatic anion
        // RDKit❌❌:             atomType = 76;
        // RDKit❌❌:             break;
        // RDKit❌❌:           }
        // RDKit❌❌:           // ... remaining 5-ring nitrogen heteroatom-neighborhood branches ...
        // RDKit❌❌:           break;
        // RDKit❌❌:         // Oxygen
        // RDKit❌❌:         case 8:
        // RDKit❌❌:           // OFUR
        // RDKit❌❌:           // Aromatic 5-ring oxygen with pi lone pair
        // RDKit❌❌:           atomType = 59;
        // RDKit❌❌:           break;
        // RDKit❌❌:         // Sulfur
        // RDKit❌❌:         case 16:
        // RDKit❌❌:           // STHI
        // RDKit❌❌:           // Aromatic 5-ring sulfur with pi lone pair
        // RDKit❌❌:           atomType = 44;
        // RDKit❌❌:           break;
        // RDKit❌❌:       }
        // RDKit❗✔️:     }
        if mmff_is_atom_in_aromatic_ring_of_size(mol, ring_info, atom.id(), 5) {
            let mut alpha_het = Vec::<AtomId>::new();
            let mut beta_het = Vec::<AtomId>::new();
            if matches!(atom.atomic_number(), 6 | 7) {
                for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                    let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                    if !mmff_is_atom_in_aromatic_ring_of_size(mol, ring_info, nbr_atom.id(), 5) {
                        continue;
                    }
                    if ring_info.are_atoms_in_same_ring_of_size(atom.id(), nbr_atom.id(), 5)
                        && (matches!(nbr_atom.atomic_number(), 8 | 16)
                            || (nbr_atom.atomic_number() == 7
                                && mmff_total_degree(
                                    mol,
                                    implicit_hydrogens,
                                    neighbor.atom_index,
                                )? == 3
                                && !is_mmff_atom_n_oxide(
                                    mol,
                                    implicit_hydrogens,
                                    neighbor.atom_index,
                                )?))
                    {
                        alpha_het.push(nbr_atom.id());
                    }
                    for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                        if neighbor2.atom_index == atom.id().index() {
                            continue;
                        }
                        let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                        if !mmff_is_atom_in_aromatic_ring_of_size(mol, ring_info, nbr2_atom.id(), 5)
                        {
                            continue;
                        }
                        if ring_info.are_atoms_in_same_ring_of_size(atom.id(), nbr2_atom.id(), 5)
                            && (matches!(nbr2_atom.atomic_number(), 8 | 16)
                                || (nbr2_atom.atomic_number() == 7
                                    && mmff_total_degree(
                                        mol,
                                        implicit_hydrogens,
                                        neighbor2.atom_index,
                                    )? == 3
                                    && !is_mmff_atom_n_oxide(
                                        mol,
                                        implicit_hydrogens,
                                        neighbor2.atom_index,
                                    )?))
                        {
                            beta_het.push(nbr2_atom.id());
                        }
                    }
                }
            }
            let is_alpha_os = alpha_het.iter().any(|atom_id| {
                matches!(
                    mol.atoms.as_slice()[atom_id.index()].atomic_number(),
                    8 | 16
                )
            });
            let is_beta_os = beta_het.iter().any(|atom_id| {
                matches!(
                    mol.atoms.as_slice()[atom_id.index()].atomic_number(),
                    8 | 16
                )
            });
            let alpha_or_beta_in_same_ring = alpha_het.iter().any(|alpha| {
                beta_het
                    .iter()
                    .any(|beta| ring_info.are_atoms_in_same_ring_of_size(*alpha, *beta, 5))
            });
            match atom.atomic_number() {
                6 => {
                    if beta_het.is_empty() {
                        let mut n_n = 0_u32;
                        let mut n_formal_charge = 0_u32;
                        let mut n_in_aromatic_5_ring = 0_u32;
                        let mut n_in_aromatic_6_ring = 0_u32;
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            if nbr_atom.atomic_number() == 7
                                && mmff_total_degree(mol, implicit_hydrogens, neighbor.atom_index)?
                                    == 3
                            {
                                n_n += 1;
                                if nbr_atom.formal_charge() > 0
                                    && !is_mmff_atom_n_oxide(
                                        mol,
                                        implicit_hydrogens,
                                        neighbor.atom_index,
                                    )?
                                {
                                    n_formal_charge += 1;
                                }
                                if mmff_is_atom_in_aromatic_ring_of_size(
                                    mol,
                                    ring_info,
                                    nbr_atom.id(),
                                    5,
                                ) {
                                    n_in_aromatic_5_ring += 1;
                                }
                                if mmff_is_atom_in_aromatic_ring_of_size(
                                    mol,
                                    ring_info,
                                    nbr_atom.id(),
                                    6,
                                ) {
                                    n_in_aromatic_6_ring += 1;
                                }
                            }
                        }
                        if (((n_n == 2) && n_in_aromatic_5_ring != 0)
                            || ((n_n == 3) && (n_in_aromatic_5_ring == 2)))
                            && n_formal_charge != 0
                            && n_in_aromatic_6_ring == 0
                        {
                            atom_type = 80;
                        }
                    }
                    if atom_type == 0 && alpha_het.len() == beta_het.len() {
                        let mut surrounded_by_benzene_c = true;
                        let mut surrounded_by_arom = true;
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            if nbr_atom.atomic_number() != 6
                                || !ring_info.is_atom_in_ring_of_size(nbr_atom.id(), 6)
                            {
                                surrounded_by_benzene_c = false;
                            }
                            if ring_info.are_atoms_in_same_ring_of_size(atom.id(), nbr_atom.id(), 5)
                                && !nbr_atom.is_aromatic()
                            {
                                surrounded_by_arom = false;
                            }
                        }
                        if (alpha_het.is_empty()
                            && beta_het.is_empty()
                            && !surrounded_by_benzene_c
                            && surrounded_by_arom)
                            || (!alpha_het.is_empty()
                                && !beta_het.is_empty()
                                && (!alpha_or_beta_in_same_ring || (!is_alpha_os && !is_beta_os)))
                        {
                            atom_type = 78;
                        }
                    }
                    if atom_type == 0
                        && !alpha_het.is_empty()
                        && (beta_het.is_empty() || is_alpha_os)
                    {
                        atom_type = 63;
                    }
                    if atom_type == 0
                        && !beta_het.is_empty()
                        && (alpha_het.is_empty() || is_beta_os)
                    {
                        atom_type = 64;
                    }
                }
                7 => {
                    let total_degree =
                        mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                    if is_mmff_atom_n_oxide(mol, implicit_hydrogens, atom.id().index())? {
                        atom_type = 82;
                    } else if alpha_het.is_empty() && beta_het.is_empty() {
                        if total_degree == 3 {
                            atom_type = 39;
                        } else {
                            atom_type = 76;
                        }
                    } else if total_degree == 3 && alpha_het.len() != beta_het.len() {
                        atom_type = 81;
                    } else if !alpha_het.is_empty() && (beta_het.is_empty() || is_alpha_os) {
                        atom_type = 65;
                    } else if !beta_het.is_empty() && (alpha_het.is_empty() || is_beta_os) {
                        atom_type = 66;
                    } else if !alpha_het.is_empty() && !beta_het.is_empty() {
                        atom_type = 79;
                    }
                }
                8 => {
                    // RDKit❗✔️:         // Oxygen
                    // RDKit❗✔️:         case 8:
                    // RDKit❗✔️:           // OFUR
                    // RDKit❗✔️:           // Aromatic 5-ring oxygen with pi lone pair
                    // RDKit❗✔️:           atomType = 59;
                    atom_type = 59;
                    // RDKit❗✔️:           break;
                }
                16 => {
                    // RDKit❗✔️:         // Sulfur
                    // RDKit❗✔️:         case 16:
                    // RDKit❗✔️:           // STHI
                    // RDKit❗✔️:           // Aromatic 5-ring sulfur with pi lone pair
                    // RDKit❗✔️:           atomType = 44;
                    atom_type = 44;
                    // RDKit❗✔️:           break;
                }
                // The C++ switch has no default: unlisted elements retain atomType == 0
                // and continue through the element-specific aliphatic branch below.
                _ => {}
            }
        }

        // RDKit✔️✔️:     if (!atomType && (rmSize.isAtomInAromaticRingOfSize(atom, 6))) {
        if atom_type == 0 && mmff_is_atom_in_aromatic_ring_of_size(mol, ring_info, atom.id(), 6) {
            // RDKit✔️✔️:       // 6-membered aromatic rings
            // RDKit✔️✔️:       switch (atom->getAtomicNum()) {
            match atom.atomic_number() {
                // RDKit✔️✔️:         // Carbon
                // RDKit✔️✔️:         case 6:
                // RDKit✔️✔️:           // CB
                // RDKit✔️✔️:           // Aromatic carbon, e.g., in benzene
                // RDKit✔️✔️:           atomType = 37;
                6 => atom_type = 37,
                // RDKit❗✔️:         // Nitrogen
                // RDKit❗✔️:         case 7:
                7 => {
                    // RDKit✔️✔️:           if (isAtomNOxide(atom)) {
                    if is_mmff_atom_n_oxide(mol, implicit_hydrogens, atom.id().index())? {
                        // RDKit✔️✔️:             // NPOX
                        // RDKit✔️✔️:             // Pyridinium N-oxide nitrogen
                        // RDKit✔️✔️:             atomType = 69;
                        atom_type = 69;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    let total_degree =
                        mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                    if atom_type == 0 && total_degree == 3 {
                        // RDKit✔️✔️:           if (atom->getTotalDegree() == 3) {
                        // RDKit✔️✔️:             // NPD+
                        // RDKit✔️✔️:             // Aromatic nitrogen in pyridinium
                        // RDKit✔️✔️:             atomType = 58;
                        atom_type = 58;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    if atom_type == 0 {
                        // RDKit✔️✔️:           // NPYD
                        // RDKit✔️✔️:           // Aromatic nitrogen with sigma lone pair
                        // RDKit✔️✔️:           atomType = 38;
                        atom_type = 38;
                        // RDKit✔️✔️:           break;
                    }
                }
                // The C++ switch has no default: unlisted elements retain atomType == 0
                // and continue through the element-specific aliphatic branch below.
                _ => {}
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
        }
        // RDKit❗✔️:   }
    }

    // RDKit❗✔️:   if (!atomType) {
    if atom_type == 0 {
        // RDKit❗✔️:     // Aliphatic heavy atom types
        //
        // RDKit❗✔️:     switch (atom->getAtomicNum()) {
        match atom.atomic_number() {
            // RDKit✔️✔️:       // Lithium
            // RDKit✔️✔️:       case 3:
            3 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           // LI+
                    // RDKit✔️✔️:           // Lithium cation
                    // RDKit✔️✔️:           atomType = 92;
                    atom_type = 92;
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit❗✔️:       // Carbon
            // RDKit❗✔️:       case 6:
            6 => {
                let total_degree = mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                // RDKit❗✔️:         // 4 neighbors
                // RDKit❗✔️:         if (atom->getTotalDegree() == 4) {
                if total_degree == 4 {
                    // RDKit❗✔️:           if (ringInfo->isAtomInRingOfSize(atom->getIdx(), 3)) {
                    if ring_info.is_atom_in_ring_of_size(atom.id(), 3) {
                        // RDKit❗✔️:             // CR3R
                        // RDKit❗✔️:             // Aliphatic carbon in 3-membered ring
                        // RDKit❗✔️:             atomType = 22;
                        atom_type = 22;
                    // RDKit❗✔️:             break;
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           if (ringInfo->isAtomInRingOfSize(atom->getIdx(), 4)) {
                    } else if ring_info.is_atom_in_ring_of_size(atom.id(), 4) {
                        // RDKit❗✔️:             // CR4R
                        // RDKit❗✔️:             // Aliphatic carbon in 4-membered ring
                        // RDKit❗✔️:             atomType = 20;
                        atom_type = 20;
                    // RDKit❗✔️:             break;
                    // RDKit❗✔️:           }
                    } else {
                        // RDKit❗✔️:           // CR
                        // RDKit❗✔️:           // Alkyl carbon
                        // RDKit❗✔️:           atomType = 1;
                        atom_type = 1;
                    }
                    // RDKit❗✔️:           break;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         // 3 neighbors
                    // RDKit❗✔️:         if (atom->getTotalDegree() == 3) {
                } else if total_degree == 3 {
                    // RDKit❗✔️:           unsigned int nN2 = 0;
                    let mut n_n2 = 0_u32;
                    // RDKit❗✔️:           unsigned int nN3 = 0;
                    let mut n_n3 = 0_u32;
                    // RDKit❗✔️:           unsigned int nO = 0;
                    let mut n_o = 0_u32;
                    // RDKit❗✔️:           unsigned int nS = 0;
                    let mut n_s = 0_u32;
                    // RDKit❗✔️:           unsigned int doubleBondedElement = 0;
                    let mut double_bonded_element = 0_u8;
                    // RDKit❗✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit❗✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        // RDKit❗✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                        // RDKit❗✔️:             // find if there is a double-bonded element
                        // RDKit❗✔️:             if ((mol.getBondBetweenAtoms(nbrAtom->getIdx(), atom->getIdx()))
                        // RDKit❗✔️:                     ->getBondType() == Bond::DOUBLE) {
                        if nbr_bond.order() == BondOrder::Double {
                            // RDKit❗✔️:               doubleBondedElement = nbrAtom->getAtomicNum();
                            double_bonded_element = nbr_atom.atomic_number();
                        }
                        // RDKit❗✔️:             }
                        let nbr_total_degree =
                            mmff_total_degree(mol, implicit_hydrogens, neighbor.atom_index)?;
                        // RDKit❗✔️:             // count how many terminal oxygen/sulfur atoms
                        // RDKit❗✔️:             // are bonded to ipso
                        // RDKit❗✔️:             if (nbrAtom->getTotalDegree() == 1) {
                        if nbr_total_degree == 1 {
                            // RDKit❗✔️:               if (nbrAtom->getAtomicNum() == 8) {
                            if nbr_atom.atomic_number() == 8 {
                                // RDKit❗✔️:                 ++nO;
                                n_o += 1;
                            // RDKit❗✔️:               } else if (nbrAtom->getAtomicNum() == 16) {
                            } else if nbr_atom.atomic_number() == 16 {
                                // RDKit❗✔️:                 ++nS;
                                n_s += 1;
                            }
                        // RDKit❗✔️:               }
                        // RDKit❗✔️:             } else if (nbrAtom->getAtomicNum() == 7) {
                        } else if nbr_atom.atomic_number() == 7 {
                            // RDKit❗✔️:               // count how many nitrogens with 3 neighbors
                            // RDKit❗✔️:               // are bonded to ipso
                            // RDKit❗✔️:               if (nbrAtom->getTotalDegree() == 3) {
                            if nbr_total_degree == 3 {
                                // RDKit❗✔️:                 ++nN3;
                                n_n3 += 1;
                            // RDKit❗✔️:               }
                            // RDKit❗✔️:               // count how many nitrogens with 2 neighbors
                            // RDKit❗✔️:               // are double-bonded to ipso
                            // RDKit❗✔️:               else if ((nbrAtom->getTotalDegree() == 2) &&
                            // RDKit❗✔️:                        ((mol.getBondBetweenAtoms(nbrAtom->getIdx(),
                            // RDKit❗✔️:                                                  atom->getIdx()))
                            // RDKit❗✔️:                             ->getBondType() == Bond::DOUBLE)) {
                            } else if nbr_total_degree == 2 && nbr_bond.order() == BondOrder::Double
                            {
                                // RDKit❗✔️:                 ++nN2;
                                n_n2 += 1;
                            }
                            // RDKit❗✔️:               }
                        }
                        // RDKit❗✔️:             }
                    }
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           // if there are two or more nitrogens with 3 neighbors each,
                    // RDKit❗✔️:           // and there are no nitrogens with two neighbors only,
                    // RDKit❗✔️:           // and carbon is double-bonded to nitrogen
                    // RDKit❗✔️:           if ((nN3 >= 2) && (!nN2) && (doubleBondedElement == 7)) {
                    if n_n3 >= 2 && n_n2 == 0 && double_bonded_element == 7 {
                        // RDKit❗✔️:             // CNN+
                        // RDKit❗✔️:             // Carbon in +N=C-N: resonance structures
                        // RDKit❗✔️:             // CGD+
                        // RDKit❗✔️:             // Guanidinium carbon
                        // RDKit❗✔️:             atomType = 57;
                        atom_type = 57;
                    // RDKit❗✔️:             break;
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           // if there are two terminal oxygen/sulfur atoms
                    // RDKit❗✔️:           if ((nO == 2) || (nS == 2)) {
                    } else if n_o == 2 || n_s == 2 {
                        // RDKit❗✔️:             // CO2M
                        // RDKit❗✔️:             // Carbon in carboxylate anion
                        // RDKit❗✔️:             // CS2M
                        // RDKit❗✔️:             // Carbon in thiocarboxylate anion
                        // RDKit❗✔️:             atomType = 41;
                        atom_type = 41;
                    // RDKit❗✔️:             break;
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           // if this carbon is in a 4-membered ring and
                    // RDKit❗✔️:           // is double-bonded to another carbon
                    // RDKit❗✔️:           if (ringInfo->isAtomInRingOfSize(atom->getIdx(), 4) &&
                    // RDKit❗✔️:               (doubleBondedElement == 6)) {
                    } else if ring_info.is_atom_in_ring_of_size(atom.id(), 4)
                        && double_bonded_element == 6
                    {
                        // RDKit❗✔️:             // CR4E
                        // RDKit❗✔️:             // Olefinic carbon in 4-membered ring
                        // RDKit❗✔️:             atomType = 30;
                        atom_type = 30;
                    // RDKit❗✔️:             break;
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           // if this carbon is is double-bonded to nitrogen,
                    // RDKit❗✔️:           // oxygen, phosphorus or sulfur
                    // RDKit❗✔️:           if ((doubleBondedElement == 7) || (doubleBondedElement == 8) ||
                    // RDKit❗✔️:               (doubleBondedElement == 15) || (doubleBondedElement == 16)) {
                    } else if matches!(double_bonded_element, 7 | 8 | 15 | 16) {
                        // RDKit❗✔️:             // C=N
                        // RDKit❗✔️:             // Imine-atomType carbon
                        // RDKit❗✔️:             // CGD
                        // RDKit❗✔️:             // Guanidine carbon
                        // RDKit❗✔️:             // C=O
                        // RDKit❗✔️:             // Generic carbonyl carbon
                        // RDKit❗✔️:             // C=OR
                        // RDKit❗✔️:             // Ketone or aldehyde carbonyl carbon
                        // RDKit❗✔️:             // C=ON
                        // RDKit❗✔️:             // Amide carbonyl carbon
                        // RDKit❗✔️:             // COO
                        // RDKit❗✔️:             // Carboxylic acid or ester carbonyl carbon
                        // RDKit❗✔️:             // COON
                        // RDKit❗✔️:             // Carbamate carbonyl carbon
                        // RDKit❗✔️:             // COOO
                        // RDKit❗✔️:             // Carbonic acid or ester carbonyl function
                        // RDKit❗✔️:             // C=OS
                        // RDKit❗✔️:             // Thioester carbonyl carbon, double bonded to O
                        // RDKit❗✔️:             // C=P
                        // RDKit❗✔️:             // Carbon doubly bonded to P
                        // RDKit❗✔️:             // C=S
                        // RDKit❗✔️:             // Thioester carbon, double bonded to S
                        // RDKit❗✔️:             // C=SN
                        // RDKit❗✔️:             // Thioamide carbon, double bonded to S
                        // RDKit❗✔️:             // CSO2
                        // RDKit❗✔️:             // Carbon in >C=SO2
                        // RDKit❗✔️:             // CS=O
                        // RDKit❗✔️:             // Sulfinyl carbon in >C=S=O
                        // RDKit❗✔️:             // CSS
                        // RDKit❗✔️:             // Thiocarboxylic acid or ester carbon
                        // RDKit❗✔️:             atomType = 3;
                        atom_type = 3;
                    // RDKit❗✔️:             break;
                    // RDKit❗✔️:           }
                    } else {
                        // RDKit✔️✔️:           // otherwise it must be generic sp2 carbon
                        // RDKit✔️✔️:           // C=C
                        // RDKit✔️✔️:           // Vinylic carbon
                        // RDKit✔️✔️:           // CSP2
                        // RDKit✔️✔️:           // Generic sp2 carbon
                        // RDKit✔️✔️:           atomType = 2;
                        atom_type = 2;
                    }
                    // RDKit✔️✔️:           break;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         // 2 neighbors
                    // RDKit❗✔️:         if (atom->getTotalDegree() == 2) {
                } else if total_degree == 2 {
                    // RDKit❗✔️:           // CSP
                    // RDKit❗✔️:           // Acetylenic carbon
                    // RDKit❗✔️:           // =C=
                    // RDKit❗✔️:           // Allenic carbon
                    // RDKit❗✔️:           atomType = 4;
                    atom_type = 4;
                    // RDKit❗✔️:           break;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         // 1 neighbor
                    // RDKit❗✔️:         if (atom->getTotalDegree() == 1) {
                } else if total_degree == 1 {
                    // RDKit❗✔️:           // C%-
                    // RDKit❗✔️:           // Isonitrile carbon
                    // RDKit❗✔️:           atomType = 60;
                    atom_type = 60;
                    // RDKit❗✔️:           break;
                }
            }
            // RDKit❗✔️:       // Nitrogen
            // RDKit❗✔️:       case 7:
            7 => {
                let total_degree = mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                let total_bond_order =
                    mmff_total_bond_order(explicit_valence, implicit_hydrogens, atom.id().index())?;
                // RDKit❗✔️:         // if the neighbor is phosphorus or sulfur
                // RDKit❗✔️:         // count the number of terminal oxygens bonded
                // RDKit❗✔️:         // to that phosphorus or sulfur atom
                let mut n_term_o_bonded_to_n = 0_u32;
                // RDKit❗✔️:         boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                // RDKit❗✔️:         for (; nbrIdx != endNbrs; ++nbrIdx) {
                for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                    // RDKit❗✔️:           const Atom *nbrAtom = mol[*nbrIdx];
                    let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                    let nbr_total_degree =
                        mmff_total_degree(mol, implicit_hydrogens, neighbor.atom_index)?;
                    // RDKit❗✔️:           // count how many terminal oxygen atoms
                    // RDKit❗✔️:           // are bonded to ipso
                    // RDKit❗✔️:           if ((nbrAtom->getAtomicNum() == 8) &&
                    // RDKit❗✔️:               (nbrAtom->getTotalDegree() == 1)) {
                    if nbr_atom.atomic_number() == 8 && nbr_total_degree == 1 {
                        // RDKit❗✔️:             ++nTermObondedToN;
                        n_term_o_bonded_to_n += 1;
                    }
                    // RDKit❗✔️:           }
                    // RDKit✔️✔️:           if (((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit✔️✔️:                 atom->getNumImplicitHs()) >= 3) &&
                    // RDKit✔️✔️:               ((nbrAtom->getAtomicNum() == 15) ||
                    // RDKit✔️✔️:                (nbrAtom->getAtomicNum() == 16))) {
                    if total_bond_order >= 3 && matches!(nbr_atom.atomic_number(), 15 | 16) {
                        // RDKit✔️✔️:             unsigned int nObondedToSP = 0;
                        let mut n_o_bonded_to_sp = 0_u32;
                        // RDKit✔️✔️:             boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                        // RDKit✔️✔️:             for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                        for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                            // RDKit✔️✔️:               const Atom *nbr2Atom = mol[*nbr2Idx];
                            let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                            // RDKit✔️✔️:               if ((nbr2Atom->getAtomicNum() == 8) &&
                            // RDKit✔️✔️:                   (nbr2Atom->getTotalDegree() == 1)) {
                            if nbr2_atom.atomic_number() == 8
                                && mmff_total_degree(mol, implicit_hydrogens, neighbor2.atom_index)?
                                    == 1
                            {
                                // RDKit✔️✔️:                 ++nObondedToSP;
                                n_o_bonded_to_sp += 1;
                            }
                            // RDKit✔️✔️:               }
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             // if there are two or more oxygens, ipso is a sulfonamide nitrogen
                        // RDKit✔️✔️:             if (!isNSO2orNSO3orNCN) {
                        if !is_nso2_or_nso3_or_ncn {
                            // RDKit✔️✔️:               isNSO2orNSO3orNCN = (nObondedToSP >= 2);
                            is_nso2_or_nso3_or_ncn = n_o_bonded_to_sp >= 2;
                        }
                        // RDKit✔️✔️:             }
                    }
                    // RDKit✔️✔️:           }
                }
                // RDKit❗✔️:         }
                // RDKit❗✔️:         // 4 neighbors
                // RDKit❗✔️:         if (atom->getTotalDegree() == 4) {
                if total_degree == 4 {
                    // RDKit✔️✔️:           if (isAtomNOxide(atom)) {
                    if is_mmff_atom_n_oxide(mol, implicit_hydrogens, atom.id().index())? {
                        // RDKit✔️✔️:             // N3OX
                        // RDKit✔️✔️:             // sp3-hybridized N-oxide nitrogen
                        // RDKit✔️✔️:             atomType = 68;
                        atom_type = 68;
                        // RDKit✔️✔️:             break;
                    } else {
                        // RDKit✔️✔️:           }
                        // RDKit✔️✔️:           // NR+
                        // RDKit✔️✔️:           // Quaternary nitrogen
                        // RDKit✔️✔️:           atomType = 34;
                        atom_type = 34;
                        // RDKit✔️✔️:           break;
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit❗✔️:         // 3 neighbors
                // RDKit❗✔️:         if (atom->getTotalDegree() == 3) {
                if atom_type == 0 && total_degree == 3 {
                    // RDKit❗✔️:           // total bond order >= 4
                    // RDKit❗✔️:           if ((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit❗✔️:                atom->getNumImplicitHs()) >= 4) {
                    if total_bond_order >= 4 {
                        // RDKit✔️✔️:             bool doubleBondedCN = false;
                        let mut double_bonded_cn = false;
                        // RDKit✔️✔️:             boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                        // RDKit✔️✔️:             for (; nbrIdx != endNbrs; ++nbrIdx) {
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            // RDKit✔️✔️:               const Atom *nbrAtom = mol[*nbrIdx];
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                            // RDKit✔️✔️:               // find if there is a double-bonded nitrogen,
                            // RDKit✔️✔️:               // or a carbon which is not bonded to other
                            // RDKit✔️✔️:               // nitrogen atoms with 3 neighbors
                            // RDKit✔️✔️:               if ((mol.getBondBetweenAtoms(nbrAtom->getIdx(), atom->getIdx()))
                            // RDKit✔️✔️:                       ->getBondType() == Bond::DOUBLE) {
                            if nbr_bond.order() == BondOrder::Double {
                                // RDKit✔️✔️:                 doubleBondedCN = ((nbrAtom->getAtomicNum() == 7) ||
                                // RDKit✔️✔️:                                   (nbrAtom->getAtomicNum() == 6));
                                double_bonded_cn =
                                    nbr_atom.atomic_number() == 7 || nbr_atom.atomic_number() == 6;
                                // RDKit✔️✔️:                 if (nbrAtom->getAtomicNum() == 6) {
                                if nbr_atom.atomic_number() == 6 {
                                    // RDKit✔️✔️:                   boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                                    // RDKit✔️✔️:                   for (; doubleBondedCN && (nbr2Idx != end2Nbrs); ++nbr2Idx) {
                                    for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index)
                                    {
                                        if !double_bonded_cn {
                                            break;
                                        }
                                        // RDKit✔️✔️:                     const Atom *nbr2Atom = mol[*nbr2Idx];
                                        let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                                        // RDKit✔️✔️:                     if (nbr2Atom->getIdx() == atom->getIdx()) {
                                        if neighbor2.atom_index == atom.id().index() {
                                            // RDKit✔️✔️:                       continue;
                                            continue;
                                        }
                                        // RDKit✔️✔️:                     }
                                        // RDKit✔️✔️:                     doubleBondedCN = (!((nbr2Atom->getAtomicNum() == 7) &&
                                        // RDKit✔️✔️:                                         (nbr2Atom->getTotalDegree() == 3)));
                                        double_bonded_cn = !(nbr2_atom.atomic_number() == 7
                                            && mmff_total_degree(
                                                mol,
                                                implicit_hydrogens,
                                                neighbor2.atom_index,
                                            )? == 3);
                                    }
                                    // RDKit✔️✔️:                   }
                                }
                                // RDKit✔️✔️:                 }
                            }
                            // RDKit✔️✔️:               }
                        }
                        // RDKit✔️✔️:             }
                        // RDKit❗✔️:             // if there is a single terminal oxygen
                        // RDKit❗✔️:             if (nTermObondedToN == 1) {
                        if n_term_o_bonded_to_n == 1 {
                            // RDKit❗✔️:               // N2OX
                            // RDKit❗✔️:               // sp2-hybridized N-oxide nitrogen
                            // RDKit❗✔️:               atomType = 67;
                            atom_type = 67;
                            // RDKit❗✔️:               break;
                            // RDKit❗✔️:             }
                            // RDKit❗✔️:             // if there are two or more terminal oxygens
                            // RDKit❗✔️:             if (nTermObondedToN >= 2) {
                        } else if n_term_o_bonded_to_n >= 2 {
                            // RDKit❗✔️:               // NO2
                            // RDKit❗✔️:               // Nitrogen in nitro group
                            // RDKit❗✔️:               // NO3
                            // RDKit❗✔️:               // Nitrogen in nitrate group
                            // RDKit❗✔️:               atomType = 45;
                            atom_type = 45;
                            // RDKit❗✔️:               break;
                        }
                        // RDKit❗✔️:             }
                        // RDKit✔️✔️:             // if the carbon bonded to ipso is bonded to 1 nitrogen
                        // RDKit✔️✔️:             // with 3 neighbors, that nitrogen is ipso (>N+=C)
                        // RDKit✔️✔️:             // alternatively, if there is no carbon but ipso is
                        // RDKit✔️✔️:             // double bonded to nitrogen, we have >N+=N
                        // RDKit✔️✔️:             if (doubleBondedCN) {
                        if atom_type == 0 && double_bonded_cn {
                            // RDKit✔️✔️:               // N+=C
                            // RDKit✔️✔️:               // Iminium nitrogen
                            // RDKit✔️✔️:               // N+=N
                            // RDKit✔️✔️:               // Positively charged nitrogen doubly bonded to N
                            // RDKit✔️✔️:               atomType = 54;
                            atom_type = 54;
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             }
                    }
                    // RDKit❗✔️:           }
                    // The following total-bond-order >= 3 delocalized-lone-pair
                    // nitrogen branch is ported only for the source-visible
                    // carbon-neighbor path below. Nitrogen-neighbor,
                    // amidinium, guanidinium, sulfur, phosphorus, and aromatic
                    // ring-cardinality paths remain unsupported until their
                    // full RDKit source blocks are ported in place.
                    // RDKit❗✔️:           // total bond order >= 3
                    // RDKit❗✔️:           if ((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit❗✔️:                atom->getNumImplicitHs()) >= 3) {
                    if atom_type == 0 && total_bond_order >= 3 {
                        // RDKit❗✔️:             bool isNCOorNCS = false;
                        let mut is_nco_or_ncs = false;
                        // RDKit✔️✔️:             bool isNCNplus = false;
                        let mut is_ncn_plus = false;
                        // RDKit✔️✔️:             bool isNGDplus = false;
                        let mut is_ngd_plus = false;
                        // RDKit✔️✔️:             bool isNNNorNNC = false;
                        let mut is_nnn_or_nnc = false;
                        // RDKit❗✔️:             bool isNbrC = false;
                        let mut is_nbr_c = false;
                        // RDKit❗✔️:             bool isNbrBenzeneC = false;
                        let mut is_nbr_benzene_c = false;
                        // RDKit❗✔️:             unsigned int elementDoubleBondedToC = 0;
                        let mut element_double_bonded_to_c = 0_u8;
                        // RDKit❗✔️:             unsigned int elementTripleBondedToC = 0;
                        let mut element_triple_bonded_to_c = 0_u8;
                        // RDKit❌❌:             unsigned int nN2bondedToC = 0;
                        // RDKit❌❌:             unsigned int nN3bondedToC = 0;
                        // RDKit❗✔️:             unsigned int nObondedToC = 0;
                        let mut n_o_bonded_to_c = 0_u32;
                        // RDKit❗✔️:             unsigned int nSbondedToC = 0;
                        let mut n_s_bonded_to_c = 0_u32;
                        // RDKit❗✔️:             // loop over neighbors
                        // RDKit❗✔️:             boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                        // RDKit❗✔️:             for (; nbrIdx != endNbrs; ++nbrIdx) {
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            // RDKit❗✔️:               const Atom *nbrAtom = mol[*nbrIdx];
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            // RDKit❗✔️:               // if the neighbor is carbon
                            // RDKit❗✔️:               if (nbrAtom->getAtomicNum() == 6) {
                            if nbr_atom.atomic_number() == 6 {
                                // RDKit❗✔️:                 isNbrC = true;
                                is_nbr_c = true;
                                // RDKit❗✔️:                 // check if we have a benzene carbon close to ipso
                                // RDKit❗✔️:                 if (nbrAtom->getIsAromatic() &&
                                // RDKit❗✔️:                     ringInfo->isAtomInRingOfSize(nbrAtom->getIdx(), 6)) {
                                if nbr_atom.is_aromatic()
                                    && ring_info.is_atom_in_ring_of_size(
                                        AtomId::new(neighbor.atom_index),
                                        6,
                                    )
                                {
                                    // RDKit❗✔️:                   isNbrBenzeneC = true;
                                    is_nbr_benzene_c = true;
                                    // RDKit❗✔️:                 }
                                }
                                // RDKit✔️✔️:                 nN2bondedToC = 0;
                                let mut n_n2_bonded_to_c = 0_u32;
                                // RDKit✔️✔️:                 nN3bondedToC = 0;
                                let mut n_n3_bonded_to_c = 0_u32;
                                // RDKit❗✔️:                 nObondedToC = 0;
                                n_o_bonded_to_c = 0;
                                // RDKit❗✔️:                 nSbondedToC = 0;
                                n_s_bonded_to_c = 0;
                                // RDKit✔️✔️:                 unsigned int nFormalCharge = 0;
                                let mut n_formal_charge = 0_u32;
                                // RDKit✔️✔️:                 unsigned int nInAromatic6Ring = 0;
                                let mut n_in_aromatic_6_ring = 0_u32;
                                // RDKit❗✔️:                 // loop over carbon neighbors
                                // RDKit❗✔️:                 boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                                // RDKit❗✔️:                 for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                                for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                                    // RDKit❗✔️:                   const Atom *nbr2Atom = mol[*nbr2Idx];
                                    let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                                    // RDKit❗✔️:                   const Bond *bond = mol.getBondBetweenAtoms(
                                    // RDKit❗✔️:                       nbrAtom->getIdx(), nbr2Atom->getIdx());
                                    let bond = &mol.bonds.as_slice()[neighbor2.bond.index()];
                                    // RDKit❗✔️:                   // check if we have oxygen or sulfur double-bonded to this
                                    // RDKit❗✔️:                   // carbon
                                    // RDKit❗✔️:                   if ((bond->getBondType() == Bond::DOUBLE) &&
                                    // RDKit❗✔️:                       ((nbr2Atom->getAtomicNum() == 8) ||
                                    // RDKit❗✔️:                        (nbr2Atom->getAtomicNum() == 16))) {
                                    if bond.order() == BondOrder::Double
                                        && matches!(nbr2_atom.atomic_number(), 8 | 16)
                                    {
                                        // RDKit❗✔️:                     isNCOorNCS = true;
                                        is_nco_or_ncs = true;
                                        // RDKit❗✔️:                   }
                                    }
                                    // RDKit❗✔️:                   // check if there is an atom double-bonded to this carbon,
                                    // RDKit❗✔️:                   // and if so find which element; if it is carbon or
                                    // RDKit❗✔️:                   // nitrogen (provided that the latter does not belong to
                                    // RDKit❗✔️:                   // multiple rings), also an aromatic bond is acceptable
                                    // RDKit❗✔️:                   if ((bond->getBondType() == Bond::DOUBLE) ||
                                    // RDKit✔️✔️:                       ((bond->getBondType() == Bond::AROMATIC) &&
                                    // RDKit✔️✔️:                        ((nbr2Atom->getAtomicNum() == 6) ||
                                    // RDKit✔️✔️:                         ((nbr2Atom->getAtomicNum() == 7) &&
                                    // RDKit✔️✔️:                          (queryIsAtomInNRings(nbr2Atom) == 1))))) {
                                    if bond.order() == BondOrder::Double
                                        || (bond.is_aromatic()
                                            && (nbr2_atom.atomic_number() == 6
                                                || (nbr2_atom.atomic_number() == 7
                                                    && ring_info.num_atom_rings(AtomId::new(
                                                        neighbor2.atom_index,
                                                    )) == 1)))
                                    {
                                        // RDKit✔️✔️:                     elementDoubleBondedToC = nbr2Atom->getAtomicNum();
                                        element_double_bonded_to_c = nbr2_atom.atomic_number();
                                        // RDKit❗✔️:                   }
                                    }
                                    // RDKit❗✔️:                   // check there is an atom triple-bonded to this carbon,
                                    // RDKit❗✔️:                   // and if so find which element
                                    // RDKit❗✔️:                   if (bond->getBondType() == Bond::TRIPLE) {
                                    if bond.order() == BondOrder::Triple {
                                        // RDKit❗✔️:                     elementTripleBondedToC = nbr2Atom->getAtomicNum();
                                        element_triple_bonded_to_c = nbr2_atom.atomic_number();
                                        // RDKit❗✔️:                   }
                                    }
                                    // RDKit✔️✔️:                   // if this carbon is bonded to a nitrogen with 3 neighbors
                                    // RDKit✔️✔️:                   if ((nbr2Atom->getAtomicNum() == 7) &&
                                    // RDKit✔️✔️:                       (nbr2Atom->getTotalDegree() == 3)) {
                                    if nbr2_atom.atomic_number() == 7
                                        && mmff_total_degree(
                                            mol,
                                            implicit_hydrogens,
                                            neighbor2.atom_index,
                                        )? == 3
                                    {
                                        // RDKit✔️✔️:                     // count the number of +1 formal charges that we have
                                        // RDKit✔️✔️:                     if (nbr2Atom->getFormalCharge() == 1) {
                                        if nbr2_atom.formal_charge() == 1 {
                                            // RDKit✔️✔️:                       ++nFormalCharge;
                                            n_formal_charge += 1;
                                            // RDKit✔️✔️:                     }
                                        }
                                        // RDKit✔️✔️:                     if (rmSize.isAtomInAromaticRingOfSize(nbrAtom, 6)) {
                                        if mmff_is_atom_in_aromatic_ring_of_size(
                                            mol,
                                            ring_info,
                                            AtomId::new(neighbor.atom_index),
                                            6,
                                        ) {
                                            // RDKit✔️✔️:                       ++nInAromatic6Ring;
                                            n_in_aromatic_6_ring += 1;
                                            // RDKit✔️✔️:                     }
                                        }
                                        // RDKit✔️✔️:                     // count how many oxygens are bonded to this nitrogen
                                        // RDKit✔️✔️:                     // with 3 neighbors
                                        // RDKit✔️✔️:                     unsigned int nObondedToN3 = 0;
                                        let mut n_o_bonded_to_n3 = 0_u32;
                                        // RDKit✔️✔️:                     boost::tie(nbr3Idx, end3Nbrs) =
                                        // RDKit✔️✔️:                         mol.getAtomNeighbors(nbr2Atom);
                                        // RDKit✔️✔️:                     for (; nbr3Idx != end3Nbrs; ++nbr3Idx) {
                                        for neighbor3 in
                                            mol.adjacency.neighbors_of(neighbor2.atom_index)
                                        {
                                            // RDKit✔️✔️:                       const Atom *nbr3Atom = mol[*nbr3Idx];
                                            let nbr3_atom =
                                                &mol.atoms.as_slice()[neighbor3.atom_index];
                                            // RDKit✔️✔️:                       if (nbr3Atom->getAtomicNum() == 8) {
                                            if nbr3_atom.atomic_number() == 8 {
                                                // RDKit✔️✔️:                         ++nObondedToN3;
                                                n_o_bonded_to_n3 += 1;
                                                // RDKit✔️✔️:                       }
                                            }
                                            // RDKit✔️✔️:                     }
                                        }
                                        // RDKit✔️✔️:                     // if there are less than 2 oxygens, this is neither
                                        // RDKit✔️✔️:                     // a nitro group nor a nitrate, so increment the counter
                                        // RDKit✔️✔️:                     // of nitrogens with 3 neighbors bonded to this carbon
                                        // RDKit✔️✔️:                     // (C-N<)
                                        // RDKit✔️✔️:                     if (nObondedToN3 < 2) {
                                        if n_o_bonded_to_n3 < 2 {
                                            // RDKit✔️✔️:                       ++nN3bondedToC;
                                            n_n3_bonded_to_c += 1;
                                            // RDKit✔️✔️:                     }
                                        }
                                        // RDKit✔️✔️:                   }
                                    }
                                    // RDKit✔️✔️:                   // if this carbon is bonded to a nitrogen with 2 neighbors
                                    // RDKit✔️✔️:                   // via a double or aromatic bond, increment the counter
                                    // RDKit✔️✔️:                   // of nitrogens with 2 neighbors bonded to this carbon
                                    // RDKit✔️✔️:                   // via a double or aromatic bond (C=N-)
                                    // RDKit✔️✔️:                   if ((nbr2Atom->getAtomicNum() == 7) &&
                                    // RDKit✔️✔️:                       (nbr2Atom->getTotalDegree() == 2) &&
                                    // RDKit✔️✔️:                       ((bond->getBondType() == Bond::DOUBLE) ||
                                    // RDKit✔️✔️:                        (bond->getBondType() == Bond::AROMATIC))) {
                                    if nbr2_atom.atomic_number() == 7
                                        && mmff_total_degree(
                                            mol,
                                            implicit_hydrogens,
                                            neighbor2.atom_index,
                                        )? == 2
                                        && (bond.order() == BondOrder::Double || bond.is_aromatic())
                                    {
                                        // RDKit✔️✔️:                     ++nN2bondedToC;
                                        n_n2_bonded_to_c += 1;
                                        // RDKit✔️✔️:                   }
                                    }
                                    // RDKit❗✔️:                   // if this carbon is bonded to an aromatic atom
                                    // RDKit❗✔️:                   if (nbr2Atom->getIsAromatic()) {
                                    if nbr2_atom.is_aromatic() {
                                        // RDKit❗✔️:                     // if it is oxygen, increment the counter of
                                        // RDKit❗✔️:                     // aromatic oxygen atoms bonded to this carbon
                                        // RDKit❗✔️:                     if (nbr2Atom->getAtomicNum() == 8) {
                                        if nbr2_atom.atomic_number() == 8 {
                                            // RDKit❗✔️:                       ++nObondedToC;
                                            n_o_bonded_to_c += 1;
                                            // RDKit❗✔️:                     }
                                        }
                                        // RDKit❗✔️:                     // if it is sulfur, increment the counter of
                                        // RDKit❗✔️:                     // aromatic sulfur atoms bonded to this carbon
                                        // RDKit❗✔️:                     if (nbr2Atom->getAtomicNum() == 16) {
                                        if nbr2_atom.atomic_number() == 16 {
                                            // RDKit❗✔️:                       ++nSbondedToC;
                                            n_s_bonded_to_c += 1;
                                            // RDKit❗✔️:                     }
                                        }
                                        // RDKit❗✔️:                   }
                                    }
                                }
                                // RDKit❗✔️:                 }
                                // RDKit✔️✔️:                 // if nitrogen is bonded to this carbon via a double or aromatic
                                // RDKit✔️✔️:                 // bond
                                // RDKit✔️✔️:                 if (elementDoubleBondedToC == 7) {
                                if element_double_bonded_to_c == 7 {
                                    // RDKit✔️✔️:                   // if 2 nitrogens with 3 neighbors and no nitrogens with 2
                                    // RDKit✔️✔️:                   // neighbors
                                    // RDKit✔️✔️:                   // are bonded to this carbon, and we have a formal charge,
                                    // RDKit✔️✔️:                   // but not a 6-membered aromatic ring, and the carbon atom
                                    // RDKit✔️✔️:                   // is not sp3, then this is an amidinium nitrogen (>N-C=N+<)
                                    // RDKit✔️✔️:                   if ((nN3bondedToC == 2) && (!nN2bondedToC) && nFormalCharge &&
                                    // RDKit✔️✔️:                       (!nInAromatic6Ring) && (nbrAtom->getTotalDegree() < 4)) {
                                    if n_n3_bonded_to_c == 2
                                        && n_n2_bonded_to_c == 0
                                        && n_formal_charge != 0
                                        && n_in_aromatic_6_ring == 0
                                        && mmff_total_degree(
                                            mol,
                                            implicit_hydrogens,
                                            neighbor.atom_index,
                                        )? < 4
                                    {
                                        // RDKit✔️✔️:                     isNCNplus = true;
                                        is_ncn_plus = true;
                                        // RDKit✔️✔️:                   }
                                    }
                                    // RDKit✔️✔️:                   // if 3 nitrogens with 3 neighbors are bonded
                                    // RDKit✔️✔️:                   // to this carbon, then this is a guanidinium nitrogen
                                    // RDKit✔️✔️:                   // ((>N-)2-C=N+<)
                                    // RDKit✔️✔️:                   if (nN3bondedToC == 3) {
                                    if n_n3_bonded_to_c == 3 {
                                        // RDKit✔️✔️:                     isNGDplus = true;
                                        is_ngd_plus = true;
                                        // RDKit✔️✔️:                   }
                                    }
                                    // RDKit✔️✔️:                 }
                                }
                            }
                            // RDKit❗✔️:               }
                            // RDKit✔️✔️:               // if the neighbor is nitrogen
                            // RDKit✔️✔️:               if (nbrAtom->getAtomicNum() == 7) {
                            if nbr_atom.atomic_number() == 7 {
                                // RDKit✔️✔️:                 unsigned int nNbondedToN = 0;
                                let mut n_n_bonded_to_n = 0_u32;
                                // RDKit✔️✔️:                 unsigned int nObondedToN = 0;
                                let mut n_o_bonded_to_n = 0_u32;
                                // RDKit✔️✔️:                 unsigned int nSbondedToN = 0;
                                let mut n_s_bonded_to_n = 0_u32;
                                // RDKit✔️✔️:                 // loop over nitrogen neighbors
                                // RDKit✔️✔️:                 boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                                // RDKit✔️✔️:                 for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                                for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                                    // RDKit✔️✔️:                   const Atom *nbr2Atom = mol[*nbr2Idx];
                                    let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                                    // RDKit✔️✔️:                   const Bond *bond = mol.getBondBetweenAtoms(
                                    // RDKit✔️✔️:                       nbrAtom->getIdx(), nbr2Atom->getIdx());
                                    let bond = &mol.bonds.as_slice()[neighbor2.bond.index()];
                                    // RDKit✔️✔️:                   // if the bond to nitrogen is double
                                    // RDKit✔️✔️:                   if (bond->getBondType() == Bond::DOUBLE) {
                                    if bond.order() == BondOrder::Double {
                                        // RDKit✔️✔️:                     // if the neighbor is carbon (N=N-C)
                                        // RDKit✔️✔️:                     if (nbr2Atom->getAtomicNum() == 6) {
                                        if nbr2_atom.atomic_number() == 6 {
                                            // RDKit✔️✔️:                       // loop over carbon neighbors
                                            // RDKit✔️✔️:                       boost::tie(nbr3Idx, end3Nbrs) =
                                            // RDKit✔️✔️:                           mol.getAtomNeighbors(nbr2Atom);
                                            // RDKit✔️✔️:                       for (; nbr3Idx != end3Nbrs; ++nbr3Idx) {
                                            for neighbor3 in
                                                mol.adjacency.neighbors_of(neighbor2.atom_index)
                                            {
                                                // RDKit✔️✔️:                         const Atom *nbr3Atom = mol[*nbr3Idx];
                                                let nbr3_atom =
                                                    &mol.atoms.as_slice()[neighbor3.atom_index];
                                                // RDKit✔️✔️:                         // if the nitrogen neighbor to ipso is met, move on
                                                // RDKit✔️✔️:                         if (nbr3Atom->getIdx() == nbrAtom->getIdx()) {
                                                if neighbor3.atom_index == neighbor.atom_index {
                                                    // RDKit✔️✔️:                           continue;
                                                    continue;
                                                }
                                                // RDKit✔️✔️:                         }
                                                // RDKit✔️✔️:                         // count how many nitrogen, oxygen, sulfur atoms
                                                // RDKit✔️✔️:                         // are bonded to this carbon
                                                // RDKit✔️✔️:                         switch (nbr3Atom->getAtomicNum()) {
                                                match nbr3_atom.atomic_number() {
                                                    // RDKit✔️✔️:                           case 7:
                                                    7 => {
                                                        // RDKit✔️✔️:                             ++nNbondedToN;
                                                        n_n_bonded_to_n += 1;
                                                    }
                                                    // RDKit✔️✔️:                             break;
                                                    // RDKit✔️✔️:                           case 8:
                                                    8 => {
                                                        // RDKit✔️✔️:                             ++nObondedToN;
                                                        n_o_bonded_to_n += 1;
                                                    }
                                                    // RDKit✔️✔️:                             break;
                                                    // RDKit✔️✔️:                           case 16:
                                                    16 => {
                                                        // RDKit✔️✔️:                             ++nSbondedToN;
                                                        n_s_bonded_to_n += 1;
                                                    }
                                                    // RDKit✔️✔️:                             break;
                                                    _ => {}
                                                }
                                                // RDKit✔️✔️:                         }
                                            }
                                            // RDKit✔️✔️:                       }
                                            // RDKit✔️✔️:                       // if there are no more nitrogens, no oxygen, no sulfur
                                            // RDKit✔️✔️:                       // bonded
                                            // RDKit✔️✔️:                       // to carbon, and the latter is not a benzene carbon
                                            // RDKit✔️✔️:                       // then it is N=N-C
                                            // RDKit✔️✔️:                       if ((!nObondedToN) && (!nSbondedToN) && (!nNbondedToN) &&
                                            // RDKit✔️✔️:                           (!isNbrBenzeneC)) {
                                            if n_o_bonded_to_n == 0
                                                && n_s_bonded_to_n == 0
                                                && n_n_bonded_to_n == 0
                                                && !is_nbr_benzene_c
                                            {
                                                // RDKit✔️✔️:                         isNNNorNNC = true;
                                                is_nnn_or_nnc = true;
                                            }
                                            // RDKit✔️✔️:                       }
                                        }
                                        // RDKit✔️✔️:                     }
                                        // RDKit✔️✔️:                     // if the neighbor is nitrogen (N=N-N) and ipso is not
                                        // RDKit✔️✔️:                     // bonded
                                        // RDKit✔️✔️:                     // to benzene carbon then it is N=N-N
                                        // RDKit✔️✔️:                     if ((nbr2Atom->getAtomicNum() == 7) && (!isNbrBenzeneC)) {
                                        if nbr2_atom.atomic_number() == 7 && !is_nbr_benzene_c {
                                            // RDKit✔️✔️:                       isNNNorNNC = true;
                                            is_nnn_or_nnc = true;
                                        }
                                        // RDKit✔️✔️:                     }
                                    }
                                    // RDKit✔️✔️:                   }
                                }
                                // RDKit✔️✔️:                 }
                            }
                        }
                        // RDKit❗✔️:             }
                        // RDKit❗✔️:             // if ipso nitrogen is bonded to carbon
                        // RDKit❗✔️:             if (isNbrC) {
                        if is_nbr_c {
                            // RDKit❗✔️:               // if neighbor carbon is triple-bonded to N, then ipso is N-C%N
                            // RDKit❗✔️:               if (elementTripleBondedToC == 7) {
                            // RDKit❗✔️:                 isNSO2orNSO3orNCN = true;
                            if element_triple_bonded_to_c == 7 {
                                is_nso2_or_nso3_or_ncn = true;
                            }
                            // RDKit❗✔️:               }
                            // RDKit✔️✔️:               // if neighbor carbon is amidinium
                            // RDKit✔️✔️:               if (isNCNplus) {
                            if is_ncn_plus {
                                // RDKit✔️✔️:                 // NCN+
                                // RDKit✔️✔️:                 // Either nitrogen in N+=C-N
                                // RDKit✔️✔️:                 atomType = 55;
                                atom_type = 55;
                                // RDKit✔️✔️:                 break;
                            }
                            // RDKit✔️✔️:               }
                            // RDKit✔️✔️:               // if neighbor carbon is guanidinium
                            // RDKit✔️✔️:               if (isNGDplus) {
                            if atom_type == 0 && is_ngd_plus {
                                // RDKit✔️✔️:                 // NGD+
                                // RDKit✔️✔️:                 // Guanidinium nitrogen
                                // RDKit✔️✔️:                 atomType = 56;
                                atom_type = 56;
                                // RDKit✔️✔️:                 break;
                            }
                            // RDKit✔️✔️:               }
                            // RDKit❗✔️:               // if neighbor carbon is not bonded to oxygen or sulfur
                            // RDKit❗✔️:               // and is not cyano, there two possibilities:
                            // RDKit❗✔️:               // 1) ipso nitrogen is bonded to benzene carbon while no oxygen
                            // RDKit❗✔️:               //    or sulfur are bonded to the latter: ipso is aniline nitrogen
                            // RDKit❗✔️:               // 2) ipso nitrogen is bonded to a carbon which is double-bonded
                            // RDKit❗✔️:               // to
                            // RDKit❗✔️:               //    carbon, nitrogen or phosphorus, or triple-bonded to carbon
                            // RDKit❗✔️:               if (((!isNCOorNCS) && (!isNSO2orNSO3orNCN)) &&
                            // RDKit❗✔️:                   (((!nObondedToC) && (!nSbondedToC) && isNbrBenzeneC) ||
                            // RDKit❗✔️:                    ((elementDoubleBondedToC == 6) ||
                            // RDKit❗✔️:                     (elementDoubleBondedToC == 7) ||
                            // RDKit❗✔️:                     (elementDoubleBondedToC == 15) ||
                            // RDKit❗✔️:                     (elementTripleBondedToC == 6)))) {
                            if atom_type == 0
                                && ((!is_nco_or_ncs) && (!is_nso2_or_nso3_or_ncn))
                                && (((n_o_bonded_to_c == 0)
                                    && (n_s_bonded_to_c == 0)
                                    && is_nbr_benzene_c)
                                    || matches!(element_double_bonded_to_c, 6 | 7 | 15)
                                    || element_triple_bonded_to_c == 6)
                            {
                                // RDKit❗✔️:                 // NC=C
                                // RDKit❗✔️:                 // Enamine or aniline nitrogen, deloc. lp
                                // RDKit❗✔️:                 // NC=N
                                // RDKit❗✔️:                 // Nitrogen in N-C=N with deloc. lp
                                // RDKit❗✔️:                 // NC=P
                                // RDKit❗✔️:                 // Nitrogen in N-C=P with deloc. lp
                                // RDKit❗✔️:                 // NC%C
                                // RDKit❗✔️:                 // Nitrogen attached to C-C triple bond
                                // RDKit❗✔️:                 atomType = 40;
                                atom_type = 40;
                                // RDKit❗✔️:                 break;
                                // RDKit❗✔️:               }
                            }
                        }
                        // RDKit❗✔️:             }
                        // RDKit❗✔️:             // if ipso is not sulfonamide while it is either amide/thioamide
                        // RDKit❗✔️:             // or >N-N=N-/>N-N=C<
                        // RDKit❗✔️:             if ((!isNSO2orNSO3orNCN) && (isNCOorNCS || isNNNorNNC)) {
                        if atom_type == 0
                            && !is_nso2_or_nso3_or_ncn
                            && (is_nco_or_ncs || is_nnn_or_nnc)
                        {
                            // RDKit❗✔️:               // NC=O
                            // RDKit❗✔️:               // Amide nitrogen
                            // RDKit❗✔️:               // NC=S
                            // RDKit❗✔️:               // Thioamide nitrogen
                            // RDKit❗✔️:               // NN=C
                            // RDKit❗✔️:               // Nitrogen in N-N=C moiety with deloc. lp
                            // RDKit❗✔️:               // NN=N
                            // RDKit❗✔️:               // Nitrogen in N-N=N moiety with deloc. lp
                            // RDKit❗✔️:               atomType = 10;
                            atom_type = 10;
                            // RDKit❗✔️:               break;
                        }
                        // RDKit❗✔️:             }
                        // RDKit❗✔️:           }
                    }
                    // RDKit❗✔️:         }
                }
                // RDKit❗✔️:         // 2 neighbors
                // RDKit❗✔️:         if (atom->getTotalDegree() == 2) {
                if atom_type == 0 && total_degree == 2 {
                    // RDKit❗✔️:           // total bond order = 4
                    // RDKit❗✔️:           if ((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit❗✔️:                atom->getNumImplicitHs()) == 4) {
                    if total_bond_order == 4 {
                        // RDKit❗✔️:             // loop over neighbors
                        // RDKit❗✔️:             bool isIsonitrile = false;
                        let mut is_isonitrile = false;
                        // RDKit❗✔️:             boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                        // RDKit❗✔️:             for (; (!isIsonitrile) && (nbrIdx != endNbrs); ++nbrIdx) {
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            if is_isonitrile {
                                break;
                            }
                            let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                            // RDKit❗✔️:               const Atom *nbrAtom = mol[*nbrIdx];
                            // RDKit❗✔️:               // if neighbor is triple-bonded
                            // RDKit❗✔️:               isIsonitrile =
                            // RDKit❗✔️:                   ((mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx()))
                            // RDKit❗✔️:                        ->getBondType() == Bond::TRIPLE);
                            is_isonitrile = nbr_bond.order() == BondOrder::Triple;
                        }
                        // RDKit❗✔️:             }
                        // RDKit❗✔️:             if (isIsonitrile) {
                        if is_isonitrile {
                            // RDKit❗✔️:               // NR%
                            // RDKit❗✔️:               // Isonitrile nitrogen
                            // RDKit❗✔️:               atomType = 61;
                            atom_type = 61;
                            // RDKit❗✔️:               break;
                        } else {
                            // RDKit❗✔️:             }
                            // RDKit❗✔️:             // =N=
                            // RDKit❗✔️:             // Central nitrogen in C=N=N or N=N=N
                            // RDKit❗✔️:             atomType = 53;
                            atom_type = 53;
                            // RDKit❗✔️:             break;
                        }
                        // RDKit❗✔️:           }
                    }
                    // RDKit❗✔️:           // total bond order = 3
                    // RDKit❗✔️:           if ((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit❗✔️:                atom->getNumImplicitHs()) == 3) {
                    if atom_type == 0 && total_bond_order == 3 {
                        // RDKit✔️✔️:             // loop over neighbors
                        // RDKit✔️✔️:             bool isNitroso = false;
                        let mut is_nitroso = false;
                        // RDKit✔️✔️:             bool isImineOrAzo = false;
                        let mut is_imine_or_azo = false;
                        // RDKit✔️✔️:             boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                        // RDKit✔️✔️:             for (; nbrIdx != endNbrs; ++nbrIdx) {
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            // RDKit✔️✔️:               const Atom *nbrAtom = mol[*nbrIdx];
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                            // RDKit✔️✔️:               // if the neighbor is double bonded (-N=)
                            // RDKit✔️✔️:               if ((mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx()))
                            // RDKit✔️✔️:                       ->getBondType() == Bond::DOUBLE) {
                            if nbr_bond.order() == BondOrder::Double {
                                // RDKit✔️✔️:                 // if it is terminal oxygen (-N=O)
                                // RDKit✔️✔️:                 isNitroso =
                                // RDKit✔️✔️:                     ((nbrAtom->getAtomicNum() == 8) && (nTermObondedToN == 1));
                                is_nitroso =
                                    nbr_atom.atomic_number() == 8 && n_term_o_bonded_to_n == 1;
                                // RDKit✔️✔️:                 // if it is carbon or nitrogen (-N=N-, -N=C<),
                                // RDKit✔️✔️:                 // ipso is imine or azo
                                // RDKit✔️✔️:                 isImineOrAzo = ((nbrAtom->getAtomicNum() == 6) ||
                                // RDKit✔️✔️:                                 (nbrAtom->getAtomicNum() == 7));
                                is_imine_or_azo = matches!(nbr_atom.atomic_number(), 6 | 7);
                            }
                            // RDKit✔️✔️:               }
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             if (isNitroso && (!isImineOrAzo)) {
                        if is_nitroso && !is_imine_or_azo {
                            // RDKit✔️✔️:               // N=O
                            // RDKit✔️✔️:               // Nitrogen in nitroso group
                            // RDKit✔️✔️:               atomType = 46;
                            atom_type = 46;
                            // RDKit✔️✔️:               break;
                            // RDKit✔️✔️:             }
                        } else if is_imine_or_azo {
                            // RDKit✔️✔️:             if (isImineOrAzo) {
                            // RDKit✔️✔️:               // N=C
                            // RDKit✔️✔️:               // Imine nitrogen
                            // RDKit✔️✔️:               // N=N
                            // RDKit✔️✔️:               // Azo-group nitrogen
                            // RDKit✔️✔️:               atomType = 9;
                            atom_type = 9;
                            // RDKit✔️✔️:               break;
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           // total bond order >= 2
                    // RDKit✔️✔️:           if ((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit✔️✔️:                atom->getNumImplicitHs()) >= 2) {
                    if atom_type == 0 && total_bond_order >= 2 {
                        // RDKit✔️✔️:             // loop over neighbors
                        // RDKit✔️✔️:             bool isNSO = false;
                        let mut is_nso = false;
                        // RDKit✔️✔️:             boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                        // RDKit✔️✔️:             for (; (!isNSO) && (nbrIdx != endNbrs); ++nbrIdx) {
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            if is_nso {
                                break;
                            }
                            // RDKit✔️✔️:               const Atom *nbrAtom = mol[*nbrIdx];
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            // RDKit✔️✔️:               // if the neighbor is sulfur bonded to a single terminal oxygen
                            // RDKit✔️✔️:               if (nbrAtom->getAtomicNum() == 16) {
                            if nbr_atom.atomic_number() == 16 {
                                // RDKit✔️✔️:                 // loop over neighbors and count how many
                                // RDKit✔️✔️:                 // terminal oxygens are bonded to sulfur
                                // RDKit✔️✔️:                 unsigned int nTermObondedToS = 0;
                                let mut n_term_o_bonded_to_s = 0_u32;
                                // RDKit✔️✔️:                 boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                                // RDKit✔️✔️:                 for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                                for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                                    // RDKit✔️✔️:                   const Atom *nbr2Atom = mol[*nbr2Idx];
                                    let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                                    // RDKit✔️✔️:                   if ((nbr2Atom->getAtomicNum() == 8) &&
                                    // RDKit✔️✔️:                       (nbr2Atom->getTotalDegree() == 1)) {
                                    if nbr2_atom.atomic_number() == 8
                                        && mmff_total_degree(
                                            mol,
                                            implicit_hydrogens,
                                            neighbor2.atom_index,
                                        )? == 1
                                    {
                                        // RDKit✔️✔️:                     ++nTermObondedToS;
                                        n_term_o_bonded_to_s += 1;
                                    }
                                    // RDKit✔️✔️:                   }
                                }
                                // RDKit✔️✔️:                 }
                                // RDKit✔️✔️:                 isNSO = (nTermObondedToS == 1);
                                is_nso = n_term_o_bonded_to_s == 1;
                            }
                            // RDKit✔️✔️:               }
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             if (isNSO) {
                        if is_nso {
                            // RDKit✔️✔️:               // NSO
                            // RDKit✔️✔️:               // Divalent nitrogen replacing monovalent O in SO2 group
                            // RDKit✔️✔️:               atomType = 48;
                            atom_type = 48;
                            // RDKit✔️✔️:               break;
                            // RDKit✔️✔️:             }
                        } else if !is_nso2_or_nso3_or_ncn {
                            // RDKit✔️✔️:             if (!isNSO2orNSO3orNCN) {
                            // RDKit✔️✔️:               // If it is not sulfonamide deprotonated nitrogen,
                            // RDKit✔️✔️:               // it is anionic nitrogen (>N::-)
                            // RDKit✔️✔️:               // NM
                            // RDKit✔️✔️:               // Anionic divalent nitrogen
                            // RDKit✔️✔️:               atomType = 62;
                            atom_type = 62;
                            // RDKit✔️✔️:               break;
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit❗✔️:         }
                }
                // RDKit✔️✔️:         // if it is sulfonamide (3 neighbors) or cyano (2 neighbors)
                // RDKit✔️✔️:         if (isNSO2orNSO3orNCN) {
                if atom_type == 0 && is_nso2_or_nso3_or_ncn {
                    // RDKit✔️✔️:           // NSO2
                    // RDKit✔️✔️:           // Sulfonamide nitrogen
                    // RDKit✔️✔️:           // NSO3
                    // RDKit✔️✔️:           // Sulfonamide nitrogen
                    // RDKit✔️✔️:           // NC%N
                    // RDKit✔️✔️:           // Nitrogen attached to cyano group
                    // RDKit✔️✔️:           atomType = 43;
                    atom_type = 43;
                    // RDKit✔️✔️:           break;
                    // RDKit✔️✔️:         }
                }
                // RDKit❗✔️:         // 1 neighbor
                // RDKit❗✔️:         if (atom->getTotalDegree() == 1) {
                if atom_type == 0 && total_degree == 1 {
                    // RDKit❗✔️:           bool isNSP = false;
                    let mut is_nsp = false;
                    // RDKit❗✔️:           bool isNAZT = false;
                    let mut is_nazt = false;
                    // RDKit❗✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit❗✔️:           for (; (!isNSP) && (!isNAZT) && (nbrIdx != endNbrs); ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        if is_nsp || is_nazt {
                            break;
                        }
                        // RDKit❗✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                        // RDKit❗✔️:             // if ipso is triple-bonded to its only neighbor
                        // RDKit❗✔️:             isNSP =
                        // RDKit❗✔️:                 ((mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx()))
                        // RDKit❗✔️:                      ->getBondType() == Bond::TRIPLE);
                        is_nsp = nbr_bond.order() == BondOrder::Triple;
                        // RDKit❗✔️:             // ipso is bonded to a nitrogen atom with 2 neighbors
                        // RDKit❗✔️:             if ((nbrAtom->getAtomicNum() == 7) &&
                        // RDKit❗✔️:                 (nbrAtom->getTotalDegree() == 2)) {
                        if nbr_atom.atomic_number() == 7
                            && mmff_total_degree(mol, implicit_hydrogens, neighbor.atom_index)? == 2
                        {
                            // RDKit❗✔️:               // loop over nitrogen neighbors
                            // RDKit❗✔️:               boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                            // RDKit❗✔️:               for (; (!isNAZT) && (nbr2Idx != end2Nbrs); ++nbr2Idx) {
                            for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                                if is_nazt {
                                    break;
                                }
                                // RDKit❗✔️:                 const Atom *nbr2Atom = mol[*nbr2Idx];
                                let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                                let nbr2_total_degree = mmff_total_degree(
                                    mol,
                                    implicit_hydrogens,
                                    neighbor2.atom_index,
                                )?;
                                // RDKit❗✔️:                 // if another nitrogen with 2 neighbors, or a carbon
                                // RDKit❗✔️:                 // with 3 neighbors is found, ipso is NAZT
                                // RDKit❗✔️:                 isNAZT = (((nbr2Atom->getAtomicNum() == 7) &&
                                // RDKit❗✔️:                            (nbr2Atom->getTotalDegree() == 2)) ||
                                // RDKit❗✔️:                           ((nbr2Atom->getAtomicNum() == 6) &&
                                // RDKit❗✔️:                            (nbr2Atom->getTotalDegree() == 3)));
                                is_nazt = (nbr2_atom.atomic_number() == 7
                                    && nbr2_total_degree == 2)
                                    || (nbr2_atom.atomic_number() == 6 && nbr2_total_degree == 3);
                            }
                            // RDKit❗✔️:               }
                        }
                        // RDKit❗✔️:             }
                    }
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           if (isNSP) {
                    if is_nsp {
                        // RDKit❗✔️:             // NSP
                        // RDKit❗✔️:             // Triply bonded nitrogen
                        // RDKit❗✔️:             atomType = 42;
                        atom_type = 42;
                        // RDKit❗✔️:             break;
                    }
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           if (isNAZT) {
                    if atom_type == 0 && is_nazt {
                        // RDKit❗✔️:             // NAZT
                        // RDKit❗✔️:             // Terminal nitrogen in azido or diazo group
                        // RDKit❗✔️:             atomType = 47;
                        atom_type = 47;
                        // RDKit❗✔️:             break;
                    }
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:         }
                }
                // RDKit✔️✔️:         // if nothing else was found
                // RDKit✔️✔️:         // NR
                // RDKit✔️✔️:         // Amine nitrogen
                // RDKit✔️✔️:         atomType = 8;
                if atom_type == 0 {
                    atom_type = 8;
                }
                // RDKit✔️✔️:         break;
            }
            // RDKit❗✔️:       // Oxygen
            // RDKit❗✔️:       case 8:
            8 => {
                let total_degree = mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                let total_bond_order =
                    mmff_total_bond_order(explicit_valence, implicit_hydrogens, atom.id().index())?;
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         // 3 neighbors
                // RDKit✔️✔️:         if (atom->getTotalDegree() == 3) {
                if total_degree == 3 {
                    // RDKit✔️✔️:           // O+
                    // RDKit✔️✔️:           // Oxonium oxygen
                    // RDKit✔️✔️:           atomType = 49;
                    atom_type = 49;
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         }
                // RDKit❗✔️:         // 2 neighbors
                // RDKit❗✔️:         if (atom->getTotalDegree() == 2) {
                if atom_type == 0 && total_degree == 2 {
                    // RDKit❗✔️:           if ((atom->getValence(Atom::ValenceType::EXPLICIT) +
                    // RDKit❗✔️:                atom->getNumImplicitHs()) == 3) {
                    if total_bond_order == 3 {
                        // RDKit❗✔️:             // O=+
                        // RDKit❗✔️:             // Oxenium oxygen
                        // RDKit❗✔️:             atomType = 51;
                        atom_type = 51;
                        // RDKit❗✔️:             break;
                    } else {
                        // RDKit❗✔️:           }
                        // RDKit❗✔️:           // count how many hydrogens are bound to ipso
                        // RDKit❗✔️:           unsigned int nHbondedToO = 0;
                        let mut n_h_bonded_to_o = 0_i32;
                        // RDKit❗✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                        // RDKit❗✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
                        for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                            // RDKit❗✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                            let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                            // RDKit❗✔️:             if (nbrAtom->getAtomicNum() == 1) {
                            if nbr_atom.atomic_number() == 1 {
                                // RDKit❗✔️:               ++nHbondedToO;
                                n_h_bonded_to_o += 1;
                            }
                            // RDKit❗✔️:             }
                        }
                        // RDKit❗✔️:           }
                        let oxygen_implicit_hydrogens = implicit_hydrogens
                            .get(atom.id().index())
                            .copied()
                            .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                                atom_index: atom.id().index(),
                                atoms: implicit_hydrogens.len(),
                            })?;
                        // RDKit❗✔️:           if ((nHbondedToO + atom->getNumImplicitHs()) == 2) {
                        if n_h_bonded_to_o + oxygen_implicit_hydrogens == 2 {
                            // RDKit❗✔️:             // OH2
                            // RDKit❗✔️:             // Oxygen in water
                            // RDKit❗✔️:             atomType = 70;
                            atom_type = 70;
                            // RDKit❗✔️:             break;
                        } else {
                            // RDKit❗✔️:           }
                            // RDKit❗✔️:           // otherwise, ipso must be one of the following
                            // RDKit❗✔️:           // OC=O
                            // RDKit❗✔️:           // Carboxylic acid or ester oxygen
                            // RDKit❗✔️:           // OC=C
                            // RDKit❗✔️:           // Enolic or phenolic oxygen
                            // RDKit❗✔️:           // OC=N
                            // RDKit❗✔️:           // Oxygen in -O-C=N- moiety
                            // RDKit❗✔️:           // OC=S
                            // RDKit❗✔️:           // Divalent oxygen in thioacid or ester
                            // RDKit❗✔️:           // ONO2
                            // RDKit❗✔️:           // Divalent nitrate "ether" oxygen
                            // RDKit❗✔️:           // ON=O
                            // RDKit❗✔️:           // Divalent nitrate "ether" oxygen
                            // RDKit❗✔️:           // OSO3
                            // RDKit❗✔️:           // Divalent oxygen in sulfate group
                            // RDKit❗✔️:           // OSO2
                            // RDKit❗✔️:           // Divalent oxygen in sulfite group
                            // RDKit❗✔️:           // OSO
                            // RDKit❗✔️:           // One of two divalent oxygens attached to sulfur
                            // RDKit❗✔️:           // OS=O
                            // RDKit❗✔️:           // Divalent oxygen in R(RO)S=O
                            // RDKit❗✔️:           // -OS
                            // RDKit❗✔️:           // Other divalent oxygen attached to sulfur
                            // RDKit❗✔️:           // OPO3
                            // RDKit❗✔️:           // Divalent oxygen in phosphate group
                            // RDKit❗✔️:           // OPO2
                            // RDKit❗✔️:           // Divalent oxygen in phosphite group
                            // RDKit❗✔️:           // OPO
                            // RDKit❗✔️:           // Divalent oxygen, one of two oxygens attached to P
                            // RDKit❗✔️:           // -OP
                            // RDKit❗✔️:           // Other divalent oxygen attached to phosphorus
                            // RDKit❗✔️:           atomType = 6;
                            atom_type = 6;
                            // RDKit❗✔️:           break;
                        }
                    }
                    // RDKit❗✔️:         }
                }
                // RDKit❗✔️:         // 1 neighbor
                // RDKit❗✔️:         if (atom->getDegree() <= 1) {
                if atom_type == 0 && graph_degree <= 1 {
                    // RDKit✔️✔️:           unsigned int nNbondedToCorNorS = 0;
                    // RDKit✔️✔️:           unsigned int nObondedToCorNorS = 0;
                    // RDKit✔️✔️:           unsigned int nSbondedToCorNorS = 0;
                    let mut n_n_bonded_to_c_or_n_or_s = 0_u32;
                    let mut n_o_bonded_to_c_or_n_or_s = 0_u32;
                    let mut n_s_bonded_to_c_or_n_or_s = 0_u32;
                    let oxygen_implicit_hydrogens = implicit_hydrogens
                        .get(atom.id().index())
                        .copied()
                        .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                            atom_index: atom.id().index(),
                            atoms: implicit_hydrogens.len(),
                        })?;
                    // RDKit✔️✔️:           bool isOxideOBondedToH =
                    // RDKit✔️✔️:               atom->getNumExplicitHs() + atom->getNumImplicitHs() > 0;
                    let mut is_oxide_o_bonded_to_h =
                        i32::from(atom.explicit_hydrogens()) + oxygen_implicit_hydrogens > 0;
                    // RDKit✔️✔️:           bool isCarboxylateO = false;
                    let mut is_carboxylate_o = false;
                    // RDKit✔️✔️:           bool isCarbonylO = false;
                    let mut is_carbonyl_o = false;
                    // RDKit✔️✔️:           bool isOxideOBondedToC = false;
                    let mut is_oxide_o_bonded_to_c = false;
                    // RDKit✔️✔️:           bool isNitrosoO = false;
                    let mut is_nitroso_o = false;
                    // RDKit✔️✔️:           bool isOxideOBondedToN = false;
                    let mut is_oxide_o_bonded_to_n = false;
                    // RDKit✔️✔️:           bool isNOxideO = false;
                    let mut is_n_oxide_o = false;
                    // RDKit✔️✔️:           bool isNitroO = false;
                    let mut is_nitro_o = false;
                    // RDKit✔️✔️:           bool isThioSulfinateO = false;
                    let mut is_thio_sulfinate_o = false;
                    // RDKit✔️✔️:           bool isSulfateO = false;
                    let mut is_sulfate_o = false;
                    // RDKit✔️✔️:           bool isSulfoxideO = false;
                    let mut is_sulfoxide_o = false;
                    // RDKit✔️✔️:           bool isPhosphateOrPerchlorateO = false;
                    let mut is_phosphate_or_perchlorate_o = false;
                    // RDKit✔️✔️:           // loop over neighbors
                    // RDKit✔️✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit✔️✔️:           for (; (nbrIdx != endNbrs) && (!isOxideOBondedToC) &&
                    // RDKit✔️✔️:                  (!isOxideOBondedToN) && (!isOxideOBondedToH) &&
                    // RDKit✔️✔️:                  (!isCarboxylateO) && (!isNitroO) && (!isNOxideO) &&
                    // RDKit✔️✔️:                  (!isThioSulfinateO) && (!isSulfateO) &&
                    // RDKit✔️✔️:                  (!isPhosphateOrPerchlorateO) && (!isCarbonylO) &&
                    // RDKit✔️✔️:                  (!isNitrosoO) && (!isSulfoxideO);
                    // RDKit✔️✔️:                ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        if is_oxide_o_bonded_to_c
                            || is_oxide_o_bonded_to_n
                            || is_oxide_o_bonded_to_h
                            || is_carboxylate_o
                            || is_nitro_o
                            || is_n_oxide_o
                            || is_thio_sulfinate_o
                            || is_sulfate_o
                            || is_phosphate_or_perchlorate_o
                            || is_carbonyl_o
                            || is_nitroso_o
                            || is_sulfoxide_o
                        {
                            break;
                        }
                        // RDKit✔️✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        // RDKit✔️✔️:             const Bond *bond =
                        // RDKit✔️✔️:                 mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx());
                        let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                        // RDKit✔️✔️:             // if the neighbor is carbon, nitrogen or sulfur
                        // RDKit✔️✔️:             if ((nbrAtom->getAtomicNum() == 6) ||
                        // RDKit✔️✔️:                 (nbrAtom->getAtomicNum() == 7) ||
                        // RDKit✔️✔️:                 (nbrAtom->getAtomicNum() == 16)) {
                        if matches!(nbr_atom.atomic_number(), 6 | 7 | 16) {
                            // RDKit✔️✔️:               // count how many terminal oxygen/sulfur atoms
                            // RDKit✔️✔️:               // or secondary nitrogens
                            // RDKit✔️✔️:               // are bonded to the carbon or nitrogen neighbor of ipso
                            // RDKit✔️✔️:               boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                            // RDKit✔️✔️:               for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                            for nbr2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                                // RDKit✔️✔️:                 const Atom *nbr2Atom = mol[*nbr2Idx];
                                let nbr2_atom = &mol.atoms.as_slice()[nbr2.atom_index];
                                let nbr2_total_degree =
                                    mmff_total_degree(mol, implicit_hydrogens, nbr2.atom_index)?;
                                // RDKit✔️✔️:                 if ((nbr2Atom->getAtomicNum() == 7) &&
                                // RDKit✔️✔️:                     (nbr2Atom->getTotalDegree() == 2)) {
                                if nbr2_atom.atomic_number() == 7 && nbr2_total_degree == 2 {
                                    // RDKit✔️✔️:                   ++nNbondedToCorNorS;
                                    n_n_bonded_to_c_or_n_or_s += 1;
                                }
                                // RDKit✔️✔️:                 }
                                // RDKit✔️✔️:                 if ((nbr2Atom->getAtomicNum() == 8) &&
                                // RDKit✔️✔️:                     (nbr2Atom->getTotalDegree() == 1)) {
                                if nbr2_atom.atomic_number() == 8 && nbr2_total_degree == 1 {
                                    // RDKit✔️✔️:                   ++nObondedToCorNorS;
                                    n_o_bonded_to_c_or_n_or_s += 1;
                                }
                                // RDKit✔️✔️:                 }
                                // RDKit✔️✔️:                 if ((nbr2Atom->getAtomicNum() == 16) &&
                                // RDKit✔️✔️:                     (nbr2Atom->getTotalDegree() == 1)) {
                                if nbr2_atom.atomic_number() == 16 && nbr2_total_degree == 1 {
                                    // RDKit✔️✔️:                   ++nSbondedToCorNorS;
                                    n_s_bonded_to_c_or_n_or_s += 1;
                                }
                                // RDKit✔️✔️:                 }
                            }
                            // RDKit✔️✔️:               }
                        }
                        // RDKit✔️✔️:             // if ipso neighbor is hydrogen
                        // RDKit✔️✔️:             isOxideOBondedToH = (nbrAtom->getAtomicNum() == 1);
                        is_oxide_o_bonded_to_h = nbr_atom.atomic_number() == 1;
                        // RDKit✔️✔️:             // if ipso neighbor is carbon
                        // RDKit✔️✔️:             if (nbrAtom->getAtomicNum() == 6) {
                        if nbr_atom.atomic_number() == 6 {
                            // RDKit✔️✔️:               // if carbon neighbor is bonded to 2 oxygens,
                            // RDKit✔️✔️:               // ipso is carboxylate oxygen
                            // RDKit✔️✔️:               isCarboxylateO = (nObondedToCorNorS == 2);
                            is_carboxylate_o = n_o_bonded_to_c_or_n_or_s == 2;
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to carbon
                            // RDKit✔️✔️:               // via a double bond, ipso is carbonyl oxygen
                            // RDKit✔️✔️:               isCarbonylO = (bond->getBondType() == Bond::DOUBLE);
                            is_carbonyl_o = nbr_bond.order() == BondOrder::Double;
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to carbon via a
                            // RDKit✔️✔️:               // single bond, and there are no other bonded oxygens,
                            // RDKit✔️✔️:               // ipso is oxide oxygen
                            // RDKit✔️✔️:               isOxideOBondedToC = ((bond->getBondType() == Bond::SINGLE) &&
                            // RDKit✔️✔️:                                    (nObondedToCorNorS == 1));
                            is_oxide_o_bonded_to_c = nbr_bond.order() == BondOrder::Single
                                && n_o_bonded_to_c_or_n_or_s == 1;
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             // if ipso neighbor is nitrogen
                        // RDKit✔️✔️:             if (nbrAtom->getAtomicNum() == 7) {
                        if nbr_atom.atomic_number() == 7 {
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to nitrogen
                            // RDKit✔️✔️:               // via a double bond, ipso is nitroso oxygen
                            // RDKit✔️✔️:               isNitrosoO = (bond->getBondType() == Bond::DOUBLE);
                            is_nitroso_o = nbr_bond.order() == BondOrder::Double;
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to nitrogen via a single bond
                            // RDKit✔️✔️:               // and there are no other bonded oxygens
                            // RDKit✔️✔️:               if ((bond->getBondType() == Bond::SINGLE) &&
                            // RDKit✔️✔️:                   (nObondedToCorNorS == 1)) {
                            if nbr_bond.order() == BondOrder::Single
                                && n_o_bonded_to_c_or_n_or_s == 1
                            {
                                let nbr_total_bond_order = mmff_total_bond_order(
                                    explicit_valence,
                                    implicit_hydrogens,
                                    neighbor.atom_index,
                                )?;
                                let nbr_total_degree = mmff_total_degree(
                                    mol,
                                    implicit_hydrogens,
                                    neighbor.atom_index,
                                )?;
                                // RDKit✔️✔️:                 // if nitrogen has 2 neighbors or, if the neighbors are 3,
                                // RDKit✔️✔️:                 // the total bond order on nitrogen is 3, ipso is oxide oxygen
                                // RDKit✔️✔️:                 isOxideOBondedToN =
                                // RDKit✔️✔️:                     ((nbrAtom->getTotalDegree() == 2) ||
                                // RDKit✔️✔️:                      ((nbrAtom->getValence(Atom::ValenceType::EXPLICIT) +
                                // RDKit✔️✔️:                        nbrAtom->getNumImplicitHs()) == 3));
                                is_oxide_o_bonded_to_n =
                                    nbr_total_degree == 2 || nbr_total_bond_order == 3;
                                // RDKit✔️✔️:                 // if the total bond order on nitrogen is 4, ipso is N-oxide
                                // RDKit✔️✔️:                 // oxygen
                                // RDKit✔️✔️:                 isNOxideO = ((nbrAtom->getValence(Atom::ValenceType::EXPLICIT) +
                                // RDKit✔️✔️:                               nbrAtom->getNumImplicitHs()) == 4);
                                is_n_oxide_o = nbr_total_bond_order == 4;
                            }
                            // RDKit✔️✔️:               }
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to nitrogen which is bonded
                            // RDKit✔️✔️:               // to multiple oxygens, ipso is nitro/nitrate oxygen
                            // RDKit✔️✔️:               isNitroO = (nObondedToCorNorS >= 2);
                            is_nitro_o = n_o_bonded_to_c_or_n_or_s >= 2;
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             // if ipso neighbor is sulfur
                        // RDKit✔️✔️:             if (nbrAtom->getAtomicNum() == 16) {
                        if nbr_atom.atomic_number() == 16 {
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to sulfur and
                            // RDKit✔️✔️:               // the latter is bonded to another sulfur,
                            // RDKit✔️✔️:               // ipso is thiosulfinate oxygen
                            // RDKit✔️✔️:               isThioSulfinateO = (nSbondedToCorNorS == 1);
                            is_thio_sulfinate_o = n_s_bonded_to_c_or_n_or_s == 1;
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to sulfur via a single
                            // RDKit✔️✔️:               // bond or, if the bond is double, there are multiple
                            // RDKit✔️✔️:               // oxygen/nitrogen atoms bonded to that sulfur,
                            // RDKit✔️✔️:               // ipso is sulfate oxygen
                            // RDKit✔️✔️:               isSulfateO = ((bond->getBondType() == Bond::SINGLE) ||
                            // RDKit✔️✔️:                             ((bond->getBondType() == Bond::DOUBLE) &&
                            // RDKit✔️✔️:                              ((nObondedToCorNorS + nNbondedToCorNorS) > 1)));
                            is_sulfate_o = nbr_bond.order() == BondOrder::Single
                                || (nbr_bond.order() == BondOrder::Double
                                    && n_o_bonded_to_c_or_n_or_s + n_n_bonded_to_c_or_n_or_s > 1);
                            // RDKit✔️✔️:               // if ipso oxygen is bonded to sulfur via a double
                            // RDKit✔️✔️:               // bond and the sum of oxygen/nitrogen atoms bonded
                            // RDKit✔️✔️:               // to that sulfur is 1, ipso is sulfoxide oxygen
                            // RDKit✔️✔️:               isSulfoxideO = ((bond->getBondType() == Bond::DOUBLE) &&
                            // RDKit✔️✔️:                               ((nObondedToCorNorS + nNbondedToCorNorS) == 1));
                            is_sulfoxide_o = nbr_bond.order() == BondOrder::Double
                                && n_o_bonded_to_c_or_n_or_s + n_n_bonded_to_c_or_n_or_s == 1;
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             // if ipso neighbor is phosphorus or chlorine
                        // RDKit✔️✔️:             isPhosphateOrPerchlorateO = ((nbrAtom->getAtomicNum() == 15) ||
                        // RDKit✔️✔️:                                          (nbrAtom->getAtomicNum() == 17));
                        is_phosphate_or_perchlorate_o = matches!(nbr_atom.atomic_number(), 15 | 17);
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           if (isOxideOBondedToC || isOxideOBondedToN || isOxideOBondedToH) {
                    if is_oxide_o_bonded_to_c || is_oxide_o_bonded_to_n || is_oxide_o_bonded_to_h {
                        // RDKit✔️✔️:             // OM
                        // RDKit✔️✔️:             // Oxide oxygen on sp3 carbon
                        // RDKit✔️✔️:             // OM2
                        // RDKit✔️✔️:             // Oxide oxygen on sp2 carbon
                        // RDKit✔️✔️:             // OM
                        // RDKit✔️✔️:             // Oxide oxygen on sp3 nitrogen (not in original MMFF.I Table III)
                        // RDKit✔️✔️:             // OM2
                        // RDKit✔️✔️:             // Oxide oxygen on sp2 nitrogen (not in original MMFF.I Table III)
                        // RDKit✔️✔️:             atomType = 35;
                        atom_type = 35;
                        // RDKit✔️✔️:             break;
                        // RDKit✔️✔️:           }
                    } else if is_carboxylate_o
                        || is_nitro_o
                        || is_n_oxide_o
                        || is_thio_sulfinate_o
                        || is_sulfate_o
                        || is_phosphate_or_perchlorate_o
                    {
                        // RDKit✔️✔️:           if (isCarboxylateO || isNitroO || isNOxideO || isThioSulfinateO ||
                        // RDKit✔️✔️:               isSulfateO || isPhosphateOrPerchlorateO) {
                        // RDKit✔️✔️:             // O2CM
                        // RDKit✔️✔️:             // Oxygen in carboxylate group
                        // RDKit✔️✔️:             // ONX
                        // RDKit✔️✔️:             // Oxygen in N-oxides
                        // RDKit✔️✔️:             // O2N
                        // RDKit✔️✔️:             // Oxygen in nitro group
                        // RDKit✔️✔️:             // O2NO
                        // RDKit✔️✔️:             // Nitro-group oxygen in nitrate
                        // RDKit✔️✔️:             // O3N
                        // RDKit✔️✔️:             // Nitrate anion oxygen
                        // RDKit✔️✔️:             // OSMS
                        // RDKit✔️✔️:             // Terminal oxygen in thiosulfinate anion
                        // RDKit✔️✔️:             // O-S
                        // RDKit✔️✔️:             // Single terminal O on tetracoordinate sulfur
                        // RDKit✔️✔️:             // O2S
                        // RDKit✔️✔️:             // One of 2 terminal O's on sulfur
                        // RDKit✔️✔️:             // O3S
                        // RDKit✔️✔️:             // One of 3 terminal O's on sulfur
                        // RDKit✔️✔️:             // O4S
                        // RDKit✔️✔️:             // Terminal O in sulfate anion
                        // RDKit✔️✔️:             // OP
                        // RDKit✔️✔️:             // Oxygen in phosphine oxide
                        // RDKit✔️✔️:             // O2P
                        // RDKit✔️✔️:             // One of 2 terminal O's on P
                        // RDKit✔️✔️:             // O3P
                        // RDKit✔️✔️:             // One of 3 terminal O's on P
                        // RDKit✔️✔️:             // O4P
                        // RDKit✔️✔️:             // One of 4 terminal O's on P
                        // RDKit✔️✔️:             // O4Cl
                        // RDKit✔️✔️:             // Oxygen in perchlorate anion
                        // RDKit✔️✔️:             atomType = 32;
                        atom_type = 32;
                        // RDKit✔️✔️:             break;
                        // RDKit✔️✔️:           }
                    } else if is_carbonyl_o || is_nitroso_o || is_sulfoxide_o {
                        // RDKit✔️✔️:           if (isCarbonylO || isNitrosoO || isSulfoxideO) {
                        // RDKit✔️✔️:             // O=C
                        // RDKit✔️✔️:             // Generic carbonyl oxygen
                        // RDKit✔️✔️:             // O=CN
                        // RDKit✔️✔️:             // Carbonyl oxygen in amides
                        // RDKit✔️✔️:             // O=CR
                        // RDKit✔️✔️:             // Carbonyl oxygen in aldehydes and ketones
                        // RDKit✔️✔️:             // O=CO
                        // RDKit✔️✔️:             // Carbonyl oxygen in acids and esters
                        // RDKit✔️✔️:             // O=N
                        // RDKit✔️✔️:             // Nitroso oxygen
                        // RDKit✔️✔️:             // O=S
                        // RDKit✔️✔️:             // Doubly bonded sulfoxide oxygen
                        // RDKit✔️✔️:             atomType = 7;
                        atom_type = 7;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:         }
                }
                // RDKit❌❌:         break;
            }
            // RDKit❗✔️:       // Fluorine
            // RDKit❗✔️:       case 9:
            9 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit❗✔️:         // 1 neighbor
                // RDKit❗✔️:         if (atom->getDegree() == 1) {
                if graph_degree == 1 {
                    // RDKit❗✔️:           // F
                    // RDKit❗✔️:           // Fluorine
                    // RDKit❗✔️:           atomType = 11;
                    atom_type = 11;
                    // RDKit❗✔️:           break;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         if (atom->getDegree() == 0) {
                } else if graph_degree == 0 {
                    // RDKit❗✔️:           // F-
                    // RDKit❗✔️:           // Fluoride anion
                    // RDKit❗✔️:           atomType = 89;
                    atom_type = 89;
                    // RDKit❗✔️:           break;
                }
                // RDKit❗✔️:         }
                // RDKit❌❌:         break;
            }
            // RDKit❗✔️:       // Sodium
            // RDKit❗✔️:       case 11:
            11 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit❗✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit❗✔️:           // NA+
                    // RDKit❗✔️:           // Sodium cation
                    // RDKit❗✔️:           atomType = 93;
                    atom_type = 93;
                    // RDKit❗✔️:           break;
                }
                // RDKit❗✔️:         }
                // RDKit❌❌:         break;
            }
            // RDKit✔️✔️:       // Magnesium
            // RDKit✔️✔️:       case 12:
            12 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           // MG+2
                    // RDKit✔️✔️:           // Dipositive magnesium cation
                    // RDKit✔️✔️:           atomType = 99;
                    atom_type = 99;
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       // Silicon
            // RDKit✔️✔️:       case 14:
            14 => {
                // RDKit✔️✔️:         // SI
                // RDKit✔️✔️:         // Silicon
                // RDKit✔️✔️:         atomType = 19;
                atom_type = 19;
                // RDKit✔️✔️:         break;
            }
            // RDKit❗✔️:       // Phosphorus
            // RDKit❗✔️:       case 15:
            15 => {
                let total_degree = mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                // RDKit❗✔️:         if (atom->getTotalDegree() == 4) {
                if total_degree == 4 {
                    // RDKit❗✔️:           // PO4
                    // RDKit❗✔️:           // Phosphate group phosphorus
                    // RDKit❗✔️:           // PO3
                    // RDKit❗✔️:           // Phosphorus with 3 attached oxygens
                    // RDKit❗✔️:           // PO2
                    // RDKit❗✔️:           // Phosphorus with 2 attached oxygens
                    // RDKit❗✔️:           // PO
                    // RDKit❗✔️:           // Phosphine oxide phosphorus
                    // RDKit❗✔️:           // PTET
                    // RDKit❗✔️:           // General tetracoordinate phosphorus
                    // RDKit❗✔️:           atomType = 25;
                    atom_type = 25;
                    // RDKit❗✔️:           break;
                } else if total_degree == 3 {
                    // RDKit✔️✔️:         if (atom->getTotalDegree() == 3) {
                    // RDKit✔️✔️:           // P
                    // RDKit✔️✔️:           // Phosphorus in phosphines
                    // RDKit✔️✔️:           atomType = 26;
                    atom_type = 26;
                    // RDKit✔️✔️:           break;
                } else if total_degree == 2 {
                    // RDKit✔️✔️:         if (atom->getTotalDegree() == 2) {
                    // RDKit✔️✔️:           // -P=C
                    // RDKit✔️✔️:           // Phosphorus doubly bonded to C
                    // RDKit✔️✔️:           atomType = 75;
                    atom_type = 75;
                    // RDKit✔️✔️:           break;
                }
                // RDKit❌❌:         break;
            }
            // RDKit❗✔️:       // Sulfur
            // RDKit❗✔️:       case 16:
            16 => {
                let total_degree = mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit❗✔️:         // 3  or 4 neighbors
                // RDKit❗✔️:         if ((atom->getTotalDegree() == 3) || (atom->getTotalDegree() == 4)) {
                if total_degree == 3 || total_degree == 4 {
                    // RDKit✔️✔️:           unsigned int nOorNbondedToS = 0;
                    let mut n_o_or_n_bonded_to_s = 0_u32;
                    // RDKit✔️✔️:           unsigned int nSbondedToS = 0;
                    let mut n_s_bonded_to_s = 0_u32;
                    // RDKit✔️✔️:           bool isCDoubleBondedToS = false;
                    let mut is_c_double_bonded_to_s = false;
                    // RDKit✔️✔️:           // loop over neighbors
                    // RDKit✔️✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit✔️✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        // RDKit✔️✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                        let nbr_graph_degree =
                            mol.adjacency.neighbors_of(neighbor.atom_index).len();
                        let nbr_total_degree =
                            mmff_total_degree(mol, implicit_hydrogens, neighbor.atom_index)?;
                        // RDKit✔️✔️:             // check if ipso sulfur is double-bonded to carbon
                        // RDKit✔️✔️:             if ((nbrAtom->getAtomicNum() == 6) &&
                        // RDKit✔️✔️:                 ((mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx()))
                        // RDKit✔️✔️:                      ->getBondType() == Bond::DOUBLE)) {
                        if nbr_atom.atomic_number() == 6 && nbr_bond.order() == BondOrder::Double {
                            // RDKit✔️✔️:               isCDoubleBondedToS = true;
                            is_c_double_bonded_to_s = true;
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:             // if the neighbor is terminal oxygen/sulfur
                        // RDKit✔️✔️:             // or secondary nitrogen, increment the respective counter
                        // RDKit✔️✔️:             if (((nbrAtom->getDegree() == 1) &&
                        // RDKit✔️✔️:                  (nbrAtom->getAtomicNum() == 8)) ||
                        // RDKit✔️✔️:                 ((nbrAtom->getTotalDegree() == 2) &&
                        // RDKit✔️✔️:                  (nbrAtom->getAtomicNum() == 7))) {
                        if (nbr_graph_degree == 1 && nbr_atom.atomic_number() == 8)
                            || (nbr_total_degree == 2 && nbr_atom.atomic_number() == 7)
                        {
                            // RDKit✔️✔️:               ++nOorNbondedToS;
                            n_o_or_n_bonded_to_s += 1;
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:             if ((nbrAtom->getDegree() == 1) &&
                        // RDKit✔️✔️:                 (nbrAtom->getAtomicNum() == 16)) {
                        if nbr_graph_degree == 1 && nbr_atom.atomic_number() == 16 {
                            // RDKit✔️✔️:               ++nSbondedToS;
                            n_s_bonded_to_s += 1;
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit❗✔️:           // if ipso sulfur has 3 neighbors and is bonded to
                    // RDKit❗✔️:           // two atoms of oxygen/nitrogen and double-bonded
                    // RDKit❗✔️:           // to carbon, or if it has 4 neighbors
                    // RDKit❗✔️:           if (((atom->getTotalDegree() == 3) && (nOorNbondedToS == 2) &&
                    // RDKit✔️✔️:                (isCDoubleBondedToS)) ||
                    // RDKit✔️✔️:               (atom->getTotalDegree() == 4)) {
                    if (total_degree == 3 && n_o_or_n_bonded_to_s == 2 && is_c_double_bonded_to_s)
                        || total_degree == 4
                    {
                        // RDKit✔️✔️:             // =SO2
                        // RDKit✔️✔️:             // Sulfone sulfur, doubly bonded to carbon
                        // RDKit✔️✔️:             atomType = 18;
                        atom_type = 18;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           // if ipso sulfur is bonded to both oxygen/nitrogen and sulfur
                    // RDKit✔️✔️:           if ((nOorNbondedToS && nSbondedToS) ||
                    // RDKit✔️✔️:               ((nOorNbondedToS == 2) && (!isCDoubleBondedToS))) {
                    if atom_type == 0
                        && ((n_o_or_n_bonded_to_s != 0 && n_s_bonded_to_s != 0)
                            || (n_o_or_n_bonded_to_s == 2 && !is_c_double_bonded_to_s))
                    {
                        // RDKit✔️✔️:             // SSOM
                        // RDKit✔️✔️:             // Tricoordinate sulfur in anionic thiosulfinate group
                        // RDKit✔️✔️:             atomType = 73;
                        atom_type = 73;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           // otherwise ipso sulfur is double bonded to oxygen or nitrogen
                    // RDKit✔️✔️:           // S=O
                    // RDKit✔️✔️:           // Sulfoxide sulfur
                    // RDKit✔️✔️:           // >S=N
                    // RDKit✔️✔️:           // Tricoordinate sulfur doubly bonded to N
                    // RDKit✔️✔️:           atomType = 17;
                    if atom_type == 0 {
                        atom_type = 17;
                    }
                    // RDKit✔️✔️:           break;
                    // RDKit❌❌:         }
                }
                // RDKit❗✔️:         // 2 neighbors
                // RDKit❗✔️:         if (atom->getTotalDegree() == 2) {
                if atom_type == 0 && total_degree == 2 {
                    // RDKit❗✔️:           // loop over neighbors
                    // RDKit❗✔️:           bool isODoubleBondedToS = false;
                    let mut is_o_double_bonded_to_s = false;
                    // RDKit❗✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit❗✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        // RDKit❗✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                        // RDKit❗✔️:             // check if ipso sulfur is double-bonded to oxygen
                        // RDKit❗✔️:             if ((nbrAtom->getAtomicNum() == 8) &&
                        // RDKit❗✔️:                 ((mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx()))
                        // RDKit❗✔️:                      ->getBondType() == Bond::DOUBLE)) {
                        if nbr_atom.atomic_number() == 8 && nbr_bond.order() == BondOrder::Double {
                            // RDKit❗✔️:               isODoubleBondedToS = true;
                            is_o_double_bonded_to_s = true;
                            // RDKit❗✔️:             }
                        }
                    }
                    // RDKit❗✔️:           }
                    // RDKit❗✔️:           // if ipso sulfur is double-bonded to oxygen
                    // RDKit❗✔️:           if (isODoubleBondedToS) {
                    if is_o_double_bonded_to_s {
                        // RDKit❗✔️:             // =S=O
                        // RDKit❗✔️:             // Sulfinyl sulfur, e.g., in C=S=O
                        // RDKit❗✔️:             atomType = 74;
                        atom_type = 74;
                        // RDKit❗✔️:             break;
                    } else {
                        // RDKit❗✔️:           }
                        // RDKit❗✔️:           // otherwise it is a thiol, sulfide or disulfide
                        // RDKit❗✔️:           // S
                        // RDKit❗✔️:           // Thiol, sulfide, or disulfide sulfur
                        // RDKit❗✔️:           atomType = 15;
                        atom_type = 15;
                        // RDKit❗✔️:           break;
                    }
                    // RDKit❗✔️:         }
                }
                // RDKit✔️✔️:         // 1 neighbor
                // RDKit✔️✔️:         if (atom->getDegree() == 1) {
                if atom_type == 0 && graph_degree == 1 {
                    // RDKit✔️✔️:           unsigned int nTermSbondedToNbr = 0;
                    let mut n_term_s_bonded_to_nbr = 0_u32;
                    // RDKit✔️✔️:           bool isCDoubleBondedToS = false;
                    let mut is_c_double_bonded_to_s = false;
                    // RDKit✔️✔️:           // find the neighbor and count how many terminal sulfur
                    // RDKit✔️✔️:           // atoms are there, including ipso
                    // RDKit✔️✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit✔️✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        // RDKit✔️✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        // RDKit✔️✔️:             boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                        // RDKit✔️✔️:             for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                        for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                            // RDKit✔️✔️:               const Atom *nbr2Atom = mol[*nbr2Idx];
                            let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                            // RDKit✔️✔️:               if ((nbr2Atom->getAtomicNum() == 16) &&
                            // RDKit✔️✔️:                   (nbr2Atom->getTotalDegree() == 1)) {
                            if nbr2_atom.atomic_number() == 16
                                && mmff_total_degree(mol, implicit_hydrogens, neighbor2.atom_index)?
                                    == 1
                            {
                                // RDKit✔️✔️:                 ++nTermSbondedToNbr;
                                n_term_s_bonded_to_nbr += 1;
                            }
                            // RDKit✔️✔️:               }
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:             // check if ipso sulfur is double-bonded to carbon
                        // RDKit✔️✔️:             if ((nbrAtom->getAtomicNum() == 6) &&
                        // RDKit✔️✔️:                 ((mol.getBondBetweenAtoms(atom->getIdx(), nbrAtom->getIdx()))
                        // RDKit✔️✔️:                      ->getBondType() == Bond::DOUBLE)) {
                        let nbr_bond = &mol.bonds.as_slice()[neighbor.bond.index()];
                        if nbr_atom.atomic_number() == 6 && nbr_bond.order() == BondOrder::Double {
                            // RDKit✔️✔️:               isCDoubleBondedToS = true;
                            is_c_double_bonded_to_s = true;
                        }
                        // RDKit✔️✔️:             }
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           // if ipso sulfur is double bonded to carbon and the latter
                    // RDKit✔️✔️:           // is not bonded to other terminal sulfur atoms, then it is
                    // RDKit✔️✔️:           // not a dithiocarboxylate, but a thioketone, etc.
                    // RDKit✔️✔️:           if (isCDoubleBondedToS && (nTermSbondedToNbr != 2)) {
                    if is_c_double_bonded_to_s && n_term_s_bonded_to_nbr != 2 {
                        // RDKit✔️✔️:             // S=C
                        // RDKit✔️✔️:             // Sulfur doubly bonded to carbon
                        // RDKit✔️✔️:             atomType = 16;
                        atom_type = 16;
                        // RDKit✔️✔️:             break;
                    } else {
                        // RDKit✔️✔️:           }
                        // RDKit✔️✔️:           // otherwise ipso must be one of these
                        // RDKit✔️✔️:           // S-P
                        // RDKit✔️✔️:           // Terminal sulfur bonded to P
                        // RDKit✔️✔️:           // SM
                        // RDKit✔️✔️:           // Anionic terminal sulfur
                        // RDKit✔️✔️:           // SSMO
                        // RDKit✔️✔️:           // Terminal sulfur in thiosulfinate group
                        // RDKit✔️✔️:           atomType = 72;
                        atom_type = 72;
                        // RDKit✔️✔️:           break;
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit❌❌:         break;
            }
            // RDKit❗✔️:       // Chlorine
            // RDKit❗✔️:       case 17:
            17 => {
                let total_degree = mmff_total_degree(mol, implicit_hydrogens, atom.id().index())?;
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         // 4 neighbors
                // RDKit✔️✔️:         if (atom->getTotalDegree() == 4) {
                if total_degree == 4 {
                    // RDKit✔️✔️:           // loop over neighbors and count the number
                    // RDKit✔️✔️:           // of bonded oxygens
                    // RDKit✔️✔️:           unsigned int nObondedToCl = 0;
                    let mut n_o_bonded_to_cl = 0_u32;
                    // RDKit✔️✔️:           boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
                    // RDKit✔️✔️:           for (; nbrIdx != endNbrs; ++nbrIdx) {
                    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
                        // RDKit✔️✔️:             const Atom *nbrAtom = mol[*nbrIdx];
                        let nbr_atom = &mol.atoms.as_slice()[neighbor.atom_index];
                        // RDKit✔️✔️:             if (nbrAtom->getAtomicNum() == 8) {
                        if nbr_atom.atomic_number() == 8 {
                            // RDKit✔️✔️:               ++nObondedToCl;
                            n_o_bonded_to_cl += 1;
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           // if there are 4 oxygens
                    // RDKit✔️✔️:           if (nObondedToCl == 4) {
                    if n_o_bonded_to_cl == 4 {
                        // RDKit✔️✔️:             // CLO4
                        // RDKit✔️✔️:             // Perchlorate anione chlorine
                        // RDKit✔️✔️:             atomType = 77;
                        atom_type = 77;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:         }
                }
                // RDKit❗✔️:         // 1 neighbor
                // RDKit❗✔️:         if (atom->getTotalDegree() == 1) {
                if atom_type == 0 && total_degree == 1 {
                    // RDKit❗✔️:           // Cl
                    // RDKit❗✔️:           // Chlorine
                    // RDKit❗✔️:           atomType = 12;
                    atom_type = 12;
                    // RDKit❗✔️:           break;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         // 0 neighbors
                    // RDKit❗✔️:         if (atom->getDegree() == 0) {
                } else if graph_degree == 0 {
                    // RDKit❗✔️:           // Cl-
                    // RDKit❗✔️:           // Chloride anion
                    // RDKit❗✔️:           atomType = 90;
                    atom_type = 90;
                    // RDKit❗✔️:           break;
                }
                // RDKit❗✔️:         }
                // RDKit❌❌:         break;
            }
            // RDKit✔️✔️:       // Potassium
            // RDKit✔️✔️:       case 19:
            19 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           // K+
                    // RDKit✔️✔️:           // Potassium cation
                    // RDKit✔️✔️:           atomType = 94;
                    atom_type = 94;
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       // Calcium
            // RDKit✔️✔️:       case 20:
            20 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           // CA+2
                    // RDKit✔️✔️:           // Dipositive calcium cation
                    // RDKit✔️✔️:           atomType = 96;
                    atom_type = 96;
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       // Iron
            // RDKit✔️✔️:       case 26:
            26 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           if (atom->getFormalCharge() == 2) {
                    if atom.formal_charge() == 2 {
                        // RDKit✔️✔️:             // FE+2
                        // RDKit✔️✔️:             // Dipositive iron cation
                        // RDKit✔️✔️:             atomType = 87;
                        atom_type = 87;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           if (atom->getFormalCharge() == 3) {
                    if atom_type == 0 && atom.formal_charge() == 3 {
                        // RDKit✔️✔️:             // FE+3
                        // RDKit✔️✔️:             // Tripositive iron cation
                        // RDKit✔️✔️:             atomType = 88;
                        atom_type = 88;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       // Copper
            // RDKit✔️✔️:       case 29:
            29 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           if (atom->getFormalCharge() == 1) {
                    if atom.formal_charge() == 1 {
                        // RDKit✔️✔️:             // CU+1
                        // RDKit✔️✔️:             // Monopositive copper cation
                        // RDKit✔️✔️:             atomType = 97;
                        atom_type = 97;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           if (atom->getFormalCharge() == 2) {
                    if atom_type == 0 && atom.formal_charge() == 2 {
                        // RDKit✔️✔️:             // CU+2
                        // RDKit✔️✔️:             // Dipositive copper cation
                        // RDKit✔️✔️:             atomType = 98;
                        atom_type = 98;
                        // RDKit✔️✔️:             break;
                    }
                    // RDKit✔️✔️:           }
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       // Zinc
            // RDKit✔️✔️:       case 30:
            30 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit✔️✔️:         if (atom->getDegree() == 0) {
                if graph_degree == 0 {
                    // RDKit✔️✔️:           // ZN+2
                    // RDKit✔️✔️:           // Dipositive zinc cation
                    // RDKit✔️✔️:           atomType = 95;
                    atom_type = 95;
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit❗✔️:       // Bromine
            // RDKit❗✔️:       case 35:
            35 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit❗✔️:         if (atom->getDegree() == 1) {
                if graph_degree == 1 {
                    // RDKit❗✔️:           // Br
                    // RDKit❗✔️:           // Bromine
                    // RDKit❗✔️:           atomType = 13;
                    atom_type = 13;
                    // RDKit❗✔️:           break;
                    // RDKit❗✔️:         }
                    // RDKit❗✔️:         if (atom->getDegree() == 0) {
                } else if graph_degree == 0 {
                    // RDKit❗✔️:           // BR-
                    // RDKit❗✔️:           // Bromide anion
                    // RDKit❗✔️:           atomType = 91;
                    atom_type = 91;
                    // RDKit❗✔️:           break;
                }
                // RDKit❗✔️:         }
                // RDKit❌❌:         break;
            }
            // RDKit❗✔️:       // Iodine
            // RDKit❗✔️:       case 53:
            53 => {
                let graph_degree = mol.adjacency.neighbors_of(atom.id().index()).len();
                // RDKit❗✔️:         if (atom->getDegree() == 1) {
                if graph_degree == 1 {
                    // RDKit❗✔️:           // I
                    // RDKit❗✔️:           // Iodine
                    // RDKit❗✔️:           atomType = 14;
                    atom_type = 14;
                    // RDKit❗✔️:           break;
                }
                // RDKit❗✔️:         }
                // RDKit❌❌:         break;
            }
            _ => {}
        }
    }

    Ok(MmffAtomProperties {
        atom_type,
        formal_charge: 0.0,
        partial_charge: 0.0,
    })
    // RDKit❌❌:   this->setMMFFAtomType(atom->getIdx(), atomType);
    // RDKit❌❌:   this->setMMFFFormalCharge(atom->getIdx(), formalCharge);
    // RDKit❌❌: }
    // END RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::setMMFFHeavyAtomType
}

fn set_mmff_hydrogen_type(
    mol: &TopologyBlock,
    atom_properties: &[MmffAtomProperties],
    atom: &Atom,
) -> Result<MmffAtomProperties, MmffMolPropertiesError> {
    // BEGIN RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::setMMFFHydrogenType (AtomTyper.cpp:2083-2373)
    // RDKit✔️✔️: // finds the MMFF atomType for a hydrogen atom
    // RDKit✔️✔️: void MMFFMolProperties::setMMFFHydrogenType(const Atom *atom) {
    // RDKit✔️✔️:   unsigned int atomType = 0;
    let mut atom_type = 0_u8;
    // RDKit✔️✔️:   bool isHOCCorHOCN = false;
    let mut is_hocc_or_hocn = false;
    // RDKit✔️✔️:   bool isHOCO = false;
    let mut is_hoco = false;
    // RDKit✔️✔️:   bool isHOP = false;
    let mut is_hop = false;
    // RDKit✔️✔️:   bool isHOS = false;
    let mut is_hos = false;
    // RDKit✔️✔️:   const ROMol &mol = atom->getOwningMol();
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx;
    // RDKit✔️✔️:   ROMol::ADJ_ITER endNbrs;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbr2Idx;
    // RDKit✔️✔️:   ROMol::ADJ_ITER end2Nbrs;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbr3Idx;
    // RDKit✔️✔️:   ROMol::ADJ_ITER end3Nbrs;
    //
    // RDKit✔️✔️:   // loop over neighbors (actually there can be only one)
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   for (; nbrIdx != endNbrs; ++nbrIdx) {
    for neighbor in mol.adjacency.neighbors_of(atom.id().index()) {
        // RDKit✔️✔️:     const Atom *nbrAtom = mol[*nbrIdx];
        let nbr_atom = mol.atoms.as_slice().get(neighbor.atom_index).ok_or(
            MmffMolPropertiesError::AtomIndexOutOfRange {
                atom_index: neighbor.atom_index,
                atoms: mol.atoms.len(),
            },
        )?;
        // RDKit✔️✔️:     switch (nbrAtom->getAtomicNum()) {
        match nbr_atom.atomic_number() {
            // RDKit✔️✔️:       // carbon, silicon
            // RDKit✔️✔️:       case 6:
            // RDKit✔️✔️:       case 14:
            // RDKit✔️✔️:         // HC
            // RDKit✔️✔️:         // Hydrogen attached to carbon
            // RDKit✔️✔️:         // HSI
            // RDKit✔️✔️:         // Hydrogen attached to silicon
            // RDKit✔️✔️:         atomType = 5;
            // RDKit✔️✔️:         break;
            6 | 14 => atom_type = 5,
            //
            // RDKit✔️✔️:       // nitrogen
            // RDKit✔️✔️:       case 7:
            7 => {
                // RDKit✔️✔️:         switch (this->getMMFFAtomType(nbrAtom->getIdx())) {
                let nbr_atom_type = atom_properties
                    .get(neighbor.atom_index)
                    .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                        atom_index: neighbor.atom_index,
                        atoms: atom_properties.len(),
                    })?
                    .atom_type;
                match nbr_atom_type {
                    // RDKit✔️✔️:           case 8:
                    // RDKit✔️✔️:           // HNR
                    // RDKit✔️✔️:           // Generic hydrogen on sp3 nitrogen, e.g. in amines
                    // RDKit✔️✔️:           // H3N
                    // RDKit✔️✔️:           // Hydrogen in ammonia
                    // RDKit✔️✔️:           case 39:
                    // RDKit✔️✔️:           // HPYL
                    // RDKit✔️✔️:           // Hydrogen on nitrogen in pyrrole
                    // RDKit✔️✔️:           case 62:
                    // RDKit✔️✔️:           // HNR
                    // RDKit✔️✔️:           // Generic hydrogen on sp3 nitrogen, e.g. in amines
                    // RDKit✔️✔️:           case 67:
                    // RDKit✔️✔️:           case 68:
                    // RDKit✔️✔️:             // HNOX
                    // RDKit✔️✔️:             // Hydrogen on N in a N-oxide
                    // RDKit✔️✔️:             atomType = 23;
                    // RDKit✔️✔️:             break;
                    8 | 39 | 62 | 67 | 68 => atom_type = 23,
                    //
                    // RDKit✔️✔️:           case 34:
                    // RDKit✔️✔️:           // NR+
                    // RDKit✔️✔️:           // Quaternary nitrogen
                    // RDKit✔️✔️:           case 54:
                    // RDKit✔️✔️:           // N+=C
                    // RDKit✔️✔️:           // Iminium nitrogen
                    // RDKit✔️✔️:           // N+=N
                    // RDKit✔️✔️:           // Positively charged nitrogen doubly bonded to N
                    // RDKit✔️✔️:           case 55:
                    // RDKit✔️✔️:           // HNN+
                    // RDKit✔️✔️:           // Hydrogen on amidinium nitrogen
                    // RDKit✔️✔️:           case 56:
                    // RDKit✔️✔️:           // HGD+
                    // RDKit✔️✔️:           // Hydrogen on guanidinium nitrogen
                    // RDKit✔️✔️:           case 58:
                    // RDKit✔️✔️:           // NPD+
                    // RDKit✔️✔️:           // Aromatic nitrogen in pyridinium
                    // RDKit✔️✔️:           case 81:
                    // RDKit✔️✔️:             // HIM+
                    // RDKit✔️✔️:             // Hydrogen on imidazolium nitrogen
                    // RDKit✔️✔️:             atomType = 36;
                    // RDKit✔️✔️:             break;
                    34 | 54 | 55 | 56 | 58 | 81 => atom_type = 36,
                    //
                    // RDKit✔️✔️:           case 9:
                    // RDKit✔️✔️:             // HN=N
                    // RDKit✔️✔️:             // Hydrogen on azo nitrogen
                    // RDKit✔️✔️:             // HN=C
                    // RDKit✔️✔️:             // Hydrogen on imine nitrogen
                    // RDKit✔️✔️:             atomType = 27;
                    // RDKit✔️✔️:             break;
                    9 => atom_type = 27,
                    //
                    // RDKit✔️✔️:           default:
                    // RDKit✔️✔️:             // HNCC
                    // RDKit✔️✔️:             // Hydrogen on enamine nitrogen
                    // RDKit✔️✔️:             // HNCN
                    // RDKit✔️✔️:             // Hydrogen in H-N-C=N moiety
                    // RDKit✔️✔️:             // HNCO
                    // RDKit✔️✔️:             // Hydrogen on amide nitrogen
                    // RDKit✔️✔️:             // HNCS
                    // RDKit✔️✔️:             // Hydrogen on thioamide nitrogen
                    // RDKit✔️✔️:             // HNNC
                    // RDKit✔️✔️:             // Hydrogen in H-N-N=C moiety
                    // RDKit✔️✔️:             // HNNN
                    // RDKit✔️✔️:             // Hydrogen in H-N-N=N moiety
                    // RDKit✔️✔️:             // HNSO
                    // RDKit✔️✔️:             // Hydrogen on NSO, NSO2, or NSO3 nitrogen
                    // RDKit✔️✔️:             // HNC%
                    // RDKit✔️✔️:             // Hydrogen on N triply bonded to C
                    // RDKit✔️✔️:             // HSP2
                    // RDKit✔️✔️:             // Generic hydrogen on sp2 nitrogen
                    // RDKit✔️✔️:             atomType = 28;
                    // RDKit✔️✔️:             break;
                    _ => atom_type = 28,
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            //
            // RDKit✔️✔️:       // oxygen
            // RDKit✔️✔️:       case 8:
            8 => {
                // RDKit✔️✔️:         switch (this->getMMFFAtomType(nbrAtom->getIdx())) {
                let nbr_atom_type = atom_properties
                    .get(neighbor.atom_index)
                    .ok_or(MmffMolPropertiesError::AtomIndexOutOfRange {
                        atom_index: neighbor.atom_index,
                        atoms: atom_properties.len(),
                    })?
                    .atom_type;
                match nbr_atom_type {
                    // RDKit✔️✔️:           case 49:
                    // RDKit✔️✔️:             // HO+
                    // RDKit✔️✔️:             // Hydrogen on oxonium oxygen
                    // RDKit✔️✔️:             atomType = 50;
                    // RDKit✔️✔️:             break;
                    49 => atom_type = 50,
                    //
                    // RDKit✔️✔️:           case 51:
                    // RDKit✔️✔️:             // HO=+
                    // RDKit✔️✔️:             // Hydrogen on oxenium oxygen
                    // RDKit✔️✔️:             atomType = 52;
                    // RDKit✔️✔️:             break;
                    51 => atom_type = 52,
                    //
                    // RDKit✔️✔️:           case 70:
                    // RDKit✔️✔️:             // HOH
                    // RDKit✔️✔️:             // Hydroxyl hydrogen in water
                    // RDKit✔️✔️:             atomType = 31;
                    // RDKit✔️✔️:             break;
                    70 => atom_type = 31,
                    // RDKit✔️✔️:           case 6:
                    // RDKit✔️✔️:             // for hydrogen bonded to atomType 6 oxygen we need to distinguish
                    // RDKit✔️✔️:             // among acidic hydrogens belonging to carboxylic/phospho acids,
                    // RDKit✔️✔️:             // enolic/phenolic/hydroxamic hydrogens and hydrogens whose oxygen
                    // RDKit✔️✔️:             // partner is bonded to sulfur. If none of these is found
                    // RDKit✔️✔️:             // it is either an alcohol or a generic hydroxyl hydrogen
                    // RDKit✔️✔️:             // loop over oxygen neighbors
                    // RDKit✔️✔️:             boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(nbrAtom);
                    // RDKit✔️✔️:             for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
                    6 => {
                        for neighbor2 in mol.adjacency.neighbors_of(neighbor.atom_index) {
                            // RDKit✔️✔️:               const Atom *nbr2Atom = mol[*nbr2Idx];
                            let nbr2_atom = &mol.atoms.as_slice()[neighbor2.atom_index];
                            // RDKit✔️✔️:               // if the neighbor of oxygen is carbon, loop over the carbon
                            // RDKit✔️✔️:               // neighbors
                            // RDKit✔️✔️:               if (nbr2Atom->getAtomicNum() == 6) {
                            if nbr2_atom.atomic_number() == 6 {
                                // RDKit✔️✔️:                 boost::tie(nbr3Idx, end3Nbrs) = mol.getAtomNeighbors(nbr2Atom);
                                // RDKit✔️✔️:                 for (; nbr3Idx != end3Nbrs; ++nbr3Idx) {
                                for neighbor3 in mol.adjacency.neighbors_of(neighbor2.atom_index) {
                                    // RDKit✔️✔️:                   const Atom *nbr3Atom = mol[*nbr3Idx];
                                    let nbr3_atom = &mol.atoms.as_slice()[neighbor3.atom_index];
                                    // RDKit✔️✔️:                   const Bond *bond = mol.getBondBetweenAtoms(
                                    // RDKit✔️✔️:                       nbr2Atom->getIdx(), nbr3Atom->getIdx());
                                    let bond = &mol.bonds.as_slice()[neighbor3.bond.index()];
                                    // RDKit✔️✔️:                   // if the starting oxygen is met, move on
                                    // RDKit✔️✔️:                   if (nbr3Atom->getIdx() == nbrAtom->getIdx()) {
                                    if neighbor3.atom_index == neighbor.atom_index {
                                        // RDKit✔️✔️:                     continue;
                                        continue;
                                    }
                                    // RDKit✔️✔️:                   }
                                    // RDKit✔️✔️:                   // if the carbon neighbor is another carbon or nitrogen
                                    // RDKit✔️✔️:                   // bonded via a double or aromatic bond, ipso is HOCC/HOCN
                                    // RDKit✔️✔️:                   if (((nbr3Atom->getAtomicNum() == 6) ||
                                    // RDKit✔️✔️:                        (nbr3Atom->getAtomicNum() == 7)) &&
                                    // RDKit✔️✔️:                       ((bond->getBondType() == Bond::DOUBLE) ||
                                    // RDKit✔️✔️:                        (bond->getBondType() == Bond::AROMATIC))) {
                                    if (nbr3_atom.atomic_number() == 6
                                        || nbr3_atom.atomic_number() == 7)
                                        && (bond.order() == BondOrder::Double
                                            || bond.order() == BondOrder::Aromatic)
                                    {
                                        // RDKit✔️✔️:                     isHOCCorHOCN = true;
                                        is_hocc_or_hocn = true;
                                    }
                                    // RDKit✔️✔️:                   }
                                    // RDKit✔️✔️:                   // if the carbon neighbor is an oxygen bonded
                                    // RDKit✔️✔️:                   // via a double bond, ipso is HOCO
                                    // RDKit✔️✔️:                   if ((nbr3Atom->getAtomicNum() == 8) &&
                                    // RDKit✔️✔️:                       (bond->getBondType() == Bond::DOUBLE)) {
                                    if nbr3_atom.atomic_number() == 8
                                        && bond.order() == BondOrder::Double
                                    {
                                        // RDKit✔️✔️:                     isHOCO = true;
                                        is_hoco = true;
                                    }
                                    // RDKit✔️✔️:                   }
                                    // RDKit✔️✔️:                 }
                                }
                            }
                            // RDKit✔️✔️:               }
                            // RDKit✔️✔️:               // if the neighbor of oxygen is phosphorus, ipso is HOCO
                            // RDKit✔️✔️:               if (nbr2Atom->getAtomicNum() == 15) {
                            if nbr2_atom.atomic_number() == 15 {
                                // RDKit✔️✔️:                 isHOP = true;
                                is_hop = true;
                            }
                            // RDKit✔️✔️:               }
                            // RDKit✔️✔️:               // if the neighbor of oxygen is sulfur, ipso is HOS
                            // RDKit✔️✔️:               if (nbr2Atom->getAtomicNum() == 16) {
                            if nbr2_atom.atomic_number() == 16 {
                                // RDKit✔️✔️:                 isHOS = true;
                                is_hos = true;
                            }
                            // RDKit✔️✔️:               }
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:             if (isHOCO || isHOP) {
                        if is_hoco || is_hop {
                            // RDKit✔️✔️:               // HOCO
                            // RDKit✔️✔️:               // Hydroxyl hydrogen in carboxylic acids
                            // RDKit✔️✔️:               atomType = 24;
                            atom_type = 24;
                            // RDKit✔️✔️:               break;
                            // RDKit✔️✔️:             }
                            // RDKit✔️✔️:             if (isHOCCorHOCN) {
                        } else if is_hocc_or_hocn {
                            // RDKit✔️✔️:               // HOCC
                            // RDKit✔️✔️:               // Enolic or phenolic hydroxyl hydrogen
                            // RDKit✔️✔️:               // HOCN
                            // RDKit✔️✔️:               // Hydroxyl hydrogen in HO-C=N moiety
                            // RDKit✔️✔️:               atomType = 29;
                            atom_type = 29;
                            // RDKit✔️✔️:               break;
                            // RDKit✔️✔️:             }
                            // RDKit✔️✔️:             if (isHOS) {
                        } else if is_hos {
                            // RDKit✔️✔️:               // HOS
                            // RDKit✔️✔️:               // Hydrogen on oxygen attached to sulfur
                            // RDKit✔️✔️:               atomType = 33;
                            atom_type = 33;
                            // RDKit✔️✔️:               break;
                            // RDKit✔️✔️:             }
                            // RDKit✔️✔️:             /* FALLTHRU */
                        } else {
                            atom_type = 21;
                        }
                    }
                    // RDKit✔️✔️:           default:
                    // RDKit✔️✔️:             // HO
                    // RDKit✔️✔️:             // Generic hydroxyl hydrogen
                    // RDKit✔️✔️:             // HOR
                    // RDKit✔️✔️:             // Hydroxyl hydrogen in alcohols
                    // RDKit✔️✔️:             atomType = 21;
                    // RDKit✔️✔️:             break;
                    _ => atom_type = 21,
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       // phosphorus and sulfur
            // RDKit✔️✔️:       case 15:
            // RDKit✔️✔️:       case 16:
            // RDKit✔️✔️:         // HP
            // RDKit✔️✔️:         // Hydrogen attached to phosphorus
            // RDKit✔️✔️:         // HS
            // RDKit✔️✔️:         // Hydrogen attached to sulfur
            // RDKit✔️✔️:         // HS=N
            // RDKit✔️✔️:         // Hydrogen attached to >S= sulfur doubly bonded to N
            // RDKit✔️✔️:         atomType = 71;
            // RDKit✔️✔️:         break;
            15 | 16 => atom_type = 71,
            // RDKit✔️✔️:     }
            _ => {}
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_MMFFAtomPropertiesPtrVect[atom->getIdx()]->mmffAtomType = atomType;
    // RDKit✔️✔️:   if (!atomType) {
    // RDKit✔️✔️:     d_valid = false;
    // RDKit✔️✔️:   }
    Ok(MmffAtomProperties {
        atom_type,
        formal_charge: 0.0,
        partial_charge: 0.0,
    })
    // RDKit✔️✔️: }
    // END RDKIT CPP METHOD RDKit::MMFF::MMFFMolProperties::setMMFFHydrogenType
}

fn mmff_total_degree(
    mol: &TopologyBlock,
    implicit_hydrogens: &[i32],
    atom_index: usize,
) -> Result<i32, MmffMolPropertiesError> {
    let atom = mol.atoms.as_slice().get(atom_index).ok_or(
        MmffMolPropertiesError::AtomIndexOutOfRange {
            atom_index,
            atoms: mol.atoms.len(),
        },
    )?;
    let implicit_hydrogens = implicit_hydrogens.get(atom_index).copied().ok_or(
        MmffMolPropertiesError::AtomIndexOutOfRange {
            atom_index,
            atoms: implicit_hydrogens.len(),
        },
    )?;
    Ok(mol.adjacency.neighbors_of(atom_index).len() as i32
        + i32::from(atom.explicit_hydrogens())
        + implicit_hydrogens)
}

fn mmff_total_bond_order(
    explicit_valence: &[i32],
    implicit_hydrogens: &[i32],
    atom_index: usize,
) -> Result<i32, MmffMolPropertiesError> {
    let explicit_valence = explicit_valence.get(atom_index).copied().ok_or(
        MmffMolPropertiesError::AtomIndexOutOfRange {
            atom_index,
            atoms: explicit_valence.len(),
        },
    )?;
    let implicit_hydrogens = implicit_hydrogens.get(atom_index).copied().ok_or(
        MmffMolPropertiesError::AtomIndexOutOfRange {
            atom_index,
            atoms: implicit_hydrogens.len(),
        },
    )?;
    Ok(explicit_valence + implicit_hydrogens)
}

fn mmff_periodic_table_row(atomic_num: u8) -> u32 {
    // BEGIN RDKIT CPP HELPER RDKit::MMFF::getPeriodicTableRow (AtomTyper.cpp:251-267)
    // RDKit✔️✔️: // given the atomic num, this function returns the periodic
    // RDKit✔️✔️: // table row number, starting from 0 for hydrogen
    // RDKit✔️✔️: unsigned int getPeriodicTableRow(const int atomicNum) {
    // RDKit✔️✔️:   unsigned int periodicTableRow = 0;
    let mut periodic_table_row = 0_u32;

    // RDKit✔️✔️:   if ((atomicNum >= 3) && (atomicNum <= 10)) {
    // RDKit✔️✔️:     periodicTableRow = 1;
    if (3..=10).contains(&atomic_num) {
        periodic_table_row = 1;
    // RDKit✔️✔️:   } else if ((atomicNum >= 11) && (atomicNum <= 18)) {
    // RDKit✔️✔️:     periodicTableRow = 2;
    } else if (11..=18).contains(&atomic_num) {
        periodic_table_row = 2;
    // RDKit✔️✔️:   } else if ((atomicNum >= 19) && (atomicNum <= 36)) {
    // RDKit✔️✔️:     periodicTableRow = 3;
    } else if (19..=36).contains(&atomic_num) {
        periodic_table_row = 3;
    // RDKit✔️✔️:   } else if ((atomicNum >= 37) && (atomicNum <= 54)) {
    // RDKit✔️✔️:     periodicTableRow = 4;
    } else if (37..=54).contains(&atomic_num) {
        periodic_table_row = 4;
    }

    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return periodicTableRow;
    // RDKit✔️✔️: }
    periodic_table_row
}

fn mmff_is_aromatic_atom_type(atom_type: u8) -> bool {
    // BEGIN RDKIT CPP DATA ForceFields::MMFF::defaultMMFFArom (Params.cpp:30-31)
    // RDKit❗✔️: const std::vector<std::uint8_t> defaultMMFFArom = {
    // RDKit❗✔️:     37, 38, 39, 44, 58, 59, 63, 64, 65, 66, 69, 76, 78, 79, 80, 81, 82};
    // BEGIN RDKIT CPP METHOD ForceFields::MMFF::MMFFAromCollection::isMMFFAromatic (Params.h:156-159)
    // RDKit❗✔️:   bool isMMFFAromatic(const unsigned int atomType) const {
    // RDKit❗✔️:     return std::find(d_params.begin(), d_params.end(), atomType) !=
    // RDKit❗✔️:            d_params.end();
    // RDKit❗✔️:   }
    matches!(
        atom_type,
        37 | 38 | 39 | 44 | 58 | 59 | 63 | 64 | 65 | 66 | 69 | 76 | 78 | 79 | 80 | 81 | 82
    )
}

fn mmff_periodic_table_row_hl(atomic_num: u8) -> u32 {
    // BEGIN RDKIT CPP HELPER RDKit::MMFF::getPeriodicTableRowHL (AtomTyper.cpp:269-292)
    // RDKit❗✔️: // given the atomic num, this function returns the periodic
    // RDKit❗✔️: // table row number, starting from 1 for helium
    // RDKit❗✔️: // Hydrogen has a special row number (0), while transition
    // RDKit❗✔️: // metals have the row number multiplied by 10
    // RDKit❗✔️: unsigned int getPeriodicTableRowHL(const int atomicNum) {
    // RDKit❗✔️:   unsigned int periodicTableRow = 0;
    let mut periodic_table_row = 0_u32;

    // RDKit❗✔️:   if (atomicNum == 2) {
    // RDKit❗✔️:     periodicTableRow = 1;
    if atomic_num == 2 {
        periodic_table_row = 1;
    // RDKit❗✔️:   } else if ((atomicNum >= 3) && (atomicNum <= 10)) {
    // RDKit❗✔️:     periodicTableRow = 2;
    } else if (3..=10).contains(&atomic_num) {
        periodic_table_row = 2;
    // RDKit❗✔️:   } else if ((atomicNum >= 11) && (atomicNum <= 18)) {
    // RDKit❗✔️:     periodicTableRow = 3;
    } else if (11..=18).contains(&atomic_num) {
        periodic_table_row = 3;
    // RDKit❗✔️:   } else if ((atomicNum >= 19) && (atomicNum <= 36)) {
    // RDKit❗✔️:     periodicTableRow = 4;
    } else if (19..=36).contains(&atomic_num) {
        periodic_table_row = 4;
    // RDKit❗✔️:   } else if ((atomicNum >= 37) && (atomicNum <= 54)) {
    // RDKit❗✔️:     periodicTableRow = 5;
    } else if (37..=54).contains(&atomic_num) {
        periodic_table_row = 5;
    }
    // RDKit❗✔️:   if (((atomicNum >= 21) && (atomicNum <= 30)) ||
    // RDKit❗✔️:       ((atomicNum >= 39) && (atomicNum <= 48))) {
    // RDKit❗✔️:     periodicTableRow *= 10;
    // RDKit❗✔️:   }
    if (21..=30).contains(&atomic_num) || (39..=48).contains(&atomic_num) {
        periodic_table_row *= 10;
    }

    // RDKit❗✔️:   return periodicTableRow;
    // RDKit❗✔️: }
    periodic_table_row
    // END RDKIT CPP HELPER RDKit::MMFF::getPeriodicTableRowHL
}

fn mmff_stretch_bend_type(angle_type: u32, bond_type1: u32, bond_type2: u32) -> u32 {
    // BEGIN RDKIT CPP HELPER RDKit::MMFF::getMMFFStretchBendType (AtomTyper.cpp:2480-2513)
    // RDKit✔️✔️: // given the angle type and the two bond types of the bond
    // RDKit✔️✔️: // which compose the angle, it returns the MMFF stretch-bend
    // RDKit✔️✔️: // type of the angle
    // RDKit✔️✔️: unsigned int getMMFFStretchBendType(const unsigned int angleType,
    // RDKit✔️✔️:                                     const unsigned int bondType1,
    // RDKit✔️✔️:                                     const unsigned int bondType2) {
    // RDKit✔️✔️:   unsigned int stretchBendType = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   switch (angleType) {
    let stretch_bend_type = match angle_type {
        // RDKit✔️✔️:     case 1:
        // RDKit✔️✔️:       stretchBendType = ((bondType1 || (bondType1 == bondType2)) ? 1 : 2);
        // RDKit✔️✔️:       break;
        1 => {
            if bond_type1 != 0 || bond_type1 == bond_type2 {
                1
            } else {
                2
            }
        }
        // RDKit✔️✔️:     case 2:
        // RDKit✔️✔️:       stretchBendType = 3;
        // RDKit✔️✔️:       break;
        2 => 3,
        // RDKit✔️✔️:     case 4:
        // RDKit✔️✔️:       stretchBendType = 4;
        // RDKit✔️✔️:       break;
        4 => 4,
        // RDKit✔️✔️:     case 3:
        // RDKit✔️✔️:       stretchBendType = 5;
        // RDKit✔️✔️:       break;
        3 => 5,
        // RDKit✔️✔️:     case 5:
        // RDKit✔️✔️:       stretchBendType = ((bondType1 || (bondType1 == bondType2)) ? 6 : 7);
        // RDKit✔️✔️:       break;
        5 => {
            if bond_type1 != 0 || bond_type1 == bond_type2 {
                6
            } else {
                7
            }
        }
        // RDKit✔️✔️:     case 6:
        // RDKit✔️✔️:       stretchBendType = 8;
        // RDKit✔️✔️:       break;
        6 => 8,
        // RDKit✔️✔️:     case 7:
        // RDKit✔️✔️:       stretchBendType = ((bondType1 || (bondType1 == bondType2)) ? 9 : 10);
        // RDKit✔️✔️:       break;
        7 => {
            if bond_type1 != 0 || bond_type1 == bond_type2 {
                9
            } else {
                10
            }
        }
        // RDKit✔️✔️:     case 8:
        // RDKit✔️✔️:       stretchBendType = 11;
        // RDKit✔️✔️:       break;
        8 => 11,
        _ => 0,
    };

    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return stretchBendType;
    // RDKit✔️✔️: }
    stretch_bend_type
}

#[cfg(test)]
#[path = "properties_tests.rs"]
mod original_properties_tests;

impl MmffAtomProperties {
    pub const fn atom_type(&self) -> u8 {
        self.atom_type
    }
    pub const fn formal_charge(&self) -> f64 {
        self.formal_charge
    }
    pub const fn partial_charge(&self) -> f64 {
        self.partial_charge
    }
}
