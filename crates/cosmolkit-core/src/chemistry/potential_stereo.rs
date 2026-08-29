// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Source reproduction protocol: dev/source_reproduction_protocol.md

use crate::{
    AtomId, BondId, ChiralTag, Hybridization, KekulizeError, Molecule, StereoError, ValenceError,
    ValenceModel,
};

// RDKit✔️✔️: enum class StereoType {
// RDKit✔️✔️:   Unspecified,
// RDKit✔️✔️:   Atom_Tetrahedral,
// RDKit✔️✔️:   Atom_SquarePlanar,
// RDKit✔️✔️:   Atom_TrigonalBipyramidal,
// RDKit✔️✔️:   Atom_Octahedral,
// RDKit✔️✔️:   Bond_Double,         // single double bond and odd-numbered cumulenes
// RDKit✔️✔️:   Bond_Cumulene_Even,  // even-numbered cumulenes
// RDKit✔️✔️:   Bond_Atropisomer
// RDKit✔️✔️: };
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoType {
    Unspecified,
    AtomTetrahedral,
    AtomSquarePlanar,
    AtomTrigonalBipyramidal,
    AtomOctahedral,
    BondDouble,
    BondEvenCumulene,
    BondAtropisomer,
}

impl StereoType {
    #[must_use]
    pub const fn is_atom_centered(self) -> bool {
        matches!(
            self,
            Self::AtomTetrahedral
                | Self::AtomSquarePlanar
                | Self::AtomTrigonalBipyramidal
                | Self::AtomOctahedral
        )
    }

    #[must_use]
    pub const fn is_bond_centered(self) -> bool {
        matches!(
            self,
            Self::BondDouble | Self::BondEvenCumulene | Self::BondAtropisomer
        )
    }
}

// RDKit✔️✔️: enum class StereoDescriptor {
// RDKit✔️✔️:   None,
// RDKit✔️✔️:   Tet_CW,
// RDKit✔️✔️:   Tet_CCW,
// RDKit✔️✔️:   Bond_Cis,
// RDKit✔️✔️:   Bond_Trans,
// RDKit✔️✔️:   Bond_AtropCW,
// RDKit✔️✔️:   Bond_AtropCCW
// RDKit✔️✔️: };
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoDescriptor {
    None,
    TetrahedralClockwise,
    TetrahedralCounterclockwise,
    BondCis,
    BondTrans,
    BondAtropClockwise,
    BondAtropCounterclockwise,
}

// RDKit✔️✔️: enum class StereoSpecified {
// RDKit✔️✔️:   Unspecified,  // no information provided
// RDKit✔️✔️:   Specified,
// RDKit✔️✔️:   Unknown  // deliberately marked as unknown
// RDKit✔️✔️: };
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoSpecified {
    Unspecified,
    Specified,
    Unknown,
}

/// Typed location of a potential stereochemical element.
///
/// `Missing` represents the source `Atom::NOATOM` default without admitting a
/// magic integer into the Rust model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoCenter {
    Missing,
    Atom(AtomId),
    Bond(BondId),
}

/// One ordered controlling-atom slot in a potential-stereo record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ControllingAtom {
    Missing,
    Atom(AtomId),
}

// RDKit✔️✔️: struct RDKIT_GRAPHMOL_EXPORT StereoInfo {
// RDKit✔️✔️:   StereoType type = StereoType::Unspecified;
// RDKit✔️✔️:   StereoSpecified specified = StereoSpecified::Unspecified;
// RDKit✔️✔️:   unsigned centeredOn = Atom::NOATOM;
// RDKit✔️✔️:   StereoDescriptor descriptor = StereoDescriptor::None;
// RDKit✔️✔️:   unsigned permutation = 0;  // for the non-tetrahedral stereo cases
// RDKit✔️✔️:   std::vector<unsigned> controllingAtoms;  // all atoms around the atom or bond.
// RDKit✔️✔️:   // Order is important
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoInfo {
    stereo_type: StereoType,
    specified: StereoSpecified,
    center: StereoCenter,
    descriptor: StereoDescriptor,
    permutation: u32,
    controlling_atoms: Vec<ControllingAtom>,
}

impl StereoInfo {
    pub(crate) fn new(
        stereo_type: StereoType,
        specified: StereoSpecified,
        center: StereoCenter,
        descriptor: StereoDescriptor,
        permutation: u32,
        controlling_atoms: Vec<ControllingAtom>,
    ) -> Result<Self, PotentialStereoError> {
        let center_matches = match center {
            StereoCenter::Atom(_) => stereo_type.is_atom_centered(),
            StereoCenter::Bond(_) => stereo_type.is_bond_centered(),
            StereoCenter::Missing => stereo_type == StereoType::Unspecified,
        };
        if !center_matches {
            return Err(PotentialStereoError::CenterTypeMismatch {
                stereo_type,
                center,
            });
        }
        Ok(Self {
            stereo_type,
            specified,
            center,
            descriptor,
            permutation,
            controlling_atoms,
        })
    }

    #[must_use]
    pub const fn stereo_type(&self) -> StereoType {
        self.stereo_type
    }

    #[must_use]
    pub const fn specified(&self) -> StereoSpecified {
        self.specified
    }

    #[must_use]
    pub const fn center(&self) -> StereoCenter {
        self.center
    }

    #[must_use]
    pub const fn descriptor(&self) -> StereoDescriptor {
        self.descriptor
    }

    #[must_use]
    pub const fn permutation(&self) -> u32 {
        self.permutation
    }

    #[must_use]
    pub fn controlling_atoms(&self) -> &[ControllingAtom] {
        &self.controlling_atoms
    }

    pub(crate) fn validate_indices(&self, molecule: &Molecule) -> Result<(), PotentialStereoError> {
        match self.center {
            StereoCenter::Missing => {}
            StereoCenter::Atom(atom) if atom.index() >= molecule.num_atoms() => {
                return Err(PotentialStereoError::AtomIndexOutOfBounds { atom });
            }
            StereoCenter::Bond(bond) if bond.index() >= molecule.num_bonds() => {
                return Err(PotentialStereoError::BondIndexOutOfBounds { bond });
            }
            StereoCenter::Atom(_) | StereoCenter::Bond(_) => {}
        }
        for controlling_atom in &self.controlling_atoms {
            if let ControllingAtom::Atom(atom) = controlling_atom
                && atom.index() >= molecule.num_atoms()
            {
                return Err(PotentialStereoError::AtomIndexOutOfBounds { atom: *atom });
            }
        }
        Ok(())
    }
}

/// Options for the single potential-stereo analysis engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PotentialStereoOptions {
    pub clean: bool,
    pub flag_possible: bool,
}

impl Default for PotentialStereoOptions {
    fn default() -> Self {
        // RDKit✔️✔️: python::def(
        // RDKit✔️✔️:     "FindPotentialStereo",
        // RDKit✔️✔️:     (std::vector<Chirality::StereoInfo>(*)(ROMol &, bool, bool)) &
        // RDKit✔️✔️:         Chirality::findPotentialStereo,
        // RDKit✔️✔️:     (python::arg("mol"), python::arg("cleanIt") = false,
        // RDKit✔️✔️:      python::arg("flagPossible") = true),
        Self {
            clean: false,
            flag_possible: true,
        }
    }
}

/// Value-style result of potential-stereo analysis.
#[derive(Debug, Clone)]
pub struct PotentialStereoAnalysis {
    pub molecule: Molecule,
    pub stereo_info: Vec<StereoInfo>,
}

impl PotentialStereoAnalysis {
    pub(crate) fn new(molecule: Molecule, stereo_info: Vec<StereoInfo>) -> Self {
        Self {
            molecule,
            stereo_info,
        }
    }
}

/// Analyze potential stereochemistry on an isolated molecule value.
///
/// The input molecule is never mutated. The returned molecule owns every
/// source-defined cleanup and computed-state change requested by `options`,
/// while `stereo_info` preserves the source record order.
pub fn analyze_potential_stereo(
    molecule: &Molecule,
    options: PotentialStereoOptions,
) -> Result<PotentialStereoAnalysis, PotentialStereoError> {
    let mut workspace = molecule.clone();
    let stereo_info =
        find_potential_stereo_in_workspace(&mut workspace, options.clean, options.flag_possible)?;
    Ok(PotentialStereoAnalysis::new(workspace, stereo_info))
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PotentialStereoError {
    #[error(transparent)]
    Stereo(#[from] StereoError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    CanonicalRanking(#[from] KekulizeError),
    #[error(transparent)]
    RingFinding(#[from] crate::RingFindingError),
    #[error("potential-stereo atom index {atom} is outside the molecule")]
    AtomIndexOutOfBounds { atom: AtomId },
    #[error("potential-stereo bond index {bond} is outside the molecule")]
    BondIndexOutOfBounds { bond: BondId },
    #[error("invalid atom degree in potential-stereo analysis for atom {atom}")]
    InvalidAtomDegree { atom: AtomId },
    #[error("invalid bond endpoint degree in potential-stereo analysis for bond {bond}")]
    InvalidBondDegree { bond: BondId },
    #[error("valence state is unavailable for potential-stereo atom {atom}")]
    MissingValenceState { atom: AtomId },
    #[error("ring state is unavailable for potential-stereo atom {atom}")]
    MissingRingState { atom: AtomId },
    #[error("potential-stereo atom {atom} has invalid implicit hydrogen count {count}")]
    InvalidImplicitHydrogenCount { atom: AtomId, count: i32 },
    #[error("potential-stereo bond {bond} has invalid integer property {property}={value:?}")]
    InvalidBondIntegerProperty {
        bond: BondId,
        property: &'static str,
        value: String,
    },
    #[error("potential-stereo atom {atom} has invalid integer property {property}={value:?}")]
    InvalidAtomIntegerProperty {
        atom: AtomId,
        property: &'static str,
        value: String,
    },
    #[error("ring state is unavailable for potential-stereo bond {bond}")]
    MissingBondRingState { bond: BondId },
    #[error("potential-stereo bond {bond} requires exactly two stereo atoms")]
    InvalidStereoAtomCount { bond: BondId },
    #[error("potential-stereo ring state is unavailable")]
    MissingRingInfo,
    #[error("potential-stereo ring does not contain a bond between atoms {begin} and {end}")]
    MissingRingBond { begin: AtomId, end: AtomId },
    #[error("stereobond duplicate comparison requires two controlling atoms")]
    MissingControllingAtom,
    #[error("potential-stereo atom rank table has no row for atom {atom}")]
    MissingAtomRank { atom: AtomId },
    #[error("potential-stereo state table {table} has length {actual}, expected {expected}")]
    InvalidStateTableLength {
        table: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("stereo atom {stereo_atom} does not control the {end} end of bond {bond}")]
    StereoAtomMismatch {
        bond: BondId,
        stereo_atom: AtomId,
        end: &'static str,
    },
    #[error("stereo type {stereo_type:?} is incompatible with center {center:?}")]
    CenterTypeMismatch {
        stereo_type: StereoType,
        center: StereoCenter,
    },
    #[error("potential-stereo analysis invariant failed: {0}")]
    InvariantViolation(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PotentialStereoInitialization {
    pub(crate) known_atoms: Vec<bool>,
    pub(crate) atom_symbols: Vec<String>,
    pub(crate) possible_atoms: Vec<bool>,
    pub(crate) known_bonds: Vec<bool>,
    pub(crate) bond_symbols: Vec<String>,
    pub(crate) possible_bonds: Vec<bool>,
}

fn atom_compare_symbol(atom: &crate::Atom) -> String {
    // RDKit✔️✔️: inline std::string getAtomCompareSymbol(const Atom &atom) {
    // RDKit✔️✔️:   // we originally tried this with boost::format, but it was WAY slower
    // RDKit✔️✔️:   return std::to_string(atom.getIsotope()) + atom.getSymbol() +
    // RDKit✔️✔️:          std::to_string(atom.getFormalCharge());
    // RDKit✔️✔️: }
    format!(
        "{}{}{}",
        atom.isotope().unwrap_or(0),
        atom.element().symbol(),
        atom.formal_charge()
    )
}

fn bond_compare_symbol(bond: &crate::Bond) -> String {
    // RDKit✔️✔️: std::string getBondSymbol(const Bond *bond) {
    // RDKit✔️✔️:   // FIX: this is not complete
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:   if (bond->getIsAromatic()) {
    // RDKit✔️✔️:     res = ":";
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     switch (bond->getBondType()) {
    // RDKit✔️✔️:       case Bond::BondType::SINGLE:
    // RDKit✔️✔️:         res = "-";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Bond::BondType::DOUBLE:
    // RDKit✔️✔️:         res = "=";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Bond::BondType::TRIPLE:
    // RDKit✔️✔️:         res = "#";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Bond::BondType::AROMATIC:
    // RDKit✔️✔️:         res = ":";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         res = "?";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if bond.is_aromatic() {
        return ":".to_string();
    }
    match bond.order() {
        crate::BondOrder::Single => "-",
        crate::BondOrder::Double => "=",
        crate::BondOrder::Triple => "#",
        crate::BondOrder::Aromatic => ":",
        _ => "?",
    }
    .to_string()
}

fn prepare_non_strict_property_cache(molecule: &mut Molecule) -> Result<(), PotentialStereoError> {
    let atom_count = molecule.num_atoms();
    let needs_update = molecule
        .derived_cache()
        .valence
        .as_ref()
        .is_none_or(|valence| {
            valence.explicit_valence.len() != atom_count
                || valence.implicit_hydrogens.len() != atom_count
        });
    if needs_update {
        // RDKit✔️✔️: if (atom->needsUpdatePropertyCache()) {
        // RDKit✔️✔️:   atom->updatePropertyCache(false);
        // RDKit✔️✔️: }
        // COSMolKit stores the same per-atom values in one typed block, so a
        // missing or incomplete row requires one equivalent non-strict pass.
        let valence = crate::assign_valence_with_options(molecule, ValenceModel::RdkitLike, false)?;
        molecule.derived_cache_mut().valence = Some(valence);
    }
    Ok(())
}

pub(crate) fn initialize_potential_stereo(
    molecule: &mut Molecule,
    flag_possible: bool,
    clean: bool,
    allow_nontetrahedral_stereo: bool,
) -> Result<PotentialStereoInitialization, PotentialStereoError> {
    // RDKit✔️✔️: void initAtomInfo(ROMol &mol, bool flagPossible, bool cleanIt,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> &knownAtoms,
    // RDKit✔️✔️:                   std::vector<std::string> &atomSymbols,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> &possibleAtoms) {
    // RDKit✔️✔️:   bool allowNontetrahedralStereo = getAllowNontetrahedralChirality();
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->needsUpdatePropertyCache()) {
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto aidx = atom->getIdx();
    // RDKit✔️✔️:     atomSymbols[aidx] = getAtomCompareSymbol(*atom);
    // RDKit✔️✔️:     if (detail::isAtomPotentialStereoAtom(atom, allowNontetrahedralStereo)) {
    // RDKit✔️✔️:       auto sinfo = detail::getStereoInfo(atom);
    // RDKit✔️✔️:       switch (sinfo.specified) {
    // RDKit✔️✔️:         case Chirality::StereoSpecified::Unknown:
    // RDKit✔️✔️:           knownAtoms.set(aidx);
    // RDKit✔️✔️:           atomSymbols[aidx] += std::to_string(aidx);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Chirality::StereoSpecified::Specified:
    // RDKit✔️✔️:           knownAtoms.set(aidx);
    // RDKit✔️✔️:           if (sinfo.descriptor == StereoDescriptor::Tet_CCW) {
    // RDKit✔️✔️:             atomSymbols[aidx] += "_CCW";
    // RDKit✔️✔️:           } else if (sinfo.descriptor == StereoDescriptor::Tet_CW) {
    // RDKit✔️✔️:             atomSymbols[aidx] += "_CW";
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             atomSymbols[aidx] += "_STEREO";
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Chirality::StereoSpecified::Unspecified:
    // RDKit✔️✔️:           if (flagPossible) {
    // RDKit✔️✔️:             possibleAtoms.set(aidx);
    // RDKit✔️✔️:             if (!cleanIt) {
    // RDKit✔️✔️:               atomSymbols[aidx] += "_" + std::to_string(aidx);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           throw ValueErrorException("bad StereoInfo.specified type");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (cleanIt) {
    // RDKit✔️✔️:       atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    prepare_non_strict_property_cache(molecule)?;
    let atom_count = molecule.num_atoms();
    let bond_count = molecule.num_bonds();
    let mut state = PotentialStereoInitialization {
        known_atoms: vec![false; atom_count],
        atom_symbols: vec![String::new(); atom_count],
        possible_atoms: vec![false; atom_count],
        known_bonds: vec![false; bond_count],
        bond_symbols: vec![String::new(); bond_count],
        possible_bonds: vec![false; bond_count],
    };
    for atom_index in 0..atom_count {
        let atom_id = AtomId::new(atom_index);
        state.atom_symbols[atom_index] = atom_compare_symbol(&molecule.atoms()[atom_index]);
        if is_atom_potential_stereo(molecule, atom_id, allow_nontetrahedral_stereo)? {
            let info = get_atom_stereo_info(molecule, atom_id, allow_nontetrahedral_stereo)?;
            match info.specified() {
                StereoSpecified::Unknown => {
                    state.known_atoms[atom_index] = true;
                    state.atom_symbols[atom_index].push_str(&atom_index.to_string());
                }
                StereoSpecified::Specified => {
                    state.known_atoms[atom_index] = true;
                    state.atom_symbols[atom_index].push_str(match info.descriptor() {
                        StereoDescriptor::TetrahedralCounterclockwise => "_CCW",
                        StereoDescriptor::TetrahedralClockwise => "_CW",
                        _ => "_STEREO",
                    });
                }
                StereoSpecified::Unspecified if flag_possible => {
                    state.possible_atoms[atom_index] = true;
                    if !clean {
                        state.atom_symbols[atom_index].push('_');
                        state.atom_symbols[atom_index].push_str(&atom_index.to_string());
                    }
                }
                StereoSpecified::Unspecified => {}
            }
        } else if clean {
            molecule.topology_block_mut().atoms[atom_index].set_chiral_tag(ChiralTag::Unspecified);
        }
    }

    // RDKit✔️✔️: void initBondInfo(ROMol &mol, bool flagPossible, bool cleanIt,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> &knownBonds,
    // RDKit✔️✔️:                   std::vector<std::string> &bondSymbols,
    // RDKit✔️✔️:                   boost::dynamic_bitset<> &possibleBonds) {
    // RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto bidx = bond->getIdx();
    // RDKit✔️✔️:     bondSymbols[bidx] = getBondSymbol(bond);
    // RDKit✔️✔️:     if (detail::isBondPotentialStereoBond(bond)) {
    // RDKit✔️✔️:       auto sinfo = detail::getStereoInfo(bond);
    // RDKit✔️✔️:       switch (sinfo.specified) {
    // RDKit✔️✔️:         case Chirality::StereoSpecified::Unknown:
    // RDKit✔️✔️:           knownBonds.set(bidx);
    // RDKit✔️✔️:           bondSymbols[bidx] += "_" + std::to_string(bidx);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Chirality::StereoSpecified::Specified:
    // RDKit✔️✔️:           knownBonds.set(bidx);
    // RDKit✔️✔️:           if (sinfo.descriptor == StereoDescriptor::Bond_Cis) {
    // RDKit✔️✔️:             bondSymbols[bidx] += "_cis";
    // RDKit✔️✔️:           } else if (sinfo.descriptor == StereoDescriptor::Bond_Trans) {
    // RDKit✔️✔️:             bondSymbols[bidx] += "_trans";
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             bondSymbols[bidx] += "_STEREO";
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Chirality::StereoSpecified::Unspecified:
    // RDKit✔️✔️:           if (flagPossible) {
    // RDKit✔️✔️:             possibleBonds.set(bidx);
    // RDKit✔️✔️:             if (!cleanIt) {
    // RDKit✔️✔️:               bondSymbols[bidx] += "_" + std::to_string(bidx);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           throw ValueErrorException("bad StereoInfo.specified type");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       auto currentStereo = bond->getStereo();
    // RDKit✔️✔️:       if (currentStereo != Bond::BondStereo::STEREOATROPCW &&
    // RDKit✔️✔️:           currentStereo != Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:         if (cleanIt) {
    // RDKit✔️✔️:           bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         knownBonds.set(bidx);
    // RDKit✔️✔️:         if (currentStereo == Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️✔️:           bondSymbols[bidx] += "_atropcw";
    // RDKit✔️✔️:         } else if (currentStereo == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:           bondSymbols[bidx] += "_atropccw";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for bond_index in 0..bond_count {
        let bond_id = BondId::new(bond_index);
        state.bond_symbols[bond_index] = bond_compare_symbol(&molecule.bonds()[bond_index]);
        if is_bond_potential_stereo(molecule, bond_id)? {
            let info = get_bond_stereo_info(molecule, bond_id)?;
            match info.specified() {
                StereoSpecified::Unknown => {
                    state.known_bonds[bond_index] = true;
                    state.bond_symbols[bond_index].push('_');
                    state.bond_symbols[bond_index].push_str(&bond_index.to_string());
                }
                StereoSpecified::Specified => {
                    state.known_bonds[bond_index] = true;
                    state.bond_symbols[bond_index].push_str(match info.descriptor() {
                        StereoDescriptor::BondCis => "_cis",
                        StereoDescriptor::BondTrans => "_trans",
                        _ => "_STEREO",
                    });
                }
                StereoSpecified::Unspecified if flag_possible => {
                    state.possible_bonds[bond_index] = true;
                    if !clean {
                        state.bond_symbols[bond_index].push('_');
                        state.bond_symbols[bond_index].push_str(&bond_index.to_string());
                    }
                }
                StereoSpecified::Unspecified => {}
            }
        } else {
            let stereo = molecule.bonds()[bond_index].stereo();
            if matches!(
                stereo,
                crate::BondStereo::AtropCw | crate::BondStereo::AtropCcw
            ) {
                state.known_bonds[bond_index] = true;
                state.bond_symbols[bond_index].push_str(if stereo == crate::BondStereo::AtropCw {
                    "_atropcw"
                } else {
                    "_atropccw"
                });
            } else if clean {
                molecule.topology_block_mut().bonds[bond_index].set_stereo(crate::BondStereo::None);
            }
        }
    }
    Ok(state)
}

fn bond_between_atoms(molecule: &Molecule, begin: AtomId, end: AtomId) -> Option<BondId> {
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(begin.index())
        .iter()
        .find_map(|neighbor| (neighbor.atom_index == end.index()).then_some(neighbor.bond))
}

pub(crate) fn are_stereobond_controlling_atoms_dupes(
    molecule: &Molecule,
    bond: BondId,
    controlling_atom_1: ControllingAtom,
    controlling_atom_2: ControllingAtom,
    atom_ranks: &[u32],
    possible_atoms: &[bool],
    known_atoms: &[bool],
    possible_bonds: &[bool],
    known_bonds: &[bool],
) -> Result<bool, PotentialStereoError> {
    for (table, actual, expected) in [
        ("atom_ranks", atom_ranks.len(), molecule.num_atoms()),
        ("possible_atoms", possible_atoms.len(), molecule.num_atoms()),
        ("known_atoms", known_atoms.len(), molecule.num_atoms()),
        ("possible_bonds", possible_bonds.len(), molecule.num_bonds()),
        ("known_bonds", known_bonds.len(), molecule.num_bonds()),
    ] {
        if actual != expected {
            return Err(PotentialStereoError::InvalidStateTableLength {
                table,
                expected,
                actual,
            });
        }
    }
    // RDKit✔️✔️: bool areStereobondControllingAtomsDupes(
    // RDKit✔️✔️:     const ROMol &mol, const Bond &bond, unsigned controllingAtom1,
    // RDKit✔️✔️:     unsigned controllingAtom2, const std::vector<unsigned> &atomRanks,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &possibleAtoms,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &knownAtoms,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &possibleBonds,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &knownBonds) {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       controllingAtom1 != Atom::NOATOM && controllingAtom2 != Atom::NOATOM,
    // RDKit✔️✔️:       "Missing a controlling atom");
    let ControllingAtom::Atom(controlling_atom_1) = controlling_atom_1 else {
        return Err(PotentialStereoError::MissingControllingAtom);
    };
    let ControllingAtom::Atom(controlling_atom_2) = controlling_atom_2 else {
        return Err(PotentialStereoError::MissingControllingAtom);
    };
    let rank_1 = atom_ranks.get(controlling_atom_1.index()).ok_or(
        PotentialStereoError::MissingAtomRank {
            atom: controlling_atom_1,
        },
    )?;
    let rank_2 = atom_ranks.get(controlling_atom_2.index()).ok_or(
        PotentialStereoError::MissingAtomRank {
            atom: controlling_atom_2,
        },
    )?;

    // RDKit✔️✔️:   if (atomRanks[controllingAtom1] != atomRanks[controllingAtom2]) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if rank_1 != rank_2 {
        return Ok(false);
    }

    // RDKit✔️✔️:   // Now that we know we have 2 neighbors with the same rank, check whether
    // RDKit✔️✔️:   // there is a common even sized ring between the controlling atoms in which
    // RDKit✔️✔️:   // the atom opposite of the bond may be chiral, which would break the tie
    // RDKit✔️✔️:   auto ringInfo = mol.getRingInfo();
    // RDKit✔️✔️:   auto atom1Members = ringInfo->atomMembers(controllingAtom1);
    // RDKit✔️✔️:   auto atom2Members = ringInfo->atomMembers(controllingAtom2);
    let rings = molecule
        .derived_cache()
        .rings
        .as_ref()
        .ok_or(PotentialStereoError::MissingRingInfo)?;
    let atom_1_members = rings.atom_members(controlling_atom_1);
    let atom_2_members = rings.atom_members(controlling_atom_2);

    // RDKit✔️✔️:   auto it1 = atom1Members.begin();
    // RDKit✔️✔️:   auto it2 = atom2Members.begin();
    // RDKit✔️✔️:   while (it1 != atom1Members.end() && it2 != atom2Members.end()) {
    let mut member_1 = 0;
    let mut member_2 = 0;
    while member_1 < atom_1_members.len() && member_2 < atom_2_members.len() {
        // RDKit✔️✔️:     if (*it1 < *it2) {
        // RDKit✔️✔️:       ++it1;
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     } else if (*it1 > *it2) {
        // RDKit✔️✔️:       ++it2;
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if atom_1_members[member_1] < atom_2_members[member_2] {
            member_1 += 1;
            continue;
        } else if atom_1_members[member_1] > atom_2_members[member_2] {
            member_2 += 1;
            continue;
        }

        // RDKit✔️✔️:     auto ring = ringInfo->atomRings().at(*it1);
        // RDKit✔️✔️:     ++it1;
        // RDKit✔️✔️:     ++it2;
        let ring = &rings.atom_rings()[atom_1_members[member_1]];
        member_1 += 1;
        member_2 += 1;

        // RDKit✔️✔️:     if (ring.size() % 2) {
        // RDKit✔️✔️:       // The common ring is odd-sized, so we can't have a tie-breaking atom
        // RDKit✔️✔️:       // directly across the ring, so skip this ring.
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if ring.len() % 2 != 0 {
            continue;
        }

        let bond_state = molecule
            .bonds()
            .get(bond.index())
            .ok_or(PotentialStereoError::BondIndexOutOfBounds { bond })?;
        // RDKit✔️✔️:     for (auto bondEnd : {bond.getBeginAtomIdx(), bond.getEndAtomIdx()}) {
        for bond_end in [bond_state.begin(), bond_state.end()] {
            // RDKit✔️✔️:       auto bondEndPosItr = std::find(ring.begin(), ring.end(), bondEnd);
            // RDKit✔️✔️:       if (bondEndPosItr != ring.end()) {
            if let Some(bond_end_position) = ring.iter().position(|atom| *atom == bond_end) {
                // RDKit✔️✔️:         auto bondEndPos = bondEndPosItr - ring.begin();
                // RDKit✔️✔️:         auto oppositePos = (bondEndPos + ring.size() / 2) % ring.size();
                // RDKit✔️✔️:
                // RDKit✔️✔️:         auto oppositeIdx = ring[oppositePos];
                let opposite = ring[(bond_end_position + ring.len() / 2) % ring.len()];
                // RDKit✔️✔️:         if (possibleAtoms[oppositeIdx] || knownAtoms[oppositeIdx]) {
                // RDKit✔️✔️:           return false;
                // RDKit✔️✔️:         }
                if possible_atoms[opposite.index()] || known_atoms[opposite.index()] {
                    return Ok(false);
                }

                // RDKit✔️✔️:         // atropisomer test: find the bond from the opposite atom that in not in
                // RDKit✔️✔️:         // the ring (if there are three bonds total)
                // RDKit✔️✔️:
                // RDKit✔️✔️:         if (mol.getAtomWithIdx(oppositeIdx)->getDegree() == 3) {
                if molecule
                    .topology_block()
                    .adjacency
                    .neighbors_of(opposite.index())
                    .len()
                    == 3
                {
                    // RDKit✔️✔️:           for (const auto &nbr :
                    // RDKit✔️✔️:                mol.atomBonds(mol.getAtomWithIdx(oppositeIdx))) {
                    for neighbor in molecule
                        .topology_block()
                        .adjacency
                        .neighbors_of(opposite.index())
                    {
                        // RDKit✔️✔️:             auto outOtherAtom = nbr->getOtherAtomIdx(oppositeIdx);
                        // RDKit✔️✔️:             auto bondOutPosItr =
                        // RDKit✔️✔️:                 std::find(ring.begin(), ring.end(), outOtherAtom);
                        // RDKit✔️✔️:             if (bondOutPosItr == ring.end()) {
                        if !ring.iter().any(|atom| atom.index() == neighbor.atom_index) {
                            // RDKit✔️✔️:               if (possibleBonds[nbr->getIdx()] || knownBonds[nbr->getIdx()]) {
                            // RDKit✔️✔️:                 return false;
                            // RDKit✔️✔️:               }
                            if possible_bonds[neighbor.bond.index()]
                                || known_bonds[neighbor.bond.index()]
                            {
                                return Ok(false);
                            }
                        }
                        // RDKit✔️✔️:             }
                        // RDKit✔️✔️:           }
                        // RDKit✔️✔️:         }
                    }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    Ok(true)
}

pub(crate) fn flag_ring_stereo(
    molecule: &mut Molecule,
    possible_ring_stereo_atoms: &mut [u32],
    possible_ring_stereo_bonds: &mut [u32],
    known_atoms: &[bool],
    possible_atoms: Option<&[bool]>,
    known_bonds: &[bool],
    possible_bonds: Option<&[bool]>,
) -> Result<(), PotentialStereoError> {
    for (table, actual, expected) in [
        (
            "possible_ring_stereo_atoms",
            possible_ring_stereo_atoms.len(),
            molecule.num_atoms(),
        ),
        (
            "possible_ring_stereo_bonds",
            possible_ring_stereo_bonds.len(),
            molecule.num_bonds(),
        ),
        ("known_atoms", known_atoms.len(), molecule.num_atoms()),
        ("known_bonds", known_bonds.len(), molecule.num_bonds()),
    ]
    .into_iter()
    .chain(possible_atoms.map(|values| ("possible_atoms", values.len(), molecule.num_atoms())))
    .chain(possible_bonds.map(|values| ("possible_bonds", values.len(), molecule.num_bonds())))
    {
        if actual != expected {
            return Err(PotentialStereoError::InvalidStateTableLength {
                table,
                expected,
                actual,
            });
        }
    }
    // RDKit✔️✔️: void flagRingStereo(ROMol &mol,
    // RDKit✔️✔️:                     std::vector<unsigned int> &possibleRingStereoAtoms,
    // RDKit✔️✔️:                     std::vector<unsigned int> &possibleRingStereoBonds,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> &knownAtoms,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> *possibleAtoms,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> &knownBonds,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> *possibleBonds) {
    // RDKit✔️✔️:   // flag possible ring stereo cases. The relevant cases here are:
    // RDKit✔️✔️:   //    1) even-sized rings with possible (or specified) atoms opposite each
    // RDKit✔️✔️:   //       other, like CC1CC(C)C1 or CC1CCC(C)CC1
    // RDKit✔️✔️:   //    2) atoms sharing a bond which fuses two or more rings, like the
    // RDKit✔️✔️:   //    central
    // RDKit✔️✔️:   //       bond in C1CCC2CCCCC2C1
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto ringInfo = mol.getRingInfo();
    let rings = molecule
        .derived_cache()
        .rings
        .as_ref()
        .ok_or(PotentialStereoError::MissingRingInfo)?;
    // RDKit✔️✔️:   boost::dynamic_bitset<> possibleAtomsInRing(mol.getNumAtoms());
    let mut possible_atoms_in_ring = vec![false; molecule.num_atoms()];
    let mut ring_stereo_other_atom_writes = Vec::new();
    // RDKit✔️✔️:   for (unsigned int ridx = 0; ridx < ringInfo->atomRings().size(); ++ridx) {
    for ring_index in 0..rings.atom_rings().len() {
        // RDKit✔️✔️:     const auto &aring = ringInfo->atomRings()[ridx];
        // RDKit✔️✔️:     const auto &bring = ringInfo->bondRings()[ridx];
        // RDKit✔️✔️:     unsigned int nHere = 0;
        // RDKit✔️✔️:     auto sz = aring.size();
        // RDKit✔️✔️:     bool ringIsOddSized = sz % 2;
        // RDKit✔️✔️:     auto halfSize = sz / 2 + ringIsOddSized;
        let atom_ring = &rings.atom_rings()[ring_index];
        let bond_ring = &rings.bond_rings()[ring_index];
        let mut count_here = 0_u32;
        let ring_size = atom_ring.len();
        let ring_is_odd_sized = ring_size % 2 != 0;
        let half_size = ring_size / 2 + usize::from(ring_is_odd_sized);

        // RDKit✔️✔️:     possibleAtomsInRing.reset();
        possible_atoms_in_ring.fill(false);
        // RDKit✔️✔️:     for (unsigned int ai = 0; ai < sz; ++ai) {
        for atom_position in 0..ring_size {
            // RDKit✔️✔️:       auto aidx = aring[ai];
            // RDKit✔️✔️:       if (!knownAtoms[aidx] && (!possibleAtoms || !possibleAtoms->test(aidx))) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            let atom = atom_ring[atom_position];
            if !known_atoms[atom.index()]
                && possible_atoms.is_none_or(|possible| !possible[atom.index()])
            {
                continue;
            }

            // RDKit✔️✔️:       for (unsigned int ringDivisor : {2, 3}) {
            for ring_divisor in [2_usize, 3] {
                // RDKit✔️✔️:         bool ringIsMultipleOfDivisor = ((sz % ringDivisor) == 0);
                // RDKit✔️✔️:         auto incrementSize = sz / ringDivisor;
                let ring_is_multiple_of_divisor = ring_size % ring_divisor == 0;
                let increment_size = ring_size / ring_divisor;

                // RDKit✔️✔️:         if (ringIsMultipleOfDivisor) {
                if ring_is_multiple_of_divisor {
                    // RDKit✔️✔️:           // find the two indices of the atoms 1/3 the way around the ring
                    // RDKit✔️✔️:
                    // RDKit✔️✔️:           unsigned int otherFoundByBondCount = 0;
                    // RDKit✔️✔️:           unsigned int otherFoundByAtomCount = 0;
                    let mut other_found_by_bond_count = 0_u32;
                    let mut other_found_by_atom_count = 0_u32;
                    // RDKit✔️✔️:           for (unsigned int indexIncrement = incrementSize; indexIncrement < sz;
                    // RDKit✔️✔️:                indexIncrement += incrementSize) {
                    let mut index_increment = increment_size;
                    while index_increment < ring_size {
                        // RDKit✔️✔️:             auto otherIdx = aring[(ai + indexIncrement) % sz];
                        // RDKit✔️✔️:             auto otherAtom = mol.getAtomWithIdx(otherIdx);
                        let other = atom_ring[(atom_position + index_increment) % ring_size];

                        // RDKit✔️✔️:             for (auto bond : mol.atomBonds(otherAtom)) {
                        for neighbor in molecule
                            .topology_block()
                            .adjacency
                            .neighbors_of(other.index())
                        {
                            // RDKit✔️✔️:               auto bidx = bond->getIdx();
                            // RDKit✔️✔️:               if ((knownBonds[bidx] ||
                            // RDKit✔️✔️:                    (possibleBonds && possibleBonds->test(bidx))) &&
                            // RDKit✔️✔️:                   std::find(bring.begin(), bring.end(), bidx) == bring.end()) {
                            if (known_bonds[neighbor.bond.index()]
                                || possible_bonds
                                    .is_some_and(|possible| possible[neighbor.bond.index()]))
                                && !bond_ring.contains(&neighbor.bond)
                            {
                                // RDKit✔️✔️:                 otherFoundByBondCount++;
                                // RDKit✔️✔️:                 break;
                                other_found_by_bond_count += 1;
                                break;
                            }
                            // RDKit✔️✔️:               }
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:             if (otherFoundByBondCount == 0) {
                        if other_found_by_bond_count == 0 {
                            // RDKit✔️✔️:               if (knownAtoms[otherIdx] ||
                            // RDKit✔️✔️:                   (possibleAtoms && possibleAtoms->test(otherIdx))) {
                            if known_atoms[other.index()]
                                || possible_atoms.is_some_and(|possible| possible[other.index()])
                            {
                                // RDKit✔️✔️:                 otherFoundByAtomCount++;
                                other_found_by_atom_count += 1;
                                // RDKit✔️✔️:               }
                            }
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                        index_increment += increment_size;
                    }

                    // RDKit✔️✔️:           if (otherFoundByBondCount == ringDivisor - 1 ||
                    // RDKit✔️✔️:               otherFoundByAtomCount == ringDivisor - 1) {
                    if other_found_by_bond_count == (ring_divisor - 1) as u32
                        || other_found_by_atom_count == (ring_divisor - 1) as u32
                    {
                        // RDKit✔️✔️:             nHere += 1 + otherFoundByBondCount;
                        count_here += 1 + other_found_by_bond_count;
                        // RDKit✔️✔️:             for (unsigned int indexIncrement = 0; indexIncrement < sz;
                        // RDKit✔️✔️:                  indexIncrement += incrementSize) {
                        let mut index_increment = 0;
                        while index_increment < ring_size {
                            // RDKit✔️✔️:               possibleAtomsInRing.set(aring[(ai + indexIncrement) % sz]);
                            possible_atoms_in_ring[atom_ring
                                [(atom_position + index_increment) % ring_size]
                                .index()] = true;
                            // RDKit✔️✔️:             }
                            index_increment += increment_size;
                        }
                        // RDKit✔️✔️:             if (ringDivisor == 2) {
                        if ring_divisor == 2 {
                            // RDKit✔️✔️:               mol.getAtomWithIdx(aidx)->setProp(
                            // RDKit✔️✔️:                   common_properties::_ringStereoOtherAtom,
                            // RDKit✔️✔️:                   aring[(ai + incrementSize) % sz]);
                            ring_stereo_other_atom_writes.push((
                                atom,
                                atom_ring[(atom_position + increment_size) % ring_size],
                            ));
                            // RDKit✔️✔️:             }

                            // RDKit✔️✔️:             continue;
                            continue;
                            // RDKit✔️✔️:           }
                        }
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }

            // RDKit✔️✔️:       // if the atom is in more than one ring, explore the common edge to
            // RDKit✔️✔️:       // see if we can find another potentially chiral atom
            // RDKit✔️✔️:       if (ringInfo->numAtomRings(aidx) > 1) {
            if rings.num_atom_rings(atom) > 1 {
                // RDKit✔️✔️:         auto previousOtherIdx = aidx;
                let mut previous_other = atom;
                // RDKit✔️✔️:         for (size_t step = 1; step <= halfSize; ++step) {
                for step in 1..=half_size {
                    // RDKit✔️✔️:           auto otherIdx = aring[(ai + step) % sz];
                    // RDKit✔️✔️:           auto bnd = mol.getBondBetweenAtoms(previousOtherIdx, otherIdx);
                    let other = atom_ring[(atom_position + step) % ring_size];
                    let common_edge = bond_between_atoms(molecule, previous_other, other).ok_or(
                        PotentialStereoError::MissingRingBond {
                            begin: previous_other,
                            end: other,
                        },
                    )?;
                    // RDKit✔️✔️:           if (ringInfo->numBondRings(bnd->getIdx()) < 2) {
                    // RDKit✔️✔️:             // We reached the end of the common edge.
                    // RDKit✔️✔️:             break;
                    // RDKit✔️✔️:           }
                    if rings.num_bond_rings(common_edge) < 2 {
                        break;
                    }
                    // RDKit✔️✔️:           if (knownAtoms[otherIdx] ||
                    // RDKit✔️✔️:               (possibleAtoms && possibleAtoms->test(otherIdx))) {
                    if known_atoms[other.index()]
                        || possible_atoms.is_some_and(|possible| possible[other.index()])
                    {
                        // RDKit✔️✔️:             // We found another chiral atom, no need to keep
                        // RDKit✔️✔️:             // searching.
                        // RDKit✔️✔️:             nHere += 2;
                        // RDKit✔️✔️:             possibleAtomsInRing.set(aidx);
                        // RDKit✔️✔️:             possibleAtomsInRing.set(otherIdx);
                        count_here += 2;
                        possible_atoms_in_ring[atom.index()] = true;
                        possible_atoms_in_ring[other.index()] = true;
                        // RDKit✔️✔️:             mol.getAtomWithIdx(aidx)->setProp(
                        // RDKit✔️✔️:                 common_properties::_ringStereoOtherAtom, otherIdx);
                        // RDKit✔️✔️:             mol.getAtomWithIdx(otherIdx)->setProp(
                        // RDKit✔️✔️:                 common_properties::_ringStereoOtherAtom, aidx);
                        ring_stereo_other_atom_writes.push((atom, other));
                        ring_stereo_other_atom_writes.push((other, atom));
                        // RDKit✔️✔️:             break;
                        break;
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           previousOtherIdx = otherIdx;
                    previous_other = other;
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     // if the ring contains at least two atoms with possible stereo,
        // RDKit✔️✔️:     // then each of those possibleAtoms should be included for ring stereo
        // RDKit✔️✔️:     if (nHere > 1) {
        if count_here > 1 {
            // RDKit✔️✔️:       for (auto aidx : aring) {
            for atom in atom_ring {
                // RDKit✔️✔️:         if (possibleAtomsInRing[aidx]) {
                // RDKit✔️✔️:           ++possibleRingStereoAtoms[aidx];
                // RDKit✔️✔️:         }
                if possible_atoms_in_ring[atom.index()] {
                    possible_ring_stereo_atoms[atom.index()] += 1;
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:       for (auto bidx : bring) {
            // RDKit✔️✔️:         ++possibleRingStereoBonds[bidx];
            // RDKit✔️✔️:       }
            for bond in bond_ring {
                possible_ring_stereo_bonds[bond.index()] += 1;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }

    if !ring_stereo_other_atom_writes.is_empty() {
        let topology = molecule.topology_properties_mut_for_private_workspace();
        for (atom, other) in ring_stereo_other_atom_writes {
            topology.atoms[atom.index()]
                .set_prop("_ringStereoOtherAtom", other.index().to_string());
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn update_atoms(
    molecule: &Molecule,
    atom_ranks: &[u32],
    atom_symbols: &mut [String],
    possible_atoms: &mut [bool],
    known_atoms: &[bool],
    fixed_atoms: &mut [bool],
    possible_ring_stereo_atoms: &mut [u32],
    possible_ring_stereo_bonds: &mut [u32],
    stereo_info: &mut Vec<StereoInfo>,
    allow_nontetrahedral_stereo: bool,
) -> Result<bool, PotentialStereoError> {
    // RDKit✔️✔️: bool updateAtoms(ROMol &mol, const std::vector<unsigned int> &aranks,
    // RDKit✔️✔️:                  std::vector<std::string> &atomSymbols,
    // RDKit✔️✔️:                  boost::dynamic_bitset<> &possibleAtoms,
    // RDKit✔️✔️:                  boost::dynamic_bitset<> &knownAtoms,
    // RDKit✔️✔️:                  boost::dynamic_bitset<> &fixedAtoms,
    // RDKit✔️✔️:                  std::vector<unsigned int> &possibleRingStereoAtoms,
    // RDKit✔️✔️:                  std::vector<unsigned int> &possibleRingStereoBonds,
    // RDKit✔️✔️:                  std::vector<StereoInfo> &sinfos) {
    for (table, actual) in [
        ("atom_ranks", atom_ranks.len()),
        ("atom_symbols", atom_symbols.len()),
        ("possible_atoms", possible_atoms.len()),
        ("known_atoms", known_atoms.len()),
        ("fixed_atoms", fixed_atoms.len()),
        (
            "possible_ring_stereo_atoms",
            possible_ring_stereo_atoms.len(),
        ),
    ] {
        if actual != molecule.num_atoms() {
            return Err(PotentialStereoError::InvalidStateTableLength {
                table,
                expected: molecule.num_atoms(),
                actual,
            });
        }
    }
    if possible_ring_stereo_bonds.len() != molecule.num_bonds() {
        return Err(PotentialStereoError::InvalidStateTableLength {
            table: "possible_ring_stereo_bonds",
            expected: molecule.num_bonds(),
            actual: possible_ring_stereo_bonds.len(),
        });
    }
    let rings = molecule
        .derived_cache()
        .rings
        .as_ref()
        .ok_or(PotentialStereoError::MissingRingInfo)?;

    // RDKit✔️✔️:   bool needAnotherRound = false;
    let mut need_another_round = false;
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    for atom_index in 0..molecule.num_atoms() {
        // RDKit✔️✔️:     auto aidx = atom->getIdx();
        // RDKit✔️✔️:     if (knownAtoms[aidx] || possibleAtoms[aidx]) {
        if known_atoms[atom_index] || possible_atoms[atom_index] {
            let atom = AtomId::new(atom_index);
            // RDKit✔️✔️:       auto sinfo = detail::getStereoInfo(atom);
            let mut info = get_atom_stereo_info(molecule, atom, allow_nontetrahedral_stereo)?;
            // RDKit✔️✔️:       if (fixedAtoms[aidx]) {
            // RDKit✔️✔️:         sinfos.push_back(std::move(sinfo));
            if fixed_atoms[atom_index] {
                stereo_info.push(info);
                // RDKit✔️✔️:       } else {
            } else {
                // RDKit✔️✔️:         std::vector<unsigned int> nbrs;
                // RDKit✔️✔️:         nbrs.reserve(sinfo.controllingAtoms.size());
                // RDKit✔️✔️:         bool haveADupe = false;
                let mut neighbor_ranks = Vec::with_capacity(info.controlling_atoms().len());
                let mut have_duplicate = false;
                // RDKit✔️✔️:         if (sinfo.type == StereoType::Atom_Tetrahedral) {
                if info.stereo_type() == StereoType::AtomTetrahedral {
                    // RDKit✔️✔️:           for (auto nbrIdx : sinfo.controllingAtoms) {
                    for controller in info.controlling_atoms() {
                        let ControllingAtom::Atom(neighbor) = controller else {
                            return Err(PotentialStereoError::MissingControllingAtom);
                        };
                        // RDKit✔️✔️:             auto rnk = aranks[nbrIdx];
                        let rank = atom_ranks[neighbor.index()];
                        // RDKit✔️✔️:             if (std::find(nbrs.begin(), nbrs.end(), rnk) != nbrs.end()) {
                        if neighbor_ranks.contains(&rank) {
                            // RDKit✔️✔️:               // ok, we just hit a duplicate rank. If the atom we're
                            // RDKit✔️✔️:               // concerned about is a candidate for ring stereo and the bond
                            // RDKit✔️✔️:               // to the atom with the duplicate rank is a ring bond
                            // RDKit✔️✔️:               if (possibleRingStereoAtoms[aidx]) {
                            if possible_ring_stereo_atoms[atom_index] != 0 {
                                // RDKit✔️✔️:                 auto bnd = mol.getBondBetweenAtoms(aidx, nbrIdx);
                                // RDKit✔️✔️:                 if (!bnd || !possibleRingStereoBonds[bnd->getIdx()]) {
                                let connecting_bond = bond_between_atoms(molecule, atom, *neighbor);
                                if connecting_bond.is_none_or(|bond| {
                                    possible_ring_stereo_bonds[bond.index()] == 0
                                }) {
                                    // RDKit✔️✔️:                   haveADupe = true;
                                    // RDKit✔️✔️:                   break;
                                    have_duplicate = true;
                                    break;
                                }
                                // RDKit✔️✔️:                 }
                                // RDKit✔️✔️:               } else {
                            } else {
                                // RDKit✔️✔️:                 haveADupe = true;
                                // RDKit✔️✔️:                 break;
                                have_duplicate = true;
                                break;
                                // RDKit✔️✔️:               }
                            }
                            // RDKit✔️✔️:             } else {
                        } else {
                            // RDKit✔️✔️:               nbrs.push_back(rnk);
                            neighbor_ranks.push(rank);
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         if (!haveADupe) {
                if !have_duplicate {
                    // RDKit✔️✔️:           auto acs = atomSymbols[aidx];
                    let mut atom_compare = atom_symbols[atom_index].clone();
                    // RDKit✔️✔️:           if (!possibleAtoms[aidx]) {
                    if !possible_atoms[atom_index] {
                        // RDKit✔️✔️:             auto sortednbrs = nbrs;
                        // RDKit✔️✔️:             std::sort(sortednbrs.begin(), sortednbrs.end());
                        let mut sorted_neighbor_ranks = neighbor_ranks.clone();
                        sorted_neighbor_ranks.sort_unstable();
                        // RDKit✔️✔️:             // FIX: only works for tetrahedral at the moment
                        // RDKit✔️✔️:             if (sinfo.type == Chirality::StereoType::Atom_Tetrahedral) {
                        if info.stereo_type() == StereoType::AtomTetrahedral {
                            // RDKit✔️✔️:               auto nSwaps = countSwapsToInterconvert(nbrs, sortednbrs);
                            let swaps = crate::source_port_helpers::count_swaps_to_interconvert(
                                &neighbor_ranks,
                                &sorted_neighbor_ranks,
                            )
                            .map_err(|error| {
                                PotentialStereoError::InvariantViolation(format!(
                                    "atom rank swap count failed: {error:?}"
                                ))
                            })?;
                            // RDKit✔️✔️:               if (nSwaps % 2 &&
                            // RDKit✔️✔️:                   (sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW ||
                            // RDKit✔️✔️:                    sinfo.descriptor == Chirality::StereoDescriptor::Tet_CW)) {
                            if swaps % 2 != 0
                                && matches!(
                                    info.descriptor(),
                                    StereoDescriptor::TetrahedralCounterclockwise
                                        | StereoDescriptor::TetrahedralClockwise
                                )
                            {
                                // RDKit✔️✔️:                 sinfo.descriptor =
                                // RDKit✔️✔️:                     sinfo.descriptor == Chirality::StereoDescriptor::Tet_CCW
                                // RDKit✔️✔️:                         ? Chirality::StereoDescriptor::Tet_CW
                                // RDKit✔️✔️:                         : Chirality::StereoDescriptor::Tet_CCW;
                                info.descriptor = match info.descriptor {
                                    StereoDescriptor::TetrahedralCounterclockwise => {
                                        StereoDescriptor::TetrahedralClockwise
                                    }
                                    StereoDescriptor::TetrahedralClockwise => {
                                        StereoDescriptor::TetrahedralCounterclockwise
                                    }
                                    descriptor => descriptor,
                                };
                                // RDKit✔️✔️:               }
                            }
                            // RDKit✔️✔️:               if (sinfo.descriptor == Chirality::StereoDescriptor::Tet_CW) {
                            // RDKit✔️✔️:                 acs = getAtomCompareSymbol(*atom) + "_CW";
                            // RDKit✔️✔️:               } else if (sinfo.descriptor ==
                            // RDKit✔️✔️:                          Chirality::StereoDescriptor::Tet_CCW) {
                            // RDKit✔️✔️:                 acs = getAtomCompareSymbol(*atom) + "_CCW";
                            // RDKit✔️✔️:               }
                            atom_compare = match info.descriptor() {
                                StereoDescriptor::TetrahedralClockwise => format!(
                                    "{}_CW",
                                    atom_compare_symbol(&molecule.atoms()[atom_index])
                                ),
                                StereoDescriptor::TetrahedralCounterclockwise => format!(
                                    "{}_CCW",
                                    atom_compare_symbol(&molecule.atoms()[atom_index])
                                ),
                                _ => atom_compare,
                            };
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:             fixedAtoms.set(aidx);
                        fixed_atoms[atom_index] = true;
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           if (atomSymbols[aidx] != acs) {
                    // RDKit✔️✔️:             atomSymbols[aidx] = acs;
                    // RDKit✔️✔️:             needAnotherRound = true;
                    // RDKit✔️✔️:           }
                    if atom_symbols[atom_index] != atom_compare {
                        atom_symbols[atom_index] = atom_compare;
                        need_another_round = true;
                    }
                    // RDKit✔️✔️:           sinfos.push_back(std::move(sinfo));
                    stereo_info.push(info);
                    // RDKit✔️✔️:         } else {
                } else {
                    // RDKit✔️✔️:           // Only do another round if we change anything here
                    // RDKit✔️✔️:           needAnotherRound |= possibleAtoms[aidx];
                    // RDKit✔️✔️:           possibleAtoms[aidx] = 0;
                    // RDKit✔️✔️:           atomSymbols[aidx] = getAtomCompareSymbol(*atom);
                    need_another_round |= possible_atoms[atom_index];
                    possible_atoms[atom_index] = false;
                    atom_symbols[atom_index] = atom_compare_symbol(&molecule.atoms()[atom_index]);

                    // RDKit✔️✔️:           // if this was creating possible ring stereo, update that info now
                    // RDKit✔️✔️:           if (possibleRingStereoAtoms[aidx]) {
                    if possible_ring_stereo_atoms[atom_index] != 0 {
                        // RDKit✔️✔️:             possibleRingStereoAtoms[aidx] = 0;
                        // RDKit✔️✔️:             needAnotherRound = true;
                        possible_ring_stereo_atoms[atom_index] = 0;
                        need_another_round = true;
                        // RDKit✔️✔️:             // we're no longer in any ring with possible ring stereo. Go
                        // RDKit✔️✔️:             // update all the other atoms/bonds in rings that we're in:
                        // RDKit✔️✔️:             for (unsigned int ridx = 0;
                        // RDKit✔️✔️:                  ridx < mol.getRingInfo()->atomRings().size(); ++ridx) {
                        for ring_index in 0..rings.atom_rings().len() {
                            // RDKit✔️✔️:               const auto &aring = mol.getRingInfo()->atomRings()[ridx];
                            // RDKit✔️✔️:               unsigned int nHere = 0;
                            let atom_ring = &rings.atom_rings()[ring_index];
                            let mut count_here = 0_u32;
                            // RDKit✔️✔️:               for (auto raidx : aring) {
                            for ring_atom in atom_ring {
                                // RDKit✔️✔️:                 // Ring stereo changed, so un-fix atoms in this ring so we can
                                // RDKit✔️✔️:                 // recheck them in the next iteration, for the case that they
                                // RDKit✔️✔️:                 // are no longer after the current atom was declared
                                // RDKit✔️✔️:                 // non-chiral
                                // RDKit✔️✔️:                 fixedAtoms[raidx] = false;
                                // RDKit✔️✔️:                 nHere += (possibleRingStereoAtoms[raidx] > 0);
                                fixed_atoms[ring_atom.index()] = false;
                                count_here +=
                                    u32::from(possible_ring_stereo_atoms[ring_atom.index()] > 0);
                                // RDKit✔️✔️:               }
                            }
                            // RDKit✔️✔️:               if (nHere <= 1) {
                            if count_here <= 1 {
                                // RDKit✔️✔️:                 // if the ring can't transmit stereo anymore, update the
                                // RDKit✔️✔️:                 // counts
                                // RDKit✔️✔️:                 if (nHere == 1) {
                                if count_here == 1 {
                                    // RDKit✔️✔️:                   // update the last potential ring stereo atom in the ring,
                                    // RDKit✔️✔️:                   // since it can't have ring stereo alone.
                                    // RDKit✔️✔️:                   for (auto raidx : aring) {
                                    for ring_atom in atom_ring {
                                        // RDKit✔️✔️:                     if (possibleRingStereoAtoms[raidx]) {
                                        // RDKit✔️✔️:                       --possibleRingStereoAtoms[raidx];
                                        // RDKit✔️✔️:                       break;
                                        // RDKit✔️✔️:                     }
                                        if possible_ring_stereo_atoms[ring_atom.index()] != 0 {
                                            possible_ring_stereo_atoms[ring_atom.index()] -= 1;
                                            break;
                                        }
                                        // RDKit✔️✔️:                   }
                                    }
                                    // RDKit✔️✔️:                 }
                                }
                                // RDKit✔️✔️:                 for (auto rbidx : mol.getRingInfo()->bondRings()[ridx]) {
                                for ring_bond in &rings.bond_rings()[ring_index] {
                                    // RDKit✔️✔️:                   if (possibleRingStereoBonds[rbidx]) {
                                    // RDKit✔️✔️:                     --possibleRingStereoBonds[rbidx];
                                    // RDKit✔️✔️:                   }
                                    if possible_ring_stereo_bonds[ring_bond.index()] != 0 {
                                        possible_ring_stereo_bonds[ring_bond.index()] -= 1;
                                    }
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
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return needAnotherRound;
    // RDKit✔️✔️: }
    Ok(need_another_round)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn update_bonds(
    molecule: &Molecule,
    atom_ranks: &[u32],
    bond_symbols: &mut [String],
    possible_atoms: &[bool],
    possible_bonds: &mut [bool],
    known_atoms: &[bool],
    known_bonds: &[bool],
    fixed_bonds: &mut [bool],
    stereo_info: &mut Vec<StereoInfo>,
) -> Result<bool, PotentialStereoError> {
    // RDKit✔️✔️: bool updateBonds(ROMol &mol, const std::vector<unsigned int> &aranks,
    // RDKit✔️✔️:                  std::vector<std::string> &bondSymbols,
    // RDKit✔️✔️:                  const boost::dynamic_bitset<> &possibleAtoms,
    // RDKit✔️✔️:                  boost::dynamic_bitset<> &possibleBonds,
    // RDKit✔️✔️:                  const boost::dynamic_bitset<> &knownAtoms,
    // RDKit✔️✔️:                  const boost::dynamic_bitset<> &knownBonds,
    // RDKit✔️✔️:                  boost::dynamic_bitset<> &fixedBonds,
    // RDKit✔️✔️:                  std::vector<StereoInfo> &sinfos) {
    for (table, actual, expected) in [
        ("atom_ranks", atom_ranks.len(), molecule.num_atoms()),
        ("possible_atoms", possible_atoms.len(), molecule.num_atoms()),
        ("known_atoms", known_atoms.len(), molecule.num_atoms()),
        ("bond_symbols", bond_symbols.len(), molecule.num_bonds()),
        ("possible_bonds", possible_bonds.len(), molecule.num_bonds()),
        ("known_bonds", known_bonds.len(), molecule.num_bonds()),
        ("fixed_bonds", fixed_bonds.len(), molecule.num_bonds()),
    ] {
        if actual != expected {
            return Err(PotentialStereoError::InvalidStateTableLength {
                table,
                expected,
                actual,
            });
        }
    }

    // RDKit✔️✔️:   bool needAnotherRound = false;
    let mut need_another_round = false;
    // RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
    for bond_index in 0..molecule.num_bonds() {
        // RDKit✔️✔️:     auto bidx = bond->getIdx();
        // RDKit✔️✔️:     if (knownBonds[bidx] || possibleBonds[bidx]) {
        if known_bonds[bond_index] || possible_bonds[bond_index] {
            let bond = BondId::new(bond_index);
            // RDKit✔️✔️:       auto sinfo = detail::getStereoInfo(bond);
            let mut info = get_bond_stereo_info(molecule, bond)?;
            // RDKit✔️✔️:       if (sinfo.type == Chirality::StereoType::Unspecified) {
            // RDKit✔️✔️:         continue;  // not a double bond nor an atropisomer bond
            // RDKit✔️✔️:       }
            if info.stereo_type() == StereoType::Unspecified {
                continue;
            }
            // RDKit✔️✔️:
            // RDKit✔️✔️:       ASSERT_INVARIANT(sinfo.controllingAtoms.size() == 4,
            // RDKit✔️✔️:                        "bad controlling atoms size");
            if info.controlling_atoms().len() != 4 {
                return Err(PotentialStereoError::InvariantViolation(format!(
                    "bond {bond} has {} controlling atoms, expected 4",
                    info.controlling_atoms().len()
                )));
            }

            // RDKit✔️✔️:       if ((sinfo.controllingAtoms[0] == Atom::NOATOM &&
            // RDKit✔️✔️:            sinfo.controllingAtoms[1] == Atom::NOATOM) ||
            // RDKit✔️✔️:           (sinfo.controllingAtoms[2] == Atom::NOATOM &&
            // RDKit✔️✔️:            sinfo.controllingAtoms[3] == Atom::NOATOM)) {
            if (info.controlling_atoms()[0] == ControllingAtom::Missing
                && info.controlling_atoms()[1] == ControllingAtom::Missing)
                || (info.controlling_atoms()[2] == ControllingAtom::Missing
                    && info.controlling_atoms()[3] == ControllingAtom::Missing)
            {
                // RDKit✔️✔️:         // we have a bond with no neighbors on one side, which means it must
                // RDKit✔️✔️:         // have a single implicit H on that side. Since the H is implicit,
                // RDKit✔️✔️:         // there is no way to know whether it is cis or trans.
                // RDKit✔️✔️:         ASSERT_INVARIANT(
                // RDKit✔️✔️:             sinfo.specified != StereoSpecified::Specified,
                // RDKit✔️✔️:             "stereo bond without neighbors can only be unspecified");
                if info.specified() == StereoSpecified::Specified {
                    return Err(PotentialStereoError::InvariantViolation(format!(
                        "specified bond {bond} has no controlling neighbors on one side"
                    )));
                }
                // RDKit✔️✔️:         fixedBonds.set(bidx);
                fixed_bonds[bond_index] = true;
                // RDKit✔️✔️:       }
            }

            // RDKit✔️✔️:       if (fixedBonds[bidx]) {
            // RDKit✔️✔️:         sinfos.push_back(std::move(sinfo));
            if fixed_bonds[bond_index] {
                stereo_info.push(info);
                // RDKit✔️✔️:       } else {
            } else {
                // RDKit✔️✔️:         bool haveADupe = false;
                // RDKit✔️✔️:         bool needsSwap = false;
                let mut have_duplicate = false;
                let mut needs_swap = false;
                // RDKit✔️✔️:         if (sinfo.controllingAtoms[0] != Atom::NOATOM &&
                // RDKit✔️✔️:             sinfo.controllingAtoms[1] != Atom::NOATOM) {
                // RDKit✔️✔️:           if (areStereobondControllingAtomsDupes(
                // RDKit✔️✔️:                   mol, *bond, sinfo.controllingAtoms[0],
                // RDKit✔️✔️:                   sinfo.controllingAtoms[1], aranks, possibleAtoms, knownAtoms,
                // RDKit✔️✔️:                   possibleBonds, knownBonds)) {
                // RDKit✔️✔️:             haveADupe = true;
                // RDKit✔️✔️:           } else if (aranks[sinfo.controllingAtoms[0]] <
                // RDKit✔️✔️:                      aranks[sinfo.controllingAtoms[1]]) {
                // RDKit✔️✔️:             std::swap(sinfo.controllingAtoms[0], sinfo.controllingAtoms[1]);
                // RDKit✔️✔️:             needsSwap = !needsSwap;
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:         if (sinfo.controllingAtoms[2] != Atom::NOATOM &&
                // RDKit✔️✔️:             sinfo.controllingAtoms[3] != Atom::NOATOM) {
                // RDKit✔️✔️:           if (areStereobondControllingAtomsDupes(
                // RDKit✔️✔️:                   mol, *bond, sinfo.controllingAtoms[2],
                // RDKit✔️✔️:                   sinfo.controllingAtoms[3], aranks, possibleAtoms, knownAtoms,
                // RDKit✔️✔️:                   possibleBonds, knownBonds)) {
                // RDKit✔️✔️:             haveADupe = true;
                // RDKit✔️✔️:           } else if (aranks[sinfo.controllingAtoms[2]] <
                // RDKit✔️✔️:                      aranks[sinfo.controllingAtoms[3]]) {
                // RDKit✔️✔️:             std::swap(sinfo.controllingAtoms[2], sinfo.controllingAtoms[3]);
                // RDKit✔️✔️:             needsSwap = !needsSwap;
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                for offset in [0_usize, 2] {
                    if info.controlling_atoms()[offset] != ControllingAtom::Missing
                        && info.controlling_atoms()[offset + 1] != ControllingAtom::Missing
                    {
                        if are_stereobond_controlling_atoms_dupes(
                            molecule,
                            bond,
                            info.controlling_atoms()[offset],
                            info.controlling_atoms()[offset + 1],
                            atom_ranks,
                            possible_atoms,
                            known_atoms,
                            possible_bonds,
                            known_bonds,
                        )? {
                            have_duplicate = true;
                        } else {
                            let ControllingAtom::Atom(first) = info.controlling_atoms()[offset]
                            else {
                                unreachable!("missing controller was excluded above")
                            };
                            let ControllingAtom::Atom(second) =
                                info.controlling_atoms()[offset + 1]
                            else {
                                unreachable!("missing controller was excluded above")
                            };
                            if atom_ranks[first.index()] < atom_ranks[second.index()] {
                                info.controlling_atoms.swap(offset, offset + 1);
                                needs_swap = !needs_swap;
                            }
                        }
                    }
                }
                // RDKit✔️✔️:         if (!haveADupe) {
                if !have_duplicate {
                    // RDKit✔️✔️:           if (needsSwap && (sinfo.descriptor == StereoDescriptor::Bond_Cis ||
                    // RDKit✔️✔️:                             sinfo.descriptor == StereoDescriptor::Bond_Trans)) {
                    if needs_swap
                        && matches!(
                            info.descriptor(),
                            StereoDescriptor::BondCis | StereoDescriptor::BondTrans
                        )
                    {
                        // RDKit✔️✔️:             sinfo.descriptor = sinfo.descriptor == StereoDescriptor::Bond_Cis
                        // RDKit✔️✔️:                                    ? StereoDescriptor::Bond_Trans
                        // RDKit✔️✔️:                                    : StereoDescriptor::Bond_Cis;
                        info.descriptor = if info.descriptor == StereoDescriptor::BondCis {
                            StereoDescriptor::BondTrans
                        } else {
                            StereoDescriptor::BondCis
                        };
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           auto gbs = bondSymbols[bidx];
                    let mut bond_compare = bond_symbols[bond_index].clone();
                    // RDKit✔️✔️:           if (sinfo.specified == StereoSpecified::Specified) {
                    if info.specified() == StereoSpecified::Specified {
                        // RDKit✔️✔️:             switch (sinfo.descriptor) {
                        // RDKit✔️✔️:               case StereoDescriptor::Bond_Cis:
                        // RDKit✔️✔️:                 gbs += "_cis";
                        // RDKit✔️✔️:                 break;
                        // RDKit✔️✔️:               case StereoDescriptor::Bond_Trans:
                        // RDKit✔️✔️:                 gbs += "_trans";
                        // RDKit✔️✔️:                 break;
                        // RDKit✔️✔️:               default:
                        // RDKit✔️✔️:                 break;
                        // RDKit✔️✔️:             }
                        match info.descriptor() {
                            StereoDescriptor::BondCis => bond_compare.push_str("_cis"),
                            StereoDescriptor::BondTrans => bond_compare.push_str("_trans"),
                            _ => {}
                        }
                        // RDKit✔️✔️:           } else if (sinfo.specified == StereoSpecified::Unknown) {
                        // RDKit✔️✔️:             gbs += "_unk";
                    } else if info.specified() == StereoSpecified::Unknown {
                        bond_compare.push_str("_unk");
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           if (bondSymbols[bidx] != gbs) {
                    // RDKit✔️✔️:             bondSymbols[bidx] = gbs;
                    // RDKit✔️✔️:             needAnotherRound = true;
                    // RDKit✔️✔️:           }
                    if bond_symbols[bond_index] != bond_compare {
                        bond_symbols[bond_index] = bond_compare;
                        need_another_round = true;
                    }
                    // RDKit✔️✔️:           if (!possibleBonds[bidx]) {
                    // RDKit✔️✔️:             fixedBonds.set(bidx);
                    // RDKit✔️✔️:           }
                    if !possible_bonds[bond_index] {
                        fixed_bonds[bond_index] = true;
                    }
                    // RDKit✔️✔️:           sinfos.push_back(std::move(sinfo));
                    stereo_info.push(info);
                    // RDKit✔️✔️:         } else if (possibleBonds[bidx]) {
                } else if possible_bonds[bond_index] {
                    // RDKit✔️✔️:           possibleBonds[bidx] = 0;
                    // RDKit✔️✔️:           bondSymbols[bidx] = getBondSymbol(bond);
                    // RDKit✔️✔️:           needAnotherRound = true;
                    possible_bonds[bond_index] = false;
                    bond_symbols[bond_index] = bond_compare_symbol(&molecule.bonds()[bond_index]);
                    need_another_round = true;
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return needAnotherRound;
    // RDKit✔️✔️: }
    Ok(need_another_round)
}

fn run_cleanup_fixed_point(
    molecule: &Molecule,
    atom_symbols: &mut [String],
    bond_symbols: &mut [String],
    possible_atoms: &mut [bool],
    possible_bonds: &mut [bool],
    known_atoms: &[bool],
    known_bonds: &[bool],
    fixed_atoms: &mut [bool],
    fixed_bonds: &mut [bool],
    possible_ring_stereo_atoms: &mut [u32],
    possible_ring_stereo_bonds: &mut [u32],
    allow_nontetrahedral_stereo: bool,
) -> Result<(Vec<StereoInfo>, Vec<u32>), PotentialStereoError> {
    use crate::canon_rank::{CanonicalRankOptions, FragmentRankScope, rank_fragment_atoms};

    let atoms_in_play = vec![true; molecule.num_atoms()];
    let bonds_in_play = vec![true; molecule.num_bonds()];
    let rank_options = CanonicalRankOptions {
        break_ties: false,
        include_chirality: false,
        include_isotopes: false,
        include_atom_maps: false,
        include_chiral_presence: false,
        include_stereo_groups: false,
        use_non_stereo_ranks: false,
        include_ring_stereo: false,
        chirality_rings_use_ring_stereo: false,
    };
    let mut result = Vec::new();
    let mut atom_ranks = vec![0_u32; molecule.num_atoms()];
    let mut need_another_round = true;
    while need_another_round {
        // RDKit✔️✔️: res.clear();
        result.clear();
        // RDKit✔️✔️: const bool includeChirality = false;
        // RDKit✔️✔️: const bool includeIsotopes = false;
        // RDKit✔️✔️: const bool breakTies = false;
        // RDKit✔️✔️: const bool includeAtomMaps = false;
        // RDKit✔️✔️: const bool includeChiralPresence = false;
        // RDKit✔️✔️: const bool useRingStereo = false;
        // RDKit✔️✔️: Canon::rankFragmentAtoms(mol, aranks, atomsInPlay, bondsInPlay,
        // RDKit✔️✔️:                          &atomSymbols, &bondSymbols, breakTies,
        // RDKit✔️✔️:                          includeChirality, includeIsotopes, includeAtomMaps,
        // RDKit✔️✔️:                          includeChiralPresence, useRingStereo);
        let scope = FragmentRankScope::new(&atoms_in_play, &bonds_in_play)
            .with_atom_symbols(atom_symbols)
            .with_bond_symbols(bond_symbols);
        atom_ranks = rank_fragment_atoms(molecule, scope, rank_options)?
            .into_iter()
            .map(|rank| {
                u32::try_from(rank).map_err(|_| {
                    PotentialStereoError::InvariantViolation(format!(
                        "canonical atom rank {rank} does not fit the RDKit unsigned rank type"
                    ))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        // RDKit✔️✔️: needAnotherRound = updateAtoms(
        // RDKit✔️✔️:     mol, aranks, atomSymbols, possibleAtoms, knownAtoms, fixedAtoms,
        // RDKit✔️✔️:     possibleRingStereoAtoms, possibleRingStereoBonds, res);
        need_another_round = update_atoms(
            molecule,
            &atom_ranks,
            atom_symbols,
            possible_atoms,
            known_atoms,
            fixed_atoms,
            possible_ring_stereo_atoms,
            possible_ring_stereo_bonds,
            &mut result,
            allow_nontetrahedral_stereo,
        )?;
        // RDKit✔️✔️: needAnotherRound |=
        // RDKit✔️✔️:     updateBonds(mol, aranks, bondSymbols, possibleAtoms, possibleBonds,
        // RDKit✔️✔️:                 knownAtoms, knownBonds, fixedBonds, res);
        need_another_round |= update_bonds(
            molecule,
            &atom_ranks,
            bond_symbols,
            possible_atoms,
            possible_bonds,
            known_atoms,
            known_bonds,
            fixed_bonds,
            &mut result,
        )?;
    }
    Ok((result, atom_ranks))
}

pub(crate) fn run_cleanup(
    molecule: &mut Molecule,
    flag_possible: bool,
    clean: bool,
) -> Result<Vec<StereoInfo>, PotentialStereoError> {
    // RDKit✔️✔️: std::vector<StereoInfo> runCleanup(ROMol &mol, bool flagPossible,
    // RDKit✔️✔️:                                    bool cleanIt) {
    // RDKit✔️✔️:   // This potentially does two passes of "canonicalization" to identify
    // RDKit✔️✔️:   // stereo atoms/bonds:
    // RDKit✔️✔️:   //   - if cleanIt is true we start with a pass which ignores possible
    // RDKit✔️✔️:   //   stereo atoms/bonds and which removes the stereo spec from any
    // RDKit✔️✔️:   //   atom/bond which doesn't have unique neighbors
    // RDKit✔️✔️:   //   - if flagPossible is true we do a pass where each unspecified
    // RDKit✔️✔️:   //   possible stereocenter is treated as if it were different from all
    // RDKit✔️✔️:   //   others. This allows us to identify every possible stereo atom/bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::dynamic_bitset<> knownAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<std::string> atomSymbols(mol.getNumAtoms());
    // RDKit✔️✔️:   boost::dynamic_bitset<> possibleAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   initAtomInfo(mol, flagPossible, cleanIt, knownAtoms, atomSymbols,
    // RDKit✔️✔️:                possibleAtoms);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<std::string> bondSymbols(mol.getNumBonds());
    // RDKit✔️✔️:   boost::dynamic_bitset<> knownBonds(mol.getNumBonds());
    // RDKit✔️✔️:   boost::dynamic_bitset<> possibleBonds(mol.getNumBonds());
    // RDKit✔️✔️:   initBondInfo(mol, flagPossible, cleanIt, knownBonds, bondSymbols,
    // RDKit✔️✔️:                possibleBonds);
    let allow_nontetrahedral_stereo = crate::stereo::get_allow_nontetrahedral_chirality();
    let mut state =
        initialize_potential_stereo(molecule, flag_possible, clean, allow_nontetrahedral_stereo)?;

    // RDKit✔️✔️:   // copy the original sets of possible atoms/bonds. We need them in the
    // RDKit✔️✔️:   // second pass
    // RDKit✔️✔️:   auto origPossibleAtoms = possibleAtoms;
    // RDKit✔️✔️:   auto origPossibleBonds = possibleBonds;
    let original_possible_atoms = state.possible_atoms.clone();
    let original_possible_bonds = state.possible_bonds.clone();

    // RDKit✔️✔️:   // tracks the number of rings with possible ring stereo that the atom is
    // RDKit✔️✔️:   // in (only set for potential stereoatoms)
    // RDKit✔️✔️:   std::vector<unsigned int> possibleRingStereoAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   // tracks the number of rings with possible ring stereo that the bond is
    // RDKit✔️✔️:   // in (set for all bonds)
    // RDKit✔️✔️:   std::vector<unsigned int> possibleRingStereoBonds(mol.getNumBonds());
    let mut possible_ring_stereo_atoms = vec![0_u32; molecule.num_atoms()];
    let mut possible_ring_stereo_bonds = vec![0_u32; molecule.num_bonds()];

    // RDKit✔️✔️:   // identify atoms which can be involved in ring stereo
    // RDKit✔️✔️:   flagRingStereo(mol, possibleRingStereoAtoms, possibleRingStereoBonds,
    // RDKit✔️✔️:                  knownAtoms, cleanIt ? nullptr : &possibleAtoms, knownBonds,
    // RDKit✔️✔️:                  cleanIt ? nullptr : &possibleBonds);
    flag_ring_stereo(
        molecule,
        &mut possible_ring_stereo_atoms,
        &mut possible_ring_stereo_bonds,
        &state.known_atoms,
        (!clean).then_some(state.possible_atoms.as_slice()),
        &state.known_bonds,
        (!clean).then_some(state.possible_bonds.as_slice()),
    )?;

    // RDKit✔️✔️:   // these are used to track which atoms/bonds have been altered
    // RDKit✔️✔️:   boost::dynamic_bitset<> fixedAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   boost::dynamic_bitset<> fixedBonds(mol.getNumBonds());
    let mut fixed_atoms = vec![false; molecule.num_atoms()];
    let mut fixed_bonds = vec![false; molecule.num_bonds()];

    // RDKit✔️✔️:   // used to tell rankFragmentAtoms to use all atoms:
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    // RDKit✔️✔️:   atomsInPlay.set();
    // RDKit✔️✔️:   boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds());
    // RDKit✔️✔️:   bondsInPlay.set();
    // RDKit✔️✔️:   std::vector<unsigned int> aranks(mol.getNumAtoms());
    // RDKit✔️✔️:   bool needAnotherRound = true;
    // RDKit✔️✔️:   while (needAnotherRound) {
    // RDKit✔️✔️:     res.clear();
    // RDKit✔️✔️:   }
    let (mut result, mut atom_ranks) = run_cleanup_fixed_point(
        molecule,
        &mut state.atom_symbols,
        &mut state.bond_symbols,
        &mut state.possible_atoms,
        &mut state.possible_bonds,
        &state.known_atoms,
        &state.known_bonds,
        &mut fixed_atoms,
        &mut fixed_bonds,
        &mut possible_ring_stereo_atoms,
        &mut possible_ring_stereo_bonds,
        allow_nontetrahedral_stereo,
    )?;

    // RDKit✔️✔️:   if (cleanIt) {
    // RDKit✔️✔️:     // remove stereo specs from atoms/bonds which should not have them
    // RDKit✔️✔️:     cleanMolStereo(mol, fixedAtoms, knownAtoms, fixedBonds, knownBonds);
    // RDKit✔️✔️:   }
    if clean {
        clean_mol_stereo(
            molecule,
            &fixed_atoms,
            &state.known_atoms,
            &fixed_bonds,
            &state.known_bonds,
        )?;
    }

    // RDKit✔️✔️:   if (flagPossible && (possibleAtoms != origPossibleAtoms ||
    // RDKit✔️✔️:                        possibleBonds != origPossibleBonds)) {
    // RDKit✔️✔️:     // if we're doing "flagPossible" mode and have done some cleanup, then
    // RDKit✔️✔️:     // we need to do another iteration
    if flag_possible
        && (state.possible_atoms != original_possible_atoms
            || state.possible_bonds != original_possible_bonds)
    {
        // RDKit✔️✔️:     possibleAtoms = origPossibleAtoms;
        // RDKit✔️✔️:     // flag every center/bond where we removed stereo as possible:
        // RDKit✔️✔️:     for (auto i = 0u; i < mol.getNumAtoms(); ++i) {
        // RDKit✔️✔️:       if (!fixedAtoms[i] && knownAtoms[i]) {
        // RDKit✔️✔️:         possibleAtoms[i] = 1;
        // RDKit✔️✔️:         knownAtoms[i] = 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (possibleAtoms[i]) {
        // RDKit✔️✔️:         atomSymbols[i] += "_" + std::to_string(i);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        state.possible_atoms = original_possible_atoms;
        for atom_index in 0..molecule.num_atoms() {
            if !fixed_atoms[atom_index] && state.known_atoms[atom_index] {
                state.possible_atoms[atom_index] = true;
                state.known_atoms[atom_index] = false;
            }
            if state.possible_atoms[atom_index] {
                state.atom_symbols[atom_index].push('_');
                state.atom_symbols[atom_index].push_str(&atom_index.to_string());
            }
        }

        // RDKit✔️✔️:     possibleBonds = origPossibleBonds;
        // RDKit✔️✔️:     for (auto i = 0u; i < mol.getNumBonds(); ++i) {
        // RDKit✔️✔️:       if (!fixedBonds[i] && knownBonds[i]) {
        // RDKit✔️✔️:         possibleBonds[i] = 1;
        // RDKit✔️✔️:         knownBonds[i] = 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (possibleBonds[i]) {
        // RDKit✔️✔️:         bondSymbols[i] += "_" + std::to_string(i);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        state.possible_bonds = original_possible_bonds;
        for bond_index in 0..molecule.num_bonds() {
            if !fixed_bonds[bond_index] && state.known_bonds[bond_index] {
                state.possible_bonds[bond_index] = true;
                state.known_bonds[bond_index] = false;
            }
            if state.possible_bonds[bond_index] {
                state.bond_symbols[bond_index].push('_');
                state.bond_symbols[bond_index].push_str(&bond_index.to_string());
            }
        }

        // RDKit✔️✔️:     flagRingStereo(mol, possibleRingStereoAtoms, possibleRingStereoBonds,
        // RDKit✔️✔️:                    knownAtoms, &possibleAtoms, knownBonds, &possibleBonds);
        flag_ring_stereo(
            molecule,
            &mut possible_ring_stereo_atoms,
            &mut possible_ring_stereo_bonds,
            &state.known_atoms,
            Some(&state.possible_atoms),
            &state.known_bonds,
            Some(&state.possible_bonds),
        )?;

        // RDKit✔️✔️:     needAnotherRound = true;
        // RDKit✔️✔️:     while (needAnotherRound) {
        // RDKit✔️✔️:       res.clear();
        // RDKit✔️✔️:       needAnotherRound = updateAtoms(
        // RDKit✔️✔️:           mol, aranks, atomSymbols, possibleAtoms, knownAtoms, fixedAtoms,
        // RDKit✔️✔️:           possibleRingStereoAtoms, possibleRingStereoBonds, res);
        // RDKit✔️✔️:       needAnotherRound |=
        // RDKit✔️✔️:           updateBonds(mol, aranks, bondSymbols, possibleAtoms, possibleBonds,
        // RDKit✔️✔️:                       knownAtoms, knownBonds, fixedBonds, res);
        // RDKit✔️✔️:     }
        (result, atom_ranks) = run_cleanup_fixed_point(
            molecule,
            &mut state.atom_symbols,
            &mut state.bond_symbols,
            &mut state.possible_atoms,
            &mut state.possible_bonds,
            &state.known_atoms,
            &state.known_bonds,
            &mut fixed_atoms,
            &mut fixed_bonds,
            &mut possible_ring_stereo_atoms,
            &mut possible_ring_stereo_bonds,
            allow_nontetrahedral_stereo,
        )?;
    }

    // RDKit✔️✔️:   boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    // RDKit✔️✔️:   findChiralAtomSpecialCases(mol, possibleSpecialCases, aranks);
    crate::stereo::find_chiral_atom_special_cases(molecule, &atom_ranks)?;

    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     atom->setProp<unsigned int>(common_properties::_ChiralAtomRank,
    // RDKit✔️✔️:                                 aranks[atom->getIdx()], true);
    // RDKit✔️✔️:   }
    let topology = molecule.topology_properties_mut_for_private_workspace();
    for (atom_index, rank) in atom_ranks.into_iter().enumerate() {
        topology.atoms[atom_index].set_computed_prop("_chiralAtomRank", rank.to_string());
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(result)
}

pub(crate) fn find_potential_stereo_in_workspace(
    molecule: &mut Molecule,
    clean: bool,
    flag_possible: bool,
) -> Result<Vec<StereoInfo>, PotentialStereoError> {
    // RDKit✔️✔️: std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
    // RDKit✔️✔️:                                             bool findPossible) {
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol);
    // RDKit✔️✔️:   }
    if !molecule
        .derived_cache()
        .rings
        .as_ref()
        .is_some_and(crate::RingInfo::is_symm_sssr)
    {
        let rings = crate::symmetrize_sssr(molecule)?;
        molecule.derived_cache_mut().rings = Some(rings);
    }

    // RDKit✔️✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    prepare_non_strict_property_cache(molecule)?;
    // RDKit✔️✔️:   std::vector<StereoInfo> res = runCleanup(mol, findPossible, cleanIt);
    let result = run_cleanup(molecule, flag_possible, clean)?;
    // RDKit✔️✔️:   mol.setProp("_potentialStereo", res, true);
    molecule.set_potential_stereo_cache(std::sync::Arc::from(result.clone()));
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(result)
}

pub(crate) fn clean_mol_stereo(
    molecule: &mut Molecule,
    fixed_atoms: &[bool],
    known_atoms: &[bool],
    fixed_bonds: &[bool],
    known_bonds: &[bool],
) -> Result<(), PotentialStereoError> {
    for (table, actual, expected) in [
        ("fixed_atoms", fixed_atoms.len(), molecule.num_atoms()),
        ("known_atoms", known_atoms.len(), molecule.num_atoms()),
        ("fixed_bonds", fixed_bonds.len(), molecule.num_bonds()),
        ("known_bonds", known_bonds.len(), molecule.num_bonds()),
    ] {
        if actual != expected {
            return Err(PotentialStereoError::InvalidStateTableLength {
                table,
                expected,
                actual,
            });
        }
    }

    // RDKit✔️✔️: void cleanMolStereo(ROMol &mol, const boost::dynamic_bitset<> &fixedAtoms,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> &knownAtoms,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> &fixedBonds,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> &knownBonds) {
    let topology = molecule.topology_block_mut();
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    for atom_index in 0..topology.atoms.len() {
        // RDKit✔️✔️:     const auto i = atom->getIdx();
        // RDKit✔️✔️:     if (!fixedAtoms[i] && knownAtoms[i]) {
        if !fixed_atoms[atom_index] && known_atoms[atom_index] {
            // RDKit✔️✔️:       switch (atom->getChiralTag()) {
            match topology.atoms[atom_index].chiral_tag() {
                // RDKit✔️✔️:         case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
                // RDKit✔️✔️:         case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
                ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw => {
                    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
                    topology.atoms[atom_index].set_chiral_tag(ChiralTag::Unspecified);
                    // RDKit✔️✔️:           for (auto nbrBond : mol.atomBonds(atom)) {
                    for neighbor in topology.adjacency.neighbors_of(atom_index) {
                        // RDKit✔️✔️:             auto bondDir = nbrBond->getBondDir();
                        let direction = topology.bonds[neighbor.bond.index()].direction();
                        // RDKit✔️✔️:             if (bondDir == Bond::BondDir::BEGINDASH ||
                        // RDKit✔️✔️:                 bondDir == Bond::BondDir::BEGINWEDGE) {
                        if matches!(
                            direction,
                            crate::BondDirection::BeginDash | crate::BondDirection::BeginWedge
                        ) {
                            // RDKit✔️✔️:               nbrBond->setBondDir(Bond::BondDir::NONE);
                            topology.bonds[neighbor.bond.index()]
                                .set_direction(crate::BondDirection::None);
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         case Atom::ChiralType::CHI_TETRAHEDRAL:
                // RDKit✔️✔️:         case Atom::ChiralType::CHI_SQUAREPLANAR:
                // RDKit✔️✔️:         case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
                // RDKit✔️✔️:         case Atom::ChiralType::CHI_OCTAHEDRAL:
                ChiralTag::Tetrahedral
                | ChiralTag::SquarePlanar
                | ChiralTag::TrigonalBipyramidal
                | ChiralTag::Octahedral => {
                    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, 0);
                    topology.atoms[atom_index].set_chiral_permutation(Some(0));
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         default:
                // RDKit✔️✔️:           break;
                _ => {} // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   bool removedStereo = false;
    let mut removed_stereo = false;
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    for bond_index in 0..topology.bonds.len() {
        // RDKit✔️✔️:     const auto i = bond->getIdx();
        // RDKit✔️✔️:     if (!fixedBonds[i] && knownBonds[i]) {
        if !fixed_bonds[bond_index] && known_bonds[bond_index] {
            // RDKit✔️✔️:       bond->setStereo(Bond::BondStereo::STEREONONE);
            // RDKit✔️✔️:       bond->setBondDir(Bond::BondDir::NONE);
            // RDKit✔️✔️:       bond->getStereoAtoms().clear();
            topology.bonds[bond_index].set_stereo(crate::BondStereo::None);
            topology.bonds[bond_index].set_direction(crate::BondDirection::None);
            topology.bonds[bond_index].set_stereo_atoms(None);
            // RDKit✔️✔️:       removedStereo = true;
            removed_stereo = true;
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   // remove any slash bond dirs that do not have a stereo neighbor
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (removedStereo) {
    if removed_stereo {
        // RDKit✔️✔️:     for (auto bond : mol.bonds()) {
        for bond_index in 0..topology.bonds.len() {
            // RDKit✔️✔️:       auto bondDir = bond->getBondDir();
            let direction = topology.bonds[bond_index].direction();
            // RDKit✔️✔️:       if (bondDir == Bond::BondDir::ENDDOWNRIGHT ||
            // RDKit✔️✔️:           bondDir == Bond::BondDir::ENDUPRIGHT) {
            if matches!(
                direction,
                crate::BondDirection::EndDownRight | crate::BondDirection::EndUpRight
            ) {
                // RDKit✔️✔️:         bool dirOk = false;
                let mut direction_ok = false;
                let endpoints = [
                    topology.bonds[bond_index].begin(),
                    topology.bonds[bond_index].end(),
                ];
                // RDKit✔️✔️:         for (auto bondEnd : {bond->getBeginAtom(), bond->getEndAtom()}) {
                for endpoint in endpoints {
                    // RDKit✔️✔️:           for (auto nbrBond : mol.atomBonds(bondEnd)) {
                    for neighbor in topology.adjacency.neighbors_of(endpoint.index()) {
                        // RDKit✔️✔️:             if (nbrBond != bond &&
                        // RDKit✔️✔️:                 nbrBond->getStereo() != Bond::BondStereo::STEREONONE) {
                        if neighbor.bond.index() != bond_index
                            && topology.bonds[neighbor.bond.index()].stereo()
                                != crate::BondStereo::None
                        {
                            // RDKit✔️✔️:               dirOk = true;
                            // RDKit✔️✔️:               break;
                            direction_ok = true;
                            break;
                            // RDKit✔️✔️:             }
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:           if (!dirOk) {
                    if !direction_ok {
                        // RDKit✔️✔️:             bond->setBondDir(Bond::BondDir::NONE);
                        topology.bonds[bond_index].set_direction(crate::BondDirection::None);
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    Ok(())
}

fn atom_total_hydrogens(molecule: &Molecule, atom: AtomId) -> Result<usize, PotentialStereoError> {
    let atom_state = molecule
        .atoms()
        .get(atom.index())
        .ok_or(PotentialStereoError::AtomIndexOutOfBounds { atom })?;
    let valence = molecule
        .derived_cache()
        .valence
        .as_ref()
        .ok_or(PotentialStereoError::MissingValenceState { atom })?;
    let implicit = *valence
        .implicit_hydrogens
        .get(atom.index())
        .ok_or(PotentialStereoError::MissingValenceState { atom })?;
    let implicit = usize::try_from(implicit).map_err(|_| {
        PotentialStereoError::InvalidImplicitHydrogenCount {
            atom,
            count: implicit,
        }
    })?;
    Ok(usize::from(atom_state.explicit_hydrogens()) + implicit)
}

pub(crate) fn get_atom_nonzero_degree(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<usize, PotentialStereoError> {
    // RDKit✔️✔️: unsigned int getAtomNonzeroDegree(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad pointer");
    // RDKit✔️✔️:   PRECONDITION(atom->hasOwningMol(), "no owning molecule");
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (auto bond : atom->getOwningMol().atomBonds(atom)) {
    // RDKit✔️✔️:     if (!bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if atom.index() >= molecule.num_atoms() {
        return Err(PotentialStereoError::AtomIndexOutOfBounds { atom });
    }
    let topology = molecule.topology_block();
    crate::chemistry::stereo::atom_nonzero_degree_from_parts(
        &topology.bonds,
        &topology.adjacency,
        atom.index(),
    )
    .map_err(PotentialStereoError::from)
}

pub(crate) fn has_protium_neighbor(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    if atom.index() >= molecule.num_atoms() {
        return Err(PotentialStereoError::AtomIndexOutOfBounds { atom });
    }
    Ok(crate::chemistry::stereo::has_protium_neighbor(
        molecule,
        atom.index(),
    ))
}

pub(crate) fn is_atom_potential_nontetrahedral_center(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // RDKit✔️✔️: bool isAtomPotentialNontetrahedralCenter(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "atom is null");
    // RDKit✔️✔️:   auto nzdegree = Chirality::detail::getAtomNonzeroDegree(atom);
    // RDKit✔️✔️:   auto impHDegree = atom->getTotalNumHs();
    // RDKit✔️✔️:   auto tnzdegree = nzdegree + impHDegree;
    // RDKit✔️✔️:   auto anum = atom->getAtomicNum();
    // RDKit✔️✔️:   if (tnzdegree > 6 || tnzdegree < 2 || (anum < 12 && anum != 4)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto chiralType = atom->getChiralTag();
    // RDKit✔️✔️:   if (chiralType >= Atom::ChiralType::CHI_SQUAREPLANAR &&
    // RDKit✔️✔️:       chiralType <= Atom::ChiralType::CHI_OCTAHEDRAL) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // with at least four neighbors but nothing specified we can start to imagine
    // RDKit✔️✔️:   // that it might be enhanced stereo
    // RDKit✔️✔️:   if (chiralType == Atom::ChiralType::CHI_UNSPECIFIED && tnzdegree >= 4) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    let atom_state = molecule
        .atoms()
        .get(atom.index())
        .ok_or(PotentialStereoError::AtomIndexOutOfBounds { atom })?;
    let total_nonzero_degree =
        get_atom_nonzero_degree(molecule, atom)? + atom_total_hydrogens(molecule, atom)?;
    let atomic_number = atom_state.atomic_number();
    if total_nonzero_degree > 6
        || total_nonzero_degree < 2
        || (atomic_number < 12 && atomic_number != 4)
    {
        return Ok(false);
    }
    if matches!(
        atom_state.chiral_tag(),
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
    ) {
        return Ok(true);
    }
    Ok(atom_state.chiral_tag() == ChiralTag::Unspecified && total_nonzero_degree >= 4)
}

pub(crate) fn is_atom_potential_tetrahedral_center(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // RDKit✔️✔️: bool isAtomPotentialTetrahedralCenter(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "atom is null");
    // RDKit✔️✔️:   auto nzDegree = getAtomNonzeroDegree(atom);
    // RDKit✔️✔️:   auto tnzDegree = nzDegree + atom->getTotalNumHs();
    // RDKit✔️✔️:   if (tnzDegree > 4) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     const auto &mol = atom->getOwningMol();
    // RDKit✔️✔️:     if (nzDegree == 4) {
    // RDKit✔️✔️:       // chirality is always possible with 4 nbrs
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     } else if (nzDegree <= 1) {
    // RDKit✔️✔️:       // chirality is never possible with 0 or 1 nbr
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     } else if (nzDegree < 3 &&
    // RDKit✔️✔️:                (atom->getAtomicNum() != 15 && atom->getAtomicNum() != 33)) {
    // RDKit✔️✔️:       // less than three neighbors is never stereogenic
    // RDKit✔️✔️:       // unless it is a phosphine/arsine with implicit H
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     } else if (atom->getAtomicNum() == 15 || atom->getAtomicNum() == 33) {
    // RDKit✔️✔️:       // from logical flow: degree is 2 or 3 (implicit H)
    // RDKit✔️✔️:       // Since InChI Software v. 1.02-standard (2009), phosphines and arsines
    // RDKit✔️✔️:       // are always treated as stereogenic even with H atom neighbors.
    // RDKit✔️✔️:       // Accept automatically.
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     } else if (nzDegree == 3) {
    // RDKit✔️✔️:       // three-coordinate with a single H we'll accept automatically:
    // RDKit✔️✔️:       if (atom->getTotalNumHs() == 1) {
    // RDKit✔️✔️:         if (detail::has_protium_neighbor(mol, atom)) {
    // RDKit✔️✔️:           // more than one H is never stereogenic
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // otherwise we default to not being a legal center
    // RDKit✔️✔️:         bool legalCenter = false;
    // RDKit✔️✔️:         // but there are a few special cases we'll accept
    // RDKit✔️✔️:         // sulfur or selenium with either a positive charge or a double
    // RDKit✔️✔️:         // bond:
    // RDKit✔️✔️:         if ((atom->getAtomicNum() == 16 || atom->getAtomicNum() == 34) &&
    // RDKit✔️✔️:             (atom->getValence(Atom::ValenceType::EXPLICIT) == 4 ||
    // RDKit✔️✔️:              (atom->getValence(Atom::ValenceType::EXPLICIT) == 3 &&
    // RDKit✔️✔️:               atom->getFormalCharge() == 1))) {
    // RDKit✔️✔️:           legalCenter = true;
    // RDKit✔️✔️:         } else if (atom->getAtomicNum() == 7) {
    // RDKit✔️✔️:           // three-coordinate N additional requirements:
    // RDKit✔️✔️:           //   in a ring of size 3  (from InChI)
    // RDKit✔️✔️:           // OR
    // RDKit✔️✔️:           /// is a bridgehead atom (RDKit extension)
    // RDKit✔️✔️:           // Also: cannot be SP2 hybridized or have a conjugated bond
    // RDKit✔️✔️:           //   (this was Github #7434)
    // RDKit✔️✔️:           if (atom->getHybridization() == Atom::HybridizationType::SP3 &&
    // RDKit✔️✔️:               !MolOps::atomHasConjugatedBond(atom) &&
    // RDKit✔️✔️:               (mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3) ||
    // RDKit✔️✔️:                queryIsAtomBridgehead(atom))) {
    // RDKit✔️✔️:             legalCenter = true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         return legalCenter;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let atom_state = molecule
        .atoms()
        .get(atom.index())
        .ok_or(PotentialStereoError::AtomIndexOutOfBounds { atom })?;
    let nonzero_degree = get_atom_nonzero_degree(molecule, atom)?;
    let total_hydrogens = atom_total_hydrogens(molecule, atom)?;
    let total_nonzero_degree = nonzero_degree + total_hydrogens;
    if total_nonzero_degree > 4 {
        return Ok(false);
    }
    if nonzero_degree == 4 {
        return Ok(true);
    }
    if nonzero_degree <= 1 {
        return Ok(false);
    }
    let atomic_number = atom_state.atomic_number();
    if nonzero_degree < 3 && atomic_number != 15 && atomic_number != 33 {
        return Ok(false);
    }
    if atomic_number == 15 || atomic_number == 33 {
        return Ok(true);
    }
    if nonzero_degree != 3 {
        return Ok(false);
    }
    if total_hydrogens == 1 {
        return Ok(!has_protium_neighbor(molecule, atom)?);
    }

    if atomic_number == 16 || atomic_number == 34 {
        let explicit_valence = molecule
            .derived_cache()
            .valence
            .as_ref()
            .and_then(|valence| valence.explicit_valence.get(atom.index()))
            .copied()
            .ok_or(PotentialStereoError::MissingValenceState { atom })?;
        return Ok(
            explicit_valence == 4 || (explicit_valence == 3 && atom_state.formal_charge() == 1)
        );
    }
    if atomic_number == 7 {
        let rings = molecule
            .derived_cache()
            .rings
            .as_ref()
            .ok_or(PotentialStereoError::MissingRingState { atom })?;
        let atom_has_conjugated_bond = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .any(|neighbor| molecule.bonds()[neighbor.bond.index()].is_conjugated());
        return Ok(atom_state.hybridization() == Hybridization::Sp3
            && !atom_has_conjugated_bond
            && (rings.is_atom_in_ring_of_size(atom, 3)
                || crate::chemistry::stereo::query_is_atom_bridgehead(
                    molecule,
                    atom.index(),
                    rings,
                ) != 0));
    }
    Ok(false)
}

pub(crate) fn is_atom_potential_stereo(
    molecule: &Molecule,
    atom: AtomId,
    allow_nontetrahedral_stereo: bool,
) -> Result<bool, PotentialStereoError> {
    // RDKit✔️✔️: bool isAtomPotentialStereoAtom(const Atom *atom,
    // RDKit✔️✔️:                                bool allowNontetrahehdralStereo) {
    // RDKit✔️✔️:   return isAtomPotentialTetrahedralCenter(atom) ||
    // RDKit✔️✔️:          (allowNontetrahehdralStereo &&
    // RDKit✔️✔️:           isAtomPotentialNontetrahedralCenter(atom));
    // RDKit✔️✔️: }
    Ok(is_atom_potential_tetrahedral_center(molecule, atom)?
        || (allow_nontetrahedral_stereo
            && is_atom_potential_nontetrahedral_center(molecule, atom)?))
}

pub(crate) fn get_atom_stereo_info(
    molecule: &Molecule,
    atom: AtomId,
    allow_nontetrahedral_stereo: bool,
) -> Result<StereoInfo, PotentialStereoError> {
    // RDKit✔️✔️: StereoInfo getStereoInfo(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "atom is null");
    // RDKit✔️✔️:   StereoInfo sinfo;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   sinfo.type = StereoType::Atom_Tetrahedral;
    // RDKit✔️✔️:   sinfo.centeredOn = atom->getIdx();
    // RDKit✔️✔️:   sinfo.controllingAtoms.reserve(atom->getDegree());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &mol = atom->getOwningMol();
    // RDKit✔️✔️:   int explicitUnknownStereo = 0;
    // RDKit✔️✔️:   for (const auto &nbri : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:     const auto &bnd = mol[nbri];
    // RDKit✔️✔️:     if (bnd->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:       explicitUnknownStereo = 1;
    // RDKit✔️✔️:     } else if (!explicitUnknownStereo) {
    // RDKit✔️✔️:       bnd->getPropIfPresent<int>(common_properties::_UnknownStereo,
    // RDKit✔️✔️:                                  explicitUnknownStereo);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     sinfo.controllingAtoms.push_back(bnd->getOtherAtomIdx(atom->getIdx()));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned> origNbrOrder = sinfo.controllingAtoms;
    // RDKit✔️✔️:   std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (explicitUnknownStereo) {
    // RDKit✔️✔️:     sinfo.specified = StereoSpecified::Unknown;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     Atom::ChiralType stereo = atom->getChiralTag();
    // RDKit✔️✔️:     if (stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW ||
    // RDKit✔️✔️:         stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:       sinfo.specified = StereoSpecified::Specified;
    // RDKit✔️✔️:       unsigned nSwaps =
    // RDKit✔️✔️:           countSwapsToInterconvert(origNbrOrder, sinfo.controllingAtoms);
    // RDKit✔️✔️:       if (nSwaps % 2) {
    // RDKit✔️✔️:         stereo = (stereo == Atom::ChiralType::CHI_TETRAHEDRAL_CCW
    // RDKit✔️✔️:                       ? Atom::ChiralType::CHI_TETRAHEDRAL_CW
    // RDKit✔️✔️:                       : Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       switch (stereo) {
    // RDKit✔️✔️:         case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
    // RDKit✔️✔️:           sinfo.descriptor = StereoDescriptor::Tet_CCW;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
    // RDKit✔️✔️:           sinfo.descriptor = StereoDescriptor::Tet_CW;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           UNDER_CONSTRUCTION("unrecognized chiral flag");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (getAllowNontetrahedralChirality() &&
    // RDKit✔️✔️:                isAtomPotentialNontetrahedralCenter(atom)) {
    // RDKit✔️✔️:       if (stereo == Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:         switch (atom->getTotalDegree()) {
    // RDKit✔️✔️:           case 4:
    // RDKit✔️✔️:             // don't assume non-tetrahedral chirality
    // RDKit✔️✔️:             stereo = Atom::ChiralType::CHI_TETRAHEDRAL;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           case 5:
    // RDKit✔️✔️:             stereo = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           case 6:
    // RDKit✔️✔️:             stereo = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           default:
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       sinfo.descriptor = StereoDescriptor::None;
    // RDKit✔️✔️:       switch (stereo) {
    // RDKit✔️✔️:         case Atom::ChiralType::CHI_TETRAHEDRAL:
    // RDKit✔️✔️:           sinfo.type = StereoType::Atom_Tetrahedral;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit✔️✔️:           sinfo.type = StereoType::Atom_SquarePlanar;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit✔️✔️:           sinfo.type = StereoType::Atom_TrigonalBipyramidal;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit✔️✔️:           sinfo.type = StereoType::Atom_Octahedral;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       unsigned int permutation;
    // RDKit✔️✔️:       if (atom->getPropIfPresent(common_properties::_chiralPermutation,
    // RDKit✔️✔️:                                  permutation)) {
    // RDKit✔️✔️:         sinfo.permutation = permutation;
    // RDKit✔️✔️:         if (!permutation) {
    // RDKit✔️✔️:           // a permutation of zero is an explicit statement that the chirality
    // RDKit✔️✔️:           // is unknown
    // RDKit✔️✔️:           sinfo.specified = Chirality::StereoSpecified::Unknown;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           sinfo.specified = Chirality::StereoSpecified::Specified;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return sinfo;
    // RDKit✔️✔️: }
    let atom_state = molecule
        .atoms()
        .get(atom.index())
        .ok_or(PotentialStereoError::AtomIndexOutOfBounds { atom })?;
    let neighbors = molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.index());
    let mut explicit_unknown_stereo = false;
    let mut original_neighbor_order = Vec::with_capacity(neighbors.len());
    for neighbor in neighbors {
        let bond = &molecule.bonds()[neighbor.bond.index()];
        if bond.direction() == crate::BondDirection::Unknown {
            explicit_unknown_stereo = true;
        } else if !explicit_unknown_stereo {
            if let Some(value) = bond.prop("_UnknownStereo") {
                explicit_unknown_stereo = value.parse::<i32>().map_err(|_| {
                    PotentialStereoError::InvalidBondIntegerProperty {
                        bond: bond.id(),
                        property: "_UnknownStereo",
                        value: value.to_string(),
                    }
                })? != 0;
            } else if bond.unknown_stereo() {
                explicit_unknown_stereo = true;
            }
        }
        original_neighbor_order.push(AtomId::new(neighbor.atom_index));
    }
    let mut sorted_neighbors = original_neighbor_order.clone();
    sorted_neighbors.sort_unstable_by_key(|neighbor| neighbor.index());

    let mut stereo_type = StereoType::AtomTetrahedral;
    let mut specified = StereoSpecified::Unspecified;
    let mut descriptor = StereoDescriptor::None;
    let mut permutation = 0;
    if explicit_unknown_stereo {
        specified = StereoSpecified::Unknown;
    } else {
        let mut stereo = atom_state.chiral_tag();
        if matches!(stereo, ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw) {
            specified = StereoSpecified::Specified;
            let swaps = crate::source_port_helpers::count_swaps_to_interconvert(
                &original_neighbor_order,
                &sorted_neighbors,
            )
            .map_err(|error| {
                PotentialStereoError::InvariantViolation(format!(
                    "atom controlling-order swap count failed: {error:?}"
                ))
            })?;
            if swaps % 2 != 0 {
                stereo = match stereo {
                    ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                    ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                    _ => unreachable!("tetrahedral tag was checked above"),
                };
            }
            descriptor = match stereo {
                ChiralTag::TetrahedralCcw => StereoDescriptor::TetrahedralCounterclockwise,
                ChiralTag::TetrahedralCw => StereoDescriptor::TetrahedralClockwise,
                _ => unreachable!("tetrahedral tag was checked above"),
            };
        } else if allow_nontetrahedral_stereo
            && is_atom_potential_nontetrahedral_center(molecule, atom)?
        {
            if stereo == ChiralTag::Unspecified {
                stereo = match neighbors.len() + atom_total_hydrogens(molecule, atom)? {
                    4 => ChiralTag::Tetrahedral,
                    5 => ChiralTag::TrigonalBipyramidal,
                    6 => ChiralTag::Octahedral,
                    _ => stereo,
                };
            }
            stereo_type = match stereo {
                ChiralTag::Tetrahedral => StereoType::AtomTetrahedral,
                ChiralTag::SquarePlanar => StereoType::AtomSquarePlanar,
                ChiralTag::TrigonalBipyramidal => StereoType::AtomTrigonalBipyramidal,
                ChiralTag::Octahedral => StereoType::AtomOctahedral,
                _ => stereo_type,
            };
            if let Some(chiral_permutation) = atom_state.chiral_permutation() {
                permutation = chiral_permutation;
                specified = if permutation == 0 {
                    StereoSpecified::Unknown
                } else {
                    StereoSpecified::Specified
                };
            }
        }
    }

    StereoInfo::new(
        stereo_type,
        specified,
        StereoCenter::Atom(atom),
        descriptor,
        permutation,
        sorted_neighbors
            .into_iter()
            .map(ControllingAtom::Atom)
            .collect(),
    )
}

fn atom_total_hydrogens_including_neighbors(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<usize, PotentialStereoError> {
    let attached = molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| molecule.atoms()[neighbor.atom_index].atomic_number() == 1)
        .count();
    Ok(atom_total_hydrogens(molecule, atom)? + attached)
}

pub(crate) fn is_bond_potential_stereo(
    molecule: &Molecule,
    bond: BondId,
) -> Result<bool, PotentialStereoError> {
    // RDKit✔️✔️: bool isBondPotentialStereoBond(const Bond *bond) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bond is null");
    // RDKit✔️✔️:   if (bond->getBondType() != Bond::BondType::DOUBLE) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // at the moment the condition for being a potential stereo bond is that
    // RDKit✔️✔️:   // each of the beginning and end neighbors must have at least 2 explicit
    // RDKit✔️✔️:   // neighbors but no more than 3 total neighbors.
    // RDKit✔️✔️:   // if it's a ring bond, the smallest ring it's in must have at least 8
    // RDKit✔️✔️:   // members
    // RDKit✔️✔️:   //  (this is common with InChI)
    // RDKit✔️✔️:   const auto beginAtom = bond->getBeginAtom();
    // RDKit✔️✔️:   auto begDegree = beginAtom->getTotalDegree();
    // RDKit✔️✔️:   const auto endAtom = bond->getEndAtom();
    // RDKit✔️✔️:   auto endDegree = endAtom->getTotalDegree();
    // RDKit✔️✔️:   if (begDegree > 1 && begDegree < 4 && endDegree > 1 && endDegree < 4 &&
    // RDKit✔️✔️:       beginAtom->getTotalNumHs(true) < 2 && endAtom->getTotalNumHs(true) < 2) {
    // RDKit✔️✔️:     // check rings
    // RDKit✔️✔️:     const auto ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️✔️:     for (const auto &bring : ri->bondRings()) {
    // RDKit✔️✔️:       if (bring.size() < minRingSizeForDoubleBondStereo &&
    // RDKit✔️✔️:           std::find(bring.begin(), bring.end(), bond->getIdx()) !=
    // RDKit✔️✔️:               bring.end()) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let bond_state = molecule
        .bonds()
        .get(bond.index())
        .ok_or(PotentialStereoError::BondIndexOutOfBounds { bond })?;
    if bond_state.order() != crate::BondOrder::Double {
        return Ok(false);
    }
    let begin = bond_state.begin();
    let end = bond_state.end();
    let begin_degree = molecule
        .topology_block()
        .adjacency
        .neighbors_of(begin.index())
        .len()
        + atom_total_hydrogens(molecule, begin)?;
    let end_degree = molecule
        .topology_block()
        .adjacency
        .neighbors_of(end.index())
        .len()
        + atom_total_hydrogens(molecule, end)?;
    if !(begin_degree > 1
        && begin_degree < 4
        && end_degree > 1
        && end_degree < 4
        && atom_total_hydrogens_including_neighbors(molecule, begin)? < 2
        && atom_total_hydrogens_including_neighbors(molecule, end)? < 2)
    {
        return Ok(false);
    }
    let rings = molecule
        .derived_cache()
        .rings
        .as_ref()
        .ok_or(PotentialStereoError::MissingBondRingState { bond })?;
    Ok(!rings
        .bond_rings()
        .iter()
        .any(|ring| ring.len() < 8 && ring.contains(&bond)))
}

fn atom_unknown_stereo_property(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    let atom_state = &molecule.atoms()[atom.index()];
    if let Some(value) = atom_state.prop("_UnknownStereo") {
        return value.parse::<i32>().map(|value| value != 0).map_err(|_| {
            PotentialStereoError::InvalidAtomIntegerProperty {
                atom,
                property: "_UnknownStereo",
                value: value.to_string(),
            }
        });
    }
    Ok(atom_state.unknown_stereo())
}

fn explore_bond_end(
    molecule: &Molecule,
    bond: BondId,
    atom: AtomId,
    controlling_atoms: &mut Vec<ControllingAtom>,
    seen_squiggle_bond: &mut bool,
) {
    let neighbors = molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.index());
    for neighbor in neighbors {
        if neighbor.bond == bond {
            continue;
        }
        if molecule.bonds()[neighbor.bond.index()].direction() == crate::BondDirection::Unknown {
            *seen_squiggle_bond = true;
        }
        controlling_atoms.push(ControllingAtom::Atom(AtomId::new(neighbor.atom_index)));
    }
    for _ in neighbors.len()..3 {
        controlling_atoms.push(ControllingAtom::Missing);
    }
}

pub(crate) fn get_bond_stereo_info(
    molecule: &Molecule,
    bond: BondId,
) -> Result<StereoInfo, PotentialStereoError> {
    let bond_state = molecule
        .bonds()
        .get(bond.index())
        .ok_or(PotentialStereoError::BondIndexOutOfBounds { bond })?;
    let begin = bond_state.begin();
    let end = bond_state.end();
    let begin_degree = molecule
        .topology_block()
        .adjacency
        .neighbors_of(begin.index())
        .len();
    let end_degree = molecule
        .topology_block()
        .adjacency
        .neighbors_of(end.index())
        .len();

    if bond_state.order() == crate::BondOrder::Double {
        // RDKit✔️✔️: if (bond->getBondType() == Bond::BondType::DOUBLE) {
        // RDKit✔️✔️:   if (beginAtom->getDegree() < 1 || endAtom->getDegree() < 1 ||
        // RDKit✔️✔️:       beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
        // RDKit✔️✔️:     throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   sinfo.type = StereoType::Bond_Double;
        // RDKit✔️✔️:   sinfo.centeredOn = bond->getIdx();
        // RDKit✔️✔️:   sinfo.controllingAtoms.reserve(4);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   bool seenSquiggleBond = false;
        // RDKit✔️✔️:   const auto &mol = bond->getOwningMol();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   auto explore_bond_end = [&mol, &bond, &sinfo,
        // RDKit✔️✔️:                            &seenSquiggleBond](const Atom *atom) {
        // RDKit✔️✔️:     for (const auto nbr : mol.atomBonds(atom)) {
        // RDKit✔️✔️:       if (nbr->getIdx() != bond->getIdx()) {
        // RDKit✔️✔️:         if (nbr->getBondDir() == Bond::BondDir::UNKNOWN) {
        // RDKit✔️✔️:           seenSquiggleBond = true;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:         sinfo.controllingAtoms.push_back(
        // RDKit✔️✔️:             nbr->getOtherAtomIdx(atom->getIdx()));
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     for (unsigned i = atom->getDegree(); i < 3; ++i) {
        // RDKit✔️✔️:       sinfo.controllingAtoms.push_back(Atom::NOATOM);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   };
        // RDKit✔️✔️:
        // RDKit✔️✔️:   explore_bond_end(beginAtom);
        // RDKit✔️✔️:   explore_bond_end(endAtom);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (!seenSquiggleBond) {
        // RDKit✔️✔️:     // check to see if either the begin or end atoms has the _UnknownStereo
        // RDKit✔️✔️:     // property set. This happens if there was a squiggle bond to an H
        // RDKit✔️✔️:     int explicitUnknownStereo = 0;
        // RDKit✔️✔️:     if ((bond->getBeginAtom()->getPropIfPresent<int>(
        // RDKit✔️✔️:              common_properties::_UnknownStereo, explicitUnknownStereo) &&
        // RDKit✔️✔️:          explicitUnknownStereo) ||
        // RDKit✔️✔️:         (bond->getEndAtom()->getPropIfPresent<int>(
        // RDKit✔️✔️:              common_properties::_UnknownStereo, explicitUnknownStereo) &&
        // RDKit✔️✔️:          explicitUnknownStereo)) {
        // RDKit✔️✔️:       seenSquiggleBond = true;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   Bond::BondStereo stereo = bond->getStereo();
        // RDKit✔️✔️:   if (stereo == Bond::BondStereo::STEREOANY ||
        // RDKit✔️✔️:       bond->getBondDir() == Bond::BondDir::EITHERDOUBLE || seenSquiggleBond) {
        // RDKit✔️✔️:     sinfo.specified = Chirality::StereoSpecified::Unknown;
        // RDKit✔️✔️:   } else if (stereo != Bond::BondStereo::STEREONONE) {
        // RDKit✔️✔️:     if (stereo == Bond::BondStereo::STEREOE ||
        // RDKit✔️✔️:         stereo == Bond::BondStereo::STEREOZ) {
        // RDKit✔️✔️:       stereo = Chirality::translateEZLabelToCisTrans(stereo);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     sinfo.specified = Chirality::StereoSpecified::Specified;
        // RDKit✔️✔️:     const auto satoms = bond->getStereoAtoms();
        // RDKit✔️✔️:     if (satoms.size() != 2) {
        // RDKit✔️✔️:       throw ValueErrorException("only can support 2 stereo neighbors");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     bool firstAtBegin;
        // RDKit✔️✔️:     if (satoms[0] == static_cast<int>(sinfo.controllingAtoms[0])) {
        // RDKit✔️✔️:       firstAtBegin = true;
        // RDKit✔️✔️:     } else if (satoms[0] == static_cast<int>(sinfo.controllingAtoms[1])) {
        // RDKit✔️✔️:       firstAtBegin = false;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       throw ValueErrorException("controlling atom mismatch at begin");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     bool firstAtEnd;
        // RDKit✔️✔️:     if (satoms[1] == static_cast<int>(sinfo.controllingAtoms[2])) {
        // RDKit✔️✔️:       firstAtEnd = true;
        // RDKit✔️✔️:     } else if (satoms[1] == static_cast<int>(sinfo.controllingAtoms[3])) {
        // RDKit✔️✔️:       firstAtEnd = false;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       throw ValueErrorException("controlling atom mismatch at end");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     auto mismatch = firstAtBegin ^ firstAtEnd;
        // RDKit✔️✔️:     if (mismatch) {
        // RDKit✔️✔️:       stereo = (stereo == Bond::BondStereo::STEREOCIS
        // RDKit✔️✔️:                     ? Bond::BondStereo::STEREOTRANS
        // RDKit✔️✔️:                     : Bond::BondStereo::STEREOCIS);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     switch (stereo) {
        // RDKit✔️✔️:       case Bond::BondStereo::STEREOCIS:
        // RDKit✔️✔️:         sinfo.descriptor = Chirality::StereoDescriptor::Bond_Cis;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       case Bond::BondStereo::STEREOTRANS:
        // RDKit✔️✔️:         sinfo.descriptor = Chirality::StereoDescriptor::Bond_Trans;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       default:
        // RDKit✔️✔️:         UNDER_CONSTRUCTION("unrecognized bond stereo type");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     sinfo.specified = Chirality::StereoSpecified::Unspecified;
        // RDKit✔️✔️:   }
        if !(1..=3).contains(&begin_degree) || !(1..=3).contains(&end_degree) {
            return Err(PotentialStereoError::InvalidBondDegree { bond });
        }
        let mut controlling_atoms = Vec::with_capacity(4);
        let mut seen_squiggle_bond = false;
        explore_bond_end(
            molecule,
            bond,
            begin,
            &mut controlling_atoms,
            &mut seen_squiggle_bond,
        );
        explore_bond_end(
            molecule,
            bond,
            end,
            &mut controlling_atoms,
            &mut seen_squiggle_bond,
        );
        if !seen_squiggle_bond {
            seen_squiggle_bond = atom_unknown_stereo_property(molecule, begin)?
                || atom_unknown_stereo_property(molecule, end)?;
        }
        let mut specified = StereoSpecified::Unspecified;
        let mut descriptor = StereoDescriptor::None;
        let mut stereo = bond_state.stereo();
        if stereo == crate::BondStereo::Any
            || bond_state.direction() == crate::BondDirection::EitherDouble
            || seen_squiggle_bond
        {
            specified = StereoSpecified::Unknown;
        } else if stereo != crate::BondStereo::None {
            stereo = match stereo {
                crate::BondStereo::E => crate::BondStereo::Trans,
                crate::BondStereo::Z => crate::BondStereo::Cis,
                other => other,
            };
            specified = StereoSpecified::Specified;
            let stereo_atoms = bond_state
                .stereo_atoms()
                .ok_or(PotentialStereoError::InvalidStereoAtomCount { bond })?;
            let first_at_begin = if controlling_atoms[0] == ControllingAtom::Atom(stereo_atoms[0]) {
                true
            } else if controlling_atoms[1] == ControllingAtom::Atom(stereo_atoms[0]) {
                false
            } else {
                return Err(PotentialStereoError::StereoAtomMismatch {
                    bond,
                    stereo_atom: stereo_atoms[0],
                    end: "begin",
                });
            };
            let first_at_end = if controlling_atoms[2] == ControllingAtom::Atom(stereo_atoms[1]) {
                true
            } else if controlling_atoms[3] == ControllingAtom::Atom(stereo_atoms[1]) {
                false
            } else {
                return Err(PotentialStereoError::StereoAtomMismatch {
                    bond,
                    stereo_atom: stereo_atoms[1],
                    end: "end",
                });
            };
            if first_at_begin ^ first_at_end {
                stereo = if stereo == crate::BondStereo::Cis {
                    crate::BondStereo::Trans
                } else {
                    crate::BondStereo::Cis
                };
            }
            descriptor = match stereo {
                crate::BondStereo::Cis => StereoDescriptor::BondCis,
                crate::BondStereo::Trans => StereoDescriptor::BondTrans,
                _ => {
                    return Err(PotentialStereoError::InvariantViolation(format!(
                        "bond {bond} has unsupported specified stereo {stereo:?}"
                    )));
                }
            };
        }
        return StereoInfo::new(
            StereoType::BondDouble,
            specified,
            StereoCenter::Bond(bond),
            descriptor,
            0,
            controlling_atoms,
        );
    }

    if bond_state.order() == crate::BondOrder::Single
        && matches!(
            bond_state.stereo(),
            crate::BondStereo::AtropCcw | crate::BondStereo::AtropCw
        )
    {
        // RDKit✔️✔️: } else if (bond->getBondType() == Bond::BondType::SINGLE &&
        // RDKit✔️✔️:            (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
        // RDKit✔️✔️:             bond->getStereo() == Bond::BondStereo::STEREOATROPCW)) {
        // RDKit✔️✔️:   if (beginAtom->getDegree() < 2 || endAtom->getDegree() < 2 ||
        // RDKit✔️✔️:       beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
        // RDKit✔️✔️:     throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   sinfo.type = StereoType::Bond_Atropisomer;
        // RDKit✔️✔️:   sinfo.centeredOn = bond->getIdx();
        // RDKit✔️✔️:   sinfo.controllingAtoms.reserve(4);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const auto &mol = bond->getOwningMol();
        // RDKit✔️✔️:   for (const auto nbr : mol.atomBonds(beginAtom)) {
        // RDKit✔️✔️:     if (nbr->getIdx() != bond->getIdx()) {
        // RDKit✔️✔️:       sinfo.controllingAtoms.push_back(
        // RDKit✔️✔️:           nbr->getOtherAtomIdx(beginAtom->getIdx()));
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (beginAtom->getDegree() == 2) {
        // RDKit✔️✔️:     sinfo.controllingAtoms.push_back(Atom::NOATOM);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (const auto nbr : mol.atomBonds(endAtom)) {
        // RDKit✔️✔️:     if (nbr->getIdx() != bond->getIdx()) {
        // RDKit✔️✔️:       sinfo.controllingAtoms.push_back(
        // RDKit✔️✔️:           nbr->getOtherAtomIdx(endAtom->getIdx()));
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (endAtom->getDegree() == 2) {
        // RDKit✔️✔️:     sinfo.controllingAtoms.push_back(Atom::NOATOM);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   Bond::BondStereo stereo = bond->getStereo();
        // RDKit✔️✔️:   sinfo.specified = Chirality::StereoSpecified::Specified;
        // RDKit✔️✔️:   switch (stereo) {
        // RDKit✔️✔️:     case Bond::BondStereo::STEREOATROPCW:
        // RDKit✔️✔️:       sinfo.descriptor = Chirality::StereoDescriptor::Bond_AtropCW;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case Bond::BondStereo::STEREOATROPCCW:
        // RDKit✔️✔️:       sinfo.descriptor = Chirality::StereoDescriptor::Bond_AtropCCW;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       UNDER_CONSTRUCTION("unrecognized bond stereo type");
        // RDKit✔️✔️:   }
        if !(2..=3).contains(&begin_degree) || !(2..=3).contains(&end_degree) {
            return Err(PotentialStereoError::InvalidBondDegree { bond });
        }
        let mut controlling_atoms = Vec::with_capacity(4);
        for endpoint in [begin, end] {
            for neighbor in molecule
                .topology_block()
                .adjacency
                .neighbors_of(endpoint.index())
            {
                if neighbor.bond != bond {
                    controlling_atoms.push(ControllingAtom::Atom(AtomId::new(neighbor.atom_index)));
                }
            }
            if molecule
                .topology_block()
                .adjacency
                .neighbors_of(endpoint.index())
                .len()
                == 2
            {
                controlling_atoms.push(ControllingAtom::Missing);
            }
        }
        let descriptor = match bond_state.stereo() {
            crate::BondStereo::AtropCw => StereoDescriptor::BondAtropClockwise,
            crate::BondStereo::AtropCcw => StereoDescriptor::BondAtropCounterclockwise,
            _ => unreachable!("atrop stereo checked above"),
        };
        return StereoInfo::new(
            StereoType::BondAtropisomer,
            StereoSpecified::Specified,
            StereoCenter::Bond(bond),
            descriptor,
            0,
            controlling_atoms,
        );
    }

    StereoInfo::new(
        StereoType::Unspecified,
        StereoSpecified::Unspecified,
        StereoCenter::Missing,
        StereoDescriptor::None,
        0,
        Vec::new(),
    )
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    fn ring_test_molecule(size: usize) -> Molecule {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let atoms = (0..size)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for index in 0..size {
            builder
                .add_bond(BondSpec::new(
                    atoms[index],
                    atoms[(index + 1) % size],
                    BondOrder::Single,
                ))
                .unwrap();
        }
        let mut molecule = builder.build().unwrap();
        molecule.derived_cache_mut().rings = Some(crate::symmetrize_sssr(&molecule).unwrap());
        molecule
    }

    fn ring_stereo_counts(
        molecule: &mut Molecule,
        possible_atom_indices: &[usize],
    ) -> (Vec<u32>, Vec<u32>) {
        let known_atoms = vec![false; molecule.num_atoms()];
        let mut possible_atoms = vec![false; molecule.num_atoms()];
        for atom in possible_atom_indices {
            possible_atoms[*atom] = true;
        }
        let known_bonds = vec![false; molecule.num_bonds()];
        let possible_bonds = vec![false; molecule.num_bonds()];
        let mut atom_counts = vec![0; molecule.num_atoms()];
        let mut bond_counts = vec![0; molecule.num_bonds()];
        flag_ring_stereo(
            molecule,
            &mut atom_counts,
            &mut bond_counts,
            &known_atoms,
            Some(&possible_atoms),
            &known_bonds,
            Some(&possible_bonds),
        )
        .unwrap();
        (atom_counts, bond_counts)
    }

    struct IterativeUpdateState {
        initialization: PotentialStereoInitialization,
        fixed_atoms: Vec<bool>,
        fixed_bonds: Vec<bool>,
        possible_ring_stereo_atoms: Vec<u32>,
        possible_ring_stereo_bonds: Vec<u32>,
    }

    fn iterative_update_state(molecule: &mut Molecule) -> IterativeUpdateState {
        let initialization = initialize_potential_stereo(molecule, true, false, true).unwrap();
        let mut possible_ring_stereo_atoms = vec![0; molecule.num_atoms()];
        let mut possible_ring_stereo_bonds = vec![0; molecule.num_bonds()];
        flag_ring_stereo(
            molecule,
            &mut possible_ring_stereo_atoms,
            &mut possible_ring_stereo_bonds,
            &initialization.known_atoms,
            Some(&initialization.possible_atoms),
            &initialization.known_bonds,
            Some(&initialization.possible_bonds),
        )
        .unwrap();
        IterativeUpdateState {
            fixed_atoms: vec![false; molecule.num_atoms()],
            fixed_bonds: vec![false; molecule.num_bonds()],
            initialization,
            possible_ring_stereo_atoms,
            possible_ring_stereo_bonds,
        }
    }

    fn run_iterative_update_round(
        molecule: &Molecule,
        state: &mut IterativeUpdateState,
    ) -> (bool, Vec<StereoInfo>, Vec<u32>) {
        use crate::canon_rank::{CanonicalRankOptions, FragmentRankScope, rank_fragment_atoms};

        let atoms_in_play = vec![true; molecule.num_atoms()];
        let bonds_in_play = vec![true; molecule.num_bonds()];
        let scope = FragmentRankScope::new(&atoms_in_play, &bonds_in_play)
            .with_atom_symbols(&state.initialization.atom_symbols)
            .with_bond_symbols(&state.initialization.bond_symbols);
        let ranks = rank_fragment_atoms(
            molecule,
            scope,
            CanonicalRankOptions {
                break_ties: false,
                include_chirality: false,
                include_isotopes: false,
                include_atom_maps: false,
                include_chiral_presence: false,
                include_stereo_groups: false,
                use_non_stereo_ranks: false,
                include_ring_stereo: false,
                chirality_rings_use_ring_stereo: false,
            },
        )
        .unwrap()
        .into_iter()
        .map(|rank| u32::try_from(rank).unwrap())
        .collect::<Vec<_>>();
        let mut records = Vec::new();
        let mut another_round = update_atoms(
            molecule,
            &ranks,
            &mut state.initialization.atom_symbols,
            &mut state.initialization.possible_atoms,
            &state.initialization.known_atoms,
            &mut state.fixed_atoms,
            &mut state.possible_ring_stereo_atoms,
            &mut state.possible_ring_stereo_bonds,
            &mut records,
            true,
        )
        .unwrap();
        another_round |= update_bonds(
            molecule,
            &ranks,
            &mut state.initialization.bond_symbols,
            &state.initialization.possible_atoms,
            &mut state.initialization.possible_bonds,
            &state.initialization.known_atoms,
            &state.initialization.known_bonds,
            &mut state.fixed_bonds,
            &mut records,
        )
        .unwrap();
        (another_round, records, ranks)
    }

    #[test]
    fn stereo_type_model_covers_every_source_variant_and_center_kind() {
        let atom = AtomId::new(0);
        let bond = BondId::new(0);
        let cases = [
            (StereoType::Unspecified, StereoCenter::Missing),
            (StereoType::AtomTetrahedral, StereoCenter::Atom(atom)),
            (StereoType::AtomSquarePlanar, StereoCenter::Atom(atom)),
            (
                StereoType::AtomTrigonalBipyramidal,
                StereoCenter::Atom(atom),
            ),
            (StereoType::AtomOctahedral, StereoCenter::Atom(atom)),
            (StereoType::BondDouble, StereoCenter::Bond(bond)),
            (StereoType::BondEvenCumulene, StereoCenter::Bond(bond)),
            (StereoType::BondAtropisomer, StereoCenter::Bond(bond)),
        ];

        for (stereo_type, center) in cases {
            let info = StereoInfo::new(
                stereo_type,
                StereoSpecified::Unspecified,
                center,
                StereoDescriptor::None,
                0,
                Vec::new(),
            )
            .expect("source stereo type must accept its typed center");
            assert_eq!(info.stereo_type(), stereo_type);
            assert_eq!(info.center(), center);
        }
    }

    #[test]
    fn descriptor_and_specified_variants_have_stable_equality_and_ordering() {
        let descriptors = [
            StereoDescriptor::None,
            StereoDescriptor::TetrahedralClockwise,
            StereoDescriptor::TetrahedralCounterclockwise,
            StereoDescriptor::BondCis,
            StereoDescriptor::BondTrans,
            StereoDescriptor::BondAtropClockwise,
            StereoDescriptor::BondAtropCounterclockwise,
        ];
        let specifications = [
            StereoSpecified::Unspecified,
            StereoSpecified::Specified,
            StereoSpecified::Unknown,
        ];

        assert_eq!(descriptors.into_iter().collect::<BTreeSet<_>>().len(), 7);
        assert_eq!(specifications.into_iter().collect::<BTreeSet<_>>().len(), 3);

        let first = StereoInfo::new(
            StereoType::AtomTetrahedral,
            StereoSpecified::Unspecified,
            StereoCenter::Atom(AtomId::new(0)),
            StereoDescriptor::None,
            0,
            vec![ControllingAtom::Atom(AtomId::new(1))],
        )
        .unwrap();
        let mut second = first.clone();
        assert_eq!(first, second);
        second.permutation = 1;
        assert!(first < second);
    }

    #[test]
    fn controlling_atom_model_preserves_missing_slots_and_order() {
        let info = StereoInfo::new(
            StereoType::BondDouble,
            StereoSpecified::Unknown,
            StereoCenter::Bond(BondId::new(2)),
            StereoDescriptor::BondTrans,
            0,
            vec![
                ControllingAtom::Atom(AtomId::new(4)),
                ControllingAtom::Missing,
                ControllingAtom::Atom(AtomId::new(1)),
                ControllingAtom::Missing,
            ],
        )
        .unwrap();

        assert_eq!(
            info.controlling_atoms(),
            [
                ControllingAtom::Atom(AtomId::new(4)),
                ControllingAtom::Missing,
                ControllingAtom::Atom(AtomId::new(1)),
                ControllingAtom::Missing,
            ]
        );
        assert_eq!(info.specified(), StereoSpecified::Unknown);
        assert_eq!(info.descriptor(), StereoDescriptor::BondTrans);
        assert_eq!(info.permutation(), 0);
    }

    #[test]
    fn constructor_rejects_center_type_mismatches() {
        let error = StereoInfo::new(
            StereoType::AtomTetrahedral,
            StereoSpecified::Specified,
            StereoCenter::Bond(BondId::new(0)),
            StereoDescriptor::TetrahedralClockwise,
            0,
            Vec::new(),
        )
        .unwrap_err();
        assert_eq!(
            error,
            PotentialStereoError::CenterTypeMismatch {
                stereo_type: StereoType::AtomTetrahedral,
                center: StereoCenter::Bond(BondId::new(0)),
            }
        );
    }

    #[test]
    fn index_validation_checks_centers_and_controlling_atoms() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let valid = StereoInfo::new(
            StereoType::BondDouble,
            StereoSpecified::Unspecified,
            StereoCenter::Bond(BondId::new(0)),
            StereoDescriptor::None,
            0,
            vec![
                ControllingAtom::Atom(AtomId::new(0)),
                ControllingAtom::Missing,
            ],
        )
        .unwrap();
        assert_eq!(valid.validate_indices(&molecule), Ok(()));

        let invalid_center = StereoInfo::new(
            StereoType::AtomTetrahedral,
            StereoSpecified::Unspecified,
            StereoCenter::Atom(AtomId::new(2)),
            StereoDescriptor::None,
            0,
            Vec::new(),
        )
        .unwrap();
        assert_eq!(
            invalid_center.validate_indices(&molecule),
            Err(PotentialStereoError::AtomIndexOutOfBounds {
                atom: AtomId::new(2)
            })
        );

        let invalid_controller = StereoInfo::new(
            StereoType::BondDouble,
            StereoSpecified::Unspecified,
            StereoCenter::Bond(BondId::new(0)),
            StereoDescriptor::None,
            0,
            vec![ControllingAtom::Atom(AtomId::new(2))],
        )
        .unwrap();
        assert_eq!(
            invalid_controller.validate_indices(&molecule),
            Err(PotentialStereoError::AtomIndexOutOfBounds {
                atom: AtomId::new(2)
            })
        );

        let invalid_bond = StereoInfo::new(
            StereoType::BondDouble,
            StereoSpecified::Unspecified,
            StereoCenter::Bond(BondId::new(1)),
            StereoDescriptor::None,
            0,
            Vec::new(),
        )
        .unwrap();
        assert_eq!(
            invalid_bond.validate_indices(&molecule),
            Err(PotentialStereoError::BondIndexOutOfBounds {
                bond: BondId::new(1)
            })
        );
    }

    #[test]
    fn analysis_options_match_the_python_find_potential_stereo_defaults() {
        assert_eq!(
            PotentialStereoOptions::default(),
            PotentialStereoOptions {
                clean: false,
                flag_possible: true,
            }
        );
    }

    #[test]
    fn analysis_result_owns_the_analyzed_value_and_ordered_records() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let info = StereoInfo::new(
            StereoType::Unspecified,
            StereoSpecified::Unspecified,
            StereoCenter::Missing,
            StereoDescriptor::None,
            0,
            Vec::new(),
        )
        .unwrap();
        let analysis = PotentialStereoAnalysis::new(molecule.clone(), vec![info.clone()]);

        assert_eq!(analysis.molecule, molecule);
        assert_eq!(analysis.stereo_info, vec![info]);
    }

    #[test]
    fn public_analysis_method_and_function_share_one_typed_core() {
        let source = Molecule::from_smiles("CC(F)=CC(Cl)C").unwrap();
        let options = PotentialStereoOptions::default();

        let via_method = source.analyze_potential_stereo(options).unwrap();
        let via_function = analyze_potential_stereo(&source, options).unwrap();

        assert_eq!(via_method.molecule, via_function.molecule);
        assert_eq!(via_method.stereo_info, via_function.stereo_info);
        assert_eq!(via_method.stereo_info.len(), 2);
        assert!(via_method.stereo_info.iter().any(|info| {
            info.stereo_type() == StereoType::AtomTetrahedral
                && matches!(info.center(), StereoCenter::Atom(_))
        }));
        assert!(via_method.stereo_info.iter().any(|info| {
            info.stereo_type() == StereoType::BondDouble
                && matches!(info.center(), StereoCenter::Bond(_))
        }));
        for info in &via_method.stereo_info {
            info.validate_indices(&via_method.molecule).unwrap();
        }
    }

    #[test]
    fn public_analysis_returns_cleanup_on_owned_value_and_preserves_source() {
        let mut source = Molecule::from_smiles("CCF").unwrap();
        source.topology_block_mut().atoms[1].set_chiral_tag(ChiralTag::TetrahedralCcw);
        let before = source.clone();
        assert!(source.potential_stereo_cache().is_none());

        let analysis = source
            .analyze_potential_stereo(PotentialStereoOptions {
                clean: true,
                flag_possible: true,
            })
            .unwrap();

        assert_eq!(source, before);
        assert_eq!(source.atoms()[1].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert!(source.potential_stereo_cache().is_none());
        assert_eq!(
            analysis.molecule.atoms()[1].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(
            analysis.molecule.potential_stereo_cache().as_deref(),
            Some(analysis.stereo_info.as_slice())
        );
    }

    #[test]
    fn public_analysis_cleaning_removes_source_redundant_direction_specs() {
        let source = Molecule::from_smiles("Br.Br/C(=N\\N=c1/nn[nH][nH]1)c1ccncc1").unwrap();
        let source_before = source.clone();

        for flag_possible in [false, true] {
            let analysis = source
                .analyze_potential_stereo(PotentialStereoOptions {
                    clean: true,
                    flag_possible,
                })
                .unwrap();
            assert_eq!(
                analysis.molecule.to_smiles(true).unwrap(),
                "Br.BrC(=NN=c1nn[nH][nH]1)c1ccncc1"
            );
            assert_eq!(
                analysis.molecule.bonds()[0].direction(),
                crate::BondDirection::None
            );
            assert_eq!(
                analysis.molecule.bonds()[2].direction(),
                crate::BondDirection::EndDownRight
            );
            assert_eq!(
                analysis.molecule.bonds()[4].direction(),
                crate::BondDirection::None
            );
        }

        assert_eq!(source, source_before);
    }

    #[test]
    fn public_analysis_is_equivalent_from_cold_and_warm_cache_state() {
        let source = Molecule::from_smiles("C[C@H](F)C(C)[C@H](F)C").unwrap();
        let cold = source
            .analyze_potential_stereo(PotentialStereoOptions::default())
            .unwrap();
        let warm = cold
            .molecule
            .analyze_potential_stereo(PotentialStereoOptions::default())
            .unwrap();

        assert_eq!(cold.stereo_info, warm.stereo_info);
        assert_eq!(cold.molecule, warm.molecule);
        assert_eq!(source.to_smiles(true).unwrap(), "CC([C@H](C)F)[C@@H](C)F");
    }

    #[test]
    fn public_analysis_returns_structured_errors_without_mutating_source() {
        let mut source = Molecule::from_smiles("F[C@](Cl)(Br)I").unwrap();
        source.topology_block_mut().bonds[0].set_prop("_UnknownStereo", "not-an-integer");
        let before = source.clone();

        let error = source
            .analyze_potential_stereo(PotentialStereoOptions::default())
            .unwrap_err();

        assert_eq!(source, before);
        assert_eq!(
            error,
            PotentialStereoError::InvalidBondIntegerProperty {
                bond: BondId::new(0),
                property: "_UnknownStereo",
                value: "not-an-integer".to_string(),
            }
        );
    }

    #[test]
    fn atom_update_fixes_unique_center_inverts_descriptor_and_converges() {
        let molecule = Molecule::from_smiles("F[C@](Cl)(Br)I").unwrap();
        let center = AtomId::new(1);
        let source_info = get_atom_stereo_info(&molecule, center, true).unwrap();
        let controllers = source_info
            .controlling_atoms()
            .iter()
            .map(|controller| match controller {
                ControllingAtom::Atom(atom) => *atom,
                ControllingAtom::Missing => panic!("tetrahedral fixture has four controllers"),
            })
            .collect::<Vec<_>>();
        assert_eq!(controllers.len(), 4);

        let mut ranks = vec![0; molecule.num_atoms()];
        for (controller, rank) in controllers.iter().zip([2, 1, 3, 4]) {
            ranks[controller.index()] = rank;
        }
        let mut atom_symbols = molecule
            .atoms()
            .iter()
            .map(atom_compare_symbol)
            .collect::<Vec<_>>();
        let mut possible_atoms = vec![false; molecule.num_atoms()];
        let mut known_atoms = vec![false; molecule.num_atoms()];
        known_atoms[center.index()] = true;
        let mut fixed_atoms = vec![false; molecule.num_atoms()];
        let mut possible_ring_atoms = vec![0; molecule.num_atoms()];
        let mut possible_ring_bonds = vec![0; molecule.num_bonds()];
        let mut records = Vec::new();

        assert!(
            update_atoms(
                &molecule,
                &ranks,
                &mut atom_symbols,
                &mut possible_atoms,
                &known_atoms,
                &mut fixed_atoms,
                &mut possible_ring_atoms,
                &mut possible_ring_bonds,
                &mut records,
                true,
            )
            .unwrap()
        );
        assert!(fixed_atoms[center.index()]);
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].center(), StereoCenter::Atom(center));
        assert_eq!(
            records[0].descriptor(),
            match source_info.descriptor() {
                StereoDescriptor::TetrahedralClockwise => {
                    StereoDescriptor::TetrahedralCounterclockwise
                }
                StereoDescriptor::TetrahedralCounterclockwise => {
                    StereoDescriptor::TetrahedralClockwise
                }
                descriptor => panic!("unexpected tetrahedral descriptor {descriptor:?}"),
            }
        );

        records.clear();
        assert!(
            !update_atoms(
                &molecule,
                &ranks,
                &mut atom_symbols,
                &mut possible_atoms,
                &known_atoms,
                &mut fixed_atoms,
                &mut possible_ring_atoms,
                &mut possible_ring_bonds,
                &mut records,
                true,
            )
            .unwrap()
        );
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn bond_update_orders_both_ends_inverts_descriptor_and_converges() {
        let molecule = Molecule::from_smiles("F/C(Cl)=C(Br)/I").unwrap();
        let bond = molecule
            .bonds()
            .iter()
            .find(|bond| bond.order() == crate::BondOrder::Double)
            .unwrap()
            .id();
        let source_info = get_bond_stereo_info(&molecule, bond).unwrap();
        assert_eq!(source_info.specified(), StereoSpecified::Specified);
        let controllers = source_info
            .controlling_atoms()
            .iter()
            .map(|controller| match controller {
                ControllingAtom::Atom(atom) => *atom,
                ControllingAtom::Missing => panic!("fully substituted alkene has four controllers"),
            })
            .collect::<Vec<_>>();

        let mut ranks = vec![0; molecule.num_atoms()];
        for (controller, rank) in controllers.iter().zip([1, 2, 4, 3]) {
            ranks[controller.index()] = rank;
        }
        let mut bond_symbols = molecule
            .bonds()
            .iter()
            .map(bond_compare_symbol)
            .collect::<Vec<_>>();
        let possible_atoms = vec![false; molecule.num_atoms()];
        let known_atoms = vec![false; molecule.num_atoms()];
        let mut possible_bonds = vec![false; molecule.num_bonds()];
        let mut known_bonds = vec![false; molecule.num_bonds()];
        known_bonds[bond.index()] = true;
        let mut fixed_bonds = vec![false; molecule.num_bonds()];
        let mut records = Vec::new();

        assert!(
            update_bonds(
                &molecule,
                &ranks,
                &mut bond_symbols,
                &possible_atoms,
                &mut possible_bonds,
                &known_atoms,
                &known_bonds,
                &mut fixed_bonds,
                &mut records,
            )
            .unwrap()
        );
        assert!(fixed_bonds[bond.index()]);
        assert_eq!(records.len(), 1);
        let ordered = records[0]
            .controlling_atoms()
            .iter()
            .map(|controller| match controller {
                ControllingAtom::Atom(atom) => ranks[atom.index()],
                ControllingAtom::Missing => 0,
            })
            .collect::<Vec<_>>();
        assert!(ordered[0] > ordered[1]);
        assert!(ordered[2] > ordered[3]);
        assert_eq!(
            records[0].descriptor(),
            match source_info.descriptor() {
                StereoDescriptor::BondCis => StereoDescriptor::BondTrans,
                StereoDescriptor::BondTrans => StereoDescriptor::BondCis,
                descriptor => panic!("unexpected double-bond descriptor {descriptor:?}"),
            }
        );

        records.clear();
        assert!(
            !update_bonds(
                &molecule,
                &ranks,
                &mut bond_symbols,
                &possible_atoms,
                &mut possible_bonds,
                &known_atoms,
                &known_bonds,
                &mut fixed_bonds,
                &mut records,
            )
            .unwrap()
        );
        assert_eq!(records.len(), 1);
    }

    #[test]
    fn atom_duplicate_removal_rolls_back_ring_transmission_counts() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let ring_atoms = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for index in 0..ring_atoms.len() {
            builder
                .add_bond(BondSpec::new(
                    ring_atoms[index],
                    ring_atoms[(index + 1) % ring_atoms.len()],
                    BondOrder::Single,
                ))
                .unwrap();
        }
        for center in [ring_atoms[0], ring_atoms[3]] {
            for _ in 0..2 {
                let substituent = builder.add_atom(AtomSpec::new(Element::C));
                builder
                    .add_bond(BondSpec::new(center, substituent, BondOrder::Single))
                    .unwrap();
            }
        }
        let mut molecule = builder.build().unwrap();
        molecule.derived_cache_mut().rings = Some(crate::symmetrize_sssr(&molecule).unwrap());
        molecule.derived_cache_mut().valence = Some(
            crate::assign_valence_with_options(&molecule, ValenceModel::RdkitLike, false).unwrap(),
        );
        let mut possible_atoms = vec![false; molecule.num_atoms()];
        possible_atoms[0] = true;
        possible_atoms[3] = true;
        let known_atoms = vec![false; molecule.num_atoms()];
        let mut fixed_atoms = vec![false; molecule.num_atoms()];
        let mut atom_symbols = molecule
            .atoms()
            .iter()
            .map(atom_compare_symbol)
            .collect::<Vec<_>>();
        atom_symbols[0].push_str("_0");
        atom_symbols[3].push_str("_3");
        let ranks = vec![0; molecule.num_atoms()];
        let mut ring_atoms = vec![0; molecule.num_atoms()];
        ring_atoms[0] = 1;
        ring_atoms[3] = 1;
        let mut ring_bonds = vec![0; molecule.num_bonds()];
        ring_bonds[..6].fill(1);
        let mut records = Vec::new();

        assert!(
            update_atoms(
                &molecule,
                &ranks,
                &mut atom_symbols,
                &mut possible_atoms,
                &known_atoms,
                &mut fixed_atoms,
                &mut ring_atoms,
                &mut ring_bonds,
                &mut records,
                true,
            )
            .unwrap()
        );
        assert!(!possible_atoms[0]);
        assert!(!possible_atoms[3]);
        assert_eq!(ring_atoms, [0; 10]);
        assert_eq!(ring_bonds, [0; 10]);
        assert_eq!(fixed_atoms, [false; 10]);
        assert!(records.is_empty());
    }

    #[test]
    fn atom_bond_updates_match_source_records_order_and_round_count() {
        for (smiles, expected_atoms, expected_bonds, expected_rounds) in [
            ("CC(F)C(C)C(C)F", 3, 0, 1),
            ("CC=C1CCC(O)CC1", 1, 1, 1),
            ("CC=CC(C=CC)C(C)C(C=CC)C=CC", 3, 4, 1),
        ] {
            let mut molecule = Molecule::from_smiles(smiles).unwrap();
            let mut state = iterative_update_state(&mut molecule);
            let mut rounds = 0;
            let final_records = loop {
                rounds += 1;
                let (another_round, records, _) = run_iterative_update_round(&molecule, &mut state);
                if !another_round {
                    break records;
                }
                assert!(rounds < molecule.num_atoms() + molecule.num_bonds() + 2);
            };
            let atom_count = final_records
                .iter()
                .filter(|record| matches!(record.center(), StereoCenter::Atom(_)))
                .count();
            let bond_count = final_records
                .iter()
                .filter(|record| matches!(record.center(), StereoCenter::Bond(_)))
                .count();
            assert_eq!((atom_count, bond_count), (expected_atoms, expected_bonds));
            assert_eq!(rounds, expected_rounds, "{smiles} refinement rounds");
            assert!(
                final_records
                    .windows(2)
                    .all(|pair| pair[0].center() <= pair[1].center()),
                "{smiles} records must retain atom-then-bond source order"
            );
        }
    }

    #[test]
    fn ring_stereo_divisor_two_three_and_odd_ring_rules_match_source() {
        let mut even = ring_test_molecule(6);
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut even, &[0, 3]);
        assert_eq!(atom_counts, [1, 0, 0, 1, 0, 0]);
        assert_eq!(bond_counts, [1; 6]);
        assert_eq!(even.atoms()[0].prop("_ringStereoOtherAtom"), Some("3"));
        assert_eq!(even.atoms()[3].prop("_ringStereoOtherAtom"), Some("0"));
        assert!(even.derived_cache().rings.is_some());

        let mut thirds = ring_test_molecule(6);
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut thirds, &[0, 2, 4]);
        assert_eq!(atom_counts, [1, 0, 1, 0, 1, 0]);
        assert_eq!(bond_counts, [1; 6]);
        assert!(
            thirds
                .atoms()
                .iter()
                .all(|atom| atom.prop("_ringStereoOtherAtom").is_none())
        );

        let mut odd = ring_test_molecule(5);
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut odd, &[0, 2]);
        assert_eq!(atom_counts, [0; 5]);
        assert_eq!(bond_counts, [0; 5]);
    }

    #[test]
    fn fused_common_edge_and_spiro_ring_transmission_are_kept_separate() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let atoms = (0..10)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for (begin, end) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ] {
            builder
                .add_bond(BondSpec::new(atoms[begin], atoms[end], BondOrder::Single))
                .unwrap();
        }
        let mut fused = builder.build().unwrap();
        fused.derived_cache_mut().rings = Some(crate::symmetrize_sssr(&fused).unwrap());
        let shared_bond = bond_between_atoms(&fused, atoms[0], atoms[1]).unwrap();
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut fused, &[0, 1]);
        assert_eq!(atom_counts[0], 2);
        assert_eq!(atom_counts[1], 2);
        assert_eq!(bond_counts[shared_bond.index()], 2);
        assert_eq!(fused.atoms()[0].prop("_ringStereoOtherAtom"), Some("1"));
        assert_eq!(fused.atoms()[1].prop("_ringStereoOtherAtom"), Some("0"));

        let mut builder = MoleculeBuilder::new();
        let atoms = (0..7)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for (begin, end) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (0, 4),
            (4, 5),
            (5, 6),
            (6, 0),
        ] {
            builder
                .add_bond(BondSpec::new(atoms[begin], atoms[end], BondOrder::Single))
                .unwrap();
        }
        let mut spiro = builder.build().unwrap();
        spiro.derived_cache_mut().rings = Some(crate::symmetrize_sssr(&spiro).unwrap());
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut spiro, &[0, 2]);
        assert_eq!(atom_counts, [1, 0, 1, 0, 0, 0, 0]);
        for bond in [
            bond_between_atoms(&spiro, atoms[0], atoms[1]).unwrap(),
            bond_between_atoms(&spiro, atoms[1], atoms[2]).unwrap(),
            bond_between_atoms(&spiro, atoms[2], atoms[3]).unwrap(),
            bond_between_atoms(&spiro, atoms[3], atoms[0]).unwrap(),
        ] {
            assert_eq!(bond_counts[bond.index()], 1);
        }
        for bond in [
            bond_between_atoms(&spiro, atoms[0], atoms[4]).unwrap(),
            bond_between_atoms(&spiro, atoms[4], atoms[5]).unwrap(),
            bond_between_atoms(&spiro, atoms[5], atoms[6]).unwrap(),
            bond_between_atoms(&spiro, atoms[6], atoms[0]).unwrap(),
        ] {
            assert_eq!(bond_counts[bond.index()], 0);
        }
    }

    #[test]
    fn cage_and_adjacent_center_regression_follow_source_ring_order() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut adjacent = ring_test_molecule(6);
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut adjacent, &[0, 1]);
        assert_eq!(atom_counts, [0; 6]);
        assert_eq!(bond_counts, [0; 6]);

        let mut builder = MoleculeBuilder::new();
        let atoms = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for begin in 0..4 {
            for end in begin + 1..4 {
                builder
                    .add_bond(BondSpec::new(atoms[begin], atoms[end], BondOrder::Single))
                    .unwrap();
            }
        }
        let mut cage = builder.build().unwrap();
        cage.derived_cache_mut().rings = Some(crate::symmetrize_sssr(&cage).unwrap());
        let (atom_counts, bond_counts) = ring_stereo_counts(&mut cage, &[0, 1, 2, 3]);
        assert!(atom_counts.iter().all(|count| *count > 0));
        assert!(bond_counts.iter().all(|count| *count > 0));
    }

    #[test]
    fn controlling_atom_duplicate_ties_use_opposite_atom_and_atrop_bond_state() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let molecule = ring_test_molecule(6);
        let bond = bond_between_atoms(&molecule, AtomId::new(0), AtomId::new(1)).unwrap();
        let ranks = vec![0, 1, 7, 2, 3, 7];
        let empty_atoms = vec![false; molecule.num_atoms()];
        let empty_bonds = vec![false; molecule.num_bonds()];
        assert!(
            are_stereobond_controlling_atoms_dupes(
                &molecule,
                bond,
                ControllingAtom::Atom(AtomId::new(2)),
                ControllingAtom::Atom(AtomId::new(5)),
                &ranks,
                &empty_atoms,
                &empty_atoms,
                &empty_bonds,
                &empty_bonds,
            )
            .unwrap()
        );
        let mut possible_atoms = empty_atoms.clone();
        possible_atoms[3] = true;
        assert!(
            !are_stereobond_controlling_atoms_dupes(
                &molecule,
                bond,
                ControllingAtom::Atom(AtomId::new(2)),
                ControllingAtom::Atom(AtomId::new(5)),
                &ranks,
                &possible_atoms,
                &empty_atoms,
                &empty_bonds,
                &empty_bonds,
            )
            .unwrap()
        );

        let mut builder = MoleculeBuilder::new();
        let atoms = (0..7)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        for index in 0..6 {
            builder
                .add_bond(BondSpec::new(
                    atoms[index],
                    atoms[(index + 1) % 6],
                    BondOrder::Single,
                ))
                .unwrap();
        }
        let external = builder
            .add_bond(BondSpec::new(atoms[3], atoms[6], BondOrder::Single))
            .unwrap();
        let mut atrop_tie = builder.build().unwrap();
        atrop_tie.derived_cache_mut().rings = Some(crate::symmetrize_sssr(&atrop_tie).unwrap());
        let bond = bond_between_atoms(&atrop_tie, atoms[0], atoms[1]).unwrap();
        let mut possible_bonds = vec![false; atrop_tie.num_bonds()];
        possible_bonds[external.index()] = true;
        assert!(
            !are_stereobond_controlling_atoms_dupes(
                &atrop_tie,
                bond,
                ControllingAtom::Atom(atoms[2]),
                ControllingAtom::Atom(atoms[5]),
                &[0, 1, 7, 2, 3, 7, 4],
                &vec![false; atrop_tie.num_atoms()],
                &vec![false; atrop_tie.num_atoms()],
                &possible_bonds,
                &vec![false; atrop_tie.num_bonds()],
            )
            .unwrap()
        );
    }

    #[test]
    fn controlling_atom_duplicate_validation_and_odd_ring_behavior_are_structured() {
        let molecule = ring_test_molecule(5);
        let bond = bond_between_atoms(&molecule, AtomId::new(0), AtomId::new(1)).unwrap();
        let ranks = vec![0, 1, 7, 2, 7];
        let atoms = vec![false; molecule.num_atoms()];
        let bonds = vec![false; molecule.num_bonds()];
        assert!(
            are_stereobond_controlling_atoms_dupes(
                &molecule,
                bond,
                ControllingAtom::Atom(AtomId::new(2)),
                ControllingAtom::Atom(AtomId::new(4)),
                &ranks,
                &atoms,
                &atoms,
                &bonds,
                &bonds,
            )
            .unwrap()
        );
        assert_eq!(
            are_stereobond_controlling_atoms_dupes(
                &molecule,
                bond,
                ControllingAtom::Missing,
                ControllingAtom::Atom(AtomId::new(4)),
                &ranks,
                &atoms,
                &atoms,
                &bonds,
                &bonds,
            ),
            Err(PotentialStereoError::MissingControllingAtom)
        );
        assert_eq!(
            are_stereobond_controlling_atoms_dupes(
                &molecule,
                bond,
                ControllingAtom::Atom(AtomId::new(2)),
                ControllingAtom::Atom(AtomId::new(4)),
                &ranks[..4],
                &atoms,
                &atoms,
                &bonds,
                &bonds,
            ),
            Err(PotentialStereoError::InvalidStateTableLength {
                table: "atom_ranks",
                expected: 5,
                actual: 4,
            })
        );
    }

    #[test]
    fn atom_eligibility_matches_rdkit_element_degree_and_valence_cases() {
        let cases = [
            ("CC(F)(Cl)CNC(C)(C)C", 1, true),
            ("CC(F)(Cl)CNC(C)(C)C", 0, false),
            ("CC(F)(Cl)CNC(C)(C)C", 4, false),
            ("CC(F)(Cl)CNC(C)(C)C", 6, true),
            ("O=S(F)CC[S+]([O-])CS=O", 1, true),
            ("O=S(F)CC[S+]([O-])CS=O", 5, true),
            ("O=S(F)CC[S+]([O-])CS=O", 8, false),
            ("O=[Se](F)CC[Se+]([O-])C[Se]=O", 1, true),
            ("O=[Se](F)CC[Se+]([O-])C[Se]=O", 5, true),
            ("O=[Se](F)CC[Se+]([O-])C[Se]=O", 8, false),
            ("OP(F)CPCP", 1, true),
            ("OP(F)CPCP", 4, true),
            ("OP(F)CPCP", 6, false),
            ("O[As](F)C[As]C[As]", 1, true),
            ("O[As](F)C[As]C[As]", 4, true),
            ("O[As](F)C[As]C[As]", 6, false),
            ("O[P]([O-])(=O)OC", 1, true),
        ];

        for (smiles, atom, expected) in cases {
            let molecule = Molecule::from_smiles(smiles).unwrap();
            assert_eq!(
                is_atom_potential_tetrahedral_center(&molecule, AtomId::new(atom)),
                Ok(expected),
                "unexpected tetrahedral eligibility for atom {atom} of {smiles}"
            );
        }
    }

    #[test]
    fn atom_eligibility_distinguishes_implicit_protium_and_isotopic_hydrogen() {
        let implicit = Molecule::from_smiles("C(F)(Cl)Br").unwrap();
        assert_eq!(atom_total_hydrogens(&implicit, AtomId::new(0)), Ok(1));
        assert!(!has_protium_neighbor(&implicit, AtomId::new(0)).unwrap());
        assert!(is_atom_potential_tetrahedral_center(&implicit, AtomId::new(0)).unwrap());

        let keep_hs = crate::SmilesParseParams::default().with_remove_hs(false);
        let protium = Molecule::from_smiles_with_params("[C@H]([H])(F)Cl", &keep_hs).unwrap();
        assert!(has_protium_neighbor(&protium, AtomId::new(0)).unwrap());
        assert!(!is_atom_potential_tetrahedral_center(&protium, AtomId::new(0)).unwrap());

        let deuterium = Molecule::from_smiles_with_params("[C@H]([2H])(F)Cl", &keep_hs).unwrap();
        assert!(!has_protium_neighbor(&deuterium, AtomId::new(0)).unwrap());
        assert!(is_atom_potential_tetrahedral_center(&deuterium, AtomId::new(0)).unwrap());
    }

    #[test]
    fn atom_nonzero_degree_uses_source_dative_direction_and_zero_bond_rules() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let donor = builder.add_atom(AtomSpec::new(Element::N));
        let acceptor = builder.add_atom(AtomSpec::new(Element::B));
        let zero_neighbor = builder.add_atom(AtomSpec::new(Element::C));
        let ordinary_neighbor = builder.add_atom(AtomSpec::new(Element::F));
        builder
            .add_bond(BondSpec::new(donor, acceptor, BondOrder::Dative))
            .unwrap();
        builder
            .add_bond(BondSpec::new(acceptor, zero_neighbor, BondOrder::Zero))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                acceptor,
                ordinary_neighbor,
                BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();

        assert_eq!(get_atom_nonzero_degree(&molecule, donor), Ok(0));
        assert_eq!(get_atom_nonzero_degree(&molecule, acceptor), Ok(2));
        assert_eq!(get_atom_nonzero_degree(&molecule, zero_neighbor), Ok(0));
        assert_eq!(get_atom_nonzero_degree(&molecule, ordinary_neighbor), Ok(1));
    }

    #[test]
    fn nitrogen_eligibility_requires_ring_or_bridgehead_sp3_nonconjugated_state() {
        let ring_n = Molecule::from_smiles("CN1CC1N(F)C").unwrap();
        assert!(is_atom_potential_tetrahedral_center(&ring_n, AtomId::new(1)).unwrap());
        assert!(!is_atom_potential_tetrahedral_center(&ring_n, AtomId::new(4)).unwrap());

        let bridgehead = Molecule::from_smiles("C1CCCCC[C@@H]2CN1CCO2").unwrap();
        assert!(is_atom_potential_tetrahedral_center(&bridgehead, AtomId::new(6)).unwrap());
        assert!(is_atom_potential_tetrahedral_center(&bridgehead, AtomId::new(8)).unwrap());

        assert_eq!(bridgehead.atoms()[8].hybridization(), Hybridization::Sp3);
        assert!(
            !bridgehead
                .topology_block()
                .adjacency
                .neighbors_of(8)
                .iter()
                .any(|neighbor| bridgehead.bonds()[neighbor.bond.index()].is_conjugated())
        );

        let conjugated = Molecule::from_smiles("O=C1CCCCC[C@@H]2CN1CCO2").unwrap();
        assert!(!is_atom_potential_tetrahedral_center(&conjugated, AtomId::new(9)).unwrap());
        assert!(
            conjugated
                .topology_block()
                .adjacency
                .neighbors_of(9)
                .iter()
                .any(|neighbor| conjugated.bonds()[neighbor.bond.index()].is_conjugated())
        );
    }

    #[test]
    fn nontetrahedral_eligibility_covers_tags_elements_and_degree_boundaries() {
        let cases = [
            ("C[S+](O)F", 1, false),
            ("C[SH](O)F", 1, true),
            ("C[S@SP](O)F", 1, true),
            ("C[Pt@SP1](F)(Cl)O", 1, true),
            ("C[Co@TB1](F)(Cl)(Br)O", 1, true),
            ("C[Co@OH1](F)(Cl)(Br)(I)O", 1, true),
            ("C[Pt@SP1]F", 1, true),
            ("[Pt@SP1]F", 0, false),
            ("C(F)(Cl)Br", 0, false),
        ];

        for (smiles, atom, expected) in cases {
            let molecule = Molecule::from_smiles(smiles).unwrap();
            assert_eq!(
                is_atom_potential_nontetrahedral_center(&molecule, AtomId::new(atom)),
                Ok(expected),
                "unexpected non-tetrahedral eligibility for atom {atom} of {smiles}"
            );
        }

        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder, ValenceAssignment};
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::BE));
        for element in [Element::F, Element::CL, Element::BR, Element::I] {
            let neighbor = builder.add_atom(AtomSpec::new(element));
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .unwrap();
        }
        let mut beryllium = builder.build().unwrap();
        beryllium.derived_cache_mut().valence = Some(ValenceAssignment {
            explicit_valence: vec![4, 1, 1, 1, 1],
            implicit_hydrogens: vec![0; 5],
        });
        assert!(is_atom_potential_nontetrahedral_center(&beryllium, center).unwrap());
    }

    #[test]
    fn combined_atom_eligibility_obeys_the_nontetrahedral_feature_switch() {
        let molecule = Molecule::from_smiles("C[Co@TB1](F)(Cl)(Br)O").unwrap();
        let center = AtomId::new(1);
        assert!(!is_atom_potential_tetrahedral_center(&molecule, center).unwrap());
        assert!(is_atom_potential_nontetrahedral_center(&molecule, center).unwrap());
        assert!(!is_atom_potential_stereo(&molecule, center, false).unwrap());
        assert!(is_atom_potential_stereo(&molecule, center, true).unwrap());
    }

    #[test]
    fn query_atoms_follow_the_same_source_eligibility_rules() {
        use crate::{
            AtomQueryPredicate, AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder, QueryNode,
            ValenceAssignment,
        };

        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_query(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))),
        );
        for element in [Element::F, Element::CL, Element::BR, Element::I] {
            let neighbor = builder.add_atom(AtomSpec::new(element));
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .unwrap();
        }
        let mut molecule = builder.build().unwrap();
        molecule.derived_cache_mut().valence = Some(ValenceAssignment {
            explicit_valence: vec![4, 1, 1, 1, 1],
            implicit_hydrogens: vec![0; 5],
        });

        assert!(molecule.atoms()[center.index()].query().is_some());
        assert!(is_atom_potential_tetrahedral_center(&molecule, center).unwrap());
    }

    #[test]
    fn atom_eligibility_reports_invalid_indices_and_unprepared_cache_state() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        assert_eq!(
            get_atom_nonzero_degree(&molecule, AtomId::new(2)),
            Err(PotentialStereoError::AtomIndexOutOfBounds {
                atom: AtomId::new(2)
            })
        );
        assert_eq!(
            has_protium_neighbor(&molecule, AtomId::new(2)),
            Err(PotentialStereoError::AtomIndexOutOfBounds {
                atom: AtomId::new(2)
            })
        );

        let unprepared = Molecule::from_smiles_with_sanitize("C(F)(Cl)Br", false).unwrap();
        assert_eq!(
            is_atom_potential_tetrahedral_center(&unprepared, AtomId::new(0)),
            Err(PotentialStereoError::MissingValenceState {
                atom: AtomId::new(0)
            })
        );
    }

    #[test]
    fn atom_stereo_info_matches_rdkit_tetrahedral_records_and_isotopes() {
        let molecule = Molecule::from_smiles("C[C@@H](O)[C@H](C)[C@H](C)O").unwrap();
        let expected = [
            (1, StereoDescriptor::TetrahedralClockwise, vec![0, 2, 3]),
            (
                3,
                StereoDescriptor::TetrahedralCounterclockwise,
                vec![1, 4, 5],
            ),
            (
                5,
                StereoDescriptor::TetrahedralCounterclockwise,
                vec![3, 6, 7],
            ),
        ];
        for (center, descriptor, controllers) in expected {
            let info = get_atom_stereo_info(&molecule, AtomId::new(center), true).unwrap();
            assert_eq!(info.stereo_type(), StereoType::AtomTetrahedral);
            assert_eq!(info.center(), StereoCenter::Atom(AtomId::new(center)));
            assert_eq!(info.specified(), StereoSpecified::Specified);
            assert_eq!(info.descriptor(), descriptor);
            assert_eq!(
                info.controlling_atoms(),
                controllers
                    .into_iter()
                    .map(|index| ControllingAtom::Atom(AtomId::new(index)))
                    .collect::<Vec<_>>()
            );
        }

        let keep_hs = crate::SmilesParseParams::default().with_remove_hs(false);
        let isotopic =
            Molecule::from_smiles_with_params("[C@](F)(Cl)([2H])[3H]", &keep_hs).unwrap();
        let info = get_atom_stereo_info(&isotopic, AtomId::new(0), true).unwrap();
        assert_eq!(
            info.descriptor(),
            StereoDescriptor::TetrahedralCounterclockwise
        );
        assert_eq!(
            info.controlling_atoms(),
            [1, 2, 3, 4].map(|index| ControllingAtom::Atom(AtomId::new(index)))
        );
    }

    #[test]
    fn atom_stereo_info_normalizes_adjacency_order_with_swap_parity() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        fn build(order: [usize; 4]) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            let center = builder
                .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
            let neighbors = [Element::F, Element::CL, Element::BR, Element::I]
                .map(|element| builder.add_atom(AtomSpec::new(element)));
            for index in order {
                builder
                    .add_bond(BondSpec::new(center, neighbors[index], BondOrder::Single))
                    .unwrap();
            }
            builder.build().unwrap()
        }

        let even = get_atom_stereo_info(&build([0, 1, 2, 3]), AtomId::new(0), true).unwrap();
        let odd = get_atom_stereo_info(&build([1, 0, 2, 3]), AtomId::new(0), true).unwrap();
        assert_eq!(even.controlling_atoms(), odd.controlling_atoms());
        assert_eq!(even.descriptor(), StereoDescriptor::TetrahedralClockwise);
        assert_eq!(
            odd.descriptor(),
            StereoDescriptor::TetrahedralCounterclockwise
        );
    }

    #[test]
    fn atom_stereo_info_honors_direction_and_property_unknown_state() {
        use crate::{AtomSpec, BondDirection, BondOrder, BondSpec, Element, MoleculeBuilder};

        for unknown_bond in [
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_direction(BondDirection::Unknown),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_prop("_UnknownStereo", "1"),
        ] {
            let mut builder = MoleculeBuilder::new();
            let center = builder
                .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
            let first = builder.add_atom(AtomSpec::new(Element::F));
            let second = builder.add_atom(AtomSpec::new(Element::CL));
            let third = builder.add_atom(AtomSpec::new(Element::BR));
            assert_eq!((center, first), (unknown_bond.begin(), unknown_bond.end()));
            builder.add_bond(unknown_bond).unwrap();
            builder
                .add_bond(BondSpec::new(center, second, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(center, third, BondOrder::Single))
                .unwrap();
            let molecule = builder.build().unwrap();
            let info = get_atom_stereo_info(&molecule, center, true).unwrap();
            assert_eq!(info.specified(), StereoSpecified::Unknown);
            assert_eq!(info.descriptor(), StereoDescriptor::None);
        }
    }

    #[test]
    fn atom_stereo_info_infers_unspecified_nontetrahedral_types_by_total_degree() {
        let cases = [
            ("C[Pt](F)(Cl)O", StereoType::AtomTetrahedral),
            ("C[Co](F)(Cl)(Br)O", StereoType::AtomTrigonalBipyramidal),
            ("C[Co](F)(Cl)(Br)(I)O", StereoType::AtomOctahedral),
        ];
        for (smiles, stereo_type) in cases {
            let molecule = Molecule::from_smiles(smiles).unwrap();
            let info = get_atom_stereo_info(&molecule, AtomId::new(1), true).unwrap();
            assert_eq!(info.stereo_type(), stereo_type, "{smiles}");
            assert_eq!(info.specified(), StereoSpecified::Unspecified, "{smiles}");
            assert_eq!(info.permutation(), 0, "{smiles}");
        }
    }

    #[test]
    fn atom_stereo_info_preserves_nontetrahedral_permutations_and_missing_ligands() {
        let cases = [
            ("C[Pt@SP1](F)Cl", StereoType::AtomSquarePlanar, 2, 3),
            (
                "C[Co@TB1](F)(Cl)Br",
                StereoType::AtomTrigonalBipyramidal,
                3,
                4,
            ),
            ("C[Co@OH1](F)(Cl)(Br)I", StereoType::AtomOctahedral, 3, 5),
        ];
        for (smiles, stereo_type, permutation, controller_count) in cases {
            let molecule = Molecule::from_smiles(smiles).unwrap();
            let info = get_atom_stereo_info(&molecule, AtomId::new(1), true).unwrap();
            assert_eq!(info.stereo_type(), stereo_type, "{smiles}");
            assert_eq!(info.specified(), StereoSpecified::Specified, "{smiles}");
            assert_eq!(info.permutation(), permutation, "{smiles}");
            assert_eq!(info.controlling_atoms().len(), controller_count, "{smiles}");
            assert!(!info.controlling_atoms().contains(&ControllingAtom::Missing));
        }

        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder, ValenceAssignment};
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::PT)
                .with_chiral_tag(ChiralTag::SquarePlanar)
                .with_chiral_permutation(0),
        );
        for element in [Element::F, Element::CL, Element::BR, Element::I] {
            let neighbor = builder.add_atom(AtomSpec::new(element));
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .unwrap();
        }
        let mut molecule = builder.build().unwrap();
        molecule.derived_cache_mut().valence = Some(ValenceAssignment {
            explicit_valence: vec![4, 1, 1, 1, 1],
            implicit_hydrogens: vec![0; 5],
        });
        let info = get_atom_stereo_info(&molecule, center, true).unwrap();
        assert_eq!(info.stereo_type(), StereoType::AtomSquarePlanar);
        assert_eq!(info.specified(), StereoSpecified::Unknown);
        assert_eq!(info.permutation(), 0);
    }

    #[test]
    fn atom_stereo_info_rejects_malformed_unknown_property_and_invalid_index() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let neighbor = builder.add_atom(AtomSpec::new(Element::F));
        let bond = builder
            .add_bond(
                BondSpec::new(center, neighbor, BondOrder::Single)
                    .with_prop("_UnknownStereo", "not-an-integer"),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        assert_eq!(
            get_atom_stereo_info(&molecule, center, true),
            Err(PotentialStereoError::InvalidBondIntegerProperty {
                bond,
                property: "_UnknownStereo",
                value: "not-an-integer".to_string(),
            })
        );
        assert_eq!(
            get_atom_stereo_info(&molecule, AtomId::new(2), true),
            Err(PotentialStereoError::AtomIndexOutOfBounds {
                atom: AtomId::new(2),
            })
        );
    }

    fn substituted_double_bond(
        stereo: crate::BondStereo,
        stereo_atoms: Option<[usize; 2]>,
        double_direction: crate::BondDirection,
        first_neighbor_direction: crate::BondDirection,
        begin_spec: crate::AtomSpec,
    ) -> Molecule {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let begin = builder.add_atom(begin_spec);
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let begin_first = builder.add_atom(AtomSpec::new(Element::F));
        let begin_second = builder.add_atom(AtomSpec::new(Element::CL));
        let end_first = builder.add_atom(AtomSpec::new(Element::BR));
        let end_second = builder.add_atom(AtomSpec::new(Element::I));
        let mut double = BondSpec::new(begin, end, BondOrder::Double)
            .with_stereo(stereo)
            .with_direction(double_direction);
        if let Some([begin_ref, end_ref]) = stereo_atoms {
            double = double.with_stereo_atoms(AtomId::new(begin_ref), AtomId::new(end_ref));
        }
        builder.add_bond(double).unwrap();
        builder
            .add_bond(
                BondSpec::new(begin, begin_first, BondOrder::Single)
                    .with_direction(first_neighbor_direction),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, begin_second, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, end_first, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, end_second, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn bond_eligibility_matches_chain_ring_and_cumulene_boundaries() {
        let chain = Molecule::from_smiles("CC=CC").unwrap();
        assert_eq!(is_bond_potential_stereo(&chain, BondId::new(1)), Ok(true));

        let seven_membered = Molecule::from_smiles("C1=CCCCCC1").unwrap();
        assert_eq!(
            is_bond_potential_stereo(&seven_membered, BondId::new(0)),
            Ok(false)
        );
        let eight_membered = Molecule::from_smiles("C1=CCCCCCC1").unwrap();
        assert_eq!(
            is_bond_potential_stereo(&eight_membered, BondId::new(0)),
            Ok(true)
        );

        let terminal_allene = Molecule::from_smiles("C=C=C").unwrap();
        assert_eq!(
            is_bond_potential_stereo(&terminal_allene, BondId::new(0)),
            Ok(false)
        );
        assert_eq!(
            is_bond_potential_stereo(&terminal_allene, BondId::new(1)),
            Ok(false)
        );
        let substituted_cumulene = Molecule::from_smiles("CC=C=CC").unwrap();
        assert_eq!(
            is_bond_potential_stereo(&substituted_cumulene, BondId::new(1)),
            Ok(true)
        );
        assert_eq!(
            is_bond_potential_stereo(&substituted_cumulene, BondId::new(2)),
            Ok(true)
        );
    }

    #[test]
    fn bond_stereo_info_preserves_endpoint_controller_order_and_padding() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder, ValenceAssignment};

        let four = substituted_double_bond(
            crate::BondStereo::None,
            None,
            crate::BondDirection::None,
            crate::BondDirection::None,
            AtomSpec::new(Element::C),
        );
        let info = get_bond_stereo_info(&four, BondId::new(0)).unwrap();
        assert_eq!(info.stereo_type(), StereoType::BondDouble);
        assert_eq!(info.specified(), StereoSpecified::Unspecified);
        assert_eq!(
            info.controlling_atoms(),
            [2, 3, 4, 5].map(|index| ControllingAtom::Atom(AtomId::new(index)))
        );

        let mut builder = MoleculeBuilder::new();
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let begin_ref = builder.add_atom(AtomSpec::new(Element::F));
        let end_ref = builder.add_atom(AtomSpec::new(Element::CL));
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, begin_ref, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, end_ref, BondOrder::Single))
            .unwrap();
        let mut padded = builder.build().unwrap();
        padded.derived_cache_mut().valence = Some(ValenceAssignment {
            explicit_valence: vec![3, 3, 1, 1],
            implicit_hydrogens: vec![1, 1, 0, 0],
        });
        padded.derived_cache_mut().rings =
            Some(crate::RingInfo::new(crate::RingFindType::SymmSssr, 4, 3));
        assert_eq!(is_bond_potential_stereo(&padded, BondId::new(0)), Ok(true));
        let info = get_bond_stereo_info(&padded, BondId::new(0)).unwrap();
        assert_eq!(
            info.controlling_atoms(),
            [
                ControllingAtom::Atom(begin_ref),
                ControllingAtom::Missing,
                ControllingAtom::Atom(end_ref),
                ControllingAtom::Missing,
            ]
        );
    }

    #[test]
    fn bond_stereo_info_unknown_state_has_source_precedence() {
        use crate::{AtomSpec, Element};

        let cases = [
            substituted_double_bond(
                crate::BondStereo::Any,
                None,
                crate::BondDirection::None,
                crate::BondDirection::None,
                AtomSpec::new(Element::C),
            ),
            substituted_double_bond(
                crate::BondStereo::None,
                None,
                crate::BondDirection::EitherDouble,
                crate::BondDirection::None,
                AtomSpec::new(Element::C),
            ),
            substituted_double_bond(
                crate::BondStereo::None,
                None,
                crate::BondDirection::None,
                crate::BondDirection::Unknown,
                AtomSpec::new(Element::C),
            ),
            substituted_double_bond(
                crate::BondStereo::None,
                None,
                crate::BondDirection::None,
                crate::BondDirection::None,
                AtomSpec::new(Element::C).with_prop("_UnknownStereo", "1"),
            ),
            substituted_double_bond(
                crate::BondStereo::None,
                None,
                crate::BondDirection::None,
                crate::BondDirection::None,
                AtomSpec::new(Element::C).with_unknown_stereo(true),
            ),
        ];
        for molecule in cases {
            let info = get_bond_stereo_info(&molecule, BondId::new(0)).unwrap();
            assert_eq!(info.specified(), StereoSpecified::Unknown);
            assert_eq!(info.descriptor(), StereoDescriptor::None);
        }
    }

    #[test]
    fn bond_stereo_info_normalizes_ez_and_stereo_atom_positions() {
        use crate::{AtomSpec, Element};

        let cases = [
            (crate::BondStereo::E, [2, 4], StereoDescriptor::BondTrans),
            (crate::BondStereo::Z, [2, 4], StereoDescriptor::BondCis),
            (crate::BondStereo::Cis, [2, 5], StereoDescriptor::BondTrans),
            (
                crate::BondStereo::Trans,
                [3, 5],
                StereoDescriptor::BondTrans,
            ),
            (crate::BondStereo::Trans, [3, 4], StereoDescriptor::BondCis),
        ];
        for (stereo, stereo_atoms, expected) in cases {
            let molecule = substituted_double_bond(
                stereo,
                Some(stereo_atoms),
                crate::BondDirection::None,
                crate::BondDirection::None,
                AtomSpec::new(Element::C),
            );
            let info = get_bond_stereo_info(&molecule, BondId::new(0)).unwrap();
            assert_eq!(info.specified(), StereoSpecified::Specified);
            assert_eq!(info.descriptor(), expected);
        }
    }

    #[test]
    fn represented_atropisomer_bond_uses_ordered_controller_slots() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        for (stereo, descriptor) in [
            (
                crate::BondStereo::AtropCw,
                StereoDescriptor::BondAtropClockwise,
            ),
            (
                crate::BondStereo::AtropCcw,
                StereoDescriptor::BondAtropCounterclockwise,
            ),
        ] {
            let mut builder = MoleculeBuilder::new();
            let begin = builder.add_atom(AtomSpec::new(Element::C));
            let end = builder.add_atom(AtomSpec::new(Element::C));
            let begin_ref = builder.add_atom(AtomSpec::new(Element::F));
            let end_ref = builder.add_atom(AtomSpec::new(Element::CL));
            builder
                .add_bond(BondSpec::new(begin, end, BondOrder::Single).with_stereo(stereo))
                .unwrap();
            builder
                .add_bond(BondSpec::new(begin, begin_ref, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(end, end_ref, BondOrder::Single))
                .unwrap();
            let molecule = builder.build().unwrap();
            let info = get_bond_stereo_info(&molecule, BondId::new(0)).unwrap();
            assert_eq!(info.stereo_type(), StereoType::BondAtropisomer);
            assert_eq!(info.specified(), StereoSpecified::Specified);
            assert_eq!(info.descriptor(), descriptor);
            assert_eq!(
                info.controlling_atoms(),
                [
                    ControllingAtom::Atom(begin_ref),
                    ControllingAtom::Missing,
                    ControllingAtom::Atom(end_ref),
                    ControllingAtom::Missing,
                ]
            );
        }
    }

    #[test]
    fn bond_stereo_info_reports_source_boundary_errors_structurally() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder, ValenceAssignment};

        let invalid_index = Molecule::from_smiles("CC").unwrap();
        assert_eq!(
            get_bond_stereo_info(&invalid_index, BondId::new(1)),
            Err(PotentialStereoError::BondIndexOutOfBounds {
                bond: BondId::new(1)
            })
        );

        let missing_stereo_atoms = substituted_double_bond(
            crate::BondStereo::E,
            None,
            crate::BondDirection::None,
            crate::BondDirection::None,
            AtomSpec::new(Element::C),
        );
        assert_eq!(
            get_bond_stereo_info(&missing_stereo_atoms, BondId::new(0)),
            Err(PotentialStereoError::InvalidStereoAtomCount {
                bond: BondId::new(0)
            })
        );

        let mismatched = substituted_double_bond(
            crate::BondStereo::E,
            Some([4, 2]),
            crate::BondDirection::None,
            crate::BondDirection::None,
            AtomSpec::new(Element::C),
        );
        assert_eq!(
            get_bond_stereo_info(&mismatched, BondId::new(0)),
            Err(PotentialStereoError::StereoAtomMismatch {
                bond: BondId::new(0),
                stereo_atom: AtomId::new(4),
                end: "begin",
            })
        );

        let malformed_property = substituted_double_bond(
            crate::BondStereo::None,
            None,
            crate::BondDirection::None,
            crate::BondDirection::None,
            AtomSpec::new(Element::C).with_prop("_UnknownStereo", "bad"),
        );
        assert_eq!(
            get_bond_stereo_info(&malformed_property, BondId::new(0)),
            Err(PotentialStereoError::InvalidAtomIntegerProperty {
                atom: AtomId::new(0),
                property: "_UnknownStereo",
                value: "bad".to_string(),
            })
        );

        let mut builder = MoleculeBuilder::new();
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let one = builder.add_atom(AtomSpec::new(Element::F));
        let two = builder.add_atom(AtomSpec::new(Element::CL));
        let three = builder.add_atom(AtomSpec::new(Element::BR));
        let end_ref = builder.add_atom(AtomSpec::new(Element::I));
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Double))
            .unwrap();
        for neighbor in [one, two, three] {
            builder
                .add_bond(BondSpec::new(begin, neighbor, BondOrder::Single))
                .unwrap();
        }
        builder
            .add_bond(BondSpec::new(end, end_ref, BondOrder::Single))
            .unwrap();
        let invalid_degree = builder.build().unwrap();
        assert_eq!(
            get_bond_stereo_info(&invalid_degree, BondId::new(0)),
            Err(PotentialStereoError::InvalidBondDegree {
                bond: BondId::new(0)
            })
        );

        let mut builder = MoleculeBuilder::new();
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(
                BondSpec::new(begin, end, BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCw),
            )
            .unwrap();
        let invalid_atrop_degree = builder.build().unwrap();
        assert_eq!(
            get_bond_stereo_info(&invalid_atrop_degree, BondId::new(0)),
            Err(PotentialStereoError::InvalidBondDegree {
                bond: BondId::new(0)
            })
        );

        let mut no_rings = Molecule::from_smiles_with_sanitize("CC=CC", false).unwrap();
        no_rings.derived_cache_mut().valence = Some(ValenceAssignment {
            explicit_valence: vec![1, 3, 3, 1],
            implicit_hydrogens: vec![3, 1, 1, 3],
        });
        assert_eq!(
            is_bond_potential_stereo(&no_rings, BondId::new(1)),
            Err(PotentialStereoError::MissingBondRingState {
                bond: BondId::new(1)
            })
        );
    }

    #[test]
    fn initialization_symbols_and_bits_cover_all_specified_states() {
        let mut specified = Molecule::from_smiles("F[C@](Cl)(Br)I.C/C=C/C").unwrap();
        let state = initialize_potential_stereo(&mut specified, true, false, true).unwrap();
        assert!(state.known_atoms[1]);
        assert!(!state.possible_atoms[1]);
        assert!(matches!(
            state.atom_symbols[1].as_str(),
            "0C0_CW" | "0C0_CCW"
        ));
        let double_bond = specified
            .bonds()
            .iter()
            .find(|bond| bond.order() == crate::BondOrder::Double)
            .unwrap()
            .id();
        assert!(state.known_bonds[double_bond.index()]);
        assert!(matches!(
            state.bond_symbols[double_bond.index()].as_str(),
            "=_cis" | "=_trans"
        ));

        use crate::{AtomSpec, BondDirection, BondOrder, BondSpec, Element, MoleculeBuilder};
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_unknown_stereo(true));
        for (neighbor_index, element) in [Element::F, Element::CL, Element::BR]
            .into_iter()
            .enumerate()
        {
            let neighbor = builder.add_atom(AtomSpec::new(element));
            builder
                .add_bond(
                    BondSpec::new(center, neighbor, BondOrder::Single).with_direction(
                        if neighbor_index == 0 {
                            BondDirection::Unknown
                        } else {
                            BondDirection::None
                        },
                    ),
                )
                .unwrap();
        }
        let mut unknown_atom = builder.build().unwrap();
        let state = initialize_potential_stereo(&mut unknown_atom, true, false, true).unwrap();
        assert!(state.known_atoms[center.index()]);
        assert_eq!(state.atom_symbols[center.index()], "0C00");

        let mut unknown_bond = substituted_double_bond(
            crate::BondStereo::Any,
            None,
            BondDirection::None,
            BondDirection::None,
            AtomSpec::new(Element::C),
        );
        unknown_bond.derived_cache_mut().rings = Some(crate::RingInfo::new(
            crate::RingFindType::SymmSssr,
            unknown_bond.num_atoms(),
            unknown_bond.num_bonds(),
        ));
        let state = initialize_potential_stereo(&mut unknown_bond, true, false, true).unwrap();
        assert!(state.known_bonds[0]);
        assert_eq!(state.bond_symbols[0], "=_0");
    }

    #[test]
    fn initialization_flag_possible_and_clean_matrix_matches_source() {
        for (flag_possible, clean, possible, decorated) in [
            (false, false, false, false),
            (false, true, false, false),
            (true, false, true, true),
            (true, true, true, false),
        ] {
            let mut molecule = Molecule::from_smiles("FC(Cl)Br.CC=CC").unwrap();
            let state = initialize_potential_stereo(&mut molecule, flag_possible, clean, true)
                .expect("initialization matrix row");
            assert_eq!(state.possible_atoms[1], possible);
            assert_eq!(state.atom_symbols[1].contains("_1"), decorated);
            let double_bond = molecule
                .bonds()
                .iter()
                .find(|bond| bond.order() == crate::BondOrder::Double)
                .unwrap()
                .id();
            assert_eq!(state.possible_bonds[double_bond.index()], possible);
            assert_eq!(
                state.bond_symbols[double_bond.index()]
                    .ends_with(&format!("_{}", double_bond.index())),
                decorated
            );
        }
    }

    #[test]
    fn initialization_prepares_non_strict_cache_and_uses_exact_source_symbols() {
        let mut molecule = Molecule::from_smiles_with_sanitize("[13C-](F)(Cl)Br", false).unwrap();
        assert!(molecule.derived_cache().valence.is_none());
        let state = initialize_potential_stereo(&mut molecule, true, false, true).unwrap();
        assert!(molecule.derived_cache().valence.is_some());
        assert_eq!(state.atom_symbols[0], "13C-1");

        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};
        let mut builder = MoleculeBuilder::new();
        let atoms = [Element::C, Element::C, Element::C, Element::C]
            .map(|element| builder.add_atom(AtomSpec::new(element)));
        builder
            .add_bond(BondSpec::new(atoms[0], atoms[1], BondOrder::Aromatic).with_aromatic(false))
            .unwrap();
        builder
            .add_bond(BondSpec::new(atoms[2], atoms[3], BondOrder::Quadruple))
            .unwrap();
        let mut symbols = builder.build().unwrap();
        let state = initialize_potential_stereo(&mut symbols, false, false, true).unwrap();
        assert_eq!(state.bond_symbols, [":", "?"]);
    }

    #[test]
    fn initialization_cleanup_clears_invalid_marks_but_preserves_atrop() {
        use crate::{AtomSpec, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let terminal =
            builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
        let terminal_neighbor = builder.add_atom(AtomSpec::new(Element::C));
        let alkene_end = builder.add_atom(AtomSpec::new(Element::C));
        let alkene_terminal = builder.add_atom(AtomSpec::new(Element::C));
        let atrop_begin = builder.add_atom(AtomSpec::new(Element::C));
        let atrop_end = builder.add_atom(AtomSpec::new(Element::C));
        let atrop_begin_ref = builder.add_atom(AtomSpec::new(Element::F));
        let atrop_end_ref = builder.add_atom(AtomSpec::new(Element::CL));
        builder
            .add_bond(BondSpec::new(
                terminal,
                terminal_neighbor,
                BondOrder::Single,
            ))
            .unwrap();
        let invalid_double = builder
            .add_bond(
                BondSpec::new(alkene_end, alkene_terminal, BondOrder::Double)
                    .with_stereo(crate::BondStereo::Any),
            )
            .unwrap();
        let atrop = builder
            .add_bond(
                BondSpec::new(atrop_begin, atrop_end, BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCcw),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                atrop_begin,
                atrop_begin_ref,
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(atrop_end, atrop_end_ref, BondOrder::Single))
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let state = initialize_potential_stereo(&mut molecule, true, true, true).unwrap();

        assert_eq!(
            molecule.atoms()[terminal.index()].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(
            molecule.bonds()[invalid_double.index()].stereo(),
            crate::BondStereo::None
        );
        assert_eq!(
            molecule.bonds()[atrop.index()].stereo(),
            crate::BondStereo::AtropCcw
        );
        assert!(state.known_bonds[atrop.index()]);
        assert_eq!(state.bond_symbols[atrop.index()], "-_atropccw");
    }

    #[test]
    fn clean_mol_stereo_clears_marked_atom_and_nontetrahedral_state_only() {
        use crate::{AtomSpec, BondDirection, BondOrder, BondSpec, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let tetra = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_prop("keep", "atom"),
        );
        let tetra_neighbor = builder.add_atom(AtomSpec::new(Element::F));
        let square = builder.add_atom(
            AtomSpec::new(Element::PT)
                .with_chiral_tag(ChiralTag::SquarePlanar)
                .with_chiral_permutation(2),
        );
        let square_neighbor = builder.add_atom(AtomSpec::new(Element::CL));
        let unrelated = builder.add_atom(AtomSpec::new(Element::C));
        let unrelated_neighbor = builder.add_atom(AtomSpec::new(Element::BR));
        let wedge = builder
            .add_bond(
                BondSpec::new(tetra, tetra_neighbor, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(square, square_neighbor, BondOrder::Single))
            .unwrap();
        let preserved = builder
            .add_bond(
                BondSpec::new(unrelated, unrelated_neighbor, BondOrder::Single)
                    .with_direction(BondDirection::BeginDash),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let mut known_atoms = vec![false; molecule.num_atoms()];
        known_atoms[tetra.index()] = true;
        known_atoms[square.index()] = true;

        clean_mol_stereo(
            &mut molecule,
            &vec![false; known_atoms.len()],
            &known_atoms,
            &vec![false; 3],
            &vec![false; 3],
        )
        .unwrap();

        assert_eq!(
            molecule.atoms()[tetra.index()].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(molecule.atoms()[tetra.index()].prop("keep"), Some("atom"));
        assert_eq!(
            molecule.bonds()[wedge.index()].direction(),
            BondDirection::None
        );
        assert_eq!(
            molecule.atoms()[square.index()].chiral_tag(),
            ChiralTag::SquarePlanar
        );
        assert_eq!(
            molecule.atoms()[square.index()].chiral_permutation(),
            Some(0)
        );
        assert_eq!(
            molecule.bonds()[preserved.index()].direction(),
            BondDirection::BeginDash
        );
    }

    #[test]
    fn clean_mol_stereo_clears_invalid_bond_and_orphaned_slashes() {
        use crate::{
            AtomSpec, BondDirection, BondOrder, BondSpec, BondStereo, Element, MoleculeBuilder,
        };

        let mut builder = MoleculeBuilder::new();
        let left = builder.add_atom(AtomSpec::new(Element::C));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        let left_ref = builder.add_atom(AtomSpec::new(Element::F));
        let right_ref = builder.add_atom(AtomSpec::new(Element::CL));
        let double = builder
            .add_bond(
                BondSpec::new(left, right, BondOrder::Double)
                    .with_stereo(BondStereo::E)
                    .with_stereo_atoms(left_ref, right_ref),
            )
            .unwrap();
        let left_slash = builder
            .add_bond(
                BondSpec::new(left, left_ref, BondOrder::Single)
                    .with_direction(BondDirection::EndUpRight),
            )
            .unwrap();
        let right_slash = builder
            .add_bond(
                BondSpec::new(right, right_ref, BondOrder::Single)
                    .with_direction(BondDirection::EndDownRight),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let atom_count = molecule.num_atoms();
        let bond_count = molecule.num_bonds();
        let mut known_bonds = vec![false; molecule.num_bonds()];
        known_bonds[double.index()] = true;

        clean_mol_stereo(
            &mut molecule,
            &vec![false; atom_count],
            &vec![false; atom_count],
            &vec![false; bond_count],
            &known_bonds,
        )
        .unwrap();

        assert_eq!(molecule.bonds()[double.index()].stereo(), BondStereo::None);
        assert_eq!(molecule.bonds()[double.index()].stereo_atoms(), None);
        assert_eq!(
            molecule.bonds()[left_slash.index()].direction(),
            BondDirection::None
        );
        assert_eq!(
            molecule.bonds()[right_slash.index()].direction(),
            BondDirection::None
        );
    }

    #[test]
    fn clean_mol_stereo_preserves_slash_after_first_endpoint_finds_stereo() {
        use crate::{
            AtomSpec, BondDirection, BondOrder, BondSpec, BondStereo, Element, MoleculeBuilder,
        };

        let mut builder = MoleculeBuilder::new();
        let x = builder.add_atom(AtomSpec::new(Element::C));
        let y = builder.add_atom(AtomSpec::new(Element::C));
        let double_other = builder.add_atom(AtomSpec::new(Element::C));
        let double_ref = builder.add_atom(AtomSpec::new(Element::F));
        let removed_left = builder.add_atom(AtomSpec::new(Element::C));
        let removed_right = builder.add_atom(AtomSpec::new(Element::C));
        let slash = builder
            .add_bond(
                BondSpec::new(x, y, BondOrder::Single).with_direction(BondDirection::EndUpRight),
            )
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(x, double_other, BondOrder::Double)
                    .with_stereo(BondStereo::E)
                    .with_stereo_atoms(y, double_ref),
            )
            .unwrap();
        let removed = builder
            .add_bond(
                BondSpec::new(removed_left, removed_right, BondOrder::Double)
                    .with_stereo(BondStereo::Any),
            )
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let atom_count = molecule.num_atoms();
        let bond_count = molecule.num_bonds();
        let mut known_bonds = vec![false; molecule.num_bonds()];
        known_bonds[removed.index()] = true;

        clean_mol_stereo(
            &mut molecule,
            &vec![false; atom_count],
            &vec![false; atom_count],
            &vec![false; bond_count],
            &known_bonds,
        )
        .unwrap();

        assert_eq!(
            molecule.bonds()[slash.index()].direction(),
            BondDirection::EndUpRight
        );
    }

    #[test]
    fn find_potential_stereo_controller_matches_rdkit_option_matrix_and_cache_state() {
        for (clean, flag_possible, expected_centers) in [
            (false, false, Vec::new()),
            (false, true, vec![1, 3]),
            (true, false, Vec::new()),
            (true, true, vec![1, 3]),
        ] {
            let mut molecule = Molecule::from_smiles("FC(Cl)C(Br)I").unwrap();
            let result =
                find_potential_stereo_in_workspace(&mut molecule, clean, flag_possible).unwrap();

            assert_eq!(
                result.iter().map(StereoInfo::center).collect::<Vec<_>>(),
                expected_centers
                    .into_iter()
                    .map(|index| StereoCenter::Atom(AtomId::new(index)))
                    .collect::<Vec<_>>()
            );
            for info in &result {
                assert_eq!(info.stereo_type(), StereoType::AtomTetrahedral);
                assert_eq!(info.specified(), StereoSpecified::Unspecified);
                assert_eq!(info.descriptor(), StereoDescriptor::None);
                assert_eq!(info.permutation(), 0);
            }
            assert!(
                molecule
                    .derived_cache()
                    .rings
                    .as_ref()
                    .is_some_and(crate::RingInfo::is_symm_sssr)
            );
            assert!(molecule.derived_cache().valence.is_some());
            assert_eq!(
                molecule.potential_stereo_cache().as_deref(),
                Some(result.as_slice())
            );
            for atom in molecule.atoms() {
                assert!(atom.is_prop_computed("_chiralAtomRank"));
                assert!(atom.prop("_chiralAtomRank").is_some());
            }
        }
    }

    #[test]
    fn find_potential_stereo_controller_preserves_rdkit_record_order_and_fields() {
        let mut molecule = Molecule::from_smiles("F[C@](Cl)(Br)I").unwrap();
        let result = find_potential_stereo_in_workspace(&mut molecule, false, true).unwrap();

        assert_eq!(result.len(), 1);
        assert_eq!(result[0].stereo_type(), StereoType::AtomTetrahedral);
        assert_eq!(result[0].specified(), StereoSpecified::Specified);
        assert_eq!(result[0].center(), StereoCenter::Atom(AtomId::new(1)));
        assert_eq!(
            result[0].descriptor(),
            StereoDescriptor::TetrahedralCounterclockwise
        );
        assert_eq!(result[0].permutation(), 0);
        assert_eq!(
            result[0].controlling_atoms(),
            [0, 2, 3, 4]
                .map(|index| ControllingAtom::Atom(AtomId::new(index)))
                .as_slice()
        );
    }

    fn assert_stereo_info_matches_golden(
        case_id: &str,
        actual: &[StereoInfo],
        expected: &serde_json::Value,
    ) {
        let expected = expected.as_array().expect("golden stereo_info array");
        assert_eq!(actual.len(), expected.len(), "{case_id}: stereo row count");
        for (row_index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            let context = format!("{case_id}: stereo row {row_index}");
            let expected_type = expected["type"].as_str().unwrap();
            assert_eq!(
                actual.stereo_type(),
                match expected_type {
                    "Atom_Tetrahedral" => StereoType::AtomTetrahedral,
                    "Atom_SquarePlanar" => StereoType::AtomSquarePlanar,
                    "Atom_TrigonalBipyramidal" => StereoType::AtomTrigonalBipyramidal,
                    "Atom_Octahedral" => StereoType::AtomOctahedral,
                    "Bond_Double" => StereoType::BondDouble,
                    "Bond_Cumulene_Even" => StereoType::BondEvenCumulene,
                    "Bond_Atropisomer" => StereoType::BondAtropisomer,
                    other => panic!("{context}: unknown golden stereo type {other}"),
                },
                "{context}: type"
            );
            assert_eq!(
                actual.specified(),
                match expected["specified"].as_str().unwrap() {
                    "Unspecified" => StereoSpecified::Unspecified,
                    "Specified" => StereoSpecified::Specified,
                    "Unknown" => StereoSpecified::Unknown,
                    other => panic!("{context}: unknown golden specified state {other}"),
                },
                "{context}: specified"
            );
            let center = usize::try_from(expected["centered_on"].as_u64().unwrap()).unwrap();
            assert_eq!(
                actual.center(),
                if actual.stereo_type().is_atom_centered() {
                    StereoCenter::Atom(AtomId::new(center))
                } else {
                    StereoCenter::Bond(BondId::new(center))
                },
                "{context}: center"
            );
            assert_eq!(
                actual.descriptor(),
                match expected["descriptor"].as_str().unwrap() {
                    "NoValue" => StereoDescriptor::None,
                    "Tet_CW" => StereoDescriptor::TetrahedralClockwise,
                    "Tet_CCW" => StereoDescriptor::TetrahedralCounterclockwise,
                    "Bond_Cis" => StereoDescriptor::BondCis,
                    "Bond_Trans" => StereoDescriptor::BondTrans,
                    "Bond_AtropCW" => StereoDescriptor::BondAtropClockwise,
                    "Bond_AtropCCW" => StereoDescriptor::BondAtropCounterclockwise,
                    other => panic!("{context}: unknown golden descriptor {other}"),
                },
                "{context}: descriptor"
            );
            assert_eq!(
                actual.permutation(),
                u32::try_from(expected["permutation"].as_u64().unwrap()).unwrap(),
                "{context}: permutation"
            );
            let expected_controllers = expected["controlling_atoms"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| {
                    let index = value.as_u64().unwrap();
                    if index == u64::from(u32::MAX) {
                        ControllingAtom::Missing
                    } else {
                        ControllingAtom::Atom(AtomId::new(usize::try_from(index).unwrap()))
                    }
                })
                .collect::<Vec<_>>();
            assert_eq!(
                actual.controlling_atoms(),
                expected_controllers,
                "{context}: controlling atoms"
            );
        }
    }

    fn assert_potential_stereo_after_state(
        case_id: &str,
        molecule: &Molecule,
        expected: &serde_json::Value,
    ) {
        let expected_atoms = expected["atoms"].as_array().unwrap();
        assert_eq!(
            molecule.num_atoms(),
            expected_atoms.len(),
            "{case_id}: atoms"
        );
        for (atom_index, (atom, expected_atom)) in
            molecule.atoms().iter().zip(expected_atoms).enumerate()
        {
            let context = format!("{case_id}: atom {atom_index}");
            assert_eq!(
                atom.chiral_tag(),
                match expected_atom["chiral_tag"].as_str().unwrap() {
                    "CHI_UNSPECIFIED" => ChiralTag::Unspecified,
                    "CHI_TETRAHEDRAL_CW" => ChiralTag::TetrahedralCw,
                    "CHI_TETRAHEDRAL_CCW" => ChiralTag::TetrahedralCcw,
                    "CHI_TETRAHEDRAL" => ChiralTag::Tetrahedral,
                    "CHI_SQUAREPLANAR" => ChiralTag::SquarePlanar,
                    "CHI_TRIGONALBIPYRAMIDAL" => ChiralTag::TrigonalBipyramidal,
                    "CHI_OCTAHEDRAL" => ChiralTag::Octahedral,
                    other => panic!("{context}: unknown golden chiral tag {other}"),
                },
                "{context}: chiral tag"
            );
            let values = expected_atom["properties"]["values"].as_object().unwrap();
            let expected_rank = values["_chiralAtomRank"].as_u64().unwrap();
            assert_eq!(
                atom.prop("_chiralAtomRank"),
                Some(expected_rank.to_string().as_str()),
                "{context}: chiral rank"
            );
            assert!(atom.is_prop_computed("_chiralAtomRank"), "{context}");

            for property in [
                "_ringStereoOtherAtom",
                "_ringStereochemCand",
                "_ringStereoAtoms",
            ] {
                let expected_value = values.get(property);
                assert_eq!(
                    atom.prop(property).is_some(),
                    expected_value.is_some(),
                    "{context}: {property} presence"
                );
                let Some(expected_value) = expected_value else {
                    continue;
                };
                match property {
                    "_ringStereoOtherAtom" => assert_eq!(
                        atom.prop(property).unwrap().parse::<u64>().unwrap(),
                        expected_value.as_u64().unwrap(),
                        "{context}: {property}"
                    ),
                    "_ringStereochemCand" => assert_eq!(
                        atom.prop(property).unwrap() == "1",
                        expected_value.as_bool().unwrap(),
                        "{context}: {property}"
                    ),
                    "_ringStereoAtoms" => {
                        let actual = atom
                            .prop(property)
                            .unwrap()
                            .split(',')
                            .filter(|value| !value.is_empty())
                            .map(|value| value.parse::<i64>().unwrap())
                            .collect::<Vec<_>>();
                        let expected = expected_value
                            .as_array()
                            .unwrap()
                            .iter()
                            .map(|value| value.as_i64().unwrap())
                            .collect::<Vec<_>>();
                        assert_eq!(actual, expected, "{context}: {property}");
                    }
                    _ => unreachable!(),
                }
            }
        }

        let expected_bonds = expected["bonds"].as_array().unwrap();
        assert_eq!(
            molecule.num_bonds(),
            expected_bonds.len(),
            "{case_id}: bonds"
        );
        for (bond_index, (bond, expected_bond)) in
            molecule.bonds().iter().zip(expected_bonds).enumerate()
        {
            let context = format!("{case_id}: bond {bond_index}");
            assert_eq!(
                format!("{:?}", bond.direction()).to_ascii_uppercase(),
                expected_bond["direction"]
                    .as_str()
                    .unwrap()
                    .replace('_', "")
                    .to_ascii_uppercase(),
                "{context}: direction"
            );
            let expected_stereo = expected_bond["stereo"].as_str().unwrap();
            assert_eq!(
                bond.stereo(),
                match expected_stereo {
                    "STEREONONE" => crate::BondStereo::None,
                    "STEREOANY" => crate::BondStereo::Any,
                    "STEREOCIS" => crate::BondStereo::Cis,
                    "STEREOTRANS" => crate::BondStereo::Trans,
                    "STEREOE" => crate::BondStereo::E,
                    "STEREOZ" => crate::BondStereo::Z,
                    "STEREOATROPCW" => crate::BondStereo::AtropCw,
                    "STEREOATROPCCW" => crate::BondStereo::AtropCcw,
                    other => panic!("{context}: unknown golden bond stereo {other}"),
                },
                "{context}: stereo"
            );
            let expected_stereo_atoms = expected_bond["stereo_atoms"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| AtomId::new(usize::try_from(value.as_u64().unwrap()).unwrap()))
                .collect::<Vec<_>>();
            assert_eq!(
                bond.stereo_atoms()
                    .map(|atoms| atoms.to_vec())
                    .unwrap_or_default(),
                expected_stereo_atoms,
                "{context}: stereo atoms"
            );
        }
    }

    #[test]
    fn find_potential_stereo_matches_every_focused_golden_row_and_option_exactly() {
        let golden_path = crate::test_data::expected_path_for_profile(
            "stereo",
            "rdkit",
            "python_stereoisomer_focused",
            "python_stereoisomer.jsonl",
        );
        let golden = std::fs::read_to_string(&golden_path)
            .unwrap_or_else(|error| panic!("failed to read {}: {error}", golden_path.display()));
        let mut compared_runs = 0_usize;
        for line in golden.lines().filter(|line| !line.is_empty()) {
            let record: serde_json::Value = serde_json::from_str(line).unwrap();
            let runs = record["potential_stereo_runs"].as_array().unwrap();
            if runs.is_empty() {
                continue;
            }
            let case_id = record["case_id"].as_str().unwrap();
            assert_eq!(record["source"]["kind"].as_str(), Some("smiles"));
            let mut source =
                Molecule::from_smiles(record["source"]["value"].as_str().unwrap()).unwrap();
            for mutation in record["source"]["mutations"].as_array().unwrap() {
                match mutation["kind"].as_str().unwrap() {
                    "set_atom_chiral_tag" => {
                        let atom = usize::try_from(mutation["atom"].as_u64().unwrap()).unwrap();
                        let tag = match mutation["value"].as_str().unwrap() {
                            "TetrahedralCw" => ChiralTag::TetrahedralCw,
                            "TetrahedralCcw" => ChiralTag::TetrahedralCcw,
                            "Unspecified" => ChiralTag::Unspecified,
                            other => panic!("{case_id}: unknown mutation chiral tag {other}"),
                        };
                        source.topology_block_mut().atoms[atom].set_chiral_tag(tag);
                    }
                    other => panic!("{case_id}: unsupported focused mutation {other}"),
                }
            }
            assert_eq!(runs.len(), 4, "{case_id}: complete option matrix");
            let mut seen_options = BTreeSet::new();
            for run in runs {
                assert_eq!(run["result"]["status"].as_str(), Some("ok"));
                let clean = run["options"]["clean_it"].as_bool().unwrap();
                let flag_possible = run["options"]["flag_possible"].as_bool().unwrap();
                assert!(seen_options.insert((clean, flag_possible)), "{case_id}");
                let mut working = source.clone();
                let actual =
                    find_potential_stereo_in_workspace(&mut working, clean, flag_possible).unwrap();
                let context = format!("{case_id} clean={clean} possible={flag_possible}");
                assert_stereo_info_matches_golden(&context, &actual, &run["result"]["stereo_info"]);
                assert_potential_stereo_after_state(&context, &working, &run["after"]);
                assert_eq!(
                    working.potential_stereo_cache().as_deref(),
                    Some(actual.as_slice()),
                    "{context}: typed potential-stereo cache"
                );
                compared_runs += 1;
            }
            assert_eq!(
                seen_options,
                BTreeSet::from([(false, false), (false, true), (true, false), (true, true)]),
                "{case_id}: option coverage"
            );
        }
        assert_eq!(compared_runs, 20);
    }
}
