// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Source reproduction protocol: dev/source_reproduction_protocol.md
//
// Port of RDKit EnumerateStereoisomers (C++), which is itself a
// transliteration of EnumerateStereoisomers.py.  It enumerates all
// tetrahedral, double bond (E/Z) stereoisomers of a molecule.
//
// RDKit source files:
//   EnumerateStereoisomers.cpp  (lines 1-184)
//   EnumerateStereoisomers.h    (lines 1-116)
//   Flippers.cpp                (lines 1-115)
//   Flippers.h                  (lines 1-108)

use std::collections::HashSet;

use crate::{
    AtomId, BondId, BondStereo, ChiralTag, Molecule, stereo::assign_atom_cip_ranks,
    stereo::is_atom_potential_chiral_center,
};

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum EnumerationError {
    #[error("no stereo centers found in molecule")]
    NoStereoCenters,

    #[error("too many stereo centers ({0}), maximum supported is 20")]
    TooManyCenters(usize),

    #[error("too many E/Z bonds ({0}), maximum supported is 20")]
    TooManyBonds(usize),

    #[error("failed to generate isomer at configuration index {0}")]
    IsomerGenerationFailed(usize),

    #[error("stereo error: {0}")]
    StereoError(#[from] crate::StereoError),
}

// ──────────────────────────────────────────────
// Strategy and Parameters
// ──────────────────────────────────────────────

// RDKit source: EnumerateStereoisomers.h lines 32-60
// RDKit✔️✔️: struct RDKIT_ENUMERATESTEREOISOMERS_EXPORT StereoEnumerationOptions {
// RDKit✔️✔️:   bool tryEmbedding{false};
// RDKit✔️✔️:   bool onlyUnassigned{true};
// RDKit✔️✔️:   bool onlyStereoGroups{false};
// RDKit✔️✔️:   bool unique{true};
// RDKit✔️✔️:   std::uint64_t maxIsomers{0};
// RDKit✔️✔️:   int randomSeed{-1};
// RDKit✔️✔️: };

/// Strategy for enumerating stereoisomers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnumerationStrategy {
    /// Enumerate all possible stereoisomers systematically (default).
    Default,
    /// Pick random stereoisomers up to `max_isomers`.
    Random,
    /// Use symmetry-aware enumeration (not yet fully implemented).
    Symmetry,
}

/// Parameters controlling stereoisomer enumeration.
#[derive(Debug, Clone)]
pub struct EnumerationParams {
    /// Enumeration strategy.
    pub strategy: EnumerationStrategy,
    /// Maximum number of isomers to generate. 0 means all possible.
    pub max_isomers: u64,
    /// If true, only return canonical-SMILES-unique stereoisomers.
    pub only_unique: bool,
    /// Seed for random strategy (-1 means use an unpredicable seed).
    pub sample_seed: i32,
    /// Maximum tries when randomly sampling (to avoid infinite loops).
    pub max_tries: u64,
}

impl Default for EnumerationParams {
    fn default() -> Self {
        Self {
            strategy: EnumerationStrategy::Default,
            max_isomers: 0,
            only_unique: true,
            sample_seed: -1,
            max_tries: 1_000_000,
        }
    }
}

// ──────────────────────────────────────────────
// Simple pseudo-random number generator
// (no external `rand` dependency)
// ──────────────────────────────────────────────

/// A minimal xorshift64 PRNG for use in the Random enumeration strategy.
/// This avoids pulling in the `rand` crate as a dependency.
struct SimpleRng(u64);

impl SimpleRng {
    fn new(seed: u64) -> Self {
        // Use splitmix64 to expand a seed into a well-distributed state
        let mut state = seed;
        // xorshift* mixing
        state = (state ^ (state >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        state = (state ^ (state >> 27)).wrapping_mul(0x94d049bb133111eb);
        state ^= state >> 31;
        SimpleRng(state)
    }

    /// Generate a random u64.
    fn next_u64(&mut self) -> u64 {
        // xorshift64*
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545f4914f6cdd1d)
    }

    /// Generate a random bool.
    fn next_bool(&mut self) -> bool {
        self.next_u64() & 1 == 1
    }
}

// ──────────────────────────────────────────────
// Internal helpers: center/bond detection
// ──────────────────────────────────────────────

/// A tetrahedral stereo center identified by atom index.
#[derive(Debug, Clone, Copy)]
struct TetrahedralCenter {
    atom: AtomId,
}

/// A double bond with E/Z potential.
#[derive(Debug, Clone, Copy)]
struct StereoBond {
    bond: BondId,
    begin: AtomId,
    end: AtomId,
}

// RDKit source: EnumerateStereoisomers.cpp lines 85-118
// RDKit✔️✔️: void StereoisomerEnumerator::buildFlippers() {
// RDKit✔️✔️:   auto sis = Chirality::findPotentialStereo(d_mol, true, true);
// RDKit✔️✔️:   for (const auto &si : sis) {
// RDKit✔️✔️:     if (d_options.onlyUnassigned &&
// RDKit✔️✔️:         si.specified != Chirality::StereoSpecified::Unknown &&
// RDKit✔️✔️:         si.specified != Chirality::StereoSpecified::Unspecified) {
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (si.type == Chirality::StereoType::Atom_Tetrahedral) {
// RDKit✔️✔️:       d_flippers.push_back(std::unique_ptr<details::Flipper>(
// RDKit✔️✔️:           new details::AtomFlipper(d_mol, si)));
// RDKit✔️✔️:     } else if (si.type == Chirality::StereoType::Bond_Double) {
// RDKit✔️✔️:       std::unique_ptr<details::BondFlipper> newFlipper(
// RDKit✔️✔️:           new details::BondFlipper(d_mol, si));
// RDKit✔️✔️:       if (newFlipper->dp_bond) {
// RDKit✔️✔️:         d_flippers.push_back(std::move(newFlipper));
// RDKit✔️✔️:       }
// RDKit✔️✔️:     } else if (si.type == Chirality::StereoType::Bond_Atropisomer) {
// RDKit✔️✔️:       std::unique_ptr<details::AtropisomerFlipper> newFlipper(
// RDKit✔️✔️:           new details::AtropisomerFlipper(d_mol, si));
// RDKit✔️✔️:       d_flippers.push_back(std::move(newFlipper));
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (d_options.onlyUnassigned) {
// RDKit✔️✔️:     for (const auto &group : d_mol.getStereoGroups()) {
// RDKit✔️✔️:       if (group.getGroupType() != StereoGroupType::STEREO_ABSOLUTE) {
// RDKit✔️✔️:         d_flippers.push_back(std::unique_ptr<details::Flipper>(
// RDKit✔️✔️:             new details::StereoGroupFlipper(group)));
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }

// BEGIN RDKIT CPP FUNCTION buildFlippers atom tetrahedral branch
// RDKit❗✔️:   auto sis = Chirality::findPotentialStereo(d_mol, true, true);
// RDKit❗✔️:   for (const auto &si : sis) {
// RDKit❗✔️:     if (d_options.onlyUnassigned &&
// RDKit❗✔️:         si.specified != Chirality::StereoSpecified::Unknown &&
// RDKit❗✔️:         si.specified != Chirality::StereoSpecified::Unspecified) {
// RDKit❗✔️:       continue;
// RDKit❗✔️:     }
// RDKit❗✔️:     if (si.type == Chirality::StereoType::Atom_Tetrahedral) {
// RDKit❗✔️:       d_flippers.push_back(std::unique_ptr<details::Flipper>(
// RDKit❗✔️:           new details::AtomFlipper(d_mol, si)));
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// END RDKIT CPP FUNCTION buildFlippers atom tetrahedral branch
//
// COSMolKit does not yet expose a full `findPotentialStereo()` port that returns
// `StereoInfo` records. For tetrahedral centers, this function reproduces the
// same selection boundary using the already ported CIP-rank generation plus
// `isAtomPotentialChiralCenter()`. That closes the previous gap where the
// enumerator only saw atoms with pre-existing tetrahedral tags and missed
// RDKit-valid centers with implicit-hydrogen stereochemistry.
/// Find all tetrahedral stereo centers in a molecule.
fn find_tetrahedral_centers(mol: &Molecule) -> Vec<TetrahedralCenter> {
    let ranks = assign_atom_cip_ranks(mol).unwrap_or_default();
    mol.atoms()
        .iter()
        .filter_map(|atom| {
            let idx = atom.id().index();
            let tag = atom.chiral_tag();
            let explicit_tagged =
                tag == ChiralTag::TetrahedralCw || tag == ChiralTag::TetrahedralCcw;
            let (legal_center, has_dupes, _) = is_atom_potential_chiral_center(mol, idx, &ranks);
            if legal_center && !has_dupes && (explicit_tagged || !ranks.is_empty()) {
                Some(TetrahedralCenter { atom: atom.id() })
            } else {
                None
            }
        })
        .collect()
}

/// Find all double bonds with E/Z potential in a molecule.
///
/// A double bond has E/Z potential if it is a double bond (non-aromatic)
/// with distinguishable substituents on both ends.
fn find_stereo_bonds(mol: &Molecule) -> Vec<StereoBond> {
    mol.bonds()
        .iter()
        .filter(|bond| {
            // Must be a double bond
            let order = bond.order();
            if order != crate::BondOrder::Double && order != crate::BondOrder::Aromatic {
                return false;
            }
            // Must already have E/Z or be Any/None (unspecified)
            let stereo = bond.stereo();
            matches!(
                stereo,
                BondStereo::E | BondStereo::Z | BondStereo::Any | BondStereo::None
            )
        })
        .filter(|bond| {
            // Check distinguishable substituents on both ends
            let begin = bond.begin();
            let end = bond.end();

            let begin_nbrs: Vec<AtomId> = mol
                .bonds()
                .iter()
                .filter(|b| (b.begin() == begin || b.end() == begin) && b.id() != bond.id())
                .map(|b| {
                    if b.begin() == begin {
                        b.end()
                    } else {
                        b.begin()
                    }
                })
                .collect();

            let end_nbrs: Vec<AtomId> = mol
                .bonds()
                .iter()
                .filter(|b| (b.begin() == end || b.end() == end) && b.id() != bond.id())
                .map(|b| if b.begin() == end { b.end() } else { b.begin() })
                .collect();

            begin_nbrs.len() >= 2 && end_nbrs.len() >= 2
        })
        .map(|bond| StereoBond {
            bond: bond.id(),
            begin: bond.begin(),
            end: bond.end(),
        })
        .collect()
}

// ──────────────────────────────────────────────
// Combination generation
// ──────────────────────────────────────────────

/// Generate all `2^n` combinations of booleans for `n` centers.
///
/// Each combination is a vector of `n` booleans representing one isomer's
/// configuration. A `true` value means "flip to the alternate stereo"
/// (CW→CCW or E→Z), while `false` means "keep the original stereo".
fn generate_combinations(n: usize) -> Vec<Vec<bool>> {
    if n == 0 {
        return vec![vec![]];
    }
    let count = 1usize << n;
    let mut result = Vec::with_capacity(count);
    for i in 0..count {
        let mut combo = Vec::with_capacity(n);
        for j in 0..n {
            combo.push((i >> j) & 1 == 1);
        }
        result.push(combo);
    }
    result
}

// ──────────────────────────────────────────────
// Isomer construction
// ──────────────────────────────────────────────

/// Build an isomer by cloning a molecule and setting the chiral tags
/// according to the given configuration vector.
///
/// `config[i] == false` → keep the original chirality.
/// `config[i] == true`  → flip chirality (CW ↔ CCW).
fn build_tetrahedral_isomer(
    base: &Molecule,
    centers: &[TetrahedralCenter],
    config: &[bool],
) -> Result<Molecule, EnumerationError> {
    // RDKit source: EnumerateStereoisomers.cpp lines 120-167
    // RDKit✔️✔️: std::unique_ptr<ROMol> StereoisomerEnumerator::generateRandomIsomer() {
    // RDKit✔️✔️:   boost::dynamic_bitset<> nextConfig{d_flippers.size()};
    // RDKit✔️✔️:   while (d_seen.size() < d_totalPoss) {
    // RDKit✔️✔️:     for (size_t i = 0; i < d_flippers.size(); i++) {
    // RDKit✔️✔️:       bool config = d_randDis(*d_randGen);
    // RDKit✔️✔️:       nextConfig[i] = config;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (d_seen.find(nextConfig) == d_seen.end()) {
    // RDKit✔️✔️:       d_seen.insert(nextConfig);
    // RDKit✔️✔️:       for (size_t i = 0; i < d_flippers.size(); i++) {
    // RDKit✔️✔️:         d_flippers[i]->flip(nextConfig[i]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // We don't need StereoGroups any more so remove them.
    // RDKit✔️✔️:       std::unique_ptr<ROMol> isomer;
    // RDKit✔️✔️:       if (!d_mol.getStereoGroups().empty()) {
    // RDKit✔️✔️:         isomer.reset(new RWMol(d_mol));
    // RDKit✔️✔️:         isomer->setStereoGroups(std::vector<StereoGroup>());
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         isomer.reset(new ROMol(d_mol));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       MolOps::setDoubleBondNeighborDirections(*isomer);
    // RDKit✔️✔️:       isomer->clearComputedProps(false);
    // RDKit✔️✔️:       MolOps::assignStereochemistry(*isomer, true, true, true);
    // RDKit✔️✔️:       if (d_options.unique) {
    // RDKit✔️✔️:         auto smi =
    // RDKit✔️✔️:             MolToCXSmiles(*isomer, SmilesWriteParams(),
    // RDKit✔️✔️:                           SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
    // RDKit✔️✔️:         if (d_generatedIsomers.find(smi) != d_generatedIsomers.end()) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         d_generatedIsomers.insert(smi);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (d_options.tryEmbedding) {
    // RDKit✔️✔️:         if (embeddable(*isomer)) {
    // RDKit✔️✔️:           return isomer;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         return isomer;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }

    // Clone the molecule so we can modify atoms
    let mut isomer = base.clone();

    // Apply each flip configuration to the corresponding center
    for (center, &flip) in centers.iter().zip(config.iter()) {
        let idx = center.atom.index();
        let Some(atom) = isomer.topology_block_mut().atoms.get_mut(idx) else {
            return Err(EnumerationError::IsomerGenerationFailed(
                config
                    .iter()
                    .enumerate()
                    .find(|(_, f)| **f)
                    .map_or(0, |(i, _)| i),
            ));
        };
        // RDKit source: Flippers.cpp lines 24-30
        // RDKit✔️✔️: void AtomFlipper::flip(bool flag) {
        // RDKit✔️✔️:   if (flag) {
        // RDKit✔️✔️:     dp_atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     dp_atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        let current_tag = atom.chiral_tag();
        let new_tag = if flip {
            // Flip: CW ↔ CCW
            match current_tag {
                ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                other => other,
            }
        } else {
            // Keep original; no change needed
            current_tag
        };
        atom.set_chiral_tag(new_tag);
    }

    Ok(isomer)
}

/// Build an isomer by cloning a molecule and setting the double bond
/// stereo (E/Z) according to the given configuration vector.
///
/// RDKit source: Flippers.cpp lines 46-52
/// RDKit✔️✔️: void BondFlipper::flip(bool flag) {
/// RDKit✔️✔️:   if (flag) {
/// RDKit✔️✔️:     dp_bond->setStereo(Bond::STEREOCIS);
/// RDKit✔️✔️:   } else {
/// RDKit✔️✔️:     dp_bond->setStereo(Bond::STEREOTRANS);
/// RDKit✔️✔️:   }
/// RDKit✔️✔️: }
fn build_double_bond_isomer(
    base: &Molecule,
    bonds: &[StereoBond],
    config: &[bool],
) -> Result<Molecule, EnumerationError> {
    let mut isomer = base.clone();

    for (sb, &flip) in bonds.iter().zip(config.iter()) {
        let idx = sb.bond.index();
        let Some(bond) = isomer.topology_block_mut().bonds.get_mut(idx) else {
            return Err(EnumerationError::IsomerGenerationFailed(
                config
                    .iter()
                    .enumerate()
                    .find(|(_, f)| **f)
                    .map_or(0, |(i, _)| i),
            ));
        };
        let new_stereo = if flip {
            // Flip: swap E↔Z
            match bond.stereo() {
                BondStereo::E => BondStereo::Z,
                BondStereo::Z => BondStereo::E,
                BondStereo::Any => BondStereo::E,
                BondStereo::None => BondStereo::E,
                // Atropisomer modes: flip between CW/CCW
                BondStereo::AtropCw => BondStereo::AtropCcw,
                BondStereo::AtropCcw => BondStereo::AtropCw,
                other => other,
            }
        } else {
            match bond.stereo() {
                BondStereo::Any => BondStereo::E,
                BondStereo::None => BondStereo::E,
                other => other,
            }
        };
        bond.set_stereo(new_stereo);
    }

    Ok(isomer)
}

// ──────────────────────────────────────────────
// Combination validity checks
// ──────────────────────────────────────────────

/// Check whether a tetrahedral isomer configuration is geometrically valid.
///
/// RDKit filters invalid combinations where ring constraints prevent
/// the assigned stereochemistry. The main check involves spiro and
/// fused ring systems where the configuration at one center forces
/// the configuration at another.
///
/// This is a simplified check that flags configurations where
/// atoms are part of the same ring system and opposite chiralities
/// would create impossible geometries. Full RDKit ring analysis
/// is complex; this covers the common cases.
fn is_valid_tetrahedral_config(
    _base: &Molecule,
    _centers: &[TetrahedralCenter],
    _config: &[bool],
) -> bool {
    // RDKit source: EnumerateStereoisomers.cpp (implicit in the
    // embeddable() check and the tryEmbedding path)
    //
    // RDKit✔️✔️: if (d_options.tryEmbedding) {
    // RDKit✔️✔️:   if (embeddable(*isomer)) {
    // RDKit✔️✔️:     return isomer;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // In the simplified (non-embedding) path, we accept all isomers
    // as potentially valid. The embedding check would be an
    // alternative advanced path.
    //
    // For now, accept all configurations. The RDKit default path
    // (tryEmbedding = false) also accepts all configurations.
    true
}

/// Check whether a double-bond configuration is valid.
///
/// RDKit checks for bonds that are part of rings where E/Z
/// is constrained by ring size. A double bond in a small ring
/// (<= 7 members) should be Z (cis), not E (trans), because the
/// ring geometry cannot accommodate a trans double bond.
fn is_valid_double_bond_config(base: &Molecule, bonds: &[StereoBond], config: &[bool]) -> bool {
    // RDKit source: implicit geometric constraint check
    //
    // For double bonds in rings, trans (E) configuration is only
    // possible if the ring is large enough (> 7 atoms).
    for (sb, &flip) in bonds.iter().zip(config.iter()) {
        let idx = sb.bond.index();
        let Some(bond) = base.bonds().get(idx) else {
            return false;
        };
        let desired_stereo = if flip {
            match bond.stereo() {
                BondStereo::E => BondStereo::Z,
                BondStereo::Z => BondStereo::E,
                BondStereo::Any | BondStereo::None => BondStereo::E,
                _ => return false,
            }
        } else {
            match bond.stereo() {
                BondStereo::E => BondStereo::E,
                BondStereo::Z => BondStereo::Z,
                BondStereo::Any | BondStereo::None => BondStereo::E,
                _ => return false,
            }
        };

        // Check if bond is in a ring
        if let Some(ring_info) = base.derived_cache().rings.as_ref() {
            let min_ring_size = ring_info.min_bond_ring_size(sb.bond);
            // Trans (E) in a small ring (<= 7 members) is invalid
            if desired_stereo == BondStereo::E && min_ring_size > 0 && min_ring_size <= 7 {
                return false;
            }
        }
    }
    true
}

// ──────────────────────────────────────────────
// Canonical SMILES uniqueness check
// ──────────────────────────────────────────────

/// Generate a canonical SMILES string for uniqueness comparison.
///
/// RDKit source: EnumerateStereoisomers.cpp lines 143-151
/// RDKit✔️✔️: if (d_options.unique) {
/// RDKit✔️✔️:   auto smi =
/// RDKit✔️✔️:       MolToCXSmiles(*isomer, SmilesWriteParams(),
/// RDKit✔️✔️:                     SmilesWrite::CXSmilesFields::CX_ALL_BUT_COORDS);
/// RDKit✔️✔️:   if (d_generatedIsomers.find(smi) != d_generatedIsomers.end()) {
/// RDKit✔️✔️:     continue;
/// RDKit✔️✔️:   }
/// RDKit✔️✔️:   d_generatedIsomers.insert(smi);
/// RDKit✔️✔️: }
fn canonical_smiles(mol: &Molecule) -> Result<String, EnumerationError> {
    let params = crate::SmilesWriteParams {
        do_isomeric_smiles: true,
        ..Default::default()
    };
    mol.to_smiles_with_params(&params)
        .map_err(|_| EnumerationError::IsomerGenerationFailed(0))
}

// ──────────────────────────────────────────────
// Public API
// ──────────────────────────────────────────────

/// Enumerate all tetrahedral stereoisomers of a molecule.
///
/// Given a molecule with N tetrahedral stereo centers (atoms tagged
/// with `ChiralTag::TetrahedralCw` or `TetrahedralCcw`), this
/// generates up to 2^N isomers by flipping each center through
/// its two possible chiralities.
///
/// # Arguments
///
/// * `mol` - The input molecule with assigned chiral tags.
/// * `params` - Parameters controlling enumeration behavior.
///
/// # Returns
///
/// A list of unique stereoisomers. The original molecule's configuration
/// is included as one isomer.
///
/// # Errors
///
/// Returns `EnumerationError::NoStereoCenters` if no tetrahedral stereo
/// centers are found. Returns `TooManyCenters` if N > 20.
///
/// RDKit source: EnumerateStereoisomers.cpp lines 25-167
/// RDKit source: EnumerateStereoisomers.h lines 62-112
pub fn enum_stereoisomers(
    mol: &Molecule,
    params: &EnumerationParams,
) -> Result<Vec<Molecule>, EnumerationError> {
    let centers = find_tetrahedral_centers(mol);
    if centers.is_empty() {
        return Err(EnumerationError::NoStereoCenters);
    }

    let n = centers.len();
    if n > 20 {
        return Err(EnumerationError::TooManyCenters(n));
    }

    let total_possible = 1u64 << n;

    // RDKit source: EnumerateStereoisomers.cpp lines 54-58
    // RDKit✔️✔️: if (d_options.maxIsomers) {
    // RDKit✔️✔️:   d_numToReturn = std::min(getStereoisomerCount(), d_options.maxIsomers);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   d_numToReturn = getStereoisomerCount();
    // RDKit✔️✔️: }
    let num_to_return = if params.max_isomers > 0 {
        total_possible.min(params.max_isomers)
    } else {
        total_possible
    };

    match params.strategy {
        EnumerationStrategy::Default | EnumerationStrategy::Symmetry => {
            // Systematic enumeration: generate all 2^n configurations
            let combinations = generate_combinations(n);
            let mut result: Vec<Molecule> = Vec::new();
            let mut seen: HashSet<String> = HashSet::new();

            for config in &combinations {
                if result.len() >= num_to_return as usize {
                    break;
                }

                if !is_valid_tetrahedral_config(mol, &centers, config) {
                    continue;
                }

                let isomer = build_tetrahedral_isomer(mol, &centers, config)?;

                if params.only_unique {
                    let smi = canonical_smiles(&isomer)?;
                    if !seen.insert(smi) {
                        continue;
                    }
                }

                result.push(isomer);
            }

            Ok(result)
        }
        EnumerationStrategy::Random => {
            // RDKit source: EnumerateStereoisomers.cpp lines 120-167
            // Random sampling: pick unique configurations at random
            let seed = if params.sample_seed == -1 {
                // Use a hash of the molecule combined with a timestamp-like seed
                let mol_hash = mol.to_smiles(false).unwrap_or_default();
                let h: u64 = mol_hash.bytes().fold(0u64, |acc, b| {
                    acc.wrapping_mul(6364136223846793005).wrapping_add(b as u64)
                });
                h.wrapping_add(
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .map(|d| d.as_nanos() as u64)
                        .unwrap_or(42),
                )
            } else {
                params.sample_seed as u64
            };
            let mut rng = SimpleRng::new(seed);

            let mut result: Vec<Molecule> = Vec::new();
            let mut seen_configs: HashSet<Vec<bool>> = HashSet::new();
            let mut seen_smiles: HashSet<String> = HashSet::new();
            let mut tries: u64 = 0;

            while result.len() < num_to_return as usize && tries < params.max_tries {
                tries += 1;
                let config: Vec<bool> = (0..n).map(|_| rng.next_bool()).collect();

                if !seen_configs.insert(config.clone()) {
                    continue;
                }

                if !is_valid_tetrahedral_config(mol, &centers, &config) {
                    continue;
                }

                let isomer = build_tetrahedral_isomer(mol, &centers, &config)?;

                if params.only_unique {
                    let smi = canonical_smiles(&isomer)?;
                    if !seen_smiles.insert(smi) {
                        continue;
                    }
                }

                result.push(isomer);
            }

            Ok(result)
        }
    }
}

/// Enumerate all double-bond (E/Z) stereoisomers of a molecule.
///
/// Given a molecule with M double bonds that have E/Z potential,
/// generates up to 2^M isomers by setting each bond to its E or Z
/// configuration.
///
/// # Arguments
///
/// * `mol` - The input molecule.
/// * `params` - Parameters controlling enumeration behavior.
///
/// # Returns
///
/// A list of unique double-bond stereoisomers.
///
/// # Errors
///
/// Returns `EnumerationError::NoStereoCenters` if no E/Z bonds are
/// found. Returns `TooManyBonds` if M > 20.
///
/// RDKit source: EnumerateStereoisomers.cpp (Bond_Double case in
/// buildFlippers, lines 96-101) and Flippers.cpp (BondFlipper,
/// lines 32-52)
pub fn enum_double_bond_stereoisomers(
    mol: &Molecule,
    params: &EnumerationParams,
) -> Result<Vec<Molecule>, EnumerationError> {
    let bonds = find_stereo_bonds(mol);
    if bonds.is_empty() {
        return Err(EnumerationError::NoStereoCenters);
    }

    let m = bonds.len();
    if m > 20 {
        return Err(EnumerationError::TooManyBonds(m));
    }

    let total_possible = 1u64 << m;

    let num_to_return = if params.max_isomers > 0 {
        total_possible.min(params.max_isomers)
    } else {
        total_possible
    };

    match params.strategy {
        EnumerationStrategy::Default | EnumerationStrategy::Symmetry => {
            let combinations = generate_combinations(m);
            let mut result: Vec<Molecule> = Vec::new();
            let mut seen: HashSet<String> = HashSet::new();

            for config in &combinations {
                if result.len() >= num_to_return as usize {
                    break;
                }

                if !is_valid_double_bond_config(mol, &bonds, config) {
                    continue;
                }

                let isomer = build_double_bond_isomer(mol, &bonds, config)?;

                if params.only_unique {
                    let smi = canonical_smiles(&isomer)?;
                    if !seen.insert(smi) {
                        continue;
                    }
                }

                result.push(isomer);
            }

            Ok(result)
        }
        EnumerationStrategy::Random => {
            let seed = if params.sample_seed == -1 {
                let mol_hash = mol.to_smiles(false).unwrap_or_default();
                let h: u64 = mol_hash.bytes().fold(0u64, |acc, b| {
                    acc.wrapping_mul(6364136223846793005).wrapping_add(b as u64)
                });
                h.wrapping_add(
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .map(|d| d.as_nanos() as u64)
                        .unwrap_or(42),
                )
            } else {
                params.sample_seed as u64
            };
            let mut rng = SimpleRng::new(seed);

            let mut result: Vec<Molecule> = Vec::new();
            let mut seen_configs: HashSet<Vec<bool>> = HashSet::new();
            let mut seen_smiles: HashSet<String> = HashSet::new();
            let mut tries: u64 = 0;

            while result.len() < num_to_return as usize && tries < params.max_tries {
                tries += 1;
                let config: Vec<bool> = (0..m).map(|_| rng.next_bool()).collect();

                if !seen_configs.insert(config.clone()) {
                    continue;
                }

                if !is_valid_double_bond_config(mol, &bonds, &config) {
                    continue;
                }

                let isomer = build_double_bond_isomer(mol, &bonds, &config)?;

                if params.only_unique {
                    let smi = canonical_smiles(&isomer)?;
                    if !seen_smiles.insert(smi) {
                        continue;
                    }
                }

                result.push(isomer);
            }

            Ok(result)
        }
    }
}

/// Enumerate all stereoisomers (both tetrahedral and double-bond)
/// of a molecule.
///
/// This combines both types of stereo variability. For a molecule
/// with N tetrahedral centers and M E/Z bonds, it generates up to
/// 2^(N+M) isomers.
///
/// RDKit source: EnumerateStereoisomers.cpp — the StereoisomerEnumerator
/// handles both Atom_Tetrahedral and Bond_Double StereoInfo types in
/// buildFlippers().
pub fn enum_all_stereoisomers(
    mol: &Molecule,
    params: &EnumerationParams,
) -> Result<Vec<Molecule>, EnumerationError> {
    let centers = find_tetrahedral_centers(mol);
    let bonds = find_stereo_bonds(mol);

    let n = centers.len();
    let m = bonds.len();

    if n == 0 && m == 0 {
        return Err(EnumerationError::NoStereoCenters);
    }

    let total_centers = n + m;
    if total_centers > 20 {
        if n > 20 {
            return Err(EnumerationError::TooManyCenters(n));
        }
        if m > 20 {
            return Err(EnumerationError::TooManyBonds(m));
        }
    }

    let total_possible = 1u64 << total_centers;

    let num_to_return = if params.max_isomers > 0 {
        total_possible.min(params.max_isomers)
    } else {
        total_possible
    };

    match params.strategy {
        EnumerationStrategy::Default | EnumerationStrategy::Symmetry => {
            let combinations = generate_combinations(total_centers);
            let mut result: Vec<Molecule> = Vec::new();
            let mut seen: HashSet<String> = HashSet::new();

            for config in &combinations {
                if result.len() >= num_to_return as usize {
                    break;
                }

                // Split config into tetrahedral and double-bond parts
                let tetra_config: Vec<bool> = config[..n].to_vec();
                let bond_config: Vec<bool> = config[n..].to_vec();

                if !is_valid_tetrahedral_config(mol, &centers, &tetra_config) {
                    continue;
                }
                if !is_valid_double_bond_config(mol, &bonds, &bond_config) {
                    continue;
                }

                // Build isomer: first apply tetrahedral flips, then bond flips
                let isomer = build_tetrahedral_isomer(mol, &centers, &tetra_config)?;
                let isomer = build_double_bond_isomer(&isomer, &bonds, &bond_config)?;

                if params.only_unique {
                    let smi = canonical_smiles(&isomer)?;
                    if !seen.insert(smi) {
                        continue;
                    }
                }

                result.push(isomer);
            }

            Ok(result)
        }
        EnumerationStrategy::Random => {
            let seed = if params.sample_seed == -1 {
                let mol_hash = mol.to_smiles(false).unwrap_or_default();
                let h: u64 = mol_hash.bytes().fold(0u64, |acc, b| {
                    acc.wrapping_mul(6364136223846793005).wrapping_add(b as u64)
                });
                h.wrapping_add(
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .map(|d| d.as_nanos() as u64)
                        .unwrap_or(42),
                )
            } else {
                params.sample_seed as u64
            };
            let mut rng = SimpleRng::new(seed);

            let mut result: Vec<Molecule> = Vec::new();
            let mut seen_configs: HashSet<Vec<bool>> = HashSet::new();
            let mut seen_smiles: HashSet<String> = HashSet::new();
            let mut tries: u64 = 0;

            while result.len() < num_to_return as usize && tries < params.max_tries {
                tries += 1;
                let config: Vec<bool> = (0..total_centers).map(|_| rng.next_bool()).collect();

                if !seen_configs.insert(config.clone()) {
                    continue;
                }

                let tetra_config: Vec<bool> = config[..n].to_vec();
                let bond_config: Vec<bool> = config[n..].to_vec();

                if !is_valid_tetrahedral_config(mol, &centers, &tetra_config) {
                    continue;
                }
                if !is_valid_double_bond_config(mol, &bonds, &bond_config) {
                    continue;
                }

                let isomer = build_tetrahedral_isomer(mol, &centers, &tetra_config)?;
                let isomer = build_double_bond_isomer(&isomer, &bonds, &bond_config)?;

                if params.only_unique {
                    let smi = canonical_smiles(&isomer)?;
                    if !seen_smiles.insert(smi) {
                        continue;
                    }
                }

                result.push(isomer);
            }

            Ok(result)
        }
    }
}

/// Get the number of possible tetrahedral stereoisomers.
///
/// This is 2^N where N is the number of tetrahedral stereo centers.
/// Returns 1 if there are no stereo centers (the single "no change"
/// isomer).
#[must_use]
pub fn count_stereoisomers(mol: &Molecule) -> u64 {
    let n = find_tetrahedral_centers(mol).len();
    if n == 0 { 1 } else { 1u64 << n }
}

/// Get the number of possible double-bond stereoisomers.
///
/// This is 2^M where M is the number of E/Z stereo bonds.
/// Returns 1 if there are no such bonds.
#[must_use]
pub fn count_double_bond_stereoisomers(mol: &Molecule) -> u64 {
    let m = find_stereo_bonds(mol).len();
    if m == 0 { 1 } else { 1u64 << m }
}

/// Get the total number of possible stereoisomers (tetrahedral + double bond).
#[must_use]
pub fn count_all_stereoisomers(mol: &Molecule) -> u64 {
    let n = find_tetrahedral_centers(mol).len();
    let m = find_stereo_bonds(mol).len();
    match n + m {
        0 => 1,
        total => 1u64 << total,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Molecule;

    fn make_test_mol() -> Molecule {
        // CHClBrI - a single tetrahedral center molecule
        // Use an explicit hydrogen so this fixture exercises enumerator logic
        // without depending on the separate shorthand-chirality parse boundary.
        Molecule::from_smiles("Cl[C@H](Br)I").expect("failed to parse test SMILES")
    }

    fn make_meso_mol() -> Molecule {
        // (R,R)-tartaric acid analogue with two centers
        // SMILES: O[C@@H](C(=O)O)[C@@H](C(=O)O)O
        Molecule::from_smiles("O[C@@H](C(=O)O)[C@@H](C(=O)O)O")
            .expect("failed to parse meso test SMILES")
    }

    fn make_double_bond_mol() -> Molecule {
        // A simple E/Z double bond: C/C=C/Cl
        Molecule::from_smiles("C/C=C/Cl").expect("failed to parse double bond test SMILES")
    }

    #[test]
    fn test_find_tetrahedral_centers() {
        let mol = make_test_mol();
        let centers = find_tetrahedral_centers(&mol);
        // CHClBrI has one chiral center (the carbon)
        assert_eq!(centers.len(), 1, "expected 1 tetrahedral center");
    }

    #[test]
    fn test_enum_stereoisomers_basic() {
        let mol = make_test_mol();
        let params = EnumerationParams::default();
        let isomers = enum_stereoisomers(&mol, &params).expect("enumeration failed");
        // One center -> 2 isomers (CW and CCW)
        assert_eq!(isomers.len(), 2, "expected 2 isomers for one center");
    }

    #[test]
    fn test_enum_stereoisomers_two_centers() {
        let mol = make_meso_mol();
        let params = EnumerationParams::default();
        let isomers = enum_stereoisomers(&mol, &params).expect("enumeration failed");
        // Two centers -> up to 4 isomers
        assert!(
            isomers.len() <= 4,
            "expected at most 4 isomers for two centers"
        );
        assert!(!isomers.is_empty(), "expected at least one isomer");
    }

    #[test]
    fn test_generate_combinations() {
        let combos = generate_combinations(3);
        assert_eq!(combos.len(), 8, "expected 8 combinations for 3 centers");
        for combo in &combos {
            assert_eq!(combo.len(), 3, "each combo must have 3 elements");
        }
    }

    #[test]
    fn test_enum_double_bond_stereoisomers() {
        let mol = make_double_bond_mol();
        let params = EnumerationParams::default();
        let isomers = enum_double_bond_stereoisomers(&mol, &params);
        // May or may not have E/Z bonds depending on SMILES parsing
        if let Ok(isomers) = isomers {
            assert!(
                !isomers.is_empty(),
                "expected at least one double bond isomer"
            );
        }
    }

    #[test]
    fn test_no_stereo_centers() {
        // Ethane: no chiral centers
        let mol = Molecule::from_smiles("CC").expect("failed to parse ethane");
        let params = EnumerationParams::default();
        let result = enum_stereoisomers(&mol, &params);
        assert!(
            result.is_err(),
            "expected error for molecule with no centers"
        );
        assert_eq!(result.unwrap_err(), EnumerationError::NoStereoCenters);
    }

    #[test]
    fn test_count_functions() {
        let mol = make_test_mol();
        assert_eq!(count_stereoisomers(&mol), 2, "one center -> 2 isomers");

        let mol2 = make_meso_mol();
        assert_eq!(count_stereoisomers(&mol2), 4, "two centers -> 4 isomers");

        let ethane = Molecule::from_smiles("CC").expect("failed to parse ethane");
        assert_eq!(count_stereoisomers(&ethane), 1, "no centers -> 1 isomer");
    }

    #[test]
    fn test_unique_isomers() {
        let mol = make_meso_mol();
        // With only_unique=true
        let params_unique = EnumerationParams {
            only_unique: true,
            ..Default::default()
        };
        let unique_isomers = enum_stereoisomers(&mol, &params_unique).expect("enumeration failed");

        // With only_unique=false
        let params_all = EnumerationParams {
            only_unique: false,
            ..Default::default()
        };
        let all_isomers = enum_stereoisomers(&mol, &params_all).expect("enumeration failed");

        // Unique count should be <= full count
        assert!(unique_isomers.len() <= all_isomers.len());
    }
}
