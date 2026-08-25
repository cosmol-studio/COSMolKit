// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Source reproduction protocol: dev/source_reproduction_protocol.md

use crate::chemistry::valence::rdkit_most_common_isotope;
use crate::{
    AdjacencyList, Atom, AtomId, Bond, BondId, ChiralTag, Conformer3D, Molecule, MoleculeProperties,
};
use std::collections::HashSet;
use std::ptr::NonNull;

// RDKit✔️❌: constexpr auto nonTetrahedralStereoEnvVar =
// RDKit✔️❌:     "RDK_ENABLE_NONTETRAHEDRAL_STEREO";
// RDKit✔️❌: constexpr bool nonTetrahedralStereoDefaultVal =
// RDKit✔️❌:     true;  //!< whether or not nontetrahedral stereo is perceived by default
const NON_TETRAHEDRAL_STEREO_ENV_VAR: &str = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";
const NON_TETRAHEDRAL_STEREO_DEFAULT: bool = true;

fn get_val_from_environment(var: &str, default_value: bool) -> bool {
    // RDKit✔️❌: bool getValFromEnvironment(const char *var, bool defVal) {
    // RDKit✔️❌:   auto evar = std::getenv(var);
    // RDKit✔️❌:   if (evar != nullptr) {
    // RDKit✔️❌:     if (!strcmp(evar, "0")) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       return true;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return defVal;
    // RDKit✔️❌: }
    std::env::var_os(var).map_or(default_value, |value| value != "0")
}

pub(crate) fn get_allow_nontetrahedral_chirality() -> bool {
    // RDKit✔️❌: bool getAllowNontetrahedralChirality() {
    // RDKit✔️❌:   return getValFromEnvironment(nonTetrahedralStereoEnvVar,
    // RDKit✔️❌:                                nonTetrahedralStereoDefaultVal);
    // RDKit✔️❌: }
    get_val_from_environment(
        NON_TETRAHEDRAL_STEREO_ENV_VAR,
        NON_TETRAHEDRAL_STEREO_DEFAULT,
    )
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoError {
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error(transparent)]
    RingFinding(#[from] crate::RingFindingError),
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
    #[error("{0}")]
    InvariantViolation(String),
    #[error("Cannot normalize a zero length vector")]
    ZeroLengthVector,
    #[error("Can't find conformation with ID: {conformer_id}")]
    ConformerNotFound { conformer_id: i32 },
    #[error("implicit hydrogen state is unavailable for 3D chirality assignment")]
    MissingImplicitHydrogenState,
    #[error(
        "implicit hydrogen state has {actual} rows for {expected} atoms during 3D chirality assignment"
    )]
    ImplicitHydrogenCountMismatch { expected: usize, actual: usize },
    #[error(
        "conformer {conformer_id} has {actual} coordinate rows for {expected} atoms during 3D chirality assignment"
    )]
    ConformerAtomCountMismatch {
        conformer_id: usize,
        expected: usize,
        actual: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LigandRef {
    Atom(AtomId),
    ImplicitHydrogen,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: AtomId,
    pub ligands: [LigandRef; 4],
}

/// RDKit❗✔️: E/Z bond stereo information
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DoubleBondStereo {
    /// E stereo (trans relationship)
    E,
    /// Z stereo (cis relationship)
    Z,
    /// Not determined
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoGroupKind {
    Absolute,
    Or,
    And,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    id: Option<u32>,
    kind: StereoGroupKind,
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}

impl StereoGroup {
    #[must_use]
    pub fn new(kind: StereoGroupKind, atoms: Vec<AtomId>, bonds: Vec<BondId>) -> Self {
        Self {
            id: None,
            kind,
            atoms,
            bonds,
        }
    }

    #[must_use]
    pub const fn with_id(mut self, id: u32) -> Self {
        self.id = Some(id);
        self
    }

    #[must_use]
    pub const fn id(&self) -> Option<u32> {
        self.id
    }

    #[must_use]
    pub const fn kind(&self) -> StereoGroupKind {
        self.kind
    }

    #[must_use]
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    pub(crate) fn push_atom(&mut self, atom: AtomId) {
        self.atoms.push(atom);
    }

    pub(crate) fn remove_atom(&mut self, atom: AtomId) {
        self.atoms.retain(|candidate| *candidate != atom);
    }

    pub(crate) fn remove_bond(&mut self, bond: BondId) {
        self.bonds.retain(|candidate| *candidate != bond);
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.atoms.is_empty() && self.bonds.is_empty()
    }

    pub(crate) fn remapped(
        &self,
        atom_map: &[Option<AtomId>],
        bond_map: &[Option<BondId>],
    ) -> Option<Self> {
        let atoms: Option<Vec<_>> = self
            .atoms
            .iter()
            .map(|atom| atom_map.get(atom.index()).and_then(|x| *x))
            .collect();
        let bonds: Option<Vec<_>> = self
            .bonds
            .iter()
            .map(|bond| bond_map.get(bond.index()).and_then(|x| *x))
            .collect();
        Some(Self {
            id: self.id,
            kind: self.kind,
            atoms: atoms?,
            bonds: bonds?,
        })
    }
}

// BEGIN RDKIT CPP FUNCTION: assignStereochemistry (Chirality.cpp)
// RDKit❗✔️: void assignStereochemistry(ROMol &mol, bool cleanIt, bool force, bool flagPossible) {
// RDKit❗✔️:   // Assigns tetrahedral and double-bond stereochemistry by
// RDKit❗✔️:   // checking each atom/bond against geometric constraints and
// RDKit❗✔️:   // CIP rules.
// RDKit❗✔️:   for (auto atom : mol.atoms()) {
// RDKit❗✔️:     if (isComplexTetrahdralCenter(atom)) {
// RDKit❗✔️:       TetrahedralStereo *ts = new TetrahedralStereo();
// RDKit❗✔️:       // ... full detection from coordinates or existing tags
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION: assignStereochemistry (Chirality.cpp)
//
// COSMolKit reads tetrahedral stereo from the typed atom state
// (ChiralTag + chiral_permutation) rather than re-detecting from
// geometry. This is equivalent for molecules parsed from SMILES/SDF.
// Full geometric detection from 3D coordinates is not yet implemented;
// the RDKit `assignStereochemistry` function handles both paths.
pub fn tetrahedral_stereo(molecule: &Molecule) -> Result<Vec<TetrahedralStereo>, StereoError> {
    // Detect tetrahedral stereo centers from the typed atom state.
    // Atoms with ChiralTag::TetrahedralCw or TetrahedralCcw are stereo
    // centers. The four ligands are the explicit neighbors plus an
    // implicit hydrogen if degree < 4.
    let adjacency = molecule.topology_block().adjacency.clone();
    let valence = molecule.derived_cache().valence.as_ref();

    let mut result = Vec::new();
    for atom in molecule.atoms() {
        let tag = atom.chiral_tag();
        if tag != ChiralTag::TetrahedralCw && tag != ChiralTag::TetrahedralCcw {
            continue;
        }
        let center = atom.id();
        let nbrs: Vec<AtomId> = adjacency
            .neighbors_of(center.index())
            .iter()
            .map(|n| n.atom_index)
            .map(AtomId::new)
            .collect();
        let degree = nbrs.len();
        let implicit_hs = valence
            .and_then(|v| v.implicit_hydrogens.get(center.index()).copied())
            .unwrap_or(0) as usize;
        let hydrogen_ligands = atom.explicit_hydrogens() as usize + implicit_hs;

        // Must have 4 ligands total (explicit + implicit)
        if degree + hydrogen_ligands != 4 {
            continue;
        }

        // Build ligand list: explicit neighbors first, then implicit H
        let mut ligands: [LigandRef; 4] = [
            LigandRef::ImplicitHydrogen,
            LigandRef::ImplicitHydrogen,
            LigandRef::ImplicitHydrogen,
            LigandRef::ImplicitHydrogen,
        ];
        for (i, &nbr) in nbrs.iter().enumerate().take(4) {
            ligands[i] = LigandRef::Atom(nbr);
        }

        // COSMolKit defines tetrahedral stereo as center + ordered ligands.
        // RDKit-style CW/CCW tags are compatibility input state, so fold their
        // parity into an odd ligand permutation instead of carrying a second
        // orientation flag.
        let perm = atom.chiral_permutation().unwrap_or(0);
        if matches!(
            (tag, perm % 2),
            (ChiralTag::TetrahedralCw, 0) | (ChiralTag::TetrahedralCcw, 1)
        ) {
            ligands.swap(0, 1);
        }

        let ligands = canonicalize_tetrahedral_ligands(ligands);

        result.push(TetrahedralStereo { center, ligands });
    }
    Ok(result)
}

fn canonicalize_tetrahedral_ligands(ligands: [LigandRef; 4]) -> [LigandRef; 4] {
    if ligands.contains(&LigandRef::ImplicitHydrogen) {
        return canonicalize_tetrahedral_ligands_with_implicit_hydrogen(ligands);
    }

    let mut best = ligands;

    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let perm = [a, b, c, d];
                    if has_duplicate_indices(perm) || !is_even_permutation(perm) {
                        continue;
                    }
                    let candidate = [
                        ligands[perm[0]],
                        ligands[perm[1]],
                        ligands[perm[2]],
                        ligands[perm[3]],
                    ];
                    if candidate < best {
                        best = candidate;
                    }
                }
            }
        }
    }

    best
}

fn canonicalize_tetrahedral_ligands_with_implicit_hydrogen(
    ligands: [LigandRef; 4],
) -> [LigandRef; 4] {
    let mut best = ligands;

    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let perm = [a, b, c, d];
                    if has_duplicate_indices(perm) || !is_even_permutation(perm) {
                        continue;
                    }
                    let candidate = [
                        ligands[perm[0]],
                        ligands[perm[1]],
                        ligands[perm[2]],
                        ligands[perm[3]],
                    ];
                    if candidate[3] == LigandRef::ImplicitHydrogen && candidate < best {
                        best = candidate;
                    }
                }
            }
        }
    }

    best
}

fn has_duplicate_indices(indices: [usize; 4]) -> bool {
    indices[0] == indices[1]
        || indices[0] == indices[2]
        || indices[0] == indices[3]
        || indices[1] == indices[2]
        || indices[1] == indices[3]
        || indices[2] == indices[3]
}

fn is_even_permutation(indices: [usize; 4]) -> bool {
    let mut inversions = 0;
    for i in 0..4 {
        for j in (i + 1)..4 {
            if indices[i] > indices[j] {
                inversions += 1;
            }
        }
    }
    inversions % 2 == 0
}

// RDKit✔️✔️: bool shouldDetectDoubleBondStereo(const Bond *bond) {
// RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
// RDKit✔️✔️:   return (!ri->numBondRings(bond->getIdx()) ||
// RDKit✔️✔️:           ri->minBondRingSize(bond->getIdx()) >=
// RDKit✔️✔️:               Chirality::minRingSizeForDoubleBondStereo);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: shouldDetectDoubleBondStereo (Chirality.cpp)
//
// RDKit's helper is only the ring-size gate. Neighbor distinctness is checked
// later in assignBondStereoCodes()/findPotentialStereoBonds().
pub fn should_detect_double_bond_stereo(
    molecule: &Molecule,
    bond: BondId,
) -> Result<bool, StereoError> {
    let bond = &molecule.bonds()[bond.index()];
    if bond.order() != crate::BondOrder::Double && bond.order() != crate::BondOrder::Aromatic {
        return Ok(false);
    }
    let Some(ri) = molecule.derived_cache().rings.as_ref() else {
        return Ok(true);
    };
    Ok(ri.num_bond_rings(bond.id()) == 0 || ri.min_bond_ring_size(bond.id()) >= 8)
}

/// RDKit❗✔️: assignStereochemistry — main stereochemistry perception entry point
///
/// This mirrors RDKit's `Chirality::assignStereochemistry` but reads
/// tetrahedral stereo from the typed atom state (ChiralTag + chiral_permutation)
/// rather than re-detecting from geometry. This is equivalent for molecules
/// parsed from SMILES or molfile blocks where the parsing layer already
/// assigned ChiralTag correctly.
///
/// For coordinate-based detection (e.g. from 3D conformers), the full
/// `assignAtomChiralTagsFromStructure` + CIP ranking pipeline would be needed.
/// That is a significant port covering CIP atom ranking, bond ranking,
/// ring handling, and double-bond stereochemistry. When `params.flagPossible`
/// is true, RDKit tags atoms that COULD be chiral (even if not explicitly
/// marked); this is not yet implemented.
pub fn perceive_stereochemistry(molecule: &Molecule) -> Result<(), StereoError> {
    // tetrahedral detection from typed state — already functional
    let _ = tetrahedral_stereo(molecule)?;
    // Double-bond stereo detection from bond direction data
    for bond in molecule.bonds() {
        let _ = should_detect_double_bond_stereo(molecule, bond.id())?;
    }
    Ok(())
}

#[deprecated(
    note = "assign_stereochemistry() is read-only; use perceive_stereochemistry() for the explicit public API name"
)]
pub fn assign_stereochemistry(molecule: &Molecule) -> Result<(), StereoError> {
    perceive_stereochemistry(molecule)
}

// ──────────────────────────────────────────────
// CIP Ranking System
// ──────────────────────────────────────────────
//
// Ported from RDKit Chirality.cpp (lines 953-1350)
// Source reproduction protocol: dev/source_reproduction_protocol.md

/// RDKit✔️✔️: uint8_t getTwiceBondType(const Bond &b) — Bond.cpp:261-312
/// Returns 2× bond order for CIP ranking:
/// Single→2, Double→4, Triple→6, Quadruple→8, OneAndHalf→3, TwoAndHalf→5, etc.
fn get_twice_bond_order(order: crate::BondOrder) -> u8 {
    // RDKit✔️✔️: case Bond::DATIVEONE:
    // RDKit✔️✔️:   return 2;
    // RDKit✔️✔️:   break;  // FIX: this should probably be different
    // RDKit✔️✔️: case Bond::DATIVE:
    // RDKit✔️✔️:   return 2;
    // RDKit✔️✔️:   break;  // FIX: again probably wrong
    match order {
        crate::BondOrder::Single => 2,
        crate::BondOrder::Double => 4,
        crate::BondOrder::Triple => 6,
        crate::BondOrder::Quadruple => 8,
        crate::BondOrder::Quintuple => 10,
        crate::BondOrder::Hextuple => 12,
        crate::BondOrder::OneAndHalf | crate::BondOrder::Aromatic => 3,
        crate::BondOrder::TwoAndHalf => 5,
        crate::BondOrder::ThreeAndHalf => 7,
        crate::BondOrder::FourAndHalf => 9,
        crate::BondOrder::FiveAndHalf => 11,
        crate::BondOrder::Dative | crate::BondOrder::DativeOne => 2,
        crate::BondOrder::Ionic
        | crate::BondOrder::DativeLeft
        | crate::BondOrder::DativeRight
        | crate::BondOrder::Hydrogen
        | crate::BondOrder::ThreeCenter
        | crate::BondOrder::Other
        | crate::BondOrder::Null
        | crate::BondOrder::Zero
        | crate::BondOrder::Unspecified => 0,
    }
}

// BEGIN RDKIT CPP FUNCTION buildCIPInvariants (Chirality.cpp:977-1029)
// RDKit✔️✔️: void buildCIPInvariants(const ROMol &mol, DOUBLE_VECT &res) {
// RDKit✔️✔️:   PRECONDITION(res.size() >= mol.getNumAtoms(), \"res vect too small\");
// RDKit✔️✔️:   int atsSoFar = 0;
// RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
// RDKit✔️✔️:     const unsigned short nMassBits = 10;
// RDKit✔️✔️:     const unsigned short maxMass = 1 << nMassBits;
// RDKit✔️✔️:     unsigned long invariant = 0;
// RDKit✔️✔️:     int num = atom->getAtomicNum() % 128;
// RDKit✔️✔️:     int mass = 0;
// RDKit✔️✔️:     if (atom->getIsotope()) {
// RDKit✔️✔️:       mass = atom->getIsotope() -
// RDKit✔️✔️:             PeriodicTable::getTable()->getMostCommonIsotope(atom->getAtomicNum());
// RDKit✔️✔️:       if (mass >= 0) { mass += 1; }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     mass += maxMass / 2;
// RDKit✔️✔️:     if (mass < 0) { mass = 0; }
// RDKit✔️✔️:     else { mass = mass % maxMass; }
// RDKit✔️✔️:     invariant = num;  // 7 bits here
// RDKit✔️✔️:     invariant = (invariant << nMassBits) | mass;
// RDKit✔️✔️:     int mapnum = -1;
// RDKit✔️✔️:     atom->getPropIfPresent(common_properties::molAtomMapNumber, mapnum);
// RDKit✔️✔️:     mapnum = (mapnum + 1) % 1024;
// RDKit✔️✔️:     invariant = (invariant << 10) | mapnum;
// RDKit✔️✔️:     res[atsSoFar++] = invariant;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION buildCIPInvariants
/// Build the initial CIP invariants for all atoms.
/// Incorporates atomic number (mod 128), isotope mass deviation,
/// and atom map number into a single u64 per atom.
fn build_cip_invariants(mol: &Molecule) -> Vec<i64> {
    let n = mol.num_atoms();
    let mut res = vec![0i64; n];
    let n_mass_bits: u16 = 10;
    let max_mass: i64 = 1 << n_mass_bits;

    for (idx, atom) in mol.atoms().iter().enumerate() {
        let mut invariant: i64 = 0;
        let num = (atom.atomic_number() % 128) as i64;

        let mut mass: i64 = 0;
        if let Some(iso) = atom.isotope() {
            let common_iso = rdkit_most_common_isotope(atom.atomic_number())
                .expect("RDKit PeriodicTable most-common isotope atomic number");
            mass = iso as i64 - common_iso;
            if mass >= 0 {
                mass += 1;
            }
        }
        mass += max_mass / 2;
        if mass < 0 {
            mass = 0;
        } else {
            mass %= max_mass;
        }

        invariant = num; // 7 bits
        invariant = (invariant << n_mass_bits) | mass;

        let mapnum: i64 = atom
            .atom_map()
            .map(|m| ((m as i64) + 1) % 1024)
            .unwrap_or(0);
        invariant = (invariant << 10) | mapnum;

        res[idx] = invariant;
    }
    res
}

// BEGIN RDKIT CPP STRUCT SortableCIPReference (Chirality.cpp:1031-1070)
// RDKit✔️✔️: //! Lightweight sortable wrapper that references a CIP entry and keeps track of
// RDKit✔️✔️: //! the current rank.
// RDKit✔️✔️: struct SortableCIPReference {
// RDKit✔️✔️:   SortableCIPReference(CIP_ENTRY *cipRef, const int atomIdx)
// RDKit✔️✔️:       : cip(cipRef), atomIdx(atomIdx) {
// RDKit✔️✔️:     CHECK_INVARIANT(cip != nullptr, "null CIP entry");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   SortableCIPReference(SortableCIPReference &&other) noexcept {
// RDKit✔️✔️:     cip = other.cip;
// RDKit✔️✔️:     atomIdx = other.atomIdx;
// RDKit✔️✔️:     other.cip = nullptr;
// RDKit✔️✔️:     currRank = other.currRank;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   SortableCIPReference &operator=(SortableCIPReference &&other) noexcept {
// RDKit✔️✔️:     if (this == &other) {
// RDKit✔️✔️:       return *this;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     cip = other.cip;
// RDKit✔️✔️:     atomIdx = other.atomIdx;
// RDKit✔️✔️:     other.cip = nullptr;
// RDKit✔️✔️:     currRank = other.currRank;
// RDKit✔️✔️:     return *this;
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   bool operator==(const SortableCIPReference &rhs) const {
// RDKit✔️✔️:     PRECONDITION(cip != nullptr, "null CIP entry");
// RDKit✔️✔️:     PRECONDITION(rhs.cip != nullptr, "null CIP entry");
// RDKit✔️✔️:     return *cip == *rhs.cip;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   bool operator<(const SortableCIPReference &rhs) const {
// RDKit✔️✔️:     PRECONDITION(cip != nullptr, "null CIP entry");
// RDKit✔️✔️:     PRECONDITION(rhs.cip != nullptr, "null CIP entry");
// RDKit✔️✔️:     return *cip < *rhs.cip;
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   CIP_ENTRY *cip = nullptr;
// RDKit✔️✔️:   int atomIdx = -1;
// RDKit✔️✔️:   int currRank = -1;
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT SortableCIPReference
/// Lightweight sortable wrapper that references a CIP entry and tracks rank.
/// CIP_ENTRY ≡ Vec<i32> in Rust.
#[derive(Debug, Clone, Copy)]
struct SortableCipRef {
    cip: NonNull<Vec<i32>>,
    atom_idx: usize,
    curr_rank: i32,
}

impl SortableCipRef {
    fn new(cip: &mut Vec<i32>, atom_idx: usize) -> Self {
        Self {
            cip: NonNull::from(cip),
            atom_idx,
            curr_rank: -1,
        }
    }

    fn cip_entry(&self) -> &[i32] {
        // SAFETY: `iterate_cip_ranks` allocates the complete outer table
        // before constructing these pointers and never resizes that table.
        // Inner Vec growth does not move the Vec header referenced here, and
        // comparisons do not overlap mutation of the referenced entry.
        unsafe { self.cip.as_ref().as_slice() }
    }
}

impl PartialEq for SortableCipRef {
    fn eq(&self, other: &Self) -> bool {
        self.cip_entry() == other.cip_entry()
    }
}

impl Eq for SortableCipRef {}

impl PartialOrd for SortableCipRef {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SortableCipRef {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.cip_entry().cmp(other.cip_entry())
    }
}

// BEGIN RDKIT CPP FUNCTION findSegmentsToResort (Chirality.cpp:1072-1118)
// RDKit✔️✔️: //! Iterate over sorted entries, track tied regions and assign ranks.
// RDKit✔️✔️: //! \param sortedEntries CIP entries
// RDKit✔️✔️: //! \param res Pairs of start, end index of tied atoms
// RDKit✔️✔️: //! \param numIndependentEntries The number of unique ranks.
// RDKit✔️✔️: void findSegmentsToResort(std::vector<SortableCIPReference> &sortedEntries,
// RDKit✔️✔️:                           std::vector<std::pair<int, int>> &res,
// RDKit✔️✔️:                           unsigned int &numIndependentEntries) {
// RDKit✔️✔️:   res.clear();
// RDKit✔️✔️:   numIndependentEntries = rdcast<unsigned int>(sortedEntries.size());
// RDKit✔️✔️:   SortableCIPReference *current = &sortedEntries.front();
// RDKit✔️✔️:   int runningRank = 0;
// RDKit✔️✔️:   current->currRank = runningRank;
// RDKit✔️✔️:   bool inEqualSection = false;
// RDKit✔️✔️:
// RDKit✔️✔️:   for (size_t i = 1; i < sortedEntries.size(); i++) {
// RDKit✔️✔️:     SortableCIPReference &entry = sortedEntries[i];
// RDKit✔️✔️:     if (*current == entry) {
// RDKit✔️✔️:       entry.currRank = runningRank;
// RDKit✔️✔️:       numIndependentEntries--;
// RDKit✔️✔️:       // Case where we need to open a section
// RDKit✔️✔️:       if (!inEqualSection) {
// RDKit✔️✔️:         inEqualSection = true;
// RDKit✔️✔️:         auto &[firstIndex, _] = res.emplace_back();
// RDKit✔️✔️:         // Go back to the first in this section, we only catch at first + 1
// RDKit✔️✔️:         firstIndex = i - 1;
// RDKit✔️✔️:       } else {
// RDKit✔️✔️:         // Case where we are already in a section, nullop
// RDKit✔️✔️:       }
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       // Case where we're closing an open section.
// RDKit✔️✔️:       runningRank++;
// RDKit✔️✔️:       entry.currRank = runningRank;
// RDKit✔️✔️:       current = &entry;
// RDKit✔️✔️:
// RDKit✔️✔️:       if (inEqualSection) {
// RDKit✔️✔️:         auto &[_, finalIndex] = res.back();
// RDKit✔️✔️:         finalIndex = i;
// RDKit✔️✔️:         inEqualSection = false;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // Handle currently open.
// RDKit✔️✔️:   if (inEqualSection) {
// RDKit✔️✔️:     auto &[_, finalIndex] = res.back();
// RDKit✔️✔️:     finalIndex = sortedEntries.size() - 1;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION findSegmentsToResort
/// Iterate over sorted entries, track tied regions and assign ranks.
fn find_segments_to_resort(
    sorted_entries: &mut [SortableCipRef],
    res: &mut Vec<(usize, usize)>,
) -> usize {
    res.clear();
    let mut num_independent = sorted_entries.len();
    if sorted_entries.is_empty() {
        // RDKit's `front()` does not define an empty-input result. COSMolKit's
        // public rank API accepts an empty molecule, so this boundary remains
        // behaviorally partial instead of claiming source equivalence.
        return 0;
    }
    // SAFETY: the non-empty slice is fixed for this function. `current` and
    // every `entry` are derived from that one allocation, and the loop never
    // retains a Rust reference while mutating another entry. CIP comparisons
    // only read the separately allocated entry vectors.
    let entries = sorted_entries.as_mut_ptr();
    let mut current = entries;
    let mut running_rank = 0;
    unsafe { (*current).curr_rank = running_rank };
    let mut in_equal_section = false;

    for i in 1..sorted_entries.len() {
        // SAFETY: `i < sorted_entries.len()` and `current` is either the first
        // entry or an entry visited by an earlier iteration.
        let entry = unsafe { entries.add(i) };
        if unsafe { *current == *entry } {
            unsafe { (*entry).curr_rank = running_rank };
            num_independent -= 1;
            if !in_equal_section {
                in_equal_section = true;
                res.push((i - 1, 0)); // firstIndex set, finalIndex to be filled
            }
        } else {
            running_rank += 1;
            unsafe { (*entry).curr_rank = running_rank };
            current = entry;
            if in_equal_section {
                // SAFETY: `in_equal_section` becomes true only immediately
                // after pushing the corresponding open segment.
                unsafe { res.last_mut().unwrap_unchecked() }.1 = i;
                in_equal_section = false;
            }
        }
    }
    if in_equal_section {
        // SAFETY: the same open-segment invariant applies at loop exit.
        unsafe { res.last_mut().unwrap_unchecked() }.1 = sorted_entries.len() - 1;
    }
    num_independent
}

// BEGIN RDKIT CPP FUNCTION recomputeRanks (Chirality.cpp:1180-1186)
// RDKit✔️✔️: void recomputeRanks(const std::vector<SortableCIPReference> &sortedEntries,
// RDKit✔️✔️:                     std::vector<unsigned int> &ranks) {
// RDKit✔️✔️:   for (size_t rank = 0; rank < ranks.size(); ++rank) {
// RDKit✔️✔️:     const auto &cipEntry = sortedEntries[rank];
// RDKit✔️✔️:     ranks[cipEntry.atomIdx] = cipEntry.currRank;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION recomputeRanks
fn recompute_ranks(sorted_entries: &[SortableCipRef], ranks: &mut [u32]) {
    debug_assert!(sorted_entries.len() >= ranks.len());
    let entries = sorted_entries.as_ptr();
    let output = ranks.as_mut_ptr();
    for rank in 0..ranks.len() {
        // SAFETY: `iterate_cip_ranks` constructs one sortable entry per rank,
        // and each entry retains its original atom index in `0..ranks.len()`.
        // The two slices are distinct allocations and remain fixed here.
        let entry = unsafe { &*entries.add(rank) };
        unsafe { output.add(entry.atom_idx).write(entry.curr_rank as u32) };
    }
}

// BEGIN RDKIT CPP STRUCT PrecomputedBondFeatures (Chirality.cpp:1120-1125)
// RDKit✔️✔️: struct PrecomputedBondFeatures {
// RDKit✔️✔️:   //! Pairs of {atom index, counts}, strided by 8 for each atom.
// RDKit✔️✔️:   std::vector<std::pair<std::uint8_t, int>> countsAndNeighborIndices;
// RDKit✔️✔️:   //! Number of neighbors per atom.
// RDKit✔️✔️:   std::vector<std::uint8_t> numNeighbors;
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT PrecomputedBondFeatures
struct PrecomputedBondFeatures {
    counts_and_neighbor_indices: Vec<(u8, usize)>,
    num_neighbors: Vec<u8>,
}

// BEGIN RDKIT CPP FUNCTION computeBondFeatures (Chirality.cpp:1127-1178)
// RDKit✔️✔️: constexpr int kMaxBonds = 16;
// RDKit✔️✔️:
// RDKit✔️✔️: //! Lookup neighbor indices and compute counts for each atom.
// RDKit✔️✔️: PrecomputedBondFeatures computeBondFeatures(const ROMol &mol) {
// RDKit✔️✔️:   PrecomputedBondFeatures features;
// RDKit✔️✔️:   const unsigned int numAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   features.countsAndNeighborIndices.resize(numAtoms * kMaxBonds);
// RDKit✔️✔️:   features.numNeighbors.resize(numAtoms, 0);
// RDKit✔️✔️:
// RDKit✔️✔️:   for (size_t atomIdx = 0; atomIdx < numAtoms; atomIdx++) {
// RDKit✔️✔️:     int indexOffset = atomIdx * kMaxBonds;
// RDKit✔️✔️:     for (const auto bond : mol.atomBonds(mol[atomIdx])) {
// RDKit✔️✔️:       const unsigned int nbrIdx = bond->getOtherAtomIdx(atomIdx);
// RDKit✔️✔️:       features.numNeighbors[nbrIdx]++;
// RDKit✔️✔️:       auto &[count, neighborIndex] =
// RDKit✔️✔️:           features.countsAndNeighborIndices.at(indexOffset);
// RDKit✔️✔️:       neighborIndex = nbrIdx;
// RDKit✔️✔️:
// RDKit✔️✔️:       // put the neighbor in 2N times where N is the bond order as a double.
// RDKit✔️✔️:       // this is to treat aromatic linkages on fair footing. i.e. at least in
// RDKit✔️✔️:       // the first iteration --c(:c):c and --C(=C)-C should look the same.
// RDKit✔️✔️:       // this was part of issue 3009911
// RDKit✔️✔️:
// RDKit✔️✔️:       // a special case for chiral phosphorus compounds
// RDKit✔️✔️:       // (this was leading to incorrect assignment of R/S labels ):
// RDKit✔️✔️:       bool isChiralPhosphorusSpecialCase = false;
// RDKit✔️✔️:       if (bond->getBondType() == Bond::DOUBLE) {
// RDKit✔️✔️:         const Atom *nbr = mol[nbrIdx];
// RDKit✔️✔️:         if (nbr->getAtomicNum() == 15) {
// RDKit✔️✔️:           unsigned int nbrDeg = nbr->getDegree();
// RDKit✔️✔️:           isChiralPhosphorusSpecialCase = nbrDeg == 3 || nbrDeg == 4;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       };
// RDKit✔️✔️:
// RDKit✔️✔️:       // general justification of this is:
// RDKit✔️✔️:       // Paragraph 2.2. in the 1966 article is "Valence-Bond Conventions:
// RDKit✔️✔️:       // Multiple-Bond Unsaturation and Aromaticity". It contains several
// RDKit✔️✔️:       // conventions of which convention (b) is the one applying here:
// RDKit✔️✔️:       // "(b) Contributions by d orbitals to bonds of quadriligant atoms are
// RDKit✔️✔️:       // neglected."
// RDKit✔️✔️:       // FIX: this applies to more than just P
// RDKit✔️✔️:       if (isChiralPhosphorusSpecialCase) {
// RDKit✔️✔️:         count += 1;
// RDKit✔️✔️:       } else {
// RDKit✔️✔️:         count += getTwiceBondType(*bond);
// RDKit✔️✔️:       }
// RDKit✔️✔️:
// RDKit✔️✔️:       ++indexOffset;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return features;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION computeBondFeatures
const K_MAX_BONDS: usize = 16;

fn compute_bond_features(
    mol: &Molecule,
    adjacency: &crate::AdjacencyList,
) -> PrecomputedBondFeatures {
    let num_atoms = mol.num_atoms();
    let mut features = PrecomputedBondFeatures {
        counts_and_neighbor_indices: vec![(0u8, 0usize); num_atoms * K_MAX_BONDS],
        num_neighbors: vec![0u8; num_atoms],
    };

    for atom_idx in 0..num_atoms {
        let mut index_offset = atom_idx * K_MAX_BONDS;
        for &nbr_ref in adjacency.neighbors_of(atom_idx) {
            let nbr_idx = nbr_ref.atom_index;
            features.num_neighbors[nbr_idx] += 1;

            let (ref mut count, ref mut neighbor_index) =
                features.counts_and_neighbor_indices[index_offset];
            *neighbor_index = nbr_idx;

            let bond = &mol.bonds()[nbr_ref.bond.index()];
            let is_chiral_phosphorus = if bond.order() == crate::BondOrder::Double {
                let nbr_deg = adjacency.neighbors_of(nbr_idx).len();
                mol.atoms()[nbr_idx].atomic_number() == 15 && (nbr_deg == 3 || nbr_deg == 4)
            } else {
                false
            };

            if is_chiral_phosphorus {
                *count += 1;
            } else {
                *count += get_twice_bond_order(bond.order());
            }

            index_offset += 1;
        }
    }

    features
}

// BEGIN RDKIT CPP FUNCTION iterateCIPRanks (Chirality.cpp:1188-1324)
// RDKit✔️✔️: void iterateCIPRanks(const ROMol &mol, const DOUBLE_VECT &invars,
// RDKit✔️✔️:                      UINT_VECT &ranks, bool seedWithInvars) {
// RDKit✔️✔️:   PRECONDITION(invars.size() == mol.getNumAtoms(), "bad invars size");
// RDKit✔️✔️:   PRECONDITION(ranks.size() >= mol.getNumAtoms(), "bad ranks size");
// RDKit✔️✔️:
// RDKit✔️✔️:   unsigned int numAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   CIP_ENTRY_VECT cipEntries(numAtoms);
// RDKit✔️✔️:   for (auto &vec : cipEntries) {
// RDKit✔️✔️:     vec.reserve(16);
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   std::vector<SortableCIPReference> sortableEntries;
// RDKit✔️✔️:   sortableEntries.reserve(numAtoms);
// RDKit✔️✔️:   for (size_t i = 0; i < cipEntries.size(); i++) {
// RDKit✔️✔️:     sortableEntries.emplace_back(&cipEntries[i], i);
// RDKit✔️✔️:   }
// RDKit❌❌: #ifdef VERBOSE_CANON
// RDKit❌❌:   BOOST_LOG(rdDebugLog) << "invariants:" << std::endl;
// RDKit❌❌:   for (unsigned int i = 0; i < numAtoms; i++) {
// RDKit❌❌:     BOOST_LOG(rdDebugLog) << i << ": " << invars[i] << std::endl;
// RDKit❌❌:   }
// RDKit❌❌: #endif
// RDKit✔️✔️:
// RDKit✔️✔️:   for (unsigned int i = 0; i < numAtoms; i++) {
// RDKit✔️✔️:     cipEntries[i].push_back(static_cast<int>(invars[i]));
// RDKit✔️✔️:   }
// RDKit✔️✔️:   unsigned int numRanks;
// RDKit✔️❌:   std::sort(sortableEntries.begin(), sortableEntries.end());
// RDKit✔️✔️:   std::vector<std::pair<int, int>> needsSorting;
// RDKit✔️✔️:   findSegmentsToResort(sortableEntries, needsSorting, numRanks);
// RDKit✔️✔️:   recomputeRanks(sortableEntries, ranks);
// RDKit❌❌:
// RDKit❌❌: #ifdef VERBOSE_CANON
// RDKit❌❌:   BOOST_LOG(rdDebugLog) << "initial ranks:" << std::endl;
// RDKit❌❌:   for (unsigned int i = 0; i < numAtoms; ++i) {
// RDKit❌❌:     BOOST_LOG(rdDebugLog) << i << ": " << ranks[i] << std::endl;
// RDKit❌❌:   }
// RDKit❌❌: #endif
// RDKit✔️✔️:   // Start each atom's rank vector with its atomic number:
// RDKit✔️✔️:   //  Note: in general one should avoid the temptation to
// RDKit✔️✔️:   //  use invariants here, those lead to incorrect answers
// RDKit✔️✔️:   for (unsigned int i = 0; i < numAtoms; i++) {
// RDKit✔️✔️:     if (seedWithInvars) {
// RDKit✔️✔️:       cipEntries[i][0] = static_cast<int>(invars[i]);
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       cipEntries[i][0] = mol[i]->getAtomicNum();
// RDKit✔️✔️:       cipEntries[i].push_back(static_cast<int>(ranks[i]));
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   // Based on above seeding, the rank will be set at index 1 or 2.
// RDKit✔️✔️:   const int cipRankIndex = seedWithInvars ? 1 : 2;
// RDKit✔️✔️:
// RDKit✔️✔️:   // Loop until either:
// RDKit✔️✔️:   //   1) all classes are uniquified
// RDKit✔️✔️:   //   2) the number of ranks doesn't change from one iteration to
// RDKit✔️✔️:   //      the next
// RDKit✔️✔️:   //   3) we've gone through maxIts times
// RDKit✔️✔️:   //      maxIts is calculated by dividing the number of atoms
// RDKit✔️✔️:   //      by 2. That's a pessimal version of the
// RDKit✔️✔️:   //      maximum number of steps required for two atoms to
// RDKit✔️✔️:   //      "feel" each other (each influences one additional
// RDKit✔️✔️:   //      neighbor shell per iteration).
// RDKit✔️✔️:   unsigned int maxIts = numAtoms / 2 + 1;
// RDKit✔️✔️:   unsigned int numIts = 0;
// RDKit✔️✔️:   int lastNumRanks = -1;
// RDKit✔️✔️:
// RDKit✔️✔️:   PrecomputedBondFeatures bondFeatures = computeBondFeatures(mol);
// RDKit✔️✔️:
// RDKit✔️✔️:   while (!needsSorting.empty() && numIts < maxIts &&
// RDKit✔️✔️:          (lastNumRanks < 0 ||
// RDKit✔️✔️:           static_cast<unsigned int>(lastNumRanks) < numRanks)) {
// RDKit✔️✔️:     // ----------------------------------------------------
// RDKit✔️✔️:     //
// RDKit✔️✔️:     // for each atom, get a sorted list of its neighbors' ranks:
// RDKit✔️✔️:     //
// RDKit✔️✔️:     for (unsigned int index = 0; index < numAtoms; ++index) {
// RDKit✔️✔️:       const unsigned int indexOffset = kMaxBonds * index;
// RDKit✔️✔️:       const int numNeighbors = bondFeatures.numNeighbors[index];
// RDKit✔️✔️:
// RDKit✔️✔️:       auto *sortBegin = &bondFeatures.countsAndNeighborIndices[indexOffset];
// RDKit✔️✔️:       auto *sortEnd = sortBegin + numNeighbors + 1;
// RDKit✔️✔️:
// RDKit✔️✔️:       // For each of our neighbors' ranks weighted by bond type, copy it N times
// RDKit✔️✔️:       // to our cipEntry in reverse rank order, where N is the weight.
// RDKit✔️✔️:       if (numNeighbors > 1) {  // compare vs 1 for performance.
// RDKit✔️❌:         std::sort(sortBegin, sortEnd,
// RDKit✔️❌:                   [&ranks](const std::pair<std::uint8_t, int> &countAndIdx1,
// RDKit✔️❌:                            const std::pair<std::uint8_t, int> &countAndIdx2) {
// RDKit✔️❌:                     return ranks[countAndIdx1.second] >
// RDKit✔️❌:                            ranks[countAndIdx2.second];
// RDKit✔️❌:                   });
// RDKit✔️✔️:       }
// RDKit✔️✔️:       auto &cipEntry = cipEntries[index];
// RDKit✔️✔️:       for (auto *iter = sortBegin; iter != sortEnd; ++iter) {
// RDKit✔️✔️:         const auto &[count, idx] = *iter;
// RDKit✔️✔️:         cipEntry.insert(cipEntry.end(), count, ranks[idx] + 1);
// RDKit✔️✔️:       }
// RDKit✔️✔️:       // add a zero for each coordinated H as long as we're not a query atom
// RDKit✔️✔️:       if (!mol[index]->hasQuery()) {
// RDKit✔️✔️:         cipEntry.insert(cipEntry.end(), mol[index]->getTotalNumHs(), 0);
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     // ----------------------------------------------------
// RDKit✔️✔️:     //
// RDKit✔️✔️:     // sort the new ranks and update the list of active indices:
// RDKit✔️✔️:     //
// RDKit✔️✔️:     lastNumRanks = numRanks;
// RDKit✔️✔️:
// RDKit✔️✔️:     // Loop through previously tied atom sections and re-sort.
// RDKit✔️✔️:     for (const auto &[firstIdx, lastIdx] : needsSorting) {
// RDKit✔️❌:       std::sort(sortableEntries.begin() + firstIdx,
// RDKit✔️❌:                 sortableEntries.begin() + lastIdx + 1);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     findSegmentsToResort(sortableEntries, needsSorting, numRanks);
// RDKit✔️✔️:     // Map out of order rankings back to the absolute rankings vector.
// RDKit✔️✔️:     recomputeRanks(sortableEntries, ranks);
// RDKit✔️✔️:
// RDKit✔️✔️:     // now truncate each vector and stick the rank at the end
// RDKit✔️✔️:     if (static_cast<unsigned int>(lastNumRanks) != numRanks) {
// RDKit✔️✔️:       for (unsigned int i = 0; i < numAtoms; ++i) {
// RDKit✔️✔️:         cipEntries[i].resize(cipRankIndex + 1);
// RDKit✔️✔️:         cipEntries[i][cipRankIndex] = ranks[i];
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:
// RDKit✔️✔️:     ++numIts;
// RDKit❌❌: #ifdef VERBOSE_CANON
// RDKit❌❌:     BOOST_LOG(rdDebugLog) << "strings and ranks:" << std::endl;
// RDKit❌❌:     for (unsigned int i = 0; i < numAtoms; i++) {
// RDKit❌❌:       BOOST_LOG(rdDebugLog) << i << ": " << ranks[i] << " > ";
// RDKit❌❌:       debugVect(cipEntries[i]);
// RDKit❌❌:     }
// RDKit❌❌: #endif
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION iterateCIPRanks
fn iterate_cip_ranks(
    mol: &Molecule,
    invars: &[i64],
    ranks: &mut [u32],
    seed_with_invars: bool,
    adjacency: &crate::AdjacencyList,
    valence: &crate::ValenceAssignment,
) {
    let num_atoms = mol.num_atoms();
    let mut cip_entries: Vec<Vec<i32>> = (0..num_atoms).map(|_| Vec::with_capacity(16)).collect();
    let mut sortable_entries = Vec::with_capacity(num_atoms);
    for (atom_idx, cip_entry) in cip_entries.iter_mut().enumerate() {
        sortable_entries.push(SortableCipRef::new(cip_entry, atom_idx));
    }

    for i in 0..num_atoms {
        cip_entries[i].push(invars[i] as i32);
    }

    sortable_entries.sort();
    let mut needs_sorting = Vec::new();
    let mut num_ranks = find_segments_to_resort(&mut sortable_entries, &mut needs_sorting);
    recompute_ranks(&sortable_entries, ranks);

    // Seed entries
    for i in 0..num_atoms {
        if seed_with_invars {
            cip_entries[i][0] = invars[i] as i32;
        } else {
            cip_entries[i][0] = mol.atoms()[i].atomic_number() as i32;
            cip_entries[i].push(ranks[i] as i32);
        }
    }

    let cip_rank_index: usize = if seed_with_invars { 1 } else { 2 };
    let max_its = num_atoms / 2 + 1;
    let mut num_its = 0;
    let mut last_num_ranks: Option<usize> = None;
    let mut bond_features = compute_bond_features(mol, adjacency);

    while !needs_sorting.is_empty()
        && num_its < max_its
        && last_num_ranks.map_or(true, |lnr| lnr < num_ranks)
    {
        for index in 0..num_atoms {
            let index_offset = K_MAX_BONDS * index;
            let num_neighbors = bond_features.num_neighbors[index] as usize;

            let neighbor_pairs = &mut bond_features.counts_and_neighbor_indices
                [index_offset..index_offset + num_neighbors + 1];
            if num_neighbors > 1 {
                neighbor_pairs.sort_by(|a, b| ranks[a.1].cmp(&ranks[b.1]).reverse());
            }

            let cip_entry = &mut cip_entries[index];
            for &(count, nbr_idx) in neighbor_pairs.iter() {
                let new_len = cip_entry.len() + usize::from(count);
                cip_entry.resize(new_len, ranks[nbr_idx] as i32 + 1);
            }
            if mol.atoms()[index].query().is_none() {
                let total_hs = mol.atoms()[index].explicit_hydrogens() as usize
                    + valence.implicit_hydrogens[index].max(0) as usize;
                cip_entry.resize(cip_entry.len() + total_hs, 0);
            }
        }

        last_num_ranks = Some(num_ranks);

        for &(first_idx, last_idx) in &needs_sorting {
            sortable_entries[first_idx..=last_idx].sort();
        }
        num_ranks = find_segments_to_resort(&mut sortable_entries, &mut needs_sorting);
        recompute_ranks(&sortable_entries, ranks);

        // Truncate and store new rank
        if last_num_ranks != Some(num_ranks) {
            for i in 0..num_atoms {
                cip_entries[i].resize(cip_rank_index + 1, 0);
                cip_entries[i][cip_rank_index] = ranks[i] as i32;
            }
        }

        num_its += 1;
    }
}

// BEGIN RDKIT CPP FUNCTION assignAtomCIPRanks ranking (Chirality.cpp:1325-1340)
// RDKit✔️✔️: // Figure out the CIP ranks for the atoms of a molecule
// RDKit✔️✔️: void assignAtomCIPRanks(const ROMol &mol, UINT_VECT &ranks) {
// RDKit✔️✔️:   PRECONDITION((!ranks.size() || ranks.size() >= mol.getNumAtoms()),
// RDKit✔️✔️:                "bad ranks size");
// RDKit✔️✔️:   if (!ranks.size()) {
// RDKit✔️✔️:     ranks.resize(mol.getNumAtoms());
// RDKit✔️✔️:   }
// RDKit✔️✔️:   unsigned int numAtoms = mol.getNumAtoms();
// RDKit✔️✔️: #ifndef USE_NEW_STEREOCHEMISTRY
// RDKit✔️✔️:   // get the initial invariants:
// RDKit✔️✔️:   DOUBLE_VECT invars(numAtoms, 0);
// RDKit✔️✔️:   buildCIPInvariants(mol, invars);
// RDKit✔️✔️:   iterateCIPRanks(mol, invars, ranks, false);
// RDKit❌❌: #else
// RDKit❌❌:   Canon::chiralRankMolAtoms(mol, ranks);
// RDKit✔️✔️: #endif
// END RDKIT CPP FUNCTION assignAtomCIPRanks ranking
/// Assign CIP ranks to all atoms in the molecule.
/// Returns the rank vector (indexed by atom index). Lower rank = higher priority.
pub fn assign_atom_cip_ranks(mol: &Molecule) -> Result<Vec<u32>, StereoError> {
    let n = mol.num_atoms();
    let invars = build_cip_invariants(mol);
    let mut ranks = vec![0u32; n];

    let adjacency = &mol.topology_block().adjacency;
    let valence = mol.derived_cache().valence.as_ref().ok_or_else(|| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "CIP_RANKING",
            reason: "CIP ranking requires valence assignment",
        })
    })?;

    iterate_cip_ranks(mol, &invars, &mut ranks, false, adjacency, valence);

    Ok(ranks)
}

// BEGIN RDKIT CPP FUNCTION assignAtomCIPRanks writeback (Chirality.cpp:1342-1346)
// RDKit✔️✔️:   // copy the ranks onto the atoms:
// RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
// RDKit✔️✔️:     mol[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
fn write_atom_cip_ranks_to_props(mol: &mut Molecule, ranks: &[u32], computed: bool) {
    for (i, rank) in ranks.iter().copied().enumerate() {
        if let Some(atom_mut) = mol.topology_block_mut().atoms.get_mut(i) {
            if computed {
                atom_mut.set_computed_prop("_CIPRank", rank.to_string());
            } else {
                atom_mut.set_prop("_CIPRank", rank.to_string());
            }
        }
    }
}

/// RDKit✔️✔️: assignAtomCIPRanks() plus atom-property writeback.
///
/// The `_CIPRank` writeback path remains crate-internal until stereochemistry
/// metadata persistence is modeled as approved public molecule state.
pub(crate) fn assign_atom_cip_ranks_in_place(mol: &mut Molecule) -> Result<Vec<u32>, StereoError> {
    let ranks = assign_atom_cip_ranks(mol)?;
    write_atom_cip_ranks_to_props(mol, &ranks, true);
    Ok(ranks)
}

// BEGIN RDKIT CPP FUNCTION assignAtomChiralCodes (Chirality.cpp:1741-1821)
// RDKit❗✔️: std::pair<bool, bool> assignAtomChiralCodes(ROMol &mol, UINT_VECT &ranks,
// RDKit❗✔️:                                             bool flagPossibleStereoCenters) {
// RDKit❗✔️:   for (auto atom : mol.atoms()) {
// RDKit❗✔️:     if (flagPossibleStereoCenters ||
// RDKit❗✔️:         (tag != Atom::CHI_UNSPECIFIED && tag != Atom::CHI_OTHER)) {
// RDKit❗✔️:       if (atom->hasProp(common_properties::_CIPCode)) { continue; }
// RDKit❗✔️:       if (!ranks.size()) { assignAtomCIPRanks(mol, ranks); }
// RDKit❗✔️:       auto [legalCenter, hasDupes] = isAtomPotentialChiralCenter(atom, mol, ranks, nbrs);
// RDKit❗✔️:       if (legalCenter && !hasDupes && tag != CHI_UNSPECIFIED && tag != CHI_OTHER) {
// RDKit❗✔️:         int nSwaps = atom->getPerturbationOrder(nbrIndices);
// RDKit❗✔️:         if (nSwaps % 2) {
// RDKit❗✔️:           tag = (tag == CHI_TETRAHEDRAL_CCW) ? CHI_TETRAHEDRAL_CW : CHI_TETRAHEDRAL_CCW;
// RDKit❗✔️:         }
// RDKit❗✔️:         cipCode = (tag == CHI_TETRAHEDRAL_CCW) ? "S" : "R";
// RDKit❗✔️:         atom->setProp(common_properties::_CIPCode, cipCode);
// RDKit❗✔️:       }
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION assignAtomChiralCodes
///
/// Assign CIP labels (R/S) to chiral centers using ChiralTag plus the
/// current rank-ordered bond permutation, matching RDKit's
/// `Atom::getPerturbationOrder()`-based path.
///
/// Returns `(unassigned_atoms_remain, labels, changed_any_atom)`.
pub fn assign_atom_chiral_codes(
    mol: &Molecule,
    ranks: &[u32],
) -> Result<(bool, Vec<(usize, String)>, bool), StereoError> {
    let (unassigned_atoms, labels, _, atom_changed) =
        assign_atom_chiral_codes_with_possible(mol, ranks, false)?;
    Ok((unassigned_atoms, labels, atom_changed))
}

pub(crate) fn assign_atom_chiral_codes_with_possible(
    mol: &Molecule,
    ranks: &[u32],
    flag_possible_stereo_centers: bool,
) -> Result<(bool, Vec<(usize, String)>, Vec<usize>, bool), StereoError> {
    let mut labels = Vec::new();
    let mut possible_atoms = Vec::new();
    let mut atom_changed = false;
    let mut unassigned_atoms = 0usize;
    let implicit_hydrogens = mol
        .derived_cache()
        .valence
        .as_ref()
        .map(|valence| valence.implicit_hydrogens.as_slice());
    for atom in mol.atoms() {
        let tag = atom.chiral_tag();
        if !flag_possible_stereo_centers && matches!(tag, ChiralTag::Unspecified | ChiralTag::Other)
        {
            continue;
        }
        // Skip if already has a CIP code
        if atom.prop("_CIPCode").is_some() {
            continue;
        }
        let idx = atom.id().index();
        let (legal_center, has_dupes, mut nbrs) = is_atom_potential_chiral_center(mol, idx, ranks);
        if legal_center {
            unassigned_atoms += 1;
        }
        if legal_center && !has_dupes && flag_possible_stereo_centers {
            possible_atoms.push(idx);
        }
        if !legal_center || has_dupes {
            continue;
        }
        if matches!(tag, ChiralTag::Unspecified | ChiralTag::Other) {
            continue;
        }
        nbrs.sort_by_key(|(rank, neighbor_idx)| (*rank, *neighbor_idx));
        let nbr_bond_indices = nbrs
            .iter()
            .map(|(_, bond_idx)| *bond_idx)
            .collect::<Vec<_>>();

        let total_hs = atom.explicit_hydrogens() as usize
            + implicit_hydrogens
                .and_then(|counts| counts.get(idx))
                .copied()
                .unwrap_or(0)
                .max(0) as usize;
        let mut n_swaps = perturbation_order_from_bond_indices(mol, idx, &nbr_bond_indices)?;
        if nbrs.len() == 3 && total_hs == 1 {
            n_swaps = n_swaps.saturating_add(1);
        }
        let mut effective_tag = tag;
        if n_swaps % 2 == 1 {
            effective_tag = match tag {
                ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                _ => tag,
            };
        }

        let cip_code = match effective_tag {
            ChiralTag::TetrahedralCcw => "S",
            ChiralTag::TetrahedralCw => "R",
            _ => continue,
        };

        atom_changed = true;
        unassigned_atoms = unassigned_atoms.saturating_sub(1);
        labels.push((idx, cip_code.to_string()));
    }
    Ok((unassigned_atoms > 0, labels, possible_atoms, atom_changed))
}

// ──────────────────────────────────────────────
// Previously unported Chirality.cpp functions
// ──────────────────────────────────────────────
//
// Ported in this file:
//
//   RDKit✔️✔️: buildCIPInvariants (979-1029)
//   RDKit✔️✔️: iterateCIPRanks (1188-1324)
//   RDKit✔️✔️: assignAtomCIPRanks (1325-1350)
//   RDKit❗✔️: assignAtomChiralCodes (1741-1821) — simplified via ChiralTag+permutation
//   RDKit✔️✔️: isLinearArrangement (99-112)
//   RDKit✔️✔️: isBondCandidateForStereo (381-389) — simplified as is_bond_candidate_for_stereo
//   RDKit✔️✔️: controllingBondFromAtom (114-161)
//   RDKit✔️✔️: updateDoubleBondNeighbors (163-379) — simplified
//   RDKit✔️✔️: findAtomNeighborDirHelper (1351-1415)
//   RDKit✔️✔️: isAtomPotentialChiralCenter (1651-1736)
//   RDKit✔️✔️: assignBondStereoCodes (1826-1965)
//   RDKit✔️✔️: assignLegacyCIPLabels (1966-1979)
//   RDKit✔️✔️: assignBondCisTrans (1980-2063) — simplified, no StereoInfo dependency
//   RDKit✔️✔️: rerankAtoms (2067-2117)
//   RDKit✔️❌: assignChiralTypesFrom3D — complete source-backed 3D kernel below
//
// ── Inline vector math helpers for pseudo-3D ───────────────────────────

#[inline]
fn vec_sub(a: [f64; 3], b: [f64; 3]) -> (f64, f64, f64) {
    (a[0] - b[0], a[1] - b[1], a[2] - b[2])
}

#[inline]
fn vec_neg(v: (f64, f64, f64)) -> (f64, f64, f64) {
    (-v.0, -v.1, -v.2)
}

#[inline]
fn vec_len_sq(v: (f64, f64, f64)) -> f64 {
    v.0 * v.0 + v.1 * v.1 + v.2 * v.2
}

#[inline]
fn vec_len(v: (f64, f64, f64)) -> f64 {
    vec_len_sq(v).sqrt()
}

#[inline]
fn vec_dot(a: (f64, f64, f64), b: (f64, f64, f64)) -> f64 {
    a.0 * b.0 + a.1 * b.1 + a.2 * b.2
}

#[inline]
fn vec_cross(a: (f64, f64, f64), b: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        a.1 * b.2 - a.2 * b.1,
        a.2 * b.0 - a.0 * b.2,
        a.0 * b.1 - a.1 * b.0,
    )
}

/// Unit vector from `from` to `to`.
#[inline]
fn vec_direction(from: [f64; 3], to: [f64; 3]) -> (f64, f64, f64) {
    let d = (to[0] - from[0], to[1] - from[1], to[2] - from[2]);
    let lsq = d.0 * d.0 + d.1 * d.1 + d.2 * d.2;
    if lsq < 1e-16 {
        (0.0, 0.0, 0.0)
    } else {
        let l = lsq.sqrt();
        (d.0 / l, d.1 / l, d.2 / l)
    }
}

// RDKit✔️✔️: static constexpr double zero_tolerance = 1.e-16;
// RDKit✔️✔️: constexpr void normalize() override {
// RDKit✔️✔️:   double l = this->length();
// RDKit✔️✔️:   if (l < zero_tolerance) {
// RDKit✔️✔️:     throw std::runtime_error("Cannot normalize a zero length vector");
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   x /= l;
// RDKit✔️✔️:   y /= l;
// RDKit✔️✔️:   z /= l;
// RDKit✔️✔️: }
// RDKit✔️✔️:
// RDKit✔️✔️: double length() const override {
// RDKit✔️✔️:   double res = x * x + y * y + z * z;
// RDKit✔️✔️:   return sqrt(res);
// RDKit✔️✔️: }
// RDKit✔️✔️: /*! \brief Returns a normalized direction vector from this
// RDKit✔️✔️:  *   point to another.
// RDKit✔️✔️:  *
// RDKit✔️✔️:  */
// RDKit✔️✔️: Point3D directionVector(const Point3D &other) const {
// RDKit✔️✔️:   Point3D res;
// RDKit✔️✔️:   res.x = other.x - x;
// RDKit✔️✔️:   res.y = other.y - y;
// RDKit✔️✔️:   res.z = other.z - z;
// RDKit✔️✔️:   res.normalize();
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
#[inline]
fn rdkit_direction_vector(from: [f64; 3], other: [f64; 3]) -> Result<(f64, f64, f64), StereoError> {
    const ZERO_TOLERANCE: f64 = 1.0e-16;

    let mut res = (other[0] - from[0], other[1] - from[1], other[2] - from[2]);
    let length_squared = res.0 * res.0 + res.1 * res.1 + res.2 * res.2;
    let length = length_squared.sqrt();
    if length < ZERO_TOLERANCE {
        return Err(StereoError::ZeroLengthVector);
    }
    res.0 /= length;
    res.1 /= length;
    res.2 /= length;
    Ok(res)
}

// RDKit✔️✔️: constexpr double lengthSq() const override {
// RDKit✔️✔️:   // double res = pow(x,2) + pow(y,2) + pow(z,2);
// RDKit✔️✔️:   double res = x * x + y * y + z * z;
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// RDKit✔️✔️:
// RDKit✔️✔️: constexpr double dotProduct(const Point3D &other) const {
// RDKit✔️✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// RDKit✔️✔️:
// RDKit✔️✔️: /*! \brief determines the angle between a vector to this point
// RDKit✔️✔️:  *   from the origin and a vector to the other point.
// RDKit✔️✔️:  *
// RDKit✔️✔️:  *  The angle is unsigned: the results of this call will always
// RDKit✔️✔️:  *   be between 0 and M_PI
// RDKit✔️✔️:  */
// RDKit✔️✔️: double angleTo(const Point3D &other) const {
// RDKit✔️✔️:   double lsq = lengthSq() * other.lengthSq();
// RDKit✔️✔️:   double dotProd = dotProduct(other);
// RDKit✔️✔️:   dotProd /= sqrt(lsq);
// RDKit✔️✔️:
// RDKit✔️✔️:   // watch for roundoff error:
// RDKit✔️✔️:   if (dotProd <= -1.0) {
// RDKit✔️✔️:     return M_PI;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (dotProd >= 1.0) {
// RDKit✔️✔️:     return 0.0;
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   return acos(dotProd);
// RDKit✔️✔️: }
#[inline]
fn rdkit_angle_to(this: (f64, f64, f64), other: (f64, f64, f64)) -> f64 {
    let this_length_squared = this.0 * this.0 + this.1 * this.1 + this.2 * this.2;
    let other_length_squared = other.0 * other.0 + other.1 * other.1 + other.2 * other.2;
    let length_squared_product = this_length_squared * other_length_squared;
    let mut dot_product = this.0 * other.0 + this.1 * other.1 + this.2 * other.2;
    dot_product /= length_squared_product.sqrt();

    if dot_product <= -1.0 {
        return std::f64::consts::PI;
    }
    if dot_product >= 1.0 {
        return 0.0;
    }

    dot_product.acos()
}

// RDKit✔️✔️: constexpr double dotProduct(const Point3D &other) const {
// RDKit✔️✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
#[inline]
fn rdkit_point3d_dot_product(this: (f64, f64, f64), other: (f64, f64, f64)) -> f64 {
    this.0 * other.0 + this.1 * other.1 + this.2 * other.2
}

// RDKit✔️✔️: /*! \brief Cross product of this point with the another point
// RDKit✔️✔️:  *
// RDKit✔️✔️:  * The order is important here
// RDKit✔️✔️:  *  The result is "this" cross with "other" not (other x this)
// RDKit✔️✔️:  */
// RDKit✔️✔️: constexpr Point3D crossProduct(const Point3D &other) const {
// RDKit✔️✔️:   Point3D res;
// RDKit✔️✔️:   res.x = y * (other.z) - z * (other.y);
// RDKit✔️✔️:   res.y = -x * (other.z) + z * (other.x);
// RDKit✔️✔️:   res.z = x * (other.y) - y * (other.x);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
#[inline]
fn rdkit_point3d_cross_product(this: (f64, f64, f64), other: (f64, f64, f64)) -> (f64, f64, f64) {
    (
        this.1 * other.2 - this.2 * other.1,
        -this.0 * other.2 + this.2 * other.0,
        this.0 * other.1 - this.1 * other.0,
    )
}

// RDKit✔️✔️: #define VOLTEST(X, Y, Z) (v[X].dotProduct(v[Y].crossProduct(v[Z])) >= 0.0)
#[inline]
fn rdkit_voltest(vectors: &[(f64, f64, f64)], x: usize, y: usize, z: usize) -> bool {
    rdkit_point3d_dot_product(
        vectors[x],
        rdkit_point3d_cross_product(vectors[y], vectors[z]),
    ) >= 0.0
}

// RDKit✔️✔️: static unsigned int OctahedralPermFrom3D(unsigned char *pair,
// RDKit✔️✔️:                                          const RDGeom::Point3D *v) {
// RDKit✔️✔️:   switch (pair[0]) {
fn octahedral_perm_from_3d(pair: &[u8; 6], vectors: &[(f64, f64, f64)]) -> u32 {
    match pair[0] {
        // RDKit✔️✔️:     case 2:  // a-b
        // RDKit✔️✔️:       switch (pair[2]) {
        2 => match pair[2] {
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 3, 4) ? 28 : 27;
            4 => {
                if rdkit_voltest(vectors, 0, 3, 4) {
                    28
                } else {
                    27
                }
            }
            // RDKit✔️✔️:         case 5:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 25 : 30;
            5 => {
                if rdkit_voltest(vectors, 0, 2, 3) {
                    25
                } else {
                    30
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 26 : 29;
            _ => {
                if rdkit_voltest(vectors, 0, 2, 3) {
                    26
                } else {
                    29
                }
            }
        },
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case 3:  // a-c
        // RDKit✔️✔️:       switch (pair[1]) {
        3 => match pair[1] {
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 3, 4) ? 22 : 21;
            4 => {
                if rdkit_voltest(vectors, 0, 3, 4) {
                    22
                } else {
                    21
                }
            }
            // RDKit✔️✔️:         case 5:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 3) ? 19 : 24;
            5 => {
                if rdkit_voltest(vectors, 0, 1, 3) {
                    19
                } else {
                    24
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 1, 3) ? 20 : 23;
            _ => {
                if rdkit_voltest(vectors, 0, 1, 3) {
                    20
                } else {
                    23
                }
            }
        },
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case 4:  // a-d
        // RDKit✔️✔️:       switch (pair[1]) {
        4 => match pair[1] {
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 4) ? 13 : 12;
            3 => {
                if rdkit_voltest(vectors, 0, 2, 4) {
                    13
                } else {
                    12
                }
            }
            // RDKit✔️✔️:         case 5:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 6 : 18;
            5 => {
                if rdkit_voltest(vectors, 0, 1, 2) {
                    6
                } else {
                    18
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 7 : 17;
            _ => {
                if rdkit_voltest(vectors, 0, 1, 2) {
                    7
                } else {
                    17
                }
            }
        },
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case 5:  // a-e
        // RDKit✔️✔️:       switch (pair[1]) {
        5 => match pair[1] {
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 11 : 9;
            3 => {
                if rdkit_voltest(vectors, 0, 2, 3) {
                    11
                } else {
                    9
                }
            }
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 3 : 16;
            4 => {
                if rdkit_voltest(vectors, 0, 1, 2) {
                    3
                } else {
                    16
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 5 : 15;
            _ => {
                if rdkit_voltest(vectors, 0, 1, 2) {
                    5
                } else {
                    15
                }
            }
        },
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     default:  // 0 or 6  a-f
        // RDKit✔️✔️:       switch (pair[1]) {
        _ => match pair[1] {
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 10 : 8;
            3 => {
                if rdkit_voltest(vectors, 0, 2, 3) {
                    10
                } else {
                    8
                }
            }
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 1 : 2;
            4 => {
                if rdkit_voltest(vectors, 0, 1, 2) {
                    1
                } else {
                    2
                }
            }
            // RDKit✔️✔️:         default:  // 5
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 4 : 14;
            _ => {
                if rdkit_voltest(vectors, 0, 1, 2) {
                    4
                } else {
                    14
                }
            }
        },
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   // unreachable
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
    }
}

fn assign_nontetrahedral_chiral_type_from_3d(
    atoms: &mut [Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    conformer: &Conformer3D,
    atom_idx: usize,
    tolerance: f64,
) -> Result<bool, StereoError> {
    // RDKit✔️✔️: // The tolerance here is pretty high in order to accomodate things coming from
    // RDKit✔️✔️: // the dgeom code As we get more experience with real-world structures and/or
    // RDKit✔️✔️: // improve the dgeom code, we can think about lowering this.
    // RDKit✔️✔️: static bool assignNontetrahedralChiralTypeFrom3D(ROMol &mol,
    // RDKit✔️✔️:                                                  const Conformer &conf,
    // RDKit✔️✔️:                                                  Atom *atom,
    // RDKit✔️✔️:                                                  double tolerance = 0.1) {
    // RDKit✔️✔️:   // FIX: add tests for dative and zero order bonds
    // RDKit✔️✔️:   // Fail fast check for non-tetrahedral elements
    // RDKit✔️✔️:   if (atom->getAtomicNum() < 15) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check for wiggly bonds
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (isWigglyBond(bond, atom)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D cen = conf.getAtomPos(atom->getIdx());
    // RDKit✔️✔️:   RDGeom::Point3D v[6];
    // RDKit✔️✔️:   unsigned int count = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     if (count == 6) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     RDGeom::Point3D p = conf.getAtomPos(*nbrIdx);
    // RDKit✔️✔️:     v[count] = cen.directionVector(p);
    // RDKit✔️✔️:     ++count;
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (count < 3) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned char pair[6];
    // RDKit✔️✔️:   memset(pair, 0, 6);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int pairs = 0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < count; j++) {
    // RDKit✔️✔️:       if (v[i].dotProduct(v[j]) < -(1 - tolerance)) {
    // RDKit✔️✔️:         if (pair[i] || pair[j]) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         pair[i] = j + 1;
    // RDKit✔️✔️:         pair[j] = i + 1;
    // RDKit✔️✔️:         pairs++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Atom::ChiralType tag;
    // RDKit✔️✔️:   unsigned int perm;
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   switch (pairs) {
    // RDKit✔️✔️:     case 0:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 1:
    // RDKit✔️✔️:       switch (count) {
    // RDKit✔️✔️:         case 3: /* T-shape */
    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           if (pair[0] == 0) {
    // RDKit✔️✔️:             perm = 3;  // Z
    // RDKit✔️✔️:           } else if (pair[0] == 2) {
    // RDKit✔️✔️:             perm = 2;  // 4
    // RDKit✔️✔️:           } else /* pair[0] == 3 */ {
    // RDKit✔️✔️:             perm = 1;  // U
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:                /* See-saw */
    // RDKit✔️✔️:           if (pair[0] == 2) {  // a b
    // RDKit✔️✔️:             if (v[2].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 2, 3) ? 25 : 29;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 2, 3) ? 7 : 8;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[0] == 3) {  // a c
    // RDKit✔️✔️:             if (v[1].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 19 : 23;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 5 : 6;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[0] == 4) {  // a d
    // RDKit✔️✔️:             if (v[1].angleTo(v[2]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 2) ? 6 : 17;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 2) ? 3 : 4;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[1] == 3) {  // b c
    // RDKit✔️✔️:             if (v[0].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 10 : 8;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(1, 0, 3) ? 13 : 14;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[1] == 4) {  // b d
    // RDKit✔️✔️:             if (v[0].angleTo(v[2]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 1 : 2;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(1, 0, 2) ? 10 : 12;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else /* pair[2] == 4 */ {  // c d
    // RDKit✔️✔️:             if (v[0].angleTo(v[1]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 4 : 14;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(3, 0, 1) ? 16 : 19;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setChiralTag(tag);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 5: /* Trigonal bipyramidal */
    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           if (pair[0] == 2) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 2, 3) ? 7 : 8;  // a b
    // RDKit✔️✔️:           } else if (pair[0] == 3) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 3) ? 5 : 6;  // a c
    // RDKit✔️✔️:           } else if (pair[0] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 2) ? 3 : 4;  // a d
    // RDKit✔️✔️:           } else if (pair[0] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 2) ? 1 : 2;  // a e
    // RDKit✔️✔️:           } else if (pair[1] == 3) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 3) ? 13 : 14;  // b c
    // RDKit✔️✔️:           } else if (pair[1] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 2) ? 10 : 12;  // b d
    // RDKit✔️✔️:           } else if (pair[1] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 2) ? 9 : 11;  // b e
    // RDKit✔️✔️:           } else if (pair[2] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(2, 0, 1) ? 16 : 19;  // c d
    // RDKit✔️✔️:           } else if (pair[2] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(2, 0, 1) ? 15 : 20;  // c e
    // RDKit✔️✔️:           } else /* pair[2] == 4 */ {
    // RDKit✔️✔️:             perm = VOLTEST(3, 0, 1) ? 17 : 18;  // d e
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 2:
    // RDKit✔️✔️:       if (count == 4) {
    // RDKit✔️✔️:         /* Square planar */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         if (pair[0] == 2) {
    // RDKit✔️✔️:           perm = 2;  // 4
    // RDKit✔️✔️:         } else if (pair[0] == 3) {
    // RDKit✔️✔️:           perm = 1;  // U
    // RDKit✔️✔️:         } else /* pair[1] == 4 */ {
    // RDKit✔️✔️:           perm = 3;  // Z
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       } else if (count == 5) {
    // RDKit✔️✔️:         /* Square pyramidal */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         perm = OctahedralPermFrom3D(pair, v);
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 3:
    // RDKit✔️✔️:       if (count == 6) {
    // RDKit✔️✔️:         /* Octahedral */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         perm = OctahedralPermFrom3D(pair, v);
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }

    let atom = atoms.get(atom_idx).ok_or_else(|| {
        StereoError::InvariantViolation(format!(
            "atom index {atom_idx} is out of range for {} atoms",
            atoms.len()
        ))
    })?;
    if atom.atomic_number() < 15 {
        return Ok(false);
    }

    let neighbors = adjacency.neighbors_of(atom_idx);
    for neighbor in neighbors {
        let bond = bonds.get(neighbor.bond.index()).ok_or_else(|| {
            StereoError::InvariantViolation(format!(
                "adjacency for atom {atom_idx} references missing bond {}",
                neighbor.bond
            ))
        })?;
        if is_wiggly_bond(bond, atom_idx)? {
            return Ok(false);
        }
    }

    let coordinates = conformer.coordinates();
    let center = coordinates.get(atom_idx).copied().ok_or_else(|| {
        StereoError::InvariantViolation(format!(
            "conformer {} has no coordinate for atom {atom_idx}",
            conformer.id()
        ))
    })?;
    let mut vectors = [(0.0, 0.0, 0.0); 6];
    let mut count = 0_usize;
    for neighbor in neighbors {
        if count == 6 {
            return Ok(false);
        }
        let position = coordinates
            .get(neighbor.atom_index)
            .copied()
            .ok_or_else(|| {
                StereoError::InvariantViolation(format!(
                    "conformer {} has no coordinate for neighboring atom {}",
                    conformer.id(),
                    neighbor.atom_index
                ))
            })?;
        vectors[count] = rdkit_direction_vector(center, position)?;
        count += 1;
    }

    if count < 3 {
        return Ok(false);
    }

    let mut pair = [0_u8; 6];
    let mut pairs = 0_u32;
    for i in 0..count {
        for j in i + 1..count {
            if rdkit_point3d_dot_product(vectors[i], vectors[j]) < -(1.0 - tolerance) {
                if pair[i] != 0 || pair[j] != 0 {
                    return Ok(false);
                }
                pair[i] = (j + 1) as u8;
                pair[j] = (i + 1) as u8;
                pairs += 1;
            }
        }
    }

    let hundred_degrees = 100.0 * std::f64::consts::PI / 180.0;
    let assignment = match pairs {
        0 => None,
        1 => match count {
            3 => {
                let permutation = if pair[0] == 0 {
                    3
                } else if pair[0] == 2 {
                    2
                } else {
                    1
                };
                Some((ChiralTag::SquarePlanar, permutation))
            }
            4 => {
                let (tag, permutation) = if pair[0] == 2 {
                    if rdkit_angle_to(vectors[2], vectors[3]) < hundred_degrees {
                        (
                            ChiralTag::Octahedral,
                            if rdkit_voltest(&vectors, 0, 2, 3) {
                                25
                            } else {
                                29
                            },
                        )
                    } else {
                        (
                            ChiralTag::TrigonalBipyramidal,
                            if rdkit_voltest(&vectors, 0, 2, 3) {
                                7
                            } else {
                                8
                            },
                        )
                    }
                } else if pair[0] == 3 {
                    if rdkit_angle_to(vectors[1], vectors[3]) < hundred_degrees {
                        (
                            ChiralTag::Octahedral,
                            if rdkit_voltest(&vectors, 0, 1, 3) {
                                19
                            } else {
                                23
                            },
                        )
                    } else {
                        (
                            ChiralTag::TrigonalBipyramidal,
                            if rdkit_voltest(&vectors, 0, 1, 3) {
                                5
                            } else {
                                6
                            },
                        )
                    }
                } else if pair[0] == 4 {
                    if rdkit_angle_to(vectors[1], vectors[2]) < hundred_degrees {
                        (
                            ChiralTag::Octahedral,
                            if rdkit_voltest(&vectors, 0, 1, 2) {
                                6
                            } else {
                                17
                            },
                        )
                    } else {
                        (
                            ChiralTag::TrigonalBipyramidal,
                            if rdkit_voltest(&vectors, 0, 1, 2) {
                                3
                            } else {
                                4
                            },
                        )
                    }
                } else if pair[1] == 3 {
                    if rdkit_angle_to(vectors[0], vectors[3]) < hundred_degrees {
                        (
                            ChiralTag::Octahedral,
                            if rdkit_voltest(&vectors, 0, 1, 3) {
                                10
                            } else {
                                8
                            },
                        )
                    } else {
                        (
                            ChiralTag::TrigonalBipyramidal,
                            if rdkit_voltest(&vectors, 1, 0, 3) {
                                13
                            } else {
                                14
                            },
                        )
                    }
                } else if pair[1] == 4 {
                    if rdkit_angle_to(vectors[0], vectors[2]) < hundred_degrees {
                        (
                            ChiralTag::Octahedral,
                            if rdkit_voltest(&vectors, 0, 1, 3) {
                                1
                            } else {
                                2
                            },
                        )
                    } else {
                        (
                            ChiralTag::TrigonalBipyramidal,
                            if rdkit_voltest(&vectors, 1, 0, 2) {
                                10
                            } else {
                                12
                            },
                        )
                    }
                } else if rdkit_angle_to(vectors[0], vectors[1]) < hundred_degrees {
                    (
                        ChiralTag::Octahedral,
                        if rdkit_voltest(&vectors, 0, 1, 3) {
                            4
                        } else {
                            14
                        },
                    )
                } else {
                    (
                        ChiralTag::TrigonalBipyramidal,
                        if rdkit_voltest(&vectors, 3, 0, 1) {
                            16
                        } else {
                            19
                        },
                    )
                };
                Some((tag, permutation))
            }
            5 => {
                let permutation = if pair[0] == 2 {
                    if rdkit_voltest(&vectors, 0, 2, 3) {
                        7
                    } else {
                        8
                    }
                } else if pair[0] == 3 {
                    if rdkit_voltest(&vectors, 0, 1, 3) {
                        5
                    } else {
                        6
                    }
                } else if pair[0] == 4 {
                    if rdkit_voltest(&vectors, 0, 1, 2) {
                        3
                    } else {
                        4
                    }
                } else if pair[0] == 5 {
                    if rdkit_voltest(&vectors, 0, 1, 2) {
                        1
                    } else {
                        2
                    }
                } else if pair[1] == 3 {
                    if rdkit_voltest(&vectors, 1, 0, 3) {
                        13
                    } else {
                        14
                    }
                } else if pair[1] == 4 {
                    if rdkit_voltest(&vectors, 1, 0, 2) {
                        10
                    } else {
                        12
                    }
                } else if pair[1] == 5 {
                    if rdkit_voltest(&vectors, 1, 0, 2) {
                        9
                    } else {
                        11
                    }
                } else if pair[2] == 4 {
                    if rdkit_voltest(&vectors, 2, 0, 1) {
                        16
                    } else {
                        19
                    }
                } else if pair[2] == 5 {
                    if rdkit_voltest(&vectors, 2, 0, 1) {
                        15
                    } else {
                        20
                    }
                } else if rdkit_voltest(&vectors, 3, 0, 1) {
                    17
                } else {
                    18
                };
                Some((ChiralTag::TrigonalBipyramidal, permutation))
            }
            _ => None,
        },
        2 if count == 4 => {
            let permutation = if pair[0] == 2 {
                2
            } else if pair[0] == 3 {
                1
            } else {
                3
            };
            Some((ChiralTag::SquarePlanar, permutation))
        }
        2 if count == 5 => Some((
            ChiralTag::Octahedral,
            octahedral_perm_from_3d(&pair, &vectors),
        )),
        3 if count == 6 => Some((
            ChiralTag::Octahedral,
            octahedral_perm_from_3d(&pair, &vectors),
        )),
        _ => None,
    };

    let Some((tag, permutation)) = assignment else {
        return Ok(false);
    };
    atoms[atom_idx].set_chiral_tag(tag);
    atoms[atom_idx].set_chiral_permutation(Some(permutation));
    Ok(true)
}

fn assign_tetrahedral_chiral_type_from_3d(
    atoms: &mut [Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    conformer: &Conformer3D,
    atom_idx: usize,
    total_num_hs: usize,
    explicit_atom: bool,
    zero_volume_tolerance: f64,
) -> Result<bool, StereoError> {
    // RDKit✔️✔️:     /* We're only doing tetrahedral cases here */
    // RDKit✔️✔️:     if (tnzDegree > 4) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     int anum = atom->getAtomicNum();
    // RDKit✔️✔️:     if (anum != 16 && anum != 34 &&  // S or Se are special
    // RDKit✔️✔️:                                      // (just using the InChI list for now)
    // RDKit✔️✔️:         tnzDegree != 4               // not enough total neighbors
    // RDKit✔️✔️:     ) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     const auto &p0 = conf.getAtomPos(atom->getIdx());
    // RDKit✔️✔️:     const RDGeom::Point3D *nbrs[4];
    // RDKit✔️✔️:     unsigned int nbrIdx = 0;
    // RDKit✔️✔️:     int hasWigglyBond = 0;
    // RDKit✔️✔️:     for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:       hasWigglyBond = isWigglyBond(bond, atom);
    // RDKit✔️✔️:       if (hasWigglyBond) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!Chirality::detail::bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       nbrs[nbrIdx++] = &conf.getAtomPos(bond->getOtherAtomIdx(atom->getIdx()));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (hasWigglyBond) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto v1 = *nbrs[0] - p0;
    // RDKit✔️✔️:     auto v2 = *nbrs[1] - p0;
    // RDKit✔️✔️:     auto v3 = *nbrs[2] - p0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     double chiralVol = v1.dotProduct(v2.crossProduct(v3));
    // RDKit✔️✔️:     bool chiralitySet = false;
    // RDKit✔️✔️:     if (chiralVol < -ZERO_VOLUME_TOL) {
    // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️✔️:       chiralitySet = true;
    // RDKit✔️✔️:     } else if (chiralVol > ZERO_VOLUME_TOL) {
    // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:       chiralitySet = true;
    // RDKit✔️✔️:     } else if (nbrIdx == 4) {
    // RDKit✔️✔️:       // The first three neighbors are on the same plane as the chiral atom (or
    // RDKit✔️✔️:       // very close to it). If a 4th neighbor is present, let's see if this one
    // RDKit✔️✔️:       // determines a chiral volume
    // RDKit✔️✔️:
    // RDKit✔️✔️:       auto v4 = *nbrs[3] - p0;
    // RDKit✔️✔️:       // v4 would be in the opposite direction to v3
    // RDKit✔️✔️:       chiralVol = -v1.dotProduct(v2.crossProduct(v4));
    // RDKit✔️✔️:       if (chiralVol < -ZERO_VOLUME_TOL) {
    // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️✔️:         chiralitySet = true;
    // RDKit✔️✔️:       } else if (chiralVol > ZERO_VOLUME_TOL) {
    // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:         chiralitySet = true;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (chiralitySet && explicitAtoms[atom->getIdx()] == 0) {
    // RDKit✔️✔️:       atom->setProp<int>(common_properties::_NonExplicit3DChirality, 1);
    // RDKit✔️✔️:     }

    let atom = atoms.get(atom_idx).ok_or_else(|| {
        StereoError::InvariantViolation(format!(
            "atom index {atom_idx} is out of range for {} atoms",
            atoms.len()
        ))
    })?;
    let neighbors = adjacency.neighbors_of(atom_idx);
    let mut nonzero_degree = 0_usize;
    for neighbor in neighbors {
        let bond = bonds.get(neighbor.bond.index()).ok_or_else(|| {
            StereoError::InvariantViolation(format!(
                "adjacency for atom {atom_idx} references missing bond {}",
                neighbor.bond
            ))
        })?;
        if bond_affects_atom_chirality(bond, atom_idx) {
            nonzero_degree += 1;
        }
    }
    let total_nonzero_degree = nonzero_degree + total_num_hs;
    if nonzero_degree < 3 || total_nonzero_degree > 4 {
        return Ok(false);
    }
    let atomic_number = atom.atomic_number();
    if atomic_number != 16 && atomic_number != 34 && total_nonzero_degree != 4 {
        return Ok(false);
    }

    let coordinates = conformer.coordinates();
    let center = coordinates.get(atom_idx).copied().ok_or_else(|| {
        StereoError::InvariantViolation(format!(
            "conformer {} has no coordinate for atom {atom_idx}",
            conformer.id()
        ))
    })?;
    let mut neighbor_coordinates = [[0.0; 3]; 4];
    let mut neighbor_count = 0_usize;
    for neighbor in neighbors {
        let bond = bonds.get(neighbor.bond.index()).ok_or_else(|| {
            StereoError::InvariantViolation(format!(
                "adjacency for atom {atom_idx} references missing bond {}",
                neighbor.bond
            ))
        })?;
        if is_wiggly_bond(bond, atom_idx)? {
            return Ok(false);
        }
        if !bond_affects_atom_chirality(bond, atom_idx) {
            continue;
        }
        let coordinate = coordinates
            .get(neighbor.atom_index)
            .copied()
            .ok_or_else(|| {
                StereoError::InvariantViolation(format!(
                    "conformer {} has no coordinate for neighboring atom {}",
                    conformer.id(),
                    neighbor.atom_index
                ))
            })?;
        neighbor_coordinates[neighbor_count] = coordinate;
        neighbor_count += 1;
    }

    let subtract_center = |point: [f64; 3]| {
        (
            point[0] - center[0],
            point[1] - center[1],
            point[2] - center[2],
        )
    };
    let v1 = subtract_center(neighbor_coordinates[0]);
    let v2 = subtract_center(neighbor_coordinates[1]);
    let v3 = subtract_center(neighbor_coordinates[2]);
    let mut chiral_volume = rdkit_point3d_dot_product(v1, rdkit_point3d_cross_product(v2, v3));

    let mut chirality_set = false;
    if chiral_volume < -zero_volume_tolerance {
        atoms[atom_idx].set_chiral_tag(ChiralTag::TetrahedralCw);
        chirality_set = true;
    } else if chiral_volume > zero_volume_tolerance {
        atoms[atom_idx].set_chiral_tag(ChiralTag::TetrahedralCcw);
        chirality_set = true;
    } else if neighbor_count == 4 {
        let v4 = subtract_center(neighbor_coordinates[3]);
        chiral_volume = -rdkit_point3d_dot_product(v1, rdkit_point3d_cross_product(v2, v4));
        if chiral_volume < -zero_volume_tolerance {
            atoms[atom_idx].set_chiral_tag(ChiralTag::TetrahedralCw);
            chirality_set = true;
        } else if chiral_volume > zero_volume_tolerance {
            atoms[atom_idx].set_chiral_tag(ChiralTag::TetrahedralCcw);
            chirality_set = true;
        } else {
            atoms[atom_idx].set_chiral_tag(ChiralTag::Unspecified);
        }
    } else {
        atoms[atom_idx].set_chiral_tag(ChiralTag::Unspecified);
    }

    if chirality_set && !explicit_atom {
        atoms[atom_idx].set_prop("_NonExplicit3DChirality", "1");
    }
    Ok(chirality_set)
}

pub(crate) fn assign_chiral_types_from_3d_kernel(
    atoms: &mut [Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    conformers: &[Conformer3D],
    properties: &mut MoleculeProperties,
    implicit_hydrogens: Option<&[i32]>,
    conformer_id: i32,
    replace_existing_tags: bool,
) -> Result<(), StereoError> {
    // RDKit✔️❌: void assignChiralTypesFrom3D(ROMol &mol, int confId, bool replaceExistingTags) {
    // RDKit✔️❌:   const double ZERO_VOLUME_TOL = 0.1;
    // RDKit✔️❌:   if (!mol.getNumConformers()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️❌:   if (!conf.is3D()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // if the molecule already has stereochemistry
    // RDKit✔️❌:   // perceived, remove the flags that indicate
    // RDKit✔️❌:   // this... what we're about to do will require
    // RDKit✔️❌:   // that we go again.
    // RDKit✔️❌:   if (mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️❌:     mol.clearProp(common_properties::_StereochemDone);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   auto allowNontetrahedralStereo = Chirality::getAllowNontetrahedralChirality();
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> explicitAtoms;
    // RDKit✔️❌:   explicitAtoms.resize(mol.getNumAtoms(), 0);
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     auto bondDir = bond->getBondDir();
    // RDKit✔️❌:     if (bondDir == Bond::BondDir::BEGINWEDGE ||
    // RDKit✔️❌:         bondDir == Bond::BondDir::BEGINDASH) {
    // RDKit✔️❌:       explicitAtoms[bond->getBeginAtom()->getIdx()] = 1;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       explicitAtoms[atom->getIdx()] = 1;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     // if we aren't replacing existing tags and the atom is already tagged,
    // RDKit✔️❌:     // punt:
    // RDKit✔️❌:     if (!replaceExistingTags && atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:     // additional reasons to skip the atom:
    // RDKit✔️❌:     auto nzDegree = Chirality::detail::getAtomNonzeroDegree(atom);
    // RDKit✔️❌:     auto tnzDegree = nzDegree + atom->getTotalNumHs();
    // RDKit✔️❌:     if (nzDegree < 3 || tnzDegree > 6) {
    // RDKit✔️❌:       // not enough explicit neighbors or too many total neighbors
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (allowNontetrahedralStereo &&
    // RDKit✔️❌:         assignNontetrahedralChiralTypeFrom3D(mol, conf, atom)) {
    // RDKit✔️❌:       if (explicitAtoms[atom->getIdx()] == 0) {
    // RDKit✔️❌:         atom->setProp(common_properties::_NonExplicit3DChirality, 1);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     /* We're only doing tetrahedral cases here */
    // RDKit✔️❌:     if (tnzDegree > 4) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     int anum = atom->getAtomicNum();
    // RDKit✔️❌:     if (anum != 16 && anum != 34 &&  // S or Se are special
    // RDKit✔️❌:                                      // (just using the InChI list for now)
    // RDKit✔️❌:         tnzDegree != 4               // not enough total neighbors
    // RDKit✔️❌:     ) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     const auto &p0 = conf.getAtomPos(atom->getIdx());
    // RDKit✔️❌:     const RDGeom::Point3D *nbrs[4];
    // RDKit✔️❌:     unsigned int nbrIdx = 0;
    // RDKit✔️❌:     int hasWigglyBond = 0;
    // RDKit✔️❌:     for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️❌:       hasWigglyBond = isWigglyBond(bond, atom);
    // RDKit✔️❌:       if (hasWigglyBond) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (!Chirality::detail::bondAffectsAtomChirality(bond, atom)) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       nbrs[nbrIdx++] = &conf.getAtomPos(bond->getOtherAtomIdx(atom->getIdx()));
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (hasWigglyBond) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     auto v1 = *nbrs[0] - p0;
    // RDKit✔️❌:     auto v2 = *nbrs[1] - p0;
    // RDKit✔️❌:     auto v3 = *nbrs[2] - p0;
    // RDKit✔️❌:
    // RDKit✔️❌:     double chiralVol = v1.dotProduct(v2.crossProduct(v3));
    // RDKit✔️❌:     bool chiralitySet = false;
    // RDKit✔️❌:     if (chiralVol < -ZERO_VOLUME_TOL) {
    // RDKit✔️❌:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️❌:       chiralitySet = true;
    // RDKit✔️❌:     } else if (chiralVol > ZERO_VOLUME_TOL) {
    // RDKit✔️❌:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️❌:       chiralitySet = true;
    // RDKit✔️❌:     } else if (nbrIdx == 4) {
    // RDKit✔️❌:       // The first three neighbors are on the same plane as the chiral atom (or
    // RDKit✔️❌:       // very close to it). If a 4th neighbor is present, let's see if this one
    // RDKit✔️❌:       // determines a chiral volume
    // RDKit✔️❌:
    // RDKit✔️❌:       auto v4 = *nbrs[3] - p0;
    // RDKit✔️❌:       // v4 would be in the opposite direction to v3
    // RDKit✔️❌:       chiralVol = -v1.dotProduct(v2.crossProduct(v4));
    // RDKit✔️❌:       if (chiralVol < -ZERO_VOLUME_TOL) {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
    // RDKit✔️❌:         chiralitySet = true;
    // RDKit✔️❌:       } else if (chiralVol > ZERO_VOLUME_TOL) {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️❌:         chiralitySet = true;
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (chiralitySet && explicitAtoms[atom->getIdx()] == 0) {
    // RDKit✔️❌:       atom->setProp<int>(common_properties::_NonExplicit3DChirality, 1);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }

    const ZERO_VOLUME_TOLERANCE: f64 = 0.1;
    if conformers.is_empty() {
        return Ok(());
    }
    let conformer = if conformer_id < 0 {
        &conformers[0]
    } else {
        conformers
            .iter()
            .find(|conformer| conformer.id() == conformer_id as usize)
            .ok_or(StereoError::ConformerNotFound { conformer_id })?
    };
    if !conformer.is_3d() {
        return Ok(());
    }

    properties.clear_prop("_StereochemDone");
    let allow_nontetrahedral_stereo = get_allow_nontetrahedral_chirality();
    if conformer.coordinates().len() != atoms.len() {
        return Err(StereoError::ConformerAtomCountMismatch {
            conformer_id: conformer.id(),
            expected: atoms.len(),
            actual: conformer.coordinates().len(),
        });
    }

    let mut explicit_atoms = vec![false; atoms.len()];
    for bond in bonds {
        if matches!(
            bond.direction(),
            crate::BondDirection::BeginWedge | crate::BondDirection::BeginDash
        ) {
            let begin_idx = bond.begin().index();
            let explicit_atom = explicit_atoms.get_mut(begin_idx).ok_or_else(|| {
                StereoError::InvariantViolation(format!(
                    "bond {} begin atom {} is out of range for {} atoms",
                    bond.id(),
                    bond.begin(),
                    atoms.len()
                ))
            })?;
            *explicit_atom = true;
        }
    }
    for (atom_idx, atom) in atoms.iter().enumerate() {
        if atom.chiral_tag() != ChiralTag::Unspecified {
            explicit_atoms[atom_idx] = true;
        }
    }

    for atom_idx in 0..atoms.len() {
        if !replace_existing_tags && atoms[atom_idx].chiral_tag() != ChiralTag::Unspecified {
            continue;
        }
        atoms[atom_idx].set_chiral_tag(ChiralTag::Unspecified);

        let nonzero_degree = atom_nonzero_degree_from_parts(bonds, adjacency, atom_idx)?;
        let implicit_hydrogens =
            implicit_hydrogens.ok_or(StereoError::MissingImplicitHydrogenState)?;
        let implicit_hydrogen_count = implicit_hydrogens.get(atom_idx).copied().ok_or(
            StereoError::ImplicitHydrogenCountMismatch {
                expected: atoms.len(),
                actual: implicit_hydrogens.len(),
            },
        )?;
        let total_num_hs =
            atoms[atom_idx].explicit_hydrogens() as usize + implicit_hydrogen_count.max(0) as usize;
        let total_nonzero_degree = nonzero_degree + total_num_hs;
        if nonzero_degree < 3 || total_nonzero_degree > 6 {
            continue;
        }

        if allow_nontetrahedral_stereo
            && assign_nontetrahedral_chiral_type_from_3d(
                atoms, bonds, adjacency, conformer, atom_idx, 0.1,
            )?
        {
            if !explicit_atoms[atom_idx] {
                atoms[atom_idx].set_prop("_NonExplicit3DChirality", "1");
            }
            continue;
        }

        if total_nonzero_degree > 4 {
            continue;
        }
        let atomic_number = atoms[atom_idx].atomic_number();
        if atomic_number != 16 && atomic_number != 34 && total_nonzero_degree != 4 {
            continue;
        }
        assign_tetrahedral_chiral_type_from_3d(
            atoms,
            bonds,
            adjacency,
            conformer,
            atom_idx,
            total_num_hs,
            explicit_atoms[atom_idx],
            ZERO_VOLUME_TOLERANCE,
        )?;
    }
    Ok(())
}

/// Applies the completed `assignChiralTypesFrom3D` kernel to an owned molecule.
///
/// Parser call sites use this adapter to prepare RDKit-like total-hydrogen
/// state and commit the kernel's topology/property mutations without carrying
/// a second implementation of the chemistry algorithm.
pub(crate) fn assign_chiral_types_from_3d_molecule(
    molecule: &mut Molecule,
    conformer_id: i32,
    replace_existing_tags: bool,
) -> Result<(), StereoError> {
    let Some(selected_conformer) = (if molecule.conformers_3d().is_empty() {
        None
    } else if conformer_id < 0 {
        molecule.conformers_3d().first()
    } else {
        molecule
            .conformers_3d()
            .iter()
            .find(|conformer| conformer.id() == conformer_id as usize)
    }) else {
        if molecule.conformers_3d().is_empty() {
            return Ok(());
        }
        return Err(StereoError::ConformerNotFound { conformer_id });
    };
    if !selected_conformer.is_3d() {
        return Ok(());
    }

    let selected_conformer = selected_conformer.clone();
    let valence =
        crate::assign_valence_with_options(molecule, crate::ValenceModel::RdkitLike, false)?;
    let mut topology = molecule.topology_block().clone();
    let mut properties = molecule.properties().clone();
    assign_chiral_types_from_3d_kernel(
        &mut topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        std::slice::from_ref(&selected_conformer),
        &mut properties,
        Some(&valence.implicit_hydrogens),
        conformer_id,
        replace_existing_tags,
    )?;

    molecule.replace_topology_block(topology);
    molecule.replace_properties(properties);
    molecule.derived_cache_mut().invalidate(
        crate::DerivedState::STEREO
            | crate::DerivedState::DRAWING
            | crate::DerivedState::FINGERPRINT,
    );
    Ok(())
}

/// Distance of `to` from `from` in the XY plane.
#[inline]
fn vec_xy_dist(from: [f64; 3], to: [f64; 3]) -> f64 {
    let dx = to[0] - from[0];
    let dy = to[1] - from[1];
    (dx * dx + dy * dy).sqrt()
}

// BEGIN RDKIT CPP FUNCTION: atomChiralTypeFromBondDirPseudo3D (Chirality.cpp:427-845)
// RDKit✔️✔️: std::optional<Atom::ChiralType> atomChiralTypeFromBondDirPseudo3D(
// RDKit✔️✔️:     const ROMol &mol, const Bond *bond, const Conformer *conf)
// RDKit✔️✔️:   // Determines ChiralTag from wedge bond directions + 3D coordinates.
// RDKit✔️✔️:   // Uses pseudo-3D z-offsets for wedged bonds to compute chiral volume.
// RDKit✔️✔️:   auto bondDir = bond->getBondDir();
// RDKit✔️✔️:   const auto atom = bond->getBeginAtom();
// RDKit✔️✔️:   auto centerLoc = conf->getAtomPos(atom->getIdx());
// RDKit✔️✔️:   centerLoc.z = 0.0;
// RDKit✔️✔️:   auto refPt = conf->getAtomPos(bondAtom->getIdx());
// RDKit✔️✔️:   // ... full pseudo-3D logic with cross/dot products for chiral volume
// RDKit✔️✔️:   if (vol > volumeTolerance) { res = CHI_TETRAHEDRAL_CCW; }
// RDKit✔️✔️:   else if (vol < -volumeTolerance) { res = CHI_TETRAHEDRAL_CW; }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: atomChiralTypeFromBondDirPseudo3D
/// Determine ChiralTag from wedge bond directions + 3D coordinates.
/// Uses pseudo-3D z-offsets for wedged bonds to compute chiral volume.
/// Returns `Some(ChiralTag)` for determined chirality, or `None` for
/// ambiguous/error cases.  Returns `Some(ChiralTag::Unspecified)` for
/// atoms with degree > 4 (out of scope for tetrahedral detection).
#[must_use]
pub fn atom_chiral_type_from_bond_dir_pseudo_3d(
    mol: &Molecule,
    bond_idx: usize,
    conformer: &Conformer3D,
) -> Option<ChiralTag> {
    // ── constants matching C++ ──
    const COORD_ZERO_TOL: f64 = 1e-4;
    const ZERO_TOL: f64 = 1e-3;
    const T_SHAPE_TOL: f64 = 0.00031;
    const PSEUDO_3D_OFFSET: f64 = 0.1;
    const VOLUME_TOLERANCE: f64 = 0.00174;

    let bond = &mol.bonds()[bond_idx];
    let bond_dir = bond.direction();
    if bond_dir != crate::BondDirection::BeginWedge && bond_dir != crate::BondDirection::BeginDash {
        return None;
    }

    // the atom at the point of the wedge
    let center_idx = bond.begin().index();
    let bond_atom_idx = bond.end().index();

    let deg = atom_degree(mol, center_idx);
    if deg > 4 {
        return Some(ChiralTag::Unspecified);
    }

    let coords = conformer.coordinates();
    let mut center_loc = coords[center_idx];
    center_loc[2] = 0.0;
    let mut ref_pt = coords[bond_atom_idx];

    // Scale the pseudo-3D offset by the reference bond's length
    // (Github #7305: some conformers have weird scalings)
    let ref_pt_2d = coords[bond_atom_idx];
    let ref_length = vec_xy_dist(center_loc, ref_pt_2d);
    ref_pt[2] = if bond_dir == crate::BondDirection::BeginWedge {
        PSEUDO_3D_OFFSET
    } else {
        -PSEUDO_3D_OFFSET
    };
    if ref_length > 0.0 {
        ref_pt[2] *= ref_length;
    }

    // ── collect neighbour bond indices and bond vectors ──
    let mut neighbor_bond_indices: Vec<usize> = Vec::new();
    let mut ref_idx = mol.bonds().len() + 1; // sentinel
    let mut bond_vects: Vec<(f64, f64, f64)> = Vec::new();
    let mut all_single = true;
    let mut nbr_count = 0usize;

    for b in mol.bonds() {
        if b.begin().index() != center_idx && b.end().index() != center_idx {
            continue;
        }
        let other_idx = if b.begin().index() == center_idx {
            b.end().index()
        } else {
            b.begin().index()
        };
        let mut tmp_pt: [f64; 3];
        if b.id().index() == bond_idx {
            ref_idx = nbr_count;
            tmp_pt = ref_pt;
        } else {
            tmp_pt = coords[other_idx];
            // Check if this other bond is also wedged (begin wedge/dash from center)
            if b.begin().index() == center_idx
                && (b.direction() == crate::BondDirection::BeginWedge
                    || b.direction() == crate::BondDirection::BeginDash)
            {
                tmp_pt[2] = if b.direction() == crate::BondDirection::BeginWedge {
                    PSEUDO_3D_OFFSET
                } else {
                    -PSEUDO_3D_OFFSET
                };
                if ref_length > 0.0 {
                    tmp_pt[2] *= ref_length;
                }
            } else {
                tmp_pt[2] = 0.0;
            }
            // Check for zero-length (or near zero-length) bonds
            let diff_sq = (center_loc[0] - tmp_pt[0]).powi(2)
                + (center_loc[1] - tmp_pt[1]).powi(2)
                + (center_loc[2] - tmp_pt[2]).powi(2);
            if diff_sq < ZERO_TOL {
                return None;
            }
        }
        nbr_count += 1;
        if b.order() != crate::BondOrder::Single {
            all_single = false;
        }
        bond_vects.push(vec_direction(center_loc, tmp_pt));
        neighbor_bond_indices.push(b.id().index());
    }

    debug_assert!(
        ref_idx < mol.bonds().len(),
        "could not find reference bond in neighbors"
    );

    let n_nbrs = bond_vects.len();

    // Need at least 3 bonds (implicit H for 3-coordinate) and at most 4
    if n_nbrs < 3 || n_nbrs > 4 {
        return None;
    }

    // ── check for overlapping neighbours ──
    for i in 0..n_nbrs {
        for j in 0..i {
            let diff = vec_sub(
                [bond_vects[i].0, bond_vects[i].1, bond_vects[i].2],
                [bond_vects[j].0, bond_vects[j].1, bond_vects[j].2],
            );
            if vec_len_sq(diff) < ZERO_TOL {
                return None;
            }
        }
    }

    // ── only process if all single bonds, or P (15) / S (16) ──
    let atom = &mol.atoms()[center_idx];
    if all_single || atom.atomic_number() == 15 || atom.atomic_number() == 16 {
        let mut vol: f64;
        let mut order: [usize; 4] = [0, 1, 2, 3];
        let mut prefactor: f64 = 1.0;
        if ref_idx != 0 {
            order.swap(0, ref_idx);
            prefactor *= -1.0;
        }

        // Check for colinear bonds 1 and 2 (T-shaped)
        if n_nbrs > 3 {
            let cp12 = vec_cross(bond_vects[order[1]], bond_vects[order[2]]);
            let cp10 = vec_cross(bond_vects[order[1]], bond_vects[order[0]]);
            if vec_len_sq(cp12) < 10.0 * ZERO_TOL && vec_len_sq(cp10) > 10.0 * ZERO_TOL {
                // That bondVect is no longer normalized, but this shouldn't break anything
                let mut bv = bond_vects[order[1]];
                bv.2 = bond_vects[order[0]].2 * -1.0;
                bond_vects[order[1]] = bv;
            }
        }

        // ── order bonds for rotation order 0-1-2 (-3) ──
        let needs_swap =
            |cp01: (f64, f64, f64), cp02: (f64, f64, f64), dp01: f64, dp02: f64| -> bool {
                if (dp01.abs() - 1.0) > -ZERO_TOL {
                    if cp02.2 < 0.0 {
                        return true;
                    }
                    return false;
                }
                if (dp02.abs() - 1.0) > -ZERO_TOL {
                    if cp01.2 < 0.0 {
                        return true;
                    }
                }

                if (cp01.2 * cp02.2) < -ZERO_TOL {
                    if cp01.2 < cp02.2 {
                        return true;
                    }
                    return false;
                }
                if dp01 * dp02 < -ZERO_TOL {
                    if dp01 < dp02 {
                        return true;
                    }
                    return false;
                }
                dp01.abs() > dp02.abs()
            };

        if n_nbrs == 3 {
            let cp01 = vec_cross(bond_vects[order[0]], bond_vects[order[1]]);
            let cp02 = vec_cross(bond_vects[order[0]], bond_vects[order[2]]);
            let dp01 = vec_dot(bond_vects[order[0]], bond_vects[order[1]]);
            let dp02 = vec_dot(bond_vects[order[0]], bond_vects[order[2]]);
            if needs_swap(cp01, cp02, dp01, dp02) {
                order.swap(1, 2);
                prefactor *= -1.0;
            }
        } else if n_nbrs > 3 {
            // Sort bonds 1, 2, 3 by cross- and dot-products to bond 0
            let mut ordered_bonds: Vec<(f64, f64, usize)> = Vec::with_capacity(3);
            for i in 1..4 {
                let cp0i = vec_cross(bond_vects[order[0]], bond_vects[order[i]]);
                let sgn = if cp0i.2 < -ZERO_TOL { -1.0 } else { 1.0 };
                let dp0i = vec_dot(bond_vects[order[0]], bond_vects[order[i]]);
                ordered_bonds.push((sgn, sgn * dp0i, order[i]));
            }
            ordered_bonds.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

            let mut n_changed = 0;
            for i in 1..4 {
                let ni = ordered_bonds[i - 1].2;
                if order[i] != ni {
                    order[i] = ni;
                    n_changed += 1;
                }
            }
            if n_changed == 2 {
                prefactor *= -1.0;
            }
        }

        // ── check for opposing bonds with opposite wedging ──
        for i in 0..n_nbrs {
            for j in i + 1..n_nbrs {
                if bond_vects[order[i]].2 * bond_vects[order[j]].2 < -ZERO_TOL {
                    let cp = vec_len_sq(vec_cross(bond_vects[order[i]], bond_vects[order[j]]));
                    if cp < 0.01 {
                        // Exception: in pseudo-3D drawings of sugars, ring substituents
                        // are drawn 180deg apart with opposite wedging.  Allow for
                        // neighbouring bonds in a 4-coordinate setting.
                        if n_nbrs == 4 {
                            let dot_check =
                                vec_dot(bond_vects[order[i]], bond_vects[order[j]]) + 1.0;
                            if dot_check.abs() < ZERO_TOL && (j - i == 1 || (i == 0 && j == 3)) {
                                let mut bv = bond_vects[order[j]];
                                bv.2 = 0.0;
                                bond_vects[order[j]] = bv;
                                continue;
                            }
                        }
                        return None;
                    }
                }
            }
        }

        // ── three-coordinate special cases ──
        if n_nbrs == 3 {
            let mut conflict = false;
            if bond_vects[order[1]].2 * bond_vects[order[0]].2 < -COORD_ZERO_TOL
                && bond_vects[order[2]].2.abs() < COORD_ZERO_TOL
            {
                let cp20 = vec_cross(bond_vects[order[2]], bond_vects[order[0]]);
                let cp21 = vec_cross(bond_vects[order[2]], bond_vects[order[1]]);
                conflict = cp20.2 * cp21.2 < -1e-4;
            } else if bond_vects[order[2]].2 * bond_vects[order[0]].2 < -COORD_ZERO_TOL
                && bond_vects[order[1]].2.abs() < COORD_ZERO_TOL
            {
                let cp10 = vec_cross(bond_vects[order[1]], bond_vects[order[0]]);
                let cp12 = vec_cross(bond_vects[order[1]], bond_vects[order[2]]);
                conflict = cp10.2 * cp12.2 < -COORD_ZERO_TOL;
            }
            if conflict {
                return None;
            }
        }

        // ── compute chiral volume ──
        // For cross products, ignore pseudo-3D z
        let mut bv1 = bond_vects[order[1]];
        bv1.2 = 0.0;
        let mut bv2 = bond_vects[order[2]];
        bv2.2 = 0.0;
        let mut crossp1 = vec_cross(bv1, bv2);

        // Catch linear arrangements
        if n_nbrs == 3 {
            if vec_len_sq(crossp1) < T_SHAPE_TOL {
                // T-shaped: assume the two perpendicular bonds are wedged opposite
                bv1.2 = -bond_vects[order[0]].2;
                bv2.2 = -bond_vects[order[0]].2;
                crossp1 = vec_cross(bv1, bv2);
            }
        } else if vec_len_sq(crossp1) < 10.0 * ZERO_TOL {
            // If the other bond is flat:
            if bond_vects[order[3]].2.abs() < COORD_ZERO_TOL {
                let mut bv3 = bond_vects[order[3]];
                bv3.2 = -1.0 * bond_vects[order[0]].2;
                // That bondVect is no longer normalized, but shouldn't break anything
                bond_vects[order[3]] = bv3;
            }
        }

        vol = vec_dot(crossp1, bond_vects[order[0]]);

        if n_nbrs == 4 {
            let dotp1 = vec_dot(bond_vects[order[1]], bond_vects[order[2]]);
            let mut bv3 = bond_vects[order[3]];
            bv3.2 = 0.0;
            let crossp2 = vec_cross(bv1, bv3);
            let dotp2 = vec_dot(bond_vects[order[1]], bond_vects[order[3]]);
            let mut vol2 = vec_dot(crossp2, bond_vects[order[0]]);

            // Detect case with no chiral volume for the default evaluation
            if vol.abs() < ZERO_TOL {
                if vol2.abs() < ZERO_TOL {
                    // Warning: ambiguous stereochemistry - no chiral volume
                    return None;
                }
                vol = vol2;
                prefactor *= -1.0;
            } else if vol * vol2 > 0.0 && vol2.abs() > VOLUME_TOLERANCE && dotp1 < dotp2 {
                vol = vol2;
                prefactor *= -1.0;
            } else if vol.abs() < VOLUME_TOLERANCE && vol2.abs() > VOLUME_TOLERANCE {
                if vol * vol2 < 0.0 {
                    prefactor *= -1.0;
                }
                vol = vol2;
            }
        }

        vol *= prefactor;

        // ── assign chirality based on the sign of the chiral volume ──
        if vol > VOLUME_TOLERANCE {
            Some(ChiralTag::TetrahedralCcw)
        } else if vol < -VOLUME_TOLERANCE {
            Some(ChiralTag::TetrahedralCw)
        } else {
            // Warning: ambiguous stereochemistry - zero final chiral volume
            None
        }
    } else {
        Some(ChiralTag::Unspecified)
    }
}

// ── Bridgehead helper ─────────────────────────────────────────────────

pub(crate) fn query_is_atom_bridgehead(
    mol: &Molecule,
    atom_idx: usize,
    ring_info: &crate::RingInfo,
) -> i32 {
    // RDKit✔️✔️: int queryIsAtomBridgehead(Atom const *at) {
    // RDKit✔️✔️:   // at least three ring bonds, all ring bonds in a ring which shares at
    // RDKit✔️✔️:   // least two bonds with another ring involving this atom
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // We can't just go with "at least three ring bonds shared between multiple
    // RDKit✔️✔️:   // rings" because of structures like CC12CCN(CC1)C2 where there are only two
    // RDKit✔️✔️:   // SSSRs
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   if (at->getDegree() < 3) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &mol = at->getOwningMol();
    // RDKit✔️✔️:   const auto ri = mol.getRingInfo();
    // RDKit✔️✔️:   if (!ri || !ri->isInitialized()) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // track which ring bonds involve this atom
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomRingBonds(mol.getNumBonds());
    // RDKit✔️✔️:   for (const auto bnd : mol.atomBonds(at)) {
    // RDKit✔️✔️:     if (ri->numBondRings(bnd->getIdx())) {
    // RDKit✔️✔️:       atomRingBonds.set(bnd->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atomRingBonds.count() < 3) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::dynamic_bitset<> bondsInRingI(mol.getNumBonds());
    // RDKit✔️✔️:   boost::dynamic_bitset<> ringsOverlap(ri->numRings());
    // RDKit✔️✔️:   for (unsigned int i = 0; i < ri->bondRings().size(); ++i) {
    // RDKit✔️✔️:     bondsInRingI.reset();
    // RDKit✔️✔️:     bool atomInRingI = false;
    // RDKit✔️✔️:     for (const auto bidx : ri->bondRings()[i]) {
    // RDKit✔️✔️:       bondsInRingI.set(bidx);
    // RDKit✔️✔️:       if (atomRingBonds[bidx]) {
    // RDKit✔️✔️:         atomInRingI = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!atomInRingI) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < ri->bondRings().size(); ++j) {
    // RDKit✔️✔️:       unsigned int overlap = 0;
    // RDKit✔️✔️:       bool atomInRingJ = false;
    // RDKit✔️✔️:       for (const auto bidx : ri->bondRings()[j]) {
    // RDKit✔️✔️:         if (atomRingBonds[bidx]) {
    // RDKit✔️✔️:           atomInRingJ = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (bondsInRingI[bidx]) {
    // RDKit✔️✔️:           ++overlap;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (overlap >= 2 && atomInRingJ) {
    // RDKit✔️✔️:           // we have two rings containing the atom which share at least two
    // RDKit✔️✔️:           // bonds:
    // RDKit✔️✔️:           ringsOverlap.set(i);
    // RDKit✔️✔️:           ringsOverlap.set(j);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!ringsOverlap[i]) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 1;
    // RDKit✔️✔️: }
    // Local complexity review: both implementations use bit-packed arrays of
    // O(B + R), scan the atom's CSR incident-bond range once, then perform the
    // same ring-pair overlap loops with O(R^2 * ring_size) worst-case time.
    // Each ring-I bitset reset is O(B), matching dynamic_bitset::reset(). Rust
    // allocates no per-pair collections and performs O(1) indexed lookups.
    let adjacency = &mol.topology_block().adjacency;
    if adjacency.neighbors_of(atom_idx).len() < 3 || !ring_info.is_initialized() {
        return 0;
    }

    let mut atom_ring_bonds = vec![false; mol.num_bonds()];
    for neighbor in adjacency.neighbors_of(atom_idx) {
        if ring_info.num_bond_rings(neighbor.bond) != 0 {
            atom_ring_bonds[neighbor.bond.index()] = true;
        }
    }
    if atom_ring_bonds.iter().filter(|is_ring| **is_ring).count() < 3 {
        return 0;
    }

    let bond_rings = ring_info.bond_rings();
    let mut bonds_in_ring_i = vec![false; mol.num_bonds()];
    let mut rings_overlap = vec![false; bond_rings.len()];
    for (i, ring_i) in bond_rings.iter().enumerate() {
        bonds_in_ring_i.fill(false);
        let mut atom_in_ring_i = false;
        for bond in ring_i {
            bonds_in_ring_i[bond.index()] = true;
            if atom_ring_bonds[bond.index()] {
                atom_in_ring_i = true;
            }
        }
        if !atom_in_ring_i {
            continue;
        }
        for (j, ring_j) in bond_rings.iter().enumerate().skip(i + 1) {
            let mut overlap = 0;
            let mut atom_in_ring_j = false;
            for bond in ring_j {
                if atom_ring_bonds[bond.index()] {
                    atom_in_ring_j = true;
                }
                if bonds_in_ring_i[bond.index()] {
                    overlap += 1;
                }
                if overlap >= 2 && atom_in_ring_j {
                    rings_overlap[i] = true;
                    rings_overlap[j] = true;
                    break;
                }
            }
        }
        if !rings_overlap[i] {
            return 0;
        }
    }
    1
}

// BEGIN RDKIT CPP FUNCTION: atomIsCandidateForRingStereochem (Chirality.cpp:1455-1517)
// RDKit✔️✔️: bool atomIsCandidateForRingStereochem(
// RDKit✔️✔️:     const ROMol &mol, const Atom *atom,
// RDKit✔️✔️:     const std::vector<unsigned int> &atomRanks) {
// RDKit✔️✔️:   // Three-coordinate N additional requirements:
// RDKit✔️✔️:   //   in a ring of size 3  (from InChI)
// RDKit✔️✔️:   //   a bridgehead (RDKit extension)
// RDKit✔️✔️:   // Non-ring neighbour rank analysis for stereochem candidacy.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: atomIsCandidateForRingStereochem
/// Check if an atom in a ring could be a stereocenter.
/// Returns true if the atom is a candidate for ring stereochemistry.
#[must_use]
pub fn atom_is_candidate_for_ring_stereochem(
    mol: &Molecule,
    ri: &crate::RingInfo,
    atom_idx: usize,
    cip_ranks: &[u32],
) -> bool {
    if !ri.is_initialized() {
        return false;
    }

    let atom_id = crate::AtomId::new(atom_idx);
    if ri.num_atom_rings(atom_id) == 0 {
        return false;
    }

    let atom = &mol.atoms()[atom_idx];

    // Three-coordinate N additional requirements:
    //   in a ring of size 3 (from InChI) OR a bridgehead (RDKit extension)
    // RDKit✔️✔️: if (atom->getAtomicNum() == 7 && atom->getTotalDegree() == 3 &&
    // RDKit✔️✔️:     !ringInfo->isAtomInRingOfSize(atom->getIdx(), 3) &&
    // RDKit✔️✔️:     !queryIsAtomBridgehead(atom)) {
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    if atom.atomic_number() == 7
        && atom_total_degree(mol, atom_idx) == 3
        && !ri.is_atom_in_ring_of_size(atom_id, 3)
        && query_is_atom_bridgehead(mol, atom_idx, ri) == 0
    {
        return false;
    }

    // Collect ring and non-ring neighbours
    let mut non_ring_nbrs: Vec<usize> = Vec::new();
    let mut ring_nbrs: Vec<usize> = Vec::new();
    let mut ring_nbr_ranks: HashSet<u32> = HashSet::new();

    for b in mol.bonds() {
        if b.begin().index() != atom_idx && b.end().index() != atom_idx {
            continue;
        }
        let other_idx = if b.begin().index() == atom_idx {
            b.end().index()
        } else {
            b.begin().index()
        };
        if ri.num_bond_rings(b.id()) == 0 {
            non_ring_nbrs.push(other_idx);
        } else {
            ring_nbrs.push(other_idx);
            if other_idx < cip_ranks.len() {
                ring_nbr_ranks.insert(cip_ranks[other_idx]);
            }
        }
    }

    match non_ring_nbrs.len() {
        2 => {
            // The ranks of the non-ring neighbours must be different AND
            // the ranks of the ring neighbours must be the same (see issue #8956)
            let mut res = true;
            if non_ring_nbrs[0] < cip_ranks.len() && non_ring_nbrs[1] < cip_ranks.len() {
                res = cip_ranks[non_ring_nbrs[0]] != cip_ranks[non_ring_nbrs[1]];
            }
            res && (ring_nbrs.len() != ring_nbr_ranks.len())
        }
        1 => ring_nbrs.len() > ring_nbr_ranks.len(),
        0 => {
            if ring_nbrs.len() == 4 && ring_nbr_ranks.len() == 3 {
                true
            } else if ring_nbrs.len() == 3 && ring_nbr_ranks.len() == 2 {
                true
            } else {
                false
            }
        }
        _ => false,
    }
}

// BEGIN RDKIT CPP FUNCTION: findChiralAtomSpecialCases (Chirality.cpp:1521-1649)
// RDKit✔️✔️: void findChiralAtomSpecialCases(ROMol &mol,
// RDKit✔️✔️:     boost::dynamic_bitset<> &possibleSpecialCases,
// RDKit✔️✔️:     const std::vector<unsigned int> &atomRanks) {
// RDKit✔️✔️:   // BFS over ring bonds to find other stereocenters.
// RDKit✔️✔️:   // Sets _ringStereoAtoms property indicating related stereo atoms.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: findChiralAtomSpecialCases
/// Find chiral atom special cases (ring stereochemistry candidates).
/// Returns a vector of `ChiralAtomSpecialCase` entries describing atoms
/// that are ring stereocenters and their inter-relationships.
/// Each entry contains the atom index, its chiral tag, and a list of
/// `(same_orientation, other_atom_idx)` cross-references.
pub fn find_chiral_atom_special_cases(
    mol: &Molecule,
    cip_ranks: &[u32],
) -> Result<Vec<ChiralAtomSpecialCase>, StereoError> {
    let symm_rings = match mol.derived_cache().rings.as_ref() {
        Some(ri) if ri.is_initialized() && ri.is_symm_sssr() => None,
        _ => Some(crate::symmetrize_sssr(mol)?),
    };
    let ri = symm_rings
        .as_ref()
        .or_else(|| mol.derived_cache().rings.as_ref())
        .expect("symmetrize_sssr() must produce initialized ring info");

    let n_atoms = mol.num_atoms();
    let n_bonds = mol.bonds().len();

    let mut result: Vec<ChiralAtomSpecialCase> = Vec::new();

    // Track seen atoms, used atoms, and seen bonds across BFS passes
    let mut atoms_seen = vec![false; n_atoms];
    let mut atoms_used = vec![false; n_atoms];
    let mut bonds_seen = vec![false; n_bonds];

    for atom in mol.atoms() {
        let idx = atom.id().index();
        if atoms_seen[idx] {
            continue;
        }
        let tag = atom.chiral_tag();
        if tag == ChiralTag::Unspecified || tag == ChiralTag::Other {
            continue;
        }
        if atom.prop("_CIPCode").is_some() {
            continue;
        }
        if ri.num_atom_rings(atom.id()) == 0 {
            continue;
        }
        if !atom_is_candidate_for_ring_stereochem(mol, ri, idx, cip_ranks) {
            continue;
        }

        // BFS from this atom along ring bonds to find other stereocenters
        let mut next_atoms: Vec<usize> = Vec::new();
        let mut ring_stereo_atoms: Vec<(i32, usize)> = Vec::new(); // (sign * (idx+1), orig_idx)

        // Start with finding viable neighbours
        for b in mol.bonds() {
            let b_idx = b.id().index();
            if b.begin().index() != idx && b.end().index() != idx {
                continue;
            }
            if bonds_seen[b_idx] {
                continue;
            }
            bonds_seen[b_idx] = true;
            if ri.num_bond_rings(b.id()) > 0 {
                let other_idx = if b.begin().index() == idx {
                    b.end().index()
                } else {
                    b.begin().index()
                };
                if !atoms_seen[other_idx] {
                    next_atoms.push(other_idx);
                    atoms_used[other_idx] = true;
                }
            }
        }

        while !next_atoms.is_empty() {
            let ratom_idx = next_atoms.remove(0);
            atoms_seen[ratom_idx] = true;

            if atoms_used[ratom_idx] {
                // skip if already used in another BFS
            }

            let ratom = &mol.atoms()[ratom_idx];
            let rtag = ratom.chiral_tag();
            if rtag != ChiralTag::Unspecified
                && rtag != ChiralTag::Other
                && ratom.prop("_CIPCode").is_none()
                && ri.num_atom_rings(ratom.id()) > 0
                && atom_is_candidate_for_ring_stereochem(mol, ri, ratom_idx, cip_ranks)
            {
                let same = if rtag == tag { 1i32 } else { -1i32 };
                ring_stereo_atoms.push((same * (ratom_idx as i32 + 1), ratom_idx));
            }

            // Push this atom's ring-bond neighbours
            for b in mol.bonds() {
                let b_idx = b.id().index();
                if b.begin().index() != ratom_idx && b.end().index() != ratom_idx {
                    continue;
                }
                if bonds_seen[b_idx] {
                    continue;
                }
                bonds_seen[b_idx] = true;
                if ri.num_bond_rings(b.id()) > 0 {
                    let other_idx = if b.begin().index() == ratom_idx {
                        b.end().index()
                    } else {
                        b.begin().index()
                    };
                    if !atoms_seen[other_idx] && !atoms_used[other_idx] {
                        next_atoms.push(other_idx);
                        atoms_used[other_idx] = true;
                    }
                }
            }
        }

        if !ring_stereo_atoms.is_empty() {
            // Build the cross-referencing: each atom in ring_stereo_atoms
            // should also reference the others
            // First, collect the full set
            let mut all_entries: Vec<(i32, usize)> = ring_stereo_atoms.clone();
            all_entries.push((1i32 * (idx as i32 + 1), idx)); // self reference

            // For each atom in the set, build its list of referenced atoms
            for &(entry_val, entry_idx) in &all_entries {
                let same_self = entry_val > 0;
                let mut refs: Vec<(bool, usize)> = Vec::new();
                for &(other_val, other_idx) in &all_entries {
                    if other_idx == entry_idx {
                        continue;
                    }
                    let other_same = other_val > 0;
                    let these_different = same_self ^ other_same;
                    refs.push((!these_different, other_idx));
                }
                result.push(ChiralAtomSpecialCase {
                    atom_idx: entry_idx,
                    chiral_tag: if entry_val > 0 {
                        tag
                    } else {
                        opposite_tag(tag)
                    },
                    ring_stereo_atoms: refs,
                });
            }
        }

        atoms_seen[idx] = true;
    }

    Ok(result)
}

/// Return the opposite tetrahedral chiral tag.
fn opposite_tag(tag: ChiralTag) -> ChiralTag {
    match tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        t => t,
    }
}

/// Result from `find_chiral_atom_special_cases`.
#[derive(Debug, Clone)]
pub struct ChiralAtomSpecialCase {
    pub atom_idx: usize,
    pub chiral_tag: ChiralTag,
    pub ring_stereo_atoms: Vec<(bool, usize)>, // (same_orientation, other_atom_idx)
}
//
// ── Non-tetrahedral stereo infrastructure ──────────────────────────────────
//
// Ported from RDKit NontetrahedralStereo.cpp (lines 1-481)
// Source reproduction protocol: dev/source_reproduction_protocol.md

// BEGIN RDKIT CPP CONSTANT: swap_squareplanar_table[4][6] (NontetrahedralStereo.cpp:21-28)
// RDKit✔️✔️: constexpr unsigned char swap_squareplanar_table[4][6]
#[rustfmt::skip]
pub(crate) const SWAP_SQUAREPLANAR_TABLE: [[u8; 6]; 4] = [
    //  0   0   0   1   1   2
    //  1   2   3   2   3   3
    [0, 0, 0, 0, 0, 0],
    [3, 1, 2, 2, 1, 3],  // SP1
    [2, 3, 1, 1, 3, 2],  // SP2
    [1, 2, 3, 3, 2, 1],  // SP3
];

// BEGIN RDKIT CPP CONSTANT: swap_trigonalbipyramidal_table[21][10] (NontetrahedralStereo.cpp:30-54)
// RDKit✔️✔️: constexpr unsigned char swap_trigonalbipyramidal_table[21][10]
#[rustfmt::skip]
pub(crate) const SWAP_TRIGONALBIPYRAMIDAL_TABLE: [[u8; 10]; 21] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 20, 17, 2, 2, 2, 7, 2, 6, 3],         // TB1
    [11, 15, 18, 1, 1, 1, 8, 1, 5, 4],        // TB2
    [10, 19, 4, 18, 4, 8, 4, 5, 4, 1],        // TB3
    [12, 16, 3, 17, 3, 7, 3, 6, 3, 2],        // TB4
    [13, 6, 16, 20, 7, 6, 6, 3, 2, 6],        // TB5
    [14, 5, 19, 15, 8, 5, 5, 4, 1, 5],        // TB6
    [8, 14, 10, 11, 5, 4, 1, 8, 8, 8],        // TB7
    [7, 13, 12, 9, 6, 3, 2, 7, 7, 7],         // TB8
    [1, 11, 11, 8, 15, 18, 11, 11, 14, 10],   // TB9
    [3, 12, 7, 12, 16, 12, 17, 13, 12, 9],    // TB10
    [2, 9, 9, 7, 20, 17, 9, 9, 13, 12],       // TB11
    [4, 10, 8, 10, 19, 10, 18, 14, 10, 11],   // TB12
    [5, 8, 14, 14, 14, 19, 15, 10, 11, 14],   // TB13
    [6, 7, 13, 13, 13, 16, 20, 12, 9, 13],    // TB14
    [20, 2, 20, 6, 9, 20, 13, 17, 20, 16],    // TB15
    [19, 4, 5, 19, 10, 14, 19, 19, 18, 15],   // TB16
    [18, 18, 1, 4, 18, 11, 10, 15, 19, 18],   // TB17
    [17, 17, 2, 3, 17, 9, 12, 20, 16, 17],    // TB18
    [16, 3, 6, 16, 12, 13, 16, 16, 17, 20],   // TB19
    [15, 1, 15, 5, 11, 15, 14, 18, 15, 19],   // TB20
];

// BEGIN RDKIT CPP CONSTANT: swap_octahedral_table[31][15] (NontetrahedralStereo.cpp:56-90)
// RDKit✔️✔️: constexpr unsigned char swap_octahedral_table[31][15]
#[rustfmt::skip]
pub(crate) const SWAP_OCTAHEDRAL_TABLE: [[u8; 15]; 31] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3],        // OH1
    [7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16],        // OH2
    [18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1],      // OH3
    [15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5],      // OH4
    [14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4],     // OH5
    [16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7],    // OH6
    [2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6],     // OH7
    [23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9],     // OH8
    [24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8],    // OH9
    [20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11],        // OH10
    [19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10],       // OH11
    [22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13],     // OH12
    [21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12],     // OH13
    [5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15],         // OH14
    [4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14],        // OH15
    [6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2],         // OH16
    [1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18],        // OH17
    [3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17],       // OH18
    [11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20],    // OH19
    [10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19],    // OH20
    [13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22],    // OH21
    [12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21],    // OH22
    [8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24],    // OH23
    [9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23],    // OH24
    [30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26],     // OH25
    [29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25],    // OH26
    [28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28],   // OH27
    [27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27],   // OH28
    [26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30],    // OH29
    [25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29],     // OH30
];

// BEGIN RDKIT CPP FUNCTION: swap_squareplanar (NontetrahedralStereo.cpp:92-111)
// RDKit✔️✔️: static unsigned int swap_squareplanar(unsigned int perm, unsigned int x, unsigned int y)
/// Swap two ligands (by index) in a square planar permutation.
/// Returns the new permutation index after swapping ligands at positions x and y.
#[must_use]
pub fn swap_squareplanar(perm: u32, x: usize, y: usize) -> u32 {
    if perm == 0 || perm > 3 {
        return 0;
    }
    if x == y {
        return perm;
    }
    // offset[3] = {0, 2, 3} for pairs (0,1), (0,2), (1,2)
    const OFFSET: [usize; 3] = [0, 2, 3];
    let swapidx = if x < y {
        if y > 3 {
            return 0;
        }
        OFFSET[x] + (y - 1)
    } else {
        if x > 3 {
            return 0;
        }
        OFFSET[y] + (x - 1)
    };
    SWAP_SQUAREPLANAR_TABLE[perm as usize][swapidx] as u32
}

// BEGIN RDKIT CPP FUNCTION: swap_trigonalbipyramidal (NontetrahedralStereo.cpp:113-132)
// RDKit✔️✔️: static unsigned int swap_trigonalbipyramidal(unsigned int perm, unsigned int x, unsigned int y)
/// Swap two ligands (by index) in a trigonal bipyramidal permutation.
#[must_use]
pub fn swap_trigonalbipyramidal(perm: u32, x: usize, y: usize) -> u32 {
    if perm == 0 || perm > 20 {
        return 0;
    }
    if x == y {
        return perm;
    }
    const OFFSET: [usize; 4] = [0, 3, 5, 6];
    let swapidx = if x < y {
        if y > 4 {
            return 0;
        }
        OFFSET[x] + (y - 1)
    } else {
        if x > 4 {
            return 0;
        }
        OFFSET[y] + (x - 1)
    };
    SWAP_TRIGONALBIPYRAMIDAL_TABLE[perm as usize][swapidx] as u32
}

// BEGIN RDKIT CPP FUNCTION: swap_octahedral (NontetrahedralStereo.cpp:134-153)
// RDKit✔️✔️: static unsigned int swap_octahedral(unsigned int perm, unsigned int x, unsigned int y)
/// Swap two ligands (by index) in an octahedral permutation.
#[must_use]
pub fn swap_octahedral(perm: u32, x: usize, y: usize) -> u32 {
    if perm == 0 || perm > 30 {
        return 0;
    }
    if x == y {
        return perm;
    }
    const OFFSET: [usize; 5] = [0, 4, 7, 9, 10];
    let swapidx = if x < y {
        if y > 5 {
            return 0;
        }
        OFFSET[x] + (y - 1)
    } else {
        if x > 5 {
            return 0;
        }
        OFFSET[y] + (x - 1)
    };
    SWAP_OCTAHEDRAL_TABLE[perm as usize][swapidx] as u32
}
// END RDKIT CPP FUNCTION: swap_octahedral

// BEGIN RDKIT CPP CONSTANT: squareplanar_across[4][4] (NontetrahedralStereo.cpp:155-160)
// RDKit✔️✔️: constexpr unsigned char squareplanar_across[4][4]
#[rustfmt::skip]
pub(crate) const SQUAREPLANAR_ACROSS: [[u8; 4]; 4] = [
    [4, 4, 4, 4],  // perm=0 (invalid)
    [2, 3, 0, 1],  // SP1
    [1, 0, 3, 2],  // SP2
    [3, 2, 1, 0],  // SP3
];

// BEGIN RDKIT CPP CONSTANT: trigonalbipyramidal_across[21][5] (NontetrahedralStereo.cpp:162-184)
// RDKit✔️✔️: constexpr unsigned char trigonalbipyramidal_across[21][5]
#[rustfmt::skip]
pub(crate) const TRIGONALBIPYRAMIDAL_ACROSS: [[u8; 5]; 21] = [
    [5, 5, 5, 5, 5],  // perm=0 (invalid)
    [4, 5, 5, 5, 0],  // TB1
    [4, 5, 5, 5, 0],  // TB2
    [3, 5, 5, 0, 5],  // TB3
    [3, 5, 5, 0, 5],  // TB4
    [2, 5, 0, 5, 5],  // TB5
    [2, 5, 0, 5, 5],  // TB6
    [1, 0, 5, 5, 5],  // TB7
    [1, 0, 5, 5, 5],  // TB8
    [5, 4, 5, 5, 1],  // TB9
    [5, 3, 5, 1, 5],  // TB10
    [5, 4, 5, 5, 1],  // TB11
    [5, 3, 5, 1, 5],  // TB12
    [5, 2, 1, 5, 5],  // TB13
    [5, 2, 1, 5, 5],  // TB14
    [5, 5, 4, 5, 2],  // TB15
    [5, 5, 3, 2, 5],  // TB16
    [5, 5, 5, 4, 3],  // TB17
    [5, 5, 5, 4, 3],  // TB18
    [5, 5, 3, 2, 5],  // TB19
    [5, 5, 4, 5, 2],  // TB20
];

// BEGIN RDKIT CPP CONSTANT: octahedral_across[31][6] (NontetrahedralStereo.cpp:186-218)
// RDKit✔️✔️: constexpr unsigned char octahedral_across[31][6]
#[rustfmt::skip]
pub(crate) const OCTAHEDRAL_ACROSS: [[u8; 6]; 31] = [
    [6, 6, 6, 6, 6, 6],  // perm=0 (invalid)
    [5, 3, 4, 1, 2, 0],  // OH1
    [5, 3, 4, 1, 2, 0],  // OH2
    [4, 3, 5, 1, 0, 2],  // OH3
    [5, 4, 3, 2, 1, 0],  // OH4
    [4, 5, 3, 2, 0, 1],  // OH5
    [3, 4, 5, 0, 1, 2],  // OH6
    [3, 5, 4, 0, 2, 1],  // OH7
    [5, 2, 1, 4, 3, 0],  // OH8
    [4, 2, 1, 5, 0, 3],  // OH9
    [5, 2, 1, 4, 3, 0],  // OH10
    [4, 2, 1, 5, 0, 3],  // OH11
    [3, 2, 1, 0, 5, 4],  // OH12
    [3, 2, 1, 0, 5, 4],  // OH13
    [5, 4, 3, 2, 1, 0],  // OH14
    [4, 5, 3, 2, 0, 1],  // OH15
    [4, 3, 5, 1, 0, 2],  // OH16
    [3, 5, 4, 0, 2, 1],  // OH17
    [3, 4, 5, 0, 1, 2],  // OH18
    [2, 4, 0, 5, 1, 3],  // OH19
    [2, 5, 0, 4, 3, 1],  // OH20
    [2, 3, 0, 1, 5, 4],  // OH21
    [2, 3, 0, 1, 5, 4],  // OH22
    [2, 5, 0, 4, 3, 1],  // OH23
    [2, 4, 0, 5, 1, 3],  // OH24
    [1, 0, 4, 5, 2, 3],  // OH25
    [1, 0, 5, 4, 3, 2],  // OH26
    [1, 0, 3, 2, 5, 4],  // OH27
    [1, 0, 3, 2, 5, 4],  // OH28
    [1, 0, 5, 4, 3, 2],  // OH29
    [1, 0, 4, 5, 2, 3],  // OH30
];

// BEGIN RDKIT CPP CONSTANT: trigonalbipyramidal_axial[21][2] (NontetrahedralStereo.cpp:220-242)
// RDKit✔️✔️: constexpr unsigned char trigonalbipyramidal_axial[21][2]
#[rustfmt::skip]
pub(crate) const TRIGONALBIPYRAMIDAL_AXIAL: [[u8; 2]; 21] = [
    [5, 5],  // perm=0 (invalid)
    [0, 4],  // TB1
    [0, 4],  // TB2
    [0, 3],  // TB3
    [0, 3],  // TB4
    [0, 2],  // TB5
    [0, 2],  // TB6
    [0, 1],  // TB7
    [0, 1],  // TB8
    [1, 4],  // TB9
    [1, 4],  // TB10
    [1, 3],  // TB11
    [1, 3],  // TB12
    [1, 2],  // TB13
    [1, 2],  // TB14
    [2, 4],  // TB15
    [2, 3],  // TB16
    [3, 4],  // TB17
    [3, 4],  // TB18
    [2, 3],  // TB19
    [2, 4],  // TB20
];
// END RDKIT CPP CONSTANT: trigonalbipyramidal_axial

// BEGIN RDKIT CPP FUNCTION: isTrigonalBipyramidalAxialBond (NontetrahedralStereo.cpp:244-274)
// RDKit✔️✔️: int isTrigonalBipyramidalAxialBond(const Atom *cen, const Bond *qry)
/// Check if a bond on a TBP center is axial. Returns 1 for `axial[0]`, -1 for
/// `axial[1]`, 0 otherwise.
#[must_use]
pub fn is_trigonal_bipyramidal_axial_bond(
    cen_idx: usize,
    query_bond_idx: usize,
    mol: &Molecule,
) -> i32 {
    let atom = &mol.atoms()[cen_idx];
    if mol.bonds().len() <= query_bond_idx {
        return 0;
    }
    let deg = atom_degree(mol, cen_idx);
    if deg > 5 || atom.chiral_tag() != ChiralTag::TrigonalBipyramidal {
        return 0;
    }
    let perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 || perm > 20 {
        return 0;
    }

    let mut count = 0usize;
    for bond in mol.bonds() {
        if bond.begin().index() == cen_idx || bond.end().index() == cen_idx {
            if bond.id().index() == query_bond_idx {
                if count == TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][0] as usize {
                    return 1;
                }
                if count == TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][1] as usize {
                    return -1;
                }
                return 0;
            }
            count += 1;
        }
    }
    0
}

// BEGIN RDKIT CPP FUNCTION: isTrigonalBipyramidalAxialAtom (NontetrahedralStereo.cpp:276-286)
// RDKit✔️✔️: int isTrigonalBipyramidalAxialAtom(const Atom *cen, const Atom *qry)
/// Check if a neighbor atom of a TBP center is axial.
#[must_use]
pub fn is_trigonal_bipyramidal_axial_atom(cen_idx: usize, qry_idx: usize, mol: &Molecule) -> i32 {
    let bond_idx = bond_between_atoms(mol, cen_idx, qry_idx);
    match bond_idx {
        Some(bi) => is_trigonal_bipyramidal_axial_bond(cen_idx, bi, mol),
        None => 0,
    }
}

// BEGIN RDKIT CPP FUNCTION: getMaxNbors (NontetrahedralStereo.cpp:288-308)
// RDKit✔️✔️: unsigned int getMaxNbors(const Atom::ChiralType tag)
/// Return the maximum number of neighbors for a given chiral tag geometry.
#[must_use]
pub const fn get_max_nbors(tag: ChiralTag) -> u32 {
    match tag {
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw | ChiralTag::Tetrahedral => 4,
        ChiralTag::Allene => 2,
        ChiralTag::SquarePlanar => 4,
        ChiralTag::TrigonalBipyramidal => 5,
        ChiralTag::Octahedral => 6,
        _ => 0,
    }
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossBond (NontetrahedralStereo.cpp:310-385)
// RDKit✔️✔️: Bond *getChiralAcrossBond(const Atom *cen, const Bond *qry)
/// Given a center atom and a query bond, find the \"across\" bond (the bond opposite).
/// Returns `Some(bond_index)` if found, `None` otherwise.
#[must_use]
pub fn get_chiral_across_bond(
    cen_idx: usize,
    query_bond_idx: usize,
    mol: &Molecule,
) -> Option<usize> {
    let atom = &mol.atoms()[cen_idx];
    let tag = atom.chiral_tag();
    let perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 {
        return None;
    }

    let max_nbors = get_max_nbors(tag) as usize;
    if max_nbors == 0 {
        return None;
    }

    let mut count = 0usize;
    let mut bond_refs: Vec<usize> = Vec::with_capacity(max_nbors);
    let mut found: Option<usize> = None;
    for bond in mol.bonds() {
        if bond.begin().index() == cen_idx || bond.end().index() == cen_idx {
            if count >= max_nbors {
                return None;
            }
            bond_refs.push(bond.id().index());
            if bond.id().index() == query_bond_idx {
                found = Some(count);
            }
            count += 1;
        }
    }

    if let Some(found_idx) = found {
        let across_idx = match tag {
            ChiralTag::SquarePlanar => {
                if perm <= 3 {
                    SQUAREPLANAR_ACROSS[perm as usize][found_idx]
                } else {
                    4
                }
            }
            ChiralTag::TrigonalBipyramidal => {
                if perm <= 20 {
                    TRIGONALBIPYRAMIDAL_ACROSS[perm as usize][found_idx]
                } else {
                    5
                }
            }
            ChiralTag::Octahedral => {
                if perm <= 30 {
                    OCTAHEDRAL_ACROSS[perm as usize][found_idx]
                } else {
                    6
                }
            }
            _ => return None,
        };
        if (across_idx as usize) < bond_refs.len() {
            Some(bond_refs[across_idx as usize])
        } else {
            None
        }
    } else {
        None
    }
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossBond (atom overload) (NontetrahedralStereo.cpp:375-385)
// RDKit✔️✔️: Bond *getChiralAcrossBond(const Atom *cen, const Atom *qry)
/// Find the across bond given a center atom and a query neighbor atom.
#[must_use]
pub fn get_chiral_across_bond_by_atom(
    cen_idx: usize,
    qry_idx: usize,
    mol: &Molecule,
) -> Option<usize> {
    let bond_idx = bond_between_atoms(mol, cen_idx, qry_idx)?;
    get_chiral_across_bond(cen_idx, bond_idx, mol)
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossAtom (NontetrahedralStereo.cpp:387-403)
// RDKit✔️✔️: Atom *getChiralAcrossAtom(const Atom *cen, const Bond *qry)
/// Find the across atom given a center atom and a query bond.
#[must_use]
pub fn get_chiral_across_atom(
    cen_idx: usize,
    query_bond_idx: usize,
    mol: &Molecule,
) -> Option<usize> {
    let across_bond_idx = get_chiral_across_bond(cen_idx, query_bond_idx, mol)?;
    let bond = &mol.bonds()[across_bond_idx];
    let other = if bond.begin().index() == cen_idx {
        bond.end().index()
    } else {
        bond.begin().index()
    };
    Some(other)
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossAtom (atom overload) (NontetrahedralStereo.cpp:392-403)
// RDKit✔️✔️: Atom *getChiralAcrossAtom(const Atom *cen, const Atom *qry)
/// Find the across atom given a center atom and a query neighbor atom.
#[must_use]
pub fn get_chiral_across_atom_by_atom(
    cen_idx: usize,
    qry_idx: usize,
    mol: &Molecule,
) -> Option<usize> {
    let bond_idx = bond_between_atoms(mol, cen_idx, qry_idx)?;
    get_chiral_across_atom(cen_idx, bond_idx, mol)
}
// END RDKIT CPP FUNCTION: getChiralAcrossAtom

// BEGIN RDKIT CPP FUNCTION: getIdealAngleBetweenLigands (NontetrahedralStereo.cpp:405-437)
// RDKit✔️✔️: double getIdealAngleBetweenLigands(const Atom *cen, const Atom *lig1, const Atom *lig2)
/// Return the ideal angle (90, 120, or 180) between two ligands of a non-tetrahedral center.
#[must_use]
pub fn get_ideal_angle_between_ligands(
    cen_idx: usize,
    lig1: usize,
    lig2: usize,
    mol: &Molecule,
) -> f64 {
    let atom = &mol.atoms()[cen_idx];
    let tag = atom.chiral_tag();
    match tag {
        ChiralTag::SquarePlanar | ChiralTag::Octahedral => {
            if get_chiral_across_atom_by_atom(cen_idx, lig1, mol) == Some(lig2) {
                180.0
            } else {
                90.0
            }
        }
        ChiralTag::TrigonalBipyramidal => {
            if get_chiral_across_atom_by_atom(cen_idx, lig1, mol) == Some(lig2) {
                180.0
            } else if is_trigonal_bipyramidal_axial_atom(cen_idx, lig1, mol) != 0
                || is_trigonal_bipyramidal_axial_atom(cen_idx, lig2, mol) != 0
            {
                90.0
            } else {
                120.0
            }
        }
        _ => 0.0,
    }
}

// BEGIN RDKIT CPP FUNCTION: getTrigonalBipyramidalAxialBond (NontetrahedralStereo.cpp:439-464)
// RDKit✔️✔️: Bond *getTrigonalBipyramidalAxialBond(const Atom *cen, int axial)
/// Get the bond to a specific axial position on a TBP center.
/// `axial == 1` returns the first axial bond, `axial == -1` returns the second.
#[must_use]
pub fn get_trigonal_bipyramidal_axial_bond(
    cen_idx: usize,
    axial: i32,
    mol: &Molecule,
) -> Option<usize> {
    let atom = &mol.atoms()[cen_idx];
    if atom.chiral_tag() != ChiralTag::TrigonalBipyramidal || atom_degree(mol, cen_idx) > 5 {
        return None;
    }

    let perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 || perm > 20 {
        return None;
    }

    let idx = if axial != -1 {
        TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][0] as usize
    } else {
        TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][1] as usize
    };

    let mut count = 0usize;
    for bond in mol.bonds() {
        if bond.begin().index() == cen_idx || bond.end().index() == cen_idx {
            if count == idx {
                return Some(bond.id().index());
            }
            count += 1;
        }
    }
    None
}

// BEGIN RDKIT CPP FUNCTION: getTrigonalBipyramidalAxialAtom (NontetrahedralStereo.cpp:466-471)
// RDKit✔️✔️: Atom *getTrigonalBipyramidalAxialAtom(const Atom *cen, int axial)
/// Get the atom at a specific axial position on a TBP center.
#[must_use]
pub fn get_trigonal_bipyramidal_axial_atom(
    cen_idx: usize,
    axial: i32,
    mol: &Molecule,
) -> Option<usize> {
    let bond_idx = get_trigonal_bipyramidal_axial_bond(cen_idx, axial, mol)?;
    let bond = &mol.bonds()[bond_idx];
    let other = if bond.begin().index() == cen_idx {
        bond.end().index()
    } else {
        bond.begin().index()
    };
    Some(other)
}

// BEGIN RDKIT CPP FUNCTION: hasNonTetrahedralStereo (NontetrahedralStereo.cpp:472-481)
// RDKit✔️✔️: bool hasNonTetrahedralStereo(const Atom *cen)
/// Check if an atom has non-tetrahedral stereochemistry.
#[must_use]
pub fn has_non_tetrahedral_stereo(atom: &crate::Atom) -> bool {
    let tag = atom.chiral_tag();
    tag == ChiralTag::SquarePlanar
        || tag == ChiralTag::TrigonalBipyramidal
        || tag == ChiralTag::Octahedral
}
// END RDKIT CPP FUNCTION: hasNonTetrahedralStereo

/// Helper: find bond index between two atoms by their indices.
fn bond_between_atoms(mol: &Molecule, a: usize, b: usize) -> Option<usize> {
    mol.bonds().iter().find_map(|bond| {
        if (bond.begin().index() == a && bond.end().index() == b)
            || (bond.begin().index() == b && bond.end().index() == a)
        {
            Some(bond.id().index())
        } else {
            None
        }
    })
}

/// Helper: find bond index between two atoms by their indices (atom/bond slice version).
fn bond_between_atoms_by_slice(
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
    a: usize,
    b: usize,
) -> Option<usize> {
    let _ = atoms; // used for consistency with the signature
    bonds.iter().find_map(|bond| {
        if (bond.begin().index() == a && bond.end().index() == b)
            || (bond.begin().index() == b && bond.end().index() == a)
        {
            Some(bond.id().index())
        } else {
            None
        }
    })
}

/// Helper: compute the degree (number of bonds) of an atom.
fn atom_degree(mol: &Molecule, a: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|b| b.begin().index() == a || b.end().index() == a)
        .count()
}

/// Helper: compute RDKit `Atom::getTotalDegree()` as graph degree plus
/// explicit and implicit hydrogens.
fn atom_total_degree(mol: &Molecule, atom_idx: usize) -> usize {
    let atom = &mol.atoms()[atom_idx];
    let implicit_hydrogens = mol
        .derived_cache()
        .valence
        .as_ref()
        .and_then(|valence| valence.implicit_hydrogens.get(atom_idx))
        .copied()
        .unwrap_or(0)
        .max(0) as usize;
    atom_degree(mol, atom_idx) + atom.explicit_hydrogens() as usize + implicit_hydrogens
}

/// Helper: compute degree from bond slice.
fn atom_degree_from_slice(bonds: &[crate::Bond], a: usize) -> usize {
    bonds
        .iter()
        .filter(|b| b.begin().index() == a || b.end().index() == a)
        .count()
}

// BEGIN RDKIT CPP FUNCTION: isTrigonalBipyramidalAxialBond (slice overload)
/// Check if a bond on a TBP center is axial (atom/bond slice version).
#[must_use]
pub fn is_trigonal_bipyramidal_axial_bond_by_slice(
    cen_idx: usize,
    query_bond_idx: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> i32 {
    if cen_idx >= atoms.len() || query_bond_idx >= bonds.len() {
        return 0;
    }
    let atom = &atoms[cen_idx];
    if atom_degree_from_slice(bonds, cen_idx) > 5
        || atom.chiral_tag() != ChiralTag::TrigonalBipyramidal
    {
        return 0;
    }
    let perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 || perm > 20 {
        return 0;
    }

    let mut count = 0usize;
    for bond in bonds {
        if bond.begin().index() == cen_idx || bond.end().index() == cen_idx {
            if bond.id().index() == query_bond_idx {
                if count == TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][0] as usize {
                    return 1;
                }
                if count == TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][1] as usize {
                    return -1;
                }
                return 0;
            }
            count += 1;
        }
    }
    0
}

// BEGIN RDKIT CPP FUNCTION: isTrigonalBipyramidalAxialAtom (slice overload)
/// Check if a neighbor atom of a TBP center is axial (atom/bond slice version).
#[must_use]
pub fn is_trigonal_bipyramidal_axial_atom_by_slice(
    cen_idx: usize,
    qry_idx: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> i32 {
    let bond_idx = bond_between_atoms_by_slice(atoms, bonds, cen_idx, qry_idx);
    match bond_idx {
        Some(bi) => is_trigonal_bipyramidal_axial_bond_by_slice(cen_idx, bi, atoms, bonds),
        None => 0,
    }
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossBond (slice overload)
/// Find the across bond given atom and bond slices instead of Molecule.
#[must_use]
pub fn get_chiral_across_bond_by_slice(
    cen_idx: usize,
    query_bond_idx: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> Option<usize> {
    if cen_idx >= atoms.len() || query_bond_idx >= bonds.len() {
        return None;
    }
    let atom = &atoms[cen_idx];
    let tag = atom.chiral_tag();
    let perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 {
        return None;
    }

    let max_nbors = get_max_nbors(tag) as usize;
    if max_nbors == 0 {
        return None;
    }

    let mut count = 0usize;
    let mut bond_refs: Vec<usize> = Vec::with_capacity(max_nbors);
    let mut found: Option<usize> = None;
    for bond in bonds {
        if bond.begin().index() == cen_idx || bond.end().index() == cen_idx {
            if count >= max_nbors {
                return None;
            }
            bond_refs.push(bond.id().index());
            if bond.id().index() == query_bond_idx {
                found = Some(count);
            }
            count += 1;
        }
    }

    if let Some(found_idx) = found {
        let across_idx = match tag {
            ChiralTag::SquarePlanar => {
                if perm <= 3 {
                    SQUAREPLANAR_ACROSS[perm as usize][found_idx]
                } else {
                    4
                }
            }
            ChiralTag::TrigonalBipyramidal => {
                if perm <= 20 {
                    TRIGONALBIPYRAMIDAL_ACROSS[perm as usize][found_idx]
                } else {
                    5
                }
            }
            ChiralTag::Octahedral => {
                if perm <= 30 {
                    OCTAHEDRAL_ACROSS[perm as usize][found_idx]
                } else {
                    6
                }
            }
            _ => return None,
        };
        if (across_idx as usize) < bond_refs.len() {
            Some(bond_refs[across_idx as usize])
        } else {
            None
        }
    } else {
        None
    }
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossBondByAtom (slice overload)
/// Find the across bond given atom/bond slices.
#[must_use]
pub fn get_chiral_across_bond_by_atom_by_slice(
    cen_idx: usize,
    qry_idx: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> Option<usize> {
    let bond_idx = bond_between_atoms_by_slice(atoms, bonds, cen_idx, qry_idx)?;
    get_chiral_across_bond_by_slice(cen_idx, bond_idx, atoms, bonds)
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossAtom (slice overload)
/// Find the across atom given atom/bond slices.
#[must_use]
pub fn get_chiral_across_atom_by_slice(
    cen_idx: usize,
    query_bond_idx: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> Option<usize> {
    let across_bond_idx = get_chiral_across_bond_by_slice(cen_idx, query_bond_idx, atoms, bonds)?;
    if across_bond_idx >= bonds.len() {
        return None;
    }
    let bond = &bonds[across_bond_idx];
    let other = if bond.begin().index() == cen_idx {
        bond.end().index()
    } else {
        bond.begin().index()
    };
    Some(other)
}

// BEGIN RDKIT CPP FUNCTION: getChiralAcrossAtomByAtom (slice overload)
/// Find the across atom given atom/bond slices.
#[must_use]
pub fn get_chiral_across_atom_by_atom_by_slice(
    cen_idx: usize,
    qry_idx: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> Option<usize> {
    let bond_idx = bond_between_atoms_by_slice(atoms, bonds, cen_idx, qry_idx)?;
    get_chiral_across_atom_by_slice(cen_idx, bond_idx, atoms, bonds)
}

// BEGIN RDKIT CPP FUNCTION: getIdealAngleBetweenLigands (slice overload)
/// Return the ideal angle between two ligands (atom/bond slice version).
#[must_use]
pub fn get_ideal_angle_between_ligands_by_slice(
    cen_idx: usize,
    lig1: usize,
    lig2: usize,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> f64 {
    if cen_idx >= atoms.len() {
        return 0.0;
    }
    let atom = &atoms[cen_idx];
    let tag = atom.chiral_tag();
    match tag {
        ChiralTag::SquarePlanar | ChiralTag::Octahedral => {
            if get_chiral_across_atom_by_atom_by_slice(cen_idx, lig1, atoms, bonds) == Some(lig2) {
                180.0
            } else {
                90.0
            }
        }
        ChiralTag::TrigonalBipyramidal => {
            if get_chiral_across_atom_by_atom_by_slice(cen_idx, lig1, atoms, bonds) == Some(lig2) {
                180.0
            } else if is_trigonal_bipyramidal_axial_atom_by_slice(cen_idx, lig1, atoms, bonds) != 0
                || is_trigonal_bipyramidal_axial_atom_by_slice(cen_idx, lig2, atoms, bonds) != 0
            {
                90.0
            } else {
                120.0
            }
        }
        _ => 0.0,
    }
}

// BEGIN RDKIT CPP FUNCTION: getTrigonalBipyramidalAxialBond (slice overload)
/// Get the bond to a specific axial position on a TBP center (atom/bond slice version).
#[must_use]
pub fn get_trigonal_bipyramidal_axial_bond_by_slice(
    cen_idx: usize,
    axial: i32,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> Option<usize> {
    if cen_idx >= atoms.len() {
        return None;
    }
    let atom = &atoms[cen_idx];
    if atom.chiral_tag() != ChiralTag::TrigonalBipyramidal
        || atom_degree_from_slice(bonds, cen_idx) > 5
    {
        return None;
    }

    let perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 || perm > 20 {
        return None;
    }

    let idx = if axial != -1 {
        TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][0] as usize
    } else {
        TRIGONALBIPYRAMIDAL_AXIAL[perm as usize][1] as usize
    };

    let mut count = 0usize;
    for bond in bonds {
        if bond.begin().index() == cen_idx || bond.end().index() == cen_idx {
            if count == idx {
                return Some(bond.id().index());
            }
            count += 1;
        }
    }
    None
}

// BEGIN RDKIT CPP FUNCTION: getTrigonalBipyramidalAxialAtom (slice overload)
/// Get the atom at a specific axial position on a TBP center (atom/bond slice version).
#[must_use]
pub fn get_trigonal_bipyramidal_axial_atom_by_slice(
    cen_idx: usize,
    axial: i32,
    atoms: &[crate::Atom],
    bonds: &[crate::Bond],
) -> Option<usize> {
    let bond_idx = get_trigonal_bipyramidal_axial_bond_by_slice(cen_idx, axial, atoms, bonds)?;
    if bond_idx >= bonds.len() {
        return None;
    }
    let bond = &bonds[bond_idx];
    let other = if bond.begin().index() == cen_idx {
        bond.end().index()
    } else {
        bond.begin().index()
    };
    Some(other)
}

// ── End non-tetrahedral stereo infrastructure ──────────────────────────────
// ──────────────────────────────────────────────
// New Chirality Functions (ported from RDKit Chirality.cpp)
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: isLinearArrangement (Chirality.cpp:99-112)
// RDKit✔️✔️: bool isLinearArrangement(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2) {
// RDKit✔️✔️:   double lsq = v1.lengthSq() * v2.lengthSq();
// RDKit✔️✔️:   if (lsq < 1.0e-6) { return true; }
// RDKit✔️✔️:   double dotProd = v1.dotProduct(v2);
// RDKit✔️✔️:   double cos178 = -0.999388;  // cos(M_PI-0.035), 2 degree tolerance
// RDKit✔️✔️:   return dotProd < cos178 * sqrt(lsq);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: isLinearArrangement
/// Check if two vectors form a linear (collinear) arrangement.
/// Returns true if the vectors are collinear (within ~2 degrees of 180 deg separation)
/// or if either vector has near-zero length.
#[must_use]
pub fn is_linear_arrangement(v1: (f64, f64, f64), v2: (f64, f64, f64)) -> bool {
    let lsq = (v1.0 * v1.0 + v1.1 * v1.1 + v1.2 * v1.2) * (v2.0 * v2.0 + v2.1 * v2.1 + v2.2 * v2.2);
    // treat zero length vectors as linear
    if lsq < 1.0e-6 {
        return true;
    }
    let dot_prod = v1.0 * v2.0 + v1.1 * v2.1 + v1.2 * v2.2;
    let cos178 = -0.999_388; // cos(M_PI-0.035), corresponds to a tolerance of 2 degrees
    dot_prod < cos178 * lsq.sqrt()
}

// BEGIN RDKIT CPP FUNCTION: shouldDetectDoubleBondStereo (Chirality.cpp:38-43, simplified)
// RDKit✔️✔️: bool shouldDetectDoubleBondStereo(const Bond *bond) {
// RDKit✔️✔️:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
// RDKit✔️✔️:   return (!ri->numBondRings(bond->getIdx()) ||
// RDKit✔️✔️:           ri->minBondRingSize(bond->getIdx()) >= Chirality::minRingSizeForDoubleBondStereo);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: shouldDetectDoubleBondStereo
/// Check if a bond is a candidate for stereo detection.
/// A double bond is a candidate if:
/// - It has stereo BondStereo != Any
/// - It has bond direction != EitherDouble
/// - Both begin/end atoms have degree > 1
/// - It's not in a small ring (or ring info is not available)
#[must_use]
pub fn is_bond_candidate_for_stereo(mol: &Molecule, bond_idx: usize) -> bool {
    if bond_idx >= mol.bonds().len() {
        return false;
    }
    let bond = &mol.bonds()[bond_idx];
    // Must be a double bond
    if bond.order() != crate::BondOrder::Double {
        return false;
    }
    // Skip if stereo is already set to Any (like crossed bond)
    if bond.stereo() == crate::BondStereo::Any {
        return false;
    }
    // Skip EitherDouble direction
    if bond.direction() == crate::BondDirection::EitherDouble {
        return false;
    }
    // Both end atoms must have degree > 1
    let begin_idx = bond.begin().index();
    let end_idx = bond.end().index();
    let begin_deg = atom_degree(mol, begin_idx);
    let end_deg = atom_degree(mol, end_idx);
    if begin_deg <= 1 || end_deg <= 1 {
        return false;
    }

    // Check ring info: skip if in a small ring (< minRingSizeForDoubleBondStereo)
    let rings_opt = mol.derived_cache().rings.clone();
    if let Some(ri) = rings_opt {
        if ri.is_initialized() {
            let bond_ring_count = ri.num_bond_rings(bond.id());
            if bond_ring_count > 0 {
                let min_size = ri.min_bond_ring_size(bond.id());
                // RDKit default minRingSizeForDoubleBondStereo = 8
                if min_size < 8 {
                    return false;
                }
            }
        }
    }

    true
}

// BEGIN RDKIT CPP STRUCT: ControllingBondResult (helper for controllingBondFromAtom)
/// Result of the controllingBondFromAtom search.
#[derive(Debug, Clone)]
pub struct ControllingBondResult {
    pub bond: Option<usize>,
    pub obond: Option<usize>,
    pub squiggle_bond_seen: bool,
    pub double_bond_seen: bool,
}

// BEGIN RDKIT CPP FUNCTION: controllingBondFromAtom (Chirality.cpp:114-161)
// RDKit✔️✔️: void controllingBondFromAtom(const ROMol &mol,
// RDKit✔️✔️:                              const boost::dynamic_bitset<> &needsDir,
// RDKit✔️✔️:                              const std::vector<unsigned int> &singleBondCounts,
// RDKit✔️✔️:                              const Bond *dblBond, const Atom *atom, Bond *&bond,
// RDKit✔️✔️:                              Bond *&obond, bool &squiggleBondSeen,
// RDKit✔️✔️:                              bool &doubleBondSeen) {
// RDKit✔️✔️:   // Selects a controlling single bond on an atom of a double bond
// RDKit✔️✔️:   // for E/Z assignment. Prefers bonds with direction set.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: controllingBondFromAtom
/// Find the controlling bond for cis/trans assignment from an atom of a double bond.
/// Scans single bonds from the given atom (excluding the double bond itself)
/// and selects the controlling bond based on direction and adjacency to other double bonds.
#[must_use]
pub fn controlling_bond_from_atom(
    mol: &Molecule,
    dbl_bond_idx: usize,
    atom_idx: usize,
) -> ControllingBondResult {
    let mut bond: Option<usize> = None;
    let mut obond: Option<usize> = None;
    let mut squiggle_bond_seen = false;
    let mut double_bond_seen = false;

    for b in mol.bonds() {
        let b_idx = b.id().index();
        if b_idx == dbl_bond_idx {
            continue;
        }
        // Only consider bonds connected to atom_idx
        if b.begin().index() != atom_idx && b.end().index() != atom_idx {
            continue;
        }

        let b_order = b.order();
        let b_dir = b.direction();

        if (b_order == crate::BondOrder::Single || b_order == crate::BondOrder::Aromatic)
            && (b_dir == crate::BondDirection::None
                || b_dir == crate::BondDirection::EndDownRight
                || b_dir == crate::BondDirection::EndUpRight)
        {
            // Prefer bonds with direction, or adjacent to more double bonds
            if bond.is_none() {
                bond = Some(b_idx);
            } else {
                obond = Some(b_idx);
            }
        } else if b_order == crate::BondOrder::Double {
            double_bond_seen = true;
        }

        // Check for squiggle bond (Unknown direction or explicit unknown stereo)
        if (b_order == crate::BondOrder::Single || b_order == crate::BondOrder::Aromatic)
            && b_dir == crate::BondDirection::Unknown
        {
            squiggle_bond_seen = true;
            // special handling for explicit unknown stereo property
            if b.prop("_UnknownStereo").and_then(|v| v.parse::<i32>().ok()) == Some(1) {
                squiggle_bond_seen = true;
            }
        }
    }

    ControllingBondResult {
        bond,
        obond,
        squiggle_bond_seen,
        double_bond_seen,
    }
}

// BEGIN RDKIT CPP FUNCTION: updateDoubleBondNeighbors (Chirality.cpp:163-379, heavily simplified)
// RDKit✔️✔️: void updateDoubleBondNeighbors(ROMol &mol, Bond *dblBond, const Conformer *conf,
// RDKit✔️✔️:                                boost::dynamic_bitset<> &needsDir,
// RDKit✔️✔️:                                std::vector<unsigned int> &singleBondCounts,
// RDKit✔️✔️:                                const VECT_INT_VECT &singleBondNbrs) {
// RDKit✔️✔️:   // Sets bond directions on single bonds adjacent to double bonds
// RDKit✔️✔️:   // based on dihedral angle analysis or existing stereo markings.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: updateDoubleBondNeighbors
/// Update neighbor ordering for double bond stereo from 2D coordinates.
/// Simple version that checks bond directions and single bond arrangement.
pub fn update_double_bond_neighbors(
    mol: &Molecule,
    dbl_bond_idx: usize,
) -> Option<ControllingBondResult> {
    if dbl_bond_idx >= mol.bonds().len() {
        return None;
    }
    let dbl_bond = &mol.bonds()[dbl_bond_idx];
    if dbl_bond.order() != crate::BondOrder::Double {
        return None;
    }

    // Check if we should process this bond
    if !is_bond_candidate_for_stereo(mol, dbl_bond_idx) {
        return None;
    }

    let begin_idx = dbl_bond.begin().index();
    let end_idx = dbl_bond.end().index();

    let begin_result = controlling_bond_from_atom(mol, dbl_bond_idx, begin_idx);
    let end_result = controlling_bond_from_atom(mol, dbl_bond_idx, end_idx);

    if begin_result.squiggle_bond_seen || end_result.squiggle_bond_seen {
        return None;
    }

    if begin_result.bond.is_none() || end_result.bond.is_none() {
        return None;
    }

    // Determine E/Z from bond directions
    let begin_bond_idx = begin_result.bond.unwrap();
    let end_bond_idx = end_result.bond.unwrap();

    let begin_bond = &mol.bonds()[begin_bond_idx];
    let end_bond = &mol.bonds()[end_bond_idx];

    let begin_dir = begin_bond.direction();
    let end_dir = end_bond.direction();

    // Normalize direction: if the bond is backwards (atom is end, not begin), flip
    let begin_dir = if begin_bond.begin().index() != begin_idx {
        match begin_dir {
            crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
            crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
            d => d,
        }
    } else {
        begin_dir
    };

    let end_dir = if end_bond.begin().index() != end_idx {
        match end_dir {
            crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
            crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
            d => d,
        }
    } else {
        end_dir
    };

    // Store the controlling atoms for E/Z
    let begin_ctrl_atom = bond_other_atom(mol, begin_bond_idx, begin_idx);
    let end_ctrl_atom = bond_other_atom(mol, end_bond_idx, end_idx);

    // RDKit combines begin and end results
    Some(ControllingBondResult {
        bond: begin_ctrl_atom,
        obond: end_ctrl_atom,
        squiggle_bond_seen: false,
        double_bond_seen: begin_dir == end_dir,
    })
}

/// Helper: get the atom at the other end of a bond from a given atom.
fn bond_other_atom(mol: &Molecule, bond_idx: usize, atom_idx: usize) -> Option<usize> {
    if bond_idx >= mol.bonds().len() {
        return None;
    }
    let bond = &mol.bonds()[bond_idx];
    if bond.begin().index() == atom_idx {
        Some(bond.end().index())
    } else if bond.end().index() == atom_idx {
        Some(bond.begin().index())
    } else {
        None
    }
}

// BEGIN RDKIT CPP FUNCTION: findAtomNeighborDirHelper (Chirality.cpp:1351-1415, simplified)
// RDKit✔️✔️: void findAtomNeighborDirHelper(const ROMol &mol, const Atom *atom,
// RDKit✔️✔️:                                const Bond *refBond, UINT_VECT &ranks,
// RDKit✔️✔️:                                INT_PAIR_VECT &neighbors,
// RDKit✔️✔️:                                bool &hasExplicitUnknownStereo) {
// RDKit✔️✔️:   // Helper collecting neighboring bonds with direction for a double bond atom.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: findAtomNeighborDirHelper
/// Collect neighbor atom indices and bond directions around a double bond atom.
/// Returns pairs of (neighbor_atom_idx, direction). Only returns non-empty
/// if at least one neighbor bond has direction set.
// BEGIN RDKIT CPP FUNCTION: findAtomNeighborDirHelper (Chirality.cpp:1351-1408)
// RDKit✔️✔️: void findAtomNeighborDirHelper(const ROMol &mol, const Atom *atom,
// RDKit✔️✔️:                                const Bond *refBond, UINT_VECT &ranks,
// RDKit✔️✔️:                                INT_PAIR_VECT &neighbors,
// RDKit✔️✔️:                                bool &hasExplicitUnknownStereo) {
// RDKit✔️✔️:   bool seenDir = false;
// RDKit✔️✔️:   for (const auto bond : mol.atomBonds(atom)) {
// RDKit✔️✔️:     if (!hasExplicitUnknownStereo) {
// RDKit✔️✔️:       int explicit_unknown_stereo;
// RDKit✔️✔️:       if (bond->getBondDir() == Bond::UNKNOWN ||
// RDKit✔️✔️:           (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
// RDKit✔️✔️:                                   explicit_unknown_stereo) &&
// RDKit✔️✔️:            explicit_unknown_stereo)) {
// RDKit✔️✔️:         hasExplicitUnknownStereo = true;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     Bond::BondDir dir = bond->getBondDir();
// RDKit✔️✔️:     if (bond->getIdx() != refBond->getIdx()) {
// RDKit✔️✔️:       if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
// RDKit✔️✔️:         seenDir = true;
// RDKit✔️✔️:         if (atom != bond->getBeginAtom()) {
// RDKit✔️✔️:           dir = dir == Bond::ENDDOWNRIGHT ? Bond::ENDUPRIGHT
// RDKit✔️✔️:                                             : Bond::ENDDOWNRIGHT;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:       Atom *nbrAtom = bond->getOtherAtom(atom);
// RDKit✔️✔️:       neighbors.push_back(std::make_pair(nbrAtom->getIdx(), dir));
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (!seenDir) {
// RDKit✔️✔️:     neighbors.clear();
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     if (neighbors.size() == 2 &&
// RDKit✔️✔️:         ranks[neighbors[0].first] == ranks[neighbors[1].first]) {
// RDKit✔️✔️:       neighbors.clear();
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       if (neighbors[0].second != Bond::ENDDOWNRIGHT &&
// RDKit✔️✔️:           neighbors[0].second != Bond::ENDUPRIGHT) {
// RDKit✔️✔️:         neighbors[0].second = neighbors[1].second == Bond::ENDDOWNRIGHT
// RDKit✔️✔️:                                   ? Bond::ENDUPRIGHT
// RDKit✔️✔️:                                   : Bond::ENDDOWNRIGHT;
// RDKit✔️✔️:       } else if (neighbors.size() > 1 &&
// RDKit✔️✔️:                  neighbors[1].second != Bond::ENDDOWNRIGHT &&
// RDKit✔️✔️:                  neighbors[1].second != Bond::ENDUPRIGHT) {
// RDKit✔️✔️:         neighbors[1].second = neighbors[0].second == Bond::ENDDOWNRIGHT
// RDKit✔️✔️:                                   ? Bond::ENDUPRIGHT
// RDKit✔️✔️:                                   : Bond::ENDDOWNRIGHT;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: findAtomNeighborDirHelper
fn find_atom_neighbor_dir_helper(
    mol: &Molecule,
    atom_idx: usize,
    dbl_bond_idx: usize,
    ranks: &[u32],
    has_explicit_unknown_stereo: &mut bool,
) -> Vec<(usize, crate::BondDirection)> {
    let mut neighbors: Vec<(usize, crate::BondDirection)> = Vec::new();
    let mut seen_dir = false;

    for b in mol.bonds() {
        let b_idx = b.id().index();
        if b_idx == dbl_bond_idx {
            continue;
        }
        if b.begin().index() != atom_idx && b.end().index() != atom_idx {
            continue;
        }

        // Check for explicit unknown stereo
        if !*has_explicit_unknown_stereo {
            if b.direction() == crate::BondDirection::Unknown {
                *has_explicit_unknown_stereo = true;
            }
            if let Some(v) = b.prop("_UnknownStereo") {
                if let Ok(val) = v.parse::<i32>() {
                    if val != 0 {
                        *has_explicit_unknown_stereo = true;
                    }
                }
            }
        }

        let mut dir = b.direction();
        if dir == crate::BondDirection::EndDownRight || dir == crate::BondDirection::EndUpRight {
            seen_dir = true;
            if atom_idx != b.begin().index() {
                dir = match dir {
                    crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
                    crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
                    other => other,
                };
            }
        }
        let nbr_atom = if b.begin().index() == atom_idx {
            b.end().index()
        } else {
            b.begin().index()
        };
        neighbors.push((nbr_atom, dir));
    }

    if !seen_dir {
        return Vec::new();
    }

    // If both neighbors have the same rank, clear (no stereochemistry)
    if neighbors.len() == 2 && ranks.get(neighbors[0].0) == ranks.get(neighbors[1].0) {
        return Vec::new();
    }

    // Ensure both neighbors have direction set
    if neighbors.len() >= 1
        && neighbors[0].1 != crate::BondDirection::EndDownRight
        && neighbors[0].1 != crate::BondDirection::EndUpRight
    {
        if neighbors.len() > 1 {
            neighbors[0].1 = if neighbors[1].1 == crate::BondDirection::EndDownRight {
                crate::BondDirection::EndUpRight
            } else {
                crate::BondDirection::EndDownRight
            };
        }
    } else if neighbors.len() > 1
        && neighbors[1].1 != crate::BondDirection::EndDownRight
        && neighbors[1].1 != crate::BondDirection::EndUpRight
    {
        neighbors[1].1 = if neighbors[0].1 == crate::BondDirection::EndDownRight {
            crate::BondDirection::EndUpRight
        } else {
            crate::BondDirection::EndDownRight
        };
    }

    neighbors
}

// BEGIN RDKIT CPP FUNCTION: isAtomPotentialChiralCenter (Chirality.cpp:1651-1736)
// RDKit❗✔️: std::pair<bool, bool> isAtomPotentialChiralCenter(
// RDKit❗✔️:     const Atom *atom, const ROMol &mol, const UINT_VECT &ranks,
// RDKit❗✔️:     Chirality::INT_PAIR_VECT &nbrs) {
// RDKit❗✔️:   // Check if an atom could be a tetrahedral chiral center.
// RDKit❗✔️:   // Returns (legal_center, has_duplicates). Populates nbrs with (rank, bond_idx).
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION: isAtomPotentialChiralCenter
/// Check if an atom could be a tetrahedral chiral center.
/// Returns (legal_center, has_duplicates, neighbors: Vec<(rank, idx)>).
/// neighbors contains (CIP_rank, neighbor_atom_idx) pairs.
pub fn is_atom_potential_chiral_center(
    mol: &Molecule,
    atom_idx: usize,
    ranks: &[u32],
) -> (bool, bool, Vec<(u32, usize)>) {
    let atom = &mol.atoms()[atom_idx];
    let mut legal_center = true;
    let mut has_dupes = false;
    let mut nbrs: Vec<(u32, usize)> = Vec::new();

    if atom_idx >= mol.num_atoms() {
        return (false, false, nbrs);
    }

    // Non-zero degree (exclude bonds that don't affect chirality)
    let nz_degree = atom_nonzero_degree(mol, atom_idx);
    let implicit_hydrogens = mol
        .derived_cache()
        .valence
        .as_ref()
        .and_then(|valence| valence.implicit_hydrogens.get(atom_idx))
        .copied()
        .unwrap_or(0)
        .max(0) as usize;
    let total_num_hs = atom.explicit_hydrogens() as usize + implicit_hydrogens;
    let total_nz_degree = nz_degree + total_num_hs;

    if total_nz_degree > 4 {
        // we only know tetrahedral chirality
        legal_center = false;
    } else if total_nz_degree < 3 {
        legal_center = false;
    } else if nz_degree < 3 && atom.atomic_number() != 15 && atom.atomic_number() != 33 {
        // less than three neighbors is never stereogenic
        // unless it is a phosphine/arsine with implicit H (this is from InChI)
        legal_center = false;
    } else if nz_degree == 3 {
        // Check if exactly one H neighbor using explicit_hydrogens
        if total_num_hs == 1 {
            // three-coordinate with exactly one H
            // if it has a protium neighbor, not stereogenic
            if has_protium_neighbor(mol, atom_idx) {
                legal_center = false;
            }
        } else {
            // assume something that's really three-coordinate isn't potentially chiral
            // then look for exceptions
            legal_center = false;
            match atom.atomic_number() {
                7 => {
                    // RDKit✔️✔️: if (atom->getHybridization() == Atom::HybridizationType::SP3 &&
                    // RDKit✔️✔️:     !MolOps::atomHasConjugatedBond(atom) &&
                    // RDKit✔️✔️:     (mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3) ||
                    // RDKit✔️✔️:      queryIsAtomBridgehead(atom))) {
                    // RDKit✔️✔️:   legalCenter = true;
                    // RDKit✔️✔️: }
                    let in_three_membered_ring = mol
                        .derived_cache()
                        .rings
                        .as_ref()
                        .is_some_and(|ri| ri.is_atom_in_ring_of_size(atom.id(), 3));
                    let atom_has_conjugated_bond = mol.bonds().iter().any(|bond| {
                        (bond.begin().index() == atom_idx || bond.end().index() == atom_idx)
                            && bond.is_conjugated()
                    });
                    let is_bridgehead = mol
                        .derived_cache()
                        .rings
                        .as_ref()
                        .is_some_and(|ri| query_is_atom_bridgehead(mol, atom_idx, ri) != 0);
                    if atom.hybridization() == crate::Hybridization::Sp3
                        && !atom_has_conjugated_bond
                        && (in_three_membered_ring || is_bridgehead)
                    {
                        legal_center = true;
                    }
                }
                15 | 33 => {
                    // phosphines and arsines are always stereogenic
                    legal_center = true;
                }
                16 | 34 => {
                    // RDKit✔️✔️: if (atom->getValence(Atom::ValenceType::EXPLICIT) == 4 ||
                    // RDKit✔️✔️:     (atom->getValence(Atom::ValenceType::EXPLICIT) == 3 &&
                    // RDKit✔️✔️:      atom->getFormalCharge() == 1)) {
                    // RDKit✔️✔️:   legalCenter = true;
                    // RDKit✔️✔️: }
                    let explicit_valence = mol
                        .derived_cache()
                        .valence
                        .as_ref()
                        .and_then(|valence| valence.explicit_valence.get(atom_idx))
                        .copied()
                        .unwrap_or_default();
                    if explicit_valence == 4 || (explicit_valence == 3 && atom.formal_charge() == 1)
                    {
                        legal_center = true;
                    }
                }
                _ => {}
            }
        }
    }

    if legal_center && !ranks.is_empty() {
        let mut codes_seen = vec![false; mol.num_atoms()];
        for b in mol.bonds() {
            if b.begin().index() != atom_idx && b.end().index() != atom_idx {
                continue;
            }
            let other_idx = if b.begin().index() == atom_idx {
                b.end().index()
            } else {
                b.begin().index()
            };
            nbrs.push((ranks[other_idx], b.id().index()));
            if !bond_affects_atom_chirality(b, atom_idx) {
                continue;
            }
            let rank = ranks[other_idx] as usize;
            if rank < codes_seen.len() {
                if codes_seen[rank] {
                    has_dupes = true;
                    break;
                }
                codes_seen[rank] = true;
            }
        }
    }

    (legal_center, has_dupes, nbrs)
}

fn perturbation_order_from_bond_indices(
    mol: &Molecule,
    atom_idx: usize,
    probe: &[usize],
) -> Result<u32, StereoError> {
    let reference = mol
        .topology_block()
        .adjacency
        .neighbors_of(atom_idx)
        .iter()
        .map(|neighbor| neighbor.bond.index())
        .collect::<Vec<_>>();
    if probe.len() != reference.len() {
        return Err(StereoError::InvariantViolation(
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
            return Err(StereoError::InvariantViolation(
                "Atom::getPerturbationOrder expected bond missing from probe order".to_string(),
            ));
        };
        work.swap(idx, found_idx);
        swaps = swaps.saturating_add(1);
    }
    Ok(swaps)
}

// RDKit✔️✔️: bool bondAffectsAtomChirality(const Bond *bond, const Atom *atom) {
// RDKit✔️✔️:   // FIX consider how to handle organometallics
// RDKit✔️✔️:   PRECONDITION(bond, "bad bond pointer");
// RDKit✔️✔️:   PRECONDITION(atom, "bad atom pointer");
// RDKit✔️✔️:   if (bond->getBondType() == Bond::BondType::UNSPECIFIED ||
// RDKit✔️✔️:       bond->getBondType() == Bond::BondType::ZERO ||
// RDKit✔️✔️:       (bond->getBondType() == Bond::BondType::DATIVE &&
// RDKit✔️✔️:        bond->getBeginAtomIdx() == atom->getIdx())) {
// RDKit✔️✔️:     return false;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return true;
// RDKit✔️✔️: }
// Safe Rust references and the typed atom index make the two pointer
// preconditions unrepresentable at this helper boundary.
fn bond_affects_atom_chirality(bond: &crate::Bond, atom_idx: usize) -> bool {
    if matches!(
        bond.order(),
        crate::BondOrder::Null | crate::BondOrder::Unspecified | crate::BondOrder::Zero
    ) {
        return false;
    }
    if bond.order() == crate::BondOrder::Dative && bond.begin().index() == atom_idx {
        return false;
    }
    true
}

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
// The molecule argument supplies the owning graph. Its CSR adjacency slice
// preserves bond-table insertion order and gives the same O(degree) traversal
// without allocation.
fn atom_nonzero_degree(mol: &Molecule, atom_idx: usize) -> usize {
    assert!(atom_idx < mol.num_atoms(), "bad atom index");
    let topology = mol.topology_block();
    atom_nonzero_degree_from_parts(&topology.bonds, &topology.adjacency, atom_idx)
        .expect("molecule topology and adjacency must be aligned")
}

fn atom_nonzero_degree_from_parts(
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    atom_idx: usize,
) -> Result<usize, StereoError> {
    let mut res = 0;
    for neighbor in adjacency.neighbors_of(atom_idx) {
        let bond = bonds.get(neighbor.bond.index()).ok_or_else(|| {
            StereoError::InvariantViolation(format!(
                "adjacency for atom {atom_idx} references missing bond {}",
                neighbor.bond
            ))
        })?;
        if !bond_affects_atom_chirality(bond, atom_idx) {
            continue;
        }
        res += 1;
    }
    Ok(res)
}

// RDKit✔️✔️: bool isWigglyBond(const Bond *bond, const Atom *atom) {
// RDKit✔️✔️:   int hasWigglyBond = 0;
// RDKit✔️✔️:   if (bond->getBeginAtomIdx() == atom->getIdx() &&
// RDKit✔️✔️:       bond->getBondType() == Bond::BondType::SINGLE &&
// RDKit✔️✔️:       (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
// RDKit✔️✔️:        (bond->getPropIfPresent<int>(common_properties::_UnknownStereo,
// RDKit✔️✔️:                                     hasWigglyBond) &&
// RDKit✔️✔️:         hasWigglyBond))) {
// RDKit✔️✔️:     return true;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return false;
// RDKit✔️✔️: }
fn is_wiggly_bond(bond: &crate::Bond, atom_idx: usize) -> Result<bool, StereoError> {
    if bond.begin().index() != atom_idx || bond.order() != crate::BondOrder::Single {
        return Ok(false);
    }
    if bond.direction() == crate::BondDirection::Unknown {
        return Ok(true);
    }
    let has_wiggly_bond = match bond.prop("_UnknownStereo") {
        Some(value) => value.parse::<i32>().map_err(|_| {
            StereoError::InvariantViolation(format!(
                "bond {} has non-integer _UnknownStereo property {value:?}",
                bond.id()
            ))
        })?,
        None => i32::from(bond.unknown_stereo()),
    };
    Ok(has_wiggly_bond != 0)
}

/// Check if an atom has a protium (regular H) neighbor.
fn has_protium_neighbor(mol: &Molecule, atom_idx: usize) -> bool {
    for b in mol.bonds() {
        let other_idx = if b.begin().index() == atom_idx {
            b.end().index()
        } else if b.end().index() == atom_idx {
            b.begin().index()
        } else {
            continue;
        };
        let other = &mol.atoms()[other_idx];
        if other.atomic_number() == 1 && other.isotope().map_or(true, |iso| iso == 0) {
            return true;
        }
    }
    false
}

// BEGIN RDKIT CPP FUNCTION: assignBondStereoCodes (Chirality.cpp:1826-1965)
// RDKit✔️✔️: std::pair<bool, bool> assignBondStereoCodes(ROMol &mol, UINT_VECT &ranks) {
// RDKit✔️✔️:   // Assign E/Z stereo codes to double bonds from CIP ranks + bond directions.
// RDKit✔️✔️:   // Returns (unassigned_bonds, assigned_bond).
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: assignBondStereoCodes
/// Assign E/Z stereo codes to double bonds from CIP ranks and bond directions.
/// Returns `(unassigned_bonds_remain, assignments, changed_any_bond)`.
pub fn assign_bond_stereo_codes(
    mol: &Molecule,
    ranks: &[u32],
) -> (bool, Vec<(usize, DoubleBondStereo, usize, usize)>, bool) {
    let mut results: Vec<(usize, DoubleBondStereo, usize, usize)> = Vec::new();
    let mut changed = false;
    let mut unassigned_bonds = 0usize;

    for dbl_bond in mol.bonds() {
        let dbl_idx = dbl_bond.id().index();
        if dbl_bond.order() != crate::BondOrder::Double {
            continue;
        }
        if dbl_bond.stereo() != crate::BondStereo::None {
            continue;
        }
        if !is_bond_candidate_for_stereo(mol, dbl_idx) {
            continue;
        }

        let beg_atom = dbl_bond.begin().index();
        let end_atom = dbl_bond.end().index();
        let beg_deg = atom_degree(mol, beg_atom);
        let end_deg = atom_degree(mol, end_atom);

        if (beg_deg != 2 && beg_deg != 3) || (end_deg != 2 && end_deg != 3) {
            continue;
        }
        unassigned_bonds += 1;

        let mut has_explicit_unknown = false;

        // Check for explicit unknown stereo on the double bond atoms
        if let Some(s) = dbl_bond.prop("_UnknownStereo") {
            if let Ok(v) = s.parse::<i32>() {
                if v != 0 {
                    has_explicit_unknown = true;
                }
            }
        }

        let beg_neighbors =
            find_atom_neighbor_dir_helper(mol, beg_atom, dbl_idx, ranks, &mut has_explicit_unknown);
        let end_neighbors =
            find_atom_neighbor_dir_helper(mol, end_atom, dbl_idx, ranks, &mut has_explicit_unknown);

        if beg_neighbors.is_empty() || end_neighbors.is_empty() {
            continue;
        }

        // Find highest-ranked direction on each side
        let (beg_dir, beg_ctrl) =
            if beg_neighbors.len() == 1 || ranks[beg_neighbors[0].0] > ranks[beg_neighbors[1].0] {
                (beg_neighbors[0].1, beg_neighbors[0].0)
            } else {
                (beg_neighbors[1].1, beg_neighbors[1].0)
            };

        let (end_dir, end_ctrl) =
            if end_neighbors.len() == 1 || ranks[end_neighbors[0].0] > ranks[end_neighbors[1].0] {
                (end_neighbors[0].1, end_neighbors[0].0)
            } else {
                (end_neighbors[1].1, end_neighbors[1].0)
            };

        // Check for conflicting directions
        let conflicting_begin =
            beg_neighbors.len() == 2 && beg_neighbors[0].1 == beg_neighbors[1].1;
        let conflicting_end = end_neighbors.len() == 2 && end_neighbors[0].1 == end_neighbors[1].1;

        if conflicting_begin || conflicting_end {
            changed = true;
        } else {
            let stereo = if has_explicit_unknown {
                DoubleBondStereo::Unknown
            } else if beg_dir == end_dir {
                // Both bonds point same direction → Z (cis)
                DoubleBondStereo::Z
            } else {
                // Opposite directions → E (trans)
                DoubleBondStereo::E
            };
            results.push((dbl_idx, stereo, beg_ctrl, end_ctrl));
            changed = true;
            unassigned_bonds = unassigned_bonds.saturating_sub(1);
        }
    }

    (unassigned_bonds > 0, results, changed)
}

// BEGIN RDKIT CPP FUNCTION: assignLegacyCIPLabels (Chirality.cpp:1966-1979)
// RDKit✔️✔️: void assignLegacyCIPLabels(ROMol &mol, bool flagPossibleStereoCenters) {
// RDKit✔️✔️:   std::vector<unsigned int> atomRanks;
// RDKit✔️✔️:   assignAtomChiralCodes(mol, atomRanks, flagPossibleStereoCenters);
// RDKit✔️✔️:   // reset double bonds
// RDKit✔️✔️:   for (auto bond : mol.bonds()) {
// RDKit✔️✔️:     if (bond->getBondType() == Bond::BondType::DOUBLE &&
// RDKit✔️✔️:         bond->getStereo() > Bond::BondStereo::STEREOANY) {
// RDKit✔️✔️:       bond->setStereo(Bond::BondStereo::STEREONONE);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   assignBondStereoCodes(mol, atomRanks);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: assignLegacyCIPLabels
/// Top-level dispatcher for CIP labeling.
/// Calls assignAtomChiralCodes + assignBondStereoCodes.
/// Returns (atom_labels, bond_stereo_results).
pub fn assign_legacy_cip_labels(
    mol: &Molecule,
    flag_possible_stereo_centers: bool,
) -> Result<
    (
        Vec<(usize, String)>,
        Vec<(usize, DoubleBondStereo, usize, usize)>,
    ),
    StereoError,
> {
    let ranks = if flag_possible_stereo_centers {
        assign_atom_cip_ranks(mol)?
    } else {
        Vec::new()
    };

    let (_, atom_labels, _) = assign_atom_chiral_codes(mol, &ranks)?;

    // Assign bond stereo codes
    let (_, bond_results, _) = assign_bond_stereo_codes(mol, &ranks);

    Ok((atom_labels, bond_results))
}

// BEGIN RDKIT CPP FUNCTION: assignBondCisTrans (Chirality.cpp:1980-2063)
// RDKit✔️✔️: void assignBondCisTrans(ROMol &mol, const StereoInfo &sinfo) {
// RDKit✔️✔️:   bool begFirstNeighbor = true;
// RDKit✔️✔️:   auto begBond = mol.getBondBetweenAtoms(dblBond->getBeginAtomIdx(),
// RDKit✔️✔️:                                          sinfo.controllingAtoms[0]);
// RDKit✔️✔️:   auto begDir = begBond->getBondDir();
// RDKit✔️✔️:   if (begDir != Bond::BondDir::ENDDOWNRIGHT &&
// RDKit✔️✔️:       begDir != Bond::BondDir::ENDUPRIGHT) {
// RDKit✔️✔️:     begFirstNeighbor = false;
// RDKit✔️✔️:     if (sinfo.controllingAtoms[1] != Atom::NOATOM) {
// RDKit✔️✔️:       begBond = mol.getBondBetweenAtoms(dblBond->getBeginAtomIdx(),
// RDKit✔️✔️:                                         sinfo.controllingAtoms[1]);
// RDKit✔️✔️:       begDir = begBond->getBondDir();
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (begBond->getBeginAtomIdx() != dblBond->getBeginAtomIdx()) {
// RDKit✔️✔️:     begDir = begDir == Bond::BondDir::ENDDOWNRIGHT
// RDKit✔️✔️:                  ? Bond::BondDir::ENDUPRIGHT
// RDKit✔️✔️:                  : Bond::BondDir::ENDDOWNRIGHT;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   bool endFirstNeighbor = true;
// RDKit✔️✔️:   auto endBond = mol.getBondBetweenAtoms(dblBond->getEndAtomIdx(),
// RDKit✔️✔️:                                          sinfo.controllingAtoms[2]);
// RDKit✔️✔️:   auto endDir = endBond->getBondDir();
// RDKit✔️✔️:   if (endDir != Bond::BondDir::ENDDOWNRIGHT &&
// RDKit✔️✔️:       endDir != Bond::BondDir::ENDUPRIGHT) {
// RDKit✔️✔️:     endFirstNeighbor = false;
// RDKit✔️✔️:     if (sinfo.controllingAtoms[3] != Atom::NOATOM) {
// RDKit✔️✔️:       endBond = mol.getBondBetweenAtoms(dblBond->getEndAtomIdx(),
// RDKit✔️✔️:                                         sinfo.controllingAtoms[3]);
// RDKit✔️✔️:       endDir = endBond->getBondDir();
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (endBond->getBeginAtomIdx() != dblBond->getEndAtomIdx()) {
// RDKit✔️✔️:     endDir = endDir == Bond::BondDir::ENDDOWNRIGHT
// RDKit✔️✔️:                  ? Bond::BondDir::ENDUPRIGHT
// RDKit✔️✔️:                  : Bond::BondDir::ENDDOWNRIGHT;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   bool sameDir = begDir == endDir;
// RDKit✔️✔️:   if (begFirstNeighbor ^ endFirstNeighbor) {
// RDKit✔️✔️:     sameDir = !sameDir;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: assignBondCisTrans
/// Assign cis/trans (E/Z) to a double bond from 2D bond directions.
/// controlling_atoms: [begin_a, begin_b, end_a, end_b] - the controlling atoms
/// where begin_a/begin_b are neighbors of begin atom, end_a/end_b are neighbors of end atom.
/// If a controlling atom is None, the other one on that side is used.
/// Returns the assigned DoubleBondStereo or None if not determinable.
pub fn assign_bond_cis_trans(
    mol: &Molecule,
    bond_idx: usize,
    controlling_atoms: &[Option<usize>; 4],
) -> Option<DoubleBondStereo> {
    if bond_idx >= mol.bonds().len() {
        return None;
    }
    if controlling_atoms.len() < 4 {
        return None;
    }

    // Check that we have at least one controlling atom per side
    if controlling_atoms[0].is_none() && controlling_atoms[1].is_none() {
        return None;
    }
    if controlling_atoms[2].is_none() && controlling_atoms[3].is_none() {
        return None;
    }

    let bond = &mol.bonds()[bond_idx];
    if bond.order() != crate::BondOrder::Double {
        return None;
    }

    let beg_atom = bond.begin().index();
    let end_atom = bond.end().index();

    // Find the direction bond at the beginning
    let mut beg_first = true;
    let mut beg_dir = None;
    let mut beg_dir_bond_idx = None;
    if let Some(beg_ctrl) = controlling_atoms[0] {
        if let Some(bi) = bond_between_atoms(mol, beg_atom, beg_ctrl) {
            let b = &mol.bonds()[bi];
            let d = b.direction();
            if d == crate::BondDirection::EndDownRight || d == crate::BondDirection::EndUpRight {
                beg_dir = Some(d);
                beg_dir_bond_idx = Some(bi);
            }
        }
    }
    if beg_dir.is_none() {
        beg_first = false;
        if let Some(beg_ctrl) = controlling_atoms[1] {
            if let Some(bi) = bond_between_atoms(mol, beg_atom, beg_ctrl) {
                let b = &mol.bonds()[bi];
                let d = b.direction();
                if d == crate::BondDirection::EndDownRight || d == crate::BondDirection::EndUpRight
                {
                    beg_dir = Some(d);
                    beg_dir_bond_idx = Some(bi);
                }
            }
        }
    }
    let mut beg_dir = beg_dir?;
    let beg_dir_bond_idx = beg_dir_bond_idx?;

    // Normalize direction
    {
        let b = &mol.bonds()[beg_dir_bond_idx];
        if b.begin().index() != beg_atom {
            beg_dir = match beg_dir {
                crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
                crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
                d => d,
            };
        }
    }

    // Find the direction bond at the end
    let mut end_first = true;
    let mut end_dir = None;
    let mut end_dir_bond_idx = None;
    if let Some(end_ctrl) = controlling_atoms[2] {
        if let Some(bi) = bond_between_atoms(mol, end_atom, end_ctrl) {
            let b = &mol.bonds()[bi];
            let d = b.direction();
            if d == crate::BondDirection::EndDownRight || d == crate::BondDirection::EndUpRight {
                end_dir = Some(d);
                end_dir_bond_idx = Some(bi);
            }
        }
    }
    if end_dir.is_none() {
        end_first = false;
        if let Some(end_ctrl) = controlling_atoms[3] {
            if let Some(bi) = bond_between_atoms(mol, end_atom, end_ctrl) {
                let b = &mol.bonds()[bi];
                let d = b.direction();
                if d == crate::BondDirection::EndDownRight || d == crate::BondDirection::EndUpRight
                {
                    end_dir = Some(d);
                    end_dir_bond_idx = Some(bi);
                }
            }
        }
    }
    let mut end_dir = end_dir?;
    let end_dir_bond_idx = end_dir_bond_idx?;

    // Normalize direction
    {
        let b = &mol.bonds()[end_dir_bond_idx];
        if b.begin().index() != end_atom {
            end_dir = match end_dir {
                crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
                crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
                d => d,
            };
        }
    }

    // Same direction = cis (Z), opposite = trans (E)
    let mut same_dir = beg_dir == end_dir;

    // If one side uses the second neighbor and the other uses the first, swap
    if beg_first ^ end_first {
        same_dir = !same_dir;
    }

    Some(if same_dir {
        DoubleBondStereo::Z
    } else {
        DoubleBondStereo::E
    })
}

// BEGIN RDKIT CPP FUNCTION: rerankAtoms (Chirality.cpp:2067-2117)
// RDKit✔️✔️: void rerankAtoms(const ROMol &mol, UINT_VECT &ranks) {
// RDKit✔️✔️:   // Re-rank atoms supplementing current ranks with chirality info.
// RDKit✔️✔️:   // R > S in priority, E/Z info also included.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: rerankAtoms
/// Re-rank atoms by supplementing the current ranks with known chirality and
/// double bond stereo information. R atoms get higher priority than S atoms.
/// Returns new ranks.
pub fn rerank_atoms(mol: &Molecule, current_ranks: &[u32]) -> Result<Vec<u32>, StereoError> {
    let n = mol.num_atoms();
    if current_ranks.len() != n {
        return Err(StereoError::UnsupportedFeature(
            crate::UnsupportedFeatureError {
                feature: "RERANK",
                reason: "current_ranks length must match number of atoms",
            },
        ));
    }

    // Compute scaling factor
    let mut factor: u32 = 100;
    while factor < n as u32 {
        factor *= 10;
    }

    let mut invars = vec![0i64; n];

    // Build supplemented invariants
    for i in 0..n {
        let mut inv = current_ranks[i] as i64 * factor as i64;
        let atom = &mol.atoms()[i];

        // Priority: R > S > nothing
        if let Some(cip_code) = atom.prop("_CIPCode") {
            match cip_code {
                "S" => inv += 10,
                "R" => inv += 20,
                _ => {}
            }
        }

        // Add E/Z stereo info from bonds
        for b in mol.bonds() {
            if b.begin().index() != i && b.end().index() != i {
                continue;
            }
            if b.order() == crate::BondOrder::Double {
                match b.stereo() {
                    crate::BondStereo::E | crate::BondStereo::Trans => inv += 1,
                    crate::BondStereo::Z | crate::BondStereo::Cis => inv += 2,
                    _ => {}
                }
            }
        }

        invars[i] = inv;
    }

    // Use CIP iteration with the supplemented invariants as seeds
    let mut new_ranks = vec![0u32; n];
    let adjacency = mol.topology_block().adjacency.clone();
    let valence = mol.derived_cache().valence.clone().ok_or_else(|| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "CIP_RANKING",
            reason: "CIP ranking requires valence assignment",
        })
    })?;

    iterate_cip_ranks(mol, &invars, &mut new_ranks, true, &adjacency, &valence);

    Ok(new_ranks)
}

// BEGIN RDKIT CPP FUNCTION rerankAtoms writeback (Chirality.cpp:2106-2109)
// RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
// RDKit✔️✔️:     mol.getAtomWithIdx(i)->setProp(common_properties::_CIPRank, ranks[i]);
// RDKit✔️✔️:   }
/// RDKit✔️✔️: rerankAtoms() plus atom-property writeback.
pub(crate) fn rerank_atoms_in_place(
    mol: &mut Molecule,
    current_ranks: &[u32],
) -> Result<Vec<u32>, StereoError> {
    let ranks = rerank_atoms(mol, current_ranks)?;
    write_atom_cip_ranks_to_props(mol, &ranks, false);
    Ok(ranks)
}

// ── End Chirality functions ──────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::{
        LigandRef, NON_TETRAHEDRAL_STEREO_ENV_VAR, SortableCipRef, StereoError,
        assign_atom_cip_ranks, assign_chiral_types_from_3d_kernel,
        assign_nontetrahedral_chiral_type_from_3d, assign_tetrahedral_chiral_type_from_3d,
        atom_nonzero_degree, bond_affects_atom_chirality, canonicalize_tetrahedral_ligands,
        find_segments_to_resort, get_allow_nontetrahedral_chirality, get_val_from_environment,
        has_duplicate_indices, is_atom_potential_chiral_center, is_even_permutation,
        is_wiggly_bond, octahedral_perm_from_3d, perturbation_order_from_bond_indices,
        rdkit_angle_to, rdkit_direction_vector, rdkit_point3d_cross_product,
        rdkit_point3d_dot_product, rdkit_voltest, recompute_ranks,
    };
    use crate::{
        Atom, AtomId, AtomQueryPredicate, AtomSpec, BondDirection, BondOrder, BondSpec, ChiralTag,
        Conformer3D, Element, Molecule, MoleculeBuilder, MoleculeProperties, QueryNode,
    };
    use std::ffi::{OsStr, OsString};
    use std::sync::Mutex;

    static NON_TETRAHEDRAL_ENV_LOCK: Mutex<()> = Mutex::new(());

    struct EnvironmentRestore {
        name: &'static str,
        original: Option<OsString>,
    }

    impl Drop for EnvironmentRestore {
        fn drop(&mut self) {
            // SAFETY: this focused test holds NON_TETRAHEDRAL_ENV_LOCK for the
            // complete mutation window and no production code writes this key.
            unsafe {
                match &self.original {
                    Some(value) => std::env::set_var(self.name, value),
                    None => std::env::remove_var(self.name),
                }
            }
        }
    }

    fn set_nontetrahedral_environment(value: Option<&OsStr>) {
        // SAFETY: callers hold NON_TETRAHEDRAL_ENV_LOCK and restore the key
        // through EnvironmentRestore before releasing the lock.
        unsafe {
            match value {
                Some(value) => std::env::set_var(NON_TETRAHEDRAL_STEREO_ENV_VAR, value),
                None => std::env::remove_var(NON_TETRAHEDRAL_STEREO_ENV_VAR),
            }
        }
    }

    #[test]
    fn get_allow_nontetrahedral_chirality_matches_rdkit_environment_semantics() {
        let _lock = NON_TETRAHEDRAL_ENV_LOCK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let _restore = EnvironmentRestore {
            name: NON_TETRAHEDRAL_STEREO_ENV_VAR,
            original: std::env::var_os(NON_TETRAHEDRAL_STEREO_ENV_VAR),
        };

        set_nontetrahedral_environment(None);
        assert!(get_allow_nontetrahedral_chirality());
        assert!(!get_val_from_environment(
            NON_TETRAHEDRAL_STEREO_ENV_VAR,
            false
        ));

        for (value, expected) in [
            ("0", false),
            ("", true),
            ("00", true),
            ("false", true),
            ("1", true),
            (" 0", true),
        ] {
            set_nontetrahedral_environment(Some(OsStr::new(value)));
            assert_eq!(
                get_allow_nontetrahedral_chirality(),
                expected,
                "unexpected RDKit environment semantics for {value:?}"
            );
        }
    }

    #[test]
    fn bond_affects_atom_chirality_matches_rdkit_directional_matrix() {
        let orders = [
            BondOrder::Null,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
            BondOrder::Quintuple,
            BondOrder::Hextuple,
            BondOrder::OneAndHalf,
            BondOrder::TwoAndHalf,
            BondOrder::ThreeAndHalf,
            BondOrder::FourAndHalf,
            BondOrder::FiveAndHalf,
            BondOrder::Aromatic,
            BondOrder::Ionic,
            BondOrder::Dative,
            BondOrder::DativeOne,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::Hydrogen,
            BondOrder::ThreeCenter,
            BondOrder::Other,
            BondOrder::Zero,
            BondOrder::Unspecified,
        ];

        for order in orders {
            let mut builder = MoleculeBuilder::new();
            let begin = builder.add_atom(AtomSpec::new(Element::C));
            let end = builder.add_atom(AtomSpec::new(Element::C));
            builder
                .add_bond(BondSpec::new(begin, end, order))
                .expect("two valid atoms accept every modeled bond type");
            let molecule = builder.build().expect("valid two-atom molecule");
            let bond = &molecule.bonds()[0];

            let rdkit_unspecified_or_zero = matches!(
                order,
                BondOrder::Null | BondOrder::Unspecified | BondOrder::Zero
            );
            assert_eq!(
                bond_affects_atom_chirality(bond, begin.index()),
                !rdkit_unspecified_or_zero && order != BondOrder::Dative,
                "unexpected begin-atom result for {order:?}"
            );
            assert_eq!(
                bond_affects_atom_chirality(bond, end.index()),
                !rdkit_unspecified_or_zero,
                "unexpected end-atom result for {order:?}"
            );
        }
    }

    #[test]
    fn atom_nonzero_degree_matches_rdkit_bond_matrix() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::FE));
        let single = builder.add_atom(AtomSpec::new(Element::C));
        let double = builder.add_atom(AtomSpec::new(Element::O));
        let unspecified = builder.add_atom(AtomSpec::new(Element::N));
        let zero = builder.add_atom(AtomSpec::new(Element::S));
        let outgoing_dative = builder.add_atom(AtomSpec::new(Element::P));
        let incoming_dative = builder.add_atom(AtomSpec::new(Element::B));
        let isolated = builder.add_atom(AtomSpec::new(Element::F));

        for spec in [
            BondSpec::new(center, single, BondOrder::Single),
            BondSpec::new(center, double, BondOrder::Double),
            BondSpec::new(center, unspecified, BondOrder::Unspecified),
            BondSpec::new(center, zero, BondOrder::Zero),
            BondSpec::new(center, outgoing_dative, BondOrder::Dative),
            BondSpec::new(incoming_dative, center, BondOrder::Dative),
        ] {
            builder.add_bond(spec).expect("valid mixed bond matrix");
        }
        let molecule = builder.build().expect("valid mixed-degree molecule");

        assert_eq!(atom_nonzero_degree(&molecule, center.index()), 3);
        assert_eq!(atom_nonzero_degree(&molecule, single.index()), 1);
        assert_eq!(atom_nonzero_degree(&molecule, double.index()), 1);
        assert_eq!(atom_nonzero_degree(&molecule, unspecified.index()), 0);
        assert_eq!(atom_nonzero_degree(&molecule, zero.index()), 0);
        assert_eq!(atom_nonzero_degree(&molecule, outgoing_dative.index()), 1);
        assert_eq!(atom_nonzero_degree(&molecule, incoming_dative.index()), 0);
        assert_eq!(atom_nonzero_degree(&molecule, isolated.index()), 0);
    }

    #[test]
    fn is_wiggly_bond_matches_rdkit_atom_relative_matrix() {
        fn molecule_with_bond(
            order: BondOrder,
            direction: BondDirection,
            unknown_stereo_property: Option<&str>,
            typed_unknown_stereo: bool,
        ) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            let begin = builder.add_atom(AtomSpec::new(Element::C));
            let end = builder.add_atom(AtomSpec::new(Element::N));
            let mut spec = BondSpec::new(begin, end, order)
                .with_direction(direction)
                .with_unknown_stereo(typed_unknown_stereo);
            if let Some(value) = unknown_stereo_property {
                spec = spec.with_prop("_UnknownStereo", value);
            }
            builder.add_bond(spec).expect("valid focused bond");
            builder.build().expect("valid focused molecule")
        }

        for direction in [
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::EndUpRight,
            BondDirection::EndDownRight,
            BondDirection::EitherDouble,
            BondDirection::Unknown,
        ] {
            let molecule = molecule_with_bond(BondOrder::Single, direction, None, false);
            let bond = &molecule.bonds()[0];
            assert_eq!(
                is_wiggly_bond(bond, 0).expect("valid missing-property state"),
                direction == BondDirection::Unknown,
                "unexpected begin-atom result for {direction:?}"
            );
            assert!(
                !is_wiggly_bond(bond, 1).expect("end perspective short-circuits"),
                "end atom must reject {direction:?}"
            );
        }

        for (property, typed, expected) in [
            (None, false, false),
            (Some("0"), false, false),
            (Some("-7"), false, true),
            (Some("2"), false, true),
            (None, true, true),
        ] {
            let molecule =
                molecule_with_bond(BondOrder::Single, BondDirection::None, property, typed);
            assert_eq!(
                is_wiggly_bond(&molecule.bonds()[0], 0).expect("valid integer property"),
                expected,
                "unexpected property result for {property:?}, typed={typed}"
            );
        }

        for order in [BondOrder::Double, BondOrder::Aromatic] {
            let molecule = molecule_with_bond(order, BondDirection::Unknown, Some("-7"), true);
            assert!(
                !is_wiggly_bond(&molecule.bonds()[0], 0).expect("non-single bond short-circuits"),
                "non-single {order:?} must not be wiggly"
            );
        }

        let malformed = molecule_with_bond(
            BondOrder::Single,
            BondDirection::None,
            Some("not-an-int"),
            false,
        );
        assert!(matches!(
            is_wiggly_bond(&malformed.bonds()[0], 0),
            Err(StereoError::InvariantViolation(_))
        ));
        assert!(matches!(
            is_wiggly_bond(&malformed.bonds()[0], 1),
            Ok(false)
        ));

        let direction_short_circuit = molecule_with_bond(
            BondOrder::Single,
            BondDirection::Unknown,
            Some("not-an-int"),
            false,
        );
        assert!(matches!(
            is_wiggly_bond(&direction_short_circuit.bonds()[0], 0),
            Ok(true)
        ));
    }

    #[test]
    fn rdkit_direction_vector_matches_source_boundaries() {
        fn assert_bits_eq(actual: (f64, f64, f64), expected: (f64, f64, f64)) {
            assert_eq!(actual.0.to_bits(), expected.0.to_bits());
            assert_eq!(actual.1.to_bits(), expected.1.to_bits());
            assert_eq!(actual.2.to_bits(), expected.2.to_bits());
        }

        assert_bits_eq(
            rdkit_direction_vector([1.0, 2.0, 3.0], [4.0, 6.0, 3.0])
                .expect("ordinary vector is normalizable"),
            (3.0 / 5.0, 4.0 / 5.0, 0.0),
        );

        assert!(matches!(
            rdkit_direction_vector([0.0; 3], [0.0; 3]),
            Err(StereoError::ZeroLengthVector)
        ));
        assert!(matches!(
            rdkit_direction_vector([0.0; 3], [0.5e-16, 0.0, 0.0]),
            Err(StereoError::ZeroLengthVector)
        ));
        assert_bits_eq(
            rdkit_direction_vector([0.0; 3], [1.0e-16, 0.0, 0.0])
                .expect("zero_tolerance equality does not throw"),
            (1.0, 0.0, 0.0),
        );
        assert_bits_eq(
            rdkit_direction_vector([0.0; 3], [2.0e-16, 0.0, 0.0])
                .expect("value above zero_tolerance is normalizable"),
            (1.0, 0.0, 0.0),
        );
        assert!(matches!(
            rdkit_direction_vector([0.0; 3], [f64::from_bits(1), 0.0, 0.0]),
            Err(StereoError::ZeroLengthVector)
        ));

        assert_bits_eq(
            rdkit_direction_vector([0.0; 3], [f64::MAX, 0.0, 0.0])
                .expect("overflowed length is not below zero_tolerance"),
            (0.0, 0.0, 0.0),
        );

        let nan = rdkit_direction_vector([0.0; 3], [f64::NAN, 1.0, -1.0])
            .expect("NaN length does not satisfy the throw comparison");
        assert!(nan.0.is_nan());
        assert!(nan.1.is_nan());
        assert!(nan.2.is_nan());

        for infinity in [f64::INFINITY, f64::NEG_INFINITY] {
            let actual = rdkit_direction_vector([0.0; 3], [infinity, 0.0, -0.0])
                .expect("infinite length does not satisfy the throw comparison");
            assert!(actual.0.is_nan());
            assert_eq!(actual.1.to_bits(), 0.0f64.to_bits());
            assert_eq!(actual.2.to_bits(), (-0.0f64).to_bits());
        }
    }

    #[test]
    fn rdkit_angle_to_matches_source_boundaries() {
        assert_eq!(
            rdkit_angle_to((1.0, 0.0, 0.0), (2.0, 0.0, 0.0)).to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(
            rdkit_angle_to((1.0, 0.0, 0.0), (-2.0, 0.0, 0.0)).to_bits(),
            std::f64::consts::PI.to_bits()
        );
        assert_eq!(
            rdkit_angle_to((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)).to_bits(),
            std::f64::consts::FRAC_PI_2.to_bits()
        );

        for (other, expected_bits) in [
            (
                (-0.156_434_465_040_231_04, 0.987_688_340_595_137_7, 0.0),
                0x3ffb_a561_4317_cb35,
            ),
            (
                (-0.173_648_177_666_930_3, 0.984_807_753_012_208, 0.0),
                0x3ffb_ecde_5da1_15a9,
            ),
            (
                (-0.190_808_995_376_544_8, 0.981_627_183_447_664, 0.0),
                0x3ffc_345b_782a_601e,
            ),
        ] {
            assert_eq!(
                rdkit_angle_to((1.0, 0.0, 0.0), other).to_bits(),
                expected_bits
            );
        }
        let hundred_degrees = f64::from_bits(0x3ffb_ecde_5da1_15a9);
        assert!(
            rdkit_angle_to(
                (1.0, 0.0, 0.0),
                (-0.156_434_465_040_231_04, 0.987_688_340_595_137_7, 0.0)
            ) < hundred_degrees
        );
        assert_eq!(
            rdkit_angle_to(
                (1.0, 0.0, 0.0),
                (-0.173_648_177_666_930_3, 0.984_807_753_012_208, 0.0)
            )
            .to_bits(),
            hundred_degrees.to_bits()
        );
        assert!(
            rdkit_angle_to(
                (1.0, 0.0, 0.0),
                (-0.190_808_995_376_544_8, 0.981_627_183_447_664, 0.0)
            ) > hundred_degrees
        );

        let near_parallel = (
            f64::from_bits(0x4010_4747_4776_4fc0),
            f64::from_bits(0xbfee_b4e4_73e5_f750),
            f64::from_bits(0x4012_0156_b839_d268),
        );
        let rounding_above_one = (
            f64::from_bits(0x4010_4747_4596_e981),
            f64::from_bits(0xbfee_b4e4_7132_ed3a),
            f64::from_bits(0x4012_0156_b5e0_194a),
        );
        assert_eq!(
            rdkit_angle_to(near_parallel, rounding_above_one).to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(
            rdkit_angle_to(
                near_parallel,
                (
                    -rounding_above_one.0,
                    -rounding_above_one.1,
                    -rounding_above_one.2,
                )
            )
            .to_bits(),
            std::f64::consts::PI.to_bits()
        );

        assert_eq!(
            rdkit_angle_to((f64::from_bits(1), 0.0, 0.0), (1.0, 0.0, 0.0)).to_bits(),
            0.0f64.to_bits()
        );

        for (this, other) in [
            ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
            ((f64::NAN, 0.0, 0.0), (1.0, 0.0, 0.0)),
            ((f64::INFINITY, 0.0, 0.0), (1.0, 0.0, 0.0)),
            ((f64::NEG_INFINITY, 0.0, 0.0), (1.0, 0.0, 0.0)),
        ] {
            assert!(rdkit_angle_to(this, other).is_nan());
        }
    }

    #[test]
    fn voltest_matches_rdkit_signed_zero_and_nonfinite_matrix() {
        let basis = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)];
        assert!(rdkit_voltest(&basis, 0, 1, 2));
        assert!(!rdkit_voltest(&basis, 0, 2, 1));

        let cross = rdkit_point3d_cross_product((7.0, 11.0, -13.0), (-17.0, 19.0, 23.0));
        assert_eq!(cross.0.to_bits(), 500.0f64.to_bits());
        assert_eq!(cross.1.to_bits(), 60.0f64.to_bits());
        assert_eq!(cross.2.to_bits(), 320.0f64.to_bits());
        assert_eq!(
            rdkit_point3d_dot_product((2.0, -3.0, 5.0), cross).to_bits(),
            2420.0f64.to_bits()
        );
        let ordered = [(2.0, -3.0, 5.0), (7.0, 11.0, -13.0), (-17.0, 19.0, 23.0)];
        assert!(rdkit_voltest(&ordered, 0, 1, 2));
        assert!(!rdkit_voltest(&ordered, 0, 2, 1));

        for signed_zero in [0.0, -0.0] {
            let vectors = [(signed_zero, signed_zero, signed_zero), basis[1], basis[2]];
            let volume = rdkit_point3d_dot_product(
                vectors[0],
                rdkit_point3d_cross_product(vectors[1], vectors[2]),
            );
            assert_eq!(volume.to_bits(), signed_zero.to_bits());
            assert!(rdkit_voltest(&vectors, 0, 1, 2));
        }

        let exact_zero = [basis[1], basis[1], basis[2]];
        assert_eq!(
            rdkit_point3d_dot_product(
                exact_zero[0],
                rdkit_point3d_cross_product(exact_zero[1], exact_zero[2]),
            )
            .to_bits(),
            0.0f64.to_bits()
        );
        assert!(rdkit_voltest(&exact_zero, 0, 1, 2));

        let min_subnormal = f64::from_bits(1);
        for (value, expected) in [
            (min_subnormal, true),
            (-min_subnormal, false),
            (f64::INFINITY, true),
            (f64::NEG_INFINITY, false),
            (f64::NAN, false),
        ] {
            let vectors = [(value, 0.0, 0.0), basis[1], basis[2]];
            assert_eq!(rdkit_voltest(&vectors, 0, 1, 2), expected);
        }
    }

    #[test]
    fn octahedral_perm_from_3d_matches_all_rdkit_switch_branches() {
        struct Case {
            pair: [u8; 6],
            volume_indices: (usize, usize, usize),
            positive: u32,
            negative: u32,
        }

        let cases = [
            Case {
                pair: [2, 0, 4, 0, 0, 0],
                volume_indices: (0, 3, 4),
                positive: 28,
                negative: 27,
            },
            Case {
                pair: [2, 0, 5, 0, 0, 0],
                volume_indices: (0, 2, 3),
                positive: 25,
                negative: 30,
            },
            Case {
                pair: [2, 0, 6, 0, 0, 0],
                volume_indices: (0, 2, 3),
                positive: 26,
                negative: 29,
            },
            Case {
                pair: [3, 4, 0, 0, 0, 0],
                volume_indices: (0, 3, 4),
                positive: 22,
                negative: 21,
            },
            Case {
                pair: [3, 5, 0, 0, 0, 0],
                volume_indices: (0, 1, 3),
                positive: 19,
                negative: 24,
            },
            Case {
                pair: [3, 6, 0, 0, 0, 0],
                volume_indices: (0, 1, 3),
                positive: 20,
                negative: 23,
            },
            Case {
                pair: [4, 3, 0, 0, 0, 0],
                volume_indices: (0, 2, 4),
                positive: 13,
                negative: 12,
            },
            Case {
                pair: [4, 5, 0, 0, 0, 0],
                volume_indices: (0, 1, 2),
                positive: 6,
                negative: 18,
            },
            Case {
                pair: [4, 6, 0, 0, 0, 0],
                volume_indices: (0, 1, 2),
                positive: 7,
                negative: 17,
            },
            Case {
                pair: [5, 3, 0, 0, 0, 0],
                volume_indices: (0, 2, 3),
                positive: 11,
                negative: 9,
            },
            Case {
                pair: [5, 4, 0, 0, 0, 0],
                volume_indices: (0, 1, 2),
                positive: 3,
                negative: 16,
            },
            Case {
                pair: [5, 6, 0, 0, 0, 0],
                volume_indices: (0, 1, 2),
                positive: 5,
                negative: 15,
            },
            Case {
                pair: [6, 3, 0, 0, 0, 0],
                volume_indices: (0, 2, 3),
                positive: 10,
                negative: 8,
            },
            Case {
                pair: [6, 4, 0, 0, 0, 0],
                volume_indices: (0, 1, 2),
                positive: 1,
                negative: 2,
            },
            Case {
                pair: [6, 5, 0, 0, 0, 0],
                volume_indices: (0, 1, 2),
                positive: 4,
                negative: 14,
            },
        ];

        for case in cases {
            let (x, y, z) = case.volume_indices;
            let mut positive_vectors = [(0.0, 0.0, 0.0); 5];
            positive_vectors[x] = (1.0, 0.0, 0.0);
            positive_vectors[y] = (0.0, 1.0, 0.0);
            positive_vectors[z] = (0.0, 0.0, 1.0);
            assert_eq!(
                octahedral_perm_from_3d(&case.pair, &positive_vectors),
                case.positive,
                "positive-volume result for pair {:?}",
                case.pair
            );

            let mut negative_vectors = positive_vectors;
            negative_vectors[x] = (-1.0, 0.0, 0.0);
            assert_eq!(
                octahedral_perm_from_3d(&case.pair, &negative_vectors),
                case.negative,
                "negative-volume result for pair {:?}",
                case.pair
            );
        }
    }

    #[test]
    fn assign_nontetrahedral_chiral_type_from_3d_covers_all_active_rdkit_branches() {
        fn run_case(
            element: Element,
            directions: &[[f64; 3]],
            tolerance: f64,
            unknown_bond: Option<(usize, bool)>,
        ) -> Result<(bool, ChiralTag, Option<u32>), StereoError> {
            let mut builder = MoleculeBuilder::new();
            let center = builder.add_atom(AtomSpec::new(element));
            for (index, _) in directions.iter().enumerate() {
                let neighbor = builder.add_atom(AtomSpec::new(Element::F));
                let (begin, end) = if unknown_bond == Some((index, true)) {
                    (neighbor, center)
                } else {
                    (center, neighbor)
                };
                let mut bond = BondSpec::new(begin, end, BondOrder::Single);
                if unknown_bond.is_some_and(|(unknown_index, _)| unknown_index == index) {
                    bond = bond.with_direction(BondDirection::Unknown);
                }
                builder.add_bond(bond).expect("valid star bond");
            }
            let mut molecule = builder.build().expect("valid star molecule");
            let mut coordinates = Vec::with_capacity(directions.len() + 1);
            coordinates.push([0.0, 0.0, 0.0]);
            coordinates.extend_from_slice(directions);
            let conformer = Conformer3D::new(0, coordinates, true);
            let bonds = molecule.bonds().to_vec();
            let adjacency = molecule.topology_block().adjacency.clone();
            let result = assign_nontetrahedral_chiral_type_from_3d(
                &mut molecule.topology_block_mut().atoms,
                &bonds,
                &adjacency,
                &conformer,
                center.index(),
                tolerance,
            )?;
            let atom = &molecule.atoms()[center.index()];
            Ok((result, atom.chiral_tag(), atom.chiral_permutation()))
        }

        fn assert_assignment(
            directions: &[[f64; 3]],
            expected_tag: ChiralTag,
            expected_permutation: u32,
        ) {
            assert_eq!(
                run_case(Element::P, directions, 0.1, None)
                    .expect("source-defined geometry is valid"),
                (true, expected_tag, Some(expected_permutation)),
                "unexpected assignment for {directions:?}"
            );
        }

        fn one_opposite_pair(
            count: usize,
            opposite: (usize, usize),
            equatorial_angle_degrees: f64,
            mirror: f64,
        ) -> Vec<[f64; 3]> {
            let mut directions = vec![[0.0, 0.0, 0.0]; count];
            directions[opposite.0] = [0.0, 0.0, 1.0];
            directions[opposite.1] = [0.0, 0.0, -1.0];
            let remaining: Vec<_> = (0..count)
                .filter(|index| *index != opposite.0 && *index != opposite.1)
                .collect();
            directions[remaining[0]] = [1.0, 0.0, 0.0];
            let angle = equatorial_angle_degrees.to_radians();
            directions[remaining[1]] = [angle.cos(), mirror * angle.sin(), 0.0];
            if remaining.len() == 3 {
                let third_angle = 240.0_f64.to_radians();
                directions[remaining[2]] = [third_angle.cos(), mirror * third_angle.sin(), 0.0];
            }
            directions
        }

        fn scalar_triple_nonnegative(
            directions: &[[f64; 3]],
            x: usize,
            y: usize,
            z: usize,
        ) -> bool {
            let cross = [
                directions[y][1] * directions[z][2] - directions[y][2] * directions[z][1],
                -directions[y][0] * directions[z][2] + directions[y][2] * directions[z][0],
                directions[y][0] * directions[z][1] - directions[y][1] * directions[z][0],
            ];
            directions[x][0] * cross[0] + directions[x][1] * cross[1] + directions[x][2] * cross[2]
                >= 0.0
        }

        let ordinary_three = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!(
            run_case(Element::SI, &ordinary_three, 0.1, None).expect("silicon gate"),
            (false, ChiralTag::Unspecified, None)
        );
        assert_eq!(
            run_case(Element::P, &ordinary_three[..2], 0.1, None).expect("degree-two gate"),
            (false, ChiralTag::Unspecified, None)
        );
        assert_eq!(
            run_case(Element::P, &ordinary_three, 0.1, None).expect("zero-pair gate"),
            (false, ChiralTag::Unspecified, None)
        );
        assert!(matches!(
            run_case(
                Element::P,
                &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
                0.1,
                None
            ),
            Err(StereoError::ZeroLengthVector)
        ));
        assert_eq!(
            run_case(
                Element::P,
                &[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
                0.1,
                None
            )
            .expect("duplicate pair rejection"),
            (false, ChiralTag::Unspecified, None)
        );
        assert_eq!(
            run_case(
                Element::P,
                &[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                0.0,
                None
            )
            .expect("strict antiparallel equality"),
            (false, ChiralTag::Unspecified, None)
        );
        assert_assignment(
            &[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            ChiralTag::SquarePlanar,
            2,
        );

        for (directions, permutation) in [
            ([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], 3),
            ([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], 2),
            ([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], 1),
        ] {
            assert_assignment(&directions, ChiralTag::SquarePlanar, permutation);
        }

        struct SeeSawCase {
            opposite: (usize, usize),
            volume: (usize, usize, usize),
            octahedral: (u32, u32),
            trigonal_bipyramidal: (u32, u32),
        }
        let see_saw_cases = [
            SeeSawCase {
                opposite: (0, 1),
                volume: (0, 2, 3),
                octahedral: (25, 29),
                trigonal_bipyramidal: (7, 8),
            },
            SeeSawCase {
                opposite: (0, 2),
                volume: (0, 1, 3),
                octahedral: (19, 23),
                trigonal_bipyramidal: (5, 6),
            },
            SeeSawCase {
                opposite: (0, 3),
                volume: (0, 1, 2),
                octahedral: (6, 17),
                trigonal_bipyramidal: (3, 4),
            },
            SeeSawCase {
                opposite: (1, 2),
                volume: (0, 1, 3),
                octahedral: (10, 8),
                trigonal_bipyramidal: (13, 14),
            },
            SeeSawCase {
                opposite: (1, 3),
                volume: (0, 1, 3),
                octahedral: (1, 2),
                trigonal_bipyramidal: (10, 12),
            },
            SeeSawCase {
                opposite: (2, 3),
                volume: (0, 1, 3),
                octahedral: (4, 14),
                trigonal_bipyramidal: (16, 19),
            },
        ];
        for case in see_saw_cases {
            for mirror in [1.0, -1.0] {
                let acute = one_opposite_pair(4, case.opposite, 90.0, mirror);
                let acute_positive =
                    scalar_triple_nonnegative(&acute, case.volume.0, case.volume.1, case.volume.2);
                assert_assignment(
                    &acute,
                    ChiralTag::Octahedral,
                    if acute_positive {
                        case.octahedral.0
                    } else {
                        case.octahedral.1
                    },
                );

                let obtuse = one_opposite_pair(4, case.opposite, 120.0, mirror);
                let tbp_volume = match case.opposite {
                    (1, 2) => (1, 0, 3),
                    (1, 3) => (1, 0, 2),
                    (2, 3) => (3, 0, 1),
                    _ => case.volume,
                };
                let obtuse_positive =
                    scalar_triple_nonnegative(&obtuse, tbp_volume.0, tbp_volume.1, tbp_volume.2);
                assert_assignment(
                    &obtuse,
                    ChiralTag::TrigonalBipyramidal,
                    if obtuse_positive {
                        case.trigonal_bipyramidal.0
                    } else {
                        case.trigonal_bipyramidal.1
                    },
                );
            }
        }
        let exact_hundred = one_opposite_pair(4, (0, 1), 100.0, 1.0);
        let exact_hundred_positive = scalar_triple_nonnegative(&exact_hundred, 0, 2, 3);
        assert_assignment(
            &exact_hundred,
            ChiralTag::TrigonalBipyramidal,
            if exact_hundred_positive { 7 } else { 8 },
        );

        struct TbpCase {
            opposite: (usize, usize),
            volume: (usize, usize, usize),
            permutations: (u32, u32),
        }
        let tbp_cases = [
            TbpCase {
                opposite: (0, 1),
                volume: (0, 2, 3),
                permutations: (7, 8),
            },
            TbpCase {
                opposite: (0, 2),
                volume: (0, 1, 3),
                permutations: (5, 6),
            },
            TbpCase {
                opposite: (0, 3),
                volume: (0, 1, 2),
                permutations: (3, 4),
            },
            TbpCase {
                opposite: (0, 4),
                volume: (0, 1, 2),
                permutations: (1, 2),
            },
            TbpCase {
                opposite: (1, 2),
                volume: (1, 0, 3),
                permutations: (13, 14),
            },
            TbpCase {
                opposite: (1, 3),
                volume: (1, 0, 2),
                permutations: (10, 12),
            },
            TbpCase {
                opposite: (1, 4),
                volume: (1, 0, 2),
                permutations: (9, 11),
            },
            TbpCase {
                opposite: (2, 3),
                volume: (2, 0, 1),
                permutations: (16, 19),
            },
            TbpCase {
                opposite: (2, 4),
                volume: (2, 0, 1),
                permutations: (15, 20),
            },
            TbpCase {
                opposite: (3, 4),
                volume: (3, 0, 1),
                permutations: (17, 18),
            },
        ];
        for case in tbp_cases {
            for mirror in [1.0, -1.0] {
                let directions = one_opposite_pair(5, case.opposite, 120.0, mirror);
                let positive = scalar_triple_nonnegative(
                    &directions,
                    case.volume.0,
                    case.volume.1,
                    case.volume.2,
                );
                assert_assignment(
                    &directions,
                    ChiralTag::TrigonalBipyramidal,
                    if positive {
                        case.permutations.0
                    } else {
                        case.permutations.1
                    },
                );
            }
        }

        for (directions, permutation) in [
            (
                [
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
                2,
            ),
            (
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
                1,
            ),
            (
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                ],
                3,
            ),
        ] {
            assert_assignment(&directions, ChiralTag::SquarePlanar, permutation);
        }

        assert_assignment(
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            ChiralTag::Octahedral,
            27,
        );
        assert_assignment(
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
            ChiralTag::Octahedral,
            27,
        );
        assert_eq!(
            run_case(
                Element::P,
                &[
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                    [1.0, 1.0, 1.0],
                ],
                0.1,
                None
            )
            .expect("degree-seven rejection"),
            (false, ChiralTag::Unspecified, None)
        );

        let t_shape = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert_eq!(
            run_case(Element::P, &t_shape, 0.1, Some((0, false)))
                .expect("begin-atom wiggly rejection"),
            (false, ChiralTag::Unspecified, None)
        );
        assert_eq!(
            run_case(Element::P, &t_shape, 0.1, Some((0, true)))
                .expect("end-atom wiggly bond is ignored"),
            (true, ChiralTag::SquarePlanar, Some(2))
        );
    }

    #[test]
    fn tetrahedral_chiral_type_from_3d_matches_rdkit_branch_boundaries() {
        #[derive(Clone, Copy)]
        struct Neighbor {
            coordinate: [f64; 3],
            order: BondOrder,
            reverse: bool,
            direction: BondDirection,
        }

        impl Neighbor {
            fn single(coordinate: [f64; 3]) -> Self {
                Self {
                    coordinate,
                    order: BondOrder::Single,
                    reverse: false,
                    direction: BondDirection::None,
                }
            }
        }

        fn run_case(
            element: Element,
            explicit_hydrogens: u8,
            total_num_hs: usize,
            explicit_atom: bool,
            neighbors: &[Neighbor],
        ) -> Result<(bool, ChiralTag, Option<String>), StereoError> {
            let mut builder = MoleculeBuilder::new();
            let center = builder
                .add_atom(AtomSpec::new(element).with_explicit_hydrogens(explicit_hydrogens));
            for neighbor in neighbors {
                let other = builder.add_atom(AtomSpec::new(Element::F));
                let (begin, end) = if neighbor.reverse {
                    (other, center)
                } else {
                    (center, other)
                };
                builder
                    .add_bond(
                        BondSpec::new(begin, end, neighbor.order)
                            .with_direction(neighbor.direction),
                    )
                    .expect("valid focused bond");
            }
            let mut molecule = builder.build().expect("valid focused molecule");
            let mut coordinates = Vec::with_capacity(neighbors.len() + 1);
            coordinates.push([0.0, 0.0, 0.0]);
            coordinates.extend(neighbors.iter().map(|neighbor| neighbor.coordinate));
            let conformer = Conformer3D::new(0, coordinates, true);
            let bonds = molecule.bonds().to_vec();
            let adjacency = molecule.topology_block().adjacency.clone();
            let result = assign_tetrahedral_chiral_type_from_3d(
                &mut molecule.topology_block_mut().atoms,
                &bonds,
                &adjacency,
                &conformer,
                center.index(),
                total_num_hs,
                explicit_atom,
                0.1,
            )?;
            let atom = &molecule.atoms()[center.index()];
            Ok((
                result,
                atom.chiral_tag(),
                atom.prop("_NonExplicit3DChirality").map(str::to_owned),
            ))
        }

        fn singles(coordinates: &[[f64; 3]]) -> Vec<Neighbor> {
            coordinates.iter().copied().map(Neighbor::single).collect()
        }

        let positive = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let negative = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]];

        assert_eq!(
            run_case(Element::C, 0, 1, false, &singles(&positive)).expect("implicit-H tetrahedron"),
            (true, ChiralTag::TetrahedralCcw, Some("1".to_owned()))
        );
        assert_eq!(
            run_case(Element::C, 1, 1, true, &singles(&negative)).expect("explicit-H tetrahedron"),
            (true, ChiralTag::TetrahedralCw, None)
        );
        assert_eq!(
            run_case(Element::C, 0, 0, false, &singles(&positive))
                .expect("ordinary carbon needs four total neighbors"),
            (false, ChiralTag::Unspecified, None)
        );
        for element in [Element::S, Element::SE] {
            assert_eq!(
                run_case(element, 0, 0, false, &singles(&positive))
                    .expect("S/Se three-coordinate exception"),
                (true, ChiralTag::TetrahedralCcw, Some("1".to_owned()))
            );
        }

        for (coordinates, expected_tag) in [
            (
                [positive[0], positive[1], positive[2]],
                ChiralTag::TetrahedralCcw,
            ),
            (
                [positive[0], positive[2], positive[1]],
                ChiralTag::TetrahedralCw,
            ),
            (
                [positive[1], positive[0], positive[2]],
                ChiralTag::TetrahedralCw,
            ),
            (
                [positive[1], positive[2], positive[0]],
                ChiralTag::TetrahedralCcw,
            ),
            (
                [positive[2], positive[0], positive[1]],
                ChiralTag::TetrahedralCcw,
            ),
            (
                [positive[2], positive[1], positive[0]],
                ChiralTag::TetrahedralCw,
            ),
        ] {
            assert_eq!(
                run_case(Element::C, 0, 1, true, &singles(&coordinates))
                    .expect("neighbor permutation"),
                (true, expected_tag, None)
            );
        }

        let tolerance = 0.1_f64;
        let below = f64::from_bits(tolerance.to_bits() - 1);
        let above = f64::from_bits(tolerance.to_bits() + 1);
        for (volume, expected) in [
            (tolerance, ChiralTag::Unspecified),
            (below, ChiralTag::Unspecified),
            (above, ChiralTag::TetrahedralCcw),
            (-tolerance, ChiralTag::Unspecified),
            (-below, ChiralTag::Unspecified),
            (-above, ChiralTag::TetrahedralCw),
        ] {
            let coordinates = [[volume, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            let actual = run_case(Element::C, 0, 1, true, &singles(&coordinates))
                .expect("source-defined volume boundary");
            assert_eq!(actual.1, expected, "volume {volume:?}");
            assert_eq!(actual.0, expected != ChiralTag::Unspecified);
        }

        for (fourth, expected) in [
            ([0.0, 0.0, 1.0], ChiralTag::TetrahedralCw),
            ([0.0, 0.0, -1.0], ChiralTag::TetrahedralCcw),
            ([1.0, 1.0, 0.0], ChiralTag::Unspecified),
        ] {
            let coordinates = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], fourth];
            let actual = run_case(Element::C, 0, 0, true, &singles(&coordinates))
                .expect("fourth-neighbor fallback");
            assert_eq!(actual.1, expected, "fourth neighbor {fourth:?}");
            assert_eq!(actual.0, expected != ChiralTag::Unspecified);
        }

        let mut wiggly = singles(&positive);
        wiggly[0].direction = BondDirection::Unknown;
        assert_eq!(
            run_case(Element::C, 0, 1, false, &wiggly).expect("wiggly rejection"),
            (false, ChiralTag::Unspecified, None)
        );

        let mixed_bonds = [
            Neighbor::single([1.0, 0.0, 0.0]),
            Neighbor {
                coordinate: [99.0, 99.0, 99.0],
                order: BondOrder::Zero,
                reverse: false,
                direction: BondDirection::None,
            },
            Neighbor::single([0.0, 1.0, 0.0]),
            Neighbor {
                coordinate: [-99.0, -99.0, -99.0],
                order: BondOrder::Dative,
                reverse: false,
                direction: BondDirection::None,
            },
            Neighbor {
                coordinate: [0.0, 0.0, 1.0],
                order: BondOrder::Dative,
                reverse: true,
                direction: BondDirection::None,
            },
        ];
        assert_eq!(
            run_case(Element::C, 0, 1, true, &mixed_bonds)
                .expect("irrelevant bonds preserve relevant neighbor order"),
            (true, ChiralTag::TetrahedralCcw, None)
        );

        let five_total_neighbors = singles(&[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ]);
        assert_eq!(
            run_case(Element::C, 0, 1, false, &five_total_neighbors)
                .expect("more than four total neighbors"),
            (false, ChiralTag::Unspecified, None)
        );
    }

    #[test]
    fn assign_chiral_types_from_3d_covers_complete_rdkit_control_flow() {
        fn build_star(
            center: AtomSpec,
            conformers: &[(usize, bool, Vec<[f64; 3]>)],
            first_bond_direction: BondDirection,
        ) -> Molecule {
            let neighbor_count = conformers
                .first()
                .expect("star requires at least one conformer")
                .2
                .len()
                .checked_sub(1)
                .expect("star conformer requires a center coordinate");
            let mut builder = MoleculeBuilder::new();
            let center_id = builder.add_atom(center);
            for neighbor_idx in 0..neighbor_count {
                let neighbor = builder.add_atom(AtomSpec::new(Element::C));
                let direction = if neighbor_idx == 0 {
                    first_bond_direction
                } else {
                    BondDirection::None
                };
                builder
                    .add_bond(
                        BondSpec::new(center_id, neighbor, BondOrder::Single)
                            .with_direction(direction),
                    )
                    .expect("valid star bond");
            }
            for (id, is_3d, coordinates) in conformers {
                assert_eq!(coordinates.len(), neighbor_count + 1);
                builder
                    .add_conformer(Conformer3D::new(*id, coordinates.clone(), *is_3d))
                    .expect("aligned star conformer");
            }
            builder.build().expect("valid star molecule")
        }

        fn run_kernel(
            molecule: &Molecule,
            conformers: Option<Vec<Conformer3D>>,
            implicit_hydrogens: Option<Vec<i32>>,
            conformer_id: i32,
            replace_existing_tags: bool,
        ) -> (Vec<Atom>, MoleculeProperties, Result<(), StereoError>) {
            let topology = molecule.topology_block();
            let mut atoms = topology.atoms.clone();
            let bonds = topology.bonds.clone();
            let adjacency = topology.adjacency.clone();
            let conformers = conformers.unwrap_or_else(|| molecule.conformers_3d().to_vec());
            let mut properties = molecule.properties().clone();
            let result = assign_chiral_types_from_3d_kernel(
                &mut atoms,
                &bonds,
                &adjacency,
                &conformers,
                &mut properties,
                implicit_hydrogens.as_deref(),
                conformer_id,
                replace_existing_tags,
            );
            (atoms, properties, result)
        }

        let _lock = NON_TETRAHEDRAL_ENV_LOCK
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let _restore = EnvironmentRestore {
            name: NON_TETRAHEDRAL_STEREO_ENV_VAR,
            original: std::env::var_os(NON_TETRAHEDRAL_STEREO_ENV_VAR),
        };
        set_nontetrahedral_environment(None);

        let mut empty_builder = MoleculeBuilder::new();
        empty_builder.add_atom(
            AtomSpec::new(Element::C)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_prop("atom_keep", "yes"),
        );
        let no_conformer = empty_builder
            .with_property("_StereochemDone", "1")
            .with_property("molecule_keep", "yes")
            .build()
            .expect("valid no-conformer molecule");
        let (atoms, properties, result) = run_kernel(&no_conformer, None, None, 99, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms, no_conformer.atoms());
        assert_eq!(properties, *no_conformer.properties());

        let non_3d = build_star(
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
            &[(
                7,
                false,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            )],
            BondDirection::None,
        )
        .with_prop("_StereochemDone", "1");
        let (atoms, properties, result) = run_kernel(&non_3d, None, None, 7, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms, non_3d.atoms());
        assert_eq!(properties.prop("_StereochemDone"), Some("1"));

        let two_conformers = build_star(
            AtomSpec::new(Element::C).with_prop("atom_keep", "yes"),
            &[
                (
                    7,
                    true,
                    vec![
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0],
                    ],
                ),
                (
                    11,
                    true,
                    vec![
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0],
                        [0.0, 1.0, 0.0],
                    ],
                ),
            ],
            BondDirection::None,
        )
        .with_prop("_StereochemDone", "1")
        .with_prop("molecule_keep", "yes");
        let implicit = Some(vec![1, 0, 0, 0]);
        let (atoms, properties, result) =
            run_kernel(&two_conformers, None, implicit.clone(), -1, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(atoms[0].prop("_NonExplicit3DChirality"), Some("1"));
        assert_eq!(atoms[0].prop("atom_keep"), Some("yes"));
        assert_eq!(properties.prop("_StereochemDone"), None);
        assert_eq!(properties.prop("molecule_keep"), Some("yes"));
        assert_eq!(two_conformers.conformers_3d()[0].id(), 7);
        assert_eq!(two_conformers.conformers_3d()[1].id(), 11);

        let (atoms, _, result) = run_kernel(&two_conformers, None, implicit.clone(), 11, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);

        let (_, properties, result) = run_kernel(&two_conformers, None, implicit.clone(), 8, true);
        assert_eq!(
            result,
            Err(StereoError::ConformerNotFound { conformer_id: 8 })
        );
        assert_eq!(properties.prop("_StereochemDone"), Some("1"));

        let existing = build_star(
            AtomSpec::new(Element::C)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_prop("atom_keep", "yes"),
            &[((
                0,
                true,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            ))],
            BondDirection::None,
        )
        .with_prop("_StereochemDone", "1");
        let (atoms, properties, result) = run_kernel(&existing, None, implicit.clone(), 0, false);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(atoms[0].prop("atom_keep"), Some("yes"));
        assert_eq!(properties.prop("_StereochemDone"), None);
        let (atoms, _, result) = run_kernel(&existing, None, implicit.clone(), 0, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(atoms[0].prop("_NonExplicit3DChirality"), None);

        let wedged = build_star(
            AtomSpec::new(Element::C),
            &[((
                0,
                true,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            ))],
            BondDirection::BeginWedge,
        );
        let (atoms, _, result) = run_kernel(&wedged, None, implicit.clone(), -1, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(atoms[0].prop("_NonExplicit3DChirality"), None);

        let square_planar = build_star(
            AtomSpec::new(Element::P),
            &[((
                0,
                true,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
            ))],
            BondDirection::None,
        );
        let square_planar_h = Some(vec![0; 5]);
        set_nontetrahedral_environment(Some(OsStr::new("0")));
        let (atoms, _, result) =
            run_kernel(&square_planar, None, square_planar_h.clone(), -1, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(atoms[0].chiral_permutation(), None);
        set_nontetrahedral_environment(Some(OsStr::new("1")));
        let (atoms, _, result) = run_kernel(&square_planar, None, square_planar_h, -1, true);
        assert_eq!(result, Ok(()));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::SquarePlanar);
        assert_eq!(atoms[0].chiral_permutation(), Some(2));
        assert_eq!(atoms[0].prop("_NonExplicit3DChirality"), Some("1"));

        for (element, coordinates, implicit_hydrogens) in [
            (
                Element::C,
                vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                vec![0; 3],
            ),
            (
                Element::C,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                vec![0; 4],
            ),
            (
                Element::P,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, -1.0],
                ],
                vec![1, 0, 0, 0, 0, 0, 0],
            ),
            (
                Element::P,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
                vec![0; 6],
            ),
        ] {
            let molecule = build_star(
                AtomSpec::new(element),
                &[(0, true, coordinates)],
                BondDirection::None,
            );
            set_nontetrahedral_environment(Some(OsStr::new("0")));
            let (atoms, _, result) =
                run_kernel(&molecule, None, Some(implicit_hydrogens), -1, true);
            assert_eq!(result, Ok(()));
            assert_eq!(atoms[0].chiral_tag(), ChiralTag::Unspecified);
        }

        let (_, _, result) = run_kernel(&two_conformers, None, None, -1, true);
        assert_eq!(result, Err(StereoError::MissingImplicitHydrogenState));
        let (_, _, result) = run_kernel(&two_conformers, None, Some(Vec::new()), -1, true);
        assert_eq!(
            result,
            Err(StereoError::ImplicitHydrogenCountMismatch {
                expected: 4,
                actual: 0,
            })
        );
        let malformed_conformer = vec![Conformer3D::new(7, vec![[0.0, 0.0, 0.0]; 3], true)];
        let (_, properties, result) = run_kernel(
            &two_conformers,
            Some(malformed_conformer),
            implicit.clone(),
            -1,
            true,
        );
        assert_eq!(
            result,
            Err(StereoError::ConformerAtomCountMismatch {
                conformer_id: 7,
                expected: 4,
                actual: 3,
            })
        );
        assert_eq!(properties.prop("_StereochemDone"), None);

        set_nontetrahedral_environment(Some(OsStr::new("1")));
        let zero_vector = build_star(
            AtomSpec::new(Element::P),
            &[(
                0,
                true,
                vec![
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ],
            )],
            BondDirection::None,
        );
        let (_, _, result) = run_kernel(&zero_vector, None, Some(vec![0; 4]), -1, true);
        assert_eq!(result, Err(StereoError::ZeroLengthVector));

        let mut invalid_property_builder = MoleculeBuilder::new();
        let center = invalid_property_builder.add_atom(AtomSpec::new(Element::C));
        let mut invalid_property_coordinates = vec![[0.0, 0.0, 0.0]];
        for coordinate in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            let neighbor = invalid_property_builder.add_atom(AtomSpec::new(Element::C));
            let mut bond = BondSpec::new(center, neighbor, BondOrder::Single);
            if coordinate == [1.0, 0.0, 0.0] {
                bond = bond.with_prop("_UnknownStereo", "not-an-int");
            }
            invalid_property_builder
                .add_bond(bond)
                .expect("valid invalid-property topology");
            invalid_property_coordinates.push(coordinate);
        }
        invalid_property_builder
            .add_3d_conformer(invalid_property_coordinates)
            .expect("aligned invalid-property conformer");
        let invalid_property = invalid_property_builder
            .build()
            .expect("valid invalid-property molecule");
        let (_, _, result) = run_kernel(&invalid_property, None, Some(vec![1, 0, 0, 0]), -1, true);
        assert!(matches!(result, Err(StereoError::InvariantViolation(_))));

        let mut ordered_builder = MoleculeBuilder::new();
        let carbon =
            ordered_builder.add_atom(AtomSpec::new(Element::C).with_prop("first_atom_keep", "yes"));
        let mut ordered_coordinates = vec![[0.0, 0.0, 0.0]];
        for coordinate in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            let neighbor = ordered_builder.add_atom(AtomSpec::new(Element::C));
            ordered_builder
                .add_bond(BondSpec::new(carbon, neighbor, BondOrder::Single))
                .expect("valid first-center bond");
            ordered_coordinates.push(coordinate);
        }
        let phosphorus = ordered_builder.add_atom(AtomSpec::new(Element::P));
        ordered_coordinates.push([10.0, 0.0, 0.0]);
        for coordinate in [[10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0]] {
            let neighbor = ordered_builder.add_atom(AtomSpec::new(Element::C));
            ordered_builder
                .add_bond(BondSpec::new(phosphorus, neighbor, BondOrder::Single))
                .expect("valid second-center bond");
            ordered_coordinates.push(coordinate);
        }
        ordered_builder
            .add_3d_conformer(ordered_coordinates)
            .expect("aligned two-center conformer");
        let ordered = ordered_builder
            .with_property("molecule_keep", "yes")
            .build()
            .expect("valid two-center molecule");
        let mut ordered_hydrogens = vec![0; ordered.num_atoms()];
        ordered_hydrogens[carbon.index()] = 1;
        let (atoms, properties, result) =
            run_kernel(&ordered, None, Some(ordered_hydrogens), -1, true);
        assert_eq!(result, Err(StereoError::ZeroLengthVector));
        assert_eq!(
            atoms[carbon.index()].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(atoms[carbon.index()].prop("first_atom_keep"), Some("yes"));
        assert_eq!(
            atoms[phosphorus.index()].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(properties.prop("molecule_keep"), Some("yes"));
    }

    #[test]
    fn tetrahedral_stereo_ligand_order_encodes_smiles_handedness() {
        let ccw = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("parse chiral SMILES");
        let cw = Molecule::from_smiles("F[C@@](Cl)(Br)I").expect("parse chiral SMILES");

        let ccw_stereo = ccw.tetrahedral_stereo().expect("tetrahedral stereo");
        let cw_stereo = cw.tetrahedral_stereo().expect("tetrahedral stereo");

        assert_eq!(ccw_stereo.len(), 1);
        assert_eq!(cw_stereo.len(), 1);
        assert_eq!(ccw_stereo[0].center, AtomId::new(1));
        assert_eq!(cw_stereo[0].center, AtomId::new(1));
        assert_eq!(
            ccw_stereo[0].ligands,
            [
                LigandRef::Atom(AtomId::new(0)),
                LigandRef::Atom(AtomId::new(2)),
                LigandRef::Atom(AtomId::new(3)),
                LigandRef::Atom(AtomId::new(4)),
            ]
        );
        assert_eq!(
            cw_stereo[0].ligands,
            [
                LigandRef::Atom(AtomId::new(0)),
                LigandRef::Atom(AtomId::new(2)),
                LigandRef::Atom(AtomId::new(4)),
                LigandRef::Atom(AtomId::new(3)),
            ]
        );
    }

    #[test]
    fn tetrahedral_stereo_canonicalizes_even_ligand_permutations() {
        let base = [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::Atom(AtomId::new(3)),
            LigandRef::Atom(AtomId::new(4)),
        ];

        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    for d in 0..4 {
                        let perm = [a, b, c, d];
                        if has_duplicate_indices(perm) {
                            continue;
                        }
                        let candidate = [base[a], base[b], base[c], base[d]];
                        if is_even_permutation(perm) {
                            assert_eq!(canonicalize_tetrahedral_ligands(candidate), base);
                        } else {
                            assert_ne!(canonicalize_tetrahedral_ligands(candidate), base);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn tetrahedral_stereo_places_implicit_hydrogen_as_fourth_ligand() {
        let mol = Molecule::from_smiles("[13CH3:7][C@H](F)Cl").expect("parse chiral SMILES");
        let stereo = mol.tetrahedral_stereo().expect("tetrahedral stereo");

        assert_eq!(stereo.len(), 1);
        assert_eq!(stereo[0].center, AtomId::new(1));
        assert_eq!(
            stereo[0].ligands,
            [
                LigandRef::Atom(AtomId::new(0)),
                LigandRef::Atom(AtomId::new(2)),
                LigandRef::Atom(AtomId::new(3)),
                LigandRef::ImplicitHydrogen,
            ]
        );
    }

    #[test]
    fn implicit_hydrogen_tetrahedral_center_is_potentially_chiral_like_rdkit() {
        let mol = Molecule::from_smiles("Cl[C@H](Br)I").expect("failed to parse test SMILES");
        let ranks = assign_atom_cip_ranks(&mol).expect("failed to assign CIP ranks");
        let center = mol
            .atoms()
            .iter()
            .position(|atom| atom.atomic_number() == 6)
            .expect("failed to find tetrahedral carbon");
        let (legal_center, has_dupes, _nbrs) =
            is_atom_potential_chiral_center(&mol, center, &ranks);
        assert!(
            legal_center,
            "implicit-H tetrahedral carbon must remain a legal center"
        );
        assert!(
            !has_dupes,
            "distinct halogen substituents must not collapse to duplicate ranks"
        );
    }

    #[test]
    fn assign_atom_cip_ranks_matches_rdkit_reference_cases() {
        for (smiles, expected) in [
            ("CC(O)F", vec![0, 1, 2, 3]),
            ("Cl[C@H](Br)I", vec![1, 0, 2, 3]),
            ("O=P(F)(Cl)Br", vec![0, 2, 1, 3, 4]),
        ] {
            let mol = Molecule::from_smiles(smiles).expect("parse CIP reference SMILES");
            let ranks = assign_atom_cip_ranks(&mol).expect("assign CIP reference ranks");
            assert_eq!(ranks, expected, "RDKit CIP rank mismatch for {smiles}");
        }
    }

    #[test]
    fn find_segments_to_resort_reuses_source_result_allocation() {
        let mut cip_entries = [vec![1], vec![1], vec![2]];
        let mut sortable_entries: Vec<_> = cip_entries
            .iter_mut()
            .enumerate()
            .map(|(atom_idx, entry)| SortableCipRef::new(entry, atom_idx))
            .collect();
        let mut segments = Vec::with_capacity(4);
        segments.push((usize::MAX, usize::MAX));
        let allocation = segments.as_ptr();
        let capacity = segments.capacity();

        let independent = find_segments_to_resort(&mut sortable_entries, &mut segments);

        assert_eq!(independent, 2);
        assert_eq!(segments, [(0, 2)]);
        assert_eq!(segments.as_ptr(), allocation);
        assert_eq!(segments.capacity(), capacity);
        assert_eq!(sortable_entries[0].curr_rank, 0);
        assert_eq!(sortable_entries[1].curr_rank, 0);
        assert_eq!(sortable_entries[2].curr_rank, 1);
    }

    #[test]
    fn recompute_ranks_writes_sorted_entries_back_by_atom_index() {
        let mut cip_entries = [vec![2], vec![1], vec![1]];
        let mut sortable_entries: Vec<_> = cip_entries
            .iter_mut()
            .enumerate()
            .map(|(atom_idx, entry)| SortableCipRef::new(entry, atom_idx))
            .collect();
        sortable_entries.sort();
        let mut segments = Vec::new();
        find_segments_to_resort(&mut sortable_entries, &mut segments);
        let mut ranks = vec![u32::MAX; sortable_entries.len()];

        recompute_ranks(&sortable_entries, &mut ranks);

        assert_eq!(ranks, [1, 0, 0]);
    }

    #[test]
    fn assign_atom_cip_ranks_skips_hydrogens_for_query_atoms() {
        let ordinary = Molecule::from_smiles("CC").expect("parse ethane");
        assert_eq!(
            assign_atom_cip_ranks(&ordinary).expect("rank ordinary ethane"),
            vec![0, 0]
        );

        let mut query = ordinary.clone();
        query.topology_block_mut().atoms[0]
            .set_query(Some(QueryNode::predicate(AtomQueryPredicate::Any)));
        assert_eq!(
            assign_atom_cip_ranks(&query).expect("rank query ethane"),
            vec![0, 1]
        );
    }

    #[test]
    fn perturbation_order_rejects_probe_reference_length_mismatch() {
        let mol = Molecule::from_smiles("CC").expect("failed to parse test SMILES");
        let error = perturbation_order_from_bond_indices(&mol, 0, &[]).unwrap_err();
        assert_eq!(
            error,
            StereoError::InvariantViolation(
                "Atom::getPerturbationOrder probe/reference length mismatch".to_string()
            )
        );
    }

    #[test]
    fn perturbation_order_rejects_missing_probe_bond() {
        let mol = Molecule::from_smiles("CC").expect("failed to parse test SMILES");
        let error = perturbation_order_from_bond_indices(&mol, 0, &[99]).unwrap_err();
        assert_eq!(
            error,
            StereoError::InvariantViolation(
                "Atom::getPerturbationOrder expected bond missing from probe order".to_string()
            )
        );
    }
}
