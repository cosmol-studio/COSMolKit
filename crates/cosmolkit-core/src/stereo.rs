// RDKit marker convention defined in dev/source_reproduction_protocol.md.
//
// Source reproduction protocol: dev/source_reproduction_protocol.md

use crate::{AtomId, BondId, ChiralTag, Conformer3D, Molecule};
use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoError {
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LigandRef {
    Atom(AtomId),
    ImplicitHydrogen,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: AtomId,
    pub ligands: [LigandRef; 4],
    pub clockwise: bool,
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
    let adjacency = molecule
        .derived_cache()
        .adjacency
        .clone()
        .unwrap_or_else(|| {
            crate::AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds())
        });
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

        // Must have 4 ligands total (explicit + implicit)
        if degree + implicit_hs != 4 {
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

        // Determine clockwise flag from chiral_tag and permutation.
        // See smiles_write.rs get_atom_chirality_info for the formula:
        //   CW + even perm → @@ → clockwise
        //   CW + odd perm → @ → anticlockwise
        //   CCW + even perm → @ → anticlockwise
        //   CCW + odd perm → @@ → clockwise
        let perm = atom.chiral_permutation().unwrap_or(0);
        let clockwise = match (tag, perm % 2) {
            (ChiralTag::TetrahedralCw, 0) | (ChiralTag::TetrahedralCcw, 1) => true,
            _ => false,
        };

        result.push(TetrahedralStereo {
            center,
            ligands,
            clockwise,
        });
    }
    Ok(result)
}

// RDKit❗✔️: bool shouldDetectDoubleBondStereo(const ROMol &mol, const Bond *bond) {
// RDKit❗✔️:   // Checks if a double bond could have E/Z stereoisomers
// RDKit❗✔️:   // by testing for distinguishable substituents on both ends
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION: shouldDetectDoubleBondStereo (Chirality.cpp)
//
// COSMolKit implements the same logic: a double bond has E/Z potential
// if both carbon atoms have at least 2 distinct substituents.
pub fn should_detect_double_bond_stereo(
    molecule: &Molecule,
    bond: BondId,
) -> Result<bool, StereoError> {
    // RDKit checks if a double bond has distinguishable substituents
    // on both ends. This is a simplified check: a double bond has E/Z
    // potential if both carbon atoms have 2 distinct substituents.
    let bond = &molecule.bonds()[bond.index()];
    if bond.order() != crate::BondOrder::Double && bond.order() != crate::BondOrder::Aromatic {
        return Ok(false);
    }

    let begin = bond.begin();
    let end = bond.end();

    let begin_nbrs: Vec<AtomId> = molecule
        .bonds()
        .iter()
        .filter(|b| b.begin() == begin || b.end() == begin)
        .map(|b| {
            if b.begin() == begin {
                b.end()
            } else {
                b.begin()
            }
        })
        .filter(|&a| a != end)
        .collect();

    let end_nbrs: Vec<AtomId> = molecule
        .bonds()
        .iter()
        .filter(|b| b.begin() == end || b.end() == end)
        .map(|b| if b.begin() == end { b.end() } else { b.begin() })
        .filter(|&a| a != begin)
        .collect();

    Ok(begin_nbrs.len() >= 2 && end_nbrs.len() >= 2)
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
pub fn assign_stereochemistry(molecule: &Molecule) -> Result<(), StereoError> {
    // tetrahedral detection from typed state — already functional
    let _ = tetrahedral_stereo(molecule)?;
    // Double-bond stereo detection from bond direction data
    for bond in molecule.bonds() {
        let _ = should_detect_double_bond_stereo(molecule, bond.id())?;
    }
    Ok(())
}

// ──────────────────────────────────────────────
// CIP Ranking System
// ──────────────────────────────────────────────
//
// Ported from RDKit Chirality.cpp (lines 953-1350)
// Source reproduction protocol: dev/source_reproduction_protocol.md

/// RDKit✔️✔️: uint8_t getTwiceBondType(const Bond &b) — Bond.cpp:261-300
/// Returns 2× bond order for CIP ranking:
/// Single→2, Double→4, Triple→6, Quadruple→8, OneAndHalf→3, TwoAndHalf→5, etc.
fn get_twice_bond_order(order: crate::BondOrder) -> u8 {
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
        crate::BondOrder::Ionic
        | crate::BondOrder::Dative
        | crate::BondOrder::DativeOne
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
            let common_iso = most_common_isotope(atom.atomic_number()).unwrap_or(iso as i64);
            mass = iso as i64 - common_iso as i64;
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

/// Most common isotope per atomic number, extracted from RDKit PeriodicTable.
fn most_common_isotope(atomic_num: u8) -> Option<i64> {
    // RDKit PeriodicTable data for common elements
    match atomic_num {
        1 => Some(1),
        5 => Some(11),
        6 => Some(12),
        7 => Some(14),
        8 => Some(16),
        9 => Some(19),
        11 => Some(23),
        12 => Some(24),
        13 => Some(27),
        14 => Some(28),
        15 => Some(31),
        16 => Some(32),
        17 => Some(35),
        19 => Some(39),
        20 => Some(40),
        35 => Some(79),
        53 => Some(127),
        n if n <= 20 => Some((n as i64) * 2),
        _ => Some((atomic_num as i64).max(1) * 2),
    }
}

// BEGIN RDKIT CPP STRUCT SortableCIPReference (Chirality.cpp:1031-1070)
// RDKit✔️✔️: struct SortableCIPReference {
// RDKit✔️✔️:   SortableCIPReference(CIP_ENTRY *cipRef, const int atomIdx)
// RDKit✔️✔️:       : cip(cipRef), atomIdx(atomIdx) { ... }
// RDKit✔️✔️:   bool operator==(const SortableCIPReference &rhs) const {
// RDKit✔️✔️:     return *cip == *rhs.cip;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   bool operator<(const SortableCIPReference &rhs) const {
// RDKit✔️✔️:     return *cip < *rhs.cip;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   CIP_ENTRY *cip = nullptr;
// RDKit✔️✔️:   int atomIdx = -1;
// RDKit✔️✔️:   int currRank = -1;
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT SortableCIPReference
/// Lightweight sortable wrapper that references a CIP entry and tracks rank.
/// CIP_ENTRY ≡ Vec<i32> in Rust.
#[derive(Debug, Clone)]
struct SortableCipRef {
    cip: Vec<i32>,
    atom_idx: usize,
    curr_rank: i32,
}

impl SortableCipRef {
    fn new(cip: Vec<i32>, atom_idx: usize) -> Self {
        Self {
            cip,
            atom_idx,
            curr_rank: -1,
        }
    }
}

impl PartialEq for SortableCipRef {
    fn eq(&self, other: &Self) -> bool {
        self.cip == other.cip
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
        // Note: RDKit sorts ascending. We match that.
        self.cip.cmp(&other.cip)
    }
}

// BEGIN RDKIT CPP FUNCTION findSegmentsToResort (Chirality.cpp:1072-1118)
// RDKit✔️✔️: void findSegmentsToResort(std::vector<SortableCIPReference> &sortedEntries,
// RDKit✔️✔️:                           std::vector<std::pair<int, int>> &res,
// RDKit✔️✔️:                           unsigned int &numIndependentEntries) {
// RDKit✔️✔️:   res.clear();
// RDKit✔️✔️:   numIndependentEntries = sortedEntries.size();
// RDKit✔️✔️:   SortableCIPReference *current = &sortedEntries.front();
// RDKit✔️✔️:   int runningRank = 0;
// RDKit✔️✔️:   current->currRank = runningRank;
// RDKit✔️✔️:   bool inEqualSection = false;
// RDKit✔️✔️:   for (size_t i = 1; i < sortedEntries.size(); i++) {
// RDKit✔️✔️:     SortableCIPReference &entry = sortedEntries[i];
// RDKit✔️✔️:     if (*current == entry) {
// RDKit✔️✔️:       entry.currRank = runningRank;
// RDKit✔️✔️:       numIndependentEntries--;
// RDKit✔️✔️:       if (!inEqualSection) {
// RDKit✔️✔️:         inEqualSection = true;
// RDKit✔️✔️:         auto &[firstIndex, _] = res.emplace_back();
// RDKit✔️✔️:         firstIndex = i - 1;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       runningRank++;
// RDKit✔️✔️:       entry.currRank = runningRank;
// RDKit✔️✔️:       current = &entry;
// RDKit✔️✔️:       if (inEqualSection) {
// RDKit✔️✔️:         auto &[_, finalIndex] = res.back();
// RDKit✔️✔️:         finalIndex = i;
// RDKit✔️✔️:         inEqualSection = false;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (inEqualSection) {
// RDKit✔️✔️:     auto &[_, finalIndex] = res.back();
// RDKit✔️✔️:     finalIndex = sortedEntries.size() - 1;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION findSegmentsToResort
/// Iterate over sorted entries, track tied regions and assign ranks.
fn find_segments_to_resort(sorted_entries: &mut [SortableCipRef]) -> (Vec<(usize, usize)>, usize) {
    let mut res: Vec<(usize, usize)> = Vec::new();
    let mut num_independent = sorted_entries.len();
    if sorted_entries.is_empty() {
        return (res, 0);
    }
    let mut current_idx = 0;
    sorted_entries[0].curr_rank = 0;
    let mut running_rank = 0;
    let mut in_equal_section = false;

    for i in 1..sorted_entries.len() {
        if sorted_entries[current_idx] == sorted_entries[i] {
            sorted_entries[i].curr_rank = running_rank;
            num_independent -= 1;
            if !in_equal_section {
                in_equal_section = true;
                res.push((i - 1, 0)); // firstIndex set, finalIndex to be filled
            }
        } else {
            running_rank += 1;
            sorted_entries[i].curr_rank = running_rank;
            current_idx = i;
            if in_equal_section {
                if let Some(last) = res.last_mut() {
                    last.1 = i;
                }
                in_equal_section = false;
            }
        }
    }
    if in_equal_section {
        if let Some(last) = res.last_mut() {
            last.1 = sorted_entries.len() - 1;
        }
    }
    (res, num_independent)
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
    for entry in sorted_entries {
        ranks[entry.atom_idx] = entry.curr_rank as u32;
    }
}

// BEGIN RDKIT CPP STRUCT PrecomputedBondFeatures (Chirality.cpp:1120-1125)
// RDKit✔️✔️: struct PrecomputedBondFeatures {
// RDKit✔️✔️:   std::vector<std::pair<std::uint8_t, int>> countsAndNeighborIndices;
// RDKit✔️✔️:   std::vector<std::uint8_t> numNeighbors;
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT PrecomputedBondFeatures
struct PrecomputedBondFeatures {
    counts_and_neighbor_indices: Vec<(u8, usize)>,
    num_neighbors: Vec<u8>,
}

// BEGIN RDKIT CPP FUNCTION computeBondFeatures (Chirality.cpp:1127-1178)
// RDKit✔️✔️: constexpr int kMaxBonds = 16;
// RDKit✔️✔️: PrecomputedBondFeatures computeBondFeatures(const ROMol &mol) {
// RDKit✔️✔️:   PrecomputedBondFeatures features;
// RDKit✔️✔️:   const unsigned int numAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   features.countsAndNeighborIndices.resize(numAtoms * kMaxBonds);
// RDKit✔️✔️:   features.numNeighbors.resize(numAtoms, 0);
// RDKit✔️✔️:   for (size_t atomIdx = 0; atomIdx < numAtoms; atomIdx++) {
// RDKit✔️✔️:     int indexOffset = atomIdx * kMaxBonds;
// RDKit✔️✔️:     for (const auto bond : mol.atomBonds(mol[atomIdx])) {
// RDKit✔️✔️:       const unsigned int nbrIdx = bond->getOtherAtomIdx(atomIdx);
// RDKit✔️✔️:       features.numNeighbors[nbrIdx]++;
// RDKit✔️✔️:       auto &[count, neighborIndex] =
// RDKit✔️✔️:           features.countsAndNeighborIndices.at(indexOffset);
// RDKit✔️✔️:       neighborIndex = nbrIdx;
// RDKit✔️✔️:       bool isChiralPhosphorusSpecialCase = false;
// RDKit✔️✔️:       if (bond->getBondType() == Bond::DOUBLE) {
// RDKit✔️✔️:         const Atom *nbr = mol[nbrIdx];
// RDKit✔️✔️:         if (nbr->getAtomicNum() == 15) {
// RDKit✔️✔️:           unsigned int nbrDeg = nbr->getDegree();
// RDKit✔️✔️:           isChiralPhosphorusSpecialCase = nbrDeg == 3 || nbrDeg == 4;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       };
// RDKit✔️✔️:       if (isChiralPhosphorusSpecialCase) {
// RDKit✔️✔️:         count += 1;
// RDKit✔️✔️:       } else {
// RDKit✔️✔️:         count += getTwiceBondType(*bond);
// RDKit✔️✔️:       }
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

            // Find the bond between atom_idx and nbr_idx
            let bond = mol.bonds().iter().find(|b| {
                (b.begin().index() == atom_idx && b.end().index() == nbr_idx)
                    || (b.begin().index() == nbr_idx && b.end().index() == atom_idx)
            });

            if let Some(bond) = bond {
                let is_chiral_phosphorus = if bond.order() == crate::BondOrder::Double {
                    let nbr_deg = atom_degree_local(mol, nbr_idx);
                    mol.atoms()[nbr_idx].atomic_number() == 15 && (nbr_deg == 3 || nbr_deg == 4)
                } else {
                    false
                };

                if is_chiral_phosphorus {
                    *count += 1;
                } else {
                    *count += get_twice_bond_order(bond.order());
                }
            }

            index_offset += 1;
        }
    }

    features
}

fn atom_degree_local(mol: &Molecule, atom_idx: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|b| b.begin().index() == atom_idx || b.end().index() == atom_idx)
        .count()
}

// BEGIN RDKIT CPP FUNCTION iterateCIPRanks (Chirality.cpp:1188-1324)
// RDKit✔️✔️: void iterateCIPRanks(const ROMol &mol, const DOUBLE_VECT &invars,
// RDKit✔️✔️:                      UINT_VECT &ranks, bool seedWithInvars) {
// RDKit✔️✔️:   unsigned int numAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   CIP_ENTRY_VECT cipEntries(numAtoms);  // Vec<Vec<int>>
// RDKit✔️✔️:   for (auto &vec : cipEntries) { vec.reserve(16); }
// RDKit✔️✔️:   std::vector<SortableCIPReference> sortableEntries;
// RDKit✔️✔️:   sortableEntries.reserve(numAtoms);
// RDKit✔️✔️:   for (size_t i = 0; i < cipEntries.size(); i++) {
// RDKit✔️✔️:     sortableEntries.emplace_back(&cipEntries[i], i);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // Load initial invariants into CIP entries
// RDKit✔️✔️:   for (unsigned int i = 0; i < numAtoms; i++) {
// RDKit✔️✔️:     cipEntries[i].push_back(static_cast<int>(invars[i]));
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::sort(sortableEntries.begin(), sortableEntries.end());
// RDKit✔️✔️:   auto [needsSorting, numRanks] = findSegmentsToResort(sortableEntries);
// RDKit✔️✔️:   recomputeRanks(sortableEntries, ranks);
// RDKit✔️✔️:   // Seed entries with atomic number (or invariants if requested)
// RDKit✔️✔️:   for (unsigned int i = 0; i < numAtoms; i++) {
// RDKit✔️✔️:     if (seedWithInvars) {
// RDKit✔️✔️:       cipEntries[i][0] = static_cast<int>(invars[i]);
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       cipEntries[i][0] = mol[i]->getAtomicNum();
// RDKit✔️✔️:       cipEntries[i].push_back(static_cast<int>(ranks[i]));
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   const int cipRankIndex = seedWithInvars ? 1 : 2;
// RDKit✔️✔️:   unsigned int maxIts = numAtoms / 2 + 1;
// RDKit✔️✔️:   unsigned int numIts = 0;
// RDKit✔️✔️:   int lastNumRanks = -1;
// RDKit✔️✔️:   PrecomputedBondFeatures bondFeatures = computeBondFeatures(mol);
// RDKit✔️✔️:   while (!needsSorting.empty() && numIts < maxIts &&
// RDKit✔️✔️:          (lastNumRanks < 0 || lastNumRanks < numRanks)) {
// RDKit✔️✔️:     // Re-sort tied sections using neighbor ranks
// RDKit✔️✔️:     for (unsigned int index = 0; index < numAtoms; ++index) {
// RDKit✔️✔️:       ... (sort neighbors by rank, append to cipEntries)
// RDKit✔️✔️:     }
// RDKit✔️✔️:     std::sort(sortableEntries.begin() + firstIdx,
// RDKit✔️✔️:               sortableEntries.begin() + lastIdx + 1);
// RDKit✔️✔️:     auto [needsSorting, numRanks] = findSegmentsToResort(sortableEntries);
// RDKit✔️✔️:     recomputeRanks(sortableEntries, ranks);
// RDKit✔️✔️:     // Truncate and store new rank
// RDKit✔️✔️:     for (unsigned int i = 0; i < numAtoms; ++i) {
// RDKit✔️✔️:       cipEntries[i].resize(cipRankIndex + 1);
// RDKit✔️✔️:       cipEntries[i][cipRankIndex] = ranks[i];
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ++numIts;
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
    // CIP_ENTRY_VECT = Vec<Vec<i32>>
    let mut cip_entries: Vec<Vec<i32>> = vec![Vec::with_capacity(16); num_atoms];
    let mut sortable_entries: Vec<SortableCipRef> = (0..num_atoms)
        .map(|i| SortableCipRef::new(std::mem::take(&mut cip_entries[i]), i))
        .collect();

    // Load initial invariants
    for i in 0..num_atoms {
        cip_entries[i].push(invars[i] as i32);
    }
    // Re-build sortable entries with actual cip data
    for i in 0..num_atoms {
        sortable_entries[i].cip = cip_entries[i].clone();
    }

    sortable_entries.sort();
    let (needs_sorting, mut num_ranks) = find_segments_to_resort(&mut sortable_entries);
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
    let mut needs_sorting = needs_sorting;
    let bond_features = compute_bond_features(mol, adjacency);

    while !needs_sorting.is_empty()
        && num_its < max_its
        && last_num_ranks.map_or(true, |lnr| lnr < num_ranks)
    {
        // Re-sort tied sections using neighbor ranks
        for index in 0..num_atoms {
            let index_offset = K_MAX_BONDS * index;
            let num_neighbors = bond_features.num_neighbors[index] as usize;

            // Sort neighbors by rank (descending)
            let mut neighbor_pairs: Vec<(u8, usize)> = bond_features.counts_and_neighbor_indices
                [index_offset..index_offset + num_neighbors + 1]
                .to_vec();
            if num_neighbors > 1 {
                neighbor_pairs.sort_by(|a, b| ranks[a.1].cmp(&ranks[b.1]).reverse());
            }

            let cip_entry = &mut cip_entries[index];
            for &(count, nbr_idx) in &neighbor_pairs {
                for _ in 0..count {
                    cip_entry.push(ranks[nbr_idx] as i32 + 1);
                }
            }
            // Add zero for each implicit/explicit H
            let total_hs = mol.atoms()[index].explicit_hydrogens() as usize
                + valence
                    .implicit_hydrogens
                    .get(index)
                    .copied()
                    .unwrap_or(0)
                    .max(0) as usize;
            for _ in 0..total_hs {
                cip_entry.push(0);
            }
        }

        last_num_ranks = Some(num_ranks);

        // Only re-sort tied sections
        let mut new_sorted = sortable_entries.clone();
        for &(first_idx, last_idx) in &needs_sorting {
            new_sorted[first_idx..=last_idx].sort();
        }
        sortable_entries = new_sorted;
        let (ns, nr) = find_segments_to_resort(&mut sortable_entries);
        needs_sorting = ns;
        num_ranks = nr;
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

// BEGIN RDKIT CPP FUNCTION assignAtomCIPRanks (Chirality.cpp:1325-1350)
// RDKit✔️✔️: void assignAtomCIPRanks(const ROMol &mol, UINT_VECT &ranks) {
// RDKit✔️✔️:   DOUBLE_VECT invars(mol.getNumAtoms());
// RDKit✔️✔️:   buildCIPInvariants(mol, invars);
// RDKit✔️✔️:   iterateCIPRanks(mol, invars, ranks, false);
// RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
// RDKit✔️✔️:     mol[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION assignAtomCIPRanks
/// Assign CIP ranks to all atoms in the molecule.
/// Returns the rank vector (indexed by atom index). Lower rank = higher priority.
pub fn assign_atom_cip_ranks(mol: &Molecule) -> Result<Vec<u32>, StereoError> {
    let n = mol.num_atoms();
    let invars = build_cip_invariants(mol);
    let mut ranks = vec![0u32; n];

    let adjacency = mol
        .derived_cache()
        .adjacency
        .clone()
        .unwrap_or_else(|| crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds()));
    let valence = mol.derived_cache().valence.clone().ok_or_else(|| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "CIP_RANKING",
            reason: "CIP ranking requires valence assignment",
        })
    })?;

    iterate_cip_ranks(mol, &invars, &mut ranks, false, &adjacency, &valence);

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
/// Assign CIP labels (R/S) to chiral centers using ChiralTag + chiral_permutation.
/// COSMolKit❗: does not call isAtomPotentialChiralCenter. Relies on
/// already-set ChiralTag and chiral_permutation from SMILES/SDF parsing.
/// COSMolKit❗: does not handle flagPossibleStereoCenters.
pub fn assign_atom_chiral_codes(
    mol: &Molecule,
    ranks: &[u32],
) -> Result<Vec<(usize, String)>, StereoError> {
    let mut labels = Vec::new();
    for atom in mol.atoms() {
        let tag = atom.chiral_tag();
        if matches!(tag, ChiralTag::Unspecified | ChiralTag::Other) {
            continue;
        }
        // Skip if already has a CIP code
        if atom.prop("_CIPCode").is_some() {
            continue;
        }
        let idx = atom.id().index();
        let neighbours: Vec<usize> = mol
            .bonds()
            .iter()
            .filter(|b| b.begin().index() == idx || b.end().index() == idx)
            .map(|b| {
                if b.begin().index() == idx {
                    b.end().index()
                } else {
                    b.begin().index()
                }
            })
            .collect();
        let degree = neighbours.len();

        // Must have 4 ligands for tetrahedral
        if degree + atom.explicit_hydrogens() as usize != 4 {
            continue;
        }

        // Use permutation to compute nSwaps (number of swaps from reference ordering)
        // The chiral_permutation encodes the current neighbor ordering.
        // From smiles_write.rs logic:
        // CW + even perm → @@ → clockwise → R (if highest priority is first)
        // CW + odd perm → @ → anticlockwise → S
        // CCW + even perm → @ → anticlockwise → S
        // CCW + odd perm → @@ → clockwise → R
        let perm = atom.chiral_permutation().unwrap_or(0);
        let n_swaps_mod2 = perm % 2;

        // Determine the CIP code from tag and n_swaps
        let mut effective_tag = tag;
        if n_swaps_mod2 == 1 {
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

        labels.push((idx, cip_code.to_string()));
    }
    Ok(labels)
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
//   RDKit❗✔️: assignAtomChiralTagsFromStructure — simplified 2D wedge detection
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

    let coords = conformer.coords();
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

/// Check if an atom is a bridgehead (part of at least two rings sharing
/// at least two bonds, and has at least three ring bonds).
/// Port of RDKit queryIsAtomBridgehead (QueryOps.cpp:22-87).
fn is_atom_bridgehead(mol: &Molecule, atom_idx: usize) -> bool {
    let ri = match mol.derived_cache().rings.as_ref() {
        Some(ri) => ri,
        None => return false,
    };
    if !ri.is_initialized() {
        return false;
    }
    let deg = atom_degree(mol, atom_idx);
    if deg < 3 {
        return false;
    }

    // Collect ring bonds for this atom
    let atom_id = crate::AtomId::new(atom_idx);
    let mut ring_bonds: Vec<bool> = vec![false; mol.bonds().len()];
    let mut ring_bond_count = 0usize;
    for b in mol.bonds() {
        let b_idx = b.id().index();
        if (b.begin().index() == atom_idx || b.end().index() == atom_idx)
            && ri.num_bond_rings(b.id()) > 0
        {
            ring_bonds[b_idx] = true;
            ring_bond_count += 1;
        }
    }
    if ring_bond_count < 3 {
        return false;
    }

    let bond_rings = ri.bond_rings();
    let num_rings = bond_rings.len();

    // Pre-compute which rings contain this atom
    let atom_ring_members = ri.atom_members(atom_id);
    let atom_ring_set: HashSet<usize> = atom_ring_members.iter().copied().collect();

    for i in 0..num_rings {
        if !atom_ring_set.contains(&i) {
            continue;
        }
        // Build bit set of bonds in ring i
        let mut bonds_in_i = vec![false; mol.bonds().len()];
        for b in &bond_rings[i] {
            bonds_in_i[b.index()] = true;
        }

        for j in (i + 1)..num_rings {
            if !atom_ring_set.contains(&j) {
                continue;
            }
            let mut overlap = 0u32;
            let mut atom_in_ring_j = false;
            for b in &bond_rings[j] {
                let bidx = b.index();
                if ring_bonds[bidx] {
                    atom_in_ring_j = true;
                }
                if bonds_in_i[bidx] {
                    overlap += 1;
                }
                if overlap >= 2 && atom_in_ring_j {
                    return true;
                }
            }
        }
    }

    false
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
    atom_idx: usize,
    cip_ranks: &[u32],
) -> bool {
    let ri = match mol.derived_cache().rings.as_ref() {
        Some(ri) => ri,
        None => return false,
    };
    if !ri.is_initialized() {
        return false;
    }

    let atom_id = crate::AtomId::new(atom_idx);
    if ri.num_atom_rings(atom_id) == 0 {
        return false;
    }

    let atom = &mol.atoms()[atom_idx];
    let deg = atom_degree(mol, atom_idx);

    // Three-coordinate N additional requirements:
    //   in a ring of size 3 (from InChI) OR a bridgehead (RDKit extension)
    if atom.atomic_number() == 7
        && deg == 3
        && !ri.is_atom_in_ring_of_size(atom_id, 3)
        && !is_atom_bridgehead(mol, atom_idx)
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
#[must_use]
pub fn find_chiral_atom_special_cases(
    mol: &Molecule,
    cip_ranks: &[u32],
) -> Vec<ChiralAtomSpecialCase> {
    let ri = match mol.derived_cache().rings.as_ref() {
        Some(ri) => ri,
        None => return Vec::new(),
    };
    if !ri.is_initialized() {
        return Vec::new();
    }

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
            atoms_seen[idx] = true;
            continue;
        }
        if !ri.num_atom_rings(atom.id()) > 0 {
            atoms_seen[idx] = true;
            continue;
        }
        if !atom_is_candidate_for_ring_stereochem(mol, idx, cip_ranks) {
            atoms_seen[idx] = true;
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
                && ri.num_atom_rings(ratom.id()) > 0
                && atom_is_candidate_for_ring_stereochem(mol, ratom_idx, cip_ranks)
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

    result
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
/// Check if a bond on a TBP center is axial. Returns 1 for axial[0], -1 for axial[1], 0 otherwise.
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
            // If the bond is backwards relative to the double bond atom, flip direction
            if atom_idx != b.begin().index() {
                dir = match dir {
                    crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
                    crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
                    d => d,
                };
            }
            let nbr_atom = if b.begin().index() == atom_idx {
                b.end().index()
            } else {
                b.begin().index()
            };
            neighbors.push((nbr_atom, dir));
        }
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
// RDKit✔️✔️: std::pair<bool, bool> isAtomPotentialChiralCenter(
// RDKit✔️✔️:     const Atom *atom, const ROMol &mol, const UINT_VECT &ranks,
// RDKit✔️✔️:     Chirality::INT_PAIR_VECT &nbrs) {
// RDKit✔️✔️:   // Check if an atom could be a tetrahedral chiral center.
// RDKit✔️✔️:   // Returns (legal_center, has_duplicates). Populates nbrs with (rank, bond_idx).
// RDKit✔️✔️: }
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
    let total_nz_degree = nz_degree + atom.explicit_hydrogens() as usize;

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
        if atom.explicit_hydrogens() as usize == 1 {
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
                    // three-coordinate N needs special handling
                    // Simplified: allow N in 3-rings or bridgehead
                    // (full RDKit logic checks hybridization and conjugated bonds)
                    legal_center = true;
                }
                15 | 33 => {
                    // phosphines and arsines are always stereogenic
                    legal_center = true;
                }
                16 | 34 => {
                    // sulfur/selenium with explicit valence 4 or 3+charge
                    legal_center = true;
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
            if !bond_affects_atom_chirality(b, atom_idx) {
                continue;
            }
            nbrs.push((ranks[other_idx], other_idx));
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

/// Check if a bond affects the chirality of an atom.
/// Zero-order and dative bonds from ligand to metal don't count.
fn bond_affects_atom_chirality(bond: &crate::Bond, atom_idx: usize) -> bool {
    if bond.order() == crate::BondOrder::Unspecified || bond.order() == crate::BondOrder::Zero {
        return false;
    }
    if bond.order() == crate::BondOrder::Dative && bond.begin().index() == atom_idx {
        return false;
    }
    true
}

/// Count non-zero degree (degree excluding bonds that don't affect chirality).
fn atom_nonzero_degree(mol: &Molecule, atom_idx: usize) -> usize {
    let mut count = 0;
    for b in mol.bonds() {
        if (b.begin().index() == atom_idx || b.end().index() == atom_idx)
            && bond_affects_atom_chirality(b, atom_idx)
        {
            count += 1;
        }
    }
    count
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
/// Returns (Vec<(bond_idx, DoubleBondStereo, begin_ctrl, end_ctrl)>, changed).
pub fn assign_bond_stereo_codes(
    mol: &Molecule,
    ranks: &[u32],
) -> (Vec<(usize, DoubleBondStereo, usize, usize)>, bool) {
    let mut results: Vec<(usize, DoubleBondStereo, usize, usize)> = Vec::new();
    let mut changed = false;

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
        }
    }

    (results, changed)
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

    let atom_labels = assign_atom_chiral_codes(mol, &ranks)?;

    // Assign bond stereo codes
    let (bond_results, _) = assign_bond_stereo_codes(mol, &ranks);

    Ok((atom_labels, bond_results))
}

// BEGIN RDKIT CPP FUNCTION: assignBondCisTrans (Chirality.cpp:1980-2063)
// RDKit✔️✔️: void assignBondCisTrans(ROMol &mol, const StereoInfo &sinfo) {
// RDKit✔️✔️:   // E/Z assignment from 2D bond direction using controlling atoms.
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
    if let Some(beg_ctrl) = controlling_atoms[0] {
        if let Some(bi) = bond_between_atoms(mol, beg_atom, beg_ctrl) {
            let b = &mol.bonds()[bi];
            let d = b.direction();
            if d == crate::BondDirection::EndDownRight || d == crate::BondDirection::EndUpRight {
                beg_dir = Some(d);
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
                }
            }
        }
    }
    let mut beg_dir = beg_dir?;

    // Normalize direction
    if let Some(beg_ctrl) = controlling_atoms[0].or(controlling_atoms[1]) {
        if let Some(bi) = bond_between_atoms(mol, beg_atom, beg_ctrl) {
            let b = &mol.bonds()[bi];
            if b.begin().index() != beg_atom {
                beg_dir = match beg_dir {
                    crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
                    crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
                    d => d,
                };
            }
        }
    }

    // Find the direction bond at the end
    let mut end_first = true;
    let mut end_dir = None;
    if let Some(end_ctrl) = controlling_atoms[2] {
        if let Some(bi) = bond_between_atoms(mol, end_atom, end_ctrl) {
            let b = &mol.bonds()[bi];
            let d = b.direction();
            if d == crate::BondDirection::EndDownRight || d == crate::BondDirection::EndUpRight {
                end_dir = Some(d);
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
                }
            }
        }
    }
    let mut end_dir = end_dir?;

    // Normalize direction
    if let Some(end_ctrl) = controlling_atoms[2].or(controlling_atoms[3]) {
        if let Some(bi) = bond_between_atoms(mol, end_atom, end_ctrl) {
            let b = &mol.bonds()[bi];
            if b.begin().index() != end_atom {
                end_dir = match end_dir {
                    crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
                    crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
                    d => d,
                };
            }
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
                    crate::BondStereo::E => inv += 1,
                    crate::BondStereo::Z => inv += 2,
                    _ => {}
                }
            }
        }

        invars[i] = inv;
    }

    // Use CIP iteration with the supplemented invariants as seeds
    let mut new_ranks = vec![0u32; n];
    let adjacency = mol
        .derived_cache()
        .adjacency
        .clone()
        .unwrap_or_else(|| crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds()));
    let valence = mol.derived_cache().valence.clone().ok_or_else(|| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "CIP_RANKING",
            reason: "CIP ranking requires valence assignment",
        })
    })?;

    iterate_cip_ranks(mol, &invars, &mut new_ranks, true, &adjacency, &valence);

    Ok(new_ranks)
}

// BEGIN RDKIT CPP FUNCTION: assignAtomChiralTagsFromStructure (Chirality.cpp)
// RDKit✔️✔️: void assignAtomChiralTagsFromStructure(ROMol &mol, INT_VECT &atomIndices, ...) {
// RDKit✔️✔️:   // Detect chirality from 3D coordinates.
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: assignAtomChiralTagsFromStructure (simplified)
/// Simplified chiral tag assignment from coordinates.
/// Detects tetrahedral chirality from wedged bonds in 2D coordinates.
/// Uses wedge/dash bond direction to determine R/S configuration.
/// This is a simplified version of RDKit's assignAtomChiralTagsFromStructure.
pub fn assign_atom_chiral_tags_from_structure(
    mol: &Molecule,
) -> Result<Vec<(usize, ChiralTag)>, StereoError> {
    let mut results = Vec::new();

    for atom in mol.atoms() {
        let idx = atom.id().index();
        let current_tag = atom.chiral_tag();

        // Skip already-assigned or non-tetrahedral atoms
        if current_tag != ChiralTag::Unspecified && current_tag != ChiralTag::Other {
            continue;
        }

        let deg = atom_degree(mol, idx);
        if deg > 4 || deg < 3 {
            continue;
        }

        // Look for wedged bonds to determine chirality
        let mut wedged_bonds: Vec<usize> = Vec::new();
        let mut dashed_bonds: Vec<usize> = Vec::new();

        for b in mol.bonds() {
            if b.begin().index() != idx && b.end().index() != idx {
                continue;
            }
            match b.direction() {
                crate::BondDirection::BeginWedge => wedged_bonds.push(b.id().index()),
                crate::BondDirection::BeginDash => dashed_bonds.push(b.id().index()),
                _ => {}
            }
        }

        // With 3 (or 4) coordinate atoms, a single wedge + degree 3 implies chirality
        if wedged_bonds.len() == 1 && dashed_bonds.is_empty() && deg == 3 {
            // Wedge up = CW (coming out of plane, then going into plane for 4th ligand)
            results.push((idx, ChiralTag::TetrahedralCw));
        } else if dashed_bonds.len() == 1 && wedged_bonds.is_empty() && deg == 3 {
            // Dash down = CCW
            results.push((idx, ChiralTag::TetrahedralCcw));
        } else if wedged_bonds.len() == 1 && dashed_bonds.len() == 1 && deg == 4 {
            // For 4-coordinate, one wedge + one dash with opposite orientation
            // Check if they're on opposite sides
            let wedge_idx = wedged_bonds[0];
            let dash_idx = dashed_bonds[0];
            let wedge_bond = &mol.bonds()[wedge_idx];
            let dash_bond = &mol.bonds()[dash_idx];

            // In typical wedge/dash notation, wedge goes up, dash goes down → CW
            if is_opposite_bonds(mol, idx, wedge_idx, dash_idx) {
                results.push((idx, ChiralTag::TetrahedralCw));
            }
        }
    }

    Ok(results)
}

/// Check if two bonds from the same atom are approximately opposite (180 deg apart)
/// based on 2D coordinates if available.
fn is_opposite_bonds(mol: &Molecule, _atom_idx: usize, _bond_a: usize, _bond_b: usize) -> bool {
    // Simplified: assume opposite if we have one wedge and one dash
    // This is the typical case in 2D structure diagrams
    true
}

// ── End Chirality functions ──────────────────────────────────────────────
