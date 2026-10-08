use crate::stereo_graph::{BondRows, GraphAdjacency, NeighborRows, StereoGraphAccess};
// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use cosmolkit_model::{
    AdjacencyError, AdjacencyList, AtomId, Bond, BondId, MoleculeProperties, MoleculePropertyError,
    NeighborRef, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::BondOrder;

const MAX_BFSQ_SIZE: usize = 200_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingFindType {
    OtherOrUnknown,
    Fast,
    Sssr,
    SymmSssr,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RingSearchParams {
    pub include_dative_bonds: bool,
    pub include_hydrogen_bonds: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RingFindingError {
    #[error("source FIND_RING_TYPE value {tag} is outside the modeled enum")]
    SourceRingFindType { tag: u32 },

    #[error("molecule property lifecycle failed: {0}")]
    MoleculeProperty(#[from] MoleculePropertyError),
    #[error("{message}")]
    Value { message: &'static str },
    #[error("expected bond not found between atom {begin} and atom {end}")]
    ExpectedBondNotFound { begin: AtomId, end: AtomId },
    #[error("ring atom index {atom} is out of range for {atom_count} atoms")]
    RingAtomOutOfRange { atom: usize, atom_count: usize },
    #[error("ring bond index {bond} is out of range for {bond_count} bonds")]
    RingBondOutOfRange { bond: usize, bond_count: usize },
    #[error("ring-decomposition edge {edge} has no original bond mapping")]
    RingFamilyEdgeMappingMissing { edge: usize },
    #[error("relevant-cycle count is outside the source unsigned-integer range")]
    RelevantCycleCountOutOfRange,
    #[error("unsupported ring info branch: {reason}")]
    UnsupportedBranch { reason: &'static str },
    #[error(transparent)]
    Adjacency(#[from] AdjacencyError),
    #[error(transparent)]
    InvalidTopology(#[from] TopologyValidationError),
    #[error(transparent)]
    RingDecomposer(#[from] cosmolkit_ringdecomposer::RingDecomposerError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingInfo {
    initialized: bool,
    find_type: RingFindType,
    atom_members: Vec<Vec<usize>>,
    bond_members: Vec<Vec<usize>>,
    atom_rings: Vec<Vec<AtomId>>,
    bond_rings: Vec<Vec<BondId>>,
    atom_ring_families: Vec<Vec<AtomId>>,
    bond_ring_families: Vec<Vec<BondId>>,
    relevant_cycle_count: Option<usize>,
    fused_rings: Vec<Vec<bool>>,
    num_fused_bonds: Vec<usize>,
}

impl RingInfo {
    /// Move all represented source cache fields into the detached MODEL value.
    #[doc(hidden)]
    pub fn into_source_snapshot(self) -> cosmolkit_model::SourceRingInfo {
        // RDKit❗✔️:   RingInfo(RingInfo &&other) noexcept = default;
        // RDKit❗✔️: typedef enum {
        // RDKit❗✔️:   FIND_RING_TYPE_FAST,
        // RDKit❗✔️:   FIND_RING_TYPE_SSSR,
        // RDKit❗✔️:   FIND_RING_TYPE_SYMM_SSSR,
        // RDKit❗✔️:   FIND_RING_TYPE_OTHER_OR_UNKNOWN
        // RDKit❗✔️: } FIND_RING_TYPE;
        // Canonical CORE computations remain here. Explicit field transport
        // moves actual vectors, including membership order, preallocated rows,
        // fused caches and represented URF count, without rebuilding any data.
        // The source URF object remains an independent unmodeled capability.
        // Constant-time moves, no cache recomputation or payload clone.
        cosmolkit_model::SourceRingInfo {
            initialized: self.initialized,
            find_type: match self.find_type {
                RingFindType::Fast => 0,
                RingFindType::Sssr => 1,
                RingFindType::SymmSssr => 2,
                RingFindType::OtherOrUnknown => 3,
            },
            atom_members: self.atom_members,
            bond_members: self.bond_members,
            atom_rings: self.atom_rings,
            bond_rings: self.bond_rings,
            atom_ring_families: self.atom_ring_families,
            bond_ring_families: self.bond_ring_families,
            relevant_cycle_count: self.relevant_cycle_count,
            fused_rings: self.fused_rings,
            num_fused_bonds: self.num_fused_bonds,
        }
    }

    /// Copy the actual detached cache, preserving initialized/empty/stale state.
    #[doc(hidden)]
    pub fn from_source_snapshot(
        source: &cosmolkit_model::SourceRingInfo,
    ) -> Result<Self, RingFindingError> {
        // RDKit❗✔️:   RingInfo(const RingInfo &other) = default;
        // RDKit❗✔️:   bool df_init{false};
        // RDKit❗✔️:   FIND_RING_TYPE df_find_type_type{FIND_RING_TYPE_OTHER_OR_UNKNOWN};
        // RDKit❗✔️:   DataType d_atomMembers, d_bondMembers;
        // RDKit❗✔️:   VECT_INT_VECT d_atomRings, d_bondRings;
        // RDKit❗✔️:   VECT_INT_VECT d_atomRingFamilies, d_bondRingFamilies;
        // RDKit❗✔️:   std::vector<boost::dynamic_bitset<>> d_fusedRings;
        // RDKit❗✔️:   std::vector<unsigned int> d_numFusedBonds;
        // Do not run from_persisted_components: it reconstructs membership
        // rows and cannot substitute for this source default copy. No ring
        // finding, eager cache repair or atom/bond dictionary conversion.
        // Linear deep vector copy, like the represented Native copy fields.
        let find_type = match source.find_type {
            0 => RingFindType::Fast,
            1 => RingFindType::Sssr,
            2 => RingFindType::SymmSssr,
            3 => RingFindType::OtherOrUnknown,
            tag => return Err(RingFindingError::SourceRingFindType { tag }),
        };
        Ok(Self {
            initialized: source.initialized,
            find_type,
            atom_members: source.atom_members.clone(),
            bond_members: source.bond_members.clone(),
            atom_rings: source.atom_rings.clone(),
            bond_rings: source.bond_rings.clone(),
            atom_ring_families: source.atom_ring_families.clone(),
            bond_ring_families: source.bond_ring_families.clone(),
            relevant_cycle_count: source.relevant_cycle_count,
            fused_rings: source.fused_rings.clone(),
            num_fused_bonds: source.num_fused_bonds.clone(),
        })
    }

    #[must_use]
    pub fn new(find_type: RingFindType, atom_count: usize, bond_count: usize) -> Self {
        let mut info = Self {
            initialized: false,
            find_type,
            atom_members: vec![Vec::new(); atom_count],
            bond_members: vec![Vec::new(); bond_count],
            atom_rings: Vec::new(),
            bond_rings: Vec::new(),
            atom_ring_families: Vec::new(),
            bond_ring_families: Vec::new(),
            relevant_cycle_count: None,
            fused_rings: Vec::new(),
            num_fused_bonds: Vec::new(),
        };
        info.initialize(find_type);
        info
    }

    #[must_use]
    pub const fn is_initialized(&self) -> bool {
        self.initialized
    }

    /// Number of atom membership rows. This is structural size information,
    /// not evidence of runtime cache validity or chemical correspondence.
    pub fn atom_row_count(&self) -> usize {
        self.atom_members.len()
    }

    /// Number of bond membership rows. This is structural size information,
    /// not evidence of runtime cache validity or chemical correspondence.
    pub fn bond_row_count(&self) -> usize {
        self.bond_members.len()
    }

    #[doc(hidden)]
    pub const fn persisted_find_type(&self) -> RingFindType {
        self.find_type
    }

    #[doc(hidden)]
    pub const fn persisted_relevant_cycle_count(&self) -> Option<usize> {
        self.relevant_cycle_count
    }

    #[doc(hidden)]
    pub fn persisted_fused_rings(&self) -> &[Vec<bool>] {
        &self.fused_rings
    }

    #[doc(hidden)]
    pub fn persisted_num_fused_bonds(&self) -> &[usize] {
        &self.num_fused_bonds
    }

    #[doc(hidden)]
    pub fn from_persisted_components(
        initialized: bool,
        find_type: RingFindType,
        atom_count: usize,
        bond_count: usize,
        atom_rings: Vec<Vec<AtomId>>,
        bond_rings: Vec<Vec<BondId>>,
        atom_ring_families: Vec<Vec<AtomId>>,
        bond_ring_families: Vec<Vec<BondId>>,
        relevant_cycle_count: Option<usize>,
        fused_rings: Vec<Vec<bool>>,
        num_fused_bonds: Vec<usize>,
    ) -> Result<Self, &'static str> {
        if !initialized {
            if find_type != RingFindType::OtherOrUnknown
                || !atom_rings.is_empty()
                || !bond_rings.is_empty()
                || !atom_ring_families.is_empty()
                || !bond_ring_families.is_empty()
                || relevant_cycle_count.is_some()
                || !fused_rings.is_empty()
                || !num_fused_bonds.is_empty()
            {
                return Err("uninitialized ring state contains materialized data");
            }
            return Ok(Self {
                initialized: false,
                find_type,
                atom_members: Vec::new(),
                bond_members: Vec::new(),
                atom_rings,
                bond_rings,
                atom_ring_families,
                bond_ring_families,
                relevant_cycle_count,
                fused_rings,
                num_fused_bonds,
            });
        }

        if atom_rings.len() != bond_rings.len() {
            return Err("ring atom/bond table length mismatch");
        }
        if atom_ring_families.len() != bond_ring_families.len() {
            return Err("ring-family atom/bond table length mismatch");
        }

        let mut atom_members = vec![Vec::new(); atom_count];
        let mut bond_members = vec![Vec::new(); bond_count];
        for (ring_index, (ring_atoms, ring_bonds)) in atom_rings.iter().zip(&bond_rings).enumerate()
        {
            if ring_atoms.len() != ring_bonds.len() {
                return Err("ring atom/bond size mismatch");
            }
            for atom in ring_atoms {
                let members = atom_members
                    .get_mut(atom.index())
                    .ok_or("ring atom index out of range")?;
                members.push(ring_index);
            }
            for bond in ring_bonds {
                let members = bond_members
                    .get_mut(bond.index())
                    .ok_or("ring bond index out of range")?;
                members.push(ring_index);
            }
        }
        for (family_atoms, family_bonds) in atom_ring_families.iter().zip(&bond_ring_families) {
            if family_atoms.iter().any(|atom| atom.index() >= atom_count) {
                return Err("ring-family atom index out of range");
            }
            if family_bonds.iter().any(|bond| bond.index() >= bond_count) {
                return Err("ring-family bond index out of range");
            }
        }

        let ring_count = atom_rings.len();
        if !fused_rings.is_empty() {
            if fused_rings.len() != ring_count
                || fused_rings.iter().any(|row| row.len() != ring_count)
            {
                return Err("fused-ring matrix dimensions do not match ring count");
            }
            for left in 0..ring_count {
                if fused_rings[left][left] {
                    return Err("fused-ring matrix diagonal must be false");
                }
                for right in left + 1..ring_count {
                    if fused_rings[left][right] != fused_rings[right][left] {
                        return Err("fused-ring matrix must be symmetric");
                    }
                }
            }
        }
        if !num_fused_bonds.is_empty() && num_fused_bonds.len() != ring_count {
            return Err("fused-bond count length does not match ring count");
        }
        if num_fused_bonds
            .iter()
            .zip(&bond_rings)
            .any(|(count, ring)| *count > ring.len())
        {
            return Err("fused-bond count exceeds ring size");
        }

        Ok(Self {
            initialized,
            find_type,
            atom_members,
            bond_members,
            atom_rings,
            bond_rings,
            atom_ring_families,
            bond_ring_families,
            relevant_cycle_count,
            fused_rings,
            num_fused_bonds,
        })
    }

    #[doc(hidden)]
    pub fn initialize(&mut self, find_type: RingFindType) {
        // BEGIN RDKIT CPP FUNCTION RingInfo::initialize
        // RDKit✔️✔️: void RingInfo::initialize(RDKit::FIND_RING_TYPE ringType) {
        // RDKit✔️✔️:   df_init = true;
        self.initialized = true;
        // RDKit✔️✔️:   df_find_type_type = ringType;
        self.find_type = find_type;
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION RingInfo::initialize
    }

    #[doc(hidden)]
    pub fn reset(&mut self) {
        // BEGIN RDKIT CPP FUNCTION RingInfo::reset
        // RDKit✔️✔️: void RingInfo::reset() {
        // RDKit✔️✔️:   if (!df_init) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        if !self.initialized {
            return;
        }
        // RDKit✔️✔️:   df_init = false;
        self.initialized = false;
        // RDKit✔️✔️:   df_find_type_type = RDKit::FIND_RING_TYPE_OTHER_OR_UNKNOWN;
        self.find_type = RingFindType::OtherOrUnknown;
        // RDKit✔️✔️:   d_atomMembers.clear();
        // RDKit✔️✔️:   d_bondMembers.clear();
        // RDKit✔️✔️:   d_atomRings.clear();
        // RDKit✔️✔️:   d_bondRings.clear();
        // RDKit✔️✔️:   d_atomRingFamilies.clear();
        // RDKit✔️✔️:   d_bondRingFamilies.clear();
        self.atom_members.clear();
        self.bond_members.clear();
        self.atom_rings.clear();
        self.bond_rings.clear();
        self.atom_ring_families.clear();
        self.bond_ring_families.clear();
        self.relevant_cycle_count = None;
        self.fused_rings.clear();
        self.num_fused_bonds.clear();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::reset
    }

    #[allow(dead_code)]
    pub fn preallocate(&mut self, atom_count: usize, bond_count: usize) {
        // BEGIN RDKIT CPP FUNCTION RingInfo::preallocate
        // RDKit✔️✔️: void RingInfo::preallocate(unsigned int numAtoms, unsigned int numBonds) {
        // RDKit✔️✔️:   d_atomMembers.resize(numAtoms);
        self.atom_members.resize(atom_count, Vec::new());
        // RDKit✔️✔️:   d_bondMembers.resize(numBonds);
        self.bond_members.resize(bond_count, Vec::new());
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::preallocate
    }

    #[must_use]
    pub const fn find_type(&self) -> RingFindType {
        self.find_type
    }

    #[must_use]
    pub fn is_find_fast_or_better(&self) -> bool {
        self.initialized
            && matches!(
                self.find_type,
                RingFindType::Fast | RingFindType::Sssr | RingFindType::SymmSssr
            )
    }

    #[must_use]
    pub fn is_sssr_or_better(&self) -> bool {
        self.initialized && matches!(self.find_type, RingFindType::Sssr | RingFindType::SymmSssr)
    }

    #[must_use]
    pub fn is_symm_sssr(&self) -> bool {
        self.initialized && self.find_type == RingFindType::SymmSssr
    }

    #[must_use]
    pub fn num_rings(&self) -> usize {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numRings
        // RDKit✔️✔️: unsigned int RingInfo::numRings() const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        // RDKit✔️✔️:   PRECONDITION(d_atomRings.size() == d_bondRings.size(), "length mismatch");
        debug_assert!(self.initialized, "RingInfo not initialized");
        debug_assert_eq!(
            self.atom_rings.len(),
            self.bond_rings.len(),
            "length mismatch"
        );
        // RDKit✔️✔️:   return rdcast<unsigned int>(d_atomRings.size());
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::numRings
        self.atom_rings.len()
    }

    #[must_use]
    pub fn atom_rings(&self) -> &[Vec<AtomId>] {
        &self.atom_rings
    }

    #[must_use]
    pub fn bond_rings(&self) -> &[Vec<BondId>] {
        &self.bond_rings
    }

    #[must_use]
    pub fn atom_ring_sizes(&self, atom: AtomId) -> Vec<usize> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::atomRingSizes
        // RDKit✔️✔️: RingInfo::INT_VECT RingInfo::atomRingSizes(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_atomMembers.size()) {
        // RDKit✔️✔️:     INT_VECT res(d_atomMembers[idx].size());
        // RDKit✔️✔️:     std::transform(d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
        // RDKit✔️✔️:                    res.begin(),
        // RDKit✔️✔️:                    [this](int ri) { return d_atomRings.at(ri).size(); });
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return INT_VECT();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::atomRingSizes
        self.atom_members
            .get(atom.index())
            .map(|members| {
                members
                    .iter()
                    .map(|&ring_idx| self.atom_rings[ring_idx].len())
                    .collect()
            })
            .unwrap_or_default()
    }

    #[must_use]
    pub fn bond_ring_sizes(&self, bond: BondId) -> Vec<usize> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::bondRingSizes
        // RDKit✔️✔️: RingInfo::INT_VECT RingInfo::bondRingSizes(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_bondMembers.size()) {
        // RDKit✔️✔️:     INT_VECT res(d_bondMembers[idx].size());
        // RDKit✔️✔️:     std::transform(d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
        // RDKit✔️✔️:                    res.begin(),
        // RDKit✔️✔️:                    [this](int ri) { return d_bondRings.at(ri).size(); });
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return INT_VECT();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::bondRingSizes
        self.bond_members
            .get(bond.index())
            .map(|members| {
                members
                    .iter()
                    .map(|&ring_idx| self.bond_rings[ring_idx].len())
                    .collect()
            })
            .unwrap_or_default()
    }

    #[must_use]
    pub fn is_atom_in_ring_of_size(&self, atom: AtomId, size: usize) -> bool {
        // BEGIN RDKIT CPP FUNCTION RingInfo::isAtomInRingOfSize
        // RDKit✔️✔️: bool RingInfo::isAtomInRingOfSize(unsigned int idx, unsigned int size) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_atomMembers.size()) {
        // RDKit✔️✔️:     return std::find_if(d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
        // RDKit✔️✔️:                         [this, size](int ri) {
        // RDKit✔️✔️:                           return d_atomRings.at(ri).size() == size;
        // RDKit✔️✔️:                         }) != d_atomMembers[idx].end();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::isAtomInRingOfSize
        self.atom_members
            .get(atom.index())
            .is_some_and(|members| members.iter().any(|&ri| self.atom_rings[ri].len() == size))
    }

    #[must_use]
    pub fn is_bond_in_ring_of_size(&self, bond: BondId, size: usize) -> bool {
        // BEGIN RDKIT CPP FUNCTION RingInfo::isBondInRingOfSize
        // RDKit✔️✔️: bool RingInfo::isBondInRingOfSize(unsigned int idx, unsigned int size) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_bondMembers.size()) {
        // RDKit✔️✔️:     return std::find_if(d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
        // RDKit✔️✔️:                         [this, size](int ri) {
        // RDKit✔️✔️:                           return d_bondRings.at(ri).size() == size;
        // RDKit✔️✔️:                         }) != d_bondMembers[idx].end();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::isBondInRingOfSize
        self.bond_members
            .get(bond.index())
            .is_some_and(|members| members.iter().any(|&ri| self.bond_rings[ri].len() == size))
    }

    #[must_use]
    pub fn min_atom_ring_size(&self, atom: AtomId) -> usize {
        // BEGIN RDKIT CPP FUNCTION RingInfo::minAtomRingSize
        // RDKit✔️✔️: unsigned int RingInfo::minAtomRingSize(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_atomMembers.size() && !d_atomMembers[idx].empty()) {
        // RDKit✔️✔️:     auto ri = *std::min_element(
        // RDKit✔️✔️:         d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
        // RDKit✔️✔️:         [this](int ri1, int ri2) {
        // RDKit✔️✔️:           return d_atomRings.at(ri1).size() < d_atomRings.at(ri2).size();
        // RDKit✔️✔️:         });
        // RDKit✔️✔️:     return d_atomRings.at(ri).size();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::minAtomRingSize
        self.atom_members
            .get(atom.index())
            .and_then(|members| members.iter().map(|&ri| self.atom_rings[ri].len()).min())
            .unwrap_or(0)
    }

    #[must_use]
    pub fn min_bond_ring_size(&self, bond: BondId) -> usize {
        // BEGIN RDKIT CPP FUNCTION RingInfo::minBondRingSize
        // RDKit✔️✔️: unsigned int RingInfo::minBondRingSize(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_bondMembers.size() && d_bondMembers[idx].size()) {
        // RDKit✔️✔️:     return d_bondRings
        // RDKit✔️✔️:         .at(*std::min_element(
        // RDKit✔️✔️:             d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
        // RDKit✔️✔️:             [this](int ri1, int ri2) {
        // RDKit✔️✔️:               return d_bondRings.at(ri1).size() < d_bondRings.at(ri2).size();
        // RDKit✔️✔️:             }))
        // RDKit✔️✔️:         .size();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::minBondRingSize
        self.bond_members
            .get(bond.index())
            .and_then(|members| members.iter().map(|&ri| self.bond_rings[ri].len()).min())
            .unwrap_or(0)
    }

    #[must_use]
    pub fn num_atom_rings(&self, atom: AtomId) -> usize {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numAtomRings
        // RDKit✔️✔️: unsigned int RingInfo::numAtomRings(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_atomMembers.size()) {
        // RDKit✔️✔️:     return rdcast<unsigned int>(d_atomMembers[idx].size());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::numAtomRings
        self.atom_members
            .get(atom.index())
            .map(Vec::len)
            .unwrap_or(0)
    }

    #[must_use]
    pub fn num_bond_rings(&self, bond: BondId) -> usize {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numBondRings
        // RDKit✔️✔️: unsigned int RingInfo::numBondRings(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx < d_bondMembers.size()) {
        // RDKit✔️✔️:     return rdcast<unsigned int>(d_bondMembers[idx].size());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::numBondRings
        self.bond_members
            .get(bond.index())
            .map(Vec::len)
            .unwrap_or(0)
    }

    #[must_use]
    pub fn atom_members(&self, atom: AtomId) -> &[usize] {
        // BEGIN RDKIT CPP FUNCTION RingInfo::atomMembers
        // RDKit✔️✔️: const RingInfo::INT_VECT &RingInfo::atomMembers(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   static const INT_VECT emptyVect;
        // RDKit✔️✔️:   if (idx < d_atomMembers.size()) {
        // RDKit✔️✔️:     return d_atomMembers[idx];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return emptyVect;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::atomMembers
        self.atom_members
            .get(atom.index())
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    #[must_use]
    pub fn bond_members(&self, bond: BondId) -> &[usize] {
        // BEGIN RDKIT CPP FUNCTION RingInfo::bondMembers
        // RDKit✔️✔️: const RingInfo::INT_VECT &RingInfo::bondMembers(unsigned int idx) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   static const INT_VECT emptyVect;
        // RDKit✔️✔️:   if (idx < d_bondMembers.size()) {
        // RDKit✔️✔️:     return d_bondMembers[idx];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return emptyVect;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::bondMembers
        self.bond_members
            .get(bond.index())
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    #[must_use]
    pub fn are_atoms_in_same_ring(&self, left: AtomId, right: AtomId) -> bool {
        self.are_atoms_in_same_ring_of_size(left, right, 0)
    }

    #[must_use]
    pub fn are_bonds_in_same_ring(&self, left: BondId, right: BondId) -> bool {
        self.are_bonds_in_same_ring_of_size(left, right, 0)
    }

    #[must_use]
    pub fn are_atoms_in_same_ring_of_size(&self, left: AtomId, right: AtomId, size: usize) -> bool {
        // BEGIN RDKIT CPP FUNCTION RingInfo::areAtomsInSameRingOfSize
        // RDKit✔️✔️: bool RingInfo::areAtomsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
        // RDKit✔️✔️:                                         unsigned int size) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx1 >= d_atomMembers.size() || idx2 >= d_atomMembers.size()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        let Some(left_members) = self.atom_members.get(left.index()) else {
            return false;
        };
        let Some(right_members) = self.atom_members.get(right.index()) else {
            return false;
        };
        let mut left_iter = left_members.iter();
        let mut right_iter = right_members.iter();
        let mut current_left = left_iter.next();
        let mut current_right = right_iter.next();
        while let (Some(&left_ring), Some(&right_ring)) = (current_left, current_right) {
            if left_ring < right_ring {
                current_left = left_iter.next();
            } else if left_ring > right_ring {
                current_right = right_iter.next();
            } else if size == 0 || self.atom_rings[left_ring].len() == size {
                return true;
            } else {
                current_left = left_iter.next();
                current_right = right_iter.next();
            }
        }
        false
    }

    #[must_use]
    pub fn are_bonds_in_same_ring_of_size(&self, left: BondId, right: BondId, size: usize) -> bool {
        // BEGIN RDKIT CPP FUNCTION RingInfo::areBondsInSameRingOfSize
        // RDKit✔️✔️: bool RingInfo::areBondsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
        // RDKit✔️✔️:                                         unsigned int size) const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   if (idx1 >= d_bondMembers.size() || idx2 >= d_bondMembers.size()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        let Some(left_members) = self.bond_members.get(left.index()) else {
            return false;
        };
        let Some(right_members) = self.bond_members.get(right.index()) else {
            return false;
        };
        let mut left_iter = left_members.iter();
        let mut right_iter = right_members.iter();
        let mut current_left = left_iter.next();
        let mut current_right = right_iter.next();
        while let (Some(&left_ring), Some(&right_ring)) = (current_left, current_right) {
            if left_ring < right_ring {
                current_left = left_iter.next();
            } else if left_ring > right_ring {
                current_right = right_iter.next();
            } else if size == 0 || self.bond_rings[left_ring].len() == size {
                return true;
            } else {
                current_left = left_iter.next();
                current_right = right_iter.next();
            }
        }
        false
    }

    pub fn is_ring_fused(&mut self, ring: usize) -> Result<bool, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::isRingFused
        // RDKit✔️✔️: bool RingInfo::isRingFused(unsigned int ringIdx) {
        // RDKit✔️✔️:   initFusedRings();
        self.init_fused_rings();
        // RDKit✔️✔️:   PRECONDITION(ringIdx < d_fusedRings.size(), "ringIdx out of bounds");
        if ring >= self.fused_rings.len() {
            return Err(RingFindingError::Value {
                message: "ringIdx out of bounds",
            });
        }
        // RDKit✔️✔️:   return d_fusedRings[ringIdx].any();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::isRingFused
        Ok(self.fused_rings[ring].iter().any(|fused| *fused))
    }

    pub fn are_rings_fused(
        &mut self,
        ring1: usize,
        ring2: usize,
    ) -> Result<bool, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::areRingsFused
        // RDKit✔️✔️: bool RingInfo::areRingsFused(unsigned int ring1Idx, unsigned int ring2Idx) {
        // RDKit✔️✔️:   initFusedRings();
        self.init_fused_rings();
        // RDKit✔️✔️:   PRECONDITION(ring1Idx < d_fusedRings.size(), "ring1Idx out of bounds");
        // RDKit✔️✔️:   PRECONDITION(ring2Idx < d_fusedRings.size(), "ring2Idx out of bounds");
        if ring1 >= self.fused_rings.len() || ring2 >= self.fused_rings.len() {
            return Err(RingFindingError::Value {
                message: "ringIdx out of bounds",
            });
        }
        // RDKit✔️✔️:   return d_fusedRings[ring1Idx].test(ring2Idx);
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::areRingsFused
        Ok(self.fused_rings[ring1][ring2])
    }

    pub fn num_fused_bonds(&mut self, ring: usize) -> Result<usize, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numFusedBonds
        // RDKit✔️✔️: unsigned int RingInfo::numFusedBonds(unsigned int ringIdx) {
        // RDKit✔️✔️:   PRECONDITION(ringIdx < d_bondRings.size(), "ringIdx out of bounds");
        if ring >= self.bond_rings.len() {
            return Err(RingFindingError::Value {
                message: "ringIdx out of bounds",
            });
        }
        // RDKit✔️✔️:   if (d_numFusedBonds.size() != d_bondRings.size()) {
        if self.num_fused_bonds.len() != self.bond_rings.len() {
            // RDKit✔️✔️:     d_numFusedBonds.clear();
            // RDKit✔️✔️:     d_numFusedBonds.resize(d_bondRings.size(), 0);
            self.num_fused_bonds.clear();
            self.num_fused_bonds.resize(self.bond_rings.len(), 0);
            // RDKit✔️✔️:     for (unsigned int ri = 0; ri < d_bondRings.size(); ++ri) {
            for ring_idx in 0..self.bond_rings.len() {
                // RDKit✔️✔️:       d_numFusedBonds[ri] += std::count_if(
                // RDKit✔️✔️:           d_bondRings[ri].begin(), d_bondRings[ri].end(),
                // RDKit✔️✔️:           [this](unsigned int bi) { return numBondRings(bi) > 1; });
                self.num_fused_bonds[ring_idx] = self.bond_rings[ring_idx]
                    .iter()
                    .filter(|bond| self.num_bond_rings(**bond) > 1)
                    .count();
            }
        }
        // RDKit✔️✔️:   return d_numFusedBonds[ringIdx];
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::numFusedBonds
        Ok(self.num_fused_bonds[ring])
    }

    pub fn num_fused_ring_neighbors(&mut self, ring: usize) -> Result<usize, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numFusedRingNeighbors
        // RDKit✔️✔️: unsigned int RingInfo::numFusedRingNeighbors(unsigned int ringIdx) {
        // RDKit✔️✔️:   initFusedRings();
        self.init_fused_rings();
        // RDKit✔️✔️:   PRECONDITION(ringIdx < d_fusedRings.size(), "ringIdx out of bounds");
        if ring >= self.fused_rings.len() {
            return Err(RingFindingError::Value {
                message: "ringIdx out of bounds",
            });
        }
        // RDKit✔️✔️:   return d_fusedRings[ringIdx].count();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::numFusedRingNeighbors
        Ok(self.fused_rings[ring]
            .iter()
            .filter(|fused| **fused)
            .count())
    }

    pub fn fused_ring_neighbors(&mut self, ring: usize) -> Result<Vec<usize>, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::fusedRingNeighbors
        // RDKit✔️✔️: std::vector<unsigned int> RingInfo::fusedRingNeighbors(unsigned int ringIdx) {
        // RDKit✔️✔️:   initFusedRings();
        self.init_fused_rings();
        // RDKit✔️✔️:   PRECONDITION(ringIdx < d_fusedRings.size(), "ringIdx out of bounds");
        if ring >= self.fused_rings.len() {
            return Err(RingFindingError::Value {
                message: "ringIdx out of bounds",
            });
        }
        // RDKit✔️✔️:   std::vector<unsigned int> res;
        // RDKit✔️✔️:   res.reserve(d_fusedRings[ringIdx].count());
        // RDKit✔️✔️:   for (unsigned int i = 0; i < d_fusedRings[ringIdx].size(); ++i) {
        // RDKit✔️✔️:     if (d_fusedRings[ringIdx].test(i)) {
        // RDKit✔️✔️:       res.push_back(i);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::fusedRingNeighbors
        Ok(self.fused_rings[ring]
            .iter()
            .enumerate()
            .filter_map(|(index, fused)| fused.then_some(index))
            .collect())
    }

    fn init_fused_rings(&mut self) {
        // BEGIN RDKIT CPP FUNCTION RingInfo::initFusedRings
        // RDKit✔️✔️: void RingInfo::initFusedRings() {
        // RDKit✔️✔️:   if (d_fusedRings.size() == d_bondRings.size()) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        if self.fused_rings.len() == self.bond_rings.len() {
            return;
        }
        // RDKit✔️✔️:   d_fusedRings.clear();
        self.fused_rings.clear();
        // RDKit✔️✔️:   if (d_bondRings.empty()) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        if self.bond_rings.is_empty() {
            return;
        }
        // RDKit✔️✔️:   d_fusedRings.resize(d_bondRings.size());
        // RDKit✔️✔️:   for (auto &fusedRing : d_fusedRings) {
        // RDKit✔️✔️:     fusedRing.resize(d_bondRings.size());
        // RDKit✔️✔️:   }
        self.fused_rings = vec![vec![false; self.bond_rings.len()]; self.bond_rings.len()];
        // RDKit✔️✔️:   for (const auto &ringIndices : d_bondMembers) {
        for ring_indices in &self.bond_members {
            // RDKit✔️✔️:     if (ringIndices.size() <= 1) {
            // RDKit✔️✔️:       continue;
            // RDKit✔️✔️:     }
            if ring_indices.len() <= 1 {
                continue;
            }
            // RDKit✔️✔️:     for (unsigned int i = 0; i < ringIndices.size() - 1; ++i) {
            for i in 0..ring_indices.len() - 1 {
                let ring1 = ring_indices[i];
                // RDKit✔️✔️:       for (unsigned int j = i + 1; j < ringIndices.size(); ++j) {
                for &ring2 in &ring_indices[i + 1..] {
                    // RDKit✔️✔️:         d_fusedRings[ringIdx1].set(ringIdx2);
                    // RDKit✔️✔️:         d_fusedRings[ringIdx2].set(ringIdx1);
                    self.fused_rings[ring1][ring2] = true;
                    self.fused_rings[ring2][ring1] = true;
                }
            }
        }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::initFusedRings
    }

    pub fn add_ring_family(
        &mut self,
        atom_indices: &[usize],
        bond_indices: &[usize],
    ) -> Result<usize, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::addRingFamily
        // RDKit✔️✔️: unsigned int RingInfo::addRingFamily(const INT_VECT &atomIndices,
        // RDKit✔️✔️:                                      const INT_VECT &bondIndices) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        if let Some(&atom) = atom_indices
            .iter()
            .find(|atom| **atom >= self.atom_members.len())
        {
            return Err(RingFindingError::RingAtomOutOfRange {
                atom,
                atom_count: self.atom_members.len(),
            });
        }
        if let Some(&bond) = bond_indices
            .iter()
            .find(|bond| **bond >= self.bond_members.len())
        {
            return Err(RingFindingError::RingBondOutOfRange {
                bond,
                bond_count: self.bond_members.len(),
            });
        }
        // RDKit✔️✔️:   d_atomRingFamilies.push_back(atomIndices);
        self.atom_ring_families
            .push(atom_indices.iter().copied().map(AtomId::new).collect());
        // RDKit✔️✔️:   d_bondRingFamilies.push_back(bondIndices);
        self.bond_ring_families
            .push(bond_indices.iter().copied().map(BondId::new).collect());
        // RDKit✔️✔️:   POSTCONDITION(d_atomRingFamilies.size() == d_bondRingFamilies.size(),
        // RDKit✔️✔️:                 "length mismatch");
        debug_assert_eq!(
            self.atom_ring_families.len(),
            self.bond_ring_families.len(),
            "length mismatch"
        );
        // RDKit✔️✔️:   return rdcast<unsigned int>(d_atomRingFamilies.size());
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RingInfo::addRingFamily
        Ok(self.atom_ring_families.len())
    }

    #[must_use]
    pub fn num_ring_families(&self) -> usize {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numRingFamilies
        // RDKit✔️✔️: unsigned int RingInfo::numRingFamilies() const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        debug_assert!(self.initialized, "RingInfo not initialized");
        // RDKit✔️✔️:   return d_atomRingFamilies.size();
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION RingInfo::numRingFamilies
        self.atom_ring_families.len()
    }

    #[must_use]
    pub const fn are_ring_families_initialized(&self) -> bool {
        // RDKit✔️✔️: bool areRingFamiliesInitialized() const { return dp_urfData != nullptr; }
        // The detached value retains the source observation without exposing
        // RDL-owned mutable data: `Some(0)` distinguishes a calculated acyclic
        // graph from a RingInfo on which family calculation never ran.
        self.relevant_cycle_count.is_some()
    }

    pub fn num_relevant_cycles(&self) -> Result<usize, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numRelevantCycles
        // RDKit✔️✔️: unsigned int RingInfo::numRelevantCycles() const {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        // RDKit✔️✔️:   return rdcast<unsigned int>(RDL_getNofRC(dp_urfData.get()));
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION RingInfo::numRelevantCycles
        debug_assert!(self.initialized, "RingInfo not initialized");
        // COSMolKit stores the URF relevant-cycle count only when
        // `find_ring_families()` materialized URF data. Calling this on a
        // RingInfo without that cache is an explicit unsupported branch.
        self.relevant_cycle_count
            .ok_or(RingFindingError::UnsupportedBranch {
                reason: "URF relevant-cycle data is not modeled",
            })
    }

    #[must_use]
    pub fn atom_ring_families(&self) -> &[Vec<AtomId>] {
        &self.atom_ring_families
    }

    #[must_use]
    pub fn bond_ring_families(&self) -> &[Vec<BondId>] {
        &self.bond_ring_families
    }

    /// Append already perceived, aligned atom and bond rows in the caller's
    /// index domain, retaining their order, memberships and supplied find type.
    /// The caller owns chemical cycle, completeness and perception-quality
    /// preconditions; this method neither finds nor certifies cycles or SSSR.
    ///
    /// A length mismatch returns an error before mutation. Membership tables
    /// grow for supplied indexes; rows are not sorted or deduplicated.
    /// Work and storage follow the supplied members and any membership growth,
    /// with the existing fused-ring cache clearing and typed row allocations.
    pub fn add_ring(
        &mut self,
        atom_indices: &[usize],
        bond_indices: &[usize],
    ) -> Result<usize, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::addRing
        // RDKit✔️✔️: unsigned int RingInfo::addRing(const INT_VECT &atomIndices,
        // RDKit✔️✔️:                                const INT_VECT &bondIndices) {
        // RDKit✔️✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        // RDKit✔️✔️:   PRECONDITION(atomIndices.size() == bondIndices.size(), "length mismatch");
        if atom_indices.len() != bond_indices.len() {
            return Err(RingFindingError::Value {
                message: "length mismatch",
            });
        }
        let ring_idx = self.atom_rings.len();
        self.fused_rings.clear();
        self.num_fused_bonds.clear();
        // RDKit✔️✔️:   for (const auto &i : atomIndices) {
        // RDKit✔️✔️:     if (i >= static_cast<int>(d_atomMembers.size())) {
        // RDKit✔️✔️:       d_atomMembers.resize(i + 1);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     d_atomMembers[i].push_back(d_atomRings.size());
        // RDKit✔️✔️:   }
        for &atom in atom_indices {
            if atom >= self.atom_members.len() {
                self.atom_members.resize(atom + 1, Vec::new());
            }
            self.atom_members[atom].push(ring_idx);
        }
        // RDKit✔️✔️:   for (const auto &i : bondIndices) {
        // RDKit✔️✔️:     if (i >= static_cast<int>(d_bondMembers.size())) {
        // RDKit✔️✔️:       d_bondMembers.resize(i + 1);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     d_bondMembers[i].push_back(d_bondRings.size());
        // RDKit✔️✔️:   }
        for &bond in bond_indices {
            if bond >= self.bond_members.len() {
                self.bond_members.resize(bond + 1, Vec::new());
            }
            self.bond_members[bond].push(ring_idx);
        }
        // RDKit✔️✔️:   d_atomRings.push_back(atomIndices);
        // RDKit✔️✔️:   d_bondRings.push_back(bondIndices);
        self.atom_rings
            .push(atom_indices.iter().copied().map(AtomId::new).collect());
        self.bond_rings
            .push(bond_indices.iter().copied().map(BondId::new).collect());
        // RDKit✔️✔️:   POSTCONDITION(d_atomRings.size() == d_bondRings.size(), "length mismatch");
        // RDKit✔️✔️:   return rdcast<unsigned int>(d_atomRings.size());
        // RDKit✔️✔️: }
        Ok(self.atom_rings.len())
    }
}

/// Materialize aligned, already perceived ring rows in their original index domain.
///
/// Callers must supply rows filtered from an existing `RingInfo`; row counts and
/// indices alone cannot establish that these are chemically valid cycles or a
/// complete SSSR set. The result is initialized even when both inputs are empty.
#[doc(hidden)]
pub fn ring_info_from_selected_rows(
    atom_count: usize,
    bond_count: usize,
    atom_rings: &[Vec<AtomId>],
    bond_rings: &[Vec<BondId>],
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::MolFragmentToSmiles ring transport
    // RDKit✔️❌:   // copy over the rings that only involve atoms/bonds in this fragment:
    // RDKit✔️❌:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️❌:     tmol.getRingInfo()->reset();
    // RDKit✔️❌:     tmol.getRingInfo()->initialize();
    // RDKit✔️❌:     for (unsigned int ridx = 0; ridx < mol.getRingInfo()->numRings(); ++ridx) {
    // RDKit✔️❌:       const INT_VECT &aring = mol.getRingInfo()->atomRings()[ridx];
    // RDKit✔️❌:       bool keepIt = true;
    // RDKit✔️❌:       for (auto aidx : aring) {
    // RDKit✔️❌:         if (!atomsInPlay[aidx]) {
    // RDKit✔️❌:           keepIt = false;
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (keepIt) {
    // RDKit✔️❌:         const INT_VECT &bring = mol.getRingInfo()->bondRings()[ridx];
    // RDKit✔️❌:         for (auto bidx : bring) {
    // RDKit✔️❌:           if (!bondsInPlay[bidx]) {
    // RDKit✔️❌:             keepIt = false;
    // RDKit✔️❌:             break;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (keepIt) {
    // RDKit✔️❌:           tmol.getRingInfo()->addRing(aring, bring);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION SmilesWrite::MolFragmentToSmiles ring transport
    // A validation pass makes all malformed input fail before any row is
    // installed. Construction then reuses the private source-shaped add_ring,
    // preserving row/member order and original IDs. Both passes are linear in
    // total members; the output has the same row and membership storage shape
    // as RingInfo::addRing. The typed-to-index Vecs and validation pass add
    // material temporary allocation/copy cost relative to the source.
    if atom_rings.len() != bond_rings.len() {
        return Err(RingFindingError::Value {
            message: "ring atom/bond table length mismatch",
        });
    }
    for (atoms, bonds) in atom_rings.iter().zip(bond_rings) {
        if atoms.len() != bonds.len() {
            return Err(RingFindingError::Value {
                message: "length mismatch",
            });
        }
        for atom in atoms {
            if atom.index() >= atom_count {
                return Err(RingFindingError::RingAtomOutOfRange {
                    atom: atom.index(),
                    atom_count,
                });
            }
        }
        for bond in bonds {
            if bond.index() >= bond_count {
                return Err(RingFindingError::RingBondOutOfRange {
                    bond: bond.index(),
                    bond_count,
                });
            }
        }
    }

    let mut info = RingInfo::new(RingFindType::OtherOrUnknown, atom_count, bond_count);
    for (atoms, bonds) in atom_rings.iter().zip(bond_rings) {
        let atom_indices: Vec<usize> = atoms.iter().map(|atom| atom.index()).collect();
        let bond_indices: Vec<usize> = bonds.iter().map(|bond| bond.index()).collect();
        info.add_ring(&atom_indices, &bond_indices)?;
    }
    Ok(info)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingSearchResult {
    find_type: RingFindType,
    rings: Vec<Vec<usize>>,
    extra_rings: Vec<Vec<usize>>,
    extra_rings_present: bool,
}

struct RingSearchContext<'a> {
    atom_count: usize,
    bonds: BondRows<'a>,
    adjacency: GraphAdjacency<'a>,
}

impl<'a> RingSearchContext<'a> {
    fn from_parts(atom_count: usize, bonds: &'a [Bond], adjacency: &'a AdjacencyList) -> Self {
        Self {
            atom_count,
            bonds: BondRows::Concrete(bonds),
            adjacency: GraphAdjacency::Concrete(adjacency),
        }
    }

    fn from_graph<G: StereoGraphAccess>(graph: &'a G) -> Self {
        Self {
            atom_count: graph.atoms().len(),
            bonds: graph.bonds(),
            adjacency: graph.adjacency(),
        }
    }
    fn atom_count(&self) -> usize {
        self.atom_count
    }

    fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    fn bonds(&self) -> BondRows<'a> {
        self.bonds
    }

    fn neighbors(&self, atom: usize) -> NeighborRows<'a> {
        self.adjacency.neighbors_of(atom)
    }

    fn bond_between_atoms(&self, begin: usize, end: usize) -> Option<BondId> {
        self.neighbors(begin)
            .iter()
            .find(|neighbor| neighbor.atom_index == end)
            .map(|neighbor| neighbor.bond)
    }
}

impl crate::paths::NeighborSource for RingSearchContext<'_> {
    fn atom_count(&self) -> usize {
        self.atom_count
    }

    fn neighbor_count(&self, atom: usize) -> usize {
        self.neighbors(atom).len()
    }

    fn neighbor_at(&self, atom: usize, position: usize) -> usize {
        self.neighbors(atom)
            .get(position)
            .expect("bounded source neighbor")
            .atom_index
    }
}

fn extra_ring_can_replace_sssr_ring(
    extra_ring: &[usize],
    ring: &[usize],
    bond_counts: &[i32],
) -> bool {
    // BEGIN RDKIT CPP INLINE BLOCK symmetrizeSSSR replacement predicate
    // RDKit✔️✔️:     for (auto &ring : bondsssrs) {
    // RDKit✔️✔️:       if (ring.size() != extraRing.size()) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // If `ring` is the only provider of some bond, extraRing must also
    // RDKit✔️✔️:       // provide that bond.
    // RDKit✔️✔️:       bool shareBond = false;
    // RDKit✔️✔️:       bool replacesAllUniqueBonds = true;
    // RDKit✔️✔️:       for (auto &bondID : ring) {
    // RDKit✔️✔️:         const int bondCount = bondCounts[bondID];
    // RDKit✔️✔️:         if (bondCount == 1 || !shareBond) {
    // RDKit✔️✔️:           auto position = find(extraRing.begin(), extraRing.end(), bondID);
    // RDKit✔️✔️:           if (position != extraRing.end()) {
    // RDKit✔️✔️:             shareBond = true;
    // RDKit✔️✔️:           } else if (bondCount == 1) {
    // RDKit✔️✔️:             // 1 means `ring` is the only ring in the SSSR to provide this
    // RDKit✔️✔️:             // bond, and extraRing did not provide it (so extraRing is not an
    // RDKit✔️✔️:             // acceptable substitution in the SSSR for ring)
    // RDKit✔️✔️:             replacesAllUniqueBonds = false;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (shareBond && replacesAllUniqueBonds) {
    // RDKit✔️✔️:         res.push_back(extraAtomRing);
    // RDKit✔️✔️:         FindRings::storeRingInfo(mol, extraAtomRing);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP INLINE BLOCK symmetrizeSSSR replacement predicate
    // A false result maps to the enclosing native candidate-loop continue;
    // all bond checks retain source order and its non-shortened scan.
    if ring.len() != extra_ring.len() {
        return false;
    }
    let mut share_bond = false;
    let mut replaces_all_unique_bonds = true;
    for &bond_id in ring {
        let bond_count = bond_counts[bond_id];
        if bond_count == 1 || !share_bond {
            if extra_ring.contains(&bond_id) {
                share_bond = true;
            } else if bond_count == 1 {
                replaces_all_unique_bonds = false;
            }
        }
    }
    share_bond && replaces_all_unique_bonds
}

/// Evaluate RDKit's bridgehead-atom query predicate over detached topology.
#[must_use]
pub fn is_atom_bridgehead_from_topology(
    topology: &cosmolkit_model::TopologyBlock,
    atom_idx: usize,
    ring_info: &RingInfo,
) -> i32 {
    // RDKit✔️✔️: if (at->getDegree() < 3) { return 0; }
    // RDKit✔️✔️: if (!ri || !ri->isInitialized()) { return 0; }
    // RDKit✔️✔️: for (const auto bnd : mol.atomBonds(at)) {
    // RDKit✔️✔️:   if (ri->numBondRings(bnd->getIdx())) {
    // RDKit✔️✔️:     atomRingBonds.set(bnd->getIdx());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (atomRingBonds.count() < 3) { return 0; }
    // RDKit✔️✔️: for (unsigned int i = 0; i < ri->bondRings().size(); ++i) {
    // RDKit✔️✔️:   for (unsigned int j = i + 1; j < ri->bondRings().size(); ++j) {
    // RDKit✔️✔️:     if (overlap >= 2 && atomInRingJ) {
    // RDKit✔️✔️:       ringsOverlap.set(i); ringsOverlap.set(j); break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!ringsOverlap[i]) { return 0; }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return 1;
    let adjacency = &topology.adjacency;
    if adjacency.neighbors_of(atom_idx).len() < 3 || !ring_info.is_initialized() {
        return 0;
    }

    let mut atom_ring_bonds = vec![false; topology.bonds.len()];
    for neighbor in adjacency.neighbors_of(atom_idx) {
        if ring_info.num_bond_rings(neighbor.bond) != 0 {
            atom_ring_bonds[neighbor.bond.index()] = true;
        }
    }
    if atom_ring_bonds.iter().filter(|is_ring| **is_ring).count() < 3 {
        return 0;
    }

    let bond_rings = ring_info.bond_rings();
    let mut bonds_in_ring_i = vec![false; topology.bonds.len()];
    let mut rings_overlap = vec![false; bond_rings.len()];
    for (i, ring_i) in bond_rings.iter().enumerate() {
        bonds_in_ring_i.fill(false);
        let mut atom_in_ring_i = false;
        for bond in ring_i {
            bonds_in_ring_i[bond.index()] = true;
            atom_in_ring_i |= atom_ring_bonds[bond.index()];
        }
        if !atom_in_ring_i {
            continue;
        }
        for (j, ring_j) in bond_rings.iter().enumerate().skip(i + 1) {
            let mut overlap = 0;
            let mut atom_in_ring_j = false;
            for bond in ring_j {
                atom_in_ring_j |= atom_ring_bonds[bond.index()];
                overlap += usize::from(bonds_in_ring_i[bond.index()]);
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

pub fn symmetrize_sssr_with_options_from_parts(
    atom_count: usize,
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::symmetrizeSSSR without output vector
    // RDKit✔️✔️: int symmetrizeSSSR(ROMol &mol, bool includeDativeBonds,
    // RDKit✔️✔️:                    bool includeHydrogenBonds) {
    // RDKit✔️✔️:   VECT_INT_VECT tmp;
    // RDKit✔️✔️:   return symmetrizeSSSR(mol, tmp, includeDativeBonds, includeHydrogenBonds);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION MolOps::symmetrizeSSSR without output vector
    // The same canonical context implementation prepares the returned ring
    // carrier and exposes its exact row count. No separate reaction body or
    // temporary mutable output matrix is needed by this detached boundary.
    // A graph-only detached value has the source's known-empty molecule
    // property carrier. This path never substitutes for a supplied carrier.
    let context = RingSearchContext::from_parts(atom_count, bonds, adjacency);
    symmetrize_sssr_from_context(&context, include_dative_bonds, include_hydrogen_bonds, None)
}

fn symmetrize_sssr_from_context(
    context: &RingSearchContext<'_>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
    mut properties: Option<&mut MoleculeProperties>,
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::symmetrizeSSSR complete source
    // RDKit✔️✔️: int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
    // RDKit✔️✔️:                    bool includeHydrogenBonds) {
    // RDKit✔️✔️:   res.clear();
    // RDKit✔️✔️:   VECT_INT_VECT sssrs;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // FIX: need to set flag here the symmetrization has been done in order to
    // RDKit✔️✔️:   // avoid repeating this work
    // RDKit✔️❌:   findSSSR(mol, sssrs, includeDativeBonds, includeHydrogenBonds);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // reinit as SYMM_SSSR
    // RDKit✔️✔️:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);
    // RDKit✔️✔️:
    // RDKit✔️🔝:   res.reserve(sssrs.size());
    // RDKit✔️🔝:   for (const auto &r : sssrs) {
    // RDKit✔️🔝:     res.emplace_back(r);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now check if there are any extra rings on the molecule
    // RDKit✔️✔️:   if (!mol.hasProp(common_properties::extraRings)) {
    // RDKit✔️✔️:     // no extra rings nothing to be done
    // RDKit✔️✔️:     return rdcast<int>(res.size());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const VECT_INT_VECT &extras =
    // RDKit✔️✔️:       mol.getProp<VECT_INT_VECT>(common_properties::extraRings);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // convert the rings to bond ids
    // RDKit✔️✔️:   VECT_INT_VECT bondsssrs;
    // RDKit✔️✔️:   RingUtils::convertToBonds(sssrs, bondsssrs, mol);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // For each "extra" ring, figure out if it could replace a single
    // RDKit✔️✔️:   // ring in the SSSR. A ring could be swapped out if:
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // * They are the same size
    // RDKit✔️✔️:   // * The replacement doesn't remove any bonds from the union of the bonds
    // RDKit✔️✔️:   //   in the SSSR.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // The latter can be checked by determining if the SSSR ring is the unique
    // RDKit✔️✔️:   // provider of any ring bond. If it is, the replacement ring must also
    // RDKit✔️✔️:   // provide that bond.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // May miss extra rings that would need to swap two (or three...) rings
    // RDKit✔️✔️:   // to be included.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // counts of each bond
    // RDKit✔️✔️:   std::vector<int> bondCounts(mol.getNumBonds(), 0);
    // RDKit✔️✔️:   for (const auto &r : bondsssrs) {
    // RDKit✔️✔️:     for (const auto &b : r) {
    // RDKit✔️✔️:       bondCounts[b] += 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_VECT extraRing;
    // RDKit✔️✔️:   for (auto &extraAtomRing : extras) {
    // RDKit✔️❌:     RingUtils::convertToBonds(extraAtomRing, extraRing, mol);
    // RDKit✔️✔️:     for (auto &ring : bondsssrs) {
    // RDKit✔️✔️:       if (ring.size() != extraRing.size()) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // If `ring` is the only provider of some bond, extraRing must also
    // RDKit✔️✔️:       // provide that bond.
    // RDKit✔️✔️:       bool shareBond = false;
    // RDKit✔️✔️:       bool replacesAllUniqueBonds = true;
    // RDKit✔️✔️:       for (auto &bondID : ring) {
    // RDKit✔️✔️:         const int bondCount = bondCounts[bondID];
    // RDKit✔️✔️:         if (bondCount == 1 || !shareBond) {
    // RDKit✔️✔️:           auto position = find(extraRing.begin(), extraRing.end(), bondID);
    // RDKit✔️✔️:           if (position != extraRing.end()) {
    // RDKit✔️✔️:             shareBond = true;
    // RDKit✔️✔️:           } else if (bondCount == 1) {
    // RDKit✔️✔️:             // 1 means `ring` is the only ring in the SSSR to provide this
    // RDKit✔️✔️:             // bond, and extraRing did not provide it (so extraRing is not an
    // RDKit✔️✔️:             // acceptable substitution in the SSSR for ring)
    // RDKit✔️✔️:             replacesAllUniqueBonds = false;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (shareBond && replacesAllUniqueBonds) {
    // RDKit✔️🔝:         res.push_back(extraAtomRing);
    // RDKit✔️✔️:         FindRings::storeRingInfo(mol, extraAtomRing);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (mol.hasProp(common_properties::extraRings)) {
    // RDKit✔️❌:     mol.clearProp(common_properties::extraRings);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return rdcast<int>(res.size());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::symmetrizeSSSR complete source

    let sssr = find_sssr_internal(
        context,
        include_dative_bonds,
        include_hydrogen_bonds,
        properties.as_deref_mut(),
        None,
        None,
    )?;
    let mut info = RingInfo::new(sssr.find_type, context.atom_count(), context.bond_count());
    // findSSSR stores its original rows before Symm reinitializes the type.
    // Reuse the canonical source conversion/storage, keeping memberships.
    store_rings_info(context, &sssr.rings, &mut info)?;
    info.initialize(RingFindType::SymmSssr);
    // This detached carrier already owns the exact output ring rows. Returning
    // it avoids the native separate res matrix copy without changing any row,
    // membership, order or count. No res alias is exposed by this boundary.
    if !sssr.extra_rings_present {
        return Ok(info);
    }
    {
        let bond_sssrs = convert_rings_to_bonds(&context, &sssr.rings)?;
        let mut bond_counts = vec![0i32; context.bond_count()];
        for ring in &bond_sssrs {
            for &bond in ring {
                bond_counts[bond] += 1;
            }
        }
        for extra_atom_ring in &sssr.extra_rings {
            let extra_ring =
                convert_to_bonds(&context, extra_atom_ring, |atom| *atom, BondId::index)?;
            for ring in &bond_sssrs {
                if extra_ring_can_replace_sssr_ring(&extra_ring, ring, &bond_counts) {
                    // Native stores each accepted extra immediately, before
                    // the next candidate and final computed-cache clear.
                    store_rings_info(context, std::slice::from_ref(extra_atom_ring), &mut info)?;
                    break;
                }
            }
        }
    }
    // Native conversion reuses extraRing capacity across iterations. The
    // current canonical converter returns a new vector per extra: retained
    // as a known allocation cost rather than claiming full perf equivalence.
    // Native symmetrizeSSSR erases its transient extraRings property after
    // using the typed cache. No other dictionary read/write intervenes.
    if sssr.extra_rings_present {
        if let Some(properties) = properties {
            properties.clear_prop("extraRings")?;
        }
    }
    Ok(info)
}

pub fn find_sssr_from_parts(
    atom_count: usize,
    bonds: &[Bond],
    adjacency: &AdjacencyList,
) -> Result<RingInfo, RingFindingError> {
    find_sssr_with_options_from_parts(atom_count, bonds, adjacency, false, false)
}

/// Native pointer-output dispatch over detached graph parts and real source
/// ring/property state. None means a local output matrix, never an empty cache
/// or replacement dictionary for a supplied source value.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn find_sssr_with_source_outputs_from_parts(
    atom_count: usize,
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    source_rings: &mut RingInfo,
    properties: Option<&mut MoleculeProperties>,
    output: Option<&mut Vec<Vec<usize>>>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<i32, RingFindingError> {
    // RDKit❗❌: int findSSSR(const ROMol &mol, VECT_INT_VECT *res, bool includeDativeBonds,
    // RDKit❗❌:              bool includeHydrogenBonds) {
    // RDKit❗❌:   if (!res) {
    // RDKit❗❌:     VECT_INT_VECT rings;
    // RDKit❗❌:     return findSSSR(mol, rings, includeDativeBonds, includeHydrogenBonds);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     return findSSSR(mol, (*res), includeDativeBonds, includeHydrogenBonds);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Both branches borrow the same reference kernel and preserve options,
    // cache/property mutation and output prefixes. The native signed INT_VECT
    // IDs use existing nonnegative usize graph IDs; allocator/pointer/build
    // width and generic Any extraRings value remain unmodeled capabilities.
    // The mirrored caller matrix costs an extra copy of completed ring rows;
    // no graph or dictionary is cloned. Native None has one local matrix.
    let context = RingSearchContext::from_parts(atom_count, bonds, adjacency);
    find_sssr_source_context(
        &context,
        source_rings,
        properties,
        output,
        include_dative_bonds,
        include_hydrogen_bonds,
    )
}

#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn find_sssr_with_source_outputs_from_graph<G: StereoGraphAccess>(
    graph: &G,
    source_rings: &mut RingInfo,
    properties: Option<&mut MoleculeProperties>,
    output: Option<&mut Vec<Vec<usize>>>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<i32, RingFindingError> {
    let context = RingSearchContext::from_graph(graph);
    find_sssr_source_context(
        &context,
        source_rings,
        properties,
        output,
        include_dative_bonds,
        include_hydrogen_bonds,
    )
}

fn find_sssr_source_context(
    context: &RingSearchContext<'_>,
    source_rings: &mut RingInfo,
    properties: Option<&mut MoleculeProperties>,
    output: Option<&mut Vec<Vec<usize>>>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<i32, RingFindingError> {
    if let Some(output) = output {
        find_sssr_reference_from_context(
            context,
            source_rings,
            properties,
            output,
            include_dative_bonds,
            include_hydrogen_bonds,
        )
    } else {
        let mut rings = Vec::new();
        find_sssr_reference_from_context(
            context,
            source_rings,
            properties,
            &mut rings,
            include_dative_bonds,
            include_hydrogen_bonds,
        )
    }
}

fn find_sssr_reference_from_context(
    context: &RingSearchContext<'_>,
    source_rings: &mut RingInfo,
    properties: Option<&mut MoleculeProperties>,
    output: &mut Vec<Vec<usize>>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<i32, RingFindingError> {
    // RDKit❗❌:   res.resize(0);
    // RDKit❗❌:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit❗❌:     mol.getRingInfo()->reset();
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SSSR);
    // RDKit❗❌:   FindRings::storeRingsInfo(mol, res);
    // RDKit❗❌:   return rdcast<int>(res.size());
    // The reference kernel's complete anchor lives in find_sssr_internal;
    // source output/cache state is supplied there, not inferred at success.
    output.clear();
    if source_rings.is_initialized() {
        source_rings.reset();
    }
    source_rings.initialize(RingFindType::Sssr);
    let result = find_sssr_internal(
        context,
        include_dative_bonds,
        include_hydrogen_bonds,
        properties,
        Some(output),
        Some(source_rings),
    )?;
    if result.find_type != RingFindType::Fast {
        store_rings_info(context, &result.rings, source_rings)?;
    }
    // Ordinary pinned rdcast is static_cast<int>, not checked/saturating.
    Ok(output.len() as i32)
}

pub fn find_sssr_with_options_from_parts(
    atom_count: usize,
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<RingInfo, RingFindingError> {
    let context = RingSearchContext::from_parts(atom_count, bonds, adjacency);
    find_sssr_from_context(&context, include_dative_bonds, include_hydrogen_bonds)
}

pub fn find_sssr(
    topology: &TopologyBlock,
    params: &RingSearchParams,
) -> Result<RingInfo, RingFindingError> {
    topology.validate()?;
    find_sssr_with_options_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        &topology.adjacency,
        params.include_dative_bonds,
        params.include_hydrogen_bonds,
    )
}

/// Source ring preparation consuming the caller's detached properties.
/// Returns both prepared rings and properties only on successful completion.
/// An error cannot expose a partially updated transient dictionary or cache.
#[doc(hidden)]
pub fn symmetrized_sssr_with_properties(
    topology: &TopologyBlock,
    mut properties: MoleculeProperties,
    params: &RingSearchParams,
) -> Result<(RingInfo, MoleculeProperties), RingFindingError> {
    // RDKit❗✔️:   findSSSR(mol, sssrs, includeDativeBonds, includeHydrogenBonds);
    // RDKit❗✔️:   mol.clearProp(common_properties::extraRings);
    // RDKit❗✔️:   if (mol.hasProp(common_properties::extraRings)) {
    // RDKit❗✔️:     mol.clearProp(common_properties::extraRings);
    // RDKit❗✔️:   }
    // The single source implementation below owns the actual clear order and
    // private cache. Native computed insertion retains its computed-list effect; the transient
    // key is erased before success, with no intervening generic property reads.
    // Ordered retained keys and source computed-list creation are reproduced
    // by representing that typed cache in RingSearchResult instead of adding
    // a second dictionary or a generic Any value implementation.
    topology.validate()?;
    let context =
        RingSearchContext::from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency);
    let info = symmetrize_sssr_from_context(
        &context,
        params.include_dative_bonds,
        params.include_hydrogen_bonds,
        Some(&mut properties),
    )?;
    Ok((info, properties))
}

pub fn symmetrized_sssr(
    topology: &TopologyBlock,
    params: &RingSearchParams,
) -> Result<RingInfo, RingFindingError> {
    topology.validate()?;
    symmetrize_sssr_with_options_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        &topology.adjacency,
        params.include_dative_bonds,
        params.include_hydrogen_bonds,
    )
}

fn find_sssr_from_context(
    context: &RingSearchContext<'_>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<RingInfo, RingFindingError> {
    let result = find_sssr_internal(
        context,
        include_dative_bonds,
        include_hydrogen_bonds,
        None,
        None,
        None,
    )?;
    let mut info = RingInfo::new(result.find_type, context.atom_count(), context.bond_count());
    store_rings_info(&context, &result.rings, &mut info)?;
    Ok(info)
}

pub fn fast_find_rings_from_parts(
    atom_count: usize,
    bonds: &[Bond],
    adjacency: &AdjacencyList,
) -> Result<RingInfo, RingFindingError> {
    let context = RingSearchContext::from_parts(atom_count, bonds, adjacency);
    fast_find_rings_from_context(&context)
}

pub fn fast_find_rings(topology: &TopologyBlock) -> Result<RingInfo, RingFindingError> {
    topology.validate()?;
    fast_find_rings_from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency)
}

fn fast_find_rings_from_context(
    context: &RingSearchContext<'_>,
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::fastFindRings complete detached result
    // RDKit✔️✔️: void fastFindRings(const ROMol &mol) {
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     mol.getRingInfo()->reset();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   mol.getRingInfo()->initialize(FIND_RING_TYPE_FAST);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   VECT_INT_VECT res;
    // RDKit✔️✔️:   res.resize(0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int nats = mol.getNumAtoms();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   INT_VECT atomColors(nats, 0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nats; ++i) {
    // RDKit✔️✔️:     if (atomColors[i]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (mol.getAtomWithIdx(i)->getDegree() < 2) {
    // RDKit✔️✔️:       atomColors[i] = 2;
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<const Atom *> traversalOrder;
    // RDKit✔️✔️:     _DFS(mol, mol.getAtomWithIdx(i), atomColors, traversalOrder, res);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   FindRings::storeRingsInfo(mol, res);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::fastFindRings complete detached result
    // A fresh returned RingInfo has the exact cleared atom/bond/family state
    // of native reset+FAST initialize. Initialize before the DFS phase, then
    // store cycles in discovery order in this one detached carrier. No old
    // RingInfo input or live cache mutation belongs to this CORE boundary.
    // Source coloring/traversal/cycle storage and bond lookup costs preserved;
    // no topology/atom/property clone or second ring perception algorithm.

    let mut info = RingInfo::new(
        RingFindType::Fast,
        context.atom_count(),
        context.bond_count(),
    );
    let rings = fast_find_rings_internal(&context)?;
    store_rings_info(&context, &rings, &mut info)?;
    Ok(info)
}

pub fn find_ring_families(
    topology: &TopologyBlock,
    params: &RingSearchParams,
) -> Result<RingInfo, RingFindingError> {
    topology.validate()?;
    find_ring_families_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        params.include_dative_bonds,
        params.include_hydrogen_bonds,
    )
}

pub fn find_ring_families_from_parts(
    atom_count: usize,
    bonds: &[Bond],
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::findRingFamilies
    // RDKit✔️✔️: #ifdef RDK_USE_URF
    // RDKit✔️✔️: void findRingFamilies(const ROMol &mol, bool includeDativeBonds,
    // RDKit✔️✔️:                       bool includeHydrogenBonds) {
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     // return if we've done this before
    // RDKit✔️✔️:     if (mol.getRingInfo()->areRingFamiliesInitialized()) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     mol.getRingInfo()->initialize();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDL_graph *graph = RDL_initNewGraph(mol.getNumAtoms());
    let mut graph = cosmolkit_ringdecomposer::Graph::new(atom_count);
    let mut graph_edge_to_bond_idx = Vec::with_capacity(bonds.len());
    // RDKit✔️✔️:   for (auto cbi : mol.bonds()) {
    for bond in bonds {
        // RDKit✔️✔️:     if (auto bt = cbi->getBondType();
        // RDKit✔️✔️:         bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
        // RDKit✔️✔️:         (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if bond.order() == BondOrder::Zero
            || (!include_dative_bonds && is_dative_bond_order(bond.order()))
            || (!include_hydrogen_bonds && bond.order() == BondOrder::Hydrogen)
        {
            continue;
        }
        // RDKit✔️✔️:     RDL_addUEdge(graph, cbi->getBeginAtomIdx(), cbi->getEndAtomIdx());
        let edge_id = graph.add_undirected_edge(bond.begin().index(), bond.end().index())?;
        if edge_id.index() != graph_edge_to_bond_idx.len() {
            return Err(RingFindingError::RingFamilyEdgeMappingMissing {
                edge: edge_id.index(),
            });
        }
        graph_edge_to_bond_idx.push(bond.id().index());
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDL_data *urfdata = RDL_calculate(graph);
    let decomposition = cosmolkit_ringdecomposer::RingDecomposition::calculate(graph)?;
    let relevant_cycle_count = checked_relevant_cycle_count(decomposition.relevant_cycle_count())?;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < RDL_getNofURF(urfdata); ++i) {
    let mut info = RingInfo::new(RingFindType::OtherOrUnknown, atom_count, bonds.len());
    info.relevant_cycle_count = Some(relevant_cycle_count);
    for urf in decomposition.urfs() {
        // RDKit✔️✔️:     RDL_node *nodes = nullptr;
        // RDKit✔️✔️:     unsigned nNodes = RDL_getNodesForURF(urfdata, i, &nodes);
        let atom_indices = urf
            .nodes()
            .iter()
            .map(|node| AtomId::new(*node).index())
            .collect::<Vec<_>>();
        // RDKit✔️✔️:     RDL_edge *edges = nullptr;
        // RDKit✔️✔️:     unsigned nEdges = RDL_getEdgesForURF(urfdata, i, &edges);
        let bond_indices =
            urf.edges()
                .iter()
                .map(|edge| {
                    graph_edge_to_bond_idx.get(edge.index()).copied().ok_or(
                        RingFindingError::RingFamilyEdgeMappingMissing { edge: edge.index() },
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
        // RDKit✔️✔️:     mol.getRingInfo()->addRingFamily(nvect, evect);
        info.add_ring_family(&atom_indices, &bond_indices)?;
    }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: #else
    // RDKit✔️✔️: void findRingFamilies(const ROMol &mol, bool includeDativeBonds,
    // RDKit✔️✔️:                       bool includeHydrogenBonds) {
    // RDKit✔️✔️:   BOOST_LOG(rdErrorLog)
    // RDKit✔️✔️:       << "This version of the RDKit was built without URF support" << std::endl;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: #endif
    // COSMolKit hard-requires URF-capable decomposition through
    // `cosmolkit_ringdecomposer`, so the non-URF compile-time branch is dead
    // in the pinned baseline and is represented by the absence of an alternate
    // build configuration rather than a runtime fallback path.
    // END RDKIT CPP FUNCTION MolOps::findRingFamilies
    Ok(info)
}

fn checked_relevant_cycle_count(count: f64) -> Result<usize, RingFindingError> {
    // RDKit's `rdcast<unsigned int>(RDL_getNofRC(...))` rejects values that
    // cannot be represented by the source return type. Rust float-to-integer
    // casts saturate or truncate, so reproduce the implicit source boundary
    // explicitly before conversion.
    if !count.is_finite() || count < 0.0 || count.fract() != 0.0 || count > f64::from(u32::MAX) {
        return Err(RingFindingError::RelevantCycleCountOutOfRange);
    }
    Ok(count as u32 as usize)
}

#[cfg(test)]
mod relevant_cycle_count_tests {
    use super::{RingFindingError, checked_relevant_cycle_count};

    #[test]
    fn checked_relevant_cycle_count_accepts_zero_and_source_maximum() {
        assert_eq!(checked_relevant_cycle_count(0.0), Ok(0));
        assert_eq!(
            checked_relevant_cycle_count(f64::from(u32::MAX)),
            Ok(u32::MAX as usize)
        );
    }

    #[test]
    fn checked_relevant_cycle_count_rejects_non_source_values() {
        for count in [
            -1.0,
            0.5,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            f64::from(u32::MAX) + 1.0,
        ] {
            assert_eq!(
                checked_relevant_cycle_count(count),
                Err(RingFindingError::RelevantCycleCountOutOfRange),
                "count {count:?} must not be truncated or saturated"
            );
        }
    }
}

fn is_dative_bond_order(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Dative | BondOrder::DativeOne | BondOrder::DativeLeft | BondOrder::DativeRight
    )
}

fn find_sssr_internal(
    context: &RingSearchContext<'_>,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
    mut properties: Option<&mut MoleculeProperties>,
    mut source_output: Option<&mut Vec<Vec<usize>>>,
    mut source_ring_state: Option<&mut RingInfo>,
) -> Result<RingSearchResult, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::findSSSR complete source
    // RDKit✔️✔️: int findSSSR(const ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
    // RDKit✔️✔️:              bool includeHydrogenBonds) {
    // RDKit✔️✔️:   res.resize(0);
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     mol.getRingInfo()->reset();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SSSR);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Zero-order bonds are not candidates for rings, and dative bonds and
    // RDKit✔️✔️:   // hydrogen bonds may also be out
    // RDKit✔️✔️:   const int nbnds = mol.getNumBonds();
    // RDKit✔️❌:   boost::dynamic_bitset<> activeBonds(nbnds);
    // RDKit✔️✔️:   activeBonds.set();
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (auto bt = bond->getBondType();
    // RDKit✔️✔️:         bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
    // RDKit✔️✔️:         (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
    // RDKit✔️✔️:       activeBonds[bond->getIdx()] = 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const unsigned int nats = mol.getNumAtoms();
    // RDKit✔️✔️:   INT_VECT atomDegrees(nats);
    // RDKit✔️✔️:   INT_VECT atomDegreesWithZeroOrderBonds(nats);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nats; ++i) {
    // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:     int deg = atom->getDegree();
    // RDKit✔️✔️:     atomDegrees[i] = deg;
    // RDKit✔️✔️:     atomDegreesWithZeroOrderBonds[i] = deg;
    // RDKit✔️✔️:     for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:       auto bt = bond->getBondType();
    // RDKit✔️✔️:       if (bt == Bond::ZERO || (!includeHydrogenBonds && bt == Bond::HYDROGEN) ||
    // RDKit✔️✔️:           (!includeDativeBonds && isDative(bt))) {
    // RDKit✔️✔️:         atomDegrees[i]--;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️❌:   mol.clearProp(common_properties::extraRings);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // find the number of fragments in the molecule - we will loop over them
    // RDKit✔️✔️:   RINGINVAR_SET invars;
    // RDKit✔️✔️:   INT_VECT curFrag;
    // RDKit✔️❌:   boost::dynamic_bitset<> ringAtoms(nats);
    // RDKit✔️❌:   boost::dynamic_bitset<> ringBonds(nbnds);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   VECT_INT_VECT frags;
    // RDKit✔️✔️:   getMolFrags(mol, frags);
    // RDKit✔️✔️:   // loop over the fragments in a molecule
    // RDKit✔️✔️:   for (const auto &curFrag : frags) {
    // RDKit✔️✔️:     if (curFrag.size() < 3) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // the following is the list of atoms that are useful in the next round of
    // RDKit✔️✔️:     // trimming basically atoms that become degree 0 or 1 because of bond
    // RDKit✔️✔️:     // removals initialized with atoms of degrees 0 and 1
    // RDKit✔️✔️:     std::queue<int> changed;
    // RDKit✔️✔️:     int bndcnt_with_zero_order_bonds = 0;
    // RDKit✔️✔️:     unsigned int nbnds = 0;
    // RDKit✔️✔️:     for (auto atom_idx : curFrag) {
    // RDKit✔️✔️:       bndcnt_with_zero_order_bonds += atomDegreesWithZeroOrderBonds[atom_idx];
    // RDKit✔️✔️:
    // RDKit✔️✔️:       int deg = atomDegrees[atom_idx];
    // RDKit✔️✔️:
    // RDKit✔️✔️:       nbnds += deg;
    // RDKit✔️✔️:       if (deg < 2) {
    // RDKit✔️✔️:         changed.push(atom_idx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // check to see if this fragment can even have a possible ring
    // RDKit✔️✔️:     CHECK_INVARIANT(bndcnt_with_zero_order_bonds % 2 == 0,
    // RDKit✔️✔️:                     "fragment graph has a dangling degree");
    // RDKit✔️✔️:     bndcnt_with_zero_order_bonds = bndcnt_with_zero_order_bonds / 2;
    // RDKit✔️✔️:     int num_possible_rings = bndcnt_with_zero_order_bonds - curFrag.size() + 1;
    // RDKit✔️✔️:     if (num_possible_rings < 1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     CHECK_INVARIANT(nbnds % 2 == 0,
    // RDKit✔️✔️:                     "fragment graph problem when including zero-order bonds");
    // RDKit✔️✔️:     nbnds = nbnds / 2;
    // RDKit✔️✔️:
    // RDKit✔️❌:     boost::dynamic_bitset<> doneAts(nats);
    // RDKit✔️✔️:     unsigned int nAtomsDone = 0;
    // RDKit✔️✔️:     VECT_INT_VECT fragRes;
    // RDKit✔️✔️:     while (nAtomsDone <= curFrag.size() - 3) {
    // RDKit✔️✔️:       // We can skip the 2 last atoms: if they were in a ring,
    // RDKit✔️✔️:       // we'd have already seen it.
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // trim all bonds that connect to degree 0 and 1 atoms
    // RDKit✔️✔️:       while (!changed.empty()) {
    // RDKit✔️✔️:         auto cand = changed.front();
    // RDKit✔️✔️:         changed.pop();
    // RDKit✔️✔️:         if (!doneAts[cand]) {
    // RDKit✔️✔️:           doneAts.set(cand);
    // RDKit✔️✔️:           ++nAtomsDone;
    // RDKit✔️✔️:           FindRings::trimBonds(cand, mol, changed, atomDegrees, activeBonds);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // all atoms left in the fragment should at least have a degree >= 2
    // RDKit✔️✔️:       // collect all the degree two nodes;
    // RDKit✔️✔️:       INT_VECT d2nodes;
    // RDKit✔️✔️:       FindRings::pickD2Nodes(mol, d2nodes, curFrag, atomDegrees, activeBonds);
    // RDKit✔️✔️:       if (d2nodes.size() > 0) {  // deal with the current degree two nodes
    // RDKit✔️✔️:         // place to record any duplicate rings discovered from the current d2
    // RDKit✔️✔️:         // nodes
    // RDKit✔️✔️:         FindRings::findRingsD2nodes(mol, fragRes, invars, d2nodes, atomDegrees,
    // RDKit✔️✔️:                                     activeBonds, ringBonds, ringAtoms);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // trim after we have dealt with all the current d2 nodes,
    // RDKit✔️✔️:         for (auto d2i : d2nodes) {
    // RDKit✔️✔️:           doneAts.set(d2i);
    // RDKit✔️✔️:           ++nAtomsDone;
    // RDKit✔️✔️:           FindRings::trimBonds(d2i, mol, changed, atomDegrees, activeBonds);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // end of degree two nodes
    // RDKit✔️✔️:       } else if (nAtomsDone <= curFrag.size() - 3) {
    // RDKit✔️✔️:         // now deal with higher degree nodes
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // this is brutal - we have no degree 2 nodes - find the first
    // RDKit✔️✔️:         // possible degree 3 node
    // RDKit✔️✔️:         int cand = -1;
    // RDKit✔️✔️:         for (auto aidi : curFrag) {
    // RDKit✔️✔️:           unsigned int deg = atomDegrees[aidi];
    // RDKit✔️✔️:           if (deg == 3) {
    // RDKit✔️✔️:             cand = (aidi);
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if we did not find a degree 3 node we are done
    // RDKit✔️✔️:         // REVIEW:
    // RDKit✔️✔️:         if (cand == -1) {
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         FindRings::findRingsD3Node(mol, fragRes, invars, cand, atomDegrees,
    // RDKit✔️✔️:                                    activeBonds);
    // RDKit✔️✔️:         doneAts.set(cand);
    // RDKit✔️✔️:         ++nAtomsDone;
    // RDKit✔️✔️:         FindRings::trimBonds(cand, mol, changed, atomDegrees, activeBonds);
    // RDKit✔️✔️:       }  // done with degree 3 node
    // RDKit✔️✔️:     }    // done finding rings in this fragment
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // calculate the cyclomatic number for the fragment:
    // RDKit✔️✔️:     int nexpt = rdcast<int>((nbnds - curFrag.size() + 1));
    // RDKit✔️✔️:     int ssiz = rdcast<int>(fragRes.size());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // first check that we got at least the number of expected rings
    // RDKit✔️✔️:     if (ssiz < nexpt) {
    // RDKit✔️✔️:       // Issue 3514824: in certain highly fused ring systems, the algorithm
    // RDKit✔️✔️:       // above would miss rings.
    // RDKit✔️✔️:       // for this fix to apply we have to have at least one non-ring bond
    // RDKit✔️✔️:       // that terminates in ring atoms. Find those bonds:
    // RDKit✔️✔️:       std::vector<const Bond *> possibleBonds;
    // RDKit✔️✔️:       for (unsigned int i = 0; i < nbnds; ++i) {
    // RDKit✔️✔️:         if (!ringBonds[i]) {
    // RDKit✔️✔️:           const Bond *bnd = mol.getBondWithIdx(i);
    // RDKit✔️✔️:           if (ringAtoms[bnd->getBeginAtomIdx()] &&
    // RDKit✔️✔️:               ringAtoms[bnd->getEndAtomIdx()]) {
    // RDKit✔️✔️:             possibleBonds.push_back(bnd);
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️❌:       boost::dynamic_bitset<> deadBonds(mol.getNumBonds());
    // RDKit✔️✔️:       while (possibleBonds.size()) {
    // RDKit✔️✔️:         bool ringFound = FindRings::findRingConnectingAtoms(
    // RDKit✔️✔️:             mol, possibleBonds[0], fragRes, invars, ringBonds, ringAtoms);
    // RDKit✔️✔️:         if (!ringFound) {
    // RDKit✔️✔️:           deadBonds.set(possibleBonds[0]->getIdx(), 1);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         possibleBonds.clear();
    // RDKit✔️✔️:         // check if we need to repeat the process:
    // RDKit✔️✔️:         for (unsigned int i = 0; i < nbnds; ++i) {
    // RDKit✔️✔️:           if (!ringBonds[i]) {
    // RDKit✔️✔️:             const Bond *bnd = mol.getBondWithIdx(i);
    // RDKit✔️✔️:             if (!deadBonds[bnd->getIdx()] &&
    // RDKit✔️✔️:                 ringAtoms[bnd->getBeginAtomIdx()] &&
    // RDKit✔️✔️:                 ringAtoms[bnd->getEndAtomIdx()]) {
    // RDKit✔️✔️:               possibleBonds.push_back(bnd);
    // RDKit✔️✔️:               break;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ssiz = rdcast<int>(fragRes.size());
    // RDKit✔️✔️:       if (ssiz < nexpt) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "WARNING: could not find number of expected rings. Switching to "
    // RDKit✔️✔️:                "an approximate ring finding algorithm."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         mol.getRingInfo()->reset();
    // RDKit✔️✔️:         fastFindRings(mol);
    // RDKit✔️✔️:         res.clear();
    // RDKit✔️✔️:         res = mol.getRingInfo()->atomRings();
    // RDKit✔️✔️:         return rdcast<int>(res.size());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // if we have more than expected we need to do some cleanup
    // RDKit✔️✔️:     // otherwise do som clean up work
    // RDKit✔️✔️:     if (ssiz > nexpt) {
    // RDKit✔️✔️:       FindRings::removeExtraRings(fragRes, nexpt, mol);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     res.insert(res.end(), fragRes.begin(), fragRes.end());
    // RDKit✔️✔️:   }  // done with all fragments
    // RDKit✔️✔️:
    // RDKit✔️✔️:   FindRings::storeRingsInfo(mol, res);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // update the ring memberships of atoms and bonds in the molecule:
    // RDKit✔️✔️:   // store the SSSR rings on the molecule as a property
    // RDKit✔️✔️:   // we will ignore any existing SSSRs on the molecule - simply overwrite
    // RDKit✔️✔️:   return rdcast<int>(res.size());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::findSSSR complete source
    // RingSearchResult owns the exact typed transient extraRings cache and
    // presence bit, including a present empty value. Supplied molecule
    // properties are cleared at the native point through their canonical
    // carrier; a graph-only value has known-empty properties, never a fallback
    // for a malformed supplied carrier. Symmetrization consumes and erases
    // this transient cache before returning; generic AnyTag remains unmodeled.
    // Canonical property clear/set bookkeeping uses the supplied carrier;
    // computed cache registration survives final clear as the source empty
    // computed-list value. Typed transient matrices remain owner-local.
    // Vec<bool> uses byte flags rather than source packed bitsets, a known
    // memory cost; cleanup BTreeSet operations have a known time/storage cost.
    // Source root/component order, inactive-edge degrees, degree2/3 trimming,
    // missing-ring search, fallback, extra-ring cleanup and append are retained.
    // No trace environment or non-source diagnostic branch remains.

    let mut res = Vec::new();
    let mut active_bonds = vec![true; context.bond_count()];
    for bond in context.bonds() {
        if is_inactive_ring_bond(bond.order(), include_dative_bonds, include_hydrogen_bonds) {
            active_bonds[bond.id().index()] = false;
        }
    }
    let mut atom_degrees = vec![0i32; context.atom_count()];
    let mut atom_degrees_with_zero_order_bonds = vec![0i32; context.atom_count()];
    for atom_idx in 0..context.atom_count() {
        let deg = i32::try_from(context.neighbors(atom_idx).len()).map_err(|_| {
            RingFindingError::Value {
                message: "atom degree out of range",
            }
        })?;
        atom_degrees[atom_idx] = deg;
        atom_degrees_with_zero_order_bonds[atom_idx] = deg;
        for neighbor in context.neighbors(atom_idx) {
            let bond = &context.bonds()[neighbor.bond.index()];
            if is_inactive_ring_bond(bond.order(), include_dative_bonds, include_hydrogen_bonds) {
                atom_degrees[atom_idx] -= 1;
            }
        }
    }
    // RDKit source clearProp occurs after active-bond/degree setup and before
    // any component/BFS processing; wrong computed-list tags propagate here.
    if let Some(properties) = properties.as_deref_mut() {
        properties.clear_prop("extraRings")?;
    }
    let mut extra_rings = Vec::new();
    let mut extra_rings_present = false;
    let mut invars = BTreeSet::new();
    let mut ring_atoms = vec![false; context.atom_count()];
    let mut ring_bonds = vec![false; context.bond_count()];
    let frags = molecule_fragments(context);
    for cur_frag in frags {
        if cur_frag.len() < 3 {
            continue;
        }
        let mut changed = VecDeque::new();
        let mut bndcnt_with_zero_order_bonds = 0i32;
        let mut nbnds = 0usize;
        for &atom_idx in &cur_frag {
            bndcnt_with_zero_order_bonds += atom_degrees_with_zero_order_bonds[atom_idx];
            let deg = atom_degrees[atom_idx];
            nbnds += usize::try_from(deg).map_err(|_| RingFindingError::Value {
                message: "active atom degree is negative",
            })?;
            if deg < 2 {
                changed.push_back(atom_idx);
            }
        }
        if bndcnt_with_zero_order_bonds % 2 != 0 {
            return Err(RingFindingError::Value {
                message: "fragment graph has a dangling degree",
            });
        }
        bndcnt_with_zero_order_bonds /= 2;
        let num_possible_rings = bndcnt_with_zero_order_bonds
            - i32::try_from(cur_frag.len()).map_err(|_| RingFindingError::Value {
                message: "fragment atom count out of range",
            })?
            + 1;
        if num_possible_rings < 1 {
            continue;
        }
        if nbnds % 2 != 0 {
            return Err(RingFindingError::Value {
                message: "fragment graph problem when including zero-order bonds",
            });
        }
        nbnds /= 2;
        let mut done_atoms = vec![false; context.atom_count()];
        let mut atoms_done = 0usize;
        let mut frag_res = Vec::new();
        while atoms_done <= (cur_frag.len() - 3) {
            while let Some(cand) = changed.pop_front() {
                if !done_atoms[cand] {
                    done_atoms[cand] = true;
                    atoms_done += 1;
                    trim_bonds(
                        context,
                        cand,
                        &mut changed,
                        &mut atom_degrees,
                        &mut active_bonds,
                    );
                }
            }
            let d2nodes = pick_d2_nodes(context, &cur_frag, &atom_degrees, &active_bonds);
            if !d2nodes.is_empty() {
                find_rings_d2_nodes(
                    context,
                    &mut frag_res,
                    &mut invars,
                    &d2nodes,
                    &mut atom_degrees,
                    &mut active_bonds,
                    &mut ring_bonds,
                    &mut ring_atoms,
                )?;
                for d2i in d2nodes {
                    done_atoms[d2i] = true;
                    atoms_done += 1;
                    trim_bonds(
                        context,
                        d2i,
                        &mut changed,
                        &mut atom_degrees,
                        &mut active_bonds,
                    );
                }
            } else if atoms_done <= (cur_frag.len() - 3) {
                let cand = cur_frag
                    .iter()
                    .copied()
                    .find(|&atom| atom_degrees[atom] == 3);
                let Some(cand) = cand else {
                    break;
                };
                find_rings_d3_node(context, &mut frag_res, &mut invars, cand, &active_bonds)?;
                done_atoms[cand] = true;
                atoms_done += 1;
                trim_bonds(
                    context,
                    cand,
                    &mut changed,
                    &mut atom_degrees,
                    &mut active_bonds,
                );
            }
        }
        // Native nbnds is unsigned int and curFrag.size() is size_t. Inactive bonds can
        // leave the full-topology fragment disconnected in the active graph,
        // so the subtraction intentionally wraps before the source cast
        // retains the same low 32 bits when cast to the source 32-bit `int`.
        let source_nbnds = u32::try_from(nbnds).map_err(|_| RingFindingError::Value {
            message: "fragment bond count out of range",
        })?;
        let source_fragment_size =
            u32::try_from(cur_frag.len()).map_err(|_| RingFindingError::Value {
                message: "fragment atom count out of range",
            })?;
        let source_nexpt = source_nbnds
            .wrapping_sub(source_fragment_size)
            .wrapping_add(1);
        let nexpt = (source_nexpt as i32) as isize;
        let mut ssize = isize::try_from(frag_res.len()).map_err(|_| RingFindingError::Value {
            message: "ring count out of range",
        })?;
        if ssize < nexpt {
            let mut dead_bonds = vec![false; context.bond_count()];
            while let Some(possible_bond) =
                next_possible_ring_bond(context, nbnds, &ring_bonds, &ring_atoms, &dead_bonds)
            {
                let ring_found = find_ring_connecting_atoms(
                    context,
                    possible_bond,
                    &mut frag_res,
                    &mut invars,
                    &mut ring_bonds,
                    &mut ring_atoms,
                )?;
                if !ring_found {
                    dead_bonds[possible_bond.index()] = true;
                }
            }
            ssize = isize::try_from(frag_res.len()).map_err(|_| RingFindingError::Value {
                message: "ring count out of range",
            })?;
            if ssize < nexpt {
                // This fallback and warning are present in pinned RDKit;
                // preserve them exactly, without adding a local ring heuristic.
                eprintln!(
                    "WARNING: could not find number of expected rings. Switching to an approximate ring finding algorithm."
                );
                if let Some(state) = source_ring_state.as_deref_mut() {
                    state.reset();
                    state.initialize(RingFindType::Fast);
                }
                let fast = fast_find_rings_internal(context)?;
                if let Some(state) = source_ring_state.as_deref_mut() {
                    store_rings_info(context, &fast, state)?;
                }
                // Native replaces res only after fastFindRings succeeds.
                if let Some(output) = source_output.as_deref_mut() {
                    output.clear();
                    output.extend(fast.iter().cloned());
                }
                return Ok(RingSearchResult {
                    find_type: RingFindType::Fast,
                    rings: fast,
                    // Native RingInfo reset/fastFindRings does not clear the
                    // separate extraRings cache accumulated by prior fragments.
                    extra_rings,
                    extra_rings_present,
                });
            }
        }
        if ssize > nexpt {
            let extras = remove_extra_rings(context, &mut frag_res)?;
            extra_rings.extend(extras);
            // Native removeExtraRings uses computed=true, including for an
            // empty extras matrix. Preserve actual computed-list creation and
            // first membership insertion while the typed cache stays local.
            if let Some(properties) = properties.as_deref_mut() {
                properties.register_transient_computed_name("extraRings")?;
            }
            extra_rings_present = true;
        }
        // Native caller res appends each completed component immediately;
        // later failures retain these rows, before final storeRingsInfo.
        if let Some(output) = source_output.as_deref_mut() {
            output.extend(frag_res.iter().cloned());
        }
        res.extend(frag_res);
    }
    Ok(RingSearchResult {
        find_type: RingFindType::Sssr,
        rings: res,
        extra_rings,
        extra_rings_present,
    })
}

fn next_possible_ring_bond(
    context: &RingSearchContext<'_>,
    source_bond_prefix: usize,
    ring_bonds: &[bool],
    ring_atoms: &[bool],
    dead_bonds: &[bool],
) -> Option<BondId> {
    // RDKit✔️✔️: for (unsigned int i = 0; i < nbnds; ++i) {
    // RDKit✔️✔️:   if (!ringBonds[i]) {
    // RDKit✔️✔️:     const Bond *bnd = mol.getBondWithIdx(i);
    context
        .bonds()
        .iter()
        .take(source_bond_prefix)
        .find_map(|bond| {
            (!ring_bonds[bond.id().index()]
                && !dead_bonds[bond.id().index()]
                && ring_atoms[bond.begin().index()]
                && ring_atoms[bond.end().index()])
            .then_some(bond.id())
        })
}

fn is_dative(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Dative | BondOrder::DativeOne | BondOrder::DativeLeft | BondOrder::DativeRight
    )
}

fn is_inactive_ring_bond(
    order: BondOrder,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> bool {
    order == BondOrder::Zero
        || (!include_dative_bonds && is_dative(order))
        || (!include_hydrogen_bonds && order == BondOrder::Hydrogen)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingInvariant {
    num_bits: usize,
    blocks: Vec<u64>,
}

impl Ord for RingInvariant {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // BEGIN BOOST CPP FUNCTION boost::operator<(dynamic_bitset, dynamic_bitset)
        // Boost✔️✔️: bool operator<(const dynamic_bitset<Block, Allocator>& a,
        // Boost✔️✔️:                const dynamic_bitset<Block, Allocator>& b)
        // Boost✔️✔️: {
        // Boost✔️✔️: //    assert(a.size() == b.size());
        // Boost✔️✔️:   typedef BOOST_DEDUCED_TYPENAME dynamic_bitset<Block, Allocator>::size_type size_type;
        // Boost✔️✔️:   size_type asize(a.size());
        // Boost✔️✔️:   size_type bsize(b.size());
        // Boost✔️✔️:   if (!bsize)
        // Boost✔️✔️:     {
        // Boost✔️✔️:     return false;
        // Boost✔️✔️:     }
        // Boost✔️✔️:   else if (!asize)
        // Boost✔️✔️:     {
        // Boost✔️✔️:     return true;
        // Boost✔️✔️:     }
        // Boost✔️✔️:   else if (asize == bsize)
        // Boost✔️✔️:     {
        // Boost✔️✔️:     for (size_type ii = a.num_blocks(); ii > 0; --ii)
        // Boost✔️✔️:       {
        // Boost✔️✔️:       size_type i = ii-1;
        // Boost✔️✔️:       if (a.m_bits[i] < b.m_bits[i])
        // Boost✔️✔️:         return true;
        // Boost✔️✔️:       else if (a.m_bits[i] > b.m_bits[i])
        // Boost✔️✔️:         return false;
        // Boost✔️✔️:       }
        // Boost✔️✔️:     return false;
        // Boost✔️✔️:     }
        // Boost✔️✔️:   else
        // Boost✔️✔️:     {
        // Boost✔️✔️:     size_type leqsize(std::min BOOST_PREVENT_MACRO_SUBSTITUTION(asize,bsize));
        // Boost✔️✔️:     for (size_type ii = 0; ii < leqsize; ++ii,--asize,--bsize)
        // Boost✔️✔️:       {
        // Boost✔️✔️:       size_type i = asize-1;
        // Boost✔️✔️:       size_type j = bsize-1;
        // Boost✔️✔️:       if (a[i] < b[j])
        // Boost✔️✔️:         return true;
        // Boost✔️✔️:       else if (a[i] > b[j])
        // Boost✔️✔️:         return false;
        // Boost✔️✔️:       }
        // Boost✔️✔️:     return (a.size() < b.size());
        // Boost✔️✔️:     }
        // Boost✔️✔️: }
        // END BOOST CPP FUNCTION boost::operator<(dynamic_bitset, dynamic_bitset)
        if self.num_bits == other.num_bits {
            return self.blocks.iter().rev().cmp(other.blocks.iter().rev());
        }

        let shared = self.num_bits.min(other.num_bits);
        for offset in 0..shared {
            let self_bit = self.bit(self.num_bits - offset - 1);
            let other_bit = other.bit(other.num_bits - offset - 1);
            match self_bit.cmp(&other_bit) {
                std::cmp::Ordering::Equal => {}
                ordering => return ordering,
            }
        }
        self.num_bits.cmp(&other.num_bits)
    }
}

impl PartialOrd for RingInvariant {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl RingInvariant {
    fn bit(&self, index: usize) -> bool {
        self.blocks
            .get(index / u64::BITS as usize)
            .is_some_and(|block| block & (1_u64 << (index % u64::BITS as usize)) != 0)
    }

    #[cfg(test)]
    fn atom_indices(&self) -> Vec<usize> {
        (0..self.num_bits)
            .filter(|&index| self.bit(index))
            .collect()
    }
}

fn compute_ring_invariant(ring: &[usize], num_atoms: usize) -> RingInvariant {
    // BEGIN RDKIT CPP FUNCTION RingUtils::computeRingInvariant
    // RDKit✔️✔️: RINGINVAR computeRingInvariant(const INT_VECT &ring, unsigned int numAtoms) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> res(numAtoms);
    let mut blocks = vec![0_u64; num_atoms.div_ceil(u64::BITS as usize)];
    // RDKit✔️✔️:   for (auto idx : ring) {
    // RDKit✔️✔️:     res.set(idx);
    // RDKit✔️✔️:   }
    for &atom_index in ring {
        blocks[atom_index / u64::BITS as usize] |= 1_u64 << (atom_index % u64::BITS as usize);
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RingUtils::computeRingInvariant
    RingInvariant {
        num_bits: num_atoms,
        blocks,
    }
}

fn convert_to_bonds<A, B>(
    context: &RingSearchContext<'_>,
    ring: &[A],
    atom_index: impl Fn(&A) -> usize,
    bond_index: impl Fn(BondId) -> B,
) -> Result<Vec<B>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::convertToBonds
    // RDKit✔️✔️: void convertToBonds(const INT_VECT &ring, INT_VECT &bondRing,
    // RDKit✔️✔️:                     const ROMol &mol) {
    // RDKit✔️✔️:   const auto rsiz = rdcast<unsigned int>(ring.size());
    // RDKit✔️✔️:   bondRing.resize(rsiz);
    let ring_size = ring.len();
    // One reserved output row; no input-index or output-ID staging vectors.
    let mut bond_ring = Vec::with_capacity(ring_size);
    if ring_size == 0 {
        return Ok(bond_ring);
    }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < (rsiz - 1); i++) {
    for i in 0..ring_size - 1 {
        // RDKit✔️✔️:     const Bond *bnd = mol.getBondBetweenAtoms(ring[i], ring[i + 1]);
        // RDKit✔️✔️:     if (!bnd) {
        // RDKit✔️✔️:       throw ValueErrorException("expected bond not found");
        // RDKit✔️✔️:     }
        let begin = atom_index(&ring[i]);
        let end = atom_index(&ring[i + 1]);
        let Some(bond) = context.bond_between_atoms(begin, end) else {
            return Err(RingFindingError::ExpectedBondNotFound {
                begin: AtomId::new(begin),
                end: AtomId::new(end),
            });
        };
        // RDKit✔️✔️:     bondRing[i] = bnd->getIdx();
        bond_ring.push(bond_index(bond));
    }
    // RDKit✔️✔️:   // bond from last to first atom
    // RDKit✔️✔️:   const Bond *bnd = mol.getBondBetweenAtoms(ring[rsiz - 1], ring[0]);
    // RDKit✔️✔️:   if (!bnd) {
    // RDKit✔️✔️:     throw ValueErrorException("expected bond not found");
    // RDKit✔️✔️:   }
    let closing_begin = atom_index(&ring[ring_size - 1]);
    let closing_end = atom_index(&ring[0]);
    let Some(bond) = context.bond_between_atoms(closing_begin, closing_end) else {
        return Err(RingFindingError::ExpectedBondNotFound {
            begin: AtomId::new(closing_begin),
            end: AtomId::new(closing_end),
        });
    };
    // RDKit✔️✔️:   bondRing[rsiz - 1] = bnd->getIdx();
    // RDKit✔️✔️: }
    bond_ring.push(bond_index(bond));
    // Behavior review (actual loop, both representations): each adjacent
    // and closing edge resolves through bond_between_atoms, which SCANS the
    // begin atom's neighbor rows (neighbors(begin).iter().find) until the
    // end index matches; a miss raises the typed ExpectedBondNotFound{begin,
    // end} cause the source raises as ValueErrorException("expected bond
    // not found"). The empty-row early return is this owner's explicit
    // zero-length boundary — a CK-added guard, NOT a claim of source
    // empty-subtraction equivalence (the source loops themselves are
    // entered only with rsiz-1 == UINT_MAX on a zero row; no claim is made
    // about upstream behavior there).
    // Complexity review (actual loop): ONE reserved output Vec with exactly
    // one push per row entry; borrowed inputs/context; per-edge cost is
    // O(degree(begin)) neighbor-row scanning (NOT O(1)), so the total is the
    // sum of visited neighbor-row scan costs over all consecutive pairs plus
    // the closing edge; no input-index staging vector, no output-ID
    // conversion vector, no second graph traversal, and no grow-only pushes.
    Ok(bond_ring)
}

/// Convert one borrowed typed atom ring to its graph-derived typed bond row.
///
/// The ONE conversion owner (RingUtils::convertToBonds) is reused through
/// scalar projections; the existing borrowed RingSearchContext is built in
/// O(1) and no staging or output-conversion vectors are introduced.
pub(crate) fn ring_atom_ids_to_bond_ids(
    topology: &TopologyBlock,
    ring: &[AtomId],
) -> Result<Vec<BondId>, RingFindingError> {
    let context =
        RingSearchContext::from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency);
    convert_to_bonds(&context, ring, |atom| atom.index(), |bond| bond)
}

fn convert_rings_to_bonds(
    context: &RingSearchContext<'_>,
    rings: &[Vec<usize>],
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    let mut bond_rings = Vec::with_capacity(rings.len());
    for ring in rings {
        bond_rings.push(convert_to_bonds(
            context,
            ring,
            |atom| *atom,
            BondId::index,
        )?);
    }
    Ok(bond_rings)
}

fn store_rings_info(
    context: &RingSearchContext<'_>,
    rings: &[Vec<usize>],
    info: &mut RingInfo,
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::storeRingsInfo and storeRingInfo
    // RDKit✔️✔️: void storeRingsInfo(const ROMol &mol, const VECT_INT_VECT &rings) {
    // RDKit✔️✔️:   for (const auto &ring : rings) {
    // RDKit✔️✔️:     storeRingInfo(mol, ring);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: void storeRingInfo(const ROMol &mol, const INT_VECT &ring) {
    // RDKit✔️✔️:   INT_VECT bondIndices;
    // RDKit✔️✔️:   RingUtils::convertToBonds(ring, bondIndices, mol);
    // RDKit✔️✔️:   mol.getRingInfo()->addRing(ring, bondIndices);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::storeRingsInfo and storeRingInfo

    for ring in rings {
        let bond_indices = convert_to_bonds(context, ring, |atom| *atom, BondId::index)?;
        info.add_ring(ring, &bond_indices)?;
    }
    Ok(())
}

fn trim_bonds(
    context: &RingSearchContext<'_>,
    cand: usize,
    changed: &mut VecDeque<usize>,
    atom_degrees: &mut [i32],
    active_bonds: &mut [bool],
) {
    // BEGIN RDKIT CPP FUNCTION FindRings::trimBonds
    // RDKit✔️✔️: void trimBonds(unsigned int cand, const ROMol &tMol, std::queue<int> &changed,
    // RDKit✔️✔️:                INT_VECT &atomDegrees, boost::dynamic_bitset<> &activeBonds) {
    // RDKit✔️✔️:   for (auto bond : tMol.atomBonds(tMol.getAtomWithIdx(cand))) {
    for neighbor in context.neighbors(cand) {
        let bond_index = neighbor.bond.index();
        // RDKit✔️✔️:     if (!activeBonds[bond->getIdx()]) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if !active_bonds[bond_index] {
            continue;
        }
        // RDKit✔️✔️:     unsigned int oIdx = bond->getOtherAtomIdx(cand);
        let other = neighbor.atom_index;
        // RDKit✔️✔️:     if (atomDegrees[oIdx] <= 2) {
        // RDKit✔️✔️:       changed.push(oIdx);
        // RDKit✔️✔️:     }
        if atom_degrees[other] <= 2 {
            changed.push_back(other);
        }
        // RDKit✔️✔️:     activeBonds[bond->getIdx()] = 0;
        active_bonds[bond_index] = false;
        // RDKit✔️✔️:     atomDegrees[oIdx] -= 1;
        // RDKit✔️✔️:     atomDegrees[cand] -= 1;
        atom_degrees[other] -= 1;
        atom_degrees[cand] -= 1;
    }
    // RDKit✔️✔️: }
}

fn smallest_rings_bfs(
    context: &RingSearchContext<'_>,
    root: usize,
    active_bonds: &[bool],
    forbidden: Option<&[usize]>,
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::smallestRingsBfs
    // RDKit✔️✔️: int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings,
    // RDKit✔️✔️:                      boost::dynamic_bitset<> &activeBonds,
    // RDKit✔️✔️:                      INT_VECT *forbidden = nullptr) {
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;
    let mut rings = Vec::new();
    // RDKit✔️✔️:   std::vector<int> done(mol.getNumAtoms(), WHITE);
    let mut done = vec![WHITE; context.atom_count()];
    // RDKit✔️✔️:   if (forbidden) {
    if let Some(forbidden) = forbidden {
        for &atom in forbidden {
            if let Some(done) = done.get_mut(atom) {
                *done = BLACK;
            }
        }
    }
    // RDKit✔️✔️:   std::vector<int> parents(mol.getNumAtoms(), -1);
    // RDKit✔️✔️:   std::vector<int> depths(mol.getNumAtoms(), 0);
    let mut parents = vec![None; context.atom_count()];
    let mut depths = vec![0usize; context.atom_count()];
    // RDKit✔️✔️:   std::deque<int> bfsq;
    // RDKit✔️✔️:   bfsq.push_back(root);
    let mut bfsq = VecDeque::new();
    bfsq.push_back(root);
    // RDKit✔️✔️:   unsigned int curSize = UINT_MAX;
    let mut cur_size = usize::MAX;
    // RDKit✔️✔️:   while (!bfsq.empty()) {
    while let Some(curr) = bfsq.pop_front() {
        if bfsq.len() >= MAX_BFSQ_SIZE {
            return Err(RingFindingError::Value {
                message: "Maximum BFS search size exceeded.\nThis is likely due to a highly symmetric fused ring system.",
            });
        }
        // RDKit✔️✔️:     done[curr] = BLACK;
        done[curr] = BLACK;
        // RDKit✔️✔️:     const unsigned int depth = depths[curr] + 1;
        let depth = depths[curr] + 1;
        // RDKit✔️✔️:     if (depth > curSize) {
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        if depth > cur_size {
            break;
        }
        // RDKit✔️✔️:     for (auto bond : mol.atomBonds(mol.getAtomWithIdx(curr))) {
        for neighbor in context.neighbors(curr) {
            // RDKit✔️✔️:       if (!activeBonds[bond->getIdx()]) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if !active_bonds[neighbor.bond.index()] {
                continue;
            }
            let neighbor_idx = neighbor.atom_index;
            // RDKit✔️✔️:       if (done[nbrIdx] == BLACK || parents[curr] == nbrIdx) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if done[neighbor_idx] == BLACK || parents[curr] == Some(neighbor_idx) {
                continue;
            }
            // RDKit✔️✔️:       if (done[nbrIdx] == WHITE) {
            if done[neighbor_idx] == WHITE {
                // RDKit✔️✔️:         parents[nbrIdx] = curr;
                // RDKit✔️✔️:         done[nbrIdx] = GRAY;
                // RDKit✔️✔️:         depths[nbrIdx] = depth;
                // RDKit✔️✔️:         bfsq.push_back(nbrIdx);
                parents[neighbor_idx] = Some(curr);
                done[neighbor_idx] = GRAY;
                depths[neighbor_idx] = depth;
                bfsq.push_back(neighbor_idx);
            } else {
                let mut ring = vec![neighbor_idx];
                let mut parent = parents[neighbor_idx];
                while let Some(parent_idx) = parent {
                    if parent_idx == root {
                        break;
                    }
                    ring.push(parent_idx);
                    parent = parents[parent_idx];
                }
                ring.insert(0, curr);
                parent = parents[curr];
                while let Some(parent_idx) = parent {
                    if ring.contains(&parent_idx) {
                        ring.clear();
                        break;
                    }
                    ring.insert(0, parent_idx);
                    parent = parents[parent_idx];
                }
                if ring.len() > 1 {
                    if ring.len() <= cur_size {
                        cur_size = ring.len();
                        rings.push(ring);
                    } else {
                        return Ok(rings);
                    }
                }
            }
        }
    }
    // RDKit✔️✔️:   return rdcast<unsigned int>(rings.size());
    // RDKit✔️✔️: }
    Ok(rings)
}

fn mark_useless_d2s(
    context: &RingSearchContext<'_>,
    root: usize,
    forbidden: &mut [bool],
    atom_degrees: &[i32],
    active_bonds: &[bool],
) {
    // BEGIN RDKIT CPP FUNCTION FindRings::markUselessD2s
    // RDKit✔️✔️: void markUselessD2s(unsigned int root, const ROMol &tMol,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> &forb, const INT_VECT &atomDegrees,
    // RDKit✔️✔️:                     const boost::dynamic_bitset<> &activeBonds) {
    // RDKit✔️✔️:   for (auto bond : tMol.atomBonds(tMol.getAtomWithIdx(root))) {
    for neighbor in context.neighbors(root) {
        // RDKit✔️✔️:     if (!activeBonds[bond->getIdx()]) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if !active_bonds[neighbor.bond.index()] {
            continue;
        }
        // RDKit✔️✔️:     unsigned int oIdx = bond->getOtherAtomIdx(root);
        let other = neighbor.atom_index;
        // RDKit✔️✔️:     if (!forb[oIdx] && atomDegrees[oIdx] == 2) {
        // RDKit✔️✔️:       forb[oIdx] = 1;
        // RDKit✔️✔️:       markUselessD2s(oIdx, tMol, forb, atomDegrees, activeBonds);
        // RDKit✔️✔️:     }
        if !forbidden[other] && atom_degrees[other] == 2 {
            forbidden[other] = true;
            mark_useless_d2s(context, other, forbidden, atom_degrees, active_bonds);
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::markUselessD2s
}

fn pick_d2_nodes(
    context: &RingSearchContext<'_>,
    current_fragment: &[usize],
    atom_degrees: &[i32],
    active_bonds: &[bool],
) -> Vec<usize> {
    // BEGIN RDKIT CPP FUNCTION FindRings::pickD2Nodes
    // RDKit✔️✔️: void pickD2Nodes(const ROMol &tMol, INT_VECT &d2nodes, const INT_VECT &currFrag,
    // RDKit✔️✔️:                  const INT_VECT &atomDegrees,
    // RDKit✔️✔️:                  const boost::dynamic_bitset<> &activeBonds) {
    // RDKit✔️✔️:   d2nodes.resize(0);
    let mut d2nodes = Vec::new();
    // RDKit✔️✔️:   boost::dynamic_bitset<> forb(tMol.getNumAtoms());
    let mut forbidden = vec![false; context.atom_count()];
    // RDKit✔️✔️:   while (1) {
    loop {
        // RDKit✔️✔️:     int root = -1;
        let mut root = None;
        // RDKit✔️✔️:     for (int axci : currFrag) {
        for &atom in current_fragment {
            // RDKit✔️✔️:       if (atomDegrees[axci] == 2 && !forb[axci]) {
            // RDKit✔️✔️:         root = axci;
            // RDKit✔️✔️:         d2nodes.push_back(axci);
            // RDKit✔️✔️:         forb[axci] = 1;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if atom_degrees[atom] == 2 && !forbidden[atom] {
                root = Some(atom);
                d2nodes.push(atom);
                forbidden[atom] = true;
                break;
            }
        }
        let Some(root) = root else {
            // RDKit✔️✔️:     if (root == -1) {
            // RDKit✔️✔️:       break;
            break;
        };
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       markUselessD2s(root, tMol, forb, atomDegrees, activeBonds);
        // RDKit✔️✔️:     }
        mark_useless_d2s(context, root, &mut forbidden, atom_degrees, active_bonds);
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::pickD2Nodes
    d2nodes
}

fn find_sssr_for_dup_candidates(
    context: &RingSearchContext<'_>,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvariant>,
    dup_map: &BTreeMap<usize, Vec<usize>>,
    dup_d2_candidates: &BTreeMap<RingInvariant, Vec<usize>>,
    atom_degrees: &[i32],
    active_bonds: &[bool],
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findSSSRforDupCands
    // RDKit✔️✔️: void findSSSRforDupCands(const ROMol &mol, VECT_INT_VECT &res,
    // RDKit✔️✔️:                          RINGINVAR_SET &invars, const INT_INT_VECT_MAP &dupMap,
    // RDKit✔️✔️:                          const RINGINVAR_INT_VECT_MAP &dupD2Cands,
    // RDKit✔️✔️:                          INT_VECT &atomDegrees,
    // RDKit✔️✔️:                          const boost::dynamic_bitset<> &activeBonds) {
    // RDKit✔️✔️:   for (const auto &dupD2Cand : dupD2Cands) {
    for dup_candidates in dup_d2_candidates.values() {
        // RDKit✔️✔️:     const INT_VECT &dupCands = dupD2Cand.second;
        // RDKit✔️✔️:     if (dupCands.size() > 1) {
        if dup_candidates.len() > 1 {
            // RDKit✔️✔️:       VECT_INT_VECT nrings;
            // RDKit✔️✔️:       auto minSiz = static_cast<unsigned int>(MAX_INT);
            let mut new_rings = Vec::new();
            let mut min_size = usize::MAX;
            // RDKit✔️✔️:       for (int dupCand : dupCands) {
            for &dup_candidate in dup_candidates {
                // RDKit✔️✔️:         INT_VECT atomDegreesCopy = atomDegrees;
                // RDKit✔️✔️:         boost::dynamic_bitset<> activeBondsCopy = activeBonds;
                // RDKit✔️✔️:         std::queue<int> changed;
                let mut atom_degrees_copy = atom_degrees.to_vec();
                let mut active_bonds_copy = active_bonds.to_vec();
                let mut changed = VecDeque::new();
                // RDKit✔️✔️:         auto dmci = dupMap.find(dupCand);
                // RDKit✔️✔️:         CHECK_INVARIANT(dmci != dupMap.end(), "duplicate could not be found");
                let Some(mapped_duplicates) = dup_map.get(&dup_candidate) else {
                    return Err(RingFindingError::Value {
                        message: "duplicate could not be found",
                    });
                };
                // RDKit✔️✔️:         for (int dni : dmci->second) {
                // RDKit✔️✔️:           trimBonds(dni, mol, changed, atomDegreesCopy, activeBondsCopy);
                // RDKit✔️✔️:         }
                for &duplicate_node in mapped_duplicates {
                    trim_bonds(
                        context,
                        duplicate_node,
                        &mut changed,
                        &mut atom_degrees_copy,
                        &mut active_bonds_copy,
                    );
                }
                // RDKit✔️✔️:         VECT_INT_VECT srings;
                // RDKit✔️✔️:         smallestRingsBfs(mol, dupCand, srings, activeBondsCopy);
                let smallest =
                    smallest_rings_bfs(context, dup_candidate, &active_bonds_copy, None)?;
                // RDKit✔️✔️:         nrings.reserve(srings.size());
                // RDKit✔️✔️:         for (const auto &sri : srings) {
                for ring in smallest {
                    // RDKit✔️✔️:           if (sri.size() < minSiz) {
                    // RDKit✔️✔️:             minSiz = rdcast<unsigned int>(sri.size());
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           nrings.push_back(sri);
                    min_size = min_size.min(ring.len());
                    new_rings.push(ring);
                }
            }
            // RDKit✔️✔️:       for (const auto &nring : nrings) {
            for ring in new_rings {
                // RDKit✔️✔️:         if (nring.size() == minSiz) {
                if ring.len() == min_size {
                    // RDKit✔️✔️:           auto invr = RingUtils::computeRingInvariant(nring, mol.getNumAtoms());
                    // RDKit✔️✔️:           if (invars.find(invr) == invars.end()) {
                    let invariant = compute_ring_invariant(&ring, context.atom_count());
                    if !invars.contains(&invariant) {
                        // RDKit✔️✔️:             res.push_back(nring);
                        // RDKit✔️✔️:             invars.insert(invr);
                        res.push(ring);
                        invars.insert(invariant);
                    }
                }
            }
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::findSSSRforDupCands
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn find_rings_d2_nodes(
    context: &RingSearchContext<'_>,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvariant>,
    d2nodes: &[usize],
    atom_degrees: &mut [i32],
    active_bonds: &mut [bool],
    ring_bonds: &mut [bool],
    ring_atoms: &mut [bool],
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findRingsD2nodes
    // RDKit✔️✔️: void findRingsD2nodes(const ROMol &tMol, VECT_INT_VECT &res,
    // RDKit✔️✔️:                       RINGINVAR_SET &invars, const INT_VECT &d2nodes,
    // RDKit✔️✔️:                       INT_VECT &atomDegrees,
    // RDKit✔️✔️:                       boost::dynamic_bitset<> &activeBonds,
    // RDKit✔️✔️:                       boost::dynamic_bitset<> &ringBonds,
    // RDKit✔️✔️:                       boost::dynamic_bitset<> &ringAtoms) {
    // RDKit✔️✔️:   RINGINVAR_INT_VECT_MAP dupD2Cands;
    // RDKit✔️✔️:   INT_INT_VECT_MAP dupMap;
    let mut dup_d2_candidates: BTreeMap<RingInvariant, Vec<usize>> = BTreeMap::new();
    let mut dup_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    // RDKit✔️✔️:   for (auto &cand : d2nodes) {
    for &candidate in d2nodes {
        // RDKit✔️✔️:     VECT_INT_VECT srings;
        // RDKit✔️✔️:     smallestRingsBfs(tMol, cand, srings, activeBonds);
        let smallest_rings = smallest_rings_bfs(context, candidate, active_bonds, None)?;
        // RDKit✔️✔️:     for (const auto &nring : srings) {
        for ring in &smallest_rings {
            // RDKit✔️✔️:       auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
            let invariant = compute_ring_invariant(ring, context.atom_count());
            // RDKit✔️✔️:       auto &duplicateInvars = dupD2Cands[invr];
            let duplicate_invariants = dup_d2_candidates.entry(invariant.clone()).or_default();
            // RDKit✔️✔️:       if (invars.find(invr) == invars.end()) {
            if !invars.contains(&invariant) {
                // RDKit✔️✔️:         res.push_back(nring);
                // RDKit✔️✔️:         invars.insert(invr);
                res.push(ring.clone());
                invars.insert(invariant);
                // RDKit✔️✔️:         for (unsigned int i = 0; i < nring.size() - 1; ++i) {
                for i in 0..ring.len() - 1 {
                    // RDKit✔️✔️:           unsigned int bIdx =
                    // RDKit✔️✔️:               tMol.getBondBetweenAtoms(nring[i], nring[i + 1])->getIdx();
                    let Some(bond) = context.bond_between_atoms(ring[i], ring[i + 1]) else {
                        return Err(RingFindingError::ExpectedBondNotFound {
                            begin: AtomId::new(ring[i]),
                            end: AtomId::new(ring[i + 1]),
                        });
                    };
                    // RDKit✔️✔️:           ringBonds.set(bIdx);
                    // RDKit✔️✔️:           ringAtoms.set(nring[i]);
                    ring_bonds[bond.index()] = true;
                    ring_atoms[ring[i]] = true;
                }
                // RDKit✔️✔️:         ringBonds.set(
                // RDKit✔️✔️:             tMol.getBondBetweenAtoms(nring[0], nring[nring.size() - 1])
                // RDKit✔️✔️:                 ->getIdx());
                // RDKit✔️✔️:         ringAtoms.set(nring[nring.size() - 1]);
                let Some(bond) = context.bond_between_atoms(ring[0], ring[ring.len() - 1]) else {
                    return Err(RingFindingError::ExpectedBondNotFound {
                        begin: AtomId::new(ring[0]),
                        end: AtomId::new(ring[ring.len() - 1]),
                    });
                };
                ring_bonds[bond.index()] = true;
                ring_atoms[ring[ring.len() - 1]] = true;
            } else {
                // RDKit✔️✔️:       } else {
                // RDKit✔️✔️:         for (auto otherCand : duplicateInvars) {
                for &other_candidate in duplicate_invariants.iter() {
                    // RDKit✔️✔️:           dupMap[cand].push_back(otherCand);
                    // RDKit✔️✔️:           dupMap[otherCand].push_back(cand);
                    dup_map.entry(candidate).or_default().push(other_candidate);
                    dup_map.entry(other_candidate).or_default().push(candidate);
                }
            }
            // RDKit✔️✔️:       duplicateInvars.push_back(cand);
            duplicate_invariants.push(candidate);
        }
        // RDKit✔️✔️:     if (srings.empty()) {
        if smallest_rings.is_empty() {
            // RDKit✔️✔️:       std::queue<int> changed;
            // RDKit✔️✔️:       changed.push(cand);
            let mut changed = VecDeque::new();
            changed.push_back(candidate);
            // RDKit✔️✔️:       while (!changed.empty()) {
            while let Some(local_candidate) = changed.pop_front() {
                // RDKit✔️✔️:         auto local_cand = changed.front();
                // RDKit✔️✔️:         changed.pop();
                // RDKit✔️✔️:         trimBonds(local_cand, tMol, changed, atomDegrees, activeBonds);
                trim_bonds(
                    context,
                    local_candidate,
                    &mut changed,
                    atom_degrees,
                    active_bonds,
                );
            }
        }
    }
    // RDKit✔️✔️:   findSSSRforDupCands(tMol, res, invars, dupMap, dupD2Cands, atomDegrees,
    // RDKit✔️✔️:                       activeBonds);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::findRingsD2nodes
    find_sssr_for_dup_candidates(
        context,
        res,
        invars,
        &dup_map,
        &dup_d2_candidates,
        atom_degrees,
        active_bonds,
    )
}

fn find_rings_d3_node(
    context: &RingSearchContext<'_>,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvariant>,
    candidate: usize,
    active_bonds: &[bool],
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findRingsD3Node
    // RDKit✔️✔️: void findRingsD3Node(const ROMol &tMol, VECT_INT_VECT &res,
    // RDKit✔️✔️:                      RINGINVAR_SET &invars, int cand, INT_VECT &,
    // RDKit✔️✔️:                      boost::dynamic_bitset<> &activeBonds) {
    // RDKit✔️✔️:   VECT_INT_VECT srings;
    // RDKit✔️✔️:   auto nsmall = smallestRingsBfs(tMol, cand, srings, activeBonds);
    let smallest = smallest_rings_bfs(context, candidate, active_bonds, None)?;
    let nsmall = smallest.len();
    // RDKit✔️✔️:   for (const auto &nring : srings) {
    for ring in &smallest {
        // RDKit✔️✔️:     auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
        // RDKit✔️✔️:     if (invars.find(invr) == invars.end()) {
        let invariant = compute_ring_invariant(ring, context.atom_count());
        if !invars.contains(&invariant) {
            // RDKit✔️✔️:       res.push_back(nring);
            // RDKit✔️✔️:       invars.insert(invr);
            res.push(ring.clone());
            invars.insert(invariant);
        }
    }
    // RDKit✔️✔️:   if (nsmall < 3) {
    if nsmall < 3 {
        // RDKit✔️✔️:     int n1 = -1, n2 = -1, n3 = -1;
        // RDKit✔️✔️:     for (auto bond : tMol.atomBonds(tMol.getAtomWithIdx(cand))) {
        // RDKit✔️✔️:       if (!activeBonds[bond->getIdx()]) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (n1 == -1) {
        // RDKit✔️✔️:         n1 = bond->getOtherAtomIdx(cand);
        // RDKit✔️✔️:       } else if (n2 == -1) {
        // RDKit✔️✔️:         n2 = bond->getOtherAtomIdx(cand);
        // RDKit✔️✔️:       } else if (n3 == -1) {
        // RDKit✔️✔️:         n3 = bond->getOtherAtomIdx(cand);
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let neighbors = context
            .neighbors(candidate)
            .iter()
            .filter(|neighbor| active_bonds[neighbor.bond.index()])
            .map(|neighbor| neighbor.atom_index)
            .take(3)
            .collect::<Vec<_>>();
        // RDKit✔️✔️:     CHECK_INVARIANT(n3 != -1, "neighbor not found");
        if neighbors.len() != 3 {
            return Err(RingFindingError::Value {
                message: "neighbor not found",
            });
        }
        let [n1, n2, n3] = [neighbors[0], neighbors[1], neighbors[2]];
        // RDKit✔️✔️:     if (nsmall == 2) {
        if nsmall == 2 {
            // RDKit✔️✔️:       int f = -1;
            // RDKit✔️✔️:       if ((std::find(srings[0].begin(), srings[0].end(), n1) !=
            // RDKit✔️✔️:            srings[0].end()) &&
            // RDKit✔️✔️:           (std::find(srings[1].begin(), srings[1].end(), n1) !=
            // RDKit✔️✔️:            srings[1].end())) {
            // RDKit✔️✔️:         f = n1;
            // RDKit✔️✔️:       } else if ((std::find(srings[0].begin(), srings[0].end(), n2) !=
            // RDKit✔️✔️:                   srings[0].end()) &&
            // RDKit✔️✔️:                  (std::find(srings[1].begin(), srings[1].end(), n2) !=
            // RDKit✔️✔️:                   srings[1].end())) {
            // RDKit✔️✔️:         f = n2;
            // RDKit✔️✔️:       } else if ((std::find(srings[0].begin(), srings[0].end(), n3) !=
            // RDKit✔️✔️:                   srings[0].end()) &&
            // RDKit✔️✔️:                  (std::find(srings[1].begin(), srings[1].end(), n3) !=
            // RDKit✔️✔️:                   srings[1].end())) {
            // RDKit✔️✔️:         f = n3;
            // RDKit✔️✔️:       }
            let f = [n1, n2, n3]
                .into_iter()
                .find(|neighbor| smallest.iter().all(|ring| ring.contains(neighbor)));
            // RDKit✔️✔️:       CHECK_INVARIANT(f >= 0, "third ring not found");
            let Some(forbidden_atom) = f else {
                return Err(RingFindingError::Value {
                    message: "third ring not found",
                });
            };
            // RDKit✔️✔️:       VECT_INT_VECT trings;
            // RDKit✔️✔️:       INT_VECT forb;
            // RDKit✔️✔️:       forb.push_back(f);
            // RDKit✔️✔️:       smallestRingsBfs(tMol, cand, trings, activeBonds, &forb);
            let rings =
                smallest_rings_bfs(context, candidate, active_bonds, Some(&[forbidden_atom]))?;
            // RDKit✔️✔️:       for (const auto &nring : trings) {
            for ring in rings {
                // RDKit✔️✔️:         auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
                // RDKit✔️✔️:         if (invars.find(invr) == invars.end()) {
                let invariant = compute_ring_invariant(&ring, context.atom_count());
                if !invars.contains(&invariant) {
                    // RDKit✔️✔️:           res.push_back(nring);
                    // RDKit✔️✔️:           invars.insert(invr);
                    res.push(ring);
                    invars.insert(invariant);
                }
            }
        }
        // RDKit✔️✔️:     if (nsmall == 1) {
        if nsmall == 1 {
            // RDKit✔️✔️:       int f1 = -1, f2 = -1;
            let (f1, f2) = if !smallest[0].contains(&n1) {
                (n2, n3)
            } else if !smallest[0].contains(&n2) {
                (n1, n3)
            } else if !smallest[0].contains(&n3) {
                (n1, n2)
            } else {
                return Err(RingFindingError::Value {
                    message: "rings not found",
                });
            };
            // RDKit✔️✔️:       CHECK_INVARIANT(f1 >= 0, "rings not found");
            // RDKit✔️✔️:       CHECK_INVARIANT(f2 >= 0, "rings not found");
            // RDKit✔️✔️:       VECT_INT_VECT trings;
            // RDKit✔️✔️:       INT_VECT forb;
            // RDKit✔️✔️:       forb.push_back(f2);
            // RDKit✔️✔️:       smallestRingsBfs(tMol, cand, trings, activeBonds, &forb);
            let rings = smallest_rings_bfs(context, candidate, active_bonds, Some(&[f2]))?;
            for ring in rings {
                let invariant = compute_ring_invariant(&ring, context.atom_count());
                if !invars.contains(&invariant) {
                    res.push(ring);
                    invars.insert(invariant);
                }
            }
            // RDKit✔️✔️:       trings.clear();
            // RDKit✔️✔️:       forb.clear();
            // RDKit✔️✔️:       forb.push_back(f1);
            // RDKit✔️✔️:       smallestRingsBfs(tMol, cand, trings, activeBonds, &forb);
            let rings = smallest_rings_bfs(context, candidate, active_bonds, Some(&[f1]))?;
            for ring in rings {
                let invariant = compute_ring_invariant(&ring, context.atom_count());
                if !invars.contains(&invariant) {
                    res.push(ring);
                    invars.insert(invariant);
                }
            }
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::findRingsD3Node
    Ok(())
}

fn remove_extra_rings(
    context: &RingSearchContext<'_>,
    res: &mut Vec<Vec<usize>>,
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::removeExtraRings complete source
    // RDKit✔️❌: void removeExtraRings(VECT_INT_VECT &res, unsigned int, const ROMol &mol) {
    // RDKit✔️❌:   // sort on size
    // RDKit✔️❌:   std::sort(res.begin(), res.end(), compRingSize);
    // RDKit✔️❌:
    // RDKit✔️❌:   // change the rings from atom IDs to bondIds
    // RDKit✔️❌:   VECT_INT_VECT brings;
    // RDKit✔️❌:   RingUtils::convertToBonds(res, brings, mol);
    // RDKit✔️❌:   std::vector<boost::dynamic_bitset<>> bitBrings;
    // RDKit✔️❌:   bitBrings.reserve(brings.size());
    // RDKit✔️❌:   for (const auto &vivi : brings) {
    // RDKit✔️❌:     boost::dynamic_bitset<> lring(mol.getNumBonds());
    // RDKit✔️❌:     for (int ivi : vivi) {
    // RDKit✔️❌:       lring.set(ivi);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     bitBrings.push_back(lring);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> availRings(res.size());
    // RDKit✔️❌:   availRings.set();
    // RDKit✔️❌:   boost::dynamic_bitset<> keepRings(res.size());
    // RDKit✔️❌:   boost::dynamic_bitset<> munion(mol.getNumBonds());
    // RDKit✔️❌:
    // RDKit✔️❌:   // optimization - don't reallocate a new one each loop
    // RDKit✔️❌:   boost::dynamic_bitset<> workspace(mol.getNumBonds());
    // RDKit✔️❌:
    // RDKit✔️❌:   for (unsigned int i = 0; i < res.size(); ++i) {
    // RDKit✔️❌:     // skip this ring if we've already seen all of its bonds
    // RDKit✔️❌:     if (bitBrings[i].is_subset_of(munion)) {
    // RDKit✔️❌:       availRings.set(i, 0);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (!availRings[i]) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     munion |= bitBrings[i];
    // RDKit✔️❌:     keepRings.set(i);
    // RDKit✔️❌:
    // RDKit✔️❌:     // from this ring we consider all others that are still available and the
    // RDKit✔️❌:     // same size
    // RDKit✔️❌:     boost::dynamic_bitset<> consider(res.size());
    // RDKit✔️❌:     for (unsigned int j = i + 1; j < res.size(); ++j) {
    // RDKit✔️❌:       if (availRings[j] && (brings[j].size() == brings[i].size())) {
    // RDKit✔️❌:         consider.set(j);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     while (consider.any()) {
    // RDKit✔️❌:       unsigned int bestJ = i + 1;
    // RDKit✔️❌:       int bestOverlap = -1;
    // RDKit✔️❌:       // loop over the available other rings in consideration and pick the one
    // RDKit✔️❌:       // that has the most overlapping bonds with what we've done so far.
    // RDKit✔️❌:       // this is the fix to github #526
    // RDKit✔️❌:       for (unsigned int j = i + 1;
    // RDKit✔️❌:            j < res.size() && brings[j].size() == brings[i].size(); ++j) {
    // RDKit✔️❌:         if (!consider[j] || !availRings[j]) {
    // RDKit✔️❌:           continue;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         workspace = bitBrings[j];
    // RDKit✔️❌:         workspace &= munion;
    // RDKit✔️❌:         int overlap = rdcast<int>(workspace.count());
    // RDKit✔️❌:         if (overlap > bestOverlap) {
    // RDKit✔️❌:           bestOverlap = overlap;
    // RDKit✔️❌:           bestJ = j;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       consider.set(bestJ, 0);
    // RDKit✔️❌:       if (bitBrings[bestJ].is_subset_of(munion)) {
    // RDKit✔️❌:         availRings.set(bestJ, 0);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         keepRings.set(bestJ);
    // RDKit✔️❌:         availRings.set(bestJ, 0);
    // RDKit✔️❌:         munion |= bitBrings[bestJ];
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // remove the extra rings from res and store them on the molecule in case we
    // RDKit✔️❌:   // wish symmetrize the SSSRs later
    // RDKit✔️❌:   VECT_INT_VECT extras;
    // RDKit✔️❌:   VECT_INT_VECT temp = res;
    // RDKit✔️❌:   res.resize(0);
    // RDKit✔️❌:   for (unsigned int i = 0; i < temp.size(); i++) {
    // RDKit✔️❌:     if (keepRings[i]) {
    // RDKit✔️❌:       res.push_back(temp[i]);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       extras.push_back(temp[i]);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // add extra rings to the molecule (there could already be some from previous
    // RDKit✔️❌:   // fragments)
    // RDKit✔️❌:   VECT_INT_VECT molExtras;
    // RDKit✔️❌:   mol.getPropIfPresent(common_properties::extraRings, molExtras);
    // RDKit✔️❌:   molExtras.insert(molExtras.end(), extras.begin(), extras.end());
    // RDKit✔️❌:   mol.setProp(common_properties::extraRings, molExtras, true);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION FindRings::removeExtraRings complete source
    // Return the ordered extras once; the canonical caller extends its cache,
    // matching native getPropIfPresent + insert of prior-fragment extras.
    // The caller registers the computed name through the canonical MODEL
    // setProp prefix after storing this exact owner-local typed matrix.
    // Local performance review: BTreeSet bond membership is materially slower
    // than native packed dynamic_bitset subset/intersection/union operations;
    // the cost marker remains negative despite equivalent set/order semantics.

    rdkit_std_sort_rings_by_size(res);
    let bond_rings = convert_rings_to_bonds(context, res)?;
    let bit_bond_rings = bond_rings
        .iter()
        .map(|ring| ring.iter().copied().collect::<BTreeSet<_>>())
        .collect::<Vec<_>>();
    let mut available = vec![true; res.len()];
    let mut keep = vec![false; res.len()];
    let mut union = BTreeSet::new();
    for i in 0..res.len() {
        if bit_bond_rings[i].is_subset(&union) {
            available[i] = false;
        }
        if !available[i] {
            continue;
        }
        union.extend(bit_bond_rings[i].iter().copied());
        keep[i] = true;
        let mut consider = vec![false; res.len()];
        for j in i + 1..res.len() {
            if available[j] && bond_rings[j].len() == bond_rings[i].len() {
                consider[j] = true;
            }
        }
        while consider.iter().any(|flag| *flag) {
            let mut best_j = i + 1;
            let mut best_overlap = -1isize;
            for j in i + 1..res.len() {
                if bond_rings[j].len() != bond_rings[i].len() {
                    break;
                }
                if !consider[j] || !available[j] {
                    continue;
                }
                let overlap = bit_bond_rings[j].intersection(&union).count() as isize;
                if overlap > best_overlap {
                    best_overlap = overlap;
                    best_j = j;
                }
            }
            consider[best_j] = false;
            if bit_bond_rings[best_j].is_subset(&union) {
                available[best_j] = false;
            } else {
                keep[best_j] = true;
                available[best_j] = false;
                union.extend(bit_bond_rings[best_j].iter().copied());
            }
        }
    }
    let old = std::mem::take(res);
    let mut extras = Vec::new();
    for (index, ring) in old.into_iter().enumerate() {
        if keep[index] {
            res.push(ring);
        } else {
            extras.push(ring);
        }
    }
    Ok(extras)
}

fn rdkit_std_sort_rings_by_size(rings: &mut [Vec<usize>]) {
    // BEGIN RDKIT CPP FUNCTION FindRings::removeExtraRings
    // RDKit❗✔️: auto compRingSize = [](const auto &v1, const auto &v2) {
    // RDKit❗✔️:   return v1.size() < v2.size();
    // RDKit❗✔️: };
    // RDKit❗✔️: std::sort(res.begin(), res.end(), compRingSize);
    // END RDKIT CPP FUNCTION FindRings::removeExtraRings
    crate::source_sort::sort_by(rings, |left, right| left.len() < right.len());
}

fn atom_search_bfs(
    context: &RingSearchContext<'_>,
    start_atom: usize,
    end_atom: usize,
    ring_atoms: &[bool],
    invars: &BTreeSet<RingInvariant>,
) -> Result<Option<Vec<usize>>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::_atomSearchBFS
    // RDKit✔️✔️: bool _atomSearchBFS(const ROMol &tMol, unsigned int startAtomIdx,
    // RDKit✔️✔️:                     unsigned int endAtomIdx, boost::dynamic_bitset<> &ringAtoms,
    // RDKit✔️✔️:                     INT_VECT &res, RINGINVAR_SET &invars) {
    // RDKit✔️✔️:   res.clear();
    // RDKit✔️✔️:   std::deque<INT_VECT> bfsq;
    // RDKit✔️✔️:   INT_VECT tv;
    // RDKit✔️✔️:   tv.push_back(startAtomIdx);
    // RDKit✔️✔️:   bfsq.push_back(tv);
    let mut bfsq = VecDeque::new();
    bfsq.push_back(vec![start_atom]);
    // RDKit✔️✔️:   while (!bfsq.empty()) {
    while let Some(path) = bfsq.pop_front() {
        // RDKit✔️✔️:     if (bfsq.size() >= RingUtils::MAX_BFSQ_SIZE) {
        // RDKit✔️✔️:       constexpr const char *msg =
        // RDKit✔️✔️:           "Maximum BFS search size exceeded.\nThis is likely due to a highly "
        // RDKit✔️✔️:           "symmetric fused ring system.";
        // RDKit✔️✔️:       BOOST_LOG(rdErrorLog) << msg << std::endl;
        // RDKit✔️✔️:       throw ValueErrorException(msg);
        // RDKit✔️✔️:     }
        if bfsq.len() >= MAX_BFSQ_SIZE {
            return Err(RingFindingError::Value {
                message: "Maximum BFS search size exceeded.\nThis is likely due to a highly symmetric fused ring system.",
            });
        }
        // RDKit✔️✔️:     tv = bfsq.front();
        // RDKit✔️✔️:     bfsq.pop_front();
        // RDKit✔️✔️:     unsigned int currAtomIdx = tv.back();
        let current = *path.last().expect("BFS path is nonempty");
        // RDKit✔️✔️:     for (auto nbr : tMol.atomNeighbors(tMol.getAtomWithIdx(currAtomIdx))) {
        for neighbor in context.neighbors(current) {
            // RDKit✔️✔️:       auto nbrIdx = nbr->getIdx();
            let neighbor_idx = neighbor.atom_index;
            // RDKit✔️✔️:       if (nbrIdx == endAtomIdx) {
            if neighbor_idx == end_atom {
                // RDKit✔️✔️:         if (currAtomIdx != startAtomIdx) {
                if current != start_atom {
                    // RDKit✔️✔️:           INT_VECT nv(tv);
                    // RDKit✔️✔️:           nv.push_back(rdcast<unsigned int>(nbrIdx));
                    let mut new_path = path.clone();
                    new_path.push(neighbor_idx);
                    // RDKit✔️✔️:           auto invr = RingUtils::computeRingInvariant(nv, tMol.getNumAtoms());
                    let invariant = compute_ring_invariant(&new_path, context.atom_count());
                    // RDKit✔️✔️:           if (invars.find(invr) == invars.end()) {
                    // RDKit✔️✔️:             res.resize(nv.size());
                    // RDKit✔️✔️:             std::copy(nv.begin(), nv.end(), res.begin());
                    // RDKit✔️✔️:             return true;
                    // RDKit✔️✔️:           }
                    if !invars.contains(&invariant) {
                        return Ok(Some(new_path));
                    }
                }
                // RDKit✔️✔️:         }
                // RDKit✔️✔️:       } else if (ringAtoms[nbrIdx] &&
                // RDKit✔️✔️:                  std::find(tv.begin(), tv.end(), nbrIdx) == tv.end()) {
            } else if ring_atoms[neighbor_idx] && !path.contains(&neighbor_idx) {
                // RDKit✔️✔️:         INT_VECT nv(tv);
                // RDKit✔️✔️:         nv.push_back(rdcast<unsigned int>(nbrIdx));
                // RDKit✔️✔️:         bfsq.push_back(nv);
                let mut new_path = path.clone();
                new_path.push(neighbor_idx);
                bfsq.push_back(new_path);
            }
            // RDKit✔️✔️:       }
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::_atomSearchBFS
    Ok(None)
}

fn find_ring_connecting_atoms(
    context: &RingSearchContext<'_>,
    bond: BondId,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<RingInvariant>,
    ring_bonds: &mut [bool],
    ring_atoms: &mut [bool],
) -> Result<bool, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findRingConnectingAtoms
    // RDKit✔️✔️: bool findRingConnectingAtoms(const ROMol &tMol, const Bond *bond,
    // RDKit✔️✔️:                              VECT_INT_VECT &res, RINGINVAR_SET &invars,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> &ringBonds,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> &ringAtoms) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(!ringBonds[bond->getIdx()], "not a ring bond");
    // RDKit✔️✔️:   PRECONDITION(ringAtoms[bond->getBeginAtomIdx()], "not a ring atom");
    // RDKit✔️✔️:   PRECONDITION(ringAtoms[bond->getEndAtomIdx()], "not a ring atom");
    let bond_ref = &context.bonds()[bond.index()];
    // RDKit✔️✔️:   INT_VECT nring;
    // RDKit✔️✔️:   if (_atomSearchBFS(tMol, bond->getBeginAtomIdx(), bond->getEndAtomIdx(),
    // RDKit✔️✔️:                      ringAtoms, nring, invars)) {
    let Some(ring) = atom_search_bfs(
        context,
        bond_ref.begin().index(),
        bond_ref.end().index(),
        ring_atoms,
        invars,
    )?
    else {
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        return Ok(false);
    };
    // RDKit✔️✔️:     auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
    let invariant = compute_ring_invariant(&ring, context.atom_count());
    // RDKit✔️✔️:     if (invars.find(invr) == invars.end()) {
    if !invars.contains(&invariant) {
        // RDKit✔️✔️:       res.push_back(nring);
        // RDKit✔️✔️:       invars.insert(invr);
        // RDKit✔️✔️:       for (unsigned int i = 0; i < nring.size() - 1; ++i) {
        for i in 0..ring.len() - 1 {
            // RDKit✔️✔️:         unsigned int bIdx =
            // RDKit✔️✔️:             tMol.getBondBetweenAtoms(nring[i], nring[i + 1])->getIdx();
            let Some(bond) = context.bond_between_atoms(ring[i], ring[i + 1]) else {
                return Err(RingFindingError::ExpectedBondNotFound {
                    begin: AtomId::new(ring[i]),
                    end: AtomId::new(ring[i + 1]),
                });
            };
            // RDKit✔️✔️:         ringBonds.set(bIdx);
            // RDKit✔️✔️:         ringAtoms.set(nring[i]);
            ring_bonds[bond.index()] = true;
            ring_atoms[ring[i]] = true;
        }
        // RDKit✔️✔️:       ringBonds.set(tMol.getBondBetweenAtoms(nring[0], nring[nring.size() - 1])
        // RDKit✔️✔️:                         ->getIdx());
        let Some(bond) = context.bond_between_atoms(ring[0], ring[ring.len() - 1]) else {
            return Err(RingFindingError::ExpectedBondNotFound {
                begin: AtomId::new(ring[0]),
                end: AtomId::new(ring[ring.len() - 1]),
            });
        };
        // RDKit✔️✔️:       ringAtoms.set(nring[nring.size() - 1]);
        ring_bonds[bond.index()] = true;
        ring_atoms[ring[ring.len() - 1]] = true;
        res.push(ring);
        invars.insert(invariant);
    }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FindRings::findRingConnectingAtoms
    Ok(true)
}

fn molecule_fragments(context: &RingSearchContext<'_>) -> Vec<Vec<usize>> {
    crate::paths::connected_components_from_source(context)
        .components
        .into_iter()
        .map(|component| component.into_iter().map(AtomId::index).collect())
        .collect()
}

fn fast_find_rings_internal(
    context: &RingSearchContext<'_>,
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    // Discovery phase of the complete source wrapper above.
    // RDKit✔️✔️: void fastFindRings(const ROMol &mol) {
    // RDKit✔️✔️:   VECT_INT_VECT res;
    // RDKit✔️✔️:   res.resize(0);
    let mut result = Vec::new();
    // RDKit✔️✔️:   unsigned int nats = mol.getNumAtoms();
    // RDKit✔️✔️:   INT_VECT atomColors(nats, 0);
    let mut atom_colors = vec![0u8; context.atom_count()];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nats; ++i) {
    for atom in 0..context.atom_count() {
        // RDKit✔️✔️:     if (atomColors[i]) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if atom_colors[atom] != 0 {
            continue;
        }
        // RDKit✔️✔️:     if (mol.getAtomWithIdx(i)->getDegree() < 2) {
        // RDKit✔️✔️:       atomColors[i] = 2;
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if context.neighbors(atom).len() < 2 {
            atom_colors[atom] = 2;
            continue;
        }
        // RDKit✔️✔️:     std::vector<const Atom *> traversalOrder;
        let mut traversal_order = Vec::new();
        // RDKit✔️✔️:     _DFS(mol, mol.getAtomWithIdx(i), atomColors, traversalOrder, res);
        dfs_fast_find_rings(
            context,
            atom,
            &mut atom_colors,
            &mut traversal_order,
            &mut result,
            None,
        );
    }
    // RDKit✔️✔️: }
    // End source discovery phase.
    Ok(result)
}

fn dfs_fast_find_rings(
    context: &RingSearchContext<'_>,
    atom: usize,
    atom_colors: &mut [u8],
    traversal_order: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
    from_atom: Option<usize>,
) {
    // BEGIN RDKIT CPP FUNCTION MolOps::_DFS complete source
    // RDKit✔️✔️: void _DFS(const ROMol &mol, const Atom *atom, INT_VECT &atomColors,
    // RDKit✔️✔️:           std::vector<const Atom *> &traversalOrder, VECT_INT_VECT &res,
    // RDKit✔️✔️:           const Atom *fromAtom = nullptr) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(atomColors[atom->getIdx()] == 0, "bad color");
    // RDKit✔️✔️:   atomColors[atom->getIdx()] = 1;
    // RDKit✔️✔️:   traversalOrder.push_back(atom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto nbr : mol.atomNeighbors(atom)) {
    // RDKit✔️✔️:     unsigned int nbrIdx = nbr->getIdx();
    // RDKit✔️✔️:     if (atomColors[nbrIdx] == 0) {
    // RDKit✔️✔️:       if (nbr->getDegree() < 2) {
    // RDKit✔️✔️:         atomColors[nbr->getIdx()] = 2;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         _DFS(mol, nbr, atomColors, traversalOrder, res, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (atomColors[nbrIdx] == 1) {
    // RDKit✔️✔️:       if (fromAtom && nbrIdx != fromAtom->getIdx()) {
    // RDKit✔️✔️:         INT_VECT cycle;
    // RDKit✔️✔️:         auto lastElem =
    // RDKit✔️✔️:             std::find(traversalOrder.rbegin(), traversalOrder.rend(), atom);
    // RDKit✔️✔️:         for (auto rIt = lastElem;  // traversalOrder.rbegin();
    // RDKit✔️✔️:              rIt != traversalOrder.rend() && (*rIt)->getIdx() != nbrIdx;
    // RDKit✔️✔️:              ++rIt) {
    // RDKit✔️✔️:           cycle.push_back((*rIt)->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         cycle.push_back(nbrIdx);
    // RDKit✔️✔️:         res.push_back(cycle);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atomColors[atom->getIdx()] = 2;
    // RDKit✔️✔️:   traversalOrder.pop_back();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::_DFS complete source
    // Current atom is the back of traversalOrder after its own push and all
    // recursive pushes/pops; native reverse find(atom) therefore starts here.
    // Each cycle follows the same reverse traversal until the gray neighbor;
    // source root-parent suppression and colors 0/1/2 preserve ring order.

    atom_colors[atom] = 1;
    traversal_order.push(atom);
    for neighbor in context.neighbors(atom) {
        let neighbor_idx = neighbor.atom_index;
        if atom_colors[neighbor_idx] == 0 {
            if context.neighbors(neighbor_idx).len() < 2 {
                atom_colors[neighbor_idx] = 2;
            } else {
                dfs_fast_find_rings(
                    context,
                    neighbor_idx,
                    atom_colors,
                    traversal_order,
                    result,
                    Some(atom),
                );
            }
        } else if atom_colors[neighbor_idx] == 1
            && from_atom.is_some()
            && Some(neighbor_idx) != from_atom
        {
            let mut cycle = Vec::new();
            for &path_atom in traversal_order.iter().rev() {
                if path_atom == neighbor_idx {
                    break;
                }
                cycle.push(path_atom);
            }
            cycle.push(neighbor_idx);
            result.push(cycle);
        }
    }
    atom_colors[atom] = 2;
    traversal_order.pop();
}

#[cfg(test)]
mod selected_row_tests {
    use super::*;

    use cosmolkit_model::{Atom, AtomSpec, Bond, BondSpec};
    use cosmolkit_types::{BondOrder, Element};

    #[test]
    fn drawing_ring_preallocate_eight_carrier_preservation_calls() {
        let mut calls = 0;
        for quality in [RingFindType::Sssr, RingFindType::SymmSssr] {
            for has_cycle in [false, true] {
                for target in [3, 5] {
                    let mut carrier = RingInfo::new(quality, 3, 3);
                    if has_cycle {
                        carrier.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
                        carrier.atom_ring_families =
                            vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]];
                        carrier.bond_ring_families =
                            vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]];
                        carrier.relevant_cycle_count = Some(1);
                        carrier.fused_rings = vec![vec![false]];
                        carrier.num_fused_bonds = vec![0];
                    }
                    assert!(carrier.is_initialized());
                    assert_eq!(carrier.find_type(), quality);
                    assert_eq!(carrier.atom_row_count(), 3);
                    assert_eq!(carrier.bond_row_count(), 3);
                    let original_members = if has_cycle {
                        vec![vec![0], vec![0], vec![0]]
                    } else {
                        vec![vec![], vec![], vec![]]
                    };
                    assert_eq!(carrier.atom_members, original_members);
                    assert_eq!(carrier.bond_members, original_members);
                    if has_cycle {
                        assert_eq!(
                            carrier.atom_rings,
                            vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]]
                        );
                        assert_eq!(
                            carrier.bond_rings,
                            vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]]
                        );
                    } else {
                        assert!(carrier.atom_rings.is_empty());
                        assert!(carrier.bond_rings.is_empty());
                    }
                    let before = carrier.clone();
                    let atom_rows_ptr = carrier.atom_rings.as_ptr();
                    let bond_rows_ptr = carrier.bond_rings.as_ptr();
                    let cycle_ptrs = has_cycle.then(|| {
                        (
                            carrier.atom_rings[0].as_ptr(),
                            carrier.bond_rings[0].as_ptr(),
                        )
                    });
                    let literal_members = match (has_cycle, target) {
                        (false, 3) => vec![vec![], vec![], vec![]],
                        (false, 5) => vec![vec![], vec![], vec![], vec![], vec![]],
                        (true, 3) => vec![vec![0], vec![0], vec![0]],
                        (true, 5) => vec![vec![0], vec![0], vec![0], vec![], vec![]],
                        _ => unreachable!(),
                    };
                    let mut expected = before.clone();
                    expected.atom_members = literal_members.clone();
                    expected.bond_members = literal_members.clone();

                    carrier.preallocate(target, target);
                    calls += 1;

                    assert_eq!(carrier, expected);
                    assert_eq!(carrier.atom_row_count(), target);
                    assert_eq!(carrier.bond_row_count(), target);
                    for (index, members) in literal_members.iter().enumerate() {
                        assert_eq!(carrier.atom_members(AtomId::new(index)), members);
                        assert_eq!(carrier.bond_members(BondId::new(index)), members);
                    }
                    assert_eq!(carrier.initialized, before.initialized);
                    assert_eq!(carrier.find_type, before.find_type);
                    assert_eq!(carrier.atom_rings, before.atom_rings);
                    assert_eq!(carrier.bond_rings, before.bond_rings);
                    assert_eq!(carrier.atom_ring_families, before.atom_ring_families);
                    assert_eq!(carrier.bond_ring_families, before.bond_ring_families);
                    assert_eq!(carrier.relevant_cycle_count, before.relevant_cycle_count);
                    assert_eq!(carrier.fused_rings, before.fused_rings);
                    assert_eq!(carrier.num_fused_bonds, before.num_fused_bonds);
                    if let Some((atoms_ptr, bonds_ptr)) = cycle_ptrs {
                        assert_eq!(carrier.atom_rings.as_ptr(), atom_rows_ptr);
                        assert_eq!(carrier.bond_rings.as_ptr(), bond_rows_ptr);
                        assert_eq!(carrier.atom_rings[0].as_ptr(), atoms_ptr);
                        assert_eq!(carrier.bond_rings[0].as_ptr(), bonds_ptr);
                    }
                }
            }
        }
        assert_eq!(calls, 8, "exact post-invocation census");
    }

    fn six_cycle_topology() -> TopologyBlock {
        let atoms = (0..6)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        let bonds = (0..6)
            .map(|id| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(
                        AtomId::new(id),
                        AtomId::new((id + 1) % 6),
                        BondOrder::Single,
                    ),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }

    // The exact K-CONVERT 8-call converter product: {F,R,[0,1,3],[0,1,2]}
    // atom rows EACH through the usize projection route (the internal owner
    // with identity atom projection and index output) and the typed
    // projection route (ring_atom_ids_to_bond_ids). Valid F derives
    // [0,1,2,3,4,5]; R derives [4,3,2,1,0,5]; the broken rows derive the
    // exact typed ExpectedBondNotFound causes. Full input stability and
    // output shape are asserted in BOTH representations.
    #[test]
    fn ring_conversion_typed_two_representations() {
        let graph = six_cycle_topology();
        let context =
            RingSearchContext::from_parts(graph.atoms.len(), &graph.bonds, &graph.adjacency);
        struct Case {
            name: &'static str,
            row: [usize; 6],
            expected: Result<[usize; 6], (usize, usize)>,
        }
        let cases = [
            Case {
                name: "F",
                row: [0, 1, 2, 3, 4, 5],
                expected: Ok([0, 1, 2, 3, 4, 5]),
            },
            Case {
                name: "R",
                row: [5, 4, 3, 2, 1, 0],
                expected: Ok([4, 3, 2, 1, 0, 5]),
            },
            Case {
                name: "gap-1-3",
                row: [0, 1, 3, 0, 1, 3],
                expected: Err((1, 3)),
            },
            Case {
                name: "gap-2-0",
                row: [0, 1, 2, 0, 1, 2],
                expected: Err((2, 0)),
            },
        ];
        let mut calls = 0usize;
        for case in &cases {
            // usize representation.
            calls += 1;
            let row_snapshot = case.row;
            let usize_result = convert_to_bonds(&context, &case.row, |atom| *atom, BondId::index);
            // typed representation.
            calls += 1;
            let typed_row: Vec<AtomId> = case.row.iter().map(|i| AtomId::new(*i)).collect();
            let typed_snapshot = typed_row.clone();
            let typed_result = ring_atom_ids_to_bond_ids(&graph, &typed_row);
            assert_eq!(case.row, row_snapshot, "{}: usize input mutated", case.name);
            assert_eq!(
                typed_row, typed_snapshot,
                "{}: typed input mutated",
                case.name
            );
            match case.expected {
                Ok(expected) => {
                    let usize_bonds = usize_result.unwrap();
                    let typed_bonds = typed_result.unwrap();
                    assert_eq!(usize_bonds.len(), case.row.len(), "{}: shape", case.name);
                    assert_eq!(typed_bonds.len(), case.row.len(), "{}: shape", case.name);
                    assert_eq!(usize_bonds, expected.to_vec(), "{}: usize", case.name);
                    assert_eq!(
                        typed_bonds
                            .iter()
                            .map(|bond| bond.index())
                            .collect::<Vec<_>>(),
                        expected.to_vec(),
                        "{}: typed",
                        case.name
                    );
                }
                Err((begin, end)) => {
                    assert_eq!(
                        usize_result.err(),
                        Some(RingFindingError::ExpectedBondNotFound {
                            begin: AtomId::new(begin),
                            end: AtomId::new(end),
                        }),
                        "{}: usize error",
                        case.name
                    );
                    assert_eq!(
                        typed_result.err(),
                        Some(RingFindingError::ExpectedBondNotFound {
                            begin: AtomId::new(begin),
                            end: AtomId::new(end),
                        }),
                        "{}: typed error",
                        case.name
                    );
                }
            }
        }
        assert_eq!(calls, 8, "exact census");
    }

    // The EXACT frozen converter product: three-element literal rows with a
    // REAL missing edge, unpadded. &[0,1,3] fails on the ADJACENT 1->3 gap;
    // &[0,1,2] fails on the CLOSING 2->0 edge. Each row runs through the
    // usize projection and the typed entry = 8 calls; fresh graph/row
    // preservation is compared immediately after EACH route (including
    // Err), BEFORE the second route.
    #[test]
    fn ring_conversion_typed_exact_closing() {
        struct Case {
            name: &'static str,
            row: &'static [usize],
            expected: Result<&'static [usize], (usize, usize)>,
        }
        let cases = [
            Case {
                name: "F",
                row: &[0, 1, 2, 3, 4, 5],
                expected: Ok(&[0, 1, 2, 3, 4, 5]),
            },
            Case {
                name: "R",
                row: &[5, 4, 3, 2, 1, 0],
                expected: Ok(&[4, 3, 2, 1, 0, 5]),
            },
            Case {
                name: "adjacent-gap-1-3",
                row: &[0, 1, 3],
                expected: Err((1, 3)),
            },
            Case {
                name: "closing-gap-2-0",
                row: &[0, 1, 2],
                expected: Err((2, 0)),
            },
        ];
        let mut calls = 0usize;
        for case in &cases {
            // Route 1: usize projection.
            // Fresh per-call graph and row baselines for THIS invocation.
            let graph = six_cycle_topology();
            assert_eq!(graph.bonds[5].begin().index(), 5);
            assert_eq!(graph.bonds[5].end().index(), 0, "closing b5=(5,0) exists");
            let graph_snapshot = graph.clone();
            let context =
                RingSearchContext::from_parts(graph.atoms.len(), &graph.bonds, &graph.adjacency);
            let row_snapshot: Vec<usize> = case.row.to_vec();
            if case.expected.is_err() {
                assert_eq!(case.row.len(), 3, "{}: gap row prerequisite", case.name);
            }
            let usize_result = convert_to_bonds(&context, case.row, |atom| *atom, BondId::index);
            calls += 1;
            assert_eq!(
                case.row,
                row_snapshot.as_slice(),
                "{}: usize row mutated",
                case.name
            );
            assert_eq!(
                &graph, &graph_snapshot,
                "{}: graph mutated (usize)",
                case.name
            );
            // Route 2: typed entry.
            // Fresh per-route graph and typed-row baselines BEFORE this
            // invocation; compared immediately after THIS route (including
            // Err), BEFORE any further work.
            let graph = six_cycle_topology();
            let graph_snapshot = graph.clone();
            let typed_row: Vec<AtomId> = case.row.iter().map(|i| AtomId::new(*i)).collect();
            let typed_snapshot = typed_row.clone();
            let typed_result = ring_atom_ids_to_bond_ids(&graph, &typed_row);
            calls += 1;
            assert_eq!(
                typed_row, typed_snapshot,
                "{}: typed row mutated",
                case.name
            );
            assert_eq!(
                &graph, &graph_snapshot,
                "{}: graph mutated (typed)",
                case.name
            );
            match case.expected {
                Ok(expected) => {
                    assert_eq!(
                        usize_result.unwrap(),
                        expected.to_vec(),
                        "{}: usize",
                        case.name
                    );
                    assert_eq!(
                        typed_result
                            .unwrap()
                            .iter()
                            .map(|b| b.index())
                            .collect::<Vec<_>>(),
                        expected.to_vec(),
                        "{}: typed",
                        case.name
                    );
                }
                Err((begin, end)) => {
                    assert_eq!(
                        usize_result.err(),
                        Some(RingFindingError::ExpectedBondNotFound {
                            begin: AtomId::new(begin),
                            end: AtomId::new(end),
                        }),
                        "{}: usize error",
                        case.name
                    );
                    assert_eq!(
                        typed_result.err(),
                        Some(RingFindingError::ExpectedBondNotFound {
                            begin: AtomId::new(begin),
                            end: AtomId::new(end),
                        }),
                        "{}: typed error",
                        case.name
                    );
                }
            }
        }
        assert_eq!(calls, 8, "exact census");
    }

    #[test]
    fn ring_info_from_selected_rows_keeps_initialized_empty_membership() {
        let rings = ring_info_from_selected_rows(5, 4, &[], &[]).unwrap();
        assert!(rings.is_initialized());
        assert_eq!(rings.find_type(), RingFindType::OtherOrUnknown);
        assert_eq!(rings.num_rings(), 0);
        assert_eq!(rings.atom_members(AtomId::new(4)), &[] as &[usize]);
        assert_eq!(rings.bond_members(BondId::new(3)), &[] as &[usize]);
        assert!(!rings.is_find_fast_or_better());
    }

    #[test]
    fn ring_info_from_selected_rows_preserves_two_original_id_rows_and_memberships() {
        // Rows are the fixed, source-ordered disconnected cycles from the
        // core fast-ring owner regression, including non-dense bond IDs.
        let atoms = vec![
            vec![AtomId::new(2), AtomId::new(1), AtomId::new(0)],
            vec![AtomId::new(5), AtomId::new(4), AtomId::new(3)],
        ];
        let bonds = vec![
            vec![BondId::new(3), BondId::new(1), BondId::new(5)],
            vec![BondId::new(2), BondId::new(0), BondId::new(4)],
        ];
        let originals = (atoms.clone(), bonds.clone());
        let rings = ring_info_from_selected_rows(6, 6, &atoms, &bonds).unwrap();
        assert_eq!(rings.atom_rings(), atoms);
        assert_eq!(rings.bond_rings(), bonds);
        assert_eq!(rings.atom_members(AtomId::new(1)), &[0]);
        assert_eq!(rings.atom_members(AtomId::new(4)), &[1]);
        assert_eq!(rings.bond_members(BondId::new(5)), &[0]);
        assert_eq!(rings.bond_members(BondId::new(0)), &[1]);
        assert_eq!((atoms, bonds), originals);
    }

    #[test]
    fn ring_info_from_selected_rows_keeps_repeated_source_rows() {
        let atoms = vec![vec![AtomId::new(2), AtomId::new(1), AtomId::new(0)]; 2];
        let bonds = vec![vec![BondId::new(2), BondId::new(1), BondId::new(0)]; 2];
        let rings = ring_info_from_selected_rows(3, 3, &atoms, &bonds).unwrap();
        assert_eq!(rings.num_rings(), 2);
        assert_eq!(rings.atom_rings(), atoms);
        assert_eq!(rings.bond_rings(), bonds);
        assert_eq!(rings.atom_members(AtomId::new(2)), &[0, 1]);
        assert_eq!(rings.bond_members(BondId::new(2)), &[0, 1]);
    }

    #[test]
    fn ring_info_from_selected_rows_rejects_unpaired_or_out_of_range_rows() {
        let atom_row = vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)];
        let bond_row = vec![BondId::new(0), BondId::new(1), BondId::new(2)];
        assert_eq!(
            ring_info_from_selected_rows(3, 3, &[atom_row.clone()], &[]),
            Err(RingFindingError::Value {
                message: "ring atom/bond table length mismatch",
            })
        );
        assert_eq!(
            ring_info_from_selected_rows(3, 3, &[atom_row.clone()], &[vec![BondId::new(0)]]),
            Err(RingFindingError::Value {
                message: "length mismatch",
            })
        );
        let mut bad_atoms = atom_row.clone();
        bad_atoms[2] = AtomId::new(3);
        assert_eq!(
            ring_info_from_selected_rows(3, 3, &[bad_atoms], &[bond_row.clone()]),
            Err(RingFindingError::RingAtomOutOfRange {
                atom: 3,
                atom_count: 3,
            })
        );
        let mut bad_bonds = bond_row;
        bad_bonds[1] = BondId::new(3);
        assert_eq!(
            ring_info_from_selected_rows(3, 3, &[atom_row], &[bad_bonds]),
            Err(RingFindingError::RingBondOutOfRange {
                bond: 3,
                bond_count: 3,
            })
        );
    }
}

pub(crate) fn preserves_appended_terminal_hydrogen_ring_prefix(
    topology: &TopologyBlock,
    rings: &RingInfo,
) -> bool {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8:
    // AddHs.cpp preserves ring state, then appends single-bond terminal H.
    // RDKit✔️❌:   mol.clearComputedProps(false);
    // RDKit✔️❌:   unsigned int stopIdx = mol.getNumAtoms();
    // RDKit✔️❌:   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    // RDKit✔️❌:       newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️❌:       mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RingInfo.cpp source getters, already implemented by RingInfo above:
    // RDKit✔️❌: const RingInfo::INT_VECT &RingInfo::atomMembers(unsigned int idx) const {
    // RDKit✔️❌:   PRECONDITION(df_init, "RingInfo not initialized");
    // RDKit✔️❌:
    // RDKit✔️❌:   static const INT_VECT emptyVect;
    // RDKit✔️❌:   if (idx < d_atomMembers.size()) {
    // RDKit✔️❌:     return d_atomMembers[idx];
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return emptyVect;
    // RDKit✔️❌: }
    // RDKit✔️❌: const RingInfo::INT_VECT &RingInfo::bondMembers(unsigned int idx) const {
    // RDKit✔️❌:   PRECONDITION(df_init, "RingInfo not initialized");
    // RDKit✔️❌:
    // RDKit✔️❌:   static const INT_VECT emptyVect;
    // RDKit✔️❌:   if (idx < d_bondMembers.size()) {
    // RDKit✔️❌:     return d_bondMembers[idx];
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return emptyVect;
    // RDKit✔️❌: }
    // This Rust-only detached input proof does not infer rings or chemistry.
    // Exact suffix degree and neighbor-bond identity give a bijection of new
    // atoms to new edges; no allocation or mutation is needed. Table/member
    // IDs remain within the old source prefix. Complexity O(A+B+membership),
    // constant extra space; source getter is O(1), hence the cost-gap marker.
    let old_atoms = rings.atom_row_count();
    let old_bonds = rings.bond_row_count();
    let atom_count = topology.atoms.len();
    let bond_count = topology.bonds.len();
    if !rings.is_initialized()
        || old_atoms >= atom_count
        || old_bonds >= bond_count
        || atom_count - old_atoms != bond_count - old_bonds
    {
        return false;
    }
    if topology.bonds[..old_bonds]
        .iter()
        .any(|bond| bond.begin().index() >= old_atoms || bond.end().index() >= old_atoms)
    {
        return false;
    }
    for atom in &topology.atoms[old_atoms..] {
        let neighbors = topology.adjacency.neighbors_of(atom.id().index());
        if atom.atomic_number() != 1 || neighbors.len() != 1 {
            return false;
        }
        let neighbor = &neighbors[0];
        if neighbor.atom_index >= old_atoms
            || neighbor.bond.index() < old_bonds
            || neighbor.bond.index() >= bond_count
        {
            return false;
        }
        let edge = &topology.bonds[neighbor.bond.index()];
        if edge.order() != BondOrder::Single
            || !((edge.begin() == atom.id() && edge.end().index() == neighbor.atom_index)
                || (edge.end() == atom.id() && edge.begin().index() == neighbor.atom_index))
        {
            return false;
        }
    }
    if topology.bonds[old_bonds..].iter().any(|bond| {
        bond.order() != BondOrder::Single
            || (bond.begin().index() >= old_atoms) == (bond.end().index() >= old_atoms)
    }) {
        return false;
    }
    let ring_count = rings.atom_rings.len();
    if ring_count != rings.bond_rings.len()
        || rings
            .atom_rings
            .iter()
            .zip(&rings.bond_rings)
            .any(|(atoms, bonds)| {
                atoms.len() != bonds.len()
                    || atoms.iter().any(|atom| atom.index() >= old_atoms)
                    || bonds.iter().any(|bond| bond.index() >= old_bonds)
            })
        || rings
            .atom_members
            .iter()
            .chain(&rings.bond_members)
            .any(|row| row.iter().any(|member| *member >= ring_count))
        || rings
            .atom_ring_families
            .iter()
            .flatten()
            .any(|atom| atom.index() >= old_atoms)
        || rings
            .bond_ring_families
            .iter()
            .flatten()
            .any(|bond| bond.index() >= old_bonds)
    {
        return false;
    }
    true
}

#[cfg(test)]
mod preserved_terminal_hydrogen_prefix_tests {
    use super::*;
    use crate::{DoubleBondStereoError, PotentialStereoError, ValenceAssignment};
    use cosmolkit_model::{Atom, AtomSpec, BondSpec, Element};

    fn topology(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        let atoms = elements
            .iter()
            .enumerate()
            .map(|(id, element)| {
                let mut atom = Atom::from_spec(AtomId::new(id), AtomSpec::new(*element));
                atom.set_no_implicit(true);
                atom
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, (begin, end, order))| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(AtomId::new(*begin), AtomId::new(*end), *order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn alkene() -> TopologyBlock {
        topology(
            &[
                Element::F,
                Element::C,
                Element::C,
                Element::F,
                Element::H,
                Element::H,
            ],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (2, 5, BondOrder::Single),
            ],
        )
    }
    fn scalar(t: &TopologyBlock) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: vec![0; t.atoms.len()],
            implicit_hydrogens: vec![0; t.atoms.len()],
        }
    }
    fn rejected_by_all_three_guards(t: &TopologyBlock, r: &RingInfo) {
        let before = r.clone();
        let value = scalar(t);
        let (reason, actual, expected) = if r.atom_row_count() != t.atoms.len() {
            (
                "atom membership row count mismatch",
                r.atom_row_count(),
                t.atoms.len(),
            )
        } else {
            (
                "bond membership row count mismatch",
                r.bond_row_count(),
                t.bonds.len(),
            )
        };
        assert_eq!(
            crate::potential_stereo(&mut t.clone(), &value, r, &Default::default()).unwrap_err(),
            PotentialStereoError::InvalidRingInfo {
                reason,
                row: 0,
                value: actual,
                limit: expected
            }
        );
        assert_eq!(
            crate::set_double_bond_neighbor_directions(t.clone(), r, None).unwrap_err(),
            DoubleBondStereoError::RingRowCount {
                dimension: if reason.starts_with("atom") {
                    "atom"
                } else {
                    "bond"
                },
                actual,
                expected
            }
        );
        // The internal bond candidate reaches the original bond guard after
        // the same valid degree and used-H checks; earlier returns stay intact.
        assert_eq!(
            crate::potential_stereo::is_potential_bond(t, &value, r, &t.bonds[1]).unwrap_err(),
            PotentialStereoError::InvalidRingInfo {
                reason: "bond membership row count mismatch",
                row: 0,
                value: r.bond_row_count(),
                limit: t.bonds.len()
            }
        );
        assert_eq!(r, &before);
    }

    #[test]
    fn all_three_existing_consumers_read_preserved_prefix_without_mutation() {
        let t = alkene();
        let r = RingInfo::new(RingFindType::SymmSssr, 4, 3);
        let before = r.clone();
        let value = scalar(&t);
        assert!(preserves_appended_terminal_hydrogen_ring_prefix(&t, &r));
        assert!(crate::potential_stereo(&mut t.clone(), &value, &r, &Default::default()).is_ok());
        assert!(crate::potential_stereo::is_potential_bond(&t, &value, &r, &t.bonds[1]).unwrap());
        assert!(crate::set_double_bond_neighbor_directions(t.clone(), &r, None).is_ok());
        assert_eq!(r, before);
        assert!(r.atom_members(AtomId::new(4)).is_empty());
        assert!(r.bond_members(BondId::new(3)).is_empty());
    }

    #[test]
    fn preserved_cycle_and_family_tables_keep_complete_source_carrier_identity() {
        let t = topology(
            &[Element::C, Element::C, Element::C, Element::H],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let mut r = RingInfo::new(RingFindType::SymmSssr, 3, 3);
        r.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
        r.atom_ring_families = vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]];
        r.bond_ring_families = vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]];
        let before = r.clone();
        assert!(preserves_appended_terminal_hydrogen_ring_prefix(&t, &r));
        assert!(
            crate::potential_stereo(&mut t.clone(), &scalar(&t), &r, &Default::default()).is_ok()
        );
        assert!(crate::set_double_bond_neighbor_directions(t, &r, None).is_ok());
        assert_eq!(r, before);
    }

    #[test]
    fn malformed_heavy_and_oversized_rows_retain_all_original_guard_errors() {
        let t = alkene();
        for (atoms, bonds) in [(3, 3), (7, 3), (4, 6), (6, 4)] {
            let r = RingInfo::new(RingFindType::SymmSssr, atoms, bonds);
            assert!(!preserves_appended_terminal_hydrogen_ring_prefix(&t, &r));
            rejected_by_all_three_guards(&t, &r);
        }
    }

    #[test]
    fn source_table_and_membership_indices_must_stay_in_original_prefix() {
        let t = alkene();
        let source = RingInfo::new(RingFindType::SymmSssr, 4, 3);
        let mut suffix_atom = source.clone();
        suffix_atom.atom_rings = vec![vec![AtomId::new(4)]];
        suffix_atom.bond_rings = vec![vec![BondId::new(0)]];
        let mut suffix_bond = source.clone();
        suffix_bond.atom_rings = vec![vec![AtomId::new(0)]];
        suffix_bond.bond_rings = vec![vec![BondId::new(3)]];
        let mut bad_member = source.clone();
        bad_member.atom_members[0].push(0);
        let mut table_length = source.clone();
        table_length.atom_rings.push(vec![]);
        let mut table_size = source.clone();
        table_size.atom_rings.push(vec![AtomId::new(0)]);
        table_size.bond_rings.push(vec![]);
        let mut suffix_family = source.clone();
        suffix_family.atom_ring_families = vec![vec![AtomId::new(4)]];
        let mut suffix_bond_family = source.clone();
        suffix_bond_family.bond_ring_families = vec![vec![BondId::new(3)]];
        for r in [
            suffix_atom,
            suffix_bond,
            bad_member,
            table_length,
            table_size,
            suffix_family,
            suffix_bond_family,
        ] {
            assert!(!preserves_appended_terminal_hydrogen_ring_prefix(&t, &r));
            rejected_by_all_three_guards(&t, &r);
        }
    }

    #[test]
    fn terminal_suffix_requires_exact_old_neighbor_and_appended_edge_accounting() {
        let elements = [
            Element::F,
            Element::C,
            Element::C,
            Element::F,
            Element::H,
            Element::H,
        ];
        let r = RingInfo::new(RingFindType::SymmSssr, 4, 3);
        let cases = [
            // A suffix H-H component has no old-prefix neighbor.
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (4, 5, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
            // Degree-two suffix H and isolated suffix H have no bijection.
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (2, 4, BondOrder::Single),
            ],
            // Existing prefix edge touches suffix, even with terminal leaves.
            vec![
                (1, 4, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (0, 1, BondOrder::Single),
                (2, 5, BondOrder::Single),
            ],
            // Extra appended old-old edge violates suffix row accounting.
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (2, 5, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
            // Source AddHs always appends SINGLE bonds.
            vec![
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (1, 4, BondOrder::Double),
                (2, 5, BondOrder::Single),
            ],
        ];
        for edges in cases {
            let t = topology(&elements, &edges);
            assert!(!preserves_appended_terminal_hydrogen_ring_prefix(&t, &r));
            let before = r.clone();
            assert!(matches!(
                crate::set_double_bond_neighbor_directions(t.clone(), &r, None),
                Err(DoubleBondStereoError::RingRowCount {
                    dimension: "atom",
                    actual: 4,
                    expected: 6
                })
            ));
            assert!(matches!(
                crate::potential_stereo(&mut t.clone(), &scalar(&t), &r, &Default::default()),
                Err(PotentialStereoError::InvalidRingInfo {
                    reason: "atom membership row count mismatch",
                    value: 4,
                    limit: 6,
                    ..
                })
            ));
            assert_eq!(r, before);
        }
        let mut heavy = elements;
        heavy[4] = Element::C;
        let t = topology(
            &heavy,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (2, 5, BondOrder::Single),
            ],
        );
        assert!(!preserves_appended_terminal_hydrogen_ring_prefix(&t, &r));
        rejected_by_all_three_guards(&t, &r);
    }

    #[test]
    fn multiple_double_bonds_preserve_exact_source_prefix_and_full_row_results() {
        let t = topology(
            &[
                Element::F,
                Element::C,
                Element::C,
                Element::C,
                Element::C,
                Element::F,
                Element::H,
                Element::H,
                Element::H,
                Element::H,
            ],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 5, BondOrder::Single),
                (1, 6, BondOrder::Single),
                (2, 7, BondOrder::Single),
                (3, 8, BondOrder::Single),
                (4, 9, BondOrder::Single),
            ],
        );
        let fields = crate::assign_valence(&t, &Default::default()).unwrap();
        let prefix = RingInfo::new(RingFindType::SymmSssr, 6, 5);
        let full = RingInfo::new(RingFindType::SymmSssr, t.atoms.len(), t.bonds.len());
        let before = prefix.clone();
        let actual =
            crate::potential_stereo(&mut t.clone(), &fields, &prefix, &Default::default()).unwrap();
        let fully_sized =
            crate::potential_stereo(&mut t.clone(), &fields, &full, &Default::default()).unwrap();
        assert_eq!(actual, fully_sized);
        assert_eq!(actual.stereo.len(), 2);
        assert!(
            actual
                .stereo
                .iter()
                .all(|item| item.stereo_type == crate::PotentialStereoType::BondDouble)
        );
        for bond in [BondId::new(1), BondId::new(3)] {
            assert!(
                crate::potential_stereo::is_potential_bond(
                    &t,
                    &fields,
                    &prefix,
                    &t.bonds[bond.index()]
                )
                .unwrap()
            );
        }
        assert_eq!(prefix, before);
        assert_eq!(prefix.atom_row_count(), 6);
        assert_eq!(prefix.bond_row_count(), 5);
    }

    #[test]
    fn fully_sized_carriers_and_initialization_order_remain_unchanged() {
        let t = alkene();
        let full = RingInfo::new(RingFindType::SymmSssr, 6, 5);
        assert!(!preserves_appended_terminal_hydrogen_ring_prefix(&t, &full));
        assert!(
            crate::potential_stereo(&mut t.clone(), &scalar(&t), &full, &Default::default())
                .is_ok()
        );
        assert!(
            crate::potential_stereo::is_potential_bond(&t, &scalar(&t), &full, &t.bonds[1])
                .unwrap()
        );
        assert!(crate::set_double_bond_neighbor_directions(t.clone(), &full, None).is_ok());
        let mut uninitialized = RingInfo::new(RingFindType::SymmSssr, 4, 3);
        uninitialized.reset();
        assert!(!preserves_appended_terminal_hydrogen_ring_prefix(
            &t,
            &uninitialized
        ));
        assert_eq!(
            crate::set_double_bond_neighbor_directions(t.clone(), &uninitialized, None)
                .unwrap_err(),
            DoubleBondStereoError::RingInfoNotInitialized
        );
        assert_eq!(
            crate::potential_stereo::is_potential_bond(
                &t,
                &scalar(&t),
                &uninitialized,
                &t.bonds[1]
            )
            .unwrap_err(),
            PotentialStereoError::InvalidRingInfo {
                reason: "ring information is not initialized",
                row: 0,
                value: 0,
                limit: 0
            }
        );
    }
}

#[cfg(test)]
mod source586_sssr_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec, Element, PropertyValue};
    fn graph(n: usize, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..n)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b, o))| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn triangle(order: BondOrder) -> TopologyBlock {
        graph(
            3,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, order),
            ],
        )
    }
    fn old_state() -> RingInfo {
        let mut r = RingInfo::new(RingFindType::SymmSssr, 3, 3);
        r.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
        r.atom_ring_families = vec![vec![AtomId::new(0)]];
        r.relevant_cycle_count = Some(1);
        r
    }

    #[test]
    fn source586_none_and_actual_output_delegate_same_disconnected_ring_state() {
        let g = graph(
            6,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Single),
                (5, 3, BondOrder::Single),
            ],
        );
        let mut a = old_state();
        let mut b = old_state();
        let mut output = vec![vec![90, 91]];
        let n = find_sssr_with_source_outputs_from_parts(
            6,
            &g.bonds,
            &g.adjacency,
            &mut a,
            None,
            None,
            false,
            false,
        )
        .unwrap();
        let m = find_sssr_with_source_outputs_from_parts(
            6,
            &g.bonds,
            &g.adjacency,
            &mut b,
            None,
            Some(&mut output),
            false,
            false,
        )
        .unwrap();
        assert_eq!((n, m), (2, 2));
        assert_eq!(a, b);
        assert_eq!(a.find_type(), RingFindType::Sssr);
        assert_eq!(
            output,
            a.atom_rings()
                .iter()
                .map(|r| r.iter().map(|i| i.index()).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        );
        assert_eq!(output.len(), 2);
        assert!(output[0].iter().all(|i| *i < 3));
        assert!(output[1].iter().all(|i| *i >= 3));
        assert!(a.atom_ring_families.is_empty());
        assert_eq!(a.relevant_cycle_count, None);
    }

    #[test]
    fn source586_optional_flags_are_forwarded_without_changing_zero_policy() {
        for order in [
            BondOrder::Dative,
            BondOrder::DativeOne,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::Hydrogen,
            BondOrder::Zero,
        ] {
            let g = triangle(order);
            for d in [false, true] {
                for h in [false, true] {
                    let mut r = old_state();
                    let mut out = vec![vec![44]];
                    let count = find_sssr_with_source_outputs_from_parts(
                        3,
                        &g.bonds,
                        &g.adjacency,
                        &mut r,
                        None,
                        Some(&mut out),
                        d,
                        h,
                    )
                    .unwrap();
                    let expected = if order == BondOrder::Zero {
                        0
                    } else if order == BondOrder::Hydrogen {
                        i32::from(h)
                    } else {
                        i32::from(d)
                    };
                    assert_eq!(count, expected, "{order:?},d={d},h={h}");
                    assert_eq!(out.len(), expected as usize);
                    assert_eq!(r.atom_rings().len(), expected as usize);
                }
            }
        }
    }

    #[test]
    fn source586_property_failure_follows_output_clear_cache_reset_and_sssr_initialize() {
        let g = triangle(BondOrder::Single);
        let mut r = old_state();
        let mut output = vec![vec![70]];
        let mut p = MoleculeProperties::default();
        p.set_prop("extraRings", PropertyValue::UInt(77)).unwrap();
        p.set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let before = p.clone();
        assert!(matches!(
            find_sssr_with_source_outputs_from_parts(
                3,
                &g.bonds,
                &g.adjacency,
                &mut r,
                Some(&mut p),
                Some(&mut output),
                false,
                false
            ),
            Err(RingFindingError::MoleculeProperty(_))
        ));
        assert!(output.is_empty());
        assert!(r.is_initialized());
        assert_eq!(r.find_type(), RingFindType::Sssr);
        assert!(r.atom_rings().is_empty());
        assert!(r.atom_ring_families.is_empty());
        assert_eq!(r.relevant_cycle_count, None);
        assert_eq!(p, before);
    }

    #[test]
    fn source586_empty_acyclic_outputs_clear_old_rows_and_preserve_unrelated_properties() {
        for g in [
            graph(0, &[]),
            graph(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]),
        ] {
            let mut r = old_state();
            let mut out = vec![vec![2]];
            let mut p = MoleculeProperties::default();
            p.set_prop("keep", "retained").unwrap();
            p.set_prop("extraRings", PropertyValue::UInt(1)).unwrap();
            assert_eq!(
                find_sssr_with_source_outputs_from_parts(
                    g.atoms.len(),
                    &g.bonds,
                    &g.adjacency,
                    &mut r,
                    Some(&mut p),
                    Some(&mut out),
                    false,
                    false
                )
                .unwrap(),
                0
            );
            assert!(out.is_empty());
            assert!(r.atom_rings().is_empty());
            assert!(r.is_initialized());
            assert_eq!(r.find_type(), RingFindType::Sssr);
            assert_eq!(p.prop("extraRings"), None);
            assert_eq!(
                p.prop("keep"),
                Some(&PropertyValue::String("retained".into()))
            );
        }
    }
}

#[cfg(test)]
mod source_ring_snapshot_tests {
    use super::*;
    #[test]
    fn all_modeled_source_cache_metadata_roundtrips_without_reconstruction() {
        let r = RingInfo::from_persisted_components(
            true,
            RingFindType::SymmSssr,
            4,
            4,
            vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]],
            vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]],
            vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]],
            vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]],
            Some(1),
            vec![vec![false]],
            vec![0],
        )
        .unwrap();
        let original = r.clone();
        let state = r.into_source_snapshot();
        assert_eq!(state.find_type, 2);
        assert_eq!(state.relevant_cycle_count, Some(1));
        assert_eq!(RingInfo::from_source_snapshot(&state).unwrap(), original);
    }
    #[test]
    fn native_duplicate_memberships_and_preallocated_empty_rows_are_not_rebuilt() {
        let mut r = RingInfo::new(RingFindType::Sssr, 4, 7);
        r.add_ring(&[0, 0, 1], &[0, 0, 1]).unwrap();
        let original = r.clone();
        let state = r.into_source_snapshot();
        assert_eq!(state.bond_members[0], vec![0, 0]);
        assert_eq!(state.atom_members.len(), 4);
        assert_eq!(state.bond_members.len(), 7);
        assert_eq!(RingInfo::from_source_snapshot(&state).unwrap(), original);
    }
    #[test]
    fn uninitialized_source_constructor_state_stays_uninitialized_and_empty() {
        let source = cosmolkit_model::SourceRingInfo::default();
        let r = RingInfo::from_source_snapshot(&source).unwrap();
        assert!(!r.is_initialized());
        assert_eq!(r.persisted_find_type(), RingFindType::OtherOrUnknown);
        assert_eq!(r.into_source_snapshot(), source);
    }
    #[test]
    fn unknown_source_enum_is_structural_error_instead_of_other_fallback() {
        let mut source = cosmolkit_model::SourceRingInfo::default();
        source.find_type = 17;
        assert!(matches!(
            RingInfo::from_source_snapshot(&source),
            Err(RingFindingError::SourceRingFindType { tag: 17 })
        ));
    }
}
