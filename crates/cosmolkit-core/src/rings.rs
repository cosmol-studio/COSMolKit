use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::{AdjacencyList, AtomId, Bond, BondId, BondOrder, Molecule};

const MAX_BFSQ_SIZE: usize = 200_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingFindType {
    OtherOrUnknown,
    Fast,
    Sssr,
    SymmSssr,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RingFindingError {
    #[error("{message}")]
    Value { message: &'static str },
    #[error("expected bond not found between atom {begin} and atom {end}")]
    ExpectedBondNotFound { begin: AtomId, end: AtomId },
    #[error("ring atom index {atom} is out of range for {atom_count} atoms")]
    RingAtomOutOfRange { atom: usize, atom_count: usize },
    #[error("ring bond index {bond} is out of range for {bond_count} bonds")]
    RingBondOutOfRange { bond: usize, bond_count: usize },
    #[error("unsupported ring info branch: {reason}")]
    UnsupportedBranch { reason: &'static str },
    #[error(transparent)]
    Adjacency(#[from] crate::AdjacencyError),
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

    pub(crate) fn initialize(&mut self, find_type: RingFindType) {
        // BEGIN RDKIT CPP FUNCTION RingInfo::initialize
        // RDKit✔️✔️: void RingInfo::initialize(RDKit::FIND_RING_TYPE ringType) {
        // RDKit✔️✔️:   df_init = true;
        self.initialized = true;
        // RDKit✔️✔️:   df_find_type_type = ringType;
        self.find_type = find_type;
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION RingInfo::initialize
    }

    #[allow(dead_code)]
    pub(crate) fn reset(&mut self) {
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
    pub(crate) fn preallocate(&mut self, atom_count: usize, bond_count: usize) {
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
        if atom_indices.len() != bond_indices.len() {
            return Err(RingFindingError::Value {
                message: "length mismatch",
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

    pub fn num_relevant_cycles(&self) -> Result<usize, RingFindingError> {
        // BEGIN RDKIT CPP FUNCTION RingInfo::numRelevantCycles
        // RDKit❗✔️: unsigned int RingInfo::numRelevantCycles() const {
        // RDKit❗✔️:   PRECONDITION(df_init, "RingInfo not initialized");
        // RDKit❗✔️:   return rdcast<unsigned int>(RDL_getNofRC(dp_urfData.get()));
        // RDKit❗✔️: };
        // END RDKIT CPP FUNCTION RingInfo::numRelevantCycles
        debug_assert!(self.initialized, "RingInfo not initialized");
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

    pub(crate) fn add_ring(
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

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingSearchResult {
    rings: Vec<Vec<usize>>,
    extra_rings: Vec<Vec<usize>>,
}

struct RingSearchContext<'a> {
    molecule: &'a Molecule,
    adjacency: AdjacencyList,
}

impl<'a> RingSearchContext<'a> {
    fn new(molecule: &'a Molecule) -> Result<Self, RingFindingError> {
        Ok(Self {
            molecule,
            adjacency: AdjacencyList::try_from_topology(molecule.num_atoms(), molecule.bonds())?,
        })
    }

    fn atom_count(&self) -> usize {
        self.molecule.num_atoms()
    }

    fn bond_count(&self) -> usize {
        self.molecule.num_bonds()
    }

    fn bonds(&self) -> &'a [Bond] {
        self.molecule.bonds()
    }

    fn neighbors(&self, atom: usize) -> &[crate::NeighborRef] {
        self.adjacency.neighbors_of(atom)
    }

    fn bond_between_atoms(&self, begin: usize, end: usize) -> Option<BondId> {
        self.neighbors(begin)
            .iter()
            .find(|neighbor| neighbor.atom_index == end)
            .map(|neighbor| neighbor.bond)
    }
}

pub fn symmetrize_sssr(molecule: &Molecule) -> Result<RingInfo, RingFindingError> {
    symmetrize_sssr_with_options(molecule, false, false)
}

pub fn symmetrize_sssr_with_options(
    molecule: &Molecule,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::symmetrizeSSSR
    // RDKit❗✔️: int symmetrizeSSSR(ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
    // RDKit❗✔️:                    bool includeHydrogenBonds) {
    // RDKit❗✔️:   res.clear();
    // RDKit❗✔️:   VECT_INT_VECT sssrs;
    let context = RingSearchContext::new(molecule)?;
    // RDKit❗✔️:   findSSSR(mol, sssrs, includeDativeBonds, includeHydrogenBonds);
    let sssr = find_sssr_internal(&context, include_dative_bonds, include_hydrogen_bonds)?;
    // RDKit❗✔️:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);
    let mut info = RingInfo::new(
        RingFindType::SymmSssr,
        context.atom_count(),
        context.bond_count(),
    );
    // RDKit❗✔️:   res.reserve(sssrs.size());
    // RDKit❗✔️:   for (const auto &r : sssrs) {
    // RDKit❗✔️:     res.emplace_back(r);
    // RDKit❗✔️:   }
    let mut rings = sssr.rings.clone();
    // RDKit❗✔️:   if (!mol.hasProp(common_properties::extraRings)) {
    // RDKit❗✔️:     return rdcast<int>(res.size());
    // RDKit❗✔️:   }
    if !sssr.extra_rings.is_empty() {
        // RDKit❗✔️:   VECT_INT_VECT bondsssrs;
        // RDKit❗✔️:   RingUtils::convertToBonds(sssrs, bondsssrs, mol);
        let bond_sssrs = convert_rings_to_bonds(&context, &sssr.rings)?;
        // RDKit❗✔️:   std::vector<int> bondCounts(mol.getNumBonds(), 0);
        // RDKit❗✔️:   for (const auto &r : bondsssrs) {
        // RDKit❗✔️:     for (const auto &b : r) {
        // RDKit❗✔️:       bondCounts[b] += 1;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        let mut bond_counts = vec![0i32; context.bond_count()];
        for ring in &bond_sssrs {
            for &bond in ring {
                bond_counts[bond] += 1;
            }
        }
        // RDKit❗✔️:   for (auto &extraAtomRing : extras) {
        for extra_atom_ring in &sssr.extra_rings {
            // RDKit❗✔️:     RingUtils::convertToBonds(extraAtomRing, extraRing, mol);
            let extra_ring = convert_to_bonds(&context, extra_atom_ring)?;
            // RDKit❗✔️:     for (auto &ring : bondsssrs) {
            for ring in &bond_sssrs {
                // RDKit❗✔️:       if (ring.size() != extraRing.size()) {
                // RDKit❗✔️:         continue;
                // RDKit❗✔️:       }
                if ring.len() != extra_ring.len() {
                    continue;
                }
                let mut share_bond = false;
                let mut replaces_all_unique_bonds = true;
                // RDKit❗✔️:       for (auto &bondID : ring) {
                for &bond_id in ring {
                    // RDKit❗✔️:         const int bondCount = bondCounts[bondID];
                    let bond_count = bond_counts[bond_id];
                    // RDKit❗✔️:         if (bondCount == 1 || !shareBond) {
                    if bond_count == 1 || !share_bond {
                        // RDKit❗✔️:           auto position = find(extraRing.begin(), extraRing.end(), bondID);
                        if extra_ring.contains(&bond_id) {
                            // RDKit❗✔️:           if (position != extraRing.end()) {
                            // RDKit❗✔️:             shareBond = true;
                            share_bond = true;
                        } else if bond_count == 1 {
                            // RDKit❗✔️:           } else if (bondCount == 1) {
                            // RDKit❗✔️:             replacesAllUniqueBonds = false;
                            replaces_all_unique_bonds = false;
                        }
                    }
                }
                // RDKit❗✔️:       if (shareBond && replacesAllUniqueBonds) {
                // RDKit❗✔️:         res.push_back(extraAtomRing);
                // RDKit❗✔️:         FindRings::storeRingInfo(mol, extraAtomRing);
                // RDKit❗✔️:         break;
                // RDKit❗✔️:       }
                if share_bond && replaces_all_unique_bonds {
                    rings.push(extra_atom_ring.clone());
                    break;
                }
            }
        }
    }
    store_rings_info(&context, &rings, &mut info)?;
    // RDKit❗✔️:   return rdcast<int>(res.size());
    // RDKit❗✔️: }
    Ok(info)
}

pub fn find_sssr(molecule: &Molecule) -> Result<RingInfo, RingFindingError> {
    let context = RingSearchContext::new(molecule)?;
    let result = find_sssr_internal(&context, false, false)?;
    let mut info = RingInfo::new(
        RingFindType::Sssr,
        context.atom_count(),
        context.bond_count(),
    );
    store_rings_info(&context, &result.rings, &mut info)?;
    Ok(info)
}

pub fn find_ring_families(
    molecule: &Molecule,
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
) -> Result<RingInfo, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::findRingFamilies
    // RDKit❗✔️: #ifdef RDK_USE_URF
    // RDKit❗✔️: void findRingFamilies(const ROMol &mol, bool includeDativeBonds,
    // RDKit❗✔️:                       bool includeHydrogenBonds) {
    // RDKit❗✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit❗✔️:     // return if we've done this before
    // RDKit❗✔️:     if (mol.getRingInfo()->areRingFamiliesInitialized()) {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     mol.getRingInfo()->initialize();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDL_graph *graph = RDL_initNewGraph(mol.getNumAtoms());
    let mut graph = cosmolkit_ringdecomposer::Graph::new(molecule.num_atoms());
    // RDKit❗✔️:   for (auto cbi : mol.bonds()) {
    for bond in molecule.bonds() {
        // RDKit❗✔️:     if (auto bt = cbi->getBondType();
        // RDKit❗✔️:         bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
        // RDKit❗✔️:         (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
        if bond.order() == BondOrder::Zero
            || (!include_dative_bonds && is_dative_bond_order(bond.order()))
            || (!include_hydrogen_bonds && bond.order() == BondOrder::Hydrogen)
        {
            continue;
        }
        // RDKit❗✔️:     RDL_addUEdge(graph, cbi->getBeginAtomIdx(), cbi->getEndAtomIdx());
        graph.add_undirected_edge(bond.begin().index(), bond.end().index())?;
    }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDL_data *urfdata = RDL_calculate(graph);
    let decomposition = cosmolkit_ringdecomposer::RingDecomposition::calculate(graph)?;
    // RDKit❗✔️:   for (unsigned int i = 0; i < RDL_getNofURF(urfdata); ++i) {
    let mut info = RingInfo::new(
        RingFindType::SymmSssr,
        molecule.num_atoms(),
        molecule.num_bonds(),
    );
    info.relevant_cycle_count = Some(decomposition.relevant_cycle_count() as usize);
    for urf in decomposition.urfs() {
        // RDKit❗✔️:     RDL_node *nodes = nullptr;
        // RDKit❗✔️:     unsigned nNodes = RDL_getNodesForURF(urfdata, i, &nodes);
        let atom_indices = urf
            .nodes()
            .iter()
            .map(|node| AtomId::new(*node).index())
            .collect::<Vec<_>>();
        // RDKit❗✔️:     RDL_edge *edges = nullptr;
        // RDKit❗✔️:     unsigned nEdges = RDL_getEdgesForURF(urfdata, i, &edges);
        let bond_indices = urf
            .edges()
            .iter()
            .map(|edge| edge.index())
            .collect::<Vec<_>>();
        // RDKit❗✔️:     mol.getRingInfo()->addRingFamily(nvect, evect);
        info.add_ring_family(&atom_indices, &bond_indices)?;
    }
    // RDKit❗✔️: }
    // RDKit❗✔️: #else
    // RDKit❗✔️: void findRingFamilies(const ROMol &mol, bool includeDativeBonds,
    // RDKit❗✔️:                       bool includeHydrogenBonds) {
    // RDKit❗✔️:   BOOST_LOG(rdErrorLog)
    // RDKit❗✔️:       << "This version of the RDKit was built without URF support" << std::endl;
    // RDKit❗✔️: }
    // RDKit❗✔️: #endif
    // END RDKIT CPP FUNCTION MolOps::findRingFamilies
    Ok(info)
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
) -> Result<RingSearchResult, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::findSSSR
    // RDKit❗✔️: int findSSSR(const ROMol &mol, VECT_INT_VECT &res, bool includeDativeBonds,
    // RDKit❗✔️:              bool includeHydrogenBonds) {
    // RDKit❗✔️:   res.resize(0);
    let mut res = Vec::new();
    // RDKit❗✔️:   boost::dynamic_bitset<> activeBonds(nbnds);
    // RDKit❗✔️:   activeBonds.set();
    let mut active_bonds = vec![true; context.bond_count()];
    // RDKit❗✔️:   for (auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (auto bt = bond->getBondType();
    // RDKit❗✔️:         bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
    // RDKit❗✔️:         (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
    // RDKit❗✔️:       activeBonds[bond->getIdx()] = 0;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    for bond in context.bonds() {
        if is_inactive_ring_bond(bond.order(), include_dative_bonds, include_hydrogen_bonds) {
            active_bonds[bond.id().index()] = false;
        }
    }
    // RDKit❗✔️:   INT_VECT atomDegrees(nats);
    // RDKit❗✔️:   INT_VECT atomDegreesWithZeroOrderBonds(nats);
    let mut atom_degrees = vec![0i32; context.atom_count()];
    let mut atom_degrees_with_zero_order_bonds = vec![0i32; context.atom_count()];
    // RDKit❗✔️:   for (unsigned int i = 0; i < nats; ++i) {
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
    // RDKit❗✔️:   mol.clearProp(common_properties::extraRings);
    let mut extra_rings = Vec::new();
    // RDKit❗✔️:   RINGINVAR_SET invars;
    // RDKit❗✔️:   boost::dynamic_bitset<> ringAtoms(nats);
    // RDKit❗✔️:   boost::dynamic_bitset<> ringBonds(nbnds);
    let mut invars = BTreeSet::new();
    let mut ring_atoms = vec![false; context.atom_count()];
    let mut ring_bonds = vec![false; context.bond_count()];
    // RDKit❗✔️:   VECT_INT_VECT frags;
    // RDKit❗✔️:   getMolFrags(mol, frags);
    let frags = molecule_fragments(context);
    // RDKit❗✔️:   for (const auto &curFrag : frags) {
    for cur_frag in frags {
        // RDKit❗✔️:     if (curFrag.size() < 3) {
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
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
        // RDKit❗✔️:     while (nAtomsDone <= curFrag.size() - 3) {
        while atoms_done <= cur_frag.len().saturating_sub(3) {
            // RDKit❗✔️:       while (!changed.empty()) {
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
            } else if atoms_done <= cur_frag.len().saturating_sub(3) {
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
        let nexpt = isize::try_from(nbnds).map_err(|_| RingFindingError::Value {
            message: "bond count out of range",
        })? - isize::try_from(cur_frag.len()).map_err(|_| RingFindingError::Value {
            message: "fragment atom count out of range",
        })? + 1;
        let mut ssize = isize::try_from(frag_res.len()).map_err(|_| RingFindingError::Value {
            message: "ring count out of range",
        })?;
        if ssize < nexpt {
            let mut possible_bonds = Vec::new();
            for i in 0..nbnds {
                if !ring_bonds[i] {
                    let bond = &context.bonds()[i];
                    if ring_atoms[bond.begin().index()] && ring_atoms[bond.end().index()] {
                        possible_bonds.push(bond.id());
                        break;
                    }
                }
            }
            let mut dead_bonds = vec![false; context.bond_count()];
            while let Some(possible_bond) = possible_bonds.first().copied() {
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
                possible_bonds.clear();
                for i in 0..nbnds {
                    if !ring_bonds[i] {
                        let bond = &context.bonds()[i];
                        if !dead_bonds[bond.id().index()]
                            && ring_atoms[bond.begin().index()]
                            && ring_atoms[bond.end().index()]
                        {
                            possible_bonds.push(bond.id());
                            break;
                        }
                    }
                }
            }
            ssize = isize::try_from(frag_res.len()).map_err(|_| RingFindingError::Value {
                message: "ring count out of range",
            })?;
            if ssize < nexpt {
                let fast = fast_find_rings_internal(context)?;
                return Ok(RingSearchResult {
                    rings: fast,
                    extra_rings: Vec::new(),
                });
            }
        }
        if ssize > nexpt {
            let extras = remove_extra_rings(context, &mut frag_res)?;
            extra_rings.extend(extras);
        }
        res.extend(frag_res);
    }
    // RDKit❗✔️:   FindRings::storeRingsInfo(mol, res);
    // RDKit❗✔️:   return rdcast<int>(res.size());
    // RDKit❗✔️: }
    Ok(RingSearchResult {
        rings: res,
        extra_rings,
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

fn compute_ring_invariant(ring: &[usize]) -> Vec<usize> {
    let mut invariant = ring.to_vec();
    invariant.sort_unstable();
    invariant.dedup();
    invariant
}

fn convert_to_bonds(
    context: &RingSearchContext<'_>,
    ring: &[usize],
) -> Result<Vec<usize>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION RingUtils::convertToBonds
    // RDKit✔️✔️: void convertToBonds(const INT_VECT &ring, INT_VECT &bondRing,
    // RDKit✔️✔️:                     const ROMol &mol) {
    // RDKit✔️✔️:   const auto rsiz = rdcast<unsigned int>(ring.size());
    // RDKit✔️✔️:   bondRing.resize(rsiz);
    let ring_size = ring.len();
    let mut bond_ring = vec![0usize; ring_size];
    if ring_size == 0 {
        return Ok(bond_ring);
    }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < (rsiz - 1); i++) {
    for i in 0..ring_size - 1 {
        // RDKit✔️✔️:     const Bond *bnd = mol.getBondBetweenAtoms(ring[i], ring[i + 1]);
        // RDKit✔️✔️:     if (!bnd) {
        // RDKit✔️✔️:       throw ValueErrorException("expected bond not found");
        // RDKit✔️✔️:     }
        let Some(bond) = context.bond_between_atoms(ring[i], ring[i + 1]) else {
            return Err(RingFindingError::ExpectedBondNotFound {
                begin: AtomId::new(ring[i]),
                end: AtomId::new(ring[i + 1]),
            });
        };
        // RDKit✔️✔️:     bondRing[i] = bnd->getIdx();
        bond_ring[i] = bond.index();
    }
    // RDKit✔️✔️:   const Bond *bnd = mol.getBondBetweenAtoms(ring[rsiz - 1], ring[0]);
    // RDKit✔️✔️:   if (!bnd) {
    // RDKit✔️✔️:     throw ValueErrorException("expected bond not found");
    // RDKit✔️✔️:   }
    let Some(bond) = context.bond_between_atoms(ring[ring_size - 1], ring[0]) else {
        return Err(RingFindingError::ExpectedBondNotFound {
            begin: AtomId::new(ring[ring_size - 1]),
            end: AtomId::new(ring[0]),
        });
    };
    // RDKit✔️✔️:   bondRing[rsiz - 1] = bnd->getIdx();
    // RDKit✔️✔️: }
    bond_ring[ring_size - 1] = bond.index();
    Ok(bond_ring)
}

fn convert_rings_to_bonds(
    context: &RingSearchContext<'_>,
    rings: &[Vec<usize>],
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    let mut bond_rings = Vec::with_capacity(rings.len());
    for ring in rings {
        bond_rings.push(convert_to_bonds(context, ring)?);
    }
    Ok(bond_rings)
}

fn store_rings_info(
    context: &RingSearchContext<'_>,
    rings: &[Vec<usize>],
    info: &mut RingInfo,
) -> Result<(), RingFindingError> {
    for ring in rings {
        let bond_indices = convert_to_bonds(context, ring)?;
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
    // RDKit❗✔️: int smallestRingsBfs(const ROMol &mol, int root, VECT_INT_VECT &rings,
    // RDKit❗✔️:                      boost::dynamic_bitset<> &activeBonds,
    // RDKit❗✔️:                      INT_VECT *forbidden = nullptr) {
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
    // RDKit❗✔️:   return rdcast<unsigned int>(rings.size());
    // RDKit❗✔️: }
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
    invars: &mut BTreeSet<Vec<usize>>,
    dup_map: &BTreeMap<usize, Vec<usize>>,
    dup_d2_candidates: &BTreeMap<Vec<usize>, Vec<usize>>,
    atom_degrees: &[i32],
    active_bonds: &[bool],
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findSSSRforDupCands
    // RDKit❗✔️: void findSSSRforDupCands(const ROMol &mol, VECT_INT_VECT &res,
    // RDKit❗✔️:                          RINGINVAR_SET &invars, const INT_INT_VECT_MAP &dupMap,
    // RDKit❗✔️:                          const RINGINVAR_INT_VECT_MAP &dupD2Cands,
    // RDKit❗✔️:                          INT_VECT &atomDegrees,
    // RDKit❗✔️:                          const boost::dynamic_bitset<> &activeBonds) {
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
                    let invariant = compute_ring_invariant(&ring);
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
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FindRings::findSSSRforDupCands
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn find_rings_d2_nodes(
    context: &RingSearchContext<'_>,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<Vec<usize>>,
    d2nodes: &[usize],
    atom_degrees: &mut [i32],
    active_bonds: &mut [bool],
    ring_bonds: &mut [bool],
    ring_atoms: &mut [bool],
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findRingsD2nodes
    // RDKit❗✔️: void findRingsD2nodes(const ROMol &tMol, VECT_INT_VECT &res,
    // RDKit❗✔️:                       RINGINVAR_SET &invars, const INT_VECT &d2nodes,
    // RDKit❗✔️:                       INT_VECT &atomDegrees,
    // RDKit❗✔️:                       boost::dynamic_bitset<> &activeBonds,
    // RDKit❗✔️:                       boost::dynamic_bitset<> &ringBonds,
    // RDKit❗✔️:                       boost::dynamic_bitset<> &ringAtoms) {
    // RDKit✔️✔️:   RINGINVAR_INT_VECT_MAP dupD2Cands;
    // RDKit✔️✔️:   INT_INT_VECT_MAP dupMap;
    let mut dup_d2_candidates: BTreeMap<Vec<usize>, Vec<usize>> = BTreeMap::new();
    let mut dup_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    // RDKit✔️✔️:   for (auto &cand : d2nodes) {
    for &candidate in d2nodes {
        // RDKit✔️✔️:     VECT_INT_VECT srings;
        // RDKit✔️✔️:     smallestRingsBfs(tMol, cand, srings, activeBonds);
        let smallest_rings = smallest_rings_bfs(context, candidate, active_bonds, None)?;
        // RDKit✔️✔️:     for (const auto &nring : srings) {
        for ring in &smallest_rings {
            // RDKit✔️✔️:       auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
            let invariant = compute_ring_invariant(ring);
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
    // RDKit❗✔️: }
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
    invars: &mut BTreeSet<Vec<usize>>,
    candidate: usize,
    active_bonds: &[bool],
) -> Result<(), RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findRingsD3Node
    // RDKit❗✔️: void findRingsD3Node(const ROMol &tMol, VECT_INT_VECT &res,
    // RDKit❗✔️:                      RINGINVAR_SET &invars, int cand, INT_VECT &,
    // RDKit❗✔️:                      boost::dynamic_bitset<> &activeBonds) {
    // RDKit✔️✔️:   VECT_INT_VECT srings;
    // RDKit✔️✔️:   auto nsmall = smallestRingsBfs(tMol, cand, srings, activeBonds);
    let smallest = smallest_rings_bfs(context, candidate, active_bonds, None)?;
    let nsmall = smallest.len();
    // RDKit✔️✔️:   for (const auto &nring : srings) {
    for ring in &smallest {
        // RDKit✔️✔️:     auto invr = RingUtils::computeRingInvariant(nring, tMol.getNumAtoms());
        // RDKit✔️✔️:     if (invars.find(invr) == invars.end()) {
        let invariant = compute_ring_invariant(ring);
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
                let invariant = compute_ring_invariant(&ring);
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
                let invariant = compute_ring_invariant(&ring);
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
                let invariant = compute_ring_invariant(&ring);
                if !invars.contains(&invariant) {
                    res.push(ring);
                    invars.insert(invariant);
                }
            }
        }
    }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FindRings::findRingsD3Node
    Ok(())
}

fn remove_extra_rings(
    context: &RingSearchContext<'_>,
    res: &mut Vec<Vec<usize>>,
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::removeExtraRings
    // RDKit❗✔️: void removeExtraRings(VECT_INT_VECT &res, unsigned int, const ROMol &mol) {
    // RDKit✔️✔️:   std::sort(res.begin(), res.end(), compRingSize);
    res.sort_by_key(Vec::len);
    // RDKit✔️✔️:   VECT_INT_VECT brings;
    // RDKit✔️✔️:   RingUtils::convertToBonds(res, brings, mol);
    let bond_rings = convert_rings_to_bonds(context, res)?;
    // RDKit✔️✔️:   std::vector<boost::dynamic_bitset<>> bitBrings;
    // RDKit✔️✔️:   bitBrings.reserve(brings.size());
    let bit_bond_rings = bond_rings
        .iter()
        .map(|ring| ring.iter().copied().collect::<BTreeSet<_>>())
        .collect::<Vec<_>>();
    // RDKit✔️✔️:   boost::dynamic_bitset<> availRings(res.size());
    // RDKit✔️✔️:   availRings.set();
    // RDKit✔️✔️:   boost::dynamic_bitset<> keepRings(res.size());
    // RDKit✔️✔️:   boost::dynamic_bitset<> munion(mol.getNumBonds());
    let mut available = vec![true; res.len()];
    let mut keep = vec![false; res.len()];
    let mut union = BTreeSet::new();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < res.size(); ++i) {
    for i in 0..res.len() {
        // RDKit✔️✔️:     if (bitBrings[i].is_subset_of(munion)) {
        // RDKit✔️✔️:       availRings.set(i, 0);
        // RDKit✔️✔️:     }
        if bit_bond_rings[i].is_subset(&union) {
            available[i] = false;
        }
        // RDKit✔️✔️:     if (!availRings[i]) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if !available[i] {
            continue;
        }
        // RDKit✔️✔️:     munion |= bitBrings[i];
        // RDKit✔️✔️:     keepRings.set(i);
        union.extend(bit_bond_rings[i].iter().copied());
        keep[i] = true;
        // RDKit✔️✔️:     boost::dynamic_bitset<> consider(res.size());
        let mut consider = vec![false; res.len()];
        // RDKit✔️✔️:     for (unsigned int j = i + 1; j < res.size(); ++j) {
        for j in i + 1..res.len() {
            // RDKit✔️✔️:       if (availRings[j] && (brings[j].size() == brings[i].size())) {
            // RDKit✔️✔️:         consider.set(j);
            // RDKit✔️✔️:       }
            if available[j] && bond_rings[j].len() == bond_rings[i].len() {
                consider[j] = true;
            }
        }
        // RDKit✔️✔️:     while (consider.any()) {
        while consider.iter().any(|flag| *flag) {
            // RDKit✔️✔️:       unsigned int bestJ = i + 1;
            // RDKit✔️✔️:       int bestOverlap = -1;
            let mut best_j = i + 1;
            let mut best_overlap = -1isize;
            // RDKit✔️✔️:       for (unsigned int j = i + 1;
            // RDKit✔️✔️:            j < res.size() && brings[j].size() == brings[i].size(); ++j) {
            for j in i + 1..res.len() {
                if bond_rings[j].len() != bond_rings[i].len() {
                    break;
                }
                if !consider[j] || !available[j] {
                    continue;
                }
                // RDKit✔️✔️:         workspace = bitBrings[j];
                // RDKit✔️✔️:         workspace &= munion;
                // RDKit✔️✔️:         int overlap = rdcast<int>(workspace.count());
                let overlap = bit_bond_rings[j].intersection(&union).count() as isize;
                // RDKit✔️✔️:         if (overlap > bestOverlap) {
                // RDKit✔️✔️:           bestOverlap = overlap;
                // RDKit✔️✔️:           bestJ = j;
                // RDKit✔️✔️:         }
                if overlap > best_overlap {
                    best_overlap = overlap;
                    best_j = j;
                }
            }
            // RDKit✔️✔️:       consider.set(bestJ, 0);
            consider[best_j] = false;
            // RDKit✔️✔️:       if (bitBrings[bestJ].is_subset_of(munion)) {
            if bit_bond_rings[best_j].is_subset(&union) {
                // RDKit✔️✔️:         availRings.set(bestJ, 0);
                available[best_j] = false;
            } else {
                // RDKit✔️✔️:       } else {
                // RDKit✔️✔️:         keepRings.set(bestJ);
                // RDKit✔️✔️:         availRings.set(bestJ, 0);
                // RDKit✔️✔️:         munion |= bitBrings[bestJ];
                keep[best_j] = true;
                available[best_j] = false;
                union.extend(bit_bond_rings[best_j].iter().copied());
            }
        }
    }
    // RDKit✔️✔️:   VECT_INT_VECT extras;
    // RDKit✔️✔️:   VECT_INT_VECT temp = res;
    // RDKit✔️✔️:   res.resize(0);
    let old = std::mem::take(res);
    let mut extras = Vec::new();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < temp.size(); i++) {
    for (index, ring) in old.into_iter().enumerate() {
        // RDKit✔️✔️:     if (keepRings[i]) {
        if keep[index] {
            // RDKit✔️✔️:       res.push_back(temp[i]);
            res.push(ring);
        } else {
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       extras.push_back(temp[i]);
            extras.push(ring);
        }
    }
    // RDKit✔️✔️:   molExtras.insert(molExtras.end(), extras.begin(), extras.end());
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FindRings::removeExtraRings
    Ok(extras)
}

fn atom_search_bfs(
    context: &RingSearchContext<'_>,
    start_atom: usize,
    end_atom: usize,
    ring_atoms: &[bool],
    invars: &BTreeSet<Vec<usize>>,
) -> Result<Option<Vec<usize>>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::_atomSearchBFS
    // RDKit❗✔️: bool _atomSearchBFS(const ROMol &tMol, unsigned int startAtomIdx,
    // RDKit❗✔️:                     unsigned int endAtomIdx, boost::dynamic_bitset<> &ringAtoms,
    // RDKit❗✔️:                     INT_VECT &res, RINGINVAR_SET &invars) {
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
                    let invariant = compute_ring_invariant(&new_path);
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
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FindRings::_atomSearchBFS
    Ok(None)
}

fn find_ring_connecting_atoms(
    context: &RingSearchContext<'_>,
    bond: BondId,
    res: &mut Vec<Vec<usize>>,
    invars: &mut BTreeSet<Vec<usize>>,
    ring_bonds: &mut [bool],
    ring_atoms: &mut [bool],
) -> Result<bool, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION FindRings::findRingConnectingAtoms
    // RDKit❗✔️: bool findRingConnectingAtoms(const ROMol &tMol, const Bond *bond,
    // RDKit❗✔️:                              VECT_INT_VECT &res, RINGINVAR_SET &invars,
    // RDKit❗✔️:                              boost::dynamic_bitset<> &ringBonds,
    // RDKit❗✔️:                              boost::dynamic_bitset<> &ringAtoms) {
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
    let invariant = compute_ring_invariant(&ring);
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
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FindRings::findRingConnectingAtoms
    Ok(true)
}

fn molecule_fragments(context: &RingSearchContext<'_>) -> Vec<Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags
    // RDKit❗✔️: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
    // RDKit❗✔️:   unsigned int natms = mol.getNumAtoms();
    // RDKit❗✔️:   mapping.resize(natms);
    // RDKit❗✔️:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
    // RDKit❗✔️:                : 0;
    // RDKit❗✔️: };
    // RDKit❗✔️: unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
    // RDKit✔️✔️:   frags.clear();
    // RDKit✔️✔️:   INT_VECT mapping;
    // RDKit❗✔️:   getMolFrags(mol, mapping);
    let mut seen = vec![false; context.atom_count()];
    let mut fragments = Vec::new();
    // RDKit✔️✔️:   INT_INT_VECT_MAP comMap;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    for start in 0..context.atom_count() {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        let mut fragment = Vec::new();
        while let Some(atom) = queue.pop_front() {
            fragment.push(atom);
            for neighbor in context.neighbors(atom) {
                if !seen[neighbor.atom_index] {
                    seen[neighbor.atom_index] = true;
                    queue.push_back(neighbor.atom_index);
                }
            }
        }
        fragments.push(fragment);
    }
    // RDKit✔️✔️:     int mi = mapping[i];
    // RDKit✔️✔️:     if (comMap.find(mi) == comMap.end()) {
    // RDKit✔️✔️:       INT_VECT comp;
    // RDKit✔️✔️:       comMap[mi] = comp;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     comMap[mi].push_back(i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    // RDKit✔️✔️:     frags.push_back((*mci).second);
    // RDKit✔️✔️:   }
    // RDKit❗✔️:   return rdcast<unsigned int>(frags.size());
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags
    fragments
}

fn fast_find_rings_internal(
    context: &RingSearchContext<'_>,
) -> Result<Vec<Vec<usize>>, RingFindingError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::fastFindRings
    // RDKit❗✔️: void fastFindRings(const ROMol &mol) {
    // RDKit❌❌:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit❌❌:     mol.getRingInfo()->reset();
    // RDKit❌❌:   }
    // RDKit❌❌:   mol.getRingInfo()->initialize(FIND_RING_TYPE_FAST);
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
    // RDKit❌❌:   FindRings::storeRingsInfo(mol, res);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::fastFindRings
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
    // BEGIN RDKIT CPP FUNCTION MolOps::_DFS
    // RDKit❗✔️: void _DFS(const ROMol &mol, const Atom *atom, INT_VECT &atomColors,
    // RDKit❗✔️:           std::vector<const Atom *> &traversalOrder, VECT_INT_VECT &res,
    // RDKit❗✔️:           const Atom *fromAtom = nullptr) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(atomColors[atom->getIdx()] == 0, "bad color");
    // RDKit✔️✔️:   atomColors[atom->getIdx()] = 1;
    // RDKit✔️✔️:   traversalOrder.push_back(atom);
    atom_colors[atom] = 1;
    traversal_order.push(atom);
    // RDKit✔️✔️:   for (const auto nbr : mol.atomNeighbors(atom)) {
    for neighbor in context.neighbors(atom) {
        // RDKit✔️✔️:     unsigned int nbrIdx = nbr->getIdx();
        let neighbor_idx = neighbor.atom_index;
        // RDKit✔️✔️:     if (atomColors[nbrIdx] == 0) {
        if atom_colors[neighbor_idx] == 0 {
            // RDKit✔️✔️:       if (nbr->getDegree() < 2) {
            // RDKit✔️✔️:         atomColors[nbr->getIdx()] = 2;
            // RDKit✔️✔️:       } else {
            if context.neighbors(neighbor_idx).len() < 2 {
                atom_colors[neighbor_idx] = 2;
            } else {
                // RDKit✔️✔️:         _DFS(mol, nbr, atomColors, traversalOrder, res, atom);
                dfs_fast_find_rings(
                    context,
                    neighbor_idx,
                    atom_colors,
                    traversal_order,
                    result,
                    Some(atom),
                );
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     } else if (atomColors[nbrIdx] == 1) {
        } else if atom_colors[neighbor_idx] == 1
            && from_atom.is_some()
            && Some(neighbor_idx) != from_atom
        {
            // RDKit✔️✔️:       if (fromAtom && nbrIdx != fromAtom->getIdx()) {
            // RDKit✔️✔️:         INT_VECT cycle;
            // RDKit✔️✔️:         auto lastElem =
            // RDKit✔️✔️:             std::find(traversalOrder.rbegin(), traversalOrder.rend(), atom);
            // RDKit✔️✔️:         for (auto rIt = lastElem;
            // RDKit✔️✔️:              rIt != traversalOrder.rend() && (*rIt)->getIdx() != nbrIdx;
            // RDKit✔️✔️:              ++rIt) {
            // RDKit✔️✔️:           cycle.push_back((*rIt)->getIdx());
            // RDKit✔️✔️:         }
            let mut cycle = Vec::new();
            for &path_atom in traversal_order.iter().rev() {
                cycle.push(path_atom);
                if path_atom == neighbor_idx {
                    break;
                }
            }
            // RDKit✔️✔️:         cycle.push_back(nbrIdx);
            // RDKit✔️✔️:         res.push_back(cycle);
            result.push(cycle);
            // RDKit✔️✔️:       }
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   atomColors[atom->getIdx()] = 2;
    // RDKit✔️✔️:   traversalOrder.pop_back();
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::_DFS
    atom_colors[atom] = 2;
    traversal_order.pop();
}

#[cfg(test)]
mod tests {
    use crate::{
        AtomId, BondId, Molecule, RingFindType, RingInfo, find_ring_families, symmetrize_sssr,
    };

    #[test]
    fn symmetrize_sssr_finds_simple_cycle_and_membership() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CCCCC1", false).unwrap();
        let rings = symmetrize_sssr(&molecule).unwrap();

        assert_eq!(rings.find_type(), RingFindType::SymmSssr);
        assert_eq!(rings.num_rings(), 1);
        assert_eq!(rings.atom_rings()[0].len(), 6);
        assert_eq!(rings.bond_rings()[0].len(), 6);
        assert_eq!(rings.num_atom_rings(AtomId::new(0)), 1);
        assert_eq!(rings.num_bond_rings(BondId::new(0)), 1);
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(0), 6));
        assert!(rings.is_bond_in_ring_of_size(BondId::new(0), 6));
    }

    #[test]
    fn symmetrize_sssr_handles_acyclic_molecule() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let rings = symmetrize_sssr(&molecule).unwrap();

        assert_eq!(rings.num_rings(), 0);
        assert_eq!(rings.min_atom_ring_size(AtomId::new(0)), 0);
        assert_eq!(rings.min_bond_ring_size(BondId::new(0)), 0);
    }

    #[test]
    fn with_assigned_rings_updates_derived_cache_through_registered_operation() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
        let original = molecule.clone();
        let result = molecule.with_assigned_rings().unwrap();

        assert_eq!(molecule, original);
        assert!(result.derived_cache().adjacency.is_some());
        let rings = result
            .derived_cache()
            .rings
            .as_ref()
            .expect("assigned rings operation should materialize ring cache");
        assert_eq!(rings.num_rings(), 1);
        assert!(rings.is_symm_sssr());
        assert!(rings.is_atom_in_ring_of_size(AtomId::new(0), 3));
        assert!(rings.is_bond_in_ring_of_size(BondId::new(0), 3));
    }

    #[test]
    fn find_ring_families_materializes_urf_data() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
        let rings = find_ring_families(&molecule, false, false).unwrap();

        assert_eq!(rings.num_ring_families(), 1);
        assert_eq!(rings.num_relevant_cycles().unwrap(), 1);
        assert_eq!(
            rings.atom_ring_families()[0],
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
        );
        assert_eq!(
            rings.bond_ring_families()[0],
            vec![BondId::new(0), BondId::new(1), BondId::new(2)]
        );
    }

    #[test]
    fn with_assigned_ring_families_updates_separate_urf_cache_through_registered_operation() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
        let original = molecule.clone();

        let result = molecule.with_assigned_ring_families().unwrap();

        assert_eq!(molecule, original);
        assert!(result.derived_cache().adjacency.is_some());
        assert!(result.derived_cache().rings.is_none());
        let ring_families = result
            .derived_cache()
            .ring_families
            .as_ref()
            .expect("assigned ring families operation should materialize URF cache");
        assert_eq!(ring_families.num_ring_families(), 1);
        assert_eq!(ring_families.num_relevant_cycles().unwrap(), 1);
        assert_eq!(
            ring_families.atom_ring_families()[0],
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
        );
    }

    #[test]
    fn with_assigned_rings_does_not_materialize_urf_cache() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();

        let result = molecule.with_assigned_rings().unwrap();

        assert!(result.derived_cache().rings.is_some());
        assert!(result.derived_cache().ring_families.is_none());
    }

    #[test]
    fn ring_info_reports_fused_ring_neighbors_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CCC2CCCCC2C1", false).unwrap();
        let mut rings = symmetrize_sssr(&molecule).unwrap();

        assert_eq!(rings.num_rings(), 2);
        assert!(rings.are_rings_fused(0, 1).unwrap());
        assert!(rings.is_ring_fused(0).unwrap());
        assert!(rings.is_ring_fused(1).unwrap());
        assert_eq!(rings.num_fused_ring_neighbors(0).unwrap(), 1);
        assert_eq!(rings.num_fused_ring_neighbors(1).unwrap(), 1);
        assert_eq!(rings.fused_ring_neighbors(0).unwrap(), vec![1]);
        assert_eq!(rings.fused_ring_neighbors(1).unwrap(), vec![0]);
        assert_eq!(rings.num_fused_bonds(0).unwrap(), 1);
        assert_eq!(rings.num_fused_bonds(1).unwrap(), 1);
    }

    #[test]
    fn ring_info_ring_family_storage_matches_rdkit_container_semantics() {
        let mut rings = RingInfo::new(RingFindType::OtherOrUnknown, 3, 3);

        assert_eq!(rings.add_ring_family(&[0, 1, 2], &[0, 1, 2]).unwrap(), 1);
        assert_eq!(rings.num_ring_families(), 1);
        assert_eq!(
            rings.atom_ring_families(),
            &[vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]]
        );
        assert_eq!(
            rings.bond_ring_families(),
            &[vec![BondId::new(0), BondId::new(1), BondId::new(2)]]
        );
        assert!(rings.num_relevant_cycles().is_err());
    }
}
