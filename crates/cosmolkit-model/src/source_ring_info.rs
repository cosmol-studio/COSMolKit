//! Explicit detached source RingInfo state. Ring algorithms remain in CORE.

use crate::{AtomId, BondId};

/// Lossless transport of the currently modeled source ring-cache fields.
/// This is actual state supplied by its owner, never a ring-finding request.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceRingInfo {
    pub initialized: bool,
    /// Pinned FIND_RING_TYPE ordinal: FAST=0, SSSR=1, SYMM_SSSR=2, OTHER=3.
    pub find_type: u32,
    pub atom_members: Vec<Vec<usize>>,
    pub bond_members: Vec<Vec<usize>>,
    pub atom_rings: Vec<Vec<AtomId>>,
    pub bond_rings: Vec<Vec<BondId>>,
    pub atom_ring_families: Vec<Vec<AtomId>>,
    pub bond_ring_families: Vec<Vec<BondId>>,
    pub relevant_cycle_count: Option<usize>,
    pub fused_rings: Vec<Vec<bool>>,
    pub num_fused_bonds: Vec<usize>,
}

impl Default for SourceRingInfo {
    fn default() -> Self {
        // RDKit❗✔️:   RingInfo() {}
        // RDKit❗✔️:   bool df_init{false};
        // RDKit❗✔️:   FIND_RING_TYPE df_find_type_type{FIND_RING_TYPE_OTHER_OR_UNKNOWN};
        // RDKit❗✔️:   DataType d_atomMembers, d_bondMembers;
        // RDKit❗✔️:   VECT_INT_VECT d_atomRings, d_bondRings;
        // RDKit❗✔️:   VECT_INT_VECT d_atomRingFamilies, d_bondRingFamilies;
        // RDKit❗✔️:   std::vector<boost::dynamic_bitset<>> d_fusedRings;
        // RDKit❗✔️:   std::vector<unsigned int> d_numFusedBonds;
        // RDKit❗✔️: typedef enum {
        // RDKit❗✔️:   FIND_RING_TYPE_FAST,
        // RDKit❗✔️:   FIND_RING_TYPE_SSSR,
        // RDKit❗✔️:   FIND_RING_TYPE_SYMM_SSSR,
        // RDKit❗✔️:   FIND_RING_TYPE_OTHER_OR_UNKNOWN
        // RDKit❗✔️: } FIND_RING_TYPE;
        // A newly constructed source ROMol owns this uninitialized cache;
        // this constructor supplies no cache for a copied/extracted molecule.
        // The existing CORE URF object remains independently unmodeled here;
        // its represented optional relevant-cycle scalar is transported only.
        // O(1), allocation-free empty vectors, like Native default members.
        Self {
            initialized: false,
            find_type: 3,
            atom_members: Vec::new(),
            bond_members: Vec::new(),
            atom_rings: Vec::new(),
            bond_rings: Vec::new(),
            atom_ring_families: Vec::new(),
            bond_ring_families: Vec::new(),
            relevant_cycle_count: None,
            fused_rings: Vec::new(),
            num_fused_bonds: Vec::new(),
        }
    }
}
