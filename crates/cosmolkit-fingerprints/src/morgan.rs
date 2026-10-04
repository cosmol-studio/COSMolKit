//! Detached Morgan generator arguments and default-generator selection.
//!
//! Source attribution: pinned RDKit `GraphMol/Fingerprints/MorganGenerator.cpp/.h`,
//! `FingerprintGenerator.cpp/.h`, and `MorganFingerprints.cpp`, plus the
//! copied `Code/RDGeneral` property/hash helpers; Boost dynamic_bitset 1.85.0
//! is pinned separately. Exact paths, notices and licenses are in
//! `THIRD_PARTY_NOTICES.md`.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};
use std::hash::{Hash, Hasher};

use crate::additional_output::AdditionalOutput;
use crate::generator::{
    FingerprintArguments, FingerprintFuncArguments, accumulate_morgan_sparse_counts_into,
    get_count_fingerprint, get_count_fingerprint_with_atom_invariants, get_fingerprint,
    get_fingerprint_with_atom_invariants, get_sparse_count_fingerprint_with_atom_invariants,
    get_sparse_fingerprint_with_atom_invariants, with_morgan_environment_inputs_and_output,
};
use crate::hash::{hash_combine, hash_value_i32, hash_value_u32};
use crate::invariants::{
    MorganAtomInvGenerator, MorganBondInvGenerator, MorganFeatureAtomInvGenerator,
};
use crate::prepared::MorganPreparedInput;
use crate::sparse_counts::SparseCountFingerprint32;
use crate::{Fingerprint, FingerprintError, MorganError};
use cosmolkit_core::{
    DenseMatrix, MatrixError, RingInfo, TopologicalDistanceMatrixParams, ValenceAssignment,
    topological_distance_matrix,
};
use cosmolkit_model::{BondOrder, ChiralTag, MoleculeProperties, PropertyValue, TopologyBlock};
use cosmolkit_search::{SearchTarget, build_prepared_query_match_context};

/// The two output widths explicitly instantiated by the pinned generator.
pub(super) trait MorganOutput: Copy {
    fn from_environment_code(code: u32) -> Self;
    fn into_u32(self) -> u32;
    fn into_u64(self) -> u64;
}

impl MorganOutput for u32 {
    fn from_environment_code(code: u32) -> Self {
        code
    }

    fn into_u32(self) -> u32 {
        self
    }

    fn into_u64(self) -> u64 {
        u64::from(self)
    }
}

impl MorganOutput for u64 {
    fn from_environment_code(code: u32) -> Self {
        u64::from(code)
    }

    fn into_u32(self) -> u32 {
        self as u32
    }

    fn into_u64(self) -> u64 {
        self
    }
}

/// Packed Boost bitset storage used for Morgan neighborhoods and atom state.
#[derive(Debug, Clone)]
struct MorganBondEnvironment {
    bit_count: usize,
    blocks: Vec<u64>,
}

impl PartialEq for MorganBondEnvironment {
    fn eq(&self, other: &Self) -> bool {
        // BEGIN BOOST CPP FUNCTION operator==(dynamic_bitset const&, dynamic_bitset const&)
        // Boost❗✔️: template <typename Block, typename Allocator>
        // Boost❗✔️: bool operator==(const dynamic_bitset<Block, Allocator>& a,
        // Boost❗✔️:                 const dynamic_bitset<Block, Allocator>& b)
        // Boost❗✔️: {
        // Boost❗✔️:     return (a.m_num_bits == b.m_num_bits)
        // Boost❗✔️:            && (a.m_bits == b.m_bits);
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION operator==(dynamic_bitset const&, dynamic_bitset const&)
        self.bit_count == other.bit_count && self.blocks == other.blocks
    }
}

impl Eq for MorganBondEnvironment {}

impl Hash for MorganBondEnvironment {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // BEGIN BOOST CPP FUNCTION hash_value(dynamic_bitset const&)
        // Boost❗✔️: template <typename Block, typename Allocator>
        // Boost❗✔️: inline std::size_t hash_value(const dynamic_bitset<Block, Allocator>& a)
        // Boost❗✔️: {
        // Boost❗✔️:     std::size_t res = hash_value(a.m_num_bits);
        // Boost❗✔️:     boost::hash_combine(res, a.m_bits);
        // Boost❗✔️:     return res;
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION hash_value(dynamic_bitset const&)
        // Behavior review: the set observes only key equality/membership; this
        // hash consumes the same logical width and packed words. Its concrete
        // bucket hash is not emitted or iterated by the source algorithm.
        // Complexity review: both source and Rust visit O(B/64) packed words.
        self.bit_count.hash(state);
        self.blocks.hash(state);
    }
}

impl Ord for MorganBondEnvironment {
    fn cmp(&self, other: &Self) -> Ordering {
        // BEGIN BOOST CPP FUNCTION operator<(dynamic_bitset const&, dynamic_bitset const&)
        // Boost❗✔️: template <typename Block, typename Allocator>
        // Boost❗✔️: bool operator<(const dynamic_bitset<Block, Allocator>& a,
        // Boost❗✔️:                const dynamic_bitset<Block, Allocator>& b)
        // Boost❗✔️: {
        // Boost❗✔️:     typedef BOOST_DEDUCED_TYPENAME dynamic_bitset<Block, Allocator>::size_type size_type;
        // Boost❗✔️:     size_type asize(a.size());
        // Boost❗✔️:     size_type bsize(b.size());
        // Boost❗✔️:     if (!bsize)
        // Boost❗✔️:         {
        // Boost❗✔️:         return false;
        // Boost❗✔️:         }
        // Boost❗✔️:     else if (!asize)
        // Boost❗✔️:         {
        // Boost❗✔️:         return true;
        // Boost❗✔️:         }
        // Boost❗✔️:     else if (asize == bsize)
        // Boost❗✔️:         {
        // Boost❗✔️:         for (size_type ii = a.num_blocks(); ii > 0; --ii)
        // Boost❗✔️:             {
        // Boost❗✔️:             size_type i = ii-1;
        // Boost❗✔️:             if (a.m_bits[i] < b.m_bits[i])
        // Boost❗✔️:                 return true;
        // Boost❗✔️:             else if (a.m_bits[i] > b.m_bits[i])
        // Boost❗✔️:                 return false;
        // Boost❗✔️:             }
        // Boost❗✔️:         return false;
        // Boost❗✔️:         }
        // Boost❗✔️:     else
        // Boost❗✔️:         {
        // Boost❗✔️:         size_type leqsize(std::min BOOST_PREVENT_MACRO_SUBSTITUTION(asize,bsize));
        // Boost❗✔️:         for (size_type ii = 0; ii < leqsize; ++ii,--asize,--bsize)
        // Boost❗✔️:             {
        // Boost❗✔️:             size_type i = asize-1;
        // Boost❗✔️:             size_type j = bsize-1;
        // Boost❗✔️:             if (a[i] < b[j])
        // Boost❗✔️:                 return true;
        // Boost❗✔️:             else if (a[i] > b[j])
        // Boost❗✔️:                 return false;
        // Boost❗✔️:             }
        // Boost❗✔️:         return (a.size() < b.size());
        // Boost❗✔️:         }
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION operator<(dynamic_bitset const&, dynamic_bitset const&)
        // Behavior review: source-width masks use unsigned packed-word comparison
        // from the most-significant word down. The separate unequal-width path
        // compares each logical top bit, then compares bit counts.
        // Complexity review: equal-width comparison is O(B/64); the source's
        // unequal-width path and this translation are O(min(B1, B2)).
        if other.bit_count == 0 {
            return if self.bit_count == 0 {
                Ordering::Equal
            } else {
                Ordering::Greater
            };
        }
        if self.bit_count == 0 {
            return Ordering::Less;
        }

        if self.bit_count == other.bit_count {
            for (&left_word, &right_word) in self.blocks.iter().rev().zip(other.blocks.iter().rev())
            {
                match left_word.cmp(&right_word) {
                    Ordering::Equal => {}
                    ordering => return ordering,
                }
            }
            return Ordering::Equal;
        }

        for offset in 0..self.bit_count.min(other.bit_count) {
            let left_index = self.bit_count - offset - 1;
            let right_index = other.bit_count - offset - 1;
            let left_bit = self.blocks[left_index / u64::BITS as usize]
                & (1u64 << (left_index % u64::BITS as usize))
                != 0;
            let right_bit = other.blocks[right_index / u64::BITS as usize]
                & (1u64 << (right_index % u64::BITS as usize))
                != 0;
            match left_bit.cmp(&right_bit) {
                Ordering::Equal => {}
                ordering => return ordering,
            }
        }
        self.bit_count.cmp(&other.bit_count)
    }
}

impl PartialOrd for MorganBondEnvironment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl MorganBondEnvironment {
    /// Create a zero-filled dynamic-bitset equivalent for the pinned LP64 source.
    fn new(bit_count: usize) -> Self {
        // BEGIN BOOST CPP FUNCTION dynamic_bitset(size_type, unsigned long)
        // Boost❗✔️: explicit dynamic_bitset(size_type num_bits, unsigned long value = 0,
        // Boost❗✔️:                    const Allocator& alloc = Allocator());
        // Boost❗✔️: m_bits.resize(calc_num_blocks(num_bits));
        // Boost❗✔️: m_num_bits = num_bits;
        // END BOOST CPP FUNCTION dynamic_bitset(size_type, unsigned long)
        // The fixed RDKit reference ABI is LP64: dynamic_bitset's default
        // unsigned-long block is 64 bits. Vec resize value-initializes the
        // same zero words; the bit count remains part of value identity.
        const BITS_PER_BLOCK: usize = u64::BITS as usize;
        let block_count = bit_count.div_ceil(BITS_PER_BLOCK);
        Self {
            bit_count,
            blocks: vec![0; block_count],
        }
    }

    fn set_all(&mut self) {
        // BEGIN BOOST CPP FUNCTION dynamic_bitset::set()
        // Boost❗✔️: dynamic_bitset<Block, Allocator>&
        // Boost❗✔️: dynamic_bitset<Block, Allocator>::set()
        // Boost❗✔️: {
        // Boost❗✔️:   std::fill(m_bits.begin(), m_bits.end(),
        // Boost❗✔️:              detail::dynamic_bitset_impl::max_limit<Block>::value);
        // Boost❗✔️:   m_zero_unused_bits();
        // Boost❗✔️:   return *this;
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION dynamic_bitset::set()
        // Behavior: fill all logical bits and clear unused high bits, including
        // the empty-bitset case. Complexity: both paths touch O(B/64) words.
        self.blocks.fill(u64::MAX);
        let used_bits = self.bit_count % u64::BITS as usize;
        if used_bits != 0 {
            let last = self.blocks.len() - 1;
            self.blocks[last] &= (1_u64 << used_bits) - 1;
        }
    }

    fn set(&mut self, bit: usize) {
        // BEGIN BOOST CPP FUNCTION dynamic_bitset::set(size_type, bool)
        // Boost❗✔️: dynamic_bitset<Block, Allocator>::set(size_type pos, bool val)
        // Boost❗✔️: {
        // Boost❗✔️:     assert(pos < m_num_bits);
        // Boost❗✔️:     if (val)
        // Boost❗✔️:         m_bits[block_index(pos)] |= bit_mask(pos);
        // Boost❗✔️:     else
        // Boost❗✔️:         reset(pos);
        // Boost❗✔️:     return *this;
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION dynamic_bitset::set(size_type, bool)
        // Behavior: the caller sets only a valid source bond or atom ID. The
        // assertion retains Boost's invalid-index boundary. The pinned
        // default block uses quotient/remainder indexing.
        assert!(
            bit < self.bit_count,
            "source dynamic_bitset index is in range"
        );
        const BITS_PER_BLOCK: usize = u64::BITS as usize;
        let block_index = bit / BITS_PER_BLOCK;
        let bit_index = bit % BITS_PER_BLOCK;
        self.blocks[block_index] |= 1u64 << bit_index;
    }

    fn contains(&self, bit: usize) -> bool {
        // BEGIN BOOST CPP FUNCTION dynamic_bitset::operator[] const
        // Boost❗✔️: bool operator[](size_type pos) const { return test(pos); }
        // Boost❗✔️: bool dynamic_bitset<Block, Allocator>::test(size_type pos) const
        // Boost❗✔️: {
        // Boost❗✔️:     assert(pos < m_num_bits);
        // Boost❗✔️:     return m_unchecked_test(pos);
        // Boost❗✔️: }
        // Boost❗✔️: bool dynamic_bitset<Block, Allocator>::m_unchecked_test(size_type pos) const
        // Boost❗✔️: {
        // Boost❗✔️:     return (m_bits[block_index(pos)] & bit_mask(pos)) != 0;
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION dynamic_bitset::operator[] const
        assert!(
            bit < self.bit_count,
            "source dynamic_bitset index is in range"
        );
        const BITS_PER_BLOCK: usize = u64::BITS as usize;
        let block_index = bit / BITS_PER_BLOCK;
        let bit_index = bit % BITS_PER_BLOCK;
        self.blocks[block_index] & (1u64 << bit_index) != 0
    }

    fn union_with(&mut self, other: &Self) {
        // BEGIN BOOST CPP FUNCTION dynamic_bitset::operator|=
        // Boost❗✔️: dynamic_bitset<Block, Allocator>::operator|=(const dynamic_bitset& rhs)
        // Boost❗✔️: {
        // Boost❗✔️:     assert(size() == rhs.size());
        // Boost❗✔️:     for (size_type i = 0; i < num_blocks(); ++i)
        // Boost❗✔️:         m_bits[i] |= rhs.m_bits[i];
        // Boost❗✔️:     return *this;
        // Boost❗✔️: }
        // END BOOST CPP FUNCTION dynamic_bitset::operator|=
        assert_eq!(self.bit_count, other.bit_count, "source bitset sizes match");
        for (block, other_block) in self.blocks.iter_mut().zip(&other.blocks) {
            *block |= *other_block;
        }
    }
}

/// The source allocates this packed atom-row mask once and reuses it by radius.
type MorganChiralAtoms = MorganBondEnvironment;

/// One lazily prepared distance matrix shared by all environments in a call.
#[derive(Default)]
pub(super) struct MorganDistanceMatrixCache<'a> {
    topology: Option<&'a TopologyBlock>,
    matrix: Option<DenseMatrix>,
}

impl<'a> MorganDistanceMatrixCache<'a> {
    fn get_or_prepare(&mut self, topology: &'a TopologyBlock) -> Result<&DenseMatrix, MatrixError> {
        if !self
            .topology
            .is_some_and(|cached| std::ptr::eq(cached, topology))
        {
            self.topology = None;
            self.matrix = None;
        }
        if self.matrix.is_none() {
            let matrix =
                topological_distance_matrix(topology, &TopologicalDistanceMatrixParams::default())?;
            self.topology = Some(topology);
            self.matrix = Some(matrix);
        }
        Ok(self
            .matrix
            .as_ref()
            .expect("the matrix is assigned when the cache is empty"))
    }
}

/// Stored environment identity used by the scheduled M02 and M03 owners.
pub(super) struct MorganAtomEnvironment<'a, OutputType: MorganOutput> {
    d_code: OutputType,
    d_atom_id: u32,
    d_layer: u32,
    topology: &'a TopologyBlock,
}

impl<'a, OutputType: MorganOutput> MorganAtomEnvironment<'a, OutputType> {
    pub(super) fn new(code: u32, atom_id: u32, layer: u32, topology: &'a TopologyBlock) -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganAtomEnv::MorganAtomEnv
        // RDKit❗✔️: MorganAtomEnv(const std::uint32_t code, const unsigned int atomId,
        // RDKit❗✔️:               const unsigned int layer, const ROMol *mol)
        // RDKit❗✔️:     : d_code(code), d_atomId(atomId), d_layer(layer), d_mol(mol) {}
        // END RDKIT CPP FUNCTION MorganAtomEnv::MorganAtomEnv
        // Behavior: store the source code, atom center, radius layer, and the
        // same detached topology supplied to the environment generator.
        Self {
            d_code: OutputType::from_environment_code(code),
            d_atom_id: atom_id,
            d_layer: layer,
            topology,
        }
    }

    pub(super) fn get_bit_id(
        &self,
        _arguments: Option<&FingerprintArguments>,
        _atom_invariants: Option<&[u32]>,
        _bond_invariants: Option<&[u32]>,
        _additional_output: Option<&mut AdditionalOutput>,
        _hash_results: bool,
        _fp_size: u64,
    ) -> OutputType {
        // BEGIN RDKIT CPP FUNCTION MorganAtomEnv<OutputType>::getBitId
        // RDKit❗✔️: OutputType MorganAtomEnv<OutputType>::getBitId(
        // RDKit❗✔️:     FingerprintArguments *,              // arguments
        // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // atomInvariants
        // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit❗✔️:     AdditionalOutput *,                  // additional Output
        // RDKit❗✔️:     const bool,                          // hashResults
        // RDKit❗✔️:     const std::uint64_t                  // fpSize
        // RDKit❗✔️: ) const {
        // RDKit❗✔️:   return d_code;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION MorganAtomEnv<OutputType>::getBitId
        // Behavior: return the stored output-width identity; every source
        // argument is deliberately ignored, including the hash and size flags.
        // Complexity: one copied scalar field, no allocation, branch, or scan.
        self.d_code
    }

    pub(super) fn update_additional_output(
        &self,
        additional_output: &mut AdditionalOutput,
        bit_id: u64,
        distance_matrix_cache: &mut MorganDistanceMatrixCache<'a>,
    ) -> Result<(), MatrixError> {
        // BEGIN RDKIT CPP FUNCTION MorganAtomEnv<OutputType>::updateAdditionalOutput
        // RDKit❗✔️: void MorganAtomEnv<OutputType>::updateAdditionalOutput(
        // RDKit❗✔️:     AdditionalOutput *additionalOutput, size_t bitId) const {
        // RDKit❗✔️:   PRECONDITION(additionalOutput, "bad output pointer");
        // RDKit❗✔️:   PRECONDITION(d_mol, "bad mol pointer");
        // RDKit❗✔️:   if (additionalOutput->bitInfoMap) {
        // RDKit❗✔️:     (*additionalOutput->bitInfoMap)[bitId].emplace_back(d_atomId, d_layer);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (additionalOutput->atomCounts) {
        // RDKit❗✔️:     (*additionalOutput->atomCounts)[d_atomId]++;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (additionalOutput->atomToBits) {
        // RDKit❗✔️:     (*additionalOutput->atomToBits)[d_atomId].push_back(bitId);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (additionalOutput->atomsPerBit) {
        // RDKit❗✔️:     std::vector<int> atomsInvolved;
        // RDKit❗✔️:     atomsInvolved.push_back(d_atomId);
        // RDKit❗✔️:     if (d_layer > 0) {
        // RDKit❗✔️:       const auto dm = MolOps::getDistanceMat(*d_mol);
        // RDKit❗✔️:       for (unsigned int i = 0; i < d_mol->getNumAtoms(); ++i) {
        // RDKit❗✔️:         if (static_cast<unsigned int>(dm[d_atomId * d_mol->getNumAtoms() + i] +
        // RDKit❗✔️:                                       .1) <= d_layer &&
        // RDKit❗✔️:             i != d_atomId) {
        // RDKit❗✔️:           atomsInvolved.push_back(i);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     (*additionalOutput->atomsPerBit)[bitId].push_back(std::move(atomsInvolved));
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION MorganAtomEnv<OutputType>::updateAdditionalOutput
        // Behavior review: borrowed mutable output and topology references encode
        // the source non-null preconditions; all five field-presence states are
        // independent, and source mutation order and repeated appends are kept.
        // The source's unsigned atomCounts increment wraps at u32 width here.
        // Complexity review: count/index writes are O(1), ordered map insertion
        // is O(log B), vector appends are amortized O(1), and positive-layer
        // atomsPerBit scans O(A) after one shared O(A^3)/O(A^2) matrix build.
        if let Some(bit_info_map) = &mut additional_output.bit_info_map {
            bit_info_map
                .entry(bit_id)
                .or_default()
                .push((self.d_atom_id, self.d_layer));
        }
        if let Some(atom_counts) = &mut additional_output.atom_counts {
            let count = &mut atom_counts[self.d_atom_id as usize];
            *count = count.wrapping_add(1);
        }
        if let Some(atom_to_bits) = &mut additional_output.atom_to_bits {
            atom_to_bits[self.d_atom_id as usize].push(bit_id);
        }
        if let Some(atoms_per_bit) = &mut additional_output.atoms_per_bit {
            let center = self.d_atom_id as usize;
            let mut atoms_involved = Vec::new();
            atoms_involved.push(self.d_atom_id as i32);
            if self.d_layer > 0 {
                let distance_matrix = distance_matrix_cache.get_or_prepare(self.topology)?;
                for atom_index in 0..self.topology.atoms.len() {
                    let distance = distance_matrix
                        .get(center, atom_index)
                        .expect("the environment center and atom rows match the matrix");
                    if (distance + 0.1) as u32 <= self.d_layer && atom_index != center {
                        atoms_involved.push(atom_index as i32);
                    }
                }
            }
            atoms_per_bit
                .entry(bit_id)
                .or_default()
                .push(atoms_involved);
        }

        Ok(())
    }
}

/// The source AccumTuple order: exact neighborhood mask, uint32 code, atom ID.
type MorganLayerCandidate = (MorganBondEnvironment, u32, u32);

/// Sort and emit one completed source layer, preserving the source's shared
/// neighborhood identity set and sticky duplicate-center state.
fn collect_morgan_layer<'a, OutputType: MorganOutput>(
    candidates: &mut [MorganLayerCandidate],
    include_redundant_environments: bool,
    only_nonzero_invariants: bool,
    atom_invariants: &[u32],
    include_atoms: &MorganBondEnvironment,
    neighborhoods: &mut HashSet<MorganBondEnvironment>,
    dead_atoms: &mut MorganBondEnvironment,
    layer: u32,
    topology: &'a TopologyBlock,
    result: &mut Vec<MorganAtomEnvironment<'a, OutputType>>,
) {
    // BEGIN RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments layer sort and dedup
    // RDKit❗✔️: typedef std::tuple<boost::dynamic_bitset<>, uint32_t, unsigned int> AccumTuple;
    // RDKit❗✔️:     std::sort(allNeighborhoodsThisRound.begin(),
    // RDKit❗✔️:               allNeighborhoodsThisRound.end());
    // RDKit❗✔️:     for (std::vector<AccumTuple>::const_iterator iter =
    // RDKit❗✔️:              allNeighborhoodsThisRound.begin();
    // RDKit❗✔️:          iter != allNeighborhoodsThisRound.end(); ++iter) {
    // RDKit❗✔️:       if (morganArguments->df_includeRedundantEnvironments ||
    // RDKit❗✔️:           neighborhoods.count(std::get<0>(*iter)) == 0) {
    // RDKit❗✔️:         if (!morganArguments->df_onlyNonzeroInvariants ||
    // RDKit❗✔️:             (*atomInvariants)[std::get<2>(*iter)]) {
    // RDKit❗✔️:           if (includeAtoms[std::get<2>(*iter)]) {
    // RDKit❗✔️:             result.push_back(new MorganAtomEnv<OutputType>(
    // RDKit❗✔️:                 std::get<1>(*iter), std::get<2>(*iter), layer + 1, &mol));
    // RDKit❗✔️:             neighborhoods.insert(std::get<0>(*iter));
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         deadAtoms[std::get<2>(*iter)] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // END RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments layer sort and dedup
    // Behavior review: the source tuple's three fields are compared in order;
    // membership is mask-only and spans layers. Redundant=true short-circuits
    // membership but still inserts emitted masks. Otherwise a duplicate makes
    // its center permanently dead; an unseen but zero/excluded row does not
    // reserve its mask. Filters use original atom invariants and selected IDs.
    // Complexity review: sorting is O(A log A) with O(B/64) equal-width mask
    // comparisons, set membership is expected O(1) with O(B/64) key work, and
    // each emitted mask is copied once into the persistent set.
    candidates.sort_unstable();
    assert_eq!(include_atoms.bit_count, topology.atoms.len());
    assert_eq!(dead_atoms.bit_count, topology.atoms.len());
    assert!(atom_invariants.len() >= topology.atoms.len());
    for (neighborhood, code, atom_id) in candidates {
        assert_eq!(neighborhood.bit_count, topology.bonds.len());
        let atom_index =
            usize::try_from(*atom_id).expect("source atom index fits the platform size type");
        assert!(atom_index < topology.atoms.len());

        if include_redundant_environments || !neighborhoods.contains(neighborhood) {
            if !only_nonzero_invariants || atom_invariants[atom_index] != 0 {
                if include_atoms.contains(atom_index) {
                    result.push(MorganAtomEnvironment::new(
                        *code,
                        *atom_id,
                        layer + 1,
                        topology,
                    ));
                    neighborhoods.insert(neighborhood.clone());
                }
            }
        } else {
            dead_atoms.set(atom_index);
        }
    }
}

/// Build source includeAtoms once so radius zero and every later layer share it.
fn selected_atom_mask(atom_count: usize, from_atoms: Option<&[u32]>) -> MorganBondEnvironment {
    // BEGIN RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments includeAtoms
    // RDKit❗✔️:   boost::dynamic_bitset<> includeAtoms(nAtoms);
    // RDKit❗✔️:   if (fromAtoms) {
    // RDKit❗✔️:     for (auto idx : *fromAtoms) {
    // RDKit❗✔️:       includeAtoms.set(idx, 1);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     includeAtoms.set();
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments includeAtoms
    // Behavior: duplicates remain idempotent; invalid IDs retain the pinned
    // bitset assertion boundary; absent selects every logical atom bit.
    // Complexity: O(F) indexed writes or O(A/64) packed-word fill, then the
    // same O(1) membership read per atom as the source bitset.
    let mut include_atoms = MorganBondEnvironment::new(atom_count);
    if let Some(from_atoms) = from_atoms {
        for &source_index in from_atoms {
            let atom_index = usize::try_from(source_index)
                .expect("source atom index fits the platform size type");
            include_atoms.set(atom_index);
        }
    } else {
        include_atoms.set_all();
    }
    include_atoms
}

/// Append source radius-zero environments in ascending atom-index order.
fn radius_zero_environments<'a, OutputType: MorganOutput>(
    topology: &'a TopologyBlock,
    current_invariants: &[u32],
    include_atoms: &MorganBondEnvironment,
    _ignore_atoms: Option<&[u32]>,
    only_nonzero_invariants: bool,
    environments: &mut Vec<MorganAtomEnvironment<'a, OutputType>>,
) {
    // BEGIN RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments radius-zero selection
    // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // ignoreAtoms
    // RDKit❗✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:     if (includeAtoms[i]) {
    // RDKit❗✔️:       if (!morganArguments->df_onlyNonzeroInvariants ||
    // RDKit❗✔️:           currentInvariants[i]) {
    // RDKit❗✔️:         result.push_back(
    // RDKit❗✔️:             new MorganAtomEnv<OutputType>(currentInvariants[i], i, 0, &mol));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments radius-zero selection
    // Behavior: source atom order, selected bits and original nonzero values
    // are retained; ignoreAtoms is deliberately not inspected. Complexity:
    // O(A) indexed scan with no second result allocation or environment heap
    // object per row; the caller already reserved source result capacity.
    let atom_count = topology.atoms.len();
    assert!(current_invariants.len() >= atom_count);
    assert_eq!(include_atoms.bit_count, atom_count);
    for atom_index in 0..atom_count {
        if include_atoms.contains(atom_index)
            && (!only_nonzero_invariants || current_invariants[atom_index] != 0)
        {
            environments.push(MorganAtomEnvironment::new(
                current_invariants[atom_index],
                atom_index as u32,
                0,
                topology,
            ));
        }
    }
}

/// Enumerate the source Morgan environments for one prepared topology.
///
/// Atom and bond invariant slices are supplied by the caller; canonical
/// integration computes them on the original topology while passing the
/// source-selected prepared topology here.
pub(super) fn generate_morgan_environments<'a, OutputType: MorganOutput>(
    topology: &'a TopologyBlock,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    atom_invariants: &[u32],
    bond_invariants: &[u32],
) -> Result<Vec<MorganAtomEnvironment<'a, OutputType>>, MorganError> {
    // BEGIN RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments
    // RDKit❗✔️: template <typename OutputType>
    // RDKit❗✔️: std::vector<AtomEnvironment<OutputType> *>
    // RDKit❗✔️: MorganEnvGenerator<OutputType>::getEnvironments(
    // RDKit❗✔️:     const ROMol &mol, FingerprintArguments *arguments,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // ignoreAtoms
    // RDKit❗✔️:     const int,                           // confId
    // RDKit❗✔️:     const AdditionalOutput *,            // additionalOutput
    // RDKit❗✔️:     const std::vector<std::uint32_t> *atomInvariants,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *bondInvariants,
    // RDKit❗✔️:     const bool  // hashResults
    // RDKit❗✔️: ) const {
    // RDKit❗✔️:   PRECONDITION(atomInvariants && (atomInvariants->size() >= mol.getNumAtoms()),
    // RDKit❗✔️:                "bad atom invariants size");
    // RDKit❗✔️:   PRECONDITION(bondInvariants && (bondInvariants->size() >= mol.getNumBonds()),
    // RDKit❗✔️:                "bad bond invariants size");
    // RDKit❗✔️:   auto *morganArguments = dynamic_cast<MorganArguments *>(arguments);
    // RDKit❗✔️:   PRECONDITION(morganArguments, "bad arguments type");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗✔️:   const unsigned int maxNumResults = (morganArguments->d_radius + 1) * nAtoms;
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<AtomEnvironment<OutputType> *> result =
    // RDKit❗✔️:       std::vector<AtomEnvironment<OutputType> *>();
    // RDKit❗✔️:   result.reserve(maxNumResults);
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we are using chirality, we need to make sure the atoms have R/S labels
    // RDKit❌❌:   if (morganArguments->df_includeChirality &&
    // RDKit❌❌:       !Chirality::getUseLegacyStereoPerception() &&
    // RDKit❌❌:       !mol.hasProp(common_properties::_CIPComputed)) {
    // RDKit❌❌:     CIPLabeler::assignCIPLabels(const_cast<ROMol &>(mol));
    // RDKit❌❌:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<OutputType> currentInvariants(atomInvariants->size());
    // RDKit❗✔️:   std::copy(atomInvariants->begin(), atomInvariants->end(),
    // RDKit❗✔️:             currentInvariants.begin());
    // RDKit❗✔️:   // will hold bit ids calculated this round to be used as invariants next
    // RDKit❗✔️:   // round
    // RDKit❗✔️:   std::vector<OutputType> nextLayerInvariants(nAtoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   // will hold up to date invariants of neighboring atoms with bond
    // RDKit❗✔️:   // types, these invariants hold information from atoms around radius
    // RDKit❗✔️:   // as big as current layer around the current atom
    // RDKit❗✔️:   std::vector<std::pair<int32_t, uint32_t>> neighborhoodInvariants;
    // RDKit❗✔️:   // Max number of neighbors expected.
    // RDKit❗✔️:   neighborhoodInvariants.reserve(8);
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::dynamic_bitset<> includeAtoms(nAtoms);
    // RDKit❗✔️:   if (fromAtoms) {
    // RDKit❗✔️:     for (auto idx : *fromAtoms) {
    // RDKit❗✔️:       includeAtoms.set(idx, 1);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     includeAtoms.set();
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::dynamic_bitset<> chiralAtoms(nAtoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   // these are the neighborhoods that have already been added
    // RDKit❗✔️:   // to the fingerprint
    // RDKit❗✔️:   std::unordered_set<boost::dynamic_bitset<>> neighborhoods;
    // RDKit❗✔️:   neighborhoods.reserve(maxNumResults);
    // RDKit❗✔️:   // these are the environments around each atom:
    // RDKit❗✔️:   std::vector<boost::dynamic_bitset<>> atomNeighborhoods(
    // RDKit❗✔️:       nAtoms, boost::dynamic_bitset<>(mol.getNumBonds()));
    // RDKit❗✔️:   // holds atoms in the environment (neighborhood) for the current layer for
    // RDKit❗✔️:   // each atom, starts with the immediate neighbors of atoms and expands
    // RDKit❗✔️:   // with every iteration
    // RDKit❗✔️:   std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods =
    // RDKit❗✔️:       atomNeighborhoods;
    // RDKit❗✔️:   boost::dynamic_bitset<> deadAtoms(nAtoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   // if df_onlyNonzeroInvariants is set order the atoms to make sure atoms
    // RDKit❗✔️:   // with zero invariants are processed last so that in case of duplicate
    // RDKit❗✔️:   // environments atoms with non-zero invariants are used
    // RDKit❗✔️:   std::vector<unsigned int> atomOrder(nAtoms);
    // RDKit❗✔️:   if (morganArguments->df_onlyNonzeroInvariants) {
    // RDKit❗✔️:     std::vector<std::pair<int32_t, uint32_t>> ordering;
    // RDKit❗✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:       if (!currentInvariants[i]) {
    // RDKit❗✔️:         ordering.emplace_back(1, i);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         ordering.emplace_back(0, i);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     std::sort(ordering.begin(), ordering.end());
    // RDKit❗✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:       atomOrder[i] = ordering[i].second;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:       atomOrder[i] = i;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // add the round 0 invariants to the result
    // RDKit❗✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:     if (includeAtoms[i]) {
    // RDKit❗✔️:       if (!morganArguments->df_onlyNonzeroInvariants || currentInvariants[i]) {
    // RDKit❗✔️:         result.push_back(
    // RDKit❗✔️:             new MorganAtomEnv<OutputType>(currentInvariants[i], i, 0, &mol));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // now do our subsequent rounds:
    // RDKit❗✔️:   for (unsigned int layer = 0; layer < morganArguments->d_radius; ++layer) {
    // RDKit❗✔️:     std::vector<AccumTuple> allNeighborhoodsThisRound;
    // RDKit❗✔️:     for (auto atomIdx : atomOrder) {
    // RDKit❗✔️:       // skip atoms which will not generate unique environments
    // RDKit❗✔️:       // (neighborhoods) anymore
    // RDKit❗✔️:       if (!deadAtoms[atomIdx]) {
    // RDKit❗✔️:         const Atom *tAtom = mol.getAtomWithIdx(atomIdx);
    // RDKit❗✔️:         if (!tAtom->getDegree()) {
    // RDKit❗✔️:           deadAtoms.set(atomIdx, 1);
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         ROMol::OEDGE_ITER beg, end;
    // RDKit❗✔️:         boost::tie(beg, end) = mol.getAtomBonds(tAtom);
    // RDKit❗✔️:
    // RDKit❗✔️:         // add up to date invariants of neighbors
    // RDKit❗✔️:         // This should keep capacity, so reallocation only triggers if we
    // RDKit❗✔️:         // haven't seen a molecule of this size.
    // RDKit❗✔️:         neighborhoodInvariants.clear();
    // RDKit❗✔️:
    // RDKit❗✔️:         while (beg != end) {
    // RDKit❗✔️:           const Bond *bond = mol[*beg];
    // RDKit❗✔️:           roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;
    // RDKit❗✔️:
    // RDKit❗✔️:           unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:           roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];
    // RDKit❗✔️:
    // RDKit❗✔️:           auto bt = static_cast<int32_t>((*bondInvariants)[bond->getIdx()]);
    // RDKit❗✔️:           neighborhoodInvariants.push_back(
    // RDKit❗✔️:               std::make_pair(bt, currentInvariants[oIdx]));
    // RDKit❗✔️:
    // RDKit❗✔️:           ++beg;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         // sort the neighbor list:
    // RDKit❗✔️:         std::sort(neighborhoodInvariants.begin(), neighborhoodInvariants.end());
    // RDKit❗✔️:         // and now calculate the new invariant and test if the atom is newly
    // RDKit❗✔️:         // "chiral"
    // RDKit❗✔️:         std::uint32_t invar = layer;
    // RDKit❗✔️:         gboost::hash_combine(invar, currentInvariants[atomIdx]);
    // RDKit❗✔️:         bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:         for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
    // RDKit❗✔️:                  neighborhoodInvariants.begin();
    // RDKit❗✔️:              it != neighborhoodInvariants.end(); ++it) {
    // RDKit❗✔️:           // add the contribution to the new invariant:
    // RDKit❗✔️:           gboost::hash_combine(invar, *it);
    // RDKit❗✔️:
    // RDKit❗✔️:           // check our "chirality":
    // RDKit❗✔️:           if (morganArguments->df_includeChirality && looksChiral &&
    // RDKit❗✔️:               !chiralAtoms[atomIdx]) {
    // RDKit❗✔️:             if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
    // RDKit❗✔️:               looksChiral = false;
    // RDKit❗✔️:             } else if (it != neighborhoodInvariants.begin() &&
    // RDKit❗✔️:                        it->second == (it - 1)->second) {
    // RDKit❗✔️:               looksChiral = false;
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         if (morganArguments->df_includeChirality && looksChiral) {
    // RDKit❗✔️:           chiralAtoms[atomIdx] = 1;
    // RDKit❗✔️:           // add an extra value to the invariant to reflect chirality:
    // RDKit❗✔️:           std::string cip = "";
    // RDKit❗✔️:           tAtom->getPropIfPresent(common_properties::_CIPCode, cip);
    // RDKit❗✔️:           if (cip == "R") {
    // RDKit❗✔️:             gboost::hash_combine(invar, 3);
    // RDKit❗✔️:           } else if (cip == "S") {
    // RDKit❗✔️:             gboost::hash_combine(invar, 2);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             gboost::hash_combine(invar, 1);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         // this rounds bit id will be next rounds atom invariant, so we save
    // RDKit❗✔️:         // it here
    // RDKit❗✔️:         nextLayerInvariants[atomIdx] = static_cast<OutputType>(invar);
    // RDKit❗✔️:
    // RDKit❗✔️:         // store the environment that generated this bit id along with the bit
    // RDKit❗✔️:         // id and the atom id
    // RDKit❗✔️:         allNeighborhoodsThisRound.push_back(
    // RDKit❗✔️:             std::make_tuple(roundAtomNeighborhoods[atomIdx],
    // RDKit❗✔️:                             static_cast<OutputType>(invar), atomIdx));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     std::sort(allNeighborhoodsThisRound.begin(),
    // RDKit❗✔️:               allNeighborhoodsThisRound.end());
    // RDKit❗✔️:     for (std::vector<AccumTuple>::const_iterator iter =
    // RDKit❗✔️:              allNeighborhoodsThisRound.begin();
    // RDKit❗✔️:          iter != allNeighborhoodsThisRound.end(); ++iter) {
    // RDKit❗✔️:       // if we haven't seen this exact environment before, add it to the
    // RDKit❗✔️:       // result
    // RDKit❗✔️:       if (morganArguments->df_includeRedundantEnvironments ||
    // RDKit❗✔️:           neighborhoods.count(std::get<0>(*iter)) == 0) {
    // RDKit❗✔️:         if (!morganArguments->df_onlyNonzeroInvariants ||
    // RDKit❗✔️:             (*atomInvariants)[std::get<2>(*iter)]) {
    // RDKit❗✔️:           if (includeAtoms[std::get<2>(*iter)]) {
    // RDKit❗✔️:             result.push_back(new MorganAtomEnv<OutputType>(
    // RDKit❗✔️:                 std::get<1>(*iter), std::get<2>(*iter), layer + 1, &mol));
    // RDKit❗✔️:             neighborhoods.insert(std::get<0>(*iter));
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // we have seen this exact environment before, this atom
    // RDKit❗✔️:         // is now out of consideration:
    // RDKit❗✔️:         deadAtoms[std::get<2>(*iter)] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // the invariants from this round become the next round invariants:
    // RDKit❗✔️:     currentInvariants.swap(nextLayerInvariants);
    // RDKit❗✔️:     std::fill(nextLayerInvariants.begin(), nextLayerInvariants.end(), 0);
    // RDKit❗✔️:
    // RDKit❗✔️:     // this rounds calculated neighbors will be next rounds initial neighbors,
    // RDKit❗✔️:     // so the radius can grow every iteration
    // RDKit❗✔️:     atomNeighborhoods = roundAtomNeighborhoods;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments
    // Behavior review: Rust checks invariant slice lengths in source order; its
    // typed generator makes the source dynamic-cast wrong-type state impossible.
    // The fixed profile uses legacy stereo perception, so the copied modern
    // CIPLabeler branch remains explicitly unsupported.
    // Rust delegates the source fromAtoms mask to selected_atom_mask, radius-zero
    // emission to radius_zero_environments, per-center expansion to
    // update_neighbor_layer, and sorted round collection/dedup/filter decisions
    // to collect_morgan_layer. Those helpers own these exact source phases; this
    // complete enclosing source anchor retains their original control-flow context.
    // Current/next invariants, persistent seen neighborhoods, dead/chiral masks,
    // fixed initial atom order and cumulative row masks preserve the described
    // phase order. Original invariants remain separate for nonzero filtering.
    // Output codes widen only when OutputType requires it; hash work stays u32.
    // Complexity review: reserve uses source u32 wrapping arithmetic; result
    // and seen set reserve O((r+1)A), row masks use O(A·B/64), atom state and
    // order use O(A), each layer creates one up-to-A candidate vector, the
    // degree scratch is reused, tuple/set costs match the audited helpers, and
    // row handoff remains O(A·B/64) per radius. Contiguous environment values
    // avoid source's separate allocation per emitted environment.
    let atom_count = topology.atoms.len();
    if atom_invariants.len() < atom_count {
        return Err(FingerprintError::PreconditionViolation {
            what: "bad atom invariants size",
        }
        .into());
    }
    if bond_invariants.len() < topology.bonds.len() {
        return Err(FingerprintError::PreconditionViolation {
            what: "bad bond invariants size",
        }
        .into());
    }

    let atom_count_source = atom_count as u32;
    let max_num_results = generator
        .radius
        .wrapping_add(1)
        .wrapping_mul(atom_count_source) as usize;
    let mut result = Vec::with_capacity(max_num_results);

    // The pinned packet profile uses legacy stereo perception. The source's
    // modern CIPLabeler branch above remains explicitly outside this profile.
    let mut current_invariants = atom_invariants.to_vec();
    let mut next_layer_invariants = vec![0; atom_count];
    let mut neighborhood_invariants = Vec::with_capacity(8);
    let include_atoms = selected_atom_mask(atom_count, arguments.from_atoms);
    let mut chiral_atoms = MorganBondEnvironment::new(atom_count);
    let mut neighborhoods = HashSet::with_capacity(max_num_results);

    let bond_count = topology.bonds.len();
    let mut atom_neighborhoods = (0..atom_count)
        .map(|_| MorganBondEnvironment::new(bond_count))
        .collect::<Vec<_>>();
    let mut round_atom_neighborhoods = atom_neighborhoods.clone();
    let mut dead_atoms = MorganBondEnvironment::new(atom_count);

    let mut atom_order = Vec::with_capacity(atom_count);
    if generator.only_nonzero_invariants {
        let mut ordering = Vec::with_capacity(atom_count);
        for atom_index in 0..atom_count_source {
            let is_zero = current_invariants[atom_index as usize] == 0;
            ordering.push((if is_zero { 1 } else { 0 }, atom_index));
        }
        ordering.sort_unstable();
        atom_order.extend(ordering.into_iter().map(|(_, atom_index)| atom_index));
    } else {
        atom_order.extend(0..atom_count_source);
    }

    radius_zero_environments::<OutputType>(
        topology,
        &current_invariants,
        &include_atoms,
        arguments.ignore_atoms,
        generator.only_nonzero_invariants,
        &mut result,
    );

    for layer in 0..generator.radius {
        let mut all_neighborhoods_this_round = Vec::new();
        for &source_atom_index in &atom_order {
            let atom_index = source_atom_index as usize;
            if !dead_atoms.contains(atom_index) {
                let Some(invariant) = update_neighbor_layer(
                    topology,
                    atom_index,
                    layer,
                    generator.fingerprint_arguments.include_chirality,
                    &mut chiral_atoms,
                    &current_invariants,
                    bond_invariants,
                    &atom_neighborhoods,
                    &mut round_atom_neighborhoods,
                    &mut neighborhood_invariants,
                ) else {
                    dead_atoms.set(atom_index);
                    continue;
                };

                next_layer_invariants[atom_index] = invariant;
                all_neighborhoods_this_round.push((
                    round_atom_neighborhoods[atom_index].clone(),
                    invariant,
                    source_atom_index,
                ));
            }
        }

        collect_morgan_layer::<OutputType>(
            &mut all_neighborhoods_this_round,
            generator.include_redundant_environments,
            generator.only_nonzero_invariants,
            atom_invariants,
            &include_atoms,
            &mut neighborhoods,
            &mut dead_atoms,
            layer,
            topology,
            &mut result,
        );

        std::mem::swap(&mut current_invariants, &mut next_layer_invariants);
        next_layer_invariants.fill(0);
        atom_neighborhoods.clone_from(&round_atom_neighborhoods);
    }

    Ok(result)
}

/// Advance one live atom's Morgan neighborhood and compute its base layer code.
/// The caller owns and reuses `neighborhood_invariants` across atom rows and
/// copies the previous neighborhood slice into `round_atom_neighborhoods`
/// once per source radius. `None` is the source isolated-atom dead/continue
/// branch; its caller leaves the pre-zeroed next-layer invariant untouched.
fn update_neighbor_layer(
    topology: &TopologyBlock,
    atom_index: usize,
    layer: u32,
    include_chirality: bool,
    chiral_atoms: &mut MorganChiralAtoms,
    current_invariants: &[u32],
    bond_invariants: &[u32],
    atom_neighborhoods: &[MorganBondEnvironment],
    round_atom_neighborhoods: &mut [MorganBondEnvironment],
    neighborhood_invariants: &mut Vec<(i32, u32)>,
) -> Option<u32> {
    // BEGIN RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments neighbor layer
    // RDKit❗✔️:         const Atom *tAtom = mol.getAtomWithIdx(atomIdx);
    // RDKit❗✔️:         if (!tAtom->getDegree()) {
    // RDKit❗✔️:           deadAtoms.set(atomIdx, 1);
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         ROMol::OEDGE_ITER beg, end;
    // RDKit❗✔️:         boost::tie(beg, end) = mol.getAtomBonds(tAtom);
    // RDKit❗✔️:
    // RDKit❗✔️:         // add up to date invariants of neighbors
    // RDKit❗✔️:         // This should keep capacity, so reallocation only triggers if we
    // RDKit❗✔️:         // haven't seen a molecule of this size.
    // RDKit❗✔️:         neighborhoodInvariants.clear();
    // RDKit❗✔️:
    // RDKit❗✔️:         while (beg != end) {
    // RDKit❗✔️:           const Bond *bond = mol[*beg];
    // RDKit❗✔️:           roundAtomNeighborhoods[atomIdx][bond->getIdx()] = 1;
    // RDKit❗✔️:
    // RDKit❗✔️:           unsigned int oIdx = bond->getOtherAtomIdx(atomIdx);
    // RDKit❗✔️:           roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];
    // RDKit❗✔️:
    // RDKit❗✔️:           auto bt = static_cast<int32_t>((*bondInvariants)[bond->getIdx()]);
    // RDKit❗✔️:           neighborhoodInvariants.push_back(
    // RDKit❗✔️:               std::make_pair(bt, currentInvariants[oIdx]));
    // RDKit❗✔️:
    // RDKit❗✔️:           ++beg;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         // sort the neighbor list:
    // RDKit❗✔️:         std::sort(neighborhoodInvariants.begin(), neighborhoodInvariants.end());
    // RDKit❗✔️:         // and now calculate the new invariant and test if the atom is newly
    // RDKit❗✔️:         // "chiral"
    // RDKit❗✔️:         std::uint32_t invar = layer;
    // RDKit❗✔️:         gboost::hash_combine(invar, currentInvariants[atomIdx]);
    // RDKit❗✔️:         bool looksChiral = (tAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:         for (std::vector<std::pair<int32_t, uint32_t>>::const_iterator it =
    // RDKit❗✔️:                  neighborhoodInvariants.begin();
    // RDKit❗✔️:              it != neighborhoodInvariants.end(); ++it) {
    // RDKit❗✔️:           // add the contribution to the new invariant:
    // RDKit❗✔️:           gboost::hash_combine(invar, *it);
    // RDKit❗✔️:
    // RDKit❗✔️:           // check our "chirality":
    // RDKit❗✔️:           if (morganArguments->df_includeChirality && looksChiral &&
    // RDKit❗✔️:               !chiralAtoms[atomIdx]) {
    // RDKit❗✔️:             if (it->first != static_cast<int32_t>(Bond::SINGLE)) {
    // RDKit❗✔️:               looksChiral = false;
    // RDKit❗✔️:             } else if (it != neighborhoodInvariants.begin() &&
    // RDKit❗✔️:                        it->second == (it - 1)->second) {
    // RDKit❗✔️:               looksChiral = false;
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         if (morganArguments->df_includeChirality && looksChiral) {
    // RDKit❗✔️:           chiralAtoms[atomIdx] = 1;
    // RDKit❗✔️:           // add an extra value to the invariant to reflect chirality:
    // RDKit❗✔️:           std::string cip = "";
    // RDKit❗✔️:           tAtom->getPropIfPresent(common_properties::_CIPCode, cip);
    // RDKit❗✔️:           if (cip == "R") {
    // RDKit❗✔️:             gboost::hash_combine(invar, 3);
    // RDKit❗✔️:           } else if (cip == "S") {
    // RDKit❗✔️:             gboost::hash_combine(invar, 2);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             gboost::hash_combine(invar, 1);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // END RDKIT CPP FUNCTION MorganEnvGenerator::getEnvironments neighbor layer
    // Source helper closure from RDGeneral/hash/hash.hpp:
    // BEGIN RDKIT CPP FUNCTION RDGeneral/hash/hash_combine
    // RDKit❗✔️: #if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: inline void hash_combine(std::hash_result_t& seed, T& v)
    // RDKit❗✔️: #else
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: inline void hash_combine(std::hash_result_t& seed, T const& v)
    // RDKit❗✔️: #endif
    // RDKit❗✔️: {
    // RDKit❗✔️:   gboost::hash<T> hasher;
    // RDKit❗✔️:   seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDGeneral/hash/hash_combine
    // BEGIN RDKIT CPP FUNCTION RDGeneral/hash/hash_value<pair>
    // RDKit❗✔️: template <class A, class B>
    // RDKit❗✔️: std::hash_result_t hash_value(std::pair<A, B> const& v) {
    // RDKit❗✔️:   std::hash_result_t seed = 0;
    // RDKit❗✔️:   hash_combine(seed, v.first);
    // RDKit❗✔️:   hash_combine(seed, v.second);
    // RDKit❗✔️:   return seed;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDGeneral/hash/hash_value<pair>
    // BEGIN RDKIT CPP FUNCTION RDProps::getPropIfPresent<string>
    // RDKit❗🔝:   template <typename T>
    // RDKit❗🔝:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit❗🔝:     return d_props.getValIfPresent(key, res);
    // RDKit❗🔝:   }
    // END RDKIT CPP FUNCTION RDProps::getPropIfPresent<string>
    // BEGIN RDKIT CPP FUNCTION Dict::getValIfPresent<string>
    // RDKit❗🔝:   bool getValIfPresent(const std::string_view what, std::string &res) const {
    // RDKit❗🔝:     for (const auto &i : _data) {
    // RDKit❗🔝:       if (i.key == what) {
    // RDKit❗🔝:         rdvalue_tostring(i.val, res);
    // RDKit❗🔝:         return true;
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:     return false;
    // RDKit❗🔝:   }
    // END RDKIT CPP FUNCTION Dict::getValIfPresent<string>
    // BEGIN RDKIT CPP FUNCTION RDValue::rdvalue_tostring modeled and unsupported tags
    // RDKit❗🔝: inline bool rdvalue_tostring(RDValue_cast_t val, std::string &res) {
    // RDKit❗🔝:   switch (val.getTag()) {
    // RDKit❗🔝:     case RDTypeTag::StringTag:
    // RDKit❗🔝:       res = rdvalue_cast<std::string>(val);
    // RDKit❗🔝:       break;
    // RDKit❗🔝:     case RDTypeTag::IntTag:
    // RDKit❗🔝:       res = boost::lexical_cast<std::string>(rdvalue_cast<int>(val));
    // RDKit❗🔝:       break;
    // RDKit❗🔝:     case RDTypeTag::DoubleTag: {
    // RDKit❗🔝:       Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit❗🔝:       res = boost::lexical_cast<std::string>(rdvalue_cast<double>(val));
    // RDKit❗🔝:       break;
    // RDKit❗🔝:     }
    // RDKit❌❌:     case RDTypeTag::UnsignedIntTag:
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<unsigned int>(val));
    // RDKit❌❌:       break;
    // RDKit❗🔝: #ifdef RDVALUE_HASBOOL
    // RDKit❗🔝:     case RDTypeTag::BoolTag:
    // RDKit❗🔝:       res = boost::lexical_cast<std::string>(rdvalue_cast<bool>(val));
    // RDKit❗🔝:       break;
    // RDKit❗🔝: #endif
    // RDKit❌❌:     case RDTypeTag::FloatTag: {
    // RDKit❌❌:       Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit❌❌:       res = boost::lexical_cast<std::string>(rdvalue_cast<float>(val));
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::VecDoubleTag: {
    // RDKit❌❌:       // vectToString uses std::imbue for locale
    // RDKit❌❌:       res = vectToString<double>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::VecFloatTag: {
    // RDKit❌❌:       // vectToString uses std::imbue for locale
    // RDKit❌❌:       res = vectToString<float>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     case RDTypeTag::VecIntTag:
    // RDKit❌❌:       res = vectToString<int>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::VecUnsignedIntTag:
    // RDKit❌❌:       res = vectToString<unsigned int>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::VecStringTag:
    // RDKit❌❌:       res = vectToString<std::string>(val);
    // RDKit❌❌:       break;
    // RDKit❌❌:     case RDTypeTag::AnyTag: {
    // RDKit❌❌:       Utils::LocaleSwitcher ls;  // for lexical cast...
    // RDKit❌❌:       try {
    // RDKit❌❌:         res = std::any_cast<std::string>(rdvalue_cast<std::any &>(val));
    // RDKit❌❌:       } catch (const std::bad_any_cast &) {
    // RDKit❌❌:         auto &rdtype = rdvalue_cast<std::any &>(val).type();
    // RDKit❌❌:         if (rdtype == typeid(long)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<long>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else if (rdtype == typeid(int64_t)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<int64_t>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else if (rdtype == typeid(uint64_t)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<uint64_t>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else if (rdtype == typeid(unsigned long)) {
    // RDKit❌❌:           res = boost::lexical_cast<std::string>(
    // RDKit❌❌:               std::any_cast<unsigned long>(rdvalue_cast<std::any &>(val)));
    // RDKit❌❌:         } else {
    // RDKit❌❌:           throw;
    // RDKit❌❌:           return false;
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:       break;
    // RDKit❌❌:     }
    // RDKit❌❌:     default:
    // RDKit❌❌:       res = "";
    // RDKit❌❌:   }
    // RDKit❗🔝:   return true;
    // RDKit❗🔝: }
    // END RDKIT CPP FUNCTION RDValue::rdvalue_tostring modeled and unsupported tags
    // This typed match preserves the final R/S/default salt without building
    // the source's temporary string: in the frozen model, Int, Double, and
    // Bool can never compare equal to R or S. Adding a new PropertyValue
    // variant makes this exhaustive match fail to compile until its source
    // conversion is reviewed. The upstream-only unsigned, float, vector, and
    // Any tags stay visibly unsupported rather than defaulting silently.
    // Complexity review: Atom::prop performs one O(log p) BTreeMap lookup
    // instead of the source Dict's O(p) vector scan. Direct tag inspection is
    // O(1) and allocation-free, avoiding lexical conversion while preserving
    // the only values observed by this branch.
    // Behavior review: this helper maps the pinned per-center source slice
    // 390–455 through the chirality salt. It reads incident bonds, unions prior
    // neighborhood rows, sorts signed-bond/u32-atom pairs, hashes in source
    // order, and applies source chirality eligibility plus R/S/default salt.
    // The Rust caller owns later source phases: next-layer assignment, candidate
    // tuple append, round sorting and redundant/seen/nonzero/fromAtoms result
    // decisions, invariant swap/clear, and cumulative neighborhood handoff.
    // None is the source isolated-center return; the caller marks that center
    // dead and leaves its already-zero next-layer invariant in place.
    // For modeled property kinds, only String can equal R or S; Int, Double
    // and Bool retain default salt 1. Unsigned, float, vector and Any tags remain
    // explicitly unsupported under the current model.
    // Complexity review: source and Rust use one reusable O(d) scratch row,
    // O(d log d) pair sorting, O(d) hashing, and O(E/64) packed-word union per
    // incident neighbor. The fixed LP64 block width matches Boost's default
    // `unsigned long`; no per-neighbor hash/vector allocation or topology scan
    // is introduced.
    assert!(
        atom_index < topology.atoms.len(),
        "source atom row is in range"
    );
    assert!(
        current_invariants.len() >= topology.atoms.len(),
        "source atom invariant precondition is satisfied"
    );
    assert!(
        bond_invariants.len() >= topology.bonds.len(),
        "source bond invariant precondition is satisfied"
    );
    assert_eq!(atom_neighborhoods.len(), topology.atoms.len());
    assert_eq!(round_atom_neighborhoods.len(), topology.atoms.len());
    assert_eq!(chiral_atoms.bit_count, topology.atoms.len());
    let neighbors = topology.adjacency.neighbors_of(atom_index);
    if neighbors.is_empty() {
        return None;
    }

    neighborhood_invariants.clear();
    if neighborhood_invariants.capacity() < 8 {
        neighborhood_invariants.reserve(8);
    }
    for neighbor in neighbors {
        let bond_index = neighbor.bond.index();
        round_atom_neighborhoods[atom_index].set(bond_index);
        round_atom_neighborhoods[atom_index].union_with(&atom_neighborhoods[neighbor.atom_index]);
        neighborhood_invariants.push((
            bond_invariants[bond_index] as i32,
            current_invariants[neighbor.atom_index],
        ));
    }

    neighborhood_invariants.sort_unstable();
    let mut invariant = layer;
    hash_combine(&mut invariant, current_invariants[atom_index]);
    let atom = &topology.atoms[atom_index];
    let mut looks_chiral = atom.chiral_tag() != ChiralTag::Unspecified;
    for (pair_index, &(bond_invariant, neighbor_invariant)) in
        neighborhood_invariants.iter().enumerate()
    {
        let mut pair_hash = 0;
        hash_combine(&mut pair_hash, hash_value_i32(bond_invariant));
        hash_combine(&mut pair_hash, hash_value_u32(neighbor_invariant));
        hash_combine(&mut invariant, pair_hash);

        if include_chirality && looks_chiral && !chiral_atoms.contains(atom_index) {
            if bond_invariant != BondOrder::Single.rdkit_code() as i32 {
                looks_chiral = false;
            } else if pair_index != 0
                && neighbor_invariant == neighborhood_invariants[pair_index - 1].1
            {
                looks_chiral = false;
            }
        }
    }

    if include_chirality && looks_chiral {
        chiral_atoms.set(atom_index);
        let chirality_salt = match atom.prop("_CIPCode") {
            Some(PropertyValue::String(cip)) => match cip.as_str() {
                "R" => 3,
                "S" => 2,
                _ => 1,
            },
            Some(PropertyValue::Int(_))
            | Some(PropertyValue::Double(_))
            | Some(PropertyValue::Bool(_))
            | None => 1,
        };
        hash_combine(&mut invariant, chirality_salt);
    }

    Some(invariant)
}

/// Morgan-specific and shared fingerprint generator options.
///
/// This type is defined inside the private `morgan` module. Public visibility
/// here supports the frozen detached owner signature without exposing a new
/// root, runtime, or binding API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganParams {
    pub radius: u32,
    pub include_chirality: bool,
    pub use_bond_types: bool,
    pub include_ring_membership: bool,
    pub only_nonzero_invariants: bool,
    pub include_redundant_environments: bool,
    pub fp_size: u32,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub bits_per_feature: u32,
}

/// Per-call Morgan selection and optional source argument overrides.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MorganCall<'a> {
    pub from_atoms: Option<&'a [u32]>,
    pub ignore_atoms: Option<&'a [u32]>,
    pub custom_atom_invariants: Option<&'a [u32]>,
    pub custom_bond_invariants: Option<&'a [u32]>,
    pub conformer_id: i32,
}

impl Default for MorganCall<'_> {
    fn default() -> Self {
        // RDKit✔️✔️: const std::vector<std::uint32_t> *fromAtoms = nullptr;
        // RDKit✔️✔️: const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
        // RDKit✔️✔️: int confId = -1;
        // RDKit✔️✔️: const std::vector<std::uint32_t> *customAtomInvariants = nullptr;
        // RDKit✔️✔️: const std::vector<std::uint32_t> *customBondInvariants = nullptr;
        Self {
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            conformer_id: -1,
        }
    }
}

impl<'a> MorganCall<'a> {
    fn fingerprint_func_arguments(self) -> FingerprintFuncArguments<'a> {
        FingerprintFuncArguments::new(
            self.from_atoms,
            self.ignore_atoms,
            self.custom_atom_invariants,
            self.custom_bond_invariants,
            self.conformer_id,
        )
    }
}

/// Select the source atom-invariant generator for a canonical Morgan call.
#[derive(Debug, Clone, Copy)]
pub enum MorganAtomInvariants<'a> {
    Connectivity,
    Features,
    FeaturePatterns(&'a [cosmolkit_search::QueryGraph]),
}

impl Default for MorganParams {
    fn default() -> Self {
        // BEGIN RDKIT CPP FUNCTION MorganArguments::MorganArguments source defaults
        // RDKit❗✔️: MorganArguments(unsigned int radius = 3, bool countSimulation = false,
        // RDKit❗✔️:                 bool includeChirality = false,
        // RDKit❗✔️:                 bool onlyNonzeroInvariants = false,
        // RDKit❗✔️:                 std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
        // RDKit❗✔️:                 std::uint32_t fpSize = 2048,
        // RDKit❗✔️:                 bool includeRedundantEnvironments = false,
        // RDKit❗✔️:                 bool useBondTypes = true)
        // RDKit❗✔️:     : FingerprintArguments(countSimulation, countBounds, fpSize, 1,
        // RDKit❗✔️:                            includeChirality),
        // RDKit❗✔️:       df_onlyNonzeroInvariants(onlyNonzeroInvariants),
        // RDKit❗✔️:       d_radius(radius),
        // RDKit❗✔️:       df_includeRedundantEnvironments(includeRedundantEnvironments),
        // RDKit❗✔️:       df_useBondTypes(useBondTypes) {};
        // END RDKIT CPP FUNCTION MorganArguments::MorganArguments source defaults
        // The independent default atom generator contributes ring membership
        // by default; bits-per-feature is fixed to one by the source base call.
        Self {
            radius: 3,
            include_chirality: false,
            use_bond_types: true,
            include_ring_membership: true,
            only_nonzero_invariants: false,
            include_redundant_environments: false,
            fp_size: 2048,
            count_simulation: false,
            count_bounds: vec![1, 2, 4, 8],
            bits_per_feature: 1,
        }
    }
}

/// A configured detached Morgan generator with source-selected defaults.
///
/// The common arguments are stored in their shared owner. The remaining
/// source fields select one Morgan environment generator and its nonzero and
/// redundant-environment policies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct MorganGenerator {
    pub(crate) radius: u32,
    pub(crate) only_nonzero_invariants: bool,
    pub(crate) include_redundant_environments: bool,
    pub(crate) fingerprint_arguments: FingerprintArguments,
    pub(crate) atom_invariants: MorganAtomInvGenerator,
    pub(crate) bond_invariants: MorganBondInvGenerator,
}

/// Construct the Morgan generator configuration and its default invariant
/// owners. Per-call custom invariant slices remain on `MorganCall`.
pub(crate) fn get_morgan_generator(params: &MorganParams) -> Result<MorganGenerator, MorganError> {
    // BEGIN RDKIT CPP FUNCTION getMorganGenerator(const MorganArguments &, ...)
    // RDKit❗✔️:   AtomEnvironmentGenerator<OutputType> *morganEnvGenerator =
    // RDKit❗✔️:       new MorganEnvGenerator<OutputType>();
    // RDKit❗✔️:   bool ownsAtomInvGenerator = ownsAtomInvGen;
    // RDKit❗✔️:   if (!atomInvariantsGenerator) {
    // RDKit❗✔️:     atomInvariantsGenerator = new MorganAtomInvGenerator();
    // RDKit❗✔️:     ownsAtomInvGenerator = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool ownsBondInvGenerator = ownsBondInvGen;
    // RDKit❗✔️:   if (!bondInvariantsGenerator) {
    // RDKit❗✔️:     bondInvariantsGenerator = new MorganBondInvGenerator(
    // RDKit❗✔️:         args.df_useBondTypes, args.df_includeChirality);
    // RDKit❗✔️:     ownsBondInvGenerator = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return new FingerprintGenerator<OutputType>(
    // RDKit❗✔️:       morganEnvGenerator, new MorganArguments(args),
    // RDKit❗✔️:       atomInvariantsGenerator, bondInvariantsGenerator,
    // RDKit❗✔️:       ownsAtomInvGenerator, ownsBondInvGenerator);
    // END RDKIT CPP FUNCTION getMorganGenerator(const MorganArguments &, ...)
    // Behavior review: this detached entry always selects the source default
    // invariant owners because the frozen API has no polymorphic generator
    // pointer inputs; its custom per-call invariant slices are preserved by
    // `MorganCall`. The copied common arguments reproduce the source's
    // MorganArguments copy and re-run its two base preconditions even if a
    // caller mutates public parameter fields after `Default` construction.
    // Complexity review: one count-bounds clone models the source argument
    // copy. Both invariant configurations are stored inline, avoiding the
    // source's two default polymorphic allocations without a graph scan or
    // per-molecule work. Radius/flags are copied in constant time.
    let fingerprint_arguments = FingerprintArguments::new(
        params.count_simulation,
        params.count_bounds.clone(),
        params.fp_size,
        params.bits_per_feature,
        params.include_chirality,
    )?;

    Ok(MorganGenerator {
        radius: params.radius,
        only_nonzero_invariants: params.only_nonzero_invariants,
        include_redundant_environments: params.include_redundant_environments,
        fingerprint_arguments,
        atom_invariants: MorganAtomInvGenerator::new(params.include_ring_membership),
        bond_invariants: MorganBondInvGenerator::new(
            params.use_bond_types,
            params.include_chirality,
        ),
    })
}

fn selected_atom_invariants(
    topology: &TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    selection: MorganAtomInvariants<'_>,
) -> Result<Vec<u32>, MorganError> {
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper atom invariant selection
    // RDKit❗✔️:   if (args.customAtomInvariants) {
    // RDKit❗✔️:     atomInvariants.reset(
    // RDKit❗✔️:         new std::vector<std::uint32_t>(*args.customAtomInvariants));
    // RDKit❗✔️:   } else if (dp_atomInvariantsGenerator) {
    // RDKit❗✔️:     atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper atom invariant selection
    // Behavior review: the shared caller invokes this only when a custom
    // atom-invariant slice is absent. Connectivity uses the configured
    // source owner and supplied assignments; feature modes use the retained
    // default six SMARTS or the present custom query slice on the original
    // topology, with the supplied final rings and valence borrowed by Search.
    // Complexity review: connectivity remains one O(A) owner pass. Feature
    // defaults reuse cached parsed queries; prepared query context borrows the
    // provided assignments without recalculation, and matching is delegated
    // to the one existing Search owner.
    match selection {
        MorganAtomInvariants::Connectivity => generator
            .atom_invariants
            .get_atom_invariants(topology, valence, rings),
        MorganAtomInvariants::Features | MorganAtomInvariants::FeaturePatterns(_) => {
            let patterns = match selection {
                MorganAtomInvariants::Connectivity => unreachable!(),
                MorganAtomInvariants::Features => None,
                MorganAtomInvariants::FeaturePatterns(patterns) => Some(patterns),
            };
            let target = SearchTarget::new(
                topology,
                coordinates,
                &topology.stereo_groups,
                Some(rings),
                Some(valence),
            );
            let query_context = build_prepared_query_match_context(topology, rings, valence)?;
            MorganFeatureAtomInvGenerator::new(patterns)
                .get_atom_invariants(&target, &query_context)
        }
    }
}

fn create_morgan_call<'call>(
    params: &MorganParams,
    call: &MorganCall<'call>,
) -> Result<(MorganGenerator, FingerprintFuncArguments<'call>), MorganError> {
    let generator = get_morgan_generator(params)?;
    let arguments = call.fingerprint_func_arguments();
    Ok((generator, arguments))
}

/// Generate raw sparse Morgan counts from final detached molecule state.
pub fn morgan_sparse_count(
    input: &MorganPreparedInput<'_>,
    params: &MorganParams,
    call: &MorganCall<'_>,
    invariants: MorganAtomInvariants<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<crate::SparseCountFingerprint, MorganError> {
    let (generator, arguments) = create_morgan_call(params, call)?;
    get_sparse_count_fingerprint_with_atom_invariants(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        &generator,
        &arguments,
        output,
        |topology, _properties, valence, rings, generator| {
            selected_atom_invariants(
                topology,
                input.coordinates,
                valence,
                rings,
                generator,
                invariants,
            )
        },
    )
}

/// Generate raw sparse Morgan presence bits from final detached molecule state.
pub fn morgan_sparse_bits(
    input: &MorganPreparedInput<'_>,
    params: &MorganParams,
    call: &MorganCall<'_>,
    invariants: MorganAtomInvariants<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<crate::SparseBitFingerprint, MorganError> {
    let (generator, arguments) = create_morgan_call(params, call)?;
    get_sparse_fingerprint_with_atom_invariants(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        &generator,
        &arguments,
        output,
        |topology, _properties, valence, rings, generator| {
            selected_atom_invariants(
                topology,
                input.coordinates,
                valence,
                rings,
                generator,
                invariants,
            )
        },
    )
}

/// Generate a hashed 32-bit sparse Morgan count fingerprint.
pub fn morgan_count(
    input: &MorganPreparedInput<'_>,
    params: &MorganParams,
    call: &MorganCall<'_>,
    invariants: MorganAtomInvariants<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<SparseCountFingerprint32, MorganError> {
    let (generator, arguments) = create_morgan_call(params, call)?;
    get_count_fingerprint_with_atom_invariants(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        &generator,
        &arguments,
        output,
        |topology, _properties, valence, rings, generator| {
            selected_atom_invariants(
                topology,
                input.coordinates,
                valence,
                rings,
                generator,
                invariants,
            )
        },
    )
}

/// Generate a dense Morgan bit fingerprint.
pub fn morgan_bits(
    input: &MorganPreparedInput<'_>,
    params: &MorganParams,
    call: &MorganCall<'_>,
    invariants: MorganAtomInvariants<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<Fingerprint, MorganError> {
    let (generator, arguments) = create_morgan_call(params, call)?;
    get_fingerprint_with_atom_invariants(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        &generator,
        &arguments,
        output,
        |topology, _properties, valence, rings, generator| {
            selected_atom_invariants(
                topology,
                input.coordinates,
                valence,
                rings,
                generator,
                invariants,
            )
        },
    )
}

/// Detached legacy `MorganFingerprints::getFingerprint` projection.
///
/// This private owner preserves the old wrapper's two sparse result branches;
/// it is not the later public `MorganParams` projection API.
pub(crate) fn get_legacy_morgan_fingerprint(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    radius: u32,
    use_chirality: bool,
    use_bond_types: bool,
    custom_atom_invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    use_counts: bool,
    only_nonzero_invariants: bool,
    atoms_setting_bits: Option<&mut BTreeMap<u32, Vec<(u32, u32)>>>,
    include_redundant_environments: bool,
) -> Result<SparseCountFingerprint32, MorganError> {
    // BEGIN RDKIT CPP FUNCTION MorganFingerprints::getFingerprint
    // RDKit❗✔️:   bool countSimulation = false;
    // RDKit❗✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit❗✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(
    // RDKit❗✔️:           radius, countSimulation, useChirality, useBondTypes,
    // RDKit❗✔️:           onlyNonzeroInvariants, includeRedundantEnvironments));
    // RDKit❗✔️:   RDKit::FingerprintFuncArguments args;
    // RDKit❗✔️:   args.fromAtoms = fromAtoms;
    // RDKit❗✔️:   args.customAtomInvariants = invariants;
    // RDKit❗✔️:   AdditionalOutput ao;
    // RDKit❗✔️:   if (atomsSettingBits) {
    // RDKit❗✔️:     args.additionalOutput = &ao;
    // RDKit❗✔️:     ao.allocateBitInfoMap();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   SparseIntVect<uint32_t> *res;
    // RDKit❗✔️:   if (!useCounts) {
    // RDKit❗✔️:     auto tmp = fpgen->getSparseFingerprint(mol, args);
    // RDKit❗✔️:     res = new SparseIntVect<uint32_t>(std::numeric_limits<uint32_t>::max());
    // RDKit❗✔️:     for (auto idx : *(tmp->dp_bits)) {
    // RDKit❗✔️:       res->setVal(idx, 1);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     res = fpgen->getSparseCountFingerprint(mol, args).release();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (atomsSettingBits) {
    // RDKit❗✔️:     atomsSettingBits->clear();
    // RDKit❗✔️:     for (const auto &pr : *(ao.bitInfoMap)) {
    // RDKit❗✔️:       (*atomsSettingBits)[pr.first] = pr.second;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // END RDKIT CPP FUNCTION MorganFingerprints::getFingerprint
    // Behavior review: the legacy wrapper fixes count simulation off and the
    // source default bits-per-feature to one, selects the full u32 output
    // domain, preserves its raw-count versus folded-bit branches, and commits
    // only the local bitInfo map after fingerprint generation succeeds.
    // The false branch writes value one for every folded nonzero ID. Since a
    // sparse count map has exactly those IDs, projecting its keys directly
    // avoids constructing an intermediate sparse-bit set without changing
    // collision, ordering, or output values. The true branch accumulates
    // directly in the source u32 result map rather than copying a u64 map.
    // Complexity review: both source branches and these projections remain
    // O(E log K) with ordered sparse storage; the direct key projection removes
    // one temporary sparse structure, and neither branch rescans the molecule.
    let mut generator_params = MorganParams::default();
    generator_params.radius = radius;
    generator_params.include_chirality = use_chirality;
    generator_params.use_bond_types = use_bond_types;
    generator_params.only_nonzero_invariants = only_nonzero_invariants;
    generator_params.include_redundant_environments = include_redundant_environments;
    let generator = get_morgan_generator(&generator_params)?;
    let arguments =
        FingerprintFuncArguments::new(from_atoms, None, custom_atom_invariants, None, -1);

    let mut staged_output = atoms_setting_bits.as_ref().map(|_| {
        let mut output = AdditionalOutput::default();
        output.allocate_bit_info_map();
        output
    });
    let output = staged_output.as_mut();
    let fingerprint = with_morgan_environment_inputs_and_output(
        topology,
        properties,
        valence,
        rings,
        &generator,
        &arguments,
        output,
        |prepared_topology, _prepared_properties, atom_invariants, bond_invariants, output| {
            let environments = generate_morgan_environments::<u32>(
                prepared_topology,
                &generator,
                &arguments,
                atom_invariants,
                bond_invariants,
            )?;
            let mut accumulated = SparseCountFingerprint32::new(u32::MAX);
            let fp_size = if use_counts { 0 } else { u64::from(u32::MAX) };
            accumulate_morgan_sparse_counts_into(
                environments,
                &generator.fingerprint_arguments,
                atom_invariants,
                bond_invariants,
                fp_size,
                output,
                &mut accumulated,
            )?;

            if use_counts {
                return Ok(accumulated);
            }

            let mut result = SparseCountFingerprint32::new(u32::MAX);
            for &bit_id in accumulated.nonzero_elements().keys() {
                result.set_value(bit_id, 1)?;
            }
            Ok(result)
        },
    )?;

    if let Some(atoms_setting_bits) = atoms_setting_bits {
        atoms_setting_bits.clear();
        let bit_info_map = staged_output
            .expect("a present source BitInfoMap allocates local AdditionalOutput")
            .bit_info_map
            .expect("the legacy wrapper allocates its local bitInfoMap");
        for (bit_id, provenance) in bit_info_map {
            atoms_setting_bits.insert(bit_id as u32, provenance);
        }
    }

    Ok(fingerprint)
}

/// Detached legacy `MorganFingerprints::getHashedFingerprint` projection.
///
/// The private wrapper preserves the legacy positive-size contract while
/// delegating hashing, accumulation, and sparse result construction to the
/// canonical count-fingerprint owner.
pub(crate) fn get_legacy_morgan_hashed_fingerprint(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    radius: u32,
    n_bits: u32,
    custom_atom_invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    use_chirality: bool,
    use_bond_types: bool,
    only_nonzero_invariants: bool,
    atoms_setting_bits: Option<&mut BTreeMap<u32, Vec<(u32, u32)>>>,
    include_redundant_environments: bool,
) -> Result<SparseCountFingerprint32, MorganError> {
    // BEGIN RDKIT CPP FUNCTION MorganFingerprints::getHashedFingerprint
    // RDKit❗✔️:   if (nBits == 0) {
    // RDKit❗✔️:     throw ValueErrorException("nBits can not be zero");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool countSimulation = false;
    // RDKit❗✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit❗✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(
    // RDKit❗✔️:           radius, countSimulation, useChirality, useBondTypes,
    // RDKit❗✔️:           onlyNonzeroInvariants, includeRedundantEnvironments, nullptr, nullptr,
    // RDKit❗✔️:           nBits));
    // RDKit❗✔️:   RDKit::FingerprintFuncArguments args;
    // RDKit❗✔️:   args.fromAtoms = fromAtoms;
    // RDKit❗✔️:   args.customAtomInvariants = invariants;
    // RDKit❗✔️:   AdditionalOutput ao;
    // RDKit❗✔️:   if (atomsSettingBits) {
    // RDKit❗✔️:     args.additionalOutput = &ao;
    // RDKit❗✔️:     ao.allocateBitInfoMap();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto res = fpgen->getCountFingerprint(mol, args).release();
    // RDKit❗✔️:   if (atomsSettingBits) {
    // RDKit❗✔️:     atomsSettingBits->clear();
    // RDKit❗✔️:     for (const auto &pr : *(ao.bitInfoMap)) {
    // RDKit❗✔️:       (*atomsSettingBits)[pr.first] = pr.second;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // END RDKIT CPP FUNCTION MorganFingerprints::getHashedFingerprint
    // Behavior review: reject nBits zero before generator or output state is
    // created, mapping the source ValueError to the existing typed
    // InvalidArguments reason. For positive nBits, configure exactly the
    // legacy flags and fpSize, then reuse the source-shaped count projection;
    // folded IDs are counted and reported to bitInfo by that shared owner.
    // Like the wrapper, stage only bitInfo when requested and replace the
    // caller map only after fingerprint generation succeeds.
    // Complexity review: option setup is O(B) for the fixed four default
    // count bounds. Generation keeps the canonical single environment pass,
    // ordered count accumulation, and source-shaped temporary-to-u32 sparse
    // map projection; this adapter adds no environment scan or count map.
    // Optional bitInfo costs one local ordered map and the source-equivalent
    // O(K log K) success-only copy to the caller map.
    if n_bits == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "nBits can not be zero",
        }
        .into());
    }

    let mut generator_params = MorganParams::default();
    generator_params.radius = radius;
    generator_params.include_chirality = use_chirality;
    generator_params.use_bond_types = use_bond_types;
    generator_params.only_nonzero_invariants = only_nonzero_invariants;
    generator_params.include_redundant_environments = include_redundant_environments;
    generator_params.fp_size = n_bits;
    let generator = get_morgan_generator(&generator_params)?;
    let arguments =
        FingerprintFuncArguments::new(from_atoms, None, custom_atom_invariants, None, -1);

    let mut staged_output = atoms_setting_bits.as_ref().map(|_| {
        let mut output = AdditionalOutput::default();
        output.allocate_bit_info_map();
        output
    });
    let fingerprint = get_count_fingerprint(
        topology,
        properties,
        valence,
        rings,
        &generator,
        &arguments,
        staged_output.as_mut(),
    )?;

    if let Some(atoms_setting_bits) = atoms_setting_bits {
        atoms_setting_bits.clear();
        let bit_info_map = staged_output
            .expect("a present source BitInfoMap allocates local AdditionalOutput")
            .bit_info_map
            .expect("the legacy wrapper allocates its local bitInfoMap");
        for (bit_id, provenance) in bit_info_map {
            atoms_setting_bits.insert(bit_id as u32, provenance);
        }
    }

    Ok(fingerprint)
}

/// Detached legacy MorganFingerprints::getFingerprintAsBitVect projection.
///
/// The legacy wrapper uses one configured Morgan generator and the canonical
/// dense fingerprint projection; this adapter preserves its size/error and
/// optional bitInfo commit boundary.
pub(crate) fn get_legacy_morgan_fingerprint_as_bit_vector(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    radius: u32,
    n_bits: u32,
    custom_atom_invariants: Option<&[u32]>,
    from_atoms: Option<&[u32]>,
    use_chirality: bool,
    use_bond_types: bool,
    only_nonzero_invariants: bool,
    atoms_setting_bits: Option<&mut BTreeMap<u32, Vec<(u32, u32)>>>,
    include_redundant_environments: bool,
) -> Result<Fingerprint, MorganError> {
    // BEGIN RDKIT CPP FUNCTION MorganFingerprints::getFingerprintAsBitVect
    // RDKit❗✔️: ExplicitBitVect *getFingerprintAsBitVect(
    // RDKit❗✔️:     const ROMol &mol, unsigned int radius, unsigned int nBits,
    // RDKit❗✔️:     std::vector<uint32_t> *invariants, const std::vector<uint32_t> *fromAtoms,
    // RDKit❗✔️:     bool useChirality, bool useBondTypes, bool onlyNonzeroInvariants,
    // RDKit❗✔️:     BitInfoMap *atomsSettingBits, bool includeRedundantEnvironments) {
    // RDKit❗✔️:   if (nBits == 0) {
    // RDKit❗✔️:     throw ValueErrorException("nBits can not be zero");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool countSimulation = false;
    // RDKit❗✔️:   std::unique_ptr<FingerprintGenerator<std::uint32_t>> fpgen(
    // RDKit❗✔️:       MorganFingerprint::getMorganGenerator<std::uint32_t>(
    // RDKit❗✔️:           radius, countSimulation, useChirality, useBondTypes,
    // RDKit❗✔️:           onlyNonzeroInvariants, includeRedundantEnvironments, nullptr, nullptr,
    // RDKit❗✔️:           nBits));
    // RDKit❗✔️:   RDKit::FingerprintFuncArguments args;
    // RDKit❗✔️:   args.fromAtoms = fromAtoms;
    // RDKit❗✔️:   args.customAtomInvariants = invariants;
    // RDKit❗✔️:   AdditionalOutput ao;
    // RDKit❗✔️:   if (atomsSettingBits) {
    // RDKit❗✔️:     args.additionalOutput = &ao;
    // RDKit❗✔️:     ao.allocateBitInfoMap();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto res = fpgen->getFingerprint(mol, args).release();
    // RDKit❗✔️:   if (atomsSettingBits) {
    // RDKit❗✔️:     atomsSettingBits->clear();
    // RDKit❗✔️:     for (const auto &pr : *(ao.bitInfoMap)) {
    // RDKit❗✔️:       (*atomsSettingBits)[pr.first] = pr.second;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MorganFingerprints::getFingerprintAsBitVect
    // Behavior review: reject nBits zero before generator/output construction;
    // fix count simulation off and forward every legacy option to the same
    // Morgan generator. Only a supplied caller map allocates local bitInfo.
    // Reuse the canonical dense projection, then replace the caller map only
    // after successful generation, matching the wrapper's commit boundary.
    // Complexity review: option setup and generator construction are O(1)
    // aside from the fixed four default bounds; dense creation is O(nBits/64)
    // and projection visits K ordered nonzero folded IDs once. Environment
    // accumulation retains the shared O(E log K) pass; optional bitInfo adds
    // one O(K log K + P) success-only map transfer and no second graph scan.
    if n_bits == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "nBits can not be zero",
        }
        .into());
    }

    let mut generator_params = MorganParams::default();
    generator_params.radius = radius;
    generator_params.include_chirality = use_chirality;
    generator_params.use_bond_types = use_bond_types;
    generator_params.only_nonzero_invariants = only_nonzero_invariants;
    generator_params.include_redundant_environments = include_redundant_environments;
    generator_params.fp_size = n_bits;
    let generator = get_morgan_generator(&generator_params)?;
    let arguments =
        FingerprintFuncArguments::new(from_atoms, None, custom_atom_invariants, None, -1);

    let mut staged_output = atoms_setting_bits.as_ref().map(|_| {
        let mut output = AdditionalOutput::default();
        output.allocate_bit_info_map();
        output
    });
    let fingerprint = get_fingerprint(
        topology,
        properties,
        valence,
        rings,
        &generator,
        &arguments,
        staged_output.as_mut(),
    )?;

    if let Some(atoms_setting_bits) = atoms_setting_bits {
        atoms_setting_bits.clear();
        let bit_info_map = staged_output
            .expect("a present source BitInfoMap allocates local AdditionalOutput")
            .bit_info_map
            .expect("the legacy wrapper allocates its local bitInfoMap");
        for (bit_id, provenance) in bit_info_map {
            atoms_setting_bits.insert(bit_id as u32, provenance);
        }
    }

    Ok(fingerprint)
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;
    use std::collections::{BTreeMap, HashSet};
    use std::panic::{AssertUnwindSafe, catch_unwind};

    use super::{
        MorganAtomEnvironment, MorganAtomInvariants, MorganBondEnvironment, MorganCall,
        MorganChiralAtoms, MorganDistanceMatrixCache, MorganLayerCandidate, MorganParams,
        collect_morgan_layer, generate_morgan_environments, get_legacy_morgan_fingerprint,
        get_legacy_morgan_fingerprint_as_bit_vector, get_legacy_morgan_hashed_fingerprint,
        get_morgan_generator, morgan_bits, morgan_count, morgan_sparse_bits, morgan_sparse_count,
        radius_zero_environments, selected_atom_mask, update_neighbor_layer,
    };
    use crate::additional_output::AdditionalOutput;
    use crate::generator::{
        FingerprintArguments, FingerprintFuncArguments, with_morgan_environment_inputs,
    };
    use crate::hash::{hash_combine, hash_value_i32, hash_value_u32};
    use crate::invariants::MorganFeatureAtomInvGenerator;
    use crate::{FingerprintError, MorganError, MorganPreparedInput};
    use cosmolkit_core::{
        RemoveHsParams, SanitizeParams, ValenceAssignment, fast_find_rings,
        remove_hydrogens_with_params, sanitize_topology,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, ChiralTag, Element,
        PropertyValue, TopologyBlock,
    };
    use cosmolkit_search::{
        SearchTarget, SmartsParseParams, build_prepared_query_match_context, parse_smarts,
    };
    use cosmolkit_smiles::{SmilesParseParams, SmilesRecord, finalize_smiles_stereo, parse_smiles};

    fn fully_prepared_record(smiles: &str) -> SmilesRecord {
        fully_prepared_record_with_valence(smiles).0
    }

    fn fully_prepared_record_with_valence(smiles: &str) -> (SmilesRecord, ValenceAssignment) {
        let parse_params = SmilesParseParams::default();
        let parsed = parse_smiles(smiles, &parse_params).expect("fixed Morgan SMILES parses");
        let sanitized = sanitize_topology(&parsed.topology, &SanitizeParams::default())
            .expect("fixed Morgan topology sanitizes");
        let removed = remove_hydrogens_with_params(
            sanitized.topology,
            parsed.coordinates,
            parsed.properties,
            &RemoveHsParams::default(),
        )
        .expect("fixed Morgan hydrogens are removed");
        let mut prepared_valence = removed.final_valence;
        let mut prepared_rings = removed.final_rings;
        let record = finalize_smiles_stereo(
            SmilesRecord {
                topology: removed.topology,
                coordinates: removed.coordinates,
                properties: removed.properties,
            },
            &parse_params,
            &mut prepared_valence,
            &mut prepared_rings,
        )
        .expect("fixed Morgan stereo is finalized");
        let prepared_valence = prepared_valence
            .expect("canonical sanitize/RemoveHs/stereo preparation yields final valence");
        (record, prepared_valence)
    }

    #[test]
    fn fingerprint_morgan_a01_canonical_output_families_match_fixed_source_product() {
        const BASE: &[(u32, i32)] = &[(864_662_311, 1), (2_245_384_272, 1), (2_246_728_737, 1)];
        const RADIUS_ONE: &[(u32, i32)] = &[
            (864_662_311, 1),
            (1_535_166_686, 1),
            (2_245_384_272, 1),
            (2_246_728_737, 1),
            (3_542_456_614, 1),
            (4_018_048_386, 1),
        ];
        const RADIUS_TWO_REDUNDANT: &[(u32, i32)] = &[
            (407_009_239, 1),
            (864_662_311, 1),
            (1_535_166_686, 1),
            (2_245_384_272, 1),
            (2_246_728_737, 1),
            (2_679_480_461, 1),
            (3_542_456_614, 1),
            (3_732_090_711, 1),
            (4_018_048_386, 1),
        ];
        const RADIUS_THREE_REDUNDANT: &[(u32, i32)] = &[
            (407_009_239, 1),
            (589_800_426, 1),
            (864_662_311, 1),
            (1_535_166_686, 1),
            (2_245_384_272, 1),
            (2_246_728_737, 1),
            (2_679_480_461, 1),
            (2_708_961_779, 1),
            (2_962_971_722, 1),
            (3_542_456_614, 1),
            (3_732_090_711, 1),
            (4_018_048_386, 1),
        ];

        fn fixed_row(radius: u32, redundant: bool) -> &'static [(u32, i32)] {
            match (radius, redundant) {
                (0, _) => BASE,
                (1, _) => RADIUS_ONE,
                (2, false) | (3, false) => RADIUS_ONE,
                (2, true) => RADIUS_TWO_REDUNDANT,
                (3, true) => RADIUS_THREE_REDUNDANT,
                _ => unreachable!("the test enumerates the four frozen radii"),
            }
        }

        let (record, valence) = fully_prepared_record_with_valence("CCO");
        let rings = fast_find_rings(&record.topology).expect("fixed CCO ring state");
        let input = MorganPreparedInput {
            topology: &record.topology,
            coordinates: &record.coordinates,
            properties: &record.properties,
            valence: &valence,
            rings: &rings,
        };
        let mut cases = 0;

        // CCO has no ring atoms, stereocenters, or zero default connectivity
        // invariants. The pinned M09 table therefore fixes the full product
        // of these source no-op flags, all four radii, and redundancy.
        for radius in 0..=3 {
            for include_chirality in [false, true] {
                for include_ring_membership in [false, true] {
                    for only_nonzero_invariants in [false, true] {
                        for include_redundant_environments in [false, true] {
                            let params = MorganParams {
                                radius,
                                include_chirality,
                                include_ring_membership,
                                only_nonzero_invariants,
                                include_redundant_environments,
                                ..MorganParams::default()
                            };
                            let call = MorganCall::default();
                            let expected_row = fixed_row(radius, include_redundant_environments);
                            let expected_raw_counts = expected_row
                                .iter()
                                .map(|&(index, count)| (u64::from(index), count))
                                .collect::<BTreeMap<_, _>>();

                            let sparse_count = morgan_sparse_count(
                                &input,
                                &params,
                                &call,
                                MorganAtomInvariants::Connectivity,
                                None,
                            )
                            .expect("fixed raw sparse Morgan count");
                            assert_eq!(
                                sparse_count.nonzero_elements(),
                                &expected_raw_counts,
                                "raw counts: radius={radius}, chirality={include_chirality}, ring_membership={include_ring_membership}, only_nonzero={only_nonzero_invariants}, redundant={include_redundant_environments}"
                            );

                            let sparse_bits = morgan_sparse_bits(
                                &input,
                                &params,
                                &call,
                                MorganAtomInvariants::Connectivity,
                                None,
                            )
                            .expect("fixed raw sparse Morgan bits");
                            let mut expected_raw_bits = expected_row
                                .iter()
                                .map(|&(index, _)| index as i32)
                                .collect::<Vec<_>>();
                            expected_raw_bits.sort_unstable();
                            assert_eq!(
                                sparse_bits.on_bits(),
                                expected_raw_bits,
                                "raw bits: radius={radius}, chirality={include_chirality}, ring_membership={include_ring_membership}, only_nonzero={only_nonzero_invariants}, redundant={include_redundant_environments}"
                            );

                            // FingerprintGenerator::getFingerprintHelper folds
                            // each fixed source ID with `% fpSize`, then the
                            // count projection adds the source count at that
                            // projected key. Expectations start only from
                            // the literal M09 row above, never from a sibling
                            // output family or the implementation under test.
                            let mut expected_hashed_counts = BTreeMap::new();
                            for &(index, count) in expected_row {
                                *expected_hashed_counts
                                    .entry(index % params.fp_size)
                                    .or_insert(0) += count;
                            }
                            let hashed_count = morgan_count(
                                &input,
                                &params,
                                &call,
                                MorganAtomInvariants::Connectivity,
                                None,
                            )
                            .expect("fixed hashed Morgan count");
                            assert_eq!(
                                hashed_count.nonzero_elements(),
                                &expected_hashed_counts,
                                "hashed counts: radius={radius}, chirality={include_chirality}, ring_membership={include_ring_membership}, only_nonzero={only_nonzero_invariants}, redundant={include_redundant_environments}"
                            );

                            let dense_bits = morgan_bits(
                                &input,
                                &params,
                                &call,
                                MorganAtomInvariants::Connectivity,
                                None,
                            )
                            .expect("fixed dense Morgan bits");
                            assert_eq!(
                                dense_bits.on_bits(),
                                expected_hashed_counts.keys().copied().collect::<Vec<_>>(),
                                "dense bits: radius={radius}, chirality={include_chirality}, ring_membership={include_ring_membership}, only_nonzero={only_nonzero_invariants}, redundant={include_redundant_environments}"
                            );
                            cases += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(cases, 64);
    }

    #[test]
    fn fingerprint_morgan_a01_canonical_feature_queries_and_metadata_are_retained() {
        let (record, valence) = fully_prepared_record_with_valence("C");
        let rings = fast_find_rings(&record.topology).expect("fixed carbon ring state");
        let input = MorganPreparedInput {
            topology: &record.topology,
            coordinates: &record.coordinates,
            properties: &record.properties,
            valence: &valence,
            rings: &rings,
        };
        let params = MorganParams {
            radius: 0,
            ..MorganParams::default()
        };
        let call = MorganCall::default();
        let mut output = AdditionalOutput {
            atom_counts: Some(Vec::new()),
            atom_to_bits: Some(Vec::new()),
            bit_info_map: Some(BTreeMap::new()),
            bit_paths: Some(BTreeMap::new()),
            atoms_per_bit: Some(BTreeMap::new()),
        };
        let default_features = morgan_sparse_count(
            &input,
            &params,
            &call,
            MorganAtomInvariants::Features,
            Some(&mut output),
        )
        .expect("canonical call uses the retained six feature patterns");
        assert_eq!(
            default_features.nonzero_elements(),
            &BTreeMap::from([(0, 1)])
        );
        assert_eq!(output.atom_counts, Some(vec![1]));
        assert_eq!(output.atom_to_bits, Some(vec![vec![0]]));
        assert_eq!(
            output.bit_info_map,
            Some(BTreeMap::from([(0, vec![(0, 0)])]))
        );
        assert_eq!(output.bit_paths, Some(BTreeMap::new()));
        assert_eq!(
            output.atoms_per_bit,
            Some(BTreeMap::from([(0, vec![vec![0]])]))
        );
        let repeated_default_features =
            morgan_sparse_count(&input, &params, &call, MorganAtomInvariants::Features, None)
                .expect("repeated canonical calls reuse the retained default query owner");
        assert_eq!(
            repeated_default_features.nonzero_elements(),
            &BTreeMap::from([(0, 1)])
        );

        let parsed_carbon = parse_smarts("[#6]", &SmartsParseParams::default())
            .expect("fixed caller-provided feature SMARTS parses");
        let feature_patterns = [parsed_carbon];
        let custom_feature = morgan_sparse_count(
            &input,
            &params,
            &call,
            MorganAtomInvariants::FeaturePatterns(&feature_patterns),
            None,
        )
        .expect("canonical call consumes the parsed custom query");
        assert_eq!(custom_feature.nonzero_elements(), &BTreeMap::from([(1, 1)]));

        let custom_atom_invariants = [7];
        let overriding_call = MorganCall {
            custom_atom_invariants: Some(&custom_atom_invariants),
            ..MorganCall::default()
        };
        let overridden_feature = morgan_sparse_count(
            &input,
            &params,
            &overriding_call,
            MorganAtomInvariants::FeaturePatterns(&feature_patterns),
            None,
        )
        .expect("source custom atom invariants override the selected query provider");
        assert_eq!(
            overridden_feature.nonzero_elements(),
            &BTreeMap::from([(7, 1)])
        );

        let empty_patterns: [cosmolkit_search::QueryGraph; 0] = [];
        let empty_feature = morgan_sparse_count(
            &input,
            &params,
            &call,
            MorganAtomInvariants::FeaturePatterns(&empty_patterns),
            None,
        )
        .expect("a present empty pattern slice remains distinct from None");
        assert_eq!(empty_feature.nonzero_elements(), &BTreeMap::from([(0, 1)]));
    }

    #[test]
    fn fingerprint_morgan_mask_fix_canonical_product_preserves_errors_and_metadata() {
        const CASES: [(usize, u32, bool); 7] = [
            (0, 0, false),
            (1, 1, false),
            (30, 0x3fff_ffff, false),
            (31, 0x7fff_ffff, false),
            (32, 0xffff_ffff, false),
            (33, 0xffff_ffff, true),
            (34, 0xffff_ffff, true),
        ];

        let carbon = parse_smarts("[#6]", &SmartsParseParams::default())
            .expect("fixed carbon feature SMARTS parses");
        let params = MorganParams {
            radius: 0,
            ..MorganParams::default()
        };
        let call = MorganCall::default();
        let mut actual_calls = 0;

        for smiles in ["C", "CCO"] {
            let (record, valence) = fully_prepared_record_with_valence(smiles);
            let rings = fast_find_rings(&record.topology)
                .expect("fixed canonical mask target has ring information");
            let input = MorganPreparedInput {
                topology: &record.topology,
                coordinates: &record.coordinates,
                properties: &record.properties,
                valence: &valence,
                rings: &rings,
            };
            let atom_count = record.topology.atoms.len();

            for &(pattern_count, expected_mask, should_error) in &CASES {
                let patterns = (0..pattern_count)
                    .map(|_| carbon.clone())
                    .collect::<Vec<_>>();
                let mut output = AdditionalOutput {
                    atom_counts: Some(Vec::new()),
                    atom_to_bits: Some(Vec::new()),
                    bit_info_map: Some(BTreeMap::new()),
                    bit_paths: Some(BTreeMap::new()),
                    atoms_per_bit: Some(BTreeMap::new()),
                };
                let topology_before = record.topology.clone();
                let coordinates_before = record.coordinates.clone();
                let properties_before = record.properties.clone();
                let valence_before = valence.clone();
                let rings_before = rings.clone();

                let result = morgan_sparse_count(
                    &input,
                    &params,
                    &call,
                    MorganAtomInvariants::FeaturePatterns(&patterns),
                    Some(&mut output),
                );
                if should_error {
                    let error = result.expect_err("ordinal32 is the first undefined shift");
                    assert!(matches!(
                        error,
                        MorganError::Fingerprint(FingerprintError::UndefinedArithmetic {
                            site: "FingerprintUtil.cpp::getFeatureInvariants (1 << i)"
                        })
                    ));
                    assert_eq!(
                        output,
                        AdditionalOutput {
                            atom_counts: Some(vec![0; atom_count]),
                            atom_to_bits: Some(vec![Vec::new(); atom_count]),
                            bit_info_map: Some(BTreeMap::new()),
                            bit_paths: Some(BTreeMap::new()),
                            atoms_per_bit: Some(BTreeMap::new()),
                        },
                        "canonical error follows AO reset and precedes accumulation"
                    );
                } else {
                    let fingerprint = result.expect("all shifts through ordinal31 are defined");
                    let expected = if expected_mask == 0 {
                        if smiles == "C" {
                            BTreeMap::from([(0, 1)])
                        } else {
                            BTreeMap::from([(0, 3)])
                        }
                    } else if smiles == "C" {
                        BTreeMap::from([(u64::from(expected_mask), 1)])
                    } else {
                        BTreeMap::from([(0, 1), (u64::from(expected_mask), 2)])
                    };
                    assert_eq!(fingerprint.nonzero_elements(), &expected);
                    assert!(output.atom_counts.is_some());
                    assert!(output.atom_to_bits.is_some());
                    assert!(output.bit_info_map.is_some());
                    assert!(output.bit_paths.is_some());
                    assert!(output.atoms_per_bit.is_some());
                }

                assert_eq!(record.topology, topology_before);
                assert_eq!(record.coordinates, coordinates_before);
                assert_eq!(record.properties, properties_before);
                assert_eq!(valence, valence_before);
                assert_eq!(rings, rings_before);
                actual_calls += 1;
            }
        }

        assert_eq!(
            actual_calls, 14,
            "the frozen canonical product is 14 real calls"
        );
    }

    #[test]
    fn fingerprint_morgan_w01_legacy_projection_complete_product_and_bit_info() {
        const RADII: [u32; 4] = [0, 1, 2, 3];
        const BENZENE_COUNTS: [&[(u32, i32)]; 4] = [
            &[(3_218_693_969, 6)],
            &[(98_513_984, 6), (3_218_693_969, 6)],
            &[(98_513_984, 6), (2_763_854_213, 6), (3_218_693_969, 6)],
            &[
                (98_513_984, 6),
                (2_763_854_213, 6),
                (3_218_693_969, 6),
                (3_741_631_696, 1),
            ],
        ];
        const BENZENE_BITS: [&[(u32, i32)]; 4] = [
            &[(3_218_693_969, 1)],
            &[(98_513_984, 1), (3_218_693_969, 1)],
            &[(98_513_984, 1), (2_763_854_213, 1), (3_218_693_969, 1)],
            &[
                (98_513_984, 1),
                (2_763_854_213, 1),
                (3_218_693_969, 1),
                (3_741_631_696, 1),
            ],
        ];
        const MAX_INDEX: [u32; 1] = [u32::MAX];

        #[derive(Clone, Copy)]
        enum InvariantProfile {
            Connectivity,
            Features,
            CustomMaximum,
        }

        let (benzene, benzene_valence) = fully_prepared_record_with_valence("c1ccccc1");
        let benzene_rings = fast_find_rings(&benzene.topology).expect("benzene rings");
        let (carbon, carbon_valence) = fully_prepared_record_with_valence("C");
        let carbon_rings = fast_find_rings(&carbon.topology).expect("carbon ring state");
        let feature_context =
            build_prepared_query_match_context(&carbon.topology, &carbon_rings, &carbon_valence)
                .expect("the isolated carbon has a prepared feature-search context");
        let feature_target = SearchTarget::new(
            &carbon.topology,
            &carbon.coordinates,
            &carbon.topology.stereo_groups,
            Some(&carbon_rings),
            Some(&carbon_valence),
        );
        let feature_invariants = MorganFeatureAtomInvGenerator::default()
            .get_atom_invariants(&feature_target, &feature_context)
            .expect("the six pinned feature patterns accept the carbon target");
        assert_eq!(feature_invariants, [0]);

        let profiles = [
            (
                "connectivity",
                InvariantProfile::Connectivity,
                &benzene,
                &benzene_valence,
                &benzene_rings,
                None,
            ),
            (
                "features",
                InvariantProfile::Features,
                &carbon,
                &carbon_valence,
                &carbon_rings,
                Some(feature_invariants.as_slice()),
            ),
            (
                "custom",
                InvariantProfile::CustomMaximum,
                &carbon,
                &carbon_valence,
                &carbon_rings,
                Some(&MAX_INDEX[..]),
            ),
        ];

        let mut calls = 0usize;
        for use_counts in [false, true] {
            for radius in RADII {
                for (name, profile, record, valence, rings, custom) in profiles {
                    let actual = get_legacy_morgan_fingerprint(
                        &record.topology,
                        &record.properties,
                        valence,
                        rings,
                        radius,
                        false,
                        true,
                        custom,
                        None,
                        use_counts,
                        false,
                        None,
                        false,
                    )
                    .unwrap_or_else(|error| {
                        panic!(
                            "fixed W01 {name}, radius={radius}, useCounts={use_counts}: {error:?}"
                        )
                    });
                    assert_eq!(actual.length(), u32::MAX);
                    let rows = actual
                        .nonzero_elements()
                        .iter()
                        .map(|(&bit_id, &count)| (bit_id, count))
                        .collect::<Vec<_>>();
                    let expected = match profile {
                        InvariantProfile::Connectivity => {
                            if use_counts {
                                BENZENE_COUNTS[radius as usize]
                            } else {
                                BENZENE_BITS[radius as usize]
                            }
                        }
                        InvariantProfile::Features => &[(0, 1)],
                        InvariantProfile::CustomMaximum => {
                            if use_counts {
                                &[(u32::MAX, 1)]
                            } else {
                                &[(0, 1)]
                            }
                        }
                    };
                    assert_eq!(
                        rows, expected,
                        "{name}, radius={radius}, useCounts={use_counts}"
                    );
                    calls += 1;
                }
            }
        }
        assert_eq!(
            calls, 24,
            "2 useCounts values x 4 radii x 3 invariant profiles"
        );

        let expected_bit_info = BTreeMap::from([(
            3_218_693_969,
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
        )]);
        let mut bit_info_calls = 0usize;
        for use_counts in [false, true] {
            let mut bit_info = BTreeMap::from([(7, vec![(99, 99)])]);
            let actual = get_legacy_morgan_fingerprint(
                &benzene.topology,
                &benzene.properties,
                &benzene_valence,
                &benzene_rings,
                0,
                false,
                true,
                None,
                None,
                use_counts,
                false,
                Some(&mut bit_info),
                false,
            )
            .expect("fixed W01 bitInfo projection succeeds");
            let expected_count = if use_counts { 6 } else { 1 };
            assert_eq!(
                actual.nonzero_elements(),
                &BTreeMap::from([(3_218_693_969, expected_count)])
            );
            assert_eq!(
                bit_info, expected_bit_info,
                "source center order, useCounts={use_counts}"
            );
            bit_info_calls += 1;
        }
        assert_eq!(bit_info_calls, 2);

        let mut unchanged_bit_info = BTreeMap::from([(7, vec![(99, 99)])]);
        let error = get_legacy_morgan_fingerprint(
            &benzene.topology,
            &benzene.properties,
            &benzene_valence,
            &benzene_rings,
            0,
            false,
            true,
            Some(&[11]),
            None,
            true,
            false,
            Some(&mut unchanged_bit_info),
            false,
        )
        .expect_err("short source atom invariants retain the source precondition");
        assert!(matches!(
            error,
            MorganError::Fingerprint(FingerprintError::PreconditionViolation {
                what: "bad atom invariants size"
            })
        ));
        assert_eq!(unchanged_bit_info, BTreeMap::from([(7, vec![(99, 99)])]));
    }

    #[test]
    fn fingerprint_morgan_w02_hashed_counts_sizes_collisions_and_bit_info() {
        const N_BITS: [u32; 3] = [1, 128, 2048];
        const RADII: [u32; 4] = [0, 1, 2, 3];
        const EXPECTED: [[&[(u32, i32)]; 4]; 3] = [
            [&[(0, 3)], &[(0, 6)], &[(0, 6)], &[(0, 6)]],
            [
                &[(33, 1), (39, 1), (80, 1)],
                &[(2, 1), (33, 1), (38, 1), (39, 1), (80, 1), (94, 1)],
                &[(2, 1), (33, 1), (38, 1), (39, 1), (80, 1), (94, 1)],
                &[(2, 1), (33, 1), (38, 1), (39, 1), (80, 1), (94, 1)],
            ],
            [
                &[(80, 1), (807, 1), (1057, 1)],
                &[(80, 1), (222, 1), (294, 1), (807, 1), (1057, 1), (1410, 1)],
                &[(80, 1), (222, 1), (294, 1), (807, 1), (1057, 1), (1410, 1)],
                &[(80, 1), (222, 1), (294, 1), (807, 1), (1057, 1), (1410, 1)],
            ],
        ];
        const CENTER_ZERO: [u32; 1] = [0];

        let (cco, cco_valence) = fully_prepared_record_with_valence("CCO");
        let cco_rings = fast_find_rings(&cco.topology).expect("CCO ring state");
        let mut projection_calls = 0usize;
        for (size_index, n_bits) in N_BITS.into_iter().enumerate() {
            for radius in RADII {
                let actual = get_legacy_morgan_hashed_fingerprint(
                    &cco.topology,
                    &cco.properties,
                    &cco_valence,
                    &cco_rings,
                    radius,
                    n_bits,
                    None,
                    None,
                    false,
                    true,
                    false,
                    None,
                    false,
                )
                .unwrap_or_else(|error| {
                    panic!("fixed W02 CCO nBits={n_bits}, radius={radius}: {error:?}")
                });
                assert_eq!(actual.length(), n_bits);
                let rows = actual
                    .nonzero_elements()
                    .iter()
                    .map(|(&bit_id, &count)| (bit_id, count))
                    .collect::<Vec<_>>();
                assert_eq!(
                    rows, EXPECTED[size_index][radius as usize],
                    "CCO nBits={n_bits}, radius={radius}"
                );
                projection_calls += 1;
            }
        }
        assert_eq!(projection_calls, 12, "3 sizes x 4 radii");

        let (benzene, benzene_valence) = fully_prepared_record_with_valence("c1ccccc1");
        let benzene_rings = fast_find_rings(&benzene.topology).expect("benzene ring state");
        let mut bit_info_calls = 0usize;
        for n_bits in N_BITS {
            let mut bit_info = BTreeMap::from([(7, vec![(99, 99)])]);
            let actual = get_legacy_morgan_hashed_fingerprint(
                &benzene.topology,
                &benzene.properties,
                &benzene_valence,
                &benzene_rings,
                1,
                n_bits,
                None,
                Some(&CENTER_ZERO),
                false,
                true,
                false,
                Some(&mut bit_info),
                false,
            )
            .expect("fixed W02 selected-center bitInfo projection succeeds");

            let (expected_rows, expected_bit_info) = match n_bits {
                1 => (&[(0, 2)][..], BTreeMap::from([(0, vec![(0, 0), (0, 1)])])),
                128 => (
                    &[(64, 1), (81, 1)][..],
                    BTreeMap::from([(64, vec![(0, 1)]), (81, vec![(0, 0)])]),
                ),
                2048 => (
                    &[(1088, 1), (1873, 1)][..],
                    BTreeMap::from([(1088, vec![(0, 1)]), (1873, vec![(0, 0)])]),
                ),
                _ => unreachable!("the fixed W02 sizes are enumerated above"),
            };
            assert_eq!(actual.length(), n_bits);
            assert_eq!(
                actual.nonzero_elements(),
                &expected_rows.iter().copied().collect(),
                "benzene selected center nBits={n_bits}"
            );
            assert_eq!(bit_info, expected_bit_info, "nBits={n_bits}");
            bit_info_calls += 1;
        }
        assert_eq!(bit_info_calls, 3, "three positive sizes with bitInfo");

        let mut unchanged_bit_info = BTreeMap::from([(7, vec![(99, 99)])]);
        let mut zero_size_calls = 0usize;
        let error = get_legacy_morgan_hashed_fingerprint(
            &cco.topology,
            &cco.properties,
            &cco_valence,
            &cco_rings,
            0,
            0,
            Some(&[11]),
            None,
            false,
            true,
            false,
            Some(&mut unchanged_bit_info),
            false,
        )
        .expect_err("zero nBits rejects before custom-invariant validation");
        zero_size_calls += 1;
        assert!(matches!(
            error,
            MorganError::Fingerprint(FingerprintError::InvalidArguments {
                reason: "nBits can not be zero"
            })
        ));
        assert_eq!(unchanged_bit_info, BTreeMap::from([(7, vec![(99, 99)])]));
        assert_eq!(zero_size_calls, 1, "one zero-size precedence call");
    }

    #[test]
    fn fingerprint_morgan_w03_bit_vector_profiles_options_bit_info_and_errors() {
        // Literal raw environment IDs are the pinned M09 CCO rows. The only
        // expected-output transform below is the source modulo fold and
        // idempotent dense set projection.
        const CCO_SOURCE_IDS: [(u32, bool, bool, &[u64]); 16] = [
            (0, false, false, &[864662311, 2245384272, 2246728737]),
            (0, false, true, &[864662311, 2245384272, 2246728737]),
            (0, true, false, &[864662311, 2245384272, 2246728737]),
            (0, true, true, &[864662311, 2245384272, 2246728737]),
            (
                1,
                false,
                false,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                1,
                false,
                true,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                1,
                true,
                false,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                1,
                true,
                true,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                2,
                false,
                false,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                2,
                false,
                true,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                2,
                true,
                false,
                &[
                    407009239, 864662311, 1535166686, 2245384272, 2246728737, 2679480461,
                    3542456614, 3732090711, 4018048386,
                ],
            ),
            (
                2,
                true,
                true,
                &[
                    407009239, 864662311, 1535166686, 2245384272, 2246728737, 2679480461,
                    3542456614, 3732090711, 4018048386,
                ],
            ),
            (
                3,
                false,
                false,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                3,
                false,
                true,
                &[
                    864662311, 1535166686, 2245384272, 2246728737, 3542456614, 4018048386,
                ],
            ),
            (
                3,
                true,
                false,
                &[
                    407009239, 589800426, 864662311, 1535166686, 2245384272, 2246728737,
                    2679480461, 2708961779, 2962971722, 3542456614, 3732090711, 4018048386,
                ],
            ),
            (
                3,
                true,
                true,
                &[
                    407009239, 589800426, 864662311, 1535166686, 2245384272, 2246728737,
                    2679480461, 2708961779, 2962971722, 3542456614, 3732090711, 4018048386,
                ],
            ),
        ];
        const BIT_SIZES: [u32; 3] = [1, 128, 2048];

        let (cco, cco_valence) = fully_prepared_record_with_valence("CCO");
        let cco_rings = fast_find_rings(&cco.topology).expect("CCO rings");
        let mut profile_calls = 0usize;
        let mut total_calls = 0usize;
        for (radius, include_redundant, only_nonzero, source_ids) in CCO_SOURCE_IDS {
            for n_bits in BIT_SIZES {
                let actual = get_legacy_morgan_fingerprint_as_bit_vector(
                    &cco.topology,
                    &cco.properties,
                    &cco_valence,
                    &cco_rings,
                    radius,
                    n_bits,
                    None,
                    None,
                    false,
                    true,
                    only_nonzero,
                    None,
                    include_redundant,
                )
                .unwrap_or_else(|error| {
                    panic!(
                        "fixed W03 CCO nBits={n_bits}, radius={radius}, redundant={include_redundant}, \
                         onlyNonzero={only_nonzero}: {error:?}"
                    )
                });
                let mut expected_bits = source_ids
                    .iter()
                    .map(|&source_id| (source_id % u64::from(n_bits)) as u32)
                    .collect::<Vec<_>>();
                expected_bits.sort_unstable();
                expected_bits.dedup();
                assert_eq!(actual.n_bits(), n_bits);
                assert_eq!(
                    actual.on_bits(),
                    expected_bits,
                    "CCO nBits={n_bits}, radius={radius}, redundant={include_redundant}, onlyNonzero={only_nonzero}"
                );
                profile_calls += 1;
                total_calls += 1;
            }
        }
        assert_eq!(profile_calls, 48, "16 literal profiles x 3 positive sizes");

        // Full Cartesian product of source booleans and positive sizes, on a
        // source-minimal zero-custom-invariant case.
        let (isolated_carbon, carbon_valence) = fully_prepared_record_with_valence("C");
        let carbon_rings =
            fast_find_rings(&isolated_carbon.topology).expect("isolated carbon rings");
        const ZERO_INVARIANT: [u32; 1] = [0];
        let mut option_product_calls = 0usize;
        for use_chirality in [false, true] {
            for use_bond_types in [false, true] {
                for only_nonzero in [false, true] {
                    for include_redundant in [false, true] {
                        for n_bits in BIT_SIZES {
                            let actual = get_legacy_morgan_fingerprint_as_bit_vector(
                                &isolated_carbon.topology,
                                &isolated_carbon.properties,
                                &carbon_valence,
                                &carbon_rings,
                                0,
                                n_bits,
                                Some(&ZERO_INVARIANT),
                                None,
                                use_chirality,
                                use_bond_types,
                                only_nonzero,
                                None,
                                include_redundant,
                            )
                            .expect("fixed W03 zero-invariant option product succeeds");
                            let expected_bits = if only_nonzero { vec![] } else { vec![0] };
                            assert_eq!(
                                actual.n_bits(),
                                n_bits,
                                "nBits={n_bits}, chirality={use_chirality}, bondTypes={use_bond_types}, onlyNonzero={only_nonzero}, redundant={include_redundant}"
                            );
                            assert_eq!(
                                actual.on_bits(),
                                expected_bits,
                                "nBits={n_bits}, chirality={use_chirality}, bondTypes={use_bond_types}, onlyNonzero={only_nonzero}, redundant={include_redundant}"
                            );
                            option_product_calls += 1;
                            total_calls += 1;
                        }
                    }
                }
            }
        }
        assert_eq!(
            option_product_calls, 48,
            "2^4 source booleans x 3 positive sizes"
        );

        // Pinned test1.cpp:1219-1279 checks these source relations for the
        // legacy wrapper, including custom atom invariants and chirality.
        let chiral_records = [
            fully_prepared_record_with_valence("C[C@H](F)Cl"),
            fully_prepared_record_with_valence("C[C@@H](F)Cl"),
            fully_prepared_record_with_valence("CC(F)Cl"),
        ];
        let chiral_rings = chiral_records
            .iter()
            .map(|(record, _)| fast_find_rings(&record.topology).expect("fixed chiral rings"))
            .collect::<Vec<_>>();
        let mut chirality_calls = 0usize;
        for use_chirality in [false, true] {
            let mut fingerprints = Vec::with_capacity(chiral_records.len());
            for ((record, valence), rings) in chiral_records.iter().zip(&chiral_rings) {
                fingerprints.push(
                    get_legacy_morgan_fingerprint_as_bit_vector(
                        &record.topology,
                        &record.properties,
                        valence,
                        rings,
                        2,
                        2048,
                        None,
                        None,
                        use_chirality,
                        true,
                        false,
                        None,
                        false,
                    )
                    .expect("fixed W03 source chirality profile succeeds"),
                );
                chirality_calls += 1;
                total_calls += 1;
            }
            if use_chirality {
                assert_ne!(fingerprints[0], fingerprints[1]);
                assert_ne!(fingerprints[0], fingerprints[2]);
                assert_ne!(fingerprints[1], fingerprints[2]);
            } else {
                assert_eq!(fingerprints[0], fingerprints[1]);
                assert_eq!(fingerprints[0], fingerprints[2]);
                assert_eq!(fingerprints[1], fingerprints[2]);
            }
        }
        assert_eq!(chirality_calls, 6, "three pinned structures x two flags");

        // Pinned test1.cpp:1219-1240 checks this custom-invariant bond-type
        // relation: CCC and CC=C differ only when bond types are enabled.
        const CUSTOM_INVARIANTS: [u32; 3] = [1, 1, 1];
        let (ccc, ccc_valence) = fully_prepared_record_with_valence("CCC");
        let ccc_rings = fast_find_rings(&ccc.topology).expect("CCC rings");
        let (cc_double, cc_double_valence) = fully_prepared_record_with_valence("CC=C");
        let cc_double_rings = fast_find_rings(&cc_double.topology).expect("CC=C rings");
        let mut bond_type_calls = 0usize;
        for use_bond_types in [true, false] {
            let single = get_legacy_morgan_fingerprint_as_bit_vector(
                &ccc.topology,
                &ccc.properties,
                &ccc_valence,
                &ccc_rings,
                2,
                2048,
                Some(&CUSTOM_INVARIANTS),
                None,
                false,
                use_bond_types,
                false,
                None,
                false,
            )
            .expect("fixed W03 single-bond custom-invariant output succeeds");
            let double = get_legacy_morgan_fingerprint_as_bit_vector(
                &cc_double.topology,
                &cc_double.properties,
                &cc_double_valence,
                &cc_double_rings,
                2,
                2048,
                Some(&CUSTOM_INVARIANTS),
                None,
                false,
                use_bond_types,
                false,
                None,
                false,
            )
            .expect("fixed W03 double-bond custom-invariant output succeeds");
            if use_bond_types {
                assert_ne!(single, double);
            } else {
                assert_eq!(single, double);
            }
            bond_type_calls += 2;
            total_calls += 2;
        }
        assert_eq!(bond_type_calls, 4, "two molecules x two bond-type flags");

        // Selected benzene center zero emits the two M09 source codes at
        // radii zero and one; size1 folds them onto one bit but keeps both
        // bitInfo provenance entries in source radius order.
        const CENTER_ZERO: [u32; 1] = [0];
        let (benzene, benzene_valence) = fully_prepared_record_with_valence("c1ccccc1");
        let benzene_rings = fast_find_rings(&benzene.topology).expect("benzene rings");
        let mut bit_info_calls = 0usize;
        for n_bits in BIT_SIZES {
            let mut bit_info = BTreeMap::from([(7, vec![(99, 99)])]);
            let actual = get_legacy_morgan_fingerprint_as_bit_vector(
                &benzene.topology,
                &benzene.properties,
                &benzene_valence,
                &benzene_rings,
                1,
                n_bits,
                None,
                Some(&CENTER_ZERO),
                false,
                true,
                false,
                Some(&mut bit_info),
                false,
            )
            .expect("fixed W03 selected-center dense bitInfo succeeds");
            let (expected_bits, expected_info) = match n_bits {
                1 => (&[0][..], BTreeMap::from([(0, vec![(0, 0), (0, 1)])])),
                128 => (
                    &[64, 81][..],
                    BTreeMap::from([(64, vec![(0, 1)]), (81, vec![(0, 0)])]),
                ),
                2048 => (
                    &[1088, 1873][..],
                    BTreeMap::from([(1088, vec![(0, 1)]), (1873, vec![(0, 0)])]),
                ),
                _ => unreachable!("fixed W03 bit sizes are enumerated above"),
            };
            assert_eq!(actual.n_bits(), n_bits);
            assert_eq!(actual.on_bits().as_slice(), expected_bits);
            assert_eq!(bit_info, expected_info, "nBits={n_bits}");
            bit_info_calls += 1;
            total_calls += 1;
        }
        assert_eq!(bit_info_calls, 3, "selected-center bitInfo at all sizes");

        // Zero size wins over invalid fromAtoms/invariants and leaves the
        // caller's existing map untouched.
        const INVALID_CENTER: [u32; 1] = [u32::MAX];
        const SHORT_INVARIANTS: [u32; 1] = [11];
        let mut unchanged_bit_info = BTreeMap::from([(7, vec![(99, 99)])]);
        let zero_error = get_legacy_morgan_fingerprint_as_bit_vector(
            &cco.topology,
            &cco.properties,
            &cco_valence,
            &cco_rings,
            0,
            0,
            Some(&SHORT_INVARIANTS),
            Some(&INVALID_CENTER),
            false,
            true,
            false,
            Some(&mut unchanged_bit_info),
            false,
        )
        .expect_err("zero nBits rejects before invalid source options");
        assert!(matches!(
            zero_error,
            MorganError::Fingerprint(FingerprintError::InvalidArguments {
                reason: "nBits can not be zero"
            })
        ));
        assert_eq!(unchanged_bit_info, BTreeMap::from([(7, vec![(99, 99)])]));
        total_calls += 1;

        let short_invariant_error = get_legacy_morgan_fingerprint_as_bit_vector(
            &cco.topology,
            &cco.properties,
            &cco_valence,
            &cco_rings,
            0,
            128,
            Some(&SHORT_INVARIANTS),
            None,
            false,
            true,
            false,
            Some(&mut unchanged_bit_info),
            false,
        )
        .expect_err("short source custom atom invariants fail");
        assert!(matches!(
            short_invariant_error,
            MorganError::Fingerprint(FingerprintError::PreconditionViolation {
                what: "bad atom invariants size"
            })
        ));
        assert_eq!(unchanged_bit_info, BTreeMap::from([(7, vec![(99, 99)])]));
        total_calls += 1;

        assert_eq!(total_calls, 111, "all W03 projection and option calls");
    }

    #[derive(Clone, Copy)]
    enum ExpectedRadiusZero {
        Rows(&'static [(u32, u32)]),
        SourceAssertion,
    }

    #[test]
    fn fingerprint_morgan_m02_get_bit_id_returns_stored_code_for_all_arguments() {
        const CODES: [(u32, u32, u64); 4] = [
            (0x0000_0000, 0x0000_0000, 0x0000_0000_0000_0000),
            (0x8000_0000, 0x8000_0000, 0x0000_0000_8000_0000),
            (0xFEDC_BA98, 0xFEDC_BA98, 0x0000_0000_FEDC_BA98),
            (0xFFFF_FFFF, 0xFFFF_FFFF, 0x0000_0000_FFFF_FFFF),
        ];
        let arguments = FingerprintArguments::new(false, vec![1, 2, 4, 8], 2048, 1, false)
            .expect("the fixed source argument values satisfy the constructor");
        let topology = TopologyBlock::default();
        let atom_invariants = [0x1234_5678, 0x9ABC_DEF0];
        let bond_invariants = [0x0BAD_F00D, 0xCAFE_BABE];
        let mut tuple_count = 0usize;
        let mut output_width_count = 0usize;

        for (code, expected_u32, expected_u64) in CODES {
            let environment_u32 = MorganAtomEnvironment::<u32>::new(code, 0, 0, &topology);
            let environment_u64 = MorganAtomEnvironment::<u64>::new(code, 0, 0, &topology);
            for hash_results in [false, true] {
                for fp_size in [0, 2048] {
                    for optional_mask in 0u8..16 {
                        let arguments_input = (optional_mask & 1 != 0).then_some(&arguments);
                        let atom_input = (optional_mask & 2 != 0).then_some(&atom_invariants[..]);
                        let bond_input = (optional_mask & 4 != 0).then_some(&bond_invariants[..]);
                        let mut additional_output_u32 = AdditionalOutput::default();
                        let actual_u32 = environment_u32.get_bit_id(
                            arguments_input,
                            atom_input,
                            bond_input,
                            if optional_mask & 8 != 0 {
                                Some(&mut additional_output_u32)
                            } else {
                                None
                            },
                            hash_results,
                            fp_size,
                        );
                        assert_eq!(actual_u32, expected_u32);
                        assert_eq!(additional_output_u32, AdditionalOutput::default());

                        let mut additional_output_u64 = AdditionalOutput::default();
                        let actual_u64 = environment_u64.get_bit_id(
                            arguments_input,
                            atom_input,
                            bond_input,
                            if optional_mask & 8 != 0 {
                                Some(&mut additional_output_u64)
                            } else {
                                None
                            },
                            hash_results,
                            fp_size,
                        );
                        assert_eq!(actual_u64, expected_u64);
                        assert_eq!(additional_output_u64, AdditionalOutput::default());

                        tuple_count += 1;
                        output_width_count += 2;
                    }
                }
            }
        }

        assert_eq!(tuple_count, 256);
        assert_eq!(output_width_count, 512);
    }

    #[test]
    fn fingerprint_morgan_m03_updates_all_optional_outputs_in_source_order() {
        let record =
            cosmolkit_smiles::parse_smiles("CCC", &cosmolkit_smiles::SmilesParseParams::default())
                .expect("the fixed three-carbon chain is valid");
        let topology = record.topology;
        assert_eq!(topology.atoms.len(), 3);
        assert_eq!(topology.bonds.len(), 2);

        let environments = [
            MorganAtomEnvironment::<u32>::new(11, 1, 0, &topology),
            MorganAtomEnvironment::<u32>::new(12, 0, 1, &topology),
            MorganAtomEnvironment::<u32>::new(13, 1, 1, &topology),
        ];
        let mut distance_matrix_cache = MorganDistanceMatrixCache::default();

        for mask in 0u8..32 {
            let mut output = AdditionalOutput {
                atom_counts: (mask & 0b00001 != 0).then(|| vec![0, 0, 0]),
                atom_to_bits: (mask & 0b00010 != 0).then(|| vec![vec![], vec![], vec![]]),
                bit_info_map: (mask & 0b00100 != 0).then(BTreeMap::new),
                bit_paths: (mask & 0b01000 != 0).then(|| BTreeMap::from([(999, vec![vec![42]])])),
                atoms_per_bit: (mask & 0b10000 != 0).then(BTreeMap::new),
            };

            for environment in &environments {
                environment
                    .update_additional_output(&mut output, 77, &mut distance_matrix_cache)
                    .expect("the valid fixed topology has a distance matrix");
            }

            let expected = AdditionalOutput {
                atom_counts: (mask & 0b00001 != 0).then(|| vec![1, 2, 0]),
                atom_to_bits: (mask & 0b00010 != 0).then(|| vec![vec![77], vec![77, 77], vec![]]),
                bit_info_map: (mask & 0b00100 != 0)
                    .then(|| BTreeMap::from([(77, vec![(1, 0), (0, 1), (1, 1)])])),
                bit_paths: (mask & 0b01000 != 0).then(|| BTreeMap::from([(999, vec![vec![42]])])),
                atoms_per_bit: (mask & 0b10000 != 0)
                    .then(|| BTreeMap::from([(77, vec![vec![1], vec![0, 1], vec![1, 0, 2]])])),
            };
            assert_eq!(output, expected, "optional-output mask {mask:#07b}");
        }
        assert_eq!(distance_matrix_cache.matrix.is_some(), true);
    }

    #[test]
    fn fingerprint_morgan_m03_atom_count_wraps_at_source_unsigned_width() {
        let record =
            cosmolkit_smiles::parse_smiles("C", &cosmolkit_smiles::SmilesParseParams::default())
                .expect("the fixed carbon atom is valid");
        let topology = record.topology;
        let environment = MorganAtomEnvironment::<u32>::new(0, 0, 0, &topology);
        let mut output = AdditionalOutput {
            atom_counts: Some(vec![u32::MAX]),
            ..AdditionalOutput::default()
        };
        let mut distance_matrix_cache = MorganDistanceMatrixCache::default();

        environment
            .update_additional_output(&mut output, 77, &mut distance_matrix_cache)
            .expect("layer zero does not require a distance matrix");
        assert_eq!(output.atom_counts, Some(vec![0]));
        assert!(distance_matrix_cache.matrix.is_none());
    }

    #[test]
    fn fingerprint_morgan_m04_atoms_per_bit_matches_line_ring_and_disconnected_rows() {
        let cases: [(&str, u32, usize, usize, [&[i32]; 3]); 3] = [
            ("CCCC", 1, 4, 3, [&[1], &[1, 0, 2], &[1, 0, 2, 3]]),
            ("C1CCCCC1", 0, 6, 6, [&[0], &[0, 1, 5], &[0, 1, 2, 4, 5]]),
            ("CC.O", 0, 3, 1, [&[0], &[0, 1], &[0, 1]]),
        ];

        for (smiles, center, atom_count, bond_count, expected_by_radius) in cases {
            let record = cosmolkit_smiles::parse_smiles(
                smiles,
                &cosmolkit_smiles::SmilesParseParams::default(),
            )
            .expect("each fixed topology is valid");
            let topology = record.topology;
            assert_eq!(topology.atoms.len(), atom_count, "{smiles}");
            assert_eq!(topology.bonds.len(), bond_count, "{smiles}");

            let mut distance_matrix_cache = MorganDistanceMatrixCache::default();
            for (layer, expected_atoms) in expected_by_radius.into_iter().enumerate() {
                let environment = MorganAtomEnvironment::<u32>::new(
                    101 + layer as u32,
                    center,
                    layer as u32,
                    &topology,
                );
                let mut output = AdditionalOutput {
                    atoms_per_bit: Some(BTreeMap::new()),
                    ..AdditionalOutput::default()
                };
                environment
                    .update_additional_output(&mut output, 91, &mut distance_matrix_cache)
                    .expect("the valid fixed topology has a distance matrix");

                assert_eq!(
                    output.atoms_per_bit,
                    Some(BTreeMap::from([(91, vec![expected_atoms.to_vec()])])),
                    "{smiles}, center {center}, radius {layer}"
                );
                assert_eq!(distance_matrix_cache.matrix.is_some(), layer > 0);
            }
        }
    }

    #[test]
    fn fingerprint_morgan_m05_radius_zero_center_selection_source_product() {
        use ExpectedRadiusZero::{Rows, SourceAssertion};

        const EMPTY: [[ExpectedRadiusZero; 5]; 2] = [
            [
                Rows(&[]),
                Rows(&[]),
                SourceAssertion,
                SourceAssertion,
                Rows(&[]),
            ],
            [
                Rows(&[]),
                Rows(&[]),
                SourceAssertion,
                SourceAssertion,
                Rows(&[]),
            ],
        ];
        const ETHANOL: [[ExpectedRadiusZero; 5]; 2] = [
            [
                Rows(&[(0, 10), (1, 0), (2, 30)]),
                Rows(&[]),
                Rows(&[(2, 30)]),
                Rows(&[(0, 10), (2, 30)]),
                Rows(&[(0, 10), (1, 0), (2, 30)]),
            ],
            [
                Rows(&[(0, 10), (2, 30)]),
                Rows(&[]),
                Rows(&[(2, 30)]),
                Rows(&[(0, 10), (2, 30)]),
                Rows(&[(0, 10), (2, 30)]),
            ],
        ];
        const PYRROLE: [[ExpectedRadiusZero; 5]; 2] = [
            [
                Rows(&[(0, 0), (1, 21), (2, 22), (3, 0), (4, 25)]),
                Rows(&[]),
                Rows(&[(4, 25)]),
                Rows(&[(0, 0), (4, 25)]),
                Rows(&[(0, 0), (1, 21), (2, 22), (3, 0), (4, 25)]),
            ],
            [
                Rows(&[(1, 21), (2, 22), (4, 25)]),
                Rows(&[]),
                Rows(&[(4, 25)]),
                Rows(&[(4, 25)]),
                Rows(&[(1, 21), (2, 22), (4, 25)]),
            ],
        ];
        const BENZENE: [[ExpectedRadiusZero; 5]; 2] = [
            [
                Rows(&[(0, 31), (1, 0), (2, 33), (3, 34), (4, 0), (5, 36)]),
                Rows(&[]),
                Rows(&[(5, 36)]),
                Rows(&[(0, 31), (5, 36)]),
                Rows(&[(0, 31), (1, 0), (2, 33), (3, 34), (4, 0), (5, 36)]),
            ],
            [
                Rows(&[(0, 31), (2, 33), (3, 34), (5, 36)]),
                Rows(&[]),
                Rows(&[(5, 36)]),
                Rows(&[(0, 31), (5, 36)]),
                Rows(&[(0, 31), (2, 33), (3, 34), (5, 36)]),
            ],
        ];
        const SALT: [[ExpectedRadiusZero; 5]; 2] = [
            [
                Rows(&[(0, 41), (1, 42)]),
                Rows(&[]),
                Rows(&[(1, 42)]),
                Rows(&[(0, 41), (1, 42)]),
                Rows(&[(0, 41), (1, 42)]),
            ],
            [
                Rows(&[(0, 41), (1, 42)]),
                Rows(&[]),
                Rows(&[(1, 42)]),
                Rows(&[(0, 41), (1, 42)]),
                Rows(&[(0, 41), (1, 42)]),
            ],
        ];

        const CASES: [(
            &str,
            usize,
            &[u32],
            &[u32],
            &[u32],
            &[u32],
            [[ExpectedRadiusZero; 5]; 2],
        ); 5] = [
            ("", 0, &[], &[0], &[0, 0], &[], EMPTY),
            (
                "CCO",
                3,
                &[10, 0, 30],
                &[2],
                &[2, 0, 2],
                &[0, 1, 2],
                ETHANOL,
            ),
            (
                "c1cc[nH]c1",
                5,
                &[0, 21, 22, 0, 25],
                &[4],
                &[4, 0, 4],
                &[0, 1, 2, 3, 4],
                PYRROLE,
            ),
            (
                "c1ccccc1",
                6,
                &[31, 0, 33, 34, 0, 36],
                &[5],
                &[5, 0, 5],
                &[0, 1, 2, 3, 4, 5],
                BENZENE,
            ),
            ("[Na+].[Cl-]", 2, &[41, 42], &[1], &[1, 0, 1], &[0, 1], SALT),
        ];
        const ONLY_NONZERO: [bool; 2] = [false, true];
        const IGNORE_ATOMS: [(&str, Option<&[u32]>); 2] = [
            ("absent", None),
            ("present_out_of_range", Some(&[u32::MAX])),
        ];

        let mut calls = 0usize;
        let mut source_assertion_cases = 0usize;
        for (smiles, atom_count, invariants, single, duplicate, all, expected_by_options) in CASES {
            let record = fully_prepared_record(smiles);
            assert_eq!(record.topology.atoms.len(), atom_count, "{smiles}");
            let selectors: [(&str, Option<&[u32]>); 5] = [
                ("none", None),
                ("empty", Some(&[])),
                ("single", Some(single)),
                ("duplicate", Some(duplicate)),
                ("all", Some(all)),
            ];

            for (only_nonzero_index, only_nonzero) in ONLY_NONZERO.into_iter().enumerate() {
                for (selector_index, (selector_name, from_atoms)) in
                    selectors.iter().copied().enumerate()
                {
                    let expected = expected_by_options[only_nonzero_index][selector_index];
                    for (ignore_name, ignore_atoms) in IGNORE_ATOMS {
                        calls += 1;
                        let actual = catch_unwind(AssertUnwindSafe(|| {
                            let include_atoms =
                                selected_atom_mask(record.topology.atoms.len(), from_atoms);
                            let mut environments = Vec::with_capacity(atom_count);
                            radius_zero_environments::<u32>(
                                &record.topology,
                                invariants,
                                &include_atoms,
                                ignore_atoms,
                                only_nonzero,
                                &mut environments,
                            );
                            environments
                        }));

                        match expected {
                            ExpectedRadiusZero::Rows(expected_rows) => {
                                let environments = actual.unwrap_or_else(|_| {
                                    panic!(
                                        "unexpected assertion for {smiles}, onlyNonzero={only_nonzero}, \
                                         fromAtoms={selector_name}, ignoreAtoms={ignore_name}"
                                    )
                                });
                                let actual_rows = environments
                                    .iter()
                                    .map(|environment| {
                                        assert_eq!(environment.d_layer, 0);
                                        (environment.d_atom_id, environment.d_code)
                                    })
                                    .collect::<Vec<_>>();
                                assert_eq!(
                                    actual_rows, expected_rows,
                                    "{smiles}, onlyNonzero={only_nonzero}, fromAtoms={selector_name}, ignoreAtoms={ignore_name}"
                                );
                            }
                            ExpectedRadiusZero::SourceAssertion => {
                                assert!(
                                    actual.is_err(),
                                    "expected the source assertion boundary for {smiles}, \
                                     onlyNonzero={only_nonzero}, fromAtoms={selector_name}, \
                                     ignoreAtoms={ignore_name}"
                                );
                                source_assertion_cases += 1;
                            }
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 100);
        assert_eq!(source_assertion_cases, 8);
    }

    fn m06_topology(atom_count: usize, bonds: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = bonds
            .iter()
            .enumerate()
            .map(|(index, &(begin, end, order))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("the fixed M06 topology rows are valid")
    }

    fn m06_remap_mask(mask: u64, old_bond_to_new: &[usize; 9]) -> u64 {
        let mut remapped = 0u64;
        for (old_bond, &new_bond) in old_bond_to_new.iter().enumerate() {
            if mask & (1u64 << old_bond) != 0 {
                remapped |= 1u64 << new_bond;
            }
        }
        remapped
    }

    #[test]
    fn fingerprint_morgan_m06_neighbor_hash_is_row_order_independent_and_source_sorted() {
        const SOURCE_BONDS: [(usize, usize, BondOrder); 9] = [
            (0, 1, BondOrder::Single),
            (0, 2, BondOrder::Double),
            (0, 3, BondOrder::Triple),
            (0, 4, BondOrder::Aromatic),
            (1, 5, BondOrder::Single),
            (2, 6, BondOrder::Double),
            (3, 7, BondOrder::Triple),
            (4, 8, BondOrder::Single),
            (9, 10, BondOrder::Single),
        ];
        const SOURCE_ATOM_INVARIANTS: [u32; 11] = [
            0xfedc_ba98,
            0x8000_0005,
            4,
            5,
            0x7fff_ffff,
            0xeeee_0005,
            0xeeee_0006,
            0xeeee_0007,
            0xeeee_0008,
            0xeeee_0009,
            0xeeee_000a,
        ];
        const SOURCE_BOND_INVARIANTS: [u32; 9] =
            [u32::MAX, 2, 0x8000_0000, 2, 101, 102, 103, 104, 105];
        const LAYER_ONE_INPUTS: [u64; 11] = [
            0x0000_000f,
            0x0000_0011,
            0x0000_0022,
            0x0000_0044,
            0x0000_0088,
            0x0000_0010,
            0x0000_0020,
            0x0000_0040,
            0x0000_0080,
            0x0000_0100,
            0x0000_0100,
        ];
        const LAYER_TWO_INPUTS: [u64; 11] = [
            0x0000_00ff,
            0x0000_001f,
            0x0000_003f,
            0x0000_007f,
            0x0000_00ff,
            0x0000_0011,
            0x0000_0022,
            0x0000_0044,
            0x0000_0088,
            0x0000_0100,
            0x0000_0100,
        ];
        const SATURATED_INPUTS: [u64; 11] = [
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_00ff,
            0x0000_0100,
            0x0000_0100,
        ];
        const IDENTITY_ATOM_ROWS: [usize; 11] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        const PERMUTED_ATOM_ROWS: [usize; 11] = [3, 0, 7, 2, 8, 4, 9, 5, 10, 1, 6];
        const IDENTITY_BOND_ROWS: [usize; 9] = [0, 1, 2, 3, 4, 5, 6, 7, 8];
        const PERMUTED_BOND_ROWS: [usize; 9] = [8, 3, 6, 1, 5, 2, 7, 0, 4];
        const LAYERS: [(u32, &[u64; 11], u32); 4] = [
            (0, &[0; 11], 982_951_719),
            (1, &LAYER_ONE_INPUTS, 677_145_243),
            (2, &LAYER_TWO_INPUTS, 3_130_291_034),
            (u32::MAX - 1, &SATURATED_INPUTS, 1_516_261_133),
        ];
        const EXPECTED_PAIR_HASHES: [u32; 4] =
            [764_723_157, 1_301_594_004, 3_449_077_584, 1_301_593_949];
        const EXPECTED_SORTED_PAIRS: [(i32, u32); 4] =
            [(i32::MIN, 5), (-1, 0x8000_0005), (2, 4), (2, 0x7fff_ffff)];

        // Independent literal pair folds from RDGeneral/hash.hpp exercise the
        // pinned signed ordering and 32-bit wrapping before the full-row cases.
        for (pair, expected_hash) in EXPECTED_SORTED_PAIRS.into_iter().zip(EXPECTED_PAIR_HASHES) {
            let mut pair_hash = 0;
            hash_combine(&mut pair_hash, hash_value_i32(pair.0));
            hash_combine(&mut pair_hash, hash_value_u32(pair.1));
            assert_eq!(pair_hash, expected_hash, "source pair {pair:?}");
        }

        let atom_orders = [IDENTITY_ATOM_ROWS, PERMUTED_ATOM_ROWS];
        let bond_orders = [IDENTITY_BOND_ROWS, PERMUTED_BOND_ROWS];
        let expected_center_masks = [
            [0x0000_000f, 0x0000_00ff, 0x0000_00ff, 0x0000_00ff],
            [0x0000_00aa, 0x0000_01fe, 0x0000_01fe, 0x0000_01fe],
        ];
        let mut calls = 0usize;

        for (atom_order_index, old_atom_to_new) in atom_orders.into_iter().enumerate() {
            let mut current_invariants = vec![0u32; SOURCE_ATOM_INVARIANTS.len()];
            for (old_atom, &new_atom) in old_atom_to_new.iter().enumerate() {
                current_invariants[new_atom] = SOURCE_ATOM_INVARIANTS[old_atom];
            }

            for (bond_order_index, new_row_to_old_bond) in bond_orders.into_iter().enumerate() {
                let mut old_bond_to_new = [0usize; 9];
                let mut bond_rows = Vec::with_capacity(SOURCE_BONDS.len());
                let mut bond_invariants = vec![0u32; SOURCE_BONDS.len()];
                for (new_bond, old_bond) in new_row_to_old_bond.into_iter().enumerate() {
                    old_bond_to_new[old_bond] = new_bond;
                    let (old_begin, old_end, order) = SOURCE_BONDS[old_bond];
                    bond_rows.push((old_atom_to_new[old_begin], old_atom_to_new[old_end], order));
                    bond_invariants[new_bond] = SOURCE_BOND_INVARIANTS[old_bond];
                }
                let topology = m06_topology(SOURCE_ATOM_INVARIANTS.len(), &bond_rows);
                let center = old_atom_to_new[0];

                for (layer_index, (layer, source_inputs, expected_hash)) in
                    LAYERS.into_iter().enumerate()
                {
                    let mut atom_neighborhoods = vec![
                        MorganBondEnvironment {
                            bit_count: SOURCE_BONDS.len(),
                            blocks: vec![0],
                        };
                        SOURCE_ATOM_INVARIANTS.len()
                    ];
                    for (old_atom, &source_mask) in source_inputs.iter().enumerate() {
                        atom_neighborhoods[old_atom_to_new[old_atom]].blocks[0] =
                            m06_remap_mask(source_mask, &old_bond_to_new);
                    }
                    let mut round_atom_neighborhoods = atom_neighborhoods.clone();
                    let mut neighborhood_invariants = Vec::with_capacity(8);
                    let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());

                    let actual_hash = update_neighbor_layer(
                        &topology,
                        center,
                        layer,
                        false,
                        &mut chiral_atoms,
                        &current_invariants,
                        &bond_invariants,
                        &atom_neighborhoods,
                        &mut round_atom_neighborhoods,
                        &mut neighborhood_invariants,
                    )
                    .expect("the fixed center has four source neighbors");

                    assert_eq!(actual_hash, expected_hash);
                    assert_eq!(neighborhood_invariants, EXPECTED_SORTED_PAIRS);
                    assert_eq!(
                        round_atom_neighborhoods[center].blocks,
                        vec![expected_center_masks[bond_order_index][layer_index]],
                        "atom row order {atom_order_index}, bond row order {bond_order_index}, layer {layer}"
                    );
                    calls += 1;
                }
            }
        }

        assert_eq!(calls, 16);
    }

    #[test]
    fn fingerprint_morgan_m06_ring_and_disconnected_rows_keep_source_bond_bits() {
        const BONDS: [(usize, usize, BondOrder); 7] = [
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 5, BondOrder::Triple),
            (5, 0, BondOrder::Single),
            (6, 7, BondOrder::Single),
        ];
        const BOND_INVARIANTS: [u32; 7] = [0x8000_0001, 20, 30, 40, 50, 5, 99];
        const CURRENT_INVARIANTS: [u32; 8] = [0x1357_9bdf, 42, 10, 11, 12, 7, 90, 91];
        const LAYER_ONE_INPUTS: [u64; 8] = [0x21, 0x03, 0x06, 0x0c, 0x18, 0x30, 0x40, 0x40];
        const LAYER_TWO_INPUTS: [u64; 8] = [0x33, 0x27, 0x0f, 0x1e, 0x3c, 0x39, 0x40, 0x40];
        const SATURATED_INPUTS: [u64; 8] = [0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x40, 0x40];
        const LAYERS: [(u32, &[u64; 8], u32, u64); 4] = [
            (0, &[0; 8], 2_944_231_836, 0x21),
            (1, &LAYER_ONE_INPUTS, 2_943_971_785, 0x33),
            (2, &LAYER_TWO_INPUTS, 2_943_695_782, 0x3f),
            (u32::MAX - 1, &SATURATED_INPUTS, 2_682_364_464, 0x3f),
        ];
        const EXPECTED_SORTED_PAIRS: [(i32, u32); 2] = [(i32::MIN + 1, 42), (5, 7)];

        let topology = m06_topology(8, &BONDS);
        let mut calls = 0usize;
        for (layer, source_inputs, expected_hash, expected_mask) in LAYERS {
            let atom_neighborhoods = source_inputs
                .iter()
                .map(|&mask| MorganBondEnvironment {
                    bit_count: BONDS.len(),
                    blocks: vec![mask],
                })
                .collect::<Vec<_>>();
            let mut round_atom_neighborhoods = atom_neighborhoods.clone();
            let mut neighborhood_invariants = Vec::with_capacity(8);
            let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());

            let actual_hash = update_neighbor_layer(
                &topology,
                0,
                layer,
                false,
                &mut chiral_atoms,
                &CURRENT_INVARIANTS,
                &BOND_INVARIANTS,
                &atom_neighborhoods,
                &mut round_atom_neighborhoods,
                &mut neighborhood_invariants,
            )
            .expect("the fixed ring center has two source neighbors");

            assert_eq!(actual_hash, expected_hash);
            assert_eq!(neighborhood_invariants, EXPECTED_SORTED_PAIRS);
            assert_eq!(round_atom_neighborhoods[0].blocks, vec![expected_mask]);
            assert_eq!(round_atom_neighborhoods[0].bit_count, 7);
            calls += 1;
        }
        assert_eq!(calls, 4);
    }

    #[test]
    fn fingerprint_morgan_m06_isolated_atom_takes_source_skip_before_scratch_clear() {
        let topology = m06_topology(1, &[]);
        let neighborhoods = [MorganBondEnvironment::new(0)];
        let mut round_neighborhoods = neighborhoods.clone();
        let mut neighborhood_invariants = vec![(17, 19)];
        let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());

        let actual = update_neighbor_layer(
            &topology,
            0,
            0,
            false,
            &mut chiral_atoms,
            &[23],
            &[],
            &neighborhoods,
            &mut round_neighborhoods,
            &mut neighborhood_invariants,
        );

        assert_eq!(actual, None);
        assert_eq!(round_neighborhoods, neighborhoods);
        assert_eq!(neighborhood_invariants, vec![(17, 19)]);
    }

    fn m07_topology(chiral_tag: ChiralTag, cip_code: Option<PropertyValue>) -> TopologyBlock {
        let mut center_spec = AtomSpec::new(Element::C).with_chiral_tag(chiral_tag);
        if let Some(cip_code) = cip_code {
            center_spec = center_spec
                .with_prop("_CIPCode", cip_code)
                .expect("the fixed CIP property key is valid");
        }

        let mut atoms = vec![Atom::from_spec(AtomId::new(0), center_spec)];
        atoms.extend(
            (1..5).map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C))),
        );
        let bonds = (1..5)
            .enumerate()
            .map(|(bond_index, atom_index)| {
                Bond::from_spec(
                    BondId::new(bond_index),
                    BondSpec::new(AtomId::new(0), AtomId::new(atom_index), BondOrder::Single),
                )
            })
            .collect();

        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("the fixed M07 four-neighbor topology is valid")
    }

    fn m07_layer_code(
        topology: &TopologyBlock,
        layer: u32,
        include_chirality: bool,
        current_invariants: &[u32],
        bond_invariants: &[u32],
        chiral_atoms: &mut MorganChiralAtoms,
    ) -> u32 {
        let atom_neighborhoods =
            vec![MorganBondEnvironment::new(topology.bonds.len()); topology.atoms.len()];
        let mut round_atom_neighborhoods = atom_neighborhoods.clone();
        let mut neighborhood_invariants = Vec::with_capacity(8);
        update_neighbor_layer(
            topology,
            0,
            layer,
            include_chirality,
            chiral_atoms,
            current_invariants,
            bond_invariants,
            &atom_neighborhoods,
            &mut round_atom_neighborhoods,
            &mut neighborhood_invariants,
        )
        .expect("the fixed M07 center has four source neighbors")
    }

    #[test]
    fn fingerprint_morgan_m07_cip_tag_and_radius_product_uses_fixed_hashes() {
        const LAYERS: [u32; 3] = [1, 2, 3];
        const BASE_HASHES: [u32; 3] = [1_169_060_574, 3_197_373_020, 3_364_079_190];
        const R_HASHES: [u32; 3] = [1_586_609_965, 3_427_473_679, 969_524_871];
        const S_HASHES: [u32; 3] = [1_586_609_964, 3_427_473_678, 969_524_870];
        const UNLABELED_HASHES: [u32; 3] = [1_586_609_967, 3_427_473_677, 969_524_889];
        const CURRENT_INVARIANTS: [u32; 5] = [0x1234_5678, 11, 22, 33, 44];
        const BOND_INVARIANTS: [u32; 4] = [1; 4];
        const CASES: [(ChiralTag, Option<&str>, [u32; 3], bool); 4] = [
            (ChiralTag::TetrahedralCw, Some("R"), R_HASHES, true),
            (ChiralTag::TetrahedralCw, Some("S"), S_HASHES, true),
            (ChiralTag::TetrahedralCw, None, UNLABELED_HASHES, true),
            (ChiralTag::Unspecified, Some("R"), BASE_HASHES, false),
        ];

        let mut calls = 0usize;
        for include_chirality in [false, true] {
            for (chiral_tag, cip_code, salted_hashes, should_be_chiral) in CASES {
                let cip_property = cip_code.map(|code| PropertyValue::String(code.to_owned()));
                let topology = m07_topology(chiral_tag, cip_property);
                for (layer_index, layer) in LAYERS.into_iter().enumerate() {
                    let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());
                    let actual = m07_layer_code(
                        &topology,
                        layer,
                        include_chirality,
                        &CURRENT_INVARIANTS,
                        &BOND_INVARIANTS,
                        &mut chiral_atoms,
                    );

                    let expected = if include_chirality {
                        salted_hashes[layer_index]
                    } else {
                        BASE_HASHES[layer_index]
                    };
                    assert_eq!(
                        actual, expected,
                        "tag={chiral_tag:?}, CIP={cip_code:?}, layer={layer}, include={include_chirality}"
                    );
                    assert_eq!(
                        chiral_atoms.contains(0),
                        include_chirality && should_be_chiral,
                        "chiralAtoms tag={chiral_tag:?}, CIP={cip_code:?}, layer={layer}, include={include_chirality}"
                    );
                    calls += 1;
                }
            }
        }

        assert_eq!(calls, 24);
    }

    #[test]
    fn fingerprint_morgan_m07_modeled_non_string_cip_properties_take_default_salt() {
        const EXPECTED_HASH: u32 = 3_427_473_677;
        const CURRENT_INVARIANTS: [u32; 5] = [0x1234_5678, 11, 22, 33, 44];
        const BOND_INVARIANTS: [u32; 4] = [1; 4];
        let cip_values = [
            PropertyValue::Int(17),
            PropertyValue::Double(1.25),
            PropertyValue::Bool(true),
            PropertyValue::String("other".to_owned()),
        ];
        let mut calls = 0usize;

        for cip_value in cip_values {
            let topology = m07_topology(ChiralTag::TetrahedralCw, Some(cip_value));
            let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());
            let actual = m07_layer_code(
                &topology,
                2,
                true,
                &CURRENT_INVARIANTS,
                &BOND_INVARIANTS,
                &mut chiral_atoms,
            );

            assert_eq!(actual, EXPECTED_HASH);
            assert!(chiral_atoms.contains(0));
            calls += 1;
        }

        assert_eq!(calls, 4);
    }

    #[test]
    fn fingerprint_morgan_m07_chiral_mask_rechecks_then_stays_set_across_layers() {
        const LAYER_ONE_DUPLICATE_BASE: u32 = 1_174_354_585;
        const LAYER_TWO_UNIQUE_R: u32 = 70_505_560;
        const LAYER_THREE_DUPLICATE_R: u32 = 3_736_007_148;
        const BOND_INVARIANTS: [u32; 4] = [1; 4];
        let topology = m07_topology(
            ChiralTag::TetrahedralCw,
            Some(PropertyValue::String("R".to_owned())),
        );
        let duplicate_neighbors = [0x1234_5678, 5, 5, 8, 11];
        let unique_neighbors = [0x1234_5678, 5, 6, 8, 11];
        let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());
        let mut calls = 0usize;

        let first = m07_layer_code(
            &topology,
            1,
            true,
            &duplicate_neighbors,
            &BOND_INVARIANTS,
            &mut chiral_atoms,
        );
        assert_eq!(first, LAYER_ONE_DUPLICATE_BASE);
        assert!(!chiral_atoms.contains(0));
        calls += 1;

        let second = m07_layer_code(
            &topology,
            2,
            true,
            &unique_neighbors,
            &BOND_INVARIANTS,
            &mut chiral_atoms,
        );
        assert_eq!(second, LAYER_TWO_UNIQUE_R);
        assert!(chiral_atoms.contains(0));
        calls += 1;

        let third = m07_layer_code(
            &topology,
            3,
            true,
            &duplicate_neighbors,
            &BOND_INVARIANTS,
            &mut chiral_atoms,
        );
        assert_eq!(third, LAYER_THREE_DUPLICATE_R);
        assert!(chiral_atoms.contains(0));
        calls += 1;
        assert_eq!(calls, 3);
    }

    #[test]
    fn fingerprint_morgan_m07_custom_non_single_bond_invariant_blocks_chirality() {
        const EXPECTED_BASE_HASH: u32 = 2_244_279_023;
        const CURRENT_INVARIANTS: [u32; 5] = [0x1234_5678, 5, 6, 8, 11];
        const CUSTOM_BOND_INVARIANTS: [u32; 4] = [2, 1, 1, 1];
        let topology = m07_topology(
            ChiralTag::TetrahedralCw,
            Some(PropertyValue::String("R".to_owned())),
        );
        let mut chiral_atoms = MorganChiralAtoms::new(topology.atoms.len());

        let actual = m07_layer_code(
            &topology,
            2,
            true,
            &CURRENT_INVARIANTS,
            &CUSTOM_BOND_INVARIANTS,
            &mut chiral_atoms,
        );

        assert_eq!(actual, EXPECTED_BASE_HASH);
        assert!(!chiral_atoms.contains(0));
    }

    fn m08_mask(bit_count: usize, bits: &[usize]) -> MorganBondEnvironment {
        let mut mask = MorganBondEnvironment::new(bit_count);
        for &bit in bits {
            mask.set(bit);
        }
        mask
    }

    fn m08_word_mask(bit_count: usize, word: u64) -> MorganBondEnvironment {
        let mut mask = MorganBondEnvironment::new(bit_count);
        assert_eq!(mask.blocks.len(), 1);
        mask.blocks[0] = word;
        mask
    }

    #[test]
    fn fingerprint_morgan_m08_dynamic_bitset_order_matches_boost_logical_bits() {
        // In set-bit-index-vector order [0, 5] would precede [1]; the pinned
        // equal-width Boost comparison examines the highest differing bit.
        let low_index = m08_mask(8, &[1]);
        let high_index = m08_mask(8, &[0, 5]);
        assert_eq!(low_index.cmp(&high_index), Ordering::Less);

        // The equal-width path compares packed words from the highest word.
        let low_word = m08_mask(130, &[64]);
        let high_word = m08_mask(130, &[127]);
        assert_eq!(low_word.cmp(&high_word), Ordering::Less);

        // Boost's different-width path compares top logical bits before size.
        let longer_low = m08_mask(4, &[0]);
        let shorter_high = m08_mask(2, &[1]);
        assert_eq!(longer_low.cmp(&shorter_high), Ordering::Less);
        let shorter_prefix = m08_mask(2, &[1]);
        let longer_prefix = m08_mask(3, &[2]);
        assert_eq!(shorter_prefix.cmp(&longer_prefix), Ordering::Less);
        assert_eq!(
            MorganBondEnvironment::new(0).cmp(&m08_mask(1, &[0])),
            Ordering::Less
        );

        let mut seen = HashSet::new();
        seen.insert(low_index.clone());
        assert!(seen.contains(&low_index));
        assert!(!seen.contains(&m08_mask(9, &[1])));
        assert_eq!(seen.len(), 1);
    }

    #[test]
    fn fingerprint_morgan_m08_collector_orders_all_ring_and_line_modes() {
        const SYMMETRIC_LINE_INPUT: [(u64, u32, u32); 2] = [(1, 11, 1), (1, 11, 0)];
        const SYMMETRIC_LINE_SORTED: [(u64, u32, u32); 2] = [(1, 11, 0), (1, 11, 1)];
        const SYMMETRIC_LINE_FALSE: [(u32, u32, u32); 1] = [(11, 0, 2)];
        const SYMMETRIC_LINE_TRUE: [(u32, u32, u32); 2] = [(11, 0, 2), (11, 1, 2)];
        const SYMMETRIC_LINE_UNIQUE: [u64; 1] = [1];
        const ASYMMETRIC_LINE_INPUT: [(u64, u32, u32); 4] =
            [(6, 50, 2), (1, 70, 0), (4, 40, 3), (3, 10, 1)];
        const ASYMMETRIC_LINE_SORTED: [(u64, u32, u32); 4] =
            [(1, 70, 0), (3, 10, 1), (4, 40, 3), (6, 50, 2)];
        const ASYMMETRIC_LINE_EMITTED: [(u32, u32, u32); 4] =
            [(70, 0, 2), (10, 1, 2), (40, 3, 2), (50, 2, 2)];
        const ASYMMETRIC_LINE_UNIQUE: [u64; 4] = [1, 3, 4, 6];
        const SYMMETRIC_RING_INPUT: [(u64, u32, u32); 3] = [(7, 44, 2), (7, 44, 0), (7, 44, 1)];
        const SYMMETRIC_RING_SORTED: [(u64, u32, u32); 3] = [(7, 44, 0), (7, 44, 1), (7, 44, 2)];
        const SYMMETRIC_RING_FALSE: [(u32, u32, u32); 1] = [(44, 0, 2)];
        const SYMMETRIC_RING_TRUE: [(u32, u32, u32); 3] = [(44, 0, 2), (44, 1, 2), (44, 2, 2)];
        const SYMMETRIC_RING_UNIQUE: [u64; 1] = [7];
        const ASYMMETRIC_RING_INPUT: [(u64, u32, u32); 4] =
            [(15, 30, 1), (14, 40, 3), (15, 10, 2), (15, 20, 0)];
        const ASYMMETRIC_RING_SORTED: [(u64, u32, u32); 4] =
            [(14, 40, 3), (15, 10, 2), (15, 20, 0), (15, 30, 1)];
        const ASYMMETRIC_RING_FALSE: [(u32, u32, u32); 2] = [(40, 3, 2), (10, 2, 2)];
        const ASYMMETRIC_RING_TRUE: [(u32, u32, u32); 4] =
            [(40, 3, 2), (10, 2, 2), (20, 0, 2), (30, 1, 2)];
        const ASYMMETRIC_RING_UNIQUE: [u64; 2] = [14, 15];

        let cases: [(
            &str,
            TopologyBlock,
            &[(u64, u32, u32)],
            &[(u64, u32, u32)],
            &[(u32, u32, u32)],
            &[(u32, u32, u32)],
            &[u64],
            &[usize],
        ); 4] = [
            (
                "symmetric line",
                m06_topology(2, &[(0, 1, BondOrder::Single)]),
                &SYMMETRIC_LINE_INPUT,
                &SYMMETRIC_LINE_SORTED,
                &SYMMETRIC_LINE_FALSE,
                &SYMMETRIC_LINE_TRUE,
                &SYMMETRIC_LINE_UNIQUE,
                &[1],
            ),
            (
                "asymmetric line",
                m06_topology(
                    4,
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Double),
                        (2, 3, BondOrder::Single),
                    ],
                ),
                &ASYMMETRIC_LINE_INPUT,
                &ASYMMETRIC_LINE_SORTED,
                &ASYMMETRIC_LINE_EMITTED,
                &ASYMMETRIC_LINE_EMITTED,
                &ASYMMETRIC_LINE_UNIQUE,
                &[],
            ),
            (
                "symmetric ring",
                m06_topology(
                    3,
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Single),
                        (2, 0, BondOrder::Single),
                    ],
                ),
                &SYMMETRIC_RING_INPUT,
                &SYMMETRIC_RING_SORTED,
                &SYMMETRIC_RING_FALSE,
                &SYMMETRIC_RING_TRUE,
                &SYMMETRIC_RING_UNIQUE,
                &[1, 2],
            ),
            (
                "asymmetric ring",
                m06_topology(
                    4,
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Single),
                        (2, 0, BondOrder::Double),
                        (2, 3, BondOrder::Single),
                    ],
                ),
                &ASYMMETRIC_RING_INPUT,
                &ASYMMETRIC_RING_SORTED,
                &ASYMMETRIC_RING_FALSE,
                &ASYMMETRIC_RING_TRUE,
                &ASYMMETRIC_RING_UNIQUE,
                &[0, 1],
            ),
        ];

        let mut calls = 0usize;
        for (
            name,
            topology,
            input_rows,
            expected_sorted_rows,
            expected_false,
            expected_true,
            expected_unique_masks,
            expected_dead_atoms,
        ) in cases
        {
            let atom_invariants = vec![1; topology.atoms.len()];
            let mut include_atoms = MorganBondEnvironment::new(topology.atoms.len());
            for atom_id in 0..topology.atoms.len() {
                include_atoms.set(atom_id);
            }

            for include_redundant in [false, true] {
                let mut candidates = input_rows
                    .iter()
                    .map(|&(mask, code, atom_id)| {
                        (m08_word_mask(topology.bonds.len(), mask), code, atom_id)
                    })
                    .collect::<Vec<_>>();
                let mut neighborhoods = HashSet::new();
                let mut dead_atoms = MorganBondEnvironment::new(topology.atoms.len());
                let mut result = Vec::new();

                collect_morgan_layer::<u32>(
                    &mut candidates,
                    include_redundant,
                    false,
                    &atom_invariants,
                    &include_atoms,
                    &mut neighborhoods,
                    &mut dead_atoms,
                    1,
                    &topology,
                    &mut result,
                );

                let actual_sorted = candidates
                    .iter()
                    .map(|(mask, code, atom_id)| (mask.clone(), *code, *atom_id))
                    .collect::<Vec<_>>();
                let expected_sorted = expected_sorted_rows
                    .iter()
                    .map(|&(mask, code, atom_id)| {
                        (m08_word_mask(topology.bonds.len(), mask), code, atom_id)
                    })
                    .collect::<Vec<_>>();
                assert_eq!(actual_sorted, expected_sorted, "{name}");

                let expected_emitted = if include_redundant {
                    expected_true
                } else {
                    expected_false
                };
                let actual_emitted = result
                    .iter()
                    .map(|environment| {
                        (
                            environment.d_code,
                            environment.d_atom_id,
                            environment.d_layer,
                        )
                    })
                    .collect::<Vec<_>>();
                assert_eq!(
                    actual_emitted, expected_emitted,
                    "{name}, redundant={include_redundant}"
                );
                assert!(
                    result
                        .iter()
                        .all(|environment| { std::ptr::eq(environment.topology, &topology) })
                );

                let actual_dead = (0..topology.atoms.len())
                    .filter(|&atom_id| dead_atoms.contains(atom_id))
                    .collect::<Vec<_>>();
                let expected_dead = if include_redundant {
                    &[][..]
                } else {
                    expected_dead_atoms
                };
                assert_eq!(
                    actual_dead, expected_dead,
                    "{name}, redundant={include_redundant}"
                );

                assert_eq!(neighborhoods.len(), expected_unique_masks.len(), "{name}");
                for &mask in expected_unique_masks {
                    assert!(neighborhoods.contains(&m08_word_mask(topology.bonds.len(), mask)));
                }
                calls += 1;
            }
        }
        assert_eq!(calls, 8);
    }

    #[test]
    fn fingerprint_morgan_m08_source_filter_and_neighborhood_state_order() {
        let topology = m06_topology(
            4,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let atom_invariants = [0, 7, 8, 9];
        let include_atoms = m08_mask(4, &[0, 1]);
        let mut neighborhoods = HashSet::new();
        let mut dead_atoms = MorganBondEnvironment::new(4);
        let mut result = Vec::new();
        let mut calls = 0usize;

        let mut first_layer: Vec<MorganLayerCandidate> = vec![
            (m08_word_mask(3, 2), 40, 3),
            (m08_word_mask(3, 1), 30, 2),
            (m08_word_mask(3, 1), 10, 0),
            (m08_word_mask(3, 1), 20, 1),
        ];
        collect_morgan_layer::<u32>(
            &mut first_layer,
            false,
            true,
            &atom_invariants,
            &include_atoms,
            &mut neighborhoods,
            &mut dead_atoms,
            0,
            &topology,
            &mut result,
        );
        calls += 1;

        assert_eq!(
            result
                .iter()
                .map(|environment| (
                    environment.d_code,
                    environment.d_atom_id,
                    environment.d_layer
                ))
                .collect::<Vec<_>>(),
            [(20, 1, 1)]
        );
        assert!(neighborhoods.contains(&m08_word_mask(3, 1)));
        assert!(!neighborhoods.contains(&m08_word_mask(3, 2)));
        assert!(dead_atoms.contains(2));
        assert!(!dead_atoms.contains(0));
        assert!(!dead_atoms.contains(3));

        // A later layer sees the same mask first, so source marks this center
        // dead before reaching its false includeAtoms test.
        let mut second_layer = vec![(m08_word_mask(3, 1), 50, 3)];
        collect_morgan_layer::<u32>(
            &mut second_layer,
            false,
            true,
            &atom_invariants,
            &include_atoms,
            &mut neighborhoods,
            &mut dead_atoms,
            1,
            &topology,
            &mut result,
        );
        calls += 1;

        assert!(dead_atoms.contains(2));
        assert!(dead_atoms.contains(3));
        assert_eq!(result.len(), 1);
        assert_eq!(neighborhoods.len(), 1);
        assert_eq!(calls, 2);
    }

    #[test]
    fn fingerprint_morgan_m09_radius_loop_matches_all_frozen_source_rows() {
        const MOLECULES: [&str; 5] = ["", "CCO", "c1cc[nH]c1", "c1ccccc1", "[Na+].[Cl-]"];
        const RADII: [u32; 4] = [0, 1, 2, 3];
        const BOOLEAN_VALUES: [bool; 2] = [false, true];
        const EXPECTED: [(&str, u32, bool, bool, &[(u64, i32)]); 80] = [
            ("", 0, false, false, &[]),
            ("", 0, false, true, &[]),
            ("", 0, true, false, &[]),
            ("", 0, true, true, &[]),
            ("", 1, false, false, &[]),
            ("", 1, false, true, &[]),
            ("", 1, true, false, &[]),
            ("", 1, true, true, &[]),
            ("", 2, false, false, &[]),
            ("", 2, false, true, &[]),
            ("", 2, true, false, &[]),
            ("", 2, true, true, &[]),
            ("", 3, false, false, &[]),
            ("", 3, false, true, &[]),
            ("", 3, true, false, &[]),
            ("", 3, true, true, &[]),
            (
                "CCO",
                0,
                false,
                false,
                &[(864662311, 1), (2245384272, 1), (2246728737, 1)],
            ),
            (
                "CCO",
                0,
                false,
                true,
                &[(864662311, 1), (2245384272, 1), (2246728737, 1)],
            ),
            (
                "CCO",
                0,
                true,
                false,
                &[(864662311, 1), (2245384272, 1), (2246728737, 1)],
            ),
            (
                "CCO",
                0,
                true,
                true,
                &[(864662311, 1), (2245384272, 1), (2246728737, 1)],
            ),
            (
                "CCO",
                1,
                false,
                false,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                1,
                false,
                true,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                1,
                true,
                false,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                1,
                true,
                true,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                2,
                false,
                false,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                2,
                false,
                true,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                2,
                true,
                false,
                &[
                    (407009239, 1),
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (2679480461, 1),
                    (3542456614, 1),
                    (3732090711, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                2,
                true,
                true,
                &[
                    (407009239, 1),
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (2679480461, 1),
                    (3542456614, 1),
                    (3732090711, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                3,
                false,
                false,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                3,
                false,
                true,
                &[
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (3542456614, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                3,
                true,
                false,
                &[
                    (407009239, 1),
                    (589800426, 1),
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (2679480461, 1),
                    (2708961779, 1),
                    (2962971722, 1),
                    (3542456614, 1),
                    (3732090711, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "CCO",
                3,
                true,
                true,
                &[
                    (407009239, 1),
                    (589800426, 1),
                    (864662311, 1),
                    (1535166686, 1),
                    (2245384272, 1),
                    (2246728737, 1),
                    (2679480461, 1),
                    (2708961779, 1),
                    (2962971722, 1),
                    (3542456614, 1),
                    (3732090711, 1),
                    (4018048386, 1),
                ],
            ),
            (
                "c1cc[nH]c1",
                0,
                false,
                false,
                &[(2132511834, 1), (3218693969, 4)],
            ),
            (
                "c1cc[nH]c1",
                0,
                false,
                true,
                &[(2132511834, 1), (3218693969, 4)],
            ),
            (
                "c1cc[nH]c1",
                0,
                true,
                false,
                &[(2132511834, 1), (3218693969, 4)],
            ),
            (
                "c1cc[nH]c1",
                0,
                true,
                true,
                &[(2132511834, 1), (3218693969, 4)],
            ),
            (
                "c1cc[nH]c1",
                1,
                false,
                false,
                &[
                    (98513984, 2),
                    (2132511834, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                1,
                false,
                true,
                &[
                    (98513984, 2),
                    (2132511834, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                1,
                true,
                false,
                &[
                    (98513984, 2),
                    (2132511834, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                1,
                true,
                true,
                &[
                    (98513984, 2),
                    (2132511834, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                2,
                false,
                false,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                2,
                false,
                true,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                2,
                true,
                false,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                2,
                true,
                true,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                3,
                false,
                false,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2266426494, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                3,
                false,
                true,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2266426494, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                ],
            ),
            (
                "c1cc[nH]c1",
                3,
                true,
                false,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2266426494, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                    (3541106057, 2),
                    (4054084389, 2),
                ],
            ),
            (
                "c1cc[nH]c1",
                3,
                true,
                true,
                &[
                    (98513984, 2),
                    (116898731, 2),
                    (1482649460, 2),
                    (2132511834, 1),
                    (2266426494, 1),
                    (2293755984, 1),
                    (2654043257, 1),
                    (2753863138, 2),
                    (3218693969, 4),
                    (3541106057, 2),
                    (4054084389, 2),
                ],
            ),
            ("c1ccccc1", 0, false, false, &[(3218693969, 6)]),
            ("c1ccccc1", 0, false, true, &[(3218693969, 6)]),
            ("c1ccccc1", 0, true, false, &[(3218693969, 6)]),
            ("c1ccccc1", 0, true, true, &[(3218693969, 6)]),
            (
                "c1ccccc1",
                1,
                false,
                false,
                &[(98513984, 6), (3218693969, 6)],
            ),
            (
                "c1ccccc1",
                1,
                false,
                true,
                &[(98513984, 6), (3218693969, 6)],
            ),
            (
                "c1ccccc1",
                1,
                true,
                false,
                &[(98513984, 6), (3218693969, 6)],
            ),
            ("c1ccccc1", 1, true, true, &[(98513984, 6), (3218693969, 6)]),
            (
                "c1ccccc1",
                2,
                false,
                false,
                &[(98513984, 6), (2763854213, 6), (3218693969, 6)],
            ),
            (
                "c1ccccc1",
                2,
                false,
                true,
                &[(98513984, 6), (2763854213, 6), (3218693969, 6)],
            ),
            (
                "c1ccccc1",
                2,
                true,
                false,
                &[(98513984, 6), (2763854213, 6), (3218693969, 6)],
            ),
            (
                "c1ccccc1",
                2,
                true,
                true,
                &[(98513984, 6), (2763854213, 6), (3218693969, 6)],
            ),
            (
                "c1ccccc1",
                3,
                false,
                false,
                &[
                    (98513984, 6),
                    (2763854213, 6),
                    (3218693969, 6),
                    (3741631696, 1),
                ],
            ),
            (
                "c1ccccc1",
                3,
                false,
                true,
                &[
                    (98513984, 6),
                    (2763854213, 6),
                    (3218693969, 6),
                    (3741631696, 1),
                ],
            ),
            (
                "c1ccccc1",
                3,
                true,
                false,
                &[
                    (98513984, 6),
                    (2763854213, 6),
                    (3218693969, 6),
                    (3741631696, 6),
                ],
            ),
            (
                "c1ccccc1",
                3,
                true,
                true,
                &[
                    (98513984, 6),
                    (2763854213, 6),
                    (3218693969, 6),
                    (3741631696, 6),
                ],
            ),
            (
                "[Na+].[Cl-]",
                0,
                false,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                0,
                false,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                0,
                true,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                0,
                true,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                1,
                false,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                1,
                false,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                1,
                true,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                1,
                true,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                2,
                false,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                2,
                false,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                2,
                true,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                2,
                true,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                3,
                false,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                3,
                false,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                3,
                true,
                false,
                &[(3737048253, 1), (3855292234, 1)],
            ),
            (
                "[Na+].[Cl-]",
                3,
                true,
                true,
                &[(3737048253, 1), (3855292234, 1)],
            ),
        ];

        let arguments = FingerprintFuncArguments::default();
        let mut calls = 0usize;
        for smiles in MOLECULES {
            let (record, prepared_valence) = fully_prepared_record_with_valence(smiles);
            let rings = fast_find_rings(&record.topology)
                .expect("fixed Morgan topology has source ring state");
            let mut generator = get_morgan_generator(&MorganParams::default())
                .expect("source-default Morgan arguments are valid");
            let atom_invariants = generator
                .atom_invariants
                .get_atom_invariants(&record.topology, &prepared_valence, &rings)
                .expect("fixed prepared topology yields default Morgan atom invariants");
            let bond_invariants = generator
                .bond_invariants
                .get_bond_invariants(&record.topology);

            for radius in RADII {
                for include_redundant in BOOLEAN_VALUES {
                    for only_nonzero in BOOLEAN_VALUES {
                        let (
                            expected_smiles,
                            expected_radius,
                            expected_redundant,
                            expected_nonzero,
                            expected_rows,
                        ) = EXPECTED[calls];
                        assert_eq!(
                            (
                                expected_smiles,
                                expected_radius,
                                expected_redundant,
                                expected_nonzero
                            ),
                            (smiles, radius, include_redundant, only_nonzero),
                            "the fixed table follows the complete declared product"
                        );

                        generator.radius = radius;
                        generator.include_redundant_environments = include_redundant;
                        generator.only_nonzero_invariants = only_nonzero;
                        let environments = generate_morgan_environments::<u64>(
                            &record.topology,
                            &generator,
                            &arguments,
                            &atom_invariants,
                            &bond_invariants,
                        )
                        .expect("fixed default invariant rows satisfy source preconditions");
                        let mut actual_counts = BTreeMap::<u64, i32>::new();
                        for environment in environments {
                            *actual_counts.entry(environment.d_code).or_default() += 1;
                        }
                        let actual_counts = actual_counts.into_iter().collect::<Vec<_>>();
                        assert_eq!(
                            actual_counts, expected_rows,
                            "{smiles}, radius={radius}, redundant={include_redundant}, only_nonzero={only_nonzero}"
                        );
                        calls += 1;
                    }
                }
            }
        }
        assert_eq!(calls, 80);
        assert_eq!(calls, EXPECTED.len());
    }

    #[test]
    fn fingerprint_morgan_g07_original_invariants_and_prepared_environment_product() {
        const INPUT: &str = "F[C@H](Cl)Br";
        const DONE_VALUES: [Option<&str>; 3] = [None, Some("0"), Some("1")];
        const DEFAULT_ATOM_INVARIANTS: [u32; 4] =
            [882_399_112, 2_245_273_601, 1_016_841_875, 3_612_926_680];
        const DEFAULT_BOND_INVARIANTS: [u32; 3] = [1, 1, 1];
        const CUSTOM_ATOM_INVARIANTS: [u32; 4] = [101, 202, 303, 404];
        const CUSTOM_BOND_INVARIANTS: [u32; 3] = [7, 11, 13];

        let original = parse_smiles(INPUT, &SmilesParseParams::default())
            .expect("the fixed raw chiral input parses");
        assert_eq!(
            original
                .topology
                .atoms
                .iter()
                .map(Atom::atomic_number)
                .collect::<Vec<_>>(),
            [9, 6, 17, 35]
        );
        assert_eq!(original.topology.bonds.len(), 3);
        let valence = cosmolkit_core::assign_valence(
            &original.topology,
            &cosmolkit_core::ValenceParams::default(),
        )
        .expect("the fixed source-valid raw topology has prepared valence");
        let rings = fast_find_rings(&original.topology)
            .expect("the fixed raw topology has initialized ring state");
        assert!(rings.is_initialized());
        let valence_before = valence.clone();
        let rings_before = rings.clone();
        let original_topology_before = original.topology.clone();
        let original_properties_before = original.properties.clone();
        let mut calls = 0usize;

        for include_chirality in [false, true] {
            for done_value in DONE_VALUES {
                for use_custom_atom_invariants in [false, true] {
                    for use_custom_bond_invariants in [false, true] {
                        let mut record = original.clone();
                        if let Some(done_value) = done_value {
                            record
                                .properties
                                .set_prop("_StereochemDone", done_value)
                                .expect("the fixed source-state marker is valid");
                        }
                        let record_topology_before = record.topology.clone();
                        let record_properties_before = record.properties.clone();
                        let should_prepare = include_chirality && done_value.is_none();

                        let params = MorganParams {
                            radius: 0,
                            include_chirality,
                            ..MorganParams::default()
                        };
                        let generator =
                            get_morgan_generator(&params).expect("fixed Morgan params are valid");
                        let custom_atom =
                            use_custom_atom_invariants.then_some(&CUSTOM_ATOM_INVARIANTS[..]);
                        let custom_bond =
                            use_custom_bond_invariants.then_some(&CUSTOM_BOND_INVARIANTS[..]);
                        let arguments =
                            FingerprintFuncArguments::new(None, None, custom_atom, custom_bond, -1);
                        let expected_atoms = if use_custom_atom_invariants {
                            CUSTOM_ATOM_INVARIANTS
                        } else {
                            DEFAULT_ATOM_INVARIANTS
                        };
                        let expected_bonds = if use_custom_bond_invariants {
                            CUSTOM_BOND_INVARIANTS
                        } else {
                            DEFAULT_BOND_INVARIANTS
                        };
                        let original_topology_address = &record.topology as *const TopologyBlock;
                        let original_properties_address =
                            &record.properties as *const cosmolkit_model::MoleculeProperties;

                        let (actual_atoms, actual_bonds, actual_rows, selected_is_original) =
                            with_morgan_environment_inputs(
                                &record.topology,
                                &record.properties,
                                &valence,
                                &rings,
                                &generator,
                                &arguments,
                                |selected_topology,
                                 selected_properties,
                                 atom_invariants,
                                 bond_invariants| {
                                    assert_eq!(
                                        atom_invariants, expected_atoms,
                                        "atom invariants come from the original molecule"
                                    );
                                    assert_eq!(
                                        bond_invariants, expected_bonds,
                                        "bond invariants come from the original molecule"
                                    );
                                    assert_eq!(
                                        std::ptr::eq(selected_topology, original_topology_address),
                                        !should_prepare,
                                        "source-selected environment topology ownership"
                                    );
                                    assert_eq!(
                                        std::ptr::eq(
                                            selected_properties,
                                            original_properties_address
                                        ),
                                        !should_prepare,
                                        "source-selected environment property ownership"
                                    );
                                    assert_eq!(
                                        selected_properties.prop("_StereochemDone"),
                                        if should_prepare {
                                            Some("1")
                                        } else {
                                            done_value
                                        }
                                    );
                                    assert_eq!(
                                        selected_properties.is_prop_computed("_StereochemDone"),
                                        should_prepare
                                    );
                                    assert_eq!(
                                        selected_topology.atoms[1]
                                            .prop("_CIPCode")
                                            .and_then(|value| value.as_string().ok()),
                                        should_prepare.then_some("R"),
                                        "only the source-selected prepared copy has assigned CIP"
                                    );

                                    let environments = generate_morgan_environments::<u32>(
                                        selected_topology,
                                        &generator,
                                        &arguments,
                                        atom_invariants,
                                        bond_invariants,
                                    )?;
                                    assert_eq!(environments.len(), 4);
                                    assert!(
                                        environments.iter().all(|environment| std::ptr::eq(
                                            environment.topology,
                                            selected_topology
                                        )),
                                        "each environment retains the source-selected molecule"
                                    );
                                    let rows = environments
                                        .iter()
                                        .map(|environment| {
                                            (
                                                environment.d_code,
                                                environment.d_atom_id,
                                                environment.d_layer,
                                            )
                                        })
                                        .collect::<Vec<_>>();
                                    Ok::<_, crate::MorganError>((
                                        atom_invariants.to_vec(),
                                        bond_invariants.to_vec(),
                                        rows,
                                        std::ptr::eq(selected_topology, original_topology_address),
                                    ))
                                },
                            )
                            .expect("the source-valid G07 composition succeeds");

                        assert_eq!(actual_atoms, expected_atoms);
                        assert_eq!(actual_bonds, expected_bonds);
                        assert_eq!(
                            actual_rows,
                            expected_atoms
                                .iter()
                                .enumerate()
                                .map(|(atom_index, &invariant)| {
                                    (invariant, atom_index as u32, 0)
                                })
                                .collect::<Vec<_>>(),
                            "radius-zero environments preserve the complete literal atom rows"
                        );
                        assert_eq!(selected_is_original, !should_prepare);
                        assert_eq!(record.topology, record_topology_before);
                        assert_eq!(record.properties, record_properties_before);
                        calls += 1;
                    }
                }
            }
        }

        assert_eq!(calls, 24, "2 × 3 × 2 × 2 complete G07 source product");
        assert_eq!(
            valence, valence_before,
            "the supplied valence is reused unchanged"
        );
        assert_eq!(
            rings, rings_before,
            "the supplied ring state is reused unchanged"
        );
        assert_eq!(original.topology, original_topology_before);
        assert_eq!(original.properties, original_properties_before);
    }
}
