//! Packed-code helpers for supplied integer atom codes.
//! Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8.
//! Atom-code generation, CIP and molecular preparation have separate owners.
//! Empty paths and out-of-width shifts are undefined in C++; InvalidArguments
//! is the explicit ROOT-approved COSMolKit safety policy. Reversed hash paths
//! longer than UINTMAX would not terminate upstream and are rejected here.
//! Huge forward paths retain upstream u32 comparison indexes and termination.
use crate::{FingerprintError, hash::hash_combine};

const CODE_BITS: u32 = crate::AtomPairsParameters::code_size();
const CHIRAL_BITS: u32 = crate::AtomPairsParameters::num_chiral_bits();
const PATH_BITS: u32 = crate::AtomPairsParameters::num_path_bits();
const MAX_PATH_LENGTH: u32 = (1 << PATH_BITS) - 1;

// Callers establish nonempty paths. This helper is the only paired-end
// comparison; source indexes are unsigned32 even with a64bit vector size.
// No sort/clone/alloc; <=floor(UINTMAX/2) paired comparisons for huge paths.
fn reverse_packed_path(path_codes: &[u32]) -> bool {
    // RDKit✔️✔️:   bool reverseIt = false;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   unsigned int j = pathCodes.size() - 1;
    // RDKit✔️✔️:   while (i < j) {
    // RDKit✔️✔️:     if (pathCodes[i] > pathCodes[j]) {
    // RDKit✔️✔️:       reverseIt = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else if (pathCodes[i] < pathCodes[j]) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++i;
    // RDKit✔️✔️:     --j;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool reverseIt = false;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   unsigned int j = pathCodes.size() - 1;
    // RDKit✔️✔️:   while (i < j) {
    // RDKit✔️✔️:     if (pathCodes[i] > pathCodes[j]) {
    // RDKit✔️✔️:       reverseIt = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else if (pathCodes[i] < pathCodes[j]) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++i;
    // RDKit✔️✔️:     --j;
    // RDKit✔️✔️:   }
    let mut i = 0u32;
    let mut j = (path_codes.len() - 1) as u32;
    while i < j {
        if path_codes[i as usize] > path_codes[j as usize] {
            return true;
        } else if path_codes[i as usize] < path_codes[j as usize] {
            break;
        }
        i += 1;
        j -= 1;
    }
    false
}

/// Source u32 min/max pair packing; dist must be strictly<31.
/// Explicit boolean, no atom generation/preparation/cache or Molecule.
/// Full u32 codes including highbits are accepted, naturally truncated.
pub fn atom_pair_code(
    code_i: u32,
    code_j: u32,
    distance: u32,
    include_chirality: bool,
) -> Result<u32, FingerprintError> {
    // RDKit✔️✔️: const unsigned int numTypeBits = 4;
    // RDKit✔️✔️: const unsigned int numPiBits = 2;
    // RDKit❌❌: const unsigned int maxNumPi = (1 << numPiBits) - 1;
    // RDKit✔️✔️: const unsigned int numBranchBits = 3;
    // RDKit❌❌: const unsigned int maxNumBranches = (1 << numBranchBits) - 1;
    // RDKit✔️✔️: const unsigned int numChiralBits = 2;
    // RDKit✔️✔️: const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
    // RDKit✔️✔️: const unsigned int numPathBits = 5;
    // RDKit✔️✔️: const unsigned int maxPathLen = (1 << numPathBits) - 1;
    // RDKit✔️✔️: std::uint32_t getAtomPairCode(std::uint32_t codeI, std::uint32_t codeJ,
    // RDKit✔️✔️:                               unsigned int dist, bool includeChirality) {
    // RDKit✔️✔️:   PRECONDITION(dist < maxPathLen, "dist too long");
    // RDKit✔️✔️:   std::uint32_t res = dist;
    // RDKit✔️✔️:   res |= std::min(codeI, codeJ) << numPathBits;
    // RDKit✔️✔️:   res |= std::max(codeI, codeJ)
    // RDKit✔️✔️:          << (numPathBits + codeSize + (includeChirality ? numChiralBits : 0));
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Source Invar::Invariant PRECONDITION mapping preserves its text.
    if distance >= MAX_PATH_LENGTH {
        return Err(FingerprintError::PreconditionViolation {
            what: "dist too long",
        });
    }
    let mut result = distance;
    result |= code_i.min(code_j) << PATH_BITS;
    result |= code_i.max(code_j)
        << (PATH_BITS + CODE_BITS + if include_chirality { CHIRAL_BITS } else { 0 });
    Ok(result)
}

/// Canonical source OR packing in u64; full u32 elements, no masks.
/// Empty or any executed shift>=64 has no C++defined number.
/// InvalidArguments is the approved COSMolKit safety policy, not an RDKit error.
pub fn topological_torsion_code(
    path_codes: &[u32],
    include_chirality: bool,
) -> Result<u64, FingerprintError> {
    // RDKit✔️✔️: const unsigned int numTypeBits = 4;
    // RDKit✔️✔️: const unsigned int numPiBits = 2;
    // RDKit❌❌: const unsigned int maxNumPi = (1 << numPiBits) - 1;
    // RDKit✔️✔️: const unsigned int numBranchBits = 3;
    // RDKit❌❌: const unsigned int maxNumBranches = (1 << numBranchBits) - 1;
    // RDKit✔️✔️: const unsigned int numChiralBits = 2;
    // RDKit✔️✔️: const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
    // RDKit✔️✔️: const unsigned int numPathBits = 5;
    // RDKit✔️✔️: const unsigned int maxPathLen = (1 << numPathBits) - 1;
    // RDKit✔️✔️: std::uint64_t getTopologicalTorsionCode(
    // RDKit✔️✔️:     const std::vector<std::uint32_t> &pathCodes, bool includeChirality) {
    // RDKit✔️✔️:   bool reverseIt = false;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   unsigned int j = pathCodes.size() - 1;
    // RDKit✔️✔️:   while (i < j) {
    // RDKit✔️✔️:     if (pathCodes[i] > pathCodes[j]) {
    // RDKit✔️✔️:       reverseIt = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else if (pathCodes[i] < pathCodes[j]) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++i;
    // RDKit✔️✔️:     --j;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int shiftSize = codeSize + (includeChirality ? numChiralBits : 0);
    // RDKit✔️✔️:   std::uint64_t res = 0;
    // RDKit✔️✔️:   if (reverseIt) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pathCodes.size(); ++i) {
    // RDKit✔️✔️:       res |= static_cast<std::uint64_t>(pathCodes[pathCodes.size() - i - 1])
    // RDKit✔️✔️:              << (shiftSize * i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pathCodes.size(); ++i) {
    // RDKit✔️✔️:       res |= static_cast<std::uint64_t>(pathCodes[i]) << (shiftSize * i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if path_codes.is_empty() {
        return Err(FingerprintError::InvalidArguments {
            reason: "empty topological torsion code path is undefined in pinned source",
        });
    }
    let shift_size = CODE_BITS + if include_chirality { CHIRAL_BITS } else { 0 };
    // Last executed index is len-1; limit=8false/6true even if all codes0.
    if path_codes.len() > ((u64::BITS - 1) / shift_size + 1) as usize {
        return Err(FingerprintError::InvalidArguments {
            reason: "topological torsion code executes a shift outside source uint64 width",
        });
    }
    let reverse_it = reverse_packed_path(path_codes);
    let mut result = 0u64;
    for i in 0..path_codes.len() {
        let index = if reverse_it {
            path_codes.len() - i - 1
        } else {
            i
        };
        result |= u64::from(path_codes[index]) << (shift_size * i as u32);
    }
    Ok(result)
}

/// Canonical u32 source hash; all nonempty terminating lengths are supported.
/// Hash has no chirality parameter or packed-code maximum length.
/// Empty is source UB; enormous reversed paths source-nonterminate.
/// Reversed-path nontermination is rejected under the approved safety policy.
pub fn topological_torsion_hash(path_codes: &[u32]) -> Result<u32, FingerprintError> {
    // RDKit✔️✔️: std::uint32_t getTopologicalTorsionHash(
    // RDKit✔️✔️:     const std::vector<std::uint32_t> &pathCodes) {
    // RDKit✔️✔️:   bool reverseIt = false;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   unsigned int j = pathCodes.size() - 1;
    // RDKit✔️✔️:   while (i < j) {
    // RDKit✔️✔️:     if (pathCodes[i] > pathCodes[j]) {
    // RDKit✔️✔️:       reverseIt = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else if (pathCodes[i] < pathCodes[j]) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++i;
    // RDKit✔️✔️:     --j;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::uint32_t res = 0;
    // RDKit✔️✔️:   if (reverseIt) {
    // RDKit❗✔️:     for (unsigned int i = 0; i < pathCodes.size(); ++i) {
    // RDKit✔️✔️:       gboost::hash_combine(res, pathCodes[pathCodes.size() - i - 1]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned int pathCode : pathCodes) {
    // RDKit✔️✔️:       gboost::hash_combine(res, pathCode);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: namespace std {
    // RDKit✔️✔️: typedef std::uint32_t hash_result_t;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: inline std::hash_result_t hash_value(unsigned int v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: #if BOOST_WORKAROUND(BOOST_MSVC, < 1300)
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: inline void hash_combine(std::hash_result_t& seed, T& v)
    // RDKit✔️✔️: #else
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: inline void hash_combine(std::hash_result_t& seed, T const& v)
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️: {
    // RDKit✔️✔️:   gboost::hash<T> hasher;
    // RDKit✔️✔️:   seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️: }
    if path_codes.is_empty() {
        return Err(FingerprintError::InvalidArguments {
            reason: "empty topological torsion hash path is undefined in pinned source",
        });
    }
    let reverse_it = reverse_packed_path(path_codes);
    let mut result = 0u32;
    if reverse_it {
        // Source unsigned32 i wraps forever when len>UINTMAX. This rejects
        // only that nonterminating source branch, not huge forward vectors.
        // ROOT approved this explicit finite safety error for the source nontermination.
        if path_codes.len() > u32::MAX as usize {
            return Err(FingerprintError::InvalidArguments {
                reason: "reversed topological torsion hash exceeds source terminating u32 index range",
            });
        }
        for i in 0..path_codes.len() {
            hash_combine(&mut result, path_codes[path_codes.len() - i - 1]);
        }
    } else {
        for &path_code in path_codes {
            hash_combine(&mut result, path_code);
        }
    }
    Ok(result)
}
