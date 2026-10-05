//! Pinned RDKit351f8f378f8ad6bbd517980c38896e66bf907af8 AtomPairs vocabulary.
//! Values are exposed through the canonical facade; the existing atom code table
//! remains the unique owner of atom-type order and its implicit trailing zero.

/// Namespace for source-defined atom-pair encoding constants, without live state.
#[derive(Debug)]
pub struct AtomPairsParameters {
    _private: (),
}
impl AtomPairsParameters {
    pub const fn version() -> &'static str {
        // RDKit✔️✔️: const std::string atomPairsVersion = "1.1.0";
        // Constant evaluation, O(1), no allocation or state access.
        "1.1.0"
    }
    pub const fn num_type_bits() -> u32 {
        // RDKit✔️✔️: const unsigned int numTypeBits = 4;
        // Constant evaluation, O(1), no allocation or state access.
        4
    }
    pub const fn num_pi_bits() -> u32 {
        // RDKit✔️✔️: const unsigned int numPiBits = 2;
        // Constant evaluation, O(1), no allocation or state access.
        2
    }
    pub const fn num_branch_bits() -> u32 {
        // RDKit✔️✔️: const unsigned int numBranchBits = 3;
        // Constant evaluation, O(1), no allocation or state access.
        3
    }
    pub const fn num_chiral_bits() -> u32 {
        // RDKit✔️✔️: const unsigned int numChiralBits = 2;
        // Constant evaluation, O(1), no allocation or state access.
        2
    }
    pub const fn code_size() -> u32 {
        // RDKit✔️✔️: const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
        // Constant evaluation, O(1), no allocation or state access.
        Self::num_type_bits() + Self::num_pi_bits() + Self::num_branch_bits()
    }
    pub const fn num_path_bits() -> u32 {
        // RDKit✔️✔️: const unsigned int numPathBits = 5;
        // Constant evaluation, O(1), no allocation or state access.
        5
    }
    pub const fn max_path_length() -> u32 {
        // RDKit✔️✔️: const unsigned int maxPathLen = (1 << numPathBits) - 1;
        // Constant evaluation, O(1), no allocation or state access.
        (1 << Self::num_path_bits()) - 1
    }
    pub const fn num_atom_pair_fingerprint_bits() -> u32 {
        // RDKit✔️✔️: const unsigned int numAtomPairFingerprintBits =
        // RDKit✔️✔️:     numPathBits + 2 * codeSize;  // note that this is only accurate if chirality
        // RDKit✔️✔️:                                  // is not included
        // Constant evaluation, O(1), no allocation or state access.
        Self::num_path_bits() + 2 * Self::code_size()
    }
    pub fn atom_types() -> Vec<u32> {
        // RDKit✔️✔️: std::vector<unsigned int> atomPairTypes(
        // RDKit✔️✔️:     RDKit::AtomPairs::atomNumberTypes,
        // RDKit✔️✔️:     RDKit::AtomPairs::atomNumberTypes +
        // RDKit✔️✔️:         sizeof(RDKit::AtomPairs::atomNumberTypes) / sizeof(unsigned int));
        // Copies the unique complete 16-entry array, preserving implicit zero.
        // Like the source vector this is one O(16) owned allocation per read.
        crate::atom_code::ATOM_NUMBER_TYPES.to_vec()
    }
}
