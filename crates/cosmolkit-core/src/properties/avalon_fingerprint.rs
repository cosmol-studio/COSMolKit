//! Fail-closed boundary for Avalon fingerprints.
//!
//! COSMolKit does not currently contain a source-level port of the Avalon
//! fingerprint implementation. Returning a local path hash under this API
//! would produce plausible-looking but non-equivalent bits, so the public
//! entry point reports a structured unsupported error until the full source
//! dependency and exact-bit parity surface are implemented.

use crate::{Fingerprint, FingerprintError, Molecule};

/// Parameters reserved for the future source-backed Avalon fingerprint port.
#[derive(Debug, Clone, PartialEq)]
pub struct AvalonFingerprintParams {
    pub min_path: u32,
    pub max_path: u32,
    pub n_bits: usize,
    pub n_bits_per_hash: u32,
    pub use_bond_order: bool,
    pub use_hs: bool,
    pub tautomeric_fingerprint: bool,
    pub from_atoms: Option<Vec<usize>>,
}

impl Default for AvalonFingerprintParams {
    fn default() -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            n_bits: 2048,
            n_bits_per_hash: 1,
            use_bond_order: true,
            use_hs: false,
            tautomeric_fingerprint: false,
            from_atoms: None,
        }
    }
}

pub fn avalon_fingerprint(
    _molecule: &Molecule,
    _params: &AvalonFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    // RDKit source: External/AvalonTools/AvalonTools.h
    // RDKit❌❌: RDKIT_AVALONLIB_EXPORT void getAvalonFP(const RDKit::ROMol &mol,
    // RDKit❌❌:                                         ExplicitBitVect &res,
    // RDKit❌❌:                                         unsigned int nBits = 512,
    // RDKit❌❌:                                         bool isQuery = false,
    // RDKit❌❌:                                         bool resetVect = true,
    // RDKit❌❌:                                         unsigned int bitFlags = avalonSSSBits);
    // RDKit source: External/AvalonTools/AvalonTools.cpp
    // RDKit❌❌: unsigned int nBytes = nBits / 8;
    // RDKit❌❌: struct reaccs_molecule_t *mp = molToReaccs(mol);
    // RDKit❌❌: reaccsToFingerprint(mp, res, bitFlags, isQuery, resetVect, nBytes);
    // RDKit❌❌: FreeMolecule(mp);
    Err(FingerprintError::UnsupportedOption {
        option: "avalon_fingerprint",
        reason: "RDKit Avalon exact-bit source port is not implemented",
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn avalon_fingerprint_fails_closed_until_exact_port_exists() {
        let err = avalon_fingerprint(&Molecule::new(), &AvalonFingerprintParams::default())
            .expect_err("unfinished Avalon fingerprint must not return approximate bits");
        assert!(matches!(
            err,
            FingerprintError::UnsupportedOption {
                option: "avalon_fingerprint",
                ..
            }
        ));
    }
}
