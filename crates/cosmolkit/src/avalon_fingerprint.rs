//! RDKit's molecule-to-Avalon adapter; the fingerprint engine owns all hashing.
use crate::{AvalonFingerprintParams, Fingerprint, MolecularIoError, Molecule};

#[derive(Debug)]
pub enum AvalonFingerprintError {
    Input(MolecularIoError),
    Engine(cosmolkit_fingerprints::AvalonError),
}
impl std::fmt::Display for AvalonFingerprintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Input(e) => e.fmt(f),
            Self::Engine(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for AvalonFingerprintError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Input(e) => e,
            Self::Engine(e) => e,
        })
    }
}
impl Molecule {
    /// Avalon fingerprint using the C++ default 512 bits and 0x007fff flags.
    /// Input conversion uses stereo-enabled, kekulized MOL, without writing back.
    /// Missing coordinates are generated only in the temporary MOL conversion;
    /// this does not require the public `depict` capability.
    pub fn fingerprint_avalon(&self) -> Result<Fingerprint, AvalonFingerprintError> {
        self.fingerprint_avalon_with_params(&AvalonFingerprintParams::default())
    }
    pub fn fingerprint_avalon_with_params(
        &self,
        params: &AvalonFingerprintParams,
    ) -> Result<Fingerprint, AvalonFingerprintError> {
        // RDKit❗✔️:   std::string molB = MolToMolBlock(mol, true);
        // RDKit❗✔️:   Utils::LocaleSwitcher ls;
        // RDKit❗✔️:   struct reaccs_molecule_t *res = MolStr2Mol((char *)molB.c_str());
        // RDKit❗✔️:   POSTCONDITION(res, "could not build a molecule");
        // RDKit❗✔️:   return res;
        // Behavior: the existing writer owns notation conversion, including
        // source-defined temporary coordinate generation. cap-fingerprints
        // supplies that internal dependency without exposing public depiction.
        // Only detached text reaches the fingerprint owner. Locale-independent
        // Rust formatting replaces ls.
        // Complexity: one MOL conversion and one detached engine evaluation;
        // no source block checkout, writeback or whole Molecule clone.
        params.validate().map_err(AvalonFingerprintError::Engine)?;
        let mut input = self.mol_write_input();
        input.allow_coordinate_generation = true;
        let block = cosmolkit_io::write_mol_block_with_params(
            input,
            &crate::MolBlockWriteParams::default(),
        )
        .map_err(MolecularIoError::MolWrite)
        .and_then(|text| crate::molecular_io::output_string(text, "MOL"))
        .map_err(AvalonFingerprintError::Input)?;
        cosmolkit_fingerprints::avalon_fingerprint(&block, params)
            .map_err(AvalonFingerprintError::Engine)
    }
}
