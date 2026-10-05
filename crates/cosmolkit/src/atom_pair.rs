//! Canonical AtomPair call transport; chemistry remains in the fingerprint owner.
use crate::{
    Fingerprint, FingerprintAdditionalOutput, Molecule, MorganReadError, SparseBitFingerprint,
    SparseCountFingerprint, SparseCountFingerprint32,
};
use cosmolkit_fingerprints::{
    AtomPairAtomInvariantsGenerator, AtomPairCall, AtomPairError, AtomPairParams,
    AtomPairPreparedInput, atom_pair_bits, atom_pair_count, atom_pair_sparse_bits,
    atom_pair_sparse_count,
};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPairFingerprintParams {
    pub generator: AtomPairParams,
    pub from_atoms: Option<Vec<u32>>,
    pub ignore_atoms: Option<Vec<u32>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub conformer_id: i32,
    pub atom_invariants_generator: Option<AtomPairAtomInvariantsGenerator>,
    pub use_legacy_stereo_perception: bool,
}
impl Default for AtomPairFingerprintParams {
    fn default() -> Self {
        Self {
            generator: AtomPairParams::default(),
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            conformer_id: -1,
            atom_invariants_generator: None,
            use_legacy_stereo_perception: true,
        }
    }
}
#[derive(Debug)]
pub enum AtomPairReadError {
    Preparation(MorganReadError),
    Generator(AtomPairError),
}
impl fmt::Display for AtomPairReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Preparation(e) => write!(f, "AtomPair preparation failed: {e}"),
            Self::Generator(e) => write!(f, "AtomPair generation failed: {e}"),
        }
    }
}
impl std::error::Error for AtomPairReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Preparation(e) => Some(e),
            Self::Generator(e) => Some(e),
        }
    }
}
impl Molecule {
    pub fn atom_pair_fingerprint(&self) -> Result<Fingerprint, AtomPairReadError> {
        self.atom_pair_fingerprint_with_params(&AtomPairFingerprintParams::default(), None)
    }
    pub fn atom_pair_fingerprint_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<Fingerprint, AtomPairReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(AtomPairReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = AtomPairCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        atom_pair_bits(&input, &params.generator, &call, output)
            .map_err(AtomPairReadError::Generator)
    }
    pub fn atom_pair_sparse_fingerprint(&self) -> Result<SparseBitFingerprint, AtomPairReadError> {
        self.atom_pair_sparse_fingerprint_with_params(&AtomPairFingerprintParams::default(), None)
    }
    pub fn atom_pair_sparse_fingerprint_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseBitFingerprint, AtomPairReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(AtomPairReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = AtomPairCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        atom_pair_sparse_bits(&input, &params.generator, &call, output)
            .map_err(AtomPairReadError::Generator)
    }
    pub fn atom_pair_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint32, AtomPairReadError> {
        self.atom_pair_count_fingerprint_with_params(&AtomPairFingerprintParams::default(), None)
    }
    pub fn atom_pair_count_fingerprint_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, AtomPairReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(AtomPairReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = AtomPairCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        atom_pair_count(&input, &params.generator, &call, output)
            .map_err(AtomPairReadError::Generator)
    }
    pub fn atom_pair_sparse_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint, AtomPairReadError> {
        self.atom_pair_sparse_count_fingerprint_with_params(
            &AtomPairFingerprintParams::default(),
            None,
        )
    }
    pub fn atom_pair_sparse_count_fingerprint_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint, AtomPairReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(AtomPairReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = AtomPairCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        atom_pair_sparse_count(&input, &params.generator, &call, output)
            .map_err(AtomPairReadError::Generator)
    }
}
