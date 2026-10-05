//! Canonical Topological Torsion call transport; chemistry remains in the fingerprint owner.
use crate::{
    AdditionalOutput, Fingerprint, Molecule, MorganReadError, SparseBitFingerprint,
    SparseCountFingerprint, SparseCountFingerprint32,
};
use cosmolkit_fingerprints::{
    AtomPairAtomInvariantsGenerator, AtomPairPreparedInput, TopologicalTorsionCall,
    TopologicalTorsionError, TopologicalTorsionParams, topological_torsion_bits,
    topological_torsion_count, topological_torsion_sparse_bits, topological_torsion_sparse_count,
};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionFingerprintParams {
    pub generator: TopologicalTorsionParams,
    pub from_atoms: Option<Vec<u32>>,
    pub ignore_atoms: Option<Vec<u32>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub conformer_id: i32,
    pub atom_invariants_generator: Option<AtomPairAtomInvariantsGenerator>,
    pub use_legacy_stereo_perception: bool,
}
impl Default for TopologicalTorsionFingerprintParams {
    fn default() -> Self {
        Self {
            generator: TopologicalTorsionParams::default(),
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
pub enum TopologicalTorsionReadError {
    Preparation(MorganReadError),
    Generator(TopologicalTorsionError),
}
impl fmt::Display for TopologicalTorsionReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Preparation(e) => write!(f, "Topological Torsion preparation failed: {e}"),
            Self::Generator(e) => write!(f, "Topological Torsion generation failed: {e}"),
        }
    }
}
impl std::error::Error for TopologicalTorsionReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Preparation(e) => Some(e),
            Self::Generator(e) => Some(e),
        }
    }
}
impl Molecule {
    pub fn topological_torsion_fingerprint(
        &self,
    ) -> Result<Fingerprint, TopologicalTorsionReadError> {
        self.topological_torsion_fingerprint_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn topological_torsion_fingerprint_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut AdditionalOutput>,
    ) -> Result<Fingerprint, TopologicalTorsionReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = TopologicalTorsionCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        topological_torsion_bits(&input, &params.generator, &call, output)
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn topological_torsion_sparse_fingerprint(
        &self,
    ) -> Result<SparseBitFingerprint, TopologicalTorsionReadError> {
        self.topological_torsion_sparse_fingerprint_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn topological_torsion_sparse_fingerprint_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut AdditionalOutput>,
    ) -> Result<SparseBitFingerprint, TopologicalTorsionReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = TopologicalTorsionCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        topological_torsion_sparse_bits(&input, &params.generator, &call, output)
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn topological_torsion_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint32, TopologicalTorsionReadError> {
        self.topological_torsion_count_fingerprint_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn topological_torsion_count_fingerprint_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut AdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, TopologicalTorsionReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = TopologicalTorsionCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        topological_torsion_count(&input, &params.generator, &call, output)
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn topological_torsion_sparse_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        self.topological_torsion_sparse_count_fingerprint_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn topological_torsion_sparse_count_fingerprint_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut AdditionalOutput>,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        let prepared = crate::morgan::prepare_morgan_read_input(self)
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let base = prepared.owner_input();
        let input = AtomPairPreparedInput {
            topology: base.topology,
            properties: base.properties,
            coordinates: base.coordinates,
            valence: base.valence,
            rings: base.rings,
            use_legacy_stereo_perception: params.use_legacy_stereo_perception,
        };
        let call = TopologicalTorsionCall {
            from_atoms: params.from_atoms.as_deref(),
            ignore_atoms: params.ignore_atoms.as_deref(),
            custom_atom_invariants: params.custom_atom_invariants.as_deref(),
            custom_bond_invariants: params.custom_bond_invariants.as_deref(),
            conformer_id: params.conformer_id,
            atom_invariants_generator: params.atom_invariants_generator,
        };
        topological_torsion_sparse_count(&input, &params.generator, &call, output)
            .map_err(TopologicalTorsionReadError::Generator)
    }
}
