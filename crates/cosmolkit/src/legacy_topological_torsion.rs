//! Thin canonical transport for the three original legacy torsion algorithms.
use crate::{Fingerprint, Molecule, SparseCountFingerprint, TopologicalTorsionReadError};
use cosmolkit_fingerprints::{
    AtomPairPreparedInput, TopologicalTorsionCall, TopologicalTorsionError,
};

/// Immutable call configuration preserving source optional-vector semantics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyTopologicalTorsionParams {
    pub torsion_atom_count: u32,
    pub include_chirality: bool,
    pub fp_size: u32,
    pub bits_per_entry: u32,
    pub from_atoms: Option<Vec<u32>>,
    pub ignore_atoms: Option<Vec<u32>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
}
impl Default for LegacyTopologicalTorsionParams {
    fn default() -> Self {
        Self::new(4, false, 2048, 4, None, None, None)
    }
}
impl LegacyTopologicalTorsionParams {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        torsion_atom_count: u32,
        include_chirality: bool,
        fp_size: u32,
        bits_per_entry: u32,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> Self {
        Self {
            torsion_atom_count,
            include_chirality,
            fp_size,
            bits_per_entry,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
        }
    }
    fn owner_params(&self) -> cosmolkit_fingerprints::LegacyTopologicalTorsionParams {
        cosmolkit_fingerprints::LegacyTopologicalTorsionParams {
            torsion_atom_count: self.torsion_atom_count,
            include_chirality: self.include_chirality,
            fp_size: self.fp_size,
            bits_per_entry: self.bits_per_entry,
        }
    }
    fn owner_call(&self) -> TopologicalTorsionCall<'_> {
        TopologicalTorsionCall {
            from_atoms: self.from_atoms.as_deref(),
            ignore_atoms: self.ignore_atoms.as_deref(),
            custom_atom_invariants: self.custom_atom_invariants.as_deref(),
            ..Default::default()
        }
    }
}
fn with_legacy_input<T>(
    molecule: &Molecule,
    params: &LegacyTopologicalTorsionParams,
    consume: impl FnOnce(
        &AtomPairPreparedInput<'_>,
        &cosmolkit_fingerprints::LegacyTopologicalTorsionParams,
        &TopologicalTorsionCall<'_>,
    ) -> Result<T, TopologicalTorsionError>,
) -> Result<T, TopologicalTorsionReadError> {
    let prepared = crate::morgan::prepare_morgan_read_input(molecule)
        .map_err(TopologicalTorsionReadError::Preparation)?;
    let base = prepared.owner_input();
    let input = AtomPairPreparedInput {
        topology: base.topology,
        properties: base.properties,
        coordinates: base.coordinates,
        valence: base.valence,
        rings: base.rings,
        use_legacy_stereo_perception: true,
    };
    consume(&input, &params.owner_params(), &params.owner_call())
        .map_err(TopologicalTorsionReadError::Generator)
}
impl Molecule {
    pub fn legacy_topological_torsion_sparse_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        self.legacy_topological_torsion_sparse_count_fingerprint_with_params(&Default::default())
    }
    pub fn legacy_topological_torsion_sparse_count_fingerprint_with_params(
        &self,
        params: &LegacyTopologicalTorsionParams,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        with_legacy_input(
            self,
            params,
            cosmolkit_fingerprints::legacy_topological_torsion_sparse_count,
        )
    }
    pub fn legacy_topological_torsion_count_fingerprint(
        &self,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        self.legacy_topological_torsion_count_fingerprint_with_params(&Default::default())
    }
    pub fn legacy_topological_torsion_count_fingerprint_with_params(
        &self,
        params: &LegacyTopologicalTorsionParams,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        with_legacy_input(
            self,
            params,
            cosmolkit_fingerprints::legacy_topological_torsion_count,
        )
    }
    pub fn legacy_topological_torsion_fingerprint(
        &self,
    ) -> Result<Fingerprint, TopologicalTorsionReadError> {
        self.legacy_topological_torsion_fingerprint_with_params(&Default::default())
    }
    pub fn legacy_topological_torsion_fingerprint_with_params(
        &self,
        params: &LegacyTopologicalTorsionParams,
    ) -> Result<Fingerprint, TopologicalTorsionReadError> {
        with_legacy_input(
            self,
            params,
            cosmolkit_fingerprints::legacy_topological_torsion_bits,
        )
    }
}
