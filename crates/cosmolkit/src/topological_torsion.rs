//! Canonical Topological Torsion call transport; chemistry remains in the fingerprint owner.
use crate::{
    Fingerprint, FingerprintAdditionalOutput, FingerprintPreparationError, Molecule,
    SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
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
    Preparation(FingerprintPreparationError),
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
    pub fn fingerprint_topological_torsion(
        &self,
    ) -> Result<Fingerprint, TopologicalTorsionReadError> {
        self.fingerprint_topological_torsion_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn fingerprint_topological_torsion_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
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
    pub fn fingerprint_topological_torsion_sparse(
        &self,
    ) -> Result<SparseBitFingerprint, TopologicalTorsionReadError> {
        self.fingerprint_topological_torsion_sparse_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn fingerprint_topological_torsion_sparse_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
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
    pub fn fingerprint_topological_torsion_count(
        &self,
    ) -> Result<SparseCountFingerprint32, TopologicalTorsionReadError> {
        self.fingerprint_topological_torsion_count_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn fingerprint_topological_torsion_count_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
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
    pub fn fingerprint_topological_torsion_sparse_count(
        &self,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        self.fingerprint_topological_torsion_sparse_count_with_params(
            &TopologicalTorsionFingerprintParams::default(),
            None,
        )
    }
    pub fn fingerprint_topological_torsion_sparse_count_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        output: Option<&mut FingerprintAdditionalOutput>,
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

/// Reusable canonical operator; the fingerprint owner holds its sole live state.
#[derive(Debug, Clone)]
pub struct TopologicalTorsionFingerprintGenerator {
    inner: cosmolkit_fingerprints::TopologicalTorsionGenerator,
}
/// Bound view of the reusable generator's mutable source options.
#[derive(Debug, Clone)]
pub struct TopologicalTorsionSettings {
    inner: cosmolkit_fingerprints::TopologicalTorsionSettings,
}
/// Immutable call selections independent of reusable generator configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologicalTorsionCallParams {
    pub from_atoms: Option<Vec<u32>>,
    pub ignore_atoms: Option<Vec<u32>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub conformer_id: i32,
    pub use_legacy_stereo_perception: bool,
}
impl Default for TopologicalTorsionCallParams {
    fn default() -> Self {
        Self {
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            conformer_id: -1,
            use_legacy_stereo_perception: true,
        }
    }
}
impl TopologicalTorsionCallParams {
    fn owner_call(&self) -> TopologicalTorsionCall<'_> {
        TopologicalTorsionCall {
            from_atoms: self.from_atoms.as_deref(),
            ignore_atoms: self.ignore_atoms.as_deref(),
            custom_atom_invariants: self.custom_atom_invariants.as_deref(),
            custom_bond_invariants: self.custom_bond_invariants.as_deref(),
            conformer_id: self.conformer_id,
            atom_invariants_generator: None,
        }
    }
}
impl TopologicalTorsionFingerprintGenerator {
    pub fn new(
        params: Option<&TopologicalTorsionParams>,
        atom_invariants: Option<AtomPairAtomInvariantsGenerator>,
    ) -> Result<Self, TopologicalTorsionReadError> {
        let defaults = TopologicalTorsionParams::default();
        let params = params.unwrap_or(&defaults);
        cosmolkit_fingerprints::TopologicalTorsionGenerator::new(params, atom_invariants)
            .map(|inner| Self { inner })
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn from_json(json: &str) -> Result<Self, TopologicalTorsionReadError> {
        cosmolkit_fingerprints::TopologicalTorsionGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn settings(&self) -> TopologicalTorsionSettings {
        TopologicalTorsionSettings {
            inner: self.inner.settings(),
        }
    }
    pub fn info_string(&self) -> Result<String, TopologicalTorsionReadError> {
        self.inner
            .info_string()
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn to_json(&self) -> Result<String, TopologicalTorsionReadError> {
        self.inner
            .to_json()
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn fingerprints(
        &self,
        molecules: &[Option<&Molecule>],
        num_threads: i32,
    ) -> Result<Vec<Option<Fingerprint>>, TopologicalTorsionReadError> {
        let prepared = molecules
            .iter()
            .map(|molecule| {
                molecule
                    .map(crate::morgan::prepare_morgan_read_input)
                    .transpose()
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let inputs = prepared
            .iter()
            .map(|p| {
                p.as_ref().map(|p| {
                    let base = p.owner_input();
                    AtomPairPreparedInput {
                        topology: base.topology,
                        properties: base.properties,
                        coordinates: base.coordinates,
                        valence: base.valence,
                        rings: base.rings,
                        use_legacy_stereo_perception: true,
                    }
                })
            })
            .collect::<Vec<_>>();
        self.inner
            .fingerprints(&inputs, num_threads)
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn sparse_fingerprints(
        &self,
        molecules: &[Option<&Molecule>],
        num_threads: i32,
    ) -> Result<Vec<Option<SparseBitFingerprint>>, TopologicalTorsionReadError> {
        let prepared = molecules
            .iter()
            .map(|molecule| {
                molecule
                    .map(crate::morgan::prepare_morgan_read_input)
                    .transpose()
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let inputs = prepared
            .iter()
            .map(|p| {
                p.as_ref().map(|p| {
                    let base = p.owner_input();
                    AtomPairPreparedInput {
                        topology: base.topology,
                        properties: base.properties,
                        coordinates: base.coordinates,
                        valence: base.valence,
                        rings: base.rings,
                        use_legacy_stereo_perception: true,
                    }
                })
            })
            .collect::<Vec<_>>();
        self.inner
            .sparse_fingerprints(&inputs, num_threads)
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn counts(
        &self,
        molecules: &[Option<&Molecule>],
        num_threads: i32,
    ) -> Result<Vec<Option<SparseCountFingerprint32>>, TopologicalTorsionReadError> {
        let prepared = molecules
            .iter()
            .map(|molecule| {
                molecule
                    .map(crate::morgan::prepare_morgan_read_input)
                    .transpose()
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let inputs = prepared
            .iter()
            .map(|p| {
                p.as_ref().map(|p| {
                    let base = p.owner_input();
                    AtomPairPreparedInput {
                        topology: base.topology,
                        properties: base.properties,
                        coordinates: base.coordinates,
                        valence: base.valence,
                        rings: base.rings,
                        use_legacy_stereo_perception: true,
                    }
                })
            })
            .collect::<Vec<_>>();
        self.inner
            .counts(&inputs, num_threads)
            .map_err(TopologicalTorsionReadError::Generator)
    }
    pub fn sparse_counts(
        &self,
        molecules: &[Option<&Molecule>],
        num_threads: i32,
    ) -> Result<Vec<Option<SparseCountFingerprint>>, TopologicalTorsionReadError> {
        let prepared = molecules
            .iter()
            .map(|molecule| {
                molecule
                    .map(crate::morgan::prepare_morgan_read_input)
                    .transpose()
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(TopologicalTorsionReadError::Preparation)?;
        let inputs = prepared
            .iter()
            .map(|p| {
                p.as_ref().map(|p| {
                    let base = p.owner_input();
                    AtomPairPreparedInput {
                        topology: base.topology,
                        properties: base.properties,
                        coordinates: base.coordinates,
                        valence: base.valence,
                        rings: base.rings,
                        use_legacy_stereo_perception: true,
                    }
                })
            })
            .collect::<Vec<_>>();
        self.inner
            .sparse_counts(&inputs, num_threads)
            .map_err(TopologicalTorsionReadError::Generator)
    }
}
macro_rules! settings_field {
    ($get:ident,$set:ident,$ty:ty) => {
        pub fn $get(&self) -> Result<$ty, TopologicalTorsionReadError> {
            self.inner
                .$get()
                .map_err(TopologicalTorsionReadError::Generator)
        }
        pub fn $set(&mut self, value: $ty) -> Result<(), TopologicalTorsionReadError> {
            self.inner
                .$set(value)
                .map_err(TopologicalTorsionReadError::Generator)
        }
    };
}
impl TopologicalTorsionSettings {
    settings_field!(torsion_atom_count, set_torsion_atom_count, u32);
    settings_field!(only_shortest_paths, set_only_shortest_paths, bool);
    settings_field!(include_chirality, set_include_chirality, bool);
    settings_field!(count_simulation, set_count_simulation, bool);
    settings_field!(fp_size, set_fp_size, u32);
    settings_field!(bits_per_feature, set_bits_per_feature, u32);
    settings_field!(count_bounds, set_count_bounds, Vec<u32>);
    pub fn params(&self) -> Result<TopologicalTorsionParams, TopologicalTorsionReadError> {
        self.inner
            .snapshot()
            .map_err(TopologicalTorsionReadError::Generator)
    }
}
fn with_reusable_input<T>(
    molecule: &Molecule,
    params: &TopologicalTorsionCallParams,
    consumer: impl FnOnce(
        &AtomPairPreparedInput<'_>,
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
        use_legacy_stereo_perception: params.use_legacy_stereo_perception,
    };
    consumer(&input, &params.owner_call()).map_err(TopologicalTorsionReadError::Generator)
}
impl Molecule {
    pub fn fingerprint_topological_torsion_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<Fingerprint, TopologicalTorsionReadError> {
        let defaults = TopologicalTorsionCallParams::default();
        let params = params.unwrap_or(&defaults);
        with_reusable_input(self, params, |input, call| {
            generator.inner.bits(input, call, output)
        })
    }
    pub fn fingerprint_topological_torsion_sparse_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseBitFingerprint, TopologicalTorsionReadError> {
        let defaults = TopologicalTorsionCallParams::default();
        let params = params.unwrap_or(&defaults);
        with_reusable_input(self, params, |input, call| {
            generator.inner.sparse_bits(input, call, output)
        })
    }
    pub fn fingerprint_topological_torsion_count_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, TopologicalTorsionReadError> {
        let defaults = TopologicalTorsionCallParams::default();
        let params = params.unwrap_or(&defaults);
        with_reusable_input(self, params, |input, call| {
            generator.inner.count(input, call, output)
        })
    }
    pub fn fingerprint_topological_torsion_sparse_count_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint, TopologicalTorsionReadError> {
        let defaults = TopologicalTorsionCallParams::default();
        let params = params.unwrap_or(&defaults);
        with_reusable_input(self, params, |input, call| {
            generator.inner.sparse_count(input, call, output)
        })
    }
}
