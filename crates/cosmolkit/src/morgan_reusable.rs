//! Thin canonical reusable Morgan state and call transport.
//! Fingerprint chemistry and worker scheduling remain in their sole domain owner.
use crate::{
    AtomPairAtomInvariantsGenerator, Fingerprint, FingerprintAdditionalOutput, Molecule,
    MorganGeneratorParams, MorganReadError, QueryGraph, SparseBitFingerprint,
    SparseCountFingerprint, SparseCountFingerprint32,
};
use cosmolkit_fingerprints::{MorganAtomProvider, MorganBondProvider, MorganCall, MorganOperator};

/// Immutable, independently captured atom-invariant provider configuration.
#[derive(Debug, Clone)]
pub struct MorganAtomInvariantsGenerator {
    inner: MorganAtomProvider,
}
impl MorganAtomInvariantsGenerator {
    #[must_use]
    pub fn connectivity(include_ring_membership: bool) -> Self {
        Self {
            inner: MorganAtomProvider::Connectivity {
                include_ring_membership,
            },
        }
    }
    #[must_use]
    pub fn features(patterns: Option<Vec<QueryGraph>>) -> Self {
        Self {
            inner: MorganAtomProvider::Features { patterns },
        }
    }
    #[must_use]
    pub fn atom_pair(generator: AtomPairAtomInvariantsGenerator) -> Self {
        Self {
            inner: MorganAtomProvider::AtomPair(generator),
        }
    }
}
/// Explicit bond-provider inputs, captured independently of live settings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar, factory)]
pub struct MorganBondInvariantsGenerator {
    pub use_bond_types: bool,
    pub include_chirality: bool,
}
impl MorganBondInvariantsGenerator {
    #[must_use]
    pub fn new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self {
            use_bond_types,
            include_chirality,
        }
    }
    #[must_use]
    pub fn use_bond_types(&self) -> bool {
        self.use_bond_types
    }
    #[must_use]
    pub fn include_chirality(&self) -> bool {
        self.include_chirality
    }
}
impl Default for MorganBondInvariantsGenerator {
    fn default() -> Self {
        Self::new(true, false)
    }
}
/// Shared persistent generator. Clones and settings views alias the same state.
#[derive(Debug, Clone)]
pub struct MorganFingerprintGenerator {
    inner: MorganOperator,
}
/// Live bound source argument view; this view retains its generator's lifetime.
#[derive(Debug, Clone)]
pub struct MorganSettings {
    inner: cosmolkit_fingerprints::MorganSettings,
}
/// Owned call selections. None and a present empty vector remain distinct.
#[derive(Debug, Clone, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(factory)]
pub struct MorganCallParams {
    pub from_atoms: Option<Vec<u32>>,
    pub ignore_atoms: Option<Vec<u32>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
}
impl MorganCallParams {
    #[must_use]
    pub fn new(
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
    ) -> Self {
        Self {
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
        }
    }
    fn owner_call(&self) -> MorganCall<'_> {
        MorganCall {
            from_atoms: self.from_atoms.as_deref(),
            ignore_atoms: self.ignore_atoms.as_deref(),
            custom_atom_invariants: self.custom_atom_invariants.as_deref(),
            custom_bond_invariants: self.custom_bond_invariants.as_deref(),
            conformer_id: -1,
        }
    }
}
impl Default for MorganCallParams {
    fn default() -> Self {
        Self::new(None, None, None, None)
    }
}
impl MorganFingerprintGenerator {
    pub fn new(
        params: Option<&MorganGeneratorParams>,
        atom_invariants: Option<&MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&MorganBondInvariantsGenerator>,
    ) -> Result<Self, MorganReadError> {
        let defaults = MorganGeneratorParams::default();
        // One creation-time source getCopy equivalent for an explicit provider.
        // The owner consumes this independent vector; calculations only borrow it.
        MorganOperator::new(
            params.unwrap_or(&defaults),
            atom_invariants.map(|provider| provider.inner.clone()),
            bond_invariants.map(|provider| MorganBondProvider {
                use_bond_types: provider.use_bond_types,
                include_chirality: provider.include_chirality,
            }),
        )
        .map(|inner| Self { inner })
        .map_err(MorganReadError::Generator)
    }
    pub fn from_json(json: &str) -> Result<Self, MorganReadError> {
        MorganOperator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(MorganReadError::Generator)
    }
    #[must_use]
    pub fn settings(&self) -> MorganSettings {
        MorganSettings {
            inner: self.inner.settings(),
        }
    }
    pub fn info_string(&self) -> Result<String, MorganReadError> {
        self.inner.info_string().map_err(MorganReadError::Generator)
    }
    pub fn to_json(&self) -> Result<cosmolkit_model::PropertyText, MorganReadError> {
        self.inner.to_json().map_err(MorganReadError::Generator)
    }
}
macro_rules! live_field {
    ($get:ident,$set:ident,$ty:ty) => {
        pub fn $get(&self) -> Result<$ty, MorganReadError> {
            self.inner.$get().map_err(MorganReadError::Generator)
        }
        pub fn $set(&mut self, value: $ty) -> Result<(), MorganReadError> {
            self.inner.$set(value).map_err(MorganReadError::Generator)
        }
    };
}
impl MorganSettings {
    live_field!(radius, set_radius, u32);
    live_field!(only_nonzero_invariants, set_only_nonzero_invariants, bool);
    live_field!(
        include_redundant_environments,
        set_include_redundant_environments,
        bool
    );
    live_field!(include_chirality, set_include_chirality, bool);
    live_field!(count_simulation, set_count_simulation, bool);
    live_field!(fp_size, set_fp_size, u32);
    live_field!(bits_per_feature, set_bits_per_feature, u32);
    live_field!(count_bounds, set_count_bounds, Vec<u32>);
    pub fn params(&self) -> Result<MorganGeneratorParams, MorganReadError> {
        self.inner.snapshot().map_err(MorganReadError::Generator)
    }
}
macro_rules! calculations {
    ($bulk:ident,$method:ident,$owner:ident,$result:ty,$batch:ty) => {
        impl MorganFingerprintGenerator {
            pub fn $bulk(
                &self,
                molecules: &[Option<&Molecule>],
                num_threads: i32,
            ) -> Result<$batch, MorganReadError> {
                // Keep one owner lock across needs decisions and calculation.
                // On poison, preserve the original ordered preparation errors
                // before surfacing the owner error; never prescan valence.
                let calculation = match self.inner.lock_for_calculation() {
                    Ok(calculation) => calculation,
                    Err(error) => {
                        for molecule in molecules.iter().flatten() {
                            crate::morgan::prepare_morgan_read_input(molecule)
                                .map_err(MorganReadError::Preparation)?;
                        }
                        return Err(MorganReadError::Generator(error));
                    }
                };
                let call = MorganCall::default();
                let prepared = molecules
                    .iter()
                    .map(|m| {
                        m.map(|molecule| {
                            crate::morgan::prepare_morgan_read_input_with_rings(
                                molecule,
                                calculation.needs_ring_info(molecule.properties(), &call),
                            )
                        })
                        .transpose()
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(MorganReadError::Preparation)?;
                let inputs = prepared
                    .iter()
                    .map(|p| p.as_ref().map(|p| p.owner_input()))
                    .collect::<Vec<_>>();
                calculation
                    .$bulk(&inputs, num_threads)
                    .map(<$batch>::from)
                    .map_err(MorganReadError::Generator)
            }
        }
        impl Molecule {
            pub fn $method(
                &self,
                generator: &MorganFingerprintGenerator,
                params: Option<&MorganCallParams>,
                additional_output: Option<&mut FingerprintAdditionalOutput>,
            ) -> Result<$result, MorganReadError> {
                let defaults = MorganCallParams::default();
                let params = params.unwrap_or(&defaults);
                let call = params.owner_call();
                let calculation = match generator.inner.lock_for_calculation() {
                    Ok(calculation) => calculation,
                    Err(error) => {
                        crate::morgan::prepare_morgan_read_input(self)
                            .map_err(MorganReadError::Preparation)?;
                        return Err(MorganReadError::Generator(error));
                    }
                };
                let prepared = crate::morgan::prepare_morgan_read_input_with_rings(
                    self,
                    calculation.needs_ring_info(self.properties(), &call),
                )
                .map_err(MorganReadError::Preparation)?;
                calculation
                    .$owner(&prepared.owner_input(), &call, additional_output)
                    .map_err(MorganReadError::Generator)
            }
        }
    };
}
calculations!(
    fingerprints,
    fingerprint_morgan_with_generator,
    bits,
    Fingerprint,
    crate::FingerprintBatch
);
calculations!(
    counts,
    fingerprint_morgan_count_with_generator,
    count,
    SparseCountFingerprint32,
    Vec<Option<SparseCountFingerprint32>>
);
calculations!(
    sparse_fingerprints,
    fingerprint_morgan_sparse_with_generator,
    sparse_bits,
    SparseBitFingerprint,
    Vec<Option<SparseBitFingerprint>>
);
calculations!(
    sparse_counts,
    fingerprint_morgan_sparse_count_with_generator,
    sparse_count,
    SparseCountFingerprint,
    Vec<Option<SparseCountFingerprint>>
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FingerprintPreparationError, MorganFingerprintParams, SmilesParseParams};

    fn cold(smiles: &str, done: Option<i32>) -> Molecule {
        let mut record = cosmolkit_smiles::parse_smiles(
            smiles,
            &SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        if let Some(value) = done {
            record
                .properties
                .set_prop("_StereochemDone", value)
                .unwrap();
        }
        Molecule::from_parts(record.topology, record.coordinates, record.properties)
            .unwrap()
            .with_assigned_valence()
            .unwrap()
    }

    fn acquisitions() -> usize {
        crate::morgan::ring_preparation_probe::acquisitions()
    }

    #[test]
    fn morgan_locked_ring_free_calls_and_bulk_keep_cold_cache() {
        let params = MorganGeneratorParams {
            radius: 2,
            // Explicit captured provider below must override this factory flag.
            include_ring_membership: true,
            ..Default::default()
        };
        let provider = MorganAtomInvariantsGenerator::connectivity(false);
        let generator =
            MorganFingerprintGenerator::new(Some(&params), Some(&provider), None).unwrap();
        let cold = cold("C1CCCCC1", None);
        let before = acquisitions();
        let configured = MorganFingerprintParams {
            generator: MorganGeneratorParams {
                radius: 2,
                include_ring_membership: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let expected = cold
            .fingerprint_morgan_with_params(&configured, None)
            .unwrap();
        for _ in 0..2 {
            assert_eq!(
                cold.fingerprint_morgan_with_generator(&generator, None, None)
                    .unwrap(),
                expected
            );
            assert_eq!(
                generator
                    .fingerprints(&[Some(&cold), None, Some(&cold)], 2)
                    .unwrap(),
                crate::FingerprintBatch::new(vec![
                    Some(expected.clone()),
                    None,
                    Some(expected.clone())
                ])
            );
        }
        assert_eq!(acquisitions(), before);
        assert!(cold.derived_cache_runtime().valid_ring_info().is_none());

        // Custom invariants bypass even a feature provider, whose own context
        // would require rings. A present zero done marker suppresses chirality
        // preparation just as source hasProp(), despite live settings changes.
        let features = MorganAtomInvariantsGenerator::features(Some(vec![]));
        let generator =
            MorganFingerprintGenerator::new(Some(&params), Some(&features), None).unwrap();
        let mut settings = generator.settings();
        settings.set_include_chirality(true).unwrap();
        let done_zero = self::cold("CCO", Some(0));
        let custom = MorganCallParams {
            custom_atom_invariants: Some(vec![11, 12, 13]),
            ..Default::default()
        };
        let baseline = acquisitions();
        for _ in 0..2 {
            let result = done_zero
                .fingerprint_morgan_sparse_count_with_generator(&generator, Some(&custom), None)
                .unwrap();
            assert!(result.nonzero_elements().contains_key(&11));
        }
        assert_eq!(acquisitions(), baseline);
        assert!(
            done_zero
                .derived_cache_runtime()
                .valid_ring_info()
                .is_none()
        );
    }

    #[test]
    fn morgan_locked_required_rings_reuse_warm_and_prepare_cold() {
        let generator = MorganFingerprintGenerator::new(
            Some(&MorganGeneratorParams {
                radius: 2,
                ..Default::default()
            }),
            None,
            None,
        )
        .unwrap();
        let cold = cold("C1CCCCC1", None);
        let warm = Molecule::from_smiles("C1CCCCC1").unwrap();
        let before = acquisitions();
        let expected = warm
            .fingerprint_morgan_with_generator(&generator, None, None)
            .unwrap();
        assert_eq!(acquisitions(), before);
        for _ in 0..2 {
            assert_eq!(
                cold.fingerprint_morgan_with_generator(&generator, None, None)
                    .unwrap(),
                expected
            );
        }
        assert_eq!(acquisitions() - before, 2);
        assert!(cold.derived_cache_runtime().valid_ring_info().is_none());
    }

    #[test]
    fn morgan_locked_preparation_errors_precede_poison_and_output_reset() {
        let generator = MorganFingerprintGenerator::new(None, None, None).unwrap();
        let valued = cold("CCO", None);
        let raw = Molecule::from_smiles_with_params(
            "CCO",
            &SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        let mut output = FingerprintAdditionalOutput::new();
        output.allocate_atom_counts();
        valued
            .fingerprint_morgan_with_generator(&generator, None, Some(&mut output))
            .unwrap();
        let previous = output.atom_counts().unwrap().to_vec();
        let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _calculation = generator.inner.lock_for_calculation().unwrap();
            panic!("fixed test poison");
        }));
        assert!(caught.is_err());
        assert!(matches!(
            raw.fingerprint_morgan_with_generator(&generator, None, Some(&mut output)),
            Err(MorganReadError::Preparation(
                FingerprintPreparationError::MissingPreparedValence
            ))
        ));
        assert_eq!(output.atom_counts().unwrap(), previous);
        let before = acquisitions();
        assert!(matches!(
            generator.fingerprints(&[Some(&valued), None, Some(&raw)], 0),
            Err(MorganReadError::Preparation(
                FingerprintPreparationError::MissingPreparedValence
            ))
        ));
        assert_eq!(
            acquisitions() - before,
            1,
            "prepare ordered items before poison"
        );
        assert!(matches!(
            valued.fingerprint_morgan_with_generator(&generator, None, Some(&mut output)),
            Err(MorganReadError::Generator(
                cosmolkit_fingerprints::MorganError::StatePoisoned
            ))
        ));
        assert_eq!(output.atom_counts().unwrap(), previous);
    }
}
