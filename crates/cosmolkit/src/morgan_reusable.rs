//! Thin canonical reusable Morgan state and call transport.
//! Fingerprint chemistry and worker scheduling remain in their sole domain owner.
use crate::{
    AtomPairAtomInvariantsGenerator, Fingerprint, FingerprintAdditionalOutput, Molecule,
    MorganParams, MorganReadError, QueryGraph, SparseBitFingerprint, SparseCountFingerprint,
    SparseCountFingerprint32,
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
/// Immutable explicit bond-provider flags, captured independently of live settings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MorganBondInvariantsGenerator {
    inner: MorganBondProvider,
}
impl MorganBondInvariantsGenerator {
    #[must_use]
    pub fn new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self {
            inner: MorganBondProvider {
                use_bond_types,
                include_chirality,
            },
        }
    }
    #[must_use]
    pub fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types
    }
    #[must_use]
    pub fn include_chirality(&self) -> bool {
        self.inner.include_chirality
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganCallParams {
    pub from_atoms: Option<Vec<u32>>,
    pub ignore_atoms: Option<Vec<u32>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub conformer_id: i32,
}
impl MorganCallParams {
    #[must_use]
    pub fn new(
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
    ) -> Self {
        Self {
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id,
        }
    }
    fn owner_call(&self) -> MorganCall<'_> {
        MorganCall {
            from_atoms: self.from_atoms.as_deref(),
            ignore_atoms: self.ignore_atoms.as_deref(),
            custom_atom_invariants: self.custom_atom_invariants.as_deref(),
            custom_bond_invariants: self.custom_bond_invariants.as_deref(),
            conformer_id: self.conformer_id,
        }
    }
}
impl Default for MorganCallParams {
    fn default() -> Self {
        Self::new(None, None, None, None, -1)
    }
}
impl MorganFingerprintGenerator {
    pub fn new(
        params: Option<&MorganParams>,
        atom_invariants: Option<&MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&MorganBondInvariantsGenerator>,
    ) -> Result<Self, MorganReadError> {
        let defaults = MorganParams::default();
        // One creation-time source getCopy equivalent for an explicit provider.
        // The owner consumes this independent vector; calculations only borrow it.
        MorganOperator::new(
            params.unwrap_or(&defaults),
            atom_invariants.map(|provider| provider.inner.clone()),
            bond_invariants.map(|provider| provider.inner),
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
    pub fn params(&self) -> Result<MorganParams, MorganReadError> {
        self.inner.snapshot().map_err(MorganReadError::Generator)
    }
}
macro_rules! calculations {
    ($bulk:ident,$method:ident,$owner:ident,$result:ty) => {
        impl MorganFingerprintGenerator {
            pub fn $bulk(
                &self,
                molecules: &[Option<&Molecule>],
                num_threads: i32,
            ) -> Result<Vec<Option<$result>>, MorganReadError> {
                // Shared ROOT preparation boundary, with no molecule/state clones.
                let prepared = molecules
                    .iter()
                    .map(|m| m.map(crate::morgan::prepare_morgan_read_input).transpose())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(MorganReadError::Preparation)?;
                let inputs = prepared
                    .iter()
                    .map(|p| p.as_ref().map(|p| p.owner_input()))
                    .collect::<Vec<_>>();
                self.inner
                    .$bulk(&inputs, num_threads)
                    .map_err(MorganReadError::Generator)
            }
        }
        impl Molecule {
            pub fn $method(
                &self,
                generator: &MorganFingerprintGenerator,
                params: Option<&MorganCallParams>,
                output: Option<&mut FingerprintAdditionalOutput>,
            ) -> Result<$result, MorganReadError> {
                let defaults = MorganCallParams::default();
                let params = params.unwrap_or(&defaults);
                let prepared = crate::morgan::prepare_morgan_read_input(self)
                    .map_err(MorganReadError::Preparation)?;
                generator
                    .inner
                    .$owner(&prepared.owner_input(), &params.owner_call(), output)
                    .map_err(MorganReadError::Generator)
            }
        }
    };
}
calculations!(
    fingerprints,
    morgan_fingerprint_with_generator,
    bits,
    Fingerprint
);
calculations!(
    counts,
    morgan_count_fingerprint_with_generator,
    count,
    SparseCountFingerprint32
);
calculations!(
    sparse_fingerprints,
    morgan_sparse_fingerprint_with_generator,
    sparse_bits,
    SparseBitFingerprint
);
calculations!(
    sparse_counts,
    morgan_sparse_count_fingerprint_with_generator,
    sparse_count,
    SparseCountFingerprint
);
