//! Canonical scalar and persistent Morgan transport, without molecule copies.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn fingerprint_morgan_with_generator(
        &self,
        generator: &ck::MorganFingerprintGenerator,
        params: Option<&ck::MorganCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::Fingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_fingerprint_with_generator(
        self.inner
            .borrow()
            .fingerprint_morgan_with_generator(generator, params, output)
    }
    pub fn fingerprint_morgan_count_with_generator(
        &self,
        generator: &ck::MorganFingerprintGenerator,
        params: Option<&ck::MorganCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint32, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_count_fingerprint_with_generator(
        self.inner
            .borrow()
            .fingerprint_morgan_count_with_generator(generator, params, output)
    }
    pub fn fingerprint_morgan_sparse_with_generator(
        &self,
        generator: &ck::MorganFingerprintGenerator,
        params: Option<&ck::MorganCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseBitFingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_sparse_fingerprint_with_generator(
        self.inner
            .borrow()
            .fingerprint_morgan_sparse_with_generator(generator, params, output)
    }
    pub fn fingerprint_morgan_sparse_count_with_generator(
        &self,
        generator: &ck::MorganFingerprintGenerator,
        params: Option<&ck::MorganCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_sparse_count_fingerprint_with_generator(
        self.inner
            .borrow()
            .fingerprint_morgan_sparse_count_with_generator(generator, params, output)
    }
    pub fn fingerprint_morgan_sparse_count(
        &self,
    ) -> Result<ck::SparseCountFingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_sparse_count_fingerprint(
        self.inner.borrow().fingerprint_morgan_sparse_count()
    }
    pub fn fingerprint_morgan_sparse_count_with_params(
        &self,
        params: &ck::MorganFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_sparse_count_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_morgan_sparse_count_with_params(params, output)
    }
    pub fn fingerprint_morgan_sparse(
        &self,
    ) -> Result<ck::SparseBitFingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_sparse_fingerprint(
        self.inner.borrow().fingerprint_morgan_sparse()
    }
    pub fn fingerprint_morgan_sparse_with_params(
        &self,
        params: &ck::MorganFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseBitFingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_sparse_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_morgan_sparse_with_params(params, output)
    }
    pub fn fingerprint_morgan_count(
        &self,
    ) -> Result<ck::SparseCountFingerprint32, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_count_fingerprint(
        self.inner.borrow().fingerprint_morgan_count()
    }
    pub fn fingerprint_morgan_count_with_params(
        &self,
        params: &ck::MorganFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint32, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_count_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_morgan_count_with_params(params, output)
    }
    pub fn fingerprint_morgan(&self) -> Result<ck::Fingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_fingerprint(
        self.inner.borrow().fingerprint_morgan()
    }
    pub fn fingerprint_morgan_with_params(
        &self,
        params: &ck::MorganFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::Fingerprint, ck::MorganReadError> {
        // COSMolKit❗✔️: .morgan_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_morgan_with_params(params, output)
    }
}
pub fn morgan_generator_fingerprints(
    generator: &ck::MorganFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::Fingerprint>>, ck::MorganReadError> {
    // COSMolKit❗✔️: self.inner.fingerprints(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.fingerprints(&rows, num_threads)
}
pub fn morgan_generator_counts(
    generator: &ck::MorganFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::SparseCountFingerprint32>>, ck::MorganReadError> {
    // COSMolKit❗✔️: self.inner.counts(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.counts(&rows, num_threads)
}
pub fn morgan_generator_sparse_fingerprints(
    generator: &ck::MorganFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::SparseBitFingerprint>>, ck::MorganReadError> {
    // COSMolKit❗✔️: self.inner.sparse_fingerprints(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.sparse_fingerprints(&rows, num_threads)
}
pub fn morgan_generator_sparse_counts(
    generator: &ck::MorganFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::SparseCountFingerprint>>, ck::MorganReadError> {
    // COSMolKit❗✔️: self.inner.sparse_counts(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.sparse_counts(&rows, num_threads)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn morgan_scalar_generator_collections_settings_and_errors() {
        let generator = ck::MorganFingerprintGenerator::new(None, None, None).unwrap();
        for text in ["CCO", "c1ccncc1", "F[C@](Cl)(Br)I"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            let params = ck::MorganFingerprintParams::default();
            let call = ck::MorganCallParams::default();
            let before = m.to_smiles().unwrap();
            assert_eq!(
                m.fingerprint_morgan_with_generator(&generator, Some(&call), None)
                    .unwrap(),
                core.fingerprint_morgan_with_generator(&generator, Some(&call), None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_count_with_generator(&generator, Some(&call), None)
                    .unwrap(),
                core.fingerprint_morgan_count_with_generator(&generator, Some(&call), None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_sparse_with_generator(&generator, Some(&call), None)
                    .unwrap(),
                core.fingerprint_morgan_sparse_with_generator(&generator, Some(&call), None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_sparse_count_with_generator(&generator, Some(&call), None)
                    .unwrap(),
                core.fingerprint_morgan_sparse_count_with_generator(&generator, Some(&call), None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_sparse_count().unwrap(),
                core.fingerprint_morgan_sparse_count().unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_sparse_count_with_params(&params, None)
                    .unwrap(),
                core.fingerprint_morgan_sparse_count_with_params(&params, None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_sparse().unwrap(),
                core.fingerprint_morgan_sparse().unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_sparse_with_params(&params, None)
                    .unwrap(),
                core.fingerprint_morgan_sparse_with_params(&params, None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_count().unwrap(),
                core.fingerprint_morgan_count().unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_count_with_params(&params, None)
                    .unwrap(),
                core.fingerprint_morgan_count_with_params(&params, None)
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan().unwrap(),
                core.fingerprint_morgan().unwrap()
            );
            assert_eq!(
                m.fingerprint_morgan_with_params(&params, None).unwrap(),
                core.fingerprint_morgan_with_params(&params, None).unwrap()
            );
            let mut a = ck::FingerprintAdditionalOutput::new();
            let mut b = ck::FingerprintAdditionalOutput::new();
            a.allocate_atom_counts();
            b.allocate_atom_counts();
            assert_eq!(
                m.fingerprint_morgan_with_generator(&generator, None, Some(&mut a))
                    .unwrap(),
                core.fingerprint_morgan_with_generator(&generator, None, Some(&mut b))
                    .unwrap()
            );
            assert_eq!(a, b);
            assert_eq!(m.to_smiles().unwrap(), before);
            for threads in [1, 2] {
                assert_eq!(
                    morgan_generator_fingerprints(&generator, &[Some(&m), None, Some(&m)], threads)
                        .unwrap(),
                    generator
                        .fingerprints(&[Some(&core), None, Some(&core)], threads)
                        .unwrap()
                );
            }
            for threads in [1, 2] {
                assert_eq!(
                    morgan_generator_counts(&generator, &[Some(&m), None, Some(&m)], threads)
                        .unwrap(),
                    generator
                        .counts(&[Some(&core), None, Some(&core)], threads)
                        .unwrap()
                );
            }
            for threads in [1, 2] {
                assert_eq!(
                    morgan_generator_sparse_fingerprints(
                        &generator,
                        &[Some(&m), None, Some(&m)],
                        threads
                    )
                    .unwrap(),
                    generator
                        .sparse_fingerprints(&[Some(&core), None, Some(&core)], threads)
                        .unwrap()
                );
            }
            for threads in [1, 2] {
                assert_eq!(
                    morgan_generator_sparse_counts(
                        &generator,
                        &[Some(&m), None, Some(&m)],
                        threads
                    )
                    .unwrap(),
                    generator
                        .sparse_counts(&[Some(&core), None, Some(&core)], threads)
                        .unwrap()
                );
            }
        }
        let mut settings = generator.settings();
        let retained = generator.settings();
        settings.set_radius(1).unwrap();
        assert_eq!(retained.radius().unwrap(), 1);
        settings.set_fp_size(128).unwrap();
        assert_eq!(generator.settings().fp_size().unwrap(), 128);
        let snapshot = settings.params().unwrap();
        settings.set_radius(3).unwrap();
        assert_eq!(snapshot.radius, 1);
        drop(generator);
        assert_eq!(settings.radius().unwrap(), 3);
        let raw = Molecule::from_smiles_with_params(
            "CCO",
            &ck::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(matches!(
            raw.fingerprint_morgan(),
            Err(ck::MorganReadError::Preparation(_))
        ));
    }
}
