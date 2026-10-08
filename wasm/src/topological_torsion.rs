//! Canonical scalar and persistent TopologicalTorsion transport, without molecule copies.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn topological_torsion_fingerprint_with_generator(
        &self,
        generator: &ck::TopologicalTorsionFingerprintGenerator,
        params: Option<&ck::TopologicalTorsionCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::Fingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_fingerprint_with_generator(
        self.inner
            .borrow()
            .topological_torsion_fingerprint_with_generator(generator, params, output)
    }
    pub fn topological_torsion_count_fingerprint_with_generator(
        &self,
        generator: &ck::TopologicalTorsionFingerprintGenerator,
        params: Option<&ck::TopologicalTorsionCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint32, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_count_fingerprint_with_generator(
        self.inner
            .borrow()
            .topological_torsion_count_fingerprint_with_generator(generator, params, output)
    }
    pub fn topological_torsion_sparse_fingerprint_with_generator(
        &self,
        generator: &ck::TopologicalTorsionFingerprintGenerator,
        params: Option<&ck::TopologicalTorsionCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseBitFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_sparse_fingerprint_with_generator(
        self.inner
            .borrow()
            .topological_torsion_sparse_fingerprint_with_generator(generator, params, output)
    }
    pub fn topological_torsion_sparse_count_fingerprint_with_generator(
        &self,
        generator: &ck::TopologicalTorsionFingerprintGenerator,
        params: Option<&ck::TopologicalTorsionCallParams>,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_sparse_count_fingerprint_with_generator(
        self.inner
            .borrow()
            .topological_torsion_sparse_count_fingerprint_with_generator(generator, params, output)
    }
    pub fn topological_torsion_sparse_count_fingerprint(
        &self,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_sparse_count_fingerprint(
        self.inner
            .borrow()
            .topological_torsion_sparse_count_fingerprint()
    }
    pub fn topological_torsion_sparse_count_fingerprint_with_params(
        &self,
        params: &ck::TopologicalTorsionFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_sparse_count_fingerprint_with_params(
        self.inner
            .borrow()
            .topological_torsion_sparse_count_fingerprint_with_params(params, output)
    }
    pub fn topological_torsion_sparse_fingerprint(
        &self,
    ) -> Result<ck::SparseBitFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_sparse_fingerprint(
        self.inner.borrow().topological_torsion_sparse_fingerprint()
    }
    pub fn topological_torsion_sparse_fingerprint_with_params(
        &self,
        params: &ck::TopologicalTorsionFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseBitFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_sparse_fingerprint_with_params(
        self.inner
            .borrow()
            .topological_torsion_sparse_fingerprint_with_params(params, output)
    }
    pub fn topological_torsion_count_fingerprint(
        &self,
    ) -> Result<ck::SparseCountFingerprint32, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_count_fingerprint(
        self.inner.borrow().topological_torsion_count_fingerprint()
    }
    pub fn topological_torsion_count_fingerprint_with_params(
        &self,
        params: &ck::TopologicalTorsionFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint32, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_count_fingerprint_with_params(
        self.inner
            .borrow()
            .topological_torsion_count_fingerprint_with_params(params, output)
    }
    pub fn topological_torsion_fingerprint(
        &self,
    ) -> Result<ck::Fingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_fingerprint(
        self.inner.borrow().topological_torsion_fingerprint()
    }
    pub fn topological_torsion_fingerprint_with_params(
        &self,
        params: &ck::TopologicalTorsionFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::Fingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_fingerprint_with_params(
        self.inner
            .borrow()
            .topological_torsion_fingerprint_with_params(params, output)
    }
}
pub fn topological_torsion_generator_fingerprints(
    generator: &ck::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::Fingerprint>>, ck::TopologicalTorsionReadError> {
    // COSMolKit❗✔️: self.inner.fingerprints(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.fingerprints(&rows, num_threads)
}
pub fn topological_torsion_generator_counts(
    generator: &ck::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::SparseCountFingerprint32>>, ck::TopologicalTorsionReadError> {
    // COSMolKit❗✔️: self.inner.counts(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.counts(&rows, num_threads)
}
pub fn topological_torsion_generator_sparse_fingerprints(
    generator: &ck::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::SparseBitFingerprint>>, ck::TopologicalTorsionReadError> {
    // COSMolKit❗✔️: self.inner.sparse_fingerprints(&rows, num_threads)
    let guards = molecules
        .iter()
        .map(|m| m.map(|m| m.inner.borrow()))
        .collect::<Vec<_>>();
    let rows = guards.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
    generator.sparse_fingerprints(&rows, num_threads)
}
pub fn topological_torsion_generator_sparse_counts(
    generator: &ck::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
    num_threads: i32,
) -> Result<Vec<Option<ck::SparseCountFingerprint>>, ck::TopologicalTorsionReadError> {
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
    fn topological_torsion_scalar_generator_collections_settings_and_errors() {
        let generator = ck::TopologicalTorsionFingerprintGenerator::new(None, None).unwrap();
        for text in ["CCO", "c1ccncc1", "F[C@](Cl)(Br)I"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            let params = ck::TopologicalTorsionFingerprintParams::default();
            let call = ck::TopologicalTorsionCallParams::default();
            let before = m.to_smiles().unwrap();
            assert_eq!(
                m.topological_torsion_fingerprint_with_generator(&generator, Some(&call), None)
                    .unwrap(),
                core.topological_torsion_fingerprint_with_generator(&generator, Some(&call), None)
                    .unwrap()
            );
            assert_eq!(
                m.topological_torsion_count_fingerprint_with_generator(
                    &generator,
                    Some(&call),
                    None
                )
                .unwrap(),
                core.topological_torsion_count_fingerprint_with_generator(
                    &generator,
                    Some(&call),
                    None
                )
                .unwrap()
            );
            assert_eq!(
                m.topological_torsion_sparse_fingerprint_with_generator(
                    &generator,
                    Some(&call),
                    None
                )
                .unwrap(),
                core.topological_torsion_sparse_fingerprint_with_generator(
                    &generator,
                    Some(&call),
                    None
                )
                .unwrap()
            );
            assert_eq!(
                m.topological_torsion_sparse_count_fingerprint_with_generator(
                    &generator,
                    Some(&call),
                    None
                )
                .unwrap(),
                core.topological_torsion_sparse_count_fingerprint_with_generator(
                    &generator,
                    Some(&call),
                    None
                )
                .unwrap()
            );
            assert_eq!(
                m.topological_torsion_sparse_count_fingerprint().unwrap(),
                core.topological_torsion_sparse_count_fingerprint().unwrap()
            );
            assert_eq!(
                m.topological_torsion_sparse_count_fingerprint_with_params(&params, None)
                    .unwrap(),
                core.topological_torsion_sparse_count_fingerprint_with_params(&params, None)
                    .unwrap()
            );
            assert_eq!(
                m.topological_torsion_sparse_fingerprint().unwrap(),
                core.topological_torsion_sparse_fingerprint().unwrap()
            );
            assert_eq!(
                m.topological_torsion_sparse_fingerprint_with_params(&params, None)
                    .unwrap(),
                core.topological_torsion_sparse_fingerprint_with_params(&params, None)
                    .unwrap()
            );
            assert_eq!(
                m.topological_torsion_count_fingerprint().unwrap(),
                core.topological_torsion_count_fingerprint().unwrap()
            );
            assert_eq!(
                m.topological_torsion_count_fingerprint_with_params(&params, None)
                    .unwrap(),
                core.topological_torsion_count_fingerprint_with_params(&params, None)
                    .unwrap()
            );
            assert_eq!(
                m.topological_torsion_fingerprint().unwrap(),
                core.topological_torsion_fingerprint().unwrap()
            );
            assert_eq!(
                m.topological_torsion_fingerprint_with_params(&params, None)
                    .unwrap(),
                core.topological_torsion_fingerprint_with_params(&params, None)
                    .unwrap()
            );
            let mut a = ck::FingerprintAdditionalOutput::new();
            let mut b = ck::FingerprintAdditionalOutput::new();
            a.allocate_atom_counts();
            b.allocate_atom_counts();
            assert_eq!(
                m.topological_torsion_fingerprint_with_generator(&generator, None, Some(&mut a))
                    .unwrap(),
                core.topological_torsion_fingerprint_with_generator(&generator, None, Some(&mut b))
                    .unwrap()
            );
            assert_eq!(a, b);
            assert_eq!(m.to_smiles().unwrap(), before);
            for threads in [1, 2] {
                assert_eq!(
                    topological_torsion_generator_fingerprints(
                        &generator,
                        &[Some(&m), None, Some(&m)],
                        threads
                    )
                    .unwrap(),
                    generator
                        .fingerprints(&[Some(&core), None, Some(&core)], threads)
                        .unwrap()
                );
            }
            for threads in [1, 2] {
                assert_eq!(
                    topological_torsion_generator_counts(
                        &generator,
                        &[Some(&m), None, Some(&m)],
                        threads
                    )
                    .unwrap(),
                    generator
                        .counts(&[Some(&core), None, Some(&core)], threads)
                        .unwrap()
                );
            }
            for threads in [1, 2] {
                assert_eq!(
                    topological_torsion_generator_sparse_fingerprints(
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
                    topological_torsion_generator_sparse_counts(
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
        settings.set_torsion_atom_count(1).unwrap();
        assert_eq!(retained.torsion_atom_count().unwrap(), 1);
        settings.set_fp_size(128).unwrap();
        assert_eq!(generator.settings().fp_size().unwrap(), 128);
        let snapshot = settings.params().unwrap();
        settings.set_torsion_atom_count(3).unwrap();
        assert_eq!(snapshot.torsion_atom_count, 1);
        drop(generator);
        assert_eq!(settings.torsion_atom_count().unwrap(), 3);
        let raw = Molecule::from_smiles_with_params(
            "CCO",
            &ck::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(matches!(
            raw.topological_torsion_fingerprint(),
            Err(ck::TopologicalTorsionReadError::Preparation(_))
        ));
    }
}

impl Molecule {
    pub fn legacy_topological_torsion_sparse_count_fingerprint(
        &self,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .legacy_topological_torsion_sparse_count_fingerprint(
        self.inner
            .borrow()
            .legacy_topological_torsion_sparse_count_fingerprint()
    }
    pub fn legacy_topological_torsion_sparse_count_fingerprint_with_params(
        &self,
        params: &ck::LegacyTopologicalTorsionParams,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .legacy_topological_torsion_sparse_count_fingerprint_with_params(
        self.inner
            .borrow()
            .legacy_topological_torsion_sparse_count_fingerprint_with_params(params)
    }
    pub fn legacy_topological_torsion_count_fingerprint(
        &self,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .legacy_topological_torsion_count_fingerprint(
        self.inner
            .borrow()
            .legacy_topological_torsion_count_fingerprint()
    }
    pub fn legacy_topological_torsion_count_fingerprint_with_params(
        &self,
        params: &ck::LegacyTopologicalTorsionParams,
    ) -> Result<ck::SparseCountFingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .legacy_topological_torsion_count_fingerprint_with_params(
        self.inner
            .borrow()
            .legacy_topological_torsion_count_fingerprint_with_params(params)
    }
    pub fn legacy_topological_torsion_fingerprint(
        &self,
    ) -> Result<ck::Fingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .legacy_topological_torsion_fingerprint(
        self.inner.borrow().legacy_topological_torsion_fingerprint()
    }
    pub fn legacy_topological_torsion_fingerprint_with_params(
        &self,
        params: &ck::LegacyTopologicalTorsionParams,
    ) -> Result<ck::Fingerprint, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .legacy_topological_torsion_fingerprint_with_params(
        self.inner
            .borrow()
            .legacy_topological_torsion_fingerprint_with_params(params)
    }
    pub fn topological_torsion_ids(&self) -> Result<Vec<u64>, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_ids(
        self.inner.borrow().topological_torsion_ids()
    }
    pub fn topological_torsion_ids_with_params(
        &self,
        torsion_atom_count: u32,
    ) -> Result<Vec<u64>, ck::TopologicalTorsionReadError> {
        // COSMolKit❗✔️: .topological_torsion_ids_with_params(
        self.inner
            .borrow()
            .topological_torsion_ids_with_params(torsion_atom_count)
    }
}

#[cfg(test)]
mod legacy_tests {
    use super::*;
    #[test]
    fn torsion_legacy_ids_and_parameter_variants() {
        for text in ["CCCCO", "c1ccncc1", "CC[C@H](F)CO"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            assert_eq!(
                m.legacy_topological_torsion_fingerprint().unwrap(),
                core.legacy_topological_torsion_fingerprint().unwrap()
            );
            for count in [3, 4, 5] {
                let p = ck::LegacyTopologicalTorsionParams {
                    torsion_atom_count: count,
                    ..Default::default()
                };
                assert_eq!(
                    m.legacy_topological_torsion_fingerprint_with_params(&p)
                        .unwrap(),
                    core.legacy_topological_torsion_fingerprint_with_params(&p)
                        .unwrap()
                );
            }
            assert_eq!(
                m.legacy_topological_torsion_count_fingerprint().unwrap(),
                core.legacy_topological_torsion_count_fingerprint().unwrap()
            );
            for count in [3, 4, 5] {
                let p = ck::LegacyTopologicalTorsionParams {
                    torsion_atom_count: count,
                    ..Default::default()
                };
                assert_eq!(
                    m.legacy_topological_torsion_count_fingerprint_with_params(&p)
                        .unwrap(),
                    core.legacy_topological_torsion_count_fingerprint_with_params(&p)
                        .unwrap()
                );
            }
            assert_eq!(
                m.legacy_topological_torsion_sparse_count_fingerprint()
                    .unwrap(),
                core.legacy_topological_torsion_sparse_count_fingerprint()
                    .unwrap()
            );
            for count in [3, 4, 5] {
                let p = ck::LegacyTopologicalTorsionParams {
                    torsion_atom_count: count,
                    ..Default::default()
                };
                assert_eq!(
                    m.legacy_topological_torsion_sparse_count_fingerprint_with_params(&p)
                        .unwrap(),
                    core.legacy_topological_torsion_sparse_count_fingerprint_with_params(&p)
                        .unwrap()
                );
            }
            assert_eq!(
                m.topological_torsion_ids().unwrap(),
                core.topological_torsion_ids().unwrap()
            );
            for count in [3, 4, 5] {
                assert_eq!(
                    m.topological_torsion_ids_with_params(count).unwrap(),
                    core.topological_torsion_ids_with_params(count).unwrap()
                );
            }
            for count in [3, 5] {
                let p = ck::TopologicalTorsionFingerprintParams {
                    generator: ck::TopologicalTorsionParams {
                        torsion_atom_count: count,
                        only_shortest_paths: true,
                        ..Default::default()
                    },
                    ..Default::default()
                };
                assert_eq!(
                    m.topological_torsion_fingerprint_with_params(&p, None)
                        .unwrap(),
                    core.topological_torsion_fingerprint_with_params(&p, None)
                        .unwrap()
                );
                assert_eq!(
                    m.topological_torsion_count_fingerprint_with_params(&p, None)
                        .unwrap(),
                    core.topological_torsion_count_fingerprint_with_params(&p, None)
                        .unwrap()
                );
                assert_eq!(
                    m.topological_torsion_sparse_fingerprint_with_params(&p, None)
                        .unwrap(),
                    core.topological_torsion_sparse_fingerprint_with_params(&p, None)
                        .unwrap()
                );
                assert_eq!(
                    m.topological_torsion_sparse_count_fingerprint_with_params(&p, None)
                        .unwrap(),
                    core.topological_torsion_sparse_count_fingerprint_with_params(&p, None)
                        .unwrap()
                );
            }
        }
    }
}
