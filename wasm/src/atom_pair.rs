//! Scalar AtomPair projection through the public facade only.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn fingerprint_atom_pair(&self) -> Result<ck::Fingerprint, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_fingerprint(
        self.inner.borrow().fingerprint_atom_pair()
    }
    pub fn fingerprint_atom_pair_with_params(
        &self,
        params: &ck::AtomPairFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::Fingerprint, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_atom_pair_with_params(params, output)
    }
    pub fn fingerprint_atom_pair_sparse(
        &self,
    ) -> Result<ck::SparseBitFingerprint, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_sparse_fingerprint(
        self.inner.borrow().fingerprint_atom_pair_sparse()
    }
    pub fn fingerprint_atom_pair_sparse_with_params(
        &self,
        params: &ck::AtomPairFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseBitFingerprint, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_sparse_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_atom_pair_sparse_with_params(params, output)
    }
    pub fn fingerprint_atom_pair_count(
        &self,
    ) -> Result<ck::SparseCountFingerprint32, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_count_fingerprint(
        self.inner.borrow().fingerprint_atom_pair_count()
    }
    pub fn fingerprint_atom_pair_count_with_params(
        &self,
        params: &ck::AtomPairFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint32, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_count_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_atom_pair_count_with_params(params, output)
    }
    pub fn fingerprint_atom_pair_sparse_count(
        &self,
    ) -> Result<ck::SparseCountFingerprint, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_sparse_count_fingerprint(
        self.inner.borrow().fingerprint_atom_pair_sparse_count()
    }
    pub fn fingerprint_atom_pair_sparse_count_with_params(
        &self,
        params: &ck::AtomPairFingerprintParams,
        output: Option<&mut ck::FingerprintAdditionalOutput>,
    ) -> Result<ck::SparseCountFingerprint, ck::AtomPairReadError> {
        // COSMolKit❗✔️: .atom_pair_sparse_count_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_atom_pair_sparse_count_with_params(params, output)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn atom_pair_scalar_values_outputs_options_and_source_failures() {
        for text in ["CCO", "c1ccncc1", "F[C@](Cl)(Br)I"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            let before = m.to_smiles().unwrap();
            assert_eq!(
                m.fingerprint_atom_pair().unwrap(),
                core.fingerprint_atom_pair().unwrap()
            );
            assert_eq!(
                m.fingerprint_atom_pair_sparse().unwrap(),
                core.fingerprint_atom_pair_sparse().unwrap()
            );
            assert_eq!(
                m.fingerprint_atom_pair_count().unwrap(),
                core.fingerprint_atom_pair_count().unwrap()
            );
            assert_eq!(
                m.fingerprint_atom_pair_sparse_count().unwrap(),
                core.fingerprint_atom_pair_sparse_count().unwrap()
            );
            for params in [
                ck::AtomPairFingerprintParams::default(),
                ck::AtomPairFingerprintParams {
                    from_atoms: Some(vec![]),
                    ..Default::default()
                },
                ck::AtomPairFingerprintParams {
                    ignore_atoms: Some(vec![0]),
                    ..Default::default()
                },
                ck::AtomPairFingerprintParams {
                    custom_atom_invariants: Some(vec![33; usize::try_from(m.num_atoms()).unwrap()]),
                    ..Default::default()
                },
            ] {
                let mut a = ck::FingerprintAdditionalOutput::new();
                a.allocate_atom_counts();
                a.allocate_atom_to_bits();
                a.allocate_bit_info_map();
                a.allocate_bit_paths();
                a.allocate_atoms_per_bit();
                let mut b = ck::FingerprintAdditionalOutput::new();
                b.allocate_atom_counts();
                b.allocate_atom_to_bits();
                b.allocate_bit_info_map();
                b.allocate_bit_paths();
                b.allocate_atoms_per_bit();
                assert_eq!(
                    m.fingerprint_atom_pair_with_params(&params, Some(&mut a))
                        .unwrap(),
                    core.fingerprint_atom_pair_with_params(&params, Some(&mut b))
                        .unwrap()
                );
                assert_eq!(a, b);
                assert_eq!(
                    m.fingerprint_atom_pair_sparse_with_params(&params, Some(&mut a))
                        .unwrap(),
                    core.fingerprint_atom_pair_sparse_with_params(&params, Some(&mut b))
                        .unwrap()
                );
                assert_eq!(a, b);
                assert_eq!(
                    m.fingerprint_atom_pair_count_with_params(&params, Some(&mut a))
                        .unwrap(),
                    core.fingerprint_atom_pair_count_with_params(&params, Some(&mut b))
                        .unwrap()
                );
                assert_eq!(a, b);
                assert_eq!(
                    m.fingerprint_atom_pair_sparse_count_with_params(&params, Some(&mut a))
                        .unwrap(),
                    core.fingerprint_atom_pair_sparse_count_with_params(&params, Some(&mut b))
                        .unwrap()
                );
                assert_eq!(a, b);
                assert_eq!(m.to_smiles().unwrap(), before);
            }
        }
        let raw = Molecule::from_smiles_with_params(
            "CC",
            &ck::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(matches!(
            raw.fingerprint_atom_pair(),
            Err(ck::AtomPairReadError::Preparation(
                ck::FingerprintPreparationError::MissingPreparedValence
            ))
        ));
        let m = Molecule::from_smiles("CCO").unwrap();
        for p in [ck::AtomPairFingerprintParams {
            generator: ck::AtomPairParams {
                fp_size: 0,
                ..Default::default()
            },
            ..Default::default()
        }] {
            assert!(matches!(
                m.fingerprint_atom_pair_with_params(&p, None),
                Err(ck::AtomPairReadError::Generator(_))
            ));
        }
        assert!(
            m.fingerprint_atom_pair_with_params(
                &ck::AtomPairFingerprintParams {
                    from_atoms: Some(vec![99]),
                    ..Default::default()
                },
                None
            )
            .unwrap()
            .on_bits()
            .is_empty()
        );
        assert_eq!(m.to_smiles().unwrap(), "CCO".into());
    }
}
