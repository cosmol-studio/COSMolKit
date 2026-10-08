//! Read-only concrete Layered/Pattern transport.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn layered_fingerprint(&self) -> Result<ck::Fingerprint, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint(
        self.inner.borrow().layered_fingerprint()
    }
    pub fn layered_fingerprint_with_params(
        &self,
        params: &ck::LayeredFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint_with_params(
        self.inner.borrow().layered_fingerprint_with_params(params)
    }
    pub fn layered_fingerprint_with_output(
        &self,
    ) -> Result<ck::LayeredFingerprintResult, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output(
        self.inner.borrow().layered_fingerprint_with_output()
    }
    pub fn layered_fingerprint_with_output_with_params(
        &self,
        params: &ck::LayeredFingerprintParams,
    ) -> Result<ck::LayeredFingerprintResult, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output_with_params(
        self.inner
            .borrow()
            .layered_fingerprint_with_output_with_params(params)
    }
    pub fn pattern_fingerprint(&self) -> Result<ck::Fingerprint, ck::PatternFingerprintError> {
        // COSMolKit❗✔️: .pattern_fingerprint(
        self.inner.borrow().pattern_fingerprint()
    }
    pub fn pattern_fingerprint_with_params(
        &self,
        params: &ck::PatternFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::PatternFingerprintError> {
        // COSMolKit❗✔️: .pattern_fingerprint_with_params(
        self.inner.borrow().pattern_fingerprint_with_params(params)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn layered_pattern_concrete_defaults_parameters_outputs_and_errors() {
        for text in ["CCO", "c1ccncc1", "CC(=O)N"] {
            let m = Molecule::from_smiles(text).unwrap();
            let canonical = ck::Molecule::from_smiles(text).unwrap();
            assert_eq!(
                m.layered_fingerprint().unwrap(),
                canonical.layered_fingerprint().unwrap()
            );
            assert_eq!(
                m.layered_fingerprint_with_params(&ck::LayeredFingerprintParams::default())
                    .unwrap(),
                canonical
                    .layered_fingerprint_with_params(&ck::LayeredFingerprintParams::default())
                    .unwrap()
            );
            assert_eq!(
                m.layered_fingerprint_with_output().unwrap(),
                canonical.layered_fingerprint_with_output().unwrap()
            );
            assert_eq!(
                m.layered_fingerprint_with_output_with_params(
                    &ck::LayeredFingerprintParams::default()
                )
                .unwrap(),
                canonical
                    .layered_fingerprint_with_output_with_params(
                        &ck::LayeredFingerprintParams::default()
                    )
                    .unwrap()
            );
            assert_eq!(
                m.pattern_fingerprint().unwrap(),
                canonical.pattern_fingerprint().unwrap()
            );
            assert_eq!(
                m.pattern_fingerprint_with_params(&ck::PatternFingerprintParams::default())
                    .unwrap(),
                canonical
                    .pattern_fingerprint_with_params(&ck::PatternFingerprintParams::default())
                    .unwrap()
            );
            let p = ck::LayeredFingerprintParams {
                atom_counts: Some(vec![0; canonical.num_atoms()]),
                ..Default::default()
            };
            assert_eq!(
                m.layered_fingerprint_with_output_with_params(&p).unwrap(),
                canonical
                    .layered_fingerprint_with_output_with_params(&p)
                    .unwrap()
            );
            assert!(matches!(
                m.layered_fingerprint_with_params(&ck::LayeredFingerprintParams {
                    min_path: 0,
                    ..Default::default()
                }),
                Err(ck::LayeredFingerprintError::InvalidArguments { .. })
            ));
            assert!(matches!(
                m.pattern_fingerprint_with_params(&ck::PatternFingerprintParams {
                    n_bits: 0,
                    ..Default::default()
                }),
                Err(ck::PatternFingerprintError::EmptyFingerprint)
            ));
        }
    }
}
