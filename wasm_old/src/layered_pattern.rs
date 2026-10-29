//! Read-only concrete Layered/Pattern transport.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn fingerprint_layered(&self) -> Result<ck::Fingerprint, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint(
        self.inner.borrow().fingerprint_layered()
    }
    pub fn fingerprint_layered_with_params(
        &self,
        params: &ck::LayeredFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint_with_params(
        self.inner.borrow().fingerprint_layered_with_params(params)
    }
    pub fn fingerprint_layered_with_output(
        &self,
    ) -> Result<ck::LayeredFingerprintResult, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output(
        self.inner.borrow().fingerprint_layered_with_output()
    }
    pub fn fingerprint_layered_with_output_with_params(
        &self,
        params: &ck::LayeredFingerprintParams,
    ) -> Result<ck::LayeredFingerprintResult, ck::LayeredFingerprintError> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output_with_params(
        self.inner
            .borrow()
            .fingerprint_layered_with_output_with_params(params)
    }
    pub fn fingerprint_pattern(&self) -> Result<ck::Fingerprint, ck::PatternFingerprintError> {
        // COSMolKit❗✔️: .pattern_fingerprint(
        self.inner.borrow().fingerprint_pattern()
    }
    pub fn fingerprint_pattern_with_params(
        &self,
        params: &ck::PatternFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::PatternFingerprintError> {
        // COSMolKit❗✔️: .pattern_fingerprint_with_params(
        self.inner.borrow().fingerprint_pattern_with_params(params)
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
                m.fingerprint_layered().unwrap(),
                canonical.fingerprint_layered().unwrap()
            );
            assert_eq!(
                m.fingerprint_layered_with_params(&ck::LayeredFingerprintParams::default())
                    .unwrap(),
                canonical
                    .fingerprint_layered_with_params(&ck::LayeredFingerprintParams::default())
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_layered_with_output().unwrap(),
                canonical.fingerprint_layered_with_output().unwrap()
            );
            assert_eq!(
                m.fingerprint_layered_with_output_with_params(
                    &ck::LayeredFingerprintParams::default()
                )
                .unwrap(),
                canonical
                    .fingerprint_layered_with_output_with_params(
                        &ck::LayeredFingerprintParams::default()
                    )
                    .unwrap()
            );
            assert_eq!(
                m.fingerprint_pattern().unwrap(),
                canonical.fingerprint_pattern().unwrap()
            );
            assert_eq!(
                m.fingerprint_pattern_with_params(&ck::PatternFingerprintParams::default())
                    .unwrap(),
                canonical
                    .fingerprint_pattern_with_params(&ck::PatternFingerprintParams::default())
                    .unwrap()
            );
            let p = ck::LayeredFingerprintParams {
                atom_counts: Some(vec![0; canonical.num_atoms()]),
                ..Default::default()
            };
            assert_eq!(
                m.fingerprint_layered_with_output_with_params(&p).unwrap(),
                canonical
                    .fingerprint_layered_with_output_with_params(&p)
                    .unwrap()
            );
            assert!(matches!(
                m.fingerprint_layered_with_params(&ck::LayeredFingerprintParams {
                    min_path: 0,
                    ..Default::default()
                }),
                Err(ck::LayeredFingerprintError::InvalidArguments { .. })
            ));
            assert!(matches!(
                m.fingerprint_pattern_with_params(&ck::PatternFingerprintParams {
                    n_bits: 0,
                    ..Default::default()
                }),
                Err(ck::PatternFingerprintError::EmptyFingerprint)
            ));
        }
    }
}
