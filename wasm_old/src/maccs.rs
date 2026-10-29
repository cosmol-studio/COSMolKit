//! Thin MACCS transport through the public facade.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn fingerprint_maccs(&self) -> Result<ck::Fingerprint, ck::MaccsFingerprintError> {
        // COSMolKit❗✔️: .maccs_fingerprint()
        self.inner.borrow().fingerprint_maccs()
    }
    pub fn fingerprint_maccs_raw(&self) -> Result<ck::Fingerprint, ck::MaccsFingerprintError> {
        // COSMolKit❗✔️: .maccs_fingerprint_raw()
        self.inner.borrow().fingerprint_maccs_raw()
    }
    pub fn fingerprint_maccs_with_params(
        &self,
        params: &ck::MaccsFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::MaccsFingerprintError> {
        // COSMolKit❗✔️: .maccs_fingerprint_with_params(&params.inner)
        self.inner.borrow().fingerprint_maccs_with_params(params)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn maccs_public_raw_parameters_and_errors() {
        for text in ["CCO", "c1ccncc1", "CC(=O)N", "CC.O"] {
            let m = Molecule::from_smiles(text).unwrap();
            let canonical = ck::Molecule::from_smiles(text).unwrap();
            let p = ck::MaccsFingerprintParams::default();
            assert_eq!(
                m.fingerprint_maccs().unwrap(),
                canonical.fingerprint_maccs().unwrap()
            );
            assert_eq!(
                m.fingerprint_maccs_raw().unwrap(),
                canonical.fingerprint_maccs_raw().unwrap()
            );
            assert_eq!(
                m.fingerprint_maccs_with_params(&p).unwrap(),
                canonical.fingerprint_maccs_with_params(&p).unwrap()
            );
            for width in [0, 64, 167] {
                let p = ck::MaccsFingerprintParams { n_bits: width };
                assert_eq!(
                    m.fingerprint_maccs_with_params(&p).unwrap_err(),
                    canonical.fingerprint_maccs_with_params(&p).unwrap_err()
                );
            }
        }
    }
}
