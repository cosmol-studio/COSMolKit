//! Canonical topological scalar transport.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn fingerprint_topological(
        &self,
    ) -> Result<ck::Fingerprint, ck::TopologicalFingerprintError> {
        // COSMolKit❗✔️: .topological_fingerprint(
        self.inner.borrow().fingerprint_topological()
    }
    pub fn fingerprint_topological_with_params(
        &self,
        params: &ck::TopologicalFingerprintParams,
    ) -> Result<ck::Fingerprint, ck::TopologicalFingerprintError> {
        // COSMolKit❗✔️: .topological_fingerprint_with_params(
        self.inner
            .borrow()
            .fingerprint_topological_with_params(params)
    }
    pub fn fingerprint_topological_with_output(
        &self,
    ) -> Result<ck::TopologicalFingerprintResult, ck::TopologicalFingerprintError> {
        // COSMolKit❗✔️: .topological_fingerprint_with_output(
        self.inner.borrow().fingerprint_topological_with_output()
    }
    pub fn fingerprint_topological_with_output_with_params(
        &self,
        params: &ck::TopologicalFingerprintParams,
        request: ck::TopologicalFingerprintOutputRequest,
    ) -> Result<ck::TopologicalFingerprintResult, ck::TopologicalFingerprintError> {
        // COSMolKit❗✔️: .topological_fingerprint_with_output_with_params(
        self.inner
            .borrow()
            .fingerprint_topological_with_output_with_params(params, request)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn topological_all_calls_and_output_requests_match_canonical() {
        for text in ["CCO", "c1ccncc1", "CC(=O)N"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            let params = ck::TopologicalFingerprintParams::default();
            assert_eq!(
                m.fingerprint_topological().unwrap(),
                core.fingerprint_topological().unwrap()
            );
            assert_eq!(
                m.fingerprint_topological_with_params(&params).unwrap(),
                core.fingerprint_topological_with_params(&params).unwrap()
            );
            assert_eq!(
                m.fingerprint_topological_with_output().unwrap(),
                core.fingerprint_topological_with_output().unwrap()
            );
            for atom_bits in [false, true] {
                for bit_info in [false, true] {
                    let request = ck::TopologicalFingerprintOutputRequest {
                        atom_bits,
                        bit_info,
                    };
                    let a = m
                        .fingerprint_topological_with_output_with_params(&params, request)
                        .unwrap();
                    let b = core
                        .fingerprint_topological_with_output_with_params(&params, request)
                        .unwrap();
                    assert_eq!(a, b);
                    assert_eq!(a.atom_bits().is_ok(), atom_bits);
                    assert_eq!(a.bit_info().is_ok(), bit_info);
                }
            }
        }
    }
}
