//! Canonical sanitization and read-only chemistry problem detection.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn sanitize(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.sanitize()
        self.inner.borrow().sanitize().map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn sanitize_with_params(
        &self,
        params: &ck::SanitizeParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.sanitize_with_params(params)
        self.inner
            .borrow()
            .sanitize_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn sanitize_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.sanitize_()
        self.inner.borrow_mut().sanitize_()
    }
    pub fn sanitize_with_params_(
        &self,
        params: &ck::SanitizeParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.sanitize_with_params_(params)
        self.inner.borrow_mut().sanitize_with_params_(params)
    }
    pub fn detect_chemistry_problems(
        &self,
    ) -> Result<ck::ChemistryProblemReport, ck::SanitizeError> {
        // COSMolKit❗✔️: self.inner.detect_chemistry_problems()
        self.inner.borrow().detect_chemistry_problems()
    }
    pub fn detect_chemistry_problems_with_params(
        &self,
        params: &ck::SanitizeParams,
    ) -> Result<ck::ChemistryProblemReport, ck::SanitizeError> {
        // COSMolKit❗✔️: self.inner.detect_chemistry_problems_with_params(params)
        self.inner
            .borrow()
            .detect_chemistry_problems_with_params(params)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn compare(
        a: Result<Molecule, ck::OperationError>,
        b: Result<ck::Molecule, ck::OperationError>,
    ) {
        match (a, b) {
            (Ok(a), Ok(b)) => assert_eq!(*a.inner.borrow(), b),
            (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
            _ => panic!("changed sanitize result"),
        }
    }
    #[test]
    fn sanitize_operations_and_ordered_reports_match_owner_for_all_stage_selections() {
        let selections = [
            ck::SanitizeOperations::NONE,
            ck::SanitizeOperations::ALL,
            ck::SanitizeOperations::CLEANUP,
            ck::SanitizeOperations::PROPERTIES,
            ck::SanitizeOperations::SYMM_RINGS,
            ck::SanitizeOperations::KEKULIZE,
            ck::SanitizeOperations::FIND_RADICALS,
            ck::SanitizeOperations::SET_AROMATICITY,
            ck::SanitizeOperations::SET_CONJUGATION,
            ck::SanitizeOperations::SET_HYBRIDIZATION,
            ck::SanitizeOperations::CLEANUP_CHIRALITY,
            ck::SanitizeOperations::ADJUST_HS,
            ck::SanitizeOperations::CLEANUP_ORGANOMETALLICS,
            ck::SanitizeOperations::CLEANUP_ATROPISOMERS,
            ck::SanitizeOperations::PROPERTIES | ck::SanitizeOperations::KEKULIZE,
        ];
        for text in [
            "",
            "CCO",
            "c1ccccc1",
            "[CH3]",
            "C(F)(F)(F)(F)F",
            "c",
            "c1cccc1",
            "C(F)(F)(F)(F)F.C(F)(F)(F)(F)F",
        ] {
            let parse = ck::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..Default::default()
            };
            let owner = ck::Molecule::from_smiles_with_params(text, &parse).unwrap();
            let m = Molecule::from_smiles_with_params(text, &parse).unwrap();
            let before = m.inner.borrow().clone();
            compare(m.sanitize(), owner.sanitize());
            assert_eq!(
                m.detect_chemistry_problems(),
                owner.detect_chemistry_problems()
            );
            let fresh = Molecule::from_smiles_with_params(text, &parse).unwrap();
            let mut expected = owner.clone();
            assert_eq!(
                fresh.sanitize_().map_err(|e| e.to_string()),
                expected.sanitize_().map_err(|e| e.to_string())
            );
            assert_eq!(*fresh.inner.borrow(), expected);
            for operations in selections {
                let params = ck::SanitizeParams { operations };
                compare(
                    m.sanitize_with_params(&params),
                    owner.sanitize_with_params(&params),
                );
                assert_eq!(
                    m.detect_chemistry_problems_with_params(&params),
                    owner.detect_chemistry_problems_with_params(&params)
                );
                let fresh = Molecule::from_smiles_with_params(text, &parse).unwrap();
                let mut expected = owner.clone();
                let result = fresh.sanitize_with_params_(&params);
                let failed = result.is_err();
                assert_eq!(
                    result.map_err(|e| e.to_string()),
                    expected
                        .sanitize_with_params_(&params)
                        .map_err(|e| e.to_string())
                );
                assert_eq!(*fresh.inner.borrow(), expected);
                if failed {
                    assert_eq!(*fresh.inner.borrow(), before);
                }
            }
            assert_eq!(*m.inner.borrow(), before);
        }
    }
}
