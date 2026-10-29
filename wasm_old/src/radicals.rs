//! Thin radical-assignment projection of the canonical public API.
use crate::Molecule;
impl Molecule {
    pub fn with_assigned_radicals(&self) -> Result<Self, cosmolkit::OperationError> {
        // COSMolKit❗✔️: self.inner.with_assigned_radicals()
        self.inner
            .borrow()
            .with_assigned_radicals()
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn assign_radicals_(&self) -> Result<(), cosmolkit::OperationError> {
        // COSMolKit❗✔️: self.inner.assign_radicals_()
        self.inner.borrow_mut().assign_radicals_()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit as ck;
    #[test]
    fn radical_projection_matches_owner_on_unsanitized_carbon_hydrogen_metal_and_charged_rows() {
        for text in [
            "", "C", "[C]", "[CH]", "[CH2]", "[CH3]", "[O]", "[OH]", "[H]", "[Fe]", "[Fe+3]",
            "[Na+]", "[Cl-]", "*", "C[CH2]", "[O][O]",
        ] {
            let owner = ck::Molecule::from_smiles_with_params(
                text,
                &ck::SmilesParseParams {
                    sanitize: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let m = Molecule::from_smiles_with_sanitize(text, false).unwrap();
            let before = m.inner.borrow().clone();
            let expected = owner.with_assigned_radicals();
            let actual = m.with_assigned_radicals();
            match (actual, expected) {
                (Ok(a), Ok(b)) => {
                    assert_eq!(*a.inner.borrow(), b);
                    if text == "[C]" {
                        assert_eq!(b.atom(ck::AtomId::new(0)).unwrap().radical_electrons(), 4);
                    }
                }
                (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
                _ => panic!("changed radical outcome"),
            }
            assert_eq!(*m.inner.borrow(), before);
            let mut expected = owner;
            let result = expected.assign_radicals_();
            assert_eq!(
                m.assign_radicals_().map_err(|e| e.to_string()),
                result.map_err(|e| e.to_string())
            );
            assert_eq!(*m.inner.borrow(), expected);
        }
    }
}
