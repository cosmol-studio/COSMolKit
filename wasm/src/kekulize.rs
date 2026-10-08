//! Thin typed Kekule operation projections.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn with_kekulized_bonds(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_kekulized_bonds()
        self.inner
            .borrow()
            .with_kekulized_bonds()
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn with_kekulized_bonds_with_params(
        &self,
        params: &ck::KekulizeParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_kekulized_bonds_with_params(params)
        self.inner
            .borrow()
            .with_kekulized_bonds_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn kekulize_bonds_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.kekulize_bonds_()
        self.inner.borrow_mut().kekulize_bonds_()
    }
    pub fn kekulize_bonds_with_params_(
        &self,
        params: &ck::KekulizeParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.kekulize_bonds_with_params_(params)
        self.inner.borrow_mut().kekulize_bonds_with_params_(params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn every_parameter_and_four_operations_match_canonical_and_preserve_failure_atomicity() {
        for text in ["c1ccccc1", "c1ccncc1", "c1ccc2ccccc2c1", "CCO"] {
            let source = ck::Molecule::from_smiles(text).unwrap();
            let m = Molecule::from_smiles(text).unwrap();
            let before = m.inner.borrow().clone();
            assert_eq!(
                *m.with_kekulized_bonds().unwrap().inner.borrow(),
                source.with_kekulized_bonds().unwrap()
            );
            let mut owner = source.clone();
            owner.kekulize_bonds_().unwrap();
            let mutating = Molecule::from_smiles(text).unwrap();
            mutating.kekulize_bonds_().unwrap();
            assert_eq!(*mutating.inner.borrow(), owner);
            for mark_atoms_bonds in [false, true] {
                for canonical in [false, true] {
                    for max_backtracks in [0, 1, 100, u32::MAX] {
                        let p = ck::KekulizeParams {
                            mark_atoms_bonds,
                            canonical,
                            max_backtracks,
                        };
                        let expected = source.with_kekulized_bonds_with_params(&p);
                        let result = m.with_kekulized_bonds_with_params(&p);
                        match (result, expected) {
                            (Ok(a), Ok(b)) => assert_eq!(*a.inner.borrow(), b),
                            (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
                            _ => panic!("projection changed outcome"),
                        };
                        let mut owner = source.clone();
                        let expected = owner.kekulize_bonds_with_params_(&p);
                        let mutating = Molecule::from_smiles(text).unwrap();
                        let result = mutating.kekulize_bonds_with_params_(&p);
                        assert_eq!(
                            result.map_err(|e| e.to_string()),
                            expected.map_err(|e| e.to_string())
                        );
                        assert_eq!(*mutating.inner.borrow(), owner);
                    }
                }
            }
            assert_eq!(*m.inner.borrow(), before);
        }
        for (text, kind) in [
            ("c", "AromaticAtomOutsideRing"),
            ("c1cccc1", "NotKekulizable"),
        ] {
            let m = Molecule::from_smiles_with_sanitize(text, false).unwrap();
            let before = m.inner.borrow().clone();
            for default in [false, true] {
                let error = if default {
                    m.kekulize_bonds_()
                } else {
                    m.kekulize_bonds_with_params_(&Default::default())
                }
                .unwrap_err();
                match error {
                    ck::OperationError::Kekulize(e) => {
                        assert!(format!("{e:?}").starts_with(kind));
                        if let ck::KekulizeError::NotKekulizable { problem_atoms } = e {
                            assert_eq!(
                                problem_atoms
                                    .iter()
                                    .map(|id| id.index())
                                    .collect::<Vec<_>>(),
                                [0, 1, 2, 3, 4]
                            );
                        }
                    }
                    other => panic!("wrong error {other:?}"),
                };
                assert_eq!(*m.inner.borrow(), before);
                assert!(m.with_kekulized_bonds().is_err());
                assert!(
                    m.with_kekulized_bonds_with_params(&Default::default())
                        .is_err()
                );
            }
        }
    }
}
