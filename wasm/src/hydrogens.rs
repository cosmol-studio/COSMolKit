//! Hydrogen operations project only the canonical runtime boundary.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn with_hydrogens(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_hydrogens()
        self.inner.borrow().with_hydrogens().map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn with_hydrogens_with_params(
        &self,
        params: &ck::AddHsParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_hydrogens_with_params(params)
        self.inner
            .borrow()
            .with_hydrogens_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn without_hydrogens(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.without_hydrogens()
        self.inner.borrow().without_hydrogens().map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn without_hydrogens_with_params(
        &self,
        params: &ck::RemoveHsParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.without_hydrogens_with_params(params)
        self.inner
            .borrow()
            .without_hydrogens_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn add_hydrogens_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.add_hydrogens_()
        self.inner.borrow_mut().add_hydrogens_()
    }
    pub fn add_hydrogens_with_params_(
        &self,
        params: &ck::AddHsParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.add_hydrogens_with_params_(params)
        self.inner.borrow_mut().add_hydrogens_with_params_(params)
    }
    pub fn remove_hydrogens_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.remove_hydrogens_()
        self.inner.borrow_mut().remove_hydrogens_()
    }
    pub fn remove_hydrogens_with_params_(
        &self,
        params: &ck::RemoveHsParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.remove_hydrogens_with_params_(params)
        self.inner
            .borrow_mut()
            .remove_hydrogens_with_params_(params)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hydrogens_all_eight_value_and_inplace_calls_match_canonical_and_fail_atomically() {
        for text in ["CCO", "c1ccccc1", "[NH4+]", "[2H]OC"] {
            let m = Molecule::from_smiles(text).unwrap();
            let mut core = ck::Molecule::from_smiles(text).unwrap();
            let original = core.clone();
            assert_eq!(
                *m.with_hydrogens().unwrap().inner.borrow(),
                core.with_hydrogens().unwrap()
            );
            for p in [
                ck::AddHsParams::default(),
                ck::AddHsParams {
                    explicit_only: true,
                    ..Default::default()
                },
                ck::AddHsParams {
                    only_on_atoms: Some(vec![ck::AtomId::new(0)]),
                    ..Default::default()
                },
            ] {
                assert_eq!(
                    *m.with_hydrogens_with_params(&p).unwrap().inner.borrow(),
                    core.with_hydrogens_with_params(&p).unwrap()
                );
            }
            assert_eq!(
                *m.without_hydrogens().unwrap().inner.borrow(),
                core.without_hydrogens().unwrap()
            );
            for p in [
                ck::RemoveHsParams::default(),
                ck::RemoveHsParams {
                    remove_isotopes: true,
                    sanitize: false,
                    ..Default::default()
                },
            ] {
                assert_eq!(
                    *m.without_hydrogens_with_params(&p).unwrap().inner.borrow(),
                    core.without_hydrogens_with_params(&p).unwrap()
                );
            }
            assert_eq!(*m.inner.borrow(), original);
            m.add_hydrogens_().unwrap();
            core.add_hydrogens_().unwrap();
            assert_eq!(*m.inner.borrow(), core);
            m.remove_hydrogens_().unwrap();
            core.remove_hydrogens_().unwrap();
            assert_eq!(*m.inner.borrow(), core);
            let p = ck::AddHsParams {
                explicit_only: true,
                ..Default::default()
            };
            m.add_hydrogens_with_params_(&p).unwrap();
            core.add_hydrogens_with_params_(&p).unwrap();
            assert_eq!(*m.inner.borrow(), core);
            let p = ck::RemoveHsParams {
                remove_isotopes: true,
                ..Default::default()
            };
            m.remove_hydrogens_with_params_(&p).unwrap();
            core.remove_hydrogens_with_params_(&p).unwrap();
            assert_eq!(*m.inner.borrow(), core);
        }
        let m = Molecule::from_smiles("CCO").unwrap();
        let before = m.inner.borrow().clone();
        let p = ck::AddHsParams {
            only_on_atoms: Some(vec![ck::AtomId::new(99)]),
            ..Default::default()
        };
        assert!(matches!(
            m.add_hydrogens_with_params_(&p),
            Err(ck::OperationError::Hydrogen(
                ck::HydrogenError::OnlyOnAtomOutOfRange { .. }
            ))
        ));
        assert_eq!(*m.inner.borrow(), before);
        assert!(m.with_hydrogens_with_params(&p).is_err());
        assert_eq!(*m.inner.borrow(), before);
    }
}
#[cfg(test)]
mod coordinate_tests {
    use super::*;
    #[test]
    fn hydrogens_added_3d_coordinates_transport_source_rows() {
        let sdf = r#"ethanol
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7000    0.2000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    1.2000    0.4000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
"#;
        let m = Molecule::from_sdf(sdf).unwrap();
        let core = ck::Molecule::from_sdf(sdf).unwrap();
        let p = ck::AddHsParams {
            add_coords: true,
            add_residue_info: true,
            skip_queries: true,
            ..Default::default()
        };
        let out = m.with_hydrogens_with_params(&p).unwrap();
        assert_eq!(
            *out.inner.borrow(),
            core.with_hydrogens_with_params(&p).unwrap()
        );
        assert_eq!(out.coordinates_3d(0).len(), 27);
        assert_eq!(m.coordinates_3d(0).len(), 9);
    }
}
