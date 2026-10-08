//! Thin canonical ring-cache and family-cache operation projections.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn with_assigned_rings(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_assigned_rings()
        self.inner.borrow().with_assigned_rings().map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn assign_rings_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.assign_rings_()
        self.inner.borrow_mut().assign_rings_()
    }
    pub fn with_assigned_ring_families(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_assigned_ring_families()
        self.inner
            .borrow()
            .with_assigned_ring_families()
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn assign_ring_families_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.assign_ring_families_()
        self.inner.borrow_mut().assign_ring_families_()
    }
    pub fn with_assigned_ring_families_with_params(
        &self,
        params: &ck::RingSearchParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: self.inner.with_assigned_ring_families_with_params(params)
        self.inner
            .borrow()
            .with_assigned_ring_families_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn assign_ring_families_with_params_(
        &self,
        params: &ck::RingSearchParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: self.inner.assign_ring_families_with_params_(params)
        self.inner
            .borrow_mut()
            .assign_ring_families_with_params_(params)
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
            _ => panic!("changed ring result"),
        }
    }
    #[test]
    fn all_ring_and_family_methods_match_canonical_full_state_and_preserve_input() {
        for text in [
            "",
            "CCO",
            "C1CCCCC1",
            "c1ccccc1",
            "c1ccc2ccccc2c1",
            "C1CC2CCC1C2",
            "C1CC1.C1CCC1",
        ] {
            for remove_hydrogens in [false, true] {
                let parse = ck::SmilesParseParams {
                    sanitize: false,
                    remove_hydrogens,
                    ..Default::default()
                };
                let owner = ck::Molecule::from_smiles_with_params(text, &parse).unwrap();
                let m = Molecule::from_smiles_with_params(text, &parse).unwrap();
                let before = m.inner.borrow().clone();
                for family in [false, true] {
                    compare(
                        if family {
                            m.with_assigned_ring_families()
                        } else {
                            m.with_assigned_rings()
                        },
                        if family {
                            owner.with_assigned_ring_families()
                        } else {
                            owner.with_assigned_rings()
                        },
                    );
                    let mut expected = owner.clone();
                    let fresh = Molecule::from_smiles_with_params(text, &parse).unwrap();
                    let (a, b) = if family {
                        (
                            fresh.assign_ring_families_(),
                            expected.assign_ring_families_(),
                        )
                    } else {
                        (fresh.assign_rings_(), expected.assign_rings_())
                    };
                    assert_eq!(a.map_err(|e| e.to_string()), b.map_err(|e| e.to_string()));
                    assert_eq!(*fresh.inner.borrow(), expected);
                }
                for include_dative_bonds in [false, true] {
                    for include_hydrogen_bonds in [false, true] {
                        let p = ck::RingSearchParams {
                            include_dative_bonds,
                            include_hydrogen_bonds,
                        };
                        compare(
                            m.with_assigned_ring_families_with_params(&p),
                            owner.with_assigned_ring_families_with_params(&p),
                        );
                        let fresh = Molecule::from_smiles_with_params(text, &parse).unwrap();
                        let mut expected = owner.clone();
                        assert_eq!(
                            fresh
                                .assign_ring_families_with_params_(&p)
                                .map_err(|e| e.to_string()),
                            expected
                                .assign_ring_families_with_params_(&p)
                                .map_err(|e| e.to_string())
                        );
                        assert_eq!(*fresh.inner.borrow(), expected);
                    }
                }
                assert_eq!(*m.inner.borrow(), before);
            }
        }
    }
}
