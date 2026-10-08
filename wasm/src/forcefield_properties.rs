//! Parameter coverage and detached MMFF property transport only.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn uff_has_all_molecule_params(&self) -> Result<bool, ck::UffParameterQueryError> {
        // COSMolKit❗✔️: pub fn uff_has_all_molecule_params(&self) -> Result<bool, UffParameterQueryError> {
        self.inner.borrow().uff_has_all_molecule_params()
    }
    pub fn mmff_has_all_molecule_params(&self) -> Result<bool, ck::MmffMolPropertiesError> {
        // COSMolKit❗✔️: pub fn mmff_has_all_molecule_params(&self) -> Result<bool, MmffMolPropertiesError> {
        self.inner.borrow().mmff_has_all_molecule_params()
    }
    pub fn mmff_properties(&self) -> Result<ck::MmffProperties, ck::MmffMolPropertiesError> {
        // COSMolKit❗✔️: pub fn mmff_properties(&self) -> Result<MmffProperties, MmffMolPropertiesError> {
        self.inner.borrow().mmff_properties()
    }
    pub fn mmff_properties_with_params(
        &self,
        params: &ck::MmffPropertiesParams,
    ) -> Result<ck::MmffProperties, ck::MmffMolPropertiesError> {
        // COSMolKit❗✔️: fn mmff_properties_with_params(
        self.inner.borrow().mmff_properties_with_params(params)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn forcefield_properties_all_queries_variants_atoms_and_errors() {
        for text in ["CCO", "c1ccncc1", "CC(=O)N"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            assert_eq!(
                m.uff_has_all_molecule_params().unwrap(),
                core.uff_has_all_molecule_params().unwrap()
            );
            assert_eq!(
                m.mmff_has_all_molecule_params().unwrap(),
                core.mmff_has_all_molecule_params().unwrap()
            );
            assert_eq!(
                m.mmff_properties().unwrap(),
                core.mmff_properties().unwrap()
            );
            for variant in ["MMFF94", "MMFF94s", "source-fallback"] {
                let p = ck::MmffPropertiesParams {
                    mmff_variant: variant.into(),
                };
                let a = m.mmff_properties_with_params(&p).unwrap();
                let b = core.mmff_properties_with_params(&p).unwrap();
                assert_eq!(a, b);
                for (i, row) in a.atoms().iter().enumerate() {
                    assert_eq!(a.atom_type(i).unwrap(), row.atom_type());
                    assert_eq!(a.formal_charge(i).unwrap(), row.formal_charge());
                    assert_eq!(a.partial_charge(i).unwrap(), row.partial_charge());
                }
                assert!(matches!(
                    a.atom_type(a.atoms().len()),
                    Err(ck::MmffMolPropertiesError::AtomIndexOutOfRange { .. })
                ));
            }
        }
        let raw = Molecule::from_smiles_with_params(
            "CCO",
            &ck::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(matches!(
            raw.uff_has_all_molecule_params(),
            Err(ck::UffParameterQueryError::Cache(_))
        ));
        assert!(
            Molecule::from_smiles("*")
                .unwrap()
                .mmff_properties()
                .is_ok()
        );
    }
}
