//! Thin canonical native 64-bit molecular hash projection.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn murcko_scaffold(&self) -> Result<Self, ck::OperationError> {
        self.inner.borrow().murcko_scaffold().map(|inner| Self {
            inner: std::cell::RefCell::new(inner),
        })
    }
    pub fn net_scaffold(&self) -> Result<Self, ck::OperationError> {
        self.inner.borrow().net_scaffold().map(|inner| Self {
            inner: std::cell::RefCell::new(inner),
        })
    }
    pub fn murcko_decompose(&self) -> Result<Self, ck::OperationError> {
        self.inner.borrow().murcko_decompose().map(|inner| Self {
            inner: std::cell::RefCell::new(inner),
        })
    }
    pub fn molecular_hash(&self) -> Result<u64, ck::MoleculeHashError> {
        // COSMolKit❗✔️: pub fn molecular_hash(&self) -> Result<u64, MoleculeHashError> {
        self.inner.borrow().molecular_hash()
    }
    pub fn molecular_hash_with_ranks(&self, ranks: &[u32]) -> Result<u64, ck::MoleculeHashError> {
        // COSMolKit❗✔️: pub fn molecular_hash_with_ranks(&self, ranks: &[u32]) -> Result<u64, MoleculeHashError> {
        self.inner.borrow().molecular_hash_with_ranks(ranks)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hashing_exact_u64_supplied_ranks_and_structured_errors() {
        for text in ["CCO", "c1ccccc1", "[NH4+]", "C[C@H](O)F"] {
            let m = Molecule::from_smiles(text).unwrap();
            let core = ck::Molecule::from_smiles(text).unwrap();
            assert_eq!(m.molecular_hash().unwrap(), core.molecular_hash().unwrap());
            for ranks in [
                vec![0; core.num_atoms()],
                (0..core.num_atoms()).map(|i| u32::MAX - i as u32).collect(),
            ] {
                assert_eq!(
                    m.molecular_hash_with_ranks(&ranks).unwrap(),
                    core.molecular_hash_with_ranks(&ranks).unwrap()
                );
            }
            assert_eq!(*m.inner.borrow(), core);
        }
        let m = Molecule::from_smiles("CCO").unwrap();
        for n in [0, 2, 4] {
            assert!(
                matches!(m.molecular_hash_with_ranks(&vec![0;n]),Err(ck::MoleculeHashError::RankCount{actual,atom_count:3}) if actual==n)
            );
        }
        assert!(matches!(
            Molecule::new().molecular_hash(),
            Err(ck::MoleculeHashError::EmptyMolecule)
        ));
        let raw = Molecule::from_smiles_with_params(
            "CCO",
            &ck::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(matches!(
            raw.molecular_hash(),
            Err(ck::MoleculeHashError::MissingPreparedValence)
        ));
        assert!(raw.molecular_hash_with_ranks(&[0, 1, 2]).is_ok());
        let mut builder = ck::MoleculeBuilder::new();
        builder.add_atom(ck::AtomSpec::new(ck::Element::C).with_atom_map(2147483648));
        let mapped = Molecule {
            inner: std::cell::RefCell::new(
                builder.build().unwrap().with_assigned_valence().unwrap(),
            ),
        };
        assert!(matches!(
            mapped.molecular_hash(),
            Err(ck::MoleculeHashError::CipRanks(
                ck::CipRankError::AtomMapOutOfRange {
                    map_number: 2147483648,
                    ..
                }
            ))
        ));
    }
}
