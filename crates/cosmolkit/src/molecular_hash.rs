//! Thin read-only projection of the sole native hashing owner.
use crate::{Molecule, MoleculeHashError};
impl Molecule {
    /// Compute the original native 64-bit hash in legacy CIP rank order.
    /// Requires already prepared valence; absent ordinary rings remain absent.
    pub fn molecular_hash(&self) -> Result<u64, MoleculeHashError> {
        let cache = self.derived_cache_runtime();
        cosmolkit_fingerprints::molecule_hash(
            self.topology(),
            cache.valence_assignment(),
            cache.valid_ring_info(),
        )
    }
    /// Compute the same native hash with exactly one supplied rank per atom.
    pub fn molecular_hash_with_ranks(&self, ranks: &[u32]) -> Result<u64, MoleculeHashError> {
        cosmolkit_fingerprints::molecule_hash_with_ranks(
            self.topology(),
            self.derived_cache_runtime().valid_ring_info(),
            ranks,
        )
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod tests {
    use super::*;
    fn sanitized_mol(smiles: &str) -> Molecule {
        // Canonical construction already runs sanitization and retains rings.
        Molecule::from_smiles(smiles).expect("valid SMILES")
    }
    #[test]
    fn test_mol_hash_benzene() {
        let mol = sanitized_mol("c1ccccc1");
        let hash = mol.molecular_hash().expect("hash");
        assert!(hash != 0);
    }

    #[test]
    fn test_mol_hash_deterministic() {
        let mol = sanitized_mol("c1ccccc1");
        let hash1 = mol.molecular_hash().expect("hash1");
        let hash2 = mol.molecular_hash().expect("hash2");
        assert_eq!(hash1, hash2);
    }

    #[test]
    fn test_mol_hash_different_molecules() {
        let mol1 = sanitized_mol("c1ccccc1");
        let mol2 = sanitized_mol("CCO");
        let hash1 = mol1.molecular_hash().expect("hash1");
        let hash2 = mol2.molecular_hash().expect("hash2");
        assert_ne!(hash1, hash2);
    }

    #[test]
    fn test_mol_hash_empty_error() {
        let mol = Molecule::new();
        assert!(mol.molecular_hash().is_err());
    }

    #[test]
    fn test_mol_hash_with_ranks() {
        let mol = sanitized_mol("c1ccccc1");
        let ranks = cosmolkit_core::assign_atom_cip_ranks(
            mol.topology(),
            mol.derived_cache_runtime().valence_assignment().unwrap(),
        )
        .expect("ranks");
        let hash = mol.molecular_hash_with_ranks(&ranks).expect("hash");
        assert!(hash != 0);
    }

    #[test]
    fn rank_length_and_missing_prepared_valence_are_structured() {
        let mol = sanitized_mol("CCO");
        for actual in [0, 2, 4] {
            assert!(
                matches!(mol.molecular_hash_with_ranks(&vec![0; actual]), Err(MoleculeHashError::RankCount {actual: n, atom_count: 3}) if n == actual)
            );
        }
        let mut builder = crate::MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let raw = builder.build().unwrap();
        assert!(matches!(
            raw.molecular_hash(),
            Err(MoleculeHashError::MissingPreparedValence)
        ));
        assert!(raw.molecular_hash_with_ranks(&[0]).is_ok());
    }
}
