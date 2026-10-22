//! Read-only ring searches; these queries never install molecular cache state.
use crate::{AtomId, BondId, Molecule, OperationError, RingSearchParams};

/// An owned ring-search result. Indices refer to the queried molecular version.
/// Symmetrized SSSR is not an enumeration of all simple cycles.
#[derive(Debug, Clone)]
pub struct RingSearchResult {
    inner: cosmolkit_core::RingInfo,
}

impl RingSearchResult {
    /// Ordered atom indices for each ring.
    pub fn atom_rings(&self) -> &[Vec<AtomId>] {
        self.inner.atom_rings()
    }
    /// Corresponding bond indices for each ring.
    pub fn bond_rings(&self) -> &[Vec<BondId>] {
        self.inner.bond_rings()
    }
    /// Number of rings in this result.
    pub fn num_rings(&self) -> usize {
        self.inner.num_rings()
    }
}

impl Molecule {
    /// Compute symmetrized SSSR without changing this molecule or its ring cache.
    pub fn symmetrized_sssr(&self) -> Result<RingSearchResult, OperationError> {
        self.symmetrized_sssr_with_params(&RingSearchParams::default())
    }

    /// Compute symmetrized SSSR with explicit bond inclusion options.
    pub fn symmetrized_sssr_with_params(
        &self,
        params: &RingSearchParams,
    ) -> Result<RingSearchResult, OperationError> {
        cosmolkit_core::symmetrized_sssr(self.topology(), params)
            .map(|inner| RingSearchResult { inner })
            .map_err(OperationError::Rings)
    }
}

#[cfg(all(test, feature = "cap-smiles"))]
mod tests {
    use super::*;

    #[test]
    fn symmetrized_sssr_is_an_owned_read_only_query() {
        let mol = Molecule::from_smiles("C12C3C4C1C5C2C3C45").unwrap();
        let before = mol.clone();
        let result = mol.symmetrized_sssr().unwrap();
        assert_eq!(result.num_rings(), 6); // Cubane has six symmetry-equivalent faces.
        assert!(result.atom_rings().iter().all(|ring| ring.len() == 4));
        assert!(result.bond_rings().iter().all(|ring| ring.len() == 4));
        let configured = mol
            .symmetrized_sssr_with_params(&RingSearchParams::default())
            .unwrap();
        assert_eq!(result.atom_rings(), configured.atom_rings());
        assert_eq!(result.bond_rings(), configured.bond_rings());
        assert_eq!(mol, before);
        drop(mol);
        assert_eq!(result.num_rings(), 6);
    }
}
