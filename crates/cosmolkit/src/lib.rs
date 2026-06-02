//! Facade crate that re-exports COSMolKit core modules.
//!
//! # Examples
//!
//! ```
//! use cosmolkit::Molecule;
//!
//! let mol = Molecule::from_smiles("CCO").unwrap();
//! let with_hydrogens = mol.with_hydrogens().unwrap();
//!
//! assert_eq!(mol.num_atoms(), 3);
//! assert!(with_hydrogens.num_atoms() > mol.num_atoms());
//! ```
//!
//! ```
//! use cosmolkit::Molecule;
//!
//! let mut mol = Molecule::from_smiles("CCO").unwrap();
//! mol.add_hydrogens_().unwrap();
//! mol.compute_2d_coordinates_().unwrap();
//!
//! assert!(mol.num_atoms() > 3);
//! assert_eq!(mol.coordinates_2d().unwrap().len(), mol.num_atoms());
//! ```

pub use cosmolkit_core as core;
pub use cosmolkit_core::bio;
pub use cosmolkit_core::io;
pub use cosmolkit_core::*;

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conformer_generation_public_api_is_reexported_from_facade() {
        let molecule = Molecule::from_smiles("CC").expect("ethane");

        let mut params = EmbedParameters::etkdg_v3();
        params.random_seed = 0xf00d;
        params.num_threads = 1;
        let generated = molecule
            .with_3d_conformer_with_params(params.clone())
            .expect("value-style conformer generation");
        assert!(molecule.conformers_3d().is_empty());
        assert_eq!(generated.conformers_3d().len(), 1);

        let direct_result =
            embed_molecule_result(&molecule, &mut params).expect("embed_molecule_result re-export");
        assert!(direct_result.ok());
        assert_eq!(direct_result.conf_id, 0);
        assert_eq!(direct_result.molecule.conformers_3d().len(), 1);

        let (embedded, conf_id) =
            embed_molecule(&molecule, &mut params).expect("direct embed_molecule re-export");
        assert_eq!(conf_id, 0);
        assert_eq!(embedded.conformers_3d().len(), 1);

        let mut multi_params = EmbedParameters::etkdg();
        multi_params.random_seed = 0xf00d;
        multi_params.num_threads = 1;
        let generated_multi = molecule
            .with_3d_conformers_with_params(2, multi_params.clone())
            .expect("value-style multi-conformer generation");
        assert_eq!(generated_multi.conformers_3d().len(), 2);

        let multi_result = embed_multiple_confs_result(&molecule, 2, &mut multi_params)
            .expect("embed_multiple_confs_result re-export");
        assert_eq!(multi_result.conf_ids, vec![0, 1]);
        assert_eq!(multi_result.generated_count(), 2);
        assert_eq!(multi_result.molecule.conformers_3d().len(), 2);

        let (embedded_multi, ids) = embed_multiple_confs(&molecule, 2, &mut multi_params)
            .expect("direct embed_multiple_confs re-export");
        assert_eq!(ids, vec![0, 1]);
        assert_eq!(embedded_multi.conformers_3d().len(), 2);
    }
}
