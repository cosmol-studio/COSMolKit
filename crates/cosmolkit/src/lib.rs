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
//! assert_eq!(mol.coords_2d().unwrap().len(), mol.num_atoms());
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
