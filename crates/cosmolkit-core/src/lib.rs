//! Core molecular graph and chemistry perception primitives.

pub mod bio;
pub mod io;
pub mod adjacency;
pub mod atom;
pub mod bond;
pub mod hydrogens;
pub mod kekulize;
pub mod molecule;
mod smiles;
pub mod valence;

pub use adjacency::{AdjacencyList, NeighborRef};
pub use atom::{Atom, ChiralTag};
pub use bond::{Bond, BondDirection, BondOrder, BondStereo};
pub use hydrogens::{AddHydrogensError, add_hydrogens_in_place};
pub use molecule::{Molecule, SmilesParseError};
pub use valence::{
    ValenceAssignment, ValenceError, ValenceModel, assign_radicals_rdkit_2025, assign_valence,
    rdkit_valence_list,
};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::{Atom, Bond, BondOrder, Molecule, version};

    #[test]
    fn version_is_not_empty() {
        assert!(!version().is_empty());
    }

    #[test]
    fn molecule_from_smiles_parses_basic_chain() {
        let mol = Molecule::from_smiles("CCO").expect("basic SMILES should parse");
        assert_eq!(mol.atomic_numbers(), vec![6, 6, 8]);
    }

    #[test]
    fn api_surface_types_exist() {
        let _atom = Atom {
            index: 0,
            atomic_num: 6,
            is_aromatic: false,
            formal_charge: 0,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: super::ChiralTag::Unspecified,
            isotope: None,
        };
        let _bond = Bond {
            index: 0,
            begin_atom: 0,
            end_atom: 1,
            order: BondOrder::Single,
            direction: super::BondDirection::None,
            stereo: super::BondStereo::None,
            stereo_atoms: Vec::new(),
        };
        let _mol = Molecule::default();
    }

    #[test]
    fn molecule_from_smiles_parses_bracket_charge() {
        let mol = Molecule::from_smiles("[N-]=[N+]=N").expect("charged atoms should parse");
        assert_eq!(mol.atomic_numbers(), vec![7, 7, 7]);
        assert_eq!(mol.atoms[0].formal_charge, -1);
        assert_eq!(mol.atoms[1].formal_charge, 1);
    }

    #[test]
    fn molecule_from_smiles_parses_aromatic_ring_order() {
        let mol = Molecule::from_smiles("O=c1cc(O)c2ccccc2o1").expect("aromatic ring should parse");
        assert_eq!(
            mol.atomic_numbers(),
            vec![8, 6, 6, 6, 8, 6, 6, 6, 6, 6, 6, 8]
        );
    }

    #[test]
    fn molecule_from_smiles_parses_fragment_dummy_and_dative_subset() {
        let mol = Molecule::from_smiles("[*:1]C.[NH3]->[Cu+2]<-[NH3]")
            .expect("representative special subset should parse");
        assert_eq!(mol.atomic_numbers(), vec![0, 6, 7, 29, 7]);
    }
}
