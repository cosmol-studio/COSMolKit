//! Core molecular graph and chemistry perception primitives.

pub mod adjacency;
pub mod atom;
pub mod batch;
pub mod bio;
pub mod bond;
pub mod canon_smiles;
pub mod distgeom;
pub mod draw;
pub mod hydrogens;
pub mod io;
pub mod kekulize;
pub mod molecule;
mod periodic_table;
mod smiles;
pub mod smiles_write;
pub mod stereo;
pub mod valence;

pub use adjacency::{AdjacencyList, NeighborRef};
pub use atom::{Atom, ChiralTag};
pub use batch::{
    BatchErrorMode, BatchExportReport, BatchRecord, BatchRecordError, BatchValidationError,
    MoleculeBatch,
};
pub use bond::{Bond, BondDirection, BondOrder, BondStereo};
pub use distgeom::DgBoundsError;
pub use draw::{PreparedDrawAtom, PreparedDrawBond, PreparedDrawMolecule, SvgDrawError};
pub use hydrogens::{
    AddHydrogensError, RemoveHydrogensError, add_hydrogens_in_place, remove_hydrogens_in_place,
};
pub use molecule::{
    ConformerStore, CoordinateDimension, Molecule, PropertyStore, SmilesParseError,
    SmilesWriteError, TopologyData,
};
pub use smiles::assign_double_bond_stereo_from_directions;
pub use smiles_write::SmilesWriteParams;
pub use stereo::{LigandRef, TetrahedralStereo};
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
            atom_map_num: None,
            props: Default::default(),
            rdkit_cip_rank: None,
        };
        let _bond = Bond {
            index: 0,
            begin_atom: 0,
            end_atom: 1,
            order: BondOrder::Single,
            is_aromatic: false,
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
        assert_eq!(mol.atoms()[0].formal_charge, -1);
        assert_eq!(mol.atoms()[1].formal_charge, 1);
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

    #[test]
    fn molecule_from_smiles_does_not_aromatize_saturated_n_heterocycle() {
        let mol = Molecule::from_smiles("CN1CCCC1").expect("saturated ring should parse");
        assert!(
            mol.atoms().iter().all(|atom| !atom.is_aromatic),
            "saturated N-heterocycle atoms should not be aromatic"
        );
        assert!(
            mol.bonds().iter().all(|bond| !bond.is_aromatic),
            "saturated N-heterocycle bonds should not be aromatic"
        );
    }

    #[test]
    fn kekulize_in_place_with_clear_aromatic_flags_true_clears_aromatic_marks() {
        let mut mol = Molecule::from_smiles("c1ccccc1").expect("benzene should parse");
        super::kekulize::kekulize_in_place(&mut mol, true).expect("kekulize should succeed");
        assert!(
            mol.atoms().iter().all(|atom| !atom.is_aromatic),
            "clearAromaticFlags=true should clear aromatic atom flags"
        );
        assert!(
            mol.bonds().iter().all(|bond| !bond.is_aromatic),
            "clearAromaticFlags=true should remove aromatic bond order annotations"
        );
    }

    #[test]
    fn kekulize_in_place_with_clear_aromatic_flags_false_preserves_aromatic_marks() {
        let mut mol = Molecule::from_smiles("c1ccccc1").expect("benzene should parse");
        super::kekulize::kekulize_in_place(&mut mol, false).expect("kekulize should succeed");
        assert!(
            mol.atoms().iter().all(|atom| atom.is_aromatic),
            "clearAromaticFlags=false should preserve aromatic atom flags"
        );
        assert!(
            mol.bonds().iter().all(|bond| bond.is_aromatic),
            "clearAromaticFlags=false should preserve aromatic bond flags"
        );
        assert!(
            mol.bonds()
                .iter()
                .all(|bond| !matches!(bond.order, super::BondOrder::Aromatic)),
            "clearAromaticFlags=false should still rewrite aromatic bond types to single/double"
        );
    }

    #[test]
    fn molecule_to_smiles_writes_basic_chain() {
        let mol = Molecule::from_smiles("CCO").expect("SMILES parser should parse CCO");
        let smiles = mol
            .to_smiles(true)
            .expect("SMILES writer should serialize CCO");
        assert_eq!(smiles, "CCO");
    }
}
