//! Shared source-backed atom-property helpers.

use crate::chemistry::valence::{bond_valence_contrib, cached_valence_assignment};
use crate::{AtomId, Hybridization, Molecule};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub(crate) enum AtomPropertyError {
    #[error("atom index {atom} is outside the molecule")]
    AtomIndex { atom: AtomId },
    #[error("explicit valence has not been computed for atom {atom}")]
    MissingExplicitValence { atom: AtomId },
    #[error("explicit valence {explicit_valence} is below physical bond count {physical_bonds} for atom {atom}")]
    ExplicitValenceBelowPhysicalBonds {
        atom: AtomId,
        explicit_valence: i32,
        physical_bonds: u32,
    },
    #[error(transparent)]
    Valence(#[from] crate::ValenceError),
}

/// Return RDKit's AtomPair pi-electron count for an atom whose explicit
/// valence property cache has already been computed.
pub(crate) fn num_pi_electrons(molecule: &Molecule, atom_id: AtomId) -> Result<u32, AtomPropertyError> {
    // RDKit✔️✔️: unsigned int numPiElectrons(const Atom &atom) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   if (atom.getIsAromatic()) {
    // RDKit✔️✔️:     res = 1;
    // RDKit✔️✔️:   } else if (atom.getHybridization() != Atom::SP3) {
    // RDKit✔️✔️:     auto val =
    // RDKit✔️✔️:         static_cast<unsigned int>(atom.getValence(Atom::ValenceType::EXPLICIT));
    // RDKit✔️✔️:     unsigned int physical_bonds = atom.getNumExplicitHs();
    // RDKit✔️✔️:     const auto &mol = atom.getOwningMol();
    // RDKit✔️✔️:     for (const auto bond : mol.atomBonds(&atom)) {
    // RDKit✔️✔️:       if (bond->getValenceContrib(&atom) != 0.0) {
    // RDKit✔️✔️:         ++physical_bonds;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     CHECK_INVARIANT(val >= physical_bonds,
    // RDKit✔️✔️:                     "explicit valence exceeds atom degree");
    // RDKit✔️✔️:     res = val - physical_bonds;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let atom = molecule
        .atoms()
        .get(atom_id.index())
        .ok_or(AtomPropertyError::AtomIndex { atom: atom_id })?;
    if atom.is_aromatic() {
        return Ok(1);
    }
    if atom.hybridization() == Hybridization::Sp3 {
        return Ok(0);
    }

    let assignment =
        cached_valence_assignment(molecule).ok_or(AtomPropertyError::MissingExplicitValence { atom: atom_id })?;
    let explicit_valence = assignment
        .explicit_valence
        .get(atom_id.index())
        .copied()
        .ok_or(AtomPropertyError::MissingExplicitValence { atom: atom_id })?;
    if explicit_valence < 0 {
        return Err(AtomPropertyError::MissingExplicitValence { atom: atom_id });
    }

    num_pi_electrons_from_explicit_valence(molecule, atom_id, explicit_valence)
}

fn num_pi_electrons_from_explicit_valence(
    molecule: &Molecule,
    atom_id: AtomId,
    explicit_valence: i32,
) -> Result<u32, AtomPropertyError> {
    let atom = molecule
        .atoms()
        .get(atom_id.index())
        .ok_or(AtomPropertyError::AtomIndex { atom: atom_id })?;
    let mut physical_bonds = u32::from(atom.explicit_hydrogens());
    for neighbor in molecule.topology_block().adjacency.neighbors_of(atom_id.index()) {
        let bond = &molecule.bonds()[neighbor.bond.index()];
        if bond_valence_contrib(bond, atom_id)? != 0.0 {
            physical_bonds += 1;
        }
    }

    let explicit_valence_u32 = explicit_valence as u32;
    if explicit_valence_u32 < physical_bonds {
        return Err(AtomPropertyError::ExplicitValenceBelowPhysicalBonds {
            atom: atom_id,
            explicit_valence,
            physical_bonds,
        });
    }
    Ok(explicit_valence_u32 - physical_bonds)
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]

    use super::*;
    use crate::{AtomSpec, BondOrder, BondSpec, Element, Hybridization};

    #[test]
    fn source_port__atom__num_pi_electrons__line_934_aromatic_and_sp3_branches() {
        let aromatic = Molecule::from_smiles("c1ccccc1").unwrap();
        for atom in aromatic.atoms() {
            assert_eq!(num_pi_electrons(&aromatic, atom.id()), Ok(1));
        }

        let sp3 = Molecule::from_smiles("CC").unwrap();
        for atom in sp3.atoms() {
            assert_eq!(atom.hybridization(), Hybridization::Sp3);
            assert_eq!(num_pi_electrons(&sp3, atom.id()), Ok(0));
        }
    }

    #[test]
    fn source_port__atom__num_pi_electrons__line_934_multiple_bonds_and_explicit_hydrogens() {
        let carbon_dioxide = Molecule::from_smiles("O=C=O").unwrap();
        assert_eq!(num_pi_electrons(&carbon_dioxide, AtomId::new(1)), Ok(2));
        assert_eq!(num_pi_electrons(&carbon_dioxide, AtomId::new(0)), Ok(1));

        let methylene = Molecule::from_smiles("[CH2]=C").unwrap();
        assert_eq!(methylene.atoms()[0].explicit_hydrogens(), 2);
        assert_eq!(num_pi_electrons(&methylene, AtomId::new(0)), Ok(1));
    }

    #[test]
    fn source_port__atom__num_pi_electrons__line_934_zero_and_dative_valence_contributions() {
        let mut zero_builder = Molecule::builder();
        let zero_left = zero_builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
        let zero_right = zero_builder.add_atom(AtomSpec::new(Element::C));
        zero_builder
            .add_bond(BondSpec::new(zero_left, zero_right, BondOrder::Zero))
            .unwrap();
        let zero = zero_builder.build().unwrap();
        assert_eq!(num_pi_electrons_from_explicit_valence(&zero, zero_left, 0), Ok(0));

        let mut dative_builder = Molecule::builder();
        let donor = dative_builder.add_atom(AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2));
        let acceptor = dative_builder.add_atom(AtomSpec::new(Element::FE).with_hybridization(Hybridization::Sp2));
        dative_builder
            .add_bond(BondSpec::new(donor, acceptor, BondOrder::Dative))
            .unwrap();
        let dative = dative_builder.build().unwrap();
        assert_eq!(num_pi_electrons_from_explicit_valence(&dative, donor, 0), Ok(0));
        assert_eq!(num_pi_electrons_from_explicit_valence(&dative, acceptor, 1), Ok(0));
    }

    #[test]
    fn source_port__atom__num_pi_electrons__line_934_requires_valid_explicit_valence_state() {
        let mut builder = Molecule::builder();
        let left = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
        let right = builder.add_atom(AtomSpec::new(Element::C));
        builder.add_bond(BondSpec::new(left, right, BondOrder::Single)).unwrap();
        let molecule = builder.build().unwrap();

        assert_eq!(
            num_pi_electrons(&molecule, left),
            Err(AtomPropertyError::MissingExplicitValence { atom: left })
        );
        assert_eq!(
            num_pi_electrons_from_explicit_valence(&molecule, left, 0),
            Err(AtomPropertyError::ExplicitValenceBelowPhysicalBonds {
                atom: left,
                explicit_valence: 0,
                physical_bonds: 1,
            })
        );
    }
}
