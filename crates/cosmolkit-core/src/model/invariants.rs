use crate::{InvariantError, Molecule};

/// Check structural invariants that every `Molecule` must satisfy.
///
/// Future agents must add checks here when adding new molecule state. Adding a
/// state field without an invariant check is not allowed unless the human
/// author explicitly approves the exception.
#[cfg(feature = "runtime-invariants")]
pub(crate) fn check_molecule_invariants(molecule: &Molecule) -> Result<(), InvariantError> {
    for (row, atom) in molecule.atoms().iter().enumerate() {
        if atom.id().index() != row {
            return Err(InvariantError::AtomIndexMismatch { row, actual: atom.id() });
        }
    }

    let atom_count = molecule.num_atoms();
    for (row, bond) in molecule.bonds().iter().enumerate() {
        if bond.id().index() != row {
            return Err(InvariantError::BondIndexMismatch { row, actual: bond.id() });
        }
        if bond.begin() == bond.end() {
            return Err(InvariantError::SelfLoopBond {
                bond: bond.id(),
                atom: bond.begin(),
            });
        }
        if bond.begin().index() >= atom_count || bond.end().index() >= atom_count {
            return Err(InvariantError::InvalidBondEndpoint {
                bond: bond.id(),
                begin: bond.begin(),
                end: bond.end(),
                atom_count,
            });
        }
        if let Some([begin_ref, end_ref]) = bond.stereo_atoms()
            && (begin_ref.index() >= atom_count || end_ref.index() >= atom_count)
        {
            return Err(InvariantError::InvalidBondEndpoint {
                bond: bond.id(),
                begin: begin_ref,
                end: end_ref,
                atom_count,
            });
        }
        if matches!(bond.stereo(), crate::BondStereo::Cis | crate::BondStereo::Trans) && bond.stereo_atoms().is_none() {
            return Err(InvariantError::BondStereoAtomsRequired {
                bond: bond.id(),
                stereo: bond.stereo(),
            });
        }
    }

    for conformer in molecule.conformers_2d() {
        if conformer.coordinates().len() != atom_count {
            return Err(InvariantError::CoordinateRowCount {
                rows: conformer.coordinates().len(),
                atom_count,
            });
        }
    }

    for conformer in molecule.conformers_3d() {
        if conformer.coordinates().len() != atom_count {
            return Err(InvariantError::ConformerRowCount {
                conformer: conformer.id(),
                rows: conformer.coordinates().len(),
                atom_count,
            });
        }
    }

    let bond_count = molecule.num_bonds();
    let sgroup_count = molecule.substance_groups().len();
    for (row, substance_group) in molecule.substance_groups().iter().enumerate() {
        if substance_group.id().index() != row {
            return Err(InvariantError::SubstanceGroupIndexMismatch {
                row,
                actual: substance_group.id(),
            });
        }
        for atom in substance_group.atoms() {
            if atom.index() >= atom_count {
                return Err(InvariantError::InvalidSubstanceGroupAtom {
                    sgroup: substance_group.id(),
                    atom: *atom,
                    atom_count,
                });
            }
        }
        for atom in substance_group.parent_atoms() {
            if atom.index() >= atom_count {
                return Err(InvariantError::InvalidSubstanceGroupAtom {
                    sgroup: substance_group.id(),
                    atom: *atom,
                    atom_count,
                });
            }
        }
        for bond in substance_group.bonds() {
            if bond.index() >= bond_count {
                return Err(InvariantError::InvalidSubstanceGroupBond {
                    sgroup: substance_group.id(),
                    bond: *bond,
                    bond_count,
                });
            }
        }
        for attach_point in substance_group.attach_points() {
            if attach_point.atom.index() >= atom_count {
                return Err(InvariantError::InvalidSubstanceGroupAtom {
                    sgroup: substance_group.id(),
                    atom: attach_point.atom,
                    atom_count,
                });
            }
            if let Some(leaving_atom) = attach_point.leaving_atom
                && leaving_atom.index() >= atom_count
            {
                return Err(InvariantError::InvalidSubstanceGroupAtom {
                    sgroup: substance_group.id(),
                    atom: leaving_atom,
                    atom_count,
                });
            }
        }
        for cstate in substance_group.cstates() {
            if cstate.bond.index() >= bond_count {
                return Err(InvariantError::InvalidSubstanceGroupBond {
                    sgroup: substance_group.id(),
                    bond: cstate.bond,
                    bond_count,
                });
            }
        }
        if let Some(parent) = substance_group.parent()
            && parent.index() >= sgroup_count
        {
            return Err(InvariantError::InvalidSubstanceGroupParent {
                sgroup: substance_group.id(),
                parent,
            });
        }
    }

    for stereo_group in molecule.stereo_groups() {
        for atom in stereo_group.atoms() {
            if atom.index() >= atom_count {
                return Err(InvariantError::InvalidBondEndpoint {
                    bond: crate::BondId::new(0),
                    begin: *atom,
                    end: *atom,
                    atom_count,
                });
            }
        }
        for bond in stereo_group.bonds() {
            if bond.index() >= bond_count {
                return Err(InvariantError::InvalidSubstanceGroupBond {
                    sgroup: crate::SubstanceGroupId::new(0),
                    bond: *bond,
                    bond_count,
                });
            }
        }
    }

    Ok(())
}

#[cfg(feature = "runtime-invariants")]
pub(crate) fn enforce_molecule_invariants(molecule: &Molecule) -> Result<(), InvariantError> {
    check_molecule_invariants(molecule)
}

#[cfg(not(feature = "runtime-invariants"))]
pub(crate) fn enforce_molecule_invariants(_molecule: &Molecule) -> Result<(), InvariantError> {
    Ok(())
}
