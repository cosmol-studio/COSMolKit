//! Shared topology value container.
//!
//! This type stores concrete graph rows and local adjacency/annotation data.
//! Topology remapping policies and live-molecule installation remain owned by
//! the runtime crate.

use crate::{AdjacencyList, Atom, Bond, StereoGroup, SubstanceGroup};
use crate::{AtomId, BondId};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum TopologyValidationError {
    #[error("atom at position {position} has id {id}, expected {position}")]
    AtomIdMismatch { position: usize, id: AtomId },
    #[error("bond at position {position} has id {id}, expected {position}")]
    BondIdMismatch { position: usize, id: BondId },
    #[error(
        "invalid bond endpoint in bond {bond}: {endpoint} atom {atom} is out of range for {atom_count} atoms"
    )]
    BondEndpointOutOfRange {
        bond: BondId,
        endpoint: &'static str,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("bond {bond} is a self-loop on atom {atom}")]
    SelfLoopBond { bond: BondId, atom: AtomId },
    #[error("bond {bond} stereo atoms {begin}-{end} are out of range for {atom_count} atoms")]
    StereoAtomOutOfRange {
        bond: BondId,
        begin: AtomId,
        end: AtomId,
        atom_count: usize,
    },
    #[error("bond {bond} stereo {stereo:?} requires stereo atom references")]
    StereoAtomsRequired {
        bond: BondId,
        stereo: crate::BondStereo,
    },
    #[error("substance group at position {position} has id {id:?}, expected {position}")]
    SubstanceGroupIdMismatch {
        position: usize,
        id: crate::SubstanceGroupId,
    },
    #[error(
        "substance group {sgroup:?} references atom {atom}, out of range for {atom_count} atoms"
    )]
    SubstanceGroupAtomOutOfRange {
        sgroup: crate::SubstanceGroupId,
        atom: AtomId,
        atom_count: usize,
    },
    #[error(
        "substance group {sgroup:?} references bond {bond}, out of range for {bond_count} bonds"
    )]
    SubstanceGroupBondOutOfRange {
        sgroup: crate::SubstanceGroupId,
        bond: BondId,
        bond_count: usize,
    },
    #[error("substance group {sgroup:?} has parent {parent:?} out of range")]
    SubstanceGroupParentOutOfRange {
        sgroup: crate::SubstanceGroupId,
        parent: crate::SubstanceGroupId,
    },
    #[error("stereo group references atom {atom}, out of range for {atom_count} atoms")]
    StereoGroupAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("stereo group references bond {bond}, out of range for {bond_count} bonds")]
    StereoGroupBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("adjacency does not match topology")]
    AdjacencyMismatch,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomMapping {
    pub old_to_new: Vec<Option<AtomId>>,
    pub new_to_old: Vec<Option<AtomId>>,
}

impl AtomMapping {
    #[must_use]
    pub fn old_to_new(&self) -> &[Option<AtomId>] {
        &self.old_to_new
    }
    #[must_use]
    pub fn new_to_old(&self) -> &[Option<AtomId>] {
        &self.new_to_old
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondMapping {
    pub old_to_new: Vec<Option<BondId>>,
    pub new_to_old: Vec<Option<BondId>>,
}

impl BondMapping {
    #[must_use]
    pub fn old_to_new(&self) -> &[Option<BondId>] {
        &self.old_to_new
    }
    #[must_use]
    pub fn new_to_old(&self) -> &[Option<BondId>] {
        &self.new_to_old
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyMapping {
    pub atoms: AtomMapping,
    pub bonds: BondMapping,
}

impl TopologyMapping {
    #[must_use]
    pub fn atoms(&self) -> &AtomMapping {
        &self.atoms
    }
    #[must_use]
    pub fn bonds(&self) -> &BondMapping {
        &self.bonds
    }
    #[must_use]
    pub fn retained_atom_indices(&self) -> Vec<usize> {
        self.atoms
            .new_to_old
            .iter()
            .filter_map(|id| id.map(AtomId::index))
            .collect()
    }
    pub fn identity(atom_count: usize, bond_count: usize) -> Self {
        Self {
            atoms: AtomMapping {
                old_to_new: (0..atom_count).map(|i| Some(AtomId::new(i))).collect(),
                new_to_old: (0..atom_count).map(|i| Some(AtomId::new(i))).collect(),
            },
            bonds: BondMapping {
                old_to_new: (0..bond_count).map(|i| Some(BondId::new(i))).collect(),
                new_to_old: (0..bond_count).map(|i| Some(BondId::new(i))).collect(),
            },
        }
    }
    pub fn with_appended(
        old_atoms: usize,
        old_bonds: usize,
        added_atoms: usize,
        added_bonds: usize,
    ) -> Self {
        let mut atom_new_to_old: Vec<_> = (0..old_atoms).map(|i| Some(AtomId::new(i))).collect();
        atom_new_to_old.extend((0..added_atoms).map(|_| None));
        let mut bond_new_to_old: Vec<_> = (0..old_bonds).map(|i| Some(BondId::new(i))).collect();
        bond_new_to_old.extend((0..added_bonds).map(|_| None));
        Self {
            atoms: AtomMapping {
                old_to_new: (0..old_atoms).map(|i| Some(AtomId::new(i))).collect(),
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: (0..old_bonds).map(|i| Some(BondId::new(i))).collect(),
                new_to_old: bond_new_to_old,
            },
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TopologyBlock {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub adjacency: AdjacencyList,
    pub substance_groups: Vec<SubstanceGroup>,
    pub stereo_groups: Vec<StereoGroup>,
}

impl Default for TopologyBlock {
    fn default() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: AdjacencyList::from_topology(0, &[]),
            substance_groups: Vec::new(),
            stereo_groups: Vec::new(),
        }
    }
}

impl TopologyBlock {
    pub fn validate(&self) -> Result<(), TopologyValidationError> {
        for (position, atom) in self.atoms.iter().enumerate() {
            if atom.id() != AtomId::new(position) {
                return Err(TopologyValidationError::AtomIdMismatch {
                    position,
                    id: atom.id(),
                });
            }
        }
        for (position, bond) in self.bonds.iter().enumerate() {
            if bond.id() != BondId::new(position) {
                return Err(TopologyValidationError::BondIdMismatch {
                    position,
                    id: bond.id(),
                });
            }
            for (endpoint, atom) in [("begin", bond.begin()), ("end", bond.end())] {
                if atom.index() >= self.atoms.len() {
                    return Err(TopologyValidationError::BondEndpointOutOfRange {
                        bond: bond.id(),
                        endpoint,
                        atom,
                        atom_count: self.atoms.len(),
                    });
                }
            }
            if bond.begin() == bond.end() {
                return Err(TopologyValidationError::SelfLoopBond {
                    bond: bond.id(),
                    atom: bond.begin(),
                });
            }
            if let Some([begin, end]) = bond.stereo_atoms()
                && (begin.index() >= self.atoms.len() || end.index() >= self.atoms.len())
            {
                return Err(TopologyValidationError::StereoAtomOutOfRange {
                    bond: bond.id(),
                    begin,
                    end,
                    atom_count: self.atoms.len(),
                });
            }
            if matches!(
                bond.stereo(),
                crate::BondStereo::Cis | crate::BondStereo::Trans
            ) && bond.stereo_atoms().is_none()
            {
                return Err(TopologyValidationError::StereoAtomsRequired {
                    bond: bond.id(),
                    stereo: bond.stereo(),
                });
            }
        }
        let atom_count = self.atoms.len();
        let bond_count = self.bonds.len();
        let sgroup_count = self.substance_groups.len();
        for (position, group) in self.substance_groups.iter().enumerate() {
            if group.id() != crate::SubstanceGroupId::new(position) {
                return Err(TopologyValidationError::SubstanceGroupIdMismatch {
                    position,
                    id: group.id(),
                });
            }
            for atom in group.atoms().iter().chain(group.parent_atoms()) {
                if atom.index() >= atom_count {
                    return Err(TopologyValidationError::SubstanceGroupAtomOutOfRange {
                        sgroup: group.id(),
                        atom: *atom,
                        atom_count,
                    });
                }
            }
            for point in group.attach_points() {
                for atom in std::iter::once(point.atom).chain(point.leaving_atom) {
                    if atom.index() >= atom_count {
                        return Err(TopologyValidationError::SubstanceGroupAtomOutOfRange {
                            sgroup: group.id(),
                            atom,
                            atom_count,
                        });
                    }
                }
            }
            for bond in group
                .bonds()
                .iter()
                .chain(group.cstates().iter().map(|state| &state.bond))
            {
                if bond.index() >= bond_count {
                    return Err(TopologyValidationError::SubstanceGroupBondOutOfRange {
                        sgroup: group.id(),
                        bond: *bond,
                        bond_count,
                    });
                }
            }
            if let Some(parent) = group.parent()
                && parent.index() >= sgroup_count
            {
                return Err(TopologyValidationError::SubstanceGroupParentOutOfRange {
                    sgroup: group.id(),
                    parent,
                });
            }
        }
        for group in &self.stereo_groups {
            for atom in group.atoms() {
                if atom.index() >= atom_count {
                    return Err(TopologyValidationError::StereoGroupAtomOutOfRange {
                        atom: *atom,
                        atom_count,
                    });
                }
            }
            for bond in group.bonds() {
                if bond.index() >= bond_count {
                    return Err(TopologyValidationError::StereoGroupBondOutOfRange {
                        bond: *bond,
                        bond_count,
                    });
                }
            }
        }
        let expected =
            AdjacencyList::try_from_topology(self.atoms.len(), &self.bonds).map_err(|error| {
                match error {
                    crate::AdjacencyError::BondAtomOutOfRange {
                        bond,
                        endpoint,
                        atom,
                        atom_count,
                    } => TopologyValidationError::BondEndpointOutOfRange {
                        bond,
                        endpoint,
                        atom,
                        atom_count,
                    },
                }
            })?;
        if self.adjacency != expected {
            return Err(TopologyValidationError::AdjacencyMismatch);
        }
        Ok(())
    }

    pub fn remove_atoms_with_mapping(&mut self, atoms_to_remove: &[AtomId]) -> TopologyMapping {
        let mut remove_atom = vec![false; self.atoms.len()];
        for atom in atoms_to_remove {
            if let Some(slot) = remove_atom.get_mut(atom.index()) {
                *slot = true;
            }
        }
        let mut atom_old_to_new = vec![None; self.atoms.len()];
        let mut atom_new_to_old = Vec::new();
        let mut atoms = Vec::with_capacity(self.atoms.len().saturating_sub(atoms_to_remove.len()));
        for atom in &self.atoms {
            if remove_atom[atom.id().index()] {
                continue;
            }
            let new_id = AtomId::new(atoms.len());
            atom_old_to_new[atom.id().index()] = Some(new_id);
            atom_new_to_old.push(Some(atom.id()));
            atoms.push(atom.clone().with_id(new_id));
        }
        let mut bond_old_to_new = vec![None; self.bonds.len()];
        let mut bond_new_to_old = Vec::new();
        let mut bonds = Vec::new();
        for bond in &self.bonds {
            let Some(begin) = atom_old_to_new.get(bond.begin().index()).and_then(|x| *x) else {
                continue;
            };
            let Some(end) = atom_old_to_new.get(bond.end().index()).and_then(|x| *x) else {
                continue;
            };
            let stereo_atoms = bond.stereo_atoms().and_then(|[left, right]| {
                Some([
                    atom_old_to_new.get(left.index()).and_then(|x| *x)?,
                    atom_old_to_new.get(right.index()).and_then(|x| *x)?,
                ])
            });
            let new_id = BondId::new(bonds.len());
            bond_old_to_new[bond.id().index()] = Some(new_id);
            bond_new_to_old.push(Some(bond.id()));
            bonds.push(bond.clone().remapped(new_id, begin, end, stereo_atoms));
        }
        let mut survives: Vec<_> = self
            .substance_groups
            .iter()
            .map(|sg| sg.can_remap_without_parent(&atom_old_to_new, &bond_old_to_new))
            .collect();
        loop {
            let mut changed = false;
            for idx in 0..self.substance_groups.len() {
                if survives[idx]
                    && self.substance_groups[idx]
                        .parent()
                        .is_some_and(|p| !survives.get(p.index()).copied().unwrap_or(false))
                {
                    survives[idx] = false;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        let mut sgroup_map = vec![None; self.substance_groups.len()];
        let mut next_sgroup_index = 0usize;
        for (idx, keep) in survives.iter().copied().enumerate() {
            if keep {
                sgroup_map[idx] = Some(crate::SubstanceGroupId::new(next_sgroup_index));
                next_sgroup_index += 1;
            }
        }
        self.substance_groups = self
            .substance_groups
            .iter()
            .enumerate()
            .filter_map(|(idx, sg)| {
                sgroup_map[idx]
                    .and_then(|id| sg.remapped(id, &atom_old_to_new, &bond_old_to_new, &sgroup_map))
            })
            .collect();
        self.stereo_groups = self
            .stereo_groups
            .iter()
            .filter_map(|g| g.remapped(&atom_old_to_new, &bond_old_to_new))
            .collect();
        self.atoms = atoms;
        self.bonds = bonds;
        self.adjacency = AdjacencyList::from_topology(self.atoms.len(), &self.bonds);
        TopologyMapping {
            atoms: AtomMapping {
                old_to_new: atom_old_to_new,
                new_to_old: atom_new_to_old,
            },
            bonds: BondMapping {
                old_to_new: bond_old_to_new,
                new_to_old: bond_new_to_old,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomSpec, BondOrder, BondSpec};

    fn atom(id: usize) -> Atom {
        Atom::from_spec(AtomId::new(id), AtomSpec::new(crate::Element::C))
    }

    #[test]
    fn validate_rejects_misaligned_ids_and_endpoints() {
        let block = TopologyBlock {
            atoms: vec![atom(1)],
            ..Default::default()
        };
        assert!(matches!(
            block.validate(),
            Err(TopologyValidationError::AtomIdMismatch { .. })
        ));

        let begin = AtomId::new(0);
        let end = AtomId::new(2);
        let bond = Bond::from_spec(BondId::new(0), BondSpec::new(begin, end, BondOrder::Single));
        let block = TopologyBlock {
            atoms: vec![atom(0)],
            bonds: vec![bond],
            adjacency: AdjacencyList::default(),
            ..Default::default()
        };
        assert!(matches!(
            block.validate(),
            Err(TopologyValidationError::BondEndpointOutOfRange { .. })
        ));
    }

    #[test]
    fn validate_rejects_stale_adjacency() {
        let begin = AtomId::new(0);
        let end = AtomId::new(1);
        let bond = Bond::from_spec(BondId::new(0), BondSpec::new(begin, end, BondOrder::Single));
        let block = TopologyBlock {
            atoms: vec![atom(0), atom(1)],
            bonds: vec![bond],
            adjacency: AdjacencyList::default(),
            ..Default::default()
        };
        assert_eq!(
            block.validate(),
            Err(TopologyValidationError::AdjacencyMismatch)
        );
    }
}
