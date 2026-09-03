//! Shared topology value container.
//!
//! This type stores concrete graph rows and local adjacency/annotation data.
//! Topology remapping policies and live-molecule installation remain owned by
//! the runtime crate.

use crate::{AdjacencyList, Atom, Bond, StereoGroup, SubstanceGroup};
use crate::{AtomId, BondId};

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
pub struct TopologyBlock<QA = (), QB = ()> {
    pub atoms: Vec<Atom<QA>>,
    pub bonds: Vec<Bond<QB>>,
    pub adjacency: AdjacencyList,
    pub substance_groups: Vec<SubstanceGroup>,
    pub stereo_groups: Vec<StereoGroup>,
}

impl<QA, QB> Default for TopologyBlock<QA, QB> {
    fn default() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: AdjacencyList::from_topology::<()>(0, &[]),
            substance_groups: Vec::new(),
            stereo_groups: Vec::new(),
        }
    }
}

impl<QA: Clone, QB: Clone> TopologyBlock<QA, QB> {
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
